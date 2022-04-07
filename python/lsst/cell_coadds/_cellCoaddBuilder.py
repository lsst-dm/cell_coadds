# This file is part of cell_coadds.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import annotations

from abc import ABCMeta, abstractmethod
from typing import Iterable, Mapping, Tuple

import lsst.geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT
import lsst.sphgeom
import lsst.utils
from lsst.daf.butler import DeferredDatasetHandle
from lsst.pipe.tasks.coaddBase import makeSkyInfo

from ._cell_coadds import UniformGrid
from ._common_components import *
from ._identifiers import *
from ._multiple_cell_coadd import MultipleCellCoadd
from ._single_cell_coadd import SingleCellCoadd

#from .singleCellCoaddBuilder import SCCBuilder

__all__ = ("CoaddInCellsConnections", "CoaddInCellsConfig", "CoaddInCellsTask", "SingleCellCoaddBuilder")


class SingleCellCoaddBuilder(pipeBase.Task, metaclass=ABCMeta):
    ConfigClass = pexConfig.Config
    _DefaultName = "singleCellCoaddBuilder"

    def __init__(self, name: str, config: pexConfig.Config, **kwargs):
        self.config = config
        super().__init__(config=config)

    @abstractmethod
    def run(
        self,
        inputs: Mapping[ObservationIdentifiers, Tuple[DeferredDatasetHandle, lsst.geom.Box2I]],
        cellInfo: pipeBase.Struct,
    ) -> pipeBase.Struct[psf, image_planes, inputs]:

        """Build a single-cell coadd

        The inner and outer bounding boxes of the cell could be obtained as
        cellInfo = skyInfo.patchInfo.getCellInfo(cellIndex)
        innerBBox = cellInfo.getInnerBBox()
        outerBBox = cellInfo.getOuterBBox()

        Parameters
        ----------
        calExpList : Iterable[`lsst.afw.image.Exposure`]
            A list of input images to coadd.
        skyInfo : `lsst.pipe.base.Struct`
            Struct with geommetric information about the patches and cells.
        cellIndex : `int`
            Index to identify the cell

        Returns
        -------
        A `SingleCellCoadd` object
        """
        raise NotImplementedError()

    registry = pexConfig.makeRegistry(doc="Internal registry")

singleCellCoaddBuilderRegistry = pexConfig.makeRegistry(doc="Registry of single cell coadd builders")
#singleCellCoaddBuilderRegistry.register("singleCellCoaddBuilder", SCCBuilder)

class CoaddInCellsConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("tract", "patch", "band", "skymap"),
    defaultTemplates={"inputCoaddName": "deep", "outputCoaddName": "deep", "warpType": "direct"},
):
    # Since we defer loading of images, we could take in both calexp and
    # warps as inputs. Unless, we don't want a dependency on warps existing.
    # The type of image will be specified in CoaddInCellsConfig.
    calexps = cT.Input(
        doc="Input exposures to be resampled and optionally PSF-matched onto a SkyMap projection/patch",
        name="calexp",
        storageClass="ExposureF",
        dimensions=("instrument", "visit", "detector"),
        deferLoad=True,
        multiple=True,
    )
    # inputWarps = cT.Input(
    #    doc=(
    #        "Input list of warps to be assemebled i.e. stacked."
    #        "WarpType (e.g. direct, psfMatched) is controlled by the warpType config parameter"
    #    ),
    #    name="{inputCoaddName}Coadd_{warpType}Warp",
    #    storageClass="ExposureF",
    #    dimensions=("tract", "patch", "skymap", "visit", "instrument"),
    #    deferLoad=True,
    #    multiple=True,
    # )
    skyMap = cT.Input(
        doc="Input definition of geometry/box and projection/wcs for coadded exposures",
        name="skyMap",
        storageClass="SkyMap",
        dimensions=("skymap",),
    )
    #cellCoadd = cT.Output(
    #    doc="Coadded image",
    #    name="{outputCoaddName}CellCoadd",
    #    storageClass="MultipleCellCoadd",
    #    dimensions=("tract", "patch", "skymap", "band", "instrument"),
    #)


class CoaddInCellsConfig(pipeBase.PipelineTaskConfig, pipelineConnections=CoaddInCellsConnections):
    """Configuration parameters for the `CoaddInCellsTask`."""

    cellIndices = pexConfig.ListField(
        dtype=int,
        doc="Cells to coadd; if set to an empty list, all cells are processed",
        default=[],
    )
    inputType = pexConfig.ChoiceField(
        doc="Type of input dataset",
        dtype=str,
        allowed={"inputWarps": "Use warps", "calexps": "Use calexps"},
        default="calexps",
    )

    singleCellCoaddBuilder = singleCellCoaddBuilderRegistry.makeField(doc="", default="sccBuilder", optional=True)
    #singleCellCoaddBuilder = pexConfig.ConfigurableField(doc="Concrete", target=SingleCellCoaddBuilder,
    #                                                     ConfigClass=SingleCellCoaddBuilder.ConfigClass)

class CoaddInCellsTask(pipeBase.PipelineTask):
    """Perform coaddition"""

    ConfigClass = CoaddInCellsConfig
    _DefaultName = "cellCoadd"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        # singleCellCoaddBuilder: pipeBase.Task
        self.makeSubtask(name="singleCellCoaddBuilder")

    def runQuantum(
        self,
        butlerQC: pipeBase.ButlerQuantumContext,
        inputRefs: pipeBase.InputQuantizedConnection,
        outputRefs: pipeBase.OutputQuantizedConnection,
    ):
        """Construct warps and then coadds

        Notes
        -----

        PipelineTask (Gen3) entry point to warp. This method is analogous to
        `runDataRef`. See lsst.pipe.tasks.makeCoaddTempExp.py for comparison.
        """
        # Read in the inputs via the butler
        inputs = butlerQC.get(inputRefs)
        # Process the skyMap to WCS and other useful sky information
        skyMap0 = inputs["skyMap"]  # skyInfo below will contain this skyMap
        if not skyMap0.config.tractBuilder.name == "cells":
            # The method below is a hack
            # It should throw an exception
            skyMap = self._convert_skyMaps(skyMap0)
        else:
            skyMap = skyMap0
        quantumDataId = butlerQC.quantum.dataId
        skyInfo = makeSkyInfo(
            skyMap,
            tractId=quantumDataId["tract"],
            patchId=quantumDataId["patch"],
        )
        skyInfo0 = makeSkyInfo(
            skyMap0,
            tractId=quantumDataId["tract"],
            patchId=quantumDataId["patch"],
        )
        packId = quantumDataId.pack("tract_patch_band")

        patchIdentifier = PatchIdentifiers(skymap=skyMap, tract=quantumDataId["tract"], patch=quantumDataId["patch"])
        # Run the warp and coaddition code
        multipleCellCoadd = self.run(inputs, skyInfo=skyInfo, quantumDataId=quantumDataId, skyInfo0=skyInfo0)

        # Persist the results via the butler
        # TODO: We cannot persist this until DM-32691 is done.
        # butlerQC.put(multipleCellCoadd, outputRefs.cellCoadd)
        return multipleCellCoadd

    def run(self, inputs, skyInfo: pipeBase.Struct, quantumDataId, skyInfo0=None) -> MultipleCellCoadd:
        """Make coadd for all the cells"""
        # Should all of these computation happen in runQuantum method itself
        # and let concrete classes implement a `run` method?
        if self.config.cellIndices:
            nx, ny = self.config.cellIndices
        else:
            nx, ny = skyInfo.patchInfo.getNumCells()

        cellIndices = range(nx)

        expList = inputs[self.config.inputType]

        cellCoadds = []
        common = CommonComponents(
            units=CoaddUnits.nJy,
            wcs=skyInfo.patchInfo.wcs,
            band=quantumDataId["band"],
            identifiers=PatchIdentifiers.from_data_id(quantumDataId),
        )

        for cellInfo in skyInfo.patchInfo:
            # Select calexps that completely overlap with the
            try:
                bbox_list = self._select_overlaps(inputs["calexps"], cellInfo=cellInfo, skyInfo=skyInfo0)
            except KeyError:
                import pdb; pdb.set_trace()
                continue

            if not bbox_list:
                continue
                # raise pipeBase.NoWorkFound("No exposures that completely overlap are found")
            scc_inputs = {
                ObservationIdentifiers.from_data_id(handle.ref.dataId): (handle, bbox)
                for handle, bbox in zip(inputs["calexps"], bbox_list)
            }
            if len(scc_inputs) == 0:
                continue

            result = self.singleCellCoaddBuilder.run(scc_inputs, cellInfo)
            cellCoadd = SingleCellCoadd(
                outer=result.image_planes,
                psf=result.psf,
                inner_bbox=cellInfo.inner_bbox,
                inputs=result.inputs,
                common=common,
                identifiers=CellIdentifiers(
                    cell=GridIdentifiers.from_info(cellInfo),
                    skymap=common.identifiers.skymap,
                    tract=common.identifiers.tract,
                    patch=common.identifiers.patch,
                    #  **common.identifiers.asdict()  # Desired
                ),
            )
            cellCoadds.append(cellCoadd)

        inner_bbox = None  # Placeholder for now
        grid = UniformGrid(skyInfo.patchInfo.inner_bbox, cellInfo.outer_bbox.getDimensions())
        return MultipleCellCoadd(cellCoadds, grid=grid, outer_cell_size=skyInfo.bbox,
                                 inner_bbox=inner_bbox, common=common, psf_image_size=41)

    @staticmethod
    def _select_overlaps(explist, cellInfo, skyInfo=None) -> List[lsst.geom.Box2I]:
        """Filter exposures for cell-based coadds.

        This methods selects from a list of exposures/warps those images that
        completely overlap with the cell, thus enabling edgeless coadds.

        Parameters
        ----------
        explist: `list` [`ExposureF`]
            List of exposures to be coadded
        cellInfo: `dict`
            The cellInfo dict, must have .wcs and .outerBBox.

        Returns
        -------
        edgeless_explist: `list` of exposures/warps to use
        """
        cell_bbox = cellInfo.outer_bbox
        cell_wcs = cellInfo.wcs
        cell_corners = [cell_wcs.pixelToSky(corner.x, corner.y) for corner in cell_bbox.getCorners()]

        overlapping_bbox = []
        for exp in explist:
            calexp_bbox = exp.get(component="bbox")
            calexp_wcs = exp.get(component="wcs")
            # TODO: Use makeSkyPolygonFromBBox function
            calexp_corners = [
                calexp_wcs.pixelToSky(corner.x, corner.y) for corner in calexp_bbox.getCorners()
            ]
            skyCell = lsst.sphgeom.ConvexPolygon([corner.getVector() for corner in cell_corners])
            skyCalexp = lsst.sphgeom.ConvexPolygon([corner.getVector() for corner in calexp_corners])

            if skyInfo:
                if skyInfo.tractInfo.outer_sky_polygon.contains(skyCalexp):
                    pass
            if skyCell.isWithin(skyCalexp):
                tiny_bbox_min_corner = calexp_wcs.skyToPixel(
                    cell_wcs.pixelToSky(cell_bbox.minX, cell_bbox.minY)
                )
                tiny_bbox_max_corner = calexp_wcs.skyToPixel(
                    cell_wcs.pixelToSky(cell_bbox.maxX, cell_bbox.maxY)
                )
                tiny_bbox = lsst.geom.Box2D(minimum=tiny_bbox_min_corner, maximum=tiny_bbox_max_corner)
                tiny_bbox = lsst.geom.Box2I(tiny_bbox)
                overlapping_bbox.append(tiny_bbox)

        return overlapping_bbox

    @staticmethod
    def _convert_skyMaps(skyMap):
        """Hacky method to change the tractBuilder"""
        import copy
        newConfig = copy.copy(skyMap.config)
        newConfig.tractBuilder.name = "cells"
        newSkyMap = skyMap.__class__(newConfig)
        return newSkyMap
