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
from typing import Iterable, Mapping, Tuple, Type

from lsst.afw.image import Exposure
from lsst.daf.butler import DataCoordinate
import lsst.geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT
import lsst.sphgeom
import lsst.utils
from lsst.daf.butler import DeferredDatasetHandle
from lsst.pipe.tasks.coaddBase import makeSkyInfo

from ._cell_coadds import GridContainerBuilder, UniformGrid
from ._common_components import *
from ._identifiers import *
from ._multiple_cell_coadd import MultipleCellCoadd
from ._single_cell_coadd import SingleCellCoadd


__all__ = ("MultipleCellsCoaddBuilderConfig", "MultipleCellsCoaddBuilderTask",
           "SingleCellCoaddBuilderTask")  # "SingleCellCoaddBuilderConfig",


class SingleCellCoaddBuilderConfig(pexConfig.Config):
    """Configuration for a single-cell coadd builder."""

    psf_dimensions = pexConfig.Field(
        doc="Dimensions of the PSF image",
        dtype=int,
        default=41,
    )


class SingleCellCoaddBuilderTask(pipeBase.Task, metaclass=ABCMeta):
    ConfigClass = SingleCellCoaddBuilderConfig
    _DefaultName = "singleCellCoaddBuilder"

    def __init__(self, name: str, config: SingleCellCoaddBuilderConfig, **kwargs):
        self.config = config
        super().__init__(name=name, config=config, **kwargs)

    @abstractmethod
    def run(
        self,
        inputs: Mapping[ObservationIdentifiers, Tuple[DeferredDatasetHandle, lsst.geom.Box2I]],
        cellInfo: pipeBase.Struct,
    ) -> pipeBase.Struct['psf', 'image_planes', 'inputs']:
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
        A pipeBase.Struct object with the following fields:
        psf, image_planes and inputs.
        """
        raise NotImplementedError()

    registry = pexConfig.makeRegistry(doc="Internal registry")

singleCellCoaddBuilderRegistry = pexConfig.makeRegistry(doc="Registry of single cell coadd builders")


class MultipleCellsCoaddBuilderConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("tract", "patch", "band", "skymap"),
    defaultTemplates={"inputCoaddName": "deep", "outputCoaddName": "deep", "warpType": "direct"},
):
    # Since we defer loading of images, we could take in both calexp and
    # warps as inputs. Unless, we don't want a dependency on warps existing.
    # The type of image will be specified in MultipleCellsCoaddConfig.
    calexps = cT.Input(
        doc="Input exposures to be resampled and optionally PSF-matched onto a SkyMap projection/patch",
        name="calexp",
        storageClass="ExposureF",
        dimensions=("instrument", "visit", "detector"),
        deferLoad=True,
        multiple=True,
    )

    skyMap = cT.Input(
        doc="Input definition of geometry/box and projection/wcs for coadded exposures",
        name="skyMap",
        storageClass="SkyMap",
        dimensions=("skymap",),
    )

    # TODO: DM-32691 Uncomment this after we can serialize cell-based coadds.
    #cellCoadd = cT.Output(
    #    doc="Coadded image",
    #    name="{outputCoaddName}CellCoadd",
    #    storageClass="MultipleCellCoadd",
    #    dimensions=("tract", "patch", "skymap", "band", "instrument"),
    #)


class MultipleCellsCoaddBuilderConfig(pipeBase.PipelineTaskConfig,
                                      pipelineConnections=MultipleCellsCoaddBuilderConnections):
    """Configuration parameters for the `MultipleCellsCoaddTask`."""

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

    psf_dimensions = pexConfig.Field(
        doc="Dimensions of the PSF image",
        dtype=int,
        default=41,
    )

    singleCellCoaddBuilder = singleCellCoaddBuilderRegistry.makeField(doc="", default="sccBuilder", optional=True)


class MultipleCellsCoaddBuilderTask(pipeBase.PipelineTask):
    """Perform coaddition"""

    ConfigClass = MultipleCellsCoaddBuilderConfig
    _DefaultName = "multipleCellCoaddBuilder"

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.makeSubtask(name="singleCellCoaddBuilder")
        self.singleCellCoaddBuilder: SingleCellCoaddBuilderTask

    def runQuantum(
        self,
        butlerQC: pipeBase.ButlerQuantumContext,
        inputRefs: pipeBase.InputQuantizedConnection,
        outputRefs: pipeBase.OutputQuantizedConnection,
    ):
        # Docstring inherited.
        inputs: dict = butlerQC.get(inputRefs)
        skyMap = inputs.pop("skyMap")  # skyInfo below will contain this skyMap
        # Ideally, we should do this check early on.
        # But skyMap is not available during config time.
        if not skyMap.config.tractBuilder.name == "cells":
            raise TypeError("skyMap is not a cells skyMap")

        quantumDataId = butlerQC.quantum.dataId
        skyInfo = makeSkyInfo(
            skyMap,
            tractId=quantumDataId["tract"],
            patchId=quantumDataId["patch"],
        )

        # Run the (warp and) coaddition code
        multipleCellCoadd = self.run(inputs[self.config.inputType],
                                     skyInfo=skyInfo,
                                     quantumDataId=quantumDataId
                                     )

        # Persist the results via the butler
        # TODO: We cannot persist this until DM-32691 is done.
        # butlerQC.put(multipleCellCoadd, outputRefs.cellCoadd)
        return multipleCellCoadd

    def run(self, expList: Iterable[DeferredDatasetHandle], skyInfo: pipeBase.Struct,
            quantumDataId: DataCoordinate) -> MultipleCellCoadd:
        """Run coaddition algorithm for all the cells.

        Parameters
        ----------
        expList: `list` of `lsst.daf.butler.DeferredDatasetHandle`
            An iterable of `lsst.daf.butler.DeferredDatasetHandle` objects,
            where the objects can be either calexp or warp images.
        skyInfo: `pipeBase.Struct`
            Struct with geommetric information about the patches and cells.
        quantumDataId: `DataCoordinate`
            An immutable dataID dictionary that uniquely refers to the dataset.

        Returns
        -------
        multipleCellCoadd: `MultipleCellCoadd`
            Cell-based coadded image.
        """

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
                bbox_list = self._select_overlaps(expList,
                                                  cellInfo=cellInfo,
                                                  skyInfo=skyInfo)
            except KeyError:
                continue

            if not bbox_list:
                continue
                # raise pipeBase.NoWorkFound("No exposures that completely overlap are found")
            scc_inputs = {
                ObservationIdentifiers.from_data_id(handle.ref.dataId): (handle, bbox)
                for handle, bbox in zip(expList, bbox_list)
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
                    band=common.identifiers.band,
                ),
            )
            cellCoadds.append(cellCoadd)

        inner_bbox = None  # Placeholder for now
        innerBBox = skyInfo.patchInfo.inner_bbox.dilatedBy(cellInfo.inner_bbox.getDimensions())  # dilatedBy is HACK
        #grid = UniformGrid(cellInfo.inner_bbox, cellInfo.inner_bbox.getDimensions())
        grid = UniformGrid(innerBBox, innerBBox.getDimensions())
        builder = GridContainerBuilder(grid.shape)
        # import pdb; pdb.set_trace()
        #builder[Index2D(x=0, y=0)] = cellCoadd
        #_mCellCoadd = builder.finish()
        #return _mCellCoadd
        #innerBBox = skyInfo.patchInfo.inner_bbox.dilatedBy(cellInfo.inner_bbox.getDimensions())  # dilatedBy is HACK
        grid = UniformGrid(innerBBox, cellInfo.inner_bbox.getDimensions()) ## Original, outer-> inner
        _mCellCoadd = MultipleCellCoadd(cellCoadds,
                                        grid=grid,
                                        outer_cell_size=cellInfo.outer_bbox.getDimensions(),
                                        inner_bbox=inner_bbox,
                                        common=common,
                                        psf_image_size=lsst.geom.Extent2I(self.config.psf_dimensions,
                                                                          self.config.psf_dimensions),
                                        )
        return _mCellCoadd

    @staticmethod
    def _select_overlaps(explist: Iterable[DeferredDatasetHandle],
                         cellInfo,
                         skyInfo=None) -> Iterable[type[lsst.geom.Box2I]]:
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
        overlapping_bbox: `list` of bounding boxes for each image in `explist`
            that overalps the given cell.
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
