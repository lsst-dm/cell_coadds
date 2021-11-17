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
from abc import abstractmethod
from typing import Iterable, List

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT
from lsst.pipe.tasks.coaddBase import makeSkyInfo
import lsst.sphgeom
import lsst.utils

from ._single_cell_coadd import SingleCellCoadd
from ._multiple_cell_coadd import MultipleCellCoadd
from ._identifiers import *

__all__ = ("CoaddInCellsConnections", "CoaddInCellsConfig", "CoaddInCellsTask")


class CoaddInCellsConnections(pipeBase.PipelineTaskConnections,
                              dimensions=("tract", "patch", "band", "skymap"),
                              defaultTemplates={"inputCoaddName": "deep",
                                                "outputCoaddName": "deep",
                                                "warpType": "direct"}):
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
    inputWarps = cT.Input(
        doc=("Input list of warps to be assemebled i.e. stacked."
             "WarpType (e.g. direct, psfMatched) is controlled by the warpType config parameter"),
        name="{inputCoaddName}Coadd_{warpType}Warp",
        storageClass="ExposureF",
        dimensions=("tract", "patch", "skymap", "visit", "instrument"),
        deferLoad=True,
        multiple=True
    )
    skyMap = cT.Input(
        doc="Input definition of geometry/box and projection/wcs for coadded exposures",
        name="skyMap",
        storageClass="SkyMap",
        dimensions=("skymap",),
    )
    cellCoadd = cT.Output(
        doc="Coadded image",
        name="{outputCoaddName}CellCoadd",
        storageClass="MultipleCellCoadd",
        dimensions=("tract", "patch", "skymap", "band", "instrument")
    )


class CoaddInCellsConfig(pipeBase.PipelineTaskConfig,
                         pipelineConnections=CoaddInCellsConnections):
    """Configuration parameters for the `CoaddInCellsTask`.
    """
    cellIndices = pexConfig.ListField(
        dtype=int,
        doc="Cells to coadd; if set to an empty list, all cells are processed",
        default=[],
    )
    inputType = pexConfig.ChoiceField(
        doc="Type of input dataset",
        dtype=str,
        allowed={"inputWarps": "Use warps",
                 "calexps": "Use calexps"},
        default="calexps"
    )


class CoaddInCellsTask(pipeBase.PipelineTask):
    """Perform coaddition
    """

    ConfigClass = CoaddInCellsConfig
    _DefaultName = "cellCoadd"

    def runQuantum(self, butlerQC: pipeBase.ButlerQuantumContext,
                   inputRefs: pipeBase.InputQuantizedConnection,
                   outputRefs: pipeBase.OutputQuantizedConnection):
        """Construct warps and then coadds

        Notes
        -----

        PipelineTask (Gen3) entry point to warp. This method is analogous to
        `runDataRef`. See lsst.pipe.tasks.makeCoaddTempExp.py for comparison.
        """
        # Read in the inputs via the butler
        inputs = butlerQC.get(inputRefs)

        # Process the skyMap to WCS and other useful sky information
        skyMap = inputs["skyMap"]  # skyInfo below will contain this skyMap
        quantumDataId = butlerQC.quantum.dataId
        skyInfo = makeSkyInfo(
            skyMap, tractId=quantumDataId["tract"],
            patchId=quantumDataId["patch"],
        )
        packId = quantumDataId.pack("tract_patch_band")

        patchIdentifier = PatchIdentifiers(skymap=skyMap, tract=quantumDataId["tract"])
        # Run the warp and coaddition code
        multipleCellCoadd = self.run(inputs, skyInfo=skyInfo)

        # Persist the results via the butler
        butlerQC.put(multipleCellCoadd, outputRefs.cellCoadd)

    def run(self, inputs,
            skyInfo: pipeBase.Struct) -> MultipleCellCoadd:
        """Make coadd for all the cells
        """
        # Should all of these computation happen in runQuantum method itself
        # and let concrete classes implement a `run` method?
        if self.config.cellIndices:
            cellIndices = self.config.cellIndices
        else:
            cellIndices = range(skyInfo.patchInfo.getNumCells())

        expList = inputs[inputs["inputType"]]

        cellCoadds = []
        for cellInfo in skyInfo.patchInfo:
            # Select calexps that completely overlap with the
            bbox_list = self.select_overlaps(inputs["calexps"], cellInfo=cellInfo)

            if not bbox_list:
                continue
                # raise pipeBase.NoWorkFound("No exposures that completely overlap are found")
            cellCoadd = self.singleCellCoaddBuilder(expList, bbox_list, cellInfo)
            cellCoadds.append(cellCoadd)

        inner_bbox = None  # Placeholder for now
        return MultipleCellCoadd(cellCoadds, inner_bbox=inner_bbox)

    @abstractmethod
    def singleCellCoaddBuilder(self, calExpList: Iterable[lsst.afw.image.ExposureF],
                               bboxList: Iterable[lsst.geom.Box2I],
                               cellInfo: pipeBase.Struct) -> SingleCellCoadd:
        """Build a single-cell coadd

        The inner and outer bounding boxes of the cell could be obtained as
        cellInfo = skyInfo.patchInfo.getCellInfo(cellIndex)
        innerBBox = cellInfo.getInnerBBox()
        outerBBox = cellInfo.getOuterBBox()

        Parameters
        ----------
        calExpList : List[`lsst.afw.image.Exposure`]
            A list of input images to coadd.
        skyInfo : `lsst.pipe.base.Struct`
            Struct with geommetric information about the patches and cells.
        cellIndex : `int`
            Index to identify the cell

        Returns
        -------
        A `SingleCellCoadd` object
        """
        # raise NotImplementedError()
        ret = SingleCellCoadd(psf=coadd_psf, inner_bbox=cellInfo.getInnerBBox(), )

    @staticmethod
    def select_overlaps(explist, cellInfo) -> List[lsst.geom.Box2I]:
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
        cell_bbox = cellInfo.outerBBox  # TBI
        cell_wcs = cellInfo.wcs
        cell_corners = [cell_wcs.pixelToSky(corner.x, corner.y) for corner in cell_bbox.getCorners()]

        overlapping_bbox = []
        for exp in explist:
            calexp_bbox = exp.get(component="bbox")
            calexp_wcs = exp.get(component="wcs")
            # TODO: Use makeSkyPolygonFromBBox function
            calexp_corners = [calexp_wcs.pixelToSky(corner.x, corner.y) for corner in calexp_bbox.getCorners()]
            skyCell = lsst.sphgeom.ConvexPolygon([corner.getVector() for corner in cell_corners])
            skyCalexp = lsst.sphgeom.ConvexPolygon([corner.getVector() for corner in calexp_corners])

            if skyCell.isWithin(skyCalexp):
                tiny_bbox_min_corner = calexp_wcs.skyToPixel(cell_wcs.pixelToSky(cell_bbox.minX, cell_bbox.minY))
                tiny_bbox_max_corner = calexp_wcs.skyToPixel(cell_wcs.pixelToSky(cell_bbox.maxX, cell_bbox.maxY))
                tiny_bbox = lsst.geom.Box2D(minimum=tiny_bbox_min_corner, maximum=tiny_bbox_max_corner)
                tiny_bbox = lsst.geom.Box2I(tiny_bbox)
                overlapping_bbox.append(tiny_bbox)

        return overlapping_bbox
