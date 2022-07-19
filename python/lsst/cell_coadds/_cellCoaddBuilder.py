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
from typing import Any, Iterable, Mapping, Optional, Tuple

import lsst.geom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.base.connectionTypes as cT
import lsst.sphgeom
import lsst.utils
from lsst.daf.butler import DataCoordinate, DeferredDatasetHandle
from lsst.pipe.tasks.coaddBase import makeSkyInfo
from lsst.skymap import CellInfo, PatchInfo

from ._cell_coadds import UniformGrid
from ._common_components import CoaddUnits, CommonComponents
from ._identifiers import CellIdentifiers, GridIdentifiers, ObservationIdentifiers, PatchIdentifiers
from ._multiple_cell_coadd import MultipleCellCoadd
from ._single_cell_coadd import SingleCellCoadd

__all__ = (
    "MultipleCellCoaddBuilderConfig",
    "MultipleCellCoaddBuilderConnections",
    "MultipleCellCoaddBuilderTask",
    "SingleCellCoaddBuilderConfig",
    "SingleCellCoaddBuilderTask",
    "singleCellCoaddBuilderRegistry",
)


class SingleCellCoaddBuilderConfig(pexConfig.Config):
    """Configuration for a single-cell coadd builder."""

    psf_dimensions = pexConfig.Field[int](
        doc="Dimensions of the PSF image",
        dtype=int,
        default=41,
    )


class SingleCellCoaddBuilderTask(pipeBase.Task, metaclass=ABCMeta):
    """An abstract interface for tasks building coadds in cells.

    `SingleCellCoaddBuilderTask` is intended to serve as an abstract interface
    for various concrete implementation of coaddition algorithms via its `run`
    method. `MultipleCellCoaddBuilderTask` is the corresponding pipeline task
    that needs to be called from the pipeline. `MultipleCellCoaddBuilderTask`
    must contain a concrete implementation inherited from
    `SingleCellCoaddBuilderTask` as its base class and registered with the
    `singleCellCoaddBuilderTaskRegistry`, say using ``@registerConfigurable``
    decorator.

    See also
    --------
    MultipleCellCoaddBuilderTask

    """

    ConfigClass = SingleCellCoaddBuilderConfig
    _DefaultName = "singleCellCoaddBuilder"

    def __init__(self, name: str, config: SingleCellCoaddBuilderConfig, **kwargs: Any):
        self.config = config
        super().__init__(name=name, config=config, **kwargs)

    @abstractmethod
    def run(
        self,
        inputs: Mapping[ObservationIdentifiers, Tuple[DeferredDatasetHandle, lsst.geom.Box2I]],
        cellInfo: CellInfo,
    ) -> pipeBase.Struct:
        """Build a single-cell coadd

        The images passed in from `MultipleCellCoaddBuilderTask` are guaranteed
        to completely overlap the outer bounding box of the cells. Any further
        exclusion of images based on quality assessment or other criteria
        should be dome in this method.

        Parameters
        ----------
        inputs: `Mapping[ObservationIdentifiers, Tuple[DeferredDatasetHandle,
                                                       lsst.geom.Box2I]]`
            A mapping from `lsst.cell_coadds.ObservationIdentifiers`` to a
            tuple containing a `DeferredDatasetHandle` pointing to the input
            image (calexp or warps) and a minimal bounding box that can be read
            without loading the entire image.
        cellInfo: `lsst.skymap.CellInfo`
            An `lsst.skymap.CellInfo` object with the following
            attributes:
            - wcs: `lsst.afw.geom.SkyWcs`
            - outer_bbox: `lsst.geom.Box2I`

        Returns
        -------
        A pipeBase.Struct object with the following fields:
        image_planes, psf and inputs.
        """
        raise NotImplementedError()

    registry = pexConfig.makeRegistry(doc="Internal registry")


singleCellCoaddBuilderRegistry = pexConfig.makeRegistry(doc="Registry of single cell coadd builders")


class MultipleCellCoaddBuilderConnections(
    pipeBase.PipelineTaskConnections,
    dimensions=("tract", "patch", "band", "skymap"),
    defaultTemplates={"inputCoaddName": "deep", "outputCoaddName": "deep", "warpType": "direct"},
):
    # Since we defer loading of images, we could take in both calexp and
    # warps as inputs. Unless, we don't want a dependency on warps existing.
    # The type of image will be specified in MultipleCellCoaddConfig.
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

    cellCoadd = cT.Output(
        doc="Coadded image",
        name="{outputCoaddName}CellCoaddPickled",
        storageClass="MultipleCellCoadd",
        dimensions=("tract", "patch", "skymap", "band", "instrument"),
    )


class MultipleCellCoaddBuilderConfig(
    pipeBase.PipelineTaskConfig, pipelineConnections=MultipleCellCoaddBuilderConnections
):
    """Configuration parameters for the `MultipleCellCoaddBuilderTask`."""

    inputType = pexConfig.ChoiceField(
        doc="Type of input dataset",
        dtype=str,
        allowed={"inputWarps": "Use warps", "calexps": "Use calexps"},
        default="calexps",
    )

    psf_dimensions: pexConfig.Field = pexConfig.Field(
        doc="Dimensions of the PSF image",
        dtype=int,
        default=41,
    )

    singleCellCoaddBuilder = singleCellCoaddBuilderRegistry.makeField(
        doc="Coaddition algorithm to use. See `SingleCellCoaddBuilderTask` for details",
        optional=False,
    )


class MultipleCellCoaddBuilderTask(pipeBase.PipelineTask):
    """Task to build cell-based coadded images.

    This is the pipeline task that needs to be called from the pipeline. It
    contains ``singleCellCoaddBuilder`` as a subtask which must implement a
    concrete coaddition algorithm in its ``run`` method. The images to be
    coadded can either be of type ``calexp`` or ``warp``, as specified in the
    ``inputType`` configuration field. A ``skymap`` with ``cells`` tract
    builder must be passed. For each cell to be coadded, the task will query
    the butler for all the input images that completely overlap the cell's
    outer bounding box and passed them to the ``run`` method of the
    ``singleCellCoaddBuilder``.

    See also
    --------
    SingleCellCoaddBuilderTask
    """

    ConfigClass = MultipleCellCoaddBuilderConfig
    _DefaultName = "multipleCellCoaddBuilder"

    def __init__(self, **kwargs: Any):
        super().__init__(**kwargs)
        self.makeSubtask(name="singleCellCoaddBuilder")
        self.singleCellCoaddBuilder: SingleCellCoaddBuilderTask

    def runQuantum(  # type: ignore[override]
        self,
        butlerQC: pipeBase.ButlerQuantumContext,
        inputRefs: pipeBase.InputQuantizedConnection,
        outputRefs: pipeBase.OutputQuantizedConnection,
    ) -> MultipleCellCoadd:
        # Docstring inherited.
        self.config: MultipleCellCoaddBuilderConfig
        inputs: dict = butlerQC.get(inputRefs)
        skyMap = inputs.pop("skyMap")  # skyInfo below will contain this skyMap
        # Ideally, we should do this check early on.
        # But skyMap is not available during config time.
        if not skyMap.config.tractBuilder.name == "cells":
            raise TypeError("skyMap is not a cells skyMap")

        quantumDataId = butlerQC.quantum.dataId
        if quantumDataId is None:
            raise ValueError("quantumDataId is None")

        skyInfo = makeSkyInfo(
            skyMap,
            tractId=quantumDataId["tract"],
            patchId=quantumDataId["patch"],
        )

        # Run the (warp and) coaddition code
        multipleCellCoadd = self.run(
            inputs[self.config.inputType],
            skyInfo=skyInfo,
            quantumDataId=quantumDataId,
        )
        # import pickle
        # with open("/project/kannawad/multipleCellCoadd_tract3828_patch19_bandi.pkl", "rb") as fp:  # noqa: W505, E501
        #     multipleCellCoadd = pickle.load(fp)

        # Persist the results via the butler
        butlerQC.put(multipleCellCoadd, outputRefs.cellCoadd)
        return multipleCellCoadd

    def run(  # type: ignore[override]  # TODO: Remove after DM-34696 is fixed.
        self,
        expList: Iterable[DeferredDatasetHandle],
        skyInfo: pipeBase.Struct,
        quantumDataId: DataCoordinate,
    ) -> MultipleCellCoadd:
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

        cellCoadds: list[SingleCellCoadd] = []
        patchInfo: PatchInfo = skyInfo.patchInfo  # type: ignore[attr-defined]
        common = CommonComponents(
            units=CoaddUnits.nJy,
            wcs=patchInfo.wcs,
            band=quantumDataId.get("band", None),
            identifiers=PatchIdentifiers.from_data_id(quantumDataId),
        )

        # This for loop could/should be parallelized, may be with asyncio?
        for cellInfo in patchInfo:
            # Select calexps that completely overlap with the cell
            try:
                bbox_list = self._select_overlaps(expList, cellInfo=cellInfo, skyInfo=skyInfo)
            except KeyError:
                continue  # Exception handling is likely a relic, should be removed.

            if not bbox_list:
                continue  # We should raise exception here, since the task would fail anyway eventually.

            scc_inputs = {
                ObservationIdentifiers.from_data_id(handle.ref.dataId): (handle, bbox)
                for handle, bbox in zip(expList, bbox_list)
            }
            if len(scc_inputs) == 0:
                continue

            result = self.singleCellCoaddBuilder.run(scc_inputs, cellInfo)
            identifiers = CellIdentifiers(
                cell=GridIdentifiers.from_info(cellInfo),
                skymap=common.identifiers.skymap,
                tract=common.identifiers.tract,
                patch=common.identifiers.patch,
                band=common.identifiers.band,
            )
            # TODO: singleCellCoaddBuilder.run should return a SingleCellCoadd
            cellCoadd = SingleCellCoadd(
                outer=result.image_planes,  # type: ignore[attr-defined]
                psf=result.psf,  # type: ignore[attr-defined]
                inner_bbox=cellInfo.inner_bbox,
                inputs=result.inputs,  # type: ignore[attr-defined]
                common=common,
                identifiers=identifiers,
            )
            cellCoadds.append(cellCoadd)

        # grid has no notion about border or inner/outer boundaries.
        # So we have to clip the outermost border when constructing the grid.
        grid_bbox = patchInfo.outer_bbox.erodedBy(patchInfo.getCellBorder())
        grid = UniformGrid(grid_bbox, patchInfo.getCellInnerDimensions())

        _mCellCoadd = MultipleCellCoadd(
            cellCoadds,
            grid=grid,
            outer_cell_size=cellInfo.outer_bbox.getDimensions(),
            inner_bbox=None,
            common=common,
            psf_image_size=lsst.geom.Extent2I(
                self.config.psf_dimensions,
                self.config.psf_dimensions,
            ),
        )
        return _mCellCoadd

    @staticmethod
    def _select_overlaps(
        explist: Iterable[DeferredDatasetHandle],
        cellInfo: CellInfo,
        skyInfo: Optional[pipeBase.Struct] = None,
    ) -> Iterable[lsst.geom.Box2I]:
        """Filter exposures for cell-based coadds.

        This methods selects from a list of exposures/warps those images that
        completely overlap with the cell, thus enabling edgeless coadds.

        Parameters
        ----------
        explist: `list` [`~lsst.daf.butler.DeferredDatasetHandle`]
            List of handles for exposures to be coadded
        cellInfo: `pipeBase.Struct` or `collections.namedtuple`
            The cellInfo dict, must have .wcs and .outerBBox.

        Returns
        -------
        overlapping_bbox: `list` [`lsst.geom.Box2I`]
            List of bounding boxes for each image in `explist` that overlaps
            the given cell.
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
                if skyInfo.tractInfo.outer_sky_polygon.contains(skyCalexp):  # type: ignore[attr-defined]
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
