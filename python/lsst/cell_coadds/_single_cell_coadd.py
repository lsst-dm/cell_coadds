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

__all__ = ("SingleCellCoadd",)

from typing import TYPE_CHECKING, FrozenSet

from lsst.afw.image import ImageD, ImageF
from lsst.geom import Box2I

from ._common_components import CommonComponents, CommonComponentsProperties
from ._image_planes import ImagePlanes, OwnedImagePlanes, ViewImagePlanes
from .typing_helpers import ImageLike

if TYPE_CHECKING:
    from ._identifiers import CellIdentifiers, ObservationIdentifiers


class SingleCellCoadd(CommonComponentsProperties):
    """A single coadd cell, built only from input images that completely
    overlap that cell.

    Parameters
    ----------
    outer: `OwnedImagePlanes`
        The actual coadded images.
    psf : `ImageD`
        The coadded PSF image.
    inner_bbox: `Box2I`
        The bounding box of the inner region of this cell; must be disjoint
        with but adjacent to all other cell inner regions.
    inputs : `frozenset` of `ObservationIdentifiers`
        Identifiers of observations that contributed to this cell.
    common: `CommonComponents`
        Image attributes common to all cells in a patch.
    identifiers: `CellIdentifiers`
        Struct of identifiers for this cell.

    Notes
    -----
    At present we assume a single PSF image per cell is sufficient to capture
    spatial variability, which seems adequate given the results we have so far
    and the cell sizes we intend to use.
    """

    def __init__(
        self,
        outer: OwnedImagePlanes,
        *,
        psf: ImageD,
        inner_bbox: Box2I,
        inputs: FrozenSet[ObservationIdentifiers],
        common: CommonComponents,
        identifiers: CellIdentifiers,
    ):
        assert outer.getBBox().contains(
            inner_bbox
        ), f"Cell inner bbox {inner_bbox} is not contained by outer bbox {outer.getBBox()}."
        self._outer = outer
        self._psf = psf
        self._inner_bbox = inner_bbox
        self._inner = ViewImagePlanes(outer, bbox=inner_bbox, make_view=self._make_view)
        self._common = common
        self._inputs = inputs
        self._identifiers = identifiers

    @property
    def inner(self) -> ImagePlanes:
        """Image planes within the inner region of this cell that is disjoint
        with all other cell inner regions.
        """
        return self._inner

    @property
    def outer(self) -> ImagePlanes:
        """Image planes within the full outer region of this cell."""
        return self._outer

    @property
    def psf_image(self) -> ImageF:
        """The coadded PSF image."""
        return self._psf

    @property
    def inputs(self) -> FrozenSet[ObservationIdentifiers]:
        """Identifiers for the input images that contributed to this cell."""
        return self._inputs

    @property
    def identifiers(self) -> CellIdentifiers:
        """Struct of unique identifiers for this cell."""
        # This overrides the patch-level property from
        # CommonComponentsProperties to provide cell-level information.
        return self._identifiers

    @property
    def common(self) -> CommonComponents:
        # Docstring inherited.
        return self._common

    def _make_view(self, image: ImageLike) -> ImageLike:
        return image[self._inner_bbox]
