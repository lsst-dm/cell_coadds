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

__all__ = (
    "GridIdentifiers",
    "PatchIdentifiers",
    "CellIdentifiers",
    "ObservationIdentifiers",
)


from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional, Union

from lsst.skymap import CellInfo, Index2D, PatchInfo

if TYPE_CHECKING:
    from lsst.daf.butler import DataCoordinate

# TODO: Implement a way to unpack these as dicts


@dataclass(frozen=True)
class GridIdentifiers:
    """Struct of identifiers that identify an element in a grid.

    This may be used to represent patches within a tract, or cells within a
    patch.
    """

    sequential: int
    """Unique sequential integer ID for this element in the grid."""

    x: int
    """The column of this element in its grid.
    """

    y: int
    """The row of this element in its grid.
    """

    @property
    def index(self) -> Index2D:
        """The index as a 2-d (x, y) tuple."""
        return Index2D(x=self.x, y=self.y)

    @classmethod
    def from_data_id(cls, data_id: DataCoordinate) -> GridIdentifiers:
        """Construct a `GridIdentifier` for a patch from a data ID.

        Parameters
        ----------
        data_id : `lsst.daf.butler.DataCoordinate`
            Fully-expanded data ID that includes the 'patch' dimension.

        Returns
        -------
        grid_identifiers : `GridIdentifiers`
            Identifiers struct for the patch.
        """
        # The cell_x, cell_y below do not refer to cell as it is used
        # throughout the rest of this package; it's a historical butler thing.
        # The 'type: ignore' directives below are present because we require a
        # a fully-expanded data ID, but this isn't embedded in the types.
        return cls(
            sequential=data_id["patch"],  # type: ignore
            x=data_id.records["patch"].cell_x,  # type: ignore
            y=data_id.records["patch"].cell_y,  # type: ignore
        )

    @classmethod
    def from_info(cls, info: Union[PatchInfo, CellInfo]) -> GridIdentifiers:
        """Construct from a skymap `PatchInfo` or `CellInfo`.

        Parameters
        ----------
        info : `PatchInfo` or `CellInfo`
            Structure describing the patch or cell.

        Returns
        -------
        grid_identifiers : `GridIdentifiers`
            Identifiers struct for the patch or cell.
        """
        return cls(sequential=info.sequential_index, x=info.index.x, y=info.index.y)


@dataclass(frozen=True)
class PatchIdentifiers:
    """Struct of identifiers for a coadd patch."""

    skymap: str
    """The name of the skymap this patch belongs to.
    """

    tract: int
    """The name of the tract this patch belongs to.
    """

    patch: GridIdentifiers
    """Identifiers for the patch itself.
    """

    band: Optional[str]
    """Name of the band, if any.
    """

    @classmethod
    def from_data_id(cls, data_id: DataCoordinate) -> PatchIdentifiers:
        """Construct from a data ID.

        Parameters
        ----------
        data_id : `lsst.daf.butler.DataCoordinate`
            Fully-expanded data ID that includes the 'patch' dimension and
            optionally the `band` dimension.

        Returns
        -------
        identifiers : `PatchIdentifiers`
            Struct of identifiers for this patch.
        """
        # The 'type: ignore' directives below are present because we require a
        # a fully-expanded data ID, but this isn't embedded in the types.
        return cls(
            skymap=data_id["skymap"],  # type: ignore
            tract=data_id["tract"],  # type: ignore
            patch=GridIdentifiers.from_data_id(data_id),
            band=data_id.get("band"),
        )


@dataclass(frozen=True)
class CellIdentifiers(PatchIdentifiers):
    """Struct of identifiers for a coadd cell."""

    cell: GridIdentifiers
    """Identifiers for the cell itself."""


@dataclass(frozen=True)
class ObservationIdentifiers:
    """Struct of identifiers for an observation that contributed to a coadd
    cell."""

    instrument: str
    """Name of the instrument that this observation was taken with.
    """

    packed: int
    """ID that uniquely identifies both the visit and detector by packing
    together their IDs.
    """

    visit: int
    """Unique identifier for the visit.

    A visit may be comprised of more than one exposure only if all were
    observed back-to-back with no dithers, allowing them to be combined early
    the processing with no resampling.  All detector-level images in a visit
    share the same visit ID.
    """

    detector: int
    """Unique identifier for the detector.
    """

    @classmethod
    def from_data_id(cls, data_id: DataCoordinate) -> ObservationIdentifiers:
        """Construct from a data ID.

        Parameters
        ----------
        data_id : `lsst.daf.butler.DataCoordinate`
            Fully-expanded data ID that includes the 'visit' and 'detector'
            dimensions.

        Returns
        -------
        identifiers : `ObservationIdentifiers`
            Struct of identifiers for this observation.
        """
        # The 'type: ignore' directives below are present because we require a
        # a fully-expanded data ID, but this isn't embedded in the types.
        return cls(
            instrument=data_id["instrument"],  # type: ignore
            packed=data_id.pack("visit_detector"),  # type: ignore
            visit=data_id["visit"],  # type: ignore
            detector=data_id["detector"],  # type: ignore
        )
