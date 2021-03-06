// -*- LSST-C++ -*-
/*
 * This file is part of cell_coadds.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "pybind11/pybind11.h"
#include "lsst/cpputils/python.h"
#include "lsst/cell_coadds/python.h"

namespace lsst {
namespace cell_coadds {

void wrapUniformGrid(utils::python::WrapperCollection &);
void wrapStitchedPsf(utils::python::WrapperCollection &);
void wrapGridContainer(utils::python::WrapperCollection &);

PYBIND11_MODULE(_cell_coadds, mod) {
    utils::python::WrapperCollection wrappers(mod, "lsst.cell_coadds");
    wrappers.addInheritanceDependency("lsst.meas.algorithms");
    wrapUniformGrid(wrappers);
    wrapGridContainer(wrappers);
    wrapStitchedPsf(wrappers);
    wrappers.finish();
}

}  // namespace cell_coadds
}  // namespace lsst
