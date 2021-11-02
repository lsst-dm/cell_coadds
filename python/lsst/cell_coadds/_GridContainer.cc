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
#include "pybind11/stl.h"

#include "boost/format.hpp"

#include "lsst/cpputils/python.h"
#include "lsst/cell_coadds/python.h"
#include "lsst/cell_coadds/GridContainer.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace cell_coadds {

namespace {

void check_index(GridIndex const& index, _GridContainerCommon const& container) {
    if (index.x < container.get_offset().x ||
        (index.x - container.get_offset().x) >= container.get_shape().x) {
        throw py::index_error((boost::format("x index %s out of range; expected a value between %s and %s") %
                               index.x % container.get_offset().x % container.get_shape().x)
                                  .str());
    }
    if (index.y < container.get_offset().y ||
        (index.y - container.get_offset().y) >= container.get_shape().y) {
        throw py::index_error((boost::format("x index %s out of range; expected a value between %s and %s") %
                               index.y % container.get_offset().y % container.get_shape().y)
                                  .str());
    }
}

}  // namespace

void wrapGridContainer(utils::python::WrapperCollection& wrappers) {
    wrappers.wrapType(
        py::class_<GridContainerBuilder<py::object>>(wrappers.module, "GridContainerBuilder"),
        [](auto& mod, auto& cls) {
            cls.def(py::init<GridIndex const&>(), "shape"_a);
            cls.def(py::init<GridIndex const&, GridIndex const&>(), "shape"_a, "offset"_a);
            cls.def_property_readonly(
                "shape", &GridContainerBuilder<py::object>::get_shape, py::return_value_policy::copy);
            cls.def_property_readonly(
                "offset", &GridContainerBuilder<py::object>::get_offset, py::return_value_policy::copy);
            cls.def("__len__", &GridContainerBuilder<py::object>::size);
            cls.def(
                "__setitem__",
                [](GridContainerBuilder<py::object>& self, GridIndex const& index, py::object value) {
                    check_index(index, self);
                    self.set(index, std::move(value));
                });
            cls.def("finish", [](GridContainerBuilder<py::object> self) { return std::move(self).finish(); });
        });
    wrappers.wrapType(
        py::class_<GridContainer<py::object>>(wrappers.module, "GridContainer"), [](auto& mod, auto& cls) {
            cls.def(py::init<GridContainerBuilder<py::object>>(), "builder"_a);
            cls.def_property_readonly(
                "shape", &GridContainer<py::object>::get_shape, py::return_value_policy::copy);
            cls.def_property_readonly(
                "offset", &GridContainer<py::object>::get_offset, py::return_value_policy::copy);
            cls.def("__len__", &GridContainer<py::object>::size);
            cls.def_property_readonly("first", &GridContainer<py::object>::get_first);
            cls.def_property_readonly("last", &GridContainer<py::object>::get_last);
            cls.def(
                "__getitem__",
                [](GridContainer<py::object> const& self, GridIndex const& index) -> py::object {
                    check_index(index, self);
                    return self[index];
                });
            cls.def("__iter__", [](GridContainer<py::object> const& self) {
                return py::make_iterator(self.begin(), self.end());
            });
            cls.def("rebuild_empty", &GridContainer<py::object>::rebuild_empty<py::object>);
            cls.def("rebuild_transformed", [](GridContainer<py::object> self, py::object callable) {
                return std::move(self).rebuild_transformed(callable);
            });
            cls.def("__copy__", [](GridContainer<py::object> pass_by_value_copies) {
                return pass_by_value_copies;
            });
            cls.def("__deepcopy__", [](GridContainer<py::object> self, py::object memo) {
                py::object deepcopy = py::module::import("copy").attr("deepcopy");
                auto builder = std::move(self).rebuild_transformed(
                    [memo, deepcopy](py::object original) { return deepcopy(original, memo); });
                return std::move(builder).finish();
            });
        });
}

}  // namespace cell_coadds
}  // namespace lsst
