#pragma once

#include <pybind11/pybind11.h>

#include "electrolytemodel.h"

namespace py = pybind11;

class PyElectrolyteModel : public ElectrolyteModel  {
public:
    /* Inherit the constructors */
    using ElectrolyteModel::ElectrolyteModel;

    /* Trampoline (need one for each virtual function) */
    std::vector<float> get_y(const std::vector<float> &conc, float T) override {
        PYBIND11_OVERRIDE_PURE(
            std::vector<float>, /* Return type */
            ElectrolyteModel,      /* Parent class */
            get_y,          /* Name of function in C++ (must match Python name) */
            conc, T      /* Argument(s) */
        );
    }

    std::vector<float> get_a(const std::vector<float> &conc, float T) override {
        PYBIND11_OVERRIDE_PURE(
            std::vector<float>, /* Return type */
            ElectrolyteModel,      /* Parent class */
            get_a,          /* Name of function in C++ (must match Python name) */
            conc, T      /* Argument(s) */
        );
    }
};

void linkElectrolyteModel(py::module& m)
{
    py::class_<ElectrolyteModel, PyElectrolyteModel>(m, "ElectrolyteModel")
        .def(py::init<>())
        .def("get_y", &ElectrolyteModel::get_y)
        .def("get_a", &ElectrolyteModel::get_a);
}