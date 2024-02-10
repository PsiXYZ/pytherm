#pragma once

#include <pybind11/pybind11.h>
#include <pytherm/activity/activitymodel.h>

namespace py = pybind11;

class PyActivityModel : public ActivityModel  {
public:
    /* Inherit the constructors */
    using ActivityModel::ActivityModel;

    /* Trampoline (need one for each virtual function) */
    std::vector<float> get_y(const std::vector<float> &conc, float T) override {
        PYBIND11_OVERRIDE_PURE(
            std::vector<float>, /* Return type */
            ActivityModel,      /* Parent class */
            get_y,          /* Name of function in C++ (must match Python name) */
            conc, T      /* Argument(s) */
        );
    }

    std::vector<float> get_a(const std::vector<float> &conc, float T) override {
        PYBIND11_OVERRIDE_PURE(
            std::vector<float>, /* Return type */
            ActivityModel,      /* Parent class */
            get_a,          /* Name of function in C++ (must match Python name) */
            conc, T      /* Argument(s) */
        );
    }
};

void linkActivityModel(py::module& m)
{
    py::class_<ActivityModel, PyActivityModel>(m, "ActivityModel")
        .def(py::init<>())
        .def("get_y", &ActivityModel::get_y)
        .def("get_a", &ActivityModel::get_a);
}