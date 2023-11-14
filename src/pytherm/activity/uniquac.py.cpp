#pragma once

#include "activitymodel.py.cpp"
#include "uniquac.cpp"


void linkUNIQUAC(py::module& m)
{
    py::class_<UNIQUAC, ActivityModel>(m, "UNIQUAC")
        .def(py::init<std::vector<float> &, std::vector<float> &, std::vector<std::vector<std::vector<float>>> &>())
        .def("get_y", &UNIQUAC::get_y, "conc"_a, "T"_a);
}