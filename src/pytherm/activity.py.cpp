#pragma once

#include <pytherm/pybind11.h>

#include "activity/activitymodel.py.cpp"
#include "activity/unifac.py.cpp"
#include "activity/uniquac.py.cpp"


void linkActivityModel(py::module& m);
void linkUNIFAC(py::module& m);
void linkUNIQUAC(py::module& m);

PYBIND11_MODULE(_activity, m)
{
    linkActivityModel(m);
    linkUNIFAC(m);
    linkUNIQUAC(m);
}
