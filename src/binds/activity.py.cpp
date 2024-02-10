#pragma once

#include <binds/pybind11.h>

#include "activity/activitymodel.py.cpp"
#include "activity/unifac.py.cpp"
#include "activity/uniquac.py.cpp"


void linkActivityModel(py::module& m);
void linkUNIFAC(py::module& m);
void linkUNIQUAC(py::module& m);

void linkActivity(py::module& m)
{
    linkActivityModel(m);
    linkUNIFAC(m);
    linkUNIQUAC(m);
}
