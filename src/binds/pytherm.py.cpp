#pragma once

#include "pybind11.h"
#include "activity.py.cpp"

void linkActivity(py::module& m);

PYBIND11_MODULE(cpp, m)
{
    linkActivity(m);
}