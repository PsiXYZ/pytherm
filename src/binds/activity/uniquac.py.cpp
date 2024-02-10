#pragma once

#include "activitymodel.py.cpp"
#include <pytherm/activity/uniquac.cpp>


void linkUNIQUAC(py::module& m)
{
    py::class_<UNIQUAC, ActivityModel>(m, "UNIQUAC")
        .def(py::init<std::vector<float> &, std::vector<float> &, std::vector<std::vector<std::vector<float>>> &>(), "r"_a, "q"_a, "res_matrix"_a)
        
        .def("get_y", &UNIQUAC::get_y, "conc"_a, "T"_a, R"pbdoc(
        Calculate activity coefficients for conc array

        Concentrations must be in molar fractions

        .. math::
            \gamma_i =  \exp\left(\ln \gamma_i^c + \ln \gamma_i^r \right)

        Parameters
        ----------
        conc : list[float]
            Input concentration array, [molar fraction]
        T : float
            Temperature, [K]

        Returns
        -------
        list[float]
            Activity coefficients

        Examples
        --------
        >>> UNIFAC.get_y([0.5, 0.5], T=298)
        )pbdoc")
        .def("get_a", &UNIQUAC::get_a, "conc"_a, "T"_a, R"pbdoc(
        Calculate activities for conc array

        Concentrations must be in molar fractions

        Parameters
        ----------
        conc : list[float]
            Input concentration array, [molar fraction]
        T : float
            Temperature, [K]

        Returns
        -------
        list[float]
            Activities

        Examples
        --------
        >>> UNIQUAC.get_a([0.5, 0.5], T=298)
        )pbdoc")
        .doc() = R"pbdoc(
        Implementation of the UNIQUAC model

        Parameters
        ----------
        r : list[float]
            list with substances' r
        q : list[float]
            list with substances' q
        res_matrix
            interaction parameters matrix
        )pbdoc";
}