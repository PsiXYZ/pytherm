#pragma once

#include "activitymodel.py.cpp"
#include "unifac.cpp"


using namespace pybind11::literals;

class PyUNIFAC : public UNIFAC
{
public:
    using UNIFAC::UNIFAC; // Inherit constructors
    vector<float> get_y(const vector<float> &conc, float T) override { PYBIND11_OVERRIDE_PURE(vector<float>, UNIFAC, get_y, conc, T); }
};

class PyUNIFAC_W : public UNIFAC_W
{
public:
    using UNIFAC_W::UNIFAC_W; // Inherit constructors
    vector<float> get_y(const vector<float> &conc, float T) override { PYBIND11_OVERRIDE(vector<float>, UNIFAC_W, get_y, conc, T); }
};

void linkUNIFAC(py::module& m)
{
    py::class_<ParametersUNIFAC>(m, "ParametersUNIFAC")
        .def(py::init<std::string>(), R"pbdoc(
            Load parameters form .dat file
            )pbdoc",
            "path"_a)
        
        .doc() = "A special class that holds UNIFAC parameters";

    py::class_<SubstancesUNIFAC>(m, "SubstancesUNIFAC")
        .def(py::init<>())

        .def("get_from_dict", &SubstancesUNIFAC::get_from_dict, R"pbdoc(
        Create :obj:`.SubstancesUNIFAC` from dict

        Parameters
        ----------
        subs_dict : dict
            Substances dict

        Examples
        --------
            >>> subs = {
            ...    "n-hexane": "2*CH3 4*CH2",
            ...    "butanone-2": "1*CH3 1*CH2 1*CH3CO",
            ... }
            >>> substances = SubstancesUNIFAC()
            >>> substances.get_from_dict(subs)
        )pbdoc",
             "subs_dict"_a)

        .doc() = "A special class that holds UNIFAC substances";

    py::class_<UNIFAC, ActivityModel, PyUNIFAC>(m, "UNIFAC")
        .def(py::init<ParametersUNIFAC &, SubstancesUNIFAC &>())

        .def("get_a", &UNIFAC::get_a, R"pbdoc(
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
        >>> UNIFAC.get_a([0.5, 0.5], T=298)
        )pbdoc",
             "conc"_a, "T"_a)

        .def("get_y", &UNIFAC::get_y, R"pbdoc(
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
        )pbdoc",
             "conc"_a, "T"_a)

        .doc() = R"pbdoc(
        Implementation of the UNIFAC model to work with molar fractions

        UNIFAC type (classic or modified) depends on :obj:`.ParametersUNIFAC`

        Parameters
        ----------
        dataset : ParametersUNIFAC
            ParametersUNIFAC object with interaction parameters
        substances : SubstancesUNIFAC
            Substances UNIFAC object with substance's group representation
        )pbdoc";

    py::class_<UNIFAC_W, UNIFAC, PyUNIFAC_W>(m, "UNIFAC_W")
        .def(py::init<ParametersUNIFAC &, SubstancesUNIFAC &, vector<float> &>())

        .def("get_a", &UNIFAC::get_a, R"pbdoc(
        Calculate activities for conc array

        Concentrations must be in weight fractions

        Parameters
        ----------
        conc : list[float]
            Input concentration array, [weight fractions]
        T : float
            Temperature, [K]

        Returns
        -------
        list[float]
            Activities

        Examples
        --------
        >>> UNIFAC.get_a([0.5, 0.5], T=298)
        )pbdoc",
             "conc"_a, "T"_a)

        .def("get_y", &UNIFAC_W::get_y, R"pbdoc(
        Calculate activity coefficients for conc array

        Concentrations must be in weight fractions

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
      )pbdoc",
             "conc"_a, "T"_a)
             
         .doc() = R"pbdoc(
        Implementation of the UNIFAC model to work with weight fractions

        UNIFAC type (classic or modified) depends on :obj:`.ParametersUNIFAC`

        Parameters
        ----------
        dataset : ParametersUNIFAC
            ParametersUNIFAC object with interaction parameters
        substances : SubstancesUNIFAC
            Substances UNIFAC object with substance's group representation
        )pbdoc";
}