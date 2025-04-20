#pragma once

#include "unifac.h"

class UNIFAC_VISCO: public UNIFAC
{
public:
    UNIFAC_VISCO(ParametersUNIFAC &parameters, SubstancesUNIFAC &substances);
    double get_GE_RT(const vector<float> &conc, double T);
};