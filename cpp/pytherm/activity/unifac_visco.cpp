#include "unifac_visco.h"


UNIFAC_VISCO::UNIFAC_VISCO(ParametersUNIFAC &parameters, SubstancesUNIFAC &substances): UNIFAC(parameters, substances)
{
    
}

double UNIFAC_VISCO::get_GE_RT(const vector<float> &conc, double T)
{   
    if (this->T != T)
    {
        this->T = T;
        calculate_psi(T);
        calculate_gamma_pure();
    }

    auto lny_comb = (this->*UNIFAC::get_lny_comb)(conc);
    auto lny_res = this->UNIFAC::get_lny_res(conc);

    double GE_RT_comb = 0;
    double DE_RT_res = 0;

    for (int i = 0; i < lny_comb.size(); ++i)
    {
        GE_RT_comb += conc[i] * lny_comb[i];
        DE_RT_res += - conc[i] * lny_res[i];
    }

    return GE_RT_comb + DE_RT_res;
}
