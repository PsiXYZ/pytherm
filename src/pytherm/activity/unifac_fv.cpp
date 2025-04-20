#include "unifac_fv.h"

void UNIFAC_FV::calculate_vdw()
{
    for (int i = 0; i < this->comp_names.size(); i++)
    {   
        float r = 0;
        float q = 0;
        for (int j = 0; j < this->groups_names.size(); j++)
        {
            r += this->R_v[j] * this->group_comp[i][j];
            q += this->Q_v[j] * this->group_comp[i][j];
        }
        r /= this->M[i];
        q /= this->M[i];

        this->r_v.push_back(r);
        this->q_v.push_back(q);
    }
}

vector<float> UNIFAC_FV::get_lny_comb_classic(const vector<float> &conc)
{
    vector<float> ln_y_comb(this->n_comps, 0);

    vector<float> phi(this->n_comps, 0);
    vector<float> theta(this->n_comps, 0);

    float s1 = 0;
    float s2 = 0;
    float s3 = 0;
    for(int i = 0; i < this->n_comps; ++i)
    {
        s1 += this->r_v[i] * conc[i]; // phi
        s2 += this->q_v[i] * conc[i]; // theta
    }

    for(int i = 0; i < this->n_comps; ++i)
    {
        phi[i] = this->r_v[i] / s1;
        theta[i] = this->q_v[i] / s2;

        ln_y_comb[i] = 1 - phi[i] + log(phi[i]) - 5 * this->q_v[i] * (1 - phi[i] / theta[i] + log(phi[i] / theta[i]));
        
    }

    return ln_y_comb;
}