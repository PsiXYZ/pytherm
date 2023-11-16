#pragma once

#include "uniquac.h"

std::vector<float> get_lny_SH(const std::vector<float> &conc, const std::vector<float> &r, const std::vector<float> &q)
{
    int n_comps = conc.size();
    std::vector<float> ln_y_comb(n_comps, 0);

    std::vector<float> phi(n_comps, 0);
    std::vector<float> theta(n_comps, 0);

    float s1 = 0;
    float s2 = 0;
    float s3 = 0;
    for (int i = 0; i < n_comps; ++i)
    {
        s1 += r[i] * conc[i]; // phi
        s2 += q[i] * conc[i]; // theta
    }

    for (int i = 0; i < n_comps; ++i)
    {
        phi[i] = r[i] / s1;
        theta[i] = q[i] / s2;

        ln_y_comb[i] = 1 - phi[i] + log(phi[i]) - 5 * q[i] * (1 - phi[i] / theta[i] + log(phi[i] / theta[i]));
    }

    return ln_y_comb;
}

UNIQUAC::UNIQUAC(std::vector<float> &r, std::vector<float> &q, std::vector<std::vector<std::vector<float>>> &res_matrix)
{
    this->r = r;
    this->q = q;
    this->T = -1;
    this->n_comp = r.size();
    this->res_matrix = res_matrix;
    this->get_lny_comb = get_lny_SH;

    // std::vector<std::vector<std::vector<float>>> m;
    // for(int i = 0; i < this->n_comp; ++i)
    // {
    //     std::vector<std::vector<float>> b;
    //     for(int j = 0; j < this->n_comp; ++j)
    //     {
    //         b.push_back(
    //             std::vector<float> {0, 0}
    //         );
    //     }
    //     m.push_back(b);
    // }
    // this->res_matrix = m;

    std::vector<std::vector<float>> m;
    for (int i = 0; i < this->n_comp; ++i)
    {
        std::vector<float> b;
        for (int j = 0; j < this->n_comp; ++j)
        {
            b.push_back(0.0);
        }
        m.push_back(b);
    }
    this->t_matrix = m;
}

std::vector<float> UNIQUAC::get_y(const std::vector<float> &conc, float T)
{
    if (this->T != T)
    {
        this->T = T;
        update_t_matrix(T);
    }

    std::vector<float> y(this->n_comp, 0);

    std::vector<float> lny_comb = this->get_lny_comb(conc, this->r, this->q);
    std::vector<float> lny_res = get_lny_res(conc);

    for (int i = 0; i < this->n_comp; ++i)
    {
        y[i] = exp(lny_comb[i] + lny_res[i]);
    }

    return y;
}

std::vector<float> UNIQUAC::get_lny_res(const std::vector<float> &conc)
{
    std::vector<float> lny_res(this->n_comp);
    for (int i = 0; i < this->n_comp; ++i)
    {
        float acc1 = 0;
        float acc2 = 0;
        for (int j = 0; j < this->n_comp; ++j)
        {
            acc1 += this->q[j] * conc[j] * this->t_matrix[j][i];
            acc2 += this->q[j] * conc[j];
        }
        float s1 = log(acc1 / acc2);

        float s2 = 0;
        for (int j = 0; j < this->n_comp; ++j)
        {
            float acc1 = 0;
            for (int k = 0; k < this->n_comp; ++k)
            {
                acc1 += this->q[k] * conc[k] * this->t_matrix[k][j];
            }
            s2 += this->q[j] * conc[j] * this->t_matrix[i][j] / acc1;
        }

        lny_res[i] = this->q[i] * (1 - s1 - s2);
    }

    return lny_res;
}

void UNIQUAC::update_t_matrix(float T)
{
    std::vector<float> temps{1, 1 / T};
    for (int i = 0; i < this->n_comp; ++i)
    {
        for (int j = 0; j < this->n_comp; ++j)
        {
            float s = 0;
            for (int k = 0; k < 2; ++k)
            {
                s += - this->res_matrix[i][j][k] * temps[k];
            }
            this->t_matrix[i][j] = exp(s);
        }
    }
}
