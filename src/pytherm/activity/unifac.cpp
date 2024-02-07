#pragma once

#include <pytherm/pybind11.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

#include "unifac.h"

namespace py = pybind11;

void tokenize(std::string const &str, const char delim, std::vector<std::string> &out)
{
    std::stringstream ss(str);
 
    std::string s;
    while (std::getline(ss, s, delim)) {
        out.push_back(s);
    }
}

void removeSpaceAtStart(std::string &str)
{
    str.erase(0, 2);
}

void tokenizeBySpace(std::string const &str, std::vector<std::string> &out)
{   
    std::string v;
    for (char c:str)
    {
        if ((c != ' '))
        {
            v += c;
        }
        else
        {
            if (!v.empty()) out.push_back(v);
            v.clear();
        }
    }
    if (!v.empty()) out.push_back(v);
}

bool contains(std::vector<std::string> &v, std::string &key)
{
    for (auto s : v)
    {
        if (s == key) return true;
    }
    return false;
}


ParametersUNIFAC::ParametersUNIFAC(std::string path)
{   
    std::ifstream in(path);
    readData(in);
    in.close();

}

void ParametersUNIFAC::readData(std::ifstream &in)
{
    std::string line;
    if (in.is_open())
    {   
        std::getline(in, line);
        if (line == "UNIFAC_parameters")
        {
            // std::cout << "NORMES" << std::endl;
        }
        else
        {
            // std::cout << "WRONG FILE TYPE" << std::endl;
        }

        while (std::getline(in, line))
        {
            if (line == "-type")
            {
                readType(in);
            }
            if (line == "-list_of_main_groups")
            {
                readMainGroups(in);
            }
            if (line == "-list_of_sub_groups")
            {
                readSubGroups(in);
            }
            if (line == "-list_of_interaction_parameters")
            {
                readInter(in);
            }
        }
    }
}

void ParametersUNIFAC::readType(std::ifstream &in)
{   
    std::string line;
    std::getline(in, line);
    removeSpaceAtStart(line);
    if (line == "modified") this->unifacType = true;
    // std::cout << line << std::endl;
}

void ParametersUNIFAC::readMainGroups(std::ifstream &in)
{   
    std::vector<std::string> mainGroups;
    std::string line;
    while (std::getline(in, line))
    {   
        removeSpaceAtStart(line);
        if (line =="") break;

        std::vector<std::string> tokens;
        tokenizeBySpace(line, tokens);
        mainGroups.push_back(tokens[1]);
    }
    this->mainGroups = mainGroups;
}

void ParametersUNIFAC::readSubGroups(std::ifstream &in)
{
    std::string line;
    while (std::getline(in, line))
    {   
        removeSpaceAtStart(line);
        if (line =="") break;

        std::vector<std::string> tokens;
        tokenizeBySpace(line, tokens);

        this->subGroups.push_back(tokens[1]);
        this->subToMain.push_back(stoi(tokens[2]));
        this->R.push_back(stof(tokens[4]));
        this->Q.push_back(stof(tokens[5]));
        // std::cout << tokens[1] << std::endl;
    }
   
}

void ParametersUNIFAC::readInter(std::ifstream &in)
{
    std::string line;
    while (std::getline(in, line))
    {
        removeSpaceAtStart(line);
        if (line =="") break;

        std::vector<std::string> tokens;
        tokenizeBySpace(line, tokens);

        int i = stoi(tokens[0]);
        int j = stoi(tokens[1]);

        std::vector<float> v1 {stof(tokens[2]), stof(tokens[3]), stof(tokens[4]),};
        this->resParams[i][j] = v1;
        std::vector<float> v2 {stof(tokens[5]), stof(tokens[6]), stof(tokens[7]),};
        this->resParams[j][i] = v2;
    }
}

int ParametersUNIFAC::get_sub_id(std::string gr_name)
{
    for (int i = 0; i <= subGroups.size(); i++)
    {
        if (subGroups[i] == gr_name) return i;
    }
}


SubstancesUNIFAC::SubstancesUNIFAC()
{

}

void SubstancesUNIFAC::get_from_dict(py::dict s_dict)
{
    for (const auto& element : s_dict)
    {   
        std::vector<std::string> tokens;
        std::string gr = element.second.cast<std::string>();
        std::string name = element.first.cast<std::string>();

        tokenize(gr, ' ', tokens);
        for (auto el : tokens)
        {   
            std::vector<std::string> tokens2;
            tokenize(el, '*', tokens2);
            this->subs[name][tokens2[1]] = stof(tokens2[0]);
        }
        this->subs_names.push_back(name);
    }    
}


UNIFAC::UNIFAC(ParametersUNIFAC &parameters, SubstancesUNIFAC &substances): ActivityModel()
{   
    this->comp_names = substances.subs_names;

    for (auto s : substances.subs_names)
    {   
        for (auto el : substances.subs[s])
        {
            std::string gr_s = el.first;
            if (!contains(this->groups_names, gr_s))
            {
                this->groups_names.push_back(el.first);
                this->sub_id_global.push_back(parameters.get_sub_id(gr_s));
            }
        }
    }

    for (int i = 0; i < groups_names.size(); i++)
    {
        this->R_v.push_back(parameters.R[this->sub_id_global[i]]);
        this->Q_v.push_back(parameters.Q[this->sub_id_global[i]]);
    }

    for (auto i : this->sub_id_global)
    {
        this->id_global.push_back(parameters.subToMain[i]);
    }

    
    for (int i = 0; i < comp_names.size(); i++)
    {
        vector<float> b;
        for (int j = 0; j < groups_names.size(); j++)
        {
            auto comp = this->comp_names[i];
            auto gr = this->groups_names[j];
            auto n = substances.subs[comp][gr];
            b.push_back(n);
        }
        this->group_comp.push_back(b);
    }

    for(int i = 0; i < groups_names.size(); i++)
    {   
        vector<vector<float>> b;
        for(int j = 0; j < groups_names.size(); j++)
        {
            int i_global = this->id_global[i];
            int j_global = this->id_global[j];
            vector<float> a_ij {
                    parameters.resParams[i_global][j_global][0], 
                    parameters.resParams[i_global][j_global][1],  
                    parameters.resParams[i_global][j_global][2]
            };
            b.push_back(a_ij);
        }
        this->res_params.push_back(b);
    }

    for(int i = 0; i < groups_names.size(); i++)
    {   
        vector<float> b;
        for(int j = 0; j < groups_names.size(); j++)
        {
            b.push_back(0.0);
        }
        this->psi.push_back(b);
    }  

    for (int comp_i = 0; comp_i < this->comp_names.size(); comp_i++)
    {   
        vector<float> b;
        for (int gr_i = 0; gr_i < this->groups_names.size(); gr_i++)
        {
            b.push_back(0.0);
        }
        this->ln_gamma_pure.push_back(b);
    }

    this->T = -1;
    this->n_comps = this->comp_names.size();
    this->n_groups = this->groups_names.size();
    calculate_vdw();

    if (parameters.unifacType) 
    {
        this->get_lny_comb = &UNIFAC::get_lny_comb_modified;
    }
    else
    {
        this->get_lny_comb = &UNIFAC::get_lny_comb_classic;
    }
}

void UNIFAC::calculate_vdw()
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
        this->r_v.push_back(r);
        this->q_v.push_back(q);
    }
}

vector<float> UNIFAC::get_y(const vector<float> &conc, float T)
{   
    if (this->T != T)
    {
        this->T = T;
        calculate_psi(T);
        calculate_gamma_pure();
    }

    vector<float> y(this->n_comps, 0);

    vector<float> lny_comb = (this->*get_lny_comb)(conc);
    vector<float> lny_res = get_lny_res(conc);

    for(int i = 0; i < this->n_comps; ++i)
    {
        y[i] = exp(lny_comb[i] + lny_res[i]);
    }

    return y;
}

vector<float> UNIFAC::get_a(const vector<float> &conc, float T)
{   
    vector<float> y = get_y(conc, T);
    vector<float> a(conc.size());
    for(int i = 0; i < this->n_comps; ++i)
    {
        a[i] = y[i] * conc[i];
    }

    return a;
}

vector<float> UNIFAC::get_lny_comb_modified(const vector<float> &conc)
{
    vector<float> ln_y_comb(this->n_comps, 0);

    vector<float> phi_m(this->n_comps, 0);
    vector<float> phi(this->n_comps, 0);
    vector<float> theta(this->n_comps, 0);

    float s1 = 0;
    float s2 = 0;
    float s3 = 0;
    for(int i = 0; i < this->n_comps; ++i)
    {
        s1 += pow(this->r_v[i], 3.0 / 4.0) * conc[i]; // phi_m
        s2 += this->r_v[i] * conc[i]; // phi
        s3 += this->q_v[i] * conc[i]; // theta
    }

    for(int i = 0; i < this->n_comps; ++i)
    {
        phi_m[i] = pow(this->r_v[i], 3.0 / 4.0) / s1;
        phi[i] = this->r_v[i] / s2;
        theta[i] = this->q_v[i] / s3;

        ln_y_comb[i] = 1 - phi_m[i] + log(phi_m[i]) - 5 * this->q_v[i] * (1 - phi[i] / theta[i] + log(phi[i] / theta[i]));
    }

    return ln_y_comb;
}

vector<float> UNIFAC::get_lny_comb_classic(const vector<float> &conc)
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

vector<float> UNIFAC::get_lny_res(const vector<float> &conc)
{
    vector<float> ln_y_res(this->n_comps, 0);
    vector<float> ln_gamma_gr = get_gamma_gr(conc);

    for(int i = 0; i < this->n_comps; ++i)
    {   
        float s = 0;
        for(int k = 0; k < this->n_groups; ++k)
        {   
            s += this->group_comp[i][k] * (ln_gamma_gr[k] - this->ln_gamma_pure[i][k]);
        }
        ln_y_res[i] = s;
    }

    return ln_y_res;
}

void UNIFAC::calculate_psi(float T)
{
    vector<double> temps {1, T, pow(T,  2)};

    for(int i = 0; i < this->groups_names.size(); i++)
    {
        for(int j = 0; j < this->groups_names.size(); j++)
        {   
            float a = 0;
            for (int k = 0; k < res_params[i][j].size(); k++)
            {
                a += res_params[i][j][k] * temps[k];
            }
            this->psi[i][j] = exp(- a / T);
        }
    }

}

vector<float> UNIFAC::get_gamma_gr(const vector<float> &conc)
{
    vector<float> gamma_gr(this->n_groups, 0);

    vector<float> X(this->n_groups, 0.0);
    vector<float> THETA(this->n_groups, 0.0);

    float total_X = 0;
    for (int i = 0; i < this->n_comps; ++i)
    {
        for (int j = 0; j < this->n_groups; ++j)
        {
            total_X += this->group_comp[i][j] * conc[i];
        }
    }
    for (int i = 0; i < this->n_groups; ++i)
    {   
        float s = 0;
        for (int j = 0; j < this->n_comps; ++j)
        {
            s += group_comp[j][i] * conc[j];
        }
        X[i] = s;
    }

    float total_THETA = 0;
    for (int i = 0; i < this->n_groups; ++i)
    {
        total_THETA += this->Q_v[i] * X[i];
    }
    for (int i = 0; i < this->n_groups; ++i)
    {
        THETA[i] = (this->Q_v[i] * X[i]) / total_THETA;
    }

    for (int k = 0; k < this->n_groups; ++k)
    {
        float s1 = 0;
        for (int m = 0; m < this->n_groups; ++m)
        {
            s1 += THETA[m] * this->psi[m][k];
        }

        float s2 = 0;
        for (int m = 0; m < this->n_groups; ++m)
        {
            float s3 = 0;
            for (int n = 0; n < this->n_groups; ++n)
            {
                s3 += THETA[n] * this->psi[n][m];
            }
            s2 += (THETA[m] * psi[k][m]) / s3;
        }
        gamma_gr[k] = this->Q_v[k] * (1 - log(s1) - s2);
    } 

    return gamma_gr;
}

void UNIFAC::calculate_gamma_pure()
{
    
    for (int i = 0; i < this->comp_names.size(); i++)
    {
        vector<float> conc(this->comp_names.size(), 0);
        conc[i] = 1;
        vector<float> ln_gamma_pure = get_gamma_gr(conc);
        for (int j = 0; j < this->groups_names.size(); j++)
        {
            this->ln_gamma_pure[i][j] = ln_gamma_pure[j];
        }
    }

}


UNIFAC_W::UNIFAC_W(ParametersUNIFAC &parameters, SubstancesUNIFAC &substances, vector<float> &Mw): UNIFAC(parameters, substances)
{
    this->Mw = Mw;
}

vector<float> UNIFAC_W::get_y(const vector<float> &conc, float T)
{   
    float w_M = 0;
    for(int i = 0; i < this->n_comps; ++i)
    {
        w_M += conc[i] / this->Mw[i];
    }

    vector<float> conc_x(this->n_comps, 0);
    for(int i = 0; i < this->n_comps; ++i)
    {
        conc_x[i] = (conc[i] / this->Mw[i]) / w_M;
    }

    vector<float> y_x(this->n_comps, 0); 
    y_x = this->UNIFAC::get_y(conc_x, T);

    vector<float> y_w(this->n_comps, 0);
    for(int i = 0; i < this->n_comps; ++i)
    {
        y_w[i] = y_x[i] / (this->Mw[i] * w_M);
    }

    return y_w;
}
