#pragma once

#include <string>
#include <map>
#include <vector>
#include <binds/pybind11.h>

#include "activitymodel.h"


class UNIFAC_FV: public ActivityModel
{
protected:
    int n_groups;
    int n_comps;
    vector<vector<float>> group_comp;
    vector<vector<float>> psi;
    vector<float> M;
    float T;
    vector<float> q_v;
    vector<float> r_v;
    vector<vector<float>> ln_gamma_pure;
    vector<string> comp_names;
    vector<string> groups_names;
    vector<int> sub_id_global;
    vector<float> Q_v;
    vector<float> R_v;
    vector<int> id_global;
    vector<vector<vector<float>>> res_params;
    void calculate_vdw();
    vector<float> get_lny_comb_classic(const vector<float> &conc);
    vector<float> get_lny_comb_modified(const vector<float> &conc);
    vector<float> (UNIFAC::*get_lny_comb)(const vector<float> &);
    vector<float> get_lny_res(const vector<float> &conc);
    void calculate_psi(float T);``
    void calculate_gamma_pure();
    vector<float> get_gamma_gr(const vector<float> &conc);

public:
    UNIFAC_FV(ParametersUNIFAC &parameters, SubstancesUNIFAC &substances);
    vector<float> get_y(const vector<float> &conc, float T) override;
    vector<float> get_a(const vector<float> &conc, float T) override;
};