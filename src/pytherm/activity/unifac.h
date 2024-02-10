#pragma once

#include <string>
#include <map>
#include <vector>
#include <binds/pybind11.h>

#include "activitymodel.h"

using std::string;
using std::vector;

namespace py = pybind11;

class ParametersUNIFAC
{
private:
    void readData(std::ifstream &in);
    void readType(std::ifstream &in);
    void readMainGroups(std::ifstream &in);
    void readSubGroups(std::ifstream &in);
    void readInter(std::ifstream &in);
public:
    bool unifacType;
    std::vector<std::string> mainGroups;
    std::vector<std::string> subGroups;
    std::vector<int> subToMain;
    std::vector<float> R;
    std::vector<float> Q;
    std::map<int, std::map<int, std::vector<float>>> resParams;
    ParametersUNIFAC(std::string path);
    int get_sub_id(std::string gr_name);
    int get_R(int id);
    int get_Q(int id);
};

class SubstancesUNIFAC
{
private:

public:
    SubstancesUNIFAC();
    std::vector<std::string> subs_names;
    std::map<std::string, std::map<std::string, float>> subs;
    void get_from_dict(py::dict s_dict);
};


class UNIFAC: public ActivityModel
{
protected:
    int n_groups;
    int n_comps;
    vector<vector<float>> group_comp;
    vector<vector<float>> psi;
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
    void calculate_psi(float T);
    void calculate_gamma_pure();
    vector<float> get_gamma_gr(const vector<float> &conc);

public:
    UNIFAC(ParametersUNIFAC &parameters, SubstancesUNIFAC &substances);
    vector<float> get_y(const vector<float> &conc, float T) override;
    vector<float> get_a(const vector<float> &conc, float T) override;
};

class UNIFAC_W: public UNIFAC
{
protected:
    vector<float> Mw;
public:
    UNIFAC_W(ParametersUNIFAC &parameters, SubstancesUNIFAC &substances, vector<float> &Mw);
    vector<float> get_y(const vector<float> &conc, float T) override;
};