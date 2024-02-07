#pragma once

#include <vector>

#include "activitymodel.h"


std::vector<float> get_lny_SH(const std::vector<float> &conc, const std::vector<float> &r, const std::vector<float> &q);

class UNIQUAC: public ActivityModel
{
private:
    std::vector<float> r;
    std::vector<float> q;
    std::vector<std::vector<std::vector<float>>> res_matrix;
    float T;
    int n_comp;
    std::vector<std::vector<float>> t_matrix;

    std::vector<float> get_lny_res(const std::vector<float> &conc);
    std::vector<float> (*get_lny_comb)(const std::vector<float> &, const std::vector<float> &,const std::vector<float> &);
    void update_t_matrix(float T);
public:
    UNIQUAC(std::vector<float> &r, std::vector<float> &q, std::vector<std::vector<std::vector<float>>> &res_matrix);
    std::vector<float> get_y(const std::vector<float> &conc, float T) override;
    std::vector<float> get_a(const std::vector<float> &conc, float T) override;
};