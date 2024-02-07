#pragma once

#include <vector>

class ActivityModel
{   
public:
    virtual std::vector<float> get_y(const std::vector<float> &conc, float T) = 0;
    virtual std::vector<float> get_a(const std::vector<float> &conc, float T) = 0;
    // ActivityModel() {}
};
