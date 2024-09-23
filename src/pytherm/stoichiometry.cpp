#include "stoichiometry.h"

int getChargeFromString(std::string &s)
{
    auto const sign_i = s.find_last_of('+-');
    if (sign_i == s.npos){
        return 0;
    }
    else {
        return std::stoi(s.substr(sign_i + 1));
    }
}
