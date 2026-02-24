#pragma once

#include <vector>
#include <unordered_map>
#include <string>
#include "includes/define.h"

namespace Kratos
{

class SectionPropertiesUtility
{
public:

    static std::unordered_map<std::string, double>
    ComputeRectangularSection(const std::vector<double>& rParameters);
    //ComputeRectangularSection(const std::vector<double>& rParameters);

    static std::unordered_map<std::string, double>
    ComputeISection(const std::vector<double>& rParameters);
};

} // namespace Kratos