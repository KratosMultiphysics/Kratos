#include "utilities/section_properties_utility.h"
#include <cmath>
#include <stdexcept>

namespace Kratos
{

std::unordered_map<std::string, double>
SectionPropertiesUtility::ComputeRectangularSection(const std::vector<double>& rParameters)
{
    if (rParameters.size() != 2) {
        throw std::runtime_error("RECT section requires [DIM1, DIM2]");
    }

    const double b = rParameters[0];
    const double h = rParameters[1];

    std::unordered_map<std::string, double> values;
    values["AREA"] = b * h;
    values["I33"]  = b * std::pow(h,3) / 12.0;
    values["I22"]  = h * std::pow(b,3) / 12.0;
    values["J"]   = b * h * (b*b + h*h) / 12.0;

    return values;
}

std::unordered_map<std::string, double>
SectionPropertiesUtility::ComputeISection(const std::vector<double>& rParameters)
{
    if (rParameters.size() != 6) {
        throw std::runtime_error("I section requires [DIM1, DIM2, DIM3, DIM4, DIM5, DIM6]");
    }

    const double h = rParameters[0];
    const double a = rParameters[1];
    const double b = rParameters[2];
    const double t_w = rParameters[3];
    const double t_a = rParameters[4];
    const double t_b = rParameters[5];
    const double h_w = h - (t_a+t_b);
    const double h_f = h - 0.5*(t_a+t_b);
    const double A = t_a*a + h_w*t_w + b*t_b;
    const double y_c = (1/2 *h_w*(h_w+t_a)*t_w + h_f*t_b*b)/A;
    const double y_s = (t_b * h_f * std::pow(b,3)) / (t_b * std::pow(b,3) + t_a * std::pow(a,3));
    //const double y_na = y_c - y_s;

    std::unordered_map<std::string, double> values;
    
    values["AREA"] = A;
    values["I22"]  = (b*std::pow(t_b,3) / 12) + (a*std::pow(t_a,3) / 12) + (t_w*std::pow(h_w,3) / 12) + std::pow((h_f - y_c),2)*b*t_b + std::pow(y_c,2)*a*t_a + std::pow((y_c - (0.5*(h_w+t_a))),2)*h_w*t_w;
    values["I33"]  = (std::pow(b,3)*t_b / 12) + (std::pow(a,3)*t_a / 12) + (std::pow(t_w,3)*h_w / 12);
    values["J"]   = 1/3 * (std::pow(t_b,3)*b + std::pow(t_a,3)*a + std::pow(t_w,3)*h_f);
    values["K1"] = t_w * h_w / A;
    values["K2"] = (5*(a*t_a + b*t_b)) / (6*A);

    return values;
}

} // namespace Kratos