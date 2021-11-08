#include <cmath>

#include "shock_capturing_entropy_viscosity_process.h"
#include "utilities/parallel_utilities.h"

namespace Kratos {

double ShockCapturingEntropyViscosityProcess::ComputeEntropy(
    const double Density,
    const double Pressure,
    const double Gamma)
{
    return Density / (Gamma - 1.0) * std::log(Pressure / std::pow(Density, Gamma));
}


double ShockCapturingEntropyViscosityProcess::ComputeH(const Element& rElement)
{
    double h_squared = 0.0;
    const auto& r_geometry = rElement.GetGeometry();

    // H is the shortest edge of the element
    for(unsigned int i=0; i<r_geometry.size(); ++i)
    {
        const unsigned j = (i + 1) % r_geometry.size();
        const auto edge = r_geometry[j] - r_geometry[i];
        h_squared = std::min(h_squared, inner_prod(edge, edge));
    }

    return std::sqrt(h_squared);
}

double ShockCapturingEntropyViscosityProcess::ComputeElementalEntropy()
{
    block_for_each(mrModelPart.Elements(), [](Element &rElement)
    {

    });
}


} // namespace Kratos 