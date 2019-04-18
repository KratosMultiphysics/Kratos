#include "custom_utilities/potential_flow_utilities.h"

namespace Kratos
{
namespace PotentialFlowUtilities
{
template <int Dim, int NumNodes>
array_1d<double, NumNodes> GetPotentialOnNormalElement(const Element& rElement)
{
    const int kutta = rElement.GetValue(KUTTA);
    array_1d<double, NumNodes> potentials;

    if (kutta == 0)
        for (unsigned int i = 0; i < NumNodes; i++)
            potentials[i] = rElement.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
    else
        for (unsigned int i = 0; i < NumNodes; i++)
            if (!rElement.GetGeometry()[i].GetValue(TRAILING_EDGE))
                potentials[i] =
                    rElement.GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
            else
                potentials[i] = rElement.GetGeometry()[i].FastGetSolutionStepValue(
                    AUXILIARY_VELOCITY_POTENTIAL);

    return potentials;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template instantiation

template array_1d<double, 3> GetPotentialOnNormalElement<2, 3>(const Element& element);
} // namespace PotentialFlow
} // namespace Kratos
