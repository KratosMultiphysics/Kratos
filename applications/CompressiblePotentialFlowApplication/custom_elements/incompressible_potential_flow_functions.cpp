#include "custom_elements/incompressible_potential_flow_functions.h"

namespace Kratos
{

template <int Dim, int NumNodes>
void ComputeLHSGaussPointContribution(const double weight,
                                      Matrix& lhs,
                                      const ElementalData<NumNodes, Dim>& data)
{
    noalias(lhs) += weight * prod(data.DN_DX, trans(data.DN_DX));
}

template <int Dim, int NumNodes>
void GetPotential(const IncompressiblePotentialFlowElement<Dim, NumNodes>* element,
                  array_1d<double, NumNodes>& phis)
{
    for (unsigned int i = 0; i < NumNodes; i++)
        phis[i] = element->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
}

template <int Dim, int NumNodes>
void ComputeVelocity(const IncompressiblePotentialFlowElement<Dim, NumNodes>* element,
                     array_1d<double, Dim>& velocity)
{
    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(element->GetGeometry(), data.DN_DX, data.N, data.vol);

    GetPotential<2, 3>(element, data.phis);

    noalias(velocity) = prod(trans(data.DN_DX), data.phis);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template instantiation

template void ComputeLHSGaussPointContribution<2, 3>(const double weight,
                                                     Matrix& lhs,
                                                     const ElementalData<3, 2>& data);

template void GetPotential<2, 3>(const IncompressiblePotentialFlowElement<2, 3>* element,
                                 array_1d<double, 3>& phis);

template void ComputeVelocity<2, 3>(const IncompressiblePotentialFlowElement<2, 3>* element,
                                    array_1d<double, 2>& velocity);
} // namespace Kratos