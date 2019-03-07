# include "custom_elements/incompressible_potential_flow_element.h"

namespace Kratos
{
template <int Dim, int NumNodes>
class IncompressiblePotentialFlowElement;

template <int Dim, int NumNodes>
void ComputeLHSGaussPointContribution(const double weight,
                                      Matrix& lhs,
                                      const ElementalData<NumNodes, Dim>& data);

template <int Dim, int NumNodes>
void GetPotential(const IncompressiblePotentialFlowElement<Dim, NumNodes>* element,
                  array_1d<double, NumNodes>& phis);

template <int Dim, int NumNodes>
void ComputeVelocity(const IncompressiblePotentialFlowElement<Dim, NumNodes>* element,
                     array_1d<double, Dim>& velocity);
}