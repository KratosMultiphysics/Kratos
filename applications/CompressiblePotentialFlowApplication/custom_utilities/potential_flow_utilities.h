# include "custom_elements/incompressible_potential_flow_element.h"

namespace Kratos
{
namespace PotentialFlow
{


//class Element;

template <int Dim, int NumNodes>
array_1d<double, NumNodes> GetPotentialOnNormalElement(const Element& rElement);

} // namespace PotentialFlow
} // namespace Kratos