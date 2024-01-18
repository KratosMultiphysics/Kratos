// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

// Application includes
#include "custom_elements/U_Pw_small_strain_axisymmetric_FIC_element.hpp"

#include "element_strategies/axisymmetric_stress_state.h"

namespace Kratos
{

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwSmallStrainAxisymmetricFICElement<TDim, TNumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwSmallStrainAxisymmetricFICElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------
template <unsigned int TDim, unsigned int TNumNodes>
Element::Pointer UPwSmallStrainAxisymmetricFICElement<TDim, TNumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new UPwSmallStrainAxisymmetricFICElement(NewId, pGeom, pProperties));
}

//----------------------------------------------------------------------------------------------------

template class UPwSmallStrainAxisymmetricFICElement<2, 3>;
template class UPwSmallStrainAxisymmetricFICElement<2, 4>;

} // Namespace Kratos
