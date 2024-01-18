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
#include "custom_elements/small_strain_U_Pw_diff_order_axisymmetric_element.hpp"
#include "element_strategies/axisymmetric_stress_state.h"

namespace Kratos
{

//----------------------------------------------------------------------------------------
Element::Pointer SmallStrainUPwDiffOrderAxisymmetricElement::Create(IndexType NewId,
                                                                    NodesArrayType const& ThisNodes,
                                                                    PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SmallStrainUPwDiffOrderAxisymmetricElement(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------
Element::Pointer SmallStrainUPwDiffOrderAxisymmetricElement::Create(IndexType NewId,
                                                                    GeometryType::Pointer pGeom,
                                                                    PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new SmallStrainUPwDiffOrderAxisymmetricElement(NewId, pGeom, pProperties));
}

//----------------------------------------------------------------------------------------------------

} // Namespace Kratos
