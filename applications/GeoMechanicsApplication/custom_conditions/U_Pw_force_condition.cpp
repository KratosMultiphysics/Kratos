// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Ignasi de Pouplana,
//                   Vahid Galavi
//

// Application includes
#include "custom_conditions/U_Pw_force_condition.hpp"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer UPwForceCondition<TDim, TNumNodes>::Create(IndexType             NewId,
                                                              NodesArrayType const& ThisNodes,
                                                              PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new UPwForceCondition(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
}

template <unsigned int TDim, unsigned int TNumNodes>
void UPwForceCondition<TDim, TNumNodes>::CalculateRHS(Vector& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    array_1d<double, 3> force_vector = this->GetGeometry()[0].FastGetSolutionStepValue(POINT_LOAD);
    std::copy_n(force_vector.begin(), TDim, rRightHandSideVector.begin());
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string UPwForceCondition<TDim, TNumNodes>::Info() const
{
    return "UPwForceCondition";
}

template class UPwForceCondition<2, 1>;
template class UPwForceCondition<3, 1>;

} // Namespace Kratos.
