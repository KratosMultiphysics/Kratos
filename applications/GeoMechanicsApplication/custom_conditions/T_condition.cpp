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
//  Main authors:    Mohamed Nabi
//                   John van Esch
//

#include "custom_conditions/T_condition.h"
#include "custom_utilities/dof_utilities.h"

namespace Kratos {

template <unsigned int TDim, unsigned int TNumNodes>
GeoTCondition<TDim, TNumNodes>::GeoTCondition() : Condition()
{
}

template <unsigned int TDim, unsigned int TNumNodes>
GeoTCondition<TDim, TNumNodes>::GeoTCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
GeoTCondition<TDim, TNumNodes>::GeoTCondition(IndexType NewId,
                                              GeometryType::Pointer pGeometry,
                                              PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
GeoTCondition<TDim, TNumNodes>::~GeoTCondition() = default;

template <unsigned int TDim, unsigned int TNumNodes>
void GeoTCondition<TDim, TNumNodes>::GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo&) const
{
    rConditionDofList = GetDofs();
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoTCondition<TDim, TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                                          VectorType& rRightHandSideVector,
                                                          const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rLeftHandSideMatrix = ZeroMatrix(TNumNodes, TNumNodes);
    rRightHandSideVector = ZeroVector(TNumNodes);

    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoTCondition<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoTCondition<TDim, TNumNodes>::CalculateAll(MatrixType& rLeftHandSideMatrix,
                                                  VectorType& rRightHandSideVector,
                                                  const ProcessInfo& rCurrentProcessInfo)
{
    CalculateRHS(rRightHandSideVector, rCurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
void GeoTCondition<TDim, TNumNodes>::CalculateRHS(VectorType& rRightHandSideVector,
                                                  const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateRHS method for a particular "
                    "condition ... illegal operation!!"
                 << std::endl;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::DofsVectorType GeoTCondition<TDim, TNumNodes>::GetDofs() const
{
    return Geo::DofUtilities::ExtractDofsFromNodes(this->GetGeometry(), TEMPERATURE);
}

template class GeoTCondition<2, 2>;
template class GeoTCondition<2, 3>;
template class GeoTCondition<2, 4>;
template class GeoTCondition<2, 5>;
template class GeoTCondition<3, 3>;
template class GeoTCondition<3, 4>;
template class GeoTCondition<3, 6>;
template class GeoTCondition<3, 8>;
template class GeoTCondition<3, 9>;

} // namespace Kratos
