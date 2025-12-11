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
//  Main authors:    Aron Noordam
//

// Application includes
#include "custom_conditions/Pw_condition.h"
#include "custom_utilities/dof_utilities.h"

namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes>
PwCondition<TDim, TNumNodes>::PwCondition() : PwCondition(0, nullptr, nullptr)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
PwCondition<TDim, TNumNodes>::PwCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : PwCondition(NewId, pGeometry, nullptr)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
PwCondition<TDim, TNumNodes>::PwCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer PwCondition<TDim, TNumNodes>::Create(IndexType               NewId,
                                                        NodesArrayType const&   ThisNodes,
                                                        PropertiesType::Pointer pProperties) const
{
    return Create(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::Pointer PwCondition<TDim, TNumNodes>::Create(IndexType NewId,
                                                        Geometry<GeometricalObject::NodeType>::Pointer pGeom,
                                                        Properties::Pointer pProperties) const
{
    return Condition::Pointer(new PwCondition(NewId, pGeom, pProperties));
}

template <unsigned int TDim, unsigned int TNumNodes>
void PwCondition<TDim, TNumNodes>::GetDofList(DofsVectorType& rConditionDofList, const ProcessInfo&) const
{
    rConditionDofList = GetDofs();
}

template <unsigned int TDim, unsigned int TNumNodes>
void PwCondition<TDim, TNumNodes>::CalculateLocalSystem(Matrix&            rLeftHandSideMatrix,
                                                        Vector&            rRightHandSideVector,
                                                        const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int conditionSize = TNumNodes;

    // Resetting the LHS
    if (rLeftHandSideMatrix.size1() != conditionSize)
        rLeftHandSideMatrix.resize(conditionSize, conditionSize, false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(conditionSize, conditionSize);

    // Resetting the RHS
    if (rRightHandSideVector.size() != conditionSize)
        rRightHandSideVector.resize(conditionSize, false);
    noalias(rRightHandSideVector) = ZeroVector(conditionSize);

    this->CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void PwCondition<TDim, TNumNodes>::CalculateLeftHandSide(Matrix& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "PwCondition::CalculateLeftHandSide is not implemented" << std::endl;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void PwCondition<TDim, TNumNodes>::CalculateRightHandSide(Vector&            rRightHandSideVector,
                                                          const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int conditionSize = TNumNodes;

    // Resetting the RHS
    if (rRightHandSideVector.size() != conditionSize)
        rRightHandSideVector.resize(conditionSize, false);
    noalias(rRightHandSideVector) = ZeroVector(conditionSize);

    this->CalculateRHS(rRightHandSideVector, rCurrentProcessInfo);

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
void PwCondition<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo&) const
{
    rResult = Geo::DofUtilities::ExtractEquationIdsFrom(GetDofs());
}

template <unsigned int TDim, unsigned int TNumNodes>
void PwCondition<TDim, TNumNodes>::CalculateAll(Matrix&            rLeftHandSideMatrix,
                                                Vector&            rRightHandSideVector,
                                                const ProcessInfo& CurrentProcessInfo)
{
    this->CalculateRHS(rRightHandSideVector, CurrentProcessInfo);
}

template <unsigned int TDim, unsigned int TNumNodes>
void PwCondition<TDim, TNumNodes>::CalculateRHS(Vector& rRightHandSideVector, const ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "calling the default CalculateRHS method for a particular condition ... "
                    "illegal operation!!"
                 << std::endl;

    KRATOS_CATCH("")
}

template <unsigned int TDim, unsigned int TNumNodes>
Condition::DofsVectorType PwCondition<TDim, TNumNodes>::GetDofs() const
{
    return Geo::DofUtilities::ExtractDofsFromNodes(GetGeometry(), WATER_PRESSURE);
}

template <unsigned int TDim, unsigned int TNumNodes>
void PwCondition<TDim, TNumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition)
}

template <unsigned int TDim, unsigned int TNumNodes>
void PwCondition<TDim, TNumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition)
}

template <unsigned int TDim, unsigned int TNumNodes>
std::string PwCondition<TDim, TNumNodes>::Info() const
{
    return "PwCondition";
}

template class PwCondition<2, 1>;
template class PwCondition<2, 2>;
template class PwCondition<2, 3>;
template class PwCondition<2, 4>;
template class PwCondition<2, 5>;
template class PwCondition<3, 1>;
template class PwCondition<3, 3>;
template class PwCondition<3, 4>;

} // Namespace Kratos.
