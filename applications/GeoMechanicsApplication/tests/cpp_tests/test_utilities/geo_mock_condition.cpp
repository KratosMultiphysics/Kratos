// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Aron Noordam
//
#include "includes/variables.h"
#include "geo_mock_condition.h"

namespace Kratos::Testing
{

void GeoMockCondition::SetLeftHandSide(const MatrixType& rLeftHandSideMatrix)
{
    mStiffnessMatrix = rLeftHandSideMatrix;
}

void GeoMockCondition::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    rLeftHandSideMatrix = mStiffnessMatrix;
    ++mCountCalculateLeftHandSideCalled;
}

void GeoMockCondition::SetMassMatrix(const MatrixType& rMassMatrix) { mMassMatrix = rMassMatrix; }

void GeoMockCondition::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    rMassMatrix = mMassMatrix;
    ++mCountCalculateMassMatrixCalled;
}

void GeoMockCondition::SetDampingMatrix(const MatrixType& rDampingMatrix)
{
    mDampingMatrix = rDampingMatrix;
}

void GeoMockCondition::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    rDampingMatrix = mDampingMatrix;
    ++mCalculateDampingMatrixCalled;
}

void GeoMockCondition::SetRightHandSide(const VectorType& rRightHandSideVector)
{
    mRhs = rRightHandSideVector;
}

void GeoMockCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    rRightHandSideVector = mRhs;
    ++mCountCalculateRightHandSideCalled;
}

void GeoMockCondition::EquationIdVector(Condition::EquationIdVectorType& rResult,
                                        const ProcessInfo&               rCurrentProcessInfo) const
{
    rResult.resize(2);

    rResult[0] = GetGeometry()[0].GetDof(DISPLACEMENT_X).EquationId();
    rResult[1] = GetGeometry()[0].GetDof(DISPLACEMENT_Y).EquationId();
}

void GeoMockCondition::GetDofList(Condition::DofsVectorType& rElementalDofList,
                                  const ProcessInfo&         rCurrentProcessInfo) const
{
    rElementalDofList.resize(2);

    rElementalDofList[0] = GetGeometry()[0].pGetDof(DISPLACEMENT_X);
    rElementalDofList[1] = GetGeometry()[0].pGetDof(DISPLACEMENT_Y);
}

std::size_t GeoMockCondition::GetCountCalculateLeftHandSideCalled() const
{
    return mCountCalculateLeftHandSideCalled;
}

std::size_t GeoMockCondition::GetCountCalculateMassMatrixCalled() const
{
    return mCountCalculateMassMatrixCalled;
}

std::size_t GeoMockCondition::GetCalculateDampingMatrixCalled() const
{
    return mCalculateDampingMatrixCalled;
}

std::size_t GeoMockCondition::GetCountCalculateRightHandSideCalled() const
{
    return mCountCalculateRightHandSideCalled;
}
} // namespace Kratos::Testing
