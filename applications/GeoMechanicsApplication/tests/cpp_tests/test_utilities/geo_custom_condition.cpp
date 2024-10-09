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
#include "custom_utilities/dof_utilities.h"
#include "geo_custom_condition.h"

namespace Kratos::Testing
{

    void GeoMockCondition::SetLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
    {
        mStiffnessMatrix = rLeftHandSideMatrix;
    }

    void GeoMockCondition::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
    {
        rLeftHandSideMatrix = mStiffnessMatrix;
        ++mCountCalculateLeftHandSideCalled;
    }

    void GeoMockCondition::SetMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
    {
        mMassMatrix = rMassMatrix;
    }

    void GeoMockCondition::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
    {
        rMassMatrix = mMassMatrix;
        ++mCountCalculateMassMatrixCalled;
    }

    void GeoMockCondition::SetDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
    {
        mDampingMatrix = rDampingMatrix;
    }

    void GeoMockCondition::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
    {
        rDampingMatrix = mDampingMatrix;
        ++mCalculateDampingMatrixCalled;
    }

    void GeoMockCondition::SetRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
    {
        mRhs = rRightHandSideVector;
    }

    void GeoMockCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
    {
        rRightHandSideVector = mRhs;
        ++mCountCalculateRightHandSideCalled;
    }

void GeoMockCondition::EquationIdVector(Condition::EquationIdVectorType& rResult,
                                          const ProcessInfo& rCurrentProcessInfo) const
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
