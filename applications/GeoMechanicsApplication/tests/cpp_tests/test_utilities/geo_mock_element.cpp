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
#include "geo_mock_element.h"

namespace Kratos::Testing
{

void GeoMockElement::SetLeftHandSide(const MatrixType& rLeftHandSideMatrix)
{
    mStiffnessMatrix = rLeftHandSideMatrix;
}

void GeoMockElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    rLeftHandSideMatrix = mStiffnessMatrix;
    ++mCountCalculateLeftHandSideCalled;
}

void GeoMockElement::SetMassMatrix(const MatrixType& rMassMatrix) { mMassMatrix = rMassMatrix; }

void GeoMockElement::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    rMassMatrix = mMassMatrix;
    ++mCountCalculateMassMatrixCalled;
}

void GeoMockElement::SetDampingMatrix(const MatrixType& rDampingMatrix)
{
    mDampingMatrix = rDampingMatrix;
}

void GeoMockElement::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    rDampingMatrix = mDampingMatrix;
    ++mCalculateDampingMatrixCalled;
}

void GeoMockElement::SetRightHandSide(const VectorType& rRightHandSideVector)
{
    mRhs = rRightHandSideVector;
}

void GeoMockElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    rRightHandSideVector = mRhs;
    ++mCountCalculateRightHandSideCalled;
}

void GeoMockElement::EquationIdVector(Element::EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    rResult.resize(2);

    rResult[0] = GetGeometry()[0].GetDof(DISPLACEMENT_X).EquationId();
    rResult[1] = GetGeometry()[0].GetDof(DISPLACEMENT_Y).EquationId();
}

void GeoMockElement::GetDofList(Element::DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    rElementalDofList.resize(2);

    rElementalDofList[0] = GetGeometry()[0].pGetDof(DISPLACEMENT_X);
    rElementalDofList[1] = GetGeometry()[0].pGetDof(DISPLACEMENT_Y);
}

std::size_t GeoMockElement::GetCountCalculateLeftHandSideCalled() const
{
    return mCountCalculateLeftHandSideCalled;
}

std::size_t GeoMockElement::GetCountCalculateMassMatrixCalled() const
{
    return mCountCalculateMassMatrixCalled;
}

std::size_t GeoMockElement::GetCalculateDampingMatrixCalled() const
{
    return mCalculateDampingMatrixCalled;
}

std::size_t GeoMockElement::GetCountCalculateRightHandSideCalled() const
{
    return mCountCalculateRightHandSideCalled;
}
} // namespace Kratos::Testing
