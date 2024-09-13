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
#include "geo_custom_element.h"

namespace Kratos::Testing
{

void GeoCustomElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    if (!mStiffnessMatrixSet) {
        mStiffnessMatrix.resize(rLeftHandSideMatrix.size1(), rLeftHandSideMatrix.size2(), false);
        for (unsigned int i = 0; i < rLeftHandSideMatrix.size1(); i++) {
            for (unsigned int j = 0; j < rLeftHandSideMatrix.size2(); j++) {
                mStiffnessMatrix(i, j) = rLeftHandSideMatrix(i, j);
            }
        }

        mStiffnessMatrixSet = true;

    }

    else {
        rLeftHandSideMatrix.resize(mStiffnessMatrix.size1(), mStiffnessMatrix.size2(), false);
        for (unsigned int i = 0; i < mStiffnessMatrix.size1(); i++) {
            for (unsigned int j = 0; j < mStiffnessMatrix.size2(); j++) {
                rLeftHandSideMatrix(i, j) = mStiffnessMatrix(i, j);
            }
        }
    }
}

void GeoCustomElement::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    if (!mMassMatrixSet) {
        mMassMatrix.resize(rMassMatrix.size1(), rMassMatrix.size2(), false);
        for (unsigned int i = 0; i < rMassMatrix.size1(); i++) {
            for (unsigned int j = 0; j < rMassMatrix.size2(); j++) {
                mMassMatrix(i, j) = rMassMatrix(i, j);
            }
        }

        mMassMatrixSet = true;
    }

    else {
        rMassMatrix.resize(mMassMatrix.size1(), mMassMatrix.size2(), false);
        for (unsigned int i = 0; i < mMassMatrix.size1(); i++) {
            for (unsigned int j = 0; j < mMassMatrix.size2(); j++) {
                rMassMatrix(i, j) = mMassMatrix(i, j);
            }
        }
    }
}

void GeoCustomElement::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    if (!mDampingMatrixSet) {
        mDampingMatrix    = rDampingMatrix;
        mDampingMatrixSet = true;
    } else rDampingMatrix = mDampingMatrix;
}

void GeoCustomElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    if (!mRhsSet) {
        mRhs    = rRightHandSideVector;
        mRhsSet = true;
    } else rRightHandSideVector = mRhs;
}

void GeoCustomElement::EquationIdVector(Element::EquationIdVectorType& rResult,
                                        const ProcessInfo&             rCurrentProcessInfo) const
{
    rResult.resize(2);

    rResult[0] = GetGeometry()[0].GetDof(DISPLACEMENT_X).EquationId();
    rResult[1] = GetGeometry()[0].GetDof(DISPLACEMENT_Y).EquationId();
}

void GeoCustomElement::GetDofList(Element::DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const
{
    rElementalDofList.resize(2);

    rElementalDofList[0] = GetGeometry()[0].pGetDof(DISPLACEMENT_X);
    rElementalDofList[1] = GetGeometry()[0].pGetDof(DISPLACEMENT_Y);
}

} // namespace Kratos::Testing
