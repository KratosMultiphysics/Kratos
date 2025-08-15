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
#pragma once

#include "includes/condition.h"

namespace Kratos::Testing
{

class GeoMockCondition : public Condition
{
public:
    GeoMockCondition() : Condition() {};

    /// Constructor using Geometry
    GeoMockCondition(IndexType NewId, GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry) {};

    ///@}
    ///@name Operators
    ///@{

    void SetLeftHandSide(const MatrixType& rLeftHandSideMatrix);

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void SetMassMatrix(const MatrixType& rMassMatrix);

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void SetDampingMatrix(const MatrixType& rDampingMatrix);

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void SetRightHandSide(const VectorType& rRightHandSideVector);

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(Condition::DofsVectorType& rDofList, const ProcessInfo& rCurrentProcessInfo) const override;

    std::size_t GetCountCalculateLeftHandSideCalled() const;

    std::size_t GetCountCalculateMassMatrixCalled() const;

    std::size_t GetCalculateDampingMatrixCalled() const;

    std::size_t GetCountCalculateRightHandSideCalled() const;

private:
    MatrixType mMassMatrix;
    MatrixType mDampingMatrix;
    MatrixType mStiffnessMatrix;
    VectorType mRhs;

    std::size_t mCountCalculateLeftHandSideCalled  = 0;
    std::size_t mCountCalculateMassMatrixCalled    = 0;
    std::size_t mCalculateDampingMatrixCalled      = 0;
    std::size_t mCountCalculateRightHandSideCalled = 0;
};

} // namespace Kratos::Testing
