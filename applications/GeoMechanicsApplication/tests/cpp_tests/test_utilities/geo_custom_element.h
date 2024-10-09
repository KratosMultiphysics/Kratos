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

#include "includes/element.h"

namespace Kratos::Testing
{

class GeoMockElement : public Element
{
public:
    GeoMockElement() : Element() {};

    /// Constructor using Geometry
    GeoMockElement(IndexType NewId, GeometryType::Pointer pGeometry) : Element(NewId, pGeometry) {};

    ///@}
    ///@name Operators
    ///@{

    void SetLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo);

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void SetMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo);

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void SetDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo);

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    void SetRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo);

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(Element::DofsVectorType& rDofList, const ProcessInfo& rCurrentProcessInfo) const override;

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
