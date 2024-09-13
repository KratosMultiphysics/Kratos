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

class GeoCustomElement : public Element
{
public:
    /**
     * Constructor using Geometry
     */
    GeoCustomElement() : Element() {};
    /// Constructor using Geometry
    GeoCustomElement(IndexType NewId, GeometryType::Pointer pGeometry)
        : Element(NewId, pGeometry) {};

    GeoCustomElement(IndexType NewId, const NodesArrayType& ThisNodes)
        : Element(NewId, std::make_shared<GeometryType>(ThisNodes))
    {
    }

    ///@}
    ///@name Operators
    ///@{

    Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GeoCustomElement>(NewId, GetGeometry().Create(ThisNodes));
    }

    Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<GeoCustomElement>(NewId, pGeom);
    }

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(Element::DofsVectorType& rDofList, const ProcessInfo& rCurrentProcessInfo) const override;

private:
    MatrixType mMassMatrix;
    MatrixType mDampingMatrix;
    MatrixType mStiffnessMatrix;
    VectorType mRhs;

    bool mStiffnessMatrixSet = false;
    bool mMassMatrixSet      = false;
    bool mDampingMatrixSet   = false;
    bool mRhsSet             = false;
};

} // namespace Kratos::Testing
