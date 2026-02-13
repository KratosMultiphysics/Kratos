// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:  Alejandro Cornejo
//

#pragma once

// Project includes
#include "timoshenko_beam_element_2D2N.h"

namespace Kratos
{

/**
 * @class NonLinearTimoshenkoBeamElement2D2N
 * @brief Derived from LinearTimoshenkoBeamElement2D2N. Currently behaves identically.
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) NonLinearTimoshenkoBeamElement2D2N
    : public LinearTimoshenkoBeamElement2D2N
{
public:
    ///@name Type Definitions
    ///@{
    using BaseType = LinearTimoshenkoBeamElement2D2N;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(NonLinearTimoshenkoBeamElement2D2N);

    ///@}
    ///@name Life Cycle
    ///@{

    NonLinearTimoshenkoBeamElement2D2N() {}

    NonLinearTimoshenkoBeamElement2D2N(IndexType NewId, GeometryType::Pointer pGeometry)
        : BaseType(NewId, pGeometry)
    {}

    NonLinearTimoshenkoBeamElement2D2N(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : BaseType(NewId, pGeometry, pProperties)
    {}

    NonLinearTimoshenkoBeamElement2D2N(NonLinearTimoshenkoBeamElement2D2N const& rOther)
        : BaseType(rOther)
    {}

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<NonLinearTimoshenkoBeamElement2D2N>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<NonLinearTimoshenkoBeamElement2D2N>(NewId, pGeom, pProperties);
    }

    Element::Pointer Clone(IndexType NewId, NodesArrayType const& rThisNodes) const override;

    /**
     * @brief Override shape function methods to return global-sized vectors automatically
     */
    void GetShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) const override;
    void GetFirstDerivativesShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) const override;
    void GetSecondDerivativesShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) const override;
    void GetThirdDerivativesShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) const override;
    void GetFourthDerivativesShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) const override;

    void GetNThetaShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) const override;
    void GetFirstDerivativesNThetaShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) const override;

    void GetNu0ShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) const override;
    void GetFirstDerivativesNu0ShapeFunctionsValues(VectorType& rN, const double Length, const double Phi, const double xi) const override;

    ///@}
    ///@name Serialization
    ///@{
    // void save(Serializer &rSerializer) const override { BaseType::save(rSerializer); }
    // void load(Serializer &rSerializer) override { BaseType::load(rSerializer); }

private:
    // Nothing else for now - behaves like the linear element
};

} // namespace Kratos
