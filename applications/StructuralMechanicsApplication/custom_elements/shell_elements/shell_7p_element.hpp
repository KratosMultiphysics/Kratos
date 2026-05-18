// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/element.h"

namespace Kratos
{

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) Shell7pElement : public Element
{
public:

    using SizeType = Element::SizeType;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(Shell7pElement);

    Shell7pElement(IndexType NewId, GeometryType::Pointer pGeometry);

    Shell7pElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    ~Shell7pElement() override = default;

    Element::Pointer Create(
        IndexType NewId,
        GeometryType::Pointer pGeom,
        PropertiesType::Pointer pProperties
    ) const override;

    Element::Pointer Create(
        IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties
    ) const override;

    void GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const override;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo
    ) const override;

    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
    ) override;
    
    void CalculateLeftHandSide(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

    void CalculateRightHandSide(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

private:

    enum class ConstitutiveLawType {
        gStVenantKirchhoff,
        oStVenantKirchhoff
        };

    enum class ConfigurationType {
        Current,
        Reference
    };

    void CovariantBaseVectors(array_1d<Vector,3>& rBaseVectors, const Matrix& rShapeFunctionGradientValues, const Vector& rNshape, 
    const ConfigurationType& rConfiguration, const double& thickness) const;

    void DirectorDerivatives(array_1d<Vector,2>& rDirectorDerivatives,const array_1d<Vector,3>& rBaseVectorCovariant,
    const Matrix& rShapeFunctionGradientValues, const double& thickness) const;

    void CovariantMetric(Matrix& rMetric,const array_1d<Vector,3>& rBaseVectorCovariant);

    void ContraVariantBaseVectors(array_1d<Vector,3>& rBaseVectors,const Matrix& rContraVariantMetric,
    const array_1d<Vector,3> rCovariantBaseVectors);

    void ContravariantMetric(Matrix& rMetric,const Matrix& rCovariantMetric);

    void JacobiDeterminante(double& rDetJacobi, const array_1d<Vector,3>& rReferenceBaseVectors) const;

    void CalculateMaterialLaw(BoundedMatrix<double, 12, 12>& CL, const Matrix& Gmkon, const double& thickness,
    const ConstitutiveLawType& option);

    void CalculatelinearBOperator(Matrix& bop, const array_1d<Vector,3>& CovariantBaseVectors, const array_1d<Vector,2>& DirectorDerivatives, 
    const Matrix& ShapeFunctionGradientValues, const Vector& Nshape, const SizeType& number_of_nodes);

    friend class Serializer;

    // A private default constructor necessary for serialization
    Shell7pElement() = default;
};

} // namespace Kratos