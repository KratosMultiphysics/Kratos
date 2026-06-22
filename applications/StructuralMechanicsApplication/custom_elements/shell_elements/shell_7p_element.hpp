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

    void CovariantBaseVectorsMidsurface(array_1d<Vector,3>& akovr, const Matrix& rShapeFunctionGradientValues, const Vector& rNshape, 
    const ConfigurationType& rConfiguration, const double& thickness) const;

    void DirectorDerivatives(array_1d<Vector,2>& a3kvp,
    const Matrix& rShapeFunctionGradientValues, const double& thickness) const;

    void CovariantMetric(Matrix& rMetric,const array_1d<Vector,3>& rBaseVectorCovariant) const;

    void ContraVariantBaseVectors(array_1d<Vector,3>& rBaseVectors,const Matrix& rContraVariantMetric,
    const array_1d<Vector,3> rCovariantBaseVectors) const;

    void ContravariantMetric(Matrix& rMetric,const Matrix& rCovariantMetric, double& detMetric_body) const;

    void JacobiDeterminante(double& DetJ, const array_1d<Vector,3>& akovr) const;

    void CovariantBaseVectorsShellBody(array_1d<Vector,3>& gkovr, const Matrix& rShapeFunctionGradientValues, 
    const Vector& rNshape, const ConfigurationType& rConfiguration, const double& Theta3, const double& thickness) const;

    void CalculateMaterialLaw(BoundedMatrix<double, 12, 12>& CL, const Matrix& gmkonr, const double& thickness,
    const ConstitutiveLawType& option, const double& Theta3, const double& fact) const;

    void CalculatelinearBOperator(Matrix& bop, const array_1d<Vector,3>& CovariantBaseVectors, const array_1d<Vector,2>& DirectorDerivatives, 
    const Matrix& ShapeFunctionGradientValues, const Vector& Nshape, const SizeType& number_of_nodes) const;

    void s8_ansqshapefunctions(array_1d<double,2>& frq, array_1d<double,2>& fsq,const double xi, const double eta) const;

    void BOperatorANSTransverseShearmodification(Matrix& Bop, const array_1d<double,2>& frq, const array_1d<double,2>& fsq,
    const array_1d<array_1d<Vector,3>,4>& akovr_ans,const array_1d<Vector,2>& a3kvp,const array_1d<Matrix,4>& DN_ans,
    const Matrix& N_ans, const SizeType& number_of_nodes) const;

    void BOperatorANSCurvatureThicknessModification(Matrix& Bop, const array_1d<array_1d<Vector,3>,4>& akovr_ct_ans, 
    const Matrix& N_ct_ans, const double r, const double s, const Vector& Np, const SizeType& number_of_nodes) const;

    void CalculateEASShapeFunctions(Matrix& M0_eas, const double r, const double s,
    const array_1d<SizeType,3>& eas_modes_per_kinematic_variable_set, const SizeType& num_eas_modes) const;

    void BasisTransformationEASShapeFunctions(Matrix& T, const Matrix& M0_eas, Matrix& M_eas,
    const array_1d<Vector,3>& akonr0_eas, const array_1d<Vector,3>& akovr, const double detJ0_surface, const double detJ_surface) const;

    friend class Serializer;

    // A private default constructor necessary for serialization
    Shell7pElement() = default;
};

} // namespace Kratos