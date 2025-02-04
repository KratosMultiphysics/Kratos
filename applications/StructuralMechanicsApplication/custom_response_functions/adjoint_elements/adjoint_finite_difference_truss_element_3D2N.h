// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//


#pragma once

#include "adjoint_finite_difference_base_element.h"

namespace Kratos
{

/** \brief AdjointFiniteDifferencingBaseElement
 *
 * This element is a wrapper for a primal truss element. It is responsible to deliver local stresses and
 * the stress displacement derivative. It is designed to be used in adjoint
 * sensitivity analysis.
 */
template <typename TPrimalElement>
class KRATOS_API(KRATOS_STRUCTURAL_MECHANICS_APPLICATION) AdjointFiniteDifferenceTrussElement
    : public AdjointFiniteDifferencingBaseElement<TPrimalElement>
{
public:

    // redefine the typedefs because of templated base class
    typedef AdjointFiniteDifferencingBaseElement<TPrimalElement> BaseType;
    typedef typename BaseType::SizeType SizeType;
    typedef typename BaseType::IndexType IndexType;
    typedef typename BaseType::GeometryType GeometryType;
    typedef typename BaseType::PropertiesType PropertiesType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::VectorType VectorType;
    typedef typename BaseType::MatrixType MatrixType;
    typedef typename BaseType::EquationIdVectorType EquationIdVectorType;
    typedef typename BaseType::DofsVectorType DofsVectorType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::IntegrationMethod IntegrationMethod;
    typedef typename BaseType::GeometryDataType GeometryDataType;

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AdjointFiniteDifferenceTrussElement);

    AdjointFiniteDifferenceTrussElement(IndexType NewId = 0)
    : BaseType(NewId)
    {
    }

    AdjointFiniteDifferenceTrussElement(IndexType NewId, typename GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry)
    {
    }

    AdjointFiniteDifferenceTrussElement(IndexType NewId,
                        typename GeometryType::Pointer pGeometry,
                        typename PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties)
    {
    }

    Element::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AdjointFiniteDifferenceTrussElement<TPrimalElement>>(
            NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(IndexType NewId,
                              typename GeometryType::Pointer pGeometry,
                              typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AdjointFiniteDifferenceTrussElement<TPrimalElement>>(
            NewId, pGeometry, pProperties);
    }

    /**
     * Calculates the derivative of stresses/stress resultants w.r.t primal displacement. The calculation is done analytically.
     * The derivative consists of two parts: The analytic derivative of the current length w.r.t. displacement
     * and an individual pre-factor for the different stresses/stress resultants.
     */
    void CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

private:
    void CheckDofs() const;

    void CheckProperties(const ProcessInfo& rCurrentProcessInfo) const;

    /**
     * Calculates the derivative of the current length w.r.t. primal displacements.
     */
    void CalculateCurrentLengthDisplacementDerivative(Vector& rDerivativeVector);

    /**
    * Calculates the stress displacement derivative pre-factor. This pre-factor gives together with the current length displacement
    * derivative the complete stress displacement derivative.
    */
    void GetDerivativePreFactor(double& rDerivativePreFactor, const ProcessInfo& rCurrentProcessInfo);

    /**
    * Calculates the individual stress displacement derivative pre-factor for the normal force.
    */
    double CalculateDerivativePreFactorFX(const ProcessInfo& rCurrentProcessInfo);

    /**
    * Calculates the individual stress displacement derivative pre-factor for the 2nd Piola-Kirchhoff stress.
    */
    double CalculateDerivativePreFactorPK2(const ProcessInfo& rCurrentProcessInfo);

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};


}
