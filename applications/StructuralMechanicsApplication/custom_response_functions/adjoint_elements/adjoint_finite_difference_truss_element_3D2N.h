// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder, https://github.com/MFusseder
//


#if !defined(ADJOINT_FINITE_DIFFERENCE_TRUSS_ELEMENT_H_INCLUDED )
#define  ADJOINT_FINITE_DIFFERENCE_TRUSS_ELEMENT_H_INCLUDED

#include "adjoint_finite_difference_base_element.h"

namespace Kratos
{

/** \brief AdjointFiniteDifferencingBaseElement
 *
 * This element is a wrapper for a primal truss element. It is responsible to deliver local stresses and
 * the stress displacement derivative. It is designed to be used in adjoint
 * sensitivity analysis.
 */
class AdjointFiniteDifferenceTrussElement : public AdjointFiniteDifferencingBaseElement
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(AdjointFiniteDifferenceTrussElement);

    AdjointFiniteDifferenceTrussElement(): AdjointFiniteDifferencingBaseElement()
    {
    }

    AdjointFiniteDifferenceTrussElement(Element::Pointer pPrimalElement);

    ~AdjointFiniteDifferenceTrussElement() override;

    void Calculate(const Variable<Vector >& rVariable,
                        Vector& rOutput,
                        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Calculates the derivative of stresses/stress resultants w.r.t primal displacement. The calculation is done analytically.
     * The derivative consists of two parts: The analytic derivative of the current length w.r.t. displacement
     * and an individual pre-factor for the different stresses/stress resultants.
     */
    void CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) override;

private:
    void CheckVariables();

    void CheckDofs();

    void CheckProperties(const ProcessInfo& rCurrentProcessInfo);

    double CalculateReferenceLength();

    double CalculateCurrentLength();

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

#endif
