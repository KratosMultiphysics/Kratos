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

    void CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) override;

private:
    void CheckVariables();

    void CheckDofs();

    void CheckProperties(const ProcessInfo& rCurrentProcessInfo);

    double GetPerturbationSizeModificationFactor(const Variable<array_1d<double,3>>& rDesignVariable) override;

    double CalculateReferenceLength();

    double CalculateCurrentLength();

    void CalculateCurrentLengthDisplacementDerivative(Vector& rDerivativeVector);

    void GetDerivativePreFactor(double& rDerivativePreFactor, const ProcessInfo& rCurrentProcessInfo);

    double CalculateDerivativePreFactorFX(const ProcessInfo& rCurrentProcessInfo);

    double CalculateDerivativePreFactorPK2(const ProcessInfo& rCurrentProcessInfo);

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};


}

#endif
