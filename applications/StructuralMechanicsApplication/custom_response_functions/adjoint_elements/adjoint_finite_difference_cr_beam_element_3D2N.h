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


#if !defined(ADJOINT_FINITE_DIFFERENCE_CR_BEAM_ELEMENT_H_INCLUDED )
#define  ADJOINT_FINITE_DIFFERENCE_CR_BEAM_ELEMENT_H_INCLUDED

#include "adjoint_finite_difference_base_element.h"
#include "custom_elements/cr_beam_element_3D2N.hpp"

namespace Kratos
{

class AdjointFiniteDifferenceCrBeamElement : public AdjointFiniteDifferencingBaseElement
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(AdjointFiniteDifferenceCrBeamElement);

    typedef Element BaseType;
    typedef BaseType::GeometryType GeometryType;
    typedef BaseType::NodesArrayType NodesArrayType;
    typedef BaseType::PropertiesType PropertiesType;
    typedef BaseType::IndexType IndexType;
    typedef BaseType::SizeType SizeType;
    typedef BaseType::MatrixType MatrixType;
    typedef BaseType::VectorType VectorType;
    typedef BaseType::EquationIdVectorType EquationIdVectorType;
    typedef BaseType::DofsVectorType DofsVectorType;

    AdjointFiniteDifferenceCrBeamElement(IndexType NewId,
                            GeometryType::Pointer pGeometry);
    AdjointFiniteDifferenceCrBeamElement(IndexType NewId,
                            GeometryType::Pointer pGeometry,
                            PropertiesType::Pointer pProperties,
                            Element::Pointer pPrimalElement);

    ~AdjointFiniteDifferenceCrBeamElement() override;

    BaseType::Pointer Create(
        IndexType NewId,
        NodesArrayType const& rThisNodes,
        PropertiesType::Pointer pProperties,
        Element::Pointer pPrimalElement) const override;

    void EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo) override;

    void GetValuesVector(Vector& rValues, int Step = 0) override;

    void Calculate(const Variable<Vector >& rVariable,
                        Vector& rOutput,
                        const ProcessInfo& rCurrentProcessInfo) override;

    void Calculate(const Variable<Matrix >& rVariable,
                        Matrix& rOutput,
                        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
                        std::vector<double>& rOutput,
                        const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
                        std::vector<double>& rValues,
                        const ProcessInfo& rCurrentProcessInfo) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) override;

protected:
    AdjointFiniteDifferenceCrBeamElement(): AdjointFiniteDifferencingBaseElement()
    {
    }


private:
    double GetDisturbanceMeasureCorrectionFactor(const Variable<array_1d<double,3>>& rDesignVariable) override;

    /**
     * pointer to the primal element
     */
    CrBeamElement3D2N::Pointer mpPrimalBeamElement;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};


}

#endif
