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

#include "utilities/adjoint_extensions.h"
#include "adjoint_finite_difference_base_element.h"

namespace Kratos
{

template <typename TPrimalElement>
class AdjointFiniteDifferenceCrBeamElement
    : public AdjointFiniteDifferencingBaseElement<TPrimalElement>
{
    class ThisExtensions : public AdjointExtensions
    {
        Element* mpElement;

    public:
        explicit ThisExtensions(Element* pElement);

        void GetFirstDerivativesVector(std::size_t NodeId,
                                       std::vector<IndirectScalar<double>>& rVector,
                                       std::size_t Step) override;

        void GetSecondDerivativesVector(std::size_t NodeId,
                                        std::vector<IndirectScalar<double>>& rVector,
                                        std::size_t Step) override;

        void GetAuxiliaryVector(std::size_t NodeId,
                                std::vector<IndirectScalar<double>>& rVector,
                                std::size_t Step) override;

        void GetFirstDerivativesVariables(std::vector<VariableData const*>& rVariables) const override;

        void GetSecondDerivativesVariables(std::vector<VariableData const*>& rVariables) const override;

        void GetAuxiliaryVariables(std::vector<VariableData const*>& rVariables) const override;
    };     
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

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AdjointFiniteDifferenceCrBeamElement);

    AdjointFiniteDifferenceCrBeamElement(IndexType NewId = 0)
    : BaseType(NewId, true)
    {
    }

    AdjointFiniteDifferenceCrBeamElement(IndexType NewId, typename GeometryType::Pointer pGeometry)
    : BaseType(NewId, pGeometry, true)
    {
    }

    AdjointFiniteDifferenceCrBeamElement(IndexType NewId,
                        typename GeometryType::Pointer pGeometry,
                        typename PropertiesType::Pointer pProperties)
    : BaseType(NewId, pGeometry, pProperties, true)
    {
    }

    Element::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AdjointFiniteDifferenceCrBeamElement<TPrimalElement>>(
            NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(IndexType NewId,
                              typename GeometryType::Pointer pGeometry,
                              typename PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<AdjointFiniteDifferenceCrBeamElement<TPrimalElement>>(
            NewId, pGeometry, pProperties);
    }

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                 std::vector< array_1d<double, 3 > >& rOutput,
                const ProcessInfo& rCurrentProcessInfo) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override
    {
        BaseType::Initialize(rCurrentProcessInfo);
        this->SetValue(ADJOINT_EXTENSIONS, Kratos::make_shared<ThisExtensions>(this));
    }
protected:


private:
    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
};


}
