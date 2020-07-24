// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:
//

#if !defined(KRATOS_ADJOINT_SOLID_ELEMENT_H_INCLUDED )
#define  KRATOS_ADJOINT_SOLID_ELEMENT_H_INCLUDED

// System includes

// External include

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "utilities/adjoint_extensions.h"
#include "custom_elements/total_lagrangian.h"

namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class AdjointSolidElement
 * @ingroup StructuralMechanicsApplication
 * @brief A template class for creating adjoint elements for solids.
 */
template <class TPrimalElement>
class AdjointSolidElement
    : public Element
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
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(AdjointSolidElement);

    ///@}
    ///@name Life Cycle
    ///@{

    AdjointSolidElement(IndexType NewId = 0);

    AdjointSolidElement(IndexType NewId, GeometryType::Pointer pGeometry);

    AdjointSolidElement(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties);

    ///@}
    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId,
                            NodesArrayType const& ThisNodes,
                            PropertiesType::Pointer pProperties) const override;

    Element::Pointer Create(IndexType NewId,
                            GeometryType::Pointer pGeom,
                            PropertiesType::Pointer pProperties) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override;

    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override;

    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                               ProcessInfo& rCurrentProcessInfo) override;

    void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                      ProcessInfo& rCurrentProcessInfo) override;

    void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) override;

    void GetValuesVector(Vector& rValues, int Step = 0) const override;

    void EquationIdVector(EquationIdVectorType& rResult,
                          const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& rCurrentProcessInfo) const override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    void CalculateSensitivityMatrix(const Variable<array_1d<double, 3>>& rDesignVariable,
                                    Matrix& rOutput,
                                    const ProcessInfo& rCurrentProcessInfo) override;

    ///@}

protected:

private:
    ///@name Member Variables
    ///@{

    TPrimalElement mPrimalElement;

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
        rSerializer.save("mPrimalElement", mPrimalElement);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
        rSerializer.load("mPrimalElement", mPrimalElement);
    }
    ///@}
};

KRATOS_API_EXTERN template class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) AdjointSolidElement<TotalLagrangian>;

///@}

} // namespace Kratos.
#endif // KRATOS_ADJOINT_SOLID_ELEMENT_H_INCLUDED  defined
