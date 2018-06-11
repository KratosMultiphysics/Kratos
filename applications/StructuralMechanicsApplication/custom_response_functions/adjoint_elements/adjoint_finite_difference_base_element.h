// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//

#if !defined(ADJOINT_FINITE_DIFFERENCE_BASE_ELEMENT_H_INCLUDED )
#define  ADJOINT_FINITE_DIFFERENCE_BASE_ELEMENT_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/element.h"

namespace Kratos
{

///@name Kratos Globals
///@{
///@}

///@name Type Definitions
///@{
///@}

///@name  Enum's
///@{
///@}

///@name  Functions
///@{
///@}

///@name Kratos Classes
///@{

/** \brief AdjointFiniteDifferencingBaseElement
 *
 * This element a wrapper for a primal element to calculate element derivatives using
 * finite differncing. It is designed to be used in semi-analytic adjoint
 * sensitivity analysis
 */
class AdjointFiniteDifferencingBaseElement : public Element
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(AdjointFiniteDifferencingBaseElement);

     typedef Element::PropertiesType PropertiesType;

     typedef Element::DofsArrayType DofsArrayType;

    ///@}

    ///@name Classes
    ///@{
    ///@}

    ///@name Life Cycle
    ///@{

    AdjointFiniteDifferencingBaseElement(IndexType NewId,
                         GeometryType::Pointer pGeometry);

    AdjointFiniteDifferencingBaseElement(IndexType NewId,
                         GeometryType::Pointer pGeometry,
                         PropertiesType::Pointer pProperties,
                         Element::Pointer pPrimalElement);

    ~AdjointFiniteDifferencingBaseElement() override;

    ///@}

    ///@name Operations
    ///@{

    // Basic

    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties,
                Element::Pointer pPrimalElement) const;

    // TODO Element::Pointer Clone (IndexType NewId, NodesArrayType const& ThisNodes) const override;

    void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

    void GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo) override;

    void GetValuesVector(Vector& values, int Step = 0) override;

    void Initialize() override
    {
        mpPrimalElement->Initialize();
    }

    void ResetConstitutiveLaw() override
    {
        mpPrimalElement->ResetConstitutiveLaw();
    }

    void CleanMemory() override
    {
        mpPrimalElement->CleanMemory();
    }

    void InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->InitializeSolutionStep(rCurrentProcessInfo);
    }

    void InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->InitializeNonLinearIteration(rCurrentProcessInfo);
    }

    void FinalizeNonLinearIteration(ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->FinalizeNonLinearIteration(rCurrentProcessInfo);
    }

    void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->FinalizeSolutionStep(rCurrentProcessInfo);
    }

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateLocalSystem(rLeftHandSideMatrix,
                                              rRightHandSideVector,
                                              rCurrentProcessInfo);
    }

    void CalculateLocalSystem(std::vector< MatrixType >& rLeftHandSideMatrices,
                                      const std::vector< Variable< MatrixType > >& rLHSVariables,
                                      std::vector< VectorType >& rRightHandSideVectors,
                                      const std::vector< Variable< VectorType > >& rRHSVariables,
                                      ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateLocalSystem(rLeftHandSideMatrices,
                                              rLHSVariables,
                                              rRightHandSideVectors,
                                              rRHSVariables,
                                              rCurrentProcessInfo);
    }

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateLeftHandSide(rLeftHandSideMatrix,
                                               rCurrentProcessInfo);
    }

    void CalculateLeftHandSide(std::vector< MatrixType >& rLeftHandSideMatrices,
					const std::vector< Variable< MatrixType > >& rLHSVariables,
					ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateLeftHandSide(rLeftHandSideMatrices,
                                               rLHSVariables,
                                               rCurrentProcessInfo);
    }

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateRightHandSide(rRightHandSideVector,
                                                rCurrentProcessInfo);
    }

    void CalculateRightHandSide(std::vector< VectorType >& rRightHandSideVectors,
					const std::vector< Variable< VectorType > >& rRHSVariables,
					ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateRightHandSide(rRightHandSideVectors,
                                                rRHSVariables,
                                                rCurrentProcessInfo);
    }

    void CalculateFirstDerivativesContributions(MatrixType& rLeftHandSideMatrix,
							VectorType& rRightHandSideVector,
							ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateFirstDerivativesContributions(rLeftHandSideMatrix,
							rRightHandSideVector,
							rCurrentProcessInfo);
    }

    void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
					        ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateFirstDerivativesLHS(rLeftHandSideMatrix,
					        rCurrentProcessInfo);
    }

    void CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector,
					      ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateFirstDerivativesRHS(rRightHandSideVector,
					        rCurrentProcessInfo);
    }

    void CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix,
							 VectorType& rRightHandSideVector,
							 ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateSecondDerivativesContributions(rLeftHandSideMatrix,
							    rRightHandSideVector,
							    rCurrentProcessInfo);
    }

    void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
					       ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateSecondDerivativesLHS(rLeftHandSideMatrix,
					        rCurrentProcessInfo);
    }

    void CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector,
					       ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateSecondDerivativesRHS(rRightHandSideVector,
					        rCurrentProcessInfo);
    }

    void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateMassMatrix(rMassMatrix,rCurrentProcessInfo);
    }

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateDampingMatrix(rDampingMatrix, rCurrentProcessInfo);
    }

    void AddExplicitContribution(ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->AddExplicitContribution(rCurrentProcessInfo);
    }

    void AddExplicitContribution(const VectorType& rRHSVector,
                                const Variable<VectorType>& rRHSVariable,
                                Variable<double >& rDestinationVariable,
                                const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->AddExplicitContribution(rRHSVector,
                                    rRHSVariable,
                                    rDestinationVariable,
                                    rCurrentProcessInfo);
    }

    void AddExplicitContribution(const VectorType& rRHSVector,
                                const Variable<VectorType>& rRHSVariable,
                                Variable<array_1d<double,3> >& rDestinationVariable,
                                const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->AddExplicitContribution(rRHSVector,
                                rRHSVariable,
                                rDestinationVariable,
                                rCurrentProcessInfo);
    }

    void AddExplicitContribution(const MatrixType& rLHSMatrix,
                                const Variable<MatrixType>& rLHSVariable,
                                Variable<Matrix>& rDestinationVariable,
                                const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->AddExplicitContribution(rLHSMatrix,
                                rLHSVariable,
                                rDestinationVariable,
                                rCurrentProcessInfo);
    }

    void Calculate(const Variable<Vector >& rVariable, Vector& rOutput,
                           const ProcessInfo& rCurrentProcessInfo) override;

    void Calculate(const Variable<Matrix >& rVariable, Matrix& rOutput,
                           const ProcessInfo& rCurrentProcessInfo) override;

    // TODO evaluate if other Calculate functions are necessary

    // TODO add functions from element.h line 744 - 882
    // Results calculation on integration points

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput,
                                        const ProcessInfo& rCurrentProcessInfo) override;

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;


    int Check(const ProcessInfo& rCurrentProcessInfo) override;

    // Sensitivity functions

    void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateSensitivityMatrix(const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo);

    void CalculateStressDesignVariableDerivative(const Variable<double>& rDesignVariable, const Variable<Vector>& rStressVariable,
                                        Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo);

    void CalculateStressDesignVariableDerivative(const Variable<array_1d<double,3>>& rDesignVariable,
                                            const Variable<Vector>& rStressVariable,
                                             Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo);

    PropertiesType::Pointer pGetProperties()
    {
        return mpPrimalElement->pGetProperties();
    }

    const PropertiesType::Pointer pGetProperties() const
    {
        return mpPrimalElement->pGetProperties();
    }

    PropertiesType& GetProperties()
    {
        return mpPrimalElement->GetProperties();
    }

    PropertiesType const& GetProperties() const
    {
        return mpPrimalElement->GetProperties();
    }

    void SetProperties(PropertiesType::Pointer pProperties)
    {
        mpPrimalElement->SetProperties(pProperties);
    }

    // TODO add functions from element.h line 1050 - 1157

    ///@}

    ///@name Public specialized Access - Temporary
    ///@{
    ///@}

protected:

    ///@name Protected Lyfe Cycle
    ///@{

    /**
     * Protected empty constructor TODO needed?
     */
    AdjointFiniteDifferencingBaseElement() : Element()
    {
    }

    ///@}
    ///@name Operations
    ///@{

    virtual void AfterPerturbation(const Variable<double>& rDesignVariable,
                            const ProcessInfo& rCurrentProcessInfo)
    {
    }

    virtual void AfterPerturbation(const Variable<array_1d<double,3>>& rDesignVariable,
                            const ProcessInfo& rCurrentProcessInfo)
    {
    }

    ///@}
    ///@name Member Variables
    ///@{

    /**
     * pointer to the primal element
     */
    Element::Pointer mpPrimalElement;

    ///@}

private:

    ///@name Private Classes
    ///@{
    ///@}

    ///@name Private Operations
    ///@{

    virtual double GetDisturbanceMeasureCorrectionFactor(const Variable<double>& rVariable);

    virtual double GetDisturbanceMeasureCorrectionFactor(const Variable<array_1d<double,3>>& rDesignVariable);

    ///@}

    ///@name Static Member Variables
    ///@{
    ///@}

    ///@name Member Variables
    ///@{
    ///@}

    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

    ///@}

    ///@name Private  Access
    ///@{
    ///@}

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Un accessible methods
    ///@{
    ///@}

};

}
#endif // ADJOINT_FINITE_DIFFERENCE_BASE_ELEMENT_H_INCLUDED
