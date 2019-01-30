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
 * finite differencing  (adjoint semi analytic approach). It is designed to be used in adjoint
 * sensitivity analysis
 */
class AdjointFiniteDifferencingBaseElement : public Element
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(AdjointFiniteDifferencingBaseElement);
    ///@}

    ///@name Classes
    ///@{
    ///@}

    ///@name Life Cycle
    ///@{
    AdjointFiniteDifferencingBaseElement() : Element()
    {}

    AdjointFiniteDifferencingBaseElement(Element::Pointer pPrimalElement);

    AdjointFiniteDifferencingBaseElement(Element::Pointer pPrimalElement, bool HasRotationDofs);

    ~AdjointFiniteDifferencingBaseElement() override;

    ///@}

    ///@}
    ///@name Informations
    ///@{

    ///@name Operations
    ///@{

    // Basic

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo) const override;


    IntegrationMethod GetIntegrationMethod() const override
    {
        return mpPrimalElement->GetIntegrationMethod();
    }

    void GetValuesVector(Vector& values, int Step = 0) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->Initialize(rCurrentProcessInfo);
    }

    void ResetConstitutiveLaw() override
    {
        mpPrimalElement->ResetConstitutiveLaw();
    }

    void CleanMemory() override
    {
        mpPrimalElement->CleanMemory();
    }

    void InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->InitializeSolutionStep(rCurrentProcessInfo);
    }

    void InitializeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->InitializeNonLinearIteration(rCurrentProcessInfo);
    }

    void FinalizeNonLinearIteration(const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->FinalizeNonLinearIteration(rCurrentProcessInfo);
    }

    void FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->FinalizeSolutionStep(rCurrentProcessInfo);
    }

    void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateLocalSystem(rLeftHandSideMatrix,
                                              rRightHandSideVector,
                                              rCurrentProcessInfo);
    }

    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                       const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateLeftHandSide(rLeftHandSideMatrix,
                                               rCurrentProcessInfo);
    }

    void CalculateLeftHandSide(std::vector< MatrixType >& rLeftHandSideMatrices,
					const std::vector< Variable< MatrixType > >& rLHSVariables,
					const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateLeftHandSide(rLeftHandSideMatrices,
                                               rLHSVariables,
                                               rCurrentProcessInfo);
    }

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateRightHandSide(rRightHandSideVector,
                                                rCurrentProcessInfo);
    }

    void CalculateRightHandSide(std::vector< VectorType >& rRightHandSideVectors,
					const std::vector< Variable< VectorType > >& rRHSVariables,
					const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateRightHandSide(rRightHandSideVectors,
                                                rRHSVariables,
                                                rCurrentProcessInfo);
    }

    void CalculateFirstDerivativesContributions(MatrixType& rLeftHandSideMatrix,
							VectorType& rRightHandSideVector,
							const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateFirstDerivativesContributions(rLeftHandSideMatrix,
							rRightHandSideVector,
							rCurrentProcessInfo);
    }

    void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
					        const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateFirstDerivativesLHS(rLeftHandSideMatrix,
					        rCurrentProcessInfo);
    }

    void CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector,
					      const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateFirstDerivativesRHS(rRightHandSideVector,
					        rCurrentProcessInfo);
    }

    void CalculateSecondDerivativesContributions(MatrixType& rLeftHandSideMatrix,
							 VectorType& rRightHandSideVector,
							 const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateSecondDerivativesContributions(rLeftHandSideMatrix,
							    rRightHandSideVector,
							    rCurrentProcessInfo);
    }

    void CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
					       const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateSecondDerivativesLHS(rLeftHandSideMatrix,
					        rCurrentProcessInfo);
    }

    void CalculateSecondDerivativesRHS(VectorType& rRightHandSideVector,
					       const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateSecondDerivativesRHS(rRightHandSideVector,
					        rCurrentProcessInfo);
    }

    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateMassMatrix(rMassMatrix,rCurrentProcessInfo);
    }

    void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateDampingMatrix(rDampingMatrix, rCurrentProcessInfo);
    }

    void AddExplicitContribution(const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->AddExplicitContribution(rCurrentProcessInfo);
    }

    void AddExplicitContribution(const VectorType& rRHSVector,
                                const Variable<VectorType>& rRHSVariable,
                                const Variable<double >& rDestinationVariable,
                                const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->AddExplicitContribution(rRHSVector,
                                    rRHSVariable,
                                    rDestinationVariable,
                                    rCurrentProcessInfo);
    }

    void AddExplicitContribution(const VectorType& rRHSVector,
                                const Variable<VectorType>& rRHSVariable,
                                const Variable<array_1d<double,3> >& rDestinationVariable,
                                const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->AddExplicitContribution(rRHSVector,
                                rRHSVariable,
                                rDestinationVariable,
                                rCurrentProcessInfo);
    }

    void AddExplicitContribution(const MatrixType& rLHSMatrix,
                                const Variable<MatrixType>& rLHSVariable,
                                const Variable<Matrix>& rDestinationVariable,
                                const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->AddExplicitContribution(rLHSMatrix,
                                rLHSVariable,
                                rDestinationVariable,
                                rCurrentProcessInfo);
    }

    void Calculate(const Variable<double >& rVariable,
			   double& Output,
			   const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "Calculate of the adjoint base element is called!" << std::endl;
    }

    void Calculate(const Variable< array_1d<double,3> >& rVariable,
			   array_1d<double,3>& Output,
			   const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "Calculate of the adjoint base element is called!" << std::endl;
    }

    void Calculate(const Variable<Vector >& rVariable,
			   Vector& Output,
			   const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "Calculate of the adjoint base element is called!" << std::endl;
    }

    void Calculate(const Variable<Matrix >& rVariable,
			   Matrix& Output,
			   const ProcessInfo& rCurrentProcessInfo) override;

    // Results calculation on integration points
    void CalculateOnIntegrationPoints(const Variable<bool>& rVariable,
					      std::vector<bool>& rOutput,
					      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "CalculateOnIntegrationPoints of the adjoint base condition is called!" << std::endl;
    }

    void CalculateOnIntegrationPoints(const Variable<double>& rVariable,
					      std::vector<double>& rOutput,
					      const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
					      std::vector< array_1d<double, 3 > >& rOutput,
					      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "CalculateOnIntegrationPoints of the adjoint base condition is called!" << std::endl;
    }

    void CalculateOnIntegrationPoints(const Variable<array_1d<double, 6 > >& rVariable,
					      std::vector< array_1d<double, 6 > >& rOutput,
					      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "CalculateOnIntegrationPoints of the adjoint base condition is called!" << std::endl;
    }

    void CalculateOnIntegrationPoints(const Variable<Vector >& rVariable,
					      std::vector< Vector >& rOutput,
					      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "CalculateOnIntegrationPoints of the adjoint base condition is called!" << std::endl;
    }

    void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable,
					      std::vector< Matrix >& rOutput,
					      const ProcessInfo& rCurrentProcessInfo) override
    {
        KRATOS_ERROR << "CalculateOnIntegrationPoints of the adjoint base condition is called!" << std::endl;
    }

    void GetValueOnIntegrationPoints(const Variable<double>& rVariable,
					     std::vector<double>& rValues,
					     const ProcessInfo& rCurrentProcessInfo) override
    {
        this->CalculateOnIntegrationPoints(rVariable, rValues, rCurrentProcessInfo);
    }

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    // Sensitivity functions

    /**
     * Calculates the pseudo-load contribution of the element w.r.t. all properties which are available at the element.
     * This is done by finite differencing of the RHS of the primal element when perturbing a property value.
     * This operation is thread-save, because the property pointer is exchanged by a local one before pertubation.
     */
    void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Calculates the pseudo-load of the design variable SHAPE (coordinates of nodes) contribution of the element.
     * This is done by finite differencing of the RHS of the primal element when perturbing a nodal coordinate.
     * This operation is currently NOT thread-save!
     */
    void CalculateSensitivityMatrix(const Variable<array_1d<double,3>>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Calculates the stress-displacement derivative of the given rStressVariable.
     * For the linear case this happens by calculating the stress using a unit-displacement for each dof.
     */
    virtual void CalculateStressDisplacementDerivative(const Variable<Vector>& rStressVariable,
                                    Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo);
    /**
     * Calculates the stress-design variable derivative of the given rStressVariable.
     * this is done by finite differencing of the Calculate function of the primal element
     */
    void CalculateStressDesignVariableDerivative(const Variable<double>& rDesignVariable, const Variable<Vector>& rStressVariable,
                                        Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo);
    /**
     * Calculates the stress-design variable derivative of the given rStressVariable.
     * this is done by finite differencing of the Calculate function of the primal element
     */
    void CalculateStressDesignVariableDerivative(const Variable<array_1d<double,3>>& rDesignVariable,
                                            const Variable<Vector>& rStressVariable,
                                             Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo);
    /**
     * Gets the pointer to the primal element.
     */
    Element::Pointer pGetPrimalElement()
    {
        return mpPrimalElement;
    }
    ///@}

    ///@name Public specialized Access - Temporary
    ///@{
    ///@}

protected:

    ///@name Protected Lyfe Cycle
    ///@{


    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    /**
     * pointer to the primal element
     */
    Element::Pointer mpPrimalElement;
    bool mHasRotationDofs = false;

    ///@}

private:

    ///@name Private Classes
    ///@{
    ///@}

    ///@name Private Operations
    ///@{

    /**
     * Get the perturbation size for a scalar variable
     */
    double GetPerturbationSize(const Variable<double>& rDesignVariable);

    /**
     * Get the perturbation size for a vector variable
     */
    double GetPerturbationSize(const Variable<array_1d<double,3>>& rDesignVariable);

    /**
     * Get the perturbation size modification factor for a scalar variable.
     * The computed factor reflects the current value of the property (design variable).
     * Note: This approach is only based on experience.
     * This can be overwritten by derived classes.
     */
    virtual double GetPerturbationSizeModificationFactor(const Variable<double>& rVariable);

    /**
     * Get the perturbation size modification factor for a vector variable.
     * The computed factor reflects the size of the element.
     * Note: This approach is only based on experience.
     * This can be overwritten by derived classes.
     */
    virtual double GetPerturbationSizeModificationFactor(const Variable<array_1d<double,3>>& rDesignVariable);

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
