// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Armin Geiser, https://github.com/armingeiser
//

#pragma once


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/element.h"
#include "structural_mechanics_application_variables.h"
#include "utilities/openmp_utils.h"
#include "custom_response_functions/adjoint_elements/adjoint_finite_difference_base_element.h"
#include "custom_elements/active_shell_3p_element.h"

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

/**
 * @brief Finite-difference adjoint wrapper for active-shell designs (ActiveAdjointFiniteDifferencingBaseElement)
 *
 * @details This element wraps a primal element of type `TPrimalElement` and computes design
 *          sensitivities via finite differences (semi-analytic adjoint approach). Most element
 *          operations are delegated to the wrapped primal element.
 *
 * @par Key extensions w.r.t. @ref AdjointFiniteDifferencingBaseElement
 * - Extends the DOF layout by appending adjoint actuation DOFs stored on the parent-geometry
 *   actuation node (`ACTIVE_SHELL_NODE_GP`) in `EquationIdVector`, `GetDofList`, and `GetValuesVector`.
 * - Adds finite-difference pseudo-load contributions w.r.t. the active shell actuation variables
 *   (`ACTIVE_SHELL_ALPHA/BETA/GAMMA/KAPPA_1/KAPPA_2/KAPPA_12`) in `CalculateSensitivityMatrix`.
 *
 * @tparam TPrimalElement Primal element type being wrapped (e.g. `ActiveShell3pElement`).
 *
 * @note This implementation is intended for Kirchhoff-Love shells and, in the active-shell use case,
 *       typically uses translational DOFs only.
 *            
 */

template <typename TPrimalElement>
class ActiveAdjointFiniteDifferencingBaseElement : public Element
{
public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>);
    ///@}

    ///@name Classes
    ///@{
    ///@}

    ///@name Life Cycle
    ///@{

    ActiveAdjointFiniteDifferencingBaseElement(IndexType NewId = 0,
                        bool HasRotationDofs = false)
    : Element(NewId),
      mpPrimalElement(Kratos::make_intrusive<TPrimalElement>(NewId, pGetGeometry())),
      mHasRotationDofs(HasRotationDofs)
    {
    }

    ActiveAdjointFiniteDifferencingBaseElement(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        bool HasRotationDofs = false)
    : Element(NewId, pGeometry),
      mpPrimalElement(Kratos::make_intrusive<TPrimalElement>(NewId, pGeometry)),
      mHasRotationDofs(HasRotationDofs)
    {
    }

    ActiveAdjointFiniteDifferencingBaseElement(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties,
                        bool HasRotationDofs = false)
    : Element(NewId, pGeometry, pProperties),
      mpPrimalElement(Kratos::make_intrusive<TPrimalElement>(NewId, pGeometry, pProperties)),
      mHasRotationDofs(HasRotationDofs)
    {       
            //CECKLEO was davon war schon da davor? lomplett neu nur klammern lassen
            KRATOS_WATCH(pGeometry);
            if (pGeometry) {
                KRATOS_WATCH(pGeometry->PointsNumber());
            } else {
                KRATOS_INFO("ActiveAdjointFiniteDifferencingBaseElement") << "pGeometry is nullptr!" << std::endl;
            }
            if (mpPrimalElement) {
                KRATOS_WATCH(mpPrimalElement->pGetGeometry());
                if (mpPrimalElement->pGetGeometry()) {
                    KRATOS_WATCH(mpPrimalElement->pGetGeometry()->PointsNumber());
                } else {
                    KRATOS_INFO("ActiveAdjointFiniteDifferencingBaseElement") << "mpPrimalElement->pGetGeometry() is nullptr!" << std::endl;
                }
            }
            KRATOS_WATCH(pProperties); //CHECKLEO -> Outbut untersuchung, weil Primal element keine Properies zum initialisieren hat 
            if (pProperties) {
                KRATOS_WATCH(pProperties->Id());
            } else {
                KRATOS_INFO("ActiveAdjointFiniteDifferencingBaseElement") << "pProperties is nullptr!" << std::endl;
            }
            if (mpPrimalElement) {
                KRATOS_WATCH(mpPrimalElement->pGetProperties());
                KRATOS_WATCH(mpPrimalElement->GetProperties().Id());
            } else {
                KRATOS_INFO("ActiveAdjointFiniteDifferencingBaseElement") << "mpPrimalElement is nullptr!" << std::endl;
            }
    }

    ///@}

    ///@}
    ///@name Information
    ///@{

    ///@name Operations
    ///@{

    Element::Pointer Create(IndexType NewId,
                              NodesArrayType const& ThisNodes,
                              PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
    }

    Element::Pointer Create(IndexType NewId,
                              GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties) const override
    {
        return Kratos::make_intrusive<ActiveAdjointFiniteDifferencingBaseElement<TPrimalElement>>(
            NewId, pGeometry, pProperties);
    }

    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const override;

    void GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& CurrentProcessInfo) const override;


    IntegrationMethod GetIntegrationMethod() const override
    {
        return mpPrimalElement->GetIntegrationMethod();
    }

    void GetValuesVector(Vector& values, int Step = 0) const override;

    void Initialize(const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->SetProperties(pGetProperties()); //CHECKLEO das ist neu -braucht man das?
        mpPrimalElement->Initialize(rCurrentProcessInfo);
    }

    void ResetConstitutiveLaw() override
    {
        mpPrimalElement->ResetConstitutiveLaw();
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

    void CalculateRightHandSide(VectorType& rRightHandSideVector,
                                        const ProcessInfo& rCurrentProcessInfo) override
    {
        mpPrimalElement->CalculateRightHandSide(rRightHandSideVector,
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
    void CalculateOnIntegrationPoints(
        const Variable<bool>& rVariable,
        std::vector<bool>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<double>& rVariable,
        std::vector<double>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 3>>& rVariable,
        std::vector<array_1d<double, 3>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 4>>& rVariable,
        std::vector<array_1d<double, 4>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 6>>& rVariable,
        std::vector<array_1d<double, 6>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<array_1d<double, 9>>& rVariable,
        std::vector<array_1d<double, 9>>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<Vector>& rVariable,
        std::vector<Vector>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    void CalculateOnIntegrationPoints(
        const Variable<Matrix>& rVariable,
        std::vector<Matrix>& rOutput,
        const ProcessInfo& rCurrentProcessInfo) override;

    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    // Sensitivity functions

    /**
     * Calculates the pseudo-load contribution of the element w.r.t. all properties which are available at the element.
     * This is done by finite differencing of the RHS of the primal element when perturbing a property value.
        * This operation is thread-safe, because the property pointer is exchanged by a local one before perturbation.
     */
    void CalculateSensitivityMatrix(const Variable<double>& rDesignVariable, Matrix& rOutput,
                                            const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * Calculates the pseudo-load of the design variable SHAPE_SENSITIVITY (coordinates of nodes) contribution of the element.
     * This is done by finite differencing of the RHS of the primal element when perturbing a nodal coordinate.
        * This operation is currently NOT thread-safe!
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
        * This is done by finite differencing of the Calculate function of the primal element.
     */
    void CalculateStressDesignVariableDerivative(const Variable<double>& rDesignVariable, const Variable<Vector>& rStressVariable,
                                        Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo);
    /**
     * Calculates the stress-design variable derivative of the given rStressVariable.
        * This is done by finite differencing of the Calculate function of the primal element.
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
    const Element::Pointer pGetPrimalElement() const
    {
        return mpPrimalElement;
    }
    ///@}

    ///@name Public specialized Access - Temporary
    ///@{
    ///@}

protected:
    ///@name Property/Material Checks (SMA-like) //CHECKLEO das ist neu
    void CheckProperties(const ProcessInfo& rCurrentProcessInfo) const; //CHECKLEO das ist neu

    ///@name Protected Lyfe Cycle
    ///@{


    ///@}
    ///@name Operations
    ///@{

    template <typename TDataType>
    void CalculateAdjointFieldOnIntegrationPoints(const Variable<TDataType>& rVariable, std::vector< TDataType >& rOutput, const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_WARNING_IF("CalculateAdjointFieldOnIntegrationPoints", OpenMPUtils::IsInParallel() != 0)
                << "The call of this non shared-memory-parallelized function within a parallel section should be avoided for efficiency reasons!" << std::endl;

        const SizeType num_nodes = mpPrimalElement->GetGeometry().PointsNumber();
        const SizeType dimension = mpPrimalElement->GetGeometry().WorkingSpaceDimension();
        const SizeType num_dofs_per_node = (mHasRotationDofs) ?  2 * dimension : dimension;
        const SizeType num_dofs = num_nodes * num_dofs_per_node;
        Vector initial_state_variables;
        initial_state_variables.resize(num_dofs);

        Vector particular_solution = ZeroVector(num_dofs);
        if(this->Has(ADJOINT_PARTICULAR_DISPLACEMENT)) {
            particular_solution = this->GetValue(ADJOINT_PARTICULAR_DISPLACEMENT);
        }

        // Build vector of variables containing the DOF-variables of the primal problem
        std::vector<Variable<double>*> primal_solution_variable_list;
        (mHasRotationDofs) ? primal_solution_variable_list = {&DISPLACEMENT_X, &DISPLACEMENT_Y, &DISPLACEMENT_Z, &ROTATION_X, &ROTATION_Y, &ROTATION_Z} :
                             primal_solution_variable_list = {&DISPLACEMENT_X, &DISPLACEMENT_Y, &DISPLACEMENT_Z};

        // Build vector of variables containing the DOF-variables of the adjoint problem
        std::vector<Variable<double>*> adjoint_solution_variable_list;
        (mHasRotationDofs) ? adjoint_solution_variable_list = {&ADJOINT_DISPLACEMENT_X, &ADJOINT_DISPLACEMENT_Y, &ADJOINT_DISPLACEMENT_Z, &ADJOINT_ROTATION_X, &ADJOINT_ROTATION_Y, &ADJOINT_ROTATION_Z} :
                             adjoint_solution_variable_list = {&ADJOINT_DISPLACEMENT_X, &ADJOINT_DISPLACEMENT_Y, &ADJOINT_DISPLACEMENT_Z};

        for (IndexType i = 0; i < num_nodes; ++i) {
            const IndexType index = i * num_dofs_per_node;
            for(IndexType j = 0; j < primal_solution_variable_list.size(); ++j) {
                initial_state_variables[index + j] = mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(*primal_solution_variable_list[j]);
                mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(*primal_solution_variable_list[j]) =
                                            this->GetGeometry()[i].FastGetSolutionStepValue(*adjoint_solution_variable_list[j]) +
                                            particular_solution[index + j];
            }
        }

        mpPrimalElement->CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);

        for (IndexType i = 0; i < num_nodes; ++i) {
            const IndexType index = i * num_dofs_per_node;
            for(IndexType j = 0; j < primal_solution_variable_list.size(); ++j) {
                mpPrimalElement->GetGeometry()[i].FastGetSolutionStepValue(*primal_solution_variable_list[j]) = initial_state_variables[index + j];
            }
        }

        KRATOS_CATCH("")
    }

    //CHECKLEO das ist neu
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info() << std::endl;
    }    

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "ActiveAdjointFiniteDifferenceElement #" << Id();
        return buffer.str();
    }

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
    double GetPerturbationSize(const Variable<double>& rDesignVariable, const ProcessInfo& rCurrentProcessInfo) const;

    /**
     * Get the perturbation size for a vector variable
     */
    double GetPerturbationSize(const Variable<array_1d<double,3>>& rDesignVariable, const ProcessInfo& rCurrentProcessInfo) const;

    /**
     * Get the perturbation size modification factor for a scalar variable.
     * The computed factor reflects the current value of the property (design variable).
     * Note: This approach is only based on experience.
     * This can be overwritten by derived classes.
     */
    virtual double GetPerturbationSizeModificationFactor(const Variable<double>& rVariable) const;

    /**
     * Get the perturbation size modification factor for a vector variable.
     * The computed factor reflects the size of the element.
     * Note: This approach is only based on experience.
     * This can be overwritten by derived classes.
     */
    virtual double GetPerturbationSizeModificationFactor(const Variable<array_1d<double,3>>& rDesignVariable) const;

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
