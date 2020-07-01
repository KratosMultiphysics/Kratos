//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Altug Emiroglu, https://github.com/emiroglu
//


#if !defined(KRATOS_MODAL_DERIVATIVE_SCHEME )
#define  KRATOS_MODAL_DERIVATIVE_SCHEME

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/dof.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/model_part.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "solving_strategies/schemes/scheme.h"

/* Application includes */
#include "rom_application_variables.h"
#include "custom_utilities/rom_finite_difference_utility.h"

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
 * @class ModalDerivativeScheme
 * @ingroup KratosCore
 * @brief This class provides the implementation of the basic tasks that are needed by the solution strategy.
 * @details It is intended to be the place for tailoring the solution strategies to problem specific tasks.
 * @tparam TSparseSpace The sparse space considered
 * @tparam TDenseSpace The dense space considered
 * @author Altug Emiroglu
 */
template<class TSparseSpace,
         class TDenseSpace //= DenseSpace<double>
         >
class ModalDerivativeScheme : public Scheme<TSparseSpace,TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ModalDerivativeScheme
    KRATOS_CLASS_POINTER_DEFINITION(ModalDerivativeScheme);

    // Base type definition
    typedef Scheme<TSparseSpace,TDenseSpace> BaseType;

    /// Data type definition
    typedef typename BaseType::TDataType TDataType;
    /// Matrix type definition
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    /// Vector type definition
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    /// Local system matrix type definition
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    /// Local system vector type definition
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    /// DoF type definition
    typedef Dof<double> TDofType;
    /// DoF array type definition
    typedef typename BaseType::DofsArrayType DofsArrayType;
    /// DoF iterator type definition
    typedef typename BaseType::DofIterator DofIterator;
    /// DoF constant iterator type definition
    typedef typename BaseType::DofConstantIterator DofConstantIterator;

    /// Elements containers definition
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    /// Conditions containers definition
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default Constructor
     * @details Initiliazes the flags
     */
    /// Constructor.
    ModalDerivativeScheme(double FiniteDifferenceStepSize, bool FiniteDifferenceTypeFlag=false)
    : 
    Scheme<TSparseSpace,TDenseSpace>(), 
    mFiniteDifferenceStepSize(FiniteDifferenceStepSize), 
    mFiniteDifferenceTypeFlag(FiniteDifferenceTypeFlag)
    {
    }

    /// Destructor.
    ~ModalDerivativeScheme() override 
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief unction to be called when it is needed to initialize an iteration. It is designed to be called at the beginning of each non linear iteration
     * @note Take care: the elemental function with the same name is NOT called here.
     * @warning Must be defined in derived classes
     * @details The function is called in the builder for memory efficiency
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    void InitializeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        KRATOS_TRY

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Definition of the first element iterator
        const auto it_elem_begin = rModelPart.ElementsBegin();

        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); ++i) {
            auto it_elem = it_elem_begin + i;
            it_elem->InitializeNonLinearIteration(r_current_process_info);
        }

        // Definition of the first condition iterator
        const auto it_cond_begin = rModelPart.ConditionsBegin();

        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); ++i) {
            auto it_cond = it_cond_begin + i;
            it_cond->InitializeNonLinearIteration(r_current_process_info);
        }

        // Definition of the first constraint iterator
        const auto it_const_begin = rModelPart.MasterSlaveConstraintsBegin();

        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.MasterSlaveConstraints().size()); ++i) {
            auto it_const = it_const_begin + i;
            it_const->InitializeNonLinearIteration(r_current_process_info);
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief It initializes a non-linear iteration (for an individual condition)
     * @warning Must be defined in derived classes
     * @param pCurrentElement The element to compute
     * @param rCurrentProcessInfo The current process info instance
     */
    void InitializeNonLinearIteration(
        Element::Pointer pCurrentElement,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        (pCurrentElement)->InitializeNonLinearIteration(rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    /**
     * @brief It initializes a non-linear iteration (for an individual condition)
     * @warning Must be defined in derived classes
     * @param rCurrentCondition The condition to compute
     * @param rCurrentProcessInfo The current process info instance
     */
    void InitializeNonLinearIteration(
        Condition::Pointer pCurrentCondition,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        (pCurrentCondition)->InitializeNonLinearIteration(rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    // Element contributions
    void CalculateSystemContributions(
        Element& rElement,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        this->CalculateLHSContribution(rElement, rLHS_Contribution, rEquationId, rCurrentProcessInfo);

        this->CalculateRHSContribution(rElement, rRHS_Contribution, rEquationId, rCurrentProcessInfo);

        rElement.EquationIdVector(rEquationId,rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateLHSContribution(
        Element& rElement,
        LocalSystemMatrixType& rLHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        if (rCurrentProcessInfo[BUILD_LEVEL] == 1 ) // Mass matrix
            rElement.CalculateMassMatrix(rLHS_Contribution, rCurrentProcessInfo);
        else if (rCurrentProcessInfo[BUILD_LEVEL] == 2) // Stiffness matrix
            rElement.CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);
        else 
            KRATOS_ERROR <<"Invalid BUILD_LEVEL: " << rCurrentProcessInfo[BUILD_LEVEL] << std::endl;
        
        rElement.EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void CalculateRHSContribution(
        Element& rElement,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY
        
        // Derivative of basis_i wrt basis_j
        std::size_t basis_i = rCurrentProcessInfo[BASIS_I];
        std::size_t basis_j = rCurrentProcessInfo[BASIS_J];
        
        // Create PhiElemental
        Vector PhiElemental;
        std::size_t num_element_nodal_dofs = 0;
        for (auto& node_i : rElement.GetGeometry())
            num_element_nodal_dofs += node_i.GetDofs().size();
        PhiElemental.resize(num_element_nodal_dofs);
        
        // Get PhiElemental
        std::size_t dof_ctr = 0;
        for (auto& node_i : rElement.GetGeometry()) {
            auto& node_i_dofs = node_i.GetDofs();
            
            const Matrix *pPhiNodal = &(node_i.GetValue(ROM_BASIS));
            for (std::size_t dof_idx = 0; dof_idx < node_i_dofs.size(); dof_idx++){
                PhiElemental[dof_ctr + dof_idx] = (*pPhiNodal)(dof_idx, basis_i);
            }
            dof_ctr += node_i_dofs.size();
        }
        
        // Build RHS contribution
        rRHS_Contribution.clear();
        std::vector<Dof<double>::Pointer> rElementalDofList;
        rElement.GetDofList(rElementalDofList, rCurrentProcessInfo);
        rRHS_Contribution.resize(rElementalDofList.size());
        for (std::size_t iRHS = 0; iRHS < rRHS_Contribution.size(); iRHS++)
            rRHS_Contribution[iRHS] = 0.0;

        // Compute element LHS derivative
        Matrix element_LHS_derivative;
        element_LHS_derivative.resize(rElementalDofList.size(),rElementalDofList.size(),false);

        // Positive perturbation
        Matrix LHS_p_perturbed;
        LHS_p_perturbed.resize(rElementalDofList.size(),rElementalDofList.size(),false);
        this->PerturbElement(rElement, 1.0, basis_j, rCurrentProcessInfo);
        rElement.CalculateLeftHandSide(LHS_p_perturbed, rCurrentProcessInfo);

        if (mFiniteDifferenceTypeFlag)
        {   // Central Difference
            
            // Negative perturbation
            Matrix LHS_n_perturbed;
            LHS_n_perturbed.resize(rElementalDofList.size(),rElementalDofList.size(),false);
            this->PerturbElement(rElement, -2.0, basis_j, rCurrentProcessInfo);
            rElement.CalculateLeftHandSide(LHS_n_perturbed, rCurrentProcessInfo);

            // Reset perturbation
            this->PerturbElement(rElement, 1.0, basis_j, rCurrentProcessInfo);

            // Compute LHS derivative
            element_LHS_derivative = (LHS_p_perturbed - LHS_n_perturbed) / (2.0*mFiniteDifferenceStepSize);
        }
        else
        {   // Forward difference

            // Reset perturbation
            this->PerturbElement(rElement, -1.0, basis_j, rCurrentProcessInfo);

            // Neutral state
            Matrix LHS;
            LHS.resize(rElementalDofList.size(),rElementalDofList.size(),false);
            rElement.CalculateLeftHandSide(LHS, rCurrentProcessInfo);

            // Compute LHS derivative
            element_LHS_derivative = (LHS_p_perturbed - LHS) / mFiniteDifferenceStepSize;

        }

        // Compute RHS contribution
        if (element_LHS_derivative.size1() != PhiElemental.size()){
            // TODO: this is a workaround. use a map as in ROM analysis
            // retrieve only displacement dofs from PhiElemental
            unsigned int disp_shifter = 3;
            Vector tmpPhiElemental;
            tmpPhiElemental.resize(PhiElemental.size() / 2);
            for (std::size_t iNode = 0; iNode < rElement.GetGeometry().size(); iNode++)
            {
                for (std::size_t iXYZ = 0; iXYZ < 3; iXYZ++)
                {
                    tmpPhiElemental[iNode * 3 + iXYZ] = PhiElemental[iNode * 6 + disp_shifter + iXYZ];
                }
            }
            rRHS_Contribution += prec_prod(element_LHS_derivative, tmpPhiElemental);
        }
        else
        {
            rRHS_Contribution += prec_prod(element_LHS_derivative, PhiElemental);
        }

        // Negate the RHS
        rRHS_Contribution *= -1.0;

        rElement.EquationIdVector(EquationId,rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    // This function perturbs the element with the vector with given eigenvector index basis_j
    void PerturbElement(
        Element& rElement,
        const double step,
        const std::size_t basis_j,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        // Initialize is necessary for updating the section properties of shell elements
        rElement.InitializeNonLinearIteration(rCurrentProcessInfo);

        // Loop over element nodes
        for (auto& node_i : rElement.GetGeometry()) {
            auto& node_i_dofs = node_i.GetDofs();
            // Loop over nodal DOFs
            auto it_dof_i = node_i_dofs.begin();
            for (std::size_t dof_idx = 0; dof_idx < node_i_dofs.size(); dof_idx++){

                // Compute and assign the perturbation
                const double dof_perturbation = step*mFiniteDifferenceStepSize*node_i.GetValue(ROM_BASIS)(dof_idx, basis_j);
                (*it_dof_i)->GetSolutionStepValue() += dof_perturbation;
                
                // Increment the dof iterator
                ++it_dof_i;
            }
        }
        
        // Finalize is necessary for updating the section properties of shell elements
        rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
    }

    // Condition contributions
    void CalculateSystemContributions(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        //rCurrentCondition.CalculateLocalSystem(LHS_Contribution, RHS_Contribution, rCurrentProcessInfo);
    }

    void CalculateLHSContribution(
        Condition& rCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        //rCurrentCondition.CalculateLeftHandSide(LHS_Contribution, rCurrentProcessInfo);
    }

    void CalculateRHSContribution(
        Condition& rCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        //rCurrentCondition.CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);
    }

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    bool mSchemeIsInitialized;

    double mFiniteDifferenceStepSize;

    bool mFiniteDifferenceTypeFlag;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

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

}; // Class ModalDerivativeScheme

} // namespace Kratos.

#endif /* KRATOS_MODAL_DERIVATIVE_SCHEME  defined */


