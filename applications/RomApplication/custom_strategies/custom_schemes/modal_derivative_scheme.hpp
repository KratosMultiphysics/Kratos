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


#if !defined(KRATOS_ROM_MODAL_DERIVATIVE_SCHEME )
#define  KRATOS_ROM_MODAL_DERIVATIVE_SCHEME

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "utilities/entities_utilities.h"
#include "solving_strategies/schemes/scheme.h"

/* Application includes */
#include "rom_application_variables.h"

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
    /// Dof pointers vector type
    typedef std::vector<TDofType::Pointer> TElementDofPointersVectorType;

    /// Node type definition
    typedef Node<3>::NodeType TNodeType;
    /// Node pointer type definition
    typedef TNodeType::Pointer TNodePointerType;
    /// Node pointers vector type
    typedef std::vector<TNodePointerType> TNodePointerVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

        /**
     * @brief Default Constructor
     * @details Initiliazes the flags
     */
    /// Constructor.
    explicit ModalDerivativeScheme(Parameters InputParameters)
    : 
    Scheme<TSparseSpace,TDenseSpace>()
    {
        KRATOS_TRY

        std::string finite_difference_type = InputParameters["finite_difference_type"].GetString();
        if (finite_difference_type == "forward")
            mFiniteDifferenceType = FiniteDifferenceType::Forward;            
        else if (finite_difference_type == "central")
            mFiniteDifferenceType = FiniteDifferenceType::Central;
        else
            KRATOS_ERROR << "\"finite_difference_type\" can only be \"forward\" or \"central\""  << std::endl;

        mFiniteDifferenceStepSize = InputParameters["finite_difference_step_size"].GetDouble();

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~ModalDerivativeScheme() override 
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Functions totally analogous to the precedent but applied to the "condition" objects
     * @param rCondition The condition to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param rEquationIdVector The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
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

    /**
     * @brief This function is designed to calculate just the LHS contribution
     * @param rElement The element to compute
     * @param LHS_Contribution The RHS vector contribution
     * @param rEquationIdVector The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateLHSContribution(
        Element& rElement,
        LocalSystemMatrixType& rLHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        switch(rCurrentProcessInfo[BUILD_LEVEL])
        {
        case 1: // Mass matrix
            rElement.CalculateMassMatrix(rLHS_Contribution, rCurrentProcessInfo);
            break;
        case 2: // Stiffness matrix
            rElement.CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);
            // Symmetrization due to corotational elements
            rLHS_Contribution += trans(rLHS_Contribution);
            rLHS_Contribution *= 0.5;
            break;
        default:
            KRATOS_ERROR << "Invalid BUILD_LEVEL: " << rCurrentProcessInfo[BUILD_LEVEL] << std::endl;
        }
        
        rElement.EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    /**
     * @brief This function is designed to calculate just the RHS contribution
     * @param rElement The element to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param rEquationIdVector The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateRHSContribution(
        Element& rElement,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        // Derivative of basis_i wrt basis_j
        const auto basis_i = rCurrentProcessInfo[BASIS_I];
        const auto basis_j = rCurrentProcessInfo[BASIS_J];

        // Get PhiElemental
        LocalSystemVectorType phi_elemental;
        this->GetPhiElemental(phi_elemental, basis_i, rElement, rCurrentProcessInfo);
        const std::size_t num_element_dofs = phi_elemental.size();

        // Compute element LHS derivative
        // Lock element nodes for OMP parallelism
        this->LockElementNodes(rElement);
        
        // Perform FD
        Matrix element_matrix_derivative(num_element_dofs, num_element_dofs);
        element_matrix_derivative.clear();        
        switch (mFiniteDifferenceType) 
        {
            case FiniteDifferenceType::Forward:
                this->ForwardDifferencingWithBasis(element_matrix_derivative, basis_j, rElement, rCurrentProcessInfo);
                break;
            case FiniteDifferenceType::Central:
                this->CentralDifferencingWithBasis(element_matrix_derivative, basis_j, rElement, rCurrentProcessInfo);
                break;
        }
        // Unlock element nodes
        this->UnlockElementNodes(rElement);

        // Compute RHS contribution
        if (rRHS_Contribution.size() != num_element_dofs)
            rRHS_Contribution.resize(num_element_dofs);
        rRHS_Contribution.clear();
        noalias(rRHS_Contribution) -= prod(element_matrix_derivative, phi_elemental);

        rElement.EquationIdVector(EquationId,rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    /**
     * @brief Functions totally analogous to the precedent but applied to the "condition" objects
     * @param rCondition The condition to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param rEquationIdVector The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateSystemContributions(
        Condition& rCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        this->CalculateLHSContribution(rCondition, rLHS_Contribution, rEquationIdVector, rCurrentProcessInfo);

        this->CalculateRHSContribution(rCondition, rRHS_Contribution, rEquationIdVector, rCurrentProcessInfo);

        rCondition.EquationIdVector(rEquationIdVector, rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    /**
     * @brief This function is designed to calculate just the LHS contribution
     * @param rCondition The condition to compute
     * @param LHS_Contribution The RHS vector contribution
     * @param rEquationIdVector The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateLHSContribution(
        Condition& rCondition,
        LocalSystemMatrixType& rLHS_Contribution,
        Condition::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        rCondition.EquationIdVector(rEquationIdVector, rCurrentProcessInfo);

        const std::size_t size = rEquationIdVector.size();
        rLHS_Contribution.resize(size, size);
        rLHS_Contribution.clear();

        KRATOS_CATCH("")
    }

    /**
     * @brief This function is designed to calculate just the RHS contribution
     * @param rCondition The condition to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param rEquationIdVector The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateRHSContribution(
        Condition& rCondition,
        LocalSystemVectorType& rRHS_Contribution,
        Condition::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY
        
        rCondition.EquationIdVector(rEquationIdVector, rCurrentProcessInfo);

        const std::size_t size = rEquationIdVector.size();
        rRHS_Contribution.resize(size, size);
        rRHS_Contribution.clear();
        
        KRATOS_CATCH("")
    }

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

    /**
     * @brief This function locks element nodes for finite differencing in parallel.
     * @details 
     * @param rElement The element to lock nodes
     */
    void LockElementNodes(Element& rElement)
    {
        KRATOS_TRY

        auto& element_nodes = rElement.GetGeometry();

        bool all_nodes_locked{false};
        auto it_node = element_nodes.begin();
        while (!all_nodes_locked)
        {   
            auto& node_lock = it_node->GetLock();
            // try to get a lock
            if (node_lock.try_lock())
            {
                // if successful: increment the node iterator
                ++it_node;

                // check if all the nodes are locked
                if (it_node == element_nodes.end())
                    all_nodes_locked = true;
            } 
            else
            {
                // if not successful: unlock all locked nodes and empty locked nodes vector
                for (auto it_locked_node = element_nodes.begin(); it_locked_node != it_node; ++it_locked_node)
                    it_locked_node->UnSetLock();
                    
                // reset the iterator to the first node
                it_node = element_nodes.begin();
            }
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This function unlocks element nodes after finite differencing in parallel.
     * @details 
     * @param rElement The element to unlock nodes
     */
    void UnlockElementNodes(Element& rElement)
    {
        KRATOS_TRY
        // Loop over element nodes
        for (auto& node_i : rElement.GetGeometry()) 
            node_i.UnSetLock();
        KRATOS_CATCH("")
    }

    /**
     * @brief This function retrieves the basis corresponding to element DOFs
     * @details 
     * @param rPhiElemental The element to unlock nodes
     * @param BasisIndex The index of the basis to retrieve
     * @param rElement The element to retrieve basis
     * @param rElementDofList The elemental dof list
     * @param rCurrentProcessInfo The current process info instance
     */
    void GetPhiElemental(
        LocalSystemVectorType& rPhiElemental,
        const std::size_t BasisIndex,
        const Element& rElement,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY

        // Get element DOF list
        TElementDofPointersVectorType r_element_dof_list;
        rElement.GetDofList(r_element_dof_list, rCurrentProcessInfo);

        // Get elemental and nodal DOFs size
        const std::size_t num_element_dofs = r_element_dof_list.size();
        const std::size_t num_element_nodes = rElement.GetGeometry().size();
        const std::size_t num_nodal_dofs = num_element_dofs/num_element_nodes;
        
        // Get PhiElemental
        if (rPhiElemental.size() != num_element_dofs)
            rPhiElemental.resize(num_element_dofs, false);
        for(std::size_t i_node = 0; i_node < num_element_nodes; ++i_node)
        {
            const auto& r_node = rElement.GetGeometry()[i_node];
            const auto& r_phi_nodal = r_node.GetValue(ROM_BASIS);
            for (std::size_t i_dof = 0; i_dof < num_nodal_dofs; ++i_dof)
            {
                const auto& rp_dof = r_element_dof_list[i_node*num_nodal_dofs+i_dof];
                rPhiElemental[i_node*num_nodal_dofs+i_dof] = r_phi_nodal(rCurrentProcessInfo[MAP_PHI].at(rp_dof->GetVariable().Key()), BasisIndex);
            }
        }
        
        KRATOS_CATCH("")
    }

    /**
     * @brief This function performs forward differencing on the element LHS matrix
     * @details 
     * @param rElementMatrixDerivative The element to unlock nodes
     * @param BasisIndex The index of the basis for perturbation
     * @param rElement The element to compute finite differencing
     * @param rCurrentProcessInfo The current process info instance
     */
    void ForwardDifferencingWithBasis(
         LocalSystemMatrixType& rElementMatrixDerivative,
         const std::size_t BasisIndex,
         Element&rElement,
         const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const std::size_t matrix_size = rElementMatrixDerivative.size1();

        // Positive perturbation
        LocalSystemMatrixType element_matrix_p_perturbed(matrix_size, matrix_size);
        this->PerturbElementWithBasis(1.0, BasisIndex, rElement, rCurrentProcessInfo);
        rElement.CalculateLeftHandSide(element_matrix_p_perturbed, rCurrentProcessInfo);
        // Symmetrization due to corotational elements
        element_matrix_p_perturbed += trans(element_matrix_p_perturbed);
        element_matrix_p_perturbed *= 0.5;

        // Reset perturbation
        this->PerturbElementWithBasis(-1.0, BasisIndex, rElement, rCurrentProcessInfo);

        // Neutral state
        LocalSystemMatrixType element_matrix(matrix_size, matrix_size);
        rElement.CalculateLeftHandSide(element_matrix, rCurrentProcessInfo);
        // Symmetrization due to corotational elements
        element_matrix += trans(element_matrix);
        element_matrix *= 0.5;

        // Compute ElementMatrixDerivative
        noalias(rElementMatrixDerivative) = (element_matrix_p_perturbed - element_matrix) / mFiniteDifferenceStepSize;

        KRATOS_CATCH("")
    }

    /**
     * @brief This function performs central differencing on the element LHS matrix
     * @details 
     * @param rElementMatrixDerivative The element to unlock nodes
     * @param BasisIndex The index of the basis for perturbation
     * @param rElement The element to compute finite differencing
     * @param rCurrentProcessInfo The current process info instance
     */
    void CentralDifferencingWithBasis(
         LocalSystemMatrixType& rElementMatrixDerivative,
         const std::size_t BasisIndex,
         Element& rElement,
         const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        const std::size_t matrix_size = rElementMatrixDerivative.size1();

        // Positive perturbation
        LocalSystemMatrixType element_matrix_p_perturbed(matrix_size, matrix_size);
        this->PerturbElementWithBasis(1.0, BasisIndex, rElement, rCurrentProcessInfo);
        rElement.CalculateLeftHandSide(element_matrix_p_perturbed, rCurrentProcessInfo);
        // Symmetrization due to corotational elements
        element_matrix_p_perturbed += trans(element_matrix_p_perturbed);
        element_matrix_p_perturbed *= 0.5;

        // Negative perturbation
        LocalSystemMatrixType element_matrix_n_perturbed(matrix_size, matrix_size);
        this->PerturbElementWithBasis(-2.0, BasisIndex, rElement, rCurrentProcessInfo);
        rElement.CalculateLeftHandSide(element_matrix_n_perturbed, rCurrentProcessInfo);
        // Symmetrization due to corotational elements
        element_matrix_n_perturbed += trans(element_matrix_n_perturbed);
        element_matrix_n_perturbed *= 0.5;

        // Reset perturbation
        this->PerturbElementWithBasis(1.0, BasisIndex, rElement, rCurrentProcessInfo);

        // Compute LHS derivative
        noalias(rElementMatrixDerivative) = (element_matrix_p_perturbed - element_matrix_n_perturbed) / (2.0 * mFiniteDifferenceStepSize);

        KRATOS_CATCH("")
    }

    /**
     * @brief This function perturbs the element with the given BasisIndex
     * @details 
     * @param Step The step direction and size
     * @param BasisIndex The index of the basis for perturbation
     * @param rElement The element to perturb
     * @param rCurrentProcessInfo The current process info instance
     */
    // This function perturbs the element with the vector with given eigenvector index BasisIndex
    void PerturbElementWithBasis(
        const double Step,
        const std::size_t BasisIndex,
        Element& rElement,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY

        // Get DOF list
        TElementDofPointersVectorType r_element_dof_list;
        rElement.GetDofList(r_element_dof_list, rCurrentProcessInfo);

        // Get PhiElemental
        LocalSystemVectorType phi_elemental;
        this->GetPhiElemental(phi_elemental, BasisIndex, rElement, rCurrentProcessInfo);

        // Apply perturbation
        const std::size_t num_element_nodes = rElement.GetGeometry().size();
        const std::size_t num_nodal_dofs = phi_elemental.size()/num_element_nodes;
        for(std::size_t i_node = 0; i_node < num_element_nodes; ++i_node)
        {
            for (std::size_t i_dof = 0; i_dof < num_nodal_dofs; ++i_dof)
            {
                const double dof_perturbation = Step*mFiniteDifferenceStepSize*phi_elemental[i_node*num_nodal_dofs+i_dof];
                auto& rp_dof = r_element_dof_list[i_node*num_nodal_dofs+i_dof];

                // Some elements need solution step value while others need the current coordinate. Thus perturb all!
                // Update solution step value
                rp_dof->GetSolutionStepValue() += dof_perturbation;
                // Update current nodal coordinates
                if (rp_dof->GetVariable().GetSourceVariable() == DISPLACEMENT){
                    auto& r_node = rElement.GetGeometry()[i_node];
                    r_node.Coordinates()[i_dof] += dof_perturbation;
                }                
            }
        }

        // Finalize and initialize is necessary to update corotational element internals. 
        // The element is in the new configuration only after initialize
        rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
        rElement.InitializeNonLinearIteration(rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

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

    double mFiniteDifferenceStepSize;

    enum FiniteDifferenceType {Forward, Central};

    FiniteDifferenceType mFiniteDifferenceType;

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


