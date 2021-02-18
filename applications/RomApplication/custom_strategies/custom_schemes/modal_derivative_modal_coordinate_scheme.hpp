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


#if !defined(KRATOS_ROM_MODAL_DERIVATIVE_MODAL_COORDINATE_SCHEME )
#define  KRATOS_ROM_MODAL_DERIVATIVE_MODAL_COORDINATE_SCHEME

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
#include "modal_derivative_scheme.hpp"

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
 * @class ModalDerivativeModalCoordinateScheme
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
class ModalDerivativeModalCoordinateScheme : public ModalDerivativeScheme<TSparseSpace,TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ModalDerivativeModalCoordinateScheme
    KRATOS_CLASS_POINTER_DEFINITION(ModalDerivativeModalCoordinateScheme);

    // Base type definition
    typedef ModalDerivativeScheme<TSparseSpace,TDenseSpace> BaseType;

    /// Matrix type definition
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    /// Vector type definition
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    /// Local system matrix type definition
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    /// Local system vector type definition
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    /// DoF type definition
    typedef typename BaseType::TDofType TDofType;
    /// Dof pointers vector type
    typedef typename BaseType::TElementDofPointersVectorType TElementDofPointersVectorType;

    /// Node pointers vector type
    typedef typename BaseType::TNodePointerVectorType TNodePointerVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

        /**
     * @brief Default Constructor
     * @details Initiliazes the flags
     */
    /// Constructor.
    explicit ModalDerivativeModalCoordinateScheme(Parameters InputParameters)
    : 
    ModalDerivativeScheme<TSparseSpace,TDenseSpace>(InputParameters)
    {
        KRATOS_TRY

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~ModalDerivativeModalCoordinateScheme() override 
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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
        const std::size_t element_dofs_size = phi_elemental.size();

        // Initialize rRHS_Contribution
        if (rRHS_Contribution.size() != element_dofs_size)
            rRHS_Contribution.resize(element_dofs_size,false);
        rRHS_Contribution.clear();

        // Compute element LHS derivative
        // Lock element nodes for OMP parallelism
        this->LockElementNodes(rElement);
        
        Matrix element_matrix_derivative(element_dofs_size, element_dofs_size);
        // Perform FD
        switch (BaseType::mFiniteDifferenceType)
        {
        case BaseType::FiniteDifferenceType::Forward:
            this->ForwardDifferencingWithBasis_LHS(element_matrix_derivative, basis_j, rElement, rCurrentProcessInfo);
            break;
        case BaseType::FiniteDifferenceType::Central:
            this->CentralDifferencingWithBasis_LHS(element_matrix_derivative, basis_j, rElement, rCurrentProcessInfo);
            break;
        }
        // Unlock element nodes
        this->UnlockElementNodes(rElement);

        // Compute RHS contribution
        noalias(rRHS_Contribution) = -prod(element_matrix_derivative, phi_elemental);

        rElement.EquationIdVector(EquationId,rCurrentProcessInfo);

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
     * @brief This function performs forward differencing on the element LHS wrt modal coordinate
     * @details 
     * @param rElementLHSDerivative The element to unlock nodes
     * @param BasisIndex The index of the basis for perturbation
     * @param rElement The element to compute finite differencing
     * @param rCurrentProcessInfo The current process info instance
     */
    void ForwardDifferencingWithBasis_LHS(
         LocalSystemMatrixType& rElementLHSDerivative,
         const std::size_t BasisIndex,
         Element&rElement,
         const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Neutral state
        LocalSystemMatrixType element_LHS_initial;
        rElement.CalculateLeftHandSide(element_LHS_initial, rCurrentProcessInfo);

        // Positive perturbation
        LocalSystemMatrixType element_LHS_p_perturbed;
        this->PerturbElementWithBasis(1.0, BasisIndex, rElement, rCurrentProcessInfo);
        rElement.CalculateLeftHandSide(element_LHS_p_perturbed, rCurrentProcessInfo);

        // Reset perturbation
        this->PerturbElementWithBasis(-1.0, BasisIndex, rElement, rCurrentProcessInfo);

        // Compute ElementMatrixDerivative
        noalias(rElementLHSDerivative) = (element_LHS_p_perturbed - element_LHS_initial) / BaseType::mFiniteDifferenceStepSize;

        KRATOS_CATCH("")
    }

    /**
     * @brief This function performs central differencing on the element LHS wrt modal coordinate
     * @details 
     * @param rElementMatrixDerivative The element to unlock nodes
     * @param BasisIndex The index of the basis for perturbation
     * @param rElement The element to compute finite differencing
     * @param rCurrentProcessInfo The current process info instance
     */
    void CentralDifferencingWithBasis_LHS(
         LocalSystemMatrixType& rElementMatrixDerivative,
         const std::size_t BasisIndex,
         Element& rElement,
         const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Positive perturbation
        LocalSystemMatrixType element_LHS_p_perturbed;
        this->PerturbElementWithBasis(1.0, BasisIndex, rElement, rCurrentProcessInfo);
        rElement.CalculateLeftHandSide(element_LHS_p_perturbed, rCurrentProcessInfo);

        // Negative perturbation
        LocalSystemMatrixType element_LHS_n_perturbed;
        this->PerturbElementWithBasis(-2.0, BasisIndex, rElement, rCurrentProcessInfo);
        rElement.CalculateLeftHandSide(element_LHS_n_perturbed, rCurrentProcessInfo);

        // Reset perturbation
        this->PerturbElementWithBasis(1.0, BasisIndex, rElement, rCurrentProcessInfo);

        // Compute LHS derivative
        noalias(rElementMatrixDerivative) = (element_LHS_p_perturbed - element_LHS_n_perturbed) / (2.0 * BaseType::mFiniteDifferenceStepSize);

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

        // Initialize is necessary for updating the section properties of shell elements
        rElement.InitializeNonLinearIteration(rCurrentProcessInfo);

        // Get DOF list
        TElementDofPointersVectorType r_element_dof_list;
        rElement.GetDofList(r_element_dof_list, rCurrentProcessInfo);

        // Get PhiElemental
        LocalSystemVectorType phi_elemental;
        this->GetPhiElemental(phi_elemental, BasisIndex, rElement, rCurrentProcessInfo);

        // Apply perturbation
        const std::size_t nodal_dof_size = phi_elemental.size()/rElement.GetGeometry().size();
        for(std::size_t i_node = 0; i_node < rElement.GetGeometry().size(); ++i_node)
        {
            for (std::size_t i_dof = 0; i_dof < nodal_dof_size; ++i_dof)
            {
                const double dof_perturbation = Step*BaseType::mFiniteDifferenceStepSize*phi_elemental[i_node*nodal_dof_size+i_dof];
                auto& rp_dof = r_element_dof_list[i_node*nodal_dof_size+i_dof];

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

        // Finalize is necessary for updating the section properties of shell elements
        rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);

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

}; // Class ModalDerivativeModalCoordinateScheme

} // namespace Kratos.

#endif /* KRATOS_ROM_MODAL_DERIVATIVE_MODAL_COORDINATE_SCHEME  defined */


