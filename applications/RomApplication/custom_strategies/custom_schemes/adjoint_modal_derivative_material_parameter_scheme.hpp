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


#if !defined(KRATOS_ROM_ADJOINT_MODAL_DERIVATIVE_MATERIAL_PARAMETER_SCHEME )
#define  KRATOS_ROM_ADJOINT_MODAL_DERIVATIVE_MATERIAL_PARAMETER_SCHEME

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
 * @class AdjointModalDerivativeMaterialParameterScheme
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
class AdjointModalDerivativeMaterialParameterScheme : public ModalDerivativeMaterialParameterScheme<TSparseSpace,TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AdjointModalDerivativeMaterialParameterScheme
    KRATOS_CLASS_POINTER_DEFINITION(AdjointModalDerivativeMaterialParameterScheme);

    // Base type definition
    typedef ModalDerivativeMaterialParameterScheme<TSparseSpace,TDenseSpace> BaseType;

    /// Local system matrix type definition
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
    /// Local system vector type definition
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    /// Dof pointers vector type
    typedef typename BaseType::TElementDofPointersVectorType TElementDofPointersVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

        /**
     * @brief Default Constructor
     * @details Initiliazes the flags
     */
    /// Constructor.
    explicit AdjointModalDerivativeMaterialParameterScheme(Parameters InputParameters)
    : 
    ModalDerivativeMaterialParameterScheme<TSparseSpace,TDenseSpace>(InputParameters)
    {
    }

    /// Destructor.
    ~AdjointModalDerivativeMaterialParameterScheme() override 
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function is designed to calculate just the LHS contribution
     * @param rElement The element to compute
     * @param rLHS_Contribution The LHS matrix contribution
     * @param rEquationIdVector The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateLHSContribution(
        Element& rElement,
        LocalSystemMatrixType& rLHS_Contribution,
        Element::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        switch(rCurrentProcessInfo[BUILD_LEVEL_ROM])
        {
        case 1: // Mass matrix
            rElement.CalculateMassMatrix(rLHS_Contribution, rCurrentProcessInfo);
            break;
        case 2: // Stiffness matrix
            rElement.CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);
            
            break;
        case 3: // Adjoint LHS matrix
        {    
            rElement.CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);
            break;
        }
        default:
            KRATOS_ERROR << "Invalid BUILD_LEVEL_ROM: " << rCurrentProcessInfo[BUILD_LEVEL_ROM] << std::endl;
        }

        // Symmetrization due to corotational elements
        rLHS_Contribution += trans(rLHS_Contribution);
        rLHS_Contribution *= 0.5;

        rElement.EquationIdVector(rEquationIdVector, rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    /**
     * @brief This function is designed to calculate just the RHS contribution
     * @param rElement The element to compute
     * @param rRHS_Contribution The RHS vector contribution
     * @param rEquationIdVector The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateRHSContribution(
        Element& rElement,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        if (rCurrentProcessInfo[BUILD_LEVEL_ROM] == 4 || rCurrentProcessInfo[BUILD_LEVEL_ROM] == 6)
        {
            rElement.EquationIdVector(rEquationIdVector,rCurrentProcessInfo);
            const std::size_t element_dofs_size = rEquationIdVector.size();

            // Initialize rRHS_contribution
            if (rRHS_Contribution.size() != element_dofs_size)
                rRHS_Contribution.resize(element_dofs_size, false);
            rRHS_Contribution.clear();

            if (rCurrentProcessInfo[BUILD_LEVEL_ROM] == 4)
            {
                // Adjoint eigenvalue derivative RHS
                CalculateAdjointRHSContribution(rElement, rRHS_Contribution, rCurrentProcessInfo);
            }
            else if (rCurrentProcessInfo[BUILD_LEVEL_ROM] == 6)
            {
                // Adjoint sensitivity contribution
                CalculateAdjointSensitivityContribution(rElement, rRHS_Contribution, rCurrentProcessInfo);
            }
        }
        else if (rCurrentProcessInfo[BUILD_LEVEL_ROM] == 5)
        {
            // Basis derivative RHS
            BaseType::CalculateRHSContribution(rElement, rRHS_Contribution, rEquationIdVector, rCurrentProcessInfo);
        }
        else
            KRATOS_ERROR << "Invalid BUILD_LEVEL_ROM: " << rCurrentProcessInfo[BUILD_LEVEL_ROM] << std::endl;

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

    /**
     * @brief This function calculates the adjoint RHS contribution
     * @param rElement The element to compute
     * @param rAdjointRHS_Contribution The adjoint RHS vector contribution
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateAdjointRHSContribution(
        Element& rElement,
        LocalSystemVectorType& rAdjointRHS_Contribution,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY

        // Lock element nodes for OMP parallelism
        this->LockElementNodes(rElement);

        // Derivative of basis_i
        const std::size_t basis_i = rCurrentProcessInfo[BASIS_I];

        // Get element DOF list
        TElementDofPointersVectorType r_element_dof_list;
        rElement.GetDofList(r_element_dof_list, rCurrentProcessInfo);
        const std::size_t num_element_dofs = r_element_dof_list.size();

        // Get PhiElemental
        LocalSystemVectorType phi_elemental;
        this->GetPhiElemental(phi_elemental, basis_i, rElement, rCurrentProcessInfo);

        // Initialize element matrices
        LocalSystemMatrixType element_LHS_p_perturbed;
        LocalSystemMatrixType element_LHS_m_perturbed;
        LocalSystemMatrixType element_LHS_derivative;
        element_LHS_p_perturbed.resize(num_element_dofs,num_element_dofs,false);
        element_LHS_m_perturbed.resize(num_element_dofs,num_element_dofs,false);
        element_LHS_derivative.resize(num_element_dofs,num_element_dofs,false);

        // Loop over element nodes
        const std::size_t num_element_nodes = rElement.GetGeometry().size();
        auto& r_element_nodes = rElement.GetGeometry();
        const std::size_t num_nodal_dofs = num_element_dofs/num_element_nodes;
        std::size_t dof_component_index;
        std::size_t dof_index;
        for (std::size_t i_node = 0; i_node < num_element_nodes; ++i_node)
        {
            auto& r_node = r_element_nodes[i_node];
            
            // Loop over nodal dofs of this element
            for (std::size_t i_dof = 0; i_dof < num_nodal_dofs; ++i_dof)
            {
                dof_index = i_node*num_nodal_dofs+i_dof;
                auto& rp_dof = r_element_dof_list[dof_index];
                
                if (rp_dof->IsFree())
                {
                    dof_component_index = rp_dof->GetVariable().GetComponentIndex();
                    const double dof_value = rp_dof->GetSolutionStepValue();
                    
                    // Positive perturbation
                    rp_dof->GetSolutionStepValue() += BaseType::mFiniteDifferenceStepSize;
                    if (rp_dof->GetVariable().GetSourceVariable() == DISPLACEMENT)
                        r_node.Coordinates()[dof_component_index] = r_node.GetInitialPosition()[dof_component_index] + rp_dof->GetSolutionStepValue();
                    rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
                    rElement.InitializeNonLinearIteration(rCurrentProcessInfo);
                    element_LHS_p_perturbed.clear();
                    rElement.CalculateLeftHandSide(element_LHS_p_perturbed, rCurrentProcessInfo);
                    // Symmetrization due to corotational elements
                    element_LHS_p_perturbed += trans(element_LHS_p_perturbed);
                    element_LHS_p_perturbed *= 0.5;

                    // Negative perturbation                    
                    rp_dof->GetSolutionStepValue() -= 2.0*BaseType::mFiniteDifferenceStepSize;
                    if (rp_dof->GetVariable().GetSourceVariable() == DISPLACEMENT)
                        r_node.Coordinates()[dof_component_index] = r_node.GetInitialPosition()[dof_component_index] + rp_dof->GetSolutionStepValue();
                    rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
                    rElement.InitializeNonLinearIteration(rCurrentProcessInfo);
                    element_LHS_m_perturbed.clear();
                    rElement.CalculateLeftHandSide(element_LHS_m_perturbed, rCurrentProcessInfo);
                    // Symmetrization due to corotational elements
                    element_LHS_m_perturbed += trans(element_LHS_m_perturbed);
                    element_LHS_m_perturbed *= 0.5;
                    
                    // Reset Perturbation
                    rp_dof->GetSolutionStepValue() = dof_value;
                    if (rp_dof->GetVariable().GetSourceVariable() == DISPLACEMENT)
                        r_node.Coordinates()[dof_component_index] = r_node.GetInitialPosition()[dof_component_index] + rp_dof->GetSolutionStepValue();
                    rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
                    rElement.InitializeNonLinearIteration(rCurrentProcessInfo);
                    
                    // Compute LHS derivative
                    element_LHS_derivative.clear();
                    noalias(element_LHS_derivative) = ((element_LHS_p_perturbed - element_LHS_m_perturbed) / (2.0*BaseType::mFiniteDifferenceStepSize));
                    
                    // Compute adjoint RHS contribution
                    rAdjointRHS_Contribution[dof_index] = -inner_prod(prod(element_LHS_derivative, phi_elemental), phi_elemental);
                }

            }
        }

        this->UnlockElementNodes(rElement);

        KRATOS_CATCH("")
    }

    /**
     * @brief This function calculates the adjoint sensitivity contribution
     * @param rElement The element to compute
     * @param rAdjointSensitivityContribution The adjoint sensitivity contribution
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateAdjointSensitivityContribution(
        Element& rElement,
        LocalSystemVectorType& rAdjointSensitivityContribution,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY

        // Compute element RHS derivative
        if ( rElement.GetProperties().Has(*BaseType::mpDerivativeParameter) )
        {            
            switch (BaseType::mFiniteDifferenceType)
            {
            case BaseType::FiniteDifferenceType::Forward:
                this->ForwardDifferencingWithMaterialParameter_RHS(rElement, rAdjointSensitivityContribution, rCurrentProcessInfo);
                break;
            case BaseType::FiniteDifferenceType::Central:
                this->CentralDifferencingWithMaterialParameter_RHS(rElement, rAdjointSensitivityContribution, rCurrentProcessInfo);
                break;
            }
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This function performs forward differencing on the element RHS wrt material parameter
     * @param rElement The element to compute
     * @param rElementRHSDerivative The element RHS derivative
     * @param rCurrentProcessInfo The current process info instance
     */
    void ForwardDifferencingWithMaterialParameter_RHS(
        Element& rElement, 
        LocalSystemVectorType& rElementRHSDerivative, 
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Compute initial RHS
        LocalSystemVectorType element_RHS_initial;
        rElement.CalculateRightHandSide(element_RHS_initial, rCurrentProcessInfo);

        // Save property pointer
        Properties::Pointer p_global_properties = rElement.pGetProperties();
        const double initial_property_value = p_global_properties->GetValue(*BaseType::mpDerivativeParameter);
        const double property_value_step_size = initial_property_value * BaseType::mFiniteDifferenceStepSize;

        // Create new property and assign it to the element
        Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
        rElement.SetProperties(p_local_property);

        // Positive perturbation
        LocalSystemVectorType element_RHS_p_perturbed;
        p_local_property->SetValue(*BaseType::mpDerivativeParameter, (initial_property_value + property_value_step_size));
        rElement.CalculateRightHandSide(element_RHS_p_perturbed, rCurrentProcessInfo);

        // Reset perturbationby giving original properties back
        rElement.SetProperties(p_global_properties);

        // Compute element RHS derivative
        noalias(rElementRHSDerivative) = (-(element_RHS_p_perturbed - element_RHS_initial) / property_value_step_size);

        KRATOS_CATCH("")
    }

    /**
     * @brief This function performs central differencing on the element RHS wrt material parameter
     * @param rElement The element to compute
     * @param rElementRHSDerivative The element RHS derivative
     * @param rCurrentProcessInfo The current process info instance
     */
    void CentralDifferencingWithMaterialParameter_RHS(
        Element& rElement, 
        LocalSystemVectorType& rElementRHSDerivative, 
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Save property pointer
        Properties::Pointer p_global_properties = rElement.pGetProperties();
        const double initial_property_value = p_global_properties->GetValue(*BaseType::mpDerivativeParameter);
        const double property_value_step_size = initial_property_value * BaseType::mFiniteDifferenceStepSize;

        // Create new property and assign it to the element
        Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
        rElement.SetProperties(p_local_property);

        // Positive perturbation
        LocalSystemVectorType element_RHS_p_perturbed;
        p_local_property->SetValue(*BaseType::mpDerivativeParameter, (initial_property_value + property_value_step_size));
        rElement.CalculateRightHandSide(element_RHS_p_perturbed, rCurrentProcessInfo);

        // Negative perturbation
        LocalSystemVectorType element_RHS_m_perturbed;
        p_local_property->SetValue(*BaseType::mpDerivativeParameter, (initial_property_value - property_value_step_size));
        rElement.CalculateRightHandSide(element_RHS_m_perturbed, rCurrentProcessInfo);

        // Reset perturbationby giving original properties back
        rElement.SetProperties(p_global_properties);

        // Compute element RHS derivative
        noalias(rElementRHSDerivative) = (-(element_RHS_p_perturbed - element_RHS_m_perturbed) / (2.0 * property_value_step_size));

        KRATOS_CATCH("")
    }

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

}; // Class AdjointModalDerivativeMaterialParameterScheme

} // namespace Kratos.

#endif /* KRATOS_ROM_ADJOINT_MODAL_DERIVATIVE_MATERIAL_PARAMETER_SCHEME  defined */


