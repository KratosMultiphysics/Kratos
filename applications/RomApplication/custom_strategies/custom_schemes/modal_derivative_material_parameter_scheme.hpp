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


#if !defined(KRATOS_ROM_MODAL_DERIVATIVE_MATERIAL_PARAMETER_SCHEME )
#define  KRATOS_ROM_MODAL_DERIVATIVE_MATERIAL_PARAMETER_SCHEME

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
 * @class ModalDerivativeMaterialParameterScheme
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
class ModalDerivativeMaterialParameterScheme : public ModalDerivativeScheme<TSparseSpace,TDenseSpace>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ModalDerivativeMaterialParameterScheme
    KRATOS_CLASS_POINTER_DEFINITION(ModalDerivativeMaterialParameterScheme);

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
    explicit ModalDerivativeMaterialParameterScheme(Parameters InputParameters)
    : 
    ModalDerivativeScheme<TSparseSpace,TDenseSpace>(InputParameters)
    {
        KRATOS_TRY

        std::string derivative_parameter = InputParameters["derivative_parameter"].GetString();
        if ( derivative_parameter == "density" )
        {
            mDerivativeParameter = DerivativeParameter::Density;
            mDerivativeMatrix = DerivativeMatrix::Mass;
        } 
        else if ( derivative_parameter == "poisson_ratio" )
        {
            mDerivativeParameter = DerivativeParameter::Poisson_Ratio;
            mDerivativeMatrix = DerivativeMatrix::Stiffness;
        }
        else if ( derivative_parameter == "young_modulus")
        {
            mDerivativeParameter = DerivativeParameter::Young_Modulus;
            mDerivativeMatrix = DerivativeMatrix::Stiffness;
        }
        else
        {
            KRATOS_ERROR << "Unknown derivative parameter : " << derivative_parameter << std::endl;
        }        
        
        KRATOS_CATCH("")
    }

    /// Destructor.
    ~ModalDerivativeMaterialParameterScheme() override 
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

        // Get element DOF list and resize rRHS_Contribution
        TElementDofPointersVectorType r_element_dof_list;
        rElement.GetDofList(r_element_dof_list, rCurrentProcessInfo);
        const std::size_t element_dofs_size = r_element_dof_list.size();

        // Initialize rRHS_contribution
        if (rRHS_Contribution.size() != element_dofs_size)
            rRHS_Contribution.resize(element_dofs_size, false);
        rRHS_Contribution.clear();

        switch(rCurrentProcessInfo[BUILD_LEVEL])
        {
        case 3:
            // Basis derivative RHS
            CalculateModalDerivativeRHSContribution(rElement, rRHS_Contribution, rCurrentProcessInfo);
            break;
        case 4:
            // Adjoint eigenvalue derivative RHS
            CalculateAdjointRHSContribution(rElement, rRHS_Contribution, rCurrentProcessInfo);
            break;
        case 5:
            // Adjoint sensitivity contribution
            CalculateAdjointSensitivityContribution(rElement, rRHS_Contribution, rCurrentProcessInfo);
            break;
        default:
            KRATOS_ERROR << "Invalid BUILD_LEVEL: " << rCurrentProcessInfo[BUILD_LEVEL] << std::endl;
        }

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
     * @brief This function calculates the partial RHS contribution for mode shape derivatives
     * @param rElement The element to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateModalDerivativeRHSContribution(
        Element& rElement,
        LocalSystemVectorType& rRHS_Contribution,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY

        // Derivative of basis_i
        const std::size_t basis_i = rCurrentProcessInfo[BASIS_I];

        // Get element DOF list
        const std::size_t element_dofs_size = rRHS_Contribution.size();

        // Get PhiElemental
        LocalSystemVectorType phi_elemental;
        this->GetPhiElemental(phi_elemental, basis_i, rElement, rCurrentProcessInfo);

        // Compute element matrix derivative
        LocalSystemMatrixType element_matrix_derivative(element_dofs_size, element_dofs_size);
        switch (mDerivativeParameter)
        {
        case DerivativeParameter::Density:
            // Linear dependency
            rElement.CalculateMassMatrix(element_matrix_derivative, rCurrentProcessInfo);
            element_matrix_derivative *= (-rCurrentProcessInfo[EIGENVALUE_VECTOR][basis_i] / rElement.GetProperties()(DENSITY));
            break;
        case DerivativeParameter::Poisson_Ratio:
            this->FiniteDifferencingWithMaterialParameter_LHS(rElement, POISSON_RATIO, element_matrix_derivative, rCurrentProcessInfo);
            break;
        case DerivativeParameter::Young_Modulus:
            this->FiniteDifferencingWithMaterialParameter_LHS(rElement, YOUNG_MODULUS, element_matrix_derivative, rCurrentProcessInfo);
            break;
        }

        // Compute RHS contribution
        noalias(rRHS_Contribution) -= prod(element_matrix_derivative, phi_elemental);

        KRATOS_CATCH("")
    }

    /**
     * @brief This function calculates the adjoint RHS contribution
     * @param rElement The element to compute
     * @param RHS_Contribution The RHS vector contribution
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

        // Get PhiElemental
        LocalSystemVectorType phi_elemental;
        this->GetPhiElemental(phi_elemental, basis_i, rElement, rCurrentProcessInfo);

        // Neutral state
        LocalSystemMatrixType element_LHS_initial;
        rElement.CalculateLeftHandSide(element_LHS_initial, rCurrentProcessInfo);

        // Initialize element matrices
        LocalSystemMatrixType element_LHS_p_perturbed;
        LocalSystemMatrixType element_LHS_derivative;
        
        // Loop over element DOFs
        std::size_t dof_idx = 0;
        for (auto& dof_i : r_element_dof_list)
        {
            // Positive perturbation
            // Initialize element internal quantities
            rElement.InitializeNonLinearIteration(rCurrentProcessInfo);
            dof_i->GetSolutionStepValue() += BaseType::mFiniteDifferenceStepSize;            
            // Some elements use current coordinates
            if (dof_i->GetVariable().GetSourceVariable() == DISPLACEMENT)
            {
                // auto& node_i = std::find_if(rElement.GetGeometry().begin(), 
                //                   rElement.GetGeometry().end(), 
                //                   [&dof_i](const auto& rNode )
                //     {
                //         return (rNode->GetId() == dof_i->GetId());
                //     }
                // );

                // node_i.Coordinates()[dof_i->GetVariable().GetComponentIndex()] += BaseType::mFiniteDifferenceStepSize;

                for (auto& node_i : rElement.GetGeometry())
                {
                    if (dof_i->GetId() == node_i.GetId())
                    {
                        node_i.Coordinates()[dof_i->GetVariable().GetComponentIndex()] += BaseType::mFiniteDifferenceStepSize;
                        break;
                    }
                }
            }
            // Finalize element internal quantities
            rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
            // Calculate perturbed LHS
            rElement.CalculateLeftHandSide(element_LHS_p_perturbed, rCurrentProcessInfo);

            // Reset perturbation
            // Initialize element internal quantities
            rElement.InitializeNonLinearIteration(rCurrentProcessInfo);            
            dof_i->GetSolutionStepValue() -= BaseType::mFiniteDifferenceStepSize;
            // Some elements use current coordinates
            if (dof_i->GetVariable().GetSourceVariable() == DISPLACEMENT)
            {
                for (auto& node_i : rElement.GetGeometry())
                {
                    if (dof_i->GetId() == node_i.GetId())
                    {
                        node_i.Coordinates()[dof_i->GetVariable().GetComponentIndex()] -= BaseType::mFiniteDifferenceStepSize;
                        break;
                    }
                }
            }
            // Finalize element internal quantities
            rElement.FinalizeNonLinearIteration(rCurrentProcessInfo);
            
            // Calculate element LHS derivative
            element_LHS_derivative = (element_LHS_p_perturbed - element_LHS_initial) / (BaseType::mFiniteDifferenceStepSize);
            
            // Calculate element DOF adjoint RHS contribution
            rAdjointRHS_Contribution[dof_idx] = -inner_prod(prod(element_LHS_derivative, phi_elemental), phi_elemental);
            
            dof_idx++;

        }

        this->UnlockElementNodes(rElement);

        KRATOS_CATCH("")
    }

    /**
     * @brief This function calculates the adjoint sensitivity contribution
     * @param rElement The element to compute
     * @param RHS_Contribution The RHS vector contribution
     * @param rCurrentProcessInfo The current process info instance
     */
    void CalculateAdjointSensitivityContribution(
        Element& rElement,
        LocalSystemVectorType& rAdjointSensitivityContribution,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY

        // Compute element LHS derivative
        switch (mDerivativeParameter)
        {
        case DerivativeParameter::Density:
            this->FiniteDifferencingWithMaterialParameter_RHS(rElement, DENSITY, rAdjointSensitivityContribution, rCurrentProcessInfo);
            break;
        case DerivativeParameter::Poisson_Ratio:
            this->FiniteDifferencingWithMaterialParameter_RHS(rElement, POISSON_RATIO, rAdjointSensitivityContribution, rCurrentProcessInfo);
            break;
        case DerivativeParameter::Young_Modulus:
            this->FiniteDifferencingWithMaterialParameter_RHS(rElement, YOUNG_MODULUS, rAdjointSensitivityContribution, rCurrentProcessInfo);
            break;
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This function performs finite differencing on the element RHS wrt material parameter
     * @param rElement The element to compute
     * @param rDerivativeParameter The derivative parameter
     * @param rElementRHSDerivative The element RHS derivative
     * @param rCurrentProcessInfo The current process info instance
     */
    void FiniteDifferencingWithMaterialParameter_RHS(
        Element& rElement, 
        const Variable<double>& rDerivativeParameter, 
        LocalSystemVectorType& rElementRHSDerivative,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Compute element RHS derivative
        switch (BaseType::mFiniteDifferenceType)
        {
        case BaseType::FiniteDifferenceType::Forward:
            this->ForwardDifferencingWithMaterialParameter_RHS(rElement, rDerivativeParameter, rElementRHSDerivative, rCurrentProcessInfo);
            break;
        case BaseType::FiniteDifferenceType::Central:
            this->CentralDifferencingWithMaterialParameter_RHS(rElement, rDerivativeParameter, rElementRHSDerivative, rCurrentProcessInfo);
            break;
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This function performs forward differencing on the element RHS wrt material parameter
     * @param rElement The element to compute
     * @param rDerivativeParameter The derivative parameter
     * @param rElementRHSDerivative The element RHS derivative
     * @param rCurrentProcessInfo The current process info instance
     */
    void ForwardDifferencingWithMaterialParameter_RHS(
        Element& rElement, 
        const Variable<double>& rDerivativeParameter, 
        LocalSystemVectorType& rElementRHSDerivative, 
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        if ( rElement.GetProperties().Has(rDerivativeParameter) )
        {
            // Compute initial RHS
            LocalSystemVectorType element_RHS_initial;
            rElement.CalculateRightHandSide(element_RHS_initial, rCurrentProcessInfo);

            // Save property pointer
            Properties::Pointer p_global_properties = rElement.pGetProperties();
            const double initial_property_value = p_global_properties->GetValue(rDerivativeParameter);
            const double property_value_step_size = initial_property_value * BaseType::mFiniteDifferenceStepSize;

            // Create new property and assign it to the element
            Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
            rElement.SetProperties(p_local_property);

            // Positive perturbation
            LocalSystemVectorType element_RHS_p_perturbed;
            p_local_property->SetValue(rDerivativeParameter, (initial_property_value + property_value_step_size));
            rElement.CalculateRightHandSide(element_RHS_p_perturbed, rCurrentProcessInfo);

            // Reset perturbationby giving original properties back
            rElement.SetProperties(p_global_properties);

            // Compute element RHS derivative
            noalias(rElementRHSDerivative) = -(element_RHS_p_perturbed - element_RHS_initial) / property_value_step_size;

        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This function performs central differencing on the element RHS wrt material parameter
     * @param rElement The element to compute
     * @param rDerivativeParameter The derivative parameter
     * @param rElementRHSDerivative The element RHS derivative
     * @param rCurrentProcessInfo The current process info instance
     */
    void CentralDifferencingWithMaterialParameter_RHS(
        Element& rElement, 
        const Variable<double>& rDerivativeParameter, 
        LocalSystemVectorType& rElementRHSDerivative, 
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        if ( rElement.GetProperties().Has(rDerivativeParameter) )
        {
            // Save property pointer
            Properties::Pointer p_global_properties = rElement.pGetProperties();
            const double initial_property_value = p_global_properties->GetValue(rDerivativeParameter);
            const double property_value_step_size = initial_property_value*BaseType::mFiniteDifferenceStepSize;
            
            // Create new property and assign it to the element
            Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
            rElement.SetProperties(p_local_property);

            // Positive perturbation
            LocalSystemVectorType element_RHS_p_perturbed;
            p_local_property->SetValue(rDerivativeParameter, (initial_property_value + property_value_step_size));
            rElement.CalculateRightHandSide(element_RHS_p_perturbed, rCurrentProcessInfo);

            // Negative perturbation
            LocalSystemVectorType element_RHS_m_perturbed;
            p_local_property->SetValue(rDerivativeParameter, (initial_property_value - property_value_step_size));
            rElement.CalculateRightHandSide(element_RHS_m_perturbed, rCurrentProcessInfo);

            // Reset perturbationby giving original properties back
            rElement.SetProperties(p_global_properties);

            // Compute element RHS derivative
            noalias(rElementRHSDerivative) = -(element_RHS_p_perturbed - element_RHS_m_perturbed) / (2.0*property_value_step_size);

        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This function performs finite differencing on the element LHS wrt material parameter
     * @param rElement The element to compute
     * @param rDerivativeParameter The derivative parameter
     * @param rElementLHSDerivative The element LHS derivative
     * @param rCurrentProcessInfo The current process info instance
     */
    void FiniteDifferencingWithMaterialParameter_LHS(
        Element& rElement, 
        const Variable<double>& rDerivativeParameter, 
        LocalSystemMatrixType& rElementLHSDerivative, 
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        switch (BaseType::mFiniteDifferenceType)
        {
        case BaseType::FiniteDifferenceType::Forward:
            this->ForwardDifferencingWithMaterialParameter_LHS(rElement, rDerivativeParameter, rElementLHSDerivative, rCurrentProcessInfo);
            break;
        case BaseType::FiniteDifferenceType::Central:
            this->CentralDifferencingWithMaterialParameter_LHS(rElement, rDerivativeParameter, rElementLHSDerivative, rCurrentProcessInfo);
            break;
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This function performs forward differencing on the element LHS wrt material parameter
     * @param rElement The element to compute
     * @param rDerivativeParameter The derivative parameter
     * @param rElementLHSDerivative The element LHS derivative
     * @param rCurrentProcessInfo The current process info instance
     */
    void ForwardDifferencingWithMaterialParameter_LHS(
        Element& rElement, 
        const Variable<double>& rDerivativeParameter, 
        LocalSystemMatrixType& rElementLHSDerivative, 
        const ProcessInfo& rCurrentProcessInfo)
    {

        KRATOS_TRY

        if ( rElement.GetProperties().Has(rDerivativeParameter) )
        {
            // Neutral state
            LocalSystemMatrixType element_LHS_initial;
            rElement.CalculateLeftHandSide(element_LHS_initial, rCurrentProcessInfo);
            
            // Save property pointer
            Properties::Pointer p_global_properties = rElement.pGetProperties();
            const double initial_property_value = p_global_properties->GetValue(rDerivativeParameter);
            const double property_value_step_size = initial_property_value*BaseType::mFiniteDifferenceStepSize;

            // Create new property and assign it to the element
            Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
            rElement.SetProperties(p_local_property);
            
            // Positive perturbation
            LocalSystemMatrixType element_LHS_p_perturbed;
            p_local_property->SetValue(rDerivativeParameter, (initial_property_value + property_value_step_size));
            rElement.CalculateLeftHandSide(element_LHS_p_perturbed, rCurrentProcessInfo);
            
            // Reset perturbationby giving original properties back
            rElement.SetProperties(p_global_properties);

            // Compute element matrix derivative
            noalias(rElementLHSDerivative) = (element_LHS_p_perturbed - element_LHS_initial) / property_value_step_size;
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This function performs central differencing on the element LHS wrt material parameter
     * @param rElement The element to compute
     * @param rDerivativeParameter The derivative parameter
     * @param rElementLHSDerivative The element LHS derivative
     * @param rCurrentProcessInfo The current process info instance
     */
    void CentralDifferencingWithMaterialParameter_LHS(
        Element& rElement, 
        const Variable<double>& rDerivativeParameter, 
        LocalSystemMatrixType& rElementLHSDerivative, 
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        if ( rElement.GetProperties().Has(rDerivativeParameter) )
        {
            // Save property pointer
            Properties::Pointer p_global_properties = rElement.pGetProperties();
            const double initial_property_value = p_global_properties->GetValue(rDerivativeParameter);
            const double property_value_step_size = initial_property_value*BaseType::mFiniteDifferenceStepSize;

            // Create new property and assign it to the element
            Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
            rElement.SetProperties(p_local_property);
            
            // Positive perturbation
            LocalSystemMatrixType element_LHS_p_perturbed;
            p_local_property->SetValue(rDerivativeParameter, (initial_property_value + property_value_step_size));
            rElement.CalculateLeftHandSide(element_LHS_p_perturbed, rCurrentProcessInfo);

            // Negative perturbation
            LocalSystemMatrixType element_LHS_n_perturbed;
            p_local_property->SetValue(rDerivativeParameter, (initial_property_value - property_value_step_size));
            rElement.CalculateLeftHandSide(element_LHS_n_perturbed, rCurrentProcessInfo);

            // Reset perturbationby giving original properties back
            rElement.SetProperties(p_global_properties);

            // Compute element matrix derivative
            noalias(rElementLHSDerivative) = (element_LHS_p_perturbed - element_LHS_n_perturbed) / (2.0*property_value_step_size);
        }

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

    enum DerivativeParameter {Density, Poisson_Ratio, Young_Modulus};

    DerivativeParameter mDerivativeParameter;

    enum DerivativeMatrix {Mass, Stiffness};

    DerivativeMatrix mDerivativeMatrix;

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

}; // Class ModalDerivativeMaterialParameterScheme

} // namespace Kratos.

#endif /* KRATOS_ROM_MODAL_DERIVATIVE_MATERIAL_PARAMETER_SCHEME  defined */


