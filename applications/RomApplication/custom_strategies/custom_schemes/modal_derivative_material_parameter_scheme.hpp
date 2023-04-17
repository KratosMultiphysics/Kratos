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
    explicit ModalDerivativeMaterialParameterScheme(Parameters InputParameters)
    : 
    ModalDerivativeScheme<TSparseSpace,TDenseSpace>(InputParameters)
    {
        KRATOS_TRY

        std::string derivative_parameter = InputParameters["derivative_parameter"].GetString();
        if ( derivative_parameter == "density" )
        {
            mpDerivativeParameter = &DENSITY;
            mDerivativeMatrixType = DerivativeMatrixType::Mass;
        } 
        else if ( derivative_parameter == "poisson_ratio" )
        {
            mpDerivativeParameter = &POISSON_RATIO;
            mDerivativeMatrixType = DerivativeMatrixType::Stiffness;
        }
        else if ( derivative_parameter == "young_modulus")
        {
            mpDerivativeParameter = &YOUNG_MODULUS;
            mDerivativeMatrixType = DerivativeMatrixType::Stiffness;
        }
        else if (derivative_parameter == "truss_prestress_pk2")
        {
            mpDerivativeParameter = &TRUSS_PRESTRESS_PK2;
            mDerivativeMatrixType = DerivativeMatrixType::Stiffness;
        }
        else
            KRATOS_ERROR << "Unknown derivative parameter : " << derivative_parameter << std::endl;
        
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
        Element::EquationIdVectorType& rEquationIdVector,
        const ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        rElement.EquationIdVector(rEquationIdVector,rCurrentProcessInfo);

        // Resize RHS contribution
        const std::size_t num_element_dofs = rEquationIdVector.size();
        if (rRHS_Contribution.size() != num_element_dofs)
            rRHS_Contribution.resize(num_element_dofs);
        rRHS_Contribution.clear();

        // Basis derivative RHS contribution
        CalculateModalDerivativeRHSContribution(rElement, rRHS_Contribution, rCurrentProcessInfo);

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
    
    Variable<double>* mpDerivativeParameter = nullptr;

    enum DerivativeMatrixType {Mass, Stiffness};

    DerivativeMatrixType mDerivativeMatrixType;

    ///@}
    ///@name Protected Operators
    ///@{

    /**
     * @brief This function calculates the partial RHS contribution for mode shape derivatives
     * @param rElement The element to compute
     * @param rRHS_Contribution The RHS vector contribution
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
        const std::size_t num_element_dofs = rRHS_Contribution.size();

        // Get PhiElemental
        LocalSystemVectorType phi_elemental;
        this->GetPhiElemental(phi_elemental, basis_i, rElement, rCurrentProcessInfo);

        // Compute element matrix derivative
        LocalSystemMatrixType element_matrix_derivative(num_element_dofs, num_element_dofs);
        element_matrix_derivative.clear();

        if (*mpDerivativeParameter == DENSITY)
        {
            // Linear dependency
            rElement.CalculateMassMatrix(element_matrix_derivative, rCurrentProcessInfo);
            // Symmetrization due to corotational elements
            element_matrix_derivative += trans(element_matrix_derivative);
            element_matrix_derivative *= 0.5;
            element_matrix_derivative *= (-rCurrentProcessInfo[EIGENVALUE_VECTOR][basis_i] / rElement.GetProperties()(DENSITY));
        }
        else if (*mpDerivativeParameter == YOUNG_MODULUS  || *mpDerivativeParameter == POISSON_RATIO || *mpDerivativeParameter == TRUSS_PRESTRESS_PK2)
        {
            this->FiniteDifferencingWithMaterialParameter_LHS(rElement, element_matrix_derivative, rCurrentProcessInfo);
        }
        else
            KRATOS_ERROR << "Unknown derivative parameter : " << *mpDerivativeParameter << std::endl;

        // Compute RHS contribution
        noalias(rRHS_Contribution) = -prod(element_matrix_derivative, phi_elemental);

        KRATOS_CATCH("")
    }

    /**
     * @brief This function performs finite differencing on the element LHS wrt material parameter
     * @param rElement The element to compute
     * @param rElementLHSDerivative The element LHS derivative
     * @param rCurrentProcessInfo The current process info instance
     */
    void FiniteDifferencingWithMaterialParameter_LHS(
        Element& rElement, 
        LocalSystemMatrixType& rElementLHSDerivative, 
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Compute element LHS derivative
        if ( rElement.GetProperties().Has(*mpDerivativeParameter) )
        {
            switch (BaseType::mFiniteDifferenceType)
            {
            case BaseType::FiniteDifferenceType::Forward:
                this->ForwardDifferencingWithMaterialParameter_LHS(rElement, rElementLHSDerivative, rCurrentProcessInfo);
                break;
            case BaseType::FiniteDifferenceType::Central:
                this->CentralDifferencingWithMaterialParameter_LHS(rElement, rElementLHSDerivative, rCurrentProcessInfo);
                break;
            }        
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This function performs forward differencing on the element LHS wrt material parameter
     * @param rElement The element to compute
     * @param rElementLHSDerivative The element LHS derivative
     * @param rCurrentProcessInfo The current process info instance
     */
    void ForwardDifferencingWithMaterialParameter_LHS(
        Element& rElement, 
        LocalSystemMatrixType& rElementLHSDerivative, 
        const ProcessInfo& rCurrentProcessInfo)
    {

        KRATOS_TRY

        // Neutral state
        LocalSystemMatrixType element_LHS_initial;
        rElement.CalculateLeftHandSide(element_LHS_initial, rCurrentProcessInfo);
        // Symmetrization due to corotational elements
        element_LHS_initial += trans(element_LHS_initial);
        element_LHS_initial *= 0.5;

        // Save property pointer
        Properties::Pointer p_global_properties = rElement.pGetProperties();
        const double initial_property_value = p_global_properties->GetValue(*mpDerivativeParameter);
        const double property_value_step_size = initial_property_value * BaseType::mFiniteDifferenceStepSize;

        // Create new property and assign it to the element
        Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
        rElement.SetProperties(p_local_property);

        // Positive perturbation
        LocalSystemMatrixType element_LHS_p_perturbed;
        p_local_property->SetValue(*mpDerivativeParameter, (initial_property_value + property_value_step_size));
        rElement.CalculateLeftHandSide(element_LHS_p_perturbed, rCurrentProcessInfo);
        // Symmetrization due to corotational elements
        element_LHS_p_perturbed += trans(element_LHS_p_perturbed);
        element_LHS_p_perturbed *= 0.5;

        // Reset perturbationby giving original properties back
        rElement.SetProperties(p_global_properties);

        // Compute element matrix derivative
        noalias(rElementLHSDerivative) = (element_LHS_p_perturbed - element_LHS_initial) / property_value_step_size;

        KRATOS_CATCH("")
    }

    /**
     * @brief This function performs central differencing on the element LHS wrt material parameter
     * @param rElement The element to compute
     * @param rElementLHSDerivative The element LHS derivative
     * @param rCurrentProcessInfo The current process info instance
     */
    void CentralDifferencingWithMaterialParameter_LHS(
        Element& rElement, 
        LocalSystemMatrixType& rElementLHSDerivative, 
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Save property pointer
        Properties::Pointer p_global_properties = rElement.pGetProperties();
        const double initial_property_value = p_global_properties->GetValue(*mpDerivativeParameter);
        const double property_value_step_size = initial_property_value * BaseType::mFiniteDifferenceStepSize;

        // Create new property and assign it to the element
        Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
        rElement.SetProperties(p_local_property);

        // Positive perturbation
        LocalSystemMatrixType element_LHS_p_perturbed;
        p_local_property->SetValue(*mpDerivativeParameter, (initial_property_value + property_value_step_size));
        rElement.CalculateLeftHandSide(element_LHS_p_perturbed, rCurrentProcessInfo);
        // Symmetrization due to corotational elements
        element_LHS_p_perturbed += trans(element_LHS_p_perturbed);
        element_LHS_p_perturbed *= 0.5;

        // Negative perturbation
        LocalSystemMatrixType element_LHS_n_perturbed;
        p_local_property->SetValue(*mpDerivativeParameter, (initial_property_value - property_value_step_size));
        rElement.CalculateLeftHandSide(element_LHS_n_perturbed, rCurrentProcessInfo);
        // Symmetrization due to corotational elements
        element_LHS_n_perturbed += trans(element_LHS_n_perturbed);
        element_LHS_n_perturbed *= 0.5;

        // Reset perturbationby giving original properties back
        rElement.SetProperties(p_global_properties);

        // Compute element matrix derivative
        noalias(rElementLHSDerivative) = (element_LHS_p_perturbed - element_LHS_n_perturbed) / (2.0 * property_value_step_size);

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

}; // Class ModalDerivativeMaterialParameterScheme

} // namespace Kratos.

#endif /* KRATOS_ROM_MODAL_DERIVATIVE_MATERIAL_PARAMETER_SCHEME  defined */


