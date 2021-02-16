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

        // Derivative of basis_i
        const auto basis_i = rCurrentProcessInfo[BASIS_I];

        // Get element DOF list
        TElementDofPointersVectorType r_element_dof_list;
        rElement.GetDofList(r_element_dof_list, rCurrentProcessInfo);
        const std::size_t element_dofs_size = r_element_dof_list.size();

        // Get PhiElemental
        LocalSystemVectorType phi_elemental(element_dofs_size);
        this->GetPhiElemental(phi_elemental, basis_i, rElement, r_element_dof_list, rCurrentProcessInfo);

        // Compute element LHS derivative
        Matrix element_matrix_derivative(element_dofs_size, element_dofs_size);
        
        switch (mDerivativeParameter) {
            case DerivativeParameter::Density:
                // Compute element matrix
                rElement.CalculateMassMatrix(element_matrix_derivative, rCurrentProcessInfo);
                // Linear dependency
                element_matrix_derivative *= (-rCurrentProcessInfo[EIGENVALUE_VECTOR][basis_i]/rElement.GetProperties()(DENSITY));
                break;
            case DerivativeParameter::Poisson_Ratio:
                this->FiniteDifferencingWithMaterialParameter(rElement, POISSON_RATIO, element_matrix_derivative, rCurrentProcessInfo);
                break;
            case DerivativeParameter::Young_Modulus:
                this->FiniteDifferencingWithMaterialParameter(rElement, YOUNG_MODULUS, element_matrix_derivative, rCurrentProcessInfo);
                break;
        }       

        // Compute RHS contribution
        if (rRHS_Contribution.size() != element_dofs_size)
            rRHS_Contribution.resize(element_dofs_size);
        rRHS_Contribution.clear();
        noalias(rRHS_Contribution) -= prod(element_matrix_derivative, phi_elemental);

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

    void FiniteDifferencingWithMaterialParameter(
        Element& rElement, 
        const Variable<double>& rDerivativeParameter, 
        Matrix& rElementMatrixDerivative, 
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        
        switch (BaseType::mFiniteDifferenceType) {
            case BaseType::FiniteDifferenceType::Forward:
                this->ForwardDifferencingWithMaterialParameter(rElement, rDerivativeParameter, rElementMatrixDerivative, rCurrentProcessInfo);
                break;
            case BaseType::FiniteDifferenceType::Central:
                this->CentralDifferencingWithMaterialParameter(rElement, rDerivativeParameter, rElementMatrixDerivative, rCurrentProcessInfo);
                break;
        }
        
        KRATOS_CATCH("")
    }

    void ForwardDifferencingWithMaterialParameter(
        Element& rElement, 
        const Variable<double>& rDerivativeParameter, 
        Matrix& rElementMatrixDerivative, 
        const ProcessInfo& rCurrentProcessInfo)
    {

        KRATOS_TRY

        if ( rElement.GetProperties().Has(rDerivativeParameter) )
        {
            LocalSystemMatrixType element_matrix_initial;
            LocalSystemMatrixType element_matrix_p_perturbed;

            // Compute initial matrix
            rElement.CalculateLeftHandSide(element_matrix_initial, rCurrentProcessInfo);
            
            // Save property pointer
            Properties::Pointer p_global_properties = rElement.pGetProperties();
            const double initial_property_value = p_global_properties->GetValue(rDerivativeParameter);
            const double property_value_step_size = initial_property_value*BaseType::mFiniteDifferenceStepSize;

            // Create new property and assign it to the element
            Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
            rElement.SetProperties(p_local_property);
            
            // Positive perturbation
            // perturb the variable
            p_local_property->SetValue(rDerivativeParameter, (initial_property_value + property_value_step_size));
            // Compute element matrix after perturbation
            rElement.CalculateLeftHandSide(element_matrix_p_perturbed, rCurrentProcessInfo);
            
            // Give element original properties back
            rElement.SetProperties(p_global_properties);

            // Compute element matrix derivative
            noalias(rElementMatrixDerivative) = (element_matrix_p_perturbed - element_matrix_initial) / property_value_step_size;
        }

        KRATOS_CATCH("")
    }

    void CentralDifferencingWithMaterialParameter(
        Element& rElement, 
        const Variable<double>& rDerivativeParameter, 
        Matrix& rElementMatrixDerivative, 
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        if ( rElement.GetProperties().Has(rDerivativeParameter) )
        {

            LocalSystemMatrixType element_matrix_p_perturbed;
            LocalSystemMatrixType element_matrix_n_perturbed;

            // Save property pointer
            Properties::Pointer p_global_properties = rElement.pGetProperties();
            const double initial_property_value = p_global_properties->GetValue(rDerivativeParameter);
            const double property_value_step_size = initial_property_value*BaseType::mFiniteDifferenceStepSize;

            // Create new property and assign it to the element
            Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
            rElement.SetProperties(p_local_property);
            
            // Positive perturbation
            // perturb the variable
            p_local_property->SetValue(rDerivativeParameter, (initial_property_value + property_value_step_size));
            // Compute element matrix after perturbation
            rElement.CalculateLeftHandSide(element_matrix_p_perturbed, rCurrentProcessInfo);

            // Negative perturbation
            // perturb the variable
            p_local_property->SetValue(rDerivativeParameter, (initial_property_value - property_value_step_size));
            // Compute element matrix after perturbation
            rElement.CalculateLeftHandSide(element_matrix_n_perturbed, rCurrentProcessInfo);

            // Give element original properties back
            rElement.SetProperties(p_global_properties);

            // Compute element matrix derivative
            noalias(rElementMatrixDerivative) = (element_matrix_p_perturbed - element_matrix_n_perturbed) / (2.0*property_value_step_size);
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


