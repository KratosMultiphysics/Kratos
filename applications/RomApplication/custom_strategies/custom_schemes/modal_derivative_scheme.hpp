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
    ModalDerivativeScheme(Variable<double>& DerivativeParameter, Parameters InputParameters)
    : 
    Scheme<TSparseSpace,TDenseSpace>(),
    mDerivativeParameter(DerivativeParameter)
    {
        KRATOS_TRY

        std::string finite_difference_type = InputParameters["finite_difference_type"].GetString();
        if (finite_difference_type == "central")
            mFiniteDifferenceTypeFlag = true;
        else if (finite_difference_type == "forward")
            mFiniteDifferenceTypeFlag = false;
        else
            KRATOS_ERROR << "\"finite_difference_type\" can only be \"forward\" or \"central\""  << std::endl;

        if (mDerivativeParameter == DENSITY)
            mDerivativeMatrixType = true;
        else if (mDerivativeParameter == MODAL_COORDINATE || mDerivativeParameter == YOUNG_MODULUS || mDerivativeParameter == POISSON_RATIO)
            mDerivativeMatrixType = false;
        else
            KRATOS_ERROR;
        
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
        
        // Derivative of basis_i
        const std::size_t basis_i = rCurrentProcessInfo[BASIS_I];

        // Get elementalGlobalDofSize
        std::size_t elementalGlobalDofSize = 0;
        for (auto& node_i : rElement.GetGeometry())
            elementalGlobalDofSize += node_i.GetDofs().size();
        
        // Get elementalLocalDofSize
        std::vector<Dof<double>::Pointer> rElementalDofList;
        rElement.GetDofList(rElementalDofList, rCurrentProcessInfo);
        std::size_t elementalLocalDofSize = rElementalDofList.size();

        // Compute element LHS derivative
        Matrix element_matrix_derivative;
        element_matrix_derivative.resize(elementalLocalDofSize,elementalLocalDofSize,false);
        if (mDerivativeParameter == MODAL_COORDINATE)
        {   // Modal parameter            

            // Derivative wrt basis_j
            const std::size_t basis_j = rCurrentProcessInfo[BASIS_J];

            if (mFiniteDifferenceTypeFlag) // Central Difference
                this->CentralDifferencingWithBasis(rElement, element_matrix_derivative, basis_j, rCurrentProcessInfo);
            else // Forward difference
                this->ForwardDifferencingWithBasis(rElement, element_matrix_derivative, basis_j, rCurrentProcessInfo);

        }
        else if (mDerivativeParameter == DENSITY || mDerivativeParameter == YOUNG_MODULUS || mDerivativeParameter == POISSON_RATIO)
        {   // Material parameter
            
            if (mFiniteDifferenceTypeFlag) // Central Difference
                this->CentralDifferencingWithMaterialParameter(rElement, element_matrix_derivative, mDerivativeParameter, rCurrentProcessInfo);
            else // Forward difference
                this->ForwardDifferencingWithMaterialParameter(rElement, element_matrix_derivative, mDerivativeParameter, rCurrentProcessInfo);

        }

        // Create PhiElementalGlobal 
        LocalSystemVectorType PhiElementalGlobal;
        PhiElementalGlobal.resize(elementalGlobalDofSize);
        // Get PhiElemental
        std::size_t elem_glob_dof_ctr = 0;
        // Loop over nodes
        for (auto& node_i : rElement.GetGeometry()) {
            auto& node_i_dofs = node_i.GetDofs();
            
            const Matrix *pPhiNodal = &(node_i.GetValue(ROM_BASIS));
            for (std::size_t node_i_glob_dof_ctr = 0; node_i_glob_dof_ctr < pPhiNodal->size1(); node_i_glob_dof_ctr++)
            {
                PhiElementalGlobal[elem_glob_dof_ctr + node_i_glob_dof_ctr] = (*pPhiNodal)(node_i_glob_dof_ctr, basis_i);
            }

            elem_glob_dof_ctr += node_i_dofs.size();
        }

        // Create PhiElementalLocal
        LocalSystemVectorType PhiElementalLocal;
        PhiElementalLocal.resize(elementalLocalDofSize);        

        // Initialize RHS contribution
        rRHS_Contribution.resize(elementalLocalDofSize);
        rRHS_Contribution.clear();

        // Build RHS contribution
        if (elementalGlobalDofSize == elementalLocalDofSize && elementalGlobalDofSize / rElement.GetGeometry().size() == 3)
        {   // there are only disp dofs both in global and in local dof sets
            rRHS_Contribution -= prod(element_matrix_derivative, PhiElementalGlobal);
        } 
        else if (elementalGlobalDofSize == elementalLocalDofSize && elementalGlobalDofSize / rElement.GetGeometry().size() == 6)
        {   // there are disp and rot dofs both in global and in local dof sets. reorder PhiElementalGlobal into PhiElementalLocal
            const std::size_t disp_shifter = 3;
            std::size_t dof_shifter = 0;
            for (auto& node_i : rElement.GetGeometry()) {
                auto& node_i_dofs = node_i.GetDofs();
                
                for (auto& dof_i : node_i_dofs) {
                    if (dof_i->GetVariable().GetSourceVariable() == DISPLACEMENT){
                        PhiElementalLocal[dof_shifter + dof_i->GetVariable().GetComponentIndex()] = PhiElementalGlobal(dof_shifter + disp_shifter + dof_i->GetVariable().GetComponentIndex());
                    } else {
                        PhiElementalLocal[dof_shifter + disp_shifter + dof_i->GetVariable().GetComponentIndex()] = PhiElementalGlobal(dof_shifter + dof_i->GetVariable().GetComponentIndex());
                    }
                }
                dof_shifter += node_i_dofs.size();
            }
            rRHS_Contribution -= prod(element_matrix_derivative, PhiElementalLocal);

        } else if (elementalGlobalDofSize != elementalLocalDofSize) 
        {   // there are disp dofs in element dof set and both disp and rot dofs in the global dof set. filter out only disp dofs
            const std::size_t disp_shifter = 3;
            const std::size_t nodal_dof_size = elementalGlobalDofSize / rElement.GetGeometry().size();
            for (std::size_t iNode = 0; iNode < rElement.GetGeometry().size(); iNode++)
            {
                for (std::size_t iXYZ = 0; iXYZ < 3; iXYZ++)
                {
                    PhiElementalLocal[iNode * 3 + iXYZ] = PhiElementalGlobal[iNode * nodal_dof_size + disp_shifter + iXYZ];
                }
            }

            rRHS_Contribution -= prod(element_matrix_derivative, PhiElementalLocal);
        }

        if (mDerivativeMatrixType)
            rRHS_Contribution *= -rCurrentProcessInfo[EIGENVALUE_VECTOR][basis_i];

        rElement.EquationIdVector(EquationId,rCurrentProcessInfo);

        KRATOS_CATCH("")
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

    // This function locks element nodes for finite differencing
    void LockElementNodes(Element& rElement)
    {
        KRATOS_TRY

        // Loop over element nodes
        for (auto& node_i : rElement.GetGeometry()) 
            node_i.SetLock();

        KRATOS_CATCH("")
    }

    // This function locks element nodes for finite differencing
    void UnlockElementNodes(Element& rElement)
    {
        KRATOS_TRY
        // Loop over element nodes
        for (auto& node_i : rElement.GetGeometry()) 
            node_i.UnSetLock();
        KRATOS_CATCH("")
    }

    void ForwardDifferencingWithMaterialParameter(
        Element& rElement,
        LocalSystemMatrixType& rElementMatrixDerivative,
        const Variable<double>& rMaterialParameter,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        if ( rElement.GetProperties().Has(rMaterialParameter) )
        {

            LocalSystemMatrixType element_matrix_initial;
            LocalSystemMatrixType element_matrix_p_perturbed;

            // Compute initial matrix
            if (mDerivativeMatrixType)
                rElement.CalculateMassMatrix(element_matrix_initial, rCurrentProcessInfo);
            else
                rElement.CalculateLeftHandSide(element_matrix_initial, rCurrentProcessInfo);
            
            // Save property pointer
            Properties::Pointer p_global_properties = rElement.pGetProperties();
            const double initial_property_value = p_global_properties->GetValue(rMaterialParameter);
            const double property_value_step_size = initial_property_value*mFiniteDifferenceStepSize;

            // Create new property and assign it to the element
            Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
            rElement.SetProperties(p_local_property);
            
            // Positive perturbation
            // perturb the variable
            p_local_property->SetValue(rMaterialParameter, (initial_property_value + property_value_step_size));
            // Compute element matrix after perturbation
            if (mDerivativeMatrixType)
                rElement.CalculateMassMatrix(element_matrix_p_perturbed, rCurrentProcessInfo);
            else
                rElement.CalculateLeftHandSide(element_matrix_p_perturbed, rCurrentProcessInfo);
            
            // Give element original properties back
            rElement.SetProperties(p_global_properties);

            noalias(rElementMatrixDerivative) = (element_matrix_p_perturbed - element_matrix_initial) / property_value_step_size;
        }
        KRATOS_CATCH("")
    }

    void CentralDifferencingWithMaterialParameter(
        Element& rElement,
        LocalSystemMatrixType& rElementMatrixDerivative,
        const Variable<double>& rMaterialParameter,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        if ( rElement.GetProperties().Has(rMaterialParameter) )
        {

            LocalSystemMatrixType element_matrix_p_perturbed;
            LocalSystemMatrixType element_matrix_n_perturbed;

            // Save property pointer
            Properties::Pointer p_global_properties = rElement.pGetProperties();
            const double initial_property_value = rElement.GetProperties()[rMaterialParameter];
            const double property_value_step_size = initial_property_value*mFiniteDifferenceStepSize;

            // Create new property and assign it to the element
            Properties::Pointer p_local_property(Kratos::make_shared<Properties>(Properties(*p_global_properties)));
            rElement.SetProperties(p_local_property);
            
            // Positive perturbation
            // perturb the variable
            p_local_property->SetValue(rMaterialParameter, (initial_property_value + property_value_step_size));
            // Compute element matrix after perturbation
            if (mDerivativeMatrixType)
                rElement.CalculateMassMatrix(element_matrix_p_perturbed, rCurrentProcessInfo);
            else
                rElement.CalculateLeftHandSide(element_matrix_p_perturbed, rCurrentProcessInfo);

            // Negative perturbation
            // perturb the variable
            p_local_property->SetValue(rMaterialParameter, (initial_property_value - property_value_step_size));
            // Compute element matrix after perturbation
            if (mDerivativeMatrixType)
                rElement.CalculateMassMatrix(element_matrix_n_perturbed, rCurrentProcessInfo);
            else
                rElement.CalculateLeftHandSide(element_matrix_n_perturbed, rCurrentProcessInfo);

            // Give element original properties back
            rElement.SetProperties(p_global_properties);

            noalias(rElementMatrixDerivative) = (element_matrix_p_perturbed - element_matrix_n_perturbed) / (2.0*property_value_step_size);

        }

        KRATOS_CATCH("")
        
    }

    void ForwardDifferencingWithBasis(
         Element&rElement,
         LocalSystemMatrixType& rElementMatrixDerivative,
         const std::size_t basis_j,
         const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Lock element nodes for OMP parallelism
        this->LockElementNodes(rElement);

        std::size_t matrix_size = rElementMatrixDerivative.size1();

        // Positive perturbation
        LocalSystemMatrixType element_matrix_p_perturbed;
        element_matrix_p_perturbed.resize(matrix_size, matrix_size, false);
        this->PerturbElementWithBasis(rElement, 1.0, basis_j, rCurrentProcessInfo);
        rElement.CalculateLeftHandSide(element_matrix_p_perturbed, rCurrentProcessInfo);

        // Reset perturbation
        this->PerturbElementWithBasis(rElement, -1.0, basis_j, rCurrentProcessInfo);

        // Neutral state
        LocalSystemMatrixType element_matrix;
        element_matrix.resize(matrix_size, matrix_size, false);
        rElement.CalculateLeftHandSide(element_matrix, rCurrentProcessInfo);

        // Compute ElementMatrixDerivative
        noalias(rElementMatrixDerivative) = (element_matrix_p_perturbed - element_matrix) / mFiniteDifferenceStepSize;

        // Unlock element nodes
        this->UnlockElementNodes(rElement);

        KRATOS_CATCH("")
    }

    void CentralDifferencingWithBasis(
         Element&rElement,
         LocalSystemMatrixType& rElementMatrixDerivative,
         const std::size_t basis_j,
         const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        // Lock element nodes for OMP parallelism
        this->LockElementNodes(rElement);

        std::size_t matrix_size = rElementMatrixDerivative.size1();

        // Positive perturbation
        LocalSystemMatrixType element_matrix_p_perturbed;
        element_matrix_p_perturbed.resize(matrix_size, matrix_size, false);
        this->PerturbElementWithBasis(rElement, 1.0, basis_j, rCurrentProcessInfo);
        rElement.CalculateLeftHandSide(element_matrix_p_perturbed, rCurrentProcessInfo);

        // Negative perturbation
        LocalSystemMatrixType element_matrix_n_perturbed;
        element_matrix_n_perturbed.resize(matrix_size, matrix_size, false);
        this->PerturbElementWithBasis(rElement, -2.0, basis_j, rCurrentProcessInfo);
        rElement.CalculateLeftHandSide(element_matrix_n_perturbed, rCurrentProcessInfo);

        // Reset perturbation
        this->PerturbElementWithBasis(rElement, 1.0, basis_j, rCurrentProcessInfo);

        // Compute LHS derivative
        noalias(rElementMatrixDerivative) = (element_matrix_p_perturbed - element_matrix_n_perturbed) / (2.0 * mFiniteDifferenceStepSize);

        this->UnlockElementNodes(rElement);

        KRATOS_CATCH("")
    }

    // This function perturbs the element with the vector with given eigenvector index basis_j
    void PerturbElementWithBasis(
        Element& rElement,
        const double Step,
        const std::size_t basis_j,
        const ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY

        // Initialize is necessary for updating the section properties of shell elements
        rElement.InitializeNonLinearIteration(rCurrentProcessInfo);
        std::size_t disp_shifter = 3;

        // Loop over element nodes
        for (auto& node_i : rElement.GetGeometry()) {
            auto& node_i_dofs = node_i.GetDofs();
            // Loop over nodal DOFs
            auto it_dof_i = node_i_dofs.begin();
            const int node_i_num_dofs = node_i_dofs.size();
            for (std::size_t dof_idx = 0; dof_idx < node_i_dofs.size(); dof_idx++){

                // Compute and assign the perturbation
                const double dof_perturbation = Step*mFiniteDifferenceStepSize*node_i.GetValue(ROM_BASIS)(dof_idx, basis_j);

                // Some elements need solution step value while others need the current coordinate.
                // Thus perturb all
                // Update solution step value
                (*it_dof_i)->GetSolutionStepValue() += dof_perturbation;

                // Update current nodal coordinates
                if (node_i_num_dofs > 3 && dof_idx > 2) // disp dofs when rots included
                    node_i.Coordinates()[dof_idx-disp_shifter] += dof_perturbation;
                else if (node_i_num_dofs < 4 ) // disp dofs when rots excluded
                    node_i.Coordinates()[dof_idx] += dof_perturbation;
                
                // Increment the dof iterator
                ++it_dof_i;
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

    bool mSchemeIsInitialized;

    bool mDerivativeMatrixType;

    Variable<double> mDerivativeParameter;

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


