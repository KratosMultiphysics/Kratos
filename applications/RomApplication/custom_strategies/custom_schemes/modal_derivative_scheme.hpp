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
    ModalDerivativeScheme() : Scheme<TSparseSpace,TDenseSpace>() 
    {
        // KRATOS_WATCH("ModalDerivativeScheme::ModalDerivativeScheme")
    }

    /// Destructor.
    ~ModalDerivativeScheme() override 
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief This is the place to initialize the Scheme.
     * @details This is intended to be called just once when the strategy is initialized
     * @param rModelPart The model part of the problem to solve
     */
    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY
        // KRATOS_WATCH("ModalDerivativeScheme::Initialize")
        mSchemeIsInitialized = true;
        KRATOS_CATCH("")
    }

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

        // KRATOS_WATCH("ModalDerivativeScheme::InitializeNonLinIteration")
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

        // KRATOS_WATCH("ModalDerivativeScheme::InitializeNonLinearIteration Element")

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

        // KRATOS_WATCH("ModalDerivativeScheme::InitializeNonLinearIteration Condition")

        (pCurrentCondition)->InitializeNonLinearIteration(rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    /**
     * @brief Function to be called when it is needed to finalize an iteration. It is designed to be called at the end of each non linear iteration
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    void FinalizeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        KRATOS_TRY

        // KRATOS_WATCH("ModalDerivativeScheme::FinalizeNonLinIteration")

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();

        // Definition of the first element iterator
        const auto it_elem_begin = rModelPart.ElementsBegin();

        // Finalizes non-linear iteration for all of the elements
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Elements().size()); ++i) {
            auto it_elem = it_elem_begin + i;
            it_elem->FinalizeNonLinearIteration(r_current_process_info);
        }

        // Definition of the first condition iterator
        const auto it_cond_begin = rModelPart.ConditionsBegin();

        // Finalizes non-linear iteration  for all of the conditions
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.Conditions().size()); ++i) {
            auto it_cond = it_cond_begin + i;
            it_cond->FinalizeNonLinearIteration(r_current_process_info);
        }

        // Definition of the first constraint iterator
        const auto it_const_begin = rModelPart.MasterSlaveConstraintsBegin();

        // Finalizes non-linear iteration for all of the constraints
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rModelPart.MasterSlaveConstraints().size()); ++i) {
            auto it_const = it_const_begin + i;
            it_const->FinalizeNonLinearIteration(r_current_process_info);
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Performing the prediction of the solution.
     * @warning Must be defined in derived classes
     * @param rModelPart The model part of the problem to solve
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    void Predict(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    /**
     * @brief Performing the update of the solution.
     * @warning Must be defined in derived classes
     * @param rModelPart The model part of the problem to solve
     * @param rDofSet Set of all primary variables
     * @param A LHS matrix
     * @param Dx Incremental update of primary variables
     * @param b RHS Vector
     */
    void Update(
        ModelPart& rModelPart,
        DofsArrayType& rDofSet,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        KRATOS_TRY
        KRATOS_CATCH("")
    }

    // Element contributions
    void CalculateSystemContributions(
        Element::Pointer pCurrentElement,
        LocalSystemMatrixType& rLHS_Contribution,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        // KRATOS_WATCH("ModalDerivativeScheme::CalculateSystemContributions")

        this->Calculate_LHS_Contribution(pCurrentElement, rLHS_Contribution, rEquationId, rCurrentProcessInfo);

        this->Calculate_RHS_Contribution(pCurrentElement, rRHS_Contribution, rEquationId, rCurrentProcessInfo);

        pCurrentElement->EquationIdVector(rEquationId,rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void Calculate_LHS_Contribution(
        Element::Pointer pCurrentElement,
        LocalSystemMatrixType& rLHS_Contribution,
        Element::EquationIdVectorType& rEquationId,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        // KRATOS_WATCH("ModalDerivativeScheme::Calculate_LHS_Contribution Element")

        if (rCurrentProcessInfo[BUILD_LEVEL] == 1)
        {   
            // Stiffness matrix contribution   
            pCurrentElement->CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);
        } 
        else if (rCurrentProcessInfo[BUILD_LEVEL] == 2) 
        {   
            // Mass matrix contribution is going to be implemented here
            KRATOS_ERROR <<"Invalid BUILD_LEVEL: " << rCurrentProcessInfo[BUILD_LEVEL] << "\nDynamic derivatives not implemented yet!" << std::endl;
        } 
        else 
        {
            KRATOS_ERROR <<"Invalid BUILD_LEVEL: " << rCurrentProcessInfo[BUILD_LEVEL] << std::endl;
        }

        pCurrentElement->EquationIdVector(rEquationId, rCurrentProcessInfo);

        KRATOS_CATCH("")
    }

    void Calculate_RHS_Contribution(
        Element::Pointer pCurrentElement,
        LocalSystemVectorType& rRHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        KRATOS_TRY

        // KRATOS_WATCH("ModalDerivativeScheme::Calculate_RHS_Contribution Element")
        // KRATOS_WATCH(pCurrentElement->Id())
        
        Matrix element_LHS_derivative;
        int eigenvalue_i = rCurrentProcessInfo[EIGENVALUE_I];
        int eigenvalue_j = rCurrentProcessInfo[EIGENVALUE_J];
        //const double eigenvalue = rCurrentProcessInfo[EIGENVALUE_VECTOR](eigenvalue_i); // This will be used for dynamic derivatives

        // KRATOS_WATCH(eigenvalue_i)
        // KRATOS_WATCH(eigenvalue_j)

        // Get PhiElemental
        Vector PhiElemental;
        int num_element_dofs = 0;
        for (auto& node_i : pCurrentElement->GetGeometry())
            num_element_dofs += node_i.GetDofs().size();
        PhiElemental.resize(num_element_dofs);
        // element_LHS_derivative.resize(num_element_dofs, num_element_dofs);
        
        unsigned int dof_ctr = 0;
        for (auto& node_i : pCurrentElement->GetGeometry()) {
            auto& node_i_dofs = node_i.GetDofs();
            
            const Matrix *pPhiNodal = &node_i.GetValue(ROM_BASIS);
            auto dof_i = node_i_dofs.begin();
            for (unsigned int dof_idx = 0; dof_idx < node_i_dofs.size(); dof_idx++){
                dof_i = dof_i + dof_idx;

                PhiElemental[dof_ctr + dof_idx] = (*pPhiNodal)(dof_idx, eigenvalue_i);
            }
            dof_ctr += node_i_dofs.size();
        }

        // KRATOS_WATCH(PhiElemental)

        rRHS_Contribution.resize(num_element_dofs);
        for (int i = 0; i < rRHS_Contribution.size(); i++)
            rRHS_Contribution[i] = 0.0;
        // Build RHS contributions
        if (rCurrentProcessInfo[BUILD_LEVEL] == 1)
        {   
            // Perturb each nodal DOF
            // Loop over element nodes
            for (auto& node_i : pCurrentElement->GetGeometry()) {
                auto& node_i_dofs = node_i.GetDofs();

                // std::cout << node_i.Id() << "\n =============================== " << std::endl;
                
                // Loop over nodal DOFs
                auto dof_i = node_i_dofs.begin();
                for (unsigned int dof_idx = 0; dof_idx < node_i_dofs.size(); dof_idx++){
                    // std::cout << dof_idx << "\n ------------------------------- " << std::endl;
                    const double perturbationMag = node_i.GetValue(ROM_BASIS)(dof_idx, eigenvalue_j);
                    if (abs(perturbationMag) > 0.0)
                    {
                        RomFiniteDifferenceUtility::CalculateLeftHandSideDOFDerivative(*pCurrentElement,
                                                                            *(*dof_i),
                                                                            perturbationMag,
                                                                            element_LHS_derivative,
                                                                            rCurrentProcessInfo);
                        // KRATOS_WATCH(element_LHS_derivative)

                        // KRATOS_WATCH(rRHS_Contribution)
                        if (rRHS_Contribution.size() != element_LHS_derivative.size1()){
                            rRHS_Contribution.resize(element_LHS_derivative.size1());
                            for (int i = 0; i < rRHS_Contribution.size(); i++)
                                rRHS_Contribution[i] = 0.0;
                        }   
                        // KRATOS_WATCH(rRHS_Contribution)
                        if (element_LHS_derivative.size1() != PhiElemental.size()){
                            // retrieve only relevant dofs from PhiElemental
                            Vector tmpPhiElemental;
                            tmpPhiElemental.resize(PhiElemental.size()/2);
                            for (int iNode = 0; iNode < pCurrentElement->GetGeometry().size(); iNode++){
                                for (int iXYZ = 0; iXYZ < 3; iXYZ++){
                                    tmpPhiElemental[iNode*3+iXYZ] = PhiElemental[iNode*6+3+iXYZ];
                                }                                
                            }
                            // KRATOS_WATCH(tmpPhiElemental)
                            rRHS_Contribution += prod(element_LHS_derivative, tmpPhiElemental);
                        }
                        else {
                            // KRATOS_WATCH(PhiElemental)
                            rRHS_Contribution += prod(element_LHS_derivative, PhiElemental);
                        }
                            
                        // KRATOS_WATCH(rRHS_Contribution)
                    }
                    dof_i++;
                }
            }
        }
        else if (rCurrentProcessInfo[BUILD_LEVEL] == 2) 
        {   
            // Mass matrix contribution is going to be implemented here
            KRATOS_ERROR <<"Invalid BUILD_LEVEL: " << rCurrentProcessInfo[BUILD_LEVEL] << "\nDynamic derivatives not implemented yet!" << std::endl;
        } 
        else 
        {
            KRATOS_ERROR <<"Invalid BUILD_LEVEL: " << rCurrentProcessInfo[BUILD_LEVEL] << std::endl;
        }

        pCurrentElement->EquationIdVector(EquationId,rCurrentProcessInfo);

        // // HACK for DEBUG
        // if (pCurrentElement->Id() == 1)
        //     std::exit(EXIT_FAILURE);

        KRATOS_CATCH("")
    }

    // Condition contributions
    void Condition_CalculateSystemContributions(
        Condition::Pointer pCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        // pCurrentCondition->CalculateLocalSystem(LHS_Contribution, RHS_Contribution, rCurrentProcessInfo);
    }

    void Condition_Calculate_LHS_Contribution(
        Condition::Pointer pCurrentCondition,
        LocalSystemMatrixType& LHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        // pCurrentCondition->CalculateLeftHandSide(LHS_Contribution, rCurrentProcessInfo);
    }

    void Condition_Calculate_RHS_Contribution(
        Condition::Pointer pCurrentCondition,
        LocalSystemVectorType& RHS_Contribution,
        Element::EquationIdVectorType& EquationId,
        ProcessInfo& rCurrentProcessInfo
        ) override
    {
        // pCurrentCondition->CalculateRightHandSide(RHS_Contribution, rCurrentProcessInfo);
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


