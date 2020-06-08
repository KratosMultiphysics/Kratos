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

        if (rCurrentProcessInfo[BUILD_LEVEL] == 1)
        {   
            // Stiffness matrix contribution   
            rElement.CalculateLeftHandSide(rLHS_Contribution, rCurrentProcessInfo);
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
        
        Matrix element_LHS_derivative;
        
        int eigenvalue_i = rCurrentProcessInfo[EIGENVALUE_I];
        int eigenvalue_j = rCurrentProcessInfo[EIGENVALUE_J];
        // KRATOS_WATCH(eigenvalue_i)
        // KRATOS_WATCH(eigenvalue_j)
        //const double eigenvalue = rCurrentProcessInfo[EIGENVALUE_VECTOR](eigenvalue_i); // This will be used for dynamic derivatives

        // Get PhiElemental
        Vector PhiElemental;
        std::size_t num_element_dofs = 0;
        for (auto& node_i : rElement.GetGeometry())
            num_element_dofs += node_i.GetDofs().size();
        PhiElemental.resize(num_element_dofs);
        
        std::size_t dof_ctr = 0;
        for (auto& node_i : rElement.GetGeometry()) {
            auto& node_i_dofs = node_i.GetDofs();
            
            const Matrix *pPhiNodal = &node_i.GetValue(ROM_BASIS);
            auto dof_i = node_i_dofs.begin();
            for (std::size_t dof_idx = 0; dof_idx < node_i_dofs.size(); dof_idx++){
                dof_i = dof_i + dof_idx;

                PhiElemental[dof_ctr + dof_idx] = (*pPhiNodal)(dof_idx, eigenvalue_i);
            }
            dof_ctr += node_i_dofs.size();
        }
        
        rRHS_Contribution.clear();
        rRHS_Contribution.resize(0);

        // Build RHS contributions
        if (rCurrentProcessInfo[BUILD_LEVEL] == 1)
        {   
            // Perturb each nodal DOF
            // Loop over element nodes
            for (auto& node_i : rElement.GetGeometry()) {
                auto& node_i_dofs = node_i.GetDofs();
                // Loop over nodal DOFs
                auto dof_i = node_i_dofs.begin();
                for (std::size_t dof_idx = 0; dof_idx < node_i_dofs.size(); dof_idx++){
                    const double perturbationMag = node_i.GetValue(ROM_BASIS)(dof_idx, eigenvalue_j)*1e-3;
                    // std::cout << "EID: " << rElement.Id() << " NID: " << node_i.Id() << " DOF: " << (*dof_i)->GetVariable().Name() << std::endl;
                    if (abs(perturbationMag) > 0.0)
                    {
                        RomFiniteDifferenceUtility::CalculateLeftHandSideDOFDerivative(rElement,
                                                                            *(*dof_i),
                                                                            perturbationMag,
                                                                            element_LHS_derivative,
                                                                            rCurrentProcessInfo);
                        
                        if (rRHS_Contribution.size() == 0){
                            // std::cout << "Resizing rRHS_Contribution" << std::endl;
                            rRHS_Contribution.resize(element_LHS_derivative.size1());
                            for (std::size_t iRHS = 0; iRHS < rRHS_Contribution.size(); iRHS++)
                                rRHS_Contribution[iRHS] = 0.0;
                        }

                        // KRATOS_WATCH(element_LHS_derivative)
                        
                        if (element_LHS_derivative.size1() != PhiElemental.size()){
                            // retrieve only relevant dofs from PhiElemental
                            Vector tmpPhiElemental;
                            tmpPhiElemental.resize(PhiElemental.size()/2);
                            for (std::size_t iNode = 0; iNode < rElement.GetGeometry().size(); iNode++){
                                for (std::size_t iXYZ = 0; iXYZ < 3; iXYZ++){
                                    tmpPhiElemental[iNode*3+iXYZ] = PhiElemental[iNode*6+3+iXYZ];
                                }                                
                            }
                            // KRATOS_WATCH(tmpPhiElemental)
                            rRHS_Contribution -= prod(element_LHS_derivative, tmpPhiElemental);
                        }
                        else{
                            // KRATOS_WATCH(PhiElemental)
                            rRHS_Contribution -= prod(element_LHS_derivative, PhiElemental);
                        }   
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

        // KRATOS_WATCH(rRHS_Contribution)

        rElement.EquationIdVector(EquationId,rCurrentProcessInfo);

        KRATOS_CATCH("")
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


