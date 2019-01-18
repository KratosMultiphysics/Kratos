//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//
#if !defined(KRATOS_CONTACT_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER )
#define  KRATOS_CONTACT_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER

/* System includes */

/* External includes */

/* Project includes */
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"

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
 * @class ContactResidualBasedBlockBuilderAndSolver
 * @ingroup ContactStructuralMechanicsApplication
 * @brief Current class provides an implementation for contact builder and solving operations.
 * @details The RHS is constituted by the unbalanced loads (residual). Degrees of freedom are reordered putting the restrained degrees of freedom at the end of the system ordered in reverse order with respect to the DofSet. Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
this information. Calculation of the reactions involves a cost very similiar to the calculation of the total residual
 * @author Vicente Mataix Ferrandiz
 * @tparam TSparseSpace The sparse matrix system considered
 * @tparam TDenseSpace The dense matrix system
 * @tparam TLinearSolver The type of linear solver considered
 * @tparam TBuilderAndSolver The builder and solver considered as base
 */
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver, //= LinearSolver<TSparseSpace,TDenseSpace>
         class TBuilderAndSolver = ResidualBasedBlockBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
         >
class ContactResidualBasedBlockBuilderAndSolver
    : public TBuilderAndSolver
{
public:
    ///@name Type Definitions
    ///@{
    
    /// Pointer definition of ContactResidualBasedBlockBuilderAndSolver
    KRATOS_CLASS_POINTER_DEFINITION(ContactResidualBasedBlockBuilderAndSolver);

    /// Definitions dependent of the base class
    typedef TBuilderAndSolver BaseType;

    typedef typename BaseType::TSchemeType TSchemeType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor.
     */
    ContactResidualBasedBlockBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BaseType(pNewLinearSystemSolver)
    {
    }

    /** Destructor.
     */
    ~ContactResidualBasedBlockBuilderAndSolver() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief This method imposses the BC of Dirichlet. It will fill with 0 the corresponding DoF
     * @param pScheme The pointer to the scheme considered
     * @param rModelPart The model part of the problem to solve 
     * @param A The LHS of the system
     * @param Dx The current solution increment
     * @param b The RHS of the system
     */
    void ApplyDirichletConditions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        FixIsolatedNodes(rModelPart);
        
        BaseType::ApplyDirichletConditions(pScheme, rModelPart, A, Dx, b);
        
        FreeIsolatedNodes(rModelPart);
    }
    
    /**
     * @brief This method buils the RHS of the system of equations
     * @param pScheme The pointer to the scheme considered
     * @param rModelPart The model part of the problem to solve 
     * @param b The RHS of the system
     */
    void BuildRHS(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemVectorType& b
        ) override
    {
        FixIsolatedNodes(rModelPart);
        
        BaseType::BuildRHS(pScheme, rModelPart, b);
        
        FreeIsolatedNodes(rModelPart);
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

    /**
     * @brief This method check the ISOLATED nodes and it fixes
     * @param rModelPart The model part to compute
     */
    void FixIsolatedNodes(ModelPart& rModelPart)
    {
        KRATOS_ERROR_IF_NOT(rModelPart.HasSubModelPart("Contact")) << "CONTACT MODEL PART NOT CREATED" << std::endl;
        KRATOS_ERROR_IF_NOT(rModelPart.HasSubModelPart("ComputingContact")) << "CONTACT COMPUTING MODEL PART NOT CREATED" << std::endl;
        ModelPart& contact_model_part = rModelPart.GetSubModelPart("Contact"); 
        ModelPart& computing_contact_model_part = rModelPart.GetSubModelPart("ComputingContact"); 
        
        // We reset the flag
        auto& nodes_array = contact_model_part.Nodes();
        #pragma omp parallel for 
        for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
            (nodes_array.begin() + i)->Set(VISITED, false);
            (nodes_array.begin() + i)->Set(ISOLATED, false);
        }
        
        // Now we set the flag in the nodes
        auto& conditions_array = computing_contact_model_part.Conditions();
        
        #pragma omp parallel for 
        for(int i = 0; i < static_cast<int>(conditions_array.size()); ++i) {
            auto it_cond = conditions_array.begin() + i;
            auto& geom = it_cond->GetGeometry();
            for (std::size_t i_node = 0; i_node < geom.size(); ++i_node) {
                geom[i_node].SetLock();
                if (geom[i_node].Is(VISITED) == false) {
                    geom[i_node].Set(ISOLATED, it_cond->Is(ISOLATED));
                    geom[i_node].Set(VISITED, true);
                } else {
                    geom[i_node].Set(ISOLATED, geom[i_node].Is(ISOLATED) && it_cond->Is(ISOLATED));
                }
                geom[i_node].UnSetLock();
            }
        }
        
        // We fix the LM
        #pragma omp parallel for 
        for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
            auto it_node = nodes_array.begin() + i;
            if (it_node->Is(ISOLATED) == true) {
                if (it_node->SolutionStepsDataHas(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE))
                    it_node->Fix(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE);
                else if (it_node->SolutionStepsDataHas(VECTOR_LAGRANGE_MULTIPLIER_X)) {
                    it_node->Fix(VECTOR_LAGRANGE_MULTIPLIER_X);
                    it_node->Fix(VECTOR_LAGRANGE_MULTIPLIER_Y);
                    it_node->Fix(VECTOR_LAGRANGE_MULTIPLIER_Z);
                }
            }
        }
    }
    
    /**
     * @brief This method releases the ISOLATED nodes 
     * @param rModelPart The model part to compute
     */
    void FreeIsolatedNodes(ModelPart& rModelPart)
    {
        KRATOS_ERROR_IF_NOT(rModelPart.HasSubModelPart("Contact")) << "CONTACT MODEL PART NOT CREATED" << std::endl;
        ModelPart& contact_model_part = rModelPart.GetSubModelPart("Contact");
        
        // We release the LM
        auto& nodes_array = contact_model_part.Nodes();
        #pragma omp parallel for 
        for(int i = 0; i < static_cast<int>(nodes_array.size()); ++i) {
            auto it_node = nodes_array.begin() + i;
            if (it_node->Is(ISOLATED) == true) {
                if (it_node->SolutionStepsDataHas(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE))
                    it_node->Free(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE);
                else if (it_node->SolutionStepsDataHas(VECTOR_LAGRANGE_MULTIPLIER_X)) {
                    it_node->Free(VECTOR_LAGRANGE_MULTIPLIER_X);
                    it_node->Free(VECTOR_LAGRANGE_MULTIPLIER_Y);
                    it_node->Free(VECTOR_LAGRANGE_MULTIPLIER_Z);
                }
            }
        }
    }
    
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

}; /* Class ContactResidualBasedBlockBuilderAndSolver */

///@}

///@name Type Definitions */
///@{


///@}

} /* namespace Kratos.*/

#endif /* KRATOS_CONTACT_RESIDUAL_BASED_BLOCK_BUILDER_AND_SOLVER  defined */
