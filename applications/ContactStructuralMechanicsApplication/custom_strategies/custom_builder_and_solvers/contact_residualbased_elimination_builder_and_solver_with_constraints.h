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
#if !defined(KRATOS_CONTACT_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER_WITH_CONSTRAINTS )
#define  KRATOS_CONTACT_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER_WITH_CONSTRAINTS

/* System includes */
#include <unordered_set>
#include <unordered_map>

/* External includes */

/* Project includes */
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints.h"

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
 * @class ContactResidualBasedEliminationBuilderAndSolverWithConstraints
 * @ingroup ContactStructuralMechanicsApplication
 * @brief Current class provides an implementation for contact builder and solving operations. (elimination)
 * @details The RHS is constituted by the unbalanced loads (residual). Degrees of freedom are reordered putting the restrained degrees of freedom at the end of the system ordered in reverse order with respect to the DofSet and not considered the inactive ones. Imposition of the dirichlet conditions is naturally dealt with as the residual already contains this information. Calculation of the reactions involves a cost very similiar to the calculation of the total residual
 * @author Vicente Mataix Ferrandiz
 * @tparam TSparseSpace The sparse matrix system considered
 * @tparam TDenseSpace The dense matrix system
 * @tparam TLinearSolver The type of linear solver considered
 */
template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class ContactResidualBasedEliminationBuilderAndSolverWithConstraints
    : public ResidualBasedEliminationBuilderAndSolverWithConstraints< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    ///@name Type Definitions
    ///@{
    
    /// Pointer definition of ContactResidualBasedEliminationBuilderAndSolverWithConstraints
    KRATOS_CLASS_POINTER_DEFINITION(ContactResidualBasedEliminationBuilderAndSolverWithConstraints);

    /// Definitions dependent of the base class
    typedef ResidualBasedEliminationBuilderAndSolverWithConstraints< TSparseSpace, TDenseSpace, TLinearSolver > BaseType;

    /// Base types definitions
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::TDataType TDataType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    /// The node type
    typedef Node<3> NodeType;

    /// The definition of the dof type
    typedef typename ModelPart::DofType DofType;

    /// The size type
    typedef std::size_t SizeType;

    /// The index type
    typedef std::size_t IndexType;

    /// Index set definition
    typedef std::unordered_set<IndexType> SetIndexType;

    ///@}
    ///@name Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor.
     */
    ContactResidualBasedEliminationBuilderAndSolverWithConstraints(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BaseType(pNewLinearSystemSolver)
    {
    }

    /** Destructor.
     */
    ~ContactResidualBasedEliminationBuilderAndSolverWithConstraints() override
    {
    }

    ///@}
    ///@name Operators
    ///@{
    
    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief It organises the dofset in order to speed up the building phase
     * @param rModelPart The model part to compute
     */
    void SetUpSystem(
        ModelPart& rModelPart
        ) override
    {
        if(rModelPart.MasterSlaveConstraints().size() > 0)
            SetUpSystemWithConstraints(rModelPart);
        else
            BaseSetUpSystem(rModelPart);
    }

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
    
    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method computes the equivalent coounter part of the SetUpSystem when using constraints
     * @param rModelPart The model part of the problem to solve
     */
    void SetUpSystemWithConstraints(ModelPart& rModelPart)
    {
        KRATOS_TRY

        // First we set up the system of equations without constraints
        BaseSetUpSystem(rModelPart);

        // Add the computation of the global ids of the solvable dofs
        IndexType counter = 0;
        for (auto& dof : BaseType::mDofSet) {
            if (dof.EquationId() < BaseType::mEquationSystemSize) {
                auto it = BaseType::mDoFSlaveSet.find(dof);
                if (it == BaseType::mDoFSlaveSet.end()) {
                    ++counter;
                }
            }
        }

        // The total system of equations to be solved
        BaseType::mDoFToSolveSystemSize = counter;

        KRATOS_CATCH("ContactResidualBasedEliminationBuilderAndSolverWithConstraints::FormulateGlobalMasterSlaveRelations failed ..");
    }

    /**
     * @brief It organises the dofset in order to speed up the building phase (base one)
     * @param rModelPart The model part to compute
     */
    void BaseSetUpSystem(ModelPart& rModelPart)
    {
        /**
         * Idem to the not contact version, except that if we fix the displacement in one slave node we should fix the corresponding LM for consistency
         */

        // We create a set of dofs of the displacement slave dofs with LM associated
        std::unordered_map<IndexType, SetIndexType> set_nodes_with_lm_associated;
        if (rModelPart.HasSubModelPart("Contact"))
            set_nodes_with_lm_associated.reserve(rModelPart.GetSubModelPart("Contact").NumberOfNodes());
        // Allocating auxiliar parameters
        IndexType node_id;
        // We start the dof loop
        for (auto& i_dof : BaseType::mDofSet) {
            node_id = i_dof.Id();
            if (IsLMDof(i_dof))
                set_nodes_with_lm_associated.insert({node_id, SetIndexType({})});
        }

        // Auxiliar keys
        const IndexType key_lm_x = VECTOR_LAGRANGE_MULTIPLIER_X.Key();
        const IndexType key_lm_y = VECTOR_LAGRANGE_MULTIPLIER_Y.Key();
        const IndexType key_lm_z = VECTOR_LAGRANGE_MULTIPLIER_Z.Key();

        // We see which LM block
        for (auto& i_dof : BaseType::mDofSet) {
            node_id = i_dof.Id();
            auto it = set_nodes_with_lm_associated.find(node_id);
            if ( it != set_nodes_with_lm_associated.end()) {
                if (i_dof.IsFixed()) {
                    const auto& r_variable = i_dof.GetVariable();
                    auto& aux_set = (it->second);
                    if (r_variable == DISPLACEMENT_X) {
                        aux_set.insert(key_lm_x);
                    } else if (r_variable == DISPLACEMENT_Y) {
                        aux_set.insert(key_lm_y);
                    } else if (r_variable == DISPLACEMENT_Z) {
                        aux_set.insert(key_lm_z);
                    }
                }
            }
        }

        // We do now the loop over the dofs
        int free_id = 0;
        int fix_id = BaseType::mDofSet.size();

        for (auto& i_dof : BaseType::mDofSet) {
            if (i_dof.IsFixed())
                i_dof.SetEquationId(--fix_id);
            else {
                node_id = i_dof.Id();
                auto it = set_nodes_with_lm_associated.find(node_id);
                if (it != set_nodes_with_lm_associated.end()) {
                    auto& aux_set = it->second;
                    if (aux_set.find((i_dof.GetVariable()).Key()) != aux_set.end()) {
                        i_dof.SetEquationId(--fix_id);
                    } else {
                        i_dof.SetEquationId(free_id++);
                    }
                } else {
                    i_dof.SetEquationId(free_id++);
                }
            }
        }

        BaseType::mEquationSystemSize = fix_id;
    }

    /**
     * @brief Checks if the degree of freedom belongs to a displacement DoF
     * @param rDoF The degree of freedom
     * @return True if the DoF corresponds with a displacement dof
     */
    static inline bool IsDisplacementDof(const DofType& rDoF)
    {
        const auto& r_variable = rDoF.GetVariable();
        if (r_variable == DISPLACEMENT_X ||
            r_variable == DISPLACEMENT_Y ||
            r_variable == DISPLACEMENT_Z) {
                return true;
        }

        return false;
    }

    /**
     * @brief Checks if the degree of freedom belongs to a LM DoF
     * @param rDoF The degree of freedom
     * @return True if the DoF corresponds with a LM dof
     */
    static inline bool IsLMDof(const DofType& rDoF)
    {
        const auto& r_variable = rDoF.GetVariable();
        if (r_variable == VECTOR_LAGRANGE_MULTIPLIER_X ||
            r_variable == VECTOR_LAGRANGE_MULTIPLIER_Y ||
            r_variable == VECTOR_LAGRANGE_MULTIPLIER_Z) {
                return true;
        }

        return false;
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

}; /* Class ContactResidualBasedEliminationBuilderAndSolverWithConstraints */

///@}

///@name Type Definitions */
///@{


///@}

} /* namespace Kratos.*/

#endif /* KRATOS_CONTACT_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER_WITH_CONSTRAINTS  defined */
