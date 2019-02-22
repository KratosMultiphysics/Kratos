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
#if !defined(KRATOS_CONTACT_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER )
#define  KRATOS_CONTACT_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER

/* System includes */
#include <unordered_set>
#include <unordered_map>

/* External includes */

/* Project includes */
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"

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
 * @class ContactResidualBasedEliminationBuilderAndSolver
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
class ContactResidualBasedEliminationBuilderAndSolver
    : public ResidualBasedEliminationBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    ///@name Type Definitions
    ///@{
    
    /// Pointer definition of ContactResidualBasedEliminationBuilderAndSolver
    KRATOS_CLASS_POINTER_DEFINITION(ContactResidualBasedEliminationBuilderAndSolver);

    /// Definitions dependent of the base class
    typedef ResidualBasedEliminationBuilderAndSolver< TSparseSpace, TDenseSpace, TLinearSolver > BaseType;

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
    typedef std::unordered_set<IndexType> IndexSetType;

    ///@}
    ///@name Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor.
     */
    ContactResidualBasedEliminationBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver)
        : BaseType(pNewLinearSystemSolver)
    {
    }

    /** Destructor.
     */
    ~ContactResidualBasedEliminationBuilderAndSolver() override
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
        /**
         * Idem to the not contact version, except that if we fix the displacement in one slave node we should fix the corresponding LM for consistency
         */

        // We create a set of dofs of the displacement slave dofs with LM associated
        std::unordered_map<IndexType, IndexSetType> set_nodes_with_lm_associated;
        if (rModelPart.HasSubModelPart("Contact"))
            set_nodes_with_lm_associated.reserve(rModelPart.GetSubModelPart("Contact").NumberOfNodes());
        // Allocating auxiliar parameters
        IndexType node_id;
        // We start the dof loop
        for (auto& i_dof : BaseType::mDofSet) {
            node_id = i_dof.Id();
            if (IsLMDof(i_dof))
                set_nodes_with_lm_associated.insert({node_id, IndexSetType({})});
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
                const auto& r_variable = i_dof.GetVariable();
                auto& aux_set = (it->second);
                if (i_dof.IsFixed()) {
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
        for (auto& i_dof : BaseType::mDofSet) {
            if (i_dof.IsFree()) {
                node_id = i_dof.Id();
                auto it = set_nodes_with_lm_associated.find(node_id);
                if (it != set_nodes_with_lm_associated.end()) {
                    auto& aux_set = it->second;
                    if (aux_set.find((i_dof.GetVariable()).Key()) != aux_set.end()) {
                        i_dof.FixDof();
                    }
                }
            }
        }

        BaseType::SetUpSystem(rModelPart);
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

}; /* Class ContactResidualBasedEliminationBuilderAndSolver */

///@}

///@name Type Definitions */
///@{


///@}

} /* namespace Kratos.*/

#endif /* KRATOS_CONTACT_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER  defined */
