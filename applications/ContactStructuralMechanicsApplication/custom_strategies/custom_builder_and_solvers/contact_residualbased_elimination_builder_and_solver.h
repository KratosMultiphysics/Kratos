// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

#pragma once

// System includes
#include <unordered_set>
#include <unordered_map>

// External includes

// Project includes
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
 * @details The RHS is constituted by the unbalanced loads (residual). Degrees of freedom are reordered putting the restrained degrees of freedom at the end of the system ordered in reverse order with respect to the DofSet and not considered the inactive ones. Imposition of the dirichlet conditions is naturally dealt with as the residual already contains this information. Calculation of the reactions involves a cost very similar to the calculation of the total residual
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

    /// Builder and solver base class
    using BaseBuilderAndSolverType = BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>;

    /// Definitions dependent on the base class
    using BaseType = ResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>;

    /// The definition of the current class
    using ClassType = ContactResidualBasedEliminationBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>;

    /// Base types definitions
    using TSchemeType = typename BaseType::TSchemeType;
    using TDataType = typename BaseType::TDataType;
    using DofsArrayType = typename BaseType::DofsArrayType;
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;
    using TSystemVectorType = typename BaseType::TSystemVectorType;

    /// The definition of the dof type
    using DofType = typename ModelPart::DofType;

    /// The size type
    using SizeType = std::size_t;

    /// The index type
    using IndexType= std::size_t;

    /// Index set definition
    using IndexSetType = std::unordered_set<IndexType>;

    ///@}
    ///@name Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{


    /**
     * @brief Default constructor
     */
    explicit ContactResidualBasedEliminationBuilderAndSolver() : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     */
    explicit ContactResidualBasedEliminationBuilderAndSolver(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters
        ) : BaseType(pNewLinearSystemSolver)
    {
        // Validate and assign defaults
        ThisParameters = this->ValidateAndAssignParameters(ThisParameters, this->GetDefaultParameters());
        this->AssignSettings(ThisParameters);
    }

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
     * @brief Create method
     * @param pNewLinearSystemSolver The linear solver for the system of equations
     * @param ThisParameters The configuration parameters
     */
    typename BaseBuilderAndSolverType::Pointer Create(
        typename TLinearSolver::Pointer pNewLinearSystemSolver,
        Parameters ThisParameters
        ) const override
    {
        return Kratos::make_shared<ClassType>(pNewLinearSystemSolver,ThisParameters);
    }

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
        // Allocating auxiliary parameters
        IndexType node_id;
        // We start the dof loop
        for (auto& i_dof : BaseType::mDofSet) {
            node_id = i_dof.Id();
            if (IsLMDof(i_dof))
                set_nodes_with_lm_associated.insert({node_id, IndexSetType({})});
        }

        // Auxiliary keys
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

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name" : "contact_residual_elimination_builder_and_solver"
        })");

        // Getting base class default parameters
        const Parameters base_default_parameters = BaseType::GetDefaultParameters();
        default_parameters.RecursivelyAddMissingParameters(base_default_parameters);
        return default_parameters;
    }

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "contact_residual_elimination_builder_and_solver";
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

    /**
     * @brief This method assigns settings to member variables
     * @param ThisParameters Parameters that are assigned to the member variables
     */
    void AssignSettings(const Parameters ThisParameters) override
    {
        BaseType::AssignSettings(ThisParameters);
    }

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