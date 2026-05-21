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

// External includes

// Project includes
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

    /// Builder and solver base class
    using BaseBuilderAndSolverType = BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver>;

    /// Definitions dependent on the base class
    using BaseType = TBuilderAndSolver;

    /// The definition of the current class
    using ClassType = ContactResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver, TBuilderAndSolver>;

    /// Type alias for the scheme used in the base class
    using TSchemeType = typename BaseType::TSchemeType;

    /// Type alias for the data type used in the base class
    using TDataType = typename BaseType::TDataType;

    /// Type alias for the collection of degrees of freedom used in the base class
    using DofsArrayType = typename BaseType::DofsArrayType;

    /// Type alias for the system matrix used in the base class
    using TSystemMatrixType = typename BaseType::TSystemMatrixType;

    /// Type alias for the system vector used in the base class
    using TSystemVectorType = typename BaseType::TSystemVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    explicit ContactResidualBasedBlockBuilderAndSolver() : BaseType()
    {
    }

    /**
     * @brief Default constructor. (with parameters)
     */
    explicit ContactResidualBasedBlockBuilderAndSolver(
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

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"(
        {
            "name" : "contact_block_builder_and_solver"
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
        return "contact_block_builder_and_solver";
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

    /**
     * @brief This method check the ISOLATED nodes and it fixes
     * @param rModelPart The model part to compute
     */
    void FixIsolatedNodes(ModelPart& rModelPart)
    {
        KRATOS_ERROR_IF_NOT(rModelPart.HasSubModelPart("Contact")) << "CONTACT MODEL PART NOT CREATED" << std::endl;
        KRATOS_ERROR_IF_NOT(rModelPart.HasSubModelPart("ComputingContact")) << "CONTACT COMPUTING MODEL PART NOT CREATED" << std::endl;
        ModelPart& r_contact_model_part = rModelPart.GetSubModelPart("Contact");
        ModelPart& r_computing_contact_model_part = rModelPart.GetSubModelPart("ComputingContact");

        // We reset the flag
        auto& r_nodes_array = r_contact_model_part.Nodes();
        block_for_each(r_nodes_array, [&](Node& rNode) {
            rNode.Set(VISITED, false);
            rNode.Set(ISOLATED, false);
        });

        // Now we set the flag in the nodes
        auto& r_conditions_array = r_computing_contact_model_part.Conditions();
        block_for_each(r_conditions_array, [&](Condition& rCond) {
            auto& r_parent_geometry = rCond.GetGeometry().GetGeometryPart(0);
            for (std::size_t i_node = 0; i_node < r_parent_geometry.size(); ++i_node) {
                r_parent_geometry[i_node].SetLock();
                if (r_parent_geometry[i_node].Is(VISITED) == false) {
                    r_parent_geometry[i_node].Set(ISOLATED, rCond.Is(ISOLATED));
                    r_parent_geometry[i_node].Set(VISITED, true);
                } else {
                    r_parent_geometry[i_node].Set(ISOLATED, r_parent_geometry[i_node].Is(ISOLATED) && rCond.Is(ISOLATED));
                }
                r_parent_geometry[i_node].UnSetLock();
            }
        });

        // We fix the LM
        block_for_each(r_nodes_array, [&](Node& rNode) {
            if (rNode.Is(ISOLATED)) {
                if (rNode.SolutionStepsDataHas(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE)) {
                    rNode.Fix(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE);
                } else if (rNode.SolutionStepsDataHas(VECTOR_LAGRANGE_MULTIPLIER_X)) {
                    rNode.Fix(VECTOR_LAGRANGE_MULTIPLIER_X);
                    rNode.Fix(VECTOR_LAGRANGE_MULTIPLIER_Y);
                    rNode.Fix(VECTOR_LAGRANGE_MULTIPLIER_Z);
                }
            }
        });
    }

    /**
     * @brief This method releases the ISOLATED nodes
     * @param rModelPart The model part to compute
     */
    void FreeIsolatedNodes(ModelPart& rModelPart)
    {
        KRATOS_ERROR_IF_NOT(rModelPart.HasSubModelPart("Contact")) << "CONTACT MODEL PART NOT CREATED" << std::endl;
        ModelPart& r_contact_model_part = rModelPart.GetSubModelPart("Contact");

        // We release the LM
        auto& r_nodes_array = r_contact_model_part.Nodes();
        block_for_each(r_nodes_array, [&](Node& rNode) {
            if (rNode.Is(ISOLATED)) {
                if (rNode.SolutionStepsDataHas(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE)) {
                    rNode.Free(LAGRANGE_MULTIPLIER_CONTACT_PRESSURE);
                } else if (rNode.SolutionStepsDataHas(VECTOR_LAGRANGE_MULTIPLIER_X)) {
                    rNode.Free(VECTOR_LAGRANGE_MULTIPLIER_X);
                    rNode.Free(VECTOR_LAGRANGE_MULTIPLIER_Y);
                    rNode.Free(VECTOR_LAGRANGE_MULTIPLIER_Z);
                }
            }
        });
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
