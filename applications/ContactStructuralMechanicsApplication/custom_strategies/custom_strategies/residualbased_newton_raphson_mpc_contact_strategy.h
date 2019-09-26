// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_RESIDUALBASED_NEWTON_RAPHSON_MPC_CONTACT_STRATEGY)
#define KRATOS_RESIDUALBASED_NEWTON_RAPHSON_MPC_CONTACT_STRATEGY

/* System Includes */

/* External Includes */

/* Project includes */
#include "contact_structural_mechanics_application_variables.h"
#include "includes/kratos_parameters.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"

// Strategies
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"

// Contact criteria
#include "custom_strategies/custom_convergencecriterias/mpc_contact_criteria.h"

// Utilities
#include "utilities/variable_utils.h"
#include "utilities/color_utilities.h"
#include "utilities/math_utils.h"

// // Processes
// #include "processes/fast_transfer_between_model_parts_process.h"

namespace Kratos {

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
 * @class ResidualBasedNewtonRaphsonMPCContactStrategy
 * @ingroup ContactStructuralMechanicsApplication
 * @brief  Contact Newton Raphson class
 * @details This class is a specialization of the Newton Raphson strategy with some custom modifications for contact problems
 * @author Vicente Mataix Ferrandiz
*/
template<class TSparseSpace,
         class TDenseSpace, // = DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >

class ResidualBasedNewtonRaphsonMPCContactStrategy :
    public ResidualBasedNewtonRaphsonStrategy< TSparseSpace, TDenseSpace, TLinearSolver >
{
public:
    ///@name Type Definitions
    ///@{

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedNewtonRaphsonMPCContactStrategy );

    typedef SolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>            StrategyBaseType;

    typedef ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    typedef ConvergenceCriteria<TSparseSpace, TDenseSpace>               TConvergenceCriteriaType;

    typedef MPCContactCriteria<TSparseSpace, TDenseSpace>                 TMPCContactCriteriaType;

    typedef typename BaseType::TBuilderAndSolverType                        TBuilderAndSolverType;

    typedef typename BaseType::TDataType                                                TDataType;

    typedef TSparseSpace                                                          SparseSpaceType;

    typedef typename BaseType::TSchemeType                                            TSchemeType;

    typedef typename BaseType::DofsArrayType                                        DofsArrayType;

    typedef typename BaseType::TSystemMatrixType                                TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType                                TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType                        LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType                        LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType                  TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType                  TSystemVectorPointerType;

    typedef ModelPart::NodesContainerType                                          NodesArrayType;

    typedef ModelPart::ElementsContainerType                                    ElementsArrayType;

    typedef ModelPart::ConditionsContainerType                                ConditionsArrayType;

    typedef ModelPart::MasterSlaveConstraintContainerType                     ConstraintArrayType;

    typedef std::size_t                                                                 IndexType;

    typedef std::size_t                                                                  SizeType;

    /**
     * @brief Default constructor
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param MaxIterations The maximum number of iterations
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    ResidualBasedNewtonRaphsonMPCContactStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        IndexType MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false,
        Parameters ThisParameters =  Parameters(R"({})")
    )
        : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver, pNewConvergenceCriteria, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag),
        mThisParameters(ThisParameters)
    {
        KRATOS_TRY;

        // We create the contact criteria
        mpMPCContactCriteria = Kratos::make_shared<TMPCContactCriteriaType>();

        Parameters default_parameters = GetDefaultParameters();
        mThisParameters.ValidateAndAssignDefaults(default_parameters);

        KRATOS_CATCH("");
    }

    /**
     * @brief Default constructor
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param pNewLinearSolver The linear solver employed
     * @param pNewConvergenceCriteria The convergence criteria employed
     * @param MaxIterations The maximum number of iterations
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    ResidualBasedNewtonRaphsonMPCContactStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver,
        IndexType MaxIterations = 30,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false,
        Parameters ThisParameters =  Parameters(R"({})")
        )
        : ResidualBasedNewtonRaphsonStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, pScheme, pNewLinearSolver, pNewConvergenceCriteria, pNewBuilderAndSolver, MaxIterations, CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag ),
        mThisParameters(ThisParameters)
    {
        KRATOS_TRY;

        // We create the contact criteria
        mpMPCContactCriteria = Kratos::make_shared<TMPCContactCriteriaType>();

        Parameters default_parameters = GetDefaultParameters();
        mThisParameters.ValidateAndAssignDefaults(default_parameters);

        KRATOS_CATCH("");
    }

    /**
     * Destructor.
     */

    ~ResidualBasedNewtonRaphsonMPCContactStrategy() override
    = default;

    //******************** OPERATIONS ACCESSIBLE FROM THE INPUT: ************************//
    //***********************************************************************************//

    /**
     * @brief Operation to predict the solution ... if it is not called a trivial predictor is used in which the
     * values of the solution step of interest are assumed equal to the old values
     */
    void Predict() override
    {
        KRATOS_TRY

        BaseType::Predict();

        // Getting model part
        ModelPart& r_model_part = StrategyBaseType::GetModelPart();

        // We get the system
        TSystemMatrixType& rA = *BaseType::mpA;
        TSystemVectorType& rDx = *BaseType::mpDx;
        TSystemVectorType& rb = *BaseType::mpb;

        // We solve the system in order to check the active set once
        TSparseSpace::SetToZero(rA);
        TSparseSpace::SetToZero(rDx);
        TSparseSpace::SetToZero(rb);

        typename TSchemeType::Pointer p_scheme = BaseType::GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = BaseType::GetBuilderAndSolver();
        p_builder_and_solver->BuildAndSolve(p_scheme, BaseType::GetModelPart(), rA, rDx, rb);

        // Check active set
        const SizeType echo_level_convergence_criteria = BaseType::mpConvergenceCriteria->GetEchoLevel();
        BaseType::mpConvergenceCriteria->SetEchoLevel(0);
        mpMPCContactCriteria->PostCriteria(r_model_part, BaseType::GetBuilderAndSolver()->GetDofSet(), rA, rDx, rb);
        BaseType::mpConvergenceCriteria->SetEchoLevel(echo_level_convergence_criteria);

        KRATOS_CATCH("")
    }

    /**
     * @brief Initialization of member variables and prior operations
     */
    void Initialize() override
    {
        KRATOS_TRY;

        // Computing nodal weights
        ComputeNodalWeights();

        BaseType::Initialize();

        KRATOS_CATCH("");
    }

    /**
     * @brief The problem of interest is solved.
     * @details This function calls sequentially: Initialize(), InitializeSolutionStep(), Predict(),
     * SolveSolutionStep() and FinalizeSolutionStep().
     * All those functions can otherwise be called separately.
     */
    double Solve() override
    {
        this->Initialize();
        this->InitializeSolutionStep();
        this->Predict();
        this->SolveSolutionStep();
        this->FinalizeSolutionStep(); // TODO: Comment for proper work of interaction

        return 0.0;
    }

    /**
     * @brief Performs all the required operations that should be done (for each step)
     * before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override
    {
        // Computing nodal weights
        ComputeNodalWeights();

        BaseType::InitializeSolutionStep();

//         // If enforcing NTN
//         const bool enforce_ntn = mThisParameters["enforce_ntn"].GetBool();
//         if (enforce_ntn) {
//             EnforcingNTN();
//         }
    }

    /**
     * @brief Performs all the required operations that should be done (for each step)
     * after solving the solution step.
     */
    void FinalizeSolutionStep() override
    {
        KRATOS_TRY;

        BaseType::FinalizeSolutionStep();

        KRATOS_CATCH("");
    }

    /**
     * @brief Solves the current step.
     * @details This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        KRATOS_TRY;

        bool is_converged = false;

        // Getting model part
        ModelPart& r_model_part = StrategyBaseType::GetModelPart();

        // We get the process info
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();

        if (r_process_info.Is(INTERACTION)) {
            // We get the system
            TSystemMatrixType& rA = *BaseType::mpA;
            TSystemVectorType& rDx = *BaseType::mpDx;
            TSystemVectorType& rb = *BaseType::mpb;

            int inner_iteration = 0;
            const SizeType echo_level_convergence_criteria = BaseType::mpConvergenceCriteria->GetEchoLevel();
            while (!is_converged && inner_iteration < mThisParameters["inner_loop_iterations"].GetInt()) {
                ++inner_iteration;

                if (echo_level_convergence_criteria > 0 && r_model_part.GetCommunicator().MyPID() == 0 ) {
                    KRATOS_INFO("Simplified semi-smooth strategy") << BOLDFONT("INNER ITERATION: ") << inner_iteration << std::endl;
                }

                // We solve one loop
                r_process_info[NL_ITERATION_NUMBER] = 1;
                is_converged = AuxiliarSolveSolutionStep();

                // We check the convergence
                if (r_process_info[NL_ITERATION_NUMBER] == 1) r_process_info[NL_ITERATION_NUMBER] = 2; // Trigger check
                is_converged = mpMPCContactCriteria->PostCriteria(r_model_part, BaseType::GetBuilderAndSolver()->GetDofSet(), rA, rDx, rb);

                if (echo_level_convergence_criteria > 0 && r_model_part.GetCommunicator().MyPID() == 0 ) {
                    if (is_converged) KRATOS_INFO("Simplified semi-smooth strategy") << BOLDFONT("Simplified semi-smooth strategy. INNER ITERATION: ") << BOLDFONT(FGRN("CONVERGED")) << std::endl;
                    else KRATOS_INFO("Simplified semi-smooth strategy") << BOLDFONT("INNER ITERATION: ") << BOLDFONT(FRED("NOT CONVERGED")) << std::endl;
                }
            }
        } else {
            is_converged = AuxiliarSolveSolutionStep();
        }

        return is_converged;

        KRATOS_CATCH("");
    }


    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise. (auxiliar method)
     */
    bool AuxiliarSolveSolutionStep()
    {
        // Getting flag INTERACTION
        ModelPart& r_model_part = StrategyBaseType::GetModelPart();
        const bool update_each_nl_iteration = mThisParameters["update_each_nl_iteration"].GetBool();
        VariableUtils().SetFlag(INTERACTION, update_each_nl_iteration, r_model_part.GetSubModelPart("ComputingContact").Conditions());

        // Pointers needed in the solution
        typename TSchemeType::Pointer p_scheme = this->GetScheme();
        typename TBuilderAndSolverType::Pointer p_builder_and_solver = this->GetBuilderAndSolver();
        auto& r_dof_set = p_builder_and_solver->GetDofSet();

        TSystemMatrixType& rA  = *BaseType::mpA;
        TSystemVectorType& rDx = *BaseType::mpDx;
        TSystemVectorType& rb  = *BaseType::mpb;

        // Initializing the parameters of the Newton-Raphson cycle
        unsigned int iteration_number = 1;
        r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
        bool is_converged = false;
        bool residual_is_updated = false;

        // Computing nodal weights
        ComputeNodalWeights();

        p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
        BaseType::mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);
        is_converged = BaseType::mpConvergenceCriteria->PreCriteria(r_model_part, p_builder_and_solver->GetDofSet(), rA, rDx, rb);

//         // If enforcing NTN
//         const bool enforce_ntn = mThisParameters["enforce_ntn"].GetBool();
//         if (enforce_ntn) {
//             EnforcingNTN();
//         }

        // Function to perform the building and the solving phase.
        if (StrategyBaseType::mRebuildLevel > 0 || StrategyBaseType::mStiffnessMatrixIsBuilt == false) {
            TSparseSpace::SetToZero(rA);
            TSparseSpace::SetToZero(rDx);
            TSparseSpace::SetToZero(rb);

            p_builder_and_solver->BuildAndSolve(p_scheme, r_model_part, rA, rDx, rb);
        } else {
            TSparseSpace::SetToZero(rDx); //Dx=0.00;
            TSparseSpace::SetToZero(rb);

            p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
        }

        // Debugging info
        BaseType::EchoInfo(iteration_number);

        // Updating the results stored in the database
        BaseType::UpdateDatabase(rA, rDx, rb, StrategyBaseType::MoveMeshFlag());

        p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
        BaseType::mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

        // Calculate reactions if required
        if (BaseType::mCalculateReactionsFlag)
            p_builder_and_solver->CalculateReactions(p_scheme, r_model_part, rA, rDx, rb);

        if (is_converged) {
            if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag()) {
                TSparseSpace::SetToZero(rb);

                p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
            }

            is_converged = BaseType::mpConvergenceCriteria->PostCriteria(r_model_part, p_builder_and_solver->GetDofSet(), rA, rDx, rb);
        }

        // Iteration Cycle... performed only for NonLinearProblems
        while (!is_converged && iteration_number++ < BaseType::mMaxIterationNumber) {
            // Setting the number of iteration
            r_model_part.GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;

            // Computing nodal weights
            ComputeNodalWeights();

            // Calling InitializeNonLinIteration
            p_scheme->InitializeNonLinIteration(r_model_part, rA, rDx, rb);
            BaseType::mpConvergenceCriteria->InitializeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            // Shaping correctly the system
            if (update_each_nl_iteration) {
                p_builder_and_solver->SetUpDofSet(p_scheme, r_model_part);
                p_builder_and_solver->SetUpSystem(r_model_part);
                p_builder_and_solver->ResizeAndInitializeVectors(p_scheme, BaseType::mpA, BaseType::mpDx, BaseType::mpb, r_model_part);
            }

            is_converged = BaseType::mpConvergenceCriteria->PreCriteria(r_model_part, p_builder_and_solver->GetDofSet(), rA, rDx, rb);

            // Call the linear system solver to find the correction mDx for the it is not called if there is no system to solve
            if (SparseSpaceType::Size(rDx) != 0) {
                if (StrategyBaseType::mRebuildLevel > 1 || !StrategyBaseType::mStiffnessMatrixIsBuilt) {
                    if (!BaseType::GetKeepSystemConstantDuringIterations()) {
                        //A = 0.00;
                        TSparseSpace::SetToZero(rA);
                        TSparseSpace::SetToZero(rDx);
                        TSparseSpace::SetToZero(rb);

                        p_builder_and_solver->BuildAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                    } else {
                        TSparseSpace::SetToZero(rDx);
                        TSparseSpace::SetToZero(rb);

                        p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                    }
                } else {
                    TSparseSpace::SetToZero(rDx);
                    TSparseSpace::SetToZero(rb);

                    p_builder_and_solver->BuildRHSAndSolve(p_scheme, r_model_part, rA, rDx, rb);
                }
            } else {
                KRATOS_WARNING("NO DOFS") << "ATTENTION: no free DOFs!! " << std::endl;
            }

            // Debugging info
            BaseType::EchoInfo(iteration_number);

            // Updating the results stored in the database
            BaseType::UpdateDatabase(rA, rDx, rb, StrategyBaseType::MoveMeshFlag());

            p_scheme->FinalizeNonLinIteration(r_model_part, rA, rDx, rb);
            BaseType::mpConvergenceCriteria->FinalizeNonLinearIteration(r_model_part, r_dof_set, rA, rDx, rb);

            residual_is_updated = false;

            // Calculate reactions if required
            if (BaseType::mCalculateReactionsFlag)
                p_builder_and_solver->CalculateReactions(p_scheme, r_model_part, rA, rDx, rb);

            if (is_converged) {
                if (BaseType::mpConvergenceCriteria->GetActualizeRHSflag()) {
                    TSparseSpace::SetToZero(rb);

                    p_builder_and_solver->BuildRHS(p_scheme, r_model_part, rb);
                    residual_is_updated = true;
                }

                is_converged = BaseType::mpConvergenceCriteria->PostCriteria(r_model_part, p_builder_and_solver->GetDofSet(), rA, rDx, rb);
            }
        }

        // Plots a warning if the maximum number of iterations is exceeded
        if (iteration_number >= BaseType::mMaxIterationNumber) {
            BaseType::MaxIterationsExceeded();
        } else {
            KRATOS_INFO_IF("NR-Strategy", this->GetEchoLevel() > 0)  << "Convergence achieved after " << iteration_number << " / "  << BaseType::mMaxIterationNumber << " iterations" << std::endl;
        }

        // Recalculate residual if needed (note that some convergence criteria need it to be recalculated)
        if (!residual_is_updated) {
            // NOTE:
            // The following part will be commented because it is time consuming
            // and there is no obvious reason to be here. If someone need this
            // part please notify the community via mailing list before uncommenting it.
            // Pooyan.

//             TSparseSpace::SetToZero(mb);
//             p_builder_and_solver->BuildRHS(p_scheme, r_model_part, mb);
        }

        // Calculate reactions if required
        if (BaseType::mCalculateReactionsFlag)
            p_builder_and_solver->CalculateReactions(p_scheme, r_model_part, rA, rDx, rb);

        return is_converged;
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

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    Parameters mThisParameters;                                   /// The configuration parameters
    typename TConvergenceCriteriaType::Pointer mpMPCContactCriteria; /// The contact criteria

    ///@}
    ///@name Protected Operators
    ///@{

    /**
     * @brief This method returns the defaulr parameters in order to avoid code duplication
     * @return Returns the default parameters
     */

    Parameters GetDefaultParameters()
    {
        Parameters default_parameters = Parameters(R"(
        {
            "inner_loop_iterations"    : 5,
            "update_each_nl_iteration" : false,
            "enforce_ntn"              : false
        })" );

        return default_parameters;
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
    ///@{

    /**
     * Copy constructor.
     */

    ResidualBasedNewtonRaphsonMPCContactStrategy(const ResidualBasedNewtonRaphsonMPCContactStrategy& Other)
    {
    };

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

//     /**
//      * @brief This inforces NTN formulation
//      */
//     void EnforcingNTN()
//     {
//         //  List of enforced nodes to not repeat
//         std::unordered_set<IndexType> enforced_nodes;
//
//         // Getting contact model part
//         ModelPart& r_root_model_part = StrategyBaseType::GetModelPart().GetRootModelPart();
//         ModelPart& r_computing_contact_model_part = StrategyBaseType::GetModelPart().GetSubModelPart("ComputingContact");
//
//         // The process info
//         const auto& r_process_info = r_root_model_part.GetProcessInfo();
//
//         // Reset the pointers of the conditions
//         for (auto& r_cond : r_computing_contact_model_part.Conditions()) {
//             if (r_cond.Has(CONSTRAINT_POINTER)) {
//                 r_cond.SetValue(CONSTRAINT_POINTER, nullptr);
//             }
//         }
//
//         // Iterate over the constraints
//         IndexType counter = 1;
//         for (auto& r_const : r_root_model_part.MasterSlaveConstraints()) {
//             r_const.SetId(counter);
//             ++counter;
//         }
//
//         // Auxiliar classes
//         Matrix original_relation_matrix, relation_matrix;
//         Vector original_constant_vector, constant_vector;
//         ModelPart::DofsVectorType original_master_dofs, master_dofs, original_slave_dofs, slave_dofs;
//
//         // Iterate over the constraints
//         for (auto& r_const : r_computing_contact_model_part.MasterSlaveConstraints()) {
//             // Getting original system
//             r_const.GetLocalSystem(original_relation_matrix, original_constant_vector, r_process_info);
//             r_const.GetDofList(original_slave_dofs, original_master_dofs, r_process_info);
//
//             // TODO: Finish rebuild
//
//             // Creating new constraint
//             r_root_model_part.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", counter, master_dofs, slave_dofs, relation_matrix, constant_vector);
//
//             // Setting to remove the old constraints
//             r_const.Set(TO_ERASE, true);
//
//             ++counter;
//         }
//
//         // Remove old constraints
//         r_root_model_part.RemoveMasterSlaveConstraintsFromAllLevels(TO_ERASE);
//
//         // Transfer constraints from the root to the computing model part
//         FastTransferBetweenModelPartsProcess(r_computing_contact_model_part, r_root_model_part, FastTransferBetweenModelPartsProcess::EntityTransfered::CONSTRAINTS).Execute();
//
//         // Reorder ids
//         counter = 1;
//         for (auto& r_const : r_root_model_part.MasterSlaveConstraints()) {
//             r_const.SetId(counter);
//             ++counter;
//         }
//     }

    /**
     * @brief This computes the nodal weights
     */
    void ComputeNodalWeights()
    {
        // Getting contact model part
        ModelPart& r_contact_model_part = StrategyBaseType::GetModelPart().GetSubModelPart("Contact");

        // Reset the NODAL_PAUX and NODAL_MAUX
        auto& r_nodes_array = r_contact_model_part.Nodes();
        VariableUtils().SetNonHistoricalVariableToZero(NODAL_PAUX, r_nodes_array);
        VariableUtils().SetNonHistoricalVariableToZero(NODAL_MAUX, r_nodes_array);

        // We set the constraints active and inactive in function of the active set
        auto& r_conditions_array = r_contact_model_part.Conditions();
        auto it_cond_begin = r_conditions_array.begin();

        // If enforcing NTN
        const bool enforce_ntn = false;
//         const bool enforce_ntn = mThisParameters["enforce_ntn"].GetBool();
//         if (enforce_ntn) {
//             VariableUtils().SetNonHistoricalVariable(NODAL_PAUX, 1.0, r_nodes_array);
//         }

        #pragma omp parallel for
        for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i) {
            auto it_cond = it_cond_begin + i;

            // Only slave conditions
            if (it_cond->Is(SLAVE)) {
                auto& r_geometry = it_cond->GetGeometry();
                Vector lumping_factor;
                lumping_factor = r_geometry.LumpingFactors(lumping_factor);
                const double domain_size = r_geometry.DomainSize();
                for (IndexType i_node = 0; i_node < r_geometry.size(); ++i_node) {
                    auto& r_node = r_geometry[i_node];
                    if (!enforce_ntn) {
                        #pragma omp atomic
                        r_node.GetValue(NODAL_PAUX) += 1.0;
                    }
                    #pragma omp atomic
                    r_node.GetValue(NODAL_MAUX) += lumping_factor[i_node] * domain_size;
                }
            }
        }
    }

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; /* Class ResidualBasedNewtonRaphsonMPCContactStrategy */
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}
}  // namespace Kratos

#endif /* KRATOS_RESIDUALBASED_NEWTON_RAPHSON_MPC_CONTACT_STRATEGY */
