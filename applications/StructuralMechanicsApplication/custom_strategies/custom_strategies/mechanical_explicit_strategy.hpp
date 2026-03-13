// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Klaus B Sautter (based on the work of JMCarbonel)
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "solving_strategies/strategies/implicit_solving_strategy.h"
#include "structural_mechanics_application_variables.h"
#include "utilities/variable_utils.h"
#include "utilities/constraint_utilities.h"
#include "utilities/parallel_utilities.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

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
 * @class MechanicalExplicitStrategy
 * @ingroup StructuralMechanicsApplciation
 * @brief This strategy is used for the explicit time integration
 * @author Klauss B Sautter (based on the work of JMCarbonel)
 */
template <class TSparseSpace,
          class TDenseSpace,  // = DenseSpace<double>,
          class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
          >
class MechanicalExplicitStrategy
    : public ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> {
public:
    ///@name Type Definitions
    ///@{

    // Base class definition
    typedef ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;

    /// Some definitions from the base class
    typedef typename BaseType::TSchemeType TSchemeType;
    typedef typename BaseType::DofsArrayType DofsArrayType;
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
    typedef typename BaseType::NodesArrayType NodesArrayType;
    typedef typename BaseType::ElementsArrayType ElementsArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    /// DoF types definition
    typedef typename Node::DofType DofType;
    typedef typename DofType::Pointer DofPointerType;

    /// Counted pointer of MechanicalExplicitStrategy
    KRATOS_CLASS_POINTER_DEFINITION(MechanicalExplicitStrategy);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param rModelPart The model part of the problem
     * @param pScheme The integration scheme
     * @param CalculateReactions The flag for the reaction calculation
     * @param ReformDofSetAtEachStep The flag that allows to compute the modification of the DOF
     * @param MoveMeshFlag The flag that allows to move the mesh
     */
    MechanicalExplicitStrategy(
        ModelPart& rModelPart,
        typename TSchemeType::Pointer pScheme,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = true)
        : ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(rModelPart, MoveMeshFlag),
          mpScheme(pScheme),
          mReformDofSetAtEachStep(ReformDofSetAtEachStep),
          mCalculateReactionsFlag(CalculateReactions)
    {
        KRATOS_TRY

        // Set EchoLevel to the default value (only time is displayed)
        BaseType::SetEchoLevel(1);

        // Set RebuildLevel to the default value
        BaseType::SetRebuildLevel(0);

        KRATOS_CATCH("")
    }

    /** Destructor.
    */
    virtual ~MechanicalExplicitStrategy()
    {
        Clear();
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Set method for the time scheme
     * @param pScheme The pointer to the time scheme considered
     */
    void SetScheme(typename TSchemeType::Pointer pScheme)
    {
        mpScheme = pScheme;
    };

    /**
     * @brief Get method for the time scheme
     * @return mpScheme: The pointer to the time scheme considered
     */
    typename TSchemeType::Pointer GetScheme()
    {
        return mpScheme;
    };

    // Set and Get Flags

    /**
     * @brief This method sets the flag mInitializeWasPerformed
     * @param InitializePerformedFlag The flag that tells if the initialize has been computed
     */
    void SetInitializePerformedFlag(bool InitializePerformedFlag = true)
    {
        mInitializeWasPerformed = InitializePerformedFlag;
    }

    /**
     * @brief This method gets the flag mInitializeWasPerformed
     * @return mInitializeWasPerformed: The flag that tells if the initialize has been computed
     */
    bool GetInitializePerformedFlag()
    {
        return mInitializeWasPerformed;
    }

    /**
     * @brief This method sets the flag mCalculateReactionsFlag
     * @param CalculateReactionsFlag The flag that tells if the reactions are computed
     */
    void SetCalculateReactionsFlag(bool CalculateReactionsFlag)
    {
        mCalculateReactionsFlag = CalculateReactionsFlag;
    }

    /**
     * @brief This method returns the flag mCalculateReactionsFlag
     * @return The flag that tells if the reactions are computed
     */
    bool GetCalculateReactionsFlag()
    {
        return mCalculateReactionsFlag;
    }

    /**
     * @brief This method sets the flag mReformDofSetAtEachStep
     * @param Flag The flag that tells if each time step the system is rebuilt
     */
    void SetReformDofSetAtEachStepFlag(bool Flag)
    {
        mReformDofSetAtEachStep = Flag;
    }

    /**
     * @brief This method returns the flag mReformDofSetAtEachStep
     * @return The flag that tells if each time step the system is rebuilt
     */
    bool GetReformDofSetAtEachStepFlag()
    {
        return mReformDofSetAtEachStep;
    }
//
    /**
     * @brief Initialization of member variables and prior operations
     */
    void Initialize() override
    {
        KRATOS_TRY

        if (!this->mInitializeWasPerformed){
            // Pointers needed in the solution
            typename TSchemeType::Pointer pScheme = GetScheme();
            ModelPart& r_model_part = BaseType::GetModelPart();

            TSystemMatrixType matrix_a_dummy = TSystemMatrixType();

            // Initialize The Scheme - OPERATIONS TO BE DONE ONCE
            if (!pScheme->SchemeIsInitialized())pScheme->Initialize(r_model_part);

            // Initialize The Elements - OPERATIONS TO BE DONE ONCE
            if (!pScheme->ElementsAreInitialized())pScheme->InitializeElements(r_model_part);

            // Initialize The Conditions- OPERATIONS TO BE DONE ONCE
            if (!pScheme->ConditionsAreInitialized())pScheme->InitializeConditions(r_model_part);

            // Set Nodal Mass and Damping to zero
            NodesArrayType& r_nodes = r_model_part.Nodes();
            VariableUtils().SetNonHistoricalVariable(NODAL_MASS, 0.0, r_nodes);
            VariableUtils().SetNonHistoricalVariable(NODAL_DISPLACEMENT_DAMPING, 0.0, r_nodes);

            // Iterate over the elements
            ElementsArrayType& r_elements = r_model_part.Elements();
            ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();

            Vector dummy_vector;
            // If we consider the rotation DoF
            const bool has_dof_for_rot_z = !r_nodes.empty() && r_nodes.begin()->HasDofFor(ROTATION_Z);
            if (has_dof_for_rot_z) {
                const array_1d<double, 3> zero_array = ZeroVector(3);
                VariableUtils().SetNonHistoricalVariable(NODAL_INERTIA, zero_array, r_nodes);
                VariableUtils().SetNonHistoricalVariable(NODAL_ROTATION_DAMPING, zero_array, r_nodes);

                block_for_each(r_elements, dummy_vector, [&r_current_process_info](Element& r_element, Vector& r_dummy_vector){
                    // Getting nodal mass and inertia from element
                    // this function needs to be implemented in the respective
                    // element to provide inertias and nodal masses
                    r_element.AddExplicitContribution(r_dummy_vector, RESIDUAL_VECTOR, NODAL_INERTIA, r_current_process_info);
                });
            } else { // Only NODAL_MASS is needed
                block_for_each(r_elements, dummy_vector, [&r_current_process_info](Element& r_element, Vector& r_dummy_vector){
                    // Getting nodal mass and inertia from element
                    // this function needs to be implemented in the respective
                    // element to provide inertias and nodal masses
                    r_element.AddExplicitContribution(r_dummy_vector, RESIDUAL_VECTOR, NODAL_MASS, r_current_process_info);
                });
            }

            // Precompute for masses and inertias
            if(r_model_part.MasterSlaveConstraints().size() > 0) {
                ConstraintUtilities::PreComputeExplicitConstraintMassAndInertia(r_model_part);
            }

            this->mInitializeWasPerformed = true;
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override {
        KRATOS_TRY

        typename TSchemeType::Pointer pScheme = GetScheme();
        ModelPart& r_model_part = BaseType::GetModelPart();

        TSystemMatrixType matrix_a_dummy = TSystemMatrixType();
        TSystemVectorType rDx = TSystemVectorType();
        TSystemVectorType rb = TSystemVectorType();

        // Initial operations ... things that are constant over the Solution Step
        pScheme->InitializeSolutionStep(r_model_part, matrix_a_dummy, rDx, rb);

        if (BaseType::mRebuildLevel > 0) { // TODO: Right now is computed in the Initialize() because is always zero, the option to set the RebuildLevel should be added in the constructor or in some place
            ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
            ElementsArrayType& r_elements = r_model_part.Elements();

            // Set Nodal Mass and Damping to zero
            NodesArrayType& r_nodes = r_model_part.Nodes();
            VariableUtils().SetNonHistoricalVariable(NODAL_MASS, 0.0, r_nodes);
            VariableUtils().SetNonHistoricalVariable(NODAL_DISPLACEMENT_DAMPING, 0.0, r_nodes);

            Vector dummy_vector;
            // If we consider the rotation DoF
            const bool has_dof_for_rot_z = r_model_part.Nodes().begin()->HasDofFor(ROTATION_Z);
            if (has_dof_for_rot_z) {
                const array_1d<double, 3> zero_array = ZeroVector(3);
                VariableUtils().SetNonHistoricalVariable(NODAL_INERTIA, zero_array, r_nodes);
                VariableUtils().SetNonHistoricalVariable(NODAL_ROTATION_DAMPING, zero_array, r_nodes);

                block_for_each(r_elements, dummy_vector, [&r_current_process_info](Element& r_element, Vector& r_dummy_vector){
                    // Getting nodal mass and inertia from element
                    // this function needs to be implemented in the respective
                    // element to provide inertias and nodal masses
                    r_element.AddExplicitContribution(r_dummy_vector, RESIDUAL_VECTOR, NODAL_INERTIA, r_current_process_info);
                });
            } else { // Only NODAL_MASS and NODAL_DISPLACEMENT_DAMPING are needed
                block_for_each(r_elements, dummy_vector, [&r_current_process_info](Element& r_element, Vector& r_dummy_vector){
                    // Getting nodal mass and inertia from element
                    // this function needs to be implemented in the respective
                    // element to provide inertias and nodal masses
                    r_element.AddExplicitContribution(r_dummy_vector, RESIDUAL_VECTOR, NODAL_MASS, r_current_process_info);
                });
            }

            // Precompute for masses and inertias
            if(r_model_part.MasterSlaveConstraints().size() > 0) {
                ConstraintUtilities::PreComputeExplicitConstraintMassAndInertia(r_model_part);
            }
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        typename TSchemeType::Pointer pScheme = GetScheme();
        ModelPart& r_model_part = BaseType::GetModelPart();

        // Some dummy sets and matrices
        DofsArrayType dof_set_dummy;
        TSystemMatrixType rA = TSystemMatrixType();
        TSystemVectorType rDx = TSystemVectorType();
        TSystemVectorType rb = TSystemVectorType();

        // Initialize the non linear iteration
        pScheme->InitializeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);

        pScheme->Predict(r_model_part, dof_set_dummy, rA, rDx, rb);

        // Pre-compute MPC contributions
        if(r_model_part.MasterSlaveConstraints().size() > 0) {
            std::vector<std::string> dof_variable_names(2);
            dof_variable_names[0] = "DISPLACEMENT";
            dof_variable_names[1] = "ROTATION";
            std::vector<std::string> residual_variable_names(2);
            residual_variable_names[0] = "FORCE_RESIDUAL";
            residual_variable_names[1] = "MOMENT_RESIDUAL";
            ConstraintUtilities::PreComputeExplicitConstraintConstribution(r_model_part, dof_variable_names, residual_variable_names);
        }

        // Explicitly integrates the equation of motion.
        pScheme->Update(r_model_part, dof_set_dummy, rA, rDx, rb);

        // Finalize the non linear iteration
        pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);

        // We compute the MPC contributions
        if(r_model_part.MasterSlaveConstraints().size() > 0) {
            ComputeExplicitConstraintConstribution(pScheme, r_model_part);
        }

        // Calculate reactions if required
        if (mCalculateReactionsFlag) {
            CalculateReactions(pScheme, r_model_part, rA, rDx, rb);
        }

        return true;
    }

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void FinalizeSolutionStep() override
    {
        typename TSchemeType::Pointer pScheme = GetScheme();
        ModelPart& r_model_part = BaseType::GetModelPart();
        TSystemMatrixType rA = TSystemMatrixType();
        TSystemVectorType rDx = TSystemVectorType();
        TSystemVectorType rb = TSystemVectorType();
        // Finalisation of the solution step,
        // operations to be done after achieving convergence, for example the
        // Final Residual Vector (rb) has to be saved in there
        // to avoid error accumulation
        pScheme->FinalizeSolutionStep(r_model_part, rA, rDx, rb);

        // Move the mesh if needed
        if (BaseType::MoveMeshFlag())
            BaseType::MoveMesh();

        // Cleaning memory after the solution
        pScheme->Clean();
    }

    /**
     * @brief Clears the internal storage
     */
    void Clear() override
    {
        KRATOS_TRY

        KRATOS_INFO("MechanicalExplicitStrategy") << "Clear function used" << std::endl;

        GetScheme()->Clear();
        mInitializeWasPerformed = false;

        KRATOS_CATCH("")
    }

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * on the input provided.
     * @details Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rModelPart The model of the problem to solve
     * @return Zero means  all ok
     */
    int Check() override
    {
        KRATOS_TRY

        BaseType::Check();

        GetScheme()->Check(BaseType::GetModelPart());

        return 0;

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
     * @brief This method computes the explicit contribution of the constraints
     * @param pScheme The pointer to the integration scheme used
     * @param rModelPart The model part which defines the problem
     */
    void ComputeExplicitConstraintConstribution(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart
        )
    {
        // First we reset the slave dofs
        ConstraintUtilities::ResetSlaveDofs(rModelPart);

        // Now we apply the constraints
        ConstraintUtilities::ApplyConstraints(rModelPart);
    }

    /**
     * @brief This method computes the reactions of the problem
     * @param pScheme The pointer to the integration scheme used
     * @param rModelPart The model part which defines the problem
     * @param rA The LHS of the system (empty)
     * @param rDx The solution of the system (empty)
     * @param rb The RHS of the system (empty)
     */
    void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart,
        TSystemMatrixType& rA,
        TSystemVectorType& rDx,
        TSystemVectorType& rb
        )
    {
        // We iterate over the nodes
        auto& r_nodes = rModelPart.Nodes();

        if (!r_nodes.empty()) {
            // If we consider rotation dofs
            const bool has_dof_for_rot_z = (r_nodes.begin())->HasDofFor(ROTATION_Z);

            // Auxiliary values
            const array_1d<double,3> zero_array = ZeroVector(3);

            // Getting
            const auto it_node_begin = r_nodes.begin();
            const IndexType disppos = it_node_begin->GetDofPosition(DISPLACEMENT_X);
            const IndexType rotppos = it_node_begin->GetDofPosition(ROTATION_X);

            // Construct loop lambda depending on whether nodes have rot_z
            std::function<void(Node&)> loop_base, loop;

            loop_base = [&disppos](Node& rNode){
                const auto force_residual = rNode.FastGetSolutionStepValue(FORCE_RESIDUAL);

                if (rNode.GetDof(DISPLACEMENT_X, disppos).IsFixed()) {
                    double& r_reaction = rNode.FastGetSolutionStepValue(REACTION_X);
                    r_reaction = -force_residual[0];
                }
                if (rNode.GetDof(DISPLACEMENT_Y, disppos + 1).IsFixed()) {
                    double& r_reaction = rNode.FastGetSolutionStepValue(REACTION_Y);
                    r_reaction = -force_residual[1];
                }
                if (rNode.GetDof(DISPLACEMENT_Z, disppos + 2).IsFixed()) {
                    double& r_reaction = rNode.FastGetSolutionStepValue(REACTION_Z);
                    r_reaction = -force_residual[2];
                }
            };

            if (has_dof_for_rot_z) {
                loop = [&rotppos, &loop_base](Node& rNode){
                    loop_base(rNode);
                    const auto moment_residual = rNode.FastGetSolutionStepValue(MOMENT_RESIDUAL);
                    if (rNode.GetDof(ROTATION_X, rotppos).IsFixed()) {
                        double& r_reaction = rNode.FastGetSolutionStepValue(REACTION_MOMENT_X);
                        r_reaction = -moment_residual[0];
                    }
                    if (rNode.GetDof(ROTATION_Y, rotppos + 1).IsFixed()) {
                        double& r_reaction = rNode.FastGetSolutionStepValue(REACTION_MOMENT_Y);
                        r_reaction = -moment_residual[1];
                    }
                    if (rNode.GetDof(ROTATION_Z, rotppos + 2).IsFixed()) {
                        double& r_reaction = rNode.FastGetSolutionStepValue(REACTION_MOMENT_Z);
                        r_reaction = -moment_residual[2];
                    }
                };
            } else {
                loop = loop_base;
            }

            // Compute on nodes
            block_for_each(r_nodes, loop);
        } // if not nodes.empty()
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

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    typename TSchemeType::Pointer mpScheme; /// The pointer to the integration scheme

    /**
     * @brief Flag telling if it is needed to reform the DofSet at each solution step or if it is possible to form it just once
     * @details
        - true  => reforme at each time step
        - false => form just one (more efficient)
        Default = false
    */
    bool mReformDofSetAtEachStep = false;

    /**
     * @brief Flag telling if it is needed or not to compute the reactions
     * @details default = true
     */
    bool mCalculateReactionsFlag = true;

    /**
     * @brief Flag telling if the initialize was performed
     * @details default = false
     */
    bool mInitializeWasPerformed = false;

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

    /** Copy constructor.
    */
    MechanicalExplicitStrategy(const MechanicalExplicitStrategy& Other){};

    ///@}

}; /* Class MechanicalExplicitStrategy */

///@}

} /* namespace Kratos.*/
