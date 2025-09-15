
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//


#if !defined(KRATOS_POROMECHANICS_EXPLICIT_STRATEGY)
#define KRATOS_POROMECHANICS_EXPLICIT_STRATEGY

/* System includes */
// #include <fstream>

// Project includes
#include "custom_strategies/custom_strategies/mechanical_explicit_strategy.hpp"

// Application includes
#include "poromechanics_application_variables.h"

namespace Kratos {

template <class TSparseSpace,
          class TDenseSpace,
          class TLinearSolver
          >
class PoromechanicsExplicitStrategy
    : public MechanicalExplicitStrategy<TSparseSpace, TDenseSpace, TLinearSolver> {
public:

    KRATOS_CLASS_POINTER_DEFINITION(PoromechanicsExplicitStrategy);

    typedef ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
    typedef MechanicalExplicitStrategy<TSparseSpace, TDenseSpace, TLinearSolver> MotherType;
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

    using MotherType::mInitializeWasPerformed;
    using MotherType::mCalculateReactionsFlag;
    using MotherType::mpScheme;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ///Constructor
    PoromechanicsExplicitStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        Parameters& rParameters,
        bool CalculateReactions = false,
        bool ReformDofSetAtEachStep = false,
        bool MoveMeshFlag = false
        ) : MechanicalExplicitStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part, pScheme,
                CalculateReactions, ReformDofSetAtEachStep, MoveMeshFlag)
        {
            Parameters default_parameters( R"(
                {
                    "rebuild_level": 0,
                    "max_radius_factor": 20.0,
                    "min_radius_factor": 0.5,
                    "initial_radius": 1.0e-12,
                    "characteristic_length": 0.05,
                    "search_neighbours_step": false,
                    "body_domain_sub_model_part_list": [],
                    "loads_sub_model_part_list": [],
                    "loads_variable_list" : []
                }  )" );

            // Validate agains defaults -- this also ensures no type mismatch
            rParameters.ValidateAndAssignDefaults(default_parameters);

            mpParameters = &rParameters;

            // Set Load SubModelParts and Variable names
            if(rParameters["loads_sub_model_part_list"].size() > 0)
            {
                mSubModelPartList.resize(rParameters["loads_sub_model_part_list"].size());
                mVariableNames.resize(rParameters["loads_variable_list"].size());

                if( mSubModelPartList.size() != mVariableNames.size() )
                    KRATOS_THROW_ERROR( std::logic_error, "For each SubModelPart there must be a corresponding nodal Variable", "" )

                for(unsigned int i = 0; i < mVariableNames.size(); i++)
                {
                    mSubModelPartList[i] = &( model_part.GetSubModelPart(rParameters["loads_sub_model_part_list"][i].GetString()) );
                    mVariableNames[i] = rParameters["loads_variable_list"][i].GetString();
                }
            }

            BaseType::SetRebuildLevel(rParameters["rebuild_level"].GetInt());
        }

    //------------------------------------------------------------------------------------

    ///Destructor
    ~PoromechanicsExplicitStrategy() override {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * @brief Initialization of member variables and prior operations
     */
    void Initialize() override
    {
        KRATOS_TRY

        if (mInitializeWasPerformed == false) {

            ModelPart& r_model_part = BaseType::GetModelPart();

            TSystemMatrixType matrix_a_dummy = TSystemMatrixType();

            // Initialize The Scheme - OPERATIONS TO BE DONE ONCE
            if (!mpScheme->SchemeIsInitialized())mpScheme->Initialize(r_model_part);

            // Initialize The Elements - OPERATIONS TO BE DONE ONCE
            if (!mpScheme->ElementsAreInitialized())mpScheme->InitializeElements(r_model_part);

            // Initialize The Conditions- OPERATIONS TO BE DONE ONCE
            if (!mpScheme->ConditionsAreInitialized())mpScheme->InitializeConditions(r_model_part);

            // Set Nodal Mass to zero
            NodesArrayType& r_nodes = r_model_part.Nodes();
            VariableUtils().SetNonHistoricalVariable(NODAL_MASS, 0.0, r_nodes);
            // TODO: Set Nodal AntiCompressibility to zero for mass-balance equation (C=1/Q, with Q being the compressibility coeff.)

            // Iterate over the elements
            ElementsArrayType& r_elements = r_model_part.Elements();
            const auto it_elem_begin = r_elements.begin();
            ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();

            Vector dummy_vector;
            #pragma omp parallel for firstprivate(dummy_vector), schedule(guided,512)
            for (int i = 0; i < static_cast<int>(r_elements.size()); ++i) {
                // Getting nodal mass and inertia from element
                // this function needs to be implemented in the respective
                // element to provide nodal masses
                auto it_elem = it_elem_begin + i;
                it_elem->AddExplicitContribution(dummy_vector, RESIDUAL_VECTOR, NODAL_MASS, r_current_process_info);
            }

            //Initialize ProcessInfo variables
            r_current_process_info[IS_CONVERGED] = true;

            mInitializeWasPerformed = true;
        }

        KRATOS_CATCH( "" )
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void InitializeSolutionStep() override {
        KRATOS_TRY

        ModelPart& r_model_part = BaseType::GetModelPart();

        TSystemMatrixType matrix_a_dummy = TSystemMatrixType();
        TSystemVectorType rDx = TSystemVectorType();
        TSystemVectorType rb = TSystemVectorType();

        // Initial operations ... things that are constant over the Solution Step
        mpScheme->InitializeSolutionStep(r_model_part, matrix_a_dummy, rDx, rb);

        if (BaseType::mRebuildLevel > 0) {
            ProcessInfo& r_current_process_info = r_model_part.GetProcessInfo();
            ElementsArrayType& r_elements = r_model_part.Elements();
            const auto it_elem_begin = r_elements.begin();

            // Set Nodal Mass and Damping to zero
            NodesArrayType& r_nodes = r_model_part.Nodes();
            VariableUtils().SetNonHistoricalVariable(NODAL_MASS, 0.0, r_nodes);

            Vector dummy_vector;
            #pragma omp parallel for firstprivate(dummy_vector), schedule(guided,512)
            for (int i = 0; i < static_cast<int>(r_elements.size()); ++i) {
                // Getting nodal mass and inertia from element
                // this function needs to be implemented in the respective
                // element to provide nodal masses
                auto it_elem = it_elem_begin + i;
                it_elem->AddExplicitContribution(dummy_vector, RESIDUAL_VECTOR, NODAL_MASS, r_current_process_info);
            }
        }

        KRATOS_CATCH("")
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    bool SolveSolutionStep() override
    {
        ModelPart& r_model_part = BaseType::GetModelPart();

        // Some dummy sets and matrices
        DofsArrayType dof_set_dummy;
        TSystemMatrixType rA = TSystemMatrixType();
        TSystemVectorType rDx = TSystemVectorType();
        TSystemVectorType rb = TSystemVectorType();

        // Initialize the non linear iteration
        mpScheme->InitializeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);

        mpScheme->Predict(r_model_part, dof_set_dummy, rA, rDx, rb);

        // Move the mesh if needed
        if (BaseType::MoveMeshFlag())
            BaseType::MoveMesh();

        // Explicitly integrates the equation of motion.
        mpScheme->Update(r_model_part, dof_set_dummy, rA, rDx, rb);

        // CONVERGENCE CHECK
        this->CheckConvergence(r_model_part);

        // Move the mesh if needed
        if (BaseType::MoveMeshFlag())
            BaseType::MoveMesh();

        // Finalize the non linear iteration
        mpScheme->FinalizeNonLinIteration(BaseType::GetModelPart(), rA, rDx, rb);

        // Calculate reactions if required
        if (mCalculateReactionsFlag) {
            this->CalculateReactions(mpScheme, r_model_part);
        }

        return true;
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    void FinalizeSolutionStep() override
    {
        ModelPart& r_model_part = BaseType::GetModelPart();
        TSystemMatrixType rA = TSystemMatrixType();
        TSystemVectorType rDx = TSystemVectorType();
        TSystemVectorType rb = TSystemVectorType();
        // Finalisation of the solution step,
        // operations to be done after achieving convergence, for example the
        // Final Residual Vector (rb) has to be saved in there
        // to avoid error accumulation
        mpScheme->FinalizeSolutionStep(r_model_part, rA, rDx, rb);

        // Cleaning memory after the solution
        mpScheme->Clean();
    }

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

protected:

    /// Member Variables
    Parameters* mpParameters;
    std::vector<ModelPart*> mSubModelPartList; /// List of every SubModelPart associated to an external load
    std::vector<std::string> mVariableNames; /// Name of the nodal variable associated to every SubModelPart

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /**
     * @brief This method computes the reactions of the problem
     * @param pScheme The pointer to the integration scheme used
     * @param rModelPart The model part which defines the problem
     * @param rA The LHS of the system (empty)
     * @param rDx The solution of the system (empty)
     * @param rb The RHS of the system (empty)
     */
    virtual void CalculateReactions(
        typename TSchemeType::Pointer pScheme,
        ModelPart& rModelPart
        )
    {
        // We iterate over the nodes
        auto& r_nodes = rModelPart.Nodes();

        // Auxiliar values
        array_1d<double, 3> force_residual = ZeroVector(3);
        double flux_residual = 0.0;

        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        const unsigned int dim = r_current_process_info[DOMAIN_SIZE];

        // Getting
        const auto it_node_begin = r_nodes.begin();

        // Iterating nodes
        #pragma omp parallel for firstprivate(force_residual,flux_residual), schedule(guided,512)
        for(int i=0; i<static_cast<int>(r_nodes.size()); ++i) {
            auto it_node = it_node_begin + i;

            noalias(force_residual) = it_node->FastGetSolutionStepValue(FORCE_RESIDUAL);
            flux_residual = it_node->FastGetSolutionStepValue(FLUX_RESIDUAL);

            if( it_node->IsFixed(DISPLACEMENT_X) == true ) {
                double& r_reaction = it_node->FastGetSolutionStepValue(REACTION_X);
                r_reaction = -force_residual[0];
            }
            if( it_node->IsFixed(DISPLACEMENT_Y) == true ) {
                double& r_reaction = it_node->FastGetSolutionStepValue(REACTION_Y);
                r_reaction = -force_residual[1];
            }
            if( it_node->IsFixed(LIQUID_PRESSURE) == true ) {
                double& r_reaction = it_node->FastGetSolutionStepValue(REACTION_LIQUID_PRESSURE);
                r_reaction = -flux_residual;
            }
            if(dim==3) {
                if( it_node->IsFixed(DISPLACEMENT_Z) == true ) {
                    double& r_reaction = it_node->FastGetSolutionStepValue(REACTION_Z);
                    r_reaction = -force_residual[2];
                }
            }
        }
    }

    virtual void CheckConvergence(ModelPart& rModelPart)
    {
        const ProcessInfo& r_current_process_info = rModelPart.GetProcessInfo();
        NodesArrayType& r_nodes = rModelPart.Nodes();
        const auto it_node_begin = rModelPart.NodesBegin();

        double l2_numerator = 0.0;
        double l2_denominator = 0.0;
        #pragma omp parallel for reduction(+:l2_numerator,l2_denominator)
        for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
            auto itCurrentNode = it_node_begin + i;
            const array_1d<double, 3>& r_current_displacement = itCurrentNode->FastGetSolutionStepValue(DISPLACEMENT);
            const array_1d<double, 3>& r_previous_displacement = itCurrentNode->FastGetSolutionStepValue(DISPLACEMENT,1);
            const double& r_current_liquid_pressure = itCurrentNode->FastGetSolutionStepValue(LIQUID_PRESSURE);
            const double& r_previous_liquid_pressure = itCurrentNode->FastGetSolutionStepValue(LIQUID_PRESSURE,1);
            array_1d<double, 3> delta_displacement;
            noalias(delta_displacement) = r_current_displacement - r_previous_displacement;
            const double delta_liquid_pressure = r_current_liquid_pressure - r_previous_liquid_pressure;
            const double norm_2_du = inner_prod(delta_displacement,delta_displacement) + delta_liquid_pressure*delta_liquid_pressure;
            const double norm_2_u_old = inner_prod(r_previous_displacement,r_previous_displacement) + r_previous_liquid_pressure*r_previous_liquid_pressure;

            l2_numerator += norm_2_du;
            l2_denominator += norm_2_u_old;
        }
        if (l2_denominator > 1.0e-12) {
            double l2_abs_error = std::sqrt(l2_numerator);
            double l2_rel_error = l2_abs_error/std::sqrt(l2_denominator);

            // std::fstream l2_error_file;
            // l2_error_file.open ("l2_error_time.txt", std::fstream::out | std::fstream::app);
            // l2_error_file.precision(12);
            // l2_error_file << r_current_process_info[TIME] << " " << l2_rel_error << std::endl;
            // l2_error_file.close();

            if (l2_rel_error <= r_current_process_info[ERROR_RATIO] && l2_abs_error <= r_current_process_info[ERROR_INTEGRATION_POINT]) {
                KRATOS_INFO("EXPLICIT CONVERGENCE CHECK") << "The simulation is converging at step: " << r_current_process_info[STEP] << std::endl;
                KRATOS_INFO("EXPLICIT CONVERGENCE CHECK") << "L2 Relative Error is: " << l2_rel_error << std::endl;
                KRATOS_INFO("EXPLICIT CONVERGENCE CHECK") << "L2 Absolute Error is: " << l2_abs_error << std::endl;
            }
        }
    }

}; // Class PoromechanicsExplicitStrategy

} // namespace Kratos

#endif // KRATOS_POROMECHANICS_EXPLICIT_STRATEGY  defined
