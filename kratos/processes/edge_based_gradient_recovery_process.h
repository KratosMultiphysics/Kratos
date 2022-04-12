//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_EDGE_BASED_GRADIENT_RECOVERY_PROCESS_INCLUDED)
#define  KRATOS_EDGE_BASED_GRADIENT_RECOVERY_PROCESS_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/global_pointer_variables.h"
#include "includes/kratos_flags.h"
#include "elements/levelset_convection_element_simplex.h"
#include "elements/levelset_convection_element_simplex_algebraic_stabilization.h"
#include "geometries/geometry_data.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"
#include "utilities/pointer_communicator.h"
#include "utilities/pointer_map_communicator.h"

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
 * @brief Edge-based gradient recovery process
 * This process implements the edge-based gradient recovery process technique
 * described in https://onlinelibrary.wiley.com/doi/epdf/10.1002/nme.4374.
 * @tparam TDim
 * @tparam TSparseSpace
 * @tparam TDenseSpace
 * @tparam TLinearSolver
 */
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver >
class EdgeBasedGradientRecoveryProcess : public Process
{
public:

    ///@name Type Definitions
    ///@{

    using SolvingStrategyType = ImplicitSolvingStrategy< TSparseSpace, TDenseSpace, TLinearSolver >;
    using BuilderAndSolverPointerType = typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer;

    ///@}
    ///@name Pointer Definitions
    ///@{

    /// Pointer definition of EdgeBasedGradientRecoveryProcess
    KRATOS_CLASS_POINTER_DEFINITION(EdgeBasedGradientRecoveryProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new Edge Based Gradient Recovery Process object
     * Level set convection proces model constructor
     * @param rModel Model container
     * @param pLinearSolver Linear solver to be used in the level set convection problem
     * @param ThisParameters Json settings encapsulating the process configuration (see also GetDefaultParameters)
     */
    EdgeBasedGradientRecoveryProcess(
        Model& rModel,
        typename TLinearSolver::Pointer pLinearSolver,
        Parameters ThisParameters)
        : EdgeBasedGradientRecoveryProcess(
            rModel.GetModelPart(ThisParameters["model_part_name"].GetString()),
            pLinearSolver,
            ThisParameters)
    {
    }

    /**
     * @brief Construct a new Edge Based Gradient Recovery Process object
     * Level set convection proces model part constructor
     * @param rBaseModelPart Origin model part
     * @param pLinearSolver Linear solver to be used in the level set convection problem
     * @param ThisParameters Json settings encapsulating the process configuration (see also GetDefaultParameters)
     */
    EdgeBasedGradientRecoveryProcess(
        ModelPart& rBaseModelPart,
        typename TLinearSolver::Pointer pLinearSolver,
        Parameters ThisParameters)
        : EdgeBasedGradientRecoveryProcess(
            rBaseModelPart,
            ThisParameters)
    {
        KRATOS_TRY

        auto p_builder_solver = Kratos::make_shared< ResidualBasedBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>>(pLinearSolver);
        InitializeGradientRecoveryStrategy(p_builder_solver);

        KRATOS_CATCH("")
    }

    /// Copy constructor.
    EdgeBasedGradientRecoveryProcess(EdgeBasedGradientRecoveryProcess const& rOther) = delete;

    /// Destructor.
    ~EdgeBasedGradientRecoveryProcess() override
    {
        mrModel.DeleteModelPart(mGradientModelPartName);
    }

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    int Check() override
    {
        //FIXME: We should check according to the origin variable type
        // Check that the level set and convection variables are in the nodal database
        VariableUtils().CheckVariableExists<Variable<double>>(*mpOriginVar, mrBaseModelPart.Nodes());
        VariableUtils().CheckVariableExists<Variable<array_1d<double,3>>>(*mpConvectVar, mrBaseModelPart.Nodes());

        return 0;
    }

    /**
     * @brief Perform the level-set convection
     * This solver provides a stabilized convection solver based on [Codina, R., 1993. Comput. Methods Appl. Mech. Engrg., 110(3-4), pp.325-342.]
     * It uses the sub-stepping approach to comply with the user defined maximum CFL number.
     * The error compensation is done according to the BFECC algorithm, which requires forward, backward, and the final forward solution steps (that triplicates the computational cost).
     * The error compensation severely disturbs the monotonicity of the results that is compensated for by implementing a limited BFECC algorithm.
     * The limiter relies on the nodal gradient of LevelSetVar (non-historical variable LevelSetGradientVar). For more info see [Kuzmin et al., Comput. Methods Appl. Mech. Engrg., 322 (2017) 23–41].
     */
    void Execute() override
    {
        KRATOS_TRY;

        // Fill the auxiliary convection model part if not done yet
        if(!mModelPartIsInitialized){
            InitializeGradientRecoveryModelPart(mrBaseModelPart);
        }

        // Evaluate steps needed to achieve target max_cfl
        const auto n_substep = EvaluateNumberOfSubsteps();

        auto& rCurrentProcessInfo = mpGradientRecoveryModelPart->GetProcessInfo();

        // Save the variables to be employed so that they can be restored after the solution
        const auto& r_previous_var = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS)->GetUnknownVariable();
        const double previous_delta_time = rCurrentProcessInfo.GetValue(DELTA_TIME);

        // Save current level set value and current and previous step velocity values
        IndexPartition<int>(mpGradientRecoveryModelPart->NumberOfNodes()).for_each(
        [&](int i_node){
            const auto it_node = mpGradientRecoveryModelPart->NodesBegin() + i_node;
            mVelocity[i_node] = it_node->FastGetSolutionStepValue(*mpConvectVar);
            mVelocityOld[i_node] = it_node->FastGetSolutionStepValue(*mpConvectVar,1);
            mOldDistance[i_node] = it_node->FastGetSolutionStepValue(*mpOriginVar,1);
            }
        );

        const double dt =  previous_delta_time / static_cast<double>(n_substep);
        rCurrentProcessInfo.SetValue(DELTA_TIME, dt);
        rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS)->SetUnknownVariable(*mpOriginVar);

        // We set these values at every time step as other processes/solvers also use them
        auto fill_process_info_function = GetFillProcessInfoFunction();
        fill_process_info_function(mrBaseModelPart);

        for(unsigned int step = 1; step <= n_substep; ++step){
            KRATOS_INFO_IF("EdgeBasedGradientRecoveryProcess", mLevelSetConvectionSettings["echo_level"].GetInt() > 0) <<
                "Doing step "<< step << " of " << n_substep << std::endl;

            if (mIsBfecc || mElementRequiresLevelSetGradient){
                mpGradientCalculator->Execute();
            }

            // Compute shape functions of old and new step
            const double Nold = 1.0 - static_cast<double>(step) / static_cast<double>(n_substep);
            const double Nnew = 1.0 - Nold;

            const double Nold_before = 1.0 - static_cast<double>(step-1) / static_cast<double>(n_substep);
            const double Nnew_before = 1.0 - Nold_before;

            // Emulate clone time step by copying the new distance onto the old one
            IndexPartition<int>(mpGradientRecoveryModelPart->NumberOfNodes()).for_each(
            [&](int i_node){
                auto it_node = mpGradientRecoveryModelPart->NodesBegin() + i_node;

                const array_1d<double,3>& r_v = mVelocity[i_node];
                const array_1d<double,3>& r_v_old = mVelocityOld[i_node];

                noalias(it_node->FastGetSolutionStepValue(*mpConvectVar)) = Nold * r_v_old + Nnew * r_v;
                noalias(it_node->FastGetSolutionStepValue(*mpConvectVar, 1)) = Nold_before * r_v_old + Nnew_before * r_v;
                it_node->FastGetSolutionStepValue(*mpOriginVar, 1) = it_node->FastGetSolutionStepValue(*mpOriginVar);
            }
            );

            // Storing the levelset variable for error calculation and Evaluating the limiter
            if (mEvaluateLimiter) {
                EvaluateLimiter();
            }

            mpSolvingStrategy->InitializeSolutionStep();
            mpSolvingStrategy->Predict();
            mpSolvingStrategy->SolveSolutionStep(); // forward convection to reach phi_n+1
            mpSolvingStrategy->FinalizeSolutionStep();

            // Error Compensation and Correction
            if (mIsBfecc) {
                ErrorCalculationAndCorrection();
            }
        }

        // Reset the processinfo to the original settings
        rCurrentProcessInfo.SetValue(DELTA_TIME, previous_delta_time);
        rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS)->SetUnknownVariable(r_previous_var);

        // Reset the velocities and levelset values to the one saved before the solution process
        IndexPartition<int>(mpGradientRecoveryModelPart->NumberOfNodes()).for_each(
        [&](int i_node){
            auto it_node = mpGradientRecoveryModelPart->NodesBegin() + i_node;
            it_node->FastGetSolutionStepValue(*mpConvectVar) = mVelocity[i_node];
            it_node->FastGetSolutionStepValue(*mpConvectVar,1) = mVelocityOld[i_node];
            it_node->FastGetSolutionStepValue(*mpOriginVar,1) = mOldDistance[i_node];
        }
        );

        KRATOS_CATCH("")
    }

    void Clear() override
    {
        mpGradientRecoveryModelPart->Nodes().clear();
        mpGradientRecoveryModelPart->Conditions().clear();
        mpGradientRecoveryModelPart->Elements().clear();
        // mpGradientRecoveryModelPart->GetProcessInfo().clear();
        mModelPartIsInitialized = false;

        mpSolvingStrategy->Clear();

        mVelocity.clear();
        mVelocityOld.clear();
        mOldDistance.clear();

        mSigmaPlus.clear();
        mSigmaMinus.clear();
        mLimiter.clear();

        mError.clear();
    }

    const Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"({
            "model_part_name" : "",
            "echo_level" : 0,
            "convection_model_part_name" : "",
            "levelset_variable_name" : "DISTANCE",
            "levelset_convection_variable_name" : "VELOCITY",
            "levelset_gradient_variable_name" : "DISTANCE_GRADIENT",
            "max_CFL" : 1.0,
            "max_substeps" : 0,
            "eulerian_error_compensation" : false,
            "BFECC_limiter_acuteness" : 2.0,
            "element_type" : "levelset_convection_supg",
            "element_settings" : {}
        })");

        return default_parameters;
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

    /// Turn back information as a string.
    std::string Info() const override {
        return "EdgeBasedGradientRecoveryProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {
        rOStream << "EdgeBasedGradientRecoveryProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
    }

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

    ModelPart& mrBaseModelPart;

    Model& mrModel;

    ModelPart* mpGradientRecoveryModelPart = nullptr;

    const Variable<double>* mpOriginVar = nullptr;

    const Variable<array_1d<double,3>>* mpConvectVar = nullptr;

    const Variable<array_1d<double,3>>* mpLevelSetGradientVar = nullptr;

    double mMaxAllowedCFL = 1.0;

	unsigned int mMaxSubsteps = 0;

    bool mIsBfecc;

    bool mElementRequiresLimiter;

    bool mElementRequiresLevelSetGradient;

    bool mEvaluateLimiter;

    double mPowerBfeccLimiter = 2.0;

    double mPowerElementalLimiter = 4.0;

    Vector mError;

    Vector mOldDistance;

    Vector mSigmaPlus;

    Vector mSigmaMinus;

    Vector mLimiter;

    std::vector<array_1d<double,3>> mVelocity;

    std::vector<array_1d<double,3>> mVelocityOld;

    bool mModelPartIsInitialized = false;

    typename SolvingStrategyType::UniquePointer mpSolvingStrategy;

    std::string mGradientModelPartName;

    std::string mConvectionElementType;

    const Element* mpConvectionFactoryElement = nullptr;

    Parameters mLevelSetConvectionSettings;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    EdgeBasedGradientRecoveryProcess(
        ModelPart& rModelPart,
        Parameters ThisParameters)
        : mrBaseModelPart(rModelPart)
        , mrModel(rModelPart.GetModel())
    {
        // Validate the common settings as well as the element formulation specific ones
        ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());
        ThisParameters["element_settings"].ValidateAndAssignDefaults(GetConvectionElementDefaultParameters(ThisParameters["element_type"].GetString()));

        // Checks and assign all the required member variables
        CheckAndAssignSettings(ThisParameters);
    }

    virtual void InitializeGradientRecoveryModelPart(ModelPart& rBaseModelPart)
    {
        KRATOS_TRY

        // Set the model part for the gradient recovery
        if (mrModel.HasModelPart(mGradientModelPartName)) {
            mrModel.DeleteModelPart(mGradientModelPartName);
        }
        mpGradientRecoveryModelPart = &(mrModel.CreateModelPart(mGradientModelPartName));

        // Add NODAL_VAUX variable to calculate the gradient with it
        mpGradientRecoveryModelPart->AddNodalSolutionStepVariable(NODAL_VAUX)

        // Emulate the origin model part nodes in the gradient model part
        auto& r_gradient_mp = *mpGradientRecoveryModelPart;
        for (const auto& r_orig_node : rBaseModelPart.Nodes()) {
            r_gradient_mp.CreateNewNode(r_orig_node.Id());
        }

        // Ensure that the nodes have the auxiliary gradient variable as a DOF
        VariableUtils().AddDof<Variable<array_1d<double,3>>>(NODAL_VAUX, rBaseModelPart);

        // Generating the edge elements
        //TODO: Think about some sort of reserve()/shrink_to_fit() for the Elements() container of the model part




        for (auto it_elem = rBaseModelPart.ElementsBegin(); it_elem != rBaseModelPart.ElementsEnd(); ++it_elem){
            // Create the new element from the factory registered one
            auto p_element = mpConvectionFactoryElement->Create(
                it_elem->Id(),
                it_elem->pGetGeometry(),
                it_elem->pGetProperties());

            mpGradientRecoveryModelPart->Elements().push_back(p_element);
        }

        mModelPartIsInitialized = true;

        KRATOS_CATCH("")
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
     * @brief Checks and assign the required member variables
     * This function checks the provided parameters, which need to have been already validated and sets the member variables
     * @param ThisParameters Json string containing the already validated process and formulation settings
     */
    void CheckAndAssignSettings(const Parameters ThisParameters)
    {
        mLevelSetConvectionSettings = ThisParameters;

        // Convection element formulation settings
        std::string element_type = ThisParameters["element_type"].GetString();
        const auto element_list = GetConvectionElementsList();
        KRATOS_ERROR_IF(std::find(element_list.begin(), element_list.end(), element_type) == element_list.end()) << "Specified \'" << element_type << "\' is not in the available elements list." << std::endl;
        mConvectionElementType = GetConvectionElementName(element_type);
        std::string element_register_name = mConvectionElementType + std::to_string(TDim) + "D" + std::to_string(TDim + 1) + "N";
        mpConvectionFactoryElement = &KratosComponents<Element>::Get(element_register_name);

        mpOriginVar = &KratosComponents<Variable<double>>::Get(ThisParameters["levelset_variable_name"].GetString());
        mpConvectVar = &KratosComponents<Variable<array_1d<double,3>>>::Get(ThisParameters["levelset_convection_variable_name"].GetString());
        if (ThisParameters["convection_model_part_name"].GetString() == "") {
            mGradientModelPartName = mrBaseModelPart.Name() + "_DistanceConvectionPart";
        } else {
            mGradientModelPartName = ThisParameters["convection_model_part_name"].GetString();
        }
    }

    void InitializeGradientRecoveryStrategy(BuilderAndSolverPointerType pBuilderAndSolver)
    {
        // Check that there is at least one element and node in the model
        KRATOS_ERROR_IF(mrBaseModelPart.NumberOfNodes() == 0) << "The model has no nodes." << std::endl;

        // Check the base model part element family (only simplex elements are supported)
        if(TDim == 2){
            KRATOS_ERROR_IF(mrBaseModelPart.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Triangle) <<
                "In 2D the element type is expected to be a triangle. Quadrilateral elements require extra artificial edges (not implemented yet)." << std::endl;
        } else if(TDim == 3) {
            KRATOS_ERROR_IF(mrBaseModelPart.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Tetrahedra) <<
                "In 3D the element type is expected to be a tetrahedra. Hexahedral elements require extra artificial edges (not implemented yet)." << std::endl;
        }

        // Generate an auxilary model part and populate it by elements of type DistanceCalculationElementSimplex
        InitializeGradientRecoveryModelPart(mrBaseModelPart);

        // Create and initialize the linear strategy
        bool calculate_norm_dx = false;
        bool calculate_reactions = false;
        bool reform_dof_at_each_iteration = false;
        auto p_scheme = Kratos::make_shared<ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>>();
        mpSolvingStrategy = Kratos::make_unique< ResidualBasedLinearStrategy<TSparseSpace,TDenseSpace,TLinearSolver > >(
            *mpGradientRecoveryModelPart,
            p_scheme,
            pBuilderAndSolver,
            calculate_reactions,
            reform_dof_at_each_iteration,
            calculate_norm_dx);

        mpSolvingStrategy->SetEchoLevel(mLevelSetConvectionSettings["echo_level"].GetInt());
        mpSolvingStrategy->Check();
        mpSolvingStrategy->Initialize();
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

    /// Assignment operator.
    EdgeBasedGradientRecoveryProcess& operator=(EdgeBasedGradientRecoveryProcess const& rOther);

    ///@}
}; // Class EdgeBasedGradientRecoveryProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// Input stream function
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::istream& operator >> (
    std::istream& rIStream,
    EdgeBasedGradientRecoveryProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>& rThis);

/// Output stream function
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const EdgeBasedGradientRecoveryProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>& rThis){

    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_EDGE_BASED_GRADIENT_RECOVERY_PROCESS_INCLUDED  defined
