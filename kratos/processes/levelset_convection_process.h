//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Ruben Zorrilla
//                   Mohammad Reza Hashemi
//

#if !defined(KRATOS_LEVELSET_CONVECTION_PROCESS_INCLUDED )
#define  KRATOS_LEVELSET_CONVECTION_PROCESS_INCLUDED

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes

// Project includes
#include "includes/convection_diffusion_settings.h"
#include "includes/define.h"
#include "includes/global_pointer_variables.h"
#include "includes/kratos_flags.h"
#include "elements/levelset_convection_element_simplex.h"
#include "elements/levelset_convection_element_simplex_algebraic_stabilization.h"
#include "geometries/geometry_data.h"
#include "processes/compute_nodal_gradient_process.h"
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

/// Short class definition.
/**takes a model part full of SIMPLICIAL ELEMENTS (triangles and tetras) and convects a level set distance
* on the top of it
*/
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver >
class LevelSetConvectionProcess
    : public Process
{
public:

    ///@name Type Definitions
    ///@{

    typedef SolvingStrategy< TSparseSpace, TDenseSpace, TLinearSolver > SolvingStrategyType;
    typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderAndSolverPointerType;
    typedef ComputeNodalGradientProcess<ComputeNodalGradientProcessSettings::SaveAsNonHistoricalVariable> ComputeGradientProcessType;
    typedef ComputeGradientProcessType::Pointer ComputeGradientProcessPointerType;

    ///@}
    ///@name Pointer Definitions
    ///@{

    /// Pointer definition of LevelSetConvectionProcess
    KRATOS_CLASS_POINTER_DEFINITION(LevelSetConvectionProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new Level Set Convection Process object
     * Level set convection proces model constructor
     * @param rModel Model container
     * @param pLinearSolver Linear solver to be used in the level set convection problem
     * @param ThisParameters Json settings encapsulating the process configuration (see also GetDefaultParameters)
     */
    LevelSetConvectionProcess(
        Model& rModel,
        typename TLinearSolver::Pointer pLinearSolver,
        Parameters ThisParameters)
        : LevelSetConvectionProcess(
            rModel.GetModelPart(ThisParameters["model_part_name"].GetString()),
            pLinearSolver,
            ThisParameters)
    {
    }

    /**
     * @brief Construct a new Level Set Convection Process object
     * Level set convection proces model part constructor
     * @param rBaseModelPart Origin model part
     * @param pLinearSolver Linear solver to be used in the level set convection problem
     * @param ThisParameters Json settings encapsulating the process configuration (see also GetDefaultParameters)
     */
    LevelSetConvectionProcess(
        ModelPart& rBaseModelPart,
        typename TLinearSolver::Pointer pLinearSolver,
        Parameters ThisParameters)
        : LevelSetConvectionProcess(
            rBaseModelPart,
            ThisParameters)
    {
        KRATOS_TRY

        auto p_builder_solver = Kratos::make_shared< ResidualBasedBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>>(pLinearSolver);
        InitializeConvectionStrategy(p_builder_solver);

        KRATOS_CATCH("")
    }

    /// Copy constructor.
    LevelSetConvectionProcess(LevelSetConvectionProcess const& rOther) = delete;

    /// Destructor.
    ~LevelSetConvectionProcess() override
    {
        mrModel.DeleteModelPart(mAuxModelPartName);
    }

    ///@}
    ///@name Operators
    ///@{

    void operator()(){
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

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
        if(mDistancePartIsInitialized == false){
            ReGenerateConvectionModelPart(mrBaseModelPart);
        }

        // Evaluate steps needed to achieve target max_cfl
        const auto n_substep = EvaluateNumberOfSubsteps();

        auto& rCurrentProcessInfo = mpDistanceModelPart->GetProcessInfo();

        // Save the variables to be employed so that they can be restored after the solution
        const auto& r_previous_var = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS)->GetUnknownVariable();
        const double previous_delta_time = rCurrentProcessInfo.GetValue(DELTA_TIME);

        // Save current level set value and current and previous step velocity values
        IndexPartition<int>(mpDistanceModelPart->NumberOfNodes()).for_each(
        [&](int i_node){
            const auto it_node = mpDistanceModelPart->NodesBegin() + i_node;
            mVelocity[i_node] = it_node->FastGetSolutionStepValue(*mpConvectVar);
            mVelocityOld[i_node] = it_node->FastGetSolutionStepValue(*mpConvectVar,1);
            mOldDistance[i_node] = it_node->FastGetSolutionStepValue(*mpLevelSetVar,1);
            }
        );

        const double dt =  previous_delta_time / static_cast<double>(n_substep);
        rCurrentProcessInfo.SetValue(DELTA_TIME, dt);
        rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS)->SetUnknownVariable(*mpLevelSetVar);

        for(unsigned int step = 1; step <= n_substep; ++step){
            KRATOS_INFO_IF("LevelSetConvectionProcess", mpSolvingStrategy->GetEchoLevel() > 0) <<
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
            IndexPartition<int>(mpDistanceModelPart->NumberOfNodes()).for_each(
            [&](int i_node){
                auto it_node = mpDistanceModelPart->NodesBegin() + i_node;

                const array_1d<double,3>& r_v = mVelocity[i_node];
                const array_1d<double,3>& r_v_old = mVelocityOld[i_node];

                noalias(it_node->FastGetSolutionStepValue(*mpConvectVar)) = Nold * r_v_old + Nnew * r_v;
                noalias(it_node->FastGetSolutionStepValue(*mpConvectVar, 1)) = Nold_before * r_v_old + Nnew_before * r_v;
                it_node->FastGetSolutionStepValue(*mpLevelSetVar, 1) = it_node->FastGetSolutionStepValue(*mpLevelSetVar);
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
        IndexPartition<int>(mpDistanceModelPart->NumberOfNodes()).for_each(
        [&](int i_node){
            auto it_node = mpDistanceModelPart->NodesBegin() + i_node;
            it_node->FastGetSolutionStepValue(*mpConvectVar) = mVelocity[i_node];
            it_node->FastGetSolutionStepValue(*mpConvectVar,1) = mVelocityOld[i_node];
            it_node->FastGetSolutionStepValue(*mpLevelSetVar,1) = mOldDistance[i_node];
        }
        );

        KRATOS_CATCH("")
    }

    void Clear() override
    {
        mpDistanceModelPart->Nodes().clear();
        mpDistanceModelPart->Conditions().clear();
        mpDistanceModelPart->Elements().clear();
        // mpDistanceModelPart->GetProcessInfo().clear();
        mDistancePartIsInitialized = false;

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
            "convection_model_part_name" : "",
            "levelset_variable_name" : "DISTANCE",
            "levelset_convection_variable_name" : "VELOCITY",
            "levelset_gradient_variable_name" : "DISTANCE_GRADIENT",
            "max_CFL" : 1.0,
            "max_substeps" : 0,
            "eulerian_error_compensation" : false,
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
        return "LevelSetConvectionProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {
        rOStream << "LevelSetConvectionProcess";
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

    ModelPart* mpDistanceModelPart = nullptr;

    const Variable<double>* mpLevelSetVar = nullptr;

    const Variable<array_1d<double,3>>* mpConvectVar = nullptr;

    const Variable<array_1d<double,3>>* mpLevelSetGradientVar = nullptr;

    double mMaxAllowedCFL = 1.0;

	unsigned int mMaxSubsteps = 0;

    bool mIsBfecc;

    bool mElementRequiresLimiter;

    bool mElementRequiresLevelSetGradient;

    bool mEvaluateLimiter;

    Vector mError;

    Vector mOldDistance;

    Vector mSigmaPlus;

    Vector mSigmaMinus;

    Vector mLimiter;

    std::vector<array_1d<double,3>> mVelocity;

    std::vector<array_1d<double,3>> mVelocityOld;

    bool mDistancePartIsInitialized = false;

    typename SolvingStrategyType::UniquePointer mpSolvingStrategy;

    std::string mAuxModelPartName;

    std::string mConvectionElementType;

    const Element* mpConvectionFactoryElement = nullptr;

    Parameters mLevelSetConvectionSettings;

    ComputeGradientProcessPointerType mpGradientCalculator = nullptr;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    LevelSetConvectionProcess(
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

        // Sets the convection diffusion problem settings
        SetConvectionProblemSettings();

        if (mIsBfecc || mElementRequiresLevelSetGradient){
            mpGradientCalculator = Kratos::make_unique<ComputeGradientProcessType>(
            mrBaseModelPart,
            *mpLevelSetVar,
            *mpLevelSetGradientVar,
            NODAL_AREA,
            false);
        }
    }

    /**
     * @brief Set the level set convection formulation settings
     * This method sets the convection diffusion settings specifying the variable to be convect, its gradient, and the convection variable
     * Additionally, it also sets the required ProcessInfo variables
     */
    void SetConvectionProblemSettings()
    {
        // Get the base model part process info
        // Note that this will be shared with the auxiliary model part used in the convection resolution
        auto& r_process_info = mrBaseModelPart.GetProcessInfo();

        // Allocate if needed the variable CONVECTION_DIFFUSION_SETTINGS of the process info, and create it if it does not exist
        if(!r_process_info.Has(CONVECTION_DIFFUSION_SETTINGS)){
            auto p_conv_diff_settings = Kratos::make_shared<ConvectionDiffusionSettings>();
            r_process_info.SetValue(CONVECTION_DIFFUSION_SETTINGS, p_conv_diff_settings);
            p_conv_diff_settings->SetUnknownVariable(*mpLevelSetVar);
            p_conv_diff_settings->SetConvectionVariable(*mpConvectVar);
            p_conv_diff_settings->SetGradientVariable(*mpLevelSetGradientVar);
        }

        // This call returns a function pointer with the ProcessInfo filling directives
        // If the user-defined level set convection requires nothing to be set, the function does nothing
        auto fill_process_info_function = GetFillProcessInfoFunction();
        fill_process_info_function(mrBaseModelPart);
    }

    virtual void ReGenerateConvectionModelPart(ModelPart& rBaseModelPart){

        KRATOS_TRY

        if (mrModel.HasModelPart(mAuxModelPartName)) {
            mrModel.DeleteModelPart(mAuxModelPartName);
        }

        mpDistanceModelPart= &(mrModel.CreateModelPart(mAuxModelPartName));


        // Check buffer size
        const auto base_buffer_size = rBaseModelPart.GetBufferSize();
        KRATOS_ERROR_IF(base_buffer_size < 2) <<
            "Base model part buffer size is " << base_buffer_size << ". Set it to a minimum value of 2." << std::endl;

        // Generate
        mpDistanceModelPart->Nodes().clear();
        mpDistanceModelPart->Conditions().clear();
        mpDistanceModelPart->Elements().clear();

        mpDistanceModelPart->SetProcessInfo(rBaseModelPart.pGetProcessInfo());
        mpDistanceModelPart->SetBufferSize(base_buffer_size);
        for(auto it_properties = rBaseModelPart.PropertiesBegin() ; it_properties != rBaseModelPart.PropertiesEnd() ; ++it_properties){
            mpDistanceModelPart->AddProperties(*(it_properties).base());
        }
        mpDistanceModelPart->Tables() = rBaseModelPart.Tables();

        // Assigning the nodes to the new model part
        mpDistanceModelPart->Nodes() = rBaseModelPart.Nodes();

        // Ensure that the nodes have distance as a DOF
        VariableUtils().AddDof<Variable<double>>(*mpLevelSetVar, rBaseModelPart);

        // Generating the elements
        mpDistanceModelPart->Elements().reserve(rBaseModelPart.NumberOfElements());
        KRATOS_ERROR_IF(mpConvectionFactoryElement == nullptr) << "Convection factory element has not been set yet." << std::endl;
        for (auto it_elem = rBaseModelPart.ElementsBegin(); it_elem != rBaseModelPart.ElementsEnd(); ++it_elem){
            // Create the new element from the factory registered one
            auto p_element = mpConvectionFactoryElement->Create(
                it_elem->Id(),
                it_elem->pGetGeometry(),
                it_elem->pGetProperties());

            mpDistanceModelPart->Elements().push_back(p_element);
        }

        // Initialize the nodal and elemental databases
        InitializeDistanceModelPartDatabases();

        // Resize the arrays
        const auto n_nodes = mpDistanceModelPart->NumberOfNodes();
        mVelocity.resize(n_nodes);
        mVelocityOld.resize(n_nodes);
        mOldDistance.resize(n_nodes);

        if (mEvaluateLimiter){
            mSigmaPlus.resize(n_nodes);
            mSigmaMinus.resize(n_nodes);
            mLimiter.resize(n_nodes);
        }

        if (mIsBfecc){
            mError.resize(n_nodes);
        }

        mDistancePartIsInitialized = true;

        KRATOS_CATCH("")
    }

    /**
     * @brief Initializes the databases values
     * This function initializes is intended to collect all the database initializations
     */
    void InitializeDistanceModelPartDatabases()
    {
        // If required, initialize the limiter elemental and nodal databases
        const array_1d<double, 3> aux_zero_vector = ZeroVector(3);
        if (mConvectionElementType == "LevelSetConvectionElementSimplexAlgebraicStabilization") {
                block_for_each(mpDistanceModelPart->Elements(), [&](Element &rElement) {rElement.SetValue(LIMITER_COEFFICIENT, 0.0);});
        }

        if (mElementRequiresLimiter){
                block_for_each(mpDistanceModelPart->Nodes(), [&](Node<3>& rNode){rNode.SetValue(LIMITER_COEFFICIENT, 0.0);});
        }
    }

    unsigned int EvaluateNumberOfSubsteps(){
        // First of all compute the maximum local CFL number
        const double dt = mpDistanceModelPart->GetProcessInfo()[DELTA_TIME];
        double max_cfl_found = block_for_each<MaxReduction<double>>(mpDistanceModelPart->Elements(), [&](Element& rElement){
            double vol;
            array_1d<double, TDim+1 > N;
            BoundedMatrix<double, TDim+1, TDim > DN_DX;
            auto& r_geom = rElement.GetGeometry();
            GeometryUtils::CalculateGeometryData(r_geom, DN_DX, N, vol);

            // Compute h
            double h=0.0;
            for(unsigned int i=0; i<TDim+1; i++){
                double h_inv = 0.0;
                for(unsigned int k=0; k<TDim; k++){
                    h_inv += DN_DX(i,k)*DN_DX(i,k);
                }
                h += 1.0/h_inv;
            }
            h = sqrt(h)/static_cast<double>(TDim+1);

            // Get average velocity at the nodes
            array_1d<double, 3 > vgauss = ZeroVector(3);
            for(unsigned int i=0; i<TDim+1; i++){
                vgauss += N[i]* r_geom[i].FastGetSolutionStepValue(*mpConvectVar);
            }

            double cfl_local = norm_2(vgauss) / h;
            return cfl_local;
        });
        max_cfl_found *= dt;

        // Synchronize maximum CFL between processes
        max_cfl_found = mpDistanceModelPart->GetCommunicator().GetDataCommunicator().MaxAll(max_cfl_found);

        unsigned int n_steps = static_cast<unsigned int>(max_cfl_found / mMaxAllowedCFL);
        if(n_steps < 1){
            n_steps = 1;
        }

		// Now we compare with the maximum set
		if (mMaxSubsteps > 0 && mMaxSubsteps < n_steps){
            n_steps = mMaxSubsteps;
        }

        return n_steps;
    }

    /**
     * @brief Convection limiter evaluation
     * This function implements the limiter evaluation
     * Note that both the standard and the high order limiter (with nodal projections contributions) are implemented
     */
    void EvaluateLimiter()
    {
        const double epsilon = 1.0e-15;
        const double power_bfecc = 2.0;
        const double power_elemental_limiter = 4.0;

        auto& r_default_comm = mpDistanceModelPart->GetCommunicator().GetDataCommunicator();
        GlobalPointersVector< Node<3 > > gp_list;

        for (int i_node = 0; i_node < static_cast<int>(mpDistanceModelPart->NumberOfNodes()); ++i_node){
            auto it_node = mpDistanceModelPart->NodesBegin() + i_node;
            GlobalPointersVector< Node<3 > >& global_pointer_list = it_node->GetValue(NEIGHBOUR_NODES);

            for (unsigned int j = 0; j< global_pointer_list.size(); ++j)
            {
                auto& global_pointer = global_pointer_list(j);
                gp_list.push_back(global_pointer);
            }
        }

        GlobalPointerCommunicator< Node<3 > > pointer_comm(r_default_comm, gp_list);

        auto coordinate_proxy = pointer_comm.Apply(
            [](GlobalPointer<Node<3> >& global_pointer) -> Point::CoordinatesArrayType
            {
                return global_pointer->Coordinates();
            }
        );

        IndexPartition<int>(mpDistanceModelPart->NumberOfNodes()).for_each(
        [&](int i_node){
            auto it_node = mpDistanceModelPart->NodesBegin() + i_node;

            it_node->SetValue(*mpLevelSetVar, it_node->FastGetSolutionStepValue(*mpLevelSetVar)); //Store mpLevelSetVar

            const auto& X_i = it_node->Coordinates();
            const auto& grad_i = it_node->GetValue(*mpLevelSetGradientVar);

            double S_plus = 0.0;
            double S_minus = 0.0;

            GlobalPointersVector< Node<3 > >& global_pointer_list = it_node->GetValue(NEIGHBOUR_NODES);

            for (unsigned int j = 0; j< global_pointer_list.size(); ++j)
            {
                // if (it_node->Id() == j_node->Id())
                //     continue;

                auto& global_pointer = global_pointer_list(j);
                auto X_j = coordinate_proxy.Get(global_pointer);

                S_plus += std::max(0.0, inner_prod(grad_i, X_i-X_j));
                S_minus += std::min(0.0, inner_prod(grad_i, X_i-X_j));
            }

            mSigmaPlus[i_node] = std::min(1.0, (std::abs(S_minus)+epsilon)/(S_plus+epsilon));
            mSigmaMinus[i_node] = std::min(1.0, (S_plus+epsilon)/(std::abs(S_minus)+epsilon));
        }
        );

        auto combined_proxy = pointer_comm.Apply(
            [&](GlobalPointer<Node<3>> &global_pointer) -> std::pair<double, array_1d<double,3>> {
                return std::make_pair(
                    global_pointer->FastGetSolutionStepValue(*mpLevelSetVar),
                    global_pointer->Coordinates());
            });

        //Calculating beta_ij in a way that the linearity is preserved on non-symmetrical meshes
        IndexPartition<int>(mpDistanceModelPart->NumberOfNodes()).for_each(
        [&](int i_node){
            auto it_node = mpDistanceModelPart->NodesBegin() + i_node;
            const double distance_i = it_node->FastGetSolutionStepValue(*mpLevelSetVar);
            const auto& X_i = it_node->Coordinates();
            const auto& grad_i = it_node->GetValue(*mpLevelSetGradientVar);

            double numerator = 0.0;
            double denominator = 0.0;

            GlobalPointersVector< Node<3 > >& global_pointer_list = it_node->GetValue(NEIGHBOUR_NODES);

            for (unsigned int j = 0; j< global_pointer_list.size(); ++j)
            {
                auto& global_pointer = global_pointer_list(j);
                auto result = combined_proxy.Get(global_pointer);
                const auto X_j = std::get<1>(result);
                const double distance_j = std::get<0>(result);

                double beta_ij = 1.0;

                if (inner_prod(grad_i, X_i-X_j) > 0)
                    beta_ij = mSigmaPlus[i_node];
                else if (inner_prod(grad_i, X_i-X_j) < 0)
                    beta_ij = mSigmaMinus[i_node];

                numerator += beta_ij*(distance_i - distance_j);
                denominator += beta_ij*std::abs(distance_i - distance_j);
            }

            const double fraction = std::abs(numerator) / (denominator + epsilon);

            if (mIsBfecc){
                mLimiter[i_node] = 1.0 - std::pow(fraction, power_bfecc);
            }

            if (mElementRequiresLimiter){
                it_node->GetValue(LIMITER_COEFFICIENT) = (1.0 - std::pow(fraction, power_elemental_limiter));
            }
        }
        );

        if (mElementRequiresLimiter){
            block_for_each(mpDistanceModelPart->Elements(), [&](Element& rElement){
                const auto& r_geometry = rElement.GetGeometry();
                double elemental_limiter = 1.0;
                for(unsigned int i_node=0; i_node< TDim+1; ++i_node) {
                    elemental_limiter = std::min(r_geometry[i_node].GetValue(LIMITER_COEFFICIENT), elemental_limiter);
                }
                rElement.GetValue(LIMITER_COEFFICIENT) = elemental_limiter;
            }
            );
        }
    }

    /**
     * @brief Eulerian error calculation and correction
     * This function implements the Backward Forward Error Compensation and Correction (BFECC) algorithm
     * Note that this assumes that the first forward convection to n+1 has been completed. Then we go backwards
     * to n* to calculate and apply the convection error.
     */
    void ErrorCalculationAndCorrection()
    {
        block_for_each(mpDistanceModelPart->Nodes(), [this](Node<3>& rNode){
            noalias(rNode.FastGetSolutionStepValue(*mpConvectVar)) = -1.0 * rNode.FastGetSolutionStepValue(*mpConvectVar);
            noalias(rNode.FastGetSolutionStepValue(*mpConvectVar, 1)) = -1.0 * rNode.FastGetSolutionStepValue(*mpConvectVar, 1);
            rNode.FastGetSolutionStepValue(*mpLevelSetVar, 1) = rNode.FastGetSolutionStepValue(*mpLevelSetVar);
        });

        mpSolvingStrategy->InitializeSolutionStep();
        mpSolvingStrategy->Predict();
        mpSolvingStrategy->SolveSolutionStep(); // backward convetion to obtain phi_n*
        mpSolvingStrategy->FinalizeSolutionStep();

        // Calculating the raw error without a limiter, etc.
        IndexPartition<int>(mpDistanceModelPart->NumberOfNodes()).for_each(
        [&](int i_node){
            auto it_node = mpDistanceModelPart->NodesBegin() + i_node;
            mError[i_node] =
                0.5*(it_node->GetValue(*mpLevelSetVar) - it_node->FastGetSolutionStepValue(*mpLevelSetVar));
        });

        IndexPartition<int>(mpDistanceModelPart->NumberOfNodes()).for_each(
        [&](int i_node){
            auto it_node = mpDistanceModelPart->NodesBegin() + i_node;
            noalias(it_node->FastGetSolutionStepValue(*mpConvectVar)) = -1.0 * it_node->FastGetSolutionStepValue(*mpConvectVar);
            noalias(it_node->FastGetSolutionStepValue(*mpConvectVar, 1)) = -1.0 * it_node->FastGetSolutionStepValue(*mpConvectVar, 1);
            const double phi_n_star = it_node->GetValue(*mpLevelSetVar) + mLimiter[i_node]*mError[i_node];
            it_node->FastGetSolutionStepValue(*mpLevelSetVar) = phi_n_star;
            it_node->FastGetSolutionStepValue(*mpLevelSetVar, 1) = phi_n_star;
        });

        mpSolvingStrategy->InitializeSolutionStep();
        mpSolvingStrategy->Predict();
        mpSolvingStrategy->SolveSolutionStep(); // forward convection to obtain the corrected phi_n+1
        mpSolvingStrategy->FinalizeSolutionStep();
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
        mElementRequiresLimiter =  ThisParameters["element_settings"].Has("include_anti_diffusivity_terms") ? ThisParameters["element_settings"]["include_anti_diffusivity_terms"].GetBool() : false;
        if (mConvectionElementType == "LevelSetConvectionElementSimplexAlgebraicStabilization"){
            ThisParameters["element_settings"]["requires_distance_gradient"].SetBool(mElementRequiresLimiter);
        }
        mElementRequiresLevelSetGradient = ThisParameters["element_settings"]["requires_distance_gradient"].GetBool();;

        // Convection related settings
        mMaxAllowedCFL = ThisParameters["max_CFL"].GetDouble();
        mMaxSubsteps = ThisParameters["max_substeps"].GetInt();
        mIsBfecc = ThisParameters["eulerian_error_compensation"].GetBool();
        mMaxAllowedCFL = ThisParameters["max_CFL"].GetDouble();
        mpLevelSetVar = &KratosComponents<Variable<double>>::Get(ThisParameters["levelset_variable_name"].GetString());
        mpConvectVar = &KratosComponents<Variable<array_1d<double,3>>>::Get(ThisParameters["levelset_convection_variable_name"].GetString());
        if (ThisParameters["convection_model_part_name"].GetString() == "") {
            mAuxModelPartName = mrBaseModelPart.Name() + "_DistanceConvectionPart";
        } else {
            mAuxModelPartName = ThisParameters["convection_model_part_name"].GetString();
        }

        // Limiter related settings
        mpLevelSetGradientVar = (mIsBfecc || mElementRequiresLevelSetGradient) ? &(KratosComponents<Variable<array_1d<double, 3>>>::Get(ThisParameters["levelset_gradient_variable_name"].GetString())) : nullptr;
        mEvaluateLimiter = (mIsBfecc || mElementRequiresLimiter) ? true : false;
    }

    /**
     * @brief Get the Convection Elements List object
     * This method returns a list with the available formulations for the level set convection
     * @return const std::vector<std::string> List containing the available formulations
     */
    const virtual inline std::vector<std::string> GetConvectionElementsList()
    {
        std::vector<std::string> elements_list = {
            "levelset_convection_supg",
            "levelset_convection_algebraic_stabilization"
        };
        return elements_list;
    }

    /**
     * @brief Get the Convection Element Name object
     * This method maps the user-defined element name to the Kratos class name
     * @param InputName User-defined element name
     * @return const std::string Kratos convection element class name
     */
    const virtual std::string GetConvectionElementName(std::string InputName)
    {
        const std::map<std::string, std::string> elements_name_map {
            {"levelset_convection_supg","LevelSetConvectionElementSimplex"},
            {"levelset_convection_algebraic_stabilization", "LevelSetConvectionElementSimplexAlgebraicStabilization"}
        };
        return elements_name_map.at(InputName);
    }

    /**
     * @brief Get the Convection Element Default Parameters object
     * For each of the available formulations, this method returns the corresponding settings
     * @param ElementType User-defined element type
     * @return const Parameters Json string encapsulating the input element settings
     */
    const virtual Parameters GetConvectionElementDefaultParameters(const std::string ElementType)
    {
        Parameters default_parameters;
        if (ElementType == "levelset_convection_supg") {
            default_parameters = Parameters(R"({
                "dynamic_tau" : 0.0,
                "cross_wind_stabilization_factor" : 0.7,
                "requires_distance_gradient" : false
            })");
        } else if (ElementType == "levelset_convection_algebraic_stabilization") {
            default_parameters = Parameters(R"({
                "include_anti_diffusivity_terms" : false,
                "requires_distance_gradient" : false
            })");
        } else {
            KRATOS_ERROR << "Default parameters are not implemented for the specified \'" << ElementType << "\' element. Available options are \n\t- \'levelset_convection_supg\'\n\t- \'levelset_convection_algebraic_stabilization\'" << std::endl;
        }

        return default_parameters;
    }

    /**
     * @brief Get the Fill Process Info Function object
     * This method returns a lambda function with the required operations to be performed in the process info
     * It has to be particularised for all the formulations. If not particularised a do nothing instruction is returned
     * @return const std::function<void(ModelPart&)> A function pointer to be called when setting up the distance model part
     */
    const virtual std::function<void(ModelPart&)> GetFillProcessInfoFunction()
    {
        std::function<void(ModelPart&)> fill_process_info_function;

        if (mConvectionElementType == "LevelSetConvectionElementSimplex") {
            fill_process_info_function = [this](ModelPart &rModelPart) {
                auto &r_process_info = rModelPart.GetProcessInfo();
                // If not present, set the DYNAMIC_TAU
                if (r_process_info.Has(DYNAMIC_TAU)) {
                    KRATOS_WARNING("LevelSetConvectionProcess") << "ProcessInfo container already has DYNAMIC_TAU. Using the existent one" << r_process_info.GetValue(DYNAMIC_TAU) << " for the level set convection." << std::endl;
                } else {
                    r_process_info.SetValue(DYNAMIC_TAU, mLevelSetConvectionSettings["element_settings"]["dynamic_tau"].GetDouble());
                }
                // Set CROSS_WIND_STABILIZATION_FACTOR
                r_process_info.SetValue(CROSS_WIND_STABILIZATION_FACTOR, mLevelSetConvectionSettings["element_settings"]["cross_wind_stabilization_factor"].GetDouble());
            };
        } else {
            fill_process_info_function = [](ModelPart &rModelPart) {};
        }

        return fill_process_info_function;
    }

    void InitializeConvectionStrategy(BuilderAndSolverPointerType pBuilderAndSolver)
    {
        // Check that there is at least one element and node in the model
        KRATOS_ERROR_IF(mrBaseModelPart.NumberOfNodes() == 0) << "The model has no nodes." << std::endl;
        KRATOS_ERROR_IF(mrBaseModelPart.NumberOfElements() == 0) << "The model has no elements." << std::endl;

        // Check that the level set and convection variables are in the nodal database
        VariableUtils().CheckVariableExists<Variable<double>>(*mpLevelSetVar, mrBaseModelPart.Nodes());
        VariableUtils().CheckVariableExists<Variable<array_1d<double,3>>>(*mpConvectVar, mrBaseModelPart.Nodes());

        // Check the base model part element family (only simplex elements are supported)
        if(TDim == 2){
            KRATOS_ERROR_IF(mrBaseModelPart.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::Kratos_Triangle) <<
                "In 2D the element type is expected to be a triangle" << std::endl;
        } else if(TDim == 3) {
            KRATOS_ERROR_IF(mrBaseModelPart.ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::Kratos_Tetrahedra) <<
                "In 3D the element type is expected to be a tetrahedra" << std::endl;
        }

        // Generate an auxilary model part and populate it by elements of type DistanceCalculationElementSimplex
        ReGenerateConvectionModelPart(mrBaseModelPart);

        // Generate a linear strategy
        bool CalculateReactions = false;
        bool ReformDofAtEachIteration = false;
        bool CalculateNormDxFlag = false;
        auto p_scheme = Kratos::make_shared<ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>>();
        mpSolvingStrategy = Kratos::make_unique< ResidualBasedLinearStrategy<TSparseSpace,TDenseSpace,TLinearSolver > >(
            *mpDistanceModelPart,
            p_scheme,
            pBuilderAndSolver,
            CalculateReactions,
            ReformDofAtEachIteration,
            CalculateNormDxFlag);

        mpSolvingStrategy->SetEchoLevel(0);
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
    LevelSetConvectionProcess& operator=(LevelSetConvectionProcess const& rOther);

    ///@}
}; // Class LevelSetConvectionProcess

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
    LevelSetConvectionProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>& rThis);

/// Output stream function
template< unsigned int TDim, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const LevelSetConvectionProcess<TDim, TSparseSpace, TDenseSpace, TLinearSolver>& rThis){

    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_LEVELSET_CONVECTION_PROCESS_INCLUDED  defined
