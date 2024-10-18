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

// External includes

// Project includes
#include "includes/define.h"
#include "processes/find_global_nodal_neighbours_process.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"

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

namespace
{

template<class TDataType, bool TIsHistorical>
struct GradientDataHandler final
{
    ///@name Type Definitions
    ///@{

    using NodeType = typename ModelPart::NodeType;

    static constexpr bool IsDataTypeScalar = std::is_same<TDataType, double>();

    using GradientDataType = typename std::conditional<IsDataTypeScalar, array_1d<double,3>, Matrix>::type;

    ///@}
    ///@name Operations
    ///@{

    static double GetOrigin(
        const NodeType& rNode,
        const Variable<TDataType>& rOriginVar,
        const std::size_t Component = 0);

    static void SetGradient(
        NodeType& rNode,
        const Variable<GradientDataType>& rGradientVar,
        const array_1d<double,3>& rGradientValue,
        const std::size_t Component = 0);

    ///@}
};

template<>
double GradientDataHandler<double, true>::GetOrigin(
    const NodeType& rNode,
    const Variable<double>& rOriginVar,
    const std::size_t Component)
{
    return rNode.FastGetSolutionStepValue(rOriginVar);
}

template<>
double GradientDataHandler<double, false>::GetOrigin(
    const NodeType& rNode,
    const Variable<double>& rOriginVar,
    const std::size_t Component)
{
    return rNode.GetValue(rOriginVar);
}

template<>
double GradientDataHandler<array_1d<double,3>, true>::GetOrigin(
    const NodeType& rNode,
    const Variable<array_1d<double,3>>& rOriginVar,
    const std::size_t Component)
{
    return rNode.FastGetSolutionStepValue(rOriginVar)[Component];
}

template<>
double GradientDataHandler<array_1d<double,3>, false>::GetOrigin(
    const NodeType& rNode,
    const Variable<array_1d<double,3>>& rOriginVar,
    const std::size_t Component)
{
    return rNode.GetValue(rOriginVar)[Component];
}

template<>
void GradientDataHandler<double, true>::SetGradient(
    NodeType& rNode,
    const Variable<array_1d<double,3>>& rGradientVar,
    const array_1d<double,3>& rGradientValue,
    const std::size_t Component)
{
    noalias(rNode.FastGetSolutionStepValue(rGradientVar)) = rGradientValue;
}

template<>
void GradientDataHandler<double, false>::SetGradient(
    NodeType& rNode,
    const Variable<array_1d<double,3>>& rGradientVar,
    const array_1d<double,3>& rGradientValue,
    const std::size_t Component)
{
    noalias(rNode.GetValue(rGradientVar)) = rGradientValue;
}

template<>
void GradientDataHandler<array_1d<double,3>, true>::SetGradient(
    NodeType& rNode,
    const Variable<Matrix>& rGradientVar,
    const array_1d<double,3>& rGradientComponentValue,
    const std::size_t Component)
{
    auto &r_grad = rNode.FastGetSolutionStepValue(rGradientVar);
    if (r_grad.size1() != 3 || r_grad.size2() != 3) {
        r_grad = ZeroMatrix(3,3);
    }

    for (std::size_t i = 0; i < 3; ++i) {
        r_grad(Component, i) = rGradientComponentValue[i];
    }
}

template<>
void GradientDataHandler<array_1d<double,3>, false>::SetGradient(
    NodeType& rNode,
    const Variable<Matrix>& rGradientVar,
    const array_1d<double,3>& rGradientComponentValue,
    const std::size_t Component)
{
    auto &r_grad = rNode.GetValue(rGradientVar);
    if (r_grad.size1() != 3 || r_grad.size2() != 3) {
        r_grad = ZeroMatrix(3,3);
    }

    for (std::size_t i = 0; i < 3; ++i) {
        r_grad(Component, i) = rGradientComponentValue[i];
    }
}

}

/**
 * @brief Edge-based gradient recovery process
 * This process implements the edge-based gradient recovery process technique
 * described in https://onlinelibrary.wiley.com/doi/epdf/10.1002/nme.4374.
 * @tparam TDataType
 * @tparam TSparseSpace
 * @tparam TDenseSpace
 * @tparam TLinearSolver
 */
template<class TDataType, class TSparseSpace, class TDenseSpace, class TLinearSolver>
class EdgeBasedGradientRecoveryProcess : public Process
{
public:

    ///@name Type Definitions
    ///@{

    static constexpr bool IsDataTypeScalar = std::is_same<TDataType, double>();

    using GradientDataType = typename std::conditional<IsDataTypeScalar, array_1d<double,3>, Matrix>::type;

    using NodeType = typename ModelPart::NodeType;

    using OriginGetFunctionType = std::function<double(
        const NodeType&,
        const Variable<TDataType>&,
        const std::size_t)>;

    using GradientSetFunctionType = std::function<void(
        NodeType &rNode,
        const Variable<GradientDataType>&,
        const array_1d<double, 3>&,
        const std::size_t)>;

    using SolvingStrategyType = ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>;

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
        : EdgeBasedGradientRecoveryProcess(rModel, ThisParameters)
    {
        auto p_builder_solver = Kratos::make_shared< ResidualBasedBlockBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>>(pLinearSolver);
        InitializeGradientRecoveryStrategy(p_builder_solver);
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
        if (mSettings["is_historical_origin_variable"].GetBool()) {
            VariableUtils().CheckVariableExists<Variable<TDataType>>(*mpOriginVar, mpOriginModelPart->Nodes());
        }
        if (mSettings["is_historical_gradient_variable"].GetBool()) {
            VariableUtils().CheckVariableExists<Variable<GradientDataType>>(*mpGradientVar, mpOriginModelPart->Nodes());
        }

        return 0;
    }

    void Execute() override
    {
        KRATOS_TRY;

        // Do the gradient recovery procedure
        // Note that we intentionally avoid calling the Check() and Finalize() methods in case
        // the Execute() is called many times as these two perform expensive operations.
        ExecuteInitialize();
        ExecuteInitializeSolutionStep();
        ExecuteFinalizeSolutionStep();

        KRATOS_CATCH("")
    }

    void ExecuteInitialize() override
    {
        // Fill the auxiliary convection model part if not done yet
        if (!mModelPartIsInitialized) {
            InitializeGradientRecoveryModelPart(*mpOriginModelPart);
        }

        // Set the gradient penalty coefficient in the gradient model part ProcessInfo container
        auto& r_process_info = mpGradientRecoveryModelPart->GetProcessInfo();
        r_process_info.SetValue(GRADIENT_PENALTY_COEFFICIENT, mSettings["gradient_penalty_coefficient"].GetDouble());

        // Initialize gradient variables in the non-historical case
        // Note that does the threadsafe memory allocation
        if (!mSettings["is_historical_gradient_variable"].GetBool()) {
            VariableUtils().SetNonHistoricalVariableToZero(*mpGradientVar, mpOriginModelPart->Nodes());
        }
    }

    void ExecuteInitializeSolutionStep() override
    {
        // Set the getter functions according to the origin and gradient variables databases
        OriginGetFunctionType origin_value_getter;
        if (mSettings["is_historical_origin_variable"].GetBool()) {
            origin_value_getter = GradientDataHandler<TDataType, true>::GetOrigin;
        } else {
            origin_value_getter = GradientDataHandler<TDataType, false>::GetOrigin;
        }

        GradientSetFunctionType gradient_value_setter;
        if (mSettings["is_historical_gradient_variable"].GetBool()) {
            gradient_value_setter = GradientDataHandler<TDataType, true>::SetGradient;
        } else {
            gradient_value_setter = GradientDataHandler<TDataType, false>::SetGradient;
        }

        const std::size_t n_components = IsDataTypeScalar ? 1 : mpOriginModelPart->GetProcessInfo().GetValue(DOMAIN_SIZE);
        for (std::size_t d = 0; d < n_components; ++d) {
            CalculateGradientComponent(origin_value_getter, gradient_value_setter, d);
        }
    }

    void ExecuteFinalizeSolutionStep() override
    {
        if (mSettings["reform_gradient_model_part_at_each_step"].GetBool()) {
            Clear();
        }
    }

    void ExecuteFinalize() override
    {
        Clear();
    }

    void Clear() override
    {
        // Empty the model part
        // Note that we call the Clear to avoid emptying the ProcessInfo
        mpGradientRecoveryModelPart->Clear();
        mModelPartIsInitialized = false;

        // Clear the linear strategy
        mpSolvingStrategy->Clear();
    }

    const Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters = Parameters(R"({
            "echo_level" : 0,
            "model_part_name" : "",
            "gradient_recovery_model_part_name" : "",
            "origin_variable" : "",
            "gradient_variable" : "",
            "gradient_penalty_coefficient" : 1.0e-6,
            "calculate_nodal_neighbours" : true,
            "is_historical_origin_variable" : true,
            "is_historical_gradient_variable": true,
            "reform_gradient_model_part_at_each_step" : false
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

    Model& mrModel;

    ModelPart* mpOriginModelPart;

    ModelPart* mpGradientRecoveryModelPart = nullptr;

    const Variable<TDataType>* mpOriginVar = nullptr;

    const Variable<GradientDataType>* mpGradientVar = nullptr;

    bool mModelPartIsInitialized = false;

    typename SolvingStrategyType::UniquePointer mpSolvingStrategy;

    std::string mGradientModelPartName;

    std::string mElementRegisterName;

    Parameters mSettings;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    EdgeBasedGradientRecoveryProcess(
        Model& rModel,
        Parameters ThisParameters)
        : Process()
        , mrModel(rModel)
    {
        // Validate the common settings as well as the element formulation specific ones
        ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());

        // Checks and assign all the required member variables
        CheckAndAssignSettings(ThisParameters);

        // Save validated parameters
        mSettings = ThisParameters;
    }

    virtual void InitializeGradientRecoveryModelPart(ModelPart& rOriginModelPart)
    {
        KRATOS_TRY

        // Set the model part for the gradient recovery
        if (mrModel.HasModelPart(mGradientModelPartName)) {
            mrModel.DeleteModelPart(mGradientModelPartName);
        }
        mpGradientRecoveryModelPart = &(mrModel.CreateModelPart(mGradientModelPartName));

        // Add NODAL_VAUX variable to calculate the gradient with it
        // Note that NODAL_MAUX is not used as the element retreives it from the non-historical database
        mpGradientRecoveryModelPart->AddNodalSolutionStepVariable(NODAL_VAUX);

        // Emulate the origin model part nodes in the gradient model part
        auto& r_gradient_mp = *mpGradientRecoveryModelPart;
        for (const auto& r_orig_node : rOriginModelPart.Nodes()) {
            auto p_new_node = r_gradient_mp.CreateNewNode(r_orig_node.Id(), r_orig_node);
        }

        // Ensure that the nodes have the auxiliary gradient variable as a DOF
        VariableUtils().AddDof<Variable<double>>(NODAL_VAUX_X, *mpGradientRecoveryModelPart);
        VariableUtils().AddDof<Variable<double>>(NODAL_VAUX_Y, *mpGradientRecoveryModelPart);
        if (mpOriginModelPart->GetProcessInfo().GetValue(DOMAIN_SIZE) == 3) {
            VariableUtils().AddDof<Variable<double>>(NODAL_VAUX_Z, *mpGradientRecoveryModelPart);
        }

        // Calculate the nodal neighbours in the origin model part
        if (mSettings["calculate_nodal_neighbours"].GetBool()) {
            FindGlobalNodalNeighboursProcess neigh_proc(rOriginModelPart);
            neigh_proc.Execute();
        }

        // Initialize flags in origin model part
        VariableUtils().SetFlag(VISITED, false, rOriginModelPart.Nodes());

        // Generating the edge elements
        // Note that we take advantage of the fact that the neighbour connectivities are the same
        std::size_t id = 0;
        const auto p_prop_0 = r_gradient_mp.CreateNewProperties(0);
        const std::size_t n_nodes = rOriginModelPart.NumberOfNodes();
        for (std::size_t i_node = 0; i_node < n_nodes; ++i_node) {
            // Get origin node neighbours to create the edge elements
            const auto it_node_orig_mp = rOriginModelPart.NodesBegin() + i_node;
            auto& r_orig_neigh = it_node_orig_mp->GetValue(NEIGHBOUR_NODES);
            for (auto& r_neigh : r_orig_neigh) {
                if (!r_neigh.Is(VISITED)) {
                    std::vector<std::size_t> aux_ids = {it_node_orig_mp->Id(), r_neigh.Id()};
                    r_gradient_mp.CreateNewElement(mElementRegisterName, ++id, aux_ids, p_prop_0);
                }
            }

            // Flag the current node as visited to avoid creating the same edge twice
            it_node_orig_mp->Set(VISITED, true);
        }

        // Set the edges model part initialization flag to avoid creating it twice
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
        mSettings = ThisParameters;

        // Set the origin and gradient variables pointers
        mpOriginVar = &KratosComponents<Variable<TDataType>>::Get(ThisParameters["origin_variable"].GetString());
        mpGradientVar = &KratosComponents<Variable<GradientDataType>>::Get(ThisParameters["gradient_variable"].GetString());

        // Check user-defined origin model part
        KRATOS_ERROR_IF(ThisParameters["model_part_name"].GetString() == "") << "Empty 'model_part_name'. This needs to be provided." << std::endl;
        mpOriginModelPart = &(mrModel.GetModelPart(ThisParameters["model_part_name"].GetString()));

        // Set the gradient recovery element name
        KRATOS_ERROR_IF_NOT(mpOriginModelPart->GetProcessInfo().Has(DOMAIN_SIZE)) << "No 'DOMAIN_SIZE' in origin model part ProcessInfo container." << std::endl;
        const std::size_t domain_size = mpOriginModelPart->GetProcessInfo().GetValue(DOMAIN_SIZE);
        mElementRegisterName = "EdgeBasedGradientRecoveryElement" + std::to_string(domain_size) + "D2N";

        // Check origin model part content
        KRATOS_ERROR_IF(mpOriginModelPart->GetCommunicator().GlobalNumberOfNodes() == 0) << "The origin model part has no nodes." << std::endl;
        if(domain_size == 2){
            KRATOS_ERROR_IF(mpOriginModelPart->ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Triangle) <<
                "In 2D the element type is expected to be a triangle. Quadrilateral elements require extra artificial edges (not implemented yet)." << std::endl;
        } else if(domain_size == 3) {
            KRATOS_ERROR_IF(mpOriginModelPart->ElementsBegin()->GetGeometry().GetGeometryFamily() != GeometryData::KratosGeometryFamily::Kratos_Tetrahedra) <<
                "In 3D the element type is expected to be a tetrahedra. Hexahedral elements require extra artificial edges (not implemented yet)." << std::endl;
        }

        // Check and set the gradient model part name
        if (ThisParameters["gradient_recovery_model_part_name"].GetString() == "") {
            mGradientModelPartName = mpOriginModelPart->Name() + "GradientRecoveryPart";
        } else {
            mGradientModelPartName = ThisParameters["gradient_recovery_model_part_name"].GetString();
        }
    }

    void InitializeGradientRecoveryStrategy(BuilderAndSolverPointerType pBuilderAndSolver)
    {
        // Generate an auxilary model part and populate it by elements of type DistanceCalculationElementSimplex
        InitializeGradientRecoveryModelPart(*mpOriginModelPart);

        // Create and initialize the linear strategy
        bool calculate_norm_dx = false;
        bool calculate_reactions = false;
        bool reform_dof_at_each_iteration = false;
        auto p_scheme = Kratos::make_shared<ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>>();
        mpSolvingStrategy = Kratos::make_unique<ResidualBasedLinearStrategy<TSparseSpace,TDenseSpace,TLinearSolver>>(
            *mpGradientRecoveryModelPart,
            p_scheme,
            pBuilderAndSolver,
            calculate_reactions,
            reform_dof_at_each_iteration,
            calculate_norm_dx);
        mpSolvingStrategy->SetEchoLevel(mSettings["echo_level"].GetInt());
        mpSolvingStrategy->Check();
        mpSolvingStrategy->Initialize();
    }

    void CalculateGradientComponent(
        const OriginGetFunctionType& rOriginValueGetter,
        const GradientSetFunctionType& rGradientValueSetter,
        const std::size_t Component)
    {
        // Set the origin variable in the edge-based nodes
        const std::size_t n_nodes = mpOriginModelPart->NumberOfNodes();
        IndexPartition<std::size_t>(n_nodes).for_each([&](std::size_t i_node){
            const auto it_orig_node = mpOriginModelPart->NodesBegin() + i_node;
            auto it_grad_node = mpGradientRecoveryModelPart->NodesBegin() + i_node;
            const double orig_val = rOriginValueGetter(*it_orig_node, *mpOriginVar, Component);
            it_grad_node->SetValue(NODAL_MAUX, orig_val);
        });

        // Solve the gradient recovery global problem
        mpSolvingStrategy->InitializeSolutionStep();
        mpSolvingStrategy->Predict();
        mpSolvingStrategy->SolveSolutionStep();
        mpSolvingStrategy->FinalizeSolutionStep();

        // Transfer the gradient solution to the origin mesh
        IndexPartition<std::size_t>(n_nodes).for_each([&](std::size_t i_node){
            auto it_orig_node = mpOriginModelPart->NodesBegin() + i_node;
            const auto it_grad_node = mpGradientRecoveryModelPart->NodesBegin() + i_node;
            const array_1d<double,3>& grad_val = it_grad_node->FastGetSolutionStepValue(NODAL_VAUX);
            rGradientValueSetter(*it_orig_node, *mpGradientVar, grad_val, Component);
        });
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
template <class TDataType, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::istream& operator>>(
    std::istream& rIStream,
    EdgeBasedGradientRecoveryProcess<TDataType, TSparseSpace, TDenseSpace, TLinearSolver>& rThis);

/// Output stream function
template <class TDataType, class TSparseSpace, class TDenseSpace, class TLinearSolver>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const EdgeBasedGradientRecoveryProcess<TDataType, TSparseSpace, TDenseSpace, TLinearSolver>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_EDGE_BASED_GRADIENT_RECOVERY_PROCESS_INCLUDED  defined
