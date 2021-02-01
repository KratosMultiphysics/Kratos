//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/checks.h"
#include "includes/model_part.h"
#include "includes/parallel_environment.h"
#include "utilities/variable_utils.h"
#include "utilities/divide_triangle_2d_3.h"
#include "utilities/divide_tetrahedra_3d_4.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.h"
#include "modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions.h"
#include "modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions.h"

// Application includes
#include "embedded_skin_visualization_process.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

/* Public functions *******************************************************/
Parameters EmbeddedSkinVisualizationProcess::StaticGetDefaultParameters()
{
    Parameters default_settings(R"(
    {
        "model_part_name"                       : "",
        "visualization_model_part_name"         : "EmbeddedSkinVisualizationModelPart",
        "reform_model_part_at_each_time_step"   : false,
        "level_set_type"                        : "",
        "shape_functions"                       : "",
        "distance_variable_name"                : "",
        "visualization_variables"               : ["VELOCITY","PRESSURE"],
        "visualization_nonhistorical_variables" : []
    })");

    return default_settings;
}

ModelPart& EmbeddedSkinVisualizationProcess::CreateAndPrepareVisualizationModelPart(
    Model& rModel,
    const Parameters rParameters)
{
    // Set visualization model part data
    const std::size_t buffer_size = 1;
    std::string visualization_mp_name = rParameters["visualization_model_part_name"].GetString();
    auto &r_origin_model_part = rModel.GetModelPart(rParameters["model_part_name"].GetString());
    const unsigned int domain_size = r_origin_model_part.GetProcessInfo()[DOMAIN_SIZE];

    // Check visualization model part data
    KRATOS_ERROR_IF_NOT(domain_size == 2 || domain_size == 3) << "Origin model part DOMAIN_SIZE is " << domain_size << "." << std::endl;
    if (visualization_mp_name == "") {
        KRATOS_WARNING("EmbeddedSkinVisualizationProcess") << "\'visualization_model_part_name\' is empty. Using \'EmbeddedSkinVisualizationModelPart\'." << std::endl;
        visualization_mp_name = "EmbeddedSkinVisualizationModelPart";
    }

    // Create the visualization model part
    auto& r_visualization_model_part = rModel.CreateModelPart(visualization_mp_name, buffer_size);

    // Set visualization model part ProcessInfo
    // Note that the visualization model part and the origin model part share the ProcessInfo container
    const auto p_process_info = r_origin_model_part.pGetProcessInfo();
    r_visualization_model_part.SetProcessInfo(p_process_info);

    // Check and add variables to visualization model part
    auto& r_origin_variables_list = r_origin_model_part.GetNodalSolutionStepVariablesList();
    auto& r_visualization_variables_list = r_visualization_model_part.GetNodalSolutionStepVariablesList();
    for (auto var_name_param : rParameters["visualization_variables"]) {
        // Check if variable exists in origin model part
        const std::string var_name = var_name_param.GetString();
        const auto& r_var = KratosComponents<VariableData>::Get(var_name);
        KRATOS_ERROR_IF_NOT(r_origin_variables_list.Has(r_var)) << "Requested variable " << var_name << " is not in the origin model part." << std::endl;
        // Add the variable to the visualization model part
        r_visualization_variables_list.Add(r_var);
    }

    // Set the origin model part as temporary communicator
    // This is required to perform all the IsDistributed checks that appear before creating the visualization mesh
    // Note that these checks might be required outside this process (i.e. in the creation of the visualization mesh output)
    // This will be updated by a proper one by the ParallelFillCommunicator after the creation of the visualization entities
    r_visualization_model_part.SetCommunicator(r_origin_model_part.pGetCommunicator());

    // If MPI, add the PARTITION_INDEX variable to the visualization model part variables
    if (r_visualization_model_part.IsDistributed()) {
        r_visualization_variables_list.Add(PARTITION_INDEX);
    }

    return r_visualization_model_part;
}

void EmbeddedSkinVisualizationProcess::CheckAndSetLevelSetType(
    const Parameters rParameters,
    LevelSetType& rLevelSetType)
{
    const std::string level_set_type = rParameters["level_set_type"].GetString();
    KRATOS_ERROR_IF(level_set_type == "") << "\'level_set_type\' is not prescribed. Admissible values are: \'continuous\' and \'discontinuous\'." << std::endl;

    if (level_set_type == "continuous") {
        rLevelSetType = LevelSetType::Continuous;
    } else if (level_set_type == "discontinuous") {
        rLevelSetType = LevelSetType::Discontinuous;
    } else {
        std::stringstream error_msg;
        error_msg << "Currently prescribed \'level_set_type\': " << level_set_type << std::endl;
        error_msg << "Admissible values are : \'continuous\' and \'discontinuous\'" << std::endl;
        KRATOS_ERROR << error_msg.str();
    }
}

void EmbeddedSkinVisualizationProcess::CheckAndSetShapeFunctionsType(
    const Parameters rParameters,
    ShapeFunctionsType& rShapeFunctionsType)
{
    const std::string shape_functions = rParameters["shape_functions"].GetString();
    KRATOS_ERROR_IF(shape_functions == "") << "\'shape_functions\' is not prescribed. Admissible values are: \'standard\' and \'ausas\'." << std::endl;

    if (shape_functions == "ausas") {
        rShapeFunctionsType = ShapeFunctionsType::Ausas;
    } else if (shape_functions == "standard") {
        rShapeFunctionsType = ShapeFunctionsType::Standard;
    } else {
        std::stringstream error_msg;
        error_msg << "Currently prescribed \'shape_functions\': " << shape_functions << std::endl;
        error_msg << "Admissible values are : \'standard\' and \'ausas\'" << std::endl;
        KRATOS_ERROR << error_msg.str();
    }
}

const std::string EmbeddedSkinVisualizationProcess::CheckAndReturnDistanceVariableName(
    const Parameters rParameters,
    const LevelSetType& rLevelSetType)
{
    std::string distance_variable_name = rParameters["distance_variable_name"].GetString();
    // If the distance variable name is not provided, try to deduce it from the FE space type
    if (distance_variable_name == "") {
        switch (rLevelSetType) {
            // Nodal-based level set
            case LevelSetType::Continuous:
                distance_variable_name = "DISTANCE";
                break;
            // Element-based level set
            case LevelSetType::Discontinuous:
                distance_variable_name = "ELEMENTAL_DISTANCES";
                break;
            default:
                KRATOS_ERROR << "Default \"distance_variable_name\" cannot be deduced from the shape functions type" << std::endl;
        }
        KRATOS_INFO("EmbeddedSkinVisualizationProcess") << "\'distance_variable_name\' is not prescribed. Using default " << distance_variable_name << std::endl;
    // User-defined distance variable name case
    } else {
        KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(distance_variable_name) || KratosComponents<Variable<Vector>>::Has(distance_variable_name))
           << "Provided \"distance_variable_name\" " << distance_variable_name << " is not in the KratosComponents. Please check the provided value." << std::endl;
    }

    return distance_variable_name;
}

template <class TDataType>
void EmbeddedSkinVisualizationProcess::FillVariablesList(
    const Parameters rParameters,
    std::vector<const Variable<TDataType>*>& rVariablesList)
{
    rVariablesList.clear();
    for (auto i_var_params : rParameters) {
        const Variable<TDataType>* p_aux = nullptr;
        const std::string i_var_name = i_var_params.GetString();
        if(KratosComponents<Variable<TDataType>>::Has(i_var_name)){
            const auto& r_var = KratosComponents<Variable<TDataType>>::Get(i_var_name);
            p_aux = &r_var;
            rVariablesList.push_back(p_aux);
        }
    }
}

EmbeddedSkinVisualizationProcess::EmbeddedSkinVisualizationProcess(
    ModelPart& rModelPart,
    ModelPart& rVisualizationModelPart,
    const std::vector<const Variable< double>* >& rVisualizationScalarVariables,
    const std::vector<const Variable< array_1d<double, 3> >* >& rVisualizationVectorVariables,
    const std::vector<const Variable< double>* >& rVisualizationNonHistoricalScalarVariables,
    const std::vector<const Variable< array_1d<double, 3> >* >& rVisualizationNonHistoricalVectorVariables,
    const LevelSetType& rLevelSetType,
    const ShapeFunctionsType& rShapeFunctionsType,
    const bool ReformModelPartAtEachTimeStep) :
    Process(),
    mrModelPart(rModelPart),
    mrVisualizationModelPart(rVisualizationModelPart),
    mLevelSetType(rLevelSetType),
    mShapeFunctionsType(rShapeFunctionsType),
    mReformModelPartAtEachTimeStep(ReformModelPartAtEachTimeStep),
    mVisualizationScalarVariables(rVisualizationScalarVariables),
    mVisualizationVectorVariables(rVisualizationVectorVariables),
    mVisualizationNonHistoricalScalarVariables(rVisualizationNonHistoricalScalarVariables),
    mVisualizationNonHistoricalVectorVariables(rVisualizationNonHistoricalVectorVariables)
{
}

EmbeddedSkinVisualizationProcess::EmbeddedSkinVisualizationProcess(
    ModelPart& rModelPart,
    ModelPart& rVisualizationModelPart,
    Parameters rParameters)
    : Process()
    , mrModelPart(rModelPart)
    , mrVisualizationModelPart(rVisualizationModelPart)
    , mLevelSetType(
        [&] (Parameters& x) {
            x.ValidateAndAssignDefaults(StaticGetDefaultParameters());
            LevelSetType aux_level_set_type;
            CheckAndSetLevelSetType(x, aux_level_set_type);
            return aux_level_set_type;
        } (rParameters)
    )
    , mShapeFunctionsType(
        [&] (Parameters& x) {
            x.ValidateAndAssignDefaults(StaticGetDefaultParameters());
            ShapeFunctionsType aux_shape_func_type;
            CheckAndSetShapeFunctionsType(x, aux_shape_func_type);
            return aux_shape_func_type;
        } (rParameters)
    )
    , mReformModelPartAtEachTimeStep(
        [&] (Parameters& x) {
            x.ValidateAndAssignDefaults(StaticGetDefaultParameters());
            return x["reform_model_part_at_each_time_step"].GetBool();
        } (rParameters)
    )
    , mpNodalDistanceVariable(
        [&] (Parameters& x) -> const Variable<double>* {
            const Variable<double>* p_aux;
            switch (mLevelSetType) {
                case LevelSetType::Continuous: {
                    x.ValidateAndAssignDefaults(StaticGetDefaultParameters());
                    const std::string dist_var_name(CheckAndReturnDistanceVariableName(x, mLevelSetType));
                    p_aux = &(KratosComponents<Variable<double>>::Get(dist_var_name));
                    return p_aux;
                } default: {
                    p_aux = nullptr;
                    return p_aux;
                }
            }
        } (rParameters)
    )
    , mpElementalDistanceVariable(
        [&] (Parameters& x) -> const Variable<Vector>* {
            const Variable<Vector>* p_aux;
            switch (mLevelSetType) {
                case LevelSetType::Discontinuous: {
                    x.ValidateAndAssignDefaults(StaticGetDefaultParameters());
                    const std::string dist_var_name(CheckAndReturnDistanceVariableName(x, mLevelSetType));
                    p_aux = &(KratosComponents<Variable<Vector>>::Get(dist_var_name));
                    return p_aux;
                } default: {
                    p_aux = nullptr;
                    return p_aux;
                }
            }
        } (rParameters)
    )
    , mVisualizationScalarVariables(
        [&] (Parameters& x) -> std::vector<const Variable<double>*> {
            x.ValidateAndAssignDefaults(StaticGetDefaultParameters());
            std::vector<const Variable<double>*> aux_list;
            FillVariablesList<double>(x["visualization_variables"], aux_list);
            return aux_list;
        } (rParameters)
    )
    , mVisualizationVectorVariables(
        [&] (Parameters& x) -> std::vector<const Variable<array_1d<double,3>>*> {
            x.ValidateAndAssignDefaults(StaticGetDefaultParameters());
            std::vector<const Variable<array_1d<double,3>>*> aux_list;
            FillVariablesList<array_1d<double,3>>(x["visualization_variables"], aux_list);
            return aux_list;
        } (rParameters)
    )
    , mVisualizationNonHistoricalScalarVariables(
        [&] (Parameters& x) -> std::vector<const Variable<double>*> {
            x.ValidateAndAssignDefaults(StaticGetDefaultParameters());
            std::vector<const Variable<double>*> aux_list;
            FillVariablesList<double>(x["visualization_nonhistorical_variables"], aux_list);
            return aux_list;
        } (rParameters)
    )
    , mVisualizationNonHistoricalVectorVariables(
        [&] (Parameters& x) -> std::vector<const Variable<array_1d<double,3>>*> {
            x.ValidateAndAssignDefaults(StaticGetDefaultParameters());
            std::vector<const Variable<array_1d<double,3>>*> aux_list;
            FillVariablesList<array_1d<double,3>>(x["visualization_nonhistorical_variables"], aux_list);
            return aux_list;
        } (rParameters)
    )
{
}

EmbeddedSkinVisualizationProcess::EmbeddedSkinVisualizationProcess(
    Model& rModel,
    Parameters rParameters)
    : EmbeddedSkinVisualizationProcess(
        [&] (Model& x, Parameters& y) -> ModelPart& {
            y.ValidateAndAssignDefaults(StaticGetDefaultParameters());
            KRATOS_ERROR_IF(y["model_part_name"].GetString() == "") << "\'model_part_name\' is empty. Please provide the origin model part name." << std::endl;
            return x.GetModelPart(y["model_part_name"].GetString());
        } (rModel, rParameters),
        [&] (Model& x, Parameters& y) -> ModelPart& {
            y.ValidateAndAssignDefaults(StaticGetDefaultParameters());
            return CreateAndPrepareVisualizationModelPart(x, y);
        } (rModel, rParameters),
        rParameters)
{
}

void EmbeddedSkinVisualizationProcess::ExecuteBeforeSolutionLoop()
{
    // If mesh update is not required (constant level set function) create the visualization mesh once
    if (!mReformModelPartAtEachTimeStep) {
        this->CreateVisualizationMesh();
    }
}

void EmbeddedSkinVisualizationProcess::ExecuteBeforeOutputStep()
{
    // If required, set the visualization mesh
    if (mReformModelPartAtEachTimeStep) {
        this->CreateVisualizationMesh();
    }

    // Copy the origin model part nodal values for the non-intersection nodes
    this->CopyOriginNodalValues();

    // Compute new intersection nodes interpolation
    this->ComputeNewNodesInterpolation();
}

void EmbeddedSkinVisualizationProcess::ExecuteAfterOutputStep()
{
    if (mReformModelPartAtEachTimeStep){
        // Clear the new nodes hash map
        mCutNodesMap.clear();

        // Clear the new elements pointers vector
        mNewElementsPointers.clear();

        // Mark all the nodes, elements and conditions to be removed
        VariableUtils().SetFlag<ModelPart::NodesContainerType>(TO_ERASE, true, mrVisualizationModelPart.Nodes());
        VariableUtils().SetFlag<ModelPart::ElementsContainerType>(TO_ERASE, true, mrVisualizationModelPart.Elements());
        VariableUtils().SetFlag<ModelPart::ConditionsContainerType>(TO_ERASE, true, mrVisualizationModelPart.Conditions());

        // Remove all the nodes, elements and conditions
        mrVisualizationModelPart.RemoveNodes(TO_ERASE);
        mrVisualizationModelPart.RemoveElements(TO_ERASE);
        mrVisualizationModelPart.RemoveConditions(TO_ERASE);

        // Remove the positive and negative sides properties
        RemoveVisualizationProperties();
    }
}

int EmbeddedSkinVisualizationProcess::Check()
{
    // Check that model part is not empty
    KRATOS_ERROR_IF(mrModelPart.NumberOfNodes() == 0) << "There are no nodes in the origin model part." << std::endl;
    KRATOS_ERROR_IF(mrModelPart.NumberOfElements() == 0) << "There are no elements in the origin model part." << std::endl;

    // Required variables check
    const auto &r_orig_node = *mrModelPart.NodesBegin();

    // Check visualization scalar variables
    for (unsigned int i_var = 0; i_var < mVisualizationScalarVariables.size(); ++i_var){
        const Variable<double>& r_var = *(mVisualizationScalarVariables[i_var]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_var, r_orig_node);
    }

    // Check visualization vector variables
    for (unsigned int i_var = 0; i_var < mVisualizationVectorVariables.size(); ++i_var){
        const Variable<array_1d<double,3>>& r_var = *(mVisualizationVectorVariables[i_var]);
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(r_var, r_orig_node);
    }

    return 0;
}

/* Protected functions ****************************************************/

/* Private functions ******************************************************/

template<>
void EmbeddedSkinVisualizationProcess::SetPartitionIndexFromOriginNode<true>(
    const Node<3>& rOriginNode,
    Node<3>& rVisualizationNode)
{
    rVisualizationNode.FastGetSolutionStepValue(PARTITION_INDEX) = rOriginNode.FastGetSolutionStepValue(PARTITION_INDEX);
}

template<>
void EmbeddedSkinVisualizationProcess::SetPartitionIndexFromOriginNode<false>(
    const Node<3>& rOriginNode,
    Node<3>& rVisualizationNode)
{
}

template<>
void EmbeddedSkinVisualizationProcess::SetPartitionIndex<true>(
    const int PartitionIndex,
    Node<3>& rVisualizationNode)
{
    rVisualizationNode.FastGetSolutionStepValue(PARTITION_INDEX) = PartitionIndex;
}

template<>
void EmbeddedSkinVisualizationProcess::SetPartitionIndex<false>(
    const int PartitionIndex,
    Node<3>& rVisualizationNode)
{
}

void EmbeddedSkinVisualizationProcess::ComputeNewNodesInterpolation()
{
    // For all the new elements, compute the interpolation with the proper shape functions
    // Note that this can be done in parallel since the intersection nodes are duplicated
    const int n_new_elems = mNewElementsPointers.size();

    #pragma omp parallel for
    for (int i_elem = 0; i_elem < n_new_elems; ++i_elem){
        // Get element geometry
        auto it_elem = mNewElementsPointers.begin() + i_elem;
        const auto &r_geometry = it_elem->GetGeometry();
        const unsigned int n_points = r_geometry.PointsNumber();

        // For the generated elements, set the new interpolation values
        // Bear in mind that the non-split geometries kept the parent geometry
        // values, since they have been built using the origin model part nodes.
        for (unsigned int i_node = 0; i_node < n_points; ++i_node){
            // Search for the current node in the cut nodes hash map
            auto p_node = r_geometry(i_node);
            const CutNodesMapType::iterator cut_node_info = mCutNodesMap.find(p_node);

            // If it is not find, keep the origin value
            // Otherwise, the nodal values are computed with the modified shape functions values
            // Note that there is no necessity of locking the nodes since are always repeated
            if (cut_node_info != mCutNodesMap.end()){

                // Get the cut node info from the map tuple iterator
                Node<3>::Pointer p_edge_node_i, p_edge_node_j;
                double weight_edge_node_i, weight_edge_node_j;
                std::tie(p_edge_node_i, p_edge_node_j, weight_edge_node_i, weight_edge_node_j) = std::get<1>(*cut_node_info);

                // Interpolate the visualization variables
                InterpolateVariablesListValues<double, true>(p_node, p_edge_node_i, p_edge_node_j, weight_edge_node_i, weight_edge_node_j, mVisualizationScalarVariables);
                InterpolateVariablesListValues<double, false>(p_node, p_edge_node_i, p_edge_node_j, weight_edge_node_i, weight_edge_node_j, mVisualizationNonHistoricalScalarVariables);
                InterpolateVariablesListValues<array_1d<double,3>, true>(p_node, p_edge_node_i, p_edge_node_j, weight_edge_node_i, weight_edge_node_j, mVisualizationVectorVariables);
                InterpolateVariablesListValues<array_1d<double,3>, false>(p_node, p_edge_node_i, p_edge_node_j, weight_edge_node_i, weight_edge_node_j, mVisualizationNonHistoricalVectorVariables);
            }
        }
    }
}

template<class TDataType, bool IsHistorical>
void EmbeddedSkinVisualizationProcess::CopyVariablesListValues(
    const ModelPart::NodeIterator& rItOriginNode,
    ModelPart::NodeIterator& rItVisualizationNode,
    const std::vector<const Variable<TDataType>*>& rVariablesList)
{
    for (unsigned int i_var = 0; i_var < rVariablesList.size(); ++i_var){
        const Variable<TDataType> &r_var = *(rVariablesList[i_var]);
        const TDataType &r_origin_value = AuxiliaryGetValue<IsHistorical>(*rItOriginNode, r_var);
        TDataType &r_visualization_value = AuxiliaryGetValue<IsHistorical>(*rItVisualizationNode, r_var);
        r_visualization_value = r_origin_value;
    }
}

template<class TDataType, bool IsHistorical>
void EmbeddedSkinVisualizationProcess::InterpolateVariablesListValues(
    const Node<3>::Pointer& rpNode,
    const Node<3>::Pointer& rpNodeI,
    const Node<3>::Pointer& rpNodeJ,
    const double WeightI,
    const double WeightJ,
    const std::vector<const Variable<TDataType>*>& rVariablesList)
{
    for (unsigned int i_var = 0; i_var < rVariablesList.size(); ++i_var){
        const auto &r_var = *(rVariablesList[i_var]);
        const TDataType &r_edge_node_i_value = AuxiliaryGetValue<IsHistorical>(*rpNodeI, r_var);
        const TDataType &r_edge_node_j_value = AuxiliaryGetValue<IsHistorical>(*rpNodeJ, r_var);
        TDataType& r_value = AuxiliaryGetValue<IsHistorical>(*rpNode, r_var);
        r_value = WeightI * r_edge_node_i_value + WeightJ * r_edge_node_j_value;
    }
}

template<>
double& EmbeddedSkinVisualizationProcess::AuxiliaryGetValue<true>(
    Node<3>& rNode,
    const Variable<double>& rVariable)
{
    return rNode.FastGetSolutionStepValue(rVariable);
}

template<>
double& EmbeddedSkinVisualizationProcess::AuxiliaryGetValue<false>(
    Node<3>& rNode,
    const Variable<double>& rVariable)
{
    return rNode.GetValue(rVariable);
}

template<>
array_1d<double,3>& EmbeddedSkinVisualizationProcess::AuxiliaryGetValue<true>(
    Node<3>& rNode,
    const Variable<array_1d<double,3>>& rVariable)
{
    return rNode.FastGetSolutionStepValue(rVariable);
}

template<>
array_1d<double,3>& EmbeddedSkinVisualizationProcess::AuxiliaryGetValue<false>(
    Node<3>& rNode,
    const Variable<array_1d<double,3>>& rVariable)
{
    return rNode.GetValue(rVariable);
}

void EmbeddedSkinVisualizationProcess::CreateVisualizationMesh()
{
    // Copy the original nodes to the visualization model part
    if (mrVisualizationModelPart.IsDistributed()) {
        this->CopyOriginNodes<true>();
    } else {
        this->CopyOriginNodes<false>();
    }

    // Creates the visualization model part geometrical entities (elements and conditions)
    this->CreateVisualizationGeometries();

    // If MPI, this creates the communication plan among processes
    ParallelEnvironment::CreateFillCommunicator(mrVisualizationModelPart)->Execute();

    // Initialize (allocate) non-historical variables
    InitializeNonHistoricalVariables<double>(mVisualizationNonHistoricalScalarVariables);
    InitializeNonHistoricalVariables<array_1d<double,3>>(mVisualizationNonHistoricalVectorVariables);
}

template<const bool IsDistributed>
void EmbeddedSkinVisualizationProcess::CopyOriginNodes()
{
    // Creates a copy of all the origin model part nodes to the visualization model part
    // Note that the original nodes will be reused when creating the splitting geometries
    const int n_nodes = mrModelPart.NumberOfNodes();
    ModelPart::NodeIterator orig_nodes_begin = mrModelPart.NodesBegin();
    for (int i_node = 0; i_node < n_nodes; ++i_node){
        auto it_node = orig_nodes_begin + i_node;
        auto p_vis_node = mrVisualizationModelPart.CreateNewNode(it_node->Id(), *it_node);
        SetPartitionIndexFromOriginNode<IsDistributed>(*it_node, *p_vis_node);
    }
}

void EmbeddedSkinVisualizationProcess::CopyOriginNodalValues()
{
    const unsigned int n_old_nodes = mrModelPart.NumberOfNodes();

    #pragma omp parallel for
    for (int i_node = 0; i_node < static_cast<int>(n_old_nodes); ++i_node){
        const auto it_origin_node = mrModelPart.NodesBegin() + i_node;
        auto it_visualization_node = mrVisualizationModelPart.NodesBegin() + i_node;

        CopyVariablesListValues<double, true>(it_origin_node, it_visualization_node, mVisualizationScalarVariables);
        CopyVariablesListValues<double, false>(it_origin_node, it_visualization_node, mVisualizationNonHistoricalScalarVariables);
        CopyVariablesListValues<array_1d<double,3>, true>(it_origin_node, it_visualization_node, mVisualizationVectorVariables);
        CopyVariablesListValues<array_1d<double,3>, false>(it_origin_node, it_visualization_node, mVisualizationNonHistoricalVectorVariables);
    }
}

void EmbeddedSkinVisualizationProcess::CreateVisualizationGeometries()
{
    int n_nodes = mrModelPart.NumberOfNodes();
    int n_elems = mrModelPart.NumberOfElements();
    int n_conds = mrModelPart.NumberOfConditions();

    // Get the current rank
    // New intersection nodes will be assigned this one
    // Note that these are duplicated among intersected elements so there is no sync required
    const auto& r_comm = mrModelPart.GetCommunicator();
    const int my_pyd = r_comm.MyPID();
    std::function<void(const int, Node<3>&)> set_partition_index_func;
    if (r_comm.IsDistributed()) {
        set_partition_index_func = &(this->SetPartitionIndex<true>);
    } else {
        set_partition_index_func = &(this->SetPartitionIndex<false>);
    }

    // Initialize the ids. for the new entries
    // For the temporal ids. there is no necessity of synchronizing between processors
    unsigned int temp_node_id = (n_nodes > 0) ? ((mrModelPart.NodesEnd()-1)->Id()) + 1 : 1;
    unsigned int temp_elem_id = (n_elems > 0) ? ((mrModelPart.ElementsEnd()-1)->Id()) + 1 : 1;
    unsigned int temp_cond_id = (n_conds > 0) ? ((mrModelPart.ConditionsEnd()-1)->Id()) + 1 : 1;

    // Auxiliar vectors to store pointers to the new entities
    // These vectors will be use when renumbering the entities ids.
    ModelPart::NodesContainerType new_nodes_vect;
    ModelPart::ConditionsContainerType new_conds_vect;

    // Set the new entities properties
    Properties::Pointer p_pos_prop, p_neg_prop;
    std::tie(p_pos_prop, p_neg_prop) = this->SetVisualizationProperties();

    // Add the elements to the visualization model part
    for (int i_elem = 0; i_elem < n_elems; ++i_elem){
        ModelPart::ElementIterator it_elem = mrModelPart.ElementsBegin() + i_elem;

        // Get element geometry
        const Geometry<Node<3>>::Pointer p_geometry = it_elem->pGetGeometry();
        const unsigned int n_nodes = p_geometry->PointsNumber();
        const Vector nodal_distances = this->SetDistancesVector(it_elem);

        // Check if the element is split
        const double zero_tol = 1.0e-10;
        const double nodal_distances_norm = norm_2(nodal_distances);
        if (nodal_distances_norm < zero_tol) {
            KRATOS_WARNING_FIRST_N("EmbeddedSkinVisualizationProcess", 10) << "Element: " << it_elem->Id() << " Distance vector norm: " << nodal_distances_norm << ". Please check the level set function." << std::endl;
        }
        const bool is_split = this->ElementIsSplit(p_geometry, nodal_distances);

        // If the element is split, create the new entities
        if (is_split){
            // Set the split utility and compute the splitting pattern
            ModifiedShapeFunctions::Pointer p_modified_shape_functions = this->SetModifiedShapeFunctionsUtility(p_geometry, nodal_distances);
            DivideGeometry::Pointer p_split_utility = p_modified_shape_functions->pGetSplittingUtil();

            // Create the auxiliar map that will be used to generate the skin
            std::unordered_map<std::pair<unsigned int,bool>, unsigned int, Hash, KeyEqual> new_nodes_map;

            // Get the split geometries from the splitting pattern
            const unsigned int n_pos_split_geom = (p_split_utility->mPositiveSubdivisions).size();
            const unsigned int n_neg_split_geom = (p_split_utility->mNegativeSubdivisions).size();
            std::vector<DivideGeometry::IndexedPointGeometryPointerType> split_geometries;
            split_geometries.reserve(n_pos_split_geom + n_neg_split_geom);
            split_geometries.insert(split_geometries.end(), (p_split_utility->mPositiveSubdivisions).begin(), (p_split_utility->mPositiveSubdivisions).end());
            split_geometries.insert(split_geometries.end(), (p_split_utility->mNegativeSubdivisions).begin(), (p_split_utility->mNegativeSubdivisions).end());

            // Create the split geometries in the visualization model part
            for (unsigned int i_geom = 0; i_geom < split_geometries.size(); ++i_geom){
                const bool pos_side = i_geom < n_pos_split_geom ? true : false;
                const DivideGeometry::IndexedPointGeometryPointerType p_sub_geom = split_geometries[i_geom];
                const unsigned int sub_geom_n_nodes = p_sub_geom->PointsNumber();

                // Fill the new element nodes array
                Element::NodesArrayType sub_geom_nodes_array;
                for (unsigned int i_sub_geom_node = 0; i_sub_geom_node < sub_geom_n_nodes; ++i_sub_geom_node){

                    DivideGeometry::IndexedPointType &sub_geom_node = p_sub_geom->operator[](i_sub_geom_node);
                    const unsigned int local_id = sub_geom_node.Id();

                    // Non-intersection node. Get the copied original geometry nodes
                    // Note that we do not reuse the original model part node since this forces
                    // the new intersection nodes to have the same buffer size, variables list, etc.
                    if (local_id < sub_geom_n_nodes){
                        const unsigned int aux_gl_id = (p_geometry->operator()(local_id))->Id();
                        sub_geom_nodes_array.push_back(mrVisualizationModelPart.pGetNode(aux_gl_id));
                    // Intersection node. The node is located in an intersected edge.
                    // Thus, take it from the geometry created by the splitting util.
                    } else {
                        const unsigned int intersected_edge_id = local_id - n_nodes;

                        // Get the intersected edge node_i and node_j
                        const unsigned int node_i = (p_split_utility->GetEdgeIdsI())[intersected_edge_id];
                        const unsigned int node_j = (p_split_utility->GetEdgeIdsJ())[intersected_edge_id];

                        // Create a new node based in the indexed point auxiliar node
                        const array_1d<double, 3> point_coords = sub_geom_node.Coordinates();
                        Node<3>::Pointer p_new_node = mrVisualizationModelPart.CreateNewNode(temp_node_id, point_coords[0], point_coords[1], point_coords[2]);
                        sub_geom_nodes_array.push_back(p_new_node);
                        new_nodes_vect.push_back(p_new_node);
                        set_partition_index_func(my_pyd, *p_new_node);

                        // Add the new node info to the hash map
                        const Node<3>::Pointer p_node_i = p_geometry->operator()(node_i);
                        const Node<3>::Pointer p_node_j = p_geometry->operator()(node_j);

                        Matrix edge_N_values;
                        if (i_geom < n_pos_split_geom){
                            p_modified_shape_functions->ComputeShapeFunctionsOnPositiveEdgeIntersections(edge_N_values);
                        } else {
                            p_modified_shape_functions->ComputeShapeFunctionsOnNegativeEdgeIntersections(edge_N_values);
                        }

                        const double node_i_weight = edge_N_values(intersected_edge_id, node_i);
                        const double node_j_weight = edge_N_values(intersected_edge_id, node_j);

                        auto new_node_info = std::make_tuple(p_node_i, p_node_j, node_i_weight, node_j_weight);
                        mCutNodesMap.insert(CutNodesMapType::value_type(p_new_node, new_node_info));

                        // Link the new node global id. to its local index in the splitting util. Note that the side
                        // (positive or negative) is saved as well. This will be used when setting up the skin conditions.
                        std::pair<unsigned int, bool> aux_info(local_id,pos_side);
                        std::pair<std::pair<unsigned int,bool>, unsigned int> new_pair(aux_info, temp_node_id);
                        new_nodes_map.insert(new_pair);

                        // Update the new nodes id. counter
                        temp_node_id++;
                    }
                }

                // Set the new element properties
                Properties::Pointer p_elem_prop = pos_side ? p_pos_prop : p_neg_prop;

                // Create the new element
                Element::Pointer p_new_elem = it_elem->Create(temp_elem_id, sub_geom_nodes_array, p_elem_prop);
                mNewElementsPointers.push_back(p_new_elem);

                // Update the new elements id. counter
                temp_elem_id++;
            }

            // Get the interface geometries from the splitting pattern
            const unsigned int n_pos_interface_geom = (p_split_utility->mPositiveInterfaces).size();
            const unsigned int n_neg_interface_geom = (p_split_utility->mNegativeInterfaces).size();

            std::vector<DivideGeometry::IndexedPointGeometryPointerType> split_interface_geometries;
            split_interface_geometries.reserve(n_pos_interface_geom + n_neg_interface_geom);
            split_interface_geometries.insert(split_interface_geometries.end(), (p_split_utility->mPositiveInterfaces).begin(), (p_split_utility->mPositiveInterfaces).end());
            split_interface_geometries.insert(split_interface_geometries.end(), (p_split_utility->mNegativeInterfaces).begin(), (p_split_utility->mNegativeInterfaces).end());

            // Create the split interface geometries in the visualization model part
            for (unsigned int i_int_geom = 0; i_int_geom < split_interface_geometries.size(); ++i_int_geom){
                const bool int_pos_side = (i_int_geom < n_pos_interface_geom) ? true : false;
                DivideGeometry::IndexedPointGeometryPointerType p_int_sub_geom = split_interface_geometries[i_int_geom];
                GeometryData::KratosGeometryType p_int_sub_geom_type = p_int_sub_geom->GetGeometryType();
                const unsigned int sub_int_geom_n_nodes = p_int_sub_geom->PointsNumber();

                // Fill the new condition nodes array
                Condition::NodesArrayType sub_int_geom_nodes_array;
                for (unsigned int i_node = 0; i_node < sub_int_geom_n_nodes; ++i_node){

                    DivideGeometry::IndexedPointType &sub_int_geom_node = p_int_sub_geom->operator[](i_node);
                    const unsigned int local_id = sub_int_geom_node.Id();

                    // Get the global id from the intersection nodes map
                    unsigned int global_id;
                    std::pair<unsigned int,bool> aux_int_info(local_id,int_pos_side);
                    auto got = new_nodes_map.find(aux_int_info);
                    if (got != new_nodes_map.end()){
                        global_id = got->second;
                    } else {
                        const std::string side = int_pos_side ? "positive" : "negative";
                        KRATOS_ERROR << "Local id " << std::get<0>(aux_int_info) << " in " << side << " side not found in new nodes map for element " << it_elem->Id();
                    }

                    sub_int_geom_nodes_array.push_back(mrVisualizationModelPart.pGetNode(global_id));
                }

                // Set the new condition geometry
                Geometry< Node<3> >::Pointer p_new_geom = SetNewConditionGeometry(
                    p_int_sub_geom_type,
                    sub_int_geom_nodes_array);

                // Set the new condition properties
                Properties::Pointer p_cond_prop = (i_int_geom < n_pos_interface_geom)? p_pos_prop : p_neg_prop;

                // Create the new condition
                Condition::Pointer p_new_cond = Kratos::make_intrusive<Condition>(temp_cond_id, p_new_geom, p_cond_prop);
                new_conds_vect.push_back(p_new_cond);
                mrVisualizationModelPart.AddCondition(p_new_cond);

                // Update the new elements id. counter
                temp_cond_id++;

            }
        // Otherwise add an element with the original geometry
        } else {
            // Copy the current element as it was in the origin model part according to its distance sign
            const bool is_positive = this->ElementIsPositive(p_geometry, nodal_distances);
            if (is_positive){
                mrVisualizationModelPart.AddElement(it_elem->Create(it_elem->Id(), p_geometry, p_pos_prop));
            } else {
                mrVisualizationModelPart.AddElement(it_elem->Create(it_elem->Id(), p_geometry, p_neg_prop));
            }
        }
    }

    // Once all the entities have been created, renumber the ids.
    // Created entities local number partial reduction
    const auto& r_data_comm = mrModelPart.GetCommunicator().GetDataCommunicator();
    int n_nodes_local = new_nodes_vect.size();
    int n_elems_local = mNewElementsPointers.size();
    int n_conds_local = new_conds_vect.size();

    std::vector<int> local_data{n_nodes_local, n_elems_local, n_conds_local};
    std::vector<int> reduced_data{0, 0, 0};
    r_data_comm.ScanSum(local_data, reduced_data);

    int n_nodes_local_scansum = reduced_data[0];
    int n_elems_local_scansum = reduced_data[1];
    int n_conds_local_scansum = reduced_data[2];

    // Origin model part number of entities
    int n_nodes_orig = mrModelPart.NumberOfNodes();
    int n_elems_orig = mrModelPart.NumberOfElements();
    int n_conds_orig = mrModelPart.NumberOfConditions();
    local_data = {n_nodes_orig, n_elems_orig, n_conds_orig};
    r_data_comm.SumAll(local_data, reduced_data);
    n_nodes_orig = reduced_data[0];
    n_elems_orig = reduced_data[1];
    n_conds_orig = reduced_data[2];

    // Initialize the new ids. values
    std::size_t new_node_id(n_nodes_orig + n_nodes_local_scansum - n_nodes_local + 1);
    std::size_t new_elem_id(n_elems_orig + n_elems_local_scansum - n_elems_local + 1);
    std::size_t new_cond_id(n_conds_orig + n_conds_local_scansum - n_conds_local + 1);

    // Set the new entities ids. values
    auto new_nodes_begin = new_nodes_vect.begin();
    auto new_conds_begin = new_conds_vect.begin();
    auto new_elems_begin = mNewElementsPointers.begin();

    for (int i_node = 0; i_node < static_cast<int>(new_nodes_vect.size()); ++i_node){
        auto it_node = new_nodes_begin + i_node;
        const unsigned int new_id = new_node_id + i_node;
        it_node->SetId(new_id);
    }

    for (int i_cond = 0; i_cond < static_cast<int>(new_conds_vect.size()); ++i_cond){
        auto it_cond = new_conds_begin + i_cond;
        const unsigned int new_id = new_cond_id + i_cond;
        it_cond->SetId(new_id);
    }

    for (int i_elem = 0; i_elem < static_cast<int>(mNewElementsPointers.size()); ++i_elem){
        auto it_elem = new_elems_begin + i_elem;
        const unsigned int new_id = new_elem_id + i_elem;
        it_elem->SetId(new_id);
    }

    mrVisualizationModelPart.AddElements(mNewElementsPointers.begin(), mNewElementsPointers.end());
    mrVisualizationModelPart.AddConditions(new_conds_vect.begin(), new_conds_vect.end());

    // Wait for all nodes to renumber its nodes
    r_data_comm.Barrier();
}

template<class TDataType>
void EmbeddedSkinVisualizationProcess::InitializeNonHistoricalVariables(const std::vector<const Variable<TDataType>*>& rNonHistoricalVariablesVector)
{
    const int n_nodes = mrVisualizationModelPart.NumberOfNodes();
#pragma omp parallel for
    for (int i_node = 0; i_node < n_nodes; ++i_node) {
        auto it_node = mrVisualizationModelPart.NodesBegin() + i_node;
        for (auto& r_var : rNonHistoricalVariablesVector) {
            it_node->SetValue(*r_var, r_var->Zero());
        }
    }
}

bool EmbeddedSkinVisualizationProcess::ElementIsPositive(
    Geometry<Node<3>>::Pointer pGeometry,
    const Vector &rNodalDistances)
{
    const unsigned int pts_number = pGeometry->PointsNumber();
    unsigned int n_pos (0);

    for (unsigned int i_node = 0; i_node < pts_number; ++i_node){
        if (rNodalDistances[i_node] > 0.0)
            n_pos++;
    }

    const bool is_positive = (n_pos == pts_number) ? true : false;

    return is_positive;
}

bool EmbeddedSkinVisualizationProcess::ElementIsSplit(
    Geometry<Node<3>>::Pointer pGeometry,
    const Vector &rNodalDistances)
{
    const unsigned int pts_number = pGeometry->PointsNumber();
    unsigned int n_pos (0), n_neg(0);

    for (unsigned int i_node = 0; i_node < pts_number; ++i_node){
        if (rNodalDistances[i_node] > 0.0)
            n_pos++;
        else
            n_neg++;
    }

    const bool is_split = (n_pos > 0 && n_neg > 0) ? true : false;

    return is_split;
}

const Vector EmbeddedSkinVisualizationProcess::SetDistancesVector(ModelPart::ElementIterator ItElem)
{
    const auto &r_geom = ItElem->GetGeometry();
    Vector nodal_distances(r_geom.PointsNumber());

    switch (mLevelSetType) {
        // Continuous nodal distance function case
        case LevelSetType::Continuous:
            for (unsigned int i_node = 0; i_node < r_geom.PointsNumber(); ++i_node) {
                nodal_distances[i_node] = r_geom[i_node].FastGetSolutionStepValue(*mpNodalDistanceVariable);
            }
            break;
        // Discontinuous elemental distance function case
        case LevelSetType::Discontinuous:
            nodal_distances = ItElem->GetValue(*mpElementalDistanceVariable);
            break;
        // Default error case
        default:
            KRATOS_ERROR << "Asking for a non-implemented modified shape functions type.";
    }

    return nodal_distances;
}

ModifiedShapeFunctions::Pointer EmbeddedSkinVisualizationProcess::SetModifiedShapeFunctionsUtility(
    const Geometry<Node<3>>::Pointer pGeometry,
    const Vector& rNodalDistances)
{
    // Get the geometry type
    const GeometryData::KratosGeometryType geometry_type = pGeometry->GetGeometryType();

    // Return the modified shape functions utility
    switch (mShapeFunctionsType) {
        case ShapeFunctionsType::Standard:
            switch (geometry_type) {
                case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
                    return Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(pGeometry, rNodalDistances);
                case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
                    return Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(pGeometry, rNodalDistances);
                default:
                    KRATOS_ERROR << "Asking for a non-implemented modified shape functions geometry.";
            }
        case ShapeFunctionsType::Ausas:
            switch (geometry_type) {
                case GeometryData::KratosGeometryType::Kratos_Triangle2D3:
                    return Kratos::make_shared<Triangle2D3AusasModifiedShapeFunctions>(pGeometry, rNodalDistances);
                case GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4:
                    return Kratos::make_shared<Tetrahedra3D4AusasModifiedShapeFunctions>(pGeometry, rNodalDistances);
                default:
                    KRATOS_ERROR << "Asking for a non-implemented Ausas modified shape functions geometry.";
            }
        default:
            KRATOS_ERROR << "Asking for a non-implemented modified shape functions type.";
    }
}

Geometry< Node<3> >::Pointer EmbeddedSkinVisualizationProcess::SetNewConditionGeometry(
    const GeometryData::KratosGeometryType &rOriginGeometryType,
    const Condition::NodesArrayType &rNewNodesArray)
{
    switch(rOriginGeometryType){
        case GeometryData::KratosGeometryType::Kratos_Line2D2:
            return Kratos::make_shared<Line2D2< Node<3> > >(rNewNodesArray);
        case GeometryData::KratosGeometryType::Kratos_Triangle3D3:
            return Kratos::make_shared<Triangle3D3< Node<3> > >(rNewNodesArray);
        default:
            KRATOS_ERROR << "Implement the visualization for the intersection geometry type " << rOriginGeometryType;
    }
}

std::tuple< Properties::Pointer , Properties::Pointer > EmbeddedSkinVisualizationProcess::SetVisualizationProperties()
{
    unsigned int max_prop_id = 0;
    for (auto it_prop = mrModelPart.GetRootModelPart().PropertiesBegin(); it_prop < mrModelPart.GetRootModelPart().PropertiesEnd(); ++it_prop){
        if (max_prop_id < it_prop->Id()){
            max_prop_id = it_prop->Id();
        }
    }
    Properties::Pointer p_pos_prop = Kratos::make_shared<Properties>(max_prop_id + 1);
    Properties::Pointer p_neg_prop = Kratos::make_shared<Properties>(max_prop_id + 2);
    mrVisualizationModelPart.AddProperties(p_pos_prop);
    mrVisualizationModelPart.AddProperties(p_neg_prop);

    return std::make_tuple(p_pos_prop , p_neg_prop);
}

void EmbeddedSkinVisualizationProcess::RemoveVisualizationProperties()
{
    // Search for the maximum property id in the base model part
    unsigned int max_prop_id = 0;
    for (auto it_prop = mrModelPart.GetRootModelPart().PropertiesBegin(); it_prop < mrModelPart.GetRootModelPart().PropertiesEnd(); ++it_prop){
        if (max_prop_id < it_prop->Id()){
            max_prop_id = it_prop->Id();
        }
    }

    // Check that the positive and negative sides properties exist
    KRATOS_ERROR_IF_NOT(mrVisualizationModelPart.HasProperties(max_prop_id + 1)) << "Visualization model part has no property " << max_prop_id + 1 << std::endl;
    KRATOS_ERROR_IF_NOT(mrVisualizationModelPart.HasProperties(max_prop_id + 2)) << "Visualization model part has no property " << max_prop_id + 2 << std::endl;

    // Remove the positive and negative sides properties
    mrVisualizationModelPart.RemoveProperties(max_prop_id + 1);
    mrVisualizationModelPart.RemoveProperties(max_prop_id + 2);
}

};  // namespace Kratos.
