//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez and Marc Nunez
//

// External includes
#include <queue>

// Project includes
#include "includes/model_part.h"
#include "utilities/variable_utils.h"
#include "utilities/builtin_timer.h"
#include "input_output/vtk_output.h"
#include "input_output/stl_io.h"
#include "containers/model.h"

// Application includes
#include "define_3d_wake_operation.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_utilities/potential_flow_utilities.h"
#include "custom_processes/move_model_part_process.h"
#include "processes/calculate_distance_to_skin_process.h"

namespace Kratos 
{

// Constructor for Define3DWakeOperation Operation
Define3DWakeOperation::Define3DWakeOperation(
    Model& rModel,
    Parameters mParameters)
{
    KRATOS_TRY;
    mParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
    mWakeSTLFileName          = mParameters["wake_stl_file_name"].GetString();
    mWakeNormal               = mParameters["wake_normal"].GetVector();
    mWakeTraslationDirection  = mParameters["wake_translation_direction"].GetVector();
    mVisualizeWakeVTK         = mParameters["visualize_wake_vtk"].GetBool();
    mShedWakeFromTrailingEdge = mParameters["shed_wake_from_trailing_edge"].GetBool();
    mShedWakeLength           = mParameters["shed_wake_length"].GetDouble();
    mShedWakeElementSize      = mParameters["shed_wake_element_size"].GetDouble();
    mShedGrowFactor           = mParameters["shed_wake_grow_factor"].GetDouble();
    mShedProjectionRootEdge   = mParameters["shed_wake_projection_root_edge"].GetDouble();
    mWakeDistanceTolerance    = mParameters["wake_distance_tolerance"].GetDouble();
    mEchoLevel                = mParameters["echo_level"].GetInt();
    
    KRATOS_ERROR_IF(mWakeNormal.size() != 3)
        << "The mWakeNormal should be a vector with 3 components!"
        << std::endl;

    // Saving the modelparts
    mpTrailingEdgeModelPart = &(rModel.GetModelPart(mParameters["trailing_edge_model_part_name"].GetString()));
    mpBodyModelPart         = &(rModel.GetModelPart(mParameters["body_model_part_name"].GetString()));
    mpWakeModelPart         = &(rModel.CreateModelPart("wake_model_part"));
    mpRootModelPart         = &(mpBodyModelPart->GetRootModelPart());

    if (mParameters["upper_surface_model_part_name"].GetString() != "")
        mpUpperSurfaceModelPart = &(rModel.GetModelPart(mParameters["upper_surface_model_part_name"].GetString()));

    if (mParameters["lower_surface_model_part_name"].GetString() != "")
        mpLowerSurfaceModelPart = &(rModel.GetModelPart(mParameters["lower_surface_model_part_name"].GetString()));

    if (mParameters["root_points_model_part_name"].GetString() != "")
        mpRootPointsModelPart = &(rModel.GetModelPart(mParameters["root_points_model_part_name"].GetString()));

    if (mParameters["tip_points_model_part_name"].GetString() != "")
        mpTipPointsModelPart = &(rModel.GetModelPart(mParameters["tip_points_model_part_name"].GetString()));
    
    if (mParameters["blunt_te_surface_model_part_name"].GetString() !="")
        mpBluntTESurfaceModelPart = &(rModel.GetModelPart(mParameters["blunt_te_surface_model_part_name"].GetString()));
    KRATOS_CATCH("");
}

void Define3DWakeOperation::Execute()
{
    KRATOS_TRY;

    InitializeVariables();
    
    InitializeTrailingEdgeSubModelpart();
    
    InitializeWakeSubModelpart();
    
    MarkTrailingEdgeAndBluntNodes();
    
    ComputeWingLowerSurfaceNormals();
    
    AddTrailingEdgeConditionsAndFindRootAndTipNodes();

    ComputeAndSaveLocalWakeNormal();
    
    if(mShedWakeFromTrailingEdge || mWakeSTLFileName == ""){
        ShedWakeSurfaceFromTheTrailingEdge();
    }else
    {
        LoadSTL();
    }

    if (norm_2(mWakeTraslationDirection) > std::numeric_limits<double>::epsilon())
        MoveWakeModelPart();
    
    if (mVisualizeWakeVTK)
        VisualizeWake();
    
    MarkWakeElements();
    
    RecomputeNodalDistancesToWakeOrWingLowerSurface();
    
    MarkKuttaElements();
    
    SaveLocalWakeNormalInElements();
    
    AddWakeNodesToWakeModelPart();
    KRATOS_CATCH("");
}

// This function initializes the variables
void Define3DWakeOperation::InitializeVariables() 
{
    KRATOS_TRY;
    block_for_each(mpRootModelPart->Nodes(), [&](Node& r_nodes)
    {
        r_nodes.SetValue(UPPER_SURFACE, false);
        r_nodes.SetValue(LOWER_SURFACE, false);
        r_nodes.SetValue(TRAILING_EDGE, false);
        r_nodes.SetValue(WAKE_DISTANCE, 0.0);
    });
    auto& r_elements = mpRootModelPart->Elements();
    VariableUtils().SetNonHistoricalVariable(WAKE, 0, r_elements);
    
    // Save wake normal
    mpRootModelPart->GetProcessInfo()[WAKE_NORMAL] = mWakeNormal;
    // Save free stream velocity direction as wake direction
    mWakeDirection = mpRootModelPart->GetProcessInfo()[FREE_STREAM_VELOCITY_DIRECTION];
    // Computing the norm of the free_stream_velocity_direction vector
    const double norm = std::sqrt(inner_prod(mWakeDirection, mWakeDirection));
    const double eps = std::numeric_limits<double>::epsilon();
    KRATOS_ERROR_IF(norm < eps)
    << "The norm of the free stream velocity should be different than 0."
    << std::endl;
    // Compute span direction as the cross product: mWakeNormal x mWakeDirection
    MathUtils<double>::CrossProduct(mSpanDirection, mWakeNormal, mWakeDirection);
    KRATOS_CATCH("");
}

// This function initializes the variables and removes all of the elements and
// nodes of the trailing edge submodelpart
void Define3DWakeOperation::InitializeTrailingEdgeSubModelpart() const
{
    KRATOS_TRY;
    if(mpRootModelPart->HasSubModelPart("trailing_edge_elements_model_part"))
    {
        // Clearing the variables and elements of the already existing
        // trailing_edge_sub_model_part
        ModelPart& trailing_edge_sub_model_part =
            mpRootModelPart->GetSubModelPart("trailing_edge_elements_model_part");

        for (auto& r_element : trailing_edge_sub_model_part.Elements()){
            r_element.SetValue(TRAILING_EDGE, false);
            r_element.SetValue(KUTTA, false);
            r_element.Reset(STRUCTURE);
            r_element.Set(TO_ERASE, true);
        }
        auto& r_nodes = trailing_edge_sub_model_part.Nodes();
        VariableUtils().SetFlag(TO_ERASE, false, r_nodes);
        trailing_edge_sub_model_part.RemoveElements(TO_ERASE);
        trailing_edge_sub_model_part.RemoveNodes(TO_ERASE);
    }
    else{
        // Creating the trailing_edge_sub_model_part
        mpRootModelPart->CreateSubModelPart("trailing_edge_elements_model_part");
    }
    KRATOS_CATCH("");
}

// This function initializes the variables and removes all of the elements and
// nodes of the wake submodelpart
void Define3DWakeOperation::InitializeWakeSubModelpart() const
{
    KRATOS_TRY;
    if(mpRootModelPart->HasSubModelPart("wake_elements_model_part"))
    {
        // Clearing the variables and elements of the already existing
        // wake_sub_model_part
        ModelPart& wake_sub_model_part =
            mpRootModelPart->GetSubModelPart("wake_elements_model_part");

        for (auto& r_element : wake_sub_model_part.Elements()){
            r_element.SetValue(WAKE, false);
            r_element.SetValue(WAKE_ELEMENTAL_DISTANCES, ZeroVector(3));
            r_element.Set(TO_ERASE, true);
        }
        auto& r_nodes = wake_sub_model_part.Nodes();
        VariableUtils().SetFlag(TO_ERASE, false, r_nodes);
        wake_sub_model_part.RemoveElements(TO_ERASE);
        wake_sub_model_part.RemoveNodes(TO_ERASE);
    }
    else{
        // Creating the wake_sub_model_part
        mpRootModelPart->CreateSubModelPart("wake_elements_model_part");
    }
    KRATOS_CATCH("");
}

// This function marks the trailing edge and 
// blunt nodes where the airfoil does not end in a thin edge, but in a flat termination
void Define3DWakeOperation::MarkTrailingEdgeAndBluntNodes()
{
    KRATOS_TRY;
    for (auto& r_node : mpTrailingEdgeModelPart->Nodes()) {
        r_node.SetValue(TRAILING_EDGE, true);
    }

    if (mpBluntTESurfaceModelPart)
    {
        for (auto& r_node : mpBluntTESurfaceModelPart->Nodes()) {
            mBluntIds.insert(r_node.Id());
            r_node.SetValue(TRAILING_EDGE, true);
        }
    }
    KRATOS_CATCH("");
}

// This function computes the wing lower surface normals and marks the upper
// and lower surfaces. The wing lower surface normals are used later in
// RecomputeComputeNodalDistancesToWakeOrWingLowerSurface inside the
// MarkKuttaElements function to check whether nodes are above or below the wake
// Upper and Lower surfaces model parts are used, if provided.
void Define3DWakeOperation::ComputeWingLowerSurfaceNormals() const
{
    KRATOS_TRY;
    // TO DISCUSS: Is it worth it to run these loops in parallel?
    // So far it has not been noted much difference in terms of speed.
    // Also, when ran in parallel, the result can be random.
    // Mark upper surface
    if (mpUpperSurfaceModelPart)
    { 
        KRATOS_INFO_IF("Define3DWakeOperation", mEchoLevel > 0) << "...Using upper_surface_model_part" << std::endl;
        for (auto& r_cond : mpUpperSurfaceModelPart->Conditions()) {
            auto& r_geometry = r_cond.GetGeometry();
            for (unsigned int j = 0; j < r_geometry.size(); j++) {
                r_geometry[j].SetLock();
                r_geometry[j].SetValue(UPPER_SURFACE, true);
                r_geometry[j].UnSetLock();
            }
        }
    }else
    {
        KRATOS_WARNING("Define3DWakeOperation") << "Marking Upper Surface automatically. Set \'upper_surface_model_part_name\' instead." << std::endl;
        KRATOS_INFO_IF("Define3DWakeOperation", mEchoLevel > 0) << "...Marking Upper Surface" << std::endl;
        for (auto& r_cond : mpBodyModelPart->Conditions()) {
            auto& r_geometry = r_cond.GetGeometry();
            const auto& surface_normal = r_geometry.UnitNormal(0);
            const double projection = inner_prod(surface_normal, mWakeNormal);
    
            if(!(projection > 0.0)){
                for (unsigned int j = 0; j < r_geometry.size(); j++) {
                    r_geometry[j].SetLock();
                    r_geometry[j].SetValue(UPPER_SURFACE, true);
                    r_geometry[j].UnSetLock();
                }
            }
        }
    }
    
    // Mark lower surface
    if (mpLowerSurfaceModelPart)
    {
        KRATOS_INFO_IF("Define3DWakeOperation", mEchoLevel > 0) << "...Using lower_surface_model_part" << std::endl;
        for (auto& r_cond : mpLowerSurfaceModelPart->Conditions()) {
            auto& r_geometry = r_cond.GetGeometry();
            const auto& surface_normal = r_geometry.UnitNormal(0);
            const double projection = inner_prod(surface_normal, mWakeNormal);
            const double switching_factor = projection > 0.0 ? 1.0 : -1.0;
            
            for (unsigned int j = 0; j < r_geometry.size(); j++) {
                r_geometry[j].SetLock();
                r_geometry[j].SetValue(NORMAL, surface_normal*switching_factor);
                r_geometry[j].SetValue(LOWER_SURFACE, true);
                r_geometry[j].UnSetLock();
            }
        }
    }else
    {
        KRATOS_WARNING("Define3DWakeOperation") << "Marking Lower Surface automatically. Set \'lower_surface_model_part_name\' instead." << std::endl;
        KRATOS_INFO_IF("Define3DWakeOperation", mEchoLevel > 0) << "...Marking Lower Surface" << std::endl;
        for (auto& r_cond : mpBodyModelPart->Conditions()) {
            auto& r_geometry = r_cond.GetGeometry();
            const auto& surface_normal = r_geometry.UnitNormal(0);
            const double projection = inner_prod(surface_normal, mWakeNormal);
    
            if(projection > 0.0){
                for (unsigned int j = 0; j < r_geometry.size(); j++) {
                    r_geometry[j].SetLock();
                    r_geometry[j].SetValue(NORMAL, surface_normal);
                    r_geometry[j].SetValue(LOWER_SURFACE, true);
                    r_geometry[j].UnSetLock();
                }
            }
        }
    }
    KRATOS_CATCH("");
}

// This function computes the local wake normal at each trailing edge node
// by avaraging the local wake normals of the surrounding conditions
void Define3DWakeOperation::ComputeAndSaveLocalWakeNormal() const
{
    KRATOS_TRY;   
    array_1d<double,3> coordinates1;
    array_1d<double,3> coordinates2 ;
    array_1d<double,3> local_span_direction;
    
    for (auto& r_cond : mpTrailingEdgeModelPart->Conditions()){
        auto& r_geometry = r_cond.GetGeometry();

        // To ensure that local_span has the same direction as mSpanDirection
        const auto& test_start_point = r_geometry[0].Coordinates();
        const auto& test_end_point   = r_geometry[1].Coordinates();
        const auto& test_local_span_direction = test_end_point - test_start_point;
        const double test_projection = inner_prod(test_local_span_direction, mSpanDirection);

        if (!(test_projection > 0.0))
        {
            coordinates1 = r_geometry[1].Coordinates();
            coordinates2 = r_geometry[0].Coordinates();
        }else
        {
            coordinates1 = r_geometry[0].Coordinates();
            coordinates2 = r_geometry[1].Coordinates();
        }
        local_span_direction = coordinates2 - coordinates1;

        array_1d<double, 3> local_wake_normal = ZeroVector(3);
        MathUtils<double>::CrossProduct(local_wake_normal, mWakeDirection, local_span_direction);

        for (unsigned int i = 0; i < r_geometry.size(); i++){
            r_geometry[i].GetValue(WAKE_NORMAL) += local_wake_normal;
        }
    }

    for (auto& r_node: mpTrailingEdgeModelPart->Nodes()){
        auto& local_wake_normal = r_node.GetValue(WAKE_NORMAL);
        const double norm = MathUtils<double>::Norm3(local_wake_normal);
        local_wake_normal /= norm;
        r_node.SetValue(WAKE_NORMAL, local_wake_normal);
    }
    KRATOS_CATCH("");
}

// This function creates the wake surface automatically by shedding it from the
// trailing edge/s in the direction of the free stream velocity (mWakeDirection).
// The user can decide how much distance is to be shedded, the element size
// of the wake surface in the wake direction and the grow factor of the element size. 
// Note that the element size in span direction is predetermined by the size of the 
// conditions constituting the trailing edge.
void Define3DWakeOperation::ShedWakeSurfaceFromTheTrailingEdge() const
{
    KRATOS_TRY;
    KRATOS_WARNING_IF("Define3DWakeOperation", mWakeSTLFileName == "") << "Generating wake automatically. Set \'wake_stl_file_name\' instead" << std::endl;
    KRATOS_INFO_IF("Define3DWakeOperation", mEchoLevel > 0) << "...Shedding wake from the trailing_edge_model_part" << std::endl;
    
    const auto p_elem_prop = mpBodyModelPart->pGetProperties(0);
    IndexType node_index = 0;
    IndexType element_index = 0;
    array_1d<double,3> start_point;
    array_1d<double,3> end_point  ;

    for (auto& r_cond : mpTrailingEdgeModelPart->Conditions()) {
        auto& r_geometry = r_cond.GetGeometry();

        // To ensure that local_span has the same direction as mSpanDirection
        const auto& test_start_point = r_geometry[0].Coordinates();
        const auto& test_end_point   = r_geometry[1].Coordinates();
        const auto& test_local_span_direction = test_end_point - test_start_point;
        const double test_projection = inner_prod(test_local_span_direction, mSpanDirection);

        if (!(test_projection > 0.0))
        {
            start_point = r_geometry[1].Coordinates();
            end_point   = r_geometry[0].Coordinates();
        }else
        {
            start_point = r_geometry[0].Coordinates();
            end_point   = r_geometry[1].Coordinates();
        }

        // Check if any node is on the wing root
        if (r_geometry[0].GetValue(WING_ROOT) != 0) {
            auto local_vector = end_point - start_point;
            start_point = start_point - (local_vector*mShedProjectionRootEdge/norm_2(local_vector));
        }

        // Create the wake surface
        double current_element_size = mShedWakeElementSize;
        double accumulated_distance = 0.0;

        array_1d<double,3> next_start = start_point + current_element_size * mWakeDirection;
        array_1d<double,3> next_end   = end_point   + current_element_size * mWakeDirection;

        CreateWakeSurfaceNodesAndElements(*mpWakeModelPart, node_index, start_point, end_point, next_start, next_end, element_index, p_elem_prop);

        while (accumulated_distance + current_element_size <= mShedWakeLength) {
            accumulated_distance += current_element_size;
            current_element_size *= mShedGrowFactor;

            start_point = next_start;
            end_point   = next_end;
            
            next_start = start_point + current_element_size * mWakeDirection;
            next_end   = end_point   + current_element_size * mWakeDirection;
            
            CreateWakeSurfaceNodesAndElements(*mpWakeModelPart, node_index, start_point, end_point, next_start, next_end, element_index, p_elem_prop);
        }
    }
    KRATOS_CATCH("");
}

void Define3DWakeOperation::CreateWakeSurfaceNodesAndElements(
    ModelPart& rModelPart,
    IndexType& rNode_index,
    const array_1d<double, 3>& rCoordinates1,
    const array_1d<double, 3>& rCoordinates2,
    const array_1d<double, 3>& rCoordinates3,
    const array_1d<double, 3>& rCoordinates4,
    IndexType& rElement_index,
    const Properties::Pointer pElemProp) const
{
    KRATOS_TRY;
    const auto nodes_ids = CreateWakeSurfaceNodes(rModelPart,
        rNode_index, rCoordinates1, rCoordinates2, rCoordinates3, rCoordinates4);
        
        const double normal_projection = ComputeFaceNormalProjectionToWakeNormal(
            rCoordinates1, rCoordinates2, rCoordinates3, rCoordinates4);
            
            CreateWakeSurfaceElements(rModelPart, normal_projection, rElement_index, nodes_ids, pElemProp);
    KRATOS_CATCH("");   
}

std::array<ModelPart::IndexType, 4> Define3DWakeOperation::CreateWakeSurfaceNodes(
    ModelPart& rModelPart,
    IndexType& rNode_index,
    const array_1d<double, 3>& rCoordinates1,
    const array_1d<double, 3>& rCoordinates2,
    const array_1d<double, 3>& rCoordinates3,
    const array_1d<double, 3>& rCoordinates4) const
{
    KRATOS_TRY;
    const auto p_node1 = rModelPart.CreateNewNode(++rNode_index, rCoordinates1[0], rCoordinates1[1], rCoordinates1[2]);
    const auto p_node2 = rModelPart.CreateNewNode(++rNode_index, rCoordinates2[0], rCoordinates2[1], rCoordinates2[2]);
    const auto p_node3 = rModelPart.CreateNewNode(++rNode_index, rCoordinates3[0], rCoordinates3[1], rCoordinates3[2]);
    const auto p_node4 = rModelPart.CreateNewNode(++rNode_index, rCoordinates4[0], rCoordinates4[1], rCoordinates4[2]);
    
    return {p_node1->Id(), p_node2->Id(), p_node3->Id(), p_node4->Id()};
    KRATOS_CATCH("");
}

double Define3DWakeOperation::ComputeFaceNormalProjectionToWakeNormal(
    const array_1d<double, 3>& rCoordinates1,
    const array_1d<double, 3>& rCoordinates2,
    const array_1d<double, 3>& rCoordinates3,
    const array_1d<double, 3>& rCoordinates4) const
{
    KRATOS_TRY;
    const auto& side1 = rCoordinates2 - rCoordinates1;
    const auto& side2 = rCoordinates3 - rCoordinates1;
    array_1d<double, 3> face_normal = ZeroVector(3);
    MathUtils<double>::CrossProduct(face_normal, side1, side2);
    return inner_prod(face_normal, mWakeNormal);
    KRATOS_CATCH("");
}

void Define3DWakeOperation::CreateWakeSurfaceElements(ModelPart& rModelPart,
                                                    const double normal_projection,
                                                    IndexType& rElement_index,
                                                    const std::array<ModelPart::IndexType, 4>& rNodes_ids,
                                                    const Properties::Pointer pElemProp) const
{
    KRATOS_TRY;
    if(normal_projection > 0.0){
        const std::vector<ModelPart::IndexType> elem_nodes_1{rNodes_ids[0], rNodes_ids[1], rNodes_ids[2]};
        const std::vector<ModelPart::IndexType> elem_nodes_2{rNodes_ids[1], rNodes_ids[3], rNodes_ids[2]};
        rModelPart.CreateNewElement("Element3D3N", ++rElement_index, elem_nodes_1, pElemProp);
        rModelPart.CreateNewElement("Element3D3N", ++rElement_index, elem_nodes_2, pElemProp);
    }
    else{
        const std::vector<ModelPart::IndexType> elem_nodes_1{rNodes_ids[0], rNodes_ids[2], rNodes_ids[1]};
        const std::vector<ModelPart::IndexType> elem_nodes_2{rNodes_ids[1], rNodes_ids[2], rNodes_ids[3]};
        rModelPart.CreateNewElement("Element3D3N", ++rElement_index, elem_nodes_1, pElemProp);
        rModelPart.CreateNewElement("Element3D3N", ++rElement_index, elem_nodes_2, pElemProp);
    }
    KRATOS_CATCH("");
}

// To load a stl mesh
void Define3DWakeOperation::LoadSTL() const
{
    KRATOS_TRY;
    KRATOS_INFO_IF("Define3DWakeOperation", mEchoLevel > 0) << "...Reading wake from stl file" << std::endl;
    
    // Operation parameters
    Parameters reading_parameters = Parameters(R"(
    {
        "open_mode"       : "read",
        "new_entity_type" : "element"
    })" );
    
    StlIO stl_read(mWakeSTLFileName, reading_parameters);
    stl_read.ReadModelPart(*mpWakeModelPart);
    KRATOS_CATCH("");
}

// To move the wake model part
void Define3DWakeOperation::MoveWakeModelPart() const
{   
    KRATOS_TRY;
    KRATOS_INFO_IF("Define3DWakeOperation", mEchoLevel > 0) << "...Moving wake_model_part" << std::endl;
    
    Parameters moving_parameters = Parameters(R"(
    {
        "origin": []
    })" );
        
    moving_parameters["origin"].SetVector(mWakeTraslationDirection);
            
    MoveModelPartProcess MoveModelPart(*mpWakeModelPart, moving_parameters);
    MoveModelPart.Execute();
    KRATOS_CATCH("");
}

// To visualize the wake
void Define3DWakeOperation::VisualizeWake() const
{
    KRATOS_TRY;
        KRATOS_INFO_IF("Define3DWakeOperation", mEchoLevel > 0) << "...Saving vtk wake_model_part in 'wake_output' folder" << std::endl;
        Parameters vtk_parameters = Parameters(R"(
        {
            "model_part_name"                    : "wake_model_part",
            "output_control_type"                : "step",
            "output_interval"                    : 1.0,
            "file_format"                        : "binary",
            "output_precision"                   : 7,
            "output_sub_model_parts"             : false,
            "output_path"                        : "wake_output",
            "save_output_files_in_folder"        : true,
            "write_deformed_configuration"       : true
        })" );
            
        VtkOutput vtk_oi(*mpWakeModelPart, vtk_parameters);
        vtk_oi.PrintOutput();
    KRATOS_CATCH("");
}

// This function checks which elements are cut by the wake and marks them as
// wake elements
void Define3DWakeOperation::MarkWakeElements() const
{
    KRATOS_TRY;
    KRATOS_INFO_IF("Define3DWakeOperation", mEchoLevel > 0) << "...Selecting wake elements" << std::endl;

    BuiltinTimer timer;
    CalculateDiscontinuousDistanceToSkinProcess<3> distance_calculator(*mpRootModelPart, *mpWakeModelPart);
    distance_calculator.Execute();

    KRATOS_INFO_IF("Define3DWakeOperation", mEchoLevel > 0) << "...Distance_calculator took " << timer.ElapsedSeconds() << " [sec]" << std::endl;

    // Variable to store element ids
    moodycamel::ConcurrentQueue<std::size_t> wake_elements_ordered_ids_concurrent_queue;
    moodycamel::ConcurrentQueue<std::size_t> trailing_edge_elements_ordered_ids_concurrent_queue;

    block_for_each(mpRootModelPart->Elements(), [&](Element& rElement)
    {
        // Check if the element is touching the trailing edge
        auto& r_geometry = rElement.GetGeometry();
        CheckIfTrailingEdgeElement(rElement, r_geometry, trailing_edge_elements_ordered_ids_concurrent_queue);

        // Mark wake elements, save their ids, save the elemental distances in
        // the element and in the nodes, and save wake nodes ids
        if (rElement.Is(TO_SPLIT) || rElement.GetValue(WAKE))
        {
            // Mark wake elements
            rElement.SetValue(WAKE, true);

            wake_elements_ordered_ids_concurrent_queue.enqueue(rElement.Id());

            // Save elemental distances in the element
            array_1d<double, 4> wake_elemental_distances = ZeroVector(4);
            wake_elemental_distances = rElement.GetValue(ELEMENTAL_DISTANCES);
            rElement.SetValue(WAKE_ELEMENTAL_DISTANCES, wake_elemental_distances);

            // Save elemental distances in the nodes
            for (unsigned int j = 0; j < wake_elemental_distances.size(); j++)
            {
                if (std::abs(wake_elemental_distances[j]) < mWakeDistanceTolerance)
                {
                    if (wake_elemental_distances[j] < 0.0)
                    {
                        wake_elemental_distances[j] = -mWakeDistanceTolerance;
                    }
                    else
                    {
                        wake_elemental_distances[j] = mWakeDistanceTolerance;
                    }
                }
                r_geometry[j].SetLock();
                r_geometry[j].SetValue(WAKE_DISTANCE, wake_elemental_distances[j]);
                r_geometry[j].UnSetLock();
            }
        }
    });

    // Variable to store element ids
    std::vector<std::size_t> wake_elements_ordered_ids;
    std::vector<std::size_t> trailing_edge_elements_ordered_ids;

    bool found_value;
    std::size_t id;

    while ( (found_value = wake_elements_ordered_ids_concurrent_queue.try_dequeue(id)) ) {
        wake_elements_ordered_ids.push_back(id);
    }

    while ( (found_value = trailing_edge_elements_ordered_ids_concurrent_queue.try_dequeue(id)) ) {
        trailing_edge_elements_ordered_ids.push_back(id);
    }

    // Add the trailing edge elements to the trailing_edge_sub_model_part
    AddTrailingEdgeAndWakeElements(wake_elements_ordered_ids, trailing_edge_elements_ordered_ids);

    KRATOS_INFO_IF("Define3DWakeOperation", mEchoLevel > 0) << "...Selecting wake elements finished" << std::endl;
    KRATOS_CATCH("");
}

// This function checks if the element is touching the trailing edge
void Define3DWakeOperation::CheckIfTrailingEdgeElement(
    Element& rElement,
    const Geometry<NodeType>& rGeometry,
    moodycamel::ConcurrentQueue<std::size_t>& rTrailingEdgeElementsOrderedIds) const
{
    KRATOS_TRY;
    // Loop over element nodes
    for (unsigned int i = 0; i < rGeometry.size(); i++)
    {
        // Elements touching the trailing edge are trailing edge elements
        if (rGeometry[i].GetValue(TRAILING_EDGE))
        {
            rElement.SetValue(TRAILING_EDGE, true);
            rTrailingEdgeElementsOrderedIds.enqueue(rElement.Id());
            break;
        }
    }

    // Here, the elements of the trailing edge are marked when it is blunt.
    // Think of a better way to do this.
    if (mpBluntTESurfaceModelPart && rElement.GetValue(TRAILING_EDGE))
    {
        unsigned int number_of_blunt_nodes = 0;
        std::unordered_set<std::size_t> indices_set(mBluntIds.begin(), mBluntIds.end());
        // Loop over element nodes
        for (unsigned int i = 0; i < rGeometry.size(); i++)
        {
            if (indices_set.find(rGeometry[i].Id()) != indices_set.end())
            {
                number_of_blunt_nodes += 1;
            }
        }
        if (number_of_blunt_nodes > 2)
        {   
            rElement.SetValue(WAKE, true);
        }
    }
    KRATOS_CATCH("");
}

// This function adds the trailing edge elements in the
// trailing_edge_sub_model_part
void Define3DWakeOperation::AddTrailingEdgeAndWakeElements(std::vector<std::size_t>& rWakeElementsOrderedIds,
                                                         std::vector<std::size_t>& rTrailingEdgeElementsOrderedIds) const
{
    KRATOS_TRY;
    std::sort(rWakeElementsOrderedIds.begin(),
              rWakeElementsOrderedIds.end());
    mpRootModelPart->GetSubModelPart("wake_elements_model_part").AddElements(rWakeElementsOrderedIds);

    std::sort(rTrailingEdgeElementsOrderedIds.begin(),
              rTrailingEdgeElementsOrderedIds.end());
    ModelPart& trailing_edge_sub_model_part =
        mpRootModelPart->GetSubModelPart("trailing_edge_elements_model_part");
    trailing_edge_sub_model_part.AddElements(rTrailingEdgeElementsOrderedIds);

    std::vector<std::size_t> trailing_edge_nodes_ordered_ids;
    trailing_edge_nodes_ordered_ids.reserve(trailing_edge_sub_model_part.NumberOfElements());

    // Add also the nodes of the elements touching the trailing edge
    for (auto& r_elem : trailing_edge_sub_model_part.Elements())
    {
        for (auto& r_node : r_elem.GetGeometry()){
            trailing_edge_nodes_ordered_ids.push_back(r_node.Id());
        }
    }

    std::sort(trailing_edge_nodes_ordered_ids.begin(),
    trailing_edge_nodes_ordered_ids.end());
    trailing_edge_sub_model_part.AddNodes(trailing_edge_nodes_ordered_ids);
    KRATOS_CATCH("");
}

// This function recomputes the wake distances from the nodes belonging to the
// elements. These distances are used later to decide which elements are KUTTA,
// which are WAKE (STRUCTURE), and which are NORMAL
void Define3DWakeOperation::RecomputeNodalDistancesToWakeOrWingLowerSurface() const
{
    KRATOS_TRY;
    ModelPart& trailing_edge_sub_model_part =
        mpRootModelPart->GetSubModelPart("trailing_edge_elements_model_part");

    block_for_each(trailing_edge_sub_model_part.Nodes(), [&](Node& r_node)
    {
        NodeType::Pointer p_closest_te_node = *mpTrailingEdgeModelPart->NodesBegin().base();
        FindClosestTrailingEdgeNode(p_closest_te_node, r_node);
        RecomputeDistance(p_closest_te_node, r_node);
    });
    KRATOS_CATCH("");
}

// This function finds the closest trailing edge node to the given point
void Define3DWakeOperation::FindClosestTrailingEdgeNode(NodeType::Pointer& pClosest_te_node,
                                                      const array_1d<double, 3>& rPoint) const
{
    KRATOS_TRY;
    double min_distance_to_te = std::numeric_limits<double>::max();
    for (auto& r_te_node : mpTrailingEdgeModelPart->Nodes())
    {
        const auto& distance_vector = rPoint - r_te_node;
        const double distance_to_te = inner_prod(distance_vector, distance_vector);
        if (distance_to_te < min_distance_to_te)
        {
            min_distance_to_te = distance_to_te;
            pClosest_te_node = &r_te_node;
        }
    }
    KRATOS_CATCH("");
}

// This function recomputes the distances from the given node to the wake or to
// the wing lower surface, depending whether the node is behind or in front of
// the trailing edge according to the wake direction (free stream velocity
// direction).
void Define3DWakeOperation::RecomputeDistance(NodeType::Pointer& pClosest_te_node,
                                            NodeType& rNode) const
{
    KRATOS_TRY;
    // Compute the distance vector from the closest trailing edge node to the
    // current node
    const auto& distance_vector = rNode - *pClosest_te_node;

    // Compute the distance in the free stream direction
    const double free_stream_direction_distance = inner_prod(distance_vector, mWakeDirection);

    //  For the nodes in front of the selected trailing edge node, the
    //  distance is computed according to the lower surface normal
    if(free_stream_direction_distance < 0.0){
        double distance = inner_prod(distance_vector, pClosest_te_node->GetValue(NORMAL));
        // Nodes close to the wing are given a negative distance
        if(std::abs(distance) < mWakeDistanceTolerance){
            distance = - mWakeDistanceTolerance;
        }
        rNode.SetValue(WAKE_DISTANCE, distance);
    }
    // For the nodes behind of the selected trailing edge node, the
    // distance is computed according to the wake normal
    else{
        const auto& local_wake_normal = pClosest_te_node->GetValue(WAKE_NORMAL);
        double distance = inner_prod(distance_vector, local_wake_normal);
        // Nodes slightly below and above the wake are given a positive distance (wake down)
        if(std::abs(distance) < mWakeDistanceTolerance){
            distance = mWakeDistanceTolerance;
        }
        rNode.SetValue(WAKE_DISTANCE, distance);
    }

    if (rNode.GetValue(TRAILING_EDGE))
    {
        rNode.SetValue(WAKE_DISTANCE, mWakeDistanceTolerance);
    }
    KRATOS_CATCH("");
}

// This function selects the kutta elements. Kutta elements are touching the
// trailing edge from below
void Define3DWakeOperation::MarkKuttaElements() const
{
    KRATOS_TRY;
    KRATOS_INFO_IF("Define3DWakeOperation", mEchoLevel > 0) << "...Selecting kutta elements" << std::endl;
    ModelPart& trailing_edge_sub_model_part =
    mpRootModelPart->GetSubModelPart("trailing_edge_elements_model_part");

    std::vector<std::size_t> wake_elements_ordered_ids;

    // TO DISCUSS: Is it worth it to run this loop in parallel?
    // So far not much different in terms of speed has been noted.
    block_for_each(trailing_edge_sub_model_part.Elements(), [&](Element& rElement)
    {
        auto& r_geometry = rElement.GetGeometry();
        // Counting number of trailing edge nodes in the element
        unsigned int number_of_te_nodes = CountNumberOfTrailindEdgeNodesInElement(r_geometry);

        KRATOS_ERROR_IF(number_of_te_nodes < 1)
            << "Number of trailing edge nodes must be 1 or larger. Element Id: "
            << rElement.Id() << " number_of_te_nodes = " << number_of_te_nodes
            << std::endl;

        // Checking only distances from nodes that are not trailing edge
        unsigned int number_of_nodes_with_negative_distance = 0;
        unsigned int number_of_nodes_with_positive_distance = 0;

        CountNumberOfPositiveAndNegativeDistances(
            r_geometry, number_of_nodes_with_negative_distance,
            number_of_nodes_with_positive_distance);
            
        SelectElementType(rElement, r_geometry, number_of_te_nodes,
                          number_of_nodes_with_negative_distance,
                          number_of_nodes_with_positive_distance);
    });

    // Remove elements that were wake and now are either kutta or normal.
    ModelPart& wake_sub_model_part = mpRootModelPart->GetSubModelPart("wake_elements_model_part");
    wake_sub_model_part.RemoveElements(TO_ERASE);
    KRATOS_INFO_IF("Define3DWakeOperation", mEchoLevel > 0) << "...Selecting kutta elements finished" << std::endl;
    KRATOS_CATCH("");
}

// This function returns the number of trailing edge nodes
unsigned int Define3DWakeOperation::CountNumberOfTrailindEdgeNodesInElement(
    const Geometry<NodeType>& rGeometry) const
{
    KRATOS_TRY;
    unsigned int number_of_te_nodes = 0;
    
    for (unsigned int j = 0; j < rGeometry.size(); j++)
    {
        const auto& r_node = rGeometry[j];
        if (r_node.GetValue(TRAILING_EDGE))
        {
            number_of_te_nodes += 1;
        }
    }
    
    return number_of_te_nodes;
    KRATOS_CATCH("");
}

// This function returns the number of non trailing edge nodes with positive and
// negative distance
void Define3DWakeOperation::CountNumberOfPositiveAndNegativeDistances(
    const Geometry<NodeType>& rGeometry,
    unsigned int& number_of_nodes_with_negative_distance,
    unsigned int& number_of_nodes_with_positive_distance) const
{
    KRATOS_TRY;
    for (unsigned int j = 0; j < rGeometry.size(); j++){
        const auto& r_node = rGeometry[j];
        if (!r_node.GetValue(TRAILING_EDGE)){
            const auto& distance = r_node.GetValue(WAKE_DISTANCE);
            if(distance < 0.0){
                number_of_nodes_with_negative_distance += 1;
            }
            else{
                number_of_nodes_with_positive_distance +=1;
            }
        }
    }
    KRATOS_CATCH("");
}

// This function selects the element type (wake, kutta, or normal) according to
// the number of non trailing edge nodes with positive and negative distance
void Define3DWakeOperation::SelectElementType(Element& rElement,
                                            const Geometry<NodeType>& rGeometry,
                                            const unsigned int number_of_te_nodes,
                                            const unsigned int number_of_nodes_with_negative_distance,
                                            const unsigned int number_of_nodes_with_positive_distance) const
{
    KRATOS_TRY;
    // Wake structure elements (cut)
    if ((number_of_nodes_with_positive_distance > 0 &&
        number_of_nodes_with_negative_distance > 0 && rElement.GetValue(WAKE)) || rElement.GetValue(WAKE))
        {
        rElement.Set(STRUCTURE);
        BoundedVector<double, 4> wake_elemental_distances = ZeroVector(4);
        for (unsigned int j = 0; j < rGeometry.size(); j++)
        {
            const auto& r_node = rGeometry[j];
            const auto& distance = r_node.GetValue(WAKE_DISTANCE);
            wake_elemental_distances[j] = distance;
        }
        rElement.SetValue(WAKE_ELEMENTAL_DISTANCES, wake_elemental_distances);
    }

    // Kutta elements (below). Kutta elements have all non_te_nodes with
    // negative distance:
    // 1 te node  -> 3 non te nodes with negative distance > 2 = 3 - 1
    // 2 te nodes -> 2 non te nodes with negative distance > 1 = 3 - 2
    else if (number_of_nodes_with_negative_distance > 3 - number_of_te_nodes)
    {
        rElement.SetValue(KUTTA, true);
        rElement.SetValue(WAKE, false);
        rElement.Set(TO_ERASE, true);
    }

    // Normal elements (above). Normal elements have all nodes with positive
    // distance.
    else
    {
        rElement.SetValue(WAKE, false);
        rElement.Set(TO_ERASE, true);
    }
    KRATOS_CATCH("");
}

// This function saves the local wake normal in the element to be used later
// inside the element to apply the wake conditions
void Define3DWakeOperation::SaveLocalWakeNormalInElements() const
{
    KRATOS_TRY;
    ModelPart& wake_sub_model_part = mpRootModelPart->GetSubModelPart("wake_elements_model_part");

    for (auto& r_elem : wake_sub_model_part.Elements()){
        NodeType::Pointer p_closest_te_node = *mpTrailingEdgeModelPart->NodesBegin().base();
        FindClosestTrailingEdgeNode(p_closest_te_node, r_elem.GetGeometry().Center());

        r_elem.GetValue(WAKE_NORMAL) = p_closest_te_node->GetValue(WAKE_NORMAL);
    }
    KRATOS_CATCH("");
}

void Define3DWakeOperation::AddWakeNodesToWakeModelPart() const
{
    KRATOS_TRY;
    ModelPart& wake_sub_model_part =
            mpRootModelPart->GetSubModelPart("wake_elements_model_part");

    std::vector<std::size_t> wake_nodes_ordered_ids;
    for (auto& r_element : wake_sub_model_part.Elements()){
        for (unsigned int i = 0; i < r_element.GetGeometry().size(); i++){
            r_element.GetGeometry()[i].SetValue(WAKE, true);
            wake_nodes_ordered_ids.push_back(r_element.GetGeometry()[i].Id());
        }
    }

    std::sort(wake_nodes_ordered_ids.begin(),
              wake_nodes_ordered_ids.end());
              wake_sub_model_part.AddNodes(wake_nodes_ordered_ids);
    KRATOS_CATCH("");
}

class KdTreeNode 
{
    public:
        Node* pNode = nullptr;
        KdTreeNode* pLeft = nullptr;
        KdTreeNode* pRight = nullptr;

        ~KdTreeNode() {
            delete pLeft;
            delete pRight;
        }
};

// This class implements a kd-tree to detect trailing edge groups
class TrueKdTree 
{
    public:
        KdTreeNode* mpRoot = nullptr;

        void Build(std::vector<Node*>& r_nodes) {
            mpRoot = BuildRecursive(r_nodes, 0, r_nodes.size(), 0);
        }

        ~TrueKdTree() {
            delete mpRoot;
        }

        std::vector<Node*> RadiusSearch(Node* p_node, double radius) const {
            std::vector<Node*> neighbors;
            array_1d<double, 3> target = p_node->Coordinates();
            RadiusSearchRecursive(mpRoot, target, radius * radius, 0, neighbors);
            return neighbors;
        }

    private:
        static KdTreeNode* BuildRecursive(std::vector<Node*>& nodes, int start, int end, int depth) {
            if (start >= end) return nullptr;

            int axis = depth % 3;
            int mid = start + (end - start) / 2;

            std::nth_element(nodes.begin() + start, nodes.begin() + mid, nodes.begin() + end,
                [axis](Node* a, Node* b) {
                    return a->Coordinates()[axis] < b->Coordinates()[axis];
                });

            KdTreeNode* node = new KdTreeNode();
            node->pNode = nodes[mid];
            node->pLeft = BuildRecursive(nodes, start, mid, depth + 1);
            node->pRight = BuildRecursive(nodes, mid + 1, end, depth + 1);

            return node;
        }

        static void RadiusSearchRecursive(KdTreeNode* pNode, const array_1d<double, 3>& target, double radius2, int depth, std::vector<Node*>& r_neighbors) {
            if (!pNode) return;

            const int axis = depth % 3;
            const array_1d<double, 3>& node_coords = pNode->pNode->Coordinates();

            double dx = node_coords[0] - target[0];
            double dy = node_coords[1] - target[1];
            double dz = node_coords[2] - target[2];
            double distance2 = dx*dx + dy*dy + dz*dz;

            if (distance2 <= radius2) {
                r_neighbors.push_back(pNode->pNode);
            }

            double diff = target[axis] - node_coords[axis];
            double diff2 = diff * diff;

            if (diff < 0) {
                RadiusSearchRecursive(pNode->pLeft, target, radius2, depth + 1, r_neighbors);
                if (diff2 < radius2) {
                    RadiusSearchRecursive(pNode->pRight, target, radius2, depth + 1, r_neighbors);
                }
            } else {
                RadiusSearchRecursive(pNode->pRight, target, radius2, depth + 1, r_neighbors);
                if (diff2 < radius2) {
                    RadiusSearchRecursive(pNode->pLeft, target, radius2, depth + 1, r_neighbors);
                }
            }
        }
};

// This function adds conditions to the trailing edge model part and
// finds the tip and root nodes, if necessary
void Define3DWakeOperation::AddTrailingEdgeConditionsAndFindRootAndTipNodes() const
{
    KRATOS_TRY;
    bool assign_conditions = (mpTrailingEdgeModelPart->NumberOfConditions() < 1);
    if (assign_conditions || mpTipPointsModelPart || mpRootPointsModelPart )
    {
        KRATOS_WARNING_IF("Define3DWakeOperation", assign_conditions) << "Assigning conditions automatically. Check the trailing_edge_model_part" << std::endl;
        KRATOS_WARNING_IF("Define3DWakeOperation", mpTipPointsModelPart)  << "Detecting tip nodes automatically. Set \'tip_points_model_part_name\' instead" << std::endl;
        KRATOS_WARNING_IF("Define3DWakeOperation", mpRootPointsModelPart) << "Detecting root nodes automatically. Set \'root_points_model_part_name\' instead" << std::endl;

        // distance_threshold controls how to detect more than one trailing edge. constant_threshold=3.0 usually works well.
        // When necessary, typically for small domains, a low threshold is needed, which can be configured with a low mShedWakeElementSize.
        const double constant_threshold = 3.0;
        const double distance_threshold = constant_threshold * mShedWakeElementSize;
        int reference_id = mpRootModelPart->NumberOfConditions();
        const auto p_cond_prop = mpTrailingEdgeModelPart->pGetProperties(0);
    
        std::vector<Node*> all_nodes;
        all_nodes.reserve(mpTrailingEdgeModelPart->NumberOfNodes());
        for (auto& r_node : mpTrailingEdgeModelPart->Nodes()) {
            all_nodes.push_back(&r_node);
        }
    
        TrueKdTree kd_tree;
        kd_tree.Build(all_nodes);
    
        std::vector<std::vector<Node*>> trailing_edge_groups;
        std::unordered_set<int> visited;
        visited.reserve(all_nodes.size());
    
        for (auto* p_node : all_nodes) {
            if (visited.count(p_node->Id()) != 0) {
                continue;
            }
    
            std::vector<Node*> current_group;
            std::queue<Node*> queue;
            queue.push(p_node);
            visited.insert(p_node->Id());
    
            while (!queue.empty()) {
                Node* p_current = queue.front();
                queue.pop();
                current_group.push_back(p_current);
    
                auto neighbors = kd_tree.RadiusSearch(p_current, distance_threshold);
                for (auto* p_neighbor : neighbors) {
                    if (visited.count(p_neighbor->Id()) == 0) {
                        visited.insert(p_neighbor->Id());
                        queue.push(p_neighbor);
                    }
                }
            }
    
            trailing_edge_groups.push_back(std::move(current_group));
        }
        
        int tip_nodes = 0;
        int root_nodes = 0;
        for (auto& group : trailing_edge_groups) {
    
            if (mpTipPointsModelPart || mpRootPointsModelPart)
            {
                auto [min_it, max_it] = std::minmax_element(
                    group.begin(), group.end(),
                    [this](const Node* a, const Node* b) {
                        return inner_prod(a->Coordinates(), mSpanDirection) <
                               inner_prod(b->Coordinates(), mSpanDirection);
                    }
                );
                
                if (mpTipPointsModelPart &&
                    !((*max_it)->GetValue(WING_TIP)) &&
                    !((*max_it)->GetValue(WING_ROOT)))
                {
                    (*max_it)->SetValue(WING_TIP, true);
                    tip_nodes += 1;
                }
                
                if (mpRootPointsModelPart &&
                    !((*min_it)->GetValue(WING_TIP)) &&
                    !((*min_it)->GetValue(WING_ROOT)))
                {
                    (*min_it)->SetValue(WING_ROOT, true);
                    root_nodes += 1;
                }
                
            }
            
            if (assign_conditions) {
                std::sort(group.begin(), group.end(),
                    [this](const Node* a, const Node* b) {
                        return inner_prod(a->Coordinates(), mSpanDirection) <
                               inner_prod(b->Coordinates(), mSpanDirection);
                    }
                );
    
                for (std::size_t i = 1; i < group.size(); ++i) {
                    ++reference_id;
                    std::vector<ModelPart::IndexType> cond{group[i - 1]->Id(), group[i]->Id()};
                    mpTrailingEdgeModelPart->CreateNewCondition("LineCondition3D2N", reference_id, cond, p_cond_prop);
                }
            }
        }
    
        KRATOS_INFO_IF("Define3DWakeOperation", mEchoLevel > 0) << "...Detected " << trailing_edge_groups.size() << " trailing edge groups" << std::endl;
        KRATOS_INFO_IF("Define3DWakeOperation", (mEchoLevel > 0 && (mpTipPointsModelPart))) << "...Detected " << tip_nodes << " tip nodes" << std::endl;
        KRATOS_INFO_IF("Define3DWakeOperation", (mEchoLevel > 0 && (mpRootPointsModelPart))) << "...Detected " << root_nodes << " root nodes" << std::endl;
        KRATOS_INFO_IF("Define3DWakeOperation", (mEchoLevel > 0 && assign_conditions)) << "...Assigned " << mpTrailingEdgeModelPart->NumberOfConditions() << " conditions" << std::endl;
    }

    if (mpTipPointsModelPart)
    {
        int tip_nodes = 0;
        for (auto& r_node : mpTipPointsModelPart->Nodes()) {
            r_node.SetValue(WING_TIP, true);
            tip_nodes += 1;
        }
        KRATOS_INFO_IF("Define3DWakeOperation", mEchoLevel > 0) << "...Detected " << tip_nodes << " tip nodes" << std::endl;
    }
    
    if (mpRootPointsModelPart)
    {
        int root_nodes = 0;
        for (auto& r_node : mpRootPointsModelPart->Nodes()) {
            r_node.SetValue(WING_ROOT, true);
            root_nodes += 1;
        }
        KRATOS_INFO_IF("Define3DWakeOperation", mEchoLevel > 0) << "...Detected " << root_nodes << " root nodes" << std::endl;
    }
       
    KRATOS_CATCH("");
}

const Parameters Define3DWakeOperation::GetDefaultParameters() const
{
    KRATOS_TRY;
    const Parameters default_parameters = Parameters(R"({
        "trailing_edge_model_part_name"    : "",
        "body_model_part_name"             : "",
        "upper_surface_model_part_name"    : "",
        "lower_surface_model_part_name"    : "",
        "root_points_model_part_name"      : "",
        "tip_points_model_part_name"       : "",
        "blunt_te_surface_model_part_name" : "",
        "wake_stl_file_name"               : "",
        "wake_normal"                      : [0.0,0.0,1.0],
        "wake_translation_direction"       : [0.0,0.0,0.0],
        "visualize_wake_vtk"               : false,
        "shed_wake_from_trailing_edge"     : false,
        "shed_wake_length"                 : 12.5,
        "shed_wake_element_size"           : 0.2,
        "shed_wake_grow_factor"            : 1.05,
        "shed_wake_projection_root_edge"   : 0.0,
        "wake_distance_tolerance"          : 1e-9,
        "echo_level"                       : 0
    })");
    return default_parameters;
    KRATOS_CATCH("");
}

} // namespace Kratos.
