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

// Project includes
#include "includes/model_part.h"
#include "processes/process.h"
#include "compressible_potential_flow_application_variables.h"
#include "define_3d_wake_process.h"
#include "utilities/variable_utils.h"
#include "custom_utilities/potential_flow_utilities.h"
#include "processes/calculate_distance_to_skin_process.h"
#include "utilities/builtin_timer.h"
#include "custom_processes/move_model_part_process.h"
#include "input_output/vtk_output.h"
#include "input_output/stl_io.h"
#include "containers/model.h"
#include <queue>

namespace Kratos {

// Constructor for Define3DWakeProcess Process
Define3DWakeProcess::Define3DWakeProcess(ModelPart& rTrailingEdgeModelPart,
                                         ModelPart& rBodyModelPart,
                                         Parameters ThisParameters)
    : Process(),
      mrTrailingEdgeModelPart(rTrailingEdgeModelPart),
      mrBodyModelPart(rBodyModelPart)
{
    Parameters default_parameters = Parameters(R"(
    {
        "wake_normal"                      : [0.0,0.0,1.0],
        "wake_stl_file_name"               : "",
        "wake_translation_direction"       : [0.0,0.0,0.0],
        "tolerance"                        : 1e-9,
        "visualize_wake_vtk"               : false,
        "upper_surface_model_part_name"    : "",
        "lower_surface_model_part_name"    : "",
        "root_points_model_part_name"      : "",
        "tip_points_model_part_name"       : "",
        "blunt_te_surface_model_part_name" : "",
        "shed_wake_from_trailing_edge"     : false,
        "shed_wake_length"                 : 12.5,
        "shed_wake_element_size"           : 0.2,
        "shed_wake_grow_factor"            : 1.05,
        "shed_wake_projection_root_edge"   : 0.0,
        "echo_level"                       : 0
    })" );
    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    mWakeNormal                  = ThisParameters["wake_normal"].GetVector();
    mWakeSTLFileName             = ThisParameters["wake_stl_file_name"].GetString();
    mWakedrTraslation            = ThisParameters["wake_translation_direction"].GetVector();
    mTolerance                   = ThisParameters["tolerance"].GetDouble();
    mVisualizeWakeVTK            = ThisParameters["visualize_wake_vtk"].GetBool();
    mUpperSurfaceModelPartName   = ThisParameters["upper_surface_model_part_name"].GetString();
    mLowerSurfaceModelPartName   = ThisParameters["lower_surface_model_part_name"].GetString();
    mRootPointsModelPartName     = ThisParameters["root_points_model_part_name"].GetString();
    mTipPointsModelPartName      = ThisParameters["tip_points_model_part_name"].GetString();
    mBluntTESurfaceModelPartName = ThisParameters["blunt_te_surface_model_part_name"].GetString();
    mShedWakeFromTrailingEdge    = ThisParameters["shed_wake_from_trailing_edge"].GetBool();
    mShedWakeLength              = ThisParameters["shed_wake_length"].GetDouble();
    mShedWakeElementSize         = ThisParameters["shed_wake_element_size"].GetDouble();
    mShedGrowFactor              = ThisParameters["shed_wake_grow_factor"].GetDouble();
    mShedProjectionRootEdge      = ThisParameters["shed_wake_projection_root_edge"].GetDouble();
    mEchoLevel                   = ThisParameters["echo_level"].GetInt();
    
    KRATOS_ERROR_IF(mWakeNormal.size() != 3)
        << "The mWakeNormal should be a vector with 3 components!"
        << std::endl;
}

void Define3DWakeProcess::ExecuteInitialize()
{
    KRATOS_TRY;

    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    block_for_each(root_model_part.Nodes(), [&](Node& r_nodes)
    {
        r_nodes.SetValue(UPPER_SURFACE, false);
        r_nodes.SetValue(LOWER_SURFACE, false);
        r_nodes.SetValue(TRAILING_EDGE, false);
        r_nodes.SetValue(WAKE_DISTANCE, 0.0);
    });
    auto& r_elements = root_model_part.Elements();
    VariableUtils().SetNonHistoricalVariable(WAKE, 0, r_elements);
    
    // Save wake normal
    root_model_part.GetProcessInfo()[WAKE_NORMAL] = mWakeNormal;
    // Save free stream velocity direction as wake direction
    mWakeDirection = root_model_part.GetProcessInfo()[FREE_STREAM_VELOCITY_DIRECTION];
    // Compute span direction as the cross product: mWakeNormal x mWakeDirection
    MathUtils<double>::CrossProduct(mSpanDirection, mWakeNormal, mWakeDirection);
    
    InitializeTrailingEdgeSubModelpart();
    
    InitializeWakeSubModelpart();
    
    MarkTrailingEdgeAndBluntNodes();
    
    ComputeWingLowerSurfaceNormals();
    
    AddTrailingEdgeConditionsAndFindRootAndTipNodes();

    ComputeAndSaveLocalWakeNormal();
    
    Model aux_model;
    // Auxiliar model and model part for loading a STL mesh or building an automatic mesh
    ModelPart& StlWakeModelPart = aux_model.CreateModelPart("wake_model_part");
    
    KRATOS_WARNING_IF("Define3DWakeProcess", mWakeSTLFileName == "") << "Generating wake automatically. Set \'wake_stl_file_name\' instead" << std::endl;

    if(mShedWakeFromTrailingEdge || mWakeSTLFileName == ""){
        ShedWakeSurfaceFromTheTrailingEdge(StlWakeModelPart);
    }else
    {
        LoadSTL(StlWakeModelPart);
    }

    MoveWakeModelPart(StlWakeModelPart);
    
    VisualizeWake(StlWakeModelPart);
    
    MarkWakeElements(StlWakeModelPart);
    
    RecomputeNodalDistancesToWakeOrWingLowerSurface();
    
    MarkKuttaElements();
    
    SaveLocalWakeNormalInElements();
    
    AddWakeNodesToWakeModelPart();

    KRATOS_CATCH("");
}

// This function initializes the variables and removes all of the elements and
// nodes of the trailing edge submodelpart
void Define3DWakeProcess::InitializeTrailingEdgeSubModelpart() const
{
    KRATOS_TRY;
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    if(root_model_part.HasSubModelPart("trailing_edge_elements_model_part"))
    {
        // Clearing the variables and elements of the already existing
        // trailing_edge_sub_model_part
        ModelPart& trailing_edge_sub_model_part =
            root_model_part.GetSubModelPart("trailing_edge_elements_model_part");

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
        root_model_part.CreateSubModelPart("trailing_edge_elements_model_part");
    }
    KRATOS_CATCH("");
}

// This function initializes the variables and removes all of the elements and
// nodes of the wake submodelpart
void Define3DWakeProcess::InitializeWakeSubModelpart() const
{
    KRATOS_TRY;
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    if(root_model_part.HasSubModelPart("wake_elements_model_part"))
    {
        // Clearing the variables and elements of the already existing
        // wake_sub_model_part
        ModelPart& wake_sub_model_part =
            root_model_part.GetSubModelPart("wake_elements_model_part");

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
        root_model_part.CreateSubModelPart("wake_elements_model_part");
    }
    KRATOS_CATCH("");
}

// This function marks the trailing edge and blunt nodes
void Define3DWakeProcess::MarkTrailingEdgeAndBluntNodes()
{
    KRATOS_TRY;
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    
    for (auto& r_node : mrTrailingEdgeModelPart.Nodes()) {
        r_node.SetValue(TRAILING_EDGE, true);
    }

    if (mBluntTESurfaceModelPartName !="")
    {
        ModelPart& blunt_te_surface_model_part = root_model_part.GetSubModelPart(mBluntTESurfaceModelPartName);
        for (auto& r_node : blunt_te_surface_model_part.Nodes()) {
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
void Define3DWakeProcess::ComputeWingLowerSurfaceNormals() const
{
    KRATOS_TRY;
    // Mark upper surface
    if (mUpperSurfaceModelPartName != "")
    { 
        KRATOS_INFO_IF("Define3DWakeProcess", mEchoLevel > 0) << "...Using upper_surface_model_part" << std::endl;
        ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
        ModelPart& upper_sub_model_part = root_model_part.GetSubModelPart(mUpperSurfaceModelPartName);
        for (auto& r_cond : upper_sub_model_part.Conditions()) {
            auto& r_geometry = r_cond.GetGeometry();
            for (unsigned int j = 0; j < r_geometry.size(); j++) {
                r_geometry[j].SetLock();
                r_geometry[j].SetValue(UPPER_SURFACE, true);
                r_geometry[j].UnSetLock();
            }
        }
    }else
    {
        KRATOS_WARNING("Define3DWakeProcess") << "Marking Upper Surface automatically. Set \'upper_surface_model_part_name\' instead." << std::endl;
        KRATOS_INFO_IF("Define3DWakeProcess", mEchoLevel > 0) << "...Marking Upper Surface" << std::endl;
        for (auto& r_cond : mrBodyModelPart.Conditions()) {
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
    if (mLowerSurfaceModelPartName != "")
    {
        KRATOS_INFO_IF("Define3DWakeProcess", mEchoLevel > 0) << "...Using lower_surface_model_part" << std::endl;
        ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
        ModelPart& lower_sub_model_part = root_model_part.GetSubModelPart(mLowerSurfaceModelPartName);
        for (auto& r_cond : lower_sub_model_part.Conditions()) {
            auto& r_geometry = r_cond.GetGeometry();
            const auto& surface_normal = r_geometry.UnitNormal(0);
            const double projection = inner_prod(surface_normal, mWakeNormal);

            double switching_factor = 1.0;
            if(!(projection > 0.0)){
                switching_factor = -1.0;
            }
            
            for (unsigned int j = 0; j < r_geometry.size(); j++) {
                r_geometry[j].SetLock();
                r_geometry[j].SetValue(NORMAL, surface_normal*switching_factor);
                r_geometry[j].SetValue(LOWER_SURFACE, true);
                r_geometry[j].UnSetLock();
            }
        }
    }else
    {
        KRATOS_WARNING("Define3DWakeProcess") << "Marking Lower Surface automatically. Set \'lower_surface_model_part_name\' instead." << std::endl;
        KRATOS_INFO_IF("Define3DWakeProcess", mEchoLevel > 0) << "...Marking Lower Surface" << std::endl;
        for (auto& r_cond : mrBodyModelPart.Conditions()) {
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
void Define3DWakeProcess::ComputeAndSaveLocalWakeNormal() const
{
    KRATOS_TRY;    
    for (auto& r_cond : mrTrailingEdgeModelPart.Conditions()){
        auto& r_geometry = r_cond.GetGeometry();

        array_1d<double,3> coordinates1;
        array_1d<double,3> coordinates2 ;
        array_1d<double,3> local_span_direction;

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

    for (auto& r_node: mrTrailingEdgeModelPart.Nodes()){
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
void Define3DWakeProcess::ShedWakeSurfaceFromTheTrailingEdge(ModelPart& StlWakeModelPart) const
{
    KRATOS_TRY;
    KRATOS_INFO_IF("Define3DWakeProcess", mEchoLevel > 0) << "...Shedding wake from the trailing_edge_model_part" << std::endl;
    
    const Properties::Pointer pElemProp = mrBodyModelPart.pGetProperties(0);
    const double max_distance = mShedWakeLength;
    const double initial_element_size = mShedWakeElementSize;
    const double growth_factor = mShedGrowFactor;

    IndexType node_index = 0;
    IndexType element_index = 0;

    for (auto& r_cond : mrTrailingEdgeModelPart.Conditions()) {
        auto& r_geometry = r_cond.GetGeometry();
        
        array_1d<double,3> start_point;
        array_1d<double,3> end_point  ;

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
        double current_element_size = initial_element_size;
        double accumulated_distance = 0.0;

        array_1d<double,3> next_start = start_point + current_element_size * mWakeDirection;
        array_1d<double,3> next_end   = end_point   + current_element_size * mWakeDirection;

        CreateWakeSurfaceNodesAndElements(StlWakeModelPart, node_index, start_point, end_point, next_start, next_end, element_index, pElemProp);

        while (accumulated_distance + current_element_size <= max_distance) {
            accumulated_distance += current_element_size;
            current_element_size *= growth_factor;

            start_point = next_start;
            end_point   = next_end;
            
            next_start = start_point + current_element_size * mWakeDirection;
            next_end   = end_point   + current_element_size * mWakeDirection;
            
            CreateWakeSurfaceNodesAndElements(StlWakeModelPart, node_index, start_point, end_point, next_start, next_end, element_index, pElemProp);
        }
    }
    KRATOS_CATCH("");
}

void Define3DWakeProcess::CreateWakeSurfaceNodesAndElements(
    ModelPart& ModelPart,
    IndexType& rNode_index,
    const array_1d<double, 3>& rCoordinates1,
    const array_1d<double, 3>& rCoordinates2,
    const array_1d<double, 3>& rCoordinates3,
    const array_1d<double, 3>& rCoordinates4,
    IndexType& rElement_index,
    const Properties::Pointer pElemProp) const
{
    KRATOS_TRY;
    const std::array<ModelPart::IndexType, 4> nodes_ids = CreateWakeSurfaceNodes(ModelPart,
        rNode_index, rCoordinates1, rCoordinates2, rCoordinates3, rCoordinates4);
        
        const double normal_projection = ComputeFaceNormalProjectionToWakeNormal(
            rCoordinates1, rCoordinates2, rCoordinates3, rCoordinates4);
            
            CreateWakeSurfaceElements(ModelPart, normal_projection, rElement_index, nodes_ids, pElemProp);
    KRATOS_CATCH("");   
}

std::array<ModelPart::IndexType, 4> Define3DWakeProcess::CreateWakeSurfaceNodes(
    ModelPart& ModelPart,
    IndexType& rNode_index,
    const array_1d<double, 3>& rCoordinates1,
    const array_1d<double, 3>& rCoordinates2,
    const array_1d<double, 3>& rCoordinates3,
    const array_1d<double, 3>& rCoordinates4) const
{
    KRATOS_TRY;
    const auto& p_node1 = ModelPart.CreateNewNode(++rNode_index, rCoordinates1[0], rCoordinates1[1], rCoordinates1[2]);
    const auto& p_node2 = ModelPart.CreateNewNode(++rNode_index, rCoordinates2[0], rCoordinates2[1], rCoordinates2[2]);
    const auto& p_node3 = ModelPart.CreateNewNode(++rNode_index, rCoordinates3[0], rCoordinates3[1], rCoordinates3[2]);
    const auto& p_node4 = ModelPart.CreateNewNode(++rNode_index, rCoordinates4[0], rCoordinates4[1], rCoordinates4[2]);
    
    return {p_node1->Id(), p_node2->Id(), p_node3->Id(), p_node4->Id()};
    KRATOS_CATCH("");
}

double Define3DWakeProcess::ComputeFaceNormalProjectionToWakeNormal(
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

void Define3DWakeProcess::CreateWakeSurfaceElements(ModelPart& ModelPart,
                                                    const double normal_projection,
                                                    IndexType& rElement_index,
                                                    const std::array<ModelPart::IndexType, 4>& rNodes_ids,
                                                    const Properties::Pointer pElemProp) const
{
    KRATOS_TRY;
    if(normal_projection > 0.0){
        const std::vector<ModelPart::IndexType> elem_nodes_1{rNodes_ids[0], rNodes_ids[1], rNodes_ids[2]};
        const std::vector<ModelPart::IndexType> elem_nodes_2{rNodes_ids[1], rNodes_ids[3], rNodes_ids[2]};
        ModelPart.CreateNewElement("Element3D3N", ++rElement_index, elem_nodes_1, pElemProp);
        ModelPart.CreateNewElement("Element3D3N", ++rElement_index, elem_nodes_2, pElemProp);
    }
    else{
        const std::vector<ModelPart::IndexType> elem_nodes_1{rNodes_ids[0], rNodes_ids[2], rNodes_ids[1]};
        const std::vector<ModelPart::IndexType> elem_nodes_2{rNodes_ids[1], rNodes_ids[2], rNodes_ids[3]};
        ModelPart.CreateNewElement("Element3D3N", ++rElement_index, elem_nodes_1, pElemProp);
        ModelPart.CreateNewElement("Element3D3N", ++rElement_index, elem_nodes_2, pElemProp);
    }
    KRATOS_CATCH("");
}

// To load a stl mesh
void Define3DWakeProcess::LoadSTL(ModelPart& StlWakeModelPart) const
{
    KRATOS_TRY;
    KRATOS_INFO_IF("Define3DWakeProcess", mEchoLevel > 0) << "...Reading wake from stl file" << std::endl;
    
    // Process parameters
    Parameters reading_parameters = Parameters(R"(
    {
        "open_mode"       : "read",
        "new_entity_type" : "element"
    })" );
    
    StlIO stl_read (mWakeSTLFileName, reading_parameters);
    stl_read.ReadModelPart(StlWakeModelPart);
    KRATOS_CATCH("");
}

// To move the wake model part
void Define3DWakeProcess::MoveWakeModelPart(ModelPart& StlWakeModelPart) const
{   
    KRATOS_TRY;
    if (norm_2(mWakedrTraslation) > std::numeric_limits<double>::epsilon()) {
        KRATOS_INFO_IF("Define3DWakeProcess", mEchoLevel > 0) << "...Moving wake_model_part" << std::endl;
        
        Parameters moving_parameters = Parameters(R"(
        {
            "origin": []
        })" );
            
        moving_parameters["origin"].SetVector(mWakedrTraslation);
                
        MoveModelPartProcess MoveModelPartProcess(StlWakeModelPart, moving_parameters);
        MoveModelPartProcess.Execute();
    }
    KRATOS_CATCH("");
}

// To visualize the wake
void Define3DWakeProcess::VisualizeWake(ModelPart& StlWakeModelPart) const
{
    KRATOS_TRY;
    if (mVisualizeWakeVTK)
    {
        KRATOS_INFO_IF("Define3DWakeProcess", mEchoLevel > 0) << "...Saving vtk wake_model_part in 'wake_output' folder" << std::endl;
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
            
        VtkOutput vtk_oi(StlWakeModelPart, vtk_parameters);
        vtk_oi.PrintOutput();
    }
    KRATOS_CATCH("");
}

// This function checks which elements are cut by the wake and marks them as
// wake elements
void Define3DWakeProcess::MarkWakeElements(ModelPart& StlWakeModelPart) const
{
    KRATOS_TRY;
    KRATOS_INFO_IF("Define3DWakeProcess", mEchoLevel > 0) << "...Selecting wake elements" << std::endl;
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();

    BuiltinTimer timer;

    CalculateDiscontinuousDistanceToSkinProcess<3> distance_calculator(root_model_part, StlWakeModelPart);
    distance_calculator.Execute();

    KRATOS_INFO_IF("Define3DWakeProcess", mEchoLevel > 0) << "...Distance_calculator took " << timer.ElapsedSeconds() << " [sec]" << std::endl;

    // Variable to store element ids
    moodycamel::ConcurrentQueue<std::size_t> wake_elements_ordered_ids_concurrent_queue;
    moodycamel::ConcurrentQueue<std::size_t> trailing_edge_elements_ordered_ids_concurrent_queue;

    block_for_each(root_model_part.Elements(), [&](Element& rElement)
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
                if (std::abs(wake_elemental_distances[j]) < mTolerance)
                {
                    if (wake_elemental_distances[j] < 0.0)
                    {
                        wake_elemental_distances[j] = -mTolerance;
                    }
                    else
                    {
                        wake_elemental_distances[j] = mTolerance;
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

    KRATOS_INFO_IF("Define3DWakeProcess", mEchoLevel > 0) << "...Selecting wake elements finished" << std::endl;
    KRATOS_CATCH("");
}

// This function checks if the element is touching the trailing edge
void Define3DWakeProcess::CheckIfTrailingEdgeElement(
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
    if (mBluntTESurfaceModelPartName != "" && rElement.GetValue(TRAILING_EDGE))
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
void Define3DWakeProcess::AddTrailingEdgeAndWakeElements(std::vector<std::size_t>& rWakeElementsOrderedIds,
                                                         std::vector<std::size_t>& rTrailingEdgeElementsOrderedIds) const
{
    KRATOS_TRY;
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    
    std::sort(rWakeElementsOrderedIds.begin(),
              rWakeElementsOrderedIds.end());
    root_model_part.GetSubModelPart("wake_elements_model_part").AddElements(rWakeElementsOrderedIds);

    std::sort(rTrailingEdgeElementsOrderedIds.begin(),
              rTrailingEdgeElementsOrderedIds.end());
    ModelPart& trailing_edge_sub_model_part =
        root_model_part.GetSubModelPart("trailing_edge_elements_model_part");
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
void Define3DWakeProcess::RecomputeNodalDistancesToWakeOrWingLowerSurface() const
{
    KRATOS_TRY;
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    ModelPart& trailing_edge_sub_model_part =
        root_model_part.GetSubModelPart("trailing_edge_elements_model_part");

    block_for_each(trailing_edge_sub_model_part.Nodes(), [&](Node& r_node)
    {
        NodeType::Pointer p_closest_te_node = *mrTrailingEdgeModelPart.NodesBegin().base();
        FindClosestTrailingEdgeNode(p_closest_te_node, r_node);
        RecomputeDistance(p_closest_te_node, r_node);
    });
    KRATOS_CATCH("");
}

// This function finds the closest trailing edge node to the given point
void Define3DWakeProcess::FindClosestTrailingEdgeNode(NodeType::Pointer& pClosest_te_node,
                                                      const array_1d<double, 3>& rPoint) const
{
    KRATOS_TRY;
    double min_distance_to_te = std::numeric_limits<double>::max();
    for (auto& r_te_node : mrTrailingEdgeModelPart.Nodes())
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
void Define3DWakeProcess::RecomputeDistance(NodeType::Pointer& pClosest_te_node,
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
        if(std::abs(distance) < mTolerance){
            distance = - mTolerance;
        }
        rNode.SetValue(WAKE_DISTANCE, distance);
    }
    // For the nodes behind of the selected trailing edge node, the
    // distance is computed according to the wake normal
    else{
        const auto& local_wake_normal = pClosest_te_node->GetValue(WAKE_NORMAL);
        double distance = inner_prod(distance_vector, local_wake_normal);
        // Nodes slightly below and above the wake are given a positive distance (wake down)
        if(std::abs(distance) < mTolerance){
            distance = mTolerance;
        }
        rNode.SetValue(WAKE_DISTANCE, distance);
    }

    if (rNode.GetValue(TRAILING_EDGE))
    {
        rNode.SetValue(WAKE_DISTANCE, mTolerance);
    }

    KRATOS_CATCH("");
}

// This function selects the kutta elements. Kutta elements are touching the
// trailing edge from below
void Define3DWakeProcess::MarkKuttaElements() const
{
    KRATOS_TRY;
    KRATOS_INFO_IF("Define3DWakeProcess", mEchoLevel > 0) << "...Selecting kutta elements" << std::endl;
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    ModelPart& trailing_edge_sub_model_part =
    root_model_part.GetSubModelPart("trailing_edge_elements_model_part");

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
    ModelPart& wake_sub_model_part = root_model_part.GetSubModelPart("wake_elements_model_part");
    wake_sub_model_part.RemoveElements(TO_ERASE);
    KRATOS_INFO_IF("Define3DWakeProcess", mEchoLevel > 0) << "...Selecting kutta elements finished" << std::endl;
    KRATOS_CATCH("");
}

// This function returns the number of trailing edge nodes
unsigned int Define3DWakeProcess::CountNumberOfTrailindEdgeNodesInElement(
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
void Define3DWakeProcess::CountNumberOfPositiveAndNegativeDistances(
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
void Define3DWakeProcess::SelectElementType(Element& rElement,
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
void Define3DWakeProcess::SaveLocalWakeNormalInElements() const
{
    KRATOS_TRY;
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    ModelPart& wake_sub_model_part = root_model_part.GetSubModelPart("wake_elements_model_part");

    for (auto& r_elem : wake_sub_model_part.Elements()){
        NodeType::Pointer p_closest_te_node = *mrTrailingEdgeModelPart.NodesBegin().base();
        FindClosestTrailingEdgeNode(p_closest_te_node, r_elem.GetGeometry().Center());

        r_elem.GetValue(WAKE_NORMAL) = p_closest_te_node->GetValue(WAKE_NORMAL);
    }
    KRATOS_CATCH("");
}

void Define3DWakeProcess::AddWakeNodesToWakeModelPart() const
{
    KRATOS_TRY;
    ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
    ModelPart& wake_sub_model_part =
            root_model_part.GetSubModelPart("wake_elements_model_part");

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
void Define3DWakeProcess::AddTrailingEdgeConditionsAndFindRootAndTipNodes() const
{
    KRATOS_TRY;
    bool assign_conditions = (mrTrailingEdgeModelPart.NumberOfConditions() < 1);
    if (assign_conditions ||
        mTipPointsModelPartName == "" ||
        mRootPointsModelPartName == "")
    {
        KRATOS_WARNING_IF("Define3DWakeProcess", assign_conditions) << "Assigning conditions automatically. Check the trailing_edge_model_part" << std::endl;
        KRATOS_WARNING_IF("Define3DWakeProcess", mTipPointsModelPartName == "")  << "Detecting tip nodes automatically. Set \'tip_points_model_part_name\' instead" << std::endl;
        KRATOS_WARNING_IF("Define3DWakeProcess", mRootPointsModelPartName == "") << "Detecting root nodes automatically. Set \'root_points_model_part_name\' instead" << std::endl;

        const double constant_threshold = 3.0;
        const double distance_threshold = constant_threshold * mShedWakeElementSize;
        ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
        int reference_id = root_model_part.NumberOfConditions();
        const Properties::Pointer pCondProp = mrTrailingEdgeModelPart.pGetProperties(0);
    
        std::vector<Node*> all_nodes;
        all_nodes.reserve(mrTrailingEdgeModelPart.NumberOfNodes());
        for (auto& r_node : mrTrailingEdgeModelPart.Nodes()) {
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
    
            if (mTipPointsModelPartName == "" ||
                mRootPointsModelPartName == "")
            {
                auto [min_it, max_it] = std::minmax_element(
                    group.begin(), group.end(),
                    [this](const Node* a, const Node* b) {
                        return inner_prod(a->Coordinates(), mSpanDirection) <
                               inner_prod(b->Coordinates(), mSpanDirection);
                    }
                );
                
                if (mTipPointsModelPartName == ""    &&
                    !((*max_it)->GetValue(WING_TIP)) &&
                    !((*max_it)->GetValue(WING_ROOT)))
                {
                    (*max_it)->SetValue(WING_TIP, true);
                    tip_nodes += 1;
                }
                
                if (mRootPointsModelPartName == ""   &&
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
                    mrTrailingEdgeModelPart.CreateNewCondition("LineCondition3D2N", reference_id, cond, pCondProp);
                }
            }
        }
    
        KRATOS_INFO_IF("Define3DWakeProcess", mEchoLevel > 0) << "...Detected " << trailing_edge_groups.size() << " trailing edge groups" << std::endl;
        KRATOS_INFO_IF("Define3DWakeProcess", (mEchoLevel > 0 && (mTipPointsModelPartName == ""))) << "...Detected " << tip_nodes << " tip nodes" << std::endl;
        KRATOS_INFO_IF("Define3DWakeProcess", (mEchoLevel > 0 && (mRootPointsModelPartName == ""))) << "...Detected " << root_nodes << " root nodes" << std::endl;
        KRATOS_INFO_IF("Define3DWakeProcess", (mEchoLevel > 0 && assign_conditions)) << "...Assigned " << mrTrailingEdgeModelPart.NumberOfConditions() << " conditions" << std::endl;
    }

    if (mTipPointsModelPartName != "")
    {
        ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
        ModelPart& tip_sub_model_part = root_model_part.GetSubModelPart(mTipPointsModelPartName);
        int tip_nodes = 0;
        for (auto& r_node : tip_sub_model_part.Nodes()) {
            r_node.SetValue(WING_TIP, true);
            tip_nodes += 1;
        }
        KRATOS_INFO_IF("Define3DWakeProcess", mEchoLevel > 0) << "...Detected " << tip_nodes << " tip nodes using " << mTipPointsModelPartName << " model part" << std::endl;
    }
    
    if (mRootPointsModelPartName != "")
    {
        int root_nodes = 0;
        ModelPart& root_model_part = mrBodyModelPart.GetRootModelPart();
        ModelPart& root_sub_model_part = root_model_part.GetSubModelPart(mRootPointsModelPartName);
        for (auto& r_node : root_sub_model_part.Nodes()) {
            r_node.SetValue(WING_ROOT, true);
            root_nodes += 1;
        }
        KRATOS_INFO_IF("Define3DWakeProcess", mEchoLevel > 0) << "...Detected " << root_nodes << " root nodes using " << mRootPointsModelPartName << " model part" << std::endl;
    }
       
    KRATOS_CATCH("");
}
} // namespace Kratos.
