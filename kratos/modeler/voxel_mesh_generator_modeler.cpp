//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes

// External includes

// Project includes
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "geometries/bounding_box.h"
#include "geometries/point.h"
#include "utilities/parallel_utilities.h"
#include "utilities/timer.h"
#include "spatial_containers/geometrical_objects_bins.h"
#include "utilities/model_part_utils.h"

#include "voxel_mesh_generator_modeler.h"
#include "coloring/voxel_mesher_coloring_factory.h"
#include "key_plane_generation/key_plane_generation_factory.h"
#include "entity_generation/voxel_mesher_entity_generation_factory.h"
#include "operation/voxel_mesher_operation_factory.h"

namespace Kratos
{

/// Default constructor.
VoxelMeshGeneratorModeler::VoxelMeshGeneratorModeler() : Modeler()
{
}

/// Constructor.
VoxelMeshGeneratorModeler::VoxelMeshGeneratorModeler(
    Model& rModel, Parameters rParameters)
    : Modeler(rModel, rParameters )
    , mpModel(&rModel)
{
    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());
    rParameters["key_plane_generator"].ValidateAndAssignDefaults(GetDefaultKeyPlaneGeneratorParameters());
}

///
const Parameters VoxelMeshGeneratorModeler::GetDefaultParameters() const
{
    return Parameters(R"(
    {
        "key_plane_generator": {},
        "entities_generator_list": [],
        "coloring_settings_list": [],
        "model_part_operations": [],
        "output_filename" : "",
        "mdpa_file_name" : "",
        "output_model_part_name" : "",
        "input_model_part_name" : "",
        "default_outside_color" : 1,
        "voxel size" : 0.001,
        "output_files" : []
    }  )");
}

///
const Parameters VoxelMeshGeneratorModeler::GetDefaultKeyPlaneGeneratorParameters() const
{
    return Parameters(R"(
    {
        "type": "bounding_box",
        "Parameters":{}
    }  )");
}

///
void VoxelMeshGeneratorModeler::SetupModelPart()
{
    ModelPart& main_model_part = ReadModelParts();

    KRATOS_INFO("Modeler") << "Generating Key planes..." << std::endl;
    KeyPlaneGenerationFactory factory;
    auto p_key_plane_generator = factory.Create(*this, mParameters["key_plane_generator"]);
    p_key_plane_generator->ValidateParameters();
    p_key_plane_generator->Generate();

    KRATOS_INFO("Modeler") << "Key Planes generated" << std::endl;

    KRATOS_INFO("Modeler") << "Preparing Internal Data Structure" << std::endl;
    PreparingTheInternalDataStructure();
    KRATOS_INFO("Modeler") << "Internal Data Structure prepared" << std::endl;

    if(mParameters.Has("coloring_settings_list")){
        ApplyColoring(mParameters["coloring_settings_list"], mParameters["default_outside_color"].GetInt());
    }
    if(mParameters.Has("entities_generator_list")){
        GenerateEntities(main_model_part, mParameters["entities_generator_list"]);
    }
    if(mParameters.Has("model_part_operations")){
        ApplyOperations(mParameters["model_part_operations"]);
    }

    for (const auto& single_output_settings : mParameters["output_files"]){
        const std::string type = single_output_settings["type"].GetString();
        if( type == "vtr"){
            mColors.WriteParaViewVTR("Coloring.vtr");
        }
        else if (type == "esri_ascii") {
            mColors.PrintAllVoxelColorsToAsciiFile(single_output_settings["file_name"].GetString());
        }
    }
}

ModelPart& VoxelMeshGeneratorModeler::ReadModelParts(){
    KRATOS_ERROR_IF_NOT( mParameters.Has("input_model_part_name") )
        << "Missing \"input_model_part_name\" in VoxelMeshGeneratorModeler Parameters." << std::endl;

    mpInputModelPart = &CreateAndGetModelPart(mParameters["input_model_part_name"].GetString());

    KRATOS_ERROR_IF_NOT(mParameters.Has("mdpa_file_name")) << "mdpa_file_name not defined" << std::endl;

    const std::string data_file_name =  mParameters["mdpa_file_name"].GetString();

    KRATOS_INFO_IF("::[VoxelMeshGeneratorModeler]::", mEchoLevel > 0) << "Importing Cad Model from: " << data_file_name << std::endl;

    // Load the mdpa
    if(data_file_name!=""){
        ModelPartIO(data_file_name).ReadModelPart(*mpInputModelPart);
    }

    // Create the target model
    return CreateAndGetModelPart(mParameters["output_model_part_name"].GetString());
}


void VoxelMeshGeneratorModeler::PreparingTheInternalDataStructure(){
    mColors.SetCoordinates(mKeyPlanes[0], mKeyPlanes[1], mKeyPlanes[2]);
    mMeshingData.SetNumberOfDivisions(mKeyPlanes[0].size(), mKeyPlanes[1].size(), mKeyPlanes[2].size());

    array_1d<double,3> margins = ZeroVector(3);
    for(int i = 0 ; i < 3 ; i++) {
        margins[i] = (mKeyPlanes[i].back() - mKeyPlanes[i].front()) * 1.0e-2;
    }

    const double margin = *std::max_element(margins.begin(), margins.end());
    mColors.ExtendBoundingBox(mpInputModelPart->Nodes(), margin);
}

void VoxelMeshGeneratorModeler::ApplyColoring(Parameters ColoringParameters, int OutsideColor){
    Timer::Start("MeshColoring");

    VoxelMesherColoringFactory factory("coloring strategy");
    mColors.SetAllColors(OutsideColor);

    for(auto parameters : ColoringParameters){

        if (!parameters.Has("outside_color")) {
            parameters.AddInt("outside_color", OutsideColor);
        }

        VoxelMesherColoring::Pointer p_coloring = factory.Create(*this, parameters);

        p_coloring->ValidateParameters();

        std::string model_part_name = parameters["model_part_name"].GetString();
        std::string type = parameters["type"].GetString();
        KRATOS_INFO("Modeler") << "Applying color to " << type << " for " << model_part_name  << " model part" << std::endl;

        p_coloring->Apply();
    }
    Timer::Stop("MeshColoring");
}

void VoxelMeshGeneratorModeler::SetStartIds(ModelPart& rTheVolumeModelPart)
{
    ModelPart& r_root_model_part = rTheVolumeModelPart.GetRootModelPart();
    std::size_t start_node_proposal = r_root_model_part.NodesArray().empty() ? 1 : r_root_model_part.NodesArray().back()->Id() + 1;
    mStartNodeId = std::max(start_node_proposal, mStartNodeId);

    std::size_t start_element_proposal = r_root_model_part.ElementsArray().empty() ? 1 : r_root_model_part.ElementsArray().back()->Id() + 1;
    mStartElementId = std::max(start_element_proposal, mStartElementId);

    std::size_t start_condition_proposal = r_root_model_part.ConditionsArray().empty()? 1 :r_root_model_part.ConditionsArray().back()->Id() + 1;
    mStartConditionId = std::max(start_condition_proposal, mStartConditionId);
}

void VoxelMeshGeneratorModeler::GenerateEntities(
    ModelPart& rTheVolumeModelPart,
    Parameters EntityGeneratorParameters
    )
{
    Timer::Start("EntityGeneration");
    VoxelMesherEntityGenerationFactory factory("entity generation strategy");

    for(auto parameters : EntityGeneratorParameters){

        auto p_generator = factory.Create(*this, parameters);
        p_generator->ValidateParameters();

        std::string model_part_name = parameters["model_part_name"].GetString();
        std::string type = parameters["type"].GetString();
        KRATOS_INFO("Modeler") << "Generate " << type << " in " << model_part_name  << " model part" << std::endl;

        p_generator->Generate();
    }
    Timer::Stop("EntityGeneration");
}

ModelPart& VoxelMeshGeneratorModeler::CreateAndGetModelPart(std::string const& FullName)
{
    return mpModel->HasModelPart(FullName) ? mpModel->GetModelPart(FullName) : mpModel->CreateModelPart(FullName);
}


Node::Pointer VoxelMeshGeneratorModeler::GenerateOrRetrieveNode(
    ModelPart& rTheVolumeModelPart,
    std::vector<ModelPart::NodeType::Pointer>& rThisNodes,
    const std::size_t I,
    const std::size_t J,
    const std::size_t K
    )
{
    auto& nodal_data = mMeshingData.GetNodalData(I, J, K);
    if(nodal_data.IsCreated()){
        rThisNodes.push_back(nodal_data.pGetNode());
        return nodal_data.pGetNode();
    }

    Point global_coordinates = mColors.GetPoint(I,J,K);
    Node::Pointer node = this->GenerateNode(rTheVolumeModelPart, global_coordinates);
    nodal_data.pSetNode(node);
    rThisNodes.push_back(node);
    return node;
}

Node::Pointer VoxelMeshGeneratorModeler::GenerateOrRetrieveQuadraticNode(
    ModelPart& rTheVolumeModelPart,
    std::vector<ModelPart::NodeType::Pointer>& rThisNodes,
    const std::size_t I,
    const std::size_t J,
    const std::size_t K
    )
{
    auto& nodal_data = mQuadraticNodeData.GetNodalData(I, J, K);
    if(nodal_data.IsCreated()){
        rThisNodes.push_back(nodal_data.pGetNode());
        return nodal_data.pGetNode();
    }

    Point global_coordinates = mColors.GetPoint(I/2.0,J/2.0,K/2.0);
    Node::Pointer node = this->GenerateNode(rTheVolumeModelPart, global_coordinates);
    nodal_data.pSetNode(node);
    rThisNodes.push_back(node);
    return node;
}

Node::Pointer VoxelMeshGeneratorModeler::GenerateNode(ModelPart& rTheVolumeModelPart, const Point& rCoordinates)
{
    Node::Pointer node = Kratos::make_intrusive< Node >( mStartNodeId++, rCoordinates[0], rCoordinates[1], rCoordinates[2]);
    // Giving model part's variables list to the node
    node->SetSolutionStepVariablesList(rTheVolumeModelPart.pGetNodalSolutionStepVariablesList());

    //set buffer size
    node->SetBufferSize(rTheVolumeModelPart.GetBufferSize());

    return node;
}

void VoxelMeshGeneratorModeler::ApplyOperations(Parameters ThisParameters)
{
    VoxelMesherOperationFactory factory("operation");
    int step = 0;
    for(auto parameters : ThisParameters){
        auto p_operation = factory.Create(*this, parameters);
        p_operation->ValidateParameters();

        std::stringstream stream;
        stream << "Operation " << ++step << " " << parameters["type"].GetString();
        const std::string label = stream.str();
        Timer::Start(label);

        p_operation->Execute();

        Timer::Stop(label);
    }
}


Internals::SkinIntersection& VoxelMeshGeneratorModeler::GetIntersections(int Color)
{
    auto iSkinIntersection = mSkinIntersections.find(Color);
    if (iSkinIntersection != mSkinIntersections.end()) {
        return iSkinIntersection->second;
    }
    auto emplace_out = mSkinIntersections.emplace(
        std::pair<int, Internals::SkinIntersection>(Color, mMeshingData.GetNumberOfDivisions()));
    return emplace_out.first->second;
}


const Internals::SkinIntersection& VoxelMeshGeneratorModeler::GetIntersections(int Color) const
{
    auto iSkinIntersection = mSkinIntersections.find(Color);
    if (iSkinIntersection != mSkinIntersections.end()) {
        return iSkinIntersection->second;
    }
    KRATOS_ERROR << "No intersections defined for color " << Color << std::endl;
}


bool VoxelMeshGeneratorModeler::IntersectionsGenerated(int Color) const
{
    return mSkinIntersections.find(Color) != mSkinIntersections.end();
}

} // namespace Kratos
