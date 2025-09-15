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
#include "utilities/model_part_utils.h"

#include "voxel_mesher_entity_generation.h"

#include "modeler/voxel_mesh_generator_modeler.h"

namespace Kratos {

VoxelMesherEntityGeneration::VoxelMesherEntityGeneration(VoxelMeshGeneratorModeler& rModeler, Parameters GenerationParameters):
    mrModeler(rModeler),
    mParameters(GenerationParameters)
{}


void VoxelMesherEntityGeneration::ValidateParameters()
{
    mParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
}


void VoxelMesherEntityGeneration::Generate() {
    std::string model_part_name = mParameters["model_part_name"].GetString();
    ModelPart& r_model_part = CreateAndGetModelPart(model_part_name);
    SetStartIds(r_model_part);

    Generate(r_model_part, mParameters);
}


Parameters VoxelMesherEntityGeneration::GetParameters() const
{
    return mParameters;
}


Parameters VoxelMesherEntityGeneration::GetDefaultParameters() const
{
    return Parameters(R"({
        "type" : "Undefined_entity_generation_type",
        "model_part_name": "Undefined",
        "color": -1,
        "properties_id": 1,
        "generated_entity": "Undefined"
    })");
}


ModelPart& VoxelMesherEntityGeneration::GetModelPart(const std::string& rName) const {
    return mrModeler.mpModel->GetModelPart(rName);
}


ModelPart& VoxelMesherEntityGeneration::CreateAndGetModelPart(std::string const& rFullName) const
{
    return mrModeler.CreateAndGetModelPart(rFullName);
}


Properties::Pointer VoxelMesherEntityGeneration::CreateAndGetProperty(ModelPart& rModelPart, std::size_t PropertyId) const
{
    if(!rModelPart.RecursivelyHasProperties(PropertyId)){
        rModelPart.CreateNewProperties(PropertyId);
    }
    return rModelPart.pGetProperties(PropertyId);
}


const VoxelMesherEntityGeneration::CartesianMeshColors& VoxelMesherEntityGeneration::GetMeshColors() const
{
    return mrModeler.mColors;
}


const array_1d<std::size_t, 3>& VoxelMesherEntityGeneration::GetNumberOfDivisions() const {
    return mrModeler.mMeshingData.GetNumberOfDivisions();
}


const VoxelMesherEntityGeneration::SkinIntersection& VoxelMesherEntityGeneration::GetIntersections(int Color) const {
    return mrModeler.GetIntersections(Color);
}


VoxelMesherEntityGeneration::SkinIntersection& VoxelMesherEntityGeneration::GetIntersections(int Color) {
    return mrModeler.GetIntersections(Color);
}


bool VoxelMesherEntityGeneration::IntersectionsGenerated(int Color) const {
    return mrModeler.IntersectionsGenerated(Color);
}


void VoxelMesherEntityGeneration::AddNodesToModelPart(ModelPart& rModelPart, ModelPart::NodesContainerType& rNewNodes) const
{
    rNewNodes.Unique();
    ModelPartUtils::AddNodesFromOrderedContainer(rModelPart, rNewNodes.begin(), rNewNodes.end());
}


void VoxelMesherEntityGeneration::SetStartIds(ModelPart& rModelPart)
{
    mrModeler.SetStartIds(rModelPart);
}


std::size_t VoxelMesherEntityGeneration::GetStartElementId() const
{
    return mrModeler.mStartElementId;
}


std::size_t VoxelMesherEntityGeneration::GetStartConditionId() const
{
    return mrModeler.mStartConditionId;
}

array_1d<std::size_t, 3> VoxelMesherEntityGeneration::GetNumberOfCells() const
{
    array_1d<std::size_t, 3> number_of_cells;
    for(int i = 0 ; i < 3 ; i++){
        number_of_cells[i] = mrModeler.mKeyPlanes[i].size() - 1;
    }
    return number_of_cells;
}


Node::Pointer VoxelMesherEntityGeneration::GenerateOrRetrieveNode(
    ModelPart& rTheVolumeModelPart,
    ModelPart::NodesContainerType& rThisNodes,
    const std::size_t I,
    const std::size_t J,
    const std::size_t K)
{
    return mrModeler.GenerateOrRetrieveNode(rTheVolumeModelPart, rThisNodes, I, J, K);
}


Node::Pointer VoxelMesherEntityGeneration::GenerateOrRetrieveQuadraticNode(
    ModelPart& rTheVolumeModelPart,
    ModelPart::NodesContainerType& rThisNodes,
    const std::size_t I,
    const std::size_t J,
    const std::size_t K)
{
    return mrModeler.GenerateOrRetrieveQuadraticNode(rTheVolumeModelPart, rThisNodes, I, J, K);
}


Node::Pointer VoxelMesherEntityGeneration::GenerateNode(ModelPart& rModelPart, const Point& rCoordintates)
{
    return mrModeler.GenerateNode(rModelPart, rCoordintates);
}


void VoxelMesherEntityGeneration::GetLinearCellNodes(
    Element::NodesArrayType& rCellNodes,
    ModelPart& rTheVolumeModelPart,
    ModelPart::NodesContainerType& rThisNodes,
    const std::size_t I,
    const std::size_t J,
    const std::size_t K)
{
    rCellNodes(0) = GenerateOrRetrieveNode(rTheVolumeModelPart, rThisNodes, I  , J  , K);
    rCellNodes(1) = GenerateOrRetrieveNode(rTheVolumeModelPart, rThisNodes, I+1, J  , K);
    rCellNodes(2) = GenerateOrRetrieveNode(rTheVolumeModelPart, rThisNodes, I+1, J+1, K);
    rCellNodes(3) = GenerateOrRetrieveNode(rTheVolumeModelPart, rThisNodes, I  , J+1, K);
    rCellNodes(4) = GenerateOrRetrieveNode(rTheVolumeModelPart, rThisNodes, I  , J  , K+1);
    rCellNodes(5) = GenerateOrRetrieveNode(rTheVolumeModelPart, rThisNodes, I+1, J  , K+1);
    rCellNodes(6) = GenerateOrRetrieveNode(rTheVolumeModelPart, rThisNodes, I+1, J+1, K+1);
    rCellNodes(7) = GenerateOrRetrieveNode(rTheVolumeModelPart, rThisNodes, I  , J+1, K+1);
}


void VoxelMesherEntityGeneration::GetQuadraticCellNodes(
    Element::NodesArrayType& rCellNodes,
    ModelPart& rTheVolumeModelPart,
    ModelPart::NodesContainerType& rThisNodes,
    const std::size_t I,
    const std::size_t J,
    const std::size_t K)
{
    /*      3----10----2
     *      |\         |\
     *      |15    23  | 14
     *      11  \ 20   9  \
     *      |   7----18+---6
     *      |24 |  26  | 22|
     *      0---+-8----1   |
     *       \ 19    25 \  17
     *       12 |  21    13|
     *         \|         \|
     *          4----16----5
     */
    GetLinearCellNodes(rCellNodes, rTheVolumeModelPart, rThisNodes, I, J, K);

    const std::size_t i = 2*I;
    const std::size_t j = 2*J;
    const std::size_t k = 2*K;
    rCellNodes( 8) = GenerateOrRetrieveQuadraticNode(rTheVolumeModelPart, rThisNodes, i+1, j  , k);
    rCellNodes( 9) = GenerateOrRetrieveQuadraticNode(rTheVolumeModelPart, rThisNodes, i+2, j+1, k);
    rCellNodes(10) = GenerateOrRetrieveQuadraticNode(rTheVolumeModelPart, rThisNodes, i+1, j+2, k);
    rCellNodes(11) = GenerateOrRetrieveQuadraticNode(rTheVolumeModelPart, rThisNodes, i  , j+1, k);
    rCellNodes(12) = GenerateOrRetrieveQuadraticNode(rTheVolumeModelPart, rThisNodes, i  , j  , k+1);
    rCellNodes(13) = GenerateOrRetrieveQuadraticNode(rTheVolumeModelPart, rThisNodes, i+2, j  , k+1);
    rCellNodes(14) = GenerateOrRetrieveQuadraticNode(rTheVolumeModelPart, rThisNodes, i+2, j+2, k+1);
    rCellNodes(15) = GenerateOrRetrieveQuadraticNode(rTheVolumeModelPart, rThisNodes, i  , j+2, k+1);
    rCellNodes(16) = GenerateOrRetrieveQuadraticNode(rTheVolumeModelPart, rThisNodes, i+1, j  , k+2);
    rCellNodes(17) = GenerateOrRetrieveQuadraticNode(rTheVolumeModelPart, rThisNodes, i+2, j+1, k+2);
    rCellNodes(18) = GenerateOrRetrieveQuadraticNode(rTheVolumeModelPart, rThisNodes, i+1, j+2, k+2);
    rCellNodes(19) = GenerateOrRetrieveQuadraticNode(rTheVolumeModelPart, rThisNodes, i  , j+1, k+2);
    rCellNodes(20) = GenerateOrRetrieveQuadraticNode(rTheVolumeModelPart, rThisNodes, i+1, j+1, k);
    rCellNodes(21) = GenerateOrRetrieveQuadraticNode(rTheVolumeModelPart, rThisNodes, i+1, j  , k+1);
    rCellNodes(22) = GenerateOrRetrieveQuadraticNode(rTheVolumeModelPart, rThisNodes, i+2, j+1, k+1);
    rCellNodes(23) = GenerateOrRetrieveQuadraticNode(rTheVolumeModelPart, rThisNodes, i+1, j+2, k+1);
    rCellNodes(24) = GenerateOrRetrieveQuadraticNode(rTheVolumeModelPart, rThisNodes, i  , j+1, k+1);
    rCellNodes(25) = GenerateOrRetrieveQuadraticNode(rTheVolumeModelPart, rThisNodes, i+1, j+1, k+2);
    rCellNodes(26) = GenerateOrRetrieveQuadraticNode(rTheVolumeModelPart, rThisNodes, i+1, j+1, k+1);
}


void VoxelMesherEntityGeneration::InitializeQuadraticData() {
    if (!mrModeler.mQuadraticNodeData.IsInitialized()) {
        const auto& r_divisions = mrModeler.mMeshingData.GetNumberOfDivisions();
        mrModeler.mQuadraticNodeData.SetNumberOfDivisions(2*r_divisions[0], 2*r_divisions[1], 2*r_divisions[2]);
    }
}

}