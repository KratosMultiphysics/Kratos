#include "voxel_mesher_coloring.h"
#include "modeler/voxel_mesh_generator_modeler.h"

#include "includes/checks.h"

namespace Kratos {

VoxelMesherColoring::VoxelMesherColoring(VoxelMeshGeneratorModeler& rModeler, Parameters ColoringParameters):
    mrModeler(rModeler),
    mParameters(ColoringParameters)
{}


void VoxelMesherColoring::ValidateParameters()
{
    mParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
}


Parameters VoxelMesherColoring::GetParameters() const
{
    return mParameters;
}


Parameters VoxelMesherColoring::GetDefaultParameters() const
{
    return Parameters(R"({
        "type" : "Undefined_coloring_type",
        "model_part_name": "Undefined",
        "color": -1,
        "cell_color": -1,
        "input_entities": "elements",
        "outside_color" : 0
    })");
}


void VoxelMesherColoring::CheckGeometryType(const GeometryData::KratosGeometryType &rType) const {
    KRATOS_ERROR_IF_NOT(rType == GeometryData::KratosGeometryType::Kratos_Triangle3D3) << " Input entities must be of type Triangle3D3" << std::endl;
}


ModelPart& VoxelMesherColoring::GetModelPart(const std::string& rName) const {
    return mrModeler.mpModel->GetModelPart(rName);
}


const ModelPart& VoxelMesherColoring::GetInputModelPart() const {
    return *(mrModeler.mpInputModelPart);
}

VoxelMesherColoring::CartesianMeshColors& VoxelMesherColoring::GetMeshColors() const {
    return mrModeler.mColors;
}


VoxelMesherColoring::SkinIntersection& VoxelMesherColoring::GetIntersections(int Color) const
{
    return mrModeler.GetIntersections(Color);
}


const std::vector<double>& VoxelMesherColoring::GetKeyPlanes(std::size_t Direction) const {
    return mrModeler.mKeyPlanes[Direction];
}

}