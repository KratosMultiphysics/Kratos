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
#include "generate_stl_intersection_with_cells.h"

namespace Kratos {

Parameters GenerateStlIntersectionWithCells::GetDefaultParameters() const
{
    return Parameters(R"({
        "type" : "generate_stl_intersection_with_cells",
        "model_part_name": "Undefined",
        "color": -1,
        "properties_id": 1,
        "generated_entity": "Element3D3N",
        "input_entities": "elements",
        "input_model_part_name": ""
    })");
}

void GenerateStlIntersectionWithCells::Generate(ModelPart& rModelPart, Parameters parameters)
{
    const CartesianMeshColors& r_colors = GetMeshColors();
    std::size_t element_id = GetStartElementId();

    const int inside_color = parameters["color"].GetInt();
    std::size_t properties_id = parameters["properties_id"].GetInt();
    Properties::Pointer p_properties = CreateAndGetProperty(rModelPart, properties_id);

    const std::string input_model_part_name = parameters["input_model_part_name"].GetString();
    if (input_model_part_name != "") {
        const std::string input_entities = parameters["input_entities"].GetString();
        GenerateIntersection(GetModelPart(input_model_part_name), input_entities, inside_color);
    }

    KRATOS_ERROR_IF_NOT(IntersectionsGenerated(inside_color)) <<
        "Intersections for color " << inside_color << " must be computed in advance to generate STL." << std::endl;

    auto& r_prototype_element = KratosComponents<Element>::Get(
        parameters["generated_entity"].GetString());
    const SkinIntersection& r_intersections = GetIntersections(inside_color);

    ModelPart::NodesContainerType new_nodes;
    ModelPart::ElementsContainerType new_elements;

    array_1d<std::size_t, 3> number_of_cells = GetNumberOfCells();

    for (std::size_t k = 0; k < number_of_cells[2]; k++) {
        for (std::size_t j = 0; j < number_of_cells[1]; j++) {
            for (std::size_t i = 0; i < number_of_cells[0]; i++) {
                if(std::lround(r_colors.GetElementalColor(i,j,k)) == inside_color){
                    const auto& r_geometries = r_intersections.CutGeometriesInCell(i, j, k);
                    // Note: this is a STL-like mesh for output, the nodes are not unique
                    auto cut_nodes = GenerateCutGeometryNodes(rModelPart, r_geometries);

                    for (auto i_node = cut_nodes.ptr_begin(); i_node != cut_nodes.ptr_end(); i_node += 3) {
                        auto geom = Kratos::make_shared<Triangle3D3<Node>>(*i_node, *(i_node+1), *(i_node+2));
                        new_elements.push_back(
                            r_prototype_element.Create(element_id++, geom, p_properties));

                    }

                    for (auto iter = cut_nodes.ptr_begin(); iter != cut_nodes.ptr_end(); ++iter) {
                        new_nodes.push_back(*iter);
                    }
                }
            }
        }
    }

    AddNodesToModelPart(rModelPart, new_nodes);
    rModelPart.AddElements(new_elements.begin(), new_elements.end());
}

PointerVector<Node> GenerateStlIntersectionWithCells::GenerateCutGeometryNodes(
    ModelPart& rTheCutModelPart,
    const std::vector<typename WorkGeometryType::Pointer>& rCutGeometries)
{
    PointerVector<Node> nodes;
    nodes.reserve(3*rCutGeometries.size());
    for (auto i_geom = rCutGeometries.begin(); i_geom < rCutGeometries.end(); ++i_geom) {
        for (auto& r_point: **i_geom) {
            nodes.push_back(GenerateNode(rTheCutModelPart, r_point));
        }
    }

    return nodes;
}


void GenerateStlIntersectionWithCells::CheckGeometryType(const GeometryData::KratosGeometryType &rType) const {
    KRATOS_ERROR_IF_NOT(rType == GeometryData::KratosGeometryType::Kratos_Triangle3D3) <<
        " Input entities must be of type Triangle3D3" << std::endl;
}


void GenerateStlIntersectionWithCells::GenerateIntersection(const ModelPart& rReferenceModelPart, const std::string& rInputEntities, int Color) {
    const CartesianMeshColors& r_colors = GetMeshColors();
    SkinIntersection& r_skin_intersection = GetIntersections(Color);

    if (rInputEntities == "elements") {
        for (const auto& r_elem: rReferenceModelPart.Elements()) {
            const auto& r_geometry = r_elem.GetGeometry();
            CheckGeometryType(r_geometry.GetGeometryType());
            GenerateIntersection(r_geometry, Color, r_colors, r_skin_intersection);
        }
    } else if (rInputEntities == "conditions") {
        for (const auto& r_cond: rReferenceModelPart.Conditions()) {
            const auto& r_geometry = r_cond.GetGeometry();
            CheckGeometryType(r_geometry.GetGeometryType());
            GenerateIntersection(r_geometry, Color, r_colors, r_skin_intersection);
        }
    } else if (rInputEntities == "geometries") {
        for(auto& r_geometry : rReferenceModelPart.Geometries()) {
            CheckGeometryType(r_geometry.GetGeometryType());
            GenerateIntersection(r_geometry, Color, r_colors, r_skin_intersection);
        }
    } else {
        KRATOS_ERROR << "The input_entities  " << rInputEntities << " is not supported. The supported input_entities are  elements and conditions" << std::endl;
    }
}


void GenerateStlIntersectionWithCells::GenerateIntersection(
    const Geometry<Node>& rGeometry,
    int Color,
    const CartesianMeshColors& rMeshColors,
    SkinIntersection& rSkinIntersection)
{
    auto p_work_geometry = rSkinIntersection.ConvertToWorkGeometry(rGeometry);
    array_1d<std::size_t, 3> min_position(3,0);
    array_1d<std::size_t, 3> max_position(3,0);
    rMeshColors.CalculateOuterMinMaxNodePositions(rGeometry, min_position, max_position);
    for(std::size_t k = min_position[2] ; k < max_position[2] ; k++){
        for(std::size_t j = min_position[1] ; j < max_position[1] ; j++){
            for(std::size_t i = min_position[0] ; i < max_position[0] ; i++){
                Point cell_min_point = rMeshColors.GetPoint(i,j,k);
                Point cell_max_point = rMeshColors.GetPoint(i+1,j+1,k+1);
                if(rGeometry.HasIntersection(cell_min_point,cell_max_point)){
                    rSkinIntersection.ComputeCutInsideCell(p_work_geometry, cell_min_point, cell_max_point, i, j, k);
                }
            }
        }
    }

}

}