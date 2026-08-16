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
#include "generate_triangular_conditions_with_face_color.h"

namespace Kratos {

Parameters GenerateTriangularConditionsWithFaceColor::GetDefaultParameters() const
{
    return Parameters(R"({
        "type" : "triangular_conditions_with_face_color",
        "model_part_name": "Undefined",
        "color": -1,
        "properties_id": 1,
        "generated_entity": "SurfaceCondition3D3N"
    })");
}

void GenerateTriangularConditionsWithFaceColor::Generate(ModelPart& rModelPart, Parameters parameters)
{
    const CartesianMeshColors& r_colors = GetMeshColors();
    std::size_t condition_id = GetStartConditionId();

    const int inside_color = parameters["color"].GetInt();
    const std::size_t properties_id = parameters["properties_id"].GetInt();
    Properties::Pointer p_properties = CreateAndGetProperty(rModelPart, properties_id);

    const std::size_t x_offset[24]={0,0,0,0, 0,1,1,0, 0,0,1,1, 1,1,1,1, 0,0,1,1, 0,1,1,0};
    const std::size_t y_offset[24]={0,0,1,1, 0,0,0,0, 0,1,1,0, 0,1,1,0, 1,1,1,1, 0,0,1,1};
    const std::size_t z_offset[24]={0,1,1,0, 0,0,1,1, 0,0,0,0, 0,0,1,1, 0,1,1,0, 1,1,1,1};

    array_1d<std::size_t, 3> number_of_cells = GetNumberOfCells();

    std::vector<ModelPart::NodeType::Pointer> new_nodes;
    std::vector<ModelPart::ConditionType::Pointer> new_conditions;
    auto& r_prototype_condition = KratosComponents<Condition>::Get(
        parameters["generated_entity"].GetString());

    Element::NodesArrayType face_nodes_1(3);
    Element::NodesArrayType face_nodes_2(3);
    for (std::size_t k = 0; k < number_of_cells[2]; k++) {
        for (std::size_t j = 0; j < number_of_cells[1]; j++) {
            for (std::size_t i = 0; i < number_of_cells[0]; i++) {
                auto& r_faces_color = r_colors.GetElementalFaceColor(i,j,k);
                for(std::size_t i_face = 0; i_face < 6; i_face++){
                    if(std::lround(r_faces_color[i_face]) == inside_color){
                        std::size_t base_index=i_face*4;
                        face_nodes_1(0) = GenerateOrRetrieveNode(rModelPart, new_nodes, i+x_offset[base_index]  , j+y_offset[base_index]    , k+z_offset[base_index]  );
                        face_nodes_1(1) = GenerateOrRetrieveNode(rModelPart, new_nodes, i+x_offset[base_index+1], j+y_offset[base_index+1]  , k+z_offset[base_index+1]);
                        face_nodes_1(2) = GenerateOrRetrieveNode(rModelPart, new_nodes, i+x_offset[base_index+2], j+y_offset[base_index+2]  , k+z_offset[base_index+2]);

                        face_nodes_2(0) = GenerateOrRetrieveNode(rModelPart, new_nodes, i+x_offset[base_index]  , j+y_offset[base_index]    , k+z_offset[base_index]  );
                        face_nodes_2(1) = GenerateOrRetrieveNode(rModelPart, new_nodes, i+x_offset[base_index+2], j+y_offset[base_index+2]  , k+z_offset[base_index+2]);
                        face_nodes_2(2) = GenerateOrRetrieveNode(rModelPart, new_nodes, i+x_offset[base_index+3], j+y_offset[base_index+3]  , k+z_offset[base_index+3]);
                        //create the new conditions
                        Condition::Pointer p_condition_1 = r_prototype_condition.Create(condition_id++, face_nodes_1, p_properties);
                        new_conditions.push_back(p_condition_1);
                        Condition::Pointer p_condition_2 = r_prototype_condition.Create(condition_id++, face_nodes_2, p_properties);
                        new_conditions.push_back(p_condition_2);
                    }
                }
            }
        }
    }

    rModelPart.AddNodes(std::move(new_nodes));
    rModelPart.AddConditions(std::move(new_conditions));
}

}