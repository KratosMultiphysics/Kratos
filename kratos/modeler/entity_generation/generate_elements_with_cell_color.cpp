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
#include "generate_elements_with_cell_color.h"

namespace Kratos {

Parameters GenerateElementsWithCellColor::GetDefaultParameters() const
{
    return Parameters(R"({
        "type" : "elements_with_cell_color",
        "model_part_name": "Undefined",
        "color": -1,
        "properties_id": 1,
        "generated_entity": "Element3D8N"
    })");
}

void GenerateElementsWithCellColor::Generate(ModelPart& rModelPart, Parameters parameters)
{
    const CartesianMeshColors& r_colors = GetMeshColors();
    std::size_t element_id = GetStartElementId();

    const int inside_color = parameters["color"].GetInt();
    std::size_t properties_id = parameters["properties_id"].GetInt();
    Properties::Pointer p_properties = CreateAndGetProperty(rModelPart, properties_id);

    array_1d<std::size_t, 3> number_of_cells = GetNumberOfCells();

    std::vector<ModelPart::NodeType::Pointer> new_nodes;
    std::vector<ModelPart::ElementType::Pointer> new_elements;

    auto& r_prototype_element = KratosComponents<Element>::Get(
        parameters["generated_entity"].GetString());

    Element::NodesArrayType cell_nodes(8);

    for (std::size_t k = 0; k < number_of_cells[2]; k++) {
        for (std::size_t j = 0; j < number_of_cells[1]; j++) {
            for (std::size_t i = 0; i < number_of_cells[0]; i++) {
                if(std::lround(r_colors.GetElementalColor(i,j,k)) == inside_color){
                    this->GetLinearCellNodes(cell_nodes, rModelPart, new_nodes, i, j, k);

                    //create the new element
                    Element::Pointer p_element = r_prototype_element.Create(element_id++, cell_nodes, p_properties);
                    new_elements.push_back(p_element);
                }
            }
        }
    }

    rModelPart.AddNodes(std::move(new_nodes));
    rModelPart.AddElements(std::move(new_elements));
}

}