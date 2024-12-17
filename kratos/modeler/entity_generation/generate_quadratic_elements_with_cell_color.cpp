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
#include "generate_quadratic_elements_with_cell_color.h"

namespace Kratos {

Parameters GenerateQuadraticElementsWithCellColor::GetDefaultParameters() const
{
    return Parameters(R"({
        "type" : "quadratic_elements_with_cell_color",
        "model_part_name": "Undefined",
        "color": -1,
        "properties_id": 1,
        "generated_entity": "Element3D27N"
    })");
}

void GenerateQuadraticElementsWithCellColor::Generate(ModelPart& rModelPart, Parameters parameters)
{
    InitializeQuadraticData();

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

    Element::NodesArrayType cell_nodes(27);

    for (std::size_t k = 0; k < number_of_cells[2]; k++) {
        for (std::size_t j = 0; j < number_of_cells[1]; j++) {
            for (std::size_t i = 0; i < number_of_cells[0]; i++) {
                if(std::lround(r_colors.GetElementalColor(i,j,k)) == inside_color){
                    this->GetQuadraticCellNodes(cell_nodes, rModelPart, new_nodes, i, j, k);

                    //create the new element
                    Element::Pointer p_element = r_prototype_element.Create(element_id++, cell_nodes, p_properties);
                    new_elements.push_back(p_element);
                }
            }
        }
    }

    rModelPart.AddNodes(new_nodes.begin(), new_nodes.end());
    rModelPart.AddElements(new_elements.begin(), new_elements.end());
}

}