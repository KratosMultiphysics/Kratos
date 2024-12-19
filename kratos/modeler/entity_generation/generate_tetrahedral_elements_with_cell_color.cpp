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
#include "geometries/tetrahedra_3d_4.h"
#include "generate_tetrahedral_elements_with_cell_color.h"

namespace Kratos {

Parameters GenerateTetrahedralElementsWithCellColor::GetDefaultParameters() const
{
    return Parameters(R"({
        "type" : "tetrahedral_elements_with_cell_color",
        "model_part_name": "Undefined",
        "color": -1,
        "properties_id": 1,
        "generated_entity": "Element3D4N"
    })");
}

void GenerateTetrahedralElementsWithCellColor::Generate(ModelPart& rModelPart, Parameters parameters)
{
    const CartesianMeshColors& r_colors = GetMeshColors();
    std::size_t element_id = GetStartElementId();

    const int inside_color = parameters["color"].GetInt();
    const std::size_t properties_id = parameters["properties_id"].GetInt();
    Properties::Pointer p_properties = CreateAndGetProperty(rModelPart, properties_id);

    auto& r_prototype_element = KratosComponents<Element>::Get(
        parameters["generated_entity"].GetString());

    array_1d<std::size_t, 3> number_of_cells = GetNumberOfCells();

    std::vector<ModelPart::NodeType::Pointer> new_nodes;
    std::vector<ModelPart::ElementType::Pointer> new_elements;

    Element::NodesArrayType cell_nodes(8);
    for (std::size_t k = 0; k < number_of_cells[2]; k++) {
        for (std::size_t j = 0; j < number_of_cells[1]; j++) {
            for (std::size_t i = 0; i < number_of_cells[0]; i++) {
                if(std::lround(r_colors.GetElementalColor(i,j,k)) == inside_color){
                    GetLinearCellNodes(cell_nodes, rModelPart, new_nodes, i, j, k);
                    CreateTetrahedraInCell(new_elements, cell_nodes, element_id, p_properties, r_prototype_element);
                    element_id += 6;
                }
            }
        }
    }

    rModelPart.AddNodes(std::move(new_nodes));
    rModelPart.AddElements(std::move(new_elements));
}

void GenerateTetrahedralElementsWithCellColor::CreateTetrahedraInCell(
    std::vector<ModelPart::ElementType::Pointer>& rElements,
    Element::NodesArrayType& rCellNodes,
    const std::size_t StartId,
    Properties::Pointer& pProperties,
    const Element& rPrototype)
{
    constexpr std::size_t nodes_per_tet = 4;
    constexpr std::size_t tets_per_cell = 6;
    using ConnectivityList = std::array<std::array<std::size_t, nodes_per_tet>, tets_per_cell>;
    static constexpr ConnectivityList local_connectivities {{ { 0,3,6,2 },
                                                                { 3,6,7,0 },
                                                                { 4,7,6,0 },
                                                                { 0,4,5,6 },
                                                                { 0,1,2,6 },
                                                                { 1,5,6,0 } }};

    for(std::size_t i = 0; i < tets_per_cell; i++) {
        const auto & connectivty = local_connectivities[i];

        auto geom = Kratos::make_shared<Tetrahedra3D4<Node>>(
            rCellNodes(connectivty[0]),
            rCellNodes(connectivty[1]),
            rCellNodes(connectivty[2]),
            rCellNodes(connectivty[3]));

        rElements.push_back(rPrototype.Create(StartId+i, geom, pProperties));
    }
}

}