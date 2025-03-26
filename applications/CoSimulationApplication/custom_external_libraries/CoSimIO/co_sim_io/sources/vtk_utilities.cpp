//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes
#include <map>
#include <unordered_map>
#include <iomanip>
#include <fstream>

// Project includes
#include "includes/vtk_utilities.hpp"
#include "includes/utilities.hpp"

namespace CoSimIO {
namespace VtkUtilities {

VtkCellType GetVtkCellTypeForElementType(ElementType I_ElementType)
{
    const std::map<ElementType, VtkCellType> element_type_to_cell_Type_map = {
        {ElementType::Hexahedra3D20,        VtkCellType::Quadratic_Hexahedron},
        {ElementType::Hexahedra3D8,         VtkCellType::Hexahedron},
        {ElementType::Prism3D6,             VtkCellType::Wedge},
        {ElementType::Pyramid3D5,           VtkCellType::Pyramid},
        {ElementType::Quadrilateral2D4,     VtkCellType::Quad},
        {ElementType::Quadrilateral2D8,     VtkCellType::Quadratic_Quad},
        {ElementType::Quadrilateral3D4,     VtkCellType::Quad},
        {ElementType::Quadrilateral3D8,     VtkCellType::Quadratic_Quad},
        {ElementType::Tetrahedra3D10,       VtkCellType::Quadratic_Tetra},
        {ElementType::Tetrahedra3D4,        VtkCellType::Tetra},
        {ElementType::Triangle2D3,          VtkCellType::Triangle},
        {ElementType::Triangle2D6,          VtkCellType::Quadratic_Triangle},
        {ElementType::Triangle3D3,          VtkCellType::Triangle},
        {ElementType::Triangle3D6,          VtkCellType::Quadratic_Triangle},
        {ElementType::Line2D2,              VtkCellType::Line},
        {ElementType::Line2D3,              VtkCellType::Quadratic_Edge},
        {ElementType::Line3D2,              VtkCellType::Line},
        {ElementType::Line3D3,              VtkCellType::Quadratic_Edge},
        {ElementType::Point2D,              VtkCellType::Vertex},
        {ElementType::Point3D,              VtkCellType::Vertex}
    };

    auto type_iter = element_type_to_cell_Type_map.find(I_ElementType);
    CO_SIM_IO_ERROR_IF(type_iter == element_type_to_cell_Type_map.end()) << "Vtk does not support element type: " << Utilities::GetElementName(I_ElementType) << std::endl; // TODO maybe return -1 or so here in the future. This way also types not supported by vtk/Paraview could be used with CoSomIO (but not visualized in Paraview)
    return type_iter->second;
}

void WriteVtk(const Info& I_Settings, const ModelPart& I_ModelPart)
{
    const std::string file_name = I_Settings.Get<std::string>("file_name");

    const int precision = I_Settings.Get<int>("precision", 12);

    std::ofstream output_file;
    output_file.open(file_name);
    Utilities::CheckStream(output_file, file_name);

    output_file << std::scientific << std::setprecision(precision);

    // write file header
    output_file << "# vtk DataFile Version 4.0\n";
    output_file << "CoSimIO FileCommunication\n";
    output_file << "ASCII\n";
    output_file << "DATASET UNSTRUCTURED_GRID\n\n";

    // write nodes and create Id map
    std::unordered_map<IdType, IdType> id_map;
    IdType vtk_id = 0;
    output_file << "POINTS " << I_ModelPart.NumberOfNodes() << " float\n";
    for (const auto& r_node : I_ModelPart.Nodes()) {
        output_file << r_node.X() << " " << r_node.Y() << " " << r_node.Z() << "\n";
        id_map[r_node.Id()] = vtk_id++;
    }
    output_file << "\n";

    // get cells size information
    std::size_t cell_list_size = 0;
    for (const auto& r_elem : I_ModelPart.Elements()) {
        cell_list_size += r_elem.NumberOfNodes() + 1; // +1 for size of connectivity
    }

    // write cells connectivity
    const auto const_id_map = id_map; // const reference to not accidentally modify the map
    output_file << "CELLS " << I_ModelPart.NumberOfElements() << " " << cell_list_size << "\n";
    for (const auto& r_elem : I_ModelPart.Elements()) {
        const std::size_t num_nodes_cell = r_elem.NumberOfNodes();
        output_file << num_nodes_cell << " ";
        std::size_t node_counter = 0;
        for (const auto& r_node : r_elem.Nodes()) {
            const IdType node_id = r_node.Id();
            auto id_iter = const_id_map.find(node_id);
            CO_SIM_IO_ERROR_IF(id_iter == const_id_map.end()) << "The node with Id " << node_id << " is not part of the ModelPart but used for Element with Id " << r_elem.Id() << std::endl;
            output_file << id_iter->second;
            if (node_counter++<num_nodes_cell-1) output_file << " "; // not adding a whitespace after last number
        }
        output_file << "\n";
    }
    output_file << "\n";

    // write cell types
    output_file << "CELL_TYPES " << I_ModelPart.NumberOfElements() << "\n";
    for (const auto& r_elem : I_ModelPart.Elements()) {
        output_file << static_cast<int>(GetVtkCellTypeForElementType(r_elem.Type())) << "\n";
    }
    output_file << "\n";

    // writing node Ids
    output_file << "POINT_DATA " << I_ModelPart.NumberOfNodes() << "\n";
    output_file << "FIELD FieldData 1" << "\n";
    output_file << "NODE_ID 1 " << I_ModelPart.NumberOfNodes() << " int\n";
    for (const auto& r_node : I_ModelPart.Nodes()) {
        output_file << r_node.Id() << "\n";
    }
    output_file << "\n";

    // writing element Ids
    output_file << "CELL_DATA " << I_ModelPart.NumberOfElements() << "\n";
    output_file << "FIELD FieldData 1" << "\n";
    output_file << "ELEMENT_ID 1 " << I_ModelPart.NumberOfElements() << " int\n";
    for (const auto& r_elem : I_ModelPart.Elements()) {
        output_file << r_elem.Id() << "\n";
    }
    output_file << "\n";

    // writing element types
    output_file << "CELL_DATA " << I_ModelPart.NumberOfElements() << "\n";
    output_file << "FIELD FieldData 1" << "\n";
    output_file << "ELEMENT_TYPE 1 " << I_ModelPart.NumberOfElements() << " int\n";
    for (const auto& r_elem : I_ModelPart.Elements()) {
        output_file << static_cast<int>(r_elem.Type()) << "\n";
    }

    output_file.close();
}

void ReadVtk(const Info& I_Settings, ModelPart& O_ModelPart)
{
    const std::string file_name = I_Settings.Get<std::string>("file_name");

    std::ifstream input_file(file_name);
    Utilities::CheckStream(input_file, file_name);

    // reading file
    std::string current_line;
    std::vector<double> nodal_coords;
    std::vector<IdType> nodal_ids;
    std::vector<IdType> element_ids;
    std::vector<ElementType> element_types;
    std::vector<ConnectivitiesType> element_connectivities;

    while (std::getline(input_file, current_line)) {
        // reading nodes
        if (current_line.find("POINTS") != std::string::npos) {
            std::size_t num_nodes;
            current_line = current_line.substr(current_line.find("POINTS") + 7); // removing "POINTS"
            std::istringstream line_stream(current_line);
            line_stream >> num_nodes;

            nodal_coords.resize(3*num_nodes);
            nodal_ids.resize(num_nodes);

            for (std::size_t i=0; i<num_nodes*3; ++i) {
                input_file >> nodal_coords[i];
            }
        }

        // reading connectivities
        if (current_line.find("CELLS") != std::string::npos) {
            std::size_t num_elems, num_nodes_per_elem;
            current_line = current_line.substr(current_line.find("CELLS") + 6); // removing "CELLS"
            std::istringstream line_stream(current_line);
            line_stream >> num_elems;

            element_ids.resize(num_elems);
            element_types.resize(num_elems);
            element_connectivities.resize(num_elems);

            for (std::size_t i=0; i<num_elems; ++i) {
                input_file >> num_nodes_per_elem;
                element_connectivities[i].resize(num_nodes_per_elem);
                for (std::size_t j=0; j<num_nodes_per_elem; ++j) {
                    input_file >> element_connectivities[i][j];
                }
            }
        }

        // reading node Ids
        if (current_line.find("NODE_ID") != std::string::npos) {
            for (std::size_t i=0; i<nodal_ids.size(); ++i) { // nodal_ids was resized to correct size above
                input_file >> nodal_ids[i];
            }
        }

        // reading element Ids
        if (current_line.find("ELEMENT_ID") != std::string::npos) {
            for (std::size_t i=0; i<element_ids.size(); ++i) { // element_ids was resized to correct size above
                input_file >> element_ids[i];
            }
        }

        // reading element types
        if (current_line.find("ELEMENT_TYPE") != std::string::npos) {
            int enum_temp;
            for (std::size_t i=0; i<element_types.size(); ++i) { // element_types was resized to correct size above
                input_file >> enum_temp; // using a temp variable as enums cannot be read directly
                element_types[i] = static_cast<CoSimIO::ElementType>(enum_temp);
            }
        }
    }

    // filling ModelPart with read information
    for (std::size_t i=0; i<nodal_ids.size(); ++i) {
        O_ModelPart.CreateNewNode(
            nodal_ids[i],
            nodal_coords[i*3],
            nodal_coords[i*3+1],
            nodal_coords[i*3+2]);
    }
    for (std::size_t i=0; i<element_ids.size(); ++i) {
        for (auto& conn : element_connectivities[i]) {
            conn = nodal_ids[conn]; // transforming vtk Ids back to original Ids
        }
        O_ModelPart.CreateNewElement(
            element_ids[i],
            element_types[i],
            element_connectivities[i]);
    }

    input_file.close();
}

} // namespace VtkUtilities
} // namespace CoSimIO
