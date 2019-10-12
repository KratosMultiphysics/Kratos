// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//					 license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher
//

// System includes
#include<unordered_map>

// External includes

// Project includes
#include "includes/model_part.h"
#include "custom_python/add_co_sim_io_to_python.h"

// CoSimIO
#include "co_simulation_api/co_sim_io.h"

namespace Kratos {
namespace Python {

namespace CoSimIO_Wrappers { // helpers namespace

void ExportGeometry(CoSim::CoSimIO& rCoSimIO, const ModelPart& rModelPart)
{
    KRATOS_ERROR << "This function is not yet implemented!" << std::endl;
}

void ImportGeometry(CoSim::CoSimIO& rCoSimIO, ModelPart& rModelPart)
{
    KRATOS_ERROR << "This function is not yet implemented!" << std::endl;
}

void ExportMesh(CoSim::CoSimIO& rCoSimIO, const ModelPart& rModelPart)
{
    KRATOS_ERROR << "This function is not yet implemented!" << std::endl;
}

void ImportMesh(CoSim::CoSimIO& rCoSimIO, ModelPart& rModelPart)
{
    std::vector<double> node_coords;
    std::vector<int> connectivities;
    std::vector<int> cell_types;
    CoSim::DataContainers::Mesh mesh = {node_coords, connectivities, cell_types};
    rCoSimIO.Import(mesh, rModelPart.Name());

    // fill ModelPart from received Mesh
    KRATOS_ERROR_IF(rModelPart.NumberOfNodes() > 0) << "ModelPart is not empty, it has nodes!" << std::endl;
    KRATOS_ERROR_IF(rModelPart.NumberOfProperties() > 0) << "ModelPart is not empty, it has properties!" << std::endl;
    KRATOS_ERROR_IF(rModelPart.IsDistributed()) << "ModelPart cannot be distributed!" << std::endl;

    const std::unordered_map<int, std::string> element_name_map = {
        // {1 , "Element3D1N"}, // does not yet exist
        {2 , "Element3D2N"},
        {3 , "Element3D3N"},
        {4 , "Element3D4N"}
    };
    const std::unordered_map<int, std::string> condition_name_map = {
        {1 , "PointCondition3D1N"},
        {2 , "LineCondition3D2N"},
        {3 , "SurfaceCondition3D3N"},
        {4 , "SurfaceCondition3D4N"}
    };

    // fill ModelPart with received entities

    const int num_nodes = node_coords.size()/3;
    for (int i=0; i<num_nodes; ++i) {
        rModelPart.CreateNewNode(i+1, node_coords[i*3], node_coords[i*3+1], node_coords[i*3+2]);
    }

    // auto p_props = rModelPart.CreateNewProperties(0);

    // int counter=0;
    // for (int i=0; i<numElems; ++i) {
    //     const int num_nodes_elem = (*numNodesPerElem)[i];
    //     std::vector<ModelPart::IndexType> elem_node_ids(num_nodes_elem);
    //     for (int j=0; j<num_nodes_elem; ++j) {
    //         elem_node_ids[j] = (*elem)[counter++];
    //     }
    //     if (UseConditions) {
    //         rModelPart.CreateNewCondition(condition_name_map.at(num_nodes_elem), i+1, elem_node_ids, p_props);
    //     } else {
    //         rModelPart.CreateNewElement(element_name_map.at(num_nodes_elem), i+1, elem_node_ids, p_props);
    //     }
    // }
}

void ExportData(CoSim::CoSimIO& rCoSimIO, const ModelPart& rModelPart)
{
    KRATOS_ERROR << "This function is not yet implemented!" << std::endl;
}

void ImportData(CoSim::CoSimIO& rCoSimIO, ModelPart& rModelPart)
{
    KRATOS_ERROR << "This function is not yet implemented!" << std::endl;
}

} // helpers namespace

void  AddCoSimIOToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef CoSim::CoSimIO::SettingsType SettingsType;
    typedef CoSim::CoSimIO CoSimIOType;

    py::class_<CoSimIOType>(m,"CoSimIO")
        .def(py::init<const std::string&, const std::string&>())
        .def(py::init<const std::string&, SettingsType>())

        .def("Connect",&CoSimIOType::Connect)
        .def("Disconnect",&CoSimIOType::Disconnect)

        .def("SendControlSignal",&CoSimIOType::SendControlSignal)
        // .def("RecvControlSignal",&CoSimIOType::RecvControlSignal) // not needed on CoSim side!

        .def("ImportGeometry", CoSimIO_Wrappers::ImportGeometry)
        .def("ExportGeometry", CoSimIO_Wrappers::ExportGeometry)

        .def("ImportMesh", CoSimIO_Wrappers::ImportMesh)
        .def("ExportMesh", CoSimIO_Wrappers::ExportMesh)

        .def("ImportData", CoSimIO_Wrappers::ImportData)
        .def("ExportData", CoSimIO_Wrappers::ExportData)
        ;

}

}  // namespace Python.
} // Namespace Kratos

