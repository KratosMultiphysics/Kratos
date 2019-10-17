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
// defining the macros that are used within the CoSimIO to point the the Kratos-native macros
// Note that this has to be done BEFORE including the header the first time
#define KRATOS_CO_SIM_INFO KRATOS_INFO
#define KRATOS_CO_SIM_INFO_IF KRATOS_INFO_IF
#define KRATOS_CO_SIM_ERROR KRATOS_ERROR
#define KRATOS_CO_SIM_ERROR_IF KRATOS_ERROR_IF
#define KRATOS_CO_SIM_ERROR_IF_NOT KRATOS_ERROR_IF_NOT
#include "co_simulation_api/co_sim_io.h"

namespace Kratos {
namespace Python {

namespace CoSimIO_Wrappers { // helpers namespace

enum class DataLocation { NodeHistorical, NodeNonHistorical, Element, Condition, ModelPart };

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

void ExportData_Scalar(CoSim::CoSimIO& rCoSimIO, const ModelPart& rModelPart, const Variable<double>& rVariable, const DataLocation DataLoc, const std::string& rIdentifier)
{
    std::vector<double> data_vals;

    if (DataLoc == DataLocation::NodeHistorical) {
        std::size_t counter = 0;
        data_vals.resize(rModelPart.NumberOfNodes());
        for (const auto& r_node : rModelPart.Nodes()) {
            data_vals[counter++] = r_node.FastGetSolutionStepValue(rVariable);
        }

    } else if (DataLoc == DataLocation::NodeNonHistorical) {
        std::size_t counter = 0;
        data_vals.resize(rModelPart.NumberOfNodes());
        for (const auto& r_node : rModelPart.Nodes()) {
            data_vals[counter++] = r_node.GetValue(rVariable);
        }

    } else if (DataLoc == DataLocation::Element) {
        std::size_t counter = 0;
        data_vals.resize(rModelPart.NumberOfElements());
        for (const auto& r_elem : rModelPart.Elements()) {
            data_vals[counter++] = r_elem.GetValue(rVariable);
        }

    } else if (DataLoc == DataLocation::Condition) {
        std::size_t counter = 0;
        data_vals.resize(rModelPart.NumberOfConditions());
        for (const auto& r_cond : rModelPart.Conditions()) {
            data_vals[counter++] = r_cond.GetValue(rVariable);
        }

    } else if (DataLoc == DataLocation::ModelPart) {
        data_vals.resize(1);
        data_vals[0] = rModelPart[rVariable];
    }

    CoSim::DataContainers::Data data_container = {data_vals};
    rCoSimIO.Export(data_container, rIdentifier);
}

void ImportData_Scalar(CoSim::CoSimIO& rCoSimIO, ModelPart& rModelPart, const Variable<double>& rVariable, const DataLocation DataLoc, const std::string& rIdentifier)
{
    std::vector<double> data_vals;
    CoSim::DataContainers::Data data_container = {data_vals};
    rCoSimIO.Import(data_container, rIdentifier);

    // TODO implement size-checks

    if (DataLoc == DataLocation::NodeHistorical) {
        std::size_t counter = 0;
        for (auto& r_node : rModelPart.Nodes()) {
            r_node.FastGetSolutionStepValue(rVariable) = data_vals[counter++];
        }

    } else if (DataLoc == DataLocation::NodeNonHistorical) {
        std::size_t counter = 0;
        for (auto& r_node : rModelPart.Nodes()) {
            r_node.GetValue(rVariable) = data_vals[counter++];
        }

    } else if (DataLoc == DataLocation::Element) {
        std::size_t counter = 0;
        for (auto& r_elem : rModelPart.Elements()) {
            r_elem.GetValue(rVariable) = data_vals[counter++];
        }

    } else if (DataLoc == DataLocation::Condition) {
        std::size_t counter = 0;
        for (auto& r_cond : rModelPart.Conditions()) {
            r_cond.GetValue(rVariable) = data_vals[counter++];
        }

    } else if (DataLoc == DataLocation::ModelPart) {
        rModelPart[rVariable] = data_vals[0];
    }
}

void ExportData_Vector(CoSim::CoSimIO& rCoSimIO, const ModelPart& rModelPart, const Variable< array_1d<double, 3> >& rVariable, const DataLocation DataLoc, const std::string& rIdentifier)
{
    std::vector<double> data_vals;

    if (DataLoc == DataLocation::NodeHistorical) {
        std::size_t counter = 0;
        data_vals.resize(rModelPart.NumberOfNodes()*3);
        for (const auto& r_node : rModelPart.Nodes()) {
            const array_1d<double, 3>& r_val = r_node.FastGetSolutionStepValue(rVariable);
            data_vals[counter++] = r_val[0];
            data_vals[counter++] = r_val[1];
            data_vals[counter++] = r_val[2];
        }

    } else if (DataLoc == DataLocation::NodeNonHistorical) {
        std::size_t counter = 0;
        data_vals.resize(rModelPart.NumberOfNodes()*3);
        for (const auto& r_node : rModelPart.Nodes()) {
            const array_1d<double, 3>& r_val = r_node.GetValue(rVariable);
            data_vals[counter++] = r_val[0];
            data_vals[counter++] = r_val[1];
            data_vals[counter++] = r_val[2];
        }

    } else if (DataLoc == DataLocation::Element) {
        std::size_t counter = 0;
        data_vals.resize(rModelPart.NumberOfElements()*3);
        for (const auto& r_elem : rModelPart.Elements()) {
            const array_1d<double, 3>& r_val = r_elem.GetValue(rVariable);
            data_vals[counter++] = r_val[0];
            data_vals[counter++] = r_val[1];
            data_vals[counter++] = r_val[2];
        }

    } else if (DataLoc == DataLocation::Condition) {
        std::size_t counter = 0;
        data_vals.resize(rModelPart.NumberOfConditions()*3);
        for (const auto& r_cond : rModelPart.Conditions()) {
            const array_1d<double, 3>& r_val = r_cond.GetValue(rVariable);
            data_vals[counter++] = r_val[0];
            data_vals[counter++] = r_val[1];
            data_vals[counter++] = r_val[2];
        }

    } else if (DataLoc == DataLocation::ModelPart) {
        data_vals.resize(3);
        const auto& r_val = rModelPart[rVariable];
        data_vals[0] = r_val[0];
        data_vals[1] = r_val[1];
        data_vals[2] = r_val[2];
    }

    CoSim::DataContainers::Data data_container = {data_vals};
    rCoSimIO.Export(data_container, rIdentifier);
}

void ImportData_Vector(CoSim::CoSimIO& rCoSimIO, ModelPart& rModelPart, const Variable< array_1d<double, 3> >& rVariable, const DataLocation DataLoc, const std::string& rIdentifier)
{
    std::vector<double> data_vals;
    CoSim::DataContainers::Data data_container = {data_vals};
    rCoSimIO.Import(data_container, rIdentifier);

    // TODO implement size-checks

    if (DataLoc == DataLocation::NodeHistorical) {
        std::size_t counter = 0;
        for (auto& r_node : rModelPart.Nodes()) {
            array_1d<double, 3>& r_val = r_node.FastGetSolutionStepValue(rVariable);
            r_val[0] = data_vals[counter++];
            r_val[1] = data_vals[counter++];
            r_val[2] = data_vals[counter++];
        }

    } else if (DataLoc == DataLocation::NodeNonHistorical) {
        std::size_t counter = 0;
        for (auto& r_node : rModelPart.Nodes()) {
            array_1d<double, 3>& r_val = r_node.GetValue(rVariable);
            r_val[0] = data_vals[counter++];
            r_val[1] = data_vals[counter++];
            r_val[2] = data_vals[counter++];
        }

    } else if (DataLoc == DataLocation::Element) {
        std::size_t counter = 0;
        for (auto& r_elem : rModelPart.Elements()) {
            array_1d<double, 3>& r_val = r_elem.GetValue(rVariable);
            r_val[0] = data_vals[counter++];
            r_val[1] = data_vals[counter++];
            r_val[2] = data_vals[counter++];
        }

    } else if (DataLoc == DataLocation::Condition) {
        std::size_t counter = 0;
        for (auto& r_cond : rModelPart.Conditions()) {
            array_1d<double, 3>& r_val = r_cond.GetValue(rVariable);
            r_val[0] = data_vals[counter++];
            r_val[1] = data_vals[counter++];
            r_val[2] = data_vals[counter++];
        }

    } else if (DataLoc == DataLocation::ModelPart) {
        auto& r_val = rModelPart[rVariable];
        r_val[0] = data_vals[0];
        r_val[1] = data_vals[1];
        r_val[2] = data_vals[2];
    }
}


void ImportData_PyList(CoSim::CoSimIO& rCoSimIO, pybind11::list& PyList, const std::string& rIdentifier)
{
    std::vector<double> data_vals;
    CoSim::DataContainers::Data data_container = {data_vals};
    rCoSimIO.Import(data_container, rIdentifier);

    KRATOS_ERROR_IF(PyList.size() != data_vals.size()) << "The size of the list has to be specified before, expected size of " << data_vals.size() << ", current size: " << PyList.size() << std::endl;

    // copy back the received values
    for (std::size_t i=0; i<data_vals.size(); ++i) {
        PyList[i] = data_vals[i];
    }

}
void ExportData_PyList(CoSim::CoSimIO& rCoSimIO, std::vector<double>& rValues, const std::string& rIdentifier)
{
    CoSim::DataContainers::Data data_container = {rValues};
    rCoSimIO.Export(data_container, rIdentifier);
}

} // helpers namespace

void  AddCoSimIOToPython(pybind11::module& m)
{
    namespace py = pybind11;

    typedef CoSim::CoSimIO::SettingsType SettingsType;
    typedef CoSim::CoSimIO CoSimIOType;

    py::class_<CoSimIOType>(m,"CoSimIO")
        .def(py::init<const std::string&, const std::string&, const bool>())
        .def(py::init<const std::string&, SettingsType, const bool>())

        .def("Connect",&CoSimIOType::Connect)
        .def("Disconnect",&CoSimIOType::Disconnect)

        .def("SendControlSignal",&CoSimIOType::SendControlSignal)
        // .def("RecvControlSignal",&CoSimIOType::RecvControlSignal) // not needed on CoSim side! (also would not work with direct exposure due to pass-by-ref)

        .def("ImportGeometry", CoSimIO_Wrappers::ImportGeometry)
        .def("ExportGeometry", CoSimIO_Wrappers::ExportGeometry)

        .def("ImportMesh", CoSimIO_Wrappers::ImportMesh)
        .def("ExportMesh", CoSimIO_Wrappers::ExportMesh)

        .def("ImportData", CoSimIO_Wrappers::ImportData_Scalar)
        .def("ExportData", CoSimIO_Wrappers::ExportData_Scalar)
        .def("ImportData", CoSimIO_Wrappers::ImportData_Vector)
        .def("ExportData", CoSimIO_Wrappers::ExportData_Vector)
        .def("ImportData", CoSimIO_Wrappers::ImportData_PyList)
        .def("ExportData", CoSimIO_Wrappers::ExportData_PyList)
        ;

    py::enum_<CoSimIO_Wrappers::DataLocation>(m,"DataLocation")
        .value("NodeHistorical", CoSimIO_Wrappers::DataLocation::NodeHistorical)
        .value("NodeNonHistorical", CoSimIO_Wrappers::DataLocation::NodeNonHistorical)
        .value("Element", CoSimIO_Wrappers::DataLocation::Element)
        .value("Condition", CoSimIO_Wrappers::DataLocation::Condition)
        .value("ModelPart", CoSimIO_Wrappers::DataLocation::ModelPart)
        ;

    py::enum_<CoSim::Internals::ControlSignal>(m,"ControlSignal")
        .value("Dummy", CoSim::Internals::ControlSignal::Dummy)
        .value("BreakSolutionLoop", CoSim::Internals::ControlSignal::BreakSolutionLoop)
        .value("ConvergenceAchieved", CoSim::Internals::ControlSignal::ConvergenceAchieved)

        .value("AdvanceInTime", CoSim::Internals::ControlSignal::AdvanceInTime)
        .value("InitializeSolutionStep", CoSim::Internals::ControlSignal::InitializeSolutionStep)
        .value("SolveSolutionStep", CoSim::Internals::ControlSignal::SolveSolutionStep)
        .value("FinalizeSolutionStep", CoSim::Internals::ControlSignal::FinalizeSolutionStep)

        .value("ImportGeometry", CoSim::Internals::ControlSignal::ImportGeometry)
        .value("ExportGeometry", CoSim::Internals::ControlSignal::ExportGeometry)
        .value("ImportMesh", CoSim::Internals::ControlSignal::ImportMesh)
        .value("ExportMesh", CoSim::Internals::ControlSignal::ExportMesh)
        .value("ImportData", CoSim::Internals::ControlSignal::ImportData)
        .value("ExportData", CoSim::Internals::ControlSignal::ExportData)
        ;
}

}  // namespace Python.
} // Namespace Kratos

