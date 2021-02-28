// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//					 license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// System includes
#include <unordered_map>
#include <map>
#include <limits>

// External includes

// Project includes
#include "includes/model_part.h"
#include "custom_python/add_co_sim_io_to_python.h"
#include "utilities/auxiliar_model_part_utilities.h"

// CoSimIO
// defining the macros that are used within the CoSimIO to point the the Kratos-native macros
// Note that this has to be done BEFORE including the header the first time
// see co_sim_io/impl/macros.hpp
#define CO_SIM_IO_INFO          KRATOS_INFO
#define CO_SIM_IO_INFO_IF       KRATOS_INFO_IF
#define CO_SIM_IO_ERROR         KRATOS_ERROR
#define CO_SIM_IO_ERROR_IF      KRATOS_ERROR_IF
#define CO_SIM_IO_ERROR_IF_NOT  KRATOS_ERROR_IF_NOT
#include "custom_external_libraries/co_sim_io/co_sim_io.hpp"

namespace Kratos {
namespace Python {

// BIG TODO: make OMP parallel most loops
// TOD use elements or conditions?? => how to switch?
namespace CoSimIO_Wrappers { // helpers namespace

// creating static buffers such that memory does not constantly have to be reallocated during the data-exchange
// this is not threadsafe, but the functions should never be executed by different threads at the same time (we don't do shared-memory parallelism from python)
// defining the static members
struct DataBuffers {
    static std::vector<double> vector_doubles;
    static std::vector<int> vector_connectivities;
    static std::vector<int> vector_types;
};

// declaring the static members
std::vector<double> DataBuffers::vector_doubles;
std::vector<int> DataBuffers::vector_connectivities;
std::vector<int> DataBuffers::vector_types;


void ExportMesh(
    CoSimIO::Info& rInfo,
    const ModelPart& rModelPart)
{
    // extract information from ModelPart
    const int num_nodes = rModelPart.NumberOfNodes();
    const int num_elems = rModelPart.NumberOfElements();

    DataBuffers::vector_doubles.resize(num_nodes*3);
    DataBuffers::vector_types.resize(num_elems);
    DataBuffers::vector_connectivities.clear();

    std::size_t node_counter = 0;
    for (const auto& r_node : rModelPart.Nodes()) {
        const auto& r_coords = r_node.GetInitialPosition(); // TODO or current coords?
        DataBuffers::vector_doubles[node_counter++] = r_coords[0];
        DataBuffers::vector_doubles[node_counter++] = r_coords[1];
        DataBuffers::vector_doubles[node_counter++] = r_coords[2];
    }

    // NOTE: See https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
    const std::map<GeometryData::KratosGeometryType, int> geo_type_vtk_cell_type_map = {
        { GeometryData::KratosGeometryType::Kratos_Point2D,          1 },
        { GeometryData::KratosGeometryType::Kratos_Point3D,          1 },
        { GeometryData::KratosGeometryType::Kratos_Line2D2,          3 },
        { GeometryData::KratosGeometryType::Kratos_Line3D2,          3 },
        { GeometryData::KratosGeometryType::Kratos_Triangle2D3,      5 },
        { GeometryData::KratosGeometryType::Kratos_Triangle3D3,      5 },
        { GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4, 9 },
        { GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4, 9 },
        { GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4,    10 },
        { GeometryData::KratosGeometryType::Kratos_Hexahedra3D8,     12 }
    };

    std::size_t elem_counter = 0;
    int connectivities_offset = std::numeric_limits<int>::max(); //in paraview the connectivities start from 0, hence we have to check beforehand what is the connectivities offset
    for (const auto& r_elem : rModelPart.Elements()) {
        const auto& r_geom = r_elem.GetGeometry();
        DataBuffers::vector_types[elem_counter++] = geo_type_vtk_cell_type_map.at(r_geom.GetGeometryType()); // TODO add check

        for (const auto& r_node : r_geom) {
            const int current_id = r_node.Id();
            DataBuffers::vector_connectivities.push_back(current_id);
            connectivities_offset = std::min(connectivities_offset, current_id);
        }
    }

    if (connectivities_offset != 0) {
        for(auto& r_connectivity : DataBuffers::vector_connectivities) {
            r_connectivity -= connectivities_offset;
        }
    }

    CoSimIO::ExportMesh(
        rInfo,
        DataBuffers::vector_doubles,
        DataBuffers::vector_connectivities,
        DataBuffers::vector_types);
}

void ImportMesh(
    CoSimIO::Info& rInfo,
    ModelPart& rModelPart)
{
    CoSimIO::ImportMesh(
        rInfo,
        DataBuffers::vector_doubles,
        DataBuffers::vector_connectivities,
        DataBuffers::vector_types);

    // fill ModelPart from received Mesh
    KRATOS_ERROR_IF(rModelPart.NumberOfNodes() > 0) << "ModelPart is not empty, it has nodes!" << std::endl;
    KRATOS_ERROR_IF(rModelPart.NumberOfProperties() > 0) << "ModelPart is not empty, it has properties!" << std::endl;
    KRATOS_ERROR_IF(rModelPart.IsDistributed()) << "ModelPart cannot be distributed!" << std::endl;

    const std::unordered_map<int, std::pair<int, std::string>> vtk_type_map = {
        // {1 , "Element3D1N"}, // does not yet exist
        {3 ,  {2, "Element2D2N"} }, // line
        {5 ,  {3, "Element2D3N"} }, // triangle
        {9 ,  {4, "Element2D4N"} }, // quad
        {10 , {4, "Element3D4N"} }, // tetra
        {12 , {8, "Element3D8N"} }  // hexa
    };

    // fill ModelPart with received entities // TODO do this in OMP and add only after creation!

    for (std::size_t i=0; i<DataBuffers::vector_doubles.size()/3; ++i) {
        rModelPart.CreateNewNode(
            i, // TODO check this, this will create nodes with Id=0
            DataBuffers::vector_doubles[i*3],
            DataBuffers::vector_doubles[i*3+1],
            DataBuffers::vector_doubles[i*3+2]);
    }

    auto p_props = rModelPart.CreateNewProperties(0);

    int counter=0;
    for (std::size_t i=0; i<DataBuffers::vector_types.size(); ++i) {
        const auto& vtk_type_info = vtk_type_map.at(DataBuffers::vector_types[i]); // TODO add check
        const int num_nodes_elem = vtk_type_info.first;
        std::vector<ModelPart::IndexType> elem_node_ids(num_nodes_elem);
        for (int j=0; j<num_nodes_elem; ++j) {
            elem_node_ids[j] = DataBuffers::vector_connectivities[counter++];
        }
        rModelPart.CreateNewElement(vtk_type_info.second, i+1, elem_node_ids, p_props);
    }
}

void ImportDataSizeCheck(const std::size_t ExpectedSize, const std::size_t ImportedSize)
{
    KRATOS_ERROR_IF(ExpectedSize != ImportedSize) << "Expected to import " << ExpectedSize << " values but got " << ImportedSize << " values instead!" << std::endl;
}

void ExportData_ModelPart_Scalar(
    CoSimIO::Info& rInfo,
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const DataLocation DataLoc)
{
    KRATOS_TRY

    AuxiliarModelPartUtilities(rModelPart).GetScalarData<double>(rVariable, (Kratos::DataLocation)DataLoc, DataBuffers::vector_doubles);
    CoSimIO::ExportData(rInfo, DataBuffers::vector_doubles);

    KRATOS_CATCH("")
}

void ImportData_ModelPart_Scalar(
    CoSimIO::Info& rInfo,
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const DataLocation DataLoc)
{
    KRATOS_TRY

    CoSimIO::ImportData(rInfo, DataBuffers::vector_doubles);

    AuxiliarModelPartUtilities(rModelPart).SetScalarData<double>(rVariable, (Kratos::DataLocation)DataLoc, DataBuffers::vector_doubles);
    KRATOS_CATCH("")
}

void ExportData_ModelPart_Vector(
    CoSimIO::Info& rInfo,
    ModelPart& rModelPart,
    const Variable< array_1d<double, 3> >& rVariable,
    const DataLocation DataLoc)
{
    KRATOS_TRY

    AuxiliarModelPartUtilities(rModelPart).GetVectorData< array_1d<double, 3> >(rVariable, (Kratos::DataLocation)DataLoc, DataBuffers::vector_doubles);
    CoSimIO::ExportData(rInfo, DataBuffers::vector_doubles);

    KRATOS_CATCH("")
}

void ImportData_ModelPart_Vector(
    CoSimIO::Info& rInfo,
    ModelPart& rModelPart,
    const Variable< array_1d<double, 3> >& rVariable,
    const DataLocation DataLoc)
{
    KRATOS_TRY

    CoSimIO::ImportData(rInfo, DataBuffers::vector_doubles);

    AuxiliarModelPartUtilities(rModelPart).SetVectorData< array_1d<double, 3> >(rVariable, (Kratos::DataLocation)DataLoc, DataBuffers::vector_doubles);
    KRATOS_CATCH("")
}

std::vector<double> ImportData_RawValues(
    CoSimIO::Info& rInfo)
{
    KRATOS_TRY

    CoSimIO::ImportData(rInfo, DataBuffers::vector_doubles);

    return DataBuffers::vector_doubles;

    KRATOS_CATCH("")
}

void ExportData_RawValues(
    CoSimIO::Info& rInfo,
    const std::vector<double>& rValues)
{
    KRATOS_TRY

    CoSimIO::ExportData(rInfo, rValues);

    KRATOS_CATCH("")
}

CoSimIO::Info InfoFromParameters(Parameters rSettings)
{
    CoSimIO::Info info;

    for (auto it = rSettings.begin(); it != rSettings.end(); ++it) {
        if      (it->IsString()) info.Set<std::string>(it.name(), it->GetString());
        else if (it->IsInt())    info.Set<int>(it.name(),    it->GetInt());
        else if (it->IsBool())   info.Set<bool>(it.name(),   it->GetBool());
        else if (it->IsDouble()) info.Set<double>(it.name(), it->GetDouble());
        else KRATOS_WARNING("Kratos-CoSimIO") << "Setting with name \"" << it.name() << "\" cannot be converted to CoSimIO::Info and is ignored!" << std::endl;
    }

    return info;
}

} // CoSimIO_Wrappers namespace

void  AddCoSimIOToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto m_co_sim_io = m.def_submodule("CoSimIO");

    py::class_<CoSimIO::Info>(m_co_sim_io,"Info")
        .def(py::init<>())
        .def("Has",       &CoSimIO::Info::Has)
        .def("GetInt",    &CoSimIO::Info::Get<int>)
        .def("GetDouble", &CoSimIO::Info::Get<double>)
        .def("GetBool",   &CoSimIO::Info::Get<bool>)
        .def("GetString", &CoSimIO::Info::Get<std::string>)
        .def("SetInt",    &CoSimIO::Info::Set<int>)
        .def("SetDouble", &CoSimIO::Info::Set<double>)
        .def("SetBool",   &CoSimIO::Info::Set<bool>)
        .def("SetString", &CoSimIO::Info::Set<std::string>)
        .def("Erase",     &CoSimIO::Info::Erase)
        .def("Clear",     &CoSimIO::Info::Clear)
        .def("Size",      &CoSimIO::Info::Size)
        .def("__str__",   PrintObject<CoSimIO::Info>);
        ;

    m_co_sim_io.def("Connect", &CoSimIO::Connect);
    m_co_sim_io.def("Disconnect", &CoSimIO::Disconnect);

    m_co_sim_io.def("IsConverged", &CoSimIO::IsConverged);

    m_co_sim_io.def("SendControlSignal", CoSimIO::Internals::SendControlSignal);

    m_co_sim_io.def("ImportMesh", CoSimIO_Wrappers::ImportMesh);
    m_co_sim_io.def("ExportMesh", CoSimIO_Wrappers::ExportMesh);

    m_co_sim_io.def("ImportData", CoSimIO_Wrappers::ImportData_ModelPart_Scalar);
    m_co_sim_io.def("ExportData", CoSimIO_Wrappers::ExportData_ModelPart_Scalar);
    m_co_sim_io.def("ImportData", CoSimIO_Wrappers::ImportData_ModelPart_Vector);
    m_co_sim_io.def("ExportData", CoSimIO_Wrappers::ExportData_ModelPart_Vector);
    m_co_sim_io.def("ImportData", CoSimIO_Wrappers::ImportData_RawValues);
    m_co_sim_io.def("ExportData", CoSimIO_Wrappers::ExportData_RawValues);


    m_co_sim_io.def("InfoFromParameters", CoSimIO_Wrappers::InfoFromParameters);

    // // m_co_sim_io.def("ImportGeometry", CoSimIO_Wrappers::ImportGeometry); // This is not yet implemented in the CoSimIO
    // // m_co_sim_io.def("ExportGeometry", CoSimIO_Wrappers::ExportGeometry); // This is not yet implemented in the CoSimIO

    py::enum_<DataLocation>(m_co_sim_io,"DataLocation")
        .value("NodeHistorical",    DataLocation::NodeHistorical)
        .value("NodeNonHistorical", DataLocation::NodeNonHistorical)
        .value("Element",           DataLocation::Element)
        .value("Condition",         DataLocation::Condition)
        .value("ModelPart",         DataLocation::ModelPart)
        ;

    py::enum_<CoSimIO::ConnectionStatus>(m_co_sim_io,"ConnectionStatus")
        .value("NotConnected",CoSimIO::ConnectionStatus::NotConnected)
        .value("Connected",CoSimIO::ConnectionStatus::Connected)
        .value("Disconnected", CoSimIO::ConnectionStatus::Disconnected)
        .value("ConnectionError", CoSimIO::ConnectionStatus::ConnectionError)
        .value("DisconnectionError", CoSimIO::ConnectionStatus::DisconnectionError)
        ;

    py::enum_<CoSimIO::ControlSignal>(m_co_sim_io,"ControlSignal")
        .value("Dummy", CoSimIO::ControlSignal::Dummy)
        // .value("BreakSolutionLoop", CoSimIO::ControlSignal::BreakSolutionLoop)
        .value("ConvergenceAchieved", CoSimIO::ControlSignal::ConvergenceAchieved)

    //     .value("AdvanceInTime", CoSimIO::ControlSignal::AdvanceInTime)
    //     .value("InitializeSolutionStep", CoSimIO::ControlSignal::InitializeSolutionStep)
    //     .value("SolveSolutionStep", CoSimIO::ControlSignal::SolveSolutionStep)
    //     .value("FinalizeSolutionStep", CoSimIO::ControlSignal::FinalizeSolutionStep)

    //     .value("ImportGeometry", CoSimIO::ControlSignal::ImportGeometry)
    //     .value("ExportGeometry", CoSimIO::ControlSignal::ExportGeometry)
    //     .value("ImportMesh", CoSimIO::ControlSignal::ImportMesh)
    //     .value("ExportMesh", CoSimIO::ControlSignal::ExportMesh)
    //     .value("ImportData", CoSimIO::ControlSignal::ImportData)
    //     .value("ExportData", CoSimIO::ControlSignal::ExportData)
        ;
}

}  // namespace Python.
} // Namespace Kratos
