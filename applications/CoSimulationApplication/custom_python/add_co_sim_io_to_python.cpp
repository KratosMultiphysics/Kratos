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
#include <unordered_map>
#include <map>
#include <limits>

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
#include "co_simulation_io/co_sim_io.h"

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


enum class DataLocation { NodeHistorical, NodeNonHistorical, Element, Condition, ModelPart };

void ExportMesh(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    const ModelPart& rModelPart)
{
    // extract information from ModelPart
    const int num_nodes = rModelPart.NumberOfNodes();
    const int num_elems = /*(UseConditions) ? (rModelPart.NumberOfConditions()) :*/ (rModelPart.NumberOfElements());

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
        rConnectionName,
        rIdentifier,
        DataBuffers::vector_doubles,
        DataBuffers::vector_connectivities,
        DataBuffers::vector_types);
}

void ImportMesh(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    ModelPart& rModelPart)
{
    CoSimIO::ImportMesh(
        rConnectionName,
        rIdentifier,
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
    // const std::unordered_map<int, std::string> condition_name_map = {
    //     {1 , "PointCondition3D1N"},
    //     {2 , "LineCondition3D2N"},
    //     {3 , "SurfaceCondition3D3N"},
    //     {4 , "SurfaceCondition3D4N"}
    // };

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
        // if (UseConditions) {
        //     rModelPart.CreateNewCondition(condition_name_map.at(num_nodes_elem), i+1, elem_node_ids, p_props);
        // } else {
            rModelPart.CreateNewElement(vtk_type_info.second, i+1, elem_node_ids, p_props);
        // }
    }
}

void ExportData_Scalar(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    const ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const DataLocation DataLoc)
{
    // TODO resize only if too small
    if (DataLoc == DataLocation::NodeHistorical) {
        std::size_t counter = 0;
        DataBuffers::vector_doubles.resize(rModelPart.NumberOfNodes());
        for (const auto& r_node : rModelPart.Nodes()) {
            DataBuffers::vector_doubles[counter++] = r_node.FastGetSolutionStepValue(rVariable);
        }

    } else if (DataLoc == DataLocation::NodeNonHistorical) {
        std::size_t counter = 0;
        DataBuffers::vector_doubles.resize(rModelPart.NumberOfNodes());
        for (const auto& r_node : rModelPart.Nodes()) {
            DataBuffers::vector_doubles[counter++] = r_node.GetValue(rVariable);
        }

    } else if (DataLoc == DataLocation::Element) {
        std::size_t counter = 0;
        DataBuffers::vector_doubles.resize(rModelPart.NumberOfElements());
        for (const auto& r_elem : rModelPart.Elements()) {
            DataBuffers::vector_doubles[counter++] = r_elem.GetValue(rVariable);
        }

    } else if (DataLoc == DataLocation::Condition) {
        std::size_t counter = 0;
        DataBuffers::vector_doubles.resize(rModelPart.NumberOfConditions());
        for (const auto& r_cond : rModelPart.Conditions()) {
            DataBuffers::vector_doubles[counter++] = r_cond.GetValue(rVariable);
        }

    } else if (DataLoc == DataLocation::ModelPart) {
        DataBuffers::vector_doubles.resize(1);
        DataBuffers::vector_doubles[0] = rModelPart[rVariable];
    }

    CoSimIO::ExportData(
        rConnectionName,
        rIdentifier,
        DataBuffers::vector_doubles);
}

void ImportData_Scalar(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const DataLocation DataLoc)
{
    CoSimIO::ImportData(
        rConnectionName,
        rIdentifier,
        DataBuffers::vector_doubles);

    // TODO implement size-checks

    if (DataLoc == DataLocation::NodeHistorical) {
        std::size_t counter = 0;
        for (auto& r_node : rModelPart.Nodes()) {
            r_node.FastGetSolutionStepValue(rVariable) = DataBuffers::vector_doubles[counter++];
        }

    } else if (DataLoc == DataLocation::NodeNonHistorical) {
        std::size_t counter = 0;
        for (auto& r_node : rModelPart.Nodes()) {
            r_node.GetValue(rVariable) = DataBuffers::vector_doubles[counter++];
        }

    } else if (DataLoc == DataLocation::Element) {
        std::size_t counter = 0;
        for (auto& r_elem : rModelPart.Elements()) {
            r_elem.GetValue(rVariable) = DataBuffers::vector_doubles[counter++];
        }

    } else if (DataLoc == DataLocation::Condition) {
        std::size_t counter = 0;
        for (auto& r_cond : rModelPart.Conditions()) {
            r_cond.GetValue(rVariable) = DataBuffers::vector_doubles[counter++];
        }

    } else if (DataLoc == DataLocation::ModelPart) {
        rModelPart[rVariable] = DataBuffers::vector_doubles[0];
    }
}

void ExportData_Vector(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    const ModelPart& rModelPart,
    const Variable< array_1d<double, 3> >& rVariable,
    const DataLocation DataLoc)
{
    if (DataLoc == DataLocation::NodeHistorical) {
        std::size_t counter = 0;
        DataBuffers::vector_doubles.resize(rModelPart.NumberOfNodes()*3);
        for (const auto& r_node : rModelPart.Nodes()) {
            const array_1d<double, 3>& r_val = r_node.FastGetSolutionStepValue(rVariable);
            DataBuffers::vector_doubles[counter++] = r_val[0];
            DataBuffers::vector_doubles[counter++] = r_val[1];
            DataBuffers::vector_doubles[counter++] = r_val[2];
        }

    } else if (DataLoc == DataLocation::NodeNonHistorical) {
        std::size_t counter = 0;
        DataBuffers::vector_doubles.resize(rModelPart.NumberOfNodes()*3);
        for (const auto& r_node : rModelPart.Nodes()) {
            const array_1d<double, 3>& r_val = r_node.GetValue(rVariable);
            DataBuffers::vector_doubles[counter++] = r_val[0];
            DataBuffers::vector_doubles[counter++] = r_val[1];
            DataBuffers::vector_doubles[counter++] = r_val[2];
        }

    } else if (DataLoc == DataLocation::Element) {
        std::size_t counter = 0;
        DataBuffers::vector_doubles.resize(rModelPart.NumberOfElements()*3);
        for (const auto& r_elem : rModelPart.Elements()) {
            const array_1d<double, 3>& r_val = r_elem.GetValue(rVariable);
            DataBuffers::vector_doubles[counter++] = r_val[0];
            DataBuffers::vector_doubles[counter++] = r_val[1];
            DataBuffers::vector_doubles[counter++] = r_val[2];
        }

    } else if (DataLoc == DataLocation::Condition) {
        std::size_t counter = 0;
        DataBuffers::vector_doubles.resize(rModelPart.NumberOfConditions()*3);
        for (const auto& r_cond : rModelPart.Conditions()) {
            const array_1d<double, 3>& r_val = r_cond.GetValue(rVariable);
            DataBuffers::vector_doubles[counter++] = r_val[0];
            DataBuffers::vector_doubles[counter++] = r_val[1];
            DataBuffers::vector_doubles[counter++] = r_val[2];
        }

    } else if (DataLoc == DataLocation::ModelPart) {
        DataBuffers::vector_doubles.resize(3);
        const auto& r_val = rModelPart[rVariable];
        DataBuffers::vector_doubles[0] = r_val[0];
        DataBuffers::vector_doubles[1] = r_val[1];
        DataBuffers::vector_doubles[2] = r_val[2];
    }

    CoSimIO::ExportData(
        rConnectionName,
        rIdentifier,
        DataBuffers::vector_doubles);
}

void ImportData_Vector(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    ModelPart& rModelPart,
    const Variable< array_1d<double, 3> >& rVariable,
    const DataLocation DataLoc)
{
    CoSimIO::ImportData(
        rConnectionName,
        rIdentifier,
        DataBuffers::vector_doubles);

    // TODO implement size-checks

    if (DataLoc == DataLocation::NodeHistorical) {
        std::size_t counter = 0;
        for (auto& r_node : rModelPart.Nodes()) {
            array_1d<double, 3>& r_val = r_node.FastGetSolutionStepValue(rVariable);
            r_val[0] = DataBuffers::vector_doubles[counter++];
            r_val[1] = DataBuffers::vector_doubles[counter++];
            r_val[2] = DataBuffers::vector_doubles[counter++];
        }

    } else if (DataLoc == DataLocation::NodeNonHistorical) {
        std::size_t counter = 0;
        for (auto& r_node : rModelPart.Nodes()) {
            array_1d<double, 3>& r_val = r_node.GetValue(rVariable);
            r_val[0] = DataBuffers::vector_doubles[counter++];
            r_val[1] = DataBuffers::vector_doubles[counter++];
            r_val[2] = DataBuffers::vector_doubles[counter++];
        }

    } else if (DataLoc == DataLocation::Element) {
        std::size_t counter = 0;
        for (auto& r_elem : rModelPart.Elements()) {
            array_1d<double, 3>& r_val = r_elem.GetValue(rVariable);
            r_val[0] = DataBuffers::vector_doubles[counter++];
            r_val[1] = DataBuffers::vector_doubles[counter++];
            r_val[2] = DataBuffers::vector_doubles[counter++];
        }

    } else if (DataLoc == DataLocation::Condition) {
        std::size_t counter = 0;
        for (auto& r_cond : rModelPart.Conditions()) {
            array_1d<double, 3>& r_val = r_cond.GetValue(rVariable);
            r_val[0] = DataBuffers::vector_doubles[counter++];
            r_val[1] = DataBuffers::vector_doubles[counter++];
            r_val[2] = DataBuffers::vector_doubles[counter++];
        }

    } else if (DataLoc == DataLocation::ModelPart) {
        auto& r_val = rModelPart[rVariable];
        r_val[0] = DataBuffers::vector_doubles[0];
        r_val[1] = DataBuffers::vector_doubles[1];
        r_val[2] = DataBuffers::vector_doubles[2];
    }
}

std::vector<double> ImportData_Values(
    const std::string& rConnectionName,
    const std::string& rIdentifier)
{
    CoSimIO::ImportData(
        rConnectionName,
        rIdentifier,
        DataBuffers::vector_doubles);

    return DataBuffers::vector_doubles;
}

void ExportData_Values(
    const std::string& rConnectionName,
    const std::string& rIdentifier,
    std::vector<double>& rValues)
{
    CoSimIO::ExportData(
        rConnectionName,
        rIdentifier,
        rValues);
}

bool IsConverged(const std::string& rConnectionName)
{
    int convergence_signal;
    CoSimIO::IsConverged(rConnectionName, convergence_signal);
    return convergence_signal;
}

} // helpers namespace

void  AddCoSimIOToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto mCoSimIO = m.def_submodule("CoSimIO");

    void (*ConnectWithSettings)(const std::string&, CoSimIO::SettingsType) = &CoSimIO::Connect;

    mCoSimIO.def("Connect", ConnectWithSettings);

    mCoSimIO.def("Disconnect", CoSimIO::Disconnect);

    mCoSimIO.def("IsConverged", CoSimIO_Wrappers::IsConverged);

    mCoSimIO.def("SendControlSignal", CoSimIO::SendControlSignal);

    mCoSimIO.def("ImportMesh", CoSimIO_Wrappers::ImportMesh);
    mCoSimIO.def("ExportMesh", CoSimIO_Wrappers::ExportMesh);

    mCoSimIO.def("ImportData", CoSimIO_Wrappers::ImportData_Scalar);
    mCoSimIO.def("ExportData", CoSimIO_Wrappers::ExportData_Scalar);
    mCoSimIO.def("ImportData", CoSimIO_Wrappers::ImportData_Vector);
    mCoSimIO.def("ExportData", CoSimIO_Wrappers::ExportData_Vector);
    mCoSimIO.def("ImportData", CoSimIO_Wrappers::ImportData_Values);
    mCoSimIO.def("ExportData", CoSimIO_Wrappers::ExportData_Values);

    // mCoSimIO.def("ImportGeometry", CoSimIO_Wrappers::ImportGeometry); // This is not yet implemented in the CoSimIO
    // mCoSimIO.def("ExportGeometry", CoSimIO_Wrappers::ExportGeometry); // This is not yet implemented in the CoSimIO

    py::enum_<CoSimIO_Wrappers::DataLocation>(mCoSimIO,"DataLocation")
        .value("NodeHistorical",    CoSimIO_Wrappers::DataLocation::NodeHistorical)
        .value("NodeNonHistorical", CoSimIO_Wrappers::DataLocation::NodeNonHistorical)
        .value("Element",           CoSimIO_Wrappers::DataLocation::Element)
        .value("Condition",         CoSimIO_Wrappers::DataLocation::Condition)
        .value("ModelPart",         CoSimIO_Wrappers::DataLocation::ModelPart)
        ;

    py::enum_<CoSimIO::ControlSignal>(mCoSimIO,"ControlSignal")
        .value("Dummy", CoSimIO::ControlSignal::Dummy)
        .value("BreakSolutionLoop", CoSimIO::ControlSignal::BreakSolutionLoop)
        .value("ConvergenceAchieved", CoSimIO::ControlSignal::ConvergenceAchieved)

        .value("AdvanceInTime", CoSimIO::ControlSignal::AdvanceInTime)
        .value("InitializeSolutionStep", CoSimIO::ControlSignal::InitializeSolutionStep)
        .value("SolveSolutionStep", CoSimIO::ControlSignal::SolveSolutionStep)
        .value("FinalizeSolutionStep", CoSimIO::ControlSignal::FinalizeSolutionStep)

        .value("ImportGeometry", CoSimIO::ControlSignal::ImportGeometry)
        .value("ExportGeometry", CoSimIO::ControlSignal::ExportGeometry)
        .value("ImportMesh", CoSimIO::ControlSignal::ImportMesh)
        .value("ExportMesh", CoSimIO::ControlSignal::ExportMesh)
        .value("ImportData", CoSimIO::ControlSignal::ImportData)
        .value("ExportData", CoSimIO::ControlSignal::ExportData)
        ;
}

}  // namespace Python.
} // Namespace Kratos

