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
#include <pybind11/functional.h>

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

#include "custom_external_libraries/co_sim_io/python/connection_status_to_python.hpp"
#include "custom_external_libraries/co_sim_io/python/info_to_python.hpp"
#include "custom_external_libraries/co_sim_io/python/vector_to_python.hpp"
#include "custom_external_libraries/co_sim_io/python/version_to_python.hpp"

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
};

// declaring the static members
std::vector<double> DataBuffers::vector_doubles;

const std::map<GeometryData::KratosGeometryType, CoSimIO::ElementType> elem_type_map {
    {GeometryData::KratosGeometryType::Kratos_Hexahedra3D20,    CoSimIO::ElementType::Hexahedra3D20   },
    {GeometryData::KratosGeometryType::Kratos_Hexahedra3D27,    CoSimIO::ElementType::Hexahedra3D27   },
    {GeometryData::KratosGeometryType::Kratos_Hexahedra3D8,     CoSimIO::ElementType::Hexahedra3D8     },
    {GeometryData::KratosGeometryType::Kratos_Prism3D15,        CoSimIO::ElementType::Prism3D15   },
    {GeometryData::KratosGeometryType::Kratos_Prism3D6,         CoSimIO::ElementType::Prism3D6     },
    {GeometryData::KratosGeometryType::Kratos_Quadrilateral2D4, CoSimIO::ElementType::Quadrilateral2D4     },
    {GeometryData::KratosGeometryType::Kratos_Quadrilateral2D8, CoSimIO::ElementType::Quadrilateral2D8     },
    {GeometryData::KratosGeometryType::Kratos_Quadrilateral2D9, CoSimIO::ElementType::Quadrilateral2D9     },
    {GeometryData::KratosGeometryType::Kratos_Quadrilateral3D4, CoSimIO::ElementType::Quadrilateral3D4     },
    {GeometryData::KratosGeometryType::Kratos_Quadrilateral3D8, CoSimIO::ElementType::Quadrilateral3D8     },
    {GeometryData::KratosGeometryType::Kratos_Quadrilateral3D9, CoSimIO::ElementType::Quadrilateral3D9     },
    {GeometryData::KratosGeometryType::Kratos_Tetrahedra3D10,   CoSimIO::ElementType::Tetrahedra3D10     },
    {GeometryData::KratosGeometryType::Kratos_Tetrahedra3D4,    CoSimIO::ElementType::Tetrahedra3D4   },
    {GeometryData::KratosGeometryType::Kratos_Triangle2D3,      CoSimIO::ElementType::Triangle2D3   },
    {GeometryData::KratosGeometryType::Kratos_Triangle2D6,      CoSimIO::ElementType::Triangle2D6   },
    {GeometryData::KratosGeometryType::Kratos_Triangle3D3,      CoSimIO::ElementType::Triangle3D3   },
    {GeometryData::KratosGeometryType::Kratos_Triangle3D6,      CoSimIO::ElementType::Triangle3D6   },
    {GeometryData::KratosGeometryType::Kratos_Line2D2,          CoSimIO::ElementType::Line2D2   },
    {GeometryData::KratosGeometryType::Kratos_Line2D3,          CoSimIO::ElementType::Line2D3   },
    {GeometryData::KratosGeometryType::Kratos_Line3D2,          CoSimIO::ElementType::Line3D2   },
    {GeometryData::KratosGeometryType::Kratos_Line3D3,          CoSimIO::ElementType::Line3D3   },
    {GeometryData::KratosGeometryType::Kratos_Point2D,          CoSimIO::ElementType::Point2D   },
    {GeometryData::KratosGeometryType::Kratos_Point3D,          CoSimIO::ElementType::Point3D   }
};

const std::map<CoSimIO::ElementType, std::string> elem_name_map {
    {CoSimIO::ElementType::Hexahedra3D20, "Element3D20N"},
    {CoSimIO::ElementType::Hexahedra3D27, "Element3D27N"},
    {CoSimIO::ElementType::Hexahedra3D8, "Element3D8N"},
    {CoSimIO::ElementType::Prism3D15, "Element3D15N"},
    {CoSimIO::ElementType::Prism3D6, "Element3D6N"},
    {CoSimIO::ElementType::Quadrilateral2D4, "Element2D4N"},
    {CoSimIO::ElementType::Quadrilateral2D8, "Element2D8N"},
    {CoSimIO::ElementType::Quadrilateral2D9, "Element2D9N"},
    {CoSimIO::ElementType::Quadrilateral3D8, "Element3D8N"},
    {CoSimIO::ElementType::Tetrahedra3D10, "Element3D10N"},
    {CoSimIO::ElementType::Tetrahedra3D4, "Element3D4N"},
    {CoSimIO::ElementType::Triangle2D3, "Element2D3N"},
    {CoSimIO::ElementType::Triangle2D6, "Element2D6N"},
    {CoSimIO::ElementType::Triangle3D3, "Element3D3N"},
    {CoSimIO::ElementType::Line2D2, "Element2D2N"},
    {CoSimIO::ElementType::Line3D2, "Element3D2N"},
    {CoSimIO::ElementType::Point2D, "Element2D1N"},
    {CoSimIO::ElementType::Point3D, "Element3D1N"}
};

void ExportMesh(
    CoSimIO::Info& rInfo,
    const ModelPart& rModelPart)
{
    KRATOS_TRY

    CoSimIO::ModelPart co_sim_io_model_part(rModelPart.Name());

    for (const auto& r_node : rModelPart.Nodes()) {
        co_sim_io_model_part.CreateNewNode(
            r_node.Id(),
            // TODO: use initial or current coordinates?
            r_node.X0(),
            r_node.Y0(),
            r_node.Z0()
        );
    };

    CoSimIO::ConnectivitiesType conn;
    for (const auto& r_elem : rModelPart.Elements()) {
        const auto& r_geom = r_elem.GetGeometry();
        if (conn.size() != r_geom.PointsNumber()) {
            conn.resize(r_geom.PointsNumber());
        };

        for (std::size_t i=0; i<r_geom.PointsNumber(); ++i) {
            conn[i] = r_geom[i].Id();
        }

        auto elem_type_it = elem_type_map.find(r_geom.GetGeometryType());
        KRATOS_ERROR_IF(elem_type_it == elem_type_map.end()) << "No CoSimIO element type found for this Kratos element type (" << static_cast<int>(r_geom.GetGeometryType()) << ")!" << std::endl;

        co_sim_io_model_part.CreateNewElement(
            r_elem.Id(),
            elem_type_it->second,
            conn
        );
    };

    CoSimIO::ExportMesh(
        rInfo,
        co_sim_io_model_part);

    KRATOS_CATCH("")
}

void ImportMesh(
    CoSimIO::Info& rInfo,
    ModelPart& rModelPart)
{
    KRATOS_TRY

    CoSimIO::ModelPart co_sim_io_model_part(rModelPart.Name());

    CoSimIO::ImportMesh(
        rInfo,
        co_sim_io_model_part);

    // fill ModelPart from received Mesh
    KRATOS_ERROR_IF(rModelPart.NumberOfNodes() > 0) << "ModelPart is not empty, it has nodes!" << std::endl;
    KRATOS_ERROR_IF(rModelPart.NumberOfProperties() > 0) << "ModelPart is not empty, it has properties!" << std::endl;
    KRATOS_ERROR_IF(rModelPart.IsDistributed()) << "ModelPart cannot be distributed!" << std::endl;

    // fill ModelPart with received entities
    for (auto node_it=co_sim_io_model_part.NodesBegin(); node_it!=co_sim_io_model_part.NodesEnd(); ++node_it) {
        rModelPart.CreateNewNode(
            (*node_it)->Id(),
            (*node_it)->X(),
            (*node_it)->Y(),
            (*node_it)->Z()
        );
    };

    auto p_props = rModelPart.CreateNewProperties(0);

    std::vector<IndexType> conn;
    for (auto elem_it=co_sim_io_model_part.ElementsBegin(); elem_it!=co_sim_io_model_part.ElementsEnd(); ++elem_it) {
        if (conn.size() != (*elem_it)->NumberOfNodes()) {
            conn.resize((*elem_it)->NumberOfNodes());
        };

        const auto nodes_begin = (*elem_it)->NodesBegin();
        for (std::size_t i=0; i<(*elem_it)->NumberOfNodes(); ++i) {
            conn[i] = (*(nodes_begin+i))->Id();
        };

        auto elem_name_it = elem_name_map.find((*elem_it)->Type());
        if (elem_name_it == elem_name_map.end()) {
            std::stringstream err;
            err << "No Kratos element found for this element type (" << static_cast<int>((*elem_it)->Type()) << ")!\nOnly the following types are available:";
            for (const auto& r_type_name_pair : elem_name_map) {
                err << "\n\t" << r_type_name_pair.second;
            }
            KRATOS_ERROR << err.str();
        }

        rModelPart.CreateNewElement(
            elem_name_it->second,
            (*elem_it)->Id(),
            conn,
            p_props
        );
    };

    KRATOS_CATCH("")
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

void ImportData_RawValues(
    CoSimIO::Info& rInfo,
    CoSimIO::VectorWrapper<double>& rValues)
{
    KRATOS_TRY

    CoSimIO::ImportData(rInfo, rValues.Vector());

    KRATOS_CATCH("")
}

void ExportData_RawValues(
    CoSimIO::Info& rInfo,
    const CoSimIO::VectorWrapper<double>& rValues)
{
    KRATOS_TRY

    CoSimIO::ExportData(rInfo, rValues.Vector());

    KRATOS_CATCH("")
}

CoSimIO::Info InfoFromParameters(Parameters rSettings)
{
    KRATOS_TRY

    CoSimIO::Info info;

    for (auto it = rSettings.begin(); it != rSettings.end(); ++it) {
        if      (it->IsString()) info.Set<std::string>(it.name(), it->GetString());
        else if (it->IsInt())    info.Set<int>(it.name(),    it->GetInt());
        else if (it->IsBool())   info.Set<bool>(it.name(),   it->GetBool());
        else if (it->IsDouble()) info.Set<double>(it.name(), it->GetDouble());
        else KRATOS_WARNING("Kratos-CoSimIO") << "Setting with name \"" << it.name() << "\" cannot be converted to CoSimIO::Info and is ignored!" << std::endl;
    }

    return info;

    KRATOS_CATCH("")
}

} // CoSimIO_Wrappers namespace

void  AddCoSimIOToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto m_co_sim_io = m.def_submodule("CoSimIO");

    CoSimIO::AddCoSimIOInfoToPython(m_co_sim_io);
    CoSimIO::AddCoSimIOConnectionStatusToPython(m_co_sim_io);
    CoSimIO::AddCoSimIOVectorToPython(m_co_sim_io);
    CoSimIO::AddCoSimIOVersionToPython(m_co_sim_io);

    m_co_sim_io.def("Connect",    &CoSimIO::Connect);
    m_co_sim_io.def("Disconnect", &CoSimIO::Disconnect);

    m_co_sim_io.def("ImportInfo", &CoSimIO::ImportInfo);
    m_co_sim_io.def("ExportInfo", &CoSimIO::ExportInfo);

    m_co_sim_io.def("ImportMesh", CoSimIO_Wrappers::ImportMesh);
    m_co_sim_io.def("ExportMesh", CoSimIO_Wrappers::ExportMesh);

    m_co_sim_io.def("ImportData", CoSimIO_Wrappers::ImportData_ModelPart_Scalar);
    m_co_sim_io.def("ExportData", CoSimIO_Wrappers::ExportData_ModelPart_Scalar);
    m_co_sim_io.def("ImportData", CoSimIO_Wrappers::ImportData_ModelPart_Vector);
    m_co_sim_io.def("ExportData", CoSimIO_Wrappers::ExportData_ModelPart_Vector);
    m_co_sim_io.def("ImportData", CoSimIO_Wrappers::ImportData_RawValues);
    m_co_sim_io.def("ExportData", CoSimIO_Wrappers::ExportData_RawValues);

    m_co_sim_io.def("Register", [](
        const CoSimIO::Info& I_Info,
        std::function<CoSimIO::Info(const CoSimIO::Info&)> FunctionPointer)
        { return CoSimIO::Register(I_Info, FunctionPointer); } );
    m_co_sim_io.def("Run",      &CoSimIO::Run);


    m_co_sim_io.def("InfoFromParameters", CoSimIO_Wrappers::InfoFromParameters);

    py::enum_<DataLocation>(m_co_sim_io,"DataLocation")
        .value("NodeHistorical",    DataLocation::NodeHistorical)
        .value("NodeNonHistorical", DataLocation::NodeNonHistorical)
        .value("Element",           DataLocation::Element)
        .value("Condition",         DataLocation::Condition)
        .value("ModelPart",         DataLocation::ModelPart)
        ;

}

}  // namespace Python.
} // Namespace Kratos
