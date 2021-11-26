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
#include "custom_utilities/co_sim_io_conversion_utilities.h"
#include "co_simulation_application_variables.h"

// CoSimIO
#include "custom_external_libraries/CoSimIO/co_sim_io/co_sim_io.hpp"

#include "custom_external_libraries/CoSimIO/co_sim_io/python/connection_status_to_python.hpp"
#include "custom_external_libraries/CoSimIO/co_sim_io/python/info_to_python.hpp"
#include "custom_external_libraries/CoSimIO/co_sim_io/python/vector_to_python.hpp"
#include "custom_external_libraries/CoSimIO/co_sim_io/python/version_to_python.hpp"

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

struct Accessor_Hist_Get_Scalar
{
    static void Execute(const Node<3>& rNode, std::vector<double>& rData, const std::size_t Index, const Variable<double>& rVariable) {
        rData[Index] = rNode.FastGetSolutionStepValue(rVariable);
    }
};
struct Accessor_NonHist_Get_Scalar
{
    template<class TEntityType>
    static void Execute(const TEntityType& rEntity, std::vector<double>& rData, const std::size_t Index, const Variable<double>& rVariable) {
        rData[Index] = rEntity.GetValue(rVariable);
    }
};
struct Accessor_Hist_Set_Scalar
{
    static void Execute(Node<3>& rNode, const std::vector<double>& rData, const std::size_t Index, const Variable<double>& rVariable) {
        rNode.FastGetSolutionStepValue(rVariable) = rData[Index];
    }
};
struct Accessor_NonHist_Set_Scalar
{
    template<class TEntityType>
    static void Execute(TEntityType& rEntity, const std::vector<double>& rData, const std::size_t Index, const Variable<double>& rVariable) {
        rEntity.GetValue(rVariable) = rData[Index];
    }
};

struct Accessor_Hist_Get_Vector
{
    template<std::size_t TSize>
    static void Execute(const Node<3>& rNode, std::vector<double>& rData, const std::size_t Index, const Variable<array_1d<double, TSize>>& rVariable) {
        const array_1d<double, TSize>& var = rNode.FastGetSolutionStepValue(rVariable);
        for (std::size_t i=0; i<TSize; ++i) { rData[Index+i] = var[i]; }
    }
};

struct Accessor_NonHist_Get_Vector
{
    template<class TEntityType, std::size_t TSize>
    static void Execute(const TEntityType& rNode, std::vector<double>& rData, const std::size_t Index, const Variable<array_1d<double, TSize>>& rVariable) {
        const array_1d<double, TSize>& var = rNode.GetValue(rVariable);
        for (std::size_t i=0; i<TSize; ++i) { rData[Index+i] = var[i]; }
    }
};

struct Accessor_Hist_Set_Vector
{
    template<std::size_t TSize>
    static void Execute(Node<3>& rNode, const std::vector<double>& rData, const std::size_t Index, const Variable<array_1d<double, TSize>>& rVariable) {
        array_1d<double, TSize>& var = rNode.FastGetSolutionStepValue(rVariable);
        for (std::size_t i=0; i<TSize; ++i) { var[i] = rData[Index+i]; }
    }
};

struct Accessor_NonHist_Set_Vector
{
    template<class TEntityType, std::size_t TSize>
    static void Execute(TEntityType& rNode, const std::vector<double>& rData, const std::size_t Index, const Variable<array_1d<double, TSize>>& rVariable) {
        array_1d<double, TSize>& var = rNode.GetValue(rVariable);
        for (std::size_t i=0; i<TSize; ++i) { var[i] = rData[Index+i]; }
    }
};

template<class TAccessor, class TContainerType, class TVarDataType>
void AccessDataWithOrder(
    TContainerType& rContainer,
    const Variable<TVarDataType>& rVariable,
    const std::vector<std::size_t>& rOrder,
    std::vector<double>& rData)
{
    IndexPartition<std::size_t>(rContainer.size()).for_each(
        [&rContainer, &rVariable, &rOrder, &rData]
            (const std::size_t Index) {
                const std::size_t entity_id = rOrder[Index];
                auto& r_entity = *rContainer.find(entity_id);
                TAccessor::Execute(r_entity, rData, Index, rVariable);
    });
}

void ExportMesh(
    CoSimIO::Info& rInfo,
    const Kratos::ModelPart& rModelPart)
{
    KRATOS_TRY

    CoSimIO::ModelPart co_sim_io_model_part(rModelPart.Name());

    CoSimIOConversionUtilities::KratosModelPartToCoSimIOModelPart(rModelPart, co_sim_io_model_part);

    CoSimIO::ExportMesh(
        rInfo,
        co_sim_io_model_part);

    KRATOS_CATCH("")
}

void ImportMesh(
    CoSimIO::Info& rInfo,
    ModelPart& rModelPart,
    const DataCommunicator& rDataComm)
{
    KRATOS_TRY

    CoSimIO::ModelPart co_sim_io_model_part(rModelPart.Name());

    CoSimIO::ImportMesh(
        rInfo,
        co_sim_io_model_part);

    CoSimIOConversionUtilities::CoSimIOModelPartToKratosModelPart(co_sim_io_model_part, rModelPart, rDataComm);

    KRATOS_CATCH("")
}

void ExportData_ModelPart_Scalar(
    CoSimIO::Info& rInfo,
    ModelPart& rModelPart,
    const Variable<double>& rVariable,
    const DataLocation DataLoc)
{
    KRATOS_TRY

    if (DataLoc == DataLocation::NodeHistorical && rModelPart.Has(NODES_ID_INDEX_MAP)) {
        AccessDataWithOrder<Accessor_Hist_Get_Scalar>(rModelPart.Nodes(), rVariable, rModelPart[NODES_ID_INDEX_MAP], DataBuffers::vector_doubles);
    } else if (DataLoc == DataLocation::NodeNonHistorical && rModelPart.Has(NODES_ID_INDEX_MAP)) {
        AccessDataWithOrder<Accessor_NonHist_Get_Scalar>(rModelPart.Nodes(), rVariable, rModelPart[NODES_ID_INDEX_MAP], DataBuffers::vector_doubles);
    } else if (DataLoc == DataLocation::Element && rModelPart.Has(ELEMENTS_ID_INDEX_MAP)) {
        AccessDataWithOrder<Accessor_NonHist_Get_Scalar>(rModelPart.Elements(), rVariable, rModelPart[ELEMENTS_ID_INDEX_MAP], DataBuffers::vector_doubles);
    } else {
        AuxiliarModelPartUtilities(rModelPart).GetScalarData<double>(rVariable, DataLoc, DataBuffers::vector_doubles);
    }

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

    if (DataLoc == DataLocation::NodeHistorical && rModelPart.Has(NODES_ID_INDEX_MAP)) {
        AccessDataWithOrder<Accessor_Hist_Set_Scalar>(rModelPart.Nodes(), rVariable, rModelPart[NODES_ID_INDEX_MAP], DataBuffers::vector_doubles);
    } else if (DataLoc == DataLocation::NodeNonHistorical && rModelPart.Has(NODES_ID_INDEX_MAP)) {
        AccessDataWithOrder<Accessor_NonHist_Set_Scalar>(rModelPart.Nodes(), rVariable, rModelPart[NODES_ID_INDEX_MAP], DataBuffers::vector_doubles);
    } else if (DataLoc == DataLocation::Element && rModelPart.Has(ELEMENTS_ID_INDEX_MAP)) {
        AccessDataWithOrder<Accessor_NonHist_Set_Scalar>(rModelPart.Elements(), rVariable, rModelPart[ELEMENTS_ID_INDEX_MAP], DataBuffers::vector_doubles);
    } else {
        AuxiliarModelPartUtilities(rModelPart).SetScalarData<double>(rVariable, DataLoc, DataBuffers::vector_doubles);
    }

    KRATOS_CATCH("")
}

void ExportData_ModelPart_Vector(
    CoSimIO::Info& rInfo,
    ModelPart& rModelPart,
    const Variable< array_1d<double, 3> >& rVariable,
    const DataLocation DataLoc)
{

    KRATOS_TRY

    if (DataLoc == DataLocation::NodeHistorical && rModelPart.Has(NODES_ID_INDEX_MAP)) {
        AccessDataWithOrder<Accessor_Hist_Get_Vector>(rModelPart.Nodes(), rVariable, rModelPart[NODES_ID_INDEX_MAP], DataBuffers::vector_doubles);
    } else if (DataLoc == DataLocation::NodeNonHistorical && rModelPart.Has(NODES_ID_INDEX_MAP)) {
        AccessDataWithOrder<Accessor_NonHist_Get_Vector>(rModelPart.Nodes(), rVariable, rModelPart[NODES_ID_INDEX_MAP], DataBuffers::vector_doubles);
    } else if (DataLoc == DataLocation::Element && rModelPart.Has(ELEMENTS_ID_INDEX_MAP)) {
        AccessDataWithOrder<Accessor_NonHist_Get_Vector>(rModelPart.Elements(), rVariable, rModelPart[ELEMENTS_ID_INDEX_MAP], DataBuffers::vector_doubles);
    } else {
        AuxiliarModelPartUtilities(rModelPart).GetVectorData< array_1d<double, 3> >(rVariable, DataLoc, DataBuffers::vector_doubles);
    }

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

    if (DataLoc == DataLocation::NodeHistorical && rModelPart.Has(NODES_ID_INDEX_MAP)) {
        AccessDataWithOrder<Accessor_Hist_Set_Vector>(rModelPart.Nodes(), rVariable, rModelPart[NODES_ID_INDEX_MAP], DataBuffers::vector_doubles);
    } else if (DataLoc == DataLocation::NodeNonHistorical && rModelPart.Has(NODES_ID_INDEX_MAP)) {
        AccessDataWithOrder<Accessor_NonHist_Set_Vector>(rModelPart.Nodes(), rVariable, rModelPart[NODES_ID_INDEX_MAP], DataBuffers::vector_doubles);
    } else if (DataLoc == DataLocation::Element && rModelPart.Has(ELEMENTS_ID_INDEX_MAP)) {
        AccessDataWithOrder<Accessor_NonHist_Set_Vector>(rModelPart.Elements(), rVariable, rModelPart[ELEMENTS_ID_INDEX_MAP], DataBuffers::vector_doubles);
    } else {
        AuxiliarModelPartUtilities(rModelPart).SetVectorData< array_1d<double, 3> >(rVariable, DataLoc, DataBuffers::vector_doubles);
    }

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


    m_co_sim_io.def("InfoFromParameters", CoSimIOConversionUtilities::InfoFromParameters);

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
