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

// External includes

// Project includes
#include "includes/model_part.h"
#include "custom_python/add_co_sim_io_to_python.h"

// CoSimIO
#include "co_simulation_api/co_sim_io.h"

namespace Kratos {
namespace Python {

namespace CoSimIO_Wrappers { // helpers namespace

void ExportGeometry(const ModelPart& rModelPart)
{

}

void ImportGeometry(ModelPart& rModelPart)
{

}

void ExportMesh(const ModelPart& rModelPart)
{

}

void ImportMesh(ModelPart& rModelPart)
{

}

void ExportData(const ModelPart& rModelPart)
{

}

void ImportData(ModelPart& rModelPart)
{

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

