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
#include "custom_python/add_co_sim_io_to_python.h"

// CoSimIO
#include "co_simulation_api/co_sim_io.h"

namespace Kratos {
namespace Python {

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

        .def("ImportGeometry",&CoSimIOType::Import<CoSim::DataContainers::Geometry>)
        .def("ExportGeometry",&CoSimIOType::Export<CoSim::DataContainers::Geometry>)
        ;






}

}  // namespace Python.
} // Namespace Kratos

