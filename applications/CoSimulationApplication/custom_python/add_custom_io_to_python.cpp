// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aditya Ghantasala
//                   Philipp Bucher
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_io_to_python.h"

// IO
#include "custom_io/co_sim_EMPIRE_API.h"

namespace Kratos {
namespace Python {

void  AddCustomIOToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto mEMPIREAPI = m.def_submodule("EMPIRE_API");

    mEMPIREAPI.def("EMPIRE_API_Connect", CoSimEMPIRE_API::EMPIRE_API_Connect);
    mEMPIREAPI.def("EMPIRE_API_Disconnect", CoSimEMPIRE_API::EMPIRE_API_Disconnect);

    mEMPIREAPI.def("EMPIRE_API_getUserDefinedText", CoSimEMPIRE_API::EMPIRE_API_getUserDefinedText);

    mEMPIREAPI.def("EMPIRE_API_sendMesh", CoSimEMPIRE_API::EMPIRE_API_sendMesh);
    // mEMPIREAPI.def("EMPIRE_API_recvMesh", CoSimEMPIRE_API::EMPIRE_API_recvMesh); // TODO check how to handle double**

    mEMPIREAPI.def("EMPIRE_API_sendDataField", CoSimEMPIRE_API::EMPIRE_API_sendDataField);
    mEMPIREAPI.def("EMPIRE_API_recvDataField", CoSimEMPIRE_API::EMPIRE_API_recvDataField);

    mEMPIREAPI.def("EMPIRE_API_sendSignal_double", CoSimEMPIRE_API::EMPIRE_API_sendSignal_double);
    mEMPIREAPI.def("EMPIRE_API_recvSignal_double", CoSimEMPIRE_API::EMPIRE_API_recvSignal_double);

    mEMPIREAPI.def("EMPIRE_API_recvConvergenceSignal", CoSimEMPIRE_API::EMPIRE_API_recvConvergenceSignal);
    mEMPIREAPI.def("EMPIRE_API_sendConvergenceSignal", CoSimEMPIRE_API::EMPIRE_API_sendConvergenceSignal);
}

}  // namespace Python.
} // Namespace Kratos

