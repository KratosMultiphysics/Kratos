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

template<bool TIsDataField>
void Wrapper_SendArray(char* name, int sizeOfArray, std::vector<double> signal)
{
    // Wrapper is needed bcs pybind cannot do the conversion to raw-ptr automatically
    if (TIsDataField) {
        CoSimEMPIRE_API::EMPIRE_API_sendDataField(name, sizeOfArray, &signal[0]);
    } else {
        CoSimEMPIRE_API::EMPIRE_API_sendSignal_double(name, sizeOfArray, &signal[0]);
    }
}

template<bool TIsDataField>
void Wrapper_ReceiveArray(char* name, int sizeOfArray, pybind11::list signal)
{
    KRATOS_ERROR_IF(static_cast<int>(signal.size()) != sizeOfArray) << "The size of the list has to be specified before, expected size of " << sizeOfArray << ", current size: " << signal.size() << std::endl;

    // Wrapper is needed bcs pybind cannot do the conversion to raw-ptr automatically
    // also the list can only be modified in place otherwise the references are not working
    std::vector<double> vec_signal(sizeOfArray);
    if (TIsDataField) {
        CoSimEMPIRE_API::EMPIRE_API_recvDataField(name, sizeOfArray, &vec_signal[0]);
    } else {
        CoSimEMPIRE_API::EMPIRE_API_recvSignal_double(name, sizeOfArray, &vec_signal[0]);
    }

    // copy back the received values
    for (int i=0; i<sizeOfArray; ++i) {
        signal[i] = vec_signal[i];
    }
}

void  AddCustomIOToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto mEMPIREAPI = m.def_submodule("EMPIRE_API");

    mEMPIREAPI.def("EMPIRE_API_Connect", CoSimEMPIRE_API::EMPIRE_API_Connect);
    mEMPIREAPI.def("EMPIRE_API_Disconnect", CoSimEMPIRE_API::EMPIRE_API_Disconnect);

    mEMPIREAPI.def("EMPIRE_API_getUserDefinedText", CoSimEMPIRE_API::EMPIRE_API_getUserDefinedText);

    mEMPIREAPI.def("EMPIRE_API_sendMesh", CoSimEMPIRE_API::EMPIRE_API_sendMesh);
    // mEMPIREAPI.def("EMPIRE_API_recvMesh", CoSimEMPIRE_API::EMPIRE_API_recvMesh); // TODO check how to handle double**

    mEMPIREAPI.def("EMPIRE_API_sendDataField", Wrapper_SendArray<true>);
    mEMPIREAPI.def("EMPIRE_API_recvDataField", Wrapper_ReceiveArray<true>);

    mEMPIREAPI.def("EMPIRE_API_sendSignal_double", Wrapper_SendArray<false>);
    mEMPIREAPI.def("EMPIRE_API_recvSignal_double", Wrapper_ReceiveArray<false>);

    mEMPIREAPI.def("EMPIRE_API_recvConvergenceSignal", CoSimEMPIRE_API::EMPIRE_API_recvConvergenceSignal);
    mEMPIREAPI.def("EMPIRE_API_sendConvergenceSignal", CoSimEMPIRE_API::EMPIRE_API_sendConvergenceSignal);
}

}  // namespace Python.
} // Namespace Kratos

