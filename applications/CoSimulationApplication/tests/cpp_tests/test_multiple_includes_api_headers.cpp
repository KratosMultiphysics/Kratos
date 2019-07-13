// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher
//

/*
This file includes the headers that are also included in other codes.
The reason is to test if the corresponding header can be included in multiple files without causing problems
(e.g. causing linking-erros if free functions are not static)
The "other" includ is supposed to be in the "add_xxx_to_python"

Furthermore it is tested if the header-file includes all necessary includes to compile standalone.
*/

// API includes
#include "custom_io/co_sim_EMPIRE_API.h"

namespace {
// defining aliases to the functions to avoid the "unused-functions" compiler-warning
const auto alias_EMPIRE_API_Connect = EMPIRE_API_Connect;
const auto alias_EMPIRE_API_Disconnect = EMPIRE_API_Disconnect;
const auto alias_EMPIRE_API_getUserDefinedText = EMPIRE_API_getUserDefinedText;
const auto alias_EMPIRE_API_sendMesh = EMPIRE_API_sendMesh;
const auto alias_EMPIRE_API_recvMesh = EMPIRE_API_recvMesh;
const auto alias_EMPIRE_API_sendDataField = EMPIRE_API_sendDataField;
const auto alias_EMPIRE_API_recvDataField = EMPIRE_API_recvDataField;
const auto alias_EMPIRE_API_sendSignal_double = EMPIRE_API_sendSignal_double;
const auto alias_EMPIRE_API_recvSignal_double = EMPIRE_API_recvSignal_double;
const auto alias_EMPIRE_API_EMPIRE_API_sendConvergenceSignal = EMPIRE_API_sendConvergenceSignal;
const auto alias_EMPIRE_API_EMPIRE_API_recvConvergenceSignal = EMPIRE_API_recvConvergenceSignal;
}