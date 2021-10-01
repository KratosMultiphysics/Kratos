// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

/*
This file includes the headers that are also included in other codes.
The reason is to test if the corresponding header can be included in multiple files without causing problems
(e.g. causing linking-erros if free functions are not static)
The "other" includ is supposed to be in the "add_io_to_python"

Furthermore it is tested if the header-file includes all necessary includes to compile standalone.
*/

// API includes
#include "custom_io/co_sim_EMPIRE_API.h"
