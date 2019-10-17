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

/* This file defines macros that are used inside the CoSimIO
Note that they are only defined here if they haven't been defined before.
This makes it possible to override them to use macros that are coming from
the code where the CoSimIO is included
*/

#ifndef KRATOS_CO_SIM_ERROR
#include <iostream>
#include <stdexcept>
struct err { // helper struct to mimic behavior of KRATOS_ERROR
   err() {std::cout << "Error: ";}
  ~err() { throw std::exception(); }
};
#define KRATOS_CO_SIM_ERROR (err(), std::cout)
#endif

#ifndef KRATOS_CO_SIM_INFO
#include <iostream>
#define CS_LOG std::cout << "[CoSimIO] "
#endif

#ifndef KRATOS_CO_SIM_INFO_IF
#define KRATOS_CO_SIM_INFO_IF if (conditional) KRATOS_CO_SIM_INFO
#endif
