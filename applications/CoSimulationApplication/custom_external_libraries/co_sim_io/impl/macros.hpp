//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#ifndef CO_SIM_IO_MACROS_H_INCLUDED
#define CO_SIM_IO_MACROS_H_INCLUDED

/* This file defines macros that are used inside the CoSimIO
Note that they are only defined here if they haven't been defined before.
This makes it possible to override them to use macros that are coming from
the code where the CoSimIO is included
*/

#ifndef CO_SIM_IO_ERROR
    #include <iostream>
    #include <stdexcept>
    struct err { // helper struct to mimic behavior of KRATOS_ERROR
    ~err() noexcept(false) { throw std::runtime_error("Error: "); } // destructors are noexcept by default
    };
    #define CO_SIM_IO_ERROR (err(), std::cout)
#endif

#ifndef CO_SIM_IO_ERROR_IF
    #define CO_SIM_IO_ERROR_IF(conditional) if (conditional) CO_SIM_IO_ERROR
#endif

#ifndef CO_SIM_IO_ERROR_IF_NOT
    #define CO_SIM_IO_ERROR_IF_NOT(conditional) if (!conditional) CO_SIM_IO_ERROR
#endif

#ifndef CO_SIM_IO_INFO
    #include <iostream>
    #define CO_SIM_IO_INFO(label) std::cout << label << ": "
#endif

#ifndef CO_SIM_IO_INFO_IF
    #define CO_SIM_IO_INFO_IF(label, conditional) if (conditional) CO_SIM_IO_INFO(label)
#endif

#endif // CO_SIM_IO_MACROS_H_INCLUDED