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

#ifndef KRATOS_CO_SIM_DATA_CONTAINERS_H_INCLUDED
#define KRATOS_CO_SIM_DATA_CONTAINERS_H_INCLUDED

// System includes
#include <vector>
#include <string>

namespace CoSim {

namespace DataContainers {

    struct Geometry
    {
        const std::vector<double>& patches;
        const std::vector<double>& triming;
        /* ... */
    };

    struct Mesh
    {
        const std::vector<double>& node_coords;
        const std::vector<int>& connectivities;
        const std::vector<int>& cell_types;
    };

    struct Data
    {
        const std::vector<double>& data;
    };
}

} // namespace CoSim

#endif /* KRATOS_CO_SIM_DATA_CONTAINERS_H_INCLUDED */