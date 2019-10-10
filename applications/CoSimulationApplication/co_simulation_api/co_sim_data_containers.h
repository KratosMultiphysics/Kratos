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
    struct Mesh
    {
        std::vector<double> node_coords;
        std::vector<int> connectivities;
        std::vector<int> cell_types;

        const std::string type="Mesh";
    };

    struct ConvergenceSignal
    {
        bool is_converged;

        const std::string type="ConvergenceSignal";
    };
}

} // namespace CoSim

#endif /* KRATOS_CO_SIM_DATA_CONTAINERS_H_INCLUDED */