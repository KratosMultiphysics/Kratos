/*
 * Author: Miguel Angel Celigueta
 *
 *  maceli@cimne.upc.edu
 */

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "dem_structures_coupling_application.h"

namespace Kratos {

    KratosDemStructuresCouplingApplication::KratosDemStructuresCouplingApplication() : KratosApplication("DemStructuresCouplingApplication"){}

    void KratosDemStructuresCouplingApplication::Register() {
        // Calling base class register to register Kratos components

        KratosApplication::Register();

        KRATOS_INFO("Dem-Struct") << std::endl;
        KRATOS_INFO("Dem-Struct") << "     KRATOS DEM STRUCTURES COUPLING APPLICATION " << std::endl;
        KRATOS_INFO("Dem-Struct") << std::endl;
        KRATOS_INFO("Dem-Struct") << "Importing DemStructuresCouplingApplication... ";
        KRATOS_INFO("") << " done." << std::endl;
    }
}  // namespace Kratos
