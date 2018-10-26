//
// Last Modified by:    Salva Latorre
//

// Project includes
#include "dem_structures_coupling_application.h"
#include "includes/kernel.h"
#include "includes/kratos_flags.h"
#include "containers/flags.h"
#include "geometries/point_3d.h"
#include "geometries/line_3d_2.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/sphere_3d_1.h"
#include "utilities/quaternion.h"

namespace Kratos {

    KratosDemStructuresCouplingApplication::KratosDemStructuresCouplingApplication() : KratosApplication("DemStructuresCouplingApplication"){}

    void KratosDemStructuresCouplingApplication::Register() {
        // Calling base class register to register Kratos components

        KratosApplication::Register();

        KRATOS_INFO("Dem-Struct") << std::endl;
        KRATOS_INFO("Dem-Struct") << "     KRATOS DEM STRUCTURES COUPLING APPLICATION " << std::endl;
        KRATOS_INFO("Dem-Struct") << std::endl;
        KRATOS_INFO("Dem-Struct") << "Importing DemStructuresCouplingApplication... ";

        KRATOS_INFO("") << "( compiled in mode \"" << Kernel::BuildType() << "\" )";

        KRATOS_INFO("") << " done." << std::endl;
    }
}  // namespace Kratos
