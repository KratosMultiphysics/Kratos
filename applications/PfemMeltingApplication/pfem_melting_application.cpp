// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Julio Marti
//


// System includes

// External includes
//

// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_27.h"
#include "pfem_melting_application.h"
#include "includes/variables.h"

namespace Kratos {



KratosPfemMeltingApplication::KratosPfemMeltingApplication()
    : KratosApplication("PfemMeltingApplication")

      {}

void KratosPfemMeltingApplication::Register() {
    KRATOS_INFO("") <<
    " KRATOS ___ ___  _  ___   __   ___ ___ ___ ___ " << std::endl <<
    "       / __/ _ || || | | / /__|   |_ _| __| __|" << std::endl <<
    "      | (_| (_) | .` || V /___| |) | || _|| _| " << std::endl <<
    "       |___|___/|_||_| |_/    |___/___|_| |_|  APPLICATION" << std::endl;

    // Registering variables
    KRATOS_REGISTER_VARIABLE(ACTIVATION_ENERGY)
    KRATOS_REGISTER_VARIABLE(ARRHENIUS_COEFFICIENT)
    KRATOS_REGISTER_VARIABLE(RADIOUS)

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(INITIAL_POSITION)


}

}  // namespace Kratos.
