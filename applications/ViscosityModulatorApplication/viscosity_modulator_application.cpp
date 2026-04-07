// KRATOS
//  _   _
// | | | (_)        |  \/  |         | |
// | | | |_ ___  ___| .  . | ___   __| |
// | | | | / __|/ __| |\/| |/ _ \ / _` |
// \ \_/ / \__ \ (__| |  | | (_) | (_| |
//  \___/|_|___/\___\_|  |_/\___/ \__,_|  APPLICATION
//                                      
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aniol Sala
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
#include "viscosity_modulator_application.h"
#include "includes/variables.h"
#include "includes/viscosity_modulator_settings.h"

namespace Kratos {



KratosViscosityModulatorApplication::KratosViscosityModulatorApplication()
    : KratosApplication("ViscosityModulatorApplication"),
      mViscosityModulatorElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node >(Element::GeometryType::PointsArrayType(3)))),
      mViscosityModulatorElement2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mViscosityModulatorElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node >(Element::GeometryType::PointsArrayType(4)))),
      mViscosityModulatorElement3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node >(Element::GeometryType::PointsArrayType(8)))) {}

void KratosViscosityModulatorApplication::Register() {
    KRATOS_INFO("") <<
    " KRATOS" << std::endl <<
    "      _                                    " << std::endl <<
    "     | | | (_)        |  \\/  |         | |" << std::endl <<
    "     | | | |_ ___  ___| .  . | ___   __| |" << std::endl <<
    "     | | | | / __|/ __| |\\/| |/ _ \\ / _` |" << std::endl <<
    "     \\ \\_/ / \\__ \\ (__| |  | | (_) | (_| |" << std::endl <<
    "      \\___/|_|___/\\___\\_|  |_/\\___/ \\__,_|" << "VISCOSITY MODULATOR APPLICATION\n" << std::endl;

    // Registering variables
    KRATOS_REGISTER_VARIABLE(AUX_SCALAR)
    KRATOS_REGISTER_VARIABLE(SHOCK_CAPTURING_INTENSITY)
    KRATOS_REGISTER_VARIABLE(USE_ANISOTROPIC_DISC_CAPTURING)
    KRATOS_REGISTER_VARIABLE(VISCOSITY_MODULATOR_SETTINGS)

    // Registering elements and conditions here
    KRATOS_REGISTER_ELEMENT("ViscosityModulatorElement2D3N", mViscosityModulatorElement2D3N);
    KRATOS_REGISTER_ELEMENT("ViscosityModulatorElement2D4N", mViscosityModulatorElement2D4N);
    KRATOS_REGISTER_ELEMENT("ViscosityModulatorElement3D4N", mViscosityModulatorElement3D4N);
    KRATOS_REGISTER_ELEMENT("ViscosityModulatorElement3D8N", mViscosityModulatorElement3D8N);

}

}  // namespace Kratos.
