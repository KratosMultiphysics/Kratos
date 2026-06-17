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
#include "geometries/line_2d_2.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_27.h"
#include "viscosity_modulator_application.h"
#include "includes/variables.h"

namespace Kratos {



KratosViscosityModulatorApplication::KratosViscosityModulatorApplication()
    : KratosApplication("ViscosityModulatorApplication"),
      // Shock-capturing convection-diffusion element
      mEulerianConvDiffShockCapturing2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node>(Element::GeometryType::PointsArrayType(3)))),
      mEulerianConvDiffShockCapturing2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node>(Element::GeometryType::PointsArrayType(4)))),
      mEulerianConvDiffShockCapturing3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node>(Element::GeometryType::PointsArrayType(4)))),
      mEulerianConvDiffShockCapturing3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node>(Element::GeometryType::PointsArrayType(8)))),
      // FluxCondition
      mFluxCondition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<Node>(Condition::GeometryType::PointsArrayType(2)))),
      mFluxCondition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<Node>(Condition::GeometryType::PointsArrayType(3)))),
      mFluxCondition3D4N(0, Condition::GeometryType::Pointer(new Quadrilateral3D4<Node>(Condition::GeometryType::PointsArrayType(4)))),
      // ConsistentFluxBoundaryCondition
      mConsistentFluxBoundaryCondition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<Node>(Condition::GeometryType::PointsArrayType(2)))),
      mConsistentFluxBoundaryCondition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<Node>(Condition::GeometryType::PointsArrayType(3)))),
      mConsistentFluxBoundaryCondition3D4N(0, Condition::GeometryType::Pointer(new Quadrilateral3D4<Node>(Condition::GeometryType::PointsArrayType(4))))
{}

void KratosViscosityModulatorApplication::Register() {
    KRATOS_INFO("") <<
    " KRATOS" << std::endl <<
    "      _                                    " << std::endl <<
    "     | | | (_)        |  \\/  |         | |" << std::endl <<
    "     | | | |_ ___  ___| .  . | ___   __| |" << std::endl <<
    "     | | | | / __|/ __| |\\/| |/ _ \\ / _` |" << std::endl <<
    "     \\ \\_/ / \\__ \\ (__| |  | | (_) | (_| |" << std::endl <<
    "      \\___/|_|___/\\___\\_|  |_/\\___/ \\__,_|" << "VISCOSITY MODULATOR APPLICATION\n" << std::endl;

    // Register CD variables (reuse CD's objects; no new variables defined)
    KRATOS_REGISTER_VARIABLE(AUX_TEMPERATURE)
    KRATOS_REGISTER_VARIABLE(EXACT_PRESSURE)
    KRATOS_REGISTER_VARIABLE(PROJECTED_SCALAR1)
    KRATOS_REGISTER_VARIABLE(TRANSFER_COEFFICIENT)
    KRATOS_REGISTER_VARIABLE(SHOCK_CAPTURING_INTENSITY)
    KRATOS_REGISTER_VARIABLE(USE_ANISOTROPIC_DISC_CAPTURING)
    KRATOS_REGISTER_VARIABLE(STEP_SOLUTION)
    KRATOS_REGISTER_VARIABLE(TRUNC_SOLUTION)
    KRATOS_REGISTER_VARIABLE(TRUNC_SOLUTION_ERROR)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(REFERENCE_VELOCITY)

    // Elements — registered under Vm-prefixed strings so they do not collide
    // with ConvectionDiffusionApplication's own element of the same base name
    // when both applications are imported in the same process.
    KRATOS_REGISTER_ELEMENT("VmEulerianConvDiffShockCapturing2D3N", mEulerianConvDiffShockCapturing2D3N);
    KRATOS_REGISTER_ELEMENT("VmEulerianConvDiffShockCapturing2D4N", mEulerianConvDiffShockCapturing2D4N);
    KRATOS_REGISTER_ELEMENT("VmEulerianConvDiffShockCapturing3D4N", mEulerianConvDiffShockCapturing3D4N);
    KRATOS_REGISTER_ELEMENT("VmEulerianConvDiffShockCapturing3D8N", mEulerianConvDiffShockCapturing3D8N);

    // Conditions — registered under the same strings as ConvectionDiffusionApplication
    // Conditions — Vm-prefixed so they do not collide with
    // ConvectionDiffusionApplication's own conditions of the same base name
    // when both applications are imported in the same process.
    KRATOS_REGISTER_CONDITION("VmFluxCondition2D2N", mFluxCondition2D2N);
    KRATOS_REGISTER_CONDITION("VmFluxCondition3D3N", mFluxCondition3D3N);
    KRATOS_REGISTER_CONDITION("VmFluxCondition3D4N", mFluxCondition3D4N);
    KRATOS_REGISTER_CONDITION("VmConsistentFluxBoundaryCondition2D2N", mConsistentFluxBoundaryCondition2D2N);
    KRATOS_REGISTER_CONDITION("VmConsistentFluxBoundaryCondition3D3N", mConsistentFluxBoundaryCondition3D3N);
    KRATOS_REGISTER_CONDITION("VmConsistentFluxBoundaryCondition3D4N", mConsistentFluxBoundaryCondition3D4N);

}

}  // namespace Kratos.
