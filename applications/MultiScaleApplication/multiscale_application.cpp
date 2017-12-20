//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:37:00 $
//   Revision:            $Revision: 1.00 $
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"

#include "geometries/point_3d.h"
#include "geometries/line_3d_2.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/quadrilateral_interface_2d_4.h"
#include "geometries/quadrilateral_interface_3d_4.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"
#include "geometries/prism_interface_3d_6.h"
#include "geometries/hexahedra_interface_3d_8.h"
#include "geometries/line_2d_2.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/tetrahedra_3d_4.h"

#include "includes/variables.h"
#include "includes/constitutive_law.h"

#include "multiscale_application.h"
#include "multiscale_application_variables.h"

namespace Kratos {
KratosMultiScaleApplication::KratosMultiScaleApplication()
    : KratosApplication("MultiScaleApplication"),
      mSmallDisplacementInterfaceElement2D4N(
          0, Element::GeometryType::Pointer(
                 new QuadrilateralInterface2D4<Node<3> >(
                     Element::GeometryType::PointsArrayType(4)))),
      mSmallDisplacementInterfaceElement3D6N(
          0, Element::GeometryType::Pointer(new PrismInterface3D6<Node<3> >(
                 Element::GeometryType::PointsArrayType(6)))),
      mSmallDisplacementInterfaceElement3D8N(
          0, Element::GeometryType::Pointer(new HexahedraInterface3D8<Node<3> >(
                 Element::GeometryType::PointsArrayType(8)))),
      mShellThickInterfaceElement3D4N(
          0, Element::GeometryType::Pointer(
                 new QuadrilateralInterface3D4<Node<3> >(
                     Element::GeometryType::PointsArrayType(4)))),
      mOptTriangleElement2D3N(
          0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(
                 Element::GeometryType::PointsArrayType(3)))),
      mConvDiffInterfaceElement2D4N(
          0, Element::GeometryType::Pointer(
                 new QuadrilateralInterface2D4<Node<3> >(
                     Element::GeometryType::PointsArrayType(4)))),
      mConvDiffInterfaceElement3D6N(
          0, Element::GeometryType::Pointer(new PrismInterface3D6<Node<3> >(
                 Element::GeometryType::PointsArrayType(6)))),
      mConvDiffInterfaceElement3D8N(
          0, Element::GeometryType::Pointer(new HexahedraInterface3D8<Node<3> >(
                 Element::GeometryType::PointsArrayType(8)))),
      mConvDiffElement2D3N(
          0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(
                 Element::GeometryType::PointsArrayType(3)))),
      mConvDiffElement2D4N(
          0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3> >(
                 Element::GeometryType::PointsArrayType(4)))),
      mConvDiffElement3D4N(
          0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(
                 Element::GeometryType::PointsArrayType(4)))),
      mConvDiffElement3D8N(
          0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(
                 Element::GeometryType::PointsArrayType(8)))),
      mEASQuadElementV22D4N(
          0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3> >(
                 Element::GeometryType::PointsArrayType(4)))),
      mQ4RIStabElement2D4N(
          0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3> >(
                 Element::GeometryType::PointsArrayType(4)))),
      mEBSTElement2D3N(
          0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(
                 Element::GeometryType::PointsArrayType(3)))),
      mAGQ4Element2D4N(
          0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3> >(
                 Element::GeometryType::PointsArrayType(4)))),
      mSmallDisplacementElasticLinkElement3D2(
          0, Element::GeometryType::Pointer(new Line3D2<Node<3> >(
                 Element::GeometryType::PointsArrayType(2)))),
      mPeriodicConditionLM2D2N(
          0, Condition::GeometryType::Pointer(new Line2D2<Node<3> >(
                 Condition::GeometryType::PointsArrayType(2)))) {}

void KratosMultiScaleApplication::Register() {
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosMultiScaleApplication... " << std::endl;

    // Register variables

    // for stress prediction
    KRATOS_REGISTER_VARIABLE(LCH_REF_RVE)
    KRATOS_REGISTER_VARIABLE(EQUIVALENT_DAMAGE)
    KRATOS_REGISTER_VARIABLE(PREDICTED_STRESS_TENSOR)
    KRATOS_REGISTER_VARIABLE(RVE_DAMAGE_SURFACE_FLAG)
    KRATOS_REGISTER_VARIABLE(RVE_PREDICTION_FLAG)

    // for thermal case
    KRATOS_REGISTER_VARIABLE(RVE_FULL_TEMPERATURE)
    KRATOS_REGISTER_VARIABLE(INITIAL_TEMP_GRAD)

    // for rve
    KRATOS_REGISTER_VARIABLE(RVE_NON_LINEAR_FLAG)
    KRATOS_REGISTER_VARIABLE(RVE_CONSTITUTIVE_LAW_FLAG)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(RVE_FULL_DISPLACEMENT)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(RVE_FULL_ROTATION)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(RVE_WPC_LAGRANGIAN_DOF)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(RVE_WPC_LAGRANGIAN_REACTION)
    KRATOS_REGISTER_VARIABLE(RVE_WPR_LAGRANGIAN_DOF)
    KRATOS_REGISTER_VARIABLE(RVE_WPR_LAGRANGIAN_REACTION)
    KRATOS_REGISTER_VARIABLE(RVE_SHELL_WRC_LAGRANGIAN_DOF_X)
    KRATOS_REGISTER_VARIABLE(RVE_SHELL_WRC_LAGRANGIAN_DOF_REACTION_X)
    KRATOS_REGISTER_VARIABLE(RVE_SHELL_WRC_LAGRANGIAN_DOF_Y)
    KRATOS_REGISTER_VARIABLE(RVE_SHELL_WRC_LAGRANGIAN_DOF_REACTION_Y)

    // for lagrange multipliers
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_LAGRANGE)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ROTATION_LAGRANGE)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(REACTION_DISPLACEMENT_LAGRANGE)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(REACTION_ROTATION_LAGRANGE)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_DOUBLE_LAGRANGE_1)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_DOUBLE_LAGRANGE_2)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ROTATION_DOUBLE_LAGRANGE_1)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ROTATION_DOUBLE_LAGRANGE_2)
    KRATOS_REGISTER_VARIABLE(DOUBLE_LAGRANGE_SCALE_FACTOR_ALPHA)
    KRATOS_REGISTER_VARIABLE(DOUBLE_LAGRANGE_SCALE_FACTOR_BETA)

    // for strategies
    KRATOS_REGISTER_VARIABLE(STRATEGY_SOLUTION_STEP_SOLVED)
    KRATOS_REGISTER_VARIABLE(STRATEGY_FINALIZE_SOLUTION_STEP_LEVEL)
    KRATOS_REGISTER_VARIABLE(CONSTITUTIVE_INTEGRATION_ERROR_CODE)
    KRATOS_REGISTER_VARIABLE(ITERATION_CONVERGENCE_FLAG)
    KRATOS_REGISTER_VARIABLE(SUGGESTED_TIME_STEP)

    // for damage constitutive law
    KRATOS_REGISTER_VARIABLE(DAMAGE_T)
    KRATOS_REGISTER_VARIABLE(DAMAGE_C)
    KRATOS_REGISTER_VARIABLE(FRACTURE_ENERGY_T)
    KRATOS_REGISTER_VARIABLE(FRACTURE_ENERGY_C)
    KRATOS_REGISTER_VARIABLE(YIELD_STRESS_T)
    KRATOS_REGISTER_VARIABLE(YIELD_STRESS_C)
    KRATOS_REGISTER_VARIABLE(DAMAGE_STRESS_T_0)
    KRATOS_REGISTER_VARIABLE(DAMAGE_STRESS_C_0)
    KRATOS_REGISTER_VARIABLE(DAMAGE_STRESS_C_P)
    KRATOS_REGISTER_VARIABLE(DAMAGE_STRESS_C_M)
    KRATOS_REGISTER_VARIABLE(DAMAGE_STRESS_C_R)
    KRATOS_REGISTER_VARIABLE(DAMAGE_STRAIN_C_P)
    KRATOS_REGISTER_VARIABLE(DAMAGE_STRAIN_C_M)
    KRATOS_REGISTER_VARIABLE(DAMAGE_COMPRESSIVE_LAW_C1)
    KRATOS_REGISTER_VARIABLE(DAMAGE_COMPRESSIVE_LAW_C2)
    KRATOS_REGISTER_VARIABLE(DAMAGE_COMPRESSIVE_LAW_C3)
    KRATOS_REGISTER_VARIABLE(BIAXIAL_COMPRESSION_MULTIPLIER)
    KRATOS_REGISTER_VARIABLE(SHEAR_COMPRESSION_REDUCTION)
    KRATOS_REGISTER_VARIABLE(DAMAGE_TENSILE_SURFACE_S1)
    KRATOS_REGISTER_VARIABLE(LUBLINER_SURFACE_PARAM_KC)
    KRATOS_REGISTER_VARIABLE(GENRANKINE_SURFACE_PARAM_A)
    KRATOS_REGISTER_VARIABLE(GENRANKINE_SURFACE_PARAM_B)
    KRATOS_REGISTER_VARIABLE(GENRANKINE_SURFACE_PARAM_C)
    KRATOS_REGISTER_VARIABLE(DAMAGE_SECANT_MATRIX)
    KRATOS_REGISTER_VARIABLE(DAMAGE_MODEL)
    KRATOS_REGISTER_VARIABLE(DAMAGE_TENSILE_MODEL)

    // for interface constitutive law
    KRATOS_REGISTER_VARIABLE(NORMAL_STIFFNESS)
    KRATOS_REGISTER_VARIABLE(TANGENTIAL_STIFFNESS)
    KRATOS_REGISTER_VARIABLE(NORMAL_STIFFNESS_COMPRESSION_MULTIPLIER)
    KRATOS_REGISTER_VARIABLE(FRACTURE_ENERGY_MODE_I)
    KRATOS_REGISTER_VARIABLE(FRACTURE_ENERGY_MODE_II)
    KRATOS_REGISTER_VARIABLE(FRACTURE_ENERGY_MODE_III)
    KRATOS_REGISTER_VARIABLE(EQUIVALENT_PLASTIC_DISPLACEMENT_JUMP_MODE_I)
    KRATOS_REGISTER_VARIABLE(EQUIVALENT_PLASTIC_DISPLACEMENT_JUMP_MODE_II)
    KRATOS_REGISTER_VARIABLE(EQUIVALENT_PLASTIC_DISPLACEMENT_JUMP_MODE_III)
    KRATOS_REGISTER_VARIABLE(INITIAL_COHESION)
    KRATOS_REGISTER_VARIABLE(INITIAL_FRICTION_ANGLE)
    KRATOS_REGISTER_VARIABLE(RESIDUAL_FRICTION_ANGLE)
    KRATOS_REGISTER_VARIABLE(INITIAL_DILATANCY_ANGLE)
    KRATOS_REGISTER_VARIABLE(RESIDUAL_DILATANCY_ANGLE)
    KRATOS_REGISTER_VARIABLE(INTERFACE_TENSILE_LAW_S0)
    KRATOS_REGISTER_VARIABLE(INTERFACE_COMPRESSIVE_LAW_S0)
    KRATOS_REGISTER_VARIABLE(INTERFACE_COMPRESSIVE_LAW_SP)
    KRATOS_REGISTER_VARIABLE(INTERFACE_COMPRESSIVE_LAW_SR)
    KRATOS_REGISTER_VARIABLE(INTERFACE_COMPRESSIVE_LAW_EP)
    KRATOS_REGISTER_VARIABLE(INTERFACE_COMPRESSIVE_LAW_C1)
    KRATOS_REGISTER_VARIABLE(INTERFACE_COMPRESSIVE_LAW_C2)
    KRATOS_REGISTER_VARIABLE(INTERFACE_COMPRESSIVE_LAW_C3)
    KRATOS_REGISTER_VARIABLE(INTERFACE_COMPRESSIVE_LAW_C4)
    KRATOS_REGISTER_VARIABLE(INTERFACE_PLASTIC_DAMAGE_FACTOR_T)
    KRATOS_REGISTER_VARIABLE(INTERFACE_PLASTIC_DAMAGE_FACTOR_C)
    KRATOS_REGISTER_VARIABLE(INTERFACE_CAP_VALUE)
    KRATOS_REGISTER_VARIABLE(INTERFACE_TRACTION)
    KRATOS_REGISTER_VARIABLE(INTERFACE_DISPLACEMENT_JUMP)
    KRATOS_REGISTER_VARIABLE(INTERFACE_COUPLE_TRACTION)
    KRATOS_REGISTER_VARIABLE(INTERFACE_ROTATION_JUMP)
    KRATOS_REGISTER_VARIABLE(INTERFACE_PLASTIC_DISPLACEMENT_JUMP)
    KRATOS_REGISTER_VARIABLE(YIELD_FUNCTION_VALUE)
    KRATOS_REGISTER_VARIABLE(INTERFACE_REDUCED_INTEGRATION)

    // for plots
    KRATOS_REGISTER_VARIABLE(YIELD_SURFACE_DATA_2D_X)
    KRATOS_REGISTER_VARIABLE(YIELD_SURFACE_DATA_2D_Y)
    KRATOS_REGISTER_VARIABLE(YIELD_SURFACE_DATA_3D_X)
    KRATOS_REGISTER_VARIABLE(YIELD_SURFACE_DATA_3D_Y)
    KRATOS_REGISTER_VARIABLE(YIELD_SURFACE_DATA_3D_Z)

    // for plastic constitutive law
    KRATOS_REGISTER_VARIABLE(ISOTROPIC_HARDENING)
    KRATOS_REGISTER_VARIABLE(KINEMATIC_HARDENING)
    KRATOS_REGISTER_VARIABLE(YIELD_STRESS_INFINITY)
    KRATOS_REGISTER_VARIABLE(ISOTROPIC_HARDENING_EXPONENT)
    KRATOS_REGISTER_VARIABLE(EQUIVALENT_PLASTIC_STRAIN)
    KRATOS_REGISTER_VARIABLE(PLASTIC_STRAIN_TENSOR)

    // for stabilized reduced integration
    KRATOS_REGISTER_VARIABLE(RI_STABILIZATION)
    KRATOS_REGISTER_VARIABLE(RI_STABILIZATION_RESIDUAL)

    // for enhanced strain elements
    KRATOS_REGISTER_VARIABLE(ENH_STRAIN_PARAM_1)
    KRATOS_REGISTER_VARIABLE(ENH_STRAIN_PARAM_2)
    KRATOS_REGISTER_VARIABLE(ENH_STRAIN_PARAM_3)
    KRATOS_REGISTER_VARIABLE(ENH_STRAIN_PARAM_4)
    KRATOS_REGISTER_VARIABLE(ENH_STRAIN_PARAM_5)

    // misc
    KRATOS_REGISTER_VARIABLE(RANDOM_IMPERFECTION_FACTOR)
    KRATOS_REGISTER_VARIABLE(DISCONTINUITY_DIRECTION)
    KRATOS_REGISTER_VARIABLE(LAMBDA_OUTPUT)

    // for shell load conditions
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SHELL_EDGE_FORCE)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SHELL_EDGE_TORQUE)

    // Register Elements
    KRATOS_REGISTER_ELEMENT("SmallDisplacementInterfaceElement2D4N",
        mSmallDisplacementInterfaceElement2D4N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementInterfaceElement3D6N",
        mSmallDisplacementInterfaceElement3D6N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementInterfaceElement3D8N",
        mSmallDisplacementInterfaceElement3D8N)
    KRATOS_REGISTER_ELEMENT(
        "ShellThickInterfaceElement3D4N", mShellThickInterfaceElement3D4N)
    KRATOS_REGISTER_ELEMENT("OptTriangleElement2D3N", mOptTriangleElement2D3N)
    KRATOS_REGISTER_ELEMENT("EASQuadElementV22D4N", mEASQuadElementV22D4N)
    KRATOS_REGISTER_ELEMENT("Q4RIStabElement2D4N", mQ4RIStabElement2D4N)
    KRATOS_REGISTER_ELEMENT(
        "ConvDiffInterfaceElement2D4N", mConvDiffInterfaceElement2D4N)
    KRATOS_REGISTER_ELEMENT(
        "ConvDiffInterfaceElement3D6N", mConvDiffInterfaceElement3D6N)
    KRATOS_REGISTER_ELEMENT(
        "ConvDiffInterfaceElement3D8N", mConvDiffInterfaceElement3D8N)
    KRATOS_REGISTER_ELEMENT("ConvDiffElement2D3N", mConvDiffElement2D3N)
    KRATOS_REGISTER_ELEMENT("ConvDiffElement2D4N", mConvDiffElement2D4N)
    KRATOS_REGISTER_ELEMENT("ConvDiffElement3D4N", mConvDiffElement3D4N)
    KRATOS_REGISTER_ELEMENT("ConvDiffElement3D8N", mConvDiffElement3D8N)
    KRATOS_REGISTER_ELEMENT("EBSTElement2D3N", mEBSTElement2D3N)
    KRATOS_REGISTER_ELEMENT("AGQ4Element2D4N", mAGQ4Element2D4N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElasticLinkElement3D2N",
        mSmallDisplacementElasticLinkElement3D2)

    // Register Conditions
    KRATOS_REGISTER_CONDITION(
        "PeriodicConditionLM2D2N", mPeriodicConditionLM2D2N)
}

}  // namespace Kratos.
