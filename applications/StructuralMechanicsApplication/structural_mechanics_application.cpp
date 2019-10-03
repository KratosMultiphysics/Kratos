// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//    Co-authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define.h"

#include "structural_mechanics_application_variables.h"
#include "structural_mechanics_application.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"

#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"

#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"

namespace Kratos {

    // We define the node type
    typedef Node<3> NodeType;

KratosStructuralMechanicsApplication::KratosStructuralMechanicsApplication()
    : KratosApplication("StructuralMechanicsApplication"),
      /* ELEMENTS */
      // Adding the truss elements
      mTrussElement3D2N(0, Element::GeometryType::Pointer(new Line3D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
      mTrussLinearElement3D2N(0, Element::GeometryType::Pointer(new Line3D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
      mCableElement3D2N(0, Element::GeometryType::Pointer(new Line3D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
      // Adding the beam elements
      mCrBeamElement3D2N(0, Element::GeometryType::Pointer(new Line3D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
      mCrLinearBeamElement3D2N(0, Element::GeometryType::Pointer(new Line3D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
      mCrBeamElement2D2N(0, Element::GeometryType::Pointer(new Line2D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
      mCrLinearBeamElement2D2N(0, Element::GeometryType::Pointer(new Line2D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
      // Adding the shells elements
      mIsotropicShellElement3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mShellThickElement3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4<NodeType >(Element::GeometryType::PointsArrayType(4))), false),
      mShellThickCorotationalElement3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4<NodeType >(Element::GeometryType::PointsArrayType(4))), true),
      mShellThinCorotationalElement3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4<NodeType >(Element::GeometryType::PointsArrayType(4))), true),
      mShellThinElement3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<NodeType >(Element::GeometryType::PointsArrayType(3))), false),
      mShellThinCorotationalElement3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<NodeType >(Element::GeometryType::PointsArrayType(3))), true),
      mShellThickCorotationalElement3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<NodeType >(Element::GeometryType::PointsArrayType(3))), true),
      // Adding the membrane element
      mPreStressMembraneElement3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mPreStressMembraneElement3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      // Adding the SPRISM element
      mSolidShellElementSprism3D6N(0, Element::GeometryType::Pointer(new Prism3D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      // Adding the nodal concentrated element
      mNodalConcentratedElement2D1N(0, Element::GeometryType::Pointer(new Point2D<NodeType >(Element::GeometryType::PointsArrayType(1))), true),
      mNodalConcentratedDampedElement2D1N(0, Element::GeometryType::Pointer(new Point2D<NodeType >(Element::GeometryType::PointsArrayType(1))), false),
      mNodalConcentratedElement3D1N(0, Element::GeometryType::Pointer(new Point3D<NodeType >(Element::GeometryType::PointsArrayType(1))), true),
      mNodalConcentratedDampedElement3D1N(0, Element::GeometryType::Pointer(new Point3D<NodeType >(Element::GeometryType::PointsArrayType(1))), false),
      // Adding the kinematic linear elements
      mSmallDisplacement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mSmallDisplacement2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mSmallDisplacement2D6N(0, Element::GeometryType::Pointer(new Triangle2D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mSmallDisplacement2D8N(0, Element::GeometryType::Pointer(new Quadrilateral2D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mSmallDisplacement2D9N(0, Element::GeometryType::Pointer(new Quadrilateral2D9<NodeType >(Element::GeometryType::PointsArrayType(9)))),
      mSmallDisplacement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mSmallDisplacement3D6N(0, Element::GeometryType::Pointer(new Prism3D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mSmallDisplacement3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mSmallDisplacement3D10N(0, Element::GeometryType::Pointer(new Tetrahedra3D10<NodeType >(Element::GeometryType::PointsArrayType(10)))),
      mSmallDisplacement3D15N(0, Element::GeometryType::Pointer(new Prism3D15<NodeType >(Element::GeometryType::PointsArrayType(15)))),
      mSmallDisplacement3D20N(0, Element::GeometryType::Pointer(new Hexahedra3D20<NodeType >(Element::GeometryType::PointsArrayType(20)))),
      mSmallDisplacement3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27<NodeType >(Element::GeometryType::PointsArrayType(27)))),

      mSmallDisplacementBbar2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType>(Element::GeometryType::PointsArrayType(4)))),
      mSmallDisplacementBbar3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType>(Element::GeometryType::PointsArrayType(8)))),

      mAxisymSmallDisplacement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mAxisymSmallDisplacement2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mAxisymSmallDisplacement2D6N(0, Element::GeometryType::Pointer(new Triangle2D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mAxisymSmallDisplacement2D8N(0, Element::GeometryType::Pointer(new Quadrilateral2D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mAxisymSmallDisplacement2D9N(0, Element::GeometryType::Pointer(new Quadrilateral2D9<NodeType >(Element::GeometryType::PointsArrayType(9)))),

      // Adding the Total lagrangian elements
      mTotalLagrangian2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mTotalLagrangian2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mTotalLagrangian2D6N(0, Element::GeometryType::Pointer(new Triangle2D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mTotalLagrangian2D8N(0, Element::GeometryType::Pointer(new Quadrilateral2D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mTotalLagrangian2D9N(0, Element::GeometryType::Pointer(new Quadrilateral2D9<NodeType >(Element::GeometryType::PointsArrayType(9)))),
      mTotalLagrangian3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mTotalLagrangian3D6N(0, Element::GeometryType::Pointer(new Prism3D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mTotalLagrangian3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mTotalLagrangian3D10N(0, Element::GeometryType::Pointer(new Tetrahedra3D10<NodeType >(Element::GeometryType::PointsArrayType(10)))),
      mTotalLagrangian3D15N(0, Element::GeometryType::Pointer(new Prism3D15<NodeType >(Element::GeometryType::PointsArrayType(15)))),
      mTotalLagrangian3D20N(0, Element::GeometryType::Pointer(new Hexahedra3D20<NodeType >(Element::GeometryType::PointsArrayType(20)))),
      mTotalLagrangian3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27<NodeType >(Element::GeometryType::PointsArrayType(27)))),
      mAxisymTotalLagrangian2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mAxisymTotalLagrangian2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mAxisymTotalLagrangian2D6N(0, Element::GeometryType::Pointer(new Triangle2D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mAxisymTotalLagrangian2D8N(0, Element::GeometryType::Pointer(new Quadrilateral2D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mAxisymTotalLagrangian2D9N(0, Element::GeometryType::Pointer(new Quadrilateral2D9<NodeType >(Element::GeometryType::PointsArrayType(9)))),
      // Adding the Updated lagrangian elements
      mUpdatedLagrangian2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mUpdatedLagrangian2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mUpdatedLagrangian2D6N(0, Element::GeometryType::Pointer(new Triangle2D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mUpdatedLagrangian2D8N(0, Element::GeometryType::Pointer(new Quadrilateral2D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mUpdatedLagrangian2D9N(0, Element::GeometryType::Pointer(new Quadrilateral2D9<NodeType >(Element::GeometryType::PointsArrayType(9)))),
      mUpdatedLagrangian3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mUpdatedLagrangian3D6N(0, Element::GeometryType::Pointer(new Prism3D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mUpdatedLagrangian3D8N( 0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mUpdatedLagrangian3D10N(0, Element::GeometryType::Pointer(new Tetrahedra3D10<NodeType >(Element::GeometryType::PointsArrayType(10)))),
      mUpdatedLagrangian3D15N(0, Element::GeometryType::Pointer(new Prism3D15<NodeType >(Element::GeometryType::PointsArrayType(15)))),
      mUpdatedLagrangian3D20N(0, Element::GeometryType::Pointer(new Hexahedra3D20<NodeType >(Element::GeometryType::PointsArrayType(20)))),
      mUpdatedLagrangian3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27<NodeType >(Element::GeometryType::PointsArrayType(27)))),
      mAxisymUpdatedLagrangian2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mAxisymUpdatedLagrangian2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mAxisymUpdatedLagrangian2D6N(0, Element::GeometryType::Pointer(new Triangle2D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mAxisymUpdatedLagrangian2D8N(0, Element::GeometryType::Pointer(new Quadrilateral2D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mAxisymUpdatedLagrangian2D9N(0, Element::GeometryType::Pointer(new Quadrilateral2D9<NodeType >(Element::GeometryType::PointsArrayType(9)))),
      // Adding the spring damper element
      mSpringDamperElement3D2N(0, Element::GeometryType::Pointer(new Line3D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
      // Adding the adjoint elements
      mAdjointFiniteDifferencingShellThinElement3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mAdjointFiniteDifferenceCrBeamElementLinear3D2N(0, Element::GeometryType::Pointer(new Line3D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
      mAdjointFiniteDifferenceTrussElement3D2N(0, Element::GeometryType::Pointer(new Line3D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
      mAdjointFiniteDifferenceTrussLinearElement3D2N(0, Element::GeometryType::Pointer(new Line3D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
      mTotalLagrangianAdjoint2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mTotalLagrangianAdjoint2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mTotalLagrangianAdjoint2D6N(0, Element::GeometryType::Pointer(new Triangle2D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mTotalLagrangianAdjoint3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mTotalLagrangianAdjoint3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mAdjointFiniteDifferencingSmallDisplacementElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mAdjointFiniteDifferencingSmallDisplacementElement3D6N(0, Element::GeometryType::Pointer(new Prism3D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mAdjointFiniteDifferencingSmallDisplacementElement3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),

      /* CONDITIONS */
      // Adding point load conditions
      mPointLoadCondition2D1N(0, Condition::GeometryType::Pointer(new Point2D<NodeType >(Condition::GeometryType::PointsArrayType(1)))),
      mPointLoadCondition3D1N(0, Condition::GeometryType::Pointer(new Point3D<NodeType >(Condition::GeometryType::PointsArrayType(1)))),
      mPointContactCondition2D1N(0, Condition::GeometryType::Pointer(new Point2D<NodeType >(Condition::GeometryType::PointsArrayType(1)))),
      mPointContactCondition3D1N(0, Condition::GeometryType::Pointer(new Point3D<NodeType >(Condition::GeometryType::PointsArrayType(1)))),
      mAxisymPointLoadCondition2D1N(0, Condition::GeometryType::Pointer(new Point2D<NodeType >(Condition::GeometryType::PointsArrayType(1)))),
      // Adding line load conditions
      mLineLoadCondition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mLineLoadCondition2D3N(0, Condition::GeometryType::Pointer(new Line2D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mLineLoadCondition3D2N(0, Condition::GeometryType::Pointer(new Line3D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mLineLoadCondition3D3N(0, Condition::GeometryType::Pointer(new Line3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mAxisymLineLoadCondition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mAxisymLineLoadCondition2D3N(0, Condition::GeometryType::Pointer(new Line2D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      // Adding surface load conditions
      mSurfaceLoadCondition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mSurfaceLoadCondition3D4N(0, Condition::GeometryType::Pointer(new Quadrilateral3D4<NodeType >( Condition::GeometryType::PointsArrayType(4)))),
      mSurfaceLoadCondition3D6N(0, Condition::GeometryType::Pointer(new Triangle3D6<NodeType >(Condition::GeometryType::PointsArrayType(6)))),
      mSurfaceLoadCondition3D8N(0, Condition::GeometryType::Pointer(new Quadrilateral3D8<NodeType >(Condition::GeometryType::PointsArrayType(8)))),
      mSurfaceLoadCondition3D9N(0, Condition::GeometryType::Pointer(new Quadrilateral3D9<NodeType >(Condition::GeometryType::PointsArrayType(9)))),
      // Adding point moment conditions
      mPointMomentCondition3D1N(0, Condition::GeometryType::Pointer(new Point3D<NodeType >(Condition::GeometryType::PointsArrayType(1)))),

      // Adding adjoint conditions
      mAdjointSemiAnalyticPointLoadCondition2D1N(0, Condition::GeometryType::Pointer(new Point2D<NodeType >(Condition::GeometryType::PointsArrayType(1)))),
      mAdjointSemiAnalyticPointLoadCondition3D1N(0, Condition::GeometryType::Pointer(new Point3D<NodeType >(Condition::GeometryType::PointsArrayType(1)))){}

void KratosStructuralMechanicsApplication::Register() {
    // calling base class register to register Kratos components
    KratosApplication::Register();

    KRATOS_INFO("") << "    KRATOS   ___|  |                   |                   |\n"
                    << "           \\___ \\  __|  __| |   |  __| __| |   |  __| _` | |\n"
                    << "                 | |   |    |   | (    |   |   | |   (   | |\n"
                    << "           _____/ \\__|_|   \\__,_|\\___|\\__|\\__,_|_|  \\__,_|_| MECHANICS\n"
                    << "Initializing KratosStructuralMechanicsApplication..." << std::endl;

    // General pourpose
    KRATOS_REGISTER_VARIABLE(INTEGRATION_ORDER); // The integration order considered on the element
    KRATOS_REGISTER_VARIABLE(LOCAL_MATERIAL_AXIS_1);
    KRATOS_REGISTER_VARIABLE(LOCAL_MATERIAL_AXIS_2);
    KRATOS_REGISTER_VARIABLE(LOCAL_MATERIAL_AXIS_3);
    KRATOS_REGISTER_VARIABLE(CENTER_OF_GRAVITY);
    KRATOS_REGISTER_VARIABLE(MASS_MOMENT_OF_INERTIA);
    KRATOS_REGISTER_VARIABLE(ELASTICITY_TENSOR);


    // Generalized eigenvalue problem
    KRATOS_REGISTER_VARIABLE(BUILD_LEVEL)
    KRATOS_REGISTER_VARIABLE(EIGENVALUE_VECTOR)
    KRATOS_REGISTER_VARIABLE(EIGENVECTOR_MATRIX)
    KRATOS_REGISTER_VARIABLE(MODAL_MASS_MATRIX)
    KRATOS_REGISTER_VARIABLE(MODAL_STIFFNESS_MATRIX)

    // Geometrical
    KRATOS_REGISTER_VARIABLE(AREA)
    KRATOS_REGISTER_VARIABLE(IT)
    KRATOS_REGISTER_VARIABLE(IY)
    KRATOS_REGISTER_VARIABLE(IZ)
    KRATOS_REGISTER_VARIABLE(CROSS_AREA)
    KRATOS_REGISTER_VARIABLE(MEAN_RADIUS)
    KRATOS_REGISTER_VARIABLE(SECTION_SIDES)
    KRATOS_REGISTER_VARIABLE(GEOMETRIC_STIFFNESS)
    KRATOS_REGISTER_VARIABLE(LOCAL_ELEMENT_ORIENTATION)
    KRATOS_REGISTER_VARIABLE(MATERIAL_ORIENTATION_ANGLE)
    KRATOS_REGISTER_VARIABLE(USE_CONSISTENT_MASS_MATRIX)
    KRATOS_REGISTER_VARIABLE(CONDENSED_DOF_LIST)

    // Truss generalized variables
    KRATOS_REGISTER_VARIABLE(TRUSS_PRESTRESS_PK2)
    KRATOS_REGISTER_VARIABLE(HARDENING_MODULUS_1D)
    KRATOS_REGISTER_VARIABLE(TANGENT_MODULUS)
    KRATOS_REGISTER_VARIABLE(PLASTIC_ALPHA)

    // Beam generalized variables
    KRATOS_REGISTER_VARIABLE(AREA_EFFECTIVE_Y)
    KRATOS_REGISTER_VARIABLE(AREA_EFFECTIVE_Z)
    KRATOS_REGISTER_VARIABLE(INERTIA_ROT_Y)
    KRATOS_REGISTER_VARIABLE(INERTIA_ROT_Z)
    KRATOS_REGISTER_VARIABLE(LOCAL_AXES_VECTOR)
    KRATOS_REGISTER_VARIABLE(TORSIONAL_INERTIA)
    KRATOS_REGISTER_VARIABLE(I22)
    KRATOS_REGISTER_VARIABLE(I33)
    KRATOS_REGISTER_VARIABLE(LUMPED_MASS_ROTATION_COEFFICIENT)

    //  Shell generalized variables
    KRATOS_REGISTER_VARIABLE(STENBERG_SHEAR_STABILIZATION_SUITABLE)
    KRATOS_REGISTER_VARIABLE(SHELL_OFFSET)
    KRATOS_REGISTER_VARIABLE(SHELL_STRAIN)
    KRATOS_REGISTER_VARIABLE(SHELL_FORCE)
    KRATOS_REGISTER_VARIABLE(SHELL_STRAIN_GLOBAL)
    KRATOS_REGISTER_VARIABLE(SHELL_FORCE_GLOBAL)
    KRATOS_REGISTER_VARIABLE(SHELL_CURVATURE)
    KRATOS_REGISTER_VARIABLE(SHELL_CURVATURE_GLOBAL)
    KRATOS_REGISTER_VARIABLE(SHELL_MOMENT)
    KRATOS_REGISTER_VARIABLE(SHELL_MOMENT_GLOBAL)
    KRATOS_REGISTER_VARIABLE(SHELL_STRESS_TOP_SURFACE)
    KRATOS_REGISTER_VARIABLE(SHELL_STRESS_TOP_SURFACE_GLOBAL)
    KRATOS_REGISTER_VARIABLE(SHELL_STRESS_MIDDLE_SURFACE)
    KRATOS_REGISTER_VARIABLE(SHELL_STRESS_MIDDLE_SURFACE_GLOBAL)
    KRATOS_REGISTER_VARIABLE(SHELL_STRESS_BOTTOM_SURFACE)
    KRATOS_REGISTER_VARIABLE(SHELL_STRESS_BOTTOM_SURFACE_GLOBAL)
    KRATOS_REGISTER_VARIABLE(VON_MISES_STRESS_TOP_SURFACE)
    KRATOS_REGISTER_VARIABLE(VON_MISES_STRESS_MIDDLE_SURFACE)
    KRATOS_REGISTER_VARIABLE(VON_MISES_STRESS_BOTTOM_SURFACE)
    KRATOS_REGISTER_VARIABLE(SHEAR_ANGLE)
    KRATOS_REGISTER_VARIABLE(SHELL_ORTHOTROPIC_STRESS_BOTTOM_SURFACE)
    KRATOS_REGISTER_VARIABLE(SHELL_ORTHOTROPIC_STRESS_TOP_SURFACE)
    KRATOS_REGISTER_VARIABLE(SHELL_ORTHOTROPIC_STRESS_BOTTOM_SURFACE_GLOBAL)
    KRATOS_REGISTER_VARIABLE(SHELL_ORTHOTROPIC_STRESS_TOP_SURFACE_GLOBAL)
    KRATOS_REGISTER_VARIABLE(SHELL_ORTHOTROPIC_4PLY_THROUGH_THICKNESS)
    KRATOS_REGISTER_VARIABLE(TSAI_WU_RESERVE_FACTOR)
    KRATOS_REGISTER_VARIABLE(SHELL_ORTHOTROPIC_LAMINA_STRENGTHS)

    // Shell energies
    KRATOS_REGISTER_VARIABLE(SHELL_ELEMENT_MEMBRANE_ENERGY)
    KRATOS_REGISTER_VARIABLE(SHELL_ELEMENT_BENDING_ENERGY)
    KRATOS_REGISTER_VARIABLE(SHELL_ELEMENT_SHEAR_ENERGY)
    KRATOS_REGISTER_VARIABLE(SHELL_ELEMENT_MEMBRANE_ENERGY_FRACTION)
    KRATOS_REGISTER_VARIABLE(SHELL_ELEMENT_BENDING_ENERGY_FRACTION)
    KRATOS_REGISTER_VARIABLE(SHELL_ELEMENT_SHEAR_ENERGY_FRACTION)

    // Prestresse membrane generalized vairiables
    KRATOS_REGISTER_VARIABLE( MEMBRANE_PRESTRESS )
    KRATOS_REGISTER_VARIABLE( PRESTRESS_VECTOR )
    KRATOS_REGISTER_VARIABLE( PRESTRESS_AXIS_1_GLOBAL )
    KRATOS_REGISTER_VARIABLE( PRESTRESS_AXIS_2_GLOBAL )
    KRATOS_REGISTER_VARIABLE( PRESTRESS_AXIS_1 )
    KRATOS_REGISTER_VARIABLE( PRESTRESS_AXIS_2 )
    KRATOS_REGISTER_VARIABLE( PROJECTION_TYPE_COMBO )

    // Formfinding
    KRATOS_REGISTER_VARIABLE(LAMBDA_MAX)
    KRATOS_REGISTER_VARIABLE(IS_FORMFINDING)
    KRATOS_REGISTER_VARIABLE(BASE_REF_1)
    KRATOS_REGISTER_VARIABLE(BASE_REF_2)

    // Cross section
    KRATOS_REGISTER_VARIABLE(SHELL_CROSS_SECTION)
    KRATOS_REGISTER_VARIABLE(SHELL_CROSS_SECTION_OUTPUT_PLY_ID)
    KRATOS_REGISTER_VARIABLE(SHELL_CROSS_SECTION_OUTPUT_PLY_LOCATION)
    KRATOS_REGISTER_VARIABLE(SHELL_ORTHOTROPIC_LAYERS)

    // Nodal stiffness
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NODAL_INITIAL_DISPLACEMENT)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NODAL_DISPLACEMENT_STIFFNESS)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NODAL_INITIAL_ROTATION)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NODAL_ROTATIONAL_STIFFNESS)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NODAL_DAMPING_RATIO)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NODAL_ROTATIONAL_DAMPING_RATIO)

    // For explicit central difference scheme
    KRATOS_REGISTER_VARIABLE(MASS_FACTOR)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MIDDLE_VELOCITY)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MIDDLE_ANGULAR_VELOCITY)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(FRACTIONAL_ACCELERATION)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(FRACTIONAL_ANGULAR_ACCELERATION)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NODAL_INERTIA)
    KRATOS_REGISTER_VARIABLE(NODAL_DISPLACEMENT_DAMPING)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NODAL_ROTATION_DAMPING)

    // CONDITIONS
    /* Moment condition */
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(POINT_MOMENT)

    // Adding the SPRISM EAS variables
    KRATOS_REGISTER_VARIABLE(ALPHA_EAS);
    KRATOS_REGISTER_VARIABLE(CONSIDER_IMPLICIT_EAS_SPRISM_ELEMENT);
    KRATOS_REGISTER_VARIABLE(CONSIDER_TOTAL_LAGRANGIAN_SPRISM_ELEMENT);
    KRATOS_REGISTER_VARIABLE(PURE_EXPLICIT_RHS_COMPUTATION);

    // Reset equations ids "flag"
    KRATOS_REGISTER_VARIABLE(RESET_EQUATION_IDS);

    // Adding the SPRISM additional variables
    KRATOS_REGISTER_VARIABLE(ANG_ROT);

    // Adding the SPRISM variable to deactivate the quadratic interpolation
    KRATOS_REGISTER_VARIABLE(CONSIDER_QUADRATIC_SPRISM_ELEMENT);

    // Strain measures
    KRATOS_REGISTER_VARIABLE(HENCKY_STRAIN_VECTOR);
    KRATOS_REGISTER_VARIABLE(HENCKY_STRAIN_TENSOR);

    KRATOS_REGISTER_VARIABLE(VON_MISES_STRESS)

    KRATOS_REGISTER_VARIABLE(REFERENCE_DEFORMATION_GRADIENT);
    KRATOS_REGISTER_VARIABLE(REFERENCE_DEFORMATION_GRADIENT_DETERMINANT);

    // Rayleigh variables
    KRATOS_REGISTER_VARIABLE(RAYLEIGH_ALPHA)
    KRATOS_REGISTER_VARIABLE(RAYLEIGH_BETA)

    // System damping
    KRATOS_REGISTER_VARIABLE(SYSTEM_DAMPING_RATIO)
    KRATOS_REGISTER_VARIABLE(SECOND_SYSTEM_DAMPING_RATIO)

    // Nodal load variables
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(POINT_LOAD)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(LINE_LOAD)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SURFACE_LOAD)

    // Condition load variables
    KRATOS_REGISTER_VARIABLE(POINT_LOADS_VECTOR)
    KRATOS_REGISTER_VARIABLE(LINE_LOADS_VECTOR)
    KRATOS_REGISTER_VARIABLE(SURFACE_LOADS_VECTOR)
    KRATOS_REGISTER_VARIABLE(POSITIVE_FACE_PRESSURES_VECTOR)
    KRATOS_REGISTER_VARIABLE(NEGATIVE_FACE_PRESSURES_VECTOR)

    // Constitutive laws variables
    KRATOS_REGISTER_VARIABLE(YIELD_STRESS_TENSION)
    KRATOS_REGISTER_VARIABLE(PLASTIC_STRAIN_VECTOR)
    KRATOS_REGISTER_VARIABLE(YIELD_STRESS_COMPRESSION)
    KRATOS_REGISTER_VARIABLE(DILATANCY_ANGLE)
    KRATOS_REGISTER_VARIABLE(SOFTENING_TYPE)
    KRATOS_REGISTER_VARIABLE(SOFTENING_TYPE_COMPRESSION)
    KRATOS_REGISTER_VARIABLE(HARDENING_CURVE)
    KRATOS_REGISTER_VARIABLE(MAX_NUMBER_NL_CL_ITERATIONS)
    KRATOS_REGISTER_VARIABLE(VISCOUS_PARAMETER)
    KRATOS_REGISTER_VARIABLE(DELAY_TIME)
    KRATOS_REGISTER_VARIABLE(MAXIMUM_STRESS)
    KRATOS_REGISTER_VARIABLE(MAXIMUM_STRESS_POSITION)
    KRATOS_REGISTER_VARIABLE(UNIAXIAL_STRESS)
    KRATOS_REGISTER_VARIABLE(FRICTION_ANGLE)
    KRATOS_REGISTER_VARIABLE(COHESION)
    KRATOS_REGISTER_VARIABLE(DAMAGE)
    KRATOS_REGISTER_VARIABLE(THRESHOLD)
    KRATOS_REGISTER_VARIABLE(INTEGRATED_STRESS_TENSOR)
    KRATOS_REGISTER_VARIABLE(PLASTIC_STRAIN_TENSOR)
    KRATOS_REGISTER_VARIABLE(CURVE_FITTING_PARAMETERS)
    KRATOS_REGISTER_VARIABLE(PLASTIC_STRAIN_INDICATORS)
    KRATOS_REGISTER_VARIABLE(EQUIVALENT_PLASTIC_STRAIN)
    KRATOS_REGISTER_VARIABLE(KINEMATIC_PLASTICITY_PARAMETERS)
    KRATOS_REGISTER_VARIABLE(KINEMATIC_HARDENING_TYPE)
    KRATOS_REGISTER_VARIABLE(CONSIDER_PERTURBATION_THRESHOLD)
    KRATOS_REGISTER_VARIABLE(TANGENT_OPERATOR_ESTIMATION)
    KRATOS_REGISTER_VARIABLE(TENSION_STRESS_TENSOR)
    KRATOS_REGISTER_VARIABLE(COMPRESSION_STRESS_TENSOR)
    KRATOS_REGISTER_VARIABLE(TENSION_STRESS_VECTOR)
    KRATOS_REGISTER_VARIABLE(COMPRESSION_STRESS_VECTOR)
    KRATOS_REGISTER_VARIABLE(EFFECTIVE_TENSION_STRESS_VECTOR)
    KRATOS_REGISTER_VARIABLE(EFFECTIVE_COMPRESSION_STRESS_VECTOR)
    KRATOS_REGISTER_VARIABLE(EXPONENTIAL_SATURATION_YIELD_STRESS)
    KRATOS_REGISTER_VARIABLE(ACCUMULATED_PLASTIC_STRAIN)
    KRATOS_REGISTER_VARIABLE(BACK_STRESS_VECTOR)
    KRATOS_REGISTER_VARIABLE(BACK_STRESS_TENSOR)
    KRATOS_REGISTER_VARIABLE(HARDENING_MODULI_VECTOR)

    // D+D- Damage Constitutive laws variables
    KRATOS_REGISTER_VARIABLE(DAMAGE_TENSION)
    KRATOS_REGISTER_VARIABLE(DAMAGE_COMPRESSION)
    KRATOS_REGISTER_VARIABLE(THRESHOLD_TENSION)
    KRATOS_REGISTER_VARIABLE(THRESHOLD_COMPRESSION)
    KRATOS_REGISTER_VARIABLE(UNIAXIAL_STRESS_TENSION)
    KRATOS_REGISTER_VARIABLE(UNIAXIAL_STRESS_COMPRESSION)
    KRATOS_REGISTER_VARIABLE(FRACTURE_ENERGY_COMPRESSION)
    KRATOS_REGISTER_VARIABLE(FRACTURE_ENERGY_DAMAGE_PROCESS)
    KRATOS_REGISTER_VARIABLE(HIGH_CYCLE_FATIGUE_COEFFICIENTS)
    KRATOS_REGISTER_VARIABLE(FATIGUE_REDUCTION_FACTOR)
    KRATOS_REGISTER_VARIABLE(NUMBER_OF_CYCLES)
    KRATOS_REGISTER_VARIABLE(WOHLER_STRESS)

    // D+D- Damage Constitutive laws variables, additional Masonry 2D & 3D
    KRATOS_REGISTER_VARIABLE(DAMAGE_ONSET_STRESS_COMPRESSION)
    KRATOS_REGISTER_VARIABLE(BIAXIAL_COMPRESSION_MULTIPLIER)
    KRATOS_REGISTER_VARIABLE(FRACTURE_ENERGY_TENSION)
    KRATOS_REGISTER_VARIABLE(RESIDUAL_STRESS_COMPRESSION)
    KRATOS_REGISTER_VARIABLE(BEZIER_CONTROLLER_C1)
    KRATOS_REGISTER_VARIABLE(BEZIER_CONTROLLER_C2)
    KRATOS_REGISTER_VARIABLE(BEZIER_CONTROLLER_C3)
    KRATOS_REGISTER_VARIABLE(YIELD_STRAIN_COMPRESSION)
    KRATOS_REGISTER_VARIABLE(SHEAR_COMPRESSION_REDUCTOR)
    KRATOS_REGISTER_VARIABLE(TRIAXIAL_COMPRESSION_COEFFICIENT)

    // Response function variables
    KRATOS_REGISTER_VARIABLE(RESPONSE_VALUE)
    // Adjoint variables
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ADJOINT_DISPLACEMENT)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ADJOINT_ROTATION)
    KRATOS_REGISTER_VARIABLE(PERTURBATION_SIZE)
    KRATOS_REGISTER_VARIABLE(ADAPT_PERTURBATION_SIZE)

    // Variables for output of sensitivities
    KRATOS_REGISTER_VARIABLE( CROSS_AREA_SENSITIVITY );
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( POINT_LOAD_SENSITIVITY );
    KRATOS_REGISTER_VARIABLE( I22_SENSITIVITY );
    KRATOS_REGISTER_VARIABLE( I33_SENSITIVITY );
    KRATOS_REGISTER_VARIABLE( THICKNESS_SENSITIVITY );
    KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_SENSITIVITY );
    KRATOS_REGISTER_VARIABLE( AREA_EFFECTIVE_Y_SENSITIVITY );
    KRATOS_REGISTER_VARIABLE( AREA_EFFECTIVE_Z_SENSITIVITY );
    KRATOS_REGISTER_VARIABLE( IS_ADJOINT );

    // Variables to for computing parts of sensitivity analysis
    KRATOS_REGISTER_VARIABLE( TRACED_STRESS_TYPE );
    KRATOS_REGISTER_VARIABLE( STRESS_DISP_DERIV_ON_GP );
    KRATOS_REGISTER_VARIABLE( STRESS_DISP_DERIV_ON_NODE);
    KRATOS_REGISTER_VARIABLE( STRESS_DESIGN_DERIVATIVE_ON_GP );
    KRATOS_REGISTER_VARIABLE( STRESS_DESIGN_DERIVATIVE_ON_NODE);
    KRATOS_REGISTER_VARIABLE( STRESS_ON_GP  );
    KRATOS_REGISTER_VARIABLE( STRESS_ON_NODE  );
    KRATOS_REGISTER_VARIABLE( DESIGN_VARIABLE_NAME );

    // Some variables related with CL
    KRATOS_REGISTER_VARIABLE(INELASTIC_FLAG)
    KRATOS_REGISTER_VARIABLE(INFINITY_YIELD_STRESS)

    //Register the truss element
    KRATOS_REGISTER_ELEMENT("TrussElement3D2N", mTrussElement3D2N)
    KRATOS_REGISTER_ELEMENT("TrussLinearElement3D2N", mTrussLinearElement3D2N)
    KRATOS_REGISTER_ELEMENT("CableElement3D2N", mCableElement3D2N)

    // Register the beam element
    KRATOS_REGISTER_ELEMENT("CrBeamElement3D2N", mCrBeamElement3D2N)
    KRATOS_REGISTER_ELEMENT("CrLinearBeamElement3D2N", mCrLinearBeamElement3D2N)
    KRATOS_REGISTER_ELEMENT("CrBeamElement2D2N", mCrBeamElement2D2N)
    KRATOS_REGISTER_ELEMENT("CrLinearBeamElement2D2N", mCrLinearBeamElement2D2N)

    //Register the shells elements
    KRATOS_REGISTER_ELEMENT("IsotropicShellElement3D3N", mIsotropicShellElement3D3N)
    KRATOS_REGISTER_ELEMENT("ShellThickElement3D4N", mShellThickElement3D4N)
    KRATOS_REGISTER_ELEMENT("ShellThickElementCorotational3D4N", mShellThickCorotationalElement3D4N)
    KRATOS_REGISTER_ELEMENT("ShellThinElementCorotational3D4N", mShellThinCorotationalElement3D4N)
    KRATOS_REGISTER_ELEMENT("ShellThinElement3D3N", mShellThinElement3D3N)
    KRATOS_REGISTER_ELEMENT("ShellThickElementCorotational3D3N", mShellThickCorotationalElement3D3N)
    KRATOS_REGISTER_ELEMENT("ShellThinElementCorotational3D3N", mShellThinCorotationalElement3D3N)

    // Register the membrane element
    KRATOS_REGISTER_ELEMENT("PreStressMembraneElement3D3N", mPreStressMembraneElement3D3N)
    KRATOS_REGISTER_ELEMENT("PreStressMembraneElement3D4N", mPreStressMembraneElement3D4N)

    // Register the SPRISM element
    KRATOS_REGISTER_ELEMENT("SolidShellElementSprism3D6N", mSolidShellElementSprism3D6N);

    // Register the nodal concentrated element
    KRATOS_REGISTER_ELEMENT("NodalConcentratedDampedElement3D1N", mNodalConcentratedDampedElement3D1N);
    KRATOS_REGISTER_ELEMENT("NodalConcentratedElement2D1N", mNodalConcentratedElement2D1N);
    KRATOS_REGISTER_ELEMENT("NodalConcentratedDampedElement2D1N", mNodalConcentratedDampedElement2D1N);
    KRATOS_REGISTER_ELEMENT("NodalConcentratedElement3D1N", mNodalConcentratedElement3D1N);

    // SOLID ELEMENTS
    // Small displacement elements
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement2D3N", mSmallDisplacement2D3N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement2D4N", mSmallDisplacement2D4N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement2D6N", mSmallDisplacement2D6N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement2D8N", mSmallDisplacement2D8N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement2D9N", mSmallDisplacement2D9N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D4N", mSmallDisplacement3D4N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D6N", mSmallDisplacement3D6N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D8N", mSmallDisplacement3D8N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D10N", mSmallDisplacement3D10N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D15N", mSmallDisplacement3D15N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D20N", mSmallDisplacement3D20N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D27N", mSmallDisplacement3D27N)

    KRATOS_REGISTER_ELEMENT("SmallDisplacementBbarElement2D4N", mSmallDisplacementBbar2D4N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementBbarElement3D8N", mSmallDisplacementBbar3D8N)

    KRATOS_REGISTER_ELEMENT("AxisymSmallDisplacementElement2D3N", mAxisymSmallDisplacement2D3N)
    KRATOS_REGISTER_ELEMENT("AxisymSmallDisplacementElement2D4N", mAxisymSmallDisplacement2D4N)
    KRATOS_REGISTER_ELEMENT("AxisymSmallDisplacementElement2D6N", mAxisymSmallDisplacement2D6N)
    KRATOS_REGISTER_ELEMENT("AxisymSmallDisplacementElement2D8N", mAxisymSmallDisplacement2D8N)
    KRATOS_REGISTER_ELEMENT("AxisymSmallDisplacementElement2D9N", mAxisymSmallDisplacement2D9N)

    // Total lagrangian elements
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement2D3N", mTotalLagrangian2D3N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement2D4N", mTotalLagrangian2D4N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement2D6N", mTotalLagrangian2D6N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement2D8N", mTotalLagrangian2D8N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement2D9N", mTotalLagrangian2D9N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement3D4N", mTotalLagrangian3D4N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement3D6N", mTotalLagrangian3D6N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement3D8N", mTotalLagrangian3D8N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement3D10N", mTotalLagrangian3D10N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement3D15N", mTotalLagrangian3D15N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement3D20N", mTotalLagrangian3D20N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement3D27N", mTotalLagrangian3D27N)

    KRATOS_REGISTER_ELEMENT("AxisymTotalLagrangianElement2D3N", mAxisymTotalLagrangian2D3N)
    KRATOS_REGISTER_ELEMENT("AxisymTotalLagrangianElement2D4N", mAxisymTotalLagrangian2D4N)
    KRATOS_REGISTER_ELEMENT("AxisymTotalLagrangianElement2D6N", mAxisymTotalLagrangian2D6N)
    KRATOS_REGISTER_ELEMENT("AxisymTotalLagrangianElement2D8N", mAxisymTotalLagrangian2D8N)
    KRATOS_REGISTER_ELEMENT("AxisymTotalLagrangianElement2D9N", mAxisymTotalLagrangian2D9N)

    // Updated lagrangian elements
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement2D3N", mUpdatedLagrangian2D3N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement2D4N", mUpdatedLagrangian2D4N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement2D6N", mUpdatedLagrangian2D6N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement2D8N", mUpdatedLagrangian2D8N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement2D9N", mUpdatedLagrangian2D9N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement3D4N", mUpdatedLagrangian3D4N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement3D6N", mUpdatedLagrangian3D6N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement3D8N", mUpdatedLagrangian3D8N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement3D10N", mUpdatedLagrangian3D10N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement3D15N", mUpdatedLagrangian3D15N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement3D20N", mUpdatedLagrangian3D20N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement3D27N", mUpdatedLagrangian3D27N)

    KRATOS_REGISTER_ELEMENT("AxisymUpdatedLagrangianElement2D3N", mAxisymUpdatedLagrangian2D3N)
    KRATOS_REGISTER_ELEMENT("AxisymUpdatedLagrangianElement2D4N", mAxisymUpdatedLagrangian2D4N)
    KRATOS_REGISTER_ELEMENT("AxisymUpdatedLagrangianElement2D6N", mAxisymUpdatedLagrangian2D6N)
    KRATOS_REGISTER_ELEMENT("AxisymUpdatedLagrangianElement2D8N", mAxisymUpdatedLagrangian2D8N)
    KRATOS_REGISTER_ELEMENT("AxisymUpdatedLagrangianElement2D9N", mAxisymUpdatedLagrangian2D9N)

    // Register the spring damper element
    KRATOS_REGISTER_ELEMENT("SpringDamperElement3D2N", mSpringDamperElement3D2N);

    //Register the adjoint elements
    KRATOS_REGISTER_ELEMENT("AdjointFiniteDifferencingShellThinElement3D3N", mAdjointFiniteDifferencingShellThinElement3D3N )
    KRATOS_REGISTER_ELEMENT("AdjointFiniteDifferenceCrBeamElementLinear3D2N", mAdjointFiniteDifferenceCrBeamElementLinear3D2N )
    KRATOS_REGISTER_ELEMENT("AdjointFiniteDifferenceTrussElement3D2N", mAdjointFiniteDifferenceTrussElement3D2N)
    KRATOS_REGISTER_ELEMENT("AdjointFiniteDifferenceTrussLinearElement3D2N", mAdjointFiniteDifferenceTrussLinearElement3D2N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianAdjointElement2D3N", mTotalLagrangianAdjoint2D3N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianAdjointElement2D4N", mTotalLagrangianAdjoint2D4N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianAdjointElement2D6N", mTotalLagrangianAdjoint2D6N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianAdjointElement3D4N", mTotalLagrangianAdjoint3D4N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianAdjointElement3D8N", mTotalLagrangianAdjoint3D8N)
    KRATOS_REGISTER_ELEMENT("AdjointFiniteDifferencingSmallDisplacementElement3D4N", mAdjointFiniteDifferencingSmallDisplacementElement3D4N)
    KRATOS_REGISTER_ELEMENT("AdjointFiniteDifferencingSmallDisplacementElement3D6N", mAdjointFiniteDifferencingSmallDisplacementElement3D6N)
    KRATOS_REGISTER_ELEMENT("AdjointFiniteDifferencingSmallDisplacementElement3D8N", mAdjointFiniteDifferencingSmallDisplacementElement3D8N)

    // Register the conditions
    // Point loads
    KRATOS_REGISTER_CONDITION("PointLoadCondition2D1N", mPointLoadCondition2D1N)
    KRATOS_REGISTER_CONDITION("PointLoadCondition3D1N", mPointLoadCondition3D1N)
    KRATOS_REGISTER_CONDITION("PointContactCondition2D1N", mPointContactCondition2D1N)
    KRATOS_REGISTER_CONDITION("PointContactCondition3D1N", mPointContactCondition3D1N)

    KRATOS_REGISTER_CONDITION("AxisymPointLoadCondition2D1N", mAxisymPointLoadCondition2D1N)

    // Line loads
    KRATOS_REGISTER_CONDITION("LineLoadCondition2D2N", mLineLoadCondition2D2N)
    KRATOS_REGISTER_CONDITION("LineLoadCondition2D3N", mLineLoadCondition2D3N)
    KRATOS_REGISTER_CONDITION("LineLoadCondition3D2N", mLineLoadCondition3D2N)
    KRATOS_REGISTER_CONDITION("LineLoadCondition3D3N", mLineLoadCondition3D3N)

    KRATOS_REGISTER_CONDITION("AxisymLineLoadCondition2D2N", mAxisymLineLoadCondition2D2N)
    KRATOS_REGISTER_CONDITION("AxisymLineLoadCondition2D3N", mAxisymLineLoadCondition2D3N)

    // Surface loads
    KRATOS_REGISTER_CONDITION("SurfaceLoadCondition3D3N", mSurfaceLoadCondition3D3N)
    KRATOS_REGISTER_CONDITION("SurfaceLoadCondition3D4N", mSurfaceLoadCondition3D4N)
    KRATOS_REGISTER_CONDITION("SurfaceLoadCondition3D6N", mSurfaceLoadCondition3D6N)
    KRATOS_REGISTER_CONDITION("SurfaceLoadCondition3D8N", mSurfaceLoadCondition3D8N)
    KRATOS_REGISTER_CONDITION("SurfaceLoadCondition3D9N", mSurfaceLoadCondition3D9N)

    // Point moment
    KRATOS_REGISTER_CONDITION("PointMomentCondition3D1N", mPointMomentCondition3D1N);

    // Adjoint conditions
    KRATOS_REGISTER_CONDITION("AdjointSemiAnalyticPointLoadCondition2D1N", mAdjointSemiAnalyticPointLoadCondition2D1N )
    KRATOS_REGISTER_CONDITION("AdjointSemiAnalyticPointLoadCondition3D1N", mAdjointSemiAnalyticPointLoadCondition3D1N )

    // Register linear elastics laws
    KRATOS_REGISTER_CONSTITUTIVE_LAW("TrussConstitutiveLaw", mTrussConstitutiveLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("TrussPlasticityConstitutiveLaw", mTrussPlasticityConstitutiveLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("BeamConstitutiveLaw", mBeamConstitutiveLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElastic3DLaw", mElasticIsotropic3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticPlaneStrain2DLaw", mLinearPlaneStrain);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticPlaneStress2DLaw", mLinearPlaneStress);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ElasticPlaneStressUncoupledShear2DLaw", mElasticIsotropicPlaneStressUncoupledShear);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticAxisym2DLaw", mAxisymElasticIsotropic);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticOrthotropic2DLaw", mLinearElasticOrthotropic2DLaw);
    // Register hyper elastic laws
    KRATOS_REGISTER_CONSTITUTIVE_LAW("KirchhoffSaintVenant3DLaw", mHyperElasticIsotropicKirchhoff3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("KirchhoffSaintVenantPlaneStress2DLaw", mHyperElasticIsotropicKirchhoffPlaneStress2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("KirchhoffSaintVenantPlaneStrain2DLaw", mHyperElasticIsotropicKirchhoffPlaneStrain2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElastic3DLaw", mHyperElasticIsotropicNeoHookean3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticPlaneStrain2DLaw", mHyperElasticIsotropicNeoHookeanPlaneStrain2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainJ2PlasticityPlaneStrain2DLaw", mSmallStrainJ2PlasticityPlaneStrain2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainJ2Plasticity3DLaw", mSmallStrainJ2Plasticity3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamagePlaneStrain2DLaw", mSmallStrainIsotropicDamagePlaneStrain2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DLaw", mSmallStrainIsotropicDamage3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamageTractionOnly3DLaw", mSmallStrainIsotropicDamageTractionOnly3D);

    // Damage and plasticity
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticityFactory", mSmallStrainIsotropicPlasticityFactory);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticityFactory", mSmallStrainKinematicPlasticityFactory);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamageFactory", mSmallStrainIsotropicDamageFactory);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ViscousGeneralizedKelvin3D", mViscousGeneralizedKelvin3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ViscousGeneralizedMaxwell3D", mViscousGeneralizedMaxwell3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("GenericSmallStrainViscoplasticity3D", mGenericSmallStrainViscoplasticity3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("PlasticityIsotropicKinematicJ2Law", mPlasticityIsotropicKinematicJ2);

    // Custom Constitutive laws
    /// Plasticity

    /* Small strain */
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DVonMisesVonMises", mSmallStrainIsotropicPlasticity3DVonMisesVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DVonMisesModifiedMohrCoulomb", mSmallStrainIsotropicPlasticity3DVonMisesModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DVonMisesDruckerPrager", mSmallStrainIsotropicPlasticity3DVonMisesDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DVonMisesTresca", mSmallStrainIsotropicPlasticity3DVonMisesTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DModifiedMohrCoulombVonMises", mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb", mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DModifiedMohrCoulombDruckerPrager", mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DModifiedMohrCoulombTresca", mSmallStrainIsotropicPlasticity3DModifiedMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DTrescaVonMises", mSmallStrainIsotropicPlasticity3DTrescaVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DTrescaModifiedMohrCoulomb", mSmallStrainIsotropicPlasticity3DTrescaModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DTrescaDruckerPrager", mSmallStrainIsotropicPlasticity3DTrescaDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DTrescaTresca", mSmallStrainIsotropicPlasticity3DTrescaTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DDruckerPragerVonMises", mSmallStrainIsotropicPlasticity3DDruckerPragerVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DDruckerPragerModifiedMohrCoulomb", mSmallStrainIsotropicPlasticity3DDruckerPragerModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DDruckerPragerDruckerPrager", mSmallStrainIsotropicPlasticity3DDruckerPragerDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DDruckerPragerTresca", mSmallStrainIsotropicPlasticity3DDruckerPragerTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DVonMisesMohrCoulomb", mSmallStrainIsotropicPlasticity3DVonMisesMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DMohrCoulombVonMises", mSmallStrainIsotropicPlasticity3DMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DMohrCoulombMohrCoulomb", mSmallStrainIsotropicPlasticity3DMohrCoulombMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DMohrCoulombDruckerPrager", mSmallStrainIsotropicPlasticity3DMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DMohrCoulombTresca", mSmallStrainIsotropicPlasticity3DMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DTrescaMohrCoulomb", mSmallStrainIsotropicPlasticity3DTrescaMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicPlasticity3DDruckerPragerMohrCoulomb", mSmallStrainIsotropicPlasticity3DDruckerPragerMohrCoulomb);

    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DVonMisesVonMises", mSmallStrainKinematicPlasticity3DVonMisesVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DVonMisesModifiedMohrCoulomb", mSmallStrainKinematicPlasticity3DVonMisesModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DVonMisesDruckerPrager", mSmallStrainKinematicPlasticity3DVonMisesDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DVonMisesTresca", mSmallStrainKinematicPlasticity3DVonMisesTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DModifiedMohrCoulombVonMises", mSmallStrainKinematicPlasticity3DModifiedMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb", mSmallStrainKinematicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DModifiedMohrCoulombDruckerPrager", mSmallStrainKinematicPlasticity3DModifiedMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DModifiedMohrCoulombTresca", mSmallStrainKinematicPlasticity3DModifiedMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DTrescaVonMises", mSmallStrainKinematicPlasticity3DTrescaVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DTrescaModifiedMohrCoulomb", mSmallStrainKinematicPlasticity3DTrescaModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DTrescaDruckerPrager", mSmallStrainKinematicPlasticity3DTrescaDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DTrescaTresca", mSmallStrainKinematicPlasticity3DTrescaTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DDruckerPragerVonMises", mSmallStrainKinematicPlasticity3DDruckerPragerVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DDruckerPragerModifiedMohrCoulomb", mSmallStrainKinematicPlasticity3DDruckerPragerModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DDruckerPragerDruckerPrager", mSmallStrainKinematicPlasticity3DDruckerPragerDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DDruckerPragerTresca", mSmallStrainKinematicPlasticity3DDruckerPragerTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DVonMisesMohrCoulomb", mSmallStrainKinematicPlasticity3DVonMisesMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DMohrCoulombVonMises", mSmallStrainKinematicPlasticity3DMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DMohrCoulombMohrCoulomb", mSmallStrainKinematicPlasticity3DMohrCoulombMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DMohrCoulombDruckerPrager", mSmallStrainKinematicPlasticity3DMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DMohrCoulombTresca", mSmallStrainKinematicPlasticity3DMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DTrescaMohrCoulomb", mSmallStrainKinematicPlasticity3DTrescaMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainKinematicPlasticity3DDruckerPragerMohrCoulomb", mSmallStrainKinematicPlasticity3DDruckerPragerMohrCoulomb);

    //Plastic Damage Model
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainPlasticDamageModel3DVonMisesVonMisesVonMises", mSmallStrainPlasticDamageModel3DVonMisesVonMisesVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainPlasticDamageModel3DVonMisesVonMisesDruckerPrager", mSmallStrainPlasticDamageModel3DVonMisesVonMisesDruckerPrager);


    /* Finite strain */

    // Kirchhoff
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DVonMisesVonMises", mHyperElasticIsotropicKirchhoffPlasticity3DVonMisesVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DVonMisesModifiedMohrCoulomb", mHyperElasticIsotropicKirchhoffPlasticity3DVonMisesModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DVonMisesDruckerPrager", mHyperElasticIsotropicKirchhoffPlasticity3DVonMisesDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DVonMisesTresca", mHyperElasticIsotropicKirchhoffPlasticity3DVonMisesTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombVonMises", mHyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb", mHyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombDruckerPrager", mHyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombTresca", mHyperElasticIsotropicKirchhoffPlasticity3DModifiedMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DTrescaVonMises", mHyperElasticIsotropicKirchhoffPlasticity3DTrescaVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DTrescaModifiedMohrCoulomb", mHyperElasticIsotropicKirchhoffPlasticity3DTrescaModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DTrescaDruckerPrager", mHyperElasticIsotropicKirchhoffPlasticity3DTrescaDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DTrescaTresca", mHyperElasticIsotropicKirchhoffPlasticity3DTrescaTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerVonMises", mHyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerModifiedMohrCoulomb", mHyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerDruckerPrager", mHyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerTresca", mHyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DVonMisesMohrCoulomb", mHyperElasticIsotropicKirchhoffPlasticity3DVonMisesMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DMohrCoulombVonMises", mHyperElasticIsotropicKirchhoffPlasticity3DMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DMohrCoulombMohrCoulomb", mHyperElasticIsotropicKirchhoffPlasticity3DMohrCoulombMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DMohrCoulombDruckerPrager", mHyperElasticIsotropicKirchhoffPlasticity3DMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DMohrCoulombTresca", mHyperElasticIsotropicKirchhoffPlasticity3DMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DTrescaMohrCoulomb", mHyperElasticIsotropicKirchhoffPlasticity3DTrescaMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerMohrCoulomb", mHyperElasticIsotropicKirchhoffPlasticity3DDruckerPragerMohrCoulomb);

    // Neo-Hookean
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DVonMisesVonMises", mHyperElasticIsotropicNeoHookeanPlasticity3DVonMisesVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DVonMisesModifiedMohrCoulomb", mHyperElasticIsotropicNeoHookeanPlasticity3DVonMisesModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DVonMisesDruckerPrager", mHyperElasticIsotropicNeoHookeanPlasticity3DVonMisesDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DVonMisesTresca", mHyperElasticIsotropicNeoHookeanPlasticity3DVonMisesTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombVonMises", mHyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb", mHyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombDruckerPrager", mHyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombTresca", mHyperElasticIsotropicNeoHookeanPlasticity3DModifiedMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DTrescaVonMises", mHyperElasticIsotropicNeoHookeanPlasticity3DTrescaVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DTrescaModifiedMohrCoulomb", mHyperElasticIsotropicNeoHookeanPlasticity3DTrescaModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DTrescaDruckerPrager", mHyperElasticIsotropicNeoHookeanPlasticity3DTrescaDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DTrescaTresca", mHyperElasticIsotropicNeoHookeanPlasticity3DTrescaTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerVonMises", mHyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerModifiedMohrCoulomb", mHyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerDruckerPrager", mHyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerTresca", mHyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DVonMisesMohrCoulomb", mHyperElasticIsotropicNeoHookeanPlasticity3DVonMisesMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DMohrCoulombVonMises", mHyperElasticIsotropicNeoHookeanPlasticity3DMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DMohrCoulombMohrCoulomb", mHyperElasticIsotropicNeoHookeanPlasticity3DMohrCoulombMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DMohrCoulombDruckerPrager", mHyperElasticIsotropicNeoHookeanPlasticity3DMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DMohrCoulombTresca", mHyperElasticIsotropicNeoHookeanPlasticity3DMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DTrescaMohrCoulomb", mHyperElasticIsotropicNeoHookeanPlasticity3DTrescaMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerMohrCoulomb", mHyperElasticIsotropicNeoHookeanPlasticity3DDruckerPragerMohrCoulomb);

    // Kirchhoff
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DVonMisesVonMises", mHyperElasticKirchhoffKinematicPlasticity3DVonMisesVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DVonMisesModifiedMohrCoulomb", mHyperElasticKirchhoffKinematicPlasticity3DVonMisesModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DVonMisesDruckerPrager", mHyperElasticKirchhoffKinematicPlasticity3DVonMisesDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DVonMisesTresca", mHyperElasticKirchhoffKinematicPlasticity3DVonMisesTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DModifiedMohrCoulombVonMises", mHyperElasticKirchhoffKinematicPlasticity3DModifiedMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb", mHyperElasticKirchhoffKinematicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DModifiedMohrCoulombDruckerPrager", mHyperElasticKirchhoffKinematicPlasticity3DModifiedMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DModifiedMohrCoulombTresca", mHyperElasticKirchhoffKinematicPlasticity3DModifiedMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DTrescaVonMises", mHyperElasticKirchhoffKinematicPlasticity3DTrescaVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DTrescaModifiedMohrCoulomb", mHyperElasticKirchhoffKinematicPlasticity3DTrescaModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DTrescaDruckerPrager", mHyperElasticKirchhoffKinematicPlasticity3DTrescaDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DTrescaTresca", mHyperElasticKirchhoffKinematicPlasticity3DTrescaTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DDruckerPragerVonMises", mHyperElasticKirchhoffKinematicPlasticity3DDruckerPragerVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DDruckerPragerModifiedMohrCoulomb", mHyperElasticKirchhoffKinematicPlasticity3DDruckerPragerModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DDruckerPragerDruckerPrager", mHyperElasticKirchhoffKinematicPlasticity3DDruckerPragerDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DDruckerPragerTresca", mHyperElasticKirchhoffKinematicPlasticity3DDruckerPragerTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DVonMisesMohrCoulomb", mHyperElasticKirchhoffKinematicPlasticity3DVonMisesMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DMohrCoulombVonMises", mHyperElasticKirchhoffKinematicPlasticity3DMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DMohrCoulombMohrCoulomb", mHyperElasticKirchhoffKinematicPlasticity3DMohrCoulombMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DMohrCoulombDruckerPrager", mHyperElasticKirchhoffKinematicPlasticity3DMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DMohrCoulombTresca", mHyperElasticKirchhoffKinematicPlasticity3DMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DTrescaMohrCoulomb", mHyperElasticKirchhoffKinematicPlasticity3DTrescaMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticKirchhoffKinematicPlasticity3DDruckerPragerMohrCoulomb", mHyperElasticKirchhoffKinematicPlasticity3DDruckerPragerMohrCoulomb);

    // Neo-Hookean
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DVonMisesVonMises", mHyperElasticNeoHookeanKinematicPlasticity3DVonMisesVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DVonMisesModifiedMohrCoulomb", mHyperElasticNeoHookeanKinematicPlasticity3DVonMisesModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DVonMisesDruckerPrager", mHyperElasticNeoHookeanKinematicPlasticity3DVonMisesDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DVonMisesTresca", mHyperElasticNeoHookeanKinematicPlasticity3DVonMisesTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DModifiedMohrCoulombVonMises", mHyperElasticNeoHookeanKinematicPlasticity3DModifiedMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb", mHyperElasticNeoHookeanKinematicPlasticity3DModifiedMohrCoulombModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DModifiedMohrCoulombDruckerPrager", mHyperElasticNeoHookeanKinematicPlasticity3DModifiedMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DModifiedMohrCoulombTresca", mHyperElasticNeoHookeanKinematicPlasticity3DModifiedMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DTrescaVonMises", mHyperElasticNeoHookeanKinematicPlasticity3DTrescaVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DTrescaModifiedMohrCoulomb", mHyperElasticNeoHookeanKinematicPlasticity3DTrescaModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DTrescaDruckerPrager", mHyperElasticNeoHookeanKinematicPlasticity3DTrescaDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DTrescaTresca", mHyperElasticNeoHookeanKinematicPlasticity3DTrescaTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DDruckerPragerVonMises", mHyperElasticNeoHookeanKinematicPlasticity3DDruckerPragerVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DDruckerPragerModifiedMohrCoulomb", mHyperElasticNeoHookeanKinematicPlasticity3DDruckerPragerModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DDruckerPragerDruckerPrager", mHyperElasticNeoHookeanKinematicPlasticity3DDruckerPragerDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DDruckerPragerTresca", mHyperElasticNeoHookeanKinematicPlasticity3DDruckerPragerTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DVonMisesMohrCoulomb", mHyperElasticNeoHookeanKinematicPlasticity3DVonMisesMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DMohrCoulombVonMises", mHyperElasticNeoHookeanKinematicPlasticity3DMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DMohrCoulombMohrCoulomb", mHyperElasticNeoHookeanKinematicPlasticity3DMohrCoulombMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DMohrCoulombDruckerPrager", mHyperElasticNeoHookeanKinematicPlasticity3DMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DMohrCoulombTresca", mHyperElasticNeoHookeanKinematicPlasticity3DMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DTrescaMohrCoulomb", mHyperElasticNeoHookeanKinematicPlasticity3DTrescaMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanKinematicPlasticity3DDruckerPragerMohrCoulomb", mHyperElasticNeoHookeanKinematicPlasticity3DDruckerPragerMohrCoulomb);

    /// Damage

    /* Small strain */
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DVonMisesVonMises", mSmallStrainIsotropicDamage3DVonMisesVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DVonMisesModifiedMohrCoulomb", mSmallStrainIsotropicDamage3DVonMisesModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DVonMisesDruckerPrager", mSmallStrainIsotropicDamage3DVonMisesDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DVonMisesTresca", mSmallStrainIsotropicDamage3DVonMisesTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DModifiedMohrCoulombVonMises", mSmallStrainIsotropicDamage3DModifiedMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DModifiedMohrCoulombModifiedMohrCoulomb", mSmallStrainIsotropicDamage3DModifiedMohrCoulombModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DModifiedMohrCoulombDruckerPrager", mSmallStrainIsotropicDamage3DModifiedMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DModifiedMohrCoulombTresca", mSmallStrainIsotropicDamage3DModifiedMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DTrescaVonMises", mSmallStrainIsotropicDamage3DTrescaVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DTrescaModifiedMohrCoulomb", mSmallStrainIsotropicDamage3DTrescaModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DTrescaDruckerPrager", mSmallStrainIsotropicDamage3DTrescaDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DTrescaTresca", mSmallStrainIsotropicDamage3DTrescaTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DDruckerPragerVonMises", mSmallStrainIsotropicDamage3DDruckerPragerVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DDruckerPragerModifiedMohrCoulomb", mSmallStrainIsotropicDamage3DDruckerPragerModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DDruckerPragerDruckerPrager", mSmallStrainIsotropicDamage3DDruckerPragerDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DDruckerPragerTresca", mSmallStrainIsotropicDamage3DDruckerPragerTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DRankineVonMises", mSmallStrainIsotropicDamage3DRankineVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DRankineModifiedMohrCoulomb", mSmallStrainIsotropicDamage3DRankineModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DRankineDruckerPrager", mSmallStrainIsotropicDamage3DRankineDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DRankineTresca", mSmallStrainIsotropicDamage3DRankineTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DSimoJuVonMises", mSmallStrainIsotropicDamage3DSimoJuVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DSimoJuModifiedMohrCoulomb", mSmallStrainIsotropicDamage3DSimoJuModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DSimoJuDruckerPrager", mSmallStrainIsotropicDamage3DSimoJuDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DSimoJuTresca", mSmallStrainIsotropicDamage3DSimoJuTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DVonMisesMohrCoulomb", mSmallStrainIsotropicDamage3DVonMisesMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DMohrCoulombVonMises", mSmallStrainIsotropicDamage3DMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DMohrCoulombMohrCoulomb", mSmallStrainIsotropicDamage3DMohrCoulombMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DMohrCoulombDruckerPrager", mSmallStrainIsotropicDamage3DMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DMohrCoulombTresca", mSmallStrainIsotropicDamage3DMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DTrescaMohrCoulomb", mSmallStrainIsotropicDamage3DTrescaMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DDruckerPragerMohrCoulomb", mSmallStrainIsotropicDamage3DDruckerPragerMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DRankineMohrCoulomb", mSmallStrainIsotropicDamage3DRankineMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage3DSimoJuMohrCoulomb", mSmallStrainIsotropicDamage3DSimoJuMohrCoulomb);

    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DVonMisesVonMises", mSmallStrainIsotropicDamage2DVonMisesVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DVonMisesModifiedMohrCoulomb", mSmallStrainIsotropicDamage2DVonMisesModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DVonMisesDruckerPrager", mSmallStrainIsotropicDamage2DVonMisesDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DVonMisesTresca", mSmallStrainIsotropicDamage2DVonMisesTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DModifiedMohrCoulombVonMises", mSmallStrainIsotropicDamage2DModifiedMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DModifiedMohrCoulombModifiedMohrCoulomb", mSmallStrainIsotropicDamage2DModifiedMohrCoulombModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DModifiedMohrCoulombDruckerPrager", mSmallStrainIsotropicDamage2DModifiedMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DModifiedMohrCoulombTresca", mSmallStrainIsotropicDamage2DModifiedMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DTrescaVonMises", mSmallStrainIsotropicDamage2DTrescaVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DTrescaModifiedMohrCoulomb", mSmallStrainIsotropicDamage2DTrescaModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DTrescaDruckerPrager", mSmallStrainIsotropicDamage2DTrescaDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DTrescaTresca", mSmallStrainIsotropicDamage2DTrescaTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DDruckerPragerVonMises", mSmallStrainIsotropicDamage2DDruckerPragerVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DDruckerPragerModifiedMohrCoulomb", mSmallStrainIsotropicDamage2DDruckerPragerModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DDruckerPragerDruckerPrager", mSmallStrainIsotropicDamage2DDruckerPragerDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DDruckerPragerTresca", mSmallStrainIsotropicDamage2DDruckerPragerTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DRankineVonMises", mSmallStrainIsotropicDamage2DRankineVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DRankineModifiedMohrCoulomb", mSmallStrainIsotropicDamage2DRankineModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DRankineDruckerPrager", mSmallStrainIsotropicDamage2DRankineDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DRankineTresca", mSmallStrainIsotropicDamage2DRankineTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DSimoJuVonMises", mSmallStrainIsotropicDamage2DSimoJuVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DSimoJuModifiedMohrCoulomb", mSmallStrainIsotropicDamage2DSimoJuModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DSimoJuDruckerPrager", mSmallStrainIsotropicDamage2DSimoJuDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DSimoJuTresca", mSmallStrainIsotropicDamage2DSimoJuTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DVonMisesMohrCoulomb", mSmallStrainIsotropicDamage2DVonMisesMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DMohrCoulombVonMises", mSmallStrainIsotropicDamage2DMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DMohrCoulombMohrCoulomb", mSmallStrainIsotropicDamage2DMohrCoulombMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DMohrCoulombDruckerPrager", mSmallStrainIsotropicDamage2DMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DMohrCoulombTresca", mSmallStrainIsotropicDamage2DMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DTrescaMohrCoulomb", mSmallStrainIsotropicDamage2DTrescaMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DDruckerPragerMohrCoulomb", mSmallStrainIsotropicDamage2DDruckerPragerMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DRankineMohrCoulomb", mSmallStrainIsotropicDamage2DRankineMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainIsotropicDamage2DSimoJuMohrCoulomb", mSmallStrainIsotropicDamage2DSimoJuMohrCoulomb);

    // HCF (High Cycle Fatigue)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawVonMisesVonMises", mSmallStrainHighCycleFatigue3DLawVonMisesVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawVonMisesModifiedMohrCoulomb", mSmallStrainHighCycleFatigue3DLawVonMisesModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawVonMisesDruckerPrager", mSmallStrainHighCycleFatigue3DLawVonMisesDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawVonMisesTresca", mSmallStrainHighCycleFatigue3DLawVonMisesTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawModifiedMohrCoulombVonMises", mSmallStrainHighCycleFatigue3DLawModifiedMohrCoulombVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawModifiedMohrCoulombModifiedMohrCoulomb", mSmallStrainHighCycleFatigue3DLawModifiedMohrCoulombModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawModifiedMohrCoulombDruckerPrager", mSmallStrainHighCycleFatigue3DLawModifiedMohrCoulombDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawModifiedMohrCoulombTresca", mSmallStrainHighCycleFatigue3DLawModifiedMohrCoulombTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawTrescaVonMises", mSmallStrainHighCycleFatigue3DLawTrescaVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawTrescaModifiedMohrCoulomb", mSmallStrainHighCycleFatigue3DLawTrescaModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawTrescaDruckerPrager", mSmallStrainHighCycleFatigue3DLawTrescaDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawTrescaTresca", mSmallStrainHighCycleFatigue3DLawTrescaTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawDruckerPragerVonMises", mSmallStrainHighCycleFatigue3DLawDruckerPragerVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawDruckerPragerModifiedMohrCoulomb", mSmallStrainHighCycleFatigue3DLawDruckerPragerModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawDruckerPragerDruckerPrager", mSmallStrainHighCycleFatigue3DLawDruckerPragerDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawDruckerPragerTresca", mSmallStrainHighCycleFatigue3DLawDruckerPragerTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawRankineVonMises", mSmallStrainHighCycleFatigue3DLawRankineVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawRankineModifiedMohrCoulomb", mSmallStrainHighCycleFatigue3DLawRankineModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawRankineDruckerPrager", mSmallStrainHighCycleFatigue3DLawRankineDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawRankineTresca", mSmallStrainHighCycleFatigue3DLawRankineTresca);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawSimoJuVonMises", mSmallStrainHighCycleFatigue3DLawSimoJuVonMises);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawSimoJuModifiedMohrCoulomb", mSmallStrainHighCycleFatigue3DLawSimoJuModifiedMohrCoulomb);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawSimoJuDruckerPrager", mSmallStrainHighCycleFatigue3DLawSimoJuDruckerPrager);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainHighCycleFatigue3DLawSimoJuTresca", mSmallStrainHighCycleFatigue3DLawSimoJuTresca);

    // d+d- laws
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombModifiedMohrCoulomb3D", mSmallStrainDplusDminusDamageModifiedMohrCoulombModifiedMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombRankine3D", mSmallStrainDplusDminusDamageModifiedMohrCoulombRankine3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombSimoJu3D", mSmallStrainDplusDminusDamageModifiedMohrCoulombSimoJu3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombVonMises3D", mSmallStrainDplusDminusDamageModifiedMohrCoulombVonMises3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombTresca3D", mSmallStrainDplusDminusDamageModifiedMohrCoulombTresca3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombDruckerPrager3D", mSmallStrainDplusDminusDamageModifiedMohrCoulombDruckerPrager3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineModifiedMohrCoulomb3D", mSmallStrainDplusDminusDamageRankineModifiedMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineRankine3D", mSmallStrainDplusDminusDamageRankineRankine3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineSimoJu3D", mSmallStrainDplusDminusDamageRankineSimoJu3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineVonMises3D", mSmallStrainDplusDminusDamageRankineVonMises3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineTresca3D", mSmallStrainDplusDminusDamageRankineTresca3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineDruckerPrager3D", mSmallStrainDplusDminusDamageRankineDruckerPrager3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuModifiedMohrCoulomb3D", mSmallStrainDplusDminusDamageSimoJuModifiedMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuRankine3D", mSmallStrainDplusDminusDamageSimoJuRankine3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuSimoJu3D", mSmallStrainDplusDminusDamageSimoJuSimoJu3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuVonMises3D", mSmallStrainDplusDminusDamageSimoJuVonMises3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuTresca3D", mSmallStrainDplusDminusDamageSimoJuTresca3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuDruckerPrager3D", mSmallStrainDplusDminusDamageSimoJuDruckerPrager3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesModifiedMohrCoulomb3D", mSmallStrainDplusDminusDamageVonMisesModifiedMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesRankine3D", mSmallStrainDplusDminusDamageVonMisesRankine3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesSimoJu3D", mSmallStrainDplusDminusDamageVonMisesSimoJu3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesVonMises3D", mSmallStrainDplusDminusDamageVonMisesVonMises3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesTresca3D", mSmallStrainDplusDminusDamageVonMisesTresca3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesDruckerPrager3D", mSmallStrainDplusDminusDamageVonMisesDruckerPrager3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaModifiedMohrCoulomb3D", mSmallStrainDplusDminusDamageTrescaModifiedMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaRankine3D", mSmallStrainDplusDminusDamageTrescaRankine3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaSimoJu3D", mSmallStrainDplusDminusDamageTrescaSimoJu3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaVonMises3D", mSmallStrainDplusDminusDamageTrescaVonMises3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaTresca3D", mSmallStrainDplusDminusDamageTrescaTresca3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaDruckerPrager3D", mSmallStrainDplusDminusDamageTrescaDruckerPrager3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerModifiedMohrCoulomb3D", mSmallStrainDplusDminusDamageDruckerPragerModifiedMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerRankine3D", mSmallStrainDplusDminusDamageDruckerPragerRankine3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerSimoJu3D", mSmallStrainDplusDminusDamageDruckerPragerSimoJu3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerVonMises3D", mSmallStrainDplusDminusDamageDruckerPragerVonMises3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerTresca3D", mSmallStrainDplusDminusDamageDruckerPragerTresca3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerDruckerPrager3D", mSmallStrainDplusDminusDamageDruckerPragerDruckerPrager3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombMohrCoulomb3D", mSmallStrainDplusDminusDamageMohrCoulombMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombRankine3D", mSmallStrainDplusDminusDamageMohrCoulombRankine3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombSimoJu3D", mSmallStrainDplusDminusDamageMohrCoulombSimoJu3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombVonMises3D", mSmallStrainDplusDminusDamageMohrCoulombVonMises3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombTresca3D", mSmallStrainDplusDminusDamageMohrCoulombTresca3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombDruckerPrager3D", mSmallStrainDplusDminusDamageMohrCoulombDruckerPrager3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineMohrCoulomb3D", mSmallStrainDplusDminusDamageRankineMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuMohrCoulomb3D", mSmallStrainDplusDminusDamageSimoJuMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesMohrCoulomb3D", mSmallStrainDplusDminusDamageVonMisesMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaMohrCoulomb3D", mSmallStrainDplusDminusDamageTrescaMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerMohrCoulomb3D", mSmallStrainDplusDminusDamageDruckerPragerMohrCoulomb3D);

    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombModifiedMohrCoulomb2D", mSmallStrainDplusDminusDamageModifiedMohrCoulombModifiedMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombRankine2D", mSmallStrainDplusDminusDamageModifiedMohrCoulombRankine2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombSimoJu2D", mSmallStrainDplusDminusDamageModifiedMohrCoulombSimoJu2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombVonMises2D", mSmallStrainDplusDminusDamageModifiedMohrCoulombVonMises2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombTresca2D", mSmallStrainDplusDminusDamageModifiedMohrCoulombTresca2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageModifiedMohrCoulombDruckerPrager2D", mSmallStrainDplusDminusDamageModifiedMohrCoulombDruckerPrager2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineModifiedMohrCoulomb2D", mSmallStrainDplusDminusDamageRankineModifiedMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineRankine2D", mSmallStrainDplusDminusDamageRankineRankine2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineSimoJu2D", mSmallStrainDplusDminusDamageRankineSimoJu2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineVonMises2D", mSmallStrainDplusDminusDamageRankineVonMises2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineTresca2D", mSmallStrainDplusDminusDamageRankineTresca2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineDruckerPrager2D", mSmallStrainDplusDminusDamageRankineDruckerPrager2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuModifiedMohrCoulomb2D", mSmallStrainDplusDminusDamageSimoJuModifiedMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuRankine2D", mSmallStrainDplusDminusDamageSimoJuRankine2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuSimoJu2D", mSmallStrainDplusDminusDamageSimoJuSimoJu2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuVonMises2D", mSmallStrainDplusDminusDamageSimoJuVonMises2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuTresca2D", mSmallStrainDplusDminusDamageSimoJuTresca2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuDruckerPrager2D", mSmallStrainDplusDminusDamageSimoJuDruckerPrager2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesModifiedMohrCoulomb2D", mSmallStrainDplusDminusDamageVonMisesModifiedMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesRankine2D", mSmallStrainDplusDminusDamageVonMisesRankine2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesSimoJu2D", mSmallStrainDplusDminusDamageVonMisesSimoJu2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesVonMises2D", mSmallStrainDplusDminusDamageVonMisesVonMises2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesTresca2D", mSmallStrainDplusDminusDamageVonMisesTresca2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesDruckerPrager2D", mSmallStrainDplusDminusDamageVonMisesDruckerPrager2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaModifiedMohrCoulomb2D", mSmallStrainDplusDminusDamageTrescaModifiedMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaRankine2D", mSmallStrainDplusDminusDamageTrescaRankine2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaSimoJu2D", mSmallStrainDplusDminusDamageTrescaSimoJu2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaVonMises2D", mSmallStrainDplusDminusDamageTrescaVonMises2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaTresca2D", mSmallStrainDplusDminusDamageTrescaTresca2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaDruckerPrager2D", mSmallStrainDplusDminusDamageTrescaDruckerPrager2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerModifiedMohrCoulomb2D", mSmallStrainDplusDminusDamageDruckerPragerModifiedMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerRankine2D", mSmallStrainDplusDminusDamageDruckerPragerRankine2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerSimoJu2D", mSmallStrainDplusDminusDamageDruckerPragerSimoJu2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerVonMises2D", mSmallStrainDplusDminusDamageDruckerPragerVonMises2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerTresca2D", mSmallStrainDplusDminusDamageDruckerPragerTresca2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerDruckerPrager2D", mSmallStrainDplusDminusDamageDruckerPragerDruckerPrager2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombMohrCoulomb2D", mSmallStrainDplusDminusDamageMohrCoulombMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombRankine2D", mSmallStrainDplusDminusDamageMohrCoulombRankine2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombSimoJu2D", mSmallStrainDplusDminusDamageMohrCoulombSimoJu2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombVonMises2D", mSmallStrainDplusDminusDamageMohrCoulombVonMises2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombTresca2D", mSmallStrainDplusDminusDamageMohrCoulombTresca2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageMohrCoulombDruckerPrager2D", mSmallStrainDplusDminusDamageMohrCoulombDruckerPrager2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageRankineMohrCoulomb2D", mSmallStrainDplusDminusDamageRankineMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageSimoJuMohrCoulomb2D", mSmallStrainDplusDminusDamageSimoJuMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageVonMisesMohrCoulomb2D", mSmallStrainDplusDminusDamageVonMisesMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageTrescaMohrCoulomb2D", mSmallStrainDplusDminusDamageTrescaMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainDplusDminusDamageDruckerPragerMohrCoulomb2D", mSmallStrainDplusDminusDamageDruckerPragerMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("DamageDPlusDMinusPlaneStressMasonry2DLaw", mDamageDPlusDMinusPlaneStressMasonry2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("DamageDPlusDMinusMasonry3DLaw", mDamageDPlusDMinusMasonry3DLaw);

    // Orthotropic Damage
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageRankine3D", mSmallStrainOrthotropicDamageRankine3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageVonMises3D", mSmallStrainOrthotropicDamageVonMises3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageDruckerPrager3D", mSmallStrainOrthotropicDamageDruckerPrager3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageTresca3D", mSmallStrainOrthotropicDamageTresca3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageMohrCoulomb3D", mSmallStrainOrthotropicDamageMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageModifiedMohrCoulomb3D", mSmallStrainOrthotropicDamageModifiedMohrCoulomb3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageSimoJu3D", mSmallStrainOrthotropicDamageSimoJu3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageRankine2D", mSmallStrainOrthotropicDamageRankine2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageVonMises2D", mSmallStrainOrthotropicDamageVonMises2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageDruckerPrager2D", mSmallStrainOrthotropicDamageDruckerPrager2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageTresca2D", mSmallStrainOrthotropicDamageTresca2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageMohrCoulomb2D", mSmallStrainOrthotropicDamageMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageModifiedMohrCoulomb2D", mSmallStrainOrthotropicDamageModifiedMohrCoulomb2D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainOrthotropicDamageSimoJu2D", mSmallStrainOrthotropicDamageSimoJu2D);

}
}  // namespace Kratos.
