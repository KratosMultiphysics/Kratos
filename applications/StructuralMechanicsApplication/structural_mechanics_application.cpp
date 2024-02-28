// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
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
#include "geometries/pyramid_3d_5.h"
#include "geometries/pyramid_3d_13.h"
#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"

#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/line_2d_4.h"
#include "geometries/line_2d_5.h"
#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/triangle_2d_10.h"
#include "geometries/triangle_2d_15.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"

namespace Kratos {

    // We define the node type
    typedef Node NodeType;

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
      mShellThickElement3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mShellThickCorotationalElement3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mShellThinCorotationalElement3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mShellThinElement3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mShellThinCorotationalElement3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mShellThickCorotationalElement3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      // Adding the membrane element
      mMembraneElement3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mMembraneElement3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      // Adding the SPRISM element
      mSolidShellElementSprism3D6N(0, Element::GeometryType::Pointer(new Prism3D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      // Adding the nodal concentrated element
      mNodalConcentratedElement2D1N(0, Element::GeometryType::Pointer(new Point2D<NodeType >(Element::GeometryType::PointsArrayType(1))), true),
      mNodalConcentratedDampedElement2D1N(0, Element::GeometryType::Pointer(new Point2D<NodeType >(Element::GeometryType::PointsArrayType(1))), false),
      mNodalConcentratedElement3D1N(0, Element::GeometryType::Pointer(new Point3D<NodeType >(Element::GeometryType::PointsArrayType(1))), true),
      mNodalConcentratedDampedElement3D1N(0, Element::GeometryType::Pointer(new Point3D<NodeType >(Element::GeometryType::PointsArrayType(1))), false),
     // Adding the mass elements
      mLineMassElement3D2N(0, Element::GeometryType::Pointer(new Line3D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
      mSurfaceMassElement3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mSurfaceMassElement3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      // Adding the kinematic linear elements
      mSmallDisplacement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mSmallDisplacement2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mSmallDisplacement2D6N(0, Element::GeometryType::Pointer(new Triangle2D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mSmallDisplacement2D8N(0, Element::GeometryType::Pointer(new Quadrilateral2D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mSmallDisplacement2D9N(0, Element::GeometryType::Pointer(new Quadrilateral2D9<NodeType >(Element::GeometryType::PointsArrayType(9)))),
      mSmallDisplacement2D10N(0, Element::GeometryType::Pointer(new Triangle2D10<NodeType >(Element::GeometryType::PointsArrayType(10)))),
      mSmallDisplacement2D15N(0, Element::GeometryType::Pointer(new Triangle2D15<NodeType >(Element::GeometryType::PointsArrayType(15)))),
      mSmallDisplacement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mSmallDisplacement3D5N(0, Element::GeometryType::Pointer(new Pyramid3D5<NodeType >(Element::GeometryType::PointsArrayType(5)))),
      mSmallDisplacement3D6N(0, Element::GeometryType::Pointer(new Prism3D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mSmallDisplacement3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mSmallDisplacement3D10N(0, Element::GeometryType::Pointer(new Tetrahedra3D10<NodeType >(Element::GeometryType::PointsArrayType(10)))),
      mSmallDisplacement3D13N(0, Element::GeometryType::Pointer(new Pyramid3D13<NodeType >(Element::GeometryType::PointsArrayType(13)))),
      mSmallDisplacement3D15N(0, Element::GeometryType::Pointer(new Prism3D15<NodeType >(Element::GeometryType::PointsArrayType(15)))),
      mSmallDisplacement3D20N(0, Element::GeometryType::Pointer(new Hexahedra3D20<NodeType >(Element::GeometryType::PointsArrayType(20)))),
      mSmallDisplacement3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27<NodeType >(Element::GeometryType::PointsArrayType(27)))),

      mSmallDisplacementBbar2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType>(Element::GeometryType::PointsArrayType(4)))),
      mSmallDisplacementBbar3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType>(Element::GeometryType::PointsArrayType(8)))),

      mSmallDisplacementShiftedBoundaryElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mSmallDisplacementShiftedBoundaryElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType>(Element::GeometryType::PointsArrayType(4)))),

      mSmallDisplacementMixedVolumetricStrainElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mSmallDisplacementMixedVolumetricStrainElement2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mSmallDisplacementMixedVolumetricStrainElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mSmallDisplacementMixedVolumetricStrainElement3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),

      mTotalLagrangianMixedVolumetricStrainElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mTotalLagrangianMixedVolumetricStrainElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),

      mTotalLagrangianQ1P0MixedElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mTotalLagrangianQ1P0MixedElement2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mTotalLagrangianQ1P0MixedElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mTotalLagrangianQ1P0MixedElement3D10N(0, Element::GeometryType::Pointer(new Tetrahedra3D10<NodeType >(Element::GeometryType::PointsArrayType(10)))),
      mTotalLagrangianQ1P0MixedElement3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mTotalLagrangianQ1P0MixedElement3D20N(0, Element::GeometryType::Pointer(new Hexahedra3D20<NodeType >(Element::GeometryType::PointsArrayType(20)))),

      mAxisymSmallDisplacement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mAxisymSmallDisplacement2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mAxisymSmallDisplacement2D6N(0, Element::GeometryType::Pointer(new Triangle2D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mAxisymSmallDisplacement2D8N(0, Element::GeometryType::Pointer(new Quadrilateral2D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mAxisymSmallDisplacement2D9N(0, Element::GeometryType::Pointer(new Quadrilateral2D9<NodeType >(Element::GeometryType::PointsArrayType(9)))),

      mZStrainDriven2p5DSmallDisplacement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mZStrainDriven2p5DSmallDisplacement2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mZStrainDriven2p5DSmallDisplacement2D6N(0, Element::GeometryType::Pointer(new Triangle2D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mZStrainDriven2p5DSmallDisplacement2D8N(0, Element::GeometryType::Pointer(new Quadrilateral2D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mZStrainDriven2p5DSmallDisplacement2D9N(0, Element::GeometryType::Pointer(new Quadrilateral2D9<NodeType >(Element::GeometryType::PointsArrayType(9)))),

      // Adding the Total lagrangian elements
      mTotalLagrangian2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mTotalLagrangian2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mTotalLagrangian2D6N(0, Element::GeometryType::Pointer(new Triangle2D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mTotalLagrangian2D8N(0, Element::GeometryType::Pointer(new Quadrilateral2D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mTotalLagrangian2D9N(0, Element::GeometryType::Pointer(new Quadrilateral2D9<NodeType >(Element::GeometryType::PointsArrayType(9)))),
      mTotalLagrangian2D10N(0, Element::GeometryType::Pointer(new Triangle2D10<NodeType >(Element::GeometryType::PointsArrayType(10)))),
      mTotalLagrangian2D15N(0, Element::GeometryType::Pointer(new Triangle2D15<NodeType >(Element::GeometryType::PointsArrayType(15)))),
      mTotalLagrangian3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mTotalLagrangian3D5N(0, Element::GeometryType::Pointer(new Pyramid3D5<NodeType >(Element::GeometryType::PointsArrayType(5)))),
      mTotalLagrangian3D6N(0, Element::GeometryType::Pointer(new Prism3D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mTotalLagrangian3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mTotalLagrangian3D10N(0, Element::GeometryType::Pointer(new Tetrahedra3D10<NodeType >(Element::GeometryType::PointsArrayType(10)))),
      mTotalLagrangian3D13N(0, Element::GeometryType::Pointer(new Pyramid3D13<NodeType >(Element::GeometryType::PointsArrayType(13)))),
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
      mUpdatedLagrangian2D10N(0, Element::GeometryType::Pointer(new Triangle2D10<NodeType >(Element::GeometryType::PointsArrayType(10)))),
      mUpdatedLagrangian2D15N(0, Element::GeometryType::Pointer(new Triangle2D15<NodeType >(Element::GeometryType::PointsArrayType(15)))),
      mUpdatedLagrangian3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mUpdatedLagrangian3D5N(0, Element::GeometryType::Pointer(new Pyramid3D5<NodeType >(Element::GeometryType::PointsArrayType(5)))),
      mUpdatedLagrangian3D6N(0, Element::GeometryType::Pointer(new Prism3D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mUpdatedLagrangian3D8N( 0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mUpdatedLagrangian3D10N(0, Element::GeometryType::Pointer(new Tetrahedra3D10<NodeType >(Element::GeometryType::PointsArrayType(10)))),
      mUpdatedLagrangian3D13N(0, Element::GeometryType::Pointer(new Pyramid3D13<NodeType >(Element::GeometryType::PointsArrayType(13)))),
      mUpdatedLagrangian3D15N(0, Element::GeometryType::Pointer(new Prism3D15<NodeType >(Element::GeometryType::PointsArrayType(15)))),
      mUpdatedLagrangian3D20N(0, Element::GeometryType::Pointer(new Hexahedra3D20<NodeType >(Element::GeometryType::PointsArrayType(20)))),
      mUpdatedLagrangian3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27<NodeType >(Element::GeometryType::PointsArrayType(27)))),
      mAxisymUpdatedLagrangian2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
      mAxisymUpdatedLagrangian2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
      mAxisymUpdatedLagrangian2D6N(0, Element::GeometryType::Pointer(new Triangle2D6<NodeType >(Element::GeometryType::PointsArrayType(6)))),
      mAxisymUpdatedLagrangian2D8N(0, Element::GeometryType::Pointer(new Quadrilateral2D8<NodeType >(Element::GeometryType::PointsArrayType(8)))),
      mAxisymUpdatedLagrangian2D9N(0, Element::GeometryType::Pointer(new Quadrilateral2D9<NodeType >(Element::GeometryType::PointsArrayType(9)))),
      // Adding the spring damper element
      mSpringDamperElement2D(0, Element::GeometryType::Pointer(new Line2D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
      mSpringDamperElement3D(0, Element::GeometryType::Pointer(new Line3D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
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
      mAdjointFiniteDifferenceSpringDamperElement3D2N(0, Element::GeometryType::Pointer(new Line3D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),

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
      mLineLoadCondition2D4N(0, Condition::GeometryType::Pointer(new Line2D4<NodeType >(Condition::GeometryType::PointsArrayType(4)))),
      mLineLoadCondition2D5N(0, Condition::GeometryType::Pointer(new Line2D5<NodeType >(Condition::GeometryType::PointsArrayType(5)))),
      mLineLoadCondition3D2N(0, Condition::GeometryType::Pointer(new Line3D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mLineLoadCondition3D3N(0, Condition::GeometryType::Pointer(new Line3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mSmallDisplacementLineLoadCondition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mSmallDisplacementLineLoadCondition2D3N(0, Condition::GeometryType::Pointer(new Line2D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mSmallDisplacementLineLoadCondition2D4N(0, Condition::GeometryType::Pointer(new Line2D4<NodeType >(Condition::GeometryType::PointsArrayType(4)))),
      mSmallDisplacementLineLoadCondition2D5N(0, Condition::GeometryType::Pointer(new Line2D5<NodeType >(Condition::GeometryType::PointsArrayType(5)))),
      mSmallDisplacementLineLoadCondition3D2N(0, Condition::GeometryType::Pointer(new Line3D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mSmallDisplacementLineLoadCondition3D3N(0, Condition::GeometryType::Pointer(new Line3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mAxisymLineLoadCondition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mAxisymLineLoadCondition2D3N(0, Condition::GeometryType::Pointer(new Line2D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      // Adding surface load conditions
      mSurfaceLoadCondition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mSurfaceLoadCondition3D4N(0, Condition::GeometryType::Pointer(new Quadrilateral3D4<NodeType >( Condition::GeometryType::PointsArrayType(4)))),
      mSurfaceLoadCondition3D6N(0, Condition::GeometryType::Pointer(new Triangle3D6<NodeType >(Condition::GeometryType::PointsArrayType(6)))),
      mSurfaceLoadCondition3D8N(0, Condition::GeometryType::Pointer(new Quadrilateral3D8<NodeType >(Condition::GeometryType::PointsArrayType(8)))),
      mSurfaceLoadCondition3D9N(0, Condition::GeometryType::Pointer(new Quadrilateral3D9<NodeType >(Condition::GeometryType::PointsArrayType(9)))),
      mSmallDisplacementSurfaceLoadCondition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mSmallDisplacementSurfaceLoadCondition3D4N(0, Condition::GeometryType::Pointer(new Quadrilateral3D4<NodeType >( Condition::GeometryType::PointsArrayType(4)))),
      mSmallDisplacementSurfaceLoadCondition3D6N(0, Condition::GeometryType::Pointer(new Triangle3D6<NodeType >(Condition::GeometryType::PointsArrayType(6)))),
      mSmallDisplacementSurfaceLoadCondition3D8N(0, Condition::GeometryType::Pointer(new Quadrilateral3D8<NodeType >(Condition::GeometryType::PointsArrayType(8)))),
      mSmallDisplacementSurfaceLoadCondition3D9N(0, Condition::GeometryType::Pointer(new Quadrilateral3D9<NodeType >(Condition::GeometryType::PointsArrayType(9)))),
      // Adding point moment conditions
      mPointMomentCondition3D1N(0, Condition::GeometryType::Pointer(new Point3D<NodeType >(Condition::GeometryType::PointsArrayType(1)))),

      // Adding adjoint conditions
      mAdjointSemiAnalyticPointLoadCondition2D1N(0, Condition::GeometryType::Pointer(new Point2D<NodeType >(Condition::GeometryType::PointsArrayType(1)))),
      mAdjointSemiAnalyticPointLoadCondition3D1N(0, Condition::GeometryType::Pointer(new Point3D<NodeType >(Condition::GeometryType::PointsArrayType(1)))),
      mAdjointSemiAnalyticSurfaceLoadCondition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mAdjointSemiAnalyticSurfaceLoadCondition3D4N(0, Condition::GeometryType::Pointer(new Quadrilateral3D4<NodeType >(Condition::GeometryType::PointsArrayType(4)))),
      mAdjointSemiAnalyticSmallDisplacementSurfaceLoadCondition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mAdjointSemiAnalyticSmallDisplacementSurfaceLoadCondition3D4N(0, Condition::GeometryType::Pointer(new Quadrilateral3D4<NodeType >(Condition::GeometryType::PointsArrayType(4)))),
      mAdjointSemiAnalyticLineLoadCondition3D2N(0, Condition::GeometryType::Pointer(new Line3D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mAdjointSemiAnalyticSmallDisplacementLineLoadCondition3D2N(0, Condition::GeometryType::Pointer(new Line3D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),

      // Adding the displacement-control condition
      mDisplacementControlCondition3D1N(0, Condition::GeometryType::Pointer(new Point3D<NodeType >(Condition::GeometryType::PointsArrayType(1)))),

      // Adding the SBM displacement conditions
      mDisplacementShiftedBoundaryCondition(0, Element::GeometryType::Pointer(new Geometry<Node>())),

      // Adding moving load conditions
      mMovingLoadCondition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mMovingLoadCondition2D3N(0, Condition::GeometryType::Pointer(new Line2D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))),
      mMovingLoadCondition3D2N(0, Condition::GeometryType::Pointer(new Line3D2<NodeType >(Condition::GeometryType::PointsArrayType(2)))),
      mMovingLoadCondition3D3N(0, Condition::GeometryType::Pointer(new Line3D3<NodeType >(Condition::GeometryType::PointsArrayType(3)))){}


void KratosStructuralMechanicsApplication::Register() {
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
    KRATOS_REGISTER_VARIABLE(LOCAL_PRESTRESS_AXIS_1);
    KRATOS_REGISTER_VARIABLE(LOCAL_PRESTRESS_AXIS_2);
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
    KRATOS_REGISTER_VARIABLE(BEAM_INITIAL_STRAIN_VECTOR)

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

    // Energy
    KRATOS_REGISTER_VARIABLE(ENERGY_DAMPING_DISSIPATION)

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
    KRATOS_REGISTER_VARIABLE( PRINCIPAL_CAUCHY_STRESS_VECTOR )
    KRATOS_REGISTER_VARIABLE( PRINCIPAL_PK2_STRESS_VECTOR )

    // Mixed formulations generalized variables
    KRATOS_REGISTER_VARIABLE( REACTION_STRAIN )

    // Formfinding
    KRATOS_REGISTER_VARIABLE(LAMBDA_MAX)
    KRATOS_REGISTER_VARIABLE(IS_FORMFINDING)
    KRATOS_REGISTER_VARIABLE(BASE_REF_1)
    KRATOS_REGISTER_VARIABLE(BASE_REF_2)

    // Cross section
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
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MOVING_LOAD)
    KRATOS_REGISTER_VARIABLE(MOVING_LOAD_LOCAL_DISTANCE)

    // Condition load variables
    KRATOS_REGISTER_VARIABLE(POINT_LOADS_VECTOR)
    KRATOS_REGISTER_VARIABLE(LINE_LOADS_VECTOR)
    KRATOS_REGISTER_VARIABLE(SURFACE_LOADS_VECTOR)
    KRATOS_REGISTER_VARIABLE(MOVING_LOADS_VECTOR)
    KRATOS_REGISTER_VARIABLE(POSITIVE_FACE_PRESSURES_VECTOR)
    KRATOS_REGISTER_VARIABLE(NEGATIVE_FACE_PRESSURES_VECTOR)

    // Displacement-Control variables
    KRATOS_REGISTER_VARIABLE(LOAD_FACTOR)
    KRATOS_REGISTER_VARIABLE(PRESCRIBED_DISPLACEMENT)

    // Response function variables
    KRATOS_REGISTER_VARIABLE(RESPONSE_VALUE)

    // Adjoint variables
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ADJOINT_DISPLACEMENT)
    KRATOS_REGISTER_VARIABLE(ADJOINT_PARTICULAR_DISPLACEMENT)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ADJOINT_ROTATION)
    KRATOS_REGISTER_VARIABLE(PERTURBATION_SIZE)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ADJOINT_CURVATURE)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ADJOINT_STRAIN)
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
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NODAL_DISPLACEMENT_STIFFNESS_SENSITIVITY);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NODAL_ROTATIONAL_STIFFNESS_SENSITIVITY);

    // Variables to for computing parts of sensitivity analysis
    KRATOS_REGISTER_VARIABLE( TRACED_STRESS_TYPE );
    KRATOS_REGISTER_VARIABLE( STRESS_DISP_DERIV_ON_GP );
    KRATOS_REGISTER_VARIABLE( STRESS_DISP_DERIV_ON_NODE);
    KRATOS_REGISTER_VARIABLE( STRESS_DESIGN_DERIVATIVE_ON_GP );
    KRATOS_REGISTER_VARIABLE( STRESS_DESIGN_DERIVATIVE_ON_NODE);
    KRATOS_REGISTER_VARIABLE( STRESS_ON_GP  );
    KRATOS_REGISTER_VARIABLE( STRESS_ON_NODE  );
    KRATOS_REGISTER_VARIABLE( DESIGN_VARIABLE_NAME );

    // for DEM-FEM 2D
    KRATOS_REGISTER_VARIABLE(IMPOSED_Z_STRAIN_VALUE)
    KRATOS_REGISTER_VARIABLE(IMPOSED_Z_STRAIN_OPTION)

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
    KRATOS_REGISTER_ELEMENT("MembraneElement3D4N", mMembraneElement3D4N)
    KRATOS_REGISTER_ELEMENT("MembraneElement3D3N", mMembraneElement3D3N)

    // Register the SPRISM element
    KRATOS_REGISTER_ELEMENT("SolidShellElementSprism3D6N", mSolidShellElementSprism3D6N);

    // Register the nodal concentrated element
    KRATOS_REGISTER_ELEMENT("NodalConcentratedDampedElement3D1N", mNodalConcentratedDampedElement3D1N);
    KRATOS_REGISTER_ELEMENT("NodalConcentratedElement2D1N", mNodalConcentratedElement2D1N);
    KRATOS_REGISTER_ELEMENT("NodalConcentratedDampedElement2D1N", mNodalConcentratedDampedElement2D1N);
    KRATOS_REGISTER_ELEMENT("NodalConcentratedElement3D1N", mNodalConcentratedElement3D1N);

    // Register mass elements
    KRATOS_REGISTER_ELEMENT("LineMassElement3D2N",    mLineMassElement3D2N);
    KRATOS_REGISTER_ELEMENT("SurfaceMassElement3D3N", mSurfaceMassElement3D3N);
    KRATOS_REGISTER_ELEMENT("SurfaceMassElement3D4N", mSurfaceMassElement3D4N);

    // SOLID ELEMENTS
    // Small displacement elements
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement2D3N", mSmallDisplacement2D3N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement2D4N", mSmallDisplacement2D4N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement2D6N", mSmallDisplacement2D6N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement2D8N", mSmallDisplacement2D8N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement2D9N", mSmallDisplacement2D9N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement2D10N", mSmallDisplacement2D10N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement2D15N", mSmallDisplacement2D15N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D4N", mSmallDisplacement3D4N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D5N", mSmallDisplacement3D5N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D6N", mSmallDisplacement3D6N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D8N", mSmallDisplacement3D8N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D10N", mSmallDisplacement3D10N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D13N", mSmallDisplacement3D13N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D15N", mSmallDisplacement3D15N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D20N", mSmallDisplacement3D20N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementElement3D27N", mSmallDisplacement3D27N)

    KRATOS_REGISTER_ELEMENT("SmallDisplacementShiftedBoundaryElement2D3N", mSmallDisplacementShiftedBoundaryElement2D3N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementShiftedBoundaryElement3D4N", mSmallDisplacementShiftedBoundaryElement3D4N)

    KRATOS_REGISTER_ELEMENT("SmallDisplacementMixedVolumetricStrainElement2D3N", mSmallDisplacementMixedVolumetricStrainElement2D3N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementMixedVolumetricStrainElement2D4N", mSmallDisplacementMixedVolumetricStrainElement2D4N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementMixedVolumetricStrainElement3D4N", mSmallDisplacementMixedVolumetricStrainElement3D4N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementMixedVolumetricStrainElement3D8N", mSmallDisplacementMixedVolumetricStrainElement3D8N)

    KRATOS_REGISTER_ELEMENT("TotalLagrangianMixedVolumetricStrainElement2D3N", mTotalLagrangianMixedVolumetricStrainElement2D3N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianMixedVolumetricStrainElement3D4N", mTotalLagrangianMixedVolumetricStrainElement3D4N)

    KRATOS_REGISTER_ELEMENT("TotalLagrangianQ1P0MixedElement2D3N", mTotalLagrangianQ1P0MixedElement2D3N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianQ1P0MixedElement2D4N", mTotalLagrangianQ1P0MixedElement2D4N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianQ1P0MixedElement3D4N", mTotalLagrangianQ1P0MixedElement3D4N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianQ1P0MixedElement3D10N", mTotalLagrangianQ1P0MixedElement3D10N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianQ1P0MixedElement3D8N", mTotalLagrangianQ1P0MixedElement3D8N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianQ1P0MixedElement3D20N", mTotalLagrangianQ1P0MixedElement3D20N)

    KRATOS_REGISTER_ELEMENT("SmallDisplacementBbarElement2D4N", mSmallDisplacementBbar2D4N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementBbarElement3D8N", mSmallDisplacementBbar3D8N)

    KRATOS_REGISTER_ELEMENT("AxisymSmallDisplacementElement2D3N", mAxisymSmallDisplacement2D3N)
    KRATOS_REGISTER_ELEMENT("AxisymSmallDisplacementElement2D4N", mAxisymSmallDisplacement2D4N)
    KRATOS_REGISTER_ELEMENT("AxisymSmallDisplacementElement2D6N", mAxisymSmallDisplacement2D6N)
    KRATOS_REGISTER_ELEMENT("AxisymSmallDisplacementElement2D8N", mAxisymSmallDisplacement2D8N)
    KRATOS_REGISTER_ELEMENT("AxisymSmallDisplacementElement2D9N", mAxisymSmallDisplacement2D9N)

    KRATOS_REGISTER_ELEMENT("ZStrainDriven2p5DSmallDisplacementElement2D3N", mZStrainDriven2p5DSmallDisplacement2D3N)
    KRATOS_REGISTER_ELEMENT("ZStrainDriven2p5DSmallDisplacementElement2D4N", mZStrainDriven2p5DSmallDisplacement2D4N)
    KRATOS_REGISTER_ELEMENT("ZStrainDriven2p5DSmallDisplacementElement2D6N", mZStrainDriven2p5DSmallDisplacement2D6N)
    KRATOS_REGISTER_ELEMENT("ZStrainDriven2p5DSmallDisplacementElement2D8N", mZStrainDriven2p5DSmallDisplacement2D8N)
    KRATOS_REGISTER_ELEMENT("ZStrainDriven2p5DSmallDisplacementElement2D9N", mZStrainDriven2p5DSmallDisplacement2D9N)

    // Total lagrangian elements
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement2D3N", mTotalLagrangian2D3N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement2D4N", mTotalLagrangian2D4N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement2D6N", mTotalLagrangian2D6N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement2D8N", mTotalLagrangian2D8N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement2D9N", mTotalLagrangian2D9N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement2D10N", mTotalLagrangian2D10N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement2D15N", mTotalLagrangian2D15N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement3D4N", mTotalLagrangian3D4N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement3D5N", mTotalLagrangian3D5N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement3D6N", mTotalLagrangian3D6N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement3D8N", mTotalLagrangian3D8N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement3D10N", mTotalLagrangian3D10N)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianElement3D13N", mTotalLagrangian3D13N)
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
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement2D10N", mUpdatedLagrangian2D10N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement2D15N", mUpdatedLagrangian2D15N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement3D4N", mUpdatedLagrangian3D4N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement3D5N", mUpdatedLagrangian3D5N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement3D6N", mUpdatedLagrangian3D6N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement3D8N", mUpdatedLagrangian3D8N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement3D10N", mUpdatedLagrangian3D10N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement3D13N", mUpdatedLagrangian3D13N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement3D15N", mUpdatedLagrangian3D15N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement3D20N", mUpdatedLagrangian3D20N)
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianElement3D27N", mUpdatedLagrangian3D27N)

    KRATOS_REGISTER_ELEMENT("AxisymUpdatedLagrangianElement2D3N", mAxisymUpdatedLagrangian2D3N)
    KRATOS_REGISTER_ELEMENT("AxisymUpdatedLagrangianElement2D4N", mAxisymUpdatedLagrangian2D4N)
    KRATOS_REGISTER_ELEMENT("AxisymUpdatedLagrangianElement2D6N", mAxisymUpdatedLagrangian2D6N)
    KRATOS_REGISTER_ELEMENT("AxisymUpdatedLagrangianElement2D8N", mAxisymUpdatedLagrangian2D8N)
    KRATOS_REGISTER_ELEMENT("AxisymUpdatedLagrangianElement2D9N", mAxisymUpdatedLagrangian2D9N)

    // Register the spring damper element
    KRATOS_REGISTER_ELEMENT("SpringDamperElement2D", mSpringDamperElement2D);
    KRATOS_REGISTER_ELEMENT("SpringDamperElement3D", mSpringDamperElement3D);
    KRATOS_REGISTER_ELEMENT("SpringDamperElement3D2N", mSpringDamperElement3D); // DEPRECATED, only here for backward compatibility, should be removed after some time

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
    KRATOS_REGISTER_ELEMENT("AdjointFiniteDifferenceSpringDamperElement3D2N", mAdjointFiniteDifferenceSpringDamperElement3D2N)

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
    KRATOS_REGISTER_CONDITION("LineLoadCondition2D4N", mLineLoadCondition2D4N)
    KRATOS_REGISTER_CONDITION("LineLoadCondition2D5N", mLineLoadCondition2D5N)
    KRATOS_REGISTER_CONDITION("LineLoadCondition3D2N", mLineLoadCondition3D2N)
    KRATOS_REGISTER_CONDITION("LineLoadCondition3D3N", mLineLoadCondition3D3N)

    KRATOS_REGISTER_CONDITION("SmallDisplacementLineLoadCondition2D2N", mSmallDisplacementLineLoadCondition2D2N)
    KRATOS_REGISTER_CONDITION("SmallDisplacementLineLoadCondition2D3N", mSmallDisplacementLineLoadCondition2D3N)
    KRATOS_REGISTER_CONDITION("SmallDisplacementLineLoadCondition2D4N", mSmallDisplacementLineLoadCondition2D4N)
    KRATOS_REGISTER_CONDITION("SmallDisplacementLineLoadCondition2D5N", mSmallDisplacementLineLoadCondition2D5N)
    KRATOS_REGISTER_CONDITION("SmallDisplacementLineLoadCondition3D2N", mSmallDisplacementLineLoadCondition3D2N)
    KRATOS_REGISTER_CONDITION("SmallDisplacementLineLoadCondition3D3N", mSmallDisplacementLineLoadCondition3D3N)

    KRATOS_REGISTER_CONDITION("AxisymLineLoadCondition2D2N", mAxisymLineLoadCondition2D2N)
    KRATOS_REGISTER_CONDITION("AxisymLineLoadCondition2D3N", mAxisymLineLoadCondition2D3N)

    // Surface loads
    KRATOS_REGISTER_CONDITION("SurfaceLoadCondition3D3N", mSurfaceLoadCondition3D3N)
    KRATOS_REGISTER_CONDITION("SurfaceLoadCondition3D4N", mSurfaceLoadCondition3D4N)
    KRATOS_REGISTER_CONDITION("SurfaceLoadCondition3D6N", mSurfaceLoadCondition3D6N)
    KRATOS_REGISTER_CONDITION("SurfaceLoadCondition3D8N", mSurfaceLoadCondition3D8N)
    KRATOS_REGISTER_CONDITION("SurfaceLoadCondition3D9N", mSurfaceLoadCondition3D9N)

    KRATOS_REGISTER_CONDITION("SmallDisplacementSurfaceLoadCondition3D3N", mSmallDisplacementSurfaceLoadCondition3D3N)
    KRATOS_REGISTER_CONDITION("SmallDisplacementSurfaceLoadCondition3D4N", mSmallDisplacementSurfaceLoadCondition3D4N)
    KRATOS_REGISTER_CONDITION("SmallDisplacementSurfaceLoadCondition3D6N", mSmallDisplacementSurfaceLoadCondition3D6N)
    KRATOS_REGISTER_CONDITION("SmallDisplacementSurfaceLoadCondition3D8N", mSmallDisplacementSurfaceLoadCondition3D8N)
    KRATOS_REGISTER_CONDITION("SmallDisplacementSurfaceLoadCondition3D9N", mSmallDisplacementSurfaceLoadCondition3D9N)

    // Point moment
    KRATOS_REGISTER_CONDITION("PointMomentCondition3D1N", mPointMomentCondition3D1N);

    // Adjoint conditions
    KRATOS_REGISTER_CONDITION("AdjointSemiAnalyticPointLoadCondition2D1N", mAdjointSemiAnalyticPointLoadCondition2D1N )
    KRATOS_REGISTER_CONDITION("AdjointSemiAnalyticPointLoadCondition3D1N", mAdjointSemiAnalyticPointLoadCondition3D1N )
    KRATOS_REGISTER_CONDITION("AdjointSemiAnalyticSurfaceLoadCondition3D3N", mAdjointSemiAnalyticSurfaceLoadCondition3D3N )
    KRATOS_REGISTER_CONDITION("AdjointSemiAnalyticSurfaceLoadCondition3D4N", mAdjointSemiAnalyticSurfaceLoadCondition3D4N )
    KRATOS_REGISTER_CONDITION("AdjointSemiAnalyticSmallDisplacementSurfaceLoadCondition3D3N", mAdjointSemiAnalyticSmallDisplacementSurfaceLoadCondition3D3N )
    KRATOS_REGISTER_CONDITION("AdjointSemiAnalyticSmallDisplacementSurfaceLoadCondition3D4N", mAdjointSemiAnalyticSmallDisplacementSurfaceLoadCondition3D4N )
    KRATOS_REGISTER_CONDITION("AdjointSemiAnalyticLineLoadCondition3D2N", mAdjointSemiAnalyticLineLoadCondition3D2N)
    KRATOS_REGISTER_CONDITION("AdjointSemiAnalyticSmallDisplacementLineLoadCondition3D2N", mAdjointSemiAnalyticSmallDisplacementLineLoadCondition3D2N)

    // Displacement-Control Conditions
    KRATOS_REGISTER_CONDITION("DisplacementControlCondition3D1N", mDisplacementControlCondition3D1N)

    // SBM displacement conditions
    KRATOS_REGISTER_CONDITION("DisplacementShiftedBoundaryCondition", mDisplacementShiftedBoundaryCondition)

    // Moving loads
    KRATOS_REGISTER_CONDITION("MovingLoadCondition2D2N", mMovingLoadCondition2D2N)
    KRATOS_REGISTER_CONDITION("MovingLoadCondition2D3N", mMovingLoadCondition2D3N)
    KRATOS_REGISTER_CONDITION("MovingLoadCondition3D2N", mMovingLoadCondition3D2N)
    KRATOS_REGISTER_CONDITION("MovingLoadCondition3D3N", mMovingLoadCondition3D3N)


    // Register linear elastics laws
    KRATOS_REGISTER_CONSTITUTIVE_LAW("TrussConstitutiveLaw", mTrussConstitutiveLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("BeamConstitutiveLaw", mBeamConstitutiveLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElastic3DLaw", mElasticIsotropic3D);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticPlaneStrain2DLaw", mLinearPlaneStrain);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticPlaneStress2DLaw", mLinearPlaneStress);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticAxisym2DLaw", mAxisymElasticIsotropic);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("UserProvidedLinearElastic2DLaw", mUserProvidedLinearElastic2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("UserProvidedLinearElastic3DLaw", mUserProvidedLinearElastic3DLaw);

}
}  // namespace Kratos.
