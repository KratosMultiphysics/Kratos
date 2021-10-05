// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//


// System includes


// External includes

// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"
#include "geometries/quadrilateral_interface_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/prism_interface_3d_6.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"
#include "geometries/hexahedra_interface_3d_8.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "geometries/quadrilateral_interface_3d_4.h"
#include "includes/variables.h"


// Application includes
#include "geo_mechanics_application.h"
#include "geo_mechanics_application_variables.h"


namespace Kratos {
// We define the node type
typedef Node<3> NodeType;

KratosGeoMechanicsApplication::KratosGeoMechanicsApplication():
    KratosApplication("GeoMechanicsApplication"),

    // transient one-phase flow elements:
    mTransientPwElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <NodeType >( Element::GeometryType::PointsArrayType(3)))),
    mTransientPwElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <NodeType >( Element::GeometryType::PointsArrayType(4)))),
    mTransientPwElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <NodeType >( Element::GeometryType::PointsArrayType(4)))),
    mTransientPwElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <NodeType >( Element::GeometryType::PointsArrayType(8)))),

    mTransientPwElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <NodeType >( Element::GeometryType::PointsArrayType(6)))),
    mTransientPwElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <NodeType >( Element::GeometryType::PointsArrayType(8)))),
    mTransientPwElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <NodeType >( Element::GeometryType::PointsArrayType(9)))),
    mTransientPwElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <NodeType >( Element::GeometryType::PointsArrayType(10)))),
    mTransientPwElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <NodeType >( Element::GeometryType::PointsArrayType(20)))),
    mTransientPwElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <NodeType >( Element::GeometryType::PointsArrayType(27)))),

    // Steady-State one-phase flow elements:
    mSteadyStatePwElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <NodeType >( Element::GeometryType::PointsArrayType(3)))),
    mSteadyStatePwElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <NodeType >( Element::GeometryType::PointsArrayType(4)))),
    mSteadyStatePwElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <NodeType >( Element::GeometryType::PointsArrayType(4)))),
    mSteadyStatePwElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <NodeType >( Element::GeometryType::PointsArrayType(8)))),

    mSteadyStatePwElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <NodeType >( Element::GeometryType::PointsArrayType(6)))),
    mSteadyStatePwElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <NodeType >( Element::GeometryType::PointsArrayType(8)))),
    mSteadyStatePwElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <NodeType >( Element::GeometryType::PointsArrayType(9)))),
    mSteadyStatePwElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <NodeType >( Element::GeometryType::PointsArrayType(10)))),
    mSteadyStatePwElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <NodeType >( Element::GeometryType::PointsArrayType(20)))),
    mSteadyStatePwElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <NodeType >( Element::GeometryType::PointsArrayType(27)))),


    // small strain elements
    mUPwSmallStrainElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <NodeType >( Element::GeometryType::PointsArrayType(3)))),
    mUPwSmallStrainElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <NodeType >( Element::GeometryType::PointsArrayType(4)))),
    mUPwSmallStrainElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <NodeType >( Element::GeometryType::PointsArrayType(4)))),
    mUPwSmallStrainElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <NodeType >( Element::GeometryType::PointsArrayType(8)))),

    mUPwSmallStrainElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <NodeType >( Element::GeometryType::PointsArrayType(6)))),
    mUPwSmallStrainElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <NodeType >( Element::GeometryType::PointsArrayType(8)))),
    mUPwSmallStrainElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <NodeType >( Element::GeometryType::PointsArrayType(9)))),
    mUPwSmallStrainElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <NodeType >( Element::GeometryType::PointsArrayType(10)))),
    mUPwSmallStrainElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <NodeType >( Element::GeometryType::PointsArrayType(20)))),
    mUPwSmallStrainElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <NodeType >( Element::GeometryType::PointsArrayType(27)))),

    // drained small strain elements
    mDrainedUPwSmallStrainElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <NodeType >( Element::GeometryType::PointsArrayType(3)))),
    mDrainedUPwSmallStrainElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <NodeType >( Element::GeometryType::PointsArrayType(4)))),
    mDrainedUPwSmallStrainElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <NodeType >( Element::GeometryType::PointsArrayType(4)))),
    mDrainedUPwSmallStrainElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <NodeType >( Element::GeometryType::PointsArrayType(8)))),

    // undrained small strain elements
    mUndrainedUPwSmallStrainElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <NodeType >( Element::GeometryType::PointsArrayType(3)))),
    mUndrainedUPwSmallStrainElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <NodeType >( Element::GeometryType::PointsArrayType(4)))),
    mUndrainedUPwSmallStrainElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <NodeType >( Element::GeometryType::PointsArrayType(4)))),
    mUndrainedUPwSmallStrainElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <NodeType >( Element::GeometryType::PointsArrayType(8)))),

    // small strain FIC elements
    mUPwSmallStrainFICElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <NodeType >( Element::GeometryType::PointsArrayType(3)))),
    mUPwSmallStrainFICElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <NodeType >( Element::GeometryType::PointsArrayType(4)))),
    mUPwSmallStrainFICElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <NodeType >( Element::GeometryType::PointsArrayType(4)))),
    mUPwSmallStrainFICElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <NodeType >( Element::GeometryType::PointsArrayType(8)))),

    // small strain different order elements
    mSmallStrainUPwDiffOrderElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <NodeType >( Element::GeometryType::PointsArrayType(6)))),
    mSmallStrainUPwDiffOrderElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <NodeType >( Element::GeometryType::PointsArrayType(8)))),
    mSmallStrainUPwDiffOrderElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <NodeType >( Element::GeometryType::PointsArrayType(9)))),
    mSmallStrainUPwDiffOrderElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <NodeType >( Element::GeometryType::PointsArrayType(10)))),
    mSmallStrainUPwDiffOrderElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <NodeType >( Element::GeometryType::PointsArrayType(20)))),
    mSmallStrainUPwDiffOrderElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <NodeType >( Element::GeometryType::PointsArrayType(27)))),

    // small strain axisymmtric elements:
    mUPwSmallStrainAxisymmetricElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <NodeType >( Element::GeometryType::PointsArrayType(3)))),
    mUPwSmallStrainAxisymmetricElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <NodeType >( Element::GeometryType::PointsArrayType(4)))),
    mUPwSmallStrainAxisymmetricElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <NodeType >( Element::GeometryType::PointsArrayType(6)))),
    mUPwSmallStrainAxisymmetricElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <NodeType >( Element::GeometryType::PointsArrayType(8)))),
    mUPwSmallStrainAxisymmetricElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <NodeType >( Element::GeometryType::PointsArrayType(9)))),
    mUPwSmallStrainAxisymmetricFICElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <NodeType >( Element::GeometryType::PointsArrayType(3)))),
    mUPwSmallStrainAxisymmetricFICElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <NodeType >( Element::GeometryType::PointsArrayType(4)))),
    mSmallStrainUPwDiffOrderAxisymmetricElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <NodeType >( Element::GeometryType::PointsArrayType(6)))),
    mSmallStrainUPwDiffOrderAxisymmetricElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <NodeType >( Element::GeometryType::PointsArrayType(8)))),
    mSmallStrainUPwDiffOrderAxisymmetricElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <NodeType >( Element::GeometryType::PointsArrayType(9)))),

    // small strain interface elements
    mUPwSmallStrainInterfaceElement2D4N( 0, Element::GeometryType::Pointer( new QuadrilateralInterface2D4 <NodeType >( Element::GeometryType::PointsArrayType(4)))),
    mUPwSmallStrainInterfaceElement3D6N( 0, Element::GeometryType::Pointer( new PrismInterface3D6 <NodeType >( Element::GeometryType::PointsArrayType(6)))),
    mUPwSmallStrainInterfaceElement3D8N( 0, Element::GeometryType::Pointer( new HexahedraInterface3D8 <NodeType >( Element::GeometryType::PointsArrayType(8)))),

    mUPwSmallStrainLinkInterfaceElement2D4N( 0, Element::GeometryType::Pointer( new QuadrilateralInterface2D4 <NodeType >( Element::GeometryType::PointsArrayType(4)))),
    mUPwSmallStrainLinkInterfaceElement3D6N( 0, Element::GeometryType::Pointer( new PrismInterface3D6 <NodeType >( Element::GeometryType::PointsArrayType(6)))),
    mUPwSmallStrainLinkInterfaceElement3D8N( 0, Element::GeometryType::Pointer( new HexahedraInterface3D8 <NodeType >( Element::GeometryType::PointsArrayType(8)))),

    // Updated-Lagrangian elements
    mUPwUpdatedLagrangianElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <NodeType >( Element::GeometryType::PointsArrayType(3)))),
    mUPwUpdatedLagrangianElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <NodeType >( Element::GeometryType::PointsArrayType(4)))),
    mUPwUpdatedLagrangianElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <NodeType >( Element::GeometryType::PointsArrayType(4)))),
    mUPwUpdatedLagrangianElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <NodeType >( Element::GeometryType::PointsArrayType(8)))),

    mUPwUpdatedLagrangianElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <NodeType >( Element::GeometryType::PointsArrayType(6)))),
    mUPwUpdatedLagrangianElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <NodeType >( Element::GeometryType::PointsArrayType(8)))),
    mUPwUpdatedLagrangianElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <NodeType >( Element::GeometryType::PointsArrayType(9)))),
    mUPwUpdatedLagrangianElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <NodeType >( Element::GeometryType::PointsArrayType(10)))),
    mUPwUpdatedLagrangianElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <NodeType >( Element::GeometryType::PointsArrayType(20)))),
    mUPwUpdatedLagrangianElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <NodeType >( Element::GeometryType::PointsArrayType(27)))),

    mUPwUpdatedLagrangianFICElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <NodeType >( Element::GeometryType::PointsArrayType(3)))),
    mUPwUpdatedLagrangianFICElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <NodeType >( Element::GeometryType::PointsArrayType(4)))),
    mUPwUpdatedLagrangianFICElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <NodeType >( Element::GeometryType::PointsArrayType(4)))),
    mUPwUpdatedLagrangianFICElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <NodeType >( Element::GeometryType::PointsArrayType(8)))),

    // Updated-Lagrangian different order elements
    mUpdatedLagrangianUPwDiffOrderElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <NodeType >( Element::GeometryType::PointsArrayType(6)))),
    mUpdatedLagrangianUPwDiffOrderElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <NodeType >( Element::GeometryType::PointsArrayType(8)))),
    mUpdatedLagrangianUPwDiffOrderElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <NodeType >( Element::GeometryType::PointsArrayType(9)))),
    mUpdatedLagrangianUPwDiffOrderElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <NodeType >( Element::GeometryType::PointsArrayType(10)))),
    mUpdatedLagrangianUPwDiffOrderElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <NodeType >( Element::GeometryType::PointsArrayType(20)))),
    mUpdatedLagrangianUPwDiffOrderElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <NodeType >( Element::GeometryType::PointsArrayType(27)))),

    // geo structural elements
    mGeoCrBeamElement2D2N(0, Element::GeometryType::Pointer(new Line2D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
    mGeoCrBeamElement3D2N(0, Element::GeometryType::Pointer(new Line3D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
    mGeoCrBeamElementLinear2D2N(0, Element::GeometryType::Pointer(new Line2D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
    mGeoCrBeamElementLinear3D2N(0, Element::GeometryType::Pointer(new Line3D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
    mGeoTrussElement2D2N(0, Element::GeometryType::Pointer(new Line2D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
    mGeoTrussElement3D2N(0, Element::GeometryType::Pointer(new Line3D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
    mGeoLinearTrussElement2D2N(0, Element::GeometryType::Pointer(new Line2D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
    mGeoLinearTrussElement3D2N(0, Element::GeometryType::Pointer(new Line3D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
    mGeoCableElement2D2N(0, Element::GeometryType::Pointer(new Line2D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),
    mGeoCableElement3D2N(0, Element::GeometryType::Pointer(new Line3D2<NodeType >(Element::GeometryType::PointsArrayType(2)))),

    // conditions
    mUPwForceCondition2D1N( 0, Condition::GeometryType::Pointer( new Point2D<NodeType >( Condition::GeometryType::PointsArrayType(1)))),
    mUPwForceCondition3D1N( 0, Condition::GeometryType::Pointer( new Point3D<NodeType >( Condition::GeometryType::PointsArrayType(1)))),
    mUPwFaceLoadCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<NodeType >( Condition::GeometryType::PointsArrayType(2)))),
    mUPwFaceLoadCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <NodeType >( Condition::GeometryType::PointsArrayType(3)))),
    mUPwFaceLoadCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <NodeType >( Condition::GeometryType::PointsArrayType(4)))),
    mUPwNormalFaceLoadCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<NodeType >( Condition::GeometryType::PointsArrayType(2)))),
    mUPwNormalFaceLoadCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <NodeType >( Condition::GeometryType::PointsArrayType(3)))),
    mUPwNormalFaceLoadCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <NodeType >( Condition::GeometryType::PointsArrayType(4)))),
    mUPwNormalFluxCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<NodeType >( Condition::GeometryType::PointsArrayType(2)))),
    mUPwNormalFluxCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <NodeType >( Condition::GeometryType::PointsArrayType(3)))),
    mUPwNormalFluxCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <NodeType >( Condition::GeometryType::PointsArrayType(4)))),

    mUPwFaceLoadCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3<NodeType >( Condition::GeometryType::PointsArrayType(3)))),

    mUPwFaceLoadInterfaceCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<NodeType >( Condition::GeometryType::PointsArrayType(2)))),
    mUPwFaceLoadInterfaceCondition3D4N( 0, Condition::GeometryType::Pointer( new QuadrilateralInterface3D4 <NodeType >( Condition::GeometryType::PointsArrayType(4)))),
    mUPwNormalFluxInterfaceCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<NodeType >( Condition::GeometryType::PointsArrayType(2)))),
    mUPwNormalFluxInterfaceCondition3D4N( 0, Condition::GeometryType::Pointer( new QuadrilateralInterface3D4 <NodeType >( Condition::GeometryType::PointsArrayType(4)))),

    mUPwNormalFluxFICCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<NodeType >( Condition::GeometryType::PointsArrayType(2)))),
    mUPwNormalFluxFICCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <NodeType >( Condition::GeometryType::PointsArrayType(3)))),
    mUPwNormalFluxFICCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <NodeType >( Condition::GeometryType::PointsArrayType(4)))),

    mLineLoadDiffOrderCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3<NodeType >( Condition::GeometryType::PointsArrayType(3)))),
    mLineNormalLoadDiffOrderCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3<NodeType >( Condition::GeometryType::PointsArrayType(3)))),
    mLineNormalFluidFluxDiffOrderCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3<NodeType >( Condition::GeometryType::PointsArrayType(3)))),
    mSurfaceLoadDiffOrderCondition3D6N( 0, Condition::GeometryType::Pointer( new Triangle3D6 <NodeType >( Condition::GeometryType::PointsArrayType(6)))),
    mSurfaceLoadDiffOrderCondition3D8N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D8 <NodeType >( Condition::GeometryType::PointsArrayType(8)))),
    mSurfaceLoadDiffOrderCondition3D9N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D9 <NodeType >( Condition::GeometryType::PointsArrayType(9)))),
    mSurfaceNormalLoadDiffOrderCondition3D6N( 0, Condition::GeometryType::Pointer( new Triangle3D6 <NodeType >( Condition::GeometryType::PointsArrayType(6)))),
    mSurfaceNormalLoadDiffOrderCondition3D8N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D8 <NodeType >( Condition::GeometryType::PointsArrayType(8)))),
    mSurfaceNormalLoadDiffOrderCondition3D9N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D9 <NodeType >( Condition::GeometryType::PointsArrayType(9)))),
    mSurfaceNormalFluidFluxDiffOrderCondition3D6N( 0, Condition::GeometryType::Pointer( new Triangle3D6 <NodeType >( Condition::GeometryType::PointsArrayType(6)))),
    mSurfaceNormalFluidFluxDiffOrderCondition3D8N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D8 <NodeType >( Condition::GeometryType::PointsArrayType(8)))),
    mSurfaceNormalFluidFluxDiffOrderCondition3D9N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D9 <NodeType >( Condition::GeometryType::PointsArrayType(9))))

    {}

void KratosGeoMechanicsApplication::Register() {
    // calling base class register to register Kratos components
    KratosApplication::Register();
    KRATOS_INFO("") << " KRATOS___                             \n"
                    << "     //   ) )                          \n"
                    << "    //         ___      ___            \n"
                    << "   //  ____  //___) ) //   ) )         \n"
                    << "  //    / / //       //   / /          \n"
                    << " ((____/ / ((____   ((___/ /  MECHANICS\n"
                    << " Initializing KratosGeoMechanicsApplication..." << std::endl;



    //Register Elements
    // transient one-phase flow elements:
    KRATOS_REGISTER_ELEMENT( "TransientPwElement2D3N", mTransientPwElement2D3N )
    KRATOS_REGISTER_ELEMENT( "TransientPwElement2D4N", mTransientPwElement2D4N )
    KRATOS_REGISTER_ELEMENT( "TransientPwElement3D4N", mTransientPwElement3D4N )
    KRATOS_REGISTER_ELEMENT( "TransientPwElement3D8N", mTransientPwElement3D8N )

    KRATOS_REGISTER_ELEMENT( "TransientPwElement2D6N", mTransientPwElement2D6N )
    KRATOS_REGISTER_ELEMENT( "TransientPwElement2D8N", mTransientPwElement2D8N )
    KRATOS_REGISTER_ELEMENT( "TransientPwElement2D9N", mTransientPwElement2D9N )
    KRATOS_REGISTER_ELEMENT( "TransientPwElement3D10N", mTransientPwElement3D10N )
    KRATOS_REGISTER_ELEMENT( "TransientPwElement3D20N", mTransientPwElement3D20N )
    KRATOS_REGISTER_ELEMENT( "TransientPwElement3D27N", mTransientPwElement3D27N )

    // Steady-State one-phase flow elements:
    KRATOS_REGISTER_ELEMENT( "SteadyStatePwElement2D3N", mSteadyStatePwElement2D3N )
    KRATOS_REGISTER_ELEMENT( "SteadyStatePwElement2D4N", mSteadyStatePwElement2D4N )
    KRATOS_REGISTER_ELEMENT( "SteadyStatePwElement3D4N", mSteadyStatePwElement3D4N )
    KRATOS_REGISTER_ELEMENT( "SteadyStatePwElement3D8N", mSteadyStatePwElement3D8N )

    KRATOS_REGISTER_ELEMENT( "SteadyStatePwElement2D6N", mSteadyStatePwElement2D6N )
    KRATOS_REGISTER_ELEMENT( "SteadyStatePwElement2D8N", mSteadyStatePwElement2D8N )
    KRATOS_REGISTER_ELEMENT( "SteadyStatePwElement2D9N", mSteadyStatePwElement2D9N )
    KRATOS_REGISTER_ELEMENT( "SteadyStatePwElement3D10N", mSteadyStatePwElement3D10N )
    KRATOS_REGISTER_ELEMENT( "SteadyStatePwElement3D20N", mSteadyStatePwElement3D20N )
    KRATOS_REGISTER_ELEMENT( "SteadyStatePwElement3D27N", mSteadyStatePwElement3D27N )

    // Small strain elements
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainElement2D3N", mUPwSmallStrainElement2D3N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainElement2D4N", mUPwSmallStrainElement2D4N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainElement3D4N", mUPwSmallStrainElement3D4N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainElement3D8N", mUPwSmallStrainElement3D8N )

    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainElement2D6N", mUPwSmallStrainElement2D6N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainElement2D8N", mUPwSmallStrainElement2D8N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainElement2D9N", mUPwSmallStrainElement2D9N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainElement3D10N", mUPwSmallStrainElement3D10N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainElement3D20N", mUPwSmallStrainElement3D20N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainElement3D27N", mUPwSmallStrainElement3D27N )

    // Drained small strain elements
    KRATOS_REGISTER_ELEMENT( "DrainedUPwSmallStrainElement2D3N", mDrainedUPwSmallStrainElement2D3N )
    KRATOS_REGISTER_ELEMENT( "DrainedUPwSmallStrainElement2D4N", mDrainedUPwSmallStrainElement2D4N )
    KRATOS_REGISTER_ELEMENT( "DrainedUPwSmallStrainElement3D4N", mDrainedUPwSmallStrainElement3D4N )
    KRATOS_REGISTER_ELEMENT( "DrainedUPwSmallStrainElement3D8N", mDrainedUPwSmallStrainElement3D8N )

    // Undrained small strain elements
    KRATOS_REGISTER_ELEMENT( "UndrainedUPwSmallStrainElement2D3N", mUndrainedUPwSmallStrainElement2D3N )
    KRATOS_REGISTER_ELEMENT( "UndrainedUPwSmallStrainElement2D4N", mUndrainedUPwSmallStrainElement2D4N )
    KRATOS_REGISTER_ELEMENT( "UndrainedUPwSmallStrainElement3D4N", mUndrainedUPwSmallStrainElement3D4N )
    KRATOS_REGISTER_ELEMENT( "UndrainedUPwSmallStrainElement3D8N", mUndrainedUPwSmallStrainElement3D8N )

    // Small strain FIC elements
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainFICElement2D3N", mUPwSmallStrainFICElement2D3N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainFICElement2D4N", mUPwSmallStrainFICElement2D4N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainFICElement3D4N", mUPwSmallStrainFICElement3D4N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainFICElement3D8N", mUPwSmallStrainFICElement3D8N )

    // Small strain different order elements
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwDiffOrderElement2D6N", mSmallStrainUPwDiffOrderElement2D6N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwDiffOrderElement2D8N", mSmallStrainUPwDiffOrderElement2D8N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwDiffOrderElement2D9N", mSmallStrainUPwDiffOrderElement2D9N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwDiffOrderElement3D10N", mSmallStrainUPwDiffOrderElement3D10N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwDiffOrderElement3D20N", mSmallStrainUPwDiffOrderElement3D20N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwDiffOrderElement3D27N", mSmallStrainUPwDiffOrderElement3D27N )

    // small strain axisymmtric elements:
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainAxisymmetricElement2D3N", mUPwSmallStrainAxisymmetricElement2D3N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainAxisymmetricElement2D4N", mUPwSmallStrainAxisymmetricElement2D4N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainAxisymmetricElement2D6N", mUPwSmallStrainAxisymmetricElement2D6N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainAxisymmetricElement2D8N", mUPwSmallStrainAxisymmetricElement2D8N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainAxisymmetricElement2D9N", mUPwSmallStrainAxisymmetricElement2D9N )

    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainAxisymmetricFICElement2D3N", mUPwSmallStrainAxisymmetricFICElement2D3N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainAxisymmetricFICElement2D4N", mUPwSmallStrainAxisymmetricFICElement2D4N )

    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwDiffOrderAxisymmetricElement2D6N", mSmallStrainUPwDiffOrderAxisymmetricElement2D6N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwDiffOrderAxisymmetricElement2D8N", mSmallStrainUPwDiffOrderAxisymmetricElement2D8N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwDiffOrderAxisymmetricElement2D9N", mSmallStrainUPwDiffOrderAxisymmetricElement2D9N )

    // Small strain interface elements
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainInterfaceElement2D4N", mUPwSmallStrainInterfaceElement2D4N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainInterfaceElement3D6N", mUPwSmallStrainInterfaceElement3D6N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainInterfaceElement3D8N", mUPwSmallStrainInterfaceElement3D8N )

    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainLinkInterfaceElement2D4N", mUPwSmallStrainLinkInterfaceElement2D4N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainLinkInterfaceElement3D6N", mUPwSmallStrainLinkInterfaceElement3D6N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainLinkInterfaceElement3D8N", mUPwSmallStrainLinkInterfaceElement3D8N )

    // Updated-Lagranian elements
    KRATOS_REGISTER_ELEMENT( "UPwUpdatedLagrangianElement2D3N", mUPwUpdatedLagrangianElement2D3N )
    KRATOS_REGISTER_ELEMENT( "UPwUpdatedLagrangianElement2D4N", mUPwUpdatedLagrangianElement2D4N )
    KRATOS_REGISTER_ELEMENT( "UPwUpdatedLagrangianElement3D4N", mUPwUpdatedLagrangianElement3D4N )
    KRATOS_REGISTER_ELEMENT( "UPwUpdatedLagrangianElement3D8N", mUPwUpdatedLagrangianElement3D8N )

    KRATOS_REGISTER_ELEMENT( "UPwUpdatedLagrangianElement2D6N", mUPwUpdatedLagrangianElement2D6N )
    KRATOS_REGISTER_ELEMENT( "UPwUpdatedLagrangianElement2D8N", mUPwUpdatedLagrangianElement2D8N )
    KRATOS_REGISTER_ELEMENT( "UPwUpdatedLagrangianElement2D9N", mUPwUpdatedLagrangianElement2D9N )
    KRATOS_REGISTER_ELEMENT( "UPwUpdatedLagrangianElement3D10N", mUPwUpdatedLagrangianElement3D10N )
    KRATOS_REGISTER_ELEMENT( "UPwUpdatedLagrangianElement3D20N", mUPwUpdatedLagrangianElement3D20N )
    KRATOS_REGISTER_ELEMENT( "UPwUpdatedLagrangianElement3D27N", mUPwUpdatedLagrangianElement3D27N )

    KRATOS_REGISTER_ELEMENT( "UPwUpdatedLagrangianFICElement2D3N", mUPwUpdatedLagrangianFICElement2D3N )
    KRATOS_REGISTER_ELEMENT( "UPwUpdatedLagrangianFICElement2D4N", mUPwUpdatedLagrangianFICElement2D4N )
    KRATOS_REGISTER_ELEMENT( "UPwUpdatedLagrangianFICElement3D4N", mUPwUpdatedLagrangianFICElement3D4N )
    KRATOS_REGISTER_ELEMENT( "UPwUpdatedLagrangianFICElement3D8N", mUPwUpdatedLagrangianFICElement3D8N )

    // Updated-Lagrangian different order elements
    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianUPwDiffOrderElement2D6N", mUpdatedLagrangianUPwDiffOrderElement2D6N )
    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianUPwDiffOrderElement2D8N", mUpdatedLagrangianUPwDiffOrderElement2D8N )
    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianUPwDiffOrderElement2D9N", mUpdatedLagrangianUPwDiffOrderElement2D9N )
    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianUPwDiffOrderElement3D10N", mUpdatedLagrangianUPwDiffOrderElement3D10N )
    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianUPwDiffOrderElement3D20N", mUpdatedLagrangianUPwDiffOrderElement3D20N )
    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianUPwDiffOrderElement3D27N", mUpdatedLagrangianUPwDiffOrderElement3D27N )


    // Register geo structural elements
    KRATOS_REGISTER_ELEMENT("GeoTrussElement2D2N", mGeoTrussElement2D2N)
    KRATOS_REGISTER_ELEMENT("GeoTrussElement3D2N", mGeoTrussElement3D2N)
    KRATOS_REGISTER_ELEMENT("GeoLinearTrussElement2D2N", mGeoLinearTrussElement2D2N)
    KRATOS_REGISTER_ELEMENT("GeoLinearTrussElement3D2N", mGeoLinearTrussElement3D2N)
    KRATOS_REGISTER_ELEMENT("GeoCableElement2D2N", mGeoCableElement2D2N)
    KRATOS_REGISTER_ELEMENT("GeoCableElement3D2N", mGeoCableElement3D2N)
    KRATOS_REGISTER_ELEMENT("GeoCrBeamElement3D2N", mGeoCrBeamElement3D2N)
    KRATOS_REGISTER_ELEMENT("GeoCrBeamElement2D2N", mGeoCrBeamElement2D2N)
    KRATOS_REGISTER_ELEMENT("GeoCrBeamElementLinear2D2N", mGeoCrBeamElementLinear2D2N)
    KRATOS_REGISTER_ELEMENT("GeoCrBeamElementLinear3D2N", mGeoCrBeamElementLinear3D2N)

    //Register Conditions
    KRATOS_REGISTER_CONDITION( "UPwForceCondition2D1N", mUPwForceCondition2D1N )
    KRATOS_REGISTER_CONDITION( "UPwForceCondition3D1N", mUPwForceCondition3D1N )
    KRATOS_REGISTER_CONDITION( "UPwFaceLoadCondition2D2N", mUPwFaceLoadCondition2D2N )
    KRATOS_REGISTER_CONDITION( "UPwFaceLoadCondition3D3N", mUPwFaceLoadCondition3D3N )
    KRATOS_REGISTER_CONDITION( "UPwFaceLoadCondition3D4N", mUPwFaceLoadCondition3D4N )
    KRATOS_REGISTER_CONDITION( "UPwNormalFaceLoadCondition2D2N", mUPwNormalFaceLoadCondition2D2N )
    KRATOS_REGISTER_CONDITION( "UPwNormalFaceLoadCondition3D3N", mUPwNormalFaceLoadCondition3D3N )
    KRATOS_REGISTER_CONDITION( "UpwNormalFaceLoadCondition3D4N", mUPwNormalFaceLoadCondition3D4N )
    KRATOS_REGISTER_CONDITION( "UPwNormalFluxCondition2D2N", mUPwNormalFluxCondition2D2N )
    KRATOS_REGISTER_CONDITION( "UPwNormalFluxCondition3D3N", mUPwNormalFluxCondition3D3N )
    KRATOS_REGISTER_CONDITION( "UPwNormalFluxCondition3D4N", mUPwNormalFluxCondition3D4N )

    KRATOS_REGISTER_CONDITION( "UPwFaceLoadCondition2D3N", mUPwFaceLoadCondition2D3N )

    KRATOS_REGISTER_CONDITION( "UPwFaceLoadInterfaceCondition2D2N", mUPwFaceLoadInterfaceCondition2D2N )
    KRATOS_REGISTER_CONDITION( "UPwFaceLoadInterfaceCondition3D4N", mUPwFaceLoadInterfaceCondition3D4N )
    KRATOS_REGISTER_CONDITION( "UPwNormalFluxInterfaceCondition2D2N", mUPwNormalFluxInterfaceCondition2D2N )
    KRATOS_REGISTER_CONDITION( "UPwNormalFluxInterfaceCondition3D4N", mUPwNormalFluxInterfaceCondition3D4N )

    KRATOS_REGISTER_CONDITION( "UPwNormalFluxFICCondition2D2N", mUPwNormalFluxFICCondition2D2N )
    KRATOS_REGISTER_CONDITION( "UPwNormalFluxFICCondition3D3N", mUPwNormalFluxFICCondition3D3N )
    KRATOS_REGISTER_CONDITION( "UPwNormalFluxFICCondition3D4N", mUPwNormalFluxFICCondition3D4N )

    KRATOS_REGISTER_CONDITION( "LineLoadDiffOrderCondition2D3N", mLineLoadDiffOrderCondition2D3N )
    KRATOS_REGISTER_CONDITION( "LineNormalLoadDiffOrderCondition2D3N", mLineNormalLoadDiffOrderCondition2D3N )
    KRATOS_REGISTER_CONDITION( "LineNormalFluidFluxDiffOrderCondition2D3N", mLineNormalFluidFluxDiffOrderCondition2D3N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadDiffOrderCondition3D6N", mSurfaceLoadDiffOrderCondition3D6N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadDiffOrderCondition3D8N", mSurfaceLoadDiffOrderCondition3D8N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadDiffOrderCondition3D9N", mSurfaceLoadDiffOrderCondition3D9N )
    KRATOS_REGISTER_CONDITION( "SurfaceNormalLoadDiffOrderCondition3D6N", mSurfaceNormalLoadDiffOrderCondition3D6N )
    KRATOS_REGISTER_CONDITION( "SurfaceNormalLoadDiffOrderCondition3D8N", mSurfaceNormalLoadDiffOrderCondition3D8N )
    KRATOS_REGISTER_CONDITION( "SurfaceNormalLoadDiffOrderCondition3D9N", mSurfaceNormalLoadDiffOrderCondition3D9N )
    KRATOS_REGISTER_CONDITION( "SurfaceNormalFluidFluxDiffOrderCondition3D6N", mSurfaceNormalFluidFluxDiffOrderCondition3D6N )
    KRATOS_REGISTER_CONDITION( "SurfaceNormalFluidFluxDiffOrderCondition3D8N", mSurfaceNormalFluidFluxDiffOrderCondition3D8N )
    KRATOS_REGISTER_CONDITION( "SurfaceNormalFluidFluxDiffOrderCondition3D9N", mSurfaceNormalFluidFluxDiffOrderCondition3D9N )

    //Register Constitutive Laws
    KRATOS_REGISTER_CONSTITUTIVE_LAW("BilinearCohesive3DLaw",           mBilinearCohesive3DLaw)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("BilinearCohesive2DLaw",           mBilinearCohesive2DLaw)

    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticPlaneStrainK02DLaw",  mLinearPlaneStrainK0Law)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticK03DLaw",             mElasticIsotropicK03DLaw)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("GeoLinearElasticPlaneStrain2DLaw", mLinearElasticPlaneStrain2DLaw)

    KRATOS_REGISTER_CONSTITUTIVE_LAW("GeoLinearElasticPlaneStress2DLaw", mLinearElasticPlaneStress2DLaw)

    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainUDSM3DLaw",            mSmallStrainUDSM3DLaw)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainUDSM2DPlaneStrainLaw", mSmallStrainUDSM2DPlaneStrainLaw)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainUDSM2DInterfaceLaw",   mSmallStrainUDSM2DInterfaceLaw)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainUDSM3DInterfaceLaw",   mSmallStrainUDSM3DInterfaceLaw)

    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainUMAT3DLaw",            mSmallStrainUMAT3DLaw)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainUMAT2DPlaneStrainLaw", mSmallStrainUMAT2DPlaneStrainLaw)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainUMAT2DInterfaceLaw",   mSmallStrainUMAT2DInterfaceLaw)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SmallStrainUMAT3DInterfaceLaw",   mSmallStrainUMAT3DInterfaceLaw)

    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElastic2DInterfaceLaw",   mLinearElastic2DInterfaceLaw)
    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElastic3DInterfaceLaw",   mLinearElastic3DInterfaceLaw)

    //Register Variables
    KRATOS_REGISTER_VARIABLE( VELOCITY_COEFFICIENT )
    KRATOS_REGISTER_VARIABLE( DT_PRESSURE_COEFFICIENT )

    KRATOS_REGISTER_VARIABLE( DT_WATER_PRESSURE )
    KRATOS_REGISTER_VARIABLE( NORMAL_FLUID_FLUX )

    KRATOS_REGISTER_VARIABLE( HYDRAULIC_HEAD )

    KRATOS_REGISTER_VARIABLE( HYDRAULIC_DISCHARGE )

    KRATOS_REGISTER_VARIABLE( DENSITY_SOLID )
    KRATOS_REGISTER_VARIABLE( BULK_MODULUS_SOLID )
    KRATOS_REGISTER_VARIABLE( BULK_MODULUS_FLUID )

    KRATOS_REGISTER_VARIABLE( K0_MAIN_DIRECTION )
    KRATOS_REGISTER_VARIABLE( K0_VALUE_XX )
    KRATOS_REGISTER_VARIABLE( K0_VALUE_YY )
    KRATOS_REGISTER_VARIABLE( K0_VALUE_ZZ )

    KRATOS_REGISTER_VARIABLE( PERMEABILITY_XX )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_YY )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_ZZ )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_XY )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_YZ )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_ZX )

    KRATOS_REGISTER_VARIABLE( MINIMUM_JOINT_WIDTH )
    KRATOS_REGISTER_VARIABLE( TRANSVERSAL_PERMEABILITY )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( FLUID_FLUX_VECTOR )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LOCAL_FLUID_FLUX_VECTOR )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LOCAL_STRESS_VECTOR )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LOCAL_RELATIVE_DISPLACEMENT_VECTOR )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_MATRIX )
    KRATOS_REGISTER_VARIABLE( LOCAL_PERMEABILITY_MATRIX )

    KRATOS_REGISTER_VARIABLE( CRITICAL_DISPLACEMENT )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( TOTAL_DISPLACEMENT )

    KRATOS_REGISTER_VARIABLE( IS_CONVERGED )

    KRATOS_REGISTER_VARIABLE( TOTAL_STRESS_TENSOR )
    KRATOS_REGISTER_VARIABLE( TOTAL_STRESS_VECTOR )

    KRATOS_REGISTER_VARIABLE( CAUCHY_STRAIN_TENSOR )
    KRATOS_REGISTER_VARIABLE( CAUCHY_STRAIN_VECTOR )

    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE )
    KRATOS_REGISTER_VARIABLE( ARC_LENGTH_LAMBDA )
    KRATOS_REGISTER_VARIABLE( ARC_LENGTH_RADIUS_FACTOR )

    KRATOS_REGISTER_VARIABLE( TIME_UNIT_CONVERTER )

    KRATOS_REGISTER_VARIABLE( LOCAL_EQUIVALENT_STRAIN )
    KRATOS_REGISTER_VARIABLE( NONLOCAL_EQUIVALENT_STRAIN )

    KRATOS_REGISTER_VARIABLE( JOINT_WIDTH )

    KRATOS_REGISTER_VARIABLE( NODAL_SMOOTHING )
    KRATOS_REGISTER_VARIABLE( NODAL_CAUCHY_STRESS_TENSOR )
    KRATOS_REGISTER_VARIABLE( NODAL_DAMAGE_VARIABLE )
    KRATOS_REGISTER_VARIABLE( NODAL_JOINT_AREA )
    KRATOS_REGISTER_VARIABLE( NODAL_JOINT_WIDTH )
    KRATOS_REGISTER_VARIABLE( NODAL_JOINT_DAMAGE )

    KRATOS_REGISTER_VARIABLE( BIOT_COEFFICIENT )

    KRATOS_REGISTER_VARIABLE( RESET_DISPLACEMENTS )
    KRATOS_REGISTER_VARIABLE( CONSIDER_GEOMETRIC_STIFFNESS )

    KRATOS_REGISTER_VARIABLE( CONSIDER_GAP_CLOSURE )

    KRATOS_REGISTER_VARIABLE( USE_CONSISTENT_MASS_MATRIX)

    KRATOS_REGISTER_VARIABLE( IGNORE_UNDRAINED )

    KRATOS_REGISTER_VARIABLE( SATURATED_SATURATION )
    KRATOS_REGISTER_VARIABLE( RESIDUAL_SATURATION )
    KRATOS_REGISTER_VARIABLE( VAN_GENUCHTEN_AIR_ENTRY_PRESSURE )
    KRATOS_REGISTER_VARIABLE( VAN_GENUCHTEN_GN )
    KRATOS_REGISTER_VARIABLE( VAN_GENUCHTEN_GL )
    KRATOS_REGISTER_VARIABLE( MINIMUM_RELATIVE_PERMEABILITY )

    KRATOS_REGISTER_VARIABLE( RETENTION_LAW )
    KRATOS_REGISTER_VARIABLE( DEGREE_OF_SATURATION )
    KRATOS_REGISTER_VARIABLE( EFFECTIVE_SATURATION )
    KRATOS_REGISTER_VARIABLE( BISHOP_COEFICIENT )
    KRATOS_REGISTER_VARIABLE( DERIVATIVE_OF_SATURATION )
    KRATOS_REGISTER_VARIABLE( RELATIVE_PERMEABILITY )

    // UDSM
    KRATOS_REGISTER_VARIABLE( UDSM_NAME )       // Also for UMAT
    KRATOS_REGISTER_VARIABLE( UDSM_NUMBER )
    KRATOS_REGISTER_VARIABLE( IS_FORTRAN_UDSM ) // Also for UMAT

    KRATOS_REGISTER_VARIABLE( UMAT_PARAMETERS )

    KRATOS_REGISTER_VARIABLE(NUMBER_OF_UMAT_STATE_VARIABLES)
    KRATOS_REGISTER_VARIABLE(NUMBER_OF_UMAT_PARAMETERS)

    KRATOS_REGISTER_VARIABLE( STATE_VARIABLES )

    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_1 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_2 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_3 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_4 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_5 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_6 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_7 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_8 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_9 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_10 )

    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_11 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_12 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_13 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_14 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_15 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_16 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_17 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_18 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_19 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_20 )

    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_21 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_22 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_23 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_24 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_25 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_26 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_27 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_28 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_29 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_30 )

    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_31 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_32 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_33 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_34 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_35 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_36 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_37 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_38 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_39 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_40 )

    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_41 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_42 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_43 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_44 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_45 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_46 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_47 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_48 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_49 )
    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE_50 )

   }
}  // namespace Kratos.
