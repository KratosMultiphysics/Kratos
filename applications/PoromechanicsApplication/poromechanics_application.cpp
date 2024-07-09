
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//


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
#include "poromechanics_application.h"

namespace Kratos
{

KratosPoromechanicsApplication::KratosPoromechanicsApplication()
    : KratosApplication("PoromechanicsApplication"),

    mUPlSmallStrainElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node >( Element::GeometryType::PointsArrayType(3)))),
    mUPlSmallStrainElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node >( Element::GeometryType::PointsArrayType(4)))),
    mUPlSmallStrainElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node >( Element::GeometryType::PointsArrayType(4)))),
    mUPlSmallStrainElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node >( Element::GeometryType::PointsArrayType(8)))),

    mUPlSmallStrainInterfaceElement2D4N( 0, Element::GeometryType::Pointer( new QuadrilateralInterface2D4 <Node >( Element::GeometryType::PointsArrayType(4)))),
    mUPlSmallStrainInterfaceElement3D6N( 0, Element::GeometryType::Pointer( new PrismInterface3D6 <Node >( Element::GeometryType::PointsArrayType(6)))),
    mUPlSmallStrainInterfaceElement3D8N( 0, Element::GeometryType::Pointer( new HexahedraInterface3D8 <Node >( Element::GeometryType::PointsArrayType(8)))),

    mUPlSmallStrainLinkInterfaceElement2D4N( 0, Element::GeometryType::Pointer( new QuadrilateralInterface2D4 <Node >( Element::GeometryType::PointsArrayType(4)))),
    mUPlSmallStrainLinkInterfaceElement3D6N( 0, Element::GeometryType::Pointer( new PrismInterface3D6 <Node >( Element::GeometryType::PointsArrayType(6)))),
    mUPlSmallStrainLinkInterfaceElement3D8N( 0, Element::GeometryType::Pointer( new HexahedraInterface3D8 <Node >( Element::GeometryType::PointsArrayType(8)))),

    mUPlSmallStrainFICElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node >( Element::GeometryType::PointsArrayType(3)))),
    mUPlSmallStrainFICElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node >( Element::GeometryType::PointsArrayType(4)))),
    mUPlSmallStrainFICElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node >( Element::GeometryType::PointsArrayType(4)))),
    mUPlSmallStrainFICElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node >( Element::GeometryType::PointsArrayType(8)))),

    mSmallStrainUPlDiffOrderElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node >( Element::GeometryType::PointsArrayType(6)))),
    mSmallStrainUPlDiffOrderElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node >( Element::GeometryType::PointsArrayType(8)))),
    mSmallStrainUPlDiffOrderElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node >( Element::GeometryType::PointsArrayType(9)))),
    mSmallStrainUPlDiffOrderElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node >( Element::GeometryType::PointsArrayType(10)))),
    mSmallStrainUPlDiffOrderElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node >( Element::GeometryType::PointsArrayType(20)))),
    mSmallStrainUPlDiffOrderElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node >( Element::GeometryType::PointsArrayType(27)))),

    mUPlPgSmallStrainElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node >( Element::GeometryType::PointsArrayType(3)))),
    mUPlPgSmallStrainElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node >( Element::GeometryType::PointsArrayType(4)))),
    mUPlPgSmallStrainElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node >( Element::GeometryType::PointsArrayType(4)))),
    mUPlPgSmallStrainElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node >( Element::GeometryType::PointsArrayType(8)))),


    mUPlForceCondition2D1N( 0, Condition::GeometryType::Pointer( new Point2D<Node >( Condition::GeometryType::PointsArrayType(1)))),
    mUPlForceCondition3D1N( 0, Condition::GeometryType::Pointer( new Point3D<Node >( Condition::GeometryType::PointsArrayType(1)))),
    mUPlFaceLoadCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node >( Condition::GeometryType::PointsArrayType(2)))),
    mUPlFaceLoadCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node >( Condition::GeometryType::PointsArrayType(3)))),
    mUPlFaceLoadCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node >( Condition::GeometryType::PointsArrayType(4)))),
    mUPlNormalFaceLoadCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node >( Condition::GeometryType::PointsArrayType(2)))),
    mUPlNormalFaceLoadCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node >( Condition::GeometryType::PointsArrayType(3)))),
    mUPlNormalFaceLoadCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node >( Condition::GeometryType::PointsArrayType(4)))),
    mUPlLiquidDischargeCondition2D1N( 0, Condition::GeometryType::Pointer( new Point2D<Node >( Condition::GeometryType::PointsArrayType(1)))),
    mUPlLiquidDischargeCondition3D1N( 0, Condition::GeometryType::Pointer( new Point3D<Node >( Condition::GeometryType::PointsArrayType(1)))),
    mUPlNormalLiquidFluxCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node >( Condition::GeometryType::PointsArrayType(2)))),
    mUPlNormalLiquidFluxCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node >( Condition::GeometryType::PointsArrayType(3)))),
    mUPlNormalLiquidFluxCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node >( Condition::GeometryType::PointsArrayType(4)))),

    mUPlFaceLoadInterfaceCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node >( Condition::GeometryType::PointsArrayType(2)))),
    mUPlFaceLoadInterfaceCondition3D4N( 0, Condition::GeometryType::Pointer( new QuadrilateralInterface3D4 <Node >( Condition::GeometryType::PointsArrayType(4)))),
    mUPlNormalLiquidFluxInterfaceCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node >( Condition::GeometryType::PointsArrayType(2)))),
    mUPlNormalLiquidFluxInterfaceCondition3D4N( 0, Condition::GeometryType::Pointer( new QuadrilateralInterface3D4 <Node >( Condition::GeometryType::PointsArrayType(4)))),

    mUPlNormalLiquidFluxFICCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node >( Condition::GeometryType::PointsArrayType(2)))),
    mUPlNormalLiquidFluxFICCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node >( Condition::GeometryType::PointsArrayType(3)))),
    mUPlNormalLiquidFluxFICCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node >( Condition::GeometryType::PointsArrayType(4)))),

    mLineLoadDiffOrderCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3<Node >( Condition::GeometryType::PointsArrayType(3)))),
    mLineNormalLoadDiffOrderCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3<Node >( Condition::GeometryType::PointsArrayType(3)))),
    mLineNormalLiquidFluxDiffOrderCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3<Node >( Condition::GeometryType::PointsArrayType(3)))),
    mSurfaceLoadDiffOrderCondition3D6N( 0, Condition::GeometryType::Pointer( new Triangle3D6 <Node >( Condition::GeometryType::PointsArrayType(6)))),
    mSurfaceLoadDiffOrderCondition3D8N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D8 <Node >( Condition::GeometryType::PointsArrayType(8)))),
    mSurfaceLoadDiffOrderCondition3D9N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D9 <Node >( Condition::GeometryType::PointsArrayType(9)))),
    mSurfaceNormalLoadDiffOrderCondition3D6N( 0, Condition::GeometryType::Pointer( new Triangle3D6 <Node >( Condition::GeometryType::PointsArrayType(6)))),
    mSurfaceNormalLoadDiffOrderCondition3D8N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D8 <Node >( Condition::GeometryType::PointsArrayType(8)))),
    mSurfaceNormalLoadDiffOrderCondition3D9N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D9 <Node >( Condition::GeometryType::PointsArrayType(9)))),
    mSurfaceNormalLiquidFluxDiffOrderCondition3D6N( 0, Condition::GeometryType::Pointer( new Triangle3D6 <Node >( Condition::GeometryType::PointsArrayType(6)))),
    mSurfaceNormalLiquidFluxDiffOrderCondition3D8N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D8 <Node >( Condition::GeometryType::PointsArrayType(8)))),
    mSurfaceNormalLiquidFluxDiffOrderCondition3D9N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D9 <Node >( Condition::GeometryType::PointsArrayType(9)))),

    mUPlPgForceCondition2D1N( 0, Condition::GeometryType::Pointer( new Point2D<Node >( Condition::GeometryType::PointsArrayType(1)))),
    mUPlPgForceCondition3D1N( 0, Condition::GeometryType::Pointer( new Point3D<Node >( Condition::GeometryType::PointsArrayType(1)))),
    mUPlPgFaceLoadCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node >( Condition::GeometryType::PointsArrayType(2)))),
    mUPlPgFaceLoadCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node >( Condition::GeometryType::PointsArrayType(3)))),
    mUPlPgFaceLoadCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node >( Condition::GeometryType::PointsArrayType(4)))),
    mUPlPgNormalFaceLoadCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node >( Condition::GeometryType::PointsArrayType(2)))),
    mUPlPgNormalFaceLoadCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node >( Condition::GeometryType::PointsArrayType(3)))),
    mUPlPgNormalFaceLoadCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node >( Condition::GeometryType::PointsArrayType(4)))),
    mUPlPgLiquidDischargeCondition2D1N( 0, Condition::GeometryType::Pointer( new Point2D<Node >( Condition::GeometryType::PointsArrayType(1)))),
    mUPlPgLiquidDischargeCondition3D1N( 0, Condition::GeometryType::Pointer( new Point3D<Node >( Condition::GeometryType::PointsArrayType(1)))),
    mUPlPgGasDischargeCondition2D1N( 0, Condition::GeometryType::Pointer( new Point2D<Node >( Condition::GeometryType::PointsArrayType(1)))),
    mUPlPgGasDischargeCondition3D1N( 0, Condition::GeometryType::Pointer( new Point3D<Node >( Condition::GeometryType::PointsArrayType(1)))),
    mUPlPgNormalLiquidFluxCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node >( Condition::GeometryType::PointsArrayType(2)))),
    mUPlPgNormalLiquidFluxCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node >( Condition::GeometryType::PointsArrayType(3)))),
    mUPlPgNormalLiquidFluxCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node >( Condition::GeometryType::PointsArrayType(4)))),
    mUPlPgNormalGasFluxCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node >( Condition::GeometryType::PointsArrayType(2)))),
    mUPlPgNormalGasFluxCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node >( Condition::GeometryType::PointsArrayType(3)))),
    mUPlPgNormalGasFluxCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node >( Condition::GeometryType::PointsArrayType(4))))

    {}

void KratosPoromechanicsApplication::Register()
{
    KRATOS_INFO("") << "Initializing KratosPoromechanicsApplication... " << std::endl;

    //Register Elements
    KRATOS_REGISTER_ELEMENT( "UPlSmallStrainElement2D3N", mUPlSmallStrainElement2D3N )
    KRATOS_REGISTER_ELEMENT( "UPlSmallStrainElement2D4N", mUPlSmallStrainElement2D4N )
    KRATOS_REGISTER_ELEMENT( "UPlSmallStrainElement3D4N", mUPlSmallStrainElement3D4N )
    KRATOS_REGISTER_ELEMENT( "UPlSmallStrainElement3D8N", mUPlSmallStrainElement3D8N )

    KRATOS_REGISTER_ELEMENT( "UPlSmallStrainInterfaceElement2D4N", mUPlSmallStrainInterfaceElement2D4N )
    KRATOS_REGISTER_ELEMENT( "UPlSmallStrainInterfaceElement3D6N", mUPlSmallStrainInterfaceElement3D6N )
    KRATOS_REGISTER_ELEMENT( "UPlSmallStrainInterfaceElement3D8N", mUPlSmallStrainInterfaceElement3D8N )

    KRATOS_REGISTER_ELEMENT( "UPlSmallStrainLinkInterfaceElement2D4N", mUPlSmallStrainLinkInterfaceElement2D4N )
    KRATOS_REGISTER_ELEMENT( "UPlSmallStrainLinkInterfaceElement3D6N", mUPlSmallStrainLinkInterfaceElement3D6N )
    KRATOS_REGISTER_ELEMENT( "UPlSmallStrainLinkInterfaceElement3D8N", mUPlSmallStrainLinkInterfaceElement3D8N )

    KRATOS_REGISTER_ELEMENT( "UPlSmallStrainFICElement2D3N", mUPlSmallStrainFICElement2D3N )
    KRATOS_REGISTER_ELEMENT( "UPlSmallStrainFICElement2D4N", mUPlSmallStrainFICElement2D4N )
    KRATOS_REGISTER_ELEMENT( "UPlSmallStrainFICElement3D4N", mUPlSmallStrainFICElement3D4N )
    KRATOS_REGISTER_ELEMENT( "UPlSmallStrainFICElement3D8N", mUPlSmallStrainFICElement3D8N )

    KRATOS_REGISTER_ELEMENT( "SmallStrainUPlDiffOrderElement2D6N", mSmallStrainUPlDiffOrderElement2D6N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPlDiffOrderElement2D8N", mSmallStrainUPlDiffOrderElement2D8N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPlDiffOrderElement2D9N", mSmallStrainUPlDiffOrderElement2D9N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPlDiffOrderElement3D10N", mSmallStrainUPlDiffOrderElement3D10N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPlDiffOrderElement3D20N", mSmallStrainUPlDiffOrderElement3D20N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPlDiffOrderElement3D27N", mSmallStrainUPlDiffOrderElement3D27N )

    KRATOS_REGISTER_ELEMENT( "UPlPgSmallStrainElement2D3N", mUPlPgSmallStrainElement2D3N )
    KRATOS_REGISTER_ELEMENT( "UPlPgSmallStrainElement2D4N", mUPlPgSmallStrainElement2D4N )
    KRATOS_REGISTER_ELEMENT( "UPlPgSmallStrainElement3D4N", mUPlPgSmallStrainElement3D4N )
    KRATOS_REGISTER_ELEMENT( "UPlPgSmallStrainElement3D8N", mUPlPgSmallStrainElement3D8N )

    //Register Conditions
    KRATOS_REGISTER_CONDITION( "UPlForceCondition2D1N", mUPlForceCondition2D1N )
    KRATOS_REGISTER_CONDITION( "UPlForceCondition3D1N", mUPlForceCondition3D1N )
    KRATOS_REGISTER_CONDITION( "UPlFaceLoadCondition2D2N", mUPlFaceLoadCondition2D2N )
    KRATOS_REGISTER_CONDITION( "UPlFaceLoadCondition3D3N", mUPlFaceLoadCondition3D3N )
    KRATOS_REGISTER_CONDITION( "UPlFaceLoadCondition3D4N", mUPlFaceLoadCondition3D4N )
    KRATOS_REGISTER_CONDITION( "UPlNormalFaceLoadCondition2D2N", mUPlNormalFaceLoadCondition2D2N )
    KRATOS_REGISTER_CONDITION( "UPlNormalFaceLoadCondition3D3N", mUPlNormalFaceLoadCondition3D3N )
    KRATOS_REGISTER_CONDITION( "UPlNormalFaceLoadCondition3D4N", mUPlNormalFaceLoadCondition3D4N )
    KRATOS_REGISTER_CONDITION( "UPlLiquidDischargeCondition2D1N", mUPlLiquidDischargeCondition2D1N )
    KRATOS_REGISTER_CONDITION( "UPlLiquidDischargeCondition3D1N", mUPlLiquidDischargeCondition3D1N )
    KRATOS_REGISTER_CONDITION( "UPlNormalLiquidFluxCondition2D2N", mUPlNormalLiquidFluxCondition2D2N )
    KRATOS_REGISTER_CONDITION( "UPlNormalLiquidFluxCondition3D3N", mUPlNormalLiquidFluxCondition3D3N )
    KRATOS_REGISTER_CONDITION( "UPlNormalLiquidFluxCondition3D4N", mUPlNormalLiquidFluxCondition3D4N )

    KRATOS_REGISTER_CONDITION( "UPlFaceLoadInterfaceCondition2D2N", mUPlFaceLoadInterfaceCondition2D2N )
    KRATOS_REGISTER_CONDITION( "UPlFaceLoadInterfaceCondition3D4N", mUPlFaceLoadInterfaceCondition3D4N )
    KRATOS_REGISTER_CONDITION( "UPlNormalLiquidFluxInterfaceCondition2D2N", mUPlNormalLiquidFluxInterfaceCondition2D2N )
    KRATOS_REGISTER_CONDITION( "UPlNormalLiquidFluxInterfaceCondition3D4N", mUPlNormalLiquidFluxInterfaceCondition3D4N )

    KRATOS_REGISTER_CONDITION( "UPlNormalLiquidFluxFICCondition2D2N", mUPlNormalLiquidFluxFICCondition2D2N )
    KRATOS_REGISTER_CONDITION( "UPlNormalLiquidFluxFICCondition3D3N", mUPlNormalLiquidFluxFICCondition3D3N )
    KRATOS_REGISTER_CONDITION( "UPlNormalLiquidFluxFICCondition3D4N", mUPlNormalLiquidFluxFICCondition3D4N )

    KRATOS_REGISTER_CONDITION( "LineLoadDiffOrderCondition2D3N", mLineLoadDiffOrderCondition2D3N )
    KRATOS_REGISTER_CONDITION( "LineNormalLoadDiffOrderCondition2D3N", mLineNormalLoadDiffOrderCondition2D3N )
    KRATOS_REGISTER_CONDITION( "LineNormalLiquidFluxDiffOrderCondition2D3N", mLineNormalLiquidFluxDiffOrderCondition2D3N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadDiffOrderCondition3D6N", mSurfaceLoadDiffOrderCondition3D6N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadDiffOrderCondition3D8N", mSurfaceLoadDiffOrderCondition3D8N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadDiffOrderCondition3D9N", mSurfaceLoadDiffOrderCondition3D9N )
    KRATOS_REGISTER_CONDITION( "SurfaceNormalLoadDiffOrderCondition3D6N", mSurfaceNormalLoadDiffOrderCondition3D6N )
    KRATOS_REGISTER_CONDITION( "SurfaceNormalLoadDiffOrderCondition3D8N", mSurfaceNormalLoadDiffOrderCondition3D8N )
    KRATOS_REGISTER_CONDITION( "SurfaceNormalLoadDiffOrderCondition3D9N", mSurfaceNormalLoadDiffOrderCondition3D9N )
    KRATOS_REGISTER_CONDITION( "SurfaceNormalLiquidFluxDiffOrderCondition3D6N", mSurfaceNormalLiquidFluxDiffOrderCondition3D6N )
    KRATOS_REGISTER_CONDITION( "SurfaceNormalLiquidFluxDiffOrderCondition3D8N", mSurfaceNormalLiquidFluxDiffOrderCondition3D8N )
    KRATOS_REGISTER_CONDITION( "SurfaceNormalLiquidFluxDiffOrderCondition3D9N", mSurfaceNormalLiquidFluxDiffOrderCondition3D9N )

    KRATOS_REGISTER_CONDITION( "UPlPgForceCondition2D1N", mUPlPgForceCondition2D1N )
    KRATOS_REGISTER_CONDITION( "UPlPgForceCondition3D1N", mUPlPgForceCondition3D1N )
    KRATOS_REGISTER_CONDITION( "UPlPgFaceLoadCondition2D2N", mUPlPgFaceLoadCondition2D2N )
    KRATOS_REGISTER_CONDITION( "UPlPgFaceLoadCondition3D3N", mUPlPgFaceLoadCondition3D3N )
    KRATOS_REGISTER_CONDITION( "UPlPgFaceLoadCondition3D4N", mUPlPgFaceLoadCondition3D4N )
    KRATOS_REGISTER_CONDITION( "UPlPgNormalFaceLoadCondition2D2N", mUPlPgNormalFaceLoadCondition2D2N )
    KRATOS_REGISTER_CONDITION( "UPlPgNormalFaceLoadCondition3D3N", mUPlPgNormalFaceLoadCondition3D3N )
    KRATOS_REGISTER_CONDITION( "UPlPgNormalFaceLoadCondition3D4N", mUPlPgNormalFaceLoadCondition3D4N )
    KRATOS_REGISTER_CONDITION( "UPlPgLiquidDischargeCondition2D1N", mUPlPgLiquidDischargeCondition2D1N )
    KRATOS_REGISTER_CONDITION( "UPlPgLiquidDischargeCondition3D1N", mUPlPgLiquidDischargeCondition3D1N )
    KRATOS_REGISTER_CONDITION( "UPlPgGasDischargeCondition2D1N", mUPlPgGasDischargeCondition2D1N )
    KRATOS_REGISTER_CONDITION( "UPlPgGasDischargeCondition3D1N", mUPlPgGasDischargeCondition3D1N )
    KRATOS_REGISTER_CONDITION( "UPlPgNormalLiquidFluxCondition2D2N", mUPlPgNormalLiquidFluxCondition2D2N )
    KRATOS_REGISTER_CONDITION( "UPlPgNormalLiquidFluxCondition3D3N", mUPlPgNormalLiquidFluxCondition3D3N )
    KRATOS_REGISTER_CONDITION( "UPlPgNormalLiquidFluxCondition3D4N", mUPlPgNormalLiquidFluxCondition3D4N )
    KRATOS_REGISTER_CONDITION( "UPlPgNormalGasFluxCondition2D2N", mUPlPgNormalGasFluxCondition2D2N )
    KRATOS_REGISTER_CONDITION( "UPlPgNormalGasFluxCondition3D3N", mUPlPgNormalGasFluxCondition3D3N )
    KRATOS_REGISTER_CONDITION( "UPlPgNormalGasFluxCondition3D4N", mUPlPgNormalGasFluxCondition3D4N )

    //Register Constitutive Laws
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ElastoPlasticMohrCoulombCohesive3DLaw",mElastoPlasticMohrCoulombCohesive3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ElastoPlasticMohrCoulombCohesive2DLaw",mElastoPlasticMohrCoulombCohesive2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ElastoPlasticModMohrCoulombCohesive3DLaw",mElastoPlasticModMohrCoulombCohesive3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ElastoPlasticModMohrCoulombCohesive2DLaw",mElastoPlasticModMohrCoulombCohesive2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("IsotropicDamageCohesive3DLaw",mIsotropicDamageCohesive3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("IsotropicDamageCohesive2DLaw",mIsotropicDamageCohesive2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ElasticCohesive3DLaw",mElasticCohesive3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ElasticCohesive2DLaw",mElasticCohesive2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("BilinearCohesive3DLaw",mBilinearCohesive3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("BilinearCohesive2DLaw",mBilinearCohesive2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ExponentialCohesive3DLaw",mExponentialCohesive3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ExponentialCohesive2DLaw",mExponentialCohesive2DLaw);

    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticSolid3DLaw", mLinearElastic3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticPlaneStrainSolid2DLaw", mLinearElasticPlaneStrain2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticPlaneStressSolid2DLaw", mLinearElasticPlaneStress2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticSolid3DLaw", mHyperElastic3DLaw);

    Serializer::Register("IsotropicDamageFlowRule", mIsotropicDamageFlowRule);
    Serializer::Register( "LocalDamageFlowRule", mLocalDamageFlowRule );
    Serializer::Register( "NonlocalDamageFlowRule", mNonlocalDamageFlowRule );

    Serializer::Register("SimoJuYieldCriterion", mSimoJuYieldCriterion);
    Serializer::Register("ModifiedMisesYieldCriterion", mModifiedMisesYieldCriterion);

    Serializer::Register("ExponentialDamageHardeningLaw", mExponentialDamageHardeningLaw);
    Serializer::Register("ModifiedExponentialDamageHardeningLaw", mModifiedExponentialDamageHardeningLaw);

    KRATOS_REGISTER_CONSTITUTIVE_LAW("SimoJuLocalDamage3DLaw",mSimoJuLocalDamage3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SimoJuLocalDamagePlaneStrain2DLaw",mSimoJuLocalDamagePlaneStrain2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SimoJuLocalDamagePlaneStress2DLaw",mSimoJuLocalDamagePlaneStress2DLaw);

    KRATOS_REGISTER_CONSTITUTIVE_LAW("SimoJuNonlocalDamage3DLaw",mSimoJuNonlocalDamage3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SimoJuNonlocalDamagePlaneStrain2DLaw",mSimoJuNonlocalDamagePlaneStrain2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("SimoJuNonlocalDamagePlaneStress2DLaw",mSimoJuNonlocalDamagePlaneStress2DLaw);

    KRATOS_REGISTER_CONSTITUTIVE_LAW("ModifiedMisesNonlocalDamage3DLaw",mModifiedMisesNonlocalDamage3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ModifiedMisesNonlocalDamagePlaneStrain2DLaw",mModifiedMisesNonlocalDamagePlaneStrain2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ModifiedMisesNonlocalDamagePlaneStress2DLaw",mModifiedMisesNonlocalDamagePlaneStress2DLaw);


    //Register Variables
    KRATOS_REGISTER_VARIABLE( VELOCITY_COEFFICIENT )
    KRATOS_REGISTER_VARIABLE( DT_LIQUID_PRESSURE_COEFFICIENT )
    KRATOS_REGISTER_VARIABLE( DT_GAS_PRESSURE_COEFFICIENT )

    KRATOS_REGISTER_VARIABLE( LIQUID_PRESSURE )
    KRATOS_REGISTER_VARIABLE( REACTION_LIQUID_PRESSURE )
    KRATOS_REGISTER_VARIABLE( GAS_PRESSURE )
    KRATOS_REGISTER_VARIABLE( REACTION_GAS_PRESSURE )
    KRATOS_REGISTER_VARIABLE( DT_LIQUID_PRESSURE )
    KRATOS_REGISTER_VARIABLE( DT_GAS_PRESSURE )
    KRATOS_REGISTER_VARIABLE( NORMAL_LIQUID_FLUX )
    KRATOS_REGISTER_VARIABLE( NORMAL_GAS_FLUX )
    KRATOS_REGISTER_VARIABLE( LIQUID_DISCHARGE )
    KRATOS_REGISTER_VARIABLE( GAS_DISCHARGE )
    KRATOS_REGISTER_VARIABLE( CAPILLARY_PRESSURE )
    KRATOS_REGISTER_VARIABLE( LIQUID_SATURATION_DEGREE )
    KRATOS_REGISTER_VARIABLE( GAS_SATURATION_DEGREE )
    KRATOS_REGISTER_VARIABLE( LIQUID_RELATIVE_PERMEABILITY )
    KRATOS_REGISTER_VARIABLE( GAS_RELATIVE_PERMEABILITY )

    KRATOS_REGISTER_VARIABLE( DENSITY_SOLID )
    KRATOS_REGISTER_VARIABLE( DENSITY_LIQUID )
    KRATOS_REGISTER_VARIABLE( DENSITY_GAS )
    KRATOS_REGISTER_VARIABLE( DYNAMIC_VISCOSITY_LIQUID )
    KRATOS_REGISTER_VARIABLE( DYNAMIC_VISCOSITY_GAS )
    KRATOS_REGISTER_VARIABLE( GAS_HENRY_SOLUBILITY_COEFF )
    KRATOS_REGISTER_VARIABLE( GAS_MOLAR_WEIGHT )
    KRATOS_REGISTER_VARIABLE( TORTUOSITY_COEFF )
    KRATOS_REGISTER_VARIABLE( GAS_DIFFUSION_COEFF )
    KRATOS_REGISTER_VARIABLE( LIQUID_DENSITY_REF_PRESSURE )
    KRATOS_REGISTER_VARIABLE( SATURATION_LAW_NAME )
    KRATOS_REGISTER_VARIABLE( RESIDUAL_LIQUID_SATURATION )
    KRATOS_REGISTER_VARIABLE( RESIDUAL_GAS_SATURATION )
    KRATOS_REGISTER_VARIABLE( GAS_ENTRY_PRESSURE )
    KRATOS_REGISTER_VARIABLE( PORE_SIZE_FACTOR )
    KRATOS_REGISTER_VARIABLE( MINIMUM_RELATIVE_PERMEABILITY )
    KRATOS_REGISTER_VARIABLE( BULK_MODULUS_SOLID )
    KRATOS_REGISTER_VARIABLE( BULK_MODULUS_LIQUID )
    KRATOS_REGISTER_VARIABLE( BULK_MODULUS_GAS )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_XX )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_YY )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_ZZ )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_XY )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_YZ )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_ZX )

    KRATOS_REGISTER_VARIABLE( NORMAL_STIFFNESS )
    KRATOS_REGISTER_VARIABLE( SHEAR_STIFFNESS )
    KRATOS_REGISTER_VARIABLE( PENALTY_STIFFNESS )
    KRATOS_REGISTER_VARIABLE( TENSILE_STRENGTH )
    KRATOS_REGISTER_VARIABLE( BETA_EQSTRAIN_SHEAR_FACTOR )
    KRATOS_REGISTER_VARIABLE( DAMAGE_EVOLUTION_LAW )
    KRATOS_REGISTER_VARIABLE( FRICTION_ANGLE )
    KRATOS_REGISTER_VARIABLE( DILATANCY_ANGLE )
    KRATOS_REGISTER_VARIABLE( COHESION )

    KRATOS_REGISTER_VARIABLE( INITIAL_JOINT_WIDTH )
    KRATOS_REGISTER_VARIABLE( TRANSVERSAL_PERMEABILITY_COEFFICIENT )
    KRATOS_REGISTER_VARIABLE( MID_PLANE_LIQUID_PRESSURE )
    KRATOS_REGISTER_VARIABLE( SLIP_TENDENCY )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LIQUID_FLUX_VECTOR )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LOCAL_LIQUID_FLUX_VECTOR )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( GAS_FLUX_VECTOR )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LOCAL_GAS_FLUX_VECTOR )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( CONTACT_STRESS_VECTOR )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LOCAL_STRESS_VECTOR )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LOCAL_RELATIVE_DISPLACEMENT_VECTOR )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_MATRIX )
    KRATOS_REGISTER_VARIABLE( LOCAL_PERMEABILITY_MATRIX )
    KRATOS_REGISTER_VARIABLE( LIQUID_PERMEABILITY_MATRIX )
    KRATOS_REGISTER_VARIABLE( GAS_PERMEABILITY_MATRIX )

    KRATOS_REGISTER_VARIABLE( CRITICAL_DISPLACEMENT )

    KRATOS_REGISTER_VARIABLE( IS_CONVERGED )

    KRATOS_REGISTER_VARIABLE( TOTAL_STRESS_TENSOR )

    KRATOS_REGISTER_VARIABLE( INITIAL_STRESS_TENSOR )

    KRATOS_REGISTER_VARIABLE( STATE_VARIABLE )
    KRATOS_REGISTER_VARIABLE( ARC_LENGTH_LAMBDA )
    KRATOS_REGISTER_VARIABLE( ARC_LENGTH_RADIUS_FACTOR )

    KRATOS_REGISTER_VARIABLE( TIME_UNIT_CONVERTER )

    KRATOS_REGISTER_VARIABLE( LOCAL_EQUIVALENT_STRAIN )
    KRATOS_REGISTER_VARIABLE( NONLOCAL_EQUIVALENT_STRAIN )

    KRATOS_REGISTER_VARIABLE( JOINT_WIDTH )

    KRATOS_REGISTER_VARIABLE( NODAL_SMOOTHING )
    KRATOS_REGISTER_VARIABLE( NODAL_CAUCHY_STRESS_TENSOR )
    KRATOS_REGISTER_VARIABLE( EFFECTIVE_STRESS_TENSOR )
    KRATOS_REGISTER_VARIABLE( NODAL_EFFECTIVE_STRESS_TENSOR )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LIQUID_PRESSURE_GRADIENT )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( GAS_PRESSURE_GRADIENT )
    KRATOS_REGISTER_VARIABLE( NODAL_JOINT_AREA )
    KRATOS_REGISTER_VARIABLE( NODAL_JOINT_WIDTH )
    KRATOS_REGISTER_VARIABLE( NODAL_JOINT_DAMAGE )
    KRATOS_REGISTER_VARIABLE( NODAL_MID_PLANE_LIQUID_PRESSURE )
    KRATOS_REGISTER_VARIABLE( NODAL_SLIP_TENDENCY )

    KRATOS_REGISTER_VARIABLE( SHEAR_FRACTURE_ENERGY )

    KRATOS_REGISTER_VARIABLE( BIOT_COEFFICIENT )
    KRATOS_REGISTER_VARIABLE( CURVE_FITTING_ETA )

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( DAMPING_FORCE )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_OLD )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_OLDER )
    KRATOS_REGISTER_VARIABLE( FLUX_RESIDUAL )
    KRATOS_REGISTER_VARIABLE( G_COEFFICIENT )
    KRATOS_REGISTER_VARIABLE( THETA_FACTOR )

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( TARGET_REACTION )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AVERAGE_REACTION )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LOADING_VELOCITY )

}

}// namespace Kratos.
