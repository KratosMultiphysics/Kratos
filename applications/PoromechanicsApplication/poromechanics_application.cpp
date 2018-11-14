
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

    mUPwSmallStrainElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType(3)))),
    mUPwSmallStrainElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType(4)))),
    mUPwSmallStrainElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType(4)))),
    mUPwSmallStrainElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType(8)))),

    mUPwSmallStrainInterfaceElement2D4N( 0, Element::GeometryType::Pointer( new QuadrilateralInterface2D4 <Node<3> >( Element::GeometryType::PointsArrayType(4)))),
    mUPwSmallStrainInterfaceElement3D6N( 0, Element::GeometryType::Pointer( new PrismInterface3D6 <Node<3> >( Element::GeometryType::PointsArrayType(6)))),
    mUPwSmallStrainInterfaceElement3D8N( 0, Element::GeometryType::Pointer( new HexahedraInterface3D8 <Node<3> >( Element::GeometryType::PointsArrayType(8)))),

    mUPwSmallStrainLinkInterfaceElement2D4N( 0, Element::GeometryType::Pointer( new QuadrilateralInterface2D4 <Node<3> >( Element::GeometryType::PointsArrayType(4)))),
    mUPwSmallStrainLinkInterfaceElement3D6N( 0, Element::GeometryType::Pointer( new PrismInterface3D6 <Node<3> >( Element::GeometryType::PointsArrayType(6)))),
    mUPwSmallStrainLinkInterfaceElement3D8N( 0, Element::GeometryType::Pointer( new HexahedraInterface3D8 <Node<3> >( Element::GeometryType::PointsArrayType(8)))),

    mUPwSmallStrainFICElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType(3)))),
    mUPwSmallStrainFICElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType(4)))),
    mUPwSmallStrainFICElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType(4)))),
    mUPwSmallStrainFICElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType(8)))),

    mSmallStrainUPwDiffOrderElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType(6)))),
    mSmallStrainUPwDiffOrderElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType(8)))),
    mSmallStrainUPwDiffOrderElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType(9)))),
    mSmallStrainUPwDiffOrderElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType(10)))),
    mSmallStrainUPwDiffOrderElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType(20)))),
    mSmallStrainUPwDiffOrderElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType(27)))),


    mUPwForceCondition2D1N( 0, Condition::GeometryType::Pointer( new Point2D<Node<3> >( Condition::GeometryType::PointsArrayType(1)))),
    mUPwForceCondition3D1N( 0, Condition::GeometryType::Pointer( new Point3D<Node<3> >( Condition::GeometryType::PointsArrayType(1)))),
    mUPwFaceLoadCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType(2)))),
    mUPwFaceLoadCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType(3)))),
    mUPwFaceLoadCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType(4)))),
    mUPwNormalFaceLoadCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType(2)))),
    mUPwNormalFaceLoadCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType(3)))),
    mUPwNormalFaceLoadCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType(4)))),
    mUPwNormalFluxCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType(2)))),
    mUPwNormalFluxCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType(3)))),
    mUPwNormalFluxCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType(4)))),

    mUPwFaceLoadInterfaceCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType(2)))),
    mUPwFaceLoadInterfaceCondition3D4N( 0, Condition::GeometryType::Pointer( new QuadrilateralInterface3D4 <Node<3> >( Condition::GeometryType::PointsArrayType(4)))),
    mUPwNormalFluxInterfaceCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType(2)))),
    mUPwNormalFluxInterfaceCondition3D4N( 0, Condition::GeometryType::Pointer( new QuadrilateralInterface3D4 <Node<3> >( Condition::GeometryType::PointsArrayType(4)))),

    mUPwNormalFluxFICCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType(2)))),
    mUPwNormalFluxFICCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType(3)))),
    mUPwNormalFluxFICCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType(4)))),

    mLineLoadDiffOrderCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3<Node<3> >( Condition::GeometryType::PointsArrayType(3)))),
    mLineNormalLoadDiffOrderCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3<Node<3> >( Condition::GeometryType::PointsArrayType(3)))),
    mLineNormalFluidFluxDiffOrderCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3<Node<3> >( Condition::GeometryType::PointsArrayType(3)))),
    mSurfaceLoadDiffOrderCondition3D6N( 0, Condition::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Condition::GeometryType::PointsArrayType(6)))),
    mSurfaceLoadDiffOrderCondition3D8N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Condition::GeometryType::PointsArrayType(8)))),
    mSurfaceLoadDiffOrderCondition3D9N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Condition::GeometryType::PointsArrayType(9)))),
    mSurfaceNormalLoadDiffOrderCondition3D6N( 0, Condition::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Condition::GeometryType::PointsArrayType(6)))),
    mSurfaceNormalLoadDiffOrderCondition3D8N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Condition::GeometryType::PointsArrayType(8)))),
    mSurfaceNormalLoadDiffOrderCondition3D9N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Condition::GeometryType::PointsArrayType(9)))),
    mSurfaceNormalFluidFluxDiffOrderCondition3D6N( 0, Condition::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Condition::GeometryType::PointsArrayType(6)))),
    mSurfaceNormalFluidFluxDiffOrderCondition3D8N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Condition::GeometryType::PointsArrayType(8)))),
    mSurfaceNormalFluidFluxDiffOrderCondition3D9N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Condition::GeometryType::PointsArrayType(9))))

    {}

void KratosPoromechanicsApplication::Register()
{
    //Calling base class register to register Kratos components
    KratosApplication::Register();
    KRATOS_INFO("") << "Initializing KratosPoromechanicsApplication... " << std::endl;

    //Register Elements
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainElement2D3N", mUPwSmallStrainElement2D3N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainElement2D4N", mUPwSmallStrainElement2D4N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainElement3D4N", mUPwSmallStrainElement3D4N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainElement3D8N", mUPwSmallStrainElement3D8N )

    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainInterfaceElement2D4N", mUPwSmallStrainInterfaceElement2D4N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainInterfaceElement3D6N", mUPwSmallStrainInterfaceElement3D6N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainInterfaceElement3D8N", mUPwSmallStrainInterfaceElement3D8N )

    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainLinkInterfaceElement2D4N", mUPwSmallStrainLinkInterfaceElement2D4N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainLinkInterfaceElement3D6N", mUPwSmallStrainLinkInterfaceElement3D6N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainLinkInterfaceElement3D8N", mUPwSmallStrainLinkInterfaceElement3D8N )

    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainFICElement2D3N", mUPwSmallStrainFICElement2D3N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainFICElement2D4N", mUPwSmallStrainFICElement2D4N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainFICElement3D4N", mUPwSmallStrainFICElement3D4N )
    KRATOS_REGISTER_ELEMENT( "UPwSmallStrainFICElement3D8N", mUPwSmallStrainFICElement3D8N )

    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwDiffOrderElement2D6N", mSmallStrainUPwDiffOrderElement2D6N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwDiffOrderElement2D8N", mSmallStrainUPwDiffOrderElement2D8N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwDiffOrderElement2D9N", mSmallStrainUPwDiffOrderElement2D9N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwDiffOrderElement3D10N", mSmallStrainUPwDiffOrderElement3D10N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwDiffOrderElement3D20N", mSmallStrainUPwDiffOrderElement3D20N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwDiffOrderElement3D27N", mSmallStrainUPwDiffOrderElement3D27N )

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
    KRATOS_REGISTER_CONSTITUTIVE_LAW("BilinearCohesive3DLaw",mBilinearCohesive3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("BilinearCohesive2DLaw",mBilinearCohesive2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("ExponentialCohesive3DLaw",mExponentialCohesive3DLaw);

    KRATOS_REGISTER_CONSTITUTIVE_LAW( "LocalDamageFlowRule", mLocalDamageFlowRule );
    KRATOS_REGISTER_CONSTITUTIVE_LAW( "NonlocalDamageFlowRule", mNonlocalDamageFlowRule );

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
    KRATOS_REGISTER_VARIABLE( DT_PRESSURE_COEFFICIENT )

    KRATOS_REGISTER_VARIABLE( DT_WATER_PRESSURE )
    KRATOS_REGISTER_VARIABLE( NORMAL_FLUID_FLUX )

    KRATOS_REGISTER_VARIABLE( DENSITY_SOLID )
    KRATOS_REGISTER_VARIABLE( BULK_MODULUS_SOLID )
    KRATOS_REGISTER_VARIABLE( BULK_MODULUS_FLUID )
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

    KRATOS_REGISTER_VARIABLE( IS_CONVERGED )

    KRATOS_REGISTER_VARIABLE( TOTAL_STRESS_TENSOR )

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
}

}// namespace Kratos.
