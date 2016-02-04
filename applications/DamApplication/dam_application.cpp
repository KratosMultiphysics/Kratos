//
//   Project Name:        KratosDamApplication $
//   Last Modified by:    $Author:     IPouplana $
//   Date:                $Date:    December 2015 $
//   Revision:            $Revision:         1.0 $
//

// External includes
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"

// Elements 
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"

#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"

// Conditions
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"

#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/variables.h"

#include "dam_application.h"

namespace Kratos
{

//Create Variables //Note that the application variables must not be defined if they already exist in KRATOS
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_POINT_LOAD )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_LINE_LOAD )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_SURFACE_LOAD )
KRATOS_CREATE_VARIABLE( double, IMPOSED_NORMAL_STRESS )
KRATOS_CREATE_VARIABLE( double, IMPOSED_TANGENTIAL_STRESS )
KRATOS_CREATE_VARIABLE( double, IMPOSED_TEMPERATURE )

//Bofang and Hidrostatic variables for evolution changes
KRATOS_CREATE_VARIABLE( std::string, GRAVITY_DIRECTION )
KRATOS_CREATE_VARIABLE( double, COORDINATE_BASE_DAM )
KRATOS_CREATE_VARIABLE( double, SURFACE_TEMP )
KRATOS_CREATE_VARIABLE( double, BOTTOM_TEMP )
KRATOS_CREATE_VARIABLE( double, HEIGHT_DAM )
KRATOS_CREATE_VARIABLE( double, AMPLITUDE )
KRATOS_CREATE_VARIABLE( double, FREQUENCY )
KRATOS_CREATE_VARIABLE( double, DAY_MAXIMUM )
KRATOS_CREATE_VARIABLE( double, SPECIFIC_WEIGHT )

// Thermal Variables
KRATOS_CREATE_VARIABLE( Matrix, THERMAL_STRESS_TENSOR )
KRATOS_CREATE_VARIABLE( Matrix, MECHANICAL_STRESS_TENSOR )
KRATOS_CREATE_VARIABLE( Matrix, THERMAL_STRAIN_TENSOR )

KRATOS_CREATE_VARIABLE( Vector, THERMAL_STRESS_VECTOR )
KRATOS_CREATE_VARIABLE( Vector, MECHANICAL_STRESS_VECTOR )
KRATOS_CREATE_VARIABLE( Vector, THERMAL_STRAIN_VECTOR )


KratosDamApplication::KratosDamApplication():

    mSmallDisplacementThermoMechanicElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mSmallDisplacementThermoMechanicElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mSmallDisplacementThermoMechanicElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mSmallDisplacementThermoMechanicElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mSmallDisplacementThermoMechanicElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) ),
    
    mSmallDisplacementThermoMechanicElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mSmallDisplacementThermoMechanicElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10, Node<3>() ) ) ) ),
    mSmallDisplacementThermoMechanicElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mSmallDisplacementThermoMechanicElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20, Node<3>() ) ) ) ),
    mSmallDisplacementThermoMechanicElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) ),

    mPointLoadCondition2D( 0, Condition::GeometryType::Pointer( new Point2D<Node<3> >( Condition::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) ),
    mPointLoadCondition3D( 0, Condition::GeometryType::Pointer( new Point3D<Node<3> >( Condition::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) ),

    mLineLoadCondition2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mLineLoadCondition3N( 0, Condition::GeometryType::Pointer( new Line2D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    
    mLineNormalLoadCondition2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mLineNormalLoadCondition3N( 0, Condition::GeometryType::Pointer( new Line2D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    
    mSurfaceLoadCondition3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mSurfaceLoadCondition4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mSurfaceLoadCondition6N( 0, Condition::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Condition::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mSurfaceLoadCondition8N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Condition::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mSurfaceLoadCondition9N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Condition::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) ),
    
    mSurfaceNormalLoadCondition3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mSurfaceNormalLoadCondition4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mSurfaceNormalLoadCondition6N( 0, Condition::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Condition::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mSurfaceNormalLoadCondition8N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Condition::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mSurfaceNormalLoadCondition9N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Condition::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )

{}

void KratosDamApplication::Register()
{
    //Calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosDamApplication... " << std::endl;

    //Register Elements
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement2D3N", mSmallDisplacementThermoMechanicElement2D3N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement2D4N", mSmallDisplacementThermoMechanicElement2D4N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement3D4N", mSmallDisplacementThermoMechanicElement3D4N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement3D8N", mSmallDisplacementThermoMechanicElement3D8N )

    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement2D6N", mSmallDisplacementThermoMechanicElement2D6N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement2D8N", mSmallDisplacementThermoMechanicElement2D8N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement2D9N", mSmallDisplacementThermoMechanicElement2D9N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement3D10N", mSmallDisplacementThermoMechanicElement3D10N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement3D20N", mSmallDisplacementThermoMechanicElement3D20N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement3D27N", mSmallDisplacementThermoMechanicElement3D27N )

    //Register Conditions
    KRATOS_REGISTER_CONDITION( "PointLoadCondition2D", mPointLoadCondition2D )
    KRATOS_REGISTER_CONDITION( "PointLoadCondition3D", mPointLoadCondition3D )

    KRATOS_REGISTER_CONDITION( "LineLoadCondition2N", mLineLoadCondition2N )
    KRATOS_REGISTER_CONDITION( "LineLoadCondition3N", mLineLoadCondition3N )
    
    KRATOS_REGISTER_CONDITION( "LineNormalLoadCondition2N", mLineNormalLoadCondition2N )
    KRATOS_REGISTER_CONDITION( "LineNormalLoadCondition3N", mLineNormalLoadCondition3N )
    
    KRATOS_REGISTER_CONDITION( "SurfaceLoadCondition3N", mSurfaceLoadCondition3N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadCondition4N", mSurfaceLoadCondition4N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadCondition6N", mSurfaceLoadCondition6N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadCondition8N", mSurfaceLoadCondition8N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadCondition9N", mSurfaceLoadCondition9N )
    
    KRATOS_REGISTER_CONDITION( "SurfaceNormalLoadCondition3N", mSurfaceNormalLoadCondition3N )
    KRATOS_REGISTER_CONDITION( "SurfaceNormalLoadCondition4N", mSurfaceNormalLoadCondition4N )
    KRATOS_REGISTER_CONDITION( "SurfaceNormalLoadCondition6N", mSurfaceNormalLoadCondition6N )
    KRATOS_REGISTER_CONDITION( "SurfaceNormalLoadCondition8N", mSurfaceNormalLoadCondition8N )
    KRATOS_REGISTER_CONDITION( "SurfaceNormalLoadCondition9N", mSurfaceNormalLoadCondition9N )

    //Register Constitutive Laws
    Serializer::Register("ThermalLinearElastic3DLaw",mThermalLinearElastic3DLaw);
    Serializer::Register("ThermalLinearElastic2DPlaneStress",mThermalLinearElastic2DPlaneStress);
    Serializer::Register("ThermalLinearElastic2DPlaneStrain",mThermalLinearElastic2DPlaneStrain);

    //Register Variables
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_POINT_LOAD )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_LINE_LOAD )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_SURFACE_LOAD )
    KRATOS_REGISTER_VARIABLE( IMPOSED_NORMAL_STRESS )
    KRATOS_REGISTER_VARIABLE( IMPOSED_TANGENTIAL_STRESS )
    KRATOS_REGISTER_VARIABLE( IMPOSED_TEMPERATURE )
    
    //Bofang and Hidrostatic variables for evolution changes
    KRATOS_REGISTER_VARIABLE( GRAVITY_DIRECTION )
    KRATOS_REGISTER_VARIABLE( COORDINATE_BASE_DAM )
    KRATOS_REGISTER_VARIABLE( SURFACE_TEMP )
    KRATOS_REGISTER_VARIABLE( BOTTOM_TEMP )
    KRATOS_REGISTER_VARIABLE( HEIGHT_DAM )
    KRATOS_REGISTER_VARIABLE( AMPLITUDE )
    KRATOS_REGISTER_VARIABLE( FREQUENCY )
    KRATOS_REGISTER_VARIABLE( DAY_MAXIMUM )
    KRATOS_REGISTER_VARIABLE( SPECIFIC_WEIGHT )
    
    // Thermal Variables
    KRATOS_REGISTER_VARIABLE( THERMAL_STRESS_TENSOR )
    KRATOS_REGISTER_VARIABLE( MECHANICAL_STRESS_TENSOR )
    KRATOS_REGISTER_VARIABLE( THERMAL_STRAIN_TENSOR )
    
    
    KRATOS_REGISTER_VARIABLE( THERMAL_STRESS_VECTOR )
    KRATOS_REGISTER_VARIABLE( MECHANICAL_STRESS_VECTOR )
    KRATOS_REGISTER_VARIABLE( THERMAL_STRAIN_VECTOR )
}

}// namespace Kratos.
