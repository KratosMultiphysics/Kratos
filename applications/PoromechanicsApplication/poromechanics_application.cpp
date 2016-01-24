//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// External includes
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"

#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/variables.h"

#include "poromechanics_application.h"

namespace Kratos
{

//Create Variables 

//Warning: Note that the application variables must not be defined if they already exist in "includes/variables.h" or in "includes/cfd_variables.h"

KRATOS_CREATE_VARIABLE( double, BETA_NEWMARK )
KRATOS_CREATE_VARIABLE( double, GAMMA_NEWMARK )
KRATOS_CREATE_VARIABLE( double, THETA_NEWMARK )

KRATOS_CREATE_VARIABLE( ConstitutiveLaw::Pointer, CONSTITUTIVE_LAW_POINTER )
KRATOS_CREATE_VARIABLE( std::string, CONSTITUTIVE_LAW_NAME )

KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_DISPLACEMENT )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( POINT_LOAD )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_POINT_LOAD )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( LINE_LOAD ) 
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_LINE_LOAD )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( SURFACE_LOAD )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_SURFACE_LOAD )
KRATOS_CREATE_VARIABLE( double, IMPOSED_NORMAL_STRESS )
KRATOS_CREATE_VARIABLE( double, IMPOSED_TANGENTIAL_STRESS )

KRATOS_CREATE_VARIABLE( double, DERIVATIVE_WATER_PRESSURE )
KRATOS_CREATE_VARIABLE( double, IMPOSED_FLUID_PRESSURE )
KRATOS_CREATE_VARIABLE( double, NORMAL_FLUID_FLUX )
KRATOS_CREATE_VARIABLE( double, IMPOSED_NORMAL_FLUID_FLUX )

KRATOS_CREATE_VARIABLE( double, DENSITY_SOLID )
KRATOS_CREATE_VARIABLE( double, BULK_MODULUS_SOLID )
KRATOS_CREATE_VARIABLE( double, BULK_MODULUS_FLUID )
KRATOS_CREATE_VARIABLE( double, PERMEABILITY_XX )
KRATOS_CREATE_VARIABLE( double, PERMEABILITY_YY )
KRATOS_CREATE_VARIABLE( double, PERMEABILITY_ZZ )
KRATOS_CREATE_VARIABLE( double, PERMEABILITY_XY )
KRATOS_CREATE_VARIABLE( double, PERMEABILITY_YZ )
KRATOS_CREATE_VARIABLE( double, PERMEABILITY_ZX )

KRATOS_CREATE_VARIABLE( Vector, CAUCHY_STRESS_VECTOR )
KRATOS_CREATE_VARIABLE( Vector, GREEN_LAGRANGE_STRAIN_VECTOR ) 
KRATOS_CREATE_VARIABLE( double, VON_MISES_STRESS )

KRATOS_CREATE_VARIABLE( Vector, FLUID_FLUX )


KratosPoromechanicsApplication::KratosPoromechanicsApplication():

    mSmallStrainUPwElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mSmallStrainUPwElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mSmallStrainUPwElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mSmallStrainUPwElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),

    mSmallStrainUPwDiffOrderElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mSmallStrainUPwDiffOrderElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mSmallStrainUPwDiffOrderElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) ),
    mSmallStrainUPwDiffOrderElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10, Node<3>() ) ) ) ),
    mSmallStrainUPwDiffOrderElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20, Node<3>() ) ) ) ),
    mSmallStrainUPwDiffOrderElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) ),

    mSmallStrainUPwFICElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mSmallStrainUPwFICElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mSmallStrainUPwFICElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mSmallStrainUPwFICElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),

    mPointLoadCondition2D( 0, Condition::GeometryType::Pointer( new Point2D<Node<3> >( Condition::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) ),
    mPointLoadCondition3D( 0, Condition::GeometryType::Pointer( new Point3D<Node<3> >( Condition::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) ),
    
    mLineLoadCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mLineNormalLoadCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mLineNormalFluidFluxCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mSurfaceLoadCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mSurfaceLoadCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mSurfaceNormalLoadCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mSurfaceNormalLoadCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mSurfaceNormalFluidFluxCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mSurfaceNormalFluidFluxCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),

    mLineLoadDiffOrderCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mLineNormalLoadDiffOrderCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mLineNormalFluidFluxDiffOrderCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3<Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mSurfaceLoadDiffOrderCondition3D6N( 0, Condition::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Condition::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mSurfaceLoadDiffOrderCondition3D8N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Condition::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mSurfaceLoadDiffOrderCondition3D9N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Condition::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) ),
    mSurfaceNormalLoadDiffOrderCondition3D6N( 0, Condition::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Condition::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mSurfaceNormalLoadDiffOrderCondition3D8N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Condition::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mSurfaceNormalLoadDiffOrderCondition3D9N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Condition::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) ),
    mSurfaceNormalFluidFluxDiffOrderCondition3D6N( 0, Condition::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Condition::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mSurfaceNormalFluidFluxDiffOrderCondition3D8N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Condition::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mSurfaceNormalFluidFluxDiffOrderCondition3D9N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Condition::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) ),
    
    mLineNormalFluidFluxFICCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mSurfaceNormalFluidFluxFICCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mSurfaceNormalFluidFluxFICCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )

{}
 	
void KratosPoromechanicsApplication::Register()
{
    //Calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosPoromechanicsApplication... " << std::endl;
 
    //Register Elements
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwElement2D3N", mSmallStrainUPwElement2D3N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwElement2D4N", mSmallStrainUPwElement2D4N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwElement3D4N", mSmallStrainUPwElement3D4N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwElement3D8N", mSmallStrainUPwElement3D8N )

    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwDiffOrderElement2D6N", mSmallStrainUPwDiffOrderElement2D6N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwDiffOrderElement2D8N", mSmallStrainUPwDiffOrderElement2D8N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwDiffOrderElement2D9N", mSmallStrainUPwDiffOrderElement2D9N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwDiffOrderElement3D10N", mSmallStrainUPwDiffOrderElement3D10N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwDiffOrderElement3D20N", mSmallStrainUPwDiffOrderElement3D20N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwDiffOrderElement3D27N", mSmallStrainUPwDiffOrderElement3D27N )

    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwFICElement2D3N", mSmallStrainUPwFICElement2D3N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwFICElement2D4N", mSmallStrainUPwFICElement2D4N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwFICElement3D4N", mSmallStrainUPwFICElement3D4N )
    KRATOS_REGISTER_ELEMENT( "SmallStrainUPwFICElement3D8N", mSmallStrainUPwFICElement3D8N )

    //Register Conditions
    KRATOS_REGISTER_CONDITION( "PointLoadCondition2D", mPointLoadCondition2D )
    KRATOS_REGISTER_CONDITION( "PointLoadCondition3D", mPointLoadCondition3D )
    
    KRATOS_REGISTER_CONDITION( "LineLoadCondition2D2N", mLineLoadCondition2D2N )
    KRATOS_REGISTER_CONDITION( "LineNormalLoadCondition2D2N", mLineNormalLoadCondition2D2N )
    KRATOS_REGISTER_CONDITION( "LineNormalFluidFluxCondition2D2N", mLineNormalFluidFluxCondition2D2N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadCondition3D3N", mSurfaceLoadCondition3D3N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadCondition3D4N", mSurfaceLoadCondition3D4N )
    KRATOS_REGISTER_CONDITION( "SurfaceNormalLoadCondition3D3N", mSurfaceNormalLoadCondition3D3N )
    KRATOS_REGISTER_CONDITION( "SurfaceNormalLoadCondition3D4N", mSurfaceNormalLoadCondition3D4N )
    KRATOS_REGISTER_CONDITION( "SurfaceNormalFluidFluxCondition3D3N", mSurfaceNormalFluidFluxCondition3D3N )
    KRATOS_REGISTER_CONDITION( "SurfaceNormalFluidFluxCondition3D4N", mSurfaceNormalFluidFluxCondition3D4N )

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

    KRATOS_REGISTER_CONDITION( "LineNormalFluidFluxFICCondition2D2N", mLineNormalFluidFluxFICCondition2D2N )
    KRATOS_REGISTER_CONDITION( "SurfaceNormalFluidFluxFICCondition3D3N", mSurfaceNormalFluidFluxFICCondition3D3N )
    KRATOS_REGISTER_CONDITION( "SurfaceNormalFluidFluxFICCondition3D4N", mSurfaceNormalFluidFluxFICCondition3D4N )
    
    //Register Constitutive Laws
    Serializer::Register("LinearElastic3DLaw",mLinearElastic3DLaw);
    Serializer::Register("LinearElastic2DPlaneStressLaw",mLinearElastic2DPlaneStressLaw);
    Serializer::Register("LinearElastic2DPlaneStrainLaw",mLinearElastic2DPlaneStrainLaw);


    //Register Variables
    KRATOS_REGISTER_VARIABLE( BETA_NEWMARK )
    KRATOS_REGISTER_VARIABLE( GAMMA_NEWMARK )
    KRATOS_REGISTER_VARIABLE( THETA_NEWMARK )
    
    KRATOS_REGISTER_VARIABLE( CONSTITUTIVE_LAW_POINTER )
    KRATOS_REGISTER_VARIABLE( CONSTITUTIVE_LAW_NAME )

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_DISPLACEMENT )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( POINT_LOAD )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_POINT_LOAD )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LINE_LOAD )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_LINE_LOAD )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( SURFACE_LOAD )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_SURFACE_LOAD )
    KRATOS_REGISTER_VARIABLE( IMPOSED_NORMAL_STRESS )
    KRATOS_REGISTER_VARIABLE( IMPOSED_TANGENTIAL_STRESS )

    KRATOS_REGISTER_VARIABLE( DERIVATIVE_WATER_PRESSURE )
    KRATOS_REGISTER_VARIABLE( IMPOSED_FLUID_PRESSURE )
    KRATOS_REGISTER_VARIABLE( NORMAL_FLUID_FLUX )
    KRATOS_REGISTER_VARIABLE( IMPOSED_NORMAL_FLUID_FLUX )

    KRATOS_REGISTER_VARIABLE( DENSITY_SOLID )
    KRATOS_REGISTER_VARIABLE( BULK_MODULUS_SOLID )
    KRATOS_REGISTER_VARIABLE( BULK_MODULUS_FLUID )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_XX )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_YY )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_ZZ )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_XY )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_YZ )
    KRATOS_REGISTER_VARIABLE( PERMEABILITY_ZX )

    KRATOS_REGISTER_VARIABLE( CAUCHY_STRESS_VECTOR )
    KRATOS_REGISTER_VARIABLE( GREEN_LAGRANGE_STRAIN_VECTOR )
    KRATOS_REGISTER_VARIABLE( VON_MISES_STRESS )

    KRATOS_REGISTER_VARIABLE( FLUID_FLUX )
}

}// namespace Kratos.
