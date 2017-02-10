//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  ilaria$
//   Date:                $Date:  July 2015$
//   Revision:            $Revision: 1.3 $
//
// 



// System includes


// External includes 


// Project includes
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"

#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"

#include "geometries/triangle_3d_3.h"

#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"

#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"

#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"

#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"

#include "geometries/line_2d.h"
#include "geometries/line_2d_2.h"

#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"

#include "geometries/point_2d.h"
#include "geometries/point_3d.h"

//#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/variables.h"
#include "includes/serializer.h"

#include "particle_mechanics_application.h"
namespace Kratos
{
    //Example
// KRATOS_CREATE_VARIABLE(double, AUX_MESH_VAR)
//KRATOS_CREATE_VARIABLE(double, IS_INTERFACE);
//KRATOS_CREATE_VARIABLE(double, NODAL_AREA);
    
    //element
    KRATOS_CREATE_VARIABLE( int, COUNTER )
    KRATOS_CREATE_VARIABLE( int, MP_NUMBER )
    KRATOS_CREATE_VARIABLE( int, MP_BOOL )
    KRATOS_CREATE_VARIABLE( double, WEIGHT )
    KRATOS_CREATE_VARIABLE( double, MP_MASS )
    KRATOS_CREATE_VARIABLE( double, MP_DENSITY )
    KRATOS_CREATE_VARIABLE( double, MP_VOLUME )
    KRATOS_CREATE_VARIABLE( double, MP_KINETIC_ENERGY )
    KRATOS_CREATE_VARIABLE( double, MP_STRAIN_ENERGY )
    KRATOS_CREATE_VARIABLE( double, MP_TOTAL_ENERGY )
    //KRATOS_CREATE_VARIABLE( double, NODAL_MASS )
    
    
    //constitutive law
    KRATOS_CREATE_VARIABLE( ConstitutiveLaw::Pointer, CONSTITUTIVE_LAW_POINTER )
    
    //nodal dofs
    //KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_DISPLACEMENT )
    KRATOS_CREATE_VARIABLE( double, AUX_R)
    KRATOS_CREATE_VARIABLE( double, AUX_T)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( AUX_R_VEL )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( AUX_T_VEL )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( AUX_R_ACC )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( AUX_T_ACC )
    KRATOS_CREATE_VARIABLE( double, NODAL_LUMPED_MASS)

    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( AUX_VELOCITY )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( AUX_ACCELERATION )
    //MP element variable
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( GAUSS_COORD )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MP_DISPLACEMENT )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MP_VELOCITY )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MP_ACCELERATION )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( AUX_MP_VELOCITY )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( AUX_MP_ACCELERATION )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MP_VOLUME_ACCELERATION )
    KRATOS_CREATE_VARIABLE( Vector, MP_CAUCHY_STRESS_VECTOR )
    KRATOS_CREATE_VARIABLE( Vector, MP_ALMANSI_STRAIN_VECTOR )
    KRATOS_CREATE_VARIABLE( Vector, PREVIOUS_MP_CAUCHY_STRESS_VECTOR )
    KRATOS_CREATE_VARIABLE( Vector, PREVIOUS_MP_ALMANSI_STRAIN_VECTOR )
    KRATOS_CREATE_VARIABLE( Matrix, MP_CONSTITUTIVE_MATRIX )
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_AUX)
    //grid node variable
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( NODAL_MOMENTUM)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( NODAL_INERTIA)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( NODAL_INTERNAL_FORCE )

    KratosParticleMechanicsApplication::KratosParticleMechanicsApplication():
        mUpdatedLagrangian2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
        mUpdatedLagrangian3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
        //mUpdatedLagrangianUP2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
        //mUpdatedLagrangianUP3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
        mUpdatedLagrangian2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4) ) ) )
        //mUpdatedLagrangianUP2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) )
        
        //mTotalLagrangian2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
        //mTotalLagrangian3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
    {}
    
    void KratosParticleMechanicsApplication::Register()
    {
    // calling base class register to register Kratos components
        KratosApplication::Register();
        //std::cout << "     KRATOS  _ |   \\  _ |   ||  | | _ |               " << std::endl;
        //std::cout << "              _| \  \\   | |  | (  | _|                " << std::endl;
        //std::cout << "           __|__/ \__\\|\\_|  | _| _| _| MECHANICS     " << std::endl;
        std::cout << " Initializing KratosParticleMechanicsApplication... " << std::endl;
        
        
// 	    KRATOS_REGISTER_VARIABLE( AUX_MESH_VAR )
// 	    KRATOS_REGISTER_VARIABLE(IS_INTERFACE);
// 	    KRATOS_REGISTER_VARIABLE(NODAL_AREA);

        //Registering elements and conditions here
        KRATOS_REGISTER_ELEMENT( "UpdatedLagrangian2D3N", mUpdatedLagrangian2D3N )
        KRATOS_REGISTER_ELEMENT( "UpdatedLagrangian3D4N", mUpdatedLagrangian3D4N )
        //KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianUP2D3N", mUpdatedLagrangianUP2D3N )
        //KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianUP3D4N", mUpdatedLagrangianUP3D4N )
        KRATOS_REGISTER_ELEMENT( "UpdatedLagrangian2D4N", mUpdatedLagrangian2D4N )
        //KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianUP2D4N", mUpdatedLagrangianUP2D4N )
        
        //KRATOS_REGISTER_ELEMENT( "TotalLagrangian2D3N", mTotalLagrangian2D3N )
        //KRATOS_REGISTER_ELEMENT( "TotalLagrangian3D4N", mTotalLagrangian3D4N )
        
        //element
        KRATOS_REGISTER_VARIABLE( COUNTER )
        KRATOS_REGISTER_VARIABLE( MP_NUMBER )
        KRATOS_REGISTER_VARIABLE( MP_BOOL )
        KRATOS_REGISTER_VARIABLE( WEIGHT )
        KRATOS_REGISTER_VARIABLE( MP_MASS )
        KRATOS_REGISTER_VARIABLE( MP_DENSITY )
        KRATOS_REGISTER_VARIABLE( MP_VOLUME )
        KRATOS_REGISTER_VARIABLE( MP_KINETIC_ENERGY )
        KRATOS_REGISTER_VARIABLE( MP_STRAIN_ENERGY )
        KRATOS_REGISTER_VARIABLE( MP_TOTAL_ENERGY )
        
        
              
        //consitutive law
        
        KRATOS_REGISTER_VARIABLE( CONSTITUTIVE_LAW_POINTER )
        
        ////nodal dofs
        //KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_DISPLACEMENT )
	KRATOS_REGISTER_VARIABLE( AUX_R )
        KRATOS_REGISTER_VARIABLE( AUX_T )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUX_R_VEL )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUX_T_VEL )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUX_R_ACC )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUX_T_ACC )
        KRATOS_REGISTER_VARIABLE( NODAL_LUMPED_MASS )


        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUX_VELOCITY )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUX_ACCELERATION )
        //MP element variable
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( GAUSS_COORD )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MP_DISPLACEMENT )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MP_VELOCITY )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MP_ACCELERATION )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUX_MP_VELOCITY )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUX_MP_ACCELERATION )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MP_VOLUME_ACCELERATION )
        KRATOS_REGISTER_VARIABLE( PREVIOUS_MP_CAUCHY_STRESS_VECTOR )
        KRATOS_REGISTER_VARIABLE( PREVIOUS_MP_ALMANSI_STRAIN_VECTOR )
        KRATOS_REGISTER_VARIABLE( MP_CAUCHY_STRESS_VECTOR )
        KRATOS_REGISTER_VARIABLE( MP_ALMANSI_STRAIN_VECTOR )
        KRATOS_REGISTER_VARIABLE( MP_CONSTITUTIVE_MATRIX )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_AUX )
        //grid node variable
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NODAL_MOMENTUM )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NODAL_INERTIA )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NODAL_INTERNAL_FORCE )

		//Hyperelastic ViscoPlastic laws
        Serializer::Register( "HyperElasticViscoplastic3DLaw", mHyperElasticViscoplastic3DLaw );
        Serializer::Register( "HyperElasticViscoplasticPlaneStrain2DLaw", mHyperElasticViscoplasticPlaneStrain2DLaw );
 
    }

}  // namespace Kratos.


