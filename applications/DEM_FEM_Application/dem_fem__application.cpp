//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.3 $
//
// 



// System includes


// External includes 


// Project includes
#include "dem_fem__application.h"


#include "includes/define.h"
#include "includes/variables.h"

// Cfeng:120416 from structual application below

#include "includes/serializer.h"

#include "geometries/point_2d.h"
#include "geometries/point_3d.h"


namespace Kratos
{


    KRATOS_CREATE_VARIABLE(Vector, PRESTRESS)
    KRATOS_CREATE_VARIABLE( Matrix, GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR)

    KRATOS_CREATE_VARIABLE(double,  PARTICLE_NUMBER_OF_NEIGHBOURS)
 

    KRATOS_CREATE_VARIABLE(double,  DEM_FEM_CONVERGENCE_RATIO)


    KRATOS_CREATE_VARIABLE(Vector,     PARTICLE_BLOCK_CONTACT_FAILURE_ID)
    KRATOS_CREATE_VARIABLE(Vector,     PARTICLE_BLOCK_CONTACT_FORCE)
    KRATOS_CREATE_VARIABLE(Vector,     PARTICLE_BLOCK_IF_INITIAL_CONTACT)
    KRATOS_CREATE_VARIABLE(WeakPointerVector<Condition >,     NEIGHBOUR_RIGID_FACES)

    KRATOS_CREATE_VARIABLE(int, IF_BOUNDARY_ELEMENT)
    KRATOS_CREATE_VARIABLE(Vector, IF_BOUNDARY_FACE)
	
	KRATOS_CREATE_VARIABLE(int, ROTATION_OPTION)


    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( PARTICLE_MOMENT );
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( PARTICLE_ROTATION_ANGLE );
    KRATOS_CREATE_VARIABLE(double,  PARTICLE_MOMENT_OF_INERTIA);

/*
    KratosDEM_FEM_Application::KratosDEM_FEM_Application():
            
            mParticle2D( 0, Element::GeometryType::Pointer( new Point2D<Node<3> >( Element::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) ),
            mParticle3D( 0, Element::GeometryType::Pointer( new Point3D<Node<3> >( Element::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) )
            //mParticle2D( 0, Particle::GeometryType::Pointer( new Point3D<Node<3> >( Particle::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) ),
            //mParticle3D( 0, Particle::GeometryType::Pointer( new Point3D<Node<3> >( Particle::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) )

    {}  //////Cfeng: this blanket is essential.120416
*/

    KratosDEM_FEM_Application::KratosDEM_FEM_Application()
    {
        
    }

    void KratosDEM_FEM_Application::Register()
    {
        // calling base class register to register Kratos components
        KratosApplication::Register();
        std::cout << "Initializing KratosDEM_FEM_Application... " << std::endl;


       KRATOS_REGISTER_VARIABLE(PRESTRESS)
       KRATOS_REGISTER_VARIABLE( CONSTITUTIVE_LAW )
       KRATOS_REGISTER_VARIABLE(GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR)

       KRATOS_REGISTER_VARIABLE(PARTICLE_NUMBER_OF_NEIGHBOURS)


       KRATOS_REGISTER_VARIABLE(DEM_FEM_CONVERGENCE_RATIO)


               
       KRATOS_REGISTER_VARIABLE(IF_BOUNDARY_ELEMENT)
       KRATOS_REGISTER_VARIABLE(IF_BOUNDARY_FACE)
	   
	   KRATOS_REGISTER_VARIABLE(ROTATION_OPTION)


       KRATOS_REGISTER_VARIABLE(PARTICLE_BLOCK_CONTACT_FAILURE_ID)
       KRATOS_REGISTER_VARIABLE(PARTICLE_BLOCK_CONTACT_FORCE)
       KRATOS_REGISTER_VARIABLE(PARTICLE_BLOCK_IF_INITIAL_CONTACT)
       KRATOS_REGISTER_VARIABLE(NEIGHBOUR_RIGID_FACES)



       KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(PARTICLE_MOMENT)
       KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(PARTICLE_ROTATION_ANGLE)
       KRATOS_REGISTER_VARIABLE(PARTICLE_MOMENT_OF_INERTIA)
         
             
    }

}  // namespace Kratos.


