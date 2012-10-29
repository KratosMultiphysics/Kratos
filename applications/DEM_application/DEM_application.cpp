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
#include "includes/define.h"


#include "includes/variables.h"
#include "DEM_application.h"
#include "geometries/point_3d.h"
#include "geometries/point_2d.h"
#include "geometries/line_3d_2.h"

namespace Kratos
{

    /*Define In Global variables.cpp
        KRATOS_CREATE_VARIABLE(double, DELTA_TIME)
        KRATOS_CREATE_VARIABLE(Vector,     PARTICLE_ROTATE_SPRING_FAILURE_TYPE)
        typedef  vector<array_1d<double,3> >  VectorArray3Double;
        KRATOS_CREATE_VARIABLE( VectorArray3Double, PARTICLE_ROTATE_SPRING_MOMENT )
   */
        KRATOS_CREATE_VARIABLE( int, SKIN_SPHERE )
        KRATOS_CREATE_VARIABLE( double, EXPORT_SKIN_SPHERE )
                
        KRATOS_CREATE_VARIABLE(double, EXPORT_ID)

        KRATOS_CREATE_VARIABLE (int, NEIGH_INITIALIZED)
        KRATOS_CREATE_VARIABLE( int, TOTAL_CONTACTS )

        //KRATOS_CREATE_VARIABLE (vector<int> EXPORT_NEIGHBOUR_LIST )

        KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( DAMP_FORCES )
        KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( TOTAL_FORCES )

        KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( PARTICLE_MOMENT );

        KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( PARTICLE_ROTATION_ANGLE );
        KRATOS_CREATE_VARIABLE(double,  PARTICLE_MOMENT_OF_INERTIA);

        KRATOS_CREATE_VARIABLE(Vector, INITIAL_AXES_TRACKING)
        KRATOS_CREATE_VARIABLE(int, plot_OPTIONS)

        KRATOS_CREATE_VARIABLE(int, VIRTUAL_MASS_OPTION)
        KRATOS_CREATE_VARIABLE(double, NODAL_MASS_COEFF)
        KRATOS_CREATE_VARIABLE(double, HISTORICAL_MIN_K)
        KRATOS_CREATE_VARIABLE(int, CRITICAL_TIME_OPTION)

        //M:possible future blocks (no FEM) interaction
        KRATOS_CREATE_VARIABLE(Vector,     PARTICLE_BLOCK_CONTACT_FAILURE_ID)
        KRATOS_CREATE_VARIABLE(Vector,     PARTICLE_BLOCK_CONTACT_FORCE)
        KRATOS_CREATE_VARIABLE(Vector,     PARTICLE_BLOCK_IF_INITIAL_CONTACT)
        KRATOS_CREATE_VARIABLE(WeakPointerVector<Element >,     NEIGHBOUR_PARTICLE_BLOCK_ELEMENTS)

        KRATOS_CREATE_VARIABLE(int, IF_BOUNDARY_ELEMENT)
        KRATOS_CREATE_VARIABLE(Vector, IF_BOUNDARY_FACE)


	KratosDEMApplication::KratosDEMApplication():
	mSphericParticle2D( 0, Element::GeometryType::Pointer( new Point2D<Node<3> >( Element::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) ),
        mSphericParticle3D( 0, Element::GeometryType::Pointer( new Point3D<Node<3> >( Element::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) ),
        mDEM_FEM_Particle2D( 0, Element::GeometryType::Pointer( new Point2D<Node<3> >( Element::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) ),
        mDEM_FEM_Particle3D( 0, Element::GeometryType::Pointer( new Point3D<Node<3> >( Element::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) ),
        mParticleContactElement( 0, Element::GeometryType::Pointer( new Line3D2<Node<3> >( Element::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) )
                             
	{}
        
	void KratosDEMApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosDEMApplication... " << std::endl;

                KRATOS_REGISTER_VARIABLE( SKIN_SPHERE )
                KRATOS_REGISTER_VARIABLE( EXPORT_SKIN_SPHERE )
                        
                KRATOS_REGISTER_VARIABLE( EXPORT_ID)

                KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( DAMP_FORCES )
                KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( TOTAL_FORCES )

              //  KRATOS_REGISTER_VARIABLE ( EXPORT_NEIGHBOUR_LIST )

                KRATOS_REGISTER_VARIABLE ( NEIGH_INITIALIZED )
                KRATOS_REGISTER_VARIABLE ( TOTAL_CONTACTS )

                KRATOS_REGISTER_VARIABLE(VIRTUAL_MASS_OPTION)
                KRATOS_REGISTER_VARIABLE(NODAL_MASS_COEFF)
                
                KRATOS_REGISTER_VARIABLE(HISTORICAL_MIN_K)
                KRATOS_REGISTER_VARIABLE(CRITICAL_TIME_OPTION)


                /* Define In Global variables.cpp
                KRATOS_REGISTER_VARIABLE(DELTA_TIME)
                KRATOS_REGISTER_VARIABLE(PARTICLE_ROTATE_SPRING_FAILURE_TYPE)
                KRATOS_REGISTER_VARIABLE(PARTICLE_ROTATE_SPRING_MOMENT)
                 */

                KRATOS_REGISTER_VARIABLE(PARTICLE_BLOCK_CONTACT_FAILURE_ID)
                KRATOS_REGISTER_VARIABLE(PARTICLE_BLOCK_CONTACT_FORCE)
                KRATOS_REGISTER_VARIABLE(PARTICLE_BLOCK_IF_INITIAL_CONTACT)
                KRATOS_REGISTER_VARIABLE(NEIGHBOUR_PARTICLE_BLOCK_ELEMENTS)

 

                KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(PARTICLE_MOMENT)
                KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(PARTICLE_ROTATION_ANGLE)
                KRATOS_REGISTER_VARIABLE(PARTICLE_MOMENT_OF_INERTIA)

                KRATOS_REGISTER_VARIABLE(INITIAL_AXES_TRACKING)
                KRATOS_REGISTER_VARIABLE(plot_OPTIONS)

                KRATOS_REGISTER_VARIABLE(IF_BOUNDARY_ELEMENT)
                KRATOS_REGISTER_VARIABLE(IF_BOUNDARY_FACE)

               
                KRATOS_REGISTER_ELEMENT("SphericParticle2D", mSphericParticle2D)
                KRATOS_REGISTER_ELEMENT("SphericParticle3D", mSphericParticle3D)

                KRATOS_REGISTER_ELEMENT("DEM_FEM_Particle2D", mDEM_FEM_Particle2D)
                KRATOS_REGISTER_ELEMENT("DEM_FEM_Particle3D", mDEM_FEM_Particle3D)
                        
                KRATOS_REGISTER_ELEMENT("ParticleContactElement", mParticleContactElement)
                

	
        }

}  // namespace Kratos.


