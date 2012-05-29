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

namespace Kratos
{

        KRATOS_CREATE_VARIABLE(double, DEM_DELTA_TIME)

        KRATOS_CREATE_VARIABLE(Vector,     PARTICLE_ROTATE_SPRING_FAILURE_TYPE)

        typedef  vector<array_1d<double,3> >  VectorArray3Double;

        KRATOS_CREATE_VARIABLE( VectorArray3Double, PARTICLE_ROTATE_SPRING_MOMENT )
      

        KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( PARTICLE_MOMENT );
        KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( PARTICLE_ROTATION_ANGLE );
        KRATOS_CREATE_VARIABLE(double,  PARTICLE_MOMENT_OF_INERTIA);

        KRATOS_CREATE_VARIABLE(Vector, INITIAL_AXES_TRACKING)
        KRATOS_CREATE_VARIABLE(int, plot_OPTIONS)

        //M:possible future blocks (no FEM) interaction
        KRATOS_CREATE_VARIABLE(Vector,     PARTICLE_BLOCK_CONTACT_FAILURE_ID)
        KRATOS_CREATE_VARIABLE(Vector,     PARTICLE_BLOCK_CONTACT_FORCE)
        KRATOS_CREATE_VARIABLE(Vector,     PARTICLE_BLOCK_IF_INITIAL_CONTACT)
        KRATOS_CREATE_VARIABLE(WeakPointerVector<Element >,     NEIGHBOUR_PARTICLE_BLOCK_ELEMENTS)

	KratosDEMApplication::KratosDEMApplication():
	mSphericParticle(0, Element::GeometryType::Pointer( new Point3D <Node<3> >( Element::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) )
	{}
        
	void KratosDEMApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosDEMApplication... " << std::endl;
                

                KRATOS_REGISTER_VARIABLE(DEM_DELTA_TIME)

                KRATOS_REGISTER_VARIABLE(PARTICLE_BLOCK_CONTACT_FAILURE_ID)
                KRATOS_REGISTER_VARIABLE(PARTICLE_BLOCK_CONTACT_FORCE)
                KRATOS_REGISTER_VARIABLE(PARTICLE_BLOCK_IF_INITIAL_CONTACT)
                KRATOS_REGISTER_VARIABLE(NEIGHBOUR_PARTICLE_BLOCK_ELEMENTS)

                KRATOS_REGISTER_VARIABLE(PARTICLE_ROTATE_SPRING_FAILURE_TYPE)
                KRATOS_REGISTER_VARIABLE(PARTICLE_ROTATE_SPRING_MOMENT)

                KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(PARTICLE_MOMENT)
                KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(PARTICLE_ROTATION_ANGLE)
                KRATOS_REGISTER_VARIABLE(PARTICLE_MOMENT_OF_INERTIA)

                KRATOS_REGISTER_VARIABLE(INITIAL_AXES_TRACKING)
                KRATOS_REGISTER_VARIABLE(plot_OPTIONS)
               
                KRATOS_REGISTER_ELEMENT("SphericParticle", mSphericParticle)
	
        }

}  // namespace Kratos.


