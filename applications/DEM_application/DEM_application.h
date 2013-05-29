//
//   Project Name:        Kratos
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_DEM_APPLICATION_H_INCLUDED )
#define  KRATOS_DEM_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"


#include "includes/variables.h"
#include "custom_elements/spheric_particle.h"
#include "custom_elements/spheric_continuum_particle.h"
#include "custom_elements/spheric_swimming_particle.h"
#include "custom_elements/DEM_FEM_Particle.h"
#include "custom_elements/Particle_Contact_Element.h"


const long double pi = 3.141592653589793238462643383279;  //this one is used in the forward_euler_scheme... does it have to be here?

namespace Kratos
{
  
   
  
  
    
    KRATOS_DEFINE_VARIABLE( int, BOUNDING_BOX_OPTION )
    KRATOS_DEFINE_VARIABLE( bool, DELTA_OPTION )
    KRATOS_DEFINE_VARIABLE( int, ROTATION_OPTION )
    KRATOS_DEFINE_VARIABLE( int, CRITICAL_TIME_OPTION )
    KRATOS_DEFINE_VARIABLE( int, FORCE_CALCULATION_TYPE )
    KRATOS_DEFINE_VARIABLE( int, DAMP_TYPE )
    KRATOS_DEFINE_VARIABLE( int, ROTA_DAMP_TYPE )
    KRATOS_DEFINE_VARIABLE( int, VIRTUAL_MASS_OPTION ) 
    KRATOS_DEFINE_VARIABLE( int, UNIFORM_MATERIAL_OPTION )
    KRATOS_DEFINE_VARIABLE( int, ACTIVATE_SEARCH )
    KRATOS_DEFINE_VARIABLE( int, LIMIT_SURFACE_OPTION )
    KRATOS_DEFINE_VARIABLE( int, NEIGH_INITIALIZED )
    
    //... to be reviewed
    KRATOS_DEFINE_VARIABLE( int, SURFACE_FRICC )
    KRATOS_DEFINE_VARIABLE( int, IF_BOUNDARY_ELEMENT )
    KRATOS_DEFINE_VARIABLE( Vector, IF_BOUNDARY_FACE )  
    KRATOS_DEFINE_VARIABLE( int, plot_OPTIONS )
    
    // CONTINUUM
 
    KRATOS_DEFINE_VARIABLE( bool, CONTINUUM_OPTION )
    KRATOS_DEFINE_VARIABLE( int, CASE_OPTION )
    KRATOS_DEFINE_VARIABLE( int, SKIN_SPHERE )
    KRATOS_DEFINE_VARIABLE( int, ROTATION_SPRING_OPTION ) 
    KRATOS_DEFINE_VARIABLE( int, NON_LINEAR_OPTION )
    KRATOS_DEFINE_VARIABLE( int, CONTACT_MESH_OPTION )
    KRATOS_DEFINE_VARIABLE( int, FAILURE_CRITERION_OPTION )
    KRATOS_DEFINE_VARIABLE( int, CONCRETE_TEST_OPTION )
  
  
//classify this ones:

KRATOS_DEFINE_VARIABLE(std::size_t,  POTENTIAL_NEIGHBOURS)

KRATOS_DEFINE_VARIABLE(double,  PARTICLE_MASS)
KRATOS_DEFINE_VARIABLE(double,  PARTICLE_INERTIA)
KRATOS_DEFINE_VARIABLE(double,  PARTICLE_DENSITY)
KRATOS_DEFINE_VARIABLE(double,  PARTICLE_STATIC_FRICTION_COEF)
KRATOS_DEFINE_VARIABLE(double,  PARTICLE_DYNAMIC_FRICTION_COEF)

KRATOS_DEFINE_VARIABLE(double,  RESTITUTION_COEFF)        

KRATOS_DEFINE_VARIABLE(double, PARTICLE_COHESION)
KRATOS_DEFINE_VARIABLE(int, PARTICLE_CONTINUUM)
KRATOS_DEFINE_VARIABLE(double,  PARTICLE_TENSION)
KRATOS_DEFINE_VARIABLE(double,  PARTICLE_FRICTION)
KRATOS_DEFINE_VARIABLE(double,  PARTICLE_ROTATION_DAMP_RATIO)
KRATOS_DEFINE_VARIABLE(int,  PARTICLE_FAILURE_ID)
KRATOS_DEFINE_VARIABLE(double,  EXPORT_PARTICLE_FAILURE_ID)




KRATOS_DEFINE_VARIABLE( int, GROUP_ID )

KRATOS_DEFINE_VARIABLE( double, GLOBAL_VARIABLES_OPTION )
KRATOS_DEFINE_VARIABLE( double, GLOBAL_KN )
KRATOS_DEFINE_VARIABLE( double, GLOBAL_KT )
KRATOS_DEFINE_VARIABLE( double, GLOBAL_KR )
KRATOS_DEFINE_VARIABLE( double, GLOBAL_RN )
KRATOS_DEFINE_VARIABLE( double, GLOBAL_RT )
KRATOS_DEFINE_VARIABLE( double, GLOBAL_RR )
KRATOS_DEFINE_VARIABLE( double, GLOBAL_FRI_ANG )

KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DELTA_DISPLACEMENT)
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DELTA_ROTA_DISPLACEMENT)        

        
//typedef boost::numeric::ublas::vector<array_1d<double,3> > VectorArray3Double;
typedef vector<array_1d<double,3> > VectorArray3Double;

KRATOS_DEFINE_VARIABLE( int, TRIHEDRON_OPTION )
KRATOS_DEFINE_VARIABLE( VectorArray3Double, ARROW_POINT )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( EULER_ANGLES )

//KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( VECTOR_ARROW_X )
//KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( VECTOR_ARROW_Y )
//KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( VECTOR_ARROW_Z )

KRATOS_DEFINE_VARIABLE( Vector, PARTICLE_CONTACT_DELTA )
KRATOS_DEFINE_VARIABLE( Vector, PARTICLE_INITIAL_DELTA )
KRATOS_DEFINE_VARIABLE( vector<int>, PARTICLE_CONTACT_FAILURE_ID )
//KRATOS_DEFINE_VARIABLE( Vector, PARTICLE_CONTACT_FAILURE_ID )
KRATOS_DEFINE_VARIABLE( VectorArray3Double, PARTICLE_CONTACT_FORCES )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(APPLIED_FORCE)

//dummy variables:   
KRATOS_DEFINE_VARIABLE(double,  DUMMY_FORCES)
KRATOS_DEFINE_VARIABLE(int,  DUMMY_SWITCH)
KRATOS_DEFINE_VARIABLE(double,  DUMMY_LOCAL_AXES)

  KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( PARTICLE_MOMENT );
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( FLUID_VEL_PROJECTED );



KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( DRAG_FORCE );
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( DRAG_REACTION );
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( PRESSURE_GRADIENT );
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( PRESSURE_GRAD_PROJECTED );
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( PARTICLE_ROTATION_ANGLE );
KRATOS_DEFINE_VARIABLE(double,  PARTICLE_MOMENT_OF_INERTIA);

KRATOS_DEFINE_VARIABLE(Vector,     INITIAL_AXES_TRACKING)
KRATOS_DEFINE_VARIABLE(int,     plot_OPTIONS)


  
 

  
  
  
  
  
  
  
  
  
  
  
  
  
        KRATOS_DEFINE_VARIABLE( int, NON_LINEAR_OPTION)    
        KRATOS_DEFINE_VARIABLE(double, SLOPE_FRACTION_N1)
        KRATOS_DEFINE_VARIABLE(double, SLOPE_FRACTION_N2)
		KRATOS_DEFINE_VARIABLE(double, NON_ELASTIC_STAGE)  
        KRATOS_DEFINE_VARIABLE(double, SLOPE_LIMIT_COEFF_C1)
        KRATOS_DEFINE_VARIABLE(double, SLOPE_LIMIT_COEFF_C2)
    
        KRATOS_DEFINE_VARIABLE (double, LOW_POISSON_FORCE)        
        KRATOS_DEFINE_VARIABLE (double, HIGH_POISSON_FORCE)
    
        KRATOS_DEFINE_VARIABLE (double, FIXED_VEL_BOT )
        KRATOS_DEFINE_VARIABLE (double, FIXED_VEL_TOP )

        KRATOS_DEFINE_VARIABLE(double, DEM_MAGIC_FACTOR )
        KRATOS_DEFINE_VARIABLE(int, NUM_OF_NEIGH )

        KRATOS_DEFINE_VARIABLE(double, PREDEFINED_SKIN)

        KRATOS_DEFINE_VARIABLE(double, AREA_VERTICAL_TAPA)
        KRATOS_DEFINE_VARIABLE(double, AREA_VERTICAL_CENTRE)

        KRATOS_DEFINE_VARIABLE(double, RADIAL_DISPLACEMENT)
      
		KRATOS_DEFINE_VARIABLE(double, DEM_STRESS_XX)
		KRATOS_DEFINE_VARIABLE(double, DEM_STRESS_XY)
		KRATOS_DEFINE_VARIABLE(double, DEM_STRESS_XZ)
		KRATOS_DEFINE_VARIABLE(double, DEM_STRESS_YX)
		KRATOS_DEFINE_VARIABLE(double, DEM_STRESS_YY)
		KRATOS_DEFINE_VARIABLE(double, DEM_STRESS_YZ)
		KRATOS_DEFINE_VARIABLE(double, DEM_STRESS_ZX)
		KRATOS_DEFINE_VARIABLE(double, DEM_STRESS_ZY)
		KRATOS_DEFINE_VARIABLE(double, DEM_STRESS_ZZ)	
		
		KRATOS_DEFINE_VARIABLE(double, LOCAL_CONTACT_AREA_HIGH)
		KRATOS_DEFINE_VARIABLE(double, LOCAL_CONTACT_AREA_LOW)
		KRATOS_DEFINE_VARIABLE(double, MEAN_CONTACT_AREA)
		
	
		KRATOS_DEFINE_VARIABLE(double, FINAL_SIMULATION_TIME)
		KRATOS_DEFINE_VARIABLE(double, INITIAL_PRESSURE_TIME) 
		KRATOS_DEFINE_VARIABLE(double, TIME_INCREASING_RATIO) 
		
        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(EXTERNAL_APPLIED_FORCE)
		
		KRATOS_DEFINE_VARIABLE(int, INT_DUMMY_1)
		KRATOS_DEFINE_VARIABLE(int, INT_DUMMY_2)  
		KRATOS_DEFINE_VARIABLE(int, INT_DUMMY_3)  
		KRATOS_DEFINE_VARIABLE(int, INT_DUMMY_4)  
		KRATOS_DEFINE_VARIABLE(int, INT_DUMMY_5)  
		KRATOS_DEFINE_VARIABLE(int, INT_DUMMY_6)  
		KRATOS_DEFINE_VARIABLE(int, INT_DUMMY_7)  
		KRATOS_DEFINE_VARIABLE(int, INT_DUMMY_8)  
		KRATOS_DEFINE_VARIABLE(int, INT_DUMMY_9) 
		KRATOS_DEFINE_VARIABLE(double, DOUBLE_DUMMY_1)
		KRATOS_DEFINE_VARIABLE(double, DOUBLE_DUMMY_2)
		KRATOS_DEFINE_VARIABLE(double, DOUBLE_DUMMY_3)
		KRATOS_DEFINE_VARIABLE(double, DOUBLE_DUMMY_4)
		KRATOS_DEFINE_VARIABLE(double, DOUBLE_DUMMY_5)
		KRATOS_DEFINE_VARIABLE(double, DOUBLE_DUMMY_6)
		KRATOS_DEFINE_VARIABLE(double, DOUBLE_DUMMY_7)
		KRATOS_DEFINE_VARIABLE(double, DOUBLE_DUMMY_8)
		KRATOS_DEFINE_VARIABLE(double, DOUBLE_DUMMY_9)
	
		KRATOS_DEFINE_VARIABLE(double, REPRESENTATIVE_VOLUME)
		
        KRATOS_DEFINE_VARIABLE (int, CONTACT_MESH_OPTION )
        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(LOCAL_CONTACT_FORCE_LOW)                 
        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(LOCAL_CONTACT_FORCE_HIGH)
		
	KRATOS_DEFINE_VARIABLE(int, ACTIVATE_SEARCH)
	KRATOS_DEFINE_VARIABLE(int, CONCRETE_TEST_OPTION)
                
        KRATOS_DEFINE_VARIABLE( double,  CONTACT_FAILURE)
        KRATOS_DEFINE_VARIABLE( double,  CONTACT_SIGMA)
        KRATOS_DEFINE_VARIABLE( double,  CONTACT_TAU)
        KRATOS_DEFINE_VARIABLE( double,  FAILURE_CRITERION_STATE) 
	
	KRATOS_DEFINE_VARIABLE( double, CONTACT_SIGMA_MAX)  
	KRATOS_DEFINE_VARIABLE( double, CONTACT_SIGMA_MIN)
	KRATOS_DEFINE_VARIABLE( double, CONTACT_TAU_ZERO)  
	KRATOS_DEFINE_VARIABLE( double, CONTACT_INTERNAL_FRICC)
	
        KRATOS_DEFINE_VARIABLE (int, FAILURE_CRITERION_OPTION )        
                
        KRATOS_DEFINE_VARIABLE(int, SKIN_SPHERE)
        KRATOS_DEFINE_VARIABLE(double, EXPORT_SKIN_SPHERE)              
                
        KRATOS_DEFINE_VARIABLE(double, EXPORT_ID)
        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( DAMP_FORCES )
        //KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DRAG_FORCE)
        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(BUOYANCY)
        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( TOTAL_FORCES )

        //SURFACE VARIABLES
        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( PARTICLE_SURFACE_CONTACT_FORCES )
        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( PARTICLE_SURFACE_ROTATE_SPRING_MOMENT )
        KRATOS_DEFINE_VARIABLE(int, LIMIT_SURFACE_OPTION)
        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( SURFACE_NORMAL_DIR )
        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( SURFACE_POINT_COOR )
        KRATOS_DEFINE_VARIABLE(int, SURFACE_FRICC)

        KRATOS_DEFINE_VARIABLE (int, NEIGH_INITIALIZED)
        KRATOS_DEFINE_VARIABLE (int, TOTAL_CONTACTS )
        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(AUX_VEL)
            
        KRATOS_DEFINE_VARIABLE(Vector,     PARTICLE_BLOCK_CONTACT_FAILURE_ID)
        KRATOS_DEFINE_VARIABLE(Vector,     PARTICLE_BLOCK_CONTACT_FORCE)
        KRATOS_DEFINE_VARIABLE(Vector,     PARTICLE_BLOCK_IF_INITIAL_CONTACT)
        KRATOS_DEFINE_VARIABLE(WeakPointerVector<Element >,     NEIGHBOUR_PARTICLE_BLOCK_ELEMENTS)

  
        KRATOS_DEFINE_VARIABLE(int, VIRTUAL_MASS_OPTION)
        KRATOS_DEFINE_VARIABLE(int, UNIFORM_MATERIAL_OPTION)
        KRATOS_DEFINE_VARIABLE(double, NODAL_MASS_COEFF)

        KRATOS_DEFINE_VARIABLE(double, HISTORICAL_MIN_K)
        KRATOS_DEFINE_VARIABLE(int, CRITICAL_TIME_OPTION)

        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( PARTICLE_MOMENT )        
        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( PARTICLE_ROTATION_ANGLE )
        KRATOS_DEFINE_VARIABLE(double,  PARTICLE_MOMENT_OF_INERTIA)
                
        KRATOS_DEFINE_VARIABLE(double,  ROLLING_FRICTION)

        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( MAX_ROTA_MOMENT )
        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( INITIAL_ROTA_MOMENT )

        KRATOS_DEFINE_VARIABLE(double, ORIENTATION_REAL )
        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( ORIENTATION_IMAG )                

        KRATOS_DEFINE_VARIABLE(Vector,     INITIAL_AXES_TRACKING)
        KRATOS_DEFINE_VARIABLE(int,     plot_OPTIONS)

 
        KRATOS_DEFINE_VARIABLE(int, IF_BOUNDARY_ELEMENT)
        KRATOS_DEFINE_VARIABLE(Vector, IF_BOUNDARY_FACE)


class KratosDEMApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosDEMSpheresApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosDEMApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosDEMApplication();

    /// Destructor.
    virtual ~KratosDEMApplication() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual void Register();



    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "KratosDEMApplication";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        //KRATOS_WATCH("in KratosDEMApplication");
        //KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{
    //const SphericParticle mSphericParticle2D;
    const SphericParticle mSphericParticle3D;
    const SphericContinuumParticle mSphericContinuumParticle3D;
    const SphericSwimmingParticle mSphericSwimmingParticle3D;

    //const DEM_FEM_Particle mDEM_FEM_Particle2D;
    const DEM_FEM_Particle mDEM_FEM_Particle3D;
    
    const Particle_Contact_Element mParticleContactElement;
    const VariablesList mVariablesList;

    //       static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{
// 		const Elem2D   mElem2D;
// 		const Elem3D   mElem3D;


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KratosDEMApplication& operator=(KratosDEMApplication const& rOther);

    /// Copy constructor.
    KratosDEMApplication(KratosDEMApplication const& rOther);


    ///@}

}; // Class KratosDEMApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_DEM_APPLICATION_H_INCLUDED  defined 


