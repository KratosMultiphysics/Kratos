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
#include "custom_elements/DEM_FEM_Particle.h"
#include "custom_elements/Particle_Contact_Element.h"


const long double pi = 3.141592653589793238462643383279;

namespace Kratos
{
    /* Define In Global variables.h
        KRATOS_DEFINE_VARIABLE(double,  DELTA_TIME);
        KRATOS_DEFINE_VARIABLE(Vector,     PARTICLE_ROTATE_SPRING_FAILURE_TYPE)
        typedef vector<array_1d<double,3> > VectorArray3Double;
        KRATOS_DEFINE_VARIABLE( VectorArray3Double, PARTICLE_ROTATE_SPRING_MOMENT )
     */

	
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
		
		KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( OLD_COORDINATES )
		
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
        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( TOTAL_FORCES )
    
        KRATOS_DEFINE_VARIABLE (int, NEIGH_INITIALIZED)
        KRATOS_DEFINE_VARIABLE (int, TOTAL_CONTACTS )
        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(AUX_VEL)

        KRATOS_DEFINE_VARIABLE(Vector,     PARTICLE_BLOCK_CONTACT_FAILURE_ID)
        KRATOS_DEFINE_VARIABLE(Vector,     PARTICLE_BLOCK_CONTACT_FORCE)
        KRATOS_DEFINE_VARIABLE(Vector,     PARTICLE_BLOCK_IF_INITIAL_CONTACT)
        KRATOS_DEFINE_VARIABLE(WeakPointerVector<Element >,     NEIGHBOUR_PARTICLE_BLOCK_ELEMENTS)

  
        KRATOS_DEFINE_VARIABLE(int, VIRTUAL_MASS_OPTION)
        KRATOS_DEFINE_VARIABLE(double, NODAL_MASS_COEFF)

        KRATOS_DEFINE_VARIABLE(double, HISTORICAL_MIN_K)
        KRATOS_DEFINE_VARIABLE(int, CRITICAL_TIME_OPTION)

        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( PARTICLE_MOMENT )        
        KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( PARTICLE_ROTATION_ANGLE )       
        KRATOS_DEFINE_VARIABLE(double,  PARTICLE_MOMENT_OF_INERTIA)
        KRATOS_DEFINE_VARIABLE(double,  ROLLING_FRICTION)

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


