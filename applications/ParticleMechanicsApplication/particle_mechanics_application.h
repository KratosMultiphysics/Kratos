//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_PARTICLE_MECHANICS_APPLICATION_H_INCLUDED )
#define  KRATOS_PARTICLE_MECHANICS_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 
#include "solid_mechanics_application.h"
//#include "structural_application.h"
//#include "pfem_solid_mechanics_application.h"
// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/constitutive_law.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

#include "includes/condition.h"   //#include "costum_conditions/nuova_condizione.h nel caso non sia gia stata definita
#include "includes/ublas_interface.h"

#include "containers/flags.h"



//element
#include "custom_elements/updated_lagrangian.hpp"
//#include "custom_elements/updated_lagrangian_quadrilateral.hpp"
//#include "custom_elements/total_lagrangian.hpp"
namespace Kratos
{
	///@name Type Definitions
	///@{
	typedef array_1d<double,3> Vector3;
	typedef array_1d<double,6> Vector6;
	///@}

	///@name Kratos Globals
	///@{ 

	// Variables definition 

    //solution
    KRATOS_DEFINE_VARIABLE(double, RAYLEIGH_ALPHA )
    KRATOS_DEFINE_VARIABLE(double, RAYLEIGH_BETA )
    //element
    KRATOS_DEFINE_VARIABLE(int, COUNTER )
    KRATOS_DEFINE_VARIABLE(int, MP_NUMBER )
    KRATOS_DEFINE_VARIABLE(double, WEIGHT )
    KRATOS_DEFINE_VARIABLE(double, MP_MASS )
    KRATOS_DEFINE_VARIABLE(double, MP_DENSITY )
    KRATOS_DEFINE_VARIABLE(double, MP_VOLUME )
    KRATOS_DEFINE_VARIABLE(double, MP_KINETIC_ENERGY )
    KRATOS_DEFINE_VARIABLE(double, MP_STRAIN_ENERGY )
    //KRATOS_DEFINE_VARIABLE(double, NODAL_MASS )
    
    KRATOS_DEFINE_VARIABLE(Vector, EXTERNAL_FORCES_VECTOR )
    KRATOS_DEFINE_VARIABLE(Vector, INTERNAL_FORCES_VECTOR )
    
    KRATOS_DEFINE_VARIABLE(Vector, CAUCHY_STRESS_VECTOR )
    KRATOS_DEFINE_VARIABLE(Vector, PK2_STRESS_VECTOR )   
    KRATOS_DEFINE_VARIABLE(Vector, GREEN_LAGRANGE_STRAIN_VECTOR )
    KRATOS_DEFINE_VARIABLE(Vector, ALMANSI_STRAIN_VECTOR )

    KRATOS_DEFINE_VARIABLE(Matrix, CAUCHY_STRESS_TENSOR );
    KRATOS_DEFINE_VARIABLE(Matrix, PK2_STRESS_TENSOR );
    KRATOS_DEFINE_VARIABLE(Matrix, GREEN_LAGRANGE_STRAIN_TENSOR );
    KRATOS_DEFINE_VARIABLE(Matrix, ALMANSI_STRAIN_TENSOR )
    
    KRATOS_DEFINE_VARIABLE(Matrix, MATERIAL_STIFFNESS_MATRIX )
    KRATOS_DEFINE_VARIABLE(Matrix, GEOMETRIC_STIFFNESS_MATRIX )
    
    
    //constitutive law
    KRATOS_DEFINE_VARIABLE(ConstitutiveLaw::Pointer, CONSTITUTIVE_LAW_POINTER )
    KRATOS_DEFINE_VARIABLE(Matrix, CONSTITUTIVE_MATRIX )
    KRATOS_DEFINE_VARIABLE(Matrix, DEFORMATION_GRADIENT )
    KRATOS_DEFINE_VARIABLE(double, DETERMINANT_F )
    
    ////nodal dofs
    //KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_DISPLACEMENT )
    
    //MP element variable
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( GAUSS_COORD )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( MP_DISPLACEMENT )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( MP_VELOCITY )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( MP_ACCELERATION )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( MP_VOLUME_ACCELERATION )
    KRATOS_DEFINE_VARIABLE(Vector, MP_CAUCHY_STRESS_VECTOR )
    KRATOS_DEFINE_VARIABLE(Vector, MP_ALMANSI_STRAIN_VECTOR )
    KRATOS_DEFINE_VARIABLE(Matrix, MP_CONSTITUTIVE_MATRIX )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_AUX )
    //grid node variable
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( NODAL_MOMENTUM )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( NODAL_INERTIA )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( NODAL_INTERNAL_FORCE )
	///@} 
	///@name Type Definitions
	///@{ 

	///@} 
	///@name  Enum's
	///@{

	///@}
	///@name  Functions 
	///@{

	///@}
	///@name Kratos Classes
	///@{

	/// Short class definition.
	/** Detail class definition.
	*/
	class KratosParticleMechanicsApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{
		

		/// Pointer definition of KratosParticleMechanicsApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosParticleMechanicsApplication);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		KratosParticleMechanicsApplication();

		/// Destructor.
		virtual ~KratosParticleMechanicsApplication(){}


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
			return "KratosParticleMechanicsApplication";
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
      	KRATOS_WATCH("in my application");
      	KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
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



		//       static const ApplicationCondition  msApplicationCondition; 

		///@} 
		///@name Member Variables 
		///@{ 
 		const UpdatedLagrangian mUpdatedLagrangian2D3N;
 		const UpdatedLagrangian mUpdatedLagrangian3D4N;
 		//const UpdatedLagrangianQuad mUpdatedLagrangianQuad2D4N;
 		//const TotalLagrangian mTotalLagrangian2D3N;
 		//const TotalLagrangian mTotalLagrangian3D4N;


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
		KratosParticleMechanicsApplication& operator=(KratosParticleMechanicsApplication const& rOther);

		/// Copy constructor.
		KratosParticleMechanicsApplication(KratosParticleMechanicsApplication const& rOther);


		///@}    

	}; // Class KratosParticleMechanicsApplication 

	///@} 


	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 

	///@} 


}  // namespace Kratos.

#endif // KRATOS_PARTICLE_MECHANICS_APPLICATION_H_INCLUDED  defined 


