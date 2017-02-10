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
//#include "custom_elements/updated_lagrangian_UP.hpp"
#include "custom_elements/updated_lagrangian_quadrilateral.hpp"
//#include "custom_elements/updated_lagrangian_UP_quadrilateral.hpp"
//#include "custom_elements/total_lagrangian.hpp"

//constitutive laws
#include "custom_constitutive/hyperelastic_viscoplastic_3D_law.hpp"
#include "custom_constitutive/hyperelastic_viscoplastic_2D_plain_strain_law.hpp"
//flow rules
#include "custom_constitutive/flow_rules/viscoplastic_flow_rule.hpp"
#include "custom_constitutive/flow_rules/bingham_viscoplastic_flow_rule.hpp"
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

    //element
    KRATOS_DEFINE_VARIABLE(int, COUNTER )
    KRATOS_DEFINE_VARIABLE(int, MP_NUMBER )
    KRATOS_DEFINE_VARIABLE(int, MP_BOOL )
    KRATOS_DEFINE_VARIABLE(double, WEIGHT )
    KRATOS_DEFINE_VARIABLE(double, MP_MASS )
    KRATOS_DEFINE_VARIABLE(double, MP_DENSITY )
    KRATOS_DEFINE_VARIABLE(double, MP_VOLUME )
    KRATOS_DEFINE_VARIABLE(double, MP_KINETIC_ENERGY )
    KRATOS_DEFINE_VARIABLE(double, MP_STRAIN_ENERGY )
    KRATOS_DEFINE_VARIABLE(double, MP_TOTAL_ENERGY )
    
    
    
    //constitutive law
	
    KRATOS_DEFINE_VARIABLE(ConstitutiveLaw::Pointer, CONSTITUTIVE_LAW_POINTER )
    
    
    ////nodal dofs
    //KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_DISPLACEMENT )
    KRATOS_DEFINE_VARIABLE(double, AUX_R)
    KRATOS_DEFINE_VARIABLE(double, AUX_T)
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( AUX_R_VEL )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( AUX_T_VEL )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( AUX_R_ACC )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( AUX_T_ACC )
    KRATOS_DEFINE_VARIABLE(double, NODAL_LUMPED_MASS)

    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( AUX_VELOCITY )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( AUX_ACCELERATION )
    //MP element variable
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( GAUSS_COORD )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( MP_DISPLACEMENT )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( MP_VELOCITY )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( MP_ACCELERATION )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( AUX_MP_VELOCITY )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( AUX_MP_ACCELERATION )
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( MP_VOLUME_ACCELERATION )
    KRATOS_DEFINE_VARIABLE(Vector, MP_CAUCHY_STRESS_VECTOR )
    KRATOS_DEFINE_VARIABLE(Vector, MP_ALMANSI_STRAIN_VECTOR )
    KRATOS_DEFINE_VARIABLE(Vector, PREVIOUS_MP_CAUCHY_STRESS_VECTOR )
    KRATOS_DEFINE_VARIABLE(Vector, PREVIOUS_MP_ALMANSI_STRAIN_VECTOR )
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
 		//const UpdatedLagrangianUP mUpdatedLagrangianUP2D3N;
 		//const UpdatedLagrangianUP mUpdatedLagrangianUP3D4N;
 		const UpdatedLagrangianQuadrilateral mUpdatedLagrangian2D4N;
 		//const UpdatedLagrangianUPQuadrilateral mUpdatedLagrangianUP2D4N;
 		
 		//const TotalLagrangian mTotalLagrangian2D3N;
 		//const TotalLagrangian mTotalLagrangian3D4N;
		
        //constitutive laws

        const HyperElasticViscoplastic3DLaw                mHyperElasticViscoplastic3DLaw;
        const HyperElasticViscoplasticPlaneStrain2DLaw     mHyperElasticViscoplasticPlaneStrain2DLaw;

		
        //Flow Rules
        //const NonLinearAssociativePlasticFlowRule     mNonLinearAssociativePlasticFlowRule;
        //const LinearAssociativePlasticFlowRule        mLinearAssociativePlasticFlowRule;
        //const IsotropicDamageFlowRule                 mIsotropicDamageFlowRule;
        const ViscoplasticFlowRule                    mViscoplasticFlowRule;
        const BinghamViscoplasticFlowRule             mBinghamViscoplasticFlowRule;
        
        //Yield Criteria
        //const MisesHuberYieldCriterion                mMisesHuberYieldCriterion;
        //const SimoJuYieldCriterion                    mSimoJuYieldCriterion;
        
        //Hardening Laws
        //const NonLinearIsotropicKinematicHardeningLaw mNonLinearIsotropicKinematicHardeningLaw;
        //const LinearIsotropicKinematicHardeningLaw    mLinearIsotropicKinematicHardeningLaw;
        //const ExponentialDamageHardeningLaw           mExponentialDamageHardeningLaw;

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


