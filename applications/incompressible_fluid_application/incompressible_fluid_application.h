//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2009-01-15 11:11:35 $
//   Revision:            $Revision: 1.16 $
//
//


#if !defined(KRATOS_KRATOS_RICCARDOS_APPLICATION_H_INCLUDED )
#define  KRATOS_KRATOS_RICCARDOS_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

#include "custom_elements/fluid_3d.h"
#include "custom_elements/fluid_2d.h"

#include "custom_conditions/fluid3d_neumann.h"
#include "custom_conditions/monolithic2d_neumann.h"
//#include "custom_conditions/monolithic3d_neumann.h"

#include "custom_conditions/no_slip_condition_2d.h"
//#include "custom_conditions/no_slip_fract_step.h"

#include "includes/variables.h"
#include "includes/condition.h"

#include "custom_elements/fluid_2dcoupled.h"
#include "custom_elements/fluid_3dcoupled.h"

#include "custom_elements/NDfluid_2d.h"
#include "custom_elements/NDfluid_3d.h"

#include "custom_elements/NDfluid_2d_CrankNicolson.h"

#include "custom_elements/asgs_2d.h"
#include "custom_elements/asgs_pr_dc_2d.h"
#include "custom_elements/asgs_compressible_2d.h"
 #include "custom_elements/asgs_comp_pr_dc_2d.h"

#include "custom_elements/asgs_3d.h"
#include "custom_elements/asgs_compressible_3d.h"
 #include "custom_elements/asgs_comp_pr_dc_3d.h"

#include "custom_elements/fluid_2dGLS_expl.h"
#include "custom_elements/fluid_3dGLS_expl.h"

#include "custom_elements/fluid_2dGLS_expl_comp.h"
#include "custom_elements/fluid_3dGLS_expl_comp.h"

#include "custom_conditions/proj_dirichlet_cond.h"

#include "custom_elements/fluid_2d_split.h"


//#include "custom_conditions/abc_2d.h"


#include "custom_elements/nonewtonian_asgs_2d.h"
#include "custom_elements/couette_nonnewtonian_asgs_2d.h"
#include "custom_elements/nonewtonian_asgs_3d.h"

//#include "custom_elements/asgs_2d_particle.h"

namespace Kratos
{

	///@name Kratos Globals
	///@{ 

	// Variables definition 
//	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(FRACT_VEL)
//	//KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(PRESS_PROJ)
//	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(CONV_PROJ)

	KRATOS_DEFINE_VARIABLE( double, MACH_NUMBER )
	KRATOS_DEFINE_VARIABLE( double, PRESSURE_COEFFICIENT )
//	KRATOS_DEFINE_VARIABLE( double, PRESSURE_OLD_IT )
//	KRATOS_DEFINE_VARIABLE( Vector, BDF_COEFFICIENTS )
//	KRATOS_DEFINE_VARIABLE(double, NODAL_MASS)
//	KRATOS_DEFINE_VARIABLE(int, AUX_INDEX)
//	KRATOS_DEFINE_VARIABLE(double, EXTERNAL_PRESSURE)
// 	KRATOS_DEFINE_VARIABLE(double, DIAMETER)
	KRATOS_DEFINE_VARIABLE(double, PERMEABILITY_INV)
	
	//for disabling elements (e.g in the interaction)
	KRATOS_DEFINE_VARIABLE(int, DISABLE)
	
	KRATOS_DEFINE_VARIABLE(double, DENSITY_AIR )

	KRATOS_DEFINE_VARIABLE(double, AIR_SOUND_VELOCITY )
	KRATOS_DEFINE_VARIABLE(double, SOUND_VELOCITY )



	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(RHS_VECTOR)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(AUX_VECTOR)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(AUX_VEL)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(AUX_VEL1)
	
			
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
	class KratosIncompressibleFluidApplication : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{
		

		/// Pointer definition of KratosIncompressibleFluidApplication
		KRATOS_CLASS_POINTER_DEFINITION(KratosIncompressibleFluidApplication);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		KratosIncompressibleFluidApplication();

		/// Destructor.
		virtual ~KratosIncompressibleFluidApplication(){}


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
			return "KratosIncompressibleFluidApplication";
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
      	KRATOS_WATCH("in KratosIncompressibleFluidApplication");
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


		///@} 
		///@name Private Operators
		///@{ 
		const Fluid3D  mFluid3D; 
		const Fluid2D  mFluid2D; 

		const Fluid2DCoupled  mFluid2DCoupled; 
		const Fluid3DCoupled  mFluid3DCoupled; 

		const Fluid3DNeumann  mFluid3DNeumann;  



		const NDFluid2D  mNDFluid2D; 
		const NDFluid3D	 mNDFluid3D;
		const NDFluid2DCrankNicolson	 mNDFluid2DCrankNicolson;

		const ASGS2D  mASGS2D; 
		const ASGSPRDC2D  mASGSPRDC2D;

		const ASGSCompressible2D  mASGSCompressible2D; 
		      
		const ASGS3D  mASGS3D;
		const ASGSCompressible3D  mASGSCompressible3D; 



		//const ASGS2DPARTICLE  mASGS2DPARTICLE; 

 		const ASGSCOMPPRDC2D mASGSCOMPPRDC2D;
 		const ASGSCOMPPRDC3D mASGSCOMPPRDC3D;

 		const Monolithic2DNeumann  mMonolithic2DNeumann; 
 		//const Monolithic3DNeumann  mMonolithic3DNeumann; 


		const Fluid2DGLS_expl  mFluid2DGLS_expl;
		const Fluid3DGLS_expl  mFluid3DGLS_expl;

		const Fluid2DGLS_expl_comp  mFluid2DGLS_expl_comp;  
		const Fluid3DGLS_expl_comp  mFluid3DGLS_expl_comp;  

		const ProjDirichletCond mProjDirichletCond;

		const Fluid2DSplit mFluid2DSplit;

		const NoNewtonianASGS2D mNoNewtonianASGS2D;
		const CouetteNonNewtonianASGS2D mCouetteNonNewtonianASGS2D;

		const NoNewtonianASGS3D mNoNewtonianASGS3D;


		const NoSlipCondition2D mNoSlipCondition2D;
// 		const NoSlipFractStep mNoSlipFractStep;

		//const ABC2D mABC2D;

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
		KratosIncompressibleFluidApplication& operator=(KratosIncompressibleFluidApplication const& rOther);

		/// Copy constructor.
		KratosIncompressibleFluidApplication(KratosIncompressibleFluidApplication const& rOther);


		///@}    

	}; // Class KratosIncompressibleFluidApplication 

	///@} 


	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 

	///@} 


}  // namespace Kratos.

#endif // KRATOS_KRATOS_RICCARDOS_APPLICATION_H_INCLUDED  defined 


