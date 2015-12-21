//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_PFEM_2_APPLICATION_H_INCLUDED )
#define  KRATOS_PFEM_2_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "custom_utilities/pfem_particle.h" //hahaha
#include "custom_utilities/pfem_particle_fluidonly.h" //hahaha

//#include "custom_elements/2fluid_2d.h" //including the file for the element
//#include "custom_elements/2fluid_3d.h" //including the file for the element
//#include "custom_elements/fluid_phase_2d.h" //including the file for the element
#include "custom_elements/fractional_step_pfem_2_2d.h" //including the file for the element
#include "custom_elements/fractional_step_pfem_2_3d.h" //including the file for the element
#include "custom_elements/monolithic_2fluid_2d.h" //including the file for the element
#include "custom_elements/monolithic_2fluid_3d.h" //including the file for the element
#include "custom_elements/nonewtonian_2fluid_2d.h" //including the file for the element
#include "custom_elements/nonewtonian_2fluid_3d.h" //including the file for the element
//#include "custom_elements/monolithic_3fluid_2d.h" //including the file for the element
//#include "custom_elements/monolithic_3fluid_3d.h" //including the file for the element
#include "custom_elements/monolithic_2fluid_2d_partintegration.h" //including the file for the element
#include "custom_elements/monolithic_2fluid_3d_partintegration.h" //including the file for the element
//#include "custom_elements/fsi_2d.h" //including the file for the element
//#include "custom_elements/fsi_3d.h" //including the file for the element
//#include "custom_elements/no_particles_solid_only_2d.h" //including the file for the element

#include "custom_conditions/fixed_velocity_2d.h" //the condition
#include "custom_conditions/fixed_velocity_3d.h" //the condition
#include "custom_conditions/fixed_pressure_2d.h" //the condition
#include "custom_conditions/fixed_pressure_3d.h" //the condition
#include "custom_conditions/autoslip_inlet_3d.h" //the condition
//#include "custom_conditions/water_fixed_velocity_2d.h" //the condition


namespace Kratos
{

	///@name Kratos Globals
	///@{ 

	// Variables definition 
	KRATOS_DEFINE_VARIABLE(double, PRESS_GRADIENT_JUMP)
	KRATOS_DEFINE_VARIABLE(double, PRESS_DISCONTINUITY)
	KRATOS_DEFINE_VARIABLE(double, INV_LAPLACIAN_ENRICH)
	KRATOS_DEFINE_VARIABLE(double, ENRICH_RHS)
	KRATOS_DEFINE_VARIABLE(double, G_VALUE)
	KRATOS_DEFINE_VARIABLE(double, GRADIENT_DISCONTINUITY)
	KRATOS_DEFINE_VARIABLE(double, PREVIOUS_ITERATION_PRESSURE)
	KRATOS_DEFINE_VARIABLE(double, FIRST_ITERATION_PRESSURE)
	KRATOS_DEFINE_VARIABLE(double, VELOCITY_OVER_ELEM_SIZE)
	KRATOS_DEFINE_VARIABLE(double, MEAN_SIZE)
	KRATOS_DEFINE_VARIABLE(double, MEAN_VELOCITY_DIFFERENCE)
	KRATOS_DEFINE_VARIABLE(double, SPECIFIC_HEAT_CAPACITY_WATER)
	KRATOS_DEFINE_VARIABLE(double, SPECIFIC_HEAT_CAPACITY_AIR)
	KRATOS_DEFINE_VARIABLE(double, DELTA_TEMPERATURE)
	KRATOS_DEFINE_VARIABLE(double, AVAILABLE_AIR_VOLUME)
	KRATOS_DEFINE_VARIABLE(double, AVAILABLE_UNBURNED_AIR_VOLUME)
	KRATOS_DEFINE_VARIABLE(double, OXYGEN_FRACTION)
	KRATOS_DEFINE_VARIABLE(double, CORRECTED_DISTANCE)
	KRATOS_DEFINE_VARIABLE(double, SOLID_PRESSURE)
	KRATOS_DEFINE_VARIABLE(double, SOLID_YP)
	KRATOS_DEFINE_VARIABLE(double, WATER_DISTANCE)
	KRATOS_DEFINE_VARIABLE(double, VOLUMETRIC_STRAIN)
	KRATOS_DEFINE_VARIABLE(double, ELASTIC_PRESSURE)
	
	KRATOS_DEFINE_VARIABLE(bool, USEFUL_ELEMENT_FOR_COMBUSTION)
	KRATOS_DEFINE_VARIABLE(Vector, ENRICH_LHS_ROW_3D)
	KRATOS_DEFINE_VARIABLE(Vector, WATER_GAUSS_POINT)
	KRATOS_DEFINE_VARIABLE(double, WATER_VOLUME)
	KRATOS_DEFINE_VARIABLE(Vector, ELEMENT_MEAN_STRESS)
	typedef PointerVector< PFEM_Particle, PFEM_Particle*, std::vector<PFEM_Particle*> > ParticlePointerVector;
	KRATOS_DEFINE_VARIABLE( ParticlePointerVector , PARTICLE_POINTERS)	
	typedef PointerVector< PFEM_Particle_Fluid, PFEM_Particle_Fluid*, std::vector<PFEM_Particle_Fluid*> > FluidParticlePointerVector;
	KRATOS_DEFINE_VARIABLE( FluidParticlePointerVector , FLUID_PARTICLE_POINTERS)	
	KRATOS_DEFINE_VARIABLE(int, NUMBER_OF_PARTICLES)
	KRATOS_DEFINE_VARIABLE(int, NUMBER_OF_PARTICLES_AUX)
	KRATOS_DEFINE_VARIABLE(int, NUMBER_OF_WATER_PARTICLES)
	KRATOS_DEFINE_VARIABLE(int, NUMBER_OF_FLUID_PARTICLES)
	KRATOS_DEFINE_VARIABLE(int, PARTICLE_POINTERS_OFFSET)	
	KRATOS_DEFINE_VARIABLE(int, WATER_PARTICLE_POINTERS_OFFSET)
	KRATOS_DEFINE_VARIABLE(int, USE_PRESS_PROJ)
	//KRATOS_DEFINE_VARIABLE(double, IS_AIR)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(ENRICH_LHS_ROW)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(ENRICH_PRESS_PROJ_NEGATIVE)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(ENRICH_PRESS_PROJ_POSITIVE)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(SURFACE_NORMAL)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(SURFACE_COORDINATES)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(PRESS_PROJ_NO_RO)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DELTA_VELOCITY)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(WATER_VELOCITY)
	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(WATER_MESH_VELOCITY)

	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(PROJECTED_VELOCITY)
	KRATOS_DEFINE_VARIABLE(double, VOLUME_CORRECTION)
	KRATOS_DEFINE_VARIABLE(double, INLET_VELOCITY)

	
//	KRATOS_DEFINE_VARIABLE(double, NODAL_AREA)


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
	class KratosPFEM2Application : public KratosApplication
	{
	public:
		///@name Type Definitions
		///@{
		

		/// Pointer definition of KratosPFEM2Application
		KRATOS_CLASS_POINTER_DEFINITION(KratosPFEM2Application);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		KratosPFEM2Application();

		/// Destructor.
		virtual ~KratosPFEM2Application(){}


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
			return "KratosPFEM2Application";
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
 	
 		//const FluidPhasePFEM22D   mFluidPhasePFEM22D; 
  		const FractionalStepPFEM22D   mFractionalStepPFEM22D; 
  		const FractionalStepPFEM23D   mFractionalStepPFEM23D; 
		const MonolithicPFEM22D   mMonolithicPFEM22D; 
		const MonolithicPFEM23D   mMonolithicPFEM23D; 
		const NoNewtonianMonolithicPFEM22D   mNoNewtonianMonolithicPFEM22D; 
		const NoNewtonianMonolithicPFEM23D   mNoNewtonianMonolithicPFEM23D; 
 		//const Monolithic3FluidPFEM22D mMonolithic3FluidPFEM22D;
 		//const Monolithic3FluidPFEM23D mMonolithic3FluidPFEM23D;
 		const MonolithicAutoSlipPFEM22D   mMonolithicAutoSlipPFEM22D; 
 		const MonolithicAutoSlipPFEM23D   mMonolithicAutoSlipPFEM23D; 
 		//const FsiPFEM22D   mFsiPFEM22D; 
 		//const FsiPFEM23D   mFsiPFEM23D; 
 		//const NoParticlesSolidOnlyPFEM22D   mNoParticlesSolidOnlyPFEM22D; 
 		
		const FixedVelocity2D   mFixedVelocity2D; 
		const FixedVelocity3D   mFixedVelocity3D;
 		const FixedPressure2D   mFixedPressure2D; 
		const FixedPressure3D   mFixedPressure3D; 
		const MonolithicAutoSlipInlet3D   mMonolithicAutoSlipInlet3D; 
		//const WaterFixedVelocity2D   mWaterFixedVelocity2D; 
 		

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
		KratosPFEM2Application& operator=(KratosPFEM2Application const& rOther);

		/// Copy constructor.
		KratosPFEM2Application(KratosPFEM2Application const& rOther);


		///@}    

	}; // Class KratosPFEM2Application 

	///@} 


	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 

	///@} 


}  // namespace Kratos.

#endif // KRATOS_PFEM_2_APPLICATION_H_INCLUDED  defined 


