//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-05-27 12:29:23 $
//   Revision:            $Revision: 1.1 $ 
//
//  this process save structural elements in a separate list 

#if !defined(KRATOS_CFL_PROCESS_INCLUDED )
#define  KRATOS_CFL_PROCESS_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// External includes 


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/condition.h"
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h" 


namespace Kratos
{

	///@name Kratos Globals
	///@{ 

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
		Update the PRESSURE_FORCE on the nodes

		
	*/
	template<unsigned int TDim>
	class CFLProcess 
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of PushStructureProcess
		KRATOS_CLASS_POINTER_DEFINITION(CFLProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		CFLProcess(ModelPart& model_part)
			: mr_model_part(model_part)
		{
		}

		/// Destructor.
		virtual ~CFLProcess()
		{
		}


		///@}
		///@name Operators 
		///@{

//		void operator()()
//		{
//			SaveStructure();
//		}


		///@}
		///@name Operations
		///@{

		double EstimateTime(const double CFLnumber, const double max_dt)
		{
		KRATOS_TRY
		//initializee dt with max dt
		//initialize dt with incredible value
		double h, dt, glob_min_dt, nu, dummy;  
		if (TDim==2)
		{
			array_1d<double,3> N = ZeroVector(3); 
			array_1d<double,3> aux = ZeroVector(3); //dimension = number of nodes
			array_1d<double,3> vel = ZeroVector(3); //dimension = number of nodes
			boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX = ZeroMatrix(3,2);
			
			//initialize it with given value
			glob_min_dt=max_dt;


			dt=0.0;
			for(ModelPart::ElementsContainerType::iterator im = mr_model_part.ElementsBegin() ; 
					im != mr_model_part.ElementsEnd() ; ++im)
			{	  
				GeometryUtils::CalculateGeometryData(im->GetGeometry(),DN_DX,N,dummy);
				//direction of the height is stored in the auxilliary vector
				for (unsigned int i=0; i<3;i++)			
				{
				aux[0]=DN_DX(i,0);
				aux[1]=DN_DX(i,1);
				aux[2]=0.0;
				//and the value of the height: hi=1/norm(nablaNi)*(NablaNi/norm(nablaNi))
				h=1.0/norm_2(aux);
				//and now the velocity and viscosity
				vel=im->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
				nu=im->GetGeometry()[i].FastGetSolutionStepValue(VISCOSITY);
				//reuse dummy to temporarily store the scalar product of vel and height
				dt=1.0/(2.0*fabs(inner_prod(vel,aux)) + 4.0*nu/h*h);
				if(dt<glob_min_dt) 
					glob_min_dt=dt;
				}
			}
		}
		if (TDim==3)
		{
			array_1d<double,4> N = ZeroVector(4); //dimension = number of nodes
			array_1d<double,3> aux = ZeroVector(3); //dimension = space dim
			array_1d<double,3> vel = ZeroVector(3); //dimension = space dim
			boost::numeric::ublas::bounded_matrix<double,4,3> DN_DX = ZeroMatrix(4,3);
			
			//initialize it with given value
			glob_min_dt=max_dt;


			dt=0.0;
			for(ModelPart::ElementsContainerType::iterator im = mr_model_part.ElementsBegin() ; 
					im != mr_model_part.ElementsEnd() ; ++im)
			{	  
				GeometryUtils::CalculateGeometryData(im->GetGeometry(),DN_DX,N,dummy);
				//direction of the height is stored in the auxilliary vector
				for (unsigned int i=0; i<4;i++)			
				{
				aux[0]=DN_DX(i,0);
				aux[1]=DN_DX(i,1);
				aux[2]=DN_DX(i,2);
				//and the value of the height: hi=1/norm(nablaNi)*(NablaNi/norm(nablaNi))
				h=1.0/norm_2(aux);
				//and now the velocity and viscosity
				vel=im->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
				nu=im->GetGeometry()[i].FastGetSolutionStepValue(VISCOSITY);
				//reuse dummy to temporarily store the scalar product of vel and height
				dt=1.0/(fabs(inner_prod(vel,aux)) + 2.0*nu/h*h);
				if(dt<glob_min_dt) 
					glob_min_dt=dt;
				}
			}
		}
		if (dt<0.0)			
		KRATOS_ERROR(std::logic_error,  "NEGATIVE VALUE OF Time step estimated" , "");
				
		KRATOS_WATCH (glob_min_dt)
		return  ( CFLnumber*glob_min_dt);
		
		KRATOS_CATCH("")
		}


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
			return "CFLProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "CFLProcess";
		}

		/// Print object's data.
		virtual void PrintData(std::ostream& rOStream) const
		{
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


		///@} 
		///@name Member Variables 
		///@{ 
		//ModelPart& mr_fluid_model_part;
		//ModelPart& mr_structure_model_part;
		ModelPart& mr_model_part;
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
//		CFLProcess& operator=(CFLProcess const& rOther);

		/// Copy constructor.
//		CFLProcess(CFLProcess const& rOther);


		///@}    

	}; // Class CFLProcess 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	/*
	inline std::istream& operator >> (std::istream& rIStream,   
		CFLProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const CFLProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	*/
	///@} 


}  // namespace Kratos.

#endif // KRATOS_CFL_PROCESS_INCLUDED  defined 


