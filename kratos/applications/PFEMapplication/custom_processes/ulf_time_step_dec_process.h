//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-01-16 17:25:31 $
//   Revision:            $Revision: 1.1 $ 
//
//  NOW FOR 2D ONLY!!! WRITE a 3D version of it!!!!

#if !defined(KRATOS_ULF_TIME_STEP_DEC_PROCESS_INCLUDED )
#define  KRATOS_ULF_TIME_STEP_DEC_PROCESS_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// External includes 

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "utilities/geometry_utilities.h"

//to compute the determinant
#include "utilities/math_utils.h"

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

	class UlfTimeStepDecProcess 
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of PushStructureProcess
		KRATOS_CLASS_POINTER_DEFINITION(UlfTimeStepDecProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		UlfTimeStepDecProcess(ModelPart& model_part)
			: mr_model_part(model_part)
		{
		}

		/// Destructor.
		virtual ~UlfTimeStepDecProcess()
		{
		}


		///@}
		///@name Operators 
		///@{

		void operator()()
		{
			
		}


		///@}
		///@name Operations
		///@{

		double EstimateDeltaTime(const double dt_max, const double domain_size)
		{
			double estimated_dt = 0.00;
			if(domain_size == 2)
				estimated_dt = EstimateDeltaTimeTemplated<2>(dt_max);
			else
				estimated_dt = EstimateDeltaTimeTemplated<3>(dt_max);
			return estimated_dt;
				
		}
		
		template< int TDim>
		double EstimateDeltaTimeTemplated(const double dt_max)
		{
		KRATOS_TRY
			
		double deltatime_new;
		double deltatime = dt_max;
		boost::numeric::ublas::bounded_matrix<double,TDim,TDim+1> aux;
		//this is for 2D only
		boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim> DN_DX;
		boost::numeric::ublas::bounded_matrix<double,TDim,TDim> Dv_dx;
		boost::numeric::ublas::bounded_matrix<double,TDim,TDim> J;
		boost::numeric::ublas::bounded_matrix<double,TDim,TDim> I = IdentityMatrix(TDim,TDim);
		boost::numeric::ublas::bounded_matrix<double,TDim,TDim> Zero = ZeroMatrix(TDim,TDim);
		array_1d<double,TDim+1> N;
		double Area;
		double detJ;	
		
		
		for(ModelPart::ElementsContainerType::iterator i = mr_model_part.ElementsBegin(); 
				i!=mr_model_part.ElementsEnd(); i++)
		{	
			//calculating the actual Jacobian
			for(unsigned int iii = 0; iii < i->GetGeometry().size(); iii++)
			{
				const array_1d<double,3>& v = i->GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY);
				for(unsigned int j=0; j <TDim; j++)
					aux(j,iii) = v[j];
			}
				
//std::cout<<"AUX from PK2 stress FCT"<<aux<<std::endl;
			GeometryUtils::CalculateGeometryData(i->GetGeometry(),DN_DX, N, Area);
				
			if (Area<=0) 
				KRATOS_ERROR(std::logic_error,"negative area at the moment of estimating the time step","");
			noalias(Dv_dx) = prod(aux,DN_DX);
			noalias(J) = I + deltatime*Dv_dx;
			
			detJ = MathUtils<double>::Det(J);
			std::cout<<"DETJ is "<<detJ<<std::endl;

			if (detJ<=0)				
			{
				deltatime_new=deltatime/(1-detJ); //x=-b/k
				deltatime=deltatime_new;
				noalias(J) = I + deltatime*Dv_dx;
				detJ = MathUtils<double>::Det(J);
				if (detJ<=0)				
				{
					deltatime_new=deltatime/(1-detJ); //x=-b/k
					deltatime=deltatime_new;				
				}
				
			}
		}
			
		return deltatime;
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
			return "UlfTimeStepDecProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "UlfTimeStepDecProcess";
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
//		UlfTimeStepDecProcess& operator=(UlfTimeStepDecProcess const& rOther);

		/// Copy constructor.
//		UlfTimeStepDecProcess(UlfTimeStepDecProcess const& rOther);


		///@}    

	}; // Class UlfTimeStepDecProcess 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		UlfTimeStepDecProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const UlfTimeStepDecProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_ULF_TIME_STEP_DEC_PROCESS_INCLUDED  defined 


