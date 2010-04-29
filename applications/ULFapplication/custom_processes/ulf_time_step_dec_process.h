/*
==============================================================================
KratosULFApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Pawel Ryzhakov
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
 
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2008-10-02 10:47:22 $
//   Revision:            $Revision: 1.7 $ 
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
		boost::numeric::ublas::bounded_matrix<double,TDim,TDim+1> aux_ac;
		//this is for 2D only
		boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim> DN_DX;
		boost::numeric::ublas::bounded_matrix<double,TDim,TDim> Dv_dx;
		boost::numeric::ublas::bounded_matrix<double,TDim,TDim> Da_dx;
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
			
			for(unsigned int iii = 0; iii < i->GetGeometry().size(); iii++)
			{
				const array_1d<double,3>& a = i->GetGeometry()[iii].FastGetSolutionStepValue(ACCELERATION);
				for(unsigned int j=0; j <TDim; j++)
					aux_ac(j,iii) = a[j];
			}
				
			//std::cout<<"AUX from PK2 stress FCT"<<aux<<std::endl;
			GeometryUtils::CalculateGeometryData(i->GetGeometry(),DN_DX, N, Area);
				
			if (Area<=0) 
				KRATOS_ERROR(std::logic_error,"negative area at the moment of estimating the time step","");
			
			noalias(Dv_dx) = prod(aux,DN_DX);
			noalias(Da_dx) = prod(aux_ac,DN_DX);

			noalias(J) = I + deltatime*Dv_dx +0.5*deltatime*deltatime*Da_dx;
			
			detJ = MathUtils<double>::Det(J);
//			std::cout<<"DETJ is "<<detJ<<std::endl;

			if (detJ<=0)				
			{
		
				deltatime_new=deltatime/(1-detJ); //x=-b/k
				deltatime=deltatime_new; 
				noalias(J) = I + deltatime*Dv_dx+0.5*deltatime*deltatime*Da_dx;
				detJ = MathUtils<double>::Det(J);
				if (detJ<=0)				
				{
					deltatime_new=deltatime/(1-detJ); //x=-b/k
					deltatime=deltatime_new;	
					noalias(J) = I + deltatime*Dv_dx+0.5*deltatime*deltatime*Da_dx;
					detJ = MathUtils<double>::Det(J);
		

						
					if (detJ<=0)				
						{
						for(unsigned int iii = 0; iii < i->GetGeometry().size(); iii++)
							
							deltatime_new=deltatime/(1-detJ); //x=-b/k
							deltatime=deltatime_new;	
							
						}
				}
				
			}
		}
		KRATOS_WATCH("Estimated delta time is")
		KRATOS_WATCH(deltatime)
		
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


