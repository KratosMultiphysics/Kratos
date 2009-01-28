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
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-05-16 13:59:01 $
//   Revision:            $Revision: 1.4 $ 12 November 2007 - 3D added
//
//		
// THIS PROCESS is INVENTED in order to CALCULATE THE PRESSURE and
//STORE it node-wise. 
//THIS IS FOR THE METHOD, where we do not calculate pressure force...

#if !defined(KRATOS_PRESSURE_CALCULATE_PROCESS_INCLUDED )
#define  KRATOS_PRESSURE_CALCULATE_PROCESS_INCLUDED



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
#include "ULF_application.h"


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

	class PressureCalculateProcess 
		: public Process
	{
	public:
		///@name Type Definitions
		///@{

		/// Pointer definition of PressureCalculateProcess
		KRATOS_CLASS_POINTER_DEFINITION(PressureCalculateProcess);

		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		PressureCalculateProcess(ModelPart& model_part, unsigned int domain_size)
			: mr_model_part(model_part),mdomain_size(domain_size)
		{
		}

		/// Destructor.
		virtual ~PressureCalculateProcess()
		{
		}


		///@}
		///@name Operators 
		///@{

		void operator()()
		{
			Execute();
		}


		///@}
		///@name Operations
		///@{

		virtual void Execute()
		{
			KRATOS_TRY
					
			ProcessInfo& proc_info = mr_model_part.GetProcessInfo();
			double dummy;
			//first initialize the pressure force to the old value

			//set the pressure to the old value
			for(ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin() ; 
				in != mr_model_part.NodesEnd() ; ++in)
			{  
				in->FastGetSolutionStepValue(PRESSURE)= 0.0;
			}

			for(ModelPart::ElementsContainerType::iterator im = mr_model_part.ElementsBegin() ; 
					im != mr_model_part.ElementsEnd() ; ++im)
			{  
				im->Calculate(PRESSURE,dummy,proc_info);
			}
KRATOS_WATCH("Execute of Pressure Calulate Process");
	/*		//
			if(mdomain_size == 2)
			{
				boost::numeric::ublas::bounded_matrix<double,3,2> DN_Dx;
				array_1d<double,3> N; 
				
				for(ModelPart::ElementsContainerType::iterator im = mr_model_part.ElementsBegin() ; 
					im != mr_model_part.ElementsEnd() ; ++im)
				{   
					//get the geometry
					Geometry< Node<3> >& geom = im->GetGeometry();

					//calculate derivatives
					double Area;
					GeometryUtils::CalculateGeometryData(geom, DN_Dx, N, Area);

					//calculate the divergence 
					const array_1d<double,3>& v0 = geom[0].FastGetSolutionStepValue(VELOCITY);
					const array_1d<double,3>& v1 = geom[1].FastGetSolutionStepValue(VELOCITY);
					const array_1d<double,3>& v2 = geom[2].FastGetSolutionStepValue(VELOCITY);

					double div_v = 	  DN_Dx(0,0)*v0[0] + DN_Dx(0,1)*v0[1]
									+ DN_Dx(1,0)*v1[0] + DN_Dx(1,1)*v1[1]
									+ DN_Dx(2,0)*v2[0] + DN_Dx(2,1)*v2[1];	
					double dp_el = K * dt * div_v * Area;

					geom[0].FastGetSolutionStepValue(PRESSURE) += dp_el*0.333333333333333333;
					geom[1].FastGetSolutionStepValue(PRESSURE) += dp_el*0.333333333333333333;
					geom[2].FastGetSolutionStepValue(PRESSURE) += dp_el*0.333333333333333333;			
				}
			}
			else
			{
				KRATOS_ERROR(std::logic_error,"pressure calculation 3D not implemented","");
	*/
					
			//divide by nodal area
			for(ModelPart::NodesContainerType::iterator in = mr_model_part.NodesBegin() ; 
				in != mr_model_part.NodesEnd() ; ++in)
			{  
				if( (in->GetValue(NEIGHBOUR_ELEMENTS)).size() != 0)
				{
					const double& ar=in->FastGetSolutionStepValue(NODAL_AREA);
					in->FastGetSolutionStepValue(PRESSURE)/=ar;
				
//					in->FastGetSolutionStepValue(PRESSURE)+=in->FastGetSolutionStepValue(PRESSURE,1);
				//here we set 0 -pressure to the lonely nodes
				}
				else
					in->FastGetSolutionStepValue(PRESSURE)=0.0;
			}
			
							
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
			return "PressureCalculateProcess";
		}

		/// Print information about this object.
		virtual void PrintInfo(std::ostream& rOStream) const
		{
			rOStream << "PressureCalculateProcess";
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
		double m_min_h;
		unsigned int mdomain_size;


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
//		PressureCalculateProcess& operator=(PressureCalculateProcess const& rOther);

		/// Copy constructor.
//		PressureCalculateProcess(PressureCalculateProcess const& rOther);


		///@}    

	}; // Class PressureCalculateProcess 

	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	inline std::istream& operator >> (std::istream& rIStream, 
		PressureCalculateProcess& rThis);

	/// output stream function
	inline std::ostream& operator << (std::ostream& rOStream, 
		const PressureCalculateProcess& rThis)
	{
		rThis.PrintInfo(rOStream);
		rOStream << std::endl;
		rThis.PrintData(rOStream);

		return rOStream;
	}
	///@} 


}  // namespace Kratos.

#endif // KRATOS_PRESSURE_CONTRIBUTION_PROCESS_INCLUDED  defined 


