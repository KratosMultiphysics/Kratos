/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
//   Date:                $Date: 2007-03-06 10:30:34 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_SOLUTION_STEP_VARIABLE_INDEXING_PYTHON_H_INCLUDED )
#define  KRATOS_SOLUTION_STEP_VARIABLE_INDEXING_PYTHON_H_INCLUDED



// System includes
//#include <string>
//#include <iostream>
//#include <sstream> 


// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"


namespace Kratos
{

	namespace Python
	{

  using namespace boost::python;
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
		*/
		template<class TContainerType, class TVariableType>
		class SolutionStepVariableIndexingPython : public def_visitor<SolutionStepVariableIndexingPython<TContainerType, TVariableType> >
		{
		public:
			///@name Type Definitions
			///@{

			/// Pointer definition of SolutionStepVariableIndexingPython
			KRATOS_CLASS_POINTER_DEFINITION(SolutionStepVariableIndexingPython);

			///@}
			///@name Life Cycle 
			///@{ 

			/// Default constructor.
			SolutionStepVariableIndexingPython(){}

			/// Copy constructor.
			SolutionStepVariableIndexingPython(const SolutionStepVariableIndexingPython& rOther);

			/// Destructor.
			virtual ~SolutionStepVariableIndexingPython(){}


			///@}
			///@name Operators 
			///@{


			///@}
			///@name Operations
			///@{

			template <class TClassType>
				void visit(TClassType& ThisClass) const
			{
				ThisClass
					//.def("__contains__", &DataValueContainerHas)
					//.def("__setitem__", &SetBufferValue)
					//.def("__getitem__", &GetBufferValue)
					.def("Has", &DataValueContainerHas)
					.def("SetSolutionStepValue", &SetValue)
					.def("SetSolutionStepValue", &SetBufferValue)
					.def("GetSolutionStepValue", &GetValue)
					.def("GetSolutionStepValue", &GetBufferValue)
					//.def("__delitem__", &DataValueContainerDeleteValue)
					//.def("Erase", &DataValueContainerDeleteValue)
					;
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


			///@}      
			///@name Friends
			///@{


			///@}


		private:
			///@name Static Member Variables 
			///@{ 


			///@} 
			///@name Member Variables 
			///@{ 


			///@} 
			///@name Private Operators
			///@{ 


			///@} 
			///@name Private Operations
			///@{

  static void SetValue(TContainerType&  rData, TVariableType const& rV, typename TVariableType::Type const& rValue)
  {
    rData.GetSolutionStepValue(rV) = rValue;
  }

  static void SetBufferValue(TContainerType&  rData, TVariableType const& rV, typename TContainerType::IndexType SolutionStepIndex,typename TVariableType::Type const& rValue)
  {
    rData.GetSolutionStepValue(rV, SolutionStepIndex) = rValue;
  }

  static typename TVariableType::Type GetValue(TContainerType const& rData, TVariableType const& rV)
  {
    return rData.GetSolutionStepValue(rV);
  }
	
inline
  static typename TVariableType::Type GetBufferValue(TContainerType const& rData, TVariableType const& rV, typename TContainerType::IndexType SolutionStepIndex)
  {
    return rData.GetSolutionStepValue(rV, SolutionStepIndex);
  }
	
  //static void DataValueContainerDeleteValue(TContainerType& rData, TVariableType const& rV)
  //{
  //  rData.Erase(rV);
  //}
	
  static bool DataValueContainerHas(TContainerType const& rData, TVariableType const& rV, typename TContainerType::IndexType SolutionStepIndex)
  {
    return rData.SolutionStepsDataHas(rV);
  }
			



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
			SolutionStepVariableIndexingPython& operator=(const SolutionStepVariableIndexingPython& rOther);



			///@}    

		}; // Class SolutionStepVariableIndexingPython 

		///@} 

		///@name Type Definitions       
		///@{ 


		///@} 
		///@name Input and output 
		///@{ 

		///@} 

	}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_SOLUTION_STEP_VARIABLE_INDEXING_PYTHON_H_INCLUDED  defined 


