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
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:34 $
//   Revision:            $Revision: 1.2 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/process_info.h"
#include "python/add_process_info_to_python.h"
#include "containers/data_value_container.h"
#include "includes/convection_diffusion_settings.h"
#include "includes/radiation_settings.h"
namespace Kratos
{
	
namespace Python
{
// 	void SetModelPartName(ModelPart& rModelPart, std::string const& NewName)
// 	{
// 		rModelPart.Name() = NewName;
// 	}
// 	std::string GetModelPartName(ModelPart const& rModelPart)
// 	{
// 		return rModelPart.Name();
// 	}
// 	ProcessInfo& GetProcessInfo(ModelPart& rModelPart )
// 	{	  return rModelPart.GetProcessInfo();	 }
// 	void SetProcessInfo(ModelPart& rModelPart, ProcessInfo& NewProcessInfo)
// 	{	  rModelPart.SetProcessInfo(NewProcessInfo);	  }

	template< class TContainerType, class TVariableType > void MySetValueHelperFunction1(
                TContainerType& el, 
                const TVariableType& rVar,
                const typename TVariableType::Type& Data)
        {
            el.SetValue(rVar,Data);
        }
        
        template< class TContainerType, class TVariableType > 
                typename TVariableType::Type MyGetValueHelperFunction1( TContainerType& el, 
                        const TVariableType& rVar )
        {
            return el.GetValue(rVar);
        }
        
//
	void  AddProcessInfoToPython()
	{
		using namespace boost::python;
		
		class_<ProcessInfo, ProcessInfo::Pointer, bases<DataValueContainer>, boost::noncopyable>("ProcessInfo")
				.def(init<>())
				.def("CreateSolutionStepInfo", &ProcessInfo::CreateSolutionStepInfo)
// 				.def("CreateTimeStepInfo",(void (ProcessInfo::*)(std::size_t)) &ProcessInfo::CreateTimeStepInfo)
// 				.def("CreateTimeStepInfo",&ProcessInfo::CreateTimeStepInfo)
// 				.def("CloneTimeStepInfo",(void (ProcessInfo::*)(std::size_t) )&ProcessInfo::CloneTimeStepInfo)
// 				.def("SetAsTimeStepInfo",(void (ProcessInfo::*)()) &ProcessInfo::SetAsTimeStepInfo)
				.def(self_ns::str(self))
				;
	}
}  // namespace Python.

} // Namespace Kratos

