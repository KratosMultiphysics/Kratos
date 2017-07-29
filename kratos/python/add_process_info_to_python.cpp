//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
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

ProcessInfo::Pointer ProcessInfoGetPreviousSolutionStepInfo(ProcessInfo & rProcessInfo)
{
	return rProcessInfo.pGetPreviousSolutionStepInfo();
}

//
void  AddProcessInfoToPython()
{
    using namespace boost::python;

    class_<ProcessInfo, ProcessInfo::Pointer, bases<DataValueContainer, Flags>, boost::noncopyable>("ProcessInfo")
    .def(init<>())
    .def("CreateSolutionStepInfo", &ProcessInfo::CreateSolutionStepInfo)
	.def("GetPreviousSolutionStepInfo", ProcessInfoGetPreviousSolutionStepInfo)
// 				.def("CreateTimeStepInfo",(void (ProcessInfo::*)(std::size_t)) &ProcessInfo::CreateTimeStepInfo)
// 				.def("CreateTimeStepInfo",&ProcessInfo::CreateTimeStepInfo)
// 				.def("CloneTimeStepInfo",(void (ProcessInfo::*)(std::size_t) )&ProcessInfo::CloneTimeStepInfo)
// 				.def("SetAsTimeStepInfo",(void (ProcessInfo::*)()) &ProcessInfo::SetAsTimeStepInfo)
    .def(self_ns::str(self))
    ;
}
}  // namespace Python.

} // Namespace Kratos

