// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


// External includes
#include <boost/python.hpp>
#include "boost/python/suite/indexing/map_indexing_suite.hpp"

// System includes
#include <sstream>

// Project includes
#include "includes/define.h"
#include "includes/kernel.h"
#include "python/add_kernel_to_python.h"

namespace Kratos
{
namespace Python
{
using namespace boost::python;

template< class TVariableType >
bool HasVariable(Kernel& rKernel, const std::string& variable_name)
{
    return KratosComponents<TVariableType>::Has(variable_name);
}

template< class TVariableType >
const TVariableType& GetVariable(Kernel& rKernel, const std::string& variable_name)
{
    if(KratosComponents<TVariableType>::Has(variable_name))
    {
        return KratosComponents<TVariableType>::Get(variable_name);
    }

    return TVariableType::StaticObject();
}

template< class TVariableType >
void PrintVariablesName(Kernel& rKernel)
{
	KratosComponents<TVariableType> kratos_components;
	kratos_components.PrintData(std::cout);
}
template< class TVariableType >
std::string GetVariableNames(Kernel& rKernel)
{
	KratosComponents<TVariableType> kratos_components;
	std::stringstream buffer;
	kratos_components.PrintData(buffer);
	return buffer.str();
}

void  AddKernelToPython()
{
//class_<std::map<std::string, const VariableData*> >("VariableDataMap")
// .def(map_indexing_suite<std::map<std::string, const VariableData> >())
//  ;
    class_<Kernel, Kernel::Pointer, boost::noncopyable >("Kernel")
    .def("Initialize",&Kernel::Initialize)
    .def("AddApplication",&Kernel::AddApplication,with_custodian_and_ward<1,2>()) // Note: custodian and ward to be checked. Pooyan.
    .def("InitializeApplication",&Kernel::InitializeApplication,with_custodian_and_ward<1,2>()) // Note: custodian and ward to be checked. Pooyan.
    //.def("",&Kernel::Initialize)
    .def("HasBoolVariable",HasVariable< Variable<bool> >)
    .def("GetBoolVariable",GetVariable< Variable<bool> >,return_internal_reference<>())
    .def("HasIntVariable",HasVariable< Variable<int> >)
    .def("GetIntVariable",GetVariable< Variable<int> >,return_internal_reference<>())
    .def("HasUnsignedIntVariable",HasVariable< Variable<unsigned int> >)
    .def("GetUnsignedIntVariable",GetVariable< Variable<unsigned int> >,return_internal_reference<>())
    .def("HasDoubleVariable",HasVariable< Variable<double> >)
    .def("GetDoubleVariable",GetVariable< Variable<double> >,return_internal_reference<>())
    .def("HasArrayVariable",HasVariable< Variable< array_1d<double,3> > >)
    .def("GetArrayVariable",GetVariable< Variable< array_1d<double,3> > >,return_internal_reference<>())
    .def("HasVectorVariable",HasVariable< Variable<Vector> >)
    .def("GetVectorVariable",GetVariable< Variable<Vector> >,return_internal_reference<>())
    .def("HasMatrixVariable",HasVariable< Variable<Matrix> >)
    .def("GetMatrixVariable",GetVariable< Variable<Matrix> >,return_internal_reference<>())
    .def("HasStringVariable",HasVariable< Variable<std::string> >)
    .def("GetStringVariable",GetVariable< Variable<std::string> >,return_internal_reference<>())
    .def("HasVariableComponent",HasVariable< VariableComponent< VectorComponentAdaptor< array_1d<double,3> > > > )
    .def("GetVariableComponent",GetVariable< VariableComponent< VectorComponentAdaptor< array_1d<double,3> > > > ,return_internal_reference<>())
    .def("HasFlagsVariable",HasVariable< Variable<Flags> >)
    .def("GetFlagsVariable",GetVariable< Variable<Flags> >,return_internal_reference<>())
    .def("HasVariableData",HasVariable< VariableData >)
	.def("PrintAllVariables", PrintVariablesName<VariableData>)
	.def("PrintBoolVariables", PrintVariablesName<Variable<bool> >)
	.def("PrintIntVariables", PrintVariablesName<Variable<int> >)
	.def("PrintUnsignedIntVariables", PrintVariablesName<Variable<int> >)
	.def("PrintDoubleVariables", PrintVariablesName<Variable <double> >)
	.def("PrintArrayVariables", PrintVariablesName<Variable <array_1d<double,3> > >)
	.def("PrintVectorVariables", PrintVariablesName<Variable <Vector> >)
	.def("PrintMatrixVariables", PrintVariablesName<Variable <Matrix> >)
	.def("PrintStringVariables", PrintVariablesName<Variable <std::string> >)
	.def("PrintFlagsVariables", PrintVariablesName<Variable <Flags> >)
	.def("PrintVariableComponentVariables", PrintVariablesName<VariableComponent< VectorComponentAdaptor< array_1d<double,3> > > >)
	.def("GetAllVariableNames", GetVariableNames<VariableData>)
	.def("GetBoolVariableNames", GetVariableNames<Variable<bool> >)
	.def("GetIntVariableNames", GetVariableNames<Variable<int> >)
	.def("GetUnsignedIntVariableNames", GetVariableNames<Variable<int> >)
	.def("GetDoubleVariableNames", GetVariableNames<Variable <double> >)
	.def("GetArrayVariableNames", GetVariableNames<Variable <array_1d<double,3> > >)
	.def("GetVectorVariableNames", GetVariableNames<Variable <Vector> >)
	.def("GetMatrixVariableNames", GetVariableNames<Variable <Matrix> >)
	.def("GetStringVariableNames", GetVariableNames<Variable <std::string> >)
	.def("GetFlagsVariableNames", GetVariableNames<Variable <Flags> >)
	.def("GetVariableComponentVariableNames", GetVariableNames<VariableComponent< VectorComponentAdaptor< array_1d<double,3> > > >)
    .def(self_ns::str(self))
    ;
}

}  // namespace Python.

} // Namespace Kratos

