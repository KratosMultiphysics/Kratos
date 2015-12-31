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



// System includes

// External includes
#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>


// Project includes
#include "includes/define.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/properties.h"
#include "includes/constitutive_law.h"
//#include "python/add_mesh_to_python.h"
#include "python/pointer_vector_set_python_interface.h"
//#include "python/variable_indexing_python.h"

namespace Kratos
{
namespace Python
{
using namespace boost::python;

typedef Mesh<Node<3>, Properties, Element, Condition> MeshType;
typedef ConstitutiveLaw ConstitutiveLawBaseType;

template< class TContainerType, class TVariableType > void SetValueHelperFunction1(
    TContainerType& el,
    const TVariableType& rVar,
    const typename TVariableType::Type& Data)
{
    el.SetValue(rVar,Data);
}

template< class TContainerType, class TVariableType >
typename TVariableType::Type GetValueHelperFunction1( TContainerType& el,
        const TVariableType& rVar )
{
    return el.GetValue(rVar);
}



template< class TContainerType, class XVariableType, class YVariableType> void SetTableHelperFunction1(
    TContainerType& el,
    const XVariableType& XVar,
    const YVariableType& YVar,
	const typename Properties::TableType& Data)
{
    el.SetTable(XVar, YVar, Data);
}

template< class TContainerType, class XVariableType, class YVariableType>
typename Properties::TableType& GetTableHelperFunction1( TContainerType& el,
        const XVariableType& XVar,
    const YVariableType& YVar )
{
    return el.GetTable(XVar, YVar);
}

void  AddPropertiesToPython()
{
    class_<Properties, Properties::Pointer, bases<Properties::BaseType > >("Properties", init<int>())
    .def("__setitem__", SetValueHelperFunction1< Element, Variable< array_1d<double, 6> > >)
    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< array_1d<double, 6> > >)
    .def("SetValue", SetValueHelperFunction1< Properties, Variable< array_1d<double, 6> > >)
    .def("GetValue", GetValueHelperFunction1< Properties, Variable< array_1d<double, 6> > >)
	
    .def("__setitem__", SetValueHelperFunction1< Element, Variable< array_1d<double, 3> > >)
    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< array_1d<double, 3> > >)
    .def("SetValue", SetValueHelperFunction1< Properties, Variable< array_1d<double, 3> > >)
    .def("GetValue", GetValueHelperFunction1< Properties, Variable< array_1d<double, 3> > >)

    .def("__setitem__", SetValueHelperFunction1< Properties, Variable< Vector > >)
    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< Vector > >)
    .def("SetValue", SetValueHelperFunction1< Properties, Variable< Vector > >)
    .def("GetValue", GetValueHelperFunction1< Properties, Variable< Vector > >)

    .def("__setitem__", SetValueHelperFunction1< Properties, Variable< Matrix > >)
    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< Matrix > >)
    .def("SetValue", SetValueHelperFunction1< Properties, Variable< Matrix > >)
    .def("GetValue", GetValueHelperFunction1< Properties, Variable< Matrix > >)

    .def("__setitem__", SetValueHelperFunction1< Properties, Variable< std::string > >)
    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< std::string > >)
    .def("SetValue", SetValueHelperFunction1< Properties, Variable< std::string > >)
    .def("GetValue", GetValueHelperFunction1< Properties, Variable< std::string > >)

    .def("__setitem__", SetValueHelperFunction1< Properties, Variable< bool > >)
    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< bool > >)
    .def("SetValue", SetValueHelperFunction1< Properties, Variable< bool > >)
    .def("GetValue", GetValueHelperFunction1< Properties, Variable< bool > >)

    .def("__setitem__", SetValueHelperFunction1< Properties, Variable< int > >)
    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< int > >)
    .def("SetValue", SetValueHelperFunction1< Properties, Variable< int > >)
    .def("GetValue", GetValueHelperFunction1< Properties, Variable< int > >)

    .def("__setitem__", SetValueHelperFunction1< Properties, Variable< double > >)
    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< double > >)
    .def("SetValue", SetValueHelperFunction1< Properties, Variable< double > >)
    .def("GetValue", GetValueHelperFunction1< Properties, Variable< double > >)

    .def("__setitem__", SetValueHelperFunction1< Properties, Variable< ConstitutiveLawBaseType::Pointer > >)
    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< ConstitutiveLawBaseType::Pointer > >)
    .def("SetValue", SetValueHelperFunction1< Properties, Variable< ConstitutiveLawBaseType::Pointer > >)
    .def("GetValue", GetValueHelperFunction1< Properties, Variable< ConstitutiveLawBaseType::Pointer > >)

	.def("GetTable", GetTableHelperFunction1< Properties, Variable< double > , Variable<double> >, return_internal_reference<>())
    .def("GetTable", GetTableHelperFunction1< Properties, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > , Variable<double> >, return_internal_reference<>())
    .def("GetTable", GetTableHelperFunction1< Properties, Variable<double>, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >, return_internal_reference<>())
    .def("GetTable", GetTableHelperFunction1< Properties, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > , VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >, return_internal_reference<>())
    .def("SetTable", SetTableHelperFunction1< Properties, Variable< double > , Variable<double> >)
    .def("SetTable", SetTableHelperFunction1< Properties, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > , Variable<double> >)
    .def("SetTable", SetTableHelperFunction1< Properties, Variable<double>, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >)
    .def("SetTable", SetTableHelperFunction1< Properties, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > , VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >)

	.def(self_ns::str(self))
    ;

    PointerVectorSetPythonInterface<MeshType::PropertiesContainerType>::CreateInterface("PropertiesArray")
    ;
}
}  // namespace Python.
} // Namespace Kratos
