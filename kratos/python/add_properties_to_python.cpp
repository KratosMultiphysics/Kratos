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


template< class TContainerType, class TVariableType >
bool HasHelperFunction_Element(TContainerType& el, const TVariableType& rVar)
{
    return el.Has(rVar);
}


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
    .def("__setitem__", SetValueHelperFunction1< Properties, Variable< array_1d<double, 6> > >)
    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< array_1d<double, 6> > >)
    .def("Has", HasHelperFunction_Element< Properties, Variable< array_1d<double, 6> > >)
    .def("SetValue", SetValueHelperFunction1< Properties, Variable< array_1d<double, 6> > >)
    .def("GetValue", GetValueHelperFunction1< Properties, Variable< array_1d<double, 6> > >)

    .def("__setitem__", SetValueHelperFunction1< Properties, Variable< array_1d<double, 3> > >)
    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< array_1d<double, 3> > >)
    .def("Has", HasHelperFunction_Element< Properties, Variable< array_1d<double, 3> > >)
    .def("SetValue", SetValueHelperFunction1< Properties, Variable< array_1d<double, 3> > >)
    .def("GetValue", GetValueHelperFunction1< Properties, Variable< array_1d<double, 3> > >)

    .def("__setitem__", SetValueHelperFunction1< Properties, Variable< Vector > >)
    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< Vector > >)
    .def("Has", HasHelperFunction_Element< Properties, Variable< Vector > >)
    .def("SetValue", SetValueHelperFunction1< Properties, Variable< Vector > >)
    .def("GetValue", GetValueHelperFunction1< Properties, Variable< Vector > >)

    .def("__setitem__", SetValueHelperFunction1< Properties, Variable< Matrix > >)
    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< Matrix > >)
    .def("Has", HasHelperFunction_Element< Properties, Variable< Matrix > >)
    .def("SetValue", SetValueHelperFunction1< Properties, Variable< Matrix > >)
    .def("GetValue", GetValueHelperFunction1< Properties, Variable< Matrix > >)

    .def("__setitem__", SetValueHelperFunction1< Properties, Variable< std::string > >)
    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< std::string > >)
    .def("Has", HasHelperFunction_Element< Properties, Variable< std::string > >)
    .def("SetValue", SetValueHelperFunction1< Properties, Variable< std::string > >)
    .def("GetValue", GetValueHelperFunction1< Properties, Variable< std::string > >)

    .def("__setitem__", SetValueHelperFunction1< Properties, Variable< bool > >)
    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< bool > >)
    .def("Has", HasHelperFunction_Element< Properties, Variable< bool > >)
    .def("SetValue", SetValueHelperFunction1< Properties, Variable< bool > >)
    .def("GetValue", GetValueHelperFunction1< Properties, Variable< bool > >)

    .def("__setitem__", SetValueHelperFunction1< Properties, Variable< int > >)
    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< int > >)
    .def("Has", HasHelperFunction_Element< Properties, Variable< int > >)
    .def("SetValue", SetValueHelperFunction1< Properties, Variable< int > >)
    .def("GetValue", GetValueHelperFunction1< Properties, Variable< int > >)

    .def("__setitem__", SetValueHelperFunction1< Properties, Variable< double > >)
    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< double > >)
    .def("Has", HasHelperFunction_Element< Properties, Variable< double > >)
    .def("SetValue", SetValueHelperFunction1< Properties, Variable< double > >)
    .def("GetValue", GetValueHelperFunction1< Properties, Variable< double > >)

    .def("__setitem__", SetValueHelperFunction1< Properties, Variable< ConstitutiveLawBaseType::Pointer > >)
    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< ConstitutiveLawBaseType::Pointer > >)
    .def("Has", HasHelperFunction_Element< Properties, Variable< ConstitutiveLawBaseType::Pointer > >)
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

    .def("HasVariables", &Properties::HasVariables)
    .def("HasTables", &Properties::HasTables)
    .def("IsEmpty", &Properties::IsEmpty)

	.def(self_ns::str(self))
    ;

    PointerVectorSetPythonInterface<MeshType::PropertiesContainerType>::CreateInterface("PropertiesArray")
    ;
}
}  // namespace Python.
} // Namespace Kratos
