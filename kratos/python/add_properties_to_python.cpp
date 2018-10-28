//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//



// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/define_python.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/properties.h"
#include "includes/constitutive_law.h"
#include "python/add_mesh_to_python.h"
#include "python/containers_interface.h"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;

typedef Node<3> NodeType;
typedef Mesh<NodeType, Properties, Element, Condition> MeshType;
typedef ConstitutiveLaw ConstitutiveLawBaseType;
typedef std::size_t IndexType;
typedef PointerVectorSet<Properties, IndexedObject> PropertiesContainerType;

template< class TContainerType, class TVariableType >
bool HasHelperFunction_Element(TContainerType& rProperties, const TVariableType& rVar)
{
    return rProperties.Has(rVar);
}


template< class TContainerType, class TVariableType >
void SetValueHelperFunction1(
    TContainerType& rProperties,
    const TVariableType& rVar,
    const typename TVariableType::Type& rData)
{
    rProperties.SetValue(rVar,rData);
}

PropertiesContainerType& GetSubProperties1(Properties& rProperties)
{
    return rProperties.GetSubProperties();
}

PropertiesContainerType const& GetSubProperties2(const Properties& rProperties)
{
    return rProperties.GetSubProperties();
}

Properties::Pointer GetSubProperty1(
    Properties& rProperties,
    IndexType Index
    )
{
    return rProperties.pGetSubProperty(Index);
}

void SetArrayValue(
    Properties& rProperties,
    const Variable<array_1d<double,3>>& rVar,
    const std::vector<double>& rData)
{
    if(rData.size() != 3)
        KRATOS_ERROR << "attempting to construct an array<double,3> by passing a list with wrong size. Input size is " << rData.size() << std::endl;

    array_1d<double,3> tmp;
    for(unsigned int i=0;i<3; ++i)
        tmp[i] = rData[i];

    rProperties.SetValue(rVar,tmp);
}

void SetVectorValue(
    Properties& rProperties,
    const Variable<Vector>& rVar,
    const std::vector<double>& rData)
{
    Vector tmp(rData.size());
    for(unsigned int i=0;i<tmp.size(); ++i)
        tmp[i] = rData[i];

    rProperties.SetValue(rVar,tmp);
}

template< class TContainerType, class TVariableType >
typename TVariableType::Type GetValueHelperFunction1( TContainerType& rContainer,
        const TVariableType& rVar )
{
    return rContainer.GetValue(rVar);
}

template< class TContainerType, class XVariableType, class YVariableType>
void SetTableHelperFunction1(
    TContainerType& rContainer,
    const XVariableType& XVar,
    const YVariableType& YVar,
    const typename Properties::TableType& rData)
{
    rContainer.SetTable(XVar, YVar, rData);
}

template< class TContainerType, class XVariableType, class YVariableType>
typename Properties::TableType& GetTableHelperFunction1( TContainerType& rContainer,
        const XVariableType& XVar,
    const YVariableType& YVar )
{
    return rContainer.GetTable(XVar, YVar);
}

void  AddPropertiesToPython(pybind11::module& m)
{
    py::class_<Properties, Properties::Pointer, Properties::BaseType >(m,"Properties")
    .def(py::init<Kratos::Properties::IndexType>())
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
//     .def("SetValue", SetArrayValue)

    .def("__setitem__", SetValueHelperFunction1< Properties, Variable< Vector > >)
    .def("__getitem__", GetValueHelperFunction1< Properties, Variable< Vector > >)
    .def("Has", HasHelperFunction_Element< Properties, Variable< Vector > >)
    .def("SetValue", SetValueHelperFunction1< Properties, Variable< Vector > >)
//     .def("SetValue", SetVectorValue)
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

	  .def("GetTable", GetTableHelperFunction1< Properties, Variable< double > , Variable<double> >, py::return_value_policy::reference_internal)
    .def("GetTable", GetTableHelperFunction1< Properties, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > , Variable<double> >, py::return_value_policy::reference_internal)
    .def("GetTable", GetTableHelperFunction1< Properties, Variable<double>, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >, py::return_value_policy::reference_internal)
    .def("GetTable", GetTableHelperFunction1< Properties, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > , VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >, py::return_value_policy::reference_internal)
    .def("SetTable", SetTableHelperFunction1< Properties, Variable< double > , Variable<double> >)
    .def("SetTable", SetTableHelperFunction1< Properties, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > , Variable<double> >)
    .def("SetTable", SetTableHelperFunction1< Properties, Variable<double>, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >)
    .def("SetTable", SetTableHelperFunction1< Properties, VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > , VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >)
    .def("__repr__", &Properties::Info) //self_ns::str(self))
    .def("HasVariables", &Properties::HasVariables)
    .def("HasTables", &Properties::HasTables)
    .def("IsEmpty", &Properties::IsEmpty)
    .def("NumberOfSubproperties", &Properties::NumberOfSubproperties)
    .def("AddSubProperty", &Properties::AddSubProperty)
    .def("GetSubProperty", GetSubProperty1)
    .def("GetSubProperties", GetSubProperties1)
    .def("GetSubProperties", GetSubProperties2)
    .def("SetSubProperties", &Properties::SetSubProperties)
    .def("__str__", PrintObject<Properties>)
    ;

    PointerVectorSetPythonInterface<MeshType::PropertiesContainerType>().CreateInterface(m,"PropertiesArray");
}
}  // namespace Python.
} // Namespace Kratos
