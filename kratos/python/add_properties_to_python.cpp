//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/mesh.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/properties.h"
#include "includes/constitutive_law.h"
#include "python/add_properties_to_python.h"
#include "python/containers_interface.h"

using AccessorBindType = std::unique_ptr<Kratos::Accessor>;

PYBIND11_MAKE_OPAQUE(AccessorBindType);

namespace Kratos::Python
{

namespace py = pybind11;

using MeshType = Mesh<Node, Properties, Element, Condition>;
using ConstitutiveLawBaseType = ConstitutiveLaw;
using IndexType = std::size_t;
using PropertiesContainerType = PointerVectorSet<Properties, IndexedObject>;
using GeometryType = Geometry<Node>;

template< class TContainerType, class TVariableType >
bool PropertiesHasHelperFunction(TContainerType& rProperties, const TVariableType& rVar)
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

PropertiesContainerType& GetSubPropertiesArray1(Properties& rProperties)
{
    return rProperties.GetSubProperties();
}

PropertiesContainerType const& GetSubPropertiesArray2(const Properties& rProperties)
{
    return rProperties.GetSubProperties();
}

bool HasSubProperties1(
    Properties& rProperties,
    IndexType Index
    )
{
    return rProperties.HasSubProperties(Index);
}

Properties::Pointer GetSubProperties1(
    Properties& rProperties,
    IndexType Index
    )
{
    return rProperties.pGetSubProperties(Index);
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

template< class TContainerType, class TVariableType >
void EraseHelperFunction1( TContainerType& rContainer,
        const TVariableType& rVar )
{
    rContainer.Erase(rVar);
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

template< class TContainerType, class XVariableType, class YVariableType>
bool HasTableHelperFunction1( TContainerType& rContainer,
        const XVariableType& XVar,
    const YVariableType& YVar )
{
    return rContainer.HasTable(XVar, YVar);
}

template<typename TVariableType>
void AddInterfaceToAccessorFold(pybind11::class_<Properties, Properties::Pointer, Properties::BaseType>& pyProperties) {

    pyProperties
    .def("GetAccessor", [](Properties &rProperties, Variable<TVariableType> &rVariable) { 
            auto accessor = &rProperties.pGetAccessor(rVariable);
            
            KRATOS_ERROR_IF(!*accessor) << "Trying to get a consumed or invalid Accessor." << std::endl;
                
            return accessor; 
        }, py::return_value_policy::reference_internal)
    .def("SetAccessor", [](Properties &rProperties, Variable<TVariableType> &rVariable, std::unique_ptr<Accessor>& rpAccessor) {
            KRATOS_ERROR_IF(!rpAccessor) << "Trying to set a consumed or invalid Accessor. Accessors are unique. Please create a different one." << std::endl;

            rProperties.SetAccessor(rVariable, std::move(rpAccessor));
        }, py::return_value_policy::reference_internal)
    .def("HasAccessor", [](Properties &rProperties, Variable<TVariableType> &rVariable) { 
            return rProperties.HasAccessor(rVariable);
        })
    ; 
}

template<typename... TArgs>
void AddInterfaceToAccessor(pybind11::class_<Properties, Properties::Pointer, Properties::BaseType>& pyProperties) {
    (AddInterfaceToAccessorFold<TArgs>(pyProperties), ...);
}

void  AddPropertiesToPython(pybind11::module& m)
{
    auto properties_module = py::class_<Properties, Properties::Pointer, Properties::BaseType>(m, "Properties")
        .def(py::init<Kratos::Properties::IndexType>())
        .def(py::init<const Properties &>())
        .def("__setitem__", SetValueHelperFunction1<Properties, Variable<array_1d<double, 6>>>)
        .def("__getitem__", GetValueHelperFunction1<Properties, Variable<array_1d<double, 6>>>)
        .def("Has", PropertiesHasHelperFunction<Properties, Variable<array_1d<double, 6>>>)
        .def("SetValue", SetValueHelperFunction1<Properties, Variable<array_1d<double, 6>>>)
        .def("GetValue", GetValueHelperFunction1<Properties, Variable<array_1d<double, 6>>>)
        .def("Erase", EraseHelperFunction1<Properties, Variable<array_1d<double, 6>>>)

        .def("__setitem__", SetValueHelperFunction1<Properties, Variable<array_1d<double, 3>>>)
        .def("__getitem__", GetValueHelperFunction1<Properties, Variable<array_1d<double, 3>>>)
        .def("Has", PropertiesHasHelperFunction<Properties, Variable<array_1d<double, 3>>>)
        .def("SetValue", SetValueHelperFunction1<Properties, Variable<array_1d<double, 3>>>)
        .def("GetValue", GetValueHelperFunction1<Properties, Variable<array_1d<double, 3>>>)
        .def("Erase", EraseHelperFunction1<Properties, Variable<array_1d<double, 3>>>)
        //     .def("SetValue", SetArrayValue)

        .def("__setitem__", SetValueHelperFunction1<Properties, Variable<Vector>>)
        .def("__getitem__", GetValueHelperFunction1<Properties, Variable<Vector>>)
        .def("Has", PropertiesHasHelperFunction<Properties, Variable<Vector>>)
        .def("SetValue", SetValueHelperFunction1<Properties, Variable<Vector>>)
        //     .def("SetValue", SetVectorValue)
        .def("GetValue", GetValueHelperFunction1<Properties, Variable<Vector>>)
        .def("Erase", EraseHelperFunction1<Properties, Variable<Vector>>)

        .def("__setitem__", SetValueHelperFunction1<Properties, Variable<Matrix>>)
        .def("__getitem__", GetValueHelperFunction1<Properties, Variable<Matrix>>)
        .def("Has", PropertiesHasHelperFunction<Properties, Variable<Matrix>>)
        .def("SetValue", SetValueHelperFunction1<Properties, Variable<Matrix>>)
        .def("GetValue", GetValueHelperFunction1<Properties, Variable<Matrix>>)
        .def("Erase", EraseHelperFunction1<Properties, Variable<Matrix>>)

        .def("__setitem__", SetValueHelperFunction1<Properties, Variable<std::string>>)
        .def("__getitem__", GetValueHelperFunction1<Properties, Variable<std::string>>)
        .def("Has", PropertiesHasHelperFunction<Properties, Variable<std::string>>)
        .def("SetValue", SetValueHelperFunction1<Properties, Variable<std::string>>)
        .def("GetValue", GetValueHelperFunction1<Properties, Variable<std::string>>)
        .def("Erase", EraseHelperFunction1<Properties, Variable<std::string>>)

        .def("__setitem__", SetValueHelperFunction1<Properties, Variable<bool>>)
        .def("__getitem__", GetValueHelperFunction1<Properties, Variable<bool>>)
        .def("Has", PropertiesHasHelperFunction<Properties, Variable<bool>>)
        .def("SetValue", SetValueHelperFunction1<Properties, Variable<bool>>)
        .def("GetValue", GetValueHelperFunction1<Properties, Variable<bool>>)
        .def("Erase", EraseHelperFunction1<Properties, Variable<bool>>)

        .def("__setitem__", SetValueHelperFunction1<Properties, Variable<int>>)
        .def("__getitem__", GetValueHelperFunction1<Properties, Variable<int>>)
        .def("Has", PropertiesHasHelperFunction<Properties, Variable<int>>)
        .def("SetValue", SetValueHelperFunction1<Properties, Variable<int>>)
        .def("GetValue", GetValueHelperFunction1<Properties, Variable<int>>)
        .def("Erase", EraseHelperFunction1<Properties, Variable<int>>)

        .def("__setitem__", SetValueHelperFunction1<Properties, Variable<double>>)
        .def("__getitem__", GetValueHelperFunction1<Properties, Variable<double>>)
        .def("Has", PropertiesHasHelperFunction<Properties, Variable<double>>)
        .def("SetValue", SetValueHelperFunction1<Properties, Variable<double>>)
        .def("GetValue", GetValueHelperFunction1<Properties, Variable<double>>)
        .def("Erase", EraseHelperFunction1<Properties, Variable<double>>)

        .def("__setitem__", SetValueHelperFunction1<Properties, Variable<ConstitutiveLawBaseType::Pointer>>)
        .def("__getitem__", GetValueHelperFunction1<Properties, Variable<ConstitutiveLawBaseType::Pointer>>)
        .def("Has", PropertiesHasHelperFunction<Properties, Variable<ConstitutiveLawBaseType::Pointer>>)
        .def("SetValue", SetValueHelperFunction1<Properties, Variable<ConstitutiveLawBaseType::Pointer>>)
        .def("GetValue", GetValueHelperFunction1<Properties, Variable<ConstitutiveLawBaseType::Pointer>>)
        .def("Erase", EraseHelperFunction1<Properties, Variable<ConstitutiveLawBaseType::Pointer>>)

        .def("GetTable", GetTableHelperFunction1<Properties, Variable<double>, Variable<double>>, py::return_value_policy::reference_internal)
        .def("SetTable", SetTableHelperFunction1<Properties, Variable<double>, Variable<double>>)
        .def("HasTable", HasTableHelperFunction1<Properties, Variable<double>, Variable<double>>)

        .def("HasVariables", &Properties::HasVariables)
        .def("HasTables", &Properties::HasTables)
        .def("IsEmpty", &Properties::IsEmpty)
        .def("NumberOfSubproperties", &Properties::NumberOfSubproperties)
        .def("AddSubProperties", &Properties::AddSubProperties)
        .def("HasSubProperties", HasSubProperties1)
        .def("GetSubProperties", GetSubProperties1)
        .def("GetSubProperties", GetSubPropertiesArray1)
        .def("GetSubProperties", GetSubPropertiesArray2)
        .def("SetSubProperties", &Properties::SetSubProperties)
        .def("__str__", PrintObject<Properties>)

        .def("GetValue", [](Properties &rProperties, const Variable<bool> &rVariable, const GeometryType &rGeometry, const Vector &rShapeFunctionVector, const ProcessInfo &rProcessInfo)
             { return rProperties.GetValue(rVariable, rGeometry, rShapeFunctionVector, rProcessInfo); })
        .def("GetValue", [](Properties &rProperties, const Variable<int> &rVariable, const GeometryType &rGeometry, const Vector &rShapeFunctionVector, const ProcessInfo &rProcessInfo)
             { return rProperties.GetValue(rVariable, rGeometry, rShapeFunctionVector, rProcessInfo); })
        .def("GetValue", [](Properties &rProperties, const Variable<double> &rVariable, const GeometryType &rGeometry, const Vector &rShapeFunctionVector, const ProcessInfo &rProcessInfo)
             { return rProperties.GetValue(rVariable, rGeometry, rShapeFunctionVector, rProcessInfo); })
        .def("GetValue", [](Properties &rProperties, const Variable<Vector> &rVariable, const GeometryType &rGeometry, const Vector &rShapeFunctionVector, const ProcessInfo &rProcessInfo)
             { return rProperties.GetValue(rVariable, rGeometry, rShapeFunctionVector, rProcessInfo); })
        .def("GetValue", [](Properties &rProperties, const Variable<Matrix> &rVariable, const GeometryType &rGeometry, const Vector &rShapeFunctionVector, const ProcessInfo &rProcessInfo)
             { return rProperties.GetValue(rVariable, rGeometry, rShapeFunctionVector, rProcessInfo); })
        .def("GetValue", [](Properties &rProperties, const Variable<array_1d<double, 3>> &rVariable, const GeometryType &rGeometry, const Vector &rShapeFunctionVector, const ProcessInfo &rProcessInfo)
             { return rProperties.GetValue(rVariable, rGeometry, rShapeFunctionVector, rProcessInfo); })
        .def("GetValue", [](Properties &rProperties, const Variable<array_1d<double, 6>> &rVariable, const GeometryType &rGeometry, const Vector &rShapeFunctionVector, const ProcessInfo &rProcessInfo)
             { return rProperties.GetValue(rVariable, rGeometry, rShapeFunctionVector, rProcessInfo); })
        .def("GetValue", [](Properties &rProperties, const Variable<array_1d<double, 9>> &rVariable, const GeometryType &rGeometry, const Vector &rShapeFunctionVector, const ProcessInfo &rProcessInfo)
             { return rProperties.GetValue(rVariable, rGeometry, rShapeFunctionVector, rProcessInfo); })
        .def("GetValue", [](Properties &rProperties, const Variable<array_1d<double, 4>> &rVariable, const GeometryType &rGeometry, const Vector &rShapeFunctionVector, const ProcessInfo &rProcessInfo)
             { return rProperties.GetValue(rVariable, rGeometry, rShapeFunctionVector, rProcessInfo); })
        .def("GetValue", [](Properties &rProperties, const Variable<std::string> &rVariable, const GeometryType &rGeometry, const Vector &rShapeFunctionVector, const ProcessInfo &rProcessInfo)
             { return rProperties.GetValue(rVariable, rGeometry, rShapeFunctionVector, rProcessInfo); })
        ;

        AddInterfaceToAccessor<
            bool,
            int,
            double,
            Vector,
            Matrix,
            array_1d<double, 3>,
            array_1d<double, 6>,
            array_1d<double, 9>,
            array_1d<double, 4>,
            std::string
        >(properties_module);

    PointerVectorSetPythonInterface<MeshType::PropertiesContainerType>().CreateInterface(m,"PropertiesArray");
}

} // Namespace Kratos::Python
