//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
//kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/kratos_parameters.h"

#include "add_kratos_parameters_to_python.h"

namespace Kratos {

namespace Python {


pybind11::list items(Parameters const& self)
{
    pybind11::list t;
    for(Parameters::const_iterator it=self.begin(); it!=self.end(); ++it)
        t.append( std::make_tuple(it.name(), *it) );
    return t;
}

pybind11::list keys(Parameters const& self)
{
    pybind11::list t;
    for(Parameters::const_iterator it=self.begin(); it!=self.end(); ++it)
        t.append(it.name());
    return t;
}

pybind11::list values(Parameters const& self)
{
    pybind11::list t;
    for(Parameters::const_iterator it=self.begin(); it!=self.end(); ++it)
        t.append(*it);
    return t;
}

template< class T >
void Append(Parameters &rParameters, const T& obj) {
    rParameters.Append(obj);
}


void  AddKratosParametersToPython(pybind11::module& m)
{
    namespace py = pybind11;



    py::class_<Parameters, Parameters::Pointer >(m,"Parameters")
    .def(py::init<>())
    .def(py::init<const std::string&>())
    .def(py::init<Parameters const&>())
    .def("WriteJsonString", &Parameters::WriteJsonString)
    .def("PrettyPrintJsonString", &Parameters::PrettyPrintJsonString)
    .def("Has", &Parameters::Has)
    .def("Clone", &Parameters::Clone)
    .def("AddValue", &Parameters::AddValue)
    .def("AddEmptyValue", &Parameters::AddEmptyValue)
    .def("RemoveValue", &Parameters::RemoveValue)
    .def("ValidateAndAssignDefaults",&Parameters::ValidateAndAssignDefaults)
    .def("RecursivelyValidateAndAssignDefaults",&Parameters::RecursivelyValidateAndAssignDefaults)
    .def("IsEquivalentTo",&Parameters::IsEquivalentTo)
    .def("HasSameKeysAndTypeOfValuesAs",&Parameters::HasSameKeysAndTypeOfValuesAs)
    //.def("GetValue", &Parameters::GetValue) //Do not export this method. users shall adopt the operator [] syntax
    .def("IsNull", &Parameters::IsNull)
    .def("IsNumber", &Parameters::IsNumber)
    .def("IsDouble", &Parameters::IsDouble)
    .def("IsInt", &Parameters::IsInt)
    .def("IsBool", &Parameters::IsBool)
    .def("IsString", &Parameters::IsString)
    .def("IsArray", &Parameters::IsArray)
    .def("IsVector", &Parameters::IsVector)
    .def("IsMatrix", &Parameters::IsMatrix)
    .def("IsSubParameter", &Parameters::IsSubParameter)
    .def("GetDouble", &Parameters::GetDouble)
    .def("GetInt", &Parameters::GetInt)
    .def("GetBool", &Parameters::GetBool)
    .def("GetString", &Parameters::GetString)
    .def("GetVector", &Parameters::GetVector)
    .def("GetMatrix", &Parameters::GetMatrix)
    .def("SetDouble", &Parameters::SetDouble)
    .def("SetInt", &Parameters::SetInt)
    .def("SetBool", &Parameters::SetBool)
    .def("SetString", &Parameters::SetString)
    .def("SetVector", &Parameters::SetVector)
    .def("SetMatrix", &Parameters::SetMatrix)
    .def("size", &Parameters::size)
    //.def("GetArrayItem", &Parameters::GetArrayItem) //Do not export this method. users shall adopt the operator [] syntax
    .def("__setitem__", &Parameters::SetValue)
    .def("__getitem__", &Parameters::GetValue)
    .def("__setitem__", &Parameters::SetArrayItem)
    .def("__getitem__", &Parameters::GetArrayItem)
    .def("__iter__", [](Parameters& self){ return py::make_iterator(self.begin(), self.end()); } , py::keep_alive<0,1>())
    .def("items", &items )
    .def("keys", &keys )
    .def("values", &values )
    .def("__str__", PrintObject<Parameters>)
    .def("AddEmptyList", &Parameters::AddEmptyArray)
    .def("Append", Append<int>) // created due to ambiguous overload int/bool...
    .def("Append", Append<bool>) // created due to ambiguous overload int/bool...
    .def("Append", Append<double>) // created due to ambiguous overload int/bool...
    .def("Append", Append<Vector>) // created due to ambiguous overload int/bool...
    .def("Append", Append<Matrix>) // created due to ambiguous overload int/bool...
    .def("Append", Append<std::string>) // created due to ambiguous overload int/bool...
    .def("Append", Append<Parameters>) // created due to ambiguous overload int/bool...
    ;

}


} // namespace Python.

} // Namespace Kratos
