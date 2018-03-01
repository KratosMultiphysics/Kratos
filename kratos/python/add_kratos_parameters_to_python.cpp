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
<<<<<<< HEAD
#include "add_kratos_parameters_to_python.h"
=======
>>>>>>> master

namespace Kratos {

namespace Python {

<<<<<<< HEAD

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

void  AddKratosParametersToPython(pybind11::module& m)
{
    using namespace pybind11;



    class_<Parameters, Parameters::Pointer >(m,"Parameters")
    .def(init<const std::string>()) //init<rapidjson::Value& >())
    .def(init<Parameters const&>())
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
    .def("__iter__", [](Parameters& self){ return make_iterator(self.begin(), self.end()); } , keep_alive<0,1>()) 
    .def("items", &items )
    .def("keys", &keys )
    .def("values", &values )
    .def("__repr__",&Parameters::Info)
    ;

=======
Parameters::iterator NonConstBegin(Parameters &el) { return el.begin(); }
Parameters::iterator NonConstEnd(Parameters &el) { return el.end(); }

boost::python::list items(Parameters const &self) {
  boost::python::list t;
  for (Parameters::const_iterator it = self.begin(); it != self.end(); ++it)
    t.append(boost::python::make_tuple(it.name(), *it));
  return t;
}

boost::python::list keys(Parameters const &self) {
  boost::python::list t;
  for (Parameters::const_iterator it = self.begin(); it != self.end(); ++it)
    t.append(it.name());
  return t;
}

boost::python::list values(Parameters const &self) {
  boost::python::list t;
  for (Parameters::const_iterator it = self.begin(); it != self.end(); ++it)
    t.append(*it);
  return t;
}

void Append(Parameters &rParameters, PyObject *type) {
  if (PyLong_CheckExact(type))
    rParameters.Append(boost::python::extract<int>(type));
  else if (PyBool_Check(type))
    rParameters.Append(boost::python::extract<bool>(type));
  else if (PyUnicode_Check(type))
    rParameters.Append(boost::python::extract<std::string>(type));
  else if (boost::python::arg_from_python<Vector>(type).convertible())
    rParameters.Append(boost::python::extract<Vector>(type));
  else if (boost::python::arg_from_python<Matrix>(type).convertible())
    rParameters.Append(boost::python::extract<Matrix>(type));
  else if (boost::python::arg_from_python<Parameters>(type).convertible())
    rParameters.Append(boost::python::extract<Parameters>(type));
  else if (boost::python::arg_from_python<double>(type).convertible())
    rParameters.Append(boost::python::extract<double>(type));
  else
    KRATOS_ERROR
        << " python object type not accepted to append in parameters Array "
        << std::endl;
}
>>>>>>> master

void AddKratosParametersToPython() {
  using namespace boost::python;

  class_<Parameters, Parameters::Pointer>(
      "Parameters", init<const std::string>()) // init<rapidjson::Value& >())
      .def(init<Parameters const &>())
      .def("WriteJsonString", &Parameters::WriteJsonString)
      .def("PrettyPrintJsonString", &Parameters::PrettyPrintJsonString)
      .def("Has", &Parameters::Has)
      .def("Clone", &Parameters::Clone)
      .def("AddValue", &Parameters::AddValue)
      .def("AddEmptyValue", &Parameters::AddEmptyValue)
      .def("RemoveValue", &Parameters::RemoveValue)
      .def("ValidateAndAssignDefaults", &Parameters::ValidateAndAssignDefaults)
      .def("RecursivelyValidateAndAssignDefaults",
           &Parameters::RecursivelyValidateAndAssignDefaults)
      .def("IsEquivalentTo", &Parameters::IsEquivalentTo)
      .def("HasSameKeysAndTypeOfValuesAs",
           &Parameters::HasSameKeysAndTypeOfValuesAs)
      //.def("GetValue", &Parameters::GetValue) //Do not export this method.
      //users shall adopt the operator [] syntax
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
      //.def("GetArrayItem", &Parameters::GetArrayItem) //Do not export this
      //method. users shall adopt the operator [] syntax
      .def("__setitem__", &Parameters::SetValue)
      .def("__getitem__", &Parameters::GetValue)
      .def("__setitem__", &Parameters::SetArrayItem)
      .def("__getitem__", &Parameters::GetArrayItem)
      .def("__iter__", boost::python::range(&NonConstBegin, &NonConstEnd))
      .def("items", &items)
      .def("keys", &keys)
      .def("values", &values)
      .def("AddEmptyList", &Parameters::AddEmptyArray)
      .def("Append", Append) // created due to ambiguous overload int/bool...
      .def(self_ns::str(self));
}

} // namespace Python.

} // Namespace Kratos
