//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include <boost/python.hpp>

namespace Kratos
{

namespace Python
{

    
void  AddKratosParametersToPython()
{
    using namespace boost::python;

    class_<Parameters, Parameters::Pointer >("Parameters", init<std::string>()) //init<rapidjson::Value& >())
        .def(init<Parameters const&>())
        .def("WriteJsonString", &Parameters::WriteJsonString)
        .def("PrettyPrintJsonString", &Parameters::PrettyPrintJsonString)
        .def("Has", &Parameters::Has)
        .def("Clone", &Parameters::Clone)
        .def("ValidateAndAssignDefaults",&Parameters::ValidateAndAssignDefaults)
        //.def("GetValue", &Parameters::GetValue) //Do not export this method. users shall adopt the operator [] syntax
        .def("IsNumber", &Parameters::IsNumber)
        .def("IsDouble", &Parameters::IsDouble)
        .def("IsInt", &Parameters::IsInt)
        .def("IsBool", &Parameters::IsBool)
        .def("IsString", &Parameters::IsString)
        .def("IsArray", &Parameters::IsArray)
        .def("IsSubParameter", &Parameters::IsSubParameter)
        .def("GetDouble", &Parameters::GetDouble)
        .def("GetInt", &Parameters::GetInt)
        .def("GetBool", &Parameters::GetBool)
        .def("GetString", &Parameters::GetString)
        .def("SetDouble", &Parameters::SetDouble)
        .def("SetInt", &Parameters::SetInt)
        .def("SetBool", &Parameters::SetBool)
        .def("SetString", &Parameters::SetString)
        .def("size", &Parameters::size)
        //.def("GetArrayItem", &Parameters::GetArrayItem) //Do not export this method. users shall adopt the operator [] syntax
        .def("__setitem__", &Parameters::SetValue)
        .def("__getitem__", &Parameters::GetValue)
        .def("__setitem__", &Parameters::SetArrayItem)
        .def("__getitem__", &Parameters::GetArrayItem)
        ; 

//     class_<KratosParameters, KratosParameters::Pointer>("KratosParameters", init<std::string>())
//         .def(init<KratosParameters const&>())
//         .def("WriteJsonString", &KratosParameters::WriteJsonString)
//         .def("PrettyPrintJsonString", &KratosParameters::PrettyPrintJsonString)
//         .def("GetValue", &KratosParameters::GetValue)
//         .def("Has", &KratosParameters::Has)
//         .def("__setitem__", &KratosParameters::SetValue)
//         .def("__getitem__", &KratosParameters::GetValue)
//         .def("__setitem__", &KratosParameters::SetArrayItem)
//         .def("__getitem__", &KratosParameters::GetArrayItem)        
//         ;   

  


}

}  // namespace Python.

} // Namespace Kratos

