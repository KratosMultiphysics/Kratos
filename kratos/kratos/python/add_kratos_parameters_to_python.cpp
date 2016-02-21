//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \
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

    class_<ParameterValue, ParameterValue::Pointer>("ParameterValue", init<rapidjson::Value& >())
        .def("Has", &ParameterValue::Has)
        .def("GetValue", &ParameterValue::GetValue)
        .def("IsNumber", &ParameterValue::IsNumber)
        .def("IsDouble", &ParameterValue::IsDouble)
        .def("IsInt", &ParameterValue::IsInt)
        .def("IsBool", &ParameterValue::IsBool)
        .def("IsString", &ParameterValue::IsString)
        .def("IsArray", &ParameterValue::IsArray)
        .def("IsSubParameter", &ParameterValue::IsSubParameter)
        .def("GetDouble", &ParameterValue::GetDouble)
        .def("GetInt", &ParameterValue::GetInt)
        .def("GetBool", &ParameterValue::GetBool)
        .def("GetString", &ParameterValue::GetString)
        .def("SetDouble", &ParameterValue::SetDouble)
        .def("SetInt", &ParameterValue::SetInt)
        .def("SetBool", &ParameterValue::SetBool)
        .def("SetString", &ParameterValue::SetString)
        .def("size", &ParameterValue::size)
        .def("GetArrayItem", &ParameterValue::GetArrayItem)
        ; 

    class_<KratosParameters, KratosParameters::Pointer>("KratosParameters", init<std::string>())
        .def("WriteJsonString", &KratosParameters::WriteJsonString)
        .def("PrettyPrintJsonString", &KratosParameters::PrettyPrintJsonString)
        .def("GetValue", &KratosParameters::GetValue)
        .def("Has", &KratosParameters::Has)
        ;   

  


}

}  // namespace Python.

} // Namespace Kratos

