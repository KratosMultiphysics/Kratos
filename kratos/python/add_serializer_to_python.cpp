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
//

// System includes

// External includes


// Project includes
#include "includes/define_python.h"
#include "includes/serializer.h"
#include "python/add_serializer_to_python.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;

template< class TObjectType >
void SerializerSave(Serializer& rSerializer, std::string const & rName, TObjectType& rObject)
{
    return rSerializer.save(rName, rObject);
}

template< class TObjectType >
void SerializerLoad(Serializer& rSerializer, std::string const & rName, TObjectType& rObject)
{
    return rSerializer.load(rName, rObject);
}

void SerializerPrint(Serializer& rSerializer)
{
    std::cout << "Serializer buffer:";
    std::cout << ((std::stringstream*)(rSerializer.pGetBuffer()))->str();
}

void  AddSerializerToPython(pybind11::module& m)
{
    py::class_<Serializer, Serializer::Pointer >(m,"Serializer")
    .def(py::init<>())
    .def(py::init<std::string const&>())
    .def(py::init<Serializer::TraceType>())
    .def(py::init<std::string const&, Serializer::TraceType>())
    .def("Load",SerializerLoad<ModelPart>)
    .def("Save",SerializerSave<ModelPart>)
    .def("Load",SerializerLoad<Parameters>)
    .def("Save",SerializerSave<Parameters>)
    .def("Print", SerializerPrint)
    ;

    py::enum_<Serializer::TraceType>(m,"SerializerTraceType")
    .value("SERIALIZER_NO_TRACE", Serializer::SERIALIZER_NO_TRACE)
    .value("SERIALIZER_TRACE_ERROR", Serializer::SERIALIZER_TRACE_ERROR)
    .value("SERIALIZER_TRACE_ALL", Serializer::SERIALIZER_TRACE_ALL)
    .export_values()
    ;

}

}  // namespace Python.

} // Namespace Kratos

