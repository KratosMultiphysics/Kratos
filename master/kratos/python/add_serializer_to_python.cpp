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
//

// System includes

// External includes


// Project includes
#include "includes/define_python.h"
#include "includes/serializer.h"
#include "includes/file_serializer.h"
#include "includes/stream_serializer.h"
#include "includes/mpi_serializer.h"
#include "python/add_serializer_to_python.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "containers/model.h"

namespace Kratos::Python
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

template< class TObjectType >
void SerializerLoadFromBeginning(Serializer& rSerializer, std::string const & rName, TObjectType& rObject)
{
    rSerializer.SetLoadState();
    return rSerializer.load(rName, rObject);
}

void SerializerPrint(Serializer& rSerializer)
{
    std::cout << "Serializer buffer:";
    std::cout << ((std::stringstream*)(rSerializer.pGetBuffer()))->str();
}

void  AddSerializerToPython(pybind11::module& m)
{

    auto serializer_py_interface = py::class_<Serializer, Serializer::Pointer >(m,"Serializer")
    .def(py::init([](const std::string& FileName) {
                    KRATOS_WARNING("DEPRECATION") << "Please use FileSerializer(FileName) instead of Serializer(FileName)" << std::endl;
                    return std::make_shared<FileSerializer>(FileName);
                }
            )
        )
    .def(py::init([](const std::string& FileName, Serializer::TraceType& rTraceType) {
                    KRATOS_WARNING("DEPRECATION") << "Please use FileSerializer(FileName,TraceType) instead of Serializer(FileName,TraceType)" << std::endl;
                    return std::make_shared<FileSerializer>(FileName,rTraceType);
                }
            )
        )
    .def("Load",SerializerLoad<ModelPart>)
    .def("LoadFromBeginning",SerializerLoadFromBeginning<ModelPart>)
    .def("Save",SerializerSave<ModelPart>)

    .def("Load",SerializerLoad<Parameters>)
    .def("LoadFromBeginning",SerializerLoadFromBeginning<Parameters>)
    .def("Save",SerializerSave<Parameters>)

    .def("Load",SerializerLoad<Model>)
    .def("LoadFromBeginning",SerializerLoadFromBeginning<Model>)
    .def("Save",SerializerSave<Model>)

    .def("Load",SerializerLoad<Flags>)
    .def("LoadFromBeginning",SerializerLoadFromBeginning<Flags>)
    .def("Save",SerializerSave<Flags>)

    .def("Set",   &Serializer::Set)
    .def("Print", SerializerPrint)
    ;

    serializer_py_interface.attr("MPI") = Serializer::MPI;
    serializer_py_interface.attr("SHALLOW_GLOBAL_POINTERS_SERIALIZATION") = Serializer::SHALLOW_GLOBAL_POINTERS_SERIALIZATION;

    py::class_<FileSerializer, FileSerializer::Pointer, Serializer >(m,"FileSerializer")
    .def(py::init<std::string const&>())
    .def(py::init<std::string const&, Serializer::TraceType>())
    ;

    py::class_<StreamSerializer, StreamSerializer::Pointer, Serializer >(m,"StreamSerializer")
    .def(py::init<>())
    .def(py::init<std::string const&>())
    .def(py::init<Serializer::TraceType>())
    .def(py::init<std::string const&, Serializer::TraceType>())
    .def(py::pickle(
            [](StreamSerializer &self) { // __getstate__
                return py::make_tuple(py::bytes(self.GetStringRepresentation()),self.GetTraceType());
            },
            [](py::tuple t) { // __setstate__, note: no `self` argument
                return Kratos::make_shared<StreamSerializer>(t[0].cast<std::string>(), t[1].cast<Serializer::TraceType>());
            }
        )
    )
    ;

    py::class_<MpiSerializer, MpiSerializer::Pointer, StreamSerializer >(m,"MpiSerializer")
    .def(py::init<>())
    .def(py::init<std::string const&>())
    .def(py::init<Serializer::TraceType>())
    .def(py::init<std::string const&, Serializer::TraceType>())
    .def(py::pickle(
            [](MpiSerializer &self) { // __getstate__
                return py::make_tuple(py::bytes(self.GetStringRepresentation()),self.GetTraceType());
            },
            [](py::tuple t) { // __setstate__, note: no `self` argument
                return Kratos::make_shared<MpiSerializer>(t[0].cast<std::string>(), t[1].cast<Serializer::TraceType>());
            }
        )
    )
    ;

    py::enum_<Serializer::TraceType>(m,"SerializerTraceType")
    .value("SERIALIZER_NO_TRACE", Serializer::SERIALIZER_NO_TRACE)
    .value("SERIALIZER_TRACE_ERROR", Serializer::SERIALIZER_TRACE_ERROR)
    .value("SERIALIZER_TRACE_ALL", Serializer::SERIALIZER_TRACE_ALL)
    .export_values()
    ;

}

}  // namespace Kratos::Python.

