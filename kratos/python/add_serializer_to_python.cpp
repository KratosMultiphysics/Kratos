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
#include "containers/model.h"

namespace Kratos
{

        //this object is needed because data may contain binary data so the treatment of strings 
    //in python makes impossible to use them (reasons related to locale)
    class StreamSerializer : public Serializer 
    {
        public:
            KRATOS_CLASS_POINTER_DEFINITION(StreamSerializer); 

            StreamSerializer(TraceType const& rTrace=SERIALIZER_NO_TRACE)
                : Serializer(rTrace)
            {}

            StreamSerializer(const std::string& data,TraceType const& rTrace=SERIALIZER_NO_TRACE)
                : Serializer(rTrace)
            {
                *(this->pGetBuffer()) << data << std::endl;  
            }

            std::string GetStringRepresentation() {
                return ((std::stringstream*)(this->pGetBuffer()))->str();
             }

            virtual ~StreamSerializer(){}

        private:

            /// Assignment operator.
            StreamSerializer& operator=(StreamSerializer const& rOther) = delete;

            /// Copy constructor.
            StreamSerializer(StreamSerializer const& rOther) = delete;
    };
    
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
    .def("Load",SerializerLoad<Model>)
    .def("Save",SerializerSave<Model>)
    .def("Print", SerializerPrint)
    ;
    
    py::class_<StreamSerializer, StreamSerializer::Pointer, Serializer >(m,"StreamSerializer")
    .def(py::init<>())
    .def(py::init<std::string const&>())
    .def(py::init<Serializer::TraceType>())
    .def(py::init<std::string const&, Serializer::TraceType>())
    .def("__getstate__", [](StreamSerializer &self) { //METHOD NEEDED FOR PICKLE
        /* Return a tuple that fully encodes the state of the object */
        return py::make_tuple(py::bytes(self.GetStringRepresentation()));
    })
    .def("__setstate__", [](StreamSerializer &self, py::tuple t) { //METHOD NEEDED FOR PICKLE
        if (t.size() != 1)
            throw std::runtime_error("Invalid state!");

        /* Invoke the in-place constructor. Note that this is needed even
           when the object just has a trivial default constructor */
        new (&self) StreamSerializer(t[0].cast<std::string>());
    })
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

