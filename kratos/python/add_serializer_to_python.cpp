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
    class RawSerializedObject 
    {
        public:
            KRATOS_CLASS_POINTER_DEFINITION(RawSerializedObject); 

            RawSerializedObject(const std::string& data)
                : mdata(data)
            {}

            virtual ~RawSerializedObject(){}

            const std::string& GetData() const {
                return mdata;
                }

        private:
            std::string mdata;

            /// Assignment operator.
            RawSerializedObject& operator=(RawSerializedObject const& rOther) = delete;

            /// Copy constructor.
            RawSerializedObject(RawSerializedObject const& rOther) = delete;
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

    py::class_<RawSerializedObject, RawSerializedObject::Pointer >(m,"RawSerializedObject")
    .def("GetData",&RawSerializedObject::GetData)
    .def("__getstate__", [](const RawSerializedObject &self) { //METHOD NEEDED FOR PICKLE
        /* Return a tuple that fully encodes the state of the object */
        return py::make_tuple(py::bytes(self.GetData()));
    })
    .def("__setstate__", [](RawSerializedObject &self, py::tuple t) {
        if (t.size() != 1)
            throw std::runtime_error("Invalid state!");

        /* Invoke the in-place constructor. Note that this is needed even
           when the object just has a trivial default constructor */
        new (&self) RawSerializedObject(t[0].cast<std::string>());
    })
    ;


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
    .def("GetStringRepresentation", [](Serializer& self)
        {
            return Kratos::make_shared<RawSerializedObject>(((std::stringstream*)(self.pGetBuffer()))->str());
        })
    .def("FillBuffer", [](Serializer& self, RawSerializedObject& input_data)
        {
            (*self.pGetBuffer()) << input_data.GetData() << std::endl;  
        } 
    ) 
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

