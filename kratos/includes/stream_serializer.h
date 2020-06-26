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
//

#if !defined(KRATOS_STREAM_SERIALIZER_H_INCLUDED )
#define  KRATOS_STREAM_SERIALIZER_H_INCLUDED

// System includes
#include <string>
#include <cstring>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/serializer.h"

//This class provides a simpler interface for serialization to a string instead of to a file
// Note that you may not override any load or save method of the Serializer. They are not virtual.
namespace Kratos
{
    class StreamSerializer : public Serializer 
    {
        public:
            KRATOS_CLASS_POINTER_DEFINITION(StreamSerializer); 

            ///this constructor simply wraps the standard Serializer and defines output to basic_iostream
            ///@param rTrace type of serialization to be employed
            explicit StreamSerializer(TraceType const& rTrace=SERIALIZER_NO_TRACE)
                : Serializer(new std::stringstream(std::ios::binary|std::ios::in|std::ios::out),rTrace)
            {}

            //this constructor generates a standard Serializer AND fills the buffer with the data contained in "data"
            ///@param data a string contained the data to be used in filling the buffer
            ///@param rTrace type of serialization to be employed
            StreamSerializer(const std::string& data,TraceType const& rTrace=SERIALIZER_NO_TRACE)
                : Serializer(new std::stringstream(std::ios::binary|std::ios::in|std::ios::out),rTrace)
            {
                *(this->pGetBuffer()) << data << std::endl;  
            }

            //get a string representation of the serialized data
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
}
#endif // KRATOS_STREAM_SERIALIZER_H_INCLUDED  defined
