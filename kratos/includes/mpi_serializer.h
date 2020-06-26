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

#if !defined(KRATOS_MPI_SERIALIZER_H_INCLUDED )
#define  KRATOS_MPI_SERIALIZER_H_INCLUDED

// System includes
#include <string>
#include <cstring>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/stream_serializer.h"

//This class provides a simpler interface for serialization to a string instead of to a file
// Note that you may not override any load or save method of the Serializer. They are not virtual.
namespace Kratos
{
    class MpiSerializer : public StreamSerializer 
    {
        public:
            KRATOS_CLASS_POINTER_DEFINITION(MpiSerializer); 

            ///this constructor simply wraps the standard Serializer and defines output to basic_iostream
            ///@param rTrace type of serialization to be employed
            explicit MpiSerializer(TraceType const& rTrace=SERIALIZER_NO_TRACE)
                : StreamSerializer(rTrace)
            {
                Set(Serializer::MPI);
                Set(Serializer::SHALLOW_GLOBAL_POINTERS_SERIALIZATION);
            }

            //this constructor generates a standard Serializer AND fills the buffer with the data contained in "data"
            ///@param data a string contained the data to be used in filling the buffer
            ///@param rTrace type of serialization to be employed
            MpiSerializer(const std::string& data,TraceType const& rTrace=SERIALIZER_NO_TRACE)
                : StreamSerializer(data,rTrace)
            {
                Set(Serializer::MPI);
                Set(Serializer::SHALLOW_GLOBAL_POINTERS_SERIALIZATION);                  
            }

            virtual ~MpiSerializer(){}

        private:

            /// Assignment operator.
            MpiSerializer& operator=(MpiSerializer const& rOther) = delete;

            /// Copy constructor.
            MpiSerializer(MpiSerializer const& rOther) = delete;
    };
}
#endif // KRATOS_MPI_SERIALIZER_H_INCLUDED  defined
