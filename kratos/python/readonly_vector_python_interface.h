// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



#if !defined(KRATOS_READONLY_VECTOR_PYTHON_INTERFACE_H_INCLUDED )
#define KRATOS_READONLY_VECTOR_PYTHON_INTERFACE_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>


// Project includes
#include "includes/define.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// A sequence converter from python.
/**
*/
template<class TContainerType>
class ReadonlyVectorPythonInterface : public vector_indexing_suite<TContainerType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ReadonlyVectorPythonInterface
    KRATOS_CLASS_POINTER_DEFINITION(ReadonlyVectorPythonInterface);

    typedef typename TContainerType::value_type data_type;
    typedef typename TContainerType::value_type key_type;
    typedef typename TContainerType::size_type index_type;
    typedef typename TContainerType::size_type size_type;
    typedef typename TContainerType::difference_type difference_type;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ReadonlyVectorPythonInterface()
    {
    }

    /// Destructor.
    virtual ~ReadonlyVectorPythonInterface() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    template <class Class>
    static void
    extension_def(Class& cl)
    {
    }

    static
    data_type const& get_item(TContainerType& container, index_type i)
    {
        return container[i];
    }

    static void
    set_item(TContainerType& container, index_type i, data_type const& v)
    {
    }

    static object
    get_slice(TContainerType& container, index_type from, index_type to)
    {
        index_type size = to - from;
        vector<data_type> result(size);
        for(index_type i = 0 ; i < size ; i++)
            result(i) = container(i + from);
        return object(result);
    }


    static void
    set_slice(TContainerType& container, index_type from,
              index_type to, data_type const& v)
    {
    }

    template <class Iter>
    static void
    set_slice(TContainerType& container, index_type from,
              index_type to, Iter first, Iter last)
    {
    }

    static void
    delete_item(TContainerType& container, index_type i)
    {
    }

    static void
    delete_slice(TContainerType& container, index_type from, index_type to)
    {
    }


    static void
    append(TContainerType& container, data_type const& v)
    {
    }

    template <class Iter>
    static void
    extend(TContainerType& container, Iter first, Iter last)
    {
    }

    static class_<TContainerType> CreateInterface(std::string const& Name)
    {
        return class_<TContainerType>(Name.c_str())
               .def(init<TContainerType>())
// 	  .def("Resize", &TContainerType::resize)
               .def("Size", &TContainerType::size)
               .def(indexing_suite<TContainerType, ReadonlyVectorPythonInterface>())
               .def(self_ns::str(self))
               ;
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{



    ///@}
    ///@name Friends
    ///@{


    ///@}

}; // Class ReadonlyVectorPythonInterface

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_READONLY_VECTOR_PYTHON_INTERFACE_H_INCLUDED defined 


