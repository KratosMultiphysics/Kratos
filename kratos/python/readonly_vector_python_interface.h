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
//                   Riccardo Rossi
//




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


