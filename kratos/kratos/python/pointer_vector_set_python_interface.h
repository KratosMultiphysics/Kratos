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



#if !defined(KRATOS_POINTER_VECTOR_SET_PYTHON_INTERFACE_H_INCLUDED )
#define KRATOS_POINTER_VECTOR_SET_PYTHON_INTERFACE_H_INCLUDED



// System includes


// External includes
/* # include <boost/python/suite/indexing/indexing_suite.hpp> */
/* # include <boost/python/suite/indexing/container_utils.hpp> */
/* # include <boost/python/iterator.hpp> */


// Project includes
#include "includes/define.h"
#include "python/pointer_indexing_python_interface.h"

namespace Kratos
{

namespace Python
{
using namespace boost::python;


/// A sequence converter from python.
/** Basically this class is a modified copy from scitbx/include/scitbx/boost_python/container_conversions.h.
*/
template<class TContainerType>
class PointerVectorSetPythonInterface
    : public pointer_indexing_suite< TContainerType
    , PointerVectorSetPythonInterface< TContainerType >
    , true
    , true
    , typename TContainerType::pointer //value_type
    , typename TContainerType::key_type
    , typename TContainerType::key_type>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PointerVectorSetPythonInterface
    KRATOS_CLASS_POINTER_DEFINITION(PointerVectorSetPythonInterface);

    typedef typename TContainerType::pointer /*value_type*/ data_type;
    typedef typename TContainerType::pointer pointer;
    typedef typename TContainerType::key_type key_type;
    typedef typename TContainerType::key_type index_type;
    typedef typename TContainerType::size_type size_type;
    typedef typename TContainerType::difference_type difference_type;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PointerVectorSetPythonInterface()
    {
    }

    /// Destructor.
    virtual ~PointerVectorSetPythonInterface() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    static void* convertible(PyObject* obj_ptr)
    {
        using namespace boost::python;
        using boost::python::allow_null; // works around gcc 2.96 bug
        {
            // Restriction to list, tuple, iter, xrange until
            // Boost.Python overload resolution is enhanced.
            if (!(   PyList_Check(obj_ptr)
                     || PyTuple_Check(obj_ptr)
                     || PyIter_Check(obj_ptr)
                     || PyRange_Check(obj_ptr))) return 0;
        }
        handle<> obj_iter(allow_null(PyObject_GetIter(obj_ptr)));
        if (!obj_iter.get())   // must be convertible to an iterator
        {
            PyErr_Clear();
            return 0;
        }
        if (true)
        {
            int obj_size = PyObject_Length(obj_ptr);
            if (obj_size < 0)   // must be a measurable sequence
            {
                PyErr_Clear();
                return 0;
            }
            /*         if (!TConversionPolicyType::check_size( */
            /*           boost::type<TTContainerTypeType>(), obj_size)) return 0; */
            bool is_range = PyRange_Check(obj_ptr);
            std::size_t i=0;
            for(;; i++)
            {
                handle<> py_elem_hdl(allow_null(PyIter_Next(obj_iter.get())));
                if (PyErr_Occurred())
                {
                    PyErr_Clear();
                    return 0;
                }
                if (!py_elem_hdl.get()) break; // end of iteration
                object py_elem_obj(py_elem_hdl);
                extract<pointer> elem_proxy(py_elem_obj);
                if (!elem_proxy.check()) return 0;
                if (is_range) break; // in a range all elements are of the same type
            }
            if (!is_range) assert(i == static_cast<std::size_t>(obj_size));
        }
        return obj_ptr;
    }

    static void construct(
        PyObject* obj_ptr,
        boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        using namespace boost::python;
        using boost::python::allow_null; // works around gcc 2.96 bug
        using boost::python::converter::rvalue_from_python_storage; // dito
        using boost::python::throw_error_already_set; // dito
        handle<> obj_iter(PyObject_GetIter(obj_ptr));
        void* storage = (
                            (rvalue_from_python_storage<TContainerType>*)
                            data)->storage.bytes;
        new (storage) TContainerType();
        data->convertible = storage;
        TContainerType& result = *((TContainerType*)storage);
        std::size_t i=0;
        for(;; i++)
        {
            handle<> py_elem_hdl(allow_null(PyIter_Next(obj_iter.get())));
            if (PyErr_Occurred()) throw_error_already_set();
            if (!py_elem_hdl.get()) break; // end of iteration
            object py_elem_obj(py_elem_hdl);
            extract<pointer> elem_proxy(py_elem_obj);
            append(result, elem_proxy());
        }
        /*       TConversionPolicyType::assert_size(boost::type<TContainerType>(), i); */
    }

    template <class Class>
    static void
    extension_def(Class& cl)
    {
        cl
            .def("append", &append)
        //	.def("extend", &base_extend)
            ;
    }

    static data_type get_item(TContainerType& container, index_type i)
    {
        return container(i);
    }

//static pointer get_item(TContainerType& container, index_type i)
//        {
//            return container(i);
//        }


    static object
    get_slice(TContainerType& container, index_type from, index_type to)
    {
        return object(TContainerType(container.find(from), container.find(to)));
    }


    //static void
    //	set_item(TContainerType& container, index_type i, data_type const& v)
    //{
    //	if(v == data_type())
    //		container.erase(i);
    //	else
    //		container[i] = v;
    //}

    static void
    set_item(TContainerType& container, index_type i, data_type /*const&*/ v)
    {
//			 container(i) = v;
        container.insert(container.begin(), v);
    }

    static void
    set_slice(TContainerType& container, index_type from,
              index_type to, data_type const& v)
    {
        if (from > to)
        {
            return;
        }
        else
        {
            container.erase(container.begin()+from, container.begin()+to);
            container.insert(container.begin()+from, v);
        }
    }
    template <class Iter>
    static void
    set_slice(TContainerType& container, index_type from,
              index_type to, Iter first, Iter last)
    {
        if (from > to)
        {
            container.insert(first, last);
        }
        else
        {
            container.erase(container.begin()+from, container.begin()+to);
            container.insert(first, last);
        }
    }

    static void
    delete_item(TContainerType& container, index_type i)
    {
        container.erase(i);
    }

    static void
    delete_slice(TContainerType& container, index_type from, index_type to)
    {
        container.erase(container.find(from), container.find(to));
    }

    static size_t
    size(TContainerType& container)
    {
        return container.size();
    }

    static void
    append(TContainerType& container, pointer v)
    {
        container.insert(container.end(),v);
    }

    //template <class Iter>
    //	static void
    //	extend(TContainerType& container, Iter first, Iter last)
    //{
    //	index_type i = container.size();
    //	TModifierType::Resize(container, i + std::distance(first,last));
    //	for( ; first != last ; first++)
    //		container.insert(i++, *first);
    //}

    static bool
    compare_index(TContainerType& container, index_type a, index_type b)
    {
        return typename TContainerType::key_compare()(a, b);
    }

    static index_type
    convert_index(TContainerType& container, PyObject* i_)
    {
        extract<key_type const&> i(i_);
        if (i.check())
        {
            return i();
        }
        else
        {
            extract<key_type> i1(i_);
            if (i1.check())
                return i1();
        }

        PyErr_SetString(PyExc_TypeError, "Invalid index type");
        throw_error_already_set();
        return index_type();
    }

    static class_<TContainerType, typename TContainerType::Pointer> CreateInterface(std::string const& Name)
    {
        boost::python::converter::registry::push_back(
            &convertible,
            &construct,
            boost::python::type_id<TContainerType>());

        return class_<TContainerType, typename TContainerType::Pointer> (Name.c_str())
               .def(init<TContainerType>())
               .def(PointerVectorSetPythonInterface<TContainerType>())
// 					.def(vector_indexing_suite<typename TContainerType::ContainerType, true>())
               .def("clear", &TContainerType::clear)
               .def("Size", &TContainerType::size)
               .def(self_ns::str(self))
               ;
    }

    static bool
    contains(TContainerType& container, key_type const& key)
    {
        return container.find(key) != container.end();
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

private:

}; // Class PointerVectorSetPythonInterface

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_POINTER_VECTOR_SET_PYTHON_INTERFACE_H_INCLUDED defined 


