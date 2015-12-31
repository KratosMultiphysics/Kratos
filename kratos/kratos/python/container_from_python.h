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



#if !defined(KRATOS_CONTAINER_FROM_PYTHON_H_INCLUDED )
#define KRATOS_CONTAINER_FROM_PYTHON_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"
#include "boost/python.hpp"
#include "boost/python/suite/indexing/vector_indexing_suite.hpp"
#include "boost/python/handle.hpp"

// Project includes
#include "includes/define.h"


namespace Kratos
{

namespace Python
{

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
/** Basically this class is a modified copy from scitbx/include/scitbx/boost_python/container_conversions.h.
    It's changed to work with index_suit method instead of scitbx conversion policies.
*/
template<class TContainerType>
class ContainerFromPython
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ContainerFromPython
    KRATOS_CLASS_POINTER_DEFINITION(ContainerFromPython);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ContainerFromPython()
    {
        boost::python::converter::registry::push_back(
            &convertible,
            &construct,
            boost::python::type_id<TContainerType>()
        );
    }

    /// Destructor.
    virtual ~ContainerFromPython() {}


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
        using boost::python::handle; // works around intel 9.1
        {
            // Restriction to list, tuple, iter, xrange until
            // Boost.Python overload resolution is enhanced.
            if (!(   PyList_Check(obj_ptr)
                     || PyTuple_Check(obj_ptr)
                     || PyIter_Check(obj_ptr)
                     || PyRange_Check(obj_ptr)))
            {
                return 0;
            }
        }

        boost::python::handle<> obj_iter( allow_null(PyObject_GetIter(obj_ptr) ) );

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
            /*           boost::type<TContainerType>(), obj_size)) return 0; */
            bool is_range = PyRange_Check(obj_ptr);
            std::size_t i=0;
            for(;; i++)
            {

                boost::python::handle<> py_elem_hdl(allow_null(PyIter_Next(obj_iter.get())));

                if (PyErr_Occurred())
                {
                    PyErr_Clear();
                    return 0;
                }
                if (!py_elem_hdl.get()) break; // end of iteration
                boost::python::object py_elem_obj(py_elem_hdl);
                boost::python::extract<typename TContainerType::value_type> elem_proxy(py_elem_obj);
                if (!elem_proxy.check()) return 0;
                if (is_range) break; // in a range all elements are of the same type
            }
            if (!is_range) assert(i == obj_size);
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
        boost::python::handle<> obj_iter(PyObject_GetIter(obj_ptr));
        void* storage = (
                            (rvalue_from_python_storage<TContainerType>*)
                            data)->storage.bytes;
        new (storage) TContainerType();
        data->convertible = storage;
        TContainerType& result = *((TContainerType*)storage);
        std::size_t i=0;
        for(;; i++)
        {

            boost::python::handle<> py_elem_hdl(allow_null(PyIter_Next(obj_iter.get())));


            if (PyErr_Occurred()) throw_error_already_set();
            if (!py_elem_hdl.get()) break; // end of iteration
            boost::python::object py_elem_obj(py_elem_hdl);
            boost::python::extract<typename TContainerType::value_type> elem_proxy(py_elem_obj);
            result.push_back(elem_proxy());
            /*         TConversionPolicyType::append(result, elem_proxy()); */
        }
        /*       TConversionPolicyType::assert_size(boost::type<TContainerType>(), i); */
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

}; // Class ContainerFromPython

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_CONTAINER_FROM_PYTHON_H_INCLUDED  defined 


