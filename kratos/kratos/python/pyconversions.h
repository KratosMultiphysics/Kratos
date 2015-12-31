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


/* -*- mode:c++ -*- */

/** @file

Boost.Python conversion functions from STL containers to Python
sequences and vice versa.

$Id: pyconversions.h,v 1.2 2007-03-06 10:30:34 rrossi Exp $

*/

#ifndef PY_CONVERSIONS_H
#define PY_CONVERSIONS_H

#include <vector>

#include <boost/python.hpp>

using namespace boost::python;

/** A wrapper of a conversion function to convert a STL vector to a
    Python tuple.  This class satisfies the requirements of the
    boost::python::to_python_converter conversion template argument.

    Copied from
    scitbx/include/scitbx/boost_python/container_conversions.h that is
    described in the <a
    href="http://www.boost.org/libs/python/doc/v2/faq.html">
    Boost.Python FAQ. </a>

    @author Ralf W. Grosse-Kunstleve <rwgk@yahoo.com> of
    <a href="http://www.lbl.gov/">Lawrence Berkeley National Laboratory</a>
*/
template < typename ContainerType >
struct to_tuple
{
    /** Creates and returns a Python @c tuple from the elements copied
    from a STL container. The @c ContainerType must be a vector,
    but may contain any type of object supported by the
    boost::python::object constructor. */
    static PyObject* convert (ContainerType const& c)
    {
        using boost::python::incref; // works around gcc 2.96 bug
        using boost::python::list; // dito
        list result;
        typename ContainerType::const_iterator i = c.begin();
        for( ; i != c.end(); ++i)
        {
            result.append(*i);
        }
        return incref(tuple(result).ptr());
    }
};

/** Converts an STL vector of T objects to Python tuple.

    Copied from
    scitbx/include/scitbx/boost_python/container_conversions.h that is
    described in the <a
    href="http://www.boost.org/libs/python/doc/v2/faq.html">
    Boost.Python FAQ. </a>

    @author Ralf W. Grosse-Kunstleve <rwgk@yahoo.com> of
    <a href="http://www.lbl.gov/">Lawrence Berkeley National Laboratory</a>
*/
template < typename T >
struct std_vector_to_tuple
{
    std_vector_to_tuple ()
    {
        to_python_converter < std::vector < T >,
                            to_tuple < std::vector < T > >  > ();
    }
};

// /** Interface to AxesType enumeration. */
// struct AxesOwner
// {
//   public:
//     AxesOwner(){};
//     typedef AxesType enum_type;
// };

// /** Conversion of AxesType enumeration. */
// struct EnumTypeConverters
//   : python::enum_as_int_converters <AxesOwner::enum_type>
// {
// };

/** Default operations on all containers for conversion from Python
    container to C++ one.

    Copied from
    scitbx/include/scitbx/boost_python/container_conversions.h that is
    described in the <a
    href="http://www.boost.org/libs/python/doc/v2/faq.html">
    Boost.Python FAQ. </a>

    @author Ralf W. Grosse-Kunstleve <rwgk@yahoo.com> of
    <a href="http://www.lbl.gov/">Lawrence Berkeley National Laboratory</a>
*/
struct default_policy
{
    static bool check_convertibility_per_element()
    {
        return false;
    }

    template <typename ContainerType>
    static bool check_size(boost::type<ContainerType>, std::size_t sz)
    {
        return true;
    }

    template <typename ContainerType>
    static void assert_size(boost::type<ContainerType>, std::size_t sz) {}

    template <typename ContainerType>
    static void reserve(ContainerType& a, std::size_t sz) {}
};

/** Operations on containers that have variable capacity for
    conversion from Python container to C++ one.

    Copied from
    scitbx/include/scitbx/boost_python/container_conversions.h that is
    described in the <a
    href="http://www.boost.org/libs/python/doc/v2/faq.html">
    Boost.Python FAQ. </a>

    @author Ralf W. Grosse-Kunstleve <rwgk@yahoo.com> of
    <a href="http://www.lbl.gov/">Lawrence Berkeley National Laboratory</a>
*/
struct variable_capacity_policy : default_policy
{
    template <typename ContainerType>
    static void reserve(ContainerType& a, std::size_t sz)
    {
        a.reserve(sz);
    }

    template <typename ContainerType, typename ValueType>
    static void set_value(ContainerType& a, std::size_t i, ValueType const& v)
    {
        assert(a.size() == i);
        a.push_back(v);
    }
};

/** Conversion of Python sequence to C++ container.

    Copied from
    scitbx/include/scitbx/boost_python/container_conversions.h that is
    described in the <a
    href="http://www.boost.org/libs/python/doc/v2/faq.html">
    Boost.Python FAQ. </a>

    @author Ralf W. Grosse-Kunstleve <rwgk@yahoo.com> of
    <a href="http://www.lbl.gov/">Lawrence Berkeley National Laboratory</a>
*/
template <typename ContainerType, typename ConversionPolicy>
struct from_python_sequence
{
    typedef typename ContainerType::value_type container_element_type;

    from_python_sequence()
    {
        boost::python::converter::registry::push_back(
            &convertible,
            &construct,
            boost::python::type_id<ContainerType>());
    }

    /** Appears to return @a obj_ptr if it is type of Python sequence
    that can be convertible to C++ container. */
    static void* convertible(PyObject* obj_ptr)
    {
        using namespace boost::python;
        using boost::python::allow_null; // works around gcc 2.96 bug
        {
            // Restriction to list, tuple, iter, xrange until
            // Boost.Python overload resolution is enhanced.
            //
            // add PySequence_Check() for numarray.
            //
            if (!(   PyList_Check(obj_ptr)
                     || PyTuple_Check(obj_ptr)
                     || PyIter_Check(obj_ptr)
                     || PyRange_Check(obj_ptr)
                     || PySequence_Check(obj_ptr) )) return 0;
        }
        handle<> obj_iter(allow_null(PyObject_GetIter(obj_ptr)));
        if (!obj_iter.get())   // must be convertible to an iterator
        {
            PyErr_Clear();
            return 0;
        }
        if (ConversionPolicy::check_convertibility_per_element())
        {
            int obj_size = PyObject_Length(obj_ptr);
            if (obj_size < 0)   // must be a measurable sequence
            {
                PyErr_Clear();
                return 0;
            }
            if (!ConversionPolicy::check_size(
                        boost::type<ContainerType>(), obj_size)) return 0;
            bool is_range = PyRange_Check(obj_ptr);
            //std::size_t i=0;
            int i = 0;
#ifndef _MSC_VER // because it causes c1001: internal compiler error
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
                extract<container_element_type> elem_proxy(py_elem_obj);
                if (!elem_proxy.check()) return 0;
                if (is_range) break; // in a range all elements are of the same type
            }
            if (!is_range) assert(i == obj_size );
#endif
        }
        return obj_ptr;
    }

    /** Constructs a C++ container from a Python sequence. */
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
                            (rvalue_from_python_storage<ContainerType>*)
                            data)->storage.bytes;
        new (storage) ContainerType();
        data->convertible = storage;
        ContainerType& result = *((ContainerType*)storage);
        std::size_t i=0;
        for(;; i++)
        {
            handle<> py_elem_hdl(allow_null(PyIter_Next(obj_iter.get())));
            if (PyErr_Occurred()) throw_error_already_set();
            if (!py_elem_hdl.get()) break; // end of iteration
            object py_elem_obj(py_elem_hdl);
            extract<container_element_type> elem_proxy(py_elem_obj);
            ConversionPolicy::set_value(result, i, elem_proxy());
        }
        ConversionPolicy::assert_size(boost::type<ContainerType>(), i);
    }

};

#endif // PY_CONVERSIONS_H
