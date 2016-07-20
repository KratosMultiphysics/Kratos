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



#if !defined(KRATOS_BOUNDED_VECTOR_PYTHON_INTERFACE_H_INCLUDED )
#define KRATOS_BOUNDED_VECTOR_PYTHON_INTERFACE_H_INCLUDED



// System includes
#include <cstddef>


// External includes


// Project includes
#include "includes/define.h"
#include "python/readonly_vector_python_interface.h"
#include "python/vector_vector_operator_python.h"
#include "python/vector_scalar_assignment_operator_python.h"
#include "python/bounded_vector_vector_assignment_operator_python.h"


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
*/
template<class TContainerType, std::size_t TSize>
class BoundedVectorPythonInterface : public ReadonlyVectorPythonInterface<TContainerType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BoundedVectorPythonInterface
    KRATOS_CLASS_POINTER_DEFINITION(BoundedVectorPythonInterface);

    typedef typename TContainerType::value_type data_type;
    typedef typename TContainerType::value_type key_type;
    typedef typename TContainerType::size_type index_type;
    typedef typename TContainerType::size_type size_type;
    typedef typename TContainerType::difference_type difference_type;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BoundedVectorPythonInterface()
    {
    }

    /// Destructor.
    virtual ~BoundedVectorPythonInterface() {}


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
            if (TSize < static_cast<std::size_t>(obj_size)) return 0;
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
                extract<typename TContainerType::value_type> elem_proxy(py_elem_obj);
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
            extract<typename TContainerType::value_type> elem_proxy(py_elem_obj);
            set_item(result, i, elem_proxy());
        }
        /*       TConversionPolicyType::assert_size(boost::type<TContainerType>(), i); */
    }

    template <class Class>
    static void
    extension_def(Class& cl)
    {
    }


    //static object
    //get_slice(TContainerType& container, index_type from, index_type to)
    //{
    //    return object((vector_range<TContainerType>(container, boost::numeric::ublas::range(from, to))));
    //}


    static void
    set_item(TContainerType& container, index_type i, data_type const& v)
    {
        if(v == data_type())
            container.erase_element(i);
        else
            container[i] = v;
    }

    static void
    set_slice(TContainerType& container, index_type from,
              index_type to, data_type const& v)
    {
        if(v == data_type())
            for(index_type i = from ; i < to ; i++)
                container.erase_element(i);
        else
            for(index_type i = from ; i < to ; i++)
                container.insert_element(i, v);
    }

    template <class Iter>
    static void
    set_slice(TContainerType& container, index_type from,
              index_type to, Iter first, Iter last)
    {
        index_type dist = std::distance(first,last) + from;
        if(to < dist)
        {
            dist -= to;
            index_type size = container.size();
            /* 	       TModifierType::Resize(container, size + dist); */
            MoveSlice(container, to + dist, to, size);
        }

        for( ; first != last ; first++)
            container(from++) = *first;

        delete_slice(container, from, to);

    }

    static void
    delete_item(TContainerType& container, index_type i)
    {
        delete_slice(container, i, i+1);
    }

    static void
    delete_slice(TContainerType& container, index_type from, index_type to)
    {
        if(from < to)
        {
            MoveSlice(container, from, to, container.size());
        }
    }


    static class_<TContainerType> CreateInterface(std::string const& Name)
    {
        boost::python::converter::registry::push_back(
            &convertible,
            &construct,
            boost::python::type_id<TContainerType>());

        return class_<TContainerType> (Name.c_str())
               .def(init<TContainerType>())
               .def(indexing_suite<TContainerType, BoundedVectorPythonInterface>())
               .def("Size", &TContainerType::size)
               .def(VectorVectorOperatorPython<TContainerType, TContainerType, TContainerType>())
               .def(VectorScalarAssignmentOperatorPython<TContainerType, double>())
               .def(BoundedVectorVectorAssignmentOperatorPython<TContainerType, zero_vector<double> >())
               .def(BoundedVectorVectorAssignmentOperatorPython<TContainerType, unit_vector<double> >())
               .def(BoundedVectorVectorAssignmentOperatorPython<TContainerType, scalar_vector<double> >())
               .def(BoundedVectorVectorAssignmentOperatorPython<TContainerType, vector<double> >())
               //	.def(BoundedVectorVectorAssignmentOperatorPython<TContainerType, mapped_vector<double> >())
               //	.def(BoundedVectorVectorAssignmentOperatorPython<TContainerType, compressed_vector<double> >())
               //	.def(BoundedVectorVectorAssignmentOperatorPython<TContainerType, coordinate_vector<double> >())
               .def(self_ns::str(self))
               ;
    }

    template <class TBaseType>
    static class_<TContainerType, bases<TBaseType> > CreateInterfaceWithBase(std::string const& Name, TBaseType Dummy)
    {
        boost::python::converter::registry::push_back(
            &convertible,
            &construct,
            boost::python::type_id<TContainerType>());

        return class_<TContainerType, bases<TBaseType> > (Name.c_str())
               .def(init<TContainerType>())
               //.def(indexing_suite<TContainerType, BoundedVectorPythonInterface>())
               .def("Size", &TContainerType::size)
               .def(VectorVectorOperatorPython<TContainerType, TContainerType, TContainerType>())
               .def(VectorScalarAssignmentOperatorPython<TContainerType, double>())
               .def(BoundedVectorVectorAssignmentOperatorPython<TContainerType, zero_vector<double> >())
               .def(BoundedVectorVectorAssignmentOperatorPython<TContainerType, unit_vector<double> >())
               .def(BoundedVectorVectorAssignmentOperatorPython<TContainerType, scalar_vector<double> >())
               .def(BoundedVectorVectorAssignmentOperatorPython<TContainerType, vector<double> >())
               //.def(BoundedVectorVectorAssignmentOperatorPython<TContainerType, mapped_vector<double> >())
               //	.def(BoundedVectorVectorAssignmentOperatorPython<TContainerType, compressed_vector<double> >())
               //	.def(BoundedVectorVectorAssignmentOperatorPython<TContainerType, coordinate_vector<double> >())
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

private:

    static TContainerType add_self(TContainerType& ThisContainer, TContainerType const& OtherContainer)
    {
        return ThisContainer + OtherContainer;
    }

    static TContainerType sub_self(TContainerType& ThisContainer, TContainerType const& OtherContainer)
    {
        return ThisContainer - OtherContainer;
    }

    static data_type mul_self(TContainerType& ThisContainer, TContainerType const& OtherContainer)
    {
        return inner_prod(ThisContainer, OtherContainer);
    }


    static void MoveSlice(TContainerType& ThisContainer, index_type Index, index_type From, index_type To)
    {
// 			    if(Index > From)
// 			      {
// 				ThisContainer.resize(ThisContainer.size() + Index - From, true);
// 				std::copy_backward(ThisContainer.begin() + From, ThisContainer.begin() + To, ThisContainer.begin() + Index + To - From);
// 			      }
// 			    else
// 			      {
// 				std::copy(ThisContainer.begin() + From, ThisContainer.begin() + To, ThisContainer.begin() + Index);
// 				ThisContainer.resize(ThisContainer.size() + Index - From, true);
// 			      }
    }


}; // Class BoundedVectorPythonInterface

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_BOUNDED_VECTOR_PYTHON_INTERFACE_H_INCLUDED defined 


