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



#if !defined(KRATOS_MATRIX_PYTHON_INTERFACE_H_INCLUDED )
#define KRATOS_MATRIX_PYTHON_INTERFACE_H_INCLUDED



// System includes


// External includes
#include "boost/smart_ptr.hpp"
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include "python/readonly_matrix_python_interface.h"
#include "python/matrix_vector_operator_python.h"
#include "python/matrix_matrix_operator_python.h"


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
template<class TMatrixType, class TFunctorType = row_major>
class MatrixPythonInterface : public ReadonlyMatrixPythonInterface<TMatrixType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MatrixPythonInterface
    KRATOS_CLASS_POINTER_DEFINITION(MatrixPythonInterface);

    typedef typename TMatrixType::value_type data_type;
    typedef typename TMatrixType::value_type key_type;
    typedef typename TMatrixType::size_type index_type;
    typedef typename TMatrixType::size_type size_type;
    typedef typename TMatrixType::difference_type difference_type;
    typedef ReadonlyMatrixPythonInterface<TMatrixType> BaseType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Destructor.
    virtual ~MatrixPythonInterface() {}


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
                extract<vector<data_type> > elem_proxy(py_elem_obj);
                if (!elem_proxy.check()) return 0;
                if (is_range) break; // in a range all elements are of the same type
            }
            if (!is_range) assert(i == obj_size);
        }
        return obj_ptr;
    }

    static void construct1(
        PyObject* obj_ptr,
        boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        using namespace boost::python;
        using boost::python::allow_null; // works around gcc 2.96 bug
        using boost::python::converter::rvalue_from_python_storage; // dito
        using boost::python::throw_error_already_set; // dito
        handle<> obj_iter(PyObject_GetIter(obj_ptr));
        void* storage = (
                            (rvalue_from_python_storage<TMatrixType>*)
                            data)->storage.bytes;
        new (storage) TMatrixType();
        data->convertible = storage;
        TMatrixType& result = *((TMatrixType*)storage);
        int size2 = PyObject_Length(obj_ptr);
        std::size_t i=0;
        for(;; i++)
        {
            handle<> py_elem_hdl(allow_null(PyIter_Next(obj_iter.get())));
            if (PyErr_Occurred()) throw_error_already_set();
            if (!py_elem_hdl.get()) break; // end of iteration
            object py_elem_obj(py_elem_hdl);
            extract<vector<data_type> > elem_proxy(py_elem_obj);
            if(i == 0)
                result.resize(elem_proxy().size(), size2, false);
            set_row(result, i, elem_proxy());
        }
    }

    static void set_item(TMatrixType& ThisMatrix, tuple index, data_type Value)
    {
        unsigned int i = extract<index_type>(index[0]);
        unsigned int j = extract<index_type>(index[1]);
        if ((i >= ThisMatrix.size1()) || (j >= ThisMatrix.size2()))
        {
            PyErr_SetString(PyExc_IndexError, "index out of range");
            throw_error_already_set();
        }
        ThisMatrix(i,j) = Value;
    }

    static void set_row(TMatrixType& ThisMatrix, index_type i, vector<data_type> const& Value)
    {
        fill_row(ThisMatrix, i, Value, TFunctorType());
    }

// 			static class_<TMatrixType, boost::shared_ptr<TMatrixType> > CreateInterface(std::string const& Name)
    static class_<TMatrixType > CreateInterface(std::string const& Name)
    {
// 				boost::python::converter::registry::push_back(
// 					&convertible,
// 					&construct1,
// 					boost::python::type_id<TMatrixType>());

//  				return class_<TMatrixType, boost::shared_ptr<TMatrixType> >(Name.c_str())
        return class_<TMatrixType >(Name.c_str())
               .def(init<TMatrixType>())
               /* 					.def("Resize", &BaseType::resize1) */
               /* 					.def("Resize", &BaseType::resize2) */
               .def("Size1", &TMatrixType::size1)
               .def("Size2", &TMatrixType::size2)
               .def("__setitem__", &set_item)
               .def("__getitem__", &BaseType::get_item)
               .def(MatrixVectorOperatorPython<TMatrixType, vector<double> >())
// 					.def(MatrixVectorOperatorPython<TMatrixType, zero_vector<double>, vector<double> >())
// 					.def(MatrixVectorOperatorPython<TMatrixType, unit_vector<double>, vector<double> >())
// 					.def(MatrixVectorOperatorPython<TMatrixType, scalar_vector<double>, vector<double> >())
// 					.def(MatrixVectorOperatorPython<TMatrixType, mapped_vector<double> >())
// 					.def(MatrixVectorOperatorPython<TMatrixType, compressed_vector<double> >())
// 					.def(MatrixVectorOperatorPython<TMatrixType, coordinate_vector<double> >())
               .def(MatrixMatrixOperatorPython<TMatrixType, TMatrixType, TMatrixType>())
               .def(self_ns::str(self))

               ;

    }


private:

    static void fill_row(TMatrixType& ThisMatrix, index_type i, vector<data_type> const& Value, row_major const& dummy)
    {
        if (i >= ThisMatrix.size1())
        {
            PyErr_SetString(PyExc_IndexError, "index out of range");
            throw_error_already_set();
        }
        size_type size = std::min(Value.size(), ThisMatrix.size2());
        for(index_type j = 0 ; j < size ; j++)
            ThisMatrix(i,j) = Value[j];
    }

    static void fill_row(TMatrixType& ThisMatrix, index_type i, vector<data_type> const& Value, upper const& dummy)
    {
        if (i >= ThisMatrix.size1())
        {
            PyErr_SetString(PyExc_IndexError, "index out of range");
            throw_error_already_set();
        }
        size_type size = std::min(Value.size(), ThisMatrix.size2());
        for(index_type j = i ; j < size ; j++)
            ThisMatrix(i,j) = Value[j];
    }

    static void fill_row(TMatrixType& ThisMatrix, index_type i, vector<data_type> const& Value, lower const& dummy)
    {
        if (i >= ThisMatrix.size1())
        {
            PyErr_SetString(PyExc_IndexError, "index out of range");
            throw_error_already_set();
        }
        size_type size = std::min(Value.size(), i + 1);
        for(index_type j = 0 ; j < size ; j++)
            ThisMatrix(i,j) = Value[j];
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

}; // Class MatrixPythonInterface

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_MATRIX_PYTHON_INTERFACE_H_INCLUDED defined 


