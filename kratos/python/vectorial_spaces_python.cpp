/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:34 $
//   Revision:            $Revision: 1.2 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "python/vector_python_interface.h"
#include "python/vector_scalar_operator_python.h"

namespace Kratos
{
	
namespace Python
{
		
  using namespace boost::python;
	
  template<class TContainerType> 
  struct UblasVectorModifier
  {
    typedef typename TContainerType::size_type index_type;
    static void Resize(TContainerType& ThisContainer, typename TContainerType::size_type NewSize)
    {
      ThisContainer.resize(NewSize, true);
    }
    static void MoveSlice(TContainerType& ThisContainer, index_type Index, index_type From, index_type To)
    {
      if(Index > From)
	{
	  ThisContainer.resize(ThisContainer.size() + Index - From, true);
	  std::copy_backward(ThisContainer.begin() + From, ThisContainer.begin() + To, ThisContainer.begin() + Index + To - From);
	}
      else
	{
	  std::copy(ThisContainer.begin() + From, ThisContainer.begin() + To, ThisContainer.begin() + Index);
	  ThisContainer.resize(ThisContainer.size() + Index - From, true);
	}
    }
  };
	
  template<class TContainerType> 
  struct UblasSparseVectorModifier
  {
    typedef typename TContainerType::size_type index_type;
    typedef typename TContainerType::value_type data_type;

    static void Resize(TContainerType& ThisContainer, typename TContainerType::size_type NewSize)
    {
      index_type size = std::min(ThisContainer.size(), NewSize);
      std::vector<std::pair<index_type, data_type> > temp;
      data_type value;

      for(index_type i = 0 ; i < size ; i++)
	if((value = ThisContainer(i)) != data_type())
	    temp.push_back(std::pair<index_type, data_type>(i,value));

      ThisContainer.clear(); // There is no way to know which resize hold the data and which not. So better to make it certain! :-)
      ThisContainer.resize(NewSize, false);

      for(typename std::vector<std::pair<index_type, data_type> >::iterator j = temp.begin() ; j != temp.end() ; j++)
	ThisContainer.insert_element(j->first, j->second);

    }

    static void MoveSlice(TContainerType& ThisContainer, index_type Index, index_type From, index_type To)
    {
      index_type i;
      index_type size = std::min(Index, From);
      index_type new_size = ThisContainer.size() + Index - From;
      std::vector<std::pair<index_type, data_type> > temp;
      data_type value;

      for(i = 0 ; i < size ; i++)
	if((value = ThisContainer(i)) != data_type())
	    temp.push_back(std::pair<index_type, data_type>(i,value));

      for(i = From ; i < To ; i++)
	if((value = ThisContainer(i)) != data_type())
	    temp.push_back(std::pair<index_type, data_type>(i + Index - From,value));

      ThisContainer.clear(); // There is no way to know which resize hold the data and which not. So better to make it certain! :-)
      ThisContainer.resize(new_size, false);

      for(typename std::vector<std::pair<index_type, data_type> >::iterator j = temp.begin() ; j != temp.end() ; j++)
	ThisContainer.insert_element(j->first, j->second);
     }
  };
	
	
  void  AddVectorialSpacesToPython()
  {

    ReadonlyVectorPythonInterface<zero_vector<double> >::CreateInterface("ZeroVector")
      .def(init<zero_vector<double>::size_type>())
      .def(VectorScalarOperatorPython<zero_vector<double>, double, vector<double> >())
      .def(VectorVectorOperatorPython<zero_vector<double>, unit_vector<double>, vector<double> >())
      .def(VectorVectorOperatorPython<zero_vector<double>, scalar_vector<double>, vector<double> >())
      .def(VectorVectorOperatorPython<zero_vector<double>, vector<double>, vector<double> >())
      //      .def(VectorVectorOperatorPython<zero_vector<double>, mapped_vector<double>, mapped_vector<double> >())
      //.def(VectorVectorOperatorPython<zero_vector<double>, compressed_vector<double>, compressed_vector<double> >())
      //.def(VectorVectorOperatorPython<zero_vector<double>, coordinate_vector<double>, coordinate_vector<double> >())
       ;

    ReadonlyVectorPythonInterface<unit_vector<double> >::CreateInterface("UnitVector")
      .def(init<unit_vector<double>::size_type, vector<double>::size_type>())
      .def(VectorScalarOperatorPython<unit_vector<double>, double, vector<double> >())
      .def(VectorVectorOperatorPython<unit_vector<double>, zero_vector<double>, vector<double> >())
      .def(VectorVectorOperatorPython<unit_vector<double>, scalar_vector<double>, vector<double> >())
      .def(VectorVectorOperatorPython<unit_vector<double>, vector<double>, vector<double> >())
      //.def(VectorVectorOperatorPython<unit_vector<double>, mapped_vector<double>, mapped_vector<double> >())
      //.def(VectorVectorOperatorPython<unit_vector<double>, compressed_vector<double>, compressed_vector<double> >())
      //.def(VectorVectorOperatorPython<unit_vector<double>, coordinate_vector<double>, coordinate_vector<double> >())
       ;

    ReadonlyVectorPythonInterface<scalar_vector<double> >::CreateInterface("ScalarVector")
      .def(init<scalar_vector<double>::size_type, scalar_vector<double>::value_type>())
      .def(VectorScalarOperatorPython<scalar_vector<double>, double, vector<double> >())
      .def(VectorVectorOperatorPython<scalar_vector<double>, zero_vector<double>, vector<double> >())
      .def(VectorVectorOperatorPython<scalar_vector<double>, unit_vector<double>, vector<double> >())
      .def(VectorVectorOperatorPython<scalar_vector<double>, vector<double>, vector<double> >())
      //.def(VectorVectorOperatorPython<scalar_vector<double>, mapped_vector<double>, mapped_vector<double> >())
      //.def(VectorVectorOperatorPython<scalar_vector<double>, compressed_vector<double>, compressed_vector<double> >())
      //.def(VectorVectorOperatorPython<scalar_vector<double>, coordinate_vector<double>, coordinate_vector<double> >())
       ;

    VectorPythonInterface<vector<double>, UblasVectorModifier<vector<double> > >::CreateInterface("Vector")
      .def(init<vector<double>::size_type>())
      .def(init<vector_expression<vector<double> > >())
      .def(VectorScalarOperatorPython<vector<double>, double, vector<double> >())
      .def(VectorVectorOperatorPython<vector<double>, zero_vector<double>, vector<double> >())
      .def(VectorVectorOperatorPython<vector<double>, unit_vector<double>, vector<double> >())
      .def(VectorVectorOperatorPython<vector<double>, scalar_vector<double>, vector<double> >())
      //.def(VectorVectorOperatorPython<vector<double>, mapped_vector<double>, vector<double> >())
       ;


    //VectorPythonInterface<mapped_vector<double>, UblasSparseVectorModifier<mapped_vector<double> > >::CreateInterface("SparseVector")
    //.def(init<mapped_vector<double>::size_type>())
    //.def("NonZeros", &mapped_vector<double>::nnz)
    //.def(VectorScalarOperatorPython<mapped_vector<double>, double, mapped_vector<double> >())
    //.def(VectorVectorOperatorPython<mapped_vector<double>, zero_vector<double>, mapped_vector<double> >())
    //.def(VectorVectorOperatorPython<mapped_vector<double>, unit_vector<double>, mapped_vector<double> >())
    //.def(VectorVectorOperatorPython<mapped_vector<double>, scalar_vector<double>, mapped_vector<double> >())
    //.def(VectorVectorOperatorPython<mapped_vector<double>, vector<double>, vector<double> >())
    //;
		
    //VectorPythonInterface<compressed_vector<double>, UblasSparseVectorModifier<compressed_vector<double> > >::CreateInterface("CompressedVector")
    //.def(init<compressed_vector<double>::size_type>())
    //.def("NonZeros", &compressed_vector<double>::nnz)
    //.def(VectorScalarOperatorPython<compressed_vector<double>, double, compressed_vector<double> >())
    //.def(VectorVectorOperatorPython<compressed_vector<double>, zero_vector<double>, compressed_vector<double> >())
    //.def(VectorVectorOperatorPython<compressed_vector<double>, unit_vector<double>, compressed_vector<double> >())
    //.def(VectorVectorOperatorPython<compressed_vector<double>, scalar_vector<double>, compressed_vector<double> >())
    //.def(VectorVectorOperatorPython<compressed_vector<double>, vector<double>, vector<double> >())
    //.def(VectorVectorOperatorPython<compressed_vector<double>, mapped_vector<double>, compressed_vector<double> >())
    //.def(VectorVectorOperatorPython<compressed_vector<double>, coordinate_vector<double>, compressed_vector<double> >())
    //;
		
    //    VectorPythonInterface<coordinate_vector<double>, UblasSparseVectorModifier<coordinate_vector<double> > >::CreateInterface("CoordinateVector")
    //.def(init<coordinate_vector<double>::size_type>())
    //.def("NonZeros", &coordinate_vector<double>::nnz)
    //.def(VectorScalarOperatorPython<coordinate_vector<double>, double, coordinate_vector<double> >())
      //.def(VectorVectorOperatorPython<coordinate_vector<double>, zero_vector<double>, coordinate_vector<double> >())
      //.def(VectorVectorOperatorPython<coordinate_vector<double>, unit_vector<double>, coordinate_vector<double> >())
      //.def(VectorVectorOperatorPython<coordinate_vector<double>, scalar_vector<double>, coordinate_vector<double> >())
      //.def(VectorVectorOperatorPython<coordinate_vector<double>, vector<double>, vector<double> >())
      //;
		


// // 		MatrixFromPython<Matrix<double> >();
		
// 		class_<Matrix<double> >("DenseMatrix", init<unsigned int, unsigned int>())
// 			.def(init<Matrix<double> >())
// 			.def("Size1", &Matrix<double>::RowsNumber)
// 			.def("Size2", &Matrix<double>::ColumnsNumber)
// 			.def("__len__", &Matrix<double>::RowsNumber)
// 			.def("__getitem__", &matrix_getrow<Matrix<double> >)
// 			.def("__getitem__", &matrix_getitem<Matrix<double> >)
// 			.def("__setitem__", &matrix_setrow<Matrix<double> >)
// 			.def("__setitem__", &matrix_setitem<Matrix<double> >)
// 			.def(self *= double())          
// 			.def(self /= double())          
// 			.def(self * Vector<double>())          
// 			.def(self += self)          
// 			.def(self -= self)          
// 			.def(self + self)          
// 			.def(self - self)          
// 			.def(self * self)          
// 			.def(self_ns::str(self))
// 			;
		
// 		class_<CSRMatrix<double> >("CSRMatrix", init<unsigned int>())
// 			.def("Size1", &CSRMatrix<double>::RowsNumber)
// 			.def("Size2", &CSRMatrix<double>::RowsNumber)
// 			.def("__len__", &CSRMatrix<double>::RowsNumber)
// 			.def("__getitem__", &matrix_getitem<CSRMatrix<double> >)
// 			.def("__setitem__", &matrix_setitem<CSRMatrix<double> >)
// 			.def("Pattern", &CSRMatrix<double>::PrintPattern)
// 			// 	 .def(self * Vector<double>())          
// 			.def(self_ns::str(self))
// 			;
 	}
	
}  // namespace Python.

} // Namespace Kratos

//   template<class TContainerType> struct MatrixFromPython
//   {
// 			MatrixFromPython()
// 			{
// 				boost::python::converter::registry::push_back(
// 					&convertible,
// 					&construct,
// 					boost::python::type_id<TContainerType>());
// 			}
			
// 			static void* convertible(PyObject* obj_ptr)
// 			{
// 				using namespace boost::python;
// 				using boost::python::allow_null; // works around gcc 2.96 bug
// 				{
// 					// Restriction to list, tuple, iter, xrange until
// 					// Boost.Python overload resolution is enhanced.
// 					if (!(   PyList_Check(obj_ptr)
// 						|| PyTuple_Check(obj_ptr)
// 						|| PyIter_Check(obj_ptr)
// 						|| PyRange_Check(obj_ptr))) return 0;
// 				}
// 				handle<> obj_iter(allow_null(PyObject_GetIter(obj_ptr)));
// 				if (!obj_iter.get()) { // must be convertible to an iterator
// 					PyErr_Clear();
// 					return 0;
// 				}
// 				int obj_size = PyObject_Length(obj_ptr);
// 				if (obj_size < 0) { // must be a measurable sequence
// 					PyErr_Clear();
// 					return 0;
// 				}
				
// 				bool is_range = PyRange_Check(obj_ptr);
// 				std::size_t i=0;
// 				for(;;i++) {
// 					handle<> py_elem_hdl(allow_null(PyIter_Next(obj_iter.get())));
// 					if (PyErr_Occurred()) {
// 						PyErr_Clear();
// 						return 0;
// 					}
// 					if (!py_elem_hdl.get()) break; // end of iteration
// 					object py_elem_obj(py_elem_hdl);
// 					extract<Vector<double> > elem_proxy(py_elem_obj);
// 					if (!elem_proxy.check()) return 0;
// 					if (is_range) break; // in a range all elements are of the same type
// 				}
// 				if (!is_range) assert(i == obj_size);
// 			return obj_ptr;
// 		}
		
// 		static void construct(
// 			PyObject* obj_ptr,
// 			boost::python::converter::rvalue_from_python_stage1_data* data)
// 		{
// 			using namespace boost::python;
// 			using boost::python::allow_null; // works around gcc 2.96 bug
// 			using boost::python::converter::rvalue_from_python_storage; // dito
// 			using boost::python::throw_error_already_set; // dito
// 			handle<> obj_iter(PyObject_GetIter(obj_ptr));
// 			void* storage = (
// 				(rvalue_from_python_storage<TContainerType>*)
// 				data)->storage.bytes;
// 			new (storage) TContainerType();
// 			data->convertible = storage;
// 			TContainerType& result = *((TContainerType*)storage);
// 			Vector<Vector<double> > temp;
// 			for(unsigned int i = 0;;i++) {
// 				handle<> py_elem_hdl(allow_null(PyIter_Next(obj_iter.get())));
// 				if (PyErr_Occurred()) throw_error_already_set();
// 				if (!py_elem_hdl.get()) break; // end of iteration
// 				object py_elem_obj(py_elem_hdl);
// 				extract<Vector<double> > elem_proxy(py_elem_obj);
// 				temp.push_back(elem_proxy());
// 			}
// 			result = Matrix<double>(temp);
// 		}
//   };
	
// 	template <class TMatrix>
// 		void matrix_setrow(TMatrix& m, unsigned int i, const Vector<double>& value)
// 	{
// 		if (i< m.RowsNumber())
// 			m[i] = value;
// 		else {
// 			PyErr_SetString(PyExc_IndexError, "index out of range");
// 			throw_error_already_set();
// 		}
// 	}
	
// 	template <class TMatrix>
// 		Vector<double> matrix_getrow(TMatrix &m, unsigned int i)
// 	{
// 		if (i < m.RowsNumber())
// 			return m[i];
// 		else {
// 			PyErr_SetString(PyExc_IndexError, "index out of range");
// 			throw_error_already_set();
// 		}
		
// 		return Vector<double>();
// 	}
	
// 	template <class TMatrix>
// 		void matrix_setitem(TMatrix& m, tuple &index, double value)
// 	{
// 		unsigned int i = extract<double>(index[0]);
// 		unsigned int j = extract<double>(index[1]);
// 		if (i< m.RowsNumber())
// 			if (j < m.ColumnsNumber())
// 				m(i,j) = value;
// 			else {
// 				PyErr_SetString(PyExc_IndexError, "index out of range");
// 				throw_error_already_set();
// 			}
// 	}
	
// 	template <class TMatrix>
// 		double matrix_getitem(TMatrix &m, tuple &index)
// 	{
// 		unsigned int i = extract<double>(index[0]);
// 		unsigned int j = extract<double>(index[1]);
// 		if (i < m.RowsNumber())
// 			if (j < m.ColumnsNumber())
// 				return m(i,j);
// 			else {
// 				PyErr_SetString(PyExc_IndexError, "index out of range");
// 				throw_error_already_set();
// 			}
			
// 			return double();
// 	}
// 	template<class TDataType>
// 		Vector<TDataType> operator*(const CSRMatrix<TDataType>& rMatrix, const Vector<TDataType>& rVector)
// 	{
// 		Vector<TDataType> result;
// 		CSRSpace<TDataType>::Mult(rMatrix, rVector, result);
// 		return result;
// 	}
