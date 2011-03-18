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
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:34 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_VECTOR_PYTHON_INTERFACE_H_INCLUDED )
#define KRATOS_VECTOR_PYTHON_INTERFACE_H_INCLUDED



// System includes 


// External includes 


// Project includes
#include "includes/define.h"
#include "python/readonly_vector_python_interface.h"
#include "python/vector_vector_operator_python.h"
#include "python/vector_scalar_assignment_operator_python.h"
#include "python/vector_vector_assignment_operator_python.h"


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
		template<class TContainerType, class TModifierType>
		class VectorPythonInterface : public ReadonlyVectorPythonInterface<TContainerType>
		{
		public:
			///@name Type Definitions
			///@{

			/// Pointer definition of VectorPythonInterface
			KRATOS_CLASS_POINTER_DEFINITION(VectorPythonInterface);

			typedef typename TContainerType::value_type data_type;
			typedef typename TContainerType::value_type key_type;
			typedef typename TContainerType::size_type index_type;
			typedef typename TContainerType::size_type size_type;
			typedef typename TContainerType::difference_type difference_type;

			///@}
			///@name Life Cycle 
			///@{ 

			/// Default constructor.
			VectorPythonInterface()
			{
			}

			/// Destructor.
			virtual ~VectorPythonInterface(){}


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
				if (!obj_iter.get()) { // must be convertible to an iterator
					PyErr_Clear();
					return 0;
				}
				if (true) {
					int obj_size = PyObject_Length(obj_ptr);
					if (obj_size < 0) { // must be a measurable sequence
						PyErr_Clear();
						return 0;
					}
					/*         if (!TConversionPolicyType::check_size( */
					/*           boost::type<TTContainerTypeType>(), obj_size)) return 0; */
					bool is_range = PyRange_Check(obj_ptr);
					std::size_t i=0;
					for(;;i++) {
						handle<> py_elem_hdl(allow_null(PyIter_Next(obj_iter.get())));
						if (PyErr_Occurred()) {
							PyErr_Clear();
							return 0;
						}
						if (!py_elem_hdl.get()) break; // end of iteration
						object py_elem_obj(py_elem_hdl);
						extract<typename TContainerType::value_type> elem_proxy(py_elem_obj);
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
				handle<> obj_iter(PyObject_GetIter(obj_ptr));
				void* storage = (
					(rvalue_from_python_storage<TContainerType>*)
					data)->storage.bytes;
				new (storage) TContainerType();
				data->convertible = storage;
				TContainerType& result = *((TContainerType*)storage);
				std::size_t i=0;
				for(;;i++) {
					handle<> py_elem_hdl(allow_null(PyIter_Next(obj_iter.get())));
					if (PyErr_Occurred()) throw_error_already_set();
					if (!py_elem_hdl.get()) break; // end of iteration
					object py_elem_obj(py_elem_hdl);
					extract<typename TContainerType::value_type> elem_proxy(py_elem_obj);
					append(result, elem_proxy());
				}
				/*       TConversionPolicyType::assert_size(boost::type<TContainerType>(), i); */
			}

			template <class Class>
				static void 
				extension_def(Class& cl)
			{
				cl
					.def("append", &base_append)
					.def("extend", &base_extend)
					;
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
					TModifierType::MoveSlice(container, to + dist, to, size);
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
					TModifierType::MoveSlice(container, from, to, container.size());
				}
			}


			static void 
				append(TContainerType& container, data_type const& v)
			{ 
				index_type i = container.size();
				TModifierType::Resize(container, i + 1);
				container.insert_element(i, v);
			}

			template <class Iter>
				static void 
				extend(TContainerType& container, Iter first, Iter last)
			{ 
				index_type i = container.size();
				TModifierType::Resize(container, i + std::distance(first,last));
				for( ; first != last ; first++)
					container.insert_element(i++, *first);
			}


			static class_<TContainerType> CreateInterface(std::string const& Name)
			{
				boost::python::converter::registry::push_back(
					&convertible,
					&construct,
					boost::python::type_id<TContainerType>());

				return class_<TContainerType> (Name.c_str())
					.def(init<TContainerType>())
					.def(indexing_suite<TContainerType, VectorPythonInterface>())
				  //.def("Resize", &TContainerType::resize)
					.def("Size", &TContainerType::size)
					.def(VectorVectorOperatorPython<TContainerType, TContainerType, TContainerType>())
					.def(VectorScalarAssignmentOperatorPython<TContainerType, double>())
					.def(VectorVectorAssignmentOperatorPython<TContainerType, zero_vector<double> >())
					.def(VectorVectorAssignmentOperatorPython<TContainerType, unit_vector<double> >())
					.def(VectorVectorAssignmentOperatorPython<TContainerType, scalar_vector<double> >())
					.def(VectorVectorAssignmentOperatorPython<TContainerType, vector<double> >())
				  //.def(VectorVectorAssignmentOperatorPython<TContainerType, mapped_vector<double> >())
				  //	.def(VectorVectorAssignmentOperatorPython<TContainerType, compressed_vector<double> >())
				  //	.def(VectorVectorAssignmentOperatorPython<TContainerType, coordinate_vector<double> >())
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
			{return ThisContainer + OtherContainer;}

			static TContainerType sub_self(TContainerType& ThisContainer, TContainerType const& OtherContainer)
			{return ThisContainer - OtherContainer;}

			static data_type mul_self(TContainerType& ThisContainer, TContainerType const& OtherContainer)
			{return inner_prod(ThisContainer, OtherContainer);}


			static void
				base_append(TContainerType& container, object v)
			{
				extract<data_type&> elem(v);
				// try if elem is an exact Data
				if (elem.check())
				{
					append(container, elem());
				}
				else
				{
					//  try to convert elem to data_type
					extract<data_type> elem(v);
					if (elem.check())
					{
						append(container, elem());
					}
					else
					{
						PyErr_SetString(PyExc_TypeError, 
							"Attempting to append an invalid type");
						throw_error_already_set();
					}
				}
			}

			static void
				base_extend(TContainerType& container, object v)
			{
				std::vector<data_type> temp;
				container_utils::extend_container(temp, v);
				extend(container, temp.begin(), temp.end());
			}
		}; // Class VectorPythonInterface 

		///@} 

		///@name Type Definitions       
		///@{ 


		///@} 
		///@name Input and output 
		///@{ 


		///@} 


	}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_VECTOR_PYTHON_INTERFACE_H_INCLUDED defined 


