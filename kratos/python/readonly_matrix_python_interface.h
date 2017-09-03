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




#if !defined(KRATOS_READONLY_MATRIX_PYTHON_INTERFACE_H_INCLUDED )
#define KRATOS_READONLY_MATRIX_PYTHON_INTERFACE_H_INCLUDED



// System includes


// External includes


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
template<class TMatrixType>
class ReadonlyMatrixPythonInterface
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ReadonlyMatrixPythonInterface
    KRATOS_CLASS_POINTER_DEFINITION(ReadonlyMatrixPythonInterface);

    typedef typename TMatrixType::value_type data_type;
    typedef typename TMatrixType::value_type key_type;
    typedef typename TMatrixType::size_type index_type;
    typedef typename TMatrixType::size_type size_type;
    typedef typename TMatrixType::difference_type difference_type;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Destructor.
    virtual ~ReadonlyMatrixPythonInterface() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    static
    data_type get_item(TMatrixType const& ThisMatrix, tuple index)
    {
        unsigned int i = extract<index_type>(index[0]);
        unsigned int j = extract<index_type>(index[1]);
        if ((i >= ThisMatrix.size1()) || (j >= ThisMatrix.size2()))
        {
            PyErr_SetString(PyExc_IndexError, "index out of range");
            throw_error_already_set();
        }

        return ThisMatrix(i,j);
    }

    /*       static  */
    /*       void resize1(TMatrixType& ThisMatrix, size_type NewSize) */
    /*       {	 */
    /*           ThisMatrix.resize(NewSize, NewSize); */
    /*       } */

    /*       static  */
    /*       void resize2(TMatrixType& ThisMatrix, size_type NewSize1, size_type NewSize2) */
    /*       {	 */
    /*           ThisMatrix.resize(NewSize1, NewSize2); */
    /*       } */



    static class_<TMatrixType> CreateInterface(std::string const& Name)
    {
        return class_<TMatrixType>(Name.c_str())
               .def(init<TMatrixType>())
               //.def("Resize", &resize1)
               //.def("Resize", &resize2)
               .def("Size1", &TMatrixType::size1)
               .def("Size2", &TMatrixType::size2)
               .def("__getitem__", &get_item)
//	  .def(self_ns::str(self))
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

}; // Class ReadonlyMatrixPythonInterface

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_READONLY_MATRIX_PYTHON_INTERFACE_H_INCLUDED defined 


