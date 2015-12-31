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


