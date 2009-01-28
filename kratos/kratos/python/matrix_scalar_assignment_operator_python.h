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


#if !defined(KRATOS_MATRIX_SCALAR_ASSIGNMENT_OPERATOR_PYTHON_H_INCLUDED )
#define  KRATOS_MATRIX_SCALAR_ASSIGNMENT_OPERATOR_PYTHON_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


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

		/// Short class definition.
		/** Detail class definition.
		*/
		template<class TMatrixType, class TScalarType>
		class MatrixScalarAssignmentOperatorPython : public def_visitor<MatrixScalarAssignmentOperatorPython<TMatrixType, TScalarType> >
		{
		public:
			///@name Type Definitions
			///@{

			/// Pointer definition of MatrixScalarAssignmentOperatorPython
			KRATOS_CLASS_POINTER_DEFINITION(MatrixScalarAssignmentOperatorPython);

			///@}
			///@name Life Cycle 
			///@{ 

			/// Default constructor.
			MatrixScalarAssignmentOperatorPython(){}

			/// Copy constructor.
			MatrixScalarAssignmentOperatorPython(const MatrixScalarAssignmentOperatorPython& rOther);

			/// Destructor.
			virtual ~MatrixScalarAssignmentOperatorPython(){}


			///@}
			///@name Operators 
			///@{


			///@}
			///@name Operations
			///@{

			template <class TClassType>
				void visit(TClassType& ThisClass) const
			{
				ThisClass
					.def("__iadd__", &iadd, return_internal_reference<1>())           
					.def("__isub__", &isub, return_internal_reference<1>())          
					.def("__imul__", &imul, return_internal_reference<1>())           
					.def("__idiv__", &idiv, return_internal_reference<1>())          
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
			///@name Static Member Variables 
			///@{ 


			///@} 
			///@name Member Variables 
			///@{ 


			///@} 
			///@name Private Operators
			///@{ 


			///@} 
			///@name Private Operations
			///@{

			static 
				TMatrixType&
				iadd(TMatrixType& ThisMatrix, TScalarType ThisScalar)
			{
				return ThisMatrix += scalar_matrix<TScalarType>(ThisMatrix.size1(), ThisMatrix.size2(), ThisScalar);
			}

			static
				TMatrixType&
				isub(TMatrixType& ThisMatrix, TScalarType ThisScalar)
			{
				return ThisMatrix += scalar_matrix<TScalarType>(ThisMatrix.size1(), ThisMatrix.size2(), ThisScalar);
			}

			static 
				TMatrixType&
				imul(TMatrixType& ThisMatrix, TScalarType ThisScalar)
			{
				return ThisMatrix *= ThisScalar;
			}

			static 
				TMatrixType&
				idiv(TMatrixType& ThisMatrix, TScalarType ThisScalar)
			{
				return ThisMatrix /= ThisScalar;
			}



			///@} 
			///@name Private  Access 
			///@{ 


			///@}    
			///@name Private Inquiry 
			///@{ 


			///@}    
			///@name Un accessible methods 
			///@{ 

			/// Assignment operator.
			MatrixScalarAssignmentOperatorPython& operator=(const MatrixScalarAssignmentOperatorPython& rOther);


			///@}    

		}; // Class MatrixScalarAssignmentOperatorPython 

		///@} 

		///@name Type Definitions       
		///@{ 


		///@} 
		///@name Input and output 
		///@{ 

		///@} 

	}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_MATRIX_SCALAR_ASSIGNMENT_OPERATOR_PYTHON_H_INCLUDED  defined 


