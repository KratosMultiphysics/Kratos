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


#if !defined(KRATOS_MATRIX_SCALAR_OPERATOR_PYTHON_H_INCLUDED )
#define  KRATOS_MATRIX_SCALAR_OPERATOR_PYTHON_H_INCLUDED



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
		template<class TMatrixType, class TScalarType, class TAddResultType = TMatrixType, class TMultResultType = TMatrixType>
		class MatrixScalarOperatorPython : public def_visitor<MatrixScalarOperatorPython<TMatrixType, TScalarType, TAddResultType, TMultResultType> >
		{
		public:
			///@name Type Definitions
			///@{

			/// Pointer definition of MatrixScalarOperatorPython
			KRATOS_CLASS_POINTER_DEFINITION(MatrixScalarOperatorPython);

			///@}
			///@name Life Cycle 
			///@{ 

			/// Default constructor.
			MatrixScalarOperatorPython(){}

			/// Copy constructor.
			MatrixScalarOperatorPython(const MatrixScalarOperatorPython& rOther){}

			/// Destructor.
			virtual ~MatrixScalarOperatorPython(){}


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
					.def("__add__", &add)           
					.def("__sub__", &sub)           
					.def("__mul__", &mul)          
					.def("__div__", &div)          
					.def("__radd__", &radd)           
					.def("__rsub__", &rsub)           
					.def("__rmul__", &rmul)          
					//.def("__rdiv__", &rdiv)          
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
				TAddResultType
				add(TMatrixType& ThisMatrix, TScalarType ThisScalar)
			{
				return ThisMatrix + scalar_matrix<TScalarType>(ThisMatrix.size1(), ThisMatrix.size2(), ThisScalar);
			}

			static
				TAddResultType
				sub(TMatrixType& ThisMatrix, TScalarType ThisScalar)
			{
				return ThisMatrix - scalar_matrix<TScalarType>(ThisMatrix.size1(), ThisMatrix.size2(), ThisScalar);
			}

			static 
				TMultResultType
				mul(TMatrixType& ThisMatrix, TScalarType ThisScalar)
			{return ThisMatrix * ThisScalar;}

			static 
				TMultResultType
				div(TMatrixType& ThisMatrix, TScalarType ThisScalar)
			{return ThisMatrix / ThisScalar;}

			static 
				TAddResultType
				radd(TMatrixType& ThisMatrix, TScalarType ThisScalar)
			{
				return scalar_matrix<TScalarType>(ThisMatrix.size1(), ThisMatrix.size2(), ThisScalar) + ThisMatrix;
			}

			static 
				TAddResultType
				rsub(TMatrixType& ThisMatrix, TScalarType ThisScalar)
			{
				return scalar_matrix<TScalarType>(ThisMatrix.size1(), ThisMatrix.size2(), ThisScalar) - ThisMatrix;
			}

			static 
				TMultResultType
				rmul(TMatrixType& ThisMatrix, TScalarType ThisScalar)
			{return ThisScalar * ThisMatrix;}

			//static 
			//	TMultResultType
			//	rdiv(TMatrixType& ThisMatrix, TScalarType ThisScalar)
			//{
			//	return scalar_matrix<TScalarType>(ThisMatrix.size1(), ThisMatrix.size2(), ThisScalar) / ThisMatrix;
			//}



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
			MatrixScalarOperatorPython& operator=(const MatrixScalarOperatorPython& rOther);


			///@}    

		}; // Class MatrixScalarOperatorPython 

		///@} 

		///@name Type Definitions       
		///@{ 


		///@} 
		///@name Input and output 
		///@{ 

		///@} 

	}  // namespace Python.

}  // namespace Kratos.

#endif // KRATOS_MATRIX_SCALAR_OPERATOR_PYTHON_H_INCLUDED  defined 


