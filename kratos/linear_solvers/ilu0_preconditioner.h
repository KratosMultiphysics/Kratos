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
//   Date:                $Date: 2007-03-06 10:30:33 $
//   Revision:            $Revision: 1.2 $
//
//



#if !defined(KRATOS_ILU0_PRECONDITIONER_H_INCLUDED )
#define  KRATOS_ILU0_PRECONDITIONER_H_INCLUDED




// System includes
#include <algorithm> 



// External includes 
#include "boost/smart_ptr.hpp"



// Project includes
#include "includes/define.h"
#include "linear_solvers/ilu_preconditioner.h"



namespace Kratos
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

	///@name  Preconditioners 
	///@{ 

	/// ILU0Preconditioner class. 
	/**   */
	template<class TSparseSpaceType, class TDenseSpaceType>
	class ILU0Preconditioner : public ILUPreconditioner<TSparseSpaceType, TDenseSpaceType>
	{
	public:
		///@name Type Definitions
		///@{

		/// Counted pointer of ILU0Preconditioner
		typedef boost::shared_ptr<ILU0Preconditioner> Pointer;


		typedef ILUPreconditioner<TSparseSpaceType, TDenseSpaceType> BaseType;


		typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

		typedef typename TSparseSpaceType::VectorType VectorType;

		typedef typename TDenseSpaceType::MatrixType DenseMatrixType;


		///@}
		///@name Life Cycle 
		///@{ 

		/// Default constructor.
		ILU0Preconditioner(){
			BaseType::L = BaseType::U = NULL;
			BaseType::iL = BaseType::jL = BaseType::iU = BaseType::jU = NULL;
		}


		/// Copy constructor.
// 		ILU0Preconditioner(const ILU0Preconditioner& Other){}


		/// Destructor.
		virtual ~ILU0Preconditioner()
		{ 
			if ( BaseType::L!=NULL) delete[]  BaseType::L;
			if (BaseType::iL!=NULL) delete[] BaseType::iL;
			if (BaseType::jL!=NULL) delete[] BaseType::jL;
			if ( BaseType::U!=NULL) delete[]  BaseType::U;
			if (BaseType::iU!=NULL) delete[] BaseType::iU;
			if (BaseType::jU!=NULL) delete[] BaseType::jU;
			
			BaseType::L = NULL;
			BaseType::iL = NULL;
			BaseType::jL = NULL;
			BaseType::U = NULL;
			BaseType::iU = NULL;
			BaseType::jU = NULL;
		}



		///@}
		///@name Operators 
		///@{

		/// Assignment operator.
// 		ILU0Preconditioner& operator=(const ILU0Preconditioner& Other)
// 		{
// 			BaseType::operator=(Other);
// 			return *this;
// 		}



		///@}
		///@name Operations
		///@{

		/** ILU0Preconditioner Initialize
		Initialize preconditioner for linear system rA*rX=rB
		@param rA  system matrix.
		@param rX Unknows vector
		@param rB Right side linear system of equations.
		*/
		virtual void Initialize(SparseMatrixType& rA, VectorType& rX, VectorType& rB) 
		{
			// ILU(0) preconditioner
			// Incomplete LU factorization with same sparcity pattern as original matrix.
			// The diagonal is included in BaseType::U. BaseType::L has non-written 1's on its diagonal.
			// See pg 274 Iterative methods for linear systems, Yousef Saad 
			// We assume that, within a row, the entries in A are sorted by increasing j
			const int size = TSparseSpaceType::Size(rX);
			int i, j, indexj, oldCountL, oldCountU, newCountL, newCountU, fillL, fillU;
			int k, indexk, indexkj, indexkjlim, jkj, indexjlim;
			double aik;
			bool diagFound;

			BaseType::mILUSize= size;
			int n = size;
			// in case a preconditioner is mistakenly initialized twice;
			if ( BaseType::L!=NULL) delete[]  BaseType::L;
			if (BaseType::iL!=NULL) delete[] BaseType::iL;
			if (BaseType::jL!=NULL) delete[] BaseType::jL;
			if ( BaseType::U!=NULL) delete[]  BaseType::U;
			if (BaseType::iU!=NULL) delete[] BaseType::iU;
			if (BaseType::jU!=NULL) delete[] BaseType::jU;
			
			BaseType::L = NULL;
			BaseType::iL = NULL;
			BaseType::jL = NULL;
			BaseType::U = NULL;
			BaseType::iU = NULL;
			BaseType::jU = NULL;


			// Create copy of matrix split in its BaseType::L and BaseType::U parts
			// Traverse matrix to count elements in rows of BaseType::L, BaseType::U
			// If there is no element in the diagonal, make room for a zero
			BaseType::iL=new int[n+1];
			BaseType::iU=new int[n+1];
			for (i=0; i<n+1; i++) 
			{	
				BaseType::iL[i]=0; BaseType::iU[i]=0;
			}

			for (typename SparseMatrixType::iterator1 a_iterator = rA.begin1(); a_iterator != rA.end1(); a_iterator++) {
				i = a_iterator.index1();
				diagFound=false;
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
				for (typename SparseMatrixType::iterator2 row_iterator = a_iterator.begin() ; row_iterator != a_iterator.end() ; ++row_iterator) {
#else
				for ( SparseMatrixType::iterator2 row_iterator = begin(a_iterator, boost::numeric::ublas::iterator1_tag()); row_iterator != end(a_iterator, boost::numeric::ublas::iterator1_tag()); ++row_iterator ){
#endif
					//j= row_iterator->first; // Changed by Pooyan!!!
					j= row_iterator.index2(); // Changed by Pooyan!!!
					if (i<=j) BaseType::iU[i]=BaseType::iU[i]+1;
					if (i>j)  BaseType::iL[i]=BaseType::iL[i]+1;
					if (i==j) diagFound=true;
				}
				if (!diagFound) BaseType::iU[i]=BaseType::iU[i]+1;
			}

			// BaseType::iL[i] is now the number of nonzero entries in the 
			// i-th row of BaseType::L. Transform to CSR indexes

			oldCountL=0; oldCountU=0; 
			for (i=0; i<n+1; i++) {
				newCountL=oldCountL+BaseType::iL[i];
				BaseType::iL[i]=oldCountL;
				oldCountL=newCountL;
				newCountU=oldCountU+BaseType::iU[i];
				BaseType::iU[i]=oldCountU;
				oldCountU=newCountU;
			}

			// Traverse again the matrix copying its entries
			BaseType::L =new double[BaseType::iL[n]];
			BaseType::jL=new int   [BaseType::iL[n]];
			BaseType::U =new double[BaseType::iU[n]];
			BaseType::jU=new int   [BaseType::iU[n]];
			fillL=0;
			fillU=0;
			//for (i=0; i<n; i++) {
			for (typename SparseMatrixType::iterator1 a_iterator = rA.begin1(); a_iterator != rA.end1(); a_iterator++) {
				i = a_iterator.index1();
				diagFound=false;
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
				for (typename SparseMatrixType::iterator2 row_iterator = a_iterator.begin() ; row_iterator != a_iterator.end() ; ++row_iterator) {
#else
				for ( SparseMatrixType::iterator2 row_iterator = begin(a_iterator, boost::numeric::ublas::iterator1_tag()); row_iterator != end(a_iterator, boost::numeric::ublas::iterator1_tag()); ++row_iterator ){
#endif
					j= row_iterator.index2(); // Changed by Pooyan!!!
					if (i==j) diagFound=true;
					if ( (j>i) && (!diagFound) ) {
						// This row does not have a diagonal entry. Make it and put a zero.
						BaseType::jU[fillU]=i;
						BaseType::U[fillU]=0.00;
						fillU++;
						diagFound=true;
					}
					if (i<=j) {
						BaseType::jU[fillU]=j; 
						BaseType::U[fillU]= *row_iterator; 
						fillU++;
					}
					if (i>j) {
						BaseType::jL[fillL]=j; 
						BaseType::L[fillL]= *row_iterator;
						fillL++;
					}
				}
			}


			// Now comes the real factorization:
			// for i=2, ... ,n
			//    for k=1, ... ,i-1    and (i,k) in nonzero pattern
			//       a[i,k]=a[i,k]/a[k,k]
			//       for j=k+1, ... ,n     and (i,j) in nonzero pattern
			//          a[i,j]=a[i,j]-a[i,k]*a[k,j]
			//       end do
			//    end do
			// end do


			for (i=1; i<n; i++) {
				for (indexk=BaseType::iL[i]; indexk<BaseType::iL[i+1]; indexk++) {
					k=BaseType::jL[indexk];

					BaseType::L[indexk]=BaseType::L[indexk]/BaseType::U[BaseType::iU[k]];
					aik=BaseType::L[indexk];

					indexkj   = BaseType::iU[k  ];  // traverses row k of BaseType::U
					indexkjlim= BaseType::iU[k+1];
					jkj       = BaseType::jU[indexkj];
					while (jkj<=k) jkj=BaseType::jU[++indexkj]; // add rows only for j>k
//****************************************************************************************
//****************************************************************************************+
//by Riccardo ... here an array bound is broken
//between this line
					indexj    = indexk+1; // traverses row i of BaseType::L beyond k
					indexjlim = BaseType::iL[i+1];
					j         = BaseType::jL[indexj];
//and this line
//the memory is broken in accessing to BaseType::jL[indexj]; when i=n-1
//when it works j is read to be 0 ... but it could be whatever value
//****************************************************************************************+
//****************************************************************************************+

					while ( (indexkj<indexkjlim) && (indexj<indexjlim) ) {
						if (j==jkj) {
							BaseType::L[indexj]= BaseType::L[indexj]-aik*BaseType::U[indexkj];
							indexj++;
							j=BaseType::jL[indexj];
							indexkj++;
							jkj=BaseType::jU[indexkj];

						} else {
							if (j<jkj) {
								indexj++;
								j=BaseType::jL[indexj];
							} else {
								indexkj++;
								jkj=BaseType::jU[indexkj];
							}
						}
					}

					indexj    = BaseType::iU[i]; // traverses row i of BaseType::U
					indexjlim = BaseType::iU[i+1];
					j         = BaseType::jU[indexj];
					while ( (indexkj<indexkjlim) && (indexj<indexjlim) ) {
						if (j==jkj) {
							BaseType::U[indexj]= BaseType::U[indexj]-aik*BaseType::U[indexkj];
							indexj++;
							j=BaseType::jU[indexj];
							indexkj++;
							jkj=BaseType::jU[indexkj];
						} else {
							if (j<jkj) {
								indexj++;
								j=BaseType::jU[indexj];
							} else {
								indexkj++;
								jkj=BaseType::jU[indexkj];
							}
						}
					}
				}
			}

			// and that's it...
			// for debugging only; write out matrices BaseType::L, D, BaseType::U
			/*      std::cout << "BaseType::L, BaseType::U:\n"; */
			/*      for (i=0; i<n; i++) { */
			/*              for (indexj=BaseType::iL[i]; indexj<BaseType::iL[i+1]; indexj++) { */
			/*                      j=BaseType::jL[indexj]; */
			/*                      std::cout << "BaseType::L[" << i << "," << j << "]=" << BaseType::L[indexj] << "\n"; */
			/*              } */
			/*      } */
			/*      for (i=0; i<n; i++) { */
			/*              for (indexj=BaseType::iU[i]; indexj<BaseType::iU[i+1]; indexj++) { */
			/*                      j=BaseType::jU[indexj]; */
			/*                      std::cout << "BaseType::U[" << i << "," << j << "]=" << BaseType::U[indexj] << "\n"; */
			/*              } */
			/*      } */



			for (i=0; i<n; i++) if (BaseType::U[BaseType::iU[i]]==0.00) 
			{
				KRATOS_ERROR(std::runtime_error, "Zero in BaseType::U diagonal found!!", "")
			}
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

		/// Return information about this object.
		virtual std::string Info() const
		{
			return "ILU0Preconditioner";
		}


		/// Print information about this object.
		virtual void  PrintInfo(std::ostream& OStream) const
		{
			OStream << "ILU0Preconditioner";
		}


		virtual void PrintData(std::ostream& OStream) const
		{
		}



		///@}      
		///@name Friends
		///@{


		///@}

	protected:
		///@name Protected static Member Variables 
		///@{ 


		///@} 
		///@name Protected member Variables 
		///@{ 


		///@} 
		///@name Protected Operators
		///@{ 


		///@} 
		///@name Protected Operations
		///@{ 


		///@} 
		///@name Protected  Access 
		///@{ 


		///@}      
		///@name Protected Inquiry 
		///@{ 


		///@}    
		///@name Protected LifeCycle 
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


		///@} 
		///@name Private  Access 
		///@{ 


		///@}    
		///@name Private Inquiry 
		///@{ 


		///@}    
		///@name Un accessible methods 
		///@{ 


		///@}    

	}; // Class ILU0Preconditioner 

	///@} 


	///@} 

	///@name Type Definitions       
	///@{ 


	///@} 
	///@name Input and output 
	///@{ 


	/// input stream function
	template<class TSparseSpaceType, class TDenseSpaceType>
		inline std::istream& operator >> (std::istream& IStream, 
		ILU0Preconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
	{
		return IStream;
	}


	/// output stream function
	template<class TSparseSpaceType, class TDenseSpaceType>
		inline std::ostream& operator << (std::ostream& OStream, 
		const ILU0Preconditioner<TSparseSpaceType, TDenseSpaceType>& rThis)
	{
		rThis.PrintInfo(OStream);
		OStream << std::endl;
		rThis.PrintData(OStream);


		return OStream;
	}
	///@} 


}  // namespace Kratos.


#endif // KRATOS_ILU0_PRECONDITIONER_H_INCLUDED  defined 

