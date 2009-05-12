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
//   Date:                $Date: 2008-11-10 14:23:33 $
//   Revision:            $Revision: 1.7 $
//
//


#if !defined(KRATOS_UBLAS_SPACE_H_INCLUDED )
#define  KRATOS_UBLAS_SPACE_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <cstddef>


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"


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

/// Short class definition.
/** Detail class definition.
*/
template<class TDataType, class TMatrixType, class TVectorType>
class UblasSpace
{
	public:
	///@name Type Definitions
	///@{
	
	/// Pointer definition of UblasSpace
	KRATOS_CLASS_POINTER_DEFINITION(UblasSpace);
	
	typedef TDataType DataType;
	
	typedef TMatrixType MatrixType;
	
	typedef TVectorType VectorType;
	
	typedef std::size_t IndexType;
	
	typedef std::size_t SizeType;
	
	typedef typename boost::shared_ptr< TMatrixType > MatrixPointerType;
	typedef typename boost::shared_ptr< TVectorType > VectorPointerType;
	
	///@}
	///@name Life Cycle 
	///@{ 
	
	/// Default constructor.
	UblasSpace(){}
	
	/// Destructor.
	virtual ~UblasSpace(){}
	
	
	///@}
	///@name Operators 
	///@{
	
	
	///@}
	///@name Operations
	///@{
	static MatrixPointerType CreateEmptyMatrixPointer(){return MatrixPointerType( new TMatrixType() ); }
	static VectorPointerType CreateEmptyVectorPointer(){return VectorPointerType( new TVectorType() ); }
	
	/// return size of vector rV
	static IndexType Size(VectorType const& rV){return rV.size();} 
	
	/// return number of rows of rM 
	static IndexType Size1(MatrixType const& rM){return rM.size1();}
	
	/// return number of columns of rM
	static IndexType Size2(MatrixType const& rM){return rM.size2();} 
	
	/// rXi = rMij
	static void GetColumn(unsigned int j, MatrixType& rM, VectorType& rX){rX = column(rM, j);} 
	
	
	///////////////////////////////// TODO: Take a close look to this method!!!!!!!!!!!!!!!!!!!!!!!!!
	/// rMij = rXi
	static void SetColumn(unsigned int j, MatrixType& rM, VectorType& rX){rX = row(rM, j);}  
	
	/// rY = rX
	static void Copy(MatrixType const& rX, MatrixType& rY){rY.assign(rX);} 
	
	/// rY = rX
	static void Copy(VectorType const& rX, VectorType& rY){rY.assign(rX);} 
	
	/// rX * rY
	static double Dot(VectorType& rX, VectorType& rY){return inner_prod(rX, rY);} 
	
	/// ||rX||2
	static double TwoNorm(VectorType const& rX){return norm_2(rX);} 
	
	
	static void Mult(MatrixType& rA, VectorType& rX, VectorType& rY)
		{
		axpy_prod(rA, rX, rY, true);
		}
	
	static void TransposeMult(MatrixType& rA, VectorType& rX, VectorType& rY)
		{
		axpy_prod(rX, rA, rY, true);
		} // rY = rAT * rX
	
		static inline SizeType GraphDegree( IndexType i, TMatrixType& A)
		{
		typename MatrixType::iterator1 a_iterator = A.begin1();
		std::advance(a_iterator,i);
		#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
		return( std::distance( a_iterator.begin(), a_iterator.end() ) );
		#else
		return( std::distance( begin(a_iterator, boost::numeric::ublas::iterator1_tag()),
				end(a_iterator, boost::numeric::ublas::iterator1_tag()) ) );
		#endif
		}
		
		static inline void GraphNeighbors( IndexType i, TMatrixType& A, std::vector<IndexType>& neighbors)
		{
		neighbors.clear();
		typename MatrixType::iterator1 a_iterator = A.begin1();
		std::advance(a_iterator,i);
		#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
		for (typename MatrixType::iterator2 row_iterator = a_iterator.begin() ;
		row_iterator != a_iterator.end() ; ++row_iterator) 
		{
		#else
		for ( typename MatrixType::iterator2 row_iterator = begin(a_iterator,
			boost::numeric::ublas::iterator1_tag()); 
		row_iterator != end(a_iterator, 
					boost::numeric::ublas::iterator1_tag()); ++row_iterator )
		{
		#endif
			neighbors.push_back( row_iterator.index2() );
		}
		}
	
	
		//********************************************************************
		//checks if a multiplication is needed and tries to do otherwise
		static void InplaceMult(VectorType& rX, const double A)
		{
	
		if( A == 1.00) 
		{}
		else if( A == -1.00)
		{
			//rX *= A;
			typename VectorType::iterator x_iterator = rX.begin();
			typename VectorType::iterator end_iterator = rX.end();
			while(x_iterator != end_iterator)
			{
				*x_iterator = -*x_iterator;
				x_iterator++;
	
			}
		}
		else
		{
			rX *= A;
		}
		}
	
		//********************************************************************
		//checks if a multiplication is needed and tries to do otherwise
		//ATTENTION it is assumed no aliasing between rX and rY
		// X = A*y;
		static void Assign(VectorType& rX, const double A, const VectorType& rY)
		{
		if( A == 1.00) 
			noalias(rX) = rY;
		else if( A == -1.00)
			noalias(rX) = -rY;
		else
			noalias(rX) = A*rY;
		}
	
		//********************************************************************
		//checks if a multiplication is needed and tries to do otherwise
		//ATTENTION it is assumed no aliasing between rX and rY
		// X += A*y;
		static void UnaliasedAdd(VectorType& rX, const double A, const VectorType& rY)
		{
		if( A == 1.00) 
			noalias(rX) += rY;
		else if( A == -1.00)
			noalias(rX) -= rY;
		else
			noalias(rX) += A*rY;
		}
	
		//********************************************************************
	static void ScaleAndAdd(const double A, const VectorType& rX, const double B, const VectorType& rY, VectorType& rZ)  // rZ = (A * rX) + (B * rY)
		{
			Assign(rZ,A,rX); //rZ = A*rX
			UnaliasedAdd(rZ,B,rY); //rZ += B*rY
	//KRATOS_WATCH(rZ);
		//typename VectorType::const_iterator x_iterator = rX.begin();
		//typename VectorType::const_iterator y_iterator = rY.begin();
		//typename VectorType::iterator z_iterator = rZ.begin();
		//typename VectorType::const_iterator end_iterator = rX.end();
	
		//while(x_iterator != end_iterator)
		//  *z_iterator++ = (A * *x_iterator++) + (B * *y_iterator++);
	
		} 
	
	static void ScaleAndAdd(const double A,const  VectorType& rX, const double B, VectorType& rY) // rY = (A * rX) + (B * rY) 
		{
			InplaceMult(rY,B);
			UnaliasedAdd(rY,A,rX);
	//KRATOS_WATCH(rY);
		//typename VectorType::const_iterator x_iterator = rX.begin();
		//typename VectorType::iterator y_iterator = rY.begin();
		//typename VectorType::const_iterator end_iterator = rX.end();
		//
		//double c = B - double(1.00);
	
		//while(x_iterator != end_iterator)
		//  {
		//    *y_iterator += (A * *x_iterator++) + (c * *y_iterator);
		//    y_iterator++;
		//  }
		} 
	
	
	/// rA[i] * rX
	static double RowDot(unsigned int i, MatrixType& rA, VectorType& rX)
		{
		return inner_prod(row(rA, i), rX);
		} 
	
	
	/// rX = A
	static void Set(VectorType& rX, TDataType A)
		{std::fill(rX.begin(), rX.end(), A);} 
	
	static void Resize(MatrixType& rA, SizeType m, SizeType n){rA.resize(m,n,false);}
	
	static void Resize(VectorType& rX, SizeType n) {rX.resize(n,false);}
	
	static void Clear(MatrixPointerType& pA){pA->clear();}
	static void Clear(VectorPointerType& pX){pX->clear();} 
/*	static void Clear(MatrixType& rA)
		{rA.clear();}
	
	static void Clear(VectorType& rX) {rX.clear();}*/
	
	template<class TOtherMatrixType>
	inline static void ClearData(TOtherMatrixType& rA){rA.clear();}
	
	inline static void ClearData(compressed_matrix<TDataType>& rA)
	{
		rA.clear();
	//    	rA.value_data() = unbounded_array<TDataType>();
		//if(rA.non_zeros() != 0) rA.value_data() = unbounded_array<TDataType>();
	}
	
	inline static void ClearData(VectorType& rX) {rX = VectorType();}
	
	template<class TOtherMatrixType>
	inline static void ResizeData(TOtherMatrixType& rA, SizeType m){rA.resize(m,false);
	std::fill(rA.begin(), rA.end(), TDataType());}
	
	inline static void ResizeData(compressed_matrix<TDataType>& rA, SizeType m) {rA.value_data().resize(m);
	std::fill(rA.value_data().begin(), rA.value_data().end(), TDataType());}
	
	inline static void ResizeData(VectorType& rX, SizeType m) {rX.resize(m);
	std::fill(rX.begin(), rX.end(), TDataType());}
	
	template<class TOtherMatrixType>
	inline static void SetToZero(TOtherMatrixType& rA){std::fill(rA.begin(), rA.end(), TDataType());}
	
	inline static void SetToZero(compressed_matrix<TDataType>& rA)
		{
			KRATOS_TRY
			typedef  unsigned int size_type;
			typedef  double value_type;
	
			size_type begin = rA.index1_data () [0];
			size_type end = rA.index1_data () [rA.size1()];
	
			for (size_type j = begin; j < end; ++ j)
			{
				rA.value_data()[j] = TDataType();
			}
			KRATOS_CATCH("");		
		}
		//{std::fill(rA.value_data().begin(), rA.value_data().end(), TDataType());}
	
	inline static void SetToZero(VectorType& rX) {std::fill(rX.begin(), rX.end(), TDataType());}
	
	template<class TOtherMatrixType, class TEquationIdVectorType>
	inline static void AssembleLHS(
		MatrixType& A,
		TOtherMatrixType& LHS_Contribution,
		TEquationIdVectorType& EquationId
		)
		{
		unsigned int system_size = A.size1();
		unsigned int local_size = LHS_Contribution.size1();
		
		for (unsigned int i_local=0; i_local<local_size; i_local++)
		{
			unsigned int i_global=EquationId[i_local];
			if ( i_global < system_size )
			{
			for (unsigned int j_local=0; j_local<local_size; j_local++)
			{
				unsigned int j_global=EquationId[j_local];
				if ( j_global < system_size )
				A(i_global,j_global) += LHS_Contribution(i_local,j_local);
			}
			}
		}
		}


//        static void GatherLocalValues(Vector& global_indices, )
//		{
//		axpy_prod(rA, rX, rY, true);
//		}
	
	
	
	
	///@}
	///@name Access
	///@{ 
	
	
	///@}
	///@name Inquiry
	///@{
	
	
	///@}      
	///@name Input and output
	///@{
	
	/// Turn back information as a string.
	virtual std::string Info() const
		{
		return "UBlasSpace";
		}
	
	/// Print information about this object.
	virtual void PrintInfo(std::ostream& rOStream) const
		{
		rOStream << "UBlasSpace";
		}
	
	/// Print object's data.
	virtual void PrintData(std::ostream& rOStream) const
		{
		}
	
        //***********************************************************************
        inline static bool IsDistributed() {
            return false;
        }

        //***********************************************************************

        inline static double GetValue(const VectorType& x, std::size_t I) {
            return x[I];
        }
        //***********************************************************************

        static void GatherValues(const VectorType& x, const std::vector<std::size_t>& IndexArray, double* pValues) {
            KRATOS_TRY

            for(std::size_t i = 0; i<IndexArray.size(); i++)
                pValues[i] = x[IndexArray[i]];

            KRATOS_CATCH("")
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
	
	/// Assignment operator.
	UblasSpace& operator=(UblasSpace const& rOther);
	
	/// Copy constructor.
	UblasSpace(UblasSpace const& rOther);
	
		
	///@}    
		
}; // Class UblasSpace 



///@} 

///@name Type Definitions       
///@{ 


///@} 
///@name Input and output 
///@{ 
	

/// input stream function
//   inline std::istream& operator >> (std::istream& rIStream, 
// 				    UblasSpace& rThis);

//   /// output stream function
//   inline std::ostream& operator << (std::ostream& rOStream, 
// 				    const UblasSpace& rThis)
//     {
//       rThis.PrintInfo(rOStream);
//       rOStream << std::endl;
//       rThis.PrintData(rOStream);

//       return rOStream;
//     }
///@} 


}  // namespace Kratos.

#endif // KRATOS_UBLAS_SPACE_H_INCLUDED  defined 


