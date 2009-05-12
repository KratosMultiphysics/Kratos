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
//   Revision:            $Revision: 1.5 $
//
//


#if !defined(KRATOS_PARALLEL_UBLAS_SPACE_H_INCLUDED )
#define  KRATOS_PARALLEL_UBLAS_SPACE_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <cstddef>
#include "omptl"
#include "omptl_algorithm"
#include "omptl_numeric"

// External includes 


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

    // The function object multiplies an element by a Factor
    template <class Type>
    class MultValue
    {
        private:
            Type Factor;   // The value to multiply by
        public:
            // Constructor initializes the value to multiply by
            MultValue ( const Type& _Val ) : Factor ( _Val ) {
            }

            // The function call for the element to be multiplied
            Type operator ( ) ( Type& elem ) const 
            {
                return elem * Factor;
            }
    };


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
  class ParallelUblasSpace
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of ParallelUblasSpace
      KRATOS_CLASS_POINTER_DEFINITION(ParallelUblasSpace);
  
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
      ParallelUblasSpace(){}

      /// Destructor.
      virtual ~ParallelUblasSpace(){}
      

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

      /// rMij = rXi
      static void SetColumn(unsigned int j, MatrixType& rM, VectorType& rX){rX = row(rM, j);} 

      /// rY = rX
      static void Copy(MatrixType const& rX, MatrixType& rY)
        {
// :TODO: Parallelize
            rY.assign(rX);
//             omptl::copy( rX.begin(), rX.end(), rY.begin() );
        } 

      /// rY = rX
      static void Copy(VectorType const& rX, VectorType& rY)
        {
//             rY.assign(rX);
            omptl::copy( rX.begin(), rX.end(), rY.begin() );
        } 

      /// rX * rY
      static TDataType Dot(VectorType const& rX, VectorType const& rY)
        {
            vector<unsigned int> partition;
            int number_of_threads = omp_get_max_threads();
            CreatePartition(number_of_threads, rX.size(), partition);

	    vector< TDataType > partial_results(number_of_threads);

	  int i;
	  #pragma omp parallel for default(shared) private(i)
	   for(i = 0; i<number_of_threads; i++)
	   {
		partial_results[i] = std::inner_product( rX.data().begin()+partition[i],
							 rX.data().begin()+partition[i+1], 
							 rY.data().begin()+partition[i],
							 TDataType() );
	   }

	   double total = TDataType();
	   for(int i = 0; i<number_of_threads; i++)
	    total += partial_results[i];
	
//              return inner_prod(rX, rY);
           return total;
        } 

      /// ||rX||2
      static double TwoNorm(VectorType const& rX)
        {
            return sqrt( Dot( rX, rX ) );
//             return norm_2(rX);
        } 
	

      static void Mult(MatrixType& rA, VectorType& rX, VectorType& rY)
	  {
          	ParallelProductNoAdd( rA, rX, rY );
// 	  axpy_prod(rA, rX, rY, true);
	  }// rY = rA * rX

      static void TransposeMult(MatrixType& rA, VectorType& rX, VectorType& rY)
	  {
// :TODO: Parallelize
	      axpy_prod(rX, rA, rY, true);
	  }// rY = rAT * rX
 
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
         omptl::transform( rX.begin(), rX.end(), rX.begin(), std::negate<double>() );
		  //rX *= A;
// 		  typename VectorType::iterator x_iterator = rX.begin();
// 		  typename VectorType::iterator end_iterator = rX.end();
// 		  while(x_iterator != end_iterator)
// 		  {
// 			  *x_iterator = -*x_iterator;
// 			  x_iterator++;
// 
// 		  }
	  }
	  else
	  {
// 		  rX *= A;
         omptl::transform( rX.begin(), rX.end(), rX.begin(), MultValue<double>( A ) );
	  }
	}


	//********************************************************************
	//checks if a multiplication is needed and tries to do otherwise
	//ATTENTION it is assumed no aliasing between rX and rY
	// X = A*y;
	static void Assign(VectorType& rX, const double A, const VectorType& rY)
	{
	  if( A == 1.00) 
        omptl::copy( rY.begin(), rY.end(), rX.begin() ); 
// 		  noalias(rX) = rY;
	  else if( A == -1.00)
      {
        omptl::transform( rY.begin(), rY.end(), rX.begin(), std::negate<double>() );
      }
// 		  noalias(rX) = -rY;
	  else
      {
         omptl::transform( rY.begin(), rY.end(), rX.begin(), MultValue<double>( A ) );   
//TODO .. parallelize
// 		  noalias(rX) = A*rY;
      }
	}

	//********************************************************************
	//checks if a multiplication is needed and tries to do otherwise
	//ATTENTION it is assumed no aliasing between rX and rY
	// X += A*y;
	static void UnaliasedAdd(VectorType& rX, const double A, const VectorType& rY)
	{
	  if( A == 1.00) 
      {
// 		  noalias(rX) += rY;
            omptl::transform( rY.data().begin(), rY.data().end(), rX.data().begin(), rX.data().begin(), std::plus<double>() );
      }
	  else if( A == -1.00)
      {
 		  noalias(rX) -= rY;
            omptl::transform( rY.data().begin(), rY.data().end(), rX.data().begin(), rX.data().begin(), std::minus<double>() );
//        omptl::transform( rY.data().begin(), rY.data().end(), rX.data().begin(), std::minus<double>() );
      }
	  else
      {
//TODO: parallelize!!!
		  noalias(rX) += A*rY;
      }
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
        //will be most probably faster in serial as the rows are short
      static double RowDot(unsigned int i, MatrixType& rA, VectorType& rX)
	{
//         return omptl::inner_product( row(rA, i).begin(), row(rA, i).end(), rX.begin() );
        return inner_prod(row(rA, i), rX);
	}

      
      /// rX = A
      static void Set(VectorType& rX, TDataType A)
      {
          //std::fill(rX.begin(), rX.end(), A);
          omptl::fill(rX.begin(), rX.end(), A );
      } 

      static void Resize(MatrixType& rA, SizeType m, SizeType n){rA.resize(m,n,false);}

      static void Resize(VectorType& rX, SizeType n) {rX.resize(n,false);}
      
	static void Clear(MatrixPointerType& pA){pA->clear();}
	static void Clear(VectorPointerType& pX){pX->clear();} 
/*      static void Clear(MatrixType& rA)
	{rA.clear();}

      static void Clear(VectorType& rX) {rX.clear();}*/
      
      template<class TOtherMatrixType>
      inline static void ClearData(TOtherMatrixType& rA){rA.clear();}

      inline static void ClearData(compressed_matrix<TDataType>& rA)
      {
	        //rA.clear();
            omptl::fill( rA.value_data().begin(), rA.value_data().end(), 0.0 );
      }

      inline static void ClearData(VectorType& rX) {rX = VectorType();}
      
//      template<class TOtherMatrixType>
  //    inline static void ResizeData(TOtherMatrixType& rA, SizeType m){rA.resize(m,false);
    //  std::fill(rA.begin(), rA.end(), TDataType());}

      inline static void ResizeData(compressed_matrix<TDataType>& rA, SizeType m) {rA.value_data().resize(m);
      omptl::fill(rA.value_data().begin(), rA.value_data().end(), TDataType());}

      inline static void ResizeData(VectorType& rX, SizeType m) {rX.resize(m);
      omptl::fill(rX.begin(), rX.end(), TDataType());}
      
      template<class TOtherMatrixType>
      inline static void SetToZero(TOtherMatrixType& rA)
	{std::fill(rA.begin(), rA.end(), TDataType());}

      inline static void SetToZero(compressed_matrix<TDataType>& rA)
	{
		KRATOS_TRY
	omptl::fill(rA.value_data().begin(), rA.value_data().end(), TDataType());
/*		typedef  unsigned int size_type;
		typedef  double value_type;

		size_type begin = rA.index1_data () [0];
		size_type end = rA.index1_data () [rA.size1()];

		for (size_type j = begin; j < end; ++ j)
		{
			rA.value_data()[j] = TDataType();
		}*/
		KRATOS_CATCH("");		
	}

      inline static void SetToZero(VectorType& rX) {omptl::fill(rX.begin(), rX.end(), TDataType());}

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
	return "ParallelUblasSpace";
	}
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
	{
	rOStream << "ParallelUblasSpace";
	}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const
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
        //y += A*x in parallel
        static void ParallelProductNoAdd( const MatrixType& A, const VectorType& in, VectorType& out)
        {
//std::cout << "in function ParallelProductNoAdd" << std::endl;
            typedef  unsigned int size_type;
            typedef  double value_type;
        
            //create partition
            vector<size_type> partition;
            int number_of_threads = omp_get_max_threads();
            CreatePartition(number_of_threads, A.size1(), partition);
            //parallel loop
            size_type  processor_row_begin, processor_row_end;
            int proc_id = 0;
        
            #pragma omp parallel
            {
                int thread_id = omp_get_thread_num();
                int number_of_rows = partition[thread_id+1] - partition[thread_id];
                typename MatrixType::index_array_type::const_iterator row_iter_begin = A.index1_data().begin()+partition[thread_id];
                typename MatrixType::index_array_type::const_iterator index_2_begin = A.index2_data().begin()+*row_iter_begin;
                typename MatrixType::value_array_type::const_iterator value_begin = A.value_data().begin()+*row_iter_begin;
                typename VectorType::iterator output_vec_begin = out.begin()+partition[thread_id];
                
            
                partial_product_no_add(    number_of_rows,
                                    row_iter_begin,
                                    index_2_begin,
                                    value_begin,
                                    in,
                                    output_vec_begin
                                );
            }
        }

        static void CreatePartition(unsigned int number_of_threads,const int number_of_rows, vector<unsigned int>& partitions)
        {
            partitions.resize(number_of_threads+1);
            int partition_size = number_of_rows / number_of_threads;
            partitions[0] = 0;
            partitions[number_of_threads] = number_of_rows;
            for(int i = 1; i<number_of_threads; i++)
               partitions[i] = partitions[i-1] + partition_size ;
        }
        
	/**
	 * calculates partial product resetting to Zero the output before
	 */
        static void partial_product_no_add(
                int number_of_rows,
                typename TMatrixType::index_array_type::const_iterator row_begin,
                typename TMatrixType::index_array_type::const_iterator index2_begin,
                typename TMatrixType::value_array_type::const_iterator value_begin,
                const VectorType& input_vec,
                typename VectorType::iterator output_vec_begin
                 )
        {
            int row_size;
            typename MatrixType::index_array_type::const_iterator row_it = row_begin;
            for(int k = 0; k < number_of_rows; k++)
            {
                row_size= *(row_it+1)-*row_it;
                row_it++;
                TDataType t = TDataType();
                
                for(int i = 0; i<row_size; i++)
                    t += *value_begin++ *   ( input_vec[*index2_begin++]);
        
                *output_vec_begin++ = t;
                    
            }
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
      ParallelUblasSpace& operator=(ParallelUblasSpace const& rOther);

      /// Copy constructor.
      ParallelUblasSpace(ParallelUblasSpace const& rOther);

        
      ///@}    
        
    }; // Class ParallelUblasSpace 



  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
//   inline std::istream& operator >> (std::istream& rIStream, 
// 				    ParallelUblasSpace& rThis);

//   /// output stream function
//   inline std::ostream& operator << (std::ostream& rOStream, 
// 				    const ParallelUblasSpace& rThis)
//     {
//       rThis.PrintInfo(rOStream);
//       rOStream << std::endl;
//       rThis.PrintData(rOStream);

//       return rOStream;
//     }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_PARALLEL_UBLAS_SPACE_H_INCLUDED  defined 

