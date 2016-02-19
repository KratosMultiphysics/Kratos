//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                    
//
#if !defined(KRATOS_PARALLEL_UBLAS_SPACE_H_INCLUDED )
#define  KRATOS_PARALLEL_UBLAS_SPACE_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <cstddef>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "includes/matrix_market_interface.h"
#include "spaces/ublas_space.h"

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
    MultValue ( const Type& _Val ) : Factor ( _Val )
    {
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

    typedef UblasSpace<TDataType, TMatrixType, TVectorType> UblasSpaceType;

    typedef typename boost::shared_ptr< TMatrixType > MatrixPointerType;
    typedef typename boost::shared_ptr< TVectorType > VectorPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ParallelUblasSpace() {}

    /// Destructor.
    virtual ~ParallelUblasSpace() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{
    static MatrixPointerType CreateEmptyMatrixPointer()
    {
        return MatrixPointerType( new TMatrixType() );
    }
    static VectorPointerType CreateEmptyVectorPointer()
    {
        return VectorPointerType( new TVectorType() );
    }

    /// return size of vector rV
    static IndexType Size(VectorType const& rV)
    {
        return rV.size();
    }

    /// return number of rows of rM
    static IndexType Size1(MatrixType const& rM)
    {
        return rM.size1();
    }

    /// return number of columns of rM
    static IndexType Size2(MatrixType const& rM)
    {
        return rM.size2();
    }

    /// rXi = rMij
    static void GetColumn(unsigned int j, MatrixType& rM, VectorType& rX)
    {
        rX = column(rM, j);
    }

    /// rMij = rXi
    static void SetColumn(unsigned int j, MatrixType& rM, VectorType& rX)
    {
        rX = row(rM, j);
    }

    /// rY = rX
    static void Copy(MatrixType const& rX, MatrixType& rY)
    {
// :TODO: Parallelize
        rY.assign(rX);
    }

    /// rY = rX
    static void Copy(VectorType const& rX, VectorType& rY)
    {
//             rY.assign(rX);
        UblasSpaceType::Copy(rX,rY);
    }

    /// rX * rY
    static TDataType Dot(VectorType const& rX, VectorType const& rY)
    {
        vector<unsigned int> partition;
#ifdef _OPENMP
        int number_of_threads = omp_get_max_threads();
#else
        int number_of_threads = 1;
#endif
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
        for(int j = 0; j<number_of_threads; j++)
            total += partial_results[j];

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

        UblasSpaceType::InplaceMult(rX,A);

    }


    //********************************************************************
    //checks if a multiplication is needed and tries to do otherwise
    //ATTENTION it is assumed no aliasing between rX and rY
    // X = A*y;
    static void Assign(VectorType& rX, const double A, const VectorType& rY)
    {

        UblasSpaceType::Assign(rX,A,rY);

    }

    //********************************************************************
    //checks if a multiplication is needed and tries to do otherwise
    //ATTENTION it is assumed no aliasing between rX and rY
    // X += A*y;
    static void UnaliasedAdd(VectorType& rX, const double A, const VectorType& rY)
    {
        UblasSpaceType::UnaliasedAdd(rX, A, rY);
    }

    //********************************************************************
    static void ScaleAndAdd(const double A, const VectorType& rX, const double B, const VectorType& rY, VectorType& rZ)  // rZ = (A * rX) + (B * rY)
    {
        Assign(rZ,A,rX); //rZ = A*rX
        UnaliasedAdd(rZ,B,rY); //rZ += B*rY
    }

    static void ScaleAndAdd(const double A,const  VectorType& rX, const double B, VectorType& rY) // rY = (A * rX) + (B * rY)
    {
        InplaceMult(rY,B);
        UnaliasedAdd(rY,A,rX);
    }


    /// rA[i] * rX
    //will be most probably faster in serial as the rows are short
    static double RowDot(unsigned int i, MatrixType& rA, VectorType& rX)
    {
        return inner_prod(row(rA, i), rX);
    }


    /// rX = A
    static void Set(VectorType& rX, TDataType A)
    {
        UblasSpaceType::set(rX, A);

    }

    static void Resize(MatrixType& rA, SizeType m, SizeType n)
    {
        rA.resize(m,n,false);
    }

    static void Resize(VectorType& rX, SizeType n)
    {
        rX.resize(n,false);
    }

    static void Clear(MatrixPointerType& pA)
    {
        pA->clear();
    }
    static void Clear(VectorPointerType& pX)
    {
        pX->clear();
    }

    template<class TOtherMatrixType>
    inline static void ClearData(TOtherMatrixType& rA)
    {
        rA.clear();
    }

    inline static void ClearData(compressed_matrix<TDataType>& rA)
    {
        UblasSpaceType::ClearData(rA);

    }

    inline static void ClearData(VectorType& rX)
    {
        rX = VectorType();
    }


    inline static void ResizeData(compressed_matrix<TDataType>& rA, SizeType m)
    {
        UblasSpaceType::ResizeData(rA,m);
    }

    inline static void ResizeData(VectorType& rX, SizeType m)
    {
        UblasSpaceType::ResizeData(rX,m);
    }

    template<class TOtherMatrixType>
    inline static void SetToZero(TOtherMatrixType& rA)
    {
        std::fill(rA.begin(), rA.end(), TDataType());
    }

    inline static void SetToZero(compressed_matrix<TDataType>& rA)
    {
        KRATOS_TRY
        UblasSpaceType::SetToZero(rA);
        KRATOS_CATCH("");
    }

    inline static void SetToZero(VectorType& rX)
    {
        UblasSpaceType::SetToZero(rX);
    }

    //***********************************************************************
    inline static bool IsDistributed()
    {
        return false;
    }

    //***********************************************************************

    inline static double GetValue(const VectorType& x, std::size_t I)
    {
        return x[I];
    }
    //***********************************************************************

    static void GatherValues(const VectorType& x, const std::vector<std::size_t>& IndexArray, double* pValues)
    {
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

    template< class TOtherMatrixType >
    static bool WriteMatrixMarketMatrix(const char *FileName, TOtherMatrixType &M, bool Symmetric)
    {
        // Use full namespace in call to make sure we are not calling this function recursively
        return Kratos::WriteMatrixMarketMatrix(FileName,M,Symmetric);
    }
    
        template< class VectorType >
    static bool WriteMatrixMarketVector(const char *FileName, VectorType& V)
    {
        // Use full namespace in call to make sure we are not calling this function recursively
        return Kratos::WriteMatrixMarketVector(FileName,V);
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
//             size_type  processor_row_begin, processor_row_end;
//             int proc_id = 0;

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
        for(unsigned int i = 1; i<number_of_threads; i++)
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

