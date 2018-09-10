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
//  Collaborator:    Vicente Mataix Ferrandiz
//
//

#if !defined(KRATOS_UBLAS_SPACE_H_INCLUDED )
#define  KRATOS_UBLAS_SPACE_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <cstddef>
#include <numeric>

#ifdef _OPENMP
#include "omp.h"
#endif


// External includes


// Project includes
#include "includes/define.h"
#include "includes/ublas_interface.h"
#include "includes/matrix_market_interface.h"
#include "utilities/dof_updater.h"

namespace Kratos
{
#ifdef _OPENMP
// The function object multiplies an element by a Factor

template <class Type>
class MultValueNoAdd
{
private:
    Type Factor; // The value to multiply by
public:
    // Constructor initializes the value to multiply by

    MultValueNoAdd(const Type& _Val) : Factor(_Val)
    {
    }

    // The function call for the element to be multiplied

    inline Type operator () (const Type& elem) const
    {
        return elem * Factor;
    }
};

template <class Type>
class MultAndAddValue : public std::binary_function<Type, Type, Type>
{
private:
    Type Factor; // The value to multiply by
public:
    // Constructor initializes the value to multiply by

    MultAndAddValue(const Type& _Val) : Factor(_Val)
    {
    }

    // The function call for the element to be multiplied

    inline Type operator () (const Type& elem1, const Type& elem2) const
    {
        return elem1 * Factor + elem2;
    }
};
#endif

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

    typedef typename Kratos::shared_ptr< TMatrixType > MatrixPointerType;
    typedef typename Kratos::shared_ptr< TVectorType > VectorPointerType;

#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it 
    template<typename T> using compressed_matrix = boost::numeric::ublas::compressed_matrix<T>;
#endif // ifdef KRATOS_USE_AMATRIX

    typedef DofUpdater< UblasSpace<TDataType,TMatrixType,TVectorType> > DofUpdaterType;
    typedef typename DofUpdaterType::UniquePointer DofUpdaterPointerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    UblasSpace()
    {
    }

    /// Destructor.

    virtual ~UblasSpace()
    {
    }


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    static MatrixPointerType CreateEmptyMatrixPointer()
    {
        return MatrixPointerType(new TMatrixType(0, 0));
    }

    static VectorPointerType CreateEmptyVectorPointer()
    {
        return VectorPointerType(new TVectorType(0));
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
	// This version is needed in order to take one column of multi column solve from AMatrix matrix and pass it to an ublas vector
	template<typename TColumnType>
	static void GetColumn(unsigned int j, Matrix& rM, TColumnType& rX)
	{
		if (rX.size() != rM.size1())
			rX.resize(rM.size1(), false);

		for (std::size_t i = 0; i < rM.size1(); i++) {
			rX[i] = rM(i, j);
		}
	}

	// This version is needed in order to take one column of multi column solve from AMatrix matrix and pass it to an ublas vector
	template<typename TColumnType>
	static void SetColumn(unsigned int j, Matrix& rM, TColumnType& rX)
	{
		for (std::size_t i = 0; i < rM.size1(); i++) {
			rM(i,j) = rX[i];
		}
	}


    /// rY = rX

    static void Copy(MatrixType const& rX, MatrixType& rY)
    {
        rY.assign(rX);
    }

    /// rY = rX

    static void Copy(VectorType const& rX, VectorType& rY)
    {
#ifndef _OPENMP
        rY.assign(rX);
#else

        const int size = rX.size();
        if (rY.size() != static_cast<unsigned int>(size))
            rY.resize(size, false);

        #pragma omp parallel for
        for (int i = 0; i < size; i++)
            rY[i] = rX[i];
#endif
    }

    /// rX * rY

    static TDataType Dot(VectorType const& rX, VectorType const& rY)
    {
#ifndef _OPENMP
        return inner_prod(rX, rY);
#else
        const int size = static_cast<int>(rX.size());

        TDataType total = TDataType();
        #pragma omp parallel for reduction( +: total), firstprivate(size)
        for(int i =0; i<size; ++i)
            total += rX[i]*rY[i];

        return total;
#endif
    }


    /// ||rX||2

    static TDataType TwoNorm(VectorType const& rX)
    {
        return std::sqrt(Dot(rX, rX));
    }

    static TDataType TwoNorm(MatrixType const& rA) // Frobenious norm
    {
        TDataType aux_sum = TDataType();
#ifndef _OPENMP
        for (int i = 0; i < static_cast<int>(rA.size1()); i++)
        {
            for (int j = 0; j < static_cast<int>(rA.size2()); j++)
            {
                aux_sum += rA(i,j) * rA(i,j);
            }
        }
#else
        #pragma omp parallel reduction(+:aux_sum)
        for (int i = 0; i < static_cast<int>(rA.size1()); i++)
        {
            for (int j = 0; j < static_cast<int>(rA.size2()); j++)
            {
                aux_sum += rA(i,j) * rA(i,j);
            }
        }
#endif
        return std::sqrt(aux_sum);
    }

    /**
     * This method computes the Jacobi norm
     * @param rA The matrix to compute the Jacobi norm
     * @return aux_sum: The Jacobi norm
     */
    static TDataType JacobiNorm(MatrixType const& rA)
    {
        TDataType aux_sum = TDataType();

#ifndef _OPENMP
        for (int i = 0; i < static_cast<int>(rA.size1()); i++)
        {
            for (int j = 0; j < static_cast<int>(rA.size2()); j++)
            {
                if (i != j)
                {
                    aux_sum += std::abs(rA(i,j));
                }
            }
        }
#else
        #pragma omp parallel for reduction(+:aux_sum)
        for (int i = 0; i < static_cast<int>(rA.size1()); i++)
        {
            for (int j = 0; j < static_cast<int>(rA.size2()); j++)
            {
                if (i != j)
                {
                    aux_sum += std::abs(rA(i,j));
                }
            }
        }
#endif

        return aux_sum;
    }

    static void Mult(const Matrix& rA, VectorType& rX, VectorType& rY)
    {
        axpy_prod(rA, rX, rY, true);
    }

    static void Mult(const compressed_matrix<TDataType>& rA, VectorType& rX, VectorType& rY)
    {
#ifndef _OPENMP
        axpy_prod(rA, rX, rY, true);
#else
        ParallelProductNoAdd(rA, rX, rY);
#endif
    }

    static void TransposeMult(MatrixType& rA, VectorType& rX, VectorType& rY)
    {
		boost::numeric::ublas::axpy_prod(rX, rA, rY, true);
    } // rY = rAT * rX

    static inline SizeType GraphDegree(IndexType i, TMatrixType& A)
    {
        typename MatrixType::iterator1 a_iterator = A.begin1();
        std::advance(a_iterator, i);
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
        return ( std::distance(a_iterator.begin(), a_iterator.end()));
#else
        return ( std::distance(begin(a_iterator, boost::numeric::ublas::iterator1_tag()),
                               end(a_iterator, boost::numeric::ublas::iterator1_tag())));
#endif
    }

    static inline void GraphNeighbors(IndexType i, TMatrixType& A, std::vector<IndexType>& neighbors)
    {
        neighbors.clear();
        typename MatrixType::iterator1 a_iterator = A.begin1();
        std::advance(a_iterator, i);
#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
        for (typename MatrixType::iterator2 row_iterator = a_iterator.begin();
                row_iterator != a_iterator.end(); ++row_iterator)
        {
#else
        for (typename MatrixType::iterator2 row_iterator = begin(a_iterator,
                boost::numeric::ublas::iterator1_tag());
                row_iterator != end(a_iterator,
                                    boost::numeric::ublas::iterator1_tag()); ++row_iterator)
        {
#endif
            neighbors.push_back(row_iterator.index2());
        }
    }


    //********************************************************************
    //checks if a multiplication is needed and tries to do otherwise

    static void InplaceMult(VectorType& rX, const double A)
    {

        if (A == 1.00)
        {
        }
        else if (A == -1.00)
        {
#ifndef _OPENMP
            typename VectorType::iterator x_iterator = rX.begin();
            typename VectorType::iterator end_iterator = rX.end();
            while (x_iterator != end_iterator)
            {
                *x_iterator = -*x_iterator;
                x_iterator++;

            }
#else
             const int size = rX.size();

            #pragma omp parallel for firstprivate(size)
            for (int i = 0; i < size; i++)
                rX[i] = -rX[i];

#endif
        }
        else
        {
#ifndef _OPENMP
            rX *= A;
#else
            const int size = rX.size();

            #pragma omp parallel for firstprivate(size)
            for (int i = 0; i < size; i++)
                rX[i] *= A;
#endif
        }
    }

    //********************************************************************
    //checks if a multiplication is needed and tries to do otherwise
    //ATTENTION it is assumed no aliasing between rX and rY
    // X = A*y;

    static void Assign(VectorType& rX, const double A, const VectorType& rY)
    {
#ifndef _OPENMP
        if (A == 1.00)
            noalias(rX) = rY;
        else if (A == -1.00)
            noalias(rX) = -rY;
        else
            noalias(rX) = A*rY;
#else
        const int size = rY.size();
        if (rX.size() != static_cast<unsigned int>(size) )
            rX.resize(size, false);

        if (A == 1.00)
        {
            #pragma omp parallel for
            for (int i = 0; i < size; i++)
                rX[i] = rY[i];
        }
        else if (A == -1.00)
        {
            #pragma omp parallel for
            for (int i = 0; i < size; i++)
                rX[i] = -rY[i];
        }
        else
        {
            #pragma omp parallel for
            for (int i = 0; i < size; i++)
                rX[i] = A * rY[i];
        }
#endif
    }

    //********************************************************************
    //checks if a multiplication is needed and tries to do otherwise
    //ATTENTION it is assumed no aliasing between rX and rY
    // X += A*y;

    static void UnaliasedAdd(VectorType& rX, const double A, const VectorType& rY)
    {
#ifndef _OPENMP
        if (A == 1.00)
            noalias(rX) += rY;
        else if (A == -1.00)
            noalias(rX) -= rY;
        else
            noalias(rX) += A*rY;
#else
        const int size = rY.size();
        if (rX.size() != static_cast<unsigned int>(size) )
            rX.resize(size, false);

        if (A == 1.00)
        {
            #pragma omp parallel for
            for (int i = 0; i < size; i++)
                rX[i] += rY[i];
        }
        else if (A == -1.00)
        {
            #pragma omp parallel for
            for (int i = 0; i < size; i++)
                rX[i] -= rY[i];
        }
        else
        {
            #pragma omp parallel for
            for (int i = 0; i < size; i++)
                rX[i] += A * rY[i];
        }


#endif

    }

    //********************************************************************

    static void ScaleAndAdd(const double A, const VectorType& rX, const double B, const VectorType& rY, VectorType& rZ) // rZ = (A * rX) + (B * rY)
    {
        Assign(rZ, A, rX); //rZ = A*rX
        UnaliasedAdd(rZ, B, rY); //rZ += B*rY
    }

    static void ScaleAndAdd(const double A, const VectorType& rX, const double B, VectorType& rY) // rY = (A * rX) + (B * rY)
    {
        InplaceMult(rY, B);
        UnaliasedAdd(rY, A, rX);
    }


    /// rA[i] * rX

    static double RowDot(unsigned int i, MatrixType& rA, VectorType& rX)
    {
        return inner_prod(row(rA, i), rX);
    }


    static void SetValue(VectorType& rX, IndexType i, TDataType value)
    {
        rX[i] = value;
    }

    /// rX = A

    static void Set(VectorType& rX, TDataType A)
    {
        std::fill(rX.begin(), rX.end(), A);
    }

    static void Resize(MatrixType& rA, SizeType m, SizeType n)
    {
        rA.resize(m, n, false);
    }

    static void Resize(MatrixPointerType& pA, SizeType m, SizeType n)
    {
        pA->resize(m, n, false);
    }

    static void Resize(VectorType& rX, SizeType n)
    {
        rX.resize(n, false);
    }

    static void Resize(VectorPointerType& pX, SizeType n)
    {
        pX->resize(n, false);
    }

    static void Clear(MatrixPointerType& pA)
    {
        pA->clear();
        pA->resize(0, 0, false);
    }

    static void Clear(VectorPointerType& pX)
    {
        pX->clear();
        pX->resize(0, false);
    }

    template<class TOtherMatrixType>
    inline static void ResizeData(TOtherMatrixType& rA, SizeType m)
    {
        rA.resize(m, false);
        //            std::fill(rA.begin(), rA.end(), TDataType());
#ifndef _OPENMP
        std::fill(rA.begin(), rA.end(), TDataType());
#else
        DataType* vals = rA.value_data().begin();
        #pragma omp parallel for firstprivate(m)
        for(int i=0; i<static_cast<int>(m); ++i)
            vals[i] = TDataType();
#endif
    }

    inline static void ResizeData(compressed_matrix<TDataType>& rA, SizeType m)
    {
        rA.value_data().resize(m);
#ifndef _OPENMP
        std::fill(rA.value_data().begin(), rA.value_data().end(), TDataType());
#else
        TDataType* vals = rA.value_data().begin();
        #pragma omp parallel for firstprivate(m)
        for(int i=0; i<static_cast<int>(m); ++i)
            vals[i] = TDataType();
#endif
    }

    inline static void ResizeData(VectorType& rX, SizeType m)
    {
        rX.resize(m, false);
#ifndef _OPENMP
        std::fill(rX.begin(), rX.end(), TDataType());
#else
        const int size = rX.size();
        #pragma omp parallel for firstprivate(size)
        for(int i=0; i<size; ++i)
            rX[i] = TDataType();
#endif
    }

    template<class TOtherMatrixType>
    inline static void SetToZero(TOtherMatrixType& rA)
    {
#ifndef _OPENMP
        std::fill(rA.begin(), rA.end(), TDataType());
#else
        TDataType* vals = rA.value_data().begin();
        const int size = rA.value_data().end() - rA.value_data().begin();
        #pragma omp parallel for firstprivate(size)
        for(int i=0; i<size; ++i)
            vals[i] = TDataType();
#endif
    }

    inline static void SetToZero(compressed_matrix<TDataType>& rA)
    {
#ifndef _OPENMP
        std::fill(rA.value_data().begin(), rA.value_data().end(), TDataType());
#else
        TDataType* vals = rA.value_data().begin();
        const int size = rA.value_data().end() - rA.value_data().begin();
        #pragma omp parallel for firstprivate(size)
        for(int i=0; i<size; ++i)
            vals[i] = TDataType();
#endif
    }

    inline static void SetToZero(VectorType& rX)
    {
#ifndef _OPENMP
        std::fill(rX.begin(), rX.end(), TDataType());
#else
        const int size = rX.size();
        #pragma omp parallel for firstprivate(size)
        for(int i=0; i<size; ++i)
            rX[i] = TDataType();
#endif
    }

    template<class TOtherMatrixType, class TEquationIdVectorType>
    inline static void AssembleLHS(
        MatrixType& A,
        TOtherMatrixType& LHS_Contribution,
        TEquationIdVectorType& EquationId
    )
    {
        unsigned int system_size = A.size1();
        unsigned int local_size = LHS_Contribution.size1();

        for (unsigned int i_local = 0; i_local < local_size; i_local++)
        {
            unsigned int i_global = EquationId[i_local];
            if (i_global < system_size)
            {
                for (unsigned int j_local = 0; j_local < local_size; j_local++)
                {
                    unsigned int j_global = EquationId[j_local];
                    if (j_global < system_size)
                        A(i_global, j_global) += LHS_Contribution(i_local, j_local);
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

    inline static bool IsDistributed()
    {
        return false;
    }

    //***********************************************************************

    inline static TDataType GetValue(const VectorType& x, std::size_t I)
    {
        return x[I];
    }
    //***********************************************************************

    static void GatherValues(const VectorType& x, const std::vector<std::size_t>& IndexArray, TDataType* pValues)
    {
        KRATOS_TRY

        for (std::size_t i = 0; i < IndexArray.size(); i++)
            pValues[i] = x[IndexArray[i]];

        KRATOS_CATCH("")
    }

    template< class TOtherMatrixType >
    static bool WriteMatrixMarketMatrix(const char* pFileName, /*const*/ TOtherMatrixType& rM, const bool Symmetric)
    {
        // Use full namespace in call to make sure we are not calling this function recursively
        return Kratos::WriteMatrixMarketMatrix(pFileName, rM, Symmetric);
    }

    template< class VectorType >
    static bool WriteMatrixMarketVector(const char* pFileName, const VectorType& rV)
    {
        // Use full namespace in call to make sure we are not calling this function recursively
        return Kratos::WriteMatrixMarketVector(pFileName, rV);
    }

    static DofUpdaterPointerType CreateDofUpdater()
    {
        DofUpdaterType tmp;
        return tmp.Create();
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

#ifdef _OPENMP
    //y += A*x in parallel

    static void ParallelProductNoAdd(const MatrixType& A, const VectorType& in, VectorType& out)
    {
        //create partition
        DenseVector<unsigned int> partition;
        unsigned int number_of_threads = omp_get_max_threads();
        unsigned int number_of_initialized_rows = A.filled1() - 1;
        CreatePartition(number_of_threads, number_of_initialized_rows, partition);
        //parallel loop
        #pragma omp parallel
        {
            int thread_id = omp_get_thread_num();
            int number_of_rows = partition[thread_id + 1] - partition[thread_id];
            typename compressed_matrix<TDataType>::index_array_type::const_iterator row_iter_begin = A.index1_data().begin() + partition[thread_id];
            typename compressed_matrix<TDataType>::index_array_type::const_iterator index_2_begin = A.index2_data().begin()+*row_iter_begin;
            typename compressed_matrix<TDataType>::value_array_type::const_iterator value_begin = A.value_data().begin()+*row_iter_begin;
            //                  typename VectorType::iterator output_vec_begin = out.begin()+partition[thread_id];


            partial_product_no_add(number_of_rows,
                                   row_iter_begin,
                                   index_2_begin,
                                   value_begin,
                                   in,
                                   partition[thread_id],
                                   out
                                  );
        }
    }

    static void CreatePartition(unsigned int number_of_threads, const int number_of_rows, DenseVector<unsigned int>& partitions)
    {
        partitions.resize(number_of_threads + 1);
        int partition_size = number_of_rows / number_of_threads;
        partitions[0] = 0;
        partitions[number_of_threads] = number_of_rows;
        for (unsigned int i = 1; i < number_of_threads; i++)
            partitions[i] = partitions[i - 1] + partition_size;
    }


    /**
     * calculates partial product resetting to Zero the output before
     */
    static void partial_product_no_add(
        int number_of_rows,
        typename compressed_matrix<TDataType>::index_array_type::const_iterator row_begin,
        typename compressed_matrix<TDataType>::index_array_type::const_iterator index2_begin,
        typename compressed_matrix<TDataType>::value_array_type::const_iterator value_begin,
        const VectorType& input_vec,
        unsigned int output_begin_index,
        VectorType& output_vec
        //                 typename VectorType::iterator output_vec_begin
    )
    {
        int row_size;
        int kkk = output_begin_index;
        typename MatrixType::index_array_type::const_iterator row_it = row_begin;
        for (int k = 0; k < number_of_rows; k++)
        {
            row_size = *(row_it + 1)-*row_it;
            row_it++;
            TDataType t = TDataType();

            for (int i = 0; i < row_size; i++)
                t += *value_begin++ * (input_vec[*index2_begin++]);

            output_vec[kkk++] = t;
            //                 *output_vec_begin++ = t;

        }
    }
#endif


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
    UblasSpace & operator=(UblasSpace const& rOther);

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


} // namespace Kratos.

#endif // KRATOS_UBLAS_SPACE_H_INCLUDED  defined
