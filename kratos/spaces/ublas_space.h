//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//  Collaborator:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <numeric>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/process_info.h"
#include "includes/ublas_interface.h"
#include "includes/matrix_market_interface.h"
#include "utilities/dof_updater.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

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
class MultAndAddValue
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

template <class TDataType, class TMatrixType, class TVectorType>
class UblasSpace;

template <class TDataType>
using TUblasSparseSpace =
    UblasSpace<TDataType, boost::numeric::ublas::compressed_matrix<TDataType>, boost::numeric::ublas::vector<TDataType>>;
template <class TDataType>
using TUblasDenseSpace =
    UblasSpace<TDataType, DenseMatrix<TDataType>, DenseVector<TDataType>>;

///@}
///@name  Enum's
///@{

// Scaling enum
enum class SCALING_DIAGONAL {NO_SCALING = 0, CONSIDER_NORM_DIAGONAL = 1, CONSIDER_MAX_DIAGONAL = 2, CONSIDER_PRESCRIBED_DIAGONAL = 3};

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class UblasSpace
 * @ingroup KratosCore
 * @brief A class template for handling data types, matrices, and vectors in a Ublas space.
 * @details This class template is designed to work with different data types, matrix types, and vector types
 * within a Ublas space. It provides typedefs and utilities for managing these types effectively.
 * @tparam TDataType The data type used in the Ublas space.
 * @tparam TMatrixType The matrix type used in the Ublas space.
 * @tparam TVectorType The vector type used in the Ublas space.
 * @author Riccardo Rossi
 */
template<class TDataType, class TMatrixType, class TVectorType>
class UblasSpace
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of UblasSpace
    KRATOS_CLASS_POINTER_DEFINITION(UblasSpace);

    /// The data type considered
    using DataType = TDataType;

    /// The matrix type considered
    using MatrixType = TMatrixType;

    /// The vector type considered
    using VectorType = TVectorType;

    /// The index type considered
    using IndexType = std::size_t;

    /// The size type considered
    using SizeType = std::size_t;

    /// The pointer to the matrix type
    using MatrixPointerType = typename Kratos::shared_ptr<TMatrixType>;

    /// The pointer to the vector type
    using VectorPointerType = typename Kratos::shared_ptr<TVectorType>;

    /// The DoF updater type
    using DofUpdaterType = DofUpdater<UblasSpace<TDataType, TMatrixType, TVectorType>>;

    /// The pointer to the DoF updater type
    using DofUpdaterPointerType = typename DofUpdaterType::UniquePointer;

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

    static TDataType TwoNorm(const Matrix& rA) // Frobenious norm
    {
        TDataType aux_sum = TDataType();
        #pragma omp parallel for reduction(+:aux_sum)
        for (int i=0; i<static_cast<int>(rA.size1()); ++i) {
            for (int j=0; j<static_cast<int>(rA.size2()); ++j) {
                aux_sum += std::pow(rA(i,j),2);
            }
        }
        return std::sqrt(aux_sum);
    }

    static TDataType TwoNorm(const compressed_matrix<TDataType> & rA) // Frobenious norm
    {
        TDataType aux_sum = TDataType();

        const auto& r_values = rA.value_data();

        #pragma omp parallel for reduction(+:aux_sum)
        for (int i=0; i<static_cast<int>(r_values.size()); ++i) {
            aux_sum += std::pow(r_values[i] , 2);
        }
        return std::sqrt(aux_sum);
    }

    /**
     * This method computes the Jacobi norm
     * @param rA The matrix to compute the Jacobi norm
     * @return aux_sum: The Jacobi norm
     */
    static TDataType JacobiNorm(const Matrix& rA)
    {
        TDataType aux_sum = TDataType();
        #pragma omp parallel for reduction(+:aux_sum)
        for (int i=0; i<static_cast<int>(rA.size1()); ++i) {
            for (int j=0; j<static_cast<int>(rA.size2()); ++j) {
                if (i != j) {
                    aux_sum += std::abs(rA(i,j));
                }
            }
        }
        return aux_sum;
    }

    static TDataType JacobiNorm(const compressed_matrix<TDataType>& rA)
    {
        TDataType aux_sum = TDataType();

        typedef typename compressed_matrix<TDataType>::const_iterator1 t_it_1;
        typedef typename compressed_matrix<TDataType>::const_iterator2 t_it_2;

        for (t_it_1 it_1 = rA.begin1(); it_1 != rA.end1(); ++it_1) {
            for (t_it_2 it_2 = it_1.begin(); it_2 != it_1.end(); ++it_2) {
                if (it_2.index1() != it_2.index2()) {
                    aux_sum += std::abs(*it_2);
                }
            }
        }
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

    template< class TOtherMatrixType >
    static void TransposeMult(TOtherMatrixType& rA, VectorType& rX, VectorType& rY)
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

    /**
     * @brief This method checks and corrects the zero diagonal values
     * @details This method returns the scale norm considering for scaling the diagonal
     * @param rProcessInfo The problem process info
     * @param rA The LHS matrix
     * @param rb The RHS vector
     * @param ScalingDiagonal The type of caling diagonal considered
     * @return The scale norm
     */
    static double CheckAndCorrectZeroDiagonalValues(
        const ProcessInfo& rProcessInfo,
        MatrixType& rA,
        VectorType& rb,
        const SCALING_DIAGONAL ScalingDiagonal = SCALING_DIAGONAL::NO_SCALING
        )
    {
        const std::size_t system_size = rA.size1();

        const auto& Avalues = rA.value_data();
        const auto& Arow_indices = rA.index1_data();

        // Define  zero value tolerance
        const double zero_tolerance = std::numeric_limits<double>::epsilon();

        // The diagonal considered
        const double scale_factor = GetScaleNorm(rProcessInfo, rA, ScalingDiagonal);

        // Detect if there is a line of all zeros and set the diagonal to a 1 if this happens
        IndexPartition(system_size).for_each([&](std::size_t Index){
            bool empty = true;

            const std::size_t col_begin = Arow_indices[Index];
            const std::size_t col_end = Arow_indices[Index + 1];

            for (std::size_t j = col_begin; j < col_end; ++j) {
                if(std::abs(Avalues[j]) > zero_tolerance) {
                    empty = false;
                    break;
                }
            }

            if(empty) {
                rA(Index, Index) = scale_factor;
                rb[Index] = 0.0;
            }
        });

        return scale_factor;
    }

    /**
     * @brief This method returns the scale norm considering for scaling the diagonal
     * @param rProcessInfo The problem process info
     * @param rA The LHS matrix
     * @param ScalingDiagonal The type of caling diagonal considered
     * @return The scale norm
     */
    static double GetScaleNorm(
        const ProcessInfo& rProcessInfo,
        const MatrixType& rA,
        const SCALING_DIAGONAL ScalingDiagonal = SCALING_DIAGONAL::NO_SCALING
        )
    {
        switch (ScalingDiagonal) {
            case SCALING_DIAGONAL::NO_SCALING:
                return 1.0;
            case SCALING_DIAGONAL::CONSIDER_PRESCRIBED_DIAGONAL: {
                KRATOS_ERROR_IF_NOT(rProcessInfo.Has(BUILD_SCALE_FACTOR)) << "Scale factor not defined at process info" << std::endl;
                return rProcessInfo.GetValue(BUILD_SCALE_FACTOR);
            }
            case SCALING_DIAGONAL::CONSIDER_NORM_DIAGONAL:
                return GetDiagonalNorm(rA)/static_cast<double>(rA.size1());
            case SCALING_DIAGONAL::CONSIDER_MAX_DIAGONAL:
                return GetMaxDiagonal(rA);
            default:
                return GetMaxDiagonal(rA);
        }
    }

    /**
     * @brief This method returns the diagonal norm considering for scaling the diagonal
     * @param rA The LHS matrix
     * @return The diagonal norm
     */
    static double GetDiagonalNorm(const MatrixType& rA)
    {
        const auto& Avalues = rA.value_data();
        const auto& Arow_indices = rA.index1_data();
        const auto& Acol_indices = rA.index2_data();

        const double diagonal_norm = IndexPartition<std::size_t>(Size1(rA)).for_each<SumReduction<double>>([&](std::size_t Index){
            const std::size_t col_begin = Arow_indices[Index];
            const std::size_t col_end = Arow_indices[Index+1];
            for (std::size_t j = col_begin; j < col_end; ++j) {
                if (Acol_indices[j] == Index ) {
                    return std::pow(Avalues[j], 2);
                }
            }
            return 0.0;
        });

        return std::sqrt(diagonal_norm);
    }

    /**
     * @brief This method returns the diagonal max value
     * @param rA The LHS matrix
     * @return The diagonal  max value
     */
    static double GetAveragevalueDiagonal(const MatrixType& rA)
    {
        return 0.5 * (GetMaxDiagonal(rA) + GetMinDiagonal(rA));
    }

    /**
     * @brief This method returns the diagonal max value
     * @param rA The LHS matrix
     * @return The diagonal  max value
     */
    static double GetMaxDiagonal(const MatrixType& rA)
    {
        const auto& Avalues = rA.value_data();
        const auto& Arow_indices = rA.index1_data();
        const auto& Acol_indices = rA.index2_data();

        return IndexPartition<std::size_t>(Size1(rA)).for_each<MaxReduction<double>>([&](std::size_t Index){
            const std::size_t col_begin = Arow_indices[Index];
            const std::size_t col_end = Arow_indices[Index+1];
            for (std::size_t j = col_begin; j < col_end; ++j) {
                if (Acol_indices[j] == Index ) {
                    return std::abs(Avalues[j]);
                }
            }
            return std::numeric_limits<double>::lowest();
        });
    }

    /**
     * @brief This method returns the diagonal min value
     * @param rA The LHS matrix
     * @return The diagonal min value
     */
    static double GetMinDiagonal(const MatrixType& rA)
    {
        const auto& Avalues = rA.value_data();
        const auto& Arow_indices = rA.index1_data();
        const auto& Acol_indices = rA.index2_data();

        return IndexPartition<std::size_t>(Size1(rA)).for_each<MinReduction<double>>([&](std::size_t Index){
            const std::size_t col_begin = Arow_indices[Index];
            const std::size_t col_end = Arow_indices[Index+1];
            for (std::size_t j = col_begin; j < col_end; ++j) {
                if (Acol_indices[j] == Index ) {
                    return std::abs(Avalues[j]);
                }
            }
            return std::numeric_limits<double>::max();
        });
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

    inline static constexpr bool IsDistributed()
    {
        return false;
    }

    /**
     * @brief Returns a list of the fastest direct solvers.
     * @details This function returns a vector of strings representing the names of the fastest direct solvers. The order of the solvers in the list may need to be updated and reordered depending on the size of the equation system.
     * @return A vector of strings containing the names of the fastest direct solvers.
     */
    inline static std::vector<std::string> FastestDirectSolverList()
    {
        std::vector<std::string> faster_direct_solvers({
            "pardiso_lu",              // LinearSolversApplication (if compiled with Intel-support)
            "pardiso_ldlt",            // LinearSolversApplication (if compiled with Intel-support)
            "sparse_lu",               // LinearSolversApplication
            "skyline_lu_factorization" // In Core, always available, but slow
        });
        return faster_direct_solvers;
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
    static bool WriteMatrixMarketVector(const char* pFileName, /*const*/ VectorType& rV)
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
