//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>

// External includes
#include "amgcl/value_type/interface.hpp"

// Project includes
#include "includes/define.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/atomic_utilities.h"

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

/**
 * @class SparseMatrixMultiplicationUtility
 * @ingroup KratosCore
 * @brief An utility to multiply sparse matrix in Ublas
 * @details Taken and adapted for ublas from external_libraries/amgcl/detail/spgemm.hpp by Denis Demidov <dennis.demidov@gmail.com>
 * @todo Remove as soon as we do not depend of Ublas anymore...
 * @author Vicente Mataix Ferrandiz
 */
class SparseMatrixMultiplicationUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TreeContactSearch
    KRATOS_CLASS_POINTER_DEFINITION( SparseMatrixMultiplicationUtility );

    /// The size type
    using SizeType = std::size_t;

    /// The index type
    using IndexType = std::size_t;

    /// The signed index type
    using SignedIndexType = std::ptrdiff_t;

    /// A vector of indexes
    using IndexVectorType = DenseVector<IndexType>;

    /// A vector of indexes (signed)
    using SignedIndexVectorType = DenseVector<SignedIndexType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    SparseMatrixMultiplicationUtility() = delete;

    /// Desctructor
    virtual ~SparseMatrixMultiplicationUtility()= default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Matrix-matrix product C = A·B
     * @detail This method uses a template for each matrix
     * @param rA The first matrix
     * @param rB The second matrix
     * @param rC The resulting matrix
     */
    template <class AMatrix, class BMatrix, class CMatrix>
    static void MatrixMultiplication(
        const AMatrix& rA,
        const BMatrix& rB,
        CMatrix& rC
        )
    {
        // We check the number of threads
        const unsigned int number_of_threads = ParallelUtilities::GetNumThreads();

        if (number_of_threads > 16) {
            MatrixMultiplicationRMerge(rA, rB, rC);
        } else {
            MatrixMultiplicationSaad(rA, rB, rC);
        }
    }

    /**
     * @brief The first is an OpenMP-enabled modification of classic algorithm from Saad
     * @details It is used whenever number of OpenMP cores is 4 or less. Saad, Yousef. Iterative methods for sparse linear systems. Siam, 2003.
     * @param A The first matrix to multiply
     * @param B The second matrix to multiply
     * @param C The resulting matrix
     */
    template <class AMatrix, class BMatrix, class CMatrix>
    static void MatrixMultiplicationSaad(
        const AMatrix& A,
        const BMatrix& B,
        CMatrix& C
        )
    {
        using ValueType = typename CMatrix::value_type;

        // Auxiliary sizes
        const SizeType nrows = A.size1();
        const SizeType ncols = B.size2();

        // Exiting just in case of empty matrix
        if ((nrows == 0) || (ncols == 0))
            return void();

        // Get access to A, B and C data
        const IndexType* index1_a = A.index1_data().begin();
        const IndexType* index2_a = A.index2_data().begin();
        const typename AMatrix::value_type* values_a = A.value_data().begin();
        const IndexType* index1_b = B.index1_data().begin();
        const IndexType* index2_b = B.index2_data().begin();
        const typename BMatrix::value_type* values_b = B.value_data().begin();
        IndexType* c_ptr = new IndexType[nrows + 1];

        c_ptr[0] = 0;

        // Definition of marker TLS
        struct TLS {
            TLS(const SizeType NumberOfColumns)
                : marker(NumberOfColumns)
            {
                std::fill(marker.begin(), marker.end(), -1);
            }
            SignedIndexVectorType marker;
        };

        IndexPartition<SignedIndexType>(static_cast<SignedIndexType>(nrows)).for_each(TLS(ncols), [&](SignedIndexType ia, TLS& rTLS) {
            const IndexType row_begin_a = index1_a[ia];
            const IndexType row_end_a   = index1_a[ia+1];

            IndexType C_cols = 0;
            for(IndexType ja = row_begin_a; ja < row_end_a; ++ja) {
                const IndexType ca = index2_a[ja];
                const IndexType row_begin_b = index1_b[ca];
                const IndexType row_end_b   = index1_b[ca+1];

                for(IndexType jb = row_begin_b; jb < row_end_b; ++jb) {
                    const IndexType cb = index2_b[jb];
                    if (rTLS.marker[cb] != ia) {
                        rTLS.marker[cb]  = ia;
                        ++C_cols;
                    }
                }
            }
            c_ptr[ia + 1] = C_cols;
        });

        // We initialize the sparse matrix
        std::partial_sum(c_ptr, c_ptr + nrows + 1, c_ptr);
        const SizeType nonzero_values = c_ptr[nrows];
        IndexType* aux_index2_c = new IndexType[nonzero_values];
        ValueType* aux_val_c = new ValueType[nonzero_values];

        IndexPartition<IndexType>(nrows).for_each(TLS(ncols), [&](IndexType ia, TLS& rTLS) {
            const IndexType row_begin_a = index1_a[ia];
            const IndexType row_end_a   = index1_a[ia+1];

            const IndexType row_beg = c_ptr[ia];
            IndexType row_end = row_beg;

            for(IndexType ja = row_begin_a; ja < row_end_a; ++ja) {
                const IndexType ca = index2_a[ja];
                const ValueType va = values_a[ja];

                const IndexType row_begin_b = index1_b[ca];
                const IndexType row_end_b   = index1_b[ca+1];

                for(IndexType jb = row_begin_b; jb < row_end_b; ++jb) {
                    const IndexType cb = index2_b[jb];
                    const ValueType vb = values_b[jb];

                    if (rTLS.marker[cb] < static_cast<SignedIndexType>(row_beg)) {
                        rTLS.marker[cb] = row_end;
                        aux_index2_c[row_end] = cb;
                        aux_val_c[row_end] = va * vb;
                        ++row_end;
                    } else {
                        aux_val_c[rTLS.marker[cb]] += va * vb;
                    }
                }
            }
        });

        // We reorder the rows
        SortRows(c_ptr, nrows, ncols, aux_index2_c, aux_val_c);

        // We fill the matrix
        CreateSolutionMatrix(C, nrows, ncols, c_ptr, aux_index2_c, aux_val_c);

        // Release memory
        delete[] c_ptr;
        delete[] aux_index2_c;
        delete[] aux_val_c;
    }

    /**
     * @brief Row-merge algorithm from Rupp et al.
     * @details The algorithm  requires less memory and shows much better scalability than classic one. It is used when number of OpenMP cores is more than 4.
     * @param A The first matrix to multiply
     * @param B The second matrix to multiply
     * @param C The resulting matrix
     */
    template <class AMatrix, class BMatrix, class CMatrix>
    static void MatrixMultiplicationRMerge(
        const AMatrix &A,
        const BMatrix &B,
        CMatrix &C
        )
    {
        using AIndex = typename AMatrix::index_array_type::value_type;
        using AValue = typename AMatrix::value_type;
        using BIndex = typename BMatrix::index_array_type::value_type;
        using BValue = typename BMatrix::value_type;
        using CIndex = typename CMatrix::index_array_type::value_type;
        using CValue = typename CMatrix::value_type;

        // Auxiliary sizes
        const AIndex nrows = A.size1();
        const BIndex ncols = B.size2();

        // Exiting just in case of empty matrix
        if ((nrows == 0) || (ncols == 0))
            return void();

        // Get access to A and B data
        const AIndex* index1_a = A.index1_data().begin();
        const AIndex* index2_a = A.index2_data().begin();
        const AValue* values_a = A.value_data().begin();
        const BIndex* index1_b = B.index1_data().begin();
        const BIndex* index2_b = B.index2_data().begin();
        const BValue* values_b = B.value_data().begin();

        // Definition of TLS
        struct TLS_max {
            IndexType my_max = 0;
        };

        const IndexType max_row_width = IndexPartition<IndexType>(nrows).for_each<MaxReduction<IndexType>>(TLS_max(), [&](IndexType i, TLS_max& rTLS) {
            const IndexType row_beg = index1_a[i];
            const IndexType row_end = index1_a[i+1];

            IndexType row_width = 0;
            for(IndexType j = row_beg; j < row_end; ++j) {
                const IndexType a_col = index2_a[j];
                row_width += index1_b[a_col + 1] - index1_b[a_col];
            }
            rTLS.my_max = std::max(rTLS.my_max, row_width);
            return rTLS.my_max;
        });

        // We check the number of threads
        const unsigned int number_of_threads = ParallelUtilities::GetNumThreads();

        std::vector<std::vector<CIndex>> tmp_col(number_of_threads);
        std::vector<std::vector<CValue>> tmp_val(number_of_threads);

        for(unsigned int i = 0; i < number_of_threads; ++i) {
            tmp_col[i].resize(3 * max_row_width);
            tmp_val[i].resize(2 * max_row_width);
        }

        // We create the c_ptr auxiliary variable
        CIndex* c_ptr = new CIndex[nrows + 1];
        c_ptr[0] = 0;

        IndexPartition<IndexType>(nrows).for_each([&](IndexType i) {
            const IndexType row_beg = index1_a[i];
            const IndexType row_end = index1_a[i+1];

            int tid = 0;
            #ifdef _OPENMP
                tid = omp_get_thread_num();
            #endif
            CIndex* t_col = &tmp_col[tid][0];

            c_ptr[i+1] = ProdRowWidth(index2_a + row_beg, index2_a + row_end,
                                      index1_b, index2_b,
                                      t_col, t_col + max_row_width, t_col + 2 * max_row_width);
        });

        // We initialize the sparse matrix
        std::partial_sum(c_ptr, c_ptr + nrows + 1, c_ptr);
        const SizeType nonzero_values = c_ptr[nrows];
        CIndex* aux_index2_c = new CIndex[nonzero_values];
        CValue* aux_val_c = new CValue[nonzero_values];

        IndexPartition<IndexType>(nrows).for_each([&](IndexType i) {
            const IndexType row_beg = index1_a[i];
            const IndexType row_end = index1_a[i+1];

            int tid = 0;
            #ifdef _OPENMP
                tid = omp_get_thread_num();
            #endif

            CIndex* t_col = tmp_col[tid].data();
            CValue* t_val = tmp_val[tid].data();

            ProdRow(index2_a + row_beg, index2_a + row_end, values_a + row_beg,
                    index1_b, index2_b, values_b,
                    aux_index2_c + c_ptr[i], aux_val_c + c_ptr[i],
                    t_col, t_val,
                    t_col + max_row_width, t_val + max_row_width);
        });

        // We fill the matrix
        CreateSolutionMatrix(C, nrows, ncols, c_ptr, aux_index2_c, aux_val_c);

        // Release memory
        delete[] c_ptr;
        delete[] aux_index2_c;
        delete[] aux_val_c;
    }

    /**
     * @brief The first is a method in order to sum to sparse matrices in a efficient way
     * @param A The resulting matrix
     * @param B The second matrix to sum
     */
    template <class AMatrix, class BMatrix>
    static void MatrixAdd(
        AMatrix& A,
        const BMatrix& B,
        const double Factor = 1.0
        )
    {
        using ValueType = typename AMatrix::value_type;

        // Auxiliary sizes
        const SizeType nrows = A.size1();
        const SizeType ncols = A.size2();

        /* Some checks */
        // Exiting just in case of empty matrix
        if ((nrows == 0) || (ncols == 0))
            return void();

        KRATOS_ERROR_IF_NOT(nrows == B.size1()) << "The second matrix has a wrong number of rows" << std::endl;
        KRATOS_ERROR_IF_NOT(ncols == B.size2()) << "The second matrix has a wrong number of columns" << std::endl;

        // Get access to A and B data
        const IndexType* index1_a = A.index1_data().begin();
        const IndexType* index2_a = A.index2_data().begin();
        const double* values_a = A.value_data().begin();
        const IndexType* index1_b = B.index1_data().begin();
        const IndexType* index2_b = B.index2_data().begin();
        const double* values_b = B.value_data().begin();

        IndexType* new_a_ptr = new IndexType[nrows + 1];
        new_a_ptr[0] = 0;

        IndexPartition<IndexType>(nrows).for_each([&](IndexType ia) {
            SignedIndexVectorType marker(ncols);
            for (int i = 0; i < static_cast<int>(ncols); ++i) {
                marker[i] = -1;
            }

            // Initialize
            IndexType new_A_cols = 0;

            // Iterate over A
            const IndexType row_begin_a = index1_a[ia];
            const IndexType row_end_a   = index1_a[ia+1];
            for(IndexType ja = row_begin_a; ja < row_end_a; ++ja) {
                const IndexType ca = index2_a[ja];
                marker[ca] = 1;
                ++new_A_cols;
            }

            // Iterate over B
            const IndexType row_begin_b = index1_b[ia];
            const IndexType row_end_b   = index1_b[ia+1];
            for(IndexType jb = row_begin_b; jb < row_end_b; ++jb) {
                const IndexType cb = index2_b[jb];
                if (marker[cb] < 0) {
                    marker[cb] = 1;
                    ++new_A_cols;
                }
            }
            new_a_ptr[ia + 1] = new_A_cols;
        });

        // We initialize the sparse matrix
        std::partial_sum(new_a_ptr, new_a_ptr + nrows + 1, new_a_ptr);
        const SizeType nonzero_values = new_a_ptr[nrows];
        IndexType* aux_index2_new_a = new IndexType[nonzero_values];
        ValueType* aux_val_new_a = new ValueType[nonzero_values];

        IndexPartition<IndexType>(nrows).for_each([&](IndexType ia) {
            SignedIndexVectorType marker(ncols);
            std::fill(marker.begin(), marker.end(), -1);

            // Initialize
            const IndexType row_beg = new_a_ptr[ia];
            IndexType row_end = row_beg;

            // Iterate over A
            const IndexType row_begin_a = index1_a[ia];
            const IndexType row_end_a   = index1_a[ia+1];
            for(IndexType ja = row_begin_a; ja < row_end_a; ++ja) {
                const IndexType ca = index2_a[ja];
                const ValueType va = values_a[ja];

                marker[ca] = row_end;
                aux_index2_new_a[row_end] = ca;
                aux_val_new_a[row_end] = va;
                ++row_end;
            }

            // Iterate over B
            const IndexType row_begin_b = index1_b[ia];
            const IndexType row_end_b   = index1_b[ia+1];
            for(IndexType jb = row_begin_b; jb < row_end_b; ++jb) {
                const IndexType cb = index2_b[jb];
                const ValueType vb = values_b[jb];

                if (marker[cb] < 0) {
                    marker[cb] = row_end;
                    aux_index2_new_a[row_end] = cb;
                    aux_val_new_a[row_end] = Factor * vb;
                    ++row_end;
                } else {
                    aux_val_new_a[marker[cb]] += Factor * vb;
                }
            }
        });

        // We reorder the rows
        SortRows(new_a_ptr, nrows, ncols, aux_index2_new_a, aux_val_new_a);

        // We fill the matrix
        CreateSolutionMatrix(A, nrows, ncols, new_a_ptr, aux_index2_new_a, aux_val_new_a);

        // Release memory
        delete[] new_a_ptr;
        delete[] aux_index2_new_a;
        delete[] aux_val_new_a;
    }

    /**
     * @brief This method computes of the transpose matrix of a given matrix
     * @param rA The resulting matrix
     * @param rB The second matrix to transpose
     */
    template <class AMatrix, class BMatrix>
    static void TransposeMatrix(
        AMatrix& rA,
        const BMatrix& rB,
        const double Factor = 1.0
        )
    {
        using ValueType = typename AMatrix::value_type;

        // Get access to B data
        const IndexType* index1 = rB.index1_data().begin();
        const IndexType* index2 = rB.index2_data().begin();
        const ValueType* data = rB.value_data().begin();
        const SizeType transpose_nonzero_values = rB.value_data().end() - rB.value_data().begin();

        const SizeType size_system_1 = rB.size1();
        const SizeType size_system_2 = rB.size2();

        if (rA.size1() != size_system_2 || rA.size2() != size_system_1 ) {
            rA.resize(size_system_2, size_system_1, false);
        }

        IndexVectorType new_a_ptr(size_system_2 + 1);
        IndexPartition<IndexType>(size_system_2 + 1).for_each([&](IndexType i) {
            new_a_ptr[i] = 0;
        });
        IndexVectorType aux_index2_new_a(transpose_nonzero_values);
        DenseVector<ValueType> aux_val_new_a(transpose_nonzero_values);

        // Auxiliary index 1
        const IndexType aux_1_index = 1;

        IndexPartition<IndexType>(size_system_1).for_each([&](IndexType i) {
            IndexType row_begin = index1[i];
            IndexType row_end   = index1[i+1];

            for (IndexType j=row_begin; j<row_end; j++) {
                AtomicAdd(new_a_ptr[index2[j] + 1], aux_1_index);
            }
        });

        // We initialize the blocks sparse matrix
        std::partial_sum(new_a_ptr.begin(), new_a_ptr.end(), &new_a_ptr[0]);

        IndexVectorType aux_indexes(size_system_2);
        IndexPartition<IndexType>(size_system_2).for_each([&](IndexType i) {
            aux_indexes[i] = 0;
        });

//         #pragma omp parallel for
        for (int i=0; i<static_cast<int>(size_system_1); ++i) {
            IndexType row_begin = index1[i];
            IndexType row_end   = index1[i+1];

            for (IndexType j=row_begin; j<row_end; j++) {
                const IndexType current_row = index2[j];
                const IndexType initial_position = new_a_ptr[current_row];

                const IndexType current_index = initial_position + aux_indexes[current_row];
                aux_index2_new_a[current_index] = i;
                aux_val_new_a[current_index] = Factor * data[j];

                // AtomicAdd(aux_indexes[current_row], aux_1_index);
                aux_indexes[current_row] += 1;
            }

        }

        // We reorder the rows
        SortRows(&new_a_ptr[0], size_system_2, size_system_1, &aux_index2_new_a[0], &aux_val_new_a[0]);

        // We fill the matrix
        CreateSolutionMatrix(rA, size_system_2, size_system_1, &new_a_ptr[0], &aux_index2_new_a[0], &aux_val_new_a[0]);
    }

    /**
     * @brief This method is designed to create the final solution sparse matrix from the auxiliary values
     * @param C The matrix solution
     * @param NRows The number of rows of the matrix
     * @param NCols The number of columns of the matrix
     * @param CPtr The indexes taht indicate the number of nonzero values in each column
     * @param AuxIndex2C The indexes of the nonzero columns
     * @param AuxValC The C array containing the values of the sparse matrix
     */
    template <class CMatrix, typename TSize, typename Ptr, typename IndexType, typename ValueType>
    static inline void CreateSolutionMatrix(
        CMatrix& C,
        const TSize NRows,
        const TSize NCols,
        const Ptr* CPtr,
        const IndexType* AuxIndex2C,
        const ValueType* AuxValC
        )
    {
        // Exiting just in case of empty matrix
        if ((NRows == 0) || (NCols == 0))
            return void();

        // Auxiliary values
        const TSize nonzero_values = CPtr[NRows];

        C = CMatrix(NRows, NCols, nonzero_values);
        IndexType* index1_c = C.index1_data().begin();
        IndexType* index2_c = C.index2_data().begin();
        ValueType* values_c = C.value_data().begin();

        index1_c[0] = 0;
        for (TSize i = 0; i < NRows; i++) {
            index1_c[i+1] = index1_c[i] + (CPtr[i+1] - CPtr[i]);
        }

        IndexPartition<IndexType>(nonzero_values).for_each([&](IndexType i) {
            KRATOS_DEBUG_ERROR_IF(AuxIndex2C[i] > static_cast<IndexType>(NCols)) << "Index " << AuxIndex2C[i] <<" is greater than the number of columns " << NCols << std::endl;
            index2_c[i] = AuxIndex2C[i];
            values_c[i] = AuxValC[i];
        });

        C.set_filled(NRows+1, nonzero_values);
    }

    /**
     * @brief This method is designed to reorder the rows by columns
     * @param NRows The number of rows of the matrix
     * @param NCols The number of columns of the matrix
     * @param CPtr The indexes that indicate the number of nonzero values in each column
     * @param Columns The columns of the problem
     * @param Values The values (to be ordered with the rows)
     */
    template <typename TSize, typename Col, typename TIndexType, typename ValueType>
    static inline void SortRows(
        const TIndexType* CPtr,
        const TSize NRows,
        const TSize NCols,
        Col* Columns,
        ValueType* Values
        )
    {
        IndexPartition<TSize>(NRows).for_each([&](TSize i_row) {
            const TIndexType row_begin = CPtr[i_row];
            const TIndexType row_end = CPtr[i_row + 1];
            TSize length = row_end - row_begin;

            // Use vectors to facilitate sorting
            std::vector<std::pair<Col, ValueType>> row_data(length);
            for (TSize j = 0; j < length; ++j) {
                row_data[j] = {Columns[j + row_begin], Values[j + row_begin]};
            }

            // Sort using a more efficient parallel algorithm
            std::sort(row_data.begin(), row_data.end(), [](const auto& a, const auto& b) {
                return a.first < b.first;
            });

            // Write back sorted data
            for (TSize j = 0; j < length; ++j) {
                Columns[j + row_begin] = row_data[j].first;
                Values[j + row_begin] = row_data[j].second;
            }
        });
    }

    /**
     * @brief This method assembles several sparse matrices into one large sparse matrix
     * @param rMatricespBlocks The pointers to the matrices we are interested in assemble
     * @param ContributionCoefficients The matrix containing the coefficients to be considered (copy, so we don't need to provide it)
     * @param TransposeBlocks The matrix containing the flags telling us to transpose the blocks (copy, so we don't need to provide it)
     */
    static inline void AssembleSparseMatrixByBlocks(
        CompressedMatrix& rMatrix,
        const DenseMatrix<CompressedMatrix*>& rMatricespBlocks,
        DenseMatrix<double> ContributionCoefficients = DenseMatrix<double>(0,0),
        DenseMatrix<bool> TransposeBlocks = DenseMatrix<bool>(0,0)
        )
    {
        const SizeType number_of_rows_blocks = rMatricespBlocks.size1();
        const SizeType number_of_columns_blocks = rMatricespBlocks.size2();

        // Fill the matrices if they are empty
        if (ContributionCoefficients.size1() == 0 && ContributionCoefficients.size2() == 0) {
            ContributionCoefficients.resize(number_of_rows_blocks, number_of_columns_blocks);
            for (IndexType i = 0; i < number_of_rows_blocks; ++i) {
                for (IndexType j = 0; j < number_of_columns_blocks; ++j) {
                    ContributionCoefficients(i, j) = 1.0;
                }
            }
        } else {
            KRATOS_ERROR_IF(ContributionCoefficients.size1() != number_of_rows_blocks || ContributionCoefficients.size2() != number_of_columns_blocks) << "The ContributionCoefficients dimensions" << ContributionCoefficients.size1() << " and " << ContributionCoefficients.size2() << "do not coincide with the dimensions of rMatricespBlocks" << number_of_rows_blocks << "and " << number_of_columns_blocks << std::endl;
        }
        if (TransposeBlocks.size1() == 0 && TransposeBlocks.size2() == 0) {
            TransposeBlocks.resize(number_of_rows_blocks, number_of_columns_blocks);
            for (IndexType i = 0; i < number_of_rows_blocks; ++i) {
                for (IndexType j = 0; j < number_of_rows_blocks; ++j) {
                    TransposeBlocks(i, j) = false;
                }
            }
        } else {
            KRATOS_ERROR_IF(TransposeBlocks.size1() != number_of_rows_blocks || TransposeBlocks.size2() != number_of_columns_blocks) << "The TransposeBlocks dimensions" << TransposeBlocks.size1() << " and " << TransposeBlocks.size2() << "do not coincide with the dimensions of rMatricespBlocks" << number_of_rows_blocks << "and " << number_of_columns_blocks << std::endl;
        }

        // Compute total size and check consistency of the different blocks
        SizeType nrows = 0, ncols = 0;
        std::vector<SizeType> row_sizes(number_of_rows_blocks);
        std::vector<SizeType> column_sizes(number_of_columns_blocks);
        for (int i=0; i<static_cast<int>(number_of_rows_blocks); ++i) {
            if (TransposeBlocks(i, 0)) {
                row_sizes[i] = (*rMatricespBlocks(i, 0)).size2();
            } else {
                row_sizes[i] = (*rMatricespBlocks(i, 0)).size1();
            }
            nrows += row_sizes[i];
        }
        for (int j=0; j<static_cast<int>(number_of_columns_blocks); ++j) {
            if (TransposeBlocks(0, j)) {
                column_sizes[j] = (*rMatricespBlocks(0, j)).size1();
            } else {
                column_sizes[j] = (*rMatricespBlocks(0, j)).size2();
            }
            ncols += column_sizes[j];
        }

        // Check consistency of all blocks
        for (int i=0; i<static_cast<int>(number_of_rows_blocks); ++i) {
            for (int j=0; j<static_cast<int>(number_of_columns_blocks); ++j) {
                if (TransposeBlocks(i, j)) {
                    KRATOS_ERROR_IF((*rMatricespBlocks(i, j)).size2() != row_sizes[i] || (*rMatricespBlocks(i, j)).size1() != column_sizes[j]) << " Not consistent size in block " << i << ", " << j << ".\t" << (*rMatricespBlocks(i, j)).size2() << ", " << (*rMatricespBlocks(i, j)).size1() << " vs " <<  row_sizes[i] << ", " << row_sizes[j] << std::endl;
                } else {
                    KRATOS_ERROR_IF((*rMatricespBlocks(i, j)).size1() != row_sizes[i] || (*rMatricespBlocks(i, j)).size2() != column_sizes[j]) << " Not consistent size in block " << i << ", " << j << ".\t" << (*rMatricespBlocks(i, j)).size1() << ", " << (*rMatricespBlocks(i, j)).size2() << " vs " <<  row_sizes[i] << ", " << row_sizes[j] << std::endl;
                }
            }
        }
        // Exiting just in case of empty matrix
        if ((nrows == 0) || (ncols == 0))
            return void();

        // We will compute nonzero terms
        IndexType* matrix_ptr = new IndexType[nrows + 1];
        IndexPartition<IndexType>(nrows + 1).for_each([&](IndexType i) {
            matrix_ptr[i] = 0;
        });

    #ifdef KRATOS_DEBUG
        IndexType check_non_zero = 0;
        DenseMatrix<IndexType> check_non_zero_blocks(number_of_rows_blocks, number_of_columns_blocks);
        for (int i=0; i<static_cast<int>(number_of_rows_blocks); ++i) {
            for (int j=0; j<static_cast<int>(number_of_columns_blocks); ++j) {
                check_non_zero_blocks(i, j) = 0;
            }
        }
    #endif

        IndexPartition<IndexType>(number_of_rows_blocks).for_each([&](IndexType i) {
            for (int k=0; k<static_cast<int>(row_sizes[i]); ++k) {
                IndexType matrix_cols_aux = 0;
                for (int j=0; j<static_cast<int>(number_of_columns_blocks); ++j) {
                #ifdef KRATOS_DEBUG
                    IndexType partial_matrix_cols_aux = 0;
                #endif
                    // Skip if empty matrix
                    CompressedMatrix& r_matrix = *rMatricespBlocks(i, j);
                    if (r_matrix.nnz() > 0) {
                        if (TransposeBlocks(i, j)) {
                            // We compute the transposed matrix
                            const SizeType size_system_1 = r_matrix.size1();
                            const SizeType size_system_2 = r_matrix.size2();
                            CompressedMatrix transpose(size_system_2, size_system_1);
                            TransposeMatrix<CompressedMatrix, CompressedMatrix>(transpose, r_matrix);
                            ComputeNonZeroBlocks(transpose, k, matrix_cols_aux);
                        #ifdef KRATOS_DEBUG
                            ComputeNonZeroBlocks(transpose, k, partial_matrix_cols_aux);
                        #endif
                        } else {
                            ComputeNonZeroBlocks(r_matrix, k, matrix_cols_aux);
                        #ifdef KRATOS_DEBUG
                            ComputeNonZeroBlocks(r_matrix, k, partial_matrix_cols_aux);
                        #endif
                        }
                    }
                #ifdef KRATOS_DEBUG
                    check_non_zero_blocks(i, j) += partial_matrix_cols_aux;
                #endif
                }
                IndexType& r_matrix_ptr_value = matrix_ptr[std::accumulate(row_sizes.begin(), row_sizes.begin() + i, 0) + k + 1];
                AtomicAdd(r_matrix_ptr_value, matrix_cols_aux);

            #ifdef KRATOS_DEBUG
                AtomicAdd(check_non_zero, matrix_cols_aux);
            #endif
            }
        });

        // Auxiliary values
        std::partial_sum(matrix_ptr, matrix_ptr + nrows + 1, matrix_ptr);
        const SizeType nonzero_values = matrix_ptr[nrows];

    #ifdef KRATOS_DEBUG
        SizeType total_nnz = 0;
        for (int i=0; i<static_cast<int>(number_of_rows_blocks); ++i) {
            for (int j=0; j<static_cast<int>(number_of_columns_blocks); ++j) {
                const SizeType block_nnz = rMatricespBlocks(i, j)->nnz();
                KRATOS_ERROR_IF_NOT(check_non_zero_blocks(i, j) == block_nnz) << "Inconsistent number of non-zero values. Check 0: " << block_nnz << " vs " << check_non_zero_blocks(i, j) << ". Block: " << i << ", " << j << std::endl;
                total_nnz += block_nnz;
            }
        }
        KRATOS_ERROR_IF_NOT(check_non_zero == total_nnz) << "Inconsistent number of non-zero values. Check 1: " << total_nnz << " vs " << check_non_zero << std::endl;
        KRATOS_ERROR_IF_NOT(nonzero_values == total_nnz) << "Inconsistent number of non-zero values. Check 2: " << total_nnz << " vs " << nonzero_values << std::endl;
    #endif

        // Initialize matrix with the corresponding non-zero values
        rMatrix = CompressedMatrix(nrows, ncols, nonzero_values);

        // Fill the new matrix
        double* Matrix_values = rMatrix.value_data().begin();
        IndexType* Matrix_index1 = rMatrix.index1_data().begin();
        IndexType* Matrix_index2 = rMatrix.index2_data().begin();

        Matrix_index1[0] = 0;
        for (IndexType i = 0; i < nrows; ++i) {
            Matrix_index1[i+1] = Matrix_index1[i] + (matrix_ptr[i + 1] - matrix_ptr[i]);
        }

        IndexPartition<IndexType>(number_of_rows_blocks).for_each([&](IndexType i) {
            for (int k=0; k<static_cast<int>(row_sizes[i]); ++k) {
                const IndexType row_beg = matrix_ptr[std::accumulate(row_sizes.begin(), row_sizes.begin() + i, 0) + k];
                IndexType row_end = row_beg;
                for (int j=0; j<static_cast<int>(number_of_columns_blocks); ++j) {
                    const SizeType initial_index_column = std::accumulate(column_sizes.begin(), column_sizes.begin() + j, 0);

                    // Skip if empty matrix
                    CompressedMatrix& r_matrix = *rMatricespBlocks(i, j);
                    if (r_matrix.nnz() > 0) {
                        if (TransposeBlocks(i, j)) {
                            // We compute the transposed matrix
                            const SizeType size_system_1 = r_matrix.size1();
                            const SizeType size_system_2 = r_matrix.size2();
                            CompressedMatrix transpose(size_system_2, size_system_1);
                            TransposeMatrix<CompressedMatrix, CompressedMatrix>(transpose, r_matrix);
                            ComputeAuxiliarValuesBlocks(transpose, Matrix_index2, Matrix_values, k, row_end, initial_index_column, ContributionCoefficients(i, j));
                        } else {
                            ComputeAuxiliarValuesBlocks(r_matrix, Matrix_index2, Matrix_values, k, row_end, initial_index_column, ContributionCoefficients(i, j));
                        }
                    }
                }
            }
        });

        // Close the matrix
        rMatrix.set_filled(nrows+1, nonzero_values);

        // Release memory
        delete[] matrix_ptr;
    }

    /**
     * @brief This is a method to check the block containing nonzero values
     * @param rMatrix The auxiliary block
     * @param CurrentRow The current row computed
     * @param rNonZeroColsAux2 The nonzero rows array
     */
    static inline void ComputeNonZeroBlocks(
        const CompressedMatrix& rMatrix,
        const int CurrentRow,
        IndexType& rNonZeroColsAux2
        )
    {
        // Get access to aux_K data
        const IndexType* aux_matrix_index1 = rMatrix.index1_data().begin();

        const IndexType row_begin = aux_matrix_index1[CurrentRow];
        const IndexType row_end   = aux_matrix_index1[CurrentRow + 1];

        for (IndexType j=row_begin; j<row_end; j++) {
            ++rNonZeroColsAux2;
        }
    }

    /**
     * @brief This is a method to compute the contribution of the auxiliary blocks
     * @param AuxK The auxiliary block
     * @param AuxIndex2 The indexes of the non zero columns
     * @param AuxVals The values of the final matrix
     * @param CurrentRow The current row computed
     * @param RowEnd The last column computed
     * @param InitialIndexColumn The initial column index of the auxiliary block in the final matrix
     */
    static inline void ComputeAuxiliarValuesBlocks(
        const CompressedMatrix& rMatrix,
        IndexType* AuxIndex2,
        double* AuxVals,
        const int CurrentRow,
        IndexType& RowEnd,
        const SizeType InitialIndexColumn,
        const double ContributionCoefficient = 1.0
        )
    {
        // Get access to aux_K data
        const double* aux_values = rMatrix.value_data().begin();
        const IndexType* aux_Matrix_index1 = rMatrix.index1_data().begin();
        const IndexType* aux_Matrix_index2 = rMatrix.index2_data().begin();

        const IndexType aux_Matrix_row_begin = aux_Matrix_index1[CurrentRow];
        const IndexType aux_Matrix_row_end   = aux_Matrix_index1[CurrentRow + 1];

        for (IndexType j=aux_Matrix_row_begin; j<aux_Matrix_row_end; j++) {
            const IndexType col_index = InitialIndexColumn + aux_Matrix_index2[j];
            AuxIndex2[RowEnd] = col_index;
            AuxVals[RowEnd] = ContributionCoefficient * aux_values[j];
            ++RowEnd;
        }
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        return "SparseMatrixMultiplicationUtility";
    }

    /// Print information about this object.
    void PrintInfo (std::ostream& rOStream) const
    {
        rOStream << "SparseMatrixMultiplicationUtility";
    }

    /// Print object's data.
    void PrintData (std::ostream& rOStream) const
    {
    }

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

    /**
     * @brief This method is oriented to merge rows
     * @param Column1 The index of the first matrix column
     * @param Column1End The last index of the first matrix column
     * @param Column2 The index of the second matrix column
     * @param Column2End The last index of the second matrix column
     * @param Column3 The index of the third matrix column
     * @return The resulting row
     */
    template <bool TNeedOut, class TAIndex, class TBIndex, class TCIndex>
    static TCIndex* MergeRows(
        const TAIndex* Column1,
        const TAIndex* Column1End,
        const TBIndex* Column2,
        const TBIndex* Column2End,
        TCIndex* Column3
        )
    {
        while(Column1 != Column1End && Column2 != Column2End) {
            TAIndex c1 = *Column1;
            TBIndex c2 = *Column2;

            if (c1 < c2) {
                if constexpr (TNeedOut) *Column3 = c1;
                ++Column1;
            } else if (c1 == c2) {
                if constexpr (TNeedOut) *Column3 = c1;
                ++Column1;
                ++Column2;
            } else {
                if constexpr (TNeedOut) *Column3 = c2;
                ++Column2;
            }
            ++Column3;
        }

        if constexpr (TNeedOut) {
            if (Column1 < Column1End) {
                return std::copy(Column1, Column1End, Column3);
            } else if (Column2 < Column2End) {
                return std::copy(Column2, Column2End, Column3);
            } else {
                return Column3;
            }
        } else {
            return Column3 + (Column1End - Column1) + (Column2End - Column2);
        }
    }

    /**
     * @brief This method is oriented to merge rows
     * @param rAlpha1 The coefficient of the first matrix
     * @param Column1 The index of the first matrix column
     * @param Column1End The last index of the first matrix column
     * @param Value1 The values of the first matrix
     * @param rAlpha2 The coefficient of the second matrix
     * @param Column2 The index of the second matrix column
     * @param Column2End The last index of the second matrix column
     * @param Value2 The values of the second matrix
     * @param Column3 The index of the third matrix column
     * @param Value3 The values of the third matrix
     * @return The resulting row
     */
    template <class TAValue,
              class TBIndex, class TBValue,
              class TCIndex, class TCValue,
              class TDIndex, class TDValue>
    static TCIndex* MergeRows(
        const TAValue& rAlpha1,
        const TBIndex* Column1,
        const TBIndex* Column1End,
        const TBValue* Value1,
        const TAValue& rAlpha2,
        const TCIndex* Column2,
        const TCIndex* Column2End,
        const TCValue* Value2,
        TDIndex* Column3,
        TDValue *Value3
        )
    {
        while(Column1 != Column1End && Column2 != Column2End) {
            TBIndex c1 = *Column1;
            TBIndex c2 = *Column2;

            if (c1 < c2) {
                ++Column1;

                *Column3 = c1;
                *Value3 = rAlpha1 * (*Value1++);
            } else if (c1 == c2) {
                ++Column1;
                ++Column2;

                *Column3 = c1;
                *Value3 = rAlpha1 * (*Value1++) + rAlpha2 * (*Value2++);
            } else {
                ++Column2;

                *Column3 = c2;
                *Value3 = rAlpha2 * (*Value2++);
            }

            ++Column3;
            ++Value3;
        }

        while(Column1 < Column1End) {
            *Column3++ = *Column1++;
            *Value3++ = rAlpha1 * (*Value1++);
        }

        while(Column2 < Column2End) {
            *Column3++ = *Column2++;
            *Value3++ = rAlpha2 * (*Value2++);
        }

        return Column3;
    }

    /**
     * @brief This method is oriented to multiply rows
     * @param AColumn The index of the first matrix column
     * @param AColumnEnd The last index of the first matrix column
     * @param BPtr The array constining the nonzero values per row of the second matrix
     * @param BColumn The index of the second matrix column
     * @param Column2End The last index of the second matrix column
     * @param Tmp1Column Indexes of the columns of first matrix
     * @param Tmp2Column Indexes of the columns of second matrix
     * @param Tmp3Column Indexes of the columns of third matrix
     * @return The resulting row
     */
    template <class TIndex>
    static TIndex ProdRowWidth(
        const TIndex* AColumn,
        const TIndex* AColumnEnd,
        const TIndex* BPtr,
        const TIndex* BColumn,
        TIndex* Tmp1Column,
        TIndex* Tmp2Column,
        TIndex* Tmp3Column
        )
    {
        const TIndex nrow = AColumnEnd - AColumn;

        /* No rows to merge, nothing to do */
        if (nrow == 0) return 0;

        /* Single row, just copy it to output */
        if (nrow == 1) return BPtr[*AColumn + 1] - BPtr[*AColumn];

        /* Two rows, merge them */
        if (nrow == 2) {
            int a1 = AColumn[0];
            int a2 = AColumn[1];

            return MergeRows<false>( BColumn + BPtr[a1], BColumn + BPtr[a1+1], BColumn + BPtr[a2], BColumn + BPtr[a2+1], Tmp1Column) - Tmp1Column;
        }

        /* Generic case (more than two rows).
        *
        * Merge rows by pairs, then merge the results together.
        * When merging two rows, the result is always wider (or equal).
        * Merging by pairs allows to work with short rows as often as possible.
        */
        // Merge first two.
        TIndex a1 = *AColumn++;
        TIndex a2 = *AColumn++;
        TIndex c_col1 = MergeRows<true>(  BColumn + BPtr[a1], BColumn + BPtr[a1+1], BColumn + BPtr[a2], BColumn + BPtr[a2+1], Tmp1Column ) - Tmp1Column;

        // Go by pairs.
        while(AColumn + 1 < AColumnEnd) {
            a1 = *AColumn++;
            a2 = *AColumn++;

            TIndex c_col2 = MergeRows<true>(  BColumn + BPtr[a1], BColumn + BPtr[a1+1], BColumn + BPtr[a2], BColumn + BPtr[a2+1], Tmp2Column ) - Tmp2Column;

            if (AColumn == AColumnEnd) {
                return MergeRows<false>( Tmp1Column, Tmp1Column + c_col1, Tmp2Column, Tmp2Column + c_col2, Tmp3Column ) - Tmp3Column;
            } else {
                c_col1 = MergeRows<true>( Tmp1Column, Tmp1Column + c_col1, Tmp2Column, Tmp2Column + c_col2, Tmp3Column ) - Tmp3Column;

                std::swap(Tmp1Column, Tmp3Column);
            }
        }

        // Merge the tail.
        a2 = *AColumn;
        return MergeRows<false>( Tmp1Column, Tmp1Column + c_col1, BColumn + BPtr[a2], BColumn + BPtr[a2+1], Tmp2Column ) - Tmp2Column;
    }

    /**
     * @brief This method is oriented to multiply rows
     * @param AColumn The index of the first matrix column
     * @param AColumnEnd The last index of the first matrix column
     * @param AValue The values of the first matrix
     * @param BPtr The array constining the nonzero values per row of the second matrix
     * @param BColumn The index of the second matrix column
     * @param BValue The values of the second matrix
     * @param OutColumn Indexes of the columns of output matrix
     * @param OutValue Values of the columns of output matrix
     * @param Tmp2Column Indexes of the columns of second matrix
     * @param Tmp2Value Values of the columns of second matrix
     * @param Tmp3Column Indexes of the columns of third matrix
     * @param Tmp3Value Values of the columns of third matrix
     * @return The resulting row
     */
    template <class TAIndex,
              class TAValue,
              class TBIndex,
              class TBValue,
              class TOutIndex,
              class TOutValue>
    static void ProdRow(
        const TAIndex* AColumn,
        const TAIndex* AColumnEnd,
        const TAValue *AValue,
        const TBIndex* BPtr,
        const TBIndex* BColumn,
        const TBValue *BValue,
        TOutIndex* OutColumn,
        TOutValue *OutValue,
        TOutIndex* Tmp2Column,
        TOutValue *Tmp2Value,
        TOutIndex* Tmp3Column,
        TOutValue *Tmp3Value
        )
    {
        const TAIndex nrow = AColumnEnd - AColumn;

        /* No rows to merge, nothing to do */
        if (nrow == 0) return;

        /* Single row, just copy it to output */
        if (nrow == 1) {
            TAIndex ac = *AColumn;
            TAValue av = *AValue;

            const TBValue *bv = BValue + BPtr[ac];
            const TBIndex* bc = BColumn + BPtr[ac];
            const TBIndex* be = BColumn + BPtr[ac+1];

            while(bc != be) {
                *OutColumn++ = *bc++;
                *OutValue++ = av * (*bv++);
            }

            return;
        }

        /* Two rows, merge them */
        if (nrow == 2) {
            TAIndex ac1 = AColumn[0];
            TAIndex ac2 = AColumn[1];

            TAValue av1 = AValue[0];
            TAValue av2 = AValue[1];

            MergeRows(av1, BColumn + BPtr[ac1], BColumn + BPtr[ac1+1], BValue + BPtr[ac1],
                      av2, BColumn + BPtr[ac2], BColumn + BPtr[ac2+1], BValue + BPtr[ac2],
                      OutColumn, OutValue);
        }

        /* Generic case (more than two rows).
        *
        * Merge rows by pairs, then merge the results together.
        * When merging two rows, the result is always wider (or equal).
        * Merging by pairs allows to work with short rows as often as possible.
        */
        // Merge first two.
        TAIndex ac1 = *AColumn++;
        TAIndex ac2 = *AColumn++;

        TAValue av1 = *AValue++;
        TAValue av2 = *AValue++;

        TOutIndex* tm1_col = OutColumn;
        TOutValue *tm1_val = OutValue;

        TOutIndex c_col1 = MergeRows(av1, BColumn + BPtr[ac1], BColumn + BPtr[ac1+1], BValue + BPtr[ac1],
                                     av2, BColumn + BPtr[ac2], BColumn + BPtr[ac2+1], BValue + BPtr[ac2],
                                     tm1_col, tm1_val ) - tm1_col;

        // Go by pairs.
        while(AColumn + 1 < AColumnEnd) {
            ac1 = *AColumn++;
            ac2 = *AColumn++;

            av1 = *AValue++;
            av2 = *AValue++;

            TOutIndex c_col2 = MergeRows( av1, BColumn + BPtr[ac1], BColumn + BPtr[ac1+1], BValue + BPtr[ac1], av2, BColumn + BPtr[ac2], BColumn + BPtr[ac2+1], BValue + BPtr[ac2], Tmp2Column, Tmp2Value ) - Tmp2Column;

            c_col1 = MergeRows(amgcl::math::identity<TAValue>(), tm1_col, tm1_col + c_col1, tm1_val,
                               amgcl::math::identity<TAValue>(), Tmp2Column, Tmp2Column + c_col2, Tmp2Value,
                               Tmp3Column, Tmp3Value) - Tmp3Column;

            std::swap(Tmp3Column, tm1_col);
            std::swap(Tmp3Value, tm1_val);
        }

        // Merge the tail if there is one.
        if (AColumn < AColumnEnd) {
            ac2 = *AColumn++;
            av2 = *AValue++;

            c_col1 = MergeRows(amgcl::math::identity<TAValue>(), tm1_col, tm1_col + c_col1, tm1_val,
                               av2, BColumn + BPtr[ac2], BColumn + BPtr[ac2+1], BValue + BPtr[ac2],
                               Tmp3Column, Tmp3Value ) - Tmp3Column;

            std::swap(Tmp3Column, tm1_col);
            std::swap(Tmp3Value, tm1_val);
        }

        // If we are lucky, tm1 now points to out.
        // Otherwise, copy the results.
        if (tm1_col != OutColumn) {
            std::copy(tm1_col, tm1_col + c_col1, OutColumn);
            std::copy(tm1_val, tm1_val + c_col1, OutValue);
        }

        return;
    }

    ///@}
}; // Class SparseMatrixMultiplicationUtility

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{


///@}

}  // namespace Kratos.