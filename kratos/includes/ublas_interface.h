//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

#pragma once

// System includes
#include <algorithm>
#include <cstddef>
#include <string>
#include <iostream>

// External includes
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/hermitian.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp>

// Project includes
#include "includes/define.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
    using namespace boost::numeric::ublas;

    template <typename TDataType> using DenseMatrix = boost::numeric::ublas::matrix<TDataType>;
    template <typename TDataType> using DenseVector = boost::numeric::ublas::vector<TDataType>;

    template <typename TDataType, std::size_t TSize1, std::size_t TSize2> using BoundedMatrix = boost::numeric::ublas::bounded_matrix<TDataType, TSize1, TSize2>;
    template <typename TDataType, std::size_t TSize> using BoundedVector = boost::numeric::ublas::bounded_vector<TDataType, TSize>;


    typedef boost::numeric::ublas::vector<double> Vector;
    typedef unit_vector<double> UnitVector;
    typedef zero_vector<double> ZeroVector;
    typedef scalar_vector<double> ScalarVector;
    //typedef sparse_vector<double> SparseVector;
    typedef mapped_vector<double> SparseVector;

    typedef compressed_vector<double> CompressedVector;
    typedef coordinate_vector<double> CoordinateVector;
    typedef vector_range<Vector> VectorRange;
    typedef vector_slice<Vector> VectorSlice;

    typedef matrix<double> Matrix;
    typedef identity_matrix<double> IdentityMatrix;
    typedef zero_matrix<double> ZeroMatrix;
    typedef scalar_matrix<double> ScalarMatrix;
    typedef triangular_matrix<double> TriangularMatrix;
    typedef symmetric_matrix<double> SymmetricMatrix;
    typedef hermitian_matrix<double> HermitianMatrix;
    typedef banded_matrix<double> BandedMatrix;
    //typedef sparse_matrix<double> SparseMatrix;
    typedef mapped_matrix<double> SparseMatrix;
    typedef coordinate_matrix<double> CoordinateMatrix;
    typedef matrix_column<Matrix> MatrixColumn;
    typedef matrix_vector_range<Matrix> MatrixVectorRange;
    typedef matrix_vector_slice<Matrix> MatrixVectorSlice;
    typedef matrix_range<Matrix> MatrixRange;
    typedef matrix_slice<Matrix> MatrixSlice;

	template <typename TExpressionType> using MatrixRow = matrix_row<TExpressionType>;

    typedef boost::numeric::ublas::compressed_matrix<double> CompressedMatrix;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

    /**
     * @brief Set every entry of a dense container to a given value.
     * @details Works for ublas::vector, ublas::matrix, ublas::bounded_matrix
     *          and Kratos::array_1d. Intended to be called right after a
     *          size-only construction, e.g. Vector a(n); SetValues(a, 0.0);
     * @tparam TContainerType Dense container type.
     * @param rContainer Container whose entries are overwritten.
     * @param rValue Value assigned to every entry.
     */
    template <class TContainerType>
    void SetValues(
        TContainerType& rContainer,
        const typename TContainerType::value_type& rValue)
    {
        std::fill(rContainer.data().begin(), rContainer.data().end(), rValue);
    }

    /**
     * @brief Construct a dense container from its sizes and set every entry to a given value.
     * @details Intended for const declarations, which cannot be filled after construction:
     *          const Vector a = CreateFilled<Vector>(1.0, n);
     * @tparam TContainerType Dense container type.
     * @tparam TSizes Types of the size arguments forwarded to the size constructor.
     * @param rValue Value assigned to every entry.
     * @param rSizes Sizes forwarded to the container constructor.
     * @return The filled container.
     */
    template <class TContainerType, class... TSizes>
    TContainerType CreateFilled(
        const typename TContainerType::value_type& rValue,
        const TSizes... rSizes)
    {
        TContainerType container(rSizes...);
        SetValues(container, rValue);
        return container;
    }

///@}
///@name Kratos Classes
///@{


///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}
}  // namespace Kratos.