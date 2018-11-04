/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Authors: Thomas Oberbichler
*/

#if !defined(KRATOS_UBLAS_WRAPPER_H_INCLUDED)
#define KRATOS_UBLAS_WRAPPER_H_INCLUDED

// System includes
#include <vector>

// External includes
#include <Eigen/Core>
#include <Eigen/Sparse>

// Project includes
#include "includes/define.h"


namespace Kratos
{

template <typename scalar_t, typename TEigenSparseMatrix = Eigen::SparseMatrix<scalar_t, Eigen::RowMajor, int>>
class UblasWrapper
{
    const std::vector<int> m_index1;
    const std::vector<int> m_index2;
    const Eigen::Map<const TEigenSparseMatrix> m_map;

public:

    template <typename TUblasSparseMatrix>
    UblasWrapper(const TUblasSparseMatrix& matrix)
    : m_index1(matrix.index1_data().begin(), matrix.index1_data().end()),
      m_index2(matrix.index2_data().begin(), matrix.index2_data().end()),
      m_map(matrix.size1(), matrix.size2(), matrix.nnz(), m_index1.data(), m_index2.data(), matrix.value_data().begin())
    {
    }

    const Eigen::Map<const TEigenSparseMatrix>& matrix() const
    {
        return m_map;
    }
};

} // namespace Kratos

#endif // defined(KRATOS_UBLAS_WRAPPER_H_INCLUDED)