#pragma once

namespace Kratos
{

template <typename TEigenSparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor, int>>
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

}