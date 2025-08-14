//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// Core includes
#include "linear_solvers/amgcl_solver.h" // AMGCLSolver
#include "linear_solvers/amgcl_solver_impl.hpp" // AMGCLAdaptor
#include "spaces/ublas_space.h" // TUblasSparseSpace, TUblasDenseSpace


namespace Kratos {


template <class TValue>
struct AMGCLAdaptor<TUblasSparseSpace<TValue>>
{
    template <int BlockSize>
    auto MakeMatrixAdaptor(const typename TUblasSparseSpace<TValue>::MatrixType& rMatrix) const
    {
        if constexpr (BlockSize == 1) {
            return amgcl::adapter::zero_copy(rMatrix.size1(),
                                             rMatrix.index1_data().begin(),
                                             rMatrix.index2_data().begin(),
                                             rMatrix.value_data().begin());
        } else {
            using BackendMatrix = amgcl::static_matrix<
                TValue,
                BlockSize,
                BlockSize
            >;
            return amgcl::adapter::block_matrix<BackendMatrix>(
                std::tuple_cat(
                    std::make_tuple(rMatrix.size1()),
                    std::tie(rMatrix.index1_data(), rMatrix.index2_data(), rMatrix.value_data())
                )
            );
        }
    }

    auto MakeVectorIterator(const typename TUblasSparseSpace<TValue>::VectorType& rVector) const
    {
        KRATOS_ERROR_IF(rVector.empty());
        return &*rVector.begin();
    }

    auto MakeVectorIterator(typename TUblasSparseSpace<TValue>::VectorType& rVector) const
    {
        KRATOS_ERROR_IF(rVector.empty());
        return &*rVector.begin();
    }
}; // struct AMGCLAdaptor


template class KRATOS_API(KRATOS_CORE) AMGCLSolver<
    TUblasSparseSpace<double>,
    TUblasDenseSpace<double>
>;

template class KRATOS_API(KRATOS_CORE) AMGCLSolver<
    TUblasSparseSpace<float>,
    TUblasDenseSpace<double>
>;


} // namespace Kratos
