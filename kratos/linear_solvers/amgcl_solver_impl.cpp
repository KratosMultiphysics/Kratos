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

// STD includes
#include <optional>


namespace Kratos {


template <class TValue>
struct AMGCLAdaptor<TUblasSparseSpace<TValue>>
{
    template <int BlockSize>
    auto MakeMatrixAdaptor(const typename TUblasSparseSpace<TValue>::MatrixType& rMatrix)
    {
        if constexpr (BlockSize == 1) {
            KRATOS_TRY
            return amgcl::adapter::zero_copy(rMatrix.size1(),
                                             rMatrix.index1_data().begin(),
                                             rMatrix.index2_data().begin(),
                                             rMatrix.value_data().begin());
            KRATOS_CATCH("")
        } else {
            using BlockType = amgcl::static_matrix<
                TValue,
                BlockSize,
                BlockSize
            >;

            KRATOS_TRY
            // amgcl::adapter::block_matrix constructs a class
            // that stores a reference to the "matrix" passed
            // into it, which in this case means the tuple
            // defined below. We need to keep it alive until
            // the hierarchy construction finishes, hence the
            // convoluted member variable.
            mIntermediateAdaptor.emplace(std::tuple_cat(
                std::tuple<std::size_t>(rMatrix.size1()),
                std::tie(rMatrix.index1_data(),
                         rMatrix.index2_data(),
                         rMatrix.value_data())
            ));

            return amgcl::adapter::block_matrix<BlockType>(mIntermediateAdaptor.value());
            KRATOS_CATCH("")
        }
    }

    template <class TStaticMatrix>
    std::size_t BlockSystemSize(const typename TUblasSparseSpace<TValue>::MatrixType& rMatrix) const noexcept
    {
        return TUblasSparseSpace<TValue>::Size1(rMatrix) / AMGCLStaticVectorTraits<TStaticMatrix>::value;
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

private:
    using TMatrix = typename TUblasSparseSpace<TValue>::MatrixType;

    std::optional<std::tuple<
        std::size_t,
        const typename TMatrix::index_array_type&,
        const typename TMatrix::index_array_type&,
        const typename TMatrix::value_array_type&
    >> mIntermediateAdaptor;
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
