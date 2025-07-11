//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//

// Core includes
#include "linear_solvers/amgcl_solver.h" // AMGCLSolver
#include "linear_solvers/amgcl_solver_impl.hpp" // AMGCLAdaptor
#include "spaces/ublas_space.h" // TUblasSparseSpace, TUblasDenseSpace


namespace Kratos {


template <class TValue>
struct AMGCLAdaptor<TUblasSparseSpace<TValue>>
{
    auto MakeMatrixAdaptor(const typename TUblasSparseSpace<TValue>::MatrixType& rMatrix) const
    {
        return amgcl::adapter::zero_copy(rMatrix.size1(),
                                         rMatrix.index1_data().begin(),
                                         rMatrix.index2_data().begin(),
                                         rMatrix.value_data().begin());
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


template class AMGCLSolver<TUblasSparseSpace<double>,
                           TUblasDenseSpace<double>,
                           Reorderer<TUblasSparseSpace<double>,
                                     TUblasDenseSpace<double>>>;

template class AMGCLSolver<TUblasSparseSpace<float>,
                           TUblasDenseSpace<double>,
                           Reorderer<TUblasSparseSpace<float>,
                                     TUblasDenseSpace<double>>>;


} // namespace Kratos
