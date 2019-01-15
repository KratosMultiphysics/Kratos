#ifndef KRATOS_AMGCL_MPI_SOLVE_FUNCTIONS_H
#define KRATOS_AMGCL_MPI_SOLVE_FUNCTIONS_H

#include <boost/property_tree/ptree.hpp>

namespace Kratos
{

// Spacialization of AMGCLScalarSolve for distribued systems.
template <class TSparseSpaceType>
typename std::enable_if<TSparseSpaceType::IsDistributed(), void>::type
AMGCLScalarSolve(
    typename TSparseSpaceType::MatrixType& rA,
    typename TSparseSpaceType::VectorType& rX,
    typename TSparseSpaceType::VectorType& rB,
    typename TSparseSpaceType::IndexType& rIterationNumber,
    double& rResidual,
    const boost::property_tree::ptree &amgclParams,
    int verbosity_level,
    bool use_gpgpu
    );

// Spacialization of AMGCLBlockSolve for distribued systems.
template <int TBlockSize, class TSparseSpaceType>
typename std::enable_if<TSparseSpaceType::IsDistributed(), void>::type
AMGCLBlockSolve(
    typename TSparseSpaceType::MatrixType & rA,
    typename TSparseSpaceType::VectorType& rX,
    typename TSparseSpaceType::VectorType& rB,
    typename TSparseSpaceType::IndexType& rIterationNumber,
    double& rResidual,
    boost::property_tree::ptree amgclParams,
    int verbosity_level,
    bool use_gpgpu
    );

} // namespace Kratos

#endif
