// Project includes
#include "solving_strategies/builder_and_solvers/builder_and_solver.h" // BuilderAndSolver
#include "spaces/ublas_space.h" // TUblasSparseSpace, TUblasDenseSpace
#include "spaces/default_spaces.h" // TDefaultSparseSpace, TDefaultDenseSpace
#include "linear_solvers/linear_solver.h" // LinearSolver
#ifdef KRATOS_USE_EIGEN_BACKEND
#include "linear_solvers/linear_solver_eigen.h" // Already-exported LinearSolver<TEigenSparseSpace<double>, TUblasDenseSpace<double>>
#endif


namespace Kratos {


template class KRATOS_API(KRATOS_CORE) BuilderAndSolver<TUblasSparseSpace<double>,
                                                        TUblasDenseSpace<double>,
                                                        LinearSolver<TUblasSparseSpace<double>,
                                                                     TUblasDenseSpace<double>>>;

template class KRATOS_API(KRATOS_CORE) BuilderAndSolver<TUblasSparseSpace<float>,
                                                        TUblasDenseSpace<double>,
                                                        LinearSolver<TUblasSparseSpace<float>,
                                                                     TUblasDenseSpace<double>>>;

#ifdef KRATOS_USE_EIGEN_BACKEND
// With the Eigen backend the default sparse space (used e.g. by
// ResidualBasedBlockBuilderAndSolver<DefaultSparseSpaceType, ...> in
// add_strategies_to_python.cpp and in the uBLAS/Eigen parity tests) differs
// from the uBLAS one explicitly instantiated above, so it needs its own
// explicit instantiation to be importable from other binaries/DLLs.
template class KRATOS_API(KRATOS_CORE) BuilderAndSolver<TDefaultSparseSpace<double>,
                                                        TDefaultDenseSpace<double>,
                                                        LinearSolver<TDefaultSparseSpace<double>,
                                                                     TDefaultDenseSpace<double>>>;
#endif


} // namespace Kratos
