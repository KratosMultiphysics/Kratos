// Project includes
#include "solving_strategies/builder_and_solvers/builder_and_solver.h" // BuilderAndSolver
#include "spaces/ublas_space.h" // TUblasSparseSpace, TUblasDenseSpace
#include "linear_solvers/linear_solver.h" // LinearSolver


namespace Kratos {


template class KRATOS_API(KRATOS_CORE) BuilderAndSolver<TUblasSparseSpace<double>,
                                                        TUblasDenseSpace<double>,
                                                        LinearSolver<TUblasSparseSpace<double>,
                                                                     TUblasDenseSpace<double>>>;

template class KRATOS_API(KRATOS_CORE) BuilderAndSolver<TUblasSparseSpace<float>,
                                                        TUblasDenseSpace<double>,
                                                        LinearSolver<TUblasSparseSpace<float>,
                                                                     TUblasDenseSpace<double>>>;


} // namespace Kratos
