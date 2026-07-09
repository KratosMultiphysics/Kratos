// Project includes
#include "solving_strategies/builder_and_solvers/builder_and_solver.h" // BuilderAndSolver
#include "spaces/default_spaces.h" // TDefaultSparseSpace, TDefaultDenseSpace
#include "linear_solvers/linear_solver.h" // LinearSolver


namespace Kratos {


// Instantiated on the configure-time selected linear-algebra backend only:
// the two backends are mutually exclusive for the real sparse spaces.
template class KRATOS_API(KRATOS_CORE) BuilderAndSolver<TDefaultSparseSpace<double>,
                                                        TDefaultDenseSpace<double>,
                                                        LinearSolver<TDefaultSparseSpace<double>,
                                                                     TDefaultDenseSpace<double>>>;

template class KRATOS_API(KRATOS_CORE) BuilderAndSolver<TDefaultSparseSpace<float>,
                                                        TDefaultDenseSpace<double>,
                                                        LinearSolver<TDefaultSparseSpace<float>,
                                                                     TDefaultDenseSpace<double>>>;


} // namespace Kratos
