#ifndef KRATOS_AMGCL_MPI_SOLVE_FUNCTIONS_H
#define KRATOS_AMGCL_MPI_SOLVE_FUNCTIONS_H

#include <boost/property_tree/ptree.hpp>

namespace Kratos
{

KRATOS_API(TRILINOS_APPLICATION) void AMGCLSolve(
    int block_size,
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>::MatrixType& rA,
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>::VectorType& rX,
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>::VectorType& rB,
    TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>::IndexType& rIterationNumber,
    double& rResidual,
    boost::property_tree::ptree amgclParams,
    int verbosity_level,
    bool use_gpgpu
    );

} // namespace Kratos

#endif
