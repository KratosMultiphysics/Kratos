//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher), based on the work of Jordi Cotela
//

#ifndef KRATOS_MPI_LINEAR_SOLVER_FACTORY_H_INCLUDED
#define KRATOS_MPI_LINEAR_SOLVER_FACTORY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "mpi/spaces/amgcl_mpi_space.h"
#include "factories/linear_solver_factory.h"


namespace Kratos
{
///@name Kratos Classes
///@{

/**
 * @class MPILinearSolverFactory
 * @ingroup KratosMPICore
 * @brief Here we add the functions needed for the registration of linear solvers
 * @details Defines the standard linear solver factory
 * @author Philipp Bucher
 * @tparam TSparseSpace The sparse space definition
 * @tparam TLocalSpace The dense space definition
 * @tparam TLinearSolverType The linear solver type
 */
template <typename TSparseSpace, typename TLocalSpace, typename TLinearSolverType>
class MPILinearSolverFactory
    : public LinearSolverFactory<TSparseSpace,TLocalSpace>
{
    ///@name Type Definitions
    ///@{

    /// The definition of the preconditioner
    typedef LinearSolver<TSparseSpace,TLocalSpace> LinearSolverType;

    ///@}
protected:
    ///@name Protected Operators
    ///@{

    /**
     * @brief This method is an auxiliar method to create a new solver
     * @return The pointer to the solver of interest
     */
    typename LinearSolverType::Pointer CreateSolver(Kratos::Parameters settings) const override
    {
        return typename LinearSolverType::Pointer(new TLinearSolverType(settings));
    }
    ///@}
};

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
template <typename TSparseSpace, typename TLocalSpace, typename TLinearSolverType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MPILinearSolverFactory<TSparseSpace,TLocalSpace,TLinearSolverType>& rThis)
{
    rOStream << "MPILinearSolverFactory" << std::endl;
    return rOStream;
}

///@}
///@name Input and output

///@}

void RegisterMPILinearSolvers();

typedef amgcl::backend::builtin<double> Backend;
typedef amgcl::mpi::distributed_matrix<Backend> amgcl_mpi_matrix;
typedef typename Backend::vector amgcl_mpi_vector;

typedef AmgclMPISpace<amgcl_mpi_matrix, amgcl_mpi_vector> MPISparseSpaceType;
typedef UblasSpace<double, Matrix, Vector> MPILocalSpaceType;

typedef LinearSolverFactory<MPISparseSpaceType,  MPILocalSpaceType> MPILinearSolverFactoryType;

#ifdef KRATOS_REGISTER_MPI_LINEAR_SOLVER
#undef KRATOS_REGISTER_MPI_LINEAR_SOLVER
#endif
#define KRATOS_REGISTER_MPI_LINEAR_SOLVER(name, reference) ; \
    KratosComponents<MPILinearSolverFactoryType>::Add(name, reference);

KRATOS_API_EXTERN template class KratosComponents<MPILinearSolverFactoryType>;

}  // namespace Kratos.

#endif // KRATOS_MPI_LINEAR_SOLVER_FACTORY_H_INCLUDED  defined
