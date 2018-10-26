//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:             BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:   Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_TRILINOS_LINEAR_SOLVER_FACTORY_H_INCLUDED )
#define  KRATOS_TRILINOS_LINEAR_SOLVER_FACTORY_H_INCLUDED

// System includes

// External includes
#include "Epetra_FECrsMatrix.h"
#include "Epetra_FEVector.h"

// Project includes
#include "includes/define.h"
#include "trilinos_space.h"
#include "includes/linear_solver_factory.h"
#include "linear_solvers/linear_solver.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class TrilinosLinearSolverFactory
 * @ingroup KratosCore
 * @brief Here we add the functions needed for the registration of linear solvers
 * @details Defines the standard linear solver factory
 * @author Vicente Mataix Ferrandiz
 * @tparam TSparseSpace The sparse space definition
 * @tparam TLocalSpace The dense space definition
 * @tparam TLinearSolverType The linear solver type
 */
template <typename TSparseSpace, typename TLocalSpace, typename TLinearSolverType>
class TrilinosLinearSolverFactory
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
                                  const TrilinosLinearSolverFactory<TSparseSpace,TLocalSpace,TLinearSolverType>& rThis)
{
    rOStream << "TrilinosLinearSolverFactory" << std::endl;
    return rOStream;
}

///@}
///@name Input and output

void RegisterTrilinosLinearSolvers();

///@}

typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;

typedef LinearSolverFactory<TrilinosSparseSpaceType,  TrilinosLocalSpaceType> TrilinosLinearSolverFactoryType;

#ifdef KRATOS_REGISTER_LINEAR_SOLVER
#undef KRATOS_REGISTER_LINEAR_SOLVER
#endif
#define KRATOS_REGISTER_TRILINOS_LINEAR_SOLVER(name, reference) \
    KratosComponents<TrilinosLinearSolverFactoryType>::Add(name, reference);

}  // namespace Kratos.

#endif // KRATOS_TRILINOS_LINEAR_SOLVER_FACTORY_H_INCLUDED  defined
