//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//
//

#if !defined(KRATOS_STANDARD_LINEAR_SOLVER_FACTORY_H_INCLUDED )
#define  KRATOS_STANDARD_LINEAR_SOLVER_FACTORY_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/linear_solver_factory.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/preconditioner.h"
#include "linear_solvers/scaling_solver.h"

namespace Kratos
{


template <typename TSparseSpace, typename TLocalSpace, typename TLinearSolverType>
class StandardLinearSolverFactory : public LinearSolverFactory<TSparseSpace,TLocalSpace>
{
protected:

    typename LinearSolver<TSparseSpace,TLocalSpace>::Pointer CreateHelper(Kratos::Parameters settings) const override
    {

        if(settings.Has("scaling") && settings["scaling"].GetBool() == true)
        {
            auto pinner_solver = typename LinearSolver<TSparseSpace,TLocalSpace>::Pointer(new TLinearSolverType(settings));

            return typename LinearSolver<TSparseSpace,TLocalSpace>::Pointer(
                       new ScalingSolver<TSparseSpace,TLocalSpace>(pinner_solver, true));

        }
        else
            return typename LinearSolver<TSparseSpace,TLocalSpace>::Pointer(new TLinearSolverType(settings));
    }
};

/// output stream function
template <typename TSparseSpace, typename TLocalSpace, typename TLinearSolverType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const StandardLinearSolverFactory<TSparseSpace,TLocalSpace,TLinearSolverType>& rThis)
{
    rOStream << "StandardLinearSolverFactory" << std::endl;
    return rOStream;
}

void RegisterPreconditioners();
void RegisterLinearSolvers();

}  // namespace Kratos.

#endif // KRATOS_STANDARD_LINEAR_SOLVER_FACTORY_H_INCLUDED  defined
