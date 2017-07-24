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

#if !defined(KRATOS_LINEAR_SOLVER_FACTORY_H_INCLUDED )
#define  KRATOS_LINEAR_SOLVER_FACTORY_H_INCLUDED



// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/preconditioner.h"
#include "linear_solvers/scaling_solver.h"

namespace Kratos
{


template <typename TSparseSpace, typename TlocalSpace, typename TLinearSolverType>
class LinearSolverFactory : public LinearSolverFactoryBase<TSparseSpace,TlocalSpace>
{
protected:

    typename LinearSolver<TSparseSpace,TlocalSpace>::Pointer CreateHelper(Kratos::Parameters settings) const
    {

        if(settings.Has("scaling") && settings["scaling"].GetBool() == true)
        {
            auto pinner_solver = typename LinearSolver<TSparseSpace,TlocalSpace>::Pointer(new TLinearSolverType(settings));

            return typename LinearSolver<TSparseSpace,TlocalSpace>::Pointer(
                       new ScalingSolver<TSparseSpace,TlocalSpace>(pinner_solver, true));

        }
        else
            return typename LinearSolver<TSparseSpace,TlocalSpace>::Pointer(new TLinearSolverType(settings));
    }
};

/// output stream function
template <typename TSparseSpace, typename TlocalSpace, typename TLinearSolverType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const LinearSolverFactory<TSparseSpace,TlocalSpace,TLinearSolverType>& rThis)
{
    rOStream << "LinearSolverFactory" << std::endl;
    return rOStream;
}
















void RegisterPreconditioners();
void RegisterLinearSolvers();





}  // namespace Kratos.

#endif // KRATOS_LINEAR_SOLVER_FACTORY_H_INCLUDED  defined 
