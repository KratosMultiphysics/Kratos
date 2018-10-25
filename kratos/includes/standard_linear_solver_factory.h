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
#include "linear_solvers/scaling_solver.h"

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
 * @class StandardLinearSolverFactory
 * @ingroup KratosCore
 * @brief Here we add the functions needed for the registration of linear solvers
 * @details Defines the standard linear solver factory
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 * @tparam TSparseSpace The sparse space definition
 * @tparam TLocalSpace The dense space definition
 * @tparam TLinearSolverType The linear solver type
 */
template <typename TSparseSpace, typename TLocalSpace, typename TLinearSolverType>
class StandardLinearSolverFactory
    : public LinearSolverFactory<TSparseSpace,TLocalSpace>
{
    ///@name Type Definitions
    ///@{

    /// The definition of the linear solver
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
        if(settings.Has("scaling") && settings["scaling"].GetBool()) {
            auto pinner_solver = typename LinearSolverType::Pointer(new TLinearSolverType(settings));

            return typename LinearSolverType::Pointer(new ScalingSolver<TSparseSpace,TLocalSpace>(pinner_solver, true));

        } else
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
                                  const StandardLinearSolverFactory<TSparseSpace,TLocalSpace,TLinearSolverType>& rThis)
{
    rOStream << "StandardLinearSolverFactory" << std::endl;
    return rOStream;
}

///@}
///@name Input and output

void RegisterLinearSolvers();

///@}

}  // namespace Kratos.

#endif // KRATOS_STANDARD_LINEAR_SOLVER_FACTORY_H_INCLUDED  defined
