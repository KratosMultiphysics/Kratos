//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//
//

#if !defined(KRATOS_STANDARD_BUILDER_AND_SOLVER_FACTORY_H_INCLUDED )
#define  KRATOS_STANDARD_BUILDER_AND_SOLVER_FACTORY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "factories/base_factory.h"

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
 * @class StandardBuilderAndSolverFactory
 * @ingroup KratosCore
 * @brief Here we add the functions needed for the registration of builder and solvers
 * @details Defines the standard builder and solver factory
 * @author Vicente Mataix Ferrandiz
 * @tparam TBuilderAndSolverType The builder and solver type
 * @tparam TLinearSolver The linear solver type considered
 * @tparam TCustomBuilderAndSolverType The builder and solver type (derived class)
 */
template <typename TBuilderAndSolverType, typename TLinearSolver, typename TCustomBuilderAndSolverType>
class StandardBuilderAndSolverFactory
    : public BaseFactory<TBuilderAndSolverType, TLinearSolver>
{
    ///@name Type Definitions
    ///@{

    // The definition of the linear solver type
    typedef TLinearSolver LinearSolverType;

    /// The definition of the builder and solver
    typedef TBuilderAndSolverType BuilderAndSolverType;

    ///@}
protected:

    ///@name Operations
    ///@{

    /**
     * @brief This method is an auxiliar method to create a new builder and solver
     * @return The pointer to the builder and solver of interest
     */
    typename BuilderAndSolverType::Pointer CreateClass(typename LinearSolverType::Pointer pLinearSolver, Kratos::Parameters Settings) const override
    {
        return typename BuilderAndSolverType::Pointer(new TCustomBuilderAndSolverType(pLinearSolver, Settings));
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
template <typename TBuilderAndSolverType, typename TLinearSolver, typename TCustomBuilderAndSolverType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const StandardBuilderAndSolverFactory<TBuilderAndSolverType, TLinearSolver, TCustomBuilderAndSolverType>& rThis)
{
    rOStream << "StandardBuilderAndSolverFactory" << std::endl;
    return rOStream;
}

///@}
///@name Input and output

void RegisterBuilderAndSolvers();

///@}

}  // namespace Kratos.

#endif // KRATOS_STANDARD_BUILDER_AND_SOLVER_FACTORY_H_INCLUDED  defined
