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

#if !defined(KRATOS_BUILDER_AND_SOLVER_FACTORY_H_INCLUDED )
#define  KRATOS_BUILDER_AND_SOLVER_FACTORY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/kratos_components.h"
#include "includes/shared_pointers.h"
#include "linear_solvers/linear_solver.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "spaces/ublas_space.h"

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
 * @class BuilderAndSolverFactory
 * @ingroup KratosCore
 * @brief Here we add the functions needed for the registration of builder and solvers
 * @details Defines the base builder and solver factory
 * @author Vicente Mataix Ferrandiz
 * @tparam TSparseSpace The sparse space definition
 * @tparam TLocalSpace The dense space definition
 * @tparam TLinearSolver The linear solver type considered
 */
template< typename TSparseSpace, typename TLocalSpace, class TLinearSolver>
class BuilderAndSolverFactory
{
public:
    ///@name Type Definitions
    ///@{

    /// The definition of the custom class
    typedef BuilderAndSolverFactory<TSparseSpace,TLocalSpace, TLinearSolver> FactoryType;

    /// The definition of the builder and solver
    typedef BuilderAndSolver<TSparseSpace,TLocalSpace, TLinearSolver> BuilderAndSolverType;

    /// Pointer definition of BuilderAndSolverFactory
    KRATOS_CLASS_POINTER_DEFINITION(BuilderAndSolverFactory );

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method checks if the linear builder and solver is registered
     * @return True if registered, false otherwise
     */
    virtual bool Has(const std::string& rSolverType)
    {
        return KratosComponents< FactoryType >::Has( rSolverType );
    }

    /**
     * @brief This method creates a new builder and solver
     * @return The pointer to the builder and solver of interest
     */
    virtual typename BuilderAndSolverType::Pointer Create(Kratos::Parameters Settings)
    {
        const std::string& builder_and_solver_type = Settings["builder_and_solver_type"].GetString();
        if(Has( builder_and_solver_type ) == false) {
            KRATOS_ERROR << "Trying to construct a builder and solver with type builder_and_solver_type= " << builder_and_solver_type << std::endl <<
                            "Which does not exist. The list of available options (for currently loaded applications) is: " << std::endl <<
                         KratosComponents< FactoryType >() << std::endl;
        }
        const auto& aux = KratosComponents< FactoryType >::Get( builder_and_solver_type );
        return aux.CreateBuilderAndSolver(Settings);
    }
    ///@}
protected:
    ///@name Protected Operators
    ///@{

    /**
     * @brief This method is an auxiliar method to create a new builder and solver
     * @return The pointer to the builder and solver of interest
     */
    virtual typename BuilderAndSolverType::Pointer CreateBuilderAndSolver(Kratos::Parameters Settings)  const
    {
        KRATOS_ERROR << "calling the base class BuilderAndSolverFactory" << std::endl;
    }
};

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// output stream function
template< typename TSparseSpace, typename TLocalSpace, class TLinearSolver>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const BuilderAndSolverFactory<TSparseSpace, TLocalSpace, TLinearSolver>& rThis)
{
    rOStream << "BuilderAndSolverFactory" << std::endl;

    return rOStream;
}

///@}

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSparseSpaceType;
typedef LinearSolver<SparseSpaceType,LocalSparseSpaceType> LinearSolverType;

typedef BuilderAndSolverFactory<SparseSpaceType, LocalSparseSpaceType, LinearSolverType> BuilderAndSolverFactoryType;

#ifdef KRATOS_REGISTER_BUILDER_AND_SOLVER
#undef KRATOS_REGISTER_BUILDER_AND_SOLVER
#endif
#define KRATOS_REGISTER_BUILDER_AND_SOLVER(name, reference) \
    KratosComponents<BuilderAndSolverFactoryType>::Add(name, reference);

}  // namespace Kratos.

#endif // KRATOS_BUILDER_AND_SOLVER_FACTORY_H_INCLUDED  defined
