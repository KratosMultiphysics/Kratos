//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_LINEAR_SOLVER_FACTORY_H_INCLUDED )
#define  KRATOS_LINEAR_SOLVER_FACTORY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/kratos_components.h"
#include "includes/shared_pointers.h"
#include "linear_solvers/linear_solver.h"
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
 * @class LinearSolverFactory
 * @ingroup KratosCore
 * @brief Here we add the functions needed for the registration of linear solvers
 * @details Defines the base linear solver factory
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 * @tparam TSparseSpace The sparse space definition
 * @tparam TLocalSpace The dense space definition
 */
template< typename TSparseSpace, typename TLocalSpace>
class LinearSolverFactory
{
public:
    ///@name Type Definitions
    ///@{

    /// The definition of the custom class
    typedef LinearSolverFactory<TSparseSpace,TLocalSpace> FactoryType;

    /// Pointer definition of LinearSolverFactory
    KRATOS_CLASS_POINTER_DEFINITION(LinearSolverFactory );

    ///@}
    ///@name Life Cycle
    ///@{

    virtual ~LinearSolverFactory(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method checks if the linear solver is registered
     * @return True if registered, false otherwise
     */
    virtual bool Has(const std::string SolverType)
    {
        return KratosComponents< FactoryType >::Has( SolverType );
    }

    /**
     * @brief This method creates a new solver
     * @return The pointer to the solver of interest
     */
    virtual typename LinearSolver<TSparseSpace,TLocalSpace>::Pointer Create(Kratos::Parameters Settings)
    {
        if(KratosComponents< FactoryType >::Has( Settings["solver_type"].GetString())== false) {
            KRATOS_ERROR << "Trying to construct a Linear solver with solver_type = " << Settings["solver_type"].GetString() << std::endl << "which does not exist. The list of available options (for currently loaded applications) is: " << std::endl <<
                         KratosComponents< FactoryType >() << std::endl;
        }
        const auto& aux = KratosComponents< FactoryType >::Get( Settings["solver_type"].GetString()  );
        return aux.CreateSolver( Settings );
    }

    ///@}
protected:
    ///@name Protected Operators
    ///@{

    /**
     * @brief This method is an auxiliar method to create a new solver
     * @return The pointer to the solver of interest
     */
    virtual typename LinearSolver<TSparseSpace,TLocalSpace>::Pointer CreateSolver(Kratos::Parameters Settings)  const
    {
        KRATOS_ERROR << "Calling the base class LinearSolverFactory" << std::endl;
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
template< typename TSparseSpace, typename TLocalSpace>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const LinearSolverFactory<TSparseSpace, TLocalSpace>& rThis)
{
    rOStream << "LinearSolverFactory" << std::endl;

    return rOStream;
}
///@}

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSparseSpaceType;

typedef LinearSolverFactory<SparseSpaceType,  LocalSparseSpaceType> LinearSolverFactoryType;

#ifdef KRATOS_REGISTER_LINEAR_SOLVER
#undef KRATOS_REGISTER_LINEAR_SOLVER
#endif
#define KRATOS_REGISTER_LINEAR_SOLVER(name, reference) \
    KratosComponents<LinearSolverFactoryType>::Add(name, reference);

typedef TUblasSparseSpace<std::complex<double>> ComplexSparseSpaceType;
typedef TUblasDenseSpace<std::complex<double>> ComplexLocalSparseSpaceType;

typedef LinearSolverFactory<ComplexSparseSpaceType,  ComplexLocalSparseSpaceType> ComplexLinearSolverFactoryType;

#ifdef KRATOS_REGISTER_COMPLEX_LINEAR_SOLVER
#undef KRATOS_REGISTER_COMPLEX_LINEAR_SOLVER
#endif
#define KRATOS_REGISTER_COMPLEX_LINEAR_SOLVER(name, reference) \
    KratosComponents<ComplexLinearSolverFactoryType>::Add(name, reference);

}  // namespace Kratos.

#endif // KRATOS_LINEAR_SOLVER_FACTORY_H_INCLUDED  defined
