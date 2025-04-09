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
    virtual bool Has(const std::string SolverType) const
    {
        return KratosComponents< FactoryType >::Has( SolverType );
    }

    /**
     * @brief This method creates a new solver
     * @return The pointer to the solver of interest
     */
    virtual typename LinearSolver<TSparseSpace,TLocalSpace>::Pointer Create(Kratos::Parameters Settings) const
    {
        std::string solver_name = Settings["solver_type"].GetString();
        // remove name of the application (if passed)
        // e.g. "LinearSolversApplication.sparse_lu" => "sparse_lu"
        solver_name = solver_name.substr(solver_name.find('.') + 1);

        KRATOS_ERROR_IF_NOT(Has(solver_name))
            << "Trying to construct a Linear solver with solver_type:\n\""
            << solver_name << "\" which does not exist.\n"
            << "The list of available options (for currently loaded applications) is:\n"
            << KratosComponents< FactoryType >() << std::endl;

        const auto& aux = KratosComponents< FactoryType >::Get(solver_name);
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

KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<LinearSolverFactoryType>;
KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<LinearSolverFactory<
    TUblasSparseSpace<float>,
    TUblasDenseSpace<double>
>>;

#ifdef KRATOS_REGISTER_LINEAR_SOLVER
#undef KRATOS_REGISTER_LINEAR_SOLVER
#endif
#define KRATOS_REGISTER_LINEAR_SOLVER(name, reference) \
    KratosComponents<LinearSolverFactoryType>::Add(name, reference);

typedef TUblasSparseSpace<std::complex<double>> ComplexSparseSpaceType;
typedef TUblasDenseSpace<std::complex<double>> ComplexLocalSparseSpaceType;

typedef LinearSolverFactory<ComplexSparseSpaceType,  ComplexLocalSparseSpaceType> ComplexLinearSolverFactoryType;

KRATOS_API_EXTERN template class KRATOS_API(KRATOS_CORE) KratosComponents<ComplexLinearSolverFactoryType>;

#ifdef KRATOS_REGISTER_COMPLEX_LINEAR_SOLVER
#undef KRATOS_REGISTER_COMPLEX_LINEAR_SOLVER
#endif
#define KRATOS_REGISTER_COMPLEX_LINEAR_SOLVER(name, reference) \
    KratosComponents<ComplexLinearSolverFactoryType>::Add(name, reference);


}  // namespace Kratos.

#endif // KRATOS_LINEAR_SOLVER_FACTORY_H_INCLUDED  defined
