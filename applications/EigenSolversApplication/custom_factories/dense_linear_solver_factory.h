/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Author: Quirin Aumann
*/

#if !defined(KRATOS_DENSE_LINEAR_SOLVER_FACTORY_H_INCLUDED )
#define  KRATOS_DENSE_LINEAR_SOLVER_FACTORY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "factories/linear_solver_factory.h"
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
 * @class DenseLinearSolverFactory
 * @ingroup EigenSolversApplication
 * @brief Here we add the functions needed for the registration of dense linear solvers
 * @details Defines the dense linear solver factory
 * @author Quirin Aumann
 * @tparam TGlobalSpace The global space definition
 * @tparam TLocalSpace The local space definition
 * @tparam TLinearSolverType The linear solver type
 */
template <typename TGlobalSpace, typename TLocalSpace, typename TLinearSolverType>
class KRATOS_API(EIGENSOLVERS_APPLICATION) DenseLinearSolverFactory
    : public LinearSolverFactory<TGlobalSpace,TLocalSpace>
{
    ///@name Type Definitions
    ///@{

    /// The definition of the preconditioner
    typedef LinearSolver<TGlobalSpace,TLocalSpace> LinearSolverType;

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
template <typename TGlobalSpace, typename TLocalSpace, typename TLinearSolverType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DenseLinearSolverFactory<TGlobalSpace,TLocalSpace,TLinearSolverType>& rThis)
{
    rOStream << "DenseLinearSolverFactory" << std::endl;
    return rOStream;
}

///@}
///@name Input and output

void RegisterDenseLinearSolvers();

///@}

typedef LinearSolverFactory<LocalSparseSpaceType, LocalSparseSpaceType> DenseLinearSolverFactoryType;

#ifdef KRATOS_REGISTER_DENSE_LINEAR_SOLVER
#undef KRATOS_REGISTER_DENSE_LINEAR_SOLVER
#endif
#define KRATOS_REGISTER_DENSE_LINEAR_SOLVER(name, reference) ; \
    KratosComponents<DenseLinearSolverFactoryType>::Add(name, reference);

KRATOS_API_EXTERN template class KRATOS_API(EIGENSOLVERS_APPLICATION) KratosComponents<DenseLinearSolverFactoryType>;

typedef LinearSolverFactory<ComplexLocalSparseSpaceType, ComplexLocalSparseSpaceType> ComplexDenseLinearSolverFactoryType;

#ifdef KRATOS_REGISTER_COMPLEX_DENSE_LINEAR_SOLVER
#undef KRATOS_REGISTER_COMPLEX_DENSE_LINEAR_SOLVER
#endif
#define KRATOS_REGISTER_COMPLEX_DENSE_LINEAR_SOLVER(name, reference) ; \
    KratosComponents<ComplexDenseLinearSolverFactoryType>::Add(name, reference);

KRATOS_API_EXTERN template class KRATOS_API(EIGENSOLVERS_APPLICATION) KratosComponents<ComplexDenseLinearSolverFactoryType>;

}  // namespace Kratos.

#endif // KRATOS_DENSE_LINEAR_SOLVER_FACTORY_H_INCLUDED  defined
