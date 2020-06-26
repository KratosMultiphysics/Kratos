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
//                   Vicente Mataix Ferrandiz
//
//

#if !defined(KRATOS_STANDARD_PRECONDITIONER_FACTORY_H_INCLUDED )
#define  KRATOS_STANDARD_PRECONDITIONER_FACTORY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "factories/preconditioner_factory.h"
#include "linear_solvers/preconditioner.h"

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
 * @class StandardPreconditionerFactory
 * @ingroup KratosCore
 * @brief Here we add the functions needed for the registration of preconditioners
 * @details Defines the standard preconditioner factory
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 * @tparam TSparseSpace The sparse space definition
 * @tparam TLocalSpace The dense space definition
 * @tparam TPreconditionerType The precondioner type
 */
template <typename TSparseSpace, typename TLocalSpace, typename TPreconditionerType>
class StandardPreconditionerFactory
    : public PreconditionerFactory<TSparseSpace,TLocalSpace>
{
    ///@name Type Definitions
    ///@{

    /// The definition of the preconditioner
    typedef Preconditioner<TSparseSpace,TLocalSpace> PreconditionerType;

    ///@}
protected:

    ///@name Operations
    ///@{

    /**
     * @brief This method is an auxiliar method to create a new solver
     * @return The pointer to the solver of interest
     */
    typename PreconditionerType::Pointer CreatePreconditioner() const override
    {
        return typename PreconditionerType::Pointer(new TPreconditionerType());
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
template <typename TSparseSpace, typename TLocalSpace, typename TPreconditionerType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const StandardPreconditionerFactory<TSparseSpace,TLocalSpace,TPreconditionerType>& rThis)
{
    rOStream << "StandardPreconditionerFactory" << std::endl;
    return rOStream;
}

///@}
///@name Input and output

void RegisterPreconditioners();

///@}

}  // namespace Kratos.

#endif // KRATOS_STANDARD_PRECONDITIONER_FACTORY_H_INCLUDED  defined
