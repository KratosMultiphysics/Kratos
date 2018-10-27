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

#if !defined(KRATOS_CONVERGENCE_CRITERIA_FACTORY_H_INCLUDED )
#define  KRATOS_CONVERGENCE_CRITERIA_FACTORY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/kratos_components.h"
#include "includes/shared_pointers.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
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
 * @class ConvergenceCriteriaFactory
 * @ingroup KratosCore
 * @brief Here we add the functions needed for the registration of convergence criterias
 * @details Defines the base convergence criteria factory
 * @author Vicente Mataix Ferrandiz
 * @tparam TSparseSpace The sparse space definition
 * @tparam TLocalSpace The dense space definition
 */
template< typename TSparseSpace, typename TLocalSpace>
class ConvergenceCriteriaFactory
{
public:
    ///@name Type Definitions
    ///@{

    /// The definition of the custom class
    typedef ConvergenceCriteriaFactory<TSparseSpace,TLocalSpace> FactoryType;

    /// The definition of the convergence criteria
    typedef ConvergenceCriteria<TSparseSpace,TLocalSpace> ConvergenceCriteriaType;

    /// Pointer definition of ConvergenceCriteriaFactory
    KRATOS_CLASS_POINTER_DEFINITION(ConvergenceCriteriaFactory );

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
     * @brief This method checks if the linear convergence criteria is registered
     * @return True if registered, false otherwise
     */
    virtual bool Has(const std::string& rSolverType)
    {
        return KratosComponents< FactoryType >::Has( rSolverType );
    }

    /**
     * @brief This method creates a new convergence criteria
     * @return The pointer to the convergence criteria of interest
     */
    virtual typename ConvergenceCriteriaType::Pointer Create(Kratos::Parameters Settings)
    {
        const std::string& convergence_criterion = Settings["convergence_criterion"].GetString();
        if(Has( convergence_criterion ) == false) {
            KRATOS_ERROR << "Trying to construct a convergence criteria with type convergence_criterion= " << convergence_criterion << std::endl <<
                            "Which does not exist. The list of available options (for currently loaded applications) is: " << std::endl <<
                         KratosComponents< FactoryType >() << std::endl;
        }
        const auto& aux = KratosComponents< FactoryType >::Get( convergence_criterion );
        return aux.CreateConvergenceCriteria(Settings);
    }
    ///@}
protected:
    ///@name Protected Operators
    ///@{

    /**
     * @brief This method is an auxiliar method to create a new convergence criteria
     * @return The pointer to the convergence criteria of interest
     */
    virtual typename ConvergenceCriteriaType::Pointer CreateConvergenceCriteria(Kratos::Parameters Settings)  const
    {
        KRATOS_ERROR << "calling the base class ConvergenceCriteriaFactory" << std::endl;
    }
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
                                  const ConvergenceCriteriaFactory<TSparseSpace, TLocalSpace>& rThis)
{
    rOStream << "ConvergenceCriteriaFactory" << std::endl;

    return rOStream;
}

///@}

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSparseSpaceType;

typedef ConvergenceCriteriaFactory<SparseSpaceType, LocalSparseSpaceType> ConvergenceCriteriaFactoryType;

#ifdef KRATOS_REGISTER_CONVERGENCE_CRITERIA
#undef KRATOS_REGISTER_CONVERGENCE_CRITERIA
#endif
#define KRATOS_REGISTER_CONVERGENCE_CRITERIA(name, reference) \
    KratosComponents<ConvergenceCriteriaFactoryType>::Add(name, reference);

}  // namespace Kratos.

#endif // KRATOS_CONVERGENCE_CRITERIA_FACTORY_H_INCLUDED  defined
