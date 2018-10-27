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

#if !defined(KRATOS_SCHEME_FACTORY_H_INCLUDED )
#define  KRATOS_SCHEME_FACTORY_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/kratos_components.h"
#include "includes/shared_pointers.h"
#include "solving_strategies/schemes/scheme.h"
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
 * @class SchemeFactory
 * @ingroup KratosCore
 * @brief Here we add the functions needed for the registration of schemes
 * @details Defines the base scheme factory
 * @author Vicente Mataix Ferrandiz
 * @tparam TSparseSpace The sparse space definition
 * @tparam TLocalSpace The dense space definition
 */
template< typename TSparseSpace, typename TLocalSpace>
class SchemeFactory
{
public:
    ///@name Type Definitions
    ///@{

    /// The definition of the custom class
    typedef SchemeFactory<TSparseSpace,TLocalSpace> FactoryType;

    /// The definition of the scheme
    typedef Scheme<TSparseSpace,TLocalSpace> SchemeType;

    /// Pointer definition of SchemeFactory
    KRATOS_CLASS_POINTER_DEFINITION(SchemeFactory );

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
     * @brief This method checks if the linear scheme is registered
     * @return True if registered, false otherwise
     */
    virtual bool Has(const std::string& rSolverType)
    {
        return KratosComponents< FactoryType >::Has( rSolverType );
    }

    /**
     * @brief This method creates a new scheme
     * @return The pointer to the scheme of interest
     */
    virtual typename SchemeType::Pointer Create(Kratos::Parameters Settings)
    {
        const std::string& scheme_type = Settings["scheme_type"].GetString();
        if(Has( scheme_type ) == false) {
            KRATOS_ERROR << "Trying to construct a scheme with type scheme_type= " << scheme_type << std::endl <<
                            "Which does not exist. The list of available options (for currently loaded applications) is: " << std::endl <<
                         KratosComponents< FactoryType >() << std::endl;
        }
        const auto& aux = KratosComponents< FactoryType >::Get( scheme_type );
        return aux.CreateScheme(Settings);
    }
    ///@}
protected:
    ///@name Protected Operators
    ///@{

    /**
     * @brief This method is an auxiliar method to create a new scheme
     * @return The pointer to the scheme of interest
     */
    virtual typename SchemeType::Pointer CreateScheme(Kratos::Parameters Settings)  const
    {
        KRATOS_ERROR << "calling the base class SchemeFactory" << std::endl;
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
                                  const SchemeFactory<TSparseSpace, TLocalSpace>& rThis)
{
    rOStream << "SchemeFactory" << std::endl;

    return rOStream;
}

///@}

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSparseSpaceType;

typedef SchemeFactory<SparseSpaceType, LocalSparseSpaceType> SchemeFactoryType;

#ifdef KRATOS_REGISTER_SCHEME
#undef KRATOS_REGISTER_SCHEME
#endif
#define KRATOS_REGISTER_SCHEME(name, reference) \
    KratosComponents<SchemeFactoryType>::Add(name, reference);

}  // namespace Kratos.

#endif // KRATOS_SCHEME_FACTORY_H_INCLUDED  defined
