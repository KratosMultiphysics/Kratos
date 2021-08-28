// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Klaus B Sautter (based on the work of MSantasusana)
//
//

/* System includes */

/* External includes */

/* Project includes */
#include "spaces/ublas_space.h"
#include "custom_strategies/custom_schemes/explicit_central_differences_scheme.hpp"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;

typedef Scheme<SparseSpaceType, LocalSpaceType> SchemeType;
typedef ExplicitCentralDifferencesScheme<SparseSpaceType,  LocalSpaceType> ExplicitCentralDifferencesSchemeType;

//NOTE: here we must create persisting objects for the strategies
static ExplicitCentralDifferencesSchemeType msExplicitCentralDifferencesScheme;

template<>
std::vector<Internals::RegisteredPrototypeBase<SchemeType>> ExplicitCentralDifferencesSchemeType::msPrototypes{
    Internals::RegisteredPrototype<ExplicitCentralDifferencesSchemeType, SchemeType>(ExplicitCentralDifferencesSchemeType::Name(), msExplicitCentralDifferencesScheme)};

///@}

} /* namespace Kratos.*/