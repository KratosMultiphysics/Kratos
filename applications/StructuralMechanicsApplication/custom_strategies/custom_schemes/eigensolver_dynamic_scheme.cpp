// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//   Project Name:        $StructuralMechanicsApplication $
//   Last modified by:    $Author: michael.andre@tum.de   $
//   Date:                $Date:         September 2016   $
//   Revision:            $Revision:                0.0   $
//

/* System includes */

/* External includes */

/* Project includes */
#include "spaces/ublas_space.h"
#include "custom_strategies/custom_schemes/eigensolver_dynamic_scheme.hpp"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;

typedef Scheme<SparseSpaceType, LocalSpaceType> SchemeType;
typedef EigensolverDynamicScheme<SparseSpaceType,  LocalSpaceType> EigensolverDynamicSchemeType;

//NOTE: here we must create persisting objects for the strategies
static EigensolverDynamicSchemeType msEigensolverDynamicScheme;

template<>
std::vector<Internals::RegisteredPrototypeBase<SchemeType>> EigensolverDynamicSchemeType::msPrototypes{
    Internals::RegisteredPrototype<EigensolverDynamicSchemeType, SchemeType>(EigensolverDynamicSchemeType::Name(), msEigensolverDynamicScheme)};

///@}

} /* namespace Kratos.*/