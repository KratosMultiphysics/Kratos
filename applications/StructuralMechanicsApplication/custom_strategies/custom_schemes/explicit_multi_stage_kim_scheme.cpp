// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Klaus B Sautter
//  Based on : "An accurate two‐stage explicit time integration scheme for structural dynamics and various dynamic problems" - Wooram Kim
//

/* System includes */

/* External includes */

/* Project includes */
#include "spaces/ublas_space.h"
#include "custom_strategies/custom_schemes/explicit_multi_stage_kim_scheme.hpp"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;

typedef Scheme<SparseSpaceType, LocalSpaceType> SchemeType;
typedef ExplicitMultiStageKimScheme<SparseSpaceType,  LocalSpaceType> ExplicitMultiStageKimSchemeType;

//NOTE: here we must create persisting objects for the strategies
static ExplicitMultiStageKimSchemeType msExplicitMultiStageKimScheme;

template<>
std::vector<Internals::RegisteredPrototypeBase<SchemeType>> ExplicitMultiStageKimSchemeType::msPrototypes{
    Internals::RegisteredPrototype<ExplicitMultiStageKimSchemeType, SchemeType>(ExplicitMultiStageKimSchemeType::Name(), msExplicitMultiStageKimScheme)};

///@}

} /* namespace Kratos.*/