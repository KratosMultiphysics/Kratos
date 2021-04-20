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
//
//

/* System includes */

/* External includes */

/* Project includes */
#include "spaces/ublas_space.h"
#include "solving_strategies/schemes/residual_based_newmark_displacement_scheme.hpp"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;

typedef Scheme<SparseSpaceType, LocalSpaceType> SchemeType;
typedef ResidualBasedNewmarkDisplacementScheme<SparseSpaceType,  LocalSpaceType> ResidualBasedNewmarkDisplacementSchemeType;

//NOTE: here we must create persisting objects for the strategies
static ResidualBasedNewmarkDisplacementSchemeType msResidualBasedNewmarkDisplacementScheme;

template<>
std::vector<Internals::RegisteredPrototypeBase<SchemeType>> ResidualBasedNewmarkDisplacementSchemeType::msPrototypes{
    Internals::RegisteredPrototype<ResidualBasedNewmarkDisplacementSchemeType, SchemeType>(ResidualBasedNewmarkDisplacementSchemeType::Name(), msResidualBasedNewmarkDisplacementScheme)};

///@}

} /* namespace Kratos.*/
