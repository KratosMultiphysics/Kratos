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
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip.h"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;

typedef Scheme<SparseSpaceType, LocalSpaceType> SchemeType;
typedef ResidualBasedIncrementalUpdateStaticSchemeSlip<SparseSpaceType,  LocalSpaceType> ResidualBasedIncrementalUpdateStaticSchemeSlipType;

//NOTE: here we must create persisting objects for the strategies
static ResidualBasedIncrementalUpdateStaticSchemeSlipType msResidualBasedIncrementalUpdateStaticSchemeSlip;

template<>
std::vector<Internals::RegisteredPrototypeBase<SchemeType>> ResidualBasedIncrementalUpdateStaticSchemeSlipType::msPrototypes{
    Internals::RegisteredPrototype<ResidualBasedIncrementalUpdateStaticSchemeSlipType, SchemeType>(ResidualBasedIncrementalUpdateStaticSchemeSlipType::Name(), msResidualBasedIncrementalUpdateStaticSchemeSlip)};

///@}

} /* namespace Kratos.*/
