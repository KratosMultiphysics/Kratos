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
#include "solving_strategies/schemes/residual_based_bdf_custom_scheme.h"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;

typedef Scheme<SparseSpaceType, LocalSpaceType> SchemeType;
typedef ResidualBasedBDFCustomScheme<SparseSpaceType,  LocalSpaceType> ResidualBasedBDFCustomSchemeType;

//NOTE: here we must create persisting objects for the strategies
static ResidualBasedBDFCustomSchemeType msResidualBasedBDFCustomScheme;

template<>
std::vector<Internals::RegisteredPrototypeBase<SchemeType>> ResidualBasedBDFCustomSchemeType::msPrototypes{
    Internals::RegisteredPrototype<ResidualBasedBDFCustomSchemeType, SchemeType>(ResidualBasedBDFCustomSchemeType::Name(), msResidualBasedBDFCustomScheme)};

///@}

} /* namespace Kratos.*/
