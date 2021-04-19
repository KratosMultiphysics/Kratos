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
#include "solving_strategies/schemes/residual_based_bdf_displacement_scheme.h"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;

typedef Scheme<SparseSpaceType, LocalSpaceType> SchemeType;
typedef ResidualBasedBDFDisplacementScheme<SparseSpaceType,  LocalSpaceType> ResidualBasedBDFDisplacementSchemeType;

//NOTE: here we must create persisting objects for the strategies
static ResidualBasedBDFDisplacementSchemeType msResidualBasedBDFDisplacementScheme;

template<>
std::vector<Internals::RegisteredPrototypeBase<SchemeType>> ResidualBasedBDFDisplacementSchemeType::msPrototypes{
    Internals::RegisteredPrototype<ResidualBasedBDFDisplacementSchemeType, SchemeType>(ResidualBasedBDFDisplacementSchemeType::Name(), msResidualBasedBDFDisplacementScheme)};

///@}

} /* namespace Kratos.*/
