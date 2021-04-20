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
#include "solving_strategies/convergencecriterias/residual_criteria.h"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;

typedef ConvergenceCriteria<SparseSpaceType, LocalSpaceType> ConvergenceCriteriaType;
typedef ResidualCriteria<SparseSpaceType,  LocalSpaceType> ResidualCriteriaType;

//NOTE: here we must create persisting objects for the strategies
static ResidualCriteriaType msResidualCriteria;

template<>
std::vector<Internals::RegisteredPrototypeBase<ConvergenceCriteriaType>> ResidualCriteriaType::msPrototypes{
    Internals::RegisteredPrototype<ResidualCriteriaType, ConvergenceCriteriaType>(ResidualCriteriaType::Name(), msResidualCriteria)};

///@}

} /* namespace Kratos.*/
