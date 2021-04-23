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
#include "solving_strategies/convergencecriterias/and_criteria.h"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;

typedef ConvergenceCriteria<SparseSpaceType, LocalSpaceType> ConvergenceCriteriaType;
typedef And_Criteria<SparseSpaceType,  LocalSpaceType> And_CriteriaType;

//NOTE: here we must create persisting objects for the strategies
static And_CriteriaType msAnd_Criteria;

template<>
std::vector<Internals::RegisteredPrototypeBase<ConvergenceCriteriaType>> And_CriteriaType::msPrototypes{
    Internals::RegisteredPrototype<And_CriteriaType, ConvergenceCriteriaType>(And_CriteriaType::Name(), msAnd_Criteria)};

///@}

} /* namespace Kratos.*/
