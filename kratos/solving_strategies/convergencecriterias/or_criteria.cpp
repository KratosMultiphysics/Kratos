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
#include "solving_strategies/convergencecriterias/or_criteria.h"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;

typedef ConvergenceCriteria<SparseSpaceType, LocalSpaceType> ConvergenceCriteriaType;
typedef Or_Criteria<SparseSpaceType,  LocalSpaceType> Or_CriteriaType;

//NOTE: here we must create persisting objects for the strategies
static Or_CriteriaType msOr_Criteria;

template<>
std::vector<Internals::RegisteredPrototypeBase<ConvergenceCriteriaType>> Or_CriteriaType::msPrototypes{
    Internals::RegisteredPrototype<Or_CriteriaType, ConvergenceCriteriaType>(Or_CriteriaType::Name(), msOr_Criteria)};

///@}

} /* namespace Kratos.*/
