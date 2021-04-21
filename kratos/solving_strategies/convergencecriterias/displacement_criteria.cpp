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
#include "solving_strategies/convergencecriterias/displacement_criteria.h"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;

typedef ConvergenceCriteria<SparseSpaceType, LocalSpaceType> ConvergenceCriteriaType;
typedef DisplacementCriteria<SparseSpaceType,  LocalSpaceType> DisplacementCriteriaType;

//NOTE: here we must create persisting objects for the strategies
static DisplacementCriteriaType msDisplacementCriteria;

template<>
std::vector<Internals::RegisteredPrototypeBase<ConvergenceCriteriaType>> DisplacementCriteriaType::msPrototypes{
    Internals::RegisteredPrototype<DisplacementCriteriaType, ConvergenceCriteriaType>(DisplacementCriteriaType::Name(), msDisplacementCriteria)};

///@}

} /* namespace Kratos.*/
