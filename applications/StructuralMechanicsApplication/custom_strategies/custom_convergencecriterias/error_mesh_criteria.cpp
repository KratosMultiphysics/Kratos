// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Anna Rehr
//                   Vicente Mataix Ferrandiz
//

/* System includes */

/* External includes */

/* Project includes */
#include "spaces/ublas_space.h"
#include "custom_strategies/custom_convergencecriterias/error_mesh_criteria.h"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;

typedef ConvergenceCriteria<SparseSpaceType, LocalSpaceType> ConvergenceCriteriaType;
typedef ErrorMeshCriteria<SparseSpaceType,  LocalSpaceType> ErrorMeshCriteriaType;

//NOTE: here we must create persisting objects for the strategies
static ErrorMeshCriteriaType msErrorMeshCriteria;

template<>
std::vector<Internals::RegisteredPrototypeBase<ConvergenceCriteriaType>> ErrorMeshCriteriaType::msPrototypes{
    Internals::RegisteredPrototype<ErrorMeshCriteriaType, ConvergenceCriteriaType>(ErrorMeshCriteriaType::Name(), msErrorMeshCriteria)};

///@}

} /* namespace Kratos.*/
