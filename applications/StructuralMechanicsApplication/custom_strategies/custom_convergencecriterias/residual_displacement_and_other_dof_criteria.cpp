// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

/* System includes */

/* External includes */

/* Project includes */
#include "spaces/ublas_space.h"
#include "custom_strategies/custom_convergencecriterias/residual_displacement_and_other_dof_criteria.h"

namespace Kratos
{

///@name Type Definitions
///@{

typedef TUblasSparseSpace<double> SparseSpaceType;
typedef TUblasDenseSpace<double> LocalSpaceType;

typedef ConvergenceCriteria<SparseSpaceType, LocalSpaceType> ConvergenceCriteriaType;
typedef ResidualDisplacementAndOtherDoFCriteria<SparseSpaceType,  LocalSpaceType> ResidualDisplacementAndOtherDoFCriteriaType;

//NOTE: here we must create persisting objects for the strategies
static ResidualDisplacementAndOtherDoFCriteriaType msResidualDisplacementAndOtherDoFCriteria;

template<>
std::vector<Internals::RegisteredPrototypeBase<ConvergenceCriteriaType>> ResidualDisplacementAndOtherDoFCriteriaType::msPrototypes{
    Internals::RegisteredPrototype<ResidualDisplacementAndOtherDoFCriteriaType, ConvergenceCriteriaType>(ResidualDisplacementAndOtherDoFCriteriaType::Name(), msResidualDisplacementAndOtherDoFCriteria)};

///@}

} /* namespace Kratos.*/
