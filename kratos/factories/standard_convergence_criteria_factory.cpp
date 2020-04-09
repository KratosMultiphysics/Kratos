//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "factories/standard_convergence_criteria_factory.h"
#include "spaces/ublas_space.h"

// Convergence criterias
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "solving_strategies/convergencecriterias/displacement_criteria.h"
#include "solving_strategies/convergencecriterias/residual_criteria.h"
#include "solving_strategies/convergencecriterias/and_criteria.h"
#include "solving_strategies/convergencecriterias/or_criteria.h"

namespace Kratos
{
    void RegisterConvergenceCriterias()
    {
        typedef TUblasSparseSpace<double> SpaceType;
        typedef TUblasDenseSpace<double> LocalSpaceType;

        typedef ConvergenceCriteria<SpaceType,  LocalSpaceType> ConvergenceCriteriaType;
        typedef DisplacementCriteria<SpaceType,  LocalSpaceType> DisplacementCriteriaType;
        typedef ResidualCriteria<SpaceType,  LocalSpaceType> ResidualCriteriaType;
        typedef And_Criteria<SpaceType,  LocalSpaceType> And_CriteriaType;
        typedef Or_Criteria<SpaceType,  LocalSpaceType> Or_CriteriaType;

        //NOTE: here we must create persisting objects for the linear solvers
        static auto DisplacementCriteriaFactory = StandardConvergenceCriteriaFactory<ConvergenceCriteriaType, DisplacementCriteriaType>();
        static auto ResidualCriteriaFactory= StandardConvergenceCriteriaFactory<ConvergenceCriteriaType, ResidualCriteriaType>();
        static auto And_CriteriaFactory= StandardConvergenceCriteriaFactory<ConvergenceCriteriaType, And_CriteriaType>();
        static auto Or_CriteriaFactory= StandardConvergenceCriteriaFactory<ConvergenceCriteriaType, Or_CriteriaType>();

        // Registration of convergence solvers
        KRATOS_REGISTER_CONVERGENCE_CRITERIA(DisplacementCriteriaType::Name(), DisplacementCriteriaFactory);
        KRATOS_REGISTER_CONVERGENCE_CRITERIA(ResidualCriteriaType::Name(), ResidualCriteriaFactory);
        KRATOS_REGISTER_CONVERGENCE_CRITERIA(And_CriteriaType::Name(), And_CriteriaFactory );
        KRATOS_REGISTER_CONVERGENCE_CRITERIA(Or_CriteriaType::Name(), Or_CriteriaFactory );
    };
} // Namespace Kratos

