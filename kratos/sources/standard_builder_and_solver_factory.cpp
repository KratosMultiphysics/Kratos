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
#include "includes/standard_convergence_criteria_factory.h"
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

//         typedef ConvergenceCriteria<SpaceType,  LocalSpaceType> ConvergenceCriteriaType;
        typedef DisplacementCriteria<SpaceType,  LocalSpaceType> DisplacementCriteriaType;
        typedef ResidualCriteria<SpaceType,  LocalSpaceType> ResidualCriteriaType;
        typedef And_Criteria<SpaceType,  LocalSpaceType> And_CriteriaType;
        typedef Or_Criteria<SpaceType,  LocalSpaceType> Or_CriteriaType;

        //NOTE: here we must create persisting objects for the linear solvers
//         static auto ConvergenceCriteriaFactory = StandardConvergenceCriteriaFactory<SpaceType,LocalSpaceType,ConvergenceCriteriaType>();
        static auto DisplacementCriteriaFactory = StandardConvergenceCriteriaFactory<SpaceType,LocalSpaceType,DisplacementCriteriaType>();
        static auto ResidualCriteriaFactory= StandardConvergenceCriteriaFactory<SpaceType,LocalSpaceType,ResidualCriteriaType>();
        static auto And_CriteriaFactory= StandardConvergenceCriteriaFactory<SpaceType,LocalSpaceType,And_CriteriaType>();
        static auto Or_CriteriaFactory= StandardConvergenceCriteriaFactory<SpaceType,LocalSpaceType,Or_CriteriaType>();

        // Registration of convergence solvers
//         KRATOS_REGISTER_CONVERGENCE_CRITERIA("ConvergenceCriteria", ConvergenceCriteriaFactory);
        KRATOS_REGISTER_CONVERGENCE_CRITERIA("DisplacementCriteria", DisplacementCriteriaFactory);
        KRATOS_REGISTER_CONVERGENCE_CRITERIA("displacement_criteria", DisplacementCriteriaFactory);
        KRATOS_REGISTER_CONVERGENCE_CRITERIA("ResidualCriteria", ResidualCriteriaFactory);
        KRATOS_REGISTER_CONVERGENCE_CRITERIA("residual_criteria", ResidualCriteriaFactory);
        KRATOS_REGISTER_CONVERGENCE_CRITERIA("And_Criteria",And_CriteriaFactory );
        KRATOS_REGISTER_CONVERGENCE_CRITERIA("and_criteria",And_CriteriaFactory );
        KRATOS_REGISTER_CONVERGENCE_CRITERIA("Or_Criteria",Or_CriteriaFactory );
        KRATOS_REGISTER_CONVERGENCE_CRITERIA("or_criteria",Or_CriteriaFactory );
    };
} // Namespace Kratos

