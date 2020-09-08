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
#include "factories/factory.h"
#include "spaces/ublas_space.h"

// Schemes
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip.h"
#include "solving_strategies/schemes/residual_based_bossak_displacement_scheme.hpp"
#include "solving_strategies/schemes/residual_based_newmark_displacement_scheme.hpp"
#include "solving_strategies/schemes/residual_based_pseudo_static_displacement_scheme.h"
#include "solving_strategies/schemes/residual_based_bdf_displacement_scheme.h"
#include "solving_strategies/schemes/residual_based_bdf_custom_scheme.h"

namespace Kratos
{
void RegisterSchemesFactories()
{
    typedef TUblasSparseSpace<double> SparseSpaceType;
    typedef TUblasDenseSpace<double> LocalSpaceType;

    typedef ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedIncrementalUpdateStaticSchemeType;
    typedef ResidualBasedIncrementalUpdateStaticSchemeSlip< SparseSpaceType, LocalSpaceType >  ResidualBasedIncrementalUpdateStaticSchemeSlipType;
    typedef ResidualBasedBossakDisplacementScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedBossakDisplacementSchemeType;
    typedef ResidualBasedNewmarkDisplacementScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedNewmarkDisplacementSchemeType;
    typedef ResidualBasedPseudoStaticDisplacementScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedPseudoStaticDisplacementSchemeType;
    typedef ResidualBasedBDFDisplacementScheme< SparseSpaceType, LocalSpaceType > ResidualBasedBDFDisplacementSchemeType;
    typedef ResidualBasedBDFCustomScheme< SparseSpaceType, LocalSpaceType > ResidualBasedBDFCustomSchemeType;

    //NOTE: here we must create persisting objects for the linear solvers
    static ResidualBasedIncrementalUpdateStaticSchemeType msResidualBasedIncrementalUpdateStaticScheme;
    static ResidualBasedIncrementalUpdateStaticSchemeSlipType msResidualBasedIncrementalUpdateStaticSchemeSlip;
    static ResidualBasedBossakDisplacementSchemeType msResidualBasedBossakDisplacementScheme;
    static ResidualBasedNewmarkDisplacementSchemeType msResidualBasedNewmarkDisplacementScheme;
    static ResidualBasedPseudoStaticDisplacementSchemeType msResidualBasedPseudoStaticDisplacementScheme;
    static ResidualBasedBDFDisplacementSchemeType msResidualBasedBDFDisplacementScheme;
    static ResidualBasedBDFCustomSchemeType msResidualBasedBDFCustomScheme;

    // Registration of schemes
    KRATOS_REGISTER_SCHEME(ResidualBasedIncrementalUpdateStaticSchemeType::Name(), msResidualBasedIncrementalUpdateStaticScheme);
    KRATOS_REGISTER_SCHEME(ResidualBasedIncrementalUpdateStaticSchemeSlipType::Name(), msResidualBasedIncrementalUpdateStaticSchemeSlip);
    KRATOS_REGISTER_SCHEME(ResidualBasedBossakDisplacementSchemeType::Name(), msResidualBasedBossakDisplacementScheme);
    KRATOS_REGISTER_SCHEME(ResidualBasedNewmarkDisplacementSchemeType::Name(), msResidualBasedNewmarkDisplacementScheme);
    KRATOS_REGISTER_SCHEME(ResidualBasedPseudoStaticDisplacementSchemeType::Name(), msResidualBasedPseudoStaticDisplacementScheme);
    KRATOS_REGISTER_SCHEME(ResidualBasedBDFDisplacementSchemeType::Name(), msResidualBasedBDFDisplacementScheme);
    KRATOS_REGISTER_SCHEME(ResidualBasedBDFCustomSchemeType::Name(), msResidualBasedBDFCustomScheme);
};
} // Namespace Kratos

