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
#include "factories/standard_scheme_factory.h"
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
    void RegisterSchemes()
    {
        typedef TUblasSparseSpace<double> SparseSpaceType;
        typedef TUblasDenseSpace<double> LocalSpaceType;

        typedef Scheme< SparseSpaceType, LocalSpaceType > BaseSchemeType;
        typedef ResidualBasedIncrementalUpdateStaticScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedIncrementalUpdateStaticSchemeType;
        typedef ResidualBasedIncrementalUpdateStaticSchemeSlip< SparseSpaceType, LocalSpaceType >  ResidualBasedIncrementalUpdateStaticSchemeSlipType;
        typedef ResidualBasedBossakDisplacementScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedBossakDisplacementSchemeType;
        typedef ResidualBasedNewmarkDisplacementScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedNewmarkDisplacementSchemeType;
        typedef ResidualBasedPseudoStaticDisplacementScheme< SparseSpaceType, LocalSpaceType >  ResidualBasedPseudoStaticDisplacementSchemeType;
        typedef ResidualBasedBDFDisplacementScheme< SparseSpaceType, LocalSpaceType > ResidualBasedBDFDisplacementSchemeType;
        typedef ResidualBasedBDFCustomScheme< SparseSpaceType, LocalSpaceType > ResidualBasedBDFCustomSchemeType;

        //NOTE: here we must create persisting objects for the linear solvers
//         static auto BaseSchemeFactory = StandardSchemeFactory<BaseSchemeType, BaseSchemeType>();
        static auto ResidualBasedIncrementalUpdateStaticSchemeFactory = StandardSchemeFactory<BaseSchemeType, ResidualBasedIncrementalUpdateStaticSchemeType>();
        static auto ResidualBasedIncrementalUpdateStaticSchemeSlipFactory = StandardSchemeFactory<BaseSchemeType, ResidualBasedIncrementalUpdateStaticSchemeSlipType>();
        static auto ResidualBasedBossakDisplacementSchemeFactory = StandardSchemeFactory<BaseSchemeType, ResidualBasedBossakDisplacementSchemeType>();
        static auto ResidualBasedNewmarkDisplacementSchemeFactory= StandardSchemeFactory<BaseSchemeType, ResidualBasedNewmarkDisplacementSchemeType>();
        static auto ResidualBasedPseudoStaticDisplacementSchemeFactory= StandardSchemeFactory<BaseSchemeType, ResidualBasedPseudoStaticDisplacementSchemeType>();
        static auto ResidualBasedBDFDisplacementSchemeFactory= StandardSchemeFactory<BaseSchemeType, ResidualBasedBDFDisplacementSchemeType>();
        static auto ResidualBasedBDFCustomSchemeFactory= StandardSchemeFactory<BaseSchemeType, ResidualBasedBDFCustomSchemeType>();

        // Registration of convergence solvers
//         KRATOS_REGISTER_SCHEME("Scheme", BaseSchemeFactory);
        KRATOS_REGISTER_SCHEME("static", ResidualBasedIncrementalUpdateStaticSchemeFactory);
        KRATOS_REGISTER_SCHEME("ResidualBasedIncrementalUpdateStaticScheme", ResidualBasedIncrementalUpdateStaticSchemeFactory);
        KRATOS_REGISTER_SCHEME("static_slip", ResidualBasedIncrementalUpdateStaticSchemeSlipFactory);
        KRATOS_REGISTER_SCHEME("ResidualBasedIncrementalUpdateStaticSchemeSlip", ResidualBasedIncrementalUpdateStaticSchemeSlipFactory);
        KRATOS_REGISTER_SCHEME("ResidualBasedBossakDisplacementScheme", ResidualBasedBossakDisplacementSchemeFactory);
        KRATOS_REGISTER_SCHEME("bossak", ResidualBasedBossakDisplacementSchemeFactory);
        KRATOS_REGISTER_SCHEME("ResidualBasedNewmarkDisplacementScheme", ResidualBasedNewmarkDisplacementSchemeFactory);
        KRATOS_REGISTER_SCHEME("newmark", ResidualBasedNewmarkDisplacementSchemeFactory);
        KRATOS_REGISTER_SCHEME("ResidualBasedPseudoStaticDisplacementScheme", ResidualBasedPseudoStaticDisplacementSchemeFactory);
        KRATOS_REGISTER_SCHEME("pseudo_static", ResidualBasedPseudoStaticDisplacementSchemeFactory);
        KRATOS_REGISTER_SCHEME("ResidualBasedBDFDisplacementScheme", ResidualBasedBDFDisplacementSchemeFactory);
        KRATOS_REGISTER_SCHEME("bdf_displacement", ResidualBasedBDFDisplacementSchemeFactory);
        KRATOS_REGISTER_SCHEME("ResidualBasedBDFCustomScheme", ResidualBasedBDFCustomSchemeFactory);
        KRATOS_REGISTER_SCHEME("bdf", ResidualBasedBDFCustomSchemeFactory);
    };
} // Namespace Kratos

