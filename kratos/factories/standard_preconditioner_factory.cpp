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
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "factories/standard_preconditioner_factory.h"
#include "spaces/ublas_space.h"
#include "spaces/default_spaces.h"

#include "linear_solvers/preconditioner.h"
#include "linear_solvers/diagonal_preconditioner.h"
#include "linear_solvers/ilu0_preconditioner.h"
#include "linear_solvers/ilu_preconditioner.h"

namespace Kratos
{
    namespace detail
    {
        template<class SpaceType, class LocalSpaceType>
        void RegisterPreconditionersForSpace()
        {
            typedef Preconditioner<SpaceType,  LocalSpaceType> PreconditionerType;
            typedef DiagonalPreconditioner<SpaceType,  LocalSpaceType> DiagonalPreconditionerType;
            typedef ILU0Preconditioner<SpaceType,  LocalSpaceType> ILU0PreconditionerType;
            typedef ILUPreconditioner<SpaceType,  LocalSpaceType> ILUPreconditionerType;

            // NOTE: declared before the factory objects because the local
            // variable below shadows the PreconditionerFactory class template
            using ComponentsType = KratosComponents<PreconditionerFactory<SpaceType, LocalSpaceType>>;

            //NOTE: here we must create persisting objects for the linear solvers
            static auto PreconditionerFactory = StandardPreconditionerFactory<SpaceType,LocalSpaceType,PreconditionerType>();
            static auto DiagonalPreconditionerFactory = StandardPreconditionerFactory<SpaceType,LocalSpaceType,DiagonalPreconditionerType>();
            static auto ILU0PreconditionerFactory= StandardPreconditionerFactory<SpaceType,LocalSpaceType,ILU0PreconditionerType>();
            static auto ILUPreconditionerFactory= StandardPreconditionerFactory<SpaceType,LocalSpaceType,ILUPreconditionerType>();

            //registration of linear solvers
            ComponentsType::Add("none", PreconditionerFactory);
            ComponentsType::Add("diagonal", DiagonalPreconditionerFactory);
            ComponentsType::Add("ilu0", ILU0PreconditionerFactory);
            ComponentsType::Add("ilu", ILUPreconditionerFactory);
        }
    } // namespace detail

    void RegisterPreconditioners()
    {
        detail::RegisterPreconditionersForSpace<TUblasSparseSpace<double>, TUblasDenseSpace<double>>();
#ifdef KRATOS_USE_EIGEN_BACKEND
        // The default (Eigen) sparse space additionally gets its own registrations
        detail::RegisterPreconditionersForSpace<TDefaultSparseSpace<double>, TDefaultDenseSpace<double>>();
#endif
    };
} // Namespace Kratos
