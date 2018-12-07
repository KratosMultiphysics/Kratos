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
#include "includes/standard_preconditioner_factory.h"
#include "spaces/ublas_space.h"

#include "linear_solvers/preconditioner.h"
#include "linear_solvers/diagonal_preconditioner.h"
#include "linear_solvers/ilu0_preconditioner.h"
#include "linear_solvers/ilu_preconditioner.h"

namespace Kratos
{
    void RegisterPreconditioners()
    {
        typedef TUblasSparseSpace<double> SpaceType;
        typedef TUblasDenseSpace<double> LocalSpaceType;

        typedef Preconditioner<SpaceType,  LocalSpaceType> PreconditionerType;
        typedef DiagonalPreconditioner<SpaceType,  LocalSpaceType> DiagonalPreconditionerType;
        typedef ILU0Preconditioner<SpaceType,  LocalSpaceType> ILU0PreconditionerType;
        typedef ILUPreconditioner<SpaceType,  LocalSpaceType> ILUPreconditionerType;

        //NOTE: here we must create persisting objects for the linear solvers
        static auto PreconditionerFactory = StandardPreconditionerFactory<SpaceType,LocalSpaceType,PreconditionerType>();
        static auto DiagonalPreconditionerFactory = StandardPreconditionerFactory<SpaceType,LocalSpaceType,DiagonalPreconditionerType>();
        static auto ILU0PreconditionerFactory= StandardPreconditionerFactory<SpaceType,LocalSpaceType,ILU0PreconditionerType>();
        static auto ILUPreconditionerFactory= StandardPreconditionerFactory<SpaceType,LocalSpaceType,ILUPreconditionerType>();

        //registration of linear solvers
        KRATOS_REGISTER_PRECONDITIONER("None", PreconditionerFactory);
        KRATOS_REGISTER_PRECONDITIONER("DiagonalPreconditioner", DiagonalPreconditionerFactory);
        KRATOS_REGISTER_PRECONDITIONER("ILU0Preconditioner", ILU0PreconditionerFactory);
        KRATOS_REGISTER_PRECONDITIONER("ILUPreconditioner",ILUPreconditionerFactory );
    };
} // Namespace Kratos

