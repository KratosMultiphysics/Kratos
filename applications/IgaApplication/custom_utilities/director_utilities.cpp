//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "director_utilities.h"
#include "factories/linear_solver_factory.h"


namespace Kratos
{
    DirectorUtilities::DirectorUtilities(
        ModelPart& rModelPart,
        Parameters JsonParameters)
        : mrModelPart(rModelPart)
        , mParameters(JsonParameters)
    {
    }

    void DirectorUtilities::ComputeDirectors()
    {
        Vector brep_ids = mParameters["brep_ids"].GetVector();
        for (IndexType i = 0; i < brep_ids.size(); ++i) {
            auto p_geom = mrModelPart.GetRootModelPart().pGetGeometry((IndexType)brep_ids[i]);
            KRATOS_WATCH(p_geom->Id())
        }

        SizeType system_dofs = 3;
        SparseMatrixType matrix(3, 3);

        Parameters solver_parameters(mParameters["linear_solver_settings"]);
        if (!solver_parameters.Has("solver_type")) solver_parameters.AddString("solver_type", "skyline_lu_factorization");

        DenseVectorType solution(system_dofs);
        DenseVectorType projector_transpose_column(system_dofs);
        auto solver = LinearSolverFactory<SparseSpaceType, LocalSpaceType>().Create(solver_parameters);

        solver->Solve(matrix, solution, projector_transpose_column);
        KRATOS_WATCH(solution)
    }
}  // namespace Kratos.

#endif // KRATOS_IGA_FLAGS_CPP_INCLUDED  defined