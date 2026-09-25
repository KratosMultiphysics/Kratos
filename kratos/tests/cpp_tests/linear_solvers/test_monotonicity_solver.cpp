//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "includes/kratos_filesystem.h"
#include "includes/matrix_market_interface.h"

/* Utility includes */
#include "includes/define.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "spaces/default_spaces.h"
#include "utilities/variable_utils.h"

/* Linear Solver includes */
#include "linear_solvers/monotonicity_preserving_solver.h"

namespace Kratos::Testing
{

using SpaceType = TDefaultSparseSpace<double>;
using LocalSpaceType = TDefaultDenseSpace<double>;
using SparseMatrixType = typename SpaceType::MatrixType;

KRATOS_TEST_CASE_IN_SUITE(MonotonictyPreservingSolver, KratosCoreFastSuite)
{
    Model current_model;
    ModelPart& model_part = current_model.CreateModelPart("Main");

    model_part.SetBufferSize(3);

    model_part.AddNodalSolutionStepVariable(TEMPERATURE);

    for (std::size_t i = 0; i < 5; i++) {
        model_part.CreateNewNode(i+1, 0.0, 0.0, 0.0);
    }

    VariableUtils().AddDof(TEMPERATURE, model_part);

    for (std::size_t i = 0; i < model_part.NumberOfNodes(); i++) {
        auto it_node = model_part.NodesBegin() + i;
        it_node->GetSolutionStepValue(TEMPERATURE) = 20.0*i;
    }

    Parameters settings( R"(
    {
        "solver_type"                    : "monotonicity_preserving",
        "inner_solver_settings"          : {
            "preconditioner_type"            : "amg",
            "solver_type"                    : "amgcl",
            "smoother_type"                  : "ilu0",
            "krylov_type"                    : "lgmres",
            "coarsening_type"                : "aggregation",
            "max_iteration"                  : 100,
            "provide_coordinates"            : false,
            "gmres_krylov_space_dimension"   : 100,
            "verbosity"                      : 1,
            "tolerance"                      : 1e-6,
            "scaling"                        : false,
            "block_size"                     : 1,
            "use_block_matrices_if_possible" : true,
            "coarse_enough"                  : 1000,
            "max_levels"                     : -1,
            "pre_sweeps"                     : 1,
            "post_sweeps"                    : 1,
            "use_gpgpu"                      : false
        }

    }  )" );

    MonotonicityPreservingSolver<SpaceType, LocalSpaceType> linear_solver(settings);
    SpaceType::VectorType rX(5);
    SpaceType::SetToZero(rX);
    SpaceType::VectorType rB(5);
    SpaceType::SetToZero(rB);
    ModelPart::DofsArrayType dof_set;
    for (std::size_t i = 0; i < model_part.NumberOfNodes(); i++) {
        auto it_node = model_part.NodesBegin() + i;
        auto p_dof = it_node->pGetDof(TEMPERATURE);
        p_dof->SetEquationId(i);
        dof_set.push_back(p_dof);
    }
    // Build the matrix by writing the CSR arrays directly (row-sorted entries),
    // valid for both the uBLAS and the Eigen backend matrix; the solver reads
    // the compressed storage through index1_data()/index2_data()/value_data().
    // row 0: (0,0)=1 (0,1)=1 (0,2)=2 (0,4)=-1 | row 1: (1,0)=1 (1,1)=1
    // row 2: (2,0)=2 (2,2)=1 | row 3: (3,3)=1 (3,4)=3 | row 4: (4,0)=-1 (4,3)=3 (4,4)=1
    const std::vector<std::vector<std::pair<std::size_t, double>>> entries {
        {{0, 1.0}, {1, 1.0}, {2, 2.0}, {4, -1.0}},
        {{0, 1.0}, {1, 1.0}},
        {{0, 2.0}, {2, 1.0}},
        {{3, 1.0}, {4, 3.0}},
        {{0, -1.0}, {3, 3.0}, {4, 1.0}}
    };
    std::size_t nnz = 0;
    for (const auto& r_row : entries) nnz += r_row.size();
    SpaceType::MatrixType rA(5, 5, nnz);
    {
        auto row_ptr = rA.index1_data().begin();
        auto col_ptr = rA.index2_data().begin();
        auto val_ptr = rA.value_data().begin();
        row_ptr[0] = 0;
        std::size_t k = 0;
        for (std::size_t i = 0; i < entries.size(); ++i) {
            for (const auto& r_entry : entries[i]) {
                col_ptr[k] = r_entry.first;
                val_ptr[k] = r_entry.second;
                ++k;
            }
            row_ptr[i + 1] = k;
        }
        rA.set_filled(entries.size() + 1, nnz);
    }
    if (linear_solver.AdditionalPhysicalDataIsNeeded()) {
        linear_solver.ProvideAdditionalData(rA, rX, rB, dof_set, model_part);
    }

    SpaceType::VectorType reference_b(5);
    reference_b[0] = 100.0;
    reference_b[1] = -20.0;
    reference_b[2] = -80.0;
    reference_b[3] = 60.0;
    reference_b[4] = -60.0;

    KRATOS_EXPECT_NEAR(rA(1,0), 0.0, 1e-10);
    KRATOS_EXPECT_NEAR(rA(0,1), 0.0, 1e-10);
    KRATOS_EXPECT_NEAR(rA(2,0), 0.0, 1e-10);
    KRATOS_EXPECT_NEAR(rA(0,2), 0.0, 1e-10);
    KRATOS_EXPECT_NEAR(rA(3,4), 0.0, 1e-10);
    KRATOS_EXPECT_NEAR(rA(4,3), 0.0, 1e-10);
    KRATOS_EXPECT_NEAR(rA(4,0), -1.0, 1e-10);
    KRATOS_EXPECT_NEAR(rA(0,4), -1.0, 1e-10);
    KRATOS_EXPECT_NEAR(rA(0,0), 4.0, 1e-10);
    KRATOS_EXPECT_NEAR(rA(1,1), 2.0, 1e-10);
    KRATOS_EXPECT_NEAR(rA(2,2), 3.0, 1e-10);
    KRATOS_EXPECT_NEAR(rA(3,3), 4.0, 1e-10);
    KRATOS_EXPECT_NEAR(rA(4,4), 4.0, 1e-10);
    KRATOS_EXPECT_VECTOR_NEAR(rB, reference_b, 1e-10);

}
}  // namespace Kratos::Testing.

