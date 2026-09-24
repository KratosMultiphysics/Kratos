//  KRATOS  _____     _ _ _
//         |_   _| __(_) (_)_ __   ___  ___
//           | || '__| | | | '_ \ / _ \/ __|
//           | || |  | | | | | | | (_) \__
//           |_||_|  |_|_|_|_| |_|\___/|___/ APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "trilinos_space_experimental.h"
#include "custom_utilities/trilinos_dof_updater_experimental.h"
#include "tests/cpp_tests/trilinos_fast_suite.h"
#include "containers/model.h"
#include "includes/variables.h"
#include "mpi/includes/mpi_data_communicator.h"

namespace Kratos::Testing
{

/// Basic definitions
using TrilinosSparseSpaceType = TrilinosSpaceExperimental<Tpetra::FECrsMatrix<>, Tpetra::FEMultiVector<>>;

namespace
{

/**
 * @brief Builds a distributed update vector and a DOF set including a ghost DOF owned by the neighbour rank.
 * @details Each rank owns LocalSize consecutive equation ids with dx[id] = id. The DOF set additionally
 * contains (when running on more than one rank) the first DOF of the next rank, so that UpdateDofs must
 * import out-of-process data.
 */
struct DofUpdaterTestSetup
{
    using GO = TrilinosSparseSpaceType::GO;
    using LO = TrilinosSparseSpaceType::LO;

    TrilinosSparseSpaceType::VectorPointerType pDx;
    TrilinosSparseSpaceType::DofsArrayType DofSet;

    DofUpdaterTestSetup(
        Model& rModel,
        const DataCommunicator& rComm,
        const int LocalSize
        )
    {
        auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(rComm);
        TrilinosSparseSpaceType::CommunicatorPointerType tpetra_comm =
            Teuchos::rcp(new TrilinosSparseSpaceType::CommunicatorType(raw_mpi_comm));

        const int rank = rComm.Rank();
        const int world_size = rComm.Size();
        const int first_local_id = rank * LocalSize;

        // Build a distributed vector whose map covers [first_local_id, first_local_id + LocalSize)
        std::vector<GO> local_gids(LocalSize);
        for (int i = 0; i < LocalSize; ++i) local_gids[i] = static_cast<GO>(first_local_id + i);
        TrilinosSparseSpaceType::MapPointerType map = Teuchos::rcp(new TrilinosSparseSpaceType::MapType(
            Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(),
            Teuchos::ArrayView<const GO>(local_gids.data(), LocalSize), 0, tpetra_comm));
        pDx = TrilinosSparseSpaceType::CreateVector(map);
        auto& r_dx = *pDx;
        for (int i = 0; i < LocalSize; ++i) {
            const GO gid = local_gids[i];
            const LO lid = map->getLocalElement(gid);
            r_dx.replaceLocalValue(lid, size_t(0), static_cast<TrilinosSparseSpaceType::ST>(gid));
        }
        TrilinosSparseSpaceType::GlobalAssemble(r_dx);

        // Build the DOF set: the locally owned DOFs plus one ghost DOF owned by the next rank
        ModelPart& r_model_part = rModel.CreateModelPart("TestModelPart");
        r_model_part.AddNodalSolutionStepVariable(PRESSURE);
        r_model_part.GetNodalSolutionStepVariablesList().AddDof(&PRESSURE);

        std::vector<int> equation_ids;
        for (int i = 0; i < LocalSize; ++i) {
            equation_ids.push_back(first_local_id + i);
        }
        if (world_size > 1) {
            const int ghost_id = ((rank + 1) % world_size) * LocalSize;
            equation_ids.push_back(ghost_id);
        }

        DofSet.reserve(equation_ids.size());
        for (const int equation_id : equation_ids) {
            auto p_node = r_model_part.CreateNewNode(equation_id + 1, 0.0, 0.0, 0.0);
            p_node->AddDof(PRESSURE);
            auto p_dof = p_node->pGetDof(PRESSURE);
            p_dof->SetEquationId(equation_id);
            DofSet.push_back(p_dof);
        }
        DofSet.Sort();
    }
};

} // namespace

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalDofUpdaterUpdateDofs, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    Model model;
    const int local_size = 2;
    DofUpdaterTestSetup setup(model, r_comm, local_size);

    // The updater created by the space must be the Tpetra-aware one
    auto p_dof_updater = TrilinosSparseSpaceType::CreateDofUpdater();
    KRATOS_EXPECT_NE(nullptr, p_dof_updater.get());
    KRATOS_EXPECT_EQ(p_dof_updater->Info(), "TrilinosDofUpdaterExperimental");

    // Update the DOFs: every DOF (owned and ghost) must receive dx[equation_id] = equation_id
    p_dof_updater->UpdateDofs(setup.DofSet, *setup.pDx);
    for (const auto& r_dof : setup.DofSet) {
        KRATOS_EXPECT_DOUBLE_EQ(r_dof.GetSolutionStepValue(), static_cast<double>(r_dof.EquationId()));
    }

    // A second update accumulates (value += dx)
    p_dof_updater->UpdateDofs(setup.DofSet, *setup.pDx);
    for (const auto& r_dof : setup.DofSet) {
        KRATOS_EXPECT_DOUBLE_EQ(r_dof.GetSolutionStepValue(), 2.0 * static_cast<double>(r_dof.EquationId()));
    }
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosExperimentalDofUpdaterClearAndReinitialize, KratosTrilinosApplicationMPITestSuite)
{
    const auto& r_comm = Testing::GetDefaultDataCommunicator();

    Model model;
    const int local_size = 2;
    DofUpdaterTestSetup setup(model, r_comm, local_size);

    TrilinosDofUpdaterExperimental<TrilinosSparseSpaceType> dof_updater;
    dof_updater.Initialize(setup.DofSet, *setup.pDx);
    dof_updater.UpdateDofs(setup.DofSet, *setup.pDx);
    for (const auto& r_dof : setup.DofSet) {
        KRATOS_EXPECT_DOUBLE_EQ(r_dof.GetSolutionStepValue(), static_cast<double>(r_dof.EquationId()));
    }

    // After Clear, UpdateDofs must re-initialize transparently
    dof_updater.Clear();
    dof_updater.UpdateDofs(setup.DofSet, *setup.pDx);
    for (const auto& r_dof : setup.DofSet) {
        KRATOS_EXPECT_DOUBLE_EQ(r_dof.GetSolutionStepValue(), 2.0 * static_cast<double>(r_dof.EquationId()));
    }
}

} // namespace Kratos::Testing
