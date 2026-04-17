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

/* Trilinos includes */
#include "Epetra_FEVector.h"

// Project includes
#include "tests/cpp_tests/trilinos_fast_suite.h"
#include "containers/model.h"
#include "custom_strategies/convergencecriterias/trilinos_mixed_generic_criteria.h"
#include "mpi/includes/mpi_data_communicator.h"
#include "mpi/utilities/model_part_communicator_utilities.h"
#include "mpi/utilities/parallel_fill_communicator.h"
#include "spaces/ublas_space.h"
#include "trilinos_space.h"

namespace Kratos::Testing
{

using DofsArrayType = ModelPart::DofsArrayType;
using TrilinosLocalSpaceType = UblasSpace<double, Matrix, Vector>;
using TrilinosSparseSpaceType = TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector>;
using TrilinosMixedGenericCriteriaType = TrilinosMixedGenericCriteria<TrilinosSparseSpaceType, TrilinosLocalSpaceType>;
using ConvergenceVariableListType = typename TrilinosMixedGenericCriteriaType::ConvergenceVariableListType;

void GenerateTestTrilinosMixedGenericCriteriaModelPart(
	ModelPart& rModelPart,
	const unsigned int NumberOfNodes = 10)
{
	const DataCommunicator& r_comm = Testing::GetDefaultDataCommunicator();

	KRATOS_ERROR_IF_NOT(r_comm.IsDistributed()) << "Only distributed DataCommunicators can be used!" << std::endl;

	ModelPartCommunicatorUtilities::SetMPICommunicator(rModelPart, r_comm);

	rModelPart.SetBufferSize(1);
	rModelPart.AddNodalSolutionStepVariable(PRESSURE);
	rModelPart.AddNodalSolutionStepVariable(VELOCITY);
	rModelPart.AddNodalSolutionStepVariable(PARTITION_INDEX);

	const int rank = r_comm.Rank();
	const int size = r_comm.Size();

	for (unsigned int i_node = NumberOfNodes * rank; i_node < NumberOfNodes * (rank + 1); ++i_node) {
		auto p_node = rModelPart.CreateNewNode(i_node + 1, 0.0, 0.0, 0.0);
		p_node->FastGetSolutionStepValue(PARTITION_INDEX) = rank;
		rModelPart.AddNode(p_node);
	}

	if (rank < (size - 1)) {
		auto p_node = rModelPart.CreateNewNode((NumberOfNodes * (rank + 1)) + 1, 0.0, 0.0, 0.0);
		p_node->FastGetSolutionStepValue(PARTITION_INDEX) = rank + 1;
		rModelPart.AddNode(p_node);
	}

	ParallelFillCommunicator(rModelPart, r_comm).Execute();
}

void SetupMixedCriteriaDofs(
	ModelPart& rModelPart,
	DofsArrayType& rAuxDofSet)
{
	rModelPart.GetNodalSolutionStepVariablesList().AddDof(&PRESSURE);
	rModelPart.GetNodalSolutionStepVariablesList().AddDof(&VELOCITY);

	for (auto& r_node : rModelPart.Nodes()) {
		r_node.AddDof(PRESSURE);
		r_node.AddDof(VELOCITY_X);
		r_node.AddDof(VELOCITY_Y);
	}

	rAuxDofSet.reserve(3 * rModelPart.NumberOfNodes());
	for (auto& r_node : rModelPart.Nodes()) {
		const IndexType base_id = 3 * (r_node.Id() - 1);
		auto p_pres_dof = r_node.pGetDof(PRESSURE);
		auto p_vel_x_dof = r_node.pGetDof(VELOCITY_X);
		auto p_vel_y_dof = r_node.pGetDof(VELOCITY_Y);

		p_pres_dof->SetEquationId(base_id);
		p_vel_x_dof->SetEquationId(base_id + 1);
		p_vel_y_dof->SetEquationId(base_id + 2);

		rAuxDofSet.push_back(p_pres_dof);
		rAuxDofSet.push_back(p_vel_x_dof);
		rAuxDofSet.push_back(p_vel_y_dof);
	}
	rAuxDofSet.Sort();
}

ConvergenceVariableListType CreateConvergenceSettings()
{
	ConvergenceVariableListType convergence_settings;
	convergence_settings.push_back(std::make_tuple(static_cast<VariableData*>(&PRESSURE), 1.0e-3, 1.0e-5));
	convergence_settings.push_back(std::make_tuple(static_cast<VariableData*>(&VELOCITY), 1.0e-3, 1.0e-5));
	return convergence_settings;
}

template<class TSparseSpaceType>
void SetTestValues(
	ModelPart& rModelPart,
	typename TSparseSpaceType::VectorType& rDx,
	const double ScaleFactor)
{
	const auto& r_comm = rModelPart.GetCommunicator().GetDataCommunicator();
	const int rank = r_comm.Rank();

	const double pressure_val = 1.0e-3;
	array_1d<double, 3> velocity_val = ZeroVector(3);
	velocity_val[0] = 2.0e-3;
	velocity_val[1] = 3.0e-3;

	for (auto& r_node : rModelPart.Nodes()) {
		r_node.FastGetSolutionStepValue(PRESSURE) = pressure_val;
		r_node.FastGetSolutionStepValue(VELOCITY) = velocity_val;

		if (r_node.FastGetSolutionStepValue(PARTITION_INDEX) == rank) {
			const IndexType base_id = 3 * (r_node.Id() - 1);
			TSparseSpaceType::SetValue(rDx, base_id, pressure_val * ScaleFactor);
			TSparseSpaceType::SetValue(rDx, base_id + 1, velocity_val[0] * ScaleFactor);
			TSparseSpaceType::SetValue(rDx, base_id + 2, velocity_val[1] * ScaleFactor);
		}
	}
}

KRATOS_TEST_CASE_IN_SUITE(TrilinosMixedGenericCriteria, KratosTrilinosApplicationMPITestSuite)
{
	Model current_model;
	ModelPart& r_model_part = current_model.CreateModelPart("TestModelPart");

	const unsigned int n_nodes = 10;
	GenerateTestTrilinosMixedGenericCriteriaModelPart(r_model_part, n_nodes);

	const auto& r_comm = r_model_part.GetCommunicator().GetDataCommunicator();
	const int rank = r_comm.Rank();

	DofsArrayType aux_dof_set;
	SetupMixedCriteriaDofs(r_model_part, aux_dof_set);

	auto mixed_generic_criteria = TrilinosMixedGenericCriteriaType(CreateConvergenceSettings());

	auto raw_mpi_comm = MPIDataCommunicator::GetMPICommunicator(r_comm);
	Epetra_MpiComm epetra_comm(raw_mpi_comm);

	const int n_local_dofs = static_cast<int>(3 * n_nodes);
	std::vector<int> local_dof_ids(n_local_dofs, 0);
	for (int i = 0; i < n_local_dofs; ++i) {
		local_dof_ids[i] = (3 * n_nodes * rank) + i;
	}
	Epetra_Map map(-1, n_local_dofs, local_dof_ids.data(), 0, epetra_comm);

	TrilinosSparseSpaceType::MatrixPointerType pA;
	TrilinosSparseSpaceType::MatrixType& rA = *pA;
	TrilinosSparseSpaceType::VectorType Dx(map);
	TrilinosSparseSpaceType::VectorPointerType pb;
	TrilinosSparseSpaceType::VectorType& rb = *pb;

	SetTestValues<TrilinosSparseSpaceType>(r_model_part, Dx, 1.0);
	mixed_generic_criteria.InitializeSolutionStep(r_model_part, aux_dof_set, rA, Dx, rb);

	bool convergence = mixed_generic_criteria.PostCriteria(r_model_part, aux_dof_set, rA, Dx, rb);
	KRATOS_EXPECT_FALSE(convergence);

	SetTestValues<TrilinosSparseSpaceType>(r_model_part, Dx, 1.0e-4);
	convergence = mixed_generic_criteria.PostCriteria(r_model_part, aux_dof_set, rA, Dx, rb);
	KRATOS_EXPECT_TRUE(convergence);
}

} // namespace Kratos::Testing
