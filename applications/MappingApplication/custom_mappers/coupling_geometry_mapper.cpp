//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//                   Tobias Teschemachen
//

// System includes

// External includes

// Project includes
#include "coupling_geometry_mapper.h"
#include "mapping_application_variables.h"
#include "mappers/mapper_define.h"
#include "custom_utilities/mapping_matrix_utilities.h"
#include "custom_utilities/mapper_utilities.h"
#include "utilities/variable_utils.h"
#include "factories/linear_solver_factory.h"
#include "utilities/sparse_matrix_multiplication_utility.h"

namespace Kratos {

void CouplingGeometryLocalSystem::CalculateAll(MatrixType& rLocalMappingMatrix,
                    EquationIdVectorType& rOriginIds,
                    EquationIdVectorType& rDestinationIds,
                    MapperLocalSystem::PairingStatus& rPairingStatus) const
{
    const IndexType slave_index = (mIsDestinationIsSlave) ? 1 : 0;
    const IndexType master_index = 1 - slave_index;

    const auto& r_geometry_master = (mIsProjection)
        ? mpGeom->GetGeometryPart(master_index) // set to master  - get projected 'mass' matrix
        : mpGeom->GetGeometryPart(slave_index); // set to slave - get consistent slave 'mass' matrix
    const auto& r_geometry_slave = mpGeom->GetGeometryPart(slave_index);

    const bool is_dual_mortar = (!mIsProjection && mIsDualMortar)
        ? true
        : false;

    const std::size_t number_of_nodes_master = r_geometry_master.size();
    const std::size_t number_of_nodes_slave = r_geometry_slave.size();

    rPairingStatus = MapperLocalSystem::PairingStatus::InterfaceInfoFound;

    if (rLocalMappingMatrix.size1() != number_of_nodes_slave || rLocalMappingMatrix.size2() != number_of_nodes_master) {
        rLocalMappingMatrix.resize(number_of_nodes_slave, number_of_nodes_master, false);
    }
    if (rOriginIds.size()      != number_of_nodes_master) rOriginIds.resize(number_of_nodes_master);
    if (rDestinationIds.size() != number_of_nodes_slave) rDestinationIds.resize(number_of_nodes_slave);

    auto sf_values_master = r_geometry_master.ShapeFunctionsValues();
    auto sf_values_slave = r_geometry_slave.ShapeFunctionsValues();
    Vector det_jacobian;
    r_geometry_slave.DeterminantOfJacobian(det_jacobian);
    KRATOS_ERROR_IF(det_jacobian.size() != 1)
        << "Coupling Geometry Mapper should only have 1 integration point coupling per local system" << std::endl;
    const double weight = r_geometry_slave.IntegrationPoints()[0].Weight();

    if (is_dual_mortar) {
        rLocalMappingMatrix.clear();
        for (IndexType integration_point_itr = 0; integration_point_itr < sf_values_slave.size1(); ++integration_point_itr) {
            for (IndexType i = 0; i < sf_values_slave.size2(); ++i) {
                rLocalMappingMatrix(i, i) = sf_values_slave(integration_point_itr, i)
                    * det_jacobian[integration_point_itr] * weight;
                KRATOS_DEBUG_ERROR_IF(sf_values_slave(integration_point_itr, i) < 0.0)
                    << "DESTINATION SHAPE FUNCTIONS LESS THAN ZERO" << std::endl;
            }
        }
    }
    else {
        KRATOS_DEBUG_ERROR_IF(sf_values_slave.size1() != sf_values_master.size1())
            << "Coupling Geometry Mapper | origin and destination shape functions have different first sizes!"
            << "\nOrigin shape functions =\n\t" << sf_values_master
            << "\nDestination shape functions =\n\t" << sf_values_slave << std::endl;
        for (IndexType integration_point_itr = 0; integration_point_itr < sf_values_slave.size1(); ++integration_point_itr) {
            for (IndexType i = 0; i < sf_values_slave.size2(); ++i) {

                KRATOS_DEBUG_ERROR_IF(sf_values_slave(integration_point_itr, i) < 0.0)
                    << "DESTINATION SHAPE FUNCTIONS LESS THAN ZERO\n" << sf_values_slave << std::endl;

                for (IndexType j = 0; j < sf_values_master.size2(); ++j) {
                    rLocalMappingMatrix(i, j) = sf_values_slave(integration_point_itr, i)
                        * sf_values_master(integration_point_itr, j)
                        * det_jacobian[integration_point_itr] * weight;

                    KRATOS_DEBUG_ERROR_IF(sf_values_master(integration_point_itr, j) < 0.0)
                        << "ORIGIN SHAPE FUNCTIONS LESS THAN ZERO\n" << sf_values_master << std::endl;
                }
            }
        }
    }

    for (IndexType i=0; i< sf_values_master.size2(); ++i) {
        rOriginIds[i] = r_geometry_master[i].GetValue(INTERFACE_EQUATION_ID);
    }
    for (IndexType i=0; i< sf_values_slave.size2(); ++i) {
        rDestinationIds[i] = r_geometry_slave[i].GetValue(INTERFACE_EQUATION_ID);
    }
}

template<class TSparseSpace, class TDenseSpace>
CouplingGeometryMapper<TSparseSpace, TDenseSpace>::CouplingGeometryMapper(
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    Parameters JsonParameters):
        mrModelPartOrigin(rModelPartOrigin),
        mrModelPartDestination(rModelPartDestination),
        mMapperSettings(JsonParameters)
{
    JsonParameters.ValidateAndAssignDefaults(GetMapperDefaultSettings());
    const bool destination_is_slave = mMapperSettings["destination_is_slave"].GetBool();

    mpModeler = (ModelerFactory::Create(
        mMapperSettings["modeler_name"].GetString(),
        rModelPartOrigin.GetModel(),
        mMapperSettings["modeler_parameters"]));

    // adds destination model part
    mpModeler->GenerateNodes(rModelPartDestination);
    mpModeler->SetupGeometryModel();
    mpModeler->PrepareGeometryModel();

    // here use whatever ModelPart(s) was created by the Modeler
    mpCouplingMP = &(rModelPartOrigin.GetModel().GetModelPart("coupling"));

    mpCouplingInterfaceMaster = (destination_is_slave)
        ? mpCouplingMP->pGetSubModelPart("interface_origin")
        : mpCouplingMP->pGetSubModelPart("interface_destination");
    mpCouplingInterfaceSlave = (destination_is_slave)
        ? mpCouplingMP->pGetSubModelPart("interface_destination")
        : mpCouplingMP->pGetSubModelPart("interface_origin");

    mpInterfaceVectorContainerMaster = Kratos::make_unique<InterfaceVectorContainerType>(*mpCouplingInterfaceMaster);
    mpInterfaceVectorContainerSlave = Kratos::make_unique<InterfaceVectorContainerType>(*mpCouplingInterfaceSlave);

    this->CreateLinearSolver();
    this->InitializeInterface();
}


template<class TSparseSpace, class TDenseSpace>
void CouplingGeometryMapper<TSparseSpace, TDenseSpace>::InitializeInterface(Kratos::Flags MappingOptions)
{
    // compose local element mappings
    const bool dual_mortar = mMapperSettings["dual_mortar"].GetBool();
    const bool precompute_mapping_matrix = mMapperSettings["precompute_mapping_matrix"].GetBool();
    const bool direct_map_to_destination = mMapperSettings["destination_is_slave"].GetBool();
    CouplingGeometryLocalSystem ref_projector_local_system(nullptr, true, dual_mortar, direct_map_to_destination);
    CouplingGeometryLocalSystem ref_slave_local_system(nullptr, false, dual_mortar, direct_map_to_destination);

    MapperUtilities::CreateMapperLocalSystemsFromGeometries(ref_projector_local_system,
                             mpCouplingMP->GetCommunicator(),
                             mMapperLocalSystemsProjector);

    MapperUtilities::CreateMapperLocalSystemsFromGeometries(ref_slave_local_system,
                             mpCouplingMP->GetCommunicator(),
                             mMapperLocalSystemsSlave);

    AssignInterfaceEquationIds(); // Has to be done every time in case of overlapping interfaces!

    // assemble projector interface mass matrix - interface_matrix_projector
    const std::size_t num_nodes_interface_slave = mpCouplingInterfaceSlave->NumberOfNodes();
    const std::size_t num_nodes_interface_master = mpCouplingInterfaceMaster->NumberOfNodes();
    mpMappingMatrix = Kratos::make_unique<MappingMatrixType>(num_nodes_interface_slave, num_nodes_interface_master);

    // TODO Philipp I am pretty sure we should separate the vector construction from the matrix construction, should be independent otherwise no clue what is happening
    MappingMatrixUtilities<TSparseSpace, TDenseSpace>::BuildMappingMatrix(
        mpMappingMatrixSlave,
        mpInterfaceVectorContainerSlave->pGetVector(),
        mpInterfaceVectorContainerSlave->pGetVector(),
        mpInterfaceVectorContainerSlave->GetModelPart(),
        mpInterfaceVectorContainerSlave->GetModelPart(),
        mMapperLocalSystemsSlave,
        0); // The echo-level is no longer needed here, refactor in separate PR

    MappingMatrixUtilities<TSparseSpace, TDenseSpace>::BuildMappingMatrix(
        mpMappingMatrixProjector,
        mpInterfaceVectorContainerMaster->pGetVector(),
        mpInterfaceVectorContainerSlave->pGetVector(),
        mpInterfaceVectorContainerMaster->GetModelPart(),
        mpInterfaceVectorContainerSlave->GetModelPart(),
        mMapperLocalSystemsProjector,
        0); // The echo-level is no longer needed here, refactor in separate PR

    // Perform consistency scaling if requested
    if (mMapperSettings["consistency_scaling"].GetBool()) {
        EnforceConsistencyWithScaling(*mpMappingMatrixSlave, *mpMappingMatrixProjector, 1.1);
    }

    // get total interface mapping matrix
    if (dual_mortar) {
        // invert the diagonal entries through the CSR arrays (backend-generic)
        auto& r_slave_matrix = *mpMappingMatrixSlave;
        auto row_indices = r_slave_matrix.index1_data();
        auto col_indices = r_slave_matrix.index2_data();
        auto values = r_slave_matrix.value_data();
        for (IndexType i = 0; i < r_slave_matrix.size1(); ++i) {
            for (auto k = row_indices[i]; k < row_indices[i + 1]; ++k) {
                if (static_cast<IndexType>(col_indices[k]) == i) {
                    double& val = values[k];
                    if (val > std::numeric_limits<double>::epsilon()) {
                        val = 1.0 / val;
                    } else {
                        val = 0.0;
                    }
                }
            }
        }

        SparseMatrixMultiplicationUtility::MatrixMultiplication(*mpMappingMatrixSlave, *mpMappingMatrixProjector, *mpMappingMatrix);
    }
    else {
        MappingMatrixUtilities<TSparseSpace, TDenseSpace>::InitializeSystemVector(mpTempVector, mpInterfaceVectorContainerSlave->GetModelPart().NumberOfNodes());
        if (precompute_mapping_matrix)  CalculateMappingMatrixWithSolver(*mpMappingMatrixSlave, *mpMappingMatrixProjector);
    }

    // Check row sum of pre-computed mapping matrices only
    if (precompute_mapping_matrix || dual_mortar) {
        const std::string base_file_name = "O_" + mrModelPartOrigin.Name() + "__D_" + mrModelPartDestination.Name() + ".mm";
        const double row_sum_tolerance = mMapperSettings["row_sum_tolerance"].GetDouble();
        MappingMatrixUtilities<TSparseSpace, TDenseSpace>::CheckRowSum(*mpMappingMatrix, base_file_name, true, row_sum_tolerance);
    }
}

template<class TSparseSpace, class TDenseSpace>
void CouplingGeometryMapper<TSparseSpace, TDenseSpace>::MapInternal(
    const Variable<double>& rOriginVariable,
    const Variable<double>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    const bool dual_mortar = mMapperSettings["dual_mortar"].GetBool();
    const bool precompute_mapping_matrix = mMapperSettings["precompute_mapping_matrix"].GetBool();

    mpInterfaceVectorContainerMaster->UpdateSystemVectorFromModelPart(rOriginVariable, MappingOptions);

    if (dual_mortar || precompute_mapping_matrix) {
        TSparseSpace::Mult(
            *mpMappingMatrix,
            mpInterfaceVectorContainerMaster->GetVector(),
            mpInterfaceVectorContainerSlave->GetVector()); // rQd = rMdo * rQo
    } else {
        TSparseSpace::Mult(
            *mpMappingMatrixProjector,
            mpInterfaceVectorContainerMaster->GetVector(),
            *mpTempVector); // rQd = rMdo * rQo

        mpLinearSolver->Solve(*mpMappingMatrixSlave, mpInterfaceVectorContainerSlave->GetVector(), *mpTempVector);
    }

    mpInterfaceVectorContainerSlave->UpdateModelPartFromSystemVector(rDestinationVariable, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void CouplingGeometryMapper<TSparseSpace, TDenseSpace>::MapInternalTranspose(
    const Variable<double>& rOriginVariable,
    const Variable<double>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    const bool dual_mortar = mMapperSettings["dual_mortar"].GetBool();
    const bool precompute_mapping_matrix = mMapperSettings["precompute_mapping_matrix"].GetBool();

    mpInterfaceVectorContainerSlave->UpdateSystemVectorFromModelPart(rDestinationVariable, MappingOptions);

    if (dual_mortar || precompute_mapping_matrix) {
        TSparseSpace::TransposeMult(
            *mpMappingMatrix,
            mpInterfaceVectorContainerSlave->GetVector(),
            mpInterfaceVectorContainerMaster->GetVector()); // rQo = rMdo^T * rQd
    } else {
        mpLinearSolver->Solve(*mpMappingMatrixSlave, *mpTempVector, mpInterfaceVectorContainerSlave->GetVector());

        TSparseSpace::TransposeMult(
            *mpMappingMatrixProjector,
            *mpTempVector,
            mpInterfaceVectorContainerMaster->GetVector()); // rQo = rMdo^T * rQd
    }

    mpInterfaceVectorContainerMaster->UpdateModelPartFromSystemVector(rOriginVariable, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void CouplingGeometryMapper<TSparseSpace, TDenseSpace>::MapInternal(
    const Variable<array_1d<double, 3>>& rOriginVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    for (const auto var_ext : {"_X", "_Y", "_Z"}) {
        const auto& var_origin = KratosComponents<Variable<double>>::Get(rOriginVariable.Name() + var_ext);
        const auto& var_destination = KratosComponents<Variable<double>>::Get(rDestinationVariable.Name() + var_ext);

        MapInternal(var_origin, var_destination, MappingOptions);
    }
}

template<class TSparseSpace, class TDenseSpace>
void CouplingGeometryMapper<TSparseSpace, TDenseSpace>::MapInternalTranspose(
    const Variable<array_1d<double, 3>>& rOriginVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    for (const auto var_ext : {"_X", "_Y", "_Z"}) {
        const auto& var_origin = KratosComponents<Variable<double>>::Get(rOriginVariable.Name() + var_ext);
        const auto& var_destination = KratosComponents<Variable<double>>::Get(rDestinationVariable.Name() + var_ext);

        MapInternalTranspose(var_origin, var_destination, MappingOptions);
    }
}

template<class TSparseSpace, class TDenseSpace>
void CouplingGeometryMapper<TSparseSpace, TDenseSpace>::EnforceConsistencyWithScaling(
    const MappingMatrixType& rInterfaceMatrixSlave,
    MappingMatrixType& rInterfaceMatrixProjected,
    const double scalingLimit)
{
    // Performs scaling of projected mapping entries as per eqn25 Wang2016

    // Get row sum vector of slave matrix
    SparseSpaceType::VectorType unit_vector(SparseSpaceType::Size2(rInterfaceMatrixSlave));
    SparseSpaceType::Set(unit_vector, 1.0);
    SparseSpaceType::VectorType slave_row_sums_vector(SparseSpaceType::Size1(rInterfaceMatrixSlave));
    SparseSpaceType::Mult(rInterfaceMatrixSlave, unit_vector, slave_row_sums_vector);

    // Get row sum vector of projected matrix
    unit_vector.resize(SparseSpaceType::Size2(rInterfaceMatrixProjected));
    SparseSpaceType::Set(unit_vector, 1.0);
    SparseSpaceType::VectorType projected_row_sums_vector(SparseSpaceType::Size1(rInterfaceMatrixProjected));
    SparseSpaceType::Mult(rInterfaceMatrixProjected, unit_vector, projected_row_sums_vector);

    // Loop over sparse rows of projected matrix and correct entries if
    // needed, through the CSR arrays (backend-generic)
    auto row_indices = rInterfaceMatrixProjected.index1_data();
    auto values = rInterfaceMatrixProjected.value_data();
    for (IndexType row_counter = 0; row_counter < rInterfaceMatrixProjected.size1(); ++row_counter)
    {
        if (std::abs(slave_row_sums_vector[row_counter]/ projected_row_sums_vector[row_counter] - 1.0) > 1e-15) {
            // Correct entries
            const double alpha = (slave_row_sums_vector[row_counter] / projected_row_sums_vector[row_counter] < scalingLimit)
                ? slave_row_sums_vector[row_counter] / projected_row_sums_vector[row_counter] : scalingLimit;
            for (auto k = row_indices[row_counter]; k < row_indices[row_counter + 1]; ++k) values[k] *= alpha;
        }
    }
}

template<class TSparseSpace, class TDenseSpace>
void CouplingGeometryMapper<TSparseSpace, TDenseSpace>::CalculateMappingMatrixWithSolver(
    MappingMatrixType& rConsistentInterfaceMatrix, MappingMatrixType& rProjectedInterfaceMatrix)
{
    // the solver-based mapping matrix is dense by construction, so the
    // compressed storage is allocated and structured up front (insert_element
    // is a uBLAS-only slow path)
    const size_t n_mm_rows = rConsistentInterfaceMatrix.size1();
    const size_t n_mm_cols = rProjectedInterfaceMatrix.size2();
    mpMappingMatrix = Kratos::make_unique<typename SparseSpaceType::MatrixType>(n_mm_rows, n_mm_cols, n_mm_rows * n_mm_cols);
    {
        auto row_indices = mpMappingMatrix->index1_data();
        auto col_indices = mpMappingMatrix->index2_data();
        row_indices[0] = 0;
        for (size_t j = 0; j < n_mm_rows; ++j) {
            for (size_t i = 0; i < n_mm_cols; ++i) col_indices[j * n_mm_cols + i] = i;
            row_indices[j + 1] = (j + 1) * n_mm_cols;
        }
        mpMappingMatrix->set_filled(n_mm_rows + 1, n_mm_rows * n_mm_cols);
    }

    const size_t n_rows = mpMappingMatrix->size1();
    const size_t n_cols = mpMappingMatrix->size2();
    typename TSparseSpace::VectorType solution(n_rows);
    typename TSparseSpace::VectorType projector_column(n_rows);
    auto values = mpMappingMatrix->value_data();

    for (size_t i = 0; i < n_cols; ++i)
    {
        for (size_t j = 0; j < n_rows; ++j) projector_column[j] = rProjectedInterfaceMatrix(j, i); // TODO try boost slice or project
        mpLinearSolver->Solve(rConsistentInterfaceMatrix, solution, projector_column);
        // the compressed storage was structured dense up front
        for (size_t j = 0; j < n_rows; ++j) values[j * n_cols + i] = solution[j];
    }
}

template<class TSparseSpace, class TDenseSpace>
void CouplingGeometryMapper<TSparseSpace, TDenseSpace>::CreateLinearSolver()
{
    if (mMapperSettings["linear_solver_settings"].Has("solver_type")) {
        mpLinearSolver = LinearSolverFactory<TSparseSpace, TDenseSpace>().Create(mMapperSettings["linear_solver_settings"]);
    }
    else {
        // TODO - replicate 'get fastest solver'
        mMapperSettings.AddString("solver_type", "skyline_lu_factorization");
        mpLinearSolver = LinearSolverFactory<TSparseSpace, TDenseSpace>().Create(mMapperSettings);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class CouplingGeometryMapper< MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType >;

}  // namespace Kratos.
