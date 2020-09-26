//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//                   Tobias Teschemachen
//

// System includes

// External includes

// Project includes
#include "coupling_geometry_mapper.h"
#include "mapping_application_variables.h"
#include "custom_utilities/mapper_typedefs.h"
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
    const auto& r_geometry_master = (mIsProjection)
        ? mpGeom->GetGeometryPart(0) // set to master  - get projected 'mass' matrix
        : mpGeom->GetGeometryPart(1); // set to slave - get consistent slave 'mass' matrix
    const auto& r_geometry_slave = mpGeom->GetGeometryPart(1);

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

    if (is_dual_mortar) {
        rLocalMappingMatrix.clear();
        for (IndexType integration_point_itr = 0; integration_point_itr < sf_values_slave.size1(); ++integration_point_itr) {
            for (IndexType i = 0; i < sf_values_slave.size2(); ++i) {
                rLocalMappingMatrix(i, i) = sf_values_slave(integration_point_itr, i)
                    * det_jacobian[integration_point_itr];
                KRATOS_DEBUG_ERROR_IF(sf_values_master(integration_point_itr, i) < 0.0)
                    << "SHAPE FUNCTIONS LESS THAN ZERO" << std::endl;
            }
        }
    }
    else {
        for (IndexType integration_point_itr = 0; integration_point_itr < sf_values_slave.size1(); ++integration_point_itr) {
            for (IndexType i = 0; i < sf_values_slave.size2(); ++i) {
                for (IndexType j = 0; j < sf_values_master.size2(); ++j) {
                    rLocalMappingMatrix(i, j) = sf_values_slave(integration_point_itr, i)
                        * sf_values_master(integration_point_itr, j)
                        * det_jacobian[integration_point_itr];

                    KRATOS_DEBUG_ERROR_IF(sf_values_master(integration_point_itr, i) < 0.0)
                        << "SHAPE FUNCTIONS LESS THAN ZERO\n" << sf_values_master << std::endl;
                    KRATOS_DEBUG_ERROR_IF(sf_values_slave(integration_point_itr, j) < 0.0)
                        << "SHAPE FUNCTIONS LESS THAN ZERO\n" << sf_values_slave << std::endl;
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

std::string CouplingGeometryLocalSystem::PairingInfo(const int EchoLevel) const
{
    //std::cout << "   >>> XXX : IS_PROJECTED_LOCAL_SYSTEM: " << mIsProjection << std::endl;
    // KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

    // std::stringstream buffer;
    // buffer << "CouplingGeometryLocalSystem based on " << mpNode->Info();
    // if (EchoLevel > 1) { // TODO leave here?
    //     buffer << " at Coodinates " << Coordinates()[0] << " | " << Coordinates()[1] << " | " << Coordinates()[2];
    //     if (mPairingStatus == MapperLocalSystem::PairingStatus::Approximation) {
    //         mpNode->SetValue(PAIRING_STATUS, 0);
    //     } else {
    //         mpNode->SetValue(PAIRING_STATUS, -1);
    //     }
    // }
    // return buffer.str();
    return "";
}



template<class TSparseSpace, class TDenseSpace>
void CouplingGeometryMapper<TSparseSpace, TDenseSpace>::InitializeInterface(Kratos::Flags MappingOptions)
{
    // compose local element mappings
    const bool dual_mortar = mMapperSettings["dual_mortar"].GetBool();
    CouplingGeometryLocalSystem ref_projector_local_system(nullptr, true, dual_mortar);
    CouplingGeometryLocalSystem ref_slave_local_system(nullptr, false, dual_mortar);

    MapperUtilities::CreateMapperLocalSystemsFromGeometries(ref_projector_local_system,
                             mpCouplingMP->GetCommunicator(),
                             mMapperLocalSystemsProjector);

    MapperUtilities::CreateMapperLocalSystemsFromGeometries(ref_slave_local_system,
                             mpCouplingMP->GetCommunicator(),
                             mMapperLocalSystemsSlave);

    AssignInterfaceEquationIds(); // Has to be done every time in case of overlapping interfaces!

    // assemble projector interface mass matrix - interface_matrix_projector
    const std::size_t num_nodes_interface_slave = mpCouplingInterfaceDestination->NumberOfNodes();
    const std::size_t num_nodes_interface_master = mpCouplingInterfaceOrigin->NumberOfNodes();
    mpMappingMatrix = Kratos::make_unique<MappingMatrixType>(num_nodes_interface_slave, num_nodes_interface_master);

    const int echo_level = mMapperSettings["echo_level"].GetInt();

    // TODO Philipp I am pretty sure we should separate the vector construction from the matrix construction, should be independent otherwise no clue what is happening
    MappingMatrixUtilities::BuildMappingMatrix<TSparseSpace, TDenseSpace>(
        mpMappingMatrixSlave,
        mpInterfaceVectorContainerDestination->pGetVector(),
        mpInterfaceVectorContainerDestination->pGetVector(),
        mpInterfaceVectorContainerDestination->GetModelPart(),
        mpInterfaceVectorContainerDestination->GetModelPart(),
        mMapperLocalSystemsSlave,
        echo_level);

    MappingMatrixUtilities::BuildMappingMatrix<TSparseSpace, TDenseSpace>(
        mpMappingMatrixProjector,
        mpInterfaceVectorContainerOrigin->pGetVector(),
        mpInterfaceVectorContainerDestination->pGetVector(),
        mpInterfaceVectorContainerOrigin->GetModelPart(),
        mpInterfaceVectorContainerDestination->GetModelPart(),
        mMapperLocalSystemsProjector,
        echo_level);

    // Perform consistency scaling if requested
    if (mMapperSettings["consistency_scaling"].GetBool()) {
        EnforceConsistencyWithScaling(*mpMappingMatrixSlave, *mpMappingMatrixProjector, 1.1);
    }

    // get total interface mapping matrix
    if (dual_mortar) {
        typedef typename MappingMatrixType::iterator1 t_it_1;
        typedef typename MappingMatrixType::iterator2 t_it_2;

        for (t_it_1 it_1 = mpMappingMatrixSlave->begin1(); it_1 != mpMappingMatrixSlave->end1(); ++it_1) {
            for (t_it_2 it_2 = it_1.begin(); it_2 != it_1.end(); ++it_2) {
                if (it_2.index1() == it_2.index2()) {
                    double& val = *it_2;
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
        // @Peter we should make this optional, the alternative is to solve the system each time we map
        // => this is done in Empire
        // lets discuss in the next meeting
        // CalculateMappingMatrixWithSolver(*mpMappingMatrixSlave, *mpMappingMatrixProjector);

        MappingMatrixUtilities::InitializeSystemVector<TSparseSpace, TDenseSpace>(mpTempVector, mpInterfaceVectorContainerDestination->GetModelPart().NumberOfNodes());
    }

    // CheckMappingMatrixConsistency();
}

template<class TSparseSpace, class TDenseSpace>
void CouplingGeometryMapper<TSparseSpace, TDenseSpace>::MapInternal(
    const Variable<double>& rOriginVariable,
    const Variable<double>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    const bool dual_mortar = mMapperSettings["dual_mortar"].GetBool();

    mpInterfaceVectorContainerOrigin->UpdateSystemVectorFromModelPart(rOriginVariable, MappingOptions);

    if (dual_mortar) {
        TSparseSpace::Mult(
            *mpMappingMatrix,
            mpInterfaceVectorContainerOrigin->GetVector(),
            mpInterfaceVectorContainerDestination->GetVector()); // rQd = rMdo * rQo
    } else {
        TSparseSpace::Mult(
            *mpMappingMatrixProjector,
            mpInterfaceVectorContainerOrigin->GetVector(),
            *mpTempVector); // rQd = rMdo * rQo

        mpLinearSolver->Solve(*mpMappingMatrixSlave, mpInterfaceVectorContainerDestination->GetVector(), *mpTempVector);
    }

    mpInterfaceVectorContainerDestination->UpdateModelPartFromSystemVector(rDestinationVariable, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void CouplingGeometryMapper<TSparseSpace, TDenseSpace>::MapInternalTranspose(
    const Variable<double>& rOriginVariable,
    const Variable<double>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    const bool dual_mortar = mMapperSettings["dual_mortar"].GetBool();

    mpInterfaceVectorContainerDestination->UpdateSystemVectorFromModelPart(rDestinationVariable, MappingOptions);

    if (dual_mortar) {
        TSparseSpace::TransposeMult(
            *mpMappingMatrix,
            mpInterfaceVectorContainerDestination->GetVector(),
            mpInterfaceVectorContainerOrigin->GetVector()); // rQo = rMdo^T * rQd
    } else {
        mpLinearSolver->Solve(*mpMappingMatrixSlave, *mpTempVector, mpInterfaceVectorContainerDestination->GetVector());

        TSparseSpace::TransposeMult(
            *mpMappingMatrixProjector,
            *mpTempVector,
            mpInterfaceVectorContainerOrigin->GetVector()); // rQo = rMdo^T * rQd
    }

    mpInterfaceVectorContainerOrigin->UpdateModelPartFromSystemVector(rOriginVariable, MappingOptions);
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
    // // Performs scaling of projected mapping entries as per eqn25 Wang2016
    // for (IndexType i = 0; i < rInterfaceMatrixSlave.size1(); ++i) {
    //     double row_sum_slave = 0.0;
    //     double row_sum_projector = 0.0;

    //     for (IndexType j = 0; j < rInterfaceMatrixSlave.size2(); ++j)
    //         row_sum_slave += rInterfaceMatrixSlave(i, j);

    //     for (IndexType j = 0; j < rInterfaceMatrixProjected.size2(); ++j)
    //         row_sum_projector += rInterfaceMatrixProjected(i, j);

    //     const double alpha = (row_sum_slave / row_sum_projector < scalingLimit)
    //         ? row_sum_slave / row_sum_projector : scalingLimit;
    //     for (IndexType j = 0; j < rInterfaceMatrixProjected.size2(); ++j)
    //             rInterfaceMatrixProjected(i, j) *= alpha;
    // }
}

template<class TSparseSpace, class TDenseSpace>
void CouplingGeometryMapper<TSparseSpace, TDenseSpace>::CalculateMappingMatrixWithSolver(
    MappingMatrixType& rConsistentInterfaceMatrix, MappingMatrixType& rProjectedInterfaceMatrix)
{
    // mpMappingMatrix = Kratos::make_unique<DenseMappingMatrixType>(DenseMappingMatrixType(
    //     rConsistentInterfaceMatrix.size1(), rProjectedInterfaceMatrix.size2()));

    // const size_t n_rows = mpMappingMatrix->size1();
    // #pragma omp parallel
    // {
    //     Vector solution(n_rows);
    //     Vector projector_column(n_rows);

    //     #pragma omp for
    //     for (int i = 0; i < static_cast<int>(mpMappingMatrix->size2()); ++i)
    //     {
    //         for (size_t j = 0; j < n_rows; ++j) projector_column[j] = rProjectedInterfaceMatrix(j, i);
    //         mpLinearSolver->Solve(rConsistentInterfaceMatrix, solution, projector_column);
    //         for (size_t j = 0; j < n_rows; ++j) (*mpMappingMatrix)(j, i) = solution[j];
    //     }
    // }
}

template<class TSparseSpace, class TDenseSpace>
void CouplingGeometryMapper<TSparseSpace, TDenseSpace>::CheckMappingMatrixConsistency()
{
    // @Peter this is already included in "MappingMatrixUtilities::CheckRowSum"

    // for (size_t row = 0; row < mpMappingMatrix->size1(); ++row) {
    //     double row_sum = 0.0;
    //     for (size_t col = 0; col < mpMappingMatrix->size2(); ++col) row_sum += (*mpMappingMatrix)(row, col);
    //     if (std::abs(row_sum - 1.0) > 1e-12) {
    //         KRATOS_WATCH(*mpMappingMatrix)
    //         KRATOS_WATCH(row_sum)
    //         KRATOS_ERROR << "mapping matrix is not consistent\n";
    //     }
    // }
}

template<class TSparseSpace, class TDenseSpace>
void CouplingGeometryMapper<TSparseSpace, TDenseSpace>::CreateLinearSolver()
{
    bool is_linear_solver_specified = false;
    if (mMapperSettings.Has("linear_solver_settings"))
    {
        if (mMapperSettings["linear_solver_settings"].Has("solver_type"))
        {
            is_linear_solver_specified = true;
            mpLinearSolver = LinearSolverFactory<TSparseSpace, TDenseSpace>().Create(mMapperSettings["linear_solver_settings"]);
        }
    }
    if (!is_linear_solver_specified)
    {
        // TODO - replicate 'get fastest solver'
        mMapperSettings.AddString("solver_type", "skyline_lu_factorization");
        mpLinearSolver = LinearSolverFactory<TSparseSpace, TDenseSpace>().Create(mMapperSettings);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class CouplingGeometryMapper< MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType >;

}  // namespace Kratos.
