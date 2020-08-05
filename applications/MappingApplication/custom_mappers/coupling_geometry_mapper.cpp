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

namespace Kratos {
namespace Internals {

typedef typename MapperDefinitions::SparseSpaceType SparseSpaceType;
typedef std::size_t SizeType;

void Assemble(
    const MapperLocalSystem::MatrixType& rLocalMappingMatrix,
    const MapperLocalSystem::EquationIdVectorType& rOriginIds,
    const MapperLocalSystem::EquationIdVectorType& rDestinationIds,
    Matrix& rMappingMatrix)
{
    KRATOS_DEBUG_ERROR_IF(rLocalMappingMatrix.size1() != rDestinationIds.size()) << "MappingMatrixAssembly: DestinationID vector size mismatch: LocalMappingMatrix-Size1: " << rLocalMappingMatrix.size1() << " | DestinationIDs-size: " << rDestinationIds.size() << std::endl;
    KRATOS_DEBUG_ERROR_IF(rLocalMappingMatrix.size2() != rOriginIds.size()) << "MappingMatrixAssembly: OriginID vector size mismatch: LocalMappingMatrix-Size2: " << rLocalMappingMatrix.size2() << " | OriginIDs-size: " << rOriginIds.size() << std::endl;

    for (IndexType i=0; i<rDestinationIds.size(); ++i) {
        for (IndexType j=0; j<rOriginIds.size(); ++j) {
            rMappingMatrix(rDestinationIds[i], rOriginIds[j]) += rLocalMappingMatrix(i,j);
        }
    }
}

// dirty copy-paste from "mapping_matrix_utilities.cpp"
void InitializeSystemVector(Kratos::unique_ptr<typename SparseSpaceType::VectorType>& rpVector,
                            const SizeType VectorSize)
{
    // The vectors dont have graphs, that why we don't always have to reinitialize them
    if (rpVector == nullptr || rpVector->size() != VectorSize) { //if the pointer is not initialized initialize it to an empty vector
        Kratos::unique_ptr<typename SparseSpaceType::VectorType> p_new_vector = Kratos::make_unique<typename SparseSpaceType::VectorType>(VectorSize);
        rpVector.swap(p_new_vector);
    }
    else {
        SparseSpaceType::SetToZero(*rpVector);
    }
}

} // namespace Internals

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
    auto integration_point_values_slave = r_geometry_slave.IntegrationPoints();
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
    // @tteschemachen here kann man theoretisch auch das Origin-MP nehmen
    // kommt drauf an, wo die Coupling-Geometries sind

    // compose local element mappings
    CreateMapperLocalSystems(mpCouplingMP->GetCommunicator(),
                             mMapperLocalSystems);

    AssignInterfaceEquationIds(); // Has to be done every time in case of overlapping interfaces!

    // assemble projector interface mass matrix - interface_matrix_projector
    const std::size_t num_nodes_interface_slave = mpCouplingInterfaceDestination->NumberOfNodes();
    const std::size_t num_nodes_interface_master = mpCouplingInterfaceOrigin->NumberOfNodes();
    Matrix interface_matrix_projector = ZeroMatrix(num_nodes_interface_slave, num_nodes_interface_master);

    MapperLocalSystem::MatrixType local_mapping_matrix;
    MapperLocalSystem::EquationIdVectorType origin_ids;
    MapperLocalSystem::EquationIdVectorType destination_ids;

    for (size_t local_projector_system = 0;
        local_projector_system < mMapperLocalSystems.size()/2; ++local_projector_system) {
        mMapperLocalSystems[local_projector_system]->PairingInfo(0);
        mMapperLocalSystems[local_projector_system]->CalculateLocalSystem(local_mapping_matrix, origin_ids, destination_ids);
        Internals::Assemble(local_mapping_matrix, origin_ids, destination_ids, interface_matrix_projector);
    }

    // assemble slave interface mass matrix - interface_matrix_slave
    // TODO for dual mortar this should be a vector not a matrix
    Matrix interface_matrix_slave = ZeroMatrix(num_nodes_interface_slave, num_nodes_interface_slave);
    for (size_t local_projector_system = mMapperLocalSystems.size() / 2;
        local_projector_system < mMapperLocalSystems.size(); ++local_projector_system)
    {
        mMapperLocalSystems[local_projector_system]->PairingInfo(0);
        mMapperLocalSystems[local_projector_system]->CalculateLocalSystem(local_mapping_matrix, origin_ids, destination_ids);
        Internals::Assemble(local_mapping_matrix, origin_ids, destination_ids, interface_matrix_slave);
    }

    // Perform consistency scaling if requested
    if (mMapperSettings["consistency_scaling"].GetBool())
        EnforceConsistencyWithScaling(interface_matrix_slave, interface_matrix_projector, 1.1);

    // get total interface mapping matrix
    Matrix inv_interface_matrix_slave(num_nodes_interface_slave, num_nodes_interface_slave, 0.0);
    if (mMapperSettings["dual_mortar"].GetBool()) {
        for (size_t i = 0; i < interface_matrix_slave.size1(); ++i) {
            if (interface_matrix_slave(i, i) > std::numeric_limits<double>::epsilon()) {
                inv_interface_matrix_slave(i, i) = 1.0 / interface_matrix_slave(i, i);
            }
        }
    }
    else {
        double aux_det_slave = 0;
        MathUtils<double>::InvertMatrix(interface_matrix_slave, inv_interface_matrix_slave, aux_det_slave);
    }
    mpMappingMatrix = Kratos::make_unique<DenseMappingMatrixType>(prod(inv_interface_matrix_slave, interface_matrix_projector));
    CheckMappingMatrixConsistency();

    Internals::InitializeSystemVector(mpInterfaceVectorContainerOrigin->pGetVector(), num_nodes_interface_master);
    Internals::InitializeSystemVector(mpInterfaceVectorContainerDestination->pGetVector(), num_nodes_interface_slave);
}

template<class TSparseSpace, class TDenseSpace>
void CouplingGeometryMapper<TSparseSpace, TDenseSpace>::MapInternal(
    const Variable<double>& rOriginVariable,
    const Variable<double>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    mpInterfaceVectorContainerOrigin->UpdateSystemVectorFromModelPart(rOriginVariable, MappingOptions);

    TSparseSpace::Mult(
        *mpMappingMatrix,
        mpInterfaceVectorContainerOrigin->GetVector(),
        mpInterfaceVectorContainerDestination->GetVector()); // rQd = rMdo * rQo

    mpInterfaceVectorContainerDestination->UpdateModelPartFromSystemVector(rDestinationVariable, MappingOptions);
}

template<class TSparseSpace, class TDenseSpace>
void CouplingGeometryMapper<TSparseSpace, TDenseSpace>::MapInternalTranspose(
    const Variable<double>& rOriginVariable,
    const Variable<double>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    mpInterfaceVectorContainerDestination->UpdateSystemVectorFromModelPart(rDestinationVariable, MappingOptions);

    TSparseSpace::TransposeMult(
        *mpMappingMatrix,
        mpInterfaceVectorContainerDestination->GetVector(),
        mpInterfaceVectorContainerOrigin->GetVector()); // rQo = rMdo^T * rQd

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
    const Matrix& rInterfaceMatrixSlave,
    Matrix& rInterfaceMatrixProjected,
    const double scalingLimit)
{
    // Performs scaling of projected mapping entries as per eqn25 Wang2016
    double row_sum_slave;
    double row_sum_projector;

    for (IndexType i = 0; i < rInterfaceMatrixSlave.size1(); ++i) {
        row_sum_slave = 0.0;
        row_sum_projector = 0.0;

        for (IndexType j = 0; j < rInterfaceMatrixSlave.size2(); ++j)
            row_sum_slave += rInterfaceMatrixSlave(i, j);

        for (IndexType j = 0; j < rInterfaceMatrixProjected.size2(); ++j)
            row_sum_projector += rInterfaceMatrixProjected(i, j);

        const double alpha = (row_sum_slave / row_sum_projector < scalingLimit)
            ? row_sum_slave / row_sum_projector : scalingLimit;
        for (IndexType j = 0; j < rInterfaceMatrixProjected.size2(); ++j)
                rInterfaceMatrixProjected(i, j) *= alpha;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class CouplingGeometryMapper< MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType >;

}  // namespace Kratos.
