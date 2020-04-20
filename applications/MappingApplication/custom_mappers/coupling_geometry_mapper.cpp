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

namespace Kratos
{

void CouplingGeometryLocalSystem::CalculateAll(MatrixType& rLocalMappingMatrix,
                    EquationIdVectorType& rOriginIds,
                    EquationIdVectorType& rDestinationIds,
                    MapperLocalSystem::PairingStatus& rPairingStatus) const
{
    const auto& r_geometry_master = (mpGeom->GetValue(IS_PROJECTED_LOCAL_SYSTEM))
        ? mpGeom->GetGeometryPart(0) // set to master  - get projected 'mass' matrix
        : mpGeom->GetGeometryPart(1); // set to slave - get consistent slave 'mass' matrix
    const auto& r_geometry_slave = mpGeom->GetGeometryPart(1);

    const bool is_dual_mortar = (mpGeom->GetValue(IS_PROJECTED_LOCAL_SYSTEM) &&
        mpGeom->GetValue(IS_DUAL_MORTAR))
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
        << "Coupling Geometry Mapper should only have 1 integration point coupling per local system"
        << std::endl;

    if (is_dual_mortar) {
        for (IndexType integration_point_itr = 0; integration_point_itr < sf_values_slave.size1(); ++integration_point_itr) { // This loop is probably redundant - it will always just be 1
            for (IndexType i = 0; i < sf_values_slave.size2(); ++i) {
                rLocalMappingMatrix(i, i) = sf_values_slave(integration_point_itr, i)
                    * det_jacobian[integration_point_itr];
                KRATOS_DEBUG_ERROR_IF(sf_values_master(integration_point_itr, i) < 0.0) 
                    << "SHAPE FUNCTIONS LESS THAN ZERO" << std::endl;
            }
        }
    }
    else {
        for (IndexType integration_point_itr = 0; integration_point_itr < sf_values_slave.size1(); ++integration_point_itr) { // This loop is probably redundant - it will always just be 1
            for (IndexType i = 0; i < sf_values_slave.size2(); ++i) {
                for (IndexType j = 0; j < sf_values_master.size2(); ++j) {
                    rLocalMappingMatrix(i, j) = sf_values_slave(integration_point_itr, i)
                        * sf_values_master(integration_point_itr, j)
                        * det_jacobian[integration_point_itr];

                    KRATOS_DEBUG_ERROR_IF(sf_values_master(integration_point_itr, i) < 0.0) 
                        << "SHAPE FUNCTIONS LESS THAN ZERO" << std::endl;
                    KRATOS_DEBUG_ERROR_IF(sf_values_slave(integration_point_itr, j) < 0.0) 
                        << "SHAPE FUNCTIONS LESS THAN ZERO" << std::endl;
                }
            }
        }
    }

    for (IndexType i=0; i< sf_values_master.size2(); ++i) {
        rOriginIds[i] = r_geometry_master[i].GetValue(INTERFACE_EQUATION_ID);
    }
    for (IndexType i=0; i< sf_values_slave.size2(); ++i) {
        rDestinationIds[0] = r_geometry_master[i].GetValue(INTERFACE_EQUATION_ID);
    }
}

std::string CouplingGeometryLocalSystem::PairingInfo(const int EchoLevel) const
{
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

    // assemble projector interface mass matrix - interface_matrix_projector
    const std::size_t num_nodes_interface_slave = mrModelPartDestination.GetSubModelPart("interface").NumberOfNodes(); // fix up and put in the mapper parameters
    const std::size_t num_nodes_interface_master = mrModelPartOrigin.GetSubModelPart("interface").NumberOfNodes();
    Matrix interface_matrix_projector = ZeroMatrix(num_nodes_interface_slave, num_nodes_interface_master);
    for (size_t local_projector_system = 0; 
        local_projector_system < mMapperLocalSystems.size()/2; local_projector_system++)
    {
        // assemble from mMapperLocalSystems(local_projector_system) to interface_matrix_projector
    }

    // assemble slave interface mass matrix - interface_matrix_slave
    Matrix interface_matrix_slave = ZeroMatrix(num_nodes_interface_slave, num_nodes_interface_slave);
    for (size_t local_projector_system = mMapperLocalSystems.size() / 2; 
        local_projector_system < mMapperLocalSystems.size(); local_projector_system++)
    {
        // assemble from mMapperLocalSystems(local_projector_system) to interface_matrix_slave
    }

    // Perform consistency scaling if requested
    if (mMapperSettings["consistency_scaling"].GetBool()) 
        EnforceConsistencyWithScaling(interface_matrix_slave, interface_matrix_projector, 1.1);

    // get total interface mapping matrix
    Matrix inv_interface_matrix_slave(num_nodes_interface_slave, num_nodes_interface_slave);
    double aux_det_slave = 0;
    MathUtils<double>::InvertMatrix(interface_matrix_slave, inv_interface_matrix_slave, aux_det_slave);
    Matrix interface_mapper = prod(interface_matrix_slave, interface_matrix_projector);

    


    // assemble to global
    // directly assemble from interface_mapper to global mapping matrix.
    //BuildMappingMatrix(MappingOptions); // TODO we probably may not need this


}

/* Performs operations that are needed for Initialization and when the interface is updated (All cases)
I.e. Operations that can be performed several times in the livetime of the mapper
*/
template<class TSparseSpace, class TDenseSpace>
void CouplingGeometryMapper<TSparseSpace, TDenseSpace>::BuildMappingMatrix(Kratos::Flags MappingOptions)
{
    AssignInterfaceEquationIds(); // Has to be done ever time in case of overlapping interfaces!

    const int echo_level = mMapperSettings["echo_level"].GetInt();

    MappingMatrixUtilities::BuildMappingMatrix<TSparseSpace, TDenseSpace>(
        mpMappingMatrix,
        mpInterfaceVectorContainerOrigin->pGetVector(),
        mpInterfaceVectorContainerDestination->pGetVector(),
        mpInterfaceVectorContainerOrigin->GetModelPart(),
        mpInterfaceVectorContainerDestination->GetModelPart(),
        mMapperLocalSystems,
        echo_level
    );

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
void CouplingGeometryMapper<TSparseSpace, TDenseSpace>::MapInternal(
    const Variable<array_1d<double, 3>>& rOriginVariable,
    const Variable<array_1d<double, 3>>& rDestinationVariable,
    Kratos::Flags MappingOptions)
{
    const std::vector<std::string> var_comps{"_X", "_Y", "_Z"};

    for (const auto& var_ext : var_comps) {
        const auto& var_origin = KratosComponents<ComponentVariableType>::Get(rOriginVariable.Name() + var_ext);
        const auto& var_destination = KratosComponents<ComponentVariableType>::Get(rDestinationVariable.Name() + var_ext);

        mpInterfaceVectorContainerOrigin->UpdateSystemVectorFromModelPart(var_origin, MappingOptions);

        TSparseSpace::Mult(
            *mpMappingMatrix,
            mpInterfaceVectorContainerOrigin->GetVector(),
            mpInterfaceVectorContainerDestination->GetVector()); // rQd = rMdo * rQo

        mpInterfaceVectorContainerDestination->UpdateModelPartFromSystemVector(var_destination, MappingOptions);
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

        for (IndexType j = 0; j < rInterfaceMatrixSlave.size2(); ++j) {
            row_sum_slave += rInterfaceMatrixSlave(i, j);
            row_sum_projector += rInterfaceMatrixProjected(i, j);
        }

        const double alpha = (row_sum_slave / row_sum_projector < scalingLimit)
            ? row_sum_slave / row_sum_projector
            : scalingLimit;
        for (IndexType j = 0; j < rInterfaceMatrixSlave.size2(); ++j) 
                rInterfaceMatrixProjected(i, j) *= alpha;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class CouplingGeometryMapper< MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType >;

}  // namespace Kratos.
