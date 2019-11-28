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
    const auto& r_geometry_origin = mpGeom->GetGeometryPart(0);
    const auto& r_geometry_destination = mpGeom->GetGeometryPart(1);

    const std::size_t number_of_nodes_origin = r_geometry_origin.size();
    const std::size_t number_of_nodes_destination = r_geometry_destination.size();

    const std::size_t mat_size = (number_of_nodes_origin + number_of_nodes_destination);

    rPairingStatus = MapperLocalSystem::PairingStatus::InterfaceInfoFound;

    if (rLocalMappingMatrix.size1() != mat_size || rLocalMappingMatrix.size2() != mat_size) {
        rLocalMappingMatrix.resize(mat_size, mat_size, false);
    }
    if (rOriginIds.size()      != number_of_nodes_origin) rOriginIds.resize(number_of_nodes_origin);
    if (rDestinationIds.size() != number_of_nodes_destination) rDestinationIds.resize(number_of_nodes_destination);

    auto sf_values_origin = r_geometry_origin.ShapeFunctionsValues();
    auto integration_point_values_origin = r_geometry_origin.IntegrationPoints();
    Vector det_jacobian;
    auto determinant_of_jacobian_values_origin = r_geometry_origin.DeterminantOfJacobian(det_jacobian);

    auto sf_values_destination = r_geometry_destination.ShapeFunctionsValues();
    //auto integration_point_values_destination = r_geometry_destination.IntegrationPoints();
    //auto determinant_of_jacobian_values_destination = r_geometry_destination.DeterminantOfJacobian();

    for (IndexType integration_point_itr=0; integration_point_itr<sf_values_origin.size1(); ++integration_point_itr){
        for (IndexType i=0; i<sf_values_origin.size2(); ++i) {
            for (IndexType j=0; j<sf_values_destination.size2(); ++j) {
                rLocalMappingMatrix(i,j) = sf_values_origin( integration_point_itr, i )
                 * sf_values_destination( integration_point_itr, j )
                 * integration_point_values_origin[integration_point_itr].Weight()
                 * determinant_of_jacobian_values_origin[i];
            }
        }
    }

    for (IndexType i=0; i<sf_values_origin.size2(); ++i) {
        rOriginIds[i] = r_geometry_origin[i].GetValue(INTERFACE_EQUATION_ID);
    }
    for (IndexType i=0; i<sf_values_destination.size2(); ++i) {
        rDestinationIds[0] = r_geometry_origin[i].GetValue(INTERFACE_EQUATION_ID);
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
void CouplingGeomteryMapper<TSparseSpace, TDenseSpace>::InitializeInterface(Kratos::Flags MappingOptions)
{
    // @tteschemachen here kann man theoretisch auch das Origin-MP nehmen
    // kommt drauf an, wo die Coupling-Geometries sind
    CreateMapperLocalSystems(mrModelPartDestination.GetCommunicator(),
                             mMapperLocalSystems);

    BuildMappingMatrix(MappingOptions);
}

/* Performs operations that are needed for Initialization and when the interface is updated (All cases)
I.e. Operations that can be performed several times in the livetime of the mapper
*/
template<class TSparseSpace, class TDenseSpace>
void CouplingGeomteryMapper<TSparseSpace, TDenseSpace>::BuildMappingMatrix(Kratos::Flags MappingOptions)
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
void CouplingGeomteryMapper<TSparseSpace, TDenseSpace>::MapInternal(
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
void CouplingGeomteryMapper<TSparseSpace, TDenseSpace>::MapInternal(
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

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class CouplingGeomteryMapper< MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType >;

}  // namespace Kratos.
