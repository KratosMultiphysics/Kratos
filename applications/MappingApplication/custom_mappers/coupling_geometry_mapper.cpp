//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "coupling_geometry_mapper.h"
#include "mapping_application_variables.h"

namespace Kratos
{
void CouplingGeometryLocalSystem::CalculateAll(MatrixType& rLocalMappingMatrix,
                    EquationIdVectorType& rOriginIds,
                    EquationIdVectorType& rDestinationIds,
                    MapperLocalSystem::PairingStatus& rPairingStatus) const
{
    const auto& r_geometry_origin = mpGeom->GetGeometryPart(0);
    const auto& r_geometry_destination = mpGeom->GetGeometryPart(1);

    const int number_of_nodes_origin = r_geometry_origin.size();
    const int number_of_nodes_destination = r_geometry_destination.size();

    const int mat_size = (number_of_nodes_origin + number_of_nodes_destination);

    rPairingStatus = MapperLocalSystem::PairingStatus::InterfaceInfoFound;

    if (rLocalMappingMatrix.size1() != mat_size || rLocalMappingMatrix.size2() != mat_size) {
        rLocalMappingMatrix.resize(mat_size, mat_size, false);
    }
    if (rOriginIds.size()      != number_of_nodes_origin) rOriginIds.resize(number_of_nodes_origin);
    if (rDestinationIds.size() != r_geometry_destination) rDestinationIds.resize(r_geometry_destination);

    auto sf_values_origin = r_geometry_origin.ShapeFunctionValues();
    auto integration_point_values_origin = r_geometry_origin.IntegrationPoints();
    auto determinant_of_jacobian_values_origin = r_geometry_origin.DeterminatOfJacobian();

    auto sf_values_destination = r_geometry_destination.ShapeFunctionValues();
    //auto integration_point_values_destination = r_geometry_destination.IntegrationPoints();
    //auto determinant_of_jacobian_values_destination = r_geometry_destination.DeterminatOfJacobian();

    for (IndexType integration_point_itr=0; integration_point_itr<sf_values_origin.size1(); ++integration_point_itr)
    {
        for (IndexType i=0; i<sf_values_origin.size2(); ++i)
        {
            for (IndexType j=0; j<sf_values_destination.size2(); ++j)
            {
                rLocalMappingMatrix(i,j) = sf_values_origin( integration_point_itr, i )
                 * sf_values_destination( integration_point_itr, j )
                 * integration_point_values_origin[integration_point_itr].Weight()
                 * determinant_of_jacobian_values_origin[i];
            }
        }
    }

    for (IndexType i=0; i<sf_values_origin.size2(); ++i)
    {
        rOriginIds[i] = r_geometry_origin[i]->GetValue(INTERFACE_EQUATION_ID);
    }
    for (IndexType i=0; i<sf_values_destination.size2(); ++i)
    {
        rDestinationIds[0] = sf_values_destination[i]->GetValue(INTERFACE_EQUATION_ID);
    }
}

std::string CouplingGeometryLocalSystem::PairingInfo(const int EchoLevel) const
{
    KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not intitialized!" << std::endl;

    std::stringstream buffer;
    buffer << "CouplingGeometryLocalSystem based on " << mpNode->Info();
    if (EchoLevel > 1) { // TODO leave here?
        buffer << " at Coodinates " << Coordinates()[0] << " | " << Coordinates()[1] << " | " << Coordinates()[2];
        if (mPairingStatus == MapperLocalSystem::PairingStatus::Approximation) {
            mpNode->SetValue(PAIRING_STATUS, 0);
        } else {
            mpNode->SetValue(PAIRING_STATUS, -1);
        }
    }
    return buffer.str();
}

}  // namespace Kratos.
