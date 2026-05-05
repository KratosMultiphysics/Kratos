//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:    Juan Ignacio Camarotti


// System includes

// External includes

// Project includes
#include "nearest_neighbor_mapper_iga.h"
#include "mapping_application_variables.h"

namespace Kratos {

void NearestNeighborInterfaceInfoIGA::ProcessSearchResult(const InterfaceObject& rInterfaceObject)
{
    SetLocalSearchWasSuccessful();

    const double neighbor_distance = MapperUtilities::ComputeDistance(this->Coordinates(), rInterfaceObject.Coordinates());

    const auto& shape_function_values = rInterfaceObject.pGetBaseGeometry()->ShapeFunctionsValues();
    const std::size_t num_values = shape_function_values.size2();

    if (neighbor_distance < mNearestNeighborDistance) {
        mNearestNeighborDistance = neighbor_distance;
        mNumberOfNearestNeighbors = 1;

        mShapeFunctionValues.clear();
        mNearestNeighborId.clear();

        // Reserve space for the vectors
        mShapeFunctionValues.reserve(mShapeFunctionValues.size() + num_values);
        mNearestNeighborId.reserve(mNearestNeighborId.size() + num_values);

        for (std::size_t i = 0; i < num_values; ++i) {
            mShapeFunctionValues.push_back(shape_function_values(i));
            mNearestNeighborId.push_back(rInterfaceObject.pGetBaseGeometry()->Points()[i].GetValue(INTERFACE_EQUATION_ID));
        }

    } 
    else if (neighbor_distance == mNearestNeighborDistance) {
        mNumberOfNearestNeighbors += 1;

        // Reserve additional space for the vectors
        mShapeFunctionValues.reserve(mShapeFunctionValues.size() + num_values); 
        mNearestNeighborId.reserve(mNearestNeighborId.size() + num_values);    

        for (std::size_t i = 0; i < num_values; ++i) {
            mShapeFunctionValues.push_back(shape_function_values(i));
            mNearestNeighborId.push_back(rInterfaceObject.pGetBaseGeometry()->Points()[i].GetValue(INTERFACE_EQUATION_ID));
        }
    }
}

void NearestNeighborLocalSystemIGA::CalculateAll(MatrixType& rLocalMappingMatrix,
                    EquationIdVectorType& rOriginIds,
                    EquationIdVectorType& rDestinationIds,
                    MapperLocalSystem::PairingStatus& rPairingStatus) const
{
    if (!mInterfaceInfos.empty()) {  
        rPairingStatus = MapperLocalSystem::PairingStatus::InterfaceInfoFound;

        rDestinationIds.resize(1);

        std::vector<int> nearest_neighbor_id;
        IndexType number_of_nearest_neighbors = 0;

        mInterfaceInfos.front()->GetValue(nearest_neighbor_id, MapperInterfaceInfo::InfoType::Dummy);
        mInterfaceInfos.front()->GetValue(number_of_nearest_neighbors, MapperInterfaceInfo::InfoType::Dummy);
        rOriginIds = nearest_neighbor_id;

        std::vector<double> sf_values;
        mInterfaceInfos.front()->GetValue(sf_values, MapperInterfaceInfo::InfoType::Dummy);

        rLocalMappingMatrix.resize(1, rOriginIds.size(), false);

        KRATOS_DEBUG_ERROR_IF(number_of_nearest_neighbors == 0) << "No nearest neighbors found for the node: " << mpNode->Id() << std::endl;
        const double mapping_weight = 1.0/number_of_nearest_neighbors;
        
        for (IndexType i=0; i<rOriginIds.size(); ++i) {
            rLocalMappingMatrix(0,i) = mapping_weight * sf_values[i]; 
        }

        KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not initialized!" << std::endl;
        rDestinationIds[0] = mpNode->GetValue(INTERFACE_EQUATION_ID);
    } else {
        ResizeToZero(rLocalMappingMatrix, rOriginIds, rDestinationIds, rPairingStatus);
    }
}

void NearestNeighborLocalSystemIGA::PairingInfo(std::ostream& rOStream, const int EchoLevel) const
{
    KRATOS_DEBUG_ERROR_IF_NOT(mpNode) << "Members are not initialized!" << std::endl;

    rOStream << "NearestNeighborLocalSystem based on " << mpNode->Info();
    if (EchoLevel > 3) {
        rOStream << " at Coordinates " << Coordinates()[0] << " | " << Coordinates()[1] << " | " << Coordinates()[2];
    }
}

void NearestNeighborLocalSystemIGA::SetPairingStatusForPrinting()
{
    if (mPairingStatus == MapperLocalSystem::PairingStatus::Approximation) {
        mpNode->SetValue(PAIRING_STATUS, 0);
    } else {
        mpNode->SetValue(PAIRING_STATUS, -1);
    }
}

}  // namespace Kratos.
