//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         DigitalTwinApplication/license.txt
//
//  Main authors:    Ihar Antonau
//


// System includes
#include <iterator>

// External includes

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/atomic_utilities.h"
#include "expression/literal_flat_expression.h"


// Application includes
#include "digital_twin_application_variables.h"

// Include base h
#include "sensor_isolation_response_utils.h"

namespace Kratos {

SensorIsolationResponseUtils::SensorIsolationResponseUtils(
    ModelPart& rSensorModelPart,
    const IndexType MaxNumberOfNeighbors,
    const double Radius)
    : mpSensorModelPart(&rSensorModelPart),
      mMaxNumberOfNeighbors(MaxNumberOfNeighbors),
      mRadius(Radius)
{
}

void SensorIsolationResponseUtils::Initialize()
{
    KRATOS_TRY

    mpSearchTree =  Kratos::make_shared<SensorIsolationResponseUtils::KDTree>(mpSensorModelPart->Nodes().ptr_begin(), mpSensorModelPart->Nodes().ptr_end(), 100);    

    KRATOS_CATCH("")
}

double SensorIsolationResponseUtils::CalculateValue()
{
    struct KDTreeThreadLocalStorage
    {
        explicit KDTreeThreadLocalStorage(const IndexType MaxNumberOfNeighbors)
        {
            mNeighbourEntityPoints.resize(MaxNumberOfNeighbors);
            mResultingSquaredDistances.resize(MaxNumberOfNeighbors);
        }

        std::vector<ModelPart::NodeType::Pointer> mNeighbourEntityPoints;
        std::vector<double> mResultingSquaredDistances;
    };    

    return block_for_each<SumReduction<double>>(mpSensorModelPart->Nodes(), KDTreeThreadLocalStorage(mMaxNumberOfNeighbors), [&](const auto& rNode, KDTreeThreadLocalStorage& rTLS){
        const auto number_of_neighbors = this->mpSearchTree->SearchInRadius(
                                            rNode,
                                            this->mRadius,
                                            rTLS.mNeighbourEntityPoints.begin(),
                                            rTLS.mResultingSquaredDistances.begin(),
                                            mMaxNumberOfNeighbors);

        KRATOS_ERROR_IF(number_of_neighbors >= mMaxNumberOfNeighbors)
            << "Maximum number of allowed neighbours reached when searching for neighbours in "
            << this->mpSensorModelPart->FullName() << " with radii = " << this->mRadius << " [ max number of allowed neighbours = "
            << mMaxNumberOfNeighbors << " ].\n";

        double result = 0.0;
        if (rNode.GetValue(SENSOR_STATUS) >= 0.999) {
            for (IndexType i = 1; i < number_of_neighbors; ++i) {
                result += rTLS.mNeighbourEntityPoints[i]->GetValue(SENSOR_STATUS);
            }
        }
        return result;
    });    
}

ContainerExpression<ModelPart::NodesContainerType> SensorIsolationResponseUtils::CalculateGradient()
{
    KRATOS_TRY

    struct KDTreeThreadLocalStorage
    {
        explicit KDTreeThreadLocalStorage(const IndexType MaxNumberOfNeighbors)
        {
            mNeighbourEntityPoints.resize(MaxNumberOfNeighbors);
            mResultingSquaredDistances.resize(MaxNumberOfNeighbors);
        }

        std::vector<ModelPart::NodeType::Pointer> mNeighbourEntityPoints;
        std::vector<double> mResultingSquaredDistances;
    };    

    auto p_expression = LiteralFlatExpression<double>::Create(mpSensorModelPart->NumberOfNodes(), {});
    IndexPartition<IndexType>(mpSensorModelPart->NumberOfNodes()).for_each([&p_expression](const auto Index) {
        *(p_expression->begin() + Index) = 0.0;
    });

    block_for_each(mpSensorModelPart->Nodes(), KDTreeThreadLocalStorage(mMaxNumberOfNeighbors), [&](const auto& rNode, KDTreeThreadLocalStorage& rTLS){
        const auto number_of_neighbors = this->mpSearchTree->SearchInRadius(
                                            rNode,
                                            this->mRadius,
                                            rTLS.mNeighbourEntityPoints.begin(),
                                            rTLS.mResultingSquaredDistances.begin(),
                                            mMaxNumberOfNeighbors);

        KRATOS_ERROR_IF(number_of_neighbors >= mMaxNumberOfNeighbors)
            << "Maximum number of allowed neighbours reached when searching for neighbours in "
            << this->mpSensorModelPart->FullName() << " with radii = " << this->mRadius << " [ max number of allowed neighbours = "
            << mMaxNumberOfNeighbors << " ].\n";

        if (rNode.GetValue(SENSOR_STATUS) >= 0.999) {
            for (IndexType i = 1; i < number_of_neighbors; ++i) {
                AtomicAdd<double>(*(p_expression->begin() + rTLS.mNeighbourEntityPoints[i]->Id() - 1), 1.0);
            }
        }
    });  

    ContainerExpression<ModelPart::NodesContainerType> result(*mpSensorModelPart);
    result.SetExpression(p_expression);
    return result;  

    KRATOS_CATCH("");
}

///@} // Kratos Classes

} /* namespace Kratos.*/