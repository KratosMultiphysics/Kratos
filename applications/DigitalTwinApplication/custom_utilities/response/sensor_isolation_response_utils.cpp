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
#include <limits>

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
    const double Radius,
    const double Beta)
    : mpSensorModelPart(&rSensorModelPart),
      mMaxNumberOfNeighbors(MaxNumberOfNeighbors),
      mRadius(Radius),
      mBeta(Beta)
{
}

void SensorIsolationResponseUtils::Initialize()
{
    KRATOS_TRY

    mNodes.resize(mpSensorModelPart->NumberOfNodes());
    IndexPartition<IndexType>(mpSensorModelPart->NumberOfNodes()).for_each([&](const auto Index) {
        mNodes[Index] = *(mpSensorModelPart->Nodes().ptr_begin() + Index);
    });

    mpSearchTree =  Kratos::make_shared<SensorIsolationResponseUtils::KDTree>(mNodes.begin(), mNodes.end(), 100);

    mClusterStatus.resize(mNodes.size());

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

    std::tie(mNumerator, mDenominator) = block_for_each<CombinedReduction<SumReduction<double>, SumReduction<double>>>(mpSensorModelPart->Nodes(), KDTreeThreadLocalStorage(mMaxNumberOfNeighbors), [&](const auto& rNode, KDTreeThreadLocalStorage& rTLS){
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

        double& value = mClusterStatus[rNode.Id() - 1];
        value = 0.0;

        for (IndexType i = 1; i < number_of_neighbors; ++i) {
            value += rTLS.mNeighbourEntityPoints[i]->GetValue(SENSOR_STATUS);
        }

        value /= number_of_neighbors;
        return std::make_tuple(value * std::exp(mBeta * value), std::exp(mBeta * value));
    });

    if (mDenominator > std::numeric_limits<double>::epsilon()) {
        return mNumerator / mDenominator;
    } else {
        return 0.0;
    }
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

    if (mDenominator > std::numeric_limits<double>::epsilon()) {
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

            for (IndexType i = 1; i < number_of_neighbors; ++i) {
                const IndexType node_index = rTLS.mNeighbourEntityPoints[i]->Id() - 1;
                const double cluster_status = mClusterStatus[node_index];
                const double value = std::exp(mBeta * cluster_status) * (1 + mBeta * cluster_status - mNumerator * mBeta / mDenominator) / (mDenominator * number_of_neighbors);
                AtomicAdd<double>(*(p_expression->begin() + rTLS.mNeighbourEntityPoints[i]->Id() - 1), value);
            }
        });
    }

    ContainerExpression<ModelPart::NodesContainerType> result(*mpSensorModelPart);
    result.SetExpression(p_expression);
    return result;

    KRATOS_CATCH("");
}

///@} // Kratos Classes

} /* namespace Kratos.*/