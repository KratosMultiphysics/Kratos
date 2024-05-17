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

#pragma once

// System includes
#include <vector>
#include <string>

// External includes

// Project includes
#include "expression/container_expression.h"
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "spatial_containers/spatial_containers.h"

// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

class SensorIsolationResponseUtils
{
public:
    ///@name Type definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(SensorIsolationResponseUtils);

    // Type definitions for tree-search
    using BucketType = Bucket<3, ModelPart::NodeType, std::vector<ModelPart::NodeType::Pointer>, ModelPart::NodeType::Pointer>;

    using KDTree = Tree<KDTreePartition<BucketType>>;

    ///@}
    ///@name Life cycle
    ///@{

    SensorIsolationResponseUtils(
        ModelPart& rSensorModelPart,
        const IndexType MaxNumberOfNeighbors,
        const double Radius,
        const double Beta);

    ///@}
    ///@name Public operations
    ///@{

    void Initialize();

    double CalculateValue();

    ContainerExpression<ModelPart::NodesContainerType> CalculateGradient();

    ///@}

private:
    ///@name Private member variables
    ///@{

    ModelPart* mpSensorModelPart;

    const IndexType mMaxNumberOfNeighbors;

    const double mRadius;

    const double mBeta;

    double mNumerator;

    double mDenominator;

    std::vector<ModelPart::NodeType::Pointer> mNodes;

    KDTree::Pointer mpSearchTree;

    std::vector<double> mClusterStatus;

    ///@}
};

///@} // Kratos Classes

} /* namespace Kratos.*/