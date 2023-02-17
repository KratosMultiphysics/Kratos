//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main author:     Reza Najian Asl
//

// System includes

// Project includes
#include "includes/model_part.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"

// Application includes

// Include base h
#include "simp_utils.h"

namespace Kratos
{

double SimpUtils::ProjectForward(
    const double X,
    const Vector& rYLimits,
    const double Beta,
    const double PenaltyFactor)
{
    const IndexType x1 = std::min<int>(std::max<int>(0, std::floor(X)), rYLimits.size() - 2);
    const IndexType x2 = x1 + 1;

    const double y1 = rYLimits[x1];
    const double y2 = rYLimits[x2];
    const double pow_val = -2.0 * Beta * (X - (x1 + x2) / 2.0);

    return (y2 - y1) / (std::pow(1 + std::exp(pow_val), PenaltyFactor)) + y1;
}

double SimpUtils::ProjectBackward(
    const double Y,
    const Vector& rYLimits,
    const double Beta,
    const double PenaltyFactor)
{
    if (Y >= rYLimits[rYLimits.size() - 1]) {
        return rYLimits.size() - 1;
    } else if (Y <= rYLimits[0]) {
        return  0;
    } else {
        for (IndexType i = 0; i < rYLimits.size() - 1; i++)
            if ((Y > rYLimits[i]) && (Y < rYLimits[i + 1])) {
                const double y1 = rYLimits[i];
                const double y2 = rYLimits[i + 1];
                return ((2.0 * i + 1) / 2.0) + (1.0 / (-2.0 * Beta)) * std::log(std::pow((y2 - y1) / (Y - y1), 1.0 / PenaltyFactor) - 1);
            }
            else if (Y == rYLimits[i]) {
                return i;
            }
            else if (Y == rYLimits[i + 1]) {
                return i + 1;
            }
    }
    return 0;
}

double SimpUtils::ProjectionDerivative(
    const double X,
    const Vector& rYLimits,
    const double Beta,
    const double PenaltyFactor)
{
    const IndexType x1 = std::min<int>(std::max<int>(0, std::floor(X)), rYLimits.size() - 2);
    const IndexType x2 = x1 + 1;

    const double y1 = rYLimits[x1];
    const double y2 = rYLimits[x2];
    const double pow_val = -2.0 * Beta * (X - (x1 + x2) / 2.0);

    const double dydx = (y2 - y1) *
                  (1.0 / std::pow(1 + std::exp(pow_val), PenaltyFactor + 1)) *
                  PenaltyFactor * 2.0 * Beta * std::exp(pow_val);

    return dydx;
}

template<class TContainerType>
void SimpUtils::ProjectContainerVariableDataHolderForward(
    ContainerVariableDataHolderBase<TContainerType>& rY,
    const ContainerVariableDataHolderBase<TContainerType>& rX,
    const Vector& rYLimits,
    const double Beta,
    const double PenaltyFactor)
{
    rY.ResizeDataTo(rX);
    IndexPartition<IndexType>(rX.GetData().size()).for_each([&](const IndexType Index) {
        rY.GetData()[Index] = ProjectForward(rX.GetData()[Index], rYLimits, Beta, PenaltyFactor);
    });
}

template<class TContainerType>
void SimpUtils::ProjectContainerVariableDataHolderBackward(
    ContainerVariableDataHolderBase<TContainerType>& rX,
    const ContainerVariableDataHolderBase<TContainerType>& rY,
    const Vector& rYLimits,
    const double Beta,
    const double PenaltyFactor)
{
    rX.ResizeDataTo(rY);
    IndexPartition<IndexType>(rY.GetData().size()).for_each([&](const IndexType Index) {
        rX.GetData()[Index] = ProjectBackward(rY.GetData()[Index], rYLimits, Beta, PenaltyFactor);
    });
}

template<class TContainerType>
void SimpUtils::ProjectContainerVariableDataHolderDerivative(
    ContainerVariableDataHolderBase<TContainerType>& rDerivative,
    const ContainerVariableDataHolderBase<TContainerType>& rX,
    const Vector& rYLimits,
    const double Beta,
    const double PenaltyFactor)
{
    rDerivative.ResizeDataTo(rX);
    IndexPartition<IndexType>(rX.GetData().size()).for_each([&](const IndexType Index) {
        rDerivative.GetData()[Index] = ProjectionDerivative(rX.GetData()[Index], rYLimits, Beta, PenaltyFactor);
    });
}

// template instantiations
#define HELMHOLTZ_UTILITIES_CONTAINER_INSTANTIATION(ContainerType)                                                                                                                                                               \
    template void SimpUtils::ProjectContainerVariableDataHolderForward(ContainerVariableDataHolderBase<ContainerType>&, const ContainerVariableDataHolderBase<ContainerType>& rX,  const Vector&, const double, const double);   \
    template void SimpUtils::ProjectContainerVariableDataHolderBackward(ContainerVariableDataHolderBase<ContainerType>&, const ContainerVariableDataHolderBase<ContainerType>& rX,  const Vector&, const double, const double);  \
    template void SimpUtils::ProjectContainerVariableDataHolderDerivative(ContainerVariableDataHolderBase<ContainerType>&, const ContainerVariableDataHolderBase<ContainerType>& rX,  const Vector&, const double, const double);

HELMHOLTZ_UTILITIES_CONTAINER_INSTANTIATION(ModelPart::NodesContainerType)
HELMHOLTZ_UTILITIES_CONTAINER_INSTANTIATION(ModelPart::ConditionsContainerType)
HELMHOLTZ_UTILITIES_CONTAINER_INSTANTIATION(ModelPart::ElementsContainerType)

#undef HELMHOLTZ_UTILITIES_CONTAINER_INSTANTIATION

}