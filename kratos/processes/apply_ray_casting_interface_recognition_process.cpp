//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Daniel Diez
//
#include "processes/apply_ray_casting_interface_recognition_process.h"

namespace Kratos
{
    template<std::size_t TDim>
    std::function<void(Node<3>&, const double)> ApplyRayCastingInterfaceRecognitionProcess<TDim>::CreateApplyNodalFunction() const
    {
        return [this](Node<3>& rNode, const double RayDistance) {
            double& r_node_distance = this->mDistanceGetterFunctor(rNode, *(this->mpDistanceVariable));
            if (std::abs(RayDistance) < this->mEpsilon) {
                r_node_distance = 0.0;
            } else if (RayDistance * r_node_distance < 0.0) {
                r_node_distance = -r_node_distance;
            }
        };
    }

    template class Kratos::ApplyRayCastingInterfaceRecognitionProcess<2>;
    template class Kratos::ApplyRayCastingInterfaceRecognitionProcess<3>;
}