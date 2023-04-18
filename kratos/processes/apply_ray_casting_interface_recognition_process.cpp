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
    ApplyRayCastingInterfaceRecognitionProcess<TDim>::ApplyRayCastingInterfaceRecognitionProcess(
        ModelPart& rVolumePart,
        ModelPart& rSkinPart,
        Parameters ThisParameters)
        : ApplyRayCastingProcess<TDim>(rVolumePart,
          rSkinPart,
          BaseType::GetDefaultParameters())
    {
        mSettings = ThisParameters;
        mSettings.ValidateAndAssignDefaults(this->GetDefaultParameters());
        mRelativeTolerance = mSettings["relative_tolerance"].GetDouble();
        mpDistanceVariable = &KratosComponents<Variable<double>>::Get(mSettings["distance_variable"].GetString());
        mDistanceGetterFunctor = CreateDistanceGetterFunctor();
    }

    template<std::size_t TDim>
    ApplyRayCastingInterfaceRecognitionProcess<TDim>::ApplyRayCastingInterfaceRecognitionProcess(
        FindIntersectedGeometricalObjectsProcess& TheFindIntersectedObjectsProcess,
        Parameters ThisParameters)
        : ApplyRayCastingProcess<TDim>(TheFindIntersectedObjectsProcess,
          BaseType::GetDefaultParameters())
    {
        mSettings = ThisParameters;
        mSettings.ValidateAndAssignDefaults(this->GetDefaultParameters());
        mRelativeTolerance = mSettings["relative_tolerance"].GetDouble();
        mpDistanceVariable = &KratosComponents<Variable<double>>::Get(mSettings["distance_variable"].GetString());
        mDistanceGetterFunctor = CreateDistanceGetterFunctor();
    }

    template<std::size_t TDim>
    const Parameters ApplyRayCastingInterfaceRecognitionProcess<TDim>::GetDefaultParameters() const
    {
        Parameters default_settings(R"({
            "interface_max_distance" : 1e-6
        })");
        default_settings.RecursivelyAddMissingParameters(BaseType::GetDefaultParameters());
        return default_settings;
    }

    template<std::size_t TDim>
    std::function<void(Node<3>&, const double)> ApplyRayCastingInterfaceRecognitionProcess<TDim>::CreateApplyNodalFunction() const
    {
        const double interface_tol = mSettings["interface_max_distance"].GetDouble();
        return [this,interface_tol](Node<3>& rNode, const double RayDistance) {
            double& r_node_distance = this->mDistanceGetterFunctor(rNode, *(this->mpDistanceVariable));
            if (std::abs(RayDistance) < interface_tol) {
                r_node_distance = 0.0;
            } else if (RayDistance * r_node_distance < 0.0) {
                r_node_distance = -r_node_distance;
            }
        };
    }

    template class Kratos::ApplyRayCastingInterfaceRecognitionProcess<2>;
    template class Kratos::ApplyRayCastingInterfaceRecognitionProcess<3>;
}