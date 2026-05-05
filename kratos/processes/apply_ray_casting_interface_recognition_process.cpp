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
        Model& rModel,
        Parameters ThisParameters)
        : ApplyRayCastingProcess<TDim>(
            rModel.GetModelPart(ThisParameters["volume_model_part"].GetString()),
            rModel.GetModelPart(ThisParameters["skin_model_part"].GetString()),
            BaseType::GetDefaultParameters())
    {
        this->mSettings = ThisParameters;
        this->mSettings.ValidateAndAssignDefaults(this->GetDefaultParameters());
        this->mRelativeTolerance = this->mSettings["relative_tolerance"].GetDouble();
        this->mpDistanceVariable = &KratosComponents<Variable<double>>::Get(this->mSettings["distance_variable"].GetString());
        this->mDistanceGetterFunctor = this->CreateDistanceGetterFunctor();
    }

    template<std::size_t TDim>
    ApplyRayCastingInterfaceRecognitionProcess<TDim>::ApplyRayCastingInterfaceRecognitionProcess(
        FindIntersectedGeometricalObjectsProcess& TheFindIntersectedObjectsProcess,
        Parameters ThisParameters)
        : ApplyRayCastingProcess<TDim>(TheFindIntersectedObjectsProcess,
          BaseType::GetDefaultParameters())
    {
        this->mSettings = ThisParameters;
        this->mSettings.ValidateAndAssignDefaults(this->GetDefaultParameters());
        this->mRelativeTolerance = this->mSettings["relative_tolerance"].GetDouble();
        this->mpDistanceVariable = &KratosComponents<Variable<double>>::Get(this->mSettings["distance_variable"].GetString());
        this->mDistanceGetterFunctor = this->CreateDistanceGetterFunctor();
    }

    template<std::size_t TDim>
    Process::Pointer ApplyRayCastingInterfaceRecognitionProcess<TDim>::Create(
        Model& rModel,
        Parameters ThisParameters)
    {
        return Kratos::make_shared<ApplyRayCastingInterfaceRecognitionProcess<TDim>>(rModel, ThisParameters);
    }

    template<std::size_t TDim>
    const Parameters ApplyRayCastingInterfaceRecognitionProcess<TDim>::GetDefaultParameters() const
    {
        Parameters default_settings(R"({
            "volume_model_part" : "",
            "skin_model_part" : "",
            "interface_max_distance" : 1e-6
        })");
        default_settings.RecursivelyAddMissingParameters(BaseType::GetDefaultParameters());
        return default_settings;
    }

    template<std::size_t TDim>
    std::function<void(Node&, const double)> ApplyRayCastingInterfaceRecognitionProcess<TDim>::CreateApplyNodalFunction() const
    {
        const double interface_tol = this->mSettings["interface_max_distance"].GetDouble();
        return [this,interface_tol](Node& rNode, const double RayDistance) {
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