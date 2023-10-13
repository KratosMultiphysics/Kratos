//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//
//

// System includes

// External includes

// Project includes
#include "modeler/copy_properties_modeler.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{

CopyPropertiesModeler::CopyPropertiesModeler(
    Model& rModel,
    Parameters ModelerParameters)
        : Modeler(rModel, ModelerParameters)
        , mpModel(&rModel)
{
    mParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
}

CopyPropertiesModeler::CopyPropertiesModeler(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart)
        : Modeler()
{
    mpModel = &rOriginModelPart.GetModel();
    KRATOS_ERROR_IF_NOT(mpModel == &rDestinationModelPart.GetModel()) << "CopyPropertiesModeler. The model parts belong to different models." << std::endl;
    mParameters.AddString("origin_model_part_name", rOriginModelPart.Name());
    mParameters.AddString("destination_model_part_name", rDestinationModelPart.Name());
}

Modeler::Pointer CopyPropertiesModeler::Create(Model& rModel,
    const Parameters ModelParameters) const
{
    return Kratos::make_shared<CopyPropertiesModeler>(rModel, ModelParameters);
}

const Parameters CopyPropertiesModeler::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"({
        "echo_level"                  : 0,
        "origin_model_part_name"      : "",
        "destination_model_part_name" : ""
    })");
    return default_parameters;
}

void CopyPropertiesModeler::SetupModelPart()
{
    const std::string origin_model_part_name = mParameters["origin_model_part_name"].GetString();
    const std::string destination_model_part_name = mParameters["destination_model_part_name"].GetString();
    auto& r_origin_model_part = mpModel->GetModelPart(origin_model_part_name);
    auto& r_destination_model_part = mpModel->GetModelPart(destination_model_part_name);

    // clear the properties of the destination model part
    r_destination_model_part.SetProperties(Kratos::make_shared<ModelPart::PropertiesContainerType>());

    // make copies of the properties
    RecursivelyCopyProperties(r_origin_model_part, r_destination_model_part);

    // replace the properties of the elements and conditions
    ReplaceProperties(r_destination_model_part.Elements(), r_destination_model_part);
    ReplaceProperties(r_destination_model_part.Conditions(), r_destination_model_part);
}

void CopyPropertiesModeler::RecursivelyCopyProperties(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart)
{
    for (auto& r_prop : rOriginModelPart.rProperties()) {
        rDestinationModelPart.AddProperties(Kratos::make_shared<Properties>(r_prop));
    }
    for (auto& r_orig_sub_mp : rOriginModelPart.SubModelParts()) {
        if (rDestinationModelPart.HasSubModelPart(r_orig_sub_mp.Name())) {
            auto& r_dest_sub_mp = rDestinationModelPart.GetSubModelPart(r_orig_sub_mp.Name());
            RecursivelyCopyProperties(r_orig_sub_mp, r_dest_sub_mp);
        }
    }
}

template<class TContainerType>
void CopyPropertiesModeler::ReplaceProperties(
    TContainerType& rContainer,
    const ModelPart& rModelPart)
{
    block_for_each(rContainer, [&](typename TContainerType::value_type& rEntity){
        auto properties_id = rEntity.GetProperties().Id();
        auto p_properties = rModelPart.pGetProperties(properties_id);
        rEntity.SetProperties(p_properties);
    });
}

}
