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
    const auto& r_origin_model_part = mpModel->GetModelPart(origin_model_part_name);
    auto& r_destination_model_part = mpModel->GetModelPart(destination_model_part_name);

    // clear the properties of the destination model part
    std::vector<std::size_t> properties_ids;
    for (const auto& r_prop : r_origin_model_part.rProperties()) {
        properties_ids.push_back(r_prop.Id());
    }
    for (auto prop_id : properties_ids) {
        r_destination_model_part.RemovePropertiesFromAllLevels(prop_id);
    }

    // make copies of the properties
    for (auto& r_prop : r_origin_model_part.rProperties()) {
        r_destination_model_part.AddProperties(Kratos::make_shared<Properties>(r_prop));
    }

    // replace the properties of the elements and conditions
    ReplaceProperties(r_destination_model_part.Elements(), r_destination_model_part);
    ReplaceProperties(r_destination_model_part.Conditions(), r_destination_model_part);
}

template<class TContainerType>
void CopyPropertiesModeler::ReplaceProperties(
    TContainerType& rContainer,
    ModelPart& rModelPart)
{
    block_for_each(rContainer, [&](typename TContainerType::value_type& rEntity){
        auto properties_id = rEntity.GetProperties().Id();
        auto p_properties = rModelPart.pGetProperties(properties_id);
        rEntity.SetProperties(p_properties);
    });
}

}
