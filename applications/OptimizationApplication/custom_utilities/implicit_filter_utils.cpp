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
#include <cmath>

// Project includes
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/variable_utils.h"

// Application includes
#include "optimization_application_variables.h"

// Include base h
#include "implicit_filter_utils.h"

namespace Kratos
{

void ImplicitFilterUtils::CalculateNodeNeighbourCount(ModelPart& rModelPart)
{
    // Calculate number of neighbour elements for each node.
    VariableUtils().SetNonHistoricalVariableToZero(NUMBER_OF_NEIGHBOUR_ELEMENTS, rModelPart.Nodes());
    block_for_each(rModelPart.Elements(), [&](ModelPart::ElementType& rElement) {
        auto& r_geometry = rElement.GetGeometry();
        for (unsigned j = 0; j < r_geometry.PointsNumber(); ++j) {
            double& r_num_neighbour =
                r_geometry[j].GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
            AtomicAdd(r_num_neighbour, 1.0);
        }
    });
    rModelPart.GetCommunicator().AssembleNonHistoricalData(NUMBER_OF_NEIGHBOUR_ELEMENTS);
}

void ImplicitFilterUtils::SetBulkRadiusForShapeFiltering(ModelPart& rModelPart)
{
    ProcessInfo &rCurrentProcessInfo = rModelPart.GetProcessInfo();
    rCurrentProcessInfo.SetValue(HELMHOLTZ_BULK_RADIUS_SHAPE,1.0);

    const double local_volume_strain_energy = block_for_each<SumReduction<double>>(rModelPart.Elements(), [&](auto& rElement) {
        double elem_val;
        rElement.Calculate(ELEMENT_STRAIN_ENERGY,elem_val,rCurrentProcessInfo);
        return elem_val;
    });

    const double local_surface_strain_energy = block_for_each<SumReduction<double>>(rModelPart.Conditions(), [&](auto& rCondition) {
        double cond_val;
        rCondition.Calculate(ELEMENT_STRAIN_ENERGY,cond_val,rCurrentProcessInfo);
        return cond_val;
    });

    double bulk_filter_size = rModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_surface_strain_energy);
    bulk_filter_size /= rModelPart.GetCommunicator().GetDataCommunicator().SumAll(local_volume_strain_energy);

    rCurrentProcessInfo.SetValue(HELMHOLTZ_BULK_RADIUS_SHAPE, bulk_filter_size);
}

void ImplicitFilterUtils::AssignProperties(
    ModelPart& rModelPart,
    Parameters PropertiesParams)
{
    KRATOS_TRY

    const auto& defaults = Parameters(R"(
    {
        "properties_id": 1,
        "Material": {}
    })" );

    PropertiesParams.ValidateAndAssignDefaults(defaults);

    auto p_properties = rModelPart.CreateNewProperties(PropertiesParams["properties_id"].GetInt());

    auto material_properties = PropertiesParams["Material"];

    if (material_properties.Has("constitutive_law")) {
        KRATOS_ERROR_IF_NOT(material_properties["constitutive_law"].Has("name"))
            << "Material constitutive law does not define the name.\n";

        auto p_constitutive_law = KratosComponents<ConstitutiveLaw>::Get(material_properties["constitutive_law"]["name"].GetString()).Clone();
        p_properties->SetValue(CONSTITUTIVE_LAW, p_constitutive_law);
    }

    if (material_properties.Has("Variables")) {
        const auto variables = material_properties["Variables"];
        for (auto itr = variables.begin(); itr != variables.end(); ++itr) {
            const auto& r_variable_name = itr.name();
            const auto variable_value_param = *itr;

            if (KratosComponents<Variable<int>>::Has(r_variable_name)) {
                p_properties->SetValue(KratosComponents<Variable<int>>::Get(r_variable_name), variable_value_param.GetInt());
            } else if (KratosComponents<Variable<double>>::Has(r_variable_name)) {
                p_properties->SetValue(KratosComponents<Variable<double>>::Get(r_variable_name), variable_value_param.GetDouble());
            } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(r_variable_name)) {
                p_properties->SetValue(KratosComponents<Variable<array_1d<double, 3>>>::Get(r_variable_name), variable_value_param.GetVector());
            } else if (KratosComponents<Variable<array_1d<double, 4>>>::Has(r_variable_name)) {
                p_properties->SetValue(KratosComponents<Variable<array_1d<double, 4>>>::Get(r_variable_name), variable_value_param.GetVector());
            } else if (KratosComponents<Variable<array_1d<double, 6>>>::Has(r_variable_name)) {
                p_properties->SetValue(KratosComponents<Variable<array_1d<double, 6>>>::Get(r_variable_name), variable_value_param.GetVector());
            } else if (KratosComponents<Variable<array_1d<double, 9>>>::Has(r_variable_name)) {
                p_properties->SetValue(KratosComponents<Variable<array_1d<double, 9>>>::Get(r_variable_name), variable_value_param.GetVector());
            }
        }
    }

    block_for_each(rModelPart.Elements(), [&p_properties](auto& rElement) {
        rElement.SetProperties(p_properties);
    });

    block_for_each(rModelPart.Conditions(), [&p_properties](auto& rCondition) {
        rCondition.SetProperties(p_properties);
    });

    KRATOS_CATCH("");
}

}
