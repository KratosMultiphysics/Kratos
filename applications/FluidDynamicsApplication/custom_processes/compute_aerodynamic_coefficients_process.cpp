//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Eduard Gómez
//
//

// System includes
#include <algorithm>
#include <sstream>
#include <map>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "includes/cfd_variables.h"
#include "includes/condition.h"
#include "includes/properties.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "compute_aerodynamic_coefficients_process.h"
#include "fluid_dynamics_application_variables.h"
#include "utilities/reduction_utilities.h"
#include "utilities/variable_utils.h"

namespace Kratos
{

ComputeAerodynamicCoefficientsProcess::ComputeAerodynamicCoefficientsProcess(
    ModelPart& rModelPart,
    Parameters Params)
    : Process(),
      mrModelPart(rModelPart)
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(Params);
}

ComputeAerodynamicCoefficientsProcess::ComputeAerodynamicCoefficientsProcess(
    Model &rModel,
    Parameters Params)
    : Process(),
      mrModelPart(rModel.GetModelPart(Params["model_part_name"].GetString()))
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(Params);
}


void ComputeAerodynamicCoefficientsProcess::CheckDefaultsAndProcessSettings(Parameters Params)
{
    KRATOS_TRY

    Parameters default_parameters( R"(
    {
        "model_part_name"     : "",
        "freestream_pressure" : 0,
        "density_database"    : "historical",
        "compute_force_coefficient": true
    })" );

    Params.ValidateAndAssignDefaults(default_parameters);
    mFreestreamPressure = Params["freestream_pressure"].GetDouble();
    mComputeForces = Params["compute_force_coefficient"].GetBool();
    mGetDensity = ProcessDensityDatabaseInput(Params["density_database"].GetString());

    KRATOS_CATCH("")
}

DensityGetter ComputeAerodynamicCoefficientsProcess::ProcessDensityDatabaseInput(const std::string& RequestedDatabase)
{
    KRATOS_TRY

    std::map<std::string, DensityGetter> available_databases;
    available_databases.emplace("non-historical", [](Properties const&, NodeType const& r_node)       { return r_node.GetValue(DENSITY); });
    available_databases.emplace("historical",     [](Properties const&, NodeType const& r_node)       { return r_node.GetSolutionStepValue(DENSITY); });
    available_databases.emplace("properties",     [](Properties const& r_properties, NodeType const&) { return r_properties.GetValue(DENSITY); });

    auto it = available_databases.find(RequestedDatabase);
    
    if (it == available_databases.end())
    {
        std::stringstream error_msg;
        error_msg << "Invalid database. Try any of :" << '\n';
        for(auto const& pair: available_databases)
        {
            error_msg << " - " << pair.first << '\n';
        }
        KRATOS_ERROR << error_msg.str();
    }
    else
    {
        return it->second;
    }

    KRATOS_CATCH("")
}


void ComputeAerodynamicCoefficientsProcess::ExecuteInitialize()
{
    VariableUtils().SetNonHistoricalVariableToZero(PRESSURE_COEFFICIENT, mrModelPart.Nodes());
}

void ComputeAerodynamicCoefficientsProcess::ExecuteFinalizeSolutionStep()
{
    Execute();
}

void ComputeAerodynamicCoefficientsProcess::ComputePressureCoefficient(const Properties& rProperties, NodeType& rNode) const
{
    const auto velocity = rNode.GetSolutionStepValue(VELOCITY);
    const auto pressure = rNode.GetSolutionStepValue(PRESSURE);
    const double density = mGetDensity(rProperties, rNode);

    const double velocity_2 = inner_prod(velocity, velocity);
    const double total_pressure = pressure + 0.5 * density * velocity_2;

    const double cp = (pressure - mFreestreamPressure) / (total_pressure - mFreestreamPressure);
    rNode.SetValue(PRESSURE_COEFFICIENT, cp);
}

array_1d<double, 3> ComputeAerodynamicCoefficientsProcess::IntegrateOnCondition(const Condition& rCondition)
{
    const auto& r_geometry = rCondition.GetGeometry();
    const auto& gauss_points = r_geometry.IntegrationPoints(rCondition.GetIntegrationMethod());

    auto cp = Vector(r_geometry.size());
    std::transform(r_geometry.begin(), r_geometry.end(), cp.begin(), [](NodeType const& r_node) { return r_node.GetValue(PRESSURE_COEFFICIENT); });
    
    Vector N;

    Vector force_coefficient = ZeroVector(3);

    for(auto& gauss_point: gauss_points)
    {
        r_geometry.ShapeFunctionsValues(N, gauss_point.Coordinates());
        const double w = gauss_point.Weight();
        const double detJ = r_geometry.DeterminantOfJacobian(gauss_point.Coordinates());
        const auto n = r_geometry.UnitNormal(gauss_point.Coordinates());

        force_coefficient += w * detJ * inner_prod(N, cp) * n;
    }

    return force_coefficient;
}


void ComputeAerodynamicCoefficientsProcess::Execute()
{
    auto& r_props = [&]() -> Properties&
    {
        if(mrModelPart.NumberOfElements() == 0)
        {
            return mrModelPart.GetProperties(0);
        }
        return mrModelPart.ElementsBegin()->GetProperties();
    }();

    block_for_each(mrModelPart.Nodes(), [&](NodeType& r_node) { ComputePressureCoefficient(r_props, r_node); });

    if(!mComputeForces) return;

    using Reduction = SumReduction<array_1d<double, 3>>;
    const array_1d<double, 3> c_f  = block_for_each<Reduction>(mrModelPart.Conditions(), IntegrateOnCondition);

    r_props.SetValue(FORCE_CM, c_f);
}



}