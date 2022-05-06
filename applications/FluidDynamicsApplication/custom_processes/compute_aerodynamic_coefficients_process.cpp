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
#include <boost/numeric/ublas/vector_expression.hpp>
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
        "freestream_density" : 0,
        "freestream_velocity" : 0,
        "freestream_pressure" : 0,
        "reference_area" : 1.0,
        "lift_direction" : [1.0, 0.0, 0.0],
        "drag_direction" : [0.0, 1.0, 0.0],
        "compute_force_coefficient": true
    })" );

    // WIP: ALLOW CHOICE TO EXECUTE ONLY DURING OUTPUT STEPS!!! (HOW?)

    Params.ValidateAndAssignDefaults(default_parameters);
    mComputeForces = Params["compute_force_coefficient"].GetBool();

    mFreestreamStaticPressure = Params["freestream_pressure"].GetDouble();

    const double freestream_density = Params["freestream_density"].GetDouble();
    const double free_stream_velocity = Params["freestream_velocity"].GetDouble();
    mFreestreamDynamicPressure = 0.5 * freestream_density * free_stream_velocity * free_stream_velocity;

    mLiftDirection = Params["lift_direction"].GetVector();
    mDragDirection = Params["drag_direction"].GetVector();

    mLiftDirection /= norm_2(mLiftDirection);
    mDragDirection /= norm_2(mDragDirection);

    mReferenceArea = Params["reference_area"].GetDouble();

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
    const auto pressure = rNode.GetSolutionStepValue(PRESSURE);
    const double cp = (pressure - mFreestreamStaticPressure) / mFreestreamDynamicPressure;
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
        const auto inwards_normal = - r_geometry.UnitNormal(gauss_point.Coordinates());

        force_coefficient += w * detJ * inner_prod(N, cp) * inwards_normal;
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

    KRATOS_TRY
    {
        block_for_each(mrModelPart.Nodes(), [&](NodeType& r_node) { ComputePressureCoefficient(r_props, r_node); });
    }
    KRATOS_CATCH("Error computing pressure coefficient")

    if(!mComputeForces) return;

    KRATOS_TRY
    {
        using Reduction = SumReduction<array_1d<double, 3>>;
        const array_1d<double, 3> c_f = (1/mReferenceArea) * block_for_each<Reduction>(mrModelPart.Conditions(), IntegrateOnCondition);

        r_props.SetValue(LIFT_COEFFICIENT, inner_prod(c_f, mLiftDirection));
        r_props.SetValue(DRAG_COEFFICIENT, inner_prod(c_f, mDragDirection));
    }
    KRATOS_CATCH("Error computing lift and drag coefficients")
}


}