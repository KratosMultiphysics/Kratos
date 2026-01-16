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

// External includes

// Project includes

// Application includes
#include "compute_lift_coefficient_process.h"
#include "fluid_dynamics_application_variables.h"
#include "utilities/variable_utils.h"


namespace Kratos
{

ComputeLiftCoefficientProcess::ComputeLiftCoefficientProcess(
    Model &rModel,
    Parameters Params)
    : Process(),
      mrModelPart(rModel.GetModelPart(Params["model_part_name"].GetString()))
{
    // Check default settings
    KRATOS_TRY
    
    Params.ValidateAndAssignDefaults(GetDefaultParameters());

    ReadFreestreamValues(Params);

    KRATOS_CATCH("")
}


const Parameters ComputeLiftCoefficientProcess::GetDefaultParameters() const
{
    return Parameters( R"(
    {
        "model_part_name"     : "PLEASE_PROVIDE_A_MODELPART_NAME",
        "reference_surface" : 0.0,
        "angle_of_attack" : 0.0
    })" );
}

void ComputeLiftCoefficientProcess::ReadFreestreamValues(const Parameters& rParams)
{
    constexpr double tol = 1e-12;

    KRATOS_ERROR_IF_NOT(rParams.Has("reference_surface"))
        << "Missing required parameter 'reference_surface'.";

    mReference_Surface = rParams["reference_surface"].GetDouble();
    KRATOS_ERROR_IF(std::abs(mReference_Surface) <= tol)
        << "Invalid value for 'reference_surface' = " << mReference_Surface
        << ". It must be non-zero.";

    KRATOS_ERROR_IF_NOT(rParams.Has("angle_of_attack"))
        << "Missing required parameter 'angle_of_attack'.";

    mAngle_of_Attack = rParams["angle_of_attack"].GetDouble();
}


void ComputeLiftCoefficientProcess::ExecuteBeforeOutputStep()
{
    Execute();
}


void ComputeLiftCoefficientProcess::Execute()
{
    KRATOS_TRY;

    double angle_rad = mAngle_of_Attack*3.14159/(180);
    std::vector<double> e_lift(3);          // direction of the lift
    e_lift[0] = std::sin(angle_rad)*-1;
    e_lift[1] = 0;
    e_lift[2] = std::cos(angle_rad);
    

    double C_L = 0.0;
    double S = 0.0;
    double S_w = 0.0;

    for (auto it_cond = mrModelPart.ConditionsBegin(); it_cond != mrModelPart.ConditionsEnd(); ++it_cond)
    {
        auto& r_geometry = it_cond->GetGeometry();
        const auto integration_method = it_cond->GetIntegrationMethod();
        const auto& r_integration_points = r_geometry.IntegrationPoints(integration_method);
        const auto& r_shape_funct = r_geometry.ShapeFunctionsValues(integration_method);

        double cp_integrated = 0.0;

        // Loop integration points
        for (unsigned int gp = 0; gp < r_integration_points.size(); ++gp)
        {
            const double weight = r_integration_points[gp].Weight();

            // Cp at integration point
            double cp_gp = 0.0;
            for (unsigned int i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node)
            {
                const double cp_node = r_geometry[i_node].GetValue(PRESSURE_COEFFICIENT);
                cp_gp += r_shape_funct(gp, i_node) * cp_node;
            }
            
            S += r_geometry.DeterminantOfJacobian(gp, integration_method);
            S_w += weight * r_geometry.DeterminantOfJacobian(gp, integration_method);
            cp_integrated += cp_gp * weight * r_geometry.DeterminantOfJacobian(gp, integration_method);
        }


        // Normal of the element
        array_1d<double, 3> normal_el;
        normal_el = r_geometry.Normal(r_geometry.Center()); // Non unit vector!
        double norm = std::sqrt(normal_el[0]*normal_el[0] + normal_el[1]*normal_el[1] + normal_el[2]*normal_el[2]);
        normal_el /= norm;

        const double lift_projection = normal_el[0] * e_lift[0] + normal_el[1] * e_lift[1] + normal_el[2] * e_lift[2];

        C_L -= cp_integrated * lift_projection;
    }
    
    C_L = C_L/mReference_Surface;
    mrModelPart.SetValue(LIFT_COEFFICIENT, C_L);
    std::cout << "S_w: " << S_w << std::endl;
    std::cout << "S: " << S << std::endl;
    std::cout << "C_L: " << C_L << std::endl;
    KRATOS_CATCH("");

}


}