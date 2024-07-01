//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquin Gonzalez-Usua
//
//

//These formulas are derived from "Derivatives.py" script

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/checks.h"
#include "utilities/variable_utils.h"
#include "utilities/math_utils.h"

// Application includes
#include "swimming_DEM_application.h"
#include "non_divergence_free_transient_porosity_solution_transient_body_force_process.h"
#include "swimming_dem_application_variables.h"

// Other applications includes
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

/* Public functions *******************************************************/
NonDivergenceFreeTransientPorositySolutionTransientBodyForceProcess::NonDivergenceFreeTransientPorositySolutionTransientBodyForceProcess(
    ModelPart& rModelPart)
    : Process(),
      mrModelPart(rModelPart)
{}

NonDivergenceFreeTransientPorositySolutionTransientBodyForceProcess::NonDivergenceFreeTransientPorositySolutionTransientBodyForceProcess(
    ModelPart& rModelPart,
    Parameters& rParameters)
    : Process(),
      mrModelPart(rModelPart)
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);
}

NonDivergenceFreeTransientPorositySolutionTransientBodyForceProcess::NonDivergenceFreeTransientPorositySolutionTransientBodyForceProcess(
    Model &rModel,
    Parameters &rParameters)
    : Process(),
      mrModelPart(rModel.GetModelPart(rParameters["model_part_name"].GetString()))
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);

}


void NonDivergenceFreeTransientPorositySolutionTransientBodyForceProcess::CheckDefaultsAndProcessSettings(Parameters &rParameters)
{
    const Parameters default_parameters = GetDefaultParameters();

    rParameters.ValidateAndAssignDefaults(default_parameters);

    mDensity       = rParameters["benchmark_parameters"]["density"].GetDouble();
    mViscosity     = rParameters["benchmark_parameters"]["viscosity"].GetDouble();
    mAlphaMin      = rParameters["benchmark_parameters"]["alpha_min"].GetDouble();
    mAlphaMax      = rParameters["benchmark_parameters"]["alpha_max"].GetDouble();
    mSigma         = rParameters["benchmark_parameters"]["sigma"].GetDouble();
    mOmega         = rParameters["benchmark_parameters"]["omega"].GetDouble();
    mk             = rParameters["benchmark_parameters"]["k"].GetDouble();
    mReynolds      = rParameters["benchmark_parameters"]["reynolds"].GetDouble();
    mX1Origin      = rParameters["benchmark_parameters"]["x1_origin"].GetDouble();
    mX2Origin      = rParameters["benchmark_parameters"]["x2_origin"].GetDouble();
    mBumpRadius    = rParameters["benchmark_parameters"]["bump_radius"].GetDouble();
    mPlateauRadius = rParameters["benchmark_parameters"]["plateau_radius"].GetDouble();
    mSqueezeAmplitude = rParameters["benchmark_parameters"]["squeeze_amplitude"].GetDouble();
    mInitialConditions = rParameters["benchmark_parameters"]["use_initial_conditions"].GetBool();
    mAlternativeFormulation = rParameters["benchmark_parameters"]["use_alternative_formulation"].GetBool();

}

const Parameters NonDivergenceFreeTransientPorositySolutionTransientBodyForceProcess::GetDefaultParameters() const
{
    const Parameters default_parameters( R"(
    {
                "model_part_name"          : "please_specify_model_part_name",
                "variable_name"            : "BODY_FORCE",
                "benchmark_name"           : "custom_body_force.vortex",
                "benchmark_parameters"     : {
                                                "velocity"    : 1.0,
                                                "viscosity"   : 0.1,
                                                "density"     : 1.0,
                                                "alpha_min"   : 1e-9,
                                                "alpha_max"   : 1.0,
                                                "sigma"       : 0.0,
                                                "k"           : 100.0,
                                                "Reynolds"    : 1000.0,
                                                "x1_origin"   : 0.5,
                                                "x2_origin"   : 0.5,
                                                "bump_radius" : 0.5,
                                                "plateau_radius" : 0.1,
                                                "use_initial_conditions" : true,
                                                "use_alternative_formulation" : false
                },
                "compute_nodal_error"      : true,
                "print_convergence_output" : false,
                "output_parameters"        : {}
    }  )" );

    return default_parameters;
}

void NonDivergenceFreeTransientPorositySolutionTransientBodyForceProcess::Execute()
{
    this->ExecuteInitialize();
    this->ExecuteInitializeSolutionStep();
}

void NonDivergenceFreeTransientPorositySolutionTransientBodyForceProcess::ExecuteInitialize()
{}

void NonDivergenceFreeTransientPorositySolutionTransientBodyForceProcess::ExecuteBeforeSolutionLoop()
{
    if (mInitialConditions == true)
    {
        this->SetBodyForceAndPorosityField();
    //     this->SetValuesOnIntegrationPoints();
    }
}

void NonDivergenceFreeTransientPorositySolutionTransientBodyForceProcess::ExecuteInitializeSolutionStep()
{
    this->SetBodyForceAndPorosityField();
    this->SetValuesOnIntegrationPoints();
}

void NonDivergenceFreeTransientPorositySolutionTransientBodyForceProcess::ExecuteFinalizeSolutionStep() {}

/* Protected functions ****************************************************/

void NonDivergenceFreeTransientPorositySolutionTransientBodyForceProcess::SetBodyForceAndPorosityField()
{
    const double time = mrModelPart.GetProcessInfo()[TIME];
    const int step = mrModelPart.GetProcessInfo()[STEP];
    const int Dim = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
    const double delta_time = mrModelPart.GetProcessInfo()[DELTA_TIME];
    const double alpha_min = mAlphaMin;
    const double rho = mDensity;
    const double nu = mViscosity;
    const double Re = mReynolds;
    const double Da = mSigma/(alpha_min * nu);
    const double omega = mOmega;
    const double k = mk;
    Matrix I = IdentityMatrix(Dim, Dim);
    const double max_alphat = omega*(1-alpha_min);
    const double Uc = max_alphat*std::sqrt(1/(8*std::pow(k,2))+2)/alpha_min;

    double du1dt, du2dt, du11, du12, du21, du22, du111, du112, du121, du122, du211, du212, du221, du222;
    // Computation of the BodyForce and Porosity fields
    for (auto it_node = mrModelPart.NodesBegin(); it_node != mrModelPart.NodesEnd(); it_node++){

        const double x1 = it_node->X();
        const double x2 = it_node->Y();

        double& r_mass_source = it_node->FastGetSolutionStepValue(MASS_SOURCE);

        double& r_alpha = it_node->FastGetSolutionStepValue(FLUID_FRACTION);
        double& r_dalphat = it_node->FastGetSolutionStepValue(FLUID_FRACTION_RATE);

        double& r_alpha1 = it_node->FastGetSolutionStepValue(FLUID_FRACTION_GRADIENT_X);
        double& r_alpha2 = it_node->FastGetSolutionStepValue(FLUID_FRACTION_GRADIENT_Y);

        double& r_body_force1 = it_node->FastGetSolutionStepValue(BODY_FORCE_X);
        double& r_body_force2 = it_node->FastGetSolutionStepValue(BODY_FORCE_Y);


        double r_u1;/* = it_node->FastGetSolutionStepValue(EXACT_VELOCITY)[0];*/
        double r_u2; /*= it_node->FastGetSolutionStepValue(EXACT_VELOCITY)[1];*/

        double& r_pressure = it_node->FastGetSolutionStepValue(EXACT_PRESSURE);

        double r_sigma; /*= it_node->FastGetSolutionStepValue(HYDRODYNAMIC_REACTION)[0];*/

        if (step == 0){

            double& r_u11 = it_node->FastGetSolutionStepValue(VELOCITY_X,1);
            double& r_u21 = it_node->FastGetSolutionStepValue(VELOCITY_Y,1);

            double& r_u12 = it_node->FastGetSolutionStepValue(VELOCITY_X,2);
            double& r_u22 = it_node->FastGetSolutionStepValue(VELOCITY_Y,2);

            r_u11 = -omega*(alpha_min - 1)*(x1*std::cos(omega*(time-delta_time)) + std::cos(2*k*(x1 + x2) + 2*omega*(time-delta_time))/(4*k))/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*(time-delta_time)),2));

            r_u21 = -omega*(alpha_min - 1)*(-x2*std::cos(omega*(time-delta_time)) + std::cos(2*k*(x1 + x2) + 2*omega*(time-delta_time))/(4*k))/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*(time-delta_time)),2));

            r_u12 = -omega*(alpha_min - 1)*(x1*std::cos(omega*(time-2.0*delta_time)) + std::cos(2*k*(x1 + x2) + 2*omega*(time-2.0*delta_time))/(4*k))/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*(time-2.0*delta_time)),2));

            r_u22 = -omega*(alpha_min - 1)*(-x2*std::cos(omega*(time-2.0*delta_time)) + std::cos(2*k*(x1 + x2) + 2*omega*(time-2.0*delta_time))/(4*k))/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*(time-2.0*delta_time)),2));

        }
        r_alpha = alpha_min + (1 - alpha_min)*std::pow(std::sin(omega*time + k*(x1 + x2)),2);
        r_u1 = -omega*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

        r_u2 = -omega*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

        du1dt = 2*std::pow(omega,2)*(1 - alpha_min)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) - omega*(alpha_min - 1)*(-omega*x1*std::sin(omega*time) - omega*std::sin(2*k*(x1 + x2) + 2*omega*time)/(2*k))/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

        du2dt = 2*std::pow(omega,2)*(1 - alpha_min)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) - omega*(alpha_min - 1)*(omega*x2*std::sin(omega*time) - omega*std::sin(2*k*(x1 + x2) + 2*omega*time)/(2*k))/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

        r_dalphat = 2*omega*(1 - alpha_min)*std::sin(omega*time + k*(x1 + x2))*std::cos(omega*time + k*(x1 + x2));

        r_alpha1 = 2*k*(1 - alpha_min)*std::sin(omega*time + k*(x1 + x2))*std::cos(omega*time + k*(x1 + x2));

        r_alpha2 = 2*k*(1 - alpha_min)*std::sin(omega*time + k*(x1 + x2))*std::cos(omega*time + k*(x1 + x2));

        du11 = 2*k*omega*(1 - alpha_min)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) - omega*(alpha_min - 1)*(-std::sin(2*k*(x1 + x2) + 2*omega*time)/2 + std::cos(omega*time))/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

        du12 = 2*k*omega*(1 - alpha_min)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + omega*(alpha_min - 1)*std::sin(2*k*(x1 + x2) + 2*omega*time)/(2*(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)));

        du111 = -8*std::pow(k,2)*omega*std::pow((1 - alpha_min),2)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),3) - 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 4*k*omega*(1 - alpha_min)*(alpha_min - 1)*(-std::sin(2*k*(x1 + x2) + 2*omega*time)/2 + std::cos(omega*time))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + k*omega*(alpha_min - 1)*std::cos(2*k*(x1 + x2) + 2*omega*time)/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

        du112 = -8*std::pow(k,2)*omega*std::pow((1 - alpha_min),2)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),3) - 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*k*omega*(1 - alpha_min)*(alpha_min - 1)*(-std::sin(2*k*(x1 + x2) + 2*omega*time)/2 + std::cos(omega*time))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) - k*omega*(1 - alpha_min)*(alpha_min - 1)*std::sin(k*(x1 + x2) + omega*time)*std::sin(2*k*(x1 + x2) + 2*omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + k*omega*(alpha_min - 1)*std::cos(2*k*(x1 + x2) + 2*omega*time)/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

        du121 = -8*std::pow(k,2)*omega*std::pow((1 - alpha_min),2)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),3) - 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*k*omega*(1 - alpha_min)*(alpha_min - 1)*(-std::sin(2*k*(x1 + x2) + 2*omega*time)/2 + std::cos(omega*time))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) - k*omega*(1 - alpha_min)*(alpha_min - 1)*std::sin(k*(x1 + x2) + omega*time)*std::sin(2*k*(x1 + x2) + 2*omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + k*omega*(alpha_min - 1)*std::cos(2*k*(x1 + x2) + 2*omega*time)/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

        du122 = -8*std::pow(k,2)*omega*std::pow((1 - alpha_min),2)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),3) - 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) - 2*k*omega*(1 - alpha_min)*(alpha_min - 1)*std::sin(k*(x1 + x2) + omega*time)*std::sin(2*k*(x1 + x2) + 2*omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + k*omega*(alpha_min - 1)*std::cos(2*k*(x1 + x2) + 2*omega*time)/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

        du21 = 2*k*omega*(1 - alpha_min)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + omega*(alpha_min - 1)*std::sin(2*k*(x1 + x2) + 2*omega*time)/(2*(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)));

        du22 = 2*k*omega*(1 - alpha_min)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) - omega*(alpha_min - 1)*(-std::sin(2*k*(x1 + x2) + 2*omega*time)/2 - std::cos(omega*time))/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

        du211 = -8*std::pow(k,2)*omega*std::pow((1 - alpha_min),2)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),3) - 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) - 2*k*omega*(1 - alpha_min)*(alpha_min - 1)*std::sin(k*(x1 + x2) + omega*time)*std::sin(2*k*(x1 + x2) + 2*omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + k*omega*(alpha_min - 1)*std::cos(2*k*(x1 + x2) + 2*omega*time)/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

        du212 = -8*std::pow(k,2)*omega*std::pow((1 - alpha_min),2)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),3) - 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*k*omega*(1 - alpha_min)*(alpha_min - 1)*(-std::sin(2*k*(x1 + x2) + 2*omega*time)/2 - std::cos(omega*time))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) - k*omega*(1 - alpha_min)*(alpha_min - 1)*std::sin(k*(x1 + x2) + omega*time)*std::sin(2*k*(x1 + x2) + 2*omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + k*omega*(alpha_min - 1)*std::cos(2*k*(x1 + x2) + 2*omega*time)/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

        du221 = -8*std::pow(k,2)*omega*std::pow((1 - alpha_min),2)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),3) - 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*k*omega*(1 - alpha_min)*(alpha_min - 1)*(-std::sin(2*k*(x1 + x2) + 2*omega*time)/2 - std::cos(omega*time))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) - k*omega*(1 - alpha_min)*(alpha_min - 1)*std::sin(k*(x1 + x2) + omega*time)*std::sin(2*k*(x1 + x2) + 2*omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + k*omega*(alpha_min - 1)*std::cos(2*k*(x1 + x2) + 2*omega*time)/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

        du222 = -8*std::pow(k,2)*omega*std::pow((1 - alpha_min),2)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),3) - 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 4*k*omega*(1 - alpha_min)*(alpha_min - 1)*(-std::sin(2*k*(x1 + x2) + 2*omega*time)/2 - std::cos(omega*time))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + k*omega*(alpha_min - 1)*std::cos(2*k*(x1 + x2) + 2*omega*time)/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

        r_pressure = Uc*nu*std::sin(Globals::Pi*x2)*std::cos(Globals::Pi*x1)*(1.0 + Re + Da);

        r_sigma = mSigma;

        const double convective1 = r_u1 * du11 + r_u2 * du12;
        const double convective2 = r_u1 * du21 + r_u2 * du22;

        const double div_of_sym_grad1 = (1.0/2.0) * (2.0 * du111 + du212 + du122);
        const double div_of_sym_grad2 = (1.0/2.0) * (du121 + du211 + 2.0 * du222);

        const double grad_of_div1 = du111 + du221;
        const double grad_of_div2 = du112 + du222;

        const double press_grad1 = -Uc*nu*Globals::Pi*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)*(1.0 + Re + Da);
        const double press_grad2 = Uc*nu*Globals::Pi*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)*(1.0 + Re + Da);

        if (mAlternativeFormulation){
            const double grad_alpha_sym_grad1 = (1.0/2.0) * (2 * r_alpha1 * du11 + r_alpha2 * (du21 + du12));
            const double grad_alpha_sym_grad2 = (1.0/2.0) * (r_alpha1 * (du12 + du21) + 2 * r_alpha2 * du22);

            const double grad_alpha_div1 = r_alpha1 * (du11 + du22);
            const double grad_alpha_div2 = r_alpha2 * (du11 + du22);

            r_body_force1 = r_alpha * du1dt + r_alpha * convective1 + r_alpha / rho * press_grad1 - 2.0 * nu * (r_alpha * div_of_sym_grad1 + grad_alpha_sym_grad1) + (2.0/3.0) * nu * (r_alpha * grad_of_div1 + grad_alpha_div1) + r_sigma * r_u1 + r_sigma * r_u2;

            r_body_force2 = r_alpha * du2dt + r_alpha * convective2 + r_alpha / rho * press_grad2 - 2.0 * nu * (r_alpha * div_of_sym_grad2 + grad_alpha_sym_grad2) + (2.0/3.0) * nu * (r_alpha * grad_of_div2 + grad_alpha_div2) + r_sigma * r_u1 + r_sigma * r_u2;

        }else{
            r_sigma /= r_alpha;

            r_body_force1 = du1dt + convective1 + 1.0/rho * press_grad1 - 2.0 * nu * div_of_sym_grad1 + (2.0/3.0) * nu * grad_of_div1 + r_sigma * r_u1 + r_sigma * r_u2;

            r_body_force2 = du2dt + convective2 + 1.0/rho * press_grad2 - 2.0 * nu * div_of_sym_grad2 + (2.0/3.0) * nu * grad_of_div2 + r_sigma * r_u1 + r_sigma * r_u2;
        }

        r_mass_source = r_dalphat + r_u1 * r_alpha1 + r_u2 * r_alpha2 + r_alpha * (du11 + du22);
        if (step == 0){

            double& r_U1 = it_node->FastGetSolutionStepValue(VELOCITY_X);
            double& r_U2 = it_node->FastGetSolutionStepValue(VELOCITY_Y);

            double& r_P = it_node->FastGetSolutionStepValue(PRESSURE);

            double& r_dU1dt =  it_node->FastGetSolutionStepValue(ACCELERATION_X);
            double& r_dU2dt =  it_node->FastGetSolutionStepValue(ACCELERATION_Y);

            r_dU1dt = du1dt;
            r_dU2dt = du2dt;

            r_U1 = r_u1;
            r_U2 = r_u2;
            r_P = r_pressure;

        }
        it_node->SetLock();
        it_node->FastGetSolutionStepValue(EXACT_VELOCITY) = ZeroVector(3);
        it_node->FastGetSolutionStepValue(EXACT_VELOCITY)[0] = r_u1;
        it_node->FastGetSolutionStepValue(EXACT_VELOCITY)[1] = r_u2;
        it_node->FastGetSolutionStepValue(HYDRODYNAMIC_REACTION) = ZeroVector(3);
        it_node->FastGetSolutionStepValue(HYDRODYNAMIC_REACTION)[0] = r_sigma;
        it_node->UnSetLock();
    }

}

void NonDivergenceFreeTransientPorositySolutionTransientBodyForceProcess::SetValuesOnIntegrationPoints()
{
    const double time = mrModelPart.GetProcessInfo()[TIME];
    const unsigned int Dim = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
    const double alpha_min = mAlphaMin;
    const double rho = mDensity;
    const double nu = mViscosity;
    const double Re = mReynolds;
    const double Da = mSigma/(alpha_min * nu);
    const double omega = mOmega;
    const double k = mk;
    Matrix I = IdentityMatrix(Dim, Dim);
    Matrix sigma = ZeroMatrix(Dim, Dim);
    const double max_alphat = omega*(1-alpha_min);
    const double Uc = max_alphat*std::sqrt(1/(8*std::pow(k,2))+2)/alpha_min;

    double alpha, alpha1, alpha2, u1, u2, pressure, du1dt, du2dt, du11, du12, du21, du22, du111, du112, du121, du122, du211, du212, du221, du222, dalphat, body_force1, body_force2;

    unsigned int n_elem = mrModelPart.NumberOfElements();

    // Computation of the BodyForce and Porosity fields
    for (int i_elem = 0; i_elem < static_cast<int>(mrModelPart.NumberOfElements()); ++i_elem){

        const auto it_elem = mrModelPart.ElementsBegin() + i_elem;
        int id_elem = it_elem->Id();
        const GeometryType& r_geometry = it_elem->GetGeometry();

        const GeometryData::IntegrationMethod integration_method = it_elem->GetIntegrationMethod();
        const auto& r_number_integration_points = r_geometry.IntegrationPointsNumber(integration_method);

        Matrix NContainer = r_geometry.ShapeFunctionsValues(integration_method);
        const unsigned int NumNodes = r_geometry.PointsNumber();

        std::vector<double> pressure_on_gauss_points;
        std::vector<double> fluid_fraction_on_gauss_points;
        std::vector<double> fluid_fraction_rate_on_gauss_points;
        std::vector<Vector> velocity_on_gauss_points(r_number_integration_points);
        Matrix body_force_on_gauss_points = ZeroMatrix(Dim, r_number_integration_points);
        std::vector<array_1d<double,3>> pressure_gradient_on_gauss_points(r_number_integration_points);
        Matrix fluid_fraction_gradient_on_gauss_points = ZeroMatrix(Dim, r_number_integration_points);
        std::vector<Matrix> velocity_gradient_on_gauss_points(r_number_integration_points);

        // if (mExactScalar.size() != n_elem)
        //     mExactScalar.resize(n_elem);

        // if (mExactVector.size() != n_elem)
        //     mExactVector.resize(n_elem);

        // if (mExactScalarGradient.size() != n_elem)
        //     mExactScalarGradient.resize(n_elem);

        // if (mExactPorosity.size() != n_elem)
        //     mExactPorosity.resize(n_elem);

        // if (mExactPorosityRate.size() != n_elem)
        //     mExactPorosityRate.resize(n_elem);

        // if (mExactPorosityGradient.size() != n_elem)
        //     mExactPorosityGradient.resize(n_elem);

        // if (mExactVectorGradient.size() != n_elem)
        //     mExactVectorGradient.resize(n_elem);

        // if (mExactBodyForce.size() != n_elem)
        //     mExactBodyForce.resize(n_elem);

        if (pressure_on_gauss_points.size() != r_number_integration_points)
            pressure_on_gauss_points.resize(r_number_integration_points);

        if (fluid_fraction_on_gauss_points.size() != r_number_integration_points)
            fluid_fraction_on_gauss_points.resize(r_number_integration_points);

        if (fluid_fraction_rate_on_gauss_points.size() != r_number_integration_points)
            fluid_fraction_rate_on_gauss_points.resize(r_number_integration_points);

        for (unsigned int g = 0; g < r_number_integration_points; g++){

            Matrix gauss_point_coordinates = ZeroMatrix(r_number_integration_points,Dim);
            for (unsigned int i = 0; i < NumNodes; ++i){
                const array_1d<double, 3>& r_coordinates = r_geometry[i].Coordinates();
                for (unsigned int d = 0; d < Dim; ++d)
                    gauss_point_coordinates(g,d) += NContainer(g,i) * r_coordinates[d];
                }

            velocity_gradient_on_gauss_points[g] = ZeroMatrix(3,3);
            pressure_gradient_on_gauss_points[g] = ZeroVector(3);
            velocity_on_gauss_points[g] = ZeroVector(3);

            const double x1 = gauss_point_coordinates(g,0);
            const double x2 = gauss_point_coordinates(g,1);

            alpha = alpha_min + (1 - alpha_min)*std::pow(std::sin(omega*time + k*(x1 + x2)),2);

            alpha = alpha_min + (1 - alpha_min)*std::pow(std::sin(omega*time + k*(x1 + x2)),2);

            u1 = -omega*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

            u2 = -omega*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

            du1dt = 2*std::pow(omega,2)*(1 - alpha_min)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) - omega*(alpha_min - 1)*(-omega*x1*std::sin(omega*time) - omega*std::sin(2*k*(x1 + x2) + 2*omega*time)/(2*k))/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

            du2dt = 2*std::pow(omega,2)*(1 - alpha_min)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) - omega*(alpha_min - 1)*(omega*x2*std::sin(omega*time) - omega*std::sin(2*k*(x1 + x2) + 2*omega*time)/(2*k))/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

            dalphat = 2*omega*(1 - alpha_min)*std::sin(omega*time + k*(x1 + x2))*std::cos(omega*time + k*(x1 + x2));

            alpha1 = 2*k*(1 - alpha_min)*std::sin(omega*time + k*(x1 + x2))*std::cos(omega*time + k*(x1 + x2));

            alpha2 = 2*k*(1 - alpha_min)*std::sin(omega*time + k*(x1 + x2))*std::cos(omega*time + k*(x1 + x2));

            du11 = 2*k*omega*(1 - alpha_min)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) - omega*(alpha_min - 1)*(-std::sin(2*k*(x1 + x2) + 2*omega*time)/2 + std::cos(omega*time))/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

            du12 = 2*k*omega*(1 - alpha_min)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + omega*(alpha_min - 1)*std::sin(2*k*(x1 + x2) + 2*omega*time)/(2*(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)));

            du111 = -8*std::pow(k,2)*omega*std::pow((1 - alpha_min),2)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),3) - 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 4*k*omega*(1 - alpha_min)*(alpha_min - 1)*(-std::sin(2*k*(x1 + x2) + 2*omega*time)/2 + std::cos(omega*time))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + k*omega*(alpha_min - 1)*std::cos(2*k*(x1 + x2) + 2*omega*time)/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

            du112 = -8*std::pow(k,2)*omega*std::pow((1 - alpha_min),2)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),3) - 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*k*omega*(1 - alpha_min)*(alpha_min - 1)*(-std::sin(2*k*(x1 + x2) + 2*omega*time)/2 + std::cos(omega*time))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) - k*omega*(1 - alpha_min)*(alpha_min - 1)*std::sin(k*(x1 + x2) + omega*time)*std::sin(2*k*(x1 + x2) + 2*omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + k*omega*(alpha_min - 1)*std::cos(2*k*(x1 + x2) + 2*omega*time)/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

            du121 = -8*std::pow(k,2)*omega*std::pow((1 - alpha_min),2)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),3) - 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*k*omega*(1 - alpha_min)*(alpha_min - 1)*(-std::sin(2*k*(x1 + x2) + 2*omega*time)/2 + std::cos(omega*time))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) - k*omega*(1 - alpha_min)*(alpha_min - 1)*std::sin(k*(x1 + x2) + omega*time)*std::sin(2*k*(x1 + x2) + 2*omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + k*omega*(alpha_min - 1)*std::cos(2*k*(x1 + x2) + 2*omega*time)/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

            du122 = -8*std::pow(k,2)*omega*std::pow((1 - alpha_min),2)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),3) - 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(x1*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) - 2*k*omega*(1 - alpha_min)*(alpha_min - 1)*std::sin(k*(x1 + x2) + omega*time)*std::sin(2*k*(x1 + x2) + 2*omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + k*omega*(alpha_min - 1)*std::cos(2*k*(x1 + x2) + 2*omega*time)/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

            du21 = 2*k*omega*(1 - alpha_min)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + omega*(alpha_min - 1)*std::sin(2*k*(x1 + x2) + 2*omega*time)/(2*(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)));

            du22 = 2*k*omega*(1 - alpha_min)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) - omega*(alpha_min - 1)*(-std::sin(2*k*(x1 + x2) + 2*omega*time)/2 - std::cos(omega*time))/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

            du211 = -8*std::pow(k,2)*omega*std::pow((1 - alpha_min),2)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),3) - 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) - 2*k*omega*(1 - alpha_min)*(alpha_min - 1)*std::sin(k*(x1 + x2) + omega*time)*std::sin(2*k*(x1 + x2) + 2*omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + k*omega*(alpha_min - 1)*std::cos(2*k*(x1 + x2) + 2*omega*time)/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

            du212 = -8*std::pow(k,2)*omega*std::pow((1 - alpha_min),2)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),3) - 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*k*omega*(1 - alpha_min)*(alpha_min - 1)*(-std::sin(2*k*(x1 + x2) + 2*omega*time)/2 - std::cos(omega*time))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) - k*omega*(1 - alpha_min)*(alpha_min - 1)*std::sin(k*(x1 + x2) + omega*time)*std::sin(2*k*(x1 + x2) + 2*omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + k*omega*(alpha_min - 1)*std::cos(2*k*(x1 + x2) + 2*omega*time)/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

            du221 = -8*std::pow(k,2)*omega*std::pow((1 - alpha_min),2)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),3) - 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*k*omega*(1 - alpha_min)*(alpha_min - 1)*(-std::sin(2*k*(x1 + x2) + 2*omega*time)/2 - std::cos(omega*time))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) - k*omega*(1 - alpha_min)*(alpha_min - 1)*std::sin(k*(x1 + x2) + omega*time)*std::sin(2*k*(x1 + x2) + 2*omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + k*omega*(alpha_min - 1)*std::cos(2*k*(x1 + x2) + 2*omega*time)/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

            du222 = -8*std::pow(k,2)*omega*std::pow((1 - alpha_min),2)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),3) - 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::sin(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 2*std::pow(k,2)*omega*(1 - alpha_min)*(alpha_min - 1)*(-x2*std::cos(omega*time) + std::cos(2*k*(x1 + x2) + 2*omega*time)/(4*k))*std::pow(std::cos(k*(x1 + x2) + omega*time),2)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + 4*k*omega*(1 - alpha_min)*(alpha_min - 1)*(-std::sin(2*k*(x1 + x2) + 2*omega*time)/2 - std::cos(omega*time))*std::sin(k*(x1 + x2) + omega*time)*std::cos(k*(x1 + x2) + omega*time)/std::pow((alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2)),2) + k*omega*(alpha_min - 1)*std::cos(2*k*(x1 + x2) + 2*omega*time)/(alpha_min + (1 - alpha_min)*std::pow(std::sin(k*(x1 + x2) + omega*time),2));

            pressure = Uc*nu*std::sin(Globals::Pi*x2)*std::cos(Globals::Pi*x1)*(1.0 + Re + Da);

            const double press_grad1 = -Uc*nu*Globals::Pi*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)*(1.0 + Re + Da);
            const double press_grad2 = Uc*nu*Globals::Pi*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)*(1.0 + Re + Da);

            sigma = mSigma * I;

            const double convective1 = u1 * du11 + u2 * du12;
            const double convective2 = u1 * du21 + u2 * du22;
            const double div_of_sym_grad1 = (1.0/2.0) * (2.0 * du111 + du212 + du122);
            const double div_of_sym_grad2 = (1.0/2.0) * (du121 + du211 + 2.0 * du222);

            const double grad_of_div1 = du111 + du221;
            const double grad_of_div2 = du112 + du222;

            if (mAlternativeFormulation){
                const double grad_alpha_sym_grad1 = (1.0/2.0) * (2 * alpha1 * du11 + alpha2 * (du21 + du12));
                const double grad_alpha_sym_grad2 = (1.0/2.0) * (alpha1 * (du12 + du21) + 2 * alpha2 * du22);

                const double grad_alpha_div1 = alpha1 * (du11 + du22);
                const double grad_alpha_div2 = alpha2 * (du11 + du22);

                body_force1 = alpha * du1dt + alpha * convective1 + alpha / rho * press_grad1 - 2.0 * nu * (alpha * div_of_sym_grad1 + grad_alpha_sym_grad1) + (2.0/3.0) * nu * (alpha * grad_of_div1 + grad_alpha_div1) + sigma(0,0) * u1 + sigma(0,1) * u2;

                body_force2 = alpha * du2dt + alpha * convective2 + alpha / rho * press_grad2 - 2.0 * nu * (alpha * div_of_sym_grad2 + grad_alpha_sym_grad2) + (2.0/3.0) * nu * (alpha * grad_of_div2 + grad_alpha_div2) + sigma(1,0) * u1 + sigma(1,1) * u2;

            }else{
                sigma /= alpha;

                body_force1 = du1dt + convective1 + 1.0/rho * press_grad1 - 2.0 * nu * div_of_sym_grad1 + (2.0/3.0) * nu * grad_of_div1 + sigma(0,0) * u1 + sigma(0,1) * u2;

                body_force2 = du2dt + convective2 + 1.0/rho * press_grad2 - 2.0 * nu * div_of_sym_grad2 + (2.0/3.0) * nu * grad_of_div2 + sigma(1,0) * u1 + sigma(1,1) * u2;
            }

            pressure_on_gauss_points[g] = pressure;
            fluid_fraction_on_gauss_points[g] = alpha;
            fluid_fraction_rate_on_gauss_points[g] = dalphat;
            fluid_fraction_gradient_on_gauss_points(0,g) = alpha1;
            fluid_fraction_gradient_on_gauss_points(1,g) = alpha2;
            velocity_on_gauss_points[g][0] = u1;
            velocity_on_gauss_points[g][1] = u2;
            body_force_on_gauss_points(0,g) = body_force1;
            body_force_on_gauss_points(1,g) = body_force2;
            pressure_gradient_on_gauss_points[g][0] = press_grad1;
            pressure_gradient_on_gauss_points[g][1] = press_grad2;
            velocity_gradient_on_gauss_points[g](0,0) = du11;
            velocity_gradient_on_gauss_points[g](0,1) = du21;
            velocity_gradient_on_gauss_points[g](1,0) = du12;
            velocity_gradient_on_gauss_points[g](1,1) = du22;
        }

        it_elem->SetValuesOnIntegrationPoints(EXACT_VELOCITY, velocity_on_gauss_points, mrModelPart.GetProcessInfo());
        it_elem->SetValuesOnIntegrationPoints(EXACT_PRESSURE, pressure_on_gauss_points, mrModelPart.GetProcessInfo());
        it_elem->SetValuesOnIntegrationPoints(RECOVERED_PRESSURE_GRADIENT, pressure_gradient_on_gauss_points, mrModelPart.GetProcessInfo());
        it_elem->SetValuesOnIntegrationPoints(EXACT_VELOCITY_GRADIENT,velocity_gradient_on_gauss_points, mrModelPart.GetProcessInfo());
        // mExactPorosity[id_elem-1] = fluid_fraction_on_gauss_points;
        // mExactPorosityGradient[id_elem-1] = fluid_fraction_gradient_on_gauss_points;
        // mExactPorosityRate[id_elem-1] = fluid_fraction_rate_on_gauss_points;
        // mExactBodyForce[id_elem-1] = body_force_on_gauss_points;
        // mExactScalar[id_elem-1] = pressure_on_gauss_points;
        // mExactVector[id_elem-1] = velocity_on_gauss_points;
        // mExactScalarGradient[id_elem-1] = pressure_gradient_on_gauss_points;
        // mExactVectorGradient[id_elem-1] = velocity_gradient_on_gauss_points;

    }
}

/* Private functions ****************************************************/

};  // namespace Kratos.
