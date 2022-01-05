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
#include "sinusoidal_porosity_solution_transient_body_force_process.h"
#include "swimming_dem_application_variables.h"

// Other applications includes
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

/* Public functions *******************************************************/
SinusoidalPorositySolutionTransientBodyForceProcess::SinusoidalPorositySolutionTransientBodyForceProcess(
    ModelPart& rModelPart)
    : Process(),
      mrModelPart(rModelPart)
{}

SinusoidalPorositySolutionTransientBodyForceProcess::SinusoidalPorositySolutionTransientBodyForceProcess(
    ModelPart& rModelPart,
    Parameters& rParameters)
    : Process(),
      mrModelPart(rModelPart)
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);
}

SinusoidalPorositySolutionTransientBodyForceProcess::SinusoidalPorositySolutionTransientBodyForceProcess(
    Model &rModel,
    Parameters &rParameters)
    : Process(),
      mrModelPart(rModel.GetModelPart(rParameters["model_part_name"].GetString()))
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);

}


void SinusoidalPorositySolutionTransientBodyForceProcess::CheckDefaultsAndProcessSettings(Parameters &rParameters)
{
    const Parameters default_parameters = GetDefaultParameters();

    rParameters.ValidateAndAssignDefaults(default_parameters);

    mDensity     = rParameters["benchmark_parameters"]["density"].GetDouble();
    mUchar       = rParameters["benchmark_parameters"]["u_char"].GetDouble();
    mDeltaAlpha  = rParameters["benchmark_parameters"]["delta_alpha"].GetDouble();
    mLength      = rParameters["benchmark_parameters"]["length"].GetDouble();
    mOmega       = rParameters["benchmark_parameters"]["omega"].GetDouble();
    mX1Origin    = rParameters["benchmark_parameters"]["x1_origin"].GetDouble();
    mX2Origin    = rParameters["benchmark_parameters"]["x2_origin"].GetDouble();
    mSqueezeAmplitude = rParameters["benchmark_parameters"]["squeeze_amplitude"].GetDouble();
    mNSafety     = rParameters["benchmark_parameters"]["n_safety"].GetDouble();
    mReynoldsNumber = rParameters["benchmark_parameters"]["n_reynolds"].GetDouble();
    mDamKohlerNumber = rParameters["benchmark_parameters"]["n_dam"].GetDouble();
    mMaxGradAlpha = rParameters["benchmark_parameters"]["max_grad_alpha"].GetDouble();
    mInitialConditions = rParameters["benchmark_parameters"]["use_initial_conditions"].GetBool();
    mAlternativeFormulation = rParameters["benchmark_parameters"]["use_alternative_formulation"].GetBool();

    double &r_nu = mViscosity;

    this->CalculateKinematicViscosity(mReynoldsNumber, r_nu);

    double dynamic_viscosity = r_nu * mDensity;

    double &permeability = mPermeability;

    this->CalculatePermeability(mDamKohlerNumber, dynamic_viscosity, permeability);

    this->CalculateWaveNumber(mMaxGradAlpha, mReynoldsNumber, dynamic_viscosity, mDeltaAlpha, mUchar, mWaveNumber);

}

void SinusoidalPorositySolutionTransientBodyForceProcess::CalculateKinematicViscosity(
    double &rReynoldsNumber,
    double &rviscosity)
{
    double u_characteristic = mUchar;
    double L = mLength;
    rviscosity = u_characteristic * L / rReynoldsNumber;
}

void SinusoidalPorositySolutionTransientBodyForceProcess::CalculatePermeability(
    double &rDamKohlerNumber,
    double &dynamic_viscosity,
    double &permeability)
{
    const double L = mLength;
    const double nu = mViscosity;
    const double u_char = mUchar;

    permeability = dynamic_viscosity * u_char / (rDamKohlerNumber * (2 * nu * (u_char/std::pow(L,2))));

}

void SinusoidalPorositySolutionTransientBodyForceProcess::CalculateWaveNumber(
    double &mMaxGradAlpha,
    double &mReynoldsNumber,
    double &dynamic_viscosity,
    double &mDeltaAlpha,
    double &mUchar,
    double &mWaveNumber)
{
    mWaveNumber = mMaxGradAlpha * mUchar / (mReynoldsNumber * mDeltaAlpha * dynamic_viscosity);

}

const Parameters SinusoidalPorositySolutionTransientBodyForceProcess::GetDefaultParameters() const
{
    const Parameters default_parameters( R"(
    {
                "model_part_name"          : "please_specify_model_part_name",
                "variable_name"            : "BODY_FORCE",
                "benchmark_name"           : "custom_body_force.vortex",
                "benchmark_parameters"     : {
                                                "velocity"    : 1.0,
                                                "length"      : 1.0,
                                                "density"     : 1.0,
                                                "frequency"   : 1.0,
                                                "alpha0"      : 0.7,
                                                "alpha_min"   : 0.5,
                                                "period"      : 0.1,
                                                "delta_alpha" : 0.25,
                                                "squeeze_amplitude" : 0.425,
                                                "x1_origin"   : 0.5,
                                                "x2_origin"   : 0.5,
                                                "omega"       : 5.0,
                                                "n_safety"    : 1.0,
                                                "n_reynolds"  : 1000.0,
                                                "n_dam"       : 0.0001,
                                                "u_char"      : 100.0,
                                                "use_alternative_formulation" : false
                },
                "compute_nodal_error"      : true,
                "print_convergence_output" : false,
                "output_parameters"        : {}
    }  )" );

    return default_parameters;
}

void SinusoidalPorositySolutionTransientBodyForceProcess::Execute()
{
    this->ExecuteInitialize();
    this->ExecuteInitializeSolutionStep();
}

void SinusoidalPorositySolutionTransientBodyForceProcess::ExecuteInitialize()
{}

void SinusoidalPorositySolutionTransientBodyForceProcess::ExecuteBeforeSolutionLoop()
{
    this->SetFluidProperties();
    if (mInitialConditions == true)
    {
        this->SetInitialBodyForceAndPorosityField();
    }
}

void SinusoidalPorositySolutionTransientBodyForceProcess::ExecuteInitializeSolutionStep()
{
    this->SetBodyForceAndPorosityField();
}

void SinusoidalPorositySolutionTransientBodyForceProcess::ExecuteFinalizeSolutionStep() {}

/* Protected functions ****************************************************/

void SinusoidalPorositySolutionTransientBodyForceProcess::SetInitialBodyForceAndPorosityField()
{
    const double time = mrModelPart.GetProcessInfo()[TIME];
    const double dim = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
    const double delta_alpha = mDeltaAlpha;
    const double rho = mDensity;
    const double nu = mViscosity;
    const double k = mWaveNumber;
    double u_char = mUchar;

    Matrix inv_permeability = ZeroMatrix(dim,dim);

    double du1dt, du2dt, du11, du12, du111, du112, du121, du122, du21, du22, du211, du212, du221, du222;

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

        double& r_u1 = it_node->FastGetSolutionStepValue(EXACT_VELOCITY_X);
        double& r_u2 = it_node->FastGetSolutionStepValue(EXACT_VELOCITY_Y);

        Matrix& permeability = it_node->FastGetSolutionStepValue(PERMEABILITY);

        double& r_pressure = it_node->FastGetSolutionStepValue(EXACT_PRESSURE);


        r_u1 = -u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        r_u2 = -u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du1dt = u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du2dt = u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        r_alpha = delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5;

        r_dalphat = 0.0;

        r_alpha1 = delta_alpha*k*std::cos(k*x1)*std::cos(k*x2);

        r_alpha2 = -delta_alpha*k*std::sin(k*x1)*std::sin(k*x2);

        du11 = delta_alpha*k*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - u_char*std::pow(x1,2)*(2*x1 - 2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5) - 2*u_char*x1*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du12 = -delta_alpha*k*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du111 = -2*std::pow(delta_alpha,2)*std::pow(k,2)*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::pow(std::cos(k*x1),2)*std::pow(std::cos(k*x2),2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),3) - delta_alpha*std::pow(k,2)*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) + 2*delta_alpha*k*u_char*std::pow(x1,2)*(2*x1 - 2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) + 4*delta_alpha*k*u_char*x1*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - 2*u_char*std::pow(x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5) - 4*u_char*x1*(2*x1 - 2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5) - 2*u_char*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du112 = 2*std::pow(delta_alpha,2)*std::pow(k,2)*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),3) - delta_alpha*std::pow(k,2)*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x2)*std::cos(k*x1)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) + delta_alpha*k*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow(1 - x2,2))*std::exp(-time)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - delta_alpha*k*u_char*std::pow(x1,2)*(2*x1 - 2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - 2*delta_alpha*k*u_char*x1*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - u_char*std::pow(x1,2)*(2*x1 - 2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5) - 2*u_char*x1*std::pow(1 - x1,2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du121 = 2*std::pow(delta_alpha,2)*std::pow(k,2)*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),3) - delta_alpha*std::pow(k,2)*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x2)*std::cos(k*x1)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) + delta_alpha*k*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow(1 - x2,2))*std::exp(-time)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - delta_alpha*k*u_char*std::pow(x1,2)*(2*x1 - 2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - 2*delta_alpha*k*u_char*x1*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - u_char*std::pow(x1,2)*(2*x1 - 2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5) - 2*u_char*x1*std::pow(1 - x1,2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du122 = -2*std::pow(delta_alpha,2)*std::pow(k,2)*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::pow(std::sin(k*x1),2)*std::pow(std::sin(k*x2),2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),3) - delta_alpha*std::pow(k,2)*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - 2*delta_alpha*k*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(12*u_char*x2 + 6*u_char*(2*x2 - 2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du21 = delta_alpha*k*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du22 = -delta_alpha*k*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - u_char*std::pow(x2,2)*(2*x2 - 2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5) - 2*u_char*x2*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du211 = -2*std::pow(delta_alpha,2)*std::pow(k,2)*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::pow(std::cos(k*x1),2)*std::pow(std::cos(k*x2),2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),3) - delta_alpha*std::pow(k,2)*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::sin(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) + 2*delta_alpha*k*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow(1 - x1,2))*std::exp(-time)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-12*u_char*x1 - 6*u_char*(2*x1 - 2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du212 = 2*std::pow(delta_alpha,2)*std::pow(k,2)*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),3) - delta_alpha*std::pow(k,2)*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::sin(k*x2)*std::cos(k*x1)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - delta_alpha*k*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow(1 - x1,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) + delta_alpha*k*u_char*std::pow(x2,2)*(2*x2 - 2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) + 2*delta_alpha*k*u_char*x2*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - u_char*std::pow(x2,2)*(2*x2 - 2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5) - 2*u_char*x2*std::pow(1 - x2,2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du221 = 2*std::pow(delta_alpha,2)*std::pow(k,2)*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),3) - delta_alpha*std::pow(k,2)*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::sin(k*x2)*std::cos(k*x1)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - delta_alpha*k*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow(1 - x1,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) + delta_alpha*k*u_char*std::pow(x2,2)*(2*x2 - 2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) + 2*delta_alpha*k*u_char*x2*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - u_char*std::pow(x2,2)*(2*x2 - 2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5) - 2*u_char*x2*std::pow(1 - x2,2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du222 = -2*std::pow(delta_alpha,2)*std::pow(k,2)*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::pow(std::sin(k*x1),2)*std::pow(std::sin(k*x2),2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),3) - delta_alpha*std::pow(k,2)*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::sin(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - 2*delta_alpha*k*u_char*std::pow(x2,2)*(2*x2 - 2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - 4*delta_alpha*k*u_char*x2*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - 2*u_char*std::pow(x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5) - 4*u_char*x2*(2*x2 - 2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5) - 2*u_char*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        r_pressure = 0.0;

        double det_permeability = MathUtils<double>::Det(permeability);
        MathUtils<double>::InvertMatrix(permeability, inv_permeability, det_permeability, -1.0);

        Matrix sigma = nu * rho * inv_permeability;

        const double convective1 = r_u1 * du11 + r_u2 * du12;
        const double convective2 = r_u1 * du21 + r_u2 * du22;

        const double div_of_sym_grad1 = (1.0/2.0) * (2 * du111 + du212 + du122);
        const double div_of_sym_grad2 = (1.0/2.0) * (du121 + du211 + 2 * du222);

        const double grad_of_div1 = du111 + du221;
        const double grad_of_div2 = du112 + du222;

        const double press_grad1 = 0.0;
        const double press_grad2 = 0.0;

        if (mAlternativeFormulation){

            const double grad_alpha_sym_grad1 = (1.0/2.0) * (2 * r_alpha1 * du11 + r_alpha2 * (du21 + du12));
            const double grad_alpha_sym_grad2 = (1.0/2.0) * (r_alpha1 * (du12 + du21) + 2 * r_alpha2 * du22);

            const double grad_alpha_div1 = r_alpha1 * (du11 + du22);
            const double grad_alpha_div2 = r_alpha2 * (du11 + du22);

            r_body_force1 = r_alpha * du1dt + r_alpha * convective1 + r_alpha / rho * press_grad1 - 2 * nu * (r_alpha * div_of_sym_grad1 + grad_alpha_sym_grad1) + (2.0/3.0) * nu * (r_alpha * grad_of_div1 + grad_alpha_div1) + sigma(0,0) * r_u1 + sigma(1,0) * r_u1;

            r_body_force2 = r_alpha * du2dt + r_alpha * convective2 + r_alpha / rho * press_grad2 - 2 * nu * (r_alpha * div_of_sym_grad2 + grad_alpha_sym_grad2) + (2.0/3.0) * nu * (r_alpha * grad_of_div2 + grad_alpha_div2) + sigma(0,1) * r_u2 + sigma(1,1) * r_u2;

        }else{
            r_body_force1 = du1dt + convective1 + 1.0/rho * press_grad1 - 2 * nu * div_of_sym_grad1 + (2.0/3.0) * nu * grad_of_div1 + sigma(0,0) * r_u1 + sigma(1,0) * r_u1;

            r_body_force2 = du2dt + convective2 + 1.0/rho * press_grad2 - 2 * nu * div_of_sym_grad2 + (2.0/3.0) * nu * grad_of_div2 + sigma(0,1) * r_u2 + sigma(1,1) * r_u2;
        }

        r_mass_source = (r_dalphat + r_u1 * r_alpha1 + r_u2 * r_alpha2 + r_alpha * (du11 + du22));

        it_node->FastGetSolutionStepValue(VELOCITY_X) = r_u1;
        it_node->FastGetSolutionStepValue(VELOCITY_Y) = r_u2;
        it_node->FastGetSolutionStepValue(PRESSURE) = r_pressure;
    }

}

void SinusoidalPorositySolutionTransientBodyForceProcess::SetBodyForceAndPorosityField()
{

    const double time = mrModelPart.GetProcessInfo()[TIME];
    const double dim = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
    const double delta_alpha = mDeltaAlpha;
    const double rho = mDensity;
    const double nu = mViscosity;
    const double k = mWaveNumber;
    double u_char = mUchar;

    Matrix inv_permeability = ZeroMatrix(dim,dim);

    double du1dt, du2dt, du11, du12, du111, du112, du121, du122, du21, du22, du211, du212, du221, du222;

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

        double& r_u1 = it_node->FastGetSolutionStepValue(EXACT_VELOCITY_X);
        double& r_u2 = it_node->FastGetSolutionStepValue(EXACT_VELOCITY_Y);

        double& r_pressure = it_node->FastGetSolutionStepValue(EXACT_PRESSURE);

        Matrix& permeability = it_node->FastGetSolutionStepValue(PERMEABILITY);

        r_u1 = -u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        r_u2 = -u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du1dt = u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du2dt = u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        r_alpha = delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5;

        r_dalphat = 0.0;

        r_alpha1 = delta_alpha*k*std::cos(k*x1)*std::cos(k*x2);

        r_alpha2 = -delta_alpha*k*std::sin(k*x1)*std::sin(k*x2);

        du11 = delta_alpha*k*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - u_char*std::pow(x1,2)*(2*x1 - 2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5) - 2*u_char*x1*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du12 = -delta_alpha*k*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du111 = -2*std::pow(delta_alpha,2)*std::pow(k,2)*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::pow(std::cos(k*x1),2)*std::pow(std::cos(k*x2),2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),3) - delta_alpha*std::pow(k,2)*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) + 2*delta_alpha*k*u_char*std::pow(x1,2)*(2*x1 - 2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) + 4*delta_alpha*k*u_char*x1*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - 2*u_char*std::pow(x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5) - 4*u_char*x1*(2*x1 - 2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5) - 2*u_char*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du112 = 2*std::pow(delta_alpha,2)*std::pow(k,2)*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),3) - delta_alpha*std::pow(k,2)*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x2)*std::cos(k*x1)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) + delta_alpha*k*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow(1 - x2,2))*std::exp(-time)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - delta_alpha*k*u_char*std::pow(x1,2)*(2*x1 - 2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - 2*delta_alpha*k*u_char*x1*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - u_char*std::pow(x1,2)*(2*x1 - 2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5) - 2*u_char*x1*std::pow(1 - x1,2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du121 = 2*std::pow(delta_alpha,2)*std::pow(k,2)*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),3) - delta_alpha*std::pow(k,2)*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x2)*std::cos(k*x1)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) + delta_alpha*k*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow(1 - x2,2))*std::exp(-time)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - delta_alpha*k*u_char*std::pow(x1,2)*(2*x1 - 2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - 2*delta_alpha*k*u_char*x1*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - u_char*std::pow(x1,2)*(2*x1 - 2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5) - 2*u_char*x1*std::pow(1 - x1,2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow(1 - x2,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du122 = -2*std::pow(delta_alpha,2)*std::pow(k,2)*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::pow(std::sin(k*x1),2)*std::pow(std::sin(k*x2),2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),3) - delta_alpha*std::pow(k,2)*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - 2*delta_alpha*k*u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow(1 - x2,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - u_char*std::pow(x1,2)*std::pow(1 - x1,2)*(12*u_char*x2 + 6*u_char*(2*x2 - 2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du21 = delta_alpha*k*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du22 = -delta_alpha*k*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - u_char*std::pow(x2,2)*(2*x2 - 2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5) - 2*u_char*x2*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du211 = -2*std::pow(delta_alpha,2)*std::pow(k,2)*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::pow(std::cos(k*x1),2)*std::pow(std::cos(k*x2),2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),3) - delta_alpha*std::pow(k,2)*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::sin(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) + 2*delta_alpha*k*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow(1 - x1,2))*std::exp(-time)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-12*u_char*x1 - 6*u_char*(2*x1 - 2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du212 = 2*std::pow(delta_alpha,2)*std::pow(k,2)*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),3) - delta_alpha*std::pow(k,2)*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::sin(k*x2)*std::cos(k*x1)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - delta_alpha*k*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow(1 - x1,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) + delta_alpha*k*u_char*std::pow(x2,2)*(2*x2 - 2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) + 2*delta_alpha*k*u_char*x2*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - u_char*std::pow(x2,2)*(2*x2 - 2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5) - 2*u_char*x2*std::pow(1 - x2,2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du221 = 2*std::pow(delta_alpha,2)*std::pow(k,2)*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),3) - delta_alpha*std::pow(k,2)*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::sin(k*x2)*std::cos(k*x1)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - delta_alpha*k*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow(1 - x1,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) + delta_alpha*k*u_char*std::pow(x2,2)*(2*x2 - 2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) + 2*delta_alpha*k*u_char*x2*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::cos(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - u_char*std::pow(x2,2)*(2*x2 - 2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5) - 2*u_char*x2*std::pow(1 - x2,2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        du222 = -2*std::pow(delta_alpha,2)*std::pow(k,2)*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::pow(std::sin(k*x1),2)*std::pow(std::sin(k*x2),2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),3) - delta_alpha*std::pow(k,2)*u_char*std::pow(x2,2)*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::sin(k*x1)*std::cos(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - 2*delta_alpha*k*u_char*std::pow(x2,2)*(2*x2 - 2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - 4*delta_alpha*k*u_char*x2*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)*std::sin(k*x1)*std::sin(k*x2)/std::pow((delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5),2) - 2*u_char*std::pow(x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5) - 4*u_char*x2*(2*x2 - 2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5) - 2*u_char*std::pow(1 - x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow(1 - x1,2))*std::exp(-time)/(delta_alpha*std::sin(k*x1)*std::cos(k*x2) + 0.5);

        r_pressure = 0.0;

        double det_permeability = MathUtils<double>::Det(permeability);
        MathUtils<double>::InvertMatrix(permeability, inv_permeability, det_permeability, -1.0);

        Matrix sigma = nu * rho * inv_permeability;

        const double convective1 = r_u1 * du11 + r_u2 * du12;
        const double convective2 = r_u1 * du21 + r_u2 * du22;

        const double div_of_sym_grad1 = (1.0/2.0) * (2 * du111 + du212 + du122);
        const double div_of_sym_grad2 = (1.0/2.0) * (du121 + du211 + 2 * du222);

        const double grad_of_div1 = du111 + du221;
        const double grad_of_div2 = du112 + du222;

        const double press_grad1 = 0.0;
        const double press_grad2 = 0.0;

        if (mAlternativeFormulation){
            const double grad_alpha_sym_grad1 = (1.0/2.0) * (2 * r_alpha1 * du11 + r_alpha2 * (du21 + du12));
            const double grad_alpha_sym_grad2 = (1.0/2.0) * (r_alpha1 * (du12 + du21) + 2 * r_alpha2 * du22);

            const double grad_alpha_div1 = r_alpha1 * (du11 + du22);
            const double grad_alpha_div2 = r_alpha2 * (du11 + du22);

            r_body_force1 = r_alpha * du1dt + r_alpha * convective1 + r_alpha / rho * press_grad1 - 2.0 * nu * (r_alpha * div_of_sym_grad1 + grad_alpha_sym_grad1) + (2.0/3.0) * nu * (r_alpha * grad_of_div1 + grad_alpha_div1) + sigma(0,0) * r_u1 + sigma(1,0) * r_u1;

            r_body_force2 = r_alpha * du2dt + r_alpha * convective2 + r_alpha / rho * press_grad2 - 2.0 * nu * (r_alpha * div_of_sym_grad2 + grad_alpha_sym_grad2) + (2.0/3.0) * nu * (r_alpha * grad_of_div2 + grad_alpha_div2) + sigma(0,1) * r_u2 + sigma(1,1) * r_u2;
        }else{
            r_body_force1 = du1dt + convective1 + 1.0/rho * press_grad1 - 2.0 * nu * div_of_sym_grad1 + (2.0/3.0) * nu * grad_of_div1 + sigma(0,0) * r_u1 + sigma(1,0) * r_u1;

            r_body_force2 = du2dt + convective2 + 1.0/rho * press_grad2 - 2.0 * nu * div_of_sym_grad2 + (2.0/3.0) * nu * grad_of_div2 + sigma(0,1) * r_u2 + sigma(1,1) * r_u2;
        }

        r_mass_source = (r_dalphat + r_u1 * r_alpha1 + r_u2 * r_alpha2 + r_alpha * (du11 + du22));

        if (mInitialConditions == true){
            if (mrModelPart.GetProcessInfo()[STEP] == 0)
            {
                it_node->FastGetSolutionStepValue(VELOCITY_X) = r_u1;
                it_node->FastGetSolutionStepValue(VELOCITY_Y) = r_u2;
                it_node->FastGetSolutionStepValue(PRESSURE) = r_pressure;
            }
        }
        else if(mInitialConditions == false && mrModelPart.GetProcessInfo()[STEP] == 1){
            it_node->FastGetSolutionStepValue(FLUID_FRACTION_OLD) = r_alpha;
        }
    }

}

void SinusoidalPorositySolutionTransientBodyForceProcess::SetFluidProperties()
{
    (mrModelPart.pGetProperties(1))->SetValue(DENSITY, mDensity);
    (mrModelPart.pGetProperties(1))->SetValue(DYNAMIC_VISCOSITY, mViscosity * mDensity);
    (mrModelPart.pGetProperties(1))->SetValue(VISCOSITY, mViscosity);

    block_for_each(mrModelPart.Elements(), [&](Element& rElement){
        rElement.SetProperties(mrModelPart.pGetProperties(1));
    });

    block_for_each(mrModelPart.Nodes(), [&](Node<3>& rNode){
        rNode.FastGetSolutionStepValue(VISCOSITY) = mViscosity;
        rNode.FastGetSolutionStepValue(DENSITY) = mDensity;
        rNode.FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = mViscosity * mDensity;
    });
}

/* Private functions ****************************************************/

};  // namespace Kratos.