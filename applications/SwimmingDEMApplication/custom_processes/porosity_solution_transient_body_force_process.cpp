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
#include "porosity_solution_transient_body_force_process.h"
#include "swimming_dem_application_variables.h"

// Other applications includes
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

/* Public functions *******************************************************/
PorositySolutionTransientBodyForceProcess::PorositySolutionTransientBodyForceProcess(
    ModelPart& rModelPart)
    : Process(),
      mrModelPart(rModelPart)
{}

PorositySolutionTransientBodyForceProcess::PorositySolutionTransientBodyForceProcess(
    ModelPart& rModelPart,
    Parameters& rParameters)
    : Process(),
      mrModelPart(rModelPart)
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);
}

PorositySolutionTransientBodyForceProcess::PorositySolutionTransientBodyForceProcess(
    Model &rModel,
    Parameters &rParameters)
    : Process(),
      mrModelPart(rModel.GetModelPart(rParameters["model_part_name"].GetString()))
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);

}


void PorositySolutionTransientBodyForceProcess::CheckDefaultsAndProcessSettings(Parameters &rParameters)
{
    const Parameters default_parameters = GetDefaultParameters();

    rParameters.ValidateAndAssignDefaults(default_parameters);

    mDensity     = rParameters["benchmark_parameters"]["density"].GetDouble();
    mUchar       = rParameters["benchmark_parameters"]["u_char"].GetDouble();
    mAlpha       = rParameters["benchmark_parameters"]["alpha"].GetDouble();
    mViscosity   = rParameters["benchmark_parameters"]["viscosity"].GetDouble();
    mInitialConditions = rParameters["benchmark_parameters"]["use_initial_conditions"].GetBool();
    mAlternativeFormulation = rParameters["benchmark_parameters"]["use_alternative_formulation"].GetBool();

}


const Parameters PorositySolutionTransientBodyForceProcess::GetDefaultParameters() const
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
                                                "alpha"       : 1.0,
                                                "u_char"      : 100.0,
                                                "use_alternative_formulation" : false
                },
                "compute_nodal_error"      : true,
                "print_convergence_output" : false,
                "output_parameters"        : {}
    }  )" );

    return default_parameters;
}

void PorositySolutionTransientBodyForceProcess::Execute()
{
    this->ExecuteInitialize();
    this->ExecuteInitializeSolutionStep();
}

void PorositySolutionTransientBodyForceProcess::ExecuteInitialize()
{}

void PorositySolutionTransientBodyForceProcess::ExecuteBeforeSolutionLoop()
{
    this->SetFluidProperties();
    if (mInitialConditions == true)
    {
        this->SetInitialBodyForceAndPorosityField();
    }
}

void PorositySolutionTransientBodyForceProcess::ExecuteInitializeSolutionStep()
{
    this->SetBodyForceAndPorosityField();
}

void PorositySolutionTransientBodyForceProcess::ExecuteFinalizeSolutionStep() {}

/* Protected functions ****************************************************/

void PorositySolutionTransientBodyForceProcess::SetInitialBodyForceAndPorosityField()
{
    const double time = mrModelPart.GetProcessInfo()[TIME];
    const double delta_time = mrModelPart.GetProcessInfo()[DELTA_TIME];
    const double dim = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
    const double rho = mDensity;
    const double nu = mViscosity;
    const double u_char = mUchar;
    const double alpha = mAlpha;
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

        double& r_u1 = it_node->FastGetSolutionStepValue(EXACT_VELOCITY)[0];
        double& r_u2 = it_node->FastGetSolutionStepValue(EXACT_VELOCITY)[1];

        double& r_pressure = it_node->FastGetSolutionStepValue(EXACT_PRESSURE);

        Matrix& permeability = it_node->FastGetSolutionStepValue(PERMEABILITY);

        r_u1 = u_char*std::pow(x1,2)*std::pow((1 - x1),2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        double u1_1 = u_char*std::pow(x1,2)*std::pow((1 - x1),2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow((1 - x2),2))*std::exp(-(time-delta_time))*std::cos(Globals::Pi*(time-delta_time))/alpha;

        double u1_2 = u_char*std::pow(x1,2)*std::pow((1 - x1),2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow((1 - x2),2))*std::exp(-(time-2*delta_time))*std::cos(Globals::Pi*(time-2*delta_time))/alpha;

        r_u2 = u_char*std::pow(x2,2)*std::pow((1 - x2),2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        double u2_1 = u_char*std::pow(x2,2)*std::pow((1 - x2),2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow((1 - x1),2))*std::exp(-(time-delta_time))*std::cos(Globals::Pi*(time-delta_time))/alpha;

        double u2_2 = u_char*std::pow(x2,2)*std::pow((1 - x2),2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow((1 - x1),2))*std::exp(-(time-2*delta_time))*std::cos(Globals::Pi*(time-2*delta_time))/alpha;

        du1dt = -Globals::Pi*u_char*std::pow(x1,2)*std::pow((1 - x1),2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow((1 - x2),2))*std::exp(-time)*std::sin(Globals::Pi*time)/alpha - u_char*std::pow(x1,2)*std::pow((1 - x1),2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du2dt = -Globals::Pi*u_char*std::pow(x2,2)*std::pow((1 - x2),2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow((1 - x1),2))*std::exp(-time)*std::sin(Globals::Pi*time)/alpha - u_char*std::pow(x2,2)*std::pow((1 - x2),2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        r_alpha = alpha;

        r_dalphat = 0.0;

        r_alpha1 = 0.0;

        r_alpha2 = 0.0;

        r_pressure = 0.0;

        du11 = u_char*std::pow(x1,2)*(2*x1 - 2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha + 2*u_char*x1*std::pow((1 - x1),2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du12 = u_char*std::pow(x1,2)*std::pow((1 - x1),2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du111 =  2*u_char*std::pow(x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha + 4*u_char*x1*(2*x1 - 2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha + 2*u_char*std::pow((1 - x1),2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du112 = u_char*std::pow(x1,2)*(2*x1 - 2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha + 2*u_char*x1*std::pow((1 - x1),2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du121 = u_char*std::pow(x1,2)*(2*x1 - 2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha + 2*u_char*x1*std::pow((1 - x1),2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du122 = u_char*std::pow(x1,2)*std::pow((1 - x1),2)*(12*u_char*x2 + 6*u_char*(2*x2 - 2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du21 = u_char*std::pow(x2,2)*std::pow((1 - x2),2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du22 = u_char*std::pow(x2,2)*(2*x2 - 2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha + 2*u_char*x2*std::pow((1 - x2),2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du211 = u_char*std::pow(x2,2)*std::pow((1 - x2),2)*(-12*u_char*x1 - 6*u_char*(2*x1 - 2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du212 = u_char*std::pow(x2,2)*(2*x2 - 2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha + 2*u_char*x2*std::pow((1 - x2),2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du221 = u_char*std::pow(x2,2)*(2*x2 - 2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha + 2*u_char*x2*std::pow((1 - x2),2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du222 = 2*u_char*std::pow(x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha + 4*u_char*x2*(2*x2 - 2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha + 2*u_char*std::pow((1 - x2),2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        for (unsigned int d = 0; d < dim; ++d){
            permeability(d,d) = 1.0e+30;
        }


        double det_permeability = MathUtils<double>::Det(permeability);
        MathUtils<double>::InvertMatrix(permeability, inv_permeability, det_permeability, -1.0);

        Matrix sigma = nu * rho * inv_permeability;

        const double convective1 = r_u1 * du11 + r_u2 * du12;
        const double convective2 = r_u1 * du21 + r_u2 * du22;

        const double div_of_sym_grad1 = (1.0/2.0) * (2.0 * du111 + du212 + du122);
        const double div_of_sym_grad2 = (1.0/2.0) * (du121 + du211 + 2.0 * du222);

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

        r_mass_source = r_dalphat + r_u1 * r_alpha1 + r_u2 * r_alpha2 + r_alpha * (du11 + du22);

        it_node->FastGetSolutionStepValue(VELOCITY_X) = r_u1;
        it_node->FastGetSolutionStepValue(VELOCITY_Y) = r_u2;
        it_node->FastGetSolutionStepValue(VELOCITY_X,1) = u1_1;
        it_node->FastGetSolutionStepValue(VELOCITY_Y,1) = u2_1;
        it_node->FastGetSolutionStepValue(VELOCITY_X,2) = u1_2;
        it_node->FastGetSolutionStepValue(VELOCITY_Y,2) = u2_2;
        it_node->FastGetSolutionStepValue(PRESSURE) = r_pressure;
        it_node->FastGetSolutionStepValue(FLUID_FRACTION_OLD) = r_alpha;
    }

}

void PorositySolutionTransientBodyForceProcess::SetBodyForceAndPorosityField()
{
    const double time = mrModelPart.GetProcessInfo()[TIME];
    const double dim = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
    const double rho = mDensity;
    const double nu = mViscosity;
    const double u_char = mUchar;
    const double alpha = mAlpha;
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

        double& r_u1 = it_node->FastGetSolutionStepValue(EXACT_VELOCITY)[0];
        double& r_u2 = it_node->FastGetSolutionStepValue(EXACT_VELOCITY)[1];

        double& r_pressure = it_node->FastGetSolutionStepValue(EXACT_PRESSURE);

        Matrix& permeability = it_node->FastGetSolutionStepValue(PERMEABILITY);

        r_u1 = u_char*std::pow(x1,2)*std::pow((1 - x1),2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        r_u2 = u_char*std::pow(x2,2)*std::pow((1 - x2),2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du1dt = -Globals::Pi*u_char*std::pow(x1,2)*std::pow((1 - x1),2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow((1 - x2),2))*std::exp(-time)*std::sin(Globals::Pi*time)/alpha - u_char*std::pow(x1,2)*std::pow((1 - x1),2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du2dt = -Globals::Pi*u_char*std::pow(x2,2)*std::pow((1 - x2),2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow((1 - x1),2))*std::exp(-time)*std::sin(Globals::Pi*time)/alpha - u_char*std::pow(x2,2)*std::pow((1 - x2),2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        r_alpha = alpha;

        r_dalphat = 0.0;

        r_alpha1 = 0.0;

        r_alpha2 = 0.0;

        r_pressure = 0.0;

        du11 = u_char*std::pow(x1,2)*(2*x1 - 2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha + 2*u_char*x1*std::pow((1 - x1),2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du12 = u_char*std::pow(x1,2)*std::pow((1 - x1),2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du111 =  2*u_char*std::pow(x1,2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha + 4*u_char*x1*(2*x1 - 2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha + 2*u_char*std::pow((1 - x1),2)*(u_char*std::pow(x2,2)*(2*x2 - 2) + 2*u_char*x2*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du112 = u_char*std::pow(x1,2)*(2*x1 - 2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha + 2*u_char*x1*std::pow((1 - x1),2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du121 = u_char*std::pow(x1,2)*(2*x1 - 2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha + 2*u_char*x1*std::pow((1 - x1),2)*(2*u_char*std::pow(x2,2) + 4*u_char*x2*(2*x2 - 2) + 2*u_char*std::pow((1 - x2),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du122 = u_char*std::pow(x1,2)*std::pow((1 - x1),2)*(12*u_char*x2 + 6*u_char*(2*x2 - 2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du21 = u_char*std::pow(x2,2)*std::pow((1 - x2),2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du22 = u_char*std::pow(x2,2)*(2*x2 - 2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha + 2*u_char*x2*std::pow((1 - x2),2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du211 = u_char*std::pow(x2,2)*std::pow((1 - x2),2)*(-12*u_char*x1 - 6*u_char*(2*x1 - 2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du212 = u_char*std::pow(x2,2)*(2*x2 - 2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha + 2*u_char*x2*std::pow((1 - x2),2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du221 = u_char*std::pow(x2,2)*(2*x2 - 2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha + 2*u_char*x2*std::pow((1 - x2),2)*(-2*u_char*std::pow(x1,2) - 4*u_char*x1*(2*x1 - 2) - 2*u_char*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        du222 = 2*u_char*std::pow(x2,2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha + 4*u_char*x2*(2*x2 - 2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha + 2*u_char*std::pow((1 - x2),2)*(-u_char*std::pow(x1,2)*(2*x1 - 2) - 2*u_char*x1*std::pow((1 - x1),2))*std::exp(-time)*std::cos(Globals::Pi*time)/alpha;

        for (unsigned int d = 0; d < dim; ++d){
            permeability(d,d) = 1.0e+30;
        };


        double det_permeability = MathUtils<double>::Det(permeability);
        MathUtils<double>::InvertMatrix(permeability, inv_permeability, det_permeability, -1.0);

        Matrix sigma = nu * rho * inv_permeability;

        const double convective1 = r_u1 * du11 + r_u2 * du12;
        const double convective2 = r_u1 * du21 + r_u2 * du22;

        const double div_of_sym_grad1 = (1.0/2.0) * (2.0 * du111 + du212 + du122);
        const double div_of_sym_grad2 = (1.0/2.0) * (du121 + du211 + 2.0 * du222);

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

        r_mass_source = r_dalphat + r_u1 * r_alpha1 + r_u2 * r_alpha2 + r_alpha * (du11 + du22);

    }

}

void PorositySolutionTransientBodyForceProcess::SetFluidProperties()
{
    (mrModelPart.pGetProperties(1))->SetValue(DENSITY, mDensity);
    (mrModelPart.pGetProperties(1))->SetValue(DYNAMIC_VISCOSITY, mViscosity * mDensity);
    (mrModelPart.pGetProperties(1))->SetValue(VISCOSITY, mViscosity);

    block_for_each(mrModelPart.Elements(), [&](Element& rElement){
        rElement.SetProperties(mrModelPart.pGetProperties(1));
    });

    block_for_each(mrModelPart.Nodes(), [&](Node& rNode){
        rNode.FastGetSolutionStepValue(VISCOSITY) = mViscosity;
        rNode.FastGetSolutionStepValue(DENSITY) = mDensity;
        rNode.FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = mViscosity * mDensity;
    });
}

/* Private functions ****************************************************/

};  // namespace Kratos.