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
#include "skrzypacz_porosity_solution_and_body_force_process.h"
#include "swimming_dem_application_variables.h"

// Other applications includes
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

/* Public functions *******************************************************/
SkrzypaczPorositySolutionAndBodyForceProcess::SkrzypaczPorositySolutionAndBodyForceProcess(
    ModelPart& rModelPart)
    : Process(),
      mrModelPart(rModelPart)
{}

SkrzypaczPorositySolutionAndBodyForceProcess::SkrzypaczPorositySolutionAndBodyForceProcess(
    ModelPart& rModelPart,
    Parameters& rParameters)
    : Process(),
      mrModelPart(rModelPart)
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);
}

SkrzypaczPorositySolutionAndBodyForceProcess::SkrzypaczPorositySolutionAndBodyForceProcess(
    Model &rModel,
    Parameters &rParameters)
    : Process(),
      mrModelPart(rModel.GetModelPart(rParameters["model_part_name"].GetString()))
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);

}


void SkrzypaczPorositySolutionAndBodyForceProcess::CheckDefaultsAndProcessSettings(Parameters &rParameters)
{
    const Parameters default_parameters = GetDefaultParameters();

    rParameters.ValidateAndAssignDefaults(default_parameters);

    mDensity     = rParameters["benchmark_parameters"]["density"].GetDouble();
    mViscosity   = rParameters["benchmark_parameters"]["viscosity"].GetDouble();
    mInitialConditions = rParameters["benchmark_parameters"]["use_initial_conditions"].GetBool();
    mAlternativeFormulation = rParameters["benchmark_parameters"]["use_alternative_formulation"].GetBool();

}

const Parameters SkrzypaczPorositySolutionAndBodyForceProcess::GetDefaultParameters() const
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
                                                "use_initial_conditions" : true,
                                                "use_alternative_formulation" : false
                },
                "compute_nodal_error"      : true,
                "print_convergence_output" : false,
                "output_parameters"        : {}
    }  )" );

    return default_parameters;
}

void SkrzypaczPorositySolutionAndBodyForceProcess::Execute()
{
    this->ExecuteInitialize();
    this->ExecuteInitializeSolutionStep();
}

void SkrzypaczPorositySolutionAndBodyForceProcess::ExecuteInitialize()
{}

void SkrzypaczPorositySolutionAndBodyForceProcess::ExecuteBeforeSolutionLoop()
{
    //this->SetFluidProperties();

    if (mInitialConditions == true)
    {
        this->SetInitialBodyForceAndPorosityField();
        this->SetValuesOnIntegrationPoints();
    }

}

void SkrzypaczPorositySolutionAndBodyForceProcess::ExecuteInitializeSolutionStep()
{}

void SkrzypaczPorositySolutionAndBodyForceProcess::ExecuteFinalizeSolutionStep() {}

/* Protected functions ****************************************************/

void SkrzypaczPorositySolutionAndBodyForceProcess::SetInitialBodyForceAndPorosityField()
{
    const double Dim = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
    const double rho = mDensity;
    const double nu = mViscosity;
    Matrix I = IdentityMatrix(Dim, Dim);

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

        double r_u1;
        double r_u2;

        double& r_U1 = it_node->FastGetSolutionStepValue(VELOCITY_X);
        double& r_U2 = it_node->FastGetSolutionStepValue(VELOCITY_Y);

        double& r_pressure = it_node->FastGetSolutionStepValue(EXACT_PRESSURE);

        Matrix& r_sigma = it_node->FastGetSolutionStepValue(PERMEABILITY);

        r_u1 = std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1);

        r_u2 = std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1);

        r_U1 = r_u1;

        r_U2 = r_u2;

        du1dt = 0.0;

        du2dt = 0.0;

        r_alpha = -0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1;

        r_dalphat = 0.0;

        r_alpha1 = -0.5*Globals::Pi*std::sin(Globals::Pi*x2)*std::cos(Globals::Pi*x1);

        r_alpha2 = -0.5*Globals::Pi*std::sin(Globals::Pi*x1)*std::cos(Globals::Pi*x2);

        du11 = Globals::Pi*std::sin(Globals::Pi*x2)*std::cos(Globals::Pi*x1)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) + 0.5*Globals::Pi*std::sin(Globals::Pi*x1)*std::pow(std::sin(Globals::Pi*x2),2)*std::cos(Globals::Pi*x1)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2);

        du12 = Globals::Pi*std::sin(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) + 0.5*Globals::Pi*std::pow(std::sin(Globals::Pi*x1),2)*std::sin(Globals::Pi*x2)*std::cos(Globals::Pi*x2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2);

        du111 = -std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) - 0.5*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x1),2)*std::pow(std::sin(Globals::Pi*x2),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 1.0*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x2),2)*std::pow(std::cos(Globals::Pi*x1),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 0.5*std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1)*std::pow(std::sin(Globals::Pi*x2),3)*std::pow(std::cos(Globals::Pi*x1),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),3);

        du112 = std::pow(Globals::Pi,2)*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) + 1.5*std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 0.5*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x1),2)*std::pow(std::sin(Globals::Pi*x2),2)*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),3);

        du121 = std::pow(Globals::Pi,2)*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) + 1.5*std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 0.5*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x1),2)*std::pow(std::sin(Globals::Pi*x2),2)*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),3);

        du122 = -std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) - 0.5*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x1),2)*std::pow(std::sin(Globals::Pi*x2),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 1.0*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x1),2)*std::pow(std::cos(Globals::Pi*x2),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 0.5*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x1),3)*std::sin(Globals::Pi*x2)*std::pow(std::cos(Globals::Pi*x2),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),3);

        du21 = -Globals::Pi*std::sin(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) + 0.5*Globals::Pi*std::sin(Globals::Pi*x2)*std::pow(std::cos(Globals::Pi*x1),2)*std::cos(Globals::Pi*x2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2);

        du22 = -Globals::Pi*std::sin(Globals::Pi*x2)*std::cos(Globals::Pi*x1)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) + 0.5*Globals::Pi*std::sin(Globals::Pi*x1)*std::cos(Globals::Pi*x1)*std::pow(std::cos(Globals::Pi*x2),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2);

        du211 = -std::pow(Globals::Pi,2)*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) - 1.5*std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 0.5*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x2),2)*std::pow(std::cos(Globals::Pi*x1),3)*std::cos(Globals::Pi*x2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),3);

        du212 = std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) - 0.5*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x1),2)*std::pow(std::cos(Globals::Pi*x2),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) - 0.5*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x2),2)*std::pow(std::cos(Globals::Pi*x1),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 0.5*std::pow(Globals::Pi,2)*std::pow(std::cos(Globals::Pi*x1),2)*std::pow(std::cos(Globals::Pi*x2),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 0.5*std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)*std::pow(std::cos(Globals::Pi*x1),2)*std::pow(std::cos(Globals::Pi*x2),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),3);

        du221 = std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) - 0.5*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x1),2)*std::pow(std::cos(Globals::Pi*x2),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) - 0.5*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x2),2)*std::pow(std::cos(Globals::Pi*x1),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 0.5*std::pow(Globals::Pi,2)*std::pow(std::cos(Globals::Pi*x1),2)*std::pow(std::cos(Globals::Pi*x2),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 0.5*std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)*std::pow(std::cos(Globals::Pi*x1),2)*std::pow(std::cos(Globals::Pi*x2),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),3);

        du222 = -std::pow(Globals::Pi,2)*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) - 1.5*std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 0.5*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x1),2)*std::cos(Globals::Pi*x1)*std::pow(std::cos(Globals::Pi*x2),3)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),3);

        r_pressure = 2*std::sin(Globals::Pi*x2)*std::cos(Globals::Pi*x1);

        it_node->FastGetSolutionStepValue(VELOCITY_X,1) = r_u1;
        it_node->FastGetSolutionStepValue(VELOCITY_Y,1) = r_u2;

        const double velocity_norm = std::sqrt(r_u1 * r_u1 + r_u2 * r_u2);
        const double kappa = (1 - r_alpha)/r_alpha;

        r_sigma = (nu * 150.0 * std::pow(kappa,2) + 1.75 * kappa * velocity_norm) * I;

        const double convective1 = r_u1 * du11 + r_u2 * du12;
        const double convective2 = r_u1 * du21 + r_u2 * du22;

        const double div_of_sym_grad1 = (1.0/2.0) * (2.0 * du111 + du212 + du122);
        const double div_of_sym_grad2 = (1.0/2.0) * (du121 + du211 + 2.0 * du222);

        const double grad_of_div1 = du111 + du221;
        const double grad_of_div2 = du112 + du222;

        const double press_grad1 = -2*Globals::Pi*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2);
        const double press_grad2 = 2*Globals::Pi*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2);

        if (mAlternativeFormulation){
            const double grad_alpha_sym_grad1 = (1.0/2.0) * (2 * r_alpha1 * du11 + r_alpha2 * (du21 + du12));
            const double grad_alpha_sym_grad2 = (1.0/2.0) * (r_alpha1 * (du12 + du21) + 2 * r_alpha2 * du22);

            const double grad_alpha_div1 = r_alpha1 * (du11 + du22);
            const double grad_alpha_div2 = r_alpha2 * (du11 + du22);

            r_body_force1 = r_alpha * du1dt + r_alpha * convective1 + r_alpha / rho * press_grad1 - 2.0 * nu * (r_alpha * div_of_sym_grad1 + grad_alpha_sym_grad1) + (2.0/3.0) * nu * (r_alpha * grad_of_div1 + grad_alpha_div1) + r_sigma(0,0) * r_u1 + r_sigma(0,1) * r_u2;

            r_body_force2 = r_alpha * du2dt + r_alpha * convective2 + r_alpha / rho * press_grad2 - 2.0 * nu * (r_alpha * div_of_sym_grad2 + grad_alpha_sym_grad2) + (2.0/3.0) * nu * (r_alpha * grad_of_div2 + grad_alpha_div2) + r_sigma(1,0) * r_u1 + r_sigma(1,1) * r_u2;

        }else{
            r_body_force1 = du1dt + convective1 + 1.0/rho * press_grad1 - 2.0 * nu * div_of_sym_grad1 + (2.0/3.0) * nu * grad_of_div1 + r_sigma(0,0) * r_u1 + r_sigma(0,1) * r_u2;

            r_body_force2 = du2dt + convective2 + 1.0/rho * press_grad2 - 2.0 * nu * div_of_sym_grad2 + (2.0/3.0) * nu * grad_of_div2 + r_sigma(1,0) * r_u1 + r_sigma(1,1) * r_u2;
        }

        r_mass_source = r_dalphat + r_u1 * r_alpha1 + r_u2 * r_alpha2 + r_alpha * (du11 + du22);

        it_node->SetLock();
        it_node->FastGetSolutionStepValue(EXACT_VELOCITY) = ZeroVector(3);
        it_node->FastGetSolutionStepValue(EXACT_VELOCITY)[0] = r_u1;
        it_node->FastGetSolutionStepValue(EXACT_VELOCITY)[1] = r_u2;
        it_node->UnSetLock();

    }

}

void SkrzypaczPorositySolutionAndBodyForceProcess::SetValuesOnIntegrationPoints()
{
    const unsigned int Dim = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
    const double rho = mDensity;
    const double nu = mViscosity;
    Matrix I = IdentityMatrix(Dim, Dim);
    Matrix sigma = ZeroMatrix(Dim, Dim);
    double u1, u2, pressure,alpha;
    double du1dt, du2dt, du11, du12, du21, du22, du111, du112, du121, du122, du211, du212, du221, du222, body_force1, body_force2;

    const unsigned int n_elem = mrModelPart.NumberOfElements();

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
        std::vector<array_1d<double,3>> body_force_on_gauss_points(r_number_integration_points);
        std::vector<array_1d<double,3>> pressure_gradient_on_gauss_points(r_number_integration_points);
        std::vector<array_1d<double,3>> fluid_fraction_gradient_on_gauss_points(r_number_integration_points);
        std::vector<Matrix> velocity_gradient_on_gauss_points(r_number_integration_points);
        std::vector<Matrix> resistance_tensor_on_gauss_points(r_number_integration_points);

        // if (mExactScalar.size() != n_elem)
        //     mExactScalar.resize(n_elem);

        // if (mExactVector.size() != n_elem)
        //     mExactVector.resize(n_elem);

        // if (mExactScalarGradient.size() != n_elem)
        //     mExactScalarGradient.resize(n_elem);

        // if (mExactPorosity.size() != n_elem)
        //     mExactPorosity.resize(n_elem);

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

        for (unsigned int g = 0; g < r_number_integration_points; g++){

            Matrix gauss_point_coordinates = ZeroMatrix(r_number_integration_points,Dim);
            for (unsigned int i = 0; i < NumNodes; ++i){
                const array_1d<double, 3>& r_coordinates = r_geometry[i].Coordinates();
                for (unsigned int d = 0; d < Dim; ++d)
                    gauss_point_coordinates(g,d) += NContainer(g,i) * r_coordinates[d];
                }

            velocity_gradient_on_gauss_points[g] = ZeroMatrix(3,3);
            resistance_tensor_on_gauss_points[g] = ZeroMatrix(3,3);
            pressure_gradient_on_gauss_points[g] = ZeroVector(3);
            fluid_fraction_gradient_on_gauss_points[g] = ZeroVector(3);
            body_force_on_gauss_points[g] = ZeroVector(3);
            velocity_on_gauss_points[g] = ZeroVector(3);

            const double x1 = gauss_point_coordinates(g,0);
            const double x2 = gauss_point_coordinates(g,1);

            double alpha1;
            double alpha2;

            u1 = std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1);

            u2 = std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1);

            du1dt = 0.0;

            du2dt = 0.0;

            alpha = -0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1;

            alpha1 = -0.5*Globals::Pi*std::sin(Globals::Pi*x2)*std::cos(Globals::Pi*x1);

            alpha2 = -0.5*Globals::Pi*std::sin(Globals::Pi*x1)*std::cos(Globals::Pi*x2);

            du11 = Globals::Pi*std::sin(Globals::Pi*x2)*std::cos(Globals::Pi*x1)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) + 0.5*Globals::Pi*std::sin(Globals::Pi*x1)*std::pow(std::sin(Globals::Pi*x2),2)*std::cos(Globals::Pi*x1)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2);

            du12 = Globals::Pi*std::sin(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) + 0.5*Globals::Pi*std::pow(std::sin(Globals::Pi*x1),2)*std::sin(Globals::Pi*x2)*std::cos(Globals::Pi*x2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2);

            du111 = -std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) - 0.5*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x1),2)*std::pow(std::sin(Globals::Pi*x2),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 1.0*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x2),2)*std::pow(std::cos(Globals::Pi*x1),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 0.5*std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1)*std::pow(std::sin(Globals::Pi*x2),3)*std::pow(std::cos(Globals::Pi*x1),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),3);

            du112 = std::pow(Globals::Pi,2)*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) + 1.5*std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 0.5*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x1),2)*std::pow(std::sin(Globals::Pi*x2),2)*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),3);

            du121 = std::pow(Globals::Pi,2)*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) + 1.5*std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 0.5*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x1),2)*std::pow(std::sin(Globals::Pi*x2),2)*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),3);

            du122 = -std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) - 0.5*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x1),2)*std::pow(std::sin(Globals::Pi*x2),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 1.0*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x1),2)*std::pow(std::cos(Globals::Pi*x2),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 0.5*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x1),3)*std::sin(Globals::Pi*x2)*std::pow(std::cos(Globals::Pi*x2),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),3);

            du21 = -Globals::Pi*std::sin(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) + 0.5*Globals::Pi*std::sin(Globals::Pi*x2)*std::pow(std::cos(Globals::Pi*x1),2)*std::cos(Globals::Pi*x2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2);

            du22 = -Globals::Pi*std::sin(Globals::Pi*x2)*std::cos(Globals::Pi*x1)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) + 0.5*Globals::Pi*std::sin(Globals::Pi*x1)*std::cos(Globals::Pi*x1)*std::pow(std::cos(Globals::Pi*x2),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2);

            du211 = -std::pow(Globals::Pi,2)*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) - 1.5*std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 0.5*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x2),2)*std::pow(std::cos(Globals::Pi*x1),3)*std::cos(Globals::Pi*x2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),3);

            du212 = std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) - 0.5*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x1),2)*std::pow(std::cos(Globals::Pi*x2),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) - 0.5*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x2),2)*std::pow(std::cos(Globals::Pi*x1),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 0.5*std::pow(Globals::Pi,2)*std::pow(std::cos(Globals::Pi*x1),2)*std::pow(std::cos(Globals::Pi*x2),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 0.5*std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)*std::pow(std::cos(Globals::Pi*x1),2)*std::pow(std::cos(Globals::Pi*x2),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),3);

            du221 = std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) - 0.5*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x1),2)*std::pow(std::cos(Globals::Pi*x2),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) - 0.5*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x2),2)*std::pow(std::cos(Globals::Pi*x1),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 0.5*std::pow(Globals::Pi,2)*std::pow(std::cos(Globals::Pi*x1),2)*std::pow(std::cos(Globals::Pi*x2),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 0.5*std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)*std::pow(std::cos(Globals::Pi*x1),2)*std::pow(std::cos(Globals::Pi*x2),2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),3);

            du222 = -std::pow(Globals::Pi,2)*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/(-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1) - 1.5*std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2)*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),2) + 0.5*std::pow(Globals::Pi,2)*std::pow(std::sin(Globals::Pi*x1),2)*std::cos(Globals::Pi*x1)*std::pow(std::cos(Globals::Pi*x2),3)/std::pow((-0.5*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2) + 1),3);

            pressure = 2*std::sin(Globals::Pi*x2)*std::cos(Globals::Pi*x1);

            const double press_grad1 = -2*std::sin(Globals::Pi*x2)*Globals::Pi*std::sin(Globals::Pi*x1);
            const double press_grad2 = 2*std::cos(Globals::Pi*x2)*Globals::Pi*std::cos(Globals::Pi*x1);

            const double velocity_norm = std::sqrt(u1 * u1 + u2 * u2);
            const double kappa = (1 - alpha)/alpha;

            sigma = (nu * 150.0 * std::pow(kappa,2) + 1.75 * kappa * velocity_norm) * I;

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
            fluid_fraction_gradient_on_gauss_points[g][0] = alpha1;
            fluid_fraction_gradient_on_gauss_points[g][1] = alpha2;
            velocity_on_gauss_points[g][0] = u1;
            velocity_on_gauss_points[g][1] = u2;
            body_force_on_gauss_points[g][0] = body_force1;
            body_force_on_gauss_points[g][1] = body_force2;
            pressure_gradient_on_gauss_points[g][0] = press_grad1;
            pressure_gradient_on_gauss_points[g][1] = press_grad2;
            velocity_gradient_on_gauss_points[g](0,0) = du11;
            velocity_gradient_on_gauss_points[g](0,1) = du21;
            velocity_gradient_on_gauss_points[g](1,0) = du12;
            velocity_gradient_on_gauss_points[g](1,1) = du22;
            resistance_tensor_on_gauss_points[g](0,0) = sigma(0,0);
            resistance_tensor_on_gauss_points[g](0,1) = sigma(0,1);
            resistance_tensor_on_gauss_points[g](1,0) = sigma(1,0);
            resistance_tensor_on_gauss_points[g](1,1) = sigma(1,1);
        }

        it_elem->SetValuesOnIntegrationPoints(EXACT_VELOCITY, velocity_on_gauss_points, mrModelPart.GetProcessInfo());
        it_elem->SetValuesOnIntegrationPoints(EXACT_PRESSURE, pressure_on_gauss_points, mrModelPart.GetProcessInfo());
        it_elem->SetValuesOnIntegrationPoints(RECOVERED_PRESSURE_GRADIENT, pressure_gradient_on_gauss_points, mrModelPart.GetProcessInfo());
        it_elem->SetValuesOnIntegrationPoints(EXACT_VELOCITY_GRADIENT,velocity_gradient_on_gauss_points, mrModelPart.GetProcessInfo());
        it_elem->SetValuesOnIntegrationPoints(FLUID_FRACTION_GRADIENT,fluid_fraction_gradient_on_gauss_points, mrModelPart.GetProcessInfo());
        it_elem->SetValuesOnIntegrationPoints(FLUID_FRACTION,fluid_fraction_on_gauss_points, mrModelPart.GetProcessInfo());
        it_elem->SetValuesOnIntegrationPoints(BODY_FORCE,body_force_on_gauss_points, mrModelPart.GetProcessInfo());
        it_elem->SetValuesOnIntegrationPoints(PERMEABILITY,resistance_tensor_on_gauss_points, mrModelPart.GetProcessInfo());

    }
}

void SkrzypaczPorositySolutionAndBodyForceProcess::SetFluidProperties()
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
