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
#include <iostream>
#include <fstream>
// External includes

// Project includes
#include "containers/model.h"
#include "includes/checks.h"
#include "utilities/variable_utils.h"
#include "utilities/math_utils.h"

// Application includes
#include "swimming_DEM_application.h"
#include "taylor_green_vortex_porosity_field_process.h"
#include "swimming_dem_application_variables.h"

// Other applications includes
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

/* Public functions *******************************************************/
TaylorGreenVortexPorosityFieldProcess::TaylorGreenVortexPorosityFieldProcess(
    ModelPart& rModelPart)
    : Process(),
      mrModelPart(rModelPart)
{}

TaylorGreenVortexPorosityFieldProcess::TaylorGreenVortexPorosityFieldProcess(
    ModelPart& rModelPart,
    Parameters& rParameters)
    : Process(),
      mrModelPart(rModelPart)
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);
}

TaylorGreenVortexPorosityFieldProcess::TaylorGreenVortexPorosityFieldProcess(
    Model &rModel,
    Parameters &rParameters)
    : Process(),
      mrModelPart(rModel.GetModelPart(rParameters["model_part_name"].GetString()))
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);

}


void TaylorGreenVortexPorosityFieldProcess::CheckDefaultsAndProcessSettings(Parameters &rParameters)
{
    const Parameters default_parameters = GetDefaultParameters();

    rParameters.ValidateAndAssignDefaults(default_parameters);

    mMeanDeviation  = rParameters["benchmark_parameters"]["mean_deviation"].GetDouble();
    mAlphaMax    = rParameters["benchmark_parameters"]["alpha_max"].GetDouble();
    mWaveNumber = rParameters["benchmark_parameters"]["wave_number"].GetInt();
    mAlternativeFormulation = rParameters["benchmark_parameters"]["use_alternative_formulation"].GetBool();
}

const Parameters TaylorGreenVortexPorosityFieldProcess::GetDefaultParameters() const
{
    const Parameters default_parameters( R"(
    {
                "model_part_name"          : "please_specify_model_part_name",
                "variable_name"            : "BODY_FORCE",
                "benchmark_name"           : "custom_body_force.vortex",
                "benchmark_parameters"     : {
                                                "mean_deviation"   : 0.5,
                                                "alpha_max"   : 1.0,
                                                "wave_number"   : 1
                },
                "compute_nodal_error"      : true,
                "print_convergence_output" : false,
                "output_parameters"        : {}
    }  )" );

    return default_parameters;
}

void TaylorGreenVortexPorosityFieldProcess::Execute(){}

array_1d<double,3> TaylorGreenVortexPorosityFieldProcess::ExecuteInTimeStep(){

    std::vector<array_1d<double,3>> velocity_on_gauss_points(16);
    std::vector<double> fluid_fraction_on_gauss_points(16);
    std::vector<array_1d<double,3>> vorticity_on_gauss_points(16);
    array_1d<double,3> data;
    double density = 0.0, viscosity = 0.0;

    double kinetic_energy = 0.0, dissipation_rate = 0.0, total_volume = 0.0;
    for (auto it_node = mrModelPart.NodesBegin(); it_node != mrModelPart.NodesEnd(); it_node++){
        viscosity = it_node->FastGetSolutionStepValue(VISCOSITY);
        density = it_node->FastGetSolutionStepValue(DENSITY);
        break;
    }

    for (int i_elem = 0; i_elem < static_cast<int>(mrModelPart.NumberOfElements()); ++i_elem){
        const auto it_elem = mrModelPart.ElementsBegin() + i_elem;
        it_elem->CalculateOnIntegrationPoints(VELOCITY, velocity_on_gauss_points, mrModelPart.GetProcessInfo());
        it_elem->CalculateOnIntegrationPoints(FLUID_FRACTION, fluid_fraction_on_gauss_points, mrModelPart.GetProcessInfo());
        it_elem->CalculateOnIntegrationPoints(VORTICITY, vorticity_on_gauss_points, mrModelPart.GetProcessInfo());
        const unsigned int ngauss = velocity_on_gauss_points.size();
        const double area = it_elem->GetGeometry().Area();
        const double Agauss = area/ngauss;
        for (unsigned int g = 0.0; g < ngauss; ++g){
            kinetic_energy += 0.5*Agauss*density*fluid_fraction_on_gauss_points[g]*(std::pow(velocity_on_gauss_points[g][0],2)+std::pow(velocity_on_gauss_points[g][1],2)+std::pow(velocity_on_gauss_points[g][2],2));
            dissipation_rate += Agauss*density*viscosity*(std::pow(vorticity_on_gauss_points[g][0],2)+std::pow(vorticity_on_gauss_points[g][1],2)+std::pow(vorticity_on_gauss_points[g][2],2));
            total_volume+=Agauss*fluid_fraction_on_gauss_points[g];
        }
    }

    data[0] = kinetic_energy/total_volume;
    data[1] = dissipation_rate/total_volume;
    data[2] = total_volume;

    return data;

}

void TaylorGreenVortexPorosityFieldProcess::ExecuteBeforeSolutionLoop()
{
    this->SetInitialConditions();
    this->SetValuesOnIntegrationPoints();
}

void TaylorGreenVortexPorosityFieldProcess::ExecuteInitializeSolutionStep()
{
    const int step = mrModelPart.GetProcessInfo()[STEP];
    const int buffer_size = mrModelPart.GetBufferSize();
    if (step == 1 && buffer_size > 1)
        this->SetInitialConditions();
}


void TaylorGreenVortexPorosityFieldProcess::ExecuteFinalizeSolutionStep() {}

/* Protected functions ****************************************************/

void TaylorGreenVortexPorosityFieldProcess::SetInitialConditions()
{
    const unsigned int k = mWaveNumber;
    const double alpha_max = mAlphaMax;
    const double deviation = mMeanDeviation;
    const double mean_alpha = alpha_max - deviation;
    const double c = mean_alpha;

    // Computation of the BodyForce and Porosity fields
    for (auto it_node = mrModelPart.NodesBegin(); it_node != mrModelPart.NodesEnd(); it_node++){

        const double x1 = it_node->X();
        const double x2 = it_node->Y();
        const double x3 = it_node->Z();

        double& r_alpha = it_node->FastGetSolutionStepValue(FLUID_FRACTION);
        double& r_dalphat = it_node->FastGetSolutionStepValue(FLUID_FRACTION_RATE);

        double& r_alpha1 = it_node->FastGetSolutionStepValue(FLUID_FRACTION_GRADIENT_X);
        double& r_alpha2 = it_node->FastGetSolutionStepValue(FLUID_FRACTION_GRADIENT_Y);
        double& r_alpha3 = it_node->FastGetSolutionStepValue(FLUID_FRACTION_GRADIENT_Z);

        double& r_pressure = it_node->FastGetSolutionStepValue(PRESSURE);

        double& u1 = it_node->FastGetSolutionStepValue(VELOCITY_X);
        double& u1_1 = it_node->FastGetSolutionStepValue(VELOCITY_X,1);
        double& u1_2 = it_node->FastGetSolutionStepValue(VELOCITY_X,2);
        double& u2 = it_node->FastGetSolutionStepValue(VELOCITY_Y);
        double& u2_1 = it_node->FastGetSolutionStepValue(VELOCITY_Y,1);
        double& u2_2 = it_node->FastGetSolutionStepValue(VELOCITY_Y,2);
        double& u3 = it_node->FastGetSolutionStepValue(VELOCITY_Z);
        double& u3_1 = it_node->FastGetSolutionStepValue(VELOCITY_Z,1);
        double& u3_2 = it_node->FastGetSolutionStepValue(VELOCITY_Z,2);

        r_pressure = 1.0/16.0*(std::cos(2*x1) + std::cos(2*x2))*(std::cos(2*x3) + 2);

        r_alpha = mean_alpha - deviation*std::cos(k*x1)*std::cos(k*x2)*std::sin(k*x3);

        r_alpha1 = k*deviation*std::sin(k*x1)*std::cos(k*x2)*std::sin(k*x3);
        r_alpha2 = k*deviation*std::cos(k*x1)*std::sin(k*x2)*std::sin(k*x3);
        r_alpha3 = -(k*deviation*std::cos(k*x1)*std::cos(k*x2)*std::cos(k*x3));

        u1 = c*std::cos(x1)*std::sin(x2)*std::sin(x3)/(mean_alpha - deviation*std::cos(k*x1)*std::cos(k*x2)*std::sin(k*x3));

        u1_1 = c*std::cos(x1)*std::sin(x2)*std::sin(x3)/(mean_alpha - deviation*std::cos(k*x1)*std::cos(k*x2)*std::sin(k*x3));

        u1_2 = c*std::cos(x1)*std::sin(x2)*std::sin(x3)/(mean_alpha - deviation*std::cos(k*x1)*std::cos(k*x2)*std::sin(k*x3));

        u2 = -c*std::sin(x1)*std::cos(x2)*std::sin(x3)/(mean_alpha - deviation*std::cos(k*x1)*std::cos(k*x2)*std::sin(k*x3));

        u2_1 = -c*std::sin(x1)*std::cos(x2)*std::sin(x3)/(mean_alpha - deviation*std::cos(k*x1)*std::cos(k*x2)*std::sin(k*x3));

        u2_2 = -c*std::sin(x1)*std::cos(x2)*std::sin(x3)/(mean_alpha - deviation*std::cos(k*x1)*std::cos(k*x2)*std::sin(k*x3));

        u3 = 0.0/(mean_alpha - deviation*std::cos(k*x1)*std::cos(k*x2)*std::sin(k*x3));

        u3_1 = 0.0/(mean_alpha - deviation*std::cos(k*x1)*std::cos(k*x2)*std::sin(k*x3));

        u3_2 = 0.0/(mean_alpha - deviation*std::cos(k*x1)*std::cos(k*x2)*std::sin(k*x3));
    }

}

/* Private functions ****************************************************/
void TaylorGreenVortexPorosityFieldProcess::SetValuesOnIntegrationPoints()
{
    const unsigned int Dim = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
    const unsigned int step = mrModelPart.GetProcessInfo()[STEP];
    const unsigned int k = mWaveNumber;
    const double alpha_max = mAlphaMax;
    const double deviation = mMeanDeviation;
    const double mean_alpha = alpha_max - deviation;

    const unsigned int n_elem = mrModelPart.NumberOfElements();

    // Computation of the BodyForce and Porosity fields
    for (int i_elem = 0; i_elem < static_cast<int>(mrModelPart.NumberOfElements()); ++i_elem){

        const auto it_elem = mrModelPart.ElementsBegin() + i_elem;
        const unsigned int id_elem = it_elem->Id();
        const GeometryType& r_geometry = it_elem->GetGeometry();

        const GeometryData::IntegrationMethod integration_method = it_elem->GetIntegrationMethod();
        const auto& r_number_integration_points = r_geometry.IntegrationPointsNumber(integration_method);

        Matrix NContainer = r_geometry.ShapeFunctionsValues(integration_method);
        const unsigned int NumNodes = r_geometry.PointsNumber();

        std::vector<double> fluid_fraction_on_gauss_points;
        std::vector<array_1d<double, 3>> fluid_fraction_gradient_on_gauss_points;

        if (fluid_fraction_on_gauss_points.size() != r_number_integration_points)
            fluid_fraction_on_gauss_points.resize(r_number_integration_points);

        if (fluid_fraction_gradient_on_gauss_points.size() != r_number_integration_points)
            fluid_fraction_gradient_on_gauss_points.resize(r_number_integration_points);


        for (unsigned int g = 0; g < r_number_integration_points; g++){

            Matrix gauss_point_coordinates = ZeroMatrix(r_number_integration_points,Dim);
            fluid_fraction_gradient_on_gauss_points[g] = ZeroVector(3);

            for (unsigned int i = 0; i < NumNodes; ++i){
                const array_1d<double, 3>& r_coordinates = r_geometry[i].Coordinates();
                for (unsigned int d = 0; d < Dim; ++d)
                    gauss_point_coordinates(g,d) += NContainer(g,i) * r_coordinates[d];
                }

            const double x1 = gauss_point_coordinates(g,0);
            const double x2 = gauss_point_coordinates(g,1);
            const double x3 = gauss_point_coordinates(g,2);

            const double alpha = mean_alpha - deviation*std::cos(k*x1)*std::cos(k*x2)*std::sin(k*x3);

            const double alpha1 = k*deviation*std::sin(k*x1)*std::cos(k*x2)*std::sin(k*x3);
            const double alpha2 = k*deviation*std::cos(k*x1)*std::sin(k*x2)*std::sin(k*x3);
            const double alpha3 = -(k*deviation*std::cos(k*x1)*std::cos(k*x2)*std::cos(k*x3));

            fluid_fraction_on_gauss_points[g] = alpha;
            fluid_fraction_gradient_on_gauss_points[g][0] = alpha1;
            fluid_fraction_gradient_on_gauss_points[g][1] = alpha2;
            fluid_fraction_gradient_on_gauss_points[g][2] = alpha3;

        }
        it_elem->SetValuesOnIntegrationPoints(FLUID_FRACTION, fluid_fraction_on_gauss_points, mrModelPart.GetProcessInfo());
        it_elem->SetValuesOnIntegrationPoints(FLUID_FRACTION_GRADIENT, fluid_fraction_gradient_on_gauss_points, mrModelPart.GetProcessInfo());
    }
}
};  // namespace Kratos.