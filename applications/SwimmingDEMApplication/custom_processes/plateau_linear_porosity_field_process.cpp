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
#include "plateau_linear_porosity_field_process.h"

#include "swimming_dem_application_variables.h"

// Other applications includes
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

/* Public functions *******************************************************/
PlateauLinearPorosityFieldProcess::PlateauLinearPorosityFieldProcess(
    ModelPart& rModelPart)
    : Process(),
      mrModelPart(rModelPart)
{}

PlateauLinearPorosityFieldProcess::PlateauLinearPorosityFieldProcess(
    ModelPart& rModelPart,
    Parameters& rParameters)
    : Process(),
      mrModelPart(rModelPart)
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);
}

PlateauLinearPorosityFieldProcess::PlateauLinearPorosityFieldProcess(
    Model &rModel,
    Parameters &rParameters)
    : Process(),
      mrModelPart(rModel.GetModelPart(rParameters["model_part_name"].GetString()))
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);

}


void PlateauLinearPorosityFieldProcess::CheckDefaultsAndProcessSettings(Parameters &rParameters)
{
    const Parameters default_parameters = GetDefaultParameters();

    rParameters.ValidateAndAssignDefaults(default_parameters);

    mAlphaMin    = rParameters["benchmark_parameters"]["alpha_min"].GetDouble();
    mAlphaMax    = rParameters["benchmark_parameters"]["alpha_max"].GetDouble();
    mBeginningSlope = rParameters["benchmark_parameters"]["begin_slope"].GetDouble();
    mBeginningPlateau = rParameters["benchmark_parameters"]["begin_plateau"].GetDouble();
    mEndingPlateau = rParameters["benchmark_parameters"]["ending_plateau"].GetDouble();
    mEndingSlope = rParameters["benchmark_parameters"]["ending_slope"].GetDouble();
    mInletVelocity = rParameters["benchmark_parameters"]["inlet_velocity"].GetDouble();
}


const Parameters PlateauLinearPorosityFieldProcess::GetDefaultParameters() const
{
    const Parameters default_parameters( R"(
    {
                "model_part_name"          : "please_specify_model_part_name",
                "variable_name"            : "BODY_FORCE",
                "benchmark_name"           : "custom_body_force.vortex",
                "benchmark_parameters"     : {
                                                "alpha_min"      : 0.5,
                                                "alpha_max"      : 1.0,
                                                "begin_slope"    : 1.0,
                                                "begin_plateau"  : 16.0,
                                                "ending_plateau" : 10.0,
                                                "ending_slope"   : 0.5,
                                                "inlet_velocity" : 0.1
                },
                "compute_nodal_error"      : true,
                "print_convergence_output" : false,
                "output_parameters"        : {}
    }  )" );

    return default_parameters;
}

void PlateauLinearPorosityFieldProcess::Execute(){}

void PlateauLinearPorosityFieldProcess::ExecuteInitialize(){}

void PlateauLinearPorosityFieldProcess::ExecuteBeforeSolutionLoop()
{
    this->SetPorosityField();
    this->SetValuesOnIntegrationPoints();
}

void PlateauLinearPorosityFieldProcess::ExecuteInitializeSolutionStep()
{
    this->SetPorosityField();
    this->SetValuesOnIntegrationPoints();
}


void PlateauLinearPorosityFieldProcess::ExecuteFinalizeSolutionStep() {}

/* Protected functions ****************************************************/

void PlateauLinearPorosityFieldProcess::SetPorosityField()
{
    const int Dim = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
    const double time = mrModelPart.GetProcessInfo()[TIME];
    const double step = mrModelPart.GetProcessInfo()[STEP];
    const double r1 = mBeginningSlope + mInletVelocity*time;
    const double r2 = mBeginningPlateau + mInletVelocity*time;
    const double r3 = mEndingPlateau + mInletVelocity*time;
    const double r4 = mEndingSlope + mInletVelocity*time;

    // Computation of the BodyForce and Porosity fields
    for (auto it_node = mrModelPart.NodesBegin(); it_node != mrModelPart.NodesEnd(); it_node++){

        const double x1 = it_node->X();
        const double x2 = it_node->Y();

        double& r_alpha = it_node->FastGetSolutionStepValue(FLUID_FRACTION);

        double& r_alpha1 = it_node->FastGetSolutionStepValue(FLUID_FRACTION_GRADIENT_X);
        double& r_alpha2 = it_node->FastGetSolutionStepValue(FLUID_FRACTION_GRADIENT_Y);

        double& r_dalphat = it_node->FastGetSolutionStepValue(FLUID_FRACTION_RATE);
        Vector& exact_u = it_node->FastGetSolutionStepValue(EXACT_VELOCITY);
        exact_u = ZeroVector(3);
        double& exact_u1 = it_node->FastGetSolutionStepValue(EXACT_VELOCITY)[0];
        double& exact_u2 = it_node->FastGetSolutionStepValue(EXACT_VELOCITY)[1];

        double& u1_1 = it_node->FastGetSolutionStepValue(VELOCITY_X,1);
        double& u1_2 = it_node->FastGetSolutionStepValue(VELOCITY_X,2);
        double& u1 = it_node->FastGetSolutionStepValue(VELOCITY_X);

        if (x1 > r1 && x1 < r2){
            double eta = (x1-r2)/(-r2+r1);
            r_alpha = 1.0-(1.0-mAlphaMin)/(1+std::exp((2*eta-1)/(eta*(1-eta))));
            r_alpha1 = (1.0 - mAlphaMin)*(-(-1 + 2*(-r2 + x1)/(r1 - r2))*(r1 - r2)/((1 - (-r2 + x1)/(r1 - r2))*std::pow(-r2 + x1,2)) + (-1 + 2*(-r2 + x1)/(r1 - r2))/(std::pow(1 - (-r2 + x1)/(r1 - r2),2)*(-r2 + x1)) + 2/((1 - (-r2 + x1)/(r1 - r2))*(-r2 + x1)))*std::exp((-1 + 2*(-r2 + x1)/(r1 - r2))*(r1 - r2)/((1 - (-r2 + x1)/(r1 - r2))*(-r2 + x1)))/std::pow(std::exp((-1 + 2*(-r2 + x1)/(r1 - r2))*(r1 - r2)/((1 - (-r2 + x1)/(r1 - r2))*(-r2 + x1))) + 1,2);
            if (r_alpha == 1.0)
                r_alpha1 = 0.0;
            r_alpha2 = 0.0;
        }
        else if (x1 >= r2 && x1 <= r3){
            r_alpha = mAlphaMin;
            r_alpha1 = 0.0;
            r_alpha2 = 0.0;
        }
        else if (x1 > r3 && x1 < r4){
            double eta = (x1-r3)/(-r3+r4);
            r_alpha = 1.0-(1.0-mAlphaMin)/(1+std::exp((2*eta-1)/(eta*(1-eta))));
            r_alpha1 = (1.0 - mAlphaMin)*(-(-1 + 2*(-r3 + x1)/(-r3 + r4))*(-r3 + r4)/((1 - (-r3 + x1)/(-r3 + r4))*std::pow(-r3 + x1,2)) + (-1 + 2*(-r3 + x1)/(-r3 + r4))/(std::pow(1 - (-r3 + x1)/(-r3 + r4),2)*(-r3 + x1)) + 2/((1 - (-r3 + x1)/(-r3 + r4))*(-r3 + x1)))*std::exp((-1 + 2*(-r3 + x1)/(-r3 + r4))*(-r3 + r4)/((1 - (-r3 + x1)/(-r3 + r4))*(-r3 + x1)))/std::pow(std::exp((-1 + 2*(-r3 + x1)/(-r3 + r4))*(-r3 + r4)/((1 - (-r3 + x1)/(-r3 + r4))*(-r3 + x1))) + 1,2);
            if (r_alpha == 1.0)
                r_alpha1 = 0.0;
            r_alpha2 = 0.0;
        }
        else if (x1 <= r1  || x1 >= r4) {
            r_alpha = mAlphaMax;
            r_alpha1 = 0.0;
            r_alpha2 = 0.0;
        }
        if (step <= 2){
            if (x1 > r1 && x1 < r2){
                r_dalphat = (1.0 - mAlphaMin)*(mInletVelocity*(-1 + 2*(-mBeginningPlateau - time*mInletVelocity + x1)/(mBeginningSlope - mBeginningPlateau))*(mBeginningSlope - mBeginningPlateau)/((1 - (-mBeginningPlateau - time*mInletVelocity + x1)/(mBeginningSlope - mBeginningPlateau))*std::pow(-mBeginningPlateau - time*mInletVelocity + x1,2)) - mInletVelocity*(-1 + 2*(-mBeginningPlateau - time*mInletVelocity + x1)/(mBeginningSlope - mBeginningPlateau))/(std::pow(1 - (-mBeginningPlateau - time*mInletVelocity + x1)/(mBeginningSlope - mBeginningPlateau),2)*(-mBeginningPlateau - time*mInletVelocity + x1)) - 2*mInletVelocity/((1 - (-mBeginningPlateau - time*mInletVelocity + x1)/(mBeginningSlope - mBeginningPlateau))*(-mBeginningPlateau - time*mInletVelocity + x1)))*std::exp((-1 + 2*(-mBeginningPlateau - time*mInletVelocity + x1)/(mBeginningSlope - mBeginningPlateau))*(mBeginningSlope - mBeginningPlateau)/((1 - (-mBeginningPlateau - time*mInletVelocity + x1)/(mBeginningSlope - mBeginningPlateau))*(-mBeginningPlateau - time*mInletVelocity + x1)))/std::pow(std::exp((-1 + 2*(-mBeginningPlateau - time*mInletVelocity + x1)/(mBeginningSlope - mBeginningPlateau))*(mBeginningSlope - mBeginningPlateau)/((1 - (-mBeginningPlateau - time*mInletVelocity + x1)/(mBeginningSlope - mBeginningPlateau))*(-mBeginningPlateau - time*mInletVelocity + x1))) + 1,2);
                if (r_alpha == 1.0 && r_alpha1 == 0.0)
                    r_dalphat = 0.0;
            }
            else if (x1 >= r2 && x1 <= r3){
                r_dalphat = 0.0;
            }
            else if (x1 > r3 && x1 < r4){
                r_dalphat = (1.0 - mAlphaMin)*(mInletVelocity*(-1 + 2*(-mEndingPlateau - time*mInletVelocity + x1)/(-mEndingPlateau + mEndingSlope))*(-mEndingPlateau + mEndingSlope)/((1 - (-mEndingPlateau - time*mInletVelocity + x1)/(-mEndingPlateau + mEndingSlope))*std::pow(-mEndingPlateau - time*mInletVelocity + x1,2)) - mInletVelocity*(-1 + 2*(-mEndingPlateau - time*mInletVelocity + x1)/(-mEndingPlateau + mEndingSlope))/(std::pow(1 - (-mEndingPlateau - time*mInletVelocity + x1)/(-mEndingPlateau + mEndingSlope),2)*(-mEndingPlateau - time*mInletVelocity + x1)) - 2*mInletVelocity/((1 - (-mEndingPlateau - time*mInletVelocity + x1)/(-mEndingPlateau + mEndingSlope))*(-mEndingPlateau - time*mInletVelocity + x1)))*std::exp((-1 + 2*(-mEndingPlateau - time*mInletVelocity + x1)/(-mEndingPlateau + mEndingSlope))*(-mEndingPlateau + mEndingSlope)/((1 - (-mEndingPlateau - time*mInletVelocity + x1)/(-mEndingPlateau + mEndingSlope))*(-mEndingPlateau - time*mInletVelocity + x1)))/std::pow(std::exp((-1 + 2*(-mEndingPlateau - time*mInletVelocity + x1)/(-mEndingPlateau + mEndingSlope))*(-mEndingPlateau + mEndingSlope)/((1 - (-mEndingPlateau - time*mInletVelocity + x1)/(-mEndingPlateau + mEndingSlope))*(-mEndingPlateau - time*mInletVelocity + x1))) + 1,2);
                if (r_alpha == 1.0 && r_alpha1 == 0.0)
                    r_dalphat = 0.0;
            }
            else if (x1 <= r1  || x1 >= r4) {
                r_dalphat = 0.0;
            }
        }


        exact_u1 = mInletVelocity;
        exact_u2 = 0.0;
        if (step == 0){
            u1 = mInletVelocity;
            // u1_1 = mInletVelocity;
            // u1_2 = mInletVelocity;
        }
    }

}

void PlateauLinearPorosityFieldProcess::SetValuesOnIntegrationPoints()
{
    const double Dim = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
    const double time = mrModelPart.GetProcessInfo()[TIME];
    const double step = mrModelPart.GetProcessInfo()[STEP];
    const double r1 = mBeginningSlope + mInletVelocity*time;
    const double r2 = mBeginningPlateau + mInletVelocity*time;
    const double r3 = mEndingPlateau + mInletVelocity*time;
    const double r4 = mEndingSlope + mInletVelocity*time;

    double u1,u2,alpha,alpha1,alpha2;

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

        std::vector<double> fluid_fraction_on_gauss_points;
        std::vector<Vector> velocity_on_gauss_points;
        std::vector<array_1d<double, 3>> fluid_fraction_gradient_on_gauss_points;

        if (fluid_fraction_on_gauss_points.size() != r_number_integration_points)
            fluid_fraction_on_gauss_points.resize(r_number_integration_points);

        if (velocity_on_gauss_points.size() != r_number_integration_points)
            velocity_on_gauss_points.resize(r_number_integration_points);
        if (fluid_fraction_gradient_on_gauss_points.size() != r_number_integration_points)
            fluid_fraction_gradient_on_gauss_points.resize(r_number_integration_points);

        for (unsigned int g = 0; g < r_number_integration_points; g++){

            Matrix gauss_point_coordinates = ZeroMatrix(r_number_integration_points,Dim);
            fluid_fraction_gradient_on_gauss_points[g] = ZeroVector(3);
            velocity_on_gauss_points[g] = ZeroVector(3);

            for (unsigned int i = 0; i < NumNodes; ++i){
                const array_1d<double, 3>& r_coordinates = r_geometry[i].Coordinates();
                for (unsigned int d = 0; d < Dim; ++d)
                    gauss_point_coordinates(g,d) += NContainer(g,i) * r_coordinates[d];
                }

            const double x1 = gauss_point_coordinates(g,0);
            const double x2 = gauss_point_coordinates(g,1);

            if (x1 > r1 && x1 < r2 ){
                double eta = (x1-r2)/(-r2+r1);
                alpha = 1.0-(1.0-mAlphaMin)/(1+std::exp((2*eta-1)/(eta*(1-eta))));
                alpha1 = (1.0 - mAlphaMin)*(-(-1 + 2*(-r2 + x1)/(r1 - r2))*(r1 - r2)/((1 - (-r2 + x1)/(r1 - r2))*std::pow(-r2 + x1,2)) + (-1 + 2*(-r2 + x1)/(r1 - r2))/(std::pow(1 - (-r2 + x1)/(r1 - r2),2)*(-r2 + x1)) + 2/((1 - (-r2 + x1)/(r1 - r2))*(-r2 + x1)))*std::exp((-1 + 2*(-r2 + x1)/(r1 - r2))*(r1 - r2)/((1 - (-r2 + x1)/(r1 - r2))*(-r2 + x1)))/std::pow(std::exp((-1 + 2*(-r2 + x1)/(r1 - r2))*(r1 - r2)/((1 - (-r2 + x1)/(r1 - r2))*(-r2 + x1))) + 1,2);
                if (alpha == 1.0)
                    alpha1 = 0.0;
                alpha2 = 0.0;
            }
            else if (x1 >= r2 && x1 <= r3){
                alpha = mAlphaMin;
                alpha1 = 0.0;
                alpha2 = 0.0;
            }
            else if (x1 > r3  && x1 < r4){
                double eta = (x1-r3)/(-r3+r4);
                alpha = 1.0-(1.0-mAlphaMin)/(1+std::exp((2*eta-1)/(eta*(1-eta))));
                alpha1 = (1.0 - mAlphaMin)*(-(-1 + 2*(-r3 + x1)/(-r3 + r4))*(-r3 + r4)/((1 - (-r3 + x1)/(-r3 + r4))*std::pow(-r3 + x1,2)) + (-1 + 2*(-r3 + x1)/(-r3 + r4))/(std::pow(1 - (-r3 + x1)/(-r3 + r4),2)*(-r3 + x1)) + 2/((1 - (-r3 + x1)/(-r3 + r4))*(-r3 + x1)))*std::exp((-1 + 2*(-r3 + x1)/(-r3 + r4))*(-r3 + r4)/((1 - (-r3 + x1)/(-r3 + r4))*(-r3 + x1)))/std::pow(std::exp((-1 + 2*(-r3 + x1)/(-r3 + r4))*(-r3 + r4)/((1 - (-r3 + x1)/(-r3 + r4))*(-r3 + x1))) + 1,2);
                if (alpha == 1.0)
                    alpha1 = 0.0;
                alpha2 = 0.0;
                }
            else if (x1 <= r1  || x1 >= r4) {
                alpha = mAlphaMax;
                alpha1 = 0.0;
                alpha2 = 0.0;
            }

            u1 = mInletVelocity;
            u2 = 0.0;

            fluid_fraction_on_gauss_points[g] = alpha;
            fluid_fraction_gradient_on_gauss_points[g][0] = alpha1;
            fluid_fraction_gradient_on_gauss_points[g][1] = alpha2;
            velocity_on_gauss_points[g][0] = u1;
            velocity_on_gauss_points[g][1] = u2;
        }

        it_elem->CalculateOnIntegrationPoints(FLUID_FRACTION, fluid_fraction_on_gauss_points, mrModelPart.GetProcessInfo());
        it_elem->CalculateOnIntegrationPoints(FLUID_FRACTION_GRADIENT, fluid_fraction_gradient_on_gauss_points, mrModelPart.GetProcessInfo());
        it_elem->SetValuesOnIntegrationPoints(EXACT_VELOCITY, velocity_on_gauss_points, mrModelPart.GetProcessInfo());
    }
}

/* Private functions ****************************************************/

};  // namespace Kratos.