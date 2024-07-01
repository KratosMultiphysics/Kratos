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
#include "flow_past_cylinder_porosity_field_process.h"
#include "swimming_dem_application_variables.h"

// Other applications includes
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

/* Public functions *******************************************************/
FlowPastCylinderPorosityFieldProcess::FlowPastCylinderPorosityFieldProcess(
    ModelPart& rModelPart)
    : Process(),
      mrModelPart(rModelPart)
{}

FlowPastCylinderPorosityFieldProcess::FlowPastCylinderPorosityFieldProcess(
    ModelPart& rModelPart,
    Parameters& rParameters)
    : Process(),
      mrModelPart(rModelPart)
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);
}

FlowPastCylinderPorosityFieldProcess::FlowPastCylinderPorosityFieldProcess(
    Model &rModel,
    Parameters &rParameters)
    : Process(),
      mrModelPart(rModel.GetModelPart(rParameters["model_part_name"].GetString()))
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);

}


void FlowPastCylinderPorosityFieldProcess::CheckDefaultsAndProcessSettings(Parameters &rParameters)
{
    const Parameters default_parameters = GetDefaultParameters();

    rParameters.ValidateAndAssignDefaults(default_parameters);

    mAlphaMin    = rParameters["benchmark_parameters"]["alpha_min"].GetDouble();
    mAlphaMax    = rParameters["benchmark_parameters"]["alpha_max"].GetDouble();
    mLength      = rParameters["benchmark_parameters"]["characteristic_length"].GetDouble();
    mViscosity = rParameters["benchmark_parameters"]["viscosity"].GetDouble();
    mPlateauRadius = rParameters["benchmark_parameters"]["plateau_radius"].GetDouble();
    mBumpRadius = rParameters["benchmark_parameters"]["bump_radius"].GetDouble();
    mX1Origin = rParameters["benchmark_parameters"]["x1_origin"].GetDouble();
    mX2Origin = rParameters["benchmark_parameters"]["x2_origin"].GetDouble();
    mAlternativeFormulation = rParameters["benchmark_parameters"]["use_alternative_formulation"].GetBool();
}


const Parameters FlowPastCylinderPorosityFieldProcess::GetDefaultParameters() const
{
    const Parameters default_parameters( R"(
    {
                "model_part_name"          : "please_specify_model_part_name",
                "variable_name"            : "BODY_FORCE",
                "benchmark_name"           : "custom_body_force.vortex",
                "benchmark_parameters"     : {
                                                "alpha_min"   : 0.5,
                                                "alpha_max"   : 1.0,
                                                "viscosity"   : 1.0,
                                                "characteristic_length" : 16.0,
                                                "x1_origin"   : 10.0,
                                                "x2_origin"   : 4.0,
                                                "bump_radius" : 0.5,
                                                "plateau_radius" : 0.1
                },
                "compute_nodal_error"      : true,
                "print_convergence_output" : false,
                "output_parameters"        : {}
    }  )" );

    return default_parameters;
}

void FlowPastCylinderPorosityFieldProcess::Execute(){}

void FlowPastCylinderPorosityFieldProcess::ExecuteInitialize(){}

void FlowPastCylinderPorosityFieldProcess::ExecuteBeforeSolutionLoop()
{
    this->SetPorosityField();
    this->SetValuesOnIntegrationPoints();
}

void FlowPastCylinderPorosityFieldProcess::ExecuteInitializeSolutionStep()
{
    this->SetPorosityField();
    this->SetValuesOnIntegrationPoints();
}


void FlowPastCylinderPorosityFieldProcess::ExecuteFinalizeSolutionStep() {}

/* Protected functions ****************************************************/

void FlowPastCylinderPorosityFieldProcess::SetPorosityField()
{
    double Dim = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
    const double step = mrModelPart.GetProcessInfo()[STEP];
    const double L = mLength;
    Matrix I = IdentityMatrix(Dim, Dim);
    const double r1 = mPlateauRadius;
    const double r2 = mBumpRadius;
    const double x10 = mX1Origin;
    const double x20 = mX2Origin;
    const double alpha_min = mAlphaMin;

    // Computation of the BodyForce and Porosity fields
    for (auto it_node = mrModelPart.NodesBegin(); it_node != mrModelPart.NodesEnd(); it_node++){

        const double x1 = it_node->X();
        const double x2 = it_node->Y();

        double& r_alpha = it_node->FastGetSolutionStepValue(FLUID_FRACTION);
        double& r_dalphat = it_node->FastGetSolutionStepValue(FLUID_FRACTION_RATE);

        double& r_alpha1 = it_node->FastGetSolutionStepValue(FLUID_FRACTION_GRADIENT_X);
        double& r_alpha2 = it_node->FastGetSolutionStepValue(FLUID_FRACTION_GRADIENT_Y);

        double& u1 = it_node->FastGetSolutionStepValue(VELOCITY_X);
        double& u1_1 = it_node->FastGetSolutionStepValue(VELOCITY_X,1);
        double& u1_2 = it_node->FastGetSolutionStepValue(VELOCITY_X,2);
        double& u2 = it_node->FastGetSolutionStepValue(VELOCITY_Y);

        Matrix& r_sigma = it_node->FastGetSolutionStepValue(PERMEABILITY);

        double velocity_norm = std::sqrt(std::pow(u1,2) + std::pow(u2,2));

        if (x1 >= r1 &&  x1 <= r2){
            r_alpha = (mAlphaMax-mAlphaMin)/(r2-r1) * x1 + (mAlphaMin*r2-mAlphaMax*r1)/(r2-r1);
            r_alpha1 = (mAlphaMax-mAlphaMin)/(r2-r1);
            r_alpha2 = 0.0;

        }
        else{
            if (x1 < r1) {
                r_alpha = mAlphaMin;
        } else if (x1 > r2){
                r_alpha = mAlphaMax;
            }
        }
        if (step == 0){
            if (std::sqrt(std::pow(x1-4,2) + std::pow(x2-4,2)) <= 0.51){
                u1 = 0.0;

                u1_1 = 0.0;

                u1_2 = 0.0;
            }
            else{
                u1 = 1.0/r_alpha;

                u1_1 = 1.0/r_alpha;

                u1_2 = 1.0/r_alpha;
            }

            u2 = 0.0;
        }

        if (r_alpha == 1.0){
                r_sigma = ZeroMatrix(Dim, Dim);
            }
        else{
            // double p_diameter = 0.1;
            // double kappa = std::pow(r_alpha,3)*std::pow(p_diameter,2)/(150*std::pow((1-r_alpha),2));
            // double F_e = 1.75/std::sqrt(150*std::pow(r_alpha,3));
            // if (mAlternativeFormulation)
            //     r_sigma = (r_alpha * ((mViscosity*r_alpha)/kappa + std::pow(r_alpha,2)*F_e/std::sqrt(kappa)*velocity_norm)) * I;
            // else
            //     r_sigma = ((mViscosity*r_alpha)/kappa + std::pow(r_alpha,2)*F_e/std::sqrt(kappa)*velocity_norm) * I;
            double p_diameter = 0.01;
            double drag_coef = 0.44;
            if (mAlternativeFormulation)
                r_sigma = (3.0/4.0*drag_coef/p_diameter*r_alpha*velocity_norm)*I;
            else
                r_sigma = (3.0/4.0*drag_coef/p_diameter*velocity_norm)*I;
        }

    }

}

/* Private functions ****************************************************/
void FlowPastCylinderPorosityFieldProcess::SetValuesOnIntegrationPoints()
{
    unsigned int Dim = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
    const double L = mLength;
    Matrix I = IdentityMatrix(Dim, Dim);
    const double r1 = mPlateauRadius;
    const double r2 = mBumpRadius;
    const double x10 = mX1Origin;
    const double x20 = mX2Origin;
    const double alpha_min = mAlphaMin;

    double alpha,alpha1,alpha2,reaction;

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
        std::vector<array_1d<double,3>> velocity_on_gauss_points;
        std::vector<array_1d<double, 3>> fluid_fraction_gradient_on_gauss_points;
        std::vector<Matrix> reaction_on_gauss_points;

        if (fluid_fraction_on_gauss_points.size() != r_number_integration_points)
            fluid_fraction_on_gauss_points.resize(r_number_integration_points);

        if (velocity_on_gauss_points.size() != r_number_integration_points)
            velocity_on_gauss_points.resize(r_number_integration_points);

        if (reaction_on_gauss_points.size() != r_number_integration_points)
            reaction_on_gauss_points.resize(r_number_integration_points);

        if (fluid_fraction_gradient_on_gauss_points.size() != r_number_integration_points)
            fluid_fraction_gradient_on_gauss_points.resize(r_number_integration_points);

        it_elem->CalculateOnIntegrationPoints(VELOCITY, velocity_on_gauss_points, mrModelPart.GetProcessInfo());

        for (unsigned int g = 0; g < r_number_integration_points; g++){

            Matrix gauss_point_coordinates = ZeroMatrix(r_number_integration_points,Dim);
            fluid_fraction_gradient_on_gauss_points[g] = ZeroVector(3);
            reaction_on_gauss_points[g] = ZeroMatrix(Dim, Dim);

            double velocity_modulus = 0.0;
            double velocity_norm;

            for (unsigned int d = 0; d < Dim; ++d){
                velocity_modulus += std::pow(velocity_on_gauss_points[g][d],2);
            }

            velocity_norm = std::sqrt(velocity_modulus);

            for (unsigned int i = 0; i < NumNodes; ++i){
                const array_1d<double, 3>& r_coordinates = r_geometry[i].Coordinates();
                for (unsigned int d = 0; d < Dim; ++d)
                    gauss_point_coordinates(g,d) += NContainer(g,i) * r_coordinates[d];
                }

            const double x1 = gauss_point_coordinates(g,0);
            const double x2 = gauss_point_coordinates(g,1);

            if (x1 >= r1 &&  x1 <= r2){
                alpha = (mAlphaMax-mAlphaMin)/(r2-r1) * x1 + (mAlphaMin*r2-mAlphaMax*r1)/(r2-r1);
                alpha1 = (mAlphaMax-mAlphaMin)/(r2-r1);
                alpha2 = 0.0;

            }
            else{
                if (x1 < r1) {
                    alpha = mAlphaMin;
                    alpha1 = 0.0;
                    alpha2 = 0.0;
            } else if (x1 > r2){
                    alpha = mAlphaMax;
                    alpha1 = 0.0;
                    alpha2 = 0.0;
                }
            }

            if (alpha == 1.0){
                reaction = 0.0;
            }
            else{
                double p_diameter = 0.01;
                double drag_coef = 0.44;
                // double kappa = std::pow(alpha,3)*std::pow(p_diameter,2)/(150*std::pow((1-alpha),2));
                // double F_e = 1.75/std::sqrt(150*std::pow(alpha,3));
                // if (mAlternativeFormulation)
                //     reaction = alpha * ((mViscosity*alpha)/kappa + std::pow(alpha,2)*F_e/std::sqrt(kappa)*velocity_norm);
                // else
                //     reaction = (mViscosity*alpha)/kappa + std::pow(alpha,2)*F_e/std::sqrt(kappa)*velocity_norm;
                if (mAlternativeFormulation)
                    reaction = 3.0/4.0*drag_coef/p_diameter*alpha*velocity_norm;
                else
                    reaction = 3.0/4.0*drag_coef/p_diameter*velocity_norm;
            }

            reaction_on_gauss_points[g] = reaction*I;
            fluid_fraction_on_gauss_points[g] = alpha;
            fluid_fraction_gradient_on_gauss_points[g][0] = alpha1;
            fluid_fraction_gradient_on_gauss_points[g][1] = alpha2;
        }

        it_elem->CalculateOnIntegrationPoints(PERMEABILITY, reaction_on_gauss_points, mrModelPart.GetProcessInfo());
        it_elem->CalculateOnIntegrationPoints(FLUID_FRACTION, fluid_fraction_on_gauss_points, mrModelPart.GetProcessInfo());
        it_elem->CalculateOnIntegrationPoints(FLUID_FRACTION_GRADIENT, fluid_fraction_gradient_on_gauss_points, mrModelPart.GetProcessInfo());
    }
}
};  // namespace Kratos.