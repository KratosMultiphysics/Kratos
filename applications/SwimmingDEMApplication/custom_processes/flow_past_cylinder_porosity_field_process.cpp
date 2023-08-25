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
    mSigma = rParameters["benchmark_parameters"]["sigma"].GetDouble();
    mPlateauRadius = rParameters["benchmark_parameters"]["plateau_radius"].GetDouble();
    mBumpRadius = rParameters["benchmark_parameters"]["bump_radius"].GetDouble();
    mX1Origin = rParameters["benchmark_parameters"]["x1_origin"].GetDouble();
    mX2Origin = rParameters["benchmark_parameters"]["x2_origin"].GetDouble();
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
                                                "sigma"       : 1.0,
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
}

void FlowPastCylinderPorosityFieldProcess::ExecuteInitializeSolutionStep(){}


void FlowPastCylinderPorosityFieldProcess::ExecuteFinalizeSolutionStep() {}

/* Protected functions ****************************************************/

void FlowPastCylinderPorosityFieldProcess::SetPorosityField()
{
    double Dim = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
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

        Matrix& sigma = it_node->FastGetSolutionStepValue(PERMEABILITY);

        sigma = mSigma * I;

        if (x1 >= r1 &&  x1 <= r2){
            r_alpha = (mAlphaMax-mAlphaMin)/(r2-r1) * x1 + (mAlphaMin*r2-mAlphaMax*r1)/(r2-r1);

            r_alpha1 = (mAlphaMax-mAlphaMin)/(r2-r1);

            r_alpha2 = 0.0;
        }
        else{
            if (x1 < r1) {
                r_alpha = mAlphaMin;
                r_alpha1 = 0.0;
                r_alpha2 = 0.0;
        } else if (x1 > r2){
                r_alpha = mAlphaMax;
                r_alpha1 = 0.0;
                r_alpha2 = 0.0;
            }
        }

        // r_alpha = (mAlphaMin-mAlphaMax)/L*x1 + mAlphaMax;

        // r_dalphat = 0.0;

        // r_alpha1 = (mAlphaMin-mAlphaMax)/L;

        // r_alpha2 = 0.0;

        // double r_sq = std::pow(x1-x10,2) + std::pow(x2-x20,2);

        // if ((std::sqrt(r_sq) > r1) && (std::sqrt(r_sq) < r2)){

        //     r_alpha = -(1 - alpha_min)*(1 - std::exp(-1/((-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) + std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2))))/(std::exp(-1/(1 - (-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) - std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2)))) + std::exp(-1/((-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) + std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2)))))) + 1;

        //     r_alpha1 = -(1 - alpha_min)*(-((2*x1 - 2*x10)*std::exp(-1/(1 - (-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) - std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2))))/((-std::pow(r1,2) + std::pow(r2,2))*std::pow((1 - (-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) - std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2))),2)) - (2*x1 - 2*x10)*std::exp(-1/((-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) + std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2))))/((-std::pow(r1,2) + std::pow(r2,2))*std::pow(((-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) + std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2))),2)))*std::exp(-1/((-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) + std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2))))/std::pow((std::exp(-1/(1 - (-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) - std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2)))) + std::exp(-1/((-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) + std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2))))),2) - (2*x1 - 2*x10)*std::exp(-1/((-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) + std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2))))/((-std::pow(r1,2) + std::pow(r2,2))*std::pow(((-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) + std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2))),2)*(std::exp(-1/(1 - (-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) - std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2)))) + std::exp(-1/((-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) + std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2)))))));

        //     r_alpha2 = -(1 - alpha_min)*(-((2*x2 - 2*x20)*std::exp(-1/(1 - (-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) - std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2))))/((-std::pow(r1,2) + std::pow(r2,2))*std::pow((1 - (-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) - std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2))),2)) - (2*x2 - 2*x20)*std::exp(-1/((-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) + std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2))))/((-std::pow(r1,2) + std::pow(r2,2))*std::pow(((-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) + std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2))),2)))*std::exp(-1/((-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) + std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2))))/std::pow((std::exp(-1/(1 - (-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) - std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2)))) + std::exp(-1/((-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) + std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2))))),2) - (2*x2 - 2*x20)*std::exp(-1/((-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) + std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2))))/((-std::pow(r1,2) + std::pow(r2,2))*std::pow(((-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) + std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2))),2)*(std::exp(-1/(1 - (-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) - std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2)))) + std::exp(-1/((-std::pow(r1,2) + std::pow((x1 - x10),2))/(-std::pow(r1,2) + std::pow(r2,2)) + std::pow((x2 - x20),2)/(-std::pow(r1,2) + std::pow(r2,2)))))));

        // }
        // else{
        //     if (std::sqrt(r_sq) <= r1) {
        //         r_alpha = mAlphaMin;
        //         r_alpha1 = 0.0;
        //         r_alpha2 = 0.0;
        //     } else if (std::sqrt(r_sq) >= r2){
        //         r_alpha = mAlphaMax;
        //         r_alpha1 = 0.0;
        //         r_alpha2 = 0.0;
        //     }
        // }

        if (std::sqrt(std::pow(x1-4,2) + std::pow(x2-4,2)) <= 0.51){
            u1 = 0.0;

            u1_1 = 0.0;

            u1_2 = 0.0;
        }
        else{
            u1 = 1.0;

            u1_1 = 1.0;

            u1_2 = 1.0;
        }

        u2 = 0.0;

    }

}

/* Private functions ****************************************************/

};  // namespace Kratos.