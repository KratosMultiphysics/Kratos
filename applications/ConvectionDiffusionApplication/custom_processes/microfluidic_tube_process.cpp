//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Aniol Sala Pascual
//
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/checks.h"
#include "utilities/variable_utils.h"
#include "utilities/math_utils.h"

// Application includes
#include "convection_diffusion_application.h"
#include "microfluidic_tube_process.h"


namespace Kratos{

/* Public functions *******************************************************/
MicrofluidicTubeProcess::MicrofluidicTubeProcess(
    ModelPart& rModelPart)
    : Process(),
    mrModelPart(rModelPart)
{}

MicrofluidicTubeProcess::MicrofluidicTubeProcess(
    ModelPart& rModelPart,
    Parameters& rParameters)
    : Process(),
      mrModelPart(rModelPart)
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);
}

MicrofluidicTubeProcess::MicrofluidicTubeProcess(
    Model &rModel,
    Parameters &rParameters)
    : Process(),
      mrModelPart(rModel.GetModelPart(rParameters["model_part_name"].GetString()))
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);

}

double MicrofluidicTubeProcess::GetVelocityValue(Node<3> node)
{
    return 0.0;
}

void MicrofluidicTubeProcess::CheckDefaultsAndProcessSettings(Parameters &rParameters)
{

    const Parameters default_parameters = GetDefaultParameters();
    rParameters.ValidateAndAssignDefaults(default_parameters);

    // Read settings
    mpOutletModelPart = &(mrModelPart.GetParentModelPart().GetSubModelPart(rParameters["outlet_model_part_name"].GetString()));
    mGeometryType = rParameters["geometry_type"].GetString();
    mDimensions = rParameters["dimensions"].GetVector();
    mPressureOutlet = rParameters["pressure_outlet"].GetDouble();
    mDiffusionCoefficient = rParameters["diffusion_coefficient"].GetDouble();

    // If the variable specific heat is defined:
    if (mrModelPart.NodesBegin()->SolutionStepsDataHas(SPECIFIC_HEAT) == true) {
        // Load the material properties from json file
        const std::string& r_input_filename = rParameters["materials_filename"].GetString();
        std::ifstream infile(r_input_filename);
        KRATOS_ERROR_IF_NOT(infile.good()) << "Materials file: " << r_input_filename << " cannot be found" << std::endl;
        std::stringstream buffer;
        buffer << infile.rdbuf();
        Parameters buoyancyParameters(buffer.str());
        Parameters material_parameters = buoyancyParameters["properties"][0]["Material"];

        // Read the values and change the value of the specific heat
        mSpecificHeat = material_parameters["Variables"]["SPECIFIC_HEAT"].GetDouble();
        mConductivity = material_parameters["Variables"]["CONDUCTIVITY"].GetDouble();
        mDensity = material_parameters["Variables"]["DENSITY"].GetDouble();
        mViscosity = material_parameters["Variables"]["DYNAMIC_VISCOSITY"].GetDouble();

        mSpecificHeat = mConductivity / (mDensity * mDiffusionCoefficient);
    }
}

const Parameters MicrofluidicTubeProcess::GetDefaultParameters() const
{
    const Parameters default_parameters( R"(
    {
                "model_part_name"              : "please_specify_model_part_name",
                "outlet_model_part_name"       : "outlet_model_part_name",
                "geometry_type"                : "circular",
                "pressure_outlet"              : 0.0,
                "dimensions"                   : [1.0, 1.0],
                "diffusion_coefficient"        : 1.0,
                "materials_filename"           : "json_file"
    }  )" );

    return default_parameters;
}


void MicrofluidicTubeProcess::Execute()
{
    this->ExecuteInitialize();
    this->ExecuteInitializeSolutionStep();
}

void MicrofluidicTubeProcess::ExecuteInitialize()
{
}

void MicrofluidicTubeProcess::ExecuteBeforeSolutionLoop()
{
    if (mrModelPart.NodesBegin()->SolutionStepsDataHas(SPECIFIC_HEAT) == true) {
        this->SetSystemProperties();
    }
}

void MicrofluidicTubeProcess::ExecuteInitializeSolutionStep()
{
    this->SetOutletVelocityAndPressure();
}

void MicrofluidicTubeProcess::ExecuteFinalizeSolutionStep() {
}

/* Protected functions ****************************************************/

void MicrofluidicTubeProcess::SetOutletVelocityAndPressure() {
    // To compute velocity at the outlet (parabolic profile)
    Vector r0 = this->GetSurfaceCenter();
    double r2;
    double rad = .5e-3;
    double a2 = std::pow(rad, 2);

    bool pressureImposed = false;
    for (auto it_node = mpOutletModelPart->NodesBegin(); it_node != mpOutletModelPart->NodesEnd(); it_node++){
        if (!pressureImposed)
        {
            it_node->pGetDof(PRESSURE)->FixDof();
            it_node->FastGetSolutionStepValue(PRESSURE) = mPressureOutlet;
            pressureImposed = true;
        }
        
        // For the moment assume the inlet normal is in the x direction
        it_node->pGetDof(VELOCITY_Y)->FixDof();
        it_node->pGetDof(VELOCITY_Z)->FixDof();
        it_node->FastGetSolutionStepValue(VELOCITY_Y) = 0.0;
        it_node->FastGetSolutionStepValue(VELOCITY_Z) = 0.0;


        // Assign velocity
        r2 = std::pow(r0[0] - it_node->Y(), 2) + std::pow(r0[1] - it_node->Z(), 2);
        
        it_node->pGetDof(VELOCITY_X)->FixDof();
        it_node->FastGetSolutionStepValue(VELOCITY_X) = -3.0e-4 * (1. - r2 / a2);
    }
    // for (auto it_node = mrModelPart.NodesBegin(); it_node != mrModelPart.NodesEnd(); it_node++){
    //     // Assign boundary conditions at the inlet
    //     it_node->pGetDof(TEMPERATURE)->FixDof();
    //     if (it_node->Z() > 1.5e-3)
    //     {
    //         it_node->FastGetSolutionStepValue(TEMPERATURE) = 1070.;
    //     } else {
    //         it_node->FastGetSolutionStepValue(TEMPERATURE) = 1085.;
    //     }
    // }

}

void MicrofluidicTubeProcess::SetSystemProperties()
{
    // Change the specific heat to match the desired coefficient diffusion
    (mrModelPart.pGetProperties(1))->SetValue(SPECIFIC_HEAT, mSpecificHeat);

    block_for_each(mrModelPart.Elements(), [&](Element& rElement){
        rElement.SetProperties(mrModelPart.pGetProperties(1));
    });

    block_for_each(mrModelPart.Nodes(), [&](Node<3>& rNode){
        rNode.FastGetSolutionStepValue(SPECIFIC_HEAT) = mSpecificHeat;
    });
    // std::cout << "Setting C = " << mSpecificHeat << std::endl;
    // std::cout << "##########################\n\n\n\n" << std::endl;
}
/* Private functions ****************************************************/

Vector MicrofluidicTubeProcess::GetSurfaceCenter() {
    Vector r0 = ZeroVector(2);
    
    for (auto it_node = mpOutletModelPart->NodesBegin(); it_node != mpOutletModelPart->NodesEnd(); it_node++){
        double y_node = it_node->Y();
        double z_node = it_node->Z();

        r0[0] += y_node;
        r0[1] += z_node;

    }
    const int n_nodes = mpOutletModelPart->Nodes().size();
    r0[0] /= n_nodes;
    r0[1] /= n_nodes;

    return r0;
}

};  // namespace Kratos.