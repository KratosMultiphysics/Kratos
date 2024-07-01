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
#include "manufactured_body_force_process.h"


namespace Kratos{

/* Public functions *******************************************************/
ManufacturedBodyForceProcess::ManufacturedBodyForceProcess(
    ModelPart& rModelPart)
    : Process(),
    mrModelPart(rModelPart)
{}

ManufacturedBodyForceProcess::ManufacturedBodyForceProcess(
    ModelPart& rModelPart,
    Parameters& rParameters)
    : Process(),
      mrModelPart(rModelPart)
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);
}

ManufacturedBodyForceProcess::ManufacturedBodyForceProcess(
    Model &rModel,
    Parameters &rParameters)
    : Process(),
      mrModelPart(rModel.GetModelPart(rParameters["model_part_name"].GetString()))
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);

}

double ManufacturedBodyForceProcess::GetFakeOutletPressure()
{
    // Assign inlet pressure based on the gradient and the
    // outlet position

    // Get the position of the false outlet
    double mean_x = 0.0;
    int n_nodes = 0;
    for (auto it_node = mpOutletModelPart->NodesBegin(); it_node != mpOutletModelPart->NodesEnd(); it_node++){
        mean_x += it_node->X();
    }
    mean_x /= n_nodes;

    // Impose a pressure gradient \Delta p = 1e-3, with p_in = 1e-3, p_out = 0.0
    double p_inlet = 1e-3, p_outlet = 0.0;
    double p_fake_outlet_pressure = p_inlet + (p_outlet - p_inlet) * mean_x / mLen;

    return p_fake_outlet_pressure;
}

void ManufacturedBodyForceProcess::CheckDefaultsAndProcessSettings(Parameters &rParameters)
{

    const Parameters default_parameters = GetDefaultParameters();
    rParameters.ValidateAndAssignDefaults(default_parameters);

    // Read the properties from json file
    const std::string& r_input_filename = rParameters["materials_filename"].GetString();
    std::ifstream infile(r_input_filename);
    KRATOS_ERROR_IF_NOT(infile.good()) << "Materials file: " << r_input_filename << " cannot be found" << std::endl;
    std::stringstream buffer;
    buffer << infile.rdbuf();
    Parameters buoyancyParameters(buffer.str());

    // Read settings
    mv0 = rParameters["central_velocity"].GetDouble();
    mRad = rParameters["radius"].GetDouble();
    mLen = rParameters["length"].GetDouble();
    mImposeConditionsInlet = rParameters["impose_inlet_conditions"].GetBool();
    mWriteExactSolutionOutput = rParameters["write_exact_solution_output"].GetBool();
    mCavity = rParameters["cavity"].GetBool();
    mpInletModelPart = &(mrModelPart.GetParentModelPart().GetSubModelPart(rParameters["inlet_model_part_name"].GetString()));
    mpOutletModelPart = &(mrModelPart.GetParentModelPart().GetSubModelPart(rParameters["outlet_model_part_name"].GetString()));
    mAddHeatFluxSource = rParameters["add_heat_flux_source"].GetBool();
    mAddBodyForceSource = rParameters["add_body_force_source"].GetBool();
    mInletTemperature = rParameters["temperature_inlet"].GetDouble();
    mTGrad = rParameters["temperature_gradient"].GetDouble();

    Parameters material_parameters = buoyancyParameters["properties"][0]["Material"];
    mSpecificHeat = material_parameters["Variables"]["SPECIFIC_HEAT"].GetDouble();
    mConductivity = material_parameters["Variables"]["CONDUCTIVITY"].GetDouble();
    mDensity = material_parameters["Variables"]["DENSITY"].GetDouble();
    mViscosity = material_parameters["Variables"]["DYNAMIC_VISCOSITY"].GetDouble();

    mDiffusionCoefficient = mConductivity / (mDensity * mSpecificHeat);
}

const Parameters ManufacturedBodyForceProcess::GetDefaultParameters() const
{
    const Parameters default_parameters( R"(
    {
                "model_part_name"              : "please_specify_model_part_name",
                "impose_inlet_conditions"      : true,
                "write_exact_solution_output"  : true,
                "cavity"                       : false,
                "inlet_model_part_name"        : "inlet_model_part_name",
                "outlet_model_part_name"       : "outlet_model_part_name",
                "central_velocity"             : 1.0,
                "radius"                       : 1.0,
                "length"                       : 1.0,
                "materials_filename"           : "materials_filename.json",
                "add_heat_flux_source"         : false,
                "add_body_force_source"        : true,
                "temperature_gradient"         : 1.0,
                "temperature_inlet"            : 1.0
    }  )" );

    return default_parameters;
}


void ManufacturedBodyForceProcess::Execute()
{
    this->ExecuteInitialize();
    this->ExecuteInitializeSolutionStep();
}

void ManufacturedBodyForceProcess::ExecuteInitialize()
{}

void ManufacturedBodyForceProcess::ExecuteBeforeSolutionLoop()
{
    this->SetSystemProperties();
}

void ManufacturedBodyForceProcess::ExecuteInitializeSolutionStep()
{
    this->SetBodyForceAndMassSource();
}

void ManufacturedBodyForceProcess::ExecuteFinalizeSolutionStep() {
    // Impose no slip conditions manually
    // for (auto it_node = mpNoSlipModelPart->NodesBegin(); it_node != mpNoSlipModelPart->NodesEnd(); it_node++){
    //     double vx = it_node->FastGetSolutionStepValue(VELOCITY_X);
    //     double vy = it_node->FastGetSolutionStepValue(VELOCITY_Y);
    //     double vz = it_node->FastGetSolutionStepValue(VELOCITY_Z);
    //     std::cout << "v = (" << vx << ", " << vy << ", " << vz << ")" << std::endl; 
    // }
}

/* Protected functions ****************************************************/

void ManufacturedBodyForceProcess::SetInitialBodyForceAndMassSource() {
    const double rho = mDensity;
    const double nu = mViscosity / rho;

    double p_inlet = 1e-3;
    double p_outlet = 0.0;
    double temperature_inlet = mInletTemperature;
    for (auto it_node = mrModelPart.NodesBegin(); it_node != mrModelPart.NodesEnd(); it_node++){
        double y = it_node->Y();
        double z = it_node->Z();

        // Read source terms
        double& body_force_x = it_node->FastGetSolutionStepValue(BODY_FORCE_X);
        double& body_force_y = it_node->FastGetSolutionStepValue(BODY_FORCE_Y);
        double& body_force_z = it_node->FastGetSolutionStepValue(BODY_FORCE_Z);
        double& heat_flux = it_node->FastGetSolutionStepValue(HEAT_FLUX);

        // Set body force
        double press_grad = (p_outlet - p_inlet) / mLen;

        // Define the square radius and square radial coord
        double r_coord_2 = pow(y, 2) + pow(z, 2);
        double rad_2 = pow(mRad, 2);
        double r_tilde_2 = r_coord_2 / rad_2;

        double vx, temperature, pressure;  // Exact solutions, used to compute body force and heat flux
        if (mAddBodyForceSource)
        {
            // Impose the body force term with non-trivial velocity (gaussian-like shape)
            double x = it_node->X();
            double const_vel_factor = mv0 / (exp(1.0) - 1.0);  // Factor that normalizes the velocity to [0, v0]
            double body_force_prefactor = 4.0 * nu * const_vel_factor / rad_2;

            vx = const_vel_factor * (exp(1.0 - r_tilde_2) - 1.0);
            pressure = p_inlet + x / mLen * (p_outlet - p_inlet);
            body_force_x = body_force_prefactor * (1. - r_tilde_2) * exp(1. - r_tilde_2);
            body_force_x += press_grad / rho;

            // Impose the heat flux term
            if (mAddHeatFluxSource)
            {
                // Non-null temperature
                temperature = temperature_inlet + mTGrad * (x / mLen);
                heat_flux = (mSpecificHeat * rho) * (vx / mLen) * mTGrad;
            } else {
                // Null temperature
                temperature = temperature_inlet;
                heat_flux = 0.0;
            }
        } else {
            if (mAddHeatFluxSource)
            {
                // No fluid, only temperature diffusion scenario
                double x = it_node->X();

                vx = 0.0;
                pressure = 0.0;
                body_force_x = 0.0;

                temperature = temperature_inlet + mTGrad * pow(x / mLen, 2);  // Quadratic
                heat_flux = -2.0 * mConductivity * mTGrad / pow(mLen, 2);
            } else {
                // Only fluid, no temperature scenario
                double x = it_node->X();
                p_inlet = 4. * mViscosity * mv0 / rad_2 * mLen;

                vx = mv0 * (1. - r_coord_2 / rad_2);
                pressure = p_inlet + x / mLen * (p_outlet - p_inlet);
                body_force_x = 0.0;

                temperature = temperature_inlet;
                heat_flux = 0.0;
            }
        }
        body_force_y = 0.0;
        body_force_z = 0.0;

        // Write the exact solutions
        if (mWriteExactSolutionOutput)
        {
            it_node->FastGetSolutionStepValue(EXACT_FLUX) = heat_flux;
            it_node->FastGetSolutionStepValue(EXACT_BODY_FORCE_X) = body_force_x;
            it_node->FastGetSolutionStepValue(EXACT_BODY_FORCE_Y) = 0.;
            it_node->FastGetSolutionStepValue(EXACT_BODY_FORCE_Z) = 0.;
            it_node->FastGetSolutionStepValue(EXACT_PRESSURE) = pressure;
            it_node->FastGetSolutionStepValue(EXACT_TEMPERATURE) = temperature;
            it_node->FastGetSolutionStepValue(EXACT_VELOCITY_X) = vx;
            it_node->FastGetSolutionStepValue(EXACT_VELOCITY_Y) = 0.;
            it_node->FastGetSolutionStepValue(EXACT_VELOCITY_Z) = 0.;
        }
    }

    // If cavity, impose velocity at the outlet
    if (mCavity) {
        std::cout << "Cavity is ON" << std::endl;

        // Fix the pressure
        for (auto it_node = mpOutletModelPart->NodesBegin(); it_node != mpOutletModelPart->NodesEnd(); it_node++){
            it_node->pGetDof(PRESSURE)->FixDof();
            it_node->FastGetSolutionStepValue(PRESSURE) = 0.0;
            break;
        }
        // Fix the velocity
        for (auto it_node = mpOutletModelPart->NodesBegin(); it_node != mpOutletModelPart->NodesEnd(); it_node++){
            it_node->pGetDof(VELOCITY_X)->FixDof();
            it_node->pGetDof(VELOCITY_Y)->FixDof();
            it_node->pGetDof(VELOCITY_Z)->FixDof();

            double y = it_node->Y();
            double z = it_node->Z();

            double r_coord_2 = pow(y, 2) + pow(z, 2);
            double rad_2 = pow(mRad, 2);
            double r_tilde_2 = r_coord_2 / rad_2;

            double vx_i;
            if (mAddBodyForceSource)
            {
                double const_vel_factor = mv0 / (exp(1.0) - 1.0);  // Factor that normalizes the velocity to [0, v0]
                vx_i = const_vel_factor * (exp(1.0 - r_tilde_2) - 1.0);
            } else {
                vx_i = mv0 * (1. - r_tilde_2);
            }
            
            // Impose conditions at the outlet
            it_node->FastGetSolutionStepValue(VELOCITY_X) = vx_i;
            it_node->FastGetSolutionStepValue(VELOCITY_Y) = 0.0;
            it_node->FastGetSolutionStepValue(VELOCITY_Z) = 0.0;
        }
    } else {
        std::cout << "Cavity is OFF" << std::endl;
        // Impose only the pressure at the outlet!
        bool first_iteration = true;
        for (auto it_node = mpOutletModelPart->NodesBegin(); it_node != mpOutletModelPart->NodesEnd(); it_node++){
            if (first_iteration) {
                it_node->pGetDof(PRESSURE)->FixDof();
                it_node->FastGetSolutionStepValue(PRESSURE) = p_outlet;
                first_iteration = false;
            }
            it_node->pGetDof(VELOCITY_Y)->FixDof();
            it_node->pGetDof(VELOCITY_Z)->FixDof();
            it_node->FastGetSolutionStepValue(VELOCITY_Y) = 0.0;
            it_node->FastGetSolutionStepValue(VELOCITY_Z) = 0.0;
        }
    }

    // Impose the conditions on the inlet
    if (mImposeConditionsInlet) {
        for (auto it_node = mpInletModelPart->NodesBegin(); it_node != mpInletModelPart->NodesEnd(); it_node++){
            double y = it_node->Y();
            double z = it_node->Z();

            double r_coord_2 = pow(y, 2) + pow(z, 2);
            double rad_2 = pow(mRad, 2);
            double r_tilde_2 = r_coord_2 / rad_2;

            double vx_i;
            if (mAddBodyForceSource)
            {
                double const_vel_factor = mv0 / (exp(1.0) - 1.0);  // Factor that normalizes the velocity to [0, v0]
                vx_i = const_vel_factor * (exp(1.0 - r_tilde_2) - 1.0);
            } else {
                vx_i = mv0 * (1. - r_tilde_2);
            }
            
            // Impose conditions at the inlet
            it_node->FastGetSolutionStepValue(VELOCITY_X) = vx_i;
            it_node->FastGetSolutionStepValue(VELOCITY_Y) = 0.0;
            it_node->FastGetSolutionStepValue(VELOCITY_Z) = 0.0;
            it_node->FastGetSolutionStepValue(TEMPERATURE) = temperature_inlet;
        }
    }
        
    // Impose no slip conditions manually
}

void ManufacturedBodyForceProcess::SetBodyForceAndMassSource() {
    this->SetInitialBodyForceAndMassSource();
}

void ManufacturedBodyForceProcess::SetSystemProperties()
{
    (mrModelPart.pGetProperties(1))->SetValue(DENSITY, mDensity);
    (mrModelPart.pGetProperties(1))->SetValue(VISCOSITY, mViscosity);

    block_for_each(mrModelPart.Elements(), [&](Element& rElement){
        rElement.SetProperties(mrModelPart.pGetProperties(1));
    });

    block_for_each(mrModelPart.Nodes(), [&](Node<3>& rNode){
        rNode.FastGetSolutionStepValue(VISCOSITY) = mViscosity;
        rNode.FastGetSolutionStepValue(DENSITY) = mDensity;
    });
}
/* Private functions ****************************************************/

};  // namespace Kratos.