//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Juan I. Camarotti
//
//

// System includes

// External includes

// Project includes

// Application includes
#include "compute_aerodynamic_coefficients_process.h"
#include "fluid_dynamics_application_variables.h"
#include "utilities/variable_utils.h"


namespace Kratos
{

ComputeAerodynamicCoefficientsProcess::ComputeAerodynamicCoefficientsProcess(
    Model &rModel,
    Parameters Params)
    : Process(),
      mrModelPart(rModel.GetModelPart(Params["model_part_name"].GetString()))
{
    // Check default settings
    KRATOS_TRY
    
    Params.ValidateAndAssignDefaults(GetDefaultParameters());

    ReadFreestreamValues(Params);

    KRATOS_CATCH("")
}


const Parameters ComputeAerodynamicCoefficientsProcess::GetDefaultParameters() const
{
    return Parameters( R"(
    {
        "model_part_name"     : "PLEASE_PROVIDE_A_MODELPART_NAME",
        "reference_surface"   : 0.0,
        "reference_chord"     : 1.0,
        "reference_span"      : 1.0,
        "moment_reference_point" : [0.0, 0.0, 0.0],
        "freestream_dynamic_pressure" : 1.0,
        "angle_of_attack"     : 0.0,
        "sideslip_angle"  : 0.0
    })" );
}

void ComputeAerodynamicCoefficientsProcess::ReadFreestreamValues(const Parameters& rParams)
{
    constexpr double tol = 1e-12;

    mReference_Surface = rParams["reference_surface"].GetDouble();
    KRATOS_ERROR_IF(std::abs(mReference_Surface) <= tol)
        << "Invalid value for 'reference_surface' = " << mReference_Surface << ".";

    mReference_Chord = rParams["reference_chord"].GetDouble();
    KRATOS_ERROR_IF(std::abs(mReference_Chord) <= tol)
        << "Invalid value for 'reference_chord' = " << mReference_Chord << ".";

    mReference_Span = rParams["reference_span"].GetDouble();
    KRATOS_ERROR_IF(std::abs(mReference_Span) <= tol)
        << "Invalid value for 'reference_span' = " << mReference_Span << ".";

    mQInf = rParams["freestream_dynamic_pressure"].GetDouble();
    KRATOS_ERROR_IF(std::abs(mQInf) <= tol)
        << "Invalid value for 'freestream_dynamic_pressure' = " << mQInf << ".";

    mAngleOfAttack = rParams["angle_of_attack"].GetDouble();
    mSideslipAngle  = rParams["sideslip_angle"].GetDouble();

    const auto& r_mrp = rParams["moment_reference_point"];
    KRATOS_ERROR_IF_NOT(r_mrp.IsArray() && r_mrp.size() == 3)
        << "'moment_reference_point' must be an array of size 3.";
    mMomentReferencePoint[0] = r_mrp[0].GetDouble();
    mMomentReferencePoint[1] = r_mrp[1].GetDouble();
    mMomentReferencePoint[2] = r_mrp[2].GetDouble();
}


void ComputeAerodynamicCoefficientsProcess::ExecuteBeforeOutputStep()
{
    Execute();
}


void ComputeAerodynamicCoefficientsProcess::Execute()
{
    KRATOS_TRY;

    // Build wind axes from alpha (AoA) and beta (sideslip)
    const double alpha = mAngleOfAttack * Globals::Pi / 180.0;
    const double beta  = mSideslipAngle  * Globals::Pi / 180.0;

    // Wind-x (drag axis): freestream direction (unit)
    array_1d<double,3> e_drag = ZeroVector(3);
    e_drag[0] =  std::cos(alpha) * std::cos(beta);
    e_drag[1] =  std::sin(beta);
    e_drag[2] =  std::sin(alpha) * std::cos(beta);

    const double n_drag = norm_2(e_drag);
    KRATOS_ERROR_IF(n_drag < 1e-15) << "Invalid wind axis: e_drag near zero.";
    e_drag /= n_drag;

    // Build a reference "up" vector to construct orthonormal basis
    array_1d<double,3> up = ZeroVector(3);
    up[0] = 0.0; up[1] = 0.0; up[2] = 1.0;
    if (std::abs(inner_prod(e_drag, up)) > 0.99) {
        up[0] = 0.0; up[1] = 1.0; up[2] = 0.0;
    }

    // Wind-y (side axis)
    array_1d<double,3> e_side = MathUtils<double>::CrossProduct(up, e_drag);
    const double n_side = norm_2(e_side);
    KRATOS_ERROR_IF(n_side < 1e-15) << "Invalid wind axis: e_side near zero.";
    e_side /= n_side;

    // Wind-z (lift axis)
    array_1d<double,3> e_lift = MathUtils<double>::CrossProduct(e_drag, e_side);
    const double n_lift = norm_2(e_lift);
    KRATOS_ERROR_IF(n_lift < 1e-15) << "Invalid wind axis: e_lift near zero.";
    e_lift /= n_lift;

    // Accumulate aerodynamic force and moment from nodal reactions
    array_1d<double,3> F_aero = ZeroVector(3);
    array_1d<double,3> M_aero = ZeroVector(3);

    for (auto& r_node : mrModelPart.Nodes())
    {
        KRATOS_ERROR_IF_NOT(r_node.SolutionStepsDataHas(REACTION))
            << "Node " << r_node.Id() << " does not have REACTION in SolutionStepData.";

        const array_1d<double,3>& r_R = r_node.FastGetSolutionStepValue(REACTION);
        const array_1d<double,3> f_node = -r_R; 

        F_aero += f_node;

        array_1d<double,3> r = r_node.Coordinates() - mMomentReferencePoint;
        M_aero += MathUtils<double>::CrossProduct(r, f_node);
    }

    // Projections in wind axes
    const double D = inner_prod(F_aero, e_drag);   // drag
    const double Y = inner_prod(F_aero, e_side);   // side force
    const double L = inner_prod(F_aero, e_lift);   // lift

    const double Mx = inner_prod(M_aero, e_drag);  // moment about wind-x (roll in wind frame)
    const double My = inner_prod(M_aero, e_side);  // moment about wind-y (pitch in wind frame)
    const double Mz = inner_prod(M_aero, e_lift);  // moment about wind-z (yaw in wind frame)

    // Denominator for non-dimensional coefficients
    const double denom_F = mQInf * mReference_Surface;

    KRATOS_ERROR_IF(std::abs(denom_F) < 1e-15) << "Invalid denominator q_inf*S_ref.";

    const double C_D = D / denom_F;
    const double C_Y = Y / denom_F;
    const double C_L = L / denom_F;

    const double C_l = Mx / (denom_F * mReference_Span);   
    const double C_m = My / (denom_F * mReference_Chord);  
    const double C_n = Mz / (denom_F * mReference_Span);   

    // Store on ModelPart
    mrModelPart.SetValue(LIFT_COEFFICIENT, C_L);
    mrModelPart.SetValue(DRAG_COEFFICIENT, C_D);
    mrModelPart.SetValue(LATERAL_FORCE_COEFFICIENT, C_Y);

    mrModelPart.SetValue(ROLLING_MOMENT_COEFFICIENT, C_l);
    mrModelPart.SetValue(PITCHING_MOMENT_COEFFICIENT, C_m);
    mrModelPart.SetValue(YAWING_MOMENT_COEFFICIENT, C_n);

    std::cout << "Aero coeffs from nodal reactions | "
              << "alpha=" << mAngleOfAttack << " deg, beta=" << mSideslipAngle << " deg : "
              << "CL=" << C_L << ", CD=" << C_D << ", CY=" << C_Y
              << ", Cl=" << C_l << ", Cm=" << C_m << ", Cn=" << C_n
              << std::endl;

    KRATOS_CATCH("");
}


}