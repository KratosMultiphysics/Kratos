#include "apply_compressible_navier_stokes_boundary_conditions_process.h"
#include "fluid_dynamics_application_variables.h"
#include "utilities/parallel_utilities.h"
#include "utilities/normal_calculation_utils.h"
namespace Kratos
{


ApplyCompressibleNavierStokesBoundaryConditionsProcess::ApplyCompressibleNavierStokesBoundaryConditionsProcess(
    Model& rModel,
    Parameters Parameters)
{
    KRATOS_TRY

    Parameters.ValidateAndAssignDefaults(GetDefaultParameters());
    mpModelPart = &rModel.GetModelPart(Parameters["model_part_name"].GetString());
    
    // Reading flow direction
    const auto flow_direction_name = Parameters["flow_direction_variable"].GetString();

    KRATOS_ERROR_IF_NOT(KratosComponents<VectorVariable>::Has(flow_direction_name))
        << "specified flow direction variable is not a known vector variable." << std::endl;
    
    mpFlowDirectionVariable = & KratosComponents<VectorVariable>::Get(flow_direction_name);

    mRefreshNormalsEveryTimeStep = Parameters["refresh_normals_every_time_step"].GetBool();

    // Reading boundary conditions to enforce
    mSubsonicBCs.reserve(Parameters["subsonic_boundary_conditions"].size());
    mSupersonicBCs.reserve(Parameters["supersonic_boundary_conditions"].size());

    for(auto var_settings : Parameters["subsonic_boundary_conditions"])
    {
        ReadBoundaryCondition(mSubsonicBCs, var_settings);
    }

    for(auto var_settings : Parameters["supersonic_boundary_conditions"])
    {
        ReadBoundaryCondition(mSupersonicBCs, var_settings);
    }

    KRATOS_CATCH("");
}


void ApplyCompressibleNavierStokesBoundaryConditionsProcess::ExecuteInitialize()
{
    NormalCalculationUtils().CalculateUnitNormals<Condition>(*mpModelPart);
}


void ApplyCompressibleNavierStokesBoundaryConditionsProcess::ExecuteInitializeSolutionStep()
{
    if(mRefreshNormalsEveryTimeStep)
    {
        NormalCalculationUtils().CalculateUnitNormals<Condition>(*mpModelPart);
    }

    // Enabling and disabling BC according to the time and active interval
    const double time = mpModelPart->GetProcessInfo().GetValue(TIME);

    for(auto& boundary_condition: mSupersonicBCs)
    {
        boundary_condition.ActivateIfInsideTimeInterval(time);
    }

    for(auto& boundary_condition: mSubsonicBCs)
    {
        boundary_condition.ActivateIfInsideTimeInterval(time);
    }

    // Enforcing boundary conditions
    block_for_each(mpModelPart->Nodes(), [&](NodeType& rNode)
    {
        // Computing mach projected onto normal
        const auto& r_normal = rNode.FastGetSolutionStepValue(NORMAL);
        const auto& r_flow_direction = rNode.FastGetSolutionStepValue(*mpFlowDirectionVariable);

        const double mach = rNode.GetValue(MACH);
        const double flow_direction_norm = norm_2(r_flow_direction);
        
        bool supersonic = false;
        
        if(std::abs(flow_direction_norm) > std::numeric_limits<double>::min())
        {
            const double mach_projection = mach * inner_prod(r_normal, r_flow_direction) / flow_direction_norm;
            supersonic = std::abs(mach_projection) >= 1.0;
        }
        
        // Chosing BC set to enforce according to mach
        const auto& active_bc  = supersonic ? mSupersonicBCs : mSubsonicBCs;

        for(const auto& boundary_condition: active_bc)
        {
            boundary_condition.Enforce(rNode);
        }

    });
}


void ApplyCompressibleNavierStokesBoundaryConditionsProcess::ExecuteFinalizeSolutionStep()
{
    block_for_each(mpModelPart->Nodes(), [&](NodeType& rNode)
    {
        for(const auto& util: mSupersonicBCs)
        {
            rNode.pGetDof(util.GetVariable())->FreeDof();
        }

        for(const auto& util: mSubsonicBCs)
        {
            rNode.pGetDof(util.GetVariable())->FreeDof();
        }
    });
}


ApplyCompressibleNavierStokesBoundaryConditionsProcess::BoundaryConditionUtility::BoundaryConditionUtility(
    const std::string &rVariableName,
    const double Value,
    const IntervalUtility& rIntervalUtility)
    : mValue(Value)
    , mInterval(rIntervalUtility)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(rVariableName))
        << "There is no double variable named " << rVariableName << std::endl;

    mpVariable = &KratosComponents<Variable<double>>::Get(rVariableName);

    KRATOS_CATCH("")
}


void ApplyCompressibleNavierStokesBoundaryConditionsProcess::BoundaryConditionUtility::ActivateIfInsideTimeInterval(const double Time)
{
    if(mInterval.IsInInterval(Time))
    {
        mEnforceInternal = &FixDof;
    }
    else
    {
        mEnforceInternal = nullptr;
    }
}


void ApplyCompressibleNavierStokesBoundaryConditionsProcess::BoundaryConditionUtility::Enforce(NodeType& rNode) const
{
    if(mEnforceInternal)
    {
       mEnforceInternal(*this, rNode);
    }
}


const Variable<double> & ApplyCompressibleNavierStokesBoundaryConditionsProcess::BoundaryConditionUtility::GetVariable() const
{
    return *mpVariable;
}


void ApplyCompressibleNavierStokesBoundaryConditionsProcess::BoundaryConditionUtility::FixDof(
    const BoundaryConditionUtility& rUtility,
    NodeType& rNode)
{
    rNode.GetSolutionStepValue(*rUtility.mpVariable) = rUtility.mValue;
    rNode.pGetDof(*rUtility.mpVariable)->FixDof();
}


void ApplyCompressibleNavierStokesBoundaryConditionsProcess::ReadBoundaryCondition(
    std::vector<BoundaryConditionUtility> & rBCList,
    Parameters Parameters)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(Parameters.Has("variable_name"))
        << "Missing string field 'variable_name' in boundary condition:\n" << Parameters  << std::endl;
    
    KRATOS_ERROR_IF_NOT(Parameters.Has("value"))
        << "Missing double (or vector) field 'value' in boundary condition:\n" << Parameters  << std::endl;

    // Reading interval
    IntervalUtility interval_utility{Parameters};
    
    // Reading value and acting depending on if it's vector or double
    if(Parameters["value"].IsNumber())
    {
        rBCList.emplace_back(
            Parameters["variable_name"].GetString(),
            Parameters["value"].GetDouble(),
            interval_utility
        );
        return;
    }

    if(Parameters["value"].IsVector())
    {
        KRATOS_ERROR_IF_NOT(Parameters.Has("constrained"))
            << "Missing boolean array field 'constrained' in boundary condition:\n"
            << Parameters  << std::endl;

        const auto values = Parameters["value"].GetVector();
        const auto constraints = Parameters["constrained"];

        const std::size_t n_components = values.size();
        const std::size_t n_constraints = constraints.size();

        KRATOS_ERROR_IF(n_constraints != n_components)
            << "The number of values specified must be the same as number of constraints specified:\n"
            << Parameters << std::endl;

        KRATOS_ERROR_IF(n_components > 3)
            << "Allowed vector variables are at most 3-dimensional:\n"
            << Parameters << std::endl;

        std::string variable_component_name = Parameters["variable_name"].GetString() + "_X";
        
        for(std::size_t i=0; i<n_components; ++i)
        {
            if(constraints[i].GetBool())
            {
                variable_component_name.back() = "XYZ"[i];
                rBCList.emplace_back(variable_component_name, values[i], interval_utility);
            }
        }

        return;
    }

    KRATOS_ERROR << "ApplyCompressibleNavierStokesBoundaryConditionsProcess supports only double and vector variables:\n" << Parameters << std::endl;

    KRATOS_CATCH("");
}


const Parameters ApplyCompressibleNavierStokesBoundaryConditionsProcess::GetDefaultParameters() const
{
    return Parameters(R"(
    {
        "model_part_name" : "main_model_part",
        "flow_direction_variable" : "PLEASE SPECIFY FLOW DIRECTION VARIABLE",
        "refresh_normals_every_time_step" : false,
        "subsonic_boundary_conditions" : [ ],
        "supersonic_boundary_conditions" : [ ]
    }
    )");
    /* Expected boundary condition format:
    [
        {                                     // For scalar variables
            "variable_name" : "DENSITY",   <--- The variable to fix. Must be a "double" variable
            "interval" : [0, "End"],       <--- Time interval for it to be fixed
            "value" : 0.0                  <--- The value to fix it to
        },
        {                                    // For vector variables
            "variable_name" : "MOMENTUM",         <--- Must be a "vector with components" variable.
            "interval" : [0.0, 15.3],             <--- Between 1 and 3 components.
            "value" : [0.0, 1.0, 5.0],
            "constrained" : [false, true, false]  <-- What components to fix. Others will be ignored.
                                                      Must have the same number of components as value.
        }
    ]
    */
}


std::string ApplyCompressibleNavierStokesBoundaryConditionsProcess::Info() const
{
    return "ApplyCompressibleNavierStokesBoundaryConditionsProcess";
}


void ApplyCompressibleNavierStokesBoundaryConditionsProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "ApplyCompressibleNavierStokesBoundaryConditionsProcess";
}


/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const ApplyCompressibleNavierStokesBoundaryConditionsProcess& rThis)
{
    rThis.PrintData(rOStream);
    return rOStream;
}


} // namespace Kratos
