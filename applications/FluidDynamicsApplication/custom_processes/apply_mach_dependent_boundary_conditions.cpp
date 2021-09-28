#include "apply_mach_dependent_boundary_conditions.h"
#include "fluid_dynamics_application_variables.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{


ApplyMachDependentBoundaryConditions::ApplyMachDependentBoundaryConditions(
    Model& rModel,
    Parameters Parameters)
{
    KRATOS_TRY

    Parameters.ValidateAndAssignDefaults(GetDefaultParameters());
    mpModelPart = &rModel.GetModelPart(Parameters["model_part_name"].GetString());

    // Reading boundary conditions to enforce
    mSubsonicBCs.reserve(Parameters["subsonic_boundary_conditions"].size());
    mSupersonicBCs.reserve(Parameters["supersonic_boundary_conditions"].size());

    for(auto & var_settings : Parameters["subsonic_boundary_conditions"])
    {
        ReadBoundaryCondition(mSubsonicBCs, var_settings);
    }

    for(auto & var_settings : Parameters["supersonic_boundary_conditions"])
    {
        ReadBoundaryCondition(mSupersonicBCs, var_settings);
    }

    KRATOS_CATCH("");
}


void ApplyMachDependentBoundaryConditions::ExecuteInitializeSolutionStep()
{
    // Enabling and disabling BC according to the time and active interval
    const double time = mpModelPart->GetProcessInfo().GetValue(TIME);

    for(auto & boundary_condition: mSupersonicBCs)
    {
        boundary_condition.ActivateIfInsideTimeInterval(time);
    }

    for(auto & boundary_condition: mSubsonicBCs)
    {
        boundary_condition.ActivateIfInsideTimeInterval(time);
    }

    // Enforcing boundary conditions
    block_for_each(mpModelPart->Nodes(), [&](NodeType & rNode)
    {
        const bool supersonic = rNode.GetValue(MACH) >= 1.0;
        
        const auto & active_bc  = supersonic ? mSupersonicBCs : mSubsonicBCs;
        const auto & passive_bc = supersonic ? mSubsonicBCs   : mSupersonicBCs;

        for(const auto & boundary_condition: active_bc)
        {
            boundary_condition.Enforce(rNode);
        }

        for(const auto & boundary_condition: passive_bc)
        {
            boundary_condition.Release(rNode);
        }

    });
}


ApplyMachDependentBoundaryConditions::BoundaryConditionUtility::
    BoundaryConditionUtility(const std::string & variable_name,
                            const double Value,
                            const IntervalUtility & rIntervalUtility)
    : mValue(Value), mInterval(rIntervalUtility)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(variable_name))
        << "There is no double variable named " << variable_name << std::endl;

    mpVariable = &KratosComponents<Variable<double>>::Get(variable_name);

    KRATOS_CATCH("")
}


void ApplyMachDependentBoundaryConditions::BoundaryConditionUtility::
    ActivateIfInsideTimeInterval(const double time)
{
    if(mInterval.IsInInterval(time))
    {
        mEnforceInternal = &FixDof;
    }
    else
    {
        mEnforceInternal = &FreeDof;
    }
}


void ApplyMachDependentBoundaryConditions::BoundaryConditionUtility::
    Enforce(NodeType & rNode) const
{
    mEnforceInternal(*this, rNode);
}


void ApplyMachDependentBoundaryConditions::BoundaryConditionUtility::
    Release(NodeType & rNode) const
{
    FreeDof(*this, rNode);
}


void ApplyMachDependentBoundaryConditions::BoundaryConditionUtility::FixDof(
    const BoundaryConditionUtility & rUtility,
    NodeType & rNode)
{
    rNode.GetSolutionStepValue(*rUtility.mpVariable) = rUtility.mValue;
    rNode.pGetDof(*rUtility.mpVariable)->FixDof();
}


void ApplyMachDependentBoundaryConditions::BoundaryConditionUtility::FreeDof(
    const BoundaryConditionUtility & rUtility,
    NodeType & rNode)
{
    rNode.pGetDof(*rUtility.mpVariable)->FreeDof();
}


void ApplyMachDependentBoundaryConditions::ReadBoundaryCondition(
    std::vector<BoundaryConditionUtility> & rBCList,
    Parameters Parameters)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(Parameters.Has("variable"))
        << "Missing string field 'variable' in boundary condition:\n" << Parameters  << std::endl;
    
    KRATOS_ERROR_IF_NOT(Parameters.Has("value"))
        << "Missing double (or vector) field 'value' in boundary condition:\n" << Parameters  << std::endl;

    // Reading interval
    IntervalUtility interval_utility{Parameters};
    
    // Reading value and acting depending on if it's vector or double
    if(Parameters["value"].IsDouble())
    {
        rBCList.emplace_back(
            Parameters["variable"].GetString(),
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
        const auto & constraints = Parameters["constrained"];

        const std::size_t n_components = values.size();
        const std::size_t n_constraints = constraints.size();

        KRATOS_ERROR_IF(n_constraints != n_components)
            << "The number of values specified must be the same as number of constraints specified:\n"
            << Parameters << std::endl;

        KRATOS_ERROR_IF(n_components > 3)
            << "Allowed vector variables are at most 3-dimensional:\n"
            << Parameters << std::endl;

        std::string variable_component_name = Parameters["variable"].GetString() + "_X";
        
        for(std::size_t i=0; i<n_components; i++)
        {
            if(constraints[i].GetBool())
            {
                variable_component_name.back() = "XYZ"[i];
                rBCList.emplace_back(variable_component_name, values[i], interval_utility);
            }
        }

        return;
    }

    KRATOS_ERROR << "ApplyMachDependentBoundaryConditions supports only double and vector variables:\n" << Parameters << std::endl;

    KRATOS_CATCH("");
}


const Parameters ApplyMachDependentBoundaryConditions::GetDefaultParameters() const
{
    return Parameters(R"(
    {
        "model_part_name" : "main_model_part",
        "subsonic_boundary_conditions" : [ ],
        "supersonic_boundary_conditions" : [ ]
    }
    )");
    /* Expected boundary condition format:
    [
        { // For scalar variables
            "variable" : "DENSITY",
            "value" : 0.0,
            "interval" : [0, "End"]
        },
        { // For vector variables
            "variable" : "MOMENTUM",
            "value" : [0.0, 1.0, 5.0],
            "constrained" : [false, true, false],
            "interval" : [0.0, 15.3]
        }
    ]
    */
}


std::string ApplyMachDependentBoundaryConditions::Info() const
{
    return "ApplyCompressibleInlet";
}


void ApplyMachDependentBoundaryConditions::
    PrintInfo(std::ostream& rOStream) const
{
    rOStream << "ApplyCompressibleInlet";
}


void ApplyMachDependentBoundaryConditions::
    PrintData(std::ostream& rOStream) const
{
}


/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const ApplyMachDependentBoundaryConditions& rThis)
{
    rThis.PrintData(rOStream);
    return rOStream;
}


} // namespace Kratos
