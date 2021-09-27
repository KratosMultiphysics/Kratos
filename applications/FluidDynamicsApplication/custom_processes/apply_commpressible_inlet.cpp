#include "apply_compressible_inlet_process.h"
#include "fluid_dynamics_application_variables.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{

ApplyCompressibleInlet::ApplyCompressibleInlet(Model& rModel, Parameters& rParameters)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());
    mpModelPart = &rModel.GetModelPart(rParameters["model_part_name"].GetString());

    // Reading boundary conditions to enforce
    mSubsonicBCs.reserve(rParameters["subsonic_boundary_conditions"].size());
    mSupersonicBCs.reserve(rParameters["supersonic_boundary_conditions"].size());

    for(auto & var_settings : rParameters["subsonic_boundary_conditions"])
    {
        ReadBoundaryCondition(mSubsonicBCs, var_settings);
    }

    for(auto & var_settings : rParameters["supersonic_boundary_conditions"])
    {
        ReadBoundaryCondition(mSupersonicBCs, var_settings);
    }

    KRATOS_CATCH("");
}

void ApplyCompressibleInlet::ExecuteInitializeSolutionStep()
{
    // Enabling and disabling BC according to the time and active interval
    const double time = mpModelPart->GetProcessInfo().GetValue(TIME);

    for(auto & boundary_condition: mSupersonicBCs)
    {
        boundary_condition.ActivateIfInsideInterval(time);
    }

    for(auto & boundary_condition: mSubsonicBCs)
    {
        boundary_condition.ActivateIfInsideInterval(time);
    }

    // Enforcing boundary conditions
    block_for_each(mpModelPart->Nodes(), [&](NodeType & rNode)
    {
        const bool supersonic = rNode.GetValue(MACH) >= 1.0;
        
        const auto & bc_list = supersonic ? mSupersonicBCs : mSubsonicBCs;

        for(const auto & boundary_condition: bc_list)
        {
            boundary_condition.Enforce(rNode);
        }
    });
}


ApplyCompressibleInlet::BoundaryConditionUtility::
BoundaryConditionUtility(const std::string & variable_name, const double Value, const double Start, const double End)
    : mValue(Value), mStart(Start), mEnd(End)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(variable_name));

    mpVariable = &KratosComponents<Variable<double>>::Get(variable_name);

    KRATOS_CATCH("")
}

void ApplyCompressibleInlet::BoundaryConditionUtility::
ActivateIfInsideInterval(const double time)
{
    if(time >= mStart && time < mEnd)
    {
        mEnforceInternal = &EnforceActive;
    }
    else
    {
        mEnforceInternal = &EnforcePassive;
    }
}


void ApplyCompressibleInlet::BoundaryConditionUtility::
Enforce(NodeType & rNode) const
{
    mEnforceInternal(*this, rNode);
}

void ApplyCompressibleInlet::BoundaryConditionUtility::
EnforceActive(const BoundaryConditionUtility & rUtility, NodeType & rNode)
{
    rNode.GetSolutionStepValue(*rUtility.mpVariable) = rUtility.mValue;
    rNode.pGetDof(*rUtility.mpVariable)->FixDof();
}

void ApplyCompressibleInlet::BoundaryConditionUtility::
EnforcePassive(const BoundaryConditionUtility & rUtility, NodeType & rNode)
{
}


constexpr char ApplyCompressibleInlet::IndexToAxis(const unsigned int i)
{
    return static_cast<char>(static_cast<unsigned int>('X') + i);
}


double ReadTime(const Parameters & raw_data)
{
    KRATOS_TRY
    
    if(raw_data.IsDouble())
    {
        return raw_data.GetDouble();
    }
    
    if(raw_data.IsString() && raw_data.GetString() == "End")
    {
        return std::numeric_limits<double>::max();
    }

    KRATOS_ERROR << "Interval must contain either a number of the string 'End'" << std::endl;

    KRATOS_CATCH("")
}

void ApplyCompressibleInlet::ReadBoundaryCondition(std::vector<BoundaryConditionUtility> & rBCList, Parameters Parameters)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(Parameters.Has("variable"))
        << "Missing string field 'variable' in boundary condition:\n" << Parameters  << std::endl;
    
    KRATOS_ERROR_IF_NOT(Parameters.Has("value"))
        << "Missing double field 'value' in boundary condition:\n" << Parameters  << std::endl;

    // Reading interval
    double interval_start = 0.0;
    double interval_end   = std::numeric_limits<double>::max();
    if(Parameters.Has("interval"))
    {
        KRATOS_ERROR_IF_NOT(Parameters["interval"].size() == 2)
            << "Field 'interval' is expected to contain an array with two entries [start, end]:\n" << Parameters << std::endl;
        interval_start = ReadTime(Parameters["interval"][0]);
        interval_end = ReadTime(Parameters["interval"][1]);
    }

    
    // Reading value and acting depending on if it's vector or double
    if(Parameters["value"].IsDouble())
    {
        rBCList.emplace_back(
            Parameters["variable"].GetString(),
            Parameters["value"].GetDouble(),
            interval_start,
            interval_end
        );
        return;
    }

    if(Parameters["value"].IsVector())
    {
        const auto values = Parameters["value"].GetVector();

        KRATOS_ERROR_IF_NOT(Parameters.Has("constrained"))
            << "Missing boolean array field 'constrained' in boundary condition:\n" << Parameters  << std::endl;

        KRATOS_ERROR_IF_NOT(values.size() == Parameters["constrained"].size())
            << "The number of values specified must be the same as number of constraints specified:\n" << Parameters << std::endl;

        for(std::size_t i=0; i<Parameters["constrained"].size(); i++)
        {
            if(Parameters["constrained"][i].GetBool())
            {
                std::string variable_name = Parameters["variable"].GetString();
                variable_name += "_" + IndexToAxis(i);
                rBCList.emplace_back(variable_name, values[i], interval_start, interval_end);
            }
        }

        return;
    }

    KRATOS_ERROR << "ApplyCompressibleInlet supports only double and vector variables:\n" << Parameters << std::endl;

    KRATOS_CATCH("");
}


const Parameters ApplyCompressibleInlet::GetDefaultParameters() const
{
    return Parameters(R"(
    {

        "model_part_name" : "main_model_part",
        "subsonic_boundary_conditions" : [ ]
        "supersonic_boundary_conditions"< : [ ]
    }
    )");
    /* Expected boundary condition format:
    [
        {
            "variable" : "DENSITY",
            "value" : 0.0
            "interval" : [0, "End"]
        },
        {
            "variable" : "MOMENTUM",
            "value" : [0.0, 1.0, 5.0],
            "constrained" : [false, true, false]
            "interval" : [0.0, 15.3]
        }
    ]
    */
}

std::string ApplyCompressibleInlet::Info() const
{
    return "ApplyCompressibleInlet";
}

void ApplyCompressibleInlet::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "ApplyCompressibleInlet";
}

void ApplyCompressibleInlet::PrintData(std::ostream& rOStream) const
{
}

/* External functions *****************************************************/

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const ApplyCompressibleInlet& rThis)
{
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos
