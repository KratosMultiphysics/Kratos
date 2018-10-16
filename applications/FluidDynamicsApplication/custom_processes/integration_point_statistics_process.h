//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#if !defined(KRATOS_INTEGRATION_POINT_STATISTICS_PROCESS_H_INCLUDED)
#define KRATOS_INTEGRATION_POINT_STATISTICS_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "custom_utilities/statistics_record.h"

namespace Kratos
{
///@addtogroup ApplicationNameApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
  */
class IntegrationPointStatisticsProcess: public Process
{
public:
///@name Type Definitions
///@{

/// Pointer definition of IntegrationPointStatisticsProcess
KRATOS_CLASS_POINTER_DEFINITION(IntegrationPointStatisticsProcess);

typedef std::map< std::string, StatisticsSampler::Pointer > StatisticsDictionary;

///@}
///@name Life Cycle
///@{

/// Constructor.
IntegrationPointStatisticsProcess(ModelPart& rModelPart, Kratos::Parameters Parameters):
    Process(),
    mrModelPart(rModelPart),
    mParameters(Parameters),
    mStartTime(0.0)
{}

/// Default constructor.
IntegrationPointStatisticsProcess(Model& rModel, Kratos::Parameters Parameters):
    Process(),
    mrModelPart(rModel.GetModelPart(Parameters["model_part_name"].GetString())),
    mParameters(Parameters),
    mStartTime(0.0)
{}

/// Destructor.
~IntegrationPointStatisticsProcess() override
{}

///@}
///@name Operators
///@{

///@}
///@name Operations
///@{

void ExecuteInitialize() override
{
    KRATOS_TRY;

    // Create a new container for statistics
    StatisticsRecord::Pointer p_turbulence_statistics = Kratos::make_shared<StatisticsRecord>();

    // Build statistics records
    CreateStatisticsFromInput(p_turbulence_statistics);

    // Initialize STATISTICS_CONTAINER in ProcessInfo
    p_turbulence_statistics->InitializeStorage(mrModelPart.Elements());
    mrModelPart.GetProcessInfo().SetValue(STATISTICS_CONTAINER,p_turbulence_statistics);

    KRATOS_CATCH("");
}

/// Update statistics.
void ExecuteFinalizeSolutionStep() override
{
    if (mrModelPart.GetProcessInfo()[TIME] >= mStartTime)
    {
        auto p_turbulence_statistics = mrModelPart.GetProcessInfo().GetValue(STATISTICS_CONTAINER);
        p_turbulence_statistics->SampleIntegrationPointResults(mrModelPart);
    }
}

/// Output simulation results to file.
void ExecuteFinalize() override
{
    auto p_turbulence_statistics = mrModelPart.GetProcessInfo().GetValue(STATISTICS_CONTAINER);
    p_turbulence_statistics->PrintToFile(mrModelPart, mOutputFileName);
}

///@}
///@name Access
///@{

///@}
///@name Inquiry
///@{

///@}
///@name Input and output
///@{

/// Turn back information as a string.
std::string Info() const override
{
    std::stringstream buffer;
    buffer << "IntegrationPointStatisticsProcess";
    return buffer.str();
}

/// Print information about this object.
void PrintInfo(std::ostream &rOStream) const override { rOStream << "IntegrationPointStatisticsProcess"; }

/// Print object's data.
void PrintData(std::ostream &rOStream) const override {}

///@}
///@name Friends
///@{

///@}

protected:
///@name Protected static Member Variables
///@{

///@}
///@name Protected member Variables
///@{

///@}
///@name Protected Operators
///@{

///@}
///@name Protected Operations
///@{

void CreateStatisticsFromInput(StatisticsRecord::Pointer pRecordedStatistics)
{
    // Validate parameters
    Kratos::Parameters default_parameters = Kratos::Parameters(R"({
        "statistics" : [],
        "output_file_name": "statistics",
        "model_part_name": ""
    })");

    mParameters.ValidateAndAssignDefaults(default_parameters);

    mOutputFileName = mParameters["output_file_name"].GetString();

    StatisticsDictionary defined_statistics;

    for (unsigned int i = 0; i < mParameters["statistics"].size(); i++)
    {
        Kratos::Parameters settings = mParameters["statistics"][i];
        KRATOS_ERROR_IF_NOT( settings.Has("type") && settings["type"].IsString() )
        << "Element " << i << " in the list of statistics passed to IntegrationPointStatitsticsProcess"
        << " has no \"type\" (string) attribute defined." << std::endl;

        std::string statistic_type = settings["type"].GetString();

        if ( statistic_type == "average" )
        {
            pRecordedStatistics->AddResult(CreateAverageSampler(settings,defined_statistics));
        }
        else if ( statistic_type == "variance" )
        {
            pRecordedStatistics->AddHigherOrderStatistic(CreateVarianceSampler(settings,defined_statistics));
        }
        else if ( statistic_type == "third_order_moment" )
        {
//            pRecordedStatistics->AddHigherOrderStatistic(CreateThirdOrderSampler(settings,defined_statistics));
        }
        else
        {
            KRATOS_ERROR << "Unknown string \"" << statistic_type << "\" passed as \"type\" argument in statistics definition." << std::endl;
        }
    }

}

StatisticsSampler::Pointer CreateAverageSampler(
    Kratos::Parameters Parameters,
    StatisticsDictionary& rDefinedStatistics) const
{
    Kratos::Parameters default_parameters(R"({
        "type" : "average",
        "variable": ""
    })");

    Parameters.ValidateAndAssignDefaults(default_parameters);
    std::string variable_name = Parameters["variable"].GetString();
    std::string type = Parameters["type"].GetString();
    KRATOS_ERROR_IF_NOT(type == "average") << "Trying to define an average statistic of unsupported type " << type << "." << std::endl;

    KRATOS_ERROR_IF(rDefinedStatistics.find(variable_name) != rDefinedStatistics.end())
    << "Duplicate definition of an average for " << variable_name << "." << std::endl;

    StatisticsSampler::Pointer new_statistic;
    if (KratosComponents<Variable<double>>::Has(variable_name))
    {
        // build double variable sampler
        Variable<double> variable = KratosComponents<Variable<double>>::Get(variable_name);
        auto value_getter = Kratos::Internals::MakeSamplerAtLocalCoordinate::ValueGetter(variable);
        new_statistic = Kratos::make_shared<ScalarAverageSampler>(value_getter,variable_name);
    }
    else if (KratosComponents<Variable<array_1d<double,3>>>::Has(variable_name))
    {
        // build vector sampler
        Variable<array_1d<double,3>> variable = KratosComponents<Variable<array_1d<double,3>>>::Get(variable_name);
        auto value_getter = Kratos::Internals::MakeSamplerAtLocalCoordinate::ValueGetter(variable);
        std::vector<std::string> tags;
        tags.push_back(std::string(variable_name+"_X"));
        tags.push_back(std::string(variable_name+"_Y"));
        tags.push_back(std::string(variable_name+"_Z"));
        new_statistic = Kratos::make_shared<VectorAverageSampler<array_1d<double,3>>>(value_getter,3,tags);
    }
    else
    {
        KRATOS_ERROR
        << "Trying to define an average statistic for variable " << variable_name
        << " which is not a variable of a supported type." << std::endl;
    }

    rDefinedStatistics[variable_name] = new_statistic;
    return new_statistic;
}

StatisticsSampler::Pointer CreateVarianceSampler(
    Kratos::Parameters Parameters,
    StatisticsDictionary& rDefinedStatistics) const
{
    Kratos::Parameters default_parameters(R"({
        "type" : "",
        "variables": []
    })");

    Parameters.ValidateAndAssignDefaults(default_parameters);
    KRATOS_ERROR_IF(Parameters["variables"].size() < 1 || Parameters["variables"].size() > 2)
    << "Unexpected number of arguments when reading \"variables\" list argument."
    << "Expected 1 or 2 values, got " << Parameters["variables"].size() << std::endl;
    std::string type = Parameters["type"].GetString();

    StatisticsSampler::Pointer new_statistic;
    if (Parameters["variables"].size() == 1)
    {
        // symmetric variance
        std::string variable_name = Parameters["variables"][0].GetString();
        StatisticsDictionary::iterator it_average = rDefinedStatistics.find(variable_name);
        KRATOS_ERROR_IF(it_average == rDefinedStatistics.end())
        << "Trying to define a variance for " << variable_name
        << " but no average has been defined for this variable." << std::endl;

        new_statistic = Kratos::make_shared<SymmetricVarianceSampler>(it_average->second);
    }
    else // size == 2
    {
        // complete or componentwise variance
        std::string first_argument = Parameters["variables"][0].GetString();
        std::string first_argument_base;
        unsigned int first_argument_index;
        bool first_argument_is_component = ProcessComponent(first_argument,first_argument_base,first_argument_index);

        std::string second_argument = Parameters["variables"][1].GetString();
        std::string second_argument_base;
        unsigned int second_argument_index;
        bool second_argument_is_component = ProcessComponent(second_argument,second_argument_base,second_argument_index);

        StatisticsSampler::Pointer first_statistic;
        if (first_argument_is_component)
        {
            auto found = rDefinedStatistics.find(first_argument_base);
            KRATOS_ERROR_IF(found == rDefinedStatistics.end())
            << "Trying to record variance for " << first_argument << " and " << second_argument
            << " but no average was declared for " << first_argument_base << "." << std::endl;
            first_statistic = found->second;

        }
        else
        {
            auto found = rDefinedStatistics.find(first_argument);
            KRATOS_ERROR_IF(found == rDefinedStatistics.end())
            << "Trying to record variance for " << first_argument << " and " << second_argument
            << " but no average was declared for " << first_argument << "." << std::endl;
            first_statistic = found->second;
        }

        StatisticsSampler::Pointer second_statistic;
        if (second_argument_is_component)
        {
            auto found = rDefinedStatistics.find(second_argument_base);
            KRATOS_ERROR_IF(found == rDefinedStatistics.end())
            << "Trying to record variance for " << first_argument << " and " << second_argument
            << " but no average was declared for " << second_argument_base << "." << std::endl;
            second_statistic = found->second;
        }
        else
        {
            auto found = rDefinedStatistics.find(second_argument);
            KRATOS_ERROR_IF(found == rDefinedStatistics.end())
            << "Trying to record variance for " << first_argument << " and " << second_argument
            << " but no average was declared for " << second_argument << "." << std::endl;
            second_statistic = found->second;
        }

        if (first_argument_is_component || second_argument_is_component)
        {
            // componentwise statistic
            new_statistic = Kratos::make_shared<ComponentwiseVarianceSampler>(
                first_statistic,first_argument_index,
                second_statistic,second_argument_index);
        }
        else
        {
            // full variable statistic
            new_statistic = Kratos::make_shared<VarianceSampler>(
                first_statistic, second_statistic);
        }
    }

    //rDefinedStatistics[statistic_name] = new_statistic;
    return new_statistic;
}

///@}
///@name Protected  Access
///@{

///@}
///@name Protected Inquiry
///@{

///@}
///@name Protected LifeCycle
///@{

///@}

private:
///@name Static Member Variables
///@{

///@}
///@name Member Variables
///@{

ModelPart& mrModelPart;

Kratos::Parameters mParameters;

double mStartTime;

std::string mOutputFileName;

///@}
///@name Private Operators
///@{

///@}
///@name Private Operations
///@{

bool ProcessComponent(
    const std::string& rInputName,
    std::string& rBaseVariableName,
    unsigned int& rComponentIndex) const
{
    bool is_component = false;
    rComponentIndex = 0; // can be used also if variable is not a component

    const std::string x_suffix = std::string("_X");
    const std::string y_suffix = std::string("_Y");
    const std::string z_suffix = std::string("_Z");

    if( StringEndsWith(rInputName,x_suffix) )
    {
        is_component = true;
        rBaseVariableName = rInputName.substr(0, rInputName.length() - 2);
        rComponentIndex = 0;
    }
    else if ( StringEndsWith(rInputName,y_suffix) )
    {
        is_component = true;
        rBaseVariableName = rInputName.substr(0, rInputName.length() - 2);
        rComponentIndex = 1;
    }
    else if ( StringEndsWith(rInputName,z_suffix) )
    {
        is_component = true;
        rBaseVariableName = rInputName.substr(0, rInputName.length() - 2);
        rComponentIndex = 2;
    }

    return is_component;
}

bool StringEndsWith(std::string const &rString, std::string const &rEnding) const {
    if (rString.length() >= rEnding.length()) {
        return (0 == rString.compare(rString.length() - rEnding.length(), rEnding.length(), rEnding));
    } else {
        return false;
    }
}

///@}
///@name Private  Access
///@{

///@}
///@name Private Inquiry
///@{

///@}
///@name Un accessible methods
///@{

/// Assignment operator.
IntegrationPointStatisticsProcess &operator=(IntegrationPointStatisticsProcess const &rOther) = delete;

/// Copy constructor.
IntegrationPointStatisticsProcess(IntegrationPointStatisticsProcess const &rOther) = delete;

///@}

}; // Class IntegrationPointStatisticsProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                IntegrationPointStatisticsProcess &rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const IntegrationPointStatisticsProcess &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_INTEGRATION_POINT_STATISTICS_PROCESS_H_INCLUDED  defined
