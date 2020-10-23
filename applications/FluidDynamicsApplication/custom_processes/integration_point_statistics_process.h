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
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// Helper process to record statistics on the integration points of the mesh.
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

/// Constructor using a ModelPart.
/** @param rModelPart ModelPart the statistics will be calculated on.
 *  @param pParameters Parameters object defining the desired statistics.
 */
IntegrationPointStatisticsProcess(ModelPart& rModelPart, Kratos::Parameters::Pointer pParameters):
    Process(),
    mrModelPart(rModelPart),
    mpParameters(pParameters),
    mDimension(3),
    mStartTime(0.0)
{}

/// Constructor using a Model.
/** @param rModel Model instance containing the ModelPart the statistics will be calculated on.
 *  @param pParameters Parameters object defining the desired statistics.
 */
IntegrationPointStatisticsProcess(Model& rModel, Kratos::Parameters::Pointer pParameters):
    Process(),
    mrModelPart(rModel.GetModelPart((*pParameters)["model_part_name"].GetString())),
    mpParameters(pParameters),
    mDimension(3),
    mStartTime(0.0)
{}

/// Destructor.
~IntegrationPointStatisticsProcess() override
{}

///@}
///@name Operations
///@{

/// Define the statistics and initialize internal storage.
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

protected:

///@name Protected Operations
///@{

void CreateStatisticsFromInput(StatisticsRecord::Pointer pRecordedStatistics)
{
    KRATOS_TRY;

    // Validate parameters
    Kratos::Parameters default_parameters = Kratos::Parameters(R"({
        "statistics" : [],
        "output_file_name": "statistics",
        "model_part_name": "",
        "dimension": 3,
        "start_time": 0.0
    })");

    Parameters& r_parameters = *mpParameters;
    r_parameters.ValidateAndAssignDefaults(default_parameters);

    mOutputFileName = r_parameters["output_file_name"].GetString();
    mStartTime = r_parameters["start_time"].GetDouble();

    mDimension = r_parameters["dimension"].GetInt();
    KRATOS_ERROR_IF( mDimension != 2 && mDimension != 3)
    << "Unsupported dimension " << mDimension << ". IntegrationPointStatisticsProcess only works for dimension 2 or 3." << std::endl;

    StatisticsDictionary defined_statistics;

    for (unsigned int i = 0; i < r_parameters["statistics"].size(); i++)
    {
        Kratos::Parameters settings = r_parameters["statistics"][i];
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
            pRecordedStatistics->AddHigherOrderStatistic(CreateThirdOrderSampler(settings,defined_statistics));
        }
        else
        {
            KRATOS_ERROR << "Unknown string \"" << statistic_type << "\" passed as \"type\" argument in statistics definition." << std::endl;
        }
    }

    KRATOS_CATCH("");
}

StatisticsSampler::Pointer CreateAverageSampler(
    Kratos::Parameters& rParameters,
    StatisticsDictionary& rDefinedStatistics) const
{
    KRATOS_TRY;

    Kratos::Parameters default_parameters(R"({
        "type" : "average",
        "variable": "",
        "tags": []
    })");

    rParameters.ValidateAndAssignDefaults(default_parameters);
    std::string variable_name = rParameters["variable"].GetString();
    std::string type = rParameters["type"].GetString();
    KRATOS_ERROR_IF_NOT(type == "average") << "Trying to define an average statistic of unsupported type " << type << "." << std::endl;

    KRATOS_ERROR_IF(rDefinedStatistics.find(variable_name) != rDefinedStatistics.end())
    << "Duplicate definition of an average for " << variable_name << "." << std::endl;

    StatisticsSampler::Pointer new_statistic;
    if (KratosComponents<Variable<double>>::Has(variable_name))
    {
        // build double variable sampler
        const Variable<double>& variable = KratosComponents<Variable<double>>::Get(variable_name);
        auto value_getter = Kratos::Internals::MakeSamplerAtLocalCoordinate::ValueGetter(variable);
        std::string tag = variable_name;
        if (rParameters["tags"].size() > 0)
        {
            KRATOS_ERROR_IF(rParameters["tags"].size() != 1) << "Only one tag is needed for scalar averages, but "
            << rParameters["tags"].size() << " were provided for variable " << variable_name << "." << std::endl;

            tag = rParameters["tags"][0].GetString();
        }
        new_statistic = Kratos::make_shared<ScalarAverageSampler>(value_getter,tag);
    }
    else if (KratosComponents<Variable<array_1d<double,3>>>::Has(variable_name))
    {
        // build vector sampler
        const Variable<array_1d<double,3>>& variable = KratosComponents<Variable<array_1d<double,3>>>::Get(variable_name);
        auto value_getter = Kratos::Internals::MakeSamplerAtLocalCoordinate::ValueGetter(variable);
        std::vector<std::string> tags;
        if (rParameters["tags"].size() > 0)
        {
            KRATOS_ERROR_IF(rParameters["tags"].size() != mDimension) << mDimension << " tags are needed for vector averages, but "
            << rParameters["tags"].size() << " were provided for variable " << variable_name << "." << std::endl;

            for (unsigned int i = 0; i < mDimension; i++) {
                tags.push_back(rParameters["tags"][i].GetString());
            }
        }
        else
        {
            tags.push_back(std::string(variable_name+"_X"));
            tags.push_back(std::string(variable_name+"_Y"));
            if (mDimension == 3) tags.push_back(std::string(variable_name+"_Z"));
        }
        new_statistic = Kratos::make_shared<VectorAverageSampler<array_1d<double,3>>>(value_getter,mDimension,tags);
    }
    else
    {
        KRATOS_ERROR
        << "Trying to define an average statistic for variable " << variable_name
        << " which is not a variable of a supported type." << std::endl;
    }

    rDefinedStatistics[variable_name] = new_statistic;
    return new_statistic;

    KRATOS_CATCH("");
}

StatisticsSampler::Pointer CreateVarianceSampler(
    Kratos::Parameters& rParameters,
    StatisticsDictionary& rDefinedStatistics) const
{
    KRATOS_TRY;

    Kratos::Parameters default_parameters(R"({
        "type" : "",
        "variables": []
    })");

    rParameters.ValidateAndAssignDefaults(default_parameters);
    KRATOS_ERROR_IF(rParameters["variables"].size() < 1 || rParameters["variables"].size() > 2)
    << "Unexpected number of arguments when reading \"variables\" list argument. "
    << "Expected 1 or 2 values, got " << rParameters["variables"].size() << "." << std::endl;

    StatisticsSampler::Pointer new_statistic;
    if (rParameters["variables"].size() == 1)
    {
        // symmetric variance
        std::string variable_name = rParameters["variables"][0].GetString();
        StatisticsDictionary::iterator it_average = rDefinedStatistics.find(variable_name);
        KRATOS_ERROR_IF(it_average == rDefinedStatistics.end())
        << "Trying to define a variance for " << variable_name
        << " but no average has been defined for this variable." << std::endl;

        new_statistic = Kratos::make_shared<SymmetricVarianceSampler>(it_average->second);
        rDefinedStatistics[MakeVarianceKey(variable_name,variable_name)] = new_statistic;
    }
    else // size == 2
    {
        // complete or componentwise variance
        std::string first_argument = rParameters["variables"][0].GetString();
        std::string first_argument_base;
        unsigned int first_argument_index;
        bool first_argument_is_component = ProcessComponent(first_argument,first_argument_base,first_argument_index);

        std::string second_argument = rParameters["variables"][1].GetString();
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

        if (
             (first_argument_is_component && !(second_argument_is_component || IsScalar(second_argument) )) ||
             (second_argument_is_component && !(first_argument_is_component || IsScalar(first_argument) ))
        )
        {
            KRATOS_ERROR << "Variances involving a component and a vector are currently not supported, please define them component-by-component." << std::endl;
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

        rDefinedStatistics[MakeVarianceKey(first_argument,second_argument)] = new_statistic;
    }
    return new_statistic;

    KRATOS_CATCH("");
}

StatisticsSampler::Pointer CreateThirdOrderSampler(
    Kratos::Parameters& rParameters,
    StatisticsDictionary& rDefinedStatistics) const
{
    KRATOS_TRY;

    Kratos::Parameters default_parameters(R"({
        "type" : "",
        "variables": []
    })");

    rParameters.ValidateAndAssignDefaults(default_parameters);
    KRATOS_ERROR_IF(rParameters["variables"].size() != 3 )
    << "Unexpected number of arguments when reading \"variables\" list argument. "
    << "Expected 3 values, got " << rParameters["variables"].size() << "." << std::endl;

    const std::string first_argument  = rParameters["variables"][0].GetString();
    const std::string second_argument = rParameters["variables"][1].GetString();
    const std::string third_argument  = rParameters["variables"][2].GetString();

    std::string first_argument_base;
    unsigned int first_argument_index;
    ProcessThirdOrderInputValue(first_argument, first_argument_base, first_argument_index);
    auto found = rDefinedStatistics.find(first_argument_base);
    KRATOS_ERROR_IF(found == rDefinedStatistics.end())
    << "Trying to record third order moment for " << first_argument << ", " << second_argument << " and " << third_argument
    << " but no average was declared for " << first_argument_base << "." << std::endl;
    StatisticsSampler::Pointer first_average = found->second;

    std::string second_argument_base;
    unsigned int second_argument_index;
    ProcessThirdOrderInputValue(second_argument, second_argument_base, second_argument_index);
    found = rDefinedStatistics.find(second_argument_base);
    KRATOS_ERROR_IF(found == rDefinedStatistics.end())
    << "Trying to record third order moment for " << first_argument << ", " << second_argument << " and " << third_argument
    << " but no average was declared for " << second_argument_base << "." << std::endl;
    StatisticsSampler::Pointer second_average = found->second;

    std::string third_argument_base;
    unsigned int third_argument_index;
    ProcessThirdOrderInputValue(third_argument, third_argument_base, third_argument_index);
    found = rDefinedStatistics.find(third_argument_base);
    KRATOS_ERROR_IF(found == rDefinedStatistics.end())
    << "Trying to record third order moment for " << first_argument << ", " << second_argument << " and " << third_argument
    << " but no average was declared for " << third_argument_base << "." << std::endl;
    StatisticsSampler::Pointer third_average = found->second;

    found = rDefinedStatistics.find(MakeVarianceKey(first_argument_base,second_argument_base));
    KRATOS_ERROR_IF(found == rDefinedStatistics.end())
    << "Trying to record third order moment for " << first_argument << ", " << second_argument << " and " << third_argument
    << " but no correlation was declared for " << first_argument_base << " and " << second_argument_base << "." << std::endl;
    StatisticsSampler::Pointer first_second_correlation = found->second;

    found = rDefinedStatistics.find(MakeVarianceKey(first_argument_base,third_argument_base));
    KRATOS_ERROR_IF(found == rDefinedStatistics.end())
    << "Trying to record third order moment for " << first_argument << ", " << second_argument << " and " << third_argument
    << " but no correlation was declared for " << first_argument_base << " and " << third_argument_base << "." << std::endl;
    StatisticsSampler::Pointer first_third_correlation = found->second;


    found = rDefinedStatistics.find(MakeVarianceKey(second_argument_base,third_argument_base));
    KRATOS_ERROR_IF(found == rDefinedStatistics.end())
    << "Trying to record third order moment for " << first_argument << ", " << second_argument << " and " << third_argument
    << " but no correlation was declared for " << second_argument_base << " and " << third_argument_base << "." << std::endl;
    StatisticsSampler::Pointer second_third_correlation = found->second;

    return Kratos::make_shared<ThirdOrderCorrelationSampler>(
        first_average, first_argument_index,
        second_average, second_argument_index,
        third_average, third_argument_index,
        first_second_correlation, first_second_correlation->ComponentIndex(first_argument_index,second_argument_index),
        first_third_correlation, first_third_correlation->ComponentIndex(first_argument_index,third_argument_index),
        second_third_correlation, second_third_correlation->ComponentIndex(second_argument_index, third_argument_index)
    );

    KRATOS_CATCH("");
}

///@}

private:

///@name Member Variables
///@{

ModelPart& mrModelPart;

Kratos::Parameters::Pointer mpParameters;

unsigned int mDimension;

double mStartTime;

std::string mOutputFileName;

///@}
///@name Private Operations
///@{

/// Helper function to check if the requested variable is a scalar variable.
bool IsScalar(const std::string& rInputName) const
{
    return KratosComponents< Variable<double> >::Has(rInputName);
}

/// Helper function to determine if the requested variable is a vector component (and if so, which one).
bool ProcessComponent(
    const std::string& rInputName,
    std::string& rBaseVariableName,
    unsigned int& rComponentIndex) const
{
    KRATOS_TRY

    bool is_component = false;

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
        KRATOS_ERROR_IF_NOT(mDimension == 3)
        << "Trying to record a statistic for " << rInputName
        << " but IntegrationPointStatisticsProcess was initialized with dimension 2." << std::endl;
    }
    else
    {
        rComponentIndex = 0; // can be used also if variable is not a component
        rBaseVariableName = rInputName;
    }

    return is_component;

    KRATOS_CATCH("");
}

/// Helper function to check whether a string ends with a given substring.
bool StringEndsWith(std::string const &rString, std::string const &rEnding) const {
    if (rString.length() >= rEnding.length()) {
        return (0 == rString.compare(rString.length() - rEnding.length(), rEnding.length(), rEnding));
    } else {
        return false;
    }
}

/// Helper function to ensure that variance keys in the dictionary of defined statistics are independent of the order of arguments.
std::string MakeVarianceKey(std::string const &rFirstArgument, std::string const &rSecondArgument) const
{
    if (rFirstArgument < rSecondArgument)
    {
        return std::string("Variance_" + rFirstArgument + "_" + rSecondArgument);
    }
    else
    {
        return std::string("Variance_" + rSecondArgument + "_" + rFirstArgument);
    }
}

/// Helper function to read input for third order correlations
void ProcessThirdOrderInputValue(
    const std::string& rArgument,
    std::string& rBaseName,
    unsigned int& rComponentIndex ) const
{
    bool argument_is_scalar = IsScalar(rArgument);
    bool argument_is_component = ProcessComponent(rArgument,rBaseName,rComponentIndex);
    KRATOS_ERROR_IF_NOT(argument_is_scalar || argument_is_component)
    << "Unexpected variable " << rArgument << " found while defining a third order correlation. "
    << "Only scalar or component variables are supported at the moment." << std::endl;
}

///@name Un accessible methods
///@{

/// Assignment operator.
IntegrationPointStatisticsProcess &operator=(IntegrationPointStatisticsProcess const &rOther) = delete;

/// Copy constructor.
IntegrationPointStatisticsProcess(IntegrationPointStatisticsProcess const &rOther) = delete;

///@}

}; // Class IntegrationPointStatisticsProcess

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
