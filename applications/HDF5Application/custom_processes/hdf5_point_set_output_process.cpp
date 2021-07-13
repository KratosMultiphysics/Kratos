//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:     Máté Kelemen
//

// Internal includes
#include "hdf5_point_set_output_process.h"

// Core includes
#include "utilities/parallel_utilities.h"


namespace Kratos
{
namespace HDF5
{


PointSetOutputProcess::PointSetOutputProcess(Parameters parameters, Model& rModel)
    : mrModelPart(rModel.GetModelPart(parameters["model_part_name"].GetString())),
      mInterval(parameters)
{
    KRATOS_TRY

    parameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
    this->InitializeFromParameters(parameters);

    // Populate points
    const bool is_historical = parameters["historical_value"].GetBool();
    const std::size_t number_of_points = parameters["positions"].size();
    mVertices.reserve(number_of_points);

    for (std::size_t i_point=0; i_point<number_of_points; ++i_point) {
        mVertices.push_back(std::make_shared<Detail::Vertex>(
            parameters["positions"].GetArrayItem(i_point).GetVector(),
            mrModelPart,
            is_historical));
    }

    KRATOS_CATCH("");
}


PointSetOutputProcess::PointSetOutputProcess(Parameters parameters,
                                             Model& rModel,
                                             Detail::VertexContainerType&& rVertices)
    : mrModelPart(rModel.GetModelPart(parameters["model_part_name"].GetString())),
      mInterval(parameters),
      mVertices(std::move(rVertices))
{
    KRATOS_TRY

    parameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
    this->InitializeFromParameters(parameters);

    KRATOS_CATCH("");
}


void PointSetOutputProcess::ExecuteFinalizeSolutionStep()
{
    KRATOS_TRY

    // Update vertices
    block_for_each(mVertices,
        [this](Detail::Vertex& rVertex){
            rVertex.Locate(*mpLocator);
        }
    );

    KRATOS_CATCH("");
}


const Parameters PointSetOutputProcess::GetDefaultParameters() const
{
    return Parameters(R"({
        "model_part_name"      : "",
        "enitity_type"         : "element",
        "interval"             : [0.0, 1e20],
        "positions"            : [[]],
        "output_variables"     : [],
        "historical_value"     : true,
        "search_configuration" : "initial",
        "search_tolerance"     : 1e-6,
        "output_path"          : "",
        "prefix"               : "point_set_output"
    })");
}


PointSetOutputProcess::Locator::Locator(ModelPart& rModelPart, const Globals::Configuration configuration, const double tolerance)
    : mLocator(rModelPart),
      mConfiguration(configuration),
      mTolerance(tolerance)
{
}


int PointSetOutputProcess::Locator::FindElement(const Point& rPoint, Kratos::Vector& rShapeFunctionValues) const
{
    KRATOS_TRY

    return mLocator.FindElement(
        rPoint,
        rShapeFunctionValues,
        mConfiguration,
        mTolerance);

    KRATOS_CATCH("");
}


void PointSetOutputProcess::InitializeFromParameters(Parameters parameters)
{
    KRATOS_TRY

    // Initialize variable names
    const std::size_t number_of_variables = parameters["output_variables"].size();
    mVariableNames.reserve(number_of_variables);

    for (std::size_t i_name; i_name<number_of_variables; ++i_name) {
        mVariableNames.push_back(parameters["output_variables"].GetArrayItem(i_name).GetString());
    }

    // Initialize locator
    const std::string configuration_name = parameters["search_configuration"].GetString();
    const double search_tolerance = parameters["search_tolerance"].GetDouble();

    Globals::Configuration configuration;
    if (configuration_name == "initial") {
        configuration = Globals::Configuration::Initial;
    }
    else if (configuration_name == "current") {
        configuration = Globals::Configuration::Current;
    }
    else {
        KRATOS_ERROR << "Invalid configuration '" << configuration_name << "', options are: 'initial', 'current'";
    }

    mpLocator = std::unique_ptr<PointSetOutputProcess::Locator>(new PointSetOutputProcess::Locator(
        mrModelPart,
        configuration,
        search_tolerance));

    // Initialize other members
    mOutputPath = parameters["output_path"].GetString();

    KRATOS_CATCH("");
}


File::Pointer PointSetOutputProcess::MakeHDF5File(const std::size_t stepID) const
{
    KRATOS_TRY

    Parameters file_parameters(R"({
        "file_name" : ""
    })");

    const std::string file_name = mOutputPath + "/step_" + std::to_string(stepID);
    file_parameters["file_name"].SetString(file_name);

    return make_shared<File>(file_parameters);

    KRATOS_CATCH("");
}


} // namespace HDF5
} // namespace Kratos