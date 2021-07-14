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

// Application includes
#include "custom_io/hdf5_vertex_container_io.h"
#include "custom_io/hdf5_file_parallel.h"
#include "custom_io/hdf5_file_serial.h"

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


void PointSetOutputProcess::ExecuteInitialize()
{
    KRATOS_TRY

    // Dump vertex coordinates
    Parameters io_parameters;
    io_parameters.AddString("prefix", mPrefix);
    VertexContainerIO(io_parameters, mpFile).WriteCoordinates(mVertices);

    KRATOS_CATCH("");
}


void PointSetOutputProcess::ExecuteFinalizeSolutionStep()
{
    KRATOS_TRY

    // Create IO
    const std::size_t step_id = mrModelPart.GetProcessInfo()[STEP];
    const std::string group_name = "step_" + std::to_string(step_id);
    Parameters io_parameters;
    io_parameters.AddString("prefix", mPrefix);
    io_parameters.AddString("variables_path", "/" + group_name);
    io_parameters.AddValue("list_of_variables", mVariables);

    VertexContainerIO vertex_io(
        io_parameters,
        mpFile);

    // Update vertices
    block_for_each(mVertices,
        [this](Detail::Vertex& rVertex){
            rVertex.Locate(*mpLocator);
        }
    );

    // Write
    vertex_io.WriteVariables(mVertices);

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
        "file_path"            : "",
        "prefix"               : "/point_set_output"
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
    mVariables = parameters["output_variables"].Clone();

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

    // Initialize file
    Parameters file_parameters;
    file_parameters.AddValue("file_name", parameters["file_path"]);

    const bool is_mpi_run = false; // TODO: detect serial/parallel
    if (is_mpi_run) {
        mpFile = std::make_shared<FileParallel>(file_parameters);
    }
    else {
        mpFile = std::make_shared<FileSerial>(file_parameters);
    }

    // Initialize other members
    mPrefix = parameters["prefix"].GetString();

    KRATOS_CATCH("");
}


} // namespace HDF5
} // namespace Kratos