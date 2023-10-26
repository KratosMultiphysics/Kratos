// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//

#include "write_output.h"
#include "includes/model_part.h"
#include "custom_utilities/element_utilities.hpp"


namespace
{

using namespace Kratos;

using NodalResultWriter = std::function<void(GidIO<>&, const ModelPart&)>;

template <typename VariableType>
NodalResultWriter MakeNodalResultWriterFor(const VariableType& rVariable)
{
    return [&rVariable](GidIO<>& rGidIO, const ModelPart& rModelPart)
    {
        const auto     Time  = rModelPart.GetProcessInfo()[TIME];
        constexpr auto Index = 0;
        rGidIO.WriteNodalResults(rVariable, rModelPart.Nodes(), Time, Index);
    };
}

using IntegrationPointResultWriter = std::function<void(GidIO<>&, const ModelPart&)>;

template <typename VariableType>
IntegrationPointResultWriter MakeIntegrationPointResultWriterFor(const VariableType& rVariable)
{
    return [&rVariable](GidIO<>& rGidIO, const ModelPart& rModelPart)
    {
        const auto     Time  = rModelPart.GetProcessInfo()[TIME];
        constexpr auto Index = 0;
        rGidIO.PrintOnGaussPoints(rVariable, rModelPart, Time, Index);
    };
}

}


namespace Kratos
{

void GeoOutputWriter::WriteGiDOutput(ModelPart&         rModelPart,
                                     Parameters         Settings,
                                     const std::string& rWorkingDirectory,
                                     bool WriteHydraulicHeadToNodes)
{
    auto output_parameters = Settings["output_processes"]["gid_output"].GetArrayItem(0)["Parameters"];
    auto gid_post_flags = output_parameters["postprocess_parameters"]["result_file_configuration"]["gidpost_flags"];

    // Calculate hydraulic head on the nodes
    auto gauss_outputs = output_parameters["postprocess_parameters"]["result_file_configuration"]["gauss_point_results"].GetStringArray();
    if (WriteHydraulicHeadToNodes && std::find(gauss_outputs.begin(), gauss_outputs.end(), "HYDRAULIC_HEAD") != gauss_outputs.end())
    {
        CalculateNodalHydraulicHead(mGidIO, rModelPart);
    }

    const auto nodal_outputs = output_parameters["postprocess_parameters"]["result_file_configuration"]["nodal_results"].GetStringArray();
    WriteNodalOutput(nodal_outputs, mGidIO, rModelPart);
    WriteIntegrationPointOutput(gauss_outputs, mGidIO, rModelPart);
}


void GeoOutputWriter::WriteNodalOutput(const std::vector<std::string>& rOutputItemNames,
                                       GidIO<>&                        rGidIO,
                                       const ModelPart&                rModelPart)
{
    const auto output_writer_map = std::map<std::string, NodalResultWriter, std::less<>>{
            {"DISPLACEMENT",        MakeNodalResultWriterFor(DISPLACEMENT)},
            {"TOTAL_DISPLACEMENT",  MakeNodalResultWriterFor(TOTAL_DISPLACEMENT)},
            {"WATER_PRESSURE",      MakeNodalResultWriterFor(WATER_PRESSURE)},
            {"NORMAL_FLUID_FLUX",   MakeNodalResultWriterFor(NORMAL_FLUID_FLUX)},
            {"VOLUME_ACCELERATION", MakeNodalResultWriterFor(VOLUME_ACCELERATION)},
            {"HYDRAULIC_DISCHARGE", MakeNodalResultWriterFor(HYDRAULIC_DISCHARGE)},
            {"HYDRAULIC_HEAD",      MakeNodalResultWriterFor(HYDRAULIC_HEAD)}
    };

    for (const auto& name : rOutputItemNames)
    {
        output_writer_map.at(name)(rGidIO, rModelPart);
    }
}


void GeoOutputWriter::WriteIntegrationPointOutput(const std::vector<std::string>& rOutputItemNames,
                                                  GidIO<>&                        rGidIO,
                                                  const ModelPart&                rModelPart)
{
    const auto output_writer_map = std::map<std::string, IntegrationPointResultWriter, std::less<>>{
            {"FLUID_FLUX_VECTOR",            MakeIntegrationPointResultWriterFor(FLUID_FLUX_VECTOR)},
            {"HYDRAULIC_HEAD",               MakeIntegrationPointResultWriterFor(HYDRAULIC_HEAD)},
            {"LOCAL_FLUID_FLUX_VECTOR",      MakeIntegrationPointResultWriterFor(LOCAL_FLUID_FLUX_VECTOR)},
            {"LOCAL_PERMEABILITY_MATRIX",    MakeIntegrationPointResultWriterFor(LOCAL_PERMEABILITY_MATRIX)},
            {"PERMEABILITY_MATRIX",          MakeIntegrationPointResultWriterFor(PERMEABILITY_MATRIX)},
            {"DEGREE_OF_SATURATION",         MakeIntegrationPointResultWriterFor(DEGREE_OF_SATURATION)},
            {"DERIVATIVE_OF_SATURATION",     MakeIntegrationPointResultWriterFor(DERIVATIVE_OF_SATURATION)},
            {"RELATIVE_PERMEABILITY",        MakeIntegrationPointResultWriterFor(RELATIVE_PERMEABILITY)},
            {"PIPE_ACTIVE",                  MakeIntegrationPointResultWriterFor(PIPE_ACTIVE)},
            {"PIPE_HEIGHT",                  MakeIntegrationPointResultWriterFor(PIPE_HEIGHT)},
            {"GREEN_LAGRANGE_STRAIN_TENSOR", MakeIntegrationPointResultWriterFor(GREEN_LAGRANGE_STRAIN_TENSOR)},
            {"ENGINEERING_STRAIN_TENSOR",    MakeIntegrationPointResultWriterFor(ENGINEERING_STRAIN_TENSOR)},
            {"CAUCHY_STRESS_TENSOR",         MakeIntegrationPointResultWriterFor(CAUCHY_STRESS_TENSOR)},
            {"TOTAL_STRESS_TENSOR",          MakeIntegrationPointResultWriterFor(TOTAL_STRESS_TENSOR)},
            {"VON_MISES_STRESS",             MakeIntegrationPointResultWriterFor(VON_MISES_STRESS)}
    };

    for (const auto& name : rOutputItemNames)
    {
        output_writer_map.at(name)(rGidIO, rModelPart);
    }
}

void GeoOutputWriter::CalculateNodalHydraulicHead(GidIO<>& rGidIO, ModelPart& rModelPart)
{
    const auto& element_var = KratosComponents<Variable<double>>::Get("HYDRAULIC_HEAD");

    for (Element element : rModelPart.Elements())
    {
        auto& rGeom = element.GetGeometry();
        const auto& rProp = element.GetProperties();
        const auto NodalHydraulicHead = GeoElementUtilities::CalculateNodalHydraulicHeadFromWaterPressures(rGeom, rProp);

        for (unsigned int node = 0; node < rGeom.PointsNumber(); ++node)
        {
            rGeom[node].SetValue(element_var, NodalHydraulicHead[node]);
        }
    }
    rGidIO.WriteNodalResultsNonHistorical(element_var, rModelPart.Nodes(), 0);
}

GiD_PostMode GeoOutputWriter::GetGiDPostModeFrom(const Parameters& rGiDPostFlags)
{
    const std::map<std::string, GiD_PostMode, std::less<>> to_post_mode{
        {"GiD_PostAscii",       GiD_PostAscii},
        {"GiD_PostAsciiZipped", GiD_PostAsciiZipped},
        {"GiD_PostBinary",      GiD_PostBinary},
        {"GiD_PostHDF5",        GiD_PostHDF5}
    };
    return to_post_mode.at(rGiDPostFlags["GiDPostMode"].GetString());
}

MultiFileFlag GeoOutputWriter::GetMultiFileFlagFrom(const Parameters& rGiDPostFlags)
{
    const std::map<std::string, MultiFileFlag, std::less<>> to_multi_file_flag{
        {"SingleFile",    SingleFile},
        {"MultipleFiles", MultipleFiles}
    };
    return to_multi_file_flag.at(rGiDPostFlags["MultiFileFlag"].GetString());
}

WriteDeformedMeshFlag GeoOutputWriter::GetWriteDeformedMeshFlagFrom(const Parameters& rGiDPostFlags)
{
    const std::map<std::string, WriteDeformedMeshFlag, std::less<>> to_write_deformed_flag{
        {"WriteDeformed",   WriteDeformed},
        {"WriteUndeformed", WriteUndeformed}
    };
    return to_write_deformed_flag.at(rGiDPostFlags["WriteDeformedMeshFlag"].GetString());
}

WriteConditionsFlag GeoOutputWriter::GetWriteConditionsFlagFrom(const Parameters& rGiDPostFlags)
{
    const std::map<std::string, WriteConditionsFlag, std::less<>> to_write_conditions_flag{
        {"WriteConditions",     WriteConditions},
        {"WriteElementsOnly",   WriteElementsOnly},
        {"WriteConditionsOnly", WriteConditionsOnly}
    };
    return to_write_conditions_flag.at(rGiDPostFlags["WriteConditionsFlag"].GetString());
}

GeoOutputWriter::GeoOutputWriter(                                     Parameters         Settings,
                                 const std::string& rWorkingDirectory,
                                 ModelPart&         rModelPart) :
      mGidIO{rWorkingDirectory + "/" + Settings["output_processes"]["gid_output"].GetArrayItem(0)["Parameters"]["output_name"].GetString(),
        GetGiDPostModeFrom(Settings["output_processes"]["gid_output"].GetArrayItem(0)["Parameters"]["postprocess_parameters"]["result_file_configuration"]["gidpost_flags"]),
        GetMultiFileFlagFrom(Settings["output_processes"]["gid_output"].GetArrayItem(0)["Parameters"]["postprocess_parameters"]["result_file_configuration"]["gidpost_flags"]),
        GetWriteDeformedMeshFlagFrom(Settings["output_processes"]["gid_output"].GetArrayItem(0)["Parameters"]["postprocess_parameters"]["result_file_configuration"]["gidpost_flags"]),
        GetWriteConditionsFlagFrom(Settings["output_processes"]["gid_output"].GetArrayItem(0)["Parameters"]["postprocess_parameters"]["result_file_configuration"]["gidpost_flags"])}
{
    mGidIO.InitializeMesh(0.0);
    mGidIO.WriteMesh(rModelPart.GetMesh());
    mGidIO.FinalizeMesh();
    mGidIO.InitializeResults(0, rModelPart.GetMesh());
}

void GeoOutputWriter::FinalizeResults()
{
    mGidIO.FinalizeResults();
}

}
