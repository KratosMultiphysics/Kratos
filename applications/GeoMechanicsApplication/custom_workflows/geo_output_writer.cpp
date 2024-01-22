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

#include "geo_output_writer.h"
#include "custom_utilities/element_utilities.hpp"
#include "includes/model_part.h"

namespace
{

using namespace Kratos;

using NodalResultWriter = std::function<void(GidIO<>&, const ModelPart&)>;

template <typename VariableType>
NodalResultWriter MakeNodalResultWriterFor(const VariableType& rVariable)
{
    return [&rVariable](GidIO<>& rGidIO, const ModelPart& rModelPart) {
        const auto     Time  = rModelPart.GetProcessInfo()[TIME];
        constexpr auto Index = 0;
        rGidIO.WriteNodalResults(rVariable, rModelPart.Nodes(), Time, Index);
    };
}

using IntegrationPointResultWriter = std::function<void(GidIO<>&, const ModelPart&)>;

template <typename VariableType>
IntegrationPointResultWriter MakeIntegrationPointResultWriterFor(const VariableType& rVariable)
{
    return [&rVariable](GidIO<>& rGidIO, const ModelPart& rModelPart) {
        const auto     Time  = rModelPart.GetProcessInfo()[TIME];
        constexpr auto Index = 0;
        rGidIO.PrintOnGaussPoints(rVariable, rModelPart, Time, Index);
    };
}

GiD_PostMode GetGiDPostModeFrom(const Parameters& rGiDPostFlags)
{
    const auto to_post_mode =
        std::map<std::string, GiD_PostMode, std::less<>>{{"GiD_PostAscii", GiD_PostAscii},
                                                         {"GiD_PostAsciiZipped", GiD_PostAsciiZipped},
                                                         {"GiD_PostBinary", GiD_PostBinary},
                                                         {"GiD_PostHDF5", GiD_PostHDF5}};
    return to_post_mode.at(rGiDPostFlags["GiDPostMode"].GetString());
}

MultiFileFlag GetMultiFileFlagFrom(const Parameters& rGiDPostFlags)
{
    const auto to_multi_file_flag = std::map<std::string, MultiFileFlag, std::less<>>{
        {"SingleFile", SingleFile}, {"MultipleFiles", MultipleFiles}};
    return to_multi_file_flag.at(rGiDPostFlags["MultiFileFlag"].GetString());
}

WriteDeformedMeshFlag GetWriteDeformedMeshFlagFrom(const Parameters& rGiDPostFlags)
{
    const auto to_write_deformed_flag = std::map<std::string, WriteDeformedMeshFlag, std::less<>>{
        {"WriteDeformed", WriteDeformed}, {"WriteUndeformed", WriteUndeformed}};
    return to_write_deformed_flag.at(rGiDPostFlags["WriteDeformedMeshFlag"].GetString());
}

WriteConditionsFlag GetWriteConditionsFlagFrom(const Parameters& rGiDPostFlags)
{
    const auto to_write_conditions_flag = std::map<std::string, WriteConditionsFlag, std::less<>>{
        {"WriteConditions", WriteConditions},
        {"WriteElementsOnly", WriteElementsOnly},
        {"WriteConditionsOnly", WriteConditionsOnly}};
    return to_write_conditions_flag.at(rGiDPostFlags["WriteConditionsFlag"].GetString());
}

} // namespace

namespace Kratos
{

GeoOutputWriter::GeoOutputWriter(const Parameters&  rGidOutputSettings,
                                 const std::string& rWorkingDirectory,
                                 ModelPart&         rModelPart)
    : mGidIO{MakeGidIO(rWorkingDirectory, rGidOutputSettings)}
{
    mGidIO.InitializeMesh(0.0);
    mGidIO.WriteMesh(rModelPart.GetMesh());
    mGidIO.FinalizeMesh();
    mGidIO.InitializeResults(0, rModelPart.GetMesh());
}

void GeoOutputWriter::WriteGiDOutput(ModelPart& rModelPart, Parameters Settings, bool WriteHydraulicHeadToNodes)
{
    // Calculate hydraulic head on the nodes
    const auto gauss_outputs =
        Settings["postprocess_parameters"]["result_file_configuration"]["gauss_point_results"].GetStringArray();
    if (WriteHydraulicHeadToNodes &&
        std::find(gauss_outputs.begin(), gauss_outputs.end(), "HYDRAULIC_HEAD") != gauss_outputs.end()) {
        CalculateNodalHydraulicHead(rModelPart);
    }

    const auto nodal_outputs =
        Settings["postprocess_parameters"]["result_file_configuration"]["nodal_results"].GetStringArray();
    WriteNodalOutput(nodal_outputs, rModelPart);
    WriteIntegrationPointOutput(gauss_outputs, rModelPart);
}

void GeoOutputWriter::WriteNodalOutput(const std::vector<std::string>& rOutputItemNames, const ModelPart& rModelPart)
{
    const auto output_writer_map = std::map<std::string, NodalResultWriter, std::less<>>{
        {"DISPLACEMENT", MakeNodalResultWriterFor(DISPLACEMENT)},
        {"TOTAL_DISPLACEMENT", MakeNodalResultWriterFor(TOTAL_DISPLACEMENT)},
        {"WATER_PRESSURE", MakeNodalResultWriterFor(WATER_PRESSURE)},
        {"NORMAL_FLUID_FLUX", MakeNodalResultWriterFor(NORMAL_FLUID_FLUX)},
        {"VOLUME_ACCELERATION", MakeNodalResultWriterFor(VOLUME_ACCELERATION)},
        {"HYDRAULIC_DISCHARGE", MakeNodalResultWriterFor(HYDRAULIC_DISCHARGE)},
        {"HYDRAULIC_HEAD", MakeNodalResultWriterFor(HYDRAULIC_HEAD)},
        {"LINE_LOAD_Y", MakeNodalResultWriterFor(LINE_LOAD_Y)},
        {"REACTION", MakeNodalResultWriterFor(REACTION)},
        {"NORMAL_CONTACT_STRESS", MakeNodalResultWriterFor(NORMAL_CONTACT_STRESS)},
        {"TANGENTIAL_CONTACT_STRESS", MakeNodalResultWriterFor(TANGENTIAL_CONTACT_STRESS)}};

    for (const auto& name : rOutputItemNames) {
        auto iter = output_writer_map.find(name);
        KRATOS_ERROR_IF(iter == output_writer_map.end())
            << "Output item '" << name << "' is not available for nodal output" << std::endl;

        iter->second(mGidIO, rModelPart);
    }
}

void GeoOutputWriter::WriteIntegrationPointOutput(const std::vector<std::string>& rOutputItemNames,
                                                  const ModelPart&                rModelPart)
{
    const auto output_writer_map = std::map<std::string, IntegrationPointResultWriter, std::less<>>{
        {"FLUID_FLUX_VECTOR", MakeIntegrationPointResultWriterFor(FLUID_FLUX_VECTOR)},
        {"HYDRAULIC_HEAD", MakeIntegrationPointResultWriterFor(HYDRAULIC_HEAD)},
        {"LOCAL_FLUID_FLUX_VECTOR", MakeIntegrationPointResultWriterFor(LOCAL_FLUID_FLUX_VECTOR)},
        {"LOCAL_PERMEABILITY_MATRIX", MakeIntegrationPointResultWriterFor(LOCAL_PERMEABILITY_MATRIX)},
        {"PERMEABILITY_MATRIX", MakeIntegrationPointResultWriterFor(PERMEABILITY_MATRIX)},
        {"DEGREE_OF_SATURATION", MakeIntegrationPointResultWriterFor(DEGREE_OF_SATURATION)},
        {"DERIVATIVE_OF_SATURATION", MakeIntegrationPointResultWriterFor(DERIVATIVE_OF_SATURATION)},
        {"RELATIVE_PERMEABILITY", MakeIntegrationPointResultWriterFor(RELATIVE_PERMEABILITY)},
        {"PIPE_ACTIVE", MakeIntegrationPointResultWriterFor(PIPE_ACTIVE)},
        {"PIPE_HEIGHT", MakeIntegrationPointResultWriterFor(PIPE_HEIGHT)},
        {"GREEN_LAGRANGE_STRAIN_TENSOR", MakeIntegrationPointResultWriterFor(GREEN_LAGRANGE_STRAIN_TENSOR)},
        {"ENGINEERING_STRAIN_TENSOR", MakeIntegrationPointResultWriterFor(ENGINEERING_STRAIN_TENSOR)},
        {"CAUCHY_STRESS_TENSOR", MakeIntegrationPointResultWriterFor(CAUCHY_STRESS_TENSOR)},
        {"TOTAL_STRESS_TENSOR", MakeIntegrationPointResultWriterFor(TOTAL_STRESS_TENSOR)},
        {"VON_MISES_STRESS", MakeIntegrationPointResultWriterFor(VON_MISES_STRESS)}};

    for (const auto& name : rOutputItemNames) {
        auto iter = output_writer_map.find(name);
        KRATOS_ERROR_IF(iter == output_writer_map.end())
            << "Output item '" << name << "' is not available for integration point output" << std::endl;

        iter->second(mGidIO, rModelPart);
    }
}

void GeoOutputWriter::CalculateNodalHydraulicHead(ModelPart& rModelPart)
{
    const auto& element_var = KratosComponents<Variable<double>>::Get("HYDRAULIC_HEAD");

    for (Element element : rModelPart.Elements()) {
        auto&       rGeom = element.GetGeometry();
        const auto& rProp = element.GetProperties();
        const auto  NodalHydraulicHead =
            GeoElementUtilities::CalculateNodalHydraulicHeadFromWaterPressures(rGeom, rProp);

        for (unsigned int node = 0; node < rGeom.PointsNumber(); ++node) {
            rGeom[node].SetValue(element_var, NodalHydraulicHead[node]);
        }
    }
    mGidIO.WriteNodalResultsNonHistorical(element_var, rModelPart.Nodes(), 0);
}

void GeoOutputWriter::FinalizeResults() { mGidIO.FinalizeResults(); }

GidIO<> GeoOutputWriter::MakeGidIO(const std::string& rWorkingDirectory, const Parameters& rGidOutputSettings)
{
    const auto gid_post_flags =
        rGidOutputSettings["postprocess_parameters"]["result_file_configuration"]["gidpost_flags"];
    return GidIO<>{rWorkingDirectory + "/" + rGidOutputSettings["output_name"].GetString(),
                   GetGiDPostModeFrom(gid_post_flags), GetMultiFileFlagFrom(gid_post_flags),
                   GetWriteDeformedMeshFlagFrom(gid_post_flags), GetWriteConditionsFlagFrom(gid_post_flags)};
}

} // namespace Kratos
