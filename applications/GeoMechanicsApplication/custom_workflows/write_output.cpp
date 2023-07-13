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

class NodeOperation
{
  public:
    virtual ~NodeOperation() = default;
    virtual void write(Kratos::GidIO<>& rGidIO, Kratos::ModelPart& rModelPart) = 0;
};

class NodeDISPLACEMENT : public NodeOperation
{
  public:
    void write(Kratos::GidIO<>& rGidIO, Kratos::ModelPart& rModelPart) override
    {
        rGidIO.WriteNodalResults(Kratos::DISPLACEMENT, rModelPart.Nodes(), 0, 0);
    }
};

class NodeTOTAL_DISPLACEMENT : public NodeOperation
{
  public:
    void write(Kratos::GidIO<>& rGidIO, Kratos::ModelPart& rModelPart) override
    {
        rGidIO.WriteNodalResults(Kratos::TOTAL_DISPLACEMENT, rModelPart.Nodes(), 0, 0);
    }
};

class NodeWATER_PRESSURE : public NodeOperation
{
  public:
    void write(Kratos::GidIO<>& rGidIO, Kratos::ModelPart& rModelPart) override
    {
        rGidIO.WriteNodalResults(Kratos::WATER_PRESSURE, rModelPart.Nodes(), 0, 0);
    }
};

class NodeNORMAL_FLUID_FLUX : public NodeOperation
{
  public:
    void write(Kratos::GidIO<>& rGidIO, Kratos::ModelPart& rModelPart) override
    {
        rGidIO.WriteNodalResults(Kratos::NORMAL_FLUID_FLUX, rModelPart.Nodes(), 0, 0);
    }
};

class NodeVOLUME_ACCELERATION : public NodeOperation
{
  public:
    void write(Kratos::GidIO<>& rGidIO, Kratos::ModelPart& rModelPart) override
    {
        rGidIO.WriteNodalResults(Kratos::VOLUME_ACCELERATION, rModelPart.Nodes(), 0, 0);
    }
};

class NodeHYDRAULIC_DISCHARGE : public NodeOperation
{
  public:
    void write(Kratos::GidIO<>& rGidIO, Kratos::ModelPart& rModelPart) override
    {
        rGidIO.WriteNodalResults(Kratos::HYDRAULIC_DISCHARGE, rModelPart.Nodes(), 0, 0);
    }
};

class NodeHYDRAULIC_HEAD : public NodeOperation
{
  public:
    void write(Kratos::GidIO<>& rGidIO, Kratos::ModelPart& rModelPart) override
    {
        rGidIO.WriteNodalResults(Kratos::HYDRAULIC_HEAD, rModelPart.Nodes(), 0, 0);
    }
};


class GaussOperation
{
  public:
    virtual ~GaussOperation() = default;
    virtual void write(Kratos::GidIO<>& rGidIO, Kratos::ModelPart& rModelPart) = 0;
};

class GaussFLUID_FLUX_VECTOR : public GaussOperation
{
  public:
    void write(Kratos::GidIO<>& rGidIO, Kratos::ModelPart& rModelPart) override
    {
        rGidIO.PrintOnGaussPoints(Kratos::FLUID_FLUX_VECTOR, rModelPart, 0, 0);
    }
};

class GaussHYDRAULIC_HEAD : public GaussOperation
{
  public:
    void write(Kratos::GidIO<>& rGidIO, Kratos::ModelPart& rModelPart) override
    {
        rGidIO.PrintOnGaussPoints(Kratos::HYDRAULIC_HEAD, rModelPart, 0, 0);
    }
};

class GaussLOCAL_FLUID_FLUX_VECTOR : public GaussOperation
{
  public:
    void write(Kratos::GidIO<>& rGidIO, Kratos::ModelPart& rModelPart) override
    {
        rGidIO.PrintOnGaussPoints(Kratos::LOCAL_FLUID_FLUX_VECTOR, rModelPart, 0, 0);
    }
};

class GaussLOCAL_PERMEABILITY_MATRIX : public GaussOperation
{
  public:
    void write(Kratos::GidIO<>& rGidIO, Kratos::ModelPart& rModelPart) override
    {
        rGidIO.PrintOnGaussPoints(Kratos::LOCAL_PERMEABILITY_MATRIX, rModelPart, 0, 0);
    }
};

class GaussPERMEABILITY_MATRIX : public GaussOperation
{
  public:
    void write(Kratos::GidIO<>& rGidIO, Kratos::ModelPart& rModelPart) override
    {
        rGidIO.PrintOnGaussPoints(Kratos::PERMEABILITY_MATRIX, rModelPart, 0, 0);
    }
};

class GaussDEGREE_OF_SATURATION : public GaussOperation
{
  public:
    void write(Kratos::GidIO<>& rGidIO, Kratos::ModelPart& rModelPart) override
    {
        rGidIO.PrintOnGaussPoints(Kratos::DEGREE_OF_SATURATION, rModelPart, 0, 0);
    }
};

class GaussDERIVATIVE_OF_SATURATION : public GaussOperation
{
  public:
    void write(Kratos::GidIO<>& rGidIO, Kratos::ModelPart& rModelPart) override
    {
        rGidIO.PrintOnGaussPoints(Kratos::DERIVATIVE_OF_SATURATION, rModelPart, 0, 0);
    }
};

class GaussRELATIVE_PERMEABILITY : public GaussOperation
{
  public:
    void write(Kratos::GidIO<>& rGidIO, Kratos::ModelPart& rModelPart) override
    {
        rGidIO.PrintOnGaussPoints(Kratos::RELATIVE_PERMEABILITY, rModelPart, 0, 0);
    }
};

class GaussPIPE_ACTIVE : public GaussOperation
{
  public:
    void write(Kratos::GidIO<>& rGidIO, Kratos::ModelPart& rModelPart) override
    {
        rGidIO.PrintOnGaussPoints(Kratos::PIPE_ACTIVE, rModelPart, 0, 0);
    }
};

class GaussPIPE_HEIGHT : public GaussOperation
{
  public:
    void write(Kratos::GidIO<>& rGidIO, Kratos::ModelPart& rModelPart) override
    {
        rGidIO.PrintOnGaussPoints(Kratos::PIPE_HEIGHT, rModelPart, 0, 0);
    }
};

}


namespace Kratos
{

void GeoOutputWriter::WriteGiDOutput(ModelPart&         rModelPart,
                                     Parameters         Settings,
                                     const std::string& rWorkingDirectory)
{
    auto output_parameters = Settings["output_processes"]["gid_output"].GetArrayItem(0)["Parameters"];
    auto gid_post_flags = output_parameters["postprocess_parameters"]["result_file_configuration"]["gidpost_flags"];

    GidIO<> gid_io{rWorkingDirectory + "/" + output_parameters["output_name"].GetString(),
                   GetGiDPostModeFrom(gid_post_flags),
                   GetMultiFileFlagFrom(gid_post_flags),
                   GetWriteDeformedMeshFlagFrom(gid_post_flags),
                   GetWriteConditionsFlagFrom(gid_post_flags)};

    gid_io.InitializeMesh(0.0);
    gid_io.WriteMesh(rModelPart.GetMesh());
    gid_io.FinalizeMesh();

    gid_io.InitializeResults(0, rModelPart.GetMesh());

    // Calculate hydraulic head on the nodes
    auto gauss_outputs = output_parameters["postprocess_parameters"]["result_file_configuration"]["gauss_point_results"].GetStringArray();
    if (std::find(gauss_outputs.begin(), gauss_outputs.end(), "HYDRAULIC_HEAD") != gauss_outputs.end())
    {
        CalculateNodalHydraulicHead(gid_io, rModelPart);
    }

    const auto nodal_outputs = output_parameters["postprocess_parameters"]["result_file_configuration"]["nodal_results"].GetStringArray();
    WriteNodalOutput(nodal_outputs, gid_io, rModelPart);
    WriteIntegrationPointOutput(gauss_outputs, gid_io, rModelPart);

    gid_io.FinalizeResults();
}


void GeoOutputWriter::WriteNodalOutput(const std::vector<std::string>& rOutputItemNames,
                                       GidIO<>&                        rGidIO,
                                       ModelPart&                      rModelPart)
{
    std::map<std::string, std::unique_ptr<NodeOperation>, std::less<>> output_writer_map;
    output_writer_map["DISPLACEMENT"]        = std::make_unique<NodeDISPLACEMENT>();
    output_writer_map["TOTAL_DISPLACEMENT"]  = std::make_unique<NodeTOTAL_DISPLACEMENT>();
    output_writer_map["WATER_PRESSURE"]      = std::make_unique<NodeWATER_PRESSURE>();
    output_writer_map["NORMAL_FLUID_FLUX"]   = std::make_unique<NodeNORMAL_FLUID_FLUX>();
    output_writer_map["VOLUME_ACCELERATION"] = std::make_unique<NodeVOLUME_ACCELERATION>();
    output_writer_map["HYDRAULIC_DISCHARGE"] = std::make_unique<NodeHYDRAULIC_DISCHARGE>();
    output_writer_map["HYDRAULIC_HEAD"]      = std::make_unique<NodeHYDRAULIC_HEAD>();

    for (const auto& name : rOutputItemNames)
    {
        output_writer_map.at(name)->write(rGidIO, rModelPart);
    }
}


void GeoOutputWriter::WriteIntegrationPointOutput(const std::vector<std::string>& rOutputItemNames,
                                                  GidIO<>&                        rGidIO,
                                                  ModelPart&                      rModelPart)
{
    std::map<std::string, std::unique_ptr<GaussOperation>, std::less<>> output_writer_map;
    output_writer_map["FLUID_FLUX_VECTOR"]         = std::make_unique<GaussFLUID_FLUX_VECTOR>();
    output_writer_map["HYDRAULIC_HEAD"]            = std::make_unique<GaussHYDRAULIC_HEAD>();
    output_writer_map["LOCAL_FLUID_FLUX_VECTOR"]   = std::make_unique<GaussLOCAL_FLUID_FLUX_VECTOR>();
    output_writer_map["LOCAL_PERMEABILITY_MATRIX"] = std::make_unique<GaussLOCAL_PERMEABILITY_MATRIX>();
    output_writer_map["PERMEABILITY_MATRIX"]       = std::make_unique<GaussPERMEABILITY_MATRIX>();
    output_writer_map["DEGREE_OF_SATURATION"]      = std::make_unique<GaussDEGREE_OF_SATURATION>();
    output_writer_map["DERIVATIVE_OF_SATURATION"]  = std::make_unique<GaussDERIVATIVE_OF_SATURATION>();
    output_writer_map["RELATIVE_PERMEABILITY"]     = std::make_unique<GaussRELATIVE_PERMEABILITY>();
    output_writer_map["PIPE_ACTIVE"]               = std::make_unique<GaussPIPE_ACTIVE>();
    output_writer_map["PIPE_HEIGHT"]               = std::make_unique<GaussPIPE_HEIGHT>();

    for (const auto& name : rOutputItemNames)
    {
        output_writer_map.at(name)->write(rGidIO, rModelPart);
    }
}


void GeoOutputWriter::CalculateNodalHydraulicHead(GidIO<>& rGidIO, ModelPart& rModelPart)
{
    const auto& element_var = KratosComponents<Variable<double>>::Get("HYDRAULIC_HEAD");

    for (Element element : rModelPart.Elements())
    {
        auto& rGeom = element.GetGeometry();
        const auto& rProp = element.GetProperties();

        const auto NodalHydraulicHead = GeoElementUtilities::CalculateNodalHydraulicHeadFromWaterPressures<3>(rGeom, rProp);

        for (unsigned int node = 0; node < 3; ++node)
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

}
