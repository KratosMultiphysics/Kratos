// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Jonathan Nuttall
//

#pragma once

#include <geo_mechanics_application.h>
#include "dgeooutput.h"

void NodeOperation::write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) {};

void NodeDISPLACEMENT::write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.WriteNodalResults(Kratos::DISPLACEMENT, model_part.Nodes(), 0, 0); }

void NodeTOTAL_DISPLACEMENT::write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.WriteNodalResults(Kratos::TOTAL_DISPLACEMENT, model_part.Nodes(), 0, 0); }

void NodeWATER_PRESSURE::write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.WriteNodalResults(Kratos::WATER_PRESSURE, model_part.Nodes(), 0, 0); }

void NodeNORMAL_FLUID_FLUX::write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.WriteNodalResults(Kratos::NORMAL_FLUID_FLUX, model_part.Nodes(), 0, 0); }

void NodeVOLUME_ACCELERATION::write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.WriteNodalResults(Kratos::VOLUME_ACCELERATION, model_part.Nodes(), 0, 0); }

void NodeHYDRAULIC_DISCHARGE::write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.WriteNodalResults(Kratos::HYDRAULIC_DISCHARGE, model_part.Nodes(), 0, 0); }

void NodeHYDRAULIC_HEAD::write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.WriteNodalResults(Kratos::HYDRAULIC_HEAD, model_part.Nodes(), 0, 0); }

void GaussOperation::write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) {};

void GaussFLUID_FLUX_VECTOR::write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.PrintOnGaussPoints(Kratos::FLUID_FLUX_VECTOR, model_part, 0, 0); }

void GaussHYDRAULIC_HEAD::write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.PrintOnGaussPoints(Kratos::HYDRAULIC_HEAD, model_part, 0, 0); }

void GaussLOCAL_FLUID_FLUX_VECTOR::write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.PrintOnGaussPoints(Kratos::LOCAL_FLUID_FLUX_VECTOR, model_part, 0, 0); }

void GaussLOCAL_PERMEABILITY_MATRIX::write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.PrintOnGaussPoints(Kratos::LOCAL_PERMEABILITY_MATRIX, model_part, 0, 0); }

void GaussPERMEABILITY_MATRIX::write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.PrintOnGaussPoints(Kratos::PERMEABILITY_MATRIX, model_part, 0, 0); }

void GaussDEGREE_OF_SATURATION::write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.PrintOnGaussPoints(Kratos::DEGREE_OF_SATURATION, model_part, 0, 0); }

void GaussDERIVATIVE_OF_SATURATION::write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.PrintOnGaussPoints(Kratos::DERIVATIVE_OF_SATURATION, model_part, 0, 0); }

void GaussRELATIVE_PERMEABILITY::write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.PrintOnGaussPoints(Kratos::RELATIVE_PERMEABILITY, model_part, 0, 0); }

void GaussPIPE_ACTIVE::write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.PrintOnGaussPoints(Kratos::PIPE_ACTIVE, model_part, 0, 0); }

void GaussPIPE_HEIGHT::write(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) { gid_io.PrintOnGaussPoints(Kratos::PIPE_HEIGHT, model_part, 0, 0); }

namespace
{
    void calculateNodalHydraulicHead(Kratos::GidIO<>& gid_io, Kratos::ModelPart& model_part) {
        auto element_var = &(Kratos::KratosComponents<Kratos::Variable<double>>::Get("HYDRAULIC_HEAD"));

        for (Kratos::Element element : model_part.Elements())
        {
            auto rGeom = element.GetGeometry();
            auto rProp = element.GetProperties();

            for (unsigned int node = 0; node < 3; ++node)
            {
                Kratos::array_1d<double, 3> NodeVolumeAcceleration;
                noalias(NodeVolumeAcceleration) = rGeom[node].FastGetSolutionStepValue(Kratos::VOLUME_ACCELERATION, 0);
                const double g = norm_2(NodeVolumeAcceleration);
                if (g > std::numeric_limits<double>::epsilon())
                {
                    const double FluidWeight = g * rProp[Kratos::DENSITY_WATER];

                    Kratos::array_1d<double, 3> NodeCoordinates;
                    noalias(NodeCoordinates) = rGeom[node].Coordinates();
                    Kratos::array_1d<double, 3> NodeVolumeAccelerationUnitVector;
                    noalias(NodeVolumeAccelerationUnitVector) = NodeVolumeAcceleration / g;

                    const double WaterPressure = rGeom[node].FastGetSolutionStepValue(Kratos::WATER_PRESSURE);
                    rGeom[node].SetValue(*element_var, -inner_prod(NodeCoordinates, NodeVolumeAccelerationUnitVector) - Kratos::PORE_PRESSURE_SIGN_FACTOR * WaterPressure / FluidWeight);
                }
                else
                {
                    rGeom[node].SetValue(*element_var, 0.0);
                }
            }
        }

        gid_io.WriteNodalResultsNonHistorical(*element_var, model_part.Nodes(), 0);
    }
}

namespace Kratos
{
    void KratosGeoOutput::outputGiD(Model &model, ModelPart &model_part, Parameters parameters, std::string workingDirectory)
    {
        std::map<std::string, GiD_PostMode> PostMode;
        PostMode["GiD_PostAscii"] = GiD_PostAscii;
        PostMode["GiD_PostAsciiZipped"] = GiD_PostAsciiZipped;
        PostMode["GiD_PostBinary"] = GiD_PostBinary;
        PostMode["GiD_PostHDF5"] = GiD_PostHDF5;

        std::map<std::string, Kratos::MultiFileFlag> MultiFiles;
        MultiFiles["SingleFile"] = Kratos::SingleFile;
        MultiFiles["MultipleFiles"] = Kratos::MultipleFiles;

        std::map<std::string, Kratos::WriteDeformedMeshFlag> DeformedFlag;
        DeformedFlag["WriteDeformed"] = Kratos::WriteDeformed;
        DeformedFlag["WriteUndeformed"] = Kratos::WriteUndeformed;

        std::map<std::string, Kratos::WriteConditionsFlag> ConditionFlag;
        ConditionFlag["WriteConditions"] = Kratos::WriteConditions;
        ConditionFlag["WriteElementsOnly"] = Kratos::WriteElementsOnly;
        ConditionFlag["WriteConditionsOnly"] = Kratos::WriteConditionsOnly;

        Parameters gid_out = parameters["output_processes"]["gid_output"].GetArrayItem(0);
        Parameters outputParameters = gid_out["Parameters"];
        std::string filename = outputParameters["output_name"].GetString();
        GiD_PostMode gid_output_type = PostMode[outputParameters["postprocess_parameters"]["result_file_configuration"]["gidpost_flags"]["GiDPostMode"].GetString()];
        MultiFileFlag multifiles_output = MultiFiles[outputParameters["postprocess_parameters"]["result_file_configuration"]["gidpost_flags"]["MultiFileFlag"].GetString()];
        WriteDeformedMeshFlag deformed_output = DeformedFlag[outputParameters["postprocess_parameters"]["result_file_configuration"]["gidpost_flags"]["WriteDeformedMeshFlag"].GetString()];
        WriteConditionsFlag condition_output = ConditionFlag[outputParameters["postprocess_parameters"]["result_file_configuration"]["gidpost_flags"]["WriteConditionsFlag"].GetString()];

        filename = workingDirectory + "/" + filename;
        GidIO<> gid_io(filename, gid_output_type, multifiles_output, deformed_output, condition_output);

        gid_io.InitializeMesh(0.0);
        gid_io.WriteMesh(model_part.GetMesh());
        gid_io.FinalizeMesh();
        gid_io.InitializeResults(0, model_part.GetMesh());

        std::unordered_map<std::string, unique_ptr<NodeOperation>> NodeOutput;
        NodeOutput["DISPLACEMENT"] = make_unique<NodeDISPLACEMENT>();
        NodeOutput["TOTAL_DISPLACEMENT"] = make_unique<NodeTOTAL_DISPLACEMENT>();
        NodeOutput["WATER_PRESSURE"] = make_unique<NodeWATER_PRESSURE>();
        NodeOutput["NORMAL_FLUID_FLUX"] = make_unique<NodeNORMAL_FLUID_FLUX>();
        NodeOutput["VOLUME_ACCELERATION"] = make_unique<NodeVOLUME_ACCELERATION>();
        NodeOutput["HYDRAULIC_DISCHARGE"] = make_unique<NodeHYDRAULIC_DISCHARGE>();
        NodeOutput["HYDRAULIC_HEAD"] = make_unique<NodeHYDRAULIC_HEAD>();

        // Calculate hydraulic head on the nodes
        auto gauss_outputs = outputParameters["postprocess_parameters"]["result_file_configuration"]["gauss_point_results"].GetStringArray();
        auto nodal_outputs = outputParameters["postprocess_parameters"]["result_file_configuration"]["nodal_results"].GetStringArray();

        if (std::find(gauss_outputs.begin(), gauss_outputs.end(), "HYDRAULIC_HEAD") != gauss_outputs.end())
        {
            calculateNodalHydraulicHead(gid_io, model_part);
        }

        for (std::string var : nodal_outputs)
        {
            NodeOutput[var]->write(gid_io, model_part);
        }

        std::unordered_map<std::string, std::unique_ptr<GaussOperation>> GaussOutput;
        GaussOutput["FLUID_FLUX_VECTOR"] = make_unique<GaussFLUID_FLUX_VECTOR>();
        GaussOutput["HYDRAULIC_HEAD"] = make_unique<GaussHYDRAULIC_HEAD>();
        GaussOutput["LOCAL_FLUID_FLUX_VECTOR"] = make_unique<GaussLOCAL_FLUID_FLUX_VECTOR>();
        GaussOutput["LOCAL_PERMEABILITY_MATRIX"] = make_unique<GaussLOCAL_PERMEABILITY_MATRIX>();
        GaussOutput["PERMEABILITY_MATRIX"] = make_unique<GaussPERMEABILITY_MATRIX>();
        GaussOutput["DEGREE_OF_SATURATION"] = make_unique<GaussDEGREE_OF_SATURATION>();
        GaussOutput["DERIVATIVE_OF_SATURATION"] = make_unique<GaussDERIVATIVE_OF_SATURATION>();
        GaussOutput["RELATIVE_PERMEABILITY"] = make_unique<GaussRELATIVE_PERMEABILITY>();
        GaussOutput["PIPE_ACTIVE"] = make_unique<GaussPIPE_ACTIVE>();
        GaussOutput["PIPE_HEIGHT"] = make_unique<GaussPIPE_HEIGHT>();

        // Now Output Gauss Point Results on Gauss Points
        for (std::string var : gauss_outputs)
        {
            GaussOutput[var]->write(gid_io, model_part);
        }

        gid_io.FinalizeResults();
    }

}
