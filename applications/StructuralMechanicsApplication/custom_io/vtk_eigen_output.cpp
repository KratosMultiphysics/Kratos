//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//
//

// System includes

// External includes

// Project includes
#include "vtk_eigen_output.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{

void VtkEigenOutput::PrintEigenOutput(
        const std::string& rLabel,
        const int AnimationStep,
        const std::vector<Variable<double>>& rRequestedDoubleResults,
        const std::vector<Variable<array_1d<double,3>>>& rRequestedVectorResults)
{
    std::ofstream output_file;
    const std::string output_file_name = GetEigenOutputFileName(AnimationStep);

    std::ios::openmode ios_flags = std::ios::out; // basic flag always needed

    if (mFileFormat == VtkOutput::FileFormat::VTK_BINARY) { ios_flags |= std::ios::binary; }

    if (AnimationStep > mLastWrittenAnimationStepIndex) {
        // No file exists yet for this animationstep, creating a new one and writing the mesh
        mLastWrittenAnimationStepIndex = AnimationStep;

        ios_flags |= std::ios::trunc;
        OpenOutputFile(output_file_name, ios_flags, output_file);

        Initialize(mrModelPart);
        WriteHeaderToFile(mrModelPart, output_file);
        WriteMeshToFile(mrModelPart, output_file);
        const SizeType num_eigenvalues = mrModelPart.GetProcessInfo()[EIGENVALUE_VECTOR].size();

        output_file << "POINT_DATA " << mrModelPart.NumberOfNodes() << "\n";
        output_file << "FIELD FieldData " << num_eigenvalues * (rRequestedDoubleResults.size() + rRequestedVectorResults.size()) << "\n";
    } else {
        // Appending results to existing file
        ios_flags |= std::ios::app;
        OpenOutputFile(output_file_name, ios_flags, output_file);
    }

    for (const auto& r_variable : rRequestedDoubleResults) {
        WriteScalarEigenVariable(mrModelPart.Nodes(), r_variable, rLabel, output_file);
    }

    for (const auto& r_variable : rRequestedVectorResults) {
        WriteVectorEigenVariable(mrModelPart.Nodes(), r_variable, rLabel, output_file);
    }

    output_file.close();
}

void VtkEigenOutput::OpenOutputFile(
    const std::string& rFileName,
    const std::ios::openmode OpenModeFlags,
    std::ofstream& rOutputFile) const
{
    rOutputFile.open(rFileName, OpenModeFlags);
    if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) {
        rOutputFile << std::scientific << std::setprecision(mDefaultPrecision);
    }

    KRATOS_ERROR_IF_NOT(rOutputFile.is_open()) << "The output file could not be opened: " << rFileName << std::endl;
}

std::string VtkEigenOutput::GetEigenOutputFileName(const int AnimationStep) const
{
    std::string output_file_name = mEigenOutputSettings["result_file_name"].GetString();

    if (output_file_name == "") { // use the name of the ModelPart in case nothing was assigned
        output_file_name = mrModelPart.Name();
    }

    output_file_name += "_EigenResults_";

    const std::string file_label = mEigenOutputSettings["file_label"].GetString();
    if (file_label == "step") {
        output_file_name += std::to_string(mrModelPart.GetProcessInfo()[STEP]);
    } else if (file_label == "time") {
        output_file_name += std::to_string(mrModelPart.GetProcessInfo()[TIME]);
    } else {
        KRATOS_ERROR << "\"file_label\" can only be \"step\" or \"time\"" << std::endl;
    }

    output_file_name += "_" + std::to_string(AnimationStep) + ".vtk";

    if (mEigenOutputSettings["save_output_files_in_folder"].GetBool()) {
        output_file_name = mEigenOutputSettings["folder_name"].GetString() + "/" + output_file_name;
    }

    return output_file_name;
}

void VtkEigenOutput::WriteScalarEigenVariable(
    const ModelPart::NodesContainerType& rNodes,
    const Variable<double>& rVariable,
    const std::string& rLabel,
    std::ofstream& rFileStream) const
{
    rFileStream << rLabel << "_" << rVariable.Name() << " 1 " << rNodes.size() << " float\n";

    for (const auto& r_node : rNodes) {
        const auto& r_result = r_node.FastGetSolutionStepValue(rVariable);
        VtkOutput::WriteScalarDataToFile((float)r_result, rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}

void VtkEigenOutput::WriteVectorEigenVariable(
    const ModelPart::NodesContainerType& rNodes,
    const Variable<array_1d<double, 3>>& rVariable,
    const std::string& rLabel,
    std::ofstream& rFileStream) const
{
    if (rNodes.size() == 0) {
        KRATOS_WARNING("VtkEigenOutput") << "Empty container!" << std::endl;
        return;
    }

    const int res_size = static_cast<int>((rNodes.begin()->FastGetSolutionStepValue(rVariable)).size());

    rFileStream << rLabel << "_" << rVariable.Name() << " " << res_size << " " << rNodes.size() << "  float\n";

    for (const auto& r_node : rNodes) {
        const auto& r_result = r_node.FastGetSolutionStepValue(rVariable);
        VtkOutput::WriteVectorDataToFile(r_result, rFileStream);
        if (mFileFormat == VtkOutput::FileFormat::VTK_ASCII) rFileStream <<"\n";
    }
}

} // namespace Kratos
