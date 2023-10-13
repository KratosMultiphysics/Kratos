//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes
#include <iostream>
#include <fstream>

// External includes

// Project includes
#include "output_eigen_values_process.h"

namespace Kratos
{

OutputEigenValuesProcess::OutputEigenValuesProcess(
    Model& rModel,
    Parameters ThisParameters)
    : Process()
    , mrModel(rModel)
    , mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
}

void OutputEigenValuesProcess::ExecuteFinalize()
{
    std::string model_part_name = mThisParameters["model_part_name"].GetString();
    ModelPart& r_model_part = mrModel.GetModelPart(model_part_name);
    std::string output_file_name = mThisParameters["output_file_name"].GetString();

    std::string contents = "Post Results File 1.0\n";

    auto node_1 = r_model_part.NodesBegin();
    const Variable<Matrix>& eigenvector_variable = KratosComponents<Variable<Matrix>>::Get("EIGENVECTOR_MATRIX");
    Matrix eigenvector_matrix = node_1->GetValue(eigenvector_variable);

    std::vector<std::string> eigenvector_list(eigenvector_matrix.size1());

    for (IndexType i = 0; i < eigenvector_matrix.size1(); ++i)
    {
        eigenvector_list[i] = "Result \"DISPLACEMENT\" \"Load Case\" " + std::to_string(i + 1) + " Vector OnNodes\nValues\n";
    }

    for (auto& r_node : r_model_part.Nodes()) {
        Matrix eigenvector_matrix_node = r_node.GetValue(eigenvector_variable);
        for (IndexType i = 0; i < eigenvector_matrix.size1(); ++i)
        {
            eigenvector_list[i] += std::to_string(r_node.Id())
                + "  " + std::to_string(eigenvector_matrix_node(i, 0))
                + "  " + std::to_string(eigenvector_matrix_node(i, 1))
                + "  " + std::to_string(eigenvector_matrix_node(i, 2))
                + "\n";
        }
    }

    for (IndexType i = 0; i < eigenvector_matrix.size1(); ++i)
    {
        eigenvector_list[i] += "End Values\n";

        contents += eigenvector_list[i];
    }

    std::ofstream file;
    file.open(output_file_name);
    file << contents;
    file.close();
}

const Parameters OutputEigenValuesProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "output_file_name"           : "",
        "model_part_name"            : ""
    })" );
    return default_parameters;
}

} // namespace Kratos
