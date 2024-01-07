//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// System includes

// External includes
#include <Eigen/Dense>
#include <Eigen/SVD>

// Project includes
#include "includes/define.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"

// Application includes

// Include base h
#include "max_singular_value_decomposition_process.h"

namespace Kratos
{

template<class TContainterType>
void MaxSingularValueDecompositionProcess::CalculateAndStoreMaxSingularValues(
    TContainterType& rContainer,
    const ProcessInfo& rProcessInfo) const
{
    KRATOS_TRY

    const auto& r_input_variable = KratosComponents<Variable<Matrix>>::Get(mInputVariableName);
    const auto& r_output_variable = KratosComponents<Variable<double>>::Get(mOutputVariableName);

    using TLS = std::tuple<Eigen::MatrixXd, Matrix>;

    block_for_each(rContainer, TLS(), [&](typename TContainterType::data_type& rEntity, TLS& rTLS) {
        auto& eigen_matrix = std::get<0>(rTLS);
        auto& input_matrix = std::get<1>(rTLS);

        rEntity.Calculate(r_input_variable, input_matrix, rProcessInfo);

        if (eigen_matrix.rows() != static_cast<int>(input_matrix.size1()) || eigen_matrix.cols() != static_cast<int>(input_matrix.size2())) {
            eigen_matrix.resize(input_matrix.size1(), input_matrix.size2());
        }

        for (IndexType i = 0; i < input_matrix.size1(); ++i) {
            for (IndexType j = 0; j < input_matrix.size2(); ++j) {
                eigen_matrix(i, j) = input_matrix(i, j);
            }
        }

        // compute the singular values for current equation
        Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::NoQRPreconditioner> svd(eigen_matrix);

        // get the maximum singular value
        rEntity.SetValue(r_output_variable, svd.singularValues()[0]);
    });

    KRATOS_CATCH("");
}

// template instantiations
template void MaxSingularValueDecompositionProcess::CalculateAndStoreMaxSingularValues(ModelPart::ConditionsContainerType&, const ProcessInfo&) const;
template void MaxSingularValueDecompositionProcess::CalculateAndStoreMaxSingularValues(ModelPart::ElementsContainerType&, const ProcessInfo&) const;

MaxSingularValueDecompositionProcess::MaxSingularValueDecompositionProcess(
    Model& rModel,
    Parameters rParameters)
: mrModel(rModel)
{
    KRATOS_TRY

    rParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    mInputVariableName = rParameters["input_variable_name"].GetString();
    mOutputVariableName = rParameters["output_variable_name"].GetString();
    mModelPartName = rParameters["model_part_name"].GetString();
    mContainerType = rParameters["container_type"].GetString();
    mEchoLevel = rParameters["echo_level"].GetInt();

    KRATOS_ERROR_IF_NOT(mContainerType == "elements" || mContainerType == "conditions")
        << "\"container_type\" should be either \"conditions\" or \"elements\". [ container_type = \""
        << mContainerType << " \"].\n";

    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<Matrix>>::Has(mInputVariableName))
        << "Given \"input_variable_name\" is not found in matrix variables list. [ \"input_variable_name\" = " << mInputVariableName << " ].\n";

    KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(mOutputVariableName))
        << "Given \"output_variable_name\" is not found in double variables list. [ \"output_variable_name\" = " << mInputVariableName << " ].\n";

    KRATOS_CATCH("");
}

void MaxSingularValueDecompositionProcess::Execute()
{
    KRATOS_TRY

    auto& r_model_part = mrModel.GetModelPart(mModelPartName);
    const auto& r_process_info = r_model_part.GetProcessInfo();

    if (mContainerType == "elements") {
        CalculateAndStoreMaxSingularValues(r_model_part.Elements(), r_process_info);
    } else if (mContainerType == "conditions") {
        CalculateAndStoreMaxSingularValues(r_model_part.Conditions(), r_process_info);
    }

    KRATOS_INFO_IF(this->Info(), mEchoLevel > 0) << "Stored max singular value of the matrix given by " << mInputVariableName << " in " << mContainerType << " in " << r_model_part.FullName() << " to " << mOutputVariableName << " .\n";

    KRATOS_CATCH("");
}

std::string MaxSingularValueDecompositionProcess::Info() const
{
    return std::string("MaxSingularValueDecompositionProcess");
}

void MaxSingularValueDecompositionProcess::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info();
}

void MaxSingularValueDecompositionProcess::PrintData(std::ostream& rOStream) const
{
}

const Parameters MaxSingularValueDecompositionProcess::GetDefaultParameters() const
{
    const auto default_parameters = Parameters(R"(
        {
            "model_part_name"     : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "input_variable_name" : "PLEASE_SPECIFY_MATRIX_VARIABLE",
            "output_variable_name": "PLEASE_SPECIFY_SCALAR_VARIABLE",
            "echo_level"          : 0,
            "container_type"      : "elements"
        })");

    return default_parameters;
}

} // namespace Kratos.
