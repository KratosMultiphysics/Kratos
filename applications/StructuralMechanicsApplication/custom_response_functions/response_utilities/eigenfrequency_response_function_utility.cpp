// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Fusseder Martin,
//                   Armin Geiser,
//                   Daniel Baumgaertner,
//                   Suneth Warnakulasuriya
//

// System includes
#include <iostream>
#include <string>
#include <cmath>
#include <tuple>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"
#include "utilities/variable_utils.h"

// Application includes
#include "finite_difference_utility.h"

// Include base h
#include "eigenfrequency_response_function_utility.h"

namespace Kratos
{

EigenfrequencyResponseFunctionUtility::EigenfrequencyResponseFunctionUtility(
    ModelPart& rModelPart,
    Parameters ResponseSettings)
    : mrModelPart(rModelPart)
{
    CheckSettingsForGradientAnalysis(ResponseSettings);
    DetermineTracedEigenfrequencies(ResponseSettings);
    if(AreSeveralEigenfrequenciesTraced())
    {
        KRATOS_ERROR_IF_NOT(ResponseSettings.Has("weighting_method")) << "Several eigenfrequencies are traced but no weighting method specified in the parameters!" << std::endl;
        if(ResponseSettings["weighting_method"].GetString() == "linear_scaling")
            CalculateLinearWeights(ResponseSettings);
        else
            KRATOS_ERROR  << "Specified weighting method for eigenfrequencies is not implemented! Available weighting methods are: 'linear_scaling'." << std::endl;
    }
    else
        UseDefaultWeight();
}

EigenfrequencyResponseFunctionUtility::EigenfrequencyResponseFunctionUtility(
    ModelPart& rModelPart,
    const double PerturbationSize,
    const std::vector<IndexType>& TracedEigenFrequencyIds,
    const std::vector<double>& WeightingFactors)
    : mrModelPart(rModelPart),
      mDelta(PerturbationSize),
      mTracedEigenfrequencyIds(TracedEigenFrequencyIds),
      mWeightingFactors(WeightingFactors)
{
}

double EigenfrequencyResponseFunctionUtility::CalculateValue()
{
    CheckIfAllNecessaryEigenvaluesAreComputed();

    double resp_function_value = 0.0;
    KRATOS_INFO("EigenfrequencyResponseFunctionUtility") << "CalculateValue:" << std::endl
    << "    #    Eigenfrequency [Hz]    weighting factor" <<std::endl;
    for(std::size_t i = 0; i < mTracedEigenfrequencyIds.size(); i++)
    {
        const double eigenfrequency = std::sqrt(GetEigenvalue(mTracedEigenfrequencyIds[i])) / (2*Globals::Pi);
        resp_function_value += mWeightingFactors[i] * eigenfrequency;
        std::stringstream msg;
        msg << std::setw(5) << mTracedEigenfrequencyIds[i]
            << std::setw(23) << eigenfrequency
            << std::setw(20) << mWeightingFactors[i]<< std::endl;
        KRATOS_INFO("") << msg.str();
    }

    return resp_function_value;
}

void EigenfrequencyResponseFunctionUtility::CalculateGradient()
{
    CheckIfAllNecessaryEigenvaluesAreComputed();
    VariableUtils().SetHistoricalVariableToZero(SHAPE_SENSITIVITY, mrModelPart.Nodes());
    PerformSemiAnalyticSensitivityAnalysis();
}

void EigenfrequencyResponseFunctionUtility::CalculateEigenFrequencyMaterialVariableSensitivity(
    const Variable<double>& rMaterialVariable,
    const Variable<double>& rMaterialSensitivityVariable)
{
    using tls_type = std::tuple<Matrix, Matrix, Matrix, Vector, Matrix, Matrix, std::vector<Vector>>;

    CheckIfAllNecessaryEigenvaluesAreComputed();

    // get gradient prefactors
    const IndexType num_of_traced_eigenfrequencies = mTracedEigenfrequencyIds.size();

    Vector traced_eigenvalues(num_of_traced_eigenfrequencies);
    Vector gradient_prefactors(num_of_traced_eigenfrequencies);

    for(IndexType i = 0; i < num_of_traced_eigenfrequencies; ++i) {
        traced_eigenvalues[i] = GetEigenvalue(mTracedEigenfrequencyIds[i]);
        gradient_prefactors[i] = 1.0 / (4.0 * Globals::Pi * std::sqrt(traced_eigenvalues[i]));
    }

    const auto& r_process_info = mrModelPart.GetProcessInfo();

    block_for_each(mrModelPart.Elements(), tls_type(), [&](ModelPart::ElementType& rElement, tls_type& rTLS) {
        auto& r_properties = rElement.GetProperties();

        KRATOS_ERROR_IF_NOT(r_properties.Has(rMaterialVariable))
            << rMaterialVariable.Name() << " not found in element with id "
            << rElement.Id() << " in " << mrModelPart.FullName() << ".\n";

        Matrix& r_ref_lhs = std::get<0>(rTLS);
        Matrix& r_ref_mass = std::get<1>(rTLS);
        Matrix& r_aux_matrix = std::get<2>(rTLS);
        Vector& r_aux_vector = std::get<3>(rTLS);
        Matrix& r_perturbed_lhs = std::get<4>(rTLS);
        Matrix& r_perturbed_mass = std::get<5>(rTLS);
        std::vector<Vector> r_eigenvectors_of_element = std::get<6>(rTLS);

        // calculate reference matrices
        rElement.CalculateLeftHandSide(r_ref_lhs, r_process_info);
        rElement.CalculateMassMatrix(r_ref_mass, r_process_info);

        // calculate the perturbed matrices
        r_properties[rMaterialVariable] += mDelta;

        rElement.CalculateLeftHandSide(r_perturbed_lhs, r_process_info);
        rElement.CalculateMassMatrix(r_perturbed_mass, r_process_info);

        r_properties[rMaterialVariable] -= mDelta;

        if (r_aux_matrix.size1() != r_perturbed_lhs.size1() || r_aux_matrix.size2() != r_perturbed_lhs.size2()) {
            r_aux_matrix.resize(r_perturbed_lhs.size1(), r_perturbed_lhs.size2(), false);
            r_aux_vector.resize(r_perturbed_lhs.size1(), false);
            r_eigenvectors_of_element.resize(num_of_traced_eigenfrequencies);
        }

        double sensitivity = 0.0;
        for(IndexType i = 0; i < num_of_traced_eigenfrequencies; ++i) {
            DetermineEigenvectorOfElement(rElement, mTracedEigenfrequencyIds[i], r_eigenvectors_of_element[i], r_process_info);
            noalias(r_aux_matrix) = ((r_perturbed_lhs - r_ref_lhs) - (r_perturbed_mass - r_ref_mass) * traced_eigenvalues[i]) / mDelta;
            noalias(r_aux_vector) = prod(r_aux_matrix , r_eigenvectors_of_element[i]);

            sensitivity += gradient_prefactors[i] * inner_prod(r_eigenvectors_of_element[i], r_aux_vector) * mWeightingFactors[i];
        }

        r_properties.SetValue(rMaterialSensitivityVariable, sensitivity);
    });

}

void EigenfrequencyResponseFunctionUtility::CheckSettingsForGradientAnalysis(Parameters rResponseSettings)
{
    const std::string gradient_mode = rResponseSettings["gradient_mode"].GetString();

    if (gradient_mode == "semi_analytic")
        mDelta = rResponseSettings["step_size"].GetDouble();
    else
        KRATOS_ERROR << "Specified gradient_mode '" << gradient_mode << "' not recognized. The only option is: semi_analytic" << std::endl;
}

void EigenfrequencyResponseFunctionUtility::DetermineTracedEigenfrequencies(Parameters rResponseSettings)
{
    mTracedEigenfrequencyIds.resize(rResponseSettings["traced_eigenfrequencies"].size(),false);
    for(std::size_t i = 0; i < mTracedEigenfrequencyIds.size(); i++)
        mTracedEigenfrequencyIds[i] = rResponseSettings["traced_eigenfrequencies"][i].GetInt();
}

bool EigenfrequencyResponseFunctionUtility::AreSeveralEigenfrequenciesTraced()
{
    return (mTracedEigenfrequencyIds.size()>1);
}

void EigenfrequencyResponseFunctionUtility::CalculateLinearWeights(Parameters rResponseSettings)
{
    KRATOS_ERROR_IF_NOT(rResponseSettings.Has("weighting_factors")) << "No weighting factors defined for given eigenfrequency response!" << std::endl;
    KRATOS_ERROR_IF_NOT(rResponseSettings["weighting_factors"].size() == mTracedEigenfrequencyIds.size()) << "The number of chosen eigenvalues does not fit to the number of weighting factors!" << std::endl;

    mWeightingFactors.resize(rResponseSettings["weighting_factors"].size(),false);

    // Read weighting factors and sum them up
    double test_sum_weight_facs = 0.0;
    for(std::size_t i = 0; i < mWeightingFactors.size(); i++)
    {
        mWeightingFactors[i] = rResponseSettings["weighting_factors"][i].GetDouble();
        test_sum_weight_facs += mWeightingFactors[i];
    }

    // Check the weighting factors and modify them for the case that their sum is unequal to one
    if(std::abs(test_sum_weight_facs - 1.0) > 1e-12)
    {
        for(std::size_t i = 0; i < mWeightingFactors.size(); i++)
            mWeightingFactors[i] /= test_sum_weight_facs;

        KRATOS_INFO("\n> EigenfrequencyResponseFunctionUtility") << "The sum of the chosen weighting factors is unequal one. A corresponding scaling process was exected!" << std::endl;
    }
}

void EigenfrequencyResponseFunctionUtility::UseDefaultWeight()
{
    mWeightingFactors.resize(1,false);
    mWeightingFactors[0] = 1.0;
}

void EigenfrequencyResponseFunctionUtility::CheckIfAllNecessaryEigenvaluesAreComputed()
{
    const int num_of_computed_eigenvalues = (mrModelPart.GetProcessInfo()[EIGENVALUE_VECTOR]).size();
    const int max_required_eigenvalue = *(std::max_element(mTracedEigenfrequencyIds.begin(), mTracedEigenfrequencyIds.end()));
    KRATOS_ERROR_IF(max_required_eigenvalue > num_of_computed_eigenvalues) << "The following Eigenfrequency shall be traced but their corresponding eigenvalue was not computed by the eigenvalue analysies: " << max_required_eigenvalue << std::endl;
}

void EigenfrequencyResponseFunctionUtility::PerformSemiAnalyticSensitivityAnalysis()
{
    const ProcessInfo &CurrentProcessInfo = mrModelPart.GetProcessInfo();

    // Predetermine all necessary eigenvalues and prefactors for gradient calculation
    const std::size_t num_of_traced_eigenfrequencies = mTracedEigenfrequencyIds.size();
    Vector traced_eigenvalues(num_of_traced_eigenfrequencies,0.0);
    Vector gradient_prefactors(num_of_traced_eigenfrequencies,0.0);
    for(std::size_t i = 0; i < num_of_traced_eigenfrequencies; i++)
    {
        traced_eigenvalues[i] = GetEigenvalue(mTracedEigenfrequencyIds[i]);
        gradient_prefactors[i] = 1.0 / (4.0 * Globals::Pi * std::sqrt(traced_eigenvalues[i]));
    }

    // Element-wise computation of gradients
    for (auto& elem_i : mrModelPart.Elements())
    {
        Matrix LHS;
        Matrix mass_matrix;
        elem_i.CalculateLeftHandSide(LHS, CurrentProcessInfo);
        elem_i.CalculateMassMatrix(mass_matrix, CurrentProcessInfo);

        const std::size_t num_dofs_element = mass_matrix.size1();
        const std::size_t domain_size = CurrentProcessInfo.GetValue(DOMAIN_SIZE);

        Matrix aux_matrix = Matrix(num_dofs_element,num_dofs_element);
        Vector aux_vector = Vector(num_dofs_element);

        // Predetermine the necessary eigenvectors
        std::vector<Vector> eigenvectors_of_element(num_of_traced_eigenfrequencies,Vector(0));
        for(std::size_t i = 0; i < num_of_traced_eigenfrequencies; i++)
            DetermineEigenvectorOfElement(elem_i, mTracedEigenfrequencyIds[i], eigenvectors_of_element[i], CurrentProcessInfo);

        const std::vector<const FiniteDifferenceUtility::array_1d_component_type*> coord_directions = {&SHAPE_SENSITIVITY_X, &SHAPE_SENSITIVITY_Y, &SHAPE_SENSITIVITY_Z};

        // Computation of derivative of state equation w.r.t. node coordinates
        for(auto& node_i : elem_i.GetGeometry())
        {
            Vector gradient_contribution(3, 0.0);
            Matrix derived_LHS = Matrix(num_dofs_element,num_dofs_element);
            Matrix derived_mass_matrix = Matrix(num_dofs_element,num_dofs_element);

            for(std::size_t coord_dir_i = 0; coord_dir_i < domain_size; coord_dir_i++)
            {
                FiniteDifferenceUtility::CalculateLeftHandSideDerivative(elem_i, LHS, *coord_directions[coord_dir_i], node_i, mDelta, derived_LHS, CurrentProcessInfo);
                FiniteDifferenceUtility::CalculateMassMatrixDerivative(elem_i, mass_matrix, *coord_directions[coord_dir_i], node_i, mDelta, derived_mass_matrix, CurrentProcessInfo);

                for(std::size_t i = 0; i < num_of_traced_eigenfrequencies; i++)
                {
                    aux_matrix.clear();
                    aux_vector.clear();

                    noalias(aux_matrix) = derived_LHS - derived_mass_matrix * traced_eigenvalues[i];
                    noalias(aux_vector) = prod(aux_matrix , eigenvectors_of_element[i]);

                    gradient_contribution[coord_dir_i] += gradient_prefactors[i] * inner_prod(eigenvectors_of_element[i] , aux_vector) * mWeightingFactors[i];
                }
            }
            noalias(node_i.FastGetSolutionStepValue(SHAPE_SENSITIVITY)) += gradient_contribution;
        }
    }
}

double EigenfrequencyResponseFunctionUtility::GetEigenvalue(const IndexType EigenFrequencyId)
{
    return (mrModelPart.GetProcessInfo()[EIGENVALUE_VECTOR])[EigenFrequencyId-1];
}

void EigenfrequencyResponseFunctionUtility::DetermineEigenvectorOfElement(
    const ModelPart::ElementType& rElement,
    const IndexType EigenFrequencyId,
    Vector& rEigenvectorOfElement,
    const ProcessInfo& rCurrentProcessInfo)
{
    std::vector<std::size_t> eq_ids;
    rElement.EquationIdVector(eq_ids, rCurrentProcessInfo);

    if (rEigenvectorOfElement.size() != eq_ids.size())
    {
        rEigenvectorOfElement.resize(eq_ids.size(), false);
    }

    // sort the values of the eigenvector into the rEigenvectorOfElement according to the dof ordering at the element
    for (auto& r_node_i : rElement.GetGeometry())
    {
        const auto& r_node_dofs = r_node_i.GetDofs();

        const Matrix& rNodeEigenvectors = r_node_i.GetValue(EIGENVECTOR_MATRIX);
        for (std::size_t dof_index = 0; dof_index < r_node_dofs.size(); dof_index++)
        {
            const auto& current_dof = *(std::begin(r_node_dofs) + dof_index);
            const std::size_t dof_index_at_element = std::distance(eq_ids.begin(), std::find(eq_ids.begin(), eq_ids.end(), current_dof->EquationId()));
            rEigenvectorOfElement(dof_index_at_element) = rNodeEigenvectors((EigenFrequencyId-1), dof_index);
        }
    }
}

} // namespace Kratos.

