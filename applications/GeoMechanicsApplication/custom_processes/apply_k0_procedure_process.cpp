// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//

#include "apply_k0_procedure_process.h"

#include <algorithm>
#include <ostream>

#include "containers/model.h"
#include "custom_constitutive/linear_elastic_law.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/process_utilities.h"
#include "geo_aliases.h"
#include "geo_mechanics_application_variables.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

namespace
{

using namespace Kratos;
using namespace std::string_literals;

void SetConsiderDiagonalEntriesOnlyAndNoShear(ModelPart::ElementsContainerType& rElements, bool Whether)
{
    block_for_each(rElements, [Whether](Element& rElement) {
        auto pLinearElasticLaw =
            dynamic_cast<GeoLinearElasticLaw*>(rElement.GetProperties().GetValue(CONSTITUTIVE_LAW).get());
        if (pLinearElasticLaw) pLinearElasticLaw->SetConsiderDiagonalEntriesOnlyAndNoShear(Whether);
    });
}

} // namespace

namespace Kratos
{

ApplyK0ProcedureProcess::ApplyK0ProcedureProcess(Model& rModel, Parameters K0Settings)
    : mSettings(std::move(K0Settings))
{
    mrModelParts =
        ProcessUtilities::GetModelPartsFromSettings(rModel, mSettings, ApplyK0ProcedureProcess::Info());
}

void ApplyK0ProcedureProcess::ExecuteInitialize()
{
    if (UseStandardProcedure()) {
        for (const auto& r_model_part : mrModelParts) {
            SetConsiderDiagonalEntriesOnlyAndNoShear(r_model_part.get().Elements(), true);
        }
    }
}

void ApplyK0ProcedureProcess::ExecuteFinalize()
{
    if (UseStandardProcedure()) {
        for (const auto& r_model_part : mrModelParts) {
            SetConsiderDiagonalEntriesOnlyAndNoShear(r_model_part.get().Elements(), false);
        }
    }
}

int ApplyK0ProcedureProcess::Check()
{
    for (const auto& r_model_part : mrModelParts) {
        KRATOS_ERROR_IF(r_model_part.get().Elements().empty())
            << "ApplyK0ProcedureProces has no elements in modelpart " << r_model_part.get().Name()
            << std::endl;
    }

    for (const auto& r_model_part : mrModelParts) {
        block_for_each(r_model_part.get().Elements(), [](Element& rElement) {
            const auto& r_properties = rElement.GetProperties();
            CheckK0MainDirection(r_properties, rElement.Id());
            CheckSufficientMaterialParameters(r_properties, rElement.Id());
            CheckOCRorPOP(r_properties, rElement.Id());
            CheckPoissonUnloadingReloading(r_properties, rElement.Id());
            CheckPhi(r_properties, rElement.Id());
            CheckK0(r_properties, rElement.Id());
        });
    }

    return 0;
}

void ApplyK0ProcedureProcess::CheckK0MainDirection(const Properties& rProperties, IndexType ElementId)
{
    KRATOS_ERROR_IF(!rProperties.Has(K0_MAIN_DIRECTION))
        << "K0_MAIN_DIRECTION is not defined for element " << ElementId << "." << std::endl;

    const auto dimension = rProperties.GetValue(CONSTITUTIVE_LAW).get()->WorkingSpaceDimension();
    if (dimension == 2) {
        KRATOS_ERROR_IF(rProperties[K0_MAIN_DIRECTION] < 0 || rProperties[K0_MAIN_DIRECTION] > 1)
            << "K0_MAIN_DIRECTION should be 0 or 1 for element " << ElementId << "." << std::endl;
    } else if (dimension == 3) {
        KRATOS_ERROR_IF(rProperties[K0_MAIN_DIRECTION] < 0 || rProperties[K0_MAIN_DIRECTION] > 2)
            << "K0_MAIN_DIRECTION should be 0, 1 or 2 for element " << ElementId << "." << std::endl;
    } else {
        KRATOS_ERROR << "dimension should be 2 or 3 for element " << ElementId << "." << std::endl;
    }
}

void ApplyK0ProcedureProcess::CheckK0(const Properties& rProperties, IndexType ElementId)
{
    const auto k0_variables = Geo::ConstVariableRefs{
        std::cref(K0_NC), std::cref(K0_VALUE_XX), std::cref(K0_VALUE_YY), std::cref(K0_VALUE_ZZ)};
    for (const auto& k0_variable_ref : k0_variables) {
        if (const auto& r_k0_variable = k0_variable_ref.get(); rProperties.Has(r_k0_variable)) {
            const double k0_var_value = rProperties[r_k0_variable];
            KRATOS_ERROR_IF(k0_var_value < 0.0)
                << r_k0_variable.Name() << " (" << k0_var_value
                << ") should be in the range [0.0,-> for element " << ElementId << "." << std::endl;
        }
    }
}

void ApplyK0ProcedureProcess::CheckPhi(const Properties& rProperties, IndexType ElementId)
{
    if (rProperties.Has(INDEX_OF_UMAT_PHI_PARAMETER) && rProperties.Has(UMAT_PARAMETERS)) {
        const auto phi_index = rProperties[INDEX_OF_UMAT_PHI_PARAMETER];
        const auto number_of_umat_parameters = static_cast<int>(rProperties[UMAT_PARAMETERS].size());

        KRATOS_ERROR_IF(phi_index < 1 || phi_index > number_of_umat_parameters)
            << "INDEX_OF_UMAT_PHI_PARAMETER (" << phi_index
            << ") is not in range 1, size of UMAT_PARAMETERS for element " << ElementId << "." << std::endl;

        const double phi = rProperties[UMAT_PARAMETERS][phi_index - 1];
        KRATOS_ERROR_IF(phi < 0.0 || phi > 90.0)
            << "Phi (" << phi << ") should be between 0 and 90 degrees for element " << ElementId
            << "." << std::endl;
    }
}

void ApplyK0ProcedureProcess::CheckOCRorPOP(const Properties& rProperties, IndexType ElementId)
{
    if (rProperties.Has(K0_NC) ||
        (rProperties.Has(INDEX_OF_UMAT_PHI_PARAMETER) && rProperties.Has(UMAT_PARAMETERS))) {
        if (rProperties.Has(OCR)) {
            const auto ocr = rProperties[OCR];
            KRATOS_ERROR_IF(ocr < 1.0) << "OCR (" << ocr << ") should be in the range [1.0,-> for element "
                                       << ElementId << "." << std::endl;
        }

        if (rProperties.Has(POP)) {
            const auto pop = rProperties[POP];
            KRATOS_ERROR_IF(pop < 0.0) << "POP (" << pop << ") should be in the range [0.0,-> for element "
                                       << ElementId << "." << std::endl;
        }
    }
}

void ApplyK0ProcedureProcess::CheckPoissonUnloadingReloading(const Properties& rProperties, IndexType ElementId)
{
    KRATOS_ERROR_IF(rProperties.Has(POISSON_UNLOADING_RELOADING) &&
                    (rProperties[POISSON_UNLOADING_RELOADING] < -1.0 ||
                     rProperties[POISSON_UNLOADING_RELOADING] >= 0.5))
        << "POISSON_UNLOADING_RELOADING (" << rProperties[POISSON_UNLOADING_RELOADING]
        << ") is not in range [-1.0, 0.5> for element " << ElementId << "." << std::endl;

    if (rProperties.Has(K0_VALUE_XX)) {
        KRATOS_ERROR_IF(rProperties.Has(POISSON_UNLOADING_RELOADING) || rProperties.Has(OCR) ||
                        rProperties.Has(POP))
            << "Insufficient material data for K0 procedure process for element "
            << ElementId << ". Poisson unloading-reloading, OCR and POP functionality cannot be combined with K0_VALUE_XX, _YY and _ZZ."
            << std::endl;
    }
}

void ApplyK0ProcedureProcess::CheckSufficientMaterialParameters(const Properties& rProperties, IndexType ElementId)
{
    KRATOS_ERROR_IF_NOT(
        rProperties.Has(K0_NC) ||
        (rProperties.Has(INDEX_OF_UMAT_PHI_PARAMETER) && rProperties.Has(UMAT_PARAMETERS)) ||
        (rProperties.Has(K0_VALUE_XX) && rProperties.Has(K0_VALUE_YY) && rProperties.Has(K0_VALUE_ZZ)))
        << "Insufficient material data for K0 procedure process for element " << ElementId << ". No K0_NC, "
        << "(INDEX_OF_UMAT_PHI_PARAMETER and UMAT_PARAMETERS) or (K0_VALUE_XX, _YY and _ZZ found)."
        << std::endl;
}

void ApplyK0ProcedureProcess::ExecuteFinalizeSolutionStep()
{
    KRATOS_TRY

    // K0 procedure for the model part:
    for (const auto& r_model_part : mrModelParts) {
        block_for_each(r_model_part.get().Elements(), [&r_model_part](Element& rElement) {
            CalculateK0Stresses(rElement, r_model_part.get().GetProcessInfo());
        });
    }
    KRATOS_CATCH("")
}

std::string ApplyK0ProcedureProcess::Info() const { return "ApplyK0ProcedureProcess"s; }

bool ApplyK0ProcedureProcess::UseStandardProcedure() const
{
    const auto setting_name = std::string{"use_standard_procedure"};
    return !mSettings.Has(setting_name) || mSettings[setting_name].GetBool();
}

array_1d<double, 3> ApplyK0ProcedureProcess::CreateK0Vector(const Element::PropertiesType& rProp)
{
    // Check for alternative K0 specifications
    array_1d<double, 3> k0_vector;
    if (rProp.Has(K0_NC)) {
        std::fill(k0_vector.begin(), k0_vector.end(), rProp[K0_NC]);
    } else if (rProp.Has(INDEX_OF_UMAT_PHI_PARAMETER) && rProp.Has(UMAT_PARAMETERS)) {
        const auto phi_in_radians = ConstitutiveLawUtilities::GetFrictionAngleInRadians(rProp);
        std::fill(k0_vector.begin(), k0_vector.end(), 1.0 - std::sin(phi_in_radians));
    } else {
        k0_vector[0] = rProp[K0_VALUE_XX];
        k0_vector[1] = rProp[K0_VALUE_YY];
        k0_vector[2] = rProp[K0_VALUE_ZZ];
    }

    return k0_vector;
}

void ApplyK0ProcedureProcess::CalculateK0Stresses(Element& rElement, const ProcessInfo& rProcessInfo)
{
    // Get K0 material parameters of this element ( probably there is something more efficient )
    const Element::PropertiesType& rProp             = rElement.GetProperties();
    const auto                     k0_main_direction = rProp[K0_MAIN_DIRECTION];

    auto k0_vector = CreateK0Vector(rProp);

    // Corrections on k0_vector by OCR or POP
    const auto PoissonUR = rProp.Has(POISSON_UNLOADING_RELOADING) ? rProp[POISSON_UNLOADING_RELOADING] : 0.;
    const auto PoissonURfactor = PoissonUR / (1. - PoissonUR);

    double POP_value = 0.0;
    if (rProp.Has(K0_NC) || rProp.Has(INDEX_OF_UMAT_PHI_PARAMETER)) {
        if (rProp.Has(OCR)) {
            // Determine OCR dependent K0 values ( constant per element! )
            k0_vector *= rProp[OCR];
            const array_1d<double, 3> correction(3, PoissonURfactor * (rProp[OCR] - 1.0));
            k0_vector -= correction;
        } else if (rProp.Has(POP)) {
            // POP is entered as positive value, convention here is compression negative.
            POP_value = -rProp[POP];
        }
    }
    // Get element stress vectors
    std::vector<ConstitutiveLaw::StressVectorType> rStressVectors;
    rElement.CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, rStressVectors, rProcessInfo);

    // Loop over integration point stress vectors
    for (auto& rStressVector : rStressVectors) {
        // Apply K0 procedure
        for (int i_dir = 0; i_dir <= 2; ++i_dir) {
            if (i_dir != k0_main_direction) {
                rStressVector[i_dir] = k0_vector[i_dir] * (rStressVector[k0_main_direction] + POP_value) -
                                       PoissonURfactor * POP_value;
            }
        }
        // Erase shear stresses
        std::fill(rStressVector.begin() + 3, rStressVector.end(), 0.0);
    }
    // Set element integration point stress tensors
    rElement.SetValuesOnIntegrationPoints(CAUCHY_STRESS_VECTOR, rStressVectors, rProcessInfo);
}

} // namespace Kratos
