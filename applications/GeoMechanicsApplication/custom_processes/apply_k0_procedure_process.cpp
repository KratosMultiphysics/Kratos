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

// System includes
#include <cmath>
#include <iostream>

// Project includes
#include "includes/element.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"

// Application includes
#include "custom_constitutive/linear_elastic_law.h"
#include "geo_mechanics_application_variables.h"
#include "includes/kratos_flags.h"

namespace
{

using namespace Kratos;

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

ApplyK0ProcedureProcess::ApplyK0ProcedureProcess(ModelPart& model_part, const Parameters& rK0Settings)
    : Process(Flags()), mrModelPart(model_part), mSettings(rK0Settings)
{
}

void ApplyK0ProcedureProcess::ExecuteInitialize()
{
    if (UseStandardProcedure())
        SetConsiderDiagonalEntriesOnlyAndNoShear(mrModelPart.Elements(), true);
}

void ApplyK0ProcedureProcess::ExecuteFinalize()
{
    if (UseStandardProcedure())
        SetConsiderDiagonalEntriesOnlyAndNoShear(mrModelPart.Elements(), false);
}

int ApplyK0ProcedureProcess::Check()
{
    block_for_each(mrModelPart.Elements(), [this](Element& rElement) {
        auto r_properties = rElement.GetProperties();
        KRATOS_ERROR_IF(!r_properties.Has(K0_MAIN_DIRECTION))
            << "K0_MAIN_DIRECTION is not defined for element " << rElement.Id() << "." << std::endl;
        KRATOS_ERROR_IF(r_properties[K0_MAIN_DIRECTION] < 0 || r_properties[K0_MAIN_DIRECTION] > 1)
            << "K0_MAIN_DIRECTION should be 0 or 1 for element " << rElement.Id() << "." << std::endl;

        KRATOS_ERROR_IF_NOT(
            r_properties.Has(K0_NC) ||
            (r_properties.Has(INDEX_OF_UMAT_PHI_PARAMETER) &&
             r_properties.Has(NUMBER_OF_UMAT_PARAMETERS) && r_properties.Has(UMAT_PARAMETERS)) ||
            (r_properties.Has(K0_VALUE_XX) && r_properties.Has(K0_VALUE_YY) && r_properties.Has(K0_VALUE_ZZ)))
            << "Insufficient material data for K0 procedure process for element " << rElement.Id() << ". No K0_NC, "
            << "(INDEX_OF_UMAT_PHI_PARAMETER, NUMBER_OF_UMAT_PARAMETERS and "
               "UMAT_PARAMETERS) or (K0_VALUE_XX, _YY and _ZZ found)."
            << std::endl;

        KRATOS_ERROR_IF(r_properties.Has(POISSON_UNLOADING_RELOADING) &&
                        (r_properties[POISSON_UNLOADING_RELOADING] < -1.0 ||
                         r_properties[POISSON_UNLOADING_RELOADING] >= 0.5))
            << "POISSON_UNLOADING_RELOADING (" << r_properties[POISSON_UNLOADING_RELOADING]
            << ") is not in range [-1.0, 0.5> for element " << rElement.Id() << "." << std::endl;

        if (r_properties.Has(K0_VALUE_XX)) {
            KRATOS_ERROR_IF(r_properties.Has(POISSON_UNLOADING_RELOADING) ||
                            r_properties.Has(OCR) || r_properties.Has(POP))
                << "Insufficient material data for K0 procedure process for element "
                << rElement.Id() << ". Poisson unloading-reloading, OCR and POP functionality cannot be combined with K0_VALUE_XX, _YY and _ZZ."
                << std::endl;
        }

        if (r_properties.Has(INDEX_OF_UMAT_PHI_PARAMETER) &&
            r_properties.Has(NUMBER_OF_UMAT_PARAMETERS) && r_properties.Has(UMAT_PARAMETERS)) {
            const auto phi_index                 = r_properties[INDEX_OF_UMAT_PHI_PARAMETER];
            const auto number_of_umat_parameters = r_properties[NUMBER_OF_UMAT_PARAMETERS];

            KRATOS_ERROR_IF(phi_index < 1 || phi_index > number_of_umat_parameters)
                << "INDEX_OF_UMAT_PHI_PARAMETER (" << phi_index << ") is not in range 1, NUMBER_OF_UMAT_PARAMETERS ("
                << number_of_umat_parameters << ") for element " << rElement.Id() << "." << std::endl;

            const double phi = r_properties[UMAT_PARAMETERS][phi_index - 1];
            KRATOS_ERROR_IF(phi < 0.0 || phi > 90.0)
                << "Phi should be between 0 and 90 degrees for element " << rElement.Id() << "."
                << std::endl;
        }
    });

    return 0;
}

void ApplyK0ProcedureProcess::ExecuteFinalizeSolutionStep()
{
    KRATOS_TRY

    // K0 procedure for the model part:
    block_for_each(mrModelPart.Elements(),
                   [this](Element& rElement) { CalculateK0Stresses(rElement); });

    KRATOS_CATCH("")
}

std::string ApplyK0ProcedureProcess::Info() const { return "ApplyK0ProcedureProcess"; }

bool ApplyK0ProcedureProcess::UseStandardProcedure() const
{
    const auto setting_name = std::string{"use_standard_procedure"};
    return !mSettings.Has(setting_name) || mSettings[setting_name].GetBool();
}

array_1d<double, 3> ApplyK0ProcedureProcess::CreateK0Vector(const Element::PropertiesType& rProp) const
{
    // Check for alternative K0 specifications
    array_1d<double, 3> k0_vector;
    if (rProp.Has(K0_NC)) {
        std::fill(k0_vector.begin(), k0_vector.end(), rProp[K0_NC]);
    } else if (rProp.Has(INDEX_OF_UMAT_PHI_PARAMETER) && rProp.Has(NUMBER_OF_UMAT_PARAMETERS) &&
               rProp.Has(UMAT_PARAMETERS)) {
        const auto phi = rProp[UMAT_PARAMETERS][rProp[INDEX_OF_UMAT_PHI_PARAMETER] - 1];
        std::fill(k0_vector.begin(), k0_vector.end(), 1.0 - std::sin(MathUtils<>::DegreesToRadians(phi)));
    } else {
        k0_vector[0] = rProp[K0_VALUE_XX];
        k0_vector[1] = rProp[K0_VALUE_YY];
        k0_vector[2] = rProp[K0_VALUE_ZZ];
    }

    return k0_vector;
}

void ApplyK0ProcedureProcess::CalculateK0Stresses(Element& rElement)
{
    // Get K0 material parameters of this element ( probably there is something more efficient )
    const Element::PropertiesType& rProp             = rElement.GetProperties();
    const int                      k0_main_direction = rProp[K0_MAIN_DIRECTION];

    auto k0_vector = CreateK0Vector(rProp);

    // Corrections on k0_vector by OCR or POP
    const auto PoissonUR = rProp.Has(POISSON_UNLOADING_RELOADING) ? rProp[POISSON_UNLOADING_RELOADING] : 0.;
    const auto PoissonURfactor = PoissonUR / (1. - PoissonUR);

    double POP_value = 0.0;
    if (rProp.Has(K0_NC) || rProp.Has(INDEX_OF_UMAT_PHI_PARAMETER)) {
        if (rProp.Has(OCR)) {
            // Determine OCR dependent K0 values ( constant per element! )
            k0_vector *= rProp[OCR];
            array_1d<double, 3> correction(3, PoissonURfactor * (rProp[OCR] - 1.0));
            k0_vector -= correction;
        } else if (rProp.Has(POP)) {
            // POP is entered as positive value, convention here is compression negative.
            POP_value = -rProp[POP];
        }
    }
    // Get element stress vectors
    const ProcessInfo& rCurrentProcessInfo = this->mrModelPart.GetProcessInfo();
    std::vector<ConstitutiveLaw::StressVectorType> rStressVectors;
    rElement.CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, rStressVectors, rCurrentProcessInfo);

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
    rElement.SetValuesOnIntegrationPoints(CAUCHY_STRESS_VECTOR, rStressVectors, rCurrentProcessInfo);
}

} // namespace Kratos
