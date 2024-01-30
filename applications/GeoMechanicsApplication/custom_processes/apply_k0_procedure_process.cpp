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

void ApplyK0ProcedureProcess::CalculateK0Stresses(Element& rElement)
{
    // Get K0 material parameters of this element ( probably there is something more efficient )
    const Element::PropertiesType& rProp = rElement.GetProperties();
    const int k0_main_direction          = rProp[K0_MAIN_DIRECTION];
    if (k0_main_direction < 0 || k0_main_direction > 1) {
        KRATOS_ERROR << "undefined K0_MAIN_DIRECTION in ApplyK0ProcedureProcess: " << k0_main_direction
                     << std::endl;
    }

    // Check for alternative K0 specifications
    array_1d<double, 3> k0_vector;
    if (rProp.Has(K0_NC)) {
        std::fill(k0_vector.begin(), k0_vector.end(), rProp[K0_NC]);
    } else if (rProp.Has(INDEX_OF_UMAT_PHI_PARAMETER) && rProp.Has(NUMBER_OF_UMAT_PARAMETERS) &&
               rProp.Has(UMAT_PARAMETERS)) {
        if (rProp[INDEX_OF_UMAT_PHI_PARAMETER] < 1 ||
            rProp[INDEX_OF_UMAT_PHI_PARAMETER] > rProp[NUMBER_OF_UMAT_PARAMETERS]) {
            KRATOS_ERROR << "undefined INDEX_OF_UMAT_PHI_PARAMETER in ApplyK0ProcedureProcess: "
                         << rProp[INDEX_OF_UMAT_PHI_PARAMETER] << std::endl;
        }
        // is more checking is possible and should that happen here?
        double phi = rProp[UMAT_PARAMETERS][rProp[INDEX_OF_UMAT_PHI_PARAMETER] - 1];
        if (phi < 0. || phi > 90.) {
            KRATOS_ERROR << "friction angle Phi out of range in ApplyK0ProcedureProcess: " << phi << std::endl;
        }
        std::fill(k0_vector.begin(), k0_vector.end(), 1.0 - std::sin(MathUtils<>::DegreesToRadians(phi)));
    } else if (rProp.Has(K0_VALUE_XX) && rProp.Has(K0_VALUE_YY) && rProp.Has(K0_VALUE_ZZ)) {
        k0_vector[0] = rProp[K0_VALUE_XX];
        k0_vector[1] = rProp[K0_VALUE_YY];
        k0_vector[2] = rProp[K0_VALUE_ZZ];
    } else {
        KRATOS_ERROR << "Insufficient material data for K0 procedure process: " << std::endl;
    }
    const auto PoissonUR = rProp.Has(POISSON_UNLOADING_RELOADING) ? rProp[POISSON_UNLOADING_RELOADING] : 0.;

    // Determine OCR dependent K0 values ( constant per element! )
    if ((rProp.Has(K0_NC) || rProp.Has(INDEX_OF_UMAT_PHI_PARAMETER)) && rProp.Has(OCR)) {
        // Modify for presence of OCR (or POP?) field values
        k0_vector *= rProp[OCR];
        array_1d<double, 3> correction(3, (PoissonUR / (1.0 - PoissonUR)) * (rProp[OCR] - 1.0));
        k0_vector -= correction;
    }
    // Get element stress vectorsn
    const ProcessInfo& rCurrentProcessInfo = this->mrModelPart.GetProcessInfo();
    std::vector<ConstitutiveLaw::StressVectorType> rStressVectors;
    rElement.CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, rStressVectors, rCurrentProcessInfo);

    // Loop over integration points
    const Element::GeometryType& rGeom = rElement.GetGeometry();
    const Element::GeometryType::IntegrationPointsArrayType& integration_points =
        rGeom.IntegrationPoints(rElement.GetIntegrationMethod());
    for (unsigned int g_point = 0; g_point < integration_points.size(); ++g_point) {
        // Apply K0 procedure
        for (int i_dir = 0; i_dir <= 2; ++i_dir) {
            if (i_dir != k0_main_direction) {
                rStressVectors[g_point][i_dir] = k0_vector[i_dir] * rStressVectors[g_point][k0_main_direction];
            }
        }
        // Erase shear stresses
        std::fill(rStressVectors[g_point].begin() + 3, rStressVectors[g_point].end(), 0.0);
    }
    // Set element integration point stress tensors
    rElement.SetValuesOnIntegrationPoints(CAUCHY_STRESS_VECTOR, rStressVectors, rCurrentProcessInfo);
}

} // namespace Kratos
