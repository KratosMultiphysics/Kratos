// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors: Marjan Fathian
//                Wijtze Pieter Kikstra
//                Anne van de Graaf
//

// Project includes
#include "custom_processes/apply_c_phi_reduction_process.h"
#include "containers/model.h"
#include "custom_constitutive/constitutive_law_dimension.h"
#include "custom_constitutive/mohr_coulomb_with_tension_cutoff.h"
#include "custom_utilities/check_utilities.hpp"
#include "custom_utilities/constitutive_law_utilities.h"
#include "custom_utilities/process_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"

namespace Kratos
{
using namespace std::string_literals;

namespace
{

bool UsesUmatParameters(const Properties& rProperties)
{
    KRATOS_ERROR_IF_NOT(rProperties.Has(CONSTITUTIVE_LAW))
        << "Properties do not have CONSTITUTIVE_LAW" << std::endl;

    const auto law_name = rProperties[CONSTITUTIVE_LAW]->Info();
    return law_name.find("UMAT") != std::string::npos || law_name.find("UDSM") != std::string::npos;
}

bool UsesInternalMohrCoulombModel(const Element& rElement)
{
    KRATOS_ERROR_IF_NOT(rElement.GetProperties().Has(CONSTITUTIVE_LAW))
        << "Properties do not have CONSTITUTIVE_LAW" << std::endl;

    return dynamic_cast<const MohrCoulombWithTensionCutOff*>(
               rElement.GetProperties()[CONSTITUTIVE_LAW].get()) != nullptr;
}
} // namespace

ApplyCPhiReductionProcess::ApplyCPhiReductionProcess(Model& rModel, const Parameters& rProcessSettings)
{
    mrModelParts = ProcessUtilities::GetModelPartsFromSettings(rModel, rProcessSettings,
                                                               ApplyCPhiReductionProcess::Info());
}

void ApplyCPhiReductionProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    if (IsStepRestarted()) mReductionIncrement *= 0.5;
    mReductionFactor = mPreviousReductionFactor - mReductionIncrement;
    KRATOS_ERROR_IF(mReductionFactor <= 0.01)
        << "Reduction factor should not drop below 0.01, calculation stopped." << std::endl;
    KRATOS_ERROR_IF(mReductionIncrement <= 0.001)
        << "Reduction increment should not drop below 0.001, calculation stopped. Final safety "
           "factor = "
        << 1.0 / mPreviousReductionFactor << std::endl;
    KRATOS_INFO("ApplyCPhiReductionProces::ExecuteInitializeSolutionStep")
        << "Try a c-phi reduction factor " << mReductionFactor << " (safety factor "
        << 1. / mReductionFactor << ") Previous reduction = " << mPreviousReductionFactor
        << " Reduction increment = " << mReductionIncrement << std::endl;

    double    phi                  = 0.;
    double    reduced_phi          = 0.;
    double    c                    = 0.;
    double    reduced_c            = 0.;
    IndexType previous_property_Id = -1;

    for (const auto& r_model_part : mrModelParts) {
        // Apply C/Phi Reduction procedure for the model part:
        block_for_each(r_model_part.get().Elements(), [this, &r_model_part, &phi, &reduced_phi, &c, &reduced_c,
                                                       &previous_property_Id](Element& rElement) {
            // Only compute new c and phi if the Id changes
            if (const auto element_property_Id = rElement.GetProperties().Id();
                element_property_Id != previous_property_Id) {
                const auto& r_part_properties = r_model_part.get().GetProperties(element_property_Id);
                phi                  = GetAndCheckPhi(r_part_properties, element_property_Id);
                reduced_phi          = ComputeReducedPhi(phi);
                c                    = GetAndCheckC(r_part_properties);
                reduced_c            = mReductionFactor * c;
                previous_property_Id = element_property_Id;
            }
            SetCPhiAtElement(rElement, reduced_phi, reduced_c);
        });
    }
    KRATOS_CATCH("")
}

void ApplyCPhiReductionProcess::ExecuteFinalizeSolutionStep()
{
    mPreviousReductionFactor = mReductionFactor;
    KRATOS_INFO("ApplyCPhiReductionProcess::ExecuteFinalizeSolutionStep")
        << "Reached safety factor so far = " << 1.0 / mReductionFactor << std::endl;
}

void ApplyCPhiReductionProcess::ExecuteFinalize()
{
    KRATOS_INFO("ApplyCPhiReductionProcess") << "Final safety factor = " << 1.0 / mReductionFactor << std::endl;
}

int ApplyCPhiReductionProcess::Check()
{
    KRATOS_ERROR_IF(std::ranges::all_of(mrModelParts, [](const auto& r_model_part) {
        return r_model_part.get().Elements().empty();
    })) << "None of the provided model parts contains at least one element. A c-phi reduction analysis requires at least one element.\n";

    for (const auto& r_model_part : mrModelParts) {
        for (const auto& r_element : r_model_part.get().Elements()) {
            auto                  r_properties = r_element.GetProperties();
            const CheckProperties check_properties(r_properties, "model part property",
                                                   CheckProperties::Bounds::AllInclusive);
            if (UsesUmatParameters(r_properties)) {
                check_properties.CheckAvailability(UMAT_PARAMETERS);
                check_properties.Check(INDEX_OF_UMAT_PHI_PARAMETER, 1,
                                       static_cast<int>(r_properties[UMAT_PARAMETERS].size()));
                check_properties.Check(INDEX_OF_UMAT_C_PARAMETER, 1,
                                       static_cast<int>(r_properties[UMAT_PARAMETERS].size()));
            } else {
                check_properties.Check(GEO_COHESION);
                check_properties.Check(GEO_FRICTION_ANGLE);
            }
        }
    }
    return 0;
}

double ApplyCPhiReductionProcess::GetAndCheckPhi(const Properties& rModelPartProperties, IndexType ElementPropertyId)
{
    // Get the initial properties from the model part. Recall that we create a separate
    // properties object with reduced c and phi for each and every element. Those reduced
    // properties objects are not linked to the original ones.

    const auto phi = ConstitutiveLawUtilities::GetFrictionAngleInDegrees(rModelPartProperties);
    KRATOS_ERROR_IF(phi < 0. || phi > 90.)
        << "Friction angle Phi in the model part property with Id " << ElementPropertyId
        << " has an invalid value: " << phi << " is out of range [0,90] (degrees)." << std::endl;
    return phi;
}

double ApplyCPhiReductionProcess::ComputeReducedPhi(double Phi) const
{
    // Phi converted to radians and then its tangent is reduced by the reduction factor
    const double tan_phi         = std::tan(MathUtils<>::DegreesToRadians(Phi));
    const double reduced_tan_phi = tan_phi * mReductionFactor;
    return std::atan(reduced_tan_phi) * 180.0 / Globals::Pi;
}

double ApplyCPhiReductionProcess::GetAndCheckC(const Properties& rModelPartProperties)
{
    // Get the initial properties from the model part. Recall that we create a separate
    // properties object with reduced c and phi for each and every element. Those reduced
    // properties objects are not linked to the original ones.

    const auto c = ConstitutiveLawUtilities::GetCohesion(rModelPartProperties);
    KRATOS_ERROR_IF(c < 0.) << "Cohesion C out of range: " << c << std::endl;
    return c;
}

void ApplyCPhiReductionProcess::SetCPhiAtElement(Element& rElement, double ReducedPhi, double ReducedC)
{
    const auto&         r_properties     = rElement.GetProperties();
    Properties::Pointer p_new_properties = Kratos::make_shared<Properties>(r_properties);

    if (UsesUmatParameters(r_properties)) {
        // Overwrite C and Phi in the UMAT_PARAMETERS
        auto umat_parameters = r_properties[UMAT_PARAMETERS];
        umat_parameters[r_properties[INDEX_OF_UMAT_PHI_PARAMETER] - 1] = ReducedPhi;
        umat_parameters[r_properties[INDEX_OF_UMAT_C_PARAMETER] - 1]   = ReducedC;
        p_new_properties->SetValue(UMAT_PARAMETERS, umat_parameters);
    } else {
        p_new_properties->SetValue(GEO_FRICTION_ANGLE, ReducedPhi);
        p_new_properties->SetValue(GEO_COHESION, ReducedC);
    }
    rElement.SetProperties(p_new_properties);

    if (UsesInternalMohrCoulombModel(rElement))
        InitializeParametersForInternalMohrCoulombModel(rElement);
}

void ApplyCPhiReductionProcess::InitializeParametersForInternalMohrCoulombModel(Element& rElement)
{
    KRATOS_TRY
    // Due to current implementation of the internal Mohr-Coulomb model, this function re-initializes material properties in it.
    std::vector<ConstitutiveLaw::Pointer> constitutive_laws;
    const ProcessInfo                     dummy_process_info;
    rElement.CalculateOnIntegrationPoints(CONSTITUTIVE_LAW, constitutive_laws, dummy_process_info);
    const auto& r_properties   = rElement.GetProperties();
    const auto  dummy_geometry = Geometry<Node>{};
    const auto  dummy_vector   = Vector();
    for (auto& p_law : constitutive_laws) {
        if (auto p_mohr_coulomb = dynamic_cast<MohrCoulombWithTensionCutOff*>(p_law.get())) {
            p_mohr_coulomb->InitializeMaterial(r_properties, dummy_geometry, dummy_vector);
        }
    }
    KRATOS_CATCH("")
}

bool ApplyCPhiReductionProcess::IsStepRestarted() const
{
    return mrModelParts[0].get().GetProcessInfo().GetValue(NUMBER_OF_CYCLES) > 1;
}

std::string ApplyCPhiReductionProcess::Info() const { return "ApplyCPhiReductionProcess"s; }

} // namespace Kratos