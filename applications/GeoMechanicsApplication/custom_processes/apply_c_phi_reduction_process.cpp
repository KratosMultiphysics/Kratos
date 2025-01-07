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
#include "includes/model_part.h"
#include "utilities/math_utils.h"

namespace Kratos
{

void ApplyCPhiReductionProcess::ExecuteInitializeSolutionStep()
{
    KRATOS_TRY

    if (IsStepRestarted()) mReductionIncrement *= 0.5;
    mReductionFactor = mPreviousReductionFactor - mReductionIncrement;
    KRATOS_ERROR_IF(mReductionFactor <= 0.01)
        << "Reduction factor should not drop below 0.01, calculation stopped." << std::endl;
    KRATOS_ERROR_IF(mReductionIncrement <= 0.001)
    << "Reduction increment should not drop below 0.001, calculation stopped." << std::endl;
    KRATOS_INFO("ApplyCPhiReductionProces::ExecuteInitializeSolutionStep")
        << "Try a c-phi reduction factor " << mReductionFactor << " (safety factor "
        << 1. / mReductionFactor << ") Previous reduction = " << mPreviousReductionFactor
        << " Reduction increment = " << mReductionIncrement << std::endl;

    double    phi                = 0.;
    double    reduced_phi        = 0.;
    double    c                  = 0.;
    double    reduced_c          = 0.;
    IndexType previousPropertyId = -1;
    // Apply C/Phi Reduction procedure for the model part:
    block_for_each(mrModelPart.Elements(),
                   [this, &phi, &reduced_phi, &c, &reduced_c, &previousPropertyId](Element& rElement) {
        // Only compute new c and phi if the Id changes
        if (rElement.GetProperties().Id() != previousPropertyId) {
            phi                = GetAndCheckPhi(rElement.GetProperties());
            reduced_phi        = ComputeReducedPhi(phi);
            c                  = GetAndCheckC(rElement.GetProperties());
            reduced_c          = mReductionFactor * c;
            previousPropertyId = rElement.GetProperties().Id();
        }
        SetCPhiAtElement(rElement, reduced_phi, reduced_c);
    });
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
    KRATOS_ERROR_IF(mrModelPart.Elements().empty())
        << "ApplyCPhiReductionProces has no elements in modelpart " << mrModelPart.Name() << std::endl;
    return 0;
}

double ApplyCPhiReductionProcess::GetAndCheckPhi(const Element::PropertiesType& rProp) const
{
    // Get the initial properties from the model part. Recall that we create a separate
    // properties object with reduced c and phi for each and every element. Those reduced
    // properties objects are not linked to the original ones.
    const auto& part_properties = mrModelPart.GetProperties(rProp.Id());

    KRATOS_ERROR_IF_NOT(part_properties.Has(UMAT_PARAMETERS))
        << "Missing required item UMAT_PARAMETERS" << std::endl;
    KRATOS_ERROR_IF_NOT(part_properties.Has(INDEX_OF_UMAT_PHI_PARAMETER))
        << "Missing required item INDEX_OF_UMAT_PHI_PARAMETER" << std::endl;

    KRATOS_ERROR_IF(part_properties[INDEX_OF_UMAT_PHI_PARAMETER] < 1 ||
                    part_properties[INDEX_OF_UMAT_PHI_PARAMETER] >
                        static_cast<int>(part_properties[UMAT_PARAMETERS].size()))
        << "Invalid INDEX_OF_UMAT_PHI_PARAMETER: " << part_properties[INDEX_OF_UMAT_PHI_PARAMETER]
        << " (out-of-bounds index)" << std::endl;
    const double phi = part_properties[UMAT_PARAMETERS][part_properties[INDEX_OF_UMAT_PHI_PARAMETER] - 1];
    KRATOS_ERROR_IF(phi < 0. || phi > 90.) << "Friction angle Phi out of range: " << phi << std::endl;
    return phi;
}

double ApplyCPhiReductionProcess::ComputeReducedPhi(double Phi) const
{
    // Phi converted to radians and then its tangent is reduced by the reduction factor
    const double tan_phi         = std::tan(MathUtils<>::DegreesToRadians(Phi));
    const double reduced_tan_phi = tan_phi * mReductionFactor;
    return std::atan(reduced_tan_phi) * 180.0 / Globals::Pi;
}

double ApplyCPhiReductionProcess::GetAndCheckC(const Element::PropertiesType& rProp) const
{
    // Get the initial properties from the model part. Recall that we create a separate
    // properties object with reduced c and phi for each and every element. Those reduced
    // properties objects are not linked to the original ones.
    const auto& part_properties = mrModelPart.GetProperties(rProp.Id());

    KRATOS_ERROR_IF_NOT(part_properties.Has(UMAT_PARAMETERS))
        << "Missing required item UMAT_PARAMETERS" << std::endl;
    KRATOS_ERROR_IF_NOT(part_properties.Has(INDEX_OF_UMAT_C_PARAMETER))
        << "Missing required item INDEX_OF_UMAT_C_PARAMETER" << std::endl;

    KRATOS_ERROR_IF(part_properties[INDEX_OF_UMAT_C_PARAMETER] < 1 ||
                    part_properties[INDEX_OF_UMAT_C_PARAMETER] >
                        static_cast<int>(part_properties[UMAT_PARAMETERS].size()))
        << "invalid INDEX_OF_UMAT_C_PARAMETER: " << part_properties[INDEX_OF_UMAT_C_PARAMETER]
        << " (out-of-bounds index)" << std::endl;
    const auto c = part_properties[UMAT_PARAMETERS][part_properties[INDEX_OF_UMAT_C_PARAMETER] - 1];
    KRATOS_ERROR_IF(c < 0.) << "Cohesion C out of range: " << c << std::endl;
    return c;
}

void ApplyCPhiReductionProcess::SetCPhiAtElement(Element& rElement, double ReducedPhi, double ReducedC) const
{
    // Get C/Phi material properties of this element
    const auto& r_prop = rElement.GetProperties();

    // Overwrite C and Phi in the UMAT_PARAMETERS
    auto Umat_parameters                                     = r_prop[UMAT_PARAMETERS];
    Umat_parameters[r_prop[INDEX_OF_UMAT_PHI_PARAMETER] - 1] = ReducedPhi;
    Umat_parameters[r_prop[INDEX_OF_UMAT_C_PARAMETER] - 1]   = ReducedC;

    SetValueAtElement(rElement, UMAT_PARAMETERS, Umat_parameters);
}

void ApplyCPhiReductionProcess::SetValueAtElement(Element&                rElement,
                                                  const Variable<Vector>& rVariable,
                                                  const Vector&           rValue) const
{
    // Copies properties
    Properties::Pointer p_new_prop = Kratos::make_shared<Properties>(rElement.GetProperties());

    // Adds new properties to the element
    p_new_prop->SetValue(rVariable, rValue);
    rElement.SetProperties(p_new_prop);
}

bool ApplyCPhiReductionProcess::IsStepRestarted() const
{
    return mrModelPart.GetProcessInfo().GetValue(NUMBER_OF_CYCLES) > 1;
}
} // namespace Kratos