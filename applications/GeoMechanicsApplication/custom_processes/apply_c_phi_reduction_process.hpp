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

#pragma once

// System includes
#include <cmath>
#include <iostream>

// Project includes
#include "custom_elements/small_strain_U_Pw_diff_order_element.hpp"
#include "custom_utilities/stress_strain_utilities.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "utilities/math_utils.h"

// Application includes
#include "geo_mechanics_application_variables.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{

class ApplyCPhiReductionProcess : public Process
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(ApplyCPhiReductionProcess);

    ApplyCPhiReductionProcess(ModelPart& rModelPart, const Parameters&)
        : Process(Flags()), mrModelPart(rModelPart)
    {
    }

    ~ApplyCPhiReductionProcess() override = default;

    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY

        // Recognize whether a step is restarted. If so, try again with a smaller reduction increment.
        if (mrModelPart.GetProcessInfo().GetValue(NUMBER_OF_CYCLES) > 1) mReductionIncrement *= 0.5;
        mReductionFactor = mPreviousReductionFactor - mReductionIncrement;
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

    void ExecuteFinalizeSolutionStep() override
    {
        mPreviousReductionFactor = mReductionFactor;
        KRATOS_INFO("ApplyCPhiReductionProcess::ExecuteFinalizeSolutionStep")
            << "Reached safety factor sofar = " << 1.0 / mReductionFactor << std::endl;
    }

    void ExecuteFinalize() override
    {
        KRATOS_INFO("ApplyCPhiReductionProcess")
            << "Final safety factor = " << 1.0 / mReductionFactor << std::endl;
    }

private:
    ModelPart& mrModelPart;
    double     mReductionFactor         = 1.0;
    double     mPreviousReductionFactor = 1.0;
    double     mReductionIncrement      = 0.1;

    double GetAndCheckPhi(const Element::PropertiesType& rProp)
    {
        // Get the initial properties from the model part. Recall that we create a separate
        // properties object with reduced c and phi for each and every element. Those reduced
        // properties objects are not linked to the original ones.
        const auto& part_properties = mrModelPart.GetProperties(rProp.Id());

        // Check for UMAT PHI Parameter
        KRATOS_ERROR_IF(!part_properties.Has(INDEX_OF_UMAT_PHI_PARAMETER) ||
                        !part_properties.Has(NUMBER_OF_UMAT_PARAMETERS) || !part_properties.Has(UMAT_PARAMETERS))
            << "Insufficient material data for C-Phi reduction process: " << std::endl;
        KRATOS_ERROR_IF(part_properties[INDEX_OF_UMAT_PHI_PARAMETER] < 1 || part_properties[INDEX_OF_UMAT_PHI_PARAMETER] > part_properties[NUMBER_OF_UMAT_PARAMETERS])
            << "undefined INDEX_OF_UMAT_PHI_PARAMETER: " << part_properties[INDEX_OF_UMAT_PHI_PARAMETER]
            << std::endl;
        const double phi =
            part_properties[UMAT_PARAMETERS][part_properties[INDEX_OF_UMAT_PHI_PARAMETER] - 1];
        KRATOS_ERROR_IF(phi < 0. || phi > 90.) << "Friction angle Phi out of range: " << phi << std::endl;
        return phi;
    }

    double ComputeReducedPhi(double Phi) const
    {
        // Phi converted to radians and then its tangent is reduced by the reduction factor
        const double tan_phi         = std::tan(MathUtils<>::DegreesToRadians(Phi));
        const double reduced_tan_phi = tan_phi * mReductionFactor;
        return std::atan(reduced_tan_phi) * 180.0 / Globals::Pi;
    }

    double GetAndCheckC(const Element::PropertiesType& rProp)
    {
        // Get the initial properties from the model part. Recall that we create a separate
        // properties object with reduced c and phi for each and every element. Those reduced
        // properties objects are not linked to the original ones.
        const auto& part_properties = mrModelPart.GetProperties(rProp.Id());
        KRATOS_ERROR_IF(!part_properties.Has(INDEX_OF_UMAT_C_PARAMETER) ||
                        !part_properties.Has(NUMBER_OF_UMAT_PARAMETERS) || !part_properties.Has(UMAT_PARAMETERS))
            << "Insufficient material data for C-phi reduction: " << std::endl;
        KRATOS_ERROR_IF(part_properties[INDEX_OF_UMAT_C_PARAMETER] < 1 || part_properties[INDEX_OF_UMAT_C_PARAMETER] > part_properties[NUMBER_OF_UMAT_PARAMETERS])
            << "undefined INDEX_OF_UMAT_C_PARAMETER: " << part_properties[INDEX_OF_UMAT_C_PARAMETER]
            << std::endl;
        const auto c = part_properties[UMAT_PARAMETERS][part_properties[INDEX_OF_UMAT_C_PARAMETER] - 1];
        KRATOS_ERROR_IF(c < 0.) << "Cohesion C out of range: " << c << std::endl;
        return c;
    }

    double GetAndCheckYoung(const Element::PropertiesType& rProp)
    {
        const auto young = mrModelPart.GetProperties(rProp.Id())[UMAT_PARAMETERS][0];
        KRATOS_ERROR_IF(young < 0.)
            << "Positive value expected for Youngs modulus UMAT_PARAMETERS(1) " << young << std::endl;
        return young;
    }

    double GetAndCheckPoisson(const Element::PropertiesType& rProp)
    {
        const auto nu = mrModelPart.GetProperties(rProp.Id())[UMAT_PARAMETERS][1];
        KRATOS_ERROR_IF(nu < -1. || nu >= 0.5)
            << "Value between -1.0 and 0.5 expected for Poissons ratio UMAT_PARAMETERS(2) " << nu
            << std::endl;
        return nu;
    }

    void SetCPhiAtElement(Element& rElement, double ReducedPhi, double ReducedC) const
    {
        // Get C/Phi material properties of this element
        const auto& r_prop = rElement.GetProperties();

        // Overwrite C and Phi in the UMAT_PARAMETERS
        auto Umat_parameters                                     = r_prop[UMAT_PARAMETERS];
        Umat_parameters[r_prop[INDEX_OF_UMAT_PHI_PARAMETER] - 1] = ReducedPhi;
        Umat_parameters[r_prop[INDEX_OF_UMAT_C_PARAMETER] - 1]   = ReducedC;

        SetValueAtElement(rElement, UMAT_PARAMETERS, Umat_parameters);
    }

    void SetValueAtElement(Element& rElement, const Variable<Vector>& rVariable, const Vector& rValue) const
    {
        // Copies properties
        Properties::Pointer p_new_prop = Kratos::make_shared<Properties>(rElement.GetProperties());

        // Adds new properties to the element
        p_new_prop->SetValue(rVariable, rValue);
        rElement.SetProperties(p_new_prop);
    }
};

} // namespace Kratos