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
#include "includes/model_part.h"
#include "utilities/math_utils.h"
#include "includes/element.h"

// Application includes
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

    class ApplyCPhiReductionProcess : public Process
    {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(ApplyCPhiReductionProcess);

        /// Constructor
        ApplyCPhiReductionProcess(ModelPart&  model_part,
                                  const Parameters& ) : Process(Flags()), mrModelPart(model_part)
        {
            mReductionFactor         = 1.;
            mPreviousReductionFactor = 1.;
            mReductionIncrement      = 0.1;
            KRATOS_INFO("ApplyCPhiReductionProcess") << "Constructor" << std::endl;
        }

        /// Destructor
        ~ApplyCPhiReductionProcess() override = default;

        void ExecuteInitializeSolutionStep() override
        {
            KRATOS_TRY
            KRATOS_INFO("ApplyCPhiReductionProcess") << "Start of Execute Initialize Solution Step" << std::endl;

            mReductionFactor -= mReductionIncrement;
            KRATOS_INFO("ApplyCPhiReductionProcess") << "Reduction factor: " << mReductionFactor << std::endl;

            // Apply C/Phi Reduction procedure for the model part:
            block_for_each(mrModelPart.Elements(), [this](Element& rElement) {
                set_C_Phi(rElement);
            });
            KRATOS_INFO("ApplyCPhiReductionProcess") << "End of Execute Initialize Solution Step" << std::endl;
            KRATOS_CATCH("")
        }

        void ExecuteFinalizeSolutionStep() override
        {
            bool cycle = false;
            if (cycle)
            {
                mReductionFactor = mPreviousReductionFactor;
                mReductionIncrement *= 0.5;
            } else {
                mPreviousReductionFactor = mReductionFactor;
            }
        }

    private:
        ModelPart& mrModelPart;
        double mReductionFactor;
        double mPreviousReductionFactor;
        double mReductionIncrement;

        void set_C_Phi(Element& rElement)
        {
            KRATOS_INFO("ApplyCPhiReductionProcess") << "Element ID: " << rElement.Id() << std::endl;
            // Get C/Phi material properties of this element
            Element::PropertiesType& rProp = rElement.GetProperties();
            ConstitutiveLaw::Pointer pConstitutiveLaw = rProp.GetValue(CONSTITUTIVE_LAW);
//            const Properties& r_material_properties = rValues.GetMaterialProperties();
// we hbben een constitutive law pointer, wat ik moet proberen in ConstitutiveLaw::Parameters zie linear_elastic_plane_strain_K0_law.cpp

            // Check for UMAT PHI Parameter
            double phi = GetAndCheckPhi(rProp);
            KRATOS_INFO("ApplyCPhiReductionProcess") << "Initial Phi = " << phi << std::endl;

            // Check for UMAT C Parameter
            double c = GetAndCheckC(rProp);
            KRATOS_INFO("ApplyCPhiReductionProcess") << "Initial C = " << c << std::endl;

            // Phi converted to radians and then its tangent is reduced by the reduction factor
            double phi_rad = MathUtils<>::DegreesToRadians(phi);
            double tan_phi = std::tan(phi_rad);
            double reduced_tan_phi = mReductionFactor * tan_phi;
            double reduced_phi_rad = std::atan(reduced_tan_phi);
            double reduced_phi = reduced_phi_rad * 180 / Globals::Pi; // TODO: RADIANSTODEGREES function!
            KRATOS_INFO("ApplyCPhiReductionProcess") << "Reduced Phi = " << reduced_phi << std::endl;

            // C is reduced by the reduction factor
            double reduced_c = mReductionFactor * c;
            KRATOS_INFO("ApplyCPhiReductionProcess") << "Reduced C = " << reduced_c << std::endl;

            auto newParameters = rProp[UMAT_PARAMETERS];
            newParameters[rProp[INDEX_OF_UMAT_PHI_PARAMETER]-1] = reduced_phi;
            newParameters[rProp[INDEX_OF_UMAT_C_PARAMETER]-1]   = reduced_c;

            SetValueAtElement(rElement, UMAT_PARAMETERS, newParameters);

        }

        double GetAndCheckPhi(const Element::PropertiesType& rProp)
        {
            // Get the initial properties from the model part. Recall that we create a separate properties
            // object with reduced c and phi for each and every element. Those reduced properties objects
            // are not linked to the original ones.
            const auto& part_properties = mrModelPart.GetProperties(rProp.Id());

            // Check for UMAT PHI Parameter
            double phi = 0.;
            if (part_properties.Has(INDEX_OF_UMAT_PHI_PARAMETER) &&
                part_properties.Has(NUMBER_OF_UMAT_PARAMETERS) &&
                part_properties.Has(UMAT_PARAMETERS)) {
                if (part_properties[INDEX_OF_UMAT_PHI_PARAMETER] < 1 ||
                    part_properties[INDEX_OF_UMAT_PHI_PARAMETER] > part_properties[NUMBER_OF_UMAT_PARAMETERS]) {
                    KRATOS_ERROR << "undefined INDEX_OF_UMAT_PHI_PARAMETER: " << part_properties[INDEX_OF_UMAT_PHI_PARAMETER] << std::endl;
                }
                // needs more checking?
                phi = part_properties[UMAT_PARAMETERS][part_properties[INDEX_OF_UMAT_PHI_PARAMETER] - 1];
                if (phi < 0. || phi > 90.) {
                    KRATOS_ERROR << "Friction angle Phi out of range: " << phi << std::endl;
                }
            }
            else {
                KRATOS_ERROR << "Insufficient material data for C-Phi reduction process: " << std::endl;
            }
            return phi;
        }

        double GetAndCheckC(const Element::PropertiesType& rProp)
        {
            // Get the initial properties from the model part. Recall that we create a separate properties
            // object with reduced c and phi for each and every element. Those reduced properties objects
            // are not linked to the original ones.
            const auto& part_properties = mrModelPart.GetProperties(rProp.Id());

            double c = 0.;
            if (part_properties.Has(INDEX_OF_UMAT_C_PARAMETER) &&
                part_properties.Has(NUMBER_OF_UMAT_PARAMETERS) &&
                part_properties.Has(UMAT_PARAMETERS)) {
                if (part_properties[INDEX_OF_UMAT_C_PARAMETER] < 1 ||
                    part_properties[INDEX_OF_UMAT_C_PARAMETER] > part_properties[NUMBER_OF_UMAT_PARAMETERS]) {
                    KRATOS_ERROR << "undefined INDEX_OF_UMAT_C_PARAMETER: " << part_properties[INDEX_OF_UMAT_C_PARAMETER] << std::endl;
                }
                // needs more checking?
                c = part_properties[UMAT_PARAMETERS][part_properties[INDEX_OF_UMAT_C_PARAMETER] - 1];
                if (c < 0.) {
                    KRATOS_ERROR << "Cohesion C out of range: " << c << std::endl;
                }
            }
            else {
                KRATOS_ERROR << "Insufficient material data for C-phi reduction: " << std::endl;
            }
            return c;
        }

        void SetValueAtElement(Element& rElement, const Variable<Vector>& rVar, const Vector& Value)
        {

            Properties& r_prop = rElement.GetProperties();
            KRATOS_INFO("SetValueAtElement") << "Properties ID of initial instance: " << r_prop.Id() << std::endl;

            // Copies properties
            Properties::Pointer p_new_prop = Kratos::make_shared<Properties>(r_prop);
            KRATOS_INFO("SetValueAtElement") << "Properties ID of new instance: " << p_new_prop->Id() << std::endl;

            // Adds new properties to the element
            p_new_prop->SetValue(rVar, Value);
            rElement.SetProperties(p_new_prop);

            KRATOS_INFO("Retrieval of written thing ") << std::endl;
            rElement.GetProperties().PrintData(std::cout);
        }

    }; // class ApplyCPhiReductionProcess
} // namespace Kratos
