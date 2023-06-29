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

            double phi = 0.;
            double reduced_phi = 0.;
            double c = 0.;
            double reduced_c = 0.;
            unsigned int previousPropertyId = 0;
            // Apply C/Phi Reduction procedure for the model part:
            block_for_each(mrModelPart.Elements(), [this,&phi,&reduced_phi,&c,&reduced_c,&previousPropertyId](Element& rElement) {
                if (mrModelPart.GetProperties(rElement.GetProperties().Id()).Id() != previousPropertyId)
                {
                    phi         = GetAndCheckPhi(rElement.GetProperties());
                    KRATOS_INFO("ApplyCPhiReductionProcess") << "Initial Phi = " << phi << std::endl;
                    reduced_phi = ComputeReducedPhi(phi);
                    KRATOS_INFO("ApplyCPhiReductionProcess") << "Reduced Phi = " << reduced_phi << std::endl;
                    c           = GetAndCheckC(rElement.GetProperties());
                    KRATOS_INFO("ApplyCPhiReductionProcess") << "Initial C = " << c << std::endl;
                    reduced_c   = mReductionFactor * c;
                    previousPropertyId = mrModelPart.GetProperties(rElement.GetProperties().Id()).Id();
                    KRATOS_INFO("ApplyCPhiReductionProcess") << "Reduced C = " << reduced_c << std::endl;
                }
                set_C_Phi_At_Element(rElement, reduced_phi, reduced_c);
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

        double ComputeReducedPhi(const double phi)
        {
            // Phi converted to radians and then its tangent is reduced by the reduction factor
            double phi_rad = MathUtils<>::DegreesToRadians(phi);
            double tan_phi = std::tan(phi_rad);
            double reduced_tan_phi = mReductionFactor * tan_phi;
            double reduced_phi_rad = std::atan(reduced_tan_phi);
            return reduced_phi_rad * 180. / Globals::Pi; // TODO: RADIANSTODEGREES function!
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

        void set_C_Phi_At_Element(Element& rElement, const double reduced_phi, const double reduced_c)
        {
            KRATOS_INFO("ApplyCPhiReductionProcess") << "Element ID: " << rElement.Id() << std::endl;
            // Get C/Phi material properties of this element
            Element::PropertiesType& rProp = rElement.GetProperties();

            auto newParameters = rProp[UMAT_PARAMETERS];
            newParameters[rProp[INDEX_OF_UMAT_PHI_PARAMETER]-1] = reduced_phi;
            newParameters[rProp[INDEX_OF_UMAT_C_PARAMETER]-1]   = reduced_c;

            SetValueAtElement(rElement, UMAT_PARAMETERS, newParameters);
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

    };

}