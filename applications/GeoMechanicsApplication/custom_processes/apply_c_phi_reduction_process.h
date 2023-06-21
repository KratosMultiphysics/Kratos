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

        //--------------------------------------------------------------------------------------------------------------

        /// Constructor
        ApplyCPhiReductionProcess(ModelPart&  model_part,
                                const Parameters& ) : Process(Flags()), mrModelPart(model_part)
        {
        }

        ///------------------------------------------------------------------------------------

        /// Destructor
        ~ApplyCPhiReductionProcess() override = default;

        //--------------------------------------------------------------------------------------------------------------

        void ExecuteInitialize() override
        {
            KRATOS_TRY
            // Apply C/Phi Reduction procedure for the model part:
            block_for_each(mrModelPart.Elements(), [this](Element& rElement) {
                set_C_phi(rElement);
            });

            KRATOS_CATCH("")
        }

        //--------------------------------------------------------------------------------------------------------------
    private:
        ModelPart& mrModelPart;

        void set_C_phi(Element& rElement)
        {
            // Get C/Phi material properties of this element
            Element::PropertiesType& rProp = rElement.GetProperties();
            ConstitutiveLaw::Pointer pConstitutiveLaw = rProp.GetValue(CONSTITUTIVE_LAW);

            // Check for UMAT PHI Parameter
            double phi = GetAndCheckPhi(rProp);

            // Check for UMAT C Parameter
            double c = GetAndCheckC (rProp);

            // C/Phi reduction factor
            double reductionFactor = 0.9;

            // Phi converted to radians and then its tangent is reduced by the reduction factor
            double phi_rad = MathUtils<>::DegreesToRadians(phi);
            double tan_phi = std::tan(phi_rad);
            double reduced_tan_phi = reductionFactor * tan_phi;
            double reduced_phi_rad = std::atan(reduced_tan_phi);
            double reduced_phi = MathUtils<>::RadiansToDegrees(reduced_phi_rad); // TODO: RADIANSTODEGREES function!

            // C is reduced by the reduction factor
            double reduced_c = reductionFactor * c;

            // Does this set the reduced values as the new values?
            rProp[UMAT_PARAMETERS][rProp[INDEX_OF_UMAT_PHI_PARAMETER] - 1] = reduced_phi;
            rProp[UMAT_PARAMETERS][rProp[INDEX_OF_UMAT_C_PARAMETER] - 1] = reduced_c;

            SetValueAtElement(rElement, UMAT_PARAMETERS, rProp[UMAT_PARAMETERS]);


        }

        double GetAndCheckPhi(const Element::PropertiesType& rProp)
        {
            // Check for UMAT PHI Parameter
            double phi = 0.;
            if (rProp.Has(INDEX_OF_UMAT_PHI_PARAMETER) && rProp.Has(NUMBER_OF_UMAT_PARAMETERS) && \
                rProp.Has(UMAT_PARAMETERS)) {
                if (rProp[INDEX_OF_UMAT_PHI_PARAMETER] < 1 ||
                    rProp[INDEX_OF_UMAT_PHI_PARAMETER] > rProp[NUMBER_OF_UMAT_PARAMETERS]) {
                    KRATOS_ERROR << "undefined INDEX_OF_UMAT_PHI_PARAMETER: "
                    << rProp[INDEX_OF_UMAT_PHI_PARAMETER] << std::endl;
                }
                // needs more checking?
                phi = rProp[UMAT_PARAMETERS][rProp[INDEX_OF_UMAT_PHI_PARAMETER] - 1];
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
            double c = 0.;
            if (rProp.Has(INDEX_OF_UMAT_C_PARAMETER) && rProp.Has(NUMBER_OF_UMAT_PARAMETERS) && \
                rProp.Has(UMAT_PARAMETERS)) {
                if (rProp[INDEX_OF_UMAT_C_PARAMETER] < 1 || rProp[INDEX_OF_UMAT_C_PARAMETER] > \
                    rProp[NUMBER_OF_UMAT_PARAMETERS]) {
                    KRATOS_ERROR << "undefined INDEX_OF_UMAT_C_PARAMETER: " \
                    << rProp[INDEX_OF_UMAT_C_PARAMETER] << std::endl;
                }
                // needs more checking?
                c = rProp[UMAT_PARAMETERS][rProp[INDEX_OF_UMAT_C_PARAMETER] - 1];
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

            // Copies properties
            Properties::Pointer p_new_prop = Kratos::make_shared<Properties>(r_prop);

            // Adds new properties to the element
            p_new_prop->SetValue(rVar, Value);
            rElement.SetProperties(p_new_prop);

        }

    }; // class ApplyCPhiReductionProcess
} // namespace Kratos
