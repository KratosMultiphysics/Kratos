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
            mReductionFactor -= mReductionIncrement;

            double phi                      = 0.;
            double reduced_phi              = 0.;
            double c                        = 0.;
            double reduced_c                = 0.;
            unsigned int previousPropertyId = 0;
            // Apply C/Phi Reduction procedure for the model part:
            block_for_each(mrModelPart.Elements(), [this,&phi,&reduced_phi,&c,&reduced_c,&previousPropertyId](Element& rElement) {
                // Only compute new c and phi if the Id changes
                if (mrModelPart.GetProperties(rElement.GetProperties().Id()).Id() != previousPropertyId)
                {
                    phi         = GetAndCheckPhi(rElement.GetProperties());
                    reduced_phi = ComputeReducedPhi(phi);
                    c           = GetAndCheckC(rElement.GetProperties());
                    reduced_c   = mReductionFactor * c;
                    previousPropertyId = mrModelPart.GetProperties(rElement.GetProperties().Id()).Id();
                    KRATOS_INFO("ApplyCPhiReductionProces") << "c = " << reduced_c << " (" << c << ") phi = " << reduced_phi <<" (" << phi << ")" <<std::endl;
                }
                set_C_Phi_At_Element(rElement, reduced_phi, reduced_c);

                // Back compute strain from resident element Cauchy stress, E = UMAT_PARAMETERS[0], nu = UMAT_PARAMETERS[1]
                backComputeElementStrains(rElement);

                // Compute RHS ( the unbalance ) from resident Cauchy stress and stress based on reduced C and Phi
            });
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

        void ExecuteFinalize() override
        {
            KRATOS_INFO("ApplyCPhiReductionProcess") << "Final safety factor = " << 1.0 / mReductionFactor << std::endl;
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
                if (c < 0.) KRATOS_ERROR << "Cohesion C out of range: " << c << std::endl;
            } else {
                KRATOS_ERROR << "Insufficient material data for C-phi reduction: " << std::endl;
            }
            return c;
        }

        double GetAndCheckYoung(const Element::PropertiesType& rProp)
        {
            const auto& part_properties = mrModelPart.GetProperties(rProp.Id());
            double young = part_properties[UMAT_PARAMETERS][0];
            if (young < 0.) KRATOS_ERROR << "Positive value expected for Youngs modulus UMAT_PARAMETERS(1) " << young << std::endl;
            return young;
        }

        double GetAndCheckPoisson(const Element::PropertiesType& rProp)
        {
            const auto& part_properties = mrModelPart.GetProperties(rProp.Id());
            double nu = part_properties[UMAT_PARAMETERS][0];
            if (nu < -1. || nu >= 0.5) KRATOS_ERROR << "Value between -1.0 and 0.5 expected for Poissons ratio UMAT_PARAMETERS(2) " << nu << std::endl;
            return nu;
        }

        Matrix& FormElasticConstitutiveTensor(const double young, const double poisson)
        {
            Matrix C = zero_matrix(VOIGT_SIZE_2D_PLANE_STRAIN, VOIGT_SIZE_2D_PLANE_STRAIN);
            const double c0 = young / ((1.0 + poisson)*(1.0 - 2.0 * poisson));
            const double c1 = (1.0 - poisson)*c0;
            const double c2 = c0 * poisson;
            const double c3 = (0.5 - poisson)*c0;

            C(INDEX_2D_PLANE_STRAIN_XX, INDEX_2D_PLANE_STRAIN_XX) = c1;
            C(INDEX_2D_PLANE_STRAIN_XX, INDEX_2D_PLANE_STRAIN_YY) = c2;
            C(INDEX_2D_PLANE_STRAIN_XX, INDEX_2D_PLANE_STRAIN_ZZ) = c2;

            C(INDEX_2D_PLANE_STRAIN_YY, INDEX_2D_PLANE_STRAIN_XX) = c2;
            C(INDEX_2D_PLANE_STRAIN_YY, INDEX_2D_PLANE_STRAIN_YY) = c1;
            C(INDEX_2D_PLANE_STRAIN_YY, INDEX_2D_PLANE_STRAIN_ZZ) = c2;

            C(INDEX_2D_PLANE_STRAIN_ZZ, INDEX_2D_PLANE_STRAIN_XX) = c2;
            C(INDEX_2D_PLANE_STRAIN_ZZ, INDEX_2D_PLANE_STRAIN_YY) = c2;
            C(INDEX_2D_PLANE_STRAIN_ZZ, INDEX_2D_PLANE_STRAIN_ZZ) = c1;

            C(INDEX_2D_PLANE_STRAIN_XY, INDEX_2D_PLANE_STRAIN_XY) = c3;

            return C;
        }

        void set_C_Phi_At_Element(Element& rElement, const double reduced_phi, const double reduced_c)
        {
            // Get C/Phi material properties of this element
            Element::PropertiesType& rProp = rElement.GetProperties();

            // Overwrite C and Phi in the UMAT_PARAMETERS
            auto newParameters = rProp[UMAT_PARAMETERS];
            newParameters[rProp[INDEX_OF_UMAT_PHI_PARAMETER]-1] = reduced_phi;
            newParameters[rProp[INDEX_OF_UMAT_C_PARAMETER]-1]   = reduced_c;

            // Write back to the element
            SetValueAtElement(rElement, UMAT_PARAMETERS, newParameters);
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

        void backComputeElementStrains(Element& rElement)
        {
            Element::PropertiesType& rProp = rElement.GetProperties();
            const double young = GetAndCheckYoung(rProp);
            const double poisson = GetAndCheckPoisson(rProp);
            Matrix C = FormElasticConstitutiveTensor(young, poisson);
            Matrix CInverse = inverse(C);

            //walk the integration points
            using NodeType       = Node;
            using GeometryType   = Geometry<NodeType>;

            const GeometryType& rGeom   = rElement.GetGeometry();
            const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
                rGeom.IntegrationPoints(rElement.GetIntegrationMethod());
            const IndexType NumGPoints = IntegrationPoints.size();
            for (IndexType GPoint = 0; GPoint < NumGPoints; ++GPoint) {
                Vector currentStress = getStressVector[GPoint];
                Vector backStrain=  prod(CInverse, currentStress);
                // haal B matrix en i.p. gewicht, tel op bij element RHS
            }
            // tel element RHS op bij systeem RHS
        }
    };

}