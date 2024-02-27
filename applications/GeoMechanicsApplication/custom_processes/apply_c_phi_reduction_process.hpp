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

    ApplyCPhiReductionProcess(ModelPart& model_part, const Parameters&)
        : Process(Flags()), mrModelPart(model_part)
    {
    }

    ~ApplyCPhiReductionProcess() override = default;

    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY
        mReductionFactor -= mReductionIncrement;

        double    phi                = 0.;
        double    reduced_phi        = 0.;
        double    c                  = 0.;
        double    reduced_c          = 0.;
        IndexType previousPropertyId = -1;
        // Apply C/Phi Reduction procedure for the model part:
        block_for_each(mrModelPart.Elements(),
                       [this, &phi, &reduced_phi, &c, &reduced_c, &previousPropertyId](Element& rElement) {
            // Only compute new c and phi if the Id changes
            if (mrModelPart.GetProperties(rElement.GetProperties().Id()).Id() != previousPropertyId) {
                phi                = GetAndCheckPhi(rElement.GetProperties());
                reduced_phi        = ComputeReducedPhi(phi);
                c                  = GetAndCheckC(rElement.GetProperties());
                reduced_c          = mReductionFactor * c;
                previousPropertyId = mrModelPart.GetProperties(rElement.GetProperties().Id()).Id();
                KRATOS_INFO("ApplyCPhiReductionProces")
                    << "c = " << reduced_c << " (" << c << ") phi = " << reduced_phi << " (" << phi
                    << ")" << std::endl;
            }
            set_C_Phi_At_Element(rElement, reduced_phi, reduced_c);

            // Back compute strain from resident element Cauchy stress, E = UMAT_PARAMETERS[0], nu = UMAT_PARAMETERS[1]
            // ComputeElementRHSContribution(rElement);

            // Compute RHS ( the unbalance ) from resident Cauchy stress and stress based on reduced C and Phi
        });
        KRATOS_CATCH("")
    }

    void ExecuteFinalizeSolutionStep() override
    {
        bool cycle = false;
        if (cycle) {
            mReductionFactor = mPreviousReductionFactor;
            mReductionIncrement *= 0.5;
        } else {
            mPreviousReductionFactor = mReductionFactor;
        }
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
        double phi = 0.;
        if (part_properties.Has(INDEX_OF_UMAT_PHI_PARAMETER) &&
            part_properties.Has(NUMBER_OF_UMAT_PARAMETERS) && part_properties.Has(UMAT_PARAMETERS)) {
            if (part_properties[INDEX_OF_UMAT_PHI_PARAMETER] < 1 ||
                part_properties[INDEX_OF_UMAT_PHI_PARAMETER] > part_properties[NUMBER_OF_UMAT_PARAMETERS]) {
                KRATOS_ERROR << "undefined INDEX_OF_UMAT_PHI_PARAMETER: "
                             << part_properties[INDEX_OF_UMAT_PHI_PARAMETER] << std::endl;
            }
            // needs more checking?
            phi = part_properties[UMAT_PARAMETERS][part_properties[INDEX_OF_UMAT_PHI_PARAMETER] - 1];
            if (phi < 0. || phi > 90.) {
                KRATOS_ERROR << "Friction angle Phi out of range: " << phi << std::endl;
            }
        } else {
            KRATOS_ERROR << "Insufficient material data for C-Phi reduction process: " << std::endl;
        }
        return phi;
    }

    double ComputeReducedPhi(const double phi) const
    {
        // Phi converted to radians and then its tangent is reduced by the reduction factor
        double phi_rad         = MathUtils<>::DegreesToRadians(phi);
        double tan_phi         = std::tan(phi_rad);
        double reduced_tan_phi = mReductionFactor * tan_phi;
        double reduced_phi_rad = std::atan(reduced_tan_phi);
        return reduced_phi_rad * 180. / Globals::Pi; // TODO: RADIANSTODEGREES function!
    }

    double GetAndCheckC(const Element::PropertiesType& rProp)
    {
        // Get the initial properties from the model part. Recall that we create a separate
        // properties object with reduced c and phi for each and every element. Those reduced
        // properties objects are not linked to the original ones.
        const auto& part_properties = mrModelPart.GetProperties(rProp.Id());

        double c = 0.;
        if (part_properties.Has(INDEX_OF_UMAT_C_PARAMETER) &&
            part_properties.Has(NUMBER_OF_UMAT_PARAMETERS) && part_properties.Has(UMAT_PARAMETERS)) {
            if (part_properties[INDEX_OF_UMAT_C_PARAMETER] < 1 ||
                part_properties[INDEX_OF_UMAT_C_PARAMETER] > part_properties[NUMBER_OF_UMAT_PARAMETERS]) {
                KRATOS_ERROR << "undefined INDEX_OF_UMAT_C_PARAMETER: "
                             << part_properties[INDEX_OF_UMAT_C_PARAMETER] << std::endl;
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
        double      young           = part_properties[UMAT_PARAMETERS][0];
        if (young < 0.)
            KRATOS_ERROR << "Positive value expected for Youngs modulus UMAT_PARAMETERS(1) "
                         << young << std::endl;
        return young;
    }

    double GetAndCheckPoisson(const Element::PropertiesType& rProp)
    {
        const auto& part_properties = mrModelPart.GetProperties(rProp.Id());
        double      nu              = part_properties[UMAT_PARAMETERS][1];
        if (nu < -1. || nu >= 0.5)
            KRATOS_ERROR
                << "Value between -1.0 and 0.5 expected for Poissons ratio UMAT_PARAMETERS(2) "
                << nu << std::endl;
        return nu;
    }

    Matrix FormElasticConstitutiveTensor(const double young, const double poisson) const
    {
        Matrix       C  = zero_matrix(VOIGT_SIZE_2D_PLANE_STRAIN, VOIGT_SIZE_2D_PLANE_STRAIN);
        const double c0 = young / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
        const double c1 = (1.0 - poisson) * c0;
        const double c2 = c0 * poisson;
        const double c3 = (0.5 - poisson) * c0;

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

    void set_C_Phi_At_Element(Element& rElement, const double reduced_phi, const double reduced_c) const
    {
        // Get C/Phi material properties of this element
        Element::PropertiesType& rProp = rElement.GetProperties();

        // Overwrite C and Phi in the UMAT_PARAMETERS
        auto newParameters                                    = rProp[UMAT_PARAMETERS];
        newParameters[rProp[INDEX_OF_UMAT_PHI_PARAMETER] - 1] = reduced_phi;
        newParameters[rProp[INDEX_OF_UMAT_C_PARAMETER] - 1]   = reduced_c;

        // Write back to the element
        SetValueAtElement(rElement, UMAT_PARAMETERS, newParameters);
    }

    void SetValueAtElement(Element& rElement, const Variable<Vector>& rVar, const Vector& Value) const
    {
        Properties& r_prop = rElement.GetProperties();
        // Copies properties
        Properties::Pointer p_new_prop = Kratos::make_shared<Properties>(r_prop);

        // Adds new properties to the element
        p_new_prop->SetValue(rVar, Value);
        rElement.SetProperties(p_new_prop);
    }

    void ComputeElementRHSContribution(Element& rElement)
    {
        Element::PropertiesType& rProp   = rElement.GetProperties();
        const double             young   = GetAndCheckYoung(rProp);
        const double             poisson = GetAndCheckPoisson(rProp);
        Matrix                   c       = FormElasticConstitutiveTensor(young, poisson);
        Matrix                   c_inverse;
        double                   det_c = 0.;
        MathUtils<double>::InvertMatrix(c, c_inverse, det_c);

        // Get element stress vectors
        const ProcessInfo& rCurrentProcessInfo = this->mrModelPart.GetProcessInfo();
        std::vector<ConstitutiveLaw::StressVectorType> rCurrentStressVectors;
        rElement.CalculateOnIntegrationPoints(CAUCHY_STRESS_VECTOR, rCurrentStressVectors, rCurrentProcessInfo);

        // walk the integration points
        using NodeType            = Node;
        using GeometryType        = Geometry<NodeType>;
        const GeometryType& rGeom = rElement.GetGeometry();

        // access to the constitutive law:

        // Definition of variables     dit is het hele vervelende protected struct met dingen gedefinieerd voor bijv. het smallstrainUPwDiffOrderElement
        // ElementVariables Variables;
        // rElement.InitializeElementVariables(Variables, rCurrentProcessInfo);

        // Create constitutive law parameters:
        ConstitutiveLaw::Parameters ConstitutiveParameters(rGeom, rProp, rCurrentProcessInfo);
        ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
        ConstitutiveParameters.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);

        // Loop over integration points
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints =
            rGeom.IntegrationPoints(rElement.GetIntegrationMethod());
        const IndexType number_of_integration_points = IntegrationPoints.size();
        for (IndexType integration_point = 0; integration_point < number_of_integration_points; ++integration_point) {
            Vector back_strain        = prod(c_inverse, rCurrentStressVectors[integration_point]);
            auto   p_constitutive_law = rProp[CONSTITUTIVE_LAW]->Clone();
            // Voer BackStrain aan constitutive law voor een gereduceerde spanning
            //    deze constitutive heeft de gereduceerde C en Phi en geen rekhistorie
            // set gauss points variables to constitutive law parameters
            // this->SetConstitutiveParameters(Variables, ConstitutiveParameters);
            // in het small_strain_U_Pw_diff_order element zit dit op regel 1423 als
            // compute constitutive stresses
            // ConstitutiveParameters.SetStressVector(mStressVector[integration_point]);
            // mConstitutiveLawVector[GPoint]->CalculateMaterialResponseCauchy(ConstitutiveParameters);
            // B ( Sig_current - Sig_reduced ) dV, tel op bij element RHS
            //    dit zit in CalculateAndAddStiffnessForce van de Geo continuum elementen
        }
        // tel element RHS op bij system RHS via de equationId's van het element
    }
};

} // namespace Kratos