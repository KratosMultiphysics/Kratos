// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

// System includes

// Project includes
#include "utilities/math_utils.h"
#include "constitutive_laws_application_variables.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "generic_small_strain_isotropic_viscoplasticity.h"
#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_plasticity.h"

// Yield surfaces
#include "custom_constitutive/auxiliary_files/yield_surfaces/generic_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/mohr_coulomb_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/rankine_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/auxiliary_files/plastic_potentials/generic_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/tresca_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/drucker_prager_plastic_potential.h"

namespace Kratos
{

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicViscoPlasticity<TConstLawIntegratorType>::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    // We get the constitutive tensor
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    // Integrate Stress plasticity
    const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());

    if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        BaseType::CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }
    this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        CalculateElasticMatrix(r_constitutive_matrix, rValues);

        double threshold           = GetThreshold();
        double plastic_dissipation = GetPlasticDissipation();
        Vector plastic_strain      = GetPlasticStrain();

        BoundedArrayType predictive_stress_vector, deviatoric_stress_vector;

        noalias(predictive_stress_vector) = prod(r_constitutive_matrix, r_strain_vector - plastic_strain);
        this->template AddInitialStressVectorContribution<BoundedArrayType>(predictive_stress_vector);

        double equivalent_stress;
        ConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, equivalent_stress, rValues);

        const double F = equivalent_stress - threshold;

        if (F >= 0.0) {
            BoundedArrayType deviatoric_stress_vector;
            double I1, J2;
            ConstitutiveLawUtilities<VoigtSize>::CalculateI1Invariant<BoundedArrayType>(predictive_stress_vector, I1);
            ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant<BoundedArrayType>(predictive_stress_vector, I1, deviatoric_stress_vector, J2);

            const auto& r_props = rValues.GetMaterialProperties();
            const double mu = r_props[MIU];
            const double sensitivity = r_props.Has(DP_EPSILON) ? r_props[DP_EPSILON] : 1.0; 
            // Perzyna model
            const double plastic_multiplier = (std::pow(equivalent_stress / threshold, 1.0 / sensitivity) - 1.0) / mu;

            array_1d<double, VoigtSize> g_flux;
            ConstLawIntegratorType::YieldSurfaceType::CalculatePlasticPotentialDerivative(predictive_stress_vector, deviatoric_stress_vector, J2, g_flux, rValues);

            const array_1d<double, VoigtSize> plastic_strain_increment = plastic_multiplier * g_flux;
            noalias(rValues.GetStressVector()) = predictive_stress_vector - prod(r_constitutive_matrix, plastic_strain_increment);

            if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
                CalculateAnalyticalTangentConstitutiveMatrix(rValues);
            }
        } else {
            noalias(rValues.GetStressVector()) = predictive_stress_vector;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/


template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicViscoPlasticity<TConstLawIntegratorType>::CalculateAnalyticalTangentConstitutiveMatrix(
    ConstitutiveLaw::Parameters& rValues
)
{
    const auto& r_props = rValues.GetMaterialProperties();
    const double young = r_props[YOUNG_MODULUS];
    const double nu = r_props[POISSON_RATIO];
    const double mu = r_props[MIU];
    const double sensitivity = r_props.Has(DP_EPSILON) ? r_props[DP_EPSILON] : 1.0; 
    const double threshold = GetThreshold();
    const auto& r_strain_vector = rValues.GetStrainVector();
    const auto plastic_strain = GetPlasticStrain();

    auto& r_tangent = rValues.GetConstitutiveMatrix();

    const double cr_tangent0 = nu - 1.0;
    const double cr_tangent1 = 1.0/mu;
    const double cr_tangent2 = nu + 1.0;
    const double cr_tangent3 = 1.0/(cr_tangent2*cr_tangent2);
    const double cr_tangent4 = plastic_strain[0] - r_strain_vector[0];
    const double cr_tangent5 = -cr_tangent4;
    const double cr_tangent6 = nu;
    const double cr_tangent7 = plastic_strain[1] - r_strain_vector[1];
    const double cr_tangent8 = -cr_tangent7;
    const double cr_tangent9 = plastic_strain[2] - r_strain_vector[2];
    const double cr_tangent10 = -cr_tangent9;
    const double cr_tangent11 = -cr_tangent0;
    const double cr_tangent12 = cr_tangent11*cr_tangent8;
    const double cr_tangent13 = 0.5*cr_tangent12;
    const double cr_tangent14 = cr_tangent10*cr_tangent11;
    const double cr_tangent15 = 0.5*cr_tangent14;
    const double cr_tangent16 = -0.5*cr_tangent10*nu - cr_tangent11*cr_tangent5 + cr_tangent13 + cr_tangent15 + cr_tangent5*cr_tangent6 - 0.5*cr_tangent8*nu;
    const double cr_tangent17 = nu - 0.5;
    const double cr_tangent18 = -cr_tangent17;
    const double cr_tangent19 = 1.0/(cr_tangent18*cr_tangent18);
    const double cr_tangent20 = 0.22222222222222224*cr_tangent19;
    const double cr_tangent21 = plastic_strain[3] - r_strain_vector[3];
    const double cr_tangent22 = cr_tangent21*cr_tangent21;
    const double cr_tangent23 = plastic_strain[4] - r_strain_vector[4];
    const double cr_tangent24 = cr_tangent23*cr_tangent23;
    const double cr_tangent25 = plastic_strain[5] - r_strain_vector[5];
    const double cr_tangent26 = cr_tangent25*cr_tangent25;
    const double cr_tangent27 = cr_tangent0*cr_tangent7;
    const double cr_tangent28 = 0.5*nu;
    const double cr_tangent29 = 0.5*cr_tangent0*cr_tangent4 + cr_tangent28*cr_tangent4;
    const double cr_tangent30 = cr_tangent0*cr_tangent9;
    const double cr_tangent31 = cr_tangent28*cr_tangent9 + 0.5*cr_tangent30;
    const double cr_tangent32 = -cr_tangent27 + cr_tangent29 + cr_tangent31 - cr_tangent6*cr_tangent7;
    const double cr_tangent33 = 1.0/(cr_tangent17*cr_tangent17);
    const double cr_tangent34 = 0.22222222222222224*cr_tangent33;
    const double cr_tangent35 = 0.5*cr_tangent27 + cr_tangent28*cr_tangent7;
    const double cr_tangent36 = cr_tangent29 - cr_tangent30 + cr_tangent35 - cr_tangent6*cr_tangent9;
    const double cr_tangent37 = cr_tangent22 + cr_tangent24 + cr_tangent26 + cr_tangent34*(cr_tangent32*cr_tangent32) + cr_tangent34*(cr_tangent36*cr_tangent36);
    const double cr_tangent38 = cr_tangent20*(cr_tangent16*cr_tangent16) + cr_tangent37;
    const double cr_tangent39 = young*young;
    const double cr_tangent40 = cr_tangent3*cr_tangent39;
    const double cr_tangent41 = cr_tangent38*cr_tangent40;
    const double cr_tangent42 = 0.8660254037844386*1.0/threshold;
    const double cr_tangent43 = 1.0/sensitivity;
    const double cr_tangent44 = std::pow(std::sqrt(cr_tangent41)*cr_tangent42, cr_tangent43);
    const double cr_tangent45 = cr_tangent44 - 1.0;
    const double cr_tangent46 = -cr_tangent16;
    const double cr_tangent47 = -0.5*cr_tangent11*cr_tangent5 + cr_tangent28*cr_tangent5;
    const double cr_tangent48 = cr_tangent10*cr_tangent28 + cr_tangent12 - cr_tangent15 + cr_tangent47 - cr_tangent6*cr_tangent8;
    const double cr_tangent49 = -cr_tangent10*cr_tangent6 - cr_tangent13 + cr_tangent14 + cr_tangent28*cr_tangent8 + cr_tangent47;
    const double cr_tangent50 = cr_tangent20*(cr_tangent46*cr_tangent46) + cr_tangent20*(cr_tangent48*cr_tangent48) + cr_tangent20*(cr_tangent49*cr_tangent49) + cr_tangent21*cr_tangent21 + cr_tangent23*cr_tangent23 + cr_tangent25*cr_tangent25;
    const double cr_tangent51 = pow(cr_tangent40*cr_tangent50, -1.0/2.0);
    const double cr_tangent52 = 1.0/(cr_tangent18*cr_tangent18*cr_tangent18*cr_tangent18);
    const double cr_tangent53 = 4.0*nu - 2.0;
    const double cr_tangent54 = -cr_tangent53;
    const double cr_tangent55 = 2.0*nu;
    const double cr_tangent56 = cr_tangent55 - 1.0;
    const double cr_tangent57 = cr_tangent48*cr_tangent56;
    const double cr_tangent58 = cr_tangent49*cr_tangent56;
    const double cr_tangent59 = cr_tangent46*cr_tangent54 + cr_tangent57 + cr_tangent58;
    const double cr_tangent60 = cr_tangent59*cr_tangent59;
    const double cr_tangent61 = 1.0/cr_tangent50;
    const double cr_tangent62 = cr_tangent46*cr_tangent56;
    const double cr_tangent63 = cr_tangent49*cr_tangent54 + cr_tangent57 + cr_tangent62;
    const double cr_tangent64 = 0.010691671651659738*cr_tangent52*cr_tangent61;
    const double cr_tangent65 = cr_tangent59*cr_tangent64;
    const double cr_tangent66 = cr_tangent63*cr_tangent65*nu;
    const double cr_tangent67 = cr_tangent55 - 1.0;
    const double cr_tangent68 = -1.0/cr_tangent67;
    const double cr_tangent69 = cr_tangent68*(0.88888888888888895 - 1.7777777777777779*nu);
    const double cr_tangent70 = 0.8660254037844386*cr_tangent69;
    const double cr_tangent71 = cr_tangent56*cr_tangent68;
    const double cr_tangent72 = 1.0*nu - 0.49999999999999994;
    const double cr_tangent73 = cr_tangent19*cr_tangent72;
    const double cr_tangent74 = 0.25*cr_tangent73;
    const double cr_tangent75 = cr_tangent54*cr_tangent74 + cr_tangent56*cr_tangent74 + cr_tangent71;
    const double cr_tangent76 = -0.76980035891950105*cr_tangent75*nu;
    const double cr_tangent77 = cr_tangent56*cr_tangent73;
    const double cr_tangent78 = 0.19245008972987526*cr_tangent77;
    const double cr_tangent79 = cr_tangent48*cr_tangent54 + cr_tangent58 + cr_tangent62;
    const double cr_tangent80 = cr_tangent79*nu;
    const double cr_tangent81 = cr_tangent65*cr_tangent80;
    const double cr_tangent82 = cr_tangent70 + cr_tangent76 + cr_tangent78 + cr_tangent81;
    const double cr_tangent83 = 0.010691671651659738*cr_tangent52*cr_tangent60*cr_tangent61 - cr_tangent66 - cr_tangent82;
    const double cr_tangent84 = 0.25*cr_tangent33*cr_tangent56;
    const double cr_tangent85 = cr_tangent32*cr_tangent84;
    const double cr_tangent86 = cr_tangent36*cr_tangent84;
    const double cr_tangent87 = -0.25*cr_tangent16*cr_tangent19*cr_tangent53 + cr_tangent85 + cr_tangent86;
    const double cr_tangent88 = -cr_tangent87;
    const double cr_tangent89 = 1.0/cr_tangent38;
    const double cr_tangent90 = -cr_tangent79*nu;
    const double cr_tangent91 = -cr_tangent63*nu;
    const double cr_tangent92 = -cr_tangent59 - cr_tangent90 - cr_tangent91;
    const double cr_tangent93 = -0.010691671651659738*cr_tangent52*cr_tangent59*cr_tangent61*cr_tangent63;
    const double cr_tangent94 = 0.096225044864937631*cr_tangent54*cr_tangent73 + 0.38490017945975052*cr_tangent71 - 0.38490017945975052*cr_tangent75*nu + 0.096225044864937631*cr_tangent77 - 0.8660254037844386*nu*(cr_tangent20*cr_tangent56*cr_tangent72 + cr_tangent69);
    const double cr_tangent95 = cr_tangent60*cr_tangent64*nu + cr_tangent94;
    const double cr_tangent96 = -cr_tangent81 - cr_tangent93 - cr_tangent95;
    const double cr_tangent97 = cr_tangent1*cr_tangent3*young;
    const double cr_tangent98 = cr_tangent51*cr_tangent97;
    const double cr_tangent99 = cr_tangent45*nu;
    const double cr_tangent100 = cr_tangent98*cr_tangent99;
    const double cr_tangent101 = -cr_tangent59*nu;
    const double cr_tangent102 = -cr_tangent101 - cr_tangent63 - cr_tangent90;
    const double cr_tangent103 = cr_tangent43*cr_tangent44*cr_tangent89;
    const double cr_tangent104 = cr_tangent102*cr_tangent103;
    const double cr_tangent105 = cr_tangent88*nu;
    const double cr_tangent106 = 0.042766686606638953*cr_tangent19*cr_tangent98;
    const double cr_tangent107 = cr_tangent105*cr_tangent106;
    const double cr_tangent108 = cr_tangent100*cr_tangent96 + cr_tangent104*cr_tangent107;
    const double cr_tangent109 = -0.010691671651659738*cr_tangent52*cr_tangent59*cr_tangent61*cr_tangent79;
    const double cr_tangent110 = -cr_tangent109 - cr_tangent66 - cr_tangent95;
    const double cr_tangent111 = -cr_tangent101 - cr_tangent79 - cr_tangent91;
    const double cr_tangent112 = cr_tangent103*cr_tangent111;
    const double cr_tangent113 = cr_tangent100*cr_tangent110 + cr_tangent107*cr_tangent112;
    const double cr_tangent114 = 1.0 - nu;
    const double cr_tangent115 = 1.0/cr_tangent67;
    const double cr_tangent116 = young*1.0/cr_tangent2;
    const double cr_tangent117 = cr_tangent115*cr_tangent116;
    const double cr_tangent118 = cr_tangent63*cr_tangent64*cr_tangent80;
    const double cr_tangent119 = cr_tangent79*cr_tangent79;
    const double cr_tangent120 = cr_tangent119*cr_tangent64*nu + cr_tangent94;
    const double cr_tangent121 = -cr_tangent109 - cr_tangent118 - cr_tangent120;
    const double cr_tangent122 = cr_tangent115*cr_tangent53*cr_tangent68;
    const double cr_tangent123 = -cr_tangent115*cr_tangent16*cr_tangent56*cr_tangent68;
    const double cr_tangent124 = cr_tangent122*cr_tangent32 + cr_tangent123 + cr_tangent86;
    const double cr_tangent125 = -cr_tangent124;
    const double cr_tangent126 = -0.010691671651659738*cr_tangent52*cr_tangent61*cr_tangent63*cr_tangent79;
    const double cr_tangent127 = -cr_tangent120 - cr_tangent126 - cr_tangent81;
    const double cr_tangent128 = cr_tangent125*nu;
    const double cr_tangent129 = cr_tangent106*cr_tangent128;
    const double cr_tangent130 = cr_tangent100*cr_tangent127 + cr_tangent104*cr_tangent129;
    const double cr_tangent131 = -cr_tangent118 + 0.010691671651659738*cr_tangent119*cr_tangent52*cr_tangent61 - cr_tangent82;
    const double cr_tangent132 = cr_tangent100*cr_tangent131 + cr_tangent112*cr_tangent129 + nu;
    const double cr_tangent133 = cr_tangent63*cr_tangent63;
    const double cr_tangent134 = cr_tangent133*cr_tangent64*nu + cr_tangent94;
    const double cr_tangent135 = -cr_tangent118 - cr_tangent134 - cr_tangent93;
    const double cr_tangent136 = cr_tangent122*cr_tangent36 + cr_tangent123 + cr_tangent85;
    const double cr_tangent137 = -cr_tangent136;
    const double cr_tangent138 = -cr_tangent126 - cr_tangent134 - cr_tangent66;
    const double cr_tangent139 = cr_tangent137*nu;
    const double cr_tangent140 = cr_tangent106*cr_tangent139;
    const double cr_tangent141 = cr_tangent100*cr_tangent138 + cr_tangent112*cr_tangent140;
    const double cr_tangent142 = -cr_tangent118 + 0.010691671651659738*cr_tangent133*cr_tangent52*cr_tangent61 - cr_tangent66 - cr_tangent70 - cr_tangent76 - cr_tangent78;
    const double cr_tangent143 = cr_tangent100*cr_tangent142 + cr_tangent104*cr_tangent140 + nu;
    const double cr_tangent144 = cr_tangent1*cr_tangent21;
    const double cr_tangent145 = -cr_tangent128 - cr_tangent139 - cr_tangent87;
    const double cr_tangent146 = std::pow(cr_tangent41, -3.0/2.0);
    const double cr_tangent147 = 0.38490017945975052*cr_tangent146*cr_tangent40;
    const double cr_tangent148 = cr_tangent0*cr_tangent147*cr_tangent45;
    const double cr_tangent149 = 0.096225044864937631*cr_tangent19;
    const double cr_tangent150 = cr_tangent149*cr_tangent51;
    const double cr_tangent151 = cr_tangent0*cr_tangent150;
    const double cr_tangent152 = cr_tangent103*cr_tangent92;
    const double cr_tangent153 = -cr_tangent105 - cr_tangent128 - cr_tangent136;
    const double cr_tangent154 = cr_tangent147*cr_tangent99;
    const double cr_tangent155 = cr_tangent150*nu;
    const double cr_tangent156 = cr_tangent104*cr_tangent155 + cr_tangent153*cr_tangent154;
    const double cr_tangent157 = -cr_tangent105 - cr_tangent124 - cr_tangent139;
    const double cr_tangent158 = cr_tangent103*cr_tangent155;
    const double cr_tangent159 = cr_tangent111*cr_tangent158 + cr_tangent154*cr_tangent157;
    const double cr_tangent160 = cr_tangent115*cr_tangent39*1.0/(cr_tangent2*cr_tangent2*cr_tangent2);
    const double cr_tangent161 = cr_tangent160*(-cr_tangent145*cr_tangent148 - cr_tangent151*cr_tangent152 + cr_tangent156 + cr_tangent159);
    const double cr_tangent162 = cr_tangent1*cr_tangent161;
    const double cr_tangent163 = cr_tangent106*cr_tangent152;
    const double cr_tangent164 = cr_tangent100*cr_tangent83 + cr_tangent105*cr_tangent163 + nu;
    const double cr_tangent165 = cr_tangent100*cr_tangent121 + cr_tangent128*cr_tangent163;
    const double cr_tangent166 = cr_tangent100*cr_tangent135 + cr_tangent139*cr_tangent163;
    const double cr_tangent167 = cr_tangent145*cr_tangent154 + cr_tangent158*cr_tangent92;
    const double cr_tangent168 = cr_tangent160*(-cr_tangent112*cr_tangent151 - cr_tangent148*cr_tangent157 + cr_tangent156 + cr_tangent167);
    const double cr_tangent169 = cr_tangent1*cr_tangent168;
    const double cr_tangent170 = cr_tangent160*(-cr_tangent104*cr_tangent151 - cr_tangent148*cr_tangent153 + cr_tangent159 + cr_tangent167);
    const double cr_tangent171 = cr_tangent1*cr_tangent170;
    const double cr_tangent172 = cr_tangent149*cr_tangent45;
    const double cr_tangent173 = -cr_tangent172*cr_tangent59 + 0.38490017945975052*cr_tangent43*cr_tangent44*cr_tangent88;
    const double cr_tangent174 = nu + 1.0;
    const double cr_tangent175 = cr_tangent174*1.0/(cr_tangent2*cr_tangent2*cr_tangent2*cr_tangent2*cr_tangent2)*(young*young*young*young);
    const double cr_tangent176 = cr_tangent146*cr_tangent175;
    const double cr_tangent177 = cr_tangent144*cr_tangent176;
    const double cr_tangent178 = 0.38490017945975052*cr_tangent125*cr_tangent43*cr_tangent44 - cr_tangent172*cr_tangent79;
    const double cr_tangent179 = 0.38490017945975052*cr_tangent137*cr_tangent43*cr_tangent44 - cr_tangent172*cr_tangent63;
    const double cr_tangent180 = cr_tangent40*(cr_tangent34*std::pow(cr_tangent0*cr_tangent4 - cr_tangent31 - cr_tangent35 + 0.99999999999999989*cr_tangent4*nu, 2) + cr_tangent37);
    const double cr_tangent181 = std::sqrt(cr_tangent180);
    const double cr_tangent182 = std::pow(cr_tangent181*cr_tangent42, cr_tangent43);
    const double cr_tangent183 = 0.8660254037844386*cr_tangent174;
    const double cr_tangent184 = cr_tangent183*(cr_tangent182 - 1.0);
    const double cr_tangent185 = std::pow(cr_tangent180, -3.0/2.0);
    const double cr_tangent186 = cr_tangent1*cr_tangent185;
    const double cr_tangent187 = cr_tangent186*1.0/(cr_tangent2*cr_tangent2*cr_tangent2*cr_tangent2)*(young*young*young);
    const double cr_tangent188 = cr_tangent184*cr_tangent187;
    const double cr_tangent189 = cr_tangent182*cr_tangent43;
    const double cr_tangent190 = cr_tangent183*cr_tangent187*cr_tangent189;
    const double cr_tangent191 = -cr_tangent184*cr_tangent97*1.0/cr_tangent181 + 0.5;
    const double cr_tangent192 = cr_tangent182 - cr_tangent189 - 1.0;
    const double cr_tangent193 = 0.8660254037844386*cr_tangent175*cr_tangent192*cr_tangent23;
    const double cr_tangent194 = cr_tangent144*cr_tangent185;
    const double cr_tangent195 = cr_tangent193*cr_tangent194;
    const double cr_tangent196 = 0.8660254037844386*cr_tangent175*cr_tangent192*cr_tangent194*cr_tangent25;
    const double cr_tangent197 = cr_tangent1*cr_tangent176;
    const double cr_tangent198 = cr_tangent197*cr_tangent23;
    const double cr_tangent199 = cr_tangent186*cr_tangent193*cr_tangent25;
    const double cr_tangent200 = cr_tangent197*cr_tangent25;
    r_tangent(0,0)=cr_tangent117*(0.042766686606638953*cr_tangent0*cr_tangent1*cr_tangent19*cr_tangent3*cr_tangent43*cr_tangent44*cr_tangent51*cr_tangent88*cr_tangent89*cr_tangent92*young + cr_tangent0*cr_tangent1*cr_tangent3*cr_tangent45*cr_tangent51*cr_tangent83*young - cr_tangent108 - cr_tangent113 - cr_tangent114);
    r_tangent(0,1)=cr_tangent117*(cr_tangent0*cr_tangent1*cr_tangent121*cr_tangent3*cr_tangent45*cr_tangent51*young + 0.042766686606638953*cr_tangent0*cr_tangent1*cr_tangent125*cr_tangent19*cr_tangent3*cr_tangent43*cr_tangent44*cr_tangent51*cr_tangent89*cr_tangent92*young - cr_tangent130 - cr_tangent132);
    r_tangent(0,2)=cr_tangent117*(cr_tangent0*cr_tangent1*cr_tangent135*cr_tangent3*cr_tangent45*cr_tangent51*young + 0.042766686606638953*cr_tangent0*cr_tangent1*cr_tangent137*cr_tangent19*cr_tangent3*cr_tangent43*cr_tangent44*cr_tangent51*cr_tangent89*cr_tangent92*young - cr_tangent141 - cr_tangent143);
    r_tangent(0,3)=cr_tangent144*cr_tangent161;
    r_tangent(0,4)=cr_tangent162*cr_tangent23;
    r_tangent(0,5)=cr_tangent162*cr_tangent25;
    r_tangent(1,0)=cr_tangent117*(cr_tangent0*cr_tangent1*cr_tangent110*cr_tangent3*cr_tangent45*cr_tangent51*young + 0.042766686606638953*cr_tangent0*cr_tangent1*cr_tangent111*cr_tangent19*cr_tangent3*cr_tangent43*cr_tangent44*cr_tangent51*cr_tangent88*cr_tangent89*young - cr_tangent108 - cr_tangent164);
    r_tangent(1,1)=cr_tangent117*(0.042766686606638953*cr_tangent0*cr_tangent1*cr_tangent111*cr_tangent125*cr_tangent19*cr_tangent3*cr_tangent43*cr_tangent44*cr_tangent51*cr_tangent89*young + cr_tangent0*cr_tangent1*cr_tangent131*cr_tangent3*cr_tangent45*cr_tangent51*young - cr_tangent114 - cr_tangent130 - cr_tangent165);
    r_tangent(1,2)=cr_tangent117*(0.042766686606638953*cr_tangent0*cr_tangent1*cr_tangent111*cr_tangent137*cr_tangent19*cr_tangent3*cr_tangent43*cr_tangent44*cr_tangent51*cr_tangent89*young + cr_tangent0*cr_tangent1*cr_tangent138*cr_tangent3*cr_tangent45*cr_tangent51*young - cr_tangent143 - cr_tangent166);
    r_tangent(1,3)=cr_tangent144*cr_tangent168;
    r_tangent(1,4)=cr_tangent169*cr_tangent23;
    r_tangent(1,5)=cr_tangent169*cr_tangent25;
    r_tangent(2,0)=cr_tangent117*(0.042766686606638953*cr_tangent0*cr_tangent1*cr_tangent102*cr_tangent19*cr_tangent3*cr_tangent43*cr_tangent44*cr_tangent51*cr_tangent88*cr_tangent89*young + cr_tangent0*cr_tangent1*cr_tangent3*cr_tangent45*cr_tangent51*cr_tangent96*young - cr_tangent113 - cr_tangent164);
    r_tangent(2,1)=cr_tangent117*(0.042766686606638953*cr_tangent0*cr_tangent1*cr_tangent102*cr_tangent125*cr_tangent19*cr_tangent3*cr_tangent43*cr_tangent44*cr_tangent51*cr_tangent89*young + cr_tangent0*cr_tangent1*cr_tangent127*cr_tangent3*cr_tangent45*cr_tangent51*young - cr_tangent132 - cr_tangent165);
    r_tangent(2,2)=cr_tangent117*(0.042766686606638953*cr_tangent0*cr_tangent1*cr_tangent102*cr_tangent137*cr_tangent19*cr_tangent3*cr_tangent43*cr_tangent44*cr_tangent51*cr_tangent89*young + cr_tangent0*cr_tangent1*cr_tangent142*cr_tangent3*cr_tangent45*cr_tangent51*young - cr_tangent114 - cr_tangent141 - cr_tangent166);
    r_tangent(2,3)=cr_tangent144*cr_tangent170;
    r_tangent(2,4)=cr_tangent171*cr_tangent23;
    r_tangent(2,5)=cr_tangent171*cr_tangent25;
    r_tangent(3,0)=cr_tangent173*cr_tangent177;
    r_tangent(3,1)=cr_tangent177*cr_tangent178;
    r_tangent(3,2)=cr_tangent177*cr_tangent179;
    r_tangent(3,3)=cr_tangent116*(cr_tangent188*cr_tangent22 - cr_tangent190*cr_tangent22 + cr_tangent191);
    r_tangent(3,4)=cr_tangent195;
    r_tangent(3,5)=cr_tangent196;
    r_tangent(4,0)=cr_tangent173*cr_tangent198;
    r_tangent(4,1)=cr_tangent178*cr_tangent198;
    r_tangent(4,2)=cr_tangent179*cr_tangent198;
    r_tangent(4,3)=cr_tangent195;
    r_tangent(4,4)=cr_tangent116*(cr_tangent188*cr_tangent24 - cr_tangent190*cr_tangent24 + cr_tangent191);
    r_tangent(4,5)=cr_tangent199;
    r_tangent(5,0)=cr_tangent173*cr_tangent200;
    r_tangent(5,1)=cr_tangent178*cr_tangent200;
    r_tangent(5,2)=cr_tangent179*cr_tangent200;
    r_tangent(5,3)=cr_tangent196;
    r_tangent(5,4)=cr_tangent199;
    r_tangent(5,5)=cr_tangent116*(cr_tangent188*cr_tangent26 - cr_tangent190*cr_tangent26 + cr_tangent191);
}


/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
void GenericSmallStrainIsotropicViscoPlasticity<TConstLawIntegratorType>::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    // We get the constitutive tensor
    Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    // Integrate Stress plasticity
    const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());

    if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        BaseType::CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }
    this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

    CalculateElasticMatrix(r_constitutive_matrix, rValues);

    double& r_threshold           = GetThreshold();
    double& r_plastic_dissipation = GetPlasticDissipation();
    Vector& r_plastic_strain      = GetPlasticStrain();

    BoundedArrayType predictive_stress_vector, deviatoric_stress_vector;

    noalias(predictive_stress_vector) = prod(r_constitutive_matrix, r_strain_vector - r_plastic_strain);
    this->template AddInitialStressVectorContribution<BoundedArrayType>(predictive_stress_vector);

    double equivalent_stress;
    ConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, equivalent_stress, rValues);

    const double F = equivalent_stress - r_threshold;

    if (F >= 0.0) {
        BoundedArrayType deviatoric_stress_vector;
        double I1, J2;
        ConstitutiveLawUtilities<VoigtSize>::CalculateI1Invariant<BoundedArrayType>(predictive_stress_vector, I1);
        ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant<BoundedArrayType>(predictive_stress_vector, I1, deviatoric_stress_vector, J2);

        const auto& r_props = rValues.GetMaterialProperties();
        const double mu = r_props[MIU];
        const double sensitivity = r_props[DP_EPSILON];
        const double plastic_multiplier = (std::pow(equivalent_stress / r_threshold, 1.0 / sensitivity) - 1.0) / mu;

        array_1d<double, VoigtSize> g_flux;
        ConstLawIntegratorType::YieldSurfaceType::CalculatePlasticPotentialDerivative(predictive_stress_vector, deviatoric_stress_vector, J2, g_flux, rValues);

        const array_1d<double, VoigtSize> plastic_strain_increment = plastic_multiplier * g_flux;
        noalias(rValues.GetStressVector()) = predictive_stress_vector - prod(r_constitutive_matrix, plastic_strain_increment);

        const double g = r_props[FRACTURE_ENERGY] / characteristic_length;

        r_plastic_dissipation += inner_prod(rValues.GetStressVector(), plastic_strain_increment) / g;
        noalias(r_plastic_strain) += plastic_strain_increment;

        double dummy = 0.0;
        ConstLawIntegratorType::CalculateEquivalentStressThreshold(r_plastic_dissipation, 1, 0, r_threshold, dummy, rValues, dummy, characteristic_length);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template <class TConstLawIntegratorType>
int GenericSmallStrainIsotropicViscoPlasticity<TConstLawIntegratorType>::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    const int check_base = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
    int check_visco = 0;
    if (!rMaterialProperties.Has(MIU)) {
        check_visco++;
        KRATOS_ERROR << "The MIU variable is not defined in the properties..." << std::endl;
    }
    return check_visco;
}

/***********************************************************************************/
/***********************************************************************************/

template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<6>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<6>>>>;

template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<DruckerPragerPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<TrescaPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<ModifiedMohrCoulombYieldSurface<TrescaPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<DruckerPragerPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<TrescaPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<TrescaPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<VonMisesYieldSurface<MohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<VonMisesPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<DruckerPragerPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<MohrCoulombYieldSurface<TrescaPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<TrescaYieldSurface<MohrCoulombPlasticPotential<3>>>>;
template class GenericSmallStrainIsotropicViscoPlasticity<GenericConstitutiveLawIntegratorPlasticity<DruckerPragerYieldSurface<MohrCoulombPlasticPotential<3>>>>;

} // namespace Kratos
