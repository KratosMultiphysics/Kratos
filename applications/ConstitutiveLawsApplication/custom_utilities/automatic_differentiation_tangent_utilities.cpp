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
//


// System includes

// External includes

// Project includes
#include "custom_utilities/automatic_differentiation_tangent_utilities.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"

// Yield surfaces
#include "custom_constitutive/auxiliary_files/yield_surfaces/von_mises_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/mohr_coulomb_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/rankine_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/simo_ju_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/drucker_prager_yield_surface.h"
#include "custom_constitutive/auxiliary_files/yield_surfaces/tresca_yield_surface.h"

// Plastic potentials
#include "custom_constitutive/auxiliary_files/plastic_potentials/rankine_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/von_mises_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/tresca_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/modified_mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/mohr_coulomb_plastic_potential.h"
#include "custom_constitutive/auxiliary_files/plastic_potentials/drucker_prager_plastic_potential.h"

#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_von_mises_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_tresca_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_rankine_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_simo_ju_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_modified_mohr_coulomb_yield_surface.h"
#include "custom_constitutive/thermal/auxiliary_files/thermal_yield_surfaces/thermal_drucker_prager_yield_surface.h"

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<VonMisesYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  const auto &r_props = rValues.GetMaterialProperties();
  const double Young = r_props[YOUNG_MODULUS];
  const double nu = r_props[POISSON_RATIO];
  const double Gf = r_props[FRACTURE_ENERGY];
  const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
  double threshold;
  YieldSurfaceType::GetInitialUniaxialThreshold(rValues, threshold);
  auto &r_Ct = rValues.GetConstitutiveMatrix();
  const auto &r_strain = rValues.GetStrainVector();

  const double cr_Ct0 = nu - 1.0;
  const double cr_Ct1 = -nu;
  const double cr_Ct2 = std::pow(cr_Ct1 + 0.5, -2);
  const double cr_Ct3 = nu * r_strain[0];
  const double cr_Ct4 = cr_Ct3;
  const double cr_Ct5 = nu * r_strain[1];
  const double cr_Ct6 = 0.5 * cr_Ct5;
  const double cr_Ct7 = -cr_Ct6;
  const double cr_Ct8 = nu * r_strain[2];
  const double cr_Ct9 = 0.5 * cr_Ct8;
  const double cr_Ct10 = -cr_Ct9;
  const double cr_Ct11 = cr_Ct1 + 1.0;
  const double cr_Ct12 = cr_Ct11 * r_strain[0];
  const double cr_Ct13 = cr_Ct11 * r_strain[1];
  const double cr_Ct14 = 0.5 * cr_Ct13;
  const double cr_Ct15 = cr_Ct11 * r_strain[2];
  const double cr_Ct16 = 0.5 * cr_Ct15;
  const double cr_Ct17 = cr_Ct10 - cr_Ct12 + cr_Ct14 + cr_Ct16 + cr_Ct4 + cr_Ct7;
  const double cr_Ct18 = cr_Ct0 * r_strain[1];
  const double cr_Ct19 = cr_Ct5;
  const double cr_Ct20 = cr_Ct0 * r_strain[2];
  const double cr_Ct21 = 0.5 * cr_Ct20;
  const double cr_Ct22 = 0.5 * cr_Ct3;
  const double cr_Ct23 = cr_Ct0 * r_strain[0];
  const double cr_Ct24 = 0.5 * cr_Ct23;
  const double cr_Ct25 = -cr_Ct22 - cr_Ct24;
  const double cr_Ct26 = cr_Ct10 + cr_Ct18 + cr_Ct19 - cr_Ct21 + cr_Ct25;
  const double cr_Ct27 = std::pow(nu - 0.5, -2);
  const double cr_Ct28 = 0.22222222222222224 * cr_Ct27;
  const double cr_Ct29 = cr_Ct8;
  const double cr_Ct30 = 0.5 * cr_Ct18;
  const double cr_Ct31 = cr_Ct20 + cr_Ct25 + cr_Ct29 - cr_Ct30 + cr_Ct7;
  const double cr_Ct32 = std::pow(r_strain[3], 2);
  const double cr_Ct33 = std::pow(r_strain[4], 2);
  const double cr_Ct34 = std::pow(r_strain[5], 2);
  const double cr_Ct35 = cr_Ct32 + cr_Ct33 + cr_Ct34;
  const double cr_Ct36 = 0.22222222222222224 * std::pow(cr_Ct17, 2) * cr_Ct2 + std::pow(cr_Ct26, 2) * cr_Ct28 + cr_Ct28 * std::pow(cr_Ct31, 2) + cr_Ct35;
  const double cr_Ct37 = nu + 1.0;
  const double cr_Ct38 = std::pow(Young, 2) / std::pow(cr_Ct37, 2);
  const double cr_Ct39 = cr_Ct36 * cr_Ct38;
  const double cr_Ct40 = std::sqrt(cr_Ct39);
  const double cr_Ct41 = threshold / cr_Ct40;
  const double cr_Ct42 = 1.1547005383792517 * cr_Ct41;
  const double cr_Ct43 = 2.0 * nu;
  const double cr_Ct44 = cr_Ct43 - 1.0;
  const double cr_Ct45 = cr_Ct27 * cr_Ct44;
  const double cr_Ct46 = 0.25 * cr_Ct45;
  const double cr_Ct47 = cr_Ct31 * cr_Ct46;
  const double cr_Ct48 = -cr_Ct47;
  const double cr_Ct49 = cr_Ct26 * cr_Ct46;
  const double cr_Ct50 = -cr_Ct49;
  const double cr_Ct51 = 4.0 * nu;
  const double cr_Ct52 = cr_Ct51 - 2.0;
  const double cr_Ct53 = cr_Ct17 * cr_Ct2 * cr_Ct52;
  const double cr_Ct54 = cr_Ct48 + cr_Ct50 + 0.25 * cr_Ct53;
  const double cr_Ct55 = -cr_Ct5;
  const double cr_Ct56 = -cr_Ct8;
  const double cr_Ct57 = cr_Ct23 + cr_Ct55 + cr_Ct56;
  const double cr_Ct58 = 1.0 / (Gf * Young / (characteristic_length * std::pow(threshold, 2)) - 0.5);
  const double cr_Ct59 = 1.0 / cr_Ct36;
  const double cr_Ct60 = cr_Ct58 * cr_Ct59;
  const double cr_Ct61 = 0.44444444444444448 * cr_Ct60;
  const double cr_Ct62 = cr_Ct57 * cr_Ct61;
  const double cr_Ct63 = 2.0 - cr_Ct51;
  const double cr_Ct64 = -cr_Ct4;
  const double cr_Ct65 = -cr_Ct16 + cr_Ct9;
  const double cr_Ct66 = -cr_Ct14 + cr_Ct6;
  const double cr_Ct67 = cr_Ct12 + cr_Ct64 + cr_Ct65 + cr_Ct66;
  const double cr_Ct68 = -cr_Ct19;
  const double cr_Ct69 = -0.5 * cr_Ct12 + cr_Ct22;
  const double cr_Ct70 = cr_Ct13 + cr_Ct65 + cr_Ct68 + cr_Ct69;
  const double cr_Ct71 = cr_Ct44 * cr_Ct70;
  const double cr_Ct72 = -cr_Ct29;
  const double cr_Ct73 = cr_Ct15 + cr_Ct66 + cr_Ct69 + cr_Ct72;
  const double cr_Ct74 = cr_Ct44 * cr_Ct73;
  const double cr_Ct75 = cr_Ct63 * cr_Ct67 + cr_Ct71 + cr_Ct74;
  const double cr_Ct76 = cr_Ct38 / std::pow(cr_Ct39, 3.0 / 2.0);
  const double cr_Ct77 = cr_Ct76 * threshold;
  const double cr_Ct78 = cr_Ct2 * cr_Ct77;
  const double cr_Ct79 = 0.12830005981991685 * cr_Ct78;
  const double cr_Ct80 = cr_Ct57 * cr_Ct79;
  const double cr_Ct81 = 1.0 / (cr_Ct43 - 1.0);
  const double cr_Ct82 = 0.8660254037844386 / threshold;
  const double cr_Ct83 = cr_Ct58;
  const double cr_Ct84 = Young / cr_Ct37;
  const double cr_Ct85 = cr_Ct84 * std::exp(cr_Ct83 * (-cr_Ct40 * cr_Ct82 + 1.0));
  const double cr_Ct86 = cr_Ct81 * cr_Ct85;
  const double cr_Ct87 = cr_Ct42 * nu;
  const double cr_Ct88 = cr_Ct44 * cr_Ct67;
  const double cr_Ct89 = cr_Ct63 * cr_Ct70 + cr_Ct74 + cr_Ct88;
  const double cr_Ct90 = cr_Ct81 / (1.0 - cr_Ct43);
  const double cr_Ct91 = cr_Ct52 * cr_Ct90;
  const double cr_Ct92 = cr_Ct26 * cr_Ct91;
  const double cr_Ct93 = cr_Ct17 * cr_Ct44 * cr_Ct90;
  const double cr_Ct94 = cr_Ct48 - cr_Ct92 + cr_Ct93;
  const double cr_Ct95 = cr_Ct63 * cr_Ct73 + cr_Ct71 + cr_Ct88;
  const double cr_Ct96 = cr_Ct31 * cr_Ct91;
  const double cr_Ct97 = cr_Ct50 + cr_Ct93 - cr_Ct96;
  const double cr_Ct98 = 1.1547005383792517 * threshold;
  const double cr_Ct99 = cr_Ct86 * (cr_Ct59 * cr_Ct83 + cr_Ct76 * cr_Ct98);
  const double cr_Ct100 = cr_Ct57 * cr_Ct99;
  const double cr_Ct101 = -cr_Ct3;
  const double cr_Ct102 = cr_Ct101 + cr_Ct18 + cr_Ct56;
  const double cr_Ct103 = cr_Ct102 * cr_Ct61;
  const double cr_Ct104 = cr_Ct102 * cr_Ct79;
  const double cr_Ct105 = cr_Ct21 + cr_Ct9;
  const double cr_Ct106 = cr_Ct30 + cr_Ct6;
  const double cr_Ct107 = -cr_Ct18;
  const double cr_Ct108 = cr_Ct22 + cr_Ct24;
  const double cr_Ct109 = -cr_Ct20;
  const double cr_Ct110 = cr_Ct28 * std::pow(cr_Ct105 + cr_Ct106 - cr_Ct23 + cr_Ct64, 2) + cr_Ct28 * std::pow(cr_Ct105 + cr_Ct107 + cr_Ct108 + cr_Ct68, 2) + cr_Ct28 * std::pow(cr_Ct106 + cr_Ct108 + cr_Ct109 + cr_Ct72, 2) + cr_Ct35;
  const double cr_Ct111 = cr_Ct110 * cr_Ct38;
  const double cr_Ct112 = std::sqrt(cr_Ct111);
  const double cr_Ct113 = cr_Ct0 * cr_Ct98 / cr_Ct112;
  const double cr_Ct114 = cr_Ct107 + cr_Ct3 + cr_Ct8;
  const double cr_Ct115 = -cr_Ct93;
  const double cr_Ct116 = 0.44444444444444448 * cr_Ct58 / cr_Ct110;
  const double cr_Ct117 = 0.5132002392796674 * cr_Ct38 * threshold / std::pow(cr_Ct111, 3.0 / 2.0);
  const double cr_Ct118 = cr_Ct81 * cr_Ct84 * std::exp(cr_Ct83 * (-cr_Ct112 * cr_Ct82 + 1.0));
  const double cr_Ct119 = cr_Ct102 * cr_Ct99;
  const double cr_Ct120 = cr_Ct101 + cr_Ct20 + cr_Ct55;
  const double cr_Ct121 = cr_Ct120 * cr_Ct61;
  const double cr_Ct122 = cr_Ct120 * cr_Ct79;
  const double cr_Ct123 = cr_Ct109 + cr_Ct3 + cr_Ct5;
  const double cr_Ct124 = cr_Ct120 * cr_Ct99;
  const double cr_Ct125 = 0.055555555555555559 * cr_Ct45;
  const double cr_Ct126 = -cr_Ct125 * cr_Ct31;
  const double cr_Ct127 = -cr_Ct125 * cr_Ct26;
  const double cr_Ct128 = 0.064150029909958425 * cr_Ct78;
  const double cr_Ct129 = cr_Ct128 * cr_Ct75 + cr_Ct60 * (cr_Ct126 + cr_Ct127 + 0.055555555555555559 * cr_Ct53);
  const double cr_Ct130 = cr_Ct85 * r_strain[3];
  const double cr_Ct131 = 0.22222222222222224 * cr_Ct93;
  const double cr_Ct132 = cr_Ct128 * cr_Ct89 + cr_Ct60 * (cr_Ct126 + cr_Ct131 - 0.22222222222222224 * cr_Ct92);
  const double cr_Ct133 = cr_Ct128 * cr_Ct95 + cr_Ct60 * (cr_Ct127 + cr_Ct131 - 0.22222222222222224 * cr_Ct96);
  const double cr_Ct134 = 0.57735026918962584 * cr_Ct41;
  const double cr_Ct135 = 0.5 * cr_Ct60;
  const double cr_Ct136 = 0.57735026918962584 * cr_Ct77;
  const double cr_Ct137 = cr_Ct135 + cr_Ct136;
  const double cr_Ct138 = cr_Ct130 * cr_Ct137;
  const double cr_Ct139 = -cr_Ct138 * r_strain[4];
  const double cr_Ct140 = -cr_Ct138 * r_strain[5];
  const double cr_Ct141 = cr_Ct85 * r_strain[4];
  const double cr_Ct142 = -cr_Ct137 * cr_Ct141 * r_strain[5];
  const double cr_Ct143 = cr_Ct85 * r_strain[5];
  r_Ct(0, 0) = cr_Ct86 * (cr_Ct0 * cr_Ct42 - cr_Ct54 * cr_Ct62 - cr_Ct75 * cr_Ct80);
  r_Ct(0, 1) = -cr_Ct86 * (cr_Ct62 * cr_Ct94 + cr_Ct80 * cr_Ct89 + cr_Ct87);
  r_Ct(0, 2) = -cr_Ct86 * (cr_Ct62 * cr_Ct97 + cr_Ct80 * cr_Ct95 + cr_Ct87);
  r_Ct(0, 3) = -cr_Ct100 * r_strain[3];
  r_Ct(0, 4) = -cr_Ct100 * r_strain[4];
  r_Ct(0, 5) = -cr_Ct100 * r_strain[5];
  r_Ct(1, 0) = -cr_Ct86 * (cr_Ct103 * cr_Ct54 + cr_Ct104 * cr_Ct75 + cr_Ct87);
  r_Ct(1, 1) = cr_Ct118 * (cr_Ct113 - cr_Ct114 * cr_Ct116 * (cr_Ct115 + cr_Ct47 + cr_Ct92) + cr_Ct114 * cr_Ct117 * cr_Ct94);
  r_Ct(1, 2) = -cr_Ct86 * (cr_Ct103 * cr_Ct97 + cr_Ct104 * cr_Ct95 + cr_Ct87);
  r_Ct(1, 3) = -cr_Ct119 * r_strain[3];
  r_Ct(1, 4) = -cr_Ct119 * r_strain[4];
  r_Ct(1, 5) = -cr_Ct119 * r_strain[5];
  r_Ct(2, 0) = -cr_Ct86 * (cr_Ct121 * cr_Ct54 + cr_Ct122 * cr_Ct75 + cr_Ct87);
  r_Ct(2, 1) = -cr_Ct86 * (cr_Ct121 * cr_Ct94 + cr_Ct122 * cr_Ct89 + cr_Ct87);
  r_Ct(2, 2) = cr_Ct118 * (cr_Ct113 - cr_Ct116 * cr_Ct123 * (cr_Ct115 + cr_Ct49 + cr_Ct96) + cr_Ct117 * cr_Ct123 * cr_Ct97);
  r_Ct(2, 3) = -cr_Ct124 * r_strain[3];
  r_Ct(2, 4) = -cr_Ct124 * r_strain[4];
  r_Ct(2, 5) = -cr_Ct124 * r_strain[5];
  r_Ct(3, 0) = -cr_Ct129 * cr_Ct130;
  r_Ct(3, 1) = -cr_Ct130 * cr_Ct132;
  r_Ct(3, 2) = -cr_Ct130 * cr_Ct133;
  r_Ct(3, 3) = cr_Ct85 * (cr_Ct134 - cr_Ct135 * cr_Ct32 - cr_Ct136 * cr_Ct32);
  r_Ct(3, 4) = cr_Ct139;
  r_Ct(3, 5) = cr_Ct140;
  r_Ct(4, 0) = -cr_Ct129 * cr_Ct141;
  r_Ct(4, 1) = -cr_Ct132 * cr_Ct141;
  r_Ct(4, 2) = -cr_Ct133 * cr_Ct141;
  r_Ct(4, 3) = cr_Ct139;
  r_Ct(4, 4) = cr_Ct85 * (cr_Ct134 - cr_Ct135 * cr_Ct33 - cr_Ct136 * cr_Ct33);
  r_Ct(4, 5) = cr_Ct142;
  r_Ct(5, 0) = -cr_Ct129 * cr_Ct143;
  r_Ct(5, 1) = -cr_Ct132 * cr_Ct143;
  r_Ct(5, 2) = -cr_Ct133 * cr_Ct143;
  r_Ct(5, 3) = cr_Ct140;
  r_Ct(5, 4) = cr_Ct142;
  r_Ct(5, 5) = cr_Ct85 * (cr_Ct134 - cr_Ct135 * cr_Ct34 - cr_Ct136 * cr_Ct34);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  const auto &r_props = rValues.GetMaterialProperties();
  const double Young = r_props[YOUNG_MODULUS];
  const double nu = r_props[POISSON_RATIO];
  const double Gf = r_props[FRACTURE_ENERGY];
  const double phi = r_props[FRICTION_ANGLE] * Globals::Pi / 180.0;
  const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
  auto &r_Ct = rValues.GetConstitutiveMatrix();
  const auto &r_strain = rValues.GetStrainVector();

  const bool has_symmetric_yield_stress = r_props.Has(YIELD_STRESS);
  const double threshold_compression = has_symmetric_yield_stress ? r_props[YIELD_STRESS] : r_props[YIELD_STRESS_COMPRESSION];
  const double threshold_tension = has_symmetric_yield_stress ? r_props[YIELD_STRESS] : r_props[YIELD_STRESS_TENSION];

  const double cr_Ct0 = nu - 1.0;
  const double cr_Ct1 = nu*r_strain[0];
  const double cr_Ct2 = nu*r_strain[1];
  const double cr_Ct3 = nu*r_strain[2];
  const double cr_Ct4 = cr_Ct0*r_strain[0];
  const double cr_Ct5 = cr_Ct0*r_strain[1];
  const double cr_Ct6 = cr_Ct0*r_strain[2];
  const double cr_Ct7 = 2.0*nu;
  const double cr_Ct8 = 1.0/(cr_Ct7 - 1.0);
  const double cr_Ct9 = std::abs(threshold_compression/threshold_tension);
  const double cr_Ct10 = std::tan(0.5*phi + 0.78539816339744828);
  const double cr_Ct11 = cr_Ct9/std::pow(cr_Ct10, 2);
  const double cr_Ct12 = std::sin(phi);
  const double cr_Ct13 = cr_Ct11 + 1;
  const double cr_Ct14 = cr_Ct12*cr_Ct13;
  const double cr_Ct15 = cr_Ct8*(cr_Ct11 + cr_Ct14 - 1);
  const double cr_Ct16 = nu + 1.0;
  const double cr_Ct17 = 1.0/cr_Ct16;
  const double cr_Ct18 = Young*cr_Ct17;
  const double cr_Ct19 = cr_Ct15*cr_Ct18*(-2*cr_Ct1 - 2*cr_Ct2 - 2*cr_Ct3 + cr_Ct4 + cr_Ct5 + cr_Ct6);
  const double cr_Ct20 = -nu;
  const double cr_Ct21 = std::pow(cr_Ct20 + 0.5, -2);
  const double cr_Ct22 = cr_Ct1;
  const double cr_Ct23 = 0.5*cr_Ct2;
  const double cr_Ct24 = -cr_Ct23;
  const double cr_Ct25 = 0.5*cr_Ct3;
  const double cr_Ct26 = -cr_Ct25;
  const double cr_Ct27 = cr_Ct20 + 1.0;
  const double cr_Ct28 = cr_Ct27*r_strain[0];
  const double cr_Ct29 = cr_Ct27*r_strain[1];
  const double cr_Ct30 = 0.5*cr_Ct29;
  const double cr_Ct31 = cr_Ct27*r_strain[2];
  const double cr_Ct32 = 0.5*cr_Ct31;
  const double cr_Ct33 = cr_Ct22 + cr_Ct24 + cr_Ct26 - cr_Ct28 + cr_Ct30 + cr_Ct32;
  const double cr_Ct34 = cr_Ct21*std::pow(cr_Ct33, 2);
  const double cr_Ct35 = cr_Ct2;
  const double cr_Ct36 = 0.5*cr_Ct1;
  const double cr_Ct37 = -cr_Ct36 - 0.5*cr_Ct4;
  const double cr_Ct38 = cr_Ct26 + cr_Ct35 + cr_Ct37 + cr_Ct5 - 0.5*cr_Ct6;
  const double cr_Ct39 = std::pow(cr_Ct38, 2);
  const double cr_Ct40 = std::pow(nu - 0.5, -2);
  const double cr_Ct41 = 0.22222222222222227*cr_Ct40;
  const double cr_Ct42 = cr_Ct3;
  const double cr_Ct43 = cr_Ct24 + cr_Ct37 + cr_Ct42 - 0.5*cr_Ct5 + cr_Ct6;
  const double cr_Ct44 = std::pow(cr_Ct43, 2);
  const double cr_Ct45 = std::pow(r_strain[3], 2);
  const double cr_Ct46 = std::pow(r_strain[4], 2);
  const double cr_Ct47 = std::pow(r_strain[5], 2);
  const double cr_Ct48 = cr_Ct45 + cr_Ct46 + cr_Ct47;
  const double cr_Ct49 = std::pow(Young, 2);
  const double cr_Ct50 = std::pow(cr_Ct16, -2);
  const double cr_Ct51 = cr_Ct49*cr_Ct50;
  const double cr_Ct52 = std::sqrt(cr_Ct51*(0.22222222222222227*cr_Ct34 + cr_Ct39*cr_Ct41 + cr_Ct41*cr_Ct44 + cr_Ct48));
  const double cr_Ct53 = 1.0/cr_Ct52;
  const double cr_Ct54 = 0.11111111111111113*cr_Ct40;
  const double cr_Ct55 = 0.5*cr_Ct45 + 0.5*cr_Ct46 + 0.5*cr_Ct47;
  const double cr_Ct56 = 0.25*r_strain[4];
  const double cr_Ct57 = 1.0/(1.0 - cr_Ct7);
  const double cr_Ct58 = 0.5*cr_Ct57;
  const double cr_Ct59 = 0.33333333333333337*cr_Ct1 - 0.33333333333333331*cr_Ct28;
  const double cr_Ct60 = 0.33333333333333337*cr_Ct2 - 0.33333333333333331*cr_Ct29;
  const double cr_Ct61 = -0.66666666666666663*cr_Ct3 + 0.66666666666666674*cr_Ct31 + cr_Ct59 + cr_Ct60;
  const double cr_Ct62 = cr_Ct61*r_strain[3];
  const double cr_Ct63 = r_strain[3]*(cr_Ct56*r_strain[5] - cr_Ct58*cr_Ct62);
  const double cr_Ct64 = 0.33333333333333337*cr_Ct3 - 0.33333333333333331*cr_Ct31;
  const double cr_Ct65 = -0.66666666666666663*cr_Ct2 + 0.66666666666666674*cr_Ct29 + cr_Ct59 + cr_Ct64;
  const double cr_Ct66 = cr_Ct65*r_strain[5];
  const double cr_Ct67 = r_strain[5]*(cr_Ct56*r_strain[3] - cr_Ct58*cr_Ct66);
  const double cr_Ct68 = 0.25*cr_Ct21;
  const double cr_Ct69 = cr_Ct61*cr_Ct65;
  const double cr_Ct70 = -0.66666666666666663*cr_Ct1 + 0.66666666666666674*cr_Ct28 + cr_Ct60 + cr_Ct64;
  const double cr_Ct71 = cr_Ct57*cr_Ct70;
  const double cr_Ct72 = cr_Ct71*(-0.25*cr_Ct46 + cr_Ct68*cr_Ct69);
  const double cr_Ct73 = 2.598076211353316*cr_Ct63 + 2.598076211353316*cr_Ct67 + 5.196152422706632*cr_Ct72;
  const double cr_Ct74 = 2.0*cr_Ct73;
  const double cr_Ct75 = cr_Ct18*cr_Ct74;
  const double cr_Ct76 = 0.33333333333333331*std::asin(cr_Ct53*cr_Ct75/(0.11111111111111113*cr_Ct34 + cr_Ct39*cr_Ct54 + cr_Ct44*cr_Ct54 + cr_Ct55));
  const double cr_Ct77 = std::cos(cr_Ct76);
  const double cr_Ct78 = 1 - cr_Ct11;
  const double cr_Ct79 = -cr_Ct12*cr_Ct78 + cr_Ct13;
  const double cr_Ct80 = 0.5*cr_Ct79;
  const double cr_Ct81 = std::sin(cr_Ct76);
  const double cr_Ct82 = cr_Ct12*(cr_Ct13 - cr_Ct78/cr_Ct12);
  const double cr_Ct83 = 0.28867513459481292*cr_Ct82;
  const double cr_Ct84 = cr_Ct77*cr_Ct80 + cr_Ct81*cr_Ct83;
  const double cr_Ct85 = cr_Ct52*cr_Ct84;
  const double cr_Ct86 = 0.16666666666666666*cr_Ct19 + 0.5*cr_Ct85;
  const double cr_Ct87 = 1.0/cr_Ct86;
  const double cr_Ct88 = std::cos(phi);
  const double cr_Ct89 = cr_Ct88*threshold_compression/cr_Ct10;
  const double cr_Ct90 = cr_Ct87*cr_Ct89;
  const double cr_Ct91 = 0.5*cr_Ct90;
  const double cr_Ct92 = cr_Ct0*cr_Ct91;
  const double cr_Ct93 = cr_Ct57*(0.16666666666666666*cr_Ct11 + 0.16666666666666666*cr_Ct14 - 0.16666666666666666);
  const double cr_Ct94 = 4.0*nu;
  const double cr_Ct95 = 2.0 - cr_Ct94;
  const double cr_Ct96 = cr_Ct25 - cr_Ct32;
  const double cr_Ct97 = cr_Ct23 - cr_Ct30;
  const double cr_Ct98 = -cr_Ct22 + cr_Ct28 + cr_Ct96 + cr_Ct97;
  const double cr_Ct99 = cr_Ct7 - 1.0;
  const double cr_Ct100 = -0.5*cr_Ct28 + cr_Ct36;
  const double cr_Ct101 = cr_Ct100 + cr_Ct29 - cr_Ct35 + cr_Ct96;
  const double cr_Ct102 = cr_Ct101*cr_Ct99;
  const double cr_Ct103 = cr_Ct100 + cr_Ct31 - cr_Ct42 + cr_Ct97;
  const double cr_Ct104 = cr_Ct103*cr_Ct99;
  const double cr_Ct105 = cr_Ct102 + cr_Ct104 + cr_Ct95*cr_Ct98;
  const double cr_Ct106 = cr_Ct105*cr_Ct21;
  const double cr_Ct107 = std::pow(cr_Ct98, 2);
  const double cr_Ct108 = 0.22222222222222227*cr_Ct21;
  const double cr_Ct109 = std::pow(cr_Ct101, 2);
  const double cr_Ct110 = std::pow(cr_Ct103, 2);
  const double cr_Ct111 = cr_Ct107*cr_Ct108 + cr_Ct108*cr_Ct109 + cr_Ct108*cr_Ct110 + cr_Ct48;
  const double cr_Ct112 = std::pow(cr_Ct111*cr_Ct51, -1.0/2.0);
  const double cr_Ct113 = 0.11111111111111113*cr_Ct21;
  const double cr_Ct114 = 1.0/(cr_Ct107*cr_Ct113 + cr_Ct109*cr_Ct113 + cr_Ct110*cr_Ct113 + cr_Ct55);
  const double cr_Ct115 = 0.33333333333333331*std::asin(cr_Ct112*cr_Ct114*cr_Ct75);
  const double cr_Ct116 = std::cos(cr_Ct115);
  const double cr_Ct117 = std::sin(cr_Ct115);
  const double cr_Ct118 = cr_Ct112*(cr_Ct116*cr_Ct80 + cr_Ct117*cr_Ct83);
  const double cr_Ct119 = 0.055555555555555566*Young*cr_Ct118*cr_Ct50;
  const double cr_Ct120 = 0.66666666666666674*nu;
  const double cr_Ct121 = cr_Ct120 - 0.33333333333333331;
  const double cr_Ct122 = 2.598076211353316*cr_Ct121;
  const double cr_Ct123 = cr_Ct122*cr_Ct45;
  const double cr_Ct124 = cr_Ct122*cr_Ct47;
  const double cr_Ct125 = cr_Ct21*cr_Ct70;
  const double cr_Ct126 = 0.66666666666666674 - 1.3333333333333335*nu;
  const double cr_Ct127 = -2.598076211353316*cr_Ct21*cr_Ct69 + 2.598076211353316*cr_Ct46;
  const double cr_Ct128 = cr_Ct114*cr_Ct57;
  const double cr_Ct129 = cr_Ct73/std::pow(cr_Ct111, 2);
  const double cr_Ct130 = 0.88888888888888906*cr_Ct129;
  const double cr_Ct131 = cr_Ct114/cr_Ct111;
  const double cr_Ct132 = cr_Ct108*cr_Ct131*cr_Ct73;
  const double cr_Ct133 = cr_Ct105*cr_Ct132 + cr_Ct106*cr_Ct130 + cr_Ct128*(-cr_Ct122*cr_Ct125*(cr_Ct120*r_strain[0] - 0.33333333333333326*cr_Ct2 - 0.66666666666666663*cr_Ct28 + 0.33333333333333343*cr_Ct29 - 0.33333333333333326*cr_Ct3 + 0.33333333333333343*cr_Ct31) + cr_Ct123 + cr_Ct124 + cr_Ct126*cr_Ct127);
  const double cr_Ct134 = 0.5*cr_Ct133;
  const double cr_Ct135 = std::pow(0.0023148148148148147 - std::pow(0.5*cr_Ct63 + 0.5*cr_Ct67 + cr_Ct72, 2)/std::pow(cr_Ct111, 3), -1.0/2.0);
  const double cr_Ct136 = 0.0080187537387448014*cr_Ct79;
  const double cr_Ct137 = 0.0046296296296296294*cr_Ct82;
  const double cr_Ct138 = cr_Ct135*(-cr_Ct116*cr_Ct137 + cr_Ct117*cr_Ct136);
  const double cr_Ct139 = cr_Ct138*cr_Ct17;
  const double cr_Ct140 = cr_Ct106*cr_Ct119 + cr_Ct134*cr_Ct139 + cr_Ct93;
  const double cr_Ct141 = -cr_Ct2;
  const double cr_Ct142 = -cr_Ct3;
  const double cr_Ct143 = cr_Ct141 + cr_Ct142 + cr_Ct4;
  const double cr_Ct144 = cr_Ct89/std::pow(0.33333333333333331*cr_Ct19 + cr_Ct85, 2);
  const double cr_Ct145 = 2.0*Young*cr_Ct144;
  const double cr_Ct146 = cr_Ct143*cr_Ct145;
  const double cr_Ct147 = cr_Ct15*(cr_Ct20 - 1.0);
  const double cr_Ct148 = 0.16666666666666666*cr_Ct147;
  const double cr_Ct149 = 0.25*cr_Ct40*cr_Ct99;
  const double cr_Ct150 = -cr_Ct149*cr_Ct43;
  const double cr_Ct151 = -cr_Ct149*cr_Ct38;
  const double cr_Ct152 = cr_Ct94 - 2.0;
  const double cr_Ct153 = cr_Ct150 + cr_Ct151 + cr_Ct152*cr_Ct33*cr_Ct68;
  const double cr_Ct154 = cr_Ct53*cr_Ct84;
  const double cr_Ct155 = cr_Ct154*cr_Ct18;
  const double cr_Ct156 = 0.22222222222222227*cr_Ct155;
  const double cr_Ct157 = cr_Ct112*cr_Ct135*cr_Ct52*(cr_Ct136*cr_Ct81 - cr_Ct137*cr_Ct77);
  const double cr_Ct158 = cr_Ct134*cr_Ct157 + cr_Ct148 + cr_Ct153*cr_Ct156;
  const double cr_Ct159 = 1.0/(Gf*Young*std::pow(cr_Ct9, 2)/(characteristic_length*std::pow(threshold_compression, 2)) - 0.5);
  const double cr_Ct160 = cr_Ct159;
  const double cr_Ct161 = cr_Ct160*cr_Ct18*cr_Ct87;
  const double cr_Ct162 = cr_Ct143*cr_Ct161;
  const double cr_Ct163 = std::exp(cr_Ct160*(-2.0*cr_Ct10*cr_Ct86/(cr_Ct88*threshold_compression) + 1.0));
  const double cr_Ct164 = cr_Ct163*cr_Ct18;
  const double cr_Ct165 = cr_Ct164*cr_Ct8;
  const double cr_Ct166 = cr_Ct91*nu;
  const double cr_Ct167 = cr_Ct98*cr_Ct99;
  const double cr_Ct168 = cr_Ct101*cr_Ct95 + cr_Ct104 + cr_Ct167;
  const double cr_Ct169 = cr_Ct119*cr_Ct21;
  const double cr_Ct170 = 2.598076211353316*cr_Ct126;
  const double cr_Ct171 = cr_Ct121*cr_Ct127;
  const double cr_Ct172 = 2.598076211353316*cr_Ct125;
  const double cr_Ct173 = cr_Ct130*cr_Ct21;
  const double cr_Ct174 = cr_Ct128*(cr_Ct123 + cr_Ct170*cr_Ct47 + cr_Ct171 - cr_Ct172*(cr_Ct121*cr_Ct65 + cr_Ct126*cr_Ct61)) + cr_Ct132*cr_Ct168 + cr_Ct168*cr_Ct173;
  const double cr_Ct175 = 0.5*cr_Ct139;
  const double cr_Ct176 = cr_Ct168*cr_Ct169 + cr_Ct174*cr_Ct175 + cr_Ct93;
  const double cr_Ct177 = cr_Ct57*cr_Ct8;
  const double cr_Ct178 = cr_Ct152*cr_Ct177;
  const double cr_Ct179 = cr_Ct177*cr_Ct33*cr_Ct99;
  const double cr_Ct180 = cr_Ct150 - cr_Ct178*cr_Ct38 + cr_Ct179;
  const double cr_Ct181 = 0.5*cr_Ct157;
  const double cr_Ct182 = cr_Ct148 + cr_Ct156*cr_Ct180 + cr_Ct174*cr_Ct181;
  const double cr_Ct183 = cr_Ct102 + cr_Ct103*cr_Ct95 + cr_Ct167;
  const double cr_Ct184 = cr_Ct128*(cr_Ct124 + cr_Ct170*cr_Ct45 + cr_Ct171 - cr_Ct172*(cr_Ct121*cr_Ct61 + cr_Ct126*cr_Ct65)) + cr_Ct132*cr_Ct183 + cr_Ct173*cr_Ct183;
  const double cr_Ct185 = cr_Ct169*cr_Ct183 + cr_Ct175*cr_Ct184 + cr_Ct93;
  const double cr_Ct186 = cr_Ct151 - cr_Ct178*cr_Ct43 + cr_Ct179;
  const double cr_Ct187 = cr_Ct148 + cr_Ct156*cr_Ct186 + cr_Ct181*cr_Ct184;
  const double cr_Ct188 = cr_Ct18*r_strain[3];
  const double cr_Ct189 = 2.598076211353316*r_strain[5];
  const double cr_Ct190 = 5.196152422706632*cr_Ct57;
  const double cr_Ct191 = 8.0*cr_Ct129;
  const double cr_Ct192 = cr_Ct131*cr_Ct74;
  const double cr_Ct193 = cr_Ct114*(-cr_Ct189*r_strain[4] + cr_Ct190*cr_Ct62) + cr_Ct191*r_strain[3] + cr_Ct192*r_strain[3];
  const double cr_Ct194 = cr_Ct118*cr_Ct188 + cr_Ct138*cr_Ct193;
  const double cr_Ct195 = cr_Ct144;
  const double cr_Ct196 = cr_Ct154*cr_Ct188;
  const double cr_Ct197 = cr_Ct157*cr_Ct193;
  const double cr_Ct198 = cr_Ct159*cr_Ct87;
  const double cr_Ct199 = cr_Ct194*cr_Ct195 + cr_Ct198*(0.5*cr_Ct196 + 0.5*cr_Ct197);
  const double cr_Ct200 = cr_Ct163*cr_Ct51*cr_Ct8;
  const double cr_Ct201 = cr_Ct143*cr_Ct200;
  const double cr_Ct202 = cr_Ct118*cr_Ct18;
  const double cr_Ct203 = cr_Ct114*(-cr_Ct189*r_strain[3] + 5.196152422706632*cr_Ct71*r_strain[4]) + cr_Ct191*r_strain[4] + cr_Ct192*r_strain[4];
  const double cr_Ct204 = cr_Ct138*cr_Ct203 + cr_Ct202*r_strain[4];
  const double cr_Ct205 = cr_Ct155*r_strain[4];
  const double cr_Ct206 = cr_Ct157*cr_Ct203;
  const double cr_Ct207 = cr_Ct195*cr_Ct204 + cr_Ct198*(0.5*cr_Ct205 + 0.5*cr_Ct206);
  const double cr_Ct208 = cr_Ct114*(cr_Ct190*cr_Ct66 - 2.598076211353316*r_strain[3]*r_strain[4]) + cr_Ct191*r_strain[5] + cr_Ct192*r_strain[5];
  const double cr_Ct209 = cr_Ct138*cr_Ct208 + cr_Ct202*r_strain[5];
  const double cr_Ct210 = cr_Ct155*r_strain[5];
  const double cr_Ct211 = cr_Ct157*cr_Ct208;
  const double cr_Ct212 = cr_Ct195*cr_Ct209 + cr_Ct198*(0.5*cr_Ct210 + 0.5*cr_Ct211);
  const double cr_Ct213 = -cr_Ct1;
  const double cr_Ct214 = cr_Ct142 + cr_Ct213 + cr_Ct5;
  const double cr_Ct215 = cr_Ct145*cr_Ct214;
  const double cr_Ct216 = cr_Ct161*cr_Ct214;
  const double cr_Ct217 = cr_Ct200*cr_Ct214;
  const double cr_Ct218 = cr_Ct141 + cr_Ct213 + cr_Ct6;
  const double cr_Ct219 = cr_Ct145*cr_Ct218;
  const double cr_Ct220 = cr_Ct161*cr_Ct218;
  const double cr_Ct221 = cr_Ct200*cr_Ct218;
  const double cr_Ct222 = 0.083333333333333329*cr_Ct147;
  const double cr_Ct223 = 0.11111111111111113*cr_Ct155;
  const double cr_Ct224 = 0.25*cr_Ct157;
  const double cr_Ct225 = cr_Ct17*cr_Ct198;
  const double cr_Ct226 = cr_Ct140*cr_Ct195 + cr_Ct225*(cr_Ct133*cr_Ct224 + cr_Ct153*cr_Ct223 + cr_Ct222);
  const double cr_Ct227 = cr_Ct163*r_strain[3];
  const double cr_Ct228 = cr_Ct17*cr_Ct49;
  const double cr_Ct229 = cr_Ct227*cr_Ct228;
  const double cr_Ct230 = cr_Ct176*cr_Ct195 + cr_Ct225*(cr_Ct174*cr_Ct224 + cr_Ct180*cr_Ct223 + cr_Ct222);
  const double cr_Ct231 = cr_Ct185*cr_Ct195 + cr_Ct225*(cr_Ct184*cr_Ct224 + cr_Ct186*cr_Ct223 + cr_Ct222);
  const double cr_Ct232 = 0.25*cr_Ct90;
  const double cr_Ct233 = 0.5*cr_Ct144;
  const double cr_Ct234 = cr_Ct194*cr_Ct233;
  const double cr_Ct235 = Young*cr_Ct225;
  const double cr_Ct236 = 0.25*cr_Ct235;
  const double cr_Ct237 = cr_Ct204*cr_Ct233;
  const double cr_Ct238 = cr_Ct198*(cr_Ct155*cr_Ct56 + 0.25*cr_Ct206) + cr_Ct237;
  const double cr_Ct239 = cr_Ct227*cr_Ct51;
  const double cr_Ct240 = cr_Ct209*cr_Ct233;
  const double cr_Ct241 = cr_Ct198*(0.25*cr_Ct210 + 0.25*cr_Ct211) + cr_Ct240;
  const double cr_Ct242 = cr_Ct163*r_strain[4];
  const double cr_Ct243 = cr_Ct228*cr_Ct242;
  const double cr_Ct244 = cr_Ct198*(0.25*cr_Ct196 + 0.25*cr_Ct197) + cr_Ct234;
  const double cr_Ct245 = cr_Ct242*cr_Ct51;
  const double cr_Ct246 = cr_Ct163*r_strain[5];
  const double cr_Ct247 = cr_Ct228*cr_Ct246;
  const double cr_Ct248 = cr_Ct246*cr_Ct51;
  r_Ct(0,0)=cr_Ct165*(-cr_Ct140*cr_Ct146 - cr_Ct158*cr_Ct162 + cr_Ct92);
  r_Ct(0,1)=-cr_Ct165*(cr_Ct146*cr_Ct176 + cr_Ct162*cr_Ct182 + cr_Ct166);
  r_Ct(0,2)=-cr_Ct165*(cr_Ct146*cr_Ct185 + cr_Ct162*cr_Ct187 + cr_Ct166);
  r_Ct(0,3)=-cr_Ct199*cr_Ct201;
  r_Ct(0,4)=-cr_Ct201*cr_Ct207;
  r_Ct(0,5)=-cr_Ct201*cr_Ct212;
  r_Ct(1,0)=-cr_Ct165*(cr_Ct140*cr_Ct215 + cr_Ct158*cr_Ct216 + cr_Ct166);
  r_Ct(1,1)=cr_Ct165*(-cr_Ct176*cr_Ct215 - cr_Ct182*cr_Ct216 + cr_Ct92);
  r_Ct(1,2)=-cr_Ct165*(cr_Ct166 + cr_Ct185*cr_Ct215 + cr_Ct187*cr_Ct216);
  r_Ct(1,3)=-cr_Ct199*cr_Ct217;
  r_Ct(1,4)=-cr_Ct207*cr_Ct217;
  r_Ct(1,5)=-cr_Ct212*cr_Ct217;
  r_Ct(2,0)=-cr_Ct165*(cr_Ct140*cr_Ct219 + cr_Ct158*cr_Ct220 + cr_Ct166);
  r_Ct(2,1)=-cr_Ct165*(cr_Ct166 + cr_Ct176*cr_Ct219 + cr_Ct182*cr_Ct220);
  r_Ct(2,2)=cr_Ct165*(-cr_Ct185*cr_Ct219 - cr_Ct187*cr_Ct220 + cr_Ct92);
  r_Ct(2,3)=-cr_Ct199*cr_Ct221;
  r_Ct(2,4)=-cr_Ct207*cr_Ct221;
  r_Ct(2,5)=-cr_Ct212*cr_Ct221;
  r_Ct(3,0)=-cr_Ct226*cr_Ct229;
  r_Ct(3,1)=-cr_Ct229*cr_Ct230;
  r_Ct(3,2)=-cr_Ct229*cr_Ct231;
  r_Ct(3,3)=cr_Ct164*(-cr_Ct188*cr_Ct234 + cr_Ct232 - cr_Ct236*r_strain[3]*(cr_Ct196 + cr_Ct197));
  r_Ct(3,4)=-cr_Ct238*cr_Ct239;
  r_Ct(3,5)=-cr_Ct239*cr_Ct241;
  r_Ct(4,0)=-cr_Ct226*cr_Ct243;
  r_Ct(4,1)=-cr_Ct230*cr_Ct243;
  r_Ct(4,2)=-cr_Ct231*cr_Ct243;
  r_Ct(4,3)=-cr_Ct244*cr_Ct245;
  r_Ct(4,4)=cr_Ct164*(-cr_Ct18*cr_Ct237*r_strain[4] + cr_Ct232 - cr_Ct235*cr_Ct56*(cr_Ct205 + cr_Ct206));
  r_Ct(4,5)=-cr_Ct241*cr_Ct245;
  r_Ct(5,0)=-cr_Ct226*cr_Ct247;
  r_Ct(5,1)=-cr_Ct230*cr_Ct247;
  r_Ct(5,2)=-cr_Ct231*cr_Ct247;
  r_Ct(5,3)=-cr_Ct244*cr_Ct248;
  r_Ct(5,4)=-cr_Ct238*cr_Ct248;
  r_Ct(5,5)=cr_Ct164*(-cr_Ct18*cr_Ct240*r_strain[5] + cr_Ct232 - cr_Ct236*r_strain[5]*(cr_Ct210 + cr_Ct211));
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<TrescaYieldSurface<TrescaPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  const auto &r_props = rValues.GetMaterialProperties();
  const double Young = r_props[YOUNG_MODULUS];
  const double nu = r_props[POISSON_RATIO];
  const double Gf = r_props[FRACTURE_ENERGY];
  const double phi = r_props[FRICTION_ANGLE];
  const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
  const double threshold = r_props[YIELD_STRESS];
  auto &r_Ct = rValues.GetConstitutiveMatrix();
  const auto &r_strain = rValues.GetStrainVector();
  const double sin_phi = std::sin(phi * Globals::Pi / 180.0);

  const double cr_Ct0 = nu - 1.0;
  const double cr_Ct1 = 1.7320508075688772*sin_phi;
  const double cr_Ct2 = cr_Ct1 - 5.196152422706632;
  const double cr_Ct3 = 1.0/cr_Ct2;
  const double cr_Ct4 = sin_phi - 1;
  const double cr_Ct5 = 1.0/cr_Ct4;
  const double cr_Ct6 = std::abs(cr_Ct5*threshold*(sin_phi + 3.0));
  const double cr_Ct7 = cr_Ct3*cr_Ct4*cr_Ct6;
  const double cr_Ct8 = cr_Ct0*cr_Ct7;
  const double cr_Ct9 = 1.0/(Gf*Young/(characteristic_length*std::pow(threshold, 2)) - 0.5);
  const double cr_Ct10 = 2.0*nu;
  const double cr_Ct11 = -cr_Ct10;
  const double cr_Ct12 = 1.0/(cr_Ct10 - 1.0);
  const double cr_Ct13 = cr_Ct12*cr_Ct3*sin_phi;
  const double cr_Ct14 = -cr_Ct13*(cr_Ct11 - 2.0);
  const double cr_Ct15 = cr_Ct0*r_strain[2];
  const double cr_Ct16 = nu*r_strain[1];
  const double cr_Ct17 = 0.5*cr_Ct16;
  const double cr_Ct18 = -cr_Ct17;
  const double cr_Ct19 = nu*r_strain[2];
  const double cr_Ct20 = cr_Ct19;
  const double cr_Ct21 = cr_Ct0*r_strain[1];
  const double cr_Ct22 = nu*r_strain[0];
  const double cr_Ct23 = 0.5*cr_Ct22;
  const double cr_Ct24 = cr_Ct0*r_strain[0];
  const double cr_Ct25 = -cr_Ct23 - 0.5*cr_Ct24;
  const double cr_Ct26 = cr_Ct15 + cr_Ct18 + cr_Ct20 - 0.5*cr_Ct21 + cr_Ct25;
  const double cr_Ct27 = cr_Ct10 - 1.0;
  const double cr_Ct28 = std::pow(nu - 0.5, -2);
  const double cr_Ct29 = 0.25*cr_Ct27*cr_Ct28;
  const double cr_Ct30 = -cr_Ct26*cr_Ct29;
  const double cr_Ct31 = 0.5*cr_Ct19;
  const double cr_Ct32 = -cr_Ct31;
  const double cr_Ct33 = cr_Ct16;
  const double cr_Ct34 = -0.5*cr_Ct15 + cr_Ct21 + cr_Ct25 + cr_Ct32 + cr_Ct33;
  const double cr_Ct35 = -cr_Ct29*cr_Ct34;
  const double cr_Ct36 = 4.0*nu;
  const double cr_Ct37 = cr_Ct36 - 2.0;
  const double cr_Ct38 = -nu;
  const double cr_Ct39 = std::pow(cr_Ct38 + 0.5, -2);
  const double cr_Ct40 = cr_Ct22;
  const double cr_Ct41 = cr_Ct38 + 1.0;
  const double cr_Ct42 = cr_Ct41*r_strain[0];
  const double cr_Ct43 = cr_Ct41*r_strain[1];
  const double cr_Ct44 = 0.5*cr_Ct43;
  const double cr_Ct45 = cr_Ct41*r_strain[2];
  const double cr_Ct46 = 0.5*cr_Ct45;
  const double cr_Ct47 = cr_Ct18 + cr_Ct32 + cr_Ct40 - cr_Ct42 + cr_Ct44 + cr_Ct46;
  const double cr_Ct48 = 0.22222222222222227*cr_Ct39;
  const double cr_Ct49 = 0.22222222222222227*cr_Ct28;
  const double cr_Ct50 = std::pow(r_strain[3], 2);
  const double cr_Ct51 = std::pow(r_strain[4], 2);
  const double cr_Ct52 = std::pow(r_strain[5], 2);
  const double cr_Ct53 = cr_Ct50 + cr_Ct51 + cr_Ct52;
  const double cr_Ct54 = nu + 1.0;
  const double cr_Ct55 = std::pow(Young, 2)/std::pow(cr_Ct54, 2);
  const double cr_Ct56 = std::sqrt(cr_Ct55*(std::pow(cr_Ct26, 2)*cr_Ct49 + std::pow(cr_Ct34, 2)*cr_Ct49 + std::pow(cr_Ct47, 2)*cr_Ct48 + cr_Ct53));
  const double cr_Ct57 = 1.0/cr_Ct56;
  const double cr_Ct58 = Young/cr_Ct54;
  const double cr_Ct59 = 0.22222222222222227*cr_Ct57*cr_Ct58;
  const double cr_Ct60 = cr_Ct9*(cr_Ct14 + cr_Ct59*(cr_Ct30 + cr_Ct35 + 0.25*cr_Ct37*cr_Ct39*cr_Ct47));
  const double cr_Ct61 = -cr_Ct16;
  const double cr_Ct62 = -cr_Ct19;
  const double cr_Ct63 = cr_Ct24 + cr_Ct61 + cr_Ct62;
  const double cr_Ct64 = cr_Ct58*cr_Ct63;
  const double cr_Ct65 = 1.0/(cr_Ct11 + 1.0);
  const double cr_Ct66 = cr_Ct65*sin_phi/(5.196152422706632 - cr_Ct1);
  const double cr_Ct67 = cr_Ct66*(cr_Ct10 + 2.0);
  const double cr_Ct68 = 2.0 - cr_Ct36;
  const double cr_Ct69 = cr_Ct31 - cr_Ct46;
  const double cr_Ct70 = cr_Ct17 - cr_Ct44;
  const double cr_Ct71 = -cr_Ct40 + cr_Ct42 + cr_Ct69 + cr_Ct70;
  const double cr_Ct72 = cr_Ct23 - 0.5*cr_Ct42;
  const double cr_Ct73 = -cr_Ct33 + cr_Ct43 + cr_Ct69 + cr_Ct72;
  const double cr_Ct74 = cr_Ct27*cr_Ct73;
  const double cr_Ct75 = -cr_Ct20 + cr_Ct45 + cr_Ct70 + cr_Ct72;
  const double cr_Ct76 = cr_Ct27*cr_Ct75;
  const double cr_Ct77 = std::sqrt(cr_Ct55*(cr_Ct48*std::pow(cr_Ct71, 2) + cr_Ct48*std::pow(cr_Ct73, 2) + cr_Ct48*std::pow(cr_Ct75, 2) + cr_Ct53));
  const double cr_Ct78 = 0.055555555555555566*cr_Ct39*cr_Ct58/cr_Ct77;
  const double cr_Ct79 = cr_Ct36*r_strain[0];
  const double cr_Ct80 = cr_Ct36*r_strain[1];
  const double cr_Ct81 = cr_Ct36*r_strain[2];
  const double cr_Ct82 = 1.0/(-cr_Ct13*cr_Ct58*(2.0*cr_Ct15 + 2.0*cr_Ct21 + 2.0*cr_Ct24 - cr_Ct79 - cr_Ct80 - cr_Ct81) + 0.5*cr_Ct56);
  const double cr_Ct83 = cr_Ct7*cr_Ct82;
  const double cr_Ct84 = cr_Ct83*(cr_Ct67 + cr_Ct78*(cr_Ct68*cr_Ct71 + cr_Ct74 + cr_Ct76));
  const double cr_Ct85 = cr_Ct82*std::exp(-cr_Ct9*(cr_Ct2*cr_Ct5*(cr_Ct58*cr_Ct66*(2.0*cr_Ct42 + 2.0*cr_Ct43 + 2.0*cr_Ct45 + cr_Ct79 + cr_Ct80 + cr_Ct81) + 0.5*cr_Ct77)/cr_Ct6 - 1));
  const double cr_Ct86 = cr_Ct58*cr_Ct85;
  const double cr_Ct87 = cr_Ct12*cr_Ct86;
  const double cr_Ct88 = cr_Ct7*nu;
  const double cr_Ct89 = cr_Ct12*cr_Ct65;
  const double cr_Ct90 = cr_Ct37*cr_Ct89;
  const double cr_Ct91 = cr_Ct27*cr_Ct47*cr_Ct89;
  const double cr_Ct92 = cr_Ct9*(cr_Ct14 + cr_Ct59*(cr_Ct30 - cr_Ct34*cr_Ct90 + cr_Ct91));
  const double cr_Ct93 = cr_Ct27*cr_Ct71;
  const double cr_Ct94 = cr_Ct83*(cr_Ct67 + cr_Ct78*(cr_Ct68*cr_Ct73 + cr_Ct76 + cr_Ct93));
  const double cr_Ct95 = cr_Ct9*(cr_Ct14 + cr_Ct59*(-cr_Ct26*cr_Ct90 + cr_Ct35 + cr_Ct91));
  const double cr_Ct96 = cr_Ct83*(cr_Ct67 + cr_Ct78*(cr_Ct68*cr_Ct75 + cr_Ct74 + cr_Ct93));
  const double cr_Ct97 = cr_Ct85*r_strain[3];
  const double cr_Ct98 = std::pow(Young, 3)*(cr_Ct83 + cr_Ct9)/std::pow(cr_Ct54, 3);
  const double cr_Ct99 = cr_Ct97*cr_Ct98;
  const double cr_Ct100 = 0.5*cr_Ct12*cr_Ct57;
  const double cr_Ct101 = cr_Ct100*cr_Ct63;
  const double cr_Ct102 = cr_Ct85*cr_Ct98*r_strain[4];
  const double cr_Ct103 = cr_Ct85*cr_Ct98*r_strain[5];
  const double cr_Ct104 = -cr_Ct22;
  const double cr_Ct105 = cr_Ct104 + cr_Ct21 + cr_Ct62;
  const double cr_Ct106 = cr_Ct105*cr_Ct58;
  const double cr_Ct107 = cr_Ct100*cr_Ct105;
  const double cr_Ct108 = cr_Ct104 + cr_Ct15 + cr_Ct61;
  const double cr_Ct109 = cr_Ct108*cr_Ct58;
  const double cr_Ct110 = cr_Ct100*cr_Ct108;
  const double cr_Ct111 = cr_Ct60 + cr_Ct84;
  const double cr_Ct112 = 0.5*cr_Ct55;
  const double cr_Ct113 = cr_Ct112*cr_Ct97;
  const double cr_Ct114 = cr_Ct92 + cr_Ct94;
  const double cr_Ct115 = cr_Ct95 + cr_Ct96;
  const double cr_Ct116 = 0.5*cr_Ct7;
  const double cr_Ct117 = 0.25*cr_Ct57;
  const double cr_Ct118 = cr_Ct117*cr_Ct55;
  const double cr_Ct119 = cr_Ct118*cr_Ct9;
  const double cr_Ct120 = cr_Ct118*cr_Ct83;
  const double cr_Ct121 = cr_Ct117*cr_Ct99;
  const double cr_Ct122 = -cr_Ct121*r_strain[4];
  const double cr_Ct123 = -cr_Ct121*r_strain[5];
  const double cr_Ct124 = cr_Ct112*cr_Ct85;
  const double cr_Ct125 = cr_Ct124*r_strain[4];
  const double cr_Ct126 = -cr_Ct102*cr_Ct117*r_strain[5];
  const double cr_Ct127 = cr_Ct124*r_strain[5];
  r_Ct(0,0)=cr_Ct87*(-cr_Ct60*cr_Ct64 - cr_Ct64*cr_Ct84 + cr_Ct8);
  r_Ct(0,1)=-cr_Ct87*(cr_Ct64*cr_Ct92 + cr_Ct64*cr_Ct94 + cr_Ct88);
  r_Ct(0,2)=-cr_Ct87*(cr_Ct64*cr_Ct95 + cr_Ct64*cr_Ct96 + cr_Ct88);
  r_Ct(0,3)=-cr_Ct101*cr_Ct99;
  r_Ct(0,4)=-cr_Ct101*cr_Ct102;
  r_Ct(0,5)=-cr_Ct101*cr_Ct103;
  r_Ct(1,0)=-cr_Ct87*(cr_Ct106*cr_Ct60 + cr_Ct106*cr_Ct84 + cr_Ct88);
  r_Ct(1,1)=cr_Ct87*(-cr_Ct106*cr_Ct92 - cr_Ct106*cr_Ct94 + cr_Ct8);
  r_Ct(1,2)=-cr_Ct87*(cr_Ct106*cr_Ct95 + cr_Ct106*cr_Ct96 + cr_Ct88);
  r_Ct(1,3)=-cr_Ct107*cr_Ct99;
  r_Ct(1,4)=-cr_Ct102*cr_Ct107;
  r_Ct(1,5)=-cr_Ct103*cr_Ct107;
  r_Ct(2,0)=-cr_Ct87*(cr_Ct109*cr_Ct60 + cr_Ct109*cr_Ct84 + cr_Ct88);
  r_Ct(2,1)=-cr_Ct87*(cr_Ct109*cr_Ct92 + cr_Ct109*cr_Ct94 + cr_Ct88);
  r_Ct(2,2)=cr_Ct87*(-cr_Ct109*cr_Ct95 - cr_Ct109*cr_Ct96 + cr_Ct8);
  r_Ct(2,3)=-cr_Ct110*cr_Ct99;
  r_Ct(2,4)=-cr_Ct102*cr_Ct110;
  r_Ct(2,5)=-cr_Ct103*cr_Ct110;
  r_Ct(3,0)=-cr_Ct111*cr_Ct113;
  r_Ct(3,1)=-cr_Ct113*cr_Ct114;
  r_Ct(3,2)=-cr_Ct113*cr_Ct115;
  r_Ct(3,3)=cr_Ct86*(cr_Ct116 - cr_Ct119*cr_Ct50 - cr_Ct120*cr_Ct50);
  r_Ct(3,4)=cr_Ct122;
  r_Ct(3,5)=cr_Ct123;
  r_Ct(4,0)=-cr_Ct111*cr_Ct125;
  r_Ct(4,1)=-cr_Ct114*cr_Ct125;
  r_Ct(4,2)=-cr_Ct115*cr_Ct125;
  r_Ct(4,3)=cr_Ct122;
  r_Ct(4,4)=cr_Ct86*(cr_Ct116 - cr_Ct119*cr_Ct51 - cr_Ct120*cr_Ct51);
  r_Ct(4,5)=cr_Ct126;
  r_Ct(5,0)=-cr_Ct111*cr_Ct127;
  r_Ct(5,1)=-cr_Ct114*cr_Ct127;
  r_Ct(5,2)=-cr_Ct115*cr_Ct127;
  r_Ct(5,3)=cr_Ct123;
  r_Ct(5,4)=cr_Ct126;
  r_Ct(5,5)=cr_Ct86*(cr_Ct116 - cr_Ct119*cr_Ct52 - cr_Ct120*cr_Ct52);

}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<RankineYieldSurface<RankinePlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  const auto &r_props = rValues.GetMaterialProperties();
  const double Young = r_props[YOUNG_MODULUS];
  const double nu = r_props[POISSON_RATIO];
  const double Gf = r_props[FRACTURE_ENERGY];
  const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
  const double threshold = r_props[YIELD_STRESS];
  auto &r_Ct = rValues.GetConstitutiveMatrix();
  const auto &r_strain = rValues.GetStrainVector();

  const double cr_Ct0 = nu - 1.0;
  const double cr_Ct1 = 2.0*nu;
  const double cr_Ct2 = 1.0/(cr_Ct1 - 1.0);
  const double cr_Ct3 = nu*r_strain[1];
  const double cr_Ct4 = -cr_Ct3;
  const double cr_Ct5 = nu*r_strain[0];
  const double cr_Ct6 = -cr_Ct5;
  const double cr_Ct7 = cr_Ct0*r_strain[2];
  const double cr_Ct8 = cr_Ct6 + cr_Ct7;
  const double cr_Ct9 = cr_Ct4 + cr_Ct8;
  const double cr_Ct10 = cr_Ct2*cr_Ct9;
  const double cr_Ct11 = nu - 0.5;
  const double cr_Ct12 = std::pow(cr_Ct11, -2);
  const double cr_Ct13 = cr_Ct0*r_strain[0];
  const double cr_Ct14 = nu*r_strain[2];
  const double cr_Ct15 = -cr_Ct14;
  const double cr_Ct16 = cr_Ct13 + cr_Ct15;
  const double cr_Ct17 = cr_Ct16 + cr_Ct4;
  const double cr_Ct18 = cr_Ct0*r_strain[1];
  const double cr_Ct19 = cr_Ct15 + cr_Ct18 + cr_Ct6;
  const double cr_Ct20 = std::pow(r_strain[3], 2);
  const double cr_Ct21 = std::pow(r_strain[4], 2);
  const double cr_Ct22 = std::pow(r_strain[5], 2);
  const double cr_Ct23 = cr_Ct20 + cr_Ct21 + cr_Ct22;
  const double cr_Ct24 = cr_Ct12*std::pow(cr_Ct17, 2) + cr_Ct12*std::pow(cr_Ct19, 2) + cr_Ct12*std::pow(cr_Ct9, 2) + cr_Ct23;
  const double cr_Ct25 = nu + 1.0;
  const double cr_Ct26 = std::pow(Young, 2)/std::pow(cr_Ct25, 2);
  const double cr_Ct27 = cr_Ct24*cr_Ct26;
  const double cr_Ct28 = std::sqrt(cr_Ct27);
  const double cr_Ct29 = 1.0/cr_Ct28;
  const double cr_Ct30 = Young/cr_Ct25;
  const double cr_Ct31 = cr_Ct29*cr_Ct30;
  const double cr_Ct32 = 0.66666666666666663*cr_Ct31;
  const double cr_Ct33 = cr_Ct19*cr_Ct2;
  const double cr_Ct34 = cr_Ct17*cr_Ct2;
  const double cr_Ct35 = 1.0/cr_Ct24;
  const double cr_Ct36 = cr_Ct35;
  const double cr_Ct37 = cr_Ct20*cr_Ct36;
  const double cr_Ct38 = cr_Ct21;
  const double cr_Ct39 = cr_Ct35*cr_Ct38;
  const double cr_Ct40 = cr_Ct22*cr_Ct36;
  const double cr_Ct41 = cr_Ct12*cr_Ct9;
  const double cr_Ct42 = cr_Ct17*cr_Ct41;
  const double cr_Ct43 = 2*cr_Ct5;
  const double cr_Ct44 = 2*cr_Ct3;
  const double cr_Ct45 = -cr_Ct44;
  const double cr_Ct46 = 2*cr_Ct14;
  const double cr_Ct47 = cr_Ct13 + cr_Ct18 - cr_Ct43 + cr_Ct45 - cr_Ct46 + cr_Ct7;
  const double cr_Ct48 = cr_Ct12*cr_Ct35*std::pow(cr_Ct47, 2);
  const double cr_Ct49 = cr_Ct12*cr_Ct19*(cr_Ct16 + cr_Ct45 + cr_Ct8);
  const double cr_Ct50 = std::sqrt(-cr_Ct36*cr_Ct42 - cr_Ct36*cr_Ct49 + cr_Ct37 + cr_Ct39 + cr_Ct40 + 0.33333333333333331*cr_Ct48);
  const double cr_Ct51 = -nu;
  const double cr_Ct52 = cr_Ct51 + 0.5;
  const double cr_Ct53 = std::pow(cr_Ct52, -2);
  const double cr_Ct54 = cr_Ct51 + 1.0;
  const double cr_Ct55 = cr_Ct54*r_strain[2];
  const double cr_Ct56 = cr_Ct5 + cr_Ct55;
  const double cr_Ct57 = cr_Ct3 + cr_Ct56;
  const double cr_Ct58 = cr_Ct54*r_strain[1];
  const double cr_Ct59 = cr_Ct14 + cr_Ct5 + cr_Ct58;
  const double cr_Ct60 = cr_Ct54*r_strain[0];
  const double cr_Ct61 = cr_Ct14 + cr_Ct60;
  const double cr_Ct62 = cr_Ct3 + cr_Ct61;
  const double cr_Ct63 = cr_Ct23 + cr_Ct53*std::pow(cr_Ct57, 2) + cr_Ct53*std::pow(cr_Ct59, 2) + cr_Ct53*std::pow(cr_Ct62, 2);
  const double cr_Ct64 = 1.0/cr_Ct63;
  const double cr_Ct65 = cr_Ct64;
  const double cr_Ct66 = cr_Ct20*cr_Ct65;
  const double cr_Ct67 = cr_Ct38*cr_Ct64;
  const double cr_Ct68 = cr_Ct22*cr_Ct65;
  const double cr_Ct69 = cr_Ct53*cr_Ct65;
  const double cr_Ct70 = cr_Ct57*cr_Ct62;
  const double cr_Ct71 = cr_Ct69*cr_Ct70;
  const double cr_Ct72 = cr_Ct43 + cr_Ct44 + cr_Ct46 + cr_Ct55 + cr_Ct58 + cr_Ct60;
  const double cr_Ct73 = std::pow(cr_Ct72, 2);
  const double cr_Ct74 = cr_Ct53*cr_Ct64;
  const double cr_Ct75 = cr_Ct73*cr_Ct74;
  const double cr_Ct76 = 0.33333333333333331*cr_Ct75;
  const double cr_Ct77 = cr_Ct44 + cr_Ct56 + cr_Ct61;
  const double cr_Ct78 = cr_Ct59*cr_Ct77;
  const double cr_Ct79 = cr_Ct69*cr_Ct78;
  const double cr_Ct80 = cr_Ct66 + cr_Ct67 + cr_Ct68 - cr_Ct71 + cr_Ct76 - cr_Ct79;
  const double cr_Ct81 = std::pow(cr_Ct80, 3);
  const double cr_Ct82 = std::pow(cr_Ct81, -1.0/2.0);
  const double cr_Ct83 = r_strain[4]*r_strain[5];
  const double cr_Ct84 = cr_Ct83*r_strain[3];
  const double cr_Ct85 = cr_Ct53;
  const double cr_Ct86 = cr_Ct57*cr_Ct59;
  const double cr_Ct87 = cr_Ct85*cr_Ct86;
  const double cr_Ct88 = cr_Ct34*(cr_Ct38 - cr_Ct87);
  const double cr_Ct89 = 9.0*cr_Ct64;
  const double cr_Ct90 = cr_Ct20*cr_Ct89;
  const double cr_Ct91 = cr_Ct21*cr_Ct89;
  const double cr_Ct92 = cr_Ct22*cr_Ct89;
  const double cr_Ct93 = cr_Ct53*cr_Ct89;
  const double cr_Ct94 = cr_Ct70*cr_Ct93;
  const double cr_Ct95 = cr_Ct78*cr_Ct93;
  const double cr_Ct96 = cr_Ct2*cr_Ct47;
  const double cr_Ct97 = cr_Ct96*(cr_Ct90 + cr_Ct91 + cr_Ct92 - cr_Ct94 - cr_Ct95);
  const double cr_Ct98 = 0.037037037037037035*cr_Ct97 + 0.037037037037037035*cr_Ct35*std::pow(cr_Ct47, 3)/std::pow(cr_Ct11, 3);
  const double cr_Ct99 = -cr_Ct10*cr_Ct37 - cr_Ct33*cr_Ct40 + cr_Ct36*cr_Ct84 - cr_Ct65*cr_Ct88 + cr_Ct98;
  const double cr_Ct100 = 0.33333333333333331*std::acos(5.196152422706632*cr_Ct31*cr_Ct82*cr_Ct99);
  const double cr_Ct101 = std::cos(cr_Ct100);
  const double cr_Ct102 = cr_Ct101*cr_Ct50;
  const double cr_Ct103 = cr_Ct10*cr_Ct32 + 1.1547005383792515*cr_Ct102 + cr_Ct32*cr_Ct33 + cr_Ct32*cr_Ct34;
  const double cr_Ct104 = 1.0/cr_Ct103;
  const double cr_Ct105 = cr_Ct104*threshold;
  const double cr_Ct106 = 2.0*cr_Ct105;
  const double cr_Ct107 = cr_Ct0*cr_Ct106;
  const double cr_Ct108 = cr_Ct57*nu;
  const double cr_Ct109 = cr_Ct59*nu;
  const double cr_Ct110 = cr_Ct54*cr_Ct62;
  const double cr_Ct111 = cr_Ct108 + cr_Ct109 + cr_Ct110;
  const double cr_Ct112 = cr_Ct111*cr_Ct53;
  const double cr_Ct113 = cr_Ct106*cr_Ct35;
  const double cr_Ct114 = cr_Ct113*cr_Ct17;
  const double cr_Ct115 = cr_Ct62*nu;
  const double cr_Ct116 = 0.5*cr_Ct74;
  const double cr_Ct117 = -cr_Ct115*cr_Ct116;
  const double cr_Ct118 = cr_Ct54*cr_Ct57;
  const double cr_Ct119 = cr_Ct77*nu;
  const double cr_Ct120 = -cr_Ct116*cr_Ct119;
  const double cr_Ct121 = std::pow(cr_Ct63, -2);
  const double cr_Ct122 = cr_Ct121*cr_Ct85;
  const double cr_Ct123 = cr_Ct111*cr_Ct122;
  const double cr_Ct124 = cr_Ct121*cr_Ct38;
  const double cr_Ct125 = std::pow(cr_Ct52, -4);
  const double cr_Ct126 = cr_Ct111*cr_Ct125;
  const double cr_Ct127 = cr_Ct121*cr_Ct126;
  const double cr_Ct128 = cr_Ct70;
  const double cr_Ct129 = cr_Ct127*cr_Ct78;
  const double cr_Ct130 = cr_Ct111*cr_Ct74;
  const double cr_Ct131 = cr_Ct130*cr_Ct57;
  const double cr_Ct132 = cr_Ct130*cr_Ct62;
  const double cr_Ct133 = cr_Ct59*(-cr_Ct131 - cr_Ct132 + 1.0);
  const double cr_Ct134 = 2.0*cr_Ct64;
  const double cr_Ct135 = cr_Ct112*cr_Ct134;
  const double cr_Ct136 = cr_Ct1 + 2.0;
  const double cr_Ct137 = -cr_Ct135*cr_Ct57 - cr_Ct135*cr_Ct59 - cr_Ct135*cr_Ct62 + cr_Ct136;
  const double cr_Ct138 = cr_Ct137*cr_Ct72;
  const double cr_Ct139 = 0.16666666666666666*cr_Ct74;
  const double cr_Ct140 = std::sqrt(cr_Ct80);
  const double cr_Ct141 = 1.0/cr_Ct140;
  const double cr_Ct142 = 1.0/(1.0 - cr_Ct1);
  const double cr_Ct143 = -cr_Ct38 + cr_Ct87;
  const double cr_Ct144 = cr_Ct143*cr_Ct62;
  const double cr_Ct145 = cr_Ct142*cr_Ct144;
  const double cr_Ct146 = cr_Ct142*cr_Ct57;
  const double cr_Ct147 = cr_Ct142*cr_Ct59;
  const double cr_Ct148 = std::pow(cr_Ct72, 3);
  const double cr_Ct149 = cr_Ct142*(-cr_Ct90 - cr_Ct91 - cr_Ct92 + cr_Ct94 + cr_Ct95);
  const double cr_Ct150 = cr_Ct149*cr_Ct72;
  const double cr_Ct151 = 0.037037037037037035*cr_Ct148*cr_Ct64/std::pow(cr_Ct52, 3) - 0.037037037037037035*cr_Ct150;
  const double cr_Ct152 = cr_Ct145*cr_Ct65 - cr_Ct146*cr_Ct66 - cr_Ct147*cr_Ct68 + cr_Ct151 + cr_Ct65*cr_Ct84;
  const double cr_Ct153 = -cr_Ct66 - cr_Ct67 - cr_Ct68 + cr_Ct71 - cr_Ct76 + cr_Ct79;
  const double cr_Ct154 = std::pow(cr_Ct153, 3);
  const double cr_Ct155 = cr_Ct26*cr_Ct63;
  const double cr_Ct156 = std::sqrt(cr_Ct155);
  const double cr_Ct157 = 1.0/cr_Ct156;
  const double cr_Ct158 = cr_Ct157*cr_Ct30;
  const double cr_Ct159 = cr_Ct158/std::sqrt(-cr_Ct154);
  const double cr_Ct160 = 0.33333333333333331*std::acos(5.196152422706632*cr_Ct152*cr_Ct159);
  const double cr_Ct161 = std::cos(cr_Ct160);
  const double cr_Ct162 = 1.1547005383792515*cr_Ct161;
  const double cr_Ct163 = cr_Ct141*cr_Ct162;
  const double cr_Ct164 = std::pow(cr_Ct155, -3.0/2.0);
  const double cr_Ct165 = 0.66666666666666663*std::pow(Young, 3)/std::pow(cr_Ct25, 3);
  const double cr_Ct166 = cr_Ct164*cr_Ct165;
  const double cr_Ct167 = cr_Ct112*cr_Ct166;
  const double cr_Ct168 = cr_Ct142*cr_Ct62;
  const double cr_Ct169 = cr_Ct142*cr_Ct64;
  const double cr_Ct170 = cr_Ct169*cr_Ct20;
  const double cr_Ct171 = 5.196152422706632*nu;
  const double cr_Ct172 = -cr_Ct170*cr_Ct171;
  const double cr_Ct173 = cr_Ct169*cr_Ct22;
  const double cr_Ct174 = -cr_Ct171*cr_Ct173;
  const double cr_Ct175 = 5.196152422706632*cr_Ct169;
  const double cr_Ct176 = cr_Ct143*cr_Ct175;
  const double cr_Ct177 = 15.588457268119896*cr_Ct121;
  const double cr_Ct178 = cr_Ct177*cr_Ct84;
  const double cr_Ct179 = cr_Ct146*cr_Ct177;
  const double cr_Ct180 = cr_Ct179*cr_Ct20;
  const double cr_Ct181 = cr_Ct147*cr_Ct177;
  const double cr_Ct182 = cr_Ct181*cr_Ct22;
  const double cr_Ct183 = 5.196152422706632*cr_Ct121;
  const double cr_Ct184 = cr_Ct145*cr_Ct183;
  const double cr_Ct185 = cr_Ct108*cr_Ct85;
  const double cr_Ct186 = cr_Ct109*cr_Ct85;
  const double cr_Ct187 = cr_Ct134*cr_Ct21;
  const double cr_Ct188 = cr_Ct134*cr_Ct86;
  const double cr_Ct189 = cr_Ct175*cr_Ct62;
  const double cr_Ct190 = 0.57735026918962573*cr_Ct142*cr_Ct75;
  const double cr_Ct191 = 0.037037037037037035*cr_Ct59;
  const double cr_Ct192 = -0.037037037037037035*nu - 0.037037037037037035;
  const double cr_Ct193 = 5.196152422706632*cr_Ct149;
  const double cr_Ct194 = cr_Ct115*cr_Ct93;
  const double cr_Ct195 = cr_Ct119*cr_Ct93;
  const double cr_Ct196 = 18.0*cr_Ct121;
  const double cr_Ct197 = cr_Ct196*cr_Ct20;
  const double cr_Ct198 = cr_Ct196*cr_Ct21;
  const double cr_Ct199 = cr_Ct196*cr_Ct22;
  const double cr_Ct200 = cr_Ct126*cr_Ct70;
  const double cr_Ct201 = cr_Ct142*cr_Ct72;
  const double cr_Ct202 = 0.19245008972987526*cr_Ct201;
  const double cr_Ct203 = 3.0*cr_Ct74;
  const double cr_Ct204 = cr_Ct115*cr_Ct203;
  const double cr_Ct205 = cr_Ct119*cr_Ct203;
  const double cr_Ct206 = 6.0*cr_Ct121;
  const double cr_Ct207 = cr_Ct20*cr_Ct206;
  const double cr_Ct208 = cr_Ct206*cr_Ct21;
  const double cr_Ct209 = cr_Ct206*cr_Ct22;
  const double cr_Ct210 = 2.598076211353316*cr_Ct152/cr_Ct153;
  const double cr_Ct211 = cr_Ct64*r_strain[3];
  const double cr_Ct212 = 0.074074074074074056*cr_Ct140*cr_Ct159*std::sin(cr_Ct160)/std::sqrt(0.037037037037037035 + cr_Ct64*std::pow(cr_Ct144*cr_Ct169 + cr_Ct151 - cr_Ct170*cr_Ct57 - cr_Ct173*cr_Ct59 + cr_Ct211*cr_Ct83, 2)/cr_Ct154);
  const double cr_Ct213 = cr_Ct142*cr_Ct158;
  const double cr_Ct214 = 0.66666666666666663*cr_Ct213;
  const double cr_Ct215 = 1.3333333333333333*cr_Ct213*nu + cr_Ct214*cr_Ct54;
  const double cr_Ct216 = -cr_Ct146*cr_Ct167 - cr_Ct147*cr_Ct167 + cr_Ct163*(-cr_Ct112*cr_Ct124 - cr_Ct116*cr_Ct118 - cr_Ct116*cr_Ct133 + cr_Ct117 + cr_Ct120 - cr_Ct123*cr_Ct20 - cr_Ct123*cr_Ct22 + cr_Ct127*cr_Ct128 + 0.5*cr_Ct129 + cr_Ct138*cr_Ct139) - cr_Ct167*cr_Ct168 + cr_Ct212*(-cr_Ct112*cr_Ct178 + cr_Ct112*cr_Ct180 + cr_Ct112*cr_Ct182 - cr_Ct112*cr_Ct184 + cr_Ct137*cr_Ct190 + cr_Ct172 + cr_Ct174 + cr_Ct176*cr_Ct54 + cr_Ct189*(cr_Ct112*cr_Ct187 - cr_Ct126*cr_Ct188 + cr_Ct185 + cr_Ct186) + cr_Ct193*(cr_Ct130*cr_Ct191 + 0.037037037037037035*cr_Ct131 + 0.037037037037037035*cr_Ct132 + cr_Ct192) - cr_Ct202*(cr_Ct112*cr_Ct197 + cr_Ct112*cr_Ct198 + cr_Ct112*cr_Ct199 + cr_Ct118*cr_Ct93 - 9.0*cr_Ct129 + cr_Ct133*cr_Ct93 + cr_Ct194 + cr_Ct195 - cr_Ct196*cr_Ct200) - cr_Ct210*(cr_Ct112*cr_Ct207 + cr_Ct112*cr_Ct208 + cr_Ct112*cr_Ct209 + cr_Ct118*cr_Ct203 - 3.0*cr_Ct129 + cr_Ct133*cr_Ct203 - cr_Ct138*cr_Ct69 - cr_Ct200*cr_Ct206 + cr_Ct204 + cr_Ct205)) + cr_Ct215;
  const double cr_Ct217 = 0.57735026918962573*cr_Ct31;
  const double cr_Ct218 = threshold/pow(cr_Ct10*cr_Ct217 + cr_Ct102 + cr_Ct217*cr_Ct33 + cr_Ct217*cr_Ct34, 2);
  const double cr_Ct219 = 1.5000000000000002*cr_Ct218;
  const double cr_Ct220 = cr_Ct17*cr_Ct219;
  const double cr_Ct221 = 0.66666666666666663*cr_Ct158;
  const double cr_Ct222 = cr_Ct26*(cr_Ct140*cr_Ct162 + cr_Ct146*cr_Ct221 + cr_Ct147*cr_Ct221 + cr_Ct214*cr_Ct62);
  const double cr_Ct223 = cr_Ct157*cr_Ct222;
  const double cr_Ct224 = cr_Ct112*cr_Ct223;
  const double cr_Ct225 = cr_Ct156*cr_Ct216;
  const double cr_Ct226 = cr_Ct224 + cr_Ct225;
  const double cr_Ct227 = 1.0/(Gf*Young/(characteristic_length*std::pow(threshold, 2)) - 0.5);
  const double cr_Ct228 = cr_Ct227;
  const double cr_Ct229 = cr_Ct104*cr_Ct228;
  const double cr_Ct230 = cr_Ct17*cr_Ct229;
  const double cr_Ct231 = cr_Ct31*std::exp(cr_Ct228*(-0.5*cr_Ct103*cr_Ct28/threshold + 1.0));
  const double cr_Ct232 = cr_Ct2*cr_Ct231;
  const double cr_Ct233 = cr_Ct1*cr_Ct105;
  const double cr_Ct234 = cr_Ct54*cr_Ct59;
  const double cr_Ct235 = cr_Ct108 + cr_Ct115 + cr_Ct234;
  const double cr_Ct236 = cr_Ct235*cr_Ct53;
  const double cr_Ct237 = -cr_Ct108*cr_Ct116;
  const double cr_Ct238 = cr_Ct54*cr_Ct77;
  const double cr_Ct239 = cr_Ct122*cr_Ct235;
  const double cr_Ct240 = cr_Ct125*cr_Ct235;
  const double cr_Ct241 = cr_Ct121*cr_Ct240;
  const double cr_Ct242 = cr_Ct241*cr_Ct78;
  const double cr_Ct243 = cr_Ct134*cr_Ct236;
  const double cr_Ct244 = -cr_Ct243*cr_Ct57 - cr_Ct243*cr_Ct62;
  const double cr_Ct245 = cr_Ct59*cr_Ct74*(cr_Ct244 + 4.0*nu);
  const double cr_Ct246 = cr_Ct136 - cr_Ct243*cr_Ct59 + cr_Ct244;
  const double cr_Ct247 = cr_Ct246*cr_Ct72;
  const double cr_Ct248 = cr_Ct166*cr_Ct236;
  const double cr_Ct249 = 5.196152422706632*cr_Ct54;
  const double cr_Ct250 = cr_Ct176*nu;
  const double cr_Ct251 = cr_Ct235*cr_Ct74;
  const double cr_Ct252 = 0.037037037037037035*cr_Ct251;
  const double cr_Ct253 = cr_Ct108*cr_Ct93;
  const double cr_Ct254 = cr_Ct240*cr_Ct70;
  const double cr_Ct255 = cr_Ct108*cr_Ct203;
  const double cr_Ct256 = -cr_Ct146*cr_Ct248 - cr_Ct147*cr_Ct248 + cr_Ct163*(-cr_Ct116*cr_Ct238 + cr_Ct117 - cr_Ct124*cr_Ct236 + cr_Ct128*cr_Ct241 + cr_Ct139*cr_Ct247 - cr_Ct20*cr_Ct239 - cr_Ct22*cr_Ct239 + cr_Ct237 + 0.5*cr_Ct242 - 0.25*cr_Ct245) - cr_Ct168*cr_Ct248 + cr_Ct212*(cr_Ct172 - cr_Ct173*cr_Ct249 - cr_Ct178*cr_Ct236 + cr_Ct180*cr_Ct236 + cr_Ct182*cr_Ct236 - cr_Ct184*cr_Ct236 + cr_Ct189*(cr_Ct118*cr_Ct85 + cr_Ct186 + cr_Ct187*cr_Ct236 - cr_Ct188*cr_Ct240) + cr_Ct190*cr_Ct246 + cr_Ct193*(cr_Ct191*cr_Ct251 + cr_Ct192 + cr_Ct252*cr_Ct57 + cr_Ct252*cr_Ct62) - cr_Ct202*(cr_Ct194 - cr_Ct196*cr_Ct254 + cr_Ct197*cr_Ct236 + cr_Ct198*cr_Ct236 + cr_Ct199*cr_Ct236 + cr_Ct238*cr_Ct93 - 9.0*cr_Ct242 + 4.5*cr_Ct245 + cr_Ct253) - cr_Ct210*(cr_Ct203*cr_Ct238 + cr_Ct204 - cr_Ct206*cr_Ct254 + cr_Ct207*cr_Ct236 + cr_Ct208*cr_Ct236 + cr_Ct209*cr_Ct236 - 3.0*cr_Ct242 + 1.5*cr_Ct245 - cr_Ct247*cr_Ct69 + cr_Ct255) + cr_Ct250) + cr_Ct215;
  const double cr_Ct257 = cr_Ct223*cr_Ct236;
  const double cr_Ct258 = cr_Ct156*cr_Ct256;
  const double cr_Ct259 = cr_Ct257 + cr_Ct258;
  const double cr_Ct260 = cr_Ct109 + cr_Ct115 + cr_Ct118;
  const double cr_Ct261 = cr_Ct260*cr_Ct53;
  const double cr_Ct262 = cr_Ct122*cr_Ct260;
  const double cr_Ct263 = cr_Ct125*cr_Ct260;
  const double cr_Ct264 = cr_Ct121*cr_Ct263;
  const double cr_Ct265 = cr_Ct264*cr_Ct78;
  const double cr_Ct266 = cr_Ct260*cr_Ct74;
  const double cr_Ct267 = cr_Ct266*cr_Ct57;
  const double cr_Ct268 = cr_Ct266*cr_Ct62;
  const double cr_Ct269 = cr_Ct59*(-cr_Ct267 - cr_Ct268 + 1.0);
  const double cr_Ct270 = cr_Ct134*cr_Ct261;
  const double cr_Ct271 = cr_Ct136 - cr_Ct270*cr_Ct57 - cr_Ct270*cr_Ct59 - cr_Ct270*cr_Ct62;
  const double cr_Ct272 = cr_Ct271*cr_Ct72;
  const double cr_Ct273 = cr_Ct166*cr_Ct261;
  const double cr_Ct274 = cr_Ct263*cr_Ct70;
  const double cr_Ct275 = -cr_Ct146*cr_Ct273 - cr_Ct147*cr_Ct273 + cr_Ct163*(-cr_Ct110*cr_Ct116 - cr_Ct116*cr_Ct269 + cr_Ct120 - cr_Ct124*cr_Ct261 + cr_Ct128*cr_Ct264 + cr_Ct139*cr_Ct272 - cr_Ct20*cr_Ct262 - cr_Ct22*cr_Ct262 + cr_Ct237 + 0.5*cr_Ct265) - cr_Ct168*cr_Ct273 + cr_Ct212*(-cr_Ct170*cr_Ct249 + cr_Ct174 - cr_Ct178*cr_Ct261 + cr_Ct180*cr_Ct261 + cr_Ct182*cr_Ct261 - cr_Ct184*cr_Ct261 + cr_Ct189*(cr_Ct185 + cr_Ct187*cr_Ct261 - cr_Ct188*cr_Ct263 + cr_Ct234*cr_Ct85) + cr_Ct190*cr_Ct271 + cr_Ct193*(cr_Ct191*cr_Ct266 + cr_Ct192 + 0.037037037037037035*cr_Ct267 + 0.037037037037037035*cr_Ct268) - cr_Ct202*(cr_Ct110*cr_Ct93 + cr_Ct195 - cr_Ct196*cr_Ct274 + cr_Ct197*cr_Ct261 + cr_Ct198*cr_Ct261 + cr_Ct199*cr_Ct261 + cr_Ct253 - 9.0*cr_Ct265 + cr_Ct269*cr_Ct93) - cr_Ct210*(cr_Ct110*cr_Ct203 + cr_Ct203*cr_Ct269 + cr_Ct205 - cr_Ct206*cr_Ct274 + cr_Ct207*cr_Ct261 + cr_Ct208*cr_Ct261 + cr_Ct209*cr_Ct261 + cr_Ct255 - 3.0*cr_Ct265 - cr_Ct272*cr_Ct69) + cr_Ct250) + cr_Ct215;
  const double cr_Ct276 = cr_Ct223*cr_Ct261;
  const double cr_Ct277 = cr_Ct156*cr_Ct275;
  const double cr_Ct278 = cr_Ct276 + cr_Ct277;
  const double cr_Ct279 = std::pow(cr_Ct27, -3.0/2.0);
  const double cr_Ct280 = cr_Ct165*r_strain[3];
  const double cr_Ct281 = cr_Ct279*cr_Ct280;
  const double cr_Ct282 = 0.57735026918962573*r_strain[3];
  const double cr_Ct283 = 2.0*cr_Ct121;
  const double cr_Ct284 = cr_Ct20*cr_Ct283;
  const double cr_Ct285 = 2.0*cr_Ct21;
  const double cr_Ct286 = cr_Ct121*cr_Ct285;
  const double cr_Ct287 = cr_Ct22*cr_Ct283;
  const double cr_Ct288 = cr_Ct283*cr_Ct53;
  const double cr_Ct289 = cr_Ct288*cr_Ct70;
  const double cr_Ct290 = cr_Ct53*cr_Ct73;
  const double cr_Ct291 = 0.66666666666666663*cr_Ct121*cr_Ct290;
  const double cr_Ct292 = cr_Ct288*cr_Ct78;
  const double cr_Ct293 = cr_Ct101*(-cr_Ct134 + cr_Ct284 + cr_Ct286 + cr_Ct287 - cr_Ct289 + cr_Ct291 - cr_Ct292)/cr_Ct50;
  const double cr_Ct294 = 5.196152422706632*cr_Ct83;
  const double cr_Ct295 = cr_Ct26*cr_Ct279;
  const double cr_Ct296 = cr_Ct20*cr_Ct83;
  const double cr_Ct297 = std::pow(Young, 4)/std::pow(cr_Ct25, 4);
  const double cr_Ct298 = 15.588457268119896*cr_Ct297/std::pow(cr_Ct27, 5.0/2.0);
  const double cr_Ct299 = cr_Ct295*r_strain[3];
  const double cr_Ct300 = std::pow(r_strain[3], 3);
  const double cr_Ct301 = cr_Ct298*r_strain[3];
  const double cr_Ct302 = cr_Ct19*cr_Ct22;
  const double cr_Ct303 = cr_Ct2*cr_Ct302;
  const double cr_Ct304 = cr_Ct29*r_strain[3];
  const double cr_Ct305 = std::pow(cr_Ct24, -2);
  const double cr_Ct306 = cr_Ct305*cr_Ct34*(-2.0*cr_Ct19*cr_Ct41 + cr_Ct285);
  const double cr_Ct307 = 5.196152422706632*cr_Ct211;
  const double cr_Ct308 = 1.1547005383792515*r_strain[3];
  const double cr_Ct309 = cr_Ct164*cr_Ct201;
  const double cr_Ct310 = cr_Ct26*cr_Ct309*cr_Ct48;
  const double cr_Ct311 = 0.19245008972987526*cr_Ct97;
  const double cr_Ct312 = cr_Ct157*r_strain[3];
  const double cr_Ct313 = 18.0*cr_Ct64;
  const double cr_Ct314 = cr_Ct196*cr_Ct53;
  const double cr_Ct315 = cr_Ct314*cr_Ct70;
  const double cr_Ct316 = cr_Ct314*cr_Ct78;
  const double cr_Ct317 = cr_Ct202*(-cr_Ct197 - cr_Ct198 - cr_Ct199 + cr_Ct313 + cr_Ct315 + cr_Ct316);
  const double cr_Ct318 = 6.0*cr_Ct305;
  const double cr_Ct319 = cr_Ct29*cr_Ct297;
  const double cr_Ct320 = 2.598076211353316*cr_Ct99*(-12.0*cr_Ct142*cr_Ct164*cr_Ct319*cr_Ct33*cr_Ct77 + cr_Ct20*cr_Ct318 + cr_Ct21*cr_Ct318 + cr_Ct22*cr_Ct318 - 3.0*cr_Ct305*cr_Ct49 + 8.0*cr_Ct309*cr_Ct319*cr_Ct96 - cr_Ct318*cr_Ct42 - 6.0*cr_Ct35)/cr_Ct80;
  const double cr_Ct321 = cr_Ct2*cr_Ct35;
  const double cr_Ct322 = cr_Ct20*cr_Ct9;
  const double cr_Ct323 = cr_Ct64*cr_Ct88;
  const double cr_Ct324 = 0.074074074074074056*cr_Ct30*cr_Ct50*cr_Ct82*std::sin(cr_Ct100)/std::sqrt(-cr_Ct35*std::pow(-cr_Ct302*cr_Ct321 - cr_Ct321*cr_Ct322 - cr_Ct323 + cr_Ct35*cr_Ct84 + cr_Ct98, 2)/cr_Ct81 + 0.037037037037037035);
  const double cr_Ct325 = cr_Ct10*cr_Ct281 + cr_Ct281*cr_Ct33 + cr_Ct281*cr_Ct34 + cr_Ct282*cr_Ct293 - cr_Ct324*(cr_Ct10*cr_Ct298*cr_Ct300 - 10.392304845413264*cr_Ct10*cr_Ct299 + cr_Ct294*cr_Ct295 + cr_Ct295*cr_Ct307*cr_Ct88 - cr_Ct296*cr_Ct298 - cr_Ct299*cr_Ct311 + cr_Ct301*cr_Ct303 + 5.196152422706632*cr_Ct304*cr_Ct306 + cr_Ct304*cr_Ct320 - cr_Ct308*cr_Ct310 + cr_Ct312*cr_Ct317);
  const double cr_Ct326 = cr_Ct222*cr_Ct312;
  const double cr_Ct327 = cr_Ct164*cr_Ct280;
  const double cr_Ct328 = cr_Ct141*cr_Ct161*(cr_Ct134 - cr_Ct284 - cr_Ct286 - cr_Ct287 + cr_Ct289 - cr_Ct291 + cr_Ct292);
  const double cr_Ct329 = cr_Ct177*r_strain[3];
  const double cr_Ct330 = cr_Ct22*cr_Ct329;
  const double cr_Ct331 = cr_Ct121*cr_Ct142*cr_Ct148*cr_Ct53;
  const double cr_Ct332 = cr_Ct183*r_strain[3];
  const double cr_Ct333 = cr_Ct53*cr_Ct86;
  const double cr_Ct334 = cr_Ct168*(cr_Ct285 - 2.0*cr_Ct333);
  const double cr_Ct335 = 0.19245008972987526*cr_Ct150;
  const double cr_Ct336 = cr_Ct202*(cr_Ct197 + cr_Ct198 + cr_Ct199 - cr_Ct313 - cr_Ct315 - cr_Ct316);
  const double cr_Ct337 = cr_Ct206*cr_Ct53;
  const double cr_Ct338 = cr_Ct210*(cr_Ct207 + cr_Ct208 + cr_Ct209 + cr_Ct283*cr_Ct290 - cr_Ct337*cr_Ct70 - cr_Ct337*cr_Ct78 - 6.0*cr_Ct64);
  const double cr_Ct339 = cr_Ct156*(-cr_Ct146*cr_Ct327 - cr_Ct147*cr_Ct327 - cr_Ct168*cr_Ct327 + cr_Ct212*(-cr_Ct145*cr_Ct332 - 10.392304845413264*cr_Ct146*cr_Ct211 + cr_Ct147*cr_Ct330 - cr_Ct177*cr_Ct296 + cr_Ct179*cr_Ct300 + cr_Ct211*cr_Ct335 + cr_Ct294*cr_Ct64 - cr_Ct308*cr_Ct331 + cr_Ct332*cr_Ct334 - cr_Ct336*r_strain[3] - cr_Ct338*r_strain[3]) + cr_Ct282*cr_Ct328);
  const double cr_Ct340 = cr_Ct104*cr_Ct227;
  const double cr_Ct341 = -cr_Ct113*r_strain[3] + cr_Ct219*cr_Ct325 - cr_Ct340*(cr_Ct326 + cr_Ct339);
  const double cr_Ct342 = cr_Ct231*cr_Ct34;
  const double cr_Ct343 = cr_Ct165*r_strain[4];
  const double cr_Ct344 = cr_Ct279*cr_Ct343;
  const double cr_Ct345 = 0.57735026918962573*r_strain[4];
  const double cr_Ct346 = 5.196152422706632*r_strain[5];
  const double cr_Ct347 = cr_Ct21*r_strain[5];
  const double cr_Ct348 = cr_Ct298*r_strain[4];
  const double cr_Ct349 = cr_Ct2*cr_Ct322;
  const double cr_Ct350 = 5.196152422706632*r_strain[4];
  const double cr_Ct351 = cr_Ct295*cr_Ct323;
  const double cr_Ct352 = 1.1547005383792515*cr_Ct310;
  const double cr_Ct353 = cr_Ct134*cr_Ct333;
  const double cr_Ct354 = cr_Ct295*cr_Ct311;
  const double cr_Ct355 = cr_Ct157*cr_Ct317;
  const double cr_Ct356 = cr_Ct29*cr_Ct320;
  const double cr_Ct357 = cr_Ct10*cr_Ct344 + cr_Ct293*cr_Ct345 - cr_Ct324*(-cr_Ct29*cr_Ct34*cr_Ct350*cr_Ct64*(-cr_Ct187 + cr_Ct353 + 2.0) + cr_Ct299*cr_Ct346 - cr_Ct301*cr_Ct347 + cr_Ct303*cr_Ct348 + cr_Ct348*cr_Ct349 + cr_Ct350*cr_Ct351 - cr_Ct352*r_strain[4] - cr_Ct354*r_strain[4] + cr_Ct355*r_strain[4] + cr_Ct356*r_strain[4]) + cr_Ct33*cr_Ct344 + cr_Ct34*cr_Ct344;
  const double cr_Ct358 = cr_Ct223*r_strain[4];
  const double cr_Ct359 = cr_Ct164*cr_Ct343;
  const double cr_Ct360 = 1.1547005383792515*cr_Ct331;
  const double cr_Ct361 = cr_Ct335*cr_Ct64;
  const double cr_Ct362 = cr_Ct156*(-cr_Ct146*cr_Ct359 - cr_Ct147*cr_Ct359 - cr_Ct168*cr_Ct359 + cr_Ct212*(cr_Ct180*r_strain[4] + cr_Ct182*r_strain[4] - cr_Ct184*r_strain[4] + cr_Ct189*r_strain[4]*(cr_Ct187 - cr_Ct353 - 2.0) + cr_Ct307*r_strain[5] - cr_Ct329*cr_Ct347 - cr_Ct336*r_strain[4] - cr_Ct338*r_strain[4] - cr_Ct360*r_strain[4] + cr_Ct361*r_strain[4]) + cr_Ct328*cr_Ct345);
  const double cr_Ct363 = -cr_Ct113*r_strain[4] + cr_Ct219*cr_Ct357 - cr_Ct340*(cr_Ct358 + cr_Ct362);
  const double cr_Ct364 = cr_Ct165*r_strain[5];
  const double cr_Ct365 = cr_Ct279*cr_Ct364;
  const double cr_Ct366 = 0.57735026918962573*r_strain[5];
  const double cr_Ct367 = 10.392304845413264*r_strain[5];
  const double cr_Ct368 = std::pow(r_strain[5], 3);
  const double cr_Ct369 = cr_Ct10*cr_Ct365 + cr_Ct293*cr_Ct366 - cr_Ct324*(-cr_Ct22*cr_Ct301*r_strain[4] + cr_Ct29*cr_Ct306*cr_Ct346 - cr_Ct295*cr_Ct33*cr_Ct367 + cr_Ct298*cr_Ct33*cr_Ct368 + cr_Ct298*cr_Ct349*r_strain[5] + cr_Ct299*cr_Ct350 + cr_Ct346*cr_Ct351 - cr_Ct352*r_strain[5] - cr_Ct354*r_strain[5] + cr_Ct355*r_strain[5] + cr_Ct356*r_strain[5]) + cr_Ct33*cr_Ct365 + cr_Ct34*cr_Ct365;
  const double cr_Ct370 = cr_Ct223*r_strain[5];
  const double cr_Ct371 = cr_Ct164*cr_Ct364;
  const double cr_Ct372 = cr_Ct156*(-cr_Ct146*cr_Ct371 - cr_Ct147*cr_Ct371 - cr_Ct168*cr_Ct371 + cr_Ct212*(-cr_Ct169*cr_Ct367*cr_Ct59 + cr_Ct180*r_strain[5] + cr_Ct181*cr_Ct368 + cr_Ct183*cr_Ct334*r_strain[5] - cr_Ct184*r_strain[5] + cr_Ct307*r_strain[4] - cr_Ct330*r_strain[4] - cr_Ct336*r_strain[5] - cr_Ct338*r_strain[5] - cr_Ct360*r_strain[5] + cr_Ct361*r_strain[5]) + cr_Ct328*cr_Ct366);
  const double cr_Ct373 = -cr_Ct113*r_strain[5] + cr_Ct219*cr_Ct369 - cr_Ct340*(cr_Ct370 + cr_Ct372);
  const double cr_Ct374 = cr_Ct113*cr_Ct19;
  const double cr_Ct375 = cr_Ct19*cr_Ct219;
  const double cr_Ct376 = cr_Ct19*cr_Ct229;
  const double cr_Ct377 = cr_Ct231*cr_Ct33;
  const double cr_Ct378 = cr_Ct113*cr_Ct9;
  const double cr_Ct379 = cr_Ct219*cr_Ct9;
  const double cr_Ct380 = cr_Ct229*cr_Ct9;
  const double cr_Ct381 = cr_Ct10*cr_Ct231;
  const double cr_Ct382 = cr_Ct105*cr_Ct36;
  const double cr_Ct383 = 0.75000000000000011*cr_Ct218;
  const double cr_Ct384 = cr_Ct112*cr_Ct382 + cr_Ct216*cr_Ct383 + cr_Ct340*(0.5*cr_Ct224 + 0.5*cr_Ct225);
  const double cr_Ct385 = cr_Ct231*r_strain[3];
  const double cr_Ct386 = cr_Ct236*cr_Ct382 + cr_Ct256*cr_Ct383 + cr_Ct340*(0.5*cr_Ct257 + 0.5*cr_Ct258);
  const double cr_Ct387 = cr_Ct261*cr_Ct382 + cr_Ct275*cr_Ct383 + cr_Ct340*(0.5*cr_Ct276 + 0.5*cr_Ct277);
  const double cr_Ct388 = cr_Ct105;
  const double cr_Ct389 = cr_Ct325*cr_Ct383;
  const double cr_Ct390 = 0.5*cr_Ct340;
  const double cr_Ct391 = cr_Ct357*cr_Ct383;
  const double cr_Ct392 = -cr_Ct340*(0.5*cr_Ct358 + 0.5*cr_Ct362) - cr_Ct382*r_strain[4] + cr_Ct391;
  const double cr_Ct393 = cr_Ct369*cr_Ct383;
  const double cr_Ct394 = -cr_Ct340*(0.5*cr_Ct370 + 0.5*cr_Ct372) - cr_Ct382*r_strain[5] + cr_Ct393;
  const double cr_Ct395 = cr_Ct231*r_strain[4];
  const double cr_Ct396 = -cr_Ct340*(0.5*cr_Ct326 + 0.5*cr_Ct339) - cr_Ct382*r_strain[3] + cr_Ct389;
  const double cr_Ct397 = cr_Ct231*r_strain[5];
  r_Ct(0,0)=cr_Ct232*(cr_Ct107 - cr_Ct112*cr_Ct114 - cr_Ct216*cr_Ct220 - cr_Ct226*cr_Ct230);
  r_Ct(0,1)=-cr_Ct232*(cr_Ct114*cr_Ct236 + cr_Ct220*cr_Ct256 + cr_Ct230*cr_Ct259 + cr_Ct233);
  r_Ct(0,2)=-cr_Ct232*(cr_Ct114*cr_Ct261 + cr_Ct220*cr_Ct275 + cr_Ct230*cr_Ct278 + cr_Ct233);
  r_Ct(0,3)=cr_Ct341*cr_Ct342;
  r_Ct(0,4)=cr_Ct342*cr_Ct363;
  r_Ct(0,5)=cr_Ct342*cr_Ct373;
  r_Ct(1,0)=-cr_Ct232*(cr_Ct112*cr_Ct374 + cr_Ct216*cr_Ct375 + cr_Ct226*cr_Ct376 + cr_Ct233);
  r_Ct(1,1)=cr_Ct232*(cr_Ct107 - cr_Ct236*cr_Ct374 - cr_Ct256*cr_Ct375 - cr_Ct259*cr_Ct376);
  r_Ct(1,2)=-cr_Ct232*(cr_Ct233 + cr_Ct261*cr_Ct374 + cr_Ct275*cr_Ct375 + cr_Ct278*cr_Ct376);
  r_Ct(1,3)=cr_Ct341*cr_Ct377;
  r_Ct(1,4)=cr_Ct363*cr_Ct377;
  r_Ct(1,5)=cr_Ct373*cr_Ct377;
  r_Ct(2,0)=-cr_Ct232*(cr_Ct112*cr_Ct378 + cr_Ct216*cr_Ct379 + cr_Ct226*cr_Ct380 + cr_Ct233);
  r_Ct(2,1)=-cr_Ct232*(cr_Ct233 + cr_Ct236*cr_Ct378 + cr_Ct256*cr_Ct379 + cr_Ct259*cr_Ct380);
  r_Ct(2,2)=cr_Ct232*(cr_Ct107 - cr_Ct261*cr_Ct378 - cr_Ct275*cr_Ct379 - cr_Ct278*cr_Ct380);
  r_Ct(2,3)=cr_Ct341*cr_Ct381;
  r_Ct(2,4)=cr_Ct363*cr_Ct381;
  r_Ct(2,5)=cr_Ct373*cr_Ct381;
  r_Ct(3,0)=-cr_Ct384*cr_Ct385;
  r_Ct(3,1)=-cr_Ct385*cr_Ct386;
  r_Ct(3,2)=-cr_Ct385*cr_Ct387;
  r_Ct(3,3)=cr_Ct231*(-cr_Ct105*cr_Ct37 + cr_Ct388 + cr_Ct389*r_strain[3] - cr_Ct390*r_strain[3]*(cr_Ct326 + cr_Ct339));
  r_Ct(3,4)=cr_Ct385*cr_Ct392;
  r_Ct(3,5)=cr_Ct385*cr_Ct394;
  r_Ct(4,0)=-cr_Ct384*cr_Ct395;
  r_Ct(4,1)=-cr_Ct386*cr_Ct395;
  r_Ct(4,2)=-cr_Ct387*cr_Ct395;
  r_Ct(4,3)=cr_Ct395*cr_Ct396;
  r_Ct(4,4)=cr_Ct231*(-cr_Ct105*cr_Ct39 + cr_Ct388 - cr_Ct390*r_strain[4]*(cr_Ct358 + cr_Ct362) + cr_Ct391*r_strain[4]);
  r_Ct(4,5)=cr_Ct394*cr_Ct395;
  r_Ct(5,0)=-cr_Ct384*cr_Ct397;
  r_Ct(5,1)=-cr_Ct386*cr_Ct397;
  r_Ct(5,2)=-cr_Ct387*cr_Ct397;
  r_Ct(5,3)=cr_Ct396*cr_Ct397;
  r_Ct(5,4)=cr_Ct392*cr_Ct397;
  r_Ct(5,5)=cr_Ct231*(-cr_Ct105*cr_Ct40 + cr_Ct388 - cr_Ct390*r_strain[5]*(cr_Ct370 + cr_Ct372) + cr_Ct393*r_strain[5]);


}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<SimoJuYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

template<>
void AutomaticDifferentiationTangentUtilities<VonMisesYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{

  const auto &r_props = rValues.GetMaterialProperties();
  const double Young = r_props[YOUNG_MODULUS];
  const double nu = r_props[POISSON_RATIO];
  const double Gf = r_props[FRACTURE_ENERGY];
  const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
  double threshold;
  YieldSurfaceType::GetInitialUniaxialThreshold(rValues, threshold);
  auto &r_Ct = rValues.GetConstitutiveMatrix();
  const auto &r_strain = rValues.GetStrainVector();

  const double cr_Ct0 = nu - 1.0;
  const double cr_Ct1 = 1.0/(1.0 - 0.5*characteristic_length*std::pow(threshold, 2)/(Gf*Young));
  const double cr_Ct2 = nu*r_strain[0];
  const double cr_Ct3 = cr_Ct2;
  const double cr_Ct4 = nu*r_strain[1];
  const double cr_Ct5 = 0.5*cr_Ct4;
  const double cr_Ct6 = -cr_Ct5;
  const double cr_Ct7 = nu*r_strain[2];
  const double cr_Ct8 = 0.5*cr_Ct7;
  const double cr_Ct9 = -cr_Ct8;
  const double cr_Ct10 = -nu;
  const double cr_Ct11 = cr_Ct10 + 1.0;
  const double cr_Ct12 = cr_Ct11*r_strain[0];
  const double cr_Ct13 = cr_Ct11*r_strain[1];
  const double cr_Ct14 = 0.5*cr_Ct13;
  const double cr_Ct15 = cr_Ct11*r_strain[2];
  const double cr_Ct16 = 0.5*cr_Ct15;
  const double cr_Ct17 = -cr_Ct12 + cr_Ct14 + cr_Ct16 + cr_Ct3 + cr_Ct6 + cr_Ct9;
  const double cr_Ct18 = std::pow(cr_Ct10 + 0.5, -2);
  const double cr_Ct19 = 0.22222222222222224*cr_Ct18;
  const double cr_Ct20 = cr_Ct0*r_strain[1];
  const double cr_Ct21 = cr_Ct4;
  const double cr_Ct22 = cr_Ct0*r_strain[2];
  const double cr_Ct23 = 0.5*cr_Ct22;
  const double cr_Ct24 = 0.5*cr_Ct2;
  const double cr_Ct25 = cr_Ct0*r_strain[0];
  const double cr_Ct26 = 0.5*cr_Ct25;
  const double cr_Ct27 = -cr_Ct24 - cr_Ct26;
  const double cr_Ct28 = cr_Ct20 + cr_Ct21 - cr_Ct23 + cr_Ct27 + cr_Ct9;
  const double cr_Ct29 = std::pow(nu - 0.5, -2);
  const double cr_Ct30 = 0.22222222222222224*cr_Ct29;
  const double cr_Ct31 = cr_Ct7;
  const double cr_Ct32 = 0.5*cr_Ct20;
  const double cr_Ct33 = cr_Ct22 + cr_Ct27 + cr_Ct31 - cr_Ct32 + cr_Ct6;
  const double cr_Ct34 = std::pow(r_strain[3], 2);
  const double cr_Ct35 = std::pow(r_strain[4], 2);
  const double cr_Ct36 = std::pow(r_strain[5], 2);
  const double cr_Ct37 = cr_Ct34 + cr_Ct35 + cr_Ct36;
  const double cr_Ct38 = nu + 1.0;
  const double cr_Ct39 = std::pow(Young, 2)/std::pow(cr_Ct38, 2);
  const double cr_Ct40 = cr_Ct39*(std::pow(cr_Ct17, 2)*cr_Ct19 + std::pow(cr_Ct28, 2)*cr_Ct30 + cr_Ct30*std::pow(cr_Ct33, 2) + cr_Ct37);
  const double cr_Ct41 = 1.1547005383792517*threshold;
  const double cr_Ct42 = cr_Ct0*(cr_Ct1*(-1.0 + cr_Ct41/std::sqrt(cr_Ct40)) + 1.0);
  const double cr_Ct43 = -cr_Ct25;
  const double cr_Ct44 = cr_Ct4 + cr_Ct43 + cr_Ct7;
  const double cr_Ct45 = -cr_Ct3;
  const double cr_Ct46 = cr_Ct23 + cr_Ct8;
  const double cr_Ct47 = cr_Ct32 + cr_Ct5;
  const double cr_Ct48 = -cr_Ct20;
  const double cr_Ct49 = -cr_Ct21;
  const double cr_Ct50 = cr_Ct24 + cr_Ct26;
  const double cr_Ct51 = -cr_Ct22;
  const double cr_Ct52 = -cr_Ct31;
  const double cr_Ct53 = std::pow(cr_Ct39*(cr_Ct30*std::pow(cr_Ct43 + cr_Ct45 + cr_Ct46 + cr_Ct47, 2) + cr_Ct30*std::pow(cr_Ct46 + cr_Ct48 + cr_Ct49 + cr_Ct50, 2) + cr_Ct30*std::pow(cr_Ct47 + cr_Ct50 + cr_Ct51 + cr_Ct52, 2) + cr_Ct37), -3.0/2.0);
  const double cr_Ct54 = 2.0*nu;
  const double cr_Ct55 = cr_Ct54 - 1.0;
  const double cr_Ct56 = 0.25*cr_Ct29*cr_Ct55;
  const double cr_Ct57 = -cr_Ct33*cr_Ct56;
  const double cr_Ct58 = -cr_Ct28*cr_Ct56;
  const double cr_Ct59 = 4.0*nu;
  const double cr_Ct60 = cr_Ct59 - 2.0;
  const double cr_Ct61 = cr_Ct53*(0.25*cr_Ct17*cr_Ct18*cr_Ct60 + cr_Ct57 + cr_Ct58);
  const double cr_Ct62 = cr_Ct1*threshold;
  const double cr_Ct63 = cr_Ct39*cr_Ct62;
  const double cr_Ct64 = 0.5132002392796674*cr_Ct63;
  const double cr_Ct65 = 1.0/(cr_Ct54 - 1.0);
  const double cr_Ct66 = Young/cr_Ct38;
  const double cr_Ct67 = cr_Ct65*cr_Ct66;
  const double cr_Ct68 = -cr_Ct16 + cr_Ct8;
  const double cr_Ct69 = -cr_Ct14 + cr_Ct5;
  const double cr_Ct70 = cr_Ct12 + cr_Ct45 + cr_Ct68 + cr_Ct69;
  const double cr_Ct71 = -0.5*cr_Ct12 + cr_Ct24;
  const double cr_Ct72 = cr_Ct13 + cr_Ct49 + cr_Ct68 + cr_Ct71;
  const double cr_Ct73 = cr_Ct15 + cr_Ct52 + cr_Ct69 + cr_Ct71;
  const double cr_Ct74 = cr_Ct1*(cr_Ct41/std::sqrt(cr_Ct39*(cr_Ct19*std::pow(cr_Ct70, 2) + cr_Ct19*std::pow(cr_Ct72, 2) + cr_Ct19*std::pow(cr_Ct73, 2) + cr_Ct37)) - 1.0);
  const double cr_Ct75 = nu*(cr_Ct74 + 1.0);
  const double cr_Ct76 = 2.0 - cr_Ct59;
  const double cr_Ct77 = cr_Ct55*cr_Ct70;
  const double cr_Ct78 = cr_Ct55*cr_Ct73;
  const double cr_Ct79 = cr_Ct72*cr_Ct76 + cr_Ct77 + cr_Ct78;
  const double cr_Ct80 = -cr_Ct4;
  const double cr_Ct81 = -cr_Ct7;
  const double cr_Ct82 = cr_Ct63/std::pow(cr_Ct40, 3.0/2.0);
  const double cr_Ct83 = 0.12830005981991685*cr_Ct18*cr_Ct82;
  const double cr_Ct84 = cr_Ct83*(cr_Ct25 + cr_Ct80 + cr_Ct81);
  const double cr_Ct85 = cr_Ct55*cr_Ct72;
  const double cr_Ct86 = cr_Ct73*cr_Ct76 + cr_Ct77 + cr_Ct85;
  const double cr_Ct87 = std::pow(Young, 3)/std::pow(cr_Ct38, 3);
  const double cr_Ct88 = cr_Ct1*cr_Ct41*cr_Ct65*cr_Ct87;
  const double cr_Ct89 = cr_Ct44*cr_Ct88;
  const double cr_Ct90 = cr_Ct53*r_strain[3];
  const double cr_Ct91 = cr_Ct53*r_strain[4];
  const double cr_Ct92 = cr_Ct53*r_strain[5];
  const double cr_Ct93 = cr_Ct70*cr_Ct76 + cr_Ct78 + cr_Ct85;
  const double cr_Ct94 = -cr_Ct2;
  const double cr_Ct95 = cr_Ct83*(cr_Ct20 + cr_Ct81 + cr_Ct94);
  const double cr_Ct96 = cr_Ct2 + cr_Ct48 + cr_Ct7;
  const double cr_Ct97 = cr_Ct65/(1.0 - cr_Ct54);
  const double cr_Ct98 = cr_Ct60*cr_Ct97;
  const double cr_Ct99 = cr_Ct17*cr_Ct55*cr_Ct97;
  const double cr_Ct100 = -cr_Ct28*cr_Ct98 + cr_Ct57 + cr_Ct99;
  const double cr_Ct101 = cr_Ct53*cr_Ct64;
  const double cr_Ct102 = cr_Ct88*cr_Ct96;
  const double cr_Ct103 = cr_Ct83*(cr_Ct22 + cr_Ct80 + cr_Ct94);
  const double cr_Ct104 = cr_Ct2 + cr_Ct4 + cr_Ct51;
  const double cr_Ct105 = -cr_Ct33*cr_Ct98 + cr_Ct58 + cr_Ct99;
  const double cr_Ct106 = cr_Ct104*cr_Ct88;
  const double cr_Ct107 = cr_Ct62*cr_Ct87;
  const double cr_Ct108 = cr_Ct107*r_strain[3];
  const double cr_Ct109 = 0.2566001196398337*cr_Ct108;
  const double cr_Ct110 = cr_Ct109*cr_Ct53;
  const double cr_Ct111 = 0.57735026918962584*cr_Ct82;
  const double cr_Ct112 = 0.5*cr_Ct74 + 0.5;
  const double cr_Ct113 = 0.57735026918962584*cr_Ct108;
  const double cr_Ct114 = -cr_Ct113*cr_Ct91;
  const double cr_Ct115 = -cr_Ct113*cr_Ct92;
  const double cr_Ct116 = 0.2566001196398337*cr_Ct107;
  const double cr_Ct117 = cr_Ct116*cr_Ct61;
  const double cr_Ct118 = cr_Ct116*cr_Ct91;
  const double cr_Ct119 = -0.57735026918962584*cr_Ct107*cr_Ct91*r_strain[5];
  const double cr_Ct120 = cr_Ct116*cr_Ct92;
  r_Ct(0,0)=cr_Ct67*(cr_Ct42 + cr_Ct44*cr_Ct61*cr_Ct64);
  r_Ct(0,1)=-cr_Ct67*(cr_Ct75 + cr_Ct79*cr_Ct84);
  r_Ct(0,2)=-cr_Ct67*(cr_Ct75 + cr_Ct84*cr_Ct86);
  r_Ct(0,3)=cr_Ct89*cr_Ct90;
  r_Ct(0,4)=cr_Ct89*cr_Ct91;
  r_Ct(0,5)=cr_Ct89*cr_Ct92;
  r_Ct(1,0)=-cr_Ct67*(cr_Ct75 + cr_Ct93*cr_Ct95);
  r_Ct(1,1)=cr_Ct67*(cr_Ct100*cr_Ct101*cr_Ct96 + cr_Ct42);
  r_Ct(1,2)=-cr_Ct67*(cr_Ct75 + cr_Ct86*cr_Ct95);
  r_Ct(1,3)=cr_Ct102*cr_Ct90;
  r_Ct(1,4)=cr_Ct102*cr_Ct91;
  r_Ct(1,5)=cr_Ct102*cr_Ct92;
  r_Ct(2,0)=-cr_Ct67*(cr_Ct103*cr_Ct93 + cr_Ct75);
  r_Ct(2,1)=-cr_Ct67*(cr_Ct103*cr_Ct79 + cr_Ct75);
  r_Ct(2,2)=cr_Ct67*(cr_Ct101*cr_Ct104*cr_Ct105 + cr_Ct42);
  r_Ct(2,3)=cr_Ct106*cr_Ct90;
  r_Ct(2,4)=cr_Ct106*cr_Ct91;
  r_Ct(2,5)=cr_Ct106*cr_Ct92;
  r_Ct(3,0)=-cr_Ct109*cr_Ct61;
  r_Ct(3,1)=-cr_Ct100*cr_Ct110;
  r_Ct(3,2)=-cr_Ct105*cr_Ct110;
  r_Ct(3,3)=cr_Ct66*(-cr_Ct111*cr_Ct34 + cr_Ct112);
  r_Ct(3,4)=cr_Ct114;
  r_Ct(3,5)=cr_Ct115;
  r_Ct(4,0)=-cr_Ct117*r_strain[4];
  r_Ct(4,1)=-cr_Ct100*cr_Ct118;
  r_Ct(4,2)=-cr_Ct105*cr_Ct118;
  r_Ct(4,3)=cr_Ct114;
  r_Ct(4,4)=cr_Ct66*(-cr_Ct111*cr_Ct35 + cr_Ct112);
  r_Ct(4,5)=cr_Ct119;
  r_Ct(5,0)=-cr_Ct117*r_strain[5];
  r_Ct(5,1)=-cr_Ct100*cr_Ct120;
  r_Ct(5,2)=-cr_Ct105*cr_Ct120;
  r_Ct(5,3)=cr_Ct115;
  r_Ct(5,4)=cr_Ct119;
  r_Ct(5,5)=cr_Ct66*(-cr_Ct111*cr_Ct36 + cr_Ct112);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  const auto &r_props = rValues.GetMaterialProperties();
  const double Young = r_props[YOUNG_MODULUS];
  const double nu = r_props[POISSON_RATIO];
  const double Gf = r_props[FRACTURE_ENERGY];
  const double phi = r_props[FRICTION_ANGLE] * Globals::Pi / 180.0;
  const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
  auto &r_Ct = rValues.GetConstitutiveMatrix();
  const auto &r_strain = rValues.GetStrainVector();

  const bool has_symmetric_yield_stress = r_props.Has(YIELD_STRESS);
  const double threshold_compression = has_symmetric_yield_stress ? r_props[YIELD_STRESS] : r_props[YIELD_STRESS_COMPRESSION];
  const double threshold_tension = has_symmetric_yield_stress ? r_props[YIELD_STRESS] : r_props[YIELD_STRESS_TENSION];

  const double cr_Ct0 = nu - 1.0;
  const double cr_Ct1 = fabs(threshold_compression/threshold_tension);
  const double cr_Ct2 = 1.0/(1.0 - 0.5*characteristic_length*std::pow(threshold_compression, 2)/(Gf*Young*std::pow(cr_Ct1, 2)));
  const double cr_Ct3 = 2.0*nu;
  const double cr_Ct4 = 1.0/(cr_Ct3 - 1.0);
  const double cr_Ct5 = std::tan(0.5*phi + 0.78539816339744828);
  const double cr_Ct6 = cr_Ct1/std::pow(cr_Ct5, 2);
  const double cr_Ct7 = std::sin(phi);
  const double cr_Ct8 = cr_Ct6 + 1;
  const double cr_Ct9 = cr_Ct7*cr_Ct8;
  const double cr_Ct10 = cr_Ct6 + cr_Ct9 - 1;
  const double cr_Ct11 = cr_Ct10*cr_Ct4;
  const double cr_Ct12 = 0.16666666666666666*cr_Ct11;
  const double cr_Ct13 = nu*r_strain[0];
  const double cr_Ct14 = 2*cr_Ct13;
  const double cr_Ct15 = nu*r_strain[1];
  const double cr_Ct16 = 2*cr_Ct15;
  const double cr_Ct17 = nu*r_strain[2];
  const double cr_Ct18 = 2*cr_Ct17;
  const double cr_Ct19 = cr_Ct0*r_strain[0];
  const double cr_Ct20 = cr_Ct0*r_strain[1];
  const double cr_Ct21 = cr_Ct0*r_strain[2];
  const double cr_Ct22 = nu + 1.0;
  const double cr_Ct23 = 1.0/cr_Ct22;
  const double cr_Ct24 = Young*cr_Ct23;
  const double cr_Ct25 = cr_Ct24*(-cr_Ct14 - cr_Ct16 - cr_Ct18 + cr_Ct19 + cr_Ct20 + cr_Ct21);
  const double cr_Ct26 = -nu;
  const double cr_Ct27 = std::pow(cr_Ct26 + 0.5, -2);
  const double cr_Ct28 = cr_Ct13;
  const double cr_Ct29 = 0.5*cr_Ct15;
  const double cr_Ct30 = -cr_Ct29;
  const double cr_Ct31 = 0.5*cr_Ct17;
  const double cr_Ct32 = -cr_Ct31;
  const double cr_Ct33 = cr_Ct26 + 1.0;
  const double cr_Ct34 = cr_Ct33*r_strain[0];
  const double cr_Ct35 = cr_Ct33*r_strain[1];
  const double cr_Ct36 = 0.5*cr_Ct35;
  const double cr_Ct37 = cr_Ct33*r_strain[2];
  const double cr_Ct38 = 0.5*cr_Ct37;
  const double cr_Ct39 = cr_Ct28 + cr_Ct30 + cr_Ct32 - cr_Ct34 + cr_Ct36 + cr_Ct38;
  const double cr_Ct40 = cr_Ct27*std::pow(cr_Ct39, 2);
  const double cr_Ct41 = cr_Ct15;
  const double cr_Ct42 = 0.5*cr_Ct21;
  const double cr_Ct43 = 0.5*cr_Ct13;
  const double cr_Ct44 = 0.5*cr_Ct19;
  const double cr_Ct45 = -cr_Ct43 - cr_Ct44;
  const double cr_Ct46 = cr_Ct20 + cr_Ct32 + cr_Ct41 - cr_Ct42 + cr_Ct45;
  const double cr_Ct47 = std::pow(cr_Ct46, 2);
  const double cr_Ct48 = std::pow(nu - 0.5, -2);
  const double cr_Ct49 = 0.22222222222222227*cr_Ct48;
  const double cr_Ct50 = cr_Ct17;
  const double cr_Ct51 = 0.5*cr_Ct20;
  const double cr_Ct52 = cr_Ct21 + cr_Ct30 + cr_Ct45 + cr_Ct50 - cr_Ct51;
  const double cr_Ct53 = std::pow(cr_Ct52, 2);
  const double cr_Ct54 = std::pow(r_strain[3], 2);
  const double cr_Ct55 = std::pow(r_strain[4], 2);
  const double cr_Ct56 = std::pow(r_strain[5], 2);
  const double cr_Ct57 = cr_Ct54 + cr_Ct55 + cr_Ct56;
  const double cr_Ct58 = std::pow(cr_Ct22, -2);
  const double cr_Ct59 = std::pow(Young, 2)*cr_Ct58;
  const double cr_Ct60 = std::sqrt(cr_Ct59*(0.22222222222222227*cr_Ct40 + cr_Ct47*cr_Ct49 + cr_Ct49*cr_Ct53 + cr_Ct57));
  const double cr_Ct61 = 1.0/cr_Ct60;
  const double cr_Ct62 = 0.11111111111111113*cr_Ct48;
  const double cr_Ct63 = 0.5*cr_Ct54 + 0.5*cr_Ct55 + 0.5*cr_Ct56;
  const double cr_Ct64 = 0.25*r_strain[4];
  const double cr_Ct65 = cr_Ct64*r_strain[5];
  const double cr_Ct66 = 0.5*r_strain[3];
  const double cr_Ct67 = 1.0/(1.0 - cr_Ct3);
  const double cr_Ct68 = 0.66666666666666663*cr_Ct17;
  const double cr_Ct69 = 0.33333333333333337*cr_Ct13;
  const double cr_Ct70 = -0.33333333333333331*cr_Ct34 + cr_Ct69;
  const double cr_Ct71 = 0.33333333333333337*cr_Ct15;
  const double cr_Ct72 = 0.33333333333333331*cr_Ct35;
  const double cr_Ct73 = cr_Ct71 - cr_Ct72;
  const double cr_Ct74 = 0.66666666666666674*cr_Ct37 - cr_Ct68 + cr_Ct70 + cr_Ct73;
  const double cr_Ct75 = cr_Ct67*cr_Ct74;
  const double cr_Ct76 = cr_Ct65 - cr_Ct66*cr_Ct75;
  const double cr_Ct77 = 2.598076211353316*r_strain[3];
  const double cr_Ct78 = cr_Ct64*r_strain[3];
  const double cr_Ct79 = 0.5*r_strain[5];
  const double cr_Ct80 = 0.66666666666666663*cr_Ct15;
  const double cr_Ct81 = 0.33333333333333337*cr_Ct17;
  const double cr_Ct82 = 0.33333333333333331*cr_Ct37;
  const double cr_Ct83 = cr_Ct81 - cr_Ct82;
  const double cr_Ct84 = 0.66666666666666674*cr_Ct35 + cr_Ct70 - cr_Ct80 + cr_Ct83;
  const double cr_Ct85 = cr_Ct67*cr_Ct84;
  const double cr_Ct86 = cr_Ct78 - cr_Ct79*cr_Ct85;
  const double cr_Ct87 = 2.598076211353316*r_strain[5];
  const double cr_Ct88 = 0.25*cr_Ct55;
  const double cr_Ct89 = 0.25*cr_Ct27;
  const double cr_Ct90 = cr_Ct74*cr_Ct84;
  const double cr_Ct91 = cr_Ct89*cr_Ct90;
  const double cr_Ct92 = 0.66666666666666674*cr_Ct34;
  const double cr_Ct93 = 0.66666666666666663*cr_Ct13;
  const double cr_Ct94 = cr_Ct73 + cr_Ct83 + cr_Ct92 - cr_Ct93;
  const double cr_Ct95 = cr_Ct67*cr_Ct94;
  const double cr_Ct96 = cr_Ct95*(-cr_Ct88 + cr_Ct91);
  const double cr_Ct97 = cr_Ct76*cr_Ct77 + cr_Ct86*cr_Ct87 + 5.196152422706632*cr_Ct96;
  const double cr_Ct98 = 2.0*cr_Ct97;
  const double cr_Ct99 = cr_Ct24*cr_Ct98;
  const double cr_Ct100 = 0.33333333333333331*std::asin(cr_Ct61*cr_Ct99/(0.11111111111111113*cr_Ct40 + cr_Ct47*cr_Ct62 + cr_Ct53*cr_Ct62 + cr_Ct63));
  const double cr_Ct101 = std::cos(cr_Ct100);
  const double cr_Ct102 = 1 - cr_Ct6;
  const double cr_Ct103 = -cr_Ct102*cr_Ct7 + cr_Ct8;
  const double cr_Ct104 = 0.5*cr_Ct103;
  const double cr_Ct105 = std::sin(cr_Ct100);
  const double cr_Ct106 = cr_Ct7*(-cr_Ct102/cr_Ct7 + cr_Ct8);
  const double cr_Ct107 = 0.28867513459481292*cr_Ct106;
  const double cr_Ct108 = cr_Ct101*cr_Ct104 + cr_Ct105*cr_Ct107;
  const double cr_Ct109 = cr_Ct108*cr_Ct60;
  const double cr_Ct110 = threshold_compression*std::cos(phi)/cr_Ct5;
  const double cr_Ct111 = 0.5*cr_Ct110;
  const double cr_Ct112 = cr_Ct0*(cr_Ct2*(cr_Ct111/(0.5*cr_Ct109 + cr_Ct12*cr_Ct25) - 1.0) + 1.0);
  const double cr_Ct113 = -cr_Ct19;
  const double cr_Ct114 = cr_Ct113 + cr_Ct15 + cr_Ct17;
  const double cr_Ct115 = cr_Ct12*(cr_Ct26 - 1.0);
  const double cr_Ct116 = cr_Ct3 - 1.0;
  const double cr_Ct117 = 0.25*cr_Ct116*cr_Ct48;
  const double cr_Ct118 = -cr_Ct117*cr_Ct52;
  const double cr_Ct119 = -cr_Ct117*cr_Ct46;
  const double cr_Ct120 = 4.0*nu;
  const double cr_Ct121 = cr_Ct120 - 2.0;
  const double cr_Ct122 = cr_Ct108*cr_Ct61;
  const double cr_Ct123 = cr_Ct122*cr_Ct24;
  const double cr_Ct124 = 0.22222222222222227*cr_Ct123;
  const double cr_Ct125 = 0.66666666666666674*nu;
  const double cr_Ct126 = cr_Ct125 - 0.33333333333333331;
  const double cr_Ct127 = 2.598076211353316*cr_Ct126;
  const double cr_Ct128 = cr_Ct127*cr_Ct54;
  const double cr_Ct129 = cr_Ct127*cr_Ct56;
  const double cr_Ct130 = cr_Ct27*cr_Ct94;
  const double cr_Ct131 = 0.66666666666666674 - 1.3333333333333335*nu;
  const double cr_Ct132 = -2.598076211353316*cr_Ct27*cr_Ct90 + 2.598076211353316*cr_Ct55;
  const double cr_Ct133 = -cr_Ct28;
  const double cr_Ct134 = cr_Ct31 - cr_Ct38;
  const double cr_Ct135 = cr_Ct29 - cr_Ct36;
  const double cr_Ct136 = cr_Ct133 + cr_Ct134 + cr_Ct135 + cr_Ct34;
  const double cr_Ct137 = std::pow(cr_Ct136, 2);
  const double cr_Ct138 = 0.11111111111111113*cr_Ct27;
  const double cr_Ct139 = -cr_Ct41;
  const double cr_Ct140 = -0.5*cr_Ct34 + cr_Ct43;
  const double cr_Ct141 = cr_Ct134 + cr_Ct139 + cr_Ct140 + cr_Ct35;
  const double cr_Ct142 = std::pow(cr_Ct141, 2);
  const double cr_Ct143 = -cr_Ct50;
  const double cr_Ct144 = cr_Ct135 + cr_Ct140 + cr_Ct143 + cr_Ct37;
  const double cr_Ct145 = std::pow(cr_Ct144, 2);
  const double cr_Ct146 = 1.0/(cr_Ct137*cr_Ct138 + cr_Ct138*cr_Ct142 + cr_Ct138*cr_Ct145 + cr_Ct63);
  const double cr_Ct147 = cr_Ct146*cr_Ct67;
  const double cr_Ct148 = 2.0 - cr_Ct120;
  const double cr_Ct149 = cr_Ct116*cr_Ct141;
  const double cr_Ct150 = cr_Ct116*cr_Ct144;
  const double cr_Ct151 = cr_Ct136*cr_Ct148 + cr_Ct149 + cr_Ct150;
  const double cr_Ct152 = cr_Ct151*cr_Ct27;
  const double cr_Ct153 = 0.22222222222222227*cr_Ct27;
  const double cr_Ct154 = cr_Ct137*cr_Ct153 + cr_Ct142*cr_Ct153 + cr_Ct145*cr_Ct153 + cr_Ct57;
  const double cr_Ct155 = cr_Ct97/std::pow(cr_Ct154, 2);
  const double cr_Ct156 = 0.88888888888888906*cr_Ct155;
  const double cr_Ct157 = cr_Ct146/cr_Ct154;
  const double cr_Ct158 = cr_Ct153*cr_Ct157*cr_Ct97;
  const double cr_Ct159 = cr_Ct147*(-cr_Ct127*cr_Ct130*(cr_Ct125*r_strain[0] - 0.33333333333333326*cr_Ct15 - 0.33333333333333326*cr_Ct17 - 0.66666666666666663*cr_Ct34 + 0.33333333333333343*cr_Ct35 + 0.33333333333333343*cr_Ct37) + cr_Ct128 + cr_Ct129 + cr_Ct131*cr_Ct132) + cr_Ct151*cr_Ct158 + cr_Ct152*cr_Ct156;
  const double cr_Ct160 = std::pow(0.0023148148148148147 - std::pow(cr_Ct66*cr_Ct76 + cr_Ct79*cr_Ct86 + cr_Ct96, 2)/std::pow(cr_Ct154, 3), -1.0/2.0);
  const double cr_Ct161 = std::sqrt(cr_Ct154*cr_Ct59);
  const double cr_Ct162 = 1.0/cr_Ct161;
  const double cr_Ct163 = 0.0080187537387448014*cr_Ct103;
  const double cr_Ct164 = 0.0046296296296296294*cr_Ct106;
  const double cr_Ct165 = cr_Ct160*cr_Ct162*cr_Ct60*(-cr_Ct101*cr_Ct164 + cr_Ct105*cr_Ct163);
  const double cr_Ct166 = 0.5*cr_Ct165;
  const double cr_Ct167 = cr_Ct115 + cr_Ct124*(cr_Ct118 + cr_Ct119 + cr_Ct121*cr_Ct39*cr_Ct89) + cr_Ct159*cr_Ct166;
  const double cr_Ct168 = 2.0*cr_Ct24;
  const double cr_Ct169 = -cr_Ct20;
  const double cr_Ct170 = -cr_Ct21;
  const double cr_Ct171 = cr_Ct14 + cr_Ct16 + cr_Ct18;
  const double cr_Ct172 = 0.33333333333333331*cr_Ct11;
  const double cr_Ct173 = cr_Ct31 + cr_Ct42;
  const double cr_Ct174 = cr_Ct29 + cr_Ct51;
  const double cr_Ct175 = std::pow(cr_Ct113 + cr_Ct133 + cr_Ct173 + cr_Ct174, 2);
  const double cr_Ct176 = cr_Ct43 + cr_Ct44;
  const double cr_Ct177 = std::pow(cr_Ct139 + cr_Ct169 + cr_Ct173 + cr_Ct176, 2);
  const double cr_Ct178 = std::pow(cr_Ct143 + cr_Ct170 + cr_Ct174 + cr_Ct176, 2);
  const double cr_Ct179 = std::sqrt(cr_Ct59*(cr_Ct175*cr_Ct49 + cr_Ct177*cr_Ct49 + cr_Ct178*cr_Ct49 + cr_Ct57));
  const double cr_Ct180 = -cr_Ct71;
  const double cr_Ct181 = -0.33333333333333331*cr_Ct19 - cr_Ct69;
  const double cr_Ct182 = -cr_Ct81;
  const double cr_Ct183 = 0.33333333333333331*std::asin(cr_Ct168*(5.196152422706632*cr_Ct67*(cr_Ct88 - cr_Ct91)*(cr_Ct180 + cr_Ct182 + cr_Ct72 + cr_Ct82 - cr_Ct92 + cr_Ct93) + cr_Ct77*(-cr_Ct4*cr_Ct66*(cr_Ct180 + cr_Ct181 - 0.33333333333333331*cr_Ct20 + 0.66666666666666674*cr_Ct21 + cr_Ct68) + cr_Ct65) + cr_Ct87*(-cr_Ct4*cr_Ct79*(cr_Ct181 + cr_Ct182 + 0.66666666666666674*cr_Ct20 - 0.33333333333333331*cr_Ct21 + cr_Ct80) + cr_Ct78))/(cr_Ct179*(cr_Ct175*cr_Ct62 + cr_Ct177*cr_Ct62 + cr_Ct178*cr_Ct62 + cr_Ct63)));
  const double cr_Ct184 = std::pow(-cr_Ct172*cr_Ct24*(cr_Ct113 + cr_Ct169 + cr_Ct170 + cr_Ct171) + cr_Ct179*(cr_Ct104*std::cos(cr_Ct183) + cr_Ct107*std::sin(cr_Ct183)), -2);
  const double cr_Ct185 = cr_Ct110*cr_Ct2;
  const double cr_Ct186 = cr_Ct184*cr_Ct185;
  const double cr_Ct187 = cr_Ct168*cr_Ct186;
  const double cr_Ct188 = cr_Ct24*cr_Ct4;
  const double cr_Ct189 = 0.33333333333333331*std::asin(cr_Ct146*cr_Ct162*cr_Ct99);
  const double cr_Ct190 = std::cos(cr_Ct189);
  const double cr_Ct191 = std::sin(cr_Ct189);
  const double cr_Ct192 = cr_Ct104*cr_Ct190 + cr_Ct107*cr_Ct191;
  const double cr_Ct193 = cr_Ct2*(cr_Ct111/(0.16666666666666666*cr_Ct10*cr_Ct24*cr_Ct67*(cr_Ct171 + cr_Ct34 + cr_Ct35 + cr_Ct37) + 0.5*cr_Ct161*cr_Ct192) - 1.0) + 1.0;
  const double cr_Ct194 = cr_Ct193*nu;
  const double cr_Ct195 = cr_Ct67*(0.16666666666666666*cr_Ct6 + 0.16666666666666666*cr_Ct9 - 0.16666666666666666);
  const double cr_Ct196 = cr_Ct116*cr_Ct136;
  const double cr_Ct197 = cr_Ct141*cr_Ct148 + cr_Ct150 + cr_Ct196;
  const double cr_Ct198 = cr_Ct162*cr_Ct192;
  const double cr_Ct199 = 0.055555555555555566*Young*cr_Ct198*cr_Ct58;
  const double cr_Ct200 = cr_Ct199*cr_Ct27;
  const double cr_Ct201 = 2.598076211353316*cr_Ct131;
  const double cr_Ct202 = cr_Ct126*cr_Ct132;
  const double cr_Ct203 = 2.598076211353316*cr_Ct130;
  const double cr_Ct204 = cr_Ct156*cr_Ct27;
  const double cr_Ct205 = 0.5*cr_Ct147*(cr_Ct128 + cr_Ct201*cr_Ct56 + cr_Ct202 - cr_Ct203*(cr_Ct126*cr_Ct84 + cr_Ct131*cr_Ct74)) + 0.5*cr_Ct158*cr_Ct197 + 0.5*cr_Ct197*cr_Ct204;
  const double cr_Ct206 = cr_Ct160*(cr_Ct163*cr_Ct191 - cr_Ct164*cr_Ct190);
  const double cr_Ct207 = cr_Ct206*cr_Ct23;
  const double cr_Ct208 = cr_Ct195 + cr_Ct197*cr_Ct200 + cr_Ct205*cr_Ct207;
  const double cr_Ct209 = -cr_Ct15;
  const double cr_Ct210 = -cr_Ct17;
  const double cr_Ct211 = cr_Ct185/std::pow(cr_Ct109 + cr_Ct172*cr_Ct25, 2);
  const double cr_Ct212 = 2.0*Young*cr_Ct211;
  const double cr_Ct213 = cr_Ct212*(cr_Ct19 + cr_Ct209 + cr_Ct210);
  const double cr_Ct214 = cr_Ct144*cr_Ct148 + cr_Ct149 + cr_Ct196;
  const double cr_Ct215 = cr_Ct147*(cr_Ct129 + cr_Ct201*cr_Ct54 + cr_Ct202 - cr_Ct203*(cr_Ct126*cr_Ct74 + cr_Ct131*cr_Ct84)) + cr_Ct158*cr_Ct214 + cr_Ct204*cr_Ct214;
  const double cr_Ct216 = 0.5*cr_Ct207;
  const double cr_Ct217 = cr_Ct195 + cr_Ct200*cr_Ct214 + cr_Ct215*cr_Ct216;
  const double cr_Ct218 = cr_Ct24*r_strain[3];
  const double cr_Ct219 = 8.0*cr_Ct155;
  const double cr_Ct220 = cr_Ct157*cr_Ct98;
  const double cr_Ct221 = cr_Ct146*(5.196152422706632*cr_Ct75*r_strain[3] - cr_Ct87*r_strain[4]) + cr_Ct219*r_strain[3] + cr_Ct220*r_strain[3];
  const double cr_Ct222 = cr_Ct122*cr_Ct218 + cr_Ct165*cr_Ct221;
  const double cr_Ct223 = cr_Ct186*cr_Ct59;
  const double cr_Ct224 = cr_Ct223;
  const double cr_Ct225 = cr_Ct224*cr_Ct4;
  const double cr_Ct226 = cr_Ct114*cr_Ct225;
  const double cr_Ct227 = cr_Ct146*(-cr_Ct77*r_strain[5] + 5.196152422706632*cr_Ct95*r_strain[4]) + cr_Ct219*r_strain[4] + cr_Ct220*r_strain[4];
  const double cr_Ct228 = cr_Ct123*r_strain[4] + cr_Ct165*cr_Ct227;
  const double cr_Ct229 = cr_Ct146*(-cr_Ct77*r_strain[4] + 5.196152422706632*cr_Ct85*r_strain[5]) + cr_Ct219*r_strain[5] + cr_Ct220*r_strain[5];
  const double cr_Ct230 = cr_Ct123*r_strain[5] + cr_Ct165*cr_Ct229;
  const double cr_Ct231 = cr_Ct152*cr_Ct199 + cr_Ct159*cr_Ct216 + cr_Ct195;
  const double cr_Ct232 = -cr_Ct13;
  const double cr_Ct233 = cr_Ct212*(cr_Ct20 + cr_Ct210 + cr_Ct232);
  const double cr_Ct234 = cr_Ct13 + cr_Ct169 + cr_Ct17;
  const double cr_Ct235 = cr_Ct4*cr_Ct67;
  const double cr_Ct236 = cr_Ct121*cr_Ct235;
  const double cr_Ct237 = cr_Ct116*cr_Ct235*cr_Ct39;
  const double cr_Ct238 = cr_Ct115 + cr_Ct124*(cr_Ct118 - cr_Ct236*cr_Ct46 + cr_Ct237) + cr_Ct165*cr_Ct205;
  const double cr_Ct239 = cr_Ct225*cr_Ct234;
  const double cr_Ct240 = cr_Ct212*(cr_Ct209 + cr_Ct21 + cr_Ct232);
  const double cr_Ct241 = cr_Ct13 + cr_Ct15 + cr_Ct170;
  const double cr_Ct242 = cr_Ct115 + cr_Ct124*(cr_Ct119 - cr_Ct236*cr_Ct52 + cr_Ct237) + cr_Ct166*cr_Ct215;
  const double cr_Ct243 = cr_Ct225*cr_Ct241;
  const double cr_Ct244 = cr_Ct224*r_strain[3];
  const double cr_Ct245 = 0.5*cr_Ct24;
  const double cr_Ct246 = cr_Ct223*cr_Ct66;
  const double cr_Ct247 = cr_Ct224*r_strain[4];
  const double cr_Ct248 = cr_Ct111*cr_Ct184*cr_Ct2*cr_Ct59*r_strain[4];
  const double cr_Ct249 = cr_Ct198*cr_Ct24;
  const double cr_Ct250 = cr_Ct211*cr_Ct24;
  const double cr_Ct251 = cr_Ct224*r_strain[5];
  const double cr_Ct252 = cr_Ct223*cr_Ct79;
  r_Ct(0,0)=cr_Ct188*(cr_Ct112 + cr_Ct114*cr_Ct167*cr_Ct187);
  r_Ct(0,1)=-cr_Ct188*(cr_Ct194 + cr_Ct208*cr_Ct213);
  r_Ct(0,2)=-cr_Ct188*(cr_Ct194 + cr_Ct213*cr_Ct217);
  r_Ct(0,3)=cr_Ct222*cr_Ct226;
  r_Ct(0,4)=cr_Ct226*cr_Ct228;
  r_Ct(0,5)=cr_Ct226*cr_Ct230;
  r_Ct(1,0)=-cr_Ct188*(cr_Ct194 + cr_Ct231*cr_Ct233);
  r_Ct(1,1)=cr_Ct188*(cr_Ct112 + cr_Ct187*cr_Ct234*cr_Ct238);
  r_Ct(1,2)=-cr_Ct188*(cr_Ct194 + cr_Ct217*cr_Ct233);
  r_Ct(1,3)=cr_Ct222*cr_Ct239;
  r_Ct(1,4)=cr_Ct228*cr_Ct239;
  r_Ct(1,5)=cr_Ct230*cr_Ct239;
  r_Ct(2,0)=-cr_Ct188*(cr_Ct194 + cr_Ct231*cr_Ct240);
  r_Ct(2,1)=-cr_Ct188*(cr_Ct194 + cr_Ct208*cr_Ct240);
  r_Ct(2,2)=cr_Ct188*(cr_Ct112 + cr_Ct187*cr_Ct241*cr_Ct242);
  r_Ct(2,3)=cr_Ct222*cr_Ct243;
  r_Ct(2,4)=cr_Ct228*cr_Ct243;
  r_Ct(2,5)=cr_Ct230*cr_Ct243;
  r_Ct(3,0)=-cr_Ct167*cr_Ct244;
  r_Ct(3,1)=-cr_Ct238*cr_Ct244;
  r_Ct(3,2)=-cr_Ct242*cr_Ct244;
  r_Ct(3,3)=cr_Ct245*(cr_Ct193 - cr_Ct211*cr_Ct218*(cr_Ct198*cr_Ct218 + cr_Ct206*cr_Ct221));
  r_Ct(3,4)=-cr_Ct228*cr_Ct246;
  r_Ct(3,5)=-cr_Ct230*cr_Ct246;
  r_Ct(4,0)=-cr_Ct167*cr_Ct247;
  r_Ct(4,1)=-cr_Ct238*cr_Ct247;
  r_Ct(4,2)=-cr_Ct242*cr_Ct247;
  r_Ct(4,3)=-cr_Ct222*cr_Ct248;
  r_Ct(4,4)=cr_Ct245*(cr_Ct193 - cr_Ct250*r_strain[4]*(cr_Ct206*cr_Ct227 + cr_Ct249*r_strain[4]));
  r_Ct(4,5)=-cr_Ct230*cr_Ct248;
  r_Ct(5,0)=-cr_Ct167*cr_Ct251;
  r_Ct(5,1)=-cr_Ct238*cr_Ct251;
  r_Ct(5,2)=-cr_Ct242*cr_Ct251;
  r_Ct(5,3)=-cr_Ct222*cr_Ct252;
  r_Ct(5,4)=-cr_Ct228*cr_Ct252;
  r_Ct(5,5)=cr_Ct245*(cr_Ct193 - cr_Ct250*r_strain[5]*(cr_Ct206*cr_Ct229 + cr_Ct249*r_strain[5]));
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<TrescaYieldSurface<TrescaPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  const auto &r_props = rValues.GetMaterialProperties();
  const double Young = r_props[YOUNG_MODULUS];
  const double nu = r_props[POISSON_RATIO];
  const double Gf = r_props[FRACTURE_ENERGY];
  const double phi = r_props[FRICTION_ANGLE];
  const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
  const double threshold = r_props[YIELD_STRESS];
  auto &r_Ct = rValues.GetConstitutiveMatrix();
  const auto &r_strain = rValues.GetStrainVector();
  const double sin_phi = std::sin(phi * Globals::Pi / 180.0);

  const double cr_Ct0 = nu - 1.0;
  const double cr_Ct1 = 1.0/(1.0 - 0.5*characteristic_length*std::pow(threshold, 2)/(Gf*Young));
  const double cr_Ct2 = 4.0*nu;
  const double cr_Ct3 = cr_Ct0*r_strain[0];
  const double cr_Ct4 = cr_Ct0*r_strain[1];
  const double cr_Ct5 = cr_Ct0*r_strain[2];
  const double cr_Ct6 = nu + 1.0;
  const double cr_Ct7 = Young/cr_Ct6;
  const double cr_Ct8 = 1.7320508075688772*sin_phi;
  const double cr_Ct9 = 1.0/(cr_Ct8 - 5.196152422706632);
  const double cr_Ct10 = 2.0*nu;
  const double cr_Ct11 = 1.0/(cr_Ct10 - 1.0);
  const double cr_Ct12 = cr_Ct11*cr_Ct9*sin_phi;
  const double cr_Ct13 = cr_Ct12*cr_Ct7*(-cr_Ct2*r_strain[0] - cr_Ct2*r_strain[1] - cr_Ct2*r_strain[2] + 2.0*cr_Ct3 + 2.0*cr_Ct4 + 2.0*cr_Ct5);
  const double cr_Ct14 = nu*r_strain[0];
  const double cr_Ct15 = cr_Ct14;
  const double cr_Ct16 = nu*r_strain[1];
  const double cr_Ct17 = 0.5*cr_Ct16;
  const double cr_Ct18 = -cr_Ct17;
  const double cr_Ct19 = nu*r_strain[2];
  const double cr_Ct20 = 0.5*cr_Ct19;
  const double cr_Ct21 = -cr_Ct20;
  const double cr_Ct22 = -nu;
  const double cr_Ct23 = cr_Ct22 + 1.0;
  const double cr_Ct24 = cr_Ct23*r_strain[0];
  const double cr_Ct25 = cr_Ct23*r_strain[1];
  const double cr_Ct26 = 0.5*cr_Ct25;
  const double cr_Ct27 = cr_Ct23*r_strain[2];
  const double cr_Ct28 = 0.5*cr_Ct27;
  const double cr_Ct29 = cr_Ct15 + cr_Ct18 + cr_Ct21 - cr_Ct24 + cr_Ct26 + cr_Ct28;
  const double cr_Ct30 = std::pow(cr_Ct22 + 0.5, -2);
  const double cr_Ct31 = 0.22222222222222227*cr_Ct30;
  const double cr_Ct32 = cr_Ct16;
  const double cr_Ct33 = 0.5*cr_Ct5;
  const double cr_Ct34 = 0.5*cr_Ct14;
  const double cr_Ct35 = 0.5*cr_Ct3;
  const double cr_Ct36 = -cr_Ct34 - cr_Ct35;
  const double cr_Ct37 = cr_Ct21 + cr_Ct32 - cr_Ct33 + cr_Ct36 + cr_Ct4;
  const double cr_Ct38 = std::pow(nu - 0.5, -2);
  const double cr_Ct39 = 0.22222222222222227*cr_Ct38;
  const double cr_Ct40 = cr_Ct19;
  const double cr_Ct41 = 0.5*cr_Ct4;
  const double cr_Ct42 = cr_Ct18 + cr_Ct36 + cr_Ct40 - cr_Ct41 + cr_Ct5;
  const double cr_Ct43 = std::pow(r_strain[3], 2);
  const double cr_Ct44 = std::pow(r_strain[4], 2);
  const double cr_Ct45 = std::pow(r_strain[5], 2);
  const double cr_Ct46 = cr_Ct43 + cr_Ct44 + cr_Ct45;
  const double cr_Ct47 = std::pow(Young, 2)/std::pow(cr_Ct6, 2);
  const double cr_Ct48 = std::sqrt(cr_Ct47*(std::pow(cr_Ct29, 2)*cr_Ct31 + std::pow(cr_Ct37, 2)*cr_Ct39 + cr_Ct39*std::pow(cr_Ct42, 2) + cr_Ct46));
  const double cr_Ct49 = 0.5*cr_Ct48;
  const double cr_Ct50 = cr_Ct13 - cr_Ct49;
  const double cr_Ct51 = sin_phi - 1;
  const double cr_Ct52 = cr_Ct51*cr_Ct9*std::abs(threshold*(sin_phi + 3.0)/cr_Ct51);
  const double cr_Ct53 = cr_Ct0*(-cr_Ct1*(1 + cr_Ct52/cr_Ct50) + 1);
  const double cr_Ct54 = -cr_Ct10;
  const double cr_Ct55 = -cr_Ct12*(cr_Ct54 - 2.0);
  const double cr_Ct56 = cr_Ct10 - 1.0;
  const double cr_Ct57 = 0.25*cr_Ct38*cr_Ct56;
  const double cr_Ct58 = -cr_Ct42*cr_Ct57;
  const double cr_Ct59 = -cr_Ct37*cr_Ct57;
  const double cr_Ct60 = cr_Ct2 - 2.0;
  const double cr_Ct61 = 1.0/cr_Ct48;
  const double cr_Ct62 = 0.22222222222222227*cr_Ct61*cr_Ct7;
  const double cr_Ct63 = cr_Ct55 + cr_Ct62*(0.25*cr_Ct29*cr_Ct30*cr_Ct60 + cr_Ct58 + cr_Ct59);
  const double cr_Ct64 = -cr_Ct3;
  const double cr_Ct65 = cr_Ct16 + cr_Ct19 + cr_Ct64;
  const double cr_Ct66 = std::pow(cr_Ct50, -2);
  const double cr_Ct67 = cr_Ct1*cr_Ct52;
  const double cr_Ct68 = cr_Ct67*cr_Ct7;
  const double cr_Ct69 = cr_Ct66*cr_Ct68;
  const double cr_Ct70 = cr_Ct11*cr_Ct7;
  const double cr_Ct71 = -cr_Ct13 + cr_Ct49;
  const double cr_Ct72 = cr_Ct1*(-cr_Ct52/cr_Ct71 + 1);
  const double cr_Ct73 = nu*(1 - cr_Ct72);
  const double cr_Ct74 = 1.0/(cr_Ct54 + 1.0);
  const double cr_Ct75 = cr_Ct74*sin_phi*(cr_Ct10 + 2.0)/(5.196152422706632 - cr_Ct8);
  const double cr_Ct76 = 2.0 - cr_Ct2;
  const double cr_Ct77 = -cr_Ct32;
  const double cr_Ct78 = cr_Ct20 - cr_Ct28;
  const double cr_Ct79 = -0.5*cr_Ct24 + cr_Ct34;
  const double cr_Ct80 = cr_Ct25 + cr_Ct77 + cr_Ct78 + cr_Ct79;
  const double cr_Ct81 = -cr_Ct15;
  const double cr_Ct82 = cr_Ct17 - cr_Ct26;
  const double cr_Ct83 = cr_Ct24 + cr_Ct78 + cr_Ct81 + cr_Ct82;
  const double cr_Ct84 = cr_Ct56*cr_Ct83;
  const double cr_Ct85 = -cr_Ct40;
  const double cr_Ct86 = cr_Ct27 + cr_Ct79 + cr_Ct82 + cr_Ct85;
  const double cr_Ct87 = cr_Ct56*cr_Ct86;
  const double cr_Ct88 = 0.055555555555555566*cr_Ct30*cr_Ct7/std::sqrt(cr_Ct47*(cr_Ct31*std::pow(cr_Ct80, 2) + cr_Ct31*std::pow(cr_Ct83, 2) + cr_Ct31*std::pow(cr_Ct86, 2) + cr_Ct46));
  const double cr_Ct89 = cr_Ct75 + cr_Ct88*(cr_Ct76*cr_Ct80 + cr_Ct84 + cr_Ct87);
  const double cr_Ct90 = -cr_Ct16;
  const double cr_Ct91 = -cr_Ct19;
  const double cr_Ct92 = std::pow(cr_Ct71, -2);
  const double cr_Ct93 = cr_Ct68*cr_Ct92;
  const double cr_Ct94 = cr_Ct93*(cr_Ct3 + cr_Ct90 + cr_Ct91);
  const double cr_Ct95 = cr_Ct56*cr_Ct80;
  const double cr_Ct96 = cr_Ct75 + cr_Ct88*(cr_Ct76*cr_Ct86 + cr_Ct84 + cr_Ct95);
  const double cr_Ct97 = cr_Ct66*r_strain[3];
  const double cr_Ct98 = cr_Ct20 + cr_Ct33;
  const double cr_Ct99 = cr_Ct17 + cr_Ct41;
  const double cr_Ct100 = -cr_Ct4;
  const double cr_Ct101 = cr_Ct34 + cr_Ct35;
  const double cr_Ct102 = -cr_Ct5;
  const double cr_Ct103 = std::pow(Young, 3)*cr_Ct67/(std::pow(cr_Ct6, 3)*std::sqrt(cr_Ct47*(cr_Ct39*std::pow(cr_Ct100 + cr_Ct101 + cr_Ct77 + cr_Ct98, 2) + cr_Ct39*std::pow(cr_Ct101 + cr_Ct102 + cr_Ct85 + cr_Ct99, 2) + cr_Ct39*std::pow(cr_Ct64 + cr_Ct81 + cr_Ct98 + cr_Ct99, 2) + cr_Ct46)));
  const double cr_Ct104 = cr_Ct103*cr_Ct97;
  const double cr_Ct105 = 0.5*cr_Ct11;
  const double cr_Ct106 = cr_Ct105*cr_Ct65;
  const double cr_Ct107 = cr_Ct103*cr_Ct66*r_strain[4];
  const double cr_Ct108 = cr_Ct103*cr_Ct66*r_strain[5];
  const double cr_Ct109 = cr_Ct75 + cr_Ct88*(cr_Ct76*cr_Ct83 + cr_Ct87 + cr_Ct95);
  const double cr_Ct110 = -cr_Ct14;
  const double cr_Ct111 = cr_Ct93*(cr_Ct110 + cr_Ct4 + cr_Ct91);
  const double cr_Ct112 = cr_Ct11*cr_Ct74;
  const double cr_Ct113 = cr_Ct112*cr_Ct60;
  const double cr_Ct114 = cr_Ct112*cr_Ct29*cr_Ct56;
  const double cr_Ct115 = cr_Ct55 + cr_Ct62*(-cr_Ct113*cr_Ct37 + cr_Ct114 + cr_Ct58);
  const double cr_Ct116 = cr_Ct100 + cr_Ct14 + cr_Ct19;
  const double cr_Ct117 = cr_Ct105*cr_Ct116;
  const double cr_Ct118 = cr_Ct93*(cr_Ct110 + cr_Ct5 + cr_Ct90);
  const double cr_Ct119 = cr_Ct55 + cr_Ct62*(-cr_Ct113*cr_Ct42 + cr_Ct114 + cr_Ct59);
  const double cr_Ct120 = cr_Ct102 + cr_Ct14 + cr_Ct16;
  const double cr_Ct121 = cr_Ct105*cr_Ct120;
  const double cr_Ct122 = cr_Ct47*cr_Ct67;
  const double cr_Ct123 = 0.5*cr_Ct122;
  const double cr_Ct124 = cr_Ct123*cr_Ct97;
  const double cr_Ct125 = 0.25*cr_Ct122*cr_Ct61*cr_Ct92;
  const double cr_Ct126 = 0.5 - 0.5*cr_Ct72;
  const double cr_Ct127 = 0.25*cr_Ct104;
  const double cr_Ct128 = -cr_Ct127*r_strain[4];
  const double cr_Ct129 = -cr_Ct127*r_strain[5];
  const double cr_Ct130 = cr_Ct123*cr_Ct66;
  const double cr_Ct131 = cr_Ct130*r_strain[4];
  const double cr_Ct132 = -0.25*cr_Ct107*r_strain[5];
  const double cr_Ct133 = cr_Ct130*r_strain[5];
  r_Ct(0,0)=cr_Ct70*(cr_Ct53 + cr_Ct63*cr_Ct65*cr_Ct69);
  r_Ct(0,1)=-cr_Ct70*(cr_Ct73 + cr_Ct89*cr_Ct94);
  r_Ct(0,2)=-cr_Ct70*(cr_Ct73 + cr_Ct94*cr_Ct96);
  r_Ct(0,3)=cr_Ct104*cr_Ct106;
  r_Ct(0,4)=cr_Ct106*cr_Ct107;
  r_Ct(0,5)=cr_Ct106*cr_Ct108;
  r_Ct(1,0)=-cr_Ct70*(cr_Ct109*cr_Ct111 + cr_Ct73);
  r_Ct(1,1)=cr_Ct70*(cr_Ct115*cr_Ct116*cr_Ct69 + cr_Ct53);
  r_Ct(1,2)=-cr_Ct70*(cr_Ct111*cr_Ct96 + cr_Ct73);
  r_Ct(1,3)=cr_Ct104*cr_Ct117;
  r_Ct(1,4)=cr_Ct107*cr_Ct117;
  r_Ct(1,5)=cr_Ct108*cr_Ct117;
  r_Ct(2,0)=-cr_Ct70*(cr_Ct109*cr_Ct118 + cr_Ct73);
  r_Ct(2,1)=-cr_Ct70*(cr_Ct118*cr_Ct89 + cr_Ct73);
  r_Ct(2,2)=cr_Ct70*(cr_Ct119*cr_Ct120*cr_Ct69 + cr_Ct53);
  r_Ct(2,3)=cr_Ct104*cr_Ct121;
  r_Ct(2,4)=cr_Ct107*cr_Ct121;
  r_Ct(2,5)=cr_Ct108*cr_Ct121;
  r_Ct(3,0)=-cr_Ct124*cr_Ct63;
  r_Ct(3,1)=-cr_Ct115*cr_Ct124;
  r_Ct(3,2)=-cr_Ct119*cr_Ct124;
  r_Ct(3,3)=cr_Ct7*(-cr_Ct125*cr_Ct43 + cr_Ct126);
  r_Ct(3,4)=cr_Ct128;
  r_Ct(3,5)=cr_Ct129;
  r_Ct(4,0)=-cr_Ct131*cr_Ct63;
  r_Ct(4,1)=-cr_Ct115*cr_Ct131;
  r_Ct(4,2)=-cr_Ct119*cr_Ct131;
  r_Ct(4,3)=cr_Ct128;
  r_Ct(4,4)=cr_Ct7*(-cr_Ct125*cr_Ct44 + cr_Ct126);
  r_Ct(4,5)=cr_Ct132;
  r_Ct(5,0)=-cr_Ct133*cr_Ct63;
  r_Ct(5,1)=-cr_Ct115*cr_Ct133;
  r_Ct(5,2)=-cr_Ct119*cr_Ct133;
  r_Ct(5,3)=cr_Ct129;
  r_Ct(5,4)=cr_Ct132;
  r_Ct(5,5)=cr_Ct7*(-cr_Ct125*cr_Ct45 + cr_Ct126);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<RankineYieldSurface<RankinePlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  const auto &r_props = rValues.GetMaterialProperties();
  const double Young = r_props[YOUNG_MODULUS];
  const double nu = r_props[POISSON_RATIO];
  const double Gf = r_props[FRACTURE_ENERGY];
  const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
  const double threshold = r_props[YIELD_STRESS];
  auto &r_Ct = rValues.GetConstitutiveMatrix();
  const auto &r_strain = rValues.GetStrainVector();

  const double cr_Ct0 = nu - 1.0;
  const double cr_Ct1 = 1.0/(1.0 - 0.5*characteristic_length*std::pow(threshold, 2)/(Gf*Young));
  const double cr_Ct2 = nu - 0.5;
  const double cr_Ct3 = std::pow(cr_Ct2, -2);
  const double cr_Ct4 = nu*r_strain[1];
  const double cr_Ct5 = -cr_Ct4;
  const double cr_Ct6 = cr_Ct0*r_strain[0];
  const double cr_Ct7 = nu*r_strain[2];
  const double cr_Ct8 = -cr_Ct7;
  const double cr_Ct9 = cr_Ct6 + cr_Ct8;
  const double cr_Ct10 = cr_Ct5 + cr_Ct9;
  const double cr_Ct11 = nu*r_strain[0];
  const double cr_Ct12 = -cr_Ct11;
  const double cr_Ct13 = cr_Ct0*r_strain[1];
  const double cr_Ct14 = cr_Ct12 + cr_Ct13 + cr_Ct8;
  const double cr_Ct15 = cr_Ct0*r_strain[2];
  const double cr_Ct16 = cr_Ct12 + cr_Ct15;
  const double cr_Ct17 = cr_Ct16 + cr_Ct5;
  const double cr_Ct18 = std::pow(r_strain[3], 2);
  const double cr_Ct19 = std::pow(r_strain[4], 2);
  const double cr_Ct20 = std::pow(r_strain[5], 2);
  const double cr_Ct21 = cr_Ct18 + cr_Ct19 + cr_Ct20;
  const double cr_Ct22 = std::pow(cr_Ct10, 2)*cr_Ct3 + std::pow(cr_Ct14, 2)*cr_Ct3 + std::pow(cr_Ct17, 2)*cr_Ct3 + cr_Ct21;
  const double cr_Ct23 = nu + 1.0;
  const double cr_Ct24 = std::pow(Young, 2)/std::pow(cr_Ct23, 2);
  const double cr_Ct25 = cr_Ct22*cr_Ct24;
  const double cr_Ct26 = std::pow(cr_Ct25, -1.0/2.0);
  const double cr_Ct27 = 2.0*nu;
  const double cr_Ct28 = 1.0/(cr_Ct27 - 1.0);
  const double cr_Ct29 = Young/cr_Ct23;
  const double cr_Ct30 = cr_Ct28*cr_Ct29;
  const double cr_Ct31 = cr_Ct26*cr_Ct30;
  const double cr_Ct32 = 0.66666666666666663*cr_Ct31;
  const double cr_Ct33 = 1.0/cr_Ct22;
  const double cr_Ct34 = cr_Ct33;
  const double cr_Ct35 = cr_Ct18*cr_Ct34;
  const double cr_Ct36 = cr_Ct19;
  const double cr_Ct37 = cr_Ct20*cr_Ct34;
  const double cr_Ct38 = cr_Ct17*cr_Ct3;
  const double cr_Ct39 = cr_Ct10*cr_Ct38;
  const double cr_Ct40 = 0.33333333333333331*cr_Ct3;
  const double cr_Ct41 = 2*cr_Ct11;
  const double cr_Ct42 = 2*cr_Ct4;
  const double cr_Ct43 = -cr_Ct42;
  const double cr_Ct44 = 2*cr_Ct7;
  const double cr_Ct45 = cr_Ct13 + cr_Ct15 - cr_Ct41 + cr_Ct43 - cr_Ct44 + cr_Ct6;
  const double cr_Ct46 = cr_Ct33*std::pow(cr_Ct45, 2);
  const double cr_Ct47 = cr_Ct3*(cr_Ct16 + cr_Ct43 + cr_Ct9);
  const double cr_Ct48 = cr_Ct14*cr_Ct47;
  const double cr_Ct49 = cr_Ct33*cr_Ct36 - cr_Ct34*cr_Ct39 - cr_Ct34*cr_Ct48 + cr_Ct35 + cr_Ct37 + cr_Ct40*cr_Ct46;
  const double cr_Ct50 = std::sqrt(cr_Ct49);
  const double cr_Ct51 = -nu;
  const double cr_Ct52 = cr_Ct51 + 0.5;
  const double cr_Ct53 = std::pow(cr_Ct52, -2);
  const double cr_Ct54 = cr_Ct51 + 1.0;
  const double cr_Ct55 = cr_Ct54*r_strain[2];
  const double cr_Ct56 = cr_Ct11 + cr_Ct4;
  const double cr_Ct57 = cr_Ct55 + cr_Ct56;
  const double cr_Ct58 = cr_Ct54*r_strain[1];
  const double cr_Ct59 = cr_Ct11 + cr_Ct7;
  const double cr_Ct60 = cr_Ct58 + cr_Ct59;
  const double cr_Ct61 = cr_Ct54*r_strain[0];
  const double cr_Ct62 = cr_Ct4 + cr_Ct7;
  const double cr_Ct63 = cr_Ct61 + cr_Ct62;
  const double cr_Ct64 = cr_Ct21 + cr_Ct53*std::pow(cr_Ct57, 2) + cr_Ct53*std::pow(cr_Ct60, 2) + cr_Ct53*std::pow(cr_Ct63, 2);
  const double cr_Ct65 = 1.0/cr_Ct64;
  const double cr_Ct66 = cr_Ct65;
  const double cr_Ct67 = cr_Ct18*cr_Ct66;
  const double cr_Ct68 = cr_Ct36*cr_Ct65;
  const double cr_Ct69 = cr_Ct20*cr_Ct66;
  const double cr_Ct70 = cr_Ct53*cr_Ct66;
  const double cr_Ct71 = cr_Ct57*cr_Ct63;
  const double cr_Ct72 = cr_Ct70*cr_Ct71;
  const double cr_Ct73 = cr_Ct55 + cr_Ct61;
  const double cr_Ct74 = cr_Ct41 + cr_Ct42 + cr_Ct44;
  const double cr_Ct75 = cr_Ct58 + cr_Ct73 + cr_Ct74;
  const double cr_Ct76 = std::pow(cr_Ct75, 2);
  const double cr_Ct77 = cr_Ct53*cr_Ct65;
  const double cr_Ct78 = cr_Ct76*cr_Ct77;
  const double cr_Ct79 = 0.33333333333333331*cr_Ct78;
  const double cr_Ct80 = cr_Ct42 + cr_Ct59;
  const double cr_Ct81 = cr_Ct73 + cr_Ct80;
  const double cr_Ct82 = cr_Ct60*cr_Ct81;
  const double cr_Ct83 = cr_Ct70*cr_Ct82;
  const double cr_Ct84 = cr_Ct67 + cr_Ct68 + cr_Ct69 - cr_Ct72 + cr_Ct79 - cr_Ct83;
  const double cr_Ct85 = std::pow(cr_Ct84, 3);
  const double cr_Ct86 = std::pow(cr_Ct85, -1.0/2.0);
  const double cr_Ct87 = 5.196152422706632*cr_Ct29;
  const double cr_Ct88 = r_strain[4]*r_strain[5];
  const double cr_Ct89 = cr_Ct88*r_strain[3];
  const double cr_Ct90 = cr_Ct17*cr_Ct28;
  const double cr_Ct91 = cr_Ct14*cr_Ct28;
  const double cr_Ct92 = cr_Ct53;
  const double cr_Ct93 = cr_Ct57*cr_Ct60;
  const double cr_Ct94 = cr_Ct92*cr_Ct93;
  const double cr_Ct95 = cr_Ct10*cr_Ct28;
  const double cr_Ct96 = cr_Ct95*(cr_Ct36 - cr_Ct94);
  const double cr_Ct97 = 0.037037037037037035/std::pow(cr_Ct2, 3);
  const double cr_Ct98 = 0.037037037037037035*cr_Ct28;
  const double cr_Ct99 = 9.0*cr_Ct65;
  const double cr_Ct100 = cr_Ct18*cr_Ct99;
  const double cr_Ct101 = cr_Ct19*cr_Ct99;
  const double cr_Ct102 = cr_Ct20*cr_Ct99;
  const double cr_Ct103 = cr_Ct53*cr_Ct99;
  const double cr_Ct104 = cr_Ct103*cr_Ct71;
  const double cr_Ct105 = cr_Ct103*cr_Ct82;
  const double cr_Ct106 = cr_Ct45*(cr_Ct100 + cr_Ct101 + cr_Ct102 - cr_Ct104 - cr_Ct105);
  const double cr_Ct107 = cr_Ct106*cr_Ct98 + cr_Ct33*std::pow(cr_Ct45, 3)*cr_Ct97;
  const double cr_Ct108 = cr_Ct26*(cr_Ct107 + cr_Ct34*cr_Ct89 - cr_Ct35*cr_Ct90 - cr_Ct37*cr_Ct91 - cr_Ct66*cr_Ct96);
  const double cr_Ct109 = 0.33333333333333331*std::acos(cr_Ct108*cr_Ct86*cr_Ct87);
  const double cr_Ct110 = std::cos(cr_Ct109);
  const double cr_Ct111 = cr_Ct110*cr_Ct50;
  const double cr_Ct112 = 1.0/(cr_Ct10*cr_Ct32 + 1.1547005383792515*cr_Ct111 + cr_Ct14*cr_Ct32 + cr_Ct17*cr_Ct32);
  const double cr_Ct113 = 2.0*threshold;
  const double cr_Ct114 = cr_Ct0*(cr_Ct1*(cr_Ct112*cr_Ct113*cr_Ct26 - 1.0) + 1.0);
  const double cr_Ct115 = -cr_Ct6;
  const double cr_Ct116 = cr_Ct115 + cr_Ct62;
  const double cr_Ct117 = 2.0*cr_Ct33;
  const double cr_Ct118 = cr_Ct112*cr_Ct117;
  const double cr_Ct119 = cr_Ct57*nu;
  const double cr_Ct120 = cr_Ct60*nu;
  const double cr_Ct121 = cr_Ct54*cr_Ct63;
  const double cr_Ct122 = cr_Ct119 + cr_Ct120 + cr_Ct121;
  const double cr_Ct123 = cr_Ct122*cr_Ct53;
  const double cr_Ct124 = 0.57735026918962573*cr_Ct31;
  const double cr_Ct125 = std::pow(cr_Ct10*cr_Ct124 + cr_Ct111 + cr_Ct124*cr_Ct14 + cr_Ct124*cr_Ct17, -2);
  const double cr_Ct126 = cr_Ct63*nu;
  const double cr_Ct127 = 0.5*cr_Ct77;
  const double cr_Ct128 = -cr_Ct126*cr_Ct127;
  const double cr_Ct129 = cr_Ct54*cr_Ct57;
  const double cr_Ct130 = cr_Ct81*nu;
  const double cr_Ct131 = -cr_Ct127*cr_Ct130;
  const double cr_Ct132 = cr_Ct122*cr_Ct18;
  const double cr_Ct133 = std::pow(cr_Ct64, -2);
  const double cr_Ct134 = cr_Ct133*cr_Ct92;
  const double cr_Ct135 = cr_Ct133*cr_Ct53;
  const double cr_Ct136 = cr_Ct122*cr_Ct135;
  const double cr_Ct137 = cr_Ct134*cr_Ct20;
  const double cr_Ct138 = std::pow(cr_Ct52, -4);
  const double cr_Ct139 = cr_Ct122*cr_Ct138;
  const double cr_Ct140 = cr_Ct133*cr_Ct139;
  const double cr_Ct141 = cr_Ct140*cr_Ct71;
  const double cr_Ct142 = cr_Ct140*cr_Ct82;
  const double cr_Ct143 = cr_Ct122*cr_Ct77;
  const double cr_Ct144 = cr_Ct143*cr_Ct57;
  const double cr_Ct145 = cr_Ct143*cr_Ct63;
  const double cr_Ct146 = cr_Ct60*(-cr_Ct144 - cr_Ct145 + 1.0);
  const double cr_Ct147 = 2.0*cr_Ct65;
  const double cr_Ct148 = cr_Ct123*cr_Ct147;
  const double cr_Ct149 = cr_Ct27 + 2.0;
  const double cr_Ct150 = -cr_Ct148*cr_Ct57 - cr_Ct148*cr_Ct60 - cr_Ct148*cr_Ct63 + cr_Ct149;
  const double cr_Ct151 = cr_Ct150*cr_Ct75;
  const double cr_Ct152 = 0.16666666666666666*cr_Ct77;
  const double cr_Ct153 = 1.0/(1.0 - cr_Ct27);
  const double cr_Ct154 = -cr_Ct36 + cr_Ct94;
  const double cr_Ct155 = cr_Ct154*cr_Ct63;
  const double cr_Ct156 = cr_Ct153*cr_Ct155;
  const double cr_Ct157 = cr_Ct153*cr_Ct57;
  const double cr_Ct158 = cr_Ct153*cr_Ct60;
  const double cr_Ct159 = cr_Ct153*(-cr_Ct100 - cr_Ct101 - cr_Ct102 + cr_Ct104 + cr_Ct105);
  const double cr_Ct160 = -0.037037037037037035*cr_Ct159*cr_Ct75 + 0.037037037037037035*cr_Ct65*std::pow(cr_Ct75, 3)/std::pow(cr_Ct52, 3);
  const double cr_Ct161 = cr_Ct156*cr_Ct66 - cr_Ct157*cr_Ct67 - cr_Ct158*cr_Ct69 + cr_Ct160 + cr_Ct66*cr_Ct89;
  const double cr_Ct162 = cr_Ct24*cr_Ct64;
  const double cr_Ct163 = std::pow(cr_Ct162, -1.0/2.0);
  const double cr_Ct164 = -cr_Ct67 - cr_Ct68 - cr_Ct69 + cr_Ct72 - cr_Ct79 + cr_Ct83;
  const double cr_Ct165 = std::pow(cr_Ct164, 3);
  const double cr_Ct166 = cr_Ct163/std::sqrt(-cr_Ct165);
  const double cr_Ct167 = 0.33333333333333331*std::acos(cr_Ct161*cr_Ct166*cr_Ct87);
  const double cr_Ct168 = std::cos(cr_Ct167);
  const double cr_Ct169 = std::sqrt(cr_Ct84);
  const double cr_Ct170 = 1.7320508075688774*cr_Ct168/cr_Ct169;
  const double cr_Ct171 = cr_Ct170*(-cr_Ct122*cr_Ct137 - cr_Ct127*cr_Ct129 - cr_Ct127*cr_Ct146 + cr_Ct128 + cr_Ct131 - cr_Ct132*cr_Ct134 - cr_Ct136*cr_Ct36 + cr_Ct141 + 0.5*cr_Ct142 + cr_Ct151*cr_Ct152);
  const double cr_Ct172 = std::pow(cr_Ct162, -3.0/2.0);
  const double cr_Ct173 = std::pow(Young, 3)/std::pow(cr_Ct23, 3);
  const double cr_Ct174 = cr_Ct172*cr_Ct173*cr_Ct92;
  const double cr_Ct175 = cr_Ct122*cr_Ct174;
  const double cr_Ct176 = cr_Ct157*cr_Ct175;
  const double cr_Ct177 = cr_Ct158*cr_Ct175;
  const double cr_Ct178 = cr_Ct153*cr_Ct63;
  const double cr_Ct179 = cr_Ct175*cr_Ct178;
  const double cr_Ct180 = cr_Ct153*cr_Ct65;
  const double cr_Ct181 = cr_Ct18*cr_Ct180;
  const double cr_Ct182 = 5.196152422706632*nu;
  const double cr_Ct183 = -cr_Ct181*cr_Ct182;
  const double cr_Ct184 = cr_Ct180*cr_Ct20;
  const double cr_Ct185 = -cr_Ct182*cr_Ct184;
  const double cr_Ct186 = 5.196152422706632*cr_Ct180;
  const double cr_Ct187 = cr_Ct154*cr_Ct186;
  const double cr_Ct188 = 15.588457268119896*cr_Ct136;
  const double cr_Ct189 = 15.588457268119896*cr_Ct135;
  const double cr_Ct190 = cr_Ct157*cr_Ct189;
  const double cr_Ct191 = cr_Ct158*cr_Ct20;
  const double cr_Ct192 = 5.196152422706632*cr_Ct156;
  const double cr_Ct193 = cr_Ct119*cr_Ct92;
  const double cr_Ct194 = cr_Ct120*cr_Ct92;
  const double cr_Ct195 = 2.0*cr_Ct19;
  const double cr_Ct196 = cr_Ct195*cr_Ct65;
  const double cr_Ct197 = cr_Ct147*cr_Ct93;
  const double cr_Ct198 = cr_Ct186*cr_Ct63;
  const double cr_Ct199 = 0.57735026918962573*cr_Ct153*cr_Ct78;
  const double cr_Ct200 = 0.037037037037037035*cr_Ct60;
  const double cr_Ct201 = -0.037037037037037035*nu - 0.037037037037037035;
  const double cr_Ct202 = 5.196152422706632*cr_Ct159;
  const double cr_Ct203 = cr_Ct103*cr_Ct126;
  const double cr_Ct204 = cr_Ct103*cr_Ct130;
  const double cr_Ct205 = 18.0*cr_Ct133;
  const double cr_Ct206 = cr_Ct18*cr_Ct205;
  const double cr_Ct207 = cr_Ct19*cr_Ct205;
  const double cr_Ct208 = cr_Ct20*cr_Ct205;
  const double cr_Ct209 = cr_Ct205*cr_Ct71;
  const double cr_Ct210 = cr_Ct153*cr_Ct75;
  const double cr_Ct211 = 0.19245008972987526*cr_Ct210;
  const double cr_Ct212 = 3.0*cr_Ct77;
  const double cr_Ct213 = cr_Ct126*cr_Ct212;
  const double cr_Ct214 = cr_Ct130*cr_Ct212;
  const double cr_Ct215 = 6.0*cr_Ct135;
  const double cr_Ct216 = 6.0*cr_Ct136;
  const double cr_Ct217 = 2.598076211353316*cr_Ct161/cr_Ct164;
  const double cr_Ct218 = cr_Ct65*r_strain[3];
  const double cr_Ct219 = 0.1111111111111111*cr_Ct29;
  const double cr_Ct220 = cr_Ct166*cr_Ct169*cr_Ct219*std::sin(cr_Ct167)/std::sqrt(0.037037037037037035 + cr_Ct65*std::pow(cr_Ct155*cr_Ct180 + cr_Ct160 - cr_Ct181*cr_Ct57 - cr_Ct184*cr_Ct60 + cr_Ct218*cr_Ct88, 2)/cr_Ct165);
  const double cr_Ct221 = cr_Ct220*(cr_Ct132*cr_Ct190 - cr_Ct136*cr_Ct192 + cr_Ct150*cr_Ct199 + cr_Ct183 + cr_Ct185 + cr_Ct187*cr_Ct54 + cr_Ct188*cr_Ct191 - cr_Ct188*cr_Ct89 + cr_Ct198*(cr_Ct123*cr_Ct196 - cr_Ct139*cr_Ct197 + cr_Ct193 + cr_Ct194) + cr_Ct202*(cr_Ct143*cr_Ct200 + 0.037037037037037035*cr_Ct144 + 0.037037037037037035*cr_Ct145 + cr_Ct201) - cr_Ct211*(cr_Ct103*cr_Ct129 + cr_Ct103*cr_Ct146 + cr_Ct123*cr_Ct206 + cr_Ct123*cr_Ct207 + cr_Ct123*cr_Ct208 - cr_Ct139*cr_Ct209 - 9.0*cr_Ct142 + cr_Ct203 + cr_Ct204) - cr_Ct217*(cr_Ct129*cr_Ct212 + cr_Ct132*cr_Ct215 - 6.0*cr_Ct141 - 3.0*cr_Ct142 + cr_Ct146*cr_Ct212 - cr_Ct151*cr_Ct70 + cr_Ct19*cr_Ct216 + cr_Ct20*cr_Ct216 + cr_Ct213 + cr_Ct214));
  const double cr_Ct222 = cr_Ct163*cr_Ct29;
  const double cr_Ct223 = cr_Ct153*cr_Ct222;
  const double cr_Ct224 = cr_Ct223*cr_Ct54;
  const double cr_Ct225 = cr_Ct223*cr_Ct27;
  const double cr_Ct226 = cr_Ct224 + cr_Ct225;
  const double cr_Ct227 = cr_Ct1*threshold;
  const double cr_Ct228 = cr_Ct227*cr_Ct26;
  const double cr_Ct229 = cr_Ct228*(cr_Ct118*cr_Ct123 + cr_Ct125*(cr_Ct171 - cr_Ct176 - cr_Ct177 - cr_Ct179 + cr_Ct221 + cr_Ct226));
  const double cr_Ct230 = 0.66666666666666663*cr_Ct222;
  const double cr_Ct231 = cr_Ct223*cr_Ct63;
  const double cr_Ct232 = cr_Ct168*cr_Ct169;
  const double cr_Ct233 = 1.0/(cr_Ct157*cr_Ct230 + cr_Ct158*cr_Ct230 + 0.66666666666666663*cr_Ct231 + 1.1547005383792515*cr_Ct232);
  const double cr_Ct234 = cr_Ct1*(cr_Ct113*cr_Ct163*cr_Ct233 - 1.0) + 1.0;
  const double cr_Ct235 = -cr_Ct234*nu;
  const double cr_Ct236 = cr_Ct54*cr_Ct60;
  const double cr_Ct237 = cr_Ct119 + cr_Ct126 + cr_Ct236;
  const double cr_Ct238 = cr_Ct147*cr_Ct53;
  const double cr_Ct239 = cr_Ct237*cr_Ct238;
  const double cr_Ct240 = 0.57735026918962573*cr_Ct222;
  const double cr_Ct241 = std::pow(cr_Ct157*cr_Ct240 + cr_Ct158*cr_Ct240 + 0.57735026918962573*cr_Ct231 + cr_Ct232, -2);
  const double cr_Ct242 = -cr_Ct119*cr_Ct127;
  const double cr_Ct243 = cr_Ct54*cr_Ct81;
  const double cr_Ct244 = cr_Ct18*cr_Ct237;
  const double cr_Ct245 = cr_Ct135*cr_Ct36;
  const double cr_Ct246 = cr_Ct138*cr_Ct237;
  const double cr_Ct247 = cr_Ct133*cr_Ct246;
  const double cr_Ct248 = cr_Ct247*cr_Ct71;
  const double cr_Ct249 = cr_Ct247*cr_Ct82;
  const double cr_Ct250 = cr_Ct237*cr_Ct57;
  const double cr_Ct251 = -cr_Ct238*cr_Ct250 - cr_Ct239*cr_Ct63;
  const double cr_Ct252 = cr_Ct60*cr_Ct77*(cr_Ct251 + 4.0*nu);
  const double cr_Ct253 = cr_Ct149 - cr_Ct239*cr_Ct60 + cr_Ct251;
  const double cr_Ct254 = cr_Ct253*cr_Ct75;
  const double cr_Ct255 = cr_Ct170*(-cr_Ct127*cr_Ct243 + cr_Ct128 - cr_Ct134*cr_Ct244 - cr_Ct137*cr_Ct237 + cr_Ct152*cr_Ct254 - cr_Ct237*cr_Ct245 + cr_Ct242 + cr_Ct248 + 0.5*cr_Ct249 - 0.25*cr_Ct252);
  const double cr_Ct256 = cr_Ct174*cr_Ct237;
  const double cr_Ct257 = cr_Ct157*cr_Ct256;
  const double cr_Ct258 = cr_Ct158*cr_Ct256;
  const double cr_Ct259 = cr_Ct237*cr_Ct63;
  const double cr_Ct260 = cr_Ct153*cr_Ct174*cr_Ct259;
  const double cr_Ct261 = 5.196152422706632*cr_Ct54;
  const double cr_Ct262 = cr_Ct187*nu;
  const double cr_Ct263 = cr_Ct189*cr_Ct237;
  const double cr_Ct264 = cr_Ct135*cr_Ct192;
  const double cr_Ct265 = cr_Ct237*cr_Ct53;
  const double cr_Ct266 = 0.037037037037037035*cr_Ct77;
  const double cr_Ct267 = cr_Ct103*cr_Ct119;
  const double cr_Ct268 = cr_Ct119*cr_Ct212;
  const double cr_Ct269 = cr_Ct215*cr_Ct237;
  const double cr_Ct270 = cr_Ct220*(cr_Ct183 - cr_Ct184*cr_Ct261 + cr_Ct190*cr_Ct244 + cr_Ct191*cr_Ct263 + cr_Ct198*(cr_Ct129*cr_Ct92 + cr_Ct194 + cr_Ct196*cr_Ct265 - cr_Ct197*cr_Ct246) + cr_Ct199*cr_Ct253 + cr_Ct202*(cr_Ct200*cr_Ct237*cr_Ct77 + cr_Ct201 + cr_Ct250*cr_Ct266 + cr_Ct259*cr_Ct266) - cr_Ct211*(cr_Ct103*cr_Ct243 + cr_Ct203 + cr_Ct206*cr_Ct265 + cr_Ct207*cr_Ct265 + cr_Ct208*cr_Ct265 - cr_Ct209*cr_Ct246 - 9.0*cr_Ct249 + 4.5*cr_Ct252 + cr_Ct267) - cr_Ct217*(cr_Ct19*cr_Ct269 + cr_Ct20*cr_Ct269 + cr_Ct212*cr_Ct243 + cr_Ct213 + cr_Ct215*cr_Ct244 - 6.0*cr_Ct248 - 3.0*cr_Ct249 + 1.5*cr_Ct252 - cr_Ct254*cr_Ct70 + cr_Ct268) - cr_Ct237*cr_Ct264 + cr_Ct262 - cr_Ct263*cr_Ct89);
  const double cr_Ct271 = -cr_Ct224 - cr_Ct225;
  const double cr_Ct272 = -cr_Ct233*cr_Ct239 + cr_Ct241*(-cr_Ct255 + cr_Ct257 + cr_Ct258 + cr_Ct260 - cr_Ct270 + cr_Ct271);
  const double cr_Ct273 = cr_Ct163*cr_Ct227;
  const double cr_Ct274 = cr_Ct10*cr_Ct273;
  const double cr_Ct275 = cr_Ct120 + cr_Ct126 + cr_Ct129;
  const double cr_Ct276 = cr_Ct238*cr_Ct275;
  const double cr_Ct277 = cr_Ct18*cr_Ct275;
  const double cr_Ct278 = cr_Ct138*cr_Ct275;
  const double cr_Ct279 = cr_Ct133*cr_Ct278;
  const double cr_Ct280 = cr_Ct279*cr_Ct71;
  const double cr_Ct281 = cr_Ct279*cr_Ct82;
  const double cr_Ct282 = cr_Ct275*cr_Ct77;
  const double cr_Ct283 = cr_Ct282*cr_Ct57;
  const double cr_Ct284 = cr_Ct282*cr_Ct63;
  const double cr_Ct285 = cr_Ct60*(-cr_Ct283 - cr_Ct284 + 1.0);
  const double cr_Ct286 = cr_Ct149 - cr_Ct276*cr_Ct57 - cr_Ct276*cr_Ct60 - cr_Ct276*cr_Ct63;
  const double cr_Ct287 = cr_Ct286*cr_Ct75;
  const double cr_Ct288 = cr_Ct170*(-cr_Ct121*cr_Ct127 - cr_Ct127*cr_Ct285 + cr_Ct131 - cr_Ct134*cr_Ct277 - cr_Ct137*cr_Ct275 + cr_Ct152*cr_Ct287 + cr_Ct242 - cr_Ct245*cr_Ct275 + cr_Ct280 + 0.5*cr_Ct281);
  const double cr_Ct289 = cr_Ct174*cr_Ct275;
  const double cr_Ct290 = cr_Ct157*cr_Ct289;
  const double cr_Ct291 = cr_Ct158*cr_Ct289;
  const double cr_Ct292 = cr_Ct178*cr_Ct289;
  const double cr_Ct293 = cr_Ct189*cr_Ct275;
  const double cr_Ct294 = cr_Ct275*cr_Ct53;
  const double cr_Ct295 = cr_Ct215*cr_Ct275;
  const double cr_Ct296 = cr_Ct220*(-cr_Ct181*cr_Ct261 + cr_Ct185 + cr_Ct190*cr_Ct277 + cr_Ct191*cr_Ct293 + cr_Ct198*(cr_Ct193 + cr_Ct196*cr_Ct294 - cr_Ct197*cr_Ct278 + cr_Ct236*cr_Ct92) + cr_Ct199*cr_Ct286 + cr_Ct202*(cr_Ct200*cr_Ct282 + cr_Ct201 + 0.037037037037037035*cr_Ct283 + 0.037037037037037035*cr_Ct284) - cr_Ct211*(cr_Ct103*cr_Ct121 + cr_Ct103*cr_Ct285 + cr_Ct204 + cr_Ct206*cr_Ct294 + cr_Ct207*cr_Ct294 + cr_Ct208*cr_Ct294 - cr_Ct209*cr_Ct278 + cr_Ct267 - 9.0*cr_Ct281) - cr_Ct217*(cr_Ct121*cr_Ct212 + cr_Ct19*cr_Ct295 + cr_Ct20*cr_Ct295 + cr_Ct212*cr_Ct285 + cr_Ct214 + cr_Ct215*cr_Ct277 + cr_Ct268 - 6.0*cr_Ct280 - 3.0*cr_Ct281 - cr_Ct287*cr_Ct70) + cr_Ct262 - cr_Ct264*cr_Ct275 - cr_Ct293*cr_Ct89);
  const double cr_Ct297 = -cr_Ct233*cr_Ct276 + cr_Ct241*(cr_Ct271 - cr_Ct288 + cr_Ct290 + cr_Ct291 + cr_Ct292 - cr_Ct296);
  const double cr_Ct298 = -cr_Ct15;
  const double cr_Ct299 = cr_Ct298 + cr_Ct56;
  const double cr_Ct300 = -cr_Ct13;
  const double cr_Ct301 = cr_Ct300 + cr_Ct59;
  const double cr_Ct302 = std::pow(cr_Ct116, 2)*cr_Ct3 + cr_Ct21 + std::pow(cr_Ct299, 2)*cr_Ct3 + cr_Ct3*std::pow(cr_Ct301, 2);
  const double cr_Ct303 = 1.0/cr_Ct302;
  const double cr_Ct304 = cr_Ct303*r_strain[3];
  const double cr_Ct305 = cr_Ct24*cr_Ct302;
  const double cr_Ct306 = std::pow(cr_Ct305, -1.0/2.0);
  const double cr_Ct307 = cr_Ct30*cr_Ct306;
  const double cr_Ct308 = 0.66666666666666663*cr_Ct307;
  const double cr_Ct309 = cr_Ct303;
  const double cr_Ct310 = cr_Ct18*cr_Ct309;
  const double cr_Ct311 = cr_Ct20*cr_Ct309;
  const double cr_Ct312 = cr_Ct299*cr_Ct3;
  const double cr_Ct313 = cr_Ct116*cr_Ct312;
  const double cr_Ct314 = cr_Ct115 + cr_Ct298;
  const double cr_Ct315 = cr_Ct300 + cr_Ct314 + cr_Ct74;
  const double cr_Ct316 = cr_Ct303*std::pow(cr_Ct315, 2);
  const double cr_Ct317 = cr_Ct3*cr_Ct301*(cr_Ct314 + cr_Ct80);
  const double cr_Ct318 = std::sqrt(cr_Ct303*cr_Ct36 - cr_Ct309*cr_Ct313 - cr_Ct309*cr_Ct317 + cr_Ct310 + cr_Ct311 + cr_Ct316*cr_Ct40);
  const double cr_Ct319 = std::pow(cr_Ct49, 3);
  const double cr_Ct320 = std::pow(cr_Ct319, -1.0/2.0);
  const double cr_Ct321 = cr_Ct28*cr_Ct299;
  const double cr_Ct322 = cr_Ct28*cr_Ct301;
  const double cr_Ct323 = cr_Ct14;
  const double cr_Ct324 = cr_Ct116*(-cr_Ct323*cr_Ct38 + cr_Ct36);
  const double cr_Ct325 = 9.0*cr_Ct33;
  const double cr_Ct326 = cr_Ct315*(cr_Ct18*cr_Ct325 + cr_Ct19*cr_Ct325 + cr_Ct20*cr_Ct325 - cr_Ct325*cr_Ct39 - cr_Ct325*cr_Ct48);
  const double cr_Ct327 = -cr_Ct303*std::pow(cr_Ct315, 3)*cr_Ct97 - cr_Ct326*cr_Ct98;
  const double cr_Ct328 = cr_Ct306*(cr_Ct28*cr_Ct324*cr_Ct34 + cr_Ct309*cr_Ct89 + cr_Ct310*cr_Ct321 + cr_Ct311*cr_Ct322 + cr_Ct327);
  const double cr_Ct329 = 0.33333333333333331*std::acos(cr_Ct320*cr_Ct328*cr_Ct87);
  const double cr_Ct330 = std::cos(cr_Ct329);
  const double cr_Ct331 = cr_Ct318*cr_Ct330;
  const double cr_Ct332 = 2.0/(-cr_Ct116*cr_Ct308 - cr_Ct299*cr_Ct308 - cr_Ct301*cr_Ct308 + 1.1547005383792515*cr_Ct331);
  const double cr_Ct333 = 0.57735026918962573*cr_Ct307;
  const double cr_Ct334 = std::pow(-cr_Ct116*cr_Ct333 - cr_Ct299*cr_Ct333 - cr_Ct301*cr_Ct333 + cr_Ct331, -2);
  const double cr_Ct335 = std::pow(cr_Ct305, -3.0/2.0);
  const double cr_Ct336 = cr_Ct321*cr_Ct335*r_strain[3];
  const double cr_Ct337 = cr_Ct173;
  const double cr_Ct338 = cr_Ct335*cr_Ct337;
  const double cr_Ct339 = cr_Ct322*cr_Ct338;
  const double cr_Ct340 = cr_Ct28*r_strain[3];
  const double cr_Ct341 = cr_Ct116*cr_Ct340;
  const double cr_Ct342 = 0.86602540378443871*r_strain[3];
  const double cr_Ct343 = std::pow(cr_Ct22, -2);
  const double cr_Ct344 = 2.0*cr_Ct343;
  const double cr_Ct345 = std::pow(Young, 4)/std::pow(cr_Ct23, 4);
  const double cr_Ct346 = cr_Ct172*cr_Ct345;
  const double cr_Ct347 = cr_Ct153*cr_Ct26*cr_Ct346*cr_Ct81*cr_Ct91;
  const double cr_Ct348 = cr_Ct26*cr_Ct45;
  const double cr_Ct349 = cr_Ct210*cr_Ct28*cr_Ct346*cr_Ct348;
  const double cr_Ct350 = cr_Ct330*(-cr_Ct117 + cr_Ct18*cr_Ct344 + cr_Ct195*cr_Ct343 + cr_Ct20*cr_Ct344 - cr_Ct323*cr_Ct343*cr_Ct47 - cr_Ct344*cr_Ct39 - 4.0*cr_Ct347 + 2.6666666666666665*cr_Ct349)/cr_Ct318;
  const double cr_Ct351 = 5.196152422706632*cr_Ct24;
  const double cr_Ct352 = cr_Ct335*cr_Ct351;
  const double cr_Ct353 = 15.588457268119896*cr_Ct345;
  const double cr_Ct354 = cr_Ct353/std::pow(cr_Ct305, 5.0/2.0);
  const double cr_Ct355 = cr_Ct18*cr_Ct88;
  const double cr_Ct356 = 10.392304845413264*cr_Ct24;
  const double cr_Ct357 = std::pow(r_strain[3], 3);
  const double cr_Ct358 = cr_Ct20*cr_Ct301;
  const double cr_Ct359 = cr_Ct354*cr_Ct358;
  const double cr_Ct360 = std::pow(cr_Ct302, -2);
  const double cr_Ct361 = 5.196152422706632*cr_Ct306;
  const double cr_Ct362 = cr_Ct360*cr_Ct361*(cr_Ct195 - 2.0*cr_Ct301*cr_Ct312);
  const double cr_Ct363 = cr_Ct28*cr_Ct33;
  const double cr_Ct364 = cr_Ct324*cr_Ct363;
  const double cr_Ct365 = cr_Ct352*r_strain[3];
  const double cr_Ct366 = std::pow(cr_Ct25, -3.0/2.0);
  const double cr_Ct367 = cr_Ct340*cr_Ct366;
  const double cr_Ct368 = cr_Ct3*cr_Ct45;
  const double cr_Ct369 = 1.1547005383792515*cr_Ct24;
  const double cr_Ct370 = cr_Ct316*cr_Ct368*cr_Ct369;
  const double cr_Ct371 = 0.19245008972987526*cr_Ct340;
  const double cr_Ct372 = cr_Ct24*cr_Ct371;
  const double cr_Ct373 = cr_Ct326*cr_Ct335;
  const double cr_Ct374 = 18.0*cr_Ct343;
  const double cr_Ct375 = cr_Ct343*cr_Ct48;
  const double cr_Ct376 = cr_Ct348*(-cr_Ct18*cr_Ct374 - cr_Ct19*cr_Ct374 - cr_Ct20*cr_Ct374 + 18.0*cr_Ct33 + 36.0*cr_Ct347 + cr_Ct374*cr_Ct39 + 9.0*cr_Ct375);
  const double cr_Ct377 = 2.598076211353316*r_strain[3];
  const double cr_Ct378 = 6.0*cr_Ct360;
  const double cr_Ct379 = cr_Ct306*cr_Ct345*cr_Ct366;
  const double cr_Ct380 = cr_Ct328*(cr_Ct18*cr_Ct378 + cr_Ct19*cr_Ct378 + cr_Ct20*cr_Ct378 + 3.0*cr_Ct301*cr_Ct379*cr_Ct47 - 6.0*cr_Ct303 - cr_Ct313*cr_Ct378 - 2.0*cr_Ct315*cr_Ct368*cr_Ct379 - 3.0*cr_Ct317*cr_Ct360)/cr_Ct49;
  const double cr_Ct381 = cr_Ct28*cr_Ct303;
  const double cr_Ct382 = cr_Ct18*cr_Ct299;
  const double cr_Ct383 = cr_Ct219*cr_Ct318*cr_Ct320*std::sin(cr_Ct329)/std::sqrt(-cr_Ct303*std::pow(cr_Ct304*cr_Ct88 + cr_Ct327 + cr_Ct358*cr_Ct381 + cr_Ct364 + cr_Ct381*cr_Ct382, 2)/cr_Ct319 + 0.037037037037037035);
  const double cr_Ct384 = cr_Ct304*cr_Ct332 - cr_Ct334*(-cr_Ct336*cr_Ct337 - cr_Ct338*cr_Ct341 - cr_Ct339*r_strain[3] + cr_Ct342*cr_Ct350 - cr_Ct383*(-cr_Ct321*cr_Ct354*cr_Ct357 + cr_Ct336*cr_Ct356 - cr_Ct340*cr_Ct359 - cr_Ct341*cr_Ct362 + cr_Ct352*cr_Ct88 - cr_Ct354*cr_Ct355 - cr_Ct364*cr_Ct365 - cr_Ct367*cr_Ct370 + cr_Ct371*cr_Ct376 + cr_Ct372*cr_Ct373 + cr_Ct377*cr_Ct380));
  const double cr_Ct385 = cr_Ct227*cr_Ct307;
  const double cr_Ct386 = cr_Ct116*cr_Ct385;
  const double cr_Ct387 = cr_Ct303*cr_Ct332;
  const double cr_Ct388 = cr_Ct321*cr_Ct338;
  const double cr_Ct389 = cr_Ct28*r_strain[4];
  const double cr_Ct390 = cr_Ct116*cr_Ct338;
  const double cr_Ct391 = 0.86602540378443871*cr_Ct350;
  const double cr_Ct392 = cr_Ct354*r_strain[3];
  const double cr_Ct393 = cr_Ct19*r_strain[5];
  const double cr_Ct394 = cr_Ct354*cr_Ct382;
  const double cr_Ct395 = cr_Ct352*cr_Ct364;
  const double cr_Ct396 = cr_Ct366*cr_Ct370;
  const double cr_Ct397 = cr_Ct14*cr_Ct38;
  const double cr_Ct398 = 0.19245008972987526*cr_Ct24;
  const double cr_Ct399 = cr_Ct373*cr_Ct398;
  const double cr_Ct400 = 0.19245008972987526*cr_Ct376;
  const double cr_Ct401 = 2.598076211353316*cr_Ct380;
  const double cr_Ct402 = -cr_Ct334*(-cr_Ct339*r_strain[4] - cr_Ct383*(cr_Ct116*cr_Ct361*cr_Ct363*r_strain[4]*(-cr_Ct117*cr_Ct19 + cr_Ct117*cr_Ct397 + 2.0) - cr_Ct359*cr_Ct389 + cr_Ct365*r_strain[5] - cr_Ct389*cr_Ct394 - cr_Ct389*cr_Ct396 + cr_Ct389*cr_Ct399 + cr_Ct389*cr_Ct400 - cr_Ct392*cr_Ct393 - cr_Ct395*r_strain[4] + cr_Ct401*r_strain[4]) - cr_Ct388*r_strain[4] - cr_Ct389*cr_Ct390 + cr_Ct391*r_strain[4]) + cr_Ct387*r_strain[4];
  const double cr_Ct403 = cr_Ct28*r_strain[5];
  const double cr_Ct404 = cr_Ct20*r_strain[4];
  const double cr_Ct405 = cr_Ct356*r_strain[5];
  const double cr_Ct406 = std::pow(r_strain[5], 3);
  const double cr_Ct407 = -cr_Ct334*(-cr_Ct339*r_strain[5] - cr_Ct383*(-cr_Ct116*cr_Ct362*cr_Ct403 + cr_Ct322*cr_Ct335*cr_Ct405 - cr_Ct322*cr_Ct354*cr_Ct406 + cr_Ct365*r_strain[4] - cr_Ct392*cr_Ct404 - cr_Ct394*cr_Ct403 - cr_Ct395*r_strain[5] - cr_Ct396*cr_Ct403 + cr_Ct399*cr_Ct403 + cr_Ct400*cr_Ct403 + cr_Ct401*r_strain[5]) - cr_Ct388*r_strain[5] - cr_Ct390*cr_Ct403 + cr_Ct391*r_strain[5]) + cr_Ct387*r_strain[5];
  const double cr_Ct408 = -cr_Ct148*cr_Ct233 + cr_Ct241*(-cr_Ct171 + cr_Ct176 + cr_Ct177 + cr_Ct179 - cr_Ct221 + cr_Ct271);
  const double cr_Ct409 = cr_Ct14*cr_Ct273;
  const double cr_Ct410 = cr_Ct228*(cr_Ct118*cr_Ct265 + cr_Ct125*(cr_Ct226 + cr_Ct255 - cr_Ct257 - cr_Ct258 - cr_Ct260 + cr_Ct270));
  const double cr_Ct411 = cr_Ct301*cr_Ct385;
  const double cr_Ct412 = cr_Ct17*cr_Ct273;
  const double cr_Ct413 = cr_Ct228*(cr_Ct118*cr_Ct294 + cr_Ct125*(cr_Ct226 + cr_Ct288 - cr_Ct290 - cr_Ct291 - cr_Ct292 + cr_Ct296));
  const double cr_Ct414 = cr_Ct299*cr_Ct385;
  const double cr_Ct415 = 0.5*cr_Ct29;
  const double cr_Ct416 = cr_Ct415*r_strain[3];
  const double cr_Ct417 = cr_Ct337*cr_Ct366;
  const double cr_Ct418 = cr_Ct90*r_strain[3];
  const double cr_Ct419 = cr_Ct173*cr_Ct323;
  const double cr_Ct420 = cr_Ct95*r_strain[3];
  const double cr_Ct421 = 2.0*cr_Ct133;
  const double cr_Ct422 = cr_Ct421*cr_Ct53;
  const double cr_Ct423 = cr_Ct110*(cr_Ct133*cr_Ct195 + 0.66666666666666663*cr_Ct135*cr_Ct76 - cr_Ct147 + cr_Ct18*cr_Ct421 + cr_Ct20*cr_Ct421 - cr_Ct422*cr_Ct71 - cr_Ct422*cr_Ct82)/cr_Ct50;
  const double cr_Ct424 = cr_Ct351*cr_Ct366;
  const double cr_Ct425 = cr_Ct353/std::pow(cr_Ct25, 5.0/2.0);
  const double cr_Ct426 = cr_Ct14*cr_Ct20;
  const double cr_Ct427 = cr_Ct425*cr_Ct426;
  const double cr_Ct428 = 5.196152422706632*cr_Ct26;
  const double cr_Ct429 = cr_Ct343*cr_Ct428*(cr_Ct195 - 2.0*cr_Ct397);
  const double cr_Ct430 = cr_Ct172*cr_Ct210*cr_Ct3*cr_Ct369*cr_Ct46;
  const double cr_Ct431 = cr_Ct106*cr_Ct366;
  const double cr_Ct432 = cr_Ct205*cr_Ct53;
  const double cr_Ct433 = cr_Ct163*cr_Ct211*(-cr_Ct206 - cr_Ct207 - cr_Ct208 + cr_Ct432*cr_Ct71 + cr_Ct432*cr_Ct82 + 18.0*cr_Ct65);
  const double cr_Ct434 = 6.0*cr_Ct343;
  const double cr_Ct435 = cr_Ct108*(cr_Ct18*cr_Ct434 + cr_Ct19*cr_Ct434 + cr_Ct20*cr_Ct434 - 6.0*cr_Ct33 - 12.0*cr_Ct347 + 8.0*cr_Ct349 - 3.0*cr_Ct375 - cr_Ct39*cr_Ct434)/cr_Ct84;
  const double cr_Ct436 = cr_Ct17*cr_Ct18;
  const double cr_Ct437 = cr_Ct65*cr_Ct96;
  const double cr_Ct438 = cr_Ct219*cr_Ct50*cr_Ct86*std::sin(cr_Ct109)/std::sqrt(-cr_Ct33*std::pow(cr_Ct107 + cr_Ct33*cr_Ct89 - cr_Ct363*cr_Ct426 - cr_Ct363*cr_Ct436 - cr_Ct437, 2)/cr_Ct85 + 0.037037037037037035);
  const double cr_Ct439 = cr_Ct227*cr_Ct306;
  const double cr_Ct440 = cr_Ct416*cr_Ct439;
  const double cr_Ct441 = cr_Ct415*r_strain[4];
  const double cr_Ct442 = cr_Ct439*cr_Ct441;
  const double cr_Ct443 = cr_Ct417*r_strain[4];
  const double cr_Ct444 = cr_Ct366*cr_Ct419;
  const double cr_Ct445 = 0.86602540378443871*cr_Ct423;
  const double cr_Ct446 = cr_Ct424*r_strain[3];
  const double cr_Ct447 = cr_Ct425*r_strain[3];
  const double cr_Ct448 = cr_Ct425*cr_Ct436;
  const double cr_Ct449 = cr_Ct424*cr_Ct437;
  const double cr_Ct450 = cr_Ct398*cr_Ct431;
  const double cr_Ct451 = 2.598076211353316*cr_Ct435;
  const double cr_Ct452 = cr_Ct415*r_strain[5];
  const double cr_Ct453 = cr_Ct439*cr_Ct452;
  const double cr_Ct454 = cr_Ct417*r_strain[5];
  r_Ct(0,0)=cr_Ct30*(cr_Ct114 + cr_Ct116*cr_Ct229);
  r_Ct(0,1)=cr_Ct30*(cr_Ct235 + cr_Ct272*cr_Ct274);
  r_Ct(0,2)=cr_Ct30*(cr_Ct235 + cr_Ct274*cr_Ct297);
  r_Ct(0,3)=cr_Ct384*cr_Ct386;
  r_Ct(0,4)=cr_Ct386*cr_Ct402;
  r_Ct(0,5)=cr_Ct386*cr_Ct407;
  r_Ct(1,0)=cr_Ct30*(cr_Ct235 + cr_Ct408*cr_Ct409);
  r_Ct(1,1)=cr_Ct30*(cr_Ct114 + cr_Ct301*cr_Ct410);
  r_Ct(1,2)=cr_Ct30*(cr_Ct235 + cr_Ct297*cr_Ct409);
  r_Ct(1,3)=cr_Ct384*cr_Ct411;
  r_Ct(1,4)=cr_Ct402*cr_Ct411;
  r_Ct(1,5)=cr_Ct407*cr_Ct411;
  r_Ct(2,0)=cr_Ct30*(cr_Ct235 + cr_Ct408*cr_Ct412);
  r_Ct(2,1)=cr_Ct30*(cr_Ct235 + cr_Ct272*cr_Ct412);
  r_Ct(2,2)=cr_Ct30*(cr_Ct114 + cr_Ct299*cr_Ct413);
  r_Ct(2,3)=cr_Ct384*cr_Ct414;
  r_Ct(2,4)=cr_Ct402*cr_Ct414;
  r_Ct(2,5)=cr_Ct407*cr_Ct414;
  r_Ct(3,0)=-cr_Ct229*cr_Ct416;
  r_Ct(3,1)=-cr_Ct410*cr_Ct416;
  r_Ct(3,2)=-cr_Ct413*cr_Ct416;
  r_Ct(3,3)=cr_Ct415*(-cr_Ct228*r_strain[3]*(cr_Ct118*r_strain[3] - cr_Ct125*(cr_Ct342*cr_Ct423 + cr_Ct367*cr_Ct419 + cr_Ct417*cr_Ct418 + cr_Ct417*cr_Ct420 - cr_Ct438*(cr_Ct218*cr_Ct424*cr_Ct96 + cr_Ct340*cr_Ct427 - cr_Ct355*cr_Ct425 - cr_Ct356*cr_Ct366*cr_Ct418 + cr_Ct357*cr_Ct425*cr_Ct90 - cr_Ct372*cr_Ct431 + cr_Ct377*cr_Ct435 + cr_Ct420*cr_Ct429 + cr_Ct424*cr_Ct88 - cr_Ct430*r_strain[3] + cr_Ct433*r_strain[3]))) + cr_Ct234);
  r_Ct(3,4)=-cr_Ct402*cr_Ct440;
  r_Ct(3,5)=-cr_Ct407*cr_Ct440;
  r_Ct(4,0)=-cr_Ct229*cr_Ct441;
  r_Ct(4,1)=-cr_Ct410*cr_Ct441;
  r_Ct(4,2)=-cr_Ct413*cr_Ct441;
  r_Ct(4,3)=-cr_Ct384*cr_Ct442;
  r_Ct(4,4)=cr_Ct415*(-cr_Ct228*r_strain[4]*(cr_Ct118*r_strain[4] - cr_Ct125*(cr_Ct389*cr_Ct444 - cr_Ct438*(cr_Ct389*cr_Ct427 + cr_Ct389*cr_Ct448 - cr_Ct389*cr_Ct450 - cr_Ct393*cr_Ct447 - cr_Ct428*cr_Ct65*cr_Ct95*r_strain[4]*(-cr_Ct196 + cr_Ct238*cr_Ct93 + 2.0) - cr_Ct430*r_strain[4] + cr_Ct433*r_strain[4] + cr_Ct446*r_strain[5] + cr_Ct449*r_strain[4] + cr_Ct451*r_strain[4]) + cr_Ct443*cr_Ct90 + cr_Ct443*cr_Ct95 + cr_Ct445*r_strain[4])) + cr_Ct234);
  r_Ct(4,5)=-cr_Ct407*cr_Ct442;
  r_Ct(5,0)=-cr_Ct229*cr_Ct452;
  r_Ct(5,1)=-cr_Ct410*cr_Ct452;
  r_Ct(5,2)=-cr_Ct413*cr_Ct452;
  r_Ct(5,3)=-cr_Ct384*cr_Ct453;
  r_Ct(5,4)=-cr_Ct402*cr_Ct453;
  r_Ct(5,5)=cr_Ct415*(-cr_Ct228*r_strain[5]*(cr_Ct118*r_strain[5] - cr_Ct125*(cr_Ct403*cr_Ct444 - cr_Ct438*(-cr_Ct366*cr_Ct405*cr_Ct91 + cr_Ct403*cr_Ct448 - cr_Ct403*cr_Ct450 - cr_Ct404*cr_Ct447 + cr_Ct406*cr_Ct425*cr_Ct91 + cr_Ct429*cr_Ct95*r_strain[5] - cr_Ct430*r_strain[5] + cr_Ct433*r_strain[5] + cr_Ct446*r_strain[4] + cr_Ct449*r_strain[5] + cr_Ct451*r_strain[5]) + cr_Ct445*r_strain[5] + cr_Ct454*cr_Ct90 + cr_Ct454*cr_Ct95)) + cr_Ct234);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<SimoJuYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<VonMisesYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
    const auto &r_props = rValues.GetMaterialProperties();
    const double Young = r_props[YOUNG_MODULUS];
    const double nu = r_props[POISSON_RATIO];
    const double Gf = r_props[FRACTURE_ENERGY];
    const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
    double threshold;
    YieldSurfaceType::GetInitialUniaxialThreshold(rValues, threshold);
    auto &r_Ct = rValues.GetConstitutiveMatrix();
    const auto &r_strain = rValues.GetStrainVector();

    const double cr_Ct0 = nu - 1.0;
    const double cr_Ct1 = nu - 0.5;
    const double cr_Ct2 = 2*nu;
    const double cr_Ct3 = 1.0/(cr_Ct2 - 1);
    const double cr_Ct4 = nu + 1.0;
    const double cr_Ct5 = Young/cr_Ct4;
    const double cr_Ct6 = cr_Ct3*cr_Ct5;
    const double cr_Ct7 = cr_Ct1*cr_Ct6;
    const double cr_Ct8 = std::pow(cr_Ct7*r_strain[2], 2.0);
    const double cr_Ct9 = nu*r_strain[1];
    const double cr_Ct10 = cr_Ct0*r_strain[0];
    const double cr_Ct11 = cr_Ct10 - cr_Ct9;
    const double cr_Ct12 = nu*r_strain[0];
    const double cr_Ct13 = cr_Ct0*r_strain[1] - cr_Ct12;
    const double cr_Ct14 = cr_Ct6*(cr_Ct11 + cr_Ct13);
    const double cr_Ct15 = 0.5*cr_Ct9;
    const double cr_Ct16 = cr_Ct6*(-0.5*cr_Ct10 + cr_Ct13 + cr_Ct15);
    const double cr_Ct17 = 1.0 - nu;
    const double cr_Ct18 = cr_Ct17*r_strain[1];
    const double cr_Ct19 = cr_Ct17*r_strain[0];
    const double cr_Ct20 = cr_Ct19 + cr_Ct9;
    const double cr_Ct21 = 1.0/(1 - cr_Ct2);
    const double cr_Ct22 = cr_Ct21*cr_Ct5;
    const double cr_Ct23 = cr_Ct22*(-0.5*cr_Ct12 - 0.5*cr_Ct18 + cr_Ct20);
    const double cr_Ct24 = 0.055555555555555552*std::pow(cr_Ct14, 2.0) + 0.22222222222222224*std::pow(cr_Ct16, 2.0) + 0.22222222222222224*std::pow(cr_Ct23, 2.0) + cr_Ct8;
    const double cr_Ct25 = std::sqrt(cr_Ct24);
    const double cr_Ct26 = 0.57735026918962584*threshold;
    const double cr_Ct27 = cr_Ct26/cr_Ct25;
    const double cr_Ct28 = cr_Ct0*cr_Ct27;
    const double cr_Ct29 = cr_Ct12 + cr_Ct18;
    const double cr_Ct30 = 0.055555555555555552*std::pow(cr_Ct22*(cr_Ct20 + cr_Ct29), 1.0);
    const double cr_Ct31 = 3.0*nu;
    const double cr_Ct32 = 2.0 - cr_Ct31;
    const double cr_Ct33 = 0.11111111111111112*std::pow(cr_Ct23, 1.0);
    const double cr_Ct34 = cr_Ct31 - 1.0;
    const double cr_Ct35 = 0.11111111111111112*std::pow(cr_Ct22*(-cr_Ct15 - 0.5*cr_Ct19 + cr_Ct29), 1.0);
    const double cr_Ct36 = cr_Ct30 + cr_Ct32*cr_Ct33 + cr_Ct34*cr_Ct35;
    const double cr_Ct37 = cr_Ct26/std::pow(cr_Ct24, 3.0/2.0);
    const double cr_Ct38 = cr_Ct22*cr_Ct37;
    const double cr_Ct39 = cr_Ct11*cr_Ct38;
    const double cr_Ct40 = 0.055555555555555552*std::pow(cr_Ct14, 1.0)*cr_Ct3;
    const double cr_Ct41 = 0.11111111111111112*std::pow(cr_Ct16, 1.0);
    const double cr_Ct42 = cr_Ct21*(cr_Ct31 - 2.0);
    const double cr_Ct43 = cr_Ct3*cr_Ct34*cr_Ct41 + cr_Ct33*cr_Ct42 + cr_Ct40;
    const double cr_Ct44 = 1.0/cr_Ct24;
    const double cr_Ct45 = 1.0/(Gf*Young/(characteristic_length*std::pow(threshold, 2)) - 0.5);
    const double cr_Ct46 = cr_Ct45;
    const double cr_Ct47 = cr_Ct44*cr_Ct46;
    const double cr_Ct48 = cr_Ct47*cr_Ct5;
    const double cr_Ct49 = cr_Ct11*cr_Ct48;
    const double cr_Ct50 = std::exp(cr_Ct46*(-1.7320508075688772*cr_Ct25/threshold + 1.0));
    const double cr_Ct51 = cr_Ct50*cr_Ct6;
    const double cr_Ct52 = -cr_Ct27*nu;
    const double cr_Ct53 = cr_Ct33*cr_Ct34;
    const double cr_Ct54 = cr_Ct30 + cr_Ct32*cr_Ct35 + cr_Ct53;
    const double cr_Ct55 = cr_Ct3*cr_Ct53 + cr_Ct40 + cr_Ct41*cr_Ct42;
    const double cr_Ct56 = cr_Ct51*cr_Ct8*(cr_Ct37 + cr_Ct47)/r_strain[2];
    const double cr_Ct57 = cr_Ct13*cr_Ct38;
    const double cr_Ct58 = cr_Ct13*cr_Ct48;
    const double cr_Ct59 = cr_Ct21*cr_Ct37;
    const double cr_Ct60 = cr_Ct44*cr_Ct45;
    const double cr_Ct61 = std::pow(Young, 2)*cr_Ct1*cr_Ct3*cr_Ct50*r_strain[2]/std::pow(cr_Ct4, 2);
    r_Ct(0,0)=cr_Ct51*(cr_Ct28 - cr_Ct36*cr_Ct39 + cr_Ct43*cr_Ct49);
    r_Ct(0,1)=cr_Ct51*(-cr_Ct39*cr_Ct54 + cr_Ct49*cr_Ct55 + cr_Ct52);
    r_Ct(0,2)=-cr_Ct11*cr_Ct56;
    r_Ct(1,0)=cr_Ct51*(-cr_Ct36*cr_Ct57 + cr_Ct43*cr_Ct58 + cr_Ct52);
    r_Ct(1,1)=cr_Ct51*(cr_Ct28 - cr_Ct54*cr_Ct57 + cr_Ct55*cr_Ct58);
    r_Ct(1,2)=-cr_Ct13*cr_Ct56;
    r_Ct(2,0)=cr_Ct61*(-cr_Ct36*cr_Ct59 + cr_Ct43*cr_Ct60);
    r_Ct(2,1)=cr_Ct61*(-cr_Ct54*cr_Ct59 + cr_Ct55*cr_Ct60);
    r_Ct(2,2)=cr_Ct50*cr_Ct7*(cr_Ct27 - cr_Ct37*cr_Ct8 - cr_Ct47*cr_Ct8);

}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  const auto &r_props = rValues.GetMaterialProperties();
  const double Young = r_props[YOUNG_MODULUS];
  const double nu = r_props[POISSON_RATIO];
  const double Gf = r_props[FRACTURE_ENERGY];
  const double phi = r_props[FRICTION_ANGLE] * Globals::Pi / 180.0;
  const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
  auto &r_Ct = rValues.GetConstitutiveMatrix();
  const auto &r_strain = rValues.GetStrainVector();

  const bool has_symmetric_yield_stress = r_props.Has(YIELD_STRESS);
  const double threshold_compression = has_symmetric_yield_stress ? r_props[YIELD_STRESS] : r_props[YIELD_STRESS_COMPRESSION];
  const double threshold_tension = has_symmetric_yield_stress ? r_props[YIELD_STRESS] : r_props[YIELD_STRESS_TENSION];

  const double cr_Ct0 = nu - 1.0;
  const double cr_Ct1 = std::cos(phi);
  const double cr_Ct2 = std::tan(0.5*phi + 0.78539816339744828);
  const double cr_Ct3 = 0.5*cr_Ct1*threshold_compression/cr_Ct2;
  const double cr_Ct4 = cr_Ct0*cr_Ct3;
  const double cr_Ct5 = 2*nu;
  const double cr_Ct6 = 1.0/(cr_Ct5 - 1);
  const double cr_Ct7 = std::abs(threshold_compression/threshold_tension);
  const double cr_Ct8 = cr_Ct7/std::pow(cr_Ct2, 2);
  const double cr_Ct9 = std::sin(phi);
  const double cr_Ct10 = cr_Ct8 + 1;
  const double cr_Ct11 = cr_Ct10*cr_Ct9;
  const double cr_Ct12 = 0.16666666666666666*cr_Ct11 + 0.16666666666666666*cr_Ct8 - 0.16666666666666666;
  const double cr_Ct13 = cr_Ct12*cr_Ct6;
  const double cr_Ct14 = -cr_Ct13;
  const double cr_Ct15 = nu*r_strain[1];
  const double cr_Ct16 = cr_Ct0*r_strain[0];
  const double cr_Ct17 = -cr_Ct15 + cr_Ct16;
  const double cr_Ct18 = nu*r_strain[0];
  const double cr_Ct19 = cr_Ct0*r_strain[1];
  const double cr_Ct20 = -cr_Ct18 + cr_Ct19;
  const double cr_Ct21 = nu + 1.0;
  const double cr_Ct22 = Young/cr_Ct21;
  const double cr_Ct23 = cr_Ct22*cr_Ct6;
  const double cr_Ct24 = cr_Ct23*(cr_Ct17 + cr_Ct20);
  const double cr_Ct25 = 0.055555555555555552*std::pow(cr_Ct24, 1.0)*cr_Ct6;
  const double cr_Ct26 = 3.0*nu;
  const double cr_Ct27 = cr_Ct26 - 1.0;
  const double cr_Ct28 = 0.5*cr_Ct15;
  const double cr_Ct29 = cr_Ct23*(-0.5*cr_Ct16 + cr_Ct20 + cr_Ct28);
  const double cr_Ct30 = std::pow(cr_Ct29, 1.0);
  const double cr_Ct31 = -nu;
  const double cr_Ct32 = cr_Ct31 + 1.0;
  const double cr_Ct33 = cr_Ct32*r_strain[1];
  const double cr_Ct34 = cr_Ct32*r_strain[0];
  const double cr_Ct35 = cr_Ct15 + cr_Ct34;
  const double cr_Ct36 = 1 - cr_Ct5;
  const double cr_Ct37 = 1.0/cr_Ct36;
  const double cr_Ct38 = cr_Ct22*cr_Ct37;
  const double cr_Ct39 = cr_Ct38*(-0.5*cr_Ct18 - 0.5*cr_Ct33 + cr_Ct35);
  const double cr_Ct40 = std::pow(cr_Ct39, 1.0);
  const double cr_Ct41 = 0.11111111111111113*cr_Ct37*(cr_Ct26 - 2.0);
  const double cr_Ct42 = -cr_Ct8;
  const double cr_Ct43 = cr_Ct42 + 1;
  const double cr_Ct44 = cr_Ct10 - cr_Ct43*cr_Ct9;
  const double cr_Ct45 = cr_Ct31 + 0.5;
  const double cr_Ct46 = std::pow(cr_Ct38*cr_Ct45*r_strain[2], 2.0);
  const double cr_Ct47 = cr_Ct18 + cr_Ct33;
  const double cr_Ct48 = cr_Ct38*(cr_Ct35 + cr_Ct47);
  const double cr_Ct49 = std::pow(cr_Ct48, 2.0);
  const double cr_Ct50 = std::pow(cr_Ct39, 2.0);
  const double cr_Ct51 = 0.44444444444444453*cr_Ct50;
  const double cr_Ct52 = cr_Ct38*(-cr_Ct28 - 0.5*cr_Ct34 + cr_Ct47);
  const double cr_Ct53 = std::pow(cr_Ct52, 2.0);
  const double cr_Ct54 = 1.0/(2.0*cr_Ct46 + 0.1111111111111111*cr_Ct49 + cr_Ct51 + 0.44444444444444453*cr_Ct53);
  const double cr_Ct55 = 0.22222222222222227*cr_Ct50;
  const double cr_Ct56 = cr_Ct46 + 0.055555555555555552*cr_Ct49 + 0.22222222222222227*cr_Ct53 + cr_Ct55;
  const double cr_Ct57 = std::pow(cr_Ct56, -1.0/2.0);
  const double cr_Ct58 = std::pow(r_strain[2], 2);
  const double cr_Ct59 = std::pow(cr_Ct45, 2);
  const double cr_Ct60 = 0.33333333333333331*cr_Ct18;
  const double cr_Ct61 = 0.66666666666666674*cr_Ct15;
  const double cr_Ct62 = 0.66666666666666674*cr_Ct34;
  const double cr_Ct63 = 0.33333333333333331*cr_Ct33;
  const double cr_Ct64 = -cr_Ct60 + cr_Ct61 + cr_Ct62 - cr_Ct63;
  const double cr_Ct65 = 0.66666666666666674*cr_Ct18;
  const double cr_Ct66 = 0.33333333333333331*cr_Ct15;
  const double cr_Ct67 = cr_Ct58*cr_Ct59 - cr_Ct64*(0.66666666666666674*cr_Ct33 - 0.33333333333333331*cr_Ct34 + cr_Ct65 - cr_Ct66);
  const double cr_Ct68 = 5.196152422706632*cr_Ct67;
  const double cr_Ct69 = std::pow(Young, 2)/std::pow(cr_Ct21, 2);
  const double cr_Ct70 = cr_Ct69/std::pow(cr_Ct36, 2);
  const double cr_Ct71 = 0.33333333333333331*std::asin(cr_Ct54*cr_Ct57*cr_Ct68*cr_Ct70);
  const double cr_Ct72 = std::cos(cr_Ct71);
  const double cr_Ct73 = 0.5*cr_Ct44*cr_Ct72;
  const double cr_Ct74 = std::sin(cr_Ct71);
  const double cr_Ct75 = cr_Ct43/cr_Ct9;
  const double cr_Ct76 = cr_Ct9*(cr_Ct10 - cr_Ct75);
  const double cr_Ct77 = 0.28867513459481292*cr_Ct74*cr_Ct76;
  const double cr_Ct78 = -cr_Ct73 + cr_Ct77;
  const double cr_Ct79 = nu - 0.5;
  const double cr_Ct80 = cr_Ct23*cr_Ct79;
  const double cr_Ct81 = std::pow(cr_Ct80*r_strain[2], 2.0);
  const double cr_Ct82 = std::pow(cr_Ct24, 2.0);
  const double cr_Ct83 = std::pow(cr_Ct29, 2.0);
  const double cr_Ct84 = cr_Ct55 + cr_Ct81 + 0.055555555555555552*cr_Ct82 + 0.22222222222222227*cr_Ct83;
  const double cr_Ct85 = std::sqrt(cr_Ct84);
  const double cr_Ct86 = 1.0/cr_Ct85;
  const double cr_Ct87 = cr_Ct78*cr_Ct86;
  const double cr_Ct88 = cr_Ct87*(cr_Ct25 + 0.11111111111111113*cr_Ct27*cr_Ct30*cr_Ct6 + cr_Ct40*cr_Ct41);
  const double cr_Ct89 = cr_Ct6*cr_Ct69;
  const double cr_Ct90 = cr_Ct60 - cr_Ct61 - cr_Ct62 + cr_Ct63;
  const double cr_Ct91 = cr_Ct37*cr_Ct90*(-0.33333333333333331*cr_Ct16 + 0.66666666666666674*cr_Ct19 - cr_Ct65 + cr_Ct66) + cr_Ct58*cr_Ct6*std::pow(cr_Ct79, 2);
  const double cr_Ct92 = 1.0/(cr_Ct51 + 2.0*cr_Ct81 + 0.1111111111111111*cr_Ct82 + 0.44444444444444453*cr_Ct83);
  const double cr_Ct93 = 5.196152422706632*cr_Ct92;
  const double cr_Ct94 = cr_Ct91*cr_Ct93;
  const double cr_Ct95 = 0.33333333333333331*std::asin(cr_Ct86*cr_Ct89*cr_Ct94);
  const double cr_Ct96 = 0.064150029909958411*cr_Ct44;
  const double cr_Ct97 = 0.037037037037037035*cr_Ct76*std::cos(cr_Ct95) + cr_Ct96*std::sin(cr_Ct95);
  const double cr_Ct98 = nu;
  const double cr_Ct99 = 3.4641016151377553*cr_Ct18;
  const double cr_Ct100 = 1.7320508075688772*cr_Ct15;
  const double cr_Ct101 = 1.7320508075688772*cr_Ct34;
  const double cr_Ct102 = 3.4641016151377553*cr_Ct33;
  const double cr_Ct103 = -cr_Ct100 - cr_Ct101 + cr_Ct102 + cr_Ct99;
  const double cr_Ct104 = 5.196152422706632*nu;
  const double cr_Ct105 = std::pow(cr_Ct84, -2);
  const double cr_Ct106 = std::pow(cr_Ct48, 1.0);
  const double cr_Ct107 = 0.22222222222222221*cr_Ct106;
  const double cr_Ct108 = 2.0 - cr_Ct26;
  const double cr_Ct109 = cr_Ct108*cr_Ct40;
  const double cr_Ct110 = std::pow(cr_Ct52, 1.0);
  const double cr_Ct111 = cr_Ct110*cr_Ct27;
  const double cr_Ct112 = 1.299038105676658*cr_Ct107 + 0.57735026918962584*cr_Ct109 + 0.57735026918962584*cr_Ct111;
  const double cr_Ct113 = 0.055555555555555552*cr_Ct106;
  const double cr_Ct114 = 0.11111111111111113*cr_Ct109 + 0.11111111111111113*cr_Ct111 + cr_Ct113;
  const double cr_Ct115 = 1.0/cr_Ct84;
  const double cr_Ct116 = std::pow(-std::pow(Young, 4)*std::pow(cr_Ct67, 2)/(std::pow(cr_Ct21, 4)*std::pow(cr_Ct36, 4)*std::pow(cr_Ct56, 3)) + 0.14814814814814814, -1.0/2.0);
  const double cr_Ct117 = cr_Ct116*cr_Ct38;
  const double cr_Ct118 = cr_Ct117*cr_Ct97*(-cr_Ct105*cr_Ct112*cr_Ct23*cr_Ct91 - cr_Ct114*cr_Ct115*cr_Ct23*cr_Ct94 + cr_Ct37*cr_Ct92*(cr_Ct103*(cr_Ct98 - 0.66666666666666674) + cr_Ct90*(cr_Ct104 - 1.7320508075688772)));
  const double cr_Ct119 = -cr_Ct118 + cr_Ct14 + cr_Ct88;
  const double cr_Ct120 = 1.0/(Gf*Young*std::pow(cr_Ct7, 2)/(characteristic_length*std::pow(threshold_compression, 2)) - 0.5);
  const double cr_Ct121 = cr_Ct120;
  const double cr_Ct122 = cr_Ct121*cr_Ct22;
  const double cr_Ct123 = cr_Ct122*cr_Ct17;
  const double cr_Ct124 = cr_Ct57*(cr_Ct73 - cr_Ct77);
  const double cr_Ct125 = cr_Ct100 + cr_Ct101 - cr_Ct102 - cr_Ct99;
  const double cr_Ct126 = std::pow(cr_Ct56, -2);
  const double cr_Ct127 = cr_Ct126*cr_Ct38*cr_Ct67;
  const double cr_Ct128 = cr_Ct54/cr_Ct56;
  const double cr_Ct129 = cr_Ct128*cr_Ct38*cr_Ct68;
  const double cr_Ct130 = 0.037037037037037035*cr_Ct72*cr_Ct9*(cr_Ct42 + cr_Ct75 - 1) - cr_Ct74*cr_Ct96;
  const double cr_Ct131 = cr_Ct117*cr_Ct130;
  const double cr_Ct132 = cr_Ct114*cr_Ct124 + cr_Ct12 + cr_Ct131*(-cr_Ct112*cr_Ct127 - cr_Ct114*cr_Ct129 + cr_Ct54*(cr_Ct125*(0.66666666666666674 - cr_Ct98) + cr_Ct64*(1.7320508075688772 - cr_Ct104)));
  const double cr_Ct133 = 0.16666666666666666*cr_Ct24*(cr_Ct11 + cr_Ct8 - 1) - cr_Ct78*cr_Ct85;
  const double cr_Ct134 = 1.0/cr_Ct133;
  const double cr_Ct135 = cr_Ct134*cr_Ct3;
  const double cr_Ct136 = cr_Ct135*cr_Ct38;
  const double cr_Ct137 = cr_Ct136*cr_Ct17;
  const double cr_Ct138 = cr_Ct134*std::exp(cr_Ct121*(1.0 - 2.0*cr_Ct133*cr_Ct2/(cr_Ct1*threshold_compression)));
  const double cr_Ct139 = cr_Ct138*cr_Ct23;
  const double cr_Ct140 = cr_Ct3*nu;
  const double cr_Ct141 = cr_Ct27*cr_Ct40;
  const double cr_Ct142 = 0.11111111111111113*cr_Ct141;
  const double cr_Ct143 = cr_Ct87*(cr_Ct142*cr_Ct6 + cr_Ct25 + cr_Ct30*cr_Ct41);
  const double cr_Ct144 = cr_Ct98 - 0.33333333333333331;
  const double cr_Ct145 = cr_Ct104 - 3.4641016151377553;
  const double cr_Ct146 = cr_Ct22*cr_Ct91;
  const double cr_Ct147 = cr_Ct108*cr_Ct110;
  const double cr_Ct148 = 1.299038105676658*cr_Ct107 + 0.57735026918962584*cr_Ct141 + 0.57735026918962584*cr_Ct147;
  const double cr_Ct149 = cr_Ct113 + cr_Ct142 + 0.11111111111111113*cr_Ct147;
  const double cr_Ct150 = cr_Ct116*cr_Ct23*cr_Ct37*cr_Ct97*(-cr_Ct105*cr_Ct146*cr_Ct148 - cr_Ct115*cr_Ct146*cr_Ct149*cr_Ct93 + cr_Ct92*(cr_Ct103*cr_Ct144 + cr_Ct145*cr_Ct90));
  const double cr_Ct151 = cr_Ct14 + cr_Ct143 - cr_Ct150;
  const double cr_Ct152 = cr_Ct12 + cr_Ct124*cr_Ct149 + cr_Ct131*(-cr_Ct127*cr_Ct148 - cr_Ct129*cr_Ct149 + cr_Ct54*(cr_Ct125*cr_Ct144 + cr_Ct145*cr_Ct64));
  const double cr_Ct153 = 1.0/r_strain[2];
  const double cr_Ct154 = cr_Ct153;
  const double cr_Ct155 = -cr_Ct154*cr_Ct81*cr_Ct87;
  const double cr_Ct156 = 10.392304845413264*cr_Ct54*cr_Ct59*r_strain[2];
  const double cr_Ct157 = cr_Ct153*cr_Ct46*cr_Ct68;
  const double cr_Ct158 = cr_Ct126*cr_Ct157;
  const double cr_Ct159 = cr_Ct128*cr_Ct157;
  const double cr_Ct160 = cr_Ct116*cr_Ct70;
  const double cr_Ct161 = cr_Ct160*cr_Ct57*cr_Ct85*cr_Ct97*(-cr_Ct156 + cr_Ct158 + cr_Ct159);
  const double cr_Ct162 = cr_Ct135*(cr_Ct124*cr_Ct154*cr_Ct46 + cr_Ct130*cr_Ct160*(cr_Ct156 - cr_Ct158 - cr_Ct159));
  const double cr_Ct163 = cr_Ct139*(cr_Ct120*(cr_Ct155 + cr_Ct161) + cr_Ct162);
  const double cr_Ct164 = cr_Ct122*cr_Ct20;
  const double cr_Ct165 = cr_Ct136*cr_Ct20;
  const double cr_Ct166 = -cr_Ct13;
  const double cr_Ct167 = cr_Ct135*cr_Ct37;
  const double cr_Ct168 = cr_Ct138*cr_Ct79*cr_Ct89*r_strain[2];
  r_Ct(0,0)=cr_Ct139*(-cr_Ct119*cr_Ct123 - cr_Ct132*cr_Ct137 + cr_Ct4);
  r_Ct(0,1)=-cr_Ct139*(cr_Ct123*cr_Ct151 + cr_Ct137*cr_Ct152 + cr_Ct140);
  r_Ct(0,2)=-cr_Ct163*cr_Ct17;
  r_Ct(1,0)=-cr_Ct139*(cr_Ct119*cr_Ct164 + cr_Ct132*cr_Ct165 + cr_Ct140);
  r_Ct(1,1)=cr_Ct139*(-cr_Ct151*cr_Ct164 - cr_Ct152*cr_Ct165 + cr_Ct4);
  r_Ct(1,2)=-cr_Ct163*cr_Ct20;
  r_Ct(2,0)=-cr_Ct168*(cr_Ct120*(-cr_Ct118 + cr_Ct166 + cr_Ct88) + cr_Ct132*cr_Ct167);
  r_Ct(2,1)=-cr_Ct168*(cr_Ct120*(cr_Ct143 - cr_Ct150 + cr_Ct166) + cr_Ct152*cr_Ct167);
  r_Ct(2,2)=cr_Ct138*cr_Ct80*(-cr_Ct121*r_strain[2]*(cr_Ct155 + cr_Ct161) - cr_Ct162*r_strain[2] + cr_Ct3);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<TrescaYieldSurface<TrescaPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  const auto &r_props = rValues.GetMaterialProperties();
  const double Young = r_props[YOUNG_MODULUS];
  const double nu = r_props[POISSON_RATIO];
  const double Gf = r_props[FRACTURE_ENERGY];
  const double phi = r_props[FRICTION_ANGLE];
  const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
  const double threshold = r_props[YIELD_STRESS];
  auto &r_Ct = rValues.GetConstitutiveMatrix();
  const auto &r_strain = rValues.GetStrainVector();
  const double sin_phi = std::sin(phi * Globals::Pi / 180.0);

  const double cr_Ct0 = nu - 1.0;
  const double cr_Ct1 = nu*r_strain[1];
  const double cr_Ct2 = cr_Ct0*r_strain[0];
  const double cr_Ct3 = -cr_Ct1 + cr_Ct2;
  const double cr_Ct4 = nu*r_strain[0];
  const double cr_Ct5 = cr_Ct0*r_strain[1] - cr_Ct4;
  const double cr_Ct6 = 2*nu;
  const double cr_Ct7 = 1.0/(cr_Ct6 - 1);
  const double cr_Ct8 = nu + 1.0;
  const double cr_Ct9 = Young/cr_Ct8;
  const double cr_Ct10 = cr_Ct7*cr_Ct9;
  const double cr_Ct11 = cr_Ct10*(cr_Ct3 + cr_Ct5);
  const double cr_Ct12 = 1.7320508075688772*sin_phi;
  const double cr_Ct13 = cr_Ct12 - 5.196152422706632;
  const double cr_Ct14 = 1.0/cr_Ct13;
  const double cr_Ct15 = 2.0*sin_phi;
  const double cr_Ct16 = cr_Ct14*cr_Ct15;
  const double cr_Ct17 = nu - 0.5;
  const double cr_Ct18 = cr_Ct10*cr_Ct17;
  const double cr_Ct19 = std::pow(cr_Ct18*r_strain[2], 2.0);
  const double cr_Ct20 = 0.5*cr_Ct1;
  const double cr_Ct21 = cr_Ct10*(-0.5*cr_Ct2 + cr_Ct20 + cr_Ct5);
  const double cr_Ct22 = -nu;
  const double cr_Ct23 = cr_Ct22 + 1.0;
  const double cr_Ct24 = cr_Ct23*r_strain[1];
  const double cr_Ct25 = cr_Ct23*r_strain[0];
  const double cr_Ct26 = cr_Ct1 + cr_Ct25;
  const double cr_Ct27 = 1.0/(1 - cr_Ct6);
  const double cr_Ct28 = cr_Ct27*cr_Ct9;
  const double cr_Ct29 = cr_Ct28*(-0.5*cr_Ct24 + cr_Ct26 - 0.5*cr_Ct4);
  const double cr_Ct30 = 0.22222222222222227*std::pow(cr_Ct29, 2.0);
  const double cr_Ct31 = std::sqrt(0.055555555555555552*std::pow(cr_Ct11, 2.0) + cr_Ct19 + 0.22222222222222227*std::pow(cr_Ct21, 2.0) + cr_Ct30);
  const double cr_Ct32 = 1.0/(-cr_Ct11*cr_Ct16 + cr_Ct31);
  const double cr_Ct33 = sin_phi - 1;
  const double cr_Ct34 = 1.0/cr_Ct33;
  const double cr_Ct35 = std::abs(cr_Ct34*threshold*(sin_phi + 3.0));
  const double cr_Ct36 = cr_Ct14*cr_Ct33*cr_Ct35;
  const double cr_Ct37 = cr_Ct32*cr_Ct36;
  const double cr_Ct38 = cr_Ct0*cr_Ct37;
  const double cr_Ct39 = cr_Ct16*cr_Ct7;
  const double cr_Ct40 = 1.0/cr_Ct31;
  const double cr_Ct41 = 0.055555555555555552*std::pow(cr_Ct11, 1.0)*cr_Ct7;
  const double cr_Ct42 = 3.0*nu;
  const double cr_Ct43 = cr_Ct42 - 1.0;
  const double cr_Ct44 = 0.11111111111111113*std::pow(cr_Ct21, 1.0);
  const double cr_Ct45 = 0.11111111111111113*std::pow(cr_Ct29, 1.0);
  const double cr_Ct46 = cr_Ct27*(cr_Ct42 - 2.0);
  const double cr_Ct47 = cr_Ct40*(cr_Ct41 + cr_Ct43*cr_Ct44*cr_Ct7 + cr_Ct45*cr_Ct46);
  const double cr_Ct48 = cr_Ct39 - cr_Ct47;
  const double cr_Ct49 = 1.0/(Gf*Young/(characteristic_length*std::pow(threshold, 2)) - 0.5);
  const double cr_Ct50 = cr_Ct49;
  const double cr_Ct51 = cr_Ct32*cr_Ct50;
  const double cr_Ct52 = cr_Ct51*cr_Ct9;
  const double cr_Ct53 = cr_Ct3*cr_Ct52;
  const double cr_Ct54 = cr_Ct15/(5.196152422706632 - cr_Ct12);
  const double cr_Ct55 = cr_Ct24 + cr_Ct4;
  const double cr_Ct56 = cr_Ct28*(cr_Ct26 + cr_Ct55);
  const double cr_Ct57 = 0.055555555555555552*std::pow(cr_Ct56, 1.0);
  const double cr_Ct58 = 2.0 - cr_Ct42;
  const double cr_Ct59 = cr_Ct28*(-cr_Ct20 - 0.5*cr_Ct25 + cr_Ct55);
  const double cr_Ct60 = 0.11111111111111113*std::pow(cr_Ct59, 1.0);
  const double cr_Ct61 = std::sqrt(cr_Ct30 + 0.055555555555555552*std::pow(cr_Ct56, 2.0) + 0.22222222222222227*std::pow(cr_Ct59, 2.0) + std::pow(cr_Ct28*r_strain[2]*(cr_Ct22 + 0.5), 2.0));
  const double cr_Ct62 = 1.0/cr_Ct61;
  const double cr_Ct63 = cr_Ct54 + cr_Ct62*(cr_Ct43*cr_Ct60 + cr_Ct45*cr_Ct58 + cr_Ct57);
  const double cr_Ct64 = 0.25*cr_Ct36/std::pow(-cr_Ct11*cr_Ct14*sin_phi + 0.5*cr_Ct31, 2);
  const double cr_Ct65 = cr_Ct28*cr_Ct64;
  const double cr_Ct66 = cr_Ct3*cr_Ct65;
  const double cr_Ct67 = std::exp(-cr_Ct50*(cr_Ct13*cr_Ct34*(cr_Ct54*cr_Ct56 + cr_Ct61)/cr_Ct35 - 1));
  const double cr_Ct68 = cr_Ct10*cr_Ct67;
  const double cr_Ct69 = cr_Ct37*nu;
  const double cr_Ct70 = cr_Ct43*cr_Ct45;
  const double cr_Ct71 = cr_Ct40*(cr_Ct41 + cr_Ct44*cr_Ct46 + cr_Ct7*cr_Ct70);
  const double cr_Ct72 = cr_Ct39 - cr_Ct71;
  const double cr_Ct73 = cr_Ct54 + cr_Ct62*(cr_Ct57 + cr_Ct58*cr_Ct60 + cr_Ct70);
  const double cr_Ct74 = cr_Ct19*cr_Ct40;
  const double cr_Ct75 = cr_Ct68*cr_Ct74*(cr_Ct51 + cr_Ct64)/r_strain[2];
  const double cr_Ct76 = cr_Ct5*cr_Ct52;
  const double cr_Ct77 = cr_Ct5*cr_Ct65;
  const double cr_Ct78 = cr_Ct32*cr_Ct49;
  const double cr_Ct79 = cr_Ct27*cr_Ct64;
  const double cr_Ct80 = std::pow(Young, 2)*cr_Ct17*cr_Ct67*cr_Ct7*r_strain[2]/std::pow(cr_Ct8, 2);
  r_Ct(0,0)=cr_Ct68*(cr_Ct38 - cr_Ct48*cr_Ct53 - cr_Ct63*cr_Ct66);
  r_Ct(0,1)=-cr_Ct68*(cr_Ct53*cr_Ct72 + cr_Ct66*cr_Ct73 + cr_Ct69);
  r_Ct(0,2)=-cr_Ct3*cr_Ct75;
  r_Ct(1,0)=-cr_Ct68*(cr_Ct48*cr_Ct76 + cr_Ct63*cr_Ct77 + cr_Ct69);
  r_Ct(1,1)=cr_Ct68*(cr_Ct38 - cr_Ct72*cr_Ct76 - cr_Ct73*cr_Ct77);
  r_Ct(1,2)=-cr_Ct5*cr_Ct75;
  r_Ct(2,0)=-cr_Ct80*(cr_Ct63*cr_Ct79 + cr_Ct78*(cr_Ct39 - cr_Ct47));
  r_Ct(2,1)=-cr_Ct80*(cr_Ct73*cr_Ct79 + cr_Ct78*(cr_Ct39 - cr_Ct71));
  r_Ct(2,2)=cr_Ct18*cr_Ct67*(cr_Ct37 - cr_Ct51*cr_Ct74 - cr_Ct64*cr_Ct74);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<RankineYieldSurface<RankinePlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  const auto &r_props = rValues.GetMaterialProperties();
  const double Young = r_props[YOUNG_MODULUS];
  const double nu = r_props[POISSON_RATIO];
  const double Gf = r_props[FRACTURE_ENERGY];
  const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
  const double threshold = r_props[YIELD_STRESS];
  auto &r_Ct = rValues.GetConstitutiveMatrix();
  const auto &r_strain = rValues.GetStrainVector();

	const double cr_Ct0 = nu - 1.0;
	const double cr_Ct1 = cr_Ct0*threshold;
	const double cr_Ct2 = nu*r_strain[1];
	const double cr_Ct3 = -cr_Ct2;
	const double cr_Ct4 = cr_Ct0*r_strain[0];
	const double cr_Ct5 = cr_Ct3 + cr_Ct4;
	const double cr_Ct6 = nu + 1.0;
	const double cr_Ct7 = Young/cr_Ct6;
	const double cr_Ct8 = cr_Ct5*cr_Ct7;
	const double cr_Ct9 = 0.5*nu;
	const double cr_Ct10 = 2*nu;
	const double cr_Ct11 = cr_Ct10 - 1;
	const double cr_Ct12 = 1.0/cr_Ct11;
	const double cr_Ct13 = cr_Ct12*(cr_Ct9 - 0.5);
	const double cr_Ct14 = -cr_Ct10;
	const double cr_Ct15 = cr_Ct14 + 1;
	const double cr_Ct16 = std::pow(cr_Ct15, -2);
	const double cr_Ct17 = 0.25*cr_Ct16;
	const double cr_Ct18 = std::pow(r_strain[2], 2);
	const double cr_Ct19 = nu - 0.5;
	const double cr_Ct20 = std::pow(cr_Ct19, 2);
	const double cr_Ct21 = cr_Ct18*cr_Ct20/std::pow(cr_Ct11, 2);
	const double cr_Ct22 = nu*r_strain[0];
	const double cr_Ct23 = -nu;
	const double cr_Ct24 = cr_Ct23 + 1.0;
	const double cr_Ct25 = cr_Ct24*r_strain[1];
	const double cr_Ct26 = cr_Ct24*r_strain[0];
	const double cr_Ct27 = cr_Ct22 + cr_Ct25 - cr_Ct26 + cr_Ct3;
	const double cr_Ct28 = std::pow(Young, 2)/std::pow(cr_Ct6, 2);
	const double cr_Ct29 = std::sqrt(cr_Ct28*(cr_Ct17*std::pow(cr_Ct27, 2) + cr_Ct21));
	const double cr_Ct30 = 1.0/cr_Ct29;
	const double cr_Ct31 = cr_Ct7*(cr_Ct10 - 1.0);
	const double cr_Ct32 = cr_Ct27*cr_Ct30*cr_Ct31;
	const double cr_Ct33 = -cr_Ct12*cr_Ct9 + cr_Ct17*cr_Ct32;
	const double cr_Ct34 = 1.0/(Gf*Young/(characteristic_length*std::pow(threshold, 2)) - 0.5);
	const double cr_Ct35 = cr_Ct34;
	const double cr_Ct36 = cr_Ct35*(cr_Ct13 + cr_Ct33);
	const double cr_Ct37 = -cr_Ct22;
	const double cr_Ct38 = cr_Ct2 - cr_Ct25 + cr_Ct26 + cr_Ct37;
	const double cr_Ct39 = 1.0/cr_Ct15;
	const double cr_Ct40 = 0.25*cr_Ct39;
	const double cr_Ct41 = cr_Ct38*cr_Ct40/std::sqrt(cr_Ct16*cr_Ct28*(cr_Ct18*std::pow(cr_Ct23 + 0.5, 2) + 0.25*std::pow(cr_Ct38, 2)));
	const double cr_Ct42 = cr_Ct12*cr_Ct7;
	const double cr_Ct43 = 0.5*cr_Ct42;
	const double cr_Ct44 = cr_Ct0*r_strain[1];
	const double cr_Ct45 = -cr_Ct2*cr_Ct43 - cr_Ct22*cr_Ct43 + cr_Ct29 + cr_Ct4*cr_Ct43 + cr_Ct43*cr_Ct44;
	const double cr_Ct46 = 1.0/cr_Ct45;
	const double cr_Ct47 = cr_Ct46*threshold;
	const double cr_Ct48 = cr_Ct39*cr_Ct47;
	const double cr_Ct49 = cr_Ct48*(cr_Ct41*cr_Ct7*(cr_Ct14 + 1.0) + 0.5);
	const double cr_Ct50 = cr_Ct46*std::exp(cr_Ct35*(-cr_Ct45/threshold + 1.0));
	const double cr_Ct51 = cr_Ct42*cr_Ct50;
	const double cr_Ct52 = nu*threshold;
	const double cr_Ct53 = cr_Ct12*(cr_Ct32*cr_Ct40 - 0.5);
	const double cr_Ct54 = cr_Ct35*cr_Ct53;
	const double cr_Ct55 = cr_Ct48*(cr_Ct31*cr_Ct41 + 0.5);
	const double cr_Ct56 = cr_Ct50*r_strain[2];
	const double cr_Ct57 = std::pow(Young, 3)*cr_Ct20*cr_Ct30*cr_Ct56*(cr_Ct35 + cr_Ct47)/(std::pow(cr_Ct11, 3)*std::pow(cr_Ct6, 3));
	const double cr_Ct58 = cr_Ct37 + cr_Ct44;
	const double cr_Ct59 = cr_Ct58*cr_Ct7;
	const double cr_Ct60 = cr_Ct12*cr_Ct19*cr_Ct28*cr_Ct56;
	const double cr_Ct61 = cr_Ct21*cr_Ct28*cr_Ct30;
	r_Ct(0,0)=cr_Ct51*(cr_Ct1 - cr_Ct36*cr_Ct8 - cr_Ct49*cr_Ct8);
	r_Ct(0,1)=-cr_Ct51*(cr_Ct52 + cr_Ct54*cr_Ct8 + cr_Ct55*cr_Ct8);
	r_Ct(0,2)=-cr_Ct5*cr_Ct57;
	r_Ct(1,0)=-cr_Ct51*(cr_Ct36*cr_Ct59 + cr_Ct49*cr_Ct59 + cr_Ct52);
	r_Ct(1,1)=cr_Ct51*(cr_Ct1 - cr_Ct54*cr_Ct59 - cr_Ct55*cr_Ct59);
	r_Ct(1,2)=-cr_Ct57*cr_Ct58;
	r_Ct(2,0)=-cr_Ct60*(cr_Ct34*(cr_Ct13 + cr_Ct33) + cr_Ct49);
	r_Ct(2,1)=-cr_Ct60*(cr_Ct34*cr_Ct53 + cr_Ct55);
	r_Ct(2,2)=cr_Ct19*cr_Ct51*(-cr_Ct35*cr_Ct61 - cr_Ct47*cr_Ct61 + threshold);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<SimoJuYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

template<>
void AutomaticDifferentiationTangentUtilities<VonMisesYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  const auto &r_props = rValues.GetMaterialProperties();
  const double Young = r_props[YOUNG_MODULUS];
  const double nu = r_props[POISSON_RATIO];
  const double Gf = r_props[FRACTURE_ENERGY];
  const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
  double threshold;
  YieldSurfaceType::GetInitialUniaxialThreshold(rValues, threshold);
  auto &r_Ct = rValues.GetConstitutiveMatrix();
  const auto &r_strain = rValues.GetStrainVector();

  const double cr_Ct0 = nu - 1.0;
  const double cr_Ct1 = 1.0/(1.0 - 0.5*characteristic_length*std::pow(threshold, 2)/(Gf*Young));
  const double cr_Ct2 = nu - 0.5;
  const double cr_Ct3 = 2*nu;
  const double cr_Ct4 = 1.0/(cr_Ct3 - 1);
  const double cr_Ct5 = nu + 1.0;
  const double cr_Ct6 = Young/cr_Ct5;
  const double cr_Ct7 = cr_Ct4*cr_Ct6;
  const double cr_Ct8 = cr_Ct2*cr_Ct7;
  const double cr_Ct9 = std::pow(cr_Ct8*r_strain[2], 2.0);
  const double cr_Ct10 = nu*r_strain[1];
  const double cr_Ct11 = cr_Ct0*r_strain[0];
  const double cr_Ct12 = -cr_Ct10 + cr_Ct11;
  const double cr_Ct13 = nu*r_strain[0];
  const double cr_Ct14 = cr_Ct0*r_strain[1];
  const double cr_Ct15 = -cr_Ct13 + cr_Ct14;
  const double cr_Ct16 = cr_Ct7*(cr_Ct12 + cr_Ct15);
  const double cr_Ct17 = 0.5*cr_Ct10;
  const double cr_Ct18 = cr_Ct7*(-0.5*cr_Ct11 + cr_Ct15 + cr_Ct17);
  const double cr_Ct19 = -nu;
  const double cr_Ct20 = cr_Ct19 + 1.0;
  const double cr_Ct21 = cr_Ct20*r_strain[1];
  const double cr_Ct22 = cr_Ct20*r_strain[0];
  const double cr_Ct23 = cr_Ct10 + cr_Ct22;
  const double cr_Ct24 = 1.0/(1 - cr_Ct3);
  const double cr_Ct25 = cr_Ct24*cr_Ct6;
  const double cr_Ct26 = cr_Ct25*(-0.5*cr_Ct13 - 0.5*cr_Ct21 + cr_Ct23);
  const double cr_Ct27 = 0.22222222222222224*std::pow(cr_Ct26, 2.0);
  const double cr_Ct28 = 0.055555555555555552*std::pow(cr_Ct16, 2.0) + 0.22222222222222224*std::pow(cr_Ct18, 2.0) + cr_Ct27 + cr_Ct9;
  const double cr_Ct29 = 0.57735026918962584*threshold;
  const double cr_Ct30 = cr_Ct0*(cr_Ct1*(-1.0 + cr_Ct29/std::sqrt(cr_Ct28)) + 1.0);
  const double cr_Ct31 = cr_Ct10 - cr_Ct11;
  const double cr_Ct32 = 0.055555555555555552*std::pow(cr_Ct16, 1.0)*cr_Ct4;
  const double cr_Ct33 = 3.0*nu;
  const double cr_Ct34 = cr_Ct33 -1.0;
  const double cr_Ct35 = 0.11111111111111112*(cr_Ct18);
  const double cr_Ct36 = 0.11111111111111112*(cr_Ct26);
  const double cr_Ct37 = cr_Ct24*(cr_Ct33 - 2.0);
  const double cr_Ct38 = cr_Ct32 + cr_Ct34*cr_Ct35*cr_Ct4 + cr_Ct36*cr_Ct37;
  const double cr_Ct39 = cr_Ct1*cr_Ct29/std::pow(cr_Ct28, 3.0/2.0);
  const double cr_Ct40 = cr_Ct39*cr_Ct6;
  const double cr_Ct41 = cr_Ct13 + cr_Ct21;
  const double cr_Ct42 = cr_Ct25*(cr_Ct23 + cr_Ct41);
  const double cr_Ct43 = cr_Ct25*(-cr_Ct17 - 0.5*cr_Ct22 + cr_Ct41);
  const double cr_Ct44 = cr_Ct1*(cr_Ct29/std::sqrt(cr_Ct27 + 0.055555555555555552*std::pow(cr_Ct42, 2.0) + 0.22222222222222224*std::pow(cr_Ct43, 2.0) + std::pow(cr_Ct25*r_strain[2]*(cr_Ct19 + 0.5), 2.0)) - 1.0) + 1.0;
  const double cr_Ct45 = cr_Ct44*nu;
  const double cr_Ct46 = 0.055555555555555552*(cr_Ct42);
  const double cr_Ct47 = cr_Ct34*cr_Ct36;
  const double cr_Ct48 = 2.0 - cr_Ct33;
  const double cr_Ct49 = 0.11111111111111112*(cr_Ct43);
  const double cr_Ct50 = cr_Ct25*cr_Ct39;
  const double cr_Ct51 = cr_Ct39*cr_Ct9;
  const double cr_Ct52 = cr_Ct51*cr_Ct7/r_strain[2];
  const double cr_Ct53 = cr_Ct13 - cr_Ct14;
  const double cr_Ct54 = cr_Ct32 + cr_Ct35*cr_Ct37 + cr_Ct4*cr_Ct47;
  const double cr_Ct55 = std::pow(Young, 2)*cr_Ct2*cr_Ct39*cr_Ct4*r_strain[2]/std::pow(cr_Ct5, 2);
  r_Ct(0,0)=cr_Ct7*(cr_Ct30 - cr_Ct31*cr_Ct38*cr_Ct40);
  r_Ct(0,1)=-cr_Ct7*(cr_Ct12*cr_Ct50*(cr_Ct46 + cr_Ct47 + cr_Ct48*cr_Ct49) + cr_Ct45);
  r_Ct(0,2)=cr_Ct31*cr_Ct52;
  r_Ct(1,0)=-cr_Ct7*(cr_Ct15*cr_Ct50*(cr_Ct34*cr_Ct49 + cr_Ct36*cr_Ct48 + cr_Ct46) + cr_Ct45);
  r_Ct(1,1)=cr_Ct7*(cr_Ct30 - cr_Ct40*cr_Ct53*cr_Ct54);
  r_Ct(1,2)=cr_Ct52*cr_Ct53;
  r_Ct(2,0)=cr_Ct38*cr_Ct55;
  r_Ct(2,1)=cr_Ct54*cr_Ct55;
  r_Ct(2,2)=cr_Ct8*(cr_Ct44 - cr_Ct51);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  const auto &r_props = rValues.GetMaterialProperties();
  const double Young = r_props[YOUNG_MODULUS];
  const double nu = r_props[POISSON_RATIO];
  const double Gf = r_props[FRACTURE_ENERGY];
  const double phi = r_props[FRICTION_ANGLE] * Globals::Pi / 180.0;
  const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
  auto &r_Ct = rValues.GetConstitutiveMatrix();
  const auto &r_strain = rValues.GetStrainVector();

  const bool has_symmetric_yield_stress = r_props.Has(YIELD_STRESS);
  const double threshold_compression = has_symmetric_yield_stress ? r_props[YIELD_STRESS] : r_props[YIELD_STRESS_COMPRESSION];
  const double threshold_tension = has_symmetric_yield_stress ? r_props[YIELD_STRESS] : r_props[YIELD_STRESS_TENSION];

  const double cr_Ct0 = nu - 1.0;
  const double cr_Ct1 = fabs(threshold_compression/threshold_tension);
  const double cr_Ct2 = 1.0/(1.0 - 0.5*characteristic_length*std::pow(threshold_compression, 2)/(Gf*Young*std::pow(cr_Ct1, 2)));
  const double cr_Ct3 = nu*r_strain[1];
  const double cr_Ct4 = cr_Ct0*r_strain[0];
  const double cr_Ct5 = -cr_Ct3 + cr_Ct4;
  const double cr_Ct6 = nu*r_strain[0];
  const double cr_Ct7 = cr_Ct0*r_strain[1];
  const double cr_Ct8 = -cr_Ct6 + cr_Ct7;
  const double cr_Ct9 = 2*nu;
  const double cr_Ct10 = 1.0/(cr_Ct9 - 1);
  const double cr_Ct11 = nu + 1.0;
  const double cr_Ct12 = Young/cr_Ct11;
  const double cr_Ct13 = cr_Ct10*cr_Ct12;
  const double cr_Ct14 = cr_Ct13*(cr_Ct5 + cr_Ct8);
  const double cr_Ct15 = std::tan(0.5*phi + 0.78539816339744828);
  const double cr_Ct16 = cr_Ct1/std::pow(cr_Ct15, 2);
  const double cr_Ct17 = std::sin(phi);
  const double cr_Ct18 = cr_Ct16 + 1;
  const double cr_Ct19 = cr_Ct17*cr_Ct18;
  const double cr_Ct20 = 0.16666666666666666*cr_Ct16 + 0.16666666666666666*cr_Ct19 - 0.16666666666666666;
  const double cr_Ct21 = cr_Ct14*cr_Ct20;
  const double cr_Ct22 = nu - 0.5;
  const double cr_Ct23 = cr_Ct13*cr_Ct22;
  const double cr_Ct24 = std::pow(cr_Ct23*r_strain[2], 2.0);
  const double cr_Ct25 = std::pow(cr_Ct14, 2.0);
  const double cr_Ct26 = 0.5*cr_Ct3;
  const double cr_Ct27 = cr_Ct13*(cr_Ct26 - 0.5*cr_Ct4 + cr_Ct8);
  const double cr_Ct28 = std::pow(cr_Ct27, 2.0);
  const double cr_Ct29 = -nu;
  const double cr_Ct30 = cr_Ct29 + 1.0;
  const double cr_Ct31 = cr_Ct30*r_strain[1];
  const double cr_Ct32 = cr_Ct30*r_strain[0];
  const double cr_Ct33 = cr_Ct3 + cr_Ct32;
  const double cr_Ct34 = 1 - cr_Ct9;
  const double cr_Ct35 = 1.0/cr_Ct34;
  const double cr_Ct36 = cr_Ct12*cr_Ct35;
  const double cr_Ct37 = cr_Ct36*(-0.5*cr_Ct31 + cr_Ct33 - 0.5*cr_Ct6);
  const double cr_Ct38 = std::pow(cr_Ct37, 2.0);
  const double cr_Ct39 = 0.22222222222222227*cr_Ct38;
  const double cr_Ct40 = cr_Ct24 + 0.055555555555555552*cr_Ct25 + 0.22222222222222227*cr_Ct28 + cr_Ct39;
  const double cr_Ct41 = std::sqrt(cr_Ct40);
  const double cr_Ct42 = -cr_Ct16;
  const double cr_Ct43 = cr_Ct42 + 1;
  const double cr_Ct44 = -cr_Ct17*cr_Ct43 + cr_Ct18;
  const double cr_Ct45 = cr_Ct29 + 0.5;
  const double cr_Ct46 = std::pow(cr_Ct36*cr_Ct45*r_strain[2], 2.0);
  const double cr_Ct47 = cr_Ct31 + cr_Ct6;
  const double cr_Ct48 = cr_Ct36*(cr_Ct33 + cr_Ct47);
  const double cr_Ct49 = std::pow(cr_Ct48, 2.0);
  const double cr_Ct50 = 0.44444444444444453*cr_Ct38;
  const double cr_Ct51 = cr_Ct36*(-cr_Ct26 - 0.5*cr_Ct32 + cr_Ct47);
  const double cr_Ct52 = std::pow(cr_Ct51, 2.0);
  const double cr_Ct53 = 1.0/(2.0*cr_Ct46 + 0.1111111111111111*cr_Ct49 + cr_Ct50 + 0.44444444444444453*cr_Ct52);
  const double cr_Ct54 = cr_Ct39 + cr_Ct46 + 0.055555555555555552*cr_Ct49 + 0.22222222222222227*cr_Ct52;
  const double cr_Ct55 = std::sqrt(cr_Ct54);
  const double cr_Ct56 = 1.0/cr_Ct55;
  const double cr_Ct57 = std::pow(r_strain[2], 2);
  const double cr_Ct58 = std::pow(cr_Ct45, 2);
  const double cr_Ct59 = 0.33333333333333331*cr_Ct6;
  const double cr_Ct60 = 0.66666666666666674*cr_Ct3;
  const double cr_Ct61 = 0.66666666666666674*cr_Ct32;
  const double cr_Ct62 = 0.33333333333333331*cr_Ct31;
  const double cr_Ct63 = -cr_Ct59 + cr_Ct60 + cr_Ct61 - cr_Ct62;
  const double cr_Ct64 = 0.66666666666666674*cr_Ct6;
  const double cr_Ct65 = 0.33333333333333331*cr_Ct3;
  const double cr_Ct66 = cr_Ct57*cr_Ct58 - cr_Ct63*(0.66666666666666674*cr_Ct31 - 0.33333333333333331*cr_Ct32 + cr_Ct64 - cr_Ct65);
  const double cr_Ct67 = 5.196152422706632*cr_Ct66;
  const double cr_Ct68 = std::pow(Young, 2)/std::pow(cr_Ct11, 2);
  const double cr_Ct69 = cr_Ct68/std::pow(cr_Ct34, 2);
  const double cr_Ct70 = 0.33333333333333331*std::asin(cr_Ct53*cr_Ct56*cr_Ct67*cr_Ct69);
  const double cr_Ct71 = std::cos(cr_Ct70);
  const double cr_Ct72 = 0.5*cr_Ct44*cr_Ct71;
  const double cr_Ct73 = std::sin(cr_Ct70);
  const double cr_Ct74 = cr_Ct43/cr_Ct17;
  const double cr_Ct75 = cr_Ct17*(cr_Ct18 - cr_Ct74);
  const double cr_Ct76 = 0.28867513459481292*cr_Ct73*cr_Ct75;
  const double cr_Ct77 = -cr_Ct72 + cr_Ct76;
  const double cr_Ct78 = cr_Ct41*cr_Ct77;
  const double cr_Ct79 = cr_Ct21 - cr_Ct78;
  const double cr_Ct80 = 0.5*threshold_compression*std::cos(phi)/cr_Ct15;
  const double cr_Ct81 = cr_Ct0*(cr_Ct2*(-1.0 + cr_Ct80/cr_Ct79) + 1.0);
  const double cr_Ct82 = cr_Ct3 - cr_Ct4;
  const double cr_Ct83 = 0.16666666666666666*cr_Ct16 + 0.16666666666666666*cr_Ct19 - 0.16666666666666666;
  const double cr_Ct84 = -cr_Ct10*cr_Ct83;
  const double cr_Ct85 = 0.055555555555555552*cr_Ct10*std::pow(cr_Ct14, 1.0);
  const double cr_Ct86 = 3.0*nu;
  const double cr_Ct87 = cr_Ct86 - 1.0;
  const double cr_Ct88 = std::pow(cr_Ct27, 1.0);
  const double cr_Ct89 = std::pow(cr_Ct37, 1.0);
  const double cr_Ct90 = 0.11111111111111113*cr_Ct35*(cr_Ct86 - 2.0);
  const double cr_Ct91 = 1.0/cr_Ct41;
  const double cr_Ct92 = cr_Ct77*cr_Ct91;
  const double cr_Ct93 = cr_Ct10*cr_Ct68;
  const double cr_Ct94 = cr_Ct59 - cr_Ct60 - cr_Ct61 + cr_Ct62;
  const double cr_Ct95 = cr_Ct10*std::pow(cr_Ct22, 2)*cr_Ct57 + cr_Ct35*cr_Ct94*(-0.33333333333333331*cr_Ct4 - cr_Ct64 + cr_Ct65 + 0.66666666666666674*cr_Ct7);
  const double cr_Ct96 = 1.0/(2.0*cr_Ct24 + 0.1111111111111111*cr_Ct25 + 0.44444444444444453*cr_Ct28 + cr_Ct50);
  const double cr_Ct97 = 5.196152422706632*cr_Ct96;
  const double cr_Ct98 = cr_Ct95*cr_Ct97;
  const double cr_Ct99 = 0.33333333333333331*std::asin(cr_Ct91*cr_Ct93*cr_Ct98);
  const double cr_Ct100 = 0.064150029909958411*cr_Ct44;
  const double cr_Ct101 = cr_Ct100*std::sin(cr_Ct99) + 0.037037037037037035*cr_Ct75*std::cos(cr_Ct99);
  const double cr_Ct102 = nu;
  const double cr_Ct103 = 3.4641016151377553*cr_Ct6;
  const double cr_Ct104 = 1.7320508075688772*cr_Ct3;
  const double cr_Ct105 = 1.7320508075688772*cr_Ct32;
  const double cr_Ct106 = 3.4641016151377553*cr_Ct31;
  const double cr_Ct107 = cr_Ct103 - cr_Ct104 - cr_Ct105 + cr_Ct106;
  const double cr_Ct108 = 5.196152422706632*nu;
  const double cr_Ct109 = std::pow(cr_Ct40, -2);
  const double cr_Ct110 = std::pow(cr_Ct48, 1.0);
  const double cr_Ct111 = 0.22222222222222221*cr_Ct110;
  const double cr_Ct112 = 2.0 - cr_Ct86;
  const double cr_Ct113 = cr_Ct112*cr_Ct89;
  const double cr_Ct114 = std::pow(cr_Ct51, 1.0);
  const double cr_Ct115 = cr_Ct114*cr_Ct87;
  const double cr_Ct116 = 1.299038105676658*cr_Ct111 + 0.57735026918962584*cr_Ct113 + 0.57735026918962584*cr_Ct115;
  const double cr_Ct117 = 0.055555555555555552*cr_Ct110;
  const double cr_Ct118 = 0.11111111111111113*cr_Ct113 + 0.11111111111111113*cr_Ct115 + cr_Ct117;
  const double cr_Ct119 = 1.0/cr_Ct40;
  const double cr_Ct120 = std::pow(-std::pow(Young, 4)*std::pow(cr_Ct66, 2)/(std::pow(cr_Ct11, 4)*std::pow(cr_Ct34, 4)*std::pow(cr_Ct54, 3)) + 0.14814814814814814, -1.0/2.0);
  const double cr_Ct121 = cr_Ct120*cr_Ct36;
  const double cr_Ct122 = -cr_Ct101*cr_Ct121*(-cr_Ct109*cr_Ct116*cr_Ct13*cr_Ct95 - cr_Ct118*cr_Ct119*cr_Ct13*cr_Ct98 + cr_Ct35*cr_Ct96*(cr_Ct107*(cr_Ct102 - 0.66666666666666674) + cr_Ct94*(cr_Ct108 - 1.7320508075688772))) + cr_Ct84 + cr_Ct92*(0.11111111111111113*cr_Ct10*cr_Ct87*cr_Ct88 + cr_Ct85 + cr_Ct89*cr_Ct90);
  const double cr_Ct123 = cr_Ct2*cr_Ct80;
  const double cr_Ct124 = cr_Ct123/std::pow(-cr_Ct21 + cr_Ct78, 2);
  const double cr_Ct125 = cr_Ct12*cr_Ct124;
  const double cr_Ct126 = cr_Ct72 - cr_Ct76;
  const double cr_Ct127 = cr_Ct2*(cr_Ct80/(cr_Ct126*cr_Ct55 + cr_Ct20*cr_Ct48) - 1.0) + 1.0;
  const double cr_Ct128 = cr_Ct127*nu;
  const double cr_Ct129 = cr_Ct87*cr_Ct89;
  const double cr_Ct130 = 0.11111111111111113*cr_Ct129;
  const double cr_Ct131 = cr_Ct112*cr_Ct114;
  const double cr_Ct132 = cr_Ct117 + cr_Ct130 + 0.11111111111111113*cr_Ct131;
  const double cr_Ct133 = cr_Ct126*cr_Ct56;
  const double cr_Ct134 = cr_Ct102 - 0.33333333333333331;
  const double cr_Ct135 = -cr_Ct103 + cr_Ct104 + cr_Ct105 - cr_Ct106;
  const double cr_Ct136 = cr_Ct108 - 3.4641016151377553;
  const double cr_Ct137 = 1.299038105676658*cr_Ct111 + 0.57735026918962584*cr_Ct129 + 0.57735026918962584*cr_Ct131;
  const double cr_Ct138 = std::pow(cr_Ct54, -2);
  const double cr_Ct139 = cr_Ct138*cr_Ct36*cr_Ct66;
  const double cr_Ct140 = cr_Ct53/cr_Ct54;
  const double cr_Ct141 = cr_Ct140*cr_Ct36*cr_Ct67;
  const double cr_Ct142 = -cr_Ct100*cr_Ct73 + 0.037037037037037035*cr_Ct17*cr_Ct71*(cr_Ct42 + cr_Ct74 - 1);
  const double cr_Ct143 = cr_Ct121*cr_Ct142;
  const double cr_Ct144 = cr_Ct123/std::pow(cr_Ct79, 2);
  const double cr_Ct145 = cr_Ct144*cr_Ct36;
  const double cr_Ct146 = 1.0/r_strain[2];
  const double cr_Ct147 = cr_Ct146;
  const double cr_Ct148 = 10.392304845413264*cr_Ct53*cr_Ct58*r_strain[2];
  const double cr_Ct149 = cr_Ct146*cr_Ct46*cr_Ct67;
  const double cr_Ct150 = cr_Ct138*cr_Ct149;
  const double cr_Ct151 = cr_Ct140*cr_Ct149;
  const double cr_Ct152 = cr_Ct120*cr_Ct69;
  const double cr_Ct153 = cr_Ct124*cr_Ct13*(cr_Ct101*cr_Ct152*cr_Ct41*cr_Ct56*(-cr_Ct148 + cr_Ct150 + cr_Ct151) - cr_Ct147*cr_Ct24*cr_Ct92);
  const double cr_Ct154 = cr_Ct6 - cr_Ct7;
  const double cr_Ct155 = cr_Ct12*cr_Ct95;
  const double cr_Ct156 = -cr_Ct101*cr_Ct120*cr_Ct13*cr_Ct35*(-cr_Ct109*cr_Ct137*cr_Ct155 - cr_Ct119*cr_Ct132*cr_Ct155*cr_Ct97 + cr_Ct96*(cr_Ct107*cr_Ct134 + cr_Ct136*cr_Ct94)) + cr_Ct84 + cr_Ct92*(cr_Ct10*cr_Ct130 + cr_Ct85 + cr_Ct88*cr_Ct90);
  const double cr_Ct157 = cr_Ct124*cr_Ct22*cr_Ct93*r_strain[2];
  r_Ct(0,0)=cr_Ct13*(cr_Ct122*cr_Ct125*cr_Ct82 + cr_Ct81);
  r_Ct(0,1)=-cr_Ct13*(cr_Ct128 + cr_Ct145*cr_Ct5*(cr_Ct132*cr_Ct133 + cr_Ct143*(-cr_Ct132*cr_Ct141 - cr_Ct137*cr_Ct139 + cr_Ct53*(cr_Ct134*cr_Ct135 + cr_Ct136*cr_Ct63)) + cr_Ct83));
  r_Ct(0,2)=cr_Ct153*cr_Ct82;
  r_Ct(1,0)=-cr_Ct13*(cr_Ct128 + cr_Ct145*cr_Ct8*(cr_Ct118*cr_Ct133 + cr_Ct143*(-cr_Ct116*cr_Ct139 - cr_Ct118*cr_Ct141 + cr_Ct53*(cr_Ct135*(0.66666666666666674 - cr_Ct102) + cr_Ct63*(1.7320508075688772 - cr_Ct108))) + cr_Ct83));
  r_Ct(1,1)=cr_Ct13*(cr_Ct125*cr_Ct154*cr_Ct156 + cr_Ct81);
  r_Ct(1,2)=cr_Ct153*cr_Ct154;
  r_Ct(2,0)=-cr_Ct122*cr_Ct157;
  r_Ct(2,1)=-cr_Ct156*cr_Ct157;
  r_Ct(2,2)=cr_Ct23*(cr_Ct127 - cr_Ct144*r_strain[2]*(cr_Ct133*cr_Ct147*cr_Ct46 + cr_Ct142*cr_Ct152*(cr_Ct148 - cr_Ct150 - cr_Ct151)));
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<TrescaYieldSurface<TrescaPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  const auto &r_props = rValues.GetMaterialProperties();
  const double Young = r_props[YOUNG_MODULUS];
  const double nu = r_props[POISSON_RATIO];
  const double Gf = r_props[FRACTURE_ENERGY];
  const double phi = r_props[FRICTION_ANGLE];
  const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
  const double threshold = r_props[YIELD_STRESS];
  auto &r_Ct = rValues.GetConstitutiveMatrix();
  const auto &r_strain = rValues.GetStrainVector();
  const double sin_phi = std::sin(phi * Globals::Pi / 180.0);

  const double cr_Ct0 = nu - 1.0;
  const double cr_Ct1 = 1.0/(1.0 - 0.5*characteristic_length*std::pow(threshold, 2)/(Gf*Young));
  const double cr_Ct2 = nu*r_strain[1];
  const double cr_Ct3 = cr_Ct0*r_strain[0];
  const double cr_Ct4 = -cr_Ct2 + cr_Ct3;
  const double cr_Ct5 = nu*r_strain[0];
  const double cr_Ct6 = cr_Ct0*r_strain[1];
  const double cr_Ct7 = -cr_Ct5 + cr_Ct6;
  const double cr_Ct8 = 2*nu;
  const double cr_Ct9 = 1.0/(cr_Ct8 - 1);
  const double cr_Ct10 = nu + 1.0;
  const double cr_Ct11 = Young/cr_Ct10;
  const double cr_Ct12 = cr_Ct11*cr_Ct9;
  const double cr_Ct13 = cr_Ct12*(cr_Ct4 + cr_Ct7);
  const double cr_Ct14 = 1.7320508075688772*sin_phi;
  const double cr_Ct15 = 1.0/(cr_Ct14 - 5.196152422706632);
  const double cr_Ct16 = 2.0*sin_phi;
  const double cr_Ct17 = cr_Ct15*cr_Ct16;
  const double cr_Ct18 = cr_Ct13*cr_Ct17;
  const double cr_Ct19 = nu - 0.5;
  const double cr_Ct20 = cr_Ct12*cr_Ct19;
  const double cr_Ct21 = std::pow(cr_Ct20*r_strain[2], 2.0);
  const double cr_Ct22 = 0.5*cr_Ct2;
  const double cr_Ct23 = cr_Ct12*(cr_Ct22 - 0.5*cr_Ct3 + cr_Ct7);
  const double cr_Ct24 = -nu;
  const double cr_Ct25 = cr_Ct24 + 1.0;
  const double cr_Ct26 = cr_Ct25*r_strain[1];
  const double cr_Ct27 = cr_Ct25*r_strain[0];
  const double cr_Ct28 = cr_Ct2 + cr_Ct27;
  const double cr_Ct29 = 1.0/(1 - cr_Ct8);
  const double cr_Ct30 = cr_Ct11*cr_Ct29;
  const double cr_Ct31 = cr_Ct30*(-0.5*cr_Ct26 + cr_Ct28 - 0.5*cr_Ct5);
  const double cr_Ct32 = 0.22222222222222227*std::pow(cr_Ct31, 2.0);
  const double cr_Ct33 = std::sqrt(0.055555555555555552*std::pow(cr_Ct13, 2.0) + cr_Ct21 + 0.22222222222222227*std::pow(cr_Ct23, 2.0) + cr_Ct32);
  const double cr_Ct34 = sin_phi - 1;
  const double cr_Ct35 = cr_Ct15*cr_Ct34*fabs(threshold*(sin_phi + 3.0)/cr_Ct34);
  const double cr_Ct36 = cr_Ct0*(-cr_Ct1*(cr_Ct35/(cr_Ct18 - cr_Ct33) + 1) + 1);
  const double cr_Ct37 = cr_Ct2 - cr_Ct3;
  const double cr_Ct38 = cr_Ct17*cr_Ct9;
  const double cr_Ct39 = 1.0/cr_Ct33;
  const double cr_Ct40 = 0.055555555555555552*std::pow(cr_Ct13, 1.0)*cr_Ct9;
  const double cr_Ct41 = 3.0*nu;
  const double cr_Ct42 = cr_Ct41 - 1.0;
  const double cr_Ct43 = 0.11111111111111113*std::pow(cr_Ct23, 1.0);
  const double cr_Ct44 = 0.11111111111111113*std::pow(cr_Ct31, 1.0);
  const double cr_Ct45 = cr_Ct29*(cr_Ct41 - 2.0);
  const double cr_Ct46 = cr_Ct38 - cr_Ct39*(cr_Ct40 + cr_Ct42*cr_Ct43*cr_Ct9 + cr_Ct44*cr_Ct45);
  const double cr_Ct47 = cr_Ct13*cr_Ct15*sin_phi;
  const double cr_Ct48 = 0.5*cr_Ct33;
  const double cr_Ct49 = 0.25*cr_Ct1*cr_Ct35;
  const double cr_Ct50 = cr_Ct49/std::pow(cr_Ct47 - cr_Ct48, 2);
  const double cr_Ct51 = cr_Ct11*cr_Ct50;
  const double cr_Ct52 = cr_Ct1*(-cr_Ct35/(-cr_Ct18 + cr_Ct33) + 1);
  const double cr_Ct53 = nu*(1 - cr_Ct52);
  const double cr_Ct54 = cr_Ct16/(5.196152422706632 - cr_Ct14);
  const double cr_Ct55 = cr_Ct26 + cr_Ct5;
  const double cr_Ct56 = cr_Ct30*(cr_Ct28 + cr_Ct55);
  const double cr_Ct57 = 0.055555555555555552*std::pow(cr_Ct56, 1.0);
  const double cr_Ct58 = cr_Ct42*cr_Ct44;
  const double cr_Ct59 = 2.0 - cr_Ct41;
  const double cr_Ct60 = cr_Ct30*(-cr_Ct22 - 0.5*cr_Ct27 + cr_Ct55);
  const double cr_Ct61 = 0.11111111111111113*std::pow(cr_Ct60, 1.0);
  const double cr_Ct62 = std::pow(cr_Ct32 + 0.055555555555555552*std::pow(cr_Ct56, 2.0) + 0.22222222222222227*std::pow(cr_Ct60, 2.0) + std::pow(cr_Ct30*r_strain[2]*(cr_Ct24 + 0.5), 2.0), -1.0/2.0);
  const double cr_Ct63 = cr_Ct49/std::pow(-cr_Ct47 + cr_Ct48, 2);
  const double cr_Ct64 = cr_Ct30*cr_Ct63;
  const double cr_Ct65 = cr_Ct21*cr_Ct39;
  const double cr_Ct66 = cr_Ct12*cr_Ct50*cr_Ct65/r_strain[2];
  const double cr_Ct67 = cr_Ct5 - cr_Ct6;
  const double cr_Ct68 = cr_Ct38 - cr_Ct39*(cr_Ct40 + cr_Ct43*cr_Ct45 + cr_Ct58*cr_Ct9);
  const double cr_Ct69 = std::pow(Young, 2)*cr_Ct19*cr_Ct50*cr_Ct9*r_strain[2]/std::pow(cr_Ct10, 2);
  r_Ct(0,0)=cr_Ct12*(cr_Ct36 + cr_Ct37*cr_Ct46*cr_Ct51);
  r_Ct(0,1)=-cr_Ct12*(cr_Ct4*cr_Ct64*(cr_Ct54 + cr_Ct62*(cr_Ct57 + cr_Ct58 + cr_Ct59*cr_Ct61)) + cr_Ct53);
  r_Ct(0,2)=cr_Ct37*cr_Ct66;
  r_Ct(1,0)=-cr_Ct12*(cr_Ct53 + cr_Ct64*cr_Ct7*(cr_Ct54 + cr_Ct62*(cr_Ct42*cr_Ct61 + cr_Ct44*cr_Ct59 + cr_Ct57)));
  r_Ct(1,1)=cr_Ct12*(cr_Ct36 + cr_Ct51*cr_Ct67*cr_Ct68);
  r_Ct(1,2)=cr_Ct66*cr_Ct67;
  r_Ct(2,0)=-cr_Ct46*cr_Ct69;
  r_Ct(2,1)=-cr_Ct68*cr_Ct69;
  r_Ct(2,2)=cr_Ct20*(-cr_Ct52 - cr_Ct63*cr_Ct65 + 1.0);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<RankineYieldSurface<RankinePlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  const auto &r_props = rValues.GetMaterialProperties();
  const double Young = r_props[YOUNG_MODULUS];
  const double nu = r_props[POISSON_RATIO];
  const double Gf = r_props[FRACTURE_ENERGY];
  const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
  const double threshold = r_props[YIELD_STRESS];
  auto &r_Ct = rValues.GetConstitutiveMatrix();
  const auto &r_strain = rValues.GetStrainVector();

	const double cr_Ct0 = nu - 1.0;
	const double cr_Ct1 = 1.0/(1.0 - 0.5*characteristic_length*std::pow(threshold, 2)/(Gf*Young));
	const double cr_Ct2 = 2*nu;
	const double cr_Ct3 = cr_Ct2 - 1;
	const double cr_Ct4 = std::pow(cr_Ct3, -2);
	const double cr_Ct5 = std::pow(r_strain[2], 2);
	const double cr_Ct6 = nu - 0.5;
	const double cr_Ct7 = std::pow(cr_Ct6, 2);
	const double cr_Ct8 = cr_Ct5*cr_Ct7;
	const double cr_Ct9 = cr_Ct4*cr_Ct8;
	const double cr_Ct10 = nu*r_strain[0];
	const double cr_Ct11 = nu*r_strain[1];
	const double cr_Ct12 = -cr_Ct11;
	const double cr_Ct13 = -nu;
	const double cr_Ct14 = cr_Ct13 + 1.0;
	const double cr_Ct15 = cr_Ct14*r_strain[1];
	const double cr_Ct16 = cr_Ct14*r_strain[0];
	const double cr_Ct17 = cr_Ct10 + cr_Ct12 + cr_Ct15 - cr_Ct16;
	const double cr_Ct18 = -cr_Ct2;
	const double cr_Ct19 = cr_Ct18 + 1;
	const double cr_Ct20 = std::pow(cr_Ct19, -2);
	const double cr_Ct21 = 0.25*cr_Ct20;
	const double cr_Ct22 = nu + 1.0;
	const double cr_Ct23 = std::pow(Young, 2)/std::pow(cr_Ct22, 2);
	const double cr_Ct24 = std::sqrt(cr_Ct23*(std::pow(cr_Ct17, 2)*cr_Ct21 + cr_Ct9));
	const double cr_Ct25 = cr_Ct0*r_strain[0];
	const double cr_Ct26 = 1.0/cr_Ct3;
	const double cr_Ct27 = Young/cr_Ct22;
	const double cr_Ct28 = cr_Ct26*cr_Ct27;
	const double cr_Ct29 = 0.5*cr_Ct28;
	const double cr_Ct30 = cr_Ct0*r_strain[1];
	const double cr_Ct31 = -cr_Ct10*cr_Ct29 - cr_Ct11*cr_Ct29 + cr_Ct25*cr_Ct29 + cr_Ct29*cr_Ct30;
	const double cr_Ct32 = cr_Ct24 + cr_Ct31;
	const double cr_Ct33 = cr_Ct0*(cr_Ct1*(-1.0 + threshold/cr_Ct32) + 1.0);
	const double cr_Ct34 = cr_Ct11 - cr_Ct25;
	const double cr_Ct35 = 0.5*nu;
	const double cr_Ct36 = 1.0/cr_Ct24;
	const double cr_Ct37 = cr_Ct2 - 1.0;
	const double cr_Ct38 = cr_Ct17*cr_Ct36*cr_Ct37;
	const double cr_Ct39 = -cr_Ct10;
	const double cr_Ct40 = cr_Ct30 + cr_Ct39;
	const double cr_Ct41 = cr_Ct23*cr_Ct4;
	const double cr_Ct42 = std::sqrt(cr_Ct41*(cr_Ct8 + 0.25*std::pow(cr_Ct34 + cr_Ct40, 2)));
	const double cr_Ct43 = cr_Ct1*threshold;
	const double cr_Ct44 = cr_Ct43/std::pow(cr_Ct31 + cr_Ct42, 2);
	const double cr_Ct45 = cr_Ct44*(cr_Ct21*cr_Ct27*cr_Ct38 - cr_Ct26*cr_Ct35 + cr_Ct26*(cr_Ct35 - 0.5));
	const double cr_Ct46 = cr_Ct27/cr_Ct19;
	const double cr_Ct47 = 0.5*cr_Ct46;
	const double cr_Ct48 = cr_Ct11 - cr_Ct15 + cr_Ct16 + cr_Ct39;
	const double cr_Ct49 = std::sqrt(cr_Ct20*cr_Ct23*(0.25*std::pow(cr_Ct48, 2) + cr_Ct5*std::pow(cr_Ct13 + 0.5, 2)));
	const double cr_Ct50 = cr_Ct1*(threshold/(cr_Ct10*cr_Ct47 + cr_Ct11*cr_Ct47 + cr_Ct15*cr_Ct47 + cr_Ct16*cr_Ct47 + cr_Ct49) - 1.0) + 1.0;
	const double cr_Ct51 = cr_Ct50*nu;
	const double cr_Ct52 = 0.25*cr_Ct46;
	const double cr_Ct53 = cr_Ct48*cr_Ct52/cr_Ct49;
	const double cr_Ct54 = cr_Ct43/std::pow(cr_Ct32, 2);
	const double cr_Ct55 = cr_Ct46*cr_Ct54;
	const double cr_Ct56 = std::pow(Young, 3)*cr_Ct44*cr_Ct7*r_strain[2]/(std::pow(cr_Ct22, 3)*std::pow(cr_Ct3, 3)*cr_Ct42);
	const double cr_Ct57 = cr_Ct10 - cr_Ct30;
	const double cr_Ct58 = cr_Ct44*(cr_Ct38*cr_Ct52 - 0.5);
	const double cr_Ct59 = cr_Ct6*r_strain[2];
	r_Ct(0,0)=cr_Ct28*(cr_Ct27*cr_Ct34*cr_Ct45 + cr_Ct33);
	r_Ct(0,1)=-cr_Ct28*(cr_Ct51 + cr_Ct55*(cr_Ct12 + cr_Ct25)*(cr_Ct37*cr_Ct53 + 0.5));
	r_Ct(0,2)=cr_Ct34*cr_Ct56;
	r_Ct(1,0)=-cr_Ct28*(cr_Ct40*cr_Ct55*(cr_Ct53*(cr_Ct18 + 1.0) + 0.5) + cr_Ct51);
	r_Ct(1,1)=cr_Ct28*(cr_Ct28*cr_Ct57*cr_Ct58 + cr_Ct33);
	r_Ct(1,2)=cr_Ct56*cr_Ct57;
	r_Ct(2,0)=-cr_Ct23*cr_Ct26*cr_Ct45*cr_Ct59;
	r_Ct(2,1)=-cr_Ct41*cr_Ct58*cr_Ct59;
	r_Ct(2,2)=cr_Ct28*cr_Ct6*(-cr_Ct23*cr_Ct36*cr_Ct54*cr_Ct9 + cr_Ct50);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<SimoJuYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalTrescaYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalTrescaYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalTrescaYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalTrescaYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalDruckerPragerYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalDruckerPragerYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalDruckerPragerYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalDruckerPragerYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalRankineYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalRankineYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalRankineYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalRankineYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalSimoJuYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalSimoJuYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalSimoJuYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalSimoJuYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<ThermalMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

// 3D exponential
template class AutomaticDifferentiationTangentUtilities<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<ThermalModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<ThermalTrescaYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<ThermalDruckerPragerYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<ThermalRankineYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<ThermalSimoJuYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<ThermalMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>;

template class AutomaticDifferentiationTangentUtilities<VonMisesYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<TrescaYieldSurface<TrescaPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<RankineYieldSurface<RankinePlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<SimoJuYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>;

// 3D Linear
template class AutomaticDifferentiationTangentUtilities<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<ThermalModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<ThermalTrescaYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<ThermalDruckerPragerYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<ThermalRankineYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<ThermalSimoJuYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<ThermalMohrCoulombYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>;

template class AutomaticDifferentiationTangentUtilities<VonMisesYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<TrescaYieldSurface<TrescaPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<RankineYieldSurface<RankinePlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<SimoJuYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>;

// 2D exponential
template class AutomaticDifferentiationTangentUtilities<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<ThermalModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<ThermalTrescaYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<ThermalDruckerPragerYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<ThermalRankineYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<ThermalSimoJuYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<ThermalMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>;

template class AutomaticDifferentiationTangentUtilities<VonMisesYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<TrescaYieldSurface<TrescaPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<RankineYieldSurface<RankinePlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<SimoJuYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>;
// 2D Linear
template class AutomaticDifferentiationTangentUtilities<ThermalVonMisesYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<ThermalModifiedMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<ThermalTrescaYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<ThermalDruckerPragerYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<ThermalRankineYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<ThermalSimoJuYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<ThermalMohrCoulombYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>;

template class AutomaticDifferentiationTangentUtilities<VonMisesYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<TrescaYieldSurface<TrescaPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<RankineYieldSurface<RankinePlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<SimoJuYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>;
} // namespace Kratos