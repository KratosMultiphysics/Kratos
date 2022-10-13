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
#include "custom_utilities/autom_differentiation_tangent_utilities.h"
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

namespace Kratos
{

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<VonMisesYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  const auto &r_props = rValues.GetMaterialProperties();
  const double Young = r_props[YOUNG_MODULUS];
  const double nu = r_props[POISSON_RATIO];
  const double Gf = r_props[FRACTURE_ENERGY];
  const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());
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
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<TrescaYieldSurface<TrescaPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  const auto &r_props = rValues.GetMaterialProperties();
  const double Young = r_props[YOUNG_MODULUS];
  const double nu = r_props[POISSON_RATIO];
  const double Gf = r_props[FRACTURE_ENERGY];
  const double phi = r_props[FRICTION_ANGLE];
  const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());
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
  const double cr_Ct56 = sqrt(cr_Ct55*(std::pow(cr_Ct26, 2)*cr_Ct49 + std::pow(cr_Ct34, 2)*cr_Ct49 + std::pow(cr_Ct47, 2)*cr_Ct48 + cr_Ct53));
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
  const double cr_Ct85 = cr_Ct82*std::exp(-1.0*cr_Ct9*(cr_Ct2*cr_Ct5*(cr_Ct58*cr_Ct66*(2.0*cr_Ct42 + 2.0*cr_Ct43 + 2.0*cr_Ct45 + cr_Ct79 + cr_Ct80 + cr_Ct81) + 0.5*cr_Ct77)/cr_Ct6 - 1));
  const double cr_Ct86 = cr_Ct58*cr_Ct85;
  const double cr_Ct87 = 1.0*cr_Ct12*cr_Ct86;
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
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  const auto &r_props = rValues.GetMaterialProperties();
  const double Young = r_props[YOUNG_MODULUS];
  const double nu = r_props[POISSON_RATIO];
  const double Gf = r_props[FRACTURE_ENERGY];
  const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());
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
  const double cr_Ct12 = pow(cr_Ct11, -2);
  const double cr_Ct13 = cr_Ct0*r_strain[0];
  const double cr_Ct14 = nu*r_strain[2];
  const double cr_Ct15 = -cr_Ct14;
  const double cr_Ct16 = cr_Ct13 + cr_Ct15;
  const double cr_Ct17 = cr_Ct16 + cr_Ct4;
  const double cr_Ct18 = cr_Ct0*r_strain[1];
  const double cr_Ct19 = cr_Ct15 + cr_Ct18 + cr_Ct6;
  const double cr_Ct20 = pow(r_strain[3], 2);
  const double cr_Ct21 = pow(r_strain[4], 2);
  const double cr_Ct22 = pow(r_strain[5], 2);
  const double cr_Ct23 = cr_Ct20 + cr_Ct21 + cr_Ct22;
  const double cr_Ct24 = cr_Ct12*pow(cr_Ct17, 2) + cr_Ct12*pow(cr_Ct19, 2) + cr_Ct12*pow(cr_Ct9, 2) + cr_Ct23;
  const double cr_Ct25 = nu + 1.0;
  const double cr_Ct26 = pow(Young, 2)/pow(cr_Ct25, 2);
  const double cr_Ct27 = cr_Ct24*cr_Ct26;
  const double cr_Ct28 = sqrt(cr_Ct27);
  const double cr_Ct29 = 1.0/cr_Ct28;
  const double cr_Ct30 = Young/cr_Ct25;
  const double cr_Ct31 = cr_Ct29*cr_Ct30;
  const double cr_Ct32 = 0.66666666666666663*cr_Ct31;
  const double cr_Ct33 = cr_Ct19*cr_Ct2;
  const double cr_Ct34 = cr_Ct17*cr_Ct2;
  const double cr_Ct35 = 1.0/cr_Ct24;
  const double cr_Ct36 = 1.0*cr_Ct35;
  const double cr_Ct37 = cr_Ct20*cr_Ct36;
  const double cr_Ct38 = 1.0*cr_Ct21;
  const double cr_Ct39 = cr_Ct35*cr_Ct38;
  const double cr_Ct40 = cr_Ct22*cr_Ct36;
  const double cr_Ct41 = cr_Ct12*cr_Ct9;
  const double cr_Ct42 = cr_Ct17*cr_Ct41;
  const double cr_Ct43 = 2*cr_Ct5;
  const double cr_Ct44 = 2*cr_Ct3;
  const double cr_Ct45 = -cr_Ct44;
  const double cr_Ct46 = 2*cr_Ct14;
  const double cr_Ct47 = cr_Ct13 + cr_Ct18 - cr_Ct43 + cr_Ct45 - cr_Ct46 + cr_Ct7;
  const double cr_Ct48 = cr_Ct12*cr_Ct35*pow(cr_Ct47, 2);
  const double cr_Ct49 = cr_Ct12*cr_Ct19*(cr_Ct16 + cr_Ct45 + cr_Ct8);
  const double cr_Ct50 = sqrt(-cr_Ct36*cr_Ct42 - cr_Ct36*cr_Ct49 + cr_Ct37 + cr_Ct39 + cr_Ct40 + 0.33333333333333331*cr_Ct48);
  const double cr_Ct51 = -nu;
  const double cr_Ct52 = cr_Ct51 + 0.5;
  const double cr_Ct53 = pow(cr_Ct52, -2);
  const double cr_Ct54 = cr_Ct51 + 1.0;
  const double cr_Ct55 = cr_Ct54*r_strain[2];
  const double cr_Ct56 = cr_Ct5 + cr_Ct55;
  const double cr_Ct57 = cr_Ct3 + cr_Ct56;
  const double cr_Ct58 = cr_Ct54*r_strain[1];
  const double cr_Ct59 = cr_Ct14 + cr_Ct5 + cr_Ct58;
  const double cr_Ct60 = cr_Ct54*r_strain[0];
  const double cr_Ct61 = cr_Ct14 + cr_Ct60;
  const double cr_Ct62 = cr_Ct3 + cr_Ct61;
  const double cr_Ct63 = cr_Ct23 + cr_Ct53*pow(cr_Ct57, 2) + cr_Ct53*pow(cr_Ct59, 2) + cr_Ct53*pow(cr_Ct62, 2);
  const double cr_Ct64 = 1.0/cr_Ct63;
  const double cr_Ct65 = 1.0*cr_Ct64;
  const double cr_Ct66 = cr_Ct20*cr_Ct65;
  const double cr_Ct67 = cr_Ct38*cr_Ct64;
  const double cr_Ct68 = cr_Ct22*cr_Ct65;
  const double cr_Ct69 = cr_Ct53*cr_Ct65;
  const double cr_Ct70 = cr_Ct57*cr_Ct62;
  const double cr_Ct71 = cr_Ct69*cr_Ct70;
  const double cr_Ct72 = cr_Ct43 + cr_Ct44 + cr_Ct46 + cr_Ct55 + cr_Ct58 + cr_Ct60;
  const double cr_Ct73 = pow(cr_Ct72, 2);
  const double cr_Ct74 = cr_Ct53*cr_Ct64;
  const double cr_Ct75 = cr_Ct73*cr_Ct74;
  const double cr_Ct76 = 0.33333333333333331*cr_Ct75;
  const double cr_Ct77 = cr_Ct44 + cr_Ct56 + cr_Ct61;
  const double cr_Ct78 = cr_Ct59*cr_Ct77;
  const double cr_Ct79 = cr_Ct69*cr_Ct78;
  const double cr_Ct80 = cr_Ct66 + cr_Ct67 + cr_Ct68 - cr_Ct71 + cr_Ct76 - cr_Ct79;
  const double cr_Ct81 = pow(cr_Ct80, 3);
  const double cr_Ct82 = pow(cr_Ct81, -1.0/2.0);
  const double cr_Ct83 = r_strain[4]*r_strain[5];
  const double cr_Ct84 = cr_Ct83*r_strain[3];
  const double cr_Ct85 = 1.0*cr_Ct53;
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
  const double cr_Ct98 = 0.037037037037037035*cr_Ct97 + 0.037037037037037035*cr_Ct35*pow(cr_Ct47, 3)/pow(cr_Ct11, 3);
  const double cr_Ct99 = -cr_Ct10*cr_Ct37 - cr_Ct33*cr_Ct40 + cr_Ct36*cr_Ct84 - cr_Ct65*cr_Ct88 + cr_Ct98;
  const double cr_Ct100 = 0.33333333333333331*acos(5.196152422706632*cr_Ct31*cr_Ct82*cr_Ct99);
  const double cr_Ct101 = cos(cr_Ct100);
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
  const double cr_Ct121 = pow(cr_Ct63, -2);
  const double cr_Ct122 = cr_Ct121*cr_Ct85;
  const double cr_Ct123 = cr_Ct111*cr_Ct122;
  const double cr_Ct124 = cr_Ct121*cr_Ct38;
  const double cr_Ct125 = pow(cr_Ct52, -4);
  const double cr_Ct126 = cr_Ct111*cr_Ct125;
  const double cr_Ct127 = cr_Ct121*cr_Ct126;
  const double cr_Ct128 = 1.0*cr_Ct70;
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
  const double cr_Ct140 = sqrt(cr_Ct80);
  const double cr_Ct141 = 1.0/cr_Ct140;
  const double cr_Ct142 = 1.0/(1.0 - cr_Ct1);
  const double cr_Ct143 = -cr_Ct38 + cr_Ct87;
  const double cr_Ct144 = cr_Ct143*cr_Ct62;
  const double cr_Ct145 = cr_Ct142*cr_Ct144;
  const double cr_Ct146 = cr_Ct142*cr_Ct57;
  const double cr_Ct147 = cr_Ct142*cr_Ct59;
  const double cr_Ct148 = pow(cr_Ct72, 3);
  const double cr_Ct149 = cr_Ct142*(-cr_Ct90 - cr_Ct91 - cr_Ct92 + cr_Ct94 + cr_Ct95);
  const double cr_Ct150 = cr_Ct149*cr_Ct72;
  const double cr_Ct151 = 0.037037037037037035*cr_Ct148*cr_Ct64/pow(cr_Ct52, 3) - 0.037037037037037035*cr_Ct150;
  const double cr_Ct152 = cr_Ct145*cr_Ct65 - cr_Ct146*cr_Ct66 - cr_Ct147*cr_Ct68 + cr_Ct151 + cr_Ct65*cr_Ct84;
  const double cr_Ct153 = -cr_Ct66 - cr_Ct67 - cr_Ct68 + cr_Ct71 - cr_Ct76 + cr_Ct79;
  const double cr_Ct154 = pow(cr_Ct153, 3);
  const double cr_Ct155 = cr_Ct26*cr_Ct63;
  const double cr_Ct156 = sqrt(cr_Ct155);
  const double cr_Ct157 = 1.0/cr_Ct156;
  const double cr_Ct158 = cr_Ct157*cr_Ct30;
  const double cr_Ct159 = cr_Ct158/sqrt(-cr_Ct154);
  const double cr_Ct160 = 0.33333333333333331*acos(5.196152422706632*cr_Ct152*cr_Ct159);
  const double cr_Ct161 = cos(cr_Ct160);
  const double cr_Ct162 = 1.1547005383792515*cr_Ct161;
  const double cr_Ct163 = cr_Ct141*cr_Ct162;
  const double cr_Ct164 = pow(cr_Ct155, -3.0/2.0);
  const double cr_Ct165 = 0.66666666666666663*pow(Young, 3)/pow(cr_Ct25, 3);
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
  const double cr_Ct212 = 0.074074074074074056*cr_Ct140*cr_Ct159*sin(cr_Ct160)/sqrt(0.037037037037037035 + cr_Ct64*pow(cr_Ct144*cr_Ct169 + cr_Ct151 - cr_Ct170*cr_Ct57 - cr_Ct173*cr_Ct59 + cr_Ct211*cr_Ct83, 2)/cr_Ct154);
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
  const double cr_Ct227 = 1.0/(Gf*Young/(characteristic_length*pow(threshold, 2)) - 0.5);
  const double cr_Ct228 = 1.0*cr_Ct227;
  const double cr_Ct229 = cr_Ct104*cr_Ct228;
  const double cr_Ct230 = cr_Ct17*cr_Ct229;
  const double cr_Ct231 = cr_Ct31*exp(cr_Ct228*(-0.5*cr_Ct103*cr_Ct28/threshold + 1.0));
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
  const double cr_Ct279 = pow(cr_Ct27, -3.0/2.0);
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
  const double cr_Ct297 = pow(Young, 4)/pow(cr_Ct25, 4);
  const double cr_Ct298 = 15.588457268119896*cr_Ct297/pow(cr_Ct27, 5.0/2.0);
  const double cr_Ct299 = cr_Ct295*r_strain[3];
  const double cr_Ct300 = pow(r_strain[3], 3);
  const double cr_Ct301 = cr_Ct298*r_strain[3];
  const double cr_Ct302 = cr_Ct19*cr_Ct22;
  const double cr_Ct303 = cr_Ct2*cr_Ct302;
  const double cr_Ct304 = cr_Ct29*r_strain[3];
  const double cr_Ct305 = pow(cr_Ct24, -2);
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
  const double cr_Ct324 = 0.074074074074074056*cr_Ct30*cr_Ct50*cr_Ct82*sin(cr_Ct100)/sqrt(-cr_Ct35*pow(-cr_Ct302*cr_Ct321 - cr_Ct321*cr_Ct322 - cr_Ct323 + cr_Ct35*cr_Ct84 + cr_Ct98, 2)/cr_Ct81 + 0.037037037037037035);
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
  const double cr_Ct341 = -cr_Ct113*r_strain[3] + cr_Ct219*cr_Ct325 - cr_Ct340*(1.0*cr_Ct326 + 1.0*cr_Ct339);
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
  const double cr_Ct363 = -cr_Ct113*r_strain[4] + cr_Ct219*cr_Ct357 - cr_Ct340*(1.0*cr_Ct358 + 1.0*cr_Ct362);
  const double cr_Ct364 = cr_Ct165*r_strain[5];
  const double cr_Ct365 = cr_Ct279*cr_Ct364;
  const double cr_Ct366 = 0.57735026918962573*r_strain[5];
  const double cr_Ct367 = 10.392304845413264*r_strain[5];
  const double cr_Ct368 = pow(r_strain[5], 3);
  const double cr_Ct369 = cr_Ct10*cr_Ct365 + cr_Ct293*cr_Ct366 - cr_Ct324*(-cr_Ct22*cr_Ct301*r_strain[4] + cr_Ct29*cr_Ct306*cr_Ct346 - cr_Ct295*cr_Ct33*cr_Ct367 + cr_Ct298*cr_Ct33*cr_Ct368 + cr_Ct298*cr_Ct349*r_strain[5] + cr_Ct299*cr_Ct350 + cr_Ct346*cr_Ct351 - cr_Ct352*r_strain[5] - cr_Ct354*r_strain[5] + cr_Ct355*r_strain[5] + cr_Ct356*r_strain[5]) + cr_Ct33*cr_Ct365 + cr_Ct34*cr_Ct365;
  const double cr_Ct370 = cr_Ct223*r_strain[5];
  const double cr_Ct371 = cr_Ct164*cr_Ct364;
  const double cr_Ct372 = cr_Ct156*(-cr_Ct146*cr_Ct371 - cr_Ct147*cr_Ct371 - cr_Ct168*cr_Ct371 + cr_Ct212*(-cr_Ct169*cr_Ct367*cr_Ct59 + cr_Ct180*r_strain[5] + cr_Ct181*cr_Ct368 + cr_Ct183*cr_Ct334*r_strain[5] - cr_Ct184*r_strain[5] + cr_Ct307*r_strain[4] - cr_Ct330*r_strain[4] - cr_Ct336*r_strain[5] - cr_Ct338*r_strain[5] - cr_Ct360*r_strain[5] + cr_Ct361*r_strain[5]) + cr_Ct328*cr_Ct366);
  const double cr_Ct373 = -cr_Ct113*r_strain[5] + cr_Ct219*cr_Ct369 - cr_Ct340*(1.0*cr_Ct370 + 1.0*cr_Ct372);
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
  const double cr_Ct388 = 1.0*cr_Ct105;
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
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

template<>
void AutomaticDifferentiationTangentUtilities<VonMisesYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{

  const auto &r_props = rValues.GetMaterialProperties();
  const double Young = r_props[YOUNG_MODULUS];
  const double nu = r_props[POISSON_RATIO];
  const double Gf = r_props[FRACTURE_ENERGY];
  const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());
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
  const double cr_Ct42 = cr_Ct0*(cr_Ct1*(-1.0 + cr_Ct41/sqrt(cr_Ct40)) + 1.0);
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
  const double cr_Ct74 = cr_Ct1*(cr_Ct41/sqrt(cr_Ct39*(cr_Ct19*std::pow(cr_Ct70, 2) + cr_Ct19*std::pow(cr_Ct72, 2) + cr_Ct19*std::pow(cr_Ct73, 2) + cr_Ct37)) - 1.0);
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
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<TrescaYieldSurface<TrescaPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  const auto &r_props = rValues.GetMaterialProperties();
  const double Young = r_props[YOUNG_MODULUS];
  const double nu = r_props[POISSON_RATIO];
  const double Gf = r_props[FRACTURE_ENERGY];
  const double phi = r_props[FRICTION_ANGLE];
  const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());
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
  const double cr_Ct48 = sqrt(cr_Ct47*(std::pow(cr_Ct29, 2)*cr_Ct31 + std::pow(cr_Ct37, 2)*cr_Ct39 + cr_Ct39*std::pow(cr_Ct42, 2) + cr_Ct46));
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
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<SimoJuYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<VonMisesYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
    const auto &r_props = rValues.GetMaterialProperties();
    const double Young = r_props[YOUNG_MODULUS];
    const double nu = r_props[POISSON_RATIO];
    const double Gf = r_props[FRACTURE_ENERGY];
    const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());
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
    const double cr_Ct46 = 1.0*cr_Ct45;
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
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<TrescaYieldSurface<TrescaPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<RankineYieldSurface<RankinePlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<SimoJuYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

template<>
void AutomaticDifferentiationTangentUtilities<VonMisesYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  const auto &r_props = rValues.GetMaterialProperties();
  const double Young = r_props[YOUNG_MODULUS];
  const double nu = r_props[POISSON_RATIO];
  const double Gf = r_props[FRACTURE_ENERGY];
  const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLength(rValues.GetElementGeometry());
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
  const double cr_Ct30 = cr_Ct0*(cr_Ct1*(-1.0 + cr_Ct29/sqrt(cr_Ct28)) + 1.0);
  const double cr_Ct31 = cr_Ct10 - cr_Ct11;
  const double cr_Ct32 = 0.055555555555555552*std::pow(cr_Ct16, 1.0)*cr_Ct4;
  const double cr_Ct33 = 3.0*nu;
  const double cr_Ct34 = cr_Ct33 -1.0;
  const double cr_Ct35 = 0.11111111111111112*std::pow(cr_Ct18, 1.0);
  const double cr_Ct36 = 0.11111111111111112*std::pow(cr_Ct26, 1.0);
  const double cr_Ct37 = cr_Ct24*(cr_Ct33 - 2.0);
  const double cr_Ct38 = cr_Ct32 + cr_Ct34*cr_Ct35*cr_Ct4 + cr_Ct36*cr_Ct37;
  const double cr_Ct39 = cr_Ct1*cr_Ct29/std::pow(cr_Ct28, 3.0/2.0);
  const double cr_Ct40 = cr_Ct39*cr_Ct6;
  const double cr_Ct41 = cr_Ct13 + cr_Ct21;
  const double cr_Ct42 = cr_Ct25*(cr_Ct23 + cr_Ct41);
  const double cr_Ct43 = cr_Ct25*(-cr_Ct17 - 0.5*cr_Ct22 + cr_Ct41);
  const double cr_Ct44 = cr_Ct1*(cr_Ct29/sqrt(cr_Ct27 + 0.055555555555555552*std::pow(cr_Ct42, 2.0) + 0.22222222222222224*std::pow(cr_Ct43, 2.0) + std::pow(cr_Ct25*r_strain[2]*(cr_Ct19 + 0.5), 2.0)) - 1.0) + 1.0;
  const double cr_Ct45 = cr_Ct44*nu;
  const double cr_Ct46 = 0.055555555555555552*std::pow(cr_Ct42, 1.0);
  const double cr_Ct47 = cr_Ct34*cr_Ct36;
  const double cr_Ct48 = 2.0 - cr_Ct33;
  const double cr_Ct49 = 0.11111111111111112*std::pow(cr_Ct43, 1.0);
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
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<TrescaYieldSurface<TrescaPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<RankineYieldSurface<RankinePlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<SimoJuYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>::
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

// 3D exponential
template class AutomaticDifferentiationTangentUtilities<VonMisesYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<TrescaYieldSurface<TrescaPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<RankineYieldSurface<RankinePlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<SimoJuYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>;
// 3D Linear
template class AutomaticDifferentiationTangentUtilities<VonMisesYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<TrescaYieldSurface<TrescaPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<RankineYieldSurface<RankinePlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<SimoJuYieldSurface<VonMisesPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>, static_cast<SizeType>(SofteningType::Linear)>;

// 2D exponential
template class AutomaticDifferentiationTangentUtilities<VonMisesYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<TrescaYieldSurface<TrescaPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<RankineYieldSurface<RankinePlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<SimoJuYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>;
template class AutomaticDifferentiationTangentUtilities<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Exponential)>;
// 2D Linear
template class AutomaticDifferentiationTangentUtilities<VonMisesYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<TrescaYieldSurface<TrescaPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<RankineYieldSurface<RankinePlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<SimoJuYieldSurface<VonMisesPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>;
template class AutomaticDifferentiationTangentUtilities<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<3>>, static_cast<SizeType>(SofteningType::Linear)>;
} // namespace Kratos