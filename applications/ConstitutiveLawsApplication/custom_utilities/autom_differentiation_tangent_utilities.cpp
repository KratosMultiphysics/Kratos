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
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void AutomaticDifferentiationTangentUtilities<RankineYieldSurface<RankinePlasticPotential<6>>, static_cast<SizeType>(SofteningType::Exponential)>::
CalculateTangentTensorAutomDiffIsotropicDamage(
  ConstitutiveLaw::Parameters rValues
  )
{
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
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
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
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
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
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
    const double cr_Ct2 = cr_Ct1*r_strain[2];
    const double cr_Ct3 = 2*nu;
    const double cr_Ct4 = 1.0/(cr_Ct3 - 1);
    const double cr_Ct5 = nu + 1.0;
    const double cr_Ct6 = Young/cr_Ct5;
    const double cr_Ct7 = cr_Ct4*cr_Ct6;
    const double cr_Ct8 = std::pow(cr_Ct2*cr_Ct7, 2.0);
    const double cr_Ct9 = nu*r_strain[1];
    const double cr_Ct10 = cr_Ct0*r_strain[0];
    const double cr_Ct11 = cr_Ct10 - cr_Ct9;
    const double cr_Ct12 = nu*r_strain[0];
    const double cr_Ct13 = cr_Ct0*r_strain[1];
    const double cr_Ct14 = -cr_Ct12 + cr_Ct13;
    const double cr_Ct15 = cr_Ct7*(cr_Ct11 + cr_Ct14 + cr_Ct2);
    const double cr_Ct16 = -nu;
    const double cr_Ct17 = r_strain[2]*(cr_Ct16 + 0.5);
    const double cr_Ct18 = -0.5*cr_Ct17;
    const double cr_Ct19 = cr_Ct16 + 1.0;
    const double cr_Ct20 = cr_Ct19*r_strain[1];
    const double cr_Ct21 = cr_Ct19*r_strain[0];
    const double cr_Ct22 = cr_Ct21 + cr_Ct9;
    const double cr_Ct23 = 1.0/(1 - cr_Ct3);
    const double cr_Ct24 = cr_Ct23*cr_Ct6;
    const double cr_Ct25 = cr_Ct24*(-0.5*cr_Ct12 + cr_Ct18 - 0.5*cr_Ct20 + cr_Ct22);
    const double cr_Ct26 = cr_Ct12 + cr_Ct20;
    const double cr_Ct27 = cr_Ct24*(cr_Ct18 - 0.5*cr_Ct21 + cr_Ct26 - 0.5*cr_Ct9);
    const double cr_Ct28 = 0.055555555555555552*std::pow(cr_Ct15, 2.0) + 0.22222222222222224*std::pow(cr_Ct25, 2.0) + 0.22222222222222224*std::pow(cr_Ct27, 2.0) + cr_Ct8;
    const double cr_Ct29 = std::sqrt(cr_Ct28);
    const double cr_Ct30 = 0.57735026918962584*threshold;
    const double cr_Ct31 = cr_Ct30/cr_Ct29;
    const double cr_Ct32 = cr_Ct0*cr_Ct31;
    const double cr_Ct33 = 0.055555555555555552*std::pow(cr_Ct24*(cr_Ct17 + cr_Ct22 + cr_Ct26), 1.0);
    const double cr_Ct34 = 3.0*nu;
    const double cr_Ct35 = 2.0 - cr_Ct34;
    const double cr_Ct36 = std::pow(cr_Ct25, 1.0);
    const double cr_Ct37 = 0.11111111111111112*cr_Ct36;
    const double cr_Ct38 = cr_Ct34 - 0.99999999999999989;
    const double cr_Ct39 = std::pow(cr_Ct27, 1.0);
    const double cr_Ct40 = 0.11111111111111112*cr_Ct39;
    const double cr_Ct41 = cr_Ct38*cr_Ct40;
    const double cr_Ct42 = cr_Ct33 + cr_Ct35*cr_Ct37 + cr_Ct41;
    const double cr_Ct43 = cr_Ct30/std::pow(cr_Ct28, 3.0/2.0);
    const double cr_Ct44 = cr_Ct24*cr_Ct43;
    const double cr_Ct45 = cr_Ct11*cr_Ct44;
    const double cr_Ct46 = 0.055555555555555552*std::pow(cr_Ct15, 1.0)*cr_Ct4;
    const double cr_Ct47 = cr_Ct23*(cr_Ct34 - 2.0);
    const double cr_Ct48 = cr_Ct37*cr_Ct47 + cr_Ct4*cr_Ct41 + cr_Ct46;
    const double cr_Ct49 = 1.0/cr_Ct28;
    const double cr_Ct50 = 1.0/(Gf*Young/(characteristic_length*std::pow(threshold, 2)) - 0.5);
    const double cr_Ct51 = 1.0*cr_Ct50;
    const double cr_Ct52 = cr_Ct49*cr_Ct51;
    const double cr_Ct53 = cr_Ct52*cr_Ct6;
    const double cr_Ct54 = cr_Ct11*cr_Ct53;
    const double cr_Ct55 = std::exp(cr_Ct51*(-1.7320508075688772*cr_Ct29/threshold + 1.0));
    const double cr_Ct56 = cr_Ct55*cr_Ct7;
    const double cr_Ct57 = -cr_Ct31*nu;
    const double cr_Ct58 = cr_Ct37*cr_Ct38;
    const double cr_Ct59 = cr_Ct33 + cr_Ct35*cr_Ct40 + cr_Ct58;
    const double cr_Ct60 = cr_Ct4*cr_Ct58 + cr_Ct40*cr_Ct47 + cr_Ct46;
    const double cr_Ct61 = cr_Ct8/r_strain[2];
    const double cr_Ct62 = 0.1111111111111111*cr_Ct1*cr_Ct7;
    const double cr_Ct63 = -cr_Ct1*cr_Ct46*cr_Ct6 + cr_Ct36*cr_Ct62 + cr_Ct39*cr_Ct62;
    const double cr_Ct64 = -cr_Ct61 + cr_Ct63;
    const double cr_Ct65 = cr_Ct43*cr_Ct64;
    const double cr_Ct66 = cr_Ct49*cr_Ct50;
    const double cr_Ct67 = cr_Ct56*(cr_Ct64*cr_Ct66 + cr_Ct65);
    const double cr_Ct68 = cr_Ct14*cr_Ct44;
    const double cr_Ct69 = cr_Ct14*cr_Ct53;
    const double cr_Ct70 = cr_Ct23*cr_Ct43;
    const double cr_Ct71 = std::pow(Young, 2)*cr_Ct2*cr_Ct4*cr_Ct55/std::pow(cr_Ct5, 2);
    r_Ct(0,0)=cr_Ct56*(cr_Ct32 - cr_Ct42*cr_Ct45 + cr_Ct48*cr_Ct54);
    r_Ct(0,1)=cr_Ct56*(-cr_Ct45*cr_Ct59 + cr_Ct54*cr_Ct60 + cr_Ct57);
    r_Ct(0,2)=-cr_Ct67*(-cr_Ct10 + cr_Ct9);
    r_Ct(1,0)=cr_Ct56*(-cr_Ct42*cr_Ct68 + cr_Ct48*cr_Ct69 + cr_Ct57);
    r_Ct(1,1)=cr_Ct56*(cr_Ct32 - cr_Ct59*cr_Ct68 + cr_Ct60*cr_Ct69);
    r_Ct(1,2)=-cr_Ct67*(cr_Ct12 - cr_Ct13);
    r_Ct(2,0)=cr_Ct71*(-cr_Ct42*cr_Ct70 + cr_Ct48*cr_Ct66);
    r_Ct(2,1)=cr_Ct71*(-cr_Ct59*cr_Ct70 + cr_Ct60*cr_Ct66);
    r_Ct(2,2)=cr_Ct1*cr_Ct56*(cr_Ct31 + cr_Ct52*r_strain[2]*(-cr_Ct61 + cr_Ct63) + cr_Ct65*r_strain[2]);
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
  KRATOS_ERROR << "This autom differentiation has not been implemented yet..." << std::endl;
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


// template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<VonMisesYieldSurface<VonMisesPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<ModifiedMohrCoulombYieldSurface<ModifiedMohrCoulombPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<TrescaYieldSurface<TrescaPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<DruckerPragerYieldSurface<DruckerPragerPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<RankineYieldSurface<RankinePlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<SimoJuYieldSurface<VonMisesPlasticPotential<6>>>>;
// template class GenericSmallStrainIsotropicDamage <GenericConstitutiveLawIntegratorDamage<MohrCoulombYieldSurface<MohrCoulombPlasticPotential<6>>>>;