//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 0.0 $
//
//   Implementation of the Gauss-Seidel two step Updated Lagrangian Velocity-Pressure element
//     ( There is a ScalingConstant to multiply the mass balance equation for a number because i read it somewhere)
//

// System includes

// External includes

// Project includes
#include "custom_elements/updated_lagrangian_V_implicit_solid_element.h"
#include "includes/cfd_variables.h"

namespace Kratos
{

/*
   * public UpdatedLagrangianVImplicitSolidElement<TDim> functions
   */

template <unsigned int TDim>
Element::Pointer UpdatedLagrangianVImplicitSolidElement<TDim>::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
{
  // return Element::Pointer( BaseType::Clone(NewId,rThisNodes) );
  UpdatedLagrangianVImplicitSolidElement NewElement(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());

  if (NewElement.mCurrentTotalCauchyStress.size() != this->mCurrentTotalCauchyStress.size())
    NewElement.mCurrentTotalCauchyStress.resize(this->mCurrentTotalCauchyStress.size());

  for (unsigned int i = 0; i < this->mCurrentTotalCauchyStress.size(); i++)
  {
    NewElement.mCurrentTotalCauchyStress[i] = this->mCurrentTotalCauchyStress[i];
  }

  if (NewElement.mCurrentDeviatoricCauchyStress.size() != this->mCurrentDeviatoricCauchyStress.size())
    NewElement.mCurrentDeviatoricCauchyStress.resize(this->mCurrentDeviatoricCauchyStress.size());

  for (unsigned int i = 0; i < this->mCurrentDeviatoricCauchyStress.size(); i++)
  {
    NewElement.mCurrentDeviatoricCauchyStress[i] = this->mCurrentDeviatoricCauchyStress[i];
  }

  if (NewElement.mUpdatedTotalCauchyStress.size() != this->mUpdatedTotalCauchyStress.size())
    NewElement.mUpdatedTotalCauchyStress.resize(this->mUpdatedTotalCauchyStress.size());

  for (unsigned int i = 0; i < this->mUpdatedTotalCauchyStress.size(); i++)
  {
    NewElement.mUpdatedTotalCauchyStress[i] = this->mUpdatedTotalCauchyStress[i];
  }

  if (NewElement.mUpdatedDeviatoricCauchyStress.size() != this->mUpdatedDeviatoricCauchyStress.size())
    NewElement.mUpdatedDeviatoricCauchyStress.resize(this->mUpdatedDeviatoricCauchyStress.size());

  for (unsigned int i = 0; i < this->mUpdatedDeviatoricCauchyStress.size(); i++)
  {
    NewElement.mUpdatedDeviatoricCauchyStress[i] = this->mUpdatedDeviatoricCauchyStress[i];
  }

  NewElement.SetData(this->GetData());
  NewElement.SetFlags(this->GetFlags());

  return Element::Pointer(new UpdatedLagrangianVImplicitSolidElement(NewElement));
}

template <>
void UpdatedLagrangianVImplicitSolidElement<2>::CalcElasticPlasticCauchySplitted(
    ElementalVariables &rElementalVariables, double TimeStep, unsigned int g, const ProcessInfo &rCurrentProcessInfo,
    double &Density, double &DeviatoricCoeff, double &VolumetricCoeff) {

    mpConstitutiveLaw = this->GetProperties().GetValue(CONSTITUTIVE_LAW);
    auto constitutive_law_values =
        ConstitutiveLaw::Parameters(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);

    Flags &constitutive_law_options = constitutive_law_values.GetOptions();
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    rElementalVariables.CurrentTotalCauchyStress = this->mCurrentTotalCauchyStress[g];
    rElementalVariables.CurrentDeviatoricCauchyStress = this->mCurrentDeviatoricCauchyStress[g];

    const Vector &r_shape_functions = row((this->GetGeometry()).ShapeFunctionsValues(), g);
    constitutive_law_values.SetShapeFunctionsValues(r_shape_functions);
    constitutive_law_values.SetStrainVector(rElementalVariables.SpatialDefRate);
    constitutive_law_values.SetStressVector(rElementalVariables.CurrentDeviatoricCauchyStress);

    mpConstitutiveLaw->CalculateMaterialResponseCauchy(constitutive_law_values);

    Density = mpConstitutiveLaw->CalculateValue(constitutive_law_values, DENSITY, Density);

    double poisson_ratio = mpConstitutiveLaw->CalculateValue(constitutive_law_values, POISSON_RATIO, poisson_ratio);
    double young_modulus = mpConstitutiveLaw->CalculateValue(constitutive_law_values, YOUNG_MODULUS, young_modulus);
    const double time_step = rCurrentProcessInfo[DELTA_TIME];
    DeviatoricCoeff = time_step * young_modulus / (2.0 * (1 + poisson_ratio));
    VolumetricCoeff =
        time_step * poisson_ratio * young_modulus / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio)) +
        2.0 / 3.0 * DeviatoricCoeff;

    const double current_first_lame = VolumetricCoeff - 2.0 / 3.0 * DeviatoricCoeff;

    this->mMaterialDeviatoricCoefficient = DeviatoricCoeff;
    this->mMaterialVolumetricCoefficient = VolumetricCoeff;
    this->mMaterialDensity = Density;

    rElementalVariables.UpdatedDeviatoricCauchyStress[0] = rElementalVariables.CurrentDeviatoricCauchyStress[0];
    rElementalVariables.UpdatedDeviatoricCauchyStress[1] = rElementalVariables.CurrentDeviatoricCauchyStress[1];
    rElementalVariables.UpdatedDeviatoricCauchyStress[2] = rElementalVariables.CurrentDeviatoricCauchyStress[2];

    rElementalVariables.UpdatedTotalCauchyStress[0] = +current_first_lame * rElementalVariables.VolumetricDefRate +
                                                      2.0 * DeviatoricCoeff * rElementalVariables.SpatialDefRate[0] +
                                                      rElementalVariables.CurrentTotalCauchyStress[0];
    rElementalVariables.UpdatedTotalCauchyStress[1] = current_first_lame * rElementalVariables.VolumetricDefRate +
                                                      2.0 * DeviatoricCoeff * rElementalVariables.SpatialDefRate[1] +
                                                      rElementalVariables.CurrentTotalCauchyStress[1];
    rElementalVariables.UpdatedTotalCauchyStress[2] =
        2.0 * DeviatoricCoeff * rElementalVariables.SpatialDefRate[2] + rElementalVariables.CurrentTotalCauchyStress[2];

    this->SetValue(CAUCHY_STRESS_VECTOR, rElementalVariables.UpdatedTotalCauchyStress);

    this->mUpdatedTotalCauchyStress[g] = rElementalVariables.UpdatedTotalCauchyStress;
    this->mUpdatedDeviatoricCauchyStress[g] = rElementalVariables.UpdatedDeviatoricCauchyStress;
}

template <>
void UpdatedLagrangianVImplicitSolidElement<3>::CalcElasticPlasticCauchySplitted(
    ElementalVariables &rElementalVariables, double TimeStep, unsigned int g, const ProcessInfo &rCurrentProcessInfo,
    double &Density, double &DeviatoricCoeff, double &VolumetricCoeff) {

    mpConstitutiveLaw = this->GetProperties().GetValue(CONSTITUTIVE_LAW);
    auto constitutive_law_values =
        ConstitutiveLaw::Parameters(this->GetGeometry(), this->GetProperties(), rCurrentProcessInfo);

    Flags &constitutive_law_options = constitutive_law_values.GetOptions();
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_STRESS, true);
    constitutive_law_options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

    rElementalVariables.CurrentTotalCauchyStress = this->mCurrentTotalCauchyStress[g];
    rElementalVariables.CurrentDeviatoricCauchyStress = this->mCurrentDeviatoricCauchyStress[g];

    const Vector &r_shape_functions = row((this->GetGeometry()).ShapeFunctionsValues(), g);
    constitutive_law_values.SetShapeFunctionsValues(r_shape_functions);
    constitutive_law_values.SetStrainVector(rElementalVariables.SpatialDefRate);
    constitutive_law_values.SetStressVector(rElementalVariables.CurrentDeviatoricCauchyStress);

    mpConstitutiveLaw->CalculateMaterialResponseCauchy(constitutive_law_values);

    Density = mpConstitutiveLaw->CalculateValue(constitutive_law_values, DENSITY, Density);

    double poisson_ratio = mpConstitutiveLaw->CalculateValue(constitutive_law_values, POISSON_RATIO, poisson_ratio);
    double young_modulus = mpConstitutiveLaw->CalculateValue(constitutive_law_values, YOUNG_MODULUS, young_modulus);
    const double time_step = rCurrentProcessInfo[DELTA_TIME];
    DeviatoricCoeff = time_step * young_modulus / (2.0 * (1 + poisson_ratio));
    VolumetricCoeff =
        time_step * poisson_ratio * young_modulus / ((1.0 + poisson_ratio) * (1.0 - 2.0 * poisson_ratio)) +
        2.0 / 3.0 * DeviatoricCoeff;

    const double current_first_lame = VolumetricCoeff - 2.0 / 3.0 * DeviatoricCoeff;

    this->mMaterialDeviatoricCoefficient = DeviatoricCoeff;
    this->mMaterialVolumetricCoefficient = VolumetricCoeff;
    this->mMaterialDensity = Density;

    rElementalVariables.UpdatedDeviatoricCauchyStress[0] = rElementalVariables.CurrentDeviatoricCauchyStress[0];
    rElementalVariables.UpdatedDeviatoricCauchyStress[1] = rElementalVariables.CurrentDeviatoricCauchyStress[1];
    rElementalVariables.UpdatedDeviatoricCauchyStress[2] = rElementalVariables.CurrentDeviatoricCauchyStress[2];
    rElementalVariables.UpdatedDeviatoricCauchyStress[3] = rElementalVariables.CurrentDeviatoricCauchyStress[3];
    rElementalVariables.UpdatedDeviatoricCauchyStress[4] = rElementalVariables.CurrentDeviatoricCauchyStress[4];
    rElementalVariables.UpdatedDeviatoricCauchyStress[5] = rElementalVariables.CurrentDeviatoricCauchyStress[5];

    rElementalVariables.UpdatedTotalCauchyStress[0] = +current_first_lame * rElementalVariables.VolumetricDefRate +
                                                      2.0 * DeviatoricCoeff * rElementalVariables.SpatialDefRate[0] +
                                                      rElementalVariables.CurrentTotalCauchyStress[0];
    rElementalVariables.UpdatedTotalCauchyStress[1] = current_first_lame * rElementalVariables.VolumetricDefRate +
                                                      2.0 * DeviatoricCoeff * rElementalVariables.SpatialDefRate[1] +
                                                      rElementalVariables.CurrentTotalCauchyStress[1];
    rElementalVariables.UpdatedTotalCauchyStress[2] = current_first_lame * rElementalVariables.VolumetricDefRate +
                                                      2.0 * DeviatoricCoeff * rElementalVariables.SpatialDefRate[2] +
                                                      rElementalVariables.CurrentTotalCauchyStress[2];
    rElementalVariables.UpdatedTotalCauchyStress[3] =
        2.0 * DeviatoricCoeff * rElementalVariables.SpatialDefRate[3] + rElementalVariables.CurrentTotalCauchyStress[3];
    rElementalVariables.UpdatedTotalCauchyStress[4] =
        2.0 * DeviatoricCoeff * rElementalVariables.SpatialDefRate[4] + rElementalVariables.CurrentTotalCauchyStress[4];
    rElementalVariables.UpdatedTotalCauchyStress[5] =
        2.0 * DeviatoricCoeff * rElementalVariables.SpatialDefRate[5] + rElementalVariables.CurrentTotalCauchyStress[5];

    this->SetValue(CAUCHY_STRESS_VECTOR, rElementalVariables.UpdatedTotalCauchyStress);

    this->mUpdatedTotalCauchyStress[g] = rElementalVariables.UpdatedTotalCauchyStress;
    this->mUpdatedDeviatoricCauchyStress[g] = rElementalVariables.UpdatedDeviatoricCauchyStress;
}

template class UpdatedLagrangianVImplicitSolidElement<2>;
template class UpdatedLagrangianVImplicitSolidElement<3>;

} // namespace Kratos
