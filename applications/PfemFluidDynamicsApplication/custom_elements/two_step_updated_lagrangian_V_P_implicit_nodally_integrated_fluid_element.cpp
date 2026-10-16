//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_nodally_integrated_fluid_element.h"
#include "includes/cfd_variables.h"
#include <cmath>

namespace Kratos
{

  /*
   * public TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim> functions
   */

  template <unsigned int TDim>
  Element::Pointer TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::Clone(IndexType NewId, NodesArrayType const &rThisNodes) const
  {

    TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement NewElement(NewId, this->GetGeometry().Create(rThisNodes), this->pGetProperties());

    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());

    return Element::Pointer(new TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement(NewElement));
  }

  template <unsigned int TDim>
  int TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo) const
  {
    KRATOS_TRY;

    // Base class checks for positive Jacobian and Id > 0
    int ierr = Element::Check(rCurrentProcessInfo);
    if (ierr != 0)
      return ierr;

    // Check that all required variables have been registered
    KRATOS_ERROR_IF(VELOCITY.Key() == 0) << "VELOCITY Key is 0. Check that the application was correctly registered.";
    KRATOS_ERROR_IF(ACCELERATION.Key() == 0) << "ACCELERATION Key is 0. Check that the application was correctly registered.";
    KRATOS_ERROR_IF(PRESSURE.Key() == 0) << "PRESSURE Key is 0. Check that the application was correctly registered.";
    KRATOS_ERROR_IF(DENSITY.Key() == 0) << "DENSITY Key is 0. Check that the application was correctly registered.";
    KRATOS_ERROR_IF(BODY_FORCE.Key() == 0) << "BODY_FORCE Key is 0. Check that the application was correctly registered.";
    KRATOS_ERROR_IF(DYNAMIC_VISCOSITY.Key() == 0) << "DYNAMIC_VISCOSITY Key is 0. Check that the application was correctly registered.";
    KRATOS_ERROR_IF(DELTA_TIME.Key() == 0) << "DELTA_TIME Key is 0. Check that the application was correctly registered.";

    const GeometryType &r_geom = this->GetGeometry();
    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for (unsigned int i = 0; i < r_geom.size(); ++i)
    {
      KRATOS_ERROR_IF_NOT(r_geom[i].SolutionStepsDataHas(VELOCITY)) << "missing VELOCITY variable on solution step data for node ", r_geom[i].Id();
      KRATOS_ERROR_IF_NOT(r_geom[i].SolutionStepsDataHas(PRESSURE)) << "missing PRESSURE variable on solution step data for node ", r_geom[i].Id();
      KRATOS_ERROR_IF_NOT(r_geom[i].SolutionStepsDataHas(BODY_FORCE)) << "missing BODY_FORCE variable on solution step data for node ", r_geom[i].Id();
      KRATOS_ERROR_IF_NOT(r_geom[i].SolutionStepsDataHas(DENSITY)) << "missing DENSITY variable on solution step data for node ", r_geom[i].Id();
      KRATOS_ERROR_IF_NOT(r_geom[i].SolutionStepsDataHas(DYNAMIC_VISCOSITY)) << "missing DYNAMIC_VISCOSITY variable on solution step data for node ", r_geom[i].Id();
      KRATOS_ERROR_IF_NOT(r_geom[i].HasDofFor(VELOCITY_X) || r_geom[i].HasDofFor(VELOCITY_Y)) << "missing VELOCITY component DOF in node ", r_geom[i].Id();
      if constexpr (TDim == 3)
      {
        KRATOS_ERROR_IF_NOT(r_geom[i].HasDofFor(VELOCITY_Z)) << "Missing VELOCITY_Z component DOF in node ", this->GetGeometry()[i].Id();
      }
      KRATOS_ERROR_IF_NOT(r_geom[i].HasDofFor(PRESSURE)) << "missing PRESSURE DOF in node ", r_geom[i].Id();
    }

    // If this is a 2D problem, check that nodes are in XY plane
    if (r_geom.WorkingSpaceDimension() == 2)
    {
      for (unsigned int i = 0; i < this->GetGeometry().size(); ++i)
      {
        if (r_geom[i].Z() != 0.0)
          KRATOS_THROW_ERROR(std::invalid_argument, "Node with non-zero Z coordinate found. Id: ", r_geom[i].Id());
      }
    }

    return ierr;

    KRATOS_CATCH("");
  }

  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::InitializeElementalVariables(ElementalVariables &rElementalVariables)
  {
    KRATOS_TRY;

    unsigned int voigtsize = 3;
    if constexpr (TDim == 3)
    {
      voigtsize = 6;
    }
    rElementalVariables.voigtsize = voigtsize;
    rElementalVariables.ConstitutiveMatrix = ZeroMatrix(voigtsize, voigtsize);
    rElementalVariables.DetFgrad = 1.0;
    rElementalVariables.DetFgradVel = 1.0;
    rElementalVariables.DeviatoricInvariant = 1.0;
    rElementalVariables.EquivalentStrainRate = 1.0;
    rElementalVariables.VolumetricDefRate = 1.0;
    rElementalVariables.SpatialDefRate = ZeroVector(voigtsize);
    rElementalVariables.MDGreenLagrangeMaterial.resize(voigtsize, false);
    noalias(rElementalVariables.MDGreenLagrangeMaterial) = ZeroVector(voigtsize);
    rElementalVariables.Fgrad = ZeroMatrix(TDim, TDim);
    rElementalVariables.InvFgrad = ZeroMatrix(TDim, TDim);
    rElementalVariables.FgradVel = ZeroMatrix(TDim, TDim);
    rElementalVariables.InvFgradVel = ZeroMatrix(TDim, TDim);
    rElementalVariables.SpatialVelocityGrad = ZeroMatrix(TDim, TDim);

    rElementalVariables.MeanPressure = 0;
    rElementalVariables.CurrentTotalCauchyStress = ZeroVector(voigtsize);
    rElementalVariables.UpdatedTotalCauchyStress = ZeroVector(voigtsize);
    rElementalVariables.CurrentDeviatoricCauchyStress = ZeroVector(voigtsize);
    rElementalVariables.UpdatedDeviatoricCauchyStress = ZeroVector(voigtsize);

    KRATOS_CATCH("");
  }

  template class TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<2>;
  template class TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<3>;

} // namespace Kratos
