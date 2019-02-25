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
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_nodally_integrated_fluid_element.h"
#include "includes/cfd_variables.h"
#include <math.h>

namespace Kratos {

  /* 
   * public TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim> functions
   */

  template< unsigned int TDim >
  Element::Pointer TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {
 
    TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement NewElement(NewId, this->GetGeometry().Create( rThisNodes ), this->pGetProperties() );

    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());

    return Element::Pointer( new TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement(NewElement) );

  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::Initialize()
  {
    KRATOS_TRY; 
    KRATOS_CATCH( "" );
  }
  
  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo)
  {

  }


  
  template< unsigned int TDim >
  int TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo)
  {
    KRATOS_TRY;

    // Base class checks for positive Jacobian and Id > 0
    int ierr = Element::Check(rCurrentProcessInfo);
    if(ierr != 0) return ierr;

    // Check that all required variables have been registered
    if(VELOCITY.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"VELOCITY Key is 0. Check that the application was correctly registered.","");
    if(ACCELERATION.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"ACCELERATION Key is 0. Check that the application was correctly registered.","");
    if(PRESSURE.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"PRESSURE Key is 0. Check that the application was correctly registered.","");
    if(BODY_FORCE.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"BODY_FORCE Key is 0. Check that the application was correctly registered.","");
    if(DENSITY.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY Key is 0. Check that the application was correctly registered.","");
    if(DYNAMIC_VISCOSITY.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"DYNAMIC_VISCOSITY Key is 0. Check that the application was correctly registered.","");
    if(DELTA_TIME.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"DELTA_TIME Key is 0. Check that the application was correctly registered.","");

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
      {
        if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(BODY_FORCE) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing BODY_FORCE variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(DENSITY) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing DENSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(DYNAMIC_VISCOSITY) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing DYNAMIC_VISCOSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
           this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
           this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY component degree of freedom on node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE component degree of freedom on node ",this->GetGeometry()[i].Id());
      }
    
    // If this is a 2D problem, check that nodes are in XY plane
    if (this->GetGeometry().WorkingSpaceDimension() == 2)
      {
        for (unsigned int i=0; i<this->GetGeometry().size(); ++i)
	  {
            if (this->GetGeometry()[i].Z() != 0.0)
	      KRATOS_THROW_ERROR(std::invalid_argument,"Node with non-zero Z coordinate found. Id: ",this->GetGeometry()[i].Id());
	  }
      }

    return ierr;

    KRATOS_CATCH("");
  }


  template <  unsigned int TDim> 
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>:: InitializeElementalVariables(ElementalVariables & rElementalVariables)
  {
    KRATOS_TRY;

    unsigned int voigtsize  = 3;
    if( TDim == 3 )
      {
	voigtsize  = 6;
      }
    rElementalVariables.voigtsize=voigtsize;


    rElementalVariables.DetFgrad=1.0;

    rElementalVariables.DetFgradVel=1.0;

    rElementalVariables.DeviatoricInvariant=1.0;

    rElementalVariables.EquivalentStrainRate=1.0;

    rElementalVariables.VolumetricDefRate=1.0;

    rElementalVariables.SpatialDefRate= ZeroVector(voigtsize);

    rElementalVariables.MDGreenLagrangeMaterial.resize(voigtsize,false);

    noalias(rElementalVariables.MDGreenLagrangeMaterial) = ZeroVector(voigtsize);
  
    rElementalVariables.Fgrad = ZeroMatrix(TDim,TDim);

    rElementalVariables.InvFgrad= ZeroMatrix(TDim,TDim);

    rElementalVariables.FgradVel= ZeroMatrix(TDim,TDim);

    rElementalVariables.InvFgradVel= ZeroMatrix(TDim,TDim);

    rElementalVariables.SpatialVelocityGrad= ZeroMatrix(TDim,TDim);

    rElementalVariables.MeanPressure=0;

    rElementalVariables.CurrentTotalCauchyStress= ZeroVector(voigtsize);

    rElementalVariables.UpdatedTotalCauchyStress=  ZeroVector(voigtsize);

    rElementalVariables.CurrentDeviatoricCauchyStress=  ZeroVector(voigtsize);

    rElementalVariables.UpdatedDeviatoricCauchyStress= ZeroVector(voigtsize);

    KRATOS_CATCH("");

  }





  template class TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<2>;
  template class TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<3>;

}
