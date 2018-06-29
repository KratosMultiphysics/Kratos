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

  // template< unsigned int TDim >
  // void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo)
  // {
  //   KRATOS_TRY; 
  //   KRATOS_CATCH( "" );
  // }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::ComputeMaterialParameters(double& Density,
										       double& DeviatoricCoeff,
										       double& VolumetricCoeff,
										       ProcessInfo &currentProcessInfo,
										       ElementalVariables &rElementalVariables)
  {

  }

  

  template< unsigned int TDim>
  double TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::ComputeNonLinearViscosity(double & equivalentStrainRate)
  {
    double FluidViscosity=0;
    return FluidViscosity;
  }
  

  template< unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::ComputeMaterialParametersGranularGas(double& Density,
												  double& DeviatoricCoeff,
												  double& VolumetricCoeff,
												  ProcessInfo &currentProcessInfo,
												  ElementalVariables &rElementalVariables)
  {

  }

  
  template< unsigned int TDim>
  double TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::ComputeJopMuIrheologyViscosity(ElementalVariables & rElementalVariables)
  {
    double FluidViscosity=0;

    return FluidViscosity;
  }


  template< unsigned int TDim>
  double TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::ComputeBercovierMuIrheologyViscosity(ElementalVariables & rElementalVariables)
  {
        double FluidViscosity=0;
  
    return FluidViscosity;
  }

  template< unsigned int TDim>
  double TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::ComputeBarkerMuIrheologyViscosity(ElementalVariables & rElementalVariables)
  {
     double FluidViscosity=0;
  
    return FluidViscosity;
  }

  
  template< unsigned int TDim>
  double TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::ComputeBarkerBercovierMuIrheologyViscosity(ElementalVariables & rElementalVariables)
  {

    double FluidViscosity=0;
   
 
    return FluidViscosity;
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
    if(VISCOSITY.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"VISCOSITY Key is 0. Check that the application was correctly registered.","");
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
        if(this->GetGeometry()[i].SolutionStepsDataHas(VISCOSITY) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing VISCOSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
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


  template<>
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<2>::ComputeMeanValueMaterialTangentMatrix(ElementalVariables & rElementalVariables,double& MeanValue,const ShapeFunctionDerivativesType& rDN_DX,const double secondLame,double & bulkModulus,const double Weight,double& MeanValueMass,const double TimeStep)
  {
   
  }

  template<>
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<3>::ComputeMeanValueMaterialTangentMatrix(ElementalVariables & rElementalVariables,double& MeanValue,const ShapeFunctionDerivativesType& rDN_DX,const double secondLame,double & bulkModulus,const double Weight,double& MeanValueMass,const double TimeStep)
  {

  
  }


  template<>
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<2>::ComputeBulkReductionCoefficient(MatrixType MassMatrix,
											  MatrixType StiffnessMatrix,
											  double& meanValueStiff,
											  double& bulkCoefficient,
											  double timeStep)
  {
  
  }


  template<>
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<3>::ComputeBulkReductionCoefficient(MatrixType MassMatrix,
											  MatrixType StiffnessMatrix,
											  double& meanValueStiff,
											  double& bulkCoefficient,
											  double timeStep)
  {
  
  }



  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::ComputeBulkMatrix(Matrix& BulkMatrix,
									       const ShapeFunctionsType& rN,
									       const double Weight)
  {
  
  }

  template< >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<2>::ComputeBulkMatrixConsistent(Matrix& BulkMatrix,
										      const double Weight)
  {
  
  }

  template< >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<3>::ComputeBulkMatrixConsistent(Matrix& BulkMatrix,
										      const double Weight)
  {
 
  }
  

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::ComputeBulkMatrixLump(Matrix& BulkMatrix,
										   const double Weight)
  {
    
  }



  template< >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<2>::ComputeBoundLHSMatrix(Matrix& BoundLHSMatrix,
										const ShapeFunctionsType& rN,
										const double Weight)
  {
  
  }

  template<  >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<3>::ComputeBoundLHSMatrix(Matrix& BoundLHSMatrix,
										const ShapeFunctionsType& rN,
										const double Weight)
  {
  

  }
 


  template< >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<2>::ComputeBoundRHSVector(VectorType& BoundRHSVector,
										const ShapeFunctionsType& rN,
										const double TimeStep,
										const double BoundRHSCoeffAcc,
										const double BoundRHSCoeffDev)
  {
  
  
  }



  template< >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<3>::ComputeBoundRHSVector(VectorType& BoundRHSVector,
										const ShapeFunctionsType& rN,
										const double TimeStep,
										const double BoundRHSCoeffAcc,
										const double BoundRHSCoeffDev)
  {
   
  }



  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::CalculateTauFIC(double& Tau,
									     double ElemSize,
									     const double Density,
									     const double Viscosity,
									     const ProcessInfo& rCurrentProcessInfo)
  {
    
  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::AddStabilizationMatrixLHS(MatrixType& rLeftHandSideMatrix,
										       Matrix& BulkAccMatrix,
										       const ShapeFunctionsType& rN,
										       const double Weight)
  {
   
  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::ComputeStabLaplacianMatrix(MatrixType& StabLaplacianMatrix,
											const ShapeFunctionDerivativesType& rDN_DX,
											const double Weight)
										
  {
 
  }



  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::AddStabilizationNodalTermsLHS(MatrixType& rLeftHandSideMatrix,
											   const double Tau,
											   const double Weight,
											   const ShapeFunctionDerivativesType& rDN_DX,
											   const SizeType i)
  {
 
  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::AddStabilizationNodalTermsRHS(VectorType& rRightHandSideVector,
											   const double Tau,
											   const double Density,
											   const double Weight,
											   const ShapeFunctionDerivativesType& rDN_DX,
											   const SizeType i)
  {


  }


  template<>
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<2>::GetPositions(Vector& rValues,const ProcessInfo& rCurrentProcessInfo,const double theta)
  {
    // GeometryType& rGeom = this->GetGeometry();
    // const SizeType NumNodes = rGeom.PointsNumber();
    // const SizeType LocalSize = 2*NumNodes;

    // if (rValues.size() != LocalSize) rValues.resize(LocalSize);

    // SizeType Index = 0;

    // for (SizeType i = 0; i < NumNodes; ++i)
    //   {
    // 	rValues[Index++] = rGeom[i].X();
    // 	rValues[Index++] = rGeom[i].Y();
    //   }
  }



  template<>
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<3>::GetPositions(Vector& rValues,const ProcessInfo& rCurrentProcessInfo,const double theta)
  {
    // GeometryType& rGeom = this->GetGeometry();
    // const SizeType NumNodes = rGeom.PointsNumber();
    // const SizeType LocalSize = 3*NumNodes;

    // if (rValues.size() != LocalSize) rValues.resize(LocalSize);

    // SizeType Index = 0;

    // for (SizeType i = 0; i < NumNodes; ++i)
    //   {
    // 	rValues[Index++] = rGeom[i].X();
    //     rValues[Index++] = rGeom[i].Y();
    //     rValues[Index++] = rGeom[i].Z();
    //   }
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


  template < > 
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<2>:: CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables,double TimeStep, unsigned int g)
  {

  

  }
 

  template < > 
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<3>:: CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables, double TimeStep, unsigned int g)
  {


  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::CalculateLocalContinuityEqForPressure(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {



  }
  

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::GetPressureVelocityValues(Vector& rValues,
										       const int Step)
  {
    // GeometryType& rGeom = this->GetGeometry();
    // const SizeType NumNodes = rGeom.PointsNumber();

    // if (rValues.size() != NumNodes) rValues.resize(NumNodes);

    // for (SizeType i = 0; i < NumNodes; ++i){
    //   rValues[i] = rGeom[i].FastGetSolutionStepValue(PRESSURE_VELOCITY,Step);

    // }
  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<TDim>::GetPressureAccelerationValues(Vector& rValues,
											   const int Step)
  {
    // GeometryType& rGeom = this->GetGeometry();
    // const SizeType NumNodes = rGeom.PointsNumber();

    // if (rValues.size() != NumNodes) rValues.resize(NumNodes);

    // for (SizeType i = 0; i < NumNodes; ++i){
    //   rValues[i] = rGeom[i].FastGetSolutionStepValue(PRESSURE_ACCELERATION,Step);

    // }
  }




  template class TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<2>;
  template class TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<3>;

}
