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
#include "custom_elements/two_step_updated_lagrangian_V_P_fluid_element.h"
#include "includes/cfd_variables.h"

namespace Kratos {

  /* 
   * public TwoStepUpdatedLagrangianVPFluidElement<TDim> functions
   */

template< unsigned int TDim >
Element::Pointer TwoStepUpdatedLagrangianVPFluidElement<TDim>::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
  // Element::Pointer p_clone= BaseType::Clone(NewId,rThisNodes);
  // return p_clone;

  TwoStepUpdatedLagrangianVPFluidElement NewElement(NewId, this->GetGeometry().Create( rThisNodes ), this->pGetProperties() );


  if ( NewElement.mCurrentFgrad.size() != this->mCurrentFgrad.size() )
    NewElement.mCurrentFgrad.resize(this->mCurrentFgrad.size());

  for(unsigned int i=0; i<this->mCurrentFgrad.size(); i++)
    {
      NewElement.mCurrentFgrad[i] = this->mCurrentFgrad[i];
    }

  if ( NewElement.mUpdatedFgrad.size() != this->mUpdatedFgrad.size() )
    NewElement.mUpdatedFgrad.resize(this->mUpdatedFgrad.size());

  for(unsigned int i=0; i<this->mUpdatedFgrad.size(); i++)
    {
      NewElement.mUpdatedFgrad[i] = this->mUpdatedFgrad[i];
    }

  if ( NewElement.mCurrentTotalCauchyStress.size() != this->mCurrentTotalCauchyStress.size() )
    NewElement.mCurrentTotalCauchyStress.resize(this->mCurrentTotalCauchyStress.size());

  for(unsigned int i=0; i<this->mCurrentTotalCauchyStress.size(); i++)
    {
      NewElement.mCurrentTotalCauchyStress[i] = this->mCurrentTotalCauchyStress[i];
    }


  if ( NewElement.mCurrentDeviatoricCauchyStress.size() != this->mCurrentDeviatoricCauchyStress.size() )
    NewElement.mCurrentDeviatoricCauchyStress.resize(this->mCurrentDeviatoricCauchyStress.size());

  for(unsigned int i=0; i<this->mCurrentDeviatoricCauchyStress.size(); i++)
    {
      NewElement.mCurrentDeviatoricCauchyStress[i] = this->mCurrentDeviatoricCauchyStress[i];
    }


  if ( NewElement.mUpdatedTotalCauchyStress.size() != this->mUpdatedTotalCauchyStress.size() )
    NewElement.mUpdatedTotalCauchyStress.resize(this->mUpdatedTotalCauchyStress.size());

  for(unsigned int i=0; i<this->mUpdatedTotalCauchyStress.size(); i++)
    {
      NewElement.mUpdatedTotalCauchyStress[i] = this->mUpdatedTotalCauchyStress[i];
    }


  if ( NewElement.mUpdatedDeviatoricCauchyStress.size() != this->mUpdatedDeviatoricCauchyStress.size() )
    NewElement.mUpdatedDeviatoricCauchyStress.resize(this->mUpdatedDeviatoricCauchyStress.size());

  for(unsigned int i=0; i<this->mUpdatedDeviatoricCauchyStress.size(); i++)
    {
      NewElement.mUpdatedDeviatoricCauchyStress[i] = this->mUpdatedDeviatoricCauchyStress[i];
    }

   NewElement.SetData(this->GetData());
   NewElement.SetFlags(this->GetFlags());

   return Element::Pointer( new TwoStepUpdatedLagrangianVPFluidElement(NewElement) );

}




  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPFluidElement<TDim>::Initialize()
  {
    KRATOS_TRY; 
 
    // LargeDisplacementElement::Initialize();
    const GeometryType& rGeom = this->GetGeometry();
    SizeType integration_points_number = rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_1);
    const unsigned int dimension       = rGeom.WorkingSpaceDimension();

    if ( this->mCurrentFgrad.size() != integration_points_number )
      this->mCurrentFgrad.resize( integration_points_number );

    if ( this->mUpdatedFgrad.size() != integration_points_number )
      this->mUpdatedFgrad.resize( integration_points_number );

    if ( this->mCurrentTotalCauchyStress.size() != integration_points_number )
      this->mCurrentTotalCauchyStress.resize( integration_points_number );
    
    if ( this->mCurrentDeviatoricCauchyStress.size() != integration_points_number )
      this->mCurrentDeviatoricCauchyStress.resize( integration_points_number );
    
    if ( this->mUpdatedTotalCauchyStress.size() != integration_points_number )
      this->mUpdatedTotalCauchyStress.resize( integration_points_number );
    
    if ( this->mUpdatedDeviatoricCauchyStress.size() != integration_points_number )
      this->mUpdatedDeviatoricCauchyStress.resize( integration_points_number );

    unsigned int voigtsize  = 3;
    if( TDim == 3 )
      {
        voigtsize  = 6;
      }
    for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
      {
        // this->mOldFgrad[PointNumber] = identity_matrix<double> (dimension);
        this->mCurrentFgrad[PointNumber] = identity_matrix<double> (dimension);
        this->mUpdatedFgrad[PointNumber] = identity_matrix<double> (dimension);
	this->mCurrentTotalCauchyStress[PointNumber] = ZeroVector(voigtsize);
	this->mCurrentDeviatoricCauchyStress[PointNumber] = ZeroVector(voigtsize);
	this->mUpdatedTotalCauchyStress[PointNumber] = ZeroVector(voigtsize);
	this->mUpdatedDeviatoricCauchyStress[PointNumber] = ZeroVector(voigtsize);
      }

    KRATOS_CATCH( "" );
  }
  
  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPFluidElement<TDim>::InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo)
  {

  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPFluidElement<TDim>::InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo)
  {
    KRATOS_TRY; 
    KRATOS_CATCH( "" );
  }


   template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPFluidElement<TDim>::ComputeMaterialParameters(double& Density,
									       double& DeviatoricCoeff,
									       double& VolumetricCoeff,
									       double timeStep,
									       const ShapeFunctionsType& N)
   {
     double FluidBulkModulus=0;
     double FluidViscosity=0;
     this->EvaluatePropertyFromANotRigidNode(Density,DENSITY);
     this->EvaluatePropertyFromANotRigidNode(FluidViscosity,VISCOSITY);
     this->EvaluatePropertyFromANotRigidNode(FluidBulkModulus,BULK_MODULUS);

     if(FluidBulkModulus==0){
       // std::cout<<"FluidBulkModulus was 0 !!!!!!!!"<<std::endl;
       FluidBulkModulus = 1000000000.0;
     }

     DeviatoricCoeff = FluidViscosity;
     VolumetricCoeff = FluidBulkModulus*timeStep;


  }







  template< unsigned int TDim >
  int TwoStepUpdatedLagrangianVPFluidElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo)
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
    if(NODAL_AREA.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"NODAL_AREA Key is 0. Check that the application was correctly registered.","");
    if(BDF_COEFFICIENTS.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"BDF_COEFFICIENTS Key is 0. Check that the application was correctly registered.","");
    if(DELTA_TIME.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"DELTA_TIME Key is 0. Check that the application was correctly registered.","");
    if(DYNAMIC_TAU.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"DYNAMIC_TAU Key is 0. Check that the application was correctly registered.","");

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
        if(this->GetGeometry()[i].SolutionStepsDataHas(NODAL_AREA) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing NODAL_AREA variable on solution step data for node ",this->GetGeometry()[i].Id());
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
  void TwoStepUpdatedLagrangianVPFluidElement<2>::ComputeMeanValueMaterialTangentMatrix(ElementalVariables & rElementalVariables,double& MeanValue,const ShapeFunctionDerivativesType& rDN_DX,const double secondLame,double & bulkModulus,const double Weight,double& MeanValueMass,const double TimeStep)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    MatrixType invGradDef=rElementalVariables.InvFgrad;
    double theta=0.5;
    double Count=0;

    for (SizeType j = 0; j < NumNodes; ++j)
      {
        for (SizeType i = 0; i < NumNodes; ++i)
	  {
	    double lagDNXi=rDN_DX(i,0)*invGradDef(0,0)+rDN_DX(i,1)*invGradDef(1,0);
	    double lagDNYi=rDN_DX(i,0)*invGradDef(0,1)+rDN_DX(i,1)*invGradDef(1,1);
	    double lagDNXj=rDN_DX(j,0)*invGradDef(0,0)+rDN_DX(j,1)*invGradDef(1,0);
	    double lagDNYj=rDN_DX(j,0)*invGradDef(0,1)+rDN_DX(j,1)*invGradDef(1,1);
	    // lagDNXi=rDN_DX(i,0);
	    // lagDNYi=rDN_DX(i,1);
	    // lagDNXj=rDN_DX(j,0);
	    // lagDNYj=rDN_DX(j,1);


            // First Row
            MeanValue += fabs(Weight * ( (FourThirds * secondLame + bulkModulus)*  lagDNXi * lagDNXj + lagDNYi * lagDNYj * secondLame )*theta);
            MeanValue += fabs(Weight * ( (nTwoThirds* secondLame + bulkModulus) *  lagDNXi * lagDNYj + lagDNYi * lagDNXj * secondLame )*theta);

            // Second Row
            MeanValue += fabs(Weight * ( (nTwoThirds * secondLame + bulkModulus) * lagDNYi * lagDNXj + lagDNXi * lagDNYj *  secondLame )*theta);
            MeanValue += fabs(Weight * ( (FourThirds * secondLame + bulkModulus) * lagDNYi * lagDNYj + lagDNXi * lagDNXj * secondLame )*theta);

	    Count+=4.0;

	  }

      }

    MeanValue*=1.0/Count;

    if(MeanValueMass!=0 && MeanValue!=0){
      bulkModulus*=MeanValueMass*2.0/TimeStep/MeanValue;
    }else{
      std::cout<<" DANGEROUS ELEMENT!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
      std::cout<<" MeanValueMass="<<MeanValueMass;
      std::cout<<"\t MeanValueMaterial= "<<MeanValue;
      std::cout<<"\t VolumetricCoeff= "<<bulkModulus<<std::endl;
      bulkModulus*=TimeStep;
    }
  }

  template<>
  void TwoStepUpdatedLagrangianVPFluidElement<3>::ComputeMeanValueMaterialTangentMatrix(ElementalVariables & rElementalVariables,double& MeanValue,const ShapeFunctionDerivativesType& rDN_DX,const double secondLame,double & bulkModulus,const double Weight,double& MeanValueMass,const double TimeStep)
  {

    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    MatrixType invGradDef=rElementalVariables.InvFgrad;
    double theta=0.5;
    double Count=0;
    for (SizeType j = 0; j < NumNodes; ++j)
      {
        for (SizeType i = 0; i < NumNodes; ++i)
	  {
	    double lagDNXi=rDN_DX(i,0)*invGradDef(0,0)+rDN_DX(i,1)*invGradDef(1,0)+rDN_DX(i,2)*invGradDef(2,0);
	    double lagDNYi=rDN_DX(i,0)*invGradDef(0,1)+rDN_DX(i,1)*invGradDef(1,1)+rDN_DX(i,2)*invGradDef(2,1);
	    double lagDNZi=rDN_DX(i,0)*invGradDef(0,2)+rDN_DX(i,1)*invGradDef(1,2)+rDN_DX(i,2)*invGradDef(2,2);
	    double lagDNXj=rDN_DX(j,0)*invGradDef(0,0)+rDN_DX(j,1)*invGradDef(1,0)+rDN_DX(j,2)*invGradDef(2,0);
	    double lagDNYj=rDN_DX(j,0)*invGradDef(0,1)+rDN_DX(j,1)*invGradDef(1,1)+rDN_DX(j,2)*invGradDef(2,1);
	    double lagDNZj=rDN_DX(j,0)*invGradDef(0,2)+rDN_DX(j,1)*invGradDef(1,2)+rDN_DX(j,2)*invGradDef(2,2);	  
	    // lagDNXi=rDN_DX(i,0);
	    // lagDNYi=rDN_DX(i,1);
	    // lagDNZi=rDN_DX(i,2);
	    // lagDNXj=rDN_DX(j,0);
	    // lagDNYj=rDN_DX(j,1);
	    // lagDNZj=rDN_DX(j,2);

            // First Row
            MeanValue += fabs(Weight * ( (FourThirds * secondLame + bulkModulus)*  lagDNXi * lagDNXj + (lagDNYi * lagDNYj +lagDNZi * lagDNZj) * secondLame ) *theta);
            MeanValue += fabs(Weight * ( (nTwoThirds* secondLame + bulkModulus) *  lagDNXi * lagDNYj + lagDNYi * lagDNXj * secondLame )*theta);
            MeanValue += fabs(Weight * ( (nTwoThirds* secondLame + bulkModulus) *  lagDNXi * lagDNZj + lagDNZi * lagDNXj * secondLame )*theta);

            // Second Row
            MeanValue += fabs(Weight * ( (nTwoThirds * secondLame + bulkModulus) * lagDNYi * lagDNXj + lagDNXi * lagDNYj *  secondLame )*theta);
            MeanValue += fabs(Weight * ( (FourThirds * secondLame + bulkModulus) * lagDNYi * lagDNYj + (lagDNXi * lagDNXj + lagDNZi * lagDNZj) * secondLame )*theta);
            MeanValue += fabs(Weight * ( (nTwoThirds * secondLame + bulkModulus) * lagDNYi * lagDNZj + lagDNZi * lagDNYj *  secondLame )*theta);

            // Third Row
            MeanValue += fabs(Weight * ( (nTwoThirds * secondLame + bulkModulus) * lagDNZi * lagDNXj + lagDNXi * lagDNZj *  secondLame )*theta);
            MeanValue += fabs(Weight * ( (nTwoThirds* secondLame + bulkModulus)  * lagDNZi * lagDNYj + lagDNYi * lagDNZj *  secondLame )*theta);
	    MeanValue += fabs(Weight * ( (FourThirds * secondLame + bulkModulus) * lagDNZi * lagDNZj + (lagDNXi * lagDNXj + lagDNYi * lagDNYj) * secondLame )*theta);
	    Count+=9.0;
	  }
   
      }
    MeanValue*=1.0/Count;

   if(MeanValueMass!=0 && MeanValue!=0){
      bulkModulus*=MeanValueMass*2.0/TimeStep/MeanValue;
    }else{
      std::cout<<" DANGEROUS ELEMENT!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
      std::cout<<" MeanValueMass="<<MeanValueMass;
      std::cout<<"\t MeanValueMaterial= "<<MeanValue;
      std::cout<<"\t VolumetricCoeff= "<<bulkModulus<<std::endl;
      bulkModulus*=TimeStep;
    }
  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPFluidElement<TDim>::ComputeBulkMatrixForPressureVel(Matrix& BulkVelMatrix,
										     const ShapeFunctionsType& rN,
										     const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    for (SizeType i = 0; i < NumNodes; ++i)
      {
	for (SizeType j = 0; j < NumNodes; ++j)
	  {
	    // LHS contribution
	    double Mij  = Weight*rN[i]*rN[j];
	    BulkVelMatrix(i,j) +=  Mij;
	  }

      }
  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPFluidElement<TDim>::ComputeBulkMatrixForPressureVelLump(Matrix& BulkVelMatrix,
											 const ShapeFunctionsType& rN,
											 const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    double coeff=1.0+TDim;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
    	// LHS contribution
    	double Mij  = Weight /coeff;
    	BulkVelMatrix(i,i) +=  Mij;
      }
 
  }


 
  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPFluidElement<TDim>::ComputeBulkMatrixForPressureAcc(Matrix& BulkAccMatrix,
										     const ShapeFunctionsType& rN,
										     const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    for (SizeType i = 0; i < NumNodes; ++i)
      {
    	for (SizeType j = 0; j < NumNodes; ++j)
    	  {
    	    // LHS contribution
    	    double Mij  = Weight*rN[i]*rN[j];
    	    BulkAccMatrix(i,j) +=  Mij;
    	  }

      }
  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPFluidElement<TDim>::ComputeBulkMatrixForPressureAccLump(Matrix& BulkAccMatrix,
											 const ShapeFunctionsType& rN,
											 const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    double coeff=1.0+TDim;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
    	// LHS contribution
    	double Mij  = Weight /coeff;
    	BulkAccMatrix(i,i) +=  Mij;
      }

  }

 

  template< >
  void TwoStepUpdatedLagrangianVPFluidElement<2>::ComputeBoundLHSMatrix(Matrix& BoundLHSMatrix,
									const ShapeFunctionsType& rN,
									const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    double coeff=1.0+2.0;
    double Mij  = Weight/ coeff;
    GeometryType& rGeom = this->GetGeometry();

      if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)){
	BoundLHSMatrix(0,0) +=  Mij;
	BoundLHSMatrix(1,1) +=  Mij;
      }else if(rGeom[0].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)){
	BoundLHSMatrix(0,0) +=  Mij;
	BoundLHSMatrix(2,2) +=  Mij;
      }else if(rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)){
	BoundLHSMatrix(1,1) +=  Mij;
	BoundLHSMatrix(2,2) +=  Mij;
      }
  }

 template<  >
  void TwoStepUpdatedLagrangianVPFluidElement<3>::ComputeBoundLHSMatrix(Matrix& BoundLHSMatrix,
									   const ShapeFunctionsType& rN,
									   const double Weight)
 {
   const SizeType NumNodes = this->GetGeometry().PointsNumber();
   double coeff=1.0+3.0;
   double Mij  = Weight/ coeff;
   GeometryType& rGeom = this->GetGeometry();

     if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)){
       BoundLHSMatrix(0,0) +=  Mij;
       BoundLHSMatrix(1,1) +=  Mij;
       BoundLHSMatrix(2,2) +=  Mij;
     }
     if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
       BoundLHSMatrix(0,0) +=  Mij;
       BoundLHSMatrix(1,1) +=  Mij;
       BoundLHSMatrix(3,3) +=  Mij;
     }
     if(rGeom[0].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
       BoundLHSMatrix(0,0) +=  Mij;
       BoundLHSMatrix(2,2) +=  Mij;
       BoundLHSMatrix(3,3) +=  Mij;
     }
     if(rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
       BoundLHSMatrix(1,1) +=  Mij;
       BoundLHSMatrix(2,2) +=  Mij;
       BoundLHSMatrix(3,3) +=  Mij;
     }

 }
 


template< >
void TwoStepUpdatedLagrangianVPFluidElement<2>::ComputeBoundRHSVector(VectorType& BoundRHSVector,
								      const ShapeFunctionsType& rN,
								      const double TimeStep,
								      const double BoundRHSCoeffAcc,
								      const double BoundRHSCoeffDev)
{
  const SizeType NumNodes = this->GetGeometry().PointsNumber();
  double coeff=1.0+2.0;
  GeometryType& rGeom = this->GetGeometry();
  array_1d<double, 3>  NormalA(3,0.0);
  NormalA.clear();
  array_1d<double, 3>  NormalB(3,0.0);
  NormalB.clear();
  array_1d<double, 3>  AccA(3,0.0);
  AccA.clear();
  array_1d<double, 3>  AccB(3,0.0);
  AccB.clear();

  VectorType ElementalAcceleration = ZeroVector(2);
  this->GetElementalAcceleration(ElementalAcceleration,0,TimeStep);

  if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE) ){
    AccA= 0.5/TimeStep*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
    AccB= 0.5/TimeStep*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
    NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
    NormalB    = rGeom[1].FastGetSolutionStepValue(NORMAL);
    BoundRHSVector[0] += (BoundRHSCoeffAcc*(AccA[0]*NormalA[0]+AccA[1]*NormalA[1]) + BoundRHSCoeffDev)/coeff ;
    BoundRHSVector[1] += (BoundRHSCoeffAcc*(AccB[0]*NormalB[0]+AccB[1]*NormalB[1]) + BoundRHSCoeffDev)/coeff ;
  }else  if(rGeom[0].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE) ){
    AccA= 0.5/TimeStep*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
    AccB= 0.5/TimeStep*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
    NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
    NormalB    = rGeom[2].FastGetSolutionStepValue(NORMAL);
    BoundRHSVector[0] += (BoundRHSCoeffAcc*(AccA[0]*NormalA[0]+AccA[1]*NormalA[1]) + BoundRHSCoeffDev)/coeff ;
    BoundRHSVector[2] += (BoundRHSCoeffAcc*(AccB[0]*NormalB[0]+AccB[1]*NormalB[1]) + BoundRHSCoeffDev)/coeff ;
  }else  if(rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE) ){
    AccA= 0.5/TimeStep*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
    AccB= 0.5/TimeStep*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
    NormalA    = rGeom[1].FastGetSolutionStepValue(NORMAL);
    NormalB    = rGeom[2].FastGetSolutionStepValue(NORMAL);
    BoundRHSVector[1] += (BoundRHSCoeffAcc*(AccA[0]*NormalA[0]+AccA[1]*NormalA[1]) + BoundRHSCoeffDev)/coeff ;
    BoundRHSVector[2] += (BoundRHSCoeffAcc*(AccB[0]*NormalB[0]+AccB[1]*NormalB[1]) + BoundRHSCoeffDev)/coeff ;
  }


}


template< >
void TwoStepUpdatedLagrangianVPFluidElement<3>::ComputeBoundRHSVector(VectorType& BoundRHSVector,
								      const ShapeFunctionsType& rN,
								      const double TimeStep,
								      const double BoundRHSCoeffAcc,
								      const double BoundRHSCoeffDev)
{
  const SizeType NumNodes = this->GetGeometry().PointsNumber();
  double coeff=1.0+3.0;
  array_1d<double, 3>  NormalA(3,0.0);
  NormalA.clear();
  array_1d<double, 3>  NormalB(3,0.0);
  NormalB.clear();
  array_1d<double, 3>  NormalC(3,0.0);
  NormalC.clear();
  array_1d<double, 3>  AccA(3,0.0);
  AccA.clear();
  array_1d<double, 3>  AccB(3,0.0);
  AccB.clear();
  array_1d<double, 3>  AccC(3,0.0);
  AccC.clear();
  GeometryType& rGeom = this->GetGeometry();

  VectorType ElementalAcceleration = ZeroVector(2);
  this->GetElementalAcceleration(ElementalAcceleration,0,TimeStep);

  if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)){
    AccA= 0.5/TimeStep*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
    AccB= 0.5/TimeStep*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
    AccC= 0.5/TimeStep*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
    NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
    NormalB    = rGeom[1].FastGetSolutionStepValue(NORMAL);
    NormalC    = rGeom[2].FastGetSolutionStepValue(NORMAL);
    BoundRHSVector[0] += (BoundRHSCoeffAcc*(AccA[0]*NormalA[0] + AccA[1]*NormalA[1] + AccA[2]*NormalA[2]) +
			  BoundRHSCoeffDev)/coeff;
    BoundRHSVector[1] += (BoundRHSCoeffAcc*(AccB[0]*NormalB[0] + AccB[1]*NormalB[1] + AccB[2]*NormalB[2]) +
			  BoundRHSCoeffDev)/coeff;
    BoundRHSVector[2] += (BoundRHSCoeffAcc*(AccC[0]*NormalC[0] + AccC[1]*NormalC[1] + AccC[2]*NormalC[2]) +
			  BoundRHSCoeffDev)/coeff;
  }else if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
    AccA= 0.5/TimeStep*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
    AccB= 0.5/TimeStep*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
    AccC= 0.5/TimeStep*(rGeom[3].FastGetSolutionStepValue(VELOCITY,0)-rGeom[3].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION,1); 
    NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
    NormalB    = rGeom[1].FastGetSolutionStepValue(NORMAL);
    NormalC    = rGeom[3].FastGetSolutionStepValue(NORMAL);
    BoundRHSVector[0] += (BoundRHSCoeffAcc*(AccA[0]*NormalA[0] + AccA[1]*NormalA[1] + AccA[2]*NormalA[2]) +
			  BoundRHSCoeffDev)/coeff;
    BoundRHSVector[1] += (BoundRHSCoeffAcc*(AccB[0]*NormalB[0] + AccB[1]*NormalB[1] + AccB[2]*NormalB[2]) +
			  BoundRHSCoeffDev)/coeff;
    BoundRHSVector[3] += (BoundRHSCoeffAcc*(AccC[0]*NormalC[0] + AccC[1]*NormalC[1] + AccC[2]*NormalC[2]) +
			  BoundRHSCoeffDev)/coeff;
  }else if(rGeom[0].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
    AccA= 0.5/TimeStep*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
    AccB= 0.5/TimeStep*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
    AccC= 0.5/TimeStep*(rGeom[3].FastGetSolutionStepValue(VELOCITY,0)-rGeom[3].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION,1); 
    NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
    NormalB    = rGeom[2].FastGetSolutionStepValue(NORMAL);
    NormalC    = rGeom[3].FastGetSolutionStepValue(NORMAL);
    BoundRHSVector[0] += (BoundRHSCoeffAcc*(AccA[0]*NormalA[0] + AccA[1]*NormalA[1] + AccA[2]*NormalA[2]) +
			  BoundRHSCoeffDev)/coeff;
    BoundRHSVector[2] += (BoundRHSCoeffAcc*(AccB[0]*NormalB[0] + AccB[1]*NormalB[1] + AccB[2]*NormalB[2]) +
			  BoundRHSCoeffDev)/coeff;
    BoundRHSVector[3] += (BoundRHSCoeffAcc*(AccC[0]*NormalC[0] + AccC[1]*NormalC[1] + AccC[2]*NormalC[2]) +
			  BoundRHSCoeffDev)/coeff;
  }else if(rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){      
    AccA= 0.5/TimeStep*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
    AccB= 0.5/TimeStep*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
    AccC= 0.5/TimeStep*(rGeom[3].FastGetSolutionStepValue(VELOCITY,0)-rGeom[3].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION,1); 
    NormalA    = rGeom[1].FastGetSolutionStepValue(NORMAL);
    NormalB    = rGeom[2].FastGetSolutionStepValue(NORMAL);
    NormalC    = rGeom[3].FastGetSolutionStepValue(NORMAL);
    BoundRHSVector[1] += (BoundRHSCoeffAcc*(AccA[0]*NormalA[0] + AccA[1]*NormalA[1] + AccA[2]*NormalA[2]) +
			  BoundRHSCoeffDev)/coeff;
    BoundRHSVector[2] += (BoundRHSCoeffAcc*(AccB[0]*NormalB[0] + AccB[1]*NormalB[1] + AccB[2]*NormalB[2]) +
			  BoundRHSCoeffDev)/coeff;
    BoundRHSVector[3] += (BoundRHSCoeffAcc*(AccC[0]*NormalC[0] + AccC[1]*NormalC[1] + AccC[2]*NormalC[2]) +
			  BoundRHSCoeffDev)/coeff;
  }



}



  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPFluidElement<TDim>::CalculateTauFIC(double& Tau,
								     double ElemSize,
								     const array_1d< double, 3 > & rAdvVel,
								     const double Density,
								     const double Viscosity,
								     const ProcessInfo& rCurrentProcessInfo)
  {
    double DeltaTime = rCurrentProcessInfo.GetValue(DELTA_TIME);

    double MeanVelocity=0;
    this->CalcMeanVelocity(MeanVelocity,0);

      Tau = 1.0 / (2.0 * Density *(0.5 * MeanVelocity / ElemSize + 0.5/DeltaTime) +  8.0 * Viscosity / (ElemSize * ElemSize) ); 
  
    if(Tau<0.0000001)
      Tau=0.0000001;

    
  }

template< unsigned int TDim >
 void TwoStepUpdatedLagrangianVPFluidElement<TDim>::AddStabilizationMatrixLHS(MatrixType& rLeftHandSideMatrix,
									      Matrix& BulkAccMatrix,
									      const ShapeFunctionsType& rN,
									      const double Weight)
 {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

   if( BulkAccMatrix.size1() != NumNodes )
     BulkAccMatrix.resize(NumNodes,NumNodes);

   BulkAccMatrix = ZeroMatrix(NumNodes,NumNodes);
    for (SizeType i = 0; i < NumNodes; ++i)
      {
	// LHS contribution
	for (SizeType j = 0; j < NumNodes; ++j)
	  {
	    double Mij = 0.0;
	    Mij = Weight * rN[i] * rN[j];
	    BulkAccMatrix(i,j) +=  Mij;
	  }
      }
    rLeftHandSideMatrix+=BulkAccMatrix;

  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPFluidElement<TDim>::ComputeStabLaplacianMatrix(MatrixType& StabLaplacianMatrix,
										const ShapeFunctionDerivativesType& rDN_DX,
										const double Weight)
										
  {
    // LHS contribution
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    for (SizeType i = 0; i < NumNodes; ++i)
      {
	for (SizeType j = 0; j < NumNodes; ++j)
	  {
	    double Lij = 0.0;
	    for (SizeType d = 0; d < TDim; ++d){
	      Lij += rDN_DX(i,d) * rDN_DX(j,d);
	    }
	    StabLaplacianMatrix(i,j) += Weight * Lij ;
	  }
      }
  }



  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPFluidElement<TDim>::AddStabilizationNodalTermsLHS(MatrixType& rLeftHandSideMatrix,
										   const double Tau,
										   const double Weight,
										   const ShapeFunctionDerivativesType& rDN_DX,
										   const SizeType i)
  {
    // LHS contribution
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    for (SizeType j = 0; j < NumNodes; ++j)
      {
        double Lij = 0.0;
        for (SizeType d = 0; d < TDim; ++d){
    	  Lij += rDN_DX(i,d) * rDN_DX(j,d);
    	}
        Lij *= Tau;

        rLeftHandSideMatrix(i,j) += Weight * Lij;
      }
  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPFluidElement<TDim>::AddStabilizationNodalTermsRHS(VectorType& rRightHandSideVector,
										   const double Tau,
										   const double Density,
										   const array_1d<double,3> BodyForce,
										   const double Weight,
										   const ShapeFunctionDerivativesType& rDN_DX,
										   const SizeType i)
  {

    // // Evaluate the pressure and pressure gradient at this point (for the G * P_n term)
    array_1d<double,TDim> OldPressureGradient(TDim,0.0);
    this->EvaluateGradientInPoint(OldPressureGradient,PRESSURE,rDN_DX);
    double RHSi = 0;
    if( this->GetGeometry()[i].SolutionStepsDataHas(VOLUME_ACCELERATION) ){ // it must be checked once at the begining only
      array_1d<double, 3 >& VolumeAcceleration = this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);
      for (SizeType d = 0; d < TDim; ++d)
	{
	  RHSi += - rDN_DX(i,d) * Tau * ( Density * VolumeAcceleration[d] );
	}
    }
    rRightHandSideVector[i] += Weight * RHSi;

  }



template< unsigned int TDim>
bool TwoStepUpdatedLagrangianVPFluidElement<TDim>::CalcMechanicsUpdated(ElementalVariables & rElementalVariables,
									const ProcessInfo& rCurrentProcessInfo,
									unsigned int g,
									const ShapeFunctionsType& N)
{
  bool computeElement=false;
  double theta=0.5;
  computeElement=this->CalcStrainRate(rElementalVariables,rCurrentProcessInfo,g,theta);
  const double TimeStep=0.5/rCurrentProcessInfo[BDF_COEFFICIENTS][2];
  this->CalcElasticPlasticCauchySplitted(rElementalVariables,TimeStep,g,N);
  return computeElement;
} 


  template<>
  void TwoStepUpdatedLagrangianVPFluidElement<2>::GetPositions(Vector& rValues,const ProcessInfo& rCurrentProcessInfo,const double theta)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 2*NumNodes;

    if (rValues.size() != LocalSize) rValues.resize(LocalSize);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
	rValues[Index++] = rGeom[i].X();
	rValues[Index++] = rGeom[i].Y();
      }
  }



  template<>
  void TwoStepUpdatedLagrangianVPFluidElement<3>::GetPositions(Vector& rValues,const ProcessInfo& rCurrentProcessInfo,const double theta)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3*NumNodes;

    if (rValues.size() != LocalSize) rValues.resize(LocalSize);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
 	rValues[Index++] = rGeom[i].X();
        rValues[Index++] = rGeom[i].Y();
        rValues[Index++] = rGeom[i].Z();
      }
  }



template <  unsigned int TDim> 
void TwoStepUpdatedLagrangianVPFluidElement<TDim>:: InitializeElementalVariables(ElementalVariables & rElementalVariables)
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
void TwoStepUpdatedLagrangianVPFluidElement<2>:: CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables,double TimeStep, unsigned int g,const ShapeFunctionsType& N)
{

  rElementalVariables.CurrentTotalCauchyStress=mCurrentTotalCauchyStress[g];
  rElementalVariables.CurrentDeviatoricCauchyStress=mCurrentDeviatoricCauchyStress[g];

  double Density  = 0;
  double CurrSecondLame  = 0;
  double CurrBulkModulus = 0;

  this->ComputeMaterialParameters(Density,CurrSecondLame,CurrBulkModulus,TimeStep,N);
 
  double CurrFirstLame  = 0;
  CurrFirstLame  =CurrBulkModulus - 2.0*CurrSecondLame/3.0;

  double DefX=rElementalVariables.SpatialDefRate[0];
  double DefY=rElementalVariables.SpatialDefRate[1];
  double DefXY=rElementalVariables.SpatialDefRate[2];

  double DefVol=rElementalVariables.VolumetricDefRate;

  double sigmaDev_xx= 2*CurrSecondLame*(DefX - DefVol/3.0);
  double sigmaDev_yy= 2*CurrSecondLame*(DefY - DefVol/3.0);
  double sigmaDev_xy= 2*CurrSecondLame*DefXY;

  double sigmaTot_xx= CurrFirstLame*DefVol + 2.0*CurrSecondLame*DefX;
  double sigmaTot_yy= CurrFirstLame*DefVol + 2.0*CurrSecondLame*DefY;
  double sigmaTot_xy= 2.0*CurrSecondLame*DefXY;

  // sigmaDev_xx=rElementalVariables.CurrentDeviatoricCauchyStress[0];
  // sigmaDev_yy=rElementalVariables.CurrentDeviatoricCauchyStress[1];
  // sigmaDev_xy=rElementalVariables.CurrentDeviatoricCauchyStress[2];

  // sigmaTot_xx+=rElementalVariables.CurrentTotalCauchyStress[0];
  // sigmaTot_yy+=rElementalVariables.CurrentTotalCauchyStress[1];
  // sigmaTot_xy+=rElementalVariables.CurrentTotalCauchyStress[2];

  sigmaTot_xx= sigmaDev_xx + rElementalVariables.MeanPressure;
  sigmaTot_yy= sigmaDev_yy + rElementalVariables.MeanPressure;
  sigmaTot_xy= sigmaDev_xy;

  // sigmaDev_xx= sigmaTot_xx - rElementalVariables.MeanPressure;
  // sigmaDev_yy= sigmaTot_yy - rElementalVariables.MeanPressure;
  // sigmaDev_xy= sigmaTot_xy;

  rElementalVariables.UpdatedDeviatoricCauchyStress[0]=sigmaDev_xx;
  rElementalVariables.UpdatedDeviatoricCauchyStress[1]=sigmaDev_yy;
  rElementalVariables.UpdatedDeviatoricCauchyStress[2]=sigmaDev_xy;

  rElementalVariables.UpdatedTotalCauchyStress[0]=sigmaTot_xx;
  rElementalVariables.UpdatedTotalCauchyStress[1]=sigmaTot_yy;
  rElementalVariables.UpdatedTotalCauchyStress[2]=sigmaTot_xy;

  // this->mCurrentTotalCauchyStress[g]=rElementalVariables.CurrentTotalCauchyStress;
  this->mUpdatedTotalCauchyStress[g]=rElementalVariables.UpdatedTotalCauchyStress;
  // this->mCurrentDeviatoricCauchyStress[g]=rElementalVariables.CurrentDeviatoricCauchyStress;
  this->mUpdatedDeviatoricCauchyStress[g]=rElementalVariables.UpdatedDeviatoricCauchyStress;

}


template < > 
void TwoStepUpdatedLagrangianVPFluidElement<3>:: CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables, double TimeStep, unsigned int g,const ShapeFunctionsType& N)
{


  rElementalVariables.CurrentTotalCauchyStress=mCurrentTotalCauchyStress[g];
  rElementalVariables.CurrentDeviatoricCauchyStress=mCurrentDeviatoricCauchyStress[g];

  double Density  = 0;
  double CurrSecondLame  = 0;
  double CurrBulkModulus = 0;

  this->ComputeMaterialParameters(Density,CurrSecondLame,CurrBulkModulus,TimeStep,N);
 
  double CurrFirstLame  = 0;
  CurrFirstLame  =CurrBulkModulus - 2.0*CurrSecondLame/3.0;

   
  double DefX=rElementalVariables.SpatialDefRate[0];
  double DefY=rElementalVariables.SpatialDefRate[1];
  double DefZ=rElementalVariables.SpatialDefRate[2];
  double DefXY=rElementalVariables.SpatialDefRate[3];
  double DefXZ=rElementalVariables.SpatialDefRate[4];
  double DefYZ=rElementalVariables.SpatialDefRate[5];

  double DefVol=rElementalVariables.VolumetricDefRate;

  double sigmaDev_xx= 2*CurrSecondLame*(DefX - DefVol/3.0);
  double sigmaDev_yy= 2*CurrSecondLame*(DefY - DefVol/3.0);
  double sigmaDev_zz= 2*CurrSecondLame*(DefZ - DefVol/3.0);
  double sigmaDev_xy= 2*CurrSecondLame*DefXY;
  double sigmaDev_xz= 2*CurrSecondLame*DefXZ;
  double sigmaDev_yz= 2*CurrSecondLame*DefYZ;

  double sigmaTot_xx= CurrFirstLame*DefVol + 2*CurrSecondLame*DefX;
  double sigmaTot_yy= CurrFirstLame*DefVol + 2*CurrSecondLame*DefY;
  double sigmaTot_zz= CurrFirstLame*DefVol + 2*CurrSecondLame*DefZ;
  double sigmaTot_xy= 2*CurrSecondLame*DefXY;
  double sigmaTot_xz= 2*CurrSecondLame*DefXZ;
  double sigmaTot_yz= 2*CurrSecondLame*DefYZ;


  // sigmaDev_xx+=rElementalVariables.CurrentDeviatoricCauchyStress[0];
  // sigmaDev_yy+=rElementalVariables.CurrentDeviatoricCauchyStress[1];
  // sigmaDev_zz+=rElementalVariables.CurrentDeviatoricCauchyStress[2];
  // sigmaDev_xy+=rElementalVariables.CurrentDeviatoricCauchyStress[3];
  // sigmaDev_xz+=rElementalVariables.CurrentDeviatoricCauchyStress[4];
  // sigmaDev_yz+=rElementalVariables.CurrentDeviatoricCauchyStress[5];

  sigmaTot_xx= sigmaDev_xx + rElementalVariables.MeanPressure;
  sigmaTot_yy= sigmaDev_yy + rElementalVariables.MeanPressure;
  sigmaTot_zz= sigmaDev_zz + rElementalVariables.MeanPressure;
  sigmaTot_xy= sigmaDev_xy;
  sigmaTot_xz= sigmaDev_xz;
  sigmaTot_yz= sigmaDev_yz;

  // sigmaTot_xx+=rElementalVariables.CurrentTotalCauchyStress[0];
  // sigmaTot_yy+=rElementalVariables.CurrentTotalCauchyStress[1];
  // sigmaTot_zz+=rElementalVariables.CurrentTotalCauchyStress[2];
  // sigmaTot_xy+=rElementalVariables.CurrentTotalCauchyStress[3];
  // sigmaTot_xz+=rElementalVariables.CurrentTotalCauchyStress[4];
  // sigmaTot_yz+=rElementalVariables.CurrentTotalCauchyStress[5];

  // sigmaDev_xx= sigmaTot_xx - rElementalVariables.MeanPressure;
  // sigmaDev_yy= sigmaTot_yy - rElementalVariables.MeanPressure;
  // sigmaDev_zz= sigmaTot_zz - rElementalVariables.MeanPressure;
  // sigmaDev_xy= sigmaTot_xy;
  // sigmaDev_xz= sigmaTot_xz;
  // sigmaDev_yz= sigmaTot_yz;


  rElementalVariables.UpdatedDeviatoricCauchyStress[0]=sigmaDev_xx;
  rElementalVariables.UpdatedDeviatoricCauchyStress[1]=sigmaDev_yy;
  rElementalVariables.UpdatedDeviatoricCauchyStress[2]=sigmaDev_zz;
  rElementalVariables.UpdatedDeviatoricCauchyStress[3]=sigmaDev_xy;
  rElementalVariables.UpdatedDeviatoricCauchyStress[4]=sigmaDev_xz;
  rElementalVariables.UpdatedDeviatoricCauchyStress[5]=sigmaDev_yz;

  rElementalVariables.UpdatedTotalCauchyStress[0]=sigmaTot_xx;
  rElementalVariables.UpdatedTotalCauchyStress[1]=sigmaTot_yy;
  rElementalVariables.UpdatedTotalCauchyStress[2]=sigmaTot_zz;
  rElementalVariables.UpdatedTotalCauchyStress[3]=sigmaTot_xy;
  rElementalVariables.UpdatedTotalCauchyStress[4]=sigmaTot_xz;
  rElementalVariables.UpdatedTotalCauchyStress[5]=sigmaTot_yz;

  // this->mCurrentTotalCauchyStress[g]=rElementalVariables.CurrentTotalCauchyStress;
  this->mUpdatedTotalCauchyStress[g]=rElementalVariables.UpdatedTotalCauchyStress;
  // this->mCurrentDeviatoricCauchyStress[g]=rElementalVariables.CurrentDeviatoricCauchyStress;
  this->mUpdatedDeviatoricCauchyStress[g]=rElementalVariables.UpdatedDeviatoricCauchyStress;

}


template < unsigned int TDim > 
void TwoStepUpdatedLagrangianVPFluidElement<TDim>:: UpdateCauchyStress(unsigned int g,ProcessInfo& rCurrentProcessInfo)
 {
   this->mCurrentTotalCauchyStress[g]=this->mUpdatedTotalCauchyStress[g];
   this->mCurrentDeviatoricCauchyStress[g]=this->mUpdatedDeviatoricCauchyStress[g];
 }


    template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPFluidElement<TDim>::CalculateLocalContinuityEqForPressure(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {


    GeometryType& rGeom = this->GetGeometry();
    // const SizeType NumNodes = rGeom.PointsNumber();
    const unsigned int NumNodes = rGeom.PointsNumber();

    // Check sizes and initialize
    if( rLeftHandSideMatrix.size1() != NumNodes ) 
      rLeftHandSideMatrix.resize(NumNodes,NumNodes);

    rLeftHandSideMatrix = ZeroMatrix(NumNodes,NumNodes);
 
    if( rRightHandSideVector.size() != NumNodes )
      rRightHandSideVector.resize(NumNodes);

    rRightHandSideVector = ZeroVector(NumNodes);
 

    bool computeElement=false;
    // computeElement=CheckSliverElements();
     computeElement=true;

    if(computeElement==true){
      // Shape functions and integration points
      ShapeFunctionDerivativesArrayType DN_DX;
      Matrix NContainer;
      VectorType GaussWeights;
      this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);
      const unsigned int NumGauss = GaussWeights.size();

      // Stabilization parameters
      double ElemSize = this->ElementSize();
      double Tau=0;

      MatrixType BulkVelMatrix = ZeroMatrix(NumNodes,NumNodes);
      MatrixType BulkAccMatrix = ZeroMatrix(NumNodes,NumNodes);
      MatrixType BulkVelMatrixLump = ZeroMatrix(NumNodes,NumNodes);
      MatrixType BulkAccMatrixLump = ZeroMatrix(NumNodes,NumNodes);
      MatrixType BoundLHSMatrix = ZeroMatrix(NumNodes,NumNodes);
      VectorType BoundRHSVector = ZeroVector(NumNodes);
      MatrixType StabLaplacianMatrix = ZeroMatrix(NumNodes,NumNodes);

      //   break;
      // }


      const double TimeStep=0.5/rCurrentProcessInfo[BDF_COEFFICIENTS][2];

      double Density = 0;
      double DeviatoricCoeff = 0;
      double VolumetricCoeff = 0;


      ElementalVariables rElementalVariables;
      this->InitializeElementalVariables(rElementalVariables);

      // Loop on integration points
      for (unsigned int g = 0; g < NumGauss; g++)
	{

	  double theta=0;
	  computeElement=this->CalcStrainRate(rElementalVariables,rCurrentProcessInfo,g,theta);
	  if(computeElement==true){
	    // this->UpdateCauchyStress(g);
	    // const double GaussWeight = fabs(GaussWeights[g]);
	    const double GaussWeight = GaussWeights[g];
	    const ShapeFunctionsType& N = row(NContainer,g);
	    const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];
	    this->ComputeMaterialParameters(Density,DeviatoricCoeff,VolumetricCoeff,TimeStep,N);

	    // Evaluate required variables at the integration point
	    array_1d<double,3> BodyForce(3,0.0);

	    // this->EvaluateInPoint(Density,DENSITY,N);

	    this->EvaluateInPoint(BodyForce,BODY_FORCE,N);

	    // Stabilization parameters
	    array_1d<double,3> ConvVel(3,0.0);
	    // this->EvaluateConvVelocity(ConvVel,N);
	    double Viscosity = 0;
	    this->EvaluateInPoint(Viscosity,VISCOSITY,N);


	    this->CalculateTauFIC(Tau,ElemSize,ConvVel,Density,DeviatoricCoeff,rCurrentProcessInfo);

	    double BulkCoeff =GaussWeight/(VolumetricCoeff);
	    if(rCurrentProcessInfo[STEP]>-1){
	      this->ComputeBulkMatrixForPressureVel(BulkVelMatrix,N,BulkCoeff);
	      this->ComputeBulkMatrixForPressureVelLump(BulkVelMatrixLump,N,BulkCoeff);
	      rLeftHandSideMatrix+=BulkVelMatrixLump;	
	    }


	    double BulkStabCoeff=BulkCoeff*Tau*Density/TimeStep;
	    if(rCurrentProcessInfo[STEP]>-1){
	      this->ComputeBulkMatrixForPressureAcc(BulkAccMatrix,N,BulkStabCoeff);
	      this->ComputeBulkMatrixForPressureAccLump(BulkAccMatrixLump,N,BulkStabCoeff);
	      rLeftHandSideMatrix+=BulkAccMatrixLump;
	    
	    }
	    double BoundLHSCoeff=Tau*ElemSize*ElemSize/GaussWeight;
	    if(TDim==3){
	      double surface3D=GaussWeight/(0.81649658*ElemSize);
	      BoundLHSCoeff=Tau*2*surface3D/ElemSize;
	    }
	    else{
	      double line2D=2*GaussWeight/ElemSize;
	      BoundLHSCoeff=Tau*2*line2D/ElemSize;
	    }
	    this->ComputeBoundLHSMatrix(BoundLHSMatrix,N,BoundLHSCoeff);
	    rLeftHandSideMatrix+=BoundLHSMatrix;

	    double NProjSpatialDefRate=0;
	    double NAcc=0;
	    this->CalcNormalProjectionsForBoundRHSVector(rElementalVariables.SpatialDefRate,NAcc,NProjSpatialDefRate,TimeStep);
	    double BoundRHSCoeffAcc=Tau*ElemSize*Density;
	    double BoundRHSCoeffDev=BoundLHSCoeff*2.0*NProjSpatialDefRate*DeviatoricCoeff;

	    if(TDim==3){
	      double surface3D=GaussWeight/(0.81649658*ElemSize);
	      BoundRHSCoeffAcc=Tau*surface3D*Density;
	      BoundRHSCoeffDev=Tau*surface3D*4.0*NProjSpatialDefRate*DeviatoricCoeff/ElemSize;
	    }
	    else{
	      double line2D=2*GaussWeight/ElemSize;
	      BoundRHSCoeffAcc=Tau*line2D*Density;
	      BoundRHSCoeffDev=Tau*line2D*4.0*NProjSpatialDefRate*DeviatoricCoeff/ElemSize;
	    }
	    this->ComputeBoundRHSVector(BoundRHSVector,N,TimeStep,BoundRHSCoeffAcc,BoundRHSCoeffDev);

	    VectorType UpdatedPressure;
	    VectorType CurrentPressure;
	    VectorType PreviousPressure;

	    UpdatedPressure = ZeroVector(NumNodes);
	    CurrentPressure = ZeroVector(NumNodes);
	    PreviousPressure = ZeroVector(NumNodes);


	    this->GetPressureValues(UpdatedPressure,0);
	    this->GetPressureValues(CurrentPressure,1);
	    this->GetPressureValues(PreviousPressure,2);

	    VectorType DeltaPressure = UpdatedPressure-CurrentPressure;
	    VectorType PreviousDeltaPressure = CurrentPressure-PreviousPressure;
	    VectorType DeltaPressureVariation = DeltaPressure - PreviousDeltaPressure;

	    // rRightHandSideVector -= prod(BulkVelMatrixLump,DeltaPressure);
	    rRightHandSideVector -=prod(BulkVelMatrixLump,UpdatedPressure);
	    rRightHandSideVector -=-prod(BulkVelMatrix,CurrentPressure);
	    // VectorType PressureAccStabilizationRHSterm = ZeroVector(NumNodes);
	    // PressureAccStabilizationRHSterm =prod(BulkAccMatrixLump,DeltaPressureVariation);
	    // PressureAccStabilizationRHSterm =prod(BulkAccMatrixLump,UpdatedPressure);
	    // PressureAccStabilizationRHSterm -=prod(BulkAccMatrix,CurrentPressure);
	
	    // rRightHandSideVector -= prod(BulkAccMatrixLump,DeltaPressureVariation);
	    rRightHandSideVector -=prod(BulkAccMatrixLump,DeltaPressure);
	    rRightHandSideVector -=-prod(BulkAccMatrix,PreviousDeltaPressure);

	    double DivU=0;
	    this->EvaluateDivergenceInPoint(DivU,VELOCITY,rDN_DX);
	    // std::cout<<"DivU "<<DivU;
	    DivU=rElementalVariables.VolumetricDefRate;
	    // std::cout<<"  VolDef "<<DivU<<std::endl;

	    double StabLaplacianWeight=Tau*GaussWeight;
	    this->ComputeStabLaplacianMatrix(StabLaplacianMatrix,rDN_DX,StabLaplacianWeight);//laplacian
	    rLeftHandSideMatrix+=StabLaplacianMatrix;

	    rRightHandSideVector -= prod(StabLaplacianMatrix,UpdatedPressure);

	    rRightHandSideVector -= prod(BoundLHSMatrix,UpdatedPressure);


	    // Add convection, stabilization and RHS contributions to the local system equation
	    for (SizeType i = 0; i < NumNodes; ++i)
	      {
		// this->AddStabilizationNodalTermsLHS(rLeftHandSideMatrix,Tau,GaussWeight,rDN_DX,i);//laplacian
         
		// RHS contribution
		// Velocity divergence
		double RHSi =  N[i] * DivU;
		rRightHandSideVector[i] += GaussWeight * RHSi;
		this->AddStabilizationNodalTermsRHS(rRightHandSideVector,Tau,Density,BodyForce,GaussWeight,rDN_DX,i);

		rRightHandSideVector[i] += BoundRHSVector[i];
		// std::cout<<"  B[i] "<<BoundRHSVector[i]<<" " <<rGeom[i].Coordinates()<<"\n ";
		// rRightHandSideVector[i] += PressureAccStabilizationRHSterm[i];
		// if(rGeom[i].Is(FREE_SURFACE) || rGeom[i].Is(RIGID)){
		//   rRightHandSideVector[i]=0;
		// }
	      }

	  }else
	    {
	      // const double GaussWeight = fabs(GaussWeights[g]);
	      const double GaussWeight = GaussWeights[g];
	      const ShapeFunctionsType& N = row(NContainer,g);
	      this->ComputeMaterialParameters(Density,DeviatoricCoeff,VolumetricCoeff,TimeStep,N);
	      double BulkCoeff =GaussWeight/(VolumetricCoeff);
	      // this->ComputeBulkMatrixForPressureVel(BulkVelMatrix,N,BulkCoeff);
	      this->ComputeBulkMatrixForPressureVelLump(BulkVelMatrixLump,N,BulkCoeff);
	      rLeftHandSideMatrix+=BulkVelMatrixLump;
	      for (SizeType n = 0; n < NumNodes; ++n)
		{
		  rGeom[n].Set(INTERFACE); //I set interface in order to not compute it in the continuity equation and the next non-linear iterations
		}


	  }

	}   
   

    }else{
      unsigned int numBoundary=0; 
      for (SizeType n = 0; n < NumNodes; ++n)
	{
	  if(rGeom[n].Is(BOUNDARY)){
	    numBoundary++;
	  } 
	}
      if(numBoundary==NumNodes){
	rLeftHandSideMatrix = IdentityMatrix(NumNodes,NumNodes);
	for (SizeType n = 0; n < NumNodes; ++n)
	  {
	    rGeom[n].FastGetSolutionStepValue(PRESSURE,0)=0;
	    rGeom[n].FastGetSolutionStepValue(PRESSURE,1)=0;
	  }
	std::cout<<"DO NOT COMPUTE THIS BOUNDARY ELEMENT: IMPOSE PRESSURE=0 "<<rGeom[0].Coordinates()<<std::endl;
      }
    }
  }
  




  template class TwoStepUpdatedLagrangianVPFluidElement<2>;
  template class TwoStepUpdatedLagrangianVPFluidElement<3>;

}
