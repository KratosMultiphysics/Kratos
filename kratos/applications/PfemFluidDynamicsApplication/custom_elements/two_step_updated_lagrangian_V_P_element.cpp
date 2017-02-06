//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 0.0 $
//
//   Implementation of the Gauss-Seidel two step Updated Lagrangian Velocity-Pressure element
//     ( There is a ScalingConstant to multiply the mass balance equation for a number because i read it somewhere)
//

// System includes

// External includes
 
// Project includes
#include "custom_elements/two_step_updated_lagrangian_V_P_element.h"
#include "includes/cfd_variables.h"

namespace Kratos {

  /*
   * public TwoStepUpdatedLagrangianVPElement<TDim> functions
   */


  template< unsigned int TDim >
  TwoStepUpdatedLagrangianVPElement<TDim>::TwoStepUpdatedLagrangianVPElement(TwoStepUpdatedLagrangianVPElement  const& rOther)
  :Element(rOther)
  {

    // if (mOldFgrad.size() !=   rOther.mOldFgrad.size() )
    //   mOldFgrad.resize(  rOther.mOldFgrad.size());

    // for(unsigned int i=0; i<  rOther.mOldFgrad.size(); i++)
    //   {
    //     mOldFgrad[i] =   rOther.mOldFgrad[i];
    //   }


    if ( mCurrentTotalCauchyStress.size() !=   rOther.mCurrentTotalCauchyStress.size() )
      mCurrentTotalCauchyStress.resize(  rOther.mCurrentTotalCauchyStress.size());

    for(unsigned int i=0; i<  rOther.mCurrentTotalCauchyStress.size(); i++)
      {
        mCurrentTotalCauchyStress[i] =   rOther.mCurrentTotalCauchyStress[i];
      }


    if ( mCurrentDeviatoricCauchyStress.size() !=   rOther.mCurrentDeviatoricCauchyStress.size() )
      mCurrentDeviatoricCauchyStress.resize(  rOther.mCurrentDeviatoricCauchyStress.size());

    for(unsigned int i=0; i<  rOther.mCurrentDeviatoricCauchyStress.size(); i++)
      {
        mCurrentDeviatoricCauchyStress[i] =   rOther.mCurrentDeviatoricCauchyStress[i];
      }


    if ( mUpdatedTotalCauchyStress.size() !=   rOther.mUpdatedTotalCauchyStress.size() )
      mUpdatedTotalCauchyStress.resize(  rOther.mUpdatedTotalCauchyStress.size());

    for(unsigned int i=0; i<  rOther.mUpdatedTotalCauchyStress.size(); i++)
      {
        mUpdatedTotalCauchyStress[i] =   rOther.mUpdatedTotalCauchyStress[i];
      }


    if ( mUpdatedDeviatoricCauchyStress.size() !=   rOther.mUpdatedDeviatoricCauchyStress.size() )
      mUpdatedDeviatoricCauchyStress.resize(  rOther.mUpdatedDeviatoricCauchyStress.size());

    for(unsigned int i=0; i<  rOther.mUpdatedDeviatoricCauchyStress.size(); i++)
      {
        mUpdatedDeviatoricCauchyStress[i] =   rOther.mUpdatedDeviatoricCauchyStress[i];
      }
    
  }


  template< unsigned int TDim >
  Element::Pointer TwoStepUpdatedLagrangianVPElement<TDim>::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {

    TwoStepUpdatedLagrangianVPElement NewElement(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    // if ( NewElement.mOldFgrad.size() !=  mOldFgrad.size() )
    //   NewElement.mOldFgrad.resize( mOldFgrad.size());

    // for(unsigned int i=0; i< mOldFgrad.size(); i++)
    //   {
    //     NewElement.mOldFgrad[i] =  mOldFgrad[i];
    //   }


    if ( NewElement.mCurrentTotalCauchyStress.size() !=  mCurrentTotalCauchyStress.size() )
      NewElement.mCurrentTotalCauchyStress.resize( mCurrentTotalCauchyStress.size());

    for(unsigned int i=0; i< mCurrentTotalCauchyStress.size(); i++)
      {
        NewElement.mCurrentTotalCauchyStress[i] =  mCurrentTotalCauchyStress[i];
      }


    if ( NewElement.mCurrentDeviatoricCauchyStress.size() !=  mCurrentDeviatoricCauchyStress.size() )
      NewElement.mCurrentDeviatoricCauchyStress.resize( mCurrentDeviatoricCauchyStress.size());

    for(unsigned int i=0; i< mCurrentDeviatoricCauchyStress.size(); i++)
      {
        NewElement.mCurrentDeviatoricCauchyStress[i] =  mCurrentDeviatoricCauchyStress[i];
      }


    if ( NewElement.mUpdatedTotalCauchyStress.size() !=  mUpdatedTotalCauchyStress.size() )
      NewElement.mUpdatedTotalCauchyStress.resize( mUpdatedTotalCauchyStress.size());

    for(unsigned int i=0; i< mUpdatedTotalCauchyStress.size(); i++)
      {
        NewElement.mUpdatedTotalCauchyStress[i] =  mUpdatedTotalCauchyStress[i];
      }


    if ( NewElement.mUpdatedDeviatoricCauchyStress.size() !=  mUpdatedDeviatoricCauchyStress.size() )
      NewElement.mUpdatedDeviatoricCauchyStress.resize( mUpdatedDeviatoricCauchyStress.size());

    for(unsigned int i=0; i< mUpdatedDeviatoricCauchyStress.size(); i++)
      {
        NewElement.mUpdatedDeviatoricCauchyStress[i] =  mUpdatedDeviatoricCauchyStress[i];
      }

    
    return Element::Pointer( new TwoStepUpdatedLagrangianVPElement(NewElement) );
  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
								     VectorType& rRightHandSideVector,
								     ProcessInfo& rCurrentProcessInfo)
  { 
    KRATOS_TRY;

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
      {
      case 1:
	{
	  this->CalculateLocalMomentumEquations(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
	  break;
	}
      case 5:
      	{
      	  this->CalculateLocalContinuityEqForPressure(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
      	  break;
      	}

      default:
	{
	  KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for TWO_STEP_UPDATED_LAGRANGIAN_V_P_ELEMENT index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
	}
      }

    KRATOS_CATCH("");
  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::EquationIdVector(EquationIdVectorType& rResult,
								 ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY;

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
      {
      case 1:
	{
	  this->VelocityEquationIdVector(rResult,rCurrentProcessInfo);
	  break;
	}
      case 5:
	{
	  this->PressureEquationIdVector(rResult,rCurrentProcessInfo);
	  break;
	}
      case 6:
	{
	  this->VelocityEquationIdVector(rResult,rCurrentProcessInfo);
	  break;
	}
      default:
	{
	  KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
	}
      }

    KRATOS_CATCH("");
  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::GetDofList(DofsVectorType& rElementalDofList,
							   ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY;

    switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
      {
      case 1:
	{
	  this->GetVelocityDofList(rElementalDofList,rCurrentProcessInfo);
	  break;
	}
      case 5:
	{
	  this->GetPressureDofList(rElementalDofList,rCurrentProcessInfo);
	  break;
	}
      case 6:
	{
	  this->GetVelocityDofList(rElementalDofList,rCurrentProcessInfo);
	  break;
	}
      default:
	{
	  KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
	}
      }

    KRATOS_CATCH("");
  }

  template< unsigned int TDim >
  GeometryData::IntegrationMethod TwoStepUpdatedLagrangianVPElement<TDim>::GetIntegrationMethod() const
  {
    return GeometryData::GI_GAUSS_1;
  }



  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::CalculateLocalMomentumEquations(MatrixType& rLeftHandSideMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY; 

    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = TDim * NumNodes;
    // const SizeType NumNodes = rGeom.PointsNumber();
    // const SizeType LocalSize = TDim * NumNodes;

    MatrixType MassMatrix;
    MatrixType BulkMatrix;

    MassMatrix = ZeroMatrix(LocalSize,LocalSize);
    BulkMatrix = ZeroMatrix(LocalSize,LocalSize);

    // Check sizes and initialize
    if( rLeftHandSideMatrix.size1() != LocalSize )
      rLeftHandSideMatrix.resize(LocalSize,LocalSize);

    rLeftHandSideMatrix = ZeroMatrix(LocalSize,LocalSize);

    if( rRightHandSideVector.size() != LocalSize )
      rRightHandSideVector.resize(LocalSize);

    rRightHandSideVector = ZeroVector(LocalSize);

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

    
      const double TimeStep=0.5/rCurrentProcessInfo[BDF_COEFFICIENTS][2];

      ElementalVariables rElementalVariables;
      this->InitializeElementalVariables(rElementalVariables);

      // Loop on integration points
      for (unsigned int g = 0; g < NumGauss; g++)
	{
	  // const double GaussWeight = fabs(GaussWeights[g]);
	  const double GaussWeight = GaussWeights[g];
	  const ShapeFunctionsType& N = row(NContainer,g);
	  const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];



	  double Pressure=0;
	  double OldPressure=0;


	  this->EvaluateInPoint(Pressure,PRESSURE,N,0);

	  this->EvaluateInPoint(OldPressure,PRESSURE,N,1);

	  if(rCurrentProcessInfo[STEP]==1){
	    rElementalVariables.MeanPressure=Pressure;
	    // std::cout<<"Pressure "<< Pressure<<std::endl;
	  }else{
	    rElementalVariables.MeanPressure=OldPressure*0.5+Pressure*0.5;
	  }
	  rElementalVariables.MeanPressure=OldPressure*0.5+Pressure*0.5;  
	  // rElementalVariables.MeanPressure=Pressure;
	  computeElement=this->CalcMechanicsUpdated(rElementalVariables,rCurrentProcessInfo,g,N);
	  if(computeElement==true){

	    // Evaluate required variables at the integration point
	    double Density=0.0;
	    array_1d<double,3> BodyForce(3,0.0);

	    this->EvaluateInPoint(BodyForce,BODY_FORCE,N);

	    double DeviatoricCoeff = 0;
	    double VolumetricCoeff = 0;

	    this->ComputeMaterialParameters(Density,DeviatoricCoeff,VolumetricCoeff,TimeStep,N);

	    // Add integration point contribution to the local mass matrix
	    // double massWeight=GaussWeight*Density*2.0/TimeStep;
	    double DynamicWeight=GaussWeight*Density;
	    double MeanValueMass=0;
	    this->ComputeLumpedMassMatrix(MassMatrix,DynamicWeight,MeanValueMass);
	      
	    // this->ComputeMomentumMassTerm(MassMatrix,N,dynamicWeight);
	    this->AddExternalForces(rRightHandSideVector,Density,BodyForce,N,GaussWeight);

	    this->AddInternalForces(rRightHandSideVector,rDN_DX,rElementalVariables,GaussWeight);


	    double MeanValueMaterial=0.0;
	    this->ComputeMeanValueMaterialTangentMatrix(rElementalVariables,MeanValueMaterial,rDN_DX,DeviatoricCoeff,VolumetricCoeff,GaussWeight,MeanValueMass,TimeStep);
	    // double MeanValueMaterial=0.0;
	    // this->ComputeMeanValueMaterialTangentMatrix(rElementalVariables,MeanValueMaterial,rDN_DX,DeviatoricCoeff,VolumetricCoeff,GaussWeight);
	    // if(MeanValueMass!=0 && MeanValueMaterial!=0){
	    //   VolumetricCoeff*=MeanValueMass*2/TimeStep/MeanValueMaterial;
	    // }else{
	    //   std::cout<<" DANGEROUS ELEMENT!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
	    //   std::cout<<" MeanValueMass="<<MeanValueMass;
	    //   std::cout<<"\t MeanValueMaterial= "<<MeanValueMaterial;
	    //   std::cout<<"\t GaussWeight= "<<GaussWeight<<std::endl;
	    //   std::cout<<"\t Density= "<<Density;
	    //   std::cout<<"\t DeviatoricCoeff= "<<DeviatoricCoeff;
	    //   std::cout<<"\t VolumetricCoeff= "<<VolumetricCoeff<<std::endl;
	    //   VolumetricCoeff*=TimeStep;
	    // }

	    double theta=0.5;
	    // theta=1.0;
	    // Add viscous term
	    this->AddCompleteTangentTerm(rElementalVariables,rLeftHandSideMatrix,rDN_DX,DeviatoricCoeff,VolumetricCoeff,theta,GaussWeight);
	  }else{
	    for (SizeType n = 0; n < NumNodes; ++n)
	      {
		rGeom[n].Set(INTERFACE); //I set interface in order to not compute it in the continuity equation and the next non-linear iterations
	      }
	  }
	}


    
      // Add residual of previous iteration to RHS
      VectorType VelocityValues = ZeroVector(LocalSize);
      VectorType UpdatedAccelerations = ZeroVector(LocalSize);
      VectorType LastAccValues = ZeroVector(LocalSize);

      this->GetVelocityValues(VelocityValues,0);
      UpdatedAccelerations = 2.0*VelocityValues/TimeStep;
      this->GetAccelerationValues(LastAccValues,0);
      this->GetVelocityValues(VelocityValues,1);
      UpdatedAccelerations += -2.0*VelocityValues/TimeStep - LastAccValues; 
      noalias( rRightHandSideVector ) += -prod(MassMatrix,UpdatedAccelerations);
      noalias( rLeftHandSideMatrix ) +=  MassMatrix*2/TimeStep;

    }

     KRATOS_CATCH( "" );
 
  }


  template<>
  void TwoStepUpdatedLagrangianVPElement<2>::AddCompleteTangentTerm(ElementalVariables & rElementalVariables,
								    MatrixType& rDampingMatrix,
								    const ShapeFunctionDerivativesType& rDN_DX,
								    const double secondLame,
								    const double bulkModulus,
								    const double theta,
								    const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    MatrixType invGradDef=rElementalVariables.InvFgrad;

    SizeType FirstRow=0;
    SizeType FirstCol=0;

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
            rDampingMatrix(FirstRow,FirstCol) += Weight * ( (FourThirds * secondLame + bulkModulus)*  lagDNXi * lagDNXj + lagDNYi * lagDNYj * secondLame ) *theta;
            rDampingMatrix(FirstRow,FirstCol+1) += Weight * ( (nTwoThirds* secondLame + bulkModulus) *  lagDNXi * lagDNYj + lagDNYi * lagDNXj * secondLame )*theta;

            // Second Row
            rDampingMatrix(FirstRow+1,FirstCol) += Weight * ( (nTwoThirds * secondLame + bulkModulus) * lagDNYi * lagDNXj + lagDNXi * lagDNYj *  secondLame )*theta;
            rDampingMatrix(FirstRow+1,FirstCol+1) += Weight * ( (FourThirds * secondLame + bulkModulus) * lagDNYi * lagDNYj + lagDNXi * lagDNXj * secondLame )*theta;

            // Update Counter
            FirstRow += 2;
	  }
        FirstRow = 0;
        FirstCol += 2;
      }
  }
 

  template<>
  void TwoStepUpdatedLagrangianVPElement<3>::AddCompleteTangentTerm(ElementalVariables & rElementalVariables,
								    MatrixType& rDampingMatrix,
								    const ShapeFunctionDerivativesType& rDN_DX,
								    const double secondLame,
								    const double bulkModulus,
								    const double theta,
								    const double Weight){

   const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    MatrixType invGradDef=rElementalVariables.InvFgrad;

    SizeType FirstRow=0;
    SizeType FirstCol=0;

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
            rDampingMatrix(FirstRow,FirstCol)     += Weight * ( (FourThirds * secondLame + bulkModulus)*  lagDNXi * lagDNXj + (lagDNYi * lagDNYj +lagDNZi * lagDNZj) * secondLame ) *theta;
            rDampingMatrix(FirstRow,FirstCol+1)   += Weight * ( (nTwoThirds* secondLame + bulkModulus) *  lagDNXi * lagDNYj + lagDNYi * lagDNXj * secondLame )*theta;
            rDampingMatrix(FirstRow,FirstCol+2)   += Weight * ( (nTwoThirds* secondLame + bulkModulus) *  lagDNXi * lagDNZj + lagDNZi * lagDNXj * secondLame )*theta;

            // Second Row
            rDampingMatrix(FirstRow+1,FirstCol)   += Weight * ( (nTwoThirds * secondLame + bulkModulus) * lagDNYi * lagDNXj + lagDNXi * lagDNYj *  secondLame )*theta;
            rDampingMatrix(FirstRow+1,FirstCol+1) += Weight * ( (FourThirds * secondLame + bulkModulus) * lagDNYi * lagDNYj + (lagDNXi * lagDNXj + lagDNZi * lagDNZj) * secondLame )*theta;
            rDampingMatrix(FirstRow+1,FirstCol+2) += Weight * ( (nTwoThirds * secondLame + bulkModulus) * lagDNYi * lagDNZj + lagDNZi * lagDNYj *  secondLame )*theta;

            // Third Row
            rDampingMatrix(FirstRow+2,FirstCol)   += Weight * ( (nTwoThirds * secondLame + bulkModulus) * lagDNZi * lagDNXj + lagDNXi * lagDNZj *  secondLame )*theta;
            rDampingMatrix(FirstRow+2,FirstCol+1) += Weight * ( (nTwoThirds* secondLame + bulkModulus)  * lagDNZi * lagDNYj + lagDNYi * lagDNZj *  secondLame )*theta;
            rDampingMatrix(FirstRow+2,FirstCol+2) += Weight * ( (FourThirds * secondLame + bulkModulus) * lagDNZi * lagDNZj + (lagDNXi * lagDNXj + lagDNYi * lagDNYj) * secondLame )*theta;
	   
	    // Update Counter
            FirstRow += 3;
	  }
        FirstRow = 0;
        FirstCol += 3;
      }
  }




  template< unsigned int TDim >
  int TwoStepUpdatedLagrangianVPElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo)
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
  void TwoStepUpdatedLagrangianVPElement<2>::VelocityEquationIdVector(EquationIdVectorType& rResult,
								      ProcessInfo& rCurrentProcessInfo)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = NumNodes*2;

    SizeType LocalIndex = 0;

    if (rResult.size() != LocalSize)
      rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);

    for (SizeType i = 0; i < NumNodes; ++i)
      {
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X,xpos).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
      }

  }

  template<>
  void TwoStepUpdatedLagrangianVPElement<3>::VelocityEquationIdVector(EquationIdVectorType& rResult,
								      ProcessInfo& rCurrentProcessInfo)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3*NumNodes;

    SizeType LocalIndex = 0;

    if (rResult.size() != LocalSize)
      rResult.resize(LocalSize, false);

    const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);

    for (SizeType i = 0; i < NumNodes; ++i)
      {
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X,xpos).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
        rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Z,xpos+2).EquationId();
      }
  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::PressureEquationIdVector(EquationIdVectorType& rResult,
									 ProcessInfo& rCurrentProcessInfo)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    if (rResult.size() != NumNodes)
      rResult.resize(NumNodes);

    const unsigned int pos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);

    for (SizeType i = 0; i < NumNodes; ++i)
      rResult[i] = rGeom[i].GetDof(PRESSURE,pos).EquationId();
  }

  template<>
  void TwoStepUpdatedLagrangianVPElement<2>::GetVelocityDofList(DofsVectorType& rElementalDofList,
								ProcessInfo& rCurrentProcessInfo)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 2*NumNodes;

    if (rElementalDofList.size() != LocalSize)
      rElementalDofList.resize(LocalSize);

    SizeType LocalIndex = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
      }
  }

  template<>
  void TwoStepUpdatedLagrangianVPElement<3>::GetVelocityDofList(DofsVectorType& rElementalDofList,
								ProcessInfo& rCurrentProcessInfo)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3*NumNodes;

    if (rElementalDofList.size() != LocalSize)
      rElementalDofList.resize(LocalSize);

    SizeType LocalIndex = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Z);
      }
  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::GetPressureDofList(DofsVectorType& rElementalDofList,
								   ProcessInfo& rCurrentProcessInfo)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    if (rElementalDofList.size() != NumNodes)
      rElementalDofList.resize(NumNodes);

    SizeType LocalIndex = 0;
    for (SizeType i = 0; i < NumNodes; ++i)
      {
        rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PRESSURE);
      }
	
  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::GetPressureValues(Vector& rValues,
								  const int Step)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    if (rValues.size() != NumNodes) rValues.resize(NumNodes);

    for (SizeType i = 0; i < NumNodes; ++i){
      rValues[i] = rGeom[i].FastGetSolutionStepValue(PRESSURE,Step);

     //  if(rValues[i]==0 && rGeom[i].Is(FREE_SURFACE))
     // 	std::cout<<"                     pressure = 0 for this free-surface node               "<<std::endl;
     // if(rValues[i]==0 && rGeom[i].Is(RIGID))
     // 	std::cout<<"                     pressure = 0 for this rigid node               "<<std::endl;


      if(rGeom[i].Is(BOUNDARY)){
	rGeom[i].FastGetSolutionStepValue(INTERF) = 1;

      }else{
      	rGeom[i].FastGetSolutionStepValue(INTERF) = 0;
      }

      if(rGeom[i].Is(FREE_SURFACE)){
	rGeom[i].FastGetSolutionStepValue(FREESURFACE) = 1;

      }else{
      	rGeom[i].FastGetSolutionStepValue(FREESURFACE) = 0;
      }

    }
  }





  template<  unsigned int TDim>
  void TwoStepUpdatedLagrangianVPElement<TDim>::CalcMeanPressure(double &meanPressure,
								 const int Step)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    meanPressure=0;
    double coeff=0.0;
    for (SizeType i = 0; i < NumNodes; ++i)
      {
	meanPressure+=rGeom[i].GetSolutionStepValue(PRESSURE,Step);
	coeff+=1.0;
      }
    meanPressure*=1.0/coeff;
  }



  template< >
  void TwoStepUpdatedLagrangianVPElement<2>::CalcMeanVelocity(double &meanVelocity,
							      const int Step)
  {

    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    //SizeType Index = 0;
   
    double velX=0;
    double velY=0;
    for (SizeType i = 0; i < NumNodes; ++i)
      {
        velX += rGeom[i].FastGetSolutionStepValue(VELOCITY_X,Step)/3.0;
        velY += rGeom[i].FastGetSolutionStepValue(VELOCITY_Y,Step)/3.0;
      }

    meanVelocity=velX*velX+velY*velY;
    meanVelocity=sqrt(meanVelocity);

  }



  template< >
  void TwoStepUpdatedLagrangianVPElement<3>::CalcMeanVelocity(double &meanVelocity,
							      const int Step)
  {

    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    //SizeType Index = 0;
   
    double velX=0;
    double velY=0;
    double velZ=0;
    for (SizeType i = 0; i < NumNodes; ++i)
      {
        velX += rGeom[i].FastGetSolutionStepValue(VELOCITY_X,Step)*0.25;
        velY += rGeom[i].FastGetSolutionStepValue(VELOCITY_Y,Step)*0.25;
        velZ += rGeom[i].FastGetSolutionStepValue(VELOCITY_Z,Step)*0.25;
      }

    meanVelocity=velX*velX+velY*velY+velZ*velZ;
    meanVelocity=sqrt(meanVelocity);

  }





template< unsigned int TDim>
void TwoStepUpdatedLagrangianVPElement<TDim>::CalculateDeltaPosition(Matrix & rDeltaPosition)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = this->GetGeometry().PointsNumber();
    unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

    rDeltaPosition = zero_matrix<double>( number_of_nodes , dimension);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        array_1d<double, 3 > & CurrentDisplacement  = this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        array_1d<double, 3 > & PreviousDisplacement = this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            rDeltaPosition(i,j) = CurrentDisplacement[j]-PreviousDisplacement[j];
        }
    }

    KRATOS_CATCH( "" )
}


  template<>
  void TwoStepUpdatedLagrangianVPElement<2>::GetVelocityValues(Vector& rValues,
							       const int Step)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 2*NumNodes;

    if (rValues.size() != LocalSize) rValues.resize(LocalSize);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_X,Step);
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y,Step);
      }
  }


  template<>
  void TwoStepUpdatedLagrangianVPElement<3>::GetVelocityValues(Vector& rValues,
							       const int Step)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3*NumNodes;

    if (rValues.size() != LocalSize) rValues.resize(LocalSize);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_X,Step);
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y,Step);
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Z,Step);
      }
  }



  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::GetElementalAcceleration(Vector& meanAcceleration,
									 const int Step,
									 const double TimeStep)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    double count=0;
    for (SizeType i = 0; i < NumNodes; ++i)
      {
	meanAcceleration+= 0.5/TimeStep*(rGeom[i].FastGetSolutionStepValue(VELOCITY,0)-rGeom[i].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[i].FastGetSolutionStepValue(ACCELERATION,1); 
	// meanAcceleration+=rGeom[i].FastGetSolutionStepValue(ACCELERATION,step);
	count+=1.0;
      }
    meanAcceleration*=1.0/count;
  }

  template< >
  void TwoStepUpdatedLagrangianVPElement<2>::GetAccelerationValues(Vector& rValues,
								   const int Step)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 2*NumNodes;

    if (rValues.size() != LocalSize) rValues.resize(LocalSize);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(ACCELERATION_X,Step);
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(ACCELERATION_Y,Step);
      }
  }

  template< >
  void TwoStepUpdatedLagrangianVPElement<3>::GetAccelerationValues(Vector& rValues,
								   const int Step)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = 3*NumNodes;

    if (rValues.size() != LocalSize) rValues.resize(LocalSize);

    SizeType Index = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(ACCELERATION_X,Step);
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(ACCELERATION_Y,Step);
        rValues[Index++] = rGeom[i].FastGetSolutionStepValue(ACCELERATION_Z,Step);
      }
  }





  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::CalculateGeometryData(ShapeFunctionDerivativesArrayType &rDN_DX,
								      Matrix &NContainer,
								      Vector &rGaussWeights)
  {
    const GeometryType& rGeom = this->GetGeometry();
    Vector DetJ;
    rGeom.ShapeFunctionsIntegrationPointsGradients(rDN_DX,DetJ,GeometryData::GI_GAUSS_1);
    NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_1);

    rGaussWeights.resize(rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_1),false);

    for (unsigned int g = 0; g < rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_1); g++){
      // rGaussWeights[g] = fabs(DetJ[g] * IntegrationPoints[g].Weight());
      rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
      if(rGaussWeights[g]<0)
	std::cout<<"NEGATIVE GAUSS WEIGHT "<<rGaussWeights[g]<<std::endl;
  
    }
    

  }

  template< unsigned int TDim >
  double TwoStepUpdatedLagrangianVPElement<TDim>::ElementSize(/*ShapeFunctionDerivativesType &rDN_DX*/)
  {
    const GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();


    double ElemSize =0;
    array_1d<double,3> Edge(3,0.0);
    Edge = rGeom[1].Coordinates() - rGeom[0].Coordinates();
    double Length = Edge[0]*Edge[0];
    for (SizeType d = 1; d < TDim; d++)
      Length += Edge[d]*Edge[d];
    ElemSize+=sqrt(Length);

    double count=1.0;
    for (SizeType i = 2; i < NumNodes; i++){
      for(SizeType j = 0; j < i; j++)
        {
    	  Edge = rGeom[i].Coordinates() - rGeom[j].Coordinates();
    	  Length = Edge[0]*Edge[0];
    	  for (SizeType d = 1; d < TDim; d++){
    	    Length += Edge[d]*Edge[d];
    	  }
    	  ElemSize+=sqrt(Length);
    	  count+=1.0;
        }
    }
    ElemSize*=1.0/count; 
    return ElemSize;

    // // calculate minimum element length (used in stabilization Tau)
    // array_1d<double,3> Edge(3,0.0);
    // Edge = rGeom[1].Coordinates() - rGeom[0].Coordinates();
    // double ElemSize = Edge[0]*Edge[0];
    // for (SizeType d = 1; d < TDim; d++)
    //   ElemSize += Edge[d]*Edge[d];

    // for (SizeType i = 2; i < NumNodes; i++)
    //   for(SizeType j = 0; j < i; j++)
    //     {
    // 	  Edge = rGeom[i].Coordinates() - rGeom[j].Coordinates();
    // 	  //for computing minimum distance
    // 	  double Length = Edge[0]*Edge[0];
    // 	  for (SizeType d = 1; d < TDim; d++){
    // 	    Length += Edge[d]*Edge[d];
    // 	  }
    // 	  if (Length < ElemSize) ElemSize = Length;
    //     }

    // return sqrt(ElemSize);

  }

template< unsigned int TDim>
bool TwoStepUpdatedLagrangianVPElement<TDim>::CalcStrainRate(ElementalVariables & rElementalVariables,
							     const ProcessInfo &rCurrentProcessInfo,
							     unsigned int g,
							     const double theta)
{

  ShapeFunctionDerivativesArrayType DN_DX;
  Matrix NContainer;
  VectorType GaussWeights;
  this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);

  // const GeometryType& rGeom = this->this->GetGeometry();
  const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];

  bool computeElement=true;

  this->CalcFGrad(rDN_DX,
		  rElementalVariables.Fgrad,
		  rElementalVariables.InvFgrad,
		  rElementalVariables.DetFgrad,
		  rCurrentProcessInfo,
		  theta);

  //it computes the material time derivative of the deformation gradient and its jacobian and inverse
  this->CalcVelDefGrad(rDN_DX,
		       rElementalVariables.FgradVel,
		       rElementalVariables.InvFgradVel,
		       rElementalVariables.DetFgradVel,
		       theta);

  //it computes the spatial velocity gradient tensor --> [L_ij]=dF_ik*invF_kj
  this->CalcSpatialVelocityGrad(rElementalVariables.InvFgrad,
				rElementalVariables.FgradVel,
				rElementalVariables.SpatialVelocityGrad);
  
  this->CalcVolDefRateFromSpatialVelGrad(rElementalVariables.VolumetricDefRate,
					 rElementalVariables.SpatialVelocityGrad);
  

  // this->CalcVolumetricDefRate(rDN_DX,
  // 			      rElementalVariables.VolumetricDefRate,
  // 			      rElementalVariables.InvFgrad,
  //                          theta );

  // //it checks whether tr(l) == div(v)
  // CheckStrain1(rElementalVariables.VolumetricDefRate,
  // 	       rElementalVariables.SpatialVelocityGrad);
      
  // CheckStrain2(rElementalVariables.SpatialVelocityGrad,
  // 	       rElementalVariables.Fgrad,
  // 	       rElementalVariables.FgradVel);
 
  //it computes Material time Derivative of Green Lagrange strain tensor in MATERIAL configuration --> [D(E)/Dt]
  this->CalcMDGreenLagrangeMaterial(rElementalVariables.Fgrad,
				    rElementalVariables.FgradVel,
				    rElementalVariables.MDGreenLagrangeMaterial);
      
  //it computes Material time Derivative of Green Lagrange strain tensor in SPATIAL configuration  --> [d]
  this->CalcSpatialDefRate(rElementalVariables.MDGreenLagrangeMaterial,
			   rElementalVariables.InvFgrad,
			   rElementalVariables.SpatialDefRate);

  computeElement=CheckStrain3(rElementalVariables.SpatialDefRate,
			      rElementalVariables.SpatialVelocityGrad);

  this->CalcDeviatoricInvariant(rElementalVariables.SpatialDefRate,
				rElementalVariables.DeviatoricInvariant);

  return computeElement;

}  




  

template <unsigned int TDim > 
void TwoStepUpdatedLagrangianVPElement<TDim>::CalcFGrad(const ShapeFunctionDerivativesType& rDN_DX,
							MatrixType &Fgrad,
							MatrixType &invFgrad,
							double &FJacobian,
							const ProcessInfo& rCurrentProcessInfo,
							const double theta)
{
  GeometryType& rGeom = this->GetGeometry();
  const SizeType NumNodes = rGeom.PointsNumber();
  const SizeType LocalSize = TDim*NumNodes;
  VectorType  NodePosition= ZeroVector(LocalSize);
  this->GetPositions(NodePosition,rCurrentProcessInfo,theta);

  Fgrad.resize(TDim,TDim);

  Fgrad=ZeroMatrix(TDim,TDim);
  for (SizeType i = 0; i < TDim; i++)
    {
      for (SizeType j = 0; j < TDim; j++)
	{
	  for (SizeType k = 0; k < NumNodes; k++)
	    {
	      Fgrad(i,j)+= NodePosition[TDim*k+i]*rDN_DX(k,j);
	    }
	}
    }


 //Inverse
  invFgrad.resize(TDim,TDim);
  invFgrad=ZeroMatrix(TDim,TDim);
  FJacobian=1;


  MathUtils<double>::InvertMatrix( Fgrad, invFgrad, FJacobian );

  // Fgrad.resize(2,2);

  // Fgrad(0,0)= NodePosition[0]*rDN_DX(0,0)+NodePosition[2]*rDN_DX(1,0)+NodePosition[4]*rDN_DX(2,0);
  // Fgrad(0,1)= NodePosition[0]*rDN_DX(0,1)+NodePosition[2]*rDN_DX(1,1)+NodePosition[4]*rDN_DX(2,1);
  // Fgrad(1,0)= NodePosition[1]*rDN_DX(0,0)+NodePosition[3]*rDN_DX(1,0)+NodePosition[5]*rDN_DX(2,0);
  // Fgrad(1,1)= NodePosition[1]*rDN_DX(0,1)+NodePosition[3]*rDN_DX(1,1)+NodePosition[5]*rDN_DX(2,1);

  // //Determinant of the material time derivative of the deformation gradient tensor
  // FJacobian= Fgrad(0,0)*Fgrad(1,1)-Fgrad(0,1)*Fgrad(1,0);
   
  // //Inverse
  // invFgrad.resize(2,2);
  // if(FJacobian==0){
  //   FJacobian=1;
  // }else{
  //   invFgrad(0,0)=  (Fgrad(1,1)/FJacobian);
  //   invFgrad(0,1)= -(Fgrad(0,1)/FJacobian);
  //   invFgrad(1,0)= -(Fgrad(1,0)/FJacobian);
  //   invFgrad(1,1)=  (Fgrad(0,0)/FJacobian); 
  // }

    // std::cout<<"             "<< FJacobian;
    // std::cout<<"__ "<< invFgrad(0,0);
    // std::cout<<"__ "<< invFgrad(1,1);
    // std::cout<<"__ "<< invFgrad(0,1);
    // std::cout<<"__ "<< invFgrad(1,0)<<std::endl;
    // std::cout<<":: "<< Fgrad(0,0);
    // std::cout<<":: "<< Fgrad(1,1);
    // std::cout<<":: "<< Fgrad(0,1);
    // std::cout<<":: "<< Fgrad(1,0)<<std::endl;


}

  

template < unsigned int TDim> 
void TwoStepUpdatedLagrangianVPElement<TDim>::CalcVelDefGrad(const ShapeFunctionDerivativesType& rDN_DX,
							     MatrixType &FgradVel,
							     MatrixType &invFgradVel,
							     double &FVelJacobian,
							     const double theta)
{
  GeometryType& rGeom = this->GetGeometry();
  const SizeType NumNodes = rGeom.PointsNumber();
  const SizeType LocalSize = TDim*NumNodes;
  VectorType CurrentVelocities = ZeroVector(LocalSize);
  VectorType UpdatedVelocities = ZeroVector(LocalSize);
  VectorType RHSVelocities = ZeroVector(LocalSize);
  this->GetVelocityValues(CurrentVelocities,1);
  this->GetVelocityValues(UpdatedVelocities,0); 
  RHSVelocities=CurrentVelocities*theta+UpdatedVelocities*(1.0-theta);

  FgradVel.resize(TDim,TDim);

  FgradVel=ZeroMatrix(TDim,TDim);
  for (SizeType i = 0; i < TDim; i++)
    {
      for (SizeType j = 0; j < TDim; j++)
	{
	  for (SizeType k = 0; k < NumNodes; k++)
	    {
	      FgradVel(i,j)+= RHSVelocities[TDim*k+i]*rDN_DX(k,j);
	    }
	}
    }


 //Inverse
  invFgradVel.resize(TDim,TDim);
  invFgradVel=ZeroMatrix(TDim,TDim);
  FVelJacobian=1;

  MathUtils<double>::InvertMatrix( FgradVel, invFgradVel, FVelJacobian );


}


template < > 
void TwoStepUpdatedLagrangianVPElement<2>::CalcVolumetricDefRate(const ShapeFunctionDerivativesType& rDN_DX,
								 double &volumetricDefRate,
								 MatrixType &invGradDef,
								 const double theta)
{


  GeometryType& rGeom = this->GetGeometry();
  const SizeType NumNodes = rGeom.PointsNumber();
  const SizeType LocalSize = 2*NumNodes;
  VectorType CurrentVelocities = ZeroVector(LocalSize);
  VectorType UpdatedVelocities = ZeroVector(LocalSize);
  VectorType RHSVelocities = ZeroVector(LocalSize);
  this->GetVelocityValues(CurrentVelocities,1);
  this->GetVelocityValues(UpdatedVelocities,0);
  RHSVelocities=CurrentVelocities*theta+UpdatedVelocities*(1.0-theta);

  double lagDNX0=rDN_DX(0,0)*invGradDef(0,0)+rDN_DX(0,1)*invGradDef(1,0);
  double lagDNX1=rDN_DX(1,0)*invGradDef(0,0)+rDN_DX(1,1)*invGradDef(1,0);
  double lagDNX2=rDN_DX(2,0)*invGradDef(0,0)+rDN_DX(2,1)*invGradDef(1,0);
  double lagDNY0=rDN_DX(0,0)*invGradDef(0,1)+rDN_DX(0,1)*invGradDef(1,1);
  double lagDNY1=rDN_DX(1,0)*invGradDef(0,1)+rDN_DX(1,1)*invGradDef(1,1);
  double lagDNY2=rDN_DX(2,0)*invGradDef(0,1)+rDN_DX(2,1)*invGradDef(1,1);


  volumetricDefRate= lagDNX0*RHSVelocities[0] + lagDNX1*RHSVelocities[2] + lagDNX2*RHSVelocities[4];
  volumetricDefRate+=lagDNY0*RHSVelocities[1] + lagDNY1*RHSVelocities[3] + lagDNY2*RHSVelocities[5];

}

template < unsigned int TDim > 
void TwoStepUpdatedLagrangianVPElement<TDim>::CalcSpatialVelocityGrad(MatrixType &invFgrad,
								      MatrixType &VelDefgrad,
								      MatrixType &SpatialVelocityGrad)
{
  SpatialVelocityGrad.resize(TDim,TDim);
  
  SpatialVelocityGrad=prod(VelDefgrad,invFgrad);

  // SpatialVelocityGrad(0,0)=VelDefgrad(0,0)*invFgrad(0,0) + VelDefgrad(0,1)*invFgrad(1,0);
  // SpatialVelocityGrad(0,1)=VelDefgrad(0,0)*invFgrad(0,1) + VelDefgrad(0,1)*invFgrad(1,1);
  // SpatialVelocityGrad(1,0)=VelDefgrad(1,0)*invFgrad(0,0) + VelDefgrad(1,1)*invFgrad(1,0);
  // SpatialVelocityGrad(1,1)=VelDefgrad(1,0)*invFgrad(0,1) + VelDefgrad(1,1)*invFgrad(1,1);

}


template < unsigned int TDim > 
void TwoStepUpdatedLagrangianVPElement<TDim>::CalcVolDefRateFromSpatialVelGrad(double &volumetricDefRate,
										    MatrixType &SpatialVelocityGrad)
{
  volumetricDefRate=0;
  for (SizeType i = 0; i < TDim; i++)
    {
      volumetricDefRate+=SpatialVelocityGrad(i,i);
    }
}



template < > 
void TwoStepUpdatedLagrangianVPElement<2>::CheckStrain1(double &VolumetricDefRate,
							     MatrixType &SpatialVelocityGrad)
{
  double trace_l=SpatialVelocityGrad(0,0)+SpatialVelocityGrad(1,1);
  if(fabs(trace_l-VolumetricDefRate)<0.0000001){
  }else{
    std::cout<<" ERROR IN CHECKSTRAIN(1) -> ";
    std::cout<<"trace_l= "<<trace_l<<" VolDefRate"<<VolumetricDefRate<<std::endl;
  }
}

template < > 
void TwoStepUpdatedLagrangianVPElement<2>::CalcMDGreenLagrangeMaterial(MatrixType &Fgrad,
								       MatrixType &VelDefgrad, 
								       VectorType &MDGreenLagrangeMaterial)
{

  // x-component
  MDGreenLagrangeMaterial[0]=VelDefgrad(0,0)*Fgrad(0,0) + VelDefgrad(1,0)*Fgrad(1,0);
  // y-component
  MDGreenLagrangeMaterial[1]=VelDefgrad(1,1)*Fgrad(1,1) + VelDefgrad(0,1)*Fgrad(0,1);
  // xy-component
  MDGreenLagrangeMaterial[2]=(VelDefgrad(0,0)*Fgrad(0,1) + VelDefgrad(1,0)*Fgrad(1,1) +
  			      VelDefgrad(0,1)*Fgrad(0,0) + VelDefgrad(1,1)*Fgrad(1,0))*0.5;
}



template < > 
void TwoStepUpdatedLagrangianVPElement<3>::CalcMDGreenLagrangeMaterial(MatrixType &Fgrad,
									    MatrixType &VelDefgrad, 
									    VectorType &MDGreenLagrangeMaterial)
{
  MatrixType FgradTransp(3,3);
  MatrixType VelDefgradTransp(3,3);
  MatrixType part1(3,3);
  MatrixType part2(3,3);

  FgradTransp=Fgrad;
  FgradTransp(0,1)=Fgrad(1,0);
  FgradTransp(0,2)=Fgrad(2,0);
  FgradTransp(1,0)=Fgrad(0,1);
  FgradTransp(1,2)=Fgrad(2,1);
  FgradTransp(2,0)=Fgrad(0,2);
  FgradTransp(2,1)=Fgrad(1,2);

  VelDefgradTransp=VelDefgrad;
  VelDefgradTransp(0,1)=VelDefgrad(1,0);
  VelDefgradTransp(0,2)=VelDefgrad(2,0);
  VelDefgradTransp(1,0)=VelDefgrad(0,1);
  VelDefgradTransp(1,2)=VelDefgrad(2,1);
  VelDefgradTransp(2,0)=VelDefgrad(0,2);
  VelDefgradTransp(2,1)=VelDefgrad(1,2);

  part1=prod(VelDefgradTransp,Fgrad);
  part2=prod(FgradTransp,VelDefgrad);

  MDGreenLagrangeMaterial[0]= ( part1(0,0) + part2(0,0) ) * 0.5;  //xx-component
  MDGreenLagrangeMaterial[1]= ( part1(1,1) + part2(1,1) ) * 0.5;  //yy-component
  MDGreenLagrangeMaterial[2]= ( part1(2,2) + part2(2,2) ) * 0.5;  //zz-component
  MDGreenLagrangeMaterial[3]= ( part1(0,1) + part2(0,1) ) * 0.5;  //xy-component
  MDGreenLagrangeMaterial[4]= ( part1(0,2) + part2(0,2) ) * 0.5;  //xz-component
  MDGreenLagrangeMaterial[5]= ( part1(1,2) + part2(1,2) ) * 0.5;  //yz-component

}




template < > 
void TwoStepUpdatedLagrangianVPElement<2>::CalcSpatialDefRate(VectorType &MDGreenLagrangeMaterial,
							      MatrixType &invFgrad,
							      VectorType &SpatialDefRate)
{
  // x-component
  SpatialDefRate[0]= invFgrad(0,0)*MDGreenLagrangeMaterial[0]*invFgrad(0,0) + 
    invFgrad(1,0)*MDGreenLagrangeMaterial[2]*invFgrad(0,0)*2 +
    invFgrad(1,0)*MDGreenLagrangeMaterial[1]*invFgrad(1,0);
  // y-component
  SpatialDefRate[1]= invFgrad(0,1)*MDGreenLagrangeMaterial[0]*invFgrad(0,1) + 
    invFgrad(0,1)*MDGreenLagrangeMaterial[2]*invFgrad(1,1)*2 +
    invFgrad(1,1)*MDGreenLagrangeMaterial[1]*invFgrad(1,1);
  // xy-component
  SpatialDefRate[2]=invFgrad(0,0)*MDGreenLagrangeMaterial[0]*invFgrad(0,1) + 
    invFgrad(0,0)*MDGreenLagrangeMaterial[2]*invFgrad(1,1) +
    invFgrad(1,0)*MDGreenLagrangeMaterial[2]*invFgrad(0,1) +
    invFgrad(1,0)*MDGreenLagrangeMaterial[1]*invFgrad(1,1);
}


template < > 
void TwoStepUpdatedLagrangianVPElement<3>::CalcSpatialDefRate(VectorType &MDGreenLagrangeMaterial,
							      MatrixType &invFgrad,
							      VectorType &SpatialDefRate)
{
  MatrixType MDGLM(3,3);
  MatrixType invFgradTransp(3,3);
  MatrixType part1(3,3);
  MatrixType totalMatrix(3,3);
  MDGLM=ZeroMatrix(3,3);
  invFgradTransp=ZeroMatrix(3,3);
  part1=ZeroMatrix(3,3);
  totalMatrix=ZeroMatrix(3,3);


  invFgradTransp=invFgrad;
  invFgradTransp(0,1)=invFgrad(1,0);
  invFgradTransp(0,2)=invFgrad(2,0);
  invFgradTransp(1,0)=invFgrad(0,1);
  invFgradTransp(1,2)=invFgrad(2,1);
  invFgradTransp(2,0)=invFgrad(0,2);
  invFgradTransp(2,1)=invFgrad(1,2);

  MDGLM(0,0)=MDGreenLagrangeMaterial[0];  //XX-component;
  MDGLM(1,1)=MDGreenLagrangeMaterial[1];  //YY-component;
  MDGLM(2,2)=MDGreenLagrangeMaterial[2];  //ZZ-component;
  MDGLM(0,1)=MDGreenLagrangeMaterial[3];  //XY-component;
  MDGLM(1,0)=MDGreenLagrangeMaterial[3];  //XY-component;
  MDGLM(0,2)=MDGreenLagrangeMaterial[4];  //ZX-component;
  MDGLM(2,0)=MDGreenLagrangeMaterial[4];  //ZX-component;
  MDGLM(1,2)=MDGreenLagrangeMaterial[5];  //YZ-component;
  MDGLM(2,1)=MDGreenLagrangeMaterial[5];  //YZ-component;

  part1=prod(MDGLM,invFgrad);

  totalMatrix=prod(invFgradTransp,part1);
 
  SpatialDefRate[0]=totalMatrix(0,0);
  SpatialDefRate[1]=totalMatrix(1,1);
  SpatialDefRate[2]=totalMatrix(2,2);
  SpatialDefRate[3]=totalMatrix(0,1);
  SpatialDefRate[4]=totalMatrix(0,2);
  SpatialDefRate[5]=totalMatrix(1,2);
}



template < > 
void TwoStepUpdatedLagrangianVPElement<2>::CalcDeviatoricInvariant(VectorType &SpatialDefRate,
									double &DeviatoricInvariant)
{
  double trace_d=SpatialDefRate[0]+SpatialDefRate[1];
  double dev_X=SpatialDefRate[0]-trace_d/3.0;
  double dev_Y=SpatialDefRate[1]-trace_d/3.0;
  DeviatoricInvariant=sqrt(2*(dev_X*dev_X+SpatialDefRate[2]*SpatialDefRate[2]+ dev_Y*dev_Y));

}


template < > 
void TwoStepUpdatedLagrangianVPElement<3>::CalcDeviatoricInvariant(VectorType &SpatialDefRate,
								   double &DeviatoricInvariant)
{
  double trace_d=SpatialDefRate[0]+SpatialDefRate[1]+SpatialDefRate[2];
  double dev_X=SpatialDefRate[0]-trace_d/3.0;
  double dev_Y=SpatialDefRate[1]-trace_d/3.0;
  double dev_Z=SpatialDefRate[2]-trace_d/3.0;
  DeviatoricInvariant=sqrt(2*(dev_X*dev_X+dev_Y*dev_Y+dev_Z*dev_Z+
			      SpatialDefRate[3]*SpatialDefRate[3]+
			      SpatialDefRate[4]*SpatialDefRate[4]+
			      SpatialDefRate[5]*SpatialDefRate[5]));
}


template < > 
void TwoStepUpdatedLagrangianVPElement<2>::CalcNormalProjectionsForBoundRHSVector(VectorType &SpatialDefRate,
										  double& NormalAcceleration,
										  double& NormalProjSpatialDefRate,
										  const double TimeStep)
{
  array_1d<double, 3>  NormalA(3,0.0);
  NormalA.clear();
  array_1d<double, 3>  NormalB(3,0.0);
  NormalB.clear();
  array_1d<double, 3>  NormalMean(3,0.0);
  NormalMean.clear();
  array_1d<double, 3>  AccA(3,0.0);
  AccA.clear();
  array_1d<double, 3>  AccB(3,0.0);
  AccB.clear();
  array_1d<double, 3>  AccMean(3,0.0);
  AccMean.clear();
  // VectorType ElementalAcceleration = ZeroVector(2);
  // this->GetElementalAcceleration(ElementalAcceleration,0,TimeStep);
  GeometryType& rGeom = this->GetGeometry();

  if(rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE)){
    AccA= 0.5/TimeStep*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-
rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
    AccB= 0.5/TimeStep*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-
rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
    AccMean = AccA*0.5+ AccB*0.5;
    NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
    NormalB    = rGeom[1].FastGetSolutionStepValue(NORMAL);
    NormalMean = NormalA*0.5+ NormalB*0.5;
  }else if(rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE)){
    AccA= 0.5/TimeStep*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-
rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
    AccB= 0.5/TimeStep*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-
rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
    AccMean = AccA*0.5+ AccB*0.5;
    NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
    NormalB    = rGeom[2].FastGetSolutionStepValue(NORMAL);
    NormalMean = NormalA*0.5+ NormalB*0.5;
  }else if(rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE)){
    AccA= 0.5/TimeStep*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-
rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
    AccB= 0.5/TimeStep*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-
rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
    AccMean = AccA*0.5+ AccB*0.5;
    NormalA    = rGeom[1].FastGetSolutionStepValue(NORMAL);
    NormalB    = rGeom[2].FastGetSolutionStepValue(NORMAL);
    NormalMean = NormalA*0.5+ NormalB*0.5;
  }else  if((rGeom[0].Is(FREE_SURFACE) || rGeom[0].Is(RIGID) ) && (rGeom[1].Is(FREE_SURFACE)  || rGeom[1].Is(RIGID)) && !(rGeom[0].Is(RIGID) && rGeom[1].Is(RIGID)) ){
    AccA= 0.5/TimeStep*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-
rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
    AccB= 0.5/TimeStep*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-
rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
    AccMean = AccA*0.5+ AccB*0.5;
    NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
    NormalB    = rGeom[1].FastGetSolutionStepValue(NORMAL);
    NormalMean = NormalA*0.5+ NormalB*0.5;
  }else if((rGeom[2].Is(FREE_SURFACE) || rGeom[2].Is(RIGID) ) && (rGeom[0].Is(FREE_SURFACE)  || rGeom[0].Is(RIGID)) && !(rGeom[2].Is(RIGID) && rGeom[0].Is(RIGID))){
    AccA= 0.5/TimeStep*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-
rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
    AccB= 0.5/TimeStep*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-
rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
    AccMean = AccA*0.5+ AccB*0.5;
    NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
    NormalB    = rGeom[2].FastGetSolutionStepValue(NORMAL);
    NormalMean = NormalA*0.5+ NormalB*0.5;
  }else if((rGeom[1].Is(FREE_SURFACE) || rGeom[1].Is(RIGID) ) && (rGeom[2].Is(FREE_SURFACE)  || rGeom[2].Is(RIGID)) && !(rGeom[1].Is(RIGID) && rGeom[2].Is(RIGID))){
    AccA= 0.5/TimeStep*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-
rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
    AccB= 0.5/TimeStep*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-
rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
    AccMean = AccA*0.5+ AccB*0.5;
    NormalA    = rGeom[1].FastGetSolutionStepValue(NORMAL);
    NormalB    = rGeom[2].FastGetSolutionStepValue(NORMAL);
    NormalMean = NormalA*0.5+ NormalB*0.5;
  }

  NormalAcceleration=0.5*(NormalA[0]*AccA[0]+NormalA[1]*AccA[1]) +0.5*(NormalB[0]*AccB[0]+NormalB[1]*AccB[1]);
  // NormalAcceleration=NormalMean[0]*AccMean[0]+NormalMean[1]*AccMean[1];

 NormalProjSpatialDefRate=NormalMean[0]*SpatialDefRate[0]*NormalMean[0]+
    NormalMean[1]*SpatialDefRate[1]*NormalMean[1]+
    2*NormalMean[0]*SpatialDefRate[2]*NormalMean[1];

  // if(rGeom[0].Is(FREE_SURFACE) && rGeom[1].Is(FREE_SURFACE)){
  //   NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
  //   NormalA    += rGeom[1].FastGetSolutionStepValue(NORMAL);
  //   NormalMean = NormalA*0.5;
  // }else if(rGeom[0].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE)){
  //   NormalA    = rGeom[2].FastGetSolutionStepValue(NORMAL);
  //   NormalA    += rGeom[0].FastGetSolutionStepValue(NORMAL);
  //   NormalMean = NormalA*0.5;
  // }else if(rGeom[1].Is(FREE_SURFACE) && rGeom[2].Is(FREE_SURFACE)){
  //   NormalA    = rGeom[1].FastGetSolutionStepValue(NORMAL);
  //   NormalA    += rGeom[2].FastGetSolutionStepValue(NORMAL);
  //   NormalMean = NormalA*0.5;
  // }else  if((rGeom[0].Is(FREE_SURFACE) || rGeom[0].Is(RIGID) ) && (rGeom[1].Is(FREE_SURFACE)  || rGeom[1].Is(RIGID)) && !(rGeom[0].Is(RIGID) && rGeom[1].Is(RIGID)) ){
  //   NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
  //   NormalA    += rGeom[1].FastGetSolutionStepValue(NORMAL);
  //   NormalMean = NormalA*0.5;
  // }else if((rGeom[2].Is(FREE_SURFACE) || rGeom[2].Is(RIGID) ) && (rGeom[0].Is(FREE_SURFACE)  || rGeom[0].Is(RIGID)) && !(rGeom[2].Is(RIGID) && rGeom[0].Is(RIGID))){
  //   NormalA    = rGeom[2].FastGetSolutionStepValue(NORMAL);
  //   NormalA    += rGeom[0].FastGetSolutionStepValue(NORMAL);
  //   NormalMean = NormalA*0.5;
  // }else if((rGeom[1].Is(FREE_SURFACE) || rGeom[1].Is(RIGID) ) && (rGeom[2].Is(FREE_SURFACE)  || rGeom[2].Is(RIGID)) && !(rGeom[1].Is(RIGID) && rGeom[2].Is(RIGID))){
  //   NormalA    = rGeom[1].FastGetSolutionStepValue(NORMAL);
  //   NormalA    += rGeom[2].FastGetSolutionStepValue(NORMAL);
  //   NormalMean = NormalA*0.5;
  // }

  // NormalAcceleration=NormalMean[0]*MeanAcceleration[0]+NormalMean[1]*MeanAcceleration[1];
 
  // NormalProjSpatialDefRate=NormalMean[0]*SpatialDefRate[0]*NormalMean[0]+
  //   NormalMean[1]*SpatialDefRate[1]*NormalMean[1]+
  //   2*NormalMean[0]*SpatialDefRate[2]*NormalMean[1];

}


template < > 
void TwoStepUpdatedLagrangianVPElement<3>::CalcNormalProjectionsForBoundRHSVector(VectorType &SpatialDefRate,
										  double& NormalAcceleration,
										  double& NormalProjSpatialDefRate,
										  const double TimeStep)
{
  array_1d<double, 3>  NormalA(3,0.0);
  NormalA.clear();
  array_1d<double, 3>  NormalB(3,0.0);
  NormalB.clear();
  array_1d<double, 3>  NormalC(3,0.0);
  NormalC.clear();
  array_1d<double, 3>  NormalMean(3,0.0);
  NormalMean.clear();
  array_1d<double, 3>  AccA(3,0.0);
  AccA.clear();
  array_1d<double, 3>  AccB(3,0.0);
  AccB.clear();
  array_1d<double, 3>  AccC(3,0.0);
  AccC.clear();
  array_1d<double, 3>  AccMean(3,0.0);
  AccMean.clear();
  // VectorType ElementalAcceleration = ZeroVector(3);
  // this->GetElementalAcceleration(ElementalAcceleration,0,TimeStep);
  GeometryType& rGeom = this->GetGeometry();

  const SizeType NumNodes = this->GetGeometry().PointsNumber();
  unsigned int numBoundary=0;
  unsigned int numFreeSurf=0;
  for (SizeType i = 0; i < NumNodes; ++i)
    {
      if(rGeom[i].Is(FREE_SURFACE) || rGeom[i].Is(RIGID)){
	numBoundary++;
	if(rGeom[i].Is(FREE_SURFACE)){
	  numFreeSurf++;
	}
      }
    }

  if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)){
    AccA= 0.5/TimeStep*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
    AccB= 0.5/TimeStep*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
    AccC= 0.5/TimeStep*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
    AccMean = (AccA+AccB+AccC)*0.33333333;
    NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
    NormalB    = rGeom[1].FastGetSolutionStepValue(NORMAL);
    NormalC    = rGeom[2].FastGetSolutionStepValue(NORMAL);
    NormalMean = (NormalA + NormalB + NormalC)*0.33333333;
    NormalAcceleration=0.33333333*(NormalA[0]*AccA[0] + NormalA[1]*AccA[1] + NormalA[2]*AccA[2] + 
				   NormalB[0]*AccB[0] + NormalB[1]*AccB[1] + NormalB[2]*AccB[2] + 
				   NormalC[0]*AccC[0] + NormalC[1]*AccC[1] + NormalC[2]*AccC[2]);
  }else if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
    AccA= 0.5/TimeStep*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
    AccB= 0.5/TimeStep*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
    AccC= 0.5/TimeStep*(rGeom[3].FastGetSolutionStepValue(VELOCITY,0)-rGeom[3].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION,1); 
    AccMean = (AccA+AccB+AccC)*0.33333333;
    NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
    NormalB    = rGeom[1].FastGetSolutionStepValue(NORMAL);
    NormalC    = rGeom[3].FastGetSolutionStepValue(NORMAL);
    NormalMean = (NormalA + NormalB + NormalC)*0.33333333;
    NormalAcceleration=0.33333333*(NormalA[0]*AccA[0] + NormalA[1]*AccA[1] + NormalA[2]*AccA[2] + 
				   NormalB[0]*AccB[0] + NormalB[1]*AccB[1] + NormalB[2]*AccB[2] + 
				   NormalC[0]*AccC[0] + NormalC[1]*AccC[1] + NormalC[2]*AccC[2]);
  }else if(rGeom[0].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
    AccA= 0.5/TimeStep*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
    AccB= 0.5/TimeStep*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
    AccC= 0.5/TimeStep*(rGeom[3].FastGetSolutionStepValue(VELOCITY,0)-rGeom[3].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION,1); 
    AccMean = (AccA+AccB+AccC)*0.33333333;
    NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
    NormalB    = rGeom[2].FastGetSolutionStepValue(NORMAL);
    NormalC    = rGeom[3].FastGetSolutionStepValue(NORMAL);
    NormalMean = (NormalA + NormalB + NormalC)*0.33333333;
    NormalAcceleration=0.33333333*(NormalA[0]*AccA[0] + NormalA[1]*AccA[1] + NormalA[2]*AccA[2] + 
				   NormalB[0]*AccB[0] + NormalB[1]*AccB[1] + NormalB[2]*AccB[2] + 
				   NormalC[0]*AccC[0] + NormalC[1]*AccC[1] + NormalC[2]*AccC[2]);
  }else if(rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){      
    AccA= 0.5/TimeStep*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
    AccB= 0.5/TimeStep*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
    AccC= 0.5/TimeStep*(rGeom[3].FastGetSolutionStepValue(VELOCITY,0)-rGeom[3].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION,1); 
    AccMean = (AccA+AccB+AccC)*0.33333333;
    NormalA    = rGeom[1].FastGetSolutionStepValue(NORMAL);
    NormalB    = rGeom[2].FastGetSolutionStepValue(NORMAL);
    NormalC    = rGeom[3].FastGetSolutionStepValue(NORMAL);
    NormalMean = (NormalA + NormalB + NormalC)*0.33333333;
    NormalAcceleration=0.33333333*(NormalA[0]*AccA[0] + NormalA[1]*AccA[1] + NormalA[2]*AccA[2] + 
				   NormalB[0]*AccB[0] + NormalB[1]*AccB[1] + NormalB[2]*AccB[2] + 
				   NormalC[0]*AccC[0] + NormalC[1]*AccC[1] + NormalC[2]*AccC[2]);
  }else if(numFreeSurf>0){
      double count=0;
      if(rGeom[0].Is(FREE_SURFACE) || rGeom[0].Is(RIGID)){
	count+=1.0;
	NormalA   = rGeom[0].FastGetSolutionStepValue(NORMAL);
	AccA= 0.5/TimeStep*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
	NormalAcceleration+=NormalA[0]*AccA[0] + NormalA[1]*AccA[1] + NormalA[2]*AccA[2];
      }
      if(rGeom[1].Is(FREE_SURFACE) || rGeom[1].Is(RIGID)){
	count+=1.0;
	NormalA   = rGeom[1].FastGetSolutionStepValue(NORMAL);
	AccA= 0.5/TimeStep*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
	NormalAcceleration+=NormalA[0]*AccA[0] + NormalA[1]*AccA[1] + NormalA[2]*AccA[2];
      }
      if(rGeom[2].Is(FREE_SURFACE) || rGeom[2].Is(RIGID)){
	count+=1.0;
	NormalA   = rGeom[2].FastGetSolutionStepValue(NORMAL);
	AccA= 0.5/TimeStep*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
	NormalAcceleration+=NormalA[0]*AccA[0] + NormalA[1]*AccA[1] + NormalA[2]*AccA[2];
      }
      if(rGeom[3].Is(FREE_SURFACE) || rGeom[3].Is(RIGID)){
	count+=1.0;
	NormalA   = rGeom[3].FastGetSolutionStepValue(NORMAL);
	AccA= 0.5/TimeStep*(rGeom[3].FastGetSolutionStepValue(VELOCITY,0)-rGeom[3].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION,1); 
	NormalAcceleration+=NormalA[0]*AccA[0] + NormalA[1]*AccA[1] + NormalA[2]*AccA[2];
      }
      NormalAcceleration*=1.0/count;
  }

   
    NormalProjSpatialDefRate=NormalMean[0]*SpatialDefRate[0]*NormalMean[0]+
      NormalMean[1]*SpatialDefRate[1]*NormalMean[1]+
      NormalMean[2]*SpatialDefRate[2]*NormalMean[2]+
      2*NormalMean[0]*SpatialDefRate[3]*NormalMean[1]+
      2*NormalMean[0]*SpatialDefRate[4]*NormalMean[2]+
      2*NormalMean[1]*SpatialDefRate[5]*NormalMean[2];


  // if(numBoundary>2){
  //   if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)){
  //     NormalA     = rGeom[0].FastGetSolutionStepValue(NORMAL);
  //     NormalA    += rGeom[1].FastGetSolutionStepValue(NORMAL);
  //     NormalA    += rGeom[2].FastGetSolutionStepValue(NORMAL);
  //     NormalMean = NormalA*0.333333333333;
  //   }else if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
  //     NormalA     = rGeom[0].FastGetSolutionStepValue(NORMAL);
  //     NormalA    += rGeom[1].FastGetSolutionStepValue(NORMAL);
  //     NormalA    += rGeom[3].FastGetSolutionStepValue(NORMAL);
  //     NormalMean = NormalA*0.333333333333;
  //   }else if(rGeom[0].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
  //     NormalA     = rGeom[0].FastGetSolutionStepValue(NORMAL);
  //     NormalA    += rGeom[2].FastGetSolutionStepValue(NORMAL);
  //     NormalA    += rGeom[3].FastGetSolutionStepValue(NORMAL);
  //     NormalMean = NormalA*0.333333333333;
  //   }else if(rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){      
  //     NormalA     = rGeom[1].FastGetSolutionStepValue(NORMAL);
  //     NormalA    += rGeom[2].FastGetSolutionStepValue(NORMAL);
  //     NormalA    += rGeom[3].FastGetSolutionStepValue(NORMAL);
  //     NormalMean = NormalA*0.333333333333;
  //   }else  if(numFreeSurf>0){
  //     double count=0;
  //     if(rGeom[0].Is(FREE_SURFACE)){
  // 	count+=1.0;
  // 	NormalA    += rGeom[0].FastGetSolutionStepValue(NORMAL);
  //     }
  //     if(rGeom[1].Is(FREE_SURFACE)){
  // 	count+=1.0;
  // 	NormalA    += rGeom[1].FastGetSolutionStepValue(NORMAL);
  //     }
  //     if(rGeom[2].Is(FREE_SURFACE)){
  // 	count+=1.0;
  // 	NormalA    += rGeom[2].FastGetSolutionStepValue(NORMAL);
  //     }
  //     if(rGeom[3].Is(FREE_SURFACE)){
  // 	count+=1.0;
  // 	NormalA    += rGeom[3].FastGetSolutionStepValue(NORMAL);
  //     }
  //     if(count!=0){
  // 	NormalMean = NormalA/count;
  //     }else{
  // 	NormalMean = NormalA;
  //     }

  //   }
  // }

  //   NormalAcceleration=NormalMean[0]*MeanAcc[0]+NormalMean[1]*MeanAcc[1]+NormalMean[2]*MeanAcc[2];
    
  //   NormalProjSpatialDefRate=NormalMean[0]*SpatialDefRate[0]*NormalMean[0]+
  //     NormalMean[1]*SpatialDefRate[1]*NormalMean[1]+
  //     NormalMean[2]*SpatialDefRate[2]*NormalMean[2]+
  //     2*NormalMean[0]*SpatialDefRate[3]*NormalMean[1]+
  //     2*NormalMean[0]*SpatialDefRate[4]*NormalMean[2]+
  //     2*NormalMean[1]*SpatialDefRate[5]*NormalMean[2];

}


template < > 
void TwoStepUpdatedLagrangianVPElement<2>::CheckStrain2(MatrixType &SpatialVelocityGrad,
							MatrixType &Fgrad,
							MatrixType &VelDefgrad)
{
  if(fabs(VelDefgrad(0,0)-SpatialVelocityGrad(0,0)*Fgrad(0,0)-SpatialVelocityGrad(0,1)*Fgrad(1,0))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(2a) VDg(0,0)="<<VelDefgrad(0,0)<<" SVG: "<<std::endl;
    std::cout<<SpatialVelocityGrad(0,0)<<" "<<SpatialVelocityGrad(0,1)<<" __ ";
    std::cout<<SpatialVelocityGrad(1,0)<<" "<<SpatialVelocityGrad(1,1)<<" __ "<<std::endl;
  }
  if(fabs(VelDefgrad(0,1)-SpatialVelocityGrad(0,0)*Fgrad(0,1)-SpatialVelocityGrad(0,1)*Fgrad(1,1))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(2b) VDg(0,0)="<<VelDefgrad(0,0)<<" SVG: "<<std::endl;
    std::cout<<SpatialVelocityGrad(0,0)<<" "<<SpatialVelocityGrad(0,1)<<" __ ";
    std::cout<<SpatialVelocityGrad(1,0)<<" "<<SpatialVelocityGrad(1,1)<<" __ "<<std::endl;
  }
  if(fabs(VelDefgrad(1,0)-SpatialVelocityGrad(1,0)*Fgrad(0,0)-SpatialVelocityGrad(1,1)*Fgrad(1,0))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(2c) VDg(0,0)="<<VelDefgrad(0,0)<<" SVG: "<<std::endl;
    std::cout<<SpatialVelocityGrad(0,0)<<" "<<SpatialVelocityGrad(0,1)<<" __ ";
    std::cout<<SpatialVelocityGrad(1,0)<<" "<<SpatialVelocityGrad(1,1)<<" __ "<<std::endl;
  }
  if(fabs(VelDefgrad(1,1)-SpatialVelocityGrad(1,0)*Fgrad(0,1)-SpatialVelocityGrad(1,1)*Fgrad(1,1))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(2d) VDg(0,0)="<<VelDefgrad(0,0)<<" SVG: "<<std::endl;
    std::cout<<SpatialVelocityGrad(0,0)<<" "<<SpatialVelocityGrad(0,1)<<" __ ";
    std::cout<<SpatialVelocityGrad(1,0)<<" "<<SpatialVelocityGrad(1,1)<<" __ "<<std::endl;
  }
}


template < > 
bool TwoStepUpdatedLagrangianVPElement<2>::CheckStrain3(VectorType &SpatialDefRate,
							MatrixType &SpatialVelocityGrad)
{
  bool computeElement=true;
  if(fabs(SpatialDefRate[0]-SpatialVelocityGrad(0,0))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(3a) Sdf[0]="<<SpatialDefRate[0]<<" SVG: "<<SpatialVelocityGrad(0,0)<<std::endl;
    computeElement=false;
  }
  if(fabs(SpatialDefRate[1]-SpatialVelocityGrad(1,1))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(3b) Sdf[0]="<<SpatialDefRate[1]<<" SVG: "<<SpatialVelocityGrad(1,1)<<std::endl;
    computeElement=false;
  }
  if(fabs(SpatialDefRate[2]-0.5*(SpatialVelocityGrad(1,0)+SpatialVelocityGrad(0,1)))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(3c) Sdf[0]="<<SpatialDefRate[2]<<" SVG10: "<<SpatialVelocityGrad(1,0)<<" SVG: "<<SpatialVelocityGrad(0,1)<<std::endl;
    computeElement=false;
  }
  return computeElement;
}

template < > 
bool TwoStepUpdatedLagrangianVPElement<3>::CheckStrain3(VectorType &SpatialDefRate,
							MatrixType &SpatialVelocityGrad)
{
  bool computeElement=true;
  if(fabs(SpatialDefRate[0]-SpatialVelocityGrad(0,0))<0.0000001){
   }else{
    std::cout<<"ERROR IN CHECKSTRAIN(3a)";
    computeElement=false;
  }
  if(fabs(SpatialDefRate[1]-SpatialVelocityGrad(1,1))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(3b)";
    computeElement=false;
  }
  if(fabs(SpatialDefRate[2]-SpatialVelocityGrad(2,2))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(3c)";
    computeElement=false;
  }
  if(fabs(SpatialDefRate[3]-0.5*(SpatialVelocityGrad(1,0)+SpatialVelocityGrad(0,1)))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(3d)";
    computeElement=false;
  }
 if(fabs(SpatialDefRate[4]-0.5*(SpatialVelocityGrad(2,0)+SpatialVelocityGrad(0,2)))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(3e)";
    computeElement=false;
  }
 if(fabs(SpatialDefRate[5]-0.5*(SpatialVelocityGrad(2,1)+SpatialVelocityGrad(1,2)))<0.0000001){
  }else{
    std::cout<<"ERROR IN CHECKSTRAIN(3f)";
    computeElement=false;
  }
  return computeElement;
}


template < > 
void TwoStepUpdatedLagrangianVPElement<3>::CalcVolumetricDefRate(const ShapeFunctionDerivativesType& rDN_DX,
								 double &volumetricDefRate,
								 MatrixType &invGradDef,
								 const double theta)
{
  std::cout<<"TO BE IMPLEMENTED ------- CalcVolumetricDefRate -------"<<std::endl;
  //you can compute the volumetric deformation rate using CalcVolDefRateFromSpatialVelGrad
}


template < > 
void TwoStepUpdatedLagrangianVPElement<3>::CheckStrain1(double &VolumetricDefRate,
							     MatrixType &SpatialVelocityGrad)
{
  std::cout<<"TO BE IMPLEMENTED ------- CheckStrain1 -------"<<std::endl;
}


template < > 
void TwoStepUpdatedLagrangianVPElement<3>::CheckStrain2(MatrixType &SpatialVelocityGrad,
							     MatrixType &Fgrad,
							     MatrixType &VelDefgrad)
{
  std::cout<<"TO BE IMPLEMENTED ------- CheckStrain2 -------"<<std::endl;
}


  template< unsigned int TDim >
bool TwoStepUpdatedLagrangianVPElement<TDim>::CheckSliverElements()
  {
    bool computeElement=true;
    unsigned int sliverNodes=0;
    ShapeFunctionDerivativesArrayType DN_DX;
    Matrix NContainer;
    VectorType GaussWeights;

    this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);
    const unsigned int NumGauss = GaussWeights.size();
    // bool sliver=false;

    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    for (unsigned int n = 0; n < NumNodes; n++)
      {
	if(rGeom[n].Is(INTERFACE)){
	  sliverNodes++;
	}
      }
    if(sliverNodes==NumNodes){
      computeElement=false;
    }
    for (unsigned int g = 0; g < NumGauss; g++)
      {
	const double GaussWeight = GaussWeights[g];
	if(fabs(GaussWeight)<0.000000000001){
	  std::cout<<" THIS IS A SLIVER! I WILL NOT COMPUTE IT "<<GaussWeight<<std::endl;
	computeElement=false;
	}

      }

    return computeElement;
  }


  template< unsigned int TDim >
  double TwoStepUpdatedLagrangianVPElement<TDim>::EquivalentStrainRate(const ShapeFunctionDerivativesType &rDN_DX) const
  {
    const GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    // Calculate Symetric gradient
    MatrixType S = ZeroMatrix(TDim,TDim);
    for (unsigned int n = 0; n < NumNodes; ++n)
      {
        const array_1d<double,3>& rVel = rGeom[n].FastGetSolutionStepValue(VELOCITY,1); // OLD VELOCITY (which is incompressible, unlike the fractional step one)
        for (unsigned int i = 0; i < TDim; ++i)
	  for (unsigned int j = 0; j < TDim; ++j)
	    S(i,j) += 0.5 * ( rDN_DX(n,j) * rVel[i] + rDN_DX(n,i) * rVel[j] );
      }

    // Norm of symetric gradient
    double NormS = 0.0;
    for (unsigned int i = 0; i < TDim; ++i)
      for (unsigned int j = 0; j < TDim; ++j)
	NormS += S(i,j) * S(i,j);

    return std::sqrt(2.0*NormS);
  }




  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::ComputeLumpedMassMatrix(Matrix& rMassMatrix,
									const double Weight,
									double & MeanValue)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    double Coeff=1.0+TDim;
    double Count=0;
    for (SizeType i = 0; i < NumNodes; ++i)
      {
    
	double Mij = Weight/Coeff;

        for ( unsigned int j = 0; j <  TDim; j++ )
	  {
            unsigned int index = i * TDim + j;
            rMassMatrix( index, index ) += Mij;
	    Count+=1.0;
	    MeanValue+=Mij;
	  }

      }
    MeanValue*=1.0/Count;
  }




  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::ComputeMomentumMassTerm(Matrix& rMassMatrix,
									     const ShapeFunctionsType& rN,
									     const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    IndexType FirstRow = 0;
    IndexType FirstCol = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
        for (SizeType j = 0; j < NumNodes; ++j)
	  {
            const double Mij = Weight * rN[i] * rN[j];
            for (SizeType d =  0; d < TDim; ++d)
	      rMassMatrix(FirstRow+d,FirstCol+d) += Mij;
            FirstCol += TDim;
	  }
        FirstRow += TDim;
        FirstCol = 0;
      }
  }

 
  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPElement<TDim>::AddExternalForces(Vector& rRHSVector,
								       const double Density,
								       const array_1d<double,3>& rBodyForce,
								       const ShapeFunctionsType& rN,
								       const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    SizeType FirstRow = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {

	if( this->GetGeometry()[i].SolutionStepsDataHas(VOLUME_ACCELERATION) ){ // it must be checked once at the begining only
	  array_1d<double, 3 >& VolumeAcceleration = this->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_ACCELERATION);
	  // Build RHS
	  for (SizeType d = 0; d < TDim; ++d)
	    {
	      // Body force
	      double RHSi = 0;
	       RHSi =  Density * rN[i] * VolumeAcceleration[d];

	      rRHSVector[FirstRow+d] += Weight * RHSi;
	    }

	}
        FirstRow += TDim;

      }
  }

  template< >
  void TwoStepUpdatedLagrangianVPElement<2>::AddDeviatoricInternalForces(Vector& rRHSVector,const ShapeFunctionDerivativesType& rDN_DX, ElementalVariables& rElementalVariables, const double Weight)
  {
 
  }
  template< >
  void TwoStepUpdatedLagrangianVPElement<3>::AddDeviatoricInternalForces(Vector& rRHSVector,const ShapeFunctionDerivativesType& rDN_DX, ElementalVariables& rElementalVariables, const double Weight)
  {
 
  }

  template< >
  void TwoStepUpdatedLagrangianVPElement<2>::AddInternalForces(Vector& rRHSVector,
								    const ShapeFunctionDerivativesType& rDN_DX,
								    ElementalVariables& rElementalVariables,
								    const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    SizeType FirstRow = 0;

    MatrixType invGradDef=rElementalVariables.InvFgrad;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
	double lagDNXi=rDN_DX(i,0)*invGradDef(0,0)+rDN_DX(i,1)*invGradDef(1,0);
	double lagDNYi=rDN_DX(i,0)*invGradDef(0,1)+rDN_DX(i,1)*invGradDef(1,1);
	// lagDNXi=rDN_DX(i,0);
	// lagDNYi=rDN_DX(i,1);

	rRHSVector[FirstRow]   += -Weight*(lagDNXi*rElementalVariables.UpdatedTotalCauchyStress[0]+
					   lagDNYi*rElementalVariables.UpdatedTotalCauchyStress[2]);

	rRHSVector[FirstRow+1] += -Weight*(lagDNYi*rElementalVariables.UpdatedTotalCauchyStress[1]+
					   lagDNXi*rElementalVariables.UpdatedTotalCauchyStress[2]);


	FirstRow += 2;
      }
  }

  template< >
  void TwoStepUpdatedLagrangianVPElement<3>::AddInternalForces(Vector& rRHSVector,
								    const ShapeFunctionDerivativesType& rDN_DX,
								    ElementalVariables& rElementalVariables,
								    const double Weight)
  {

    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    SizeType FirstRow = 0;

    MatrixType invGradDef=rElementalVariables.InvFgrad;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
	double lagDNXi=rDN_DX(i,0)*invGradDef(0,0)+rDN_DX(i,1)*invGradDef(1,0)+rDN_DX(i,2)*invGradDef(2,0);
	double lagDNYi=rDN_DX(i,0)*invGradDef(0,1)+rDN_DX(i,1)*invGradDef(1,1)+rDN_DX(i,2)*invGradDef(2,1);
	double lagDNZi=rDN_DX(i,0)*invGradDef(0,2)+rDN_DX(i,1)*invGradDef(1,2)+rDN_DX(i,2)*invGradDef(2,2);

	// lagDNXi=rDN_DX(i,0);
	// lagDNYi=rDN_DX(i,1);
	// lagDNZi=rDN_DX(i,2);

	rRHSVector[FirstRow]   += -Weight*(lagDNXi*rElementalVariables.UpdatedTotalCauchyStress[0]+
					   lagDNYi*rElementalVariables.UpdatedTotalCauchyStress[3]+
					   lagDNZi*rElementalVariables.UpdatedTotalCauchyStress[4]);

	rRHSVector[FirstRow+1] += -Weight*(lagDNYi*rElementalVariables.UpdatedTotalCauchyStress[1]+
					   lagDNXi*rElementalVariables.UpdatedTotalCauchyStress[3]+
					   lagDNZi*rElementalVariables.UpdatedTotalCauchyStress[5]);

	rRHSVector[FirstRow+2] += -Weight*(lagDNZi*rElementalVariables.UpdatedTotalCauchyStress[2]+
					   lagDNXi*rElementalVariables.UpdatedTotalCauchyStress[4]+
					   lagDNYi*rElementalVariables.UpdatedTotalCauchyStress[5]);

	FirstRow += 3;
      }


  }



  //template< unsigned int TDim >
  //static const TwoStepUpdatedLagrangianVPElement<TDim>::ShapeFunctionsType TwoStepUpdatedLagrangianVPElement<TDim>::msNg = TwoStepUpdatedLagrangianVPElement<TDim>::InitializeShapeFunctions();

  /*
   * Template class definition (this should allow us to compile the desired template instantiations)
   */

  template class TwoStepUpdatedLagrangianVPElement<2>;
  template class TwoStepUpdatedLagrangianVPElement<3>;

}
