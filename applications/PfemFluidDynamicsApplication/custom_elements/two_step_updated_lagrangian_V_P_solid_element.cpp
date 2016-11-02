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
#include "custom_elements/two_step_updated_lagrangian_V_P_solid_element.h"
#include "includes/cfd_variables.h" 

namespace Kratos {

  /*
   * public TwoStepUpdatedLagrangianVPSolidElement<TDim> functions
   */
  
template< unsigned int TDim >
Element::Pointer TwoStepUpdatedLagrangianVPSolidElement<TDim>::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
  // return Element::Pointer( BaseType::Clone(NewId,rThisNodes) );
  TwoStepUpdatedLagrangianVPSolidElement NewElement(NewId, this->GetGeometry().Create( rThisNodes ), this->pGetProperties() );
  

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

  return Element::Pointer( new TwoStepUpdatedLagrangianVPSolidElement(NewElement) );
}


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPSolidElement<TDim>::Initialize()
  {
    KRATOS_TRY;

    // LargeDisplacementElement::Initialize();
    const GeometryType& rGeom = this->GetGeometry();
    SizeType integration_points_number = rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_1);
    // const unsigned int dimension       = rGeom.WorkingSpaceDimension();

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
	this->mCurrentTotalCauchyStress[PointNumber] = ZeroVector(voigtsize);
	this->mCurrentDeviatoricCauchyStress[PointNumber] = ZeroVector(voigtsize);
	this->mUpdatedTotalCauchyStress[PointNumber] = ZeroVector(voigtsize);
	this->mUpdatedDeviatoricCauchyStress[PointNumber] = ZeroVector(voigtsize);
      }

    KRATOS_CATCH( "" );
  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPSolidElement<TDim>::InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo)
  {
   KRATOS_TRY;

    const GeometryType& rGeom = this->GetGeometry();
    SizeType integration_points_number = rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_1);

    for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
      {
        this->UpdateCauchyStress(PointNumber);
      }

    KRATOS_CATCH( "" );

  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPSolidElement<TDim>::InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo)
  {

  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPSolidElement<TDim>::ComputeMaterialParameters(double& DeviatoricCoeff,
									       double& VolumetricCoeff,
									       double timeStep,
									       const ShapeFunctionsType& N)
  {

  
    // double YoungModulus=this->pGetProperties()[YOUNG_MODULUS];
    // double PoissonRatio=this->pGetProperties()[POISSON_RATIO];
    std::cout<<"I am setting from inside YoungModulus = 100000000.0 ... work in progress...."<<std::endl;
    const double YoungModulus = 100000000.0;
    const double PoissonRatio = 0.495;
 
    // FirstLame = timeStep*PoissonRatio*YoungModulus/((1.0+PoissonRatio)*(1.0-2.0*PoissonRatio)); 
    DeviatoricCoeff = timeStep*YoungModulus/(1.0+PoissonRatio)*0.5;
    // BulkModulus = FirstLame + 2.0*SecondLame/3.0;
    VolumetricCoeff = timeStep*PoissonRatio*YoungModulus/((1.0+PoissonRatio)*(1.0-2.0*PoissonRatio)) + 2.0*DeviatoricCoeff/3.0;
  }



  template< unsigned int TDim >
  int TwoStepUpdatedLagrangianVPSolidElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo)
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
    if(MESH_VELOCITY.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"MESH_VELOCITY Key is 0. Check that the application was correctly registered.","");
    if(FRACT_VEL.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"FRACT_VEL Key is 0. Check that the application was correctly registered.","");
    if(PRESSURE_OLD_IT.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"PRESSURE_OLD_IT Key is 0. Check that the application was correctly registered.","");
    if(NODAL_AREA.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"NODAL_AREA Key is 0. Check that the application was correctly registered.","");
    if(CONV_PROJ.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"CONV_PROJ Key is 0. Check that the application was correctly registered.","");
    if(PRESS_PROJ.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"PRESS_PROJ Key is 0. Check that the application was correctly registered.","");
    if(DIVPROJ.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"DIVPROJ Key is 0. Check that the application was correctly registered.","");
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
        if(this->GetGeometry()[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing MESH_VELOCITY variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(FRACT_VEL) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing FRACT_VEL variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE_OLD_IT) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE_OLD_IT variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(NODAL_AREA) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing NODAL_AREA variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(CONV_PROJ) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing CONV_PROJ variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(PRESS_PROJ) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESS_PROJ variable on solution step data for node ",this->GetGeometry()[i].Id());
        if(this->GetGeometry()[i].SolutionStepsDataHas(DIVPROJ) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing DIVPROJ variable on solution step data for node ",this->GetGeometry()[i].Id());
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
  void TwoStepUpdatedLagrangianVPSolidElement<2>::AddCompleteTangentTerm(ElementalVariables & rElementalVariables,
									 MatrixType& rDampingMatrix,
									 const ShapeFunctionDerivativesType& rDN_DX,
									 const double secondLame,
									 const double bulkModulus,
									 const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    MatrixType invGradDef=rElementalVariables.InvFgrad;
    double theta=0.5;

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
 

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPSolidElement<TDim>::ComputeBulkMatrixForPressureVelLump(Matrix& BulkVelMatrix,
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
  void TwoStepUpdatedLagrangianVPSolidElement<TDim>::ComputeBulkMatrixForPressureAccLump(Matrix& BulkAccMatrix,
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


 

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPSolidElement<TDim>::ComputeBulkMatrixForPressureVel(Matrix& BulkVelMatrix,
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
  void TwoStepUpdatedLagrangianVPSolidElement<TDim>::ComputeBulkMatrixForPressureAcc(Matrix& BulkAccMatrix,
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
 void TwoStepUpdatedLagrangianVPSolidElement<TDim>::AddStabilizationMatrixLHS(MatrixType& rLeftHandSideMatrix,
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
 void TwoStepUpdatedLagrangianVPSolidElement<TDim>::ComputeStabLaplacianMatrix(MatrixType& StabLaplacianMatrix,
									       const ShapeFunctionDerivativesType& rDN_DX,									   
									       const double Weight
)
										
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
	    StabLaplacianMatrix(i,j) += Weight * Lij;
	  }
      }
  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPSolidElement<TDim>::AddStabilizationNodalTermsLHS(MatrixType& rLeftHandSideMatrix,
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
  void TwoStepUpdatedLagrangianVPSolidElement<TDim>::AddStabilizationNodalTermsRHS(VectorType& rRightHandSideVector,
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
	  // RHSi += rDN_DX(i,d) * Tau * ( - OldPressureGradient[d]  );
	  RHSi += - rDN_DX(i,d) * Tau * ( Density * VolumeAcceleration[d] );

	}
    }
    rRightHandSideVector[i] += Weight * RHSi;

  }


  template<>
  void TwoStepUpdatedLagrangianVPSolidElement<3>::AddCompleteTangentTerm(ElementalVariables & rElementalVariables,
									 MatrixType& rDampingMatrix,
									 const ShapeFunctionDerivativesType& rDN_DX,
									 const double secondLame,
									 const double bulkModulus,
									 const double Weight){

   const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    MatrixType invGradDef=rElementalVariables.InvFgrad;
    double theta=0.5;

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



  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPSolidElement<TDim>::CalculateTauFIC(double& Tau,
								     double ElemSize,
								     const array_1d< double, 3 > & rAdvVel,
								     const double Density,
								     const double Viscosity,
								     const ProcessInfo& rCurrentProcessInfo)
  {
    const double DeltaTime = rCurrentProcessInfo.GetValue(DELTA_TIME);
    const double TimeFactor = rCurrentProcessInfo.GetValue(DYNAMIC_TAU);

    // Compute mean advective velocity norm
    double AdvVelNorm = 0.0;
    for (unsigned int d = 0; d < TDim; ++d)
      AdvVelNorm += rAdvVel[d] * rAdvVel[d];

    AdvVelNorm = sqrt(AdvVelNorm);

    Tau = 1.0 / (Density * ( TimeFactor / DeltaTime + 2.0 * AdvVelNorm / ElemSize) + 4.0 * Viscosity / (ElemSize * ElemSize) ); 
    Tau = 1.0 / (2.0 * Density * AdvVelNorm / ElemSize +  8.0 * Viscosity / (ElemSize * ElemSize) ); //viscosity already includes density
    Tau=0;

  }




template< unsigned int TDim>
bool TwoStepUpdatedLagrangianVPSolidElement<TDim>::CalcMechanicsUpdated(ElementalVariables & rElementalVariables,
									const ProcessInfo& rCurrentProcessInfo,
									unsigned int g,
									const ShapeFunctionsType& N)
{

  bool computeElement=false;
  computeElement=this->CalcStrainRateUpdated(rElementalVariables,rCurrentProcessInfo,g);
  const double TimeStep=0.5/rCurrentProcessInfo[BDF_COEFFICIENTS][2];
  this->CalcElasticPlasticCauchySplitted(rElementalVariables,TimeStep,g,N);
  return computeElement;
} 




template <  unsigned int TDim> 
void TwoStepUpdatedLagrangianVPSolidElement<TDim>:: InitializeElementalVariables(ElementalVariables & rElementalVariables)
 {
    unsigned int voigtsize  = 3;
    if( TDim == 3 )
    {
        voigtsize  = 6;
    }
    rElementalVariables.voigtsize=voigtsize;

    rElementalVariables.DetFgrad=1;
    rElementalVariables.DetFgradVel=1;
    rElementalVariables.DeviatoricInvariant=1;
    rElementalVariables.VolumetricDefRate=1;
    rElementalVariables.SpatialDefRate.resize(voigtsize);
    rElementalVariables.MDGreenLagrangeMaterial.resize(voigtsize);
    rElementalVariables.Fgrad.resize(TDim,TDim);
    rElementalVariables.InvFgrad.resize(TDim,TDim);
    rElementalVariables.FgradVel.resize(TDim,TDim);
    rElementalVariables.InvFgradVel.resize(TDim,TDim);
    rElementalVariables.SpatialVelocityGrad.resize(TDim,TDim);

    rElementalVariables.MeanPressure=0;
    rElementalVariables.CurrentTotalCauchyStress.resize(voigtsize);
    rElementalVariables.UpdatedTotalCauchyStress.resize(voigtsize);
    rElementalVariables.CurrentDeviatoricCauchyStress.resize(voigtsize);
    rElementalVariables.UpdatedDeviatoricCauchyStress.resize(voigtsize);

 }


template < > 
void TwoStepUpdatedLagrangianVPSolidElement<2>:: CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables,double TimeStep, unsigned int g,const ShapeFunctionsType& N)
{

  rElementalVariables.CurrentTotalCauchyStress=this->mCurrentTotalCauchyStress[g];
  rElementalVariables.CurrentDeviatoricCauchyStress=this->mCurrentDeviatoricCauchyStress[g];

  double CurrSecondLame  = 0;
  double CurrBulkModulus = 0;

  this->ComputeMaterialParameters(CurrSecondLame,CurrBulkModulus,TimeStep,N);
 
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

  sigmaDev_xx+=rElementalVariables.CurrentDeviatoricCauchyStress[0];
  sigmaDev_yy+=rElementalVariables.CurrentDeviatoricCauchyStress[1];
  sigmaDev_xy+=rElementalVariables.CurrentDeviatoricCauchyStress[2];

  sigmaTot_xx+=rElementalVariables.CurrentTotalCauchyStress[0];
  sigmaTot_yy+=rElementalVariables.CurrentTotalCauchyStress[1];
  sigmaTot_xy+=rElementalVariables.CurrentTotalCauchyStress[2];

  // sigmaTot_xx= sigmaDev_xx + rElementalVariables.MeanPressure;
  // sigmaTot_yy= sigmaDev_yy + rElementalVariables.MeanPressure;
  // sigmaTot_xy= sigmaDev_xy;

  // sigmaDev_xx= sigmaTot_xx - rElementalVariables.MeanPressure;
  // sigmaDev_yy= sigmaTot_yy - rElementalVariables.MeanPressure;
  // sigmaDev_xy= sigmaTot_xy;

  rElementalVariables.UpdatedDeviatoricCauchyStress[0]=sigmaDev_xx;
  rElementalVariables.UpdatedDeviatoricCauchyStress[1]=sigmaDev_yy;
  rElementalVariables.UpdatedDeviatoricCauchyStress[2]=sigmaDev_xy;

  rElementalVariables.UpdatedTotalCauchyStress[0]=sigmaTot_xx;
  rElementalVariables.UpdatedTotalCauchyStress[1]=sigmaTot_yy;
  rElementalVariables.UpdatedTotalCauchyStress[2]=sigmaTot_xy;

  this->mCurrentTotalCauchyStress[g]=rElementalVariables.CurrentTotalCauchyStress;
  this->mUpdatedTotalCauchyStress[g]=rElementalVariables.UpdatedTotalCauchyStress;
  this->mCurrentDeviatoricCauchyStress[g]=rElementalVariables.CurrentDeviatoricCauchyStress;
  this->mUpdatedDeviatoricCauchyStress[g]=rElementalVariables.UpdatedDeviatoricCauchyStress;

}


template < > 
void TwoStepUpdatedLagrangianVPSolidElement<3>:: CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables, double TimeStep, unsigned int g,const ShapeFunctionsType& N)
{


  rElementalVariables.CurrentTotalCauchyStress=this->mCurrentTotalCauchyStress[g];
  rElementalVariables.CurrentDeviatoricCauchyStress=this->mCurrentDeviatoricCauchyStress[g];

  double CurrSecondLame  = 0;
  double CurrBulkModulus = 0;

  this->ComputeMaterialParameters(CurrSecondLame,CurrBulkModulus,TimeStep,N);
 
  double CurrFirstLame   = 0;
  CurrFirstLame  = CurrBulkModulus - 2.0*CurrSecondLame/3.0;
  
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


  sigmaDev_xx+=rElementalVariables.CurrentDeviatoricCauchyStress[0];
  sigmaDev_yy+=rElementalVariables.CurrentDeviatoricCauchyStress[1];
  sigmaDev_zz+=rElementalVariables.CurrentDeviatoricCauchyStress[2];
  sigmaDev_xy+=rElementalVariables.CurrentDeviatoricCauchyStress[3];
  sigmaDev_xz+=rElementalVariables.CurrentDeviatoricCauchyStress[4];
  sigmaDev_yz+=rElementalVariables.CurrentDeviatoricCauchyStress[5];

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

  this->mCurrentTotalCauchyStress[g]=rElementalVariables.CurrentTotalCauchyStress;
  this->mUpdatedTotalCauchyStress[g]=rElementalVariables.UpdatedTotalCauchyStress;
  this->mCurrentDeviatoricCauchyStress[g]=rElementalVariables.CurrentDeviatoricCauchyStress;
  this->mUpdatedDeviatoricCauchyStress[g]=rElementalVariables.UpdatedDeviatoricCauchyStress;

}



template < unsigned int TDim > 
void TwoStepUpdatedLagrangianVPSolidElement<TDim>:: UpdateCauchyStress(unsigned int g)
 {
   this->mCurrentTotalCauchyStress[g]=this->mUpdatedTotalCauchyStress[g];
   this->mCurrentDeviatoricCauchyStress[g]=this->mUpdatedDeviatoricCauchyStress[g];
 }





  template class TwoStepUpdatedLagrangianVPSolidElement<2>;
  template class TwoStepUpdatedLagrangianVPSolidElement<3>;

}
