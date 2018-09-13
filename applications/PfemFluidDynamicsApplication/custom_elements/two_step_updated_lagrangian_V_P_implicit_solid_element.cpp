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
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_solid_element.h"
#include "includes/cfd_variables.h"

namespace Kratos {

  /*
   * public TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim> functions
   */

  template< unsigned int TDim >
  Element::Pointer TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim>::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {
    // return Element::Pointer( BaseType::Clone(NewId,rThisNodes) );
    TwoStepUpdatedLagrangianVPImplicitSolidElement NewElement(NewId, this->GetGeometry().Create( rThisNodes ), this->pGetProperties() );

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

    return Element::Pointer( new TwoStepUpdatedLagrangianVPImplicitSolidElement(NewElement) );
  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim>::Initialize()
  {
    KRATOS_TRY;

    // LargeDisplacementElement::Initialize();
    const GeometryType& rGeom = this->GetGeometry();
    SizeType integration_points_number = rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_1);
    // SizeType integration_points_number = rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_4);
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
	this->mCurrentTotalCauchyStress[PointNumber] = ZeroVector(voigtsize);
	this->mCurrentDeviatoricCauchyStress[PointNumber] = ZeroVector(voigtsize);
	this->mUpdatedTotalCauchyStress[PointNumber] = ZeroVector(voigtsize);
	this->mUpdatedDeviatoricCauchyStress[PointNumber] = ZeroVector(voigtsize);
      }

    KRATOS_CATCH( "" );
  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim>::InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo)
  {
    KRATOS_TRY;

    const GeometryType& rGeom = this->GetGeometry();
    SizeType integration_points_number = rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_1);
    // SizeType integration_points_number = rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_4);

    for ( unsigned int PointNumber = 0; PointNumber < integration_points_number; PointNumber++ )
      {
        this->UpdateCauchyStress(PointNumber,rCurrentProcessInfo);
      }

    KRATOS_CATCH( "" );

  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim>::InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo)
  {

  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim>::ComputeMaterialParameters(double& Density,
										       double& DeviatoricCoeff,
										       double& VolumetricCoeff,
										       ProcessInfo &currentProcessInfo,
										       ElementalVariables& rElementalVariables)
  {

    Density=this->GetProperties()[DENSITY];
    double YoungModulus=this->GetProperties()[YOUNG_MODULUS];
    double PoissonRatio=this->GetProperties()[POISSON_RATIO];
    double timeStep=currentProcessInfo[DELTA_TIME];


    // this->EvaluateInPoint(YoungModulus,YOUNG_MODULUS,N);
    // this->EvaluateInPoint(PoissonRatio,POISSON_RATIO,N);
    // this->EvaluateInPoint(Density,DENSITY,N);


    // YoungModulus=10000000;//falling cylinder
    // PoissonRatio=0.35;//falling cylinder
    // YoungModulus=1000000;//dam break fsi
    // PoissonRatio=0;//dam break fsi
    // YoungModulus=100000000;
    // PoissonRatio=0;

    // FirstLame = timeStep*PoissonRatio*YoungModulus/((1.0+PoissonRatio)*(1.0-2.0*PoissonRatio));
    DeviatoricCoeff = timeStep*YoungModulus/(1.0+PoissonRatio)*0.5;
    // BulkModulus = FirstLame + 2.0*SecondLame/3.0;
    VolumetricCoeff = timeStep*PoissonRatio*YoungModulus/((1.0+PoissonRatio)*(1.0-2.0*PoissonRatio)) + 2.0*DeviatoricCoeff/3.0;

    this->mMaterialDeviatoricCoefficient=DeviatoricCoeff;
    this->mMaterialVolumetricCoefficient=VolumetricCoeff;
    this->mMaterialDensity=Density;

  }



  template< unsigned int TDim >
  int TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo)
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


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim>::ComputeBulkMatrixForPressureVelLump(Matrix& BulkVelMatrix,
												 const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    double coeff=1.0+TDim;
    // coeff=6.0;
    if(TDim==2 && NumNodes ==6){
      double Mij = Weight/57.0;
      double consistent=1.0;
      for (SizeType i = 0; i < NumNodes; ++i)
	{
	  if(i<3){
	    consistent=coeff;
	  }else{
	    consistent=16.0;
	  }

	  BulkVelMatrix(i,i) +=  Mij*consistent;

	}
    }else{
      for (SizeType i = 0; i < NumNodes; ++i)
	{
	  // LHS contribution
	  double Mij  = Weight /coeff;
	  BulkVelMatrix(i,i) +=  Mij;
	}
      if(NumNodes>4){
	std::cout<<"ComputeBulkMatrixForPressureVelLump 3D not yet implemented!"<<std::endl;
      }
    }


  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim>::ComputeBulkMatrixForPressureVel(Matrix& BulkVelMatrix,
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



  // template< unsigned int TDim>
  // bool TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim>::CalcMechanicsUpdated(ElementalVariables & rElementalVariables,
  // 									const ProcessInfo& rCurrentProcessInfo,
  // 									const ShapeFunctionDerivativesType& rDN_DX,
  // 									unsigned int g)
  // {

  //   double theta=this->GetThetaMomentum();
  //   bool computeElement=this->CalcStrainRate(rElementalVariables,rCurrentProcessInfo,rDN_DX,theta);
  //   const double TimeStep=rCurrentProcessInfo[DELTA_TIME];
  //   this->CalcElasticPlasticCauchySplitted(rElementalVariables,TimeStep,g);
  //   return computeElement;
  // }


  template<>
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<2>::GetPositions(Vector& rValues,const ProcessInfo& rCurrentProcessInfo,const double theta)
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
    	// rValues[Index++] = rGeom[i].X0();
    	// rValues[Index++] = rGeom[i].Y0();
      }
  }


  template<>
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<3>::GetPositions(Vector& rValues,const ProcessInfo& rCurrentProcessInfo,const double theta)
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
	// rValues[Index++] = rGeom[i].X0();
        // rValues[Index++] = rGeom[i].Y0();
        // rValues[Index++] = rGeom[i].Z0();
      }
  }



  template <  unsigned int TDim>
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim>:: InitializeElementalVariables(ElementalVariables & rElementalVariables)
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
    rElementalVariables.EquivalentStrainRate=1;
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
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<2>:: CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables, double TimeStep, unsigned int g)
  {

    rElementalVariables.CurrentTotalCauchyStress=this->mCurrentTotalCauchyStress[g];
    rElementalVariables.CurrentDeviatoricCauchyStress=this->mCurrentDeviatoricCauchyStress[g];

    double CurrSecondLame  = this->mMaterialDeviatoricCoefficient;
    // double CurrBulkModulus = this->mMaterialVolumetricCoefficient;
    // double CurrFirstLame  =CurrBulkModulus - 2.0*CurrSecondLame/3.0;

    double DefX=rElementalVariables.SpatialDefRate[0];
    double DefY=rElementalVariables.SpatialDefRate[1];
    double DefXY=rElementalVariables.SpatialDefRate[2];

    double DefVol=rElementalVariables.VolumetricDefRate;

    double sigmaDev_xx= 2*CurrSecondLame*(DefX - DefVol/3.0);
    double sigmaDev_yy= 2*CurrSecondLame*(DefY - DefVol/3.0);
    double sigmaDev_xy= 2*CurrSecondLame*DefXY;

    // double sigmaTot_xx= CurrFirstLame*DefVol + 2.0*CurrSecondLame*DefX;
    // double sigmaTot_yy= CurrFirstLame*DefVol + 2.0*CurrSecondLame*DefY;
    // double sigmaTot_xy= 2.0*CurrSecondLame*DefXY;

    // sigmaTot_xx=rElementalVariables.CurrentTotalCauchyStress[0] + rElementalVariables.MeanPressure + sigmaDev_xx;
    // sigmaTot_yy=rElementalVariables.CurrentTotalCauchyStress[1] + rElementalVariables.MeanPressure + sigmaDev_yy;
    // sigmaTot_xy=rElementalVariables.CurrentTotalCauchyStress[2] + sigmaDev_xy;

    sigmaDev_xx+=rElementalVariables.CurrentDeviatoricCauchyStress[0];
    sigmaDev_yy+=rElementalVariables.CurrentDeviatoricCauchyStress[1];
    sigmaDev_xy+=rElementalVariables.CurrentDeviatoricCauchyStress[2];

    double sigmaTot_xx= sigmaDev_xx + rElementalVariables.MeanPressure;
    double sigmaTot_yy= sigmaDev_yy + rElementalVariables.MeanPressure;
    double sigmaTot_xy= sigmaDev_xy;

    // sigmaTot_xx+=rElementalVariables.CurrentTotalCauchyStress[0];
    // sigmaTot_yy+=rElementalVariables.CurrentTotalCauchyStress[1];
    // sigmaTot_xy+=rElementalVariables.CurrentTotalCauchyStress[2];

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
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<3>:: CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables, double TimeStep, unsigned int g)
  {

    rElementalVariables.CurrentTotalCauchyStress=this->mCurrentTotalCauchyStress[g];
    rElementalVariables.CurrentDeviatoricCauchyStress=this->mCurrentDeviatoricCauchyStress[g];

    double CurrSecondLame  = this->mMaterialDeviatoricCoefficient;
    // double CurrBulkModulus = this->mMaterialVolumetricCoefficient;
    // double CurrFirstLame = CurrBulkModulus - 2.0*CurrSecondLame/3.0;

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

    // double sigmaTot_xx= CurrFirstLame*DefVol + 2*CurrSecondLame*DefX;
    // double sigmaTot_yy= CurrFirstLame*DefVol + 2*CurrSecondLame*DefY;
    // double sigmaTot_zz= CurrFirstLame*DefVol + 2*CurrSecondLame*DefZ;
    // double sigmaTot_xy= 2*CurrSecondLame*DefXY;
    // double sigmaTot_xz= 2*CurrSecondLame*DefXZ;
    // double sigmaTot_yz= 2*CurrSecondLame*DefYZ;

    sigmaDev_xx+=rElementalVariables.CurrentDeviatoricCauchyStress[0];
    sigmaDev_yy+=rElementalVariables.CurrentDeviatoricCauchyStress[1];
    sigmaDev_zz+=rElementalVariables.CurrentDeviatoricCauchyStress[2];
    sigmaDev_xy+=rElementalVariables.CurrentDeviatoricCauchyStress[3];
    sigmaDev_xz+=rElementalVariables.CurrentDeviatoricCauchyStress[4];
    sigmaDev_yz+=rElementalVariables.CurrentDeviatoricCauchyStress[5];

    double sigmaTot_xx= sigmaDev_xx + rElementalVariables.MeanPressure;
    double sigmaTot_yy= sigmaDev_yy + rElementalVariables.MeanPressure;
    double sigmaTot_zz= sigmaDev_zz + rElementalVariables.MeanPressure;
    double sigmaTot_xy= sigmaDev_xy;
    double sigmaTot_xz= sigmaDev_xz;
    double sigmaTot_yz= sigmaDev_yz;

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
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim>:: UpdateCauchyStress(unsigned int g,ProcessInfo &rCurrentProcessInfo)
  {

    double theta=this->GetThetaContinuity();
    ElementalVariables rElementalVariables;
    this->InitializeElementalVariables(rElementalVariables);
    ShapeFunctionDerivativesArrayType DN_DX;
    Matrix NContainer;
    VectorType GaussWeights;
    this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);
    const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];
    // bool computeElement=this->CalcStrainRate(rElementalVariables,rCurrentProcessInfo,rDN_DX,theta);
    bool computeElement=this->CalcCompleteStrainRate(rElementalVariables,rCurrentProcessInfo,rDN_DX,theta);
    const double TimeStep=rCurrentProcessInfo[DELTA_TIME];
    if(computeElement==true){
      this->CalcElasticPlasticCauchySplitted(rElementalVariables,TimeStep,g);
    }

    // ElementalVariables rElementalVariables;
    // this->InitializeElementalVariables(rElementalVariables);

    // double theta=0;
    // bool computeElement=this->CalcStrainRate(rElementalVariables,rCurrentProcessInfo,g,theta);

    //     // //SIGMA(2D)=(F·S·Ft)/J --> already checked -> ok!
    // double SigmaXX= (  rElementalVariables.Fgrad(0,0)*rElementalVariables.UpdatedTotalCauchyStress[0]*rElementalVariables.Fgrad(0,0) +
    // 		    2*rElementalVariables.Fgrad(0,0)*rElementalVariables.UpdatedTotalCauchyStress[2]*rElementalVariables.Fgrad(0,1) +
    // 		      rElementalVariables.Fgrad(0,1)*rElementalVariables.UpdatedTotalCauchyStress[1]*rElementalVariables.Fgrad(0,1))/rElementalVariables.DetFgrad ;

    // double SigmaYY= (  rElementalVariables.Fgrad(1,0)*rElementalVariables.UpdatedTotalCauchyStress[0]*rElementalVariables.Fgrad(1,0) +
    // 		    2*rElementalVariables.Fgrad(1,1)*rElementalVariables.UpdatedTotalCauchyStress[2]*rElementalVariables.Fgrad(1,0) +
    // 		      rElementalVariables.Fgrad(1,1)*rElementalVariables.UpdatedTotalCauchyStress[1]*rElementalVariables.Fgrad(1,1))/rElementalVariables.DetFgrad ;

    // double SigmaXY= (rElementalVariables.Fgrad(0,0)*rElementalVariables.UpdatedTotalCauchyStress[0]*rElementalVariables.Fgrad(1,0) +
    // 		    rElementalVariables.Fgrad(0,0)*rElementalVariables.UpdatedTotalCauchyStress[2]*rElementalVariables.Fgrad(1,1) +
    // 		    rElementalVariables.Fgrad(0,1)*rElementalVariables.UpdatedTotalCauchyStress[2]*rElementalVariables.Fgrad(1,0) +
    // 		    rElementalVariables.Fgrad(0,1)*rElementalVariables.UpdatedTotalCauchyStress[1]*rElementalVariables.Fgrad(1,1))/rElementalVariables.DetFgrad;

    // rElementalVariables.UpdatedTotalCauchyStress[0]=SigmaXX;
    // rElementalVariables.UpdatedTotalCauchyStress[1]=SigmaYY;
    // rElementalVariables.UpdatedTotalCauchyStress[2]=SigmaXY;

    // this->mUpdatedTotalCauchyStress[g]=rElementalVariables.UpdatedTotalCauchyStress;

    //     // //SIGMA(2D)=(F·S·Ft)/J --> already checked -> ok!
    // SigmaXX= (  rElementalVariables.Fgrad(0,0)*rElementalVariables.UpdatedDeviatoricCauchyStress[0]*rElementalVariables.Fgrad(0,0) +
    // 	       2*rElementalVariables.Fgrad(0,0)*rElementalVariables.UpdatedDeviatoricCauchyStress[2]*rElementalVariables.Fgrad(0,1) +
    // 	       rElementalVariables.Fgrad(0,1)*rElementalVariables.UpdatedDeviatoricCauchyStress[1]*rElementalVariables.Fgrad(0,1))/rElementalVariables.DetFgrad ;

    // SigmaYY= (  rElementalVariables.Fgrad(1,0)*rElementalVariables.UpdatedDeviatoricCauchyStress[0]*rElementalVariables.Fgrad(1,0) +
    // 	       2*rElementalVariables.Fgrad(1,1)*rElementalVariables.UpdatedDeviatoricCauchyStress[2]*rElementalVariables.Fgrad(1,0) +
    // 	       rElementalVariables.Fgrad(1,1)*rElementalVariables.UpdatedDeviatoricCauchyStress[1]*rElementalVariables.Fgrad(1,1))/rElementalVariables.DetFgrad ;

    // SigmaXY= (rElementalVariables.Fgrad(0,0)*rElementalVariables.UpdatedDeviatoricCauchyStress[0]*rElementalVariables.Fgrad(1,0) +
    // 	     rElementalVariables.Fgrad(0,0)*rElementalVariables.UpdatedDeviatoricCauchyStress[2]*rElementalVariables.Fgrad(1,1) +
    // 	     rElementalVariables.Fgrad(0,1)*rElementalVariables.UpdatedDeviatoricCauchyStress[2]*rElementalVariables.Fgrad(1,0) +
    // 	     rElementalVariables.Fgrad(0,1)*rElementalVariables.UpdatedDeviatoricCauchyStress[1]*rElementalVariables.Fgrad(1,1))/rElementalVariables.DetFgrad;

    // rElementalVariables.UpdatedDeviatoricCauchyStress[0]=SigmaXX;
    // rElementalVariables.UpdatedDeviatoricCauchyStress[1]=SigmaYY;
    // rElementalVariables.UpdatedDeviatoricCauchyStress[2]=SigmaXY;

    // this->mUpdatedDeviatoricCauchyStress[g]=rElementalVariables.UpdatedDeviatoricCauchyStress;

    // std::cout<<"  UpdateCauchyStress ";
    this->mCurrentTotalCauchyStress[g]=this->mUpdatedTotalCauchyStress[g];
    this->mCurrentDeviatoricCauchyStress[g]=this->mUpdatedDeviatoricCauchyStress[g];


  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitSolidElement<TDim>::CalculateLocalContinuityEqForPressure(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    // Check sizes and initialize
    if( rLeftHandSideMatrix.size1() != NumNodes )
      rLeftHandSideMatrix.resize(NumNodes,NumNodes);

    rLeftHandSideMatrix = ZeroMatrix(NumNodes,NumNodes);

    if( rRightHandSideVector.size() != NumNodes )
      rRightHandSideVector.resize(NumNodes);

    rRightHandSideVector = ZeroVector(NumNodes);

    // MatrixType BulkVelMatrix = ZeroMatrix(NumNodes,NumNodes);

    // Shape functions and integration points
    ShapeFunctionDerivativesArrayType DN_DX;
    Matrix NContainer;
    VectorType GaussWeights;
    this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);
    const unsigned int NumGauss = GaussWeights.size();

    // const double TimeStep=rCurrentProcessInfo[DELTA_TIME];
    double theta=this->GetThetaContinuity();

    ElementalVariables rElementalVariables;
    this->InitializeElementalVariables(rElementalVariables);

    // double Density  = mMaterialDensity;
    // double DeviatoricCoeff = mMaterialVolumetricCoefficient;
    double VolumetricCoeff = this->mMaterialVolumetricCoefficient;

    double totalVolume=0;

    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; g++)
      {
	const double GaussWeight = GaussWeights[g];
	totalVolume+=GaussWeight;
	const ShapeFunctionsType& N = row(NContainer,g);
	const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];
	// bool computeElement=this->CalcStrainRate(rElementalVariables,rCurrentProcessInfo,rDN_DX,theta);
	bool computeElement=this->CalcCompleteStrainRate(rElementalVariables,rCurrentProcessInfo,rDN_DX,theta);
	if(computeElement==true){
	  // double BulkCoeff =GaussWeight/(VolumetricCoeff);
	  // this->ComputeBulkMatrixForPressureVel(BulkVelMatrix,N,BulkCoeff);

	  for (SizeType i = 0; i < NumNodes; ++i)
	    {
	      // RHS contribution
	      // Velocity divergence
	      double RHSi =  N[i] * rElementalVariables.VolumetricDefRate;
	      rRightHandSideVector[i] += GaussWeight * RHSi;
	    }

	}

      }

    MatrixType BulkVelMatrixLump = ZeroMatrix(NumNodes,NumNodes);
    double lumpedBulkCoeff =totalVolume/(VolumetricCoeff);
    this->ComputeBulkMatrixForPressureVelLump(BulkVelMatrixLump,lumpedBulkCoeff);

    rLeftHandSideMatrix+=BulkVelMatrixLump;
    // rLeftHandSideMatrix+=BulkVelMatrix;

    VectorType UpdatedPressure = ZeroVector(NumNodes);
    VectorType CurrentPressure = ZeroVector(NumNodes);;

    this->GetPressureValues(UpdatedPressure,0);
    this->GetPressureValues(CurrentPressure,1);

    VectorType DeltaPressure = UpdatedPressure-CurrentPressure;

    rRightHandSideVector -= prod(BulkVelMatrixLump,DeltaPressure);
    // rRightHandSideVector -= prod(BulkVelMatrix,DeltaPressure);

  }


  template class TwoStepUpdatedLagrangianVPImplicitSolidElement<2>;
  template class TwoStepUpdatedLagrangianVPImplicitSolidElement<3>;

}
