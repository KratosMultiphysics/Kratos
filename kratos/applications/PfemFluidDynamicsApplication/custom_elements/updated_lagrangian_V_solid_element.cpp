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
#include "custom_elements/updated_lagrangian_V_solid_element.h"
#include "includes/cfd_variables.h" 

namespace Kratos {

  /*
   * public UpdatedLagrangianVSolidElement<TDim> functions
   */
  
template< unsigned int TDim >
Element::Pointer UpdatedLagrangianVSolidElement<TDim>::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
  // return Element::Pointer( BaseType::Clone(NewId,rThisNodes) );
  UpdatedLagrangianVSolidElement NewElement(NewId, this->GetGeometry().Create( rThisNodes ), this->pGetProperties() );
 
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

  return Element::Pointer( new UpdatedLagrangianVSolidElement(NewElement) );
}




template< unsigned int TDim>
bool UpdatedLagrangianVSolidElement<TDim>::CalcMechanicsUpdated(ElementalVariables & rElementalVariables,
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




template < > 
void UpdatedLagrangianVSolidElement<2>:: CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables,double TimeStep, unsigned int g,const ShapeFunctionsType& N)
{

  rElementalVariables.CurrentTotalCauchyStress=this->mCurrentTotalCauchyStress[g];
  rElementalVariables.CurrentDeviatoricCauchyStress=this->mCurrentDeviatoricCauchyStress[g];

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

  sigmaDev_xx+=rElementalVariables.CurrentDeviatoricCauchyStress[0];
  sigmaDev_yy+=rElementalVariables.CurrentDeviatoricCauchyStress[1];
  sigmaDev_xy+=rElementalVariables.CurrentDeviatoricCauchyStress[2];
  // std::cout<<"  sxx:"<<rElementalVariables.CurrentDeviatoricCauchyStress[0];
  sigmaTot_xx+=rElementalVariables.CurrentTotalCauchyStress[0];
  sigmaTot_yy+=rElementalVariables.CurrentTotalCauchyStress[1];
  sigmaTot_xy+=rElementalVariables.CurrentTotalCauchyStress[2];

  // // sigmaDev_xx= sigmaTot_xx - rElementalVariables.MeanPressure;
  // // sigmaDev_yy= sigmaTot_yy - rElementalVariables.MeanPressure;
  // // sigmaDev_xy= sigmaTot_xy;

  // sigmaTot_xx= sigmaDev_xx + rElementalVariables.MeanPressure;
  // sigmaTot_yy= sigmaDev_yy + rElementalVariables.MeanPressure;
  // sigmaTot_xy= sigmaDev_xy;

  rElementalVariables.UpdatedDeviatoricCauchyStress[0]=sigmaDev_xx;
  rElementalVariables.UpdatedDeviatoricCauchyStress[1]=sigmaDev_yy;
  rElementalVariables.UpdatedDeviatoricCauchyStress[2]=sigmaDev_xy;

  rElementalVariables.UpdatedTotalCauchyStress[0]=sigmaTot_xx;
  rElementalVariables.UpdatedTotalCauchyStress[1]=sigmaTot_yy;
  rElementalVariables.UpdatedTotalCauchyStress[2]=sigmaTot_xy;


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
   

  // this->mCurrentTotalCauchyStress[g]=rElementalVariables.CurrentTotalCauchyStress;
  this->mUpdatedTotalCauchyStress[g]=rElementalVariables.UpdatedTotalCauchyStress;
  // this->mCurrentDeviatoricCauchyStress[g]=rElementalVariables.CurrentDeviatoricCauchyStress;
  this->mUpdatedDeviatoricCauchyStress[g]=rElementalVariables.UpdatedDeviatoricCauchyStress;

}

template < > 
void UpdatedLagrangianVSolidElement<3>:: CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables, double TimeStep, unsigned int g,const ShapeFunctionsType& N)
{

  rElementalVariables.CurrentTotalCauchyStress=this->mCurrentTotalCauchyStress[g];
  rElementalVariables.CurrentDeviatoricCauchyStress=this->mCurrentDeviatoricCauchyStress[g];

  double Density  = 0;
  double CurrSecondLame  = 0;
  double CurrBulkModulus = 0;

  this->ComputeMaterialParameters(Density,CurrSecondLame,CurrBulkModulus,TimeStep,N);
  
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

  // sigmaTot_xx= sigmaDev_xx + rElementalVariables.MeanPressure;
  // sigmaTot_yy= sigmaDev_yy + rElementalVariables.MeanPressure;
  // sigmaTot_zz= sigmaDev_zz + rElementalVariables.MeanPressure;
  // sigmaTot_xy= sigmaDev_xy;
  // sigmaTot_xz= sigmaDev_xz;
  // sigmaTot_yz= sigmaDev_yz;

  sigmaTot_xx+=rElementalVariables.CurrentTotalCauchyStress[0];
  sigmaTot_yy+=rElementalVariables.CurrentTotalCauchyStress[1];
  sigmaTot_zz+=rElementalVariables.CurrentTotalCauchyStress[2];
  sigmaTot_xy+=rElementalVariables.CurrentTotalCauchyStress[3];
  sigmaTot_xz+=rElementalVariables.CurrentTotalCauchyStress[4];
  sigmaTot_yz+=rElementalVariables.CurrentTotalCauchyStress[5];

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




  template< unsigned int TDim >
  void UpdatedLagrangianVSolidElement<TDim>::CalculateLocalContinuityEqForPressure(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {


    GeometryType& rGeom = this->GetGeometry();
    // const SizeType NumNodes = rGeom.PointsNumber();
    const unsigned int NumNodes = rGeom.PointsNumber();
    // for (unsigned int n = 0; n < NumNodes; n++)
    //   {
    // 	if(rGeom[n].Is(BOUNDARY)){
    // 	  rGeom[n].Set(FLUID);	  
    // 	}
    //   }

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
      // const unsigned int NumGauss = GaussWeights.size();
      unsigned int NumGauss = GaussWeights.size();

      // Stabilization parameters

      // MatrixType BulkVelMatrix = ZeroMatrix(NumNodes,NumNodes);
      MatrixType BulkVelMatrixLump = ZeroMatrix(NumNodes,NumNodes);


      const double TimeStep=0.5/rCurrentProcessInfo[BDF_COEFFICIENTS][2];

      double Density  = 0;
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
	    const double GaussWeight = fabs(GaussWeights[g]);
	    const ShapeFunctionsType& N = row(NContainer,g);
	    // const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];
	    this->ComputeMaterialParameters(Density,DeviatoricCoeff,VolumetricCoeff,TimeStep,N);

	    double BulkCoeff =GaussWeight/(VolumetricCoeff);
	    if(rCurrentProcessInfo[STEP]>-1){
	      // this->ComputeBulkMatrixForPressureVel(BulkVelMatrix,N,BulkCoeff);
	      this->ComputeBulkMatrixForPressureVelLump(BulkVelMatrixLump,N,BulkCoeff);
	      rLeftHandSideMatrix+=BulkVelMatrixLump;	
	    }

	    VectorType UpdatedPressure;
	    VectorType CurrentPressure;

	    UpdatedPressure = ZeroVector(NumNodes);
	    CurrentPressure = ZeroVector(NumNodes);

	    this->GetPressureValues(UpdatedPressure,0);
	    this->GetPressureValues(CurrentPressure,1);

	    VectorType DeltaPressure = UpdatedPressure-CurrentPressure;

	    rRightHandSideVector -= prod(BulkVelMatrixLump,DeltaPressure);
	    double DivU=0;
	    DivU=rElementalVariables.VolumetricDefRate;

	    // Add convection, stabilization and RHS contributions to the local system equation
	    for (SizeType i = 0; i < NumNodes; ++i)
	      {
		// RHS contribution
		// Velocity divergence
		double RHSi =  N[i] * DivU;
		rRightHandSideVector[i] += GaussWeight * RHSi;
	
	      }

	  }

	}   
   

    }
  }
  

  template class UpdatedLagrangianVSolidElement<2>;
  template class UpdatedLagrangianVSolidElement<3>;

}
