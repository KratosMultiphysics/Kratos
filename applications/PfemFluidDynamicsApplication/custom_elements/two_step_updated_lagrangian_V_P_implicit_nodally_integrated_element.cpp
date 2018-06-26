//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:               June 2018 $
//   Revision:            $Revision:                 0.0 $
//
//   Implementation of the Gauss-Seidel two step Updated Lagrangian Velocity-Pressure element
//     ( There is a ScalingConstant to multiply the mass balance equation for a number because i read it somewhere)
//

// System includes

// External includes
 
// Project includes
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_nodally_integrated_element.h"
#include "includes/cfd_variables.h"

namespace Kratos {


  template< unsigned int TDim >
  Element::Pointer TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {
    KRATOS_TRY;

    TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement NewElement(NewId, this->GetGeometry().Create( rThisNodes ), this->pGetProperties() );
    return Element::Pointer( new TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement(NewElement) );
    
    KRATOS_CATCH("");

  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
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
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
											    ProcessInfo& rCurrentProcessInfo)
  { 
    KRATOS_TRY;


    
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = TDim * NumNodes;

    // Check sizes and initialize
    if( rLeftHandSideMatrix.size1() != LocalSize )
      rLeftHandSideMatrix.resize(LocalSize,LocalSize);

    rLeftHandSideMatrix = ZeroMatrix(LocalSize,LocalSize);

    // Shape functions and integration points
    ShapeFunctionDerivativesArrayType DN_DX;
    Matrix NContainer;
    VectorType GaussWeights;
    this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);
    
    const ShapeFunctionDerivativesType& rDN_DX = DN_DX[0]; 

    SizeType FirstRow=0;
    SizeType FirstCol=0;


    /// to check----é sicuramente sbagliato!! non ci fare cas0!!!!!!!!!!
    for (SizeType j = 0; j < NumNodes; ++j)
      {
	const double nodalPatchVolume=rGeom[j].FastGetSolutionStepValue(NODAL_AREA);
	for (SizeType i = 0; i < NumNodes; ++i)
	  {
          // First Row
            rLeftHandSideMatrix(FirstRow,FirstCol) += rDN_DX(j,0)*GaussWeights[0]/nodalPatchVolume;
            rLeftHandSideMatrix(FirstRow,FirstCol+1) += rDN_DX(j,1)*GaussWeights[0]/nodalPatchVolume;

            // Second Row
            rLeftHandSideMatrix(FirstRow+1,FirstCol) += rDN_DX(j,0)*GaussWeights[0]/nodalPatchVolume;
            rLeftHandSideMatrix(FirstRow+1,FirstCol+1) += rDN_DX(j,1)*GaussWeights[0]/nodalPatchVolume;

            // Update Counter
            FirstRow += 2;
	  }
        FirstRow = 0;
        FirstCol += 2;
      }


    
    KRATOS_CATCH("");
  }

    template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
    {
      KRATOS_TRY;

      // std::cout<<"InitializeNonLinearIteration "<<std::endl;
      GeometryType& rGeom = this->GetGeometry();
      const unsigned int NumNodes = rGeom.PointsNumber();
      // Shape functions and integration points
      ShapeFunctionDerivativesArrayType DN_DX;
      Matrix NContainer;
      VectorType GaussWeights;
      this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);
      ElementalVariables rElementalVariables;
      this->InitializeElementalVariables(rElementalVariables);
      // Loop on integration points
      const unsigned int NumGauss = GaussWeights.size();
      const ShapeFunctionDerivativesType& rDN_DX = DN_DX[0]; 
      if(NumGauss==1){
	this->CalcElementalStrains(rElementalVariables,rCurrentProcessInfo,rDN_DX);
      }else{
	std::cout<<"a different structure is required for more gauss points"<<std::endl;
      }
 
      double elementVolume=0;
      if(TDim==3){
	elementVolume=rGeom.Volume()*0.25;
      }else if(TDim==2){
	elementVolume=rGeom.Area()/3.0;
      }

      for (unsigned int i = 0; i < NumNodes; i++)
	{
	  double& meshSize = rGeom[i].FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE);
	  WeakPointerVector<Element >& neighb_elems = rGeom[i].GetValue(NEIGHBOUR_ELEMENTS);
	  double numberOfNeighElems=double(neighb_elems.size());
	  meshSize+=this->ElementSize()/numberOfNeighElems;

	  if(rGeom[i].Is(FREE_SURFACE) || (rGeom[i].Is(SOLID) && rGeom[i].Is(BOUNDARY))){
	    this->NodalFreeSurfaceLength(i);
	  }

	  // std::cout<<"meshSize "<<meshSize<<"  numberOfNeighElems="<<numberOfNeighElems<<std::endl;
	  
	  const double nodalVolume=rGeom[i].FastGetSolutionStepValue(NODAL_AREA);
	  rGeom[i].FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD)+=rElementalVariables.Fgrad*elementVolume/nodalVolume;
      	  rGeom[i].FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD_VEL)+=rElementalVariables.FgradVel*elementVolume/nodalVolume;
      	  rGeom[i].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE_BIS)+=rElementalVariables.SpatialDefRate*elementVolume/nodalVolume;
	  
	  VectorType nodalSFDneighbours=rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER);
	  unsigned int nodalSFDneighboursSize=nodalSFDneighbours.size();
	   
	  for (unsigned int j = 0; j< NumNodes; j++)
	    {
	      unsigned int position=rGeom[j].Id();
	      double dnDX=rDN_DX(j,0)*elementVolume/nodalVolume;
	      double dnDY=rDN_DX(j,1)*elementVolume/nodalVolume;
	      double dnDZ=rDN_DX(j,2)*elementVolume/nodalVolume;

	      unsigned int SFDposition=0;
	      for (unsigned int k = 0; k< nodalSFDneighboursSize; k++)
		{
		  if(position==nodalSFDneighbours[k]){
		    // rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_COMMON_ELEMENTS)[k]+=1.0;
		    rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_VOLUME_COMMON_ELEMENTS)[k]=elementVolume*3.0;
		    if(TDim==3){
		      rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition]+=dnDX;
		      rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition+1]+=dnDY;
		      rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition+2]+=dnDZ;		      
		      break;
		    }else if(TDim==2){
		      rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition]+=dnDX;
		      rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition+1]+=dnDY;
		      break;
		    }
		  }
		  if(TDim==3){
		    SFDposition+=3;
		  }else if(TDim==2){
		    SFDposition+=2;
		  }
		}
	      

	    }
	}
	 


      
      KRATOS_CATCH( "" );

    }



  
  template< >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<2>::NodalFreeSurfaceLength(unsigned int nodeIndex)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    array_1d<double,2> Edge(2,0.0);

    for (SizeType i = 0; i < NumNodes; i++){
      
      if((rGeom[i].Is(FREE_SURFACE)  || (rGeom[i].Is(SOLID) && rGeom[i].Is(BOUNDARY))) && i!=nodeIndex){
	noalias(Edge) = rGeom[nodeIndex].Coordinates() - rGeom[i].Coordinates();
	rGeom[i].FastGetSolutionStepValue(NODAL_FREESURFACE_AREA) += sqrt(Edge[0]*Edge[0] + Edge[1]*Edge[1])/2.0;
      }
    }
      
  }

  template< >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<3>::NodalFreeSurfaceLength(unsigned int nodeIndex)
  {
    std::cout<<"NodalFreeSurfaceLength is not yet implemented for the 3D"<<std::endl;      
  }

  
    template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::CalculateLocalMomentumEquations(MatrixType& rLeftHandSideMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY; 

    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = TDim * NumNodes;

    MatrixType MassMatrix= ZeroMatrix(LocalSize,LocalSize);
    MatrixType StiffnessMatrix= ZeroMatrix(LocalSize,LocalSize);

    // Check sizes and initialize
    if( rLeftHandSideMatrix.size1() != LocalSize )
      rLeftHandSideMatrix.resize(LocalSize,LocalSize);

    rLeftHandSideMatrix = ZeroMatrix(LocalSize,LocalSize);

    if( rRightHandSideVector.size() != LocalSize )
      rRightHandSideVector.resize(LocalSize);

    rRightHandSideVector = ZeroVector(LocalSize);

    // Shape functions and integration points
    ShapeFunctionDerivativesArrayType DN_DX;
    Matrix NContainer;
    VectorType GaussWeights;
    this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);

    const double TimeStep=rCurrentProcessInfo[DELTA_TIME];

    double totalVolume=0;
    double MeanValueMass=0;
    double Density=this->GetProperties()[DENSITY];

      if(TDim==3){
	totalVolume=rGeom.Volume();
      }else if(TDim==2){
	totalVolume=rGeom.Area();
      }

      this->ComputeExternalForces(rRightHandSideVector,Density,totalVolume);


    double lumpedDynamicWeight=totalVolume*Density;
    this->ComputeLumpedMassMatrix(MassMatrix,lumpedDynamicWeight,MeanValueMass);    

    // double BulkReductionCoefficient=1.0;
    // double MeanValueStiffness=0.0;
    // this->ComputeBulkReductionCoefficient(MassMatrix,StiffnessMatrix,MeanValueStiffness,BulkReductionCoefficient,TimeStep);
    // if(BulkReductionCoefficient!=1.0){
    //   // VolumetricCoeff*=BulkReductionCoefficient;
    //   VolumetricCoeff*=MeanValueMass*2.0/(TimeStep*MeanValueStiffness);
    //   StiffnessMatrix= ZeroMatrix(LocalSize,LocalSize);

    //   for (unsigned int g = 0; g < NumGauss; g++)
    // 	{
    // 	  const double GaussWeight = GaussWeights[g];
    // 	  const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g]; 
    // 	  this->ComputeCompleteTangentTerm(rElementalVariables,StiffnessMatrix,rDN_DX,DeviatoricCoeff,VolumetricCoeff,theta,nodalPatchVolume);
    // 	}
    // }

    // Add residual of previous iteration to RHS
    VectorType VelocityValues = ZeroVector(LocalSize);
    VectorType AccelerationValues = ZeroVector(LocalSize);

    //2nd order 
    this->GetAccelerationValues(AccelerationValues,0);
    this->GetVelocityValues(VelocityValues,0);
    noalias(AccelerationValues)+=-2.0*VelocityValues/TimeStep;
    this->GetVelocityValues(VelocityValues,1);
    noalias(AccelerationValues)+=2.0*VelocityValues/TimeStep;//these are negative accelerations

    // MassMatrix= ZeroMatrix(LocalSize,LocalSize);
    noalias( rRightHandSideVector )+= prod(MassMatrix,AccelerationValues);
    noalias( rLeftHandSideMatrix ) +=  StiffnessMatrix + MassMatrix*2/TimeStep;
  


    KRATOS_CATCH( "" );
 
  }



  template < > 
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<2>::CalcElementalStrains(ElementalVariables & rElementalVariables,
											 const ProcessInfo &rCurrentProcessInfo,
											 const ShapeFunctionDerivativesType& rDN_DX)
  {
    const double theta=this->GetThetaMomentum();
    unsigned int dimension=this->GetGeometry().WorkingSpaceDimension();
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = dimension*NumNodes;
    VectorType  NodePosition= ZeroVector(LocalSize);
    VectorType VelocityValues = ZeroVector(LocalSize);
    VectorType RHSVelocities = ZeroVector(LocalSize);
    this->GetPositions(NodePosition,rCurrentProcessInfo,theta);
    this->GetVelocityValues(RHSVelocities,0); 
    RHSVelocities*=theta;
    this->GetVelocityValues(VelocityValues,1);
    RHSVelocities+=VelocityValues*(1.0-theta);
    rElementalVariables.Fgrad=ZeroMatrix(dimension,dimension);
    rElementalVariables.FgradVel=ZeroMatrix(dimension,dimension);
    for (SizeType i = 0; i < dimension; i++)
      {
	for (SizeType j = 0; j < dimension; j++)
	  {
	    for (SizeType k = 0; k < NumNodes; k++)
	      {
		rElementalVariables.Fgrad(i,j)+= NodePosition[dimension*k+i]*rDN_DX(k,j);
		rElementalVariables.FgradVel(i,j)+= RHSVelocities[dimension*k+i]*rDN_DX(k,j);
	      }

	  }
      }

    //Inverse
    rElementalVariables.InvFgrad=ZeroMatrix(dimension,dimension);
    rElementalVariables.DetFgrad=1;
    MathUtils<double>::InvertMatrix(rElementalVariables.Fgrad, 
				    rElementalVariables.InvFgrad, 
				    rElementalVariables.DetFgrad);

    // rElementalVariables.InvFgradVel=ZeroMatrix(dimension,dimension);
    // rElementalVariables.DetFgradVel=1;
    // MathUtils<double>::InvertMatrix(rElementalVariables.FgradVel, 
    // 				    rElementalVariables.InvFgradVel, 
    // 				    rElementalVariables.DetFgradVel);

    //it computes the spatial velocity gradient tensor --> [L_ij]=dF_ik*invF_kj
    rElementalVariables.SpatialVelocityGrad.resize(dimension,dimension);
    rElementalVariables.SpatialVelocityGrad=prod(rElementalVariables.FgradVel,rElementalVariables.InvFgrad);

    rElementalVariables.VolumetricDefRate=0;
    for (SizeType i = 0; i < dimension; i++)
      {
	rElementalVariables.VolumetricDefRate+=rElementalVariables.SpatialVelocityGrad(i,i);
      }

    rElementalVariables.SpatialDefRate[0]=rElementalVariables.SpatialVelocityGrad(0,0);
    rElementalVariables.SpatialDefRate[1]=rElementalVariables.SpatialVelocityGrad(1,1);
    rElementalVariables.SpatialDefRate[2]=0.5*(rElementalVariables.SpatialVelocityGrad(1,0)+rElementalVariables.SpatialVelocityGrad(0,1));
  
    double aThird=1.0/3.0;
    double dev_X=rElementalVariables.SpatialDefRate[0]-
      (rElementalVariables.SpatialDefRate[0]+rElementalVariables.SpatialDefRate[1])*aThird;
    double dev_Y=rElementalVariables.SpatialDefRate[1]-
      (rElementalVariables.SpatialDefRate[0]+rElementalVariables.SpatialDefRate[1])*aThird;
    rElementalVariables.DeviatoricInvariant=sqrt(2*(dev_X*dev_X + dev_Y*dev_Y +
						    rElementalVariables.SpatialDefRate[2]*rElementalVariables.SpatialDefRate[2]));

    rElementalVariables.EquivalentStrainRate=sqrt((2.0*rElementalVariables.SpatialDefRate[0]*rElementalVariables.SpatialDefRate[0] +
						   2.0*rElementalVariables.SpatialDefRate[1]*rElementalVariables.SpatialDefRate[1] +
						   4.0*rElementalVariables.SpatialDefRate[2]*rElementalVariables.SpatialDefRate[2]));


  }  

  template < > 
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<3>::CalcElementalStrains(ElementalVariables & rElementalVariables,
											 const ProcessInfo &rCurrentProcessInfo,
											 const ShapeFunctionDerivativesType& rDN_DX)
  {
    const double theta=this->GetThetaMomentum();
    unsigned int dimension=this->GetGeometry().WorkingSpaceDimension();
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = dimension*NumNodes;
    VectorType  NodePosition= ZeroVector(LocalSize);
    VectorType VelocityValues = ZeroVector(LocalSize);
    VectorType RHSVelocities = ZeroVector(LocalSize);
    this->GetPositions(NodePosition,rCurrentProcessInfo,theta);
    this->GetVelocityValues(RHSVelocities,0); 
    RHSVelocities*=theta;
    this->GetVelocityValues(VelocityValues,1);
    RHSVelocities+=VelocityValues*(1.0-theta);

    rElementalVariables.Fgrad=ZeroMatrix(dimension,dimension);
    rElementalVariables.FgradVel=ZeroMatrix(dimension,dimension);
    for (SizeType i = 0; i < dimension; i++)
      {
	for (SizeType j = 0; j < dimension; j++)
	  {
	    for (SizeType k = 0; k < NumNodes; k++)
	      {
		rElementalVariables.Fgrad(i,j)+= NodePosition[dimension*k+i]*rDN_DX(k,j);
		rElementalVariables.FgradVel(i,j)+= RHSVelocities[dimension*k+i]*rDN_DX(k,j);
	      }
	  }
      }

    //Inverse
    rElementalVariables.InvFgrad=ZeroMatrix(dimension,dimension);
    rElementalVariables.DetFgrad=1;
    MathUtils<double>::InvertMatrix(rElementalVariables.Fgrad, 
				    rElementalVariables.InvFgrad, 
				    rElementalVariables.DetFgrad);

    // rElementalVariables.InvFgradVel=ZeroMatrix(dimension,dimension);
    // rElementalVariables.DetFgradVel=1;
    // MathUtils<double>::InvertMatrix(rElementalVariables.FgradVel, 
    // 				    rElementalVariables.InvFgradVel, 
    // 				    rElementalVariables.DetFgradVel);


    //it computes the spatial velocity gradient tensor --> [L_ij]=dF_ik*invF_kj
    rElementalVariables.SpatialVelocityGrad.resize(dimension,dimension);
    rElementalVariables.SpatialVelocityGrad=prod(rElementalVariables.FgradVel,rElementalVariables.InvFgrad);
  
    rElementalVariables.VolumetricDefRate=0;
    for (SizeType i = 0; i < dimension; i++)
      {
	rElementalVariables.VolumetricDefRate+=rElementalVariables.SpatialVelocityGrad(i,i);
      }
  
    rElementalVariables.SpatialDefRate[0]=rElementalVariables.SpatialVelocityGrad(0,0);
    rElementalVariables.SpatialDefRate[1]=rElementalVariables.SpatialVelocityGrad(1,1);
    rElementalVariables.SpatialDefRate[2]=rElementalVariables.SpatialVelocityGrad(2,2);
    rElementalVariables.SpatialDefRate[3]=0.5*(rElementalVariables.SpatialVelocityGrad(1,0)+rElementalVariables.SpatialVelocityGrad(0,1));
    rElementalVariables.SpatialDefRate[4]=0.5*(rElementalVariables.SpatialVelocityGrad(2,0)+rElementalVariables.SpatialVelocityGrad(0,2));
    rElementalVariables.SpatialDefRate[5]=0.5*(rElementalVariables.SpatialVelocityGrad(2,1)+rElementalVariables.SpatialVelocityGrad(1,2));
    // computeElement=CheckStrain3(rElementalVariables.SpatialDefRate,rElementalVariables.SpatialVelocityGrad);

    double aThird=1.0/3.0;
    double dev_X=rElementalVariables.SpatialDefRate[0]-
      (rElementalVariables.SpatialDefRate[0]+rElementalVariables.SpatialDefRate[1]+rElementalVariables.SpatialDefRate[2])*aThird;
    double dev_Y=rElementalVariables.SpatialDefRate[1]-
      (rElementalVariables.SpatialDefRate[0]+rElementalVariables.SpatialDefRate[1]+rElementalVariables.SpatialDefRate[2])*aThird;
    double dev_Z=rElementalVariables.SpatialDefRate[2]-
      (rElementalVariables.SpatialDefRate[0]+rElementalVariables.SpatialDefRate[1]+rElementalVariables.SpatialDefRate[2])*aThird;
    rElementalVariables.DeviatoricInvariant=sqrt(2*(dev_X*dev_X + dev_Y*dev_Y + dev_Z*dev_Z +
						    rElementalVariables.SpatialDefRate[3]*rElementalVariables.SpatialDefRate[3] +
						    rElementalVariables.SpatialDefRate[4]*rElementalVariables.SpatialDefRate[4] +
						    rElementalVariables.SpatialDefRate[5]*rElementalVariables.SpatialDefRate[5]));

    rElementalVariables.EquivalentStrainRate=sqrt(2.0*(rElementalVariables.SpatialDefRate[0]*rElementalVariables.SpatialDefRate[0] +
						       rElementalVariables.SpatialDefRate[1]*rElementalVariables.SpatialDefRate[1] +
						       rElementalVariables.SpatialDefRate[2]*rElementalVariables.SpatialDefRate[2] +
						       2.0*rElementalVariables.SpatialDefRate[3]*rElementalVariables.SpatialDefRate[3] +
						       2.0*rElementalVariables.SpatialDefRate[4]*rElementalVariables.SpatialDefRate[4] +
						       2.0*rElementalVariables.SpatialDefRate[5]*rElementalVariables.SpatialDefRate[5]));


  }
  

  
  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::CalculateGeometryData(ShapeFunctionDerivativesArrayType &rDN_DX,
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
      // if(rGaussWeights[g]<0)
      // 	std::cout<<"NEGATIVE GAUSS WEIGHT "<<rGaussWeights[g]<<std::endl;
    }
  }

  
  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
												    std::vector<double>& rValues,
												    const ProcessInfo& rCurrentProcessInfo )
  {
    if ( rVariable == YIELDED)
      {
	rValues[0]=this->GetValue(YIELDED);
      }
    if ( rVariable == FLOW_INDEX)
      {
	rValues[0]=this->GetValue(FLOW_INDEX);
      }	
  }


  template<>
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<2>::ComputeCompleteTangentTerm(ElementalVariables & rElementalVariables,
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

    SizeType FirstRow=0;
    SizeType FirstCol=0;

    for (SizeType j = 0; j < NumNodes; ++j)
      {

        for (SizeType i = 0; i < NumNodes; ++i)
	  {
	    // double lagDNXi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,0);
	    // double lagDNYi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,1);
	    // double lagDNXj=rDN_DX(j,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(j,1)*rElementalVariables.InvFgrad(1,0);
	    // double lagDNYj=rDN_DX(j,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(j,1)*rElementalVariables.InvFgrad(1,1);
	    double coeffI=this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA)/Weight;
	    double coeffJ=this->GetGeometry()[j].FastGetSolutionStepValue(NODAL_AREA)/Weight;
	    coeffI=1.0;
	    coeffJ=1.0;
	    double lagDNXi=rDN_DX(i,0)/coeffI;
	    double lagDNYi=rDN_DX(i,1)/coeffI;
	    double lagDNXj=rDN_DX(j,0)/coeffJ;
	    double lagDNYj=rDN_DX(j,1)/coeffJ;


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
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<3>::ComputeCompleteTangentTerm(ElementalVariables & rElementalVariables,
											       MatrixType& rDampingMatrix,
											       const ShapeFunctionDerivativesType& rDN_DX,
											       const double secondLame,
											       const double bulkModulus,
											       const double theta,
											       const double Weight){
    
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    SizeType FirstRow=0;
    SizeType FirstCol=0;

    for (SizeType j = 0; j < NumNodes; ++j)
      {
        for (SizeType i = 0; i < NumNodes; ++i)
	  {
	    // double lagDNXi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,0)+rDN_DX(i,2)*rElementalVariables.InvFgrad(2,0);
	    // double lagDNYi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,1)+rDN_DX(i,2)*rElementalVariables.InvFgrad(2,1);
	    // double lagDNZi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,2)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,2)+rDN_DX(i,2)*rElementalVariables.InvFgrad(2,2);
	    // double lagDNXj=rDN_DX(j,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(j,1)*rElementalVariables.InvFgrad(1,0)+rDN_DX(j,2)*rElementalVariables.InvFgrad(2,0);
	    // double lagDNYj=rDN_DX(j,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(j,1)*rElementalVariables.InvFgrad(1,1)+rDN_DX(j,2)*rElementalVariables.InvFgrad(2,1);
	    // double lagDNZj=rDN_DX(j,0)*rElementalVariables.InvFgrad(0,2)+rDN_DX(j,1)*rElementalVariables.InvFgrad(1,2)+rDN_DX(j,2)*rElementalVariables.InvFgrad(2,2);
	    
	    double lagDNXi=rDN_DX(i,0);
	    double lagDNYi=rDN_DX(i,1);
	    double lagDNZi=rDN_DX(i,2);
	    double lagDNXj=rDN_DX(j,0);
	    double lagDNYj=rDN_DX(j,1);
	    double lagDNZj=rDN_DX(j,2);

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
  int TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo)
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

  template< >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<2>::ComputeInternalForces(Vector& rRHSVector,
											  const ShapeFunctionDerivativesType& rDN_DX,
											  ElementalVariables& rElementalVariables,
											  const double Weight)
  {

    std::cout<<"DO NOT ENTER HERE ComputeInternalForces"<<std::endl;
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    SizeType FirstRow = 0;    

    for (SizeType i = 0; i < NumNodes; ++i)
      {
	// double lagDNXi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,0);
	// double lagDNYi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,1);
	double lagDNXi=rDN_DX(i,0);
	double lagDNYi=rDN_DX(i,1);

	rRHSVector[FirstRow]   += -Weight*(lagDNXi*rGeom[i].GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[0] +
					   lagDNYi*rGeom[i].GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[2]);
	
	rRHSVector[FirstRow+1] += -Weight*(lagDNYi*rGeom[i].GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[1] +
					   lagDNXi*rGeom[i].GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[2]);

	// rRHSVector[FirstRow]   += -Weight*(lagDNXi*(rGeom[i].GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,0)[0]+rGeom[i].GetSolutionStepValue(PRESSURE)) +
	// 				   lagDNYi*rGeom[i].GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[2]);
	
	// rRHSVector[FirstRow+1] += -Weight*(lagDNYi*(rGeom[i].GetSolutionStepValue(NODAL_DEVIATORIC_CAUCHY_STRESS,0)[1]+rGeom[i].GetSolutionStepValue(PRESSURE)) +
	// 				   lagDNXi*rGeom[i].GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[2]);

	// rRHSVector[FirstRow]   += -Weight*(lagDNXi*sigmaX + lagDNYi*sigmaXY);
	
	// rRHSVector[FirstRow+1] += -Weight*(lagDNYi*sigmaY + lagDNXi*sigmaXY);
	

	
	// rRHSVector[FirstRow]   += -Weight*(lagDNXi*rElementalVariables.UpdatedTotalCauchyStress[0] +
	// 				   lagDNYi*rElementalVariables.UpdatedTotalCauchyStress[2]);
	
	// rRHSVector[FirstRow+1] += -Weight*(lagDNYi*rElementalVariables.UpdatedTotalCauchyStress[1] +
	// 				   lagDNXi*rElementalVariables.UpdatedTotalCauchyStress[2]);
	
	FirstRow += 2;
      }
  }

  template< >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<3>::ComputeInternalForces(Vector& rRHSVector,
											  const ShapeFunctionDerivativesType& rDN_DX,
											  ElementalVariables& rElementalVariables,
											  const double Weight)
  {
    std::cout<<"DO NOT ENTER HERE ComputeInternalForces"<<std::endl;

    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    SizeType FirstRow = 0;

    for (SizeType i = 0; i < NumNodes; ++i)
      {
	// double lagDNXi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,0)+rDN_DX(i,2)*rElementalVariables.InvFgrad(2,0);
	// double lagDNYi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,1)+rDN_DX(i,2)*rElementalVariables.InvFgrad(2,1);
	// double lagDNZi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,2)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,2)+rDN_DX(i,2)*rElementalVariables.InvFgrad(2,2);
	double lagDNXi=rDN_DX(i,0);
	double lagDNYi=rDN_DX(i,1);
	double lagDNZi=rDN_DX(i,2);

	rRHSVector[FirstRow]   += -Weight*(lagDNXi*rGeom[i].GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[0] +
					   lagDNYi*rGeom[i].GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[3] +
					   lagDNZi*rGeom[i].GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[4]);

	rRHSVector[FirstRow+1] += -Weight*(lagDNYi*rGeom[i].GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[1] +
					   lagDNXi*rGeom[i].GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[3] +
					   lagDNZi*rGeom[i].GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[5]);

	rRHSVector[FirstRow+2] += -Weight*(lagDNZi*rGeom[i].GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[2] +
					   lagDNXi*rGeom[i].GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[4] +
					   lagDNYi*rGeom[i].GetSolutionStepValue(NODAL_CAUCHY_STRESS,0)[5]);

	FirstRow += 3;
      }


  }


  

  
  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>::ComputeExternalForces(Vector& rRHSVector,
											 const double Density,
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
	      double coeff=0;
	      if(TDim==2){
		coeff=1.0/3.0;
	      }else if(TDim==3){
		coeff=0.25;
	      }
	      rRHSVector[FirstRow+d] += Weight * Density * coeff * VolumeAcceleration[d];
	      // std::cout<<"rRHSVector[FirstRow+d] "<<rRHSVector[FirstRow+d]<<std::endl;
	    }

	}

        FirstRow += TDim;

      }
  }




  /*
   * Template class definition (this should allow us to compile the desired template instantiations)
   */

  template class TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<2>;
  template class TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<3>;

}
