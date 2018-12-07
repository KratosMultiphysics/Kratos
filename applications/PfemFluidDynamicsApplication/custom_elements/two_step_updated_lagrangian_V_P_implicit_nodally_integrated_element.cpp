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
      const ShapeFunctionDerivativesType& rDN_DX = DN_DX[0];
      
      // const unsigned int NumGauss = GaussWeights.size();
      // if(NumGauss==1){
      // 	this->CalcElementalStrains(rElementalVariables,rCurrentProcessInfo,rDN_DX);
      // }else{
      // 	std::cout<<"a different structure is required for more gauss points"<<std::endl;
      // }
 
      double elementVolume=0;
      if(TDim==3){
	elementVolume=rGeom.Volume()*0.25;
      }else if(TDim==2){
	elementVolume=rGeom.Area()/3.0;
      }

      for (unsigned int i = 0; i < NumNodes; i++)
	{

	  VectorType nodalSFDneighbours=rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS_ORDER);
	  unsigned int nodalSFDneighboursSize=nodalSFDneighbours.size();

	  if(nodalSFDneighboursSize>1){
	    double& meshSize = rGeom[i].FastGetSolutionStepValue(NODAL_MEAN_MESH_SIZE);
	    WeakPointerVector<Element >& neighb_elems = rGeom[i].GetValue(NEIGHBOUR_ELEMENTS);
	    double numberOfNeighElems=double(neighb_elems.size());
	    meshSize+=this->ElementSize()/numberOfNeighElems;

	    if(rGeom[i].Is(FREE_SURFACE)){
	      this->NodalFreeSurfaceLength(i);
	    }
	  
	    const double nodalVolume=rGeom[i].FastGetSolutionStepValue(NODAL_VOLUME);

	    // the nodal strain measure could be also be computed as follows (now they are computed in the strategy)
	    // rGeom[i].FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD)+=rElementalVariables.Fgrad*elementVolume/nodalVolume;
	    // rGeom[i].FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD_VEL)+=rElementalVariables.FgradVel*elementVolume/nodalVolume;
	    // rGeom[i].FastGetSolutionStepValue(NODAL_SPATIAL_DEF_RATE_BIS)+=rElementalVariables.SpatialDefRate*elementVolume/nodalVolume;

	      for (unsigned int j = 0; j< NumNodes; j++)
		{
		  unsigned int position=rGeom[j].Id();
	     
		  unsigned int SFDposition=0;
		  for (unsigned int k = 0; k< nodalSFDneighboursSize; k++)
		    {
		      if(position==nodalSFDneighbours[k]){
			rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition]   += rDN_DX(j,0)*elementVolume/nodalVolume;
			rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition+1] += rDN_DX(j,1)*elementVolume/nodalVolume;
			if(TDim==3){
			  rGeom[i].FastGetSolutionStepValue(NODAL_SFD_NEIGHBOURS)[SFDposition+2] += rDN_DX(j,2)*elementVolume/nodalVolume;
			}
			break;
		      }
		      SFDposition+=TDim;
		    }
	      
		}
	    
	  }
	  else{
	    std::cout<<rGeom[i].Id()<<"  this node is isolated!!! "<<std::endl;
	    for (unsigned int k = 0; k< TDim; k++)
	      {
		rGeom[i].FastGetSolutionStepValue(NODAL_DEFORMATION_GRAD)(k,k)=1.0;
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
    // unsigned int countFreeSurface=1;
    for (SizeType i = 0; i < NumNodes; i++){
      
      if((rGeom[i].Is(FREE_SURFACE)  || (rGeom[i].Is(SOLID) && rGeom[i].Is(BOUNDARY))) && i!=nodeIndex){
	noalias(Edge) = rGeom[nodeIndex].Coordinates() - rGeom[i].Coordinates();
	rGeom[nodeIndex].FastGetSolutionStepValue(NODAL_FREESURFACE_AREA) += sqrt(Edge[0]*Edge[0] + Edge[1]*Edge[1])/2.0;
	// countFreeSurface+=1;
      }
    }
    // if(countFreeSurface==NumNodes){
    //   WeakPointerVector<Element >& neighb_elemsA = rGeom[0].GetValue(NEIGHBOUR_ELEMENTS);
    //   WeakPointerVector<Element >& neighb_elemsB = rGeom[1].GetValue(NEIGHBOUR_ELEMENTS);
    //   WeakPointerVector<Element >& neighb_elemsC = rGeom[2].GetValue(NEIGHBOUR_ELEMENTS);
    //   if(neighb_elemsA.size()==1 && neighb_elemsB.size()==1 && neighb_elemsC.size()==1){
    // 	rGeom[nodeIndex].Set(ISOLATED);
    //   }
    // }
      
  }

  template< >
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<3>::NodalFreeSurfaceLength(unsigned int nodeIndex)
  {

    GeometryType& rGeom = this->GetGeometry();
    // const SizeType NumNodes = rGeom.PointsNumber();
    array_1d<double,2> Edge(3,0.0);


    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)){
      
      const double a = MathUtils<double>::Norm3(rGeom.GetPoint(0)-rGeom.GetPoint(1));
      const double b = MathUtils<double>::Norm3(rGeom.GetPoint(1)-rGeom.GetPoint(2));
      const double c = MathUtils<double>::Norm3(rGeom.GetPoint(2)-rGeom.GetPoint(0));
	  
      const double s = (a+b+c) / 2.0;
	  
      rGeom[ nodeIndex].FastGetSolutionStepValue(NODAL_FREESURFACE_AREA) += std::sqrt(s*(s-a)*(s-b)*(s-c))/3.0;
      
    }
    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
            
      const double a = MathUtils<double>::Norm3(rGeom.GetPoint(0)-rGeom.GetPoint(1));
      const double b = MathUtils<double>::Norm3(rGeom.GetPoint(1)-rGeom.GetPoint(3));
      const double c = MathUtils<double>::Norm3(rGeom.GetPoint(3)-rGeom.GetPoint(0));
	  
      const double s = (a+b+c) / 2.0;
	  
      rGeom[ nodeIndex].FastGetSolutionStepValue(NODAL_FREESURFACE_AREA) += std::sqrt(s*(s-a)*(s-b)*(s-c))/3.0;

    }
    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
      
      const double a = MathUtils<double>::Norm3(rGeom.GetPoint(0)-rGeom.GetPoint(2));
      const double b = MathUtils<double>::Norm3(rGeom.GetPoint(2)-rGeom.GetPoint(3));
      const double c = MathUtils<double>::Norm3(rGeom.GetPoint(3)-rGeom.GetPoint(0));
	  
      const double s = (a+b+c) / 2.0;
	  
      rGeom[ nodeIndex].FastGetSolutionStepValue(NODAL_FREESURFACE_AREA) += std::sqrt(s*(s-a)*(s-b)*(s-c))/3.0;
      
    }
    if(rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){      
      
      const double a = MathUtils<double>::Norm3(rGeom.GetPoint(1)-rGeom.GetPoint(2));
      const double b = MathUtils<double>::Norm3(rGeom.GetPoint(2)-rGeom.GetPoint(3));
      const double c = MathUtils<double>::Norm3(rGeom.GetPoint(3)-rGeom.GetPoint(1));
	  
      const double s = (a+b+c) / 2.0;
	  
      rGeom[ nodeIndex].FastGetSolutionStepValue(NODAL_FREESURFACE_AREA) += std::sqrt(s*(s-a)*(s-b)*(s-c))/3.0;
      
    }

  }


  template<>
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<2>::GetNodesPosition(Vector& rValues,const ProcessInfo& rCurrentProcessInfo, double theta)
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
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<3>::GetNodesPosition(Vector& rValues,const ProcessInfo& rCurrentProcessInfo, double theta)
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

  template < > 
  void TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<2>::CalcElementalStrains(ElementalVariables & rElementalVariables,
											 const ProcessInfo &rCurrentProcessInfo,
											 const ShapeFunctionDerivativesType& rDN_DX)
  {
    unsigned int dimension=this->GetGeometry().WorkingSpaceDimension();
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = dimension*NumNodes;
    VectorType  NodesPosition= ZeroVector(LocalSize);
    VectorType VelocityValues = ZeroVector(LocalSize);
    VectorType RHSVelocities = ZeroVector(LocalSize);
    double theta=0.5;
    // if(rGeom[0].Is(SOLID) && rGeom[1].Is(SOLID) && rGeom[2].Is(SOLID)){
    //   theta=1.0;
    // }
    this->GetNodesPosition(NodesPosition,rCurrentProcessInfo,theta);
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
		rElementalVariables.Fgrad(i,j)+=NodesPosition[dimension*k+i]*rDN_DX(k,j);
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
    unsigned int dimension=this->GetGeometry().WorkingSpaceDimension();
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();
    const SizeType LocalSize = dimension*NumNodes;
    VectorType  NodesPosition= ZeroVector(LocalSize);
    VectorType VelocityValues = ZeroVector(LocalSize);
    VectorType RHSVelocities = ZeroVector(LocalSize);
    double theta=0.5;
    if(rGeom[0].Is(SOLID) && rGeom[1].Is(SOLID) && rGeom[2].Is(SOLID) && rGeom[3].Is(SOLID)){
      theta=1.0;
    }
    this->GetNodesPosition(NodesPosition,rCurrentProcessInfo,theta);
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
		rElementalVariables.Fgrad(i,j)+=NodesPosition[dimension*k+i]*rDN_DX(k,j);
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

  template class TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<2>;
  template class TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<3>;

}
