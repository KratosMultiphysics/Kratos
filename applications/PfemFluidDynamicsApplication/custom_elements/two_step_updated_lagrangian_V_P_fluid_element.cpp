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
 
    TwoStepUpdatedLagrangianVPFluidElement NewElement(NewId, this->GetGeometry().Create( rThisNodes ), this->pGetProperties() );

    NewElement.SetData(this->GetData());
    NewElement.SetFlags(this->GetFlags());

    return Element::Pointer( new TwoStepUpdatedLagrangianVPFluidElement(NewElement) );

  }




  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPFluidElement<TDim>::Initialize()
  {
    KRATOS_TRY; 
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
									       double timeStep)
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
  void TwoStepUpdatedLagrangianVPFluidElement<2>::ComputeMeanValueMaterialTangentMatrix(ElementalVariables & rElementalVariables,double& MeanValue,const ShapeFunctionDerivativesType& rDN_DX,const double secondLame,double & bulkModulus,const double Weight,double& MeanValueMass,const double TimeStep)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    const double FourThirds = 4.0 / 3.0;
    const double nTwoThirds = -2.0 / 3.0;

    double theta=0.5;
    double Count=0;
    for (SizeType j = 0; j < NumNodes; ++j)
      {
        for (SizeType i = 0; i < NumNodes; ++i)
	  {
	    double lagDNXi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,0);
	    double lagDNYi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,1);
	    double lagDNXj=rDN_DX(j,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(j,1)*rElementalVariables.InvFgrad(1,0);
	    double lagDNYj=rDN_DX(j,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(j,1)*rElementalVariables.InvFgrad(1,1);
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

    double theta=0.5;
    double Count=0;
    for (SizeType j = 0; j < NumNodes; ++j)
      {
        for (SizeType i = 0; i < NumNodes; ++i)
	  {
	    double lagDNXi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,0)+rDN_DX(i,2)*rElementalVariables.InvFgrad(2,0);
	    double lagDNYi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,1)+rDN_DX(i,2)*rElementalVariables.InvFgrad(2,1);
	    double lagDNZi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,2)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,2)+rDN_DX(i,2)*rElementalVariables.InvFgrad(2,2);
	    double lagDNXj=rDN_DX(j,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(j,1)*rElementalVariables.InvFgrad(1,0)+rDN_DX(j,2)*rElementalVariables.InvFgrad(2,0);
	    double lagDNYj=rDN_DX(j,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(j,1)*rElementalVariables.InvFgrad(1,1)+rDN_DX(j,2)*rElementalVariables.InvFgrad(2,1);
	    double lagDNZj=rDN_DX(j,0)*rElementalVariables.InvFgrad(0,2)+rDN_DX(j,1)*rElementalVariables.InvFgrad(1,2)+rDN_DX(j,2)*rElementalVariables.InvFgrad(2,2);	  
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


  template<>
  void TwoStepUpdatedLagrangianVPFluidElement<2>::ComputeBulkReductionCoefficient(MatrixType MassMatrix,
										  MatrixType StiffnessMatrix,
										  double& meanValueStiff,
										  double& bulkCoefficient,
										  double timeStep)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    IndexType FirstRow = 0;
    IndexType FirstCol = 0;
    double meanValueMass=0;
    double countStiff=0;
    double countMass=0;
    for (SizeType j = 0; j < NumNodes; ++j)
      {
        for (SizeType i = 0; i < NumNodes; ++i)
	  {
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow,FirstCol));
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow,FirstCol+1));
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow+1,FirstCol));
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow+1,FirstCol+1));
            // Update Counter
	    countStiff+=4.0;

	    meanValueMass+= fabs(MassMatrix(FirstRow,FirstCol));
	    meanValueMass+= fabs(MassMatrix(FirstRow+1,FirstCol+1));
            // Update Counter
	    countMass+=2.0;

            FirstRow += 2;

	  }
	FirstRow = 0;
        FirstCol += 2;
      }

    meanValueStiff*=1.0/countStiff;
    meanValueMass*=1.0/countMass;
    
    if(meanValueMass!=0 && meanValueStiff!=0){
      bulkCoefficient=meanValueMass*4/(timeStep*meanValueStiff);
    }else{
      std::cout<<" DANGEROUS ELEMENT!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
      std::cout<<" coordinates "<<this->GetGeometry()[0].X()<<" "<<this->GetGeometry()[0].Y()<<std::endl;
      std::cout<<" MeanValueMass="<<meanValueMass;
      std::cout<<"\t MeanValueMaterial= "<<meanValueStiff;
      bulkCoefficient=timeStep;
    }
  }


  template<>
  void TwoStepUpdatedLagrangianVPFluidElement<3>::ComputeBulkReductionCoefficient(MatrixType MassMatrix,
										  MatrixType StiffnessMatrix,
										  double& meanValueStiff,
										  double& bulkCoefficient,
										  double timeStep)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();
    IndexType FirstRow = 0;
    IndexType FirstCol = 0;
    double meanValueMass=0;
    double countStiff=0;
    double countMass=0;
    for (SizeType j = 0; j < NumNodes; ++j)
      {
        for (SizeType i = 0; i < NumNodes; ++i)
	  {
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow,FirstCol));
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow,FirstCol+1));
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow,FirstCol+2));
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow+1,FirstCol));
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow+1,FirstCol+1));
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow+1,FirstCol+2));
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow+2,FirstCol));
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow+2,FirstCol+1));
	    meanValueStiff+= fabs(StiffnessMatrix(FirstRow+2,FirstCol+2));
	    countStiff+=9.0;
	    
	    meanValueMass+= fabs(MassMatrix(FirstRow,FirstCol));
	    meanValueMass+= fabs(MassMatrix(FirstRow+1,FirstCol+1));
	    meanValueMass+= fabs(MassMatrix(FirstRow+2,FirstCol+2));
	    countMass+=3.0;

	    // Update Counter
            FirstRow += 3;
	  }
	FirstRow = 0;
        FirstCol += 3;
      }

    meanValueStiff*=1.0/countStiff;
    meanValueMass*=1.0/countMass;
    
    if(meanValueMass!=0 && meanValueStiff!=0){
      bulkCoefficient=meanValueMass*2.0/timeStep/meanValueStiff;
    }else{
      std::cout<<" DANGEROUS ELEMENT!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
      std::cout<<" MeanValueMass="<<meanValueMass;
      std::cout<<"\t MeanValueMaterial= "<<meanValueStiff;
      bulkCoefficient=timeStep;
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
											 const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    double coeff=1.0+TDim;
    if((NumNodes==3 && TDim==2) || (NumNodes==4 && TDim==3)){
      for (SizeType i = 0; i < NumNodes; ++i)
	{
	  // LHS contribution
	  double Mij  = Weight /coeff;
	  BulkVelMatrix(i,i) +=  Mij;
	}
    }else{
      std::cout<<"... ComputeBulkMatrixForPressureVelLump TO IMPLEMENT"<<std::endl;
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
											 const double Weight)
  {
    const SizeType NumNodes = this->GetGeometry().PointsNumber();

    double coeff=1.0+TDim;
    if((NumNodes==3 && TDim==2) || (NumNodes==4 && TDim==3)){
      for (SizeType i = 0; i < NumNodes; ++i)
	{
	  // LHS contribution
	  double Mij  = Weight /coeff;
	  BulkAccMatrix(i,i) +=  Mij;
	}
    }else{
      std::cout<<"... ComputeBulkMatrixForPressureAccLump TO IMPLEMENT"<<std::endl;
    }
  }

 

  template< >
  void TwoStepUpdatedLagrangianVPFluidElement<2>::ComputeBoundLHSMatrix(Matrix& BoundLHSMatrix,
									const ShapeFunctionsType& rN,
									const double Weight)
  {
    GeometryType& rGeom = this->GetGeometry();
    //const SizeType NumNodes = rGeom.PointsNumber();

    // for (SizeType i = 0; i < (NumNodes-1); i++)
    //   {
    //     for (SizeType j = (i+1); j < NumNodes; j++)
    // 	  {
    // 	    if(rGeom[i].Is(FREE_SURFACE) && rGeom[j].Is(FREE_SURFACE)){
    // 	      if(rGeom[i].IsNot(INLET))
    // 		BoundLHSMatrix(i,i) +=  Weight / 3.0;
    // 	      if(rGeom[j].IsNot(INLET))
    // 		BoundLHSMatrix(j,j) +=  Weight / 3.0;
    // 	    }
    // 	  }

    //   }

    // for (SizeType i = 0; i < (NumNodes-1); i++)
    //   {
    //     for (SizeType j = (i+1); j < NumNodes; j++)
    // 	  {
    // 	    if(rGeom[i].Is(FREE_SURFACE) && rGeom[j].Is(FREE_SURFACE)){
    // 	      if(rGeom[i].IsNot(INLET)){
    // 		BoundLHSMatrix(i,i) +=  Weight*rN[i]*rN[i];
    // 		BoundLHSMatrix(i,j) +=  Weight*rN[i]*rN[j];
    // 	      }
    // 	      if(rGeom[j].IsNot(INLET)){
    // 		BoundLHSMatrix(j,j) +=  Weight*rN[j]*rN[j];
    // 		BoundLHSMatrix(j,i) +=  Weight*rN[j]*rN[i];
    // 	      }
    // 	    }
    // 	  }

    //   }

    // if(this->Is(TO_ERASE)){
    //   WeakPointerVector<Node<3> >& rN0 = rGeom[0].GetValue(NEIGHBOUR_NODES);
    //   WeakPointerVector<Node<3> >& rN1 = rGeom[1].GetValue(NEIGHBOUR_NODES);
    //   WeakPointerVector<Node<3> >& rN2 = rGeom[2].GetValue(NEIGHBOUR_NODES);
    //   if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE) &&
    // 	 (rN0.size()<NumNodes || rN1.size()<NumNodes)){
    // 	if(rGeom[0].IsNot(INLET))
    // 	  BoundLHSMatrix(0,0) +=  Weight / 3.0;
    // 	if(rGeom[1].IsNot(INLET))
    // 	  BoundLHSMatrix(1,1) +=  Weight / 3.0;
    //   }else if(rGeom[0].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE) &&
    // 	 (rN0.size()<NumNodes || rN2.size()<NumNodes)){
    // 	if(rGeom[0].IsNot(INLET))
    // 	  BoundLHSMatrix(0,0) +=  Weight / 3.0;
    // 	if(rGeom[2].IsNot(INLET))
    // 	  BoundLHSMatrix(2,2) +=  Weight / 3.0;
    //   }else if(rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE) &&
    // 	 (rN2.size()<NumNodes || rN1.size()<NumNodes)){
    // 	if(rGeom[1].IsNot(INLET))
    // 	  BoundLHSMatrix(1,1) +=  Weight / 3.0;
    // 	if(rGeom[2].IsNot(INLET))
    // 	  BoundLHSMatrix(2,2) +=  Weight / 3.0;
    //   }
    // }else{

      if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)){
	if(rGeom[0].IsNot(INLET))
	  BoundLHSMatrix(0,0) +=  Weight / 3.0;
	if(rGeom[1].IsNot(INLET))
	  BoundLHSMatrix(1,1) +=  Weight / 3.0;
      }
      if(rGeom[0].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)){
	if(rGeom[0].IsNot(INLET))
	  BoundLHSMatrix(0,0) +=  Weight / 3.0;
	if(rGeom[2].IsNot(INLET))
	  BoundLHSMatrix(2,2) +=  Weight / 3.0;
      }
      if(rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)){
	if(rGeom[1].IsNot(INLET))
	  BoundLHSMatrix(1,1) +=  Weight / 3.0;
	if(rGeom[2].IsNot(INLET))
	  BoundLHSMatrix(2,2) +=  Weight / 3.0;
      }
    // }

  }

  template<  >
  void TwoStepUpdatedLagrangianVPFluidElement<3>::ComputeBoundLHSMatrix(Matrix& BoundLHSMatrix,
									const ShapeFunctionsType& rN,
									const double Weight)
  {
    GeometryType& rGeom = this->GetGeometry(); 
    //const SizeType NumNodes = rGeom.PointsNumber();

    // for (SizeType i = 0; i < (NumNodes-2); i++)
    //   {
    //     for (SizeType j = (i+1); j < (NumNodes-1); j++)
    // 	  {
    // 	    for (SizeType k = (j+1); k < NumNodes; k++)
    // 	      {
    // 		if(rGeom[i].Is(FREE_SURFACE) && rGeom[j].Is(FREE_SURFACE) && rGeom[k].Is(FREE_SURFACE)){
    // 		  if(rGeom[i].IsNot(INLET)){
    // 		    BoundLHSMatrix(i,i) +=  Weight*rN[i]*rN[i];
    // 		    BoundLHSMatrix(i,j) +=  Weight*rN[i]*rN[j];
    // 		    BoundLHSMatrix(i,k) +=  Weight*rN[i]*rN[k];
    // 		  }
    // 		  if(rGeom[j].IsNot(INLET)){
    // 		    BoundLHSMatrix(j,i) +=  Weight*rN[j]*rN[i];
    // 		    BoundLHSMatrix(j,j) +=  Weight*rN[j]*rN[j];
    // 		    BoundLHSMatrix(j,k) +=  Weight*rN[j]*rN[k];
    // 		  }
    // 		  if(rGeom[k].IsNot(INLET)){
    // 		    BoundLHSMatrix(k,i) +=  Weight*rN[k]*rN[i];
    // 		    BoundLHSMatrix(k,j) +=  Weight*rN[k]*rN[j];
    // 		    BoundLHSMatrix(k,k) +=  Weight*rN[k]*rN[k];
    // 		  }
    // 		}
    // 	      }
    // 	  }

    //   }


    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)){
      if(rGeom[0].IsNot(INLET))
    	BoundLHSMatrix(0,0) +=  Weight / 4.0;
      if(rGeom[1].IsNot(INLET))
    	BoundLHSMatrix(1,1) +=  Weight / 4.0;
      if(rGeom[2].IsNot(INLET))
    	BoundLHSMatrix(2,2) +=  Weight / 4.0;
    }
    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
      if(rGeom[0].IsNot(INLET))
    	BoundLHSMatrix(0,0) +=  Weight / 4.0;
      if(rGeom[1].IsNot(INLET))
    	BoundLHSMatrix(1,1) +=  Weight / 4.0;
      if(rGeom[3].IsNot(INLET))
    	BoundLHSMatrix(3,3) +=  Weight / 4.0;
    }
    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
      if(rGeom[0].IsNot(INLET))
    	BoundLHSMatrix(0,0) +=  Weight / 4.0;
      if(rGeom[2].IsNot(INLET))
    	BoundLHSMatrix(2,2) +=  Weight / 4.0;
      if(rGeom[3].IsNot(INLET))
    	BoundLHSMatrix(3,3) +=  Weight / 4.0;
    }
    if(rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
      if(rGeom[1].IsNot(INLET))
    	BoundLHSMatrix(1,1) +=  Weight / 4.0;
      if(rGeom[2].IsNot(INLET))
    	BoundLHSMatrix(2,2) +=  Weight / 4.0;
      if(rGeom[3].IsNot(INLET))
    	BoundLHSMatrix(3,3) +=  Weight / 4.0;
    }

  }
 


  template< >
  void TwoStepUpdatedLagrangianVPFluidElement<2>::ComputeBoundRHSVector(VectorType& BoundRHSVector,
									const ShapeFunctionsType& rN,
									const double TimeStep,
									const double BoundRHSCoeffAcc,
									const double BoundRHSCoeffDev)
  {
    GeometryType& rGeom = this->GetGeometry();
    //const SizeType NumNodes = rGeom.PointsNumber();
    array_1d<double, 3>  AccA(3,0.0);
    array_1d<double, 3>  AccB(3,0.0);

    // for (SizeType i = 0; i < (NumNodes-1); i++)
    //   {
    // 	for (SizeType j = (i+1); j < NumNodes; j++)
    // 	  {
    // 	    if(rGeom[i].Is(FREE_SURFACE) && rGeom[j].Is(FREE_SURFACE)){
    // 	      AccA= 0.5/TimeStep*(rGeom[i].FastGetSolutionStepValue(VELOCITY,0)-rGeom[i].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[i].FastGetSolutionStepValue(ACCELERATION,1); 
    // 	      AccB= 0.5/TimeStep*(rGeom[j].FastGetSolutionStepValue(VELOCITY,0)-rGeom[j].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[j].FastGetSolutionStepValue(ACCELERATION,1); 
    // 	      const array_1d<double, 3> &NormalA    = rGeom[i].FastGetSolutionStepValue(NORMAL);
    // 	      const array_1d<double, 3> &NormalB    = rGeom[j].FastGetSolutionStepValue(NORMAL);
    // 	      double coeff=3.0;
    // 	      if(rGeom[i].IsNot(INLET)) //to change into moving wall!!!!!
    // 		BoundRHSVector[i] += (BoundRHSCoeffAcc*(AccA[0]*NormalA[0]+AccA[1]*NormalA[1]) + 
    // 				      BoundRHSCoeffDev)/coeff ;
    // 	      if(rGeom[j].IsNot(INLET))
    // 		BoundRHSVector[j] += (BoundRHSCoeffAcc*(AccB[0]*NormalB[0]+AccB[1]*NormalB[1]) + 
    // 				      BoundRHSCoeffDev)/coeff ;
    // 	    }
    // 	  }

    //   }

    // for (SizeType i = 0; i < (NumNodes-1); i++)
    //   {
    // 	for (SizeType j = (i+1); j < NumNodes; j++)
    // 	  {
    // 	    if(rGeom[i].Is(FREE_SURFACE) && rGeom[j].Is(FREE_SURFACE)){
    // 	      AccA= 0.5/TimeStep*(rGeom[i].FastGetSolutionStepValue(VELOCITY,0)-rGeom[i].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[i].FastGetSolutionStepValue(ACCELERATION,1); 
    // 	      AccB= 0.5/TimeStep*(rGeom[j].FastGetSolutionStepValue(VELOCITY,0)-rGeom[j].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[j].FastGetSolutionStepValue(ACCELERATION,1); 
    // 	      const array_1d<double, 3> &NormalA    = rGeom[i].FastGetSolutionStepValue(NORMAL);
    // 	      const array_1d<double, 3> &NormalB    = rGeom[j].FastGetSolutionStepValue(NORMAL);
    // 	      if(rGeom[i].IsNot(INLET)) 
    // 		BoundRHSVector[i] += (BoundRHSCoeffAcc*(AccA[0]*NormalA[0]+AccA[1]*NormalA[1]) + 
    // 				      BoundRHSCoeffDev) * rN[i];
    // 	      if(rGeom[j].IsNot(INLET))
    // 		BoundRHSVector[j] += (BoundRHSCoeffAcc*(AccB[0]*NormalB[0]+AccB[1]*NormalB[1]) + 
    // 				      BoundRHSCoeffDev) * rN[j] ;
    // 	    }
    // 	  }

    //   }


    // if(this->Is(TO_ERASE)){
    //   WeakPointerVector<Node<3> >& rN0 = rGeom[0].GetValue(NEIGHBOUR_NODES);
    //   WeakPointerVector<Node<3> >& rN1 = rGeom[1].GetValue(NEIGHBOUR_NODES);
    //   WeakPointerVector<Node<3> >& rN2 = rGeom[2].GetValue(NEIGHBOUR_NODES);

    //   if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE) && (rN0.size()<NumNodes || rN1.size()<NumNodes)){ 
    // 	AccA= 0.5/TimeStep*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
    // 	AccB= 0.5/TimeStep*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
    // 	const array_1d<double, 3> &NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
    // 	const array_1d<double, 3> &NormalB    = rGeom[1].FastGetSolutionStepValue(NORMAL);
    // 	if(rGeom[0].IsNot(INLET)) //to change into moving wall!!!!!
    // 	  BoundRHSVector[0] += (BoundRHSCoeffAcc*(AccA[0]*NormalA[0]+AccA[1]*NormalA[1]) + BoundRHSCoeffDev)/3.0;
    // 	if(rGeom[1].IsNot(INLET))
    // 	  BoundRHSVector[1] += (BoundRHSCoeffAcc*(AccB[0]*NormalB[0]+AccB[1]*NormalB[1]) + BoundRHSCoeffDev)/3.0;
    //   }else  if(rGeom[0].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE) && (rN0.size()<NumNodes || rN2.size()<NumNodes) ){
    // 	AccA= 0.5/TimeStep*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
    // 	AccB= 0.5/TimeStep*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
    // 	const array_1d<double, 3> &NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
    // 	const array_1d<double, 3> &NormalB    = rGeom[2].FastGetSolutionStepValue(NORMAL);
    // 	if(rGeom[0].IsNot(INLET))
    // 	  BoundRHSVector[0] += (BoundRHSCoeffAcc*(AccA[0]*NormalA[0]+AccA[1]*NormalA[1]) + BoundRHSCoeffDev)/3.0;
    // 	if(rGeom[2].IsNot(INLET))   
    // 	  BoundRHSVector[2] += (BoundRHSCoeffAcc*(AccB[0]*NormalB[0]+AccB[1]*NormalB[1]) + BoundRHSCoeffDev)/3.0;
    //   }else  if(rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && (rN1.size()<NumNodes || rN2.size()<NumNodes) ){
    // 	AccA= 0.5/TimeStep*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
    // 	AccB= 0.5/TimeStep*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
    // 	const array_1d<double, 3> &NormalA    = rGeom[1].FastGetSolutionStepValue(NORMAL);
    // 	const array_1d<double, 3> &NormalB    = rGeom[2].FastGetSolutionStepValue(NORMAL);
    // 	if(rGeom[1].IsNot(INLET))
    // 	  BoundRHSVector[1] += (BoundRHSCoeffAcc*(AccA[0]*NormalA[0]+AccA[1]*NormalA[1]) + BoundRHSCoeffDev)/3.0;
    // 	if(rGeom[2].IsNot(INLET))
    // 	  BoundRHSVector[2] += (BoundRHSCoeffAcc*(AccB[0]*NormalB[0]+AccB[1]*NormalB[1]) + BoundRHSCoeffDev)/3.0;
    //   }

    // }else{
      const double factor = 0.5/TimeStep;
      const double one_third = 1.0/3.0;

      if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE) ){ 
	noalias(AccA)= factor*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
	noalias(AccB)= factor*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
	const array_1d<double, 3> &NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
	const array_1d<double, 3> &NormalB    = rGeom[1].FastGetSolutionStepValue(NORMAL);
	if(rGeom[0].IsNot(INLET)) //to change into moving wall!!!!!
	  BoundRHSVector[0] += one_third * (BoundRHSCoeffAcc*(AccA[0]*NormalA[0]+AccA[1]*NormalA[1]) + BoundRHSCoeffDev);
	if(rGeom[1].IsNot(INLET))
	  BoundRHSVector[1] += one_third * (BoundRHSCoeffAcc*(AccB[0]*NormalB[0]+AccB[1]*NormalB[1]) + BoundRHSCoeffDev);
      }
      if(rGeom[0].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE) ){
	noalias(AccA)= factor*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
	noalias(AccB)= factor*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
	const array_1d<double, 3> &NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
	const array_1d<double, 3> &NormalB    = rGeom[2].FastGetSolutionStepValue(NORMAL);
	if(rGeom[0].IsNot(INLET))
	  BoundRHSVector[0] += one_third * (BoundRHSCoeffAcc*(AccA[0]*NormalA[0]+AccA[1]*NormalA[1]) + BoundRHSCoeffDev);
	if(rGeom[2].IsNot(INLET))   
	  BoundRHSVector[2] += one_third * (BoundRHSCoeffAcc*(AccB[0]*NormalB[0]+AccB[1]*NormalB[1]) + BoundRHSCoeffDev);
      }
      if(rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE) ){
	noalias(AccA)= factor*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
	noalias(AccB)= factor*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
	const array_1d<double, 3> &NormalA    = rGeom[1].FastGetSolutionStepValue(NORMAL);
	const array_1d<double, 3> &NormalB    = rGeom[2].FastGetSolutionStepValue(NORMAL);
	if(rGeom[1].IsNot(INLET))
	  BoundRHSVector[1] += one_third * (BoundRHSCoeffAcc*(AccA[0]*NormalA[0]+AccA[1]*NormalA[1]) + BoundRHSCoeffDev);
	if(rGeom[2].IsNot(INLET))
	  BoundRHSVector[2] += one_third * (BoundRHSCoeffAcc*(AccB[0]*NormalB[0]+AccB[1]*NormalB[1]) + BoundRHSCoeffDev);
      }
    // }

  }



  template< >
  void TwoStepUpdatedLagrangianVPFluidElement<3>::ComputeBoundRHSVector(VectorType& BoundRHSVector,
									const ShapeFunctionsType& rN,
									const double TimeStep,
									const double BoundRHSCoeffAcc,
									const double BoundRHSCoeffDev)
  {
    GeometryType& rGeom = this->GetGeometry();
    //const SizeType NumNodes = rGeom.PointsNumber();
    array_1d<double, 3>  AccA(3,0.0);
    array_1d<double, 3>  AccB(3,0.0);
    array_1d<double, 3>  AccC(3,0.0);


    // for (SizeType i = 0; i < (NumNodes-2); i++)
    //   {
    // 	for (SizeType j = (i+1); j < (NumNodes-1); j++)
    // 	  {
    // 	    for (SizeType k = (j+1); k < NumNodes; k++)
    // 	      {
    // 		if(rGeom[i].Is(FREE_SURFACE) && rGeom[j].Is(FREE_SURFACE) && rGeom[k].Is(FREE_SURFACE)){
    // 		  AccA= 0.5/TimeStep*(rGeom[i].FastGetSolutionStepValue(VELOCITY,0)-rGeom[i].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[i].FastGetSolutionStepValue(ACCELERATION,1); 
    // 		  AccB= 0.5/TimeStep*(rGeom[j].FastGetSolutionStepValue(VELOCITY,0)-rGeom[j].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[j].FastGetSolutionStepValue(ACCELERATION,1); 
    // 		  AccC= 0.5/TimeStep*(rGeom[k].FastGetSolutionStepValue(VELOCITY,0)-rGeom[k].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[k].FastGetSolutionStepValue(ACCELERATION,1); 
		  
    // 		  const array_1d<double, 3> &NormalA    = rGeom[i].FastGetSolutionStepValue(NORMAL);
    // 		  const array_1d<double, 3> &NormalB    = rGeom[j].FastGetSolutionStepValue(NORMAL);
    // 		  const array_1d<double, 3> &NormalC    = rGeom[k].FastGetSolutionStepValue(NORMAL);
    // 		  if(rGeom[i].IsNot(INLET)) 
    // 		    BoundRHSVector[i] += (BoundRHSCoeffAcc*(AccA[0]*NormalA[0] + AccA[1]*NormalA[1] + AccA[2]*NormalA[2]) + 
    // 					  BoundRHSCoeffDev) * rN[i];
    // 		  if(rGeom[j].IsNot(INLET))
    // 		    BoundRHSVector[j] += (BoundRHSCoeffAcc*(AccB[0]*NormalB[0] + AccB[1]*NormalB[1] + AccB[2]*NormalB[2]) + 
    // 					  BoundRHSCoeffDev) * rN[j] ;
    // 		  if(rGeom[k].IsNot(INLET))
    // 		    BoundRHSVector[k] += (BoundRHSCoeffAcc*(AccC[0]*NormalC[0] + AccC[1]*NormalC[1] + AccC[2]*NormalC[2]) + 
    // 					  BoundRHSCoeffDev) * rN[k] ;
    // 		}
    // 	      }
    // 	  }

    //   }
    const double factor = 0.5/TimeStep;

    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)){      
      noalias(AccA)= factor*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
      noalias(AccB)= factor*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
      noalias(AccC)= factor*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
      const array_1d<double,3> &NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
      const array_1d<double,3> &NormalB    = rGeom[1].FastGetSolutionStepValue(NORMAL);
      const array_1d<double,3> &NormalC    = rGeom[2].FastGetSolutionStepValue(NORMAL);
      if(rGeom[0].IsNot(INLET))
	BoundRHSVector[0] += 0.25 * (BoundRHSCoeffAcc*(AccA[0]*NormalA[0] + AccA[1]*NormalA[1] + AccA[2]*NormalA[2]) +
			      BoundRHSCoeffDev);
      if(rGeom[1].IsNot(INLET))
	BoundRHSVector[1] += 0.25 * (BoundRHSCoeffAcc*(AccB[0]*NormalB[0] + AccB[1]*NormalB[1] + AccB[2]*NormalB[2]) +
			      BoundRHSCoeffDev);
      if(rGeom[2].IsNot(INLET))
	BoundRHSVector[2] += 0.25 * (BoundRHSCoeffAcc*(AccC[0]*NormalC[0] + AccC[1]*NormalC[1] + AccC[2]*NormalC[2]) +
			      BoundRHSCoeffDev);
    }
    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
      noalias(AccA)= factor*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
      noalias(AccB)= factor*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
      noalias(AccC)= factor*(rGeom[3].FastGetSolutionStepValue(VELOCITY,0)-rGeom[3].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION,1); 
      const array_1d<double,3> &NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
      const array_1d<double,3> &NormalB    = rGeom[1].FastGetSolutionStepValue(NORMAL);
      const array_1d<double,3> &NormalC    = rGeom[3].FastGetSolutionStepValue(NORMAL);
      if(rGeom[0].IsNot(INLET))
	BoundRHSVector[0] += 0.25 * (BoundRHSCoeffAcc*(AccA[0]*NormalA[0] + AccA[1]*NormalA[1] + AccA[2]*NormalA[2]) +
			      BoundRHSCoeffDev);
      if(rGeom[1].IsNot(INLET))
	BoundRHSVector[1] += 0.25 * (BoundRHSCoeffAcc*(AccB[0]*NormalB[0] + AccB[1]*NormalB[1] + AccB[2]*NormalB[2]) +
			      BoundRHSCoeffDev);
      if(rGeom[3].IsNot(INLET))
	BoundRHSVector[3] += 0.25 * (BoundRHSCoeffAcc*(AccC[0]*NormalC[0] + AccC[1]*NormalC[1] + AccC[2]*NormalC[2]) +
			      BoundRHSCoeffDev);
    }
    if(rGeom[0].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
      noalias(AccA)= factor*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
      noalias(AccB)= factor*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
      noalias(AccC)= factor*(rGeom[3].FastGetSolutionStepValue(VELOCITY,0)-rGeom[3].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION,1); 
      const array_1d<double,3> &NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
      const array_1d<double,3> &NormalB    = rGeom[2].FastGetSolutionStepValue(NORMAL);
      const array_1d<double,3> &NormalC    = rGeom[3].FastGetSolutionStepValue(NORMAL);
      if(rGeom[0].IsNot(INLET))
	BoundRHSVector[0] += 0.25 * (BoundRHSCoeffAcc*(AccA[0]*NormalA[0] + AccA[1]*NormalA[1] + AccA[2]*NormalA[2]) +
			      BoundRHSCoeffDev);
      if(rGeom[2].IsNot(INLET))
	BoundRHSVector[2] += 0.25 * (BoundRHSCoeffAcc*(AccB[0]*NormalB[0] + AccB[1]*NormalB[1] + AccB[2]*NormalB[2]) +
			      BoundRHSCoeffDev);
      if(rGeom[3].IsNot(INLET))
	BoundRHSVector[3] += 0.25 * (BoundRHSCoeffAcc*(AccC[0]*NormalC[0] + AccC[1]*NormalC[1] + AccC[2]*NormalC[2]) +
			      BoundRHSCoeffDev);
    }
    if(rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){      
      noalias(AccA)= factor*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
      noalias(AccB)= factor*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
      noalias(AccC)= factor*(rGeom[3].FastGetSolutionStepValue(VELOCITY,0)-rGeom[3].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION,1); 
      const array_1d<double,3> &NormalA    = rGeom[1].FastGetSolutionStepValue(NORMAL);
      const array_1d<double,3> &NormalB    = rGeom[2].FastGetSolutionStepValue(NORMAL);
      const array_1d<double,3> &NormalC    = rGeom[3].FastGetSolutionStepValue(NORMAL);
      if(rGeom[1].IsNot(INLET))
	BoundRHSVector[1] += 0.25 * (BoundRHSCoeffAcc*(AccA[0]*NormalA[0] + AccA[1]*NormalA[1] + AccA[2]*NormalA[2]) +
			      BoundRHSCoeffDev);
      if(rGeom[2].IsNot(INLET))
	BoundRHSVector[2] += 0.25 * (BoundRHSCoeffAcc*(AccB[0]*NormalB[0] + AccB[1]*NormalB[1] + AccB[2]*NormalB[2]) +
			      BoundRHSCoeffDev);
      if(rGeom[3].IsNot(INLET))
	BoundRHSVector[3] += 0.25 * (BoundRHSCoeffAcc*(AccC[0]*NormalC[0] + AccC[1]*NormalC[1] + AccC[2]*NormalC[2]) +
			      BoundRHSCoeffDev);
    }

    // if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)){      
    //   noalias(AccA)= factor*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
    //   noalias(AccB)= factor*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
    //   noalias(AccC)= factor*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
    //   const array_1d<double,3> &NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
    //   const array_1d<double,3> &NormalB    = rGeom[1].FastGetSolutionStepValue(NORMAL);
    //   const array_1d<double,3> &NormalC    = rGeom[2].FastGetSolutionStepValue(NORMAL);
    //   if(rGeom[0].IsNot(INLET))
    // 	BoundRHSVector[0] += 0.25 * (BoundRHSCoeffAcc*(AccA[0]*NormalA[0] + AccA[1]*NormalA[1] + AccA[2]*NormalA[2]) +
    // 			      BoundRHSCoeffDev);
    //   if(rGeom[1].IsNot(INLET))
    // 	BoundRHSVector[1] += 0.25 * (BoundRHSCoeffAcc*(AccB[0]*NormalB[0] + AccB[1]*NormalB[1] + AccB[2]*NormalB[2]) +
    // 			      BoundRHSCoeffDev);
    //   if(rGeom[2].IsNot(INLET))
    // 	BoundRHSVector[2] += 0.25 * (BoundRHSCoeffAcc*(AccC[0]*NormalC[0] + AccC[1]*NormalC[1] + AccC[2]*NormalC[2]) +
    // 			      BoundRHSCoeffDev);
    // }else if(rGeom[0].Is(FREE_SURFACE)  && rGeom[1].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
    //   noalias(AccA)= factor*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
    //   noalias(AccB)= factor*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
    //   noalias(AccC)= factor*(rGeom[3].FastGetSolutionStepValue(VELOCITY,0)-rGeom[3].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION,1); 
    //   const array_1d<double,3> &NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
    //   const array_1d<double,3> &NormalB    = rGeom[1].FastGetSolutionStepValue(NORMAL);
    //   const array_1d<double,3> &NormalC    = rGeom[3].FastGetSolutionStepValue(NORMAL);
    //   if(rGeom[0].IsNot(INLET))
    // 	BoundRHSVector[0] += 0.25 * (BoundRHSCoeffAcc*(AccA[0]*NormalA[0] + AccA[1]*NormalA[1] + AccA[2]*NormalA[2]) +
    // 			      BoundRHSCoeffDev);
    //   if(rGeom[1].IsNot(INLET))
    // 	BoundRHSVector[1] += 0.25 * (BoundRHSCoeffAcc*(AccB[0]*NormalB[0] + AccB[1]*NormalB[1] + AccB[2]*NormalB[2]) +
    // 			      BoundRHSCoeffDev);
    //   if(rGeom[3].IsNot(INLET))
    // 	BoundRHSVector[3] += 0.25 * (BoundRHSCoeffAcc*(AccC[0]*NormalC[0] + AccC[1]*NormalC[1] + AccC[2]*NormalC[2]) +
    // 			      BoundRHSCoeffDev);
    // }else if(rGeom[0].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){
    //   noalias(AccA)= factor*(rGeom[0].FastGetSolutionStepValue(VELOCITY,0)-rGeom[0].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[0].FastGetSolutionStepValue(ACCELERATION,1); 
    //   noalias(AccB)= factor*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
    //   noalias(AccC)= factor*(rGeom[3].FastGetSolutionStepValue(VELOCITY,0)-rGeom[3].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION,1); 
    //   const array_1d<double,3> &NormalA    = rGeom[0].FastGetSolutionStepValue(NORMAL);
    //   const array_1d<double,3> &NormalB    = rGeom[2].FastGetSolutionStepValue(NORMAL);
    //   const array_1d<double,3> &NormalC    = rGeom[3].FastGetSolutionStepValue(NORMAL);
    //   if(rGeom[0].IsNot(INLET))
    // 	BoundRHSVector[0] += 0.25 * (BoundRHSCoeffAcc*(AccA[0]*NormalA[0] + AccA[1]*NormalA[1] + AccA[2]*NormalA[2]) +
    // 			      BoundRHSCoeffDev);
    //   if(rGeom[2].IsNot(INLET))
    // 	BoundRHSVector[2] += 0.25 * (BoundRHSCoeffAcc*(AccB[0]*NormalB[0] + AccB[1]*NormalB[1] + AccB[2]*NormalB[2]) +
    // 			      BoundRHSCoeffDev);
    //   if(rGeom[3].IsNot(INLET))
    // 	BoundRHSVector[3] += 0.25 * (BoundRHSCoeffAcc*(AccC[0]*NormalC[0] + AccC[1]*NormalC[1] + AccC[2]*NormalC[2]) +
    // 			      BoundRHSCoeffDev);
    // }else if(rGeom[1].Is(FREE_SURFACE)  && rGeom[2].Is(FREE_SURFACE)  && rGeom[3].Is(FREE_SURFACE)){      
    //   noalias(AccA)= factor*(rGeom[1].FastGetSolutionStepValue(VELOCITY,0)-rGeom[1].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[1].FastGetSolutionStepValue(ACCELERATION,1); 
    //   noalias(AccB)= factor*(rGeom[2].FastGetSolutionStepValue(VELOCITY,0)-rGeom[2].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[2].FastGetSolutionStepValue(ACCELERATION,1); 
    //   noalias(AccC)= factor*(rGeom[3].FastGetSolutionStepValue(VELOCITY,0)-rGeom[3].FastGetSolutionStepValue(VELOCITY,1)) - rGeom[3].FastGetSolutionStepValue(ACCELERATION,1); 
    //   const array_1d<double,3> &NormalA    = rGeom[1].FastGetSolutionStepValue(NORMAL);
    //   const array_1d<double,3> &NormalB    = rGeom[2].FastGetSolutionStepValue(NORMAL);
    //   const array_1d<double,3> &NormalC    = rGeom[3].FastGetSolutionStepValue(NORMAL);
    //   if(rGeom[1].IsNot(INLET))
    // 	BoundRHSVector[1] += 0.25 * (BoundRHSCoeffAcc*(AccA[0]*NormalA[0] + AccA[1]*NormalA[1] + AccA[2]*NormalA[2]) +
    // 			      BoundRHSCoeffDev);
    //   if(rGeom[2].IsNot(INLET))
    // 	BoundRHSVector[2] += 0.25 * (BoundRHSCoeffAcc*(AccB[0]*NormalB[0] + AccB[1]*NormalB[1] + AccB[2]*NormalB[2]) +
    // 			      BoundRHSCoeffDev);
    //   if(rGeom[3].IsNot(INLET))
    // 	BoundRHSVector[3] += 0.25 * (BoundRHSCoeffAcc*(AccC[0]*NormalC[0] + AccC[1]*NormalC[1] + AccC[2]*NormalC[2]) +
    // 			      BoundRHSCoeffDev);
    // }



  }



  template <unsigned int TDim>
  void TwoStepUpdatedLagrangianVPFluidElement<TDim>::CalculateTauFIC(double& Tau,
								     double ElemSize,
								     const double Density,
								     const double Viscosity,
								     const ProcessInfo& rCurrentProcessInfo)
  {
    double DeltaTime = rCurrentProcessInfo.GetValue(DELTA_TIME);
    if(rCurrentProcessInfo.GetValue(DELTA_TIME)<rCurrentProcessInfo.GetValue(PREVIOUS_DELTA_TIME)){
      DeltaTime = 0.5*rCurrentProcessInfo.GetValue(DELTA_TIME)+0.5*rCurrentProcessInfo.GetValue(PREVIOUS_DELTA_TIME);
    }

    double MeanVelocity=0;
    this->CalcMeanVelocity(MeanVelocity,0);

    Tau = 1.0 / (2.0 * Density *(0.5 * MeanVelocity / ElemSize + 0.5/DeltaTime) +  8.0 * Viscosity / (ElemSize * ElemSize) ); 
  
    if(Tau<0.0000001){
      Tau=0.0000001;
    }

    if(MeanVelocity==0){
      Tau=0;
    }
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
										   const double Weight,
										   const ShapeFunctionDerivativesType& rDN_DX,
										   const SizeType i)
  {

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



  // template< unsigned int TDim>
  // bool TwoStepUpdatedLagrangianVPFluidElement<TDim>::CalcMechanicsUpdated(ElementalVariables & rElementalVariables,
  // 									  const ProcessInfo& rCurrentProcessInfo,
  // 									  const ShapeFunctionDerivativesType& rDN_DX,
  // 									  unsigned int g)
  // {
  //   double theta=this->GetThetaMomentum();
  //   const double TimeStep=rCurrentProcessInfo[DELTA_TIME];
  //   bool computeElement=this->CalcStrainRate(rElementalVariables,rCurrentProcessInfo,rDN_DX,theta);
  //   this->CalcElasticPlasticCauchySplitted(rElementalVariables,TimeStep,g);
  //   return computeElement;
  // } 


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
  void TwoStepUpdatedLagrangianVPFluidElement<2>:: CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables,double TimeStep, unsigned int g)
  {

    double Density  = 0;
    double CurrSecondLame  = 0;
    double CurrBulkModulus = 0;

    this->ComputeMaterialParameters(Density,CurrSecondLame,CurrBulkModulus,TimeStep);
 
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

  }


  template < > 
  void TwoStepUpdatedLagrangianVPFluidElement<3>:: CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables, double TimeStep, unsigned int g)
  {

    double Density  = 0;
    double CurrSecondLame  = 0;
    double CurrBulkModulus = 0;

    this->ComputeMaterialParameters(Density,CurrSecondLame,CurrBulkModulus,TimeStep);
 
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

  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPFluidElement<TDim>::CalculateLocalContinuityEqForPressure(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
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

    // Shape functions and integration points
    ShapeFunctionDerivativesArrayType DN_DX;
    Matrix NContainer;
    VectorType GaussWeights;
    this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);
    const unsigned int NumGauss = GaussWeights.size();

    // MatrixType BulkVelMatrix = ZeroMatrix(NumNodes,NumNodes);
    // MatrixType BulkAccMatrix = ZeroMatrix(NumNodes,NumNodes);

    double TimeStep=rCurrentProcessInfo[DELTA_TIME];
    double theta=this->GetThetaContinuity();
    double ElemSize = this->ElementSize();

    ElementalVariables rElementalVariables;
    this->InitializeElementalVariables(rElementalVariables);

    double Density = 0;
    double DeviatoricCoeff = 0;
    double VolumetricCoeff = 0;
    this->ComputeMaterialParameters(Density,DeviatoricCoeff,VolumetricCoeff,TimeStep);
    
    double Tau=0;
    this->CalculateTauFIC(Tau,ElemSize,Density,DeviatoricCoeff,rCurrentProcessInfo);

    double totalVolume=0;
    bool computeElement=false;
    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; g++)
      {
	const double GaussWeight = GaussWeights[g];
	totalVolume+=GaussWeight;
	const ShapeFunctionsType& N = row(NContainer,g);
	const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];
	// computeElement=this->CalcStrainRate(rElementalVariables,rCurrentProcessInfo,rDN_DX,theta);
	computeElement=this->CalcCompleteStrainRate(rElementalVariables,rCurrentProcessInfo,rDN_DX,theta);

	if(computeElement==true){
	  // Evaluate required variables at the integration point

	  // double BulkCoeff =GaussWeight/(VolumetricCoeff);
	  // this->ComputeBulkMatrixForPressureVel(BulkVelMatrix,N,BulkCoeff);
	
	  // double BulkStabCoeff=BulkCoeff*Tau*Density/TimeStep;
	  // this->ComputeBulkMatrixForPressureAcc(BulkAccMatrix,N,BulkStabCoeff);

	  double BoundLHSCoeff=Tau*4.0*GaussWeight/(ElemSize*ElemSize);
 	  if(TDim==3){
	    BoundLHSCoeff=Tau*2*GaussWeight/(0.81649658*ElemSize*ElemSize);
	  }

	  this->ComputeBoundLHSMatrix(rLeftHandSideMatrix,N,BoundLHSCoeff);

	  double NProjSpatialDefRate=this->CalcNormalProjectionDefRate(rElementalVariables.SpatialDefRate);

	  double BoundRHSCoeffAcc=Tau*Density*2*GaussWeight/ElemSize;
	  double BoundRHSCoeffDev=Tau*8.0*NProjSpatialDefRate*DeviatoricCoeff*GaussWeight/(ElemSize*ElemSize);
	  if(TDim==3){
	    BoundRHSCoeffAcc=Tau*GaussWeight*Density/(0.81649658*ElemSize);
	    BoundRHSCoeffDev=Tau*GaussWeight*4.0*NProjSpatialDefRate*DeviatoricCoeff/(0.81649658*ElemSize*ElemSize);
	  }

	  this->ComputeBoundRHSVector(rRightHandSideVector,N,TimeStep,BoundRHSCoeffAcc,BoundRHSCoeffDev);

	  double StabLaplacianWeight=Tau*GaussWeight;
	  this->ComputeStabLaplacianMatrix(rLeftHandSideMatrix,rDN_DX,StabLaplacianWeight);

	  for (SizeType i = 0; i < NumNodes; ++i)
	    {         
	      // RHS contribution
	      // Velocity divergence   
	      rRightHandSideVector[i] += GaussWeight * N[i] * rElementalVariables.VolumetricDefRate;
	      this->AddStabilizationNodalTermsRHS(rRightHandSideVector,Tau,Density,GaussWeight,rDN_DX,i);
	    }
	}

      }   
   

    if(computeElement==true){

      VectorType PressureValues= ZeroVector(NumNodes);
      VectorType PressureValuesForRHS= ZeroVector(NumNodes);
      this->GetPressureValues(PressureValuesForRHS,0);
      //the LHS matrix up to now just contains the laplacian term and the bound term
      noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,PressureValuesForRHS);

      this->GetPressureValues(PressureValues,1);
      noalias(PressureValuesForRHS)+=-PressureValues;
      MatrixType BulkMatrix = ZeroMatrix(NumNodes,NumNodes);
      double lumpedBulkCoeff =totalVolume/(VolumetricCoeff);
      double lumpedBulkStabCoeff=lumpedBulkCoeff*Tau*Density/TimeStep;
      this->ComputeBulkMatrixForPressureVelLump(BulkMatrix,lumpedBulkCoeff);
      noalias(rLeftHandSideMatrix)+=BulkMatrix;
      noalias(rRightHandSideVector) -= prod(BulkMatrix,PressureValuesForRHS);
      // rRightHandSideVector -= prod(BulkVelMatrix,PressureValuesForRHS);

      this->GetPressureVelocityValues(PressureValues,0);
      noalias(PressureValuesForRHS)+=-PressureValues*TimeStep;
      noalias(BulkMatrix) = ZeroMatrix(NumNodes,NumNodes);
      this->ComputeBulkMatrixForPressureAccLump(BulkMatrix,lumpedBulkStabCoeff);
      noalias(rLeftHandSideMatrix)+=BulkMatrix;
      noalias(rRightHandSideVector) -=prod(BulkMatrix,PressureValuesForRHS);
      // rRightHandSideVector -=prod(BulkAccMatrix,PressureVelocityVariation);

    }else{
      double lumpedBulkCoeff =totalVolume*Tau*Density/(TimeStep*VolumetricCoeff);
      MatrixType BulkVelMatrixLump = ZeroMatrix(NumNodes,NumNodes);
      this->ComputeBulkMatrixForPressureVelLump(BulkVelMatrixLump,lumpedBulkCoeff);
      noalias(rLeftHandSideMatrix)+=BulkVelMatrixLump;
      VectorType PressureValues= ZeroVector(NumNodes);
      VectorType PressureValuesForRHS= ZeroVector(NumNodes);
      this->GetPressureValues(PressureValuesForRHS,0);
      this->GetPressureValues(PressureValues,1);
      noalias(PressureValuesForRHS)+=-PressureValues;
      noalias(rRightHandSideVector) -= prod(BulkVelMatrixLump,PressureValuesForRHS);
    }
 


  }
  

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPFluidElement<TDim>::GetPressureVelocityValues(Vector& rValues,
									       const int Step)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    if (rValues.size() != NumNodes) rValues.resize(NumNodes);

    for (SizeType i = 0; i < NumNodes; ++i){
      rValues[i] = rGeom[i].FastGetSolutionStepValue(PRESSURE_VELOCITY,Step);

    }
  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPFluidElement<TDim>::GetPressureAccelerationValues(Vector& rValues,
										   const int Step)
  {
    GeometryType& rGeom = this->GetGeometry();
    const SizeType NumNodes = rGeom.PointsNumber();

    if (rValues.size() != NumNodes) rValues.resize(NumNodes);

    for (SizeType i = 0; i < NumNodes; ++i){
      rValues[i] = rGeom[i].FastGetSolutionStepValue(PRESSURE_ACCELERATION,Step);

    }
  }




  template class TwoStepUpdatedLagrangianVPFluidElement<2>;
  template class TwoStepUpdatedLagrangianVPFluidElement<3>;

}
