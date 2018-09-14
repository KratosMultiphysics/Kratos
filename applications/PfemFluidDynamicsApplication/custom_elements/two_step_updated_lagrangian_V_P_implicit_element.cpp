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
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_element.h"
#include "includes/cfd_variables.h"

namespace Kratos {


  template< unsigned int TDim >
  Element::Pointer TwoStepUpdatedLagrangianVPImplicitElement<TDim>::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {
    KRATOS_TRY;

    TwoStepUpdatedLagrangianVPImplicitElement NewElement(NewId, this->GetGeometry().Create( rThisNodes ), this->pGetProperties() );
    return Element::Pointer( new TwoStepUpdatedLagrangianVPImplicitElement(NewElement) );

    KRATOS_CATCH("");

  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitElement<TDim>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
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
  void TwoStepUpdatedLagrangianVPImplicitElement<TDim>::CalculateLocalMomentumEquations(MatrixType& rLeftHandSideMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
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
    const unsigned int NumGauss = GaussWeights.size();
    const double TimeStep=rCurrentProcessInfo[DELTA_TIME];

    double theta=this->GetThetaMomentum();

    ElementalVariables rElementalVariables;
    this->InitializeElementalVariables(rElementalVariables);

    double totalVolume=0;
    double MeanValueMass=0;
    double Density=0.0;
    double DeviatoricCoeff = 0;
    double VolumetricCoeff = 0;
    // this->ComputeMaterialParameters(Density,DeviatoricCoeff,VolumetricCoeff,TimeStep);

    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; g++)
      {
	const double GaussWeight = GaussWeights[g];
	totalVolume+=GaussWeight;
	const ShapeFunctionsType& N = row(NContainer,g);
	const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];

	double Pressure=0;
	double OldPressure=0;

	this->EvaluateInPoint(Pressure,PRESSURE,N,0);

	this->EvaluateInPoint(OldPressure,PRESSURE,N,1);

	rElementalVariables.MeanPressure=OldPressure*(1-theta)+Pressure*theta;

	bool computeElement=this->CalcMechanicsUpdated(rElementalVariables,rCurrentProcessInfo,rDN_DX,g);

	this->ComputeMaterialParameters(Density,DeviatoricCoeff,VolumetricCoeff,rCurrentProcessInfo,rElementalVariables);

	this->CalcElasticPlasticCauchySplitted(rElementalVariables,TimeStep,g);

	// std::vector<double> rOutput;
	// this->GetElementalValueForOutput(YIELDED,rOutput);

	if(computeElement==true){
	  // Add integration point contribution to the local mass matrix
	  // double DynamicWeight=GaussWeight*Density;
	  // this->ComputeMassMatrix(MassMatrix,N,DynamicWeight,MeanValueMass);

	  this->AddExternalForces(rRightHandSideVector,Density,N,GaussWeight);

	  this->AddInternalForces(rRightHandSideVector,rDN_DX,rElementalVariables,GaussWeight);

	  // double lumpedDynamicWeight=GaussWeight*Density;
	  // this->ComputeLumpedMassMatrix(MassMatrix,lumpedDynamicWeight,MeanValueMass);

	  // double MeanValueMaterial=0.0;
	  // this->ComputeMeanValueMaterialTangentMatrix(rElementalVariables,MeanValueMaterial,rDN_DX,DeviatoricCoeff,VolumetricCoeff,GaussWeight,MeanValueMass,TimeStep);
	  // double deviatoricCoeffTemp=DeviatoricCoeff;
	  // DeviatoricCoeff=0;
	  // // Add viscous term
	  // this->ComputeCompleteTangentTerm(rElementalVariables,rLeftHandSideMatrix,rDN_DX,DeviatoricCoeff,VolumetricCoeff,theta,GaussWeight);

	  // double staticFrictionCoefficient=0.34;
	  // double dynamicFrictionCoefficient=0.6;
	  // double inertialNumberZero=0.279;
	  // double grainDiameter=0.001;
	  // double grainDensity=2500;
	  // double meanPressure=fabs(rElementalVariables.MeanPressure);
	  // double deltaFrictionCoefficient=dynamicFrictionCoefficient-staticFrictionCoefficient;
	  // double smallScale=grainDiameter/sqrt(meanPressure/grainDensity);
	  // double denominator=inertialNumberZero+rElementalVariables.EquivalentStrainRate*smallScale;
	  // if(meanPressure!=0 && rElementalVariables.EquivalentStrainRate!=0){
	  //   DeviatoricCoeff=inertialNumberZero*deltaFrictionCoefficient*meanPressure*smallScale/pow(denominator,2);
	  // }
	  this->ComputeCompleteTangentTerm(rElementalVariables,StiffnessMatrix,rDN_DX,DeviatoricCoeff,VolumetricCoeff,theta,GaussWeight);
	  // DeviatoricCoeff=deviatoricCoeffTemp;
	}
      }

    double lumpedDynamicWeight=totalVolume*Density;
    this->ComputeLumpedMassMatrix(MassMatrix,lumpedDynamicWeight,MeanValueMass);

    double BulkReductionCoefficient=1.0;
    double MeanValueStiffness=0.0;
    this->ComputeBulkReductionCoefficient(MassMatrix,StiffnessMatrix,MeanValueStiffness,BulkReductionCoefficient,TimeStep);
    if(BulkReductionCoefficient!=1.0){
      // VolumetricCoeff*=BulkReductionCoefficient;
      VolumetricCoeff*=MeanValueMass*2.0/(TimeStep*MeanValueStiffness);
      StiffnessMatrix= ZeroMatrix(LocalSize,LocalSize);

      for (unsigned int g = 0; g < NumGauss; g++)
    	{
    	  const double GaussWeight = GaussWeights[g];
    	  const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];
    	  this->ComputeCompleteTangentTerm(rElementalVariables,StiffnessMatrix,rDN_DX,DeviatoricCoeff,VolumetricCoeff,theta,GaussWeight);
    	}
    }


    // Add residual of previous iteration to RHS
    VectorType VelocityValues = ZeroVector(LocalSize);
    VectorType AccelerationValues = ZeroVector(LocalSize);

    // //1st order
    // this->GetVelocityValues(VelocityValues,0);
    // AccelerationValues = VelocityValues/TimeStep;
    // this->GetAccelerationValues(LastAccValues,0);
    // this->GetVelocityValues(VelocityValues,1);
    // AccelerationValues += -VelocityValues/TimeStep;
    // noalias( rRightHandSideVector ) += -prod(MassMatrix,AccelerationValues);
    // noalias( rLeftHandSideMatrix ) +=  MassMatrix/TimeStep;

    //2nd order
    this->GetAccelerationValues(AccelerationValues,0);
    this->GetVelocityValues(VelocityValues,0);
    noalias(AccelerationValues)+=-2.0*VelocityValues/TimeStep;
    this->GetVelocityValues(VelocityValues,1);
    noalias(AccelerationValues)+=2.0*VelocityValues/TimeStep;//these are negative accelerations
    noalias( rRightHandSideVector )+= prod(MassMatrix,AccelerationValues);
    noalias( rLeftHandSideMatrix ) +=  StiffnessMatrix + MassMatrix*2/TimeStep;


    // // Add residual of previous iteration to RHS
    // VectorType VelocityValues = ZeroVector(LocalSize);
    // VectorType UpdatedAccelerations = ZeroVector(LocalSize);
    // VectorType LastAccValues = ZeroVector(LocalSize);

    // // //1st order
    // // this->GetVelocityValues(VelocityValues,0);
    // // UpdatedAccelerations = VelocityValues/TimeStep;
    // // this->GetAccelerationValues(LastAccValues,0);
    // // this->GetVelocityValues(VelocityValues,1);
    // // UpdatedAccelerations += -VelocityValues/TimeStep;
    // // // UpdatedAccelerations =LastAccValues;
    // // noalias( rRightHandSideVector ) += -prod(MassMatrix,UpdatedAccelerations);
    // // noalias( rLeftHandSideMatrix ) +=  MassMatrix/TimeStep;

    // //2nd order
    // this->GetVelocityValues(VelocityValues,0);
    // UpdatedAccelerations = 2.0*VelocityValues/TimeStep;
    // this->GetAccelerationValues(LastAccValues,0);
    // this->GetVelocityValues(VelocityValues,1);
    // UpdatedAccelerations += -2.0*VelocityValues/TimeStep - LastAccValues;
    // noalias( rRightHandSideVector ) += -prod(MassMatrix,UpdatedAccelerations);
    // noalias( rLeftHandSideMatrix ) +=  StiffnessMatrix;
    // noalias( rLeftHandSideMatrix ) +=  MassMatrix*2/TimeStep;

    KRATOS_CATCH( "" );

  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPImplicitElement<TDim>::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
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
  void TwoStepUpdatedLagrangianVPImplicitElement<2>::ComputeCompleteTangentTerm(ElementalVariables & rElementalVariables,
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
	    double lagDNXi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,0);
	    double lagDNYi=rDN_DX(i,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(i,1)*rElementalVariables.InvFgrad(1,1);
	    double lagDNXj=rDN_DX(j,0)*rElementalVariables.InvFgrad(0,0)+rDN_DX(j,1)*rElementalVariables.InvFgrad(1,0);
	    double lagDNYj=rDN_DX(j,0)*rElementalVariables.InvFgrad(0,1)+rDN_DX(j,1)*rElementalVariables.InvFgrad(1,1);
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
  void TwoStepUpdatedLagrangianVPImplicitElement<3>::ComputeCompleteTangentTerm(ElementalVariables & rElementalVariables,
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
  int TwoStepUpdatedLagrangianVPImplicitElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo)
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




  /*
   * Template class definition (this should allow us to compile the desired template instantiations)
   */

  template class TwoStepUpdatedLagrangianVPImplicitElement<2>;
  template class TwoStepUpdatedLagrangianVPImplicitElement<3>;

}
