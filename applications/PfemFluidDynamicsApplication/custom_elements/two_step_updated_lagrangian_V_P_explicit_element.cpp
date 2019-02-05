//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:              April 2018 $
//   Revision:            $Revision:                 0.0 $
//
//   Implementation of the Gauss-Seidel two step Updated Lagrangian Velocity-Pressure element
//     ( There is a ScalingConstant to multiply the mass balance equation for a number because i read it somewhere)
//

// System includes

// External includes

// Project includes
#include "custom_elements/two_step_updated_lagrangian_V_P_explicit_element.h"
#include "includes/cfd_variables.h"

namespace Kratos {

  /*
   * public TwoStepUpdatedLagrangianVPExplicitElement<TDim> functions
   */


  // template< unsigned int TDim >
  // TwoStepUpdatedLagrangianVPExplicitElement<TDim>::TwoStepUpdatedLagrangianVPExplicitElement(TwoStepUpdatedLagrangianVPExplicitElement  const& rOther)
  // :Element(rOther)
  // {
  //   KRATOS_TRY;
  //   KRATOS_CATCH("");
  // }


  template< unsigned int TDim >
  Element::Pointer TwoStepUpdatedLagrangianVPExplicitElement<TDim>::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
  {
    KRATOS_TRY;

    TwoStepUpdatedLagrangianVPExplicitElement NewElement(NewId, this->GetGeometry().Create( rThisNodes ), this->pGetProperties() );
    return Element::Pointer( new TwoStepUpdatedLagrangianVPExplicitElement(NewElement) );

    KRATOS_CATCH("");

  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPExplicitElement<TDim>::Initialize()
  {
    KRATOS_TRY;
    // std::cout<<"INITIALIZE  !!!"<<std::endl;
    // GeometryType& rGeom = this->GetGeometry();
    // const SizeType NumNodes = rGeom.PointsNumber();
    // for (SizeType i = 0; i < NumNodes; ++i){
    //   rGeom[i].FastGetSolutionStepValue(PRESSURE)=0;
    // }
    KRATOS_CATCH( "" );
  }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPExplicitElement<TDim>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
									     VectorType& rRightHandSideVector,
									     ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY;

    // this->CalculateExplicitContinuityEquation(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
    this->CalculateMassMatrixMomentum(rLeftHandSideMatrix,rCurrentProcessInfo);
    this->CalculateRightHandSideMomentum(rRightHandSideVector,rCurrentProcessInfo);
    KRATOS_CATCH("");
  }

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPExplicitElement<TDim>::GetValueOnIntegrationPoints( const Variable<double>& rVariable,
										     std::vector<double>& rValues,
										     const ProcessInfo& rCurrentProcessInfo )
  {
    if ( rVariable == YIELDED)
      {
	rValues[0]=this->GetValue(YIELDED);
      }
  }



  template< unsigned int TDim >
  int TwoStepUpdatedLagrangianVPExplicitElement<TDim>::Check(const ProcessInfo &rCurrentProcessInfo)
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
    if(DENSITY.Key() == 0)
      KRATOS_THROW_ERROR(std::invalid_argument,"DENSITY Key is 0. Check that the application was correctly registered.","");
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
	if(this->GetGeometry()[i].SolutionStepsDataHas(DENSITY) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing DENSITY variable on solution step data for node ",this->GetGeometry()[i].Id());
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
	if(this->GetGeometry()[i].HasDofFor(DENSITY) == false)
	  KRATOS_THROW_ERROR(std::invalid_argument,"missing DENSITY component degree of freedom on node ",this->GetGeometry()[i].Id());
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




  // template< unsigned int TDim>
  // bool TwoStepUpdatedLagrangianVPExplicitElement<TDim>::CalcMechanicsUpdated(ElementalVariables & rElementalVariables,
  // 									     const ProcessInfo& rCurrentProcessInfo,
  // 									     const ShapeFunctionDerivativesType& rDN_DX,
  // 									     unsigned int g)
  // {

  //   double theta=this->GetThetaMomentum();
  //   // bool computeElement=this->CalcStrainRate(rElementalVariables,rCurrentProcessInfo,rDN_DX,theta);
  //   bool computeElement=this->CalcCompleteStrainRate(rElementalVariables,rCurrentProcessInfo,rDN_DX,theta);
  //   return computeElement;

  // }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPExplicitElement<TDim>::CalculateMassMatrix(Matrix& rMassMatrix,
									    ProcessInfo& rCurrentProcessInfo)
  {

    KRATOS_TRY

      // std::cout<<"Calculate Mass Matrix for explicit PFEM"<<std::endl;

      GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    if(rMassMatrix.size1() != NumNodes )
      rMassMatrix.resize(NumNodes,NumNodes,false);

    noalias(rMassMatrix)= ZeroMatrix(NumNodes,NumNodes);

    // ShapeFunctionDerivativesArrayType DN_DX;
    // Matrix NContainer;
    // VectorType GaussWeights;
    // this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);
    VectorType GaussWeights;
    this->CalculateGeometryData(GaussWeights);
    const unsigned int NumGauss = GaussWeights.size();
    double totalVolume=0;

    // double VolumetricCoeff = this->mMaterialDeviatoricCoefficient;

    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; g++)
      {
	const double GaussWeight = GaussWeights[g];
	totalVolume+=GaussWeight;
      }

    // double timeStep=rCurrentProcessInfo[DELTA_TIME];
    MatrixType BulkMatrix = ZeroMatrix(NumNodes,NumNodes);
    // double lumpedBulkCoeff =totalVolume/(VolumetricCoeff);
    double lumpedBulkCoeff =totalVolume;
    this->ComputeBulkMatrixLump(BulkMatrix,lumpedBulkCoeff);
    // this->ComputeBulkMatrixConsistent(BulkMatrix,lumpedBulkCoeff);
    rMassMatrix+=BulkMatrix;

    KRATOS_CATCH( "" )

      }




  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPExplicitElement<TDim>::CalculateMassMatrixMomentum(Matrix& rMassMatrix,
										    ProcessInfo& rCurrentProcessInfo)
  {

    KRATOS_TRY

      // std::cout<<"Calculate Mass Matrix for explicit PFEM"<<std::endl;

      GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = TDim * NumNodes;

    if(rMassMatrix.size1() != LocalSize )
      rMassMatrix.resize(LocalSize,LocalSize,false);

    noalias(rMassMatrix)= ZeroMatrix(LocalSize,LocalSize);

    // ShapeFunctionDerivativesArrayType DN_DX;
    // Matrix NContainer;
    VectorType GaussWeights;
    // this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);
    this->CalculateGeometryData(GaussWeights);
    const unsigned int NumGauss = GaussWeights.size();
    double totalVolume=0;

    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; g++)
      {
	const double GaussWeight = GaussWeights[g];
	totalVolume+=GaussWeight;
      }

    //the time step is not included, the rhs is defined consequently
    double  Weight=totalVolume*this->GetProperties()[DENSITY];
    //lumped
    if((NumNodes==3 && TDim==2) || (NumNodes==4 && TDim==3)){
      double Coeff=1.0+TDim;
      for (SizeType i = 0; i < NumNodes; ++i)
	{
	  double Mij = Weight/Coeff;
	  for ( unsigned int j = 0; j <  TDim; j++ )
	    {
	      unsigned int index = i * TDim + j;
	      // std::cout<<"index "<<index<<"  massTerm "<<Mij<<std::endl;
	      rMassMatrix( index, index ) += Mij;
	    }

	}
    }
    else if(NumNodes==6 && TDim==2){
      double Mij = Weight/57.0;
      double consistent=1.0;
      for (SizeType i = 0; i < NumNodes; ++i)
	{
	  if(i<3){
	    consistent=3.0;
	  }else{
	    consistent=16.0;
	  }
	  for ( unsigned int j = 0; j <  TDim; j++ )
	    {
	      unsigned int index = i * TDim + j;
	      rMassMatrix( index, index ) += Mij*consistent;
	    }

	}

    }else{
      std::cout<<"ComputeLumpedMassMatrix 3D quadratic not yet implemented!"<<std::endl;
    }


    KRATOS_CATCH( "" )

      }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPExplicitElement<TDim>::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

      // it just compute the external and internal forces terms. The inertial forces are already taken into account with the temporal scheme

      GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    MatrixType MassMatrix= ZeroMatrix(NumNodes,NumNodes);

    if( rRightHandSideVector.size() != NumNodes )
      rRightHandSideVector.resize(NumNodes);

    rRightHandSideVector = ZeroVector(NumNodes);

    // Shape functions and integration points
    // ShapeFunctionDerivativesArrayType DN_DX;
    // Matrix NContainer;
    VectorType GaussWeights;
    // this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);
    this->CalculateGeometryData(GaussWeights);
    const unsigned int NumGauss = GaussWeights.size();

    // double theta=this->GetThetaContinuity();

    ElementalVariables rElementalVariables;
    this->InitializeElementalVariables(rElementalVariables);

    // double VolumetricCoeff = this->mMaterialDeviatoricCoefficient;
    // double Density = this->mMaterialDensity;
    // double DeviatoricCoeff = this->mMaterialDeviatoricCoefficient;
    // double TimeStep=rCurrentProcessInfo[DELTA_TIME];
    // double Tau=0;
    // double ElemSize = this->ElementSize();

    double totalVolume=0;

    // MatrixType BulkMatrixConsistent = ZeroMatrix(NumNodes,NumNodes);

    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; g++)
      {
	const double GaussWeight = GaussWeights[g];
	totalVolume+=GaussWeight;
	// const ShapeFunctionsType& N = row(NContainer,g);
	// const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];
	// bool computeElement=this->CalcCompleteStrainRate(rElementalVariables,rCurrentProcessInfo,rDN_DX,theta);
      }

    MatrixType BulkMatrix = ZeroMatrix(NumNodes,NumNodes);
    // double lumpedBulkCoeff =totalVolume/(VolumetricCoeff);
    double lumpedBulkCoeff =totalVolume;
    // VectorType UpdatedPressure = ZeroVector(NumNodes);
    VectorType UpdatedDensity = ZeroVector(NumNodes);

    // this->GetPressureValues(UpdatedPressure,0);
    this->GetDensityValues(UpdatedDensity,0);
    this->ComputeBulkMatrixRHS(BulkMatrix,lumpedBulkCoeff);
    rRightHandSideVector += prod(BulkMatrix,UpdatedDensity);


    KRATOS_CATCH( "" )
      }



  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPExplicitElement<TDim>::CalculateRightHandSideMomentum( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
  {
    KRATOS_TRY

      // it just compute the external and internal forces terms. The inertial forces are already taken into account with the temporal scheme

      GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = TDim * NumNodes;

    MatrixType MassMatrix= ZeroMatrix(LocalSize,LocalSize);

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

    // double theta=this->GetThetaMomentum();

    ElementalVariables rElementalVariables;
    this->InitializeElementalVariables(rElementalVariables);

    double Density=0.0;
    double DeviatoricCoeff = 0;
    double VolumetricCoeff = 0;

    // Loop on integration points
    for (unsigned int g = 0; g < NumGauss; g++)
      {
	const double GaussWeight = GaussWeights[g];
	const ShapeFunctionsType& N = row(NContainer,g);
	const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];

	double Pressure=0;
	// double OldPressure=0;
	this->EvaluateInPoint(Pressure,PRESSURE,N,0);
	// this->EvaluateInPoint(OldPressure,PRESSURE,N,1);
	// rElementalVariables.MeanPressure=OldPressure*(1-theta)+Pressure*theta;
	rElementalVariables.MeanPressure=Pressure;

	bool computeElement=this->CalcMechanicsUpdated(rElementalVariables,rCurrentProcessInfo,rDN_DX,g);

	this->ComputeMaterialParameters(Density,DeviatoricCoeff,VolumetricCoeff,rCurrentProcessInfo,rElementalVariables);

	this->CalcElasticPlasticCauchySplitted(rElementalVariables,TimeStep,g);

	if(computeElement==true){
	  // Add integration point contribution to the local mass matrix
	  this->AddExternalForces(rRightHandSideVector,Density,N,GaussWeight);
	  this->AddInternalForces(rRightHandSideVector,rDN_DX,rElementalVariables,GaussWeight);
	}
      }

    KRATOS_CATCH( "" )
      }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPExplicitElement<TDim>::AddExplicitContribution(const VectorType& rRHSVector,
										const Variable<VectorType>& rRHSVariable,
										Variable<array_1d<double,3> >& rDestinationVariable,
										const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

      const unsigned int number_of_nodes = this->GetGeometry().PointsNumber();
    const unsigned int dimension       = this->GetGeometry().WorkingSpaceDimension();

    if( rRHSVariable == EXTERNAL_FORCES_VECTOR && rDestinationVariable == EXTERNAL_FORCE )
      {

	for(unsigned int i=0; i< number_of_nodes; i++)
	  {
	    int index = dimension * i;

	    this->GetGeometry()[i].SetLock();

	    array_1d<double, 3 > &ExternalForce = this->GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_FORCE);
	    for(unsigned int j=0; j<dimension; j++)
	      {
		ExternalForce[j] += rRHSVector[index + j];
	      }

	    this->GetGeometry()[i].UnSetLock();
	  }
      }

    if( rRHSVariable == INTERNAL_FORCES_VECTOR && rDestinationVariable == INTERNAL_FORCE )
      {

	for(unsigned int i=0; i< number_of_nodes; i++)
	  {
	    int index = dimension * i;

	    this->GetGeometry()[i].SetLock();

	    array_1d<double, 3 > &InternalForce = this->GetGeometry()[i].FastGetSolutionStepValue(INTERNAL_FORCE);
	    for(unsigned int j=0; j<dimension; j++)
	      {
		InternalForce[j] += rRHSVector[index + j];
	      }

	    this->GetGeometry()[i].UnSetLock();
	  }
      }


    if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL )
      {

	for(unsigned int i=0; i< number_of_nodes; i++)
	  {
	    int index = dimension * i;

	    this->GetGeometry()[i].SetLock();

	    array_1d<double, 3 > &ForceResidual = this->GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);

	    for(unsigned int j=0; j<dimension; j++)
	      {
		ForceResidual[j] += rRHSVector[index + j];
	      }

	    this->GetGeometry()[i].UnSetLock();
	  }
      }

    KRATOS_CATCH( "" )
      }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPExplicitElement<TDim>::AddExplicitContribution(const VectorType& rRHSVector,
										const Variable<VectorType>& rRHSVariable,
										Variable<double >& rDestinationVariable,
										const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

      const unsigned int number_of_nodes = this->GetGeometry().PointsNumber();

    if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == NODAL_ERROR )
      {

	for(unsigned int i=0; i< number_of_nodes; i++)
	  {
	    this->GetGeometry()[i].SetLock();

	    double &ForceResidual = this->GetGeometry()[i].FastGetSolutionStepValue(NODAL_ERROR);
	    ForceResidual += rRHSVector[i];

	    this->GetGeometry()[i].UnSetLock();
	  }
      }

    KRATOS_CATCH( "" )
      }


  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPExplicitElement<TDim>::CalculateExplicitContinuityEquation(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  {
    GeometryType& rGeom = this->GetGeometry();
    const unsigned int NumNodes = rGeom.PointsNumber();

    // Check sizes and initialize
    if( rLeftHandSideMatrix.size1() != NumNodes )
      rLeftHandSideMatrix.resize(NumNodes,NumNodes,false);

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

    // const double TimeStep=rCurrentProcessInfo[DELTA_TIME];
    double theta=this->GetThetaContinuity();

    ElementalVariables rElementalVariables;
    this->InitializeElementalVariables(rElementalVariables);

    double VolumetricCoeff = this->mMaterialVolumetricCoefficient;

    double totalVolume=0;

    // MatrixType BulkMatrixConsistent = ZeroMatrix(NumNodes,NumNodes);

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
	  // this->ComputeBulkMatrix(BulkMatrixConsistent,N,BulkCoeff);

	  // // if(this->Is(FLUID)){
	  // this->CalculateTauFIC(Tau,ElemSize,Density,DeviatoricCoeff,rCurrentProcessInfo);
	  // double BoundLHSCoeff=Tau*4.0*GaussWeight/(ElemSize*ElemSize);
	  // this->ComputeBoundLHSMatrix(rLeftHandSideMatrix,N,BoundLHSCoeff);
	  // double NProjSpatialDefRate=this->CalcNormalProjectionDefRate(rElementalVariables.SpatialDefRate);
	  // double BoundRHSCoeffAcc=Tau*Density*2*GaussWeight/ElemSize;
	  // double BoundRHSCoeffDev=Tau*8.0*NProjSpatialDefRate*DeviatoricCoeff*GaussWeight/(ElemSize*ElemSize);
	  // this->ComputeBoundRHSVector(rRightHandSideVector,N,TimeStep,BoundRHSCoeffAcc,BoundRHSCoeffDev);
	  // double StabLaplacianWeight=Tau*GaussWeight;
	  // this->ComputeStabLaplacianMatrix(rLeftHandSideMatrix,rDN_DX,StabLaplacianWeight);
	  // // }

	  for (SizeType i = 0; i < NumNodes; ++i)
	    {
	      // RHS contribution
	      // Velocity divergence
	      rRightHandSideVector[i] += GaussWeight * N[i] * rElementalVariables.VolumetricDefRate;
	      // // if(this->Is(FLUID)){
	      // this->AddStabilizationNodalTermsRHS(rRightHandSideVector,Tau,Density,GaussWeight,rDN_DX,i);
	      // }
	    }

	}

      }

    // double timeStep=rCurrentProcessInfo[DELTA_TIME];
    MatrixType BulkMatrix = ZeroMatrix(NumNodes,NumNodes);
    double lumpedBulkCoeff =totalVolume/(VolumetricCoeff);
    this->ComputeBulkMatrixLump(BulkMatrix,lumpedBulkCoeff);
    // this->ComputeBulkMatrixConsistent(BulkMatrix,lumpedBulkCoeff);
    rLeftHandSideMatrix+=BulkMatrix;

    VectorType UpdatedPressure = ZeroVector(NumNodes);
    this->GetPressureValues(UpdatedPressure,0);
    BulkMatrix = ZeroMatrix(NumNodes,NumNodes);
    this->ComputeBulkMatrixRHS(BulkMatrix,lumpedBulkCoeff);
    rRightHandSideVector += prod(BulkMatrix,UpdatedPressure);
    // rRightHandSideVector += prod(BulkMatrixConsistent,UpdatedPressure);
  }


  /*
   * Template class definition (this should allow us to compile the desired template instantiations)
   */

  template class TwoStepUpdatedLagrangianVPExplicitElement<2>;
  template class TwoStepUpdatedLagrangianVPExplicitElement<3>;

}
