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
  void TwoStepUpdatedLagrangianVPSolidElement<TDim>::Initialize()
  {
    KRATOS_TRY;

    // LargeDisplacementElement::Initialize();
    const GeometryType& rGeom = this->GetGeometry();
    SizeType integration_points_number = rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_1);
    const unsigned int dimension       = rGeom.WorkingSpaceDimension();
 
    //Resize historic deformation gradient
    if ( this->mOldFgrad.size() != integration_points_number )
      this->mOldFgrad.resize( integration_points_number ); 
 
    if ( this->mDetFgrad.size() != integration_points_number )
      this->mDetFgrad.resize( integration_points_number, false );

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
        this->mDetFgrad[PointNumber] = 1;
        this->mOldFgrad[PointNumber] = identity_matrix<double> (dimension);
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



  // template< unsigned int TDim >
  // void TwoStepUpdatedLagrangianVPSolidElement<TDim>::CalculateLocalMomentumEquations(MatrixType& rLeftHandSideMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
  // {

  //   const GeometryType& rGeom = this->GetGeometry();
  //   const SizeType NumNodes = rGeom.PointsNumber();
  //   const SizeType LocalSize = TDim * NumNodes;

  //   // Check sizes and initialize
  //   if( rLeftHandSideMatrix.size1() != LocalSize )
  //     rLeftHandSideMatrix.resize(LocalSize,LocalSize);

  //   rLeftHandSideMatrix = ZeroMatrix(LocalSize,LocalSize);

  //   if( rRightHandSideVector.size() != LocalSize )
  //     rRightHandSideVector.resize(LocalSize);

  //   rRightHandSideVector = ZeroVector(LocalSize);

  //   // Shape functions and integration points
  //   ShapeFunctionDerivativesArrayType DN_DX;
  //   Matrix NContainer;
  //   VectorType GaussWeights;
  //   this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);
  //   const unsigned int NumGauss = GaussWeights.size();


  //   MatrixType MassMatrix = ZeroMatrix(LocalSize,LocalSize);
  //   MatrixType BulkMatrix = ZeroMatrix(LocalSize,LocalSize);

  //   const double TimeStep=0.5/rCurrentProcessInfo[BDF_COEFFICIENTS][2];

  //   // Loop on integration points
  //   for (unsigned int g = 0; g < NumGauss; g++)
  //     {
  // 	const double GaussWeight = GaussWeights[g];
  // 	const ShapeFunctionsType& N = row(NContainer,g);
  // 	const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];


  // 	ElementalVariables rElementalVariables;
  // 	this->InitializeElementalVariables(rElementalVariables,g);
	
  // 	double Pressure=0;
  // 	double OldPressure=0;
  // 	this->EvaluateInPoint(Pressure,PRESSURE,N,0);
  // 	this->EvaluateInPoint(OldPressure,PRESSURE,N,1);

  // 	rElementalVariables.MeanPressure=OldPressure*0.5+Pressure*0.5;

  // 	this->CalcMechanicsUpdated(rElementalVariables,rCurrentProcessInfo,g);
    
  //      	// Evaluate required variables at the integration point
  // 	double Density;
  // 	array_1d<double,3> BodyForce(3,0.0);

  // 	this->EvaluateInPoint(Density,DENSITY,N);
  // 	this->EvaluateInPoint(BodyForce,BODY_FORCE,N);

  // 	// Add integration point contribution to the local mass matrix
  // 	// double massWeight=GaussWeight*Density*2.0/TimeStep;
  // 	double dynamicWeight=GaussWeight*Density;

  // 	 this->ComputeLumpedMassMatrix(MassMatrix,dynamicWeight);
  // 	// this->ComputeMomentumMassTerm(MassMatrix,N,dynamicWeight);

  // 	this->AddExternalForces(rRightHandSideVector,Density,BodyForce,N,GaussWeight);

  // 	this->AddInternalForces(rRightHandSideVector,rDN_DX,rElementalVariables,GaussWeight);

  // 	// this->AddDynamicForces(rRightHandSideVector,dynamicWeight,TimeStep);

  // 	const double YoungModulus = 100000000.0;
  // 	const double PoissonRatio = 0;
  // 	double CurrFirstLame   = 0;
  // 	double CurrSecondLame  = 0;
  // 	double CurrBulkModulus = 0;

  // 	this->ComputeMaterialParameters(YoungModulus,PoissonRatio,CurrFirstLame,CurrSecondLame,CurrBulkModulus,TimeStep);
  // 	// Add viscous term
  //      	this->AddCompleteTangentTerm(rElementalVariables,rLeftHandSideMatrix,rDN_DX,CurrSecondLame,CurrBulkModulus,GaussWeight);

  //     }
    
  //   // Add residual of previous iteration to RHS
  //   VectorType VelocityValues = ZeroVector(LocalSize);
  //   this->GetVelocityValues(VelocityValues,0);
  //   VectorType UpdatedAccelerations = ZeroVector(LocalSize);
  //   UpdatedAccelerations = 2.0*VelocityValues/TimeStep;
  //   VectorType LastAccValues = ZeroVector(LocalSize);
  //   this->GetAccelerationValues(LastAccValues,0);
  //   this->GetVelocityValues(VelocityValues,1);
  //   UpdatedAccelerations += -2.0*VelocityValues/TimeStep - LastAccValues; 
  //   noalias( rRightHandSideVector ) += -prod(MassMatrix,UpdatedAccelerations);
  //   noalias( rLeftHandSideMatrix ) +=  MassMatrix*2/TimeStep;




  // }

  // template< unsigned int TDim >
  // void TwoStepUpdatedLagrangianVPSolidElement<TDim>::CalculateLocalContinuityEqForPressure(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
  // {


  //   GeometryType& rGeom = this->GetGeometry();
  //   const SizeType NumNodes = rGeom.PointsNumber();

  //   // Check sizes and initialize
  //   if( rLeftHandSideMatrix.size1() != NumNodes )
  //     rLeftHandSideMatrix.resize(NumNodes,NumNodes);

  //   rLeftHandSideMatrix = ZeroMatrix(NumNodes,NumNodes);

  //   if( rRightHandSideVector.size() != NumNodes )
  //     rRightHandSideVector.resize(NumNodes);

  //   rRightHandSideVector = ZeroVector(NumNodes);

  //   // Shape functions and integration points
  //   ShapeFunctionDerivativesArrayType DN_DX;
  //   Matrix NContainer;
  //   VectorType GaussWeights;
  //   this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);
  //   const unsigned int NumGauss = GaussWeights.size();

  //   // Stabilization parameters
  //   double ElemSize = this->ElementSize();
  //   double Tau=0;

  //   MatrixType BulkVelMatrix = ZeroMatrix(NumNodes,NumNodes);
  //   MatrixType BulkAccMatrix = ZeroMatrix(NumNodes,NumNodes);

  //   const double YoungModulus = 100000000.0;
  //   const double PoissonRatio = 0;
  //   double CurrFirstLame   = 0;
  //   double CurrSecondLame  = 0;
  //   double CurrBulkModulus = 0;
  //   const double TimeStep=0.5/rCurrentProcessInfo[BDF_COEFFICIENTS][2];

  //   this->ComputeMaterialParameters(YoungModulus,PoissonRatio,CurrFirstLame,CurrSecondLame,CurrBulkModulus,TimeStep);

  //   // Loop on integration points
  //   for (unsigned int g = 0; g < NumGauss; g++)
  //     {
  // 	// this->UpdateCauchyStress(g);
  //       const double GaussWeight = GaussWeights[g];
  //       const ShapeFunctionsType& N = row(NContainer,g);
  //       const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];

  //       // Evaluate required variables at the integration point
  //       double Density;
  // 	//        array_1d<double,3> Velocity(3,0.0);
  // 	//        array_1d<double,3> MeshVelocity(3,0.0);
  //       array_1d<double,3> BodyForce(3,0.0);
  //       array_1d<double,3> MomentumProjection(3,0.0);

  //       this->EvaluateInPoint(Density,DENSITY,N);
  // 	//        this->EvaluateInPoint(Velocity,VELOCITY,N);
  // 	//        this->EvaluateInPoint(MeshVelocity,MESH_VELOCITY,N);
  //       this->EvaluateInPoint(BodyForce,BODY_FORCE,N);
  // 	//        this->EvaluateInPoint(MomentumProjection,ADVPROJ,N);
  //       this->EvaluateInPoint(MomentumProjection,PRESS_PROJ,N);

  // 	// Evaluate the pressure and pressure gradient at this point (for the G * P_n term)
  //       array_1d<double,TDim> OldPressureGradient(TDim,0.0);
  //       this->EvaluateGradientInPoint(OldPressureGradient,PRESSURE,rDN_DX);

  //       // Stabilization parameters
  //       array_1d<double,3> ConvVel(3,0.0);
  //       // this->EvaluateConvVelocity(ConvVel,N);
  //       double Viscosity = 0;
  // 	this->EvaluateInPoint(Viscosity,VISCOSITY,N);

  // 	if(PoissonRatio>0.499){
  // 	  this->CalculateTauFIC(Tau,ElemSize,ConvVel,Density,Viscosity,rCurrentProcessInfo);
  // 	}


  // 	double BulkCoeff =GaussWeight/(CurrBulkModulus);
  // 	this->AddBulkMatrixForPressureVel(BulkVelMatrix,N,BulkCoeff);

  // 	// double Pressure=0;
  // 	// this->EvaluateInPoint(Pressure,PRESSURE,N,1);

  // 	// BulkCoeff*=Tau*Density/TimeStep;
  // 	// this->AddBulkMatrixForPressureAcc(BulkAccMatrix,N,BulkCoeff);
  // 	rLeftHandSideMatrix+=BulkVelMatrix;
  // 	// rLeftHandSideMatrix+=BulkAccMatrix;

  // 	// VectorType NewPressure = ZeroVector(NumNodes);
  // 	// VectorType OldPressure = ZeroVector(NumNodes);
  // 	// this->GetPressureValues(NewPressure,0);
  // 	// this->GetPressureValues(OldPressure,1);
  // 	// VectorType PressureAccRHSterm = ZeroVector(NumNodes);
  // 	// VectorType DeltaPressure = NewPressure-OldPressure;
  // 	// PressureAccRHSterm =prod(BulkAccMatrix,DeltaPressure);

  //       double DivU=0;
  //       this->EvaluateDivergenceInPoint(DivU,VELOCITY,rDN_DX);

  // 	// Add convection, stabilization and RHS contributions to the local system equation
  //       for (SizeType i = 0; i < NumNodes; ++i)
  // 	  {
  //           // // LHS contribution
  //           // for (SizeType j = 0; j < NumNodes; ++j)
  // 	    //   {
  //           //     double Lij = 0.0;
  //           //     for (SizeType d = 0; d < TDim; ++d){
  // 	    // 	  Lij += rDN_DX(i,d) * rDN_DX(j,d);
  // 	    // 	}
  // 	    // 	// Lij *= (LaplacianCoeff + Tau);
  //           //     Lij *= Tau;

  //           //     rLeftHandSideMatrix(i,j) += GaussWeight * Lij;
  // 	    //   }

  //           // RHS contribution
  // 	    // Velocity divergence
  //           double RHSi =  N[i] * DivU;

  //           for (SizeType d = 0; d < TDim; ++d)
  // 	      {
  //              RHSi += rDN_DX(i,d) * Tau * ( Density  * ( BodyForce[d] ) - OldPressureGradient[d]  );
  // 		if(d==2){
  // 		  RHSi += - rDN_DX(i,d) * Tau * ( Density * 9.81 );
  // 		}
  // 	      }

  //           rRightHandSideVector[i] += GaussWeight * RHSi;
  //           // rRightHandSideVector[i] += PressureAccRHSterm[i];
  // 	  }

  // 	VectorType vectorPressure = ZeroVector(NumNodes);
  // 	VectorType nodalPressure = ZeroVector(NumNodes);
  // 	VectorType nodalOldPressure = ZeroVector(NumNodes);
  // 	this->GetPressureValues(nodalPressure,0);
  // 	this->GetPressureValues(nodalOldPressure,1);
  // 	vectorPressure=nodalPressure-nodalOldPressure;
  // 	rRightHandSideVector -= prod(rLeftHandSideMatrix,vectorPressure);

  //     }
  // }
  

  template< unsigned int TDim >
  void TwoStepUpdatedLagrangianVPSolidElement<TDim>::ComputeMaterialParameters(double& DeviatoricCoeff,
									       double& VolumetricCoeff,
									       double timeStep)
  {

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

//   template<>
//   void TwoStepUpdatedLagrangianVPSolidElement<2>::VelocityEquationIdVector(EquationIdVectorType& rResult,
// 									   ProcessInfo& rCurrentProcessInfo)
//   {
//     GeometryType& rGeom = this->GetGeometry();
//     const SizeType NumNodes = rGeom.PointsNumber();
//     const SizeType LocalSize = NumNodes*2;

//     SizeType LocalIndex = 0;

//     if (rResult.size() != LocalSize)
//       rResult.resize(LocalSize, false);

//     const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);

//     for (SizeType i = 0; i < NumNodes; ++i)
//       {
//         rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X,xpos).EquationId();
//         rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
//       }
//   }

//   template<>
//   void TwoStepUpdatedLagrangianVPSolidElement<3>::VelocityEquationIdVector(EquationIdVectorType& rResult,
// 									   ProcessInfo& rCurrentProcessInfo)
//   {
//     GeometryType& rGeom = this->GetGeometry();
//     const SizeType NumNodes = rGeom.PointsNumber();
//     const SizeType LocalSize = 3*NumNodes;

//     SizeType LocalIndex = 0;

//     if (rResult.size() != LocalSize)
//       rResult.resize(LocalSize, false);

//     const unsigned int xpos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);

//     for (SizeType i = 0; i < NumNodes; ++i)
//       {
//         rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X,xpos).EquationId();
//         rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y,xpos+1).EquationId();
//         rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Z,xpos+2).EquationId();
//       }
//   }

//   template< unsigned int TDim >
//   void TwoStepUpdatedLagrangianVPSolidElement<TDim>::PressureEquationIdVector(EquationIdVectorType& rResult,
// 									      ProcessInfo& rCurrentProcessInfo)
//   {
//     GeometryType& rGeom = this->GetGeometry();
//     const SizeType NumNodes = rGeom.PointsNumber();

//     if (rResult.size() != NumNodes)
//       rResult.resize(NumNodes);

//     const unsigned int pos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);

//     for (SizeType i = 0; i < NumNodes; ++i)
//       rResult[i] = rGeom[i].GetDof(PRESSURE,pos).EquationId();
//   }

//   template<>
//   void TwoStepUpdatedLagrangianVPSolidElement<2>::GetVelocityDofList(DofsVectorType& rElementalDofList,
// 								     ProcessInfo& rCurrentProcessInfo)
//   {
//     GeometryType& rGeom = this->GetGeometry();
//     const SizeType NumNodes = rGeom.PointsNumber();
//     const SizeType LocalSize = 2*NumNodes;

//     if (rElementalDofList.size() != LocalSize)
//       rElementalDofList.resize(LocalSize);

//     SizeType LocalIndex = 0;

//     for (SizeType i = 0; i < NumNodes; ++i)
//       {
//         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
//         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
//       }
//   }

//   template<>
//   void TwoStepUpdatedLagrangianVPSolidElement<3>::GetVelocityDofList(DofsVectorType& rElementalDofList,
// 								     ProcessInfo& rCurrentProcessInfo)
//   {
//     GeometryType& rGeom = this->GetGeometry();
//     const SizeType NumNodes = rGeom.PointsNumber();
//     const SizeType LocalSize = 3*NumNodes;

//     if (rElementalDofList.size() != LocalSize)
//       rElementalDofList.resize(LocalSize);

//     SizeType LocalIndex = 0;

//     for (SizeType i = 0; i < NumNodes; ++i)
//       {
//         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
//         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
//         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Z);
//       }
//   }

//   template< unsigned int TDim >
//   void TwoStepUpdatedLagrangianVPSolidElement<TDim>::GetPressureDofList(DofsVectorType& rElementalDofList,
// 									ProcessInfo& rCurrentProcessInfo)
//   {
//     GeometryType& rGeom = this->GetGeometry();
//     const SizeType NumNodes = rGeom.PointsNumber();

//     if (rElementalDofList.size() != NumNodes)
//       rElementalDofList.resize(NumNodes);

//     SizeType LocalIndex = 0;

//     for (SizeType i = 0; i < NumNodes; ++i)
//       {
//         rElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PRESSURE);
//       }
//   }

//   template< unsigned int TDim >
//   void TwoStepUpdatedLagrangianVPSolidElement<TDim>::GetPressureValues(Vector& rValues,
// 								       const int Step)
//   {
//     GeometryType& rGeom = this->GetGeometry();
//     const SizeType NumNodes = rGeom.PointsNumber();

//     if (rValues.size() != NumNodes) rValues.resize(NumNodes);

//     for (SizeType i = 0; i < NumNodes; ++i)
//       rValues[i] = rGeom[i].FastGetSolutionStepValue(PRESSURE,Step);
//   }

  
//   template<>
//   void TwoStepUpdatedLagrangianVPSolidElement<2>::GetUpdatedPositions(Vector& rValues,
// 								      const ProcessInfo& rCurrentProcessInfo)
//   {
//     GeometryType& rGeom = this->GetGeometry();
//     const SizeType NumNodes = rGeom.PointsNumber();
//     const SizeType LocalSize = 2*NumNodes;
//     const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

//     if (rValues.size() != LocalSize) rValues.resize(LocalSize);

//     SizeType Index = 0;

//     for (SizeType i = 0; i < NumNodes; ++i)
//       {
//         rValues[Index++] = rGeom[i].X()-(rGeom[i].FastGetSolutionStepValue(VELOCITY_X,0)+rGeom[i].FastGetSolutionStepValue(VELOCITY_X,1))/BDFcoeffs[1];
// 	rValues[Index++] = rGeom[i].Y()-(rGeom[i].FastGetSolutionStepValue(VELOCITY_Y,0)+rGeom[i].FastGetSolutionStepValue(VELOCITY_Y,1))/BDFcoeffs[1];
//       }
//   }

//   template<>
//   void TwoStepUpdatedLagrangianVPSolidElement<3>::GetUpdatedPositions(Vector& rValues,
// 								      const ProcessInfo& rCurrentProcessInfo)
//  {
//     GeometryType& rGeom = this->GetGeometry();
//     const SizeType NumNodes = rGeom.PointsNumber();
//     const SizeType LocalSize = 3*NumNodes;
//     const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

//     if (rValues.size() != LocalSize) rValues.resize(LocalSize);

//     SizeType Index = 0;

//     for (SizeType i = 0; i < NumNodes; ++i)
//       {
//         rValues[Index++] = rGeom[i].X()-(rGeom[i].FastGetSolutionStepValue(VELOCITY_X,0)+rGeom[i].FastGetSolutionStepValue(VELOCITY_X,1))/BDFcoeffs[1];
// 	rValues[Index++] = rGeom[i].Y()-(rGeom[i].FastGetSolutionStepValue(VELOCITY_Y,0)+rGeom[i].FastGetSolutionStepValue(VELOCITY_Y,1))/BDFcoeffs[1];
// 	rValues[Index++] = rGeom[i].Z()-(rGeom[i].FastGetSolutionStepValue(VELOCITY_Z,0)+rGeom[i].FastGetSolutionStepValue(VELOCITY_Z,1))/BDFcoeffs[1];
//       }
//   }
 


// template<  unsigned int TDim>
//   void TwoStepUpdatedLagrangianVPSolidElement<TDim>::CalcMeanPressure(double &meanPressure,
// 								   const int Step)
//   {
//     GeometryType& rGeom = this->GetGeometry();
//     const SizeType NumNodes = rGeom.PointsNumber();
//     meanPressure=0;
//     double coeff=0.0;
//     for (SizeType i = 0; i < NumNodes; ++i)
//       {
// 	meanPressure+=rGeom[i].GetSolutionStepValue(PRESSURE,Step);
// 	coeff+=1.0;
//       }
//     meanPressure*=1.0/coeff;
//   }



// template<>
//   void TwoStepUpdatedLagrangianVPSolidElement<2>::GetPositions(Vector& rValues)
//   {
//     GeometryType& rGeom = this->GetGeometry();
//     const SizeType NumNodes = rGeom.PointsNumber();
//     const SizeType LocalSize = 2*NumNodes;

//     if (rValues.size() != LocalSize) rValues.resize(LocalSize);

//     SizeType Index = 0;

//     for (SizeType i = 0; i < NumNodes; ++i)
//       {
//         // rValues[Index++] = rGeom[i].X();
//         // rValues[Index++] = rGeom[i].Y();
// 	rValues[Index++] = rGeom[i].X0()+rGeom[i].GetSolutionStepValue(DISPLACEMENT_X);
//         rValues[Index++] = rGeom[i].Y0()+rGeom[i].GetSolutionStepValue(DISPLACEMENT_Y);
//       }
//   }



//   template<>
//   void TwoStepUpdatedLagrangianVPSolidElement<3>::GetPositions(Vector& rValues)
//   {
//     GeometryType& rGeom = this->GetGeometry();
//     const SizeType NumNodes = rGeom.PointsNumber();
//     const SizeType LocalSize = 3*NumNodes;

//     if (rValues.size() != LocalSize) rValues.resize(LocalSize);

//     SizeType Index = 0;

//     for (SizeType i = 0; i < NumNodes; ++i)
//       {
//  	rValues[Index++] = rGeom[i].X0()+rGeom[i].GetSolutionStepValue(DISPLACEMENT_X);
//         rValues[Index++] = rGeom[i].Y0()+rGeom[i].GetSolutionStepValue(DISPLACEMENT_Y);
//         rValues[Index++] = rGeom[i].Z0()+rGeom[i].GetSolutionStepValue(DISPLACEMENT_Z);
//       }
//   }


//   template<>
//   void TwoStepUpdatedLagrangianVPSolidElement<2>::GetVelocityValues(Vector& rValues,
// 								    const int Step)
//   {
//     GeometryType& rGeom = this->GetGeometry();
//     const SizeType NumNodes = rGeom.PointsNumber();
//     const SizeType LocalSize = 2*NumNodes;

//     if (rValues.size() != LocalSize) rValues.resize(LocalSize);

//     SizeType Index = 0;

//     for (SizeType i = 0; i < NumNodes; ++i)
//       {
//         rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_X,Step);
//         rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y,Step);
//       }
//   }


//   template<>
//   void TwoStepUpdatedLagrangianVPSolidElement<3>::GetVelocityValues(Vector& rValues,
// 								    const int Step)
//   {
//     GeometryType& rGeom = this->GetGeometry();
//     const SizeType NumNodes = rGeom.PointsNumber();
//     const SizeType LocalSize = 3*NumNodes;

//     if (rValues.size() != LocalSize) rValues.resize(LocalSize);

//     SizeType Index = 0;

//     for (SizeType i = 0; i < NumNodes; ++i)
//       {
//         rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_X,Step);
//         rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Y,Step);
//         rValues[Index++] = rGeom[i].FastGetSolutionStepValue(VELOCITY_Z,Step);
//       }
//   }



//   template< >
//   void TwoStepUpdatedLagrangianVPSolidElement<2>::GetAccelerationValues(Vector& rValues,
// 									const int Step)
//   {
//     GeometryType& rGeom = this->GetGeometry();
//     const SizeType NumNodes = rGeom.PointsNumber();
//     const SizeType LocalSize = 2*NumNodes;

//     if (rValues.size() != LocalSize) rValues.resize(LocalSize);

//     SizeType Index = 0;

//     for (SizeType i = 0; i < NumNodes; ++i)
//       {
//         rValues[Index++] = rGeom[i].FastGetSolutionStepValue(ACCELERATION_X,Step);
//         rValues[Index++] = rGeom[i].FastGetSolutionStepValue(ACCELERATION_Y,Step);
//       }
//   }

//   template< >
//   void TwoStepUpdatedLagrangianVPSolidElement<3>::GetAccelerationValues(Vector& rValues,
// 									const int Step)
//   {
//     GeometryType& rGeom = this->GetGeometry();
//     const SizeType NumNodes = rGeom.PointsNumber();
//     const SizeType LocalSize = 3*NumNodes;

//     if (rValues.size() != LocalSize) rValues.resize(LocalSize);

//     SizeType Index = 0;

//     for (SizeType i = 0; i < NumNodes; ++i)
//       {
//         rValues[Index++] = rGeom[i].FastGetSolutionStepValue(ACCELERATION_X,Step);
//         rValues[Index++] = rGeom[i].FastGetSolutionStepValue(ACCELERATION_Y,Step);
//         rValues[Index++] = rGeom[i].FastGetSolutionStepValue(ACCELERATION_Z,Step);
//       }
//   }





//   template< unsigned int TDim >
//   void TwoStepUpdatedLagrangianVPSolidElement<TDim>::CalculateGeometryData(ShapeFunctionDerivativesArrayType &rDN_DX,
// 									   Matrix &NContainer,
// 									   Vector &rGaussWeights)
//   {
//     const GeometryType& rGeom = this->GetGeometry();
//     Vector DetJ;
//     rGeom.ShapeFunctionsIntegrationPointsGradients(rDN_DX,DetJ,GeometryData::GI_GAUSS_1);
//     NContainer = rGeom.ShapeFunctionsValues(GeometryData::GI_GAUSS_1);
//     const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::GI_GAUSS_1);

//     rGaussWeights.resize(rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_1),false);

//     for (unsigned int g = 0; g < rGeom.IntegrationPointsNumber(GeometryData::GI_GAUSS_1); g++){
//       rGaussWeights[g] = DetJ[g] * IntegrationPoints[g].Weight();
  
//     }
    

//   }

//   template< unsigned int TDim >
//   double TwoStepUpdatedLagrangianVPSolidElement<TDim>::ElementSize(/*ShapeFunctionDerivativesType &rDN_DX*/)
//   {
//     const GeometryType& rGeom = this->GetGeometry();
//     const SizeType NumNodes = rGeom.PointsNumber();

//     // calculate minimum element length (used in stabilization Tau)
//     array_1d<double,3> Edge(3,0.0);
//     Edge = rGeom[1].Coordinates() - rGeom[0].Coordinates();
//     double ElemSize = Edge[0]*Edge[0];
//     for (SizeType d = 1; d < TDim; d++)
//       ElemSize += Edge[d]*Edge[d];

//     for (SizeType i = 2; i < NumNodes; i++)
//       for(SizeType j = 0; j < i; j++)
//         {
// 	  Edge = rGeom[i].Coordinates() - rGeom[j].Coordinates();
// 	  double Length = Edge[0]*Edge[0];
// 	  for (SizeType d = 1; d < TDim; d++)
// 	    Length += Edge[d]*Edge[d];
// 	  if (Length < ElemSize) ElemSize = Length;
//         }
//     return sqrt(ElemSize);

//   }



  // template< unsigned int TDim >
  // double TwoStepUpdatedLagrangianVPSolidElement<TDim>::EquivalentStrainRate(const ShapeFunctionDerivativesType &rDN_DX) const
  // {
  //   const GeometryType& rGeom = this->GetGeometry();
  //   const unsigned int NumNodes = rGeom.PointsNumber();

  //   // Calculate Symetric gradient
  //   MatrixType S = ZeroMatrix(TDim,TDim);
  //   for (unsigned int n = 0; n < NumNodes; ++n)
  //     {
  //       const array_1d<double,3>& rVel = rGeom[n].FastGetSolutionStepValue(VELOCITY,1); // OLD VELOCITY (which is incompressible, unlike the fractional step one)
  //       for (unsigned int i = 0; i < TDim; ++i)
  // 	  for (unsigned int j = 0; j < TDim; ++j)
  // 	    S(i,j) += 0.5 * ( rDN_DX(n,j) * rVel[i] + rDN_DX(n,i) * rVel[j] );
  //     }

  //   // Norm of symetric gradient
  //   double NormS = 0.0;
  //   for (unsigned int i = 0; i < TDim; ++i)
  //     for (unsigned int j = 0; j < TDim; ++j)
  // 	NormS += S(i,j) * S(i,j);

  //   return std::sqrt(2.0*NormS);
  // }



  // template< unsigned int TDim >
  // void TwoStepUpdatedLagrangianVPSolidElement<TDim>::ComputeLumpedMassMatrix(Matrix& rMassMatrix,
  // 									     const double Weight)
  // {
  //   const SizeType NumNodes = this->GetGeometry().PointsNumber();

  //   double coeff=1.0+TDim;
  //   for (SizeType i = 0; i < NumNodes; ++i)
  //     {
    
  // 	double Mij = Weight/coeff;

  //       for ( unsigned int j = 0; j <  TDim; j++ )
  // 	  {
  //           unsigned int index = i * TDim + j;
  //           rMassMatrix( index, index ) += Mij;
  // 	  }

  //     }
  // }




  // template< unsigned int TDim >
  // void TwoStepUpdatedLagrangianVPSolidElement<TDim>::ComputeMomentumMassTerm(Matrix& rMassMatrix,
  // 									     const ShapeFunctionsType& rN,
  // 									     const double Weight)
  // {
  //   const SizeType NumNodes = this->GetGeometry().PointsNumber();

  //   IndexType FirstRow = 0;
  //   IndexType FirstCol = 0;

  //   for (SizeType i = 0; i < NumNodes; ++i)
  //     {
  //       for (SizeType j = 0; j < NumNodes; ++j)
  // 	  {
  //           const double Mij = Weight * rN[i] * rN[j];
  //           for (SizeType d =  0; d < TDim; ++d)
  // 	      rMassMatrix(FirstRow+d,FirstCol+d) += Mij;
  //           FirstCol += TDim;
  // 	  }
  //       FirstRow += TDim;
  //       FirstCol = 0;
  //     }
  // }

 
  // template< unsigned int TDim >
  // void TwoStepUpdatedLagrangianVPSolidElement<TDim>::AddExternalForces(Vector& rRHSVector,
  // 								       const double Density,
  // 								       const array_1d<double,3>& rBodyForce,
  // 								       const ShapeFunctionsType& rN,
  // 								       const double Weight)
  // {
  //   const SizeType NumNodes = this->GetGeometry().PointsNumber();

  //   SizeType FirstRow = 0;

  //   for (SizeType i = 0; i < NumNodes; ++i)
  //     {

  //       // Build RHS
  //       for (SizeType d = 0; d < TDim; ++d)
  // 	  {
  //           // Body force
  //           // double RHSi = Density * rN[i] * rBodyForce[d];
  //           double RHSi = 0;

  // 	    if(d==2){
  // 	      RHSi = - Density * rN[i] * 9.81;
  // 	    } 
  // 	    rRHSVector[FirstRow+d] += Weight * RHSi;
  // 	  }

  //       FirstRow += TDim;
  //     }
  // }

  // template< >
  // void TwoStepUpdatedLagrangianVPSolidElement<2>::AddDeviatoricInternalForces(Vector& rRHSVector,const ShapeFunctionDerivativesType& rDN_DX, ElementalVariables& rElementalVariables, const double Weight)
  // {
 
  // }
  // template< >
  // void TwoStepUpdatedLagrangianVPSolidElement<3>::AddDeviatoricInternalForces(Vector& rRHSVector,const ShapeFunctionDerivativesType& rDN_DX, ElementalVariables& rElementalVariables, const double Weight)
  // {
 
  // }

  // template< >
  // void TwoStepUpdatedLagrangianVPSolidElement<2>::AddInternalForces(Vector& rRHSVector,
  // 								    const ShapeFunctionDerivativesType& rDN_DX,
  // 								    ElementalVariables& rElementalVariables,
  // 								    const double Weight)
  // {
  //   const SizeType NumNodes = this->GetGeometry().PointsNumber();

  //   SizeType FirstRow = 0;

  //   MatrixType invGradDef=rElementalVariables.InvFgrad;

  //   for (SizeType i = 0; i < NumNodes; ++i)
  //     {
  // 	double lagDNXi=rDN_DX(i,0)*invGradDef(0,0)+rDN_DX(i,1)*invGradDef(1,0);
  // 	double lagDNYi=rDN_DX(i,0)*invGradDef(0,1)+rDN_DX(i,1)*invGradDef(1,1);
  // 	// lagDNXi=rDN_DX(i,0);
  // 	// lagDNYi=rDN_DX(i,1);

  // 	rRHSVector[FirstRow]   += -Weight*(lagDNXi*rElementalVariables.UpdatedTotalCauchyStress[0]+
  // 					   lagDNYi*rElementalVariables.UpdatedTotalCauchyStress[2]);

  // 	rRHSVector[FirstRow+1] += -Weight*(lagDNYi*rElementalVariables.UpdatedTotalCauchyStress[1]+
  // 					   lagDNXi*rElementalVariables.UpdatedTotalCauchyStress[2]);

  // 	FirstRow += 2;
  //     }
  // }

  // template< >
  // void TwoStepUpdatedLagrangianVPSolidElement<3>::AddInternalForces(Vector& rRHSVector,
  // 								    const ShapeFunctionDerivativesType& rDN_DX,
  // 								    ElementalVariables& rElementalVariables,
  // 								    const double Weight)
  // {

  //   const SizeType NumNodes = this->GetGeometry().PointsNumber();

  //   SizeType FirstRow = 0;

  //   MatrixType invGradDef=rElementalVariables.InvFgrad;

  //   for (SizeType i = 0; i < NumNodes; ++i)
  //     {
  // 	double lagDNXi=rDN_DX(i,0)*invGradDef(0,0)+rDN_DX(i,1)*invGradDef(1,0)+rDN_DX(i,2)*invGradDef(2,0);
  // 	double lagDNYi=rDN_DX(i,0)*invGradDef(0,1)+rDN_DX(i,1)*invGradDef(1,1)+rDN_DX(i,2)*invGradDef(2,1);
  // 	double lagDNZi=rDN_DX(i,0)*invGradDef(0,2)+rDN_DX(i,1)*invGradDef(1,2)+rDN_DX(i,2)*invGradDef(2,2);

  // 	// lagDNXi=rDN_DX(i,0);
  // 	// lagDNYi=rDN_DX(i,1);
  // 	// lagDNZi=rDN_DX(i,2);

  // 	rRHSVector[FirstRow]   += -Weight*(lagDNXi*rElementalVariables.UpdatedTotalCauchyStress[0]+
  // 					   lagDNYi*rElementalVariables.UpdatedTotalCauchyStress[3]+
  // 					   lagDNZi*rElementalVariables.UpdatedTotalCauchyStress[4]);

  // 	rRHSVector[FirstRow+1] += -Weight*(lagDNYi*rElementalVariables.UpdatedTotalCauchyStress[1]+
  // 					   lagDNXi*rElementalVariables.UpdatedTotalCauchyStress[3]+
  // 					   lagDNZi*rElementalVariables.UpdatedTotalCauchyStress[5]);

  // 	rRHSVector[FirstRow+2] += -Weight*(lagDNZi*rElementalVariables.UpdatedTotalCauchyStress[2]+
  // 					   lagDNXi*rElementalVariables.UpdatedTotalCauchyStress[4]+
  // 					   lagDNYi*rElementalVariables.UpdatedTotalCauchyStress[5]);

  // 	FirstRow += 3;
  //     }


  // }

  // template<  unsigned int TDim >
  // void TwoStepUpdatedLagrangianVPSolidElement<TDim>::AddDynamicForces(Vector& rRHSVector,
  // 								   const double Weight,
  // 								   const double TimeStep)
  // {

  //   double coeff=1.0+TDim;

  //   const SizeType NumNodes = this->GetGeometry().PointsNumber();
  //   const SizeType LocalSize = TDim * NumNodes;

  //   VectorType VelocityValues = ZeroVector(LocalSize);
  //   this->GetVelocityValues(VelocityValues,0);
  //   VectorType UpdatedAccelerations = ZeroVector(LocalSize);
  //   UpdatedAccelerations = 2.0*VelocityValues/TimeStep;
  //   VectorType LastAccValues = ZeroVector(LocalSize);
  //   this->GetAccelerationValues(LastAccValues,0);
  //   this->GetVelocityValues(VelocityValues,1);
  //   UpdatedAccelerations += -2.0*VelocityValues/TimeStep - LastAccValues; 

  //   SizeType FirstRow = 0;
  //   for (SizeType i = 0; i < NumNodes; ++i)
  //     {
  //       // for ( unsigned int j = 0; j <  TDim; j++ )
  // 	//   {
  //       //     unsigned int index = i * TDim + j;
  //       //     rRHSVector[index] += - (Weight /coeff) * UpdatedAccelerations[index];
  // 	//   }
  // 	rRHSVector[FirstRow]   += - (Weight /coeff) * UpdatedAccelerations[FirstRow] ;
  // 	rRHSVector[FirstRow+1] += - (Weight /coeff) * UpdatedAccelerations[FirstRow+1];
  // 	FirstRow += 2;
  //     }
  // }


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
  void TwoStepUpdatedLagrangianVPSolidElement<TDim>::ComputeBulkMatrixForPressureVel(Matrix& BulkVelMatrix,
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
  void TwoStepUpdatedLagrangianVPSolidElement<TDim>::ComputeBulkMatrixForPressureAcc(Matrix& BulkAccMatrix,
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

    // Evaluate the pressure and pressure gradient at this point (for the G * P_n term)
    array_1d<double,TDim> OldPressureGradient(TDim,0.0);
    this->EvaluateGradientInPoint(OldPressureGradient,PRESSURE,rDN_DX);
    double RHSi = 0;
    for (SizeType d = 0; d < TDim; ++d)
      {
	RHSi += rDN_DX(i,d) * Tau * ( Density  * ( BodyForce[d] ) - OldPressureGradient[d]  );
	
	if(d==2){
	  RHSi += - rDN_DX(i,d) * Tau * ( Density * 9.81 );
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
void TwoStepUpdatedLagrangianVPSolidElement<TDim>::CalcMechanicsUpdated(ElementalVariables & rElementalVariables,
									const ProcessInfo& rCurrentProcessInfo,
									unsigned int g)
{

  this->CalcStrainRateUpdated(rElementalVariables,rCurrentProcessInfo,g);
  const double TimeStep=0.5/rCurrentProcessInfo[BDF_COEFFICIENTS][2];
  this->CalcElasticPlasticCauchySplitted(rElementalVariables,TimeStep,g);

} 



// template< unsigned int TDim>
// void TwoStepUpdatedLagrangianVPSolidElement<TDim>::CalculateDeltaPosition(Matrix & rDeltaPosition)
// {
//     KRATOS_TRY

//     const unsigned int number_of_nodes = this->GetGeometry().PointsNumber();
//     unsigned int dimension = this->GetGeometry().WorkingSpaceDimension();

//     rDeltaPosition = zero_matrix<double>( number_of_nodes , dimension);

//     for ( unsigned int i = 0; i < number_of_nodes; i++ )
//     {
//         array_1d<double, 3 > & CurrentDisplacement  = this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
//         array_1d<double, 3 > & PreviousDisplacement = this->GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);

//         for ( unsigned int j = 0; j < dimension; j++ )
//         {
//             rDeltaPosition(i,j) = CurrentDisplacement[j]-PreviousDisplacement[j];
//         }
//     }

//     KRATOS_CATCH( "" )
// }



// template< unsigned int TDim>
// void TwoStepUpdatedLagrangianVPSolidElement<TDim>::CalcStrainRateUpdated(ElementalVariables & rElementalVariables,
// 									 const ProcessInfo &rCurrentProcessInfo,
// 									 unsigned int g)
// {

//   ShapeFunctionDerivativesArrayType DN_DX;
//   Matrix NContainer;
//   VectorType GaussWeights;
//   this->CalculateGeometryData(DN_DX,NContainer,GaussWeights);

//   // const GeometryType& rGeom = this->this->GetGeometry();
//   const ShapeFunctionDerivativesType& rDN_DX = DN_DX[g];


//   this->CalcFGrad(rDN_DX,rElementalVariables.Fgrad,
// 		  rElementalVariables.InvFgrad,
// 		  rElementalVariables.DetFgrad,
// 		  rCurrentProcessInfo);


//   //it computes the material time derivative of the deformation gradient and its jacobian and inverse
//   this->CalcVelDefGrad(rDN_DX,rElementalVariables.FgradVel,
// 		       rElementalVariables.InvFgradVel,
// 		       rElementalVariables.DetFgradVel);

//   //it computes the spatial velocity gradient tensor --> [l]
//   this->CalcSpatialVelocityGrad(rElementalVariables.InvFgrad,
// 				rElementalVariables.FgradVel,
// 				rElementalVariables.SpatialVelocityGrad);
  
//   this->CalcVolDefRateFromSpatialVelGrad(rElementalVariables.VolumetricDefRate,
// 					 rElementalVariables.SpatialVelocityGrad);
  

//   // this->CalcVolumetricDefRate(rDN_DX,
//   // 			      rElementalVariables.VolumetricDefRate,
//   // 			      rElementalVariables.InvFgrad);

//   // //it checks whether tr(l) == div(v)
//   // CheckStrain1(rElementalVariables.VolumetricDefRate,
//   // 	       rElementalVariables.SpatialVelocityGrad);
      
//   // CheckStrain2(rElementalVariables.SpatialVelocityGrad,
//   // 	       rElementalVariables.Fgrad,
//   // 	       rElementalVariables.FgradVel);
 
//   //it computes Material time Derivative of Green Lagrange strain tensor in MATERIAL configuration --> [D(E)/Dt]
//   this->CalcMDGreenLagrangeMaterial(rElementalVariables.Fgrad,
// 				    rElementalVariables.FgradVel,
// 				    rElementalVariables.MDGreenLagrangeMaterial);
      
//   //it computes Material time Derivative of Green Lagrange strain tensor in SPATIAL configuration  --> [d]
//   this->CalcSpatialDefRate(rElementalVariables.MDGreenLagrangeMaterial,
// 			   rElementalVariables.InvFgrad,
// 			   rElementalVariables.SpatialDefRate);

//   CheckStrain3(rElementalVariables.SpatialDefRate,
//   	       rElementalVariables.SpatialVelocityGrad);

//   this->CalcDeviatoricInvariant(rElementalVariables.SpatialDefRate,
// 				rElementalVariables.DeviatoricInvariant);

// }  


// template <unsigned int TDim > 
// void TwoStepUpdatedLagrangianVPSolidElement<TDim>::CalcFGrad(const ShapeFunctionDerivativesType& rDN_DX,
// 							     MatrixType &Fgrad,
// 							     MatrixType &invFgrad,
// 							     double &FJacobian,
// 							     const ProcessInfo& rCurrentProcessInfo)
// {
//   GeometryType& rGeom = this->GetGeometry();
//   const SizeType NumNodes = rGeom.PointsNumber();
//   const SizeType LocalSize = TDim*NumNodes;
//   VectorType  NodePosition= ZeroVector(LocalSize);
//   this->GetUpdatedPositions(NodePosition,rCurrentProcessInfo);
//   // this->GetPositions(NodePosition);

//   Fgrad.resize(TDim,TDim);

//   Fgrad=ZeroMatrix(TDim,TDim);
//   for (SizeType i = 0; i < TDim; i++)
//     {
//       for (SizeType j = 0; j < TDim; j++)
// 	{
// 	  for (SizeType k = 0; k < NumNodes; k++)
// 	    {
// 	      Fgrad(i,j)+= NodePosition[TDim*k+i]*rDN_DX(k,j);
// 	    }
// 	}
//     }


//  //Inverse
//   invFgrad.resize(TDim,TDim);
  
//   MathUtils<double>::InvertMatrix( Fgrad, invFgrad, FJacobian );

//   // Fgrad.resize(2,2);

//   // Fgrad(0,0)= NodePosition[0]*rDN_DX(0,0)+NodePosition[2]*rDN_DX(1,0)+NodePosition[4]*rDN_DX(2,0);
//   // Fgrad(0,1)= NodePosition[0]*rDN_DX(0,1)+NodePosition[2]*rDN_DX(1,1)+NodePosition[4]*rDN_DX(2,1);
//   // Fgrad(1,0)= NodePosition[1]*rDN_DX(0,0)+NodePosition[3]*rDN_DX(1,0)+NodePosition[5]*rDN_DX(2,0);
//   // Fgrad(1,1)= NodePosition[1]*rDN_DX(0,1)+NodePosition[3]*rDN_DX(1,1)+NodePosition[5]*rDN_DX(2,1);

//   // //Determinant of the material time derivative of the deformation gradient tensor
//   // FJacobian= Fgrad(0,0)*Fgrad(1,1)-Fgrad(0,1)*Fgrad(1,0);
   
//   // //Inverse
//   // invFgrad.resize(2,2);
//   // if(FJacobian==0){
//   //   FJacobian=1;
//   // }else{
//   //   invFgrad(0,0)=  (Fgrad(1,1)/FJacobian);
//   //   invFgrad(0,1)= -(Fgrad(0,1)/FJacobian);
//   //   invFgrad(1,0)= -(Fgrad(1,0)/FJacobian);
//   //   invFgrad(1,1)=  (Fgrad(0,0)/FJacobian); 
//   // }

// }




// template < unsigned int TDim> 
// void TwoStepUpdatedLagrangianVPSolidElement<TDim>::CalcVelDefGrad(const ShapeFunctionDerivativesType& rDN_DX,
// 								  MatrixType &FgradVel,
// 								  MatrixType &invFgradVel,
// 								  double &FVelJacobian)
// {
//   GeometryType& rGeom = this->GetGeometry();
//   const SizeType NumNodes = rGeom.PointsNumber();
//   const SizeType LocalSize = TDim*NumNodes;
//   VectorType CurrentVelocities = ZeroVector(LocalSize);
//   VectorType UpdatedVelocities = ZeroVector(LocalSize);
//   VectorType RHSVelocities = ZeroVector(LocalSize);
//   this->GetVelocityValues(CurrentVelocities,1);
//   this->GetVelocityValues(UpdatedVelocities,0);
//   RHSVelocities=CurrentVelocities*0.5+UpdatedVelocities*0.5;

//   FgradVel.resize(TDim,TDim);

//   FgradVel=ZeroMatrix(TDim,TDim);
//   for (SizeType i = 0; i < TDim; i++)
//     {
//       for (SizeType j = 0; j < TDim; j++)
// 	{
// 	  for (SizeType k = 0; k < NumNodes; k++)
// 	    {
// 	      FgradVel(i,j)+= RHSVelocities[TDim*k+i]*rDN_DX(k,j);
// 	    }
// 	}
//     }


//  //Inverse
//   invFgradVel.resize(TDim,TDim);
  
//   MathUtils<double>::InvertMatrix( FgradVel, invFgradVel, FVelJacobian );



//   // FgradVel.resize(2,2);


//   // FgradVel(0,0)= RHSVelocities[0]*rDN_DX(0,0)+RHSVelocities[2]*rDN_DX(1,0)+RHSVelocities[4]*rDN_DX(2,0);
//   // FgradVel(0,1)= RHSVelocities[0]*rDN_DX(0,1)+RHSVelocities[2]*rDN_DX(1,1)+RHSVelocities[4]*rDN_DX(2,1);
//   // FgradVel(1,0)= RHSVelocities[1]*rDN_DX(0,0)+RHSVelocities[3]*rDN_DX(1,0)+RHSVelocities[5]*rDN_DX(2,0);
//   // FgradVel(1,1)= RHSVelocities[1]*rDN_DX(0,1)+RHSVelocities[3]*rDN_DX(1,1)+RHSVelocities[5]*rDN_DX(2,1);

//   // //Determinant of the material time derivative of the deformation gradient tensor
//   // FVelJacobian= FgradVel(0,0)*FgradVel(1,1)-FgradVel(0,1)*FgradVel(1,0);
   
//   // //Inverse
//   // invFgradVel.resize(2,2);
//   // if(FVelJacobian==0){
//   //   FVelJacobian=1;
//   // }else{
//   //   invFgradVel(0,0)=  (FgradVel(1,1)/FVelJacobian);
//   //   invFgradVel(0,1)= -(FgradVel(0,1)/FVelJacobian);
//   //   invFgradVel(1,0)= -(FgradVel(1,0)/FVelJacobian);
//   //   invFgradVel(1,1)=  (FgradVel(0,0)/FVelJacobian); 
//   // }

// }


// template < > 
// void TwoStepUpdatedLagrangianVPSolidElement<2>::CalcVolumetricDefRate(const ShapeFunctionDerivativesType& rDN_DX,
// 								      double &volumetricDefRate,
// 								      MatrixType &invGradDef)
// {


//   GeometryType& rGeom = this->GetGeometry();
//   const SizeType NumNodes = rGeom.PointsNumber();
//   const SizeType LocalSize = 2*NumNodes;
//   VectorType CurrentVelocities = ZeroVector(LocalSize);
//   VectorType UpdatedVelocities = ZeroVector(LocalSize);
//   VectorType RHSVelocities = ZeroVector(LocalSize);
//   this->GetVelocityValues(CurrentVelocities,1);
//   this->GetVelocityValues(UpdatedVelocities,0);
//   RHSVelocities=CurrentVelocities*0.5+UpdatedVelocities*0.5;

//   double lagDNX0=rDN_DX(0,0)*invGradDef(0,0)+rDN_DX(0,1)*invGradDef(1,0);
//   double lagDNX1=rDN_DX(1,0)*invGradDef(0,0)+rDN_DX(1,1)*invGradDef(1,0);
//   double lagDNX2=rDN_DX(2,0)*invGradDef(0,0)+rDN_DX(2,1)*invGradDef(1,0);
//   double lagDNY0=rDN_DX(0,0)*invGradDef(0,1)+rDN_DX(0,1)*invGradDef(1,1);
//   double lagDNY1=rDN_DX(1,0)*invGradDef(0,1)+rDN_DX(1,1)*invGradDef(1,1);
//   double lagDNY2=rDN_DX(2,0)*invGradDef(0,1)+rDN_DX(2,1)*invGradDef(1,1);


//   volumetricDefRate= lagDNX0*RHSVelocities[0] + lagDNX1*RHSVelocities[2] + lagDNX2*RHSVelocities[4];
//   volumetricDefRate+=lagDNY0*RHSVelocities[1] + lagDNY1*RHSVelocities[3] + lagDNY2*RHSVelocities[5];

// }

// template < unsigned int TDim > 
// void TwoStepUpdatedLagrangianVPSolidElement<TDim>::CalcSpatialVelocityGrad(MatrixType &invFgrad,
// 									   MatrixType &VelDefgrad,
// 									   MatrixType &SpatialVelocityGrad)
// {
//   SpatialVelocityGrad.resize(TDim,TDim);
  
//   SpatialVelocityGrad=prod(VelDefgrad,invFgrad);

//   // SpatialVelocityGrad(0,0)=VelDefgrad(0,0)*invFgrad(0,0) + VelDefgrad(0,1)*invFgrad(1,0);
//   // SpatialVelocityGrad(0,1)=VelDefgrad(0,0)*invFgrad(0,1) + VelDefgrad(0,1)*invFgrad(1,1);
//   // SpatialVelocityGrad(1,0)=VelDefgrad(1,0)*invFgrad(0,0) + VelDefgrad(1,1)*invFgrad(1,0);
//   // SpatialVelocityGrad(1,1)=VelDefgrad(1,0)*invFgrad(0,1) + VelDefgrad(1,1)*invFgrad(1,1);

// }


// template < unsigned int TDim > 
// void TwoStepUpdatedLagrangianVPSolidElement<TDim>::CalcVolDefRateFromSpatialVelGrad(double &volumetricDefRate,
// 										    MatrixType &SpatialVelocityGrad)
// {
//   volumetricDefRate=0;
//   for (SizeType i = 0; i < TDim; i++)
//     {
//       volumetricDefRate+=SpatialVelocityGrad(i,i);
//     }
// }



// template < > 
// void TwoStepUpdatedLagrangianVPSolidElement<2>::CheckStrain1(double &VolumetricDefRate,
// 							     MatrixType &SpatialVelocityGrad)
// {
//   double trace_l=SpatialVelocityGrad(0,0)+SpatialVelocityGrad(1,1);
//   if(fabs(trace_l-VolumetricDefRate)<0.0000001){
//   }else{
//     std::cout<<" ERROR IN CHECKSTRAIN(1) -> ";
//     std::cout<<"trace_l= "<<trace_l<<" VolDefRate"<<VolumetricDefRate<<std::endl;
//   }
// }

// template < > 
// void TwoStepUpdatedLagrangianVPSolidElement<2>::CalcMDGreenLagrangeMaterial(MatrixType &Fgrad,
// 									    MatrixType &VelDefgrad, 
// 									    VectorType &MDGreenLagrangeMaterial)
// {
//   // x-component
//   MDGreenLagrangeMaterial[0]=VelDefgrad(0,0)*Fgrad(0,0) + VelDefgrad(1,0)*Fgrad(1,0);
//   // y-component
//   MDGreenLagrangeMaterial[1]=VelDefgrad(1,1)*Fgrad(1,1) + VelDefgrad(0,1)*Fgrad(0,1);
//   // xy-component
//   MDGreenLagrangeMaterial[2]=(VelDefgrad(0,0)*Fgrad(0,1) + VelDefgrad(1,0)*Fgrad(1,1) +
// 			      VelDefgrad(0,1)*Fgrad(0,0) + VelDefgrad(1,1)*Fgrad(1,0))*0.5;
// }



// template < > 
// void TwoStepUpdatedLagrangianVPSolidElement<3>::CalcMDGreenLagrangeMaterial(MatrixType &Fgrad,
// 									    MatrixType &VelDefgrad, 
// 									    VectorType &MDGreenLagrangeMaterial)
// {
//   MatrixType FgradTransp(3,3);
//   MatrixType VelDefgradTransp(3,3);
//   MatrixType part1(3,3);
//   MatrixType part2(3,3);

//   FgradTransp=Fgrad;
//   FgradTransp(0,1)=Fgrad(1,0);
//   FgradTransp(0,2)=Fgrad(2,0);
//   FgradTransp(1,0)=Fgrad(0,1);
//   FgradTransp(1,2)=Fgrad(2,1);
//   FgradTransp(2,0)=Fgrad(0,2);
//   FgradTransp(2,1)=Fgrad(1,2);

//   VelDefgradTransp=VelDefgrad;
//   VelDefgradTransp(0,1)=VelDefgrad(1,0);
//   VelDefgradTransp(0,2)=VelDefgrad(2,0);
//   VelDefgradTransp(1,0)=VelDefgrad(0,1);
//   VelDefgradTransp(1,2)=VelDefgrad(2,1);
//   VelDefgradTransp(2,0)=VelDefgrad(0,2);
//   VelDefgradTransp(2,1)=VelDefgrad(1,2);

//   part1=prod(VelDefgradTransp,Fgrad);
//   part2=prod(FgradTransp,VelDefgrad);

//   MDGreenLagrangeMaterial[0]= ( part1(0,0) + part2(0,0) ) * 0.5;  //xx-component
//   MDGreenLagrangeMaterial[1]= ( part1(1,1) + part2(1,1) ) * 0.5;  //yy-component
//   MDGreenLagrangeMaterial[2]= ( part1(2,2) + part2(2,2) ) * 0.5;  //zz-component
//   MDGreenLagrangeMaterial[3]= ( part1(0,1) + part2(0,1) ) * 0.5;  //xy-component
//   MDGreenLagrangeMaterial[4]= ( part1(0,2) + part2(0,2) ) * 0.5;  //xz-component
//   MDGreenLagrangeMaterial[5]= ( part1(1,2) + part2(1,2) ) * 0.5;  //yz-component

// }





// template < > 
// void TwoStepUpdatedLagrangianVPSolidElement<2>::CalcSpatialDefRate(VectorType &MDGreenLagrangeMaterial,
// 								   MatrixType &invFgrad,
// 								   VectorType &SpatialDefRate)
// {
//   // x-component
//   SpatialDefRate[0]= invFgrad(0,0)*MDGreenLagrangeMaterial[0]*invFgrad(0,0) + 
//     invFgrad(1,0)*MDGreenLagrangeMaterial[2]*invFgrad(0,0)*2 +
//     invFgrad(1,0)*MDGreenLagrangeMaterial[1]*invFgrad(1,0);
//   // y-component
//   SpatialDefRate[1]= invFgrad(0,1)*MDGreenLagrangeMaterial[0]*invFgrad(0,1) + 
//     invFgrad(0,1)*MDGreenLagrangeMaterial[2]*invFgrad(1,1)*2 +
//     invFgrad(1,1)*MDGreenLagrangeMaterial[1]*invFgrad(1,1);
//   // xy-component
//   SpatialDefRate[2]=invFgrad(0,0)*MDGreenLagrangeMaterial[0]*invFgrad(0,1) + 
//     invFgrad(0,0)*MDGreenLagrangeMaterial[2]*invFgrad(1,1) +
//     invFgrad(1,0)*MDGreenLagrangeMaterial[2]*invFgrad(0,1) +
//     invFgrad(1,0)*MDGreenLagrangeMaterial[1]*invFgrad(1,1);
// }


// template < > 
// void TwoStepUpdatedLagrangianVPSolidElement<3>::CalcSpatialDefRate(VectorType &MDGreenLagrangeMaterial,
// 								   MatrixType &invFgrad,
// 								   VectorType &SpatialDefRate)
// {
//   MatrixType MDGLM(3,3);
//   MatrixType invFgradTransp(3,3);
//   MatrixType part1(3,3);
//   MatrixType totalMatrix(3,3);

//   invFgradTransp=invFgrad;
//   invFgradTransp(0,1)=invFgrad(1,0);
//   invFgradTransp(0,2)=invFgrad(2,0);
//   invFgradTransp(1,0)=invFgrad(0,1);
//   invFgradTransp(1,2)=invFgrad(2,1);
//   invFgradTransp(2,0)=invFgrad(0,2);
//   invFgradTransp(2,1)=invFgrad(1,2);

//   MDGLM(0,0)=MDGreenLagrangeMaterial[0];  //XX-component;
//   MDGLM(1,1)=MDGreenLagrangeMaterial[1];  //YY-component;
//   MDGLM(2,2)=MDGreenLagrangeMaterial[2];  //ZZ-component;
//   MDGLM(0,1)=MDGreenLagrangeMaterial[3];  //XY-component;
//   MDGLM(1,0)=MDGreenLagrangeMaterial[3];  //XY-component;
//   MDGLM(0,2)=MDGreenLagrangeMaterial[4];  //ZX-component;
//   MDGLM(2,0)=MDGreenLagrangeMaterial[4];  //ZX-component;
//   MDGLM(1,2)=MDGreenLagrangeMaterial[5];  //YZ-component;
//   MDGLM(2,1)=MDGreenLagrangeMaterial[5];  //YZ-component;

//   part1=prod(MDGLM,invFgrad);

//   totalMatrix=prod(invFgradTransp,part1);
 
//   SpatialDefRate[0]=totalMatrix(0,0);
//   SpatialDefRate[1]=totalMatrix(1,1);
//   SpatialDefRate[2]=totalMatrix(2,2);
//   SpatialDefRate[3]=totalMatrix(0,1);
//   SpatialDefRate[4]=totalMatrix(0,2);
//   SpatialDefRate[5]=totalMatrix(1,2);
// }



// template < > 
// void TwoStepUpdatedLagrangianVPSolidElement<2>::CalcDeviatoricInvariant(VectorType &SpatialDefRate,
// 									double &DeviatoricInvariant)
// {
//   double trace_d=SpatialDefRate[0]+SpatialDefRate[1];
//   double dev_X=SpatialDefRate[0]-trace_d/3.0;
//   double dev_Y=SpatialDefRate[1]-trace_d/3.0;
//   DeviatoricInvariant=sqrt(2*(dev_X*dev_X+SpatialDefRate[2]*SpatialDefRate[2]+ dev_Y*dev_Y));

// }


// template < > 
// void TwoStepUpdatedLagrangianVPSolidElement<3>::CalcDeviatoricInvariant(VectorType &SpatialDefRate,
// 									double &DeviatoricInvariant)
// {
//   double trace_d=SpatialDefRate[0]+SpatialDefRate[1]+SpatialDefRate[2];
//   double dev_X=SpatialDefRate[0]-trace_d/3.0;
//   double dev_Y=SpatialDefRate[1]-trace_d/3.0;
//   double dev_Z=SpatialDefRate[2]-trace_d/3.0;
//   DeviatoricInvariant=sqrt(2*(dev_X*dev_X+dev_Y*dev_Y+dev_Z*dev_Z+
// 			      SpatialDefRate[3]*SpatialDefRate[3]+
// 			      SpatialDefRate[4]*SpatialDefRate[4]+
// 			      SpatialDefRate[5]*SpatialDefRate[5]));
// }


// template < > 
// void TwoStepUpdatedLagrangianVPSolidElement<2>::CheckStrain2(MatrixType &SpatialVelocityGrad,
// 							     MatrixType &Fgrad,
// 							     MatrixType &VelDefgrad)
// {
//   if(fabs(VelDefgrad(0,0)-SpatialVelocityGrad(0,0)*Fgrad(0,0)-SpatialVelocityGrad(0,1)*Fgrad(1,0))<0.0000001){
//   }else{
//     std::cout<<"ERROR IN CHECKSTRAIN(2a)";
//   }
//   if(fabs(VelDefgrad(0,1)-SpatialVelocityGrad(0,0)*Fgrad(0,1)-SpatialVelocityGrad(0,1)*Fgrad(1,1))<0.0000001){
//   }else{
//     std::cout<<"ERROR IN CHECKSTRAIN(2b)";
//   }
//   if(fabs(VelDefgrad(1,0)-SpatialVelocityGrad(1,0)*Fgrad(0,0)-SpatialVelocityGrad(1,1)*Fgrad(1,0))<0.0000001){
//   }else{
//     std::cout<<"ERROR IN CHECKSTRAIN(2c)";
//   }
//   if(fabs(VelDefgrad(1,1)-SpatialVelocityGrad(1,0)*Fgrad(0,1)-SpatialVelocityGrad(1,1)*Fgrad(1,1))<0.0000001){
//   }else{
//     std::cout<<"ERROR IN CHECKSTRAIN(2d)";
//   }
// }


// template < > 
// void TwoStepUpdatedLagrangianVPSolidElement<2>::CheckStrain3(VectorType &SpatialDefRate,
// 							     MatrixType &SpatialVelocityGrad)
// {
//   if(fabs(SpatialDefRate[0]-SpatialVelocityGrad(0,0))<0.0000001){
//   }else{
//     std::cout<<"ERROR IN CHECKSTRAIN(3a)";
//   }
//   if(fabs(SpatialDefRate[1]-SpatialVelocityGrad(1,1))<0.0000001){
//   }else{
//     std::cout<<"ERROR IN CHECKSTRAIN(3b)";
//   }
//   if(fabs(SpatialDefRate[2]-0.5*(SpatialVelocityGrad(1,0)+SpatialVelocityGrad(0,1)))<0.0000001){
//   }else{
//     std::cout<<"ERROR IN CHECKSTRAIN(3c)";
//   }
// }

// template < > 
// void TwoStepUpdatedLagrangianVPSolidElement<3>::CheckStrain3(VectorType &SpatialDefRate,
// 							     MatrixType &SpatialVelocityGrad)
// {
//   if(fabs(SpatialDefRate[0]-SpatialVelocityGrad(0,0))<0.0000001){
//    }else{
//     std::cout<<"ERROR IN CHECKSTRAIN(3a)";
//   }
//   if(fabs(SpatialDefRate[1]-SpatialVelocityGrad(1,1))<0.0000001){
//   }else{
//     std::cout<<"ERROR IN CHECKSTRAIN(3b)";
//   }
//   if(fabs(SpatialDefRate[2]-SpatialVelocityGrad(2,2))<0.0000001){
//   }else{
//     std::cout<<"ERROR IN CHECKSTRAIN(3c)";
//   }
//   if(fabs(SpatialDefRate[3]-0.5*(SpatialVelocityGrad(1,0)+SpatialVelocityGrad(0,1)))<0.0000001){
//   }else{
//     std::cout<<"ERROR IN CHECKSTRAIN(3d)";
//   }
//  if(fabs(SpatialDefRate[4]-0.5*(SpatialVelocityGrad(2,0)+SpatialVelocityGrad(0,2)))<0.0000001){
//   }else{
//     std::cout<<"ERROR IN CHECKSTRAIN(3e)";
//   }
//  if(fabs(SpatialDefRate[5]-0.5*(SpatialVelocityGrad(2,1)+SpatialVelocityGrad(1,2)))<0.0000001){
//   }else{
//     std::cout<<"ERROR IN CHECKSTRAIN(3f)";
//   }
// }


template <  unsigned int TDim> 
void TwoStepUpdatedLagrangianVPSolidElement<TDim>:: InitializeElementalVariables(ElementalVariables & rElementalVariables, unsigned int g)
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
void TwoStepUpdatedLagrangianVPSolidElement<2>:: CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables,double TimeStep, unsigned int g)
{

  rElementalVariables.CurrentTotalCauchyStress=this->mCurrentTotalCauchyStress[g];
  rElementalVariables.CurrentDeviatoricCauchyStress=this->mCurrentDeviatoricCauchyStress[g];

  double CurrSecondLame  = 0;
  double CurrBulkModulus = 0;

  this->ComputeMaterialParameters(CurrSecondLame,CurrBulkModulus,TimeStep);
 
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
void TwoStepUpdatedLagrangianVPSolidElement<3>:: CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables, double TimeStep, unsigned int g)
{


  rElementalVariables.CurrentTotalCauchyStress=this->mCurrentTotalCauchyStress[g];
  rElementalVariables.CurrentDeviatoricCauchyStress=this->mCurrentDeviatoricCauchyStress[g];

  double CurrSecondLame  = 0;
  double CurrBulkModulus = 0;

  this->ComputeMaterialParameters(CurrSecondLame,CurrBulkModulus,TimeStep);
 
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


// template < > 
// void TwoStepUpdatedLagrangianVPSolidElement<3>::CalcVolumetricDefRate(const ShapeFunctionDerivativesType& rDN_DX,
// 								      double &volumetricDefRate,
// 								      MatrixType &invGradDef)
// {
//   std::cout<<"TO BE IMPLEMENTED ------- CalcVolumetricDefRate -------"<<std::endl;
//   //you can compute the volumetric deformation rate using CalcVolDefRateFromSpatialVelGrad
// }


// template < > 
// void TwoStepUpdatedLagrangianVPSolidElement<3>::CheckStrain1(double &VolumetricDefRate,
// 							     MatrixType &SpatialVelocityGrad)
// {
//   std::cout<<"TO BE IMPLEMENTED ------- CheckStrain1 -------"<<std::endl;
// }


// template < > 
// void TwoStepUpdatedLagrangianVPSolidElement<3>::CheckStrain2(MatrixType &SpatialVelocityGrad,
// 							     MatrixType &Fgrad,
// 							     MatrixType &VelDefgrad)
// {
//   std::cout<<"TO BE IMPLEMENTED ------- CheckStrain2 -------"<<std::endl;
// }





  template class TwoStepUpdatedLagrangianVPSolidElement<2>;
  template class TwoStepUpdatedLagrangianVPSolidElement<3>;

}
