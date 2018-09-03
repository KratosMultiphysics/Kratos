//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2013 $
//   Revision:            $Revision:                    0.0 $
//
//


// System includes

// External includes

// Project includes
#include "custom_conditions/beam_point_rigid_contact_LM_3D_condition.hpp"

#include "conatct_mechanics_application_variables.h"

namespace Kratos
{
  //************************************************************************************
  //************************************************************************************
  BeamPointRigidContactLM3DCondition::BeamPointRigidContactLM3DCondition(IndexType NewId, GeometryType::Pointer
									   pGeometry)
  : BeamPointRigidContactCondition(NewId, pGeometry)
  {
    //DO NOT ADD DOFS HERE!!!

  }

  //************************************************************************************
  //************************************************************************************
  BeamPointRigidContactLM3DCondition::BeamPointRigidContactLM3DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
  : BeamPointRigidContactCondition(NewId, pGeometry, pProperties)
  {
  }


  //************************************************************************************
  //************************************************************************************
  BeamPointRigidContactLM3DCondition::BeamPointRigidContactLM3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall)
  : BeamPointRigidContactCondition(NewId, pGeometry, pProperties, pRigidWall)
  {

  }

  //************************************************************************************
  //************************************************************************************
  BeamPointRigidContactLM3DCondition::BeamPointRigidContactLM3DCondition( BeamPointRigidContactLM3DCondition const& rOther )
  : BeamPointRigidContactCondition(rOther)
  {
  }

  //************************************************************************************
  //************************************************************************************

  Condition::Pointer BeamPointRigidContactLM3DCondition::Create(IndexType NewId, NodesArrayType
								 const& ThisNodes,  PropertiesType::Pointer pProperties) const
  {
    return Kratos::make_shared<BeamPointRigidContactLM3DCondition>(NewId,GetGeometry().Create(ThisNodes), pProperties);
  }


  //************************************************************************************
  //************************************************************************************


  BeamPointRigidContactLM3DCondition::~BeamPointRigidContactLM3DCondition()
  {

  }

  //************* GETTING METHODS

  //***********************************************************************************
  //***********************************************************************************

  void BeamPointRigidContactLM3DCondition::GetDofList(DofsVectorType& rConditionDofList,
						      ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

      rConditionDofList.resize(0);
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    for (unsigned int i = 0; i < number_of_nodes; i++)
      {
        rConditionDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
        rConditionDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));

	if( dimension == 3 ){
	  rConditionDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
	  rConditionDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
	  rConditionDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Y));
	}

        rConditionDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Z));

        rConditionDofList.push_back(GetGeometry()[i].pGetDof(LAGRANGE_MULTIPLIER_NORMAL));
      }


    KRATOS_CATCH( "" )
  }

  //***********************************************************************************
  //***********************************************************************************

  void BeamPointRigidContactLM3DCondition::EquationIdVector(EquationIdVectorType& rResult,
							    ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int condition_size        = number_of_nodes * ((dimension * (dimension-1))+1);

    if (rResult.size() != condition_size)
      rResult.resize( condition_size, false );

    for (unsigned int i = 0; i < number_of_nodes; i++)
      {
        int index = i * (dimension * (dimension-1));
        rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();

	if( dimension == 3){
	  rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
	  rResult[index + 3] = GetGeometry()[i].GetDof(ROTATION_X).EquationId();
          rResult[index + 4] = GetGeometry()[i].GetDof(ROTATION_Y).EquationId();
	}

	rResult[index + 5] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();

	rResult[index + 6] = GetGeometry()[i].GetDof(LAGRANGE_MULTIPLIER_NORMAL).EquationId();
      }

    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************

  void BeamPointRigidContactLM3DCondition::GetValuesVector(Vector& rValues, int Step)
  {
    KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       condition_size  = number_of_nodes * ((dimension * (dimension-1))+1);

    if ( rValues.size() != condition_size )
      rValues.resize( condition_size, false );

    for (unsigned int i = 0; i < number_of_nodes; i++)
      {
        unsigned int index = i * (dimension * (dimension-1));
        rValues[index]     = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Y, Step );

        if ( dimension == 3 ){
	  rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( DISPLACEMENT_Z, Step );
	  rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( ROTATION_X, Step );
	  rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Y, Step );
	}

	rValues[index + 5] = GetGeometry()[i].GetSolutionStepValue( ROTATION_Z, Step );
	rValues[index + 6] = GetGeometry()[i].GetSolutionStepValue( LAGRANGE_MULTIPLIER_NORMAL, Step );
      }

    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************

  void BeamPointRigidContactLM3DCondition::GetFirstDerivativesVector( Vector& rValues, int Step )
  {
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       condition_size    = number_of_nodes * ((dimension * (dimension-1))+1);

    if ( rValues.size() != condition_size ) rValues.resize( condition_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
        unsigned int index = i * dimension;
        rValues[index]     = GetGeometry()[i].GetSolutionStepValue( VELOCITY_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Y, Step );

        if ( dimension == 3 ){
	  rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( VELOCITY_Z, Step );
	  rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_VELOCITY_X, Step );
	  rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_VELOCITY_Y, Step );
	}

	rValues[index + 5] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_VELOCITY_Z, Step );
	rValues[index + 6] = 0;
      }

    KRATOS_CATCH( "" )
   }


  //***********************************************************************************
  //***********************************************************************************

  void BeamPointRigidContactLM3DCondition::GetSecondDerivativesVector( Vector& rValues, int Step )
  {
    KRATOS_TRY

      const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();
    unsigned int       condition_size    = number_of_nodes * ((dimension * (dimension-1))+1);

    if ( rValues.size() != condition_size ) rValues.resize( condition_size, false );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
        unsigned int index = i * (dimension * (dimension-1));
        rValues[index]     = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_X, Step );
        rValues[index + 1] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Y, Step );

        if ( dimension == 3 ){
	  rValues[index + 2] = GetGeometry()[i].GetSolutionStepValue( ACCELERATION_Z, Step );
	  rValues[index + 3] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_X, Step );
	  rValues[index + 4] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_Y, Step );

	}

	rValues[index + 5] = GetGeometry()[i].GetSolutionStepValue( ANGULAR_ACCELERATION_Z, Step );
	rValues[index + 5] = 0;
      }

    KRATOS_CATCH( "" )
  }


  //************************************************************************************
  //************************************************************************************

  void BeamPointRigidContactLM3DCondition::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
  {
      KRATOS_TRY


       ConditionVariables ContactVariables;

       SpatialBoundingBox::BoundingBoxParameters BoxParameters(this->GetGeometry()[0], ContactVariables.Gap.Normal, ContactVariables.Gap.Tangent, ContactVariables.Surface.Normal, ContactVariables.Surface.Tangent, ContactVariables.RelativeDisplacement);

       //to perform contact with a tube radius must be set
       BoxParameters.SetRadius(GetGeometry()[0].GetValue(CROSS_SECTION_AREA));

       if ( this->mpRigidWall->IsInside( BoxParameters, rCurrentProcessInfo ) ) {

  	   const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
  	   array_1d<double, 3> &ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);

	   double TangentRelativeMovement = 0;
	   TangentRelativeMovement = this->CalculateTangentRelativeMovement( TangentRelativeMovement, ContactVariables );

           mTangentialVariables.PreviousTangentForceModulus = 0.0;
           for (unsigned int i = 0; i < dimension; ++i) {
              mTangentialVariables.PreviousTangentForceModulus += ContactForce[i] * ContactVariables.Surface.Tangent[i];
           }


        }
        else {
           mTangentialVariables.PreviousTangentForceModulus = 0.0;
        }

       mTangentialVariables.DeltaTime = rCurrentProcessInfo[DELTA_TIME];

       mTangentialVariables.Sign = 1;

       mTangentialVariables.FrictionCoefficient        =  GetProperties()[FRICTION_COEFFICIENT];//0.3;
       mTangentialVariables.DynamicFrictionCoefficient =  GetProperties()[MU_DYNAMIC];//0.2;
       mTangentialVariables.StaticFrictionCoefficient  =  GetProperties()[MU_STATIC];//0.3;

       //std::cout<<" Friction Coef "<<mTangentialVariables.FrictionCoefficient<<" dynamic "<<mTangentialVariables.DynamicFrictionCoefficient<<" static "<<mTangentialVariables.StaticFrictionCoefficient<<std::endl;

       ClearNodalForces();


    KRATOS_CATCH( "" )

  }


  //************************************************************************************
  //************************************************************************************

  void BeamPointRigidContactLM3DCondition::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
  {
    //added to control force evolution per step:
    // array_1d<double, 3> &ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);
    // mTangentialVariables.PreviousTangentForceModulus = norm_2(ContactForce);

    ClearNodalForces();
  }

  //************************************************************************************
  //************************************************************************************

  void BeamPointRigidContactLM3DCondition::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
								    VectorType& rRightHandSideVector,
								    Flags& rCalculationFlags)

  {
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * ((dimension * (dimension-1))+1);

    if ( rCalculationFlags.Is(BeamPointRigidContactCondition::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
      {
        if ( rLeftHandSideMatrix.size1() != MatSize )
	  rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
      }


    //resizing as needed the RHS
    if ( rCalculationFlags.Is(BeamPointRigidContactCondition::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
      {
        if ( rRightHandSideVector.size() != MatSize )
	  rRightHandSideVector.resize( MatSize, false );

	rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS

      }
  }

  //************************************************************************************
  //************************************************************************************

  void BeamPointRigidContactLM3DCondition::InitializeConditionVariables (ConditionVariables& rVariables,
									const ProcessInfo& rCurrentProcessInfo)
  {



  }

  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************

  void BeamPointRigidContactLM3DCondition::CalculateKinematics(ConditionVariables& rVariables,
							       const ProcessInfo& rCurrentProcessInfo,
							       const double& rPointNumber)
  {
    KRATOS_TRY

    SpatialBoundingBox::BoundingBoxParameters BoxParameters(this->GetGeometry()[0], rVariables.Gap.Normal, rVariables.Gap.Tangent, rVariables.Surface.Normal, rVariables.Surface.Tangent, rVariables.RelativeDisplacement);

    //to perform contact with a tube radius must be set
    BoxParameters.SetRadius(GetGeometry()[0].GetValue(CROSS_SECTION_AREA));

    if( this->mpRigidWall->IsInside( BoxParameters, rCurrentProcessInfo ) ){

      rVariables.Options.Set(ACTIVE,true);

      //double& NormalLM = GetGeometry()[0].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_NORMAL);

      //check previous normal reaction from the lagrange multiplier
      // if( NormalLM > 0 )
      // 	rVariables.Options.Set(ACTIVE,false);

      //get contact properties and parameters
      this->CalculateContactFactors( rVariables );

    }
    else{

      rVariables.Options.Set(ACTIVE,false);

    }


    KRATOS_CATCH( "" )
      }


  //************************************************************************************
  //************************************************************************************


  void BeamPointRigidContactLM3DCondition::CalculateContactFactors(ConditionVariables &rVariables)
  {

    KRATOS_TRY

    WeakPointerVector<Node<3> >& rN = GetGeometry()[0].GetValue(NEIGHBOUR_NODES);

    array_1d<double,3> Contact_Point = GetGeometry()[0].Coordinates();
    array_1d<double,3> Neighb_Point;

    double distance = 0;
    double counter = 0;

    for(unsigned int i = 0; i < rN.size(); i++)
      {
	if(rN[i].Is(BOUNDARY)){

	  Neighb_Point[0] = rN[i].X();
	  Neighb_Point[1] = rN[i].Y();
	  Neighb_Point[2] = rN[i].Z();

	  distance += norm_2(Contact_Point-Neighb_Point);

	  counter ++;
	}
      }

    if( counter != 0 )
      distance /= counter;

    if( distance == 0 )
      distance = 1;

    //get contact properties and parameters
    double PenaltyParameter = GetProperties()[PENALTY_PARAMETER];
    double ElasticModulus   = GetProperties()[YOUNG_MODULUS];

    double factor = 4;
    if( distance < 1.0 ){ //take a number bigger than 1.0 (length units)
      int order = (int)((-1) * std::log10(distance) + 1) ;
      distance *= factor * pow(10,order);
    }

    //beam reduction
    PenaltyParameter *= 1e-8;

    rVariables.Penalty.Normal  = distance * PenaltyParameter * ElasticModulus;
    rVariables.Penalty.Tangent = rVariables.Penalty.Normal;


    //std::cout<<" Node "<<GetGeometry()[0].Id()<<" Contact Factors "<<rVariables.Penalty.Normal<<" Gap Normal "<<rVariables.Gap.Normal<<" Gap Tangent "<<rVariables.Gap.Tangent<<" Surface.Normal "<<rVariables.Surface.Normal<<" Surface.Tangent "<<rVariables.Surface.Tangent<<" distance "<<distance<<" ElasticModulus "<<ElasticModulus<<" PenaltyParameter "<<PenaltyParameter<<std::endl;

    // std::cout<<" Penalty.Normal "<<rVariables.Penalty.Normal<<" Penalty.Tangent "<<rVariables.Penalty.Tangent<<std::endl;

    KRATOS_CATCH( "" )
  }


  //***********************************************************************************
  //***********************************************************************************

  void BeamPointRigidContactLM3DCondition::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
								ConditionVariables& rVariables,
								double& rIntegrationWeight)

  {
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();


    if( rVariables.Options.Is(ACTIVE)){

      unsigned int size = dimension*2+1;
      MatrixType Kuug = ZeroMatrix(size,size);

      for( unsigned int i = 0; i<dimension+1; i++ )
	{
	  Kuug(i,dimension*2) = -rVariables.Surface.Normal[i];
	  Kuug(dimension*2,i) = -rVariables.Surface.Normal[i];
	}

      //Stabilization term
      double penalty = 1;
      double stab_LM = -rVariables.Gap.Normal * rVariables.Gap.Normal * penalty;

      Kuug(dimension*2, dimension*2) = stab_LM;

      //Building the Local Stiffness Matrix
      BeamMathUtilsType::AddMatrix( rLeftHandSideMatrix, Kuug, 0, 0 );

      // std::cout<<std::endl;
      // std::cout<<" Penalty.Normal "<<rVariables.Penalty.Normal<<" rVariables.Gap.Normal "<<rVariables.Gap.Normal<<" rVariables.Surface.Normal "<<rVariables.Surface.Normal<<" rIntegrationWeight "<<rIntegrationWeight<<" nxn : "<<outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal)<<std::endl;

      // this->CalculateAndAddKuugTangent( rLeftHandSideMatrix,  rVariables, rIntegrationWeight );
      // std::cout<<std::endl;
      // std::cout<<" Kcont "<<rLeftHandSideMatrix<<std::endl;

      // KRATOS_WATCH( rLeftHandSideMatrix )

    }
    else{

      unsigned int size = dimension*2+1;
      rLeftHandSideMatrix= ZeroMatrix(size,size);

      //to avoid system with a zero row and column
      MatrixType Kuug = ZeroMatrix(size,size);

      Kuug(dimension*2,dimension*2) = 1;

      //Building the Local Stiffness Matrix
      BeamMathUtilsType::AddMatrix( rLeftHandSideMatrix, Kuug, 0, 0 );


    }


    KRATOS_CATCH( "" )
  }

  //************* Tangent Contact Force constitutive matrix      **********************
  //***********************************************************************************

   void BeamPointRigidContactLM3DCondition::CalculateAndAddKuugTangent(MatrixType& rLeftHandSideMatrix, ConditionVariables& rVariables, double& rIntegrationWeight)
   {
       KRATOS_TRY


       KRATOS_CATCH( "" )

   }

  //***********************************************************************************
  //***********************************************************************************

  void BeamPointRigidContactLM3DCondition::CalculateAndAddContactForces(VectorType& rRightHandSideVector,
									 ConditionVariables& rVariables,
									 double& rIntegrationWeight)

  {
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    if( rVariables.Options.Is(ACTIVE)){

       this->CalculateAndAddNormalContactForce( rRightHandSideVector, rVariables, rIntegrationWeight );
       //this->CalculateAndAddTangentContactForce( rRightHandSideVector, rVariables, rIntegrationWeight );

       //KRATOS_WATCH( rRightHandSideVector )

    }
    else{

      rRightHandSideVector = ZeroVector(dimension * 2 + 1 );

    }


    KRATOS_CATCH( "" )
  }

  //**************************** Calculate Normal Contact Force ***********************
  //***********************************************************************************

  void BeamPointRigidContactLM3DCondition::CalculateAndAddNormalContactForce(VectorType& rRightHandSideVector,
									      ConditionVariables& rVariables,
									      double& rIntegrationWeight)
  {

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

      double& NormalLM = GetGeometry()[0].FastGetSolutionStepValue(LAGRANGE_MULTIPLIER_NORMAL);

      double penalty = 1;

      rRightHandSideVector[dimension*2] = rVariables.Gap.Normal * (1.0 + NormalLM * rVariables.Gap.Normal * penalty);


      GetGeometry()[0].SetLock();

      array_1d<double, 3 >& ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);


      for(unsigned int j = 0; j < dimension; j++)
	{
	  ContactForce[j] = NormalLM * rVariables.Surface.Normal[j];
	  rRightHandSideVector[j] += ContactForce[j];

	}

      if( NormalLM > 0 )
	ContactForce.clear();

      GetGeometry()[0].UnSetLock();


      //std::cout<<"[ID: "<<GetGeometry()[0].Id()<<"]:  Normal "<<rVariables.Surface.Normal<<" Gap "<<rVariables.Gap.Normal<<" LM_N "<<NormalLM<<std::endl;


  }

  //**************************** Calculate Tangent Contact Force **********************
  //***********************************************************************************

  void BeamPointRigidContactLM3DCondition::CalculateAndAddTangentContactForce(VectorType& rRightHandSideVector,
									       ConditionVariables& rVariables,
									       double& rIntegrationWeight)
  {

      KRATOS_TRY


      KRATOS_CATCH( "" )


  }

  //**************************** Calculate Tangent Force Modulus **********************
  //***********************************************************************************

  double& BeamPointRigidContactLM3DCondition::CalculateTangentRelativeMovement( double& rTangentRelativeMovement, ConditionVariables& rVariables )
  {
       const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

       const array_1d<double, 3> & CurrentDisplacement  =  GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
       const array_1d<double, 3> & PreviousDisplacement =  GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT, 1);
       array_1d<double, 3 > DeltaDisplacement           =  CurrentDisplacement-PreviousDisplacement;


       //Calculate the relative rotation
       const array_1d<double, 3> & CurrentRotation      =  GetGeometry()[0].FastGetSolutionStepValue(ROTATION);
       const array_1d<double, 3> & PreviousRotation     =  GetGeometry()[0].FastGetSolutionStepValue(ROTATION, 1);
       array_1d<double, 3 > DeltaRotation               =  CurrentRotation-PreviousRotation;


       double SectionMeanRadius = GetProperties()[CROSS_SECTION_AREA];
       array_1d<double, 3 > RadiusVector;
       for(unsigned int i = 0; i < 3; i++)
	 {
	   RadiusVector[i] =  SectionMeanRadius * rVariables.Surface.Normal[i];
	 }

       array_1d<double, 3 > DeltaRotationDisplacement;
       MathUtils<double>::CrossProduct( DeltaRotationDisplacement, DeltaRotation, RadiusVector );

       // bool Regularization = false;

       // if( Regularization == true ){

       // 	 //regularization taking the contigous boundary segments
       // 	 WeakPointerVector<Node<3> >& rN = GetGeometry()[0].GetValue(NEIGHBOUR_NODES);

       // 	 array_1d<double, 3 > NeighbDeltaDisplacement;

       // 	 double counter = 0;

       // 	 for(unsigned int i = 0; i < rN.size(); i++)
       // 	   {
       // 	     if(rN[i].Is(BOUNDARY)){

       // 	       const array_1d<double, 3> & NeighbCurrentDisplacement  =  GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
       // 	       const array_1d<double, 3> & NeighbPreviousDisplacement =  GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT, 1);

       // 	       array_1d<double, 3 > NeighbDeltaDisplacement           =  NeighbCurrentDisplacement-NeighbPreviousDisplacement;

       // 	       DeltaDisplacement += NeighbDeltaDisplacement;

       // 	       counter ++;
       // 	     }
       // 	   }
       //        if( counter!= 0)
       // 	    DeltaDisplacement /= counter;

       // }

       VectorType WallDisplacement = mTangentialVariables.DeltaTime * this->mpRigidWall->GetVelocity();

       rTangentRelativeMovement = 0.0;
       VectorType TotalTangentRelativeMovement = ZeroVector(dimension);

       for (unsigned int i = 0; i < dimension; ++i)
	 {
	   TotalTangentRelativeMovement[i] = DeltaDisplacement[i] + DeltaRotationDisplacement[i] - WallDisplacement[i];
	 }

       rTangentRelativeMovement = MathUtils<double>::Norm3(TotalTangentRelativeMovement);

       if( rTangentRelativeMovement !=0 )
	 rVariables.Surface.Tangent = TotalTangentRelativeMovement/rTangentRelativeMovement;


       rVariables.Gap.Tangent = rTangentRelativeMovement;


       return rTangentRelativeMovement;

  }

  //**************************** Check Coulomb law for Tangent Contact Force **********
  //***********************************************************************************

  double BeamPointRigidContactLM3DCondition::CalculateCoulombsFrictionLaw(double & rTangentRelativeMovement, double & rNormalForceModulus , ConditionVariables& rVariables)
  {
       mTangentialVariables.FrictionCoefficient = this->CalculateFrictionCoefficient(rTangentRelativeMovement);


       double TangentForceModulus = rVariables.Penalty.Tangent * rVariables.Gap.Tangent; //+ mTangentialVariables.PreviousTangentForceModulus;



       if ( fabs(TangentForceModulus) >  mTangentialVariables.FrictionCoefficient * fabs(rNormalForceModulus) && fabs(rVariables.Gap.Tangent) > 1e-200) {

	 mTangentialVariables.Sign =  rVariables.Gap.Tangent/ fabs( rVariables.Gap.Tangent ) ;

	 TangentForceModulus =  mTangentialVariables.Sign * mTangentialVariables.FrictionCoefficient * fabs(rNormalForceModulus) ;
	 mTangentialVariables.Slip = true;

       }
       else {
	 mTangentialVariables.Slip = false;
       }


       return TangentForceModulus;
  }



  //**************************** Check friction coefficient ***************************
  //***********************************************************************************

  double BeamPointRigidContactLM3DCondition::CalculateFrictionCoefficient(double & rTangentRelativeMovement)
  {

       //---FRICTION LAW in function of the relative sliding velocity ---//

       double Velocity = 0;
       Velocity = rTangentRelativeMovement / mTangentialVariables.DeltaTime;


       //Addicional constitutive parameter  C
       //which describes how fast the static coefficient approaches the dynamic:
       double C=0.1;

       //Addicional constitutive parameter  E
       //regularization parameter (->0, classical Coulomb law)
       double E=0.01;


       double FrictionCoefficient = mTangentialVariables.DynamicFrictionCoefficient + ( mTangentialVariables.StaticFrictionCoefficient-  mTangentialVariables.DynamicFrictionCoefficient ) * exp( (-1) * C * fabs(Velocity) );


       //Square root regularization
       FrictionCoefficient *= fabs(Velocity)/sqrt( ( Velocity * Velocity ) + ( E * E ) );

       //Hyperbolic regularization
       //FrictionCoefficient *= tanh( fabs(Velocity)/E );

       return FrictionCoefficient;

  }


} // Namespace Kratos
