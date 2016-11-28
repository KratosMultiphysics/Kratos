//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_conditions/point_rigid_contact_penalty_3D_condition.hpp"

#include "contact_mechanics_application_variables.h"

#include "custom_friction/friction_law.hpp"

namespace Kratos
{
  
  //************************************************************************************
  //************************************************************************************
  
  PointRigidContactPenalty3DCondition::PointRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : PointRigidContactCondition(NewId, pGeometry)
  {
    //DO NOT ADD DOFS HERE!!!
  }

  //************************************************************************************
  //************************************************************************************
  
  PointRigidContactPenalty3DCondition::PointRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : PointRigidContactCondition(NewId, pGeometry, pProperties)
  {
    //DO NOT ADD DOFS HERE!!!    
  }


  //************************************************************************************
  //************************************************************************************
  
  PointRigidContactPenalty3DCondition::PointRigidContactPenalty3DCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall)
    : PointRigidContactCondition(NewId, pGeometry, pProperties, pRigidWall)
  {
    //DO NOT ADD DOFS HERE!!!
  }

  //************************************************************************************
  //************************************************************************************
  
  PointRigidContactPenalty3DCondition::PointRigidContactPenalty3DCondition( PointRigidContactPenalty3DCondition const& rOther )
    : PointRigidContactCondition(rOther)
  {
    //DO NOT ADD DOFS HERE!!! 
  }

  //************************************************************************************
  //************************************************************************************

  Condition::Pointer PointRigidContactPenalty3DCondition::Create(IndexType NewId, const NodesArrayType& ThisNodes, PropertiesType::Pointer pProperties) const
  {
    return Condition::Pointer(new PointRigidContactPenalty3DCondition(NewId,GetGeometry().Create(ThisNodes), pProperties));
  }

  //************************************CLONE*******************************************
  //************************************************************************************

  Condition::Pointer PointRigidContactPenalty3DCondition::Clone(IndexType NewId, const NodesArrayType& ThisNodes) const
  {
    return Condition::Pointer(new PointRigidContactPenalty3DCondition(NewId,GetGeometry().Create(ThisNodes), pGetProperties(), mpRigidWall));
  }


  //************************************************************************************
  //************************************************************************************

  PointRigidContactPenalty3DCondition::~PointRigidContactPenalty3DCondition()
  {

  }

  //************************************************************************************
  //************************************************************************************

  void PointRigidContactPenalty3DCondition::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    GeneralVariables ContactVariables;

    SpatialBoundingBox::BoundingBoxParameters BoxParameters(this->GetGeometry()[0], ContactVariables.Gap.Normal, ContactVariables.Gap.Tangent, ContactVariables.Surface.Normal, ContactVariables.Surface.Tangent, ContactVariables.RelativeDisplacement);           

    if ( this->mpRigidWall->IsInside( BoxParameters, rCurrentProcessInfo) ) {

      const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
      array_1d<double, 3> &ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);

      mTangentialVariables.PreviousTangentForceModulus = 0.0;
      for (unsigned int i = 0; i < dimension; ++i) {
    	mTangentialVariables.PreviousTangentForceModulus += ContactForce[i] * ContactVariables.Surface.Tangent[i];
      }

    }
    else {
      
      mTangentialVariables.PreviousTangentForceModulus = 0.0;
    } 

    mTangentialVariables.PreviousTangentForceModulus = 0.0;
    
    mTangentialVariables.DeltaTime = rCurrentProcessInfo[DELTA_TIME];
    mTangentialVariables.Sign = 1;

    mTangentialVariables.DynamicFrictionCoefficient = 0.0;//0.2;
    mTangentialVariables.StaticFrictionCoefficient  = 0.0;//0.3;

    if( GetProperties().Has(FRICTION_ACTIVE) ){
      if( GetProperties()[FRICTION_ACTIVE] ){

	if( GetProperties().Has(MU_DYNAMIC) )
	  mTangentialVariables.DynamicFrictionCoefficient = GetProperties()[MU_DYNAMIC];

	if( GetProperties().Has(MU_STATIC) )
	  mTangentialVariables.StaticFrictionCoefficient = GetProperties()[MU_STATIC];
      }
    }
	  
    mTangentialVariables.FrictionCoefficient = mTangentialVariables.StaticFrictionCoefficient;      

    // Compute the neighbour distance, then a stress-"like" may be computed.
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

    mTangentialVariables.Neighb_distance = distance; 

    //when implex is active --> it deletes contact forces at last implex step<--
    //std::cout<<" Initialize Solution Step ["<<this->Id()<<"]"<<std::endl;
    //ClearNodalForces();

    KRATOS_CATCH( "" )

  }


  //************************************************************************************
  //************************************************************************************
  void PointRigidContactPenalty3DCondition::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY
      
    rCurrentProcessInfo[NUMBER_OF_ACTIVE_CONTACTS] = 0;
    rCurrentProcessInfo[NUMBER_OF_STICK_CONTACTS]  = 0;
    rCurrentProcessInfo[NUMBER_OF_SLIP_CONTACTS]   = 0;
   
    ClearNodalForces();

    KRATOS_CATCH( "" )
    
  }

  //************************************************************************************
  //************************************************************************************

  void PointRigidContactPenalty3DCondition::InitializeGeneralVariables (GeneralVariables& rVariables, 
									const ProcessInfo& rCurrentProcessInfo)
  {
    KRATOS_TRY

    KRATOS_CATCH( "" )
  }

  //*********************************COMPUTE KINEMATICS*********************************
  //************************************************************************************

  void PointRigidContactPenalty3DCondition::CalculateKinematics(GeneralVariables& rVariables, const ProcessInfo& rCurrentProcessInfo, const double& rPointNumber)
  {
    KRATOS_TRY

    SpatialBoundingBox::BoundingBoxParameters BoxParameters(this->GetGeometry()[0], rVariables.Gap.Normal, rVariables.Gap.Tangent, rVariables.Surface.Normal, rVariables.Surface.Tangent, rVariables.RelativeDisplacement);

    if( this->mpRigidWall->IsInside( BoxParameters, rCurrentProcessInfo ) ){

      rVariables.Options.Set(ACTIVE,true);

      rVariables.Gap.Normal = fabs(rVariables.Gap.Normal);

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

  void PointRigidContactPenalty3DCondition::CalculateContactFactors(GeneralVariables &rVariables)
  {

    KRATOS_TRY

    WeakPointerVector<Node<3> >& rN = GetGeometry()[0].GetValue(NEIGHBOUR_NODES);

    array_1d<double,3> Contact_Point = GetGeometry()[0].Coordinates();
    array_1d<double,3> Neighb_Point;

    double distance = 0;
    double counter  = 0;

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
    double PenaltyParameter = 1;
    if( GetProperties().Has(PENALTY_PARAMETER) )
      PenaltyParameter = GetProperties()[PENALTY_PARAMETER];

    WeakPointerVector<Element >& rE = GetGeometry()[0].GetValue(NEIGHBOUR_ELEMENTS);
    double ElasticModulus = 0;
    if( GetProperties().Has(YOUNG_MODULUS) )
      ElasticModulus = GetProperties()[YOUNG_MODULUS];
    else
      ElasticModulus = rE.front().GetProperties()[YOUNG_MODULUS];
      
    // the Modified Cam Clay model does not have a constant Young modulus, so something similar to that is computed
    if (ElasticModulus <= 1.0e-5) {
      std::vector<double> mModulus;
      ProcessInfo SomeProcessInfo;
      for ( unsigned int i = 0; i < rE.size(); i++)
	{
	  rE[i].CalculateOnIntegrationPoints(EQUIVALENT_YOUNG_MODULUS, mModulus, SomeProcessInfo);
	  ElasticModulus += mModulus[0];
	}
      ElasticModulus /= double(rE.size());
    }

    double factor = 4;
    if( distance < 1.0 ){ //take a number bigger than 1.0 (length units)
      int order = (int)((-1) * std::log10(distance) + 1) ;
      distance *= factor * pow(10,order);
    }

    rVariables.Penalty.Normal  = distance * PenaltyParameter * ElasticModulus;

    double PenaltyRatio = 1;
    if( GetProperties().Has(TANGENTIAL_PENALTY_RATIO) )
      PenaltyRatio = GetProperties()[TANGENTIAL_PENALTY_RATIO];

    rVariables.Penalty.Tangent = rVariables.Penalty.Normal * PenaltyRatio ;  

    //set contact normal
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
      {
	GetGeometry()[i].SetLock();
	 
	array_1d<double, 3> &ContactNormal  = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_NORMAL);

	for(unsigned int i=0; i<3; i++)
	  ContactNormal[i] = rVariables.Surface.Normal[i];	 

	GetGeometry()[i].UnSetLock();
      }
      
    KRATOS_CATCH( "" )

 }


  //***********************************************************************************
  //***********************************************************************************

  void PointRigidContactPenalty3DCondition::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
								GeneralVariables& rVariables,
								double& rIntegrationWeight)

  {
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if( rVariables.Options.Is(ACTIVE)){

      noalias(rLeftHandSideMatrix) = rVariables.Penalty.Normal * rIntegrationWeight  * outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal);

      this->CalculateAndAddKuugTangent( rLeftHandSideMatrix,  rVariables, rIntegrationWeight);

    }
    else{

      rLeftHandSideMatrix= ZeroMatrix(dimension,dimension);
    }

    KRATOS_CATCH( "" )
  }

  //************* Tangent Contact Force constitutive matrix      **********************
  //***********************************************************************************

  void PointRigidContactPenalty3DCondition::CalculateAndAddKuugTangent(MatrixType& rLeftHandSideMatrix, GeneralVariables& rVariables, double& rIntegrationWeight)
  {
    KRATOS_TRY
      
    //const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    
    double NormalForceModulus = 0;
    NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables );
    double TangentForceModulus = this->CalculateCoulombsFrictionLaw( rVariables.Gap.Tangent, NormalForceModulus, rVariables );
    //std::cout<<" TangentForceModulus "<<TangentForceModulus<<std::endl;
   
    if( fabs(TangentForceModulus) >= 1e-25 ){
      
      if( mTangentialVariables.Slip ){
	//simpler expression:
	// Matrix Kuug = rLeftHandSideMatrix;	
	// std::cout<<" Kuug "<<Kuug<<std::endl;
	
	noalias(rLeftHandSideMatrix) += mTangentialVariables.FrictionCoefficient * rVariables.Penalty.Normal * rIntegrationWeight * ( outer_prod(rVariables.Surface.Tangent, rVariables.Surface.Normal) );
	
	noalias(rLeftHandSideMatrix) += mTangentialVariables.FrictionCoefficient * rVariables.Penalty.Normal * rIntegrationWeight * (rVariables.Gap.Normal/rVariables.Gap.Tangent) * ( IdentityMatrix(3,3) - outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal) );

	//std::cout<<" A:Kuug "<<rLeftHandSideMatrix<<std::endl;
	
      }
      else {


	noalias(rLeftHandSideMatrix) += rVariables.Penalty.Tangent * rIntegrationWeight * outer_prod(rVariables.Surface.Tangent, rVariables.Surface.Tangent);
	
	noalias(rLeftHandSideMatrix) += rVariables.Penalty.Tangent * rIntegrationWeight * ( IdentityMatrix(3,3) - outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal) );

	//std::cout<<" B:Kuug "<<rLeftHandSideMatrix<<std::endl;
      }

    }
    
    KRATOS_CATCH( "" )
  }

  void PointRigidContactPenalty3DCondition::CalculateAndAddContactForces(VectorType& rRightHandSideVector,
									 GeneralVariables& rVariables,
									 double& rIntegrationWeight)
  {
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if( rVariables.Options.Is(ACTIVE)){

      mTangentialVariables.IntegrationWeight  = rIntegrationWeight;

      this->CalculateAndAddNormalContactForce( rRightHandSideVector, rVariables, rIntegrationWeight );
      this->CalculateAndAddTangentContactForce( rRightHandSideVector, rVariables, rIntegrationWeight );

    }
    else{

      rRightHandSideVector = ZeroVector(dimension);

    }


    KRATOS_CATCH( "" )
  }

  //**************************** Calculate Normal Contact Force ***********************
  //***********************************************************************************

  void PointRigidContactPenalty3DCondition::CalculateAndAddNormalContactForce(VectorType& rRightHandSideVector,
									      GeneralVariables& rVariables,
									      double& rIntegrationWeight)
  {

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    double NormalForceModulus = 0;
    NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables );

    NormalForceModulus *= rIntegrationWeight;

    for(unsigned int j = 0; j < dimension; j++)
      {
	rRightHandSideVector[j] = NormalForceModulus * rVariables.Surface.Normal[j];
      }


    GetGeometry()[0].SetLock();

    array_1d<double, 3 >& ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);


    for(unsigned int j = 0; j < dimension; j++)
      {
	ContactForce[j] = rRightHandSideVector[j];
      }
    

    if ( GetGeometry()[0].SolutionStepsDataHas( EFFECTIVE_CONTACT_FORCE )) { 

      double EffectiveNormalModulus = CalculateEffectiveNormalForceModulus( NormalForceModulus*(-1) / rIntegrationWeight);
      EffectiveNormalModulus *= (-1) * rIntegrationWeight;
 
      array_1d<double, 3 >& EffectiveContactForce = GetGeometry()[0].FastGetSolutionStepValue(EFFECTIVE_CONTACT_FORCE);
      for (unsigned int j = 0; j < dimension; j++) {
	EffectiveContactForce[j] = EffectiveNormalModulus * rVariables.Surface.Normal[j];
      }
    }

    GetGeometry()[0].UnSetLock();

  }

  //**************************** Calculate Tangent Contact Force **********************
  //***********************************************************************************

  void PointRigidContactPenalty3DCondition::CalculateAndAddTangentContactForce(VectorType& rRightHandSideVector,
									       GeneralVariables& rVariables,
									       double& rIntegrationWeight)
  {

    mTangentialVariables.IntegrationWeight = rIntegrationWeight;
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    double NormalForceModulus = 0;
    NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables ); 

    double TangentForceModulus =  this->CalculateCoulombsFrictionLaw( rVariables.Gap.Tangent, NormalForceModulus, rVariables );

    TangentForceModulus *= rIntegrationWeight;

    GetGeometry()[0].SetLock();

    array_1d<double, 3 > & ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);

    for (unsigned int i = 0; i < dimension ; ++i) {
      rRightHandSideVector[i] += TangentForceModulus * rVariables.Surface.Tangent[i];
      ContactForce[i] += TangentForceModulus * rVariables.Surface.Tangent[i];
    }


    if ( GetGeometry()[0].SolutionStepsDataHas( EFFECTIVE_CONTACT_FORCE )) { 
      array_1d<double, 3 > & EffectiveContactForce = GetGeometry()[0].FastGetSolutionStepValue( EFFECTIVE_CONTACT_FORCE );
      for (unsigned int i = 0; i < dimension ; ++i) {
	EffectiveContactForce[i] += TangentForceModulus * rVariables.Surface.Tangent[i];
      }
    }

    GetGeometry()[0].UnSetLock();

  }

  //**************************** Calculate Normal Force Modulus ***********************
  //***********************************************************************************

  double& PointRigidContactPenalty3DCondition::CalculateNormalForceModulus ( double& rNormalForceModulus, GeneralVariables& rVariables )
  {

    rNormalForceModulus = (rVariables.Penalty.Normal * rVariables.Gap.Normal);

    return rNormalForceModulus;

  }



  //**************************** Check Coulomb law for Tangent Contact Force **********
  //***********************************************************************************

  double PointRigidContactPenalty3DCondition::CalculateCoulombsFrictionLaw(double & rTangentRelativeMovement, double & rNormalForceModulus , GeneralVariables& rVariables)
  {

    mTangentialVariables.FrictionCoefficient = this->CalculateFrictionCoefficient(rTangentRelativeMovement);

    double FrictionalForceModulus      =  rVariables.Penalty.Tangent * rVariables.Gap.Tangent;      
    double SlipFrictionalForceModulus  =  mTangentialVariables.FrictionCoefficient * fabs(rNormalForceModulus);
       
    if( fabs(FrictionalForceModulus) > fabs(SlipFrictionalForceModulus) ){

      mTangentialVariables.Slip = true;
	 
      FrictionalForceModulus = SlipFrictionalForceModulus;

    }
    else {
      mTangentialVariables.Slip = false;
	 
    }

    return FrictionalForceModulus;
       
  }



  //**************************** Check friction coefficient ***************************
  //***********************************************************************************

  double PointRigidContactPenalty3DCondition::CalculateEffectiveNormalForceModulus(const double& rNormalForceModulus)
  {

    double EffectiveForce = rNormalForceModulus;

    if ( GetGeometry()[0].HasDofFor(WATER_PRESSURE) )
      {
	double WaterForce = GetGeometry()[0].FastGetSolutionStepValue( WATER_PRESSURE );
	if (WaterForce > 0.0)
	  WaterForce = 0.0;
	WaterForce *= mTangentialVariables.Neighb_distance;

	EffectiveForce = rNormalForceModulus - WaterForce ; // due to the sign convention
      }

    return EffectiveForce; 
  }


  double PointRigidContactPenalty3DCondition::CalculateFrictionCoefficient(double & rTangentRelativeMovement)
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


    double FrictionCoefficient = mTangentialVariables.DynamicFrictionCoefficient + ( mTangentialVariables.StaticFrictionCoefficient - mTangentialVariables.DynamicFrictionCoefficient ) * exp( (-1) * C * fabs(Velocity) );


    //Square root regularization
    FrictionCoefficient *= fabs(Velocity)/sqrt( ( Velocity * Velocity ) + ( E * E ) );

    //Hyperbolic regularization
    //FrictionCoefficient *= tanh( fabs(Velocity)/E );

    return FrictionCoefficient;

  }


} // Namespace Kratos



