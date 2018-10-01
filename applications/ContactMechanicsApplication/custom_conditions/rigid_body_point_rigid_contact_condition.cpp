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
#include "custom_conditions/rigid_body_point_rigid_contact_condition.hpp"

#include "contact_mechanics_application_variables.h"

namespace Kratos
{

//***********************************************************************************
//***********************************************************************************
RigidBodyPointRigidContactCondition::RigidBodyPointRigidContactCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : PointRigidContactCondition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//***********************************************************************************
//***********************************************************************************
RigidBodyPointRigidContactCondition::RigidBodyPointRigidContactCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : PointRigidContactCondition(NewId, pGeometry, pProperties)
{

    //DO NOT ADD DOFS HERE!!!
}


//************************************************************************************
//************************************************************************************
RigidBodyPointRigidContactCondition::RigidBodyPointRigidContactCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties, SpatialBoundingBox::Pointer pRigidWall)
  : PointRigidContactCondition(NewId, pGeometry, pProperties)
{
    mpRigidWall = pRigidWall;

    const unsigned int inode = GetGeometry().PointsNumber()-1;

    mMasterElements = GetGeometry()[inode].GetValue(MASTER_ELEMENTS);

    //DO NOT ADD DOFS HERE!!!
}


//************************************************************************************
//************************************************************************************
RigidBodyPointRigidContactCondition::RigidBodyPointRigidContactCondition( RigidBodyPointRigidContactCondition const& rOther )
    : PointRigidContactCondition(rOther)
    , mMasterElements(rOther.mMasterElements)

{
}

//***********************************************************************************
//***********************************************************************************
Condition::Pointer RigidBodyPointRigidContactCondition::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
  return Kratos::make_shared<RigidBodyPointRigidContactCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}


//************************************CLONE*******************************************
//************************************************************************************

Condition::Pointer RigidBodyPointRigidContactCondition::Clone( IndexType NewId, NodesArrayType const& rThisNodes ) const
{
  return this->Create( NewId, rThisNodes, pGetProperties() );
}


//***********************************************************************************
//***********************************************************************************
RigidBodyPointRigidContactCondition::~RigidBodyPointRigidContactCondition()
{
}

//************* GETTING METHODS

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::GetDofList(DofsVectorType& rConditionDofList,
				    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rConditionDofList.resize(0);

    Element& MasterElement = mMasterElements.back();
    MasterElement.GetDofList(rConditionDofList, rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::EquationIdVector(EquationIdVectorType& rResult,
					  ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    rResult.resize( 0, false );

    Element& MasterElement = mMasterElements.back();
    MasterElement.EquationIdVector(rResult, rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::GetValuesVector(Vector& rValues, int Step)
{
    KRATOS_TRY

    Element& MasterElement = mMasterElements.back();
    MasterElement.GetValuesVector(rValues, Step);


    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::GetFirstDerivativesVector( Vector& rValues, int Step )
{
    KRATOS_TRY

    Element& MasterElement = mMasterElements.back();
    MasterElement.GetFirstDerivativesVector(rValues, Step);

    KRATOS_CATCH( "" )
}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::GetSecondDerivativesVector( Vector& rValues, int Step )
{
    KRATOS_TRY

    Element& MasterElement = mMasterElements.back();
    MasterElement.GetSecondDerivativesVector(rValues, Step);

    KRATOS_CATCH( "" )
}

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::AddExplicitContribution(const VectorType& rRHSVector,
						 const Variable<VectorType>& rRHSVariable,
						 Variable<array_1d<double,3> >& rDestinationVariable,
						 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    if( rRHSVariable == CONTACT_FORCES_VECTOR && rDestinationVariable == CONTACT_FORCE )
      {

	for(unsigned int i=0; i< number_of_nodes; i++)
	  {
	    int index = i * (dimension * (dimension-1));

	    GetGeometry()[i].SetLock();

	    array_1d<double, 3 > &ContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE);
	    for(unsigned int j=0; j<dimension; j++)
	      {
		ContactForce[j] += rRHSVector[index + j];
	      }

	    GetGeometry()[i].UnSetLock();
	  }
      }


    if( rRHSVariable == RESIDUAL_VECTOR && rDestinationVariable == FORCE_RESIDUAL )
      {

	for(unsigned int i=0; i< number_of_nodes; i++)
	  {
	    int index = i * (dimension * (dimension-1));

	    GetGeometry()[i].SetLock();

	    array_1d<double, 3 > &ForceResidual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
	    for(unsigned int j=0; j<dimension; j++)
	      {
		ForceResidual[j] += rRHSVector[index + j];
	      }

	    GetGeometry()[i].UnSetLock();
	  }
      }

    KRATOS_CATCH( "" )
}

//************* STARTING - ENDING  METHODS
//***********************************************************************************
//***********************************************************************************

//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
								   VectorType& rRightHandSideVector,
								   Flags& rCalculationFlags)

{
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * (dimension * 2);

    if ( rCalculationFlags.Is(RigidBodyPointRigidContactCondition::COMPUTE_LHS_MATRIX) ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }


    //resizing as needed the RHS
    if ( rCalculationFlags.Is(RigidBodyPointRigidContactCondition::COMPUTE_RHS_VECTOR) ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != MatSize )
	    rRightHandSideVector.resize( MatSize, false );

	rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS

    }
}



//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************

void RigidBodyPointRigidContactCondition::CalculateKinematics(ConditionVariables& rVariables, const ProcessInfo& rCurrentProcessInfo, const double& rPointNumber)
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

    rVariables.DeltaTime = rCurrentProcessInfo[DELTA_TIME];
    KRATOS_CATCH( "" )
  }

  //************************************************************************************
  //************************************************************************************


  void RigidBodyPointRigidContactCondition::CalculateContactFactors(ConditionVariables &rVariables)
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

    rVariables.ContributoryFactor = distance;

    //get contact properties and parameters
    double PenaltyParameter = GetProperties()[PENALTY_PARAMETER];
    double ElasticModulus   = mMasterElements.front().GetProperties()[YOUNG_MODULUS];

    //reduction of the penalty parameter:
    PenaltyParameter *=1e-6;

    double factor = 4;
    if( distance < 1.0 ){ //take a number bigger than 1.0 (length units)
      int order = (int)((-1) * std::log10(distance) + 1) ;
      distance *= factor * pow(10,order);
    }

    rVariables.Penalty.Normal  = distance * PenaltyParameter * ElasticModulus;
    rVariables.Penalty.Tangent = rVariables.Penalty.Normal;


    PointType CentroidPosition = GetGeometry()[0].Coordinates() - mMasterElements.front().GetGeometry()[0].Coordinates();


    //compute the skewsymmmetric tensor of the distance
    this->VectorToSkewSymmetricTensor(CentroidPosition, rVariables.SkewSymDistance);


    //std::cout<<" Node "<<GetGeometry()[0].Id()<<" Contact Factors "<<rVariables.Penalty.Normal<<" Gap Normal "<<rVariables.Gap.Normal<<" Gap Tangent "<<rVariables.Gap.Tangent<<" Surface.Normal "<<rVariables.Surface.Normal<<" Surface.Tangent "<<rVariables.Surface.Tangent<<" distance "<<distance<<" ElasticModulus "<<ElasticModulus<<" PenaltyParameter "<<PenaltyParameter<<std::endl;

    // std::cout<<" Penalty.Normal "<<rVariables.Penalty.Normal<<" Penalty.Tangent "<<rVariables.Penalty.Tangent<<std::endl;

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




//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


void RigidBodyPointRigidContactCondition::CalculateAndAddKuug(MatrixType& rLeftHandSideMatrix,
							      ConditionVariables& rVariables,
							      double& rIntegrationWeight)

{
    KRATOS_TRY


    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();


    if( rVariables.Options.Is(ACTIVE)){

      //Force
      Matrix ForceMatrix(3,3);
      noalias(ForceMatrix) = ZeroMatrix(3,3);

      noalias(ForceMatrix) = rVariables.Penalty.Normal * rIntegrationWeight  * outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal);


      for(unsigned int i=0; i<dimension; i++)
	{
	  for(unsigned int j=0; j<dimension; j++)
	    {
	      rLeftHandSideMatrix(i,j) += ForceMatrix(i,j);
	    }
	}


      //Moment
      Matrix MomentMatrix(3,3);
      noalias(MomentMatrix) = ZeroMatrix(3,3);

      noalias(MomentMatrix) = prod(ForceMatrix,rVariables.SkewSymDistance);

      for(unsigned int i=0; i<dimension; i++)
	{
	  for(unsigned int j=0; j<dimension; j++)
	    {
	      rLeftHandSideMatrix(i+dimension,j) -= MomentMatrix(i,j);
	    }
	}



      // std::cout<<std::endl;
      // std::cout<<" Penalty.Normal "<<rVariables.Penalty.Normal<<" rVariables.Gap.Normal "<<rVariables.Gap.Normal<<" rVariables.Surface.Normal "<<rVariables.Surface.Normal<<" rVariables.Surface.Tangent "<<rVariables.Surface.Tangent<<" rIntegrationWeight "<<rIntegrationWeight<<" nxn : "<<outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal)<<std::endl;

      this->CalculateAndAddKuugTangent( rLeftHandSideMatrix,  rVariables, rIntegrationWeight);
      // std::cout<<std::endl;
      // std::cout<<" Kcont "<<rLeftHandSideMatrix<<std::endl;

    }
    else{

      rLeftHandSideMatrix= ZeroMatrix(dimension*2,dimension*2);

    }

    //KRATOS_WATCH( rLeftHandSideMatrix )

    KRATOS_CATCH( "" )
}


//************* Tangent Contact Force constitutive matrix      **********************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::CalculateAndAddKuugTangent(MatrixType& rLeftHandSideMatrix, ConditionVariables& rVariables, double& rIntegrationWeight)
{

  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();


  double NormalForceModulus = 0;
  NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables );

  double TangentForceModulus = this->CalculateCoulombsFrictionLaw( rVariables.Gap.Tangent, NormalForceModulus, rVariables );

  //Force
  Matrix ForceMatrix(3,3);
  noalias(ForceMatrix) = ZeroMatrix(3,3);


  if( fabs(TangentForceModulus) >= 1e-25 ){

    if ( rVariables.Slip ) {
      //simpler expression:
      noalias(ForceMatrix) = rVariables.FrictionCoefficient * rVariables.Penalty.Normal * rIntegrationWeight * ( outer_prod(rVariables.Surface.Tangent, rVariables.Surface.Normal) );

      noalias(ForceMatrix) += rVariables.FrictionCoefficient * rVariables.Penalty.Normal * rIntegrationWeight * ( outer_prod(rVariables.Surface.Tangent, rVariables.Surface.Normal) + (rVariables.Gap.Normal/rVariables.Gap.Tangent) *( IdentityMatrix(3,3) - outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal) ));

      //extra term (2D)
      //if( dimension == 2 )
	//noalias(ForceMatrix) -= rVariables.FrictionCoefficient * rVariables.Penalty.Normal * rIntegrationWeight * (rVariables.Gap.Normal/rVariables.Gap.Tangent) * (outer_prod(rVariables.Surface.Tangent, VectorType( rVariables.Surface.Tangent - ( inner_prod(rVariables.RelativeDisplacement,rVariables.Surface.Normal) * rVariables.Surface.Tangent ) - ( inner_prod(rVariables.Surface.Normal,rVariables.RelativeDisplacement) * rVariables.Surface.Normal) ) ) );


    }
    else {

      noalias(ForceMatrix) = rVariables.Penalty.Tangent * rIntegrationWeight * outer_prod(rVariables.Surface.Tangent, rVariables.Surface.Tangent);

      noalias(ForceMatrix) += rVariables.Penalty.Tangent * rIntegrationWeight * ( IdentityMatrix(3,3) - outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal) );

      //extra term (2D)
      //if( dimension == 2 )
        //noalias(ForceMatrix) -= rVariables.Penalty.Tangent * rIntegrationWeight * (outer_prod(rVariables.Surface.Tangent, VectorType( rVariables.Surface.Tangent - ( inner_prod(rVariables.RelativeDisplacement,rVariables.Surface.Normal) * rVariables.Surface.Tangent ) - ( inner_prod(rVariables.Surface.Normal,rVariables.RelativeDisplacement) * rVariables.Surface.Normal) ) ) );

    }

  }


  //noalias(ForceMatrix) = rVariables.Penalty.Normal * rIntegrationWeight  * outer_prod(rVariables.Surface.Normal, rVariables.Surface.Normal);


  for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
	{
	  rLeftHandSideMatrix(i,j) += ForceMatrix(i,j);
	}
    }

  // std::cout<<" KuuT "<<ForceMatrix<<std::endl;

  //Moment
  Matrix MomentMatrix(3,3);
  noalias(MomentMatrix) = ZeroMatrix(3,3);

  MomentMatrix = prod(ForceMatrix,rVariables.SkewSymDistance);

  for(unsigned int i=0; i<dimension; i++)
    {
      for(unsigned int j=0; j<dimension; j++)
	{
	  rLeftHandSideMatrix(i+dimension,j) -= MomentMatrix(i,j);
	}
    }


}


//***********************************************************************************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::CalculateAndAddContactForces(VectorType& rRightHandSideVector,
							      ConditionVariables& rVariables,
							      double& rIntegrationWeight)

{
    KRATOS_TRY

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if( rVariables.Options.Is(ACTIVE)){

       this->CalculateAndAddNormalContactForce( rRightHandSideVector, rVariables, rIntegrationWeight );
       this->CalculateAndAddTangentContactForce( rRightHandSideVector, rVariables, rIntegrationWeight );

    }
    else{

      rRightHandSideVector = ZeroVector(dimension*2);

    }

    KRATOS_CATCH( "" )
}


//**************************** Calculate Normal Contact Force ***********************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::CalculateAndAddNormalContactForce(VectorType& rRightHandSideVector,
									    ConditionVariables& rVariables,
									    double& rIntegrationWeight)
{

  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();


  double NormalForceModulus = 0;
  NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables );

  NormalForceModulus *= rIntegrationWeight;

  VectorType ContactForceVector = ZeroVector(3);

  for (unsigned int i = 0; i < dimension ; ++i) {
    ContactForceVector[i]    = NormalForceModulus * rVariables.Surface.Normal[i];
    rRightHandSideVector[i] += ContactForceVector[i];
  }

  GetGeometry()[0].SetLock();

  array_1d<double, 3 >& ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);


  for(unsigned int j = 0; j < dimension; j++)
    {
      ContactForce[j] += NormalForceModulus * rVariables.Surface.Normal[j];
    }

  GetGeometry()[0].UnSetLock();

  VectorType ContactTorque = ZeroVector(3);

  //ContactTorque = MathUtils<double>::CrossProduct( rVariables.CentroidPosition, ContactForceVector);
  //std::cout<<" [ContactTorqueA]: "<<ContactTorque;

  ContactTorque = prod(rVariables.SkewSymDistance,ContactForceVector); //  = (D x F)

  // std::cout<<" [ContactTorqueB]: "<<ContactTorque;
  // std::cout<<" [ContactForce]: "<<ContactForceVector;
  // std::cout<<" [Normal]: "<<rVariables.Surface.Normal;
  // std::cout<<" [Distance]: "<<rVariables.SkewSymDistance;
  // std::cout<<std::endl;

  //Contact torque due to contact force on beam surface
  for (unsigned int i =0; i < dimension; ++i) {
    rRightHandSideVector[i+dimension] += ContactTorque[i];
  }


  //std::cout<<" Fcont "<<rRightHandSideVector<<std::endl;

}

//**************************** Calculate Tangent Contact Force **********************
//***********************************************************************************

void RigidBodyPointRigidContactCondition::CalculateAndAddTangentContactForce(VectorType& rRightHandSideVector,
									     ConditionVariables& rVariables,
									     double& rIntegrationWeight)
{

  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

  double NormalForceModulus = 0;
  NormalForceModulus = this->CalculateNormalForceModulus( NormalForceModulus, rVariables );

  double TangentForceModulus =  this->CalculateCoulombsFrictionLaw( rVariables.Gap.Tangent , NormalForceModulus, rVariables );

  TangentForceModulus *= rIntegrationWeight;


  VectorType ContactForceVector = ZeroVector(3);

  for (unsigned int i = 0; i < dimension ; ++i) {
    ContactForceVector[i]    = TangentForceModulus * rVariables.Surface.Tangent[i];
    rRightHandSideVector[i] += ContactForceVector[i];
  }


  GetGeometry()[0].SetLock();

  array_1d<double, 3 >& ContactForce = GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCE);


  for(unsigned int j = 0; j < dimension; j++)
    {
      ContactForce[j] += TangentForceModulus * rVariables.Surface.Tangent[j];
    }

  GetGeometry()[0].UnSetLock();

  VectorType ContactTorque = ZeroVector(3);

  //ContactTorque = MathUtils<double>::CrossProduct( rVariables.CentroidPosition, ContactForceVector);
  ContactTorque = prod(rVariables.SkewSymDistance,ContactForceVector); //  = (D x F)

  // std::cout<<" [ContactTorque]: "<<ContactTorque;
  // std::cout<<" [ContactForce]:  "<<ContactForceVector;
  // std::cout<<" [Normal]:  "<<rVariables.Surface.Normal;
  // std::cout<<std::endl;

  //Contact torque due to contact tangent force on beam surface
  for (unsigned int i =0; i < dimension; ++i) {
    rRightHandSideVector[i+dimension] += ContactTorque[i];
  }

}

//**************************** Calculate Normal Force Modulus ***********************
//***********************************************************************************

double& RigidBodyPointRigidContactCondition::CalculateNormalForceModulus ( double& rNormalForceModulus, ConditionVariables& rVariables )
{

  rNormalForceModulus = (rVariables.Penalty.Normal * rVariables.Gap.Normal);

  return rNormalForceModulus;

}


//**************************** Check Coulomb law for Tangent Contact Force **********
//***********************************************************************************

double RigidBodyPointRigidContactCondition::CalculateCoulombsFrictionLaw(double & rTangentRelativeMovement, double & rNormalForceModulus , ConditionVariables& rVariables)
{
  rVariables.FrictionCoefficient = this->CalculateFrictionCoefficient(rTangentRelativeMovement,rVariables.DeltaTime);


  double TangentForceModulus = rVariables.Penalty.Tangent * rVariables.Gap.Tangent; //+ rVariables.PreviousTangentForceModulus;

  //std::cout<<" Gap.Tangent "<<rVariables.Gap.Tangent<<std::endl;


  if ( fabs(TangentForceModulus) >  rVariables.FrictionCoefficient * fabs(rNormalForceModulus) && fabs(rVariables.Gap.Tangent) > 1e-200) {

    TangentForceModulus = rVariables.FrictionCoefficient * fabs(rNormalForceModulus) ;
    rVariables.Slip = true;

  }
  else {
    rVariables.Slip = false;
  }


  return TangentForceModulus;
}



//**************************** Check friction coefficient ***************************
//***********************************************************************************

double RigidBodyPointRigidContactCondition::CalculateFrictionCoefficient(const double& rTangentRelativeMovement, const double& rDeltaTime)
{
    KRATOS_TRY

    //---FRICTION LAW in function of the relative sliding velocity ---//

    double DynamicFrictionCoefficient = 0.0;//0.2;
    double StaticFrictionCoefficient  = 0.0;//0.3;

    if( GetProperties().Has(FRICTION_ACTIVE) ){
      if( GetProperties()[FRICTION_ACTIVE] ){

	if( GetProperties().Has(MU_DYNAMIC) )
	  DynamicFrictionCoefficient = GetProperties()[MU_DYNAMIC];

	if( GetProperties().Has(MU_STATIC) )
	  StaticFrictionCoefficient = GetProperties()[MU_STATIC];
      }
    }


  double Velocity = 0;
  Velocity = rTangentRelativeMovement / rDeltaTime;


  //Addicional constitutive parameter  C
  //which describes how fast the static coefficient approaches the dynamic:
  double C = 2;//0.1;

  //Addicional constitutive parameter  E
  //regularization parameter (->0, classical Coulomb law)
  double E = 0;//0.01;


  double FrictionCoefficient = DynamicFrictionCoefficient + ( StaticFrictionCoefficient - DynamicFrictionCoefficient ) * exp( (-1) * C * fabs(Velocity) );


  //Square root regularization
  FrictionCoefficient *= fabs(Velocity)/sqrt( ( Velocity * Velocity ) + ( E * E ) );

  //Hyperbolic regularization
  //FrictionCoefficient *= tanh( fabs(Velocity)/E );

  return FrictionCoefficient;

  KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************

void RigidBodyPointRigidContactCondition::VectorToSkewSymmetricTensor( const Vector& rVector,
								       Matrix& rSkewSymmetricTensor )
{
    KRATOS_TRY

    //Initialize Local Matrices
    if( rSkewSymmetricTensor.size1() != 3 )
      rSkewSymmetricTensor.resize(3, 3, false);

    rSkewSymmetricTensor = ZeroMatrix(3,3);

    rSkewSymmetricTensor( 0, 1 ) = -rVector[2];
    rSkewSymmetricTensor( 0, 2 ) =  rVector[1];
    rSkewSymmetricTensor( 1, 2 ) = -rVector[0];

    rSkewSymmetricTensor( 1, 0 ) =  rVector[2];
    rSkewSymmetricTensor( 2, 0 ) = -rVector[1];
    rSkewSymmetricTensor( 2, 1 ) =  rVector[0];


    KRATOS_CATCH( "" )

}


//***********************************************************************************
//***********************************************************************************


int RigidBodyPointRigidContactCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
  return 0;
}

//***********************************************************************************
//***********************************************************************************

} // Namespace Kratos.
