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
#include "custom_conditions/deformable_contact/contact_domain_LM_3D_condition.hpp"

#include "contact_mechanics_application_variables.h"


namespace Kratos
{

  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  ContactDomainLM3DCondition::ContactDomainLM3DCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : ContactDomainCondition( NewId, pGeometry )
  {
    //DO NOT ADD DOFS HERE!!!
  }


  //******************************CONSTRUCTOR*******************************************
  //************************************************************************************

  ContactDomainLM3DCondition::ContactDomainLM3DCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : ContactDomainCondition( NewId, pGeometry, pProperties )
  {
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
  }


  //******************************COPY CONSTRUCTOR**************************************
  //************************************************************************************

  ContactDomainLM3DCondition::ContactDomainLM3DCondition( ContactDomainLM3DCondition const& rOther)
    :ContactDomainCondition(rOther)
  {
  }


  //*******************************ASSIGMENT OPERATOR***********************************
  //************************************************************************************

  ContactDomainLM3DCondition&  ContactDomainLM3DCondition::operator=(ContactDomainLM3DCondition const& rOther)
  {
    ContactDomainCondition::operator=(rOther);

    return *this;
  }


  //*********************************OPERATIONS*****************************************
  //************************************************************************************

  Condition::Pointer ContactDomainLM3DCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
  {
    return Kratos::make_intrusive<ContactDomainLM3DCondition>(NewId, GetGeometry().Create( ThisNodes ), pProperties);
  }

  //************************************CLONE*******************************************
  //************************************************************************************

  Condition::Pointer ContactDomainLM3DCondition::Clone( IndexType NewId, NodesArrayType const& ThisNodes ) const
  {
    return this->Create(NewId, ThisNodes, pGetProperties());
  }

  //*******************************DESTRUCTOR*******************************************
  //************************************************************************************


  ContactDomainLM3DCondition::~ContactDomainLM3DCondition()
  {
  }


  //************* STARTING - ENDING  METHODS
  //************************************************************************************
  //************************************************************************************

  void ContactDomainLM3DCondition::SetMasterGeometry()
  {
    KRATOS_TRY

    Element::ElementType& rMasterElement = GetValue(MASTER_ELEMENTS).back();
    mContactVariables.SetMasterElement(rMasterElement);

    Element::NodeType&    rMasterNode  = GetValue(MASTER_NODES).front();
    mContactVariables.SetMasterNode(rMasterNode);

    int slave = -1;

    Geometry< Node<3> >& rMasterGeometry = rMasterElement.GetGeometry();
    Geometry< Node<3> >& rGeometry       = GetGeometry();


    for(unsigned int i=0; i<rMasterGeometry.PointsNumber(); i++)
      {
        if(rMasterNode.Id()==rMasterGeometry[i].Id())
	  {
	    slave=i;
	  }
      }

    //set the second slave in order
    // unsigned int node_a = 0, node_b = 0;
    // unsigned int slave_a = 0, slave_b = 0;
    // if( this->Is(SELECTED) ){
    // 	unsigned int slave_a = GetValue(MASTER_NODES).front()->Id(), slave_b = GetValue(MASTER_NODES).back()->Id();
    // 	for(unsigned int j=0; j<GetGeometry().PointsNumber(); j++)
    // 	{

    // 	    if(slave_a==rMasterGeometry[j].Id())
    // 	    {
    // 		node_a = j;
    // 		break;
    // 	    }
    // 	}
    // 	for(unsigned int j=0; j<GetGeometry().PointsNumber(); j++)
    // 	{

    // 	    if(slave_b==rMasterGeometry[j].Id())
    // 	    {
    // 		node_b = j;
    // 		break;
    // 	    }
    // 	}

    // 	mContactUtilities.GetOppositeEdge(node_a, node_b, slave_a, slave_b);

    // 	for(unsigned int j=0; j<GetGeometry().PointsNumber(); j++)
    // 	{

    // 	    if(rMasterGeometry[slave_a].Id()==GetGeometry()[j].Id())
    // 	    {
    // 		node_a = j;
    // 		break;
    // 	    }
    // 	}
    // 	for(unsigned int j=0; j<GetGeometry().PointsNumber(); j++)
    // 	{

    // 	    if(rMasterGeometry[slave_b].Id()==GetGeometry()[j].Id())
    // 	    {
    // 		node_b = j;
    // 		break;
    // 	    }
    // 	}
    // }


    //set conectivities a-b-c-d correspondence with 1-2-3-4 (local contact element):
    //nodes 1/2/3 and slave for FaceType (PATCH-A)
    //nodes 1/2 and 3/4 for EdgeType (PATCH-B)

    //order contains the relationship between contact element and body element condident nodes numeration

    if(slave>=0)
      {
	// Clear nodes and slaves before push back quantities
	mContactVariables.nodes.resize(0);
	mContactVariables.slaves.resize(0);

	// Master element positions [0m,1m,2m,3m]
	// Current contact element positions [0c,1c,2c,3c] (i.e. [0c=1m,1c=3m,2c=0m,3c=2m])
	mContactVariables.order.resize(rGeometry.PointsNumber(),false);

	bool iset = false;
        for(unsigned int i=0; i<GetGeometry().PointsNumber(); i++)
	  {
	    iset = false;
	    for(unsigned int j=0; j<GetGeometry().PointsNumber(); j++)
	      {

                if(rGeometry[i].Id()==rMasterGeometry[j].Id())
		  {
                    mContactVariables.order[i] = j;
 		    //std::cout<<" order master ["<<i<<"] = "<<j<<std::endl;
                    iset=true;
		    break;
		  }

	      }

            if(iset==false)
	      {
         	mContactVariables.order[i] = slave;
                mContactVariables.slaves.push_back(i);
	      }

	  }

	//set the second slave in order
	if( this->Is(SELECTED) ){

	  for(unsigned int i=0; i<GetGeometry().PointsNumber(); i++)
	    {
	      iset = false;
	      for(unsigned int j=0; j<mContactVariables.order.size(); j++)
	  	{
	  	  if(mContactVariables.order[j] == i)
	  	    {
	  	      iset = true;
	  	      break;
	  	    }
	  	}

	      if(iset==false)
	  	mContactVariables.order[mContactVariables.slaves.back()] = i;
	    }

	}

	//Permute
	std::vector<unsigned int> permute (7);

	permute[0]=0;
	permute[1]=1;
	permute[2]=2;
	permute[3]=3;
	permute[4]=0;
	permute[5]=1;
	permute[6]=2;


	DenseMatrix<unsigned int> lpofa; //points that define the faces
	GetGeometry().NodesInFaces(lpofa);

	//reorder counter-clock-wise
	if(mContactVariables.slaves.size() == 1){

	  for( unsigned int i=3; i>0; i-- ) //counterclockwise order to force inside normal
	    {
	      mContactVariables.nodes.push_back(lpofa(i,mContactVariables.slaves.front()));
	    }

	  mContactVariables.nodes.push_back(mContactVariables.slaves.front());

	}
	else if(mContactVariables.slaves.size() == 2){

	    // mContactUtilities.GetOppositeEdge(node_a, node_b, slave_a, slave_b);
	    // mContactVariables.nodes.push_back(node_a);
	    // mContactVariables.nodes.push_back(node_b);
	    // mContactVariables.slaves.front() = slave_a;
	    // mContactVariables.slaves.back()  = slave_b;

	    if( mContactVariables.slaves.back() < mContactVariables.slaves.front() )
		std::iter_swap(mContactVariables.slaves.begin(), mContactVariables.slaves.begin()+1);

	    //now slave.front() < slave.back(): possible order
	    if( mContactVariables.slaves.back() == permute[mContactVariables.slaves.front()+1] ){

		if( mContactVariables.slaves.front() != 1 ){
		    //slaves[3,4]: (0,1) (2,3) | masters[1,2]: (2,3) (0,1)
		    mContactVariables.nodes.push_back(permute[mContactVariables.slaves.back()+1]);
		    mContactVariables.nodes.push_back(permute[mContactVariables.slaves.back()+2]);
		}
		else{
		    //slaves[3,4]: (1,2) | masters[1,2]: (0,3)
		    mContactVariables.nodes.push_back(permute[mContactVariables.slaves.back()+2]);
		    mContactVariables.nodes.push_back(permute[mContactVariables.slaves.back()+1]);
		}
	    }
	    else if( mContactVariables.slaves.back() == permute[mContactVariables.slaves.front()+2] ){

		if( mContactVariables.slaves.front() == 0 ){
		    //slaves[3,4]: (0,2) | masters[1,2]: (3,1)
		    mContactVariables.nodes.push_back(permute[mContactVariables.slaves.back()+1]);
		    mContactVariables.nodes.push_back(permute[mContactVariables.slaves.front()+1]);
		}
		else{
		    //slaves[3,4]: (1,3) | masters[1,2]: (2,0)
		    mContactVariables.nodes.push_back(permute[mContactVariables.slaves.front()+1]);
		    mContactVariables.nodes.push_back(permute[mContactVariables.slaves.back()+1]);
		}

	    }
	    else if( mContactVariables.slaves.back() == permute[mContactVariables.slaves.front()+3] ){

		//slaves[3,4]: (0,3) | masters[1,2]: (1,2)
		mContactVariables.nodes.push_back(permute[mContactVariables.slaves.front()+1]);
		mContactVariables.nodes.push_back(permute[mContactVariables.slaves.back()-1]);
	    }
	    else{
		std::cout<<" Something is wrong with EDGE slaves distribution "<<std::endl;
	    }

	    // if( mContactVariables.slaves.back() < mContactVariables.slaves.front() ){
	    // 	std::iter_swap(mContactVariables.slaves.begin(), mContactVariables.slaves.begin()+1);
	    // }

	    // //now slave.front() < slave.back(): possible order
	    // if( mContactVariables.slaves.back() == permute[mContactVariables.slaves.front()+1] ){

	    // 	//(0,1) (1,2) (2,3) | (2,3) (3,0) (0,1)
	    // 	mContactVariables.nodes.push_back(permute[mContactVariables.slaves.back()+1]);
	    // 	mContactVariables.nodes.push_back(permute[mContactVariables.slaves.back()+2]);
	    // }
	    // else if( mContactVariables.slaves.back() == permute[mContactVariables.slaves.front()+2] ){

	    // 	//(1,3) (0,2) | (0,2) (1,3)
	    // 	mContactVariables.nodes.push_back(permute[mContactVariables.slaves.back()+1]);
	    // 	mContactVariables.nodes.push_back(permute[mContactVariables.slaves.back()-1]);
	    // }
	    // else if( mContactVariables.slaves.back() == permute[mContactVariables.slaves.front()+3] ){

	    // 	//(0,3) | (2,1)
	    // 	mContactVariables.nodes.push_back(permute[mContactVariables.slaves.back()-1]);
	    // 	mContactVariables.nodes.push_back(permute[mContactVariables.slaves.front()+1]);
	    // }
	    // else{
	    // 	std::cout<<" Something is wrong with EDGE slaves distribution "<<std::endl;
	    // }

	    mContactVariables.nodes.push_back(mContactVariables.slaves.front());
	    mContactVariables.nodes.push_back(mContactVariables.slaves.back());

	    //check correct slave edge orientation:

	    //Reference Position
	    PointType V1 = GetGeometry()[mContactVariables.nodes[1]].Coordinates()-GetGeometry()[mContactVariables.nodes[0]].Coordinates();
	    PointType V2 = GetGeometry()[mContactVariables.nodes[2]].Coordinates()-GetGeometry()[mContactVariables.nodes[3]].Coordinates();

	    //Compute Reference Normal
	    PointType Normal = mContactUtilities.CalculateSurfaceNormal(Normal,V1,V2);

	    V2 = GetGeometry()[mContactVariables.nodes[0]].FastGetSolutionStepValue(NORMAL);

	    if((inner_prod(V2,Normal))<0){ //change the slaves nodes order for consistency of element orientation
		std::cout<<" normal away "<<std::endl;
		std::iter_swap(mContactVariables.slaves.begin(), mContactVariables.slaves.begin()+1);
		std::iter_swap(mContactVariables.nodes.end()-2, mContactVariables.nodes.end()-1);
	    }


	}
	else{
	  std::cout<<" Problem in the SLAVES DETECTION "<<std::endl;
	}

	mContactVariables.SetMasterGeometry( rMasterElement.GetGeometry() );

      }
    else{
      KRATOS_THROW_ERROR( std::invalid_argument, "MASTERNODE do not belongs to MASTER ELEMENT", "" )
	}

    if( this->Is(SELECTED) )
	std::cout<<" CONTACT EDGE("<<GetGeometry()[0].Id()<<")("<<GetGeometry()[1].Id()<<")("<<GetGeometry()[2].Id()<<")("<<GetGeometry()[3].Id()<<")"<<std::endl;
    else
	std::cout<<" CONTACT ("<<GetGeometry()[0].Id()<<")("<<GetGeometry()[1].Id()<<")("<<GetGeometry()[2].Id()<<")("<<GetGeometry()[3].Id()<<")"<<std::endl;

    std::cout<<" MASTER  ("<<rMasterElement.GetGeometry()[0].Id()<<")("<<rMasterElement.GetGeometry()[1].Id()<<")("<<rMasterElement.GetGeometry()[2].Id()<<")("<<rMasterElement.GetGeometry()[3].Id()<<")"<<std::endl;
    std::cout<<" nodes:["<<mContactVariables.nodes[0]<<" "<<mContactVariables.nodes[1]<<" "<<mContactVariables.nodes[2]<<" "<<mContactVariables.nodes[3]<<"]"<<" slaves:["<<mContactVariables.slaves[0]<<" "<<mContactVariables.slaves[1]<<"]"<<std::endl;

    std::cout<<" order:["<<mContactVariables.order[0]<<" "<<mContactVariables.order[1]<<" "<<mContactVariables.order[2]<<" "<<mContactVariables.order[3]<<"]"<<std::endl;
    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************



  void ContactDomainLM3DCondition::CalculatePreviousGap() //prediction of the lagrange multiplier
  {
    //std::cout<<" Element Type EDGE "<< this->Is(SELECTED) <<std::endl;

    if( this->Is(SELECTED) )
      CalculatePreviousGapEdgeType();
    else
      CalculatePreviousGapFaceType();

  }


  //************************************************************************************
  //************************************************************************************


  void ContactDomainLM3DCondition::TransformCovariantToContravariantBase(SurfaceBase& Covariant,SurfaceBase& Contravariant) //prediction of the lagrange multiplier face type element
  {

    // Covariant metric
    Covariant.Metric.resize(3,3,false);
    array_1d<double, 3> DirectionC;
    MathUtils<double>::CrossProduct(DirectionC,Covariant.DirectionA,Covariant.DirectionB);
    if( norm_2(DirectionC)!=0 )
	DirectionC /= norm_2(DirectionC);

    Covariant.Metric(0,0) = inner_prod(Covariant.DirectionA, Covariant.DirectionA);
    Covariant.Metric(0,1) = inner_prod(Covariant.DirectionA, Covariant.DirectionB);
    Covariant.Metric(0,2) = inner_prod(Covariant.DirectionA, DirectionC);

    Covariant.Metric(1,0) = inner_prod(Covariant.DirectionB, Covariant.DirectionA);
    Covariant.Metric(1,1) = inner_prod(Covariant.DirectionB, Covariant.DirectionB);
    Covariant.Metric(1,2) = inner_prod(Covariant.DirectionB, DirectionC);

    Covariant.Metric(2,0) = inner_prod(DirectionC, Covariant.DirectionA);
    Covariant.Metric(2,1) = inner_prod(DirectionC, Covariant.DirectionB);
    Covariant.Metric(2,2) = inner_prod(DirectionC, DirectionC);

    // Contravariant vectors and contravariant metric
    Contravariant.Metric.resize(3,3,false);
    double MetricDet;
    MathUtils<double>::InvertMatrix3(Covariant.Metric,Contravariant.Metric, MetricDet);

    //transform DirectionA to the contravariant base
    Contravariant.DirectionA = Covariant.DirectionA * Contravariant.Metric(0,0) + Covariant.DirectionB * Contravariant.Metric(1,0);

    //transform DirectionB to the contravariant base
    Contravariant.DirectionB = Covariant.DirectionA * Contravariant.Metric(0,1) + Covariant.DirectionB * Contravariant.Metric(1,1);

    // check:
    // array_1d<double, 3> DirectionC;
    // MathUtils<double>::CrossProduct(DirectionC,Covariant.DirectionA,Covariant.DirectionB);
    // double V = norm_2(DirectionC);
    // Contravariant.Metric(0,0) = inner_prod(Contravariant.DirectionA, Contravariant.DirectionA);
    // Contravariant.Metric(0,1) = inner_prod(Contravariant.DirectionA, Contravariant.DirectionB);
    // Contravariant.Metric(1,0) = inner_prod(Contravariant.DirectionB, Contravariant.DirectionA);
    // Contravariant.Metric(1,1) = inner_prod(Contravariant.DirectionB, Contravariant.DirectionB);

    // std::cout<<" Check "<<inner_prod(Contravariant.DirectionA, Covariant.DirectionB)<<" "<<inner_prod(Contravariant.DirectionB, Covariant.DirectionA)<<" Det "<<1.0/MathUtils<double>::Det(Contravariant.Metric)<<" "<<MathUtils<double>::Det(Covariant.Metric)<<std::endl;

    // std::cout<<" CvMetric "<<Covariant.Metric<<" CvA "<<Covariant.DirectionA<<" CvB "<<Covariant.DirectionB<<std::endl;
    // std::cout<<" CnMetric "<<Contravariant.Metric<<" CnA "<<Contravariant.DirectionA<<" CnB "<<Contravariant.DirectionB<<std::endl;


  }

  //************************************************************************************
  //************************************************************************************


  void ContactDomainLM3DCondition::CalculatePreviousGapFaceType() //prediction of the lagrange multiplier face type element
  {

    //Contact face node1-node2-node3
    unsigned int node1=mContactVariables.nodes[0];
    unsigned int node2=mContactVariables.nodes[1];
    unsigned int node3=mContactVariables.nodes[2];
    unsigned int slave=mContactVariables.slaves.front();

    //Get Reference Normal (alternative)
    //mContactVariables.ReferenceSurface.Normal = GetValue(NORMAL);

    PointType P1  =  GetGeometry()[node1].Coordinates() - (GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,0) - GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1) );
    PointType P2  =  GetGeometry()[node2].Coordinates() - (GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,0) - GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1) );
    PointType P3  =  GetGeometry()[node3].Coordinates() - (GetGeometry()[node3].FastGetSolutionStepValue(DISPLACEMENT,0) - GetGeometry()[node3].FastGetSolutionStepValue(DISPLACEMENT,1) );
    PointType PS  =  GetGeometry()[slave].Coordinates() - (GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,0) - GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,1) );

    //Set Reference Tangent
    mContactVariables.Tangent.CovariantBase.DirectionA = P2 - P1;
    mContactVariables.Tangent.CovariantBase.DirectionB = P3 - P1;

    TransformCovariantToContravariantBase(mContactVariables.Tangent.CovariantBase, mContactVariables.Tangent.ContravariantBase);

    //Compute Current Normal
    mContactVariables.ReferenceSurface.Normal = mContactUtilities.CalculateSurfaceNormal(mContactVariables.ReferenceSurface.Normal, mContactVariables.Tangent.CovariantBase.DirectionA, mContactVariables.Tangent.CovariantBase.DirectionB);

    PointType V = GetGeometry()[node1].FastGetSolutionStepValue(NORMAL);

    if((inner_prod(V,mContactVariables.ReferenceSurface.Normal))<0) //change the slaves nodes order for consistency
      std::cout<<" IT IS A INSIDE NORMAL !!! "<<std::endl;

    //1.- Compute tension vector:  (must be updated each iteration)
    Matrix StressMatrix(3,3);
    noalias(StressMatrix) = ZeroMatrix(3,3);
    Matrix F(3,3);
    noalias(F) = ZeroMatrix(3,3);

    //a.- Assign initial 2nd Piola Kirchhoff stress:
    Condition* MasterCondition = GetValue(MASTER_CONDITION).get();

    //Get previous mechanics stored in the master node/condition
    Vector StressVector;
    StressVector = MasterCondition->GetValue(CAUCHY_STRESS_VECTOR);  //it means that has been stored
    F            = MasterCondition->GetValue(DEFORMATION_GRADIENT);  //it means that has been stored


    StressMatrix = MathUtils<double>::StressVectorToTensor( StressVector );

    //we are going to need F here from Cn-1 to Cn
    // F0 from C0 to Cn is need for the stress recovery on domain elements

    double detF =MathUtils<double>::Det(F);

    //b.- Compute the 1srt Piola Kirchhoff stress tensor  (P=J*CauchyStress*F^-T)
    ConstitutiveLaw Constitutive;
    Constitutive.TransformStresses(StressMatrix,F,detF,ConstitutiveLaw::StressMeasure_Cauchy,ConstitutiveLaw::StressMeasure_PK1);

    //c.- Compute (n-1) normals, tangents and relative displacements from historic mX on boundaries

    //Previous normal and tangent:  n_n-1,t_n-1

    //Previous Position
    P1  =  GetGeometry()[node1].Coordinates() - (GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1) - GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,2) );
    P2  =  GetGeometry()[node2].Coordinates() - (GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1) - GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,2) );
    P3  =  GetGeometry()[node3].Coordinates() - (GetGeometry()[node3].FastGetSolutionStepValue(DISPLACEMENT,1) - GetGeometry()[node3].FastGetSolutionStepValue(DISPLACEMENT,2) );
    PS  =  GetGeometry()[slave].Coordinates() - (GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,1) - GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,2) );

    //Set Previous Tangent
    PointType V1 = P2 - P1;
    PointType V2 = P3 - P1;

    //Compute Previous Normal
    mContactVariables.PreStepSurface.Normal=mContactUtilities.CalculateSurfaceNormal(mContactVariables.PreStepSurface.Normal, V1, V2);


    if((inner_prod(mContactVariables.PreStepSurface.Normal,mContactVariables.ReferenceSurface.Normal))<0){ //to give the correct direction
      std::cout<<" Previous Contact Normal Sense INVERTED Face "<<mContactVariables.PreStepSurface.Normal<<std::endl;
      mContactVariables.PreStepSurface.Normal*=-1;
    }

    if(!(norm_2(mContactVariables.PreStepSurface.Normal)))
      {
        mContactVariables.PreStepSurface.Normal=mContactVariables.ReferenceSurface.Normal;
      }


    //d.- Compute the tension (or traction) vector T=P*N (in the Reference configuration)
    mContactVariables.TractionVector=prod(StressMatrix,mContactVariables.PreStepSurface.Normal);

    //e.- Compute A_n-1,B_n-1,L_n-1

    //A_n-1, B_n-1, L_n-1:
    std::vector<BaseLengths> PreviousBase(3);
    mContactUtilities.CalculateBaseDistances(PreviousBase,P1,P2,P3,PS,mContactVariables.PreStepSurface.Normal);
    double EquivalentArea = 1;
    mContactUtilities.CalculateBaseArea(EquivalentArea,PreviousBase[0].L,PreviousBase[1].L,PreviousBase[2].L);
    double FactorArea = sqrt(EquivalentArea);

    //complete the computation of the stabilization gap
    double ContactFactor = mContactVariables.StabilizationFactor * FactorArea;
    double ContactFactorTangent = ContactFactor * GetProperties()[TANGENTIAL_PENALTY_RATIO];

    //f.-obtain the (g_N)3 and (g_T)3 for the n-1 configuration

    mContactVariables.PreviousGap.Normal = 0;
    mContactVariables.PreviousGap.Normal = inner_prod((PS-P1),mContactVariables.PreStepSurface.Normal);


    mContactVariables.Tangent.PreviousGapA.Covariant = mContactVariables.PreviousGap.Normal * inner_prod(mContactVariables.Tangent.CovariantBase.DirectionA,mContactVariables.PreStepSurface.Normal);
    mContactVariables.Tangent.PreviousGapB.Covariant = mContactVariables.PreviousGap.Normal * inner_prod(mContactVariables.Tangent.CovariantBase.DirectionB,mContactVariables.PreStepSurface.Normal);

    mContactVariables.Tangent.PreviousGapA.Contravariant = mContactVariables.PreviousGap.Normal * inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionA,mContactVariables.PreStepSurface.Normal);
    mContactVariables.Tangent.PreviousGapB.Contravariant = mContactVariables.PreviousGap.Normal * inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionB,mContactVariables.PreStepSurface.Normal);


    PointType D1  =  GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1)-GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,2);
    PointType D2  =  GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1)-GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,2);
    PointType D3  =  GetGeometry()[node3].FastGetSolutionStepValue(DISPLACEMENT,1)-GetGeometry()[node3].FastGetSolutionStepValue(DISPLACEMENT,2);
    PointType DS  =  GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,1)-GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,2);


    //(g_N)3
    mContactVariables.PreviousGap.Normal*= inner_prod(mContactVariables.ReferenceSurface.Normal,mContactVariables.PreStepSurface.Normal);
    mContactVariables.PreviousGap.Normal+=inner_prod(mContactVariables.ReferenceSurface.Normal,(D1*(-PreviousBase[0].A/PreviousBase[0].L)));
    mContactVariables.PreviousGap.Normal+=inner_prod(mContactVariables.ReferenceSurface.Normal,(D2*(-PreviousBase[1].A/PreviousBase[1].L)));
    mContactVariables.PreviousGap.Normal+=inner_prod(mContactVariables.ReferenceSurface.Normal,(D3*(-PreviousBase[2].A/PreviousBase[2].L)));
    mContactVariables.PreviousGap.Normal+=inner_prod(mContactVariables.ReferenceSurface.Normal,DS);


    //(g_T)3
    mContactVariables.Tangent.PreviousGapA.Covariant+=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionA,(D1*(-PreviousBase[0].A/PreviousBase[0].L)));
    mContactVariables.Tangent.PreviousGapA.Covariant+=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionA,(D2*(-PreviousBase[1].A/PreviousBase[1].L)));
    mContactVariables.Tangent.PreviousGapA.Covariant+=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionA,(D3*(-PreviousBase[2].A/PreviousBase[2].L)));
    mContactVariables.Tangent.PreviousGapA.Covariant+=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionA,DS);

    mContactVariables.Tangent.PreviousGapB.Covariant+=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionB,(D1*(-PreviousBase[0].A/PreviousBase[0].L)));
    mContactVariables.Tangent.PreviousGapB.Covariant+=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionB,(D2*(-PreviousBase[1].A/PreviousBase[1].L)));
    mContactVariables.Tangent.PreviousGapB.Covariant+=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionB,(D3*(-PreviousBase[2].A/PreviousBase[2].L)));
    mContactVariables.Tangent.PreviousGapB.Covariant+=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionB,DS);


    mContactVariables.Tangent.PreviousGapA.Contravariant+=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionA,(D1*(-PreviousBase[0].A/PreviousBase[0].L)));
    mContactVariables.Tangent.PreviousGapA.Contravariant+=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionA,(D2*(-PreviousBase[1].A/PreviousBase[1].L)));
    mContactVariables.Tangent.PreviousGapA.Contravariant+=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionA,(D3*(-PreviousBase[2].A/PreviousBase[2].L)));
    mContactVariables.Tangent.PreviousGapA.Contravariant+=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionA,DS);


    mContactVariables.Tangent.PreviousGapB.Contravariant+=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionB,(D1*(-PreviousBase[0].A/PreviousBase[0].L)));
    mContactVariables.Tangent.PreviousGapB.Contravariant+=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionB,(D2*(-PreviousBase[1].A/PreviousBase[1].L)));
    mContactVariables.Tangent.PreviousGapB.Contravariant+=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionB,(D3*(-PreviousBase[2].A/PreviousBase[2].L)));
    mContactVariables.Tangent.PreviousGapB.Contravariant+=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionB,DS);


    //d_n-1=X_n - X_n-1

    //g.- get total effective gap as: gap_n^eff=gap_n+(PreviousTimeStep/CurrentTimeStep)*(gap_n-gap_n-1)

    //gap_n-1 (in function of the n-1 position of hte other node) gap_n-1=(g_N)3_n-1+2*Tau*tn_n-1

    double NormalTensil=0,TangentTensilcvA=0,TangentTensilcvB=0,TangentTensilcnA=0,TangentTensilcnB=0;

    //h.- Compute normal component of the tension vector:   (tn=n·P·N)
    NormalTensil=inner_prod(mContactVariables.ReferenceSurface.Normal,mContactVariables.TractionVector);


    //i.- Compute tangent component of the tension vector:  (tt=cvt·P·N)
    TangentTensilcvA=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionA,mContactVariables.TractionVector);
    TangentTensilcvB=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionB,mContactVariables.TractionVector);


    //j.- Compute tangent component of the tension vector:  (tt=cnt·P·N)
    TangentTensilcnA=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionA,mContactVariables.TractionVector);
    TangentTensilcnB=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionB,mContactVariables.TractionVector);


    mContactVariables.PreviousGap.Normal  += 3 * ContactFactor * NormalTensil;

    mContactVariables.Tangent.PreviousGapA.Covariant += 3 * ContactFactorTangent * TangentTensilcvA;
    mContactVariables.Tangent.PreviousGapB.Covariant += 3 * ContactFactorTangent * TangentTensilcvB;

    mContactVariables.Tangent.PreviousGapA.Contravariant += 3 * ContactFactorTangent * TangentTensilcnA;
    mContactVariables.Tangent.PreviousGapB.Contravariant += 3 * ContactFactorTangent * TangentTensilcnB;

    // std::cout<<"ConditionID:  "<<this->Id()<<std::endl;
    // std::cout<<" -> Previous Tractions [tN:"<<NormalTensil<<"] "<<std::endl;
    // std::cout<<" Previous Normal Gap [gN:"<<mContactVariables.PreviousGap.Normal<<"] "<<std::endl;
  }


  //************************************************************************************
  //************************************************************************************

  void ContactDomainLM3DCondition::CalculatePreviousGapEdgeType() //prediction of the lagrange multiplier face type element
  {

    //Contact face segment node1-node2
    unsigned int node1=mContactVariables.nodes[0];
    unsigned int node2=mContactVariables.nodes[1];
    unsigned int slave1=mContactVariables.slaves[0];
    unsigned int slave2=mContactVariables.slaves[1];


    //Reference Position
    PointType P1  =  GetGeometry()[node1].Coordinates() - (GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,0) - GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1) );
    PointType P2  =  GetGeometry()[node2].Coordinates() - (GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,0) - GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1) );
    PointType PS1  =  GetGeometry()[slave1].Coordinates() - (GetGeometry()[slave1].FastGetSolutionStepValue(DISPLACEMENT,0) - GetGeometry()[slave1].FastGetSolutionStepValue(DISPLACEMENT,1) );
    PointType PS2  =  GetGeometry()[slave2].Coordinates() - (GetGeometry()[slave2].FastGetSolutionStepValue(DISPLACEMENT,0) - GetGeometry()[slave2].FastGetSolutionStepValue(DISPLACEMENT,1) );


    //Set Reference Tangent
    mContactVariables.Tangent.CovariantBase.DirectionA = P2 - P1;
    mContactVariables.Tangent.CovariantBase.DirectionB = PS1 - PS2;

    TransformCovariantToContravariantBase(mContactVariables.Tangent.CovariantBase, mContactVariables.Tangent.ContravariantBase);

    //Compute Reference Normal
    mContactVariables.ReferenceSurface.Normal=mContactUtilities.CalculateSurfaceNormal(mContactVariables.ReferenceSurface.Normal,mContactVariables.Tangent.CovariantBase.DirectionA,mContactVariables.Tangent.CovariantBase.DirectionB);


    //1.- Compute tension vector:  (must be updated each iteration)
    Matrix StressMatrix(3,3);
    noalias(StressMatrix) = ZeroMatrix(3,3);
    Matrix F(3,3);
    noalias(F) = ZeroMatrix(3,3);

    //a.- Assign initial 2nd Piola Kirchhoff stress:
    Condition* MasterCondition = GetValue(MASTER_CONDITION).get();

    //Get previous mechanics stored in the master node/condition
    Vector StressVector;
    StressVector = MasterCondition->GetValue(CAUCHY_STRESS_VECTOR);  //it means that has been stored
    F            = MasterCondition->GetValue(DEFORMATION_GRADIENT);  //it means that has been stored

    StressMatrix = MathUtils<double>::StressVectorToTensor( StressVector );

    //we are going to need F here from Cn-1 to Cn
    // F0 from C0 to Cn is need for the stress recovery on domain elements

    double detF =MathUtils<double>::Det(F);

    //b.- Compute the 1srt Piola Kirchhoff stress tensor  (P=J*CauchyStress*F^-T)
    ConstitutiveLaw Constitutive;
    Constitutive.TransformStresses(StressMatrix,F,detF,ConstitutiveLaw::StressMeasure_Cauchy,ConstitutiveLaw::StressMeasure_PK1);


    //c.- Compute (n-1) normals, tangents and relative displacements from historic mX on boundaries

    //Previous normal and tangent:  n_n-1,t_n-1

    //Previous Position
    P1  =  GetGeometry()[node1].Coordinates() - (GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1) - GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,2) );
    P2  =  GetGeometry()[node2].Coordinates() - (GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1) - GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,2) );
    PS1  =  GetGeometry()[slave1].Coordinates() - (GetGeometry()[slave1].FastGetSolutionStepValue(DISPLACEMENT,1) - GetGeometry()[slave1].FastGetSolutionStepValue(DISPLACEMENT,2) );
    PS2  =  GetGeometry()[slave2].Coordinates() - (GetGeometry()[slave2].FastGetSolutionStepValue(DISPLACEMENT,1) - GetGeometry()[slave2].FastGetSolutionStepValue(DISPLACEMENT,2) );

    //Set Previous Tangent
    PointType V1 = P2-P1;
    PointType V2 = PS1-PS2;

    //Compute Previous Normal
    mContactVariables.PreStepSurface.Normal=mContactUtilities.CalculateSurfaceNormal(mContactVariables.PreStepSurface.Normal,V1,V2);

    if((inner_prod(mContactVariables.PreStepSurface.Normal,mContactVariables.ReferenceSurface.Normal))<0){ //to give the correct direction
      std::cout<<" Previous Contact Normal Sense INVERTED Edge "<<mContactVariables.PreStepSurface.Normal<<std::endl;
      mContactVariables.PreStepSurface.Normal*=-1;
    }

    if(!(norm_2(mContactVariables.PreStepSurface.Normal)))
      {
        mContactVariables.PreStepSurface.Normal=mContactVariables.ReferenceSurface.Normal;
      }


    //d.- Compute the tension (or traction) vector T=P*N (in the Reference configuration)
    mContactVariables.TractionVector=prod(StressMatrix,mContactVariables.PreStepSurface.Normal);


    //Reference normal: n_n,t_n  -> mContactVariables.ReferenceSurface.Normal / mContactVariables.ReferenceSurface.Tangent

    //e.- Compute A_n-1,B_n-1,L_n-1

    // std::cout<<" Pre Normal ["<<this->Id()<<"] "<<mContactVariables.PreStepSurface.Normal<<std::endl;

    // std::cout<<" P1 ("<<P1[0]<<" "<<P1[1]<<" "<<P1[2]<<")"<<std::endl;
    // std::cout<<" P2 ("<<P2[0]<<" "<<P2[1]<<" "<<P2[2]<<")"<<std::endl;

    // std::cout<<" PS1 ("<<PS1[0]<<" "<<PS1[1]<<" "<<PS1[2]<<")"<<std::endl;
    // std::cout<<" PS2 ("<<PS2[0]<<" "<<PS2[1]<<" "<<PS2[2]<<")"<<std::endl;


    //A_n-1, B_n-1, L_n-1:
    std::vector<BaseLengths> PreviousBase(3);
    mContactUtilities.CalculateEdgeDistances(PreviousBase,P1,P2,PS1,PS2,mContactVariables.PreStepSurface.Normal);

    PointType NormalDirection;
    MathUtils<double>::CrossProduct(NormalDirection, mContactVariables.Tangent.CovariantBase.DirectionA,mContactVariables.Tangent.CovariantBase.DirectionB);
    double EquivalentArea = 0.5 * norm_2( NormalDirection );
    double FactorArea = 0.25 * (PreviousBase[0].L + PreviousBase[1].L) * (PreviousBase[0].L + PreviousBase[1].L);

    //complete the computation of the stabilization gap
    double ContactFactor = mContactVariables.StabilizationFactor * FactorArea;
    double ContactFactorTangent = ContactFactor * GetProperties()[TANGENTIAL_PENALTY_RATIO];

    //f.-obtain the (g_N)3 and (g_T)3 for the n-1 configuration

    mContactVariables.PreviousGap.Normal  = 0;
    V1 = 0.5 * (P1+P2);
    V2 = 0.5 * (PS1+PS2);
    mContactVariables.PreviousGap.Normal = inner_prod((V2-V1),mContactVariables.PreStepSurface.Normal);


    mContactVariables.Tangent.PreviousGapA.Covariant = mContactVariables.PreviousGap.Normal * inner_prod(mContactVariables.Tangent.CovariantBase.DirectionA,mContactVariables.PreStepSurface.Normal);
    mContactVariables.Tangent.PreviousGapB.Covariant = mContactVariables.PreviousGap.Normal * inner_prod(mContactVariables.Tangent.CovariantBase.DirectionB,mContactVariables.PreStepSurface.Normal);

    mContactVariables.Tangent.PreviousGapA.Contravariant = mContactVariables.PreviousGap.Normal * inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionA,mContactVariables.PreStepSurface.Normal);
    mContactVariables.Tangent.PreviousGapB.Contravariant = mContactVariables.PreviousGap.Normal * inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionB,mContactVariables.PreStepSurface.Normal);

    PointType D1  =  GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1)-GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,2);
    PointType D2  =  GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1)-GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,2);
    PointType DS1  =  GetGeometry()[slave1].FastGetSolutionStepValue(DISPLACEMENT,1)-GetGeometry()[slave1].FastGetSolutionStepValue(DISPLACEMENT,2);
    PointType DS2  =  GetGeometry()[slave2].FastGetSolutionStepValue(DISPLACEMENT,1)-GetGeometry()[slave2].FastGetSolutionStepValue(DISPLACEMENT,2);

    //(g_N)3
    mContactVariables.PreviousGap.Normal*= inner_prod(mContactVariables.ReferenceSurface.Normal,mContactVariables.PreStepSurface.Normal);

    mContactVariables.PreviousGap.Normal+=inner_prod(mContactVariables.ReferenceSurface.Normal,(D1*(-PreviousBase[0].A/PreviousBase[0].L)));
    mContactVariables.PreviousGap.Normal+=inner_prod(mContactVariables.ReferenceSurface.Normal,(D2*(-PreviousBase[0].B/PreviousBase[0].L)));
    mContactVariables.PreviousGap.Normal+=inner_prod(mContactVariables.ReferenceSurface.Normal,(DS1*(PreviousBase[1].A/PreviousBase[1].L)));
    mContactVariables.PreviousGap.Normal+=inner_prod(mContactVariables.ReferenceSurface.Normal,(DS1*(PreviousBase[1].B/PreviousBase[1].L)));

    //(g_T)3
    mContactVariables.Tangent.PreviousGapA.Covariant+=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionA,(D1*(-PreviousBase[0].A/PreviousBase[0].L)));
    mContactVariables.Tangent.PreviousGapA.Covariant+=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionA,(D2*(-PreviousBase[0].B/PreviousBase[0].L)));
    mContactVariables.Tangent.PreviousGapA.Covariant+=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionA,(DS1*(PreviousBase[1].A/PreviousBase[1].L)));
    mContactVariables.Tangent.PreviousGapA.Covariant+=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionA,(DS2*(PreviousBase[1].B/PreviousBase[1].L)));

    mContactVariables.Tangent.PreviousGapB.Covariant+=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionB,(D1*(-PreviousBase[0].A/PreviousBase[0].L)));
    mContactVariables.Tangent.PreviousGapB.Covariant+=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionB,(D2*(-PreviousBase[0].B/PreviousBase[0].L)));
    mContactVariables.Tangent.PreviousGapB.Covariant+=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionB,(DS1*(PreviousBase[1].A/PreviousBase[1].L)));
    mContactVariables.Tangent.PreviousGapB.Covariant+=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionB,(DS2*(PreviousBase[1].B/PreviousBase[1].L)));


    mContactVariables.Tangent.PreviousGapA.Contravariant+=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionA,(D1*(-PreviousBase[0].A/PreviousBase[0].L)));
    mContactVariables.Tangent.PreviousGapA.Contravariant+=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionA,(D2*(-PreviousBase[0].B/PreviousBase[0].L)));
    mContactVariables.Tangent.PreviousGapA.Contravariant+=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionA,(DS1*(PreviousBase[1].A/PreviousBase[1].L)));
    mContactVariables.Tangent.PreviousGapA.Contravariant+=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionA,(DS2*(PreviousBase[1].B/PreviousBase[1].L)));

    mContactVariables.Tangent.PreviousGapB.Contravariant+=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionB,(D1*(-PreviousBase[0].A/PreviousBase[0].L)));
    mContactVariables.Tangent.PreviousGapB.Contravariant+=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionB,(D2*(-PreviousBase[0].B/PreviousBase[0].L)));
    mContactVariables.Tangent.PreviousGapB.Contravariant+=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionB,(DS1*(PreviousBase[1].A/PreviousBase[1].L)));
    mContactVariables.Tangent.PreviousGapB.Contravariant+=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionB,(DS2*(PreviousBase[1].B/PreviousBase[1].L)));

    double EquivalentHeigh = EquivalentArea /PreviousBase[0].L;


    //d_n-1=X_n - X_n-1

    //g.- get total effective gap as: gap_n^eff=gap_n+(PreviousTimeStep/CurrentTimeStep)*(gap_n-gap_n-1)

    //gap_n-1 (in function of the n-1 position of hte other node) gap_n-1=(g_N)3_n-1+2*Tau*tn_n-1

    double NormalTensil=0,TangentTensilcvA=0,TangentTensilcvB=0,TangentTensilcnA=0,TangentTensilcnB=0;

    //h.- Compute normal component of the tension vector:   (tn=n·P·N)
    NormalTensil=inner_prod(mContactVariables.ReferenceSurface.Normal,mContactVariables.TractionVector);


    //i.- Compute tangent component of the tension vector:  (tt=cvt·P·N)
    TangentTensilcvA=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionA,mContactVariables.TractionVector);
    TangentTensilcvB=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionB,mContactVariables.TractionVector);


    //j.- Compute tangent component of the tension vector:  (tt=cnt·P·N)
    TangentTensilcnA=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionA,mContactVariables.TractionVector);
    TangentTensilcnB=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionB,mContactVariables.TractionVector);


    mContactVariables.PreviousGap.Normal  += 3 * ContactFactor * NormalTensil / EquivalentHeigh;

    mContactVariables.Tangent.PreviousGapA.Covariant += 3 * ContactFactorTangent * TangentTensilcvA / EquivalentHeigh;
    mContactVariables.Tangent.PreviousGapB.Covariant += 3 * ContactFactorTangent * TangentTensilcvB / EquivalentHeigh;

    mContactVariables.Tangent.PreviousGapA.Contravariant += 3 * ContactFactorTangent * TangentTensilcnA / EquivalentHeigh;
    mContactVariables.Tangent.PreviousGapB.Contravariant += 3 * ContactFactorTangent * TangentTensilcnB / EquivalentHeigh;


    // std::cout<<"ConditionID:  "<<this->Id()<<std::endl;
    // std::cout<<" -> Previous Tractions [tN:"<<NormalTensil<<"] "<<std::endl;
    // std::cout<<" Previous Normal Gap [gN:"<<mContactVariables.PreviousGap.Normal<<"] "<<std::endl;
  }


  //**********************************COMPUTE TAU STAB**********************************
  //************************************************************************************


  void ContactDomainLM3DCondition::CalculateContactFactor( ProcessInfo& rCurrentProcessInfo )
  {
    //Initilialize Tau for the stabilization
    double alpha_stab = 0.1;
    alpha_stab = GetProperties()[TAU_STAB];

    ElementType& rMasterElement = mContactVariables.GetMasterElement();

    //Look at the nodes, get the slave and get the Emin

    //Contact face segment node1-node2
    unsigned int slave = mContactVariables.slaves.front();

    const Properties& SlaveProperties  = GetGeometry()[slave].GetValue(NEIGHBOUR_ELEMENTS).front().GetProperties();
    const Properties& MasterProperties = rMasterElement.GetProperties();
    double Eslave  = 1e9;
    if( SlaveProperties.Has(YOUNG_MODULUS) ){
	Eslave  = SlaveProperties[YOUNG_MODULUS];
    }
    else if( SlaveProperties.Has(C10) ){
	Eslave = SlaveProperties[C10];
    }

    double Emaster = 1e9;
    if( MasterProperties.Has(YOUNG_MODULUS) ){
	Emaster = MasterProperties[YOUNG_MODULUS];
    }
    else if( MasterProperties.Has(C10) ){
	Emaster = MasterProperties[C10];
    }

    //STANDARD OPTION
    if(Emaster>Eslave)
      Emaster=Eslave;

    mContactVariables.StabilizationFactor = alpha_stab/Emaster;

    //std::cout<<" Emaster "<<Emaster<<" tau "<<mContactVariables.StabilizationFactor<<std::endl;

    //EXPERIMENTAL OPTION
    // const GeometryType::IntegrationPointsArrayType& integration_points = rMasterElement.GetGeometry().IntegrationPoints( mThisIntegrationMethod );

    // //Get Current ConstitutiveMatrix
    // int size = integration_points.size();
    // std::vector<Matrix> ConstitutiveMatrix(size);
    // rMasterElement.CalculateOnIntegrationPoints(CONSTITUTIVE_MATRIX,ConstitutiveMatrix,rCurrentProcessInfo);

    // //Calc Norm of the Constitutive tensor:
    // double Cnorm=0;
    // for(int i=0; i<size; i++){
    //   for(int j=0; j<size; j++){
    // 	Cnorm += ConstitutiveMatrix[0](i,j)*ConstitutiveMatrix[0](i,j);
    //   }
    // }

    // Cnorm = sqrt(Cnorm)*0.5;


    // if(Emin>Cnorm){
    //   //std::cout<<" --Tau Stab A "<<mContactVariables.StabilizationFactor<<" -- Tau Stab B "<<alpha_stab/Cnorm<<std::endl;
    //   mContactVariables.StabilizationFactor=alpha_stab/Cnorm;
    // }



  }

  //************************************************************************************
  //************************************************************************************


  //************* COMPUTING  METHODS
  //************************************************************************************
  //************************************************************************************


  //********************************CALCULATE EXPLICIT MULTIPLIERS**********************
  //************************************************************************************

  void ContactDomainLM3DCondition::CalculateExplicitFactors(ConditionVariables& rVariables, ProcessInfo& rCurrentProcessInfo)
  {

    // std::cout<<" Master Nodes "<<GetValue(MASTER_NODES).size()<<std::endl;

    Element::ElementType& rMasterElement = GetValue(MASTER_ELEMENTS).back();
    mContactVariables.SetMasterElement(rMasterElement);


    // if( this->Is(SELECTED) )
    //   std::cout<<" ELEMENT ["<<this->Id()<<"] EDGE "<<std::endl;
    // else
    //   std::cout<<" ELEMENT ["<<this->Id()<<"] FACE "<<std::endl;

    // std::cout<<" CONTACT ("<<GetGeometry()[0].Id()<<")("<<GetGeometry()[1].Id()<<")("<<GetGeometry()[2].Id()<<")("<<GetGeometry()[3].Id()<<")"<<std::endl;
    // std::cout<<" MASTER  ("<<rMasterElement.GetGeometry()[0].Id()<<")("<<rMasterElement.GetGeometry()[1].Id()<<")("<<rMasterElement.GetGeometry()[2].Id()<<")("<<rMasterElement.GetGeometry()[3].Id()<<")"<<std::endl;
    // std::cout<<" Nodes ("<<mContactVariables.nodes[0]<<" "<<mContactVariables.nodes[1]<<" "<<mContactVariables.nodes[2]<<" "<<mContactVariables.nodes[3]<<")"<<std::endl;
    //std::cout<<" Order ("<<mContactVariables.order[0]<<" "<<mContactVariables.order[1]<<" "<<mContactVariables.order[2]<<" "<<mContactVariables.order[3]<<")"<<std::endl;

    // std::cout<<" Slaves ("<<GetValue(MASTER_NODES).front()->Id()<<" "<<GetValue(MASTER_NODES).back()->Id()<<") ["<<mContactVariables.slaves[0]<<" "<<mContactVariables.slaves[1]<<"]"<<std::endl;


    if( this->Is(SELECTED) )
      CalculateExplicitFactorsEdgeType(rVariables, rCurrentProcessInfo);
    else
      CalculateExplicitFactorsFaceType(rVariables, rCurrentProcessInfo);



    // //set contact normal
    // array_1d<double, 3> &ContactNormal  = GetGeometry()[slave].FastGetSolutionStepValue(CONTACT_NORMAL);

    // for(unsigned int i=0; i<3; i++)
    //   ContactNormal[i] = rVariables.Contact.CurrentSurface.Normal[i];


    // std::cout<<" GapN "<<rVariables.Contact.CurrentGap.Normal<<std::endl;
    // std::cout<<" LmN "<<rVariables.Contact.Multiplier.Normal<<std::endl;


    // if(rVariables.Contact.Options.Is(ACTIVE))
    //   std::cout<<"["<<this->Id()<<"] ELEMENT ACTIVE "<<std::endl;
    // else
    //   std::cout<<"["<<this->Id()<<"] ELEMENT NOT ACTIVE "<<std::endl;


    if(mContactVariables.IterationCounter < 1)
      mContactVariables.IterationCounter += 1;

  }

  //************************************************************************************
  //************************************************************************************

  void ContactDomainLM3DCondition::CalculateExplicitFactorsFaceType(ConditionVariables& rVariables, ProcessInfo& rCurrentProcessInfo)
  {

    //Contact face node1-node2-node3
    unsigned int node1=mContactVariables.nodes[0];
    unsigned int node2=mContactVariables.nodes[1];
    unsigned int node3=mContactVariables.nodes[2];
    unsigned int slave=mContactVariables.slaves.front();


    //1.- Compute tension vector:  (must be updated each iteration)
    Matrix StressMatrix(3,3);
    noalias(StressMatrix) = ZeroMatrix(3,3);

    //a.- Assign initial 2nd Piola Kirchhoff stress:
    StressMatrix=MathUtils<double>::StressVectorToTensor( rVariables.StressVector );

    // UL
    //b.- Compute the 1srt Piola Kirchhoff stress tensor
    ConstitutiveLaw Constitutive;
    StressMatrix = Constitutive.TransformStresses(StressMatrix, rVariables.F, rVariables.detF, ConstitutiveLaw::StressMeasure_Cauchy, ConstitutiveLaw::StressMeasure_PK1);

    //c.- Compute the tension (or traction) vector T=P*N (in the Reference configuration)
    mContactVariables.TractionVector=prod(StressMatrix,mContactVariables.ReferenceSurface.Normal);


    // std::cout<<" F "<<rVariables.F<<std::endl;
    // std::cout<<" StressTensor "<<StressMatrix<<std::endl;
    // std::cout<<" Reference Normal "<<mContactVariables.ReferenceSurface.Normal<<std::endl;
    // std::cout<<" Traction Vector "<<mContactVariables.TractionVector<<std::endl;


    //d.- Compute the Current Normal and Tangent

    //Current Position
    PointType P1  =  GetGeometry()[node1].Coordinates();
    PointType P2  =  GetGeometry()[node2].Coordinates();
    PointType P3  =  GetGeometry()[node3].Coordinates();
    PointType PS  =  GetGeometry()[slave].Coordinates();

    //Set Current Tangent
    rVariables.Contact.Tangent.CovariantBase.DirectionA = P2 - P1;
    rVariables.Contact.Tangent.CovariantBase.DirectionB = P3 - P1;

    TransformCovariantToContravariantBase(rVariables.Contact.Tangent.CovariantBase, rVariables.Contact.Tangent.ContravariantBase);

    // std::cout<<" cvMetric "<<rVariables.Contact.Tangent.CovariantBase.Metric<<std::endl;
    // std::cout<<" cnMetric "<<rVariables.Contact.Tangent.ContravariantBase.Metric<<std::endl;

    //Compute Current Normal
    rVariables.Contact.CurrentSurface.Normal=mContactUtilities.CalculateSurfaceNormal(rVariables.Contact.CurrentSurface.Normal,rVariables.Contact.Tangent.CovariantBase.DirectionA,rVariables.Contact.Tangent.CovariantBase.DirectionB);

    if( inner_prod(rVariables.Contact.CurrentSurface.Normal,mContactVariables.ReferenceSurface.Normal)<0 ){ //to give the correct direction
      std::cout<<" Current Contact Normal Sense INVERTED Face "<<rVariables.Contact.CurrentSurface.Normal<<std::endl;
      rVariables.Contact.CurrentSurface.Normal*=-1;
    }

    if(!(norm_2(rVariables.Contact.CurrentSurface.Normal)))
      rVariables.Contact.CurrentSurface.Normal=mContactVariables.ReferenceSurface.Normal;



    // std::cout<<" Current Normal   "<<rVariables.Contact.CurrentSurface.Normal<<std::endl;
    // std::cout<<" Reference Normal "<<mContactVariables.ReferenceSurface.Normal<<std::endl;
    // std::cout<<" PreStep Normal   "<<mContactVariables.PreStepSurface.Normal<<std::endl;

    //4.- Compute Effective Gaps: (g^eff=g_n3+2*Tau*tn=2*Tau*Multiplier.Normal)

    //Reference normal: n_n,t_n  -> mContactVariables.ReferenceSurface.Normal / rVariables.Contact.Tangent
    //Current normal:   n,t      -> rVariables.Contact.CurrentSurface.Normal /  rVariables.Contact.CurrentSurface.Tangent

    //e.- Compute A_n,B_n,L_n
    rVariables.Contact.ReferenceBase.resize(3);
    rVariables.Contact.CurrentBase.resize(3);

    //a, b, l:
    mContactUtilities.CalculateBaseDistances(rVariables.Contact.CurrentBase,P1,P2,P3,PS,rVariables.Contact.CurrentSurface.Normal);
    mContactUtilities.CalculateBaseArea(rVariables.Contact.Tangent.CurrentArea,rVariables.Contact.CurrentBase[0].L,rVariables.Contact.CurrentBase[1].L,rVariables.Contact.CurrentBase[2].L);

    //A, B, L:

    //Reference Position
    P1  =  GetGeometry()[node1].Coordinates() - (GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1) );
    P2  =  GetGeometry()[node2].Coordinates() - (GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1) );
    P3  =  GetGeometry()[node3].Coordinates() - (GetGeometry()[node3].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node3].FastGetSolutionStepValue(DISPLACEMENT,1) );
    PS  =  GetGeometry()[slave].Coordinates() - (GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,1) );

    mContactUtilities.CalculateBaseDistances(rVariables.Contact.ReferenceBase,P1,P2,P3,PS,mContactVariables.ReferenceSurface.Normal);

    //reference integration weight update
    mContactUtilities.CalculateBaseArea(rVariables.Contact.Tangent.ReferenceArea,rVariables.Contact.ReferenceBase[0].L,rVariables.Contact.ReferenceBase[1].L,rVariables.Contact.ReferenceBase[2].L);

    rVariables.Contact.Tangent.FactorArea = sqrt(rVariables.Contact.Tangent.ReferenceArea);

    // std::cout<<" L1 :"<<rVariables.Contact.ReferenceBase[0].L<<" A :"<<rVariables.Contact.ReferenceBase[0].A<<" B :"<<rVariables.Contact.ReferenceBase[0].B<<std::endl;
    // std::cout<<" L2 :"<<rVariables.Contact.ReferenceBase[1].L<<" A :"<<rVariables.Contact.ReferenceBase[1].A<<" B :"<<rVariables.Contact.ReferenceBase[1].B<<std::endl;
    // std::cout<<" L3 :"<<rVariables.Contact.ReferenceBase[2].L<<" A :"<<rVariables.Contact.ReferenceBase[2].A<<" B :"<<rVariables.Contact.ReferenceBase[2].B<<std::endl;

    //complete the computation of the stabilization gap
    rVariables.Contact.ContactFactor.Normal  =  mContactVariables.StabilizationFactor * rVariables.Contact.Tangent.FactorArea;
    rVariables.Contact.ContactFactor.Tangent =  rVariables.Contact.ContactFactor.Normal * GetProperties()[TANGENTIAL_PENALTY_RATIO];


    //f.-obtain the (g_N)3 and (g_T)3 for the n configuration

    double ReferenceGapN = inner_prod((PS-P1), mContactVariables.ReferenceSurface.Normal);

    //covariant gap
    double ReferenceGapcvTA = ReferenceGapN * inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionA,mContactVariables.ReferenceSurface.Normal);
    double ReferenceGapcvTB = ReferenceGapN * inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionB,mContactVariables.ReferenceSurface.Normal);

    //contravariant gap
    double ReferenceGapcnTA = ReferenceGapN * inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionA,mContactVariables.ReferenceSurface.Normal);
    double ReferenceGapcnTB = ReferenceGapN * inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionB,mContactVariables.ReferenceSurface.Normal);

    PointType D1  =  GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT)-GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1);
    PointType D2  =  GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT)-GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1);
    PointType D3  =  GetGeometry()[node3].FastGetSolutionStepValue(DISPLACEMENT)-GetGeometry()[node3].FastGetSolutionStepValue(DISPLACEMENT,1);
    PointType DS  =  GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT)-GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,1);

    //(g_N)3
    ReferenceGapN*=inner_prod(rVariables.Contact.CurrentSurface.Normal,mContactVariables.ReferenceSurface.Normal);
    ReferenceGapN+=inner_prod(rVariables.Contact.CurrentSurface.Normal,(D1*(-rVariables.Contact.ReferenceBase[0].A/rVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapN+=inner_prod(rVariables.Contact.CurrentSurface.Normal,(D2*(-rVariables.Contact.ReferenceBase[1].A/rVariables.Contact.ReferenceBase[1].L)));
    ReferenceGapN+=inner_prod(rVariables.Contact.CurrentSurface.Normal,(D3*(-rVariables.Contact.ReferenceBase[2].A/rVariables.Contact.ReferenceBase[2].L)));
    ReferenceGapN+=inner_prod(rVariables.Contact.CurrentSurface.Normal,DS);

    //(g_cvTA)3
    ReferenceGapcvTA+=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionA,(D1*(-rVariables.Contact.ReferenceBase[0].A/rVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapcvTA+=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionA,(D2*(-rVariables.Contact.ReferenceBase[1].A/rVariables.Contact.ReferenceBase[1].L)));
    ReferenceGapcvTA+=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionA,(D3*(-rVariables.Contact.ReferenceBase[2].A/rVariables.Contact.ReferenceBase[2].L)));
    ReferenceGapcvTA+=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionA,DS);

    //(g_cvTB)3
    ReferenceGapcvTB+=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionB,(D1*(-rVariables.Contact.ReferenceBase[0].A/rVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapcvTB+=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionB,(D2*(-rVariables.Contact.ReferenceBase[1].A/rVariables.Contact.ReferenceBase[1].L)));
    ReferenceGapcvTB+=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionB,(D3*(-rVariables.Contact.ReferenceBase[2].A/rVariables.Contact.ReferenceBase[2].L)));
    ReferenceGapcvTB+=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionB,DS);

    //(g_cnTA)3
    ReferenceGapcnTA+=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionA,(D1*(-rVariables.Contact.ReferenceBase[0].A/rVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapcnTA+=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionA,(D2*(-rVariables.Contact.ReferenceBase[1].A/rVariables.Contact.ReferenceBase[1].L)));
    ReferenceGapcnTA+=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionA,(D3*(-rVariables.Contact.ReferenceBase[2].A/rVariables.Contact.ReferenceBase[2].L)));
    ReferenceGapcnTA+=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionA,DS);

    //(g_cnTB)3
    ReferenceGapcnTB+=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionB,(D1*(-rVariables.Contact.ReferenceBase[0].A/rVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapcnTB+=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionB,(D2*(-rVariables.Contact.ReferenceBase[1].A/rVariables.Contact.ReferenceBase[1].L)));
    ReferenceGapcnTB+=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionB,(D3*(-rVariables.Contact.ReferenceBase[2].A/rVariables.Contact.ReferenceBase[2].L)));
    ReferenceGapcnTB+=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionB,DS);

    rVariables.Contact.Tangent.EquivalentHeigh = 1;

    //...

    rVariables.Contact.CurrentGap.Normal=ReferenceGapN; //(g_N)3 -- needed in the Kcont1 computation

    rVariables.Contact.Tangent.A.CurrentGap.Covariant=ReferenceGapcvTA; //(g_cvTA)3 -- needed in the Kcont1 computation
    rVariables.Contact.Tangent.A.CurrentGap.Contravariant=ReferenceGapcnTA; //(g_cvTB)3 -- needed in the Kcont1 computation

    rVariables.Contact.Tangent.B.CurrentGap.Covariant=ReferenceGapcvTB; //(g_cnTA)3 -- needed in the Kcont1 computation
    rVariables.Contact.Tangent.B.CurrentGap.Contravariant=ReferenceGapcnTB; //(g_cnTB)3 -- needed in the Kcont1 computation

    //g.- get total effective gap as: gap_n^eff=gap_n+(PreviousTimeStep/CurrentTimeStep)*(gap_n-gap_n-1)

    //gap_n   (in function of the n position of the other node) gap_n=(g_N)3+2*Tau*tn_n


    //h.- Compute normal component of the tension vector:   (tn=n·P·N)
    rVariables.Contact.CurrentTensil.Normal=inner_prod(rVariables.Contact.CurrentSurface.Normal,mContactVariables.TractionVector);

    //i.- Compute tangent component of the tension vector:  (tt=cvt·P·N)
    rVariables.Contact.Tangent.A.CurrentTensil.Covariant=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionA,mContactVariables.TractionVector);
    rVariables.Contact.Tangent.B.CurrentTensil.Covariant=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionB,mContactVariables.TractionVector);


    //j.- Compute tangent component of the tension vector:  (tt=cnt·P·N)
    rVariables.Contact.Tangent.A.CurrentTensil.Contravariant=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionA,mContactVariables.TractionVector);
    rVariables.Contact.Tangent.B.CurrentTensil.Contravariant=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionB,mContactVariables.TractionVector);


    ReferenceGapN += 3 * rVariables.Contact.ContactFactor.Normal * rVariables.Contact.CurrentTensil.Normal;

    ReferenceGapcvTA += 3 * rVariables.Contact.ContactFactor.Tangent * rVariables.Contact.Tangent.A.CurrentTensil.Covariant;
    ReferenceGapcvTB += 3 * rVariables.Contact.ContactFactor.Tangent * rVariables.Contact.Tangent.B.CurrentTensil.Covariant;

    ReferenceGapcnTA += 3 * rVariables.Contact.ContactFactor.Tangent * rVariables.Contact.Tangent.A.CurrentTensil.Contravariant;
    ReferenceGapcnTB += 3 * rVariables.Contact.ContactFactor.Tangent * rVariables.Contact.Tangent.B.CurrentTensil.Contravariant;


    //5.- Compute (Lagrange) Multipliers

    //From effective gaps set active contact domain:

    double EffectiveGapN = ReferenceGapN;

    double EffectiveGapcvTA = ReferenceGapcvTA;
    double EffectiveGapcvTB = ReferenceGapcvTB;

    double EffectiveGapcnTA = ReferenceGapcnTA;
    double EffectiveGapcnTB = ReferenceGapcnTB;


    double CurrentTimeStep  = rCurrentProcessInfo[DELTA_TIME];
    ProcessInfo& rPreviousProcessInfo = rCurrentProcessInfo.GetPreviousSolutionStepInfo();
    double PreviousTimeStep = rPreviousProcessInfo[DELTA_TIME];


    if(mContactVariables.PreviousGap.Normal!=0 && mContactVariables.IterationCounter<1){
      EffectiveGapN+=(CurrentTimeStep/PreviousTimeStep)*(ReferenceGapN-mContactVariables.PreviousGap.Normal);
    }

    if(mContactVariables.Tangent.PreviousGapA.Covariant!=0 && mContactVariables.IterationCounter<1){
      EffectiveGapcvTA+=(CurrentTimeStep/PreviousTimeStep)*(ReferenceGapcvTA-mContactVariables.Tangent.PreviousGapA.Covariant);
    }

    if(mContactVariables.Tangent.PreviousGapB.Covariant!=0 && mContactVariables.IterationCounter<1){
      EffectiveGapcvTB+=(CurrentTimeStep/PreviousTimeStep)*(ReferenceGapcvTB-mContactVariables.Tangent.PreviousGapB.Covariant);
    }

    if(mContactVariables.Tangent.PreviousGapA.Contravariant!=0 && mContactVariables.IterationCounter<1){
      EffectiveGapcnTA+=(CurrentTimeStep/PreviousTimeStep)*(ReferenceGapcnTA-mContactVariables.Tangent.PreviousGapA.Contravariant);
    }

    if(mContactVariables.Tangent.PreviousGapB.Contravariant!=0 && mContactVariables.IterationCounter<1){
      EffectiveGapcnTB+=(CurrentTimeStep/PreviousTimeStep)*(ReferenceGapcnTB-mContactVariables.Tangent.PreviousGapB.Contravariant);
    }

    //CHECK IF THE ELEMENT IS ACTIVE:

    rVariables.Contact.Options.Set(SLIP,false);

    //CORRECTION: to skip tip contact elements problems:

    bool check_fictious_geometry = false;

    //Check ORTHOGONAL FACES in contact
    if(check_fictious_geometry ==true){
      PointType& SlaveNormal  =  GetGeometry()[slave].FastGetSolutionStepValue(NORMAL);
      double orthogonal = inner_prod(SlaveNormal,rVariables.Contact.CurrentSurface.Normal);

      if(EffectiveGapN<=0 && fabs(orthogonal)<=1){

	bool real_contact = CheckFictiousContacts(rVariables);

	if(!real_contact || fabs(orthogonal)<=0.25){
	  EffectiveGapN = 1; //not active element: geometrically wrong
	  std::cout<<" DISABLE ContactElement "<<this->Id()<<" real contact "<<real_contact<<" geometrically orthogonal "<<orthogonal<<std::endl;
	}
      }
    }
    //Check ORTHOGONAL FACES in contact

    //decimal correction from tension vector calculation
    // if(fabs(EffectiveGapN)<= 1e-20 && fabs(EffectiveGapN)<= 1e-20)
    //   EffectiveGapN = 0;
    //decimal correction from tension vector calculation

    // std::cout<<" PreviousGapN  "<< mContactVariables.PreviousGap.Normal <<std::endl;
    // std::cout<<" ReferenceGapN "<< ReferenceGapN <<std::endl;
    // std::cout<<" EffectiveGapN "<< EffectiveGapN <<std::endl;
    // std::cout<<" CurrentTensil "<< rVariables.Contact.CurrentTensil.Normal <<std::endl;
    if(EffectiveGapN<=0)   //if(EffectiveGap<0){
      {

        rVariables.Contact.Options.Set(ACTIVE,true);  //normal contact active

        //Initialize friction parameter
        rVariables.Contact.FrictionCoefficient = 0;

        if( (fabs(EffectiveGapcvTA+EffectiveGapcvTB)<=rVariables.Contact.FrictionCoefficient*fabs(EffectiveGapN)) ||
	    (fabs(EffectiveGapcnTA+EffectiveGapcnTB)<=rVariables.Contact.FrictionCoefficient*fabs(EffectiveGapN)) )
	  {
	    rVariables.Contact.Options.Set(SLIP,false); //contact stick case active
	  }
        else
	  {

	    rVariables.Contact.Options.Set(SLIP,true);  //contact slip  case active
	  }

      }
    else
      {
	rVariables.Contact.Options.Set(ACTIVE,false); //normal contact not active
      }


    //temporary
    rVariables.Contact.Options.Set(ACTIVE,true);//impose normal contact active
    rVariables.Contact.Options.Set(SLIP,false); //impose stick

    //From total current gap compute multipliers:

    //rVariables.Contact.Multiplier.Normal = EffectiveGap*(1./(2.0*rVariables.Contact.ContactFactor.Normal)); //posible computation of the Lagrange Multiplier
    rVariables.Contact.Multiplier.Normal =rVariables.Contact.CurrentTensil.Normal;
    rVariables.Contact.Multiplier.Normal+=rVariables.Contact.CurrentGap.Normal*(1.0/(3.0*rVariables.Contact.ContactFactor.Normal));

    rVariables.Contact.Tangent.A.Multiplier=rVariables.Contact.Tangent.A.CurrentTensil.Covariant;
    rVariables.Contact.Tangent.A.Multiplier+=rVariables.Contact.Tangent.A.CurrentGap.Covariant*(1.0/(3.0*rVariables.Contact.ContactFactor.Tangent));

    rVariables.Contact.Tangent.B.Multiplier=rVariables.Contact.Tangent.B.CurrentTensil.Covariant;
    rVariables.Contact.Tangent.B.Multiplier+=rVariables.Contact.Tangent.B.CurrentGap.Covariant*(1.0/(3.0*rVariables.Contact.ContactFactor.Tangent));


    if(rVariables.Contact.Tangent.A.Multiplier<0)  //add the sign of the Lagrange Multiplier
      {
        rVariables.Contact.Tangent.A.GapSign*=(-1);
      }

    if(rVariables.Contact.Tangent.B.Multiplier<0)  //add the sign of the Lagrange Multiplier
      {
        rVariables.Contact.Tangent.B.GapSign*=(-1);
      }

    //check for distorted patches
    if(rVariables.Contact.Options.Is(ACTIVE)){

      double distorted_0=1,distorted_1=1,distorted_2=1,distorted_3=1;
      distorted_0=fabs((-rVariables.Contact.ReferenceBase[0].A/rVariables.Contact.ReferenceBase[0].L));
      distorted_1=fabs((-rVariables.Contact.ReferenceBase[1].A/rVariables.Contact.ReferenceBase[1].L));
      distorted_2=fabs((-rVariables.Contact.ReferenceBase[2].A/rVariables.Contact.ReferenceBase[2].L));
      distorted_3=1;

      double dist = 1.01+fabs(rVariables.Contact.CurrentGap.Normal);//1e12;

      if(distorted_0>dist || distorted_1>dist || distorted_2>dist || distorted_3>dist){
	rVariables.Contact.Options.Set(ACTIVE,false);
	std::cout<<" DISTORTED FACE ELEMENT "<<this->Id()<<" : (d0="<<distorted_0<<",d1="<<distorted_1<<",d2="<<distorted_2<<",d3="<<distorted_3<<") "<<std::endl;
      }
    }


  }


  //************************************************************************************
  //************************************************************************************

  void ContactDomainLM3DCondition::CalculateExplicitFactorsEdgeType(ConditionVariables& rVariables, ProcessInfo& rCurrentProcessInfo)
  {

    //Contact face node1-node2-node3
    unsigned int node1=mContactVariables.nodes[0];
    unsigned int node2=mContactVariables.nodes[1];
    unsigned int slave1=mContactVariables.slaves[0];
    unsigned int slave2=mContactVariables.slaves[1];

    //Condition* MasterCondition = GetValue(MASTER_CONDITION);

    //1.- Compute tension vector:  (must be updated each iteration)
    Matrix StressMatrix(3,3);
    noalias(StressMatrix) = ZeroMatrix(3,3);

    //a.- Assign initial 2nd Piola Kirchhoff stress:
    StressMatrix=MathUtils<double>::StressVectorToTensor( rVariables.StressVector );

    // UL
    //b.- Compute the 1srt Piola Kirchhoff stress tensor
    ConstitutiveLaw Constitutive;
    StressMatrix = Constitutive.TransformStresses(StressMatrix, rVariables.F, rVariables.detF, ConstitutiveLaw::StressMeasure_Cauchy, ConstitutiveLaw::StressMeasure_PK1);

    //c.- Compute the tension (or traction) vector T=P*N (in the Reference configuration)
    mContactVariables.TractionVector=prod(StressMatrix,mContactVariables.ReferenceSurface.Normal);

    // std::cout<<" F "<<rVariables.F<<std::endl;
    // std::cout<<" StressTensor "<<StressMatrix<<std::endl;
    // std::cout<<" Reference Normal "<<mContactVariables.ReferenceSurface.Normal<<std::endl;
    // std::cout<<" Traction Vector "<<mContactVariables.TractionVector<<std::endl;

    //d.- Compute the Current Normal and Tangent

    //Current Position
    PointType P1  =  GetGeometry()[node1].Coordinates();
    PointType P2  =  GetGeometry()[node2].Coordinates();
    PointType PS1 =  GetGeometry()[slave1].Coordinates();
    PointType PS2 =  GetGeometry()[slave2].Coordinates();

    //Set Current Tangent
    rVariables.Contact.Tangent.CovariantBase.DirectionA = P2 - P1;
    rVariables.Contact.Tangent.CovariantBase.DirectionB = PS1 - PS2;

    TransformCovariantToContravariantBase(rVariables.Contact.Tangent.CovariantBase, rVariables.Contact.Tangent.ContravariantBase);

    // std::cout<<" cvMetric "<<rVariables.Contact.Tangent.CovariantBase.Metric<<std::endl;
    // std::cout<<" cnMetric "<<rVariables.Contact.Tangent.ContravariantBase.Metric<<std::endl;

    //Compute Current Normal
    rVariables.Contact.CurrentSurface.Normal=mContactUtilities.CalculateSurfaceNormal(rVariables.Contact.CurrentSurface.Normal,rVariables.Contact.Tangent.CovariantBase.DirectionA, rVariables.Contact.Tangent.CovariantBase.DirectionB);


    if( inner_prod(rVariables.Contact.CurrentSurface.Normal,mContactVariables.ReferenceSurface.Normal)<0){ //to give the correct direction
      std::cout<<" Current Contact Normal Sense INVERTED Edge "<<rVariables.Contact.CurrentSurface.Normal<<std::endl;
      rVariables.Contact.CurrentSurface.Normal*=-1;
    }


    if(!(norm_2(rVariables.Contact.CurrentSurface.Normal)))
      rVariables.Contact.CurrentSurface.Normal=mContactVariables.ReferenceSurface.Normal;


    // std::cout<<" Current Normal   "<<rVariables.Contact.CurrentSurface.Normal<<std::endl;
    // std::cout<<" Reference Normal "<<mContactVariables.ReferenceSurface.Normal<<std::endl;
    // std::cout<<" PreStep Normal   "<<mContactVariables.PreStepSurface.Normal<<std::endl;

    //4.- Compute Effective Gaps: (g^eff=g_n3+2*Tau*tn=2*Tau*Multiplier.Normal)

    //Reference normal: n_n,t_n  -> mContactVariables.ReferenceSurface.Normal / rVariables.Contact.Tangent
    //Current normal:   n,t      -> rVariables.Contact.CurrentSurface.Normal /  rVariables.Contact.CurrentSurface.Tangent

    //e.- Compute A_n,B_n,L_n
    rVariables.Contact.ReferenceBase.resize(2);
    rVariables.Contact.CurrentBase.resize(2);

    //a, b, l:
    mContactUtilities.CalculateEdgeDistances(rVariables.Contact.CurrentBase,P1,P2,PS1,PS2,rVariables.Contact.CurrentSurface.Normal);
    PointType NormalDirection;
    MathUtils<double>::CrossProduct( NormalDirection, rVariables.Contact.Tangent.CovariantBase.DirectionA,rVariables.Contact.Tangent.CovariantBase.DirectionB);
    rVariables.Contact.Tangent.CurrentArea = 0.5 * norm_2(NormalDirection);


    //A, B, L:

    //Reference Position
    P1  =  GetGeometry()[node1].Coordinates() - (GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1) );
    P2  =  GetGeometry()[node2].Coordinates() - (GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1) );
    PS1  =  GetGeometry()[slave1].Coordinates() - (GetGeometry()[slave1].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[slave1].FastGetSolutionStepValue(DISPLACEMENT,1) );
    PS2  =  GetGeometry()[slave2].Coordinates() - (GetGeometry()[slave2].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[slave2].FastGetSolutionStepValue(DISPLACEMENT,1) );

    mContactUtilities.CalculateEdgeDistances(rVariables.Contact.ReferenceBase,P1,P2,PS1,PS2,mContactVariables.ReferenceSurface.Normal);

    //reference integration weight update
    PointType V1 = P2 - P1;
    PointType V2 = PS1 - PS2;

    PointType V3;
    MathUtils<double>::CrossProduct(V3,V1,V2);
    rVariables.Contact.Tangent.ReferenceArea = 0.5 * norm_2(V3);


    rVariables.Contact.Tangent.FactorArea = 0.25 * (rVariables.Contact.ReferenceBase[0].L + rVariables.Contact.ReferenceBase[1].L) * (rVariables.Contact.ReferenceBase[0].L + rVariables.Contact.ReferenceBase[1].L);


    // std::cout<<" L1 :"<<rVariables.Contact.ReferenceBase[0].L<<" A :"<<rVariables.Contact.ReferenceBase[0].A<<" B :"<<rVariables.Contact.ReferenceBase[0].B<<std::endl;
    // std::cout<<" L2 :"<<rVariables.Contact.ReferenceBase[1].L<<" A :"<<rVariables.Contact.ReferenceBase[1].A<<" B :"<<rVariables.Contact.ReferenceBase[1].B<<std::endl;


    //complete the computation of the stabilization gap
    rVariables.Contact.ContactFactor.Normal  =  mContactVariables.StabilizationFactor * rVariables.Contact.Tangent.FactorArea;
    rVariables.Contact.ContactFactor.Tangent =  rVariables.Contact.ContactFactor.Normal * GetProperties()[TANGENTIAL_PENALTY_RATIO];;

    //f.-obtain the (g_N)3 and (g_T)3 for the n configuration

    PointType M1 = 0.5 * (P1+P2);
    PointType M2 = 0.5 * (PS1+PS2);
    double ReferenceGapN = inner_prod((M2-M1), mContactVariables.ReferenceSurface.Normal);

    //covariant gap
    double ReferenceGapcvTA = ReferenceGapN * inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionA,mContactVariables.ReferenceSurface.Normal);
    double ReferenceGapcvTB = ReferenceGapN * inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionB,mContactVariables.ReferenceSurface.Normal);

    //contravariant gap
    double ReferenceGapcnTA = ReferenceGapN * inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionA,mContactVariables.ReferenceSurface.Normal);
    double ReferenceGapcnTB = ReferenceGapN * inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionB,mContactVariables.ReferenceSurface.Normal);

    PointType D1  =  GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT)-GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1);
    PointType D2  =  GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT)-GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1);
    PointType DS1  =  GetGeometry()[slave1].FastGetSolutionStepValue(DISPLACEMENT)-GetGeometry()[slave1].FastGetSolutionStepValue(DISPLACEMENT,1);
    PointType DS2  =  GetGeometry()[slave2].FastGetSolutionStepValue(DISPLACEMENT)-GetGeometry()[slave2].FastGetSolutionStepValue(DISPLACEMENT,1);

    //(g_N)3
    ReferenceGapN*=inner_prod(rVariables.Contact.CurrentSurface.Normal,mContactVariables.ReferenceSurface.Normal);

    ReferenceGapN+=inner_prod(rVariables.Contact.CurrentSurface.Normal,(D1*(-rVariables.Contact.ReferenceBase[0].A/rVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapN+=inner_prod(rVariables.Contact.CurrentSurface.Normal,(D2*(-rVariables.Contact.ReferenceBase[0].B/rVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapN+=inner_prod(rVariables.Contact.CurrentSurface.Normal,(DS1*(rVariables.Contact.ReferenceBase[1].A/rVariables.Contact.ReferenceBase[1].L)));
    ReferenceGapN+=inner_prod(rVariables.Contact.CurrentSurface.Normal,(DS2*(rVariables.Contact.ReferenceBase[1].B/rVariables.Contact.ReferenceBase[1].L)));


    //(g_cvTA)3
    ReferenceGapcvTA+=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionA,(D1*(-rVariables.Contact.ReferenceBase[0].A/rVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapcvTA+=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionA,(D2*(-rVariables.Contact.ReferenceBase[0].B/rVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapcvTA+=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionA,(DS1*(rVariables.Contact.ReferenceBase[1].A/rVariables.Contact.ReferenceBase[1].L)));
    ReferenceGapcvTA+=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionA,(DS2*(rVariables.Contact.ReferenceBase[1].B/rVariables.Contact.ReferenceBase[1].L)));

    //(g_cvTB)3
    ReferenceGapcvTB+=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionB,(D1*(-rVariables.Contact.ReferenceBase[0].A/rVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapcvTB+=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionB,(D2*(-rVariables.Contact.ReferenceBase[0].B/rVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapcvTB+=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionB,(DS1*(rVariables.Contact.ReferenceBase[1].A/rVariables.Contact.ReferenceBase[1].L)));
    ReferenceGapcvTB+=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionB,(DS2*(rVariables.Contact.ReferenceBase[1].B/rVariables.Contact.ReferenceBase[1].L)));

    //(g_cnTA)3
    ReferenceGapcnTA+=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionA,(D1*(-rVariables.Contact.ReferenceBase[0].A/rVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapcnTA+=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionA,(D2*(-rVariables.Contact.ReferenceBase[0].B/rVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapcnTA+=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionA,(DS1*(rVariables.Contact.ReferenceBase[1].A/rVariables.Contact.ReferenceBase[1].L)));
    ReferenceGapcnTA+=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionA,(DS2*(rVariables.Contact.ReferenceBase[1].B/rVariables.Contact.ReferenceBase[1].L)));

    //(g_cnTB)3
    ReferenceGapcnTB+=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionB,(D1*(-rVariables.Contact.ReferenceBase[0].A/rVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapcnTB+=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionB,(D2*(-rVariables.Contact.ReferenceBase[0].B/rVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapcnTB+=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionB,(DS1*(rVariables.Contact.ReferenceBase[1].A/rVariables.Contact.ReferenceBase[1].L)));
    ReferenceGapcnTB+=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionB,(DS2*(rVariables.Contact.ReferenceBase[1].B/rVariables.Contact.ReferenceBase[1].L)));

    rVariables.Contact.Tangent.EquivalentHeigh = rVariables.Contact.Tangent.ReferenceArea/rVariables.Contact.ReferenceBase[0].L;

    //...

    rVariables.Contact.CurrentGap.Normal=ReferenceGapN; //(g_N)3 -- needed in the Kcont1 computation

    rVariables.Contact.Tangent.A.CurrentGap.Covariant=ReferenceGapcvTA; //(g_cvTA)3 -- needed in the Kcont1 computation
    rVariables.Contact.Tangent.A.CurrentGap.Contravariant=ReferenceGapcnTA; //(g_cnTA)3 -- needed in the Kcont1 computation

    rVariables.Contact.Tangent.B.CurrentGap.Covariant=ReferenceGapcvTB; //(g_cvTB)3 -- needed in the Kcont1 computation
    rVariables.Contact.Tangent.B.CurrentGap.Contravariant=ReferenceGapcnTB; //(g_nTB)3 -- needed in the Kcont1 computation

    //g.- get total effective gap as: gap_n^eff=gap_n+(PreviousTimeStep/CurrentTimeStep)*(gap_n-gap_n-1)

    //gap_n   (in function of the n position of the other node) gap_n=(g_N)3+2*Tau*tn_n


    //h.- Compute normal component of the tension vector:   (tn=n·P·N)
    rVariables.Contact.CurrentTensil.Normal=inner_prod(rVariables.Contact.CurrentSurface.Normal,mContactVariables.TractionVector);


    //i.- Compute tangent component of the tension vector:  (tt=cvt·P·N)
    rVariables.Contact.Tangent.A.CurrentTensil.Covariant=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionA,mContactVariables.TractionVector);
    rVariables.Contact.Tangent.B.CurrentTensil.Covariant=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionB,mContactVariables.TractionVector);


    //j.- Compute tangent component of the tension vector:  (tt=cnt·P·N)
    rVariables.Contact.Tangent.A.CurrentTensil.Contravariant=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionA,mContactVariables.TractionVector);
    rVariables.Contact.Tangent.B.CurrentTensil.Contravariant=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionB,mContactVariables.TractionVector);


    ReferenceGapN += 3 * rVariables.Contact.ContactFactor.Normal * rVariables.Contact.CurrentTensil.Normal/rVariables.Contact.Tangent.EquivalentHeigh;

    ReferenceGapcvTA += 3 * rVariables.Contact.ContactFactor.Tangent * rVariables.Contact.Tangent.A.CurrentTensil.Covariant/rVariables.Contact.Tangent.EquivalentHeigh;
    ReferenceGapcvTB += 3 * rVariables.Contact.ContactFactor.Tangent * rVariables.Contact.Tangent.B.CurrentTensil.Covariant/rVariables.Contact.Tangent.EquivalentHeigh;

    ReferenceGapcnTA += 3 * rVariables.Contact.ContactFactor.Tangent * rVariables.Contact.Tangent.A.CurrentTensil.Contravariant/rVariables.Contact.Tangent.EquivalentHeigh;
    ReferenceGapcnTB += 3 * rVariables.Contact.ContactFactor.Tangent * rVariables.Contact.Tangent.B.CurrentTensil.Contravariant/rVariables.Contact.Tangent.EquivalentHeigh;


    //5.- Compute (Lagrange) Multipliers

    //From effective gaps set active contact domain:

    double EffectiveGapN = ReferenceGapN;

    double EffectiveGapcvTA = ReferenceGapcvTA;
    double EffectiveGapcvTB = ReferenceGapcvTB;

    double EffectiveGapcnTA = ReferenceGapcnTA;
    double EffectiveGapcnTB = ReferenceGapcnTB;


    double CurrentTimeStep  = rCurrentProcessInfo[DELTA_TIME];
    ProcessInfo& rPreviousProcessInfo = rCurrentProcessInfo.GetPreviousSolutionStepInfo();
    double PreviousTimeStep = rPreviousProcessInfo[DELTA_TIME];


    if(mContactVariables.PreviousGap.Normal!=0 && mContactVariables.IterationCounter<1){
      EffectiveGapN+=(CurrentTimeStep/PreviousTimeStep)*(ReferenceGapN-mContactVariables.PreviousGap.Normal);
    }

    if(mContactVariables.Tangent.PreviousGapA.Covariant!=0 && mContactVariables.IterationCounter<1){
      EffectiveGapcvTA+=(CurrentTimeStep/PreviousTimeStep)*(ReferenceGapcvTA-mContactVariables.Tangent.PreviousGapA.Covariant);
    }

    if(mContactVariables.Tangent.PreviousGapB.Covariant!=0 && mContactVariables.IterationCounter<1){
      EffectiveGapcvTB+=(CurrentTimeStep/PreviousTimeStep)*(ReferenceGapcvTB-mContactVariables.Tangent.PreviousGapB.Covariant);
    }

    if(mContactVariables.Tangent.PreviousGapA.Contravariant!=0 && mContactVariables.IterationCounter<1){
      EffectiveGapcnTA+=(CurrentTimeStep/PreviousTimeStep)*(ReferenceGapcnTA-mContactVariables.Tangent.PreviousGapA.Contravariant);
    }

    if(mContactVariables.Tangent.PreviousGapB.Contravariant!=0 && mContactVariables.IterationCounter<1){
      EffectiveGapcnTB+=(CurrentTimeStep/PreviousTimeStep)*(ReferenceGapcnTB-mContactVariables.Tangent.PreviousGapB.Contravariant);
    }

    //CHECK IF THE ELEMENT IS ACTIVE:

    rVariables.Contact.Options.Set(SLIP,false);

    //CORRECTION: to skip tip contact elements problems:

    bool check_fictious_geometry = false;

    //Check ORTHOGONAL FACES in contact
    if(check_fictious_geometry ==true){
      //slave -> slave1
      PointType& SlaveNormal  =  GetGeometry()[slave1].FastGetSolutionStepValue(NORMAL);
      double orthogonal = inner_prod(SlaveNormal,rVariables.Contact.CurrentSurface.Normal);

      if(EffectiveGapN<=0 && fabs(orthogonal)<=1){

	bool real_contact = CheckFictiousContacts(rVariables);

	if(!real_contact || fabs(orthogonal)<=0.25){
	  EffectiveGapN = 1; //not active element: geometrically wrong
	  std::cout<<" DISABLE ContactElement "<<this->Id()<<" real contact "<<real_contact<<" geometrically orthogonal "<<orthogonal<<std::endl;
	}
      }
    }
    //Check ORTHOGONAL FACES in contact

    //decimal correction from tension vector calculation
    // if(fabs(EffectiveGapN)<= 1e-15 && fabs(EffectiveGapN)<= 1e-15)
    //   EffectiveGapN = 0;
    //decimal correction from tension vector calculation

    // std::cout<<" PreviousGapN  "<< mContactVariables.PreviousGap.Normal <<std::endl;
    // std::cout<<" ReferenceGapN "<< ReferenceGapN <<std::endl;
    // std::cout<<" EffectiveGapN "<< EffectiveGapN <<std::endl;
    // std::cout<<" CurrentTensil "<< rVariables.Contact.CurrentTensil.Normal <<std::endl;
    if(EffectiveGapN<=0)   //if(EffectiveGap<0){
      {
        //Initialize friction parameter
        rVariables.Contact.FrictionCoefficient = 0;

        rVariables.Contact.Options.Set(ACTIVE, true);  //normal contact active

        if( (fabs(EffectiveGapcvTA+EffectiveGapcvTB)<=rVariables.Contact.FrictionCoefficient*fabs(EffectiveGapN)) ||
	    (fabs(EffectiveGapcnTA+EffectiveGapcnTB)<=rVariables.Contact.FrictionCoefficient*fabs(EffectiveGapN)) )
	  {
	    rVariables.Contact.Options.Set(SLIP,false); //contact stick case active
	  }
        else
	  {

	    rVariables.Contact.Options.Set(SLIP,true);  //contact slip  case active
	  }

      }
    else
      {
	rVariables.Contact.Options.Set(ACTIVE,false); //normal contact not active
      }

    //temporary
    rVariables.Contact.Options.Set(SLIP,false); //impose stick

    //From total current gap compute multipliers:

    //rVariables.Contact.Multiplier.Normal = EffectiveGap*(1./(2.0*rVariables.Contact.ContactFactor.Normal)); //posible computation of the Lagrange Multiplier
    rVariables.Contact.Multiplier.Normal =rVariables.Contact.CurrentTensil.Normal;
    rVariables.Contact.Multiplier.Normal+=rVariables.Contact.CurrentGap.Normal*(rVariables.Contact.Tangent.EquivalentHeigh/(3.0*rVariables.Contact.ContactFactor.Normal));

    rVariables.Contact.Tangent.A.Multiplier=rVariables.Contact.Tangent.A.CurrentTensil.Covariant;
    rVariables.Contact.Tangent.A.Multiplier+=rVariables.Contact.Tangent.A.CurrentGap.Covariant*(rVariables.Contact.Tangent.EquivalentHeigh/(3.0*rVariables.Contact.ContactFactor.Tangent));

    rVariables.Contact.Tangent.B.Multiplier=rVariables.Contact.Tangent.B.CurrentTensil.Covariant;
    rVariables.Contact.Tangent.B.Multiplier+=rVariables.Contact.Tangent.B.CurrentGap.Covariant*(rVariables.Contact.Tangent.EquivalentHeigh/(3.0*rVariables.Contact.ContactFactor.Tangent));


    if(rVariables.Contact.Tangent.A.Multiplier<0)  //add the sign of the Lagrange Multiplier
      {
        rVariables.Contact.Tangent.A.GapSign*=(-1);
      }

    if(rVariables.Contact.Tangent.B.Multiplier<0)  //add the sign of the Lagrange Multiplier
      {
        rVariables.Contact.Tangent.B.GapSign*=(-1);
      }

    //check for distorted patches
    if(rVariables.Contact.Options.Is(ACTIVE)){

      double distorted_0=1,distorted_1=1,distorted_2=1,distorted_3=1;
      distorted_0=fabs((-rVariables.Contact.ReferenceBase[0].A/rVariables.Contact.ReferenceBase[0].L));
      distorted_1=fabs((-rVariables.Contact.ReferenceBase[0].B/rVariables.Contact.ReferenceBase[0].L));
      distorted_2=fabs((rVariables.Contact.ReferenceBase[1].A/rVariables.Contact.ReferenceBase[1].L));
      distorted_3=fabs((rVariables.Contact.ReferenceBase[1].B/rVariables.Contact.ReferenceBase[1].L));

      double dist = 1.01+fabs(rVariables.Contact.CurrentGap.Normal);//1e12;

      if(distorted_0>dist || distorted_1>dist || distorted_2>dist || distorted_3>dist){
	rVariables.Contact.Options.Set(ACTIVE,false);
	std::cout<<" DISTORTED EDGE ELEMENT "<<this->Id()<<" : (d0="<<distorted_0<<",d1="<<distorted_1<<",d2="<<distorted_2<<",d3="<<distorted_3<<") "<<std::endl;
      }
    }

  }



  //************************************************************************************
  //************************************************************************************


  void ContactDomainLM3DCondition::CalculateDomainShapeN(ConditionVariables& rVariables)
  {

    unsigned int ndi,ndj,ndk,ndl,ndm,ndn;

    if( this->Is(SELECTED) ){

      ndi=mContactVariables.nodes[0];
      ndj=mContactVariables.nodes[1];
      ndk=mContactVariables.nodes[2];
      ndl=mContactVariables.nodes[3];
      ndm=4;
      ndn=5;

    }
    else{

      ndi=mContactVariables.nodes[0];
      ndj=mContactVariables.nodes[1];
      ndk=mContactVariables.nodes[2];
      ndl=mContactVariables.nodes[3];
      ndm=4;

    }

    //Set discrete variations of the shape function on the normal and tangent directions:

    Matrix DN_DX = rVariables.DN_DX;
    for(unsigned int i=0;i<DN_DX.size1();i++)
      {
	for(unsigned int j=0;j<DN_DX.size2();j++)
	  {
	    rVariables.DN_DX(i,j) = DN_DX(mContactVariables.order[i],j);
	  }
      }

    //std::cout<<" Variables_DN_DX "<<rVariables.DN_DX<<std::endl;

    if( this->Is(SELECTED) ){

      //dN_dn:

      rVariables.Contact.dN_dn.resize(6);
      noalias(rVariables.Contact.dN_dn) = ZeroVector(6);

      rVariables.Contact.dN_dn[ndi]=(-1)*rVariables.Contact.CurrentBase[0].A/rVariables.Contact.CurrentBase[0].L;
      rVariables.Contact.dN_dn[ndj]=(-1)*rVariables.Contact.CurrentBase[0].B/rVariables.Contact.CurrentBase[0].L;
      rVariables.Contact.dN_dn[ndk]= rVariables.Contact.CurrentBase[1].A/rVariables.Contact.CurrentBase[1].L;
      rVariables.Contact.dN_dn[ndl]= rVariables.Contact.CurrentBase[1].B/rVariables.Contact.CurrentBase[1].L;

      //std::cout<<" dN_dn "<<rVariables.Contact.dN_dn<<std::endl;

      //dN_drn:

      rVariables.Contact.dN_drn.resize(6);
      noalias(rVariables.Contact.dN_drn) = ZeroVector(6);

      rVariables.Contact.dN_drn[ndi]=(-1)*rVariables.Contact.ReferenceBase[0].A/rVariables.Contact.ReferenceBase[0].L;
      rVariables.Contact.dN_drn[ndj]=(-1)*rVariables.Contact.ReferenceBase[0].B/rVariables.Contact.ReferenceBase[0].L;
      rVariables.Contact.dN_drn[ndk]= rVariables.Contact.ReferenceBase[1].A/rVariables.Contact.ReferenceBase[1].L;
      rVariables.Contact.dN_drn[ndl]= rVariables.Contact.ReferenceBase[1].B/rVariables.Contact.ReferenceBase[1].L;;

      // std::cout<<" dN_drn "<<rVariables.Contact.dN_drn<<std::endl;

      //A.dN_dt:

      rVariables.Contact.Tangent.A.dN_dt.resize(6);
      noalias(rVariables.Contact.Tangent.A.dN_dt) = ZeroVector(6);


      rVariables.Contact.Tangent.A.dN_dt[ndi]=-1.0;
      rVariables.Contact.Tangent.A.dN_dt[ndj]= 1.0;

      //B.dN_dt:

      rVariables.Contact.Tangent.B.dN_dt.resize(6);
      noalias(rVariables.Contact.Tangent.B.dN_dt) = ZeroVector(6);

      rVariables.Contact.Tangent.B.dN_dt[ndk]= 1.0;
      rVariables.Contact.Tangent.B.dN_dt[ndl]=-1.0;

      //A.TsigmaP : 1-2-5-6
      rVariables.Contact.Tangent.A.Tsigma.resize(6);

      //B.TsigmaP :
      rVariables.Contact.Tangent.B.Tsigma.resize(6);

      //NsigmaP :
      rVariables.Contact.Nsigma.resize(6);

    }
    else{

      //dN_dn:
      rVariables.Contact.dN_dn.resize(5);
      noalias(rVariables.Contact.dN_dn) = ZeroVector(5);

      rVariables.Contact.dN_dn[ndi]=(-1)*rVariables.Contact.CurrentBase[0].A/rVariables.Contact.CurrentBase[0].L;
      rVariables.Contact.dN_dn[ndj]=(-1)*rVariables.Contact.CurrentBase[1].A/rVariables.Contact.CurrentBase[1].L;
      rVariables.Contact.dN_dn[ndk]=(-1)*rVariables.Contact.CurrentBase[2].A/rVariables.Contact.CurrentBase[2].L;
      rVariables.Contact.dN_dn[ndl]= 1.0;

      //std::cout<<" dN_dn "<<rVariables.Contact.dN_dn<<std::endl;

      //dN_drn:

      rVariables.Contact.dN_drn.resize(5);
      noalias(rVariables.Contact.dN_drn) = ZeroVector(5);

      rVariables.Contact.dN_drn[ndi]=(-1)*rVariables.Contact.ReferenceBase[0].A/rVariables.Contact.ReferenceBase[0].L;
      rVariables.Contact.dN_drn[ndj]=(-1)*rVariables.Contact.ReferenceBase[1].A/rVariables.Contact.ReferenceBase[1].L;
      rVariables.Contact.dN_drn[ndk]=(-1)*rVariables.Contact.ReferenceBase[2].A/rVariables.Contact.ReferenceBase[2].L;
      rVariables.Contact.dN_drn[ndl]= 1.0;

      // std::cout<<" dN_drn "<<rVariables.Contact.dN_drn<<std::endl;

      //A.dN_dt:

      rVariables.Contact.Tangent.A.dN_dt.resize(5);
      rVariables.Contact.Tangent.A.dN_dt.clear();

      rVariables.Contact.Tangent.A.dN_dt[ndi]=-1.0;
      rVariables.Contact.Tangent.A.dN_dt[ndj]= 1.0;

      //B.dN_dt:

      rVariables.Contact.Tangent.B.dN_dt.resize(5);
      rVariables.Contact.Tangent.B.dN_dt.clear();

      rVariables.Contact.Tangent.B.dN_dt[ndi]=-1.0;
      rVariables.Contact.Tangent.B.dN_dt[ndk]= 1.0;

      //A.TsigmaP : 1-3-2-5
      rVariables.Contact.Tangent.A.Tsigma.resize(5);

      //B.TsigmaP :
      rVariables.Contact.Tangent.B.Tsigma.resize(5);

      //NsigmaP :
      rVariables.Contact.Nsigma.resize(5);

    }

    //rVariables.Contact.CurTangent=(rVariables.Contact.Tangent.cvta+rVariables.Contact.Tangent.cvtb).direction();

    // std::cout<<" Contact Factor "<< rVariables.Contact.ContactFactor.Normal <<std::endl;
    // std::cout<<" Tau "<<mContactVariables.StabilizationFactor<<std::endl;
    // std::cout<<" Area "<<rVariables.Contact.Tangent.CurrentArea<<std::endl;
    // std::cout<<" FactorArea "<<rVariables.Contact.Tangent.FactorArea<<std::endl;
    // std::cout<<" CHeigh "<<rVariables.Contact.Tangent.EquivalentHeigh<<std::endl;
    // std::cout<<" cvTab "<<rVariables.Contact.Tangent.CovariantBase.Metric<<std::endl;
    // std::cout<<" cnDirectionA "<<rVariables.Contact.Tangent.ContravariantBase.DirectionA<<std::endl;
    // std::cout<<" cvDirectionA "<<rVariables.Contact.Tangent.CovariantBase.DirectionA<<std::endl;

    //A.TsigmaP :
    FSigmaP(rVariables,rVariables.Contact.Tangent.A.Tsigma,rVariables.Contact.Tangent.CovariantBase.DirectionA,ndi,ndj,ndk,ndl,ndm,ndn);

    // std::cout<<" cnDirectionB "<<rVariables.Contact.Tangent.ContravariantBase.DirectionB<<std::endl;
    // std::cout<<" cvDirectionB "<<rVariables.Contact.Tangent.CovariantBase.DirectionB<<std::endl;

    //B.TsigmaP :
    FSigmaP(rVariables,rVariables.Contact.Tangent.B.Tsigma,rVariables.Contact.Tangent.CovariantBase.DirectionB,ndi,ndj,ndk,ndl,ndm,ndn);

    // std::cout<<" NormalDirection "<<rVariables.Contact.CurrentSurface.Normal<<std::endl;

    //NsigmaP :
    FSigmaP(rVariables,rVariables.Contact.Nsigma,rVariables.Contact.CurrentSurface.Normal,ndi,ndj,ndk,ndl,ndm,ndn);


    rVariables.DN_DX = DN_DX;


    // unsigned int size = 5;
    // unsigned int dimension = 3;

    // if( this->Is(SELECTED) )
    //   size = 6;

    // double fcont = 0;
    // for (unsigned int ndi=0; ndi<size; ndi++)
    //   {
    // 	for (unsigned int i=0; i<dimension; i++)
    // 	  {
    // 	    this->CalculateNormalForce(fcont,rVariables,ndi,i);
    // 	    std::cout<<" rfcont "<<fcont * (1.0/3.0) * rVariables.Contact.Tangent.ReferenceArea<<" ndi "<<ndi<<" i "<<i<<std::endl;
    // 	    std::cout<<" cfcont "<<fcont * (1.0/3.0) * rVariables.Contact.Tangent.CurrentArea<<" ndi "<<ndi<<" i "<<i<<std::endl;
    // 	  }
    //   }

    // double kcont = 0;
    // for (unsigned int ndi=0; ndi<size; ndi++)
    //   {
    // 	for (unsigned int ndj=0; ndj<size; ndj++)
    // 	  {
    // 	    for (unsigned int i=0; i<dimension; i++)
    // 	      {
    // 		for (unsigned int j=0; j<dimension; j++)
    // 		  {
    // 		    kcont=0;
    // 		    this->CalculateContactStiffness(kcont,rVariables,ndi,ndj,i,j);
    // 		    std::cout<<" kcont "<<kcont * (1.0/3.0) * rVariables.Contact.Tangent.CurrentArea<<" ndi "<<ndi<<" ndj "<<ndj<<" i "<<i<<" j "<<j<<std::endl;
    // 		  }
    // 	      }
    // 	  }
    //   }



  }


  //************************************************************************************
  //************************************************************************************

  void ContactDomainLM3DCondition::FSigmaP(ConditionVariables& rVariables, std::vector<Vector > &rSigmaP, PointType& rDirVector,unsigned int &ndi,unsigned int &ndj,unsigned int &ndk,unsigned int &ndl,unsigned int &ndm,unsigned int &ndn)
  {
    if( this->Is(SELECTED) ){

      //node for computation / node for assignment (nd1,nd2,nd5,nd6)

      FSigmaPnd(rVariables, rSigmaP, rDirVector, ndi, ndi);

      FSigmaPnd(rVariables, rSigmaP, rDirVector, ndj, ndj);

      rSigmaP[ndk] = ZeroVector(3);

      rSigmaP[ndl] = ZeroVector(3);

      FSigmaPnd(rVariables, rSigmaP, rDirVector, ndk, ndm);

      FSigmaPnd(rVariables, rSigmaP, rDirVector, ndl, ndn);

    }
    else{

      //node for computation / node for assignment (nd1,nd3,nd2,nd5)

      FSigmaPnd(rVariables, rSigmaP, rDirVector, ndi, ndi);

      FSigmaPnd(rVariables, rSigmaP, rDirVector, ndj, ndj);

      FSigmaPnd(rVariables, rSigmaP, rDirVector, ndk, ndk);

      rSigmaP[ndl] = ZeroVector(3);

      FSigmaPnd(rVariables, rSigmaP, rDirVector, ndl, ndm);

    }

    // if( this->Is(SELECTED) )
    // 	std::cout<<" CONTACT EDGE("<<GetGeometry()[0].Id()<<")("<<GetGeometry()[1].Id()<<")("<<GetGeometry()[2].Id()<<")("<<GetGeometry()[3].Id()<<") slaves: "<<mContactVariables.slaves[0]<<", "<<mContactVariables.slaves[1]<<std::endl;
    // else
    // 	std::cout<<" CONTACT ("<<GetGeometry()[0].Id()<<")("<<GetGeometry()[1].Id()<<")("<<GetGeometry()[2].Id()<<")("<<GetGeometry()[3].Id()<<") slave:"<<mContactVariables.slaves[0]<<std::endl;

    // std::cout<<" ORDER ("<<GetGeometry()[ndi].Id()<<")("<<GetGeometry()[ndj].Id()<<")("<<GetGeometry()[ndk].Id()<<")("<<GetGeometry()[ndl].Id()<<")"<<std::endl;

    // for(unsigned int i=0;i<rSigmaP.size();i++)
    // {
    // 	std::cout<<" position: "<<i<<" ";
    // 	std::cout<<rSigmaP[i]<<std::endl;
    // }


  }

  //************************************************************************************
  //************************************************************************************

  void ContactDomainLM3DCondition::FSigmaPnd(ConditionVariables& rVariables, std::vector<Vector >& rSigmaP, PointType& rDirVector,unsigned int &ndi,unsigned int &ndj)
  {
    //Computation with the ndi and storage to ndj
    rSigmaP[ndj].resize(3);
    noalias(rSigmaP[ndj]) = ZeroVector(3);

    //simplify nomenclature
    PointType& Normal = mContactVariables.ReferenceSurface.Normal;

    //n0N0*(S00*Ni0+S01*Ni1+S02*Ni2)+n0N1*(S11*Ni1+S01*Ni0+S12*Ni2)+n0N2*(S22*Ni2+S12*Ni1+S02*Ni0)
    //part1:
    rSigmaP[ndj][0]= rDirVector[0]*Normal[0]*(rVariables.StressVector[0]*rVariables.DN_DX(ndi,0)+rVariables.StressVector[3]*rVariables.DN_DX(ndi,1)+rVariables.StressVector[5]*rVariables.DN_DX(ndi,2))+
      rDirVector[0]*Normal[1]*(rVariables.StressVector[3]*rVariables.DN_DX(ndi,0)+rVariables.StressVector[1]*rVariables.DN_DX(ndi,1)+rVariables.StressVector[4]*rVariables.DN_DX(ndi,2))+
      rDirVector[0]*Normal[2]*(rVariables.StressVector[5]*rVariables.DN_DX(ndi,0)+rVariables.StressVector[4]*rVariables.DN_DX(ndi,1)+rVariables.StressVector[2]*rVariables.DN_DX(ndi,2));
    //n1N1*(S11*Ni1+S01*Ni0+S12*Ni2)+n1N0*(S00*Ni0+S01*Ni1+S02*Ni2)+n1N2*(S22*Ni2+S12*Ni1+S02*Ni0)
    rSigmaP[ndj][1]= rDirVector[1]*Normal[1]*(rVariables.StressVector[3]*rVariables.DN_DX(ndi,0)+rVariables.StressVector[1]*rVariables.DN_DX(ndi,1)+rVariables.StressVector[4]*rVariables.DN_DX(ndi,2))+
      rDirVector[1]*Normal[0]*(rVariables.StressVector[0]*rVariables.DN_DX(ndi,0)+rVariables.StressVector[3]*rVariables.DN_DX(ndi,1)+rVariables.StressVector[5]*rVariables.DN_DX(ndi,2))+
      rDirVector[1]*Normal[2]*(rVariables.StressVector[5]*rVariables.DN_DX(ndi,0)+rVariables.StressVector[4]*rVariables.DN_DX(ndi,1)+rVariables.StressVector[2]*rVariables.DN_DX(ndi,2));
    //n2N2*(S22*Ni2+S12*Ni1+S02*Ni0)+n2N1*(S11*Ni1+S01*Ni0+S12*Ni2)+n2N0*(S00*Ni0+S01*Ni1+S02*Ni2)
    rSigmaP[ndj][2]= rDirVector[2]*Normal[2]*(rVariables.StressVector[5]*rVariables.DN_DX(ndi,0)+rVariables.StressVector[4]*rVariables.DN_DX(ndi,1)+rVariables.StressVector[2]*rVariables.DN_DX(ndi,2))+
      rDirVector[2]*Normal[1]*(rVariables.StressVector[3]*rVariables.DN_DX(ndi,0)+rVariables.StressVector[1]*rVariables.DN_DX(ndi,1)+rVariables.StressVector[4]*rVariables.DN_DX(ndi,2))+
      rDirVector[2]*Normal[0]*(rVariables.StressVector[0]*rVariables.DN_DX(ndi,0)+rVariables.StressVector[3]*rVariables.DN_DX(ndi,1)+rVariables.StressVector[5]*rVariables.DN_DX(ndi,2));


    //part2:
    std::vector<Vector> FD(9);

    FD[0].resize(9);
    FD[0][0]=(rVariables.F(0,0)*rVariables.ConstitutiveMatrix(0,0)+rVariables.F(0,1)*rVariables.ConstitutiveMatrix(3,0)+rVariables.F(0,2)*rVariables.ConstitutiveMatrix(5,0));
    FD[0][1]=(rVariables.F(0,0)*rVariables.ConstitutiveMatrix(0,1)+rVariables.F(0,1)*rVariables.ConstitutiveMatrix(3,1)+rVariables.F(0,2)*rVariables.ConstitutiveMatrix(5,1));
    FD[0][2]=(rVariables.F(0,0)*rVariables.ConstitutiveMatrix(0,2)+rVariables.F(0,1)*rVariables.ConstitutiveMatrix(3,2)+rVariables.F(0,2)*rVariables.ConstitutiveMatrix(5,2));
    FD[0][3]=(rVariables.F(0,0)*rVariables.ConstitutiveMatrix(0,3)+rVariables.F(0,1)*rVariables.ConstitutiveMatrix(3,3)+rVariables.F(0,2)*rVariables.ConstitutiveMatrix(5,3));
    FD[0][4]=(rVariables.F(0,0)*rVariables.ConstitutiveMatrix(0,3)+rVariables.F(0,1)*rVariables.ConstitutiveMatrix(3,3)+rVariables.F(0,2)*rVariables.ConstitutiveMatrix(5,3));
    FD[0][5]=(rVariables.F(0,0)*rVariables.ConstitutiveMatrix(0,4)+rVariables.F(0,1)*rVariables.ConstitutiveMatrix(3,4)+rVariables.F(0,2)*rVariables.ConstitutiveMatrix(5,4));
    FD[0][6]=(rVariables.F(0,0)*rVariables.ConstitutiveMatrix(0,4)+rVariables.F(0,1)*rVariables.ConstitutiveMatrix(3,4)+rVariables.F(0,2)*rVariables.ConstitutiveMatrix(5,4));
    FD[0][7]=(rVariables.F(0,0)*rVariables.ConstitutiveMatrix(0,5)+rVariables.F(0,1)*rVariables.ConstitutiveMatrix(3,5)+rVariables.F(0,2)*rVariables.ConstitutiveMatrix(5,5));
    FD[0][8]=(rVariables.F(0,0)*rVariables.ConstitutiveMatrix(0,5)+rVariables.F(0,1)*rVariables.ConstitutiveMatrix(3,5)+rVariables.F(0,2)*rVariables.ConstitutiveMatrix(5,5));

    FD[1].resize(9);
    FD[1][0]=(rVariables.F(1,1)*rVariables.ConstitutiveMatrix(1,0)+rVariables.F(1,0)*rVariables.ConstitutiveMatrix(3,0)+rVariables.F(1,2)*rVariables.ConstitutiveMatrix(4,0));
    FD[1][1]=(rVariables.F(1,1)*rVariables.ConstitutiveMatrix(1,1)+rVariables.F(1,0)*rVariables.ConstitutiveMatrix(3,1)+rVariables.F(1,2)*rVariables.ConstitutiveMatrix(4,1));
    FD[1][2]=(rVariables.F(1,1)*rVariables.ConstitutiveMatrix(1,2)+rVariables.F(1,0)*rVariables.ConstitutiveMatrix(3,2)+rVariables.F(1,2)*rVariables.ConstitutiveMatrix(4,2));
    FD[1][3]=(rVariables.F(1,1)*rVariables.ConstitutiveMatrix(1,3)+rVariables.F(1,0)*rVariables.ConstitutiveMatrix(3,3)+rVariables.F(1,2)*rVariables.ConstitutiveMatrix(4,3));
    FD[1][4]=(rVariables.F(1,1)*rVariables.ConstitutiveMatrix(1,3)+rVariables.F(1,0)*rVariables.ConstitutiveMatrix(3,3)+rVariables.F(1,2)*rVariables.ConstitutiveMatrix(4,3));
    FD[1][5]=(rVariables.F(1,1)*rVariables.ConstitutiveMatrix(1,4)+rVariables.F(1,0)*rVariables.ConstitutiveMatrix(3,4)+rVariables.F(1,2)*rVariables.ConstitutiveMatrix(4,4));
    FD[1][6]=(rVariables.F(1,1)*rVariables.ConstitutiveMatrix(1,4)+rVariables.F(1,0)*rVariables.ConstitutiveMatrix(3,4)+rVariables.F(1,2)*rVariables.ConstitutiveMatrix(4,4));
    FD[1][7]=(rVariables.F(1,1)*rVariables.ConstitutiveMatrix(1,5)+rVariables.F(1,0)*rVariables.ConstitutiveMatrix(3,5)+rVariables.F(1,2)*rVariables.ConstitutiveMatrix(4,5));
    FD[1][8]=(rVariables.F(1,1)*rVariables.ConstitutiveMatrix(1,5)+rVariables.F(1,0)*rVariables.ConstitutiveMatrix(3,5)+rVariables.F(1,2)*rVariables.ConstitutiveMatrix(4,5));

    FD[2].resize(9);
    FD[2][0]=(rVariables.F(2,2)*rVariables.ConstitutiveMatrix(2,0)+rVariables.F(2,1)*rVariables.ConstitutiveMatrix(4,0)+rVariables.F(2,0)*rVariables.ConstitutiveMatrix(5,0));
    FD[2][1]=(rVariables.F(2,2)*rVariables.ConstitutiveMatrix(2,1)+rVariables.F(2,1)*rVariables.ConstitutiveMatrix(4,1)+rVariables.F(2,0)*rVariables.ConstitutiveMatrix(5,1));
    FD[2][2]=(rVariables.F(2,2)*rVariables.ConstitutiveMatrix(2,2)+rVariables.F(2,1)*rVariables.ConstitutiveMatrix(4,2)+rVariables.F(2,0)*rVariables.ConstitutiveMatrix(5,2));
    FD[2][3]=(rVariables.F(2,2)*rVariables.ConstitutiveMatrix(2,3)+rVariables.F(2,1)*rVariables.ConstitutiveMatrix(4,3)+rVariables.F(2,0)*rVariables.ConstitutiveMatrix(5,3));
    FD[2][4]=(rVariables.F(2,2)*rVariables.ConstitutiveMatrix(2,3)+rVariables.F(2,1)*rVariables.ConstitutiveMatrix(4,3)+rVariables.F(2,0)*rVariables.ConstitutiveMatrix(5,3));
    FD[2][5]=(rVariables.F(2,2)*rVariables.ConstitutiveMatrix(2,4)+rVariables.F(2,1)*rVariables.ConstitutiveMatrix(4,4)+rVariables.F(2,0)*rVariables.ConstitutiveMatrix(5,4));
    FD[2][6]=(rVariables.F(2,2)*rVariables.ConstitutiveMatrix(2,4)+rVariables.F(2,1)*rVariables.ConstitutiveMatrix(4,4)+rVariables.F(2,0)*rVariables.ConstitutiveMatrix(5,4));
    FD[2][7]=(rVariables.F(2,2)*rVariables.ConstitutiveMatrix(2,5)+rVariables.F(2,1)*rVariables.ConstitutiveMatrix(4,5)+rVariables.F(2,0)*rVariables.ConstitutiveMatrix(5,5));
    FD[2][8]=(rVariables.F(2,2)*rVariables.ConstitutiveMatrix(2,5)+rVariables.F(2,1)*rVariables.ConstitutiveMatrix(4,5)+rVariables.F(2,0)*rVariables.ConstitutiveMatrix(5,5));

    FD[3].resize(9);
    FD[3][0]=(rVariables.F(0,1)*rVariables.ConstitutiveMatrix(1,0)+rVariables.F(0,0)*rVariables.ConstitutiveMatrix(3,0)+rVariables.F(0,2)*rVariables.ConstitutiveMatrix(4,0));
    FD[3][1]=(rVariables.F(0,1)*rVariables.ConstitutiveMatrix(1,1)+rVariables.F(0,0)*rVariables.ConstitutiveMatrix(3,1)+rVariables.F(0,2)*rVariables.ConstitutiveMatrix(4,1));
    FD[3][2]=(rVariables.F(0,1)*rVariables.ConstitutiveMatrix(1,2)+rVariables.F(0,0)*rVariables.ConstitutiveMatrix(3,2)+rVariables.F(0,2)*rVariables.ConstitutiveMatrix(4,2));
    FD[3][3]=(rVariables.F(0,1)*rVariables.ConstitutiveMatrix(1,3)+rVariables.F(0,0)*rVariables.ConstitutiveMatrix(3,3)+rVariables.F(0,2)*rVariables.ConstitutiveMatrix(4,3));
    FD[3][4]=(rVariables.F(0,1)*rVariables.ConstitutiveMatrix(1,3)+rVariables.F(0,0)*rVariables.ConstitutiveMatrix(3,3)+rVariables.F(0,2)*rVariables.ConstitutiveMatrix(4,3));
    FD[3][5]=(rVariables.F(0,1)*rVariables.ConstitutiveMatrix(1,4)+rVariables.F(0,0)*rVariables.ConstitutiveMatrix(3,4)+rVariables.F(0,2)*rVariables.ConstitutiveMatrix(4,4));
    FD[3][6]=(rVariables.F(0,1)*rVariables.ConstitutiveMatrix(1,4)+rVariables.F(0,0)*rVariables.ConstitutiveMatrix(3,4)+rVariables.F(0,2)*rVariables.ConstitutiveMatrix(4,4));
    FD[3][7]=(rVariables.F(0,1)*rVariables.ConstitutiveMatrix(1,5)+rVariables.F(0,0)*rVariables.ConstitutiveMatrix(3,5)+rVariables.F(0,2)*rVariables.ConstitutiveMatrix(4,5));
    FD[3][8]=(rVariables.F(0,1)*rVariables.ConstitutiveMatrix(1,5)+rVariables.F(0,0)*rVariables.ConstitutiveMatrix(3,5)+rVariables.F(0,2)*rVariables.ConstitutiveMatrix(4,5));

    FD[4].resize(9);
    FD[4][0]=(rVariables.F(1,0)*rVariables.ConstitutiveMatrix(0,0)+rVariables.F(1,1)*rVariables.ConstitutiveMatrix(3,0)+rVariables.F(1,2)*rVariables.ConstitutiveMatrix(5,0));
    FD[4][1]=(rVariables.F(1,0)*rVariables.ConstitutiveMatrix(0,1)+rVariables.F(1,1)*rVariables.ConstitutiveMatrix(3,1)+rVariables.F(1,2)*rVariables.ConstitutiveMatrix(5,1));
    FD[4][2]=(rVariables.F(1,0)*rVariables.ConstitutiveMatrix(0,2)+rVariables.F(1,1)*rVariables.ConstitutiveMatrix(3,2)+rVariables.F(1,2)*rVariables.ConstitutiveMatrix(5,2));
    FD[4][3]=(rVariables.F(1,0)*rVariables.ConstitutiveMatrix(0,3)+rVariables.F(1,1)*rVariables.ConstitutiveMatrix(3,3)+rVariables.F(1,2)*rVariables.ConstitutiveMatrix(5,3));
    FD[4][4]=(rVariables.F(1,0)*rVariables.ConstitutiveMatrix(0,3)+rVariables.F(1,1)*rVariables.ConstitutiveMatrix(3,3)+rVariables.F(1,2)*rVariables.ConstitutiveMatrix(5,3));
    FD[4][5]=(rVariables.F(1,0)*rVariables.ConstitutiveMatrix(0,4)+rVariables.F(1,1)*rVariables.ConstitutiveMatrix(3,4)+rVariables.F(1,2)*rVariables.ConstitutiveMatrix(5,4));
    FD[4][6]=(rVariables.F(1,0)*rVariables.ConstitutiveMatrix(0,4)+rVariables.F(1,1)*rVariables.ConstitutiveMatrix(3,4)+rVariables.F(1,2)*rVariables.ConstitutiveMatrix(5,4));
    FD[4][7]=(rVariables.F(1,0)*rVariables.ConstitutiveMatrix(0,5)+rVariables.F(1,1)*rVariables.ConstitutiveMatrix(3,5)+rVariables.F(1,2)*rVariables.ConstitutiveMatrix(5,5));
    FD[4][8]=(rVariables.F(1,0)*rVariables.ConstitutiveMatrix(0,5)+rVariables.F(1,1)*rVariables.ConstitutiveMatrix(3,5)+rVariables.F(1,2)*rVariables.ConstitutiveMatrix(5,5));

    FD[5].resize(9);
    FD[5][0]=(rVariables.F(1,2)*rVariables.ConstitutiveMatrix(2,0)+rVariables.F(1,1)*rVariables.ConstitutiveMatrix(4,0)+rVariables.F(1,0)*rVariables.ConstitutiveMatrix(5,0));
    FD[5][1]=(rVariables.F(1,2)*rVariables.ConstitutiveMatrix(2,1)+rVariables.F(1,1)*rVariables.ConstitutiveMatrix(4,1)+rVariables.F(1,0)*rVariables.ConstitutiveMatrix(5,1));
    FD[5][2]=(rVariables.F(1,2)*rVariables.ConstitutiveMatrix(2,2)+rVariables.F(1,1)*rVariables.ConstitutiveMatrix(4,2)+rVariables.F(1,0)*rVariables.ConstitutiveMatrix(5,2));
    FD[5][3]=(rVariables.F(1,2)*rVariables.ConstitutiveMatrix(2,3)+rVariables.F(1,1)*rVariables.ConstitutiveMatrix(4,3)+rVariables.F(1,0)*rVariables.ConstitutiveMatrix(5,3));
    FD[5][4]=(rVariables.F(1,2)*rVariables.ConstitutiveMatrix(2,3)+rVariables.F(1,1)*rVariables.ConstitutiveMatrix(4,3)+rVariables.F(1,0)*rVariables.ConstitutiveMatrix(5,3));
    FD[5][5]=(rVariables.F(1,2)*rVariables.ConstitutiveMatrix(2,4)+rVariables.F(1,1)*rVariables.ConstitutiveMatrix(4,4)+rVariables.F(1,0)*rVariables.ConstitutiveMatrix(5,4));
    FD[5][6]=(rVariables.F(1,2)*rVariables.ConstitutiveMatrix(2,4)+rVariables.F(1,1)*rVariables.ConstitutiveMatrix(4,4)+rVariables.F(1,0)*rVariables.ConstitutiveMatrix(5,4));
    FD[5][7]=(rVariables.F(1,2)*rVariables.ConstitutiveMatrix(2,5)+rVariables.F(1,1)*rVariables.ConstitutiveMatrix(4,5)+rVariables.F(1,0)*rVariables.ConstitutiveMatrix(5,5));
    FD[5][8]=(rVariables.F(1,2)*rVariables.ConstitutiveMatrix(2,5)+rVariables.F(1,1)*rVariables.ConstitutiveMatrix(4,5)+rVariables.F(1,0)*rVariables.ConstitutiveMatrix(5,5));

    FD[6].resize(9);
    FD[6][0]=(rVariables.F(2,1)*rVariables.ConstitutiveMatrix(1,0)+rVariables.F(2,0)*rVariables.ConstitutiveMatrix(3,0)+rVariables.F(2,2)*rVariables.ConstitutiveMatrix(4,0));
    FD[6][1]=(rVariables.F(2,1)*rVariables.ConstitutiveMatrix(1,1)+rVariables.F(2,0)*rVariables.ConstitutiveMatrix(3,1)+rVariables.F(2,2)*rVariables.ConstitutiveMatrix(4,1));
    FD[6][2]=(rVariables.F(2,1)*rVariables.ConstitutiveMatrix(1,2)+rVariables.F(2,0)*rVariables.ConstitutiveMatrix(3,2)+rVariables.F(2,2)*rVariables.ConstitutiveMatrix(4,2));
    FD[6][3]=(rVariables.F(2,1)*rVariables.ConstitutiveMatrix(1,3)+rVariables.F(2,0)*rVariables.ConstitutiveMatrix(3,3)+rVariables.F(2,2)*rVariables.ConstitutiveMatrix(4,3));
    FD[6][4]=(rVariables.F(2,1)*rVariables.ConstitutiveMatrix(1,3)+rVariables.F(2,0)*rVariables.ConstitutiveMatrix(3,3)+rVariables.F(2,2)*rVariables.ConstitutiveMatrix(4,3));
    FD[6][5]=(rVariables.F(2,1)*rVariables.ConstitutiveMatrix(1,4)+rVariables.F(2,0)*rVariables.ConstitutiveMatrix(3,4)+rVariables.F(2,2)*rVariables.ConstitutiveMatrix(4,4));
    FD[6][6]=(rVariables.F(2,1)*rVariables.ConstitutiveMatrix(1,4)+rVariables.F(2,0)*rVariables.ConstitutiveMatrix(3,4)+rVariables.F(2,2)*rVariables.ConstitutiveMatrix(4,4));
    FD[6][7]=(rVariables.F(2,1)*rVariables.ConstitutiveMatrix(1,5)+rVariables.F(2,0)*rVariables.ConstitutiveMatrix(3,5)+rVariables.F(2,2)*rVariables.ConstitutiveMatrix(4,5));
    FD[6][8]=(rVariables.F(2,1)*rVariables.ConstitutiveMatrix(1,5)+rVariables.F(2,0)*rVariables.ConstitutiveMatrix(3,5)+rVariables.F(2,2)*rVariables.ConstitutiveMatrix(4,5));

    FD[7].resize(9);
    FD[7][0]=(rVariables.F(0,2)*rVariables.ConstitutiveMatrix(2,0)+rVariables.F(0,1)*rVariables.ConstitutiveMatrix(4,0)+rVariables.F(0,0)*rVariables.ConstitutiveMatrix(5,0));
    FD[7][1]=(rVariables.F(0,2)*rVariables.ConstitutiveMatrix(2,1)+rVariables.F(0,1)*rVariables.ConstitutiveMatrix(4,1)+rVariables.F(0,0)*rVariables.ConstitutiveMatrix(5,1));
    FD[7][2]=(rVariables.F(0,2)*rVariables.ConstitutiveMatrix(2,2)+rVariables.F(0,1)*rVariables.ConstitutiveMatrix(4,2)+rVariables.F(0,0)*rVariables.ConstitutiveMatrix(5,2));
    FD[7][3]=(rVariables.F(0,2)*rVariables.ConstitutiveMatrix(2,3)+rVariables.F(0,1)*rVariables.ConstitutiveMatrix(4,3)+rVariables.F(0,0)*rVariables.ConstitutiveMatrix(5,3));
    FD[7][4]=(rVariables.F(0,2)*rVariables.ConstitutiveMatrix(2,3)+rVariables.F(0,1)*rVariables.ConstitutiveMatrix(4,3)+rVariables.F(0,0)*rVariables.ConstitutiveMatrix(5,3));
    FD[7][5]=(rVariables.F(0,2)*rVariables.ConstitutiveMatrix(2,4)+rVariables.F(0,1)*rVariables.ConstitutiveMatrix(4,4)+rVariables.F(0,0)*rVariables.ConstitutiveMatrix(5,4));
    FD[7][6]=(rVariables.F(0,2)*rVariables.ConstitutiveMatrix(2,4)+rVariables.F(0,1)*rVariables.ConstitutiveMatrix(4,4)+rVariables.F(0,0)*rVariables.ConstitutiveMatrix(5,4));
    FD[7][7]=(rVariables.F(0,2)*rVariables.ConstitutiveMatrix(2,5)+rVariables.F(0,1)*rVariables.ConstitutiveMatrix(4,5)+rVariables.F(0,0)*rVariables.ConstitutiveMatrix(5,5));
    FD[7][8]=(rVariables.F(0,2)*rVariables.ConstitutiveMatrix(2,5)+rVariables.F(0,1)*rVariables.ConstitutiveMatrix(4,5)+rVariables.F(0,0)*rVariables.ConstitutiveMatrix(5,5));

    FD[8].resize(9);
    FD[8][0]=(rVariables.F(2,0)*rVariables.ConstitutiveMatrix(0,0)+rVariables.F(2,1)*rVariables.ConstitutiveMatrix(3,0)+rVariables.F(2,2)*rVariables.ConstitutiveMatrix(5,0));
    FD[8][1]=(rVariables.F(2,0)*rVariables.ConstitutiveMatrix(0,1)+rVariables.F(2,1)*rVariables.ConstitutiveMatrix(3,1)+rVariables.F(2,2)*rVariables.ConstitutiveMatrix(5,1));
    FD[8][2]=(rVariables.F(2,0)*rVariables.ConstitutiveMatrix(0,2)+rVariables.F(2,1)*rVariables.ConstitutiveMatrix(3,2)+rVariables.F(2,2)*rVariables.ConstitutiveMatrix(5,2));
    FD[8][3]=(rVariables.F(2,0)*rVariables.ConstitutiveMatrix(0,3)+rVariables.F(2,1)*rVariables.ConstitutiveMatrix(3,3)+rVariables.F(2,2)*rVariables.ConstitutiveMatrix(5,3));
    FD[8][4]=(rVariables.F(2,0)*rVariables.ConstitutiveMatrix(0,3)+rVariables.F(2,1)*rVariables.ConstitutiveMatrix(3,3)+rVariables.F(2,2)*rVariables.ConstitutiveMatrix(5,3));
    FD[8][5]=(rVariables.F(2,0)*rVariables.ConstitutiveMatrix(0,4)+rVariables.F(2,1)*rVariables.ConstitutiveMatrix(3,4)+rVariables.F(2,2)*rVariables.ConstitutiveMatrix(5,4));
    FD[8][6]=(rVariables.F(2,0)*rVariables.ConstitutiveMatrix(0,4)+rVariables.F(2,1)*rVariables.ConstitutiveMatrix(3,4)+rVariables.F(2,2)*rVariables.ConstitutiveMatrix(5,4));
    FD[8][7]=(rVariables.F(2,0)*rVariables.ConstitutiveMatrix(0,5)+rVariables.F(2,1)*rVariables.ConstitutiveMatrix(3,5)+rVariables.F(2,2)*rVariables.ConstitutiveMatrix(5,5));
    FD[8][8]=(rVariables.F(2,0)*rVariables.ConstitutiveMatrix(0,5)+rVariables.F(2,1)*rVariables.ConstitutiveMatrix(3,5)+rVariables.F(2,2)*rVariables.ConstitutiveMatrix(5,5));

    // simpler
    // FD[0].resize(9);
    // FD[0][0]=(rVariables.F(0,0)*rVariables.ConstitutiveMatrix(0,0));
    // FD[0][1]=(rVariables.F(0,0)*rVariables.ConstitutiveMatrix(0,1));
    // FD[0][2]=(rVariables.F(0,0)*rVariables.ConstitutiveMatrix(0,2));
    // FD[0][3]=(rVariables.F(0,1)*rVariables.ConstitutiveMatrix(3,3));
    // FD[0][4]=(rVariables.F(0,1)*rVariables.ConstitutiveMatrix(3,3));
    // FD[0][5]= 0.0;
    // FD[0][6]= 0.0;
    // FD[0][7]=(rVariables.F(0,2)*rVariables.ConstitutiveMatrix(5,5));
    // FD[0][8]=(rVariables.F(0,2)*rVariables.ConstitutiveMatrix(5,5));

    // FD[1].resize(9);
    // FD[1][0]=(rVariables.F(1,1)*rVariables.ConstitutiveMatrix(1,0));
    // FD[1][1]=(rVariables.F(1,1)*rVariables.ConstitutiveMatrix(1,1));
    // FD[1][2]=(rVariables.F(1,1)*rVariables.ConstitutiveMatrix(1,2));
    // FD[1][3]=(rVariables.F(1,0)*rVariables.ConstitutiveMatrix(3,3));
    // FD[1][4]=(rVariables.F(1,0)*rVariables.ConstitutiveMatrix(3,3));
    // FD[1][5]=(rVariables.F(1,2)*rVariables.ConstitutiveMatrix(4,4));
    // FD[1][6]=(rVariables.F(1,2)*rVariables.ConstitutiveMatrix(4,4));
    // FD[1][7]= 0.0;
    // FD[1][8]= 0.0;

    // FD[2].resize(9);
    // FD[2][0]=(rVariables.F(2,2)*rVariables.ConstitutiveMatrix(2,0));
    // FD[2][1]=(rVariables.F(2,2)*rVariables.ConstitutiveMatrix(2,1));
    // FD[2][2]=(rVariables.F(2,2)*rVariables.ConstitutiveMatrix(2,2));
    // FD[2][3]= 0.0;
    // FD[2][4]= 0.0;
    // FD[2][5]=(rVariables.F(2,1)*rVariables.ConstitutiveMatrix(4,4));
    // FD[2][6]=(rVariables.F(2,1)*rVariables.ConstitutiveMatrix(4,4));
    // FD[2][7]=(rVariables.F(2,0)*rVariables.ConstitutiveMatrix(5,5));
    // FD[2][8]=(rVariables.F(2,0)*rVariables.ConstitutiveMatrix(5,5));

    // FD[3].resize(9);
    // FD[3][0]=(rVariables.F(0,1)*rVariables.ConstitutiveMatrix(1,0));
    // FD[3][1]=(rVariables.F(0,1)*rVariables.ConstitutiveMatrix(1,1));
    // FD[3][2]=(rVariables.F(0,1)*rVariables.ConstitutiveMatrix(1,2));
    // FD[3][3]=(rVariables.F(0,0)*rVariables.ConstitutiveMatrix(3,3));
    // FD[3][4]=(rVariables.F(0,0)*rVariables.ConstitutiveMatrix(3,3));
    // FD[3][5]=(rVariables.F(0,2)*rVariables.ConstitutiveMatrix(4,4));
    // FD[3][6]=(rVariables.F(0,2)*rVariables.ConstitutiveMatrix(4,4));
    // FD[3][7]= 0.0;
    // FD[3][8]= 0.0;

    // FD[4].resize(9);
    // FD[4][0]=(rVariables.F(1,0)*rVariables.ConstitutiveMatrix(0,0));
    // FD[4][1]=(rVariables.F(1,0)*rVariables.ConstitutiveMatrix(0,1));
    // FD[4][2]=(rVariables.F(1,0)*rVariables.ConstitutiveMatrix(0,2));
    // FD[4][3]=(rVariables.F(1,1)*rVariables.ConstitutiveMatrix(3,3));
    // FD[4][4]=(rVariables.F(1,1)*rVariables.ConstitutiveMatrix(3,3));
    // FD[4][5]= 0.0;
    // FD[4][6]= 0.0;
    // FD[4][7]=(rVariables.F(1,2)*rVariables.ConstitutiveMatrix(5,5));
    // FD[4][8]=(rVariables.F(1,2)*rVariables.ConstitutiveMatrix(5,5));

    // FD[5].resize(9);
    // FD[5][0]=(rVariables.F(1,2)*rVariables.ConstitutiveMatrix(2,0));
    // FD[5][1]=(rVariables.F(1,2)*rVariables.ConstitutiveMatrix(2,1));
    // FD[5][2]=(rVariables.F(1,2)*rVariables.ConstitutiveMatrix(2,2));
    // FD[5][3]= 0.0;
    // FD[5][4]= 0.0;
    // FD[5][5]=(rVariables.F(1,1)*rVariables.ConstitutiveMatrix(4,4));
    // FD[5][6]=(rVariables.F(1,1)*rVariables.ConstitutiveMatrix(4,4));
    // FD[5][7]=(rVariables.F(1,0)*rVariables.ConstitutiveMatrix(5,5));
    // FD[5][8]=(rVariables.F(1,0)*rVariables.ConstitutiveMatrix(5,5));

    // FD[6].resize(9);
    // FD[6][0]=(rVariables.F(2,1)*rVariables.ConstitutiveMatrix(1,0));
    // FD[6][1]=(rVariables.F(2,1)*rVariables.ConstitutiveMatrix(1,1));
    // FD[6][2]=(rVariables.F(2,1)*rVariables.ConstitutiveMatrix(1,2));
    // FD[6][3]=(rVariables.F(2,0)*rVariables.ConstitutiveMatrix(3,3));
    // FD[6][4]=(rVariables.F(2,0)*rVariables.ConstitutiveMatrix(3,3));
    // FD[6][5]=(rVariables.F(2,2)*rVariables.ConstitutiveMatrix(4,4));
    // FD[6][6]=(rVariables.F(2,2)*rVariables.ConstitutiveMatrix(4,4));
    // FD[6][7]= 0.0;
    // FD[6][8]= 0.0;

    // FD[7].resize(9);
    // FD[7][0]=(rVariables.F(0,2)*rVariables.ConstitutiveMatrix(2,0));
    // FD[7][1]=(rVariables.F(0,2)*rVariables.ConstitutiveMatrix(2,1));
    // FD[7][2]=(rVariables.F(0,2)*rVariables.ConstitutiveMatrix(2,2));
    // FD[7][3]= 0.0;
    // FD[7][4]= 0.0;
    // FD[7][5]=(rVariables.F(0,1)*rVariables.ConstitutiveMatrix(4,4));
    // FD[7][6]=(rVariables.F(0,1)*rVariables.ConstitutiveMatrix(4,4));
    // FD[7][7]=(rVariables.F(0,0)*rVariables.ConstitutiveMatrix(5,5));
    // FD[7][8]=(rVariables.F(0,0)*rVariables.ConstitutiveMatrix(5,5));

    // FD[8].resize(9);
    // FD[8][0]=(rVariables.F(2,0)*rVariables.ConstitutiveMatrix(0,0));
    // FD[8][1]=(rVariables.F(2,0)*rVariables.ConstitutiveMatrix(0,1));
    // FD[8][2]=(rVariables.F(2,0)*rVariables.ConstitutiveMatrix(0,2));
    // FD[8][3]=(rVariables.F(2,1)*rVariables.ConstitutiveMatrix(3,3));
    // FD[8][4]=(rVariables.F(2,1)*rVariables.ConstitutiveMatrix(3,3));
    // FD[8][5]= 0.0;
    // FD[8][6]= 0.0;
    // FD[8][7]=(rVariables.F(2,2)*rVariables.ConstitutiveMatrix(5,5));
    // FD[8][8]=(rVariables.F(2,2)*rVariables.ConstitutiveMatrix(5,5));


    std::vector<Vector> FDB(9);

    FDB[0].resize(3);   FDB[1].resize(3);   FDB[2].resize(3);
    FDB[3].resize(3);   FDB[4].resize(3);   FDB[5].resize(3);
    FDB[6].resize(3);   FDB[7].resize(3);   FDB[8].resize(3);

    for(unsigned int i=0; i<9; i++)
      {
	noalias(FDB[i]) = ZeroVector(3);
      }

    for(unsigned int i=0; i<9; i++)
      {
        for(unsigned int j=0; j<3; j++)
	  {

	    FDB[i][j]= FD[i][0]*rVariables.F(j,0)*rVariables.DN_DX(ndi,0)+
	      FD[i][1]*rVariables.F(j,1)*rVariables.DN_DX(ndi,1)+
	      FD[i][2]*rVariables.F(j,2)*rVariables.DN_DX(ndi,2)+
	      (FD[i][3]+FD[i][4])*(0.5)*(rVariables.F(j,0)*rVariables.DN_DX(ndi,1)+rVariables.F(j,1)*rVariables.DN_DX(ndi,0))+
	      (FD[i][5]+FD[i][6])*(0.5)*(rVariables.F(j,1)*rVariables.DN_DX(ndi,2)+rVariables.F(j,2)*rVariables.DN_DX(ndi,1))+
	      (FD[i][7]+FD[i][8])*(0.5)*(rVariables.F(j,2)*rVariables.DN_DX(ndi,0)+rVariables.F(j,0)*rVariables.DN_DX(ndi,2));
	  }

      }

    //std::cout<<" FBD [0] "<<FDB[0]<<std::endl;
    //std::cout<<" FBD [1] "<<FDB[1]<<std::endl;
    //std::cout<<" FBD [2] "<<FDB[2]<<std::endl;
    //std::cout<<" FBD [3] "<<FDB[3]<<std::endl;


    for(unsigned int i=0; i<3; i++)
      {
    	rSigmaP[ndj][i] += rDirVector[0]*Normal[0]*(FDB[0][i])+
    	  rDirVector[1]*Normal[1]*(FDB[1][i])+
    	  rDirVector[2]*Normal[2]*(FDB[2][i])+

    	  rDirVector[0]*Normal[1]*(FDB[3][i])+
    	  rDirVector[1]*Normal[0]*(FDB[4][i])+
    	  rDirVector[1]*Normal[2]*(FDB[5][i])+

    	  rDirVector[2]*Normal[1]*(FDB[6][i])+
    	  rDirVector[0]*Normal[2]*(FDB[7][i])+
    	  rDirVector[2]*Normal[0]*(FDB[8][i]);

      }

  }


  //***********************************************************************************
  //************************************************************************************

  double& ContactDomainLM3DCondition::CalculateIntegrationWeight(double& rIntegrationWeight)
  {
    return rIntegrationWeight;
  }


  //************************************************************************************
  //************************************************************************************

  void ContactDomainLM3DCondition::CalculateNormalForce (double &F,ConditionVariables& rVariables,unsigned int& ndi,unsigned int& idir)
  {
    KRATOS_TRY

    F  = rVariables.Contact.Multiplier.Normal*rVariables.Contact.dN_dn[ndi]*rVariables.Contact.CurrentSurface.Normal[idir];
    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************


  void ContactDomainLM3DCondition::CalculateTangentStickForce (double &F,ConditionVariables& rVariables,unsigned int& ndi,unsigned int& idir)
  {
    KRATOS_TRY

    if( rVariables.Contact.Options.Is(ContactDomainUtilities::COMPUTE_FRICTION_FORCES) )
      {
	F = rVariables.Contact.CurrentGap.Normal * rVariables.Contact.CurrentSurface.Normal[idir] +
	  rVariables.Contact.Tangent.A.CurrentGap.Contravariant * rVariables.Contact.Tangent.CovariantBase.DirectionA[idir] +
	  rVariables.Contact.Tangent.B.CurrentGap.Contravariant * rVariables.Contact.Tangent.CovariantBase.DirectionB[idir];


	F *= ( rVariables.Contact.Tangent.A.Multiplier * rVariables.Contact.Tangent.A.dN_dt[ndi]+
	       rVariables.Contact.Tangent.B.Multiplier * rVariables.Contact.Tangent.B.dN_dt[ndi] );

	F += rVariables.Contact.dN_drn[ndi] * ( rVariables.Contact.Tangent.A.Multiplier * rVariables.Contact.Tangent.ContravariantBase.DirectionA[idir] + rVariables.Contact.Tangent.B.Multiplier * rVariables.Contact.Tangent.ContravariantBase.DirectionB[idir] );

      }
    else{
      F=0.0;
    }

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************


  void ContactDomainLM3DCondition::CalculateTangentSlipForce (double &F,ConditionVariables& rVariables,unsigned int& ndi,unsigned int& idir)
  {
    KRATOS_TRY

    if( rVariables.Contact.Options.Is(ContactDomainUtilities::COMPUTE_FRICTION_FORCES) )
      {
	F = (rVariables.Contact.Tangent.A.GapSign * rVariables.Contact.Tangent.CovariantBase.DirectionA[ndi] +
	     rVariables.Contact.Tangent.B.GapSign * rVariables.Contact.Tangent.CovariantBase.DirectionB[ndi] );

	F *= rVariables.Contact.CurrentGap.Normal * rVariables.Contact.CurrentSurface.Normal[idir];

	F += rVariables.Contact.dN_drn[ndi] * ( rVariables.Contact.Tangent.ContravariantBase.DirectionA[idir] +
						rVariables.Contact.Tangent.ContravariantBase.DirectionB[idir] );

	F *= ( rVariables.Contact.Multiplier.Normal * rVariables.Contact.FrictionCoefficient );
      }
    else{
      F=0.0;
    }

    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************

  void ContactDomainLM3DCondition::CalculateContactStiffness (double &Kcont,ConditionVariables& rVariables,unsigned int& ndi,unsigned int& ndj,unsigned int& idir,unsigned int& jdir)
  {
    KRATOS_TRY

    Kcont=0;

    double athird = 1.0/3.0;

    //Normal contact contribution:
    //KI:
    Kcont = (athird*rVariables.Contact.Tangent.EquivalentHeigh/(rVariables.Contact.ContactFactor.Normal)) *
        (rVariables.Contact.dN_dn[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dn[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]);

    Kcont-= (rVariables.Contact.Multiplier.Normal * rVariables.Contact.CurrentGap.Normal) *
        (rVariables.Contact.Tangent.ContravariantBase.Metric(0,0)*rVariables.Contact.Tangent.A.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.Tangent.A.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+
         rVariables.Contact.Tangent.ContravariantBase.Metric(0,1)*rVariables.Contact.Tangent.A.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.Tangent.B.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+
         rVariables.Contact.Tangent.ContravariantBase.Metric(1,0)*rVariables.Contact.Tangent.B.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.Tangent.A.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+
         rVariables.Contact.Tangent.ContravariantBase.Metric(1,1)*rVariables.Contact.Tangent.B.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.Tangent.B.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]);

    Kcont-= (rVariables.Contact.Multiplier.Normal) *
        (rVariables.Contact.Tangent.A.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dn[ndj]*rVariables.Contact.Tangent.ContravariantBase.DirectionA[jdir]+
         rVariables.Contact.Tangent.B.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dn[ndj]*rVariables.Contact.Tangent.ContravariantBase.DirectionB[jdir]+
         rVariables.Contact.dN_dn[ndi]*rVariables.Contact.Tangent.ContravariantBase.DirectionA[idir]*rVariables.Contact.Tangent.A.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+
         rVariables.Contact.dN_dn[ndi]*rVariables.Contact.Tangent.ContravariantBase.DirectionB[idir]*rVariables.Contact.Tangent.B.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]);

    Kcont-= (rVariables.Contact.dN_dn[ndi] * rVariables.Contact.CurrentSurface.Normal[idir]) *
        (rVariables.Contact.Tangent.A.CurrentTensil.Contravariant*rVariables.Contact.Tangent.A.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+
         rVariables.Contact.Tangent.B.CurrentTensil.Contravariant*rVariables.Contact.Tangent.B.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]);

    //covariant instead of contravariant in fortran implementation
    // Kcont = (athird*rVariables.Contact.Tangent.EquivalentHeigh/(rVariables.Contact.ContactFactor.Normal)) *
    //     (rVariables.Contact.dN_dn[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dn[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir] );
    // Kcont-= (rVariables.Contact.Multiplier.Normal * rVariables.Contact.CurrentGap.Normal ) *
    //    (rVariables.Contact.Tangent.CovariantBase.Metric(0,0)*rVariables.Contact.Tangent.A.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.Tangent.A.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+
    //     rVariables.Contact.Tangent.CovariantBase.Metric(0,1)*rVariables.Contact.Tangent.A.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.Tangent.B.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+
    //     rVariables.Contact.Tangent.CovariantBase.Metric(1,0)*rVariables.Contact.Tangent.B.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.Tangent.A.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+
    //     rVariables.Contact.Tangent.CovariantBase.Metric(1,1)*rVariables.Contact.Tangent.B.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.Tangent.B.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]);
    // Kcont -= (rVariables.Contact.Multiplier.Normal) *
    //   (rVariables.Contact.Tangent.A.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dn[ndj]*rVariables.Contact.Tangent.CovariantBase.DirectionA[jdir]+
    //    rVariables.Contact.Tangent.B.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dn[ndj]*rVariables.Contact.Tangent.CovariantBase.DirectionB[jdir]+
    //    rVariables.Contact.dN_dn[ndi]*rVariables.Contact.Tangent.CovariantBase.DirectionA[idir]*rVariables.Contact.Tangent.A.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+
    //    rVariables.Contact.dN_dn[ndi]*rVariables.Contact.Tangent.CovariantBase.DirectionB[idir]*rVariables.Contact.Tangent.B.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]);
    // Kcont -= (rVariables.Contact.dN_dn[ndi] * rVariables.Contact.CurrentSurface.Normal[idir]) *
    //   (rVariables.Contact.Tangent.A.CurrentTensil.Covariant*rVariables.Contact.Tangent.A.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+
    //    rVariables.Contact.Tangent.B.CurrentTensil.Covariant*rVariables.Contact.Tangent.B.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]);


    // rstiff = cte*geofct_n/tau_n*hdNn_a(in)*ne_a(idime)*hdNn_a(jn)*ne_a(jdime)

    // std::cout<<" ndi "<<ndi<<" ndj "<<ndj<<" i "<<idir<<" j "<<jdir<<std::endl;
    // std::cout<<" Kg "<<Kcont;
    // //std::cout<<" constant "<<constant<<" Kg "<<Kcont*constant;
    // double K1=Kcont;
    //std::cout<<" Ka "<<Kcont;
    //KII:

    //current configuration factor (UL-SL)
    double factor = 1.0; //rVariables.Contact.Tangent.CurrentArea/rVariables.Contact.Tangent.ReferenceArea;

    Kcont+= factor*rVariables.Contact.dN_dn[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.Nsigma[ndj][jdir];
    //std::cout<<" Kb "<<Kcont;

    rVariables.Contact.Options.Set(ContactDomainUtilities::COMPUTE_FRICTION_STIFFNESS,false); //friction needs an special treatment --> correct linearization is needed.

    //Stick contact contribution:
    if(rVariables.Contact.Options.IsNot(SLIP))
      {

	//std::cout<<" + stick ";
	if(rVariables.Contact.Options.Is(ContactDomainUtilities::COMPUTE_FRICTION_STIFFNESS))
	  {
	    //std::cout<<"(mu_on)";
	    //KI:
	    double raux=rVariables.Contact.CurrentGap.Normal*rVariables.Contact.CurrentSurface.Normal[idir]+
	      rVariables.Contact.Tangent.A.CurrentGap.Contravariant*rVariables.Contact.Tangent.CovariantBase.DirectionA[idir]+
	      rVariables.Contact.Tangent.B.CurrentGap.Contravariant*rVariables.Contact.Tangent.CovariantBase.DirectionB[idir]; //it must be GapcnTA (contravariant gap)
	    //std::cout<<" raux "<<raux;
	    //KI:
	    Kcont+= (raux*rVariables.Contact.Tangent.A.dN_dt[ndi]+rVariables.Contact.dN_drn[ndi]*rVariables.Contact.Tangent.ContravariantBase.DirectionA[idir]) *
	      ( (athird*rVariables.Contact.Tangent.EquivalentHeigh/(rVariables.Contact.ContactFactor.Normal) )*rVariables.Contact.dN_drn[ndj]*rVariables.Contact.Tangent.CovariantBase.DirectionA[jdir]+
		rVariables.Contact.Multiplier.Normal*rVariables.Contact.CurrentSurface.Normal[jdir]*(rVariables.Contact.Tangent.CovariantBase.Metric(0,0)*rVariables.Contact.Tangent.A.dN_dt[ndj]+rVariables.Contact.Tangent.CovariantBase.Metric(0,1)*rVariables.Contact.Tangent.B.dN_dt[ndj]) );

            Kcont+= (raux*rVariables.Contact.Tangent.B.dN_dt[ndi]+rVariables.Contact.dN_drn[ndi]*rVariables.Contact.Tangent.ContravariantBase.DirectionB[idir])*
	      ( (athird*rVariables.Contact.Tangent.EquivalentHeigh/(rVariables.Contact.ContactFactor.Normal) )*rVariables.Contact.dN_drn[ndj]*rVariables.Contact.Tangent.CovariantBase.DirectionB[jdir]+
		rVariables.Contact.Multiplier.Normal*rVariables.Contact.CurrentSurface.Normal[jdir]*(rVariables.Contact.Tangent.CovariantBase.Metric(1,0)*rVariables.Contact.Tangent.A.dN_dt[ndj]+rVariables.Contact.Tangent.CovariantBase.Metric(1,1)*rVariables.Contact.Tangent.B.dN_dt[ndj]) );

            Kcont+= rVariables.Contact.Tangent.A.Multiplier*(rVariables.Contact.Tangent.A.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dn[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.Tangent.A.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+(rVariables.Contact.Tangent.ContravariantBase.DirectionA[idir]*rVariables.Contact.Tangent.CovariantBase.DirectionA[jdir]+rVariables.Contact.Tangent.ContravariantBase.DirectionB[idir]*rVariables.Contact.Tangent.CovariantBase.DirectionB[jdir])*rVariables.Contact.Tangent.A.dN_dt[ndi]*rVariables.Contact.dN_drn[ndj]+(rVariables.Contact.Tangent.A.CurrentGap.Contravariant*rVariables.Contact.Tangent.A.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+ rVariables.Contact.Tangent.B.CurrentGap.Contravariant*rVariables.Contact.Tangent.B.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir])*rVariables.Contact.Tangent.A.dN_dt(ndi)*rVariables.Contact.CurrentSurface.Normal(idir)-raux*rVariables.Contact.Tangent.A.dN_dt(ndi)*rVariables.Contact.Tangent.CovariantBase.DirectionA(jdir)*rVariables.Contact.Tangent.A.dN_dt(ndj)-raux*rVariables.Contact.Tangent.B.dN_dt(ndi)*rVariables.Contact.Tangent.CovariantBase.DirectionB(jdir)*rVariables.Contact.Tangent.A.dN_dt(ndj));

            Kcont-= rVariables.Contact.Tangent.B.Multiplier*(rVariables.Contact.Tangent.B.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dn[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.Tangent.B.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+(rVariables.Contact.Tangent.ContravariantBase.DirectionA[idir]*rVariables.Contact.Tangent.CovariantBase.DirectionA[jdir]+rVariables.Contact.Tangent.ContravariantBase.DirectionB[idir]*rVariables.Contact.Tangent.CovariantBase.DirectionB[jdir])*rVariables.Contact.Tangent.B.dN_dt[ndi]*rVariables.Contact.dN_drn[ndj]+(rVariables.Contact.Tangent.A.CurrentGap.Contravariant*rVariables.Contact.Tangent.A.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+ rVariables.Contact.Tangent.B.CurrentGap.Contravariant*rVariables.Contact.Tangent.B.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir])*rVariables.Contact.Tangent.B.dN_dt(ndi)*rVariables.Contact.CurrentSurface.Normal(idir)-raux*rVariables.Contact.Tangent.A.dN_dt(ndi)*rVariables.Contact.Tangent.CovariantBase.DirectionA(jdir)*rVariables.Contact.Tangent.B.dN_dt(ndj)-raux*rVariables.Contact.Tangent.B.dN_dt(ndi)*rVariables.Contact.Tangent.CovariantBase.DirectionB(jdir)*rVariables.Contact.Tangent.B.dN_dt(ndj));

	    // rstiff = rstiff + (raux*dNt1_a(in) + t1cn(idime)*hdNn_n(in))*
	    //                    ( cte*geofct_n/tau_t*t1cv(jdime)*hdNn_n(jn) +
	    //                      lam_n*ne_a(jdime)*(gmcov(1,1)*dNt1_a(jn) + gmcov(1,2)*dNt2_a(jn)) ) +

	    //                    (raux*dNt2_a(in) + t2cn(idime)*hdNn_n(in))*
	    //                    ( cte*geofct_n/tau_t*t2cv(jdime)*hdNn_n(jn) +
	    //                      lam_n*ne_a(jdime)*(gmcov(2,1)*dNt1_a(jn) + gmcov(2,2)*dNt2_a(jn)) ) +

	    //           lam_t1*( ne_a(idime)*dNt1_a(in)*ne_a(jdime)*hdNn_a(jn) +
	    // 			  ne_a(idime)*hdNn_n(in)*ne_a(jdime)*dNt1_a(jn) +  &

	    //                    ( t1cn(idime)*t1cv(jdime) + t2cn(idime)*t2cv(jdime) )*dNt1_a(in)*hdNn_n(jn) +
	    //                    ( gt1_cv*ne_a(jdime)*dNt1_a(jn) + gt2_cv*ne_a(jdime)*dNt2_a(jn) )*ne_a(idime)*dNt1_a(in) -
	    //                    raux*dNt1_a(in)*t1cv(jdime)*dNt1_a(jn) -
	    // 			  raux*dNt2_a(in)*t2cv(jdime)*dNt1_a(jn) ) +              &

	    //           lam_t2*( ne_a(idime)*dNt2_a(in)*ne_a(jdime)*hdNn_a(jn) + ne_a(idime)*hdNn_n(in)*ne_a(jdime)*dNt2_a(jn) +
	    //                    ( t1cn(idime)*t1cv(jdime) + t2cn(idime)*t2cv(jdime) )*dNt2_a(in)*hdNn_n(jn) +
	    //                    ( gt1_cv*ne_a(jdime)*dNt1_a(jn) + gt2_cv*ne_a(jdime)*dNt2_a(jn) )*ne_a(idime)*dNt2_a(in) -
	    //                    raux*dNt1_a(in)*t1cv(jdime)*dNt2_a(jn) -
	    // 			  raux*dNt2_a(in)*t2cv(jdime)*dNt2_a(jn) )

	    //std::cout<<" Kc "<<Kcont;

            //KII:
	    Kcont+= ((rVariables.Contact.CurrentGap.Normal*rVariables.Contact.CurrentSurface.Normal[idir]+rVariables.Contact.Tangent.A.CurrentGap.Contravariant*rVariables.Contact.Tangent.CovariantBase.DirectionA[idir]+rVariables.Contact.Tangent.B.CurrentGap.Contravariant*rVariables.Contact.Tangent.CovariantBase.DirectionB[idir])*(rVariables.Contact.Tangent.A.dN_dt[ndi]*rVariables.Contact.Tangent.A.Tsigma[ndj][jdir]+rVariables.Contact.Tangent.B.dN_dt[ndi]*rVariables.Contact.Tangent.B.Tsigma[ndj][jdir])+(rVariables.Contact.dN_drn[ndi]*(rVariables.Contact.Tangent.ContravariantBase.DirectionA[idir]*rVariables.Contact.Tangent.A.Tsigma[ndj][jdir]+rVariables.Contact.Tangent.ContravariantBase.DirectionB[idir]*rVariables.Contact.Tangent.B.Tsigma[ndj][jdir])));

	    //std::cout<<" Kd "<<Kcont<<std::endl;

	    // rstiff = rstiff + ( he_a*ne_a(idime) + gt1_cn*t1cv(idime) + gt2_cn*t2cv(idime) )*
	    //                               ( dNt1_a(in)*raux1 + dNt2_a(in)*raux2 ) +                        &
	    //                               hdNn_n(in)*( t1cn(idime)*raux1 + t2cn(idime)*raux2 )
	  }
      }
    else
      {
        //Slip contact contribution:
	//std::cout<<" + slip ";
        if(rVariables.Contact.Options.Is(ContactDomainUtilities::COMPUTE_FRICTION_STIFFNESS))
	  {


	    //          raux = he_a*ne_a(idime) + gt1_cn*t1cv(idime) + gt2_cn*t2cv(idime)

	    //std::cout<<"(mu_on)";
            //KI:
	    //     rstiff = rstiff + (raux*dNt1_a(in) + t1cn(idime)*hdNn_n(in))*
	    //                                   ( cte*geofct_n/tau_t*t1cv(jdime)*hdNn_n(jn) +
	    //                                     lam_n*ne_a(jdime)*(gmcov(1,1)*dNt1_a(jn) + gmcov(1,2)*dNt2_a(jn)) ) +
	    //                       (raux*dNt2_a(in) + t2cn(idime)*hdNn_n(in))*
	    //   	                             ( cte*geofct_n/tau_t*t2cv(jdime)*hdNn_n(jn) +
	    //                                     lam_n*ne_a(jdime)*(gmcov(2,1)*dNt1_a(jn) + gmcov(2,2)*dNt2_a(jn)) ) +

	    //                          lam_t1*( ( gt1_cv*ne_a(jdime)*dNt1_a(jn) + gt2_cv*ne_a(jdime)*dNt2_a(jn) )*ne_a(idime)*dNt1_a(in) -
	    //                                   raux*dNt1_a(in)*t1cv(jdime)*dNt1_a(jn) - raux*dNt2_a(in)*t2cv(jdime)*dNt1_a(jn) -
	    //                                   t1cn(idime)*hdNn_n(in)*t1cv(jdime)*dNt1_a(jn) - t2cn(idime)*hdNn_n(in)*t2cv(jdime)*dNt1_a(jn) )+

	    //                          lam_t2*( ( gt1_cv*ne_a(jdime)*dNt1_a(jn) + gt2_cv*ne_a(jdime)*dNt2_a(jn) )*ne_a(idime)*dNt2_a(in) -
	    //                                   raux*dNt1_a(in)*t1cv(jdime)*dNt2_a(jn) - raux*dNt2_a(in)*t2cv(jdime)*dNt2_a(jn) -
	    //                                   t1cn(idime)*hdNn_n(in)*t1cv(jdime)*dNt2_a(jn) - t2cn(idime)*hdNn_n(in)*t2cv(jdime)*dNt2_a(jn)

            //KII:


	  }
      }


    // std::cout<<" ndi "<<ndi<<" ndj "<<ndj<<" i "<<idir<<" j "<<jdir<<std::endl;
    // std::cout<<" Kcont "<<Kcont;
    // std::cout<<" constant "<<athird * rVariables.Contact.Tangent.CurrentArea<<std::endl;

    //std::cout<<" Ks "<<Kcont-K1<<std::endl;


    KRATOS_CATCH(" ")
  }


  //************************************************************************************
  //************************************************************************************


  ContactDomainUtilities::PointType & ContactDomainLM3DCondition::CalculateCurrentTangent ( PointType &rTangent )
  {
    KRATOS_TRY

    unsigned int node1=mContactVariables.nodes[0];
    unsigned int node2=mContactVariables.nodes[1];

    PointType P1  =  GetGeometry()[node1].Coordinates();
    PointType P2  =  GetGeometry()[node2].Coordinates();

    //Set Reference Tangent
    rTangent=mContactUtilities.CalculateFaceTangent(rTangent,P1,P2);

    return rTangent;

    KRATOS_CATCH(" ")
  }

  //************************************************************************************
  //************************************************************************************


  inline bool ContactDomainLM3DCondition::CheckFictiousContacts(ConditionVariables& rVariables)
  {
    KRATOS_TRY

    bool real_contact = false;

    //Contact face segment node1-node2
    unsigned int node1=mContactVariables.nodes[0];
    unsigned int node2=mContactVariables.nodes[1];
    unsigned int slave=mContactVariables.slaves.front();

    double offset_factor = rVariables.Contact.CurrentGap.Normal;

    PointType PS  =  GetGeometry()[slave].Coordinates();

    PointType Normal = GetGeometry()[slave].FastGetSolutionStepValue(NORMAL);
    double  Shrink   = 1;//GetGeometry()[slave].FastGetSolutionStepValue(SHRINK_FACTOR);
    PointType Offset = GetGeometry()[slave].FastGetSolutionStepValue(OFFSET);
    offset_factor = norm_2(Offset);

    //modify slave position projection following slave normal
    double Sx1 = PS[0]+Normal[0]*Shrink*offset_factor;
    double Sy1 = PS[1]+Normal[1]*Shrink*offset_factor;

    //modify slave position projection following master normal
    double Mx1 = PS[0]-rVariables.Contact.CurrentSurface.Normal[0]*Shrink*offset_factor;
    double My1 = PS[1]-rVariables.Contact.CurrentSurface.Normal[1]*Shrink*offset_factor;

    //Domain neighbours:


    //Check slave node inside the contacting domain:

    //node1:
    ElementWeakPtrVectorType& nElements1 = GetGeometry()[node1].GetValue(NEIGHBOUR_ELEMENTS);
    //node2:
    ElementWeakPtrVectorType& nElements2 = GetGeometry()[node2].GetValue(NEIGHBOUR_ELEMENTS);

    bool is_inside_a = false;
    //following slave normal projection of the slave Sx1 and Sy1
    for(auto& i_nelem : nElements1)
    {
      GeometryType::PointsArrayType& vertices=i_nelem.GetGeometry().Points();

      is_inside_a = mContactUtilities.CalculatePosition( vertices[0].X(), vertices[0].Y(),
                                                         vertices[1].X(), vertices[1].Y(),
                                                         vertices[2].X(), vertices[2].Y(),
                                                         Sx1, Sy1);
      if(is_inside_a)
        break;
    }

    if(!is_inside_a){

      for(auto i_nelem : nElements2)
      {
        GeometryType::PointsArrayType& vertices=i_nelem.GetGeometry().Points();

        is_inside_a = mContactUtilities.CalculatePosition( vertices[0].X(), vertices[0].Y(),
                                                           vertices[1].X(), vertices[1].Y(),
                                                           vertices[2].X(), vertices[2].Y(),
                                                           Sx1, Sy1);

        if(is_inside_a)
          break;
      }

    }

    bool is_inside_b = false;
    //Check projection of the slave node inside the contacting domain:
    //following master normal projection of the slave Mx1 and My1
    for(auto& i_nelem : nElements1)
    {
      GeometryType::PointsArrayType& vertices=i_nelem.GetGeometry().Points();

      is_inside_b = mContactUtilities.CalculatePosition( vertices[0].X(), vertices[0].Y(),
                                                         vertices[1].X(), vertices[1].Y(),
                                                         vertices[2].X(), vertices[2].Y(),
                                                         Mx1, My1);

      if(is_inside_b)
        break;
    }

    if(!is_inside_b){


      for(auto i_nelem : nElements2)
      {
        GeometryType::PointsArrayType& vertices=i_nelem.GetGeometry().Points();

        is_inside_b = mContactUtilities.CalculatePosition( vertices[0].X(), vertices[0].Y(),
                                                           vertices[1].X(), vertices[1].Y(),
                                                           vertices[2].X(), vertices[2].Y(),
                                                           Mx1, My1);

        if(is_inside_b)
          break;
      }

    }

    if(is_inside_a && is_inside_b)
      real_contact = true; //if the slave node is inside of the domain --> real contact
    else
      real_contact = false;

    if( real_contact == false && (is_inside_b || is_inside_a) ){ //following the master normal is in.
      std::cout<<" THERE IS a SERIOUS DOUBT IN A FICTIOUS CONTACT "<<this->Id()<<std::endl;

      PointType P1  =  GetGeometry()[node1].Coordinates();
      PointType P2  =  GetGeometry()[node2].Coordinates();

      bool is_obtuse = mContactUtilities.CalculateObtuseAngle( P1[0], P1[1],
							       P2[0], P2[1],
							       PS[0], PS[1] );

      if(!is_obtuse){
	real_contact=true;
      }
      else{
	std::cout<<" BUT IT IS OBTUSE --> FICTIOUS "<<std::endl;
      }

    }

    double projection = inner_prod(Normal,rVariables.Contact.CurrentSurface.Normal);

    if(real_contact==false && fabs(projection)>0.707){
      real_contact = true;
      std::cout<<" NORMALS say that this is a REAL CONTACT "<<std::endl;
    }

    if(real_contact==false){
      std::cout<<" S normal ("<<Sx1<<","<<Sy1<<")"<<std::endl;
      std::cout<<" P normal ("<<Mx1<<","<<My1<<")"<<std::endl;
      std::cout<<" Current Normal "<<rVariables.Contact.CurrentSurface.Normal<<" Slave normal "<<Normal<<std::endl;
      std::cout<<" Shrink "<<Shrink<<" offset_factor "<<offset_factor<<" Shrink*offset_factor "<<Shrink*offset_factor<<std::endl;
    }

    return real_contact;

    KRATOS_CATCH(" ")
  }


  void ContactDomainLM3DCondition::save( Serializer& rSerializer ) const
  {
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ContactDomainCondition )
  }

  void ContactDomainLM3DCondition::load( Serializer& rSerializer )
  {
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ContactDomainCondition )
  }



} // Namespace Kratos
