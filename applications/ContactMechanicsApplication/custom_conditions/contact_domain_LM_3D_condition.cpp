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
#include "custom_conditions/contact_domain_LM_3D_condition.hpp"

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
    return Condition::Pointer(new ContactDomainLM3DCondition( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
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
    Element::ElementType& MasterElement = GetValue(MASTER_ELEMENTS).back();
    mContactVariables.SetMasterElement(MasterElement);

    Element::NodeType&    MasterNode  = GetValue(MASTER_NODES).back();
    mContactVariables.SetMasterNode(MasterNode);

    int  slave = -1;
    
    for(unsigned int i=0; i<MasterElement.GetGeometry().PointsNumber(); i++)
    {
        if(MasterNode.Id()==MasterElement.GetGeometry()[i].Id())
        {
	    slave=i;
	}
    }


    if(slave>=0)
    {

        NodesArrayType vertex;
	mContactVariables.order.resize(GetGeometry().PointsNumber(),false);
	
	int counter = 0;
        for(unsigned int i=0; i<GetGeometry().PointsNumber(); i++)
        {
            bool iset=false;
            for(unsigned int j=0; j<GetGeometry().PointsNumber(); j++)
            {

                if(GetGeometry()[i].Id()==MasterElement.GetGeometry()[j].Id())
                {
                    //vertex[i]=GetGeometry()[i];
		    //vertex.push_back( GetGeometry()[i] );
   		    //mContactVariables.nodes.push_back(i);
                    mContactVariables.order[i] = j;
 		    //std::cout<<" order ["<<i<<"] = "<<j<<std::endl;
                    iset=true;
		    break;
                }

            }

            if(iset==false)
            {
		//vertex[i]=GetGeometry()[i];
                //vertex[i]=MasterNode;
		//vertex.push_back( MasterNode );
         	mContactVariables.order[i] = slave;
		//std::cout<<" order ["<<i<<"] = "<<slave<<std::endl;
                mContactVariables.slaves.push_back(i);
            }

	    counter += mContactVariables.order[i];
        }


	if( counter != 6 )
	  KRATOS_THROW_ERROR( std::invalid_argument, "Something is Wrong with MASTERNODE and contact order", "" )

	//Permute
	std::vector<unsigned int> permute (7);
      
	permute[0]=0; 
	permute[1]=1;
	permute[2]=2;
	permute[3]=3;
	permute[4]=0;
	permute[5]=1;
	permute[6]=2;

	//reorder counter-clock-wise
	if(mContactVariables.slaves.size() == 1){

	  for( unsigned int i=1; i<3; i++ )
	    {
	      mContactVariables.nodes.push_back(permute[mContactVariables.slaves.back()+i]);
	    }
	  
	}
	else if(mContactVariables.slaves.size() == 2){

	  if( mContactVariables.slaves.back() == permute[mContactVariables.slaves.front()+1] ){
	    mContactVariables.nodes.push_back(permute[mContactVariables.slaves.back()+1]);
	    mContactVariables.nodes.push_back(permute[mContactVariables.slaves.back()+2]);
	  }
	  else{

	    if( mContactVariables.slaves.back() == permute[mContactVariables.slaves.front()+2] ){
	      mContactVariables.nodes.push_back(permute[mContactVariables.slaves.back()+1]);
	      mContactVariables.nodes.push_back(permute[mContactVariables.slaves.front()+1]);
	    }
	    else{

	      std::iter_swap(mContactVariables.slaves.begin(), mContactVariables.slaves.begin()+1);
		  
	      if( mContactVariables.slaves.back() == permute[mContactVariables.slaves.front()+1] ){
		mContactVariables.nodes.push_back(permute[mContactVariables.slaves.back()+1]);
		mContactVariables.nodes.push_back(permute[mContactVariables.slaves.back()+2]);
	      }
	      else{
		mContactVariables.nodes.push_back(permute[mContactVariables.slaves.back()+1]);
		mContactVariables.nodes.push_back(permute[mContactVariables.slaves.front()+1]);
	      } 
		  
	    }
	  }
	  
	}
	else{

	  std::cout<<" Problem in the SLAVES DETECTION "<<std::endl;
	}
    
	mContactVariables.SetMasterGeometry( MasterElement.GetGeometry() );
       
    }
    else
    {
        KRATOS_THROW_ERROR( std::invalid_argument, "MASTERNODE do not belongs to MASTER ELEMENT", "" )
    }

    
}

//************************************************************************************
//************************************************************************************


 
void ContactDomainLM3DCondition::CalculatePreviousGap() //prediction of the lagrange multiplier
{

  if( this->Is(SELECTED) )
    CalculatePreviousGapEdgeType(); 
  else
    CalculatePrefiousGapFaceType();

}
  

//************************************************************************************
//************************************************************************************

  
void ContactDomainLM3DCondition::TransformCovariantToContravariantBase(SurfaceBase& Covariant,SurfaceBase& Contravariant) //prediction of the lagrange multiplier face type element
{

  Covariant.Metric.resize(2,2);

  Covariant.Metric(0,0) = Covariant.DirectionA * Covariant.DirectionA;
  Covariant.Metric(0,1) = Covariant.DirectionA * Covariant.DirectionB;
  Covariant.Metric(1,0) = Covariant.DirectionB * Covariant.DirectionA;
  Covariant.Metric(1,1) = Covariant.DirectionB * Covariant.DirectionB;


  double MetricDet;
  MathUtils<double>::InvertMatrix2(Covariant.Metric,Contravariant.Metric, MetricDet);

  //transform DirectionA to the contravariant base
  Contravariant.DirectionA = Covariant.DirectionA * Contravariant.Metric(0,0) + Covariant.DirectionB * Contravariant.Metric(0,1);

  //transform DirectionB to the contravariant base
  Contravariant.DirectionB = Covariant.DirectionA * Contravariant.Metric(1,0) + Covariant.DirectionB * Contravariant.Metric(1,1);

}
  
//************************************************************************************
//************************************************************************************

  
void ContactDomainLM3DCondition::CalculatePreviousGapFaceType() //prediction of the lagrange multiplier face type element
{
 
    //Contact face node1-node2-node3
    unsigned int node1=mContactVariables.nodes[0];
    unsigned int node2=mContactVariables.nodes[1];
    unsigned int node3=mContactVariables.nodes[2];
    unsigned int slave=mContactVariables.slaves.back();

    //Get Reference Normal
    mContactVariables.ReferenceSurface.Normal=GetValue(NORMAL);

    //1.- Compute tension vector:  (must be updated each iteration)
    Matrix StressMatrix(3,3);
    noalias(StressMatrix) = ZeroMatrix(3,3);
    Matrix F(3,3);
    noalias(F) = ZeroMatrix(3,3);
    
    //a.- Assign initial 2nd Piola Kirchhoff stress:
    Condition::Pointer MasterCondition = GetValue(MASTER_CONDITION);

    //Get previous mechanics stored in the master node/condition
    Vector StressVector;
    StressVector = MasterCondition->GetValue(CAUCHY_STRESS_VECTOR);  //it means that has been stored
    F            = MasterCondition->GetValue(DEFORMATION_GRADIENT);  //it means that has been stored

    StressMatrix = MathUtils<double>::StressVectorToTensor( StressVector );

    //we are going to need F here from Cn-1 to Cn
    // F0 from C0 to Cn is need for the stress recovery on domain elements

    double detF =MathUtils<double>::Det(F);

    //b.- Compute the 1srt Piola Kirchhoff stress tensor  (P=J*CauchyStress*F^-T)
    mConstitutiveLawVector[0]->TransformStresses(StressMatrix,F,detF,ConstitutiveLaw::StressMeasure_Cauchy,ConstitutiveLaw::StressMeasure_PK1);

    //c.- Compute the tension (or traction) vector T=P*N (in the Reference configuration)   
    mContactVariables.TractionVector=prod(StressMatrix,mContactVariables.PreStepSurface.Normal);

    //d.- Compute (n-1) normals, tangents and relative displacements from historic mX on boundaries

    //Previous normal and tangent:  n_n-1,t_n-1

    //Previous Position    
    PointType P1  =  GetGeometry()[node1].Coordinates() - (GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1) - GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,2) );
    PointType P2  =  GetGeometry()[node2].Coordinates() - (GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1) - GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,2) );
    PointType P3  =  GetGeometry()[node3].Coordinates() - (GetGeometry()[node3].FastGetSolutionStepValue(DISPLACEMENT,1) - GetGeometry()[node3].FastGetSolutionStepValue(DISPLACEMENT,2) );
    PointType PS  =  GetGeometry()[slave].Coordinates() - (GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,1) - GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,2) );

    //Set Previous Tangent    
    mContactVariables.Tangent.CovariantBase.DirectionA = P2 - P1;
    mContactVariables.Tangent.CovariantBase.DirectionB = P3 - P1;

    TransformCovariantToContravariantBase(mContactVariables.Tangent.CovariantBase, mContactVariables.ContravariantBase);

    //Compute Previous Normal
    mContactVariables.PreStepSurface.Normal=mContactUtilities.CalculateSurfaceNormal(mContactVariables.PreStepSurface.Normal,mContactVariables.Tangent.CovariantBase.DirectionA,mContactVariables.Tangent.CovariantBase.DirectionB);


    if((inner_prod(mContactVariables.PreStepSurface.Normal,mContactVariables.ReferenceSurface.Normal))<0) //to give the correct direction
        mContactVariables.PreStepSurface.Normal*=-1;

    if(!(norm_2(mContactVariables.PreStepSurface.Normal)))
    {
        mContactVariables.PreStepSurface.Normal=mContactVariables.ReferenceSurface.Normal;
    }

  
    //e.- Compute A_n-1,B_n-1,L_n-1

    //A_n-1, B_n-1, L_n-1:
    std::vector<BaseLengths> PreviousBase(3);
    double EquivalentArea = 1;
    mContactUtilities.CalculateBaseDistances(PreviousBase,P1,P2,P3,PS,mContactVariables.PreStepSurface.Normal);
    mContactUtilities.CalculateBaseArea(EquivalentArea,PreviousBase[0].L,PreviousBase[1].L,PreviousBase[2].L);

    //complete the computation of the stabilization gap
    double ContactFactor = mContactVariables.StabilizationFactor * EquivalentArea;


    //f.-obtain the (g_N)3 and (g_T)3 for the n-1 configuration

    mContactVariables.PreviousGap.Normal  = 0;
    mContactVariables.PreviousGap.Normal = inner_prod((PS-P1),mContactVariables.PreStepSurface.Normal);

    
    mContactVariables.Tangent.PreviousGapA.Covariant = mContactVariables.PreviousGap.Normal * inner_prod(mContactVariables.Tangent.CovariantBase.DirectionA,mContactVariables.PreviousGap.Normal);
    mContactVariables.Tangent.PreviousGapB.Covariant = mContactVariables.PreviousGap.Normal * inner_prod(mContactVariables.Tangent.CovariantBase.DirectionB,mContactVariables.PreviousGap.Normal);
    
    mContactVariables.Tangent.PreviousGapA.Contravariant = mContactVariables.PreviousGap.Normal * inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionA,mContactVariables.PreviousGap.Normal);
    mContactVariables.Tangent.PreviousGapB.Contravariant = mContactVariables.PreviousGap.Normal * inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionB,mContactVariables.PreviousGap.Normal);


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
    NormalTensil=inner_prod(mContactVariables.PreStepSurface.Normal,mContactVariables.TractionVector);

    
    //i.- Compute tangent component of the tension vector:  (tt=cvt·P·N)
    TangentTensilcvA=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionA,mContactVariables.TractionVector);
    TangentTensilcvB=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionB,mContactVariables.TractionVector);
      
      
    //j.- Compute tangent component of the tension vector:  (tt=cnt·P·N)
    TangentTensilcnA=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionA,mContactVariables.TractionVector);
    TangentTensilcnB=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionB,mContactVariables.TractionVector);


    mContactVariables.PreviousGap.Normal  += 3 * ContactFactor * NormalTensil;

    mContactVariables.Tangent.PreviousGapA.Covariant += 3 * ContactFactor * TangentTensilcvA;
    mContactVariables.Tangent.PreviousGapB.Covariant += 3 * ContactFactor * TangentTensilcvB;

    mContactVariables.Tangent.PreviousGapA.Contravariant += 3 * ContactFactor * TangentTensilcnA;
    mContactVariables.Tangent.PreviousGapB.Contravariant += 3 * ContactFactor * TangentTensilcnB;


    //std::cout<<"ConditionID:  "<<this->Id()<<" -> Previous Tractions [tN:"<<NormalTensil<<"] "<<std::endl; 
    //std::cout<<" Previous Normal Gap [gN:"<<mContactVariables.PreviousGap.Normal<<"] "<<std::endl; 
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

    //Get Reference Normal
    mContactVariables.ReferenceSurface.Normal=GetValue(NORMAL);

    //1.- Compute tension vector:  (must be updated each iteration)
    Matrix StressMatrix(3,3);
    noalias(StressMatrix) = ZeroMatrix(3,3);
    Matrix F(3,3);
    noalias(F) = ZeroMatrix(3,3);
    
    //a.- Assign initial 2nd Piola Kirchhoff stress:
    Condition::Pointer MasterCondition = GetValue(MASTER_CONDITION);

    //Get previous mechanics stored in the master node/condition
    Vector StressVector;
    StressVector = MasterCondition->GetValue(CAUCHY_STRESS_VECTOR);  //it means that has been stored
    F            = MasterCondition->GetValue(DEFORMATION_GRADIENT);  //it means that has been stored

    StressMatrix = MathUtils<double>::StressVectorToTensor( StressVector );

    //we are going to need F here from Cn-1 to Cn
    // F0 from C0 to Cn is need for the stress recovery on domain elements

    double detF =MathUtils<double>::Det(F);

    //b.- Compute the 1srt Piola Kirchhoff stress tensor  (P=J*CauchyStress*F^-T)
    mConstitutiveLawVector[0]->TransformStresses(StressMatrix,F,detF,ConstitutiveLaw::StressMeasure_Cauchy,ConstitutiveLaw::StressMeasure_PK1);

    //c.- Compute the tension (or traction) vector T=P*N (in the Reference configuration)   
    mContactVariables.TractionVector=prod(StressMatrix,mContactVariables.PreStepSurface.Normal);

    //d.- Compute (n-1) normals, tangents and relative displacements from historic mX on boundaries

    //Previous normal and tangent:  n_n-1,t_n-1

    //Previous Position    
    PointType P1  =  GetGeometry()[node1].Coordinates() - (GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1) - GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,2) );
    PointType P2  =  GetGeometry()[node2].Coordinates() - (GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1) - GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,2) );
    PointType PS1  =  GetGeometry()[node3].Coordinates() - (GetGeometry()[slave1].FastGetSolutionStepValue(DISPLACEMENT,1) - GetGeometry()[slave1].FastGetSolutionStepValue(DISPLACEMENT,2) );

    PointType PS2  =  GetGeometry()[slave].Coordinates() - (GetGeometry()[slave2].FastGetSolutionStepValue(DISPLACEMENT,1) - GetGeometry()[slave2].FastGetSolutionStepValue(DISPLACEMENT,2) );

    //Set Previous Tangent    
    mContactVariables.Tangent.CovariantBase.DirectionA = P2 - P1;
    mContactVariables.Tangent.CovariantBase.DirectionB = PS2 - PS1;

    TransformCovariantToContravariantBase(mContactVariables.Tangent.CovariantBase, mContactVariables.ContravariantBase);

    //Take MasterCondition for computing reference normal:    
    PointType M1  =  MasterCondition.GetGeometry()[0].Coordinates() - (MasterCondition.GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT,1) - MasterCondition.GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT,2) );
    PointType M2  =  MasterCondition.GetGeometry()[1].Coordinates() - (MasterCondition.GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT,1) - MasterCondition.GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT,2) );
    PointType M3  =  MasterCondition.GetGeometry()[2].Coordinates() - (MasterCondition.GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT,1) - MasterCondition.GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT,2) );

    PointType V1 = M2-M1;
    PointType V2 = M3-M1;
    //Compute Previous Normal
    mContactVariables.PreStepSurface.Normal=mContactUtilities.CalculateSurfaceNormal(mContactVariables.PreStepSurface.Normal,V1,V2);

    //std::cout<<" Pre Normal ["<<this->Id()<<"] "<<mContactVariables.PreStepSurface.Normal<<std::endl;

    if((inner_prod(mContactVariables.PreStepSurface.Normal,mContactVariables.ReferenceSurface.Normal))<0) //to give the correct direction
        mContactVariables.PreStepSurface.Normal*=-1;

    if(!(norm_2(mContactVariables.PreStepSurface.Normal)))
    {
        mContactVariables.PreStepSurface.Normal=mContactVariables.ReferenceSurface.Normal;
    }


    //Reference normal: n_n,t_n  -> mContactVariables.ReferenceSurface.Normal / mContactVariables.ReferenceSurface.Tangent

    //e.- Compute A_n-1,B_n-1,L_n-1

    //A_n-1, B_n-1, L_n-1:
    std::vector<BaseLengths> PreviousBase(3);
    mContactUtilities.CalculateEdgeDistances(PreviousBase,P1,P2,PS1,PS2,mContactVariables.PreStepSurface.Normal);
    double EquivalentArea = 0.25 * (PreviousBase[0].L * PreviousBase[1].L) * (PreviousBase[0].L * PreviousBase[1].L);

    //complete the computation of the stabilization gap
    double ContactFactor = mContactVariables.StabilizationFactor * EquivalentArea;

    //f.-obtain the (g_N)3 and (g_T)3 for the n-1 configuration

    mContactVariables.PreviousGap.Normal  = 0;
    M1 = 0.5 * (P1+P2);
    M2 = 0.5 * (PS1+PS2);
    mContactVariables.PreviousGap.Normal = inner_prod((M2-M1),mContactVariables.PreStepSurface.Normal);

    
    mContactVariables.Tangent.PreviousGapA.Covariant = mContactVariables.PreviousGap.Normal*(mContactVariables.Tangent.CovariantBase.DirectionA*mContactVariables.PreviousGap.Normal);
    mContactVariables.Tangent.PreviousGapB.Covariant = mContactVariables.PreviousGap.Normal*(mContactVariables.Tangent.CovariantBase.DirectionB*mContactVariables.PreviousGap.Normal);
    
    mContactVariables.Tangent.PreviousGapA.Contravariant = mContactVariables.PreviousGap.Normal*(mContactVariables.Tangent.ContravariantBase.DirectionA*mContactVariables.PreviousGap.Normal);
    mContactVariables.Tangent.PreviousGapB.Contravariant = mContactVariables.PreviousGap.Normal*(mContactVariables.Tangent.ContravariantBase.DirectionB*mContactVariables.PreviousGap.Normal);

    PointType D1  =  GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1)-GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,2);
    PointType D2  =  GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1)-GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,2);
    PointType DS1  =  GetGeometry()[slave1].FastGetSolutionStepValue(DISPLACEMENT,1)-GetGeometry()[slave1].FastGetSolutionStepValue(DISPLACEMENT,2);
    PointType DS2  =  GetGeometry()[slave2].FastGetSolutionStepValue(DISPLACEMENT,1)-GetGeometry()[slave2].FastGetSolutionStepValue(DISPLACEMENT,2);

    
    //(g_N)3
    mContactVariables.PreviousGap.Normal*= inner_prod(mContactVariables.ReferenceSurface.Normal,mContactVariables.PreStepSurface.Normal);
    
    mContactVariables.PreviousGap.Normal+=inner_prod(mContactVariables.ReferenceSurface.Normal,(D1*(-PreviousBase[0].A/PreviousBase[0].L)));
    mContactVariables.PreviousGap.Normal+=inner_prod(mContactVariables.ReferenceSurface.Normal,(D2*(-PreviousBase[0].B/PreviousBase[0].L)));
    mContactVariables.PreviousGap.Normal+=inner_prod(mContactVariables.ReferenceSurface.Normal,(DS1*(PreviousBase[1].A/PreviousBase[1].L)));
    mContactVariables.PreviousGap.Normal+=inner_prod(mContactVariables.ReferenceSurface.Normal,(DS1*(PreviousBase[1].B/PreviousBase[1].L))));
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


    double EquivalentHeigh = norm_2(MathUtils<double>::CrossProduct(mContactVariables.Tangent.CovariantBase.DirectionA,mContactVariables.Tangent.CovariantBase.DirectionB) );

    EquivalentHeigh *= (0.5/PreviousBase[0].L);
  

    //d_n-1=X_n - X_n-1

    //g.- get total effective gap as: gap_n^eff=gap_n+(PreviousTimeStep/CurrentTimeStep)*(gap_n-gap_n-1)

    //gap_n-1 (in function of the n-1 position of hte other node) gap_n-1=(g_N)3_n-1+2*Tau*tn_n-1

    double NormalTensil=0,TangentTensilcvA=0,TangentTensilcvB=0,TangentTensilcnA=0,TangentTensilcnB=0;    
    
    //h.- Compute normal component of the tension vector:   (tn=n·P·N)
    NormalTensil=inner_prod(mContactVariables.PreStepSurface.Normal,mContactVariables.TractionVector);

    
    //i.- Compute tangent component of the tension vector:  (tt=cvt·P·N)
    TangentTensilcvA=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionA,mContactVariables.TractionVector);
    TangentTensilcvB=inner_prod(mContactVariables.Tangent.CovariantBase.DirectionB,mContactVariables.TractionVector);
      
      
    //j.- Compute tangent component of the tension vector:  (tt=cnt·P·N)
    TangentTensilcnA=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionA,mContactVariables.TractionVector);
    TangentTensilcnB=inner_prod(mContactVariables.Tangent.ContravariantBase.DirectionB,mContactVariables.TractionVector);


    mContactVariables.PreviousGap.Normal  += 3 * ContactFactor * NormalTensil / EquivalentHeigh;

    mContactVariables.Tangent.PreviousGapA.Covariant += 3 * ContactFactor * TangentTensilcvA / EquivalentHeigh;
    mContactVariables.Tangent.PreviousGapB.Covariant += 3 * ContactFactor * TangentTensilcvB / EquivalentHeigh;

    mContactVariables.Tangent.PreviousGapA.Contravariant += 3 * ContactFactor * TangentTensilcnA / EquivalentHeigh;
    mContactVariables.Tangent.PreviousGapB.Contravariant += 3 * ContactFactor * TangentTensilcnB / EquivalentHeigh;


    //std::cout<<"ConditionID:  "<<this->Id()<<" -> Previous Tractions [tN:"<<NormalTensil<<"] "<<std::endl; 
    //std::cout<<" Previous Normal Gap [gN:"<<mContactVariables.PreviousGap.Normal<<"] "<<std::endl; 
}


//**********************************COMPUTE TAU STAB**********************************
//************************************************************************************


void ContactDomainLM3DCondition::CalculateContactFactor( ProcessInfo& rCurrentProcessInfo )
{
    //Initilialize Tau for the stabilization
    double alpha_stab = 0.1;
    alpha_stab = GetProperties()[TAU_STAB];

    ElementType& MasterElement = mContactVariables.GetMasterElement();
  
    //Look at the nodes, get the slave and get the Emin

    //Contact face segment node1-node2
    unsigned int slave = mContactVariables.slaves.back();


    double Eslave = GetGeometry()[slave].GetValue(NEIGHBOUR_ELEMENTS)[0].GetProperties()[YOUNG_MODULUS];
    double Emin   = MasterElement.GetProperties()[YOUNG_MODULUS];

    //STANDARD OPTION
    if(Emin>Eslave)
	Emin=Eslave;

    mContactVariables.StabilizationFactor = alpha_stab/Emin;

    //EXPERIMENTAL OPTION
    // const GeometryType::IntegrationPointsArrayType& integration_points = MasterElement.GetGeometry().IntegrationPoints( mThisIntegrationMethod );
   
    // //Get Current ConstitutiveMatrix
    // int size = integration_points.size();
    // std::vector<Matrix> ConstitutiveMatrix(size);   
    // MasterElement.CalculateOnIntegrationPoints(CONSTITUTIVE_MATRIX,ConstitutiveMatrix,rCurrentProcessInfo);

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

void ContactDomainLM3DCondition::CalculateExplicitFactors(GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo)
{

  if( this->Is(SELECTED) )
    CalculateExplicitFactorsEdgeType(); 
  else
    CalculateExplicitFactorsFaceType();



  //set contact normal
  array_1d<double, 3> &ContactNormal  = GetGeometry()[slave].FastGetSolutionStepValue(CONTACT_NORMAL);
  
  for(unsigned int i=0; i<3; i++)
    ContactNormal[i] = rVariables.Contact.CurrentSurface.Normal[i];

    
  if(mContactVariables.IterationCounter < 1)
    mContactVariables.IterationCounter += 1;
  
}

//************************************************************************************
//************************************************************************************

void ContactDomainLM3DCondition::CalculateExplicitFactorsFaceType(GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo)
{

    //Contact face node1-node2-node3
    unsigned int node1=mContactVariables.nodes[0];
    unsigned int node2=mContactVariables.nodes[1];
    unsigned int node3=mContactVariables.nodes[2];
    unsigned int slave=mContactVariables.slaves.back();


    //1.- Compute tension vector:  (must be updated each iteration)
    Matrix StressMatrix(3,3);
    noalias(StressMatrix) = ZeroMatrix(3,3);
    
    //a.- Assign initial 2nd Piola Kirchhoff stress:
    StressMatrix=MathUtils<double>::StressVectorToTensor( rVariables.StressVector );

    // UL
    //b.- Compute the 1srt Piola Kirchhoff stress tensor
    StressMatrix = mConstitutiveLawVector[0]->TransformStresses(StressMatrix, rVariables.F, rVariables.detF, ConstitutiveLaw::StressMeasure_Cauchy, ConstitutiveLaw::StressMeasure_PK1);

    //c.- Compute the tension (or traction) vector T=P*N (in the Reference configuration)
    mContactVariables.TractionVector=prod(StressMatrix,mContactVariables.ReferenceSurface.Normal);


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

    //Compute Current Normal
    rVariables.Contact.CurrentSurface.Normal=mContactUtilities.CalculateSurfaceNormal(rVariables.Contact.CurrentSurface.Normal,rVariables.Contact.Tangent.CovariantBase.DirectionA,rVariables.Contact.Tangent.CovariantBase.DirectionB);


    if(double(inner_prod(rVariables.Contact.CurrentSurface.Normal,mContactVariables.ReferenceSurface.Normal))<0) //to give the correct direction
        rVariables.Contact.CurrentSurface.Normal*=-1;

    rVariables.Contact.CurrentSurface.Normal /= norm_2(rVariables.Contact.CurrentSurface.Normal);  //to be unitary

    if(!(norm_2(rVariables.Contact.CurrentSurface.Normal)))
        rVariables.Contact.CurrentSurface.Normal=mContactVariables.ReferenceSurface.Normal;



    //4.- Compute Effective Gaps: (g^eff=g_n3+2*Tau*tn=2*Tau*Multiplier.Normal)

    //Reference normal: n_n,t_n  -> mContactVariables.ReferenceSurface.Normal / rVariables.Contact.Tangent
    //Current normal:   n,t      -> rVariables.Contact.CurrentSurface.Normal /  rVariables.Contact.CurrentSurface.Tangent

    //e.- Compute A_n,B_n,L_n
    rVariables.Contact.ReferenceBase.resize(3);
    rVariables.Contact.CurrentBase.resize(3);

    //a, b, l:
    mContactUtilities.CalculateBaseDistances(rVariables.Contact.CurrentBase,P1,P2,P3,PS,rVariables.Contact.CurrentSurface.Normal);
    mContactUtilities.CalculateBaseArea(rVariables.Contact.Tangent.EquivalentArea,rVariables.Contact.CurrentBase[0].L,rVariables.Contact.CurrentBase[1].L,rVariables.Contact.CurrentBase[2].L);

    //A, B, L:

    //Reference Position
    PointType P1  =  GetGeometry()[node1].Coordinates() - (GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1) );
    PointType P2  =  GetGeometry()[node2].Coordinates() - (GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1) );
    PointType P3  =  GetGeometry()[node3].Coordinates() - (GetGeometry()[node3].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node3].FastGetSolutionStepValue(DISPLACEMENT,1) );
    PointType PS  =  GetGeometry()[slave].Coordinates() - (GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,1) );

    mContactUtilities.CalculateBaseDistances(rVariables.Contact.ReferenceBase,P1,P2,P3,PS,mContactVariables.ReferenceSurface.Normal);


    //complete the computation of the stabilization gap
    rVariables.Contact.ContactFactor.Normal  =  mContactVariables.StabilizationFactor * rVariables.Contact.Tangent.EquivalentArea;
    rVariables.Contact.ContactFactor.Tangent =  rVariables.Contact.ContactFactor.Normal;
    
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
    rVariables.Contact.CurrentTensil.Normal=inner_prod(mContactVariables.CurrentSurface.Normal,mContactVariables.TractionVector);

    
    //i.- Compute tangent component of the tension vector:  (tt=cvt·P·N)
    rVariables.Contact.Tangent.A.CurrentTensil.Covariant=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionA,mContactVariables.TractionVector);
    rVariables.Contact.Tangent.B.CurrentTensil.Covariant=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionB,mContactVariables.TractionVector);
      
      
    //j.- Compute tangent component of the tension vector:  (tt=cnt·P·N)
    rVariables.Contact.Tangent.A.CurrentTensil.Contravariant=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionA,mContactVariables.TractionVector);
    rVariables.Contact.Tangent.B.CurrentTensil.Contravariant=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionB,mContactVariables.TractionVector);

    
    ReferenceGapN+= 3 * rVariables.Contact.ContactFactor.Normal * rVariables.Contact.CurrentTensil.Normal;

    ReferenceGapcvTA+= 3 * rVariables.Contact.ContactFactor.Tangent * rVariables.Contact.Tangent.A.CurrentTensil.Covariant;
    ReferenceGapcvTB+= 3 * rVariables.Contact.ContactFactor.Tangent * rVariables.Contact.Tangent.B.CurrentTensil.Covariant;
    
    ReferenceGapcnTA+= 3 * rVariables.Contact.ContactFactor.Tangent * rVariables.Contact.Tangent.A.CurrentTensil.Contravariant;
    ReferenceGapcnTB+= 3 * rVariables.Contact.ContactFactor.Tangent * rVariables.Contact.Tangent.B.CurrentTensil.Contravariant;


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
    if(fabs(EffectiveGapN)<= 1e-15 && fabs(EffectiveGapN)<= 1e-15)
      EffectiveGapN = 0;
    //decimal correction from tension vector calculation


    if(EffectiveGapN<=0)   //if(EffectiveGap<0){
    {

        rVariables.Contact.Options.Set(ACTIVE);  //normal contact active

        if( (fabs(EffectiveGapcvTa+EffectiveGapcvTB)<=rVariables.Contact.FrictionCoefficient*fabs(EffectiveGapN)) ||
	    (fabs(EffectiveGapcnTa+EffectiveGapcnTB)<=rVariables.Contact.FrictionCoefficient*fabs(EffectiveGapN)) )
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
}


//************************************************************************************
//************************************************************************************

void ContactDomainLM3DCondition::CalculateExplicitFactorsEdgeType(GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo)
{

    //Contact face node1-node2-node3
    unsigned int node1=mContactVariables.nodes[0];
    unsigned int node2=mContactVariables.nodes[1];
    unsigned int slave1=mContactVariables.slave[0];
    unsigned int slave2=mContactVariables.slave[1];


    //1.- Compute tension vector:  (must be updated each iteration)
    Matrix StressMatrix(3,3);
    noalias(StressMatrix) = ZeroMatrix(3,3);
    
    //a.- Assign initial 2nd Piola Kirchhoff stress:
    StressMatrix=MathUtils<double>::StressVectorToTensor( rVariables.StressVector );

    // UL
    //b.- Compute the 1srt Piola Kirchhoff stress tensor
    StressMatrix = mConstitutiveLawVector[0]->TransformStresses(StressMatrix, rVariables.F, rVariables.detF, ConstitutiveLaw::StressMeasure_Cauchy, ConstitutiveLaw::StressMeasure_PK1);

    //c.- Compute the tension (or traction) vector T=P*N (in the Reference configuration)
    mContactVariables.TractionVector=prod(StressMatrix,mContactVariables.ReferenceSurface.Normal);


    //d.- Compute the Current Normal and Tangent

    //Current Position    
    PointType P1  =  GetGeometry()[node1].Coordinates();
    PointType P2  =  GetGeometry()[node2].Coordinates();
    PointType PS1 =  GetGeometry()[slave1].Coordinates();
    PointType PS2 =  GetGeometry()[slave2].Coordinates();

    //Set Current Tangent    
    rVariables.Contact.Tangent.CovariantBase.DirectionA = P2 - P1;
    rVariables.Contact.Tangent.CovariantBase.DirectionB = PS2 - PS1;

    TransformCovariantToContravariantBase(rVariables.Contact.Tangent.CovariantBase, rVariables.Contact.Tangent.ContravariantBase);


    //Take MasterCondition for computing current normal:    
    PointType M1  =  MasterCondition.GetGeometry()[0].Coordinates();
    PointType M2  =  MasterCondition.GetGeometry()[1].Coordinates();
    PointType M3  =  MasterCondition.GetGeometry()[2].Coordinates();

    PointType V1 = M2-M1;
    PointType V2 = M3-M1;
   
    //Compute Current Normal
    rVariables.Contact.CurrentSurface.Normal=mContactUtilities.CalculateSurfaceNormal(rVariables.Contact.CurrentSurface.Normal,V1,V2);


    if(double(inner_prod(rVariables.Contact.CurrentSurface.Normal,mContactVariables.ReferenceSurface.Normal))<0) //to give the correct direction
        rVariables.Contact.CurrentSurface.Normal*=-1;

    rVariables.Contact.CurrentSurface.Normal /= norm_2(rVariables.Contact.CurrentSurface.Normal);  //to be unitary

    if(!(norm_2(rVariables.Contact.CurrentSurface.Normal)))
        rVariables.Contact.CurrentSurface.Normal=mContactVariables.ReferenceSurface.Normal;



    //4.- Compute Effective Gaps: (g^eff=g_n3+2*Tau*tn=2*Tau*Multiplier.Normal)

    //Reference normal: n_n,t_n  -> mContactVariables.ReferenceSurface.Normal / rVariables.Contact.Tangent
    //Current normal:   n,t      -> rVariables.Contact.CurrentSurface.Normal /  rVariables.Contact.CurrentSurface.Tangent

    //e.- Compute A_n,B_n,L_n
    rVariables.Contact.ReferenceBase.resize(3);
    rVariables.Contact.CurrentBase.resize(3);

    //a, b, l:
    mContactUtilities.CalculateEdgeDistances(rVariables.Contact.CurrentBase,P1,P2,PS1,PS2,rVariables.Contact.CurrentSurface.Normal);
    rVariables.Contact.Tangent.EquivalentArea = 0.25 * (rVariables.Contact.CurrentBase[0].L * rVariables.Contact.CurrentBase[1].L) * (rVariables.Contact.CurrentBase[0].L * rVariables.Contact.CurrentBase[1].L);

    //A, B, L:

    //Reference Position
    PointType P1  =  GetGeometry()[node1].Coordinates() - (GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1) );
    PointType P2  =  GetGeometry()[node2].Coordinates() - (GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1) );
    PointType PS1  =  GetGeometry()[slave1].Coordinates() - (GetGeometry()[slave1].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[slave1].FastGetSolutionStepValue(DISPLACEMENT,1) );
    PointType PS2  =  GetGeometry()[slave2].Coordinates() - (GetGeometry()[slave2].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[slave2].FastGetSolutionStepValue(DISPLACEMENT,1) );

    mContactUtilities.CalculateEdgeDistances(rVariables.Contact.ReferenceBase,P1,P2,PS1,PS2,mContactVariables.ReferenceSurface.Normal);


    //complete the computation of the stabilization gap
    rVariables.Contact.ContactFactor.Normal  =  mContactVariables.StabilizationFactor * rVariables.Contact.Tangent.EquivalentArea;
    rVariables.Contact.ContactFactor.Tangent =  rVariables.Contact.ContactFactor.Normal;
    
    //f.-obtain the (g_N)3 and (g_T)3 for the n configuration

    M1 = 0.5 * (P1+P2);
    M2 = 0.5 * (PS1+PS2);
    double ReferenceGapN = inner_prod((M2-M1),mContactVariables.PreStepSurface.Normal);

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

    rVariables.Contact.Tangent.EquivalentHeigh = norm_2(MathUtils<double>::CrossProduct(rVariables.Contact.Tangent.CovariantBase.DirectionA,rVariables.Contact.Tangent.CovariantBase.DirectionB) );

    rVariables.Contact.Tangent.EquivalentHeigh *= (0.5/rVariables.Contact.ReferenceBase[0].L);
    
    //...
    
    rVariables.Contact.CurrentGap.Normal=ReferenceGapN; //(g_N)3 -- needed in the Kcont1 computation

    rVariables.Contact.Tangent.A.CurrentGap.Covariant=ReferenceGapcvTA; //(g_cvTA)3 -- needed in the Kcont1 computation
    rVariables.Contact.Tangent.A.CurrentGap.Contravariant=ReferenceGapcnTA; //(g_cnTA)3 -- needed in the Kcont1 computation
    
    rVariables.Contact.Tangent.B.CurrentGap.Covariant=ReferenceGapcvTB; //(g_cvTB)3 -- needed in the Kcont1 computation
    rVariables.Contact.Tangent.B.CurrentGap.Contravariant=ReferenceGapcnTB; //(g_nTB)3 -- needed in the Kcont1 computation

    //g.- get total effective gap as: gap_n^eff=gap_n+(PreviousTimeStep/CurrentTimeStep)*(gap_n-gap_n-1)

    //gap_n   (in function of the n position of the other node) gap_n=(g_N)3+2*Tau*tn_n

    
    //h.- Compute normal component of the tension vector:   (tn=n·P·N)
    rVariables.Contact.CurrentTensil.Normal=inner_prod(mContactVariables.CurrentSurface.Normal,mContactVariables.TractionVector);

    
    //i.- Compute tangent component of the tension vector:  (tt=cvt·P·N)
    rVariables.Contact.Tangent.A.CurrentTensil.Covariant=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionA,mContactVariables.TractionVector);
    rVariables.Contact.Tangent.B.CurrentTensil.Covariant=inner_prod(rVariables.Contact.Tangent.CovariantBase.DirectionB,mContactVariables.TractionVector);
      
      
    //j.- Compute tangent component of the tension vector:  (tt=cnt·P·N)
    rVariables.Contact.Tangent.A.CurrentTensil.Contravariant=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionA,mContactVariables.TractionVector);
    rVariables.Contact.Tangent.B.CurrentTensil.Contravariant=inner_prod(rVariables.Contact.Tangent.ContravariantBase.DirectionB,mContactVariables.TractionVector);

    
    ReferenceGapN+= 3 * rVariables.Contact.ContactFactor.Normal * rVariables.Contact.CurrentTensil.Normal;

    ReferenceGapcvTA+= 3 * rVariables.Contact.ContactFactor.Tangent * rVariables.Contact.Tangent.A.CurrentTensil.Covariant;
    ReferenceGapcvTB+= 3 * rVariables.Contact.ContactFactor.Tangent * rVariables.Contact.Tangent.B.CurrentTensil.Covariant;
    
    ReferenceGapcnTA+= 3 * rVariables.Contact.ContactFactor.Tangent * rVariables.Contact.Tangent.A.CurrentTensil.Contravariant;
    ReferenceGapcnTB+= 3 * rVariables.Contact.ContactFactor.Tangent * rVariables.Contact.Tangent.B.CurrentTensil.Contravariant;


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
    if(fabs(EffectiveGapN)<= 1e-15 && fabs(EffectiveGapN)<= 1e-15)
      EffectiveGapN = 0;
    //decimal correction from tension vector calculation


    if(EffectiveGapN<=0)   //if(EffectiveGap<0){
    {

        rVariables.Contact.Options.Set(ACTIVE);  //normal contact active

        if( (fabs(EffectiveGapcvTa+EffectiveGapcvTB)<=rVariables.Contact.FrictionCoefficient*fabs(EffectiveGapN)) ||
	    (fabs(EffectiveGapcnTa+EffectiveGapcnTB)<=rVariables.Contact.FrictionCoefficient*fabs(EffectiveGapN)) )
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
    
}



//************************************************************************************
//************************************************************************************


void ContactDomainLM3DCondition::CalculateDomainShapeN(GeneralVariables& rVariables)
{

    unsigned int ndi,ndj,ndk,ndl,ndm,ndn;

    if( this->Is(SELECTED) ){
      
      ndi=mContactVariables.nodes[0];
      ndj=mContactVariables.nodes[1];
      ndk=mContactVariables.slaves[0];
      ndl=mContactVariables.slaves[1];
      ndm=4;
      ndn=5;
      
    }
    else{

      ndi=mContactVariables.nodes[0];
      ndj=mContactVariables.nodes[1];
      ndk=mContactVariables.nodes[2];
      ndl=mContactVariables.slaves.back();
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

    // std::cout<<" DN_DX "<<DN_DX<<std::endl;
    // std::cout<<" Variables_DN_DX "<<rVariables.DN_DX<<std::endl;


    if( this->Is(SELECTED) ){

      //dN_dn:
      
      rVariables.Contact.dN_dn.resize(6);
      noalias(rVariables.Contact.dN_dn) = ZeroVector(6);
    
      rVariables.Contact.dN_dn[ndi]=(-1)*rVariables.Contact.CurrentBase[0].A/rVariables.Contact.CurrentBase[0].L;
      rVariables.Contact.dN_dn[ndj]=(-1)*rVariables.Contact.CurrentBase[0].B/rVariables.Contact.CurrentBase[0].L; 
      rVariables.Contact.dN_dn[ndk]= rVariables.Contact.CurrentBase[1].A/rVariables.Contact.CurrentBase[1].L; 
      rVariables.Contact.dN_dn[ndl]= rVariables.Contact.CurrentBase[1].B/rVariables.Contact.CurrentBase[1].L;


      //dN_dr:
    
      rVariables.Contact.dN_dr.resize(6);
      noalias(rVariables.Contact.dN_dr) = ZeroVector(6);
    
      rVariables.Contact.dN_dr[ndi]=(-1)*rVariables.Contact.ReferenceBase[0].A/rVariables.Contact.ReferenceBase[0].L;
      rVariables.Contact.dN_dr[ndj]=(-1)*rVariables.Contact.ReferenceBase[0].B/rVariables.Contact.ReferenceBase[0].L;
      rVariables.Contact.dN_dr[ndk]= rVariables.Contact.ReferenceBase[1].A/rVariables.Contact.ReferenceBase[1].L;
      rVariables.Contact.dN_dr[ndl]= rVariables.Contact.ReferenceBase[1].B/rVariables.Contact.ReferenceBase[1].L;;


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
      
      
      rVariables.Contact.dN_dn.resize(5);
      noalias(rVariables.Contact.dN_dn) = ZeroVector(5);
    
      rVariables.Contact.dN_dn[ndi]=(-1)*rVariables.Contact.CurrentBase[0].A/rVariables.Contact.CurrentBase[0].L;
      rVariables.Contact.dN_dn[ndj]=(-1)*rVariables.Contact.CurrentBase[1].A/rVariables.Contact.CurrentBase[1].L; 
      rVariables.Contact.dN_dn[ndk]=(-1)*rVariables.Contact.CurrentBase[2].A/rVariables.Contact.CurrentBase[2].L; 
      rVariables.Contact.dN_dn[ndl]= 1.0;


      //dN_dr:
    
      rVariables.Contact.dN_dr.resize(5);
      noalias(rVariables.Contact.dN_dr) = ZeroVector(5);
    
      rVariables.Contact.dN_dr[ndi]=(-1)*rVariables.Contact.ReferenceBase[0].A/rVariables.Contact.ReferenceBase[0].L;
      rVariables.Contact.dN_dr[ndj]=(-1)*rVariables.Contact.ReferenceBase[1].A/rVariables.Contact.ReferenceBase[1].L;
      rVariables.Contact.dN_dr[ndk]=(-1)*rVariables.Contact.ReferenceBase[2].A/rVariables.Contact.ReferenceBase[2].L;
      rVariables.Contact.dN_dr[ndl]= 1.0;


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

  //A.TsigmaP :
  FSigmaP(rVariables,rVariables.Contact.Tangent.A.Tsigma,rVariables.Contact.Tangent.CovariantBase.DirectionA,ndi,ndj,ndk,ndl,ndm,ndn);


  //B.TsigmaP :
  FSigmaP(rVariables,rVariables.Contact.Tangent.B.Tsigma,rVariables.Contact.Tangent.CovariantBase.DirectionB,ndi,ndj,ndk,ndl,ndm,ndn);

 
  //NsigmaP :
  FSigmaP(rVariables,rVariables.Contact.Nsigma,rVariables.Contact.CurrentSurface.Normal,ndi,ndj,ndk,ndl,ndm,ndn);


  rVariables.DN_DX = DN_DX;

}


//************************************************************************************
//************************************************************************************

void ContactDomainLM3DCondition::FSigmaP(GeneralVariables& rVariables, std::vector<Vector > &SigmaP, PointType& DirVector,unsigned int &ndi,unsigned int &ndj,unsigned int &ndk,unsigned int &ndl,unsigned int &ndm,unsigned int &ndn)
{
    if( this->Is(SELECTED) ){

      //node for computation / node for assignment (nd1,nd2,nd5,nd6)
      
      FSigmaPnd(rVariables, SigmaP, DirVector, ndi, ndi);
      
      FSigmaPnd(rVariables, SigmaP, DirVector, ndj, ndj);
      
      SigmaP[ndk] = ZeroVector(3);
      
      SigmaP[ndl] = ZeroVector(3);
      
      FSigmaPnd(rVariables, SigmaP, DirVector, ndk, ndm);

      FSigmaPnd(rVariables, SigmaP, DirVector, ndl, ndn);
      
    }
    else{
      
      //node for computation / node for assignment (nd1,nd3,nd2,nd5)

      FSigmaPnd(rVariables, SigmaP, DirVector, ndi, ndi);
  
      FSigmaPnd(rVariables, SigmaP, DirVector, ndj, ndj);
      
      FSigmaPnd(rVariables, SigmaP, DirVector, ndk, ndk);
      
      SigmaP[ndl] = ZeroVector(3);
      
      FSigmaPnd(rVariables, SigmaP, DirVector, ndl, ndm);
      
    }

}

//************************************************************************************
//************************************************************************************


void ContactDomainLM3DCondition::FSigmaPnd(GeneralVariables& rVariables, std::vector<Vector > &SigmaP, PointType& DirVector,unsigned int &ndi,unsigned int &ndj)
{
    //Computation with the ndi and storage to ndj
    SigmaP[ndj].resize(3);
    noalias(SigmaP[ndj]) = ZeroVector(3);

    //part1:
    SigmaP[ndj][0]= DirVector[0]*mContactVariables.ReferenceSurface.Normal[0]*(rVariables.StressVector[0]*rVariables.DN_DX(ndi,0)+rVariables.StressVector[3]*rVariables.DN_DX(ndi,1)+rVariables.StressVector[5]*rVariables.DN_DX(ndi,2))+
                    DirVector[0]*mContactVariables.ReferenceSurface.Normal[1]*(rVariables.StressVector[3]*rVariables.DN_DX(ndi,0)+rVariables.StressVector[1]*rVariables.DN_DX(ndi,1)+rVariables.StressVector[4]*rVariables.DN_DX(ndi,2))+
                    DirVector[0]*mContactVariables.ReferenceSurface.Normal[2]*(rVariables.StressVector[5]*rVariables.DN_DX(ndi,0)+rVariables.StressVector[4]*rVariables.DN_DX(ndi,1)+rVariables.StressVector[2]*rVariables.DN_DX(ndi,2));

    SigmaP[ndj][1]= DirVector[1]*mContactVariables.ReferenceSurface.Normal[1]*(rVariables.StressVector[3]*rVariables.DN_DX(ndi,0)+rVariables.StressVector[1]*rVariables.DN_DX(ndi,1)+rVariables.StressVector[4]*rVariables.DN_DX(ndi,2))+
                    DirVector[1]*mContactVariables.ReferenceSurface.Normal[0]*(rVariables.StressVector[0]*rVariables.DN_DX(ndi,0)+rVariables.StressVector[5]*rVariables.DN_DX(ndi,1)+rVariables.StressVector[5]*rVariables.DN_DX(ndi,2))+
                    DirVector[1]*mContactVariables.ReferenceSurface.Normal[2]*(rVariables.StressVector[5]*rVariables.DN_DX(ndi,0)+rVariables.StressVector[4]*rVariables.DN_DX(ndi,1)+rVariables.StressVector[2]*rVariables.DN_DX(ndi,2));

    SigmaP[ndj][2]= DirVector[2]*mContactVariables.ReferenceSurface.Normal[2]*(rVariables.StressVector[5]*rVariables.DN_DX(ndi,0)+rVariables.StressVector[4]*rVariables.DN_DX(ndi,1)+rVariables.StressVector[2]*rVariables.DN_DX(ndi,2))+
                    DirVector[2]*mContactVariables.ReferenceSurface.Normal[1]*(rVariables.StressVector[3]*rVariables.DN_DX(ndi,0)+rVariables.StressVector[1]*rVariables.DN_DX(ndi,1)+rVariables.StressVector[4]*rVariables.DN_DX(ndi,2))+
                    DirVector[2]*mContactVariables.ReferenceSurface.Normal[0]*(rVariables.StressVector[0]*rVariables.DN_DX(ndi,0)+rVariables.StressVector[3]*rVariables.DN_DX(ndi,1)+rVariables.StressVector[5]*rVariables.DN_DX(ndi,2));

    
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

    // std::cout<<" FBD [0] "<<FDB[0]<<std::endl;
    // std::cout<<" FBD [1] "<<FDB[1]<<std::endl;
    // std::cout<<" FBD [2] "<<FDB[2]<<std::endl;
    // std::cout<<" FBD [3] "<<FDB[3]<<std::endl;


    for(unsigned int i=0; i<3; i++)
    {
      SigmaP[ndj][i]+=  DirVector[0]*mContactVariables.ReferenceSurface.Normal[0]*(FDB[0][i])+
            	        DirVector[1]*mContactVariables.ReferenceSurface.Normal[1]*(FDB[1][i])+
                        DirVector[2]*mContactVariables.ReferenceSurface.Normal[2]*(FDB[2][i])+

                        DirVector[0]*mContactVariables.ReferenceSurface.Normal[1]*(FDB[3][i])+
                        DirVector[1]*mContactVariables.ReferenceSurface.Normal[0]*(FDB[4][i])+
                        DirVector[1]*mContactVariables.ReferenceSurface.Normal[2]*(FDB[5][i])+

                        DirVector[2]*mContactVariables.ReferenceSurface.Normal[1]*(FDB[6][i])+
                        DirVector[0]*mContactVariables.ReferenceSurface.Normal[2]*(FDB[7][i])+
                        DirVector[2]*mContactVariables.ReferenceSurface.Normal[0]*(FDB[8][i]);

    }



}

//***********************************************************************************
//************************************************************************************

double& ContactDomainLM3DCondition::CalculateIntegrationWeight(double& rIntegrationWeight)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if ( dimension == 2 ){   
      ElementType& MasterElement = mContactVariables.GetMasterElement();
      rIntegrationWeight *= MasterElement.GetProperties()[THICKNESS];
    }
	
    return rIntegrationWeight;
}


//************************************************************************************
//************************************************************************************

void ContactDomainLM3DCondition::CalculateNormalForce (double &F,GeneralVariables& rVariables,unsigned int& ndi,unsigned int& idir)
{    

    F  = rVariables.Contact.Multiplier.Normal*rVariables.Contact.dN_dn[ndi]*rVariables.Contact.CurrentSurface.Normal[idir];
    F *= (1.0/3.0) * rVariables.Contact.Tangent.EquivalentArea;

}

//************************************************************************************
//************************************************************************************


void ContactDomainLM3DCondition::CalculateTangentStickForce (double &F,GeneralVariables& rVariables,unsigned int& ndi,unsigned int& idir)
{

	if( rVariables.Contact.Options.Is(ContactDomainUtilities::COMPUTE_FRICTION_FORCES) )
	{
	        F = rVariables.Contact.CurrentGap.Normal * rVariables.Contact.CurrentSurface.Normal[idir] +
		    rVariables.Contact.Tangent.A.CurrentGap.Contravariant * rVariables.Contact.Tangent.CovariantBase.DirectionA[idir] +
 		    rVariables.Contact.Tangent.B.CurrentGap.Contravariant * rVariables.Contact.Tangent.CovariantBase.DirectionB[idir];

		  
		F *= ( rVariables.Contact.Tangent.A.Multiplier * rVariables.Contact.Tangent.A.dN_dt[ndi]+
		       rVariables.Contact.Tangent.B.Multiplier * rVariables.Contact.Tangent.B.dN_dt[ndi] );
		
		F += rVariables.Contact.dN_drn[ndi] * ( rVariables.Contact.Tangent.A.Multiplier * rVariables.Contact.Tangent.Contravariant.DirectionA[idir] + rVariables.Contact.Tangent.B.Multiplier * rVariables.Contact.Tangent.Contravariant.DirectionB[idir] );

		F *= (1.0/3.0) * rVariables.Contact.Tangent.EquivalentArea;
	}
	else{
		F=0.0;
	}
	
}

//************************************************************************************
//************************************************************************************


void ContactDomainLM3DCondition::CalculateTangentSlipForce (double &F,GeneralVariables& rVariables,unsigned int& ndi,unsigned int& idir)
{
         
	if( rVariables.Contact.Options.Is(ContactDomainUtilities::COMPUTE_FRICTION_FORCES) )
	{

	         F = (rVariables.Contact.Tangent.A.GapSign * rVariables.Contact.Tangent.CovariantBase.DirectionA[ndi] +
		      rVariables.Contact.Tangent.B.GapSign * rVariables.Contact.Tangent.CovariantBase.DirectionB[ndi] );
		  
		 F *= rVariables.Contact.CurrentGap.Normal * rVariables.Contact.CurrentSurface.Normal[idir];
		  
		 F += rVariables.Contact.dN_drn[ndi] * ( rVariables.Contact.Tangent.Contravariant.DirectionA[idir] +
							 rVariables.Contact.Tangent.Contravariant.DirectionB[idir] );

		 F *= ( rVariables.Contact.Multiplier.Normal * rVariables.Contact.FrictionCoefficient );

		 F *= (1.0/3.0) * rVariables.Contact.Tangent.EquivalentArea;
	}
	else{
		 F=0.0;
	}



}


//************************************************************************************
//************************************************************************************

void ContactDomainLM3DCondition::CalcContactStiffness (double &Kcont,GeneralVariables& rVariables,unsigned int& ndi,unsigned int& ndj,unsigned int& idir,unsigned int& jdir)
{
    Kcont=0;

    double athird = 1.0/3.0;
    
    //Normal contact contribution:
    //KI:       
    Kcont =  (athird*rVariables.Contact.Tangent.EquivalentHeigh/(rVariables.Contact.ContactFactor.Normal)) *
             (rVariables.Contact.dN_dn[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dn[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir] )-    
             (mContactDomain.LmN * rVariables.Contact.CurrentGap.Normal ) *
             (rVariables.Contact.Tangent.ContravariantBase.Metric(0,0)*rVariables.Contact.Tangent.A.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.Tangent.A.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+
	      rVariables.Contact.Tangent.ContravariantBase.Metric(0,1)*rVariables.Contact.Tangent.A.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.Tangent.B.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+
	      rVariables.Contact.Tangent.ContravariantBase.Metric(1,0)*rVariables.Contact.Tangent.B.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.Tangent.A.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+
	      rVariables.Contact.Tangent.ContravariantBase.Metric(1,1)*rVariables.Contact.Tangent.B.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.Tangent.B.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir])-
             (rVariables.Contact.Multiplier.Normal) *
             (rVariables.Contact.Tangent.A.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dn[ndj]*rVariables.Contact.Tangent.ContravariantBase.DirectionA[jdir]+
	      rVariables.Contact.Tangent.B.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dn[ndj]*rVariables.Contact.Tangent.ContravariantBase.DirectionB[jdir]+
	      rVariables.Contact.dN_dn[ndi]*rVariables.Contact.Tangent.ContravariantBase.DirectionA[idir]*rVariables.Contact.Tangent.A.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+
	      rVariables.Contact.dN_dn[ndi]*rVariables.Contact.Tangent.ContravariantBase.DirectionB[idir]*rVariables.Contact.Tangent.B.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir])-
             (rVariables.Contact.dN_dn[ndi] * rVariables.Contact.CurrentSurface.Normal[idir]) *
             (mContactDomain.CurrentTensil.TangentcnA*rVariables.Contact.Tangent.A.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+
	      mContactDomain.CurrentTensil.TangentcnB*rVariables.Contact.Tangent.B.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]);

    // rstiff = cte*geofct_n/tau_n*hdNn_a(in)*ne_a(idime)*hdNn_a(jn)*ne_a(jdime)
    
    // std::cout<<" ndi "<<ndi<<" ndj "<<ndj<<" i "<<idir<<" j "<<jdir<<std::endl;
    // std::cout<<" Kg "<<Kcont;
    // //std::cout<<" constant "<<constant<<" Kg "<<Kcont*constant;
    // double K1=Kcont;

    //KII:
    Kcont+= rVariables.Contact.dN_dn[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.Nsigma[ndj][jdir];


    //Stick contact contribution:
    if(rVariables.Contact.Options.Is(NOT_SLIP))
    {
      //std::cout<<" + stick ";
      if(rVariables.Contact.Options.Is(ContactDomainUtilities::COMPUTE_FRICTION_STIFFNESS))
        {
	  //std::cout<<"(mu_on)";
	  //KI:
	  double raux=rVariables.Contact.CurrentGap.Normal*rVariables.Contact.CurrentSurface.Normal[idir]+ 
	              rVariables.Contact.Tangent.A.CurrentGap.Contravariant*rVariables.Contact.Tangent.CovariantBase.DirectionA[idir]+ 
           	      rVariables.Contact.Tangent.B.CurrentGap.Contravariant*rVariables.Contact.Tangent.CovariantBase.DirectionB[idir]; //it must be GapcnTA (contravariant gap)
	  
	  //KI:
	  Kcont+= (raux*rVariables.Contact.Tangent.A.dN_dt[ndi]+rVariables.Contact.dN_drn[ndi]*rVariables.Contact.Tangent.ContravariantBase.DirectionA[idir]) *
	          ( (athird*rVariables.Contact.Tangent.EquivalentHeigh/(rVariables.Contact.ContactFactor.Normal) )*rVariables.Contact.dN_drn[ndj]*rVariables.Contact.Tangent.CovariantBase.DirectionA[jdir]+
		    rVariables.Contact.Multiplier.Normal*rVariables.Contact.CurrentSurface.Normal[jdir]*(rVariables.Contact.Tangent.CovariantBase.Metric(0,0)*rVariables.Contact.Tangent.A.dN_dt[ndj]+rVariables.Contact.Tangent.CovariantBase.Metric(0,1)*rVariables.Contact.Tangent.B.dN_dt[ndj]) )+
	    (raux*rVariables.Contact.Tangent.B.dN_dt[ndi]+rVariables.Contact.dN_drn[ndi]*rVariables.Contact.Tangent.ContravariantBase.DirectionB[idir])*
  	         ( (athird*rVariables.Contact.Tangent.EquivalentHeigh/(rVariables.Contact.ContactFactor.Normal) )*rVariables.Contact.dN_drn[ndj]*rVariables.Contact.Tangent.CovariantBase.DirectionB[jdir]+
		   rVariables.Contact.Multiplier.Normal*rVariables.Contact.CurrentSurface.Normal[jdir]*(rVariables.Contact.Tangent.CovariantBase.Metric(1,0)*rVariables.Contact.Tangent.A.dN_dt[ndj]+rVariables.Contact.Tangent.CovariantBase.Metric(1,1)*rVariables.Contact.Tangent.B.dN_dt[ndj]) )+
	    rVariables.Contact.Tangent.A.Multiplier*(rVariables.Contact.Tangent.A.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dn[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+
			    rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.Tangent.A.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+
			   (rVariables.Contact.Tangent.ContravariantBase.DirectionA[idir]*rVariables.Contact.Tangent.CovariantBase.DirectionA[jdir]+rVariables.Contact.Tangent.ContravariantBase.DirectionB[idir]*rVariables.Contact.Tangent.CovariantBase.DirectionB[jdir])*rVariables.Contact.Tangent.A.dN_dt[ndi]*rVariables.Contact.dN_drn[ndj]+
			   (rVariables.Contact.Tangent.A.CurrentGap.Contvariant*rVariables.Contact.Tangent.A.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+ rVariables.Contact.Tangent.B.CurrentGap.Contvariant*rVariables.Contact.Tangent.B.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir])*rVariables.Contact.Tangent.A.dN_dt(ndi)*rVariables.Contact.CurrentSurface.Normal(idir)-
			   raux*rVariables.Contact.Tangent.A.dN_dt(ndi)*rVariables.Contact.Tangent.CovariantBase.DirectionA(jdir)*rVariables.Contact.Tangent.A.dN_dt(ndj)-
			   raux*rVariables.Contact.Tangent.B.dN_dt(ndi)*rVariables.Contact.Tangent.CovariantBase.DirectionB(jdir)*rVariables.Contact.Tangent.A.dN_dt(ndj))-	    
	    rVariables.Contact.Tangent.B.Multiplier*(rVariables.Contact.Tangent.B.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dn[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+
			   rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.Tangent.B.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+
			   (rVariables.Contact.Tangent.ContravariantBase.DirectionA[idir]*rVariables.Contact.Tangent.CovariantBase.DirectionA[jdir]+rVariables.Contact.Tangent.ContravariantBase.DirectionB[idir]*rVariables.Contact.Tangent.CovariantBase.DirectionB[jdir])*rVariables.Contact.Tangent.B.dN_dt[ndi]*rVariables.Contact.dN_drn[ndj]+
			   (rVariables.Contact.Tangent.A.CurrentGap.Contvariant*rVariables.Contact.Tangent.A.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+ rVariables.Contact.Tangent.B.CurrentGap.Contvariant*rVariables.Contact.Tangent.B.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir])*rVariables.Contact.Tangent.B.dN_dt(ndi)*rVariables.Contact.CurrentSurface.Normal(idir)-
			   raux*rVariables.Contact.Tangent.A.dN_dt(ndi)*rVariables.Contact.Tangent.CovariantBase.DirectionA(jdir)*rVariables.Contact.Tangent.B.dN_dt(ndj)-
			   raux*rVariables.Contact.Tangent.B.dN_dt(ndi)*rVariables.Contact.Tangent.CovariantBase.DirectionB(jdir)*rVariables.Contact.Tangent.B.dN_dt(ndj));

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

            //KII:
	    Kcont+= ((rVariables.Contact.CurrentGap.Normal*rVariables.Contact.CurrentSurface.Normal[idir]+rVariables.Contact.Tangent.A.CurrentGap.Contravariant*rVariables.Contact.Tangent.CovariantBase.DirectionA[idir]+rVariables.Contact.Tangent.B.CurrentGap.Contravariant*rVariables.Contact.Tangent.CovariantBase.DirectionB[idir])*
	     (rVariables.Contact.Tangent.A.dN_dt[ndi]*rVariables.Contact.Tangent.A.Tsigma[ndj][jdir]+rVariables.Contact.Tangent.B.dN_dt[ndi]*rVariables.Contact.Tangent.B.Tsigma[ndj][jdir])+
	     (rVariables.Contact.dN_drn[ndi]*(rVariables.Contact.Tangent.ContravariantBase.DirectionA[idir]*rVariables.Contact.Tangent.A.Tsigma[ndj][jdir]+rVariables.Contact.Tangent.ContravariantBase.DirectionB[idir]*rVariables.Contact.Tangent.B.Tsigma[ndj][jdir])));

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


    Kcont *= (1.0/3.0) * rVariables.Contact.Tangent.EquivalentArea;
    
    //std::cout<<" Ks "<<Kcont-K1<<std::endl;

}


//************************************************************************************
//************************************************************************************


ContactDomainUtilities::PointType & ContactDomainLM3DCondition::CalculateCurrentTangent ( PointType &rTangent )
{

	unsigned int node1=mContactVariables.nodes[0];
	unsigned int node2=mContactVariables.nodes[1];

	PointType P1  =  GetGeometry()[node1].Coordinates();
	PointType P2  =  GetGeometry()[node2].Coordinates();
  
	//Set Reference Tangent
	rTangent=mContactUtilities.CalculateFaceTangent(rTangent,P1,P2);
  
	return rTangent;

}

//************************************************************************************
//************************************************************************************


inline bool ContactDomainLM3DCondition::CheckFictiousContacts(GeneralVariables& rVariables)
{

  bool real_contact = false;

  //Contact face segment node1-node2
  unsigned int node1=mContactVariables.nodes[0];
  unsigned int node2=mContactVariables.nodes[1];    
  unsigned int slave=mContactVariables.slaves.back();

  double offset_factor = rVariables.Contact.CurrentGap.Normal; 
  
  PointType PS  =  GetGeometry()[slave].Coordinates();

  PointType Normal =GetGeometry()[slave].FastGetSolutionStepValue(NORMAL); 
  double  Shrink             =1;//GetGeometry()[slave].FastGetSolutionStepValue(SHRINK_FACTOR);   
  PointType Offset =GetGeometry()[slave].FastGetSolutionStepValue(OFFSET);   
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
  WeakPointerVector<Element >& rNeighbours_n1 = GetGeometry()[node1].GetValue(NEIGHBOUR_ELEMENTS);
  //node2:
  WeakPointerVector<Element >& rNeighbours_n2 = GetGeometry()[node2].GetValue(NEIGHBOUR_ELEMENTS);

  unsigned int NumberOfNeighbours_n1 = rNeighbours_n1.size();
  unsigned int NumberOfNeighbours_n2 = rNeighbours_n2.size();

  bool is_inside_a = false;
  //following slave normal projection of the slave Sx1 and Sy1
  for(unsigned int i = 0; i < NumberOfNeighbours_n1; i++)
    {
      GeometryType::PointsArrayType& vertices=rNeighbours_n1[i].GetGeometry().Points();
  
      is_inside_a = mContactUtilities.CalculatePosition( vertices[0].X(), vertices[0].Y(),
							       vertices[1].X(), vertices[1].Y(),
							       vertices[2].X(), vertices[2].Y(),
							       Sx1, Sy1);

      if(is_inside_a)
	break;
    }

  if(!is_inside_a){
  
    for(unsigned int i = 0; i < NumberOfNeighbours_n2; i++)
      {
	GeometryType::PointsArrayType& vertices=rNeighbours_n2[i].GetGeometry().Points();
      
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
  for(unsigned int i = 0; i < NumberOfNeighbours_n1; i++)
    {
      GeometryType::PointsArrayType& vertices=rNeighbours_n1[i].GetGeometry().Points();
  
      is_inside_b = mContactUtilities.CalculatePosition( vertices[0].X(), vertices[0].Y(),
							       vertices[1].X(), vertices[1].Y(),
							       vertices[2].X(), vertices[2].Y(),
							       Mx1, My1);

      if(is_inside_b)
	break;
    }

  if(!is_inside_b){
  

    for(unsigned int i = 0; i < NumberOfNeighbours_n2; i++)
      {
	GeometryType::PointsArrayType& vertices=rNeighbours_n2[i].GetGeometry().Points();
      
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


