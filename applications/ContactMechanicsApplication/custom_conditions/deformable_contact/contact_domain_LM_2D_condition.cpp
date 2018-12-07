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
#include "custom_conditions/deformable_contact/contact_domain_LM_2D_condition.hpp"

#include "contact_mechanics_application_variables.h"


namespace Kratos
{


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ContactDomainLM2DCondition::ContactDomainLM2DCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : ContactDomainCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ContactDomainLM2DCondition::ContactDomainLM2DCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : ContactDomainCondition( NewId, pGeometry, pProperties )
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

ContactDomainLM2DCondition::ContactDomainLM2DCondition( ContactDomainLM2DCondition const& rOther)
    :ContactDomainCondition(rOther)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

ContactDomainLM2DCondition&  ContactDomainLM2DCondition::operator=(ContactDomainLM2DCondition const& rOther)
{
    ContactDomainCondition::operator=(rOther);

    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Condition::Pointer ContactDomainLM2DCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
  return Kratos::make_shared<ContactDomainLM2DCondition>(NewId, GetGeometry().Create( ThisNodes ), pProperties);
}

//************************************CLONE*******************************************
//************************************************************************************

Condition::Pointer ContactDomainLM2DCondition::Clone( IndexType NewId, NodesArrayType const& ThisNodes ) const
{
  return this->Create(NewId, ThisNodes, pGetProperties());
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

ContactDomainLM2DCondition::~ContactDomainLM2DCondition()
{
}

//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************


void ContactDomainLM2DCondition::SetMasterGeometry()

{
    KRATOS_TRY
    // std::cout<<" MASTER_ELEMENTS "<<GetValue(MASTER_ELEMENTS).size()<<" MASTER_NODES "<<GetValue(MASTER_NODES).size()<<std::endl;
    Element::ElementType& MasterElement = GetValue(MASTER_ELEMENTS).back();
    mContactVariables.SetMasterElement(MasterElement);

    Element::NodeType&    MasterNode   = GetValue(MASTER_NODES).back();
    mContactVariables.SetMasterNode(MasterNode);

    int  slave=-1;
    for(unsigned int i=0; i<MasterElement.GetGeometry().PointsNumber(); i++)
    {
        if(MasterNode.Id()==MasterElement.GetGeometry()[i].Id())
        {
	    slave=i;
	}
    }


    if(slave>=0)
    {
	// Clear nodes and slaves before push back quantities
	mContactVariables.nodes.resize(0);
	mContactVariables.slaves.resize(0);

        //NodesArrayType vertex;
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


	if( counter != 3 )
	  KRATOS_THROW_ERROR( std::invalid_argument, "Something is Wrong with MASTERNODE and contact order", "" )

	//Permute
	std::vector<unsigned int> permute (5);

	permute[0]=0;
	permute[1]=1;
	permute[2]=2;
	permute[3]=0;
	permute[4]=1;

	//reorder counter-clock-wise
	mContactVariables.nodes.push_back(permute[mContactVariables.slaves.back()+1]);
	mContactVariables.nodes.push_back(permute[mContactVariables.slaves.back()+2]);

	mContactVariables.SetMasterGeometry( MasterElement.GetGeometry() );

    }
    else
    {
        KRATOS_THROW_ERROR( std::invalid_argument, "MASTERNODE do not belongs to MASTER ELEMENT", "" )

    }

    // std::cout<<" CONTACT ("<<GetGeometry()[0].Id()<<")("<<GetGeometry()[1].Id()<<")("<<GetGeometry()[2].Id()<<")"<<std::endl;

    // std::cout<<" MASTER  ("<<MasterElement.GetGeometry()[0].Id()<<")("<<MasterElement.GetGeometry()[1].Id()<<")("<<MasterElement.GetGeometry()[2].Id()<<")"<<std::endl;
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************

void ContactDomainLM2DCondition::CalculatePreviousGap() //prediction of the lagrange multiplier
{

    //Contact face segment node1-node2
    unsigned int node1=mContactVariables.nodes[0];
    unsigned int node2=mContactVariables.nodes[1];
    unsigned int slave=mContactVariables.slaves.back();

    //Get Reference Normal
    PointType PP1  =  GetGeometry()[node1].Coordinates() - (GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT) -GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1) );
    PointType PP2  =  GetGeometry()[node2].Coordinates() - (GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1) );

    mContactVariables.ReferenceSurface.Normal=mContactUtilities.CalculateFaceNormal(mContactVariables.ReferenceSurface.Normal,PP1,PP2);

    //std::cout<<" COMPUTED NORMAL "<<mContactVariables.ReferenceSurface.Normal<<" CONDITION NORMAL "<<GetValue(NORMAL)<<std::endl;
    //mContactVariables.ReferenceSurface.Normal=GetValue(NORMAL);

    //std::cout<<" Got Normal ["<<this->Id()<<"] "<<mContactVariables.ReferenceSurface.Normal<<std::endl;

    //Set Reference Tangent
    mContactVariables.ReferenceSurface.Tangent=mContactUtilities.CalculateFaceTangent(mContactVariables.ReferenceSurface.Tangent,mContactVariables.ReferenceSurface.Normal);

    //4.- Compute Effective Gaps: (g^eff=g_n3+2*Tau*tn=2*Tau*Multiplier.Normal)

    //a.- Recover tensils from previous step: (it must be done once for each time step)

    //compare to auxiliar variables stored in the contact nodes to restore the LocalTensils
    //from the previous configuration

    // Element::NodeType&    MasterNode   = GetValue(MASTER_NODES).back();

    Condition::Pointer MasterCondition = GetValue(MASTER_CONDITION);


    //Get previous mechanics stored in the master node/condition
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    Matrix StressMatrix(dimension,dimension);
    noalias(StressMatrix)= ZeroMatrix(dimension,dimension);
    Matrix F(dimension,dimension);
    noalias(F)= ZeroMatrix(dimension,dimension);

    Vector StressVector;
    StressVector = MasterCondition->GetValue(CAUCHY_STRESS_VECTOR);  //it means that has been stored
    F            = MasterCondition->GetValue(DEFORMATION_GRADIENT);  //it means that has been stored

    // std::cout<<" Contact Condition["<<this->Id()<<"]"<<std::endl;
    // std::cout<<" Master Condition : "<<MasterCondition->Id()<<" Contact "<<std::endl;
    // std::cout<<" MC StressVector "<<StressVector<<std::endl;
    // std::cout<<" MC F "<<F<<std::endl;

    StressMatrix = MathUtils<double>::StressVectorToTensor( StressVector );

    //we are going to need F here from Cn-1 to Cn
    // F0 from C0 to Cn is need for the stress recovery on domain elements

    double detF =MathUtils<double>::Det(F);

    //b.- Compute the 1srt Piola Kirchhoff stress tensor  (P=J*CauchyStress*F^-T)
    ConstitutiveLaw Constitutive;
    Constitutive.TransformStresses(StressMatrix,F,detF,ConstitutiveLaw::StressMeasure_Cauchy,ConstitutiveLaw::StressMeasure_PK1);

    //Compute the tension (or traction) vector T=P*N (in the Reference configuration)

    //c.- Transform to 3 components
    Matrix StressMatrix3D(3,3);
    noalias(StressMatrix3D) = ZeroMatrix(3,3);
    for(unsigned int i=0; i<2; i++)
    {
	for(unsigned int j=0; j<2; j++)
	{
	    StressMatrix3D(i,j)=StressMatrix(i,j);
	}
    }



    //c.- Compute (n-1) normals, tangents and relative displacements from historic mX on boundaries

    //Previous normal and tangent:  n_n-1,t_n-1
    //Previous Position
    PointType PS  =  GetGeometry()[slave].Coordinates() - (GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,1) - GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,2) );
    PointType P1  =  GetGeometry()[node1].Coordinates() - (GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1) -GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,2) );
    PointType P2  =  GetGeometry()[node2].Coordinates() - (GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1) - GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,2) );

    //Compute Previous Normal
    mContactVariables.PreStepSurface.Normal=mContactUtilities.CalculateFaceNormal(mContactVariables.PreStepSurface.Normal,P1,P2);

    //std::cout<<" Pre Normal ["<<this->Id()<<"] "<<mContactVariables.PreStepSurface.Normal<<std::endl;

    // if((inner_prod(mContactVariables.PreStepSurface.Normal,mContactVariables.ReferenceSurface.Normal))<0) //to give the correct direction
    //     mContactVariables.PreStepSurface.Normal*=-1;

    // mContactVariables.PreStepSurface.Normal /= norm_2(mContactVariables.PreStepSurface.Normal);       //to be unitary

    // if(!(norm_2(mContactVariables.PreStepSurface.Normal)))
    // {
    //     mContactVariables.PreStepSurface.Normal=mContactVariables.ReferenceSurface.Normal;
    // }

    //Set Previous Tangent
    mContactVariables.PreStepSurface.Tangent=mContactUtilities.CalculateFaceTangent(mContactVariables.PreStepSurface.Tangent,mContactVariables.PreStepSurface.Normal);

    //Traction vector
    mContactVariables.TractionVector=prod(StressMatrix3D,mContactVariables.PreStepSurface.Normal);



    //Reference normal: n_n,t_n  -> mContactVariables.ReferenceSurface.Normal / mContactVariables.ReferenceSurface.Tangent

    //d.- Compute A_n-1,B_n-1,L_n-1

    //A_n-1, B_n-1, L_n-1:
    BaseLengths PreviousBase;
    mContactUtilities.CalculateBaseDistances(PreviousBase,P1,P2,PS,mContactVariables.PreStepSurface.Normal);

    //std::cout<<" L :"<<PreviousBase.L<<" A :"<<PreviousBase.A<<" B :"<<PreviousBase.B<<std::endl;

    //complete the computation of the stabilization gap
    double ContactFactor = mContactVariables.StabilizationFactor * PreviousBase.L;
    double ContactFactorTangent = ContactFactor * GetProperties()[TANGENTIAL_PENALTY_RATIO];

    //std::cout<<" Tau "<<ContactFactor<<std::endl;

    //e.-obtain the (g_N)3 and (g_T)3 for the n-1 configuration

    PointType DS  =  GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,1)-GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,2);
    PointType D1  =  GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1)-GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,2);
    PointType D2  =  GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1)-GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,2);


    mContactVariables.PreviousGap.Normal  = 0;
    mContactVariables.PreviousGap.Tangent = 0;

    //(g_N)3
    mContactVariables.PreviousGap.Normal = inner_prod((PS-P1),mContactVariables.PreStepSurface.Normal);
    mContactVariables.PreviousGap.Tangent = mContactVariables.PreviousGap.Normal;

    mContactVariables.PreviousGap.Normal*= inner_prod(mContactVariables.ReferenceSurface.Normal,mContactVariables.PreStepSurface.Normal);

    mContactVariables.PreviousGap.Normal+=inner_prod(mContactVariables.ReferenceSurface.Normal,(D1*(-PreviousBase.A/PreviousBase.L)));
    mContactVariables.PreviousGap.Normal+=inner_prod(mContactVariables.ReferenceSurface.Normal,(D2*(-PreviousBase.B/PreviousBase.L)));
    mContactVariables.PreviousGap.Normal+=inner_prod(mContactVariables.ReferenceSurface.Normal,DS);

    //(g_T)3
    mContactVariables.PreviousGap.Tangent*= inner_prod(mContactVariables.ReferenceSurface.Tangent,mContactVariables.PreStepSurface.Normal);

    mContactVariables.PreviousGap.Tangent+=inner_prod(mContactVariables.ReferenceSurface.Tangent,(D1*(-PreviousBase.A/PreviousBase.L)));
    mContactVariables.PreviousGap.Tangent+=inner_prod(mContactVariables.ReferenceSurface.Tangent,(D2*(-PreviousBase.B/PreviousBase.L)));
    mContactVariables.PreviousGap.Tangent+=inner_prod(mContactVariables.ReferenceSurface.Tangent,DS);


    //d_n-1=X_n - X_n-1

    //f.- get total effective gap as: gap_n^eff=gap_n+(PreviousTimeStep/CurrentTimeStep)*(gap_n-gap_n-1)

    //gap_n-1 (in function of the n-1 position of hte other node) gap_n-1=(g_N)3_n-1+2*Tau*tn_n-1

    double NormalTensil=0,TangentTensil=0;


    //g.- Compute normal component of the tension vector:   (tn=n·P·N)
    //NormalTensil=inner_prod(mContactVariables.PreStepSurface.Normal,mContactVariables.TractionVector);
    NormalTensil=inner_prod(mContactVariables.ReferenceSurface.Normal,mContactVariables.TractionVector);

    //h.- Compute tangent component of the tension vector:  (tt=t·P·N)
    //TangentTensil=inner_prod(mContactVariables.PreStepSurface.Tangent,mContactVariables.TractionVector);
    TangentTensil=inner_prod(mContactVariables.ReferenceSurface.Tangent,mContactVariables.TractionVector);

    mContactVariables.PreviousGap.Normal  += 2 * ContactFactor * NormalTensil;
    mContactVariables.PreviousGap.Tangent += 2 * ContactFactorTangent * TangentTensil;


    //std::cout<<"ConditionID:  "<<this->Id()<<" -> Previous Tractions [tN:"<<NormalTensil<<", tT:"<<TangentTensil<<"] "<<std::endl;
    //std::cout<<" Previous Gaps [gN:"<<mContactVariables.PreviousGap.Normal<<", gT:"<<mContactVariables.PreviousGap.Tangent<<"] "<<std::endl;
}


//**********************************COMPUTE TAU STAB**********************************
//************************************************************************************


void ContactDomainLM2DCondition::CalculateContactFactor( ProcessInfo& rCurrentProcessInfo )
{
    //Initilialize Tau for the stabilization
    double alpha_stab = 0.1;
    alpha_stab = GetProperties()[TAU_STAB];

    ElementType& rMasterElement = mContactVariables.GetMasterElement();

    //Look at the nodes, get the slave and get the Emin

    //Contact face segment node1-node2
    unsigned int slave = mContactVariables.slaves.back();

    const Properties& SlaveProperties  = GetGeometry()[slave].GetValue(NEIGHBOUR_ELEMENTS)[0].GetProperties();
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

void ContactDomainLM2DCondition::CalculateExplicitFactors(ConditionVariables& rVariables, ProcessInfo& rCurrentProcessInfo)
{

    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    //Contact face segment node1-node2
    unsigned int node1=mContactVariables.nodes[0];
    unsigned int node2=mContactVariables.nodes[1];
    unsigned int slave=mContactVariables.slaves.back();

    // std::cout<<" ************ CONTACT ELEMENT "<<this->Id()<<" ************* "<<std::endl;
    // std::cout<<std::endl;

    // Element::ElementType& MasterElement = GetValue(MASTER_ELEMENTS).back();

    // std::cout<<" master element "<<MasterElement.Id()<<std::endl;
    // std::cout<<" Elastic Modulus "<<MasterElement.GetProperties()[YOUNG_MODULUS]<<std::endl;

    // std::cout<<" Nodes 1,2,S "<<node1<<" "<<node2<<" "<<slave<<std::endl;

    // std::cout<<" Nodes ID  "<<GetGeometry()[0].Id()<<" "<<GetGeometry()[1].Id()<<" "<<GetGeometry()[2].Id()<<std::endl;

    // std::cout<<" Master Nodes ID  "<<(*mpMasterGeometry)[0].Id()<<" "<<(*mpMasterGeometry)[1].Id()<<" "<<(*mpMasterGeometry)[2].Id()<<std::endl;

    //1.- Compute tension vector:  (must be updated each iteration)
    Matrix StressMatrix(dimension,dimension);

    //a.- Assign initial 2nd Piola Kirchhoff stress:
    StressMatrix=MathUtils<double>::StressVectorToTensor( rVariables.StressVector );

    // UL
    //b.- Compute the 1srt Piola Kirchhoff stress tensor
    ConstitutiveLaw Constitutive;
    StressMatrix = Constitutive.TransformStresses(StressMatrix, rVariables.F, rVariables.detF, ConstitutiveLaw::StressMeasure_Cauchy, ConstitutiveLaw::StressMeasure_PK1);

    // UTL
    //b.- Compute the 1srt Piola Kirchhoff stress tensor  (P=F·S)
    //StressMatrix = Constitutive.TransformStresses(StressMatrix, rVariables.F, rVariables.detF, ConstitutiveLaw::StressMeasure_PK2, ConstitutiveLaw::StressMeasure_PK1);
    //StressMatrix=prod(rVariables.F,StressMatrix);

    //b.- Transform to 3 components
    Matrix StressMatrix3D(3,3);
    noalias(StressMatrix3D) = ZeroMatrix(3,3);
    for(unsigned int i=0; i<2; i++)
      {
	for(unsigned int j=0; j<2; j++)
	  {
	    StressMatrix3D(i,j)=StressMatrix(i,j);
	  }
      }

    //c.- Compute the tension (or traction) vector T=P*N (in the Reference configuration)
    mContactVariables.TractionVector=prod(StressMatrix3D,mContactVariables.ReferenceSurface.Normal);


    //d.- Compute the Current Normal and Tangent

    PointType PS  =  GetGeometry()[slave].Coordinates();
    PointType P1  =  GetGeometry()[node1].Coordinates();
    PointType P2  =  GetGeometry()[node2].Coordinates();

    //compute the current normal vector
    rVariables.Contact.CurrentSurface.Normal=mContactUtilities.CalculateFaceNormal(rVariables.Contact.CurrentSurface.Normal,P1,P2);

    //std::cout<<" Current Normal "<<rVariables.Contact.CurrentSurface.Normal<<std::endl;

    // if(double(inner_prod(rVariables.Contact.CurrentSurface.Normal,mContactVariables.ReferenceSurface.Normal))<0) //to give the correct direction
    //     rVariables.Contact.CurrentSurface.Normal*=-1;


    // if(!(norm_2(rVariables.Contact.CurrentSurface.Normal)))
    //   rVariables.Contact.CurrentSurface.Normal = mContactVariables.ReferenceSurface.Normal;
    // else
    //   rVariables.Contact.CurrentSurface.Normal /= norm_2(rVariables.Contact.CurrentSurface.Normal);  //to be unitary

    //compute the current tangent vector
    //rVariables.Contact.CurrentSurface.Tangent=mContactUtilities.CalculateFaceTangent(rVariables.Contact.CurrentSurface.Tangent,P1,P2);
    rVariables.Contact.CurrentSurface.Tangent=mContactUtilities.CalculateFaceTangent(rVariables.Contact.CurrentSurface.Tangent,rVariables.Contact.CurrentSurface.Normal);


    //Current normal:   mContactVariables.ReferenceSurface.Normal

    //2.- Compute normal component of the tension vector:   (tn=n·P·N)
    rVariables.Contact.CurrentTensil.Normal=inner_prod(rVariables.Contact.CurrentSurface.Normal,mContactVariables.TractionVector);

    //3.- Compute tangent component of the tension vector:  (tt=t·P·N)
    rVariables.Contact.CurrentTensil.Tangent=inner_prod(rVariables.Contact.CurrentSurface.Tangent,mContactVariables.TractionVector);

    //4.- Compute Effective Gaps: (g^eff=g_n3+2*Tau*tn=2*Tau*Multiplier.Normal)

    //Reference normal: n_n,t_n  -> mContactVariables.ReferenceSurface.Normal / rVariables.Contact.Tangent
    //Current normal:   n,t      -> rVariables.Contact.CurrentSurface.Normal /  rVariables.Contact.CurrentSurface.Tangent


    // std::cout<<" 2nd PK stress "<<rVariables.StressVector<<std::endl;
    // std::cout<<" 1st PK stress "<<StressMatrix3D<<std::endl;
    // std::cout<<" DeformationGradientF "<<rVariables.F<<std::endl;
    // std::cout<<" Traction  Vector "<<mContactVariables.TractionVector<<std::endl;


    // std::cout<<" reference face  normal  "<<mContactVariables.ReferenceSurface.Normal<<std::endl;
    // std::cout<<" reference face  tangent  "<<mContactVariables.ReferenceSurface.Tangent<<std::endl;

    // std::cout<<" current face  normal  "<<rVariables.Contact.CurrentSurface.Normal<<std::endl;
    // std::cout<<" current face  tangent  "<<rVariables.Contact.CurrentSurface.Tangent<<std::endl;



    //d.- Compute A_n,B_n,L_n
    rVariables.Contact.ReferenceBase.resize(1);
    rVariables.Contact.CurrentBase.resize(1);

    //a, b, l:
    mContactUtilities.CalculateBaseDistances (rVariables.Contact.CurrentBase[0],P1,P2,PS,rVariables.Contact.CurrentSurface.Normal);

    //Write Current Positions:
    // std::cout<<" Current position node 1 "<<P1<<std::endl;
    // std::cout<<" Current position node 2 "<<P2<<std::endl;
    // std::cout<<" Current position node s "<<PS<<std::endl;


    //A, B, L:

    PS =  GetGeometry()[slave].Coordinates() - ( GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,1) );
    P1 =  GetGeometry()[node1].Coordinates() - ( GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1) );
    P2 =  GetGeometry()[node2].Coordinates() - ( GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1) );

    mContactUtilities.CalculateBaseDistances (rVariables.Contact.ReferenceBase[0],P1,P2,PS,mContactVariables.ReferenceSurface.Normal);


    //complete the computation of the stabilization gap
    rVariables.Contact.ContactFactor.Normal  =  mContactVariables.StabilizationFactor * rVariables.Contact.ReferenceBase[0].L;
    rVariables.Contact.ContactFactor.Tangent =  rVariables.Contact.ContactFactor.Normal * GetProperties()[TANGENTIAL_PENALTY_RATIO];

    //e.-obtain the (g_N)3 and (g_T)3 for the n configuration

    //Write Reference Positions:
    // std::cout<<" Reference position node 1 "<<P1<<std::endl;
    // std::cout<<" Reference position node 2 "<<P2<<std::endl;
    // std::cout<<" Reference position node s "<<PS<<std::endl;

    double ReferenceGapN = inner_prod((PS - P1), mContactVariables.ReferenceSurface.Normal);

    // std::cout<<" Reference GAP "<<ReferenceGapN<<std::endl;

    double ReferenceGapT = ReferenceGapN;

    double H = ReferenceGapN;

    PointType DS  =  GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT)-GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,1);
    PointType D1  =  GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT)-GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1);
    PointType D2  =  GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT)-GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1);

    //Write Displacements:
    // std::cout<<" displacement node 1 "<<D1<<std::endl;
    // std::cout<<" displacement node 2 "<<D2<<std::endl;
    // std::cout<<" displacement node s "<<DS<<std::endl;

    //(g_N)3
    ReferenceGapN*=inner_prod(rVariables.Contact.CurrentSurface.Normal,mContactVariables.ReferenceSurface.Normal);

    ReferenceGapN+=inner_prod(rVariables.Contact.CurrentSurface.Normal,(D1*(-rVariables.Contact.ReferenceBase[0].A/rVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapN+=inner_prod(rVariables.Contact.CurrentSurface.Normal,(D2*(-rVariables.Contact.ReferenceBase[0].B/rVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapN+=inner_prod(rVariables.Contact.CurrentSurface.Normal,DS);

    //(g_T)3
    ReferenceGapT*=inner_prod(rVariables.Contact.CurrentSurface.Tangent,mContactVariables.ReferenceSurface.Normal);

    ReferenceGapT+=inner_prod(rVariables.Contact.CurrentSurface.Tangent,(D1*(-rVariables.Contact.ReferenceBase[0].A/rVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapT+=inner_prod(rVariables.Contact.CurrentSurface.Tangent,(D2*(-rVariables.Contact.ReferenceBase[0].B/rVariables.Contact.ReferenceBase[0].L)));
    ReferenceGapT+=inner_prod(rVariables.Contact.CurrentSurface.Tangent,DS);



    //std::cout<<" L :"<<rVariables.Contact.ReferenceBase[0].L<<" A :"<<rVariables.Contact.ReferenceBase[0].A<<" B :"<<rVariables.Contact.ReferenceBase[0].B<<std::endl;
    // std::cout<<" gN3 ref "<<ReferenceGapN<<std::endl;


    rVariables.Contact.CurrentGap.Normal=ReferenceGapN; //(g_N)3 -- needed in the Kcont1 computation
    rVariables.Contact.CurrentGap.Tangent=ReferenceGapT; //(g_T)3 -- needed in the Kcont1 computation

    //f.- get total effective gap as: gap_n^eff=gap_n+(PreviousTimeStep/CurrentTimeStep)*(gap_n-gap_n-1)

    //gap_n   (in function of the n position of the other node) gap_n=(g_N)3+2*Tau*tn_n

    ReferenceGapN+= 2 * rVariables.Contact.ContactFactor.Normal * rVariables.Contact.CurrentTensil.Normal;
    ReferenceGapT+= 2 * rVariables.Contact.ContactFactor.Tangent * rVariables.Contact.CurrentTensil.Tangent;

    //std::cout<<" ReferenceGapN "<<ReferenceGapN<<" ContactFactor "<<rVariables.Contact.ContactFactor.Normal<<" TensilNormal "<<rVariables.Contact.CurrentTensil.Normal<<std::endl;

    rVariables.Contact.TangentialGapSign=1;

    if( rVariables.Contact.CurrentGap.Tangent < 0 )
    {
        rVariables.Contact.TangentialGapSign=-1;
    }


    if(H==0) rVariables.Contact.TangentialGapSign=0;

    //CORRECTION: to skip change on diagonals problems //convergence problems !!! take care;
    if(rVariables.Contact.CurrentGap.Normal>0 && ReferenceGapN<0)
      {
	//look at the magnitud
	if(fabs(rVariables.Contact.CurrentGap.Normal) > 2*fabs(ReferenceGapN)){
	  ReferenceGapN = 0; //rVariables.Contact.CurrentGap.Normal +  rVariables.Contact.ContactFactor.Normal * rVariables.Contact.CurrentTensil.Normal;
	}

      }

    //std::cout<<" Tensil "<<rVariables.Contact.CurrentTensil.ReferenceSurface.Normal<<" Tau "<<rVariables.Contact.ContactFactor.Normal<<" product "<<2*rVariables.Contact.ContactFactor*rVariables.Contact.CurrentTensil.ReferenceSurface.Normal<<std::endl;
    //std::cout<<" gN3 ref total "<<ReferenceGapN<<std::endl;

    //5.- Compute (Lagrange) Multipliers

    //From effective gaps set active contact domain:

    double EffectiveGapN = ReferenceGapN;
    double EffectiveGapT = ReferenceGapT;

    double CurrentTimeStep  = rCurrentProcessInfo[DELTA_TIME];
    ProcessInfo& rPreviousProcessInfo = rCurrentProcessInfo.GetPreviousSolutionStepInfo();
    double PreviousTimeStep = rPreviousProcessInfo[DELTA_TIME];


    if(mContactVariables.PreviousGap.Normal!=0 && mContactVariables.IterationCounter<1){
      //std::cout<<" Effective prediction first iteration +:"<<(CurrentTimeStep/PreviousTimeStep)*(ReferenceGapN-mContactVariables.PreviousGap.Normal)<<" PreviousGap.Normal "<<mContactVariables.PreviousGap.Normal<<" ReferenceGapN "<<ReferenceGapN<<std::endl;

      EffectiveGapN+=(CurrentTimeStep/PreviousTimeStep)*(ReferenceGapN-mContactVariables.PreviousGap.Normal);
    }

    //std::cout<<" EffectiveGapN "<<EffectiveGapN<<" PreviousEffectiveGapN "<<ReferenceGapN<<std::endl;
    //only in the first iteration:
    //mContactVariables.PreviousGap.Normal=ReferenceGapN;


    if(mContactVariables.PreviousGap.Tangent!=0 && mContactVariables.IterationCounter<1){
      //std::cout<<" Effective prediction first iteration +:"<<(CurrentTimeStep/PreviousTimeStep)*(ReferenceGapT-mContactVariables.PreviousGap.Tangent)<<" PreviousGap.Tangent "<<mContactVariables.PreviousGap.Tangent<<" ReferenceGapT "<<ReferenceGapT<<std::endl;

      EffectiveGapT+=(CurrentTimeStep/PreviousTimeStep)*(ReferenceGapT-mContactVariables.PreviousGap.Tangent);
    }

    //std::cout<<" EffectiveGapT "<<EffectiveGapT<<" PreviousEffectiveGapT "<<ReferenceGapT<<std::endl;
    //only in the first iteration:
    //mContactVariables.PreviousGap.Tangent=ReferenceGapT;


    // std::cout<<" ConditionID:  "<<this->Id()<<" -> Previous Gap [gN:"<<ReferenceGapN<<", gT:"<<ReferenceGapT<<"] "<<std::endl;
    // std::cout<<" -> Effective Gap [gN:"<<EffectiveGapN<<", gT:"<<EffectiveGapT<<"] "<<std::endl;

    //std::cout<<" PreTimeStep "<<Time.PreStep<<" TimeStep "<<Time.Step<<std::endl;

    //CHECK IF THE ELEMENT IS ACTIVE:

    rVariables.Contact.Options.Set(SLIP,false);

    // if( rVariables.Contact.CurrentGap.Normal>0 && EffectiveGapN<0 && rVariables.Contact.CurrentGap.Normal>rVariables.Contact.ReferenceBase[0].L*0.001){
    //   EffectiveGapN =  rVariables.Contact.CurrentGap.Normal;
    // }

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
    // if(fabs(EffectiveGapN)<= 1e-15 && fabs(EffectiveGapN)<= 1e-15)
    //   EffectiveGapN = 0;
    //decimal correction from tension vector calculation


    if(EffectiveGapN<=0)   //if(EffectiveGap<0){
    {

        //Calculate Friction coeffitients with a regularized law

        //Initialize friction parameter
        rVariables.Contact.FrictionCoefficient = 0;

        //Tangent velocity and stablish friction parameter
        PointType TangentVelocity (3,0.0);

        //Calculate Relative Velocity
        this->CalculateRelativeVelocity(rVariables,TangentVelocity,rCurrentProcessInfo);

        //Calculate Friction Coefficient
        this->CalculateFrictionCoefficient(rVariables,TangentVelocity);

	//std::cout<<" Friction Coefficient ["<<this->Id()<<"]"<<rVariables.Contact.FrictionCoefficient<<" Sign "<<rVariables.Contact.TangentialGapSign<<std::endl;

        rVariables.Contact.Options.Set(ACTIVE,true);  //normal contact active

	//std::cout<<" EffectiveGapT ["<<this->Id()<<"] :"<<EffectiveGapT<<std::endl;
	//std::cout<<" mu * EffectiveGapN ["<<this->Id()<<"] :"<<rVariables.Contact.FrictionCoefficient*fabs(EffectiveGapN)<<", mu :"<<rVariables.Contact.FrictionCoefficient<<std::endl;

        if(fabs(EffectiveGapT)<rVariables.Contact.FrictionCoefficient*fabs(EffectiveGapN))
        {
	    rVariables.Contact.Options.Set(SLIP,false); //contact stick case active

        }
        else
        {
	    rVariables.Contact.Options.Set(SLIP,true);  //contact slip  case active
	    //std::cout<<" EffectiveGapT ["<<this->Id()<<"] :"<<EffectiveGapT<<" Reference Gap "<<ReferenceGapT<<" mu*Fn "<<(rVariables.Contact.FrictionCoefficient*fabs(EffectiveGapN))<<std::endl;
        }


	//this->Set(ACTIVE); not here, if is reset is not going to enter here anymore
    }
    else
    {
	rVariables.Contact.Options.Set(ACTIVE,false); //normal contact not active
	//this->Reset(ACTIVE); not here, if is reset is not going to enter here anymore
    }


    //rVariables.Contact.Options.Set(SLIP,false); //impose stick
    //rVariables.Contact.Options.Set(SLIP,true); //impose slip

    //From total current gap compute multipliers:

    //rVariables.Contact.Multiplier.Normal = EffectiveGap*(1./(2.0*rVariables.Contact.ContactFactor.Normal)); //posible computation of the Lagrange Multiplier
    rVariables.Contact.Multiplier.Normal =rVariables.Contact.CurrentTensil.Normal;
    rVariables.Contact.Multiplier.Normal+=rVariables.Contact.CurrentGap.Normal*(1./(2.0*rVariables.Contact.ContactFactor.Normal));

    rVariables.Contact.Multiplier.Tangent =rVariables.Contact.CurrentTensil.Tangent;
    rVariables.Contact.Multiplier.Tangent+=rVariables.Contact.CurrentGap.Tangent*(1./(2.0*rVariables.Contact.ContactFactor.Tangent));

    //std::cout<<" Mt "<<rVariables.Contact.Multiplier.Tangent<<" Mn*mu "<<rVariables.Contact.Multiplier.Normal*rVariables.Contact.FrictionCoefficient<<std::endl;

    if( rVariables.Contact.Multiplier.Normal < 0 ){

      if( rVariables.Contact.Multiplier.Tangent < 0 )  //set the sign of the Lagrange Multiplier
	{
	  rVariables.Contact.TangentialGapSign = 1;
	}
      else
	{
	  rVariables.Contact.TangentialGapSign = -1;
	}

    }

    //std::cout<<" Tangent "<<rVariables.Contact.Multiplier.Tangent<<std::endl;

    if(rVariables.Contact.Options.Is(ACTIVE) && rVariables.Contact.CurrentGap.Normal>0){
      // int active = 0;
      // if(this->Is(ACTIVE))
      // 	active = 1;
      // std::cout<<" Condition ["<<this->Id()<<"]:  Active "<<active<<" Effective GapN "<<EffectiveGapN<<" Multiplier.Normal "<<rVariables.Contact.Multiplier.Normal<<" CurrentTensil.N "<<rVariables.Contact.CurrentTensil.Normal<<" GapN "<<rVariables.Contact.CurrentGap.Normal<<" ReferenceGapN "<<ReferenceGapN<<" Tau "<<rVariables.Contact.ContactFactor.Normal<<" iteration "<<mContactVariables.IterationCounter<<std::endl;
    }

    //set contact normal
    array_1d<double, 3> &ContactNormal  = GetGeometry()[slave].FastGetSolutionStepValue(CONTACT_NORMAL);

    for(unsigned int i=0; i<3; i++)
      ContactNormal[i] = rVariables.Contact.CurrentSurface.Normal[i];


    if(mContactVariables.IterationCounter < 1)
      mContactVariables.IterationCounter += 1;

}


//************************************************************************************
//************************************************************************************


void ContactDomainLM2DCondition::CalculateDomainShapeN(ConditionVariables& rVariables)
{

    unsigned int ndi=mContactVariables.nodes[0];
    unsigned int ndj=mContactVariables.nodes[1];
    unsigned int ndk=mContactVariables.slaves[0];
    unsigned int ndr=3;

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

    //dN_dn:

    rVariables.Contact.dN_dn=ZeroVector(4);

    rVariables.Contact.dN_dn[ndi]=(-1)*rVariables.Contact.CurrentBase[0].A/rVariables.Contact.CurrentBase[0].L;
    rVariables.Contact.dN_dn[ndj]=(-1)*rVariables.Contact.CurrentBase[0].B/rVariables.Contact.CurrentBase[0].L;
    rVariables.Contact.dN_dn[ndk]= 1.0;


    //dN_drn:

    rVariables.Contact.dN_drn=ZeroVector(4);

    rVariables.Contact.dN_drn[ndi]=(-1)*rVariables.Contact.ReferenceBase[0].A/rVariables.Contact.ReferenceBase[0].L;
    rVariables.Contact.dN_drn[ndj]=(-1)*rVariables.Contact.ReferenceBase[0].B/rVariables.Contact.ReferenceBase[0].L;
    rVariables.Contact.dN_drn[ndk]= 1.0;


    //dN_dt:

    rVariables.Contact.dN_dt=ZeroVector(4);

    rVariables.Contact.dN_dt[ndi]=   1.0/rVariables.Contact.CurrentBase[0].L;
    rVariables.Contact.dN_dt[ndj]=(-1.0)/rVariables.Contact.CurrentBase[0].L;

    // std::cout<<" dN_dn "<<rVariables.Contact.dN_dn<<std::endl;
    // std::cout<<" dN_drn"<<rVariables.Contact.dN_drn<<std::endl;
    // std::cout<<" dN_dt "<<rVariables.Contact.dN_dt<<std::endl;

    //TsigmaP :
    rVariables.Contact.Tsigma.resize(4);
    FSigmaP(rVariables,rVariables.Contact.Tsigma,rVariables.Contact.CurrentSurface.Tangent,ndi,ndj,ndk,ndr);

    // for(unsigned int i=0; i<rVariables.Contact.Tsigma.size(); i++)
    //   std::cout<<" i: "<<i<<" Tsigma "<<rVariables.Contact.Tsigma[i]<<std::endl;


    //NsigmaP :
    rVariables.Contact.Nsigma.resize(4);
    FSigmaP(rVariables,rVariables.Contact.Nsigma,rVariables.Contact.CurrentSurface.Normal,ndi,ndj,ndk,ndr);

    // for(unsigned int i=0; i<rVariables.Contact.Nsigma.size(); i++)
    //   std::cout<<" i: "<<i<<" Nsigma "<<rVariables.Contact.Nsigma[i]<<std::endl;


    rVariables.DN_DX = DN_DX;

}


//************************************************************************************
//************************************************************************************

void ContactDomainLM2DCondition::FSigmaP(ConditionVariables& rVariables, std::vector<Vector > &SigmaP, PointType& DirVector,unsigned int &ndi,unsigned int &ndj,unsigned int &ndk,unsigned int &ndr)
{
    //Computation with the ndi and storage to ndj

    // std::cout<<" DirVector "<<DirVector<<std::endl;
    // std::cout<<" StressVector "<<rVariables.StressVector<<std::endl;

    FSigmaPnd(rVariables, SigmaP, DirVector, ndi, ndi);

    FSigmaPnd(rVariables, SigmaP, DirVector, ndj, ndj);

    SigmaP[ndk]=ZeroVector(2);

    FSigmaPnd(rVariables, SigmaP, DirVector, ndk, ndr);

}

//************************************************************************************
//************************************************************************************


void ContactDomainLM2DCondition::FSigmaPnd(ConditionVariables& rVariables, std::vector<Vector > &SigmaP, PointType& DirVector,unsigned int &ndi,unsigned int &ndj)
{
    //Computation with the ndi and storage to ndj
    SigmaP[ndj]=ZeroVector(2);

    //part1:
    SigmaP[ndj][0]= DirVector[0]*mContactVariables.ReferenceSurface.Normal[0]*(rVariables.StressVector[0]*rVariables.DN_DX(ndi,0)+rVariables.StressVector[2]*rVariables.DN_DX(ndi,1))+ DirVector[0]*mContactVariables.ReferenceSurface.Normal[1]*(rVariables.StressVector[1]*rVariables.DN_DX(ndi,1)+rVariables.StressVector[2]*rVariables.DN_DX(ndi,0));


    SigmaP[ndj][1]= DirVector[1]*mContactVariables.ReferenceSurface.Normal[1]*(rVariables.StressVector[1]*rVariables.DN_DX(ndi,1)+rVariables.StressVector[2]*rVariables.DN_DX(ndi,0))+ DirVector[1]*mContactVariables.ReferenceSurface.Normal[0]*(rVariables.StressVector[0]*rVariables.DN_DX(ndi,0)+rVariables.StressVector[2]*rVariables.DN_DX(ndi,1));

    //part2:
    std::vector<Vector> FD(4);

    FD[0].resize(4);
    FD[0][0]=(rVariables.F(0,0)*rVariables.ConstitutiveMatrix(0,0)+rVariables.F(0,1)*rVariables.ConstitutiveMatrix(2,0));
    FD[0][1]=(rVariables.F(0,0)*rVariables.ConstitutiveMatrix(0,1)+rVariables.F(0,1)*rVariables.ConstitutiveMatrix(2,1));
    FD[0][2]=(rVariables.F(0,0)*rVariables.ConstitutiveMatrix(0,2)+rVariables.F(0,1)*rVariables.ConstitutiveMatrix(2,2));
    FD[0][3]=(rVariables.F(0,0)*rVariables.ConstitutiveMatrix(0,2)+rVariables.F(0,1)*rVariables.ConstitutiveMatrix(2,2));

    FD[1].resize(4);
    FD[1][0]=(rVariables.F(1,1)*rVariables.ConstitutiveMatrix(1,0)+rVariables.F(1,0)*rVariables.ConstitutiveMatrix(2,0));
    FD[1][1]=(rVariables.F(1,1)*rVariables.ConstitutiveMatrix(1,1)+rVariables.F(1,0)*rVariables.ConstitutiveMatrix(2,1));
    FD[1][2]=(rVariables.F(1,1)*rVariables.ConstitutiveMatrix(1,2)+rVariables.F(1,0)*rVariables.ConstitutiveMatrix(2,2));
    FD[1][3]=(rVariables.F(1,1)*rVariables.ConstitutiveMatrix(1,2)+rVariables.F(1,0)*rVariables.ConstitutiveMatrix(2,2));

    FD[2].resize(4);
    FD[2][0]=(rVariables.F(0,1)*rVariables.ConstitutiveMatrix(1,0)+rVariables.F(0,0)*rVariables.ConstitutiveMatrix(2,0));
    FD[2][1]=(rVariables.F(0,1)*rVariables.ConstitutiveMatrix(1,1)+rVariables.F(0,0)*rVariables.ConstitutiveMatrix(2,1));
    FD[2][2]=(rVariables.F(0,1)*rVariables.ConstitutiveMatrix(1,2)+rVariables.F(0,0)*rVariables.ConstitutiveMatrix(2,2));
    FD[2][3]=(rVariables.F(0,1)*rVariables.ConstitutiveMatrix(1,2)+rVariables.F(0,0)*rVariables.ConstitutiveMatrix(2,2));

    FD[3].resize(4);
    FD[3][0]=(rVariables.F(1,0)*rVariables.ConstitutiveMatrix(0,0)+rVariables.F(1,1)*rVariables.ConstitutiveMatrix(2,0));
    FD[3][1]=(rVariables.F(1,0)*rVariables.ConstitutiveMatrix(0,1)+rVariables.F(1,1)*rVariables.ConstitutiveMatrix(2,1));
    FD[3][2]=(rVariables.F(1,0)*rVariables.ConstitutiveMatrix(0,2)+rVariables.F(1,1)*rVariables.ConstitutiveMatrix(2,2));
    FD[3][3]=(rVariables.F(1,0)*rVariables.ConstitutiveMatrix(0,2)+rVariables.F(1,1)*rVariables.ConstitutiveMatrix(2,2));

    std::vector<Vector> FDB(4);

    FDB[0]=ZeroVector(2);
    FDB[1]=ZeroVector(2);
    FDB[2]=ZeroVector(2);
    FDB[3]=ZeroVector(2);

    for(int i=0; i<4; i++)
    {
        for(int j=0; j<2; j++)
        {
            FDB[i][j]=FD[i][0]*rVariables.F(j,0)*rVariables.DN_DX(ndi,0)+FD[i][1]*rVariables.F(j,1)*rVariables.DN_DX(ndi,1)+
                      (FD[i][2]+FD[i][3])*(0.5)*(rVariables.F(j,0)*rVariables.DN_DX(ndi,1)+rVariables.F(j,1)*rVariables.DN_DX(ndi,0));

        }

    }

    // std::cout<<" FBD [0] "<<FDB[0]<<std::endl;
    // std::cout<<" FBD [1] "<<FDB[1]<<std::endl;
    // std::cout<<" FBD [2] "<<FDB[2]<<std::endl;
    // std::cout<<" FBD [3] "<<FDB[3]<<std::endl;

    for(int i=0; i<2; i++)
    {
        SigmaP[ndj][i]+=DirVector[0]*mContactVariables.ReferenceSurface.Normal[0]*(FDB[0][i])+
                        DirVector[1]*mContactVariables.ReferenceSurface.Normal[1]*(FDB[1][i])+
                        DirVector[0]*mContactVariables.ReferenceSurface.Normal[1]*(FDB[2][i])+
                        DirVector[1]*mContactVariables.ReferenceSurface.Normal[0]*(FDB[3][i]);
    }



}

//***********************************************************************************
//************************************************************************************

double& ContactDomainLM2DCondition::CalculateIntegrationWeight(double& rIntegrationWeight)
{
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if ( dimension == 2 ){
      ElementType& MasterElement = mContactVariables.GetMasterElement();
      if ( MasterElement.GetProperties().Has(THICKNESS) )
	  rIntegrationWeight *= MasterElement.GetProperties()[THICKNESS];
    }

    return rIntegrationWeight;
}


//************************************************************************************
//************************************************************************************

void ContactDomainLM2DCondition::CalculateNormalForce (double &F,ConditionVariables& rVariables,unsigned int& ndi,unsigned int& idir)
{
    F = rVariables.Contact.Multiplier.Normal*rVariables.Contact.dN_dn[ndi]*rVariables.Contact.CurrentSurface.Normal[idir];

}

//************************************************************************************
//************************************************************************************


void ContactDomainLM2DCondition::CalculateTangentStickForce (double &F,ConditionVariables& rVariables,unsigned int& ndi,unsigned int& idir)
{

	if( rVariables.Contact.Options.Is(ContactDomainUtilities::COMPUTE_FRICTION_FORCES) )
	{
		F=rVariables.Contact.Multiplier.Tangent*(rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]+rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Tangent[idir]);
	}
	else{
		F=0.0;
	}

	//std::cout<<" Multiplier Tangent "<<rVariables.Contact.Multiplier.Tangent<<" GapNormal "<<rVariables.Contact.CurrentGap.Normal<<" dN_dt "<<rVariables.Contact.dN_dt<<" Normal "<<rVariables.Contact.CurrentSurface.Normal<<" dN_drn "<<rVariables.Contact.dN_drn<<std::endl;

}

//************************************************************************************
//************************************************************************************


void ContactDomainLM2DCondition::CalculateTangentSlipForce (double &F,ConditionVariables& rVariables,unsigned int& ndi,unsigned int& idir)
{

	if( rVariables.Contact.Options.Is(ContactDomainUtilities::COMPUTE_FRICTION_FORCES) )
	{
		F=rVariables.Contact.Multiplier.Normal*(rVariables.Contact.FrictionCoefficient*rVariables.Contact.TangentialGapSign)*(rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]+rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Tangent[idir]);

	}
	else{
		F=0.0;
	}



}


//************************************************************************************
//************************************************************************************

void ContactDomainLM2DCondition::CalculateContactStiffness (double &Kcont,ConditionVariables& rVariables,unsigned int& ndi,unsigned int& ndj,unsigned int& idir,unsigned int& jdir)
{
    Kcont=0;

    //Normal contact contribution:
    //KI:
    Kcont = (0.5/rVariables.Contact.ContactFactor.Normal) * (rVariables.Contact.dN_dn[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dn[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir] )
      - rVariables.Contact.Multiplier.Normal*(rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]
					      + rVariables.Contact.dN_dn[ndi]*rVariables.Contact.CurrentSurface.Tangent[idir]*rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]
					      + rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dn[ndj]*rVariables.Contact.CurrentSurface.Tangent[jdir])
      - rVariables.Contact.CurrentTensil.Tangent*(rVariables.Contact.dN_dn[ndi]*rVariables.Contact.CurrentSurface.Normal[idir])*rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir];

    // std::cout<<" ndi "<<ndi<<" ndj "<<ndj<<" i "<<idir<<" j "<<jdir<<std::endl;
    // std::cout<<" Kg "<<Kcont;
    // //std::cout<<" constant "<<constant<<" Kg "<<Kcont*constant;
    // double K1=Kcont;

    //KII:
    Kcont+= rVariables.Contact.dN_dn[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.Nsigma[ndj][jdir];


    if(rVariables.Contact.Options.Is(ContactDomainUtilities::COMPUTE_FRICTION_STIFFNESS))
      {

	//Stick contact contribution:
	if(rVariables.Contact.Options.Is(NOT_SLIP))
	  {
	    //std::cout<<" + stick ";
	    //std::cout<<"(mu_on)";
            //KI:
            Kcont+= (0.5/rVariables.Contact.ContactFactor.Normal) * (rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]+rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Tangent[idir])*(rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+rVariables.Contact.dN_drn[ndj]*rVariables.Contact.CurrentSurface.Tangent[jdir])
	      + rVariables.Contact.Multiplier.Tangent * (rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dn[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]
							 + rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]
							 - rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Tangent[idir]*rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]
							 - rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Tangent[jdir])
	      + rVariables.Contact.CurrentTensil.Normal*(rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]+rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Tangent[idir])*(rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]);

            //KII:
            Kcont+= (rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]+rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Tangent[idir])*rVariables.Contact.Tsigma[ndj][jdir];

            //(he_a*ddN_dt_a(in)*ne_a(idime) + hddN_dn_n(in)*te_a(idime))*raux
	  }
	else{
        //Slip contact contribution:
	//std::cout<<" + slip ";
	    //std::cout<<"(mu_on)";
            //KI:
	    Kcont+= (rVariables.Contact.FrictionCoefficient*rVariables.Contact.TangentialGapSign)*
	      ((0.5/rVariables.Contact.ContactFactor.Normal)*(rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]+rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Tangent[idir])*(rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dn[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir])
		    +rVariables.Contact.Multiplier.Normal*(rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dn[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]
							 + rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]
							 - rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Tangent[idir]*rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]
							 - rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Tangent[jdir])
	            -rVariables.Contact.CurrentTensil.Tangent*(rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]+rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Tangent[idir])*(rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]));

            //KII:
            Kcont+= (rVariables.Contact.FrictionCoefficient*rVariables.Contact.TangentialGapSign)*(rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]+rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Tangent[idir])*rVariables.Contact.Tsigma[ndj][jdir];

        }
      }

    //std::cout<<" Ks "<<Kcont-K1<<std::endl;

}


//************************************************************************************
//************************************************************************************


ContactDomainUtilities::PointType & ContactDomainLM2DCondition::CalculateCurrentTangent ( PointType &rTangent )
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


inline bool ContactDomainLM2DCondition::CheckFictiousContacts(ConditionVariables& rVariables)
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


void ContactDomainLM2DCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ContactDomainCondition )
    // std::cout<<" Restart Save MASTER_ELEMENTS "<<GetValue(MASTER_ELEMENTS).size()<<" MASTER_NODES "<<GetValue(MASTER_NODES).size()<<std::endl;
}

void ContactDomainLM2DCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ContactDomainCondition )
    // std::cout<<" Restart Load MASTER_ELEMENTS "<<GetValue(MASTER_ELEMENTS).size()<<" MASTER_NODES "<<GetValue(MASTER_NODES).size()<<std::endl;
}



} // Namespace Kratos
