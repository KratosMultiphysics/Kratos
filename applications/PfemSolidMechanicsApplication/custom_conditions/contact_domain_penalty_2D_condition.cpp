//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "includes/kratos_flags.h"
#include "custom_conditions/contact_domain_penalty_2D_condition.hpp"

#include "pfem_solid_mechanics_application.h"


namespace Kratos
{


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ContactDomainPenalty2DCondition::ContactDomainPenalty2DCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : ContactDomainLM2DCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ContactDomainPenalty2DCondition::ContactDomainPenalty2DCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : ContactDomainLM2DCondition( NewId, pGeometry, pProperties )
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

ContactDomainPenalty2DCondition::ContactDomainPenalty2DCondition( ContactDomainPenalty2DCondition const& rOther)
    :ContactDomainLM2DCondition(rOther)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

ContactDomainPenalty2DCondition&  ContactDomainPenalty2DCondition::operator=(ContactDomainPenalty2DCondition const& rOther)
{
    ContactDomainLM2DCondition::operator=(rOther);

    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Condition::Pointer ContactDomainPenalty2DCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer(new ContactDomainPenalty2DCondition( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************


ContactDomainPenalty2DCondition::~ContactDomainPenalty2DCondition()
{
}


//************* STARTING - ENDING  METHODS
//************************************************************************************
//************************************************************************************

//*****************************COMPUTE PENALTY FACTOR*********************************
//************************************************************************************

void ContactDomainPenalty2DCondition::CalculateContactFactor( ProcessInfo& rCurrentProcessInfo )
{
    //Initilialize Tau for the stabilization
    double penalty_parameter = 1000;
    penalty_parameter = GetProperties()[PENALTY_PARAMETER];

    ElementType& MasterElement = mContactVariables.GetMasterElement();
  
    //Look at the nodes, get the slave and get the Emin

    //Contact face segment node1-node2
    unsigned int slave = mContactVariables.slaves.back();


    double Eslave = GetGeometry()[slave].GetValue(NEIGHBOUR_ELEMENTS)[0].GetProperties()[YOUNG_MODULUS];
    double Emax   = MasterElement.GetProperties()[YOUNG_MODULUS];

    //STANDARD OPTION
    if(Emax<Eslave)
	Emax=Eslave;


    mContactVariables.PenaltyFactor = 0.5 * penalty_parameter ;

    // mContactVariables.PenaltyParameter = 0.5 / mContactVariables.StabilizationFactor ;
   
}



//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


//****************************CALCULATE EXPLICIT PENALTY PARAMETERS*******************
//************************************************************************************

void ContactDomainPenalty2DCondition::CalculateExplicitFactors(GeneralVariables& rVariables, ProcessInfo& rCurrentProcessInfo)
{


    //Contact face segment node1-node2
    unsigned int node1=mContactVariables.nodes[0];
    unsigned int node2=mContactVariables.nodes[1];    
    unsigned int slave=mContactVariables.slaves.back();

    // std::cout<<" ************ CONTACT ELEMENT "<<this->Id()<<" ************* "<<std::endl;
    // std::cout<<std::endl;
    
    // unsigned int vsize=GetValue(MASTER_ELEMENTS).size();
    // Element::ElementType& MasterElement = GetValue(MASTER_ELEMENTS)[vsize-1];
    
    // std::cout<<" master element "<<MasterElement.Id()<<std::endl;
    // std::cout<<" Elastic Modulus "<<MasterElement.GetProperties()[YOUNG_MODULUS]<<std::endl;

    // std::cout<<" Nodes 1,2,S "<<node1<<" "<<node2<<" "<<slave<<std::endl;

    // std::cout<<" Nodes ID  "<<GetGeometry()[0].Id()<<" "<<GetGeometry()[1].Id()<<" "<<GetGeometry()[2].Id()<<std::endl;

    // std::cout<<" Master Nodes ID  "<<(*mpMasterGeometry)[0].Id()<<" "<<(*mpMasterGeometry)[1].Id()<<" "<<(*mpMasterGeometry)[2].Id()<<std::endl;


    //1.- Compute the Current Normal and Tangent

    VectorType PS  =  GetGeometry()[slave].Coordinates() + ( GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,1) );
    VectorType P1  =  GetGeometry()[node1].Coordinates() + ( GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1) );
    VectorType P2  =  GetGeometry()[node2].Coordinates() + ( GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1) );

    //compute the current normal vector
    rVariables.Contact.CurrentSurface.Normal=ComputeFaceNormal(rVariables.Contact.CurrentSurface.Normal,P1,P2);

    //std::cout<<" Current Normal "<<rVariables.Contact.CurrentSurface.Normal<<std::endl;

    if(double(inner_prod(rVariables.Contact.CurrentSurface.Normal,mContactVariables.ReferenceSurface.Normal))<0) //to give the correct direction
        rVariables.Contact.CurrentSurface.Normal*=-1;

    rVariables.Contact.CurrentSurface.Normal /= norm_2(rVariables.Contact.CurrentSurface.Normal);  //to be unitary

    if(!(norm_2(rVariables.Contact.CurrentSurface.Normal)))
        rVariables.Contact.CurrentSurface.Normal=mContactVariables.ReferenceSurface.Normal;


    //compute the current tangent vector
    //rVariables.Contact.CurrentSurface.Tangent=ComputeFaceTangent(rVariables.Contact.CurrentSurface.Tangent,P1,P2);
    rVariables.Contact.CurrentSurface.Tangent = this->ComputeFaceTangent(rVariables.Contact.CurrentSurface.Tangent, rVariables.Contact.CurrentSurface.Normal);



    //std::cout<<" reference face  normal  "<<mContactVariables.ReferenceSurface.Normal<<std::endl;
    //std::cout<<" reference face  tangent  "<<mContactVariables.ReferenceSurface.Tangent<<std::endl;


    //std::cout<<" current face  normal  "<<rVariables.Contact.CurrentSurface.Normal<<std::endl;
    //std::cout<<" current face  tangent  "<<rVariables.Contact.CurrentSurface.Tangent<<std::endl;

    //Current normal:   mContactVariables.ReferenceSurface.Normal

    //Reference normal: n_n,t_n  -> mContactVariables.ReferenceSurface.Normal / rVariables.Contact.Tangent
    //Current normal:   n,t      -> rVariables.Contact.CurrentSurface.Normal /  rVariables.Contact.CurrentSurface.Tangent

    //d.- Compute A_n,B_n,L_n
    rVariables.Contact.ReferenceBase.resize(1);
    rVariables.Contact.CurrentBase.resize(1);

    //a, b, l:
    CalcBaseDistances (rVariables.Contact.CurrentBase[0],P1,P2,PS,rVariables.Contact.CurrentSurface.Normal);

    //Write Current Positions:
    // std::cout<<" Current position node 1 "<<P1<<std::endl;
    // std::cout<<" Current position node 2 "<<P2<<std::endl;
    // std::cout<<" Current position node s "<<PS<<std::endl;


    //A, B, L:

    PS =  GetGeometry()[slave].Coordinates();
    P1 =  GetGeometry()[node1].Coordinates();
    P2 =  GetGeometry()[node2].Coordinates();

    CalcBaseDistances (rVariables.Contact.ReferenceBase[0],P1,P2,PS,mContactVariables.ReferenceSurface.Normal);


    //complete the computation of the stabilization gap
    //rVariables.Contact.ContactFactor = mContactVariables.StabilizationFactor * rVariables.Contact.ReferenceBase[0].L;
    rVariables.Contact.ContactFactor = mContactVariables.PenaltyFactor / rVariables.Contact.ReferenceBase[0].L;

    //e.-obtain the (g_N)3 and (g_T)3 for the n configuration
    //Write Current Positions:
    // std::cout<<" Reference position node 1 "<<P1<<std::endl;
    // std::cout<<" Reference position node 2 "<<P2<<std::endl;
    // std::cout<<" Reference position node s "<<PS<<std::endl;

    double ReferenceGapN = inner_prod((PS - P1),mContactVariables.ReferenceSurface.Normal);
    //std::cout<<" Reference GAP "<<ReferenceGapN<<std::endl;
    
    double ReferenceGapT = ReferenceGapN;

    double H = ReferenceGapN;

    VectorType DS  =  GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT)-GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,1);
    VectorType D1  =  GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT)-GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1);
    VectorType D2  =  GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT)-GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1);

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

 
    //Write Displacements:
    // std::cout<<" displacement node 1 "<<D1<<std::endl;
    // std::cout<<" displacement node 2 "<<D2<<std::endl;
    // std::cout<<" displacement node s "<<DS<<std::endl;

    // std::cout<<" L :"<<rVariables.Contact.ReferenceBase[0].L<<" A :"<<rVariables.Contact.ReferenceBase[0].A<<" B :"<<rVariables.Contact.ReferenceBase[0].B<<std::endl;
    // std::cout<<" gN3 ref "<<ReferenceGapN<<std::endl;


     // std::cout<<" current   normal  "<<rVariables.Contact.CurrentSurface.Normal<<std::endl;
     // std::cout<<" reference normal  "<<mContactVariables.ReferenceSurface.Normal<<std::endl;


    rVariables.Contact.CurrentGap.Normal=ReferenceGapN; //(g_N)3 -- needed in the Kcont1 computation
    rVariables.Contact.CurrentGap.Tangent=ReferenceGapT; //(g_T)3 -- needed in the Kcont1 computation

    //f.- get total effective gap as: gap_n^eff=gap_n+(PreviousTimeStep/CurrentTimeStep)*(gap_n-gap_n-1)

    //gap_n   (in function of the n position of the other node)

    rVariables.Contact.TangentialGapSign=1;

    if((rVariables.Contact.CurrentGap.Tangent)<0)
    {
        rVariables.Contact.TangentialGapSign*=(-1);
    }

    if(H==0) rVariables.Contact.TangentialGapSign=0;

    //5.- Compute (Lagrange) Multipliers


    // double CurrentTimeStep  = rCurrentProcessInfo[DELTA_TIME];
    // double PreviousTimeStep = rCurrentProcessInfo[PREVIOUS_DELTA_TIME];
        

    //rVariables.Contact.Penalty.Normal  = H * (0.5/rVariables.Contact.ContactFactor); //(1/(2*Tau)) is now KN, the penalty parameter
    rVariables.Contact.Penalty.Normal  = H * rVariables.Contact.ContactFactor; 
    rVariables.Contact.Penalty.Tangent = 0;
    //std::cout<<" PreTime "<<Time.PreStep<<" Time "<<Time.Step<<std::endl;

    //CHECK IF THE ELEMENT IS ACTIVE:

    rVariables.Contact.Options.Set(SLIP,false); 

    if(rVariables.Contact.Penalty.Normal<0)   //if(EffectiveGap<0){
    {

	 rVariables.Contact.Options.Set(ACTIVE,true); //normal contact active

        // if(fabs(EffectiveGapT)<=rVariables.Contact.FrictionCoefficient*fabs(EffectiveGapN))
        // {
        //     rVariables.Contact.Options.Set(SLIP,false);  //contact stick case active
	//     rCurrentProcessInfo[NUMBER_OF_STICK_CONTACTS] += 1;
        // }
        // else
        // {
        //     rVariables.Contact.Options.Set(SLIP,true);  //contact slip  case active
	//     rCurrentProcessInfo[NUMBER_OF_SLIP_CONTACTS] += 1;
        // }
	rVariables.Contact.Options.Set(COMPUTE_FRICTION_STIFFNESS); 
	rVariables.Contact.Options.Set(COMPUTE_FRICTION_FORCES);
	//rCurrentProcessInfo[NUMBER_OF_ACTIVE_CONTACTS] += 1;

    }
    else
    {
	rVariables.Contact.Options.Set(ACTIVE,false); //normal contact not active
    }


    rVariables.Contact.Options.Set(SLIP,false); //impose stick


    if(rVariables.Contact.Penalty.Tangent<0)  //add the sign of the Lagrange Multiplier
    {
        rVariables.Contact.TangentialGapSign*=(-1);
    }

    if(mContactVariables.IterationCounter < 1)
      mContactVariables.IterationCounter += 1;

    
    //now Multiplier.Normal is the penalty parameter:
    rVariables.Contact.Penalty.Normal = rVariables.Contact.ContactFactor;


}


//************************************************************************************
//************************************************************************************

void ContactDomainPenalty2DCondition::CalculateNormalForce (double &F,GeneralVariables& rVariables,unsigned int& ndi,unsigned int& idir)
{    
	F = rVariables.Contact.CurrentGap.Normal*rVariables.Contact.Penalty.Normal*rVariables.Contact.dN_dn[ndi]*rVariables.Contact.CurrentSurface.Normal[idir];

}

//************************************************************************************
//************************************************************************************


void ContactDomainPenalty2DCondition::CalculateTangentStickForce (double &F,GeneralVariables& rVariables,unsigned int& ndi,unsigned int& idir)
{

	if( rVariables.Contact.Options.Is(COMPUTE_FRICTION_FORCES) )
	{
		F=rVariables.Contact.Penalty.Tangent*(rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]+rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Tangent[idir]);
	}
	else{
		F=0.0;
	}

}

//************************************************************************************
//************************************************************************************


void ContactDomainPenalty2DCondition::CalculateTangentSlipForce (double &F,GeneralVariables& rVariables,unsigned int& ndi,unsigned int& idir)
{

	if( rVariables.Contact.Options.Is(COMPUTE_FRICTION_FORCES) )
	{
		F=rVariables.Contact.Penalty.Normal*(rVariables.Contact.FrictionCoefficient*rVariables.Contact.TangentialGapSign)*(rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]+rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Tangent[idir]);

	}
	else{
		F=0.0;
	}



}



//************************************************************************************
//************************************************************************************

void ContactDomainPenalty2DCondition::CalcContactStiffness (double &Kcont,GeneralVariables& rVariables,unsigned int& ndi,unsigned int& ndj,unsigned int& idir,unsigned int& jdir)
{
    Kcont=0;

    
    //Normal contact penalty contribution:
    //KI:
    double KN = rVariables.Contact.Multiplier.Normal;

    Kcont= KN*( (rVariables.Contact.dN_dn[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dn[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir])-
		(rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dn[ndj]*rVariables.Contact.CurrentSurface.Tangent[jdir])-
		(rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dn[ndi]*rVariables.Contact.CurrentSurface.Tangent[idir]*rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir])-
		((rVariables.Contact.CurrentGap.Normal*rVariables.Contact.CurrentGap.Normal)*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]));

    // std::cout<<" ndi "<<ndi<<" ndj "<<ndj<<" i "<<idir<<" j "<<jdir<<std::endl;
    // std::cout<<" Kg "<<Kcont;
    // //std::cout<<" constant "<<constant<<" Kg "<<Kcont*constant;
    // double K1=Kcont;


    //Stick contact contribution:
    // if(rVariables.Contact.Options.Is(NOT_SLIP))
    // {
    // 	//std::cout<<" + stick ";
    //     if(rVariables.Contact.Options.Is(COMPUTE_FRICTION_STIFFNESS))
    //     {
    // 	    //std::cout<<"(mu_on)";
    //         //KI:
    //         Kcont+= (0.5/(mContactVariables.StabilizationFactor*rVariables.Contact.ReferenceBase[0].L))*(rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]+rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Tangent[idir])*
    //                 (rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+rVariables.Contact.dN_drn[ndj]*rVariables.Contact.CurrentSurface.Tangent[jdir])+ rVariables.Contact.Multiplier.Tangent*(rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dn[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+ rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]- rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Tangent[idir]*rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]- rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Tangent[jdir])+ rVariables.Contact.CurrentTensil.Normal*(rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]+rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Tangent[idir])*(rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]);

    //         //KII:
    //         Kcont+= (rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]+rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Tangent[idir])*rVariables.Contact.Tsigma[ndj][jdir];

    //         //(he_a*ddN_dt_a(in)*ne_a(idime) + hddN_dn_n(in)*te_a(idime))*raux
    //     }
    // }
    // else
    // {
    //     //Slip contact contribution:
    // 	//std::cout<<" + slip ";
    //     if(rVariables.Contact.Options.Is(COMPUTE_FRICTION_STIFFNESS))
    //     {
    // 	    //std::cout<<"(mu_on)";
    //         //KI:
    //         Kcont+= (rVariables.Contact.FrictionCoefficient*rVariables.Contact.TangentialGapSign)*((0.5/(mContactVariables.StabilizationFactor*rVariables.Contact.ReferenceBase[0].L))*(rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]+rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Tangent[idir])* (rVariables.Contact.dN_dn[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir])+
    //                 rVariables.Contact.NormalMultiplier*(rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dn[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]+ rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]- rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Tangent[idir]*rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]- rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]*rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Tangent[jdir])- rVariables.Contact.CurrentTensil.Tangent*(rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]+rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Tangent[idir])*(rVariables.Contact.dN_dt[ndj]*rVariables.Contact.CurrentSurface.Normal[jdir]));


    //         //KII:
    //         Kcont+= (rVariables.Contact.FrictionCoefficient*rVariables.Contact.TangentialGapSign)*(rVariables.Contact.CurrentGap.Normal*rVariables.Contact.dN_dt[ndi]*rVariables.Contact.CurrentSurface.Normal[idir]+rVariables.Contact.dN_drn[ndi]*rVariables.Contact.CurrentSurface.Tangent[idir])*rVariables.Contact.Nsigma[ndj][jdir];

    //     }
    // }

    //std::cout<<" Ks "<<Kcont-K1<<std::endl;

}


void ContactDomainPenalty2DCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ContactDomainLM2DCondition );
}

void ContactDomainPenalty2DCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ContactDomainLM2DCondition );
}



} // Namespace Kratos


