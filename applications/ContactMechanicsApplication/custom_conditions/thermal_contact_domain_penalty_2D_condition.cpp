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
#include "includes/kratos_flags.h"
#include "custom_conditions/thermal_contact_domain_penalty_2D_condition.hpp"

#include "contact_mechanics_application_variables.h"


namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ThermalContactDomainPenalty2DCondition::ThermalContactDomainPenalty2DCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : ThermalContactDomainCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}


//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ThermalContactDomainPenalty2DCondition::ThermalContactDomainPenalty2DCondition( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
    : ThermalContactDomainCondition( NewId, pGeometry, pProperties )
{
    mThisIntegrationMethod = GetGeometry().GetDefaultIntegrationMethod();
}


//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

ThermalContactDomainPenalty2DCondition::ThermalContactDomainPenalty2DCondition( ThermalContactDomainPenalty2DCondition const& rOther)
    :ThermalContactDomainCondition(rOther)
{
}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

ThermalContactDomainPenalty2DCondition&  ThermalContactDomainPenalty2DCondition::operator=(ThermalContactDomainPenalty2DCondition const& rOther)
{
    ThermalContactDomainCondition::operator=(rOther);

    return *this;
}


//*********************************OPERATIONS*****************************************
//************************************************************************************

Condition::Pointer ThermalContactDomainPenalty2DCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
  return Kratos::make_shared<ThermalContactDomainPenalty2DCondition>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}


//*******************************DESTRUCTOR*******************************************
//************************************************************************************


ThermalContactDomainPenalty2DCondition::~ThermalContactDomainPenalty2DCondition()
{
}




//************* COMPUTING  METHODS
//************************************************************************************
//************************************************************************************


void ThermalContactDomainPenalty2DCondition::SetMasterGeometry()

{
    unsigned int vsize=GetValue(MASTER_ELEMENTS).size();
    Element::ElementType& MasterElement = GetValue(MASTER_ELEMENTS)[vsize-1];
    mContactVariables.SetMasterElement(MasterElement);


    vsize=GetValue(MASTER_NODES).size();
    Element::NodeType&    MasterNode   = GetValue(MASTER_NODES)[vsize-1];
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

        NodesArrayType vertex;
	mContactVariables.order.resize(GetGeometry().PointsNumber(),false);

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
        }



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


	// for(unsigned int i=0; i<MasterElement.GetGeometry().size(); i++)
        // {
	//   vertex.push_back(MasterElement.GetGeometry()[i]);
	// }
	// mpMasterGeometry= GetGeometry().Create(vertex);

	mContactVariables.SetMasterGeometry( MasterElement.GetGeometry() );

    }
    else
    {
        KRATOS_THROW_ERROR( std::invalid_argument, "MASTERNODE do not belongs to MASTER ELEMENT", "" );

    }


}

//*********************************COMPUTE KINEMATICS*********************************
//************************************************************************************


void ThermalContactDomainPenalty2DCondition::CalculateKinematics(GeneralVariables& rVariables,
						    ProcessInfo& rCurrentProcessInfo,
						    const unsigned int& rPointNumber)
{
    KRATOS_TRY

    //Calculate Current Contact Projections
    this->CalcProjections(rVariables,rCurrentProcessInfo);

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void ThermalContactDomainPenalty2DCondition::CalcProjections(GeneralVariables & rVariables, ProcessInfo& rCurrentProcessInfo)
{

    //Contact face segment node1-node2
    unsigned int node1=mContactVariables.nodes[0];
    unsigned int node2=mContactVariables.nodes[1];
    unsigned int slave=mContactVariables.slaves.back();

    //a.- Compute the Reference Normal and Tangent

    //Get Reference Normal
    //rVariables.ReferenceSurface.Normal=GetValue(NORMAL);

    PointType PP1  =  GetGeometry()[node1].Coordinates() - (GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT) -GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1) );
    PointType PP2  =  GetGeometry()[node2].Coordinates() - (GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1) );

    rVariables.ReferenceSurface.Normal=mContactUtilities.CalculateFaceNormal(rVariables.ReferenceSurface.Normal,PP1,PP2);


    //Set Reference Tangent
    rVariables.ReferenceSurface.Tangent=mContactUtilities.CalculateFaceTangent(rVariables.ReferenceSurface.Tangent,rVariables.ReferenceSurface.Normal);


    //b.- Compute the Current Normal and Tangent

    PointType PS  =  GetGeometry()[slave].Coordinates();
    PointType P1  =  GetGeometry()[node1].Coordinates();
    PointType P2  =  GetGeometry()[node2].Coordinates();

    //compute the current normal vector
    rVariables.CurrentSurface.Normal=mContactUtilities.CalculateFaceNormal(rVariables.CurrentSurface.Normal,P1,P2);


    // if(double(inner_prod(rVariables.CurrentSurface.Normal,rVariables.ReferenceSurface.Normal))<0) //to give the correct direction
    //     rVariables.CurrentSurface.Normal*=-1;

    // rVariables.CurrentSurface.Normal /= norm_2(rVariables.CurrentSurface.Normal);  //to be unitary

    // if(!(norm_2(rVariables.CurrentSurface.Normal)))
    //     rVariables.CurrentSurface.Normal=rVariables.ReferenceSurface.Normal;


    //compute the current tangent vector
    //rVariables.CurrentSurface.Tangent=mContactUtilities.CalculateFaceTangent(rVariables.CurrentSurface.Tangent,P1,P2);
    rVariables.CurrentSurface.Tangent=mContactUtilities.CalculateFaceTangent(rVariables.CurrentSurface.Tangent,rVariables.CurrentSurface.Normal);

    // std::cout<<" current face  normal  "<<rVariables.CurrentSurface.Normal<<std::endl;
    // std::cout<<" current face  tangent  "<<rVariables.CurrentSurface.Tangent<<std::endl;

    //Current normal:   rVariables.ReferenceSurface.Normal


    //c.- Compute A_n,B_n,L_n
    rVariables.ReferenceBase.resize(1);
    rVariables.CurrentBase.resize(1);

    //a, b, l:
    mContactUtilities.CalculateBaseDistances (rVariables.CurrentBase[0],P1,P2,PS,rVariables.CurrentSurface.Normal);

    //Write Current Positions:
    // std::cout<<" Current position node 1 "<<P1<<std::endl;
    // std::cout<<" Current position node 2 "<<P2<<std::endl;
    // std::cout<<" Current position node s "<<PS<<std::endl;

    //Set Projection Vector:
    rVariables.ProjectionsVector = ZeroVector(3);
    rVariables.ProjectionsVector[slave] = 1;
    rVariables.ProjectionsVector[node1] = -(rVariables.CurrentBase[0].A/rVariables.CurrentBase[0].L);
    rVariables.ProjectionsVector[node2] = -(rVariables.CurrentBase[0].B/rVariables.CurrentBase[0].L);


    //A, B, L:

    PS =  GetGeometry()[slave].Coordinates() - ( GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[slave].FastGetSolutionStepValue(DISPLACEMENT,1) );
    P1 =  GetGeometry()[node1].Coordinates() - ( GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node1].FastGetSolutionStepValue(DISPLACEMENT,1) );
    P2 =  GetGeometry()[node2].Coordinates() - ( GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT) - GetGeometry()[node2].FastGetSolutionStepValue(DISPLACEMENT,1) );

    mContactUtilities.CalculateBaseDistances (rVariables.ReferenceBase[0],P1,P2,PS,rVariables.ReferenceSurface.Normal);

    //d.-obtain the (thermal_gap)

    double& TS  =  GetGeometry()[slave].FastGetSolutionStepValue(TEMPERATURE);
    double& T1  =  GetGeometry()[node1].FastGetSolutionStepValue(TEMPERATURE);
    double& T2  =  GetGeometry()[node2].FastGetSolutionStepValue(TEMPERATURE);

    rVariables.ThermalGap  = TS * rVariables.ProjectionsVector[slave];
    rVariables.ThermalGap += T1 * rVariables.ProjectionsVector[node1];
    rVariables.ThermalGap += T2 * rVariables.ProjectionsVector[node2];


    //CHECK IF THE ELEMENT IS ACTIVE and THERE IS FRICTION:

    //Check if is active checking the mechanical force component
    rVariables.Options.Set(ACTIVE,false);


    const unsigned int number_of_nodes = GetGeometry().PointsNumber();


    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
	PointType & ContactForce = GetGeometry()[i].FastGetSolutionStepValue(CONTACT_FORCE);
	if(norm_2(ContactForce)>0){
	    rVariables.Options.Set(ACTIVE,true);
	    break;
	}
    }


    //Check if is stick or slip checking RelativeTangent Velocity and the Friction (Tangent) force component
    rVariables.Options.Set(SLIP,false); //impose stick

    PointType TangentVelocity (3,0.0);
    CalculateRelativeVelocity(rVariables, TangentVelocity, rCurrentProcessInfo);

    rVariables.RelativeVelocityNorm = norm_2(TangentVelocity);

    PointType & ContactForce = GetGeometry()[slave].FastGetSolutionStepValue(CONTACT_FORCE);

    PointType TangentForce = ContactForce;
    TangentForce -= (inner_prod(ContactForce,rVariables.CurrentSurface.Normal)) * ContactForce;

    rVariables.FrictionForceNorm    = norm_2(TangentForce);


    if(rVariables.RelativeVelocityNorm>0 && rVariables.FrictionForceNorm>0){
        rVariables.Options.Set(SLIP,true);
    }
}


//***********************************************************************************
//************************************************************************************

double& ThermalContactDomainPenalty2DCondition::CalculateIntegrationWeight(double& rIntegrationWeight)
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


void ThermalContactDomainPenalty2DCondition::CalculateThermalConductionForce (double &F, GeneralVariables& rVariables, unsigned int& ndi)
{
	F = mContactVariables.StabilizationFactor * rVariables.ThermalGap * rVariables.ProjectionsVector[ndi];
	//std::cout<<" Thermal Conduction Force [stab: "<<mContactVariables.StabilizationFactor<<" therma_gap: "<<rVariables.ThermalGap<<" F "<<F<<"]"<<std::endl;

}

//************************************************************************************
//************************************************************************************


void ThermalContactDomainPenalty2DCondition::CalculateThermalFrictionForce (double &F, GeneralVariables& rVariables, unsigned int& ndi)
{


    PointType & ContactForce = GetGeometry()[ndi].FastGetSolutionStepValue(CONTACT_FORCE);

    PointType TangentForce = ContactForce;
    TangentForce -= (inner_prod(ContactForce,rVariables.CurrentSurface.Normal)) * ContactForce;

    double FrictionForceNorm = norm_2(TangentForce);

    //this must be introduced as a property
    double HeatWorkFraction = 0.9;


    F=0.0;
    if( rVariables.Options.Is(SLIP) )
    {
      F = rVariables.CurrentSurface.Tangent[ndi] * FrictionForceNorm * rVariables.RelativeVelocityNorm * 0.5 * HeatWorkFraction;
    }


    //std::cout<<" Thermal Friction Force [friction_norm: "<<FrictionForceNorm<<" velocity_norm: "<<rVariables.RelativeVelocityNorm<<" F "<<F<<"]"<<std::endl;

}


//************************************************************************************
//************************************************************************************


ContactDomainUtilities::PointType & ThermalContactDomainPenalty2DCondition::CalculateCurrentTangent ( PointType &rTangent )
{

	unsigned int node1=mContactVariables.nodes[0];
	unsigned int node2=mContactVariables.nodes[1];

	PointType P1  =  GetGeometry()[node1].Coordinates();
	PointType P2  =  GetGeometry()[node2].Coordinates();

	//Set Reference Tangent
	rTangent = mContactUtilities.CalculateFaceTangent(rTangent,P1,P2);

	return rTangent;

}

//************************************************************************************
//************************************************************************************
/**
 * This function provides the place to perform checks on the completeness of the input.
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 */
int  ThermalContactDomainPenalty2DCondition::Check( const ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    return 0;

    KRATOS_CATCH( "" );
}


void ThermalContactDomainPenalty2DCondition::save( Serializer& rSerializer ) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ThermalContactDomainCondition );

}

void ThermalContactDomainPenalty2DCondition::load( Serializer& rSerializer )
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ThermalContactDomainCondition );

}



} // Namespace Kratos
