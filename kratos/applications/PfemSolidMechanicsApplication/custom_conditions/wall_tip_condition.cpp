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
#include "pfem_solid_mechanics_application.h"

namespace Kratos
{



//************************************************************************************
//************************************************************************************
  
  KRATOS_CREATE_LOCAL_FLAG( WallTipCondition, WALL_TIP, 6 );


//************************************************************************************
//************************************************************************************
WallTipCondition::WallTipCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
  //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
WallTipCondition::WallTipCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )
{

}

Condition::Pointer WallTipCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new WallTipCondition( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

WallTipCondition::~WallTipCondition()
{
}

//************************************************************************************
//************************************************************************************

void WallTipCondition::Initialize()
{
    KRATOS_TRY
      mWallTip.wall_tip_set =false;     
    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void WallTipCondition::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY

      if(!mWallTip.wall_tip_set){
  
	double & TipRadius                   = GetGeometry()[0].FastGetSolutionStepValue( WALL_TIP_RADIUS );
	array_1d<double,3> & TipPosition     = GetGeometry()[0].FastGetSolutionStepValue( WALL_REFERENCE_POINT );
	array_1d<double,3> & TipVelocity     = GetGeometry()[0].FastGetSolutionStepValue( WALL_VELOCITY );
                 
      
	mWallTip.Radius   = TipRadius;    //Constant
	mWallTip.Center   = TipPosition;  //Must be updated
	mWallTip.Velocity = TipVelocity;

	mWallTip.Position = TipPosition;

	// std::cout<<" InitializeStep [TipCenter:"<<mWallTip.Position<<", TipVelc: "<<mWallTip.Velocity<<", TipRadius: "<<TipRadius<<"] (ID: "<<GetGeometry()[0].Id()<<"; "<<GetGeometry()[1].Id()<<") "<<std::endl;
     
	if(TipRadius!=0)
	  mWallTip.wall_tip_set = true;
      }

  
    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void WallTipCondition::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
{

    KRATOS_TRY

      if(mWallTip.wall_tip_set){

	mWallTip.Position = mWallTip.Center + mWallTip.Velocity * CurrentProcessInfo[DELTA_TIME];

	array_1d<double,3> & TipPosition1    = GetGeometry()[0].FastGetSolutionStepValue( WALL_REFERENCE_POINT );	
	array_1d<double,3> & TipPosition2    = GetGeometry()[1].FastGetSolutionStepValue( WALL_REFERENCE_POINT );
	
	TipPosition1 = mWallTip.Position;
	TipPosition2 = mWallTip.Position;

	mWallTip.wall_tip_set = false;
	
	// std::cout<<" FinalizeStep [TipCenter:"<<mWallTip.Position<<", TipVelc: "<<mWallTip.Velocity<<", TipRadius: "<<mWallTip.Radius<<"] (ID: "<<GetGeometry()[0].Id()<<"; "<<GetGeometry()[1].Id()<<") "<<std::endl;
	
      }
      // else{

      // 	std::cout<<" Tool NOT SET "<<std::endl;
      // }


    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************



} // Namespace Kratos


