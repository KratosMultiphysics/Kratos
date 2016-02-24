//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes
#include "wall_condition.hpp"

// Project includes
#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
WallCondition::WallCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
  //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
WallCondition::WallCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )
{

}


//************************************CLONE*******************************************
//************************************************************************************


Condition::Pointer WallCondition::Clone( IndexType NewId, NodesArrayType const& ThisNodes ) const
{
  return this->Create(NewId, ThisNodes, pGetProperties());
}


//************************************************************************************
//************************************************************************************

Condition::Pointer WallCondition::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new WallCondition( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}


WallCondition::~WallCondition()
{
}

//************************************************************************************
//************************************************************************************

void WallCondition::Initialize()
{
    KRATOS_TRY

      mWallVariables.wall_tip_set =false;     

    KRATOS_CATCH( "" )
}



//************************************************************************************
//************************************************************************************

void WallCondition::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
{
    KRATOS_TRY

      // if(!mWallVariables.wall_tip_set){
  
      // 	double & TipRadius                   = GetGeometry()[0].FastGetSolutionStepValue( WALL_TIP_RADIUS );
      // 	array_1d<double,3> & TipPosition     = GetGeometry()[0].FastGetSolutionStepValue( WALL_REFERENCE_POINT );
      // 	array_1d<double,3> & TipVelocity     = GetGeometry()[0].FastGetSolutionStepValue( WALL_VELOCITY );
                 
      
      // 	mWallVariables.Radius   = TipRadius;    //Constant
      // 	mWallVariables.Center   = TipPosition;  //Must be updated
      // 	mWallVariables.Velocity = TipVelocity;

      // 	mWallVariables.Position = TipPosition;

      // 	// std::cout<<" InitializeStep [TipCenter:"<<mWallVariables.Position<<", TipVelc: "<<mWallVariables.Velocity<<", TipRadius: "<<TipRadius<<"] (ID: "<<GetGeometry()[0].Id()<<"; "<<GetGeometry()[1].Id()<<") "<<std::endl;
     
      // 	if(TipRadius!=0)
      // 	  mWallVariables.wall_tip_set = true;
      // }

  
    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

void WallCondition::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
{

    KRATOS_TRY

      // if(mWallVariables.wall_tip_set){

      // 	mWallVariables.Position = mWallVariables.Center + mWallVariables.Velocity * CurrentProcessInfo[DELTA_TIME];

      // 	array_1d<double,3> & TipPosition1    = GetGeometry()[0].FastGetSolutionStepValue( WALL_REFERENCE_POINT );	
      // 	array_1d<double,3> & TipPosition2    = GetGeometry()[1].FastGetSolutionStepValue( WALL_REFERENCE_POINT );
	
      // 	TipPosition1 = mWallVariables.Position;
      // 	TipPosition2 = mWallVariables.Position;

      // 	mWallVariables.wall_tip_set = false;
	
      // 	// std::cout<<" FinalizeStep [TipCenter:"<<mWallVariables.Position<<", TipVelc: "<<mWallVariables.Velocity<<", TipRadius: "<<mWallVariables.Radius<<"] (ID: "<<GetGeometry()[0].Id()<<"; "<<GetGeometry()[1].Id()<<") "<<std::endl;
	
      // }
      // // else{

      // // 	std::cout<<" Tool NOT SET "<<std::endl;
      // // }


    KRATOS_CATCH( "" )

}

//************************************************************************************
//************************************************************************************



} // Namespace Kratos


