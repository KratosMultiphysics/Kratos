//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes


// External includes


// Project includes
#include "custom_conditions/beam_point_pressure_condition.hpp"

#include "pfem_solid_mechanics_application.h"

namespace Kratos
{
  //************************************************************************************
  //************************************************************************************
  BeamPointPressureCondition::BeamPointPressureCondition(IndexType NewId, GeometryType::Pointer
									   pGeometry)
  : Condition(NewId, pGeometry)
  {
    //DO NOT ADD DOFS HERE!!!

  }

  //************************************************************************************
  //************************************************************************************
  BeamPointPressureCondition::BeamPointPressureCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
  : Condition(NewId, pGeometry, pProperties)
  {
  }

  //************************************************************************************
  //************************************************************************************
  BeamPointPressureCondition::BeamPointPressureCondition( BeamPointPressureCondition const& rOther )
  : Condition(rOther)
  {
  }

  //************************************************************************************
  //************************************************************************************

  Condition::Pointer BeamPointPressureCondition::Create(IndexType NewId, NodesArrayType
								 const& ThisNodes,  PropertiesType::Pointer pProperties) const
  {
    return Condition::Pointer(new BeamPointPressureCondition(NewId,GetGeometry().Create(ThisNodes), pProperties));
  }


  //************************************************************************************
  //************************************************************************************


  BeamPointPressureCondition::~BeamPointPressureCondition()
  {

  }

  //************************************************************************************
  //************************************************************************************
 
  void BeamPointPressureCondition::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
  {


  }


  //************************************************************************************
  //************************************************************************************
  void BeamPointPressureCondition::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
  {
   
  }
    
 


} // Namespace Kratos



