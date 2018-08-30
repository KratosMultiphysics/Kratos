//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		BSD License
//					Kratos default license: kratos/license.txt
//
//  Main authors:    Bodhinanda Chandra
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "custom_utilities/stress_invariants_utilities.hpp"
#include "custom_constitutive/yield_criteria/modified_cam_clay_yield_criterion.hpp"
#include "includes/mat_variables.h"


namespace Kratos
{

//*******************************CONSTRUCTOR******************************************
//************************************************************************************
ModifiedCamClayYieldCriterion::ModifiedCamClayYieldCriterion()
  :YieldCriterion()
{

}

//*****************************INITIALIZATION CONSTRUCTOR*****************************
//************************************************************************************

ModifiedCamClayYieldCriterion::ModifiedCamClayYieldCriterion(HardeningLawPointer pHardeningLaw)
  :YieldCriterion(pHardeningLaw)
{

}


//*******************************ASSIGMENT OPERATOR***********************************
//************************************************************************************

ModifiedCamClayYieldCriterion& ModifiedCamClayYieldCriterion::operator=(ModifiedCamClayYieldCriterion const& rOther)
{
  YieldCriterion::operator=(rOther);
  return *this;
}

//*******************************COPY CONSTRUCTOR*************************************
//************************************************************************************

ModifiedCamClayYieldCriterion::ModifiedCamClayYieldCriterion(ModifiedCamClayYieldCriterion const& rOther)
  :YieldCriterion(rOther)
{

}


//********************************DESTRUCTOR******************************************
//************************************************************************************

ModifiedCamClayYieldCriterion::~ModifiedCamClayYieldCriterion()
{
}



//************************* CALCULATE YIELD FUNCTION  ******************
//**********************************************************************

double& ModifiedCamClayYieldCriterion::CalculateYieldCondition(double& rStateFunction, const Vector& rStressVector, const double& rAlpha)
{

  return rStateFunction; 
}


//*******************************CALCULATE FIRST YIELD FUNCTION DERIVATIVE *****************
//************************************************************************************
void ModifiedCamClayYieldCriterion::CalculateYieldFunctionDerivative(const Vector& rStressVector, Vector& rFirstDerivative, const double& rAlpha)
{

}

//*******************************CALCULATE SECOND YIELD FUNCTION DERIVATIVE *****************
//************************************************************************************
void ModifiedCamClayYieldCriterion::CalculateYieldFunctionSecondDerivative(const Vector& rStressVector, Vector& rSecondDerivative)
{

}

double ModifiedCamClayYieldCriterion::GetSmoothingLodeAngle()
{
  return 27.0*GetPI()/180.0;
}

double ModifiedCamClayYieldCriterion::GetPI()
{
  return std::atan(1.0)*4.0;
}

void ModifiedCamClayYieldCriterion::save( Serializer& rSerializer ) const
{
   KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, YieldCriterion )
}

void ModifiedCamClayYieldCriterion::load( Serializer& rSerializer )
{
   KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, YieldCriterion )
}


}  // namespace Kratos.
