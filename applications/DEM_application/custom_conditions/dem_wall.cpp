/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */
/* *********************************************************
*
*   Last Modified by:    $Author: Feng Chun $
*   Date:                $Date: 2013-10-10 
*   Revision:            $Revision: 1.0
*
* ***********************************************************/


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/dem_wall.h"
#include "DEM_application.h"

#include "custom_utilities/GeometryFunctions.h"

namespace Kratos
{
	using namespace GeometryFunctions;

//***********************************************************************************
//***********************************************************************************


// Constructor

DEMWall::DEMWall()
{
}

// Constructor

DEMWall::DEMWall(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
}

// Constructor

DEMWall::DEMWall(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
    //setting up the nodal degrees of freedom
}

//***********************************************************************************
//***********************************************************************************

Condition::Pointer DEMWall::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new DEMWall(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//***********************************************************************************
//***********************************************************************************
// Destructor

DEMWall::~DEMWall()
{
}

//***********************************************************************************
//***********************************************************************************

void DEMWall::Initialize()
{
  /*
  mTgOfFrictionAngle = GetProperties()[WALL_FRICTION];
  mYoungModulus = GetProperties()[YOUNG_MODULUS];
  mPoissonRatio = GetProperties()[POISSON_RATIO];
  */
  
  this->GetGeometry()[0].FastGetSolutionStepValue(NON_DIMENSIONAL_VOLUME_WEAR) = 0.0;
  this->GetGeometry()[1].FastGetSolutionStepValue(NON_DIMENSIONAL_VOLUME_WEAR) = 0.0;
  this->GetGeometry()[2].FastGetSolutionStepValue(NON_DIMENSIONAL_VOLUME_WEAR) = 0.0;
  
  this->GetGeometry()[0].FastGetSolutionStepValue(IMPACT_WEAR) = 0.0;
  this->GetGeometry()[1].FastGetSolutionStepValue(IMPACT_WEAR) = 0.0;
  this->GetGeometry()[2].FastGetSolutionStepValue(IMPACT_WEAR) = 0.0;
}

//***********************************************************************************
//***********************************************************************************

void DEMWall::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo) {
}

void DEMWall::CalculateElasticForces(
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo) {
}

void DEMWall::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo){
}


void DEMWall::CalculateNormal(array_1d<double, 3>& rnormal){
    
    array_1d<double, 3> v1, v2;
        
    v1[0] = GetGeometry()[1].X() - GetGeometry()[0].X();
    v1[1] = GetGeometry()[1].Y() - GetGeometry()[0].Y();
    v1[2] = GetGeometry()[1].Z() - GetGeometry()[0].Z();

    v2[0] = GetGeometry()[2].X() - GetGeometry()[0].X();
    v2[1] = GetGeometry()[2].Y() - GetGeometry()[0].Y();
    v2[2] = GetGeometry()[2].Z() - GetGeometry()[0].Z();

    MathUtils<double>::CrossProduct(rnormal, v1, v2);

    rnormal /= MathUtils<double>::Norm3(rnormal);
}
  
 void DEMWall::AddExplicitContribution(const VectorType& rRHS,
                         const Variable<VectorType>& rRHSVariable,
                         Variable<array_1d<double,3> >& rDestinationVariable,
                         const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY


    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

   if (rDestinationVariable == EXTERNAL_FORCE )
      {

    for (unsigned int i=0; i< number_of_nodes; i++)
      {
        int index = dimension * i;

        GetGeometry()[i].SetLock();

        array_1d<double, 3 > &ExternalForce = GetGeometry()[i].FastGetSolutionStepValue(EXTERNAL_FORCE);
        for(unsigned int j=0; j<dimension; j++) {
            ExternalForce[j] += rRHS[index + j];
        }

        GetGeometry()[i].UnSetLock();
      }
      }

    if (rDestinationVariable == FORCE_RESIDUAL )
      {

    for (unsigned int i=0; i< number_of_nodes; i++)
      {
        int index = dimension * i;

        GetGeometry()[i].SetLock();

        array_1d<double, 3 > &ForceResidual = GetGeometry()[i].FastGetSolutionStepValue(FORCE_RESIDUAL);
        for(unsigned int j=0; j<dimension; j++)
          {
        ForceResidual[j] += rRHS[index + j];
          }

        GetGeometry()[i].UnSetLock();
      }
      }

    KRATOS_CATCH( "" )
}

 
double DEMWall::GetYoung()                                                      { return GetProperties()[YOUNG_MODULUS]; }
double DEMWall::GetTgOfFrictionAngle()                                          { return GetProperties()[WALL_FRICTION]; }
double DEMWall::GetPoisson()                                                    { return GetProperties()[POISSON_RATIO]; }


void DEMWall::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo)
{
    
}

void DEMWall::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)   
{  
  
}

//***********************************************************************************
//***********************************************************************************

} // Namespace Kratos.
