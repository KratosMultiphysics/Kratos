/*
Kratos Multi-Physics

Copyright (c) 2015, Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
		in the documentation and/or other materials provided with the distribution.
	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
			This product includes Kratos Multi-Physics technology.
	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


// Project includes
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
    KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMWall::Initialize) shouldn't be accessed, use derived class instead", 0);
}

//***********************************************************************************
//***********************************************************************************

void DEMWall::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& r_process_info) {
}

void DEMWall::CalculateElasticForces(
    VectorType& rRightHandSideVector,
    ProcessInfo& r_process_info) {
}

void DEMWall::InitializeSolutionStep(ProcessInfo& r_process_info){
}


void DEMWall::CalculateNormal(array_1d<double, 3>& rnormal){

   KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMWall::CalculateNormal) shouldn't be accessed, use derived class instead", "");
}

 void DEMWall::AddExplicitContribution(const VectorType& rRHS,
                         const Variable<VectorType>& rRHSVariable,
                         Variable<array_1d<double,3> >& rDestinationVariable,
                         const ProcessInfo& r_process_info)
{
    KRATOS_TRY


    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int dimension       = GetGeometry().WorkingSpaceDimension();

   if( rDestinationVariable == EXTERNAL_FORCE )
      {

    for(unsigned int i=0; i< number_of_nodes; i++)
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

    if( rDestinationVariable == FORCE_RESIDUAL )
      {

    for(unsigned int i=0; i< number_of_nodes; i++)
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


void DEMWall::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& r_process_info)
{

}

void DEMWall::FinalizeSolutionStep(ProcessInfo& r_process_info)
{

}

//***********************************************************************************
//***********************************************************************************

} // Namespace Kratos.
