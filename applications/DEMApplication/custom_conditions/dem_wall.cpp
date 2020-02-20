//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

// Project includes
#include "custom_conditions/dem_wall.h"
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

void DEMWall::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMWall::Initialize) shouldn't be accessed, use derived class instead", 0);
}

//***********************************************************************************
//***********************************************************************************

void DEMWall::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    ProcessInfo& r_process_info) {
      KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMWall::CalculateRightHandSide) shouldn't be accessed, use derived class instead", 0);
}

void DEMWall::CalculateElasticForces(
    VectorType& rRightHandSideVector,
    ProcessInfo& r_process_info) {
      KRATOS_THROW_ERROR(std::runtime_error, "This function (DEMWall::CalculateElasticForces) shouldn't be accessed, use derived class instead", 0);
}

void DEMWall::GetDeltaDisplacement( array_1d<double, 3> & delta_displacement, int inode)
{
  delta_displacement = this->GetGeometry()[inode].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
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
double DEMWall::GetTgOfFrictionAngle()                                          { return GetProperties()[FRICTION]; }
double DEMWall::GetPoisson()                                                    { return GetProperties()[POISSON_RATIO]; }


void DEMWall::FinalizeSolutionStep(ProcessInfo& r_process_info)
{

}

void DEMWall::GetRightHadSideVector(std::vector<array_1d <double, 3> >& rRightHandSideVector)
{

  for(unsigned int a = 0; a< mRightHandSideVector.size(); a++)
  {
    for(unsigned int b = 0; b< 3; b++)
    {
      rRightHandSideVector[a][b] = mRightHandSideVector[a][b];
    }
  }

}

void DEMWall::SetRightHadSideVector(const std::vector<array_1d <double, 3> >& rRightHandSideVector)
{
  for(unsigned int a = 0; a< mRightHandSideVector.size(); a++)
  {
    for(unsigned int b = 0; b< 3; b++)
    {
      mRightHandSideVector[a][b] = rRightHandSideVector[a][b];
    }
  }

}

void DEMWall::AddToRightHadSideVector(const std::vector<array_1d <double, 3> >& rRightHandSideVector)
{
  for(unsigned int a = 0; a< mRightHandSideVector.size(); a++)
  {
    for(unsigned int b = 0; b< 3; b++)
    {
      mRightHandSideVector[a][b] += rRightHandSideVector[a][b];
    }
  }

}

//***********************************************************************************
//***********************************************************************************

} // Namespace Kratos.
