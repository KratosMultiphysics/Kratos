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

// External includes

// Project includes
#include "includes/define.h"
#include "custom_conditions/RigidEdge.h"
#include "DEM_application.h"

#include "custom_utilities/GeometryFunctions.h"

namespace Kratos
{
	using namespace GeometryFunctions;
//************************************************************************************
//************************************************************************************
RigidEdge3D::RigidEdge3D( IndexType NewId,
        GeometryType::Pointer pGeometry)
    : Condition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//**** life cycle ********************************************************************
//************************************************************************************
RigidEdge3D::RigidEdge3D( IndexType NewId, GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties
                                            )
    : Condition( NewId, pGeometry, pProperties )
{
}

RigidEdge3D::RigidEdge3D( IndexType NewId, GeometryType::Pointer pGeometry,
        PropertiesType::Pointer pProperties,
        Condition::Pointer Master,
        Condition::Pointer Slave,
        Point<3>& MasterContactLocalPoint,
        Point<3>& SlaveContactLocalPoint,
        int SlaveIntegrationPointIndex
                                            )
    : Condition( NewId, pGeometry, pProperties )
{

}

//********************************************************
//**** Operations ****************************************
//********************************************************


Condition::Pointer RigidEdge3D::Create( IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer( new RigidEdge3D(NewId, GetGeometry().Create(ThisNodes),
                               pProperties));
}
/**
 * Destructor. Never to be called manually
 */
RigidEdge3D::~RigidEdge3D()
{
}

//************************************************************************************
//************************************************************************************

/**
 * calculates only the RHS vector (certainly to be removed due to contact algorithm)
 */
void RigidEdge3D::CalculateRightHandSide( VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
  const unsigned int number_of_nodes = GetGeometry().size();
   unsigned int               MatSize = number_of_nodes * 3;

	if (rRightHandSideVector.size() != MatSize)
	{
		rRightHandSideVector.resize(MatSize, false);
	}
	rRightHandSideVector = ZeroVector(MatSize); 
	
	
	ParticleWeakVectorType& rNeighbours    = this->GetValue(NEIGHBOUR_PARTICLE_OF_RIGID_FACE);

	
	for (ParticleWeakIteratorType neighbour_iterator = rNeighbours.begin(); neighbour_iterator != rNeighbours.end(); neighbour_iterator++)
	{
		ConditionWeakVectorType& rRFnei    = neighbour_iterator->GetValue(NEIGHBOUR_RIGID_FACES);
		
		for(unsigned int i_nei = 0; i_nei < rRFnei.size(); i_nei++)
		{
			if( rRFnei[i_nei].Id() == this->Id() )
			{
				double weight[4] = {0.0};
				double ContactForce[3] = {0.0};
				
				unsigned int ino = 15 * i_nei;
				
				weight[0] = neighbour_iterator->GetValue(NEIGHBOUR_RIGID_FACES_PRAM)[ino + 10];
				weight[1] = neighbour_iterator->GetValue(NEIGHBOUR_RIGID_FACES_PRAM)[ino + 11];
				weight[2] = neighbour_iterator->GetValue(NEIGHBOUR_RIGID_FACES_PRAM)[ino + 12];
				weight[3] = neighbour_iterator->GetValue(NEIGHBOUR_RIGID_FACES_PRAM)[ino + 13];
				
				ino = 3 * i_nei;
				
				ContactForce[0] = neighbour_iterator->GetValue(NEIGHBOUR_RIGID_FACES_CONTACT_FORCE)[ino + 0];
				ContactForce[1] = neighbour_iterator->GetValue(NEIGHBOUR_RIGID_FACES_CONTACT_FORCE)[ino + 1];
				ContactForce[2] = neighbour_iterator->GetValue(NEIGHBOUR_RIGID_FACES_CONTACT_FORCE)[ino + 2];
				
				for(unsigned int inode = 0; inode < GetGeometry().size(); inode++ )
				{
					unsigned int ino1 =  inode * 3;
					
					rRightHandSideVector[ino1 + 0] += -ContactForce[0] * weight[inode];
					rRightHandSideVector[ino1 + 1] += -ContactForce[1] * weight[inode];
					rRightHandSideVector[ino1 + 2] += -ContactForce[2] * weight[inode];
					
				}
				
				
			}
		}
	
	}
}


void RigidEdge3D::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo)
{
  if (rVariable == RIGID_FACE_COMPUTE_MOVEMENT)
  {
	const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int               MatSize = number_of_nodes * 3;

	if (Output.size() != MatSize)
	{
		Output.resize(MatSize, false);
	}
	Output = ZeroVector(MatSize); 
	  
	double delta_t     = rCurrentProcessInfo[DELTA_TIME];
	double CyclePerSec = rCurrentProcessInfo[RIGID_FACE_ROTA_SPEED];
	double NormalV     = rCurrentProcessInfo[RIGID_FACE_AXIAL_SPEED];

	double GXvel       = rCurrentProcessInfo[RIGID_FACE_ROTA_GLOBAL_VELOCITY][0];
	double GYvel       = rCurrentProcessInfo[RIGID_FACE_ROTA_GLOBAL_VELOCITY][1];
	double GZvel       = rCurrentProcessInfo[RIGID_FACE_ROTA_GLOBAL_VELOCITY][2];

	double Xnormal     = rCurrentProcessInfo[RIGID_FACE_ROTA_AXIAL_DIR][0];
	double Ynormal     = rCurrentProcessInfo[RIGID_FACE_ROTA_AXIAL_DIR][1];
	double Znormal     = rCurrentProcessInfo[RIGID_FACE_ROTA_AXIAL_DIR][2];

	double  Xorigin    = rCurrentProcessInfo[RIGID_FACE_ROTA_ORIGIN_COORD][0];   
	double  Yorigin    = rCurrentProcessInfo[RIGID_FACE_ROTA_ORIGIN_COORD][1];
	double  Zorigin    = rCurrentProcessInfo[RIGID_FACE_ROTA_ORIGIN_COORD][2]; 
	
	///movement of the original point
	int time_step           = rCurrentProcessInfo[TIME_STEPS];			
	double begin_time       = rCurrentProcessInfo[RIGID_FACE_BEGIN_TIME];
	double real_rota_time   = delta_t * time_step - begin_time;
			
	
	double n[3] = {Xnormal, Ynormal, Znormal};
	GeometryFunctions::normalize(n);

	double omiga = CyclePerSec * 2.0 * M_PI;
	
	double vel = NormalV;

	double g_v[3] = {GXvel, GYvel, GZvel};

	Xorigin += (g_v[0] + n[0] * vel) * real_rota_time; 
	Yorigin += (g_v[1] + n[1] * vel) * real_rota_time; 
	Zorigin += (g_v[2] + n[2] * vel) * real_rota_time; 

	
	double origin[3] = {Xorigin, Yorigin, Zorigin};

	double coord[3], vector1[3], vector2[3];
	double dist, dist1;

	double a[3][3];
	double local_vel[3],global_vel[3];
	
	for(unsigned int j = 0; j < number_of_nodes; j++)
	{
		array_1d<double, 3> Nodecoord;
		Nodecoord = this->GetGeometry()(j)->Coordinates();
		
		coord[0] = Nodecoord[0];
		coord[1] = Nodecoord[1];
		coord[2] = Nodecoord[2];

		vector1[0] = coord[0] - origin[0];
		vector1[1] = coord[1] - origin[1];
		vector1[2] = coord[2] - origin[2];


		dist  = fabs(GeometryFunctions::DotProduct(vector1,n));
		dist1 = GeometryFunctions::DistanceOfTwoPoint(coord,origin);

		dist = sqrt( dist1 * dist1 - dist * dist);

		if(dist < 1.0e-6)
		{
			global_vel[0] = n[0] * vel;
			global_vel[1] = n[1] * vel;
			global_vel[2] = n[2] * vel;
		}
		else
		{
			local_vel[0] = 0.0;
			local_vel[1] = dist * omiga;
			local_vel[2] = vel;

			GeometryFunctions::normalize(vector1);
			
			GeometryFunctions::CrossProduct(n,vector1,vector2);
			
			GeometryFunctions::normalize(vector2);  
			
			GeometryFunctions::CrossProduct(vector2,n,vector1);
			
			GeometryFunctions::normalize(vector1);

			a[0][0] = vector1[0];
			a[0][1] = vector1[1];
			a[0][2] = vector1[2];

			a[1][0] = vector2[0];
			a[1][1] = vector2[1];
			a[1][2] = vector2[2];

			a[2][0] = n[0];
			a[2][1] = n[1];
			a[2][2] = n[2];

			GeometryFunctions::VectorLocal2Global(a,local_vel,global_vel);	
		}
		
		Output[3 * j + 0] = (global_vel[0] + g_v[0]);
		Output[3 * j + 1] = (global_vel[1] + g_v[1]);
		Output[3 * j + 2] = (global_vel[2] + g_v[2]);			
	}
  }
}



} // Namespace Kratos
