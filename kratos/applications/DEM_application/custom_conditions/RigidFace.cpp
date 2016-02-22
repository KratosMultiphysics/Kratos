//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

// Project includes
#include "includes/define.h"
#include "custom_conditions/RigidFace.h"
#include "DEM_application.h"

#include "custom_utilities/GeometryFunctions.h"

namespace Kratos {
    
    using namespace GeometryFunctions;

//***********************************************************************************
//***********************************************************************************

// Constructor

RigidFace3D::RigidFace3D() {}

// Constructor

RigidFace3D::RigidFace3D(IndexType NewId, GeometryType::Pointer pGeometry) : DEMWall(NewId, pGeometry) {}

// Constructor

RigidFace3D::RigidFace3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
           : DEMWall(NewId, pGeometry, pProperties) {
    //setting up the nodal degrees of freedom
    }

//***********************************************************************************
//***********************************************************************************

Condition::Pointer RigidFace3D::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new RigidFace3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//***********************************************************************************
//***********************************************************************************
// Destructor

RigidFace3D::~RigidFace3D() {}


//***********************************************************************************
//***********************************************************************************

void RigidFace3D::Initialize() {
    /*
  mTgOfFrictionAngle = GetProperties()[WALL_FRICTION];
    mYoungModulus      = GetProperties()[YOUNG_MODULUS];
    mPoissonRatio      = GetProperties()[POISSON_RATIO];
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

void RigidFace3D::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                         ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * 3;
    
    if (rRightHandSideVector.size() != MatSize)
    {
        rRightHandSideVector.resize(MatSize, false);
    }
    rRightHandSideVector = ZeroVector(MatSize);
    
    std::vector<SphericParticle*>& rNeighbours = this->mNeighbourSphericParticles;
    
    for (unsigned int i=0; i<rNeighbours.size(); i++)
    {
        if(rNeighbours[i]->Is(BLOCKED)) continue; //Inlet Generator Spheres are ignored when integrating forces.
        
        std::vector<DEMWall*>& rRFnei = rNeighbours[i]->mNeighbourRigidFaces;

        for (unsigned int i_nei = 0; i_nei < rRFnei.size(); i_nei++)
        {
            if ( rRFnei[i_nei]->Id() == this->Id() )
            {
                double weight[4] = {0.0};
                double ContactForce[3] = {0.0};

                unsigned int ino = 16 * i_nei;
                const std::vector<double>& neighbour_rigid_faces_pram = rNeighbours[i]->mNeighbourRigidFacesPram;
                
                weight[0] = neighbour_rigid_faces_pram[ino + 10];
                weight[1] = neighbour_rigid_faces_pram[ino + 11];
                weight[2] = neighbour_rigid_faces_pram[ino + 12];
                weight[3] = neighbour_rigid_faces_pram[ino + 13];

                const array_1d<double, 3>& neighbour_rigid_faces_contact_force = rNeighbours[i]->mNeighbourRigidFacesTotalContactForce[i_nei];
                ContactForce[0] = neighbour_rigid_faces_contact_force[0];
                ContactForce[1] = neighbour_rigid_faces_contact_force[1];
                ContactForce[2] = neighbour_rigid_faces_contact_force[2];

                for (unsigned int inode = 0; inode < GetGeometry().size(); inode++)
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

void RigidFace3D::CalculateElasticForces(VectorType& rElasticForces,
                                         ProcessInfo& rCurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * 3;
    
    if (rElasticForces.size() != MatSize)
    {
        rElasticForces.resize(MatSize, false);
    }
    rElasticForces = ZeroVector(MatSize);
    
    std::vector<SphericParticle*>& rNeighbours = this->mNeighbourSphericParticles;
    
    for (unsigned int i=0; i<rNeighbours.size(); i++)
    {
        std::vector<DEMWall*>& rRFnei = rNeighbours[i]->mNeighbourRigidFaces;

        for (unsigned int i_nei = 0; i_nei < rRFnei.size(); i_nei++)
        {
            if ( rRFnei[i_nei]->Id() == this->Id() )
            {
                double weight[4] = {0.0};
                double ContactElasticForce[3] = {0.0};

                unsigned int ino = 16 * i_nei;
                const std::vector<double>& neighbour_rigid_faces_pram = rNeighbours[i]->mNeighbourRigidFacesPram;
                
                weight[0] = neighbour_rigid_faces_pram[ino + 10];
                weight[1] = neighbour_rigid_faces_pram[ino + 11];
                weight[2] = neighbour_rigid_faces_pram[ino + 12];
                weight[3] = neighbour_rigid_faces_pram[ino + 13];

                const array_1d<double, 3>& neighbour_rigid_faces_elastic_contact_force = rNeighbours[i]->mNeighbourRigidFacesElasticContactForce[i_nei];                    
                ContactElasticForce[0] = neighbour_rigid_faces_elastic_contact_force[0];
                ContactElasticForce[1] = neighbour_rigid_faces_elastic_contact_force[1];
                ContactElasticForce[2] = neighbour_rigid_faces_elastic_contact_force[2];

                for (unsigned int inode = 0; inode < GetGeometry().size(); inode++)
                {
                    unsigned int ino1 =  inode * 3;
                    rElasticForces[ino1 + 0] += -ContactElasticForce[0] * weight[inode];
                    rElasticForces[ino1 + 1] += -ContactElasticForce[1] * weight[inode];
                    rElasticForces[ino1 + 2] += -ContactElasticForce[2] * weight[inode];
                }
            }
        }
    }
}


void RigidFace3D::CalculateNormal(array_1d<double, 3>& rnormal){

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


void RigidFace3D::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo)
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

      double Xorigin    = rCurrentProcessInfo[RIGID_FACE_ROTA_ORIGIN_COORD][0];   
      double Yorigin    = rCurrentProcessInfo[RIGID_FACE_ROTA_ORIGIN_COORD][1];
      double Zorigin    = rCurrentProcessInfo[RIGID_FACE_ROTA_ORIGIN_COORD][2]; 
      
      ///movement of the original point
      int time_step           = rCurrentProcessInfo[TIME_STEPS];			
      double begin_time       = rCurrentProcessInfo[RIGID_FACE_BEGIN_TIME];
      double real_rota_time   = delta_t * time_step - begin_time;
          
      
      double n[3] = {Xnormal, Ynormal, Znormal};
      GeometryFunctions::normalize(n);

      double omiga = CyclePerSec * 2.0 * KRATOS_M_PI;
      
      double vel = NormalV;

      double g_v[3] = {GXvel, GYvel, GZvel};

      Xorigin += (g_v[0] + n[0] * vel) * real_rota_time; 
      Yorigin += (g_v[1] + n[1] * vel) * real_rota_time; 
      Zorigin += (g_v[2] + n[2] * vel) * real_rota_time; 

      
      double origin[3] = {Xorigin, Yorigin, Zorigin};

      double vector1[3], vector2[3];
      double dist, dist1;

      double a[3][3];
      double local_vel[3],global_vel[3];
      
        for(unsigned int j = 0; j < number_of_nodes; j++)
        {
          const array_1d<double, 3>& Nodecoord = this->GetGeometry()[j].Coordinates();

          vector1[0] = Nodecoord[0] - origin[0];
          vector1[1] = Nodecoord[1] - origin[1];
          vector1[2] = Nodecoord[2] - origin[2];

          dist  = fabs(GeometryFunctions::DotProduct(vector1,n));
          dist1 = GeometryFunctions::DistanceOfTwoPoint(Nodecoord,origin);

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

void RigidFace3D::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)   
{
  
  
}

//***********************************************************************************
//***********************************************************************************

} // Namespace Kratos.
