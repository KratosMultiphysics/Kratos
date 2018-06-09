//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
//

// Project includes
#include "custom_conditions/RigidFace.h"
#include "custom_elements/spheric_particle.h"
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

    const unsigned int number_of_nodes = GetGeometry().size();

    for (unsigned int i=0; i< number_of_nodes; i++)
    {
        this->GetGeometry()[i].FastGetSolutionStepValue(NON_DIMENSIONAL_VOLUME_WEAR) = 0.0;
        this->GetGeometry()[i].FastGetSolutionStepValue(IMPACT_WEAR) = 0.0;
    }
}

//***********************************************************************************
//***********************************************************************************

void RigidFace3D::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                         ProcessInfo& r_process_info)
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
            int Contact_Type = rNeighbours[i]->mContactConditionContactTypes[i_nei];
            
            if ( ( rRFnei[i_nei]->Id() == this->Id() ) && (Contact_Type > 0 ) )
            {
                
                array_1d<double, 4> weights_vector = rNeighbours[i]->mContactConditionWeights[i_nei];
                double weight = 0.0;
                
                double ContactForce[3] = {0.0};

                const array_1d<double, 3>& neighbour_rigid_faces_contact_force = rNeighbours[i]->mNeighbourRigidFacesTotalContactForce[i_nei];

                ContactForce[0] = neighbour_rigid_faces_contact_force[0];
                ContactForce[1] = neighbour_rigid_faces_contact_force[1];
                ContactForce[2] = neighbour_rigid_faces_contact_force[2];

                for (unsigned int k=0; k< number_of_nodes; k++)
                {
                    weight = weights_vector[k];
  
                    unsigned int w =  k * 3;

                    rRightHandSideVector[w + 0] += -ContactForce[0] * weight;
                    rRightHandSideVector[w + 1] += -ContactForce[1] * weight;
                    rRightHandSideVector[w + 2] += -ContactForce[2] * weight;
                }
                
            }//if the condition neighbour of my sphere neighbour is myself.
        }//Loop spheres neighbours (condition)
    }//Loop condition neighbours (spheres)
}//CalculateRightHandSide

void RigidFace3D::CalculateElasticForces(VectorType& rElasticForces,
                                         ProcessInfo& r_process_info)
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
        
        if(rNeighbours[i]->Is(BLOCKED)) continue; //Inlet Generator Spheres are ignored when integrating forces.
          
        std::vector<DEMWall*>& rRFnei = rNeighbours[i]->mNeighbourRigidFaces;

        for (unsigned int i_nei = 0; i_nei < rRFnei.size(); i_nei++)
        {
            int Contact_Type = rNeighbours[i]->mContactConditionContactTypes[i_nei];
            
            if ( ( rRFnei[i_nei]->Id() == this->Id() ) && (Contact_Type > 0 ) )
            {
                array_1d<double, 4> weights_vector = rNeighbours[i]->mContactConditionWeights[i_nei];
                double weight = 0.0;
                
                double ContactElasticForce[3] = {0.0};

                const array_1d<double, 3>& neighbour_rigid_faces_elastic_contact_force = rNeighbours[i]->mNeighbourRigidFacesElasticContactForce[i_nei];                    
                ContactElasticForce[0] = neighbour_rigid_faces_elastic_contact_force[0];
                ContactElasticForce[1] = neighbour_rigid_faces_elastic_contact_force[1];
                ContactElasticForce[2] = neighbour_rigid_faces_elastic_contact_force[2];

                for (unsigned int k=0; k< number_of_nodes; k++)
                {
                    weight = weights_vector[k];
  
                    unsigned int w =  k * 3;

                    rElasticForces[w + 0] += -ContactElasticForce[0] * weight;
                    rElasticForces[w + 1] += -ContactElasticForce[1] * weight;
                    rElasticForces[w + 2] += -ContactElasticForce[2] * weight;
                }
            }//if the condition neighbour of my sphere neighbour is myself.
        }//Loop spheres neighbours (condition)
    }//Loop condition neighbours (spheres)
}//CalculateRightHandSide


void RigidFace3D::CalculateNormal(array_1d<double, 3>& rnormal){

//    array_1d<double, 3> v1, v2;
//    Geometry<Node<3> >& geom = GetGeometry();

//    v1[0] = geom[1][0] - geom[0][0];
//    v1[1] = geom[1][1] - geom[0][1];
//    v1[2] = geom[1][2] - geom[0][2];

//    v2[0] = geom[2][0] - geom[0][0];
//    v2[1] = geom[2][1] - geom[0][1];
//    v2[2] = geom[2][2] - geom[0][2];

//    MathUtils<double>::CrossProduct(rnormal, v1, v2);

//    rnormal /= MathUtils<double>::Norm3(rnormal);
    Geometry<Node<3> >& geom = GetGeometry();
    if (geom.size()>2){
        double v1[3];
        double v2[3];

        double p0[3] = {geom[0][0], geom[0][1], geom[0][2]};
        double p1[3] = {geom[1][0], geom[1][1], geom[1][2]};
        double p2[3] = {geom[2][0], geom[2][1], geom[2][2]};

        v1[0] = p1[0] - p0[0];
        v1[1] = p1[1] - p0[1];
        v1[2] = p1[2] - p0[2];

        v2[0] = p2[0] - p0[0];
        v2[1] = p2[1] - p0[1];
        v2[2] = p2[2] - p0[2];
        DEM_SET_TO_CROSS_OF_FIRST_TWO_3(v1, v2, rnormal)
        const double norm_n_inv = 1.0 / DEM_MODULUS_3(rnormal) ;
        DEM_MULTIPLY_BY_SCALAR_3(rnormal, norm_n_inv);
    }
    else{
        KRATOS_ERROR << "Calculating normal direction for line or point is not possible" << std::endl;
    }
}


void RigidFace3D::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& r_process_info)
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
        
      double delta_t     = r_process_info[DELTA_TIME];
      double CyclePerSec = r_process_info[RIGID_FACE_ROTA_SPEED];
      double NormalV     = r_process_info[RIGID_FACE_AXIAL_SPEED];

      double GXvel       = r_process_info[RIGID_FACE_ROTA_GLOBAL_VELOCITY][0];
      double GYvel       = r_process_info[RIGID_FACE_ROTA_GLOBAL_VELOCITY][1];
      double GZvel       = r_process_info[RIGID_FACE_ROTA_GLOBAL_VELOCITY][2];

      double Xnormal     = r_process_info[RIGID_FACE_ROTA_AXIAL_DIR][0];
      double Ynormal     = r_process_info[RIGID_FACE_ROTA_AXIAL_DIR][1];
      double Znormal     = r_process_info[RIGID_FACE_ROTA_AXIAL_DIR][2];

      double Xorigin    = r_process_info[RIGID_FACE_ROTA_ORIGIN_COORD][0];   
      double Yorigin    = r_process_info[RIGID_FACE_ROTA_ORIGIN_COORD][1];
      double Zorigin    = r_process_info[RIGID_FACE_ROTA_ORIGIN_COORD][2]; 
      
      ///movement of the original point
      int time_step           = r_process_info[TIME_STEPS];			
      double begin_time       = r_process_info[RIGID_FACE_BEGIN_TIME];
      double real_rota_time   = delta_t * time_step - begin_time;
          
      
      double n[3] = {Xnormal, Ynormal, Znormal};
      GeometryFunctions::normalize(n);

      double omiga = CyclePerSec * 2.0 * Globals::Pi;
      
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

array_1d<double, 3> RigidFace3D::GetVelocity() {
        
    size_t FE_size = this->GetGeometry().size();
    array_1d<double, 3> rigid_face_velocity = ZeroVector(3);
    double factor = 1.0;
    
    for (std::size_t inode = 0; inode < FE_size; inode++) {
        
        DEM_ADD_SECOND_TO_FIRST(rigid_face_velocity, this->GetGeometry()[inode].FastGetSolutionStepValue(VELOCITY))
    }
    
    if (FE_size) factor /= FE_size;
    
    DEM_MULTIPLY_BY_SCALAR_3(rigid_face_velocity, factor)
    
    return rigid_face_velocity;
}

void RigidFace3D::ComputeConditionRelativeData(int rigid_neighbour_index,
                                               SphericParticle* const particle,
                                               double LocalCoordSystem[3][3],
                                               double& DistPToB,
                                               array_1d<double, 4>& Weight,
                                               array_1d<double, 3>& wall_delta_disp_at_contact_point,
                                               array_1d<double, 3>& wall_velocity_at_contact_point,
                                               int& ContactType)
{
    size_t FE_size = this->GetGeometry().size();
    
    std::vector<double> TempWeight;
    TempWeight.resize(FE_size);

    double total_weight = 0.0;
    int points = 0;
    int inode1 = 0, inode2 = 0;

    for (unsigned int inode = 0; inode < FE_size; inode++) {

        if (Weight[inode] > 1.0e-12) {
            total_weight = total_weight + Weight[inode];
            points++;
            if (points == 1) {inode1 = inode;}
            if (points == 2) {inode2 = inode;}
        }

        if (fabs(total_weight - 1.0) < 1.0e-12) {
            break;
        }
    }

    bool contact_exists = true;
    array_1d<double, 3>& node_coordinates = particle->GetGeometry()[0].Coordinates();

    const double radius = particle->GetInteractionRadius();

    if (points == 3 || points == 4)
    {
        unsigned int dummy_current_edge_index;
        contact_exists = GeometryFunctions::FacetCheck(this->GetGeometry(), node_coordinates, radius, LocalCoordSystem, DistPToB, TempWeight, dummy_current_edge_index);
        ContactType = 1;
        Weight[0]=TempWeight[0];
        Weight[1]=TempWeight[1];
        Weight[2]=TempWeight[2];
        if (points == 4)
        {
            Weight[3] = TempWeight[3];
        }
        else
        {
            Weight[3] = 0.0;
        }
    }

    else if (points == 2) {

        double eta = 0.0;
        contact_exists = GeometryFunctions::EdgeCheck(this->GetGeometry()[inode1], this->GetGeometry()[inode2], node_coordinates, radius, LocalCoordSystem, DistPToB, eta);

        Weight[inode1] = 1-eta;
        Weight[inode2] = eta;
        ContactType = 2;

    }

    else if (points == 1) {
        contact_exists = GeometryFunctions::VertexCheck(this->GetGeometry()[inode1], node_coordinates, radius, LocalCoordSystem, DistPToB);
        Weight[inode1] = 1.0;
        ContactType = 3;
    }

    if (contact_exists == false) {ContactType = -1;}

    for (std::size_t inode = 0; inode < FE_size; inode++) {
        noalias(wall_velocity_at_contact_point) += this->GetGeometry()[inode].FastGetSolutionStepValue(VELOCITY) * Weight[inode];

        array_1d<double, 3>  wall_delta_displacement = ZeroVector(3);
        this->GetDeltaDisplacement(wall_delta_displacement, inode);
        noalias(wall_delta_disp_at_contact_point) += wall_delta_displacement* Weight[inode];

    }
}//ComputeConditionRelativeData

int RigidFace3D::CheckSide(SphericParticle *p_particle){
    const array_1d<double, 3>& particle_center_coors = p_particle->GetGeometry()[0];
    const Geometry<Node<3> >& geom = this->GetGeometry();
    const array_1d<double, 3> a0 = geom[1] - geom[0];
    const array_1d<double, 3> a1 = geom[2] - geom[0];
    const array_1d<double, 3> a2 = particle_center_coors - geom[0];
    const double ball_to_vertices_determinant = DEM_DETERMINANT_3x3(a0, a1, a2);  // each side corresponds to a different sign of this determinant

    return RigidFace3D::Sign(ball_to_vertices_determinant);
}

// The following function checks if the projection of a point onto the plane defined by this face
// falls inside it or not. Algorithm taken from https://math.stackexchange.com/questions/544946/determine-if-projection-of-3d-point-onto-plane-is-within-a-triangle
bool RigidFace3D::CheckProjectionFallsInside(SphericParticle *p_particle)
{
    const array_1d<double, 3>& P = p_particle->GetGeometry()[0].Coordinates();
    const Geometry<Node<3> >& geom = GetGeometry();
    const array_1d<double, 3> w  = P - geom[0].Coordinates();
    array_1d<double, 3> u1 = geom[1].Coordinates() - geom[0].Coordinates();
    array_1d<double, 3> u2 = geom[2].Coordinates() - geom[0].Coordinates();
    array_1d<double, 3> u2_copy;
    noalias(u2_copy) = u2;
    array_1d<double, 3> n;
    GeometryFunctions::CrossProduct(u1, u2, n);
    GeometryFunctions::CrossProduct(w, u2_copy, u2);
    const double beta = DEM_INNER_PRODUCT_3(u2, n);
    GeometryFunctions::CrossProduct(u1, w, u1);
    const double gamma = DEM_INNER_PRODUCT_3(u1, n);
    const double n2 = DEM_INNER_PRODUCT_3(n, n);
    const double alpha = n2 - beta - gamma;

    const bool falls_inside = alpha >=  0 && beta >=  0 && gamma >= 0
                           && alpha <= n2 && beta <= n2 && gamma <= n2;

    return falls_inside;
}

void RigidFace3D::FinalizeSolutionStep(ProcessInfo& r_process_info)   
{
  
}

//***********************************************************************************
//***********************************************************************************

} // Namespace Kratos.
