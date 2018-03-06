//
// Author: Guillermo Casas gcasas@cimne.upc.edu
//

// System includes
#include <string>
#include <iostream>
#include <cmath>

#include <fstream>

// External includes

// Project includes
#include "analytic_spheric_particle.h"

namespace Kratos
{
// using namespace GeometryFunctions;

AnalyticSphericParticle::AnalyticSphericParticle()
    : SphericParticle()
{
    ClearImpactMemberVariables();
}

AnalyticSphericParticle::AnalyticSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry)
    : SphericParticle(NewId, pGeometry)
{
    ClearImpactMemberVariables();
}

AnalyticSphericParticle::AnalyticSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : SphericParticle(NewId, pGeometry, pProperties)
{
    ClearImpactMemberVariables();
}

AnalyticSphericParticle::AnalyticSphericParticle(IndexType NewId, NodesArrayType const& ThisNodes)
    : SphericParticle(NewId, ThisNodes)
{
    ClearImpactMemberVariables();
}

AnalyticSphericParticle::AnalyticSphericParticle(Element::Pointer p_spheric_particle)
{
    GeometryType::Pointer p_geom = p_spheric_particle->pGetGeometry();
    PropertiesType::Pointer pProperties = p_spheric_particle->pGetProperties();
    AnalyticSphericParticle(p_spheric_particle->Id(), p_geom, pProperties);
}

AnalyticSphericParticle& AnalyticSphericParticle::operator=(const AnalyticSphericParticle& rOther) {

    SphericParticle::operator=(rOther);

    NeighboursContactStatus = rOther.NeighboursContactStatus;
    mNumberOfCollidingSpheres = rOther.mNumberOfCollidingSpheres;
    mNumberOfCollidingSpheresWithFaces = rOther.mNumberOfCollidingSpheresWithFaces;
    mNumberOfCollidingSpheresWithEdges = rOther.mNumberOfCollidingSpheresWithEdges;
    mCollidingIds = rOther.mCollidingIds;
    mCollidingRadii = rOther.mCollidingRadii;
    mCollidingNormalVelocities = rOther.mCollidingNormalVelocities;
    mCollidingTangentialVelocities = rOther.mCollidingTangentialVelocities;
    mCollidingLinearImpulse = rOther.mCollidingLinearImpulse;
    mContactingNeighbourIds = rOther.mContactingNeighbourIds;
    mCollidingFaceIds = rOther.mCollidingFaceIds;
    mCollidingFaceNormalVelocities = rOther.mCollidingFaceNormalVelocities;
    mCollidingFaceTangentialVelocities = rOther.mCollidingFaceTangentialVelocities;
    mCollidingFaceSecondTangentialVelocities = rOther.mCollidingFaceSecondTangentialVelocities;
    mCollidingFaceCollisionTypes = rOther.mCollidingFaceCollisionTypes;
    mContactingFaceNeighbourIds = rOther.mContactingFaceNeighbourIds;

    //Nothing done for std::unique_ptr<ParticleDataBuffer> mpDataBuffer;

    return *this;
}


int AnalyticSphericParticle::GetNumberOfCollisions(){return mNumberOfCollidingSpheres;}
int AnalyticSphericParticle::GetNumberOfCollisionsWithFaces(){return mNumberOfCollidingSpheresWithFaces;}
int AnalyticSphericParticle::GetNumberOfCollisionsWithEdges(){return mNumberOfCollidingSpheresWithEdges;}

array_1d<int, 4> &AnalyticSphericParticle::GetCollidingIds(){return mCollidingIds;}
array_1d<int, 4> &AnalyticSphericParticle::GetCollidingFaceIds(){return mCollidingFaceIds;}
array_1d<double, 4> &AnalyticSphericParticle::GetCollidingNormalRelativeVelocity(){return mCollidingNormalVelocities;}
array_1d<double, 4> &AnalyticSphericParticle::GetCollidingFaceNormalRelativeVelocity(){return mCollidingFaceNormalVelocities;}
array_1d<double, 4> &AnalyticSphericParticle::GetCollidingTangentialRelativeVelocity(){return mCollidingTangentialVelocities;}
array_1d<double, 4> &AnalyticSphericParticle::GetCollidingFaceTangentialRelativeVelocity(){return mCollidingFaceTangentialVelocities;}
array_1d<double, 4> &AnalyticSphericParticle::GetCollidingLinearImpulse(){return mCollidingLinearImpulse;}

/*
void AnalyticSphericParticle::GetCollidingFaceIds(array_1d<int, 4>& colliding_ids_with_walls)
{
    colliding_ids_with_walls = mCollidingFaceIds;
}

void AnalyticSphericParticle::GetCollidingNormalRelativeVelocity(array_1d<double, 4>& colliding_normal_vel)
{
    colliding_normal_vel = mCollidingNormalVelocities;
}

void AnalyticSphericParticle::GetCollidingFaceNormalRelativeVelocity(array_1d<double, 4>& colliding_normal_vel)
{
    colliding_normal_vel = mCollidingFaceNormalVelocities;
}

void AnalyticSphericParticle::GetCollidingTangentialRelativeVelocity(array_1d<double, 4>& colliding_tangential_vel)
{
    colliding_tangential_vel = mCollidingTangentialVelocities;
}

void AnalyticSphericParticle::GetCollidingFaceTangentialRelativeVelocity(array_1d<double, 4>& colliding_tangential_vel)
{
    colliding_tangential_vel = mCollidingFaceTangentialVelocities;
}

void AnalyticSphericParticle::GetCollidingLinearImpulse(array_1d<double, 4>& colliding_linear_impulse)
{
    colliding_linear_impulse = mCollidingLinearImpulse;
}
*/

void AnalyticSphericParticle::ClearImpactMemberVariables()
{
    mNumberOfCollidingSpheres = 0;
    mNumberOfCollidingSpheresWithFaces = 0;
    mNumberOfCollidingSpheresWithEdges = 0;

    for (unsigned int i = 0; i < mMaxCollidingSpheres; ++i){
        mCollidingIds[i] = 0;
        mCollidingRadii[i] = 0.0;
        mCollidingNormalVelocities[i] = 0.0;
        mCollidingTangentialVelocities[i] = 0.0;
        mCollidingLinearImpulse[i] = 0.0;

        mCollidingFaceIds[i] = 0;
        mCollidingFaceNormalVelocities[i] = 0.0;
        mCollidingFaceTangentialVelocities[i] = 0.0;
        mCollidingFaceSecondTangentialVelocities[i] = 0.0;
        //mCollidingFaceCollisionTypes[i] = "f";
        
    }
}

Element::Pointer AnalyticSphericParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    GeometryType::Pointer p_geom = GetGeometry().Create(ThisNodes);

    return Element::Pointer(new AnalyticSphericParticle(NewId, p_geom, pProperties));
}

void AnalyticSphericParticle::PushBackIdToContactingNeighbours(BaseBufferType & data_buffer, int id)
{
    GetPointerToDerivedDataBuffer(data_buffer)->mCurrentContactingNeighbourIds.push_back(id);
}

void AnalyticSphericParticle::PushBackIdToContactingFaceNeighbours(BaseBufferType & data_buffer, int p_wall_id)
{
    GetPointerToDerivedDataBuffer(data_buffer)->mCurrentContactingFaceNeighbourIds.push_back(p_wall_id);
}


void AnalyticSphericParticle::ClearNeighbours(BaseBufferType & data_buffer)
{
    GetPointerToDerivedDataBuffer(data_buffer)->mCurrentContactingNeighbourIds.clear();
    GetPointerToDerivedDataBuffer(data_buffer)->mCurrentContactingFaceNeighbourIds.clear();
}


void AnalyticSphericParticle::EvaluateBallToBallForcesForPositiveIndentiations(SphericParticle::ParticleDataBuffer & data_buffer,
                                                       const ProcessInfo& r_process_info,
                                                       double LocalElasticContactForce[3],
                                                       double DeltDisp[3],
                                                       double LocalDeltDisp[3],
                                                       double RelVel[3],
                                                       const double indentation,
                                                       double ViscoDampingLocalContactForce[3],
                                                       double& cohesive_force,
                                                       SphericParticle* p_neighbour_element,
                                                       bool& sliding,
                                                       double LocalCoordSystem[3][3],
                                                       double OldLocalCoordSystem[3][3],
                                                       array_1d<double, 3>& neighbour_elastic_contact_force)
{
    SphericParticle::EvaluateBallToBallForcesForPositiveIndentiations(data_buffer,
                                                                      r_process_info,
                                                                      LocalElasticContactForce,
                                                                      DeltDisp,
                                                                      LocalDeltDisp,
                                                                      RelVel,
                                                                      indentation,
                                                                      ViscoDampingLocalContactForce,
                                                                      cohesive_force,
                                                                      p_neighbour_element,   // ALREADY INCLUDED IN DATABUFFER
                                                                      sliding,
                                                                      LocalCoordSystem,
                                                                      OldLocalCoordSystem,
                                                                      neighbour_elastic_contact_force);

    const auto id = data_buffer.mpOtherParticle->Id();
    
    if (IsNewNeighbour(id) && mNumberOfCollidingSpheres < mMaxCollidingSpheres){
        RecordNewImpact(data_buffer);
    }

    PushBackIdToContactingNeighbours(data_buffer, int(id));

}

void AnalyticSphericParticle::ComputeBallToRigidFaceContactForce(SphericParticle::ParticleDataBuffer & data_buffer,
                                                         array_1d<double, 3>& r_elastic_force,
                                                         array_1d<double, 3>& r_contact_force,
                                                         double& RollingResistance,
                                                         array_1d<double, 3>& rigid_element_force,
                                                         ProcessInfo& r_process_info,                                                        
                                                         int search_control)

{

    SphericParticle::ComputeBallToRigidFaceContactForce(data_buffer, 
                                        r_elastic_force, 
                                        r_contact_force, 
                                        RollingResistance, 
                                        rigid_element_force, 
                                        r_process_info, 
                                        search_control);

    
    //const auto face_id = data_buffer.mNeighbourRigidFaces->Id(); //definir neighbor wall id

    std::vector<DEMWall*>& rNeighbours = this->mNeighbourRigidFaces;
        
    for (unsigned int i = 0; i < rNeighbours.size(); i++) {
        DEMWall* p_wall = rNeighbours[i];        
        if(p_wall == NULL) continue;
        if(p_wall->IsPhantom()){
            p_wall->CheckSide(this);
            continue;
        }    
        int p_wall_id;
        p_wall_id = p_wall->Id();    
        if (IsNewFaceNeighbour(p_wall_id) && mNumberOfCollidingSpheresWithFaces < mMaxCollidingFaceSpheres){
           
            RecordNewFaceImpact(data_buffer);
        }
        PushBackIdToContactingFaceNeighbours(data_buffer, int(p_wall_id));
    }      
}


bool AnalyticSphericParticle::IsNewNeighbour(const int neighbour_id)
{
    const bool already_in_contact = std::find(mContactingNeighbourIds.begin(), mContactingNeighbourIds.end(), neighbour_id) != mContactingNeighbourIds.end();
    return !already_in_contact;
}

 bool AnalyticSphericParticle::IsNewFaceNeighbour(const int p_wall_id)
{
   bool const already_in_contact = std::find(mContactingFaceNeighbourIds.begin(), mContactingFaceNeighbourIds.end(), p_wall_id) != mContactingFaceNeighbourIds.end();
   return !already_in_contact;
} 

void AnalyticSphericParticle::RecordNewImpact(BaseBufferType & data_buffer)
{
    mCollidingIds[mNumberOfCollidingSpheres] = data_buffer.mpOtherParticle->Id();
    mCollidingRadii[mNumberOfCollidingSpheres] = data_buffer.mOtherRadius;
    mCollidingNormalVelocities[mNumberOfCollidingSpheres] = data_buffer.mLocalRelVel[2];
    mCollidingTangentialVelocities[mNumberOfCollidingSpheres] = std::sqrt(data_buffer.mLocalRelVel[0] * data_buffer.mLocalRelVel[0] + data_buffer.mLocalRelVel[1] * data_buffer.mLocalRelVel[1]);
    
    /*
    double mass = GetMass();
    double other_mass = data_buffer.mpOtherParticle->GetMass();

    const array_1d<double, 3>& vel = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    //array_1d<double, 3> local_vel;
    //GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, vel, local_vel);
    //double v_pre = local_vel[2];
    double v_pre = vel[2];
    double other_v_pre = v_pre - data_buffer.mLocalRelVel[2];

    double p1 = mass * v_pre;
    double p2 = other_mass * other_v_pre;
    double q1 = mass * v_pre * v_pre;
    double q2 = other_mass * other_v_pre * other_v_pre; 

    double a = (p1*p2)/other_mass;
    double b = mass/other_mass;
    double c = (q1*q2)/mass;

    double v_post_1 = (2*a + sqrt(4*a*a-4*(1+b)*(a*a/b-c))) / (2*(1+b));
    //double v_post_2 = (2*a - sqrt(4*a*a-4*(1+b)*(a*a/b-c))) / (2*(1+b));
    //criteri per decidir la solucio correcta
    double LinearImpulse = mass * (v_post_1 - v_pre);
    */
    mCollidingLinearImpulse[mNumberOfCollidingSpheres] = 0.0;
    ++mNumberOfCollidingSpheres;
}

 void AnalyticSphericParticle::RecordNewFaceImpact(BaseBufferType & data_buffer)
{

    mCollidingFaceNormalVelocities[mNumberOfCollidingSpheresWithFaces] = data_buffer.mLocalRelVel[2];
    mCollidingFaceTangentialVelocities[mNumberOfCollidingSpheresWithFaces] = std::sqrt(data_buffer.mLocalRelVel[0] * data_buffer.mLocalRelVel[0] + data_buffer.mLocalRelVel[1] * data_buffer.mLocalRelVel[1]);
    ++mNumberOfCollidingSpheresWithFaces;   


    // //char impact_type = data_buffer.mImpactType;  WORK IN PROGRESS
    // if (impact_type == "f")
    // {
    //     mCollidingFaceNormalVelocities[mNumberOfCollidingSpheresWithFaces] = data_buffer.mLocalRelVel[2];
    //     mCollidingFaceTangentialVelocities[mNumberOfCollidingSpheresWithFaces] = std::sqrt(data_buffer.mLocalRelVel[0] * data_buffer.mLocalRelVel[0] + data_buffer.mLocalRelVel[1] * data_buffer.mLocalRelVel[1]);
    //     ++mNumberOfCollidingSpheresWithFaces;        
    // }
    // else if (impact_type == "e"){ // investigate correct components
    //     // mCollidingEdgeNormalVelocities[mNumberOfCollidingSpheresWithEdges] = data_buffer.mLocalRelVel[2];
    //     // mCollidingEdgeVelocitiesAlongEdge[mNumberOfCollidingSpheresWithEdges] = data_buffer.mLocalRelVel[0];
    //     // mCollidingEdgeVelocitiesTransverseToEdge[mNumberOfCollidingSpheresWithEdges] = data_buffer.mLocalRelVel[1];
    //     // ++mNumberOfCollidingSpheresWithEdges; 
    // }
    // else if (impact_type == "v"){

    // }
    // else{

    // } // To-Do kratos error


}  

void AnalyticSphericParticle::FinalizeForceComputation(BaseType::ParticleDataBuffer & data_buffer)
{
    mContactingNeighbourIds = AnalyticSphericParticle::GetPointerToDerivedDataBuffer(data_buffer)->mCurrentContactingNeighbourIds;
    mContactingFaceNeighbourIds = AnalyticSphericParticle::GetPointerToDerivedDataBuffer(data_buffer)->mCurrentContactingFaceNeighbourIds;
    ClearNeighbours(data_buffer);
}

}  // namespace Kratos.
