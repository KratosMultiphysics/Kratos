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


int AnalyticSphericParticle::GetNumberOfCollisions()
{
    return mNumberOfCollidingSpheres;
}

void AnalyticSphericParticle::GetCollidingIds(array_1d<int, 4>& colliding_ids)
{
    colliding_ids = mCollidingIds;
}

void AnalyticSphericParticle::GetCollidingNormalRelativeVelocity(array_1d<double, 4>& colliding_normal_vel)
{
    colliding_normal_vel = mCollidingNormalVelocities;
}

void AnalyticSphericParticle::GetCollidingTangentialRelativeVelocity(array_1d<double, 4>& colliding_tangential_vel)
{
    colliding_tangential_vel = mCollidingTangentialVelocities;
}

void AnalyticSphericParticle::ClearImpactMemberVariables()
{
    mNumberOfCollidingSpheres = 0;

    for (unsigned int i = 0; i < 4; ++i){
        mCollidingIds[i] = 0;
        mCollidingRadii[i] = 0.0;
        mCollidingNormalVelocities[i] = 0.0;
        mCollidingTangentialVelocities[i] = 0.0;
    }
}

Element::Pointer AnalyticSphericParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new AnalyticSphericParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

void AnalyticSphericParticle::PushBackIdToContactingNeighbours(BaseBufferType & data_buffer, int id)
{
    GetPointerToDerivedDataBuffer(data_buffer)->mCurrentContactingNeighbourIds.push_back(id);
}

void AnalyticSphericParticle::ClearNeighbours(BaseBufferType & data_buffer)
{
    GetPointerToDerivedDataBuffer(data_buffer)->mCurrentContactingNeighbourIds.clear();
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
                                                                      p_neighbour_element,
                                                                      sliding,
                                                                      LocalCoordSystem,
                                                                      OldLocalCoordSystem,
                                                                      neighbour_elastic_contact_force);

    const auto id = data_buffer.mpOtherParticle->Id();

    if (IsNewNeighbour(id)){
        RecordNewImpact(data_buffer);
    }

    PushBackIdToContactingNeighbours(data_buffer, int(id));

}

bool AnalyticSphericParticle::IsNewNeighbour(const int nighbour_id)
{
    const bool already_in_contact = std::find(mContactingNeighbourIds.begin(), mContactingNeighbourIds.end(), nighbour_id) != mContactingNeighbourIds.end();

    return !already_in_contact;
}

void AnalyticSphericParticle::RecordNewImpact(BaseBufferType & data_buffer)
{
    mCollidingIds[mNumberOfCollidingSpheres] = data_buffer.mpOtherParticle->Id();
    mCollidingRadii[mNumberOfCollidingSpheres] = data_buffer.mOtherRadius;
    mCollidingNormalVelocities[mNumberOfCollidingSpheres] = data_buffer.mLocalRelVel[2];
    mCollidingTangentialVelocities[mNumberOfCollidingSpheres] = std::sqrt(data_buffer.mLocalRelVel[0] * data_buffer.mLocalRelVel[0] + data_buffer.mLocalRelVel[1] * data_buffer.mLocalRelVel[1]);
    ++mNumberOfCollidingSpheres;
}

void AnalyticSphericParticle::FinalizeForceComputation(BaseType::ParticleDataBuffer & data_buffer)
{
    mContactingNeighbourIds = AnalyticSphericParticle::GetPointerToDerivedDataBuffer(data_buffer)->mCurrentContactingNeighbourIds;
    ClearNeighbours(data_buffer);
}

}  // namespace Kratos.
