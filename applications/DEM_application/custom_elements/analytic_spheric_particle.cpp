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
    mNumberOfCollidingSpheres = 0;

    for (unsigned int i = 0; i < 4; ++i){
        mCollidingIds[i] = 0;
        mCollidingRadii[i] = 0.0;
        mCollidingNormalVelocities[i] = 0.0;
        mCollidingTangentialVelocities[i] = 0.0;
    }
}

AnalyticSphericParticle::AnalyticSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry)
    : SphericParticle(NewId, pGeometry)
{
    mNumberOfCollidingSpheres = 0;

    for (unsigned int i = 0; i < 4; ++i){
        mCollidingIds[i] = 0;
        mCollidingRadii[i] = 0.0;
        mCollidingNormalVelocities[i] = 0.0;
        mCollidingTangentialVelocities[i] = 0.0;
    }
}

AnalyticSphericParticle::AnalyticSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : SphericParticle(NewId, pGeometry, pProperties)
{
    mNumberOfCollidingSpheres = 0;

    for (unsigned int i = 0; i < 4; ++i){
        mCollidingIds[i] = 0;
        mCollidingRadii[i] = 0.0;
        mCollidingNormalVelocities[i] = 0.0;
        mCollidingTangentialVelocities[i] = 0.0;
    }
}

AnalyticSphericParticle::AnalyticSphericParticle(IndexType NewId, NodesArrayType const& ThisNodes)
    : SphericParticle(NewId, ThisNodes)
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

void AnalyticSphericParticle::CalculateRelativePositions(ParticleDataBuffer & data_buffer)
{
    SphericParticle::CalculateRelativePositions(data_buffer);
    const auto id = data_buffer.mpOtherParticle->Id();

    if (IsNewNeighbour(id)){
        RecordNewImpact(data_buffer);
    }

    data_buffer.mCurrentNeighbourIds.push_back(id);
}

bool AnalyticSphericParticle::IsNewNeighbour(const int nighbour_id)
{
    for (int i = 0; i < int(mContactingNeighbourIds.size()); ++i){
       if (mContactingNeighbourIds[i] == int(nighbour_id)){
           return false;
       }
    }

    return true;
}

void AnalyticSphericParticle::RecordNewImpact(ParticleDataBuffer & data_buffer)
{
    mCollidingIds[mNumberOfCollidingSpheres] = data_buffer.mpOtherParticle->Id();
    mCollidingRadii[mNumberOfCollidingSpheres] = data_buffer.mOtherRadius;
    mCollidingNormalVelocities[mNumberOfCollidingSpheres] = data_buffer.mLocalRelVel[0];
    mCollidingTangentialVelocities[mNumberOfCollidingSpheres] = std::sqrt(data_buffer.mLocalRelVel[0] * data_buffer.mLocalRelVel[0] + data_buffer.mLocalRelVel[1] * data_buffer.mLocalRelVel[1]);
    ++mNumberOfCollidingSpheres;
}

void AnalyticSphericParticle::FinalizeForceComputation(ParticleDataBuffer & data_buffer)
{
    mNumberOfCollidingSpheres = 0;
    mContactingNeighbourIds = data_buffer.mCurrentNeighbourIds;
}

}  // namespace Kratos.
