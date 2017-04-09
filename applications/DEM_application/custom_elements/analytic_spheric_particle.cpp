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
    ClearMemberVariables();
}

AnalyticSphericParticle::AnalyticSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry)
    : SphericParticle(NewId, pGeometry)
{
    ClearMemberVariables();
}

AnalyticSphericParticle::AnalyticSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : SphericParticle(NewId, pGeometry, pProperties)
{
    ClearMemberVariables();
}

AnalyticSphericParticle::AnalyticSphericParticle(IndexType NewId, NodesArrayType const& ThisNodes)
    : SphericParticle(NewId, ThisNodes)
{
    ClearMemberVariables();
}

array_1d<int, 4> AnalyticSphericParticle::GetCollidingIds()
{
    return mCollidingIds;
}


void AnalyticSphericParticle::ClearMemberVariables()
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

void AnalyticSphericParticle::CalculateRelativePositions(BaseBufferType & data_buffer)
{
    SphericParticle::CalculateRelativePositions(data_buffer);
    const auto id = data_buffer.mpOtherParticle->Id();

    if (IsNewNeighbour(id)){
        RecordNewImpact(data_buffer);
        KRATOS_WATCH(id)
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
    mCollidingNormalVelocities[mNumberOfCollidingSpheres] = data_buffer.mLocalRelVel[0];
    mCollidingTangentialVelocities[mNumberOfCollidingSpheres] = std::sqrt(data_buffer.mLocalRelVel[0] * data_buffer.mLocalRelVel[0] + data_buffer.mLocalRelVel[1] * data_buffer.mLocalRelVel[1]);
    ++mNumberOfCollidingSpheres;
}

void AnalyticSphericParticle::FinalizeForceComputation(BaseType::ParticleDataBuffer & data_buffer)
{
    mNumberOfCollidingSpheres = 0;
    mContactingNeighbourIds = AnalyticSphericParticle::GetPointerToDerivedDataBuffer(data_buffer)->mCurrentContactingNeighbourIds;
    ClearNeighbours(data_buffer);
}

}  // namespace Kratos.
