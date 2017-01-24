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

}  // namespace Kratos.
