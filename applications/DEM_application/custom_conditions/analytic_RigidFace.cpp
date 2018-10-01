// Author: Guillermo Casas (gcasas@cimne.upc.edu)

// Project includes
#include "custom_conditions/analytic_RigidFace.h"
#include "../custom_elements/spheric_particle.h"

namespace Kratos {

// Constructor

AnalyticRigidFace3D::AnalyticRigidFace3D() : mNumberThroughput(0){}

// Constructor

AnalyticRigidFace3D::AnalyticRigidFace3D(IndexType NewId, GeometryType::Pointer pGeometry) : RigidFace3D(NewId, pGeometry), mNumberThroughput(0){}

// Constructor

AnalyticRigidFace3D::AnalyticRigidFace3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
           : RigidFace3D(NewId, pGeometry, pProperties), mNumberThroughput(0){}


//***********************************************************************************
//***********************************************************************************

Condition::Pointer AnalyticRigidFace3D::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new AnalyticRigidFace3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

// Destructor

AnalyticRigidFace3D::~AnalyticRigidFace3D(){}

int AnalyticRigidFace3D::CheckSide(SphericParticle* p_particle)
{
    const int side_sign = BaseType::CheckSide(p_particle);

    const int signed_id = int(p_particle->Id()) * side_sign;

    // the particle just changed side if it can be found in the old list with the opposite sign:
    const bool just_changed_side = AnalyticRigidFace3D::IsInside(- signed_id, mOldContactingNeighbourSignedIds);
    #pragma omp critical
    {
        mContactingNeighbourSignedIds.push_back(signed_id);
        if (just_changed_side){
            const bool is_a_crosser = CheckProjectionFallsInside(p_particle);

            if (is_a_crosser || true){
                mNumberThroughput += side_sign;
                mCrossers.push_back(signed_id);
                mMasses.push_back(p_particle->GetMass());
                array_1d<double, 3> normal;
                CalculateNormal(normal);
                array_1d<double, 3> particle_vel = p_particle->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
                const double normal_vel_component = DEM_INNER_PRODUCT_3(particle_vel, normal);
                mCollidingNormalVelocities.push_back(normal_vel_component);
                noalias(particle_vel) += - normal_vel_component * normal;
                mCollidingTangentialVelocities.push_back(DEM_MODULUS_3(particle_vel));
            }
        }
    }

    return signed_id;
}

// the signed number of crossings
int AnalyticRigidFace3D::GetNumberThroughput()
{
    return mNumberThroughput;
}

int AnalyticRigidFace3D::AreThereNewCrossings()
{
    return int(mCrossers.size());
}

std::vector<int> AnalyticRigidFace3D::GetSignedCollidingIds()
{
    return mCrossers;
}

std::vector<double> AnalyticRigidFace3D::GetCollidingNormalRelativeVelocity()
{
    return mCollidingNormalVelocities;
}

std::vector<double> AnalyticRigidFace3D::GetCollidingTangentialRelativeVelocity()
{
    return mCollidingTangentialVelocities;
}

std::vector<double> AnalyticRigidFace3D::GetMasses()
{
    return mMasses;
}

void AnalyticRigidFace3D::InitializeSolutionStep(ProcessInfo& r_process_info)
{
    RigidFace3D::InitializeSolutionStep(r_process_info);
    mOldContactingNeighbourSignedIds.swap(mContactingNeighbourSignedIds);
    mContactingNeighbourSignedIds.clear();
    mCrossers.clear();
    mMasses.clear();
    mCollidingNormalVelocities.clear();
    mCollidingTangentialVelocities.clear();
    mNumberThroughput = 0;
}


} // Namespace Kratos.
