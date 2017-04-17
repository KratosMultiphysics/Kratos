// Author: Guillermo Casas (gcasas@cimne.upc.edu)

// Project includes
#include "custom_conditions/analytic_RigidFace.h"
#include "../custom_elements/spheric_particle.h"

namespace Kratos {

// Constructor

AnalyticRigidFace3D::AnalyticRigidFace3D() : mNumberOfCrossingSpheres(0){}

// Constructor

AnalyticRigidFace3D::AnalyticRigidFace3D(IndexType NewId, GeometryType::Pointer pGeometry) : RigidFace3D(NewId, pGeometry), mNumberOfCrossingSpheres(0){}

// Constructor

AnalyticRigidFace3D::AnalyticRigidFace3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
           : RigidFace3D(NewId, pGeometry, pProperties), mNumberOfCrossingSpheres(0){}


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


void AnalyticRigidFace3D::ComputeConditionRelativeData(int rigid_neighbour_index,
                                                       SphericParticle* const particle,
                                                       double LocalCoordSystem[3][3],
                                                       double& DistPToB,
                                                       array_1d<double, 4>& Weight,
                                                       array_1d<double, 3>& wall_delta_disp_at_contact_point,
                                                       array_1d<double, 3>& wall_velocity_at_contact_point,
                                                       int& ContactType)
{
    BaseType::ComputeConditionRelativeData(rigid_neighbour_index,
                                           particle,
                                           LocalCoordSystem,
                                           DistPToB,
                                           Weight,
                                           wall_delta_disp_at_contact_point,
                                           wall_velocity_at_contact_point,
                                           ContactType);

}

int AnalyticRigidFace3D::CheckSide(SphericParticle* p_particle)
{
    const int side_sign = BaseType::CheckSide(p_particle);

    const int signed_id = int(p_particle->Id() * side_sign);

    // the particle just changed side if it can be found in the old list with the opposite sign:
    const auto beginning = std::begin(mOldContactingNeighbourSignedIds);
    const auto end       = std::end(mOldContactingNeighbourSignedIds);
    const bool just_changed_side = (end != std::find(beginning, end, - signed_id));

#pragma omp critical
{
    mContactingNeighbourSignedIds.push_back(signed_id);
    if (just_changed_side){
        mAllCrossers.push_back(fabs(signed_id));
        const bool is_a_crosser = CheckProjectionFallsInSide(p_particle);
        if (is_a_crosser || true){
            ++mNumberOfCrossingSpheres;
            mCrossers.push_back(fabs(signed_id));
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

int AnalyticRigidFace3D::GetNumberOfCrossings()
{
    return mNumberOfCrossingSpheres;
}

std::vector<int> AnalyticRigidFace3D::GetIdsOfCrossers()
{
    return mCrossers;
}

std::vector<int> AnalyticRigidFace3D::GetSignedCollidingIds()
{
    return mContactingNeighbourSignedIds;
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
{   RigidFace3D::InitializeSolutionStep(r_process_info);
    mOldContactingNeighbourSignedIds.swap(mContactingNeighbourSignedIds);
    mContactingNeighbourSignedIds.clear();
    mCrossers.clear();
    mMasses.clear();
    mCollidingNormalVelocities.clear();
    mCollidingTangentialVelocities.clear();
    mNumberOfCrossingSpheres = 0;
}


} // Namespace Kratos.
