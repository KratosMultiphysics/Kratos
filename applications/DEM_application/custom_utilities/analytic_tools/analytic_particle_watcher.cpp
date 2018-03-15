//   $Main author: Guillermo Casas
//

// Project includes

// System includes
#include <limits>
#include <iostream>
#include <iomanip>
#include <list>

// External includes
#ifdef _OPENMP
#include <omp.h>
#endif

// Project includes
#include "analytic_particle_watcher.h"

namespace Kratos
{

typedef ModelPart::ElementsContainerType::iterator ElementsIteratorType;
typedef Kratos::AnalyticSphericParticle AnalyticParticle;

void AnalyticParticleWatcher::MakeMeasurements(ModelPart& analytic_model_part)
{
    const double current_time = analytic_model_part.GetProcessInfo()[TIME];
    InterParticleImpactDataOfAParticle time_step_database(current_time);
    FaceParticleImpactDataOfAParticle face_time_step_database(current_time);

    for (ElementsIteratorType i_elem = analytic_model_part.ElementsBegin(); i_elem != analytic_model_part.ElementsEnd(); ++i_elem){
        AnalyticParticle& particle = dynamic_cast<Kratos::AnalyticSphericParticle&>(*(*(i_elem.base())));

        const int n_collisions = particle.GetNumberOfCollisions();
        
        if (n_collisions){
            const int id = int(i_elem->Id());
            InterParticleImpactDataOfATimeStep& particle_database = GetParticleDataBase(id);
            const array_1d<int, 4> &colliding_ids = particle.GetCollidingIds();
            const array_1d<double, 4> &colliding_normal_vel = particle.GetCollidingNormalRelativeVelocity();
            const array_1d<double, 4> &colliding_tangential_vel = particle.GetCollidingTangentialRelativeVelocity();
            const array_1d<double, 4> &colliding_linear_impulse = particle.GetCollidingLinearImpulse();

            /*
            array_1d<double, 4> colliding_normal_vel;
            array_1d<double, 4> colliding_tangential_vel;
            array_1d<double, 4> colliding_linear_impulse;
            particle.GetCollidingIds();
            particle.GetCollidingNormalRelativeVelocity(colliding_normal_vel);
            particle.GetCollidingTangentialRelativeVelocity(colliding_tangential_vel);
            particle.GetCollidingLinearImpulse(colliding_linear_impulse);
            */

            for (int i = 0; i < n_collisions; ++i){
                time_step_database.PushBackImpacts(id, colliding_ids[i], colliding_normal_vel[i], colliding_tangential_vel[i]);
                particle_database.PushBackImpacts(current_time, colliding_ids[i], colliding_normal_vel[i], colliding_tangential_vel[i], colliding_linear_impulse[i]);
            }
        }

        const int n_collisions_with_walls = particle.GetNumberOfCollisionsWithFaces();
        if (n_collisions_with_walls){           
            const int id = int(i_elem->Id());
            FaceParticleImpactDataOfATimeStep& flat_wall_particle_database = GetParticleFaceDataBase(id);
            const array_1d<int, 4> &colliding_ids_with_walls = particle.GetCollidingFaceIds();
            const array_1d<double, 4> &colliding_normal_vel = particle.GetCollidingFaceNormalRelativeVelocity();
            const array_1d<double, 4> &colliding_tangential_vel = particle.GetCollidingFaceTangentialRelativeVelocity();

            /*
            array_1d<int, 4> colliding_ids_with_walls;
            array_1d<double, 4> colliding_normal_vel;
            array_1d<double, 4> colliding_tangential_vel;
            particle.GetCollidingFaceIds(colliding_ids_with_walls);
            particle.GetCollidingFaceNormalRelativeVelocity(colliding_normal_vel);
            particle.GetCollidingFaceTangentialRelativeVelocity(colliding_tangential_vel);
            */

            for (int i = 0; i < n_collisions_with_walls; ++i){
                face_time_step_database.PushBackImpacts(id, colliding_ids_with_walls[i], colliding_normal_vel[i], colliding_tangential_vel[i]);
                flat_wall_particle_database.PushBackImpacts(current_time, colliding_ids_with_walls[i], colliding_normal_vel[i], colliding_tangential_vel[i]);
            }
        }
    }

    mInterParticleImpactDataOfAllParticles.push_back(time_step_database);
    mFaceParticleImpactDataOfAllParticles.push_back(face_time_step_database);
    if (time_step_database.GetNumberOfImpacts()){}
}

void AnalyticParticleWatcher::SetNodalMaxImpactVelocities(ModelPart& analytic_model_part)
{
    for (ElementsIteratorType i_elem = analytic_model_part.ElementsBegin(); i_elem != analytic_model_part.ElementsEnd(); ++i_elem){
        AnalyticParticle& particle = dynamic_cast<Kratos::AnalyticSphericParticle&>(*(*(i_elem.base())));

        const int id = int(i_elem->Id());
        InterParticleImpactDataOfATimeStep& particle_database = GetParticleDataBase(id);
        double db_normal_impact_velocity = 0.0;
        double db_tangential_impact_velocity = 0.0;
        particle_database.GetMaxCollidingSpeedFromDatabase(db_normal_impact_velocity, db_tangential_impact_velocity);

        // get current nodal values
        double& current_max_normal_velocity = particle.GetGeometry()[0].FastGetSolutionStepValue(NORMAL_IMPACT_VELOCITY);
        double& current_max_tangential_velocity = particle.GetGeometry()[0].FastGetSolutionStepValue(TANGENTIAL_IMPACT_VELOCITY);

        // choose max between current and database
        current_max_normal_velocity = std::max(current_max_normal_velocity, db_normal_impact_velocity);
        current_max_tangential_velocity = std::max(current_max_tangential_velocity, db_tangential_impact_velocity);
    }

    //double avg_normal_impact_velocity = normal_velocity_sum/number_of_elements; // retrieve average impact velocity
    //double avg_tangential_impact_velocity = tangential_velocity_sum/number_of_elements;
}


void AnalyticParticleWatcher::SetNodalMaxFaceImpactVelocities(ModelPart& analytic_model_part)
{
    for (ElementsIteratorType i_elem = analytic_model_part.ElementsBegin(); i_elem != analytic_model_part.ElementsEnd(); ++i_elem){
        AnalyticParticle& particle = dynamic_cast<Kratos::AnalyticSphericParticle&>(*(*(i_elem.base())));

        const int id = int(i_elem->Id());
        FaceParticleImpactDataOfATimeStep& particle_database = GetParticleFaceDataBase(id);
        double db_normal_impact_velocity = 0.0;
        double db_tangential_impact_velocity = 0.0;
        particle_database.GetMaxCollidingSpeedFromDatabase(db_normal_impact_velocity, db_tangential_impact_velocity);

        // get current nodal values
        double& current_max_normal_velocity = particle.GetGeometry()[0].FastGetSolutionStepValue(FACE_NORMAL_IMPACT_VELOCITY);
        double& current_max_tangential_velocity = particle.GetGeometry()[0].FastGetSolutionStepValue(FACE_TANGENTIAL_IMPACT_VELOCITY);

        // choose max between current and database
        current_max_normal_velocity = std::max(current_max_normal_velocity, db_normal_impact_velocity);
        current_max_tangential_velocity = std::max(current_max_tangential_velocity, db_tangential_impact_velocity);
    }
}


void AnalyticParticleWatcher::SetNodalMaxLinearImpulse(ModelPart& analytic_model_part)
{
    for (ElementsIteratorType i_elem = analytic_model_part.ElementsBegin(); i_elem != analytic_model_part.ElementsEnd(); ++i_elem){
        AnalyticParticle& particle = dynamic_cast<Kratos::AnalyticSphericParticle&>(*(*(i_elem.base())));

        const int id = int(i_elem->Id());
        InterParticleImpactDataOfATimeStep& particle_database = GetParticleDataBase(id);
        double db_linear_impulse = 0.0;
        particle_database.GetMaxLinearImpulseFromDatabase(db_linear_impulse);

        // get current nodal values
        double& current_max_linear_impulse = particle.GetGeometry()[0].FastGetSolutionStepValue(LINEAR_IMPULSE);

        // choose max between current and database
        current_max_linear_impulse = std::max(current_max_linear_impulse, db_linear_impulse);
    }
}

void AnalyticParticleWatcher::GetParticleData(int id,
                                              std::list<double> times,
                                              std::list<int> neighbour_ids,
                                              std::list<double> normal_relative_vel,
                                              std::list<double> tangential_relative_vel)
{
    mInterParticleImpactDataOfAllTimeSteps[id].FillUpPythonLists(times, neighbour_ids, normal_relative_vel, tangential_relative_vel);
}

void AnalyticParticleWatcher::GetAllParticlesData(ModelPart& analytic_model_part,
                                                  std::list<double> times,
                                                  std::list<int> neighbour_ids,
                                                  std::list<double> normal_relative_vel,
                                                  std::list<double> tangential_relative_vel)
{
    times.clear();
    neighbour_ids.clear();
    normal_relative_vel.clear();
    tangential_relative_vel.clear();

    for (ElementsIteratorType i_elem = analytic_model_part.ElementsBegin(); i_elem != analytic_model_part.ElementsEnd(); ++i_elem){
        std::list<double> times_i;
        std::list<int> neighbour_ids_i;
        std::list<double> normal_relative_vel_i;
        std::list<double> tangential_relative_vel_i;

        const int id = int(i_elem->Id());

        GetParticleData(id, times_i, neighbour_ids_i, normal_relative_vel_i, tangential_relative_vel_i);
        times.insert(times.end(), times_i.begin(), times_i.end());
        neighbour_ids.insert(neighbour_ids.end(), neighbour_ids_i.begin(), neighbour_ids_i.end());
        normal_relative_vel.insert(normal_relative_vel.end(), normal_relative_vel_i.begin(), normal_relative_vel_i.end());
        tangential_relative_vel.insert(tangential_relative_vel.end(), tangential_relative_vel_i.begin(), tangential_relative_vel_i.end());
    }

}

void AnalyticParticleWatcher::GetTimeStepsData(std::list<int> ids,
                                               std::list<int> neighbour_ids,
                                               std::list<double> normal_relative_vel,
                                               std::list<double> tangential_relative_vel)
{
    ids.clear();
    neighbour_ids.clear();
    normal_relative_vel.clear();
    tangential_relative_vel.clear();

    const int n_time_steps = mInterParticleImpactDataOfAllParticles.size();

    std::list<int> ids_i;
    std::list<int> neighbour_ids_i;
    std::list<double> normal_relative_vel_i;
    std::list<double> tangential_relative_vel_i;
    
    for (int i = 0; i < n_time_steps; ++i){
        mInterParticleImpactDataOfAllParticles[i].FillUpPythonLists(ids_i, neighbour_ids_i, normal_relative_vel_i, tangential_relative_vel_i);
        
        ids.insert(ids.end(), ids_i.begin(), ids_i.end());
        neighbour_ids.insert(neighbour_ids.end(), neighbour_ids_i.begin(), neighbour_ids_i.end());
        normal_relative_vel.insert(normal_relative_vel.end(), normal_relative_vel_i.begin(), normal_relative_vel_i.end());
        tangential_relative_vel.insert(tangential_relative_vel.end(), tangential_relative_vel_i.begin(), tangential_relative_vel_i.end());
    }
}


AnalyticParticleWatcher::InterParticleImpactDataOfATimeStep& AnalyticParticleWatcher::GetParticleDataBase(int id)
{
    if (mInterParticleImpactDataOfAllTimeSteps.find(id) == mInterParticleImpactDataOfAllTimeSteps.end()){
        AnalyticParticleWatcher::InterParticleImpactDataOfATimeStep new_particle_database(id);
        mInterParticleImpactDataOfAllTimeSteps[id] = new_particle_database;
    }
    return mInterParticleImpactDataOfAllTimeSteps[id];
}

AnalyticParticleWatcher::FaceParticleImpactDataOfATimeStep& AnalyticParticleWatcher::GetParticleFaceDataBase(int id)
{
    if (mFaceParticleImpactDataOfAllTimeSteps.find(id) == mFaceParticleImpactDataOfAllTimeSteps.end()){
        AnalyticParticleWatcher::FaceParticleImpactDataOfATimeStep new_particle_database(id);
        mFaceParticleImpactDataOfAllTimeSteps[id] = new_particle_database;
    }
    return mFaceParticleImpactDataOfAllTimeSteps[id];
}

/// Turn back information as a string.
std::string AnalyticParticleWatcher::Info() const {
        return "AnalyticParticleWatcher";
}

/// Print information about this object.
void AnalyticParticleWatcher::PrintInfo(std::ostream& rOStream) const {}

/// Print object's data.
void AnalyticParticleWatcher::PrintData(std::ostream& rOStream) const {}


} // namespace Kratos
