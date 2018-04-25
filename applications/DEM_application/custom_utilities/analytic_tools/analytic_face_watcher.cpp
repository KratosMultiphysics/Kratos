//   $Author: Guillermo Casas

// Project includes

// System includes
#include <limits>
#include <iostream>
#include <iomanip>

// External includes
#ifdef _OPENMP
#include <omp.h>
#endif

// Project includes
#include "analytic_face_watcher.h"
#include "DEM_application.h"

namespace Kratos
{

typedef ModelPart::ConditionsContainerType::iterator ConditionsIteratorType;
typedef AnalyticRigidFace3D AnalyticFace;

void AnalyticFaceWatcher::ClearData()
{
    mSetOfIds.clear();
    mVectorOfTimeStepDatabases.clear();
    mMapOfFaceHistoryDatabases.clear();
}

void AnalyticFaceWatcher::MakeMeasurements(ModelPart& analytic_model_part)
{
    const double current_time = analytic_model_part.GetProcessInfo()[TIME];
    CrossingsTimeStepDataBase time_step_database(current_time);

    for (ConditionsIteratorType i_cond = analytic_model_part.ConditionsBegin(); i_cond != analytic_model_part.ConditionsEnd(); ++i_cond){
        AnalyticFace& face = dynamic_cast<Kratos::AnalyticRigidFace3D&>(*(*(i_cond.base())));
        const int n_crossings = abs(face.GetNumberThroughput());

        if (n_crossings){
            KRATOS_WATCH("CROSSINGSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS")
            const int id = int(i_cond->Id());
            FaceHistoryDatabase& face_database = GetFaceDataBase(id);
            std::vector<int> colliding_ids = face.GetSignedCollidingIds();
            std::vector<double> particle_masses = face.GetMasses();
            std::vector<double> colliding_normal_vel = face.GetCollidingNormalRelativeVelocity();
            std::vector<double> colliding_tangential_vel = face.GetCollidingTangentialRelativeVelocity();

            for (unsigned int i = 0; i < colliding_ids.size(); ++i){
                time_step_database.PushBackCrossings(id, colliding_ids[i], particle_masses[i], colliding_normal_vel[i], colliding_tangential_vel[i]);
                face_database.PushBackCrossings(current_time, colliding_ids[i], particle_masses[i], colliding_normal_vel[i], colliding_tangential_vel[i]);
            }
        }
    }

    mVectorOfTimeStepDatabases.push_back(time_step_database);
}

void AnalyticFaceWatcher::GetFaceData(int id,
                                      std::list<double> times,
                                      std::list<int> neighbour_ids,
                                      std::list<double> masses,
                                      std::list<double> normal_relative_vel,
                                      std::list<double> tangential_relative_vel)
{
    KRATOS_WATCH("AnalyticFaceWatcher::GetFaceData")
    mMapOfFaceHistoryDatabases[id].FillUpPythonLists(times, neighbour_ids, masses, normal_relative_vel, tangential_relative_vel);
}

void AnalyticFaceWatcher::GetAllFacesData(ModelPart& analytic_model_part,
                                          std::list<double> times,
                                          std::list<int> neighbour_ids,
                                          std::list<double> masses,
                                          std::list<double> normal_relative_vel,
                                          std::list<double> tangential_relative_vel)
{
    times.clear();
    neighbour_ids.clear();
    masses.clear();
    normal_relative_vel.clear();
    tangential_relative_vel.clear();
    
    KRATOS_WATCH("AnalyticFaceWatcher::GetAllFacesData outside loop")

    for (ConditionsIteratorType i_cond = analytic_model_part.ConditionsBegin(); i_cond != analytic_model_part.ConditionsEnd(); ++i_cond){
        std::list<double> times_i;
        std::list<int> neighbour_ids_i;
        std::list<double> masses_i;
        std::list<double> normal_relative_vel_i;
        std::list<double> tangential_relative_vel_i;
        
        const int id = int(i_cond->Id());
                
        GetFaceData(id, times_i, neighbour_ids_i, masses_i, normal_relative_vel_i, tangential_relative_vel_i);

        times.insert(times.end(), times_i.begin(), times_i.end());
        neighbour_ids.insert(neighbour_ids.end(), neighbour_ids_i.begin(), neighbour_ids_i.end());
        masses.insert(masses.end(), masses_i.begin(), masses_i.end());
        normal_relative_vel.insert(normal_relative_vel.end(), normal_relative_vel_i.begin(), normal_relative_vel_i.end());
        tangential_relative_vel.insert(tangential_relative_vel_i.end(), tangential_relative_vel_i.begin(), tangential_relative_vel_i.end());
        KRATOS_WATCH("AnalyticFaceWatcher::GetAllFacesData inside loop")
        KRATOS_WATCH(times.size())
        KRATOS_WATCH(neighbour_ids.size())
        KRATOS_WATCH(masses.size())
        KRATOS_WATCH(normal_relative_vel.size())
        KRATOS_WATCH(tangential_relative_vel.size())
    }

}

void AnalyticFaceWatcher::GetTimeStepsData(std::list<int> ids,
                                           std::list<int> neighbour_ids,
                                           std::list<double> masses,
                                           std::list<double> normal_relative_vel,
                                           std::list<double> tangential_relative_vel)
{
    ids.clear();
    neighbour_ids.clear();
    masses.clear();
    normal_relative_vel.clear();
    tangential_relative_vel.clear();

    const int n_time_steps = mVectorOfTimeStepDatabases.size();

    for (int i = 0; i < n_time_steps; ++i){
        std::list<int> ids_i;
        std::list<int> neighbour_ids_i;
        std::list<double> masses_i;
        std::list<double> normal_relative_vel_i;
        std::list<double> tangential_relative_vel_i;

        mVectorOfTimeStepDatabases[i].FillUpPythonLists(ids_i, neighbour_ids_i, masses_i, normal_relative_vel_i, tangential_relative_vel_i);
        
        ids.insert(ids.end(), ids_i.begin(), ids_i.end());
        neighbour_ids.insert(neighbour_ids.end(), neighbour_ids_i.begin(), neighbour_ids_i.end());
        masses.insert(masses.end(), masses_i.begin(), masses_i.end());
        normal_relative_vel.insert(normal_relative_vel.end(), normal_relative_vel_i.begin(), normal_relative_vel_i.end());
        tangential_relative_vel.insert(tangential_relative_vel.end(), tangential_relative_vel_i.begin(), tangential_relative_vel_i.end());
        KRATOS_WATCH("AnalyticFaceWatcher::GetTimeStepsData")
        KRATOS_WATCH(ids.size())
        KRATOS_WATCH(neighbour_ids.size())
        KRATOS_WATCH(masses.size())
        KRATOS_WATCH(normal_relative_vel.size())
        KRATOS_WATCH(tangential_relative_vel.size())
    }
}

void AnalyticFaceWatcher::GetTotalFlux(std::list<double> &times, std::list<int> &n_particles, std::list<double> &mass)
{
    times.clear();
    n_particles.clear();
    mass.clear();

    const int n_time_steps = mVectorOfTimeStepDatabases.size();

    for (int i = 0; i < n_time_steps; ++i){
        times.push_back(mVectorOfTimeStepDatabases[i].GetTime());
        n_particles.push_back(mVectorOfTimeStepDatabases[i].GetTotalThroughput());
        mass.push_back(mVectorOfTimeStepDatabases[i].GetTotalMassThroughput());
    }
    KRATOS_WATCH("AnalyticFaceWatcher::GetTotalFlux")
    std::list<int>::iterator it_int;
    std::list<double>::iterator it_double;
    std::list<double>::iterator it_double_2;
    for (it_double = times.begin(); it_double != times.end(); it_double++) {
        KRATOS_WATCH(*it_double)
    }
    for (it_int = n_particles.begin(); it_int != n_particles.end(); it_int++) {
        KRATOS_WATCH(*it_int)
    }
    for (it_double_2 = mass.begin(); it_double_2 != mass.end(); it_double_2++) {
        KRATOS_WATCH(*it_double_2)
    }
}

AnalyticFaceWatcher::FaceHistoryDatabase& AnalyticFaceWatcher::GetFaceDataBase(int id)
{
    if (mMapOfFaceHistoryDatabases.find(id) == mMapOfFaceHistoryDatabases.end()){
        AnalyticFaceWatcher::FaceHistoryDatabase new_face_database(id);
        mMapOfFaceHistoryDatabases[id] = new_face_database;
    }
    KRATOS_WATCH("AnalyticFaceWatcher::GetTotalFlux")
    return mMapOfFaceHistoryDatabases[id];
}

/// Turn back information as a string.
std::string AnalyticFaceWatcher::Info() const {
        return "AnalyticFaceWatcher";
}

/// Print information about this object.
void AnalyticFaceWatcher::PrintInfo(std::ostream& rOStream) const {}

/// Print object's data.
void AnalyticFaceWatcher::PrintData(std::ostream& rOStream) const {}


} // namespace Kratos
