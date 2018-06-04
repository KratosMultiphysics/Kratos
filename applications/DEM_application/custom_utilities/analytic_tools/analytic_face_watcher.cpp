//   $Author: Guillermo Casas

// Project includes
#include "analytic_face_watcher.h"

// System includes
#include <limits>
#include <iostream>
#include <iomanip>

// External includes
#ifdef _OPENMP
#include <omp.h>
#endif

// Project includes
#include "DEM_application.h"

// Project includes
#include "includes/define.h"


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

void AnalyticFaceWatcher::MakeMeasurements()
{
    const double current_time = mrModelPart.GetProcessInfo()[TIME];
    CrossingsTimeStepDataBase time_step_database(current_time);
    for (ConditionsIteratorType i_cond = mrModelPart.ConditionsBegin(); i_cond != mrModelPart.ConditionsEnd(); ++i_cond){
        AnalyticFace& face = dynamic_cast<Kratos::AnalyticRigidFace3D&>(*(*(i_cond.base())));
        const int n_crossings = face.AreThereNewCrossings();
        if (n_crossings){
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
                                      pybind11::list times,
                                      pybind11::list neighbour_ids,
                                      pybind11::list masses,
                                      pybind11::list normal_relative_vel,
                                      pybind11::list tangential_relative_vel)
{
    mMapOfFaceHistoryDatabases[id].FillUpPythonLists(times, neighbour_ids, masses, normal_relative_vel, tangential_relative_vel);
}

void AnalyticFaceWatcher::GetAllFacesData(ModelPart& analytic_model_part,
                                          pybind11::list& times,
                                          pybind11::list& neighbour_ids,
                                          pybind11::list& masses,
                                          pybind11::list& normal_relative_vel,
                                          pybind11::list& tangential_relative_vel)
{
    times.attr("clear")();
    neighbour_ids.attr("clear")();
    masses.attr("clear")();
    normal_relative_vel.attr("clear")();
    tangential_relative_vel.attr("clear")();

    for (ConditionsIteratorType i_cond = analytic_model_part.ConditionsBegin(); i_cond != analytic_model_part.ConditionsEnd(); ++i_cond){
        pybind11::list times_i;
        pybind11::list neighbour_ids_i;
        pybind11::list masses_i;
        pybind11::list normal_relative_vel_i;
        pybind11::list tangential_relative_vel_i;

        const int id = int(i_cond->Id());

        GetFaceData(id, times_i, neighbour_ids_i, masses_i, normal_relative_vel_i, tangential_relative_vel_i);

        times.append(times_i[id]);
        neighbour_ids.append(neighbour_ids_i[id]);
        masses.append(masses_i[id]);
        normal_relative_vel.append(normal_relative_vel_i[id]);
        tangential_relative_vel.append(tangential_relative_vel_i[id]);
        //times.insert(times.end(), times_i.begin(), times_i.end());
   }

}

void AnalyticFaceWatcher::GetTimeStepsData(pybind11::list& ids,
                                           pybind11::list& neighbour_ids,
                                           pybind11::list& masses,
                                           pybind11::list& normal_relative_vel,
                                           pybind11::list& tangential_relative_vel)
{

    ids.attr("clear")();
    neighbour_ids.attr("clear")();
    masses.attr("clear")();
    normal_relative_vel.attr("clear")();
    tangential_relative_vel.attr("clear")();

    const int n_time_steps = mVectorOfTimeStepDatabases.size();

    for (int i = 0; i < n_time_steps; ++i){
        pybind11::list ids_i;
        pybind11::list neighbour_ids_i;
        pybind11::list masses_i;
        pybind11::list normal_relative_vel_i;
        pybind11::list tangential_relative_vel_i;

        mVectorOfTimeStepDatabases[i].FillUpPythonLists(ids_i,
                                                        neighbour_ids_i,
                                                        masses_i,
                                                        normal_relative_vel_i,
                                                        tangential_relative_vel_i);

        ids.append(ids_i[i]);
        neighbour_ids.append(neighbour_ids_i[i]);
        masses.append(masses_i[i]);
        normal_relative_vel.append(normal_relative_vel_i[i]);
        tangential_relative_vel.append(tangential_relative_vel_i[i]);
        //ids.insert(ids.end(), ids_i.begin(), ids_i.end());
   }
}

void AnalyticFaceWatcher::GetTotalFlux(pybind11::list &times,
                                       pybind11::list &n_particles,
                                       pybind11::list &mass,
                                       pybind11::list &vel_nr,
                                       pybind11::list &vel_tg)
{
    const int n_time_steps = mVectorOfTimeStepDatabases.size();

    for (int i = 0; i < n_time_steps; ++i){
        times.append(mVectorOfTimeStepDatabases[i].GetTime());
        n_particles.append(mVectorOfTimeStepDatabases[i].GetTotalThroughput());
        mass.append(mVectorOfTimeStepDatabases[i].GetTotalMassThroughput());
        vel_nr.append(mVectorOfTimeStepDatabases[i].GetRelVelNormalxMass());
        vel_tg.append(mVectorOfTimeStepDatabases[i].GetRelVelTangentialxMass());
    }
}

AnalyticFaceWatcher::FaceHistoryDatabase& AnalyticFaceWatcher::GetFaceDataBase(int id)
{
    if (mMapOfFaceHistoryDatabases.find(id) == mMapOfFaceHistoryDatabases.end()){
        AnalyticFaceWatcher::FaceHistoryDatabase new_face_database(id);
        mMapOfFaceHistoryDatabases[id] = new_face_database;
    }

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
