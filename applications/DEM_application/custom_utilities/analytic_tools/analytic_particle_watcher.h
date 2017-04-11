#ifndef ANALYTIC_PARTICLE_WATCHER_H
#define ANALYTIC_PARTICLE_WATCHER_H

// System includes

#include <limits>
#include <iostream>
#include <iomanip>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "../../custom_elements/analytic_spheric_particle.h"


/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif
#include "boost/python/list.hpp"

namespace Kratos
{
class AnalyticParticleWatcher {

public:   

KRATOS_CLASS_POINTER_DEFINITION(AnalyticParticleWatcher);

/// Default constructor

AnalyticParticleWatcher(){}

/// Destructor

virtual ~AnalyticParticleWatcher(){}


class ImpactsTimeStepDataBase  // It holds the historical information gathered in a single time step
{
    public:

    ImpactsTimeStepDataBase(const double time) : mNImpacts(0), mTime(time){}
    ~ImpactsTimeStepDataBase(){}

    int GetNumberOfImpacts()
    {
        return mNImpacts;
    }

    void PushBackImpacts(const int id1, const int id2, const double normal_vel, const double tang_vel)
    {
        if (ImpactIsNew(id2)){
            ++mNImpacts;
            mId1.push_back(id1);
            mId2.push_back(id2);
            mRelVelNormal.push_back(normal_vel);
            mRelVelTangential.push_back(tang_vel);
        }
    }

    void FillUpPythonLists(boost::python::list& ids,
                           boost::python::list& neighbour_ids,
                           boost::python::list& normal_relative_vel,
                           boost::python::list& tangential_relative_vel)
    {
        for (int i = 0; i < mNImpacts; ++i){
            AnalyticParticleWatcher::ClearList(ids);
            AnalyticParticleWatcher::ClearList(neighbour_ids);
            AnalyticParticleWatcher::ClearList(normal_relative_vel);
            AnalyticParticleWatcher::ClearList(tangential_relative_vel);
            ids.append(mId1[i]);
            neighbour_ids.append(mId2[i]);
            normal_relative_vel.append(mRelVelNormal[i]);
            tangential_relative_vel.append(mRelVelTangential[i]);
        }
    }

    private:

        int mNImpacts;
        double mTime;
        std::vector<int> mId1;
        std::vector<int> mId2;
        std::vector<double> mRelVelNormal;
        std::vector<double> mRelVelTangential;

        bool ImpactIsNew(const int id_2)
        {
            return std::find(mId1.begin(), mId1.end(), id_2) != mId1.end();
        }
    };

class ParticleHistoryDatabase // It holds the historical information gathered for a single particle
    {
        public:

        ParticleHistoryDatabase(): mNImpacts(0), mId(0){}
        ParticleHistoryDatabase(const int id) : mNImpacts(0), mId(id){}
        ~ParticleHistoryDatabase(){}

        void PushBackImpacts(const double time, const int id2, const double normal_vel, const double tang_vel)
        {
            ++mNImpacts;
            mTimes.push_back(time);
            mId2.push_back(id2);
            mRelVelNormal.push_back(normal_vel);
            mRelVelTangential.push_back(tang_vel);
        }

        void FillUpPythonLists(boost::python::list& times,
                               boost::python::list& neighbour_ids,
                               boost::python::list& normal_relative_vel,
                               boost::python::list& tangential_relative_vel)
        {
            for (int i = 0; i < mNImpacts; ++i){
                AnalyticParticleWatcher::ClearList(times);
                AnalyticParticleWatcher::ClearList(neighbour_ids);
                AnalyticParticleWatcher::ClearList(normal_relative_vel);
                AnalyticParticleWatcher::ClearList(tangential_relative_vel);
                times.append(mTimes[i]);
                neighbour_ids.append(mId2[i]);
                normal_relative_vel.append(mRelVelNormal[i]);
                tangential_relative_vel.append(mRelVelTangential[i]);
            }
        }

    private:

        int mNImpacts;
        int mId;
        std::vector<double> mTimes;
        std::vector<int> mId2;
        std::vector<double> mRelVelNormal;
        std::vector<double> mRelVelTangential;
};

static void ClearList(boost::python::list& my_list); // its best to pass empty lists in the first place to avoid this operation

void GetParticleData(int id,
                     boost::python::list times,
                     boost::python::list neighbour_ids,
                     boost::python::list normal_relative_vel,
                     boost::python::list tangential_relative_vel);

void GetAllParticlesData(ModelPart& analytic_model_part,
                         boost::python::list times,
                         boost::python::list neighbour_ids,
                         boost::python::list normal_relative_vel,
                         boost::python::list tangential_relative_vel);

void GetTimeStepsData(boost::python::list ids,
                      boost::python::list neighbour_ids,
                      boost::python::list normal_relative_vel,
                      boost::python::list tangential_relative_vel);

virtual void MakeMeasurements(ModelPart& analytic_model_part);

virtual ParticleHistoryDatabase& GetParticleDataBase(int id);

/// Turn back information as a string
virtual std::string Info() const;

/// Print information about this object
virtual void PrintInfo(std::ostream& rOStream) const;

/// Print object's data
virtual void PrintData(std::ostream& rOStream) const;


private:

std::set<int> mSetOfIds;
std::vector<ImpactsTimeStepDataBase> mVectorOfTimeStepDatabases;
std::map<int, ParticleHistoryDatabase> mMapOfParticleHistoryDatabases;

/// Assignment operator
AnalyticParticleWatcher & operator=(AnalyticParticleWatcher const& rOther);

}; // Class AnalyticParticleWatcher

} // namespace Kratos.

#endif // ANALYTIC_PARTICLE_WATCHER_H
