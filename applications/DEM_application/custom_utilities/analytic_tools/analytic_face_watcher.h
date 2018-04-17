//   $Author: Guillermo Casas
#ifndef ANALYTIC_FACE_WATCHER_H
#define ANALYTIC_FACE_WATCHER_H

// System includes

#include <limits>
#include <iostream>
#include <iomanip>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "../../custom_conditions/RigidFace.h"

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

namespace Kratos
{
class KRATOS_API(DEM_APPLICATION) AnalyticFaceWatcher {

public:

KRATOS_CLASS_POINTER_DEFINITION(AnalyticFaceWatcher);

/// Default constructor

AnalyticFaceWatcher(){}

/// Destructor

virtual ~AnalyticFaceWatcher(){}

template<typename T>
static inline int Sign(T x)
{
    return (T(0) < x) - (x < T(0));
}

class CrossingsTimeStepDataBase  // It holds the historical information gathered in a single time step
{
    public:

    CrossingsTimeStepDataBase(const double time) : mNCrossings(0), mNSignedCrossings(0), mTime(time), mMass(0.0){}
    ~CrossingsTimeStepDataBase(){}

    int GetNumberOfCrossings()
    {
        return mNCrossings;
    }

    void PushBackCrossings(const int id1, const int id2, const double mass, const double normal_vel, const double tang_vel)
    {
        ++mNCrossings;
        mNSignedCrossings += Sign(id2);
        mMass += mass;
        mMasses.push_back(mass);
        mId1.push_back(id1);
        mId2.push_back(std::abs(id2));
        mRelVelNormal.push_back(normal_vel);
        mRelVelTangential.push_back(tang_vel);
    }    

    int GetTotalThroughput()
    {
        return mNSignedCrossings;
    }

    double GetTotalMassThroughput()
    {
        return mMass;
    }

    double GetTime()
    {
        return mTime;
    }

    void FillUpPythonLists(std::list<int>& ids,
                           std::list<int>& neighbour_ids,
                           std::list<double>& masses,
                           std::list<double>& normal_relative_vel,
                           std::list<double>& tangential_relative_vel)
    {
        ids.clear();
        neighbour_ids.clear();
        masses.clear();
        normal_relative_vel.clear();
        tangential_relative_vel.clear();

        for (int i = 0; i < mNCrossings; ++i){
            ids.push_back(mId1[i]);
            neighbour_ids.push_back(mId2[i]);
            masses.push_back(mMasses[i]);
            normal_relative_vel.push_back(mRelVelNormal[i]);
            tangential_relative_vel.push_back(mRelVelTangential[i]);
        }
    }

    private:

        int mNCrossings;
        int mNSignedCrossings;
        double mTime;
        double mMass;
        std::vector<double> mMasses;
        std::vector<int> mId1;
        std::vector<int> mId2;
        std::vector<double> mRelVelNormal;
        std::vector<double> mRelVelTangential;
    };

class FaceHistoryDatabase // It holds the historical information gathered for a single face
    {
    public:

    FaceHistoryDatabase(): mNCrossings(0), mNSignedCrossings(0), /*mId(0),*/ mMass(0.0){}
    FaceHistoryDatabase(const int id): mNCrossings(0), mNSignedCrossings(0), /*mId(id),*/ mMass(0.0){}
    ~FaceHistoryDatabase(){}

    void PushBackCrossings(const double time, const int id2, const double mass, const double normal_vel, const double tang_vel)
    {
        ++mNCrossings;
        mNSignedCrossings += Sign(id2);
        mMass += mass * Sign(normal_vel);
        mTimes.push_back(time);
        mId2.push_back(std::abs(id2));
        mMasses.push_back(mass * Sign(normal_vel));
        mRelVelNormal.push_back(normal_vel);
        mRelVelTangential.push_back(tang_vel);
    }

    int GetTotalThroughput()
    {
        return mNSignedCrossings;
    }

    double GetTotalMassThroughput()
    {
        return mMass;
    }

    void FillUpPythonLists(std::list<double>& times,
                           std::list<int>& neighbour_ids,
                           std::list<double>& masses,
                           std::list<double>& normal_relative_vel,
                           std::list<double>& tangential_relative_vel)
    {
        times.clear();
        neighbour_ids.clear();
        masses.clear();
        normal_relative_vel.clear();
        tangential_relative_vel.clear();

        for (int i = 0; i < mNCrossings; ++i){
            times.push_back(mTimes[i]);
            neighbour_ids.push_back(mId2[i]);
            masses.push_back(mMasses[i]);
            normal_relative_vel.push_back(mRelVelNormal[i]);
            tangential_relative_vel.push_back(mRelVelTangential[i]);
        }
    }

    private:

        int mNCrossings;
        int mNSignedCrossings;
        //int mId;        
        double mMass;
        std::vector<double> mTimes;
        std::vector<double> mMasses;
        std::vector<int> mId2;
        std::vector<double> mRelVelNormal;
        std::vector<double> mRelVelTangential;
};

void ClearData();

void GetFaceData(int id,
                 std::list<double> times,
                 std::list<int> neighbour_ids,
                 std::list<double> masses,
                 std::list<double> normal_relative_vel,
                 std::list<double> tangential_relative_vel);

void GetAllFacesData(ModelPart& analytic_model_part,
                     std::list<double> times,
                     std::list<int> neighbour_ids,
                     std::list<double> masses,
                     std::list<double> normal_relative_vel,
                     std::list<double> tangential_relative_vel);

void GetTimeStepsData(std::list<int> ids,
                      std::list<int> neighbour_ids,
                      std::list<double> masses,
                      std::list<double> normal_relative_vel,
                      std::list<double> tangential_relative_vel);

void GetTotalFlux(std::list<double> &times, std::list<int> &n_particles, std::list<double> &mass);

virtual void MakeMeasurements(ModelPart& analytic_model_part);

virtual FaceHistoryDatabase& GetFaceDataBase(int id);

/// Turn back information as a string
virtual std::string Info() const;

/// Print information about this object
virtual void PrintInfo(std::ostream& rOStream) const;

/// Print object's data
virtual void PrintData(std::ostream& rOStream) const;


private:

std::set<int> mSetOfIds;
std::vector<CrossingsTimeStepDataBase> mVectorOfTimeStepDatabases;
std::map<int, FaceHistoryDatabase> mMapOfFaceHistoryDatabases;

/// Assignment operator
AnalyticFaceWatcher & operator=(AnalyticFaceWatcher const& rOther);

}; // Class AnalyticFaceWatcher

} // namespace Kratos.

#endif // ANALYTIC_FACE_WATCHER_H
