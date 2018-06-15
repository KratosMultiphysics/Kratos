//   $Author: Guillermo Casas
#ifndef ANALYTIC_FACE_WATCHER_H
#define ANALYTIC_FACE_WATCHER_H

#include <pybind11/pybind11.h>

// System includes

#include <limits>
#include <iostream>
#include <iomanip>
#include <list>

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
class AnalyticFaceWatcher {

public:

KRATOS_CLASS_POINTER_DEFINITION(AnalyticFaceWatcher);

/// Default constructor

AnalyticFaceWatcher(ModelPart& model_part):mrModelPart(model_part){}

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
        mMass += mass*Sign(id2);
        mRelVelNormalxMass += mass*normal_vel;
        mRelVelTangentialxMass += mass*tang_vel;
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

    double GetRelVelNormalxMass()
    {
        return mRelVelNormalxMass;
    }

    double GetRelVelTangentialxMass()
    {
        return mRelVelTangentialxMass;
    }


    void FillUpPythonLists(pybind11::list& ids,
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

        for (int i = 0; i < mNCrossings; ++i){
            //ids.push_back(mId1[i]);
            ids.append(mId1[i]);
            neighbour_ids.append(mId2[i]);
            masses.append(mMasses[i]);
            normal_relative_vel.append(mRelVelNormal[i]);
            tangential_relative_vel.append(mRelVelTangential[i]);
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
        double mRelVelNormalxMass;
        double mRelVelTangentialxMass;

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

    void FillUpPythonLists(pybind11::list& times,
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

        for (int i = 0; i < mNCrossings; ++i){
            //times.push_back(mTimes[i]);
            times.append(mTimes[i]);
            neighbour_ids.append(mId2[i]);
            masses.append(mMasses[i]);
            normal_relative_vel.append(mRelVelNormal[i]);
            tangential_relative_vel.append(mRelVelTangential[i]);
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
        double mRelVelNormalxMass;
        double mRelVelTangentialxMass;
};

void ClearData();

void GetFaceData(int id,
                 pybind11::list times,
                 pybind11::list neighbour_ids,
                 pybind11::list masses,
                 pybind11::list normal_relative_vel,
                 pybind11::list tangential_relative_vel);

void GetAllFacesData(ModelPart& analytic_model_part,
                     pybind11::list& times,
                     pybind11::list& neighbour_ids,
                     pybind11::list& masses,
                     pybind11::list& normal_relative_vel,
                     pybind11::list& tangential_relative_vel);

void GetTimeStepsData(pybind11::list& ids,
                      pybind11::list& neighbour_ids,
                      pybind11::list& masses,
                      pybind11::list& normal_relative_vel,
                      pybind11::list& tangential_relative_vel);

void GetTotalFlux(pybind11::list &times,
                  pybind11::list &n_particles,
                  pybind11::list &mass,
                  pybind11::list &normal_relative_vel,
                  pybind11::list &tangential_relative_vel);

virtual void MakeMeasurements();

virtual FaceHistoryDatabase& GetFaceDataBase(int id);

/// Turn back information as a string
virtual std::string Info() const;

/// Print information about this object
virtual void PrintInfo(std::ostream& rOStream) const;

/// Print object's data
virtual void PrintData(std::ostream& rOStream) const;


private:

ModelPart& mrModelPart;
std::set<int> mSetOfIds;
std::vector<CrossingsTimeStepDataBase> mVectorOfTimeStepDatabases;
std::map<int, FaceHistoryDatabase> mMapOfFaceHistoryDatabases;

/// Assignment operator
AnalyticFaceWatcher & operator=(AnalyticFaceWatcher const& rOther);

}; // Class AnalyticFaceWatcher

} // namespace Kratos.

#endif // ANALYTIC_FACE_WATCHER_H
