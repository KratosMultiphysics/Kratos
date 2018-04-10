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

namespace Kratos
{
class AnalyticParticleWatcher {

public:   

KRATOS_CLASS_POINTER_DEFINITION(AnalyticParticleWatcher);

/// Default constructor

AnalyticParticleWatcher(){}

/// Destructor

virtual ~AnalyticParticleWatcher(){}


class InterParticleImpactDataOfAParticle  // It holds the historical information gathered in a single time step
{
    public:

    InterParticleImpactDataOfAParticle(const double time) : mNImpacts(0)/*, mTime(time)*/{}
    ~InterParticleImpactDataOfAParticle(){}

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

    void FillUpPythonLists(std::list<int>& ids,
                           std::list<int>& neighbour_ids,
                           std::list<double>& normal_relative_vel,
                           std::list<double>& tangential_relative_vel)
    {
        ids.clear();
        neighbour_ids.clear();
        normal_relative_vel.clear();
        tangential_relative_vel.clear();

        for (int i = 0; i < mNImpacts; ++i){
            ids.push_back(mId1[i]);
            neighbour_ids.push_back(mId2[i]);
            normal_relative_vel.push_back(mRelVelNormal[i]);
            tangential_relative_vel.push_back(mRelVelTangential[i]);
        }
    }

    private:

        int mNImpacts;
        /*double mTime;*/
        std::vector<int> mId1;
        std::vector<int> mId2;
        std::vector<double> mRelVelNormal;
        std::vector<double> mRelVelTangential;
       

        bool ImpactIsNew(const int id_2)
        {
            return std::find(mId1.begin(), mId1.end(), id_2) != mId1.end();
        }
    };

class InterParticleImpactDataOfATimeStep // It holds the historical information gathered for a single particle
    {
        public:

        InterParticleImpactDataOfATimeStep(): mNImpacts(0)/*, mId(0)*/{}
        InterParticleImpactDataOfATimeStep(const int id) : mNImpacts(0)/*, mId(id)*/{}
        ~InterParticleImpactDataOfATimeStep(){}

        void PushBackImpacts(const double time, const int id2, const double normal_vel, const double tang_vel, const double linear_impulse)
        {
            ++mNImpacts;
            mTimes.push_back(time);
            mId2.push_back(id2);
            mRelVelNormal.push_back(normal_vel);
            mRelVelTangential.push_back(tang_vel);
            mLinearImpulse.push_back(linear_impulse);
        }

        void FillUpPythonLists(std::list<double>& times,
                               std::list<int32_t>& neighbour_ids,
                               std::list<double>& normal_relative_vel,
                               std::list<double>& tangential_relative_vel)
        {
            times.clear();
            neighbour_ids.clear();
            normal_relative_vel.clear();
            tangential_relative_vel.clear();

            for (int i = 0; i < mNImpacts; ++i){
                times.push_back(mTimes[i]);
                neighbour_ids.push_back(mId2[i]);
                normal_relative_vel.push_back(mRelVelNormal[i]);
                tangential_relative_vel.push_back(mRelVelTangential[i]);
            }
        }

        void GetMaxCollidingSpeedFromDatabase(double& db_normal_impact_velocity, double& db_tangential_impact_velocity){
            if(mRelVelNormal.size()){
                db_normal_impact_velocity = std::abs(*(std::max_element(mRelVelNormal.begin(), mRelVelNormal.end())));
                db_tangential_impact_velocity = std::abs(*(std::max_element(mRelVelTangential.begin(), mRelVelTangential.end())));
            }
            else {
                db_normal_impact_velocity = 0.0;
                db_tangential_impact_velocity = 0.0;
            }

        }


        void GetMaxLinearImpulseFromDatabase(double& db_linear_impulse){
            if(mRelVelNormal.size()){
                db_linear_impulse = std::abs(*(std::max_element(mLinearImpulse.begin(), mLinearImpulse.end())));
                
            }
            else {
                db_linear_impulse = 0.0;
                
            }

        }

    private:

        int mNImpacts;
        //int mId;
        std::vector<double> mTimes;
        std::vector<int> mId2;
        std::vector<double> mRelVelNormal;
        std::vector<double> mRelVelTangential;
        std::vector<double> mLinearImpulse;
};




// FaceParticleImpactDataOfAParticle

class FaceParticleImpactDataOfAParticle  // It holds the historical information gathered in a single time step against flat walls
{
    public:

    FaceParticleImpactDataOfAParticle(const double time) : mNImpacts(0)/*, mTime(time)*/{}
    ~FaceParticleImpactDataOfAParticle(){}

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

    void FillUpPythonLists(std::list<int>& ids,
                           std::list<int>& neighbour_ids,
                           std::list<double>& normal_relative_vel,
                           std::list<double>& tangential_relative_vel)
    {
        ids.clear();
        neighbour_ids.clear();
        normal_relative_vel.clear();
        tangential_relative_vel.clear();

        for (int i = 0; i < mNImpacts; ++i){
            ids.push_back(mId1[i]);
            neighbour_ids.push_back(mId2[i]);
            normal_relative_vel.push_back(mRelVelNormal[i]);
            tangential_relative_vel.push_back(mRelVelTangential[i]);
        }
    }

    private:

        int mNImpacts;
        /*double mTime;*/
        std::vector<int> mId1;
        std::vector<int> mId2;
        std::vector<double> mRelVelNormal;
        std::vector<double> mRelVelTangential;

        bool ImpactIsNew(const int id_2)
        {
            return std::find(mId1.begin(), mId1.end(), id_2) != mId1.end();
        }
    };

class FaceParticleImpactDataOfATimeStep // It holds the historical information gathered for a single particle against flat walls
    {
        public:

        FaceParticleImpactDataOfATimeStep(): mNImpacts(0)/*, mId(0)*/{}
        FaceParticleImpactDataOfATimeStep(const int id) : mNImpacts(0)/*, mId(id)*/{}
        ~FaceParticleImpactDataOfATimeStep(){}

        void PushBackImpacts(const double time, const int id2, const double normal_vel, const double tang_vel)
        {
            ++mNImpacts;
            mTimes.push_back(time);
            mId2.push_back(id2);
            mRelVelNormal.push_back(normal_vel);
            mRelVelTangential.push_back(tang_vel);
        }

        void FillUpPythonLists(std::list<double>& times,
                               std::list<int>& neighbour_ids,
                               std::list<double>& normal_relative_vel,
                               std::list<double>& tangential_relative_vel)
        {
            times.clear();
            neighbour_ids.clear();
            normal_relative_vel.clear();
            tangential_relative_vel.clear();

            for (int i = 0; i < mNImpacts; ++i){
                times.push_back(mTimes[i]);
                neighbour_ids.push_back(mId2[i]);
                normal_relative_vel.push_back(mRelVelNormal[i]);
                tangential_relative_vel.push_back(mRelVelTangential[i]);
            }
        }

        void GetMaxCollidingSpeedFromDatabase(double& db_normal_impact_velocity, double& db_tangential_impact_velocity){
            if(mRelVelNormal.size()){
                db_normal_impact_velocity = std::abs(*(std::max_element(mRelVelNormal.begin(), mRelVelNormal.end())));
                db_tangential_impact_velocity = std::abs(*(std::max_element(mRelVelTangential.begin(), mRelVelTangential.end())));
            }
            else {
                db_normal_impact_velocity = 0.0;
                db_tangential_impact_velocity = 0.0;
            }

        }

    private:

        int mNImpacts;
        //int mId;
        std::vector<double> mTimes;
        std::vector<int> mId2;
        std::vector<double> mRelVelNormal;
        std::vector<double> mRelVelTangential;
};

// 

void GetParticleData(int id,
                     std::list<double> times,
                     std::list<int> neighbour_ids,
                     std::list<double> normal_relative_vel,
                     std::list<double> tangential_relative_vel);

void GetAllParticlesData(ModelPart& analytic_model_part,
                         std::list<double> times,
                         std::list<int> neighbour_ids,
                         std::list<double> normal_relative_vel,
                         std::list<double> tangential_relative_vel);

void GetTimeStepsData(std::list<int> ids,
                      std::list<int> neighbour_ids,
                      std::list<double> normal_relative_vel,
                      std::list<double> tangential_relative_vel);

virtual void MakeMeasurements(ModelPart &analytic_model_part);

virtual void SetNodalMaxImpactVelocities(ModelPart &analytic_model_part);
virtual void SetNodalMaxFaceImpactVelocities(ModelPart &analytic_model_part);
virtual void SetNodalMaxLinearImpulse(ModelPart &analytic_model_part);

virtual InterParticleImpactDataOfATimeStep& GetParticleDataBase(int id);
virtual FaceParticleImpactDataOfATimeStep& GetParticleFaceDataBase(int id);

/// Turn back information as a string
virtual std::string Info() const;

/// Print information about this object
virtual void PrintInfo(std::ostream& rOStream) const;

/// Print object's data
virtual void PrintData(std::ostream& rOStream) const;


private:

std::set<int> mSetOfIds;

std::vector<InterParticleImpactDataOfAParticle> mInterParticleImpactDataOfAllParticles;
std::map<int, InterParticleImpactDataOfATimeStep> mInterParticleImpactDataOfAllTimeSteps;

std::vector<FaceParticleImpactDataOfAParticle> mFaceParticleImpactDataOfAllParticles;
std::map<int, FaceParticleImpactDataOfATimeStep> mFaceParticleImpactDataOfAllTimeSteps;

// std::vector<EdgeParticleImpactDataOfAParticle> mEdgeParticleImpactDataOfAllParticles; inactive
// std::map<int, EdgeParticleImpactDataOfATimeStep> mEdgeParticleImpactDataOfAllTimeSteps;

// std::vector<VertexParticleImpactDataOfAParticle> mVertexParticleImpactDataOfAllParticles; inactive
// std::map<int, VertexParticleImpactDataOfATimeStep> mVertexParticleImpactDataOfAllTimeSteps;

/// Assignment operator
AnalyticParticleWatcher & operator=(AnalyticParticleWatcher const& rOther);

}; // Class AnalyticParticleWatcher

} // namespace Kratos.

#endif // ANALYTIC_PARTICLE_WATCHER_H
