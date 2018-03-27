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


class InterParticleImpactDataOfAllParticlesSingleTimeStep  // It holds the historical information gathered in a single time step
{
    public:

    InterParticleImpactDataOfAllParticlesSingleTimeStep(const double time) : mNImpacts(0)/*, mTime(time)*/{}
    ~InterParticleImpactDataOfAllParticlesSingleTimeStep(){}

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

    void PushBackImpacts(InterParticleImpactDataOfAllParticlesSingleTimeStep& other_list_of_impacts)
    {
        for (int i=0; i<other_list_of_impacts.GetNumberOfImpacts(); i++) {
            PushBackImpacts(other_list_of_impacts.mId1[i], other_list_of_impacts.mId2[i], other_list_of_impacts.mRelVelNormal[i], other_list_of_impacts.mRelVelTangential[i]);
        }
    }

    void FillUpPythonLists(boost::python::list& ids,
                           boost::python::list& neighbour_ids,
                           boost::python::list& normal_relative_vel,
                           boost::python::list& tangential_relative_vel)
    {
        AnalyticParticleWatcher::ClearList(ids);
        AnalyticParticleWatcher::ClearList(neighbour_ids);
        AnalyticParticleWatcher::ClearList(normal_relative_vel);
        AnalyticParticleWatcher::ClearList(tangential_relative_vel);

        for (int i = 0; i < mNImpacts; ++i){            
            ids.append(mId1[i]);
            neighbour_ids.append(mId2[i]);
            normal_relative_vel.append(mRelVelNormal[i]);
            tangential_relative_vel.append(mRelVelTangential[i]);
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
            return std::find(mId1.begin(), mId1.end(), id_2) == mId1.end();
        }
    };

class InterParticleImpactDataOfAllTimeStepsSingleParticle // It holds the historical information gathered for a single particle
    {
        public:

        InterParticleImpactDataOfAllTimeStepsSingleParticle(): mNImpacts(0)/*, mId(0)*/{}
        InterParticleImpactDataOfAllTimeStepsSingleParticle(const int id) : mNImpacts(0)/*, mId(id)*/{}
        ~InterParticleImpactDataOfAllTimeStepsSingleParticle(){}

        void PushBackImpacts(const double time, const int id2, const double normal_vel, const double tang_vel, const double linear_impulse)
        {
            ++mNImpacts;
            mTimes.push_back(time);
            mId2.push_back(id2);
            mRelVelNormal.push_back(normal_vel);
            mRelVelTangential.push_back(tang_vel);
            mLinearImpulse.push_back(linear_impulse);            
        }

        void FillUpPythonLists(boost::python::list& times,
                               boost::python::list& neighbour_ids,
                               boost::python::list& normal_relative_vel,
                               boost::python::list& tangential_relative_vel)
        {
            AnalyticParticleWatcher::ClearList(times);
            AnalyticParticleWatcher::ClearList(neighbour_ids);
            AnalyticParticleWatcher::ClearList(normal_relative_vel);
            AnalyticParticleWatcher::ClearList(tangential_relative_vel);

            for (int i = 0; i < mNImpacts; ++i){                
                times.append(mTimes[i]);
                neighbour_ids.append(mId2[i]);
                normal_relative_vel.append(mRelVelNormal[i]);
                tangential_relative_vel.append(mRelVelTangential[i]);
            }
        }

        void GetMaxCollidingSpeedFromDatabase(double& db_normal_impact_velocity, double& db_tangential_impact_velocity){
            if(mRelVelNormal.size()){
                for(int i=0; i<(int)mRelVelNormal.size(); i++){
                    const double abs_normal_value = std::abs(mRelVelNormal[i]);
                    db_normal_impact_velocity = std::max(db_normal_impact_velocity, abs_normal_value);
                    const double abs_tg_value = std::abs(mRelVelTangential[i]);
                    db_tangential_impact_velocity = std::max(db_tangential_impact_velocity, abs_tg_value);
                }
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




// FaceParticleImpactDataOfAllParticlesSingleTimeStep

class FaceParticleImpactDataOfAllParticlesSingleTimeStep  // It holds the historical information gathered in a single time step against flat walls
{
    public:

    FaceParticleImpactDataOfAllParticlesSingleTimeStep(const double time) : mNImpacts(0)/*, mTime(time)*/{}
    ~FaceParticleImpactDataOfAllParticlesSingleTimeStep(){}

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

    void PushBackImpacts(FaceParticleImpactDataOfAllParticlesSingleTimeStep& other_list_of_impacts)
    {
        for (int i=0; i<other_list_of_impacts.GetNumberOfImpacts(); i++) {
            PushBackImpacts(other_list_of_impacts.mId1[i], other_list_of_impacts.mId2[i], other_list_of_impacts.mRelVelNormal[i], other_list_of_impacts.mRelVelTangential[i]);
        }
    }

    void FillUpPythonLists(boost::python::list& ids,
                           boost::python::list& neighbour_ids,
                           boost::python::list& normal_relative_vel,
                           boost::python::list& tangential_relative_vel)
    {
        AnalyticParticleWatcher::ClearList(ids);
        AnalyticParticleWatcher::ClearList(neighbour_ids);
        AnalyticParticleWatcher::ClearList(normal_relative_vel);
        AnalyticParticleWatcher::ClearList(tangential_relative_vel);

        for (int i = 0; i < mNImpacts; ++i){            
            ids.append(mId1[i]);
            neighbour_ids.append(mId2[i]);
            normal_relative_vel.append(mRelVelNormal[i]);
            tangential_relative_vel.append(mRelVelTangential[i]);
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
            return std::find(mId1.begin(), mId1.end(), id_2) == mId1.end();
        }
    };

class FaceParticleImpactDataOfAllTimeStepsSingleParticle // It holds the historical information gathered for a single particle against flat walls
    {
        public:

        FaceParticleImpactDataOfAllTimeStepsSingleParticle(): mNImpacts(0)/*, mId(0)*/{}
        FaceParticleImpactDataOfAllTimeStepsSingleParticle(const int id) : mNImpacts(0)/*, mId(id)*/{}
        ~FaceParticleImpactDataOfAllTimeStepsSingleParticle(){}

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
            AnalyticParticleWatcher::ClearList(times);
            AnalyticParticleWatcher::ClearList(neighbour_ids);
            AnalyticParticleWatcher::ClearList(normal_relative_vel);
            AnalyticParticleWatcher::ClearList(tangential_relative_vel);

            for (int i = 0; i < mNImpacts; ++i){                
                times.append(mTimes[i]);
                neighbour_ids.append(mId2[i]);
                normal_relative_vel.append(mRelVelNormal[i]);
                tangential_relative_vel.append(mRelVelTangential[i]);
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

virtual void MakeMeasurements(ModelPart &analytic_model_part);

virtual void SetNodalMaxImpactVelocities(ModelPart &analytic_model_part);
virtual void SetNodalMaxFaceImpactVelocities(ModelPart &analytic_model_part);
virtual void SetNodalMaxLinearImpulse(ModelPart &analytic_model_part);

virtual InterParticleImpactDataOfAllTimeStepsSingleParticle& GetParticleDataBase(int id, std::map<int, InterParticleImpactDataOfAllTimeStepsSingleParticle>& data_base);
virtual FaceParticleImpactDataOfAllTimeStepsSingleParticle& GetParticleFaceDataBase(int id, std::map<int, FaceParticleImpactDataOfAllTimeStepsSingleParticle>& data_base);

/// Turn back information as a string
virtual std::string Info() const;

/// Print information about this object
virtual void PrintInfo(std::ostream& rOStream) const;

/// Print object's data
virtual void PrintData(std::ostream& rOStream) const;


private:

std::set<int> mSetOfIds;

std::vector<InterParticleImpactDataOfAllParticlesSingleTimeStep> mInterParticleImpactDataOfAllParticles;
std::map<int, InterParticleImpactDataOfAllTimeStepsSingleParticle> mInterParticleImpactDataOfAllTimeSteps;

std::vector<FaceParticleImpactDataOfAllParticlesSingleTimeStep> mFaceParticleImpactDataOfAllParticles;
std::map<int, FaceParticleImpactDataOfAllTimeStepsSingleParticle> mFaceParticleImpactDataOfAllTimeSteps;

// std::vector<EdgeParticleImpactDataOfAParticle> mEdgeParticleImpactDataOfAllParticles; inactive
// std::map<int, EdgeParticleImpactDataOfATimeStep> mEdgeParticleImpactDataOfAllTimeSteps;

// std::vector<VertexParticleImpactDataOfAParticle> mVertexParticleImpactDataOfAllParticles; inactive
// std::map<int, VertexParticleImpactDataOfATimeStep> mVertexParticleImpactDataOfAllTimeSteps;

/// Assignment operator
AnalyticParticleWatcher & operator=(AnalyticParticleWatcher const& rOther);

}; // Class AnalyticParticleWatcher

} // namespace Kratos.

#endif // ANALYTIC_PARTICLE_WATCHER_H
