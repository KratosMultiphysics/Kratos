
//   Project Name:        Kratos
//   Last Modified by:    $Author: Ferran Arrufat $
//   Date:                $Date: 2016-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


// System includes
#include <string>
#include <iostream> 

// External includes


// Project includes
#include "thermal_spheric_particle.h"
#include "sintering_spheric_continuum_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
//#include "DEM_application.h"
#include "utilities/openmp_utils.h"
#include "DEM_application_variables.h"


//......................
//std::cout<<print...<<std::endl;

namespace Kratos
{
	void SinteringSphericContinuumParticle::Initialize(const ProcessInfo& r_process_info) {
	
		ThermalSphericParticle<SphericContinuumParticle>::Initialize(r_process_info);
		mSinteringDisplacement = 0;
		mSinteringDrivingForce = 0;
	}

	void SinteringSphericContinuumParticle::UpdateContinuumNeighboursVector(ProcessInfo& r_process_info)
	{
            if(mNeighbourElements.size() == mContinuumInitialNeighborsSize) return;
            
            const unsigned int initial_number_of_cont_neighbours = mContinuumInitialNeighborsSize;
            
            std::vector<SphericParticle*> cont_ini_neighbour_elems;
            std::vector<int> cont_ini_ids;
            std::vector<double> cont_ini_deltas;
            std::vector<int> cont_ini_failure_ids;  
            std::vector<SphericParticle*> discont_ini_neighbour_elems; 
            std::vector<int> discont_ini_ids;
            std::vector<double> discont_ini_deltas;
            std::vector<SphericParticle*> discont_neighbour_elems;
                        
            //Assuming that the neighbours are ordered (first initial continuum, then initial discontinuum, then other discontinuum)
            for (unsigned int i = 0; i < mContinuumInitialNeighborsSize; i++) {
                cont_ini_neighbour_elems.push_back(mNeighbourElements[i]);
                cont_ini_ids.push_back(mIniNeighbourIds[i]);
                cont_ini_deltas.push_back(mIniNeighbourDelta[i]);
                cont_ini_failure_ids.push_back(mIniNeighbourFailureId[i]);                
            }            

            //// Discontinuum and initial
            for (unsigned int i = mContinuumInitialNeighborsSize; i < mInitialNeighborsSize; i++) {
                array_1d<double, 3> other_to_me_vect;
                noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - mNeighbourElements[i]->GetGeometry()[0].Coordinates();
                double distance = DEM_MODULUS_3(other_to_me_vect);
                double radius_sum = GetRadius() + mNeighbourElements[i]->GetRadius();
                double initial_delta = GetInitialDelta(i);
                double initial_dist = radius_sum - initial_delta;
                double indentation = initial_dist - distance;
                
                if(indentation>0.0){
                    cont_ini_neighbour_elems.push_back(mNeighbourElements[i]);
                    cont_ini_ids.push_back(mIniNeighbourIds[i]);
                    cont_ini_deltas.push_back(mIniNeighbourDelta[i]);
                    cont_ini_failure_ids.push_back(0);                        
                } else {
                    discont_ini_ids.push_back(mNeighbourElements[i]->Id());
                    discont_ini_deltas.push_back(initial_delta);
                    discont_ini_neighbour_elems.push_back(mNeighbourElements[i]);                    
                }            
            }
            
            
            //The rest (discontinuum not initial)
            for (unsigned int i = mInitialNeighborsSize; i < mNeighbourElements.size(); i++) {
                array_1d<double, 3> other_to_me_vect;
                noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - mNeighbourElements[i]->GetGeometry()[0].Coordinates();
                double distance = DEM_MODULUS_3(other_to_me_vect);
                double radius_sum = GetRadius() + mNeighbourElements[i]->GetRadius();
                double initial_delta = GetInitialDelta(i);
                double initial_dist = radius_sum - initial_delta;
                double indentation = initial_dist - distance;
                
                if(indentation>0.0){
                    cont_ini_neighbour_elems.push_back(mNeighbourElements[i]);
                    cont_ini_ids.push_back(mIniNeighbourIds[i]);
                    cont_ini_deltas.push_back(mIniNeighbourDelta[i]);
                    cont_ini_failure_ids.push_back(0);                        
                } else {
                    discont_neighbour_elems.push_back(mNeighbourElements[i]);                    
                }            
            }

            //Adding everything up in single vectors mNeighbourElements, mIniNeighbourIds, mIniNeighbourDelta  
            mContinuumInitialNeighborsSize = cont_ini_neighbour_elems.size();
            mInitialNeighborsSize = mContinuumInitialNeighborsSize + discont_ini_neighbour_elems.size();
            
            mIniNeighbourIds.resize(mInitialNeighborsSize);
            mIniNeighbourDelta.resize(mInitialNeighborsSize);
            mIniNeighbourFailureId.resize(mContinuumInitialNeighborsSize);
            
            //Continuum initial
            for (unsigned int i = 0; i < mContinuumInitialNeighborsSize; i++) {
                mNeighbourElements[i] = cont_ini_neighbour_elems[i];
                mIniNeighbourIds[i]   = cont_ini_ids[i];
                mIniNeighbourDelta[i] = cont_ini_deltas[i]; 
                mIniNeighbourFailureId[i] = cont_ini_failure_ids[i];
            }

            mContinuumConstitutiveLawArray.resize(mContinuumInitialNeighborsSize);

            for (unsigned int i = initial_number_of_cont_neighbours; i < mContinuumInitialNeighborsSize; i++) {
                DEMContinuumConstitutiveLaw::Pointer NewContinuumConstitutiveLaw = GetProperties()[DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER]-> Clone();
                mContinuumConstitutiveLawArray[i] = NewContinuumConstitutiveLaw;
                mContinuumConstitutiveLawArray[i]->Initialize();
            }

            for (unsigned int i = mContinuumInitialNeighborsSize; i < mInitialNeighborsSize; i++) {
                mNeighbourElements[i] = discont_ini_neighbour_elems[i];
                mIniNeighbourIds[i]   = discont_ini_ids[i];
                mIniNeighbourDelta[i] = discont_ini_deltas[i];                
            }
            
            for (unsigned int i = mInitialNeighborsSize; i < mNeighbourElements.size(); i++) {
                mNeighbourElements[i] = discont_neighbour_elems[i];                
            }                            
            
	};
        
        
        void SphericContinuumParticle::SetInitialSphereContacts(ProcessInfo& r_process_info) {
         
        std::vector<SphericContinuumParticle*> ContinuumInitialNeighborsElements;
        std::vector<SphericContinuumParticle*> DiscontinuumInitialNeighborsElements;
        std::vector<int> DiscontinuumInitialNeighborsIds;
        std::vector<double> DiscontinuumInitialNeighborsDeltas;
        mIniNeighbourFailureId.clear(); // We will have to build this vector, we still don't know its size, it applies only to continuum particles
        size_t continuum_ini_size    = 0;
        size_t discontinuum_ini_size = 0;
        unsigned int neighbours_size = mNeighbourElements.size();
        mIniNeighbourIds.resize(neighbours_size);
        mIniNeighbourDelta.resize(neighbours_size);

        for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
            
            SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            array_1d<double, 3> other_to_me_vect;
            noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();

            double distance = DEM_MODULUS_3(other_to_me_vect);
            double radius_sum = GetRadius() + neighbour_iterator->GetRadius();
            double initial_delta = radius_sum - distance;
            int r_other_continuum_group = neighbour_iterator->mContinuumGroup; // finding out neighbor's Continuum Group Id
            
            if ((r_other_continuum_group  == this->mContinuumGroup) && (this->mContinuumGroup != 0)) {
                
                mIniNeighbourIds[continuum_ini_size]   = neighbour_iterator->Id();
                mIniNeighbourDelta[continuum_ini_size] = initial_delta;
                mIniNeighbourFailureId.push_back(0);
                ContinuumInitialNeighborsElements.push_back(neighbour_iterator);
                continuum_ini_size++;

            } else {
                
                DiscontinuumInitialNeighborsIds.push_back(neighbour_iterator->Id());
                DiscontinuumInitialNeighborsDeltas.push_back(initial_delta);
                DiscontinuumInitialNeighborsElements.push_back(neighbour_iterator);
                discontinuum_ini_size++;
            }            
        }

        mContinuumInitialNeighborsSize = continuum_ini_size;
        mInitialNeighborsSize = neighbours_size;
        
        for (unsigned int j = 0; j < continuum_ini_size; j++) {
            mNeighbourElements[j] = ContinuumInitialNeighborsElements[j];
        }
            
        for (unsigned int k = 0; k < discontinuum_ini_size; k++) {
            
            mIniNeighbourIds[continuum_ini_size + k]   = DiscontinuumInitialNeighborsIds[k];
            mIniNeighbourDelta[continuum_ini_size + k] = DiscontinuumInitialNeighborsDeltas[k];
            mNeighbourElements[continuum_ini_size + k] = DiscontinuumInitialNeighborsElements[k];
        }
    }//SetInitialSphereContacts 
        
        

	void SinteringSphericContinuumParticle::SetInitialSinteringSphereContacts(ProcessInfo& r_process_info)
	{
		std::vector<SphericContinuumParticle*> ContinuumInitialNeighborsElements;
		std::vector<SphericContinuumParticle*> DiscontinuumInitialNeighborsElements;
		std::vector<int> DiscontinuumInitialNeighborsIds;
		std::vector<double> DiscontinuumInitialNeighborsDeltas;
		mIniNeighbourFailureId.clear(); 
		size_t continuum_ini_size = 0;
		size_t discontinuum_ini_size = 0;
		unsigned int neighbours_size = mNeighbourElements.size();
		mIniNeighbourIds.resize(neighbours_size);
		mIniNeighbourDelta.resize(neighbours_size);
		mActualNeighbourSinteringDisplacement.reserve(neighbours_size);
		mOldNeighbourSinteringDisplacement.reserve(neighbours_size);

		for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {

			SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
			array_1d<double, 3> other_to_me_vect;
			noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();

			double distance = DEM_MODULUS_3(other_to_me_vect);
			double radius_sum = GetRadius() + neighbour_iterator->GetRadius();
			double initial_delta = radius_sum - distance;
			//int r_other_continuum_group = neighbour_iterator->mContinuumGroup; // finding out neighbor's Continuum Group Id

			if (initial_delta > 0) {

				mIniNeighbourIds[continuum_ini_size] = neighbour_iterator->Id();
				mIniNeighbourDelta[continuum_ini_size] = initial_delta;
				mIniNeighbourFailureId.push_back(0);
				ContinuumInitialNeighborsElements.push_back(neighbour_iterator);
				continuum_ini_size++;

			}
			else {

				DiscontinuumInitialNeighborsIds.push_back(neighbour_iterator->Id());
				DiscontinuumInitialNeighborsDeltas.push_back(initial_delta);
				DiscontinuumInitialNeighborsElements.push_back(neighbour_iterator);
				discontinuum_ini_size++;
			}
		}

		mContinuumInitialNeighborsSize = continuum_ini_size;
		mInitialNeighborsSize = neighbours_size;

		for (unsigned int j = 0; j < continuum_ini_size; j++) {
			mNeighbourElements[j] = ContinuumInitialNeighborsElements[j];
		}

		for (unsigned int k = 0; k < discontinuum_ini_size; k++) {

			mIniNeighbourIds[continuum_ini_size + k] = DiscontinuumInitialNeighborsIds[k];
			mIniNeighbourDelta[continuum_ini_size + k] = DiscontinuumInitialNeighborsDeltas[k];
			mNeighbourElements[continuum_ini_size + k] = DiscontinuumInitialNeighborsElements[k];
		}

	CreateContinuumConstitutiveLaws(); //// Reorder????
	}//SetInitialSphereContacts    /*

	void SinteringSphericContinuumParticle::InitializeForceComputation(ProcessInfo& r_process_info)
	{
		double temperature = GetTemperature();
		double sintering_start_temp = GetProperties()[SINTERING_START_TEMPERATURE];
		

		if (temperature > sintering_start_temp)
		{
			if (this->Is(DEMFlags::IS_SINTERING))
			{
				UpdateContinuumNeighboursVector(r_process_info);
				mOldNeighbourSinteringDisplacement = mActualNeighbourSinteringDisplacement;
				mActualNeighbourSinteringDisplacement.clear();
			}
			else // first step of sintering procedure (T = TsintStart) during heating
			{
				SetInitialSinteringSphereContacts(r_process_info);
				this->Set(DEMFlags::IS_SINTERING, true);
			} 
		}
		else
		{
			if (this->Is(DEMFlags::IS_SINTERING)) // last step of sintering procedure (T = TsintStart) during cooling - transition from sintering model to KDEM
			{
				// Maybe we can put here SetInitialSphereContacts method to establish "the neighbours data vectors" in the typical way
				this->Set(DEMFlags::IS_SINTERING, false);
				//SetInitialSphereContacts;
			}
			else
			{
				return;// Discontiuum CL (before sintering) or Continuum CL (after sintering)
			}
		}
	}
        
        void SinteringSphericContinuumParticle::ComputeOtherBallToBallForces(array_1d<double, 3>& other_ball_to_ball_forces){
            other_ball_to_ball_forces[2] = -this->mSinteringDrivingForce;
        }

	
        double SinteringSphericContinuumParticle::GetInitialDelta(int index) {
            return 0.0;                 
        }

	 
}  // namespace Kratos.
