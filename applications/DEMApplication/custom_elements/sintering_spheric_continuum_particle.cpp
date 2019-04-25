//Authors: S. Nosewicz (IPPT PAN, Warsaw, Poland) and M.A. Celigueta (CIMNE, maceli@cimne.upc.edu)
// Date: April 2015


// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "thermal_spheric_particle.h"
#include "sintering_spheric_continuum_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "utilities/openmp_utils.h"
#include "DEM_application_variables.h"


//......................
//KRATOS_INFO("DEM") <<print...<<std::endl;

namespace Kratos
{
	void SinteringSphericContinuumParticle::Initialize(const ProcessInfo& r_process_info) {
            KRATOS_TRY
            ThermalSphericParticle<SphericContinuumParticle>::Initialize(r_process_info);
            mSinteringDisplacement = 0;
            mSinteringDrivingForce = 0;
            KRATOS_CATCH("")
	}

        void SinteringSphericContinuumParticle::InitializeSolutionStep(ProcessInfo& r_process_info) {
            KRATOS_TRY
			ThermalSphericParticle<SphericContinuumParticle>::InitializeSolutionStep(r_process_info);
            const double temperature = GetTemperature();
            const double sintering_start_temp = GetProperties()[SINTERING_START_TEMPERATURE];

            if (temperature > sintering_start_temp) { this->Set(DEMFlags::IS_SINTERING, true);  }
            else                                    { this->Set(DEMFlags::IS_SINTERING, false); }
            KRATOS_CATCH("")
        }

	void SinteringSphericContinuumParticle::UpdateContinuumNeighboursVector(ProcessInfo& r_process_info) {
            KRATOS_TRY
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
                if(mNeighbourElements[i] == NULL) {
                    discont_ini_ids.push_back(mIniNeighbourIds[i]);
                    discont_ini_deltas.push_back(mIniNeighbourDelta[i]);
                    discont_ini_neighbour_elems.push_back(mNeighbourElements[i]);
                    continue;
                }
                array_1d<double, 3> other_to_me_vect;
                noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - mNeighbourElements[i]->GetGeometry()[0].Coordinates();
                double distance = DEM_MODULUS_3(other_to_me_vect);
                double radius_sum = GetRadius() + mNeighbourElements[i]->GetRadius();
                double initial_dist = radius_sum - mIniNeighbourDelta[i];
                double indentation = initial_dist - distance;

                if(indentation>0.0 && mNeighbourElements[i]->Is(DEMFlags::IS_SINTERING)){
                    cont_ini_neighbour_elems.push_back(mNeighbourElements[i]);
                    cont_ini_ids.push_back(mIniNeighbourIds[i]);
                    cont_ini_deltas.push_back(mIniNeighbourDelta[i]);
                    cont_ini_failure_ids.push_back(0);
                } else {
                    discont_ini_ids.push_back(mIniNeighbourIds[i]);
                    discont_ini_deltas.push_back(mIniNeighbourDelta[i]);
                    discont_ini_neighbour_elems.push_back(mNeighbourElements[i]);
                }
            }

            //The rest (discontinuum not initial)
            for (unsigned int i = mInitialNeighborsSize; i < mNeighbourElements.size(); i++) {
                array_1d<double, 3> other_to_me_vect;
                noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - mNeighbourElements[i]->GetGeometry()[0].Coordinates();
                double distance = DEM_MODULUS_3(other_to_me_vect);
                double radius_sum = GetRadius() + mNeighbourElements[i]->GetRadius();
                double indentation = radius_sum - distance;

                if(indentation>0.0 && mNeighbourElements[i]->Is(DEMFlags::IS_SINTERING)){
                    cont_ini_neighbour_elems.push_back(mNeighbourElements[i]);
                    cont_ini_ids.push_back(mNeighbourElements[i]->Id());
                    cont_ini_deltas.push_back(0.0);
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
            mActualNeighbourSinteringDisplacement.resize(mContinuumInitialNeighborsSize);
            mOldNeighbourSinteringDisplacement.resize(mContinuumInitialNeighborsSize);

            //Continuum
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

                mBondElements.push_back(NULL);
                //std::string ElementName;
                //ElementName = std::string("ParticleContactElement");
                //const Element& rReferenceElement = KratosComponents<Element>::Get(ElementName);

                mActualNeighbourSinteringDisplacement[i] = 0.0;
                mOldNeighbourSinteringDisplacement[i] = 0.0;
            }

            for (unsigned int i = mContinuumInitialNeighborsSize; i < mInitialNeighborsSize; i++) {
                mNeighbourElements[i] = discont_ini_neighbour_elems[i-mContinuumInitialNeighborsSize];
                mIniNeighbourIds[i]   = discont_ini_ids[i-mContinuumInitialNeighborsSize];
                mIniNeighbourDelta[i] = discont_ini_deltas[i-mContinuumInitialNeighborsSize];
		if(discont_ini_neighbour_elems.size() != mInitialNeighborsSize - mContinuumInitialNeighborsSize) {
		  KRATOS_WATCH("ERROR 1 IN SINTERING PARTICLE WHEN CONVERTING DISCONTINUUM INTO CONTINUUM")
		}
            }

            for (unsigned int i = mInitialNeighborsSize; i < mNeighbourElements.size(); i++) {
                mNeighbourElements[i] = discont_neighbour_elems[i-mInitialNeighborsSize];
		if(discont_neighbour_elems.size() != mNeighbourElements.size() - mInitialNeighborsSize) {
		  KRATOS_WATCH("ERROR 1 IN SINTERING PARTICLE WHEN CONVERTING DISCONTINUUM INTO CONTINUUM")
		}
            }

            KRATOS_CATCH("")
	};


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
	}//SetInitialSinteringSphereContacts

	void SinteringSphericContinuumParticle::InitializeForceComputation(ProcessInfo& r_process_info)
	{
                if (this->Is(DEMFlags::IS_SINTERING))
                {
                        UpdateContinuumNeighboursVector(r_process_info);
                        mOldNeighbourSinteringDisplacement = mActualNeighbourSinteringDisplacement;
                        mActualNeighbourSinteringDisplacement.clear();
                }
	}

        void SinteringSphericContinuumParticle::ComputeOtherBallToBallForces(array_1d<double, 3>& other_ball_to_ball_forces) {
            other_ball_to_ball_forces[2] = -this->mSinteringDrivingForce;
        }

        /*double SinteringSphericContinuumParticle::GetInitialDelta(int index) {
            return 0.0;
        }*/

        void SinteringSphericContinuumParticle::ComputeContactArea(const double rmin, double indentation, double& calculation_area) { //TODO: add "thermal" word to this function
                double actual_neck_radius;
                if (this->Is(DEMFlags::IS_SINTERING)) {
                        indentation = -indentation;
                        double geo_a = rmin;
                        double geo_c = rmin;
                        double geo_b = rmin + rmin + indentation; //// ZMIANY
                        double geo_aproj = (geo_a*geo_a + geo_b*geo_b - geo_c*geo_c) / (2 * geo_b);
                        actual_neck_radius = 1.43 * std::sqrt(geo_a*geo_a - geo_aproj*geo_aproj);
                }
                else {
                        actual_neck_radius = std::sqrt(rmin*indentation);
                }

                calculation_area = Globals::Pi * actual_neck_radius * actual_neck_radius;
        }
}  // namespace Kratos.
