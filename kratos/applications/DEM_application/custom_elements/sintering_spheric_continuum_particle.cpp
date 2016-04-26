
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
#include "DEM_application.h"
#include "utilities/openmp_utils.h"
#include "DEM_application_variables.h"


//......................
//std::cout<<print...<<std::endl;

namespace Kratos
{
	void SinteringSphericContinuumParticle::Initialize(const ProcessInfo& r_process_info) {
	
		ThermalSphericParticle<SphericContinuumParticle>::Initialize(r_process_info);
		sintering_displ = 0;
		sinter_driv_force = 0;
	}

	void SinteringSphericContinuumParticle::UpdatingNeighboursVector(ProcessInfo& r_process_info)
	{
		//double h = 1;
		//double h1;double
		//h1 = GetRadius() + 1;
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
	}//SetInitialSphereContacts    /*

	void SinteringSphericContinuumParticle::InitializeForceComputation(ProcessInfo& r_process_info)
	{
		double temperature = GetTemperature();
		double sintering_start_temp = GetProperties()[SINTERING_START_TEMPERATURE];
		//UpdateTemperatureDependentRadius(r_process_info);
		

		if (temperature > sintering_start_temp)
		{
			if (this->Is(DEMFlags::IS_SINTERING))
			{
				UpdatingNeighboursVector(r_process_info);
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
	};

	void SinteringSphericContinuumParticle::AddUpForcesAndProject(double OldCoordSystem[3][3],
		double LocalCoordSystem[3][3],
		double LocalContactForce[3],
		double LocalElasticContactForce[3],
		double GlobalContactForce[3],
		double GlobalElasticContactForce[3],
		double ViscoDampingLocalContactForce[3],
		const double cohesive_force,
		double sinter_driv_force,
		array_1d<double, 3>& r_elastic_force,
		array_1d<double, 3>& r_contact_force,
		const unsigned int i_neighbour_count)
	{
		for (unsigned int index = 0; index < 3; index++) {
			LocalContactForce[index] = LocalElasticContactForce[index] + ViscoDampingLocalContactForce[index];
		}
		LocalContactForce[2] -= cohesive_force;  /// ????
		LocalContactForce[2] = LocalContactForce[2] - sinter_driv_force;

		GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalElasticContactForce, GlobalElasticContactForce);
		GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalContactForce, GlobalContactForce);

		// Saving contact forces (We need to, since tangential elastic force is history-dependent)
		DEM_COPY_SECOND_TO_FIRST_3(mNeighbourElasticContactForces[i_neighbour_count], GlobalElasticContactForce)
			DEM_COPY_SECOND_TO_FIRST_3(mNeighbourTotalContactForces[i_neighbour_count], GlobalContactForce)
			DEM_ADD_SECOND_TO_FIRST(mContactForce, GlobalContactForce)
			DEM_ADD_SECOND_TO_FIRST(r_elastic_force, GlobalElasticContactForce)
			DEM_ADD_SECOND_TO_FIRST(r_contact_force, GlobalContactForce)
	}


	 
}  // namespace Kratos.
