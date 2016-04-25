
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

	void SinteringSphericContinuumParticle::ComputeBallToBallContactForce(array_1d<double, 3>& rElasticForce,
                                                                            array_1d<double, 3 > & rContactForce,
                                                                            array_1d<double, 3>& rInitialRotaMoment,
                                                                            ProcessInfo& r_process_info,
                                                                            const double dt,
                                                                            const bool multi_stage_RHS) {
		KRATOS_TRY

                const int time_steps = r_process_info[TIME_STEPS];
		const int& search_control = r_process_info[SEARCH_CONTROL];
		vector<int>& search_control_vector = r_process_info[SEARCH_CONTROL_VECTOR];

		const array_1d<double, 3>& vel = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
		const array_1d<double, 3>& delta_displ = this->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
		const array_1d<double, 3>& ang_vel = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);

		for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {

			if (mNeighbourElements[i] == NULL) continue;
			if (this->Is(NEW_ENTITY) && mNeighbourElements[i]->Is(NEW_ENTITY)) continue;

			ThermalSphericParticle<SphericContinuumParticle>* neighbour_iterator = dynamic_cast<ThermalSphericParticle<SphericContinuumParticle>*>(mNeighbourElements[i]);

			unsigned int neighbour_iterator_id = neighbour_iterator->Id();

			array_1d<double, 3> other_to_me_vect;
			noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();

			const double& other_radius = neighbour_iterator->GetRadius();

			double distance = DEM_MODULUS_3(other_to_me_vect);
			double radius_sum = GetRadius() + other_radius;
                        
			double initial_delta = GetInitialDelta(i);
            
			double initial_dist = radius_sum - initial_delta;
			double indentation = initial_dist - distance;
			double myYoung = GetYoung();
			double myPoisson = GetPoisson();

			double kn_el = 0.0;
			double kt_el = 0.0;
			double DeltDisp[3] = { 0.0 };
			double RelVel[3] = { 0.0 };
			double LocalCoordSystem[3][3] = { { 0.0 },{ 0.0 },{ 0.0 } };
			double OldLocalCoordSystem[3][3] = { { 0.0 },{ 0.0 },{ 0.0 } };
			bool sliding = false;

			double contact_tau = 0.0;
			double contact_sigma = 0.0;
			double failure_criterion_state = 0.0;
			double acumulated_damage = 0.0;

			// Getting neighbour properties
			double other_young = neighbour_iterator->GetYoung();
			double other_poisson = neighbour_iterator->GetPoisson();
			double equiv_poisson;
			if ((myPoisson + other_poisson) != 0.0) { equiv_poisson = 2.0 * myPoisson * other_poisson / (myPoisson + other_poisson); }
			else { equiv_poisson = 0.0; }

			double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);
			double calculation_area = 0.0;

	//		if (i < mContinuumInitialNeighborsSize && this->IsNot(DEMFlags::IS_SINTERING)) {
	//			mContinuumConstitutiveLawArray[i]->GetContactArea(GetRadius(), other_radius, mContIniNeighArea, i, calculation_area);
	//			mContinuumConstitutiveLawArray[i]->CalculateElasticConstants(kn_el, kt_el, initial_dist, equiv_young, equiv_poisson, calculation_area);
	//		}

			EvaluateDeltaDisplacement(DeltDisp, RelVel, LocalCoordSystem, OldLocalCoordSystem, other_to_me_vect, vel, delta_displ, neighbour_iterator, distance);

			if (this->Is(DEMFlags::HAS_ROTATION)) {
				DisplacementDueToRotationMatrix(DeltDisp, RelVel, OldLocalCoordSystem, other_radius, dt, ang_vel, neighbour_iterator);
			}

			double LocalDeltDisp[3] = { 0.0 };
			double LocalElasticContactForce[3] = { 0.0 }; // 0: first tangential, // 1: second tangential, // 2: normal force
			double GlobalElasticContactForce[3] = { 0.0 };
			double OldLocalElasticContactForce[3] = { 0.0 };

			GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltDisp, LocalDeltDisp);

			RotateOldContactForces(OldLocalCoordSystem, LocalCoordSystem, mNeighbourElasticContactForces[i]);

			// Here we recover the old local forces projected in the new coordinates in the way they were in the old ones; Now they will be increased if necessary
			GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, mNeighbourElasticContactForces[i], OldLocalElasticContactForce);

			GlobalElasticContactForce[0] = mNeighbourElasticContactForces[i][0];
			GlobalElasticContactForce[1] = mNeighbourElasticContactForces[i][1];
			GlobalElasticContactForce[2] = mNeighbourElasticContactForces[i][2];

			GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce);
			//we recover this way the old local forces projected in the new coordinates in the way they were in the old ones; Now they will be increased if its the necessary

			double ViscoDampingLocalContactForce[3] = { 0.0 };
			double equiv_visco_damp_coeff_normal;
			double equiv_visco_damp_coeff_tangential;
            double ElasticLocalRotationalMoment[3] = {0.0};
            double ViscoLocalRotationalMoment[3] = {0.0};

			double LocalRelVel[3] = { 0.0 };
			GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, RelVel, LocalRelVel);
			double rel_vel = LocalRelVel[2];
			double thermalRelVel = 0;
			//UpdateTemperatureDependentRadii(r_process_info, this, neighbour_iterator);
			UpdateNormalRelativeVelocityDueToThermalExpansion(r_process_info, thermalRelVel, this, neighbour_iterator);
			rel_vel = rel_vel - thermalRelVel;

			double penetration; 

			if (i < mContinuumInitialNeighborsSize) {
				int failure_id = mIniNeighbourFailureId[i];
				penetration = radius_sum - distance;
				sintering_displ = mOldNeighbourSinteringDisplacement[i];

				if (this->Is(DEMFlags::IS_SINTERING))
				{
					mContinuumConstitutiveLawArray[i]->CalculateForcesOfSintering(r_process_info, OldLocalElasticContactForce, LocalElasticContactForce, rel_vel, penetration, 
					sintering_displ, sinter_driv_force, this, neighbour_iterator, ViscoDampingLocalContactForce);
					mActualNeighbourSinteringDisplacement.push_back(sintering_displ);  // adding the sintering displacements to vector (only for sintering period (continuum CL)
				}
				else 
				{ 
					mContinuumConstitutiveLawArray[i]->CalculateForces(r_process_info, LocalElasticContactForce, LocalDeltDisp, kn_el, kt_el, contact_sigma, contact_tau, failure_criterion_state, sintering_displ, penetration,
																	   calculation_area, acumulated_damage, this, neighbour_iterator, i, r_process_info[TIME_STEPS], sliding, search_control, search_control_vector);
					mContinuumConstitutiveLawArray[i]->CalculateViscoDampingCoeff(equiv_visco_damp_coeff_normal, equiv_visco_damp_coeff_tangential, this, neighbour_iterator, kn_el, kt_el);
					mContinuumConstitutiveLawArray[i]->CalculateViscoDamping(LocalRelVel, ViscoDampingLocalContactForce, penetration, equiv_visco_damp_coeff_normal, equiv_visco_damp_coeff_tangential, sliding, failure_id);
				}
				
			}
			else if (indentation > 0.0) {
				double cohesive_force = 0.0;
				const double previous_indentation = indentation + LocalDeltDisp[2];
				mDiscontinuumConstitutiveLaw->CalculateForces(r_process_info, OldLocalElasticContactForce, LocalElasticContactForce, LocalDeltDisp, LocalRelVel, indentation, previous_indentation, ViscoDampingLocalContactForce, cohesive_force, this, neighbour_iterator, sliding);
			}

			// Transforming to global forces and adding up
			double LocalContactForce[3] = { 0.0 };
			double GlobalContactForce[3] = { 0.0 };

			if (this->Is(DEMFlags::HAS_STRESS_TENSOR) && (i < mContinuumInitialNeighborsSize)) { // We leave apart the discontinuum neighbors (the same for the walls). The neighbor would not be able to do the same if we activate it. 
				mContinuumConstitutiveLawArray[i]->AddPoissonContribution(equiv_poisson, LocalCoordSystem, LocalElasticContactForce[2], calculation_area, mSymmStressTensor, this, neighbour_iterator);
			}

			AddUpForcesAndProject(OldLocalCoordSystem, LocalCoordSystem, LocalContactForce, LocalElasticContactForce, GlobalContactForce,
				GlobalElasticContactForce, ViscoDampingLocalContactForce, 0.0, sinter_driv_force, rElasticForce, rContactForce, i);


			if (this->IsNot(DEMFlags::IS_SINTERING))
			{
				if (this->Is(DEMFlags::HAS_ROTATION)) {
					if (this->Is(DEMFlags::HAS_ROLLING_FRICTION) && !multi_stage_RHS) {
						const double coeff_acc = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) / dt;
						noalias(rInitialRotaMoment) = coeff_acc * ang_vel; // the moment needed to stop the spin in a single time step
					}
					ComputeMoments(LocalElasticContactForce[2], mNeighbourTotalContactForces[i], rInitialRotaMoment, LocalCoordSystem[2], neighbour_iterator, indentation);
					if (i < mContinuumInitialNeighborsSize) {
						if (mIniNeighbourFailureId[i] == 0) {
							mContinuumConstitutiveLawArray[i]->ComputeParticleRotationalMoments(this, neighbour_iterator, equiv_young, distance, calculation_area,
                                                                                            LocalCoordSystem, ElasticLocalRotationalMoment, ViscoLocalRotationalMoment);
						}
					}
				}
                        }
            
                        AddUpMomentsAndProject(LocalCoordSystem, ElasticLocalRotationalMoment, ViscoLocalRotationalMoment);

			if (r_process_info[CONTACT_MESH_OPTION] == 1 && (i < mContinuumInitialNeighborsSize) && this->Id() < neighbour_iterator_id) {
				CalculateOnContactElements(i, LocalElasticContactForce, contact_sigma, contact_tau, failure_criterion_state, acumulated_damage, time_steps);
			}

			if (this->Is(DEMFlags::HAS_STRESS_TENSOR) && (i < mContinuumInitialNeighborsSize)) {
				AddNeighbourContributionToStressTensor(GlobalElasticContactForce, LocalCoordSystem[2], distance, radius_sum);
			}

			AddContributionToRepresentativeVolume(distance, radius_sum, calculation_area);

		} // for each neighbor
		
		KRATOS_CATCH("")
	} //  ComputeBallToBallContactForce   
	 
}  // namespace Kratos.
