//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

// System includes

// External includes

// Project includes
#include "sintering_continuum.h"

namespace Kratos
{
	DEMContinuumConstitutiveLaw::Pointer SinteringContinuum::Clone() const {
		DEMContinuumConstitutiveLaw::Pointer p_clone(new SinteringContinuum(*this));
		return p_clone;
	}

	////////////////////// CALCULATE SINTERING NORMAL FORCE ///////////////////////////////

	void SinteringContinuum::GetContactArea(const double  radius,
		                                      const double  other_radius,
		                                      const Vector& vector_of_initial_areas,
		                                      const int     neighbour_position,
		                                      double&       calculation_area) {}

	void SinteringContinuum::CalculateElasticConstants(double& kn_el,
		                                                 double& kt_el,
		                                                 double  initial_dist,
		                                                 double  equiv_young,
		                                                 double  equiv_poisson,
		                                                 double  calculation_area,
		                                                 SphericContinuumParticle* element1,
		                                                 SphericContinuumParticle* element2,
		                                                 double indentation) {}

	void SinteringContinuum::CalculateSinteringForces(const ProcessInfo& r_process_info,
		                                                const double OldLocalElasticContactForce[3],
		                                                double       LocalElasticContactForce[3],
		                                                const double rel_vel,
		                                                double       indentation,
		                                                double&      sintering_displ,
		                                                double&      sinter_driv_force,
		                                                SphericContinuumParticle* element1,
		                                                SphericContinuumParticle* element2,
		                                                double ViscoDampingLocalContactForce[3])
	{
		double sintering_displ_old           = 0;
		double damping_force                 = 0;
		double equiv_visco_damp_coeff_normal = 0;
		double minimal_radius                = 0;
		double kn                            = 0;
		double delta_time                    = r_process_info[DELTA_TIME];
		indentation                          = -indentation;

		InitializeContact(element1, element2, indentation, minimal_radius, kn, sintering_displ);

		double geo_a              = minimal_radius;
		double geo_c              = minimal_radius;
		double geo_b              = minimal_radius + minimal_radius + indentation; // ZMIANY
		double geo_aproj          = (geo_a * geo_a + geo_b * geo_b - geo_c * geo_c) / (2 * geo_b);
		double actual_neck_radius = 1.43 * std::sqrt(geo_a * geo_a - geo_aproj * geo_aproj);
		double dihedral_angle     = element1->GetProperties()[DIHEDRAL_ANGLE];
		dihedral_angle            = dihedral_angle * (Globals::Pi / 180); // ZMIANY
		//double final_neck_radius = minimal_radius * std::sin((dihedral_angle * (Globals::Pi / 180)) / 2);
		double final_neck_radius  = minimal_radius * std::sin(dihedral_angle / 2); // ZMIANY

		CalculateDampingCoeff(equiv_visco_damp_coeff_normal, kn, element1, element2);

		sintering_displ_old = sintering_displ;
		double elastic_contact_force;
		double elastic_contact_force_old = -OldLocalElasticContactForce[2];

		if (actual_neck_radius < final_neck_radius) {
			double k_const = 1.3806504e-23;  // Boltzmann constant
			double R_const = 8.3144621;      // Gas constant

			double atomic_volume  = element1->GetProperties()[ATOMIC_VOLUME];
			double surface_energy = element1->GetProperties()[SURFACE_ENERGY];

			//double relaxation_time    = element1->GetProperties()[RELAXATION_TIME];
			//double large_visco_coeff  = element1->GetProperties()[LARGE_VISCOSITY_COEFFICIENT];
			//double thermal_alpha      = element1->GetProperties()[THERMAL_EXPANSION_COEFFICIENT];
			double pre_Dgb              = element1->GetProperties()[PRE_EXP_DIFFUSION_COEFFICIENT];
			double gb_width             = element1->GetProperties()[GB_WIDTH];
			double enth_activ           = element1->GetProperties()[ENTHAPLY_ACTIVATION];
			double temperature_element1 = element1->GetGeometry()[0].GetSolutionStepValue(TEMPERATURE);
			double temperature_element2 = element2->GetGeometry()[0].GetSolutionStepValue(TEMPERATURE);
			double temperature          = (temperature_element1 + temperature_element2) / 2;
			double D_gb                 = pre_Dgb * std::exp(-enth_activ / (R_const * temperature));
			double D_eff                = (D_gb * gb_width * atomic_volume) / (k_const * temperature);
			double visco_coeff          = (Globals::Pi * actual_neck_radius * actual_neck_radius * actual_neck_radius * actual_neck_radius) / (8 * D_eff);
			sinter_driv_force           = Globals::Pi * surface_energy * (4 * minimal_radius * (1 - std::cos(dihedral_angle / 2)) + actual_neck_radius * std::sin(dihedral_angle / 2));
			elastic_contact_force       = elastic_contact_force_old * (1 / kn - delta_time * 0.5 / visco_coeff) + rel_vel * delta_time;
			elastic_contact_force       = elastic_contact_force / (1 / kn + delta_time * 0.5 / visco_coeff);
			double sinter_vel           = (elastic_contact_force + elastic_contact_force_old) / (2 * visco_coeff);
			sintering_displ             = sintering_displ_old + sinter_vel * delta_time;
			damping_force               = equiv_visco_damp_coeff_normal * (rel_vel - sinter_vel);
		}
		else { // equilibrium state achieved
			sinter_driv_force     = 0;
			damping_force         = equiv_visco_damp_coeff_normal * rel_vel;
			elastic_contact_force = kn * (indentation - sintering_displ) * 0.666666666666666;
		}

		//sintering_displ                = - sintering_displ;
		LocalElasticContactForce[2]      = -elastic_contact_force;
		ViscoDampingLocalContactForce[2] = -damping_force;
	}

	////////////////////// CALCULATE SINTERING NORMAL FORCE ///////////////////////////////

	////////////////////// CALCULATE FORCE /////////////////////////////////////

	void SinteringContinuum::CalculateForcesOfSintering(const ProcessInfo& r_process_info,
		                                                  const double OldLocalElasticContactForce[3],
		                                                  double       LocalElasticContactForce[3],
		                                                  const double rel_vel,
		                                                  const double indentation,
		                                                  double&      sintering_displ,
		                                                  double&      sinter_driv_force,
		                                                  SphericContinuumParticle* element1,
		                                                  SphericContinuumParticle* element2,
		                                                  double ViscoDampingLocalContactForce[3])
	{
		KRATOS_TRY

		//SphericParticle* element1;
		//SphericParticle* element2;

		if (indentation >= 0.0) { // COMPRESSION
			CalculateSinteringForces(r_process_info,
			                        OldLocalElasticContactForce,
			                        LocalElasticContactForce,
			                        rel_vel,
			                        indentation,
			                        sintering_displ,
			                        sinter_driv_force,
			                        element1,
			                        element2,
			                        ViscoDampingLocalContactForce);
			}
			else {
				//std::system("pause");
			}

		KRATOS_CATCH("")
	}

	////////////////////// CALCULATE FORCE /////////////////////////////////////

	void SinteringContinuum::CalculateDampingCoeff(double& equiv_visco_damp_coeff_normal,
		                                             double  kn,
		                                             SphericContinuumParticle* const element1,
		                                             SphericContinuumParticle* const element2) {

		const double my_mass          = element1->GetMass();
		const double other_mass       = element2->GetMass();
		const double equiv_mass       = 1.0 / (1.0 / my_mass + 1.0 / other_mass);
		const double& equiv_gamma     = (*mpProperties)[DAMPING_GAMMA];
		equiv_visco_damp_coeff_normal = 2.0 * equiv_gamma * sqrt(equiv_mass * kn);
	}

	// MA
	void SinteringContinuum::CalculateForces(const ProcessInfo& r_process_info,
		                                       double       OldLocalElasticContactForce[3],
		                                       double       LocalElasticContactForce[3],
		                                       double       LocalElasticExtraContactForce[3],
		                                       double       LocalCoordSystem[3][3],
		                                       double       LocalDeltDisp[3],
		                                       const double kn_el,
		                                       const double kt_el,
		                                       double&      contact_sigma,
		                                       double&      contact_tau,
		                                       double&      failure_criterion_state,
		                                       double       equiv_young,
		                                       double       equiv_shear,
		                                       double       indentation,
		                                       double       calculation_area,
		                                       double&      acumulated_damage,
		                                       SphericContinuumParticle* element1,
		                                       SphericContinuumParticle* element2,
		                                       int     i_neighbour_count,
		                                       int     time_steps,
		                                       bool&   sliding,
		                                       double& equiv_visco_damp_coeff_normal,
		                                       double& equiv_visco_damp_coeff_tangential,
		                                       double  LocalRelVel[3],
		                                       double  ViscoDampingLocalContactForce[3]) {

		KRATOS_TRY

		SinteringSphericContinuumParticle* p_sintering_element1 = dynamic_cast<SinteringSphericContinuumParticle*>(element1);
		p_sintering_element1->mSinteringDisplacement = p_sintering_element1->mOldNeighbourSinteringDisplacement[i_neighbour_count];

		if (element1->Is(DEMThermalFlags::IS_SINTERING) && element2->Is(DEMThermalFlags::IS_SINTERING)) {
			CalculateForcesOfSintering(r_process_info, OldLocalElasticContactForce, LocalElasticContactForce, LocalRelVel[2], indentation,
			p_sintering_element1->mSinteringDisplacement, p_sintering_element1->mSinteringDrivingForce,
			element1, element2, ViscoDampingLocalContactForce);
			p_sintering_element1->mActualNeighbourSinteringDisplacement.push_back(p_sintering_element1->mSinteringDisplacement);  // adding the sintering displacements to vector (only for sintering period (continuum CL)
		}
		else {
			CalculateNormalForcesAfterSintering(LocalElasticContactForce,
				                                  p_sintering_element1->mSinteringDisplacement,
				                                  indentation,
				                                  calculation_area,
				                                  element1,
				                                  element2,
				                                  i_neighbour_count);

			CalculateTangentialForces(OldLocalElasticContactForce,
				                        LocalElasticContactForce,
				                        LocalElasticExtraContactForce,
				                        ViscoDampingLocalContactForce,
				                        LocalCoordSystem,
				                        LocalDeltDisp,
				                        LocalRelVel,
				                        kt_el,
				                        equiv_shear,
				                        contact_sigma,
				                        contact_tau,
				                        indentation,
				                        calculation_area,
				                        failure_criterion_state,
				                        element1,
				                        element2,
				                        i_neighbour_count,
				                        sliding,
				                        r_process_info);
		}

		KRATOS_CATCH("")
	}

	void SinteringContinuum::CalculateNormalForcesAfterSintering(double LocalElasticContactForce[3], // HERTZIAN CL
		                                                           double sintering_displ,
		                                                           double indentation,
		                                                           double calculation_area,
		                                                           SphericContinuumParticle* element1,
		                                                           SphericContinuumParticle* element2,
		                                                           int i_neighbour_count) {

		KRATOS_TRY

		int& failure_type = element1->mIniNeighbourFailureId[i_neighbour_count];

		double minimal_radius;
		double kn = 0;
		indentation = -indentation;

		InitializeContact(element1, element2, indentation, minimal_radius, kn, sintering_displ);

		const double& tension_limit = (*mpProperties)[CONTACT_SIGMA_MIN];
		const double limit_force = tension_limit * calculation_area;

		if (indentation <= 0.0) { // COMPRESSION
			LocalElasticContactForce[2] = -kn * (indentation - sintering_displ) * 0.666666666666666;
		}
		else { // tension
			if (failure_type == 0) {
				LocalElasticContactForce[2] = -kn * (indentation - sintering_displ) * 0.666666666666666;
				if (fabs(LocalElasticContactForce[2]) > limit_force) {
					failure_type = 4; // tension failure
					LocalElasticContactForce[2] = 0.0;
				}
			}
			else LocalElasticContactForce[2] = 0.0;
		}

		KRATOS_CATCH("")
	}

	void SinteringContinuum::InitializeContact(SphericContinuumParticle* const element1,
		                                         SphericContinuumParticle* const element2,
		                                         double& indentation,
		                                         double& minimal_radius,
		                                         double& kn,
		                                         double  sintering_displ) {
		// Get equivalent and minimal Radius
		const double my_radius      = element1->GetRadius();
		const double other_radius   = element2->GetRadius();
		const double radius_sum     = my_radius + other_radius;
		const double radius_sum_inv = 1.0 / radius_sum;
		double equiv_radius         = my_radius * other_radius * radius_sum_inv;
		minimal_radius              = std::min(my_radius, other_radius);

		// Get equivalent Young's Modulus
		const double my_young      = element1->GetYoung();
		const double other_young   = element2->GetYoung();
		const double my_poisson    = element1->GetPoisson();
		const double other_poisson = element2->GetPoisson();
		const double equiv_young   = my_young * other_young / (other_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - other_poisson * other_poisson));

		//const double my_shear_modulus    = 0.5 * my_young / (1.0 + my_poisson);
		//const double other_shear_modulus = 0.5 * other_young / (1.0 + other_poisson);
		//const double equiv_shear         = 1.0 / ((2.0 - my_poisson) / my_shear_modulus + (2.0 - other_poisson) / other_shear_modulus);

		// Normal and Tangent elastic constants
		const double sqrt_equiv_radius = sqrt(equiv_radius);
		const double sqrt_indentation_with_sintering_displ = sqrt(std::abs(indentation - sintering_displ));
		kn = 2.0 * equiv_young * sqrt_equiv_radius * sqrt_indentation_with_sintering_displ;
	}

	//void SinteringContinuum::ComputeParticleRotationalMoments(SphericContinuumParticle* element,
	//                                                          SphericContinuumParticle* neighbor,
	//                                                          double equiv_young,
	//                                                          double distance,
	//                                                          double calculation_area,
	//                                                          double LocalCoordSystem[3][3],
	//                                                          double ElasticLocalRotationalMoment[3],
	//                                                          double ViscoLocalRotationalMoment[3],
	//                                                          double equiv_poisson,
	//                                                          double indentation) {

	//	if (element->Is(DEMThermalFlags::IS_SINTERING)) return;

	//	DEM_KDEM::ComputeParticleRotationalMoments(element,
	//		                                         neighbor,
	//		                                         equiv_young,
	//		                                         distance,
	//		                                         calculation_area,
	//		                                         LocalCoordSystem,
	//	                                           ElasticLocalRotationalMoment,
	//		                                         ViscoLocalRotationalMoment,
	//		                                         equiv_poisson,
	//		                                         indentation);
	//}

} // namespace Kratos
