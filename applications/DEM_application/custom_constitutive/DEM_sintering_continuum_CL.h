//Authors: S. Nosewicz (IPPT PAN, Warsaw, Poland) and M.A. Celigueta (CIMNE)
// Date: April 2015


#if !defined(DEM_sintering_continuum_CL_H_INCLUDED)
#define  DEM_sintering_continuum_CL_H_INCLUDED

/* Project includes */
#include "DEM_KDEM_CL.h"
#include "DEM_continuum_constitutive_law.h"


namespace Kratos {

	class DEM_sintering_continuum : public DEM_KDEM {
	public:

		KRATOS_CLASS_POINTER_DEFINITION(DEM_sintering_continuum);

		DEM_sintering_continuum() {
		}

		void Initialize();

		void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;

		~DEM_sintering_continuum() {
		}

		DEMContinuumConstitutiveLaw::Pointer Clone() const;

		void CalculateSinteringForces(const ProcessInfo& r_process_info,
			const double OldLocalElasticContactForce[3],
			double LocalElasticContactForce[3],
			const double rel_vel,
			const double indentation,
			double& sintering_displ,
			double& sinter_driv_force,
			SphericContinuumParticle* element1,
			SphericContinuumParticle* element2,
			double ViscoDampingLocalContactForce[3]);

		void CalculateNormalForces(const ProcessInfo& r_process_info,
			const double OldLocalElasticContactForce[3],
			double LocalElasticContactForce[3],
			const double rel_vel,
			const double indentation,
			double& sintering_displ,
			double& sinter_driv_force,
			SphericContinuumParticle* element1,
			SphericContinuumParticle* element2,
			double ViscoDampingLocalContactForce[3]);

		void CalculateForcesOfSintering(const ProcessInfo& r_process_info,
			const double OldLocalElasticContactForce[3],
			double LocalElasticContactForce[3],
			const double rel_vel,
			const double indentation,
			double& sintering_displ,
			double& sinter_driv_force,
			SphericContinuumParticle* element1,
			SphericContinuumParticle* element2,
			double ViscoDampingLocalContactForce[3]);

		void CalculateDampingCoeff(double& equiv_visco_damp_coeff_normal,
			double kn,
			SphericContinuumParticle* const element1,
			SphericContinuumParticle* const element2);

		void CalculateForces(const ProcessInfo& r_process_info, // After sintering, Hertzian Continuum CL
			double LocalElasticContactForce[3],
			double LocalDeltDisp[3],
			const double kn_el,
			double kt_el,
			double& contact_sigma,
			double& contact_tau,
			double& failure_criterion_state,
			double equiv_young,
			double indentation,
			double calculation_area,
			double& acumulated_damage,
			SphericContinuumParticle* element1,
			SphericContinuumParticle* element2,
			int i_neighbour_count,
			int time_steps,
			bool& sliding,
			int search_control,
			vector<int>& search_control_vector);

		void CalculateNormalForcesAfterSintering(double LocalElasticContactForce[3], // HERTZIAN CL
			const double kn_el,
			double indentation,
			double calculation_area,
			SphericContinuumParticle* element1,
			SphericContinuumParticle* element2,
			int i_neighbour_count);

		void InitializeContact(SphericContinuumParticle* const element1, SphericContinuumParticle* const element2, double& indentation,
		     double& minimal_radius, double& kn, double sintering_displ);


		




		

	private:

		friend class Serializer;

		virtual void save(Serializer& rSerializer) const {
			KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMContinuumConstitutiveLaw)
				//rSerializer.save("MyMemberName",myMember);
		}

		virtual void load(Serializer& rSerializer) {
			KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMContinuumConstitutiveLaw)
				//rSerializer.load("MyMemberName",myMember);
		}

	};

} /* namespace Kratos.*/
#endif /* DEM_KDEM_H_INCLUDED  defined */


