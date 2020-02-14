//Authors: S. Nosewicz (IPPT PAN, Warsaw, Poland) and M.A. Celigueta (CIMNE)
// Date: April 2015


#if !defined(DEM_sintering_continuum_CL_H_INCLUDED)
#define  DEM_sintering_continuum_CL_H_INCLUDED

/* Project includes */
#include "DEM_KDEM_CL.h"
#include "DEM_continuum_constitutive_law.h"


namespace Kratos {

	class KRATOS_API(DEM_APPLICATION) DEM_sintering_continuum : public DEM_KDEM {
	public:

		KRATOS_CLASS_POINTER_DEFINITION(DEM_sintering_continuum);

		DEM_sintering_continuum() {
		}

		void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) override;

		~DEM_sintering_continuum() {
		}

		DEMContinuumConstitutiveLaw::Pointer Clone() const override;

                void GetContactArea(const double radius,
                                    const double other_radius,
                                    const Vector& vector_of_initial_areas,
                                    const int neighbour_position,
                                    double& calculation_area) override;

                void CalculateElasticConstants(double& kn_el, double& kt_el, double initial_dist, double equiv_young,
                                             double equiv_poisson, double calculation_area, SphericContinuumParticle* element1, SphericContinuumParticle* element2) override;

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

		void CalculateForces(const ProcessInfo& r_process_info,
                                    double OldLocalElasticContactForce[3],
                                    double LocalElasticContactForce[3],
                                    double LocalElasticExtraContactForce[3],
                                    double LocalCoordSystem[3][3],
                                    double LocalDeltDisp[3],
                                    const double kn_el,
                                    const double kt_el,
                                    double& contact_sigma,
                                    double& contact_tau,
                                    double& failure_criterion_state,
                                    double equiv_young,
                                    double equiv_shear,
                                    double indentation,
                                    double calculation_area,
                                    double& acumulated_damage,
                                    SphericContinuumParticle* element1,
                                    SphericContinuumParticle* element2,
                                    int i_neighbour_count,
                                    int time_steps,
                                    bool& sliding,
                                    int search_control,
                                    DenseVector<int>& search_control_vector,
                                    double &equiv_visco_damp_coeff_normal,
                                    double &equiv_visco_damp_coeff_tangential,
                                    double LocalRelVel[3],
                                    double ViscoDampingLocalContactForce[3]) override;

		void CalculateNormalForcesAfterSintering(double LocalElasticContactForce[3], // HERTZIAN CL
			const double kn_el,
			double indentation,
			double calculation_area,
			SphericContinuumParticle* element1,
			SphericContinuumParticle* element2,
			int i_neighbour_count);

		void InitializeContact(SphericContinuumParticle* const element1, SphericContinuumParticle* const element2, double& indentation,
		     double& minimal_radius, double& kn, double sintering_displ);


                void ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                                    SphericContinuumParticle* neighbor,
                                                    double equiv_young,
                                                    double distance,
                                                    double calculation_area,
                                                    double LocalCoordSystem[3][3],
                                                    double ElasticLocalRotationalMoment[3],
                                                    double ViscoLocalRotationalMoment[3],
                                                    double equiv_poisson,
                                                    double indentation) override;



	private:

		friend class Serializer;

		virtual void save(Serializer& rSerializer) const override{
			KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMContinuumConstitutiveLaw)
				//rSerializer.save("MyMemberName",myMember);
		}

		virtual void load(Serializer& rSerializer) override{
			KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMContinuumConstitutiveLaw)
				//rSerializer.load("MyMemberName",myMember);
		}

	};

} /* namespace Kratos.*/
#endif /* DEM_KDEM_H_INCLUDED  defined */


