/////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu
// Date: Oct 2023
/////////////////////////////////////////////////

#if !defined(DEM_PARALLEL_BOND_BILINEAR_DAMAGE_CL_H_INCLUDE)
#define DEM_PARALLEL_BOND_BILINEAR_DAMAGE_CL_H_INCLUDE

// Project includes
#include "DEM_parallel_bond_CL.h"

namespace Kratos{

    class KRATOS_API(DEM_APPLICATION) DEM_parallel_bond_bilinear_damage : public DEM_parallel_bond {

        typedef DEM_parallel_bond BaseClassType;

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_parallel_bond_bilinear_damage);

        DEM_parallel_bond_bilinear_damage() {}

        void Check(Properties::Pointer pProp) const override;

        ~DEM_parallel_bond_bilinear_damage() {}

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;
        double ComputeNormalUnbondedForce(double unbonded_indentation) override;
        double GetBondKn(double bond_equiv_young, double calculation_area, double distance) override;

        void CalculateNormalForces(double LocalElasticContactForce[3],
                            const double kn_el,
                            double equiv_young,
                            double indentation,
                            double indentation_particle,
                            double calculation_area,
                            double& acumulated_damage,
                            SphericContinuumParticle* element1,
                            SphericContinuumParticle* element2,
                            int i_neighbour_count,
                            int time_steps,
                            const ProcessInfo& r_process_info,
                            double& contact_sigma) override;
            
        void CalculateViscoDampingCoeff(double &equiv_visco_damp_coeff_normal,
                                double &equiv_visco_damp_coeff_tangential,
                                SphericContinuumParticle* element1,
                                SphericContinuumParticle* element2,
                                const double kn_el,
                                const double kt_el) override;

        void CalculateTangentialForces(double OldLocalElasticContactForce[3],
                                            double LocalElasticContactForce[3],
                                            double LocalElasticExtraContactForce[3],
                                            double ViscoDampingLocalContactForce[3],
                                            double LocalCoordSystem[3][3],
                                            double LocalDeltDisp[3],
                                            double LocalRelVel[3],
                                            const double kt_el,
                                            const double equiv_shear,
                                            double& contact_sigma,
                                            double& contact_tau,
                                            double indentation_particle,
                                            double calculation_area,
                                            double& failure_criterion_state,
                                            SphericContinuumParticle* element1,
                                            SphericContinuumParticle* element2,
                                            int i_neighbour_count,
                                            bool& sliding,
                                            const ProcessInfo& r_process_info) override;

        void CheckFailure(const int i_neighbour_count, 
                            SphericContinuumParticle* element1, 
                            SphericContinuumParticle* element2,
                            double& contact_sigma,
                            double& contact_tau, 
                            double LocalElasticContactForce[3],
                            double ViscoDampingLocalContactForce[3],
                            double ElasticLocalRotationalMoment[3],
                            double ViscoLocalRotationalMoment[3]) override;     

        double mDamageNormal = 0.0;
        double mDamageTangential = 0.0;
        double mDamageMoment = 0.0;
        const double mDamageThresholdTolerance = 0.9999;
        double mDamageReal = 0.0;
        bool mDebugPrintingOption = false;
        double mInitialIndentationForBondedPart = 0.0;
    
    protected:

    private:

    };   

} // namespace Kratos

#endif /*DEM_PARALLEL_BOND_BILINEAR_DAMAGE_CL_H_INCLUDE defined*/