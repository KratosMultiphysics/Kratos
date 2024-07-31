/////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu
// Date: Oct 2023
/////////////////////////////////////////////////

#if !defined(DEM_PARALLEL_BOND_BILINEAR_DAMAGE_MIXED_CL_H_INCLUDE)
#define DEM_PARALLEL_BOND_BILINEAR_DAMAGE_MIXED_CL_H_INCLUDE

// Project includes
#include "DEM_parallel_bond_bilinear_damage_CL.h"

namespace Kratos{

    class KRATOS_API(DEM_APPLICATION) DEM_parallel_bond_bilinear_damage_mixed : public DEM_parallel_bond_bilinear_damage {

        typedef DEM_parallel_bond_bilinear_damage BaseClassType;

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_parallel_bond_bilinear_damage_mixed);

        DEM_parallel_bond_bilinear_damage_mixed() {}

        void Check(Properties::Pointer pProp) const override;

        ~DEM_parallel_bond_bilinear_damage_mixed() {}

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;

        double ComputeNormalUnbondedForce(double unbonded_indentation) override;

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
                            double indentation_particle,
                            double calculation_area,
                            double& acumulated_damage,
                            SphericContinuumParticle* element1,
                            SphericContinuumParticle* element2,
                            int i_neighbour_count,
                            int time_steps,
                            bool& sliding,
                            double &equiv_visco_damp_coeff_normal,
                            double &equiv_visco_damp_coeff_tangential,
                            double LocalRelVel[3],
                            double ViscoDampingLocalContactForce[3]) override;

        /*
        double mDamageNormal = 0.0;
        double mDamageTangential = 0.0;
        double mDamageMoment = 0.0;
        const double mDamageThresholdTolerance = 0.999;
        double mDamageReal = 0.0;
        bool mDebugPrintingOption = false;
        double mInitialIndentationForBondedPart = 0.0;*/
    
    protected:

    private:

    };   

} // namespace Kratos

#endif /*DEM_PARALLEL_BOND_BILINEAR_DAMAGE_MIXED_CL_H_INCLUDE defined*/