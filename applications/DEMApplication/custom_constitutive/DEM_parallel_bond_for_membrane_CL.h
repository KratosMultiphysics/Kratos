/////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu
// Date: Feb 2023
/////////////////////////////////////////////////

#if !defined(DEM_PARALLEL_BOND_FOR_MEMBRANE_CL_H_INCLUDE)
#define DEM_PARALLEL_BOND_FOR_MEMBRANE_CL_H_INCLUDE

// Project includes
#include "DEM_parallel_bond_CL.h"

namespace Kratos{

    class KRATOS_API(DEM_APPLICATION) DEM_parallel_bond_for_membrane : public DEM_parallel_bond {

        typedef DEMContinuumConstitutiveLaw BaseClassType;

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_parallel_bond_for_membrane);

        DEM_parallel_bond_for_membrane() {}

        ~DEM_parallel_bond_for_membrane() {}

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;

        void CalculateElasticConstants(double& kn_el, double& kt_el, double initial_dist, double equiv_young,
                                    double equiv_poisson, double calculation_area, SphericContinuumParticle* element1, SphericContinuumParticle* element2, double indentation) override;
        void InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) override;

        double ComputeNormalUnbondedForce(double unbonded_indentation) override;

        void CalculateUnbondedViscoDampingForce(double LocalRelVel[3],
                                                double UnbondedViscoDampingLocalContactForce[3],
                                                SphericParticle* const element1,
                                                SphericParticle* const element2) override; 

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

    protected:

    private:

    };   

} // namespace Kratos

#endif /*DEM_PARALLEL_BOND_FOR_MEMBRANE_CL_H_INCLUDE defined*/