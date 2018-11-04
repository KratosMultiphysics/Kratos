
#if !defined(DEM_KDEMFabric2D_CL_H_INCLUDED)
#define DEM_KDEMFabric2D_CL_H_INCLUDED

#include "dem_kdem_2d_cl.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) DEM_KDEMFabric2D : public DEM_KDEM2D {

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_KDEMFabric2D);

        DEM_KDEMFabric2D() {}

        ~DEM_KDEMFabric2D() {}

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) const override;

        void CalculateContactArea(double radius, double other_radius, double& calculation_area) override;
        
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

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMContinuumConstitutiveLaw)
        }
        
        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMContinuumConstitutiveLaw)
        }
    };
} // namespace Kratos
#endif // DEM_KDEMFabric2D_CL_H_INCLUDED defined
