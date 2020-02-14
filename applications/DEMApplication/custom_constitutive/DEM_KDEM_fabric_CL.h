
#if !defined(DEM_KDEM_FABRIC_CL_H_INCLUDED)
#define  DEM_KDEM_FABRIC_CL_H_INCLUDED

/* Project includes */
#include "DEM_continuum_constitutive_law.h"
#include "DEM_KDEM_CL.h"


namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) DEM_KDEMFabric : public DEM_KDEM {

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_KDEMFabric);

        DEM_KDEMFabric() {}

        ~DEM_KDEMFabric() {}

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) override;

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;

        void AddPoissonContribution(const double equiv_poisson,
                                    double LocalCoordSystem[3][3],
                                    double& normal_force,
                                    double calculation_area,
                                    BoundedMatrix<double, 3, 3>* mSymmStressTensor,
                                    SphericContinuumParticle* element1,
                                    SphericContinuumParticle* element2,
                                    const ProcessInfo& r_process_info,
                                    const int i_neighbor_count,
                                    const double indentation) override;

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMContinuumConstitutiveLaw)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMContinuumConstitutiveLaw)
                    //rSerializer.load("MyMemberName",myMember);
        }

    };

} /* namespace Kratos.*/
#endif /* DEM_KDEM_FABRIC_H_INCLUDED  defined */
