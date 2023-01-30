
#if !defined(DEM_KDEM_SOFT_TORQUE_H_INCLUDED)
#define  DEM_KDEM_SOFT_TORQUE_H_INCLUDED

/* Project includes */
#include "DEM_KDEM_CL.h"


namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) DEM_KDEM_soft_torque : public DEM_KDEM {

        typedef DEM_KDEM BaseClassType;

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_KDEM_soft_torque);

        DEM_KDEM_soft_torque() {}

        ~DEM_KDEM_soft_torque() {}

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;

        virtual void ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                                      SphericContinuumParticle* neighbor,
                                                      double equiv_young,
                                                      double distance,
                                                      double calculation_area,
                                                      double LocalCoordSystem[3][3],
                                                      double ElasticLocalRotationalMoment[3],
                                                      double ViscoLocalRotationalMoment[3],
                                                      double equiv_poisson,
                                                      double indentation,
                                                      double LocalElasticContactForce[3]) override;

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override{
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEM_KDEM)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) override{
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEM_KDEM)
                    //rSerializer.load("MyMemberName",myMember);
        }

    };

} /* namespace Kratos.*/
#endif /* DEM_KDEM_SOFT_TORQUE_H_INCLUDED  defined */
