
#if !defined(DEM_DEMPACK_TORQUE_CL_H_INCLUDED)
#define  DEM_DEMPACK_TORQUE_CL_H_INCLUDED

/* Project includes */
#include "DEM_Dempack_CL.h"
//#include "DEM_discontinuum_constitutive_law.h"


namespace Kratos {

    class DEM_Dempack_torque : public DEM_Dempack {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_Dempack_torque);

        DEM_Dempack_torque() {
        }

        void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;

        ~DEM_Dempack_torque() {
        }

        DEMContinuumConstitutiveLaw::Pointer Clone() const;       

        void ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                              SphericContinuumParticle* neighbor,
                                              double equiv_young,
                                              double distance,
                                              double calculation_area,
                                              double LocalCoordSystem[3][3],
                                              array_1d<double, 3>& mContactMoment);
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
#endif /* DEM_DEMPACK_TORQUE_CL_H_INCLUDED  defined */
