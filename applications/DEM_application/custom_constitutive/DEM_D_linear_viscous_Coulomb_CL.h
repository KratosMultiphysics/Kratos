
#if !defined(DEM_D_LINEAR_VISCOUS_COULOMB_CL_H_INCLUDED)
#define  DEM_D_LINEAR_VISCOUS_COULOMB_CL_H_INCLUDED

/* Project includes */
#include "../custom_elements/spheric_continuum_particle.h"
#include "DEM_discontinuum_constitutive_law.h"
#include "../custom_elements/spheric_particle.h"
namespace Kratos
{

class DEM_D_linear_viscous_Coulomb:public DEMDiscontinuumConstitutiveLaw
{
public:

    KRATOS_CLASS_POINTER_DEFINITION( DEM_D_linear_viscous_Coulomb );

    DEM_D_linear_viscous_Coulomb(){}

    DEMDiscontinuumConstitutiveLaw::Pointer Clone() const;

    void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;

    ~DEM_D_linear_viscous_Coulomb(){}

    void Initialize(const ProcessInfo& rCurrentProcessInfo);

    virtual void NormalForceCalculation(double LocalElasticContactForce[3], double kn, double indentation);

    virtual void CalculateContactForces(double LocalElasticContactForce[3],
                                        double indentation,
                                        double kn_el,
                                        double LocalDeltDisp[3],
                                        double kt_el,
                                        int& neighbour_failure_id,
                                        double equiv_tg_of_fri_ang);



private:

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, DEMDiscontinuumConstitutiveLaw )
                //rSerializer.save("MyMemberName",myMember);
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, DEMDiscontinuumConstitutiveLaw )
                //rSerializer.load("MyMemberName",myMember);
    }
};

} /* namespace Kratos.*/
#endif /* DEM_D_LINEAR_VISCOUS_COULOMB_CL_H_INCLUDED  defined */
