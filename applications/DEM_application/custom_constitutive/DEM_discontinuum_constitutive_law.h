
#if !defined(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED)
#define  DEM_DISCONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED

/* Project includes */
#include "includes/define.h"
#include "../custom_utilities/AuxiliaryFunctions.h"
#include "../custom_utilities/properties_proxies.h"
#include "includes/serializer.h"
#include "includes/properties.h"
#include "containers/flags.h"

namespace Kratos
{

/**
 * Base class of constitutive laws.
 */
  class /*__declspec( dllexport )*/ DEMDiscontinuumConstitutiveLaw : public Flags
  {

  public:
    
    KRATOS_CLASS_POINTER_DEFINITION( DEMDiscontinuumConstitutiveLaw );

      DEMDiscontinuumConstitutiveLaw();
      
      DEMDiscontinuumConstitutiveLaw( const DEMDiscontinuumConstitutiveLaw &rReferenceDiscontinuumConstitutiveLaw);
      
    void Initialize(const ProcessInfo& rCurrentProcessInfo);
      
    virtual void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;
      
      virtual ~DEMDiscontinuumConstitutiveLaw();
      
    virtual DEMDiscontinuumConstitutiveLaw::Pointer Clone() const;

    virtual void CalculateContactForces(double LocalElasticContactForce[3],
                                        double indentation,
                                        double kn_el,
                                        double LocalDeltDisp[3],
                                        double kt_el,
                                        int& neighbour_failure_id,
                                        double equiv_tg_of_fri_ang);

    virtual void PlasticityAndDamage(double LocalElasticContactForce[3],
                                     double kn_el,
                                     double equiv_young,
                                     double indentation,
                                     double calculation_area,
                                     double radius_sum_i,
                                     double& failure_criterion_state,
                                     double& acumulated_damage,
                                     int i_neighbour_count,
                                     int mapping_new_cont,
                                     int mapping_new_ini,
                                     int time_steps);

    virtual void CalculateViscoDamping(double LocalRelVel[3],
                                       double ViscoDampingLocalContactForce[3],
                                       double indentation,
                                       double equiv_visco_damp_coeff_normal,
                                       double equiv_visco_damp_coeff_tangential,
                                       bool sliding,
                                       int mDampType);
    
   
  private:
    
      friend class Serializer;

      virtual void save(Serializer& rSerializer) const
      {
          KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Flags )
          //rSerializer.save("MyMemberName",myMember);

      }

      virtual void load(Serializer& rSerializer)
      {
          KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Flags )
          //rSerializer.load("MyMemberName",myMember);

    }

};


KRATOS_DEFINE_VARIABLE(DEMDiscontinuumConstitutiveLaw::Pointer, DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER)

} /* namespace Kratos.*/
#endif /* DEM_CONSTITUTIVE_LAW_H_INCLUDED  defined */

