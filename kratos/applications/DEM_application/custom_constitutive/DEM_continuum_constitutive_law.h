
#if !defined(DEM_CONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED)
#define  DEM_CONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED

/* Project includes */
#include "includes/define.h"
#include "../custom_utilities/AuxiliaryFunctions.h"
#include "../custom_utilities/properties_proxies.h"
#include "includes/serializer.h"
#include "includes/properties.h"
#include "containers/flags.h"

#include "custom_utilities/GeometryFunctions.h"
#include "../custom_elements/discrete_element.h"
#include "../custom_elements/Particle_Contact_Element.h"
#include "containers/vector_component_adaptor.h"
#include "containers/array_1d.h"
namespace Kratos
{

/**
 * Base class of constitutive laws.
 */
  class /*__declspec( dllexport )*/ DEMContinuumConstitutiveLaw : public Flags
  {

  public:
    
    KRATOS_CLASS_POINTER_DEFINITION( DEMContinuumConstitutiveLaw );

      DEMContinuumConstitutiveLaw();
      
      DEMContinuumConstitutiveLaw( const DEMContinuumConstitutiveLaw &rReferenceContinuumConstitutiveLaw);
      
    virtual void Initialize(const ProcessInfo& rCurrentProcessInfo);
      
    virtual void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;
      
      virtual ~DEMContinuumConstitutiveLaw();
      
    virtual DEMContinuumConstitutiveLaw::Pointer Clone() const;

    virtual void CalculateContactForces(double mRadius,
                                        double mSqrtOfRealMass,
                                        double other_radius,
                                        double otherSqrtMass,
                                        double distance,
                                        double initial_delta,
                                        int& neighbour_failure_id,
                                        ProcessInfo& rCurrentProcessInfo,
                                        PropertiesProxy *myProperties,
                                        PropertiesProxy *neighbourProperties,
                                        int mapping_new_ini,
                                        int mapping_new_cont,
                                        unsigned int i_neighbour_count,
                                        double LocalElasticContactForce[3],
                                        double ViscoDampingLocalContactForce[3],
                                        double LocalDeltDisp[3],
                                        Vector mcont_ini_neigh_area,
                                        array_1d<double, 4 > &mHistory_mapping_new_cont,
                                        double mDempack_damping,
                                        int mDampType,
                                        int mIniNeighbourFailureId_mapping_new_ini,
                                        double LocalCoordSystem[3][3],
                                        double RelVel[3]); //FF


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


KRATOS_DEFINE_VARIABLE(DEMContinuumConstitutiveLaw::Pointer, DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER)

} /* namespace Kratos.*/
#endif /* DEM_CONSTITUTIVE_LAW_H_INCLUDED  defined */

