
#if !defined(DEM_DEMPACK1_CL_H_INCLUDED)
#define  DEM_DEMPACK1_CL_H_INCLUDED

/* Project includes */
#include "DEM_continuum_constitutive_law.h"
#include "../custom_elements/spheric_continuum_particle.h"
#include "DEM_discontinuum_constitutive_law.h"
#include "../custom_elements/spheric_particle.h"
namespace Kratos
{

class DEM_Dempack1:public DEMContinuumConstitutiveLaw
{
public:

    KRATOS_CLASS_POINTER_DEFINITION( DEM_Dempack1 );

    DEM_Dempack1(){}

    DEMContinuumConstitutiveLaw::Pointer Clone() const;

    void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;

    ~DEM_Dempack1(){}

    double mN1;
    double mN2;
    double mN3;
    double mC1;
    double mC2;
    double mC3;
    double mYoungPlastic;
    double mPlasticityLimit;
    double mDamageMaxDisplacementFactor;
    double mTensionLimit;
    double mTauZero;
    double mContactInternalFriccion;
    double mTanContactInternalFriccion;
    double mSinContactInternalFriccion;
    double mCosContactInternalFriccion;
    int mFailureCriterionOption;


    void Initialize(const ProcessInfo& rCurrentProcessInfo);


    void CalculateContactForces(double mRadius,
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



    void PlasticityAndDamage(double LocalElasticContactForce[3],
                             double kn,
                             double equiv_young,
                             double indentation,
                             double corrected_area,
                             double radius_sum_i,
                             double& failure_criterion_state,
                             double& acumulated_damage,
                             int& neighbour_failure_id,
                             int mapping_new_cont,
                             int mapping_new_ini,
                             int time_steps,
                             array_1d<double, 4 > &mHistory_mapping_new_cont,
                             int mIniNeighbourFailureId_mapping_new_ini);



    void EvaluateFailureCriteria(double LocalElasticContactForce[3],
                                 double ShearForceNow,
                                 double calculation_area,
                                 int i_neighbour_count,
                                 double& contact_sigma,
                                 double& contact_tau,
                                 double& failure_criterion_state,
                                 bool& sliding,
                                 int mapping_new_ini,
                                 int mFailureCriterionOption,
                                 double mTauZero,
                                 double mTanContactInternalFriccion,
                                 double mSinContactInternalFriccion,
                                 double mCosContactInternalFriccion,
                                 int mIniNeighbourFailureId_mapping_new_ini,
                                 int& neighbour_failure_id);



private:

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, DEMContinuumConstitutiveLaw )
                //rSerializer.save("MyMemberName",myMember);
    }

    virtual void load(Serializer& rSerializer)
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, DEMContinuumConstitutiveLaw )
                //rSerializer.load("MyMemberName",myMember);
    }

};

} /* namespace Kratos.*/
#endif /* DEM_DEMPACK1_H_INCLUDED  defined */
