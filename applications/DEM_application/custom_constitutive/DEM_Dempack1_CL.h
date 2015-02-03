
#if !defined(DEM_DEMPACK1_CL_H_INCLUDED)
#define  DEM_DEMPACK1_CL_H_INCLUDED

/* Project includes */
#include "DEM_continuum_constitutive_law.h"
//#include "DEM_discontinuum_constitutive_law.h"


namespace Kratos
{

class DEM_Dempack1:public DEMContinuumConstitutiveLaw
{
public:

    KRATOS_CLASS_POINTER_DEFINITION( DEM_Dempack1 );

    DEM_Dempack1(){}

    //DEMContinuumConstitutiveLaw(const DEMContinuumConstitutiveLaw &rReferenceContinuumConstitutiveLaw);

    void Initialize(const ProcessInfo& rCurrentProcessInfo);

    void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;

    ~DEM_Dempack1(){}

    DEMContinuumConstitutiveLaw::Pointer Clone() const;

    void CalculateContactForces(double mRadius,
                                double mRealMass,
                                double other_radius,
                                double otherMass,
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
                                array_1d<double, 6 > &mHistory_mapping_new_cont,
                                double mDempack_damping,
                                int mDampType,
                                int mIniNeighbourFailureId_mapping_new_ini,
                                double LocalCoordSystem[3][3],
                                double RelVel[3]); //FF



    void PlasticityAndDamage1D(double LocalElasticContactForce[3],
                const double kn_el,
                double equiv_young,
                double indentation,
                double calculation_area,
                double radius_sum_i,
                double& failure_criterion_state,
                double& acumulated_damage,
                int i_neighbour_count,
                int mapping_new_cont,
                int mapping_new_ini,
                const double mN1,
                const double mN2,
                const double mN3,
                const double mYoungPlastic,
                const double mPlasticityLimit,
                const double mC1,
                const double mC2,
                const double mC3,
                const double mTensionLimit,
                const double mDamageMaxDisplacementFactor,
                array_1d <double, 6> &mHistory_mapping_new_cont,
                int &mNeighbourFailureId_i_neighbour_count,
                int &mIniNeighbourFailureId_mapping_new_ini,
                int time_steps);



    void EvaluateFailureCriteria(
                const double contact_sigma,
                const double contact_tau,
                double& failure_criterion_state,
                bool& sliding,
                const int FailureCriterionOption,
                const double TauZero,
                const double TanContactInternalFriccion,
                const double SinContactInternalFriccion,
                const double CosContactInternalFriccion,
                int& NeighbourFailureId_i_neighbour_count,
                int& IniNeighbourFailureId_mapping_new_ini,
                const double TensionLimit);



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
