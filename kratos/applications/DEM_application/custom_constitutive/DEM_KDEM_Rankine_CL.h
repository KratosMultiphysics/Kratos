
#if !defined(DEM_KDEM_RANKINE_CL_H_INCLUDED)
#define  DEM_KDEM_RANKINE_CL_H_INCLUDED

/* Project includes */
#include "DEM_continuum_constitutive_law.h"
#include "DEM_KDEM_CL.h"


namespace Kratos {

    class DEM_KDEM_Rankine : public DEM_KDEM {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_KDEM_Rankine);

        DEM_KDEM_Rankine() {
        }    

        void SetConstitutiveLawInProperties(Properties::Pointer pProp) const override;

        ~DEM_KDEM_Rankine() {
        }

        DEMContinuumConstitutiveLaw::Pointer Clone() const override;
        
        void CheckFailure(const int i_neighbour_count, SphericContinuumParticle* element1, SphericContinuumParticle* element2) override;

        void CalculateNormalForces(double LocalElasticContactForce[3],
                const double kn_el,
                double equiv_young,
                double indentation,
                double calculation_area,
                double& acumulated_damage,
                SphericContinuumParticle* element1,
                SphericContinuumParticle* element2,
                int i_neighbour_count,
                int time_steps) override;

        void CalculateTangentialForces(double OldLocalElasticContactForce[3],
                double LocalElasticContactForce[3],
                double LocalElasticExtraContactForce[3],
                double LocalCoordSystem[3][3],             
                double LocalDeltDisp[3],                
                const double kt_el,
                const double equiv_shear,
                double& contact_sigma,
                double& contact_tau,
                double indentation,
                double calculation_area,
                double& failure_criterion_state,
                SphericContinuumParticle* element1,
                SphericContinuumParticle* element2,
                int i_neighbour_count,
                bool& sliding,
                int search_control,
                vector<int>& search_control_vector,
                const ProcessInfo& r_process_info) override;                

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
#endif /* DEM_KDEM_RANKINE_H_INCLUDED  defined */
