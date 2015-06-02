
#if !defined(DEM_D_HERTZ_VISCOUS_COULOMB_CL_H_INCLUDED)
#define  DEM_D_HERTZ_VISCOUS_COULOMB_CL_H_INCLUDED

/* Project includes */
//#include "../custom_elements/spheric_continuum_particle.h"
#include "DEM_discontinuum_constitutive_law.h"

namespace Kratos {

    class DEM_D_Hertz_viscous_Coulomb : public DEMDiscontinuumConstitutiveLaw {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_D_Hertz_viscous_Coulomb);

        DEM_D_Hertz_viscous_Coulomb() {
        }
      
        void Initialize(const ProcessInfo& rCurrentProcessInfo);         

        void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;

        ~DEM_D_Hertz_viscous_Coulomb() {
        }

        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const;

        void CalculateForces(double LocalElasticContactForce[3],
                double LocalDeltDisp[3],
                double kn_el,
                double kt_el,
                double indentation,
                double& failure_criterion_state,
                bool& sliding,
                SphericParticle* element1,
                SphericParticle* element2,
                int &mNeighbourFailureId_count,
                double mapping_new_cont);
              
        void CalculateNormalForceHertz(double LocalElasticContactForce[3], const double kn_el, const double indentation);
        
        void CalculateTangentialForceLinear(double LocalElasticContactForce[3],
                double LocalDeltDisp[3],
                const double kt_el,
                const double indentation,
                double& failure_criterion_state,
                bool& sliding,
                SphericParticle* element1,
                SphericParticle* element2,
                int &mNeighbourFailureId_count,
                double mapping_new_cont);

        void CalculateViscoDamping(double LocalRelVel[3],
                double ViscoDampingLocalContactForce[3],
                double indentation,
                double equiv_visco_damp_coeff_normal,
                double equiv_visco_damp_coeff_tangential,
                bool sliding,
                int mDampType);
        
        
        
        void InitializeContact(SphericParticle* const element1, SphericParticle* const element2);  
        void InitializeContactWithFEM(SphericParticle* const element, DEMWall* const wall, const double ini_delta=0.0);
        
        double CalculateNormalForce(const double indentation,
                                    SphericParticle* const element1,
                                    SphericParticle* const element2);
        double CalculateNormalForceWithFEM(const double indentation,
                                    SphericParticle* const element,
                                    DEMWall* const wall);
        
        double CalculateCohesiveNormalForce(SphericParticle* const element1,
                                            SphericParticle* const element2);
        double CalculateCohesiveNormalForceWithFEM(SphericParticle* const element,
                                            DEMWall* const wall);
        
        void CalculateTangentialForce(const double normal_force,
                                      double LocalElasticContactForce[3],
                                      const double LocalDeltDisp[3],            
                                      bool& sliding,
                                      SphericParticle* const element1,
                                      SphericParticle* const element2);
        void CalculateTangentialForceWithFEM(const double normal_force,
                                                    double LocalElasticContactForce[3],
                                                    const double LocalDeltDisp[3],            
                                                    bool& sliding,
                                                    SphericParticle* const element,
                                                    DEMWall* const wall,
                                                    double indentation);
        
        void CalculateViscoDampingForce(double LocalRelVel[3],
                                        double ViscoDampingLocalContactForce[3],
                                        bool sliding,
                                        SphericParticle* const element,
                                        SphericParticle* const element2);                                        
        void CalculateViscoDampingForceWithFEM(double LocalRelVel[3],
                                        double ViscoDampingLocalContactForce[3],
                                        bool sliding,
                                        SphericParticle* const element,
                                        DEMWall* const wall,
                                        double indentation);
                            
    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
                    //rSerializer.load("MyMemberName",myMember);
        }
    };

} /* namespace Kratos.*/
#endif /* DEM_D_HERTZ_VISCOUS_COULOMB_CL_H_INCLUDED  defined */
