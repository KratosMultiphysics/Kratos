
#if !defined(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED)
#define  DEM_DISCONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED

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
#include "../custom_conditions/dem_wall.h"

namespace Kratos {
    /**
     * Base class of constitutive laws.
     */

    class SphericParticle; // forward declaration of spheric cont particle

    class /*__declspec( dllexport )*/ DEMDiscontinuumConstitutiveLaw : public Flags {
    public:

        double mKn;
        double mKt;

        KRATOS_CLASS_POINTER_DEFINITION(DEMDiscontinuumConstitutiveLaw);

        DEMDiscontinuumConstitutiveLaw();

        DEMDiscontinuumConstitutiveLaw(const DEMDiscontinuumConstitutiveLaw &rReferenceDiscontinuumConstitutiveLaw);

        virtual void Initialize(const ProcessInfo& rCurrentProcessInfo);

        virtual void SetConstitutiveLawInProperties(Properties::Pointer pProp) const;
        
        virtual std::string GetTypeOfLaw();

        virtual ~DEMDiscontinuumConstitutiveLaw();

        virtual DEMDiscontinuumConstitutiveLaw::Pointer Clone() const;

        virtual void CalculateContactArea(double radius, double other_radius, double &calculation_area);

        virtual void CalculateElasticConstants(double &kn_el,
                double &kt_el,
                double initial_dist,
                double equiv_young,
                double equiv_poisson,
                double calculation_area);


        virtual void CalculateForces(ProcessInfo& rCurrentProcessInfo,
                                     double LocalElasticContactForce[3],
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
        
        virtual void CalculateElasticEnergy(double& normal_elastic_energy,
                                                                double indentation,
                                                                double& cohesive_force,
                                                                SphericParticle* element1,
                                                                SphericParticle* element2);


        virtual void CalculateNormalForceLinear(double LocalElasticContactForce[3], const double kn_el, const double indentation);
        virtual void CalculateTangentialForceLinear(double LocalElasticContactForce[3],
                double LocalDeltDisp[3],
                const double kt_el,
                const double indentation,
                double& failure_criterion_state,
                bool& sliding,
                SphericParticle* element1,
                SphericParticle* element2,
                int &mNeighbourFailureId_count,
                double mapping_new_cont);
        virtual void CalculateNormalForceHertz(double LocalElasticContactForce[3], const double kn_el, const double indentation);
        virtual void CalculateViscoDamping(double LocalRelVel[3],
                double ViscoDampingLocalContactForce[3],
                double indentation,
                double equiv_visco_damp_coeff_normal,
                double equiv_visco_damp_coeff_tangential,
                bool sliding);

        virtual void CalculateViscoDampingCoeff(double &equiv_visco_damp_coeff_normal,
                double &equiv_visco_damp_coeff_tangential,
                SphericParticle* element1,
                SphericParticle* element2,
                double kn_el,
                double kt_el);

        virtual void InitializeContact(SphericParticle * const element1, SphericParticle * const element2, const double ini_delta = 0.0);
        virtual void InitializeContactWithFEM(SphericParticle* const element, DEMWall* const wall, const double ini_delta=0.0);
        
        virtual void GetContactStiffness(SphericParticle* const element1, SphericParticle* const element2, const double ini_delta, double& kn,double& kt);
        
        virtual void CalculateForces(ProcessInfo& rCurrentProcessInfo,
                                                        const double OldLocalContactForce[3],
                                                        double LocalElasticContactForce[3],
                                                        double LocalDeltDisp[3],
                                                        double LocalRelVel[3],            
                                                        double indentation,
                                                        double previous_indentation,
                                                        double ViscoDampingLocalContactForce[3],
                                                        double& cohesive_force,
                                                        SphericParticle* element1,
                                                        SphericParticle* element2,
                                                        bool& sliding);
        
        virtual void CalculateForcesWithFEM( ProcessInfo& rCurrentProcessInfo,
                                            const double OldLocalContactForce[3],
                                            double LocalElasticContactForce[3],
                                            double LocalDeltDisp[3],
                                            double LocalRelVel[3],            
                                            double indentation,
                                            double previous_indentation,
                                            double ViscoDampingLocalContactForce[3],
                                            double& cohesive_force,
                                            SphericParticle* const element,
                                            DEMWall* const wall,
                                            bool& sliding);
        
        virtual double CalculateNormalForce(const double indentation, SphericParticle * const element1, SphericParticle * const element2);
        
        virtual double CalculateNormalForce(const double indentation);
    
        virtual double CalculateNormalForceWithFEM(const double indentation, SphericParticle* const element, DEMWall* const wall);
        
        virtual void CalculateTangentialForce(const double normal_force,
                                            double LocalElasticContactForce[3],
                                            const double LocalDeltDisp[3],
                                            bool& sliding,
                                            SphericParticle * const element1,
                                            SphericParticle * const element2);
        virtual void CalculateTangentialForceWithFEM(const double OldLocalContactForce[3],
                                                    double LocalElasticContactForce[3],
                                                    const double LocalDeltDisp[3],            
                                                    bool& sliding,
                                                    SphericParticle* const element,
                                                    DEMWall* const wall,
                                                    double indentation,
                                                    double previous_indentation);
        
        virtual double CalculateCohesiveNormalForce(SphericParticle * const element1, SphericParticle * const element2, const double indentation);
        virtual double CalculateCohesiveNormalForceWithFEM(SphericParticle* const element, DEMWall* const wall, const double indentation);

        virtual void CalculateViscoDampingForce(double LocalRelVel[3],
                                                double ViscoDampingLocalContactForce[3],
                                                SphericParticle * const element1,
                                                SphericParticle * const element2);                                                
        virtual void CalculateViscoDampingForceWithFEM(double LocalRelVel[3],
                                                    double ViscoDampingLocalContactForce[3],
                                                    bool sliding,
                                                    SphericParticle* const element,
                                                    DEMWall* const wall,
                                                    double indentation);
        
    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags)
                    //rSerializer.save("MyMemberName",myMember);

        }

        virtual void load(Serializer& rSerializer) {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags)
                    //rSerializer.load("MyMemberName",myMember);
        }
    };

    KRATOS_DEFINE_VARIABLE(DEMDiscontinuumConstitutiveLaw::Pointer, DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER)

} /* namespace Kratos.*/
#endif /* DEM_CONSTITUTIVE_LAW_H_INCLUDED  defined */

