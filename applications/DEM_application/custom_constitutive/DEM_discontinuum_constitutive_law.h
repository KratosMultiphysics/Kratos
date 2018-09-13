
#if !defined(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED)
#define  DEM_DISCONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED

/* Project includes */
#include "includes/define.h"
#include "../custom_utilities/AuxiliaryFunctions.h"
#include "includes/serializer.h"
#include "containers/flags.h"

#include "../custom_utilities/GeometryFunctions.h"
#include "../custom_elements/discrete_element.h"
#include "../custom_elements/Particle_Contact_Element.h"
#include "containers/vector_component_adaptor.h"
#include "containers/array_1d.h"


namespace Kratos {
    
    class Properties;
    class SphericParticle; // forward declaration of spheric cont particle

    class KRATOS_API(DEM_APPLICATION) DEMDiscontinuumConstitutiveLaw : public Flags {
    public:

        double mKn;
        double mKt;

        KRATOS_CLASS_POINTER_DEFINITION(DEMDiscontinuumConstitutiveLaw);

        DEMDiscontinuumConstitutiveLaw();

        DEMDiscontinuumConstitutiveLaw(const DEMDiscontinuumConstitutiveLaw &rReferenceDiscontinuumConstitutiveLaw);

        virtual void Initialize(const ProcessInfo& r_process_info);

        virtual void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) const;
        
        virtual std::string GetTypeOfLaw();

        virtual ~DEMDiscontinuumConstitutiveLaw();

        virtual DEMDiscontinuumConstitutiveLaw::Pointer Clone() const;

        virtual void CalculateContactArea(double radius, double other_radius, double &calculation_area);

        virtual void CalculateElasticConstants(double &kn_el,
                double &kt_el,
                double initial_dist,
                double equiv_young,
                double equiv_poisson,
                double calculation_area,
                SphericParticle* element1,
                SphericParticle* element2);

        
        virtual void CalculateElasticEnergy(double& normal_elastic_energy,
                                                                double indentation,
                                                                double& cohesive_force,
                                                                SphericParticle* element1,
                                                                SphericParticle* element2);

        
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
        virtual void InitializeContactWithFEM(SphericParticle* const element, Condition* const wall, const double indentation, const double ini_delta = 0.0);
        
        virtual void GetContactStiffness(SphericParticle* const element1, SphericParticle* const element2, const double ini_delta, double& kn,double& kt);
        
        virtual void CalculateForces(const ProcessInfo& r_process_info,
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
                                                        bool& sliding, double LocalCoordSystem[3][3]);
        
        virtual void CalculateForcesWithFEM( ProcessInfo& r_process_info,
                                            const double OldLocalContactForce[3],
                                            double LocalElasticContactForce[3],
                                            double LocalDeltDisp[3],
                                            double LocalRelVel[3],            
                                            double indentation,
                                            double previous_indentation,
                                            double ViscoDampingLocalContactForce[3],
                                            double& cohesive_force,
                                            SphericParticle* const element,
                                            Condition* const wall,
                                            bool& sliding);
                
        virtual double CalculateNormalForce(const double indentation);
        virtual double CalculateNormalForce(SphericParticle* const element1, SphericParticle* const element2, const double indentation, double LocalCoordSystem[3][3]);
        virtual double CalculateNormalForce(SphericParticle* const element, Condition* const wall, const double indentation); 
        virtual double CalculateCohesiveNormalForce(SphericParticle * const element1, SphericParticle * const element2, const double indentation);
        virtual double CalculateCohesiveNormalForceWithFEM(SphericParticle* const element, Condition* const wall, const double indentation);       
        virtual double LocalPeriod(const int i, SphericParticle* element1,SphericParticle* element2);


    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags)
                    //rSerializer.save("MyMemberName",myMember);

        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags)
                    //rSerializer.load("MyMemberName",myMember);
        }
    };

    //This definition is done here to avoid recursive inclusion of header files
    KRATOS_DEFINE_APPLICATION_VARIABLE(DEM_APPLICATION, DEMDiscontinuumConstitutiveLaw::Pointer, DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER)

} /* namespace Kratos.*/
#endif /* DEM_CONSTITUTIVE_LAW_H_INCLUDED  defined */

