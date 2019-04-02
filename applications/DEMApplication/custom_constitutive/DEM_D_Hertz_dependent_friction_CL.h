//Authors: J. Irazabal (CIMNE)
//   Date: November 2016

#if !defined(DEM_D_DEPENDENT_FRICTION_CL_H_INCLUDED)
#define DEM_D_DEPENDENT_FRICTION_CL_H_INCLUDED

#include <string>
#include <iostream>
#include "DEM_discontinuum_constitutive_law.h"

namespace Kratos {

    class SphericDiscreteContactFeaturesParticle;

    class KRATOS_API(DEM_APPLICATION) DEM_D_Hertz_dependent_friction : public DEMDiscontinuumConstitutiveLaw {

    public:
        using DEMDiscontinuumConstitutiveLaw::CalculateNormalForce;

        KRATOS_CLASS_POINTER_DEFINITION(DEM_D_Hertz_dependent_friction);

        DEM_D_Hertz_dependent_friction() {}

        ~DEM_D_Hertz_dependent_friction() {}

        void Initialize(const ProcessInfo& r_process_info) override;

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) const override;

        void Check(Properties::Pointer pProp) const override;

        std::string GetTypeOfLaw() override;

        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const override;

        void InitializeDependentContact(double equiv_radius,
                                        const double equiv_level_of_fouling,
                                        const double equiv_young,
                                        const double equiv_shear,
                                        const double indentation);

        void DamageContact(SphericDiscreteContactFeaturesParticle* const element1,
                           SphericDiscreteContactFeaturesParticle* const element2,
                           double& equiv_radius,
                           const double equiv_level_of_fouling,
                           const double equiv_young,
                           const double equiv_shear, double& indentation,
                           const double normal_contact_force);

        void InitializeDependentContactWithFEM(double effective_radius,
                                               const double equiv_level_of_fouling,
                                               const double equiv_young,
                                               const double equiv_shear,
                                               const double indentation);

        void DamageContactWithFEM(SphericDiscreteContactFeaturesParticle* const element,
                                  Condition* const wall,
                                  double& effective_radius,
                                  const double equiv_level_of_fouling,
                                  const double equiv_young,
                                  const double equiv_shear,
                                  double& indentation,
                                  const double normal_contact_force);

        void CalculateForces(const ProcessInfo& r_process_info,
                             const double OldLocalElasticContactForce[3],
                             double LocalElasticContactForce[3],
                             double LocalDeltDisp[3],
                             double LocalRelVel[3],
                             double indentation,
                             double previous_indentation,
                             double ViscoDampingLocalContactForce[3],
                             double& cohesive_force,
                             SphericDiscreteContactFeaturesParticle* element1,
                             SphericDiscreteContactFeaturesParticle* element2,
                             bool& sliding,
                             double LocalCoordSystem[3][3]);

        void CalculateForcesWithFEM(ProcessInfo& r_process_info,
                                    const double OldLocalElasticContactForce[3],
                                    double LocalElasticContactForce[3],
                                    double LocalDeltDisp[3],
                                    double LocalRelVel[3],
                                    double indentation,
                                    double previous_indentation,
                                    double ViscoDampingLocalContactForce[3],
                                    double& cohesive_force,
                                    SphericDiscreteContactFeaturesParticle* const element,
                                    Condition* const wall,
                                    bool& sliding);

        double CalculateNormalForce(const double indentation) override;

        double CalculateCohesiveNormalForce(SphericDiscreteContactFeaturesParticle* const element1,
                                            SphericDiscreteContactFeaturesParticle* const element2,
                                            const double indentation);

        double CalculateCohesiveNormalForceWithFEM(SphericDiscreteContactFeaturesParticle* const element,
                                                   Condition* const wall,
                                                   const double indentation);

        void CalculateTangentialForce(const double normal_contact_force,
                                      const double OldLocalElasticContactForce[3],
                                      double LocalElasticContactForce[3],
                                      double ViscoDampingLocalContactForce[3],
                                      const double LocalDeltDisp[3],
                                      bool& sliding,
                                      SphericDiscreteContactFeaturesParticle* const element1,
                                      SphericDiscreteContactFeaturesParticle* const element2,
                                      const double equiv_radius,
                                      const double equiv_young,
                                      double indentation,
                                      double previous_indentation,
                                      double& AuxElasticShearForce,
                                      double& MaximumAdmisibleShearForce);

        void CalculateTangentialForceWithFEM(const double normal_contact_force,
                                             const double OldLocalElasticContactForce[3],
                                             double LocalElasticContactForce[3],
                                             double ViscoDampingLocalContactForce[3],
                                             const double LocalDeltDisp[3],
                                             bool& sliding,
                                             SphericDiscreteContactFeaturesParticle* const element,
                                             Condition* const wall,
                                             const double equiv_radius,
                                             const double equiv_young,
                                             double indentation,
                                             double previous_indentation,
                                             double& AuxElasticShearForce,
                                             double& MaximumAdmisibleShearForce);

        void CalculateViscoDampingForce(double LocalRelVel[3],
                                        double ViscoDampingLocalContactForce[3],
                                        SphericDiscreteContactFeaturesParticle* const element1,
                                        SphericDiscreteContactFeaturesParticle* const element2);

        void CalculateViscoDampingForceWithFEM(double LocalRelVel[3],
                                        double ViscoDampingLocalContactForce[3],
                                        SphericDiscreteContactFeaturesParticle* const element,
                                        Condition* const wall);

        void CalculateElasticEnergyDEM(double& elastic_energy,
                                       double indentation,
                                       double LocalElasticContactForce[3]);

        void CalculateInelasticFrictionalEnergyDEM(double& inelastic_frictional_energy,
                                                   double& AuxElasticShearForce,
                                                   double LocalElasticContactForce[3]);

        void CalculateInelasticViscodampingEnergyDEM(double& inelastic_viscodamping_energy,
                                                     double ViscoDampingLocalContactForce[3],
                                                     double LocalDeltDisp[3]);

        void CalculateElasticEnergyFEM(double& elastic_energy,
                                       double indentation,
                                       double LocalElasticContactForce[3]);

        void CalculateInelasticFrictionalEnergyFEM(double& inelastic_frictional_energy,
                                                   double& AuxElasticShearForce,
                                                   double LocalElasticContactForce[3]);

        void CalculateInelasticViscodampingEnergyFEM(double& inelastic_viscodamping_energy,
                                                     double ViscoDampingLocalContactForce[3],
                                                     double LocalDeltDisp[3]);

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
        }

        virtual void load(Serializer& rSerializer)  override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
        }

    }; //class DEM_D_Hertz_dependent_friction

} /* namespace Kratos.*/
#endif /* DEM_D_DEPENDENT_FRICTION_CL_H_INCLUDED  defined */
