//Authors: J. Irazabal (CIMNE)
//   Date: November 2016

#if !defined(DEM_D_CONICAL_DAMAGE_CL_H_INCLUDED)
#define DEM_D_CONICAL_DAMAGE_CL_H_INCLUDED

#include <string>
#include <iostream>

// Project includes
#include "DEM_D_Hertz_viscous_Coulomb_CL.h"
#include "custom_elements/spheric_particle.h"
#include "custom_elements/contact_info_spheric_particle.h"

namespace Kratos {

    class SphericParticle;

    class KRATOS_API(DEM_APPLICATION) DEM_D_Conical_damage : public DEM_D_Hertz_viscous_Coulomb {

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_D_Conical_damage);

        DEM_D_Conical_damage() {}

        ~DEM_D_Conical_damage() {}


        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) override;

        void Check(Properties::Pointer pProp) const override;

        std::string GetTypeOfLaw() override;

        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const override;

        void InitializeDependentContact(double equiv_radius,
                                        const double equiv_level_of_fouling,
                                        const double equiv_young,
                                        const double equiv_shear,
                                        const double indentation);

        void DamageContact(ContactInfoSphericParticle* const element1,
                           ContactInfoSphericParticle* const element2,
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

        void DamageContactWithFEM(ContactInfoSphericParticle* const element,
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
                             SphericParticle* element1,
                             SphericParticle* element2,
                             bool& sliding,
                             double LocalCoordSystem[3][3]) override;

        void CalculateForcesWithFEM(ProcessInfo& r_process_info,
                                    const double OldLocalElasticContactForce[3],
                                    double LocalElasticContactForce[3],
                                    double LocalDeltDisp[3],
                                    double LocalRelVel[3],
                                    double indentation,
                                    double previous_indentation,
                                    double ViscoDampingLocalContactForce[3],
                                    double& cohesive_force,
                                    SphericParticle* const element,
                                    Condition* const wall,
                                    bool& sliding) override;

        void CalculateTangentialForce(const double normal_contact_force,
                                      const double OldLocalElasticContactForce[3],
                                      double LocalElasticContactForce[3],
                                      double ViscoDampingLocalContactForce[3],
                                      const double LocalDeltDisp[3],
                                      bool& sliding,
                                      ContactInfoSphericParticle* const element1,
                                      ContactInfoSphericParticle* const element2,
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
                                             ContactInfoSphericParticle* const element,
                                             Condition* const wall,
                                             const double equiv_radius,
                                             const double equiv_young,
                                             double indentation,
                                             double previous_indentation,
                                             double& AuxElasticShearForce,
                                             double& MaximumAdmisibleShearForce);

        void CalculateViscoDampingForce(double LocalRelVel[3],
                                        double ViscoDampingLocalContactForce[3],
                                        ContactInfoSphericParticle* const element1,
                                        ContactInfoSphericParticle* const element2);

        void CalculateViscoDampingForceWithFEM(double LocalRelVel[3],
                                               double ViscoDampingLocalContactForce[3],
                                               ContactInfoSphericParticle* const element,
                                               Condition* const wall);

        using DEM_D_Hertz_viscous_Coulomb::CalculateNormalForce;

        using DEM_D_Hertz_viscous_Coulomb::CalculateElasticEnergyDEM;
        using DEM_D_Hertz_viscous_Coulomb::CalculateInelasticFrictionalEnergyDEM;
        using DEM_D_Hertz_viscous_Coulomb::CalculateInelasticViscodampingEnergyDEM;
        using DEM_D_Hertz_viscous_Coulomb::CalculateElasticEnergyFEM;
        using DEM_D_Hertz_viscous_Coulomb::CalculateInelasticFrictionalEnergyFEM;
        using DEM_D_Hertz_viscous_Coulomb::CalculateInelasticViscodampingEnergyFEM;

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
        }

        virtual void load(Serializer& rSerializer)  override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
        }

    }; //class DEM_D_Conical_damage

} /* namespace Kratos.*/
#endif /* DEM_D_CONICAL_DAMAGE_CL_H_INCLUDED  defined */
