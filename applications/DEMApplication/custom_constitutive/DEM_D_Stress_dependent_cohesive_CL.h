//Author: J. Irazábal (CIMNE)
//   Date: March 2019

#if !defined(DEM_D_STRESS_DEPENDENT_COHESIVE_LAW_CL_H_INCLUDED)
#define DEM_D_STRESS_DEPENDENT_COHESIVE_LAW_CL_H_INCLUDED

#include <string>
#include <iostream>

// Project includes
#include "DEM_D_Linear_viscous_Coulomb_CL.h"
#include "custom_elements/spheric_particle.h"
#include "custom_elements/contact_info_spheric_particle.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) DEM_D_Stress_Dependent_Cohesive : public DEM_D_Linear_viscous_Coulomb {

    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEM_D_Stress_Dependent_Cohesive);

        DEM_D_Stress_Dependent_Cohesive() {}

        ~DEM_D_Stress_Dependent_Cohesive() {}

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) override;

        void Check(Properties::Pointer pProp) const override;

        std::string GetTypeOfLaw() override;

        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const override;

        void InitializeContact(ContactInfoSphericParticle* const element1,
                               ContactInfoSphericParticle* const element2,
                               const double indentation);

        void InitializeContactWithFEM(ContactInfoSphericParticle* const element,
                                      Condition* const wall,
                                      const double indentation,
                                      const double ini_delta = 0.0);

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

        double CalculateNormalForce(const double indentation) override;

        double CalculateCohesiveNormalForce(ContactInfoSphericParticle* const element1,
                                            ContactInfoSphericParticle* const element2,
                                            const double normal_contact_force,
                                            const double indentation);

        double CalculateCohesiveNormalForceWithFEM(ContactInfoSphericParticle* const element,
                                                   Condition* const wall,
                                                   const double normal_contact_force,
                                                   const double indentation);

        void CalculateViscoDampingForce(double LocalRelVel[3],
                                        double ViscoDampingLocalContactForce[3],
                                        ContactInfoSphericParticle* const element1,
                                        ContactInfoSphericParticle* const element2);

        void CalculateViscoDampingForceWithFEM(double LocalRelVel[3],
                                               double ViscoDampingLocalContactForce[3],
                                               ContactInfoSphericParticle* const element,
                                               Condition* const wall);

        template <class NeighbourClassType>
        void CalculateTangentialForceWithNeighbour(const double normal_contact_force,
                                                   const double OldLocalElasticContactForce[3],
                                                   double LocalElasticContactForce[3],
                                                   double ViscoDampingLocalContactForce[3],
                                                   const double LocalDeltDisp[3],
                                                   bool& sliding,
                                                   ContactInfoSphericParticle* const element,
                                                   NeighbourClassType* const neighbour,
                                                   double indentation,
                                                   double previous_indentation,
                                                   double& AuxElasticShearForce,
                                                   double& MaximumAdmisibleShearForce);

        using DEM_D_Linear_viscous_Coulomb::CalculateViscoDampingForce;
        using DEM_D_Linear_viscous_Coulomb::CalculateViscoDampingForceWithFEM;

        using DEM_D_Linear_viscous_Coulomb::CalculateElasticEnergyDEM;
        using DEM_D_Linear_viscous_Coulomb::CalculateInelasticFrictionalEnergyDEM;
        using DEM_D_Linear_viscous_Coulomb::CalculateInelasticViscodampingEnergyDEM;
        using DEM_D_Linear_viscous_Coulomb::CalculateElasticEnergyFEM;
        using DEM_D_Linear_viscous_Coulomb::CalculateInelasticFrictionalEnergyFEM;
        using DEM_D_Linear_viscous_Coulomb::CalculateInelasticViscodampingEnergyFEM;

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
        }

    }; //class DEM_D_Stress_Dependent_Cohesive

} /* namespace Kratos.*/
#endif /* DEM_D_STRESS_DEPENDENT_COHESIVE_LAW_H defined */
