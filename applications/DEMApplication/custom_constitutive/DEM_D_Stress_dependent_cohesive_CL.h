//Author: J. Irazábal (CIMNE)
//   Date: March 2019

#if !defined(DEM_D_STRESS_DEPENDENT_COHESIVE_LAW_CL_H_INCLUDED)
#define DEM_D_STRESS_DEPENDENT_COHESIVE_LAW_CL_H_INCLUDED

#include <string>
#include <iostream>
#include "DEM_discontinuum_constitutive_law.h"

namespace Kratos {

    class SphericParticle;

    class KRATOS_API(DEM_APPLICATION) DEM_D_Stress_Dependent_Cohesive : public DEMDiscontinuumConstitutiveLaw {

    public:

        using DEMDiscontinuumConstitutiveLaw::CalculateNormalForce;

        KRATOS_CLASS_POINTER_DEFINITION(DEM_D_Stress_Dependent_Cohesive);

        DEM_D_Stress_Dependent_Cohesive() {}

        ~DEM_D_Stress_Dependent_Cohesive() {}

        void Initialize(const ProcessInfo& r_process_info) override;

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) const override;

        void Check(Properties::Pointer pProp) const override;

        std::string GetTypeOfLaw() override;

        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const override;

        void InitializeContact(SphericParticle* const element1,
                               SphericParticle* const element2,
                               const double indentation) override;

        void InitializeContactWithFEM(SphericParticle* const element,
                                      Condition* const wall,
                                      const double indentation,
                                      const double ini_delta = 0.0) override;

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

        double CalculateCohesiveNormalForce(SphericParticle* const element1,
                                            SphericParticle* const element2,
                                            const double normal_contact_force,
                                            const double indentation);

        double CalculateCohesiveNormalForceWithFEM(SphericParticle* const element,
                                                   Condition* const wall,
                                                   const double normal_contact_force,
                                                   const double indentation);

        template <class NeighbourClassType>
        void CalculateTangentialForceWithNeighbour(const double normal_contact_force,
                                                   const double OldLocalElasticContactForce[3],
                                                   double LocalElasticContactForce[3],
                                                   double ViscoDampingLocalContactForce[3],
                                                   const double LocalDeltDisp[3],
                                                   bool& sliding,
                                                   SphericParticle* const element,
                                                   NeighbourClassType* const neighbour,
                                                   double indentation,
                                                   double previous_indentation,
                                                   double& AuxElasticShearForce,
                                                   double& MaximumAdmisibleShearForce);

        void CalculateViscoDampingForce(double LocalRelVel[3],
                                        double ViscoDampingLocalContactForce[3],
                                        SphericParticle* const element1,
                                        SphericParticle* const element2);

        void CalculateViscoDampingForceWithFEM(double LocalRelVel[3],
                                               double ViscoDampingLocalContactForce[3],
                                               SphericParticle* const element,
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

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
        }

    }; //class DEM_D_Stress_Dependent_Cohesive

} /* namespace Kratos.*/
#endif /* DEM_D_STRESS_DEPENDENT_COHESIVE_LAW_H defined */
