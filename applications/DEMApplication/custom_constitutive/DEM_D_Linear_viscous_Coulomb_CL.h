//Authors: M.A. Celigueta and S. Latorre (CIMNE)
//   Date: July 2015

#if !defined(DEM_D_LINEAR_VISCOUS_COULOMB_CL_H_INCLUDED)
#define DEM_D_LINEAR_VISCOUS_COULOMB_CL_H_INCLUDED

#include <string>
#include <iostream>
#include "DEM_discontinuum_constitutive_law.h"

namespace Kratos {

    class SphericParticle;

    class KRATOS_API(DEM_APPLICATION) DEM_D_Linear_viscous_Coulomb : public DEMDiscontinuumConstitutiveLaw {

    public:

        using DEMDiscontinuumConstitutiveLaw::CalculateNormalForce;

        KRATOS_CLASS_POINTER_DEFINITION(DEM_D_Linear_viscous_Coulomb);

        DEM_D_Linear_viscous_Coulomb() {}

        ~DEM_D_Linear_viscous_Coulomb() {}

        void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true) override;

        void Check(Properties::Pointer pProp) const override;

        std::string GetTypeOfLaw() override;

        DEMDiscontinuumConstitutiveLaw::Pointer Clone() const override;

        void InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) override;

        void InitializeContactWithFEM(SphericParticle* const element, Condition* const wall, const double indentation, const double ini_delta = 0.0) override;

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
                            bool& sliding, double LocalCoordSystem[3][3]) override;

        void CalculateForcesWithFEM(const ProcessInfo& r_process_info,
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
                                            const double indentation) override;

        double CalculateCohesiveNormalForceWithFEM(SphericParticle* const element,
                                                   Condition* const wall,
                                                   const double indentation) override;

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
                                                   double& MaximumAdmisibleShearForce){

            LocalElasticContactForce[0] = OldLocalElasticContactForce[0] - mKt * LocalDeltDisp[0];
            LocalElasticContactForce[1] = OldLocalElasticContactForce[1] - mKt * LocalDeltDisp[1];

            AuxElasticShearForce = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);

            const double my_tg_of_static_friction_angle        = GetTgOfStaticFrictionAngleOfElement(element);
            const double neighbour_tg_of_static_friction_angle = neighbour->GetProperties()[STATIC_FRICTION];
            const double equiv_tg_of_static_fri_ang            = 0.5 * (my_tg_of_static_friction_angle + neighbour_tg_of_static_friction_angle);

            const double my_tg_of_dynamic_friction_angle        = GetTgOfDynamicFrictionAngleOfElement(element);
            const double neighbour_tg_of_dynamic_friction_angle = neighbour->GetProperties()[DYNAMIC_FRICTION];
            const double equiv_tg_of_dynamic_fri_ang            = 0.5 * (my_tg_of_dynamic_friction_angle + neighbour_tg_of_dynamic_friction_angle);

            if(equiv_tg_of_static_fri_ang < 0.0 || equiv_tg_of_dynamic_fri_ang < 0.0) {
                KRATOS_ERROR << "The averaged friction is negative for one contact of element with Id: "<< GetElementId(element)<<std::endl;
            }

            MaximumAdmisibleShearForce = normal_contact_force * equiv_tg_of_static_fri_ang;
            if (AuxElasticShearForce > MaximumAdmisibleShearForce) MaximumAdmisibleShearForce = normal_contact_force * equiv_tg_of_dynamic_fri_ang;

            const double tangential_contact_force_0 = LocalElasticContactForce[0] + ViscoDampingLocalContactForce[0];
            const double tangential_contact_force_1 = LocalElasticContactForce[1] + ViscoDampingLocalContactForce[1];

            const double ActualTotalShearForce = sqrt(tangential_contact_force_0 * tangential_contact_force_0 + tangential_contact_force_1 * tangential_contact_force_1);

            if (ActualTotalShearForce > MaximumAdmisibleShearForce) {

                const double ActualElasticShearForce = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);

                const double dot_product = LocalElasticContactForce[0] * ViscoDampingLocalContactForce[0] + LocalElasticContactForce[1] * ViscoDampingLocalContactForce[1];
                const double ViscoDampingLocalContactForceModule = sqrt(ViscoDampingLocalContactForce[0] * ViscoDampingLocalContactForce[0] +\
                                                                        ViscoDampingLocalContactForce[1] * ViscoDampingLocalContactForce[1]);

                if (dot_product >= 0.0) {

                    if (ActualElasticShearForce > MaximumAdmisibleShearForce) {
                        const double fraction = MaximumAdmisibleShearForce / ActualElasticShearForce;
                        LocalElasticContactForce[0]      = LocalElasticContactForce[0] * fraction;
                        LocalElasticContactForce[1]      = LocalElasticContactForce[1] * fraction;
                        ViscoDampingLocalContactForce[0] = 0.0;
                        ViscoDampingLocalContactForce[1] = 0.0;
                    }
                    else {
                        const double ActualViscousShearForce = MaximumAdmisibleShearForce - ActualElasticShearForce;
                        const double fraction = ActualViscousShearForce / ViscoDampingLocalContactForceModule;
                        ViscoDampingLocalContactForce[0] *= fraction;
                        ViscoDampingLocalContactForce[1] *= fraction;
                    }
                }
                else {
                    if (ViscoDampingLocalContactForceModule >= ActualElasticShearForce) {
                        const double fraction = (MaximumAdmisibleShearForce + ActualElasticShearForce) / ViscoDampingLocalContactForceModule;
                        ViscoDampingLocalContactForce[0] *= fraction;
                        ViscoDampingLocalContactForce[1] *= fraction;
                    }
                    else {
                        const double fraction = MaximumAdmisibleShearForce / ActualElasticShearForce;
                        LocalElasticContactForce[0]      = LocalElasticContactForce[0] * fraction;
                        LocalElasticContactForce[1]      = LocalElasticContactForce[1] * fraction;
                        ViscoDampingLocalContactForce[0] = 0.0;
                        ViscoDampingLocalContactForce[1] = 0.0;
                    }
                }
                sliding = true;
            }
        }

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

    protected:
        double GetTgOfDynamicFrictionAngleOfElement(SphericParticle* element);
        double GetTgOfStaticFrictionAngleOfElement(SphericParticle* element);
        std::size_t GetElementId(SphericParticle* element);

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMDiscontinuumConstitutiveLaw)
                    //rSerializer.load("MyMemberName",myMember);
        }

    }; //class DEM_D_Linear_viscous_Coulomb

} /* namespace Kratos.*/
#endif /* DEM_D_LINEAR_VISCOUS_COULOMB_CL_H_INCLUDED  defined */
