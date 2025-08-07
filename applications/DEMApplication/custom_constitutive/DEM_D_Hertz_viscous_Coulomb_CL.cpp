// Authors: M.A. Celigueta and S. Latorre (CIMNE)
// Date: October 2015

#include "DEM_D_Hertz_viscous_Coulomb_CL.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Hertz_viscous_Coulomb::Clone() const {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Hertz_viscous_Coulomb(*this));
        return p_clone;
    }

    std::unique_ptr<DEMDiscontinuumConstitutiveLaw> DEM_D_Hertz_viscous_Coulomb::CloneUnique() {
        return Kratos::make_unique<DEM_D_Hertz_viscous_Coulomb>();
    }

    std::string DEM_D_Hertz_viscous_Coulomb::GetTypeOfLaw() {
        std::string type_of_law = "Hertz";
        return type_of_law;
    }

    void DEM_D_Hertz_viscous_Coulomb::Check(Properties::Pointer pProp) const {
        if(!pProp->Has(STATIC_FRICTION)) pProp->GetValue(STATIC_FRICTION) = 0.0;
        if(!pProp->Has(DYNAMIC_FRICTION)) pProp->GetValue(DYNAMIC_FRICTION) = 0.0;
        if(!pProp->Has(FRICTION_DECAY)) pProp->GetValue(FRICTION_DECAY) = 0.0;
        if(!pProp->Has(COEFFICIENT_OF_RESTITUTION)) pProp->GetValue(COEFFICIENT_OF_RESTITUTION) = 0.0;
    }

    void DEM_D_Hertz_viscous_Coulomb::CalculateForces(const ProcessInfo& r_process_info, const double OldLocalElasticContactForce[3], double LocalElasticContactForce[3], double LocalDeltDisp[3], double LocalRelVel[3], double indentation, double previous_indentation, double ViscoDampingLocalContactForce[3], double& cohesive_force, SphericParticle* element1, SphericParticle* element2, bool& sliding, double LocalCoordSystem[3][3]) {
        InitializeContact(element1, element2, indentation);

        LocalElasticContactForce[2] = 2.0 * mKn * indentation / 3.0;
        CalculateViscoDampingForce(LocalRelVel, ViscoDampingLocalContactForce, element1, element2);

        double normal_contact_force = LocalElasticContactForce[2] + ViscoDampingLocalContactForce[2];
        if (normal_contact_force < 0.0) {
            normal_contact_force = 0.0;
            ViscoDampingLocalContactForce[2] = -1.0 * LocalElasticContactForce[2];
        }

        double AuxElasticShearForce, MaximumAdmisibleShearForce;
        CalculateTangentialForceWithNeighbour(normal_contact_force, OldLocalElasticContactForce, LocalElasticContactForce, ViscoDampingLocalContactForce, LocalDeltDisp, LocalRelVel, sliding, element1, element2, indentation, previous_indentation, AuxElasticShearForce, MaximumAdmisibleShearForce);

        double& elastic_energy = element1->GetElasticEnergy();
        CalculateElasticEnergyDEM(elastic_energy, indentation, LocalElasticContactForce);

        if(sliding && MaximumAdmisibleShearForce != 0.0){
            double& inelastic_frictional_energy = element1->GetInelasticFrictionalEnergy();
            CalculateInelasticFrictionalEnergyDEM(inelastic_frictional_energy, AuxElasticShearForce, LocalElasticContactForce);
        }

        double& inelastic_viscodamping_energy = element1->GetInelasticViscodampingEnergy();
        CalculateInelasticViscodampingEnergyDEM(inelastic_viscodamping_energy, ViscoDampingLocalContactForce, LocalDeltDisp);

    }

    void DEM_D_Hertz_viscous_Coulomb::CalculateForcesWithFEM(const ProcessInfo& r_process_info, const double OldLocalElasticContactForce[3], double LocalElasticContactForce[3], double LocalDeltDisp[3], double LocalRelVel[3], double indentation, double previous_indentation, double ViscoDampingLocalContactForce[3], double& cohesive_force, SphericParticle* const element, Condition* const wall, bool& sliding) {
        InitializeContactWithFEM(element, wall, indentation);

        LocalElasticContactForce[2] = 2.0 * mKn * indentation / 3.0;
        CalculateViscoDampingForceWithFEM(LocalRelVel, ViscoDampingLocalContactForce, element, wall);

        double normal_contact_force = LocalElasticContactForce[2] + ViscoDampingLocalContactForce[2];
        if (normal_contact_force < 0.0) {
            normal_contact_force = 0.0;
            ViscoDampingLocalContactForce[2] = -1.0 * LocalElasticContactForce[2];
        }

        double AuxElasticShearForce, MaximumAdmisibleShearForce;
        CalculateTangentialForceWithNeighbour(normal_contact_force, OldLocalElasticContactForce, LocalElasticContactForce, ViscoDampingLocalContactForce, LocalDeltDisp, LocalRelVel, sliding, element, wall, indentation, previous_indentation, AuxElasticShearForce, MaximumAdmisibleShearForce);

        double& elastic_energy = element->GetElasticEnergy();
        CalculateElasticEnergyFEM(elastic_energy, indentation, LocalElasticContactForce);

        if(sliding && MaximumAdmisibleShearForce != 0.0){
            double& inelastic_frictional_energy = element->GetInelasticFrictionalEnergy();
            CalculateInelasticFrictionalEnergyFEM(inelastic_frictional_energy, AuxElasticShearForce, LocalElasticContactForce);
        }

        double& inelastic_viscodamping_energy = element->GetInelasticViscodampingEnergy();
        CalculateInelasticViscodampingEnergyFEM(inelastic_viscodamping_energy, ViscoDampingLocalContactForce, LocalDeltDisp);
    }

    void DEM_D_Hertz_viscous_Coulomb::InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double indentation) {
        // Equivalent radius
        const double r1 = element1->GetRadius();
        const double r2 = element2->GetRadius();
        const double r12 = r1 * r2 / (r1 + r2);

        // Equivalent Young's modulus
        const double E1 = element1->GetYoung();
        const double E2 = element2->GetYoung();
        const double v1 = element1->GetPoisson();
        const double v2 = element2->GetPoisson();
        const double E12 = E1 * E2 / (E2 * (1.0 - v1 * v1) + E1 * (1.0 - v2 * v2));

        // Equivalent Shear modulus
        const double G1 = 0.5 * E1 / (1.0 + v1);
        const double G2 = 0.5 * E2 / (1.0 + v2);
        const double G12 = 1.0 / ((2.0 - v1) / G1 + (2.0 - v2) / G2);

        // Normal and tangent elastic constants
        mKn = 2.0 * E12 * sqrt(r12 * indentation);
        mKt = 4.0 * G12 * mKn / E12;
    }

    void DEM_D_Hertz_viscous_Coulomb::InitializeContactWithFEM(SphericParticle* const element, Condition* const wall, const double indentation, const double ini_delta) {
        // Equivalent radius
        const double r12 = element->GetRadius();

        // Equivalent Young's modulus
        const double E1 = element->GetYoung();
        const double E2 = wall->GetProperties()[YOUNG_MODULUS];
        const double v1 = element->GetPoisson();
        const double v2 = wall->GetProperties()[POISSON_RATIO];
        const double E12 = E1 * E2 / (E2 * (1.0 - v1 * v1) + E1 * (1.0 - v2 * v2));

        // Equivalent Shear modulus
        const double G1 = 0.5 * E1 / (1.0 + v1);
        const double G2 = 0.5 * E2 / (1.0 + v2);
        const double G12 = 1.0 / ((2.0 - v1) / G1 + (2.0 - v2) / G2);

        //Normal and Tangent elastic constants
        mKn = 2.0 * E12 * sqrt(r12 * indentation);
        mKt = 4.0 * G12 * mKn / E12;
    }

    void DEM_D_Hertz_viscous_Coulomb::CalculateViscoDampingForce(double LocalRelVel[3], double ViscoDampingLocalContactForce[3], SphericParticle* const element1, SphericParticle* const element2) {
        const double m1 = element1->GetMass();
        const double m2 = element2->GetMass();
        const double equiv_mass = 1.0 / (1.0 / m1 + 1.0 / m2);

        Properties& properties_of_this_contact = element1->GetProperties().GetSubProperties(element2->GetProperties().Id());
        const double damping_gamma = properties_of_this_contact[DAMPING_GAMMA];
        const double equiv_visco_damp_coeff_tangential = 2.0 * damping_gamma * sqrt(equiv_mass * mKt);
        const double equiv_visco_damp_coeff_normal     = 2.0 * damping_gamma * sqrt(equiv_mass * mKn);

        ViscoDampingLocalContactForce[0] = -equiv_visco_damp_coeff_tangential * LocalRelVel[0];
        ViscoDampingLocalContactForce[1] = -equiv_visco_damp_coeff_tangential * LocalRelVel[1];
        ViscoDampingLocalContactForce[2] = -equiv_visco_damp_coeff_normal     * LocalRelVel[2];
    }

    void DEM_D_Hertz_viscous_Coulomb::CalculateViscoDampingForceWithFEM(double LocalRelVel[3], double ViscoDampingLocalContactForce[3], SphericParticle* const element, Condition* const wall) {
        const double m1 = element->GetMass();

        Properties& properties_of_this_contact = element->GetProperties().GetSubProperties(wall->GetProperties().Id());
        const double damping_gamma = properties_of_this_contact[DAMPING_GAMMA];
        const double tangential_damping_coefficient = 2.0 * damping_gamma * sqrt(m1 * mKt);
        const double normal_damping_coefficient     = 2.0 * damping_gamma * sqrt(m1 * mKn);

        ViscoDampingLocalContactForce[0] = -tangential_damping_coefficient * LocalRelVel[0];
        ViscoDampingLocalContactForce[1] = -tangential_damping_coefficient * LocalRelVel[1];
        ViscoDampingLocalContactForce[2] = -normal_damping_coefficient     * LocalRelVel[2];
    }

    template<class NeighbourClassType>
    void DEM_D_Hertz_viscous_Coulomb::CalculateTangentialForceWithNeighbour(const double normal_contact_force, const double OldLocalElasticContactForce[3], double LocalElasticContactForce[3], double ViscoDampingLocalContactForce[3], const double LocalDeltDisp[3], const double LocalRelVel[3], bool& sliding, SphericParticle* const element, NeighbourClassType* const neighbour, double indentation, double previous_indentation, double& AuxElasticShearForce, double& MaximumAdmisibleShearForce) {
        if (previous_indentation > indentation) {
            const double minoring_factor = sqrt(indentation / previous_indentation);
            LocalElasticContactForce[0] = OldLocalElasticContactForce[0] * minoring_factor - mKt * LocalDeltDisp[0];
            LocalElasticContactForce[1] = OldLocalElasticContactForce[1] * minoring_factor - mKt * LocalDeltDisp[1];
        }
        else {
            LocalElasticContactForce[0] = OldLocalElasticContactForce[0] - mKt * LocalDeltDisp[0];
            LocalElasticContactForce[1] = OldLocalElasticContactForce[1] - mKt * LocalDeltDisp[1];
        }

        const double tangential_contact_force_0 = LocalElasticContactForce[0] + ViscoDampingLocalContactForce[0];
        const double tangential_contact_force_1 = LocalElasticContactForce[1] + ViscoDampingLocalContactForce[1];

        AuxElasticShearForce = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);
        const double ActualTotalShearForce = sqrt(tangential_contact_force_0 * tangential_contact_force_0 + tangential_contact_force_1 * tangential_contact_force_1);

        Properties& properties_of_this_contact = element->GetProperties().GetSubProperties(neighbour->GetProperties().Id());
        double friction = properties_of_this_contact[DYNAMIC_FRICTION];
        MaximumAdmisibleShearForce = normal_contact_force * friction;

        if (ActualTotalShearForce > MaximumAdmisibleShearForce) {
            const double ActualElasticShearForce = AuxElasticShearForce;
            const double dot_product = LocalElasticContactForce[0] * ViscoDampingLocalContactForce[0] + LocalElasticContactForce[1] * ViscoDampingLocalContactForce[1];
            const double ViscoDampingLocalContactForceModule = sqrt(ViscoDampingLocalContactForce[0] * ViscoDampingLocalContactForce[0] + ViscoDampingLocalContactForce[1] * ViscoDampingLocalContactForce[1]);

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

    // Energy calculations
    void DEM_D_Hertz_viscous_Coulomb::CalculateElasticEnergyDEM(double& elastic_energy, double indentation, double LocalElasticContactForce[3]) {
        double normal_elastic     = 0.20*LocalElasticContactForce[2]*indentation; //each ball in a contact with another ball receives half the contact energy
        double tangential_elastic = 0.25*(LocalElasticContactForce[0]*LocalElasticContactForce[0]+LocalElasticContactForce[1]*LocalElasticContactForce[1])/mKt; //each ball in a contact with another ball receives half the contact energy
        elastic_energy += normal_elastic;
        elastic_energy += tangential_elastic;
    }

    void DEM_D_Hertz_viscous_Coulomb::CalculateInelasticFrictionalEnergyDEM(double& inelastic_frictional_energy, double& AuxElasticShearForce, double LocalElasticContactForce[3]) {
        double frictional_energy = 0.25*((AuxElasticShearForce*AuxElasticShearForce)-(LocalElasticContactForce[0]*LocalElasticContactForce[0]+LocalElasticContactForce[1]*LocalElasticContactForce[1]))/mKt;
        inelastic_frictional_energy += frictional_energy;
    }

    void DEM_D_Hertz_viscous_Coulomb::CalculateInelasticViscodampingEnergyDEM(double& inelastic_viscodamping_energy, double ViscoDampingLocalContactForce[3], double LocalDeltDisp[3]) {
        double viscodamping_energy = 0.50 * sqrt(ViscoDampingLocalContactForce[0]*ViscoDampingLocalContactForce[0]*LocalDeltDisp[0]*LocalDeltDisp[0]+ViscoDampingLocalContactForce[1]*ViscoDampingLocalContactForce[1]*LocalDeltDisp[1]*LocalDeltDisp[1]+ViscoDampingLocalContactForce[2]*ViscoDampingLocalContactForce[2]*LocalDeltDisp[2]*LocalDeltDisp[2]);
        inelastic_viscodamping_energy += viscodamping_energy;
    }

    void DEM_D_Hertz_viscous_Coulomb::CalculateElasticEnergyFEM(double& elastic_energy, double indentation, double LocalElasticContactForce[3]) {
        double normal_elastic     = 0.40*LocalElasticContactForce[2]*indentation; //each ball in a contact with a wall receives all the contact energy
        double tangential_elastic = 0.50*(LocalElasticContactForce[0]*LocalElasticContactForce[0]+LocalElasticContactForce[1]*LocalElasticContactForce[1])/mKt;
        elastic_energy += normal_elastic;
        elastic_energy += tangential_elastic;
    }

    void DEM_D_Hertz_viscous_Coulomb::CalculateInelasticFrictionalEnergyFEM(double& inelastic_frictional_energy, double& AuxElasticShearForce, double LocalElasticContactForce[3]) {
        double frictional_energy = 0.50 * ((AuxElasticShearForce*AuxElasticShearForce)-(LocalElasticContactForce[0]*LocalElasticContactForce[0]+LocalElasticContactForce[1]*LocalElasticContactForce[1]))/mKt;
        inelastic_frictional_energy += frictional_energy;
    }

    void DEM_D_Hertz_viscous_Coulomb::CalculateInelasticViscodampingEnergyFEM(double& inelastic_viscodamping_energy, double ViscoDampingLocalContactForce[3], double LocalDeltDisp[3]) {
        double viscodamping_energy = sqrt(ViscoDampingLocalContactForce[0]*ViscoDampingLocalContactForce[0]*LocalDeltDisp[0]*LocalDeltDisp[0]+ViscoDampingLocalContactForce[1]*ViscoDampingLocalContactForce[1]*LocalDeltDisp[1]*LocalDeltDisp[1]+ViscoDampingLocalContactForce[2]*ViscoDampingLocalContactForce[2]*LocalDeltDisp[2]*LocalDeltDisp[2]);
        inelastic_viscodamping_energy += viscodamping_energy;
    }

    double DEM_D_Hertz_viscous_Coulomb::GetTangentialStiffness() {
        return mKt;
    }
}
