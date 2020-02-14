// Authors: J. Irazábal (CIMNE)
// Date: November 2016

#include "DEM_D_Conical_damage_CL.h"
#include "DEM_D_Hertz_viscous_Coulomb_CL.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Conical_damage::Clone() const {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Conical_damage(*this));
        return p_clone;
    }

    void DEM_D_Conical_damage::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        if(verbose) KRATOS_INFO("DEM") << "Assigning DEM_D_Conical_damage to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
        this->Check(pProp);
    }

    void DEM_D_Conical_damage::Check(Properties::Pointer pProp) const {
        DEMDiscontinuumConstitutiveLaw::Check(pProp);
        if(!pProp->Has(CONICAL_DAMAGE_CONTACT_RADIUS)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable CONICAL_DAMAGE_CONTACT_RADIUS should be present in the properties when using DEM_D_Conical_damage. 90.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(CONICAL_DAMAGE_CONTACT_RADIUS) = 0.0;
        }
        if(!pProp->Has(CONICAL_DAMAGE_MAX_STRESS)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable CONICAL_DAMAGE_MAX_STRESS should be present in the properties when using DEM_D_Conical_damage. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(CONICAL_DAMAGE_MAX_STRESS) = 1.0e20;
        }
        if(!pProp->Has(CONICAL_DAMAGE_ALPHA)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable CONICAL_DAMAGE_ALPHA should be present in the properties when using DEM_D_Conical_damage. 90.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(CONICAL_DAMAGE_ALPHA) = 90.0;
        }
        if(!pProp->Has(CONICAL_DAMAGE_GAMMA)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable CONICAL_DAMAGE_GAMMA should be present in the properties when using DEM_D_Conical_damage. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(CONICAL_DAMAGE_GAMMA) = 0.0;
        }
    }

    std::string DEM_D_Conical_damage::GetTypeOfLaw() {
        std::string type_of_law = "Conical_damge";
        return type_of_law;
    }

    /////////////////////////
    // DEM-DEM INTERACTION //
    /////////////////////////

    void DEM_D_Conical_damage::InitializeDependentContact(double equiv_radius,
                                                          const double equiv_level_of_fouling,
                                                          const double equiv_young,
                                                          const double equiv_shear,
                                                          const double indentation) {
        //Normal and Tangent elastic constants
        const double sqrt_equiv_radius_and_indentation = sqrt(equiv_level_of_fouling * equiv_radius * indentation);
        mKn = 2.0 * equiv_young * sqrt_equiv_radius_and_indentation;
        mKt = 4.0 * equiv_shear * mKn / equiv_young;
    }

    void DEM_D_Conical_damage::DamageContact(ContactInfoSphericParticle* const element1,
                                             ContactInfoSphericParticle* const element2,
                                             double& equiv_radius,
                                             const double equiv_level_of_fouling,
                                             const double equiv_young,
                                             const double equiv_shear,
                                             double& indentation,
                                             const double normal_contact_force) {

        //Get new Equivalent Radius
        double equiv_radius_new = (equiv_young * sqrt(6 * normal_contact_force)) / (pow(Globals::Pi * element1->GetParticleConicalDamageMaxStress(),1.5));

        double offset = 0.0;
        if (equiv_radius_new > equiv_level_of_fouling * equiv_radius) {
            const double AlphaFunction = element1->GetProperties()[CONICAL_DAMAGE_ALPHA_FUNCTION];
            offset = (equiv_radius_new - equiv_radius) * AlphaFunction;

            equiv_radius = equiv_radius_new;
        }

        for (unsigned int i = 0; element1->mNeighbourElements.size(); i++) {
            if (element1->mNeighbourElements[i]->Id() == element2->Id()) {
                element1->mNeighbourContactRadius[i] = equiv_radius;
                if (indentation > offset) element1->mNeighbourIndentation[i] = indentation - offset;
                else element1->mNeighbourIndentation[i] = 0.0;
                indentation = element1->mNeighbourIndentation[i];
                break;
            }
        }

        //New Normal and Tangent elastic constants
        const double sqrt_equiv_radius_and_indentation = sqrt(equiv_level_of_fouling * equiv_radius * indentation);
        mKn = 2.0 * equiv_young * sqrt_equiv_radius_and_indentation;
        mKt = 4.0 * equiv_shear * mKn / equiv_young;
    }

    void DEM_D_Conical_damage::CalculateForces(const ProcessInfo& r_process_info,
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
                                               double LocalCoordSystem[3][3]) {

        ContactInfoSphericParticle* p_element1 = dynamic_cast<ContactInfoSphericParticle*>(element1);
        ContactInfoSphericParticle* p_element2 = dynamic_cast<ContactInfoSphericParticle*>(element2);

        //Get equivalent Radius
        const double my_radius      = p_element1->GetParticleConicalDamageContactRadius();
        const double other_radius   = p_element2->GetParticleConicalDamageContactRadius();
        const double radius_sum     = my_radius + other_radius;
        const double radius_sum_inv = 1.0 / radius_sum;
        double equiv_radius         = my_radius * other_radius * radius_sum_inv;

        double elastic_indentation = indentation;

        for (unsigned int i = 0; p_element1->mNeighbourElements.size(); i++) {
            if (p_element1->mNeighbourElements[i]->Id() == p_element2->Id()) {
                if (p_element1->mNeighbourContactRadius[i] > equiv_radius) {
                    equiv_radius = p_element1->mNeighbourContactRadius[i];
                    elastic_indentation = p_element1->mNeighbourIndentation[i] + (indentation - previous_indentation);
                    if (elastic_indentation < 0.0) p_element1->mNeighbourIndentation[i] = elastic_indentation = 0.0;
                }
                break;
            }
        }

        if (elastic_indentation > 0.0) {
            //Get equivalent Young's Modulus
            const double my_young        = p_element1->GetYoung();
            const double other_young     = p_element2->GetYoung();
            const double my_poisson      = p_element1->GetPoisson();
            const double other_poisson   = p_element2->GetPoisson();
            const double equiv_young     = my_young * other_young / (other_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - other_poisson * other_poisson));

            //Get equivalent Shear Modulus
            const double my_shear_modulus = 0.5 * my_young / (1.0 + my_poisson);
            const double other_shear_modulus = 0.5 * other_young / (1.0 + other_poisson);
            const double equiv_shear = 1.0 / ((2.0 - my_poisson)/my_shear_modulus + (2.0 - other_poisson)/other_shear_modulus);

            //Level of fouling in case it is considered
            const double equiv_level_of_fouling = 0.5 * ((1.0 + p_element1->GetLevelOfFouling()) + (1.0 + p_element2->GetLevelOfFouling()));

            InitializeDependentContact(equiv_radius, equiv_level_of_fouling, equiv_young, equiv_shear, elastic_indentation);

            LocalElasticContactForce[2] = DEM_D_Hertz_viscous_Coulomb::CalculateNormalForce(elastic_indentation);

            CalculateViscoDampingForce(LocalRelVel, ViscoDampingLocalContactForce, p_element1, p_element2);

            double normal_contact_force = LocalElasticContactForce[2] + ViscoDampingLocalContactForce[2];

            if (normal_contact_force < 0.0) {
                normal_contact_force = 0.0;
                ViscoDampingLocalContactForce[2] = -1.0 * LocalElasticContactForce[2];
            }

            double contact_stress = (3 * normal_contact_force) / (2 * Globals::Pi * equiv_level_of_fouling * equiv_radius * elastic_indentation);

            if (contact_stress > p_element1->GetParticleConicalDamageMaxStress()) {
                DamageContact(p_element1, p_element2, equiv_radius, equiv_level_of_fouling, equiv_young, equiv_shear, elastic_indentation, normal_contact_force);
                LocalElasticContactForce[2] = DEM_D_Hertz_viscous_Coulomb::CalculateNormalForce(elastic_indentation);
            }

            double AuxElasticShearForce;
            double MaximumAdmisibleShearForce;

            CalculateTangentialForce(normal_contact_force, OldLocalElasticContactForce, LocalElasticContactForce, ViscoDampingLocalContactForce, LocalDeltDisp,
                                     sliding, p_element1, p_element2, equiv_radius, equiv_young, elastic_indentation, previous_indentation, AuxElasticShearForce, MaximumAdmisibleShearForce);

            double& elastic_energy = p_element1->GetElasticEnergy();
            DEM_D_Hertz_viscous_Coulomb::CalculateElasticEnergyDEM(elastic_energy, elastic_indentation, LocalElasticContactForce);

            if(sliding && MaximumAdmisibleShearForce != 0.0){
                double& inelastic_frictional_energy = p_element1->GetInelasticFrictionalEnergy();
                DEM_D_Hertz_viscous_Coulomb::CalculateInelasticFrictionalEnergyDEM(inelastic_frictional_energy, AuxElasticShearForce, LocalElasticContactForce);
            }

            double& inelastic_viscodamping_energy = p_element1->GetInelasticViscodampingEnergy();
            DEM_D_Hertz_viscous_Coulomb::CalculateInelasticViscodampingEnergyDEM(inelastic_viscodamping_energy, ViscoDampingLocalContactForce, LocalDeltDisp);
        }
    }

    /////////////////////////
    // DEM-FEM INTERACTION //
    /////////////////////////

    void DEM_D_Conical_damage::InitializeDependentContactWithFEM(double effective_radius,
                                                                 const double equiv_level_of_fouling,
                                                                 const double equiv_young,
                                                                 const double equiv_shear,
                                                                 const double indentation) {
        //Normal and Tangent elastic constants
        const double sqrt_equiv_radius_and_indentation = sqrt(equiv_level_of_fouling * effective_radius * indentation);
        mKn = 2.0 * equiv_young * sqrt_equiv_radius_and_indentation;
        mKt = 4.0 * equiv_shear * mKn / equiv_young;
    }

    void DEM_D_Conical_damage::DamageContactWithFEM(ContactInfoSphericParticle* const element,
                                                    Condition* const wall,
                                                    double& effective_radius,
                                                    const double equiv_level_of_fouling,
                                                    const double equiv_young,
                                                    const double equiv_shear,
                                                    double& indentation,
                                                    const double normal_contact_force) {
        //Get new Equivalent Radius
        double effective_radius_new = (equiv_young * sqrt(6 * normal_contact_force)) / (pow(Globals::Pi * element->GetParticleConicalDamageMaxStress(),1.5));

        double offset = 0.0;
        if (effective_radius_new > equiv_level_of_fouling * effective_radius) {
            const double AlphaFunction = element->GetProperties()[CONICAL_DAMAGE_ALPHA_FUNCTION];
            offset = (effective_radius_new - effective_radius) * AlphaFunction;

            effective_radius = effective_radius_new;
        }

        for (unsigned int i = 0; element->mNeighbourRigidFaces.size(); i++) {
            if (element->mNeighbourRigidFaces[i]->Id() == wall->Id()) {
                element->mNeighbourRigidContactRadius[i] = effective_radius;
                if (indentation > offset) element->mNeighbourRigidIndentation[i] = indentation - offset;
                else element->mNeighbourRigidIndentation[i] = 0.0;
                element->mNeighbourRigidIndentation[i] = indentation;
                break;
            }
        }

        //New Normal and Tangent elastic constants
        const double sqrt_equiv_radius_and_indentation = sqrt(equiv_level_of_fouling * effective_radius * indentation);
        mKn = 2.0 * equiv_young * sqrt_equiv_radius_and_indentation;
        mKt = 4.0 * equiv_shear * mKn / equiv_young;
    }

    void DEM_D_Conical_damage::CalculateForcesWithFEM(ProcessInfo& r_process_info,
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
                                                      bool& sliding) {

        ContactInfoSphericParticle* p_element = dynamic_cast<ContactInfoSphericParticle*>(element);

        //Get effective Radius
        double effective_radius    = p_element->GetParticleConicalDamageContactRadius();

        double elastic_indentation = indentation;

        for (unsigned int i = 0; p_element->mNeighbourRigidFaces.size(); i++) {
            if (p_element->mNeighbourRigidFaces[i]->Id() == wall->Id()) {
                if (p_element->mNeighbourRigidContactRadius[i] > effective_radius) {
                    effective_radius = p_element->mNeighbourRigidContactRadius[i];
                    elastic_indentation = p_element->mNeighbourRigidIndentation[i] + (indentation - previous_indentation);
                    if (elastic_indentation < 0.0) p_element->mNeighbourRigidIndentation[i] = elastic_indentation = 0.0;
                }
                break;
            }
        }

        if (elastic_indentation > 0.0) {
            //Get equivalent Young's Modulus
            const double my_young            = p_element->GetYoung();
            const double walls_young         = wall->GetProperties()[YOUNG_MODULUS];
            const double my_poisson          = p_element->GetPoisson();
            const double walls_poisson       = wall->GetProperties()[POISSON_RATIO];
            const double equiv_young         = my_young * walls_young / (walls_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - walls_poisson * walls_poisson));

            //Get equivalent Shear Modulus
            const double my_shear_modulus    = 0.5 * my_young / (1.0 + my_poisson);
            const double walls_shear_modulus = 0.5 * walls_young / (1.0 + walls_poisson);
            const double equiv_shear         = 1.0 / ((2.0 - my_poisson)/my_shear_modulus + (2.0 - walls_poisson)/walls_shear_modulus);


            //Level of fouling in case it is considered
            const double equiv_level_of_fouling = 1.0 + p_element->GetLevelOfFouling();

            InitializeDependentContactWithFEM(effective_radius, equiv_level_of_fouling, equiv_young, equiv_shear, elastic_indentation);

            LocalElasticContactForce[2] = DEM_D_Hertz_viscous_Coulomb::CalculateNormalForce(elastic_indentation);

            CalculateViscoDampingForceWithFEM(LocalRelVel, ViscoDampingLocalContactForce, p_element, wall);

            double normal_contact_force = LocalElasticContactForce[2] + ViscoDampingLocalContactForce[2];

            if (normal_contact_force < 0.0) {
                normal_contact_force = 0.0;
                ViscoDampingLocalContactForce[2] = -1.0 * LocalElasticContactForce[2];
            }

            double contact_stress = (3 * normal_contact_force) / (2 * Globals::Pi * equiv_level_of_fouling * effective_radius * elastic_indentation);

            if (contact_stress > p_element->GetParticleConicalDamageMaxStress()) {
                DamageContactWithFEM(p_element, wall, effective_radius, equiv_level_of_fouling, equiv_young, equiv_shear, elastic_indentation, normal_contact_force);
                LocalElasticContactForce[2] = DEM_D_Hertz_viscous_Coulomb::CalculateNormalForce(elastic_indentation);
            }

            double AuxElasticShearForce;
            double MaximumAdmisibleShearForce;

            CalculateTangentialForceWithFEM(normal_contact_force, OldLocalElasticContactForce, LocalElasticContactForce, ViscoDampingLocalContactForce, LocalDeltDisp,
                                            sliding, p_element, wall, effective_radius, equiv_young, elastic_indentation, previous_indentation, AuxElasticShearForce, MaximumAdmisibleShearForce);

            double& elastic_energy = p_element->GetElasticEnergy();
            DEM_D_Hertz_viscous_Coulomb::CalculateElasticEnergyFEM(elastic_energy, elastic_indentation, LocalElasticContactForce);

            if(sliding && MaximumAdmisibleShearForce != 0.0){
                double& inelastic_frictional_energy = p_element->GetInelasticFrictionalEnergy();
                DEM_D_Hertz_viscous_Coulomb::CalculateInelasticFrictionalEnergyFEM(inelastic_frictional_energy, AuxElasticShearForce, LocalElasticContactForce);
            }

            double& inelastic_viscodamping_energy = p_element->GetInelasticViscodampingEnergy();
            DEM_D_Hertz_viscous_Coulomb::CalculateInelasticViscodampingEnergyFEM(inelastic_viscodamping_energy, ViscoDampingLocalContactForce, LocalDeltDisp);
        }
    }

    void DEM_D_Conical_damage::CalculateViscoDampingForce(double LocalRelVel[3],
                                                          double ViscoDampingLocalContactForce[3],
                                                          ContactInfoSphericParticle* const element1,
                                                          ContactInfoSphericParticle* const element2) {

        const double my_mass    = element1->GetMass();
        const double other_mass = element2->GetMass();

        const double equiv_mass = 1.0 / (1.0/my_mass + 1.0/other_mass);

        const double my_gamma    = element1->GetProperties()[DAMPING_GAMMA];
        const double other_gamma = element2->GetProperties()[DAMPING_GAMMA];
        const double equiv_gamma = 0.5 * (my_gamma + other_gamma);

        const double equiv_visco_damp_coeff_normal     = 2.0 * equiv_gamma * sqrt(equiv_mass * mKn);
        const double equiv_visco_damp_coeff_tangential = 2.0 * equiv_gamma * sqrt(equiv_mass * mKt);

        ViscoDampingLocalContactForce[0] = - equiv_visco_damp_coeff_tangential * LocalRelVel[0];
        ViscoDampingLocalContactForce[1] = - equiv_visco_damp_coeff_tangential * LocalRelVel[1];
        ViscoDampingLocalContactForce[2] = - equiv_visco_damp_coeff_normal     * LocalRelVel[2];
    }

    void DEM_D_Conical_damage::CalculateTangentialForce(const double normal_contact_force,
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
                                                        double& MaximumAdmisibleShearForce) {

        LocalElasticContactForce[0] = OldLocalElasticContactForce[0] - mKt * LocalDeltDisp[0];
        LocalElasticContactForce[1] = OldLocalElasticContactForce[1] - mKt * LocalDeltDisp[1];

        if (previous_indentation > indentation) {
            const double minoring_factor = sqrt (indentation / previous_indentation);
            LocalElasticContactForce[0] = OldLocalElasticContactForce[0] * minoring_factor - mKt * LocalDeltDisp[0];
            LocalElasticContactForce[1] = OldLocalElasticContactForce[1] * minoring_factor - mKt * LocalDeltDisp[1];
        }

        AuxElasticShearForce = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);

        const double my_tg_of_friction_angle        = element1->GetTgOfFrictionAngle();
        const double neighbour_tg_of_friction_angle = element2->GetTgOfFrictionAngle();
        double equiv_tg_of_fri_ang                  = 0.5 * (my_tg_of_friction_angle + neighbour_tg_of_friction_angle);

        if(equiv_tg_of_fri_ang < 0.0) {
            KRATOS_ERROR << "The averaged friction is negative for one contact of element with Id: "<< element1->Id()<<std::endl;
        }

        if (fabs(equiv_tg_of_fri_ang) > 1.0e-12) {
            double critical_force = 0.6666666666666667 * Globals::Pi * equiv_radius * indentation * element1->GetParticleConicalDamageMaxStress();
            double critical_force_inv = 1.0  / critical_force;
            equiv_tg_of_fri_ang *= pow((normal_contact_force * critical_force_inv), element1->GetParticleConicalDamageGamma());
        }

        for (unsigned int i = 0; element1->mNeighbourElements.size(); i++) {
            if (element1->mNeighbourElements[i]->Id() == element2->Id()) {
                if (element1->mNeighbourTgOfFriAng[i] <= equiv_tg_of_fri_ang) equiv_tg_of_fri_ang = element1->mNeighbourTgOfFriAng[i];
                else element1->mNeighbourTgOfFriAng[i] = equiv_tg_of_fri_ang;
                break;
            }
        }

        MaximumAdmisibleShearForce = normal_contact_force * equiv_tg_of_fri_ang;

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

    void DEM_D_Conical_damage::CalculateViscoDampingForceWithFEM(double LocalRelVel[3],
                                                                 double ViscoDampingLocalContactForce[3],
                                                                 ContactInfoSphericParticle* const element,
                                                                 Condition* const wall) {

        const double my_mass    = element->GetMass();
        const double gamma = element->GetProperties()[DAMPING_GAMMA];
        const double normal_damping_coefficient     = 2.0 * gamma * sqrt(my_mass * mKn);
        const double tangential_damping_coefficient = 2.0 * gamma * sqrt(my_mass * mKt);

        ViscoDampingLocalContactForce[0] = - tangential_damping_coefficient * LocalRelVel[0];
        ViscoDampingLocalContactForce[1] = - tangential_damping_coefficient * LocalRelVel[1];
        ViscoDampingLocalContactForce[2] = - normal_damping_coefficient     * LocalRelVel[2];

    }

    void DEM_D_Conical_damage::CalculateTangentialForceWithFEM(const double normal_contact_force,
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
                                                               double& MaximumAdmisibleShearForce) {

        LocalElasticContactForce[0] = OldLocalElasticContactForce[0] - mKt * LocalDeltDisp[0];
        LocalElasticContactForce[1] = OldLocalElasticContactForce[1] - mKt * LocalDeltDisp[1];

        if (previous_indentation > indentation) {
            const double minoring_factor = sqrt (indentation / previous_indentation);
            LocalElasticContactForce[0] = OldLocalElasticContactForce[0] * minoring_factor - mKt * LocalDeltDisp[0];
            LocalElasticContactForce[1] = OldLocalElasticContactForce[1] * minoring_factor - mKt * LocalDeltDisp[1];
        }

        AuxElasticShearForce = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);

        const double my_tg_of_friction_angle   = element->GetTgOfFrictionAngle();
        const double wall_tg_of_friction_angle = wall->GetProperties()[FRICTION];
        double equiv_tg_of_fri_ang             = 0.5 * (my_tg_of_friction_angle + wall_tg_of_friction_angle);

        if (fabs(equiv_tg_of_fri_ang) > 1.0e-12) {
            double critical_force = 0.6666666666666667 * Globals::Pi * equiv_radius * indentation * element->GetParticleConicalDamageMaxStress();
            double critical_force_inv = 1.0  / critical_force;
            equiv_tg_of_fri_ang *= pow((normal_contact_force * critical_force_inv), element->GetParticleConicalDamageGamma());
        }

        for (unsigned int i = 0; element->mNeighbourRigidFaces.size(); i++) {
            if (element->mNeighbourRigidFaces[i]->Id() == wall->Id()) {
                if (element->mNeighbourRigidTgOfFriAng[i] <= equiv_tg_of_fri_ang) equiv_tg_of_fri_ang = element->mNeighbourRigidTgOfFriAng[i];
                else element->mNeighbourRigidTgOfFriAng[i] = equiv_tg_of_fri_ang;
                break;
            }
        }

        MaximumAdmisibleShearForce = normal_contact_force * equiv_tg_of_fri_ang;

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

} // namespace Kratos
