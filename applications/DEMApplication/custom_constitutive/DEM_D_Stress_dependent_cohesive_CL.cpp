// Authors: J. Irazábal (CIMNE)
// Date: April 2019

#include "DEM_D_Stress_dependent_cohesive_CL.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Stress_Dependent_Cohesive::Clone() const {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Stress_Dependent_Cohesive(*this));
        return p_clone;
    }

    void DEM_D_Stress_Dependent_Cohesive::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        if(verbose) KRATOS_INFO("DEM") << "Assigning DEM_D_Stress_Dependent_Cohesive to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
        this->Check(pProp);
    }

    void DEM_D_Stress_Dependent_Cohesive::Check(Properties::Pointer pProp) const {
        DEMDiscontinuumConstitutiveLaw::Check(pProp);
        if(!pProp->Has(PARTICLE_INITIAL_COHESION)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable PARTICLE_INITIAL_COHESION should be present in the properties when using DEM_D_Stress_Dependent_Cohesive. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(PARTICLE_INITIAL_COHESION) = 0.0;
        }
        if(!pProp->Has(AMOUNT_OF_COHESION_FROM_STRESS)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable AMOUNT_OF_COHESION_FROM_STRESS should be present in the properties when using DEM_D_Stress_Dependent_Cohesive. 1.0e20 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(AMOUNT_OF_COHESION_FROM_STRESS) = 1.0e20;
        }
    }

    std::string DEM_D_Stress_Dependent_Cohesive::GetTypeOfLaw() {
        std::string type_of_law = "Stress_Dependent_Cohesive";
        return type_of_law;
    }

    /////////////////////////
    // DEM-DEM INTERACTION //
    /////////////////////////

    void DEM_D_Stress_Dependent_Cohesive::InitializeContact(SphericParticle* const element1,
                                                            SphericParticle* const element2,
                                                            double& contact_area,
                                                            const double indentation) {

        //Get equivalent Radius
        const double my_radius       = element1->GetRadius();
        const double other_radius    = element2->GetRadius();
        const double radius_sum      = my_radius + other_radius;
        const double equiv_radius    = 2.0 * (my_radius * other_radius) / radius_sum;

        //Get equivalent Young's Modulus
        const double my_young        = element1->GetYoung();
        const double other_young     = element2->GetYoung();
        const double my_poisson      = element1->GetPoisson();
        const double other_poisson   = element2->GetPoisson();
        const double equiv_young     = my_young * other_young / (other_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - other_poisson * other_poisson));

        //Get equivalent Shear Modulus
        const double my_shear_modulus = 0.5 * my_young / (1.0 + my_poisson);
        const double other_shear_modulus = 0.5 * other_young / (1.0 + other_poisson);
        const double equiv_shear = 1.0 / ((2.0 - my_poisson)/my_shear_modulus + (2.0 - other_poisson)/other_shear_modulus);

        const double normalize_dist = radius_sum / (radius_sum - indentation);
        contact_area = Globals::Pi * equiv_radius * equiv_radius * normalize_dist;

        //Normal and Tangent elastic constants
        mKn = contact_area * equiv_young / (radius_sum - indentation);
        mKt = 4.0 * equiv_shear * mKn / equiv_young;
    }

    void DEM_D_Stress_Dependent_Cohesive::CalculateForces(const ProcessInfo& r_process_info,
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

        double contact_area;

        InitializeContact(element1, element2, contact_area, indentation);

        LocalElasticContactForce[2]  = CalculateNormalForce(indentation);

        CalculateViscoDampingForce(LocalRelVel, ViscoDampingLocalContactForce, element1, element2);

        double normal_contact_force = LocalElasticContactForce[2] + ViscoDampingLocalContactForce[2];

        if (normal_contact_force < 0.0) {
            normal_contact_force = 0.0;
            ViscoDampingLocalContactForce[2] = -1.0 * LocalElasticContactForce[2];
        }

        bool initial_time_step = false;

        if (r_process_info[TIME_STEPS] == 0) initial_time_step = true;

        cohesive_force = CalculateCohesiveNormalForce(element1, element2, normal_contact_force, contact_area, indentation, initial_time_step);

        double AuxElasticShearForce;
        double MaximumAdmisibleShearForce;

        CalculateTangentialForceWithNeighbour(normal_contact_force, OldLocalElasticContactForce, LocalElasticContactForce, ViscoDampingLocalContactForce, LocalDeltDisp,
                                              sliding, element1, element2, contact_area, indentation, previous_indentation, AuxElasticShearForce, MaximumAdmisibleShearForce);

        double& elastic_energy = element1->GetElasticEnergy();
        DEM_D_Linear_viscous_Coulomb::CalculateElasticEnergyDEM(elastic_energy, indentation, LocalElasticContactForce);

        if (AuxElasticShearForce > MaximumAdmisibleShearForce && MaximumAdmisibleShearForce != 0.0){
            double& inelastic_frictional_energy = element1->GetInelasticFrictionalEnergy();
            DEM_D_Linear_viscous_Coulomb::CalculateInelasticFrictionalEnergyDEM(inelastic_frictional_energy, AuxElasticShearForce, LocalElasticContactForce);
        }

        double& inelastic_viscodamping_energy = element1->GetInelasticViscodampingEnergy();
        DEM_D_Linear_viscous_Coulomb::CalculateInelasticViscodampingEnergyDEM(inelastic_viscodamping_energy, ViscoDampingLocalContactForce, LocalDeltDisp);

    }

    /////////////////////////
    // DEM-FEM INTERACTION //
    /////////////////////////

    void DEM_D_Stress_Dependent_Cohesive::InitializeContactWithFEM(SphericParticle* const element,
                                                                   Condition* const wall,
                                                                   double& contact_area,
                                                                   const double indentation,
                                                                   const double ini_delta) {

        //Get effective Radius
        const double my_radius           = element->GetRadius(); //Get equivalent Radius
        const double effective_radius    = my_radius - ini_delta;

        //Get equivalent Young's Modulus
        const double my_young            = element->GetYoung();
        const double walls_young         = wall->GetProperties()[YOUNG_MODULUS];
        const double my_poisson          = element->GetPoisson();
        const double walls_poisson       = wall->GetProperties()[POISSON_RATIO];
        const double equiv_young         = my_young * walls_young / (walls_young * (1.0 - my_poisson * my_poisson) + my_young * (1.0 - walls_poisson * walls_poisson));

        //Get equivalent Shear Modulus
        const double my_shear_modulus    = 0.5 * my_young / (1.0 + my_poisson);
        const double walls_shear_modulus = 0.5 * walls_young / (1.0 + walls_poisson);
        const double equiv_shear         = 1.0 / ((2.0 - my_poisson)/my_shear_modulus + (2.0 - walls_poisson)/walls_shear_modulus);

        const double normalize_dist = my_radius / (my_radius - indentation);
        contact_area = Globals::Pi * effective_radius * effective_radius * normalize_dist;

        mKn = contact_area * equiv_young / (my_radius - indentation);
        mKt = 4.0 * equiv_shear * mKn / equiv_young;
    }

    void DEM_D_Stress_Dependent_Cohesive::CalculateForcesWithFEM(ProcessInfo& r_process_info,
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

        double contact_area;

        InitializeContactWithFEM(element, wall, contact_area, indentation);

        LocalElasticContactForce[2] = CalculateNormalForce(indentation);

        CalculateViscoDampingForceWithFEM(LocalRelVel, ViscoDampingLocalContactForce, element, wall);

        double normal_contact_force = LocalElasticContactForce[2] + ViscoDampingLocalContactForce[2];

        if (normal_contact_force < 0.0) {
            normal_contact_force = 0.0;
            ViscoDampingLocalContactForce[2] = -1.0 * LocalElasticContactForce[2];
        }

        bool initial_time_step = false;

        if (r_process_info[TIME_STEPS] == 0) initial_time_step = true;

        cohesive_force = CalculateCohesiveNormalForceWithFEM(element, wall, normal_contact_force, contact_area, indentation, initial_time_step);

        double AuxElasticShearForce;
        double MaximumAdmisibleShearForce;

        CalculateTangentialForceWithNeighbour(normal_contact_force, OldLocalElasticContactForce, LocalElasticContactForce, ViscoDampingLocalContactForce, LocalDeltDisp,
                                              sliding, element, wall, contact_area, indentation, previous_indentation, AuxElasticShearForce, MaximumAdmisibleShearForce);

        double& elastic_energy = element->GetElasticEnergy();
        DEM_D_Linear_viscous_Coulomb::CalculateElasticEnergyFEM(elastic_energy, indentation, LocalElasticContactForce);

        if (AuxElasticShearForce > MaximumAdmisibleShearForce && MaximumAdmisibleShearForce != 0.0){
            double& inelastic_frictional_energy = element->GetInelasticFrictionalEnergy();
            DEM_D_Linear_viscous_Coulomb::CalculateInelasticFrictionalEnergyFEM(inelastic_frictional_energy, AuxElasticShearForce, LocalElasticContactForce);
        }

        double& inelastic_viscodamping_energy = element->GetInelasticViscodampingEnergy();
        DEM_D_Linear_viscous_Coulomb::CalculateInelasticViscodampingEnergyFEM(inelastic_viscodamping_energy, ViscoDampingLocalContactForce, LocalDeltDisp);
    }


    double DEM_D_Stress_Dependent_Cohesive::CalculateNormalForce(const double indentation) {
        return mKn * indentation;
    }

    double DEM_D_Stress_Dependent_Cohesive::CalculateCohesiveNormalForce(SphericParticle* const element1,
                                                                         SphericParticle* const element2,
                                                                         const double normal_contact_force,
                                                                         const double contact_area,
                                                                         const double indentation,
                                                                         const bool initial_time_step) {

        ContactInfoSphericParticle* p_element1 = dynamic_cast<ContactInfoSphericParticle*>(element1);
        ContactInfoSphericParticle* p_element2 = dynamic_cast<ContactInfoSphericParticle*>(element2);

        double equiv_cohesion = 0.0;
        double equiv_amount_of_cohesion_from_stress = 0.5 * (p_element1->GetAmountOfCohesionFromStress() + p_element2->GetAmountOfCohesionFromStress());

        for (unsigned int i = 0; p_element1->mNeighbourElements.size(); i++) {
            if (p_element1->mNeighbourElements[i]->Id() == p_element2->Id()) {
                if (initial_time_step) p_element1->mNeighbourCohesion[i] = 0.5 * (p_element1->GetParticleInitialCohesion() + p_element2->GetParticleInitialCohesion());
                equiv_cohesion = std::min(0.5 * (p_element1->GetParticleCohesion() + p_element2->GetParticleCohesion()), equiv_amount_of_cohesion_from_stress * p_element1->mNeighbourContactStress[i]);
                if (p_element1->mNeighbourCohesion[i] != 0.0) equiv_cohesion = std::max(p_element1->mNeighbourCohesion[i], equiv_cohesion);

                double contact_stress = normal_contact_force / contact_area;
                p_element1->mNeighbourContactStress[i] = std::max(p_element1->mNeighbourContactStress[i], contact_stress);
                break;
            }
        }

        const double cohesive_force = equiv_cohesion * contact_area;

        return cohesive_force;
    }

    double DEM_D_Stress_Dependent_Cohesive::CalculateCohesiveNormalForceWithFEM(SphericParticle* const element,
                                                                                Condition* const wall,
                                                                                const double normal_contact_force,
                                                                                const double contact_area,
                                                                                const double indentation,
                                                                                const bool initial_time_step) {

        ContactInfoSphericParticle* p_element = dynamic_cast<ContactInfoSphericParticle*>(element);

        double equiv_cohesion = 0.0;
        double equiv_amount_of_cohesion_from_stress = p_element->GetAmountOfCohesionFromStress();

        for (unsigned int i = 0; p_element->mNeighbourRigidFaces.size(); i++) {
            if (p_element->mNeighbourRigidFaces[i]->Id() == wall->Id()) {
                if (initial_time_step) p_element->mNeighbourRigidCohesion[i] = 0.5 * (p_element->GetParticleInitialCohesion() + wall->GetProperties()[WALL_INITIAL_COHESION]);
                equiv_cohesion = std::min(0.5 * (p_element->GetParticleCohesion() + wall->GetProperties()[WALL_COHESION]), equiv_amount_of_cohesion_from_stress * p_element->mNeighbourRigidContactStress[i]);
                if (p_element->mNeighbourRigidCohesion[i] != 0.0) equiv_cohesion = std::max(p_element->mNeighbourRigidCohesion[i], equiv_cohesion);

                double contact_stress = normal_contact_force / contact_area;
                p_element->mNeighbourRigidContactStress[i] = std::max(p_element->mNeighbourRigidContactStress[i], contact_stress);
                break;
            }
        }

        const double cohesive_force   = equiv_cohesion * contact_area;

        return cohesive_force;
    }

    template<class NeighbourClassType>

    void DEM_D_Stress_Dependent_Cohesive::CalculateTangentialForceWithNeighbour(const double normal_contact_force,
                                                                                const double OldLocalElasticContactForce[3],
                                                                                double LocalElasticContactForce[3],
                                                                                double ViscoDampingLocalContactForce[3],
                                                                                const double LocalDeltDisp[3],
                                                                                bool& sliding,
                                                                                SphericParticle* const element,
                                                                                NeighbourClassType* const neighbour,
                                                                                const double contact_area,
                                                                                double indentation,
                                                                                double previous_indentation,
                                                                                double& AuxElasticShearForce,
                                                                                double& MaximumAdmisibleShearForce) {

        ContactInfoSphericParticle* p_element = dynamic_cast<ContactInfoSphericParticle*>(element);

        LocalElasticContactForce[0] = OldLocalElasticContactForce[0] - mKt * LocalDeltDisp[0];
        LocalElasticContactForce[1] = OldLocalElasticContactForce[1] - mKt * LocalDeltDisp[1];

        AuxElasticShearForce = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);

        const double my_tg_of_static_friction_angle        = p_element->GetTgOfStaticFrictionAngle();
        const double neighbour_tg_of_static_friction_angle = neighbour->GetProperties()[STATIC_FRICTION];
        const double equiv_tg_of_static_fri_ang            = 0.5 * (my_tg_of_static_friction_angle + neighbour_tg_of_static_friction_angle);

        const double my_tg_of_dynamic_friction_angle        = p_element->GetTgOfDynamicFrictionAngle();
        const double neighbour_tg_of_dynamic_friction_angle = neighbour->GetProperties()[DYNAMIC_FRICTION];
        const double equiv_tg_of_dynamic_fri_ang            = 0.5 * (my_tg_of_dynamic_friction_angle + neighbour_tg_of_dynamic_friction_angle);

        if(equiv_tg_of_static_fri_ang < 0.0 || equiv_tg_of_dynamic_fri_ang < 0.0) {
            KRATOS_ERROR << "The averaged friction is negative for one contact of element with Id: "<< p_element->Id()<<std::endl;
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
} // namespace Kratos
