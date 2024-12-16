// Authors: J. Irazábal (CIMNE)
// Date: April 2019

#include "DEM_D_Stress_dependent_cohesive_CL.h"

namespace Kratos {

    DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Stress_Dependent_Cohesive::Clone() const {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Stress_Dependent_Cohesive(*this));
        return p_clone;
    }

    std::unique_ptr<DEMDiscontinuumConstitutiveLaw> DEM_D_Stress_Dependent_Cohesive::CloneUnique() {
        return Kratos::make_unique<DEM_D_Stress_Dependent_Cohesive>();
    }

    void DEM_D_Stress_Dependent_Cohesive::Check(Properties::Pointer pProp) const {
        DEM_D_Linear_viscous_Coulomb::Check(pProp);
        if(!pProp->Has(INITIAL_COHESION)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable INITIAL_COHESION should be present in the properties when using DEM_D_Stress_Dependent_Cohesive. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(INITIAL_COHESION) = 0.0;
        }
        if(!pProp->Has(AMOUNT_OF_COHESION_FROM_STRESS)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable AMOUNT_OF_COHESION_FROM_STRESS should be present in the properties when using DEM_D_Stress_Dependent_Cohesive. 1.0e20 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(AMOUNT_OF_COHESION_FROM_STRESS) = 1.0e20;
        }
    }

    std::string DEM_D_Stress_Dependent_Cohesive::GetTypeOfLaw() {
        std::string type_of_law = "Stress_dependent";
        return type_of_law;
    }

    ///////////////////////////////////////
    //        DEM-DEM INTERACTION        //
    ///////////////////////////////////////

    void DEM_D_Stress_Dependent_Cohesive::CalculateIndentedContactArea(const double radius,
                                                                       const double other_radius,
                                                                       const double indentation,
                                                                       double &calculation_area) {

        KRATOS_TRY

        //Get equivalent Radius
        const double radius_sum      = radius + other_radius;
        const double equiv_radius    = 2.0 * (radius * other_radius) / radius_sum;
        const double normalize_dist = radius_sum / (radius_sum - indentation);
        calculation_area = Globals::Pi * equiv_radius * equiv_radius * normalize_dist;

        KRATOS_CATCH("")
    }

    void DEM_D_Stress_Dependent_Cohesive::InitializeContact(SphericParticle* const element1,
                                                            SphericParticle* const element2,
                                                            const double indentation) {

        KRATOS_TRY

        // Particles Radius Sum
        const double radius = element1->GetRadius();
        const double other_radius = element2->GetRadius();
        const double radius_sum = radius + other_radius;

        // Get equivalent Young's Modulus
        const double my_young = element1->GetYoung();
        const double other_young = element2->GetYoung();
        const double equiv_young = my_young * other_young / (other_young + my_young);

        // Get equivalent Poisson
        const double my_poisson = element1->GetPoisson();
        const double other_poisson = element2->GetPoisson();
        const double equiv_poisson = 2.0 * my_poisson * other_poisson / (my_poisson + other_poisson);

        double contact_area = 0.0;
        CalculateIndentedContactArea(radius, other_radius, indentation, contact_area);

        // Normal and Tangent elastic constants
        // Taken from 'A methodology for calibrating parameters in discrete element models based on machine learning surrogates.'
        // https://link.springer.com/article/10.1007/s40571-022-00550-1
        mKn = contact_area * equiv_young / (radius_sum - indentation);
        mKt = mKn * (2.0 * (1.0 - equiv_poisson) / (2.0 - equiv_poisson));

        KRATOS_CATCH("")
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

        KRATOS_TRY

        InitializeContact(element1, element2, indentation);

        LocalElasticContactForce[2]  = CalculateNormalForce(indentation);

        CalculateViscoDampingForce(LocalRelVel, ViscoDampingLocalContactForce, element1, element2);

        double normal_contact_force = LocalElasticContactForce[2] + ViscoDampingLocalContactForce[2];

        if (normal_contact_force < 0.0) {
            normal_contact_force = 0.0;
            ViscoDampingLocalContactForce[2] = -1.0 * LocalElasticContactForce[2];
        }

        bool initial_time_step = false;

        if (r_process_info[TIME_STEPS] == 0) initial_time_step = true;

        cohesive_force = CalculateStressDependentCohesiveNormalForce(element1, element2, normal_contact_force, indentation, initial_time_step);

        double AuxElasticShearForce;
        double MaximumAdmisibleShearForce;

        CalculateTangentialForceWithNeighbour(normal_contact_force, OldLocalElasticContactForce, LocalElasticContactForce, ViscoDampingLocalContactForce, LocalDeltDisp,
                                              LocalRelVel, sliding, element1, element2, indentation, previous_indentation, AuxElasticShearForce, MaximumAdmisibleShearForce);

        double& elastic_energy = element1->GetElasticEnergy();
        DEM_D_Linear_viscous_Coulomb::CalculateElasticEnergyDEM(elastic_energy, indentation, LocalElasticContactForce);

        if (AuxElasticShearForce > MaximumAdmisibleShearForce && MaximumAdmisibleShearForce != 0.0){
            double& inelastic_frictional_energy = element1->GetInelasticFrictionalEnergy();
            DEM_D_Linear_viscous_Coulomb::CalculateInelasticFrictionalEnergyDEM(inelastic_frictional_energy, AuxElasticShearForce, LocalElasticContactForce);
        }

        double& inelastic_viscodamping_energy = element1->GetInelasticViscodampingEnergy();
        DEM_D_Linear_viscous_Coulomb::CalculateInelasticViscodampingEnergyDEM(inelastic_viscodamping_energy, ViscoDampingLocalContactForce, LocalDeltDisp);

        KRATOS_CATCH("")
    }

    double DEM_D_Stress_Dependent_Cohesive::CalculateStressDependentCohesiveNormalForce(SphericParticle* const element1,
                                                                                        SphericParticle* const element2,
                                                                                        const double normal_contact_force,
                                                                                        const double indentation,
                                                                                        const bool initial_time_step) {

        KRATOS_TRY

        ContactInfoSphericParticle* p_element1 = dynamic_cast<ContactInfoSphericParticle*>(element1);
        ContactInfoSphericParticle* p_element2 = dynamic_cast<ContactInfoSphericParticle*>(element2);

        const double radius = p_element1->GetRadius();
        const double other_radius = p_element2->GetRadius();

        double contact_area = 0.0;
        CalculateIndentedContactArea(radius, other_radius, indentation, contact_area);

        double cohesion = 0.0;
        Properties& properties_of_this_contact = element1->GetProperties().GetSubProperties(element2->GetProperties().Id());
        double amount_of_cohesion_from_stress = properties_of_this_contact[AMOUNT_OF_COHESION_FROM_STRESS];

        for (unsigned int i = 0; p_element1->mNeighbourElements.size(); i++) {
            if (p_element1->mNeighbourElements[i]->Id() == p_element2->Id()) {
                if (initial_time_step) p_element1->mNeighbourCohesion[i] = properties_of_this_contact[INITIAL_COHESION];
                cohesion = std::min(properties_of_this_contact[PARTICLE_COHESION], amount_of_cohesion_from_stress * p_element1->mNeighbourContactStress[i]);
                if (p_element1->mNeighbourCohesion[i] != 0.0) cohesion = std::max(p_element1->mNeighbourCohesion[i], cohesion);

                double contact_stress = normal_contact_force / contact_area;
                p_element1->mNeighbourContactStress[i] = std::max(p_element1->mNeighbourContactStress[i], contact_stress);
                break;
            }
        }

        const double cohesive_force = cohesion * contact_area;

        return cohesive_force;

        KRATOS_CATCH("")
    }

    ///////////////////////////////////////
    //        DEM-FEM INTERACTION        //
    ///////////////////////////////////////

    void DEM_D_Stress_Dependent_Cohesive::CalculateIndentedContactAreaWithFEM(const double radius,
                                                                              const double indentation,
                                                                              double &calculation_area,
                                                                              const double ini_delta) {

        KRATOS_TRY

        //Get effective Radius
        const double effective_radius    = radius - ini_delta;
        const double normalize_dist = radius / (radius - indentation);
        calculation_area = Globals::Pi * effective_radius * effective_radius * normalize_dist;

        KRATOS_CATCH("")
    }

    void DEM_D_Stress_Dependent_Cohesive::InitializeContactWithFEM(SphericParticle* const element,
                                                                   Condition* const wall,
                                                                   const double indentation,
                                                                   const double ini_delta) {

        KRATOS_TRY

        //Particle Radius
        const double radius = element->GetRadius();

        //Get equivalent Young's Modulus
        const double my_young = element->GetYoung();
        const double walls_young = wall->GetProperties()[YOUNG_MODULUS];
        const double equiv_young = my_young * walls_young / (walls_young + my_young);

        //Get equivalent Poisson
        const double my_poisson = element->GetPoisson();
        const double walls_poisson = wall->GetProperties()[POISSON_RATIO];
        const double equiv_poisson = 0.5 * (my_poisson + walls_poisson);

        double contact_area = 0.0;
        CalculateIndentedContactAreaWithFEM(radius, indentation, contact_area);

        // Normal and Tangent elastic constants
        // Taken from 'A methodology for calibrating parameters in discrete element models based on machine learning surrogates.'
        // https://link.springer.com/article/10.1007/s40571-022-00550-1
        mKn = contact_area * equiv_young / (radius - indentation);
        mKt = mKn * (2.0 * (1.0 - equiv_poisson) / (2.0 - equiv_poisson));

        KRATOS_CATCH("")
    }

    void DEM_D_Stress_Dependent_Cohesive::CalculateForcesWithFEM(const ProcessInfo& r_process_info,
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

        KRATOS_TRY

        InitializeContactWithFEM(element, wall, indentation);

        LocalElasticContactForce[2] = CalculateNormalForce(indentation);

        CalculateViscoDampingForceWithFEM(LocalRelVel, ViscoDampingLocalContactForce, element, wall);

        double normal_contact_force = LocalElasticContactForce[2] + ViscoDampingLocalContactForce[2];

        if (normal_contact_force < 0.0) {
            normal_contact_force = 0.0;
            ViscoDampingLocalContactForce[2] = -1.0 * LocalElasticContactForce[2];
        }

        bool initial_time_step = false;

        if (r_process_info[TIME_STEPS] == 0) initial_time_step = true;

        cohesive_force = CalculateStressDependentCohesiveNormalForceWithFEM(element, wall, normal_contact_force, indentation, initial_time_step);

        double AuxElasticShearForce;
        double MaximumAdmisibleShearForce;

        CalculateTangentialForceWithNeighbour(normal_contact_force, OldLocalElasticContactForce, LocalElasticContactForce, ViscoDampingLocalContactForce, LocalDeltDisp,
                                              LocalRelVel, sliding, element, wall, indentation, previous_indentation, AuxElasticShearForce, MaximumAdmisibleShearForce);

        double& elastic_energy = element->GetElasticEnergy();
        DEM_D_Linear_viscous_Coulomb::CalculateElasticEnergyFEM(elastic_energy, indentation, LocalElasticContactForce);

        if (AuxElasticShearForce > MaximumAdmisibleShearForce && MaximumAdmisibleShearForce != 0.0){
            double& inelastic_frictional_energy = element->GetInelasticFrictionalEnergy();
            DEM_D_Linear_viscous_Coulomb::CalculateInelasticFrictionalEnergyFEM(inelastic_frictional_energy, AuxElasticShearForce, LocalElasticContactForce);
        }

        double& inelastic_viscodamping_energy = element->GetInelasticViscodampingEnergy();
        DEM_D_Linear_viscous_Coulomb::CalculateInelasticViscodampingEnergyFEM(inelastic_viscodamping_energy, ViscoDampingLocalContactForce, LocalDeltDisp);

        KRATOS_CATCH("")
    }

    double DEM_D_Stress_Dependent_Cohesive::CalculateStressDependentCohesiveNormalForceWithFEM(SphericParticle* const element,
                                                                                               Condition* const wall,
                                                                                               const double normal_contact_force,
                                                                                               const double indentation,
                                                                                               const bool initial_time_step) {

        KRATOS_TRY

        ContactInfoSphericParticle* p_element = dynamic_cast<ContactInfoSphericParticle*>(element);

        const double radius           = element->GetRadius();

        double contact_area = 0.0;
        CalculateIndentedContactAreaWithFEM(radius, indentation, contact_area);

        double equiv_cohesion = 0.0;
        Properties& properties_of_this_contact = element->GetProperties().GetSubProperties(wall->GetProperties().Id());
        double equiv_amount_of_cohesion_from_stress = properties_of_this_contact[AMOUNT_OF_COHESION_FROM_STRESS];

        for (unsigned int i = 0; p_element->mNeighbourRigidFaces.size(); i++) {
            if (p_element->mNeighbourRigidFaces[i]->Id() == wall->Id()) {
                if (initial_time_step) p_element->mNeighbourRigidCohesion[i] = properties_of_this_contact[INITIAL_COHESION];
                equiv_cohesion = std::min(properties_of_this_contact[PARTICLE_COHESION], equiv_amount_of_cohesion_from_stress * p_element->mNeighbourRigidContactStress[i]);
                if (p_element->mNeighbourRigidCohesion[i] != 0.0) equiv_cohesion = std::max(p_element->mNeighbourRigidCohesion[i], equiv_cohesion);

                double contact_stress = normal_contact_force / contact_area;
                p_element->mNeighbourRigidContactStress[i] = std::max(p_element->mNeighbourRigidContactStress[i], contact_stress);
                break;
            }
        }

        const double cohesive_force   = equiv_cohesion * contact_area;

        return cohesive_force;

        KRATOS_CATCH("")
    }

} // namespace Kratos
