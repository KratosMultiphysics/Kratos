// Last modified by: S. Latorre (CIMNE)
// Date: October 2015

#include "DEM_discontinuum_constitutive_law.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

    DEMDiscontinuumConstitutiveLaw::DEMDiscontinuumConstitutiveLaw() {}

    //copy constructor
    DEMDiscontinuumConstitutiveLaw::DEMDiscontinuumConstitutiveLaw(const DEMDiscontinuumConstitutiveLaw &rReferenceDiscontinuumConstitutiveLaw) {}

    void DEMDiscontinuumConstitutiveLaw::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        if (verbose) KRATOS_INFO("DEM")  << "Assigning " << pProp->GetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_NAME) << " to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
        this->Check(pProp);
    }

    std::unique_ptr<DEMDiscontinuumConstitutiveLaw> DEMDiscontinuumConstitutiveLaw::CloneUnique() {
        KRATOS_ERROR << "This function (DEMDiscontinuumConstitutiveLaw::CloneUnique) shouldn't be accessed, use derived class instead"<<std::endl;
    }

    DEMDiscontinuumConstitutiveLaw::Pointer DEMDiscontinuumConstitutiveLaw::Clone() const {
        DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEMDiscontinuumConstitutiveLaw(*this));
        return p_clone;
    }

    DEMDiscontinuumConstitutiveLaw::~DEMDiscontinuumConstitutiveLaw() {}

    void DEMDiscontinuumConstitutiveLaw::Check(Properties::Pointer pProp) const {
        KRATOS_ERROR << "This function (DEMContinuumConstitutiveLaw::Check) shouldn't be accessed, use derived class instead"<<std::endl;
    }

    void DEMDiscontinuumConstitutiveLaw::GetContactStiffness(SphericParticle* const element1, SphericParticle* const element2, const double ini_delta, double& kn,double& kt){

      InitializeContact(element1, element2, ini_delta);
      kn = mKn;
      kt = mKt;
    }

    void DEMDiscontinuumConstitutiveLaw::InitializeContact(SphericParticle* const element1, SphericParticle* const element2, const double ini_delta) {
        KRATOS_ERROR << "This function (DEMContinuumConstitutiveLaw::InitializeContact) shouldn't be accessed, use derived class instead"<<std::endl;
    }

    void DEMDiscontinuumConstitutiveLaw::InitializeContactWithFEM(SphericParticle* const element, Condition* const wall, const double indentation, const double ini_delta) {
        KRATOS_ERROR << "This function (DEMContinuumConstitutiveLaw::InitializeContactWithFEM) shouldn't be accessed, use derived class instead"<<std::endl;
    }

    std::string DEMDiscontinuumConstitutiveLaw::GetTypeOfLaw() {
        KRATOS_ERROR << "This function (DEMContinuumConstitutiveLaw::GetTypeOfLaw) shouldn't be accessed, use derived class instead"<<std::endl;
        std::string type_of_law = "";
        return type_of_law;
    }

    void DEMDiscontinuumConstitutiveLaw::CalculateForces(const ProcessInfo& r_process_info,
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
                                                        bool& sliding, double LocalCoordSystem[3][3]) {

        KRATOS_ERROR << "This function (DEMContinuumConstitutiveLaw::CalculateForces) shouldn't be accessed, use derived class instead"<<std::endl;
    }

    void DEMDiscontinuumConstitutiveLaw::CalculateForcesWithFEM(const ProcessInfo& r_process_info,
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
                                                                bool& sliding) {

        KRATOS_ERROR << "This function (DEMContinuumConstitutiveLaw::CalculateForcesWithFEM) shouldn't be accessed, use derived class instead"<<std::endl;
    }

    double DEMDiscontinuumConstitutiveLaw::CalculateNormalForce(const double indentation) {
        KRATOS_ERROR << "This function (DEMContinuumConstitutiveLaw::CalculateNormalForce) shouldn't be accessed, use derived class instead"<<std::endl;
    }

    double DEMDiscontinuumConstitutiveLaw::CalculateNormalForce(SphericParticle* const element1, SphericParticle* const element2, const double indentation,
        double LocalCoordSystem[3][3]) {
        KRATOS_ERROR << "This function (DEMDiscontinuumConstitutiveLaw::CalculateNormalForce) shouldn't be accessed, use derived class instead"<< std::endl;
        return 0.0;
    }

    double DEMDiscontinuumConstitutiveLaw::CalculateNormalForce(SphericParticle* const element, Condition* const wall, const double indentation){
        KRATOS_ERROR << "This function (DEMDiscontinuumConstitutiveLaw::CalculateNormalForce) shouldn't be accessed, use derived class instead"<< std::endl;
        return 0.0;
    }

    double DEMDiscontinuumConstitutiveLaw::CalculateCohesiveNormalForce(SphericParticle* const element1, SphericParticle* const element2, const double indentation) {
        KRATOS_ERROR << "This function (DEMContinuumConstitutiveLaw::CalculateCohesiveNormalForce) shouldn't be accessed, use derived class instead"<<std::endl;
        return 0.0;
    }

    double DEMDiscontinuumConstitutiveLaw::CalculateCohesiveNormalForceWithFEM(SphericParticle* const element, Condition* const wall, const double indentation){
        KRATOS_ERROR << "This function (DEMContinuumConstitutiveLaw::CalculateCohesiveNormalForceWithFEM) shouldn't be accessed, use derived class instead"<<std::endl;
        return 0.0;
    }


    //TODO: This function will be deleted in the near future.
    // double DEMDiscontinuumConstitutiveLaw::LocalPeriod(const int i,
    //                                                     SphericParticle* element1,
    //                                                     SphericParticle* element2) {
    //     KRATOS_ERROR << "This function (DEMContinuumConstitutiveLaw::LocalPeriod) is deprecated."<<std::endl;
    //     double myYoung = element1->GetYoung();
    //     double other_young = element2->GetYoung();
    //     double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);
    //     const double my_radius = element1->GetRadius();
    //     const double other_radius = element2->GetRadius();
    //     double calculation_area = 0;
    //     CalculateContactArea(my_radius, other_radius, calculation_area);

    //     double radius_sum = my_radius + other_radius;
    //     const double radius_sum_inv  = 1.0 / radius_sum;
    //     const double equiv_radius    = my_radius * other_radius * radius_sum_inv;
    //     const double modified_radius = equiv_radius * 0.31225;    // sqrt(alpha * (2.0 - alpha)) = 0.31225
    //     double kn = equiv_young * Globals::Pi * modified_radius;  // 2.0 * equiv_young * sqrt_equiv_radius;

    //     const double mRealMass = element1->GetMass();
    //     const double other_real_mass = element2->GetMass();
    //     double equiv_mass = (mRealMass*other_real_mass)/(mRealMass+other_real_mass);

    //     // calculation of damping gamma
    //     Properties& properties_of_this_contact = element1->GetProperties().GetSubProperties(element2->GetProperties().Id());
    //     const double damping_gamma = properties_of_this_contact[DAMPING_GAMMA];
    //     const double friction_coeff = properties_of_this_contact[STATIC_FRICTION];
    //     const double viscous_damping_coeff     = 2.0 * damping_gamma * sqrt(equiv_mass * kn);
    //     double rescaled_damping = viscous_damping_coeff/(2*equiv_mass);
    //     double sqr_period = sqrt(1+friction_coeff*friction_coeff) * kn / equiv_mass - rescaled_damping*rescaled_damping;
    //     return sqr_period;
    // }

} // KRATOS
