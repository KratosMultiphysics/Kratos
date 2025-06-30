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

    double DEMDiscontinuumConstitutiveLaw::GetTangentialStiffness(){
        KRATOS_ERROR << "This function (DEMContinuumConstitutiveLaw::GetTangentialStiffness) shouldn't be accessed, use derived class instead"<<std::endl;
        return 0.0;
    }

} // KRATOS
