#include "DEM_continuum_constitutive_law.h"
#include "custom_elements/spheric_continuum_particle.h"

namespace Kratos {

    DEMContinuumConstitutiveLaw::DEMContinuumConstitutiveLaw() {
        //KRATOS_INFO("DEM") << " DEMContinuumConstitutiveLaw constructor.." << std::endl;

    } // Class DEMContinuumConstitutiveLaw

    DEMContinuumConstitutiveLaw::DEMContinuumConstitutiveLaw(const DEMContinuumConstitutiveLaw &rReferenceContinuumConstitutiveLaw) {
        //KRATOS_INFO("DEM") << " DEMContinuumConstitutiveLaw copy constructor.." << std::endl;
    }

    std::string DEMContinuumConstitutiveLaw::GetTypeOfLaw() {
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMContinuumConstitutiveLaw::GetTypeOfLaw) should not be called.","")
        std::string type_of_law = "";
        return type_of_law;
    }

    void DEMContinuumConstitutiveLaw::Initialize(SphericContinuumParticle* owner_sphere) {
    }

    void DEMContinuumConstitutiveLaw::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        if(verbose) KRATOS_INFO("DEM") << "Assigning DEMContinuumConstitutiveLaw to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
        this->Check(pProp);
    }

    void DEMContinuumConstitutiveLaw::Check(Properties::Pointer pProp) const {
        if(!pProp->Has(FRICTION)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable FRICTION should be present in the properties when using DEMDiscontinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(FRICTION) = 0.0;
        }
        if(!pProp->Has(YOUNG_MODULUS)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable YOUNG_MODULUS should be present in the properties when using DEMDiscontinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(YOUNG_MODULUS) = 0.0;
        }
        if(!pProp->Has(POISSON_RATIO)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable POISSON_RATIO should be present in the properties when using DEMDiscontinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(POISSON_RATIO) = 0.0;
        }
        if(!pProp->Has(DAMPING_GAMMA)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable DAMPING_GAMMA should be present in the properties when using DEMDiscontinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(DAMPING_GAMMA) = 0.0;
        }
        if(!pProp->Has(COEFFICIENT_OF_RESTITUTION)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable COEFFICIENT_OF_RESTITUTION should be present in the properties when using DEMDiscontinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(COEFFICIENT_OF_RESTITUTION) = 0.0;
        }
    }

    DEMContinuumConstitutiveLaw::Pointer DEMContinuumConstitutiveLaw::Clone() const {
        DEMContinuumConstitutiveLaw::Pointer p_clone(new DEMContinuumConstitutiveLaw(*this));
        return p_clone;
    }

    DEMContinuumConstitutiveLaw::~DEMContinuumConstitutiveLaw() { //KRATOS_INFO("DEM") << "Law destructor..." ;
    }

    void DEMContinuumConstitutiveLaw::CalculateViscoDamping(double LocalRelVel[3],
                                                            double ViscoDampingLocalContactForce[3],
                                                            double indentation,
                                                            double equiv_visco_damp_coeff_normal,
                                                            double equiv_visco_damp_coeff_tangential,
                                                            bool& sliding,
                                                            int failure_id) {

        KRATOS_TRY

         if ((indentation > 0) || (failure_id == 0)) {
             ViscoDampingLocalContactForce[2] = -equiv_visco_damp_coeff_normal * LocalRelVel[2];
         }
         if (((indentation > 0) || (failure_id == 0)) && (sliding == false)) {
             ViscoDampingLocalContactForce[0] = -equiv_visco_damp_coeff_tangential * LocalRelVel[0];
             ViscoDampingLocalContactForce[1] = -equiv_visco_damp_coeff_tangential * LocalRelVel[1];
         }

        KRATOS_CATCH("")
    }

    void DEMContinuumConstitutiveLaw::ComputeParticleRotationalMoments(SphericContinuumParticle* element,
                                                                       SphericContinuumParticle* neighbor,
                                                                       double equiv_young,
                                                                       double distance,
                                                                       double calculation_area,
                                                                       double LocalCoordSystem[3][3],
                                                                       double ElasticLocalRotationalMoment[3],
                                                                       double ViscoLocalRotationalMoment[3],
                                                                       double equiv_poisson,
                                                                       double indentation) {
        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMContinuumConstitutiveLaw::ComputeParticleRotationalMoments) should not be called.","")
    }

    void DEMContinuumConstitutiveLaw::AddPoissonContribution(const double equiv_poisson, double LocalCoordSystem[3][3], double& normal_force, double calculation_area, BoundedMatrix<double, 3, 3>* mSymmStressTensor,
                                                             SphericContinuumParticle* element1, SphericContinuumParticle* element2, const ProcessInfo& r_process_info, const int i_neighbor_count, const double indentation) {
    }


    double DEMContinuumConstitutiveLaw::LocalMaxSearchDistance(const int i,
                                          SphericContinuumParticle* element1,
                                          SphericContinuumParticle* element2) {

        KRATOS_THROW_ERROR(std::runtime_error,"This function (DEMContinuumConstitutiveLaw::LocalMaxSearchDistance) should not be called.","")

    }

    double DEMContinuumConstitutiveLaw::LocalPeriod(const int i, SphericContinuumParticle* element1,
                                                                 SphericContinuumParticle* element2) {

        // calculation of equivalent young modulus
        double myYoung = element1->GetYoung();
        double other_young = element2->GetYoung();
        double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);

        const double my_radius = element1->GetRadius();
        const double other_radius = element2->GetRadius();
        double calculation_area = 0;
        CalculateContactArea(my_radius, other_radius, calculation_area);

        double radius_sum = my_radius + other_radius;
        double initial_delta = element1->GetInitialDelta(i);
        double initial_dist = radius_sum - initial_delta;

        // calculation of elastic constants
        double kn_el = equiv_young * calculation_area / initial_dist;

        const double mRealMass = element1->GetMass();  // { mRealMass = real_mass;  GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = real_mass;}
        const double &other_real_mass = element2->GetMass();

        // calculation of damping gamma
        const double my_gamma    = element1->GetProperties()[DAMPING_GAMMA];
        const double other_gamma = element2->GetProperties()[DAMPING_GAMMA];
        const double equiv_gamma = 0.5 * (my_gamma + other_gamma);
        double equiv_mass = (mRealMass*other_real_mass)/(mRealMass+other_real_mass);
        const double viscous_damping_coeff     = 2.0 * equiv_gamma * sqrt(equiv_mass * kn_el);
        //double viscous_damping_coeff = (1-equiv_coefficientOfRestitution) * 2.0 * sqrt(kn_el * equiv_mass);

        double rescaled_damping = viscous_damping_coeff/(2*equiv_mass);
        //double a = 1.4142-equiv_gamma*equiv_gamma;

        //double sqr_period = kn_el / equiv_mass - rescaled_damping*rescaled_damping;
        double sqr_period = sqrt(2.0) * kn_el / equiv_mass - rescaled_damping*rescaled_damping;   //esta es la correcta en continuu suponiendo un maximo de Kt= Kn
        return sqr_period;
    }

    bool DEMContinuumConstitutiveLaw::CheckRequirementsOfStressTensor() {

        return false;

    }

} //kratos
