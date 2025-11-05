//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//
// References using this model:
// R.L. Rangel et al. (2024). Multiscale data-driven modeling of the thermomechanical behavior of granular media with thermal expansion effects. Computers and Geotechnics, 176:106789.
// R.L. Rangel et al. (2024). A continuum-discrete multiscale methodology using machine learning for thermal analysis of granular media. Computers and Geotechnics, 168:106118.
// Guo & Zhao (2014). coupled FEM/DEM approach for hierarchical multiscale modelling of granular media. IJNME, 99(11):789-818.
//

#include "DEM_D_Linear_Simple_Coulomb_CL.h"
#include "custom_elements/spheric_particle.h"

namespace Kratos {

  DEMDiscontinuumConstitutiveLaw::Pointer DEM_D_Linear_Simple_Coulomb::Clone() const {
    DEMDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_D_Linear_Simple_Coulomb(*this));
    return p_clone;
  }

  std::unique_ptr<DEMDiscontinuumConstitutiveLaw> DEM_D_Linear_Simple_Coulomb::CloneUnique() {
    return Kratos::make_unique<DEM_D_Linear_Simple_Coulomb>();
  }

  std::string DEM_D_Linear_Simple_Coulomb::GetTypeOfLaw() {
    std::string type_of_law = "Linear";
    return type_of_law;
  }

  void DEM_D_Linear_Simple_Coulomb::Check(Properties::Pointer pProp) const {
    if (!pProp->Has(STATIC_FRICTION)) {
      if (!pProp->Has(FRICTION)) { //deprecated since April 6th, 2020
        KRATOS_WARNING("DEM") << std::endl;
        KRATOS_WARNING("DEM") << "WARNING: Variable STATIC_FRICTION or FRICTION should be present in the properties when using DEMDiscontinuumConstitutiveLaw. 0.0 value assigned by default." << std::endl;
        KRATOS_WARNING("DEM") << std::endl;
        pProp->GetValue(STATIC_FRICTION) = 0.0;
      }
      else {
        pProp->GetValue(STATIC_FRICTION) = pProp->GetValue(FRICTION);
      }
    }
  }

  void DEM_D_Linear_Simple_Coulomb::CalculateForces(const ProcessInfo& r_process_info,
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
    // Element properties
    const double r1 = element1->GetRadius();
    const double r2 = element2->GetRadius();
    const double m1 = element1->GetMass();
    const double m2 = element2->GetMass();
    const double E1 = element1->GetYoung();
    const double E2 = element2->GetYoung();
    const double v1 = element1->GetPoisson();
    const double v2 = element2->GetPoisson();

    // Effective properties
    const double reff = r1 * r2 / (r1 + r2);
    const double meff = m1 * m2 / (m1 + m2);
    const double Eeff = 1.0 / ((1.0 - v1*v1) / E1 + (1.0 - v2*v2) / E2);
    
    // Compute normal force elastic
    const double Kn = 2.0 * Eeff * sqrt(reff * indentation);
    const double Fne = (2.0/3.0) * Kn * indentation;

    // Compute normal force viscous
    Properties& properties_of_this_contact = element1->GetProperties().GetSubProperties(element2->GetProperties().Id());
    const double phi = properties_of_this_contact[DAMPING_GAMMA];
    double Fnv = -(2.0 * phi * sqrt(meff * Kn)) * LocalRelVel[2];

    // Check for artificial cohesion
    double Fn = Fne + Fnv;
    if (Fn < 0.0) {
        Fnv = -Fne;
    }

    // Store forces in their respective arrays
    LocalElasticContactForce[2] = Fne;
    ViscoDampingLocalContactForce[2] = Fnv;

    // Compute tangential force
    // TODO...

    // Calculate elastic energy (each particle in a contact with another particle receives half the contact energy)
    double& elastic_energy = element1->GetElasticEnergy();
    elastic_energy += 0.20 * LocalElasticContactForce[2] * indentation;;
  }

  void DEM_D_Linear_Simple_Coulomb::CalculateForcesWithFEM(const ProcessInfo& r_process_info,
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
    // Element properties
    const double E1 = element->GetYoung();
    const double E2 = wall->GetProperties()[YOUNG_MODULUS];
    const double v1 = element->GetPoisson();
    const double v2 = wall->GetProperties()[POISSON_RATIO];

    // Effective properties
    const double reff = element->GetRadius();
    const double meff = element->GetMass();
    const double Eeff = 1.0 / ((1.0 - v1*v1) / E1 + (1.0 - v2*v2) / E2);

    // Compute normal force elastic
    const double Kn = 2.0 * Eeff * sqrt(reff * indentation);
    const double Fne = (2.0/3.0) * Kn * indentation;

    // Compute normal force viscous
    Properties& properties_of_this_contact = element->GetProperties().GetSubProperties(wall->GetProperties().Id());
    const double phi = properties_of_this_contact[DAMPING_GAMMA];
    double Fnv = -(2.0 * phi * sqrt(meff * Kn)) * LocalRelVel[2];

    // Check for artificial cohesion
    double Fn = Fne + Fnv;
    if (Fn < 0.0) {
        Fnv = -Fne;
    }

    // Store forces in their respective arrays
    LocalElasticContactForce[2] = Fne;
    ViscoDampingLocalContactForce[2] = Fnv;

    // Compute tangential force
    // TODO...

    // Calculate elastic energy (each particle in a contact with a wall receives all the contact energy)
    double& elastic_energy = element->GetElasticEnergy();
    elastic_energy += 0.40 * LocalElasticContactForce[2] * indentation;;
    
  }

  template<class NeighbourClassType>
  void DEM_D_Linear_Simple_Coulomb::CalculateTangentialForceWithNeighbour(const double normal_contact_force,
                                                                          const double OldLocalElasticContactForce[3],
                                                                          double LocalElasticContactForce[3],
                                                                          const double LocalDeltDisp[3],
                                                                          bool& sliding,
                                                                          SphericParticle* const element,
                                                                          NeighbourClassType* const neighbour) {
    // Compute shear force
    LocalElasticContactForce[0] = OldLocalElasticContactForce[0] - this->mKt * LocalDeltDisp[0];
    LocalElasticContactForce[1] = OldLocalElasticContactForce[1] - this->mKt * LocalDeltDisp[1];
    const double tangent_contact_force = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);

    // Compute maximum admissible shear force
    Properties& properties_of_this_contact = element->GetProperties().GetSubProperties(neighbour->GetProperties().Id());
    const double friction_angle_tg = std::tan(properties_of_this_contact[STATIC_FRICTION]);
    const double MaximumAdmisibleShearForce = normal_contact_force * friction_angle_tg;

    // Check for sliding: apply Coulomb friction condition
    if (tangent_contact_force > MaximumAdmisibleShearForce) {
        sliding = true;
        const double fraction = MaximumAdmisibleShearForce / tangent_contact_force;
        LocalElasticContactForce[0] *= fraction;
        LocalElasticContactForce[1] *= fraction;
    }
  }

} // namespace Kratos
