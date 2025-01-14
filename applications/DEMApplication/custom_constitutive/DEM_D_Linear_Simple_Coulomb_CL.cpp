//  Kratos Multi-Physics - DEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//
// References using this model:
// R.L. Rangel et al. (2024). Multiscale data-driven modeling of the thermomechanical behavior of granular media with thermal expansion effects. Comput Geotech, 176:106789.
// R.L. Rangel et al. (2024). A continuum--discrete multiscale methodology using machine learning for thermal analysis of granular media. Comput Geotech, 168:106118.
// S. Zhao et al. (2020). Multiscale modeling of thermo-mechanical responses of granular materials: A hierarchical continuum--discrete coupling approach. CMAME, 367:113100.
//

#include "DEM_D_Linear_Simple_Coulomb_C.h"
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
    // Compute stiffness coefficients
    const double my_radius    = element1->GetRadius();
    const double other_radius = element2->GetRadius();
    const double equiv_radius = 2.0 * my_radius * other_radius / (my_radius + other_radius);
    const double my_young     = element1->GetYoung();
    const double my_poisson   = element1->GetPoisson();
    
    mKn = my_young * equiv_radius; // normal
    mKt = my_poisson * mKn;        // tangent

    // Compute normal and tangent forces
    const double normal_contact_force = mKn * indentation;
    LocalElasticContactForce[2] = normal_contact_force;

    CalculateTangentialForceWithNeighbour(normal_contact_force, OldLocalElasticContactForce, LocalElasticContactForce, LocalDeltDisp, sliding, element1, element2);
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
    // Compute stiffness coefficients
    const double my_radius  = element->GetRadius();
    const double my_young   = element->GetYoung();
    const double my_poisson = element->GetPoisson();

    mKn = my_young * my_radius; // normal
    mKt = my_poisson * mKn;     // tangent

    // Compute normal and tangent forces
    const double normal_contact_force = mKn * indentation;
    LocalElasticContactForce[2] = normal_contact_force;

    CalculateTangentialForceWithNeighbour(normal_contact_force, OldLocalElasticContactForce, LocalElasticContactForce, LocalDeltDisp, sliding, element, wall);
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
    LocalElasticContactForce[0] = OldLocalElasticContactForce[0] - mKt * LocalDeltDisp[0];
    LocalElasticContactForce[1] = OldLocalElasticContactForce[1] - mKt * LocalDeltDisp[1];
    const double tangent_contact_force = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);

    // Compute maximum admisible shear force
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
