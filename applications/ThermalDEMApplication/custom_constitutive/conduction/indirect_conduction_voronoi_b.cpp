//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

// System includes

// External includes

// Project includes
#include "indirect_conduction_voronoi_b.h"

namespace Kratos {
  //------------------------------------------------------------------------------------------------------------
  IndirectConductionVoronoiB::IndirectConductionVoronoiB() {}
  IndirectConductionVoronoiB::~IndirectConductionVoronoiB() {}

  //------------------------------------------------------------------------------------------------------------
  double IndirectConductionVoronoiB::GetSearchDistance(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    return particle->GetParticleRadius() * r_process_info[MAX_CONDUCTION_DISTANCE];
  }

  //------------------------------------------------------------------------------------------------------------
  double IndirectConductionVoronoiB::ComputeHeatFlux(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Check if particles are close enough
    if (!particle->CheckSurfaceDistance(r_process_info[MAX_CONDUCTION_DISTANCE]))
      return 0.0;

    // Compute heat transfer coefficient
    double h = 0.0;

    if (particle->mNeighborType & WALL_NEIGHBOR)
      h = SphereWallCoeff(r_process_info, particle);

    else if (particle->GetParticleRadius() == particle->GetNeighborRadius())
      h = SphereSphereMonoSizeCoeff(r_process_info, particle);

    else
      h = SphereSphereMultiSizeCoeff(r_process_info, particle);

    // Compute heat flux
    return h * (particle->GetNeighborTemperature() - particle->GetParticleTemperature());

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double IndirectConductionVoronoiB::SphereWallCoeff(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    const double core                  = r_process_info[ISOTHERMAL_CORE_RADIUS];
    const double fluid_conductivity    = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
    const double particle_conductivity = particle->GetParticleConductivity();
    const double distance              = particle->mNeighborDistanceAdjusted;
    const double contact_radius        = particle->mContactRadiusAdjusted;
    const double particle_radius       = particle->GetParticleRadius();

    // Get radius of voronoi cell face
    const double rij = particle->GetVoronoiCellFaceRadius(r_process_info);
    if (rij <= contact_radius)
      return 0.0;

    const double rc = core * particle_radius;
    const double a  = (1.0 / rc - 1.0 / particle_radius) / (2.0 * particle_conductivity) + 1.0 / (2 * fluid_conductivity * particle_radius);
    const double b  = 1.0 / (2 * fluid_conductivity * distance);
    const double c0 = distance / sqrt(rij * rij + distance * distance);
    const double c1 = distance / sqrt(contact_radius * contact_radius + distance * distance);
    const double f  = (a - b * c0) / (a - b * c1);
    double ln = 0.0;
    if (f > 0.0)
      ln = log(f);

    // Heat transfer coefficient
    return Globals::Pi * ln / b;
    
    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double IndirectConductionVoronoiB::SphereSphereMonoSizeCoeff(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    const double core                   = r_process_info[ISOTHERMAL_CORE_RADIUS];
    const double fluid_conductivity     = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
    const double effective_conductivity = particle->ComputeEffectiveConductivity();
    const double distance               = particle->mNeighborDistanceAdjusted / 2.0;
    const double contact_radius         = particle->mContactRadiusAdjusted;
    const double particle_radius        = particle->GetParticleRadius();

    // Get radius of voronoi cell face
    const double rij = particle->GetVoronoiCellFaceRadius(r_process_info);
    if (rij <= contact_radius)
      return 0.0;
    
    const double rc = core * particle_radius;
    const double a  = (1.0 / rc - 1.0 / particle_radius) / (2.0 * effective_conductivity) + 1.0 / (fluid_conductivity * particle_radius);
    const double b  = 1.0 / (fluid_conductivity * distance);
    const double c0 = distance / sqrt(rij * rij + distance * distance);
    const double c1 = distance / sqrt(contact_radius * contact_radius + distance * distance);
    const double f  = (a - b * c0) / (a - b * c1);
    double ln = 0.0;
    if (f > 0.0)
      ln = log(f);

    // Heat transfer coefficient
    return Globals::Pi * ln / b;

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double IndirectConductionVoronoiB::SphereSphereMultiSizeCoeff(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    const double core                  = r_process_info[ISOTHERMAL_CORE_RADIUS];
    const double fluid_conductivity    = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
    const double particle_conductivity = particle->GetParticleConductivity();
    const double neighbor_conductivity = particle->GetNeighborConductivity();
    const double distance              = particle->mNeighborDistanceAdjusted;
    const double contact_radius        = particle->mContactRadiusAdjusted;
    const double particle_radius       = particle->GetParticleRadius();
    const double neighbor_radius       = particle->GetNeighborRadius();

    // Get radius of voronoi cell face
    const double rij = particle->GetVoronoiCellFaceRadius(r_process_info);
    if (rij <= contact_radius)
      return 0.0;

    const double An = Globals::Pi * rij * rij; // area of neighboring voronoi cells

    const double gamma1 = particle_radius / distance;
    const double gamma2 = neighbor_radius / distance;
    const double dgamma = gamma2 - gamma1;

    const double A = (particle_conductivity + fluid_conductivity * (1.0 / core - 1.0)) / (particle_conductivity * gamma1);
    const double B = (neighbor_conductivity + fluid_conductivity * (1.0 / core - 1.0)) / (neighbor_conductivity * gamma2);

    const double lambda = (1.0 + dgamma * A) * (1.0 - dgamma * B);

    const double delmax = 0.5 * (sqrt((4.0 * An) / (Globals::Pi * distance * distance * (1.0 - dgamma * dgamma)) + 1.0) - dgamma);
    const double delmin = 0.5 * (sqrt((4.0 * contact_radius * contact_radius) / (distance * distance * (1.0 - dgamma * dgamma)) + 1.0) - dgamma);

    const double Xmax = ((A + B) * delmax + dgamma * B - 1.0) / sqrt(std::abs(lambda));
    const double Xmin = ((A + B) * delmin + dgamma * B - 1.0) / sqrt(std::abs(lambda));

    const double Y1 = (Xmax - Xmin) / (1.0 - Xmax * Xmin);
    const double Y2 = (Xmax - Xmin) / (1.0 + Xmax * Xmin);

    // Heat transfer coefficient
    if (lambda > 0.0)
      return Globals::Pi * fluid_conductivity * distance * (1.0 - dgamma * dgamma) * log(std::abs((1.0 - Y1) / (1.0 + Y1))) / (2.0 * sqrt(std::abs(lambda)));
    else if (lambda < 0.0)
      return Globals::Pi * fluid_conductivity * distance * (1.0 - dgamma * dgamma) * atan(Y2) / (2.0 * sqrt(std::abs(lambda)));
    else
      return Globals::Pi * fluid_conductivity * distance * (1.0 - dgamma * dgamma) * (1.0 / delmin - 1.0 / delmax) / (A + B);

    KRATOS_CATCH("")
  }

} // namespace Kratos
