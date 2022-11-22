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
#include "indirect_conduction_surround_layer.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  IndirectConductionSurroundLayer::IndirectConductionSurroundLayer() {}
  IndirectConductionSurroundLayer::~IndirectConductionSurroundLayer() {}

  //------------------------------------------------------------------------------------------------------------
  double IndirectConductionSurroundLayer::GetSearchDistance(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    return particle->GetParticleRadius() * r_process_info[FLUID_LAYER_THICKNESS];
  }

  //------------------------------------------------------------------------------------------------------------
  double IndirectConductionSurroundLayer::ComputeHeatFlux(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    // Check if particles are close enough
    if (!particle->CheckSurfaceDistance(r_process_info[FLUID_LAYER_THICKNESS]))
      return 0.0;

    // Compute heat transfer coefficient
    double h = 0.0;

    if (particle->mNeighborType & WALL_NEIGHBOR)
      h = SphereWallCoeff(r_process_info, particle);

    else if (particle->mNeighborType & PARTICLE_NEIGHBOR)
      h = SphereSphereCoeff(r_process_info, particle);

    // Compute heat flux
    return h * (particle->GetNeighborTemperature() - particle->GetParticleTemperature());

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double IndirectConductionSurroundLayer::SphereWallCoeff(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    const double min_dist           = r_process_info[MIN_CONDUCTION_DISTANCE];
    const double fluid_conductivity = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
    const double fluid_layer        = r_process_info[FLUID_LAYER_THICKNESS];
    const double distance_0         = particle->mNeighborDistance;
    const double distance           = particle->mNeighborDistanceAdjusted;
    const double particle_radius    = particle->GetParticleRadius();
    double a, b, c, r_in, r_out;

    a = (distance - particle_radius) / particle_radius;

    if (distance_0 > particle_radius + min_dist)
      r_in = 0.0;
    else
      r_in = sqrt(1.0 - pow(min_dist / particle_radius - a - 1.0, 2.0));

    if (a > sqrt(pow((particle_radius + (fluid_layer * particle_radius)) / particle_radius, 2.0) - 1.0) - 1.0)
      r_out = sqrt(pow((particle_radius + (fluid_layer * particle_radius)) / particle_radius, 2.0) - pow(a + 1.0, 2.0));
    else
      r_out = 1.0;

    b = sqrt(1.0 - r_out * r_out);
    c = sqrt(1.0 - r_in  * r_in);

    // Heat transfer coefficient from analytical solution of the integral expression
    return 2.0 * Globals::Pi * fluid_conductivity * particle_radius * ((a + 1.0) * log(fabs(b - a - 1.0) / fabs(a - c + 1.0)) + b - c);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double IndirectConductionSurroundLayer::SphereSphereCoeff(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    NumericalIntegrationMethod& integ = particle->GetNumericalIntegrationMethod();
    const double min_dist             = r_process_info[MIN_CONDUCTION_DISTANCE];
    const double fluid_conductivity   = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
    const double fluid_layer          = r_process_info[FLUID_LAYER_THICKNESS];
    const double distance             = particle->mNeighborDistanceAdjusted;
    const double contact_radius       = particle->mContactRadiusAdjusted;
    const double particle_radius      = particle->GetParticleRadius();
    const double neighbor_radius      = particle->GetNeighborRadius();
    const double r_min                = std::min(particle_radius, neighbor_radius);
    const double r_max                = std::max(particle_radius, neighbor_radius);
    double upp_lim;

    // Compute upper limit of integral
    const double param = pow((r_max + (fluid_layer * r_max)), 2.0);

    if (distance <= sqrt(param - r_min * r_min))
      upp_lim = r_min;
    else
      upp_lim = sqrt(param - pow(((param - r_min * r_min + distance * distance) / (2.0 * distance)), 2.0));

    // Fill integration parameters
    integ.CleanParameters();
    integ.mpEvalIntegrand = &IndirectConductionSurroundLayer::EvalIntegrandSurrLayer;
    integ.mLimMin         = contact_radius;
    integ.mLimMax         = upp_lim;
    integ.mParams.p1      = distance;
    integ.mParams.p2      = min_dist;
    integ.mParams.p3      = particle_radius;
    integ.mParams.p4      = neighbor_radius;

    // Heat transfer coefficient from integral expression solved numerically
    return fluid_conductivity * integ.SolveIntegral();

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double IndirectConductionSurroundLayer::EvalIntegrandSurrLayer(NumericalIntegrationMethod* method) {
    KRATOS_TRY

    const double r    = method->mCoord;
    const double d    = method->mParams.p1;
    const double dmin = method->mParams.p2;
    const double r1   = method->mParams.p3;
    const double r2   = method->mParams.p4;

    return 2.0 * Globals::Pi * r / std::max(dmin, d - sqrt(r1 * r1 - r * r) - sqrt(r2 * r2 - r * r));

    KRATOS_CATCH("")
  }

} // namespace Kratos
