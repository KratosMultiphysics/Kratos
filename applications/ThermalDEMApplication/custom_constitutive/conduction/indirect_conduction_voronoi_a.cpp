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
#include "indirect_conduction_voronoi_a.h"

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  IndirectConductionVoronoiA::IndirectConductionVoronoiA() {}
  IndirectConductionVoronoiA::~IndirectConductionVoronoiA() {}

  //------------------------------------------------------------------------------------------------------------
  double IndirectConductionVoronoiA::GetSearchDistance(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    return particle->GetParticleRadius() * r_process_info[MAX_CONDUCTION_DISTANCE];
  }

  //------------------------------------------------------------------------------------------------------------
  double IndirectConductionVoronoiA::ComputeHeatFlux(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
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
  double IndirectConductionVoronoiA::SphereWallCoeff(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    NumericalIntegrationMethod& integ  = particle->GetNumericalIntegrationMethod();
    const double fluid_conductivity    = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
    const double particle_conductivity = particle->GetParticleConductivity();
    const double distance              = particle->mNeighborDistanceAdjusted;
    const double contact_radius        = particle->mContactRadiusAdjusted;
    const double particle_radius       = particle->GetParticleRadius();
    double rij, upp_lim;

    // Get radius of voronoi cell face
    rij = particle->GetVoronoiCellFaceRadius(r_process_info);
    if (rij <= contact_radius)
      return 0.0;

    // Compute upper limit of integral
    upp_lim = particle_radius * rij / sqrt(rij * rij + distance * distance);

    // Fill integration parameters
    integ.CleanParameters();
    integ.mpEvalIntegrand = &IndirectConductionVoronoiA::EvalIntegrandVoronoiWall;
    integ.mLimMin         = contact_radius;
    integ.mLimMax         = upp_lim;
    integ.mParams.p1      = distance;
    integ.mParams.p2      = fluid_conductivity;
    integ.mParams.p3      = particle_conductivity;
    integ.mParams.p4      = particle_radius;
    integ.mParams.p5      = rij;

    // Heat transfer coefficient from integral expression solved numerically
    return integ.SolveIntegral();

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double IndirectConductionVoronoiA::SphereSphereMonoSizeCoeff(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    NumericalIntegrationMethod& integ   = particle->GetNumericalIntegrationMethod();
    const double fluid_conductivity     = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
    const double effective_conductivity = particle->ComputeEffectiveConductivity();
    const double distance               = particle->mNeighborDistanceAdjusted;
    const double contact_radius         = particle->mContactRadiusAdjusted;
    const double particle_radius        = particle->GetParticleRadius();
    double rij, upp_lim;

    // Get radius of voronoi cell face
    rij = particle->GetVoronoiCellFaceRadius(r_process_info);
    if (rij <= contact_radius)
      return 0.0;

    // Compute upper limit of integral
    upp_lim = particle_radius * rij / sqrt(rij * rij + distance * distance / 4.0);

    // Fill integration parameters
    integ.CleanParameters();
    integ.mpEvalIntegrand = &IndirectConductionVoronoiA::EvalIntegrandVoronoiMono;
    integ.mLimMin         = contact_radius;
    integ.mLimMax         = upp_lim;
    integ.mParams.p1      = distance;
    integ.mParams.p2      = fluid_conductivity;
    integ.mParams.p3      = effective_conductivity;
    integ.mParams.p4      = particle_radius;
    integ.mParams.p5      = rij;

    // Heat transfer coefficient from integral expression solved numerically
    return integ.SolveIntegral();

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double IndirectConductionVoronoiA::SphereSphereMultiSizeCoeff(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) {
    KRATOS_TRY

    NumericalIntegrationMethod& integ  = particle->GetNumericalIntegrationMethod();
    const double fluid_conductivity    = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
    const double particle_conductivity = particle->GetParticleConductivity();
    const double neighbor_conductivity = particle->GetNeighborConductivity();
    const double distance              = particle->mNeighborDistanceAdjusted;
    const double contact_radius        = particle->mContactRadiusAdjusted;
    const double particle_radius       = particle->GetParticleRadius();
    const double neighbor_radius       = particle->GetNeighborRadius();
    double rij, rij_, upp_lim, D1, D2;

    // Get radius of voronoi cell face
    rij = particle->GetVoronoiCellFaceRadius(r_process_info);
    if (rij <= contact_radius)
      return 0.0;

    if (particle->mNeighborInContact)
      D1 = sqrt(particle_radius * particle_radius - contact_radius * contact_radius);
    else
      D1 = (particle_radius * particle_radius - neighbor_radius * neighbor_radius + distance * distance) / (2 * distance);

    D2 = distance - D1;

    if (particle_radius <= neighbor_radius)
      upp_lim = particle_radius * rij / sqrt(rij * rij + D1 * D1);
    else
      upp_lim = neighbor_radius * rij / sqrt(rij * rij + D2 * D2);

    rij_ = D2 * upp_lim / sqrt(neighbor_radius * neighbor_radius - upp_lim * upp_lim);

    // Fill integration parameters
    integ.CleanParameters();
    integ.mpEvalIntegrand = &IndirectConductionVoronoiA::EvalIntegrandVoronoiMulti;
    integ.mLimMin         = contact_radius;
    integ.mLimMax         = upp_lim;
    integ.mParams.p1      = distance;
    integ.mParams.p2      = fluid_conductivity;
    integ.mParams.p3      = particle_conductivity;
    integ.mParams.p4      = neighbor_conductivity;
    integ.mParams.p5      = particle_radius;
    integ.mParams.p6      = neighbor_radius;
    integ.mParams.p7      = rij;
    integ.mParams.p8      = rij_;
    integ.mParams.p9      = D1;
    integ.mParams.p10     = D2;

    // Heat transfer coefficient from integral expression solved numerically
    return integ.SolveIntegral();

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double IndirectConductionVoronoiA::EvalIntegrandVoronoiWall(NumericalIntegrationMethod* method) {
    KRATOS_TRY

    const double r   = method->mCoord;
    const double d   = method->mParams.p1;
    const double kf  = method->mParams.p2;
    const double kp  = method->mParams.p3;
    const double rp  = method->mParams.p4;
    const double rij = method->mParams.p5;

    return 2.0 * Globals::Pi * r / ((sqrt(rp * rp - r * r) - r * d / rij) / kp + (d - sqrt(rp * rp - r * r)) / kf);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double IndirectConductionVoronoiA::EvalIntegrandVoronoiMono(NumericalIntegrationMethod* method) {
    KRATOS_TRY

    const double r    = method->mCoord;
    const double d    = method->mParams.p1;
    const double kf   = method->mParams.p2;
    const double keff = method->mParams.p3;
    const double rp   = method->mParams.p4;
    const double rij  = method->mParams.p5;

    return 2.0 * Globals::Pi * r / ((sqrt(rp * rp - r * r) - r * d / (2.0 * rij)) / keff + 2.0 * (d / 2.0 - sqrt(rp * rp - r * r)) / kf);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  double IndirectConductionVoronoiA::EvalIntegrandVoronoiMulti(NumericalIntegrationMethod* method) {
    KRATOS_TRY

    const double r    = method->mCoord;
    const double d    = method->mParams.p1;
    const double kf   = method->mParams.p2;
    const double k1   = method->mParams.p3;
    const double k2   = method->mParams.p4;
    const double r1   = method->mParams.p5;
    const double r2   = method->mParams.p6;
    const double rij  = method->mParams.p7;
    const double rij_ = method->mParams.p8;
    const double D1   = method->mParams.p9;
    const double D2   = method->mParams.p10;

    const double beta1 = sqrt(r1 * r1 - r * r);
    const double beta2 = sqrt(r2 * r2 - r * r);

    return 2.0 * Globals::Pi * r / ((beta1 - D1 * r / rij) / k1 + (beta2 - D2 * r / rij_) / k2 + (d - beta1 - beta2) / kf);

    KRATOS_CATCH("")
  }

} // namespace Kratos
