//
//   Project Name:                     ThermalDEM $
//   Last Modified by:    $Author: Ferran Arrufat $
//   Date:                $Date:    February 2015 $
//   Revision:            $Revision:      1.0.0.0 $
//

// System includes
#include <string>
#include <iostream>
#include <limits>

// Project includes
#include "thermal_spheric_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{
  // Constructor/Destructor methods

  template <class TBaseElement>
  ThermalSphericParticle<TBaseElement>::~ThermalSphericParticle() {}

  // Initialization methods

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::Initialize(const ProcessInfo& r_process_info) {
    // Initialize base class
    TBaseElement::Initialize(r_process_info);

    // Get material properties
    mSpecificHeat        = GetProperties()[SPECIFIC_HEAT];
    mThermalConductivity = GetProperties()[THERMAL_CONDUCTIVITY];

    // Set thermal flags
    if (r_process_info[DIRECT_CONDUCTION_OPTION])   this->Set(DEMFlags::HAS_DIRECT_CONDUCTION,   true);
    else                                            this->Set(DEMFlags::HAS_DIRECT_CONDUCTION,   false);

    if (r_process_info[INDIRECT_CONDUCTION_OPTION]) this->Set(DEMFlags::HAS_INDIRECT_CONDUCTION, true);
    else                                            this->Set(DEMFlags::HAS_INDIRECT_CONDUCTION, false);

    if (r_process_info[CONVECTION_OPTION])          this->Set(DEMFlags::HAS_CONVECTION,          true);
    else                                            this->Set(DEMFlags::HAS_CONVECTION,          false);

    if (r_process_info[CONVECTION_OPTION])          this->Set(DEMFlags::HAS_RADIATION,           true);
    else                                            this->Set(DEMFlags::HAS_RADIATION,           false);
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::InitializeSolutionStep(const ProcessInfo& r_process_info) {
    // Initialize base class
    TBaseElement::InitializeSolutionStep(r_process_info);
  }

  // Calculate right hand side (forces and heat fluxes)

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::CalculateRightHandSide(const ProcessInfo& r_process_info, double dt, const array_1d<double, 3>& gravity) {
    // Force components
    TBaseElement::CalculateRightHandSide(r_process_info, dt, gravity);
    
    // Heat flux components
    ComputeHeatFluxes(r_process_info);
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeHeatFluxes(const ProcessInfo& r_process_info) {
    // Initialize heat flux contributions
    mConductiveHeatFlux = 0.0;
    mConvectiveHeatFlux = 0.0;
    mRadiativeHeatFlux  = 0.0;
    mTotalHeatFlux      = 0.0;

    // Direct conduction
    if (this->Is(DEMFlags::HAS_DIRECT_CONDUCTION)) {
      ComputeBallToBallDirectConductionHeatFlux(r_process_info);
      ComputeBallToRigidFaceDirectConductionHeatFlux(r_process_info);
    }

    // Indirect conduction
    if (this->Is(DEMFlags::HAS_INDIRECT_CONDUCTION)) {
      ComputeBallToBallIndirectConductionHeatFlux(r_process_info);
      ComputeBallToRigidFaceIndirectConductionHeatFlux(r_process_info);
    }

    // Convection
    if (this->Is(DEMFlags::HAS_CONVECTION)) {
      ComputeConvectiveHeatFlux(r_process_info);
    }

    // Radiation
    if (this->Is(DEMFlags::HAS_RADIATION)) {
      ComputeRadiativeHeatFlux(r_process_info);
    }

    // Sum up contributions
    mTotalHeatFlux = mConductiveHeatFlux + mConvectiveHeatFlux + mRadiativeHeatFlux;
  }

  // Compute heat fluxes components

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeBallToBallDirectConductionHeatFlux(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Direct conduction model
    std::string model = r_process_info[DIRECT_CONDUCTION_MODEL];

    // Loop over neighbour particles
    for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
      if (mNeighbourElements[i] == NULL) continue;
      ThermalSphericParticle<TBaseElement>* neighbour_iterator = dynamic_cast<ThermalSphericParticle<TBaseElement>*>(mNeighbourElements[i]);

      // Get particles radii
      double this_radius  = GetRadius();
      double other_radius = neighbour_iterator->GetRadius();

      // Compute direction and distance between centroids
      array_1d<double, 3> direction;
      direction[0] = GetGeometry()[0].Coordinates()[0] - neighbour_iterator->GetGeometry()[0].Coordinates()[0];
      direction[1] = GetGeometry()[0].Coordinates()[1] - neighbour_iterator->GetGeometry()[0].Coordinates()[1];
      direction[2] = GetGeometry()[0].Coordinates()[2] - neighbour_iterator->GetGeometry()[0].Coordinates()[2];

      double distance = DEM_MODULUS_3(direction);

      // Check if particles are in contact
      if (distance >= this_radius + other_radius)
        continue;

      // Get particles temperatures
      double this_temp  = GetParticleTemperature();
      double other_temp = neighbour_iterator->GetParticleTemperature();
      double temp_grad  = other_temp - this_temp;

      // Get common properties
      double other_conductivity = neighbour_iterator->mThermalConductivity;

      // Compute heat flux according to selected model
      if (model.compare("batchelor_obrien") == 0) {

        // Compute effective thermal conductivity
        double eff_conductivity = mThermalConductivity * other_conductivity / (mThermalConductivity + other_conductivity);

        // Compute contact radius
        double contact_radius = sqrt(fabs(this_radius*this_radius - pow(((this_radius*this_radius - other_radius*other_radius + distance*distance) / (2*distance)),2)));

        // Compute heat flux
        mConductiveHeatFlux += 4 * eff_conductivity * contact_radius * temp_grad;
      }
      else if (model.compare("thermal_pipe") == 0) {

        // Compute average thermal conductivity
        double avg_conductivity = (this_radius + other_radius) / (this_radius/mThermalConductivity + other_radius/other_conductivity);

        // Compute contact area
        double contact_radius = sqrt(fabs(this_radius*this_radius - pow(((this_radius*this_radius - other_radius*other_radius + distance*distance) / (2*distance)),2)));
        double contact_area   = Globals::Pi*contact_radius*contact_radius;

        // Compute heat flux
        mConductiveHeatFlux += avg_conductivity * contact_area * temp_grad / distance;
      }
      else if (model.compare("collisional") == 0) {
        
        // Get properties
        double this_density       = GetDensity();
        double this_mass          = GetMass();
        double this_poisson       = GetPoisson();
        double this_young         = GetYoung();
        double other_density      = neighbour_iterator->GetDensity();
        double other_mass         = neighbour_iterator->GetMass();
        double other_poisson      = neighbour_iterator->GetPoisson();
        double other_young        = neighbour_iterator->GetYoung();
        double other_heatcapacity = neighbour_iterator->mSpecificHeat;

        // Compute effective parameters
        double eff_radius = this_radius * other_radius / (this_radius + other_radius);
        double eff_mass   = this_mass * other_mass / (this_mass + other_mass);
        double eff_young  = 1 / ((1-this_poisson*this_poisson)/this_young + (1-other_poisson*other_poisson)/other_young);

        // Compute expected collision time
        double impact_normal_velocity = 0.0;   // TODO: save impact normal velocity
        double col_time = 0.0;                 // TODO: track collision time
        double expect_col_time = 2.87 * pow(eff_mass*eff_mass / (eff_radius*eff_young*eff_young*impact_normal_velocity),1/5);

        // Check if collision time is smaller than expected value, otherwise use static model (batchelor_obrien)
        if (col_time < expect_col_time && impact_normal_velocity != 0.0) {
          // Compute max contact radius
          double contact_radius_max = pow(15*eff_radius*eff_mass*impact_normal_velocity*impact_normal_velocity / (16*eff_young),1/5);

          // Compute coefficients
          double a1 = this_density * mSpecificHeat;
          double a2 = other_density * other_heatcapacity;
          double b1 = a1 * mThermalConductivity;
          double b2 = a2 * other_conductivity;
          double c  = a1/a2;

          // Fourier number (Assumption: average of both particles)
          double fo1 = mThermalConductivity * expect_col_time / (a1*contact_radius_max*contact_radius_max);
          double fo2 = other_conductivity   * expect_col_time / (a2*contact_radius_max*contact_radius_max);
          double fo  = (fo1+fo2)/2;

          // Compute coefficient C
          double c1 = -2.300*c*c +  8.909*c -  4.235;
          double c2 =  8.169*c*c - 33.770*c + 24.885;
          double c3 = -5.758*c*c + 24.464*c - 20.511;

          double C_coeff = 0.435 * (sqrt(c2*c2-4*c1*(c3-fo)) - c2) / c1;

          // Compute heat flux
          mConductiveHeatFlux += C_coeff * Globals::Pi * contact_radius_max*contact_radius_max * pow(expect_col_time,-1/2) * temp_grad / (pow(b1,-1/2) + pow(b2,-1/2));
        }
        else {
          // Compute effective thermal conductivity
          double eff_conductivity = mThermalConductivity * other_conductivity / (mThermalConductivity + other_conductivity);

          // Compute contact radius
          double contact_radius = sqrt(fabs(this_radius*this_radius - pow(((this_radius*this_radius - other_radius*other_radius + distance*distance) / (2*distance)),2)));

          // Compute heat flux
          mConductiveHeatFlux += 4 * eff_conductivity * contact_radius * temp_grad;
        }
      }
    }

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeBallToRigidFaceDirectConductionHeatFlux(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Direct conduction model
    std::string model = r_process_info[DIRECT_CONDUCTION_MODEL];

    // Loop over neighbour walls
    std::vector<DEMWall*>& rNeighbours = mNeighbourRigidFaces;
    for (unsigned int i = 0; i < rNeighbours.size(); i++) {
      DEMWall* wall = rNeighbours[i];
      if (wall == NULL) continue;
      
      // Compute distance between particle and wall
      // (stolen from ComputeBallToRigidFaceContactForce in spheric_particle)
      double dummy1[3][3];
      DEM_SET_COMPONENTS_TO_ZERO_3x3(dummy1);
      array_1d<double, 4> dummy2 = this->mContactConditionWeights[i];
      array_1d<double, 3> dummy3 = ZeroVector(3);
      array_1d<double, 3> dummy4 = ZeroVector(3);
      double distance = 0.0;
      int ContactType = -1;
      wall->ComputeConditionRelativeData(i, this, dummy1, distance, dummy2, dummy3, dummy4, ContactType);

      // Check if particle is in contact with wall
      if (ContactType != 1 && ContactType != 2 && ContactType != 3)
        continue;

      // Get temperatures
      // (Assumption: wall temperature is the average of its nodes)
      double wall_temp = 0.0;
      double n_nodes = wall->GetGeometry().size();
      for (unsigned int i = 0; i < n_nodes; i++) {
        double node_temp = wall->GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE);
        wall_temp += node_temp;
      }
      wall_temp /= n_nodes;

      double particle_temp = GetParticleTemperature();
      double temp_grad = wall_temp - particle_temp;

      // Get common properties
      double particle_radius = GetRadius();
      double wall_conductivity = wall->GetProperties()[THERMAL_CONDUCTIVITY];

      // Compute heat flux according to selected model
      if (model.compare("batchelor_obrien") == 0) {

        // Compute effective thermal conductivity
        double eff_conductivity = mThermalConductivity * wall_conductivity / (mThermalConductivity + wall_conductivity);

        // Compute contact radius
        double indentation    = particle_radius - distance;
        double contact_radius = sqrt(indentation * (2 * particle_radius - indentation));

        // Compute heat flux
        mConductiveHeatFlux += 4 * eff_conductivity * contact_radius * temp_grad;
      }
      else if (model.compare("thermal_pipe") == 0) {
        
        // Compute average thermal conductivity
        // (Assumption: average conductivity considers particle only)
        double avg_conductivity = mThermalConductivity;

        // Compute contact area
        double indentation    = particle_radius - distance;
        double contact_radius = sqrt(indentation * (2 * particle_radius - indentation));
        double contact_area   = Globals::Pi*contact_radius*contact_radius;

        // Compute heat flux
        // (Assumption: twice the particle-wall distance as if the wall was another particle)
        mConductiveHeatFlux += avg_conductivity * contact_area * temp_grad / (2*distance);

      }
      else if (model.compare("collisional") == 0) {
        // Get properties
        double particle_density  = GetDensity();
        double particle_mass     = GetMass();
        double particle_poisson  = GetPoisson();
        double particle_young    = GetYoung();
        double wall_density      = wall->GetProperties()[DENSITY];
        double wall_poisson      = wall->GetPoisson();
        double wall_young        = wall->GetYoung();
        double wall_heatcapacity = wall->GetProperties()[SPECIFIC_HEAT];

        // Compute effective parameters
        double eff_radius = particle_radius;
        double eff_mass   = particle_mass;
        double eff_young  = 1 / ((1-particle_poisson*particle_poisson)/ particle_young + (1-wall_poisson*wall_poisson)/wall_young);

        // Compute expected collision time
        double impact_normal_velocity = 0.0;   // TODO: save impact normal velocity
        double col_time = 0.0;                 // TODO: track collision time
        double expect_col_time = 2.87 * pow(eff_mass*eff_mass / (eff_radius*eff_young*eff_young*impact_normal_velocity),1/5);

        // Check if collision time is smaller than expected value, otherwise use static model (batchelor_obrien)
        if (col_time < expect_col_time && impact_normal_velocity != 0.0) {
          // Compute max contact radius
          double contact_radius_max = pow(15*eff_radius*eff_mass*impact_normal_velocity*impact_normal_velocity / (16*eff_young),1/5);

          // Compute coefficients
          double a1 = particle_density * mSpecificHeat;
          double a2 = wall_density * wall_heatcapacity;
          double b1 = a1 * mThermalConductivity;
          double b2 = a2 * wall_conductivity;
          double c  = a1/a2;

          // Fourier number
          double fo = mThermalConductivity * expect_col_time / (a1*contact_radius_max*contact_radius_max);

          // Compute coefficient C
          double c1 = -2.300*c*c +  8.909*c -  4.235;
          double c2 =  8.169*c*c - 33.770*c + 24.885;
          double c3 = -5.758*c*c + 24.464*c - 20.511;

          double C_coeff = 0.435 * (sqrt(c2*c2-4*c1*(c3-fo)) - c2) / c1;

          // Compute heat flux
          mConductiveHeatFlux += C_coeff * Globals::Pi * contact_radius_max*contact_radius_max * pow(expect_col_time,-1/2) * temp_grad / (pow(b1,-1/2) + pow(b2,-1/2));
        }
        else {
          // Compute effective thermal conductivity
          double eff_conductivity = mThermalConductivity * wall_conductivity / (mThermalConductivity + wall_conductivity);

          // Compute contact radius
          double indentation    = particle_radius - distance;
          double contact_radius = sqrt(indentation * (2 * particle_radius - indentation));

          // Compute heat flux
          mConductiveHeatFlux += 4 * eff_conductivity * contact_radius * temp_grad;
        }
      }
    }

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeBallToBallIndirectConductionHeatFlux(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Indirect conduction model
    std::string model = r_process_info[INDIRECT_CONDUCTION_MODEL];

    // Loop over neighbour particles
    for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
      if (mNeighbourElements[i] == NULL) continue;
      ThermalSphericParticle<TBaseElement>* neighbour_iterator = dynamic_cast<ThermalSphericParticle<TBaseElement>*>(mNeighbourElements[i]);

      // Get particles temperatures
      double this_temp  = GetParticleTemperature();
      double other_temp = neighbour_iterator->GetParticleTemperature();
      double temp_grad  = other_temp - this_temp;

      // Get particles radii
      double this_radius  = GetRadius();
      double other_radius = neighbour_iterator->GetRadius();
      
      // Get interstitial fluid properties
      double fluid_conductivity = r_process_info[FLUID_THERMAL_CONDUCTIVITY];

      // Compute heat flux according to selected model
      if (model.compare("surrounding_layer") == 0) {
        
        // Get model parameters
        double layer    = r_process_info[FLUID_LAYER_THICKNESS];
        double min_dist = r_process_info[MIN_CONDUCTION_DISTANCE];
        
        // Compute direction and distance between centroids
        array_1d<double, 3> direction;
        direction[0] = GetGeometry()[0].Coordinates()[0] - neighbour_iterator->GetGeometry()[0].Coordinates()[0];
        direction[1] = GetGeometry()[0].Coordinates()[1] - neighbour_iterator->GetGeometry()[0].Coordinates()[1];
        direction[2] = GetGeometry()[0].Coordinates()[2] - neighbour_iterator->GetGeometry()[0].Coordinates()[2];

        double distance = DEM_MODULUS_3(direction);

        // Check if particles are close enough
        if (distance > this_radius + other_radius + layer * std::max(this_radius, other_radius))
          continue;

        // Compute heat transfer coefficient
        double h = 0.0;
        
        if (this_radius == other_radius) {
          
          double a = (distance - 2*this_radius) / this_radius;
          double r_in;
          double r_out;

          if (distance > 2*this_radius + min_dist)
            r_in = 0.0;
          else
            r_in = sqrt(1-pow(min_dist/this_radius-a-1,2));

          if (a > sqrt(pow((this_radius + (layer*this_radius)) / this_radius,2) - 1) - 1)
            r_out = sqrt(pow((this_radius + (layer*this_radius)) / this_radius,2) - (a+1)*(a+1));
          else
            r_out = 1.0;

          double b = sqrt(1-r_out*r_out);
          double c = sqrt(1-r_in*r_in);

          h = 2 * Globals::Pi * fluid_conductivity * this_radius * ((a+1) * log(abs((b-a-1) / (a-c+1))) + b - c);
        }
        else {

          // Compute lower limit of integral (contact radius)
          double low_lim;

          if (distance < this_radius + other_radius)
            low_lim = sqrt(fabs(this_radius*this_radius - pow(((this_radius*this_radius - other_radius*other_radius + distance*distance) / (2*distance)),2)));
          else
            low_lim = 0.0;

          // Compute upper limit of integral
          double r_min = std::min(this_radius,other_radius);
          double r_max = std::max(this_radius,other_radius);
          double param = pow((r_max+(layer*r_max)),2);
          double upp_lim;
          
          if (distance <= sqrt(param-r_min*r_min))
            upp_lim = r_min;
          else
            upp_lim = sqrt(param - pow(((param - r_min*r_min + distance*distance) / (2*distance)),2));

          // Heat transfer coefficient from integral expression solved numerically
          h = fluid_conductivity * IntegralSurrLayer(r_process_info, distance, this_radius, other_radius, low_lim, upp_lim);
        }

        // Compute heat flux
        mConductiveHeatFlux += h * temp_grad;
      }
    }

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeBallToRigidFaceIndirectConductionHeatFlux(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Indirect conduction model
    std::string model = r_process_info[INDIRECT_CONDUCTION_MODEL];

    // Loop over neighbour walls
    std::vector<DEMWall*>& rNeighbours = mNeighbourRigidFaces;
    for (unsigned int i = 0; i < rNeighbours.size(); i++) {
      DEMWall* wall = rNeighbours[i];
      if (wall == NULL) continue;
      
      // Get temperatures
      // (Assumption: wall temperature is the average of its nodes)
      double wall_temp = 0.0;
      double n_nodes = wall->GetGeometry().size();
      for (unsigned int i = 0; i < n_nodes; i++) {
        double node_temp = wall->GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE);
        wall_temp += node_temp;
      }
      wall_temp /= n_nodes;

      double particle_temp = GetParticleTemperature();
      double temp_grad = wall_temp - particle_temp;

      // Get particles radii
      double particle_radius  = GetRadius();
      
      // Get interstitial fluid properties
      double fluid_conductivity = r_process_info[FLUID_THERMAL_CONDUCTIVITY];

      // Compute heat flux according to selected model
      if (model.compare("surrounding_layer") == 0) {

        // Get model parameters
        double layer    = r_process_info[FLUID_LAYER_THICKNESS];
        double min_dist = r_process_info[MIN_CONDUCTION_DISTANCE];
        
        // Compute distance between particle and wall
        // (stolen from ComputeBallToRigidFaceContactForce in spheric_particle)
        double dummy1[3][3];
        DEM_SET_COMPONENTS_TO_ZERO_3x3(dummy1);
        array_1d<double, 4> dummy2 = this->mContactConditionWeights[i];
        array_1d<double, 3> dummy3 = ZeroVector(3);
        array_1d<double, 3> dummy4 = ZeroVector(3);
        double distance = 0.0;
        int ContactType = -1;
        wall->ComputeConditionRelativeData(i, this, dummy1, distance, dummy2, dummy3, dummy4, ContactType);

        // Check if particle is close enough to wall
        if (distance > particle_radius + layer * particle_radius)
          continue;

        // Compute heat transfer coefficient
        // (Assumption: particle-Wall is treated as 2 equal size particles)          
        double a = (distance - 2*particle_radius) / particle_radius;
        double r_in;
        double r_out;

        if (distance > 2*particle_radius + min_dist)
          r_in = 0.0;
        else
          r_in = sqrt(1-pow(min_dist/particle_radius-a-1,2));

        if (a > sqrt(pow((particle_radius + (layer*particle_radius)) / particle_radius,2) - 1) - 1)
          r_out = sqrt(pow((particle_radius + (layer*particle_radius)) / particle_radius,2) - (a+1)*(a+1));
        else
          r_out = 1.0;

        double b = sqrt(1-r_out*r_out);
        double c = sqrt(1-r_in*r_in);

        double h = 2 * Globals::Pi * fluid_conductivity * particle_radius * ((a+1) * log(abs((b-a-1) / (a-c+1))) + b - c);

        // Compute heat flux
        mConductiveHeatFlux += h * temp_grad;
      }
    }

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeConvectiveHeatFlux(const ProcessInfo& r_process_info) {
    KRATOS_TRY
    
    // Get particle properties
    double radius                    = GetRadius();
    double surface_area              = 4 * Globals::Pi * radius*radius;
    double char_length               = 2 * radius;
    double particle_temp             = GetParticleTemperature();
    array_1d<double, 3> particle_vel = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);

    // Get surrounding fluid properties
    double              fluid_density      = r_process_info[FLUID_DENSITY];
    double              fluid_viscosity    = r_process_info[FLUID_VISCOSITY];
    double              fluid_conductivity = r_process_info[FLUID_THERMAL_CONDUCTIVITY];
    double              fluid_heatcapacity = r_process_info[FLUID_HEAT_CAPACITY];
    double              fluid_temp         = r_process_info[FLUID_TEMPERATURE];
    array_1d<double, 3> fluid_velocity     = r_process_info[FLUID_VELOCITY];

    // Relative velocity between particle and fluid
    array_1d<double, 3> rel_velocity_vec = fluid_velocity;

    for (unsigned int i = 0; i < rel_velocity_vec.size(); i++)
      rel_velocity_vec[i] -= particle_vel[i];

    double rel_velocity_mod = DEM_MODULUS_3(rel_velocity_vec);

    // Compute Prandtl and Reynolds numbers
    double Pr = fluid_viscosity * fluid_heatcapacity / fluid_conductivity;
    double Re = fluid_density * char_length * rel_velocity_mod / fluid_viscosity;

    // Compute Nusselt number according to selected correlation
    std::string model = r_process_info[CONVECTION_MODEL];
    double Nu = 0.0;

    if (model.compare("sphere_hanz_marshall") == 0) {
      Nu = 2 + 0.6 * pow(Re,1/2) * pow(Pr,1/3);
    }
    else if (model.compare("sphere_whitaker") == 0) {
      Nu = 2 + (0.4 * pow(Re,1/2) + 0.06 * pow(Re,2/3)) * pow(Pr,2/5);
    }

    // Compute thermal convection coefficient 
    double convection_coeff = Nu * fluid_conductivity / char_length;

    // Compute heat flux
    mConvectiveHeatFlux += convection_coeff * surface_area * (fluid_temp - particle_temp);

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeRadiativeHeatFlux(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    // Stefan-Boltzmann constant
    const double B = 5.670374419e-8;

    // Compute heat flux according to selected model
    std::string model = r_process_info[RADIATION_MODEL];

    if (model.compare("continuum_zhou") == 0) {
      // TODO: Needs the porosity
      mRadiativeHeatFlux += 0.0;
    }
    else if (model.compare("continuum_krause") == 0) {

      // Get particle properties
      double this_radius     = GetRadius();
      double this_area       = 4 * Globals::Pi * this_radius*this_radius;
      double this_emissivity = GetProperties()[EMISSIVITY];
      double this_temp       = GetParticleTemperature();

      // Initialize parameters
      double num = 0.0;
      double den = 0.0;

      // Loop over neighbour particles
      // (Assumption: walls not considered)
      for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
        if (mNeighbourElements[i] == NULL) continue;
        ThermalSphericParticle<TBaseElement>* neighbour_iterator = dynamic_cast<ThermalSphericParticle<TBaseElement>*>(mNeighbourElements[i]);
        
        // Compute direction and distance between centroids
        array_1d<double, 3> direction;
        direction[0] = GetGeometry()[0].Coordinates()[0] - neighbour_iterator->GetGeometry()[0].Coordinates()[0];
        direction[1] = GetGeometry()[0].Coordinates()[1] - neighbour_iterator->GetGeometry()[0].Coordinates()[1];
        direction[2] = GetGeometry()[0].Coordinates()[2] - neighbour_iterator->GetGeometry()[0].Coordinates()[2];

        double distance = DEM_MODULUS_3(direction);

        // Check if particles are close enough
        // (Assumption: radiation influence factor applied to the maximum radius)
        double other_radius = neighbour_iterator->GetRadius();
        if (distance > r_process_info[RADIATION_RADIUS_FACTOR] * std::max(this_radius, other_radius))
          continue;

        // Get particle properties
        double other_area       = 4 * Globals::Pi * other_radius*other_radius;
        double other_emissivity = neighbour_iterator->GetProperties()[EMISSIVITY];
        double other_temp       = neighbour_iterator->GetParticleTemperature();

        // Update parameters
        den += B * other_emissivity * other_area / 2;
        num += den * pow(other_temp,4);
      }

      // Averaged environment temperature
      double env_temp = pow(num/den,1/4);

      // Compute heat flux
      mRadiativeHeatFlux += this_emissivity * B * this_area * (pow(env_temp,4) - pow(this_temp,4));
    }

    KRATOS_CATCH("")
  }

  // Update methods

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::UpdateTemperature(const ProcessInfo& r_process_info) {
    // Condition to avoid issue when mSpecificHeat is equal zero
    if (mSpecificHeat > 0) {
      // Compute new temperature
      double dt              = r_process_info[DELTA_TIME];
      double thermal_inertia = GetMass() * mSpecificHeat;
      double temp_increment  = mTotalHeatFlux / thermal_inertia * dt;
      double temp_new        = GetParticleTemperature() + temp_increment;
      
      // Set new temperature
      SetParticleTemperature(temp_new);
    }
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::UpdateTemperatureDependentRadius(const ProcessInfo& r_process_info) {
    double this_temp      = GetParticleTemperature();
    double fluid_temp     = r_process_info[FLUID_TEMPERATURE];
    double relative_temp  = this_temp - fluid_temp; // temp in Kelvin
    double thermal_alpha  = GetProperties()[THERMAL_EXPANSION_COEFFICIENT];
    double updated_radius = GetRadius() * (1 + thermal_alpha * relative_temp);
    SetRadius(updated_radius);
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::UpdateNormalRelativeDisplacementAndVelocityDueToThermalExpansion(const ProcessInfo& r_process_info,
                                                                                                              double& thermalDeltDisp,
                                                                                                              double& thermalRelVel,
                                                                                                              ThermalSphericParticle<TBaseElement>* element2) {
    //thermalRelVel                      = 0;
    //double temperature                 = GetParticleTemperature();
    //double other_temperature           = element2->GetParticleTemperature();
    //double previous_temperature        = mPreviousTemperature;
    //double previous_other_temperature  = element2->mPreviousTemperature;
    //double thermal_alpha               = GetProperties()[THERMAL_EXPANSION_COEFFICIENT];
    //double updated_radius              = GetRadius();
    //double updated_other_radius        = element2->GetRadius();
    //double dt                          = r_process_info[DELTA_TIME];
    //double temperature_increment_elem1 = temperature - previous_temperature;
    //double temperature_increment_elem2 = other_temperature - previous_other_temperature;
    //thermalDeltDisp                    = updated_radius * thermal_alpha * temperature_increment_elem1;
    //thermalDeltDisp                    = thermalDeltDisp + updated_other_radius * thermal_alpha * temperature_increment_elem2;
    //thermalRelVel                      = updated_radius * thermal_alpha * temperature_increment_elem1 / dt;
    //thermalRelVel                      = thermalRelVel + updated_other_radius * thermal_alpha * temperature_increment_elem2 / dt;
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::RelativeDisplacementAndVelocityOfContactPointDueToOtherReasons(const ProcessInfo& r_process_info,
                                                                                                            double DeltDisp[3], //IN GLOBAL AXES
                                                                                                            double RelVel[3],   //IN GLOBAL AXES
                                                                                                            double OldLocalCoordSystem[3][3],
                                                                                                            double LocalCoordSystem[3][3],
                                                                                                            SphericParticle* neighbour_iterator) {
    //double thermalDeltDisp = 0;
    //double thermalRelVel   = 0;
    //ThermalSphericParticle<TBaseElement>* thermal_neighbour_iterator = dynamic_cast<ThermalSphericParticle<TBaseElement>*>(neighbour_iterator);
    //UpdateNormalRelativeDisplacementAndVelocityDueToThermalExpansion(r_process_info, thermalDeltDisp, thermalRelVel, thermal_neighbour_iterator);
    //double LocalRelVel[3] = {0.0};
    //GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, RelVel, LocalRelVel); //TODO: can we do this in global axes directly?
    ////LocalRelVel[2] -= thermalRelVel;
    //GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalRelVel, RelVel);
  }

  // Finalization methods

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::FinalizeSolutionStep(const ProcessInfo& r_process_info) {
    TBaseElement::FinalizeSolutionStep(r_process_info);
    UpdateTemperature(r_process_info);
    mPreviousTemperature = GetParticleTemperature();
    GetGeometry()[0].GetSolutionStepValue(HEATFLUX) = mTotalHeatFlux;
  }

  // Auxiliary computation methods

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeContactArea(const double rmin, double indentation, double& calculation_area) {
    calculation_area = Globals::Pi*rmin*rmin;
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeAddedSearchDistance(const ProcessInfo& r_process_info, double& added_search_distance) {
    KRATOS_TRY

    if (this->Is(DEMFlags::HAS_INDIRECT_CONDUCTION)){
      std::string model = r_process_info[INDIRECT_CONDUCTION_MODEL];
      if (model.compare("surrounding_layer") == 0) {
        double model_search_distance = GetRadius() * r_process_info[FLUID_LAYER_THICKNESS];
        added_search_distance = std::max(added_search_distance, model_search_distance);
      }
    }

    if (this->Is(DEMFlags::HAS_RADIATION)) {
      std::string model = r_process_info[RADIATION_MODEL];
      if (model.compare("sphere_hanz_marshall") == 0 ||
          model.compare("sphere_whitaker") == 0) {
        double model_search_distance = GetRadius() * (r_process_info[RADIATION_RADIUS_FACTOR] - 1);
        added_search_distance = std::max(added_search_distance, model_search_distance);
      }
    }

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::IntegralSurrLayer(const ProcessInfo& r_process_info, double d, double r1, double r2, double a, double b) {
    KRATOS_TRY

    // Initialization
    double fa  = EvalIntegrandSurrLayer(r_process_info, d, r1, r2, a);
    double fb  = EvalIntegrandSurrLayer(r_process_info, d, r1, r2, b);
    double fc  = EvalIntegrandSurrLayer(r_process_info, d, r1, r2, (a+b)/2);
    double tol = r_process_info[INTEGRAL_TOLERANCE];
    constexpr double eps = std::numeric_limits<double>::epsilon();
    if (tol < 10.0*eps)
      tol = 10.0*eps;

    // Solve integral recursively with adaptive Simpson quadrature
    return SolveIntegralSurrLayer(r_process_info, d, r1, r2, a, b, tol, fa, fb, fc);

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::SolveIntegralSurrLayer(const ProcessInfo& r_process_info, double d, double r1, double r2, double a, double b, double tol, double fa, double fb, double fc) {
    KRATOS_TRY

    // TODO: in order to catch possible erros that can occur in singularities,
    //       add a min value for subdivision size (to contain machine representable point) and a max number of function evaluation (+- 10000).

    double c  = (a + b) / 2;
    double fd = EvalIntegrandSurrLayer(r_process_info, d, r1, r2, (a+c)/2);
    double fe = EvalIntegrandSurrLayer(r_process_info, d, r1, r2, (c+b)/2);
    double I1 = (b-a)/6  * (fa + 4*fc + fb);
    double I2 = (b-a)/12 * (fa + 4*fd + 2*fc + 4*fe + fb);

    if (fabs(I2-I1) <= tol) {
      return I2 + (I2 - I1) / 15;
    }
    else { // sub-divide interval recursively
      double Ia = SolveIntegralSurrLayer(r_process_info, d, r1, r2, a, c, tol, fa ,fd, fc);
      double Ib = SolveIntegralSurrLayer(r_process_info, d, r1, r2, c, b, tol, fc, fe, fb);
      return Ia + Ib;
    }

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  double ThermalSphericParticle<TBaseElement>::EvalIntegrandSurrLayer(const ProcessInfo& r_process_info, double d, double r1, double r2, double r) {
    KRATOS_TRY
    return 2 * Globals::Pi * r / std::max(r_process_info[MIN_CONDUCTION_DISTANCE], d - sqrt(r1*r1-r*r) - sqrt(r2*r2-r*r));
    KRATOS_CATCH("")
  }

  // Get/Set methods

  template <class TBaseElement>
  const double& ThermalSphericParticle<TBaseElement>::GetParticleTemperature() {
    return GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::SetParticleTemperature(const double temperature) {
    GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE) = temperature;
  }

  //Explicit Instantiation
  template class ThermalSphericParticle<SphericParticle>;
  template class ThermalSphericParticle<SphericContinuumParticle>;

} // namespace Kratos
