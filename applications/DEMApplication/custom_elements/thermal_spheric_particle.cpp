//
//   Project Name:                     ThermalDEM $
//   Last Modified by:    $Author: Ferran Arrufat $
//   Date:                $Date:    February 2015 $
//   Revision:            $Revision:      1.0.0.0 $
//

// System includes
#include <string>
#include <iostream>

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

    // Sum up contributions
    mTotalHeatFlux = mConductiveHeatFlux + mConvectiveHeatFlux;
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

      // Get properties
      double other_conductivity = neighbour_iterator->mThermalConductivity;
      double other_heatcapacity = neighbour_iterator->mSpecificHeat;

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
        double this_density  = GetDensity();
        double this_mass     = GetMass();
        double this_poisson  = GetPoisson();
        double this_young    = GetYoung();
        double other_density = neighbour_iterator->GetDensity();
        double other_mass    = neighbour_iterator->GetMass();
        double other_poisson = neighbour_iterator->GetPoisson();
        double other_young   = neighbour_iterator->GetYoung();

        // Compute effective parameters
        double eff_radius = this_radius * other_radius / (this_radius + other_radius);
        double eff_mass   = this_mass * other_mass / (this_mass + other_mass);
        double eff_young  = 1 / ((1-this_poisson*this_poisson)/this_young + (1-other_poisson*other_poisson)/other_young);

        // Compute expected collision time
        double impact_normal_velocity = 0.01;  // TODO: save impact normal velocity
        double col_time = 0.0;                 // TODO: track collision time
        double expect_col_time = 2.87 * pow(eff_mass*eff_mass / (eff_radius*eff_young*eff_young*impact_normal_velocity),1/5);

        // Check if collision time is smaller than expected value, otherwise use batchelor_obrien model
        if (col_time < expect_col_time && impact_normal_velocity != 0.0) {
          // Compute max contact radius
          double contact_radius_max = pow(15*eff_radius*eff_mass*impact_normal_velocity*impact_normal_velocity / (16*eff_young),1/5);

          // Compute coefficients
          double a1 = this_density * mSpecificHeat;
          double a2 = other_density * other_heatcapacity;
          double b1 = a1 * mThermalConductivity;
          double b2 = a2 * other_conductivity;
          double c  = a1/a2;

          // Fourier number (taken as the aerage of both particles)
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
          // Compute average thermal conductivity
          double avg_conductivity = (this_radius + other_radius) / (this_radius/mThermalConductivity + other_radius/other_conductivity);

          // Compute contact area
          double contact_radius = sqrt(fabs(this_radius*this_radius - pow(((this_radius*this_radius - other_radius*other_radius + distance*distance) / (2*distance)),2)));
          double contact_area   = Globals::Pi*contact_radius*contact_radius;

          // Compute heat flux
          mConductiveHeatFlux += avg_conductivity * contact_area * temp_grad / distance;
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
      
      // Check if particle is in contact with wall (TODO)
      // Compute distance between particle centroid and wall
      //array_1d<double,3> direction = GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();
      //double distance = DEM_MODULUS_3(direction);
      double identation = 0.001;

      // Compute heat flux according to selected model
      if (model.compare("batchelor_obrien") == 0) {
        // Get particle properties
        double particle_radius = GetRadius();
        double particle_temp   = GetParticleTemperature();

        // Get wall properties
        double wall_conductivity = wall->GetProperties()[THERMAL_CONDUCTIVITY];

        // Get wall temperature as the average of its nodes
        double wall_temp = 0.0;
        double n_nodes = wall->GetGeometry().size();
        for (unsigned int i = 0; i < n_nodes; i++) {
          double node_temp = wall->GetGeometry()[i].FastGetSolutionStepValue(TEMPERATURE);
          wall_temp += node_temp;
        }
        wall_temp /= n_nodes;

        // Compute effective thermal conductivity
        double eff_conductivity = mThermalConductivity * wall_conductivity / (mThermalConductivity + wall_conductivity);

        // Compute contact radius
        double contact_radius = sqrt(identation * (2 * particle_radius - identation));

        // Compute heat flux
        mConductiveHeatFlux += 4 * eff_conductivity * contact_radius * (wall_temp - particle_temp);
      }
      else if (model.compare("thermal_pipe") == 0) {
        mConductiveHeatFlux += 0.0;
      }
      else if (model.compare("collisional") == 0) {
        mConductiveHeatFlux += 0.0;
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

      // Compute heat flux according to selected model
      if (model.compare("surrounding_layer") == 0) {
        // Compute heat flux
        mConductiveHeatFlux += 0.0;
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
      
      // Compute heat flux according to selected model
      if (model.compare("surrounding_layer") == 0) {
        // Compute heat flux
        mConductiveHeatFlux += 0.0;
      }
    }

    KRATOS_CATCH("")
  }

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeConvectiveHeatFlux(const ProcessInfo& r_process_info) {
    KRATOS_TRY
    
    // Get particle properties
    double radius                    = GetRadius();
    double surface_area              = 4 * Globals::Pi*radius*radius;
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

  // Auxiliary computation methods

  template <class TBaseElement>
  void ThermalSphericParticle<TBaseElement>::ComputeContactArea(const double rmin, double indentation, double& calculation_area) {
    calculation_area = Globals::Pi*rmin*rmin;
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
