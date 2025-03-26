#ifndef CALCULATE_GLOBAL_PHYSICAL_PROPERTIES_H
#define CALCULATE_GLOBAL_PHYSICAL_PROPERTIES_H

// /* External includes */

// System includes

// Project includes
#include "utilities/timer.h"
#include "custom_utilities/create_and_destroy.h"
#include "custom_elements/Particle_Contact_Element.h"
#include "includes/variables.h"

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

/* Project includes */
#include "includes/define.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{

class SphericElementGlobalPhysicsCalculator
    {
     public:

     typedef ModelPart::ElementsContainerType ElementsArrayType;

     KRATOS_CLASS_POINTER_DEFINITION(SphericElementGlobalPhysicsCalculator);

      /// Default constructor.

      SphericElementGlobalPhysicsCalculator(ModelPart& r_model_part)
      {
          mInitialCenterOfMassAndMass = CalculateCenterOfMass(r_model_part);
          mInitialMass                = CalculateTotalMass(r_model_part);
      }

      /// Destructor.

      virtual ~SphericElementGlobalPhysicsCalculator(){}


      //***************************************************************************************************************
      //***************************************************************************************************************

      double CalculateTotalVolume(ModelPart& r_model_part)
      {
          OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);
          double added_volume = 0.0;

          #pragma omp parallel for reduction(+ : added_volume)
          for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){

              for (ElementsArrayType::iterator it = GetElementPartitionBegin(r_model_part, k); it != GetElementPartitionEnd(r_model_part, k); ++it){
                  if (it->GetGeometry()[0].Is(BLOCKED)) { // we exclude blocked elements from the volume calculation (e.g., inlet injectors)
                      continue;
                  }
                  if ((it)->IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)) {
                      SphericParticle& r_spheric_particle = dynamic_cast<Kratos::SphericParticle&> (*it);
                      const double particle_radius = r_spheric_particle.GetRadius();
                      added_volume += 4.0 / 3.0 * Globals::Pi * particle_radius * particle_radius * particle_radius;
                  }
            }
        }

        return added_volume;
      }

      double CalculatePorosityWithinSphere(ModelPart& r_model_part, const double radius, const array_1d<double, 3>& center)
      {
          OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);
          double sphere_volume_inside_range = 0.0;
          const double total_volume = 4.0 / 3.0 * Globals::Pi * radius * radius * radius;

          #pragma omp parallel for reduction(+ : sphere_volume_inside_range)
          for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){

              for (ElementsArrayType::iterator it = GetElementPartitionBegin(r_model_part, k); it != GetElementPartitionEnd(r_model_part, k); ++it){
                  
                  if ((it)->IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)) {
                      SphericParticle& r_spheric_particle = dynamic_cast<Kratos::SphericParticle&> (*it);
                      const double particle_radius = r_spheric_particle.GetRadius();
                      const array_1d<double, 3>& particle_coordinates = r_spheric_particle.GetGeometry()[0].Coordinates();
                      const double distance = std::sqrt(std::pow(particle_coordinates[0] - center[0], 2) + std::pow(particle_coordinates[1] - center[1], 2) + std::pow(particle_coordinates[2] - center[2], 2));
                      if (distance < radius - particle_radius) {
                        sphere_volume_inside_range += 4.0 / 3.0 * Globals::Pi * particle_radius * particle_radius * particle_radius;
                      } else if (distance <= radius + particle_radius) {
                        const double other_part_d = radius - (radius * radius + distance * distance - particle_radius * particle_radius) / (distance * 2);
                        const double my_part_d = particle_radius - (particle_radius * particle_radius + distance * distance - radius * radius) / (distance * 2);
                        const double cross_volume = Globals::Pi * other_part_d * other_part_d * (radius - 1.0 / 3.0 * other_part_d) + Globals::Pi * my_part_d * my_part_d * (particle_radius - 1.0 / 3.0 * my_part_d);
                        sphere_volume_inside_range += cross_volume;
                      }
                  }
              }
          }
          return 1.0 - sphere_volume_inside_range / total_volume;
      }

      //***************************************************************************************************************
      //***************************************************************************************************************
      // Returns the minimum value of a double variable in the model part.

    double CalculateMaxNodalVariable(ModelPart& r_model_part, const Variable<double>& r_variable) {
        ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();

        KRATOS_ERROR_IF(pElements.size() == 0) << "Cannot compute maximum of the required nodal variable. Empty model part. Could not compute the maximum of the required variable " << r_variable << std::endl;

        ElementsArrayType::iterator it_begin = pElements.ptr_begin();

        KRATOS_ERROR_IF_NOT(it_begin->GetGeometry()[0].SolutionStepsDataHas(r_variable)) << "Cannot compute maximum of the required nodal variable. Missing nodal variable " << r_variable << std::endl;

        std::vector<double> max_values;
        double max_val = - std::numeric_limits<double>::max();
        max_values.resize(ParallelUtilities::GetNumThreads());

        for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){
            max_values[k] = max_val;
        }

        OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), pElements.size(), mElementsPartition);

        unsigned int elem_counter;

        #pragma omp parallel for private(elem_counter)
        for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){
            elem_counter = mElementsPartition[k];

            for (ElementsArrayType::iterator it = GetElementPartitionBegin(r_model_part, k); it != GetElementPartitionEnd(r_model_part, k); ++it){
                max_values[k] = std::max(max_values[k], (it)->GetGeometry()[0].FastGetSolutionStepValue(r_variable));
                elem_counter++;
            }
        }

        // getting the maximum between threads:
        for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){
            max_val = std::max(max_val, max_values[k]);
        }

        return max_val;
    }

      //***************************************************************************************************************
      //***************************************************************************************************************

    double CalculateMinNodalVariable(ModelPart& r_model_part, const Variable<double>& r_variable) {
        ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();

        KRATOS_ERROR_IF(pElements.size() == 0) << "Cannot compute minimum of the required nodal variable. Empty model part. Could not compute the maximum of the required variable " << r_variable << std::endl;

        ElementsArrayType::iterator it_begin = pElements.ptr_begin();

        KRATOS_ERROR_IF_NOT(it_begin->GetGeometry()[0].SolutionStepsDataHas(r_variable)) << "Cannot compute minimum of the required nodal variable. Missing variable " << r_variable << std::endl;

        std::vector<double> min_values;
        double min_val = std::numeric_limits<double>::max();
        min_values.resize(ParallelUtilities::GetNumThreads());

        for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){
            min_values[k] = min_val;
        }

        OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), pElements.size(), mElementsPartition);

        unsigned int elem_counter;

        #pragma omp parallel for private(elem_counter)
        for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){
            elem_counter = mElementsPartition[k];

            for (ElementsArrayType::iterator it = GetElementPartitionBegin(r_model_part, k); it != GetElementPartitionEnd(r_model_part, k); ++it){
                min_values[k] = std::min(min_values[k], (it)->GetGeometry()[0].FastGetSolutionStepValue(r_variable));
                elem_counter++;
            }
        }

        // getting the minimum between threads:
        for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){
            min_val = std::min(min_val, min_values[k]);
        }

        return min_val;
    }

      //***************************************************************************************************************
      //***************************************************************************************************************

      double CalculateD50(ModelPart& r_model_part)
      {
          const unsigned int size = r_model_part.GetCommunicator().LocalMesh().Elements().size();
          OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), size, mElementsPartition);

          std::vector<double> radii;
          radii.resize(size);
          unsigned int particle_counter = 0;

          #pragma omp parallel for private(particle_counter)
          for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){
              particle_counter = mElementsPartition[k];

              for (ElementsArrayType::iterator it = GetElementPartitionBegin(r_model_part, k); it != GetElementPartitionEnd(r_model_part, k); ++it){
                  SphericParticle& r_spheric_particle = dynamic_cast<Kratos::SphericParticle&> (*it);
                  radii[particle_counter] = r_spheric_particle.GetRadius();
                  particle_counter++;
                }

            }
          if (particle_counter) {
            std::sort(radii.begin(), radii.end());
            int half   = div(size, 2).quot;
            bool even  = (size%2 == 0);
            double d50 = even ? 2 * radii[half] : radii[half] + radii[half + 1];

            return d50;
          }
          else {
            return 0.00;
          }
      }

      //***************************************************************************************************************
      //***************************************************************************************************************

      double CalculateTotalMass(ModelPart& r_model_part)
      {
          OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(),r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

          double added_mass = 0.0;

          #pragma omp parallel for reduction(+ : added_mass)
          for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){

              for (ElementsArrayType::iterator it = GetElementPartitionBegin(r_model_part, k); it != GetElementPartitionEnd(r_model_part, k); ++it){
                  if ((it)->IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)) {
                      double particle_mass = (it)->GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS);
                      added_mass += particle_mass;
                  }
              }
        }

        return added_mass;
      }

      //***************************************************************************************************************
      //***************************************************************************************************************

      array_1d<double, 3> CalculateCenterOfMass(ModelPart& r_model_part)
      {
          OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

          const double total_mass_inv         = 1 / CalculateTotalMass(r_model_part);
          double cm_x = 0.0;
          double cm_y = 0.0;
          double cm_z = 0.0;

          #pragma omp parallel for reduction(+ : cm_x, cm_y, cm_z)
          for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){

              for (ElementsArrayType::iterator it = GetElementPartitionBegin(r_model_part, k); it != GetElementPartitionEnd(r_model_part, k); ++it){
                  if ((it)->IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)) {
                      double particle_mass = (it)->GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS);
                      cm_x += particle_mass * (it)->GetGeometry()[0].Coordinates()[0];
                      cm_y += particle_mass * (it)->GetGeometry()[0].Coordinates()[1];
                      cm_z += particle_mass * (it)->GetGeometry()[0].Coordinates()[2];
                    }
              }
          }

          array_1d<double, 3> center_of_mass;
          center_of_mass[0] = total_mass_inv * cm_x;
          center_of_mass[1] = total_mass_inv * cm_y;
          center_of_mass[2] = total_mass_inv * cm_z;

          return center_of_mass;
      }

      //***************************************************************************************************************
      //***************************************************************************************************************

      double CalculateGravitationalPotentialEnergy(ModelPart& r_model_part, const array_1d<double, 3> reference_point)
      {
          double gravitational_energy;
          const double total_mass                               = CalculateTotalMass(r_model_part);
          if (total_mass == 0)  gravitational_energy = 0.0;
          else {
              const array_1d<double, 3>& gravity                    = r_model_part.GetProcessInfo()[GRAVITY];
              const array_1d<double, 3> center_of_mass              = CalculateCenterOfMass(r_model_part);
              const array_1d<double, 3> center_of_mass_to_reference = reference_point - center_of_mass;
              gravitational_energy = total_mass * (center_of_mass_to_reference[0] * gravity[0] + center_of_mass_to_reference[1] * gravity[1] + center_of_mass_to_reference[2] * gravity[2]);
          }
          return gravitational_energy;
      }

      //***************************************************************************************************************
      //***************************************************************************************************************

      double CalculateTranslationalKinematicEnergy(ModelPart& r_model_part)
      {
          OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

          double kinematic_energy = 0.0;

          #pragma omp parallel for reduction(+ : kinematic_energy)
          for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){

              for (ElementsArrayType::iterator it = GetElementPartitionBegin(r_model_part, k); it != GetElementPartitionEnd(r_model_part, k); ++it){
                  if ((it)->IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)) {
                      double particle_translational_kinematic_energy = 0.0;

                      (it)->Calculate(PARTICLE_TRANSLATIONAL_KINEMATIC_ENERGY, particle_translational_kinematic_energy, r_model_part.GetProcessInfo());

                      kinematic_energy += particle_translational_kinematic_energy;
                  }
              }

          }

          return kinematic_energy;
      }

      //***************************************************************************************************************
      //***************************************************************************************************************

      double CalculateRotationalKinematicEnergy(ModelPart& r_model_part)
      {
          OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

          double rotational_kinematic_energy = 0.0;

          #pragma omp parallel for reduction(+ : rotational_kinematic_energy)
          for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){

              for (ElementsArrayType::iterator it = GetElementPartitionBegin(r_model_part, k); it != GetElementPartitionEnd(r_model_part, k); ++it){
                  if ((it)->IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)) {
                      double particle_rotational_kinematic_energy = 0.0;

                      (it)->Calculate(PARTICLE_ROTATIONAL_KINEMATIC_ENERGY, particle_rotational_kinematic_energy, r_model_part.GetProcessInfo());

                      rotational_kinematic_energy += particle_rotational_kinematic_energy;
                  }
              }

          }

          return rotational_kinematic_energy;
      }

      //***************************************************************************************************************
      //***************************************************************************************************************

      double CalculateElasticEnergy(ModelPart& r_model_part)
      {
          OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

          double elastic_energy = 0.0;

          #pragma omp parallel for reduction(+ : elastic_energy)
          for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){

              for (ElementsArrayType::iterator it = GetElementPartitionBegin(r_model_part, k); it != GetElementPartitionEnd(r_model_part, k); ++it){
                  if ((it)->IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)) {
                      double particle_elastic_energy = 0.0;

                      (it)->Calculate(PARTICLE_ELASTIC_ENERGY, particle_elastic_energy, r_model_part.GetProcessInfo());

                      elastic_energy += particle_elastic_energy;
                  }
              }

          }

          return elastic_energy;
      }

      //***************************************************************************************************************
      //***************************************************************************************************************

      double CalculateInelasticFrictionalEnergy(ModelPart& r_model_part)
      {
          OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

          double frictional_energy = 0.0;

          #pragma omp parallel for reduction(+ : frictional_energy)
          for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){

              for (ElementsArrayType::iterator it = GetElementPartitionBegin(r_model_part, k); it != GetElementPartitionEnd(r_model_part, k); ++it){
                  if ((it)->IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)) {
                      double particle_frictional_energy = 0.0;

                      (it)->Calculate(PARTICLE_INELASTIC_FRICTIONAL_ENERGY, particle_frictional_energy, r_model_part.GetProcessInfo());

                      frictional_energy += particle_frictional_energy;
                  }
              }

          }

          return frictional_energy;
      }

      double CalculateInelasticViscodampingEnergy(ModelPart& r_model_part)
      {
          OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

          double viscodamping_energy = 0.0;

          #pragma omp parallel for reduction(+ : viscodamping_energy)
          for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){

              for (ElementsArrayType::iterator it = GetElementPartitionBegin(r_model_part, k); it != GetElementPartitionEnd(r_model_part, k); ++it){
                  if ((it)->IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)) {
                      double particle_viscodamping_energy = 0.0;

                      (it)->Calculate(PARTICLE_INELASTIC_VISCODAMPING_ENERGY, particle_viscodamping_energy, r_model_part.GetProcessInfo());

                      viscodamping_energy += particle_viscodamping_energy;
                  }
              }

          }

          return viscodamping_energy;
      }

      double CalculateInelasticRollingResistanceEnergy(ModelPart& r_model_part)
      {
          OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

          double rollingresistance_energy = 0.0;

          #pragma omp parallel for reduction(+ : rollingresistance_energy)
          for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){

              for (ElementsArrayType::iterator it = GetElementPartitionBegin(r_model_part, k); it != GetElementPartitionEnd(r_model_part, k); ++it){
                  if ((it)->IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)) {
                      double particle_rollingresistance_energy = 0.0;

                      (it)->Calculate(PARTICLE_INELASTIC_ROLLING_RESISTANCE_ENERGY, particle_rollingresistance_energy, r_model_part.GetProcessInfo());

                      rollingresistance_energy += particle_rollingresistance_energy;
                  }
              }

          }

          return rollingresistance_energy;
      }

      double CalculateParticleNumberTimesMaxNormalBallToBallForceTimesRadius(ModelPart& r_model_part)
      {

        double global_max_normal_force_times_radius = 0.0;

        ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();

        #pragma omp parallel
        {
            double local_max = 0.0;

            #pragma omp for
            for (int k = 0; k < (int)pElements.size(); k++) {

                ElementsArrayType::iterator it = pElements.ptr_begin() + k;

                if ((it)->IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)) {
                    double particle_normal_ball_to_ball_force_times_radius = 0.0;

                    (it)->Calculate(PARTICLE_MAX_NORMAL_BALL_TO_BALL_FORCE_TIMES_RADIUS, particle_normal_ball_to_ball_force_times_radius, r_model_part.GetProcessInfo());

                    if (local_max < particle_normal_ball_to_ball_force_times_radius) {
                        local_max = particle_normal_ball_to_ball_force_times_radius;
                    }
                }
            }

            #pragma omp critical
            {
                if (global_max_normal_force_times_radius < local_max) {
                    global_max_normal_force_times_radius = local_max;
                }
            }
        }

          return global_max_normal_force_times_radius * r_model_part.GetCommunicator().LocalMesh().Elements().size();
      }

      //***************************************************************************************************************
      //***************************************************************************************************************

      array_1d<double, 3> CalculateTotalMomentum(ModelPart& r_model_part)
      {
          OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

          double m_x = 0.0;
          double m_y = 0.0;
          double m_z = 0.0;

          #pragma omp parallel for reduction(+ : m_x, m_y, m_z)
          for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){

              for (ElementsArrayType::iterator it = GetElementPartitionBegin(r_model_part, k); it != GetElementPartitionEnd(r_model_part, k); ++it){
                  if ((it)->IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)) {
                      array_1d<double, 3> particle_momentum;
                      (it)->Calculate(MOMENTUM, particle_momentum, r_model_part.GetProcessInfo());
                      m_x += particle_momentum[0];
                      m_y += particle_momentum[1];
                      m_z += particle_momentum[2];
                  }
              }

          }

          array_1d<double, 3> momentum;
          momentum[0] = m_x;
          momentum[1] = m_y;
          momentum[2] = m_z;

          return momentum;
      }

      //***************************************************************************************************************
      //***************************************************************************************************************

      array_1d<double, 3> CalulateTotalAngularMomentum(ModelPart& r_model_part)
      {
          OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

          const array_1d<double, 3> center_of_mass   = CalculateCenterOfMass(r_model_part);
          double am_x = 0.0;
          double am_y = 0.0;
          double am_z = 0.0;

          #pragma omp parallel for reduction(+ : am_x, am_y, am_z)
          for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){

              for (ElementsArrayType::iterator it = GetElementPartitionBegin(r_model_part, k); it != GetElementPartitionEnd(r_model_part, k); ++it){
                  if ((it)->IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)) {
                      array_1d<double, 3> particle_momentum;
                      array_1d<double, 3> particle_local_angular_momentum;
                      array_1d<double, 3> center_of_mass_to_particle = (it)->GetGeometry()[0].Coordinates() - center_of_mass;

                      (it)->Calculate(MOMENTUM, particle_momentum, r_model_part.GetProcessInfo());
                      (it)->Calculate(ANGULAR_MOMENTUM, particle_local_angular_momentum, r_model_part.GetProcessInfo());

                      array_1d<double, 3> aux;
                      Kratos::MathUtils<double>::CrossProduct(aux, particle_momentum, center_of_mass_to_particle);

                      am_x += particle_local_angular_momentum[0] + aux[0];
                      am_y += particle_local_angular_momentum[1] + aux[1];
                      am_z += particle_local_angular_momentum[2] + aux[2];
                  }
              }
          }

          array_1d<double, 3> angular_momentum;
          angular_momentum[0] = am_x;
          angular_momentum[1] = am_y;
          angular_momentum[2] = am_z;

          return angular_momentum;
      }

      //***************************************************************************************************************
      //***************************************************************************************************************
      // Check by how much Newton's Third Law is violated
      array_1d<double, 3> CalculateSumOfInternalForces(ModelPart& r_model_part)
      {
            OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(),r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);
            double sum_of_contact_forces_x = 0.0;
            double sum_of_contact_forces_y = 0.0;
            double sum_of_contact_forces_z = 0.0;

            #pragma omp parallel for reduction(+ : sum_of_contact_forces_x, sum_of_contact_forces_y, sum_of_contact_forces_z)
            for (int k = 0; k < ParallelUtilities::GetNumThreads(); ++k){

                for (ElementsArrayType::iterator it = GetElementPartitionBegin(r_model_part, k); it != GetElementPartitionEnd(r_model_part, k); ++it){
                    if ((it)->IsNot(DEMFlags::BELONGS_TO_A_CLUSTER)){
                        const array_1d<double, 3>& contact_force = (it)->GetGeometry()[0].FastGetSolutionStepValue(CONTACT_FORCES);
                        sum_of_contact_forces_x += contact_force[0];
                        sum_of_contact_forces_y += contact_force[1];
                        sum_of_contact_forces_z += contact_force[2];
                    }
                }
            }

            array_1d<double, 3> sum_of_contact_forces;
            sum_of_contact_forces[0] = sum_of_contact_forces_x;
            sum_of_contact_forces[1] = sum_of_contact_forces_y;
            sum_of_contact_forces[2] = sum_of_contact_forces_z;
            return sum_of_contact_forces;
      }

      double CalculateSumOfParticlesWithinSphere(ModelPart& sphere_model_part, const double radius, const array_1d<double, 3>& center)
      {
        OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), sphere_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);
        double total_particle_number = 0;

        #pragma omp parallel for reduction(+ : total_particle_number)
        for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){

            for (ElementsArrayType::iterator it = GetElementPartitionBegin(sphere_model_part, k); it != GetElementPartitionEnd(sphere_model_part, k); ++it){

                double r = (it)->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
                double x = (it)->GetGeometry()[0].X();
                double y = (it)->GetGeometry()[0].Y();
                double z = (it)->GetGeometry()[0].Z();

                double center_to_sphere_distance = std::sqrt(std::pow(x - center[0], 2) + std::pow(y - center[1], 2) + std::pow(z - center[2], 2));

                if (center_to_sphere_distance < (radius - r)) {
                    total_particle_number += 1;
                }
            }
        }
        return total_particle_number;
      }

      //***************************************************************************************************************
      //***************************************************************************************************************
        ///@}
        ///@name Access
        ///@{

        array_1d<double, 3> GetInitialCenterOfMass()
        {
            return mInitialCenterOfMassAndMass;
        }

        ///@}
        ///@name Inquiry
        ///@{


        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a stemplate<class T, std::size_t dim> tring.

        virtual std::string Info() const
        {
            return "";
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const
        {
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const
        {
        }


        ///@}
        ///@name Friends
        ///@{

        std::vector<unsigned int>& GetElementPartition()
        {
          return (mElementsPartition);
        }

        ElementsArrayType::iterator GetElementPartitionBegin(ModelPart& r_model_part, unsigned int k)
        {
          ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();
          return (pElements.ptr_begin() + mElementsPartition[k]);
        }

        ElementsArrayType::iterator GetElementPartitionEnd(ModelPart& r_model_part, unsigned int k)
        {
          ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();
          return (pElements.ptr_begin() + mElementsPartition[k + 1]);
        }

        ///@}

    protected:
        ///@name Protected static Member r_variables
        ///@{


        ///@}
        ///@name Protected member r_variables
        ///@{ template<class T, std::size_t dim>


        ///@}
        ///@name Protected Operators
        ///@{


        ///@}
        ///@name Protected Operations
        ///@{


        ///@}
        ///@name Protected  Access
        ///@{
        std::vector<unsigned int> mElementsPartition;

        ///@}
        ///@name Protected Inquiry
        ///@{


        ///@}
        ///@name Protected LifeCycle
        ///@{


        ///@}

    private:


        ///@name Static Member r_variables
        ///@{


        ///@}
        ///@name Member r_variables
        ///@{

        array_1d<double, 3> mInitialCenterOfMassAndMass;
        double mInitialMass;


        ///@}
        ///@name Private Operators
        ///@{

        ///@}
        ///@name Private Operations
        ///@{


        ///@}
        ///@name Private  Access
        ///@{


        ///@}
        ///@name Private Inquiry
        ///@{


        ///@}
        ///@name Un accessible methods
        ///@{

        /// Assignment operator.
        SphericElementGlobalPhysicsCalculator & operator=(SphericElementGlobalPhysicsCalculator const& rOther);


        ///@}

    }; // Class SphericElementGlobalPhysicsCalculator

class ContactElementGlobalPhysicsCalculator
    {
     public:

     typedef ModelPart::ElementsContainerType ElementsArrayType;

     KRATOS_CLASS_POINTER_DEFINITION(ContactElementGlobalPhysicsCalculator);

    /// Default constructor.

    ContactElementGlobalPhysicsCalculator(){}

    /// Destructor.

    virtual ~ContactElementGlobalPhysicsCalculator(){}


    //***************************************************************************************************************
    //***************************************************************************************************************

    std::vector<std::vector<double>> CalculateTotalStressTensor(ModelPart& r_model_part, double Lx, double Ly, double Lz)
    {
        OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);
        std::vector<std::vector<double>> measured_total_stress_tensor(3, std::vector<double>(3));

        double s_00, s_01, s_02, s_10, s_11, s_12, s_20, s_21, s_22;
        s_00 = s_01 = s_02 = s_10 = s_11 = s_12 = s_20 = s_21 = s_22 = 0.0;

        #pragma omp parallel for reduction(+ : s_00, s_01, s_02, s_10, s_11, s_12, s_20, s_21, s_22)
        for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){

            for (ElementsArrayType::iterator it = GetElementPartitionBegin(r_model_part, k); it != GetElementPartitionEnd(r_model_part, k); ++it){

                double x_0 = (it)->GetGeometry()[0].X();
                double x_1 = (it)->GetGeometry()[1].X();
                double y_0 = (it)->GetGeometry()[0].Y();
                double y_1 = (it)->GetGeometry()[1].Y();
                double z_0 = (it)->GetGeometry()[0].Z();
                double z_1 = (it)->GetGeometry()[1].Z();

                double dx, dy, dz;
                
                /*
                dx = x_0 - x_1;
                if (dx > 0.5 * Lx){
                    dx -= Lx;
                }
                else if (dx < -0.5 * Lx){
                    dx += Lx;
                }
                
                dy = y_0 - y_1;
                if (dy > 0.5 * Ly){
                    dy -= Ly;
                }
                else if (dy < -0.5 * Ly){
                    dy += Ly;
                }

                dz = z_0 - z_1;
                if (dz > 0.5 * Lz){
                    dz -= Lz;
                }
                else if (dz < -0.5 * Lz){
                    dz += Lz;
                }*/

                dx = x_0 - x_1;
                dy = y_0 - y_1;
                dz = z_0 - z_1;

                //only consider the inner contacts
                if (std::abs(dx) > 0.5 * Lx || std::abs(dy) > 0.5 * Ly || std::abs(dz) > 0.5 * Lz){
                    continue;
                }

                const array_1d<double, 3>& contact_force = (it)->GetValue(GLOBAL_CONTACT_FORCE);

                s_00 += contact_force[0] * dx;
                s_01 += contact_force[0] * dy;
                s_02 += contact_force[0] * dz;
                s_10 += contact_force[1] * dx;
                s_11 += contact_force[1] * dy;
                s_12 += contact_force[1] * dz;
                s_20 += contact_force[2] * dx;
                s_21 += contact_force[2] * dy;
                s_22 += contact_force[2] * dz;
            }
        }

        measured_total_stress_tensor[0][0] = s_00;
        measured_total_stress_tensor[0][1] = s_01;
        measured_total_stress_tensor[0][2] = s_02;
        measured_total_stress_tensor[1][0] = s_10;
        measured_total_stress_tensor[1][1] = s_11;
        measured_total_stress_tensor[1][2] = s_12;
        measured_total_stress_tensor[2][0] = s_20;
        measured_total_stress_tensor[2][1] = s_21;
        measured_total_stress_tensor[2][2] = s_22;

        return measured_total_stress_tensor;
    }

    std::vector<std::vector<double>> CalculateTotalStressTensorWithinSphere(ModelPart& contact_model_part, const double radius, const array_1d<double, 3>& center)
    {
        OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), contact_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);
        std::vector<std::vector<double>> measured_stress_tensor(3, std::vector<double>(3));

        double s_00, s_01, s_02, s_10, s_11, s_12, s_20, s_21, s_22;
        s_00 = s_01 = s_02 = s_10 = s_11 = s_12 = s_20 = s_21 = s_22 = 0.0;

        #pragma omp parallel for reduction(+ : s_00, s_01, s_02, s_10, s_11, s_12, s_20, s_21, s_22)
        for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){

            for (ElementsArrayType::iterator it = GetElementPartitionBegin(contact_model_part, k); it != GetElementPartitionEnd(contact_model_part, k); ++it){

                double x_0 = (it)->GetGeometry()[0].X();
                double x_1 = (it)->GetGeometry()[1].X();
                double y_0 = (it)->GetGeometry()[0].Y();
                double y_1 = (it)->GetGeometry()[1].Y();
                double z_0 = (it)->GetGeometry()[0].Z();
                double z_1 = (it)->GetGeometry()[1].Z();
                double r_0 = (it)->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
                double r_1 = (it)->GetGeometry()[1].FastGetSolutionStepValue(RADIUS);
                double r   = 0.5 * (r_0 + r_1);

                double center_to_sphere_distance_0 = std::sqrt(std::pow(x_0 - center[0], 2) + std::pow(y_0 - center[1], 2) + std::pow(z_0 - center[2], 2));
                double center_to_sphere_distance_1 = std::sqrt(std::pow(x_1 - center[0], 2) + std::pow(y_1 - center[1], 2) + std::pow(z_1 - center[2], 2));

                if (center_to_sphere_distance_0 < (radius - r) || center_to_sphere_distance_1 < (radius - r)) {
                    const array_1d<double, 3>& contact_force = (it)->GetValue(GLOBAL_CONTACT_FORCE);
                    double contact_force_vector[3] = {contact_force[0], contact_force[1], contact_force[2]};
                    double vector_l[3] = {x_0 - x_1, y_0 - y_1, z_0 - z_1};
                    double tensor[3][3];

                    for (int i = 0; i < 3; i++){
                        for (int j = 0; j < 3; j++){
                            tensor[i][j] = contact_force_vector[i] * vector_l[j];
                        }
                    }

                    s_00 += tensor[0][0];
                    s_01 += tensor[0][1];
                    s_02 += tensor[0][2];
                    s_10 += tensor[1][0];
                    s_11 += tensor[1][1];
                    s_12 += tensor[1][2];
                    s_20 += tensor[2][0];
                    s_21 += tensor[2][1];
                    s_22 += tensor[2][2];
                }
            }
        }

        measured_stress_tensor[0][0] = s_00;
        measured_stress_tensor[0][1] = s_01;
        measured_stress_tensor[0][2] = s_02;
        measured_stress_tensor[1][0] = s_10;
        measured_stress_tensor[1][1] = s_11;
        measured_stress_tensor[1][2] = s_12;
        measured_stress_tensor[2][0] = s_20;
        measured_stress_tensor[2][1] = s_21;
        measured_stress_tensor[2][2] = s_22;

        return measured_stress_tensor;

    }

    std::vector<std::vector<double>> CalculateFabricTensorWithinSphere(ModelPart& contact_model_part, const double radius, const array_1d<double, 3>& center)
    {
        OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), contact_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);
        std::vector<std::vector<double>> measured_fabric_tensor(3, std::vector<double>(3));

        double total_tensor[3][3];
        total_tensor[0][0] = total_tensor[0][1] = total_tensor[0][2] = total_tensor[1][0] = total_tensor[1][1] = total_tensor[1][2] = total_tensor[2][0] = total_tensor[2][1] = total_tensor[2][2] = 0.0;
        double total_contact_number = 0;

        double s_00, s_01, s_02, s_10, s_11, s_12, s_20, s_21, s_22;
        s_00 = s_01 = s_02 = s_10 = s_11 = s_12 = s_20 = s_21 = s_22 = 0.0;

        #pragma omp parallel for reduction(+ : s_00, s_01, s_02, s_10, s_11, s_12, s_20, s_21, s_22, total_contact_number)
        for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){

            for (ElementsArrayType::iterator it = GetElementPartitionBegin(contact_model_part, k); it != GetElementPartitionEnd(contact_model_part, k); ++it){

                double x_0 = (it)->GetGeometry()[0].X();
                double x_1 = (it)->GetGeometry()[1].X();
                double y_0 = (it)->GetGeometry()[0].Y();
                double y_1 = (it)->GetGeometry()[1].Y();
                double z_0 = (it)->GetGeometry()[0].Z();
                double z_1 = (it)->GetGeometry()[1].Z();
                double r_0 = (it)->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
                double r_1 = (it)->GetGeometry()[1].FastGetSolutionStepValue(RADIUS);

                double center_to_sphere_distance_0 = std::sqrt(std::pow(x_0 - center[0], 2) + std::pow(y_0 - center[1], 2) + std::pow(z_0 - center[2], 2));
                double center_to_sphere_distance_1 = std::sqrt(std::pow(x_1 - center[0], 2) + std::pow(y_1 - center[1], 2) + std::pow(z_1 - center[2], 2));

                if (center_to_sphere_distance_0 < (radius - r_0) || center_to_sphere_distance_1 < (radius - r_1)) {
                    double vector1[3] = {x_1 - x_0, y_1 - y_0, z_1 - z_0};
                    double v1_norm = std::sqrt(std::pow(vector1[0], 2) + std::pow(vector1[1], 2) + std::pow(vector1[2], 2));
                    double vector1_unit[3] = {vector1[0] / v1_norm, vector1[1] / v1_norm, vector1[2] / v1_norm};
                    double tensor[3][3];

                    for (int i = 0; i < 3; i++){
                        for (int j = 0; j < 3; j++){
                            tensor[i][j] = vector1_unit[i] * vector1_unit[j];
                        }
                    }

                    s_00 += tensor[0][0];
                    s_01 += tensor[0][1];
                    s_02 += tensor[0][2];
                    s_10 += tensor[1][0];
                    s_11 += tensor[1][1];
                    s_12 += tensor[1][2];
                    s_20 += tensor[2][0];
                    s_21 += tensor[2][1];
                    s_22 += tensor[2][2];
                    total_contact_number += 1;
                }
            }
        }

        if (total_contact_number) {
            measured_fabric_tensor[0][0] = s_00 / total_contact_number;
            measured_fabric_tensor[0][1] = s_01 / total_contact_number;
            measured_fabric_tensor[0][2] = s_02 / total_contact_number;
            measured_fabric_tensor[1][0] = s_10 / total_contact_number;
            measured_fabric_tensor[1][1] = s_11 / total_contact_number;
            measured_fabric_tensor[1][2] = s_12 / total_contact_number;
            measured_fabric_tensor[2][0] = s_20 / total_contact_number;
            measured_fabric_tensor[2][1] = s_21 / total_contact_number;
            measured_fabric_tensor[2][2] = s_22 / total_contact_number;
        }

        return measured_fabric_tensor;
    }

    double CalculateAveragedCoordinationNumberWithinSphere(ModelPart& sphere_model_part, ModelPart& contact_model_part, const double radius, const array_1d<double, 3>& center)
    {
        OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), contact_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);
        double total_particle_number = 0.0;
        double total_contact_number = 0.0;

        #pragma omp parallel for reduction(+ : total_contact_number)
        for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){

            for (ElementsArrayType::iterator it = GetElementPartitionBegin(contact_model_part, k); it != GetElementPartitionEnd(contact_model_part, k); ++it){

                double x_0 = (it)->GetGeometry()[0].X();
                double x_1 = (it)->GetGeometry()[1].X();
                double y_0 = (it)->GetGeometry()[0].Y();
                double y_1 = (it)->GetGeometry()[1].Y();
                double z_0 = (it)->GetGeometry()[0].Z();
                double z_1 = (it)->GetGeometry()[1].Z();
                double r_0 = (it)->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
                double r_1 = (it)->GetGeometry()[1].FastGetSolutionStepValue(RADIUS);

                double center_to_sphere_distance_0 = std::sqrt(std::pow(x_0 - center[0], 2) + std::pow(y_0 - center[1], 2) + std::pow(z_0 - center[2], 2));
                double center_to_sphere_distance_1 = std::sqrt(std::pow(x_1 - center[0], 2) + std::pow(y_1 - center[1], 2) + std::pow(z_1 - center[2], 2));

                if (center_to_sphere_distance_0 < (radius - r_0)) {
                    total_contact_number += 1;
                }

                if (center_to_sphere_distance_1 < (radius - r_1)) {
                    total_contact_number += 1;
                }
            }
        }

        OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), sphere_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

        #pragma omp parallel for reduction(+ : total_particle_number)
        for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){

            for (ElementsArrayType::iterator it = GetElementPartitionBegin(sphere_model_part, k); it != GetElementPartitionEnd(sphere_model_part, k); ++it){

                double r = (it)->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
                double x = (it)->GetGeometry()[0].X();
                double y = (it)->GetGeometry()[0].Y();
                double z = (it)->GetGeometry()[0].Z();

                double center_to_sphere_distance = std::sqrt(std::pow(x - center[0], 2) + std::pow(y - center[1], 2) + std::pow(z - center[2], 2));

                if (center_to_sphere_distance < (radius - r)) {
                    total_particle_number += 1;
                }
            }
        }
        
        double measured_coordination_number = 0.0;

        if (total_particle_number) {
            measured_coordination_number = total_contact_number / total_particle_number;
        }

        return measured_coordination_number;

    }

    double CalculateUnbalancedForceWithinSphere(ModelPart& sphere_model_part, ModelPart& contact_model_part, const double radius, const array_1d<double, 3>& center)
    {
        OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), sphere_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);
        double total_particle_force_modulus_square = 0.0;
        double averaged_total_particle_force_modulus_square = 0.0;
        double particle_number_count = 0;

        #pragma omp parallel for reduction(+ : total_particle_force_modulus_square, particle_number_count)

        for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){

            for (ElementsArrayType::iterator it = GetElementPartitionBegin(sphere_model_part, k); it != GetElementPartitionEnd(sphere_model_part, k); ++it){

                double r = (it)->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
                double x = (it)->GetGeometry()[0].X();
                double y = (it)->GetGeometry()[0].Y();
                double z = (it)->GetGeometry()[0].Z();

                double center_to_sphere_distance = std::sqrt(std::pow(x - center[0], 2) + std::pow(y - center[1], 2) + std::pow(z - center[2], 2));

                if (center_to_sphere_distance < (radius - r)) {
                    const array_1d<double, 3>& total_force = (it)->GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES);
                    double total_force_vector[3] = {total_force[0], total_force[1], total_force[2]};
                    double total_force_vector_modulus = std::sqrt(std::pow(total_force_vector[0], 2) + std::pow(total_force_vector[1], 2) + std::pow(total_force_vector[2], 2));
                    total_particle_force_modulus_square += std::pow(total_force_vector_modulus, 2);
                    particle_number_count += 1;
                }
            }
        }

        if (particle_number_count) {
            averaged_total_particle_force_modulus_square = total_particle_force_modulus_square / particle_number_count;
        }

        OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), contact_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);
        double total_contact_force_modulus_square = 0.0;
        double averaged_contact_force_modulus_square = 0.0;
        double total_contact_number = 0;

        #pragma omp parallel for reduction(+ : total_contact_force_modulus_square, total_contact_number)
        
        for (int k = 0; k < ParallelUtilities::GetNumThreads(); k++){

            for (ElementsArrayType::iterator it = GetElementPartitionBegin(contact_model_part, k); it != GetElementPartitionEnd(contact_model_part, k); ++it){

                double x_0 = (it)->GetGeometry()[0].X();
                double x_1 = (it)->GetGeometry()[1].X();
                double y_0 = (it)->GetGeometry()[0].Y();
                double y_1 = (it)->GetGeometry()[1].Y();
                double z_0 = (it)->GetGeometry()[0].Z();
                double z_1 = (it)->GetGeometry()[1].Z();
                double r_0 = (it)->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
                double r_1 = (it)->GetGeometry()[1].FastGetSolutionStepValue(RADIUS);
                double r   = 0.5 * (r_0 + r_1);

                double center_to_sphere_distance_0 = std::sqrt(std::pow(x_0 - center[0], 2) + std::pow(y_0 - center[1], 2) + std::pow(z_0 - center[2], 2));
                double center_to_sphere_distance_1 = std::sqrt(std::pow(x_1 - center[0], 2) + std::pow(y_1 - center[1], 2) + std::pow(z_1 - center[2], 2));

                if (center_to_sphere_distance_0 < (radius - r) || center_to_sphere_distance_1 < (radius - r)) {
                    const array_1d<double, 3>& contact_force = (it)->GetValue(GLOBAL_CONTACT_FORCE);
                    double contact_force_vector[3] = {contact_force[0], contact_force[1], contact_force[2]};
                    double contact_force_vector_modulus = std::sqrt(std::pow(contact_force_vector[0], 2) + std::pow(contact_force_vector[1], 2) + std::pow(contact_force_vector[2], 2));
                    total_contact_force_modulus_square += std::pow(contact_force_vector_modulus, 2);
                    total_contact_number += 1;
                }
            }
        }

        if (total_contact_number) {
            averaged_contact_force_modulus_square = total_contact_force_modulus_square / total_contact_number;
        }

        double unbalanced_force = 0.0;

        if (averaged_contact_force_modulus_square) {
            unbalanced_force = std::sqrt(averaged_total_particle_force_modulus_square / averaged_contact_force_modulus_square);
        } else {
            unbalanced_force = 0.0;
        }

        return unbalanced_force;
    }

    private:

        std::vector<unsigned int> mElementsPartition;

        ElementsArrayType::iterator GetElementPartitionBegin(ModelPart& r_model_part, unsigned int k)
        {
            ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();
            return (pElements.ptr_begin() + mElementsPartition[k]);
        }
    
        ElementsArrayType::iterator GetElementPartitionEnd(ModelPart& r_model_part, unsigned int k)
        {
            ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();
            return (pElements.ptr_begin() + mElementsPartition[k + 1]);
        }

        ///@name Static Member r_variables
        ///@{


        ///@}
        ///@name Member r_variables
        ///@{


    }; // Class ContactElementGlobalPhysicsCalculator

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // CALCULATE_GLOBAL_PHYSICAL_PROPERTIES_H
