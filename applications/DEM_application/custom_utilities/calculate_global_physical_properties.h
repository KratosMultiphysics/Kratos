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
          OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);
          double added_volume = 0.0;

          #pragma omp parallel for reduction(+ : added_volume)
          for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){

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

      //***************************************************************************************************************
      //***************************************************************************************************************
      // Returns the minimum value of a double variable in the model part.

      double CalculateMaxNodalVariable(ModelPart& r_model_part, const Variable<double>& r_variable)
      {
        ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();

        if (pElements.size() == 0){
            KRATOS_THROW_ERROR(std::invalid_argument, "Cannot compute maximum of the required nodal variable. Empty model part. Could not compute the maximum of the required variable ", r_variable);
          }

        ElementsArrayType::iterator it_begin = pElements.ptr_begin();

        if (!it_begin->GetGeometry()[0].SolutionStepsDataHas(r_variable)){
            KRATOS_THROW_ERROR(std::invalid_argument, "Cannot compute maximum of the required nodal variable. Missing nodal variable ", r_variable);
          }

        Vector max_values;
        double max_val = - std::numeric_limits<double>::max();
        max_values.resize(OpenMPUtils::GetNumThreads());

        for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){
            max_values[k] = max_val;
          }

        OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), pElements.size(), mElementsPartition);

        unsigned int elem_counter;

        #pragma omp parallel for private(elem_counter)
        for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){
            elem_counter = mElementsPartition[k];

            for (ElementsArrayType::iterator it = GetElementPartitionBegin(r_model_part, k); it != GetElementPartitionEnd(r_model_part, k); ++it){
                max_values[k] = std::max(max_values[k], (it)->GetGeometry()[0].FastGetSolutionStepValue(r_variable));
                elem_counter++;

              }

          }

        // getting the maximum between threads:

        for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){
            max_val = std::max(max_val, max_values[k]);
          }

        return max_val;
      }

      //***************************************************************************************************************
      //***************************************************************************************************************

      double CalculateMinNodalVariable(ModelPart& r_model_part, const Variable<double>& r_variable)
      {
        ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();

        if (pElements.size() == 0){
            KRATOS_THROW_ERROR(std::invalid_argument, "Cannot compute minimum of the required nodal variable. Empty model part. Could not compute the maximum of the required variable ", r_variable);
          }

        ElementsArrayType::iterator it_begin = pElements.ptr_begin();

        if (!it_begin->GetGeometry()[0].SolutionStepsDataHas(r_variable)){
            KRATOS_THROW_ERROR(std::invalid_argument, "Cannot compute minimum of the required nodal variable. Missing variable ", r_variable);
          }

        Vector min_values;
        double min_val = std::numeric_limits<double>::max();
        min_values.resize(OpenMPUtils::GetNumThreads());

        for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){
            min_values[k] = min_val;
          }

        OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), pElements.size(), mElementsPartition);

        unsigned int elem_counter;

        #pragma omp parallel for private(elem_counter)
        for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){
            elem_counter = mElementsPartition[k];

            for (ElementsArrayType::iterator it = GetElementPartitionBegin(r_model_part, k); it != GetElementPartitionEnd(r_model_part, k); ++it){
                min_values[k] = std::min(min_values[k], (it)->GetGeometry()[0].FastGetSolutionStepValue(r_variable));
                elem_counter++;

              }

          }

        // getting the minimum between threads:

        for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){
            min_val = std::min(min_val, min_values[k]);
          }

        return min_val;
      }

      //***************************************************************************************************************
      //***************************************************************************************************************

      double CalculateD50(ModelPart& r_model_part)
      {
          const unsigned int size = r_model_part.GetCommunicator().LocalMesh().Elements().size();
          OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), size, mElementsPartition);

          Vector radii;
          radii.resize(size);
          unsigned int particle_counter = 0;

          #pragma omp parallel for private(particle_counter)
          for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){
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
          OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(),r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

          double added_mass = 0.0;

          #pragma omp parallel for reduction(+ : added_mass)
          for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){

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
          OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

          const double total_mass_inv         = 1 / CalculateTotalMass(r_model_part);
          double cm_x = 0.0;
          double cm_y = 0.0;
          double cm_z = 0.0;

          #pragma omp parallel for reduction(+ : cm_x, cm_y, cm_z)
          for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){

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
          OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

          double kinematic_energy = 0.0;

          #pragma omp parallel for reduction(+ : kinematic_energy)
          for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){

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
          OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

          double rotational_kinematic_energy = 0.0;

          #pragma omp parallel for reduction(+ : rotational_kinematic_energy)
          for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){

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
          OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

          double elastic_energy = 0.0;

          #pragma omp parallel for reduction(+ : elastic_energy)
          for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){

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
          OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

          double frictional_energy = 0.0;

          #pragma omp parallel for reduction(+ : frictional_energy)
          for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){

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
          OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

          double viscodamping_energy = 0.0;

          #pragma omp parallel for reduction(+ : viscodamping_energy)
          for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){

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

      //***************************************************************************************************************
      //***************************************************************************************************************

      array_1d<double, 3> CalculateTotalMomentum(ModelPart& r_model_part)
      {
          OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

          double m_x = 0.0;
          double m_y = 0.0;
          double m_z = 0.0;

          #pragma omp parallel for reduction(+ : m_x, m_y, m_z)
          for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){

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
          OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);

          const array_1d<double, 3> center_of_mass   = CalculateCenterOfMass(r_model_part);
          double am_x = 0.0;
          double am_y = 0.0;
          double am_z = 0.0;

          #pragma omp parallel for reduction(+ : am_x, am_y, am_z)
          for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){

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
            OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(),r_model_part.GetCommunicator().LocalMesh().Elements().size(), mElementsPartition);
            double sum_of_contact_forces_x = 0.0;
            double sum_of_contact_forces_y = 0.0;
            double sum_of_contact_forces_z = 0.0;

            #pragma omp parallel for reduction(+ : sum_of_contact_forces_x, sum_of_contact_forces_y, sum_of_contact_forces_z)
            for (int k = 0; k < OpenMPUtils::GetNumThreads(); ++k){

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
      
      void SetGravity(ModelPart& r_model_part) {

        static bool called_for_the_first_time = true;
        static double current_reference_time = r_model_part.GetProcessInfo()[TIME];
        std::string gravities_filename = "TimeAngle.csv";

        static std::ifstream input_file(gravities_filename);
        std::string line;
        static size_t gravity_number = 0;
        double elapsed_time = r_model_part.GetProcessInfo()[TIME] - current_reference_time;
        // Minimum time span to wait between gravity shifts to ensure that the powder realocates
        const double delta_security_time = 1.0;
        // We choose the maximum velocity admissible in order to change the gravity vector
        const double maximum_squared_velocity_module = 0.05 * 0.05; // squares will be compared
        static bool is_time_to_change_gravity_value = true;
        double node_i_squared_velocity_module = 0.0;
        // There are 113 different gravities in TimeAngle.csv
        const size_t total_number_of_gravities = 113;
        
        if (!called_for_the_first_time) is_time_to_change_gravity_value = true;

        ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();

        ElementsArrayType::iterator it_begin = pElements.ptr_begin();
        ElementsArrayType::iterator it_end = pElements.ptr_end();

        for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it) {
                   
            array_1d<double, 3> velocity = (it)->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);

            node_i_squared_velocity_module = velocity[0] * velocity[0] + velocity[1] * velocity[1] + velocity[2] * velocity[2];

            if (node_i_squared_velocity_module > maximum_squared_velocity_module) {

                is_time_to_change_gravity_value = false;
                break;
            }                
        }

        if ((called_for_the_first_time) or (is_time_to_change_gravity_value and (elapsed_time > delta_security_time))) {

            if (gravity_number > total_number_of_gravities) {
            
                std::cout << "\n\nNo more gravity vectors left... Simulation should be stopped soon...\n\n";
                is_time_to_change_gravity_value = false;
                current_reference_time = r_model_part.GetProcessInfo()[TIME];
                return;
            }

            std::getline(input_file, line);
            std::istringstream iss{line};
            std::vector<std::string> tokens;
            std::string token;

            while (std::getline(iss, token, ',')) tokens.push_back(token);
            
            array_1d<double, 3>& gravity = r_model_part.GetProcessInfo()[GRAVITY];

            const double gravity_x = std::stof(tokens[3]);
            const double gravity_y = std::stof(tokens[4]);
            const double gravity_z = std::stof(tokens[5]);

            gravity[0] = 9.81 * gravity_x;
            gravity[1] = 9.81 * gravity_y;
            gravity[2] = 9.81 * gravity_z;

            is_time_to_change_gravity_value = called_for_the_first_time = false;
            current_reference_time = r_model_part.GetProcessInfo()[TIME];
            
            ++gravity_number;
            
            std::cout << "\n**********CHANGING GRAVITY VECTOR...\nGravity number " << gravity_number << "\n";
            KRATOS_WATCH(gravity[0])
            KRATOS_WATCH(gravity[1])
            KRATOS_WATCH(gravity[2])
            std::cout << "\n";
        }
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

        DenseVector<unsigned int>& GetElementPartition()
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
        DenseVector<unsigned int> mElementsPartition;

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

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

} // namespace Kratos.

#endif // CALCULATE_GLOBAL_PHYSICAL_PROPERTIES_H
