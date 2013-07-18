#ifndef CALCULATE_GLOBAL_PHYSICAL_PROPERTIES_H
#define CALCULATE_GLOBAL_PHYSICAL_PROPERTIES_H

// /* External includes */

// System includes

// Project includes
#include "utilities/timer.h"
#include "custom_utilities/create_and_destroy.h"

#include "custom_elements/spheric_swimming_particle.h"
#include "custom_elements/Particle_Contact_Element.h"

#include "includes/variables.h"
#include "DEM_application.h"

/* System includes */
#include <limits>
#include <iostream>
#include <iomanip>
#include <iostream>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{

class SphericElementGlobalPhysicsCalculator
    {
     public:

     typedef ModelPart::ElementsContainerType                          ElementsArrayType;
     typedef ModelPart::NodesContainerType::ContainerType              NodesContainerType;


     KRATOS_CLASS_POINTER_DEFINITION(SphericElementGlobalPhysicsCalculator);

      /// Default constructor.

      SphericElementGlobalPhysicsCalculator(ModelPart& rModelPart)
      {
          mInitialCenterOfMassAndMass = CalculateCenterOfMass(rModelPart);
          mInitialMass                = CalculateTotalMass(rModelPart);
      }

      /// Destructor.

      virtual ~SphericElementGlobalPhysicsCalculator(){}

      double CalculateTotalVolume(ModelPart& rModelPart)
      {
          ElementsArrayType& pElements        = rModelPart.GetCommunicator().LocalMesh().Elements();

          OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), pElements.size(), this->GetElementPartition());
          double added_volume = 0.0;

          #pragma omp parallel for
          for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){
              ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
              ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

              for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                  const double& particle_radius = (it)->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
                  double particle_volume = 4/3 * pi * particle_radius * particle_radius * particle_radius;
                  added_volume += particle_volume;
              }

          }

        return added_volume;
      }

      double CalculateTotalMass(ModelPart& rModelPart)
      {
          ElementsArrayType& pElements        = rModelPart.GetCommunicator().LocalMesh().Elements();

          OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), pElements.size(), this->GetElementPartition());
          double added_mass = 0.0;

          #pragma omp parallel for
          for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){
              ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
              ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

              for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                  double particle_mass = (it)->GetGeometry()(0)->FastGetSolutionStepValue(NODAL_MASS);
                  added_mass += particle_mass;
              }

          }

        return added_mass;
      }

      array_1d<double, 3> CalculateCenterOfMass(ModelPart& rModelPart)
      {
          ElementsArrayType& pElements        = rModelPart.GetCommunicator().LocalMesh().Elements();

          OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), pElements.size(), this->GetElementPartition());
          array_1d<double, 3> center_of_mass;
          double total_mass = CalculateTotalMass(rModelPart);
          center_of_mass[0] = 0.0;
          center_of_mass[1] = 0.0;
          center_of_mass[2] = 0.0;

          #pragma omp parallel for
          for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){
              ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
              ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

              for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                  double particle_mass = (it)->GetGeometry()(0)->FastGetSolutionStepValue(NODAL_MASS);
                  center_of_mass += particle_mass * (it)->GetGeometry()(0)->Coordinates();
              }

          }

          center_of_mass /= total_mass;

          return center_of_mass;
      }

      double CalculateGravitationalPotentialEnergy(ModelPart& rModelPart, const array_1d<double, 3> reference_point)
      {

          ProcessInfo& rCurrentProcessInfo                = rModelPart.GetProcessInfo();
          const array_1d<double, 3>& gravity              = rCurrentProcessInfo[GRAVITY];
          array_1d<double, 3> center_of_mass              = CalculateCenterOfMass(rModelPart);
          double total_mass                               = CalculateTotalMass(rModelPart);
          array_1d<double, 3> center_of_mass_to_reference = reference_point - center_of_mass;

          double potential_energy = total_mass * (center_of_mass_to_reference[0] * gravity[0] + center_of_mass_to_reference[1] * gravity[1] + center_of_mass_to_reference[2] * gravity[2]);

          return potential_energy;
      }

      double CalculateKineticEnergy(ModelPart& rModelPart)
      {
          ElementsArrayType& pElements        = rModelPart.GetCommunicator().LocalMesh().Elements();
          ProcessInfo& rCurrentProcessInfo    = rModelPart.GetProcessInfo();
          OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), pElements.size(), this->GetElementPartition());
          double kinetic_energy   = 0.0;

          #pragma omp parallel for
          for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){
              ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
              ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

              for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                  double particle_kinetic_energy            = 0.0;

                  (it)->Calculate(KINETIC_ENERGY, particle_kinetic_energy, rCurrentProcessInfo);

                  kinetic_energy   += particle_kinetic_energy;

              }

          }

          return kinetic_energy;
      }


      double CalculateElasticEnergy(ModelPart& rModelPart)
      {
          ElementsArrayType& pElements        = rModelPart.GetCommunicator().LocalMesh().Elements();
          ProcessInfo& rCurrentProcessInfo    = rModelPart.GetProcessInfo();
          OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), pElements.size(), this->GetElementPartition());
          double potential_energy = 0.0;

          #pragma omp parallel for
          for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){
              ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
              ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

              for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                  double particle_contacts_potential_energy = 0.0;

                  (it)->Calculate(ELASTIC_ENERGY_OF_CONTACTS, particle_contacts_potential_energy, rCurrentProcessInfo);

                  potential_energy += 0.5 * particle_contacts_potential_energy;

              }

          }

          return potential_energy;
      }

      array_1d<double, 3> CalculateTotalMomentum(ModelPart& rModelPart)
      {
          ElementsArrayType& pElements        = rModelPart.GetCommunicator().LocalMesh().Elements();
          ProcessInfo& rCurrentProcessInfo    = rModelPart.GetProcessInfo();
          OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), pElements.size(), this->GetElementPartition());
          array_1d<double, 3> momentum;
          momentum[0] = 0.0;
          momentum[1] = 0.0;
          momentum[2] = 0.0;

          #pragma omp parallel for
          for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){
              ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
              ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

              for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                  array_1d<double, 3> particle_momentum;
                  particle_momentum[0] = 0.0;
                  particle_momentum[1] = 0.0;
                  particle_momentum[2] = 0.0;

                  (it)->Calculate(MOMENTUM, particle_momentum, rCurrentProcessInfo);

                  momentum += particle_momentum;

              }

          }

          return momentum;
      }

      array_1d<double, 3> CalulateTotalAngularMomentum(ModelPart& rModelPart)
      {
          ElementsArrayType& pElements        = rModelPart.GetCommunicator().LocalMesh().Elements();
          ProcessInfo& rCurrentProcessInfo    = rModelPart.GetProcessInfo();
          OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), pElements.size(), this->GetElementPartition());
          array_1d<double, 3> center_of_mass;
          array_1d<double, 3> angular_momentum;

          center_of_mass = CalculateCenterOfMass(rModelPart);
          angular_momentum[0] = 0.0;
          angular_momentum[1] = 0.0;
          angular_momentum[2] = 0.0;

          #pragma omp parallel for
          for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++){
              ElementsArrayType::iterator it_begin = pElements.ptr_begin() + this->GetElementPartition()[k];
              ElementsArrayType::iterator it_end   = pElements.ptr_begin() + this->GetElementPartition()[k + 1];

              for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it){
                  array_1d<double, 3> particle_momentum;
                  array_1d<double, 3> particle_local_angular_momentum;
                  array_1d<double, 3> particle_coordinates       = (it)->GetGeometry()(0)->Coordinates();
                  array_1d<double, 3> center_of_mass_to_particle = particle_coordinates - center_of_mass;

                  (it)->Calculate(MOMENTUM, particle_momentum, rCurrentProcessInfo);
                  (it)->Calculate(ANGULAR_MOMENTUM, particle_local_angular_momentum, rCurrentProcessInfo);

                  angular_momentum += particle_local_angular_momentum + Kratos::MathUtils<double>::CrossProduct(center_of_mass_to_particle, particle_momentum);
              }

          }

          return angular_momentum;
      }


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
        vector<unsigned int>&    GetElementPartition(){return (mElementPartition);};
        ///@}

    protected:
        ///@name Protected static Member rVariables
        ///@{


        ///@}
        ///@name Protected member rVariables
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
        vector<unsigned int>                mElementPartition;

        ///@}
        ///@name Protected Inquiry
        ///@{


        ///@}
        ///@name Protected LifeCycle
        ///@{


        ///@}

    private:


        ///@name Static Member rVariables
        ///@{


        ///@}
        ///@name Member rVariables
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




/// output stream function
// 	template<std::size_t TDim>
// 	inline std::ostream& operator << (std::ostream& rOStream)
// 	{
// 		rThis.PrintInfo(rOStream);
// 		rOStream << std::endl;
// 		rThis.PrintData(rOStream);
//
// 		return rOStream;
// 	}
///@}


} // namespace Kratos.

#endif // CALCULATE_GLOBAL_PHYSICAL_PROPERTIES_H
