//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:           MSantasusana $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:               March 2015 $
//   Revision:            $Revision:                  0.0 $
//
//


#if !defined(KRATOS_ENERGY_UTILITIES_H_INCLUDED)
#define  KRATOS_ENERGY_UTILITIES_H_INCLUDED


// System includes
//#include <cmath>
//#include <set>

// #ifdef _OPENMP
// #include <omp.h>
// #endif

// External includes
//#include "boost/smart_ptr.hpp"
#include "utilities/timer.h"

// Project includes
//#include "includes/define.h"
#include "includes/variables.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "utilities/math_utils.h"
#include "solid_mechanics_application_variables.h"

namespace Kratos
{

  /// Short class definition.

  /** Computes the energy
   */

  class EnergyUtilities
  {
  public:

    ///@name Type Definitions
    ///@{

    typedef ModelPart::ElementsContainerType                  ElementsContainerType;
    typedef ModelPart::MeshType::GeometryType                          GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    EnergyUtilities(){ mEchoLevel = 0;  mParallel = true; };

    EnergyUtilities(bool Parallel){ mEchoLevel = 0;  mParallel = Parallel; };


    /// Destructor.
    virtual ~EnergyUtilities(){};


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    //**************************************************************************
    //**************************************************************************


    double GetTotalKinematicEnergy(ModelPart& rModelPart)
    {
        KRATOS_TRY

        double KinematicEnergy = 0.0;

        array_1d<double, 3> vel;
        double mass = 0.0;
        double vel_arg = 0.0;

        ModelPart::NodesContainerType& pNodes = rModelPart.Nodes();

        #ifdef _OPENMP
                int number_of_threads = omp_get_max_threads();
        #else
                int number_of_threads = 1;
        #endif

        if( mParallel == false ){ number_of_threads = 1; }

        vector<unsigned int> node_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

        std::vector<double> KinematicEnergyPartition(number_of_threads);
        for(int i=0; i<number_of_threads; i++){
            KinematicEnergyPartition[i] = 0.0;
        }

        #pragma omp parallel private (KinematicEnergy, vel, mass, vel_arg)
        {
            int k = OpenMPUtils::ThisThread();
            ModelPart::NodesContainerType::iterator NodeBegin = pNodes.begin() + node_partition[k];
            ModelPart::NodesContainerType::iterator NodeEnd = pNodes.begin() + node_partition[k + 1];


            for (ModelPart::NodesContainerType::const_iterator in = NodeBegin; in != NodeEnd; in++){

                mass    = in->FastGetSolutionStepValue(NODAL_MASS);
                vel     = in->FastGetSolutionStepValue(VELOCITY);
                vel_arg = MathUtils<double>::Norm3(vel);
                KinematicEnergyPartition[k] += 0.5*vel_arg*vel_arg*mass;
            }
        }

        for(int i=0; i<number_of_threads; i++){
            KinematicEnergy   += KinematicEnergyPartition[i];

        }

        return KinematicEnergy;


        KRATOS_CATCH( "" )
    }



    double GetGravitationalEnergy(ModelPart& rModelPart)
    {
        KRATOS_TRY

        double PotentialEnergy = 0.0;
        double gh = 0.0;


        array_1d<double, 3> coord;
        double mass = 0.0;
        Vector VolumeAcceleration(3);
	noalias(VolumeAcceleration) = ZeroVector(3);

        ModelPart::NodesContainerType& pNodes = rModelPart.Nodes();

        #ifdef _OPENMP
                int number_of_threads = omp_get_max_threads();
        #else
                int number_of_threads = 1;
        #endif

        if( mParallel == false ){ number_of_threads = 1; }

        vector<unsigned int> node_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

        std::vector<double> PotentialEnergyPartition(number_of_threads);
        for(int i=0; i<number_of_threads; i++){
            PotentialEnergyPartition[i] = 0.0;
        }

        #pragma omp parallel private (gh, VolumeAcceleration, mass, coord)
        {
            int k = OpenMPUtils::ThisThread();
            ModelPart::NodesContainerType::iterator NodeBegin = pNodes.begin() + node_partition[k];
            ModelPart::NodesContainerType::iterator NodeEnd = pNodes.begin() + node_partition[k + 1];


            for (ModelPart::NodesContainerType::const_iterator in = NodeBegin; in != NodeEnd; in++)
            {

                mass = in->FastGetSolutionStepValue(NODAL_MASS);
                coord = in->Coordinates();

                gh = 0.0;
                VolumeAcceleration = in->FastGetSolutionStepValue(VOLUME_ACCELERATION);  //if( in->SolutionStepsDataHas(VOLUME_ACCELERATION) ){ //MSI:: should be checked once at the beginning only

                for (unsigned int i=0;i<3;i++){
                    gh += coord(i)*-VolumeAcceleration(i); //negative makes gravity introduce positive potential energy
                }

                PotentialEnergyPartition[k] += mass*gh;

            }
        }

        for(int i=0; i<number_of_threads; i++){
            PotentialEnergy += PotentialEnergyPartition[i];

        }

        return PotentialEnergy;


        KRATOS_CATCH( "" )
    }


    //Total Energy on initial configuration
    double GetTotalStrainEnergy(ModelPart& rModelPart)
    {

        KRATOS_TRY

        double StrainEnergy = 0.0;
        ProcessInfo& rCurrentProcessInfo                  = rModelPart.GetProcessInfo();

        ModelPart::ElementsContainerType& pElements = rModelPart.Elements();

        #ifdef _OPENMP
                int number_of_threads = omp_get_max_threads();
        #else
                int number_of_threads = 1;
        #endif

        if( mParallel == false ){ number_of_threads = 1; }

        vector<unsigned int> element_partition;
        OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);
        std::vector<double> StrainEnergyPartition(number_of_threads);
        for(int i=0; i<number_of_threads; i++){
            StrainEnergyPartition[i] = 0.0;
        }

        #pragma omp parallel
        {
            int k = OpenMPUtils::ThisThread();
            ModelPart::ElementsContainerType::iterator ElemBegin = pElements.begin() + element_partition[k];
            ModelPart::ElementsContainerType::iterator ElemEnd = pElements.begin() + element_partition[k + 1];

            for(ModelPart::ElementsContainerType::const_iterator ie = ElemBegin; ie != ElemEnd; ie++)
            {
                const GeometryType::IntegrationMethod IntegrationMethod = ie->GetGeometry().GetDefaultIntegrationMethod();
                const GeometryType::IntegrationPointsArrayType& integration_points = ie->GetGeometry().IntegrationPoints( IntegrationMethod );

                std::vector<double> rOutput;
                ie->CalculateOnIntegrationPoints(STRAIN_ENERGY,rOutput,rCurrentProcessInfo );

                for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
                {

                    StrainEnergyPartition[k] += rOutput[PointNumber];

                }

            }

        }

        for(int i=0; i<number_of_threads; i++){
            StrainEnergy += StrainEnergyPartition[i];

        }

        if(StrainEnergy==0){
                std::cout<<"   [ NO_STRAIN_ENERGY] "<<std::endl;
        }

        return StrainEnergy;


        KRATOS_CATCH( "" )
    }



    double GetExternallyAppliedEnergy(ModelPart& rModelPart)
    {
      KRATOS_TRY

     //work done by conditions...

     return 0.0;
      KRATOS_CATCH( "" )
    }

    void CalculateNodalMass(ModelPart::NodesContainerType& pNodes, ModelPart::ElementsContainerType& pElements, ProcessInfo& rCurrentProcessInfo)
    {

      KRATOS_TRY
      #ifdef _OPENMP
              int number_of_threads = omp_get_max_threads();
      #else
              int number_of_threads = 1;
      #endif

      vector<unsigned int> node_partition;
      OpenMPUtils::CreatePartition(number_of_threads, pNodes.size(), node_partition);

      vector<unsigned int> element_partition;
      OpenMPUtils::CreatePartition(number_of_threads, pElements.size(), element_partition);


      #pragma omp parallel
      {

        #pragma omp for

          for(int k=0; k<number_of_threads; k++)
          {
              ModelPart::NodesContainerType::iterator i_begin=pNodes.ptr_begin()+node_partition[k];
              ModelPart::NodesContainerType::iterator i_end=pNodes.ptr_begin()+node_partition[k+1];

              for(ModelPart::NodeIterator i=i_begin; i!= i_end; ++i)
              {
                  double& nodal_mass    =  i->FastGetSolutionStepValue(NODAL_MASS);
                  nodal_mass = 0.0;
              }
          }

      }

      //Calculate and assemble Mass Matrix on nodes

      unsigned int index = 0;

      #pragma omp parallel
      {
          int k = OpenMPUtils::ThisThread();
          ModelPart::ElementsContainerType::iterator ElemBegin = pElements.begin() + element_partition[k];
          ModelPart::ElementsContainerType::iterator ElemEnd = pElements.begin() + element_partition[k + 1];

          for (ModelPart::ElementsContainerType::iterator itElem = ElemBegin; itElem != ElemEnd; itElem++)  //MSI: To be parallelized
          {
              Matrix MassMatrix;

              Element::GeometryType& geometry = itElem->GetGeometry();

              (itElem)->CalculateMassMatrix(MassMatrix, rCurrentProcessInfo); //already lumped.

              const unsigned int dimension   = geometry.WorkingSpaceDimension();

              index = 0;
              for (unsigned int i = 0; i <geometry.size(); i++)
              {
                  index = i*dimension;

                  double& mass = geometry(i)->FastGetSolutionStepValue(NODAL_MASS);

                  geometry(i)->SetLock();
                  mass += MassMatrix(index,index);
                  geometry(i)->UnSetLock();
              }
          }
      }

      KRATOS_CATCH( "" )

    } //CalculateNodalMass



    //**************************************************************************
    //**************************************************************************

    /**
     * level of echo for the error calculation
     */
    virtual void SetEchoLevel(int Level)
    {
      mEchoLevel = Level;
    }

    int GetEchoLevel()
    {
      return mEchoLevel;
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{


    ///@}


  protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}
  private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    int mEchoLevel;

    bool mParallel;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    //************************************************************************************
    //************************************************************************************


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Unaccessible methods
    ///@{

    ///@}


  }; // Class EnergyUtilities

  ///@}
  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}

} // namespace Kratos.

#endif // KRATOS_ENERGY_UTILITIES_H_INCLUDED defined


