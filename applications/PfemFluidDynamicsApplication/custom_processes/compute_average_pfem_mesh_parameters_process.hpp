//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:                   AFranci $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:                October 2019 $
//   Revision:            $Revision:                     0.0 $
//
//

#if !defined(KRATOS_COMPUTE_AVERAGE_PFEM_MESH_PARAMETERS_PROCESS_H_INCLUDED)
#define KRATOS_COMPUTE_AVERAGE_PFEM_MESH_PARAMETERS_PROCESS_H_INCLUDED

// External includes

// System includes

// Project includes
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"

#include "includes/model_part.h"
#include "custom_utilities/mesh_error_calculation_utilities.hpp"
#include "custom_utilities/mesher_utilities.hpp"

#include "pfem_fluid_dynamics_application_variables.h"
#include "custom_processes/mesher_process.hpp"

/// VARIABLES used:
// Data:
// Flags:    (checked)
//           (set)
//           (modified)
//           (reset)
//(set):=(set in this process)

namespace Kratos
{

  ///@name Kratos Classes
  ///@{

  /// Refine Mesh Elements Process 2D and 3D
  /** The process labels the nodes to be refined (TO_REFINE)
      if the ThresholdVariable  is larger than a ReferenceThreshold
  */

  class ComputeAveragePfemMeshParametersProcess
      : public MesherProcess
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(ComputeAveragePfemMeshParametersProcess);

    typedef ModelPart::NodeType NodeType;
    typedef ModelPart::ConditionType ConditionType;
    typedef ModelPart::PropertiesType PropertiesType;
    typedef ConditionType::GeometryType GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ComputeAveragePfemMeshParametersProcess(ModelPart &rModelPart,
                                            MesherUtilities::MeshingParameters &rRemeshingParameters,
                                            int EchoLevel)
        : mrModelPart(rModelPart),
          mrRemesh(rRemeshingParameters)
    {
      KRATOS_INFO("ComputeAveragePfemMeshParametersProcess") << " activated " << std::endl;

      mEchoLevel = EchoLevel;
    }

    /// Destructor.
    virtual ~ComputeAveragePfemMeshParametersProcess() {}

    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
      Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /// Execute method is used to execute the Process algorithms.
    void Execute() override
    {
      KRATOS_TRY

      if (mEchoLevel > 1)
        std::cout << "  COMPUTE AVERAGE PFEM MESH PARAMETERS PROCESS ]; " << std::endl;

      const unsigned int dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
      const unsigned int numberOfRefiningBoxes = mrRemesh.UseRefiningBox.size();
      Vector inBoxesNodes(numberOfRefiningBoxes);
      noalias(inBoxesNodes) = ZeroVector(numberOfRefiningBoxes);
      Vector inBoxesMeanNodalSize(numberOfRefiningBoxes);
      noalias(inBoxesMeanNodalSize) = ZeroVector(numberOfRefiningBoxes);
      double preliminaryOutOfBoxesFluidNodes = 0;
      double preliminaryOutOfBoxesMeanNodalSize = 0;
      bool homogeneousMesh = true;
      for (ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin(); i_node != mrModelPart.NodesEnd(); i_node++)
      {
        if (i_node->Is(FLUID))
        {
          if (numberOfRefiningBoxes == 0 || (numberOfRefiningBoxes == 1 && mrRemesh.UseRefiningBox[0] == false))
          {
            preliminaryOutOfBoxesFluidNodes += 1.0;
            preliminaryOutOfBoxesMeanNodalSize += i_node->FastGetSolutionStepValue(NODAL_H);
          }
          else
          {
            unsigned outOfRefiningBoxes = true;
            for (unsigned int index = 0; index < numberOfRefiningBoxes; index++)
            {
              array_1d<double, 3> RefiningBoxMinimumPoint = mrRemesh.RefiningBoxMinimumPoint[index];
              array_1d<double, 3> RefiningBoxMaximumPoint = mrRemesh.RefiningBoxMaximumPoint[index];
              if (mrRemesh.UseRefiningBox[index] == true)
              {
                homogeneousMesh = false;
                if (dimension == 2)
                {
                  if (i_node->X() > RefiningBoxMinimumPoint[0] && i_node->Y() > RefiningBoxMinimumPoint[1] &&
                      i_node->X() < RefiningBoxMaximumPoint[0] && i_node->Y() < RefiningBoxMaximumPoint[1])
                  {
                    outOfRefiningBoxes = false;
                    inBoxesNodes[index] += 1.0;
                    inBoxesMeanNodalSize[index] += i_node->FastGetSolutionStepValue(NODAL_H); // this is a preliminary evaluation of the local mesh size
                  }
                }
                else if (dimension == 3)
                {
                  if (i_node->X() > RefiningBoxMinimumPoint[0] && i_node->Y() > RefiningBoxMinimumPoint[1] && i_node->Z() > RefiningBoxMinimumPoint[2] &&
                      i_node->X() < RefiningBoxMaximumPoint[0] && i_node->Y() < RefiningBoxMaximumPoint[1] && i_node->Z() < RefiningBoxMaximumPoint[2])
                  {
                    outOfRefiningBoxes = false;
                    inBoxesNodes[index] += 1.0;
                    inBoxesMeanNodalSize[index] += i_node->FastGetSolutionStepValue(NODAL_H); // this is a preliminary evaluation of the local mesh size
                  }
                }
              }
            }
            // CONSIDER ONLY THE NODES OUT FROM THE REFINEMENT AREAS
            if (outOfRefiningBoxes == true)
            {
              preliminaryOutOfBoxesFluidNodes += 1.0;
              preliminaryOutOfBoxesMeanNodalSize += i_node->FastGetSolutionStepValue(NODAL_H); // this is a preliminary evaluation of the local mesh size
            }
          }
        }
      }
      preliminaryOutOfBoxesMeanNodalSize *= 1.0 / preliminaryOutOfBoxesFluidNodes;

      mrRemesh.Refine->CriticalRadius = preliminaryOutOfBoxesMeanNodalSize;
      mrRemesh.Refine->InitialRadius = preliminaryOutOfBoxesMeanNodalSize;

      if (homogeneousMesh == false)
      {
        double outOfBoxesFluidNodes = 0;
        double outOfBoxesMeanNodalSize = 0;
        for (ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin(); i_node != mrModelPart.NodesEnd(); i_node++)
        {
          if (i_node->Is(FLUID))
          {
            unsigned outOfRefiningBoxes = true;
            for (unsigned int index = 0; index < numberOfRefiningBoxes; index++)
            {
              const double transitionDistanceInInputMesh = mrRemesh.RefiningBoxElementsInTransitionZone[index] * mrRemesh.Refine->CriticalRadius;
              array_1d<double, 3> RefiningBoxMinimumPoint = mrRemesh.RefiningBoxMinimumPoint[index];
              array_1d<double, 3> RefiningBoxMaximumPoint = mrRemesh.RefiningBoxMaximumPoint[index];
              if (mrRemesh.UseRefiningBox[index] == true)
              {
                if (dimension == 2)
                {
                  if (i_node->X() > (RefiningBoxMinimumPoint[0] - transitionDistanceInInputMesh) && i_node->Y() > (RefiningBoxMinimumPoint[1] - transitionDistanceInInputMesh) &&
                      i_node->X() < (RefiningBoxMaximumPoint[0] + transitionDistanceInInputMesh) && i_node->Y() < (RefiningBoxMaximumPoint[1] + transitionDistanceInInputMesh))
                  {
                    outOfRefiningBoxes = false;
                  }
                }
                else if (dimension == 3)
                {
                  if (i_node->X() > (RefiningBoxMinimumPoint[0] - transitionDistanceInInputMesh) && i_node->Y() > (RefiningBoxMinimumPoint[1] - transitionDistanceInInputMesh) && i_node->Z() > (RefiningBoxMinimumPoint[2] - transitionDistanceInInputMesh) &&
                      i_node->X() < (RefiningBoxMaximumPoint[0] + transitionDistanceInInputMesh) && i_node->Y() < (RefiningBoxMaximumPoint[1] + transitionDistanceInInputMesh) && i_node->Z() < (RefiningBoxMaximumPoint[2] + transitionDistanceInInputMesh))
                  {
                    outOfRefiningBoxes = false;
                  }
                }
              }
            }
            // CONSIDER ONLY THE NODES OUT FROM THE REFINEMENT AREAS
            if (outOfRefiningBoxes == true)
            {
              outOfBoxesFluidNodes += 1.0;
              outOfBoxesMeanNodalSize += i_node->FastGetSolutionStepValue(NODAL_H);
            }
          }
        }
        if (outOfBoxesFluidNodes == 0)
        {
          mrRemesh.Refine->CriticalRadius = preliminaryOutOfBoxesMeanNodalSize;
          mrRemesh.Refine->InitialRadius = preliminaryOutOfBoxesMeanNodalSize;
          std::cout << "the coarse zone is too thin, I'll take the preliminary mesh size estimation: " << preliminaryOutOfBoxesMeanNodalSize << std::endl;
        }
        else
        {
          outOfBoxesMeanNodalSize *= 1.0 / outOfBoxesFluidNodes;
          mrRemesh.Refine->CriticalRadius = outOfBoxesMeanNodalSize;
          mrRemesh.Refine->InitialRadius = outOfBoxesMeanNodalSize;
        }

        std::cout << "Mesh size outside the refining boxes is: " << outOfBoxesMeanNodalSize << std::endl;
        for (unsigned int index = 0; index < numberOfRefiningBoxes; index++)
        {
          double localMeshSize = inBoxesMeanNodalSize[index] / inBoxesNodes[index];
          std::cout << "Mesh size inside refining box n." << index << " is: " << localMeshSize << std::endl;

          mrRemesh.SetRefiningBoxMeshSize(index, localMeshSize);
          const double tolerance = mrRemesh.RefiningBoxMeshSize[index] * 0.01;
          const double differenceOfSize = outOfBoxesMeanNodalSize - mrRemesh.RefiningBoxMeshSize[index];

          mrRemesh.RefiningBoxMinimumPoint[index][0] += -tolerance; // the finest nodes at the frontier should not be erased
          mrRemesh.RefiningBoxMinimumPoint[index][1] += -tolerance;
          mrRemesh.RefiningBoxMinimumPoint[index][2] += -tolerance;

          mrRemesh.RefiningBoxMaximumPoint[index][0] += tolerance;
          mrRemesh.RefiningBoxMaximumPoint[index][1] += tolerance;
          mrRemesh.RefiningBoxMaximumPoint[index][2] += tolerance;

          const double transitionDistance = mrRemesh.RefiningBoxElementsInTransitionZone[index] * std::abs(differenceOfSize);

          mrRemesh.RefiningBoxShiftedMinimumPoint[index][0] = mrRemesh.RefiningBoxMinimumPoint[index][0] - transitionDistance;
          mrRemesh.RefiningBoxShiftedMinimumPoint[index][1] = mrRemesh.RefiningBoxMinimumPoint[index][1] - transitionDistance;
          mrRemesh.RefiningBoxShiftedMinimumPoint[index][2] = mrRemesh.RefiningBoxMinimumPoint[index][2] - transitionDistance;

          mrRemesh.RefiningBoxShiftedMaximumPoint[index][0] = mrRemesh.RefiningBoxMaximumPoint[index][0] + transitionDistance;
          mrRemesh.RefiningBoxShiftedMaximumPoint[index][1] = mrRemesh.RefiningBoxMaximumPoint[index][1] + transitionDistance;
          mrRemesh.RefiningBoxShiftedMaximumPoint[index][2] = mrRemesh.RefiningBoxMaximumPoint[index][2] + transitionDistance;
        }
      }

      KRATOS_CATCH(" ")
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

    /// Turn back information as a string.
    std::string Info() const override
    {
      return "ComputeAveragePfemMeshParametersProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
      rOStream << "ComputeAveragePfemMeshParametersProcess";
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

  private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Static Member Variables
    ///@{
    ModelPart &mrModelPart;

    MesherUtilities::MeshingParameters &mrRemesh;

    MesherUtilities mMesherUtilities;

    int mEchoLevel;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{

    ///@}/*  */
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ComputeAveragePfemMeshParametersProcess &operator=(ComputeAveragePfemMeshParametersProcess const &rOther);

    /// this function is a private function

    /// Copy constructor.
    // Process(Process const& rOther);

    ///@}

  }; // Class Process

  ///@}

  ///@name Type Definitions
  ///@{

  ///@}
  ///@name Input and output
  ///@{

  /// input stream function
  inline std::istream &operator>>(std::istream &rIStream,
                                  ComputeAveragePfemMeshParametersProcess &rThis);

  /// output stream function
  inline std::ostream &operator<<(std::ostream &rOStream,
                                  const ComputeAveragePfemMeshParametersProcess &rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}

} // namespace Kratos.

#endif // KRATOS_SET_INLET_PROCESS_H_INCLUDED  defined
