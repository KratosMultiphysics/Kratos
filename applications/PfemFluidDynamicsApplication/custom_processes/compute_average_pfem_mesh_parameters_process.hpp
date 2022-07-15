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
      const unsigned int numberOfRefiningBoxes = mrRemesh.UseRefiningBoxList.size();
      array_1d<double, 1> refinedFluidNodesList = ZeroVector(1); // to change!
      array_1d<double, 1> refinedMeanNodalSizeList = ZeroVector(1);
      refinedFluidNodesList.resize(numberOfRefiningBoxes);
      refinedMeanNodalSizeList.resize(numberOfRefiningBoxes);
      refinedFluidNodesList = ZeroVector(numberOfRefiningBoxes);
      refinedMeanNodalSizeList = ZeroVector(numberOfRefiningBoxes);
      double fluidNodes = 0;
      double meanNodalSize = 0;
      for (ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin(); i_node != mrModelPart.NodesEnd(); i_node++)
      {
        if (i_node->Is(FLUID))
        {
          if (numberOfRefiningBoxes == 0)
          {
            fluidNodes += 1.0;
            meanNodalSize += i_node->FastGetSolutionStepValue(NODAL_H);
          }
          else
          {
            unsigned outOfRefiningBoxes = true;
            double transitionElementsInInputMesh = 10.0;
            for (unsigned int index = 0; index < numberOfRefiningBoxes; index++)
            {
              double transitionDistanceInInputMesh = transitionElementsInInputMesh * mrRemesh.RefiningBoxMeshSizeList[index];
              array_1d<double, 3> RefiningBoxMinimumPointList = mrRemesh.RefiningBoxMinimumPointList[index];
              array_1d<double, 3> RefiningBoxMaximumPointList = mrRemesh.RefiningBoxMaximumPointList[index];
              if (dimension == 2)
              {
                if (i_node->X() > (RefiningBoxMinimumPointList[0] - transitionDistanceInInputMesh) && i_node->Y() > (RefiningBoxMinimumPointList[1] - transitionDistanceInInputMesh) &&
                    i_node->X() < (RefiningBoxMaximumPointList[0] + transitionDistanceInInputMesh) && i_node->Y() < (RefiningBoxMaximumPointList[1] + transitionDistanceInInputMesh))
                {
                  outOfRefiningBoxes = false;
                  // break;
                  double localSize = i_node->FastGetSolutionStepValue(NODAL_H) * 2.0; // local size is used to avoid the frontiere zone in the mean mesh size computation
                  if (i_node->X() > (RefiningBoxMinimumPointList[0] + localSize) && i_node->Y() > (RefiningBoxMinimumPointList[1] + localSize) &&
                      i_node->X() < (RefiningBoxMaximumPointList[0] - localSize) && i_node->Y() < (RefiningBoxMaximumPointList[1] - localSize))
                  {
                    refinedFluidNodesList[index] += 1.0;
                    refinedMeanNodalSizeList[index] += i_node->FastGetSolutionStepValue(NODAL_H);
                  }
                }
              }
              else if (dimension == 3)
              {
                if (i_node->X() > (RefiningBoxMinimumPointList[0] - transitionDistanceInInputMesh) && i_node->Y() > (RefiningBoxMinimumPointList[1] - transitionDistanceInInputMesh) && i_node->Z() > (RefiningBoxMinimumPointList[2] - transitionDistanceInInputMesh) &&
                    i_node->X() < (RefiningBoxMaximumPointList[0] + transitionDistanceInInputMesh) && i_node->Y() < (RefiningBoxMaximumPointList[1] + transitionDistanceInInputMesh) && i_node->Z() < (RefiningBoxMaximumPointList[2] + transitionDistanceInInputMesh))
                {
                  outOfRefiningBoxes = false;
                  // break;
                  double localSize = i_node->FastGetSolutionStepValue(NODAL_H) * 2.0; // local size is used to avoid the frontiere zone in the mean mesh size computation
                  if (i_node->X() > (RefiningBoxMinimumPointList[0] + localSize) && i_node->Y() > (RefiningBoxMinimumPointList[1] + localSize) && i_node->Z() > (RefiningBoxMinimumPointList[2] + localSize) &&
                      i_node->X() < (RefiningBoxMaximumPointList[0] - localSize) && i_node->Y() < (RefiningBoxMaximumPointList[1] - localSize) && i_node->Z() < (RefiningBoxMaximumPointList[2] - localSize))
                  {
                    refinedFluidNodesList[index] += 1.0;
                    refinedMeanNodalSizeList[index] += i_node->FastGetSolutionStepValue(NODAL_H);
                  }
                }
              }
            }
            // CONSIDER ONLY THE NODES OUT FROM THE REFINEMENT AREAS
            if (outOfRefiningBoxes == true)
            {
              fluidNodes += 1.0;
              meanNodalSize += i_node->FastGetSolutionStepValue(NODAL_H);
            }
          }
        }
      }
      meanNodalSize *= 1.0 / fluidNodes;

      std::cout << "new computation of meanNodalSize " << meanNodalSize << std::endl;
      std::cout << "new computation of fluidNodes " << fluidNodes << std::endl;

      mrRemesh.Refine->CriticalRadius = meanNodalSize;
      mrRemesh.Refine->InitialRadius = meanNodalSize;

      for (unsigned int index = 0; index < numberOfRefiningBoxes; index++)
      {
        // {
        //   mrRemesh.RefiningBoxMeshSizeList[index] *= 0.8;
        // }
        std::cout << " user provided fine mesh size " << mrRemesh.RefiningBoxMeshSizeList[index] << std::endl;
        double localMeshSize = refinedMeanNodalSizeList[index] / refinedFluidNodesList[index];
        mrRemesh.SetRefiningBoxMeshSizeList(index, localMeshSize);
        std::cout << " fluid nodesin refined zone " << refinedFluidNodesList[index] << std::endl;
        std::cout << " new computed fine mesh size " << mrRemesh.RefiningBoxMeshSizeList[index] << std::endl;
        double tolerance = mrRemesh.RefiningBoxMeshSizeList[index] * 0.01;
        double differenceOfSize = meanNodalSize - mrRemesh.RefiningBoxMeshSizeList[index];

        mrRemesh.RefiningBoxMinimumPointList[index][0] += -tolerance; // the finest nodes at the frontier should not be erased
        mrRemesh.RefiningBoxMinimumPointList[index][1] += -tolerance;
        mrRemesh.RefiningBoxMinimumPointList[index][2] += -tolerance;

        mrRemesh.RefiningBoxMaximumPointList[index][0] += tolerance;
        mrRemesh.RefiningBoxMaximumPointList[index][1] += tolerance;
        mrRemesh.RefiningBoxMaximumPointList[index][2] += tolerance;

        double transitionDistance = 5.0 * fabs(differenceOfSize);

        mrRemesh.RefiningBoxMinExternalPointList[index][0] = mrRemesh.RefiningBoxMinimumPointList[index][0] - transitionDistance;
        mrRemesh.RefiningBoxMinExternalPointList[index][1] = mrRemesh.RefiningBoxMinimumPointList[index][1] - transitionDistance;
        mrRemesh.RefiningBoxMinExternalPointList[index][2] = mrRemesh.RefiningBoxMinimumPointList[index][2] - transitionDistance;
        // mrRemesh.RefiningBoxMinInternalPointList[index][0] = mrRemesh.RefiningBoxMinimumPointList[index][0] + mrRemesh.RefiningBoxMeshSizeList[index];
        // mrRemesh.RefiningBoxMinInternalPointList[index][1] = mrRemesh.RefiningBoxMinimumPointList[index][1] + mrRemesh.RefiningBoxMeshSizeList[index];
        // mrRemesh.RefiningBoxMinInternalPointList[index][2] = mrRemesh.RefiningBoxMinimumPointList[index][2] + mrRemesh.RefiningBoxMeshSizeList[index];

        mrRemesh.RefiningBoxMaxExternalPointList[index][0] = mrRemesh.RefiningBoxMaximumPointList[index][0] + transitionDistance;
        mrRemesh.RefiningBoxMaxExternalPointList[index][1] = mrRemesh.RefiningBoxMaximumPointList[index][1] + transitionDistance;
        mrRemesh.RefiningBoxMaxExternalPointList[index][2] = mrRemesh.RefiningBoxMaximumPointList[index][2] + transitionDistance;
        // mrRemesh.RefiningBoxMaxInternalPointList[index][0] = mrRemesh.RefiningBoxMaximumPointList[index][0] - mrRemesh.RefiningBoxMeshSizeList[index];
        // mrRemesh.RefiningBoxMaxInternalPointList[index][1] = mrRemesh.RefiningBoxMaximumPointList[index][1] - mrRemesh.RefiningBoxMeshSizeList[index];
        // mrRemesh.RefiningBoxMaxInternalPointList[index][2] = mrRemesh.RefiningBoxMaximumPointList[index][2] - mrRemesh.RefiningBoxMeshSizeList[index];

        std::cout << "mrRemesh.RefiningBoxMaxExternalPointList[index][0] " << mrRemesh.RefiningBoxMaxExternalPointList[index][0] << std::endl;
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
