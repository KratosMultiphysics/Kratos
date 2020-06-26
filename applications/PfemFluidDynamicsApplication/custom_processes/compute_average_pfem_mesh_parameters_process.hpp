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

///VARIABLES used:
//Data:
//StepData: CONTACT_FORCE, DISPLACEMENT
//Flags:    (checked)
//          (set)
//          (modified)
//          (reset)
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
    KRATOS_INFO("ComputeAveragePfemMeshParametersProcess") << " inlet_management CONSTRUCTOR ";

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

    bool refiningBox = mrRemesh.UseRefiningBox;

    array_1d<double, 3> &minExternalPointRefiningBox = mrRemesh.RefiningBoxMinExternalPoint;
    array_1d<double, 3> &minInternalPointRefiningBox = mrRemesh.RefiningBoxMinInternalPoint;
    array_1d<double, 3> &maxExternalPointRefiningBox = mrRemesh.RefiningBoxMaxExternalPoint;
    array_1d<double, 3> &maxInternalPointRefiningBox = mrRemesh.RefiningBoxMaxInternalPoint;
    array_1d<double, 3> &RefiningBoxMinimumPoint = mrRemesh.RefiningBoxMinimumPoint;
    array_1d<double, 3> &RefiningBoxMaximumPoint = mrRemesh.RefiningBoxMaximumPoint;

    double fluidNodes = 0;
    double meanNodalSize = 0;

    // double refinedFluidNodes = 0;
    // double refinedMeanNodalSize = 0;

    const unsigned int dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
    // unsigned int count=0;
    for (ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin(); i_node != mrModelPart.NodesEnd(); i_node++)
    {
      if (refiningBox == false)
      {
        if (i_node->Is(FLUID))
        {
          fluidNodes += 1.0;
          meanNodalSize += i_node->FastGetSolutionStepValue(NODAL_H);
        }
      }
      else
      {
        if (dimension == 2)
        {
          if (i_node->X() < RefiningBoxMinimumPoint[0] || i_node->Y() < RefiningBoxMinimumPoint[1] ||
              i_node->X() > RefiningBoxMaximumPoint[0] || i_node->Y() > RefiningBoxMaximumPoint[1])
          {
            //CONSIDER ONLY THE NODES OUT FROM THE REFINEMENT AREA
            if (i_node->Is(FLUID))
            {
              fluidNodes += 1.0;
              meanNodalSize += i_node->FastGetSolutionStepValue(NODAL_H);
            }
          }
          // else{
          //   if (i_node->Is(FLUID))
          //   {
          //     refinedFluidNodes += 1.0;
          //     refinedMeanNodalSize += i_node->FastGetSolutionStepValue(NODAL_H);
          //   }
          // }
        }
        else if (dimension == 3)
        {
          if (i_node->X() < RefiningBoxMinimumPoint[0] || i_node->Y() < RefiningBoxMinimumPoint[1] || i_node->Z() < RefiningBoxMinimumPoint[2] ||
              i_node->X() > RefiningBoxMaximumPoint[0] || i_node->Y() > RefiningBoxMaximumPoint[1] || i_node->Z() > RefiningBoxMaximumPoint[2])
          {
            //CONSIDER ONLY THE NODES OUT FROM THE REFINEMENT AREA
            if (i_node->Is(FLUID))
            {
              fluidNodes += 1.0;
              meanNodalSize += i_node->FastGetSolutionStepValue(NODAL_H);
            }
          }
          // else
          // {
          //   if (i_node->Is(FLUID))
          //   {
          //     refinedFluidNodes += 1.0;
          //     refinedMeanNodalSize += i_node->FastGetSolutionStepValue(NODAL_H);
          //   }
          // }
        }
      }
    }
    meanNodalSize *= 1.0 / fluidNodes;
    // refinedMeanNodalSize *= 1.0 / refinedFluidNodes;

    mrRemesh.Refine->CriticalRadius = meanNodalSize;
    mrRemesh.Refine->InitialRadius = meanNodalSize;

    if (dimension == 3)
    {
      mrRemesh.RefiningBoxMeshSize *= 0.8;
    }

    // std::cout << fluidNodes << " nodes in Not Refined area with mean element size: " << mrRemesh.Refine->CriticalRadius << std::endl;
    // std::cout << refinedFluidNodes << " nodes in Refined area with mean element size:  " << mrRemesh.RefiningBoxMeshSize << std::endl;
    // std::cout << " othermeanNodalSize " << refinedMeanNodalSize << std::endl;

    double smallSize = meanNodalSize;

    if (meanNodalSize < mrRemesh.RefiningBoxMeshSize)
    {
      smallSize = mrRemesh.RefiningBoxMeshSize;

      RefiningBoxMinimumPoint[0] += 0.01 * smallSize; //the finest nodes at the frontier should not be erased
      RefiningBoxMinimumPoint[1] += 0.01 * smallSize;
      RefiningBoxMinimumPoint[2] += 0.01 * smallSize;

      RefiningBoxMaximumPoint[0] += -0.01 * smallSize;
      RefiningBoxMaximumPoint[1] += -0.01 * smallSize;
      RefiningBoxMaximumPoint[2] += -0.01 * smallSize;
    }
    else // the mesh is finer in the RefiningBox
    {
      RefiningBoxMinimumPoint[0] += -0.01 * smallSize; //the finest nodes at the frontier should not be erased
      RefiningBoxMinimumPoint[1] += -0.01 * smallSize;
      RefiningBoxMinimumPoint[2] += -0.01 * smallSize;

      RefiningBoxMaximumPoint[0] += 0.01 * smallSize;
      RefiningBoxMaximumPoint[1] += 0.01 * smallSize;
      RefiningBoxMaximumPoint[2] += 0.01 * smallSize;
    }

    minExternalPointRefiningBox[0] = RefiningBoxMinimumPoint[0] - mrRemesh.Refine->CriticalRadius;
    minExternalPointRefiningBox[1] = RefiningBoxMinimumPoint[1] - mrRemesh.Refine->CriticalRadius;
    minExternalPointRefiningBox[2] = RefiningBoxMinimumPoint[2] - mrRemesh.Refine->CriticalRadius;
    minInternalPointRefiningBox[0] = RefiningBoxMinimumPoint[0] + mrRemesh.RefiningBoxMeshSize;
    minInternalPointRefiningBox[1] = RefiningBoxMinimumPoint[1] + mrRemesh.RefiningBoxMeshSize;
    minInternalPointRefiningBox[2] = RefiningBoxMinimumPoint[2] + mrRemesh.RefiningBoxMeshSize;

    maxExternalPointRefiningBox[0] = RefiningBoxMaximumPoint[0] + mrRemesh.Refine->CriticalRadius;
    maxExternalPointRefiningBox[1] = RefiningBoxMaximumPoint[1] + mrRemesh.Refine->CriticalRadius;
    maxExternalPointRefiningBox[2] = RefiningBoxMaximumPoint[2] + mrRemesh.Refine->CriticalRadius;
    maxInternalPointRefiningBox[0] = RefiningBoxMaximumPoint[0] - mrRemesh.RefiningBoxMeshSize;
    maxInternalPointRefiningBox[1] = RefiningBoxMaximumPoint[1] - mrRemesh.RefiningBoxMeshSize;
    maxInternalPointRefiningBox[2] = RefiningBoxMaximumPoint[2] - mrRemesh.RefiningBoxMeshSize;

    // std::cout<<" RefiningBoxMinimumPoint "<<mrRemesh.RefiningBoxMinimumPoint <<std::endl;
    // std::cout<<" minExternalPointRefiningBox "<<mrRemesh.RefiningBoxMinExternalPoint <<std::endl;
    // std::cout<<" minInternalPointRefiningBox "<<mrRemesh.RefiningBoxMinInternalPoint <<std::endl;
    // std::cout<<"     RefiningBoxMaximumPoint "<<mrRemesh.RefiningBoxMaximumPoint <<std::endl;
    // std::cout<<"     maxExternalPointRefiningBox "<<mrRemesh.RefiningBoxMaxExternalPoint <<std::endl;
    // std::cout<<"     maxInternalPointRefiningBox "<<mrRemesh.RefiningBoxMaxInternalPoint <<std::endl;

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
  //Process(Process const& rOther);

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
