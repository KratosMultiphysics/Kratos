//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:                   AFranci $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:                October 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

#if !defined(KRATOS_INLET_MANAGEMENT_PROCESS_H_INCLUDED)
#define KRATOS_INLET_MANAGEMENT_PROCESS_H_INCLUDED

// External includes

// System includes

// Project includes
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"

#include "includes/model_part.h"
#include "custom_utilities/mesh_error_calculation_utilities.hpp"
#include "custom_utilities/mesher_utilities.hpp"
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

  class InletManagementProcess
      : public MesherProcess
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(InletManagementProcess);

    typedef ModelPart::NodeType NodeType;
    typedef ModelPart::ConditionType ConditionType;
    typedef ModelPart::PropertiesType PropertiesType;
    typedef ConditionType::GeometryType GeometryType;

    typedef GlobalPointersVector<Node> NodeWeakPtrVectorType;
    typedef GlobalPointersVector<Element> ElementWeakPtrVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    InletManagementProcess(ModelPart &rModelPart,
                           MesherUtilities::MeshingParameters &rRemeshingParameters,
                           int EchoLevel)
        : mrModelPart(rModelPart),
          mrRemesh(rRemeshingParameters)
    {
      KRATOS_INFO("InletManagementProcess") << " activated " << std::endl;

      mEchoLevel = EchoLevel;
    }

    /// Destructor.
    virtual ~InletManagementProcess() {}

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
        std::cout << " [ INLET MANAGEMENT PROCESS: " << std::endl;

      if (mrModelPart.Name() != mrRemesh.SubModelPartName)
        std::cout << " ModelPart Supplied do not corresponds to the Meshing Domain: (" << mrModelPart.Name() << " != " << mrRemesh.SubModelPartName << ")" << std::endl;

      const ProcessInfo &rCurrentProcessInfo = mrModelPart.GetProcessInfo();
      double currentTime = rCurrentProcessInfo[TIME];
      double timeInterval = rCurrentProcessInfo[DELTA_TIME];

      unsigned int numberOfLagrangianInletNodes = mrRemesh.Info->NumberOfLagrangianInletNodes;

      if (currentTime < 2 * timeInterval)
      {
        mrRemesh.Info->RemovedNodes = 0;
        if (mEchoLevel > 1)
          std::cout << " First meshes: I repare the mesh without adding new nodes" << std::endl;
        mrRemesh.Info->InitialNumberOfNodes = mrRemesh.Info->NumberOfNodes;
      }

      if (currentTime > 1.5 * timeInterval && numberOfLagrangianInletNodes > 0)
      {
        CheckAndCreateNewInletLayer();
      }
      else
      {
        CountInletNodes();
      }

      if (mEchoLevel > 1)
        std::cout << "   INLET MANAGEMENT PROCESS ]; " << std::endl;

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
      return "InletManagementProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
      rOStream << "InletManagementProcess";
    }

    /// Print object's data.s
    void PrintData(std::ostream &rOStream) const override
    {
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

    void CheckAndCreateNewInletLayer()

    {
      KRATOS_TRY

      if (mEchoLevel > 1)
        std::cout << " CheckAndCreateNewInletLayer " << std::endl;
      const unsigned int dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
      double maxSeparation = mrRemesh.Refine->CriticalRadius;

      std::vector<Node::Pointer> clonedNodes;
      clonedNodes.clear();
      clonedNodes.resize(1);
      unsigned int numberClonedNodes = 0;
      unsigned int sizeClonedNodes = 0;

      for (ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin(); i_node != mrModelPart.NodesEnd(); i_node++)
      {
        if (i_node->Is(PFEMFlags::LAGRANGIAN_INLET))
        {
          ElementWeakPtrVectorType &neighb_elems = i_node->GetValue(NEIGHBOUR_ELEMENTS);
          NodeWeakPtrVectorType &rN = i_node->GetValue(NEIGHBOUR_NODES);

          if ((neighb_elems.size() == 0 && rN.size() == 0) || i_node->Is(RIGID))
          {

            const array_1d<double, 3> &inletDisplacement = i_node->FastGetSolutionStepValue(DISPLACEMENT);
            double distanceFromOrigin = sqrt(inletDisplacement[0] * inletDisplacement[0] +
                                             inletDisplacement[1] * inletDisplacement[1]);
            if (dimension == 3)
            {
              distanceFromOrigin = sqrt(inletDisplacement[0] * inletDisplacement[0] +
                                        inletDisplacement[1] * inletDisplacement[1] +
                                        inletDisplacement[2] * inletDisplacement[2]);
            }

            if (distanceFromOrigin > maxSeparation)
            {

              if (i_node->Is(FLUID))
              {
                Node::Pointer pnode = i_node->Clone();
                sizeClonedNodes = numberClonedNodes + 1;
                clonedNodes.resize(sizeClonedNodes);
                clonedNodes[numberClonedNodes] = pnode;
                numberClonedNodes++;
              }
              // inlet nodes will be replaced at their initial position
              i_node->X() = i_node->X0();
              i_node->Y() = i_node->Y0();
              i_node->FastGetSolutionStepValue(DISPLACEMENT_X, 0) = 0;
              i_node->FastGetSolutionStepValue(DISPLACEMENT_X, 1) = 0;
              i_node->FastGetSolutionStepValue(DISPLACEMENT_Y, 0) = 0;
              i_node->FastGetSolutionStepValue(DISPLACEMENT_Y, 1) = 0;
              if (dimension == 3)
              {
                i_node->Z() = i_node->Z0();
                i_node->FastGetSolutionStepValue(DISPLACEMENT_Z, 0) = 0;
                i_node->FastGetSolutionStepValue(DISPLACEMENT_Z, 1) = 0;
              }

            } /// if maxSeparation> limit
          }
        }
      }

      for (unsigned int i = 0; i < sizeClonedNodes; i++)
      {

        Node::Pointer pnode = clonedNodes[i];
        double NodeIdParent = MesherUtilities::GetMaxNodeId(mrModelPart.GetParentModelPart());
        double NodeId = MesherUtilities::GetMaxNodeId(mrModelPart);
        unsigned int id = NodeIdParent + 1; // total model part node size

        if (NodeId > NodeIdParent)
        {
          id = NodeId + 1;
        }
        pnode->SetId(id);
        pnode->Free(VELOCITY_X);
        pnode->Free(VELOCITY_Y);
        if (dimension == 3)
        {
          pnode->Free(VELOCITY_Z);
        }
        pnode->Reset(PFEMFlags::LAGRANGIAN_INLET);
        pnode->Reset(INLET);
        pnode->Reset(RIGID);
        pnode->Reset(BOUNDARY);
        mrRemesh.NodalPreIds.push_back(pnode->Id());
        mrModelPart.AddNode(pnode);
      }

      KRATOS_CATCH("")
    }

    void CountInletNodes()

    {
      KRATOS_TRY

      unsigned int eulerianInletNodes = 0;
      unsigned int lagrangianInletNodes = 0;
      const unsigned int dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
      const double tolerance = -1e-12;

      for (ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin(); i_node != mrModelPart.NodesEnd(); i_node++)
      {
        if (i_node->Is(PFEMFlags::EULERIAN_INLET))
        {
          const array_1d<double, 3> &inletVelocity = i_node->FastGetSolutionStepValue(VELOCITY);
          const array_1d<double, 3> &inletNormal = i_node->FastGetSolutionStepValue(NORMAL);
          const double inwardX = inletVelocity[0] * inletNormal[0];
          const double inwardY = inletVelocity[1] * inletNormal[1];
          if (dimension == 2)
          {
            if (inwardX < tolerance || inwardY < tolerance)
            {
              eulerianInletNodes += 1;
            }
          }
          else if (dimension == 3)
          {
            const double inwardZ = inletVelocity[2] * inletNormal[2];
            if (inwardX < tolerance || inwardY < tolerance || inwardZ < tolerance)
            {
              eulerianInletNodes += 1;
            }
          }
        }
        else if (i_node->Is(PFEMFlags::LAGRANGIAN_INLET))
        {
          lagrangianInletNodes += 1;
        }
      }
      mrRemesh.Info->NumberOfEulerianInletNodes = eulerianInletNodes;
      mrRemesh.Info->NumberOfLagrangianInletNodes = lagrangianInletNodes;

      KRATOS_CATCH("")
    }

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
    InletManagementProcess &
    operator=(InletManagementProcess const &rOther);

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
                                  InletManagementProcess &rThis);

  /// output stream function
  inline std::ostream &operator<<(std::ostream &rOStream,
                                  const InletManagementProcess &rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}

} // namespace Kratos.

#endif // KRATOS_INLET_MANAGEMENT_PROCESS_H_INCLUDED  defined
