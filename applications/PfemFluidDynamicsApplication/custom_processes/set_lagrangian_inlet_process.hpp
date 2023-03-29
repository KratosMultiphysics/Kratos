//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:                   AFranci $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:                October 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

#if !defined(KRATOS_SET_LAGRANGIAN_INLET_PROCESS_H_INCLUDED)
#define KRATOS_SET_LAGRANGIAN_INLET_PROCESS_H_INCLUDED

// External includes

// System includes

// Project includes
#include "containers/variables_list_data_value_container.h"
#include "spatial_containers/spatial_containers.h"

#include "includes/model_part.h"
#include "custom_utilities/mesh_error_calculation_utilities.hpp"
#include "custom_utilities/mesher_utilities.hpp"

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

  class SetLagrangianInletProcess
      : public Process
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(SetLagrangianInletProcess);

    typedef ModelPart::NodeType NodeType;
    typedef ModelPart::ConditionType ConditionType;
    typedef ModelPart::PropertiesType PropertiesType;
    typedef ConditionType::GeometryType GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SetLagrangianInletProcess(ModelPart &rModelPart,
                              int EchoLevel)
        : mrModelPart(rModelPart)
    {
      KRATOS_INFO("SetLagrangianInletProcess") << " activated " << std::endl;

      mEchoLevel = EchoLevel;
    }

    /// Destructor.
    virtual ~SetLagrangianInletProcess() {}

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
        std::cout << "  SET LAGRANGIAN INLET PROCESS ]; " << std::endl;

      for (ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin(); i_node != mrModelPart.NodesEnd(); i_node++)
      {
        i_node->Set(INLET);
        i_node->Set(RIGID);
        i_node->Set(PFEMFlags::LAGRANGIAN_INLET);
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
      return "SetLagrangianInletProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
      rOStream << "SetLagrangianInletProcess";
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

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    SetLagrangianInletProcess &operator=(SetLagrangianInletProcess const &rOther);

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
                                  SetLagrangianInletProcess &rThis);

  /// output stream function
  inline std::ostream &operator<<(std::ostream &rOStream,
                                  const SetLagrangianInletProcess &rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}

} // namespace Kratos.

#endif // KRATOS_SET_LAGRANGIAN_INLET_PROCESS_H_INCLUDED  defined
