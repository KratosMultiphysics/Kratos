//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:                Navaneeth K Narayanan $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                August 2018 $
//   Revision:            $Revision:                    0.0 $
//
//

#ifndef BUILD_MASTER_CONDTIONS_FOR_HYDROSTATIC_PROCESS_H
#define BUILD_MASTER_CONDTIONS_FOR_HYDROSTATIC_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes
#include <boost/timer.hpp>

// Project includes
#include "includes/model_part.h"

#include "includes/condition.h"
#include "processes/process.h"
#include "structural_mechanics_application_variables.h"

///VARIABLES used:
//Data:      MASTER_CONDITION(set)
//StepData:
//Flags:    (checked)
//          (set)     CONTACT(set)
//          (modified)
//          (reset)
// (set):=(set in this process)

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
   */
class BuildMasterConditionsForHydrostaticLoadingProcess
    : public Process
{
public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of BuildMasterConditionsForHydrostaticLoadingProcess
  KRATOS_CLASS_POINTER_DEFINITION(BuildMasterConditionsForHydrostaticLoadingProcess);

  typedef ModelPart::ConditionType ConditionType;
  typedef ModelPart::PropertiesType PropertiesType;
  typedef ConditionType::GeometryType GeometryType;

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  BuildMasterConditionsForHydrostaticLoadingProcess(ModelPart &rModelPart,Element &rNodalELement)
      : mrModelPart(rModelPart),mrNodalElement(rNodalELement)

  {
  }

  /// Destructor.
  virtual ~BuildMasterConditionsForHydrostaticLoadingProcess()
  {
  }

  ///@}
  ///@name Operators
  ///@{

  void operator()()
  {
    Execute();
  }

  ///@}
  ///@name Operations
  ///@{

  void Execute() override
  {
    KRATOS_TRY;
    

    // Master conditions and nodes for the created nodal element
    for (ModelPart::ConditionsContainerType::iterator it_cond = mrModelPart.ConditionsBegin(); it_cond != mrModelPart.ConditionsEnd(); it_cond++)

    {

      mrNodalElement.GetValue(MASTER_CONDITIONS).push_back(Condition::WeakPointer(*(it_cond).base()));
    }
    std::cout << "Number of master conditons ::" << mrNodalElement.GetValue(MASTER_CONDITIONS).size() << std::endl;

    for (ModelPart::NodesContainerType::iterator it_node = mrModelPart.NodesBegin(); it_node != mrModelPart.NodesEnd(); it_node++)

    {

      mrNodalElement.GetValue(MASTER_NODES).push_back(Node<3>::WeakPointer(*(it_node).base()));
    }
    std::cout << "Number of master nodes ::" << mrNodalElement.GetValue(MASTER_NODES).size() << std::endl;

    KRATOS_CATCH(" ");
  }

  ///@}
  ///@name Access
  ///@{

  ///@}u
  ///@name Inquiry
  ///@{

  ///@}
  ///@name Input and output
  ///@{

  /// Turn back information as a string.
  std::string Info() const override
  {
    return "BuildMasterConditionsForHydrostaticLoadingProcess";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream &rOStream) const override
  {
    rOStream << "BuildMasterConditionsForHydrostaticLoadingProcess";
  }

  /// Print object's data.
  void PrintData(std::ostream &rOStream) const override
  {
  }

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
  
  ModelPart &mrModelPart;
  Element &mrNodalElement;

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
  BuildMasterConditionsForHydrostaticLoadingProcess &operator=(BuildMasterConditionsForHydrostaticLoadingProcess const &rOther);

  /// Copy constructor.
  //BuildMasterConditionsForHydrostaticLoadingProcess(BuildMasterConditionsForHydrostaticLoadingProcess const& rOther);

  ///@}

}; // Class BuildMasterConditionsForHydrostaticLoadingProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                BuildMasterConditionsForHydrostaticLoadingProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const BuildMasterConditionsForHydrostaticLoadingProcess &rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_BUILD_CONTACT_CONDITIONS_PROCESS_H_INCLUDED  defined
