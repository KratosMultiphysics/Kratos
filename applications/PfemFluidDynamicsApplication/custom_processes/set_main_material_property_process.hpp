//
//   Project Name:        KratosPfemFluidApplication $
//   Created by:          $Author:           AFranci $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:          April 2022 $
//   Revision:            $Revision:             0.0 $
//
//

#if !defined(KRATOS_SET_MAIN_MATERIAL_PROPERTY_PROCESS_H_INCLUDED)
#define KRATOS_SET_MAIN_MATERIAL_PROPERTY_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes

#include "spatial_containers/spatial_containers.h"

#include "custom_processes/set_main_material_property_process.hpp"
#include "custom_utilities/mesher_utilities.hpp"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "utilities/math_utils.h"
#include "custom_processes/mesher_process.hpp"

/// VARIABLES used:
// Data:
// StepData:
// Flags:    (checked)
//           (set)
//           (modified)
//           (reset)

namespace Kratos
{

  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{
  typedef ModelPart::NodesContainerType NodesContainerType;
  typedef ModelPart::ElementsContainerType ElementsContainerType;
  typedef ModelPart::MeshType::GeometryType::PointsArrayType PointsArrayType;
  typedef std::size_t SizeType;

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
  class SetMainMaterialPropertyProcess
      : public MesherProcess
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SetMainMaterialPropertyProcess
    KRATOS_CLASS_POINTER_DEFINITION(SetMainMaterialPropertyProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SetMainMaterialPropertyProcess(ModelPart &rModelPart)
        : mrModelPart(rModelPart)
    {
    }

    /// Destructor.
    virtual ~SetMainMaterialPropertyProcess()
    {
    }

    void operator()()
    {
      Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    void Execute() override
    {
      KRATOS_TRY
      const auto elem_begin = mrModelPart.ElementsBegin();
      unsigned int main_property_id = elem_begin->GetProperties().Id();
      ProcessInfo &rCurrentProcessInfo = mrModelPart.GetProcessInfo();
      rCurrentProcessInfo[MAIN_MATERIAL_PROPERTY] = main_property_id;
      KRATOS_CATCH(" ")
    }; // namespace Kratos

    ///@}
    ///@name Operators
    ///@{

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
      return "SetMainMaterialPropertyProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
      rOStream << "SetMainMaterialPropertyProcess";
    }

  protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{
    ModelPart &mrModelPart;

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
    SetMainMaterialPropertyProcess &operator=(SetMainMaterialPropertyProcess const &rOther);

    /// Copy constructor.
    // SetMainMaterialPropertyProcess(SetMainMaterialPropertyProcess const& rOther);

    ///@}
  }; // Class SetMainMaterialPropertyProcess

  ///@}

  ///@name Type Definitions
  ///@{

  ///@}
  ///@name Input and output
  ///@{

  /// input stream function
  inline std::istream &operator>>(std::istream &rIStream,
                                  SetMainMaterialPropertyProcess &rThis);

  /// output stream function
  inline std::ostream &operator<<(std::ostream &rOStream,
                                  const SetMainMaterialPropertyProcess &rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}

} // namespace Kratos.

#endif // KRATOS_SET_MAIN_MATERIAL_PROPERTY_PROCESS_H_INCLUDED  defined
