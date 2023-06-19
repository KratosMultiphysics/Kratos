//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:                   AFranci $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:                October 2019 $
//   Revision:            $Revision:                     0.0 $
//
//

#if !defined(KRATOS_FIND_NODAL_H_FOR_RIGID_WALLS_PROCESS_H_INCLUDED)
#define KRATOS_FIND_NODAL_H_FOR_RIGID_WALLS_PROCESS_H_INCLUDED

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

  class FindNodalHForRigidWallsProcess
      : public Process
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Process
    KRATOS_CLASS_POINTER_DEFINITION(FindNodalHForRigidWallsProcess);

    typedef ModelPart::NodeType NodeType;
    typedef ModelPart::ConditionType ConditionType;
    typedef ModelPart::PropertiesType PropertiesType;
    typedef ConditionType::GeometryType GeometryType;
    typedef std::size_t SizeType;
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FindNodalHForRigidWallsProcess(ModelPart &rModelPart)
        : mrModelPart(rModelPart)
    {
      KRATOS_INFO("FindNodalHForRigidWallsProcess") << " activated " << std::endl;
    }

    /// Destructor.
    virtual ~FindNodalHForRigidWallsProcess() {}

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

    void Execute() override
    {
      KRATOS_TRY
      // Check if variables are available

      // initialize to zero
      for (ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin(); i_node != mrModelPart.NodesEnd(); i_node++)
      {
        if (i_node->Is(RIGID))
        {
          i_node->FastGetSolutionStepValue(NODAL_H_WALL) = 0;
        }
      }

      for (IndexType i = 0; i < mrModelPart.Elements().size(); ++i)
      {
        auto it_element = mrModelPart.ElementsBegin() + i;
        auto &r_geom = it_element->GetGeometry();
        const SizeType number_of_nodes = r_geom.size();

        for (IndexType k = 0; k < number_of_nodes - 1; ++k)
        {
          if (r_geom[k].Is(RIGID))
          {
            double &r_h1 = r_geom[k].FastGetSolutionStepValue(NODAL_H_WALL);
            for (IndexType l = k + 1; l < number_of_nodes; ++l)
            {
              if (r_geom[l].Is(RIGID))
              {
                double hedge = norm_2(r_geom[l].Coordinates() - r_geom[k].Coordinates());
                double &r_h2 = r_geom[l].FastGetSolutionStepValue(NODAL_H_WALL);
                // Get minimum between the existent value and the considered edge length
                r_geom[k].FastGetSolutionStepValue(NODAL_H_WALL) = std::max(r_h1, hedge);
                r_geom[l].FastGetSolutionStepValue(NODAL_H_WALL) = std::max(r_h2, hedge);
              }
            }
          }
        }
      }
      // Synchronize between processes
      // mrModelPart.GetCommunicator().SynchronizeCurrentDataToMin(NODAL_H_WALL);
      KRATOS_CATCH("")
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
      return "FindNodalHForRigidWallsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
      rOStream << "FindNodalHForRigidWallsProcess";
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
    FindNodalHForRigidWallsProcess &operator=(FindNodalHForRigidWallsProcess const &rOther);

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
                                  FindNodalHForRigidWallsProcess &rThis);

  /// output stream function
  inline std::ostream &operator<<(std::ostream &rOStream,
                                  const FindNodalHForRigidWallsProcess &rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}

} // namespace Kratos.

#endif // KRATOS_SET_INLET_PROCESS_H_INCLUDED  defined
