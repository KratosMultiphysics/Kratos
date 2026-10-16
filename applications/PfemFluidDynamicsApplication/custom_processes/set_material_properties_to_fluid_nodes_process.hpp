//
//   Project Name:        KratosPfemFluidApplication $
//   Created by:          $Author:           AFranci $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:        October 2018 $
//   Revision:            $Revision:             0.0 $
//
//

#if !defined(KRATOS_SET_MATERIAL_PROPERTIES_TO_FLUID_NODES_PROCESS_H_INCLUDED)
#define KRATOS_SET_MATERIAL_PROPERTIES_TO_FLUID_NODES_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes

#include "spatial_containers/spatial_containers.h"

#include "custom_processes/set_material_properties_to_fluid_nodes_process.hpp"
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
  class SetMaterialPropertiesToFluidNodesProcess
      : public MesherProcess
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SetMaterialPropertiesToFluidNodesProcess
    KRATOS_CLASS_POINTER_DEFINITION(SetMaterialPropertiesToFluidNodesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SetMaterialPropertiesToFluidNodesProcess(ModelPart &rModelPart)
        : mrModelPart(rModelPart)
    {
    }

    /// Destructor.
    virtual ~SetMaterialPropertiesToFluidNodesProcess()
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

    const int NoPropertyId = std::numeric_limits<int>::max();

    const auto& rNodalVariables =
        mrModelPart.GetNodalSolutionStepVariablesList();

    // =========================================================================
    // Check that PROPERTY_ID is available
    // =========================================================================

    KRATOS_ERROR_IF_NOT(rNodalVariables.Has(PROPERTY_ID))
        << "PROPERTY_ID must be included in the "
        << "NodalSolutionStepVariablesList."
        << std::endl;


    // =========================================================================
    // 1. Initialise PROPERTY_ID for all nodes
    // =========================================================================

    for (ModelPart::NodeIterator itNode = mrModelPart.NodesBegin();
         itNode != mrModelPart.NodesEnd();
         ++itNode)
    {
        itNode->FastGetSolutionStepValue(PROPERTY_ID) = NoPropertyId;
    }


    // =========================================================================
    // 2. Determine the material associated with every node
    //
    // IMPORTANT:
    //
    // If a node belongs to several fluid elements with different materials,
    // the material with the smallest Properties ID is selected.
    //
    // This part is intentionally SERIAL to guarantee deterministic results.
    // =========================================================================

    for (ModelPart::ElementIterator itElem = mrModelPart.ElementsBegin();
         itElem != mrModelPart.ElementsEnd();
         ++itElem)
    {
        const int ElemPropertyId =
            static_cast<int>(itElem->GetProperties().Id());

        Geometry<Node>& rGeometry =
            itElem->GetGeometry();

        const SizeType NumNodes =
            rGeometry.PointsNumber();

        for (SizeType i = 0; i < NumNodes; ++i)
        {
            Node& rNode = rGeometry[i];

            int& rNodePropertyId =
                rNode.FastGetSolutionStepValue(PROPERTY_ID);

            if (ElemPropertyId < rNodePropertyId)
            {
                rNodePropertyId = ElemPropertyId;
            }
        }
    }

    // =========================================================================
    // 3. Copy the selected material properties to the nodes
    // =========================================================================

#pragma omp parallel
    {
        ModelPart::NodeIterator NodeBegin;
        ModelPart::NodeIterator NodeEnd;

        OpenMPUtils::PartitionedIterators(
            mrModelPart.Nodes(),
            NodeBegin,
            NodeEnd);

        for (ModelPart::NodeIterator itNode = NodeBegin;
             itNode != NodeEnd;
             ++itNode)
        {
            Node& rNode = *itNode;

            const int PropertyId =
                rNode.FastGetSolutionStepValue(PROPERTY_ID);

            if (PropertyId == NoPropertyId)
            {
                continue;
            }

            Properties& rProperties =
                mrModelPart.GetProperties(
                    static_cast<unsigned int>(PropertyId));


            // =============================================================
            // GENERAL MATERIAL PROPERTIES
            // =============================================================

            if (rNodalVariables.Has(DENSITY) &&
                rProperties.Has(DENSITY))
            {
                rNode.FastGetSolutionStepValue(DENSITY) =
                    rProperties[DENSITY];
            }

            if (rNodalVariables.Has(BULK_MODULUS) &&
                rProperties.Has(BULK_MODULUS))
            {
                rNode.FastGetSolutionStepValue(BULK_MODULUS) =
                    rProperties[BULK_MODULUS];
            }

            if (rNodalVariables.Has(DYNAMIC_VISCOSITY) &&
                rProperties.Has(DYNAMIC_VISCOSITY))
            {
                rNode.FastGetSolutionStepValue(DYNAMIC_VISCOSITY) =
                    rProperties[DYNAMIC_VISCOSITY];
            }


            // =============================================================
            // BINGHAM / HERSCHEL-BULKLEY VARIABLES
            // =============================================================

            if (rNodalVariables.Has(YIELD_SHEAR) &&
                rProperties.Has(YIELD_SHEAR))
            {
                rNode.FastGetSolutionStepValue(YIELD_SHEAR) =
                    rProperties[YIELD_SHEAR];
            }

            if (rNodalVariables.Has(FLOW_INDEX) &&
                rProperties.Has(FLOW_INDEX))
            {
                rNode.FastGetSolutionStepValue(FLOW_INDEX) =
                    rProperties[FLOW_INDEX];
            }

            if (rNodalVariables.Has(ADAPTIVE_EXPONENT) &&
                rProperties.Has(ADAPTIVE_EXPONENT))
            {
                rNode.FastGetSolutionStepValue(ADAPTIVE_EXPONENT) =
                    rProperties[ADAPTIVE_EXPONENT];
            }


            // =============================================================
            // FRICTIONAL VISCOPLASTIC MODEL
            // =============================================================

            if (rNodalVariables.Has(INTERNAL_FRICTION_ANGLE) &&
                rProperties.Has(INTERNAL_FRICTION_ANGLE))
            {
                rNode.FastGetSolutionStepValue(INTERNAL_FRICTION_ANGLE) =
                    rProperties[INTERNAL_FRICTION_ANGLE];
            }

            if (rNodalVariables.Has(COHESION) &&
                rProperties.Has(COHESION))
            {
                rNode.FastGetSolutionStepValue(COHESION) =
                    rProperties[COHESION];
            }


            // =============================================================
            // MU(I) RHEOLOGY
            // =============================================================

            if (rNodalVariables.Has(STATIC_FRICTION) &&
                rProperties.Has(STATIC_FRICTION))
            {
                rNode.FastGetSolutionStepValue(STATIC_FRICTION) =
                    rProperties[STATIC_FRICTION];
            }

            if (rNodalVariables.Has(DYNAMIC_FRICTION) &&
                rProperties.Has(DYNAMIC_FRICTION))
            {
                rNode.FastGetSolutionStepValue(DYNAMIC_FRICTION) =
                    rProperties[DYNAMIC_FRICTION];
            }

            if (rNodalVariables.Has(INERTIAL_NUMBER_ZERO) &&
                rProperties.Has(INERTIAL_NUMBER_ZERO))
            {
                rNode.FastGetSolutionStepValue(INERTIAL_NUMBER_ZERO) =
                    rProperties[INERTIAL_NUMBER_ZERO];
            }

            if (rNodalVariables.Has(GRAIN_DIAMETER) &&
                rProperties.Has(GRAIN_DIAMETER))
            {
                rNode.FastGetSolutionStepValue(GRAIN_DIAMETER) =
                    rProperties[GRAIN_DIAMETER];
            }

            if (rNodalVariables.Has(GRAIN_DENSITY) &&
                rProperties.Has(GRAIN_DENSITY))
            {
                rNode.FastGetSolutionStepValue(GRAIN_DENSITY) =
                    rProperties[GRAIN_DENSITY];
            }

            if (rNodalVariables.Has(REGULARIZATION_COEFFICIENT) &&
                rProperties.Has(REGULARIZATION_COEFFICIENT))
            {
                rNode.FastGetSolutionStepValue(REGULARIZATION_COEFFICIENT) =
                    rProperties[REGULARIZATION_COEFFICIENT];
            }
        }
    }

    KRATOS_CATCH("")
}

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
  return "SetMaterialPropertiesToFluidNodesProcess";
}

/// Print information about this object.
void PrintInfo(std::ostream &rOStream) const override
{
  rOStream << "SetMaterialPropertiesToFluidNodesProcess";
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
SetMaterialPropertiesToFluidNodesProcess &operator=(SetMaterialPropertiesToFluidNodesProcess const &rOther);

/// Copy constructor.
// SetMaterialPropertiesToFluidNodesProcess(SetMaterialPropertiesToFluidNodesProcess const& rOther);

///@}
}
; // Class SetMaterialPropertiesToFluidNodesProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                SetMaterialPropertiesToFluidNodesProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const SetMaterialPropertiesToFluidNodesProcess &rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_SET_MATERIAL_PROPERTIES_TO_FLUID_NODES_PROCESS_H_INCLUDED  defined
