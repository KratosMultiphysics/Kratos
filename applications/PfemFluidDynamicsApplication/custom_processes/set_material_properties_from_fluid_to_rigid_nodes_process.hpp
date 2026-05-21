//
//   Project Name:        KratosPfemFluidApplication $
//   Created by:          $Author:           AFranci $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:        October 2018 $
//   Revision:            $Revision:             0.0 $
//
//

#if !defined(KRATOS_SET_MATERIAL_PROPERTIES_FROM_FLUID_TO_RIGID_NODES_PROCESS_H_INCLUDED)
#define KRATOS_SET_MATERIAL_PROPERTIES_FROM_FLUID_TO_RIGID_NODES_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes

#include "spatial_containers/spatial_containers.h"

#include "custom_processes/set_material_properties_from_fluid_to_rigid_nodes_process.hpp"
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
  class SetMaterialPropertiesFromFluidToRigidNodesProcess
      : public MesherProcess
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of SetMaterialPropertiesFromFluidToRigidNodesProcess
    KRATOS_CLASS_POINTER_DEFINITION(SetMaterialPropertiesFromFluidToRigidNodesProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SetMaterialPropertiesFromFluidToRigidNodesProcess(ModelPart &rRigidModelPart, ModelPart &rFluidModelPart)
        : mrRigidModelPart(rRigidModelPart),
          mrFluidModelPart(rFluidModelPart)
    {
    }

    /// Destructor.
    virtual ~SetMaterialPropertiesFromFluidToRigidNodesProcess()
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

      double density = 0;
      double bulk_modulus = 0;
      double viscosity = 0;
      double flow_index = 1;
      double yield_shear = 0;
      double adaptive_exponent = 0;
      double static_friction = 0;
      double dynamic_friction = 0;
      double inertial_number_zero = 0;
      double grain_diameter = 0;
      double grain_density = 0;
      double regularization_coefficient = 0;
      double friction_angle = 0;
      double cohesion = 0;

#pragma omp parallel
      {

        ModelPart::ElementIterator ElemBegin;
        ModelPart::ElementIterator ElemEnd;
        OpenMPUtils::PartitionedIterators(mrFluidModelPart.Elements(), ElemBegin, ElemEnd);
        for (ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
        {
          ModelPart::PropertiesType &elemProperties = itElem->GetProperties();

          density = elemProperties[DENSITY];
          bulk_modulus = elemProperties[BULK_MODULUS];
          viscosity = elemProperties[DYNAMIC_VISCOSITY];

          if (elemProperties.Has(YIELD_SHEAR)) // Bingham model
          {
            flow_index = elemProperties[FLOW_INDEX];
            yield_shear = elemProperties[YIELD_SHEAR];
            adaptive_exponent = elemProperties[ADAPTIVE_EXPONENT];
          }
          else if (elemProperties.Has(INTERNAL_FRICTION_ANGLE)) // Frictional Viscoplastic model
          {
            friction_angle = elemProperties[INTERNAL_FRICTION_ANGLE];
            cohesion = elemProperties[COHESION];
            adaptive_exponent = elemProperties[ADAPTIVE_EXPONENT];
          }
          else if (elemProperties.Has(STATIC_FRICTION)) // Mu(I)-rheology
          {
            static_friction = elemProperties[STATIC_FRICTION];
            dynamic_friction = elemProperties[DYNAMIC_FRICTION];
            inertial_number_zero = elemProperties[INERTIAL_NUMBER_ZERO];
            grain_diameter = elemProperties[GRAIN_DIAMETER];
            grain_density = elemProperties[GRAIN_DENSITY];
            regularization_coefficient = elemProperties[REGULARIZATION_COEFFICIENT];
          }
          break;
        }
      }
#pragma omp parallel
      {

        ModelPart::NodeIterator NodeBegin;
        ModelPart::NodeIterator NodeEnd;
        OpenMPUtils::PartitionedIterators(mrRigidModelPart.Nodes(), NodeBegin, NodeEnd);
        for (ModelPart::NodeIterator iNode = NodeBegin; iNode != NodeEnd; ++iNode)
        {

          if (mrFluidModelPart.GetNodalSolutionStepVariablesList().Has(BULK_MODULUS))
            iNode->FastGetSolutionStepValue(BULK_MODULUS) = bulk_modulus;

          if (mrFluidModelPart.GetNodalSolutionStepVariablesList().Has(DENSITY))
            iNode->FastGetSolutionStepValue(DENSITY) = density;

          if (mrFluidModelPart.GetNodalSolutionStepVariablesList().Has(DYNAMIC_VISCOSITY))
            iNode->FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = viscosity;

          if (mrFluidModelPart.GetNodalSolutionStepVariablesList().Has(YIELD_SHEAR))
            iNode->FastGetSolutionStepValue(YIELD_SHEAR) = yield_shear;

          if (mrFluidModelPart.GetNodalSolutionStepVariablesList().Has(FLOW_INDEX))
            iNode->FastGetSolutionStepValue(FLOW_INDEX) = flow_index;

          if (mrFluidModelPart.GetNodalSolutionStepVariablesList().Has(ADAPTIVE_EXPONENT))
            iNode->FastGetSolutionStepValue(ADAPTIVE_EXPONENT) = adaptive_exponent;

          if (mrFluidModelPart.GetNodalSolutionStepVariablesList().Has(INTERNAL_FRICTION_ANGLE))
            iNode->FastGetSolutionStepValue(INTERNAL_FRICTION_ANGLE) = friction_angle;

          if (mrFluidModelPart.GetNodalSolutionStepVariablesList().Has(COHESION))
            iNode->FastGetSolutionStepValue(COHESION) = cohesion;

          if (mrFluidModelPart.GetNodalSolutionStepVariablesList().Has(ADAPTIVE_EXPONENT))
            iNode->FastGetSolutionStepValue(ADAPTIVE_EXPONENT) = adaptive_exponent;

          if (mrFluidModelPart.GetNodalSolutionStepVariablesList().Has(STATIC_FRICTION))
            iNode->FastGetSolutionStepValue(STATIC_FRICTION) = static_friction;

          if (mrFluidModelPart.GetNodalSolutionStepVariablesList().Has(DYNAMIC_FRICTION))
            iNode->FastGetSolutionStepValue(DYNAMIC_FRICTION) = dynamic_friction;

          if (mrFluidModelPart.GetNodalSolutionStepVariablesList().Has(INERTIAL_NUMBER_ZERO))
            iNode->FastGetSolutionStepValue(INERTIAL_NUMBER_ZERO) = inertial_number_zero;

          if (mrFluidModelPart.GetNodalSolutionStepVariablesList().Has(GRAIN_DIAMETER))
            iNode->FastGetSolutionStepValue(GRAIN_DIAMETER) = grain_diameter;

          if (mrFluidModelPart.GetNodalSolutionStepVariablesList().Has(GRAIN_DENSITY))
            iNode->FastGetSolutionStepValue(GRAIN_DENSITY) = grain_density;

          if (mrFluidModelPart.GetNodalSolutionStepVariablesList().Has(REGULARIZATION_COEFFICIENT))
            iNode->FastGetSolutionStepValue(REGULARIZATION_COEFFICIENT) = regularization_coefficient;
        }
      }

      KRATOS_CATCH(" ")
    };

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
      return "SetMaterialPropertiesFromFluidToRigidNodesProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
      rOStream << "SetMaterialPropertiesFromFluidToRigidNodesProcess";
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
    ModelPart &mrRigidModelPart;
    ModelPart &mrFluidModelPart;

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
    SetMaterialPropertiesFromFluidToRigidNodesProcess &operator=(SetMaterialPropertiesFromFluidToRigidNodesProcess const &rOther);

    /// Copy constructor.
    // SetMaterialPropertiesFromFluidToRigidNodesProcess(SetMaterialPropertiesFromFluidToRigidNodesProcess const& rOther);

    ///@}

  }; // namespace Kratos

  ///@}

  ///@name Type Definitions
  ///@{

  ///@}
  ///@name Input and output
  ///@{

  /// input stream function
  inline std::istream &operator>>(std::istream &rIStream,
                                  SetMaterialPropertiesFromFluidToRigidNodesProcess &rThis);

  /// output stream function
  inline std::ostream &operator<<(std::ostream &rOStream,
                                  const SetMaterialPropertiesFromFluidToRigidNodesProcess &rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}

} // namespace Kratos.

#endif // KRATOS_SET_MATERIAL_PROPERTIES_FROM_FLUID_TO_RIGID_NODES_PROCESS_H_INCLUDED  defined
