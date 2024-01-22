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

    void Execute() override{
        KRATOS_TRY

#pragma omp parallel
        {

            ModelPart::ElementIterator ElemBegin;
    ModelPart::ElementIterator ElemEnd;
    OpenMPUtils::PartitionedIterators(mrModelPart.Elements(), ElemBegin, ElemEnd);
    for (ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem)
    {
      ModelPart::PropertiesType &elemProperties = itElem->GetProperties();

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

      double density = elemProperties[DENSITY];
      double bulk_modulus = elemProperties[BULK_MODULUS];
      double viscosity = elemProperties[DYNAMIC_VISCOSITY];
      unsigned int elem_property_id = elemProperties.Id();

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

      Geometry<Node> &rGeom = itElem->GetGeometry();
      const SizeType NumNodes = rGeom.PointsNumber();
      for (SizeType i = 0; i < NumNodes; ++i)
      {

        if (mrModelPart.GetNodalSolutionStepVariablesList().Has(PROPERTY_ID))
        {
          rGeom[i].FastGetSolutionStepValue(PROPERTY_ID) = elem_property_id;
        }

        if (mrModelPart.GetNodalSolutionStepVariablesList().Has(BULK_MODULUS))
          rGeom[i].FastGetSolutionStepValue(BULK_MODULUS) = bulk_modulus;

        if (mrModelPart.GetNodalSolutionStepVariablesList().Has(DENSITY))
          rGeom[i].FastGetSolutionStepValue(DENSITY) = density;

        if (mrModelPart.GetNodalSolutionStepVariablesList().Has(DYNAMIC_VISCOSITY))
          rGeom[i].FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = viscosity;

        if (mrModelPart.GetNodalSolutionStepVariablesList().Has(YIELD_SHEAR))
          rGeom[i].FastGetSolutionStepValue(YIELD_SHEAR) = yield_shear;

        if (mrModelPart.GetNodalSolutionStepVariablesList().Has(FLOW_INDEX))
          rGeom[i].FastGetSolutionStepValue(FLOW_INDEX) = flow_index;

        if (mrModelPart.GetNodalSolutionStepVariablesList().Has(ADAPTIVE_EXPONENT))
          rGeom[i].FastGetSolutionStepValue(ADAPTIVE_EXPONENT) = adaptive_exponent;

        if (mrModelPart.GetNodalSolutionStepVariablesList().Has(INTERNAL_FRICTION_ANGLE))
          rGeom[i].FastGetSolutionStepValue(INTERNAL_FRICTION_ANGLE) = friction_angle;

        if (mrModelPart.GetNodalSolutionStepVariablesList().Has(COHESION))
          rGeom[i].FastGetSolutionStepValue(COHESION) = cohesion;

        if (mrModelPart.GetNodalSolutionStepVariablesList().Has(ADAPTIVE_EXPONENT))
          rGeom[i].FastGetSolutionStepValue(ADAPTIVE_EXPONENT) = adaptive_exponent;

        if (mrModelPart.GetNodalSolutionStepVariablesList().Has(STATIC_FRICTION))
          rGeom[i].FastGetSolutionStepValue(STATIC_FRICTION) = static_friction;

        if (mrModelPart.GetNodalSolutionStepVariablesList().Has(DYNAMIC_FRICTION))
          rGeom[i].FastGetSolutionStepValue(DYNAMIC_FRICTION) = dynamic_friction;

        if (mrModelPart.GetNodalSolutionStepVariablesList().Has(INERTIAL_NUMBER_ZERO))
          rGeom[i].FastGetSolutionStepValue(INERTIAL_NUMBER_ZERO) = inertial_number_zero;

        if (mrModelPart.GetNodalSolutionStepVariablesList().Has(GRAIN_DIAMETER))
          rGeom[i].FastGetSolutionStepValue(GRAIN_DIAMETER) = grain_diameter;

        if (mrModelPart.GetNodalSolutionStepVariablesList().Has(GRAIN_DENSITY))
          rGeom[i].FastGetSolutionStepValue(GRAIN_DENSITY) = grain_density;

        if (mrModelPart.GetNodalSolutionStepVariablesList().Has(REGULARIZATION_COEFFICIENT))
          rGeom[i].FastGetSolutionStepValue(REGULARIZATION_COEFFICIENT) = regularization_coefficient;
      }
    }

  }

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
