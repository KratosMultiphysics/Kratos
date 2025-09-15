//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:               June 2021 $
//   Revision:            $Revision:                 0.0 $
//
//

#if !defined(KRATOS_THREE_STEP_SECOND_ORDER_PSPG_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED)
#define KRATOS_THREE_STEP_SECOND_ORDER_PSPG_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

#include "pfem_fluid_dynamics_application_variables.h"

#include "custom_elements/three_step_second_order_updated_lagrangian_element.h"

#include "includes/model_part.h"

namespace Kratos
{

  ///@addtogroup FluidDynamicsApplication
  ///@{

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

  /// A stabilized element for the incompressible Navier-Stokes equations.
  /**
   */

  template <unsigned int TDim>
  class ThreeStepSecondOrderPspgUpdatedLagrangianElement : public ThreeStepSecondOrderUpdatedLagrangianElement<TDim>
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ThreeStepSecondOrderUpdatedLagrangianElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ThreeStepSecondOrderPspgUpdatedLagrangianElement);

    typedef ThreeStepSecondOrderUpdatedLagrangianElement<TDim> BaseType;

    /// Node type (default is: Node)
    typedef Node NodeType;

    /// Geometry type (using with given NodeType)
    typedef Geometry<NodeType> GeometryType;

    /// Definition of nodes container type, redefined from GeometryType
    typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

    /// Vector type for local contributions to the linear system
    typedef Vector VectorType;

    /// Matrix type for local contributions to the linear system
    typedef Matrix MatrixType;

    typedef std::size_t IndexType;

    typedef std::size_t SizeType;

    typedef std::vector<std::size_t> EquationIdVectorType;

    typedef std::vector<Dof<double>::Pointer> DofsVectorType;

    typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

    /// Type for shape function values container
    typedef Kratos::Vector ShapeFunctionsType;

    /// Type for a matrix containing the shape function gradients
    typedef Kratos::Matrix ShapeFunctionDerivativesType;

    /// Type for an array of shape function gradient matrices
    typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

    /// Reference type definition for constitutive laws
    typedef ConstitutiveLaw ConstitutiveLawType;

    /// Pointer type for constitutive laws
    typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;

    typedef typename BaseType::PropertiesType PropertiesType;

    typedef typename BaseType::PropertiesType::Pointer pPropertiesType;

    typedef typename BaseType::ElementalVariables ElementalVariables;

    ///@}
    ///@name Life Cycle
    ///@{

    //Constructors.

    /// Default constructor.
    /**
       * @param NewId Index number of the new element (optional)
       */
    ThreeStepSecondOrderPspgUpdatedLagrangianElement(IndexType NewId = 0) : BaseType(NewId)
    {
    }

    /// Constructor using an array of nodes.
    /**
       * @param NewId Index of the new element
       * @param ThisNodes An array containing the nodes of the new element
       */
    ThreeStepSecondOrderPspgUpdatedLagrangianElement(IndexType NewId, const NodesArrayType &ThisNodes) : BaseType(NewId, ThisNodes)
    {
    }

    /// Constructor using a geometry object.
    /**
       * @param NewId Index of the new element
       * @param pGeometry Pointer to a geometry object
       */
    ThreeStepSecondOrderPspgUpdatedLagrangianElement(IndexType NewId, GeometryType::Pointer pGeometry) : BaseType(NewId, pGeometry)
    {
    }

    /// Constructor using geometry and properties.
    /**
       * @param NewId Index of the new element
       * @param pGeometry Pointer to a geometry object
       * @param pProperties Pointer to the element's properties
       */
    ThreeStepSecondOrderPspgUpdatedLagrangianElement(IndexType NewId, GeometryType::Pointer pGeometry, pPropertiesType pProperties) : BaseType(NewId, pGeometry, pProperties)
    {
    }

    /// copy constructor

    ThreeStepSecondOrderPspgUpdatedLagrangianElement(ThreeStepSecondOrderPspgUpdatedLagrangianElement const &rOther) : BaseType(rOther){};

    /// Destructor.
    virtual ~ThreeStepSecondOrderPspgUpdatedLagrangianElement()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /// Create a new element of this type
    /**
       * Returns a pointer to a new ThreeStepSecondOrderPspgUpdatedLagrangianElement element, created using given input
       * @param NewId: the ID of the new element
       * @param ThisNodes: the nodes of the new element
       * @param pProperties: the properties assigned to the new element
       * @return a Pointer to the new element
       */
    Element::Pointer Create(IndexType NewId, NodesArrayType const &ThisNodes,
                            pPropertiesType pProperties) const override
    {
      return Element::Pointer(new ThreeStepSecondOrderPspgUpdatedLagrangianElement(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
    }

    Element::Pointer Clone(IndexType NewId, NodesArrayType const &ThisNodes) const override;

    /// Calculate the element's local contribution to the system for the current step.
    void CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
                              VectorType &rRightHandSideVector,
                              const ProcessInfo &rCurrentProcessInfo) override;

    void CalculateFSplusPSPGPressureSystem(MatrixType &rLeftHandSideMatrix,
                                           VectorType &rRightHandSideVector,
                                           const ProcessInfo &rCurrentProcessInfo);
    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Elemental Data
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
      std::stringstream buffer;
      buffer << "ThreeStepSecondOrderPspgUpdatedLagrangianElement #" << this->Id();
      return buffer.str();
    }

    //        /// Print object's data.
    //        virtual void PrintData(std::ostream& rOStream) const;

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

    ConstitutiveLaw::Pointer mpConstitutiveLaw = nullptr;

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

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer &rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
    }

    void load(Serializer &rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
    }

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
    ThreeStepSecondOrderPspgUpdatedLagrangianElement &operator=(ThreeStepSecondOrderPspgUpdatedLagrangianElement const &rOther);

    /* /// Copy constructor. */
    /* ThreeStepSecondOrderPspgUpdatedLagrangianElement(ThreeStepSecondOrderPspgUpdatedLagrangianElement const& rOther); */

    ///@}

  }; // Class ThreeStepSecondOrderPspgUpdatedLagrangianElement

  ///@}

  ///@name Type Definitions
  ///@{

  ///@}
  ///@name Input and output
  ///@{

  /// input stream function
  template <unsigned int TDim>
  inline std::istream &operator>>(std::istream &rIStream,
                                  ThreeStepSecondOrderPspgUpdatedLagrangianElement<TDim> &rThis)
  {
    return rIStream;
  }

  /// output stream function
  template <unsigned int TDim>
  inline std::ostream &operator<<(std::ostream &rOStream,
                                  const ThreeStepSecondOrderPspgUpdatedLagrangianElement<TDim> &rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }

} // namespace Kratos.

#endif // KRATOS_THREE_STEP_SECOND_ORDER_UPDATED_LAGRANGIAN_ELEMENT  defined
