//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 0.0 $
//
//

#if !defined(KRATOS_TWO_STEP_UPDATED_LAGRANGIAN_V_P_IMPLICIT_NODALLY_INTEGRATED_ELEMENT_H_INCLUDED )
#define  KRATOS_TWO_STEP_UPDATED_LAGRANGIAN_V_P_IMPLICIT_NODALLY_INTEGRATED_ELEMENT_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes


// Project includes 
#include "containers/array_1d.h"
#include "includes/define.h"
/* #include "includes/element.h" */
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

#include "pfem_fluid_dynamics_application_variables.h"

#include "custom_elements/two_step_updated_lagrangian_element.h"

#include "includes/model_part.h"
/* #include "includes/node.h" */

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

  template< unsigned int TDim >
    class TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement : public TwoStepUpdatedLagrangianElement<TDim>
    {
  
    protected:


      ///@name Protected static Member Variables
      ///@{


      ///@}
      ///@name Protected member Variables
      ///@{
 

    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement
      KRATOS_CLASS_POINTER_DEFINITION(TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement);

      ///base type:
      typedef TwoStepUpdatedLagrangianElement<TDim> BaseType;
	    
      /// Node type (default is: Node<3>)
      typedef Node <3> NodeType;

      /// Geometry type (using with given NodeType)
      typedef Geometry<NodeType> GeometryType;

      /// Definition of nodes container type, redefined from GeometryType
      typedef Geometry<NodeType>::PointsArrayType NodesArrayType;

      /// Vector type for local contributions to the linear system
      typedef Kratos::Vector VectorType;

      /// Matrix type for local contributions to the linear system
      typedef Kratos::Matrix MatrixType;

      typedef std::size_t IndexType;

      typedef std::size_t SizeType;

      typedef std::vector<std::size_t> EquationIdVectorType;

      typedef std::vector< Dof<double>::Pointer > DofsVectorType;

      typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

      typedef VectorMap<IndexType, DataValueContainer> SolutionStepsElementalDataContainerType;

      /// Type for shape function values container
      typedef Kratos::Vector ShapeFunctionsType;

      /// Type for a matrix containing the shape function gradients
      typedef Kratos::Matrix ShapeFunctionDerivativesType;

      /// Type for an array of shape function gradient matrices
      typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

      /* typedef Element::PropertiesType::Pointer PropertiesType::Pointer; */

      typedef typename BaseType::PropertiesType PropertiesType;

      typedef typename BaseType::PropertiesType::Pointer pPropertiesType;

      typedef typename BaseType::ElementalVariables ElementalVariables;


      ///@}
      ///@name Life Cycle
      ///@{

      //Constructors.

      /// Default constuctor.
      /**
       * @param NewId Index number of the new element (optional)
       */
    TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement(IndexType NewId = 0) :
      BaseType(NewId)
      {}

      /// Constructor using an array of nodes.
      /**
       * @param NewId Index of the new element
       * @param ThisNodes An array containing the nodes of the new element
       */
    TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement(IndexType NewId, const NodesArrayType& ThisNodes) :
      BaseType(NewId, ThisNodes)
        {}

      /// Constructor using a geometry object.
      /**
       * @param NewId Index of the new element
       * @param pGeometry Pointer to a geometry object
       */
    TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement(IndexType NewId, GeometryType::Pointer pGeometry) :
      BaseType(NewId, pGeometry)
        {}

      /// Constuctor using geometry and properties.
      /**
       * @param NewId Index of the new element
       * @param pGeometry Pointer to a geometry object
       * @param pProperties Pointer to the element's properties
       */
    TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement(IndexType NewId, GeometryType::Pointer pGeometry, pPropertiesType pProperties) :
      BaseType(NewId, pGeometry, pProperties)
        {}



      /// copy constructor

    TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement(TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement const& rOther):
      BaseType(rOther)
      {}
      
      /// Destructor.
      virtual ~TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement()
        {}


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      /// Create a new element of this type
      /**
       * Returns a pointer to a new TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement element, created using given input
       * @param NewId: the ID of the new element
       * @param ThisNodes: the nodes of the new element
       * @param pProperties: the properties assigned to the new element
       * @return a Pointer to the new element
       */
      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
			      pPropertiesType pProperties) const override
        {
	  return Element::Pointer(new TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement(NewId, this->GetGeometry().Create(ThisNodes), pProperties));
        }


      Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

      void Initialize() override{};

      /// Initializes the element and all geometric information required for the problem.
      void InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo) override{};

      void InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo) override;

      void CalcElementalStrains(ElementalVariables & rElementalVariables,
				const ProcessInfo &rCurrentProcessInfo,
				const ShapeFunctionDerivativesType& rDN_DX);

      void GetNodesPosition(Vector& rValues,
			    const ProcessInfo& rCurrentProcessInfo,
			    double theta);
            
      void CalculateGeometryData(ShapeFunctionDerivativesArrayType& rDN_DX,
				 Matrix& rNContainer,
				 Vector& rGaussWeights);
      
      void InitializeElementalVariables(ElementalVariables & rElementalVariables) override{
	KRATOS_TRY;
	KRATOS_THROW_ERROR(std::logic_error,"InitializeElementalVariables","");
	KRATOS_CATCH("");
      };



      /// Checks the input and that all required Kratos variables have been registered.
      /**
       * This function provides the place to perform checks on the completeness of the input.
       * It is designed to be called only once (or anyway, not often) typically at the beginning
       * of the calculations, so to verify that nothing is missing from the input
       * or that no common error is found.
       * @param rCurrentProcessInfo The ProcessInfo of the ModelPart that contains this element.
       * @return 0 if no errors were found.
       */
      int Check(const ProcessInfo& rCurrentProcessInfo) override;

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
	  buffer << "TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement #" <<  BaseType::Id();
	  return buffer.str();
	}

      /// Print information about this object.
      void PrintInfo(std::ostream& rOStream) const override
      {
	rOStream << "TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement" << TDim << "D";
      }

    protected:

      void NodalFreeSurfaceLength(unsigned int nodeIndex);

      void GetValueOnIntegrationPoints( const Variable<double>& rVariable,
					std::vector<double>& rValues,
					const ProcessInfo& rCurrentProcessInfo ) override;
   
    private:


      friend class Serializer;

      void save(Serializer& rSerializer) const override
      {
	KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
      }

      void load(Serializer& rSerializer) override
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
      TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement & operator=(TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement const& rOther);

      /* /// Copy constructor. */
      /* TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement(TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement const& rOther); */

      ///@}

    }; // Class TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  template< unsigned int TDim >
    inline std::istream& operator >>(std::istream& rIStream,
                                     TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>& rThis)
    {
      return rIStream;
    }

  /// output stream function
  template< unsigned int TDim >
    inline std::ostream& operator <<(std::ostream& rOStream,
                                     const TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<TDim>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }

} // namespace Kratos.

#endif // KRATOS_TWO_STEP_UPDATED_LAGRANGIAN_V_P_IMPLICIT_NODALLY_INTEGRATED_ELEMENT  defined
