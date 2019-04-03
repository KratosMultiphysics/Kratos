//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 0.0 $
//
//

#if !defined(KRATOS_UPDATED_LAGRANGIAN_V_EXPLICIT_SOLID_ELEMENT_H_INCLUDED )
#define  KRATOS_UPDATED_LAGRANGIAN_V_EXPLICIT_SOLID_ELEMENT_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

#include "custom_elements/two_step_updated_lagrangian_V_P_explicit_solid_element.h" 

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
  /* class UpdatedLagrangianVExplicitSolidElement : public Element */
  template< unsigned int TDim > 
    /* class UpdatedLagrangianVExplicitSolidElement : public TwoStepUpdatedLagrangianVPExplicitElement<TDim> */
    class UpdatedLagrangianVExplicitSolidElement : public TwoStepUpdatedLagrangianVPExplicitSolidElement<TDim>
    {
  
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of UpdatedLagrangianVExplicitSolidElement
      KRATOS_CLASS_POINTER_DEFINITION(UpdatedLagrangianVExplicitSolidElement);

      ///base type: 
      typedef TwoStepUpdatedLagrangianVPExplicitSolidElement<TDim> BaseType;

      /// Node type (default is: Node<3>)
      typedef Node <3> NodeType;

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

      typedef std::vector< Dof<double>::Pointer > DofsVectorType;

      typedef PointerVectorSet<Dof<double>, IndexedObject> DofsArrayType;

      /// Type for shape function values container
      typedef Kratos::Vector ShapeFunctionsType;

      /// Type for a matrix containing the shape function gradients
      typedef Kratos::Matrix ShapeFunctionDerivativesType;

      /// Type for an array of shape function gradient matrices
      typedef GeometryType::ShapeFunctionsGradientsType ShapeFunctionDerivativesArrayType;

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
    UpdatedLagrangianVExplicitSolidElement(IndexType NewId = 0) :
      BaseType(NewId)
      {}

      /// Constructor using an array of nodes.
      /**
       * @param NewId Index of the new element
       * @param ThisNodes An array containing the nodes of the new element
       */
    UpdatedLagrangianVExplicitSolidElement(IndexType NewId, const NodesArrayType& ThisNodes) :
      BaseType(NewId, ThisNodes)
	{}

      /// Constructor using a geometry object.
      /**
       * @param NewId Index of the new element
       * @param pGeometry Pointer to a geometry object
       */
    UpdatedLagrangianVExplicitSolidElement(IndexType NewId, GeometryType::Pointer pGeometry) :
      BaseType(NewId, pGeometry)
	{}

      /// Constuctor using geometry and properties.
      /**
       * @param NewId Index of the new element
       * @param pGeometry Pointer to a geometry object
       * @param pProperties Pointer to the element's properties
       */
    UpdatedLagrangianVExplicitSolidElement(IndexType NewId, GeometryType::Pointer pGeometry, pPropertiesType pProperties) : BaseType(NewId, pGeometry, pProperties)
	{}


      /// copy constructor

    UpdatedLagrangianVExplicitSolidElement(UpdatedLagrangianVExplicitSolidElement const& rOther):
      BaseType(rOther)
      {}

      /// Destructor.
      virtual ~UpdatedLagrangianVExplicitSolidElement()
	{}


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      /// Create a new element of this type
      /**
       * Returns a pointer to a new UpdatedLagrangianVExplicitSolidElement element, created using given input
       * @param NewId: the ID of the new element
       * @param ThisNodes: the nodes of the new element
       * @param pProperties: the properties assigned to the new element
       * @return a Pointer to the new element
       */
      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
			      pPropertiesType pProperties) const override
	{
	  return Element::Pointer(new UpdatedLagrangianVExplicitSolidElement(NewId, BaseType::GetGeometry().Create(ThisNodes), pProperties));
	}

      Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;


      /* virtual void Initialize(); */

      /// Initializes the element and all geometric information required for the problem.
      /* virtual void InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo); */


      /* virtual void InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo); */

      void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
				 ProcessInfo& rCurrentProcessInfo) override
      {
	KRATOS_TRY;
	KRATOS_THROW_ERROR(std::logic_error,"UpdatedLagrangianVExplicitSolidElement::CalculateLeftHandSide not implemented","");
	KRATOS_CATCH("");
      }


      // The following methods have different implementations depending on TDim
      /// Provides the global indices for each one of this element's local rows
      /**
       * this determines the elemental equation ID vector for all elemental
       * DOFs
       * @param rResult A vector containing the global Id of each row
       * @param rCurrentProcessInfo the current process info object (unused)
       */

      /// Returns a list of the element's Dofs
      /**
       * @param ElementalDofList the list of DOFs
       * @param rCurrentProcessInfo the current process info instance
       */


      ///@}
      ///@name Access
      ///@{

      ///@}
      ///@name Elemental Data
      ///@{

      /// Checks the input and that all required Kratos variables have been registered.
      /**
       * This function provides the place to perform checks on the completeness of the input.
       * It is designed to be called only once (or anyway, not often) typically at the beginning
       * of the calculations, so to verify that nothing is missing from the input
       * or that no common error is found.
       * @param rCurrentProcessInfo The ProcessInfo of the ModelPart that contains this element.
       * @return 0 if no errors were found.
       */
      /* virtual int Check(const ProcessInfo& rCurrentProcessInfo); */

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
	  buffer << "UpdatedLagrangianVExplicitSolidElement #" << BaseType::Id();
	  return buffer.str();
	}

      /// Print information about this object.
      void PrintInfo(std::ostream& rOStream) const override
      {
	rOStream << "UpdatedLagrangianVExplicitSolidElement" << TDim << "D";
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


      ///@}
      ///@name Protected Operators
      ///@{


      ///@}
      ///@name Protected Operations
      ///@{

      /// Add integration point contribution to the mass matrix.
      /**
       * A constistent mass matrix is used.
       * @param rMassMatrix The local matrix where the result will be added.
       * @param rN Elemental shape functions.
       * @param Weight Multiplication coefficient for the matrix, typically Density times integration point weight.
       */
   
      
      /* bool CalcMechanicsUpdated(ElementalVariables & rElementalVariables, */
      /* 				const ProcessInfo& rCurrentProcessInfo, */
      /* 				const ShapeFunctionDerivativesType& rDN_DX, */
      /* 				unsigned int g); */

	
      void CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables,
					    double TimeStep,
					    unsigned int g) override;


      double GetThetaMomentum () override{return 1.0;};

      double GetThetaContinuity () override{return 1.0;};


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
      UpdatedLagrangianVExplicitSolidElement & operator=(UpdatedLagrangianVExplicitSolidElement const& rOther);

      /* /// Copy constructor. */
      /* UpdatedLagrangianVExplicitSolidElement(UpdatedLagrangianVExplicitSolidElement const& rOther); */

      ///@}

    }; // Class UpdatedLagrangianVExplicitSolidElement

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  template< unsigned int TDim >
    inline std::istream& operator >>(std::istream& rIStream,
                                     UpdatedLagrangianVExplicitSolidElement<TDim>& rThis)
    {
      return rIStream;
    }

  /// output stream function
  template< unsigned int TDim >
    inline std::ostream& operator <<(std::ostream& rOStream,
                                     const UpdatedLagrangianVExplicitSolidElement<TDim>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }

} // namespace Kratos.

#endif // KRATOS_TWO_STEP_UPDATED_LAGRANGIAN_V_EXPLICIT_SOLID_ELEMENT  defined
