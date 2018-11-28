//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:              April 2018 $
//   Revision:            $Revision:                 0.0 $
//
//

#if !defined(KRATOS_TWO_STEP_UPDATED_LAGRANGIAN_V_P_EXPLICIT_ELEMENT_H_INCLUDED )
#define  KRATOS_TWO_STEP_UPDATED_LAGRANGIAN_V_P_EXPLICIT_ELEMENT_H_INCLUDED

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
    class TwoStepUpdatedLagrangianVPExplicitElement : public TwoStepUpdatedLagrangianElement<TDim>
    /* class TwoStepUpdatedLagrangianVPExplicitElement : public Element */
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

      /// Pointer definition of TwoStepUpdatedLagrangianVPExplicitElement
      KRATOS_CLASS_POINTER_DEFINITION(TwoStepUpdatedLagrangianVPExplicitElement);


      ///base type: 
      typedef TwoStepUpdatedLagrangianElement<TDim> BaseType;
      /* typedef Element BaseType; */

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
    TwoStepUpdatedLagrangianVPExplicitElement(IndexType NewId = 0) :
      BaseType(NewId)
      {}

      /// Constructor using an array of nodes.
      /**
       * @param NewId Index of the new element
       * @param ThisNodes An array containing the nodes of the new element
       */
    TwoStepUpdatedLagrangianVPExplicitElement(IndexType NewId, const NodesArrayType& ThisNodes) :
      BaseType(NewId, ThisNodes)
        {}

      /// Constructor using a geometry object.
      /**
       * @param NewId Index of the new element
       * @param pGeometry Pointer to a geometry object
       */
    TwoStepUpdatedLagrangianVPExplicitElement(IndexType NewId, GeometryType::Pointer pGeometry) :
      BaseType(NewId, pGeometry)
        {}

      /// Constuctor using geometry and properties.
      /**
       * @param NewId Index of the new element
       * @param pGeometry Pointer to a geometry object
       * @param pProperties Pointer to the element's properties
       */
    TwoStepUpdatedLagrangianVPExplicitElement(IndexType NewId, GeometryType::Pointer pGeometry, pPropertiesType pProperties) :
      BaseType(NewId, pGeometry, pProperties)
        {}



      /// copy constructor

    TwoStepUpdatedLagrangianVPExplicitElement(TwoStepUpdatedLagrangianVPExplicitElement const& rOther):
      BaseType(rOther)
      {}
      
      /// Destructor.
      virtual ~TwoStepUpdatedLagrangianVPExplicitElement()
        {}


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      /// Create a new element of this type
      /**
       * Returns a pointer to a new TwoStepUpdatedLagrangianVPExplicitElement element, created using given input
       * @param NewId: the ID of the new element
       * @param ThisNodes: the nodes of the new element
       * @param pProperties: the properties assigned to the new element
       * @return a Pointer to the new element
       */
      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
			      pPropertiesType pProperties) const override
        {
	  return Element::Pointer(new TwoStepUpdatedLagrangianVPExplicitElement(NewId, BaseType::GetGeometry().Create(ThisNodes), pProperties));
        }


      Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

      void Initialize() override;

      /// Initializes the element and all geometric information required for the problem.
      void InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo) override{};


      void InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo) override{};

      /// Calculate the element's local contribution to the system for the current step.
      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
				VectorType& rRightHandSideVector,
				ProcessInfo& rCurrentProcessInfo) override;


      void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
				 ProcessInfo& rCurrentProcessInfo) override {};

      void CalculateRightHandSide(VectorType& rRightHandSideVector,
				  ProcessInfo& rCurrentProcessInfo) override;

      
      void CalculateRightHandSideMomentum(VectorType& rRightHandSideVector,
					  ProcessInfo& rCurrentProcessInfo);
 
      /* // The following methods have different implementations depending on TDim */
      /* /// Provides the global indices for each one of this element's local rows */
      /* /\** */
      /*  * this determines the elemental equation ID vector for all elemental */
      /*  * DOFs */
      /*  * @param rResult A vector containing the global Id of each row */
      /*  * @param rCurrentProcessInfo the current process info object (unused) */
      /*  *\/ */
      /* virtual void EquationIdVector(EquationIdVectorType& rResult, */
      /* 				    ProcessInfo& rCurrentProcessInfo); */

      /* /// Returns a list of the element's Dofs */
      /* /\** */
      /*  * @param ElementalDofList the list of DOFs */
      /*  * @param rCurrentProcessInfo the current process info instance */
      /*  *\/ */
      /* virtual void GetDofList(DofsVectorType& rElementalDofList, */
      /* 			      ProcessInfo& rCurrentProcessInfo); */

    
      /* virtual GeometryData::IntegrationMethod GetIntegrationMethod() const; */
    
      void UpdateCauchyStress(unsigned int g,ProcessInfo& rCurrentProcessInfo) override{};

      /* virtual void InitializeElementalVariables(ElementalVariables & rElementalVariables){ */
      /* 	KRATOS_TRY; */
      /* 	KRATOS_THROW_ERROR(std::logic_error,"InitializeElementalVariables",""); */
      /* 	KRATOS_CATCH(""); */
      /* }; */


      void AddExplicitContribution(const VectorType& rRHSVector, 
				   const Variable<VectorType>& rRHSVariable, 
				   Variable<array_1d<double,3> >& rDestinationVariable, 
				   const ProcessInfo& rCurrentProcessInfo) override;

      void AddExplicitContribution(const VectorType& rRHSVector, 
				   const Variable<VectorType>& rRHSVariable, 
				   Variable<double > & rDestinationVariable, 
				   const ProcessInfo& rCurrentProcessInfo) override;
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
	  buffer << "TwoStepUpdatedLagrangianVPExplicitElement #" << BaseType::Id();
	  return buffer.str();
	}

      /// Print information about this object.
      void PrintInfo(std::ostream& rOStream) const override
      {
	rOStream << "TwoStepUpdatedLagrangianVPExplicitElement" << TDim << "D";
      }

      //        /// Print object's data.
      //        virtual void PrintData(std::ostream& rOStream) const;

      ///@}
      ///@name Friends
      ///@{

      ///@}
    protected:

      /* double mMaterialDeviatoricCoefficient=0; */
      /* double mMaterialVolumetricCoefficient=0; */
      /* double mMaterialDensity=0; */

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

      void GetValueOnIntegrationPoints( const Variable<double>& rVariable,
					std::vector<double>& rValues,
					const ProcessInfo& rCurrentProcessInfo ) override;
      
      void CalculateLocalMomentumEquations(MatrixType& rLeftHandSideMatrix,
					   VectorType& rRightHandSideVector,
					   ProcessInfo& rCurrentProcessInfo) override{std::cout<<"CalculateLocalMomentumEquations is not defined for explicit element"<<std::endl;};

      void CalculateLocalContinuityEqForPressure(MatrixType& rLeftHandSideMatrix,
						 VectorType& rRightHandSideVector,
						 ProcessInfo& rCurrentProcessInfo) override{
	std::cout<<"CalculateLocalContinuityEqForPressure not defined for explicit element"<<std::endl;};
      
      void CalculateExplicitContinuityEquation(MatrixType& rLeftHandSideMatrix,
					       VectorType& rRightHandSideVector,
					       ProcessInfo& rCurrentProcessInfo) override;


      void ComputeMaterialParameters (double& Density,
				      double& DeviatoricCoeff,
				      double& VolumetricCoeff,
				      ProcessInfo& rCurrentProcessInfo,
				      ElementalVariables& rElementalVariables) override{};

      void ComputeMaterialParametersGranularGas (double& Density,
						 double& DeviatoricCoeff,
						 double& VolumetricCoeff,
						 ProcessInfo& rCurrentProcessInfo,
						 ElementalVariables& rElementalVariables) override{};

      double GetThetaMomentum () override{return 1.0;};

      double GetThetaContinuity () override{return 1.0;};
     
      void GetPositions(Vector& rValues,
			const ProcessInfo& rCurrentProcessInfo,
			const double theta) override{};

      /// Add integration point contribution to the mass matrix.
      /**
       * A constistent mass matrix is used.
       * @param rMassMatrix The local matrix where the result will be added.
       * @param rN Elemental shape functions.
       * @param Weight Multiplication coefficient for the matrix, typically Density times integration point weight.
       */
      void CalculateMassMatrix(Matrix& rMassMatrix,
			       ProcessInfo& rCurrentProcessInfo) override; 
      
      void CalculateMassMatrixMomentum(Matrix& rMassMatrix,
				       ProcessInfo& rCurrentProcessInfo);
         
      void ComputeBulkMatrixRHS(MatrixType& BulkMatrix,
				const double Weight) override{};    
      
      /* bool CalcMechanicsUpdated(ElementalVariables & rElementalVariables, */
      /* 				const ProcessInfo& rCurrentProcessInfo, */
      /* 				const ShapeFunctionDerivativesType& rDN_DX, */
      /* 				unsigned int g); */

	
      void CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables,
					    double TimeStep,
					    unsigned int g)override{};


      /* /// Helper function to print results on gauss points */
      /* /\** Reads a variable from the element's database and returns it in a format */
      /*  * that can be used by GetValueOnIntegrationPoints functions. */
      /*  * @see GetValueOnIntegrationPoints */
      /*  *\/ */
      /* template<class TValueType> */
      /* 	void GetElementalValueForOutput(const Kratos::Variable<TValueType>& rVariable, */
      /* 					std::vector<TValueType>& rOutput) */
      /* 	{ */
      /* 	  unsigned int NumValues = this->GetGeometry().IntegrationPointsNumber(GeometryData::GI_GAUSS_1); */
      /* 	  /\* unsigned int NumValues = this->GetGeometry().IntegrationPointsNumber(GeometryData::GI_GAUSS_4); *\/ */
      /* 	  rOutput.resize(NumValues); */
      /* 	  /\* */
      /* 	    The cast is done to avoid modification of the element's data. Data modification */
      /* 	    would happen if rVariable is not stored now (would initialize a pointer to &rVariable */
      /* 	    with associated value of 0.0). This is catastrophic if the variable referenced */
      /* 	    goes out of scope. */
      /* 	  *\/ */
      /* 	  const TwoStepUpdatedLagrangianVPExplicitElement<TDim>* const_this = static_cast<const TwoStepUpdatedLagrangianVPExplicitElement<TDim>*> (this); */
      /* 	  const TValueType& Val = const_this->GetValue(rVariable); */

      /* 	  for (unsigned int i = 0; i < NumValues; i++){ */
      /* 	    rOutput[i] = Val; */
      /* 	  } */
      /* 	} */
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
      TwoStepUpdatedLagrangianVPExplicitElement & operator=(TwoStepUpdatedLagrangianVPExplicitElement const& rOther);

      /* /// Copy constructor. */
      /* TwoStepUpdatedLagrangianVPExplicitElement(TwoStepUpdatedLagrangianVPExplicitElement const& rOther); */

      ///@}

    }; // Class TwoStepUpdatedLagrangianVPExplicitElement

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  template< unsigned int TDim >
    inline std::istream& operator >>(std::istream& rIStream,
                                     TwoStepUpdatedLagrangianVPExplicitElement<TDim>& rThis)
    {
      return rIStream;
    }

  /// output stream function
  template< unsigned int TDim >
    inline std::ostream& operator <<(std::ostream& rOStream,
                                     const TwoStepUpdatedLagrangianVPExplicitElement<TDim>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }

} // namespace Kratos.

#endif // KRATOS_TWO_STEP_UPDATED_LAGRANGIAN_V_P_EXPLICIT_ELEMENT  defined
