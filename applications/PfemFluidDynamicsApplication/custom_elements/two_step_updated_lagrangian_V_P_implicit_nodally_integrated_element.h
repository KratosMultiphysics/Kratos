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

      /// Calculate the element's local contribution to the system for the current step.
      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
				VectorType& rRightHandSideVector,
				ProcessInfo& rCurrentProcessInfo) override;


      void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
				 ProcessInfo& rCurrentProcessInfo) override;

      void CalculateRightHandSide(VectorType& rRightHandSideVector,
				  ProcessInfo& rCurrentProcessInfo) override
      {
	KRATOS_TRY;
	KRATOS_THROW_ERROR(std::logic_error,"TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement::CalculateRightHandSide not implemented","");
	KRATOS_CATCH("");
      }

 
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

      void CalcElementalStrains(ElementalVariables & rElementalVariables,
				const ProcessInfo &rCurrentProcessInfo,
				const ShapeFunctionDerivativesType& rDN_DX);
            
      void CalculateGeometryData(ShapeFunctionDerivativesArrayType& rDN_DX,
				 Matrix& rNContainer,
				 Vector& rGaussWeights);
      
      void UpdateCauchyStress(unsigned int g,ProcessInfo& rCurrentProcessInfo) override{};

      void InitializeElementalVariables(ElementalVariables & rElementalVariables) override{
	KRATOS_TRY;
	KRATOS_THROW_ERROR(std::logic_error,"InitializeElementalVariables","");
	KRATOS_CATCH("");
      };

      void ComputeExternalForces(Vector& rRHSVector,const double Density,const double Weight);

      void ComputeInternalForces(Vector& rRHSVector,
				 const ShapeFunctionDerivativesType& rDN_DX,
				 ElementalVariables& rElementalVariables,
				 const double Weight);

      /* void CalculateDeltaPosition (Matrix & rDeltaPosition); */

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
	  buffer << "TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement #" <<  BaseType::Id();
	  return buffer.str();
	}

      /// Print information about this object.
      void PrintInfo(std::ostream& rOStream) const override
      {
	rOStream << "TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement" << TDim << "D";
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

      void NodalFreeSurfaceLength(unsigned int nodeIndex);

      void GetValueOnIntegrationPoints( const Variable<double>& rVariable,
					std::vector<double>& rValues,
					const ProcessInfo& rCurrentProcessInfo ) override;
      
      void CalculateLocalMomentumEquations(MatrixType& rLeftHandSideMatrix,
					   VectorType& rRightHandSideVector,
					   ProcessInfo& rCurrentProcessInfo) override;

      void CalculateLocalContinuityEqForPressure(MatrixType& rLeftHandSideMatrix,
						 VectorType& rRightHandSideVector,
						 ProcessInfo& rCurrentProcessInfo) override{};

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

      double GetThetaMomentum () override{
	std::cout<<"I SHOULD NOT ENTER HERE!"<<std::endl;
	return 1.0;
      };

      double GetThetaContinuity () override{
	std::cout<<"I SHOULD NOT ENTER HERE!"<<std::endl;
	return 1.0;
      };
      
  

      void GetPositions(Vector& rValues,
			const ProcessInfo& rCurrentProcessInfo,
			const double theta) override{};

      void ComputeMeanValueMaterialTangentMatrix(ElementalVariables& rElementalVariables,
						 double& MeanValue,
						 const ShapeFunctionDerivativesType& rShapeDeriv,
						 const double secondLame,
						 double& bulkModulus,
						 const double Weight,
						 double& MeanValueMass,
						 const double TimeStep){};

      virtual void ComputeBulkReductionCoefficient(MatrixType MassMatrix,
						   MatrixType StiffnessMatrix,
						   double& meanValueStiff,
						   double& bulkCoefficient,
						   double timeStep){};

      void ComputeCompleteTangentTerm(ElementalVariables& rElementalVariables,
				      MatrixType& rDampingMatrix,
				      const ShapeFunctionDerivativesType& rShapeDeriv,
				      const double secondLame,
				      const double bulkModulus,
				      const double theta,
				      const double Weight);

      virtual void ComputeBulkMatrixLump(MatrixType& BulkMatrix,
					 const double Weight){};
      
      virtual void ComputeBulkMatrixConsistent(MatrixType& BulkMatrix,
					       const double Weight){};
      
      virtual void ComputeBulkMatrix(MatrixType& BulkMatrix,
				     const ShapeFunctionsType& rN,
				     const double Weight){};
	
      /* virtual void ComputeBulkMatrixForPressureVelLump(MatrixType& BulkVelMatrix, */
      /* 						   const double Weight){}; */
      
      /* virtual void ComputeBulkMatrixForPressureAccLump(MatrixType& BulkAccMatrix, */
      /* 						   const double Weight){}; */

      /* virtual void ComputeBulkMatrixForPressureVel(MatrixType& BulkVelMatrix, */
      /* 						   const ShapeFunctionsType& rN, */
      /* 						   const double Weight){}; */
      
      /* virtual void ComputeBulkMatrixForPressureAcc(MatrixType& BulkAccMatrix, */
      /* 						   const ShapeFunctionsType& rN, */
      /* 						   const double Weight){}; */

      virtual void ComputeBoundLHSMatrix(MatrixType& BoundLHSMatrix,
					 const ShapeFunctionsType& rN,
					 const double Weight){};

      virtual void ComputeBoundRHSVector(VectorType& BoundRHSVector,
					 const ShapeFunctionsType& rN,
					 const double TimeStep,
					 const double BoundRHSCoeffAcc,
					 const double BoundRHSCoeffDev){};

      virtual void ComputeStabLaplacianMatrix(MatrixType& StabLaplacianMatrix,
					      const ShapeFunctionDerivativesType& rShapeDeriv,
					      const double Weight){};
  
      void CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables,
					    double TimeStep,
					    unsigned int g) override{};

      virtual void CalculateTauFIC(double& TauOne,
				   double ElemSize,
				   const double Density,
				   const double Viscosity,
				   const ProcessInfo& rCurrentProcessInfo){};

      virtual void AddStabilizationMatrixLHS(MatrixType& rLeftHandSideMatrix,
					     Matrix& BulkAccMatrix,
					     const ShapeFunctionsType& rN,
					     const double Weight){};

      virtual void AddStabilizationNodallyTermsLHS(MatrixType& rLeftHandSideMatrix,
						 const double Tau,
						 const double Weight,
						 const ShapeFunctionDerivativesType& rDN_DX,
						 const SizeType i){};

      virtual void AddStabilizationNodallyTermsRHS(VectorType& rRightHandSideVector,
						 const double Tau,
						 const double Density,
						 const double Weight,
						 const ShapeFunctionDerivativesType& rDN_DX,
						 const SizeType i){};

   
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
