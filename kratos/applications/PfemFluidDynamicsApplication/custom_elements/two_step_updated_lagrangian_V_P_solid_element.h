//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 0.0 $
//
//

#if !defined(KRATOS_TWO_STEP_UPDATED_LAGRANGIAN_V_P_SOLID_ELEMENT_H_INCLUDED )
#define  KRATOS_TWO_STEP_UPDATED_LAGRANGIAN_V_P_SOLID_ELEMENT_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

 
/* // Project includes */
/* #include "containers/array_1d.h" */
/* #include "includes/define.h" */
/* /\* #include "includes/element.h" *\/ */
/* #include "includes/serializer.h" */
/* #include "geometries/geometry.h" */
/* #include "utilities/math_utils.h" */

#include "custom_elements/two_step_updated_lagrangian_V_P_element.h" 

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
/* class TwoStepUpdatedLagrangianVPSolidElement : public Element */
  template< unsigned int TDim > 
  class TwoStepUpdatedLagrangianVPSolidElement : public TwoStepUpdatedLagrangianVPElement<TDim>
    {
  
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of TwoStepUpdatedLagrangianVPSolidElement
      KRATOS_CLASS_POINTER_DEFINITION(TwoStepUpdatedLagrangianVPSolidElement);

      ///base type: 
      typedef TwoStepUpdatedLagrangianVPElement<TDim> BaseType;

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

      typedef VectorMap<IndexType, DataValueContainer> SolutionStepsElementalDataContainerType;

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
    TwoStepUpdatedLagrangianVPSolidElement(IndexType NewId = 0) :
      BaseType(NewId)
      {}

      /// Constructor using an array of nodes.
      /**
       * @param NewId Index of the new element
       * @param ThisNodes An array containing the nodes of the new element
       */
    TwoStepUpdatedLagrangianVPSolidElement(IndexType NewId, const NodesArrayType& ThisNodes) :
      BaseType(NewId, ThisNodes)
	{}

      /// Constructor using a geometry object.
      /**
       * @param NewId Index of the new element
       * @param pGeometry Pointer to a geometry object
       */
    TwoStepUpdatedLagrangianVPSolidElement(IndexType NewId, GeometryType::Pointer pGeometry) :
      BaseType(NewId, pGeometry)
	{}

      /// Constuctor using geometry and properties.
      /**
       * @param NewId Index of the new element
       * @param pGeometry Pointer to a geometry object
       * @param pProperties Pointer to the element's properties
       */
    TwoStepUpdatedLagrangianVPSolidElement(IndexType NewId, GeometryType::Pointer pGeometry, pPropertiesType pProperties) : BaseType(NewId, pGeometry, pProperties)
	{}


   /// copy constructor

    TwoStepUpdatedLagrangianVPSolidElement(TwoStepUpdatedLagrangianVPSolidElement const& rOther):
      BaseType(rOther)
      {}

      /// Destructor.
      virtual ~TwoStepUpdatedLagrangianVPSolidElement()
	{}


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      /// Create a new element of this type
      /**
       * Returns a pointer to a new TwoStepUpdatedLagrangianVPSolidElement element, created using given input
       * @param NewId: the ID of the new element
       * @param ThisNodes: the nodes of the new element
       * @param pProperties: the properties assigned to the new element
       * @return a Pointer to the new element
       */
      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
			      pPropertiesType pProperties) const
	{
	  return Element::Pointer(new TwoStepUpdatedLagrangianVPSolidElement(NewId, BaseType::GetGeometry().Create(ThisNodes), pProperties));
	}

      Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const;


      virtual void Initialize();

      /// Initializes the element and all geometric information required for the problem.
      virtual void InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo);


      virtual void InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo);

      virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
					 ProcessInfo& rCurrentProcessInfo)
      {
	KRATOS_TRY;
	KRATOS_THROW_ERROR(std::logic_error,"TwoStepUpdatedLagrangianVPSolidElement::CalculateLeftHandSide not implemented","");
	KRATOS_CATCH("");
      }

      virtual void CalculateRightHandSide(VectorType& rRightHandSideVector,
					  ProcessInfo& rCurrentProcessInfo)
      {
	KRATOS_TRY;
	KRATOS_THROW_ERROR(std::logic_error,"TwoStepUpdatedLagrangianVPSolidElement::CalculateRightHandSide not implemented","");
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


      virtual void UpdateCauchyStress(unsigned int g);

      virtual void InitializeElementalVariables(ElementalVariables & rElementalVariables);

      /* virtual void CalculateDeltaPosition (Matrix & rDeltaPosition); */

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
      virtual int Check(const ProcessInfo& rCurrentProcessInfo);

      ///@}
      ///@name Inquiry
      ///@{


      ///@}
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
	{
	  std::stringstream buffer;
	  buffer << "TwoStepUpdatedLagrangianVPSolidElement #" << BaseType::Id();
	  return buffer.str();
	}

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
      {
	rOStream << "TwoStepUpdatedLagrangianVPSolidElement" << TDim << "D";
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


      void ComputeMaterialParameters (double& DeviatoricCoeff,
				      double& VolumetricCoeff,
				      double timeStep,
				      const ShapeFunctionsType& rN);


      /// Add integration point contribution to the mass matrix.
      /**
       * A constistent mass matrix is used.
       * @param rMassMatrix The local matrix where the result will be added.
       * @param rN Elemental shape functions.
       * @param Weight Multiplication coefficient for the matrix, typically Density times integration point weight.
       */
   

      void ComputeMeanValueMaterialTangentMatrix(ElementalVariables& rElementalVariables,
						 double& MeanValue,
						 const ShapeFunctionDerivativesType& rShapeDeriv,
						 const double secondLame,
						 const double bulkModulus,
						 const double Weight){};

      void AddCompleteTangentTerm(ElementalVariables& rElementalVariables,
				  MatrixType& rDampingMatrix,
				  const ShapeFunctionDerivativesType& rShapeDeriv,
				  const double secondLame,
				  const double bulkModulus,
				  const double Weight);
	
      void ComputeBulkMatrixForPressureVelLump(MatrixType& BulkVelMatrix,
					       const ShapeFunctionsType& rN,
					       const double Weight);

      void ComputeBulkMatrixForPressureAccLump(MatrixType& BulkAccMatrix,
					       const ShapeFunctionsType& rN,
					       const double Weight);


      void ComputeBulkMatrixForPressureVel(MatrixType& BulkVelMatrix,
				       const ShapeFunctionsType& rN,
				       const double Weight);

      void ComputeBulkMatrixForPressureAcc(MatrixType& BulkAccMatrix,
				       const ShapeFunctionsType& rN,
				       const double Weight);

     void ComputeBoundLHSMatrix(MatrixType& BoundLHSMatrix,
				const ShapeFunctionsType& rN,
				const double Weight){std::cout<<"ComputeBoundLHSMatrix solid"<<std::endl;};

     void ComputeBoundRHSVector(VectorType& BoundRHSVector,
				const ShapeFunctionsType& rN,
				const double Weight){std::cout<<"ComputeBoundRHSvector solid"<<std::endl;};

      void ComputeStabLaplacianMatrix(MatrixType& StabLaplacianMatrix,
				      const ShapeFunctionDerivativesType& rShapeDeriv,
				      const double Weight);
      
      bool CalcMechanicsUpdated(ElementalVariables & rElementalVariables,
				const ProcessInfo& rCurrentProcessInfo,
				unsigned int g,
				const ShapeFunctionsType& N);

	
      void CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables,
					    double TimeStep,
					    unsigned int g,
					    const ShapeFunctionsType& rN);

      virtual void CalculateTauFIC(double& TauOne,
				   double ElemSize,
				   const array_1d< double, 3 > & rAdvVel,
				   const double Density,
				   const double Viscosity,
				   const ProcessInfo& rCurrentProcessInfo);

      void AddStabilizationMatrixLHS(MatrixType& rLeftHandSideMatrix,
					     Matrix& BulkAccMatrix,
					     const ShapeFunctionsType& rN,
					     const double Weight);

      void AddStabilizationNodalTermsLHS(MatrixType& rLeftHandSideMatrix,
						 const double Tau,
						 const double Weight,
						 const ShapeFunctionDerivativesType& rDN_DX,
						 const SizeType i);

      void AddStabilizationNodalTermsRHS(VectorType& rRightHandSideVector,
					 const double Tau,
					 const double Density,
					 const array_1d<double,3> BodyForce,
					 const double Weight,
					 const ShapeFunctionDerivativesType& rDN_DX,
					 const SizeType i);

     /// Write the value of a variable at a point inside the element to a double
      /**
       * Evaluate a nodal variable in the point where the form functions take the
       * values given by rShapeFunc and write the result to rResult.
       * This is an auxiliary function used to compute values in integration points.
       * @param rResult: The variable where the value will be added to
       * @param rVariable: The nodal variable to be read
       * @param rShapeFunc: The values of the form functions in the point
       */
      template< class TVariableType >
	void EvaluateInPoint(TVariableType& rResult,
			     const Kratos::Variable<TVariableType>& Var,
			     const ShapeFunctionsType& rShapeFunc)
	{
	  GeometryType& rGeom = this->GetGeometry();
	  const SizeType NumNodes = rGeom.PointsNumber();

	  rResult = rShapeFunc[0] * rGeom[0].FastGetSolutionStepValue(Var);

	  for(SizeType i = 1; i < NumNodes; i++)
	    {
	      rResult += rShapeFunc[i] * rGeom[i].FastGetSolutionStepValue(Var);
	    }
	}

      /// Write the value of a variable at a point inside the element to a double
      /**
       * Evaluate a nodal variable in the point where the form functions take the
       * values given by rShapeFunc and write the result to rResult.
       * This is an auxiliary function used to compute values in integration points.
       * @param rResult The variable where the value will be added to
       * @param rVariable The nodal variable to be read
       * @param rShapeFunc The values of the form functions in the point
       * @param Step Number of time steps back
       */
      template< class TVariableType >
	void EvaluateInPoint(TVariableType& rResult,
			     const Kratos::Variable<TVariableType>& Var,
			     const ShapeFunctionsType& rShapeFunc,
			     const IndexType Step)
	{
	  GeometryType& rGeom = this->GetGeometry();
	  const SizeType NumNodes = rGeom.PointsNumber();

	  rResult = rShapeFunc[0] * rGeom[0].FastGetSolutionStepValue(Var,Step);

	  for(SizeType i = 1; i < NumNodes; i++)
	    {
	      rResult += rShapeFunc[i] * rGeom[i].FastGetSolutionStepValue(Var,Step);
	    }
	}

      void EvaluateGradientInPoint(array_1d<double,TDim>& rResult,
				   const Kratos::Variable<double>& Var,
				   const ShapeFunctionDerivativesType& rDN_DX)
      {
	GeometryType& rGeom = this->GetGeometry();
	const SizeType NumNodes = rGeom.PointsNumber();
	    
	const double& var = rGeom[0].FastGetSolutionStepValue(Var);
	for (SizeType d = 0; d < TDim; ++d)
	  rResult[d] = rDN_DX(0,d) * var;
			
	for(SizeType i = 1; i < NumNodes; i++)
	  {
	    const double& var = rGeom[i].FastGetSolutionStepValue(Var);
	    for (SizeType d = 0; d < TDim; ++d)
	      rResult[d] += rDN_DX(i,d) * var;
				
	  }
      }

      void EvaluateDivergenceInPoint(double& rResult,
				     const Kratos::Variable< array_1d<double,3> >& Var,
				     const ShapeFunctionDerivativesType& rDN_DX)
      {
	GeometryType& rGeom = this->GetGeometry();
	const SizeType NumNodes = rGeom.PointsNumber();

	rResult = 0.0;
	for(SizeType i = 0; i < NumNodes; i++)
	  {
	    const array_1d<double,3>& var = rGeom[i].FastGetSolutionStepValue(Var);
	    for (SizeType d = 0; d < TDim; ++d)
	      {
		rResult += rDN_DX(i,d) * var[d];
	      }
	  }
      }

      /// Helper function to print results on gauss points
      /** Reads a variable from the element's database and returns it in a format
       * that can be used by GetValueOnIntegrationPoints functions.
       * @see GetValueOnIntegrationPoints
       */
      template<class TValueType>
	void GetElementalValueForOutput(const Kratos::Variable<TValueType>& rVariable,
					std::vector<TValueType>& rOutput)
	{
	  unsigned int NumValues = this->GetGeometry().IntegrationPointsNumber(GeometryData::GI_GAUSS_1);
	  rOutput.resize(NumValues);
	  /*
	    The cast is done to avoid modification of the element's data. Data modification
	    would happen if rVariable is not stored now (would initialize a pointer to &rVariable
	    with associated value of 0.0). This is catastrophic if the variable referenced
	    goes out of scope.
	  */
	  const TwoStepUpdatedLagrangianVPSolidElement<TDim>* const_this = static_cast<const TwoStepUpdatedLagrangianVPSolidElement<TDim>*> (this);
	  const TValueType& Val = const_this->GetValue(rVariable);

	  for (unsigned int i = 0; i < NumValues; i++)
	    rOutput[i] = Val;
	}

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

      virtual void save(Serializer& rSerializer) const
      {
	KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );
      }

      virtual void load(Serializer& rSerializer)
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
      TwoStepUpdatedLagrangianVPSolidElement & operator=(TwoStepUpdatedLagrangianVPSolidElement const& rOther);

      /* /// Copy constructor. */
      /* TwoStepUpdatedLagrangianVPSolidElement(TwoStepUpdatedLagrangianVPSolidElement const& rOther); */

      ///@}

    }; // Class TwoStepUpdatedLagrangianVPSolidElement

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  template< unsigned int TDim >
    inline std::istream& operator >>(std::istream& rIStream,
                                     TwoStepUpdatedLagrangianVPSolidElement<TDim>& rThis)
    {
      return rIStream;
    }

  /// output stream function
  template< unsigned int TDim >
    inline std::ostream& operator <<(std::ostream& rOStream,
                                     const TwoStepUpdatedLagrangianVPSolidElement<TDim>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }

} // namespace Kratos.

#endif // KRATOS_TWO_STEP_UPDATED_LAGRANGIAN_V_P_SOLID_ELEMENT  defined
