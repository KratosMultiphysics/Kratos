//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 0.0 $
//
//

#if !defined(KRATOS_TWO_STEP_UPDATED_LAGRANGIAN_V_P_ELEMENT_H_INCLUDED )
#define  KRATOS_TWO_STEP_UPDATED_LAGRANGIAN_V_P_ELEMENT_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "containers/array_1d.h"
#include "includes/define.h"
#include "includes/element.h"
#include "includes/serializer.h"
#include "geometries/geometry.h"
#include "utilities/math_utils.h"


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
  class TwoStepUpdatedLagrangianVPElement : public Element
    {
  

    protected:

        typedef struct
      {
      	unsigned int voigtsize;
      	// strain state
      	double DetFgrad;
      	double DetFgradVel;
      	double DeviatoricInvariant;
      	double VolumetricDefRate;
      	VectorType SpatialDefRate;
      	VectorType MDGreenLagrangeMaterial;
      	MatrixType Fgrad;
      	MatrixType InvFgrad;
      	MatrixType FgradVel;
      	MatrixType InvFgradVel;
      	MatrixType SpatialVelocityGrad;
      	// Stress state
      	double MeanPressure;
      	VectorType CurrentTotalCauchyStress;
      	VectorType UpdatedTotalCauchyStress;
      	VectorType CurrentDeviatoricCauchyStress;
      	VectorType UpdatedDeviatoricCauchyStress;
            
      } ElementalVariables;
  

      std::vector< Matrix > mOldFgrad;
      std::vector< Vector > mCurrentTotalCauchyStress;
      std::vector< Vector > mCurrentDeviatoricCauchyStress;
      std::vector< Vector > mUpdatedTotalCauchyStress;
      std::vector< Vector > mUpdatedDeviatoricCauchyStress;
      Vector mDetFgrad;

    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of TwoStepUpdatedLagrangianVPElement
      KRATOS_CLASS_POINTER_DEFINITION(TwoStepUpdatedLagrangianVPElement);

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

      /* typedef Element::PropertiesType::Pointer PropertiesType::Pointer; */

      typedef Element::PropertiesType PropertiesType;

      ///@}
      ///@name Life Cycle
      ///@{

      //Constructors.

      /// Default constuctor.
      /**
       * @param NewId Index number of the new element (optional)
       */
    TwoStepUpdatedLagrangianVPElement(IndexType NewId = 0) :
      Element(NewId)
      {}

      /// Constructor using an array of nodes.
      /**
       * @param NewId Index of the new element
       * @param ThisNodes An array containing the nodes of the new element
       */
    TwoStepUpdatedLagrangianVPElement(IndexType NewId, const NodesArrayType& ThisNodes) :
      Element(NewId, ThisNodes)
        {}

      /// Constructor using a geometry object.
      /**
       * @param NewId Index of the new element
       * @param pGeometry Pointer to a geometry object
       */
    TwoStepUpdatedLagrangianVPElement(IndexType NewId, GeometryType::Pointer pGeometry) :
      Element(NewId, pGeometry)
        {}

      /// Constuctor using geometry and properties.
      /**
       * @param NewId Index of the new element
       * @param pGeometry Pointer to a geometry object
       * @param pProperties Pointer to the element's properties
       */
    TwoStepUpdatedLagrangianVPElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
      Element(NewId, pGeometry, pProperties)
        {}

      /// Destructor.
      virtual ~TwoStepUpdatedLagrangianVPElement()
        {}


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      /// Create a new element of this type
      /**
       * Returns a pointer to a new TwoStepUpdatedLagrangianVPElement element, created using given input
       * @param NewId: the ID of the new element
       * @param ThisNodes: the nodes of the new element
       * @param pProperties: the properties assigned to the new element
       * @return a Pointer to the new element
       */
      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
			      PropertiesType::Pointer pProperties) const
        {
	  return Element::Pointer(new TwoStepUpdatedLagrangianVPElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
        }

      virtual void Initialize(){};

      /// Initializes the element and all geometric information required for the problem.
      virtual void InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo){};


      virtual void InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo){};

      /// Calculate the element's local contribution to the system for the current step.
      virtual void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
					VectorType& rRightHandSideVector,
					ProcessInfo& rCurrentProcessInfo);


      virtual void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
					 ProcessInfo& rCurrentProcessInfo)
      {
	KRATOS_TRY;
	KRATOS_THROW_ERROR(std::logic_error,"TwoStepUpdatedLagrangianVPElement::CalculateLeftHandSide not implemented","");
	KRATOS_CATCH("");
      }

      virtual void CalculateRightHandSide(VectorType& rRightHandSideVector,
					  ProcessInfo& rCurrentProcessInfo)
      {
	KRATOS_TRY;
	KRATOS_THROW_ERROR(std::logic_error,"TwoStepUpdatedLagrangianVPElement::CalculateRightHandSide not implemented","");
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
      virtual void EquationIdVector(EquationIdVectorType& rResult,
				    ProcessInfo& rCurrentProcessInfo);

      /// Returns a list of the element's Dofs
      /**
       * @param ElementalDofList the list of DOFs
       * @param rCurrentProcessInfo the current process info instance
       */
      virtual void GetDofList(DofsVectorType& rElementalDofList,
			      ProcessInfo& rCurrentProcessInfo);

    
      virtual GeometryData::IntegrationMethod GetIntegrationMethod() const;
    
      virtual void UpdateCauchyStress(unsigned int g){};

      virtual void InitializeElementalVariables(ElementalVariables & rElementalVariables, unsigned int g){};

      void CalculateDeltaPosition (Matrix & rDeltaPosition);

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
	  buffer << "TwoStepUpdatedLagrangianVPElement #" << Id();
	  return buffer.str();
	}

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
      {
	rOStream << "TwoStepUpdatedLagrangianVPElement" << TDim << "D";
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

      void CalculateLocalMomentumEquations(MatrixType& rLeftHandSideMatrix,
						   VectorType& rRightHandSideVector,
						   ProcessInfo& rCurrentProcessInfo);

      void CalculateLocalContinuityEqForPressure(MatrixType& rLeftHandSideMatrix,
							 VectorType& rRightHandSideVector,
							 ProcessInfo& rCurrentProcessInfo);

      virtual void ComputeMaterialParameters (double& DeviatoricCoeff,
					      double& VolumetricCoeff,
					      double timeStep){};
      
      void VelocityEquationIdVector(EquationIdVectorType& rResult,
				    ProcessInfo& rCurrentProcessInfo);

      void PressureEquationIdVector(EquationIdVectorType& rResult,
				    ProcessInfo& rCurrentProcessInfo);

      void GetVelocityDofList(DofsVectorType& rElementalDofList,
			      ProcessInfo& rCurrentProcessInfo);

      void GetPressureDofList(DofsVectorType& rElementalDofList,
			      ProcessInfo& rCurrentProcessInfo);


      void CalcMeanVelocity(double& meanVelocity,
			    const int Step);
      
      void CalcMeanPressure(double& meanPressure,
			    const int Step);
      

      void GetPressureValues(Vector& rValues,
			     const int Step = 0);

      void GetVelocityValues(Vector& rValues,
			     const int Step = 0);

      void GetPositions(Vector& rValues);

      void GetUpdatedPositions(Vector& rValues,
			       const ProcessInfo& rCurrentProcessInfo);

      void GetAccelerationValues(Vector& rValues,
				 const int Step = 0);

      /// Determine integration point weights and shape funcition derivatives from the element's geometry.
      void CalculateGeometryData(ShapeFunctionDerivativesArrayType& rDN_DX,
				 Matrix& rNContainer,
				 Vector& rGaussWeights);

      double ElementSize(/*ShapeFunctionDerivativesType& rDN_DX*/);


      /**
       * @brief EquivalentStrainRate Calculate the second invariant of the strain rate tensor GammaDot = (2SijSij)^0.5.
       *
       * @note Our implementation of non-Newtonian consitutive models such as Bingham relies on this funcition being
       * defined on all fluid elements.
       *
       * @param rDN_DX Shape function derivatives at the integration point.
       * @return GammaDot = (2SijSij)^0.5.
       */
      double EquivalentStrainRate(const ShapeFunctionDerivativesType &rDN_DX) const;


      /// Add integration point contribution to the mass matrix.
      /**
       * A constistent mass matrix is used.
       * @param rMassMatrix The local matrix where the result will be added.
       * @param rN Elemental shape functions.
       * @param Weight Multiplication coefficient for the matrix, typically Density times integration point weight.
       */
      void ComputeMomentumMassTerm(Matrix& rMassMatrix,
					   const ShapeFunctionsType& rN,
					   const double Weight);

      void ComputeLumpedMassMatrix(Matrix& rMassMatrix,
					   const double Weight);

    
      void AddExternalForces( Vector& rRHSVector,
				      const double Density,
				      const array_1d<double,3>& rBodyForce,
				      const ShapeFunctionsType& rN,
				      const double Weight);
      
      void AddInternalForces( Vector& rRHSVector,
				      const ShapeFunctionDerivativesType& rDN_DX,
				      ElementalVariables& rElementalVariables,
				      const double Weight);

      void AddDynamicForces(Vector& rRHSVector,
				     const double Weight,
				     const double TimeStep);

      void AddDeviatoricInternalForces( Vector& rRHSVector,
						const ShapeFunctionDerivativesType& rDN_DX,
						ElementalVariables& rElementalVariables,
						const double Weight);

      virtual void AddCompleteTangentTerm(ElementalVariables& rElementalVariables,
					  MatrixType& rDampingMatrix,
					  const ShapeFunctionDerivativesType& rShapeDeriv,
					  const double secondLame,
					  const double bulkModulus,
					  const double Weight){};
	
      virtual void ComputeBulkMatrixForPressureVel(MatrixType& BulkVelMatrix,
						   const ShapeFunctionsType& rN,
						   const double Weight){};
      
      virtual void ComputeBulkMatrixForPressureAcc(MatrixType& BulkAccMatrix,
						   const ShapeFunctionsType& rN,
						   const double Weight){};

      virtual void ComputeStabLaplacianMatrix(MatrixType& StabLaplacianMatrix,
					      const ShapeFunctionDerivativesType& rShapeDeriv,
					      const double Weight){};
      
      virtual void CalcMechanicsUpdated(ElementalVariables & rElementalVariables,
					const ProcessInfo& rCurrentProcessInfo,
					unsigned int g){};

      void CalcStrainRateUpdated(ElementalVariables & rElementalVariables,
					 const ProcessInfo& rCurrentProcessInfo,
					 unsigned int g);

      void CalcVelDefGrad(const ShapeFunctionDerivativesType& rDN_DX,
				  MatrixType &FgradVel,
				  MatrixType &invFgradVel,
				  double &FVelJacobian);

      void CalcFGrad(const ShapeFunctionDerivativesType& rDN_DX,
			     MatrixType &Fgrad,
			     MatrixType &invFgrad,
			     double &FJacobian,
			     const ProcessInfo& rCurrentProcessInfo);

      void CalcVolumetricDefRate(const ShapeFunctionDerivativesType& rDN_DX,
					 double &volumetricDefRate,
					 MatrixType &invGradDef);

      void CalcVolDefRateFromSpatialVelGrad(double &volumetricDefRate,
						    MatrixType &SpatialVelocityGrad);


     void CalcSpatialVelocityGrad(MatrixType &invFgrad,
					   MatrixType &VelDefgrad,
					   MatrixType &SpatialVelocityGrad);

      void CalcMDGreenLagrangeMaterial(MatrixType &Fgrad,
					       MatrixType &VelDefgrad, 
					       VectorType &MDGreenLagrangeMaterial);

      void CalcSpatialDefRate(VectorType &MDGreenLagrangeMaterial,
				      MatrixType &invFgrad,
				      VectorType &SpatialDefRate);

      void CalcDeviatoricInvariant(VectorType &SpatialDefRate,
					   double &DeviatoricInvariant);

      void CheckStrain1(double &VolumetricDefRate,
				MatrixType &SpatialVelocityGrad);

      void CheckStrain2(MatrixType &SpatialVelocityGrad,
				 MatrixType &Fgrad,
				 MatrixType &VelDefgrad);
	
      void CheckStrain3(VectorType &SpatialDefRate,
				MatrixType &SpatialVelocityGrad);
	
      virtual void CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables,
						    double TimeStep,
						    unsigned int g){};

      virtual void CalculateTauFIC(double& TauOne,
				   double ElemSize,
				   const array_1d< double, 3 > & rAdvVel,
				   const double Density,
				   const double Viscosity,
				   const ProcessInfo& rCurrentProcessInfo){};

      virtual void AddStabilizationMatrixLHS(MatrixType& rLeftHandSideMatrix,
					     Matrix& BulkAccMatrix,
					     const ShapeFunctionsType& rN,
					     const double Weight){};

      virtual void AddStabilizationNodalTermsLHS(MatrixType& rLeftHandSideMatrix,
						 const double Tau,
						 const double Weight,
						 const ShapeFunctionDerivativesType& rDN_DX,
						 const SizeType i){};

      virtual void AddStabilizationNodalTermsRHS(VectorType& rRightHandSideVector,
						 const double Tau,
						 const double Density,
						 const array_1d<double,3> BodyForce,
						 const double Weight,
						 const ShapeFunctionDerivativesType& rDN_DX,
						 const SizeType i){};

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
	  const TwoStepUpdatedLagrangianVPElement<TDim>* const_this = static_cast<const TwoStepUpdatedLagrangianVPElement<TDim>*> (this);
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
      TwoStepUpdatedLagrangianVPElement & operator=(TwoStepUpdatedLagrangianVPElement const& rOther);

      /// Copy constructor.
      TwoStepUpdatedLagrangianVPElement(TwoStepUpdatedLagrangianVPElement const& rOther);

      ///@}

    }; // Class TwoStepUpdatedLagrangianVPElement

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  template< unsigned int TDim >
    inline std::istream& operator >>(std::istream& rIStream,
                                     TwoStepUpdatedLagrangianVPElement<TDim>& rThis)
    {
      return rIStream;
    }

  /// output stream function
  template< unsigned int TDim >
    inline std::ostream& operator <<(std::ostream& rOStream,
                                     const TwoStepUpdatedLagrangianVPElement<TDim>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }

} // namespace Kratos.

#endif // KRATOS_TWO_STEP_UPDATED_LAGRANGIAN_V_P_ELEMENT  defined
