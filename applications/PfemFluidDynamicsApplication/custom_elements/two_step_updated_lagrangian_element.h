//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:              April 2018 $
//   Revision:            $Revision:                 0.0 $
//
//

#if !defined(KRATOS_TWO_STEP_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED )
#define  KRATOS_TWO_STEP_UPDATED_LAGRANGIAN_ELEMENT_H_INCLUDED

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
#include "utilities/geometry_utilities.h"

#include "pfem_fluid_dynamics_application_variables.h"

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
    class TwoStepUpdatedLagrangianElement : public Element
    {
  
    protected:


      ///@name Protected static Member Variables
      ///@{


      ///@}
      ///@name Protected member Variables
      ///@{
      typedef struct
      {
      	unsigned int voigtsize;
      	// strain state
      	double DetFgrad;
      	double DetFgradVel;
      	double DeviatoricInvariant;
      	double EquivalentStrainRate;
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
  

    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of TwoStepUpdatedLagrangianElement
      KRATOS_CLASS_POINTER_DEFINITION(TwoStepUpdatedLagrangianElement);

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
    TwoStepUpdatedLagrangianElement(IndexType NewId = 0) :
      Element(NewId)
      {}

      /// Constructor using an array of nodes.
      /**
       * @param NewId Index of the new element
       * @param ThisNodes An array containing the nodes of the new element
       */
    TwoStepUpdatedLagrangianElement(IndexType NewId, const NodesArrayType& ThisNodes) :
      Element(NewId, ThisNodes)
        {}

      /// Constructor using a geometry object.
      /**
       * @param NewId Index of the new element
       * @param pGeometry Pointer to a geometry object
       */
    TwoStepUpdatedLagrangianElement(IndexType NewId, GeometryType::Pointer pGeometry) :
      Element(NewId, pGeometry)
        {}

      /// Constuctor using geometry and properties.
      /**
       * @param NewId Index of the new element
       * @param pGeometry Pointer to a geometry object
       * @param pProperties Pointer to the element's properties
       */
    TwoStepUpdatedLagrangianElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) :
      Element(NewId, pGeometry, pProperties)
        {}



      /// copy constructor

      TwoStepUpdatedLagrangianElement(TwoStepUpdatedLagrangianElement const& rOther);
      
      /// Destructor.
      virtual ~TwoStepUpdatedLagrangianElement()
        {}


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      /// Create a new element of this type
      /**
       * Returns a pointer to a new TwoStepUpdatedLagrangianElement element, created using given input
       * @param NewId: the ID of the new element
       * @param ThisNodes: the nodes of the new element
       * @param pProperties: the properties assigned to the new element
       * @return a Pointer to the new element
       */
      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
			      PropertiesType::Pointer pProperties) const override
        {
	  return Element::Pointer(new TwoStepUpdatedLagrangianElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
        }


      Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override ;

      void Initialize() override{};

      /// Initializes the element and all geometric information required for the problem.
      void InitializeSolutionStep(ProcessInfo &rCurrentProcessInfo) override{};


      void InitializeNonLinearIteration(ProcessInfo &rCurrentProcessInfo) override{};

      /// Calculate the element's local contribution to the system for the current step.
      void CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
				VectorType& rRightHandSideVector,
				ProcessInfo& rCurrentProcessInfo) override{};


      void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
				 ProcessInfo& rCurrentProcessInfo) override
      {
	KRATOS_TRY;
	KRATOS_THROW_ERROR(std::logic_error,"TwoStepUpdatedLagrangianElement::CalculateLeftHandSide not implemented","");
	KRATOS_CATCH("");
      }

      void CalculateRightHandSide(VectorType& rRightHandSideVector,
				  ProcessInfo& rCurrentProcessInfo) override{};

      
      /* virtual void CalculateRightHandSideMomentum(VectorType& rRightHandSideVector, */
      /* 						  ProcessInfo& rCurrentProcessInfo){}; */
 
      // The following methods have different implementations depending on TDim
      /// Provides the global indices for each one of this element's local rows
      /**
       * this determines the elemental equation ID vector for all elemental
       * DOFs
       * @param rResult A vector containing the global Id of each row
       * @param rCurrentProcessInfo the current process info object (unused)
       */
      void EquationIdVector(EquationIdVectorType& rResult,
			    ProcessInfo& rCurrentProcessInfo) override;

      /// Returns a list of the element's Dofs
      /**
       * @param ElementalDofList the list of DOFs
       * @param rCurrentProcessInfo the current process info instance
       */
      void GetDofList(DofsVectorType& rElementalDofList,
		      ProcessInfo& rCurrentProcessInfo) override;

    
      GeometryData::IntegrationMethod GetIntegrationMethod() const override;
    
      virtual void UpdateCauchyStress(unsigned int g,ProcessInfo& rCurrentProcessInfo){};

      virtual void InitializeElementalVariables(ElementalVariables & rElementalVariables){};

      void CalculateDeltaPosition (Matrix & rDeltaPosition);

      void AddExplicitContribution(const VectorType& rRHSVector, 
				   const Variable<VectorType>& rRHSVariable, 
				   Variable<array_1d<double,3> >& rDestinationVariable, 
				   const ProcessInfo& rCurrentProcessInfo) override{};

      void AddExplicitContribution(const VectorType& rRHSVector, 
				   const Variable<VectorType>& rRHSVariable, 
				   Variable<double > & rDestinationVariable, 
				   const ProcessInfo& rCurrentProcessInfo) override {};
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
      int Check(const ProcessInfo& rCurrentProcessInfo) override {return 0;};

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
	  buffer << "TwoStepUpdatedLagrangianElement #" << Id();
	  return buffer.str();
	}

      /// Print information about this object.
      void PrintInfo(std::ostream& rOStream) const  override{}

      //        /// Print object's data.
      //        virtual void PrintData(std::ostream& rOStream) const;

      ///@}
      ///@name Friends
      ///@{

      ///@}
    protected:

      double mMaterialDeviatoricCoefficient=0;
      double mMaterialVolumetricCoefficient=0;
      double mMaterialDensity=0;

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
					const ProcessInfo& rCurrentProcessInfo ) override{};
      
      virtual void CalculateLocalMomentumEquations(MatrixType& rLeftHandSideMatrix,
						   VectorType& rRightHandSideVector,
						   ProcessInfo& rCurrentProcessInfo){};

      virtual void CalculateLocalContinuityEqForPressure(MatrixType& rLeftHandSideMatrix,
							 VectorType& rRightHandSideVector,
							 ProcessInfo& rCurrentProcessInfo){};
      
      virtual void CalculateExplicitContinuityEquation(MatrixType& rLeftHandSideMatrix,
						       VectorType& rRightHandSideVector,
						       ProcessInfo& rCurrentProcessInfo){};


      virtual void ComputeMaterialParameters (double& Density,
					      double& DeviatoricCoeff,
					      double& VolumetricCoeff,
					      ProcessInfo& rCurrentProcessInfo,
					      ElementalVariables& rElementalVariables){};

      virtual void ComputeMaterialParametersGranularGas (double& Density,
							 double& DeviatoricCoeff,
							 double& VolumetricCoeff,
							 ProcessInfo& rCurrentProcessInfo,
							 ElementalVariables& rElementalVariables){};

      virtual double GetThetaMomentum (){
	std::cout<<"I SHOULD NOT ENTER HERE!"<<std::endl;
	return 0.0;
      };

      virtual double GetThetaContinuity (){
	std::cout<<"I SHOULD NOT ENTER HERE!"<<std::endl;
	return 0.0;
      };
      
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
      
      void GetDensityValues(Vector& rValues,
			    const int Step = 0);

      void GetVelocityValues(Vector& rValues,
			     const int Step = 0);

      void GetDisplacementValues(Vector& rValues,
				 const int Step = 0);

      virtual void GetPositions(Vector& rValues,
				const ProcessInfo& rCurrentProcessInfo,
				const double theta){};

      void GetAccelerationValues(Vector& rValues,
				 const int Step = 0);


      void GetElementalAcceleration(Vector& rValues,
				    const int Step,
				    const double TimeStep);

      /// Determine integration point weights and shape funcition derivatives from the element's geometry.
      void CalculateGeometryData(ShapeFunctionDerivativesArrayType& rDN_DX,
				 Matrix& rNContainer,
				 Vector& rGaussWeights);

      void CalculateGeometryData(Vector& rGaussWeights);
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
      void CalculateMassMatrix(Matrix& rMassMatrix,
			       ProcessInfo& rCurrentProcessInfo) override{};
      
      /* virtual void CalculateMassMatrixMomentum(Matrix& rMassMatrix, */
      /* 				       ProcessInfo& rCurrentProcessInfo){}; */
      
      /* void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix, */
      /* 				    VectorType& rRightHandSideVector, */
      /* 				    Flags& rCalculationFlags); */

      void ComputeMassMatrix(Matrix& rMassMatrix,
			     const ShapeFunctionsType& rN,
			     const double Weight,
			     double& MeanValue);
      
      void ComputeLumpedMassMatrix(Matrix& rMassMatrix,
				   const double Weight,
				   double& MeanValue);

    
      void AddExternalForces( Vector& rRHSVector,
			      const double Density,
			      const ShapeFunctionsType& rN,
			      const double Weight);
      
      void AddInternalForces( Vector& rRHSVector,
			      const ShapeFunctionDerivativesType& rDN_DX,
			      ElementalVariables& rElementalVariables,
			      const double Weight);

      void ComputeBulkMatrixLump(MatrixType& BulkMatrix,
				 const double Weight);
      
      void ComputeBulkMatrixConsistent(MatrixType& BulkMatrix,
				       const double Weight);
      
      void ComputeBulkMatrix(MatrixType& BulkMatrix,
			     const ShapeFunctionsType& rN,
			     const double Weight);
       
      virtual void ComputeBulkMatrixRHS(MatrixType& BulkMatrix,
					const double Weight){};
       
      
      bool CalcMechanicsUpdated(ElementalVariables & rElementalVariables,
				const ProcessInfo& rCurrentProcessInfo,
				const ShapeFunctionDerivativesType& rDN_DX,
				unsigned int g);

      bool CalcStrainRate(ElementalVariables & rElementalVariables,
			  const ProcessInfo& rCurrentProcessInfo,
			  const ShapeFunctionDerivativesType& rDN_DX,
			  const double theta);

      bool CalcCompleteStrainRate(ElementalVariables & rElementalVariables,
				  const ProcessInfo& rCurrentProcessInfo,
				  const ShapeFunctionDerivativesType& rDN_DX,
				  const double theta);

      void CalcVelDefGrad(const ShapeFunctionDerivativesType& rDN_DX,
			  MatrixType &FgradVel,
			  MatrixType &invFgradVel,
			  double &FVelJacobian,
			  const double theta);

      void CalcFGrad(const ShapeFunctionDerivativesType& rDN_DX,
		     MatrixType &Fgrad,
		     MatrixType &invFgrad,
		     double &FJacobian,
		     const ProcessInfo& rCurrentProcessInfo,
		     const double theta);

      void CalcVolumetricDefRate(const ShapeFunctionDerivativesType& rDN_DX,
				 double &volumetricDefRate,
				 MatrixType &invGradDef,
				 const double theta);

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

      void CalcEquivalentStrainRate(VectorType &SpatialDefRate,
				    double &EquivalentStrainRate);

      double CalcNormalProjectionDefRate(const VectorType &SpatialDefRate, 
					 const array_1d<double, 3> NormalVector); 
      
      double CalcNormalProjectionDefRate(VectorType &SpatialDefRate);

      void CheckStrain1(double &VolumetricDefRate,
			MatrixType &SpatialVelocityGrad);

      void CheckStrain2(MatrixType &SpatialVelocityGrad,
			MatrixType &Fgrad,
			MatrixType &VelDefgrad);
	
      bool CheckStrain3(VectorType &SpatialDefRate,
			MatrixType &SpatialVelocityGrad);

	
      virtual void CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables,
						    double TimeStep,
						    unsigned int g){};


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

      template< class TVariableType >
	void EvaluatePropertyFromANotRigidNode(TVariableType& rResult,
					       const Kratos::Variable<TVariableType>& Var)
	{
	  GeometryType& rGeom = this->GetGeometry();
	  const SizeType NumNodes = rGeom.PointsNumber();
	  bool found=false;
	  for(SizeType i = 0; i < NumNodes; i++)
	    {
	      if(rGeom[i].IsNot(RIGID)){
		rResult = rGeom[i].FastGetSolutionStepValue(Var);
		found=true;
		break;
	      }
	    }
	  if(found==false)
	    std::cout<<"ATTENTION! NO PROPERTIES HAVE BEEN FOUNDED FOR THIS ELEMENT! X,Y "<<rGeom[1].X()<<","<<rGeom[1].Y()<<std::endl;
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
	  /* unsigned int NumValues = this->GetGeometry().IntegrationPointsNumber(GeometryData::GI_GAUSS_4); */
	  rOutput.resize(NumValues);
	  /*
	    The cast is done to avoid modification of the element's data. Data modification
	    would happen if rVariable is not stored now (would initialize a pointer to &rVariable
	    with associated value of 0.0). This is catastrophic if the variable referenced
	    goes out of scope.
	  */
	  const TwoStepUpdatedLagrangianElement<TDim>* const_this = static_cast<const TwoStepUpdatedLagrangianElement<TDim>*> (this);
	  const TValueType& Val = const_this->GetValue(rVariable);

	  for (unsigned int i = 0; i < NumValues; i++){
	    rOutput[i] = Val;
	  }
	}


       
      void GetOutwardsUnitNormalForTwoPoints(array_1d<double, 3> &NormalVector, 
					     unsigned int idA, 
					     unsigned int idB, 
					     unsigned int otherId) 
      { 
	GeometryType& rGeom = this->GetGeometry(); 
	double deltaX= rGeom[idB].X()-rGeom[idA].X(); 
	double deltaY= rGeom[idB].Y()-rGeom[idA].Y(); 
	double elementSize=sqrt(deltaX*deltaX+deltaY*deltaY); // this is just to have an idea of the size of the element 
	if(fabs(deltaX)>fabs(deltaY)){//to avoid division by zero or small numbers 
	  NormalVector[0]=-deltaY/deltaX; 
	  NormalVector[1]=1.0; 
	  double normNormal=NormalVector[0]*NormalVector[0] + NormalVector[1]*NormalVector[1]; 
	  NormalVector*=1.0/sqrt(normNormal); 
	}else{ 
	  NormalVector[0]=1.0; 
	  NormalVector[1]=-deltaX/deltaY; 
	  double normNormal=NormalVector[0]*NormalVector[0] + NormalVector[1]*NormalVector[1]; 
	  NormalVector*=1.0/sqrt(normNormal); 
	} 
	//to determine if the computed normal outwards or inwards 
	const array_1d<double, 3> MeanPoint=(rGeom[idB].Coordinates()+rGeom[idA].Coordinates())*0.5; 
	const array_1d<double, 3> DistanceA=rGeom[otherId].Coordinates() - (MeanPoint + NormalVector*elementSize); 
	const array_1d<double, 3> DistanceB=rGeom[otherId].Coordinates() - (MeanPoint - NormalVector*elementSize); 
	const double normA=DistanceA[0]*DistanceA[0] + DistanceA[1]*DistanceA[1]; 
	const double normB=DistanceB[0]*DistanceB[0] + DistanceB[1]*DistanceB[1]; 
	if(normB>normA){ 
	  NormalVector*=-1.0; 
	}   
      } 

      void GetOutwardsUnitNormalForThreePoints(array_1d<double, 3> &NormalVector, 
					       unsigned int idA, 
					       unsigned int idB, 
					       unsigned int idC, 
					       unsigned int otherId) 
      { 
	GeometryType& rGeom = this->GetGeometry(); 
	const array_1d<double, 3> TangentXi  = rGeom[idB].Coordinates() - rGeom[idA].Coordinates(); 
        const array_1d<double, 3> TangentEta = rGeom[idC].Coordinates() - rGeom[idA].Coordinates(); 
 
        MathUtils<double>::CrossProduct(NormalVector, TangentXi, TangentEta); 
	double normNormal=NormalVector[0]*NormalVector[0] + NormalVector[1]*NormalVector[1] + NormalVector[2]*NormalVector[2]; 
	NormalVector*=1.0/sqrt(normNormal); 
   
	//to determine if the computed normal outwards or inwards 
	double deltaX= rGeom[idB].X()-rGeom[idA].X(); 
	double deltaY= rGeom[idB].Y()-rGeom[idA].Y(); 
	double deltaZ= rGeom[idB].Z()-rGeom[idA].Z(); 
	double elementSize=sqrt(deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ); // this is just to have an idea of the size of the element 
	const array_1d<double, 3> MeanPoint=(rGeom[idC].Coordinates()+rGeom[idB].Coordinates()+rGeom[idA].Coordinates())*0.3333333333; 
	const array_1d<double, 3> DistanceA=rGeom[otherId].Coordinates() - (MeanPoint + NormalVector*elementSize); 
	const array_1d<double, 3> DistanceB=rGeom[otherId].Coordinates() - (MeanPoint - NormalVector*elementSize); 
	const double normA=DistanceA[0]*DistanceA[0] + DistanceA[1]*DistanceA[1] + DistanceA[2]*DistanceA[2]; 
	const double normB=DistanceB[0]*DistanceB[0] + DistanceB[1]*DistanceB[1] + DistanceB[2]*DistanceB[2]; 
	if(normB>normA){ 
	  NormalVector*=-1.0; 
	} 
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
      TwoStepUpdatedLagrangianElement & operator=(TwoStepUpdatedLagrangianElement const& rOther);

      /* /// Copy constructor. */
      /* TwoStepUpdatedLagrangianElement(TwoStepUpdatedLagrangianElement const& rOther); */

      ///@}

    }; // Class TwoStepUpdatedLagrangianElement

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  template< unsigned int TDim >
    inline std::istream& operator >>(std::istream& rIStream,
                                     TwoStepUpdatedLagrangianElement<TDim>& rThis)
    {
      return rIStream;
    }

  /// output stream function
  template< unsigned int TDim >
    inline std::ostream& operator <<(std::ostream& rOStream,
                                     const TwoStepUpdatedLagrangianElement<TDim>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }

} // namespace Kratos.

#endif // KRATOS_TWO_STEP_UPDATED_LAGRANGIAN_ELEMENT  defined
