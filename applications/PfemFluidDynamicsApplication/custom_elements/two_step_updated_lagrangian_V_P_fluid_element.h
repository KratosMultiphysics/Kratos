//
//   Project Name:        KratosFluidDynamicsApplication $
//   Last modified by:    $Author:               AFranci $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 0.0 $
//
//

#if !defined(KRATOS_TWO_STEP_UPDATED_LAGRANGIAN_V_P_FLUID_ELEMENT_H_INCLUDED )
#define  KRATOS_TWO_STEP_UPDATED_LAGRANGIAN_V_P_FLUID_ELEMENT_H_INCLUDED

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
  template< unsigned int TDim >
    class TwoStepUpdatedLagrangianVPFluidElement : public TwoStepUpdatedLagrangianVPElement<TDim>
    {

    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of TwoStepUpdatedLagrangianVPFluidElement
      KRATOS_CLASS_POINTER_DEFINITION(TwoStepUpdatedLagrangianVPFluidElement);

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
    TwoStepUpdatedLagrangianVPFluidElement(IndexType NewId = 0) :
      BaseType(NewId)
      {}

      /// Constructor using an array of nodes.
      /**
       * @param NewId Index of the new element
       * @param ThisNodes An array containing the nodes of the new element
       */
    TwoStepUpdatedLagrangianVPFluidElement(IndexType NewId, const NodesArrayType& ThisNodes) :
      BaseType(NewId, ThisNodes)
	{}

      /// Constructor using a geometry object.
      /**
       * @param NewId Index of the new element
       * @param pGeometry Pointer to a geometry object
       */
    TwoStepUpdatedLagrangianVPFluidElement(IndexType NewId, GeometryType::Pointer pGeometry) :
      BaseType(NewId, pGeometry)
	{}

      /// Constuctor using geometry and properties.
      /**
       * @param NewId Index of the new element
       * @param pGeometry Pointer to a geometry object
       * @param pProperties Pointer to the element's properties
       */
    TwoStepUpdatedLagrangianVPFluidElement(IndexType NewId, GeometryType::Pointer pGeometry, pPropertiesType pProperties) :
      BaseType(NewId, pGeometry, pProperties)
	{}


     /// copy constructor

    TwoStepUpdatedLagrangianVPFluidElement(TwoStepUpdatedLagrangianVPFluidElement const& rOther):
        BaseType(rOther)
      {}
 

      /// Destructor.
      virtual ~TwoStepUpdatedLagrangianVPFluidElement()
	{}


      ///@}
      ///@name Operators
      ///@{


      ///@}
      ///@name Operations
      ///@{

      /// Create a new element of this type
      /**
       * Returns a pointer to a new TwoStepUpdatedLagrangianVPFluidElement element, created using given input
       * @param NewId: the ID of the new element
       * @param ThisNodes: the nodes of the new element
       * @param pProperties: the properties assigned to the new element
       * @return a Pointer to the new element
       */
      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,
			      pPropertiesType pProperties) const
	{
	  return Element::Pointer(new TwoStepUpdatedLagrangianVPFluidElement(NewId, BaseType::GetGeometry().Create(ThisNodes), pProperties));
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
	KRATOS_THROW_ERROR(std::logic_error,"TwoStepUpdatedLagrangianVPFluidElement::CalculateLeftHandSide not implemented","");
	KRATOS_CATCH("");
      }

      virtual void CalculateRightHandSide(VectorType& rRightHandSideVector,
					  ProcessInfo& rCurrentProcessInfo)
      {
	KRATOS_TRY;
	KRATOS_THROW_ERROR(std::logic_error,"TwoStepUpdatedLagrangianVPFluidElement::CalculateRightHandSide not implemented","");
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

      virtual void InitializeElementalVariables(ElementalVariables & rElementalVariables);

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
	  buffer << "TwoStepUpdatedLagrangianVPFluidElement #" << BaseType::Id();
	  return buffer.str();
	}

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const
      {
	rOStream << "TwoStepUpdatedLagrangianVPFluidElement" << TDim << "D";
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

      void ComputeMaterialParameters (double& Density,
				      double& DeviatoricCoeff,
				      double& VolumetricCoeff,
				      double timeStep,
				      ElementalVariables& rElementalVariables);

   

      /**
       * A constistent mass matrix is used.
       * @param rMassMatrix The local matrix where the result will be added.
       * @param rN Elemental shape functions.
       * @param Weight Multiplication coefficient for the matrix, typically Density times integration point weight.
       */
   
      virtual void ComputeMeanValueMaterialTangentMatrix(ElementalVariables& rElementalVariables,
						 double& MeanValue,
						 const ShapeFunctionDerivativesType& rShapeDeriv,
						 const double secondLame,
						 double& bulkModulus,
						 const double Weight,
						 double& MeanValueMass,
						 const double TimeStep);   
	
      virtual void ComputeBulkReductionCoefficient(MatrixType MassMatrix,
						   MatrixType StiffnessMatrix,
						   double& meanValueStiff,
						   double& bulkCoefficient,
						   double timeStep);
      
      double ComputeNonLinearViscosity(double & equivalentStrainRate);

     void ComputeBulkMatrixForPressureVelLump(MatrixType& BulkVelMatrix,
					      const double Weight);
	
     void ComputeBulkMatrixForPressureAccLump(MatrixType& BulkAccMatrix,
					      const double Weight);

     void ComputeBulkMatrixForPressureVel(MatrixType& BulkVelMatrix,
					  const ShapeFunctionsType& rN,
					  const double Weight);
	
     void ComputeBulkMatrixForPressureAcc(MatrixType& BulkAccMatrix,
					  const ShapeFunctionsType& rN,
					  const double Weight);

     void ComputeBoundLHSMatrix(MatrixType& BoundLHSMatrix,
				const ShapeFunctionsType& rN,
				const double Weight);

     void ComputeBoundRHSVector(VectorType& BoundRHSVector,
				const ShapeFunctionsType& rN,
				const double TimeStep,
				const double BoundRHSCoeffAcc,
				const double BoundRHSCoeffDev);

     void ComputeStabLaplacianMatrix(MatrixType& StabLaplacianMatrix,
				     const ShapeFunctionDerivativesType& rShapeDeriv,
				     const double Weight);


      /* bool CalcMechanicsUpdated(ElementalVariables & rElementalVariables, */
      /* 				const ProcessInfo& rCurrentProcessInfo, */
      /* 				const ShapeFunctionDerivativesType& rDN_DX, */
      /* 				unsigned int g); */

      void GetPositions(Vector& rValues,
			const ProcessInfo& rCurrentProcessInfo,
			const double theta);
	
      void CalcElasticPlasticCauchySplitted(ElementalVariables & rElementalVariables,
					    double TimeStep,
					    unsigned int g);

      virtual void CalculateTauFIC(double& TauOne,
				   double ElemSize,
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
					 const double Weight,
					 const ShapeFunctionDerivativesType& rDN_DX,
					 const SizeType i);

      void CalculateLocalContinuityEqForPressure(MatrixType& rLeftHandSideMatrix,
						 VectorType& rRightHandSideVector,
						 ProcessInfo& rCurrentProcessInfo); 

      void GetPressureVelocityValues(Vector& rValues,
				     const int Step);


      void GetPressureAccelerationValues(Vector& rValues,
					 const int Step);

      double GetThetaMomentum (){return 0.5;};

      double GetThetaContinuity (){return 1.0;};
      


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
      TwoStepUpdatedLagrangianVPFluidElement & operator=(TwoStepUpdatedLagrangianVPFluidElement const& rOther);

      /* /// Copy constructor. */
      /* TwoStepUpdatedLagrangianVPFluidElement(TwoStepUpdatedLagrangianVPFluidElement const& rOther); */

      ///@}

    }; // Class TwoStepUpdatedLagrangianVPFluidElement

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  template< unsigned int TDim >
    inline std::istream& operator >>(std::istream& rIStream,
                                     TwoStepUpdatedLagrangianVPFluidElement<TDim>& rThis)
    {
      return rIStream;
    }

  /// output stream function
  template< unsigned int TDim >
    inline std::ostream& operator <<(std::ostream& rOStream,
                                     const TwoStepUpdatedLagrangianVPFluidElement<TDim>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }

} // namespace Kratos.

#endif // KRATOS_TWO_STEP_UPDATED_LAGRANGIAN_V_P_FLUID_ELEMENT  defined
