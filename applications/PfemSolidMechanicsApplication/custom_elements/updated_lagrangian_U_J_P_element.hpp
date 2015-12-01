//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_UPDATED_LAGRANGIAN_U_J_P_ELEMENT_H_INCLUDED )
#define      KRATOS_UPDATED_LAGRANGIAN_U_J_P_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "solid_mechanics_application.h"
#include "custom_elements/large_displacement_element.hpp"


namespace Kratos
{
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

   /// Updated Lagrangian Large Displacement Lagrangian U-J-P Element for 3D and 2D geometries. Linear Triangles and Tetrahedra (base class)
   /// This is the element from Chapter 5, Vol 2. of Z


   class UpdatedLagrangianUJPElement
      : public LargeDisplacementElement
   {

      protected:
         typedef struct
         {
            unsigned int voigtsize;

            double Alpha;
            double Beta; 

            double NodalMeanStress;
            double ElementalMeanStress;

            double NodalJacobian;

            Vector StressVector;
            Vector StressVectorEC;

         } UJPGeneralVariables;
      
      public:

         ///@name Type Definitions
         ///@{
         ///Reference type definition for constitutive laws
         typedef ConstitutiveLaw ConstitutiveLawType;
         ///Pointer type for constitutive laws
         typedef ConstitutiveLawType::Pointer ConstitutiveLawPointerType;
         ///StressMeasure from constitutive laws
         typedef ConstitutiveLawType::StressMeasure StressMeasureType;
         ///Type definition for integration methods
         typedef GeometryData::IntegrationMethod IntegrationMethod;

         /// Counted pointer of LargeDisplacementUPElement
         KRATOS_CLASS_POINTER_DEFINITION( UpdatedLagrangianUJPElement );
         ///@}

         ///@name Life Cycle
         ///@{

         /// Empty constructor needed for serialization
         UpdatedLagrangianUJPElement();

         /// Default constructors
         UpdatedLagrangianUJPElement(IndexType NewId, GeometryType::Pointer pGeometry);

         UpdatedLagrangianUJPElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

         ///Copy constructor
         UpdatedLagrangianUJPElement(UpdatedLagrangianUJPElement const& rOther);


         /// Destructor.
         virtual ~UpdatedLagrangianUJPElement();

         ///@}
         ///@name Operators
         ///@{

         /// Assignment operator.
         UpdatedLagrangianUJPElement& operator=(UpdatedLagrangianUJPElement const& rOther);


         ///@}
         ///@name Operations
         ///@{
         /**
          * Returns the currently selected integration method
          * @return current integration method selected
          */
         /**
          * creates a new total lagrangian updated element pointer
          * @param NewId: the ID of the new element
          * @param ThisNodes: the nodes of the new element
          * @param pProperties: the properties assigned to the new element
          * @return a Pointer to the new element
          */
         Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const;

         /**
          * clones the selected element variables, creating a new one
          * @param NewId: the ID of the new element
          * @param ThisNodes: the nodes of the new element
          * @param pProperties: the properties assigned to the new element
          * @return a Pointer to the new element
          */
         Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const;

         //************* GETTING METHODS

         //SET

         /**
          * Set a double  Value on the Element Constitutive Law
          */
         void SetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);


         //GET:

         /**
          * Get on rVariable a double Value from the Element Constitutive Law
          */
         void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);

         void GetValueOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo);

         void GetValueOnIntegrationPoints( const Variable<Matrix>& rVariable, std::vector<Matrix>& rValue, const ProcessInfo& rCurrentProcessInfo);

         //************* STARTING - ENDING  METHODS

         /**
          * Called to initialize the element.
          * Must be called before any calculation is done
          */
         void Initialize();

         /**
          * Sets on rElementalDofList the degrees of freedom of the considered element geometry
          */
         void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo);

         /**
          * Sets on rResult the ID's of the element degrees of freedom
          */
         void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo);

         /**
          * Sets on rValues the nodal displacements
          */
         void GetValuesVector(Vector& rValues, int Step = 0);

         /**
          * Sets on rValues the nodal velocities
          */
         void GetFirstDerivativesVector(Vector& rValues, int Step = 0);

         /**
          * Sets on rValues the nodal accelerations
          */
         void GetSecondDerivativesVector(Vector& rValues, int Step = 0);

         /**
          * Called at the end of eahc solution step
          */
         void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo);

         //************************************************************************************
         //************************************************************************************
         /**
          * This function provides the place to perform checks on the completeness of the input.
          * It is designed to be called only once (or anyway, not often) typically at the beginning
          * of the calculations, so to verify that nothing is missing from the input
          * or that no common error is found.
          * @param rCurrentProcessInfo
          */
         int Check(const ProcessInfo& rCurrentProcessInfo);

         /**
          * Calculate Element Kinematics
          */
         virtual void CalculateKinematics(GeneralVariables& rVariables,
               const double& rPointNumber);


         /**
          * Calculation of the Deformation Gradient F
          */
         void CalculateDeformationGradient(const Matrix& rDN_DX,
               Matrix& rF,
               Matrix& rDeltaPosition);

         /**
          * Calculation of the Deformation Matrix  BL
          */
         virtual void CalculateDeformationMatrix(Matrix& rB,
               Matrix& rF,
               Matrix& rDN_DX);

         ///@}
         ///@name Access
         ///@{

         ///@}
         ///@name Inquiry
         ///@{
         ///@}
         ///@name Input and output
         ///@{
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

         /**
           * Variables to scale some equations (maybe)
           */
         double ScalingConstant1;
         double ScalingConstant2;

         /**
          * Container for historical total elastic deformation measure F0 = dx/dX
          */
         std::vector< Matrix > mDeformationGradientF0;

         /**
          * Container for the total deformation gradient determinants
          */
         Vector mDeterminantF0;

         /**** 
           the time step (requiered). It shall be somewhere else.
          ****/    
         double mTimeStep;

         /*** 
           Just to check a few things
          ***/
         //bool mCompressibleWater;

         ///@}
         ///@name Protected Operators
         ///@{

         /**
          * Calculates the elemental contributions
          * \f$ K^e = w\,B^T\,D\,B \f$ and
          * \f$ r^e \f$
          */
         virtual void CalculateElementalSystem(LocalSystemComponents& rLocalSystem,
               ProcessInfo& rCurrentProcessInfo);

         ///@}
         ///@name Protected Operations
         ///@{

         /**
          * Calculation and addition of the matrices of the LHS
          */

         virtual void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
               GeneralVariables& rVariables,
               double& rIntegrationWeight);

         /**
          * Calculation and addition of the vectors of the RHS
          */

         virtual void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
               GeneralVariables& rVariables,
               Vector& rVolumeForce,
               double& rIntegrationWeight);

         /**
          * Initialize Element General Variables
          */
         virtual void InitializeGeneralVariables(GeneralVariables & rVariables, const ProcessInfo& rCurrentProcessInfo);

         /**
          * Finalize Element Internal Variables
          */
         virtual void FinalizeStepVariables(GeneralVariables & rVariables, const double& rPointNumber);



         /**
          * Calculation of the Material Stiffness Matrix. Kuum = BT * D * B
          */
         virtual void CalculateAndAddKuum(MatrixType& rK,
               GeneralVariables & rVariables,
               UJPGeneralVariables& rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
          */
         virtual void CalculateAndAddKuug(MatrixType& rK,
               GeneralVariables & rVariables,
               UJPGeneralVariables& rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the KuJ matrix
          */
         virtual void CalculateAndAddKuJ (MatrixType& rK,
               GeneralVariables & rVariables,
               UJPGeneralVariables& rElementVariables,
               double& rIntegrationWeight
               );
         /**
          * Calculation of the Kup matrix
          */
         virtual void CalculateAndAddKup (MatrixType& rK,
               GeneralVariables & rVariables,
               UJPGeneralVariables& rElementVariables,
               double& rIntegrationWeight
               );
         /**
          * Calculation of the KJu matrix
          */
         virtual void CalculateAndAddKJu(MatrixType& rK,
               GeneralVariables & rVariables,
               UJPGeneralVariables& rElementVariables,
               double& rIntegrationWeight
               );
         /**
          * Calculation of the KJJ matrix
          */
         virtual void CalculateAndAddKJJ(MatrixType& rK,
               GeneralVariables & rVariables,
               UJPGeneralVariables& rElementVariables,
               double& rIntegrationWeight
               );
         /**
          * Calculation of the KJp matrix
          */
         virtual void CalculateAndAddKJp(MatrixType& rK,
               GeneralVariables & rVariables,
               UJPGeneralVariables& rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Kpu matrix
          */
         virtual void CalculateAndAddKpu(MatrixType& rK,
               GeneralVariables & rVariables,
               UJPGeneralVariables& rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the KpJ matrix
          */
         virtual void CalculateAndAddKpJ(MatrixType& rK,
               GeneralVariables & rVariables,
               UJPGeneralVariables& rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Kpp matrix
          */
         virtual void CalculateAndAddKpp(MatrixType& rK,
               GeneralVariables & rVariables,
               UJPGeneralVariables& rElementVariables,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the Kpp Stabilization Term matrix
          */
         virtual void CalculateAndAddKppStab(MatrixType& rK,
               GeneralVariables & rVariables,
               UJPGeneralVariables& rElementVariables,
               double& rIntegrationWeight
               );
         /**
           * Calculation of the KJJ Stabilization Term matrix
           */
         virtual void CalculateAndAddKJJStab(MatrixType& rK,
               GeneralVariables & rVariables,
               UJPGeneralVariables& rElementVariables,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the External Forces Vector. Fe = N * t + N * b
          */
         void CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
               GeneralVariables& rVariables,
               Vector& rVolumeForce,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Internal Forces due to Jacobian-Balance
          */
         virtual void CalculateAndAddJacobianForces( VectorType& rRightHandSideVector,
               GeneralVariables & rVariables,
               UJPGeneralVariables& rElementVariables,
               double& rIntegrationWeight
               );



         /**
          * Calculation of the Internal Forces due to Pressure-Balance
          */
         virtual void CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
               GeneralVariables & rVariables,
               UJPGeneralVariables& rElementVariables,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the Stabilization for the Pressure
          */
         virtual void CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
               GeneralVariables & rVariables,
               UJPGeneralVariables& rElementVariables,
               double& rIntegrationWeight
               );


         /**
           * Calculation of the Stabilization for the Jacobian
           */
         virtual void CalculateAndAddStabilizedJacobian( VectorType& rRightHandSideVector,
               GeneralVariables & rVariables,
               UJPGeneralVariables& rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Internal Forces due to sigma. Fi = B * sigma
          */
         virtual void CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
               GeneralVariables & rVariables,
               UJPGeneralVariables& rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Initialize System Matrices
          */
         void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
               VectorType& rRightHandSideVector,
               Flags& rCalculationFlags);

         //on integration points:
         /**
          * Calculate a double Variable on the Element Constitutive Law
          */
         void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rOutput, const ProcessInfo& rCurrentProcessInfo);

         void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rOutput, const ProcessInfo& rCurrentProcessInfo);

         void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rOutput, const ProcessInfo& rCurrentProcessInfo);


         virtual double GetElementSize( const Matrix& rDN_DX);

         /**
          * Get the Historical Deformation Gradient to calculate after finalize the step
          */
         virtual void GetHistoricalVariables( GeneralVariables& rVariables, 
               const double& rPointNumber );

         /**
          * Calculation of the Volume Change of the Element
          */
         virtual double& CalculateVolumeChange(double& rVolumeChange, GeneralVariables& rVariables);


         virtual void CalculateThisElementGeneralVariables( UJPGeneralVariables& rElementVariables, const GeneralVariables& rVariables);
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
         ///@name Private Operators
         ///@{


         ///@}
         ///@name Private Operations
         ///@{


         ///@}
         ///@name Private  Access
         ///@{
         ///@}

         ///@}
         ///@name Serialization
         ///@{
         friend class Serializer;

         // A private default constructor necessary for serialization

         virtual void save(Serializer& rSerializer) const;

         virtual void load(Serializer& rSerializer);


         ///@name Private Inquiry
         ///@{
         ///@}
         ///@name Un accessible methods
         ///@{
         ///@}


}; // Class UpdatedLagrangianUJPElement



} // namespace Kratos
#endif // KRATOS_UPDATED_LAGRANGIAN_U_J_P_ELEMENT_H_INCLUDED

