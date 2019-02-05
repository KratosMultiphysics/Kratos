//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:              LMonforte $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                 -0.1 $
//
//

#if !defined(KRATOS_UPDATED_LAGRANGIAN_U_P_wP_ELEMENT_H_INCLUDED )
#define  KRATOS_UPDATED_LAGRANGIAN_U_P_wP_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/updated_lagrangian_U_Pressure_element.hpp"


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

   /// Updated Lagrangian Large Displacement Lagrangian U-wP Element for 3D and 2D geometries. Linear Triangles and Tetrahedra (base class)


   class KRATOS_API(PFEM_SOLID_MECHANICS_APPLICATION) UpdatedLagrangianUPwPElement
      : public UpdatedLagrangianUPressureElement
   {
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
         KRATOS_CLASS_POINTER_DEFINITION( UpdatedLagrangianUPwPElement );
         ///@}

         ///@name Life Cycle
         ///@{

         /// Empty constructor needed for serialization
         UpdatedLagrangianUPwPElement();

         /// Default constructors
         UpdatedLagrangianUPwPElement(IndexType NewId, GeometryType::Pointer pGeometry);

         UpdatedLagrangianUPwPElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

         ///Copy constructor
         UpdatedLagrangianUPwPElement(UpdatedLagrangianUPwPElement const& rOther);


         /// Destructor.
         virtual ~UpdatedLagrangianUPwPElement();

         ///@}
         ///@name Operators
         ///@{

         /// Assignment operator.
         UpdatedLagrangianUPwPElement& operator=(UpdatedLagrangianUPwPElement const& rOther);


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
         Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

         /**
          * clones the selected element variables, creating a new one
          * @param NewId: the ID of the new element
          * @param ThisNodes: the nodes of the new element
          * @param pProperties: the properties assigned to the new element
          * @return a Pointer to the new element
          */
         Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

         //************* GETTING METHODS


         //************* STARTING - ENDING  METHODS

         /**
          * Sets on rElementalDofList the degrees of freedom of the considered element geometry
          */
         void GetDofList(DofsVectorType& rElementalDofList, ProcessInfo& rCurrentProcessInfo) override;

         /**
          * Sets on rResult the ID's of the element degrees of freedom
          */
         void EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo) override;

         /**
          * Sets on rValues the nodal displacements
          */
         void GetValuesVector(Vector& rValues, int Step = 0) override;

         /**
          * Sets on rValues the nodal velocities
          */
         void GetFirstDerivativesVector(Vector& rValues, int Step = 0) override;

         /**
          * Sets on rValues the nodal accelerations
          */
         void GetSecondDerivativesVector(Vector& rValues, int Step = 0) override;

         /**
          * Called at the end of eahc solution step
          */
         void FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) override;

         //************************************************************************************
         //************************************************************************************
         /**
          * This function provides the place to perform checks on the completeness of the input.
          * It is designed to be called only once (or anyway, not often) typically at the beginning
          * of the calculations, so to verify that nothing is missing from the input
          * or that no common error is found.
          * @param rCurrentProcessInfo
          */
         int Check(const ProcessInfo& rCurrentProcessInfo) override;

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

         double mTimeStep; // because I don't kwnow how to do it better 

         ///@}
         ///@name Protected Operators
         ///@{

         /**
          * Calculates the elemental contributions
          * \f$ K^e = w\,B^T\,D\,B \f$ and
          * \f$ r^e \f$
          */
         void CalculateElementalSystem(LocalSystemComponents& rLocalSystem,
               ProcessInfo& rCurrentProcessInfo) override;

         ///@}
         ///@name Protected Operations
         ///@{

         /**
          * Calculation and addition of the matrices of the LHS
          */

         void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
               ElementDataType& rVariables,
               double& rIntegrationWeight) override;

         /**
          * Calculation and addition of the vectors of the RHS
          */

         void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
               ElementDataType& rVariables,
               Vector& rVolumeForce,
               double& rIntegrationWeight) override;

         /**
          * Initialize Element General Variables
          */
         void InitializeElementData(ElementDataType & rVariables, const ProcessInfo& rCurrentProcessInfo) override;

         /**
          * Calculation of the geometric terms due to the water pressure 
          */
         virtual void CalculateAndAddUnconsideredKuuTerms(MatrixType& rK,
               ElementDataType & rVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Ku wP Matrix
          */
         virtual void CalculateAndAddKuwP(MatrixType& rK,
               ElementDataType & rVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the KwP U Matrix
          */
         virtual void CalculateAndAddKwPu(MatrixType& rK,
               ElementDataType & rVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the KwP P Matrix
          */
         virtual void CalculateAndAddKwPP(MatrixType& rK,
               ElementDataType & rVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the K wP wP Matrix
          */
         virtual void CalculateAndAddKwPwP(MatrixType& rK,
               ElementDataType & rVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Stabilization Tangent Matrix
          */
         virtual void CalculateAndAddKwPwPStab(MatrixType& rK,
               ElementDataType & rVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the External Forces Vector. Fe = N * t + N * b
          */
         void CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
               ElementDataType& rVariables,
               Vector& rVolumeForce,
               double& rIntegrationWeight
               ) override;


         /**
          * Calculation of the Internal Forces due to Pressure-Balance
          */
         void CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
               ElementDataType & rVariables,
               double& rIntegrationWeight
               ) override;


         /**
          * Calculation of the Internal Forces due to Pressure-Balance
          */
         void CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
               ElementDataType & rVariables,
               double& rIntegrationWeight
               ) override;

         /**
          * Calculation of the Internal Forces due to sigma. Fi = B * sigma
          */
         void CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
               ElementDataType & rVariables,
               double& rIntegrationWeight
               ) override;

         /**
          * Calculation of the Mass Balance ( ie water pressure equation)
          */
         virtual void CalculateAndAddWaterPressureForces( VectorType& rRightHandSideVector,
               ElementDataType& rVariables,
               double& rIntegrationWeight
               );
         /**
          * Stabilization of the MassBalance equation
          */
         virtual void CalculateAndAddStabilizedWaterPressure( VectorType& rRightHandSideVector, 
               ElementDataType& rVariables,
               double& rIntegartionWeight
               );

         /**
          * Initialize System Matrices
          */
         void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
               VectorType& rRightHandSideVector,
               Flags& rCalculationFlags) override;


         /**
          * Calculation of the Volume Change of the Element
          */
         double& CalculateVolumeChange(double& rVolumeChange, ElementDataType& rVariables) override;


         void GetConstants( double& rScalingConstant, double& rWaterBulk, double& rDeltaTime, double& rPermeability);

         virtual double GetElementSize( const Matrix& rDN_DX);

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

         void save(Serializer& rSerializer) const override;

         void load(Serializer& rSerializer) override;


         ///@name Private Inquiry
         ///@{
         ///@}
         ///@name Un accessible methods
         ///@{
         ///@}


}; // Class UpdatedLagrangianUPwPElement



} // namespace Kratos
#endif // KRATOS_UPDATED_LAGRANGIAN_U_P_wP_ELEMENT_H_INCLUDED

