//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:              LMonforte $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                 -0.1 $
//
//

#if !defined(KRATOS_AXISYM_UPDATED_LAGRANGIAN_U_J_wP_ELEMENT_H_INCLUDED )
#define  KRATOS_AXISYM_UPDATED_LAGRANGIAN_U_J_wP_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/axisym_updated_lagrangian_U_J_element.hpp"


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


   class AxisymUpdatedLagrangianUJwPElement
      : public AxisymUpdatedLagrangianUJElement
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
         KRATOS_CLASS_POINTER_DEFINITION( AxisymUpdatedLagrangianUJwPElement );
         ///@}

         ///@name Life Cycle
         ///@{

         /// Empty constructor needed for serialization
         AxisymUpdatedLagrangianUJwPElement();

         /// Default constructors
         AxisymUpdatedLagrangianUJwPElement(IndexType NewId, GeometryType::Pointer pGeometry);

         AxisymUpdatedLagrangianUJwPElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

         ///Copy constructor
         AxisymUpdatedLagrangianUJwPElement(AxisymUpdatedLagrangianUJwPElement const& rOther);


         /// Destructor.
         virtual ~AxisymUpdatedLagrangianUJwPElement();

         ///@}
         ///@name Operators
         ///@{

         /// Assignment operator.
         AxisymUpdatedLagrangianUJwPElement& operator=(AxisymUpdatedLagrangianUJwPElement const& rOther);


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

         /**
          * Get on rVariable a double Value from the Element Constitutive Law
          */

         void CalculateOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo);
         void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo);
         void CalculateOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);

         //************* STARTING - ENDING  METHODS

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
          * Calculation of the geometric terms due to the water pressure 
          */
         virtual void CalculateAndAddUnconsideredKuuTerms(MatrixType& rK,
               GeneralVariables & rVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Ku wP Matrix
          */
         virtual void CalculateAndAddKuwP(MatrixType& rK,
               GeneralVariables & rVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the KwP U Matrix
          */
         virtual void CalculateAndAddKwPu(MatrixType& rK,
               GeneralVariables & rVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the KwP J Matrix
          */
         virtual void CalculateAndAddKwPJ(MatrixType& rK,
               GeneralVariables & rVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the K wP wP Matrix
          */
         virtual void CalculateAndAddKwPwP(MatrixType& rK,
               GeneralVariables & rVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Stabilization Tangent Matrix
          */
         virtual void CalculateAndAddKwPwPStab(MatrixType& rK,
               GeneralVariables & rVariables,
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
          * Calculation of the Internal Forces due to Pressure-Balance
          */
         virtual void CalculateAndAddJacobianForces(VectorType& rRightHandSideVector,
               GeneralVariables & rVariables,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the Internal Forces due to Pressure-Balance
          */
         virtual void CalculateAndAddStabilizedJacobian(VectorType& rRightHandSideVector,
               GeneralVariables & rVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Internal Forces due to sigma. Fi = B * sigma
          */
         virtual void CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
               GeneralVariables & rVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Mass Balance ( ie water pressure equation)
          */
         virtual void CalculateAndAddWaterPressureForces( VectorType& rRightHandSideVector,
               GeneralVariables& rVariables,
               double& rIntegrationWeight
               );
         /**
          * Stabilization of the MassBalance equation
          */
         virtual void CalculateAndAddStabilizedWaterPressure( VectorType& rRightHandSideVector, 
               GeneralVariables& rVariables,
               double& rIntegartionWeight
               );

         /**
          * Initialize System Matrices
          */
         void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
               VectorType& rRightHandSideVector,
               Flags& rCalculationFlags);




         /**
          * Calculation of the Volume Change of the Element
          */
         virtual double& CalculateVolumeChange(double& rVolumeChange, GeneralVariables& rVariables);


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

         virtual void save(Serializer& rSerializer) const;

         virtual void load(Serializer& rSerializer);


         ///@name Private Inquiry
         ///@{
         ///@}
         ///@name Un accessible methods
         ///@{
         ///@}


}; // Class AxisymUpdatedLagrangianUJwPElement



} // namespace Kratos
#endif // KRATOS_UPDATED_LAGRANGIAN_U_J_wP_ELEMENT_H_INCLUDED

