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


   class UpdatedLagrangianUPwPElement
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
         virtual void CalculateElementalSystem(LocalSystemComponents& rLocalSystem,
               ProcessInfo& rCurrentProcessInfo);

         ///@}
         ///@name Protected Operations
         ///@{

         /**
          * Calculation and addition of the matrices of the LHS
          */

         virtual void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
               ElementVariables& rVariables,
               double& rIntegrationWeight);

         /**
          * Calculation and addition of the vectors of the RHS
          */

         virtual void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
               ElementVariables& rVariables,
               Vector& rVolumeForce,
               double& rIntegrationWeight);

         /**
          * Initialize Element General Variables
          */
         virtual void InitializeElementVariables(ElementVariables & rVariables, const ProcessInfo& rCurrentProcessInfo);

         /**
          * Calculation of the geometric terms due to the water pressure 
          */
         virtual void CalculateAndAddUnconsideredKuuTerms(MatrixType& rK,
               ElementVariables & rVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Ku wP Matrix
          */
         virtual void CalculateAndAddKuwP(MatrixType& rK,
               ElementVariables & rVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the KwP U Matrix
          */
         virtual void CalculateAndAddKwPu(MatrixType& rK,
               ElementVariables & rVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the KwP P Matrix
          */
         virtual void CalculateAndAddKwPP(MatrixType& rK,
               ElementVariables & rVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the K wP wP Matrix
          */
         virtual void CalculateAndAddKwPwP(MatrixType& rK,
               ElementVariables & rVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Stabilization Tangent Matrix
          */
         virtual void CalculateAndAddKwPwPStab(MatrixType& rK,
               ElementVariables & rVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the External Forces Vector. Fe = N * t + N * b
          */
         void CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
               ElementVariables& rVariables,
               Vector& rVolumeForce,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the Internal Forces due to Pressure-Balance
          */
         virtual void CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
               ElementVariables & rVariables,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the Internal Forces due to Pressure-Balance
          */
         virtual void CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
               ElementVariables & rVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Internal Forces due to sigma. Fi = B * sigma
          */
         virtual void CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
               ElementVariables & rVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Mass Balance ( ie water pressure equation)
          */
         virtual void CalculateAndAddWaterPressureForces( VectorType& rRightHandSideVector,
               ElementVariables& rVariables,
               double& rIntegrationWeight
               );
         /**
          * Stabilization of the MassBalance equation
          */
         virtual void CalculateAndAddStabilizedWaterPressure( VectorType& rRightHandSideVector, 
               ElementVariables& rVariables,
               double& rIntegartionWeight
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



         /**
          * Calculation of the Volume Change of the Element
          */
         virtual double& CalculateVolumeChange(double& rVolumeChange, ElementVariables& rVariables);


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


}; // Class UpdatedLagrangianUPwPElement



} // namespace Kratos
#endif // KRATOS_UPDATED_LAGRANGIAN_U_P_wP_ELEMENT_H_INCLUDED

