//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:              LMonforte $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                 -0.1 $
//
//

#if !defined(KRATOS_UPDATED_LAGRANGIAN_U_J_wP_ELEMENT_H_INCLUDED )
#define  KRATOS_UPDATED_LAGRANGIAN_U_J_wP_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/updated_lagrangian_U_J_element.hpp"


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


   class UpdatedLagrangianUJwPElement
      : public UpdatedLagrangianUJElement
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
         KRATOS_CLASS_POINTER_DEFINITION( UpdatedLagrangianUJwPElement );
         ///@}

         ///@name Life Cycle
         ///@{

         /// Empty constructor needed for serialization
         UpdatedLagrangianUJwPElement();

         /// Default constructors
         UpdatedLagrangianUJwPElement(IndexType NewId, GeometryType::Pointer pGeometry);

         UpdatedLagrangianUJwPElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

         ///Copy constructor
         UpdatedLagrangianUJwPElement(UpdatedLagrangianUJwPElement const& rOther);


         /// Destructor.
         virtual ~UpdatedLagrangianUJwPElement();

         ///@}
         ///@name Operators
         ///@{

         /// Assignment operator.
         UpdatedLagrangianUJwPElement& operator=(UpdatedLagrangianUJwPElement const& rOther);


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

         /**
          * Get on rVariable a double Value from the Element Constitutive Law
          */

         void CalculateOnIntegrationPoints(const Variable<Vector>& rVariable, std::vector<Vector>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

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


         //************* COMPUTING  METHODS

         /**
          * this is called during the assembling process in order
          * to calculate the elemental mass matrix
          * @param rMassMatrix: the elemental mass matrix
          * @param rCurrentProcessInfo: the current process info instance
          */
         void CalculateMassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo) override;

         /**
          * this is called during the assembling process in order
          * to calculate the elemental damping matrix
          * @param rDampingMatrix: the elemental damping matrix
          * @param rCurrentProcessInfo: the current process info instance
          */
         void CalculateDampingMatrix(MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo) override;


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
          * Initialize System Matrices
          */
         void InitializeSystemMatrices(MatrixType& rLeftHandSideMatrix,
               VectorType& rRightHandSideVector,
               Flags& rCalculationFlags) override;
    
    
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


}; // Class UpdatedLagrangianUJwPElement



} // namespace Kratos
#endif // KRATOS_UPDATED_LAGRANGIAN_U_J_wP_ELEMENT_H_INCLUDED

