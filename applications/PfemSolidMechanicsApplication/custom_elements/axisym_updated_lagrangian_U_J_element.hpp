//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:              LMonforte $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_AXISYM_UPDATED_LAGRANGIAN_U_J_ELEMENT_H_INCLUDED )
#define  KRATOS_AXISYM_UPDATED_LAGRANGIAN_U_J_ELEMENT_H_INCLUDED

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



   class AxisymUpdatedLagrangianUJElement
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
         KRATOS_CLASS_POINTER_DEFINITION( AxisymUpdatedLagrangianUJElement );
         ///@}

         ///@name Life Cycle
         ///@{

         /// Empty constructor needed for serialization
         AxisymUpdatedLagrangianUJElement();

         /// Default constructors
         AxisymUpdatedLagrangianUJElement(IndexType NewId, GeometryType::Pointer pGeometry);

         AxisymUpdatedLagrangianUJElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

         ///Copy constructor
         AxisymUpdatedLagrangianUJElement(AxisymUpdatedLagrangianUJElement const& rOther);


         /// Destructor.
         virtual ~AxisymUpdatedLagrangianUJElement();

         ///@}
         ///@name Operators
         ///@{

         /// Assignment operator.
         AxisymUpdatedLagrangianUJElement& operator=(AxisymUpdatedLagrangianUJElement const& rOther);


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
         Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override ;

         /**
          * clones the selected element variables, creating a new one
          * @param NewId: the ID of the new element
          * @param ThisNodes: the nodes of the new element
          * @param pProperties: the properties assigned to the new element
          * @return a Pointer to the new element
          */
         Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const  override;

         //************* GETTING METHODS


         //************* STARTING - ENDING  METHODS

         /**
          * Called to initialize the element.
          * Must be called before any calculation is done
          */
         void Initialize() override;

         //************* COMPUTING  METHODS

         /**
          * this is called during the assembling process in order
          * to calculate the elemental mass matrix
          * @param rMassMatrix: the elemental mass matrix
          * @param rCurrentProcessInfo: the current process info instance
          */
         void CalculateMassMatrix(MatrixType& rMassMatrix, 
               ProcessInfo& rCurrentProcessInfo) override;


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

         /**
          * Calculate Element Kinematics
          */
         virtual void CalculateKinematics(ElementDataType& rVariables,
               const double& rPointNumber) override;


         /**
          * Calculation of the Deformation Gradient F
          */
         void CalculateDeformationGradient(const Matrix& rDN_DX,
               Matrix& rF,
               Matrix& rDeltaPosition, 
               double & rCurrentRadius, 
               double & rReferenceRadius);

         /**
          * Calculate Radius in the current and deformed geometry
          */
         virtual void CalculateRadius(double & rCurrentRadius, double & rReferenceRadius, const Vector& rN);


         /**
          * Calculation of the Deformation Matrix  BL
          */
         void CalculateDeformationMatrix(Matrix& rB,
                                            Matrix& rF,
                                            Vector& rN,
                                            double& rCurrentRadius);

         
         /**
          * Calculation of the Green Lagrange Strain Vector
          */
         void CalculateGreenLagrangeStrain(const Matrix& rF,
               Vector& rStrainVector) override;

         /**
          * Calculation of the Almansi Strain Vector
          */
         void CalculateAlmansiStrain(const Matrix& rF,
               Vector& rStrainVector) override;

         

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
               ElementDataType& rVariables,
               double& rIntegrationWeight) override;

         /**
          * Calculation and addition of the vectors of the RHS
          */

         virtual void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
               ElementDataType& rVariables,
               Vector& rVolumeForce,
               double& rIntegrationWeight) override;

         /**
          * Initialize Element General Variables
          */
         virtual void InitializeElementData(ElementDataType & rVariables, const ProcessInfo& rCurrentProcessInfo) override;




         /**
          * Calculation of the Material Stiffness Matrix. Kuum = BT * D * B
          */
         virtual void CalculateAndAddKuum(MatrixType& rK,
               ElementDataType & rVariables,
               double& rIntegrationWeight
               ) override;

         /**
          * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
          */
         virtual void CalculateAndAddKuug(MatrixType& rK,
               ElementDataType & rVariables,
               double& rIntegrationWeight
               ) override;

         /**
          * Calculation of the Kup matrix
          */
         virtual void CalculateAndAddKuJ (MatrixType& rK,
               ElementDataType & rVariables,
               double& rIntegrationWeight
               ) override;

         /**
          * Calculation of the Kpu matrix
          */
         virtual void CalculateAndAddKJu(MatrixType& rK,
               ElementDataType & rVariables,
               double& rIntegrationWeight
               ) override;



         /**
          * Calculation of the Kpp matrix
          */
         virtual void CalculateAndAddKJJ(MatrixType& rK,
               ElementDataType & rVariables,
               double& rIntegrationWeight
               ) override;

         /**
          * Calculation of the Kpp Stabilization Term matrix
          */
         virtual void CalculateAndAddKJJStab(MatrixType& rK, ElementDataType & rVariables,
               double& rIntegrationWeight
               ) override;


         /**
          * Calculation of the Internal Forces due to Pressure-Balance
          */
         virtual void CalculateAndAddJacobianForces(VectorType& rRightHandSideVector,
               ElementDataType & rVariables,
               double& rIntegrationWeight
               ) override;

         /**
          * Calculation of the Internal Forces due to Pressure-Balance
          */
         virtual void CalculateAndAddStabilizedJacobian(VectorType& rRightHandSideVector,
               ElementDataType & rVariables,
               double& rIntegrationWeight
               ) override;




         /**
          * Get the Historical Deformation Gradient to calculate after finalize the step
          */
         virtual void GetHistoricalVariables( ElementDataType& rVariables, 
               const double& rPointNumber ) override;


         /*
          * Function to modify the deformation gradient to the constitutitve equation
          */
         virtual void ComputeConstitutiveVariables( ElementDataType& rVariables, Matrix& rFT, double& rdetFT) override; 

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

         virtual void save(Serializer& rSerializer) const override;

         virtual void load(Serializer& rSerializer) override;


         ///@name Private Inquiry
         ///@{
         ///@}
         ///@name Un accessible methods
         ///@{
         ///@}


}; // Class AxisymUpdatedLagrangianUJElement



} // namespace Kratos
#endif // KRATOS_UPDATED_LAGRANGIAN_U_J_ELEMENT_H_INCLUDED

