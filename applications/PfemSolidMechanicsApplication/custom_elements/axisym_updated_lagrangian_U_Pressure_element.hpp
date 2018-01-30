//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:              LMonforte $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                 -0.1 $
//
//

#if !defined(KRATOS_AXISYM_UPDATED_LAGRANGIAN_UPRESSURE_ELEMENT_H_INCLUDED )
#define  KRATOS_AXISYM_UPDATED_LAGRANGIAN_UPRESSURE_ELEMENT_H_INCLUDED

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

   // THIS IS THE SECOND VERSION OF THE UP-ELEMENT FOR ANY CONSTITUTIVE EQUATION. 
   // In this version, the pressure is the pressure and it is not some "2D" invented pressure,..
   // the constitutive equation should NOT be the UP version


   class AxisymUpdatedLagrangianUPressureElement
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
         KRATOS_CLASS_POINTER_DEFINITION( AxisymUpdatedLagrangianUPressureElement );
         ///@}

         ///@name Life Cycle
         ///@{

         /// Empty constructor needed for serialization
         AxisymUpdatedLagrangianUPressureElement();

         /// Default constructors
         AxisymUpdatedLagrangianUPressureElement(IndexType NewId, GeometryType::Pointer pGeometry);

         AxisymUpdatedLagrangianUPressureElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

         ///Copy constructor
         AxisymUpdatedLagrangianUPressureElement(AxisymUpdatedLagrangianUPressureElement const& rOther);


         /// Destructor.
         virtual ~AxisymUpdatedLagrangianUPressureElement();

         ///@}
         ///@name Operators
         ///@{

         /// Assignment operator.
         AxisymUpdatedLagrangianUPressureElement& operator=(AxisymUpdatedLagrangianUPressureElement const& rOther);


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
          * Called to initialize the element.
          * Must be called before any calculation is done
          */
         virtual void Initialize();

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
         virtual void InitializeElementVariables(ElementVariables & rVariables, 
               const ProcessInfo& rCurrentProcessInfo);


         /**
          * Calculation of the Material Stiffness Matrix. Kuum = BT * D * B
          */
         virtual void CalculateAndAddKuum(MatrixType& rK,
               ElementVariables & rVariables,
               ThisElementVariables & rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
          */
         virtual void CalculateAndAddKuug(MatrixType& rK,
               ElementVariables & rVariables,
               ThisElementVariables & rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Kup matrix
          */
         virtual void CalculateAndAddKup (MatrixType& rK,
               ElementVariables & rVariables,
               ThisElementVariables & rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Kpu matrix
          */
         virtual void CalculateAndAddKpu(MatrixType& rK,
               ElementVariables & rVariables,
               ThisElementVariables & rElementVariables,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the Kpp matrix
          */
         virtual void CalculateAndAddKpp(MatrixType& rK,
               ElementVariables & rVariables,
               ThisElementVariables & rElementVariables,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the Kpp Stabilization Term matrix
          */
         virtual void CalculateAndAddKppStab(MatrixType& rK,
               ElementVariables & rVariables,
               ThisElementVariables & rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Internal Forces due to Pressure-Balance
          */
         virtual void CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
               ElementVariables & rVariables,
               ThisElementVariables & rElementVariables,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the Internal Forces due to Pressure-Balance
          */
         virtual void CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
               ElementVariables & rVariables,
               ThisElementVariables & rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculate ThisElementVariables ( variables that are easy computations that are performed several times 
          */
         virtual void CalculateThisElementVariables( ThisElementVariables& rElementVariables, const ElementVariables & rVariables);



         /**
          * Calculate Element Kinematics
          */
         virtual void CalculateKinematics(ElementVariables& rVariables,
               const double& rPointNumber);


         /**
          * Calculate Radius in the current and deformed geometry
          */
         void CalculateRadius(double & rCurrentRadius,
               double & rReferenceRadius,
               const Vector& rN);

         /**
          * Calculation of the Deformation Gradient F
          */
         void CalculateDeformationGradient(const Matrix& rDN_DX,
               Matrix& rF,
               Matrix& rDeltaPosition,
               double & rCurrentRadius,
               double & rReferenceRadius);

         /**
          * Calculation of the Deformation Matrix  BL
          */
         void CalculateDeformationMatrix(Matrix& rB,
               Matrix& rDN_DX,
               Vector& rN,
               double & rCurrentRadius);


         /**
          * Calculation of the Green Lagrange Strain Vector
          */
         void CalculateGreenLagrangeStrain(const Matrix& rF,
               Vector& rStrainVector);

         /**
          * Calculation of the Almansi Strain Vector
          */
         void CalculateAlmansiStrain(const Matrix& rF,
               Vector& rStrainVector);



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


}; // Class AxisymUpdatedLagrangianUPressureElement



} // namespace Kratos
#endif // KRATOS_____

