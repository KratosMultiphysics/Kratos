//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:              LMonforte $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                 -0.1 $
//
//

#if !defined(KRATOS_UPDATED_LAGRANGIAN_UP_SECOND_ELEMENT_H_INCLUDED )
#define  KRATOS_UPDATED_LAGRANGIAN_UP_SECOND_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/updated_lagrangian_U_P_element.hpp"

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

   /// Large Displacement Lagrangian U-P Element for 3D and 2D geometries. Linear Triangles and Tetrahedra (base class)

   // THIS IS THE SECOND VERSION OF THE UP-ELEMENT FOR ANY CONSTITUTIVE EQUATION. 
   // In this version, the pressure is the pressure and it is not some "2D" invented pressure,..
   // the constitutive equation should not be the UP version


   class UpdatedLagrangianUPSecondElement
      : public UpdatedLagrangianUPElement
   {

      protected:
         typedef struct
         {
            unsigned int voigtsize; 

            double NodalMeanStress;
            double ElementalMeanStress;

            Vector StressVector;

            Matrix DeviatoricTensor; 
            
         } ThisElementGeneralVariables;
      
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
         KRATOS_CLASS_POINTER_DEFINITION( UpdatedLagrangianUPSecondElement );
         ///@}

         ///@name Life Cycle
         ///@{

         /// Empty constructor needed for serialization
         UpdatedLagrangianUPSecondElement();

         /// Default constructors
         UpdatedLagrangianUPSecondElement(IndexType NewId, GeometryType::Pointer pGeometry);

         UpdatedLagrangianUPSecondElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

         ///Copy constructor
         UpdatedLagrangianUPSecondElement(UpdatedLagrangianUPSecondElement const& rOther);


         /// Destructor.
         virtual ~UpdatedLagrangianUPSecondElement();

         ///@}
         ///@name Operators
         ///@{

         /// Assignment operator.
         UpdatedLagrangianUPSecondElement& operator=(UpdatedLagrangianUPSecondElement const& rOther);


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

         /*
          * Get on rVariable a Matrix Value from the Element Constitutive Law
          */
         virtual void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo);

         /**
          * Calculate a Matrix Variable on the Element Constitutive Law
          */
         void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo);


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

         double mElementStabilizationNumber;

         double mElementScalingNumber; 

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
          * Set Variables of the Element to the Parameters of the Constitutive Law
          */
         virtual void SetGeneralVariables(GeneralVariables& rVariables,
               ConstitutiveLaw::Parameters& rValues,
               const int & rPointNumber);

         /**
          * Initialize Element General Variables
          */
         virtual void InitializeGeneralVariables(GeneralVariables & rVariables, 
               const ProcessInfo& rCurrentProcessInfo);



         virtual double GetElementSize( const Matrix& rDN_DX);
         /**
          * Set Variables of the Element to the Parameters of the Constitutive Law
          */
         //virtual void SetGeneralVariables(GeneralVariables& rVariables,
         //                                 ConstitutiveLaw::Parameters& rValues,
         //                                 const int & rPointNumber);


         /**
          * Calculation of the Material Stiffness Matrix. Kuum = BT * D * B
          */
         virtual void CalculateAndAddKuum(MatrixType& rK,
               GeneralVariables & rVariables,
               ThisElementGeneralVariables & rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
          */
         virtual void CalculateAndAddKuug(MatrixType& rK,
               GeneralVariables & rVariables,
               ThisElementGeneralVariables & rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Kup matrix
          */
         virtual void CalculateAndAddKup (MatrixType& rK,
               GeneralVariables & rVariables,
               ThisElementGeneralVariables & rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Kpu matrix
          */
         virtual void CalculateAndAddKpu(MatrixType& rK,
               GeneralVariables & rVariables,
               ThisElementGeneralVariables & rElementVariables,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the Kpp matrix
          */
         virtual void CalculateAndAddKpp(MatrixType& rK,
               GeneralVariables & rVariables,
               ThisElementGeneralVariables & rElementVariables,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the Kpp Stabilization Term matrix
          */
         virtual void CalculateAndAddKppStab(MatrixType& rK,
               GeneralVariables & rVariables,
               ThisElementGeneralVariables & rElementVariables,
               double& rIntegrationWeight
               );
         /**
          * Calculation of the External Forces Vector. Fe = N * t + N * b
          */
         //void CalculateAndAddExternalForces(VectorType& rRightHandSideVector,
         //                                   GeneralVariables& rVariables,
         //                                   Vector& rVolumeForce,
         //                                   double& rIntegrationWeight
         //                                  );

         /*
          * Calculation of the Internal Forces due to sigma. Fi = B * sigma
          */
         void CalculateAndAddInternalForces(VectorType& rRightHandSideVector,
               GeneralVariables & rVariables,
               ThisElementGeneralVariables & rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Internal Forces due to Pressure-Balance
          */
         virtual void CalculateAndAddPressureForces(VectorType& rRightHandSideVector,
               GeneralVariables & rVariables,
               ThisElementGeneralVariables & rElementVariables,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the Internal Forces due to Pressure-Balance
          */
         virtual void CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
               GeneralVariables & rVariables,
               ThisElementGeneralVariables & rElementVariables,
               double& rIntegrationWeight
               );

         /**
           * Calculate ThisElementGeneralVariables ( variables that are easy computations that are performed several times 
           */
         virtual void CalculateThisElementGeneralVariables( ThisElementGeneralVariables& rElementGeneralVariables, const GeneralVariables & rVariables);


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


}; // Class UpdatedLagrangianUPSecondElement



} // namespace Kratos
#endif // KRATOS_____

