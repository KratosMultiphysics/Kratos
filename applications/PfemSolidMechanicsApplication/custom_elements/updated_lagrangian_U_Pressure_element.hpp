//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_UPDATED_LAGRANGIAN_UP_SECOND_ELEMENT_H_INCLUDED )
#define  KRATOS_UPDATED_LAGRANGIAN_UP_SECOND_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/solid_elements/updated_lagrangian_U_P_element.hpp"

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


   class KRATOS_API(PFEM_SOLID_MECHANICS_APPLICATION) UpdatedLagrangianUPressureElement
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
            
         } ThisElementData;
      
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
         ///Type for size
         typedef GeometryData::SizeType SizeType;
         ///Type for element variables
         typedef UpdatedLagrangianUPElement::ElementDataType ElementDataType;

         /// Counted pointer of LargeDisplacementUPElement
         KRATOS_CLASS_POINTER_DEFINITION( UpdatedLagrangianUPressureElement );
         ///@}

         ///@name Life Cycle
         ///@{

         /// Empty constructor needed for serialization
         UpdatedLagrangianUPressureElement();

         /// Default constructors
         UpdatedLagrangianUPressureElement(IndexType NewId, GeometryType::Pointer pGeometry);

         UpdatedLagrangianUPressureElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

         ///Copy constructor
         UpdatedLagrangianUPressureElement(UpdatedLagrangianUPressureElement const& rOther);


         /// Destructor.
         virtual ~UpdatedLagrangianUPressureElement();

         ///@}
         ///@name Operators
         ///@{

         /// Assignment operator.
         UpdatedLagrangianUPressureElement& operator=(UpdatedLagrangianUPressureElement const& rOther);


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

         /*
          * Get on rVariable a Matrix Value from the Element Constitutive Law
          */
         void GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable, std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo) override;
         
         
         // Get Values defined in order to avoid a clang warning (?)
         void GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) override;

         /**
          * Calculate a Matrix Variable on the Element Constitutive Law
          */
         virtual void CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& rOutput, const ProcessInfo& rCurrentProcessInfo) override;


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
         virtual void InitializeElementData(ElementDataType & rVariables, 
               const ProcessInfo& rCurrentProcessInfo) override;



         /**
          * Calculation of the Material Stiffness Matrix. Kuum = BT * D * B
          */
         virtual void CalculateAndAddKuumElemUP(MatrixType& rK,
               ElementDataType & rVariables,
               ThisElementData & rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Geometric Stiffness Matrix. Kuug = BT * S
          */
         virtual void CalculateAndAddKuugElemUP(MatrixType& rK,
               ElementDataType & rVariables,
               ThisElementData & rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Kup matrix
          */
         virtual void CalculateAndAddKupElemUP (MatrixType& rK,
               ElementDataType & rVariables,
               ThisElementData & rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Kpu matrix
          */
         virtual void CalculateAndAddKpuElemUP(MatrixType& rK,
               ElementDataType & rVariables,
               ThisElementData & rElementVariables,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the Kpp matrix
          */
         virtual void CalculateAndAddKppElemUP(MatrixType& rK,
               ElementDataType & rVariables,
               ThisElementData & rElementVariables,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the Kpp Stabilization Term matrix
          */
         virtual void CalculateAndAddKppStabElemUP(MatrixType& rK,
               ElementDataType & rVariables,
               ThisElementData & rElementVariables,
               double& rIntegrationWeight
               );

         /*
          * Calculation of the Internal Forces due to sigma. Fi = B * sigma
          */
         void CalculateAndAddInternalForcesElemUP(VectorType& rRightHandSideVector,
               ElementDataType & rVariables,
               ThisElementData & rElementVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Internal Forces due to Pressure-Balance
          */
         virtual void CalculateAndAddPressureForcesElemUP(VectorType& rRightHandSideVector,
               ElementDataType & rVariables,
               ThisElementData & rElementVariables,
               double& rIntegrationWeight
               );


         /**
          * Calculation of the Internal Forces due to Pressure-Balance
          */
         virtual void CalculateAndAddStabilizedPressureElemUP(VectorType& rRightHandSideVector,
               ElementDataType & rVariables,
               ThisElementData & rElementVariables,
               double& rIntegrationWeight
               );

         /**
           * Calculate ThisElementData ( variables that are easy computations that are performed several times 
           */
         virtual void CalculateThisElementData( ThisElementData& rElementVariables, const ElementDataType & rVariables);


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


}; // Class UpdatedLagrangianUPressureElement



} // namespace Kratos
#endif // KRATOS_____

