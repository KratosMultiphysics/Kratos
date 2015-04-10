//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_AXISYM_SPATIAL_LAGRANGIAN_U_wP_STAB_ELEMENT_H_INCLUDED )
#define  KRATOS_AXISYM_SPATIAL_LAGRANGIAN_U_wP_STAB_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/axisym_spatial_lagrangian_U_wP_element.hpp"

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

   /// Stabilization of the Axisim Spatial Lagrangian Large Displacement Lagrangian U-wP Element for 3D and 2D geometries. 

   /**
    * Implements a Large Displacement Lagrangian definition for structural analysis.
    * This works for arbitrary geometries in 3D and 2D (base class)
    */

   class AxisymSpatialLagrangianUwPStabElement
      : public AxisymSpatialLagrangianUwPElement
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
         KRATOS_CLASS_POINTER_DEFINITION( AxisymSpatialLagrangianUwPStabElement );
         ///@}

         ///@name Life Cycle
         ///@{

         /// Empty constructor needed for serialization
         AxisymSpatialLagrangianUwPStabElement();

         /// Default constructors
         AxisymSpatialLagrangianUwPStabElement(IndexType NewId, GeometryType::Pointer pGeometry);

         AxisymSpatialLagrangianUwPStabElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

         ///Copy constructor
         AxisymSpatialLagrangianUwPStabElement(AxisymSpatialLagrangianUwPStabElement const& rOther);


         /// Destructor.
         virtual ~AxisymSpatialLagrangianUwPStabElement();

         ///@}
         ///@name Operators
         ///@{

         /// Assignment operator.
         AxisymSpatialLagrangianUwPStabElement& operator=(AxisymSpatialLagrangianUwPStabElement const& rOther);


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


         /**
          * Calculation of the Kpp Stabilization Term matrix
          */
         virtual void CalculateAndAddKppStab(MatrixType& rK,
               GeneralVariables & rVariables,
               double& rIntegrationWeight
               );

         /**
          * Calculation of the Internal Forces due to Pressure-Balance
          */
         virtual void CalculateAndAddStabilizedPressure(VectorType& rRightHandSideVector,
               GeneralVariables & rVariables,
               double& rIntegrationWeight
               );



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


}; // Class AxisymSpatialLagrangianUwPStabElement



} // namespace Kratos
#endif // KRATOS_SPATIAL_LAGRANGIAN_U_wP_STAB_ELEMENT_H_INCLUDED

