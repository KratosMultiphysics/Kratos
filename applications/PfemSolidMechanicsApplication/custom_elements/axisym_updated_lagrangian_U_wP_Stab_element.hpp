//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2015 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_AXISYM_UPDATED_LAGRANGIAN_U_wP_STAB_ELEMENT_H_INCLUDED )
#define  KRATOS_AXISYM_UPDATED_LAGRANGIAN_U_wP_STAB_ELEMENT_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_elements/axisym_updated_lagrangian_U_wP_element.hpp"

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

   /// Stabilization of the Axisim Updated Lagrangian Large Displacement Lagrangian U-wP Element for 3D and 2D geometries. 

   /**
    * Implements a Large Displacement Lagrangian definition for structural analysis.
    * This works for arbitrary geometries in 3D and 2D (base class)
    */

   class AxisymUpdatedLagrangianUwPStabElement
      : public AxisymUpdatedLagrangianUwPElement
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
         KRATOS_CLASS_POINTER_DEFINITION( AxisymUpdatedLagrangianUwPStabElement );
         ///@}

         ///@name Life Cycle
         ///@{

         /// Empty constructor needed for serialization
         AxisymUpdatedLagrangianUwPStabElement();

         /// Default constructors
         AxisymUpdatedLagrangianUwPStabElement(IndexType NewId, GeometryType::Pointer pGeometry);

         AxisymUpdatedLagrangianUwPStabElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

         ///Copy constructor
         AxisymUpdatedLagrangianUwPStabElement(AxisymUpdatedLagrangianUwPStabElement const& rOther);


         /// Destructor.
         virtual ~AxisymUpdatedLagrangianUwPStabElement();

         ///@}
         ///@name Operators
         ///@{

         /// Assignment operator.
         AxisymUpdatedLagrangianUwPStabElement& operator=(AxisymUpdatedLagrangianUwPStabElement const& rOther);


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

    virtual void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
                                    ElementVariables& rVariables,
                                    double& rIntegrationWeight);

         /**
     * Calculation and addition of the vectors of the RHS
          */

    virtual void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
               ElementVariables & rVariables,
                                    Vector& rVolumeForce,
                                    double& rIntegrationWeight);




         /**
	  * Initialize Element General Variables
	  */
        virtual void InitializeElementVariables(ElementVariables & rVariables, const ProcessInfo& rCurrentProcessInfo);


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


}; // Class AxisymUpdatedLagrangianUwPStabElement



} // namespace Kratos
#endif // KRATOS_UPDATED_LAGRANGIAN_U_wP_STAB_ELEMENT_H_INCLUDED

