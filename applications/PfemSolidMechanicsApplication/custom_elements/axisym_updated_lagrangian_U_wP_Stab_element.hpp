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

   class KRATOS_API(PFEM_SOLID_MECHANICS_APPLICATION) AxisymUpdatedLagrangianUwPStabElement
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
         Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

         /**
          * clones the selected element variables, creating a new one
          * @param NewId: the ID of the new element
          * @param ThisNodes: the nodes of the new element
          * @param pProperties: the properties assigned to the new element
          * @return a Pointer to the new element
          */
         Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

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

    void CalculateAndAddLHS(LocalSystemComponents& rLocalSystem,
                                    ElementDataType& rVariables,
                                    double& rIntegrationWeight) override;

         /**
     * Calculation and addition of the vectors of the RHS
          */

    void CalculateAndAddRHS(LocalSystemComponents& rLocalSystem,
               ElementDataType & rVariables,
                                    Vector& rVolumeForce,
                                    double& rIntegrationWeight) override;




         /**
	  * Initialize Element General Variables
	  */
        void InitializeElementData(ElementDataType & rVariables, const ProcessInfo& rCurrentProcessInfo) override;


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


}; // Class AxisymUpdatedLagrangianUwPStabElement



} // namespace Kratos
#endif // KRATOS_UPDATED_LAGRANGIAN_U_wP_STAB_ELEMENT_H_INCLUDED

