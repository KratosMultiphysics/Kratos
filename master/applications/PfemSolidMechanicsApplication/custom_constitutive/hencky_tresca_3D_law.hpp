//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:            LMonforte $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2015 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined (KRATOS_HENCKY_TRESCA_3D_LAW_H_INCLUDED)
#define  KRATOS_HENCKY_TRESCA_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/non_linear_hencky_plastic_3D_law.hpp"

#include "custom_constitutive/custom_flow_rules/tresca_explicit_plastic_flow_rule.hpp"
#include "custom_constitutive/custom_yield_criteria/tresca_yield_criterion.hpp"
#include "custom_constitutive/custom_hardening_laws/cam_clay_hardening_law.hpp"

namespace Kratos
{
   /**
    * Defines a hyperelastic-plastic isotropic constitutive law J2 in 3D 
    * With stress split in an isochoric and volumetric parts
    * This material law is defined by the parameters needed by the yield criterion:

    * The functionality is limited to large displacements 
    */

   class KRATOS_API(PFEM_SOLID_MECHANICS_APPLICATION) HenckyTresca3DLaw : public NonLinearHenckyElasticPlastic3DLaw
   {
      public:
         /**
          * Type Definitions
          */
         typedef ProcessInfo      ProcessInfoType;
         typedef ConstitutiveLaw         BaseType;
         typedef std::size_t             SizeType;

         typedef FlowRule::Pointer                FlowRulePointer;
         typedef YieldCriterion::Pointer    YieldCriterionPointer;
         typedef HardeningLaw::Pointer        HardeningLawPointer;
         typedef Properties::Pointer            PropertiesPointer;

         /**
          * Counted pointer of HenckyTresca3DLaw
          */

         KRATOS_CLASS_POINTER_DEFINITION( HenckyTresca3DLaw );

         /**
          * Life Cycle
          */

         /**
          * Default constructor.
          */
         HenckyTresca3DLaw();


         HenckyTresca3DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 

         /**
          * Copy constructor.
          */
         HenckyTresca3DLaw (const HenckyTresca3DLaw& rOther);


         /**
          * Assignment operator.
          */

         //HenckyTresca3DLaw& operator=(const HenckyTresca3DLaw& rOther);

         /**
          * Clone function (has to be implemented by any derived class)
          * @return a pointer to a new instance of this constitutive law
          */
         ConstitutiveLaw::Pointer Clone() const override;

         /**
          * Destructor.
          */
         virtual ~HenckyTresca3DLaw();

         /**
          * Operators
          */

         /**
          * Operations needed by the base class:
          */


         /**
          * This function is designed to be called once to perform all the checks needed
          * on the input provided. Checks can be "expensive" as the function is designed
          * to catch user's errors.
          * @param rMaterialProperties
          * @param rElementGeometry
          * @param rCurrentProcessInfo
          * @return
          */
         //int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo);



         /**
          * Input and output
          */
         /**
          * Turn back information as a string.
          */
         //String Info() const override;
         /**
          * Print information about this object.
          */
         //void PrintInfo(std::ostream& rOStream) const override;
         /**
          * Print object's data.
          */
         //void PrintData(std::ostream& rOStream) const override;

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

         void save(Serializer& rSerializer) const override
         {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, NonLinearHenckyElasticPlastic3DLaw )
         }

         void load(Serializer& rSerializer) override
         {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, NonLinearHenckyElasticPlastic3DLaw )
         }



   }; // Class HenckyTresca3DLaw
}  // namespace Kratos.
#endif // KRATOS_HENCKY_TRESCA_3D_LAW_H_INCLUDED defined
