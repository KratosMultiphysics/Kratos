// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#if !defined(KRATOS_SMALL_STRAIN_UDSM_2D_INTERFACE_LAW_H_INCLUDED )
#define  KRATOS_SMALL_STRAIN_UDSM_2D_INTERFACE_LAW_H_INCLUDED

// System includes
#include "includes/define.h"

// External includes

// Project includes
#include "small_strain_udsm_3D_law.hpp"

namespace Kratos
{
   ///@addtogroup ConstitutiveModelsApplication
   ///@{

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

   /// Short class definition.
   /** Detail class definition.
    */
   class KRATOS_API(GEO_MECHANICS_APPLICATION) SmallStrainUDSM2DInterfaceLaw: public SmallStrainUDSM3DLaw
   {
   public:

      // The base class ConstitutiveLaw type definition
      using BaseType = ConstitutiveLaw;

      /// The size type definition
      using SizeType = std::size_t;

      /// Static definition of the dimension
      static constexpr SizeType Dimension = N_DIM_2D;

      /// Static definition of the VoigtSize
      static constexpr SizeType VoigtSize = VOIGT_SIZE_2D_INTERFACE;

      /// Pointer definition of SmallStrainUDSM2DInterfaceLaw
      KRATOS_CLASS_POINTER_DEFINITION( SmallStrainUDSM2DInterfaceLaw );


      //@}
      //@name Life Cycle
      //@{

      //----------------------------------------------------------------------------------------
      /**
       * @brief Default constructor.
       */
      SmallStrainUDSM2DInterfaceLaw();

      /**
       * @brief Clone method
       */
      ConstitutiveLaw::Pointer Clone() const override;

      /**
       * Copy constructor.
       */
      SmallStrainUDSM2DInterfaceLaw(SmallStrainUDSM2DInterfaceLaw const& rOther);

      /**
       * @brief Destructor.
       */
      virtual ~SmallStrainUDSM2DInterfaceLaw();

      // Assignment operator:
      SmallStrainUDSM2DInterfaceLaw& operator=(SmallStrainUDSM2DInterfaceLaw const& rOther);

      Vector& GetValue( const Variable<Vector> &rThisVariable, Vector &rValue ) override;

      void SetValue( const Variable<Vector>& rVariable,
                     const Vector& rValue,
                     const ProcessInfo& rCurrentProcessInfo ) override;

      //----------------------------------------------------------------------------------------
      /**
       * @brief Dimension of the law:
       */
      SizeType WorkingSpaceDimension() override
      {
         return Dimension;
      }

      /**
       * @brief Voigt tensor size:
       */
      SizeType GetStrainSize() const override
      {
         return VoigtSize;
      }

      /**
       * @brief Returns the expected strain measure of this constitutive law (by default Green-Lagrange)
       * @return the expected strain measure
       */
      StrainMeasure GetStrainMeasure() override
      {
         return StrainMeasure_Infinitesimal;
      }

      /**
       * returns the stress measure of this constitutive law (by default 1st Piola-Kirchhoff stress in voigt notation)
       * @return the expected stress measure
       */
      virtual StressMeasure GetStressMeasure() override
      {
         return StressMeasure_Cauchy;
      }


      /**
       * @brief It calculates the strain vector
       * @param rValues The internal values of the law
       * @param rStrainVector The strain vector in Voigt notation
       */
      void CalculateCauchyGreenStrain(ConstitutiveLaw::Parameters& rValues, Vector& rStrainVector) override;

      ///@}
      ///@name Inquiry
      ///@{


      ///@}
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const override
      {
         std::stringstream buffer;
         buffer << "SmallStrainUDSM2DInterfaceLaw";
         return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const override
      {
         rOStream << "SmallStrainUDSM2DInterfaceLaw";
      }

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const override
      {
         rOStream << "SmallStrainUDSM2DInterfaceLaw Data";
      }


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

      ///@}
      ///@name Protected  Access
      ///@{
      void UpdateInternalDeltaStrainVector(ConstitutiveLaw::Parameters &rValues) override;
      void SetExternalStressVector(Vector& rStressVector) override;
      void SetInternalStressVector(const Vector& rStressVector) override;
      void SetInternalStrainVector(const Vector& rStrainVector) override;
      void CopyConstitutiveMatrix(ConstitutiveLaw::Parameters &rValues, Matrix& rConstitutiveMatrix) override;


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

      indexStress3D getIndex3D(indexStress2DInterface index2D);

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
      ///@name Serialization
      ///@{
      friend class Serializer;

      virtual void save(Serializer& rSerializer) const override
      {
         KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
      }

      virtual void load(Serializer& rSerializer) override
      {
         KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
      }

      ///@}
      ///@name Private Inquiry
      ///@{


      ///@}
      ///@name Un accessible methods
      ///@{

      ///@}

   }; // Class SmallStrainUDSM3DLaw

   ///@}

   ///@name Type Definitions
   ///@{


   ///@}
   ///@name Input and output
   ///@{

   ///@}

   ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SMALL_STRAIN_UDSM_2D_INTERFACE_LAW_H_INCLUDED  defined


