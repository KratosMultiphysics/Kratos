//
//   Project Name:        KratosUmatApplication      $
//   Created by:          $Author:       LlMonforte  $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:      September 2017 $
//   Revision:            $Revision:             0.0 $
//
//

#if !defined(KRATOS_SMALL_STRAIN_UMAT_3D_LAW_H_INCLUDED )
#define  KRATOS_SMALL_STRAIN_UMAT_3D_LAW_H_INCLUDED

// System includes
#include <iostream>
#include <cmath>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/properties.h"
#include "custom_laws/umat_3D_law.hpp"


namespace Kratos
{

  class KRATOS_API(UMAT_APPLICATION) SmallStrainUmat3DLaw : public Umat3DLaw
  {
  protected:

    using VoigtIndexType = const unsigned int(*)[2];

  public:

    ///@name Type Definitions
    ///@{
    typedef ProcessInfo                                           ProcessInfoType;
    typedef ConstitutiveLaw                                              BaseType;
    typedef Umat3DLaw                                                 DerivedType;

    typedef DerivedType::MaterialTensorType                    MaterialTensorType;
    typedef DerivedType::PlaneArrayType                            PlaneArrayType;
    typedef DerivedType::SpaceArrayType                            SpaceArrayType;

    /// Pointer definition of SmallStrainUmat3DLaw
    KRATOS_CLASS_POINTER_DEFINITION( SmallStrainUmat3DLaw );
    
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    SmallStrainUmat3DLaw();

    /// Copy constructor.
    SmallStrainUmat3DLaw(const SmallStrainUmat3DLaw& rOther);

    /// Clone.
    ConstitutiveLaw::Pointer Clone() const override;

    /// Assignment operator.
    SmallStrainUmat3DLaw& operator=(const SmallStrainUmat3DLaw& rOther);

    /// Destructor.
    virtual ~SmallStrainUmat3DLaw();

    ///@}
    ///@name Operators
    ///@{

    
    ///@}
    ///@name Operations
    ///@{


    /**
     * This function is designed to be called once to check compatibility with element and the constitutive law
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures) override;


    ///@}
    ///@name Access
    ///@{

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
        buffer << "SmallStrainUmatLaw";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SmallStrainUmat3DLaw";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "SmallStrainUmat3DLaw Data";
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
    ///@name Serialization
    ///@{
    friend class Serializer;

    virtual void save ( Serializer& rSerializer ) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS ( rSerializer, Umat3DLaw )
    }

    virtual void load ( Serializer& rSerializer ) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS ( rSerializer, Umat3DLaw )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
  }; // Class SmallStrainUmat3DLaw
  
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_SMALL_STRAIN_UMAT_3D_LAW_H_INCLUDED  defined 

