//
//   Project Name:        KratosUmatApplication      $
//   Created by:          $Author:       LlMonforte  $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:      September 2017 $
//   Revision:            $Revision:             0.0 $
//
//

#if !defined(KRATOS_LARGE_STRAIN_UMAT_3D_LAW_H_INCLUDED )
#define  KRATOS_LARGE_STRAIN_UMAT_3D_LAW_H_INCLUDED

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

  class KRATOS_API(UMAT_APPLICATION) LargeStrainUmat3DLaw : public Umat3DLaw
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

    /// Pointer definition of LargeStrainUmat3DLaw
    KRATOS_CLASS_POINTER_DEFINITION( LargeStrainUmat3DLaw );
    
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    LargeStrainUmat3DLaw();

    /// Copy constructor.
    LargeStrainUmat3DLaw(const LargeStrainUmat3DLaw& rOther);

    /// Clone.
    ConstitutiveLaw::Pointer Clone() const override;

    /// Assignment operator.
    LargeStrainUmat3DLaw& operator=(const LargeStrainUmat3DLaw& rOther);

    /// Destructor.
    virtual ~LargeStrainUmat3DLaw();

    ///@}
    ///@name Operators
    ///@{

    
    ///@}
    ///@name Operations
    ///@{

    /**
     * Material parameters are inizialized
     */
    virtual void InitializeMaterial ( const Properties& props,
				      const GeometryType& geom,
				      const Vector& ShapeFunctionsValues ) override;

    /**
     * Computes the material response:
     * Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    virtual void CalculateMaterialResponseKirchhoff (Parameters & rValues) override;

    
    /**
     * Updates the material response:
     * Cauchy stresses and Internal Variables
     * @param rValues
     * @see   Parameters
     */
    virtual void FinalizeMaterialResponseCauchy(Parameters & rValues) override;


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
        buffer << "LargeStrainUmat3DLaw";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "LargeStrainUmat3DLaw";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "LargeStrainUmat3DLaw Data";
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
  }; // Class LargeStrainUmat3DLaw
  
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_LARGE_STRAIN_UMAT_3D_LAW_H_INCLUDED  defined 

