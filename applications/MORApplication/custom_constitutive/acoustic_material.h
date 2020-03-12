// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: 
//
//  Main authors:    
//                   
//

#if !defined (ACOUSTIC_MATERIAL_H_INCLUDED)
#define  ACOUSTIC_MATERIAL_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

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

/**
 * @class AcousticMaterial
 * @ingroup MORApplication
 * @brief 
 * @details 
 * @author 
 */
class KRATOS_API(MOR_APPLICATION) AcousticMaterial
    : public ConstitutiveLaw
{
public:

    ///@name Type Definitions
    ///@{

    /// The process info type definition
    typedef ProcessInfo      ProcessInfoType;

    /// The base class ConstitutiveLaw type definition
    typedef ConstitutiveLaw         BaseType;

    /// The size type definition
    typedef std::size_t             SizeType;

    // Adding the respective using to avoid overload conflicts
    using BaseType::Has;
    using BaseType::GetValue;

    /// Counted pointer of AcousticMaterial
    KRATOS_CLASS_POINTER_DEFINITION( AcousticMaterial );

    ///@}
    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    AcousticMaterial();

    /**
     * @brief Clone method
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    AcousticMaterial (const AcousticMaterial& rOther);

    /**
     * @brief Destructor.
     */
    ~AcousticMaterial() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

 

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
    * @brief It checks the size of the constitutive matrix rConstitutiveMatrix and resize it if neccessary
    * @param rConstitutiveMatrix The constitutive matrix
    */
    void CheckClearElasticMatrix(Matrix& rConstitutiveMatrix);

    /**
    * @brief It calculates the constitutive matrix rConstitutiveMatrix
    * @param rConstitutiveMatrix The constitutive matrix
    * @param rValues Parameters of the constitutive law
    */
    virtual void CalculateElasticMatrix(
        Matrix& rConstitutiveMatrix,
        ConstitutiveLaw::Parameters& rValues
        );


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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw)
    }


}; // Class AcousticMaterial
}  // namespace Kratos.
#endif // ACOUSTIC_MATERIAL_H_INCLUDED  defined
