// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined (KRATOS_AXISYM_ELASTIC_ISOTROPIC_LAW_H_INCLUDED)
#define  KRATOS_AXISYM_ELASTIC_ISOTROPIC_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/elastic_isotropic_3d.h"

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

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) AxisymElasticIsotropic : public ElasticIsotropic3D
{
public:

    ///@name Type Definitions
    ///@{

    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;

    /**
     * Counted pointer of AxisymElasticIsotropic
     */
    KRATOS_CLASS_POINTER_DEFINITION( AxisymElasticIsotropic );

    ///@name Lyfe Cycle
    ///@{

    /**
     * Default constructor.
     */
    AxisymElasticIsotropic();

    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    AxisymElasticIsotropic (const AxisymElasticIsotropic& rOther);

    /**
     * Destructor.
     */
    ~AxisymElasticIsotropic() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures: The Features of the law
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return 4;
    }

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

    /**
     * It calculates the constitutive matrix C
     * @param C: The constitutive matrix
     * @param E: The Young Modulus
     * @param NU: The poisson coefficient
     */
    void CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues) override;

    /**
     * It calculates the strain vector
     * @param rValues: The Internalvalues of the law
     * @param rStrainVector: The strain vector in Voigt notation
     */
    void CalculateCauchyGreenStrain(
        ConstitutiveLaw::Parameters& rValues,
        Vector& rStrainVector
    ) override;

    ///@}
    ///@name Private Operations
    ///@{
    ///@}


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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw)
    }


}; // Class AxisymElasticIsotropic
}  // namespace Kratos.
#endif // KRATOS_AXISYM_ELASTIC_ISOTROPIC_LAW_H_INCLUDED  defined
