// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined (KRATOS_LINEAR_PLANE_STRAIN_LAW_H_INCLUDED)
#define  KRATOS_LINEAR_PLANE_STRAIN_LAW_H_INCLUDED

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
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) LinearPlaneStrain : public ElasticIsotropic3D
{
public:
    ///@name Type Definitions
    ///@{

    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw       CLBaseType;
    typedef ElasticIsotropic3D      BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of LinearPlaneStrain
     */

    KRATOS_CLASS_POINTER_DEFINITION( LinearPlaneStrain );

    ///@name Life Cycle
    ///@{

    /**
     * Default constructor.
     */
    LinearPlaneStrain();

    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    LinearPlaneStrain (const LinearPlaneStrain& rOther);


    /**
     * Destructor.
     */
    ~LinearPlaneStrain() override;

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
     * Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override
    {
        return 2;
    };

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return 3;
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

    /**
     * returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    bool& GetValue(const Variable<bool>& rThisVariable, bool& rValue) override;

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
    * It calculates the constitutive matrix C
    * @param C: The constitutive matrix
    * @param rValues Parameters of the constitutive law
    */
    void CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues) override;

    /**
    * It calculates the stress vector
    * @param rStrainVector The strain vector in Voigt notation
    * @param rStressVector The stress vector in Voigt notation
    * @param rValues Parameters of the constitutive law
    */
    void CalculatePK2Stress(
        const Vector& rStrainVector,
        Vector& rStressVector,
        ConstitutiveLaw::Parameters& rValues
    ) override;

    /**
    * It calculates the strain vector
    * @param rValues The internal values of the law
    * @param rStrainVector The strain vector in Voigt notation
    */
    void CalculateCauchyGreenStrain(
        ConstitutiveLaw::Parameters& rValues,
        Vector& rStrainVector
    ) override;

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ElasticIsotropic3D)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ElasticIsotropic3D)
    }


}; // Class LinearPlaneStrain
}  // namespace Kratos.
#endif // KRATOS_LINEAR_PLANE_STRAIN_LAW_H_INCLUDED  defined
