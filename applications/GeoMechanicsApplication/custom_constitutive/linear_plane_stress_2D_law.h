// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi,
//                   Vahid Galavi
//

#if !defined (KRATOS_GEO_LINEAR_PLANE_STRESS_LAW_H_INCLUDED)
#define  KRATOS_GEO_LINEAR_PLANE_STRESS_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/elastic_isotropic_K0_3d_law.h"

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
 * @class LinearPlaneStress2DLaw
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a small deformation linear elastic constitutive model for plane stress cases
 * @details This class derives from the linear elastic case on 3D
 * @author Riccardo Rossi
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) LinearPlaneStress2DLaw 
    : public ElasticIsotropicK03DLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// The base class ConstitutiveLaw type definition
    typedef ConstitutiveLaw       CLBaseType;
    
    /// The base class ElasticIsotropicK03DLaw type definition
    typedef ElasticIsotropicK03DLaw      BaseType;

    // Adding the respective using to avoid overload conflicts
    using BaseType::Has;
    using BaseType::GetValue;
    
    /// The size type definition
    typedef std::size_t             SizeType;
    
    /// Static definition of the dimension
    static constexpr SizeType Dimension = 2;
    
    /// Static definition of the VoigtSize
    static constexpr SizeType VoigtSize = 3;
    
    /// Counted pointer of LinearPlaneStress2DLaw
    KRATOS_CLASS_POINTER_DEFINITION( LinearPlaneStress2DLaw );

    ///@name Life Cycle
    ///@{

    /**
     * Default constructor.
     */
    LinearPlaneStress2DLaw();

    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    LinearPlaneStress2DLaw (const LinearPlaneStress2DLaw& rOther);


    /**
     * Destructor.
     */
    ~LinearPlaneStress2DLaw() override;

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
        return Dimension;
    };

    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return VoigtSize;
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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ElasticIsotropicK03DLaw)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ElasticIsotropicK03DLaw)
    }


}; // Class LinearPlaneStress2DLaw
}  // namespace Kratos.
#endif // KRATOS_GEO_LINEAR_PLANE_STRESS_LAW_H_INCLUDED  defined
