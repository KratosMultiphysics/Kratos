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

#if !defined (KRATOS_ELASTIC_ISOTROPIC_3D_LAW_H_INCLUDED)
#define  KRATOS_ELASTIC_ISOTROPIC_3D_LAW_H_INCLUDED

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
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ElasticIsotropic3D : public ConstitutiveLaw
{
public:

    ///@name Type Definitions
    ///@{

    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of ElasticIsotropic3D
     */

    KRATOS_CLASS_POINTER_DEFINITION( ElasticIsotropic3D );

    ///@}
    ///@name Lyfe Cycle
    ///@{

    /**
     * Default constructor.
     */
    ElasticIsotropic3D();

    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    ElasticIsotropic3D (const ElasticIsotropic3D& rOther);

    /**
     * Destructor.
     */
    ~ElasticIsotropic3D() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function is designed to be called once to check compatibility with element
     * @param rFeatures The Features of the law
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
    * @brief Dimension of the law:
    */
    SizeType WorkingSpaceDimension() override
    {
        return 3;
    };

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return 6;
    };

    /**
     * Computes the material response:
     * PK1 stresses and algorithmic ConstitutiveMatrix
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters & rValues) override;

    /**
     * Computes the material response:
     * PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void CalculateMaterialResponsePK2 (ConstitutiveLaw::Parameters & rValues) override;

    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff (ConstitutiveLaw::Parameters & rValues) override;

    /**
     * Computes the material response:
     * Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters & rValues) override;

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues The internal values of the law
      * @see   Parameters
      */
    void FinalizeMaterialResponsePK1 (ConstitutiveLaw::Parameters & rValues) override;

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues The internal values of the law
      * @see   Parameters
      */
    void FinalizeMaterialResponsePK2 (ConstitutiveLaw::Parameters & rValues) override;

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues The internal values of the law
      * @see   Parameters
      */
    void FinalizeMaterialResponseKirchhoff (ConstitutiveLaw::Parameters & rValues)  override;

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues The internal values of the law
      * @see   Parameters
      */
    void FinalizeMaterialResponseCauchy (ConstitutiveLaw::Parameters & rValues) override;

    /**
     * calculates the value of a specified variable
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue) override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input or that no common error is found.
     * @param rMaterialProperties The properties of the material
     * @param rElementGeometry The geometry of the element
     * @param rCurrentProcessInfo The current process info instance
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

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
    * It checks the size of the constitutive matrix C and resize it if neccessary
    * @param C The constitutive matrix
    */
    void CheckClearElasticMatrix(
        Matrix& C
    );

    /**
    * It calculates the constitutive matrix C
    * @param C The constitutive matrix
    * @param rValues Parameters of the constitutive law
    */
    virtual void CalculateElasticMatrix(
        Matrix& C,
        ConstitutiveLaw::Parameters& rValues
    );

    /**
    * It calculates the stress vector
    * @param rStrainVector The strain vector in Voigt notation
    * @param rStressVector The stress vector in Voigt notation
    * @param rValues Parameters of the constitutive law
    */
    virtual void CalculatePK2Stress(
        const Vector& rStrainVector,
        Vector& rStressVector,
        ConstitutiveLaw::Parameters& rValues
    );

    /**
    * It calculates the strain vector
    * @param rValues The internal values of the law
    * @param rStrainVector The strain vector in Voigt notation
    */
    virtual void CalculateCauchyGreenStrain(
        ConstitutiveLaw::Parameters& rValues,
        Vector& rStrainVector
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


}; // Class ElasticIsotropic3D
}  // namespace Kratos.
#endif // KRATOS_ELASTIC_ISOTROPIC_3D_LAW_H_INCLUDED  defined
