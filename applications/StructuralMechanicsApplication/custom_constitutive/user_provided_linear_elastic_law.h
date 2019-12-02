// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Riccardo Rossi
//

#if !defined (KRATOS_USER_PROVIDED_LINEAR_ELASTIC_LAW_H_INCLUDED)
#define  KRATOS_USER_PROVIDED_LINEAR_ELASTIC_LAW_H_INCLUDED

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
 * @class UserProvidedLinearElasticLaw
 * @ingroup StructuralMechanicsApplication
 * @brief This class defines a linear elastic law with user provided constitutive tensor
 * @details The constitutive behaviour is given by the constitutive tensor defined in Matrix
 * type variable ELASTICITY_TENSOR, which is assumed to be stored in the properties
 * @author Ruben Zorrilla
 * @author Riccardo Rossi
 */
template <unsigned int TDim>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) UserProvidedLinearElasticLaw
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

    /// Static definition of the dimension
    static constexpr SizeType Dimension = TDim;

    /// Static definition of the StrainSize
    static constexpr SizeType StrainSize = (TDim * 3) - 3;

    /// Counted pointer of UserProvidedLinearElasticLaw
    KRATOS_CLASS_POINTER_DEFINITION( UserProvidedLinearElasticLaw );

    ///@}
    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    UserProvidedLinearElasticLaw();

    /**
     * @brief Clone method
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    UserProvidedLinearElasticLaw (const UserProvidedLinearElasticLaw& rOther);

    /**
     * @brief Destructor.
     */
    ~UserProvidedLinearElasticLaw();

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
        return Dimension;
    };

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return StrainSize;
    };

    /**
     * @brief Returns the expected strain measure of this constitutive law (by default Green-Lagrange)
     * @return the expected strain measure
     */
    StrainMeasure GetStrainMeasure() override
    {
        return StrainMeasure_Infinitesimal;
    }

    /**
     * @brief Returns the stress measure of this constitutive law (by default 2st Piola-Kirchhoff stress in voigt notation)
     * @return the expected stress measure
     */
    StressMeasure GetStressMeasure() override
    {
        return StressMeasure_Cauchy;
    }

    /**
     * @brief Computes the material response:
     * @details PK1 stresses and algorithmic ConstitutiveMatrix
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters & rValues) override;

    /**
     * @brief Computes the material response:
     * @details PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void CalculateMaterialResponsePK2 (ConstitutiveLaw::Parameters & rValues) override;

    /**
     * @brief Computes the material response:
     * @details Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff (ConstitutiveLaw::Parameters & rValues) override;

    /**
     * @brief Computes the material response:
     * @details Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues The internal values of the law
     * @see   Parameters
     */
    void CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters & rValues) override;

    /**
     * Initialize the material response in terms of 1st Piola-Kirchhoff stresses
     * @param rValues The internal values of the law
     * @see Parameters
     */
    void InitializeMaterialResponsePK1 (Parameters& rValues) override {};

    /**
     * Initialize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @param rValues The internal values of the law
     * @see Parameters
     */
    void InitializeMaterialResponsePK2 (Parameters& rValues) override {};

    /**
     * Initialize the material response in terms of Kirchhoff stresses
     * @param rValues The internal values of the law
     * @see Parameters
     */
    void InitializeMaterialResponseKirchhoff (Parameters& rValues) override {};

    /**
     * Initialize the material response in terms of Cauchy stresses
     * @param rValues The internal values of the law
     * @see Parameters
     */
    void InitializeMaterialResponseCauchy (Parameters& rValues) override {};

    /**
     * Finalize the material response in terms of 1st Piola-Kirchhoff stresses
     * @param rValues The internal values of the law
     * @see Parameters
     */
    void FinalizeMaterialResponsePK1 (Parameters& rValues) override {};

    /**
     * Finalize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @param rValues The internal values of the law
     * @see Parameters
     */
    void FinalizeMaterialResponsePK2 (Parameters& rValues) override {};

    /**
     * Finalize the material response in terms of Kirchhoff stresses
     * @param rValues The internal values of the law
     * @see Parameters
     */
    void FinalizeMaterialResponseKirchhoff (Parameters& rValues) override {};

    /**
     * Finalize the material response in terms of Cauchy stresses
     * @param rValues The internal values of the law
     * @see Parameters
     */

    void FinalizeMaterialResponseCauchy (Parameters& rValues) override {};

    /**
     * @brief It calculates the value of a specified variable (double case)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<double>& rThisVariable,
        double& rValue) override;

    /**
     * @brief It calculates the value of a specified variable (Vector case)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Vector& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Vector>& rThisVariable,
        Vector& rValue) override;

    /**
     * @brief It calculates the value of a specified variable (Matrix case)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Matrix& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue) override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning of the calculations, so to verify that nothing is missing from the input or that no common error is found.
     * @param rMaterialProperties The properties of the material
     * @param rElementGeometry The geometry of the element
     * @param rCurrentProcessInfo The current process info instance
     * @return 0 if OK, 1 otherwise
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo) override;

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
    * @brief It calculates the constitutive matrix rConstitutiveMatrix
    * @param rConstitutiveMatrix The constitutive matrix
    * @param rValues Parameters of the constitutive law
    */
    virtual void CalculateElasticMatrix(
        Matrix& rConstitutiveMatrix,
        ConstitutiveLaw::Parameters& rValues
        );

    /**
     * @brief It calculates the stress vector
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
     * @brief It calculates the strain vector
     * @param rValues The internal values of the law
     * @param rStrainVector The strain vector in Voigt notation
     */
    virtual void CalculateGreenLagrangeStrainVector(
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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
    }


}; // Class UserProvidedLinearElasticLaw
}  // namespace Kratos.
#endif // KRATOS_USER_PROVIDED_LINEAR_ELASTIC_LAW_H_INCLUDED  defined
