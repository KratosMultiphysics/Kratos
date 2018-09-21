// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Fernando Rastellini
//

#if !defined (KRATOS_RULE_OF_MIXTURES_LAW_H_INCLUDED)
#define  KRATOS_RULE_OF_MIXTURES_LAW_H_INCLUDED

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
 * @class RuleOfMixturesLaw
 * @ingroup StructuralMechanicsApplication
 * @brief This law defines a parallel rule of mixture (classic law of mixture)
 * @details The constitutive law show have defined a subproperties in order to work properly
 * @author Vicente Mataix Ferrandiz
 * @author Fernando Rastellini
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) RuleOfMixturesLaw
    : public ConstitutiveLaw
{
public:

    ///@name Type Definitions
    ///@{

    /// The definition of the process info
    typedef ProcessInfo ProcessInfoType;
    
    /// The definition of the CL base  class
    typedef ConstitutiveLaw    BaseType;
    
    /// The definition of the size type
    typedef std::size_t        SizeType;
    
    /// The definition of the index type
    typedef std::size_t       IndexType;
    
    /// Pointer definition of RuleOfMixturesLaw
    KRATOS_CLASS_POINTER_DEFINITION( RuleOfMixturesLaw );

    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    RuleOfMixturesLaw();

    /**
     * @brief Constructor with values
     * @param rSubPropertiesIDs The list of subproperties ids
     * @param rCombinationFactors The list of subproperties combination factors
     * @param rMaterialRotationAngles The rotation angles of the layers
     */
    RuleOfMixturesLaw(
        const std::vector<IndexType>& rSubPropertiesIDs,
        const std::vector<double>& rCombinationFactors,
        const std::vector<double>& rMaterialRotationAngles
        );

    /**
     * @brief Copy constructor.
     */
    RuleOfMixturesLaw (const RuleOfMixturesLaw& rOther);

    /**
     * @brief Destructor.
     */
    ~RuleOfMixturesLaw() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Clone operation
     */
    ConstitutiveLaw::Pointer Clone() const override;
    
    /**
     * @brief Creates a new constitutive law pointer
     * @param NewParameters The configuration parameters of the new constitutive law
     * @return a Pointer to the new constitutive law
     */
    ConstitutiveLaw::Pointer Create(Kratos::Parameters NewParameters) const override;

    /**
     * @brief This function is designed to be called once to check compatibility with element
     * @param rFeatures The Features of the law
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * @brief Dimension of the law
     * @details This is not used, so 0 is returned
     */
    SizeType WorkingSpaceDimension() override
    {
        return 0;
    };

    /**
     * @brief Voigt tensor size
     * @details This is not used, so 0 is returned
     */
    SizeType GetStrainSize() override
    {
        return 0;
    };

    /**
     * @brief This is to be called at the very beginning of the calculation (e.g. from InitializeElement) in order to initialize all relevant attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void InitializeMaterial(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues
        );

    /**
     * @brief Computes the material response: PK1 stresses and algorithmic ConstitutiveMatrix
     * @param rValues The Internalvalues of the law
     * @see   Parameters
     */
    void CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters & rValues) override;

    /**
     * @brief Computes the material response: PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues The Internalvalues of the law
     * @see   Parameters
     */
    void CalculateMaterialResponsePK2 (ConstitutiveLaw::Parameters & rValues) override;

    /**
     * @brief Computes the material response: Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues The Internalvalues of the law
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff (ConstitutiveLaw::Parameters & rValues) override;

    /**
     * @brief Computes the material response: Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues The Internalvalues of the law
     * @see   Parameters
     */
    void CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters & rValues) override;

    /**
      * @brief Updates the material response: Cauchy stresses and Internal Variables
      * @param rValues The Internalvalues of the law
      * @see   Parameters
      */
    void FinalizeMaterialResponsePK1 (ConstitutiveLaw::Parameters & rValues) override;

    /**
      * @brief Updates the material response: Cauchy stresses and Internal Variables
      * @param rValues The Internalvalues of the law
      * @see   Parameters
      */
    void FinalizeMaterialResponsePK2 (ConstitutiveLaw::Parameters & rValues) override;

    /**
      * @brief Updates the material response: Cauchy stresses and Internal Variables
      * @param rValues The Internalvalues of the law
      * @see   Parameters
      */
    void FinalizeMaterialResponseKirchhoff (ConstitutiveLaw::Parameters & rValues)  override;

    /**
      * @brief Updates the material response: Cauchy stresses and Internal Variables
      * @param rValues The Internalvalues of the law
      * @see   Parameters
      */
    void FinalizeMaterialResponseCauchy (ConstitutiveLaw::Parameters & rValues) override;

    /**
     * @brief It calculates the value of a specified variable (double)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    double& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<double>& rThisVariable, 
        double& rValue
        ) override;

    /**
     * @brief It calculates the value of a specified variable (Vector)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Vector& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Vector>& rThisVariable, 
        Vector& rValue
        ) override;

    /**
     * @brief It calculates the value of a specified variable (Matrix)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Matrix& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Matrix>& rThisVariable, 
        Matrix& rValue
        ) override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
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

    ///@}

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    std::unordered_map<IndexType, ConstitutiveLaw::Pointer> mConstitutiveLaws; /// The map containing the constitutive laws (must be cloned, the ones contained on the properties can conflict between them)
    std::unordered_map<IndexType, double> mCombinationFactors;                 /// The map containing the combination factors of the different layers of the material
    std::unordered_map<IndexType, double> mMaterialRotationAngles;             /// The map containing the rotation angles of the different layers

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
        rSerializer.save("CombinationFactors", mCombinationFactors);
        rSerializer.save("MaterialRotationAngles", mMaterialRotationAngles);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw)
        rSerializer.load("CombinationFactors", mCombinationFactors);
        rSerializer.load("MaterialRotationAngles", mMaterialRotationAngles);
    }


}; // Class RuleOfMixturesLaw
}  // namespace Kratos.
#endif // KRATOS_RULE_OF_MIXTURES_LAW_H_INCLUDED  defined
