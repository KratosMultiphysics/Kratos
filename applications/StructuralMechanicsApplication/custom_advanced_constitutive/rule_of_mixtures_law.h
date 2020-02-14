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
 * This law combines parallel CL considering the following principles:
 *  - All layer have the same strain
 *  - The total stress is the addition of the strain in each layer
 *  - The constitutive tensor is the addition of the constitutive tensor of each layer
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
     * @param rCombinationFactors The list of subproperties combination factors
     */
    RuleOfMixturesLaw(const std::vector<double>& rCombinationFactors);

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
     * @brief Dimension of the law
     * @details This is not used, so 0 is returned
     */
    SizeType WorkingSpaceDimension() override;

    /**
     * @brief Voigt tensor size
     * @details This is not used, so 0 is returned
     */
    SizeType GetStrainSize() override;

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresInitializeMaterialResponse() override
    {
        return true;
    }

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresFinalizeMaterialResponse() override
    {
        return true;
    }
    
    /**
     * @brief Returns whether this constitutive Law has specified variable (boolean)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<bool>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (integer)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<int>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (double)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<double>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (Vector)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Vector>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (Matrix)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Matrix>& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (array of 3 components)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     * @note Fixed size array of 3 doubles (e.g. for 2D stresses, plastic strains, ...)
     */
    bool Has(const Variable<array_1d<double, 3 > >& rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (array of 6 components)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     * @note Fixed size array of 6 doubles (e.g. for stresses, plastic strains, ...)
     */
    bool Has(const Variable<array_1d<double, 6 > >& rThisVariable) override;

    /**
     * @brief Returns the value of a specified variable (boolean)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    bool& GetValue(
        const Variable<bool>& rThisVariable,
        bool& rValue
        ) override;

    /**
     * @briefReturns the value of a specified variable (integer)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    int& GetValue(
        const Variable<int>& rThisVariable,
        int& rValue
        ) override;

    /**
     * @brief Returns the value of a specified variable (double)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& GetValue(
        const Variable<double>& rThisVariable,
        double& rValue
        ) override;

    /**
     * @brief Returns the value of a specified variable (Vector)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Vector& GetValue(
        const Variable<Vector>& rThisVariable,
        Vector& rValue
        ) override;

    /**
     * @brief Returns the value of a specified variable (Matrix)
     * @param rThisVariable the variable to be returned
     * @return rValue output: the value of the specified variable
     */
    Matrix& GetValue(
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue
        ) override;

    /**
     * @brief Returns the value of a specified variable (array of 3 components)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    array_1d<double, 3 > & GetValue(
        const Variable<array_1d<double, 3 > >& rThisVariable,
        array_1d<double, 3 > & rValue
        ) override;

    /**
     * @brief Returns the value of a specified variable (array of 6 components)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return the value of the specified variable
     */
    array_1d<double, 6 > & GetValue(
        const Variable<array_1d<double, 6 > >& rThisVariable,
        array_1d<double, 6 > & rValue
        ) override;

    /**
     * @brief Sets the value of a specified variable (boolean)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<bool>& rThisVariable,
        const bool& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets the value of a specified variable (integer)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<int>& rThisVariable,
        const int& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets the value of a specified variable (double)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
     void SetValue(
        const Variable<double>& rThisVariable,
        const double& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets the value of a specified variable (Vector)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<Vector >& rThisVariable,
        const Vector& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets the value of a specified variable (Matrix)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<Matrix >& rThisVariable,
        const Matrix& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets the value of a specified variable (array of 3 components)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<array_1d<double, 3 > >& rThisVariable,
        const array_1d<double, 3 > & rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Sets the value of a specified variable (array of 6 components)
     * @param rThisVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<array_1d<double, 6 > >& rThisVariable,
        const array_1d<double, 6 > & rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculates the value of a specified variable (bool)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    bool& CalculateValue(
        Parameters& rParameterValues,
        const Variable<bool>& rThisVariable,
        bool& rValue
        ) override;

    /**
     * @brief Calculates the value of a specified variable (int)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    int& CalculateValue(
        Parameters& rParameterValues,
        const Variable<int>& rThisVariable,
        int& rValue
        ) override;

    /**
     * @brief Calculates the value of a specified variable (double)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    double& CalculateValue(
        Parameters& rParameterValues,
        const Variable<double>& rThisVariable,
        double& rValue
        ) override;

    /**
     * @brief Calculates the value of a specified variable (Vector)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Vector& CalculateValue(
        Parameters& rParameterValues,
        const Variable<Vector>& rThisVariable,
        Vector& rValue
        ) override;

    /**
     * @brief Calculates the value of a specified variable (Matrix)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
     Matrix& CalculateValue(
        Parameters& rParameterValues,
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue
        ) override;

    /**
     * @brief Calculates the value of a specified variable (array of 3 components)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
     array_1d<double, 3 > & CalculateValue(
        Parameters& rParameterValues,
        const Variable<array_1d<double, 3 > >& rVariable,
        array_1d<double, 3 > & rValue
        ) override;

    /**
     * @brief Returns the value of a specified variable (array of 6 components)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return The value of the specified variable
     */
    array_1d<double, 6 > & CalculateValue(
        Parameters& rParameterValues,
        const Variable<array_1d<double, 6 > >& rVariable,
        array_1d<double, 6 > & rValue
        ) override;

     /**
      * @brief Is called to check whether the provided material parameters in the Properties match the requirements of current constitutive model.
      * @param rMaterialProperties the current Properties to be validated against.
      * @return true, if parameters are correct; false, if parameters are insufficient / faulty
      * @note  this has to be implemented by each constitutive model. Returns false in base class since no valid implementation is contained here.
      */
    bool ValidateInput(const Properties& rMaterialProperties) override;

    /**
     * @brief Returns the expected strain measure of this constitutive law (by default linear strains)
     * @return the expected strain measure
     */
    StrainMeasure GetStrainMeasure() override;

    /**
     * @brief Returns the stress measure of this constitutive law (by default 1st Piola-Kirchhoff stress in voigt notation)
     * @return the expected stress measure
     */
    StressMeasure GetStressMeasure() override;

    /**
     * @brief Returns whether this constitutive model is formulated in incremental strains/stresses
     * @note By default, all constitutive models should be formulated in total strains
     * @return true, if formulated in incremental strains/stresses, false otherwise
     */
    bool IsIncremental() override;

    /**
     * @brief This is to be called at the very beginning of the calculation
     * @details (e.g. from InitializeElement) in order to initialize all relevant attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void InitializeMaterial(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues
        ) override;

    /**
     * @brief To be called at the beginning of each solution step
     * @details (e.g. from Element::InitializeSolutionStep)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void InitializeSolutionStep(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief To be called at the end of each solution step
     * @details (e.g. from Element::FinalizeSolutionStep)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void FinalizeSolutionStep(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief  To be called at the beginning of each step iteration
     * @details (e.g. from Element::InitializeNonLinearIteration)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void InitializeNonLinearIteration(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief To be called at the end of each step iteration
     * @details (e.g. from Element::FinalizeNonLinearIteration)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void FinalizeNonLinearIteration(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK1 (Parameters& rValues) override;

    /**
     * @brief Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK2 (Parameters& rValues) override;

    /**
     * @brief Computes the material response in terms of Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseKirchhoff (Parameters& rValues) override;

    /**
     * @brief Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy (Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponsePK1 (Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponsePK2 (Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponseKirchhoff (Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void InitializeMaterialResponseCauchy (Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK1 (Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK2 (Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseKirchhoff (Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy (Parameters& rValues) override;

    /**
     * @brief This can be used in order to reset all internal variables of the
     * constitutive law (e.g. if a model should be reset to its reference state)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void ResetMaterial(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues
        ) override;

    /**
     * @brief This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return 0 if OK, 1 otherwise
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

    std::vector<ConstitutiveLaw::Pointer> mConstitutiveLaws; /// The vector containing the constitutive laws (must be cloned, the ones contained on the properties can conflict between them)
    std::vector<double> mCombinationFactors;                 /// The vector containing the combination factors of the different layers of the material

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
        rSerializer.save("ConstitutiveLaws", mConstitutiveLaws);
        rSerializer.save("CombinationFactors", mCombinationFactors);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw)
        rSerializer.load("ConstitutiveLaws", mConstitutiveLaws);
        rSerializer.load("CombinationFactors", mCombinationFactors);
    }


}; // Class RuleOfMixturesLaw
}  // namespace Kratos.
#endif // KRATOS_RULE_OF_MIXTURES_LAW_H_INCLUDED  defined
