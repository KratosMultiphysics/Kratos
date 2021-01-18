// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Klaus B. Sautter
//

#if !defined (KRATOS_WRINKLING_LINEAR_2D_LAW_H_INCLUDED)
#define  KRATOS_WRINKLING_LINEAR_2D_LAW_H_INCLUDED

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
 * @class WrinklingLinear2DLaw
 * @ingroup StructuralMechanicsApplication
 * @brief This law defines a wrinkling modification for any linear 2D claw
 * @author Klaus B. Sautter
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) WrinklingLinear2DLaw
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

    /// Pointer definition of WrinklingLinear2DLaw
    KRATOS_CLASS_POINTER_DEFINITION( WrinklingLinear2DLaw );

    enum class WrinklingType {
      Taut,
      Slack,
      Wrinkle
    };

    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    WrinklingLinear2DLaw();

    /**
     * @brief Copy constructor.
     */
    WrinklingLinear2DLaw (const WrinklingLinear2DLaw& rOther);

    /**
     * @brief Destructor.
     */
    ~WrinklingLinear2DLaw() override = default;

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
     * @brief Returns whether this constitutive Law has specified variable (bool)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<bool>& rThisVariable) override
    {
        return mpConstitutiveLaw->Has(rThisVariable);
    }

    /**
     * @brief Returns whether this constitutive Law has specified variable (int)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<int>& rThisVariable) override
    {
        return mpConstitutiveLaw->Has(rThisVariable);
    }

    /**
     * @brief Returns whether this constitutive Law has specified variable (double)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<double>& rThisVariable) override
    {
        return mpConstitutiveLaw->Has(rThisVariable);
    }

    /**
     * @brief Returns whether this constitutive Law has specified variable (Vector)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Vector>& rThisVariable) override
    {
        return mpConstitutiveLaw->Has(rThisVariable);
    }

    /**
     * @brief Returns whether this constitutive Law has specified variable (Matrix)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Matrix>& rThisVariable) override
    {
        return mpConstitutiveLaw->Has(rThisVariable);
    }

    /**
     * @brief Returns whether this constitutive Law has specified variable (array_1d<double, 3 >)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<array_1d<double, 3 > >& rThisVariable) override
    {
        return mpConstitutiveLaw->Has(rThisVariable);
    }

    /**
     * @brief Returns whether this constitutive Law has specified variable (array_1d<double, 6 >)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<array_1d<double, 6 > >& rThisVariable) override
    {
        return mpConstitutiveLaw->Has(rThisVariable);
    }

    /**
     * @brief Returns the value of a specified variable (bool)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    bool& GetValue(const Variable<bool>& rThisVariable,bool& rValue) override
    {
        mpConstitutiveLaw->GetValue(rThisVariable,rValue);
        return rValue;
    }

    /**
     * @brief Returns the value of a specified variable (int)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    int& GetValue(const Variable<int>& rThisVariable,int& rValue) override
    {
        mpConstitutiveLaw->GetValue(rThisVariable,rValue);
        return rValue;
    }

    /**
     * @brief Returns the value of a specified variable (double)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& GetValue(const Variable<double>& rThisVariable,double& rValue) override
    {
        mpConstitutiveLaw->GetValue(rThisVariable,rValue);
        return rValue;
    }

    /**
     * @brief Returns the value of a specified variable (Vector)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Vector& GetValue(const Variable<Vector>& rThisVariable,Vector& rValue) override
    {
        mpConstitutiveLaw->GetValue(rThisVariable,rValue);
        return rValue;
    }

    /**
     * @brief Returns the value of a specified variable (Matrix)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Matrix& GetValue(const Variable<Matrix>& rThisVariable,Matrix& rValue) override
    {
        mpConstitutiveLaw->GetValue(rThisVariable,rValue);
        return rValue;
    }

    /**
     * @brief Returns the value of a specified variable (array_1d<double, 3 > )
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    array_1d<double, 3 >& GetValue(const Variable<array_1d<double, 3 >>& rThisVariable,array_1d<double, 3 >& rValue) override
    {
        mpConstitutiveLaw->GetValue(rThisVariable,rValue);
        return rValue;
    }


    /**
     * @brief Returns the value of a specified variable (array_1d<double, 6 > )
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    array_1d<double, 6 >& GetValue(const Variable<array_1d<double, 6 >>& rThisVariable,array_1d<double, 6 >& rValue) override
    {
        mpConstitutiveLaw->GetValue(rThisVariable,rValue);
        return rValue;
    }

    /**
     * @brief Calculates the value of a specified variable (bool)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    bool& CalculateValue(Parameters& rParameterValues,const Variable<bool>& rThisVariable,
        bool& rValue) override
    {
        mpConstitutiveLaw->CalculateValue(rParameterValues,rThisVariable,rValue);
        return rValue;
    }

    /**
     * @brief Calculates the value of a specified variable (int)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    int& CalculateValue(Parameters& rParameterValues,
        const Variable<int>& rThisVariable,int& rValue) override
    {
        mpConstitutiveLaw->CalculateValue(rParameterValues,rThisVariable,rValue);
        return rValue;
    }

    /**
     * @brief Calculates the value of a specified variable (double)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    double& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<double>& rThisVariable,double& rValue) override
    {
        mpConstitutiveLaw->CalculateValue(rParameterValues,rThisVariable,rValue);
        return rValue;
    }

    /**
     * @brief Calculates the value of a specified variable (Vector)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Vector& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Vector>& rThisVariable,Vector& rValue) override
    {
        mpConstitutiveLaw->CalculateValue(rParameterValues,rThisVariable,rValue);
        return rValue;
    }

    /**
     * @brief Calculates the value of a specified variable (Matrix)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Matrix& CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Matrix>& rThisVariable,Matrix& rValue) override
    {
        mpConstitutiveLaw->CalculateValue(rParameterValues,rThisVariable,rValue);
        return rValue;
    }

    /**
     * @brief Calculates the value of a specified variable (array_1d<double, 3 >)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    array_1d<double, 3 >& CalculateValue(Parameters& rParameterValues,
        const Variable<array_1d<double, 3 >>& rThisVariable,array_1d<double, 3 >& rValue) override
    {
        mpConstitutiveLaw->CalculateValue(rParameterValues,rThisVariable,rValue);
        return rValue;
    }

    /**
     * @brief Calculates the value of a specified variable (array_1d<double, 6 >)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    array_1d<double, 6 >& CalculateValue(Parameters& rParameterValues,
        const Variable<array_1d<double, 6 >>& rThisVariable,array_1d<double, 6 >& rValue) override
    {
        mpConstitutiveLaw->CalculateValue(rParameterValues,rThisVariable,rValue);
        return rValue;
    }

    /**
     * @brief Sets the value of a specified variable (boolean)
     * @param rVariable the variable to be returned
     * @param Value new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(const Variable<bool>& rThisVariable,const bool& rValue,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }

    /**
     * @brief Sets the value of a specified variable (int)
     * @param rVariable the variable to be returned
     * @param Value new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(const Variable<int>& rThisVariable,const int& rValue,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }

    /**
     * @brief Sets the value of a specified variable (double)
     * @param rVariable the variable to be returned
     * @param Value new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(const Variable<double>& rThisVariable,const double& rValue,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }

    /**
     * @brief Sets the value of a specified variable (Vector)
     * @param rVariable the variable to be returned
     * @param Value new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(const Variable<Vector>& rThisVariable,const Vector& rValue,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }

    /**
     * @brief Sets the value of a specified variable (Matrix)
     * @param rVariable the variable to be returned
     * @param Value new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(const Variable<Matrix>& rThisVariable,const Matrix& rValue,
        const ProcessInfo& rCurrentProcessInfo) override
    {
        mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }

    /**
     * @brief Sets the value of a specified variable (array_1d<double, 3 > )
     * @param rVariable the variable to be returned
     * @param Value new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(const Variable<array_1d<double, 3 >>& rThisVariable,
        const array_1d<double, 3 >& rValue,const ProcessInfo& rCurrentProcessInfo) override
    {
        mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }

    /**
     * @brief Sets the value of a specified variable (array_1d<double, 6 > )
     * @param rVariable the variable to be returned
     * @param Value new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(const Variable<array_1d<double, 6 >>& rThisVariable,
        const array_1d<double, 6 >& rValue,const ProcessInfo& rCurrentProcessInfo) override
    {
        mpConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }

    /**
     * This is to be called at the very beginning of the calculation
     * (e.g. from InitializeElement) in order to initialize all relevant
     * attributes of the constitutive law
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
     * Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK2 (Parameters& rValues) override;

    /**
     * Initialize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponsePK2 (Parameters& rValues) override;


    /**
     * Finalize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK2 (Parameters& rValues) override;


    /**
     * This can be used in order to reset all internal variables of the
     * constitutive law (e.g. if a model should be reset to its reference state)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void ResetMaterial(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues
        ) override;

    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures) override;


    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rMaterialProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Calculates the principal vectors
     * @param rPrincipalVector the principal vectors
     * @param rNonPrincipalVector reference state
     */
    void PrincipalVector(Vector& rPrincipalVector, const Vector& rNonPrincipalVector);


      /**
     * @brief Checks for taunt/slack/wrinkles
     * @param rWrinklingState the current wrinkling state
     * @param rStress the stress
     * @param rStrain the strain
     */
    void CheckWrinklingState(WrinklingType& rWrinklingState, const Vector& rStress, const Vector& rStrain,
      Vector& rWrinklingDirectionVector);

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

    ConstitutiveLaw::Pointer mpConstitutiveLaw;


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
        rSerializer.save("ConstitutiveLaw", mpConstitutiveLaw);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw)
        rSerializer.load("ConstitutiveLaw", mpConstitutiveLaw);
    }


}; // Class WrinklingLinear2DLaw
}  // namespace Kratos.
#endif // KRATOS_WRINKLING_LINEAR_2D_LAW_H_INCLUDED  defined
