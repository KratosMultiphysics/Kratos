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

#if !defined (KRATOS_WRINKLING_2D_LAW_H_INCLUDED)
#define  KRATOS_WRINKLING_2D_LAW_H_INCLUDED

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
 * @class Wrinkling2DLaw
 * @ingroup StructuralMechanicsApplication
 * @brief This law defines a wrinkling modification for any 2D claw
 * @author Klaus B. Sautter
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) Wrinkling2DLaw
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

    /// Pointer definition of Wrinkling2DLaw
    KRATOS_CLASS_POINTER_DEFINITION( Wrinkling2DLaw );

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
    Wrinkling2DLaw();

    /**
     * @brief Copy constructor.
     */
    Wrinkling2DLaw (const Wrinkling2DLaw& rOther);

    /**
     * @brief Destructor.
     */
    ~Wrinkling2DLaw() override;

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


    template <class T>
    bool THas(const Variable<T>& rTemplateVariable) const;


    bool Has(const Variable<bool>& rThisVariable) override;


    bool Has(const Variable<int>& rThisVariable) override;


    bool Has(const Variable<double>& rThisVariable) override;


    bool Has(const Variable<Vector>& rThisVariable) override;


    bool Has(const Variable<Matrix>& rThisVariable) override;


    bool Has(const Variable<array_1d<double, 3 > >& rThisVariable) override;


    bool Has(const Variable<array_1d<double, 6 > >& rThisVariable) override;



    template <class T>
    T& TGetValue(const Variable<T>& rTemplateVariable, T& rTemplateValue);


    bool& GetValue(
        const Variable<bool>& rThisVariable,
        bool& rValue
        ) override;


    int& GetValue(
        const Variable<int>& rThisVariable,
        int& rValue
        ) override;


    double& GetValue(
        const Variable<double>& rThisVariable,
        double& rValue
        ) override;


    Vector& GetValue(
        const Variable<Vector>& rThisVariable,
        Vector& rValue
        ) override;


    Matrix& GetValue(
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue
        ) override;


    array_1d<double, 3 > & GetValue(
        const Variable<array_1d<double, 3 > >& rThisVariable,
        array_1d<double, 3 > & rValue
        ) override;


    array_1d<double, 6 > & GetValue(
        const Variable<array_1d<double, 6 > >& rThisVariable,
        array_1d<double, 6 > & rValue
        ) override;

    bool& CalculateValue(
        Parameters& rParameterValues,
        const Variable<bool>& rThisVariable,
        bool& rValue
        ) override;

    int& CalculateValue(
        Parameters& rParameterValues,
        const Variable<int>& rThisVariable,
        int& rValue
        ) override;

    double& CalculateValue(
        Parameters& rParameterValues,
        const Variable<double>& rThisVariable,
        double& rValue
        ) override;

    Vector& CalculateValue(
        Parameters& rParameterValues,
        const Variable<Vector>& rThisVariable,
        Vector& rValue
        ) override;

     Matrix& CalculateValue(
        Parameters& rParameterValues,
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue
        ) override;

     array_1d<double, 3 > & CalculateValue(
        Parameters& rParameterValues,
        const Variable<array_1d<double, 3 > >& rVariable,
        array_1d<double, 3 > & rValue
        ) override;

    array_1d<double, 6 > & CalculateValue(
        Parameters& rParameterValues,
        const Variable<array_1d<double, 6 > >& rVariable,
        array_1d<double, 6 > & rValue
        ) override;

    void SetValue(
        const Variable<bool>& rThisVariable,
        const bool& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    void SetValue(
        const Variable<int>& rThisVariable,
        const int& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

     void SetValue(
        const Variable<double>& rThisVariable,
        const double& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    void SetValue(
        const Variable<Vector >& rThisVariable,
        const Vector& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    void SetValue(
        const Variable<Matrix >& rThisVariable,
        const Matrix& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    void SetValue(
        const Variable<array_1d<double, 3 > >& rThisVariable,
        const array_1d<double, 3 > & rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    void SetValue(
        const Variable<array_1d<double, 6 > >& rThisVariable,
        const array_1d<double, 6 > & rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    void InitializeMaterial(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues
        ) override;

    void InitializeSolutionStep(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    void FinalizeSolutionStep(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;


    void InitializeNonLinearIteration(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;


    void FinalizeNonLinearIteration(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    void CalculateMaterialResponsePK2 (Parameters& rValues) override;


    void InitializeMaterialResponsePK2 (Parameters& rValues) override;


    void FinalizeMaterialResponsePK2 (Parameters& rValues) override;


    void ResetMaterial(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues
        ) override;


    void GetLawFeatures(Features& rFeatures) override;


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


}; // Class Wrinkling2DLaw
}  // namespace Kratos.
#endif // KRATOS_WRINKLING_2D_LAW_H_INCLUDED  defined
