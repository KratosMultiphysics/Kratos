// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#pragma once

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
 * @class ThicknessIntegratedIsotropicConstitutiveLaw
 * @ingroup StructuralMechanicsApplication
 * @brief This constitutive law is implemented to be used together with the CSDSG3ThickShellElement3D3N element
 * It computes the constitutive matrix and stress vector by integrating over the thickness using a set of
 * sub-properties, each one with its own subproperty, constitutive model (isotropic plane stress).
 * Input: the Generalized strain vector of 8 components (3 membrane + 3 bending + 2 shear)
 * Output: the Generalized stress vector of 8 components (3 membrane + 3 bending + 2 shear) and 8x8 generalized
 * constitutive matrix.
 * The number of integration points can be defined by the user, by default 5 points are used.
 * @author Alejandro Cornejo
 */
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) ThicknessIntegratedIsotropicConstitutiveLaw
    : public ConstitutiveLaw
{
public:

    ///@name Type Definitions
    ///@{

    /// The definition of the process info
    using ProcessInfoType = ProcessInfo;

    /// The definition of the CL base class
    using BaseType = ConstitutiveLaw;

    /// The definition of the size type
    using SizeType = std::size_t;

    /// The definition of the index type
    using IndexType = std::size_t;

    /// The define the working dimension size, generalized strains/stresses size
    static constexpr SizeType VoigtSize = 8;

    /// The Dimension of the CL
    static constexpr SizeType Dimension = 3;

    /// Definition of the machine precision tolerance
    static constexpr double machine_tolerance = std::numeric_limits<double>::epsilon();

    /// Pointer definition of ThicknessIntegratedIsotropicConstitutiveLaw
    KRATOS_CLASS_POINTER_DEFINITION(ThicknessIntegratedIsotropicConstitutiveLaw);

    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    ThicknessIntegratedIsotropicConstitutiveLaw();

    /**
     * @brief Constructor with values
     * @param rThicknessIntegrationPoints The amount of thickness integration points
     */
    ThicknessIntegratedIsotropicConstitutiveLaw(const IndexType rThicknessIntegrationPoints);

    /**
     * @brief Copy constructor.
     */
    ThicknessIntegratedIsotropicConstitutiveLaw(const ThicknessIntegratedIsotropicConstitutiveLaw &rOther);

    /**
     * @brief Destructor.
     */
    ~ThicknessIntegratedIsotropicConstitutiveLaw() override;

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
    SizeType GetStrainSize() const override;

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresInitializeMaterialResponse() override
    {
        return mConstitutiveLaws[0]->RequiresInitializeMaterialResponse(); // All the same;
    }

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresFinalizeMaterialResponse() override
    {
        return mConstitutiveLaws[0]->RequiresFinalizeMaterialResponse(); // All the same;
    }

    /**
     * @brief Returns whether this constitutive Law has specified variable (generic)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    template<class TDataType>
    bool Has(const Variable<TDataType>& rThisVariable)
    {
        return mConstitutiveLaws[0]->Has(rThisVariable); // All the same
    }

    /**
     * @brief Returns whether this constitutive Law has specified variable (generic)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    template<class TDataType>
    TDataType& GetValue(
        const Variable<TDataType>& rThisVariable,
        TDataType& rValue
        )
    {
        KRATOS_TRY

        IndexType number_of_laws = mConstitutiveLaws.size();
        TDataType ip_value;
        if (mConstitutiveLaws[0]->Has(rThisVariable)) {
            mConstitutiveLaws[0]->GetValue(rThisVariable, rValue);

            for (IndexType i = 1; i < number_of_laws; ++i) {
                mConstitutiveLaws[i]->GetValue(rThisVariable, ip_value);
                rValue += ip_value;
            }
            rValue /= static_cast<double>(number_of_laws);
        }
        return rValue;

        KRATOS_CATCH("Generic GetValue")
    }

    /**
     * @brief Sets the value of a specified variable (generic)
     * @param rThisVariable the variable to be checked for
     */
    template<class TDataType>
    void SetValue(
        const Variable<TDataType>& rThisVariable,
        const TDataType& rValue,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY

        if (mConstitutiveLaws[0]->Has(rThisVariable)) {
            for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
                mConstitutiveLaws[i]->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
            }
        }

        KRATOS_CATCH("Generic SetValue")
    }

    /**
     * @brief CalculateValue of a specified variable (generic)
     * @param rThisVariable the variable to be checked for
     */
    template<class TDataType>
    TDataType& CalculateValue(
        Parameters& rParameterValues,
        const Variable<TDataType>& rThisVariable,
        TDataType& rValue)
        {
            KRATOS_TRY

            const Properties& r_material_properties = rParameterValues.GetMaterialProperties();
            const auto sub_property = r_material_properties.GetSubProperties().begin();
            IndexType number_of_laws = mConstitutiveLaws.size();

            TDataType aux_value;

            rParameterValues.SetMaterialProperties(*(sub_property));
            mConstitutiveLaws[0]->CalculateValue(rParameterValues, rThisVariable, rValue);

            for (IndexType i = 1; i < number_of_laws; ++i) {
                mConstitutiveLaws[i]->CalculateValue(rParameterValues, rThisVariable, aux_value);
                rValue += aux_value;
            }
            rValue /= static_cast<double>(number_of_laws);

            rParameterValues.SetMaterialProperties(r_material_properties);
            
            return rValue;

            KRATOS_CATCH("Generic CalculateValue")
        }

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
        const Variable<Vector>& rThisVariable,
        const Vector& rValue,
        const ProcessInfo& rCurrentProcessInfo
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


    std::vector<ConstitutiveLaw::Pointer>& GetConstitutiveLaws()
    {
        return mConstitutiveLaws;
    }

    IndexType& GetNumberOfThicknessIntegrationPoints()
    {
        return mThicknessIntegrationPoints;
    }

    void SetNumberOfThicknessIntegrationPoints(const IndexType NumberOfPoints)
    {
        mThicknessIntegrationPoints = NumberOfPoints;
    }

    double GetMaxReferenceEdgeLength(const GeometryType& rGeometry) const
    {
        double max_length = 0.0;

        const auto& r_coord_1 = rGeometry[0].GetInitialPosition();
        const auto& r_coord_2 = rGeometry[1].GetInitialPosition();
        const auto& r_coord_3 = rGeometry[2].GetInitialPosition();

        const double length_12 = norm_2(r_coord_2 - r_coord_1);
        const double length_23 = norm_2(r_coord_3 - r_coord_2);
        const double length_31 = norm_2(r_coord_1 - r_coord_3);

        max_length = std::max(length_12, length_23);
        max_length = std::max(max_length, length_31);

        return max_length;
    }

    void CalculateCoordinatesAndWeights(
        std::vector<double> &rCoordinates,
        std::vector<double> &rWeights,
        const IndexType NumberOfPoints,
        const Properties& rMaterialProperties);

    /**
     * @brief Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK1 (Parameters& rValues) override;

    /**
     * @brief Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK2(Parameters& rValues) override;

    /**
     * @brief Computes the material response in terms of Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseKirchhoff(Parameters& rValues) override;

    /**
     * @brief Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy(Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponsePK1(Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponsePK2(Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponseKirchhoff(Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void InitializeMaterialResponseCauchy(Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK1(Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK2(Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseKirchhoff(Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy(Parameters& rValues) override;

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
        ) const override;




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
    IndexType mThicknessIntegrationPoints = 5;               /// The number of thickness integration points

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
        rSerializer.save("ThicknessIntegrationPoints", mThicknessIntegrationPoints);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw)
        rSerializer.load("ConstitutiveLaws", mConstitutiveLaws);
        rSerializer.load("ThicknessIntegrationPoints", mThicknessIntegrationPoints);
    }


}; // Class ThicknessIntegratedIsotropicConstitutiveLaw
}  // namespace Kratos.
