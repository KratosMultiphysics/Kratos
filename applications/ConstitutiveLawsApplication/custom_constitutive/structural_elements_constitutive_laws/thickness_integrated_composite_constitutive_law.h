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

#include "custom_utilities/advanced_constitutive_law_utilities.h"

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
 * @class ThicknessIntegratedCompositeConstitutiveLaw
 * @ingroup StructuralMechanicsApplication
 * @brief This constitutive law is implemented to be used together with the CSDSG3ThickShellElement3D3N or MITCThickShellElement3D4N element
 * It computes the constitutive matrix and stress vector by integrating over the thickness using a set of
 * sub-properties, each one with its own subproperty and constitutive model (isotropic plane stress).
 * For each layer, the user needs to provide: thickness, z coordinate of the layer w.r.t the bending axis and the Euler angle of rotation
 * w.r.t the local axes of the shell.
 * Input: the Generalized strain vector of 8 components (3 membrane + 3 bending + 2 shear)
 * Output: the Generalized stress vector of 8 components (3 membrane + 3 bending + 2 shear) and 8x8 generalized
 * constitutive matrix.
 * The number of integration points can be defined by the user, by default 5 points are used.
 * @author Alejandro Cornejo
 */
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) ThicknessIntegratedCompositeConstitutiveLaw
    : public ConstitutiveLaw
{
public:

    ///@name Type Definitions
    ///@{

    /// The definition of the process info
    using ProcessInfoType = ProcessInfo;

    /// The definition of the CL base class
    using BaseType = ConstitutiveLaw;

    /// Unhide the base class SetValue overloads
    using BaseType::SetValue;

    /// Unhide the base class GetValue overloads
    using BaseType::GetValue;

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

    /// Pointer definition of ThicknessIntegratedCompositeConstitutiveLaw
    KRATOS_CLASS_POINTER_DEFINITION(ThicknessIntegratedCompositeConstitutiveLaw);

    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    ThicknessIntegratedCompositeConstitutiveLaw();

    /**
     * @brief Constructor with values
     * @param rThicknessIntegrationPoints The amount of thickness integration points
     */
    ThicknessIntegratedCompositeConstitutiveLaw(
        const std::vector<double>& rZCoordinates,
        const std::vector<double>& rEulerAngles,
        const std::vector<double>& rThicknesses);

    /**
     * @brief Copy constructor.
     */
    ThicknessIntegratedCompositeConstitutiveLaw(const ThicknessIntegratedCompositeConstitutiveLaw &rOther);

    /**
     * @brief Destructor.
     */
    ~ThicknessIntegratedCompositeConstitutiveLaw() override;

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
        for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
            if (mConstitutiveLaws[i_layer]->RequiresInitializeMaterialResponse())
                return true;
        }
        return false;
    }

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresFinalizeMaterialResponse() override
    {
        for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
            if (mConstitutiveLaws[i_layer]->RequiresFinalizeMaterialResponse())
                return true;
        }
        return false;
    }

    /**
     * @brief Returns whether this constitutive Law has specified variable (generic)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    template<class TDataType>
    bool THas(const Variable<TDataType>& rThisVariable)
    {
        for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
            if (mConstitutiveLaws[i_layer]->Has(rThisVariable))
                return true;
        }
        return false;
    }

    /**
     * @brief Returns whether this constitutive Law has specified variable (generic)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    template<class TDataType>
    TDataType& TGetValue(
        const Variable<TDataType>& rThisVariable,
        TDataType& rValue
        )
    {
        KRATOS_TRY

        // IndexType number_of_laws = mConstitutiveLaws.size();
        // TDataType ip_value;
        // if (mConstitutiveLaws[0]->Has(rThisVariable)) {
        //     mConstitutiveLaws[0]->GetValue(rThisVariable, rValue);

        //     for (IndexType i = 1; i < number_of_laws; ++i) {
        //         mConstitutiveLaws[i]->GetValue(rThisVariable, ip_value);
        //         rValue += ip_value;
        //     }
        //     rValue /= static_cast<double>(number_of_laws);
        // }
        return rValue;

        KRATOS_CATCH("Generic GetValue")
    }

    /**
     * @brief Sets the value of a specified variable (generic)
     * @param rThisVariable the variable to be checked for
     */
    template<class TDataType>
    void TSetValue(
        const Variable<TDataType>& rThisVariable,
        const TDataType& rValue,
        const ProcessInfo& rCurrentProcessInfo
        )
    {
        KRATOS_TRY

        // if (mConstitutiveLaws[0]->Has(rThisVariable)) {
        //     for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
        //         mConstitutiveLaws[i]->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
        //     }
        // }

        KRATOS_CATCH("Generic SetValue")
    }

    /**
     * @brief CalculateValue of a specified variable (generic)
     * @param rThisVariable the variable to be checked for
     */
    template<class TDataType>
    TDataType& TCalculateValue(
        Parameters& rParameterValues,
        const Variable<TDataType>& rThisVariable,
        TDataType& rValue)
        {
            KRATOS_TRY

            // const Properties& r_material_properties = rParameterValues.GetMaterialProperties();
            // const auto sub_property = r_material_properties.GetSubProperties().begin();
            // IndexType number_of_laws = mConstitutiveLaws.size();

            // TDataType aux_value;

            // rParameterValues.SetMaterialProperties(*(sub_property));
            // mConstitutiveLaws[0]->CalculateValue(rParameterValues, rThisVariable, rValue);

            // for (IndexType i = 1; i < number_of_laws; ++i) {
            //     mConstitutiveLaws[i]->CalculateValue(rParameterValues, rThisVariable, aux_value);
            //     rValue += aux_value;
            // }
            // rValue /= static_cast<double>(number_of_laws);

            // rParameterValues.SetMaterialProperties(r_material_properties);

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

    /**
     * @brief Computes the material maximum edge length for the given geometry, used to compute the integration points in thickness
     */
    double GetMaxReferenceEdgeLength(const GeometryType& rGeometry) const
    {
        return AdvancedConstitutiveLawUtilities<3>::GetMaxReferenceEdgeLengthForShell(rGeometry);
    }

    /**
     * @brief Computes the shear modulus of the shell
     */
    void CalculateShearModulus(double &Gyz, double &Gxz, Parameters &rValues);

    /**
     * @brief Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK1(Parameters& rValues) override;

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
    std::vector<double> mZCoordinates; /// The vector containing the z-coordinate of the centroid of each layer in the thickness of the shell
    std::vector<double> mEulerAngles; /// The Euler angle of rotation of each layer w.r.t the local axes of the shell
    std::vector<double> mThicknesses; /// The thickness of each layer

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
        rSerializer.save("ZCoordinates", mZCoordinates);
        rSerializer.save("EulerAngles", mEulerAngles);
        rSerializer.save("Thicknesses", mThicknesses);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw)
        rSerializer.load("ConstitutiveLaws", mConstitutiveLaws);
        rSerializer.load("ThicknessIntegrationPoints", mZCoordinates);
        rSerializer.load("EulerAngles", mEulerAngles);
        rSerializer.load("Thicknesses", mThicknesses);
    }


}; // Class ThicknessIntegratedCompositeConstitutiveLaw
}  // namespace Kratos.
