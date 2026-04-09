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
#include "custom_utilities/constitutive_law_utilities.h"

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

        TDataType ip_value;

            for (IndexType i = 1; i < mConstitutiveLaws.size(); ++i) {
                if (mConstitutiveLaws[i]->Has(rThisVariable)) {
                    mConstitutiveLaws[i]->GetValue(rThisVariable, ip_value);
                    rValue += ip_value * mThicknesses[i]; // We integrate the value through the thickness, multiplying by the layer thickness
                }
            }

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

        for (IndexType i = 0; i < mConstitutiveLaws.size(); ++i) {
            if (mConstitutiveLaws[0]->Has(rThisVariable))
                mConstitutiveLaws[i]->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
        }

        KRATOS_CATCH("Generic SetValue")
    }

    /**
     * @brief CalculateValue of a specified variable (generic)
     * @param rThisVariable the variable to be checked for
     */
    template<class TDataType>
    TDataType& TCalculateValue(
        Parameters& rValues,
        const Variable<TDataType>& rThisVariable,
        TDataType& rValue)
        {
            KRATOS_TRY

            const IndexType number_of_laws = mConstitutiveLaws.size();
            const auto subprop_strain_size = mConstitutiveLaws[0]->GetStrainSize(); // 3
            const auto subprop_dimension = mConstitutiveLaws[0]->WorkingSpaceDimension(); // 2
            const auto& r_material_properties = rValues.GetMaterialProperties();

            const Vector generalized_strain_vector = rValues.GetStrainVector(); // size 8

            Vector& r_stress_vector = rValues.GetStressVector(); // size 3
            Vector& r_strain_vector = rValues.GetStrainVector(); // size 3
            Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix(); // size 3x3
            r_strain_vector.resize(subprop_strain_size, false);
            r_stress_vector.resize(subprop_strain_size, false);
            r_constitutive_matrix.resize(subprop_strain_size, subprop_strain_size, false);
            r_strain_vector.clear();
            r_stress_vector.clear();
            r_constitutive_matrix.clear();
            
            Matrix F(subprop_dimension, subprop_dimension); // 2x2
            double weight, z_coord, aux_weight, detF, Euler_angle;
            
            double Gyz = 0.0;
            double Gxz = 0.0;
            
            const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
            // Rotation and strain-rotation matrices
            BoundedMatrix<double, 2, 2> T;
            BoundedMatrix<double, 3, 3> Tvoigt;

            TDataType ip_value;

            // We perform the integration through the thickness
            for (IndexType i_layer = 0; i_layer < number_of_laws; ++i_layer) {
                // Assign subprops of the layer
                Properties &r_subprop = *(it_prop_begin + i_layer);
                rValues.SetMaterialProperties(r_subprop);

                CalculateShearModulus(Gyz, Gxz, rValues);

                // We retrieve the layer information
                z_coord = mZCoordinates[i_layer];
                Euler_angle = mEulerAngles[i_layer];

                r_strain_vector[0] = generalized_strain_vector[0] + z_coord * generalized_strain_vector[3]; // xx
                r_strain_vector[1] = generalized_strain_vector[1] + z_coord * generalized_strain_vector[4]; // yy
                r_strain_vector[2] = generalized_strain_vector[2] + z_coord * generalized_strain_vector[5]; // xy

                // We rotate the strain to layer local axes
                AdvancedConstitutiveLawUtilities<3>::CalculateRotationOperatorEuler1(Euler_angle, T);
                ConstitutiveLawUtilities<3>::CalculateRotationOperatorVoigt(T, Tvoigt);
                r_strain_vector = prod(Tvoigt, r_strain_vector);

                // In case the 2D Cls work in finite strain
                noalias(F) = AdvancedConstitutiveLawUtilities<3>::ComputeEquivalentSmallDeformationDeformationGradient(r_strain_vector);
                detF = MathUtils<double>::Det2(F);
                rValues.SetDeterminantF(detF);
                rValues.SetDeformationGradientF(F);

                mConstitutiveLaws[i_layer]->CalculateValue(rValues, rThisVariable, ip_value);
                rValue += ip_value * mThicknesses[i_layer]; // We integrate the value through the thickness, multiplying by the layer thickness
            }

            r_strain_vector.resize(VoigtSize, false);
            noalias(r_strain_vector) = generalized_strain_vector;

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

    /**
     * @brief This is to calculate and initialize the mShearReductionFactors at the very beginning of the calculation
     * @param rMaterialProperties the Properties instance of the current element
     * This implementation is based on E. Oñate Vol 2: Beams plates and shells. Page 399 of the document
     * Eq. (7.48a) from section 7.3  "computation of transverse shear correction parameters".
     */
    void InitializeShearReductionFactors(
        const Properties &rMaterialProperties);

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
    std::vector<double> mShearReductionFactors = std::vector<double>(2, 1.0); /// The shear reduction factors for the shear components of the constitutive matrix, initialized to 1.0 (no reduction)

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.save("ConstitutiveLaws", mConstitutiveLaws);
        rSerializer.save("ZCoordinates", mZCoordinates);
        rSerializer.save("EulerAngles", mEulerAngles);
        rSerializer.save("Thicknesses", mThicknesses);
        rSerializer.save("ShearReductionFactors", mShearReductionFactors);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.load("ConstitutiveLaws", mConstitutiveLaws);
        rSerializer.load("ThicknessIntegrationPoints", mZCoordinates);
        rSerializer.load("EulerAngles", mEulerAngles);
        rSerializer.load("Thicknesses", mThicknesses);
        rSerializer.load("ShearReductionFactors", mShearReductionFactors);
    }


}; // Class ThicknessIntegratedCompositeConstitutiveLaw
}  // namespace Kratos.
