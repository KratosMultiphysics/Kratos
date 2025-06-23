// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:                 BSD License
//                               license: structural_mechanics_application/license.txt
//
//  Main authors:    Alireza Taherzadeh-Fard
//                   Alejandro Cornejo Velazquez
//                   Sergio Jimenez Reyes
//                   Lucia Gratiela Barbu
//

#pragma once

// System includes

// External includes

// Project includes
#include "custom_constitutive/composites/rule_of_mixtures_law.h"


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
 * @class TractionSeparationLaw3D
 * @ingroup StructuralMechanicsApplication
 * @brief This law defines a parallel rule of mixture (classic law of mixture)
 * @details The constitutive law show have defined a subproperties in order to work properly
 * This law combines parallel CL considering the following principles:
 *  - All layer have the same strain
 *  - The total stress is the addition of the strain in each layer
 *  - The constitutive tensor is the addition of the constitutive tensor of each layer
 * @author Vicente Mataix Ferrandiz
 * @author Fernando Rastellini
 * @author Alejandro Cornejo Velazquez
 */
template<unsigned int TDim>
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) TractionSeparationLaw3D
    : public ParallelRuleOfMixturesLaw<TDim>
{
public:

    ///@name Type Definitions
    ///@{
    typedef ParallelRuleOfMixturesLaw<TDim> BaseType;

    typedef std::size_t SizeType;

    /// The define the working dimension size, already defined in the integrator
    static constexpr SizeType VoigtSize = (TDim == 3) ? 6 : 3;

    /// The define the Voigt size, already defined in the  integrator
    static constexpr SizeType Dimension = TDim;

    /// Definition of the machine precision tolerance
    static constexpr double machine_tolerance = std::numeric_limits<double>::epsilon();

    /// Pointer definition of TractionSeparationLaw3D
    KRATOS_CLASS_POINTER_DEFINITION(TractionSeparationLaw3D);

    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    TractionSeparationLaw3D();

    /**
     * @brief Constructor with values
     * @param rCombinationFactors The list of subproperties combination factors
     */
    TractionSeparationLaw3D(const std::vector<double> &rCombinationFactors);

    /**
     * @brief Copy constructor.
     */
    TractionSeparationLaw3D(const TractionSeparationLaw3D &rOther);

    /**
     * @brief Destructor.
     */
    ~TractionSeparationLaw3D() override;

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
     * @brief Returns whether this constitutive Law has specified variable (Vector)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Vector>& rThisVariable) override;

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
      * @brief Is called to check whether the provided material parameters in the Properties match the requirements of current constitutive model.
      * @param rMaterialProperties the current Properties to be validated against.
      * @return true, if parameters are correct; false, if parameters are insufficient / faulty
      * @note  this has to be implemented by each constitutive model. Returns false in base class since no valid implementation is contained here.
      */
    bool ValidateInput(const Properties& rMaterialProperties) override;

    /**
     * @brief This is to be called at the very beginning of the calculation
     * @details (e.g. from InitializeElement) in order to initialize all relevant attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void InitializeMaterial(
        const Properties& rMaterialProperties,
        const ConstitutiveLaw::GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues
        ) override;

    /**
     * @brief Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK2 (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Computes the material response in terms of Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK2 (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Finalize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Calculate delamination damage in different modes
     * @see Parameters
     */
    double CalculateDelaminationDamageExponentialSoftening (
        ConstitutiveLaw::Parameters& rValues,
        const double GI,
        const double E,
        const double T0,
        const double equivalent_stress);

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
        const ConstitutiveLaw::GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    /**
     * @brief This method computes the tangent tensor
     * @param rValues The constitutive law parameters and flags
     */
    void CalculateTangentTensor(
        ConstitutiveLaw::Parameters& rValues,
        const ConstitutiveLaw::StressMeasure& rStressMeasure
    );

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

    Vector mDelaminationDamageModeOne;
    Vector mDelaminationDamageModeTwo;
    Vector mThresholdModeOne;
    Vector mThresholdModeTwo;

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ParallelRuleOfMixturesLaw<TDim> )
        rSerializer.save("DelaminationDamageModeOne", mDelaminationDamageModeOne);
        rSerializer.save("DelaminationDamageModeTwo", mDelaminationDamageModeTwo);
        rSerializer.save("ThresholdModeOne", mThresholdModeOne);
        rSerializer.save("ThresholdModeTwo", mThresholdModeTwo);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ParallelRuleOfMixturesLaw<TDim>)
        rSerializer.load("DelaminationDamageModeOne", mDelaminationDamageModeOne);
        rSerializer.load("DelaminationDamageModeTwo", mDelaminationDamageModeTwo);
        rSerializer.load("ThresholdModeOne", mThresholdModeOne);
        rSerializer.load("ThresholdModeTwo", mThresholdModeTwo);
    }

}; // Class TractionSeparationLaw3D
}  // namespace Kratos.
