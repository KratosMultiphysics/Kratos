// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:        BSD License
//                  license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//                   Sergio Jimenez
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "custom_constitutive/elastic_isotropic_3d.h"
#include "custom_constitutive/linear_plane_strain.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "constitutive_laws_application_variables.h"

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
 * @class AssociativePlasticDamageModel
 * @ingroup StructuralMechanicsApplication
 * @brief This law defines a parallel rule of mixture (classic law of mixture)
 * @details This constitutive law unifies the High cycle, Ultra Low cycle and Low cycle fatigue processes
 * by means of a plastic damage model. Source: A thermodynamically consistent plastic-damage framework for
localized failure in quasi-brittle solids: Material model and strain
localization analysis (Wu and Cervera https://doi.org/10.1016/j.ijsolstr.2016.03.005)
 * @author Alejandro Cornejo
 * @author Sergio Jimenez
 */
template<class TYieldSurfaceType>
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) AssociativePlasticDamageModel
    : public std::conditional<TYieldSurfaceType::VoigtSize == 6, ElasticIsotropic3D, LinearPlaneStrain >::type
{
public:

    ///@name Type Definitions
    ///@{
    /// The define the working dimension size, already defined in the integrator
    /// The definition of the size type
    typedef std::size_t SizeType;

    static constexpr SizeType Dimension = TYieldSurfaceType::Dimension;

    /// The define the Voigt size, already defined in the  integrator
    static constexpr SizeType VoigtSize = TYieldSurfaceType::VoigtSize;

    /// The definition of the process info
    typedef ProcessInfo ProcessInfoType;

    /// The definition of the CL base  class
    typedef typename std::conditional<VoigtSize == 6, ElasticIsotropic3D, LinearPlaneStrain >::type BaseType;

    /// Definition of the machine precision tolerance
    static constexpr double machine_tolerance = std::numeric_limits<double>::epsilon();

    static constexpr double tolerance = 1.0e-8;

    /// The node definition
    typedef Node NodeType;

    /// The geometry definition
    typedef Geometry<NodeType> GeometryType;

    /// The definition of the Voigt array type
    typedef array_1d<double, VoigtSize> BoundedVectorType;

    /// The definition of the bounded matrix type
    typedef BoundedMatrix<double, VoigtSize, VoigtSize> BoundedMatrixType;

    /// Pointer definition of AssociativePlasticDamageModel
    KRATOS_CLASS_POINTER_DEFINITION(AssociativePlasticDamageModel);

    struct PlasticDamageParameters {
        BoundedMatrixType ComplianceMatrixIncrement{ZeroMatrix(VoigtSize, VoigtSize)};
        BoundedMatrixType ComplianceMatrix{ZeroMatrix(VoigtSize, VoigtSize)};
        BoundedMatrixType ComplianceMatrixCompression{ZeroMatrix(VoigtSize, VoigtSize)};
        BoundedMatrixType ConstitutiveMatrix{ZeroMatrix(VoigtSize, VoigtSize)};
        BoundedMatrixType TangentTensor{ZeroMatrix(VoigtSize, VoigtSize)};
        BoundedVectorType PlasticFlow{ZeroVector(VoigtSize)};
        BoundedVectorType PlasticStrain{ZeroVector(VoigtSize)};
        BoundedVectorType PlasticStrainIncrement{ZeroVector(VoigtSize)};
        BoundedVectorType StrainVector{ZeroVector(VoigtSize)};
        BoundedVectorType StressVector{ZeroVector(VoigtSize)};
        double NonLinearIndicator          = 0.0; // F
        double PlasticConsistencyIncrement = 0.0; // Lambda dot
        double UniaxialStress              = 0.0;
        double DamageDissipation           = 0.0; // Kappa d
        double DamageDissipationIncrement  = 0.0; // Kappa d dot
        double PlasticDissipation          = 0.0; // Kappa p
        double PlasticDissipationIncrement = 0.0; // Kappa p dot
        double TotalDissipation            = 0.0; // Kappa
        double CharacteristicLength        = 0.0;
        double Threshold                   = 0.0;
        double Slope                       = 0.0; // d(Threshold)/d(dissipation)
        double PlasticDamageProportion     = 0.5; // 0-> Plastic    1->Damage
    };

    /// The definition of the lambdas to compute implicitly the threshold
    typedef std::function<double(const double, const double, ConstitutiveLaw::Parameters& , PlasticDamageParameters &)> ResidualFunctionType;

    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    AssociativePlasticDamageModel()
    {}

    /**
     * @brief Destructor.
     */
    ~AssociativePlasticDamageModel() override
    {}

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
    * Copy constructor.
    */
    AssociativePlasticDamageModel(const AssociativePlasticDamageModel &rOther)
        : BaseType(rOther),
          mPlasticDissipation(rOther.mPlasticDissipation),
          mDamageDissipation(rOther.mDamageDissipation),
          mThreshold(rOther.mThreshold),
          mPlasticStrain(rOther.mPlasticStrain),
          mOldStrain(rOther.mOldStrain),
          mComplianceMatrix(rOther.mComplianceMatrix),
          mComplianceMatrixCompression(rOther.mComplianceMatrixCompression)
    {
    }

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresInitializeMaterialResponse() override
    {
        return false;
    }

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresFinalizeMaterialResponse() override
    {
        return true;
    }

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
     * @brief Calculates the value of a specified variable (double)
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
     * @brief Initialize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponsePK2 (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Initialize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void InitializeMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues) override;

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


    /**
     * @brief This method computes the tangent tensor
     * @param rValues The constitutive law parameters and flags
     */
    void CalculateTangentTensor(
        ConstitutiveLaw::ConstitutiveLaw::Parameters& rValues,
        PlasticDamageParameters &rPlasticDamageParameters);


    /**
     * @brief This method computes the compliance elasatic matrix
     * (https://en.wikiversity.org/wiki/Introduction_to_Elasticity/Constitutive_relations)
     * @param rValues The constitutive law parameters and flags
     */
    void CalculateElasticComplianceMatrix(
        BoundedMatrixType& rConstitutiveMatrix,
        const Properties& rMaterialProperties
        );

    /**
     * @brief This method add somehting if the increment is positive
     */
    void AddIfPositive(
        double&  rToBeAdded,
        const double Increment
        )
    {
        if (Increment > machine_tolerance)
            rToBeAdded += Increment;
    }

    /**
     * @brief This method evaluates the Macaulay brackets
     */
    double MacaullyBrackets(const double Number)
    {
        return (Number > machine_tolerance) ? Number : 0.0;
    }

    /**
     * @brief This method computes the normalized
     * plastic dissipation
     */
    void CalculatePlasticDissipationIncrement(
        const Properties &rMaterialProperties,
        PlasticDamageParameters &rPlasticDamageParameters);

    /**
     * @brief This method computes the normalized
     * plastic dissipation
     */
    void CalculateDamageDissipationIncrement(
        const Properties &rMaterialProperties,
        PlasticDamageParameters &rPlasticDamageParameters);

    /**
     * @brief This method computes the threshold
     * and the dThreshold/dKappa according to energy dissipation
     */
    void CalculateThresholdAndSlope(
        ConstitutiveLaw::Parameters& rValues,
        PlasticDamageParameters &rPlasticDamageParameters);

    /**
     * @brief This method computes the plastic flow (dF/dS)
     */
    void CalculateFlowVector(
        ConstitutiveLaw::Parameters& rValues,
        PlasticDamageParameters &rPlasticDamageParameters);

    /**
     * @brief This method computes the Plastic Strain increment
     */
    void CalculatePlasticStrainIncrement(
        ConstitutiveLaw::Parameters& rValues,
        PlasticDamageParameters &rPlasticDamageParameters);

    /**
     * @brief This method computes the Compliance matrix increment
     */
    void CalculateComplianceMatrixIncrement(
        ConstitutiveLaw::Parameters& rValues,
        PlasticDamageParameters &rPlasticDamageParameters);

    /**
     * @brief This method computes the plastic consistency increment
     */
    void CalculatePlasticConsistencyIncrement(
        ConstitutiveLaw::Parameters& rValues,
        PlasticDamageParameters &rPlasticDamageParameters);

    /**
     * @brief This method integrates the stress
     */
    void IntegrateStressPlasticDamageMechanics(
        ConstitutiveLaw::Parameters& rValues,
        PlasticDamageParameters &rPlasticDamageParameters);

    /**
     * @brief This method computes the constitutive matrix
     */
    void CalculateConstitutiveMatrix(
        ConstitutiveLaw::Parameters& rValues,
        PlasticDamageParameters &rPlasticDamageParameters);

    /**
     * @brief This method updates all the internal variables
     */
    void UpdateInternalVariables(
        PlasticDamageParameters &rPlasticDamageParameters
    );

    /**
     * @brief This method updates all the internal variables
     */
    void CheckMinimumFractureEnergy(
        ConstitutiveLaw::Parameters& rValues,
        PlasticDamageParameters &rPDParameters
    );

    /**
     * @brief This method initializes all the values
     * in the PlasticDamageParameters
     */
    void InitializePlasticDamageParameters(
        const BoundedVectorType& rStrainVector,
        const Properties& rMaterialProperties,
        const double CharateristicLength,
        PlasticDamageParameters &rPlasticDamageParameters
        )
    {
        rPlasticDamageParameters.PlasticDissipation     = mPlasticDissipation;
        rPlasticDamageParameters.DamageDissipation      = mDamageDissipation;
        rPlasticDamageParameters.TotalDissipation       = mPlasticDissipation + mDamageDissipation;
        rPlasticDamageParameters.Threshold              = mThreshold;
        noalias(rPlasticDamageParameters.PlasticStrain) = mPlasticStrain;
        noalias(rPlasticDamageParameters.ComplianceMatrix) = mComplianceMatrix;
        noalias(rPlasticDamageParameters.ComplianceMatrixCompression) = mComplianceMatrixCompression;
        noalias(rPlasticDamageParameters.StrainVector) = rStrainVector;
        rPlasticDamageParameters.CharacteristicLength  = CharateristicLength;
        rPlasticDamageParameters.PlasticDamageProportion = rMaterialProperties[PLASTIC_DAMAGE_PROPORTION];
    }

    /**
     * @brief This method computes the continuum
     * analytical tangent tensor
     */
    void CalculateAnalyticalTangentTensor(
        ConstitutiveLaw::Parameters& rValues,
        PlasticDamageParameters &rParam
        );

    /**
     * @brief This method increases the damage and plastic
     * dissipation with the increment
     */
    void AddNonLinearDissipation(
        PlasticDamageParameters &rPDParameters
        )
    {
        rPDParameters.DamageDissipation  += rPDParameters.DamageDissipationIncrement;
        rPDParameters.DamageDissipation = (rPDParameters.DamageDissipation > 0.99999) ?
            0.99999 : rPDParameters.DamageDissipation;

        rPDParameters.PlasticDissipation += rPDParameters.PlasticDissipationIncrement;
        rPDParameters.PlasticDissipation = (rPDParameters.PlasticDissipation > 0.99999) ?
            0.99999 : rPDParameters.PlasticDissipation;

        rPDParameters.TotalDissipation = (rPDParameters.PlasticDissipation +
            rPDParameters.DamageDissipation);
        rPDParameters.TotalDissipation = (rPDParameters.TotalDissipation > 0.99999) ?
            0.99999 : rPDParameters.TotalDissipation;
    }

    /**
     * @brief This method computes an averaged
     * volumetric fracture energy depending
     * if it is in tension or compression
     */
    static double CalculateVolumetricFractureEnergy( // g_F
        const Properties& rMaterialProperties,
        PlasticDamageParameters &rPDParameters
        );

    /**
     * @brief This method computes the denominator
     * of the expression for computing the
     * plastic multiplier
     */
    double CalculatePlasticDenominator(
        ConstitutiveLaw::Parameters& rValues,
        PlasticDamageParameters &rParam);

    /**
     * @brief This method solves a non-linear
     * equation that related the dissipation
     * with the threshold
     */
    double CalculateThresholdImplicitExpression(
        ResidualFunctionType &rF,
        ResidualFunctionType &rdF_dk,
        ConstitutiveLaw::Parameters &rValues,
        PlasticDamageParameters &rPDParameters,
        const double MaxThreshold = std::numeric_limits<double>::max());

    /**
     * @brief This method computes the slope or
     * d(threshold)/d(dissipation) by finite
     * differences
     */
    double CalculateSlopeFiniteDifferences(
        ResidualFunctionType &rF,
        ResidualFunctionType &rdF_dk,
        ConstitutiveLaw::Parameters &rValues,
        PlasticDamageParameters &rPDParameters,
        const double MaxThreshold = std::numeric_limits<double>::max());

    /**
     * @brief Implicit function that relates the
     * plastic-damage energy dissipation with the
     * uniaxial stress threshold
     */
    ResidualFunctionType ExponentialSofteningImplicitFunction();

    /**
     * @brief Implicit function derivative to be used
     * in the minimization of the implicit function
     */
    ResidualFunctionType ExponentialSofteningImplicitFunctionDerivative();

    /**
     * @brief Implicit function that relates the
     * plastic-damage energy dissipation with the
     * uniaxial stress threshold
     */
    ResidualFunctionType ExponentialHardeningImplicitFunction();

    /**
     * @brief Implicit function derivative to be used
     * in the minimization of the implicit function
     */
    ResidualFunctionType ExponentialHardeningImplicitFunctionDerivative();

    /**
     * @brief Implicit function that relates the
     * plastic-damage energy dissipation with the
     * uniaxial stress threshold
     */
    ResidualFunctionType CurveByPointsHardeningImplicitFunction();

    /**
     * @brief Implicit function derivative to be used
     * in the minimization of the implicit function
     */
    ResidualFunctionType CurveByPointsHardeningImplicitFunctionDerivative();
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

    double mPlasticDissipation = 0.0;
    double mDamageDissipation  = 0.0;
    double mThreshold          = 0.0;
    BoundedVectorType mPlasticStrain    = ZeroVector(VoigtSize);
    BoundedVectorType mOldStrain        = ZeroVector(VoigtSize);
    BoundedMatrixType mComplianceMatrix = ZeroMatrix(VoigtSize, VoigtSize);
    BoundedMatrixType mComplianceMatrixCompression = ZeroMatrix(VoigtSize, VoigtSize);

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

    void save(Serializer &rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.save("PlasticDissipation", mPlasticDissipation);
        rSerializer.save("DamageDissipation", mDamageDissipation);
        rSerializer.save("Threshold", mThreshold);
        rSerializer.save("PlasticStrain", mPlasticStrain);
        rSerializer.save("OldStrain", mOldStrain);
        rSerializer.save("ComplianceMatrix", mComplianceMatrix);
        rSerializer.save("ComplianceMatrixCompression", mComplianceMatrixCompression);
    }

    void load(Serializer &rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.load("PlasticDissipation", mPlasticDissipation);
        rSerializer.load("DamageDissipation", mDamageDissipation);
        rSerializer.load("Threshold", mThreshold);
        rSerializer.load("PlasticStrain", mPlasticStrain);
        rSerializer.load("OldStrain", mOldStrain);
        rSerializer.load("ComplianceMatrix", mComplianceMatrix);
        rSerializer.load("ComplianceMatrixCompression", mComplianceMatrixCompression);
    }


}; // Class AssociativePlasticDamageModel
}  // namespace Kratos.
