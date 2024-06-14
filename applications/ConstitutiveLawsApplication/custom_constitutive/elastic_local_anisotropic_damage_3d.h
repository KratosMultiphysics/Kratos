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
//  Main authors:   athira vadakkekkara
//  Collaborator:

#pragma once

// System includes

// External includes

// Project includes
#include "custom_constitutive/elastic_isotropic_3d.h"

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
 * @class ElasticAnisotropicDamage
 * @ingroup StructuralMechanicsApplication
 * @brief Defines a damage with hardening constitutive law in 3D
 * @details This material law is defined by the parameters:
 * - YOUNG_MODULUS
 * - POISSON_RATIO
 * - HARDEDNING_CURVE: Type of hardening model: 0 exponential, 1 multilinear
 * - STRESS_LIMITS: list of stress values in which the corresponding hardening
 *   parameter is valid
 * - HARDENING_PARAMETERS: List of hardening modules (max three branches considered)
 * @warning Valid for small strains
 * @note
 * @author Marcelo Raschi
 */
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) ElasticAnisotropicDamage
    : public ElasticIsotropic3D
{
public:

    ///@name Type Definitions
    ///@{
    typedef ProcessInfo ProcessInfoType;
    typedef ConstitutiveLaw BaseType;
    typedef std::size_t SizeType;
    /// Static definition of the dimension
    static constexpr SizeType Dimension = 3;

    /// Static definition of the VoigtSize
    static constexpr SizeType VoigtSize = 6;

    // The definition of the bounded vector type
    typedef BoundedVector<double, Dimension > BoundedVectorType;

    // The definition of the bounded vector type
    typedef BoundedVector<double, VoigtSize > BoundedVectorVoigtType;

    /// The definition of the bounded matrix type
    typedef BoundedMatrix<double, Dimension, Dimension> BoundedMatrixType;

    /// The definition of the bounded matrix type
    typedef BoundedMatrix<double, VoigtSize, VoigtSize> BoundedMatrixVoigtType;

    /// The definition of the bounded matrix type
    typedef BoundedMatrix<double, VoigtSize, Dimension> BoundedMatrix6x3Type;

    /// The definition of the bounded matrix type
    typedef BoundedMatrix<double, Dimension, VoigtSize> BoundedMatrix3x6Type;

    // Counted pointer
    KRATOS_CLASS_POINTER_DEFINITION(ElasticAnisotropicDamage);

    ///@}
    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    ElasticAnisotropicDamage();

    /**
     * @brief Copy constructor.
     */
    ElasticAnisotropicDamage(const ElasticAnisotropicDamage& rOther);

    /**
     * @brief Destructor.
     */
    ~ElasticAnisotropicDamage() override;

    /**
     * @brief Clone function
     * @return A pointer to a new instance of this constitutive law
     */
    ConstitutiveLaw::Pointer Clone() const override;

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
     * @brief Returns the value of a specified variable (Vector)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Vector& GetValue(
        const Variable<Vector>& rThisVariable,
        Vector& rValues
        ) override;

    /**
     * @brief Sets the value of a specified variable (Vector)
     * @param rThisVariable The variable to be returned
     * @param rValue New value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<Vector>& rThisVariable,
        const Vector& rValues,
        const ProcessInfo& rProcessInfo
        ) override;

    /**
     * @brief This is to be called at the very beginning of the calculation
     * @details (e.g. from InitializeElement) in order to initialize all relevant attributes of the constitutive law
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     */
    void InitializeMaterial(const Properties& rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const Vector& rShapeFunctionsValues) override;

    /**
     * @brief This method computes the stress and constitutive tensor
     * @param rValues The norm of the deviation stress
     * @param rStrainVariable
     */
    void CalculateStressResponse(ConstitutiveLaw::Parameters& rParametersValues,
                                 Vector& rDamageVector,
                                 Vector& rStrainVariables);

    void TensorProduct6(Matrix& rOutput,
                        const Vector& rVector1,
                        const Vector& rVector2);



    /**
     * @brief Computes the material response in terms of 2nd Piola-Kirchhoff stress
     * @param rValues The specific parameters of the current constitutive law
     * @see Parameters
     */
    void CalculateMaterialResponsePK2(Parameters& rValues) override;

    /**
     * @brief Indicates if this CL requires initialization of the material response,
     * called by the element in InitializeMaterialResponse.
     */
    bool RequiresInitializeMaterialResponse() override
    {
        return false;
    }

    /**
     * @brief Indicates if this CL requires finalization of the material response,
     * called by the element in FinalizeMaterialResponse.
     */
    bool RequiresFinalizeMaterialResponse() override
    {
        return true;
    }

    /**
     * @brief Finalize the material response in terms of Cauchy stresses
     * @param rValues The specific parameters of the current constitutive law
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy(Parameters& rValues) override;

    // /**
    //  * @brief calculates the value of a specified variable (double)
    //  * @param rValues the needed parameters for the CL calculation
    //  * @param rThisVariable the variable to be returned
    //  * @param rValue a reference to the returned value
    //  * @return rValue output: the value of the specified variable
    //  */
    // double& CalculateValue(Parameters& rValues,
    //                        const Variable<double>& rThisVariable,
    //                        double& rValue) override;

    /**
     * @brief calculates the value of a specified variable (Vector)
     * @param rValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
     Vector& CalculateValue(Parameters& rValues,
                            const Variable<Vector>& rThisVariable,
                            Vector& rValue) override;

    /**
     * @brief calculates the damage effect tensor M
     * @return damage effect tensor M
     */
    void GetDamageEffectTensor(BoundedMatrixVoigtType& DamageEffectTensor,
                               const BoundedVectorType& DamageVector
                               );

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input or that no common error is found.
     * @param rMaterialProperties The properties of the material
     * @param rElementGeometry The geometry of the element
     * @param rCurrentProcessInfo The current process info instance
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

    ///@}
    ///@name Inquiry
    ///@{
    ///@}
    ///@name Input and output
    ///@{

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
        rOStream << "Small Strain Isotropic Damage 3D constitutive law\n";
    };

protected:

    ///@name Protected static Member Variables
    ///@{
        const double eps = 1e-8;
    ///@}

    ///@name Protected member Variables
    ///@{
    Vector mDamageVector;
    Vector mStrainVariables = ZeroVector(3);
    ///@}

    ///@name Protected Operators
    ///@{
    ///@}

    ///@name Protected Operations
    ///@{
    /**

     * @brief This method computes principal values of stresses/strains
     * @param VectorForm Stresses/Strains in vector form
     * @param Pri_Values principal values in vector form
     * @param MinValue minimum of the principal values
     */
    void GetEigenValues(BoundedVectorType& Pri_Values,
                        const Variable<Vector>& rThisVariable,
                        const Vector& VectorForm);

    ///@}

    /**
     * @brief This method calculates the linearized tangent operator
     */
    void CalculateParameters(BoundedMatrixVoigtType& EffStiffnessMatrix,
                             BoundedMatrix3x6Type& dEprdE,
                             BoundedMatrixType& dkdEpr,
                             ConstitutiveLaw::Parameters& rParametersValues,
                             const Vector& DamageVector
                             );
    /**
     * @brief This method calculates the linearized tangent operator
     */

    void CalculatePartialDerivatives(array_1d<BoundedMatrix<double, 6, 6>, 3>& dHdk,
                                    const Properties& rMaterialProperties,
                                    const Vector& DamageVector,
                                    const BoundedVectorType& Kappa0,
                                    const BoundedVectorType& Beta1,
                                    const BoundedVectorType& Beta2,
                                    const BoundedVectorType& Kappa
                                    );

    /**
     * @brief This method converts stress or strain vectors to tensors
     */
    void VectorToTensor(BoundedMatrixType& TensorForm,
                        const Vector& VectorForm,
                        const Variable<Vector>& rThisVariable
                        );
    /**
     * @brief the derivatives of eigen values with respect to the matrix elements
     * @param DerivativesofEigenvalues
     * @param EigenvaluesVector
     * @param Voigtform
     */
    void CalculateDerivativesofEigenvalues(BoundedMatrix3x6Type &DerivativesofEigenvalues,
                                           BoundedVectorType &EigenvaluesVector,
                                           const BoundedVectorVoigtType &Voigtform,
                                           const Variable<Vector>& rThisVariable);


    void MultiplyTensors(BoundedMatrixVoigtType& dSdE,
                        const array_1d<BoundedMatrix<double, 6, 6>, 6>& dHdE,
                        const Vector& StrainVector);
    void GetdHdE(array_1d<BoundedMatrix<double, 6, 6>, 6>& dHdE,
                const array_1d<BoundedMatrix<double, 6, 6>, 3>& dHdk,
                const BoundedMatrix3x6Type& dkdE);
    void GetdHdk(array_1d<BoundedMatrix<double, 6, 6>, 3>& dHdk,
                const array_1d<BoundedMatrix<double, 6, 6>, 3>& dHdD,
                const BoundedMatrixType& dDdkappa);
        /**
     * @brief This method computes the tangent tensor
     * @param rValues The constitutive law parameters and flags
     */
    void CalculateTangentTensor(ConstitutiveLaw::Parameters &rValues);

    /**
     * @brief This method computes the secant tensor
     * @param rValues The constitutive law parameters and flags
     */
    void CalculateSecantTensor(ConstitutiveLaw::Parameters& rValues, Matrix& rSecantTensor);
private:

    ///@name Static Member Variables
    ///@{

    ///@}

    ///@name Member Variables
    ///@{

    ///@}

    ///@name Private Operators
    ///@}

    ///@name Private  Access
    ///@{
    ///@}

    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;

    void load(Serializer& rSerializer) override;

}; // class ElasticAnisotropicDamage
} // namespace Kratos
