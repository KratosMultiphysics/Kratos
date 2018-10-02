// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo 
//  
//

#if !defined(KRATOS_GENERIC_SMALL_STRAIN_D_PLUS_D_MINUS_DAMAGE_H_INCLUDED)
#define KRATOS_GENERIC_SMALL_STRAIN_D_PLUS_D_MINUS_DAMAGE_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_constitutive/linear_plane_strain.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    // The size type definition
    typedef std::size_t SizeType;
    
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
 * @class GenericSmallStrainDplusDminusDamage
 * @ingroup StructuralMechanicsApplication
 * @brief This class is the base class which define all the constitutive laws for damage in small deformation
 * @details This class considers a constitutive law integrator as an intermediate utility to compute the damage
 * @tparam TConstLawIntegratorType The constitutive law integrator considered
 * @author Alejandro Cornejo 
 */
template <class TConstLawIntegratorTensionType, class TConstLawIntegratorCompressionType>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) GenericSmallStrainDplusDminusDamage
    : public std::conditional<TConstLawIntegratorTensionType::VoigtSize == 6, ElasticIsotropic3D, LinearPlaneStrain >::type
{
public:
    ///@name Type Definitions
    ///@{

    /// The define the working dimension size, already defined in the integrator
    static constexpr SizeType Dimension = TConstLawIntegratorTensionType::Dimension;

    /// The define the Voigt size, already defined in the  integrator
    static constexpr SizeType VoigtSize = TConstLawIntegratorTensionType::VoigtSize;
    
    /// Definition of the base class
    typedef typename std::conditional<VoigtSize == 6, ElasticIsotropic3D, LinearPlaneStrain >::type BaseType;

    /// Counted pointer of GenericYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(GenericSmallStrainDplusDminusDamage);

    /// The node definition
    typedef Node<3> NodeType;
    
    /// The geometry definition
    typedef Geometry<NodeType> GeometryType;
    
    /// Definition of the machine precision tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

	struct DamageParameters {
		double DamageTension = 0.0;
		double DamageCompression = 0.0;
		double ThresholdTension = 0.0;
		double ThresholdCompression = 0.0;
		array_1d<double, VoigtSize> TensionStressVector;
		array_1d<double, VoigtSize> CompressionStressVector;
        double UniaxialTensionStress = 0.0;
        double UniaxialCompressionStress = 0.0;
	};
    ///@}
    ///@name Life Cycle
    ///@{

    /**
    * Default constructor.
    */
    GenericSmallStrainDplusDminusDamage()
    {
    }

    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
    {
        return Kratos::make_shared<GenericSmallStrainDplusDminusDamage<TConstLawIntegratorTensionType, TConstLawIntegratorCompressionType>>(*this);
    }

    /**
    * Copy constructor.
    */
    GenericSmallStrainDplusDminusDamage(const GenericSmallStrainDplusDminusDamage &rOther)
        : BaseType(rOther),
          mTensionDamage(rOther.mTensionDamage),
          mTensionThreshold(rOther.mTensionThreshold),
          mNonConvTensionDamage(rOther.mNonConvTensionDamage),
          mNonConvTensionThreshold(rOther.mNonConvTensionThreshold),
          mCompressionDamage(rOther.mCompressionDamage),
          mCompressionThreshold(rOther.mCompressionThreshold),
          mNonConvCompressionDamage(rOther.mNonConvCompressionDamage),
          mNonConvCompressionThreshold(rOther.mNonConvCompressionThreshold)
    {
    }

    /**
    * Destructor.
    */
    ~GenericSmallStrainDplusDminusDamage() override
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Computes the material response in terms of Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Integrates the predictive tension stress vector if necessary
     * @param F_compression = uniaxial_stress_tension - threshold
     */
    void IntegrateStressTensionIfNecessary(
        const double F_tension, 
        DamageParameters& Parameters, 
        array_1d<double, VoigtSize>& IntegratedStressVectorTension,
        ConstitutiveLaw::Parameters& rValues
        );

    /**
     * @brief Integrates the predictive tension stress vector if necessary
     * @param F_compression = uniaxial_stress_compression - threshold
     */
    void IntegrateStressCompressionIfNecessary(
        const double F_compression, 
        DamageParameters& Parameters,
        array_1d<double, VoigtSize>& IntegratedStressVectorCompression,
        ConstitutiveLaw::Parameters& rValues
        );

    /**
     * @brief Computes the inetgarted stress vector S = A:D0:A:E
     */
    void CalculateIntegratedStressVector(
        Vector& IntegratedStressVectorTension,
        const DamageParameters& Parameters,
        ConstitutiveLaw::Parameters& rValues
        );

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
     * @brief To be called at the end of each solution step
     * @details (e.g. from Element::FinalizeSolutionStep)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param rCurrentProcessInfo the current ProcessInfo instance
     */
    void FinalizeSolutionStep(
        const Properties &rMaterialProperties,
        const GeometryType &rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Finalize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Finalize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Finalize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues) override;
    /**
     * Finalize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (double)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<double> &rThisVariable) override;

    /**
     * @brief Returns whether this constitutive Law has specified variable (Vector)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Vector> &rThisVariable) override;
    
    /**
     * @brief Returns whether this constitutive Law has specified variable (Matrix)
     * @param rThisVariable the variable to be checked for
     * @return true if the variable is defined in the constitutive law
     */
    bool Has(const Variable<Matrix> &rThisVariable) override;

    /**
     * @brief Sets the value of a specified variable (double)
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<double> &rThisVariable,
        const double& rValue,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief Returns the value of a specified variable (double)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& GetValue(
        const Variable<double> &rThisVariable,
        double& rValue
        ) override;

    /**
     * @brief Returns the value of a specified variable (Vector)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Vector& GetValue(
        const Variable<Vector> &rThisVariable,
        Vector& rValue
        ) override;

    /**
     * @brief Returns the value of a specified variable (matrix)
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    Matrix& GetValue(
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue
        ) override;

    /**
     * @brief Returns the value of a specified variable (double)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    double& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<double>& rThisVariable,
        double& rValue) override;

    /**
     * @brief Returns the value of a specified variable (vector)
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
     * @brief Returns the value of a specified variable (matrix)
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
     * of the calculations, so to verify that nothing is missing from the input or that no common error is found.
     * @param rMaterialProperties The properties of the material
     * @param rElementGeometry The geometry of the element
     * @param rCurrentProcessInfo The current process info instance
     * @return 0 if OK, 1 otherwise
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * @brief This method computes the elastic tensor
     * @param rLinearElasticMatrix The elastic tensor
     * @param YoungModulus The properties of the material
     * @param PoissonCoefficient The poisson coefficient
     * @return 0 if OK, 1 otherwise
     */
    void CalculateLinearElasticMatrix(
        Matrix& rLinearElasticMatrix, 
        const double YoungModulus, 
        const double PoissonCoefficient
        );
    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

    ///@}

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
    // Tension values
    double& GetTensionThreshold() { return mTensionThreshold; }
    double& GetTensionDamage() { return mTensionDamage; }
    double& GetNonConvTensionThreshold() { return mNonConvTensionThreshold; }
    double& GetNonConvTensionDamage() { return mNonConvTensionDamage; }

    void SetTensionThreshold(const double toThreshold) { mTensionThreshold = toThreshold; }
    void SetTensionDamage(const double toDamage) { mTensionDamage = toDamage; }
    void SetNonConvTensionThreshold(const double toThreshold) { mNonConvTensionThreshold = toThreshold; }
    void SetNonConvTensionDamage(const double toDamage) { mNonConvTensionDamage = toDamage; }

    // Compression values
    double& GetCompressionThreshold() { return mCompressionThreshold; }
    double& GetCompressionDamage() { return mCompressionDamage; }
    double& GetNonConvCompressionThreshold() { return mNonConvCompressionThreshold; }
    double& GetNonConvCompressionDamage() { return mNonConvCompressionDamage; }

    void SetCompressionThreshold(const double toThreshold) { mCompressionThreshold = toThreshold; }
    void SetCompressionDamage(const double toDamage) { mCompressionDamage = toDamage; }
    void SetNonConvCompressionThreshold(const double toThreshold) { mNonConvCompressionThreshold = toThreshold; }
    void SetNonConvCompressionDamage(const double toDamage) { mNonConvCompressionDamage = toDamage; }
    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    // Converged values
    double mTensionDamage = 0.0;
    double mTensionThreshold = 0.0;
    // double mUniaxialStress = 0.0;

    // Non Converged values
    double mNonConvTensionDamage = 0.0;
    double mNonConvTensionThreshold = 0.0;

    double mCompressionDamage = 0.0;
    double mCompressionThreshold = 0.0;
    // double mUniaxialStress = 0.0;

    // Non Converged values
    double mNonConvCompressionDamage = 0.0;
    double mNonConvCompressionThreshold = 0.0;
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method computes the tangent tensor
     * @param rValues The constitutive law parameters and flags
     */
    void CalculateTangentTensor(ConstitutiveLaw::Parameters &rValues);

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    // Serialization

    friend class Serializer;

    void save(Serializer &rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.save("TensionDamage", mTensionDamage);
        rSerializer.save("TensionThreshold", mTensionThreshold);
        rSerializer.save("NonConvTensionDamage", mNonConvTensionDamage);
        rSerializer.save("NonConvTensionThreshold", mNonConvTensionThreshold);
        rSerializer.save("CompressionDamage", mCompressionDamage);
        rSerializer.save("CompressionThreshold", mCompressionThreshold);
        rSerializer.save("NonConvCompressionnDamage", mNonConvCompressionDamage);
        rSerializer.save("NonConvCompressionThreshold", mNonConvCompressionThreshold);
    }

    void load(Serializer &rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.load("TensionDamage", mTensionDamage);
        rSerializer.load("TensionThreshold", mTensionThreshold);
        rSerializer.load("NonConvTensionDamage", mNonConvTensionDamage);
        rSerializer.load("NonConvTensionThreshold", mNonConvTensionThreshold);
        rSerializer.load("CompressionDamage", mCompressionDamage);
        rSerializer.load("CompressionThreshold", mCompressionThreshold);
        rSerializer.load("NonConvCompressionnDamage", mNonConvCompressionDamage);
        rSerializer.load("NonConvCompressionThreshold", mNonConvCompressionThreshold);
    }

    ///@}

}; // Class 

} // namespace Kratos
#endif
