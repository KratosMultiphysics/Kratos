// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Philip Kalkbrenner
//                   Alejandro Cornejo
//
//
#if !defined(KRATOS_PLANE_STRESS_D_PLUS_D_MINUS_DAMAGE_MASONRY_2D_H_INCLUDED)
#define KRATOS_PLANE_STRESS_D_PLUS_D_MINUS_DAMAGE_MASONRY_2D_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/linear_plane_stress.h"

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
 * @class DamageDPlusDMinusMasonry2DLaw
 * @ingroup StructuralMechanicsApplication
 */
 class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) DamageDPlusDMinusMasonry2DLaw
    : public LinearPlaneStress
{
public:
    ///@name Type Definitions
    ///@{


    /// The define the working dimension size, already defined in the integrator
    static constexpr SizeType Dimension = 2;

    /// The define the Voigt size, already defined in the  integrator
    static constexpr SizeType VoigtSize = 3;

    /// Definition of the base class
    typedef LinearPlaneStress BaseType;

    // Adding the respective using to avoid overload conflicts
    using BaseType::Has;
    using BaseType::GetValue;

    /// Counted pointer of GenericYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(DamageDPlusDMinusMasonry2DLaw);

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
    DamageDPlusDMinusMasonry2DLaw();


    /**
    * Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override
	{
        return Kratos::make_shared<DamageDPlusDMinusMasonry2DLaw>(*this);
    }
    /**
     * @brief Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override
    {
        return Dimension;
    };

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return VoigtSize;
    };

    /**
    * Copy constructor.
    */
    DamageDPlusDMinusMasonry2DLaw(const DamageDPlusDMinusMasonry2DLaw &rOther)
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
    ~DamageDPlusDMinusMasonry2DLaw() override
    {
    }

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
    bool IntegrateStressTensionIfNecessary(
        const double F_tension,
        DamageParameters& Parameters,
        array_1d<double, VoigtSize>& IntegratedStressVectorTension,
        const array_1d<double,3> rIntegratedStressVector,
        ConstitutiveLaw::Parameters& rValues
        );

    /**
     * @brief Integrates the predictive tension stress vector if necessary
     * @param F_compression = uniaxial_stress_compression - threshold
     */
    bool IntegrateStressCompressionIfNecessary(
        const double F_compression,
        DamageParameters& Parameters,
        array_1d<double, VoigtSize>& IntegratedStressVectorCompression,
        array_1d<double, VoigtSize> rIntegratedStressVector,
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
     * @brief Sets the value of a specified variable (Vector)
     * @param rVariable the variable to be returned
     * @param rValue new value of the specified variable
     * @param rCurrentProcessInfo the process info
     */
    void SetValue(
        const Variable<Vector> &rThisVariable,
        const Vector& rValue,
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

    void SetTensionStress(const double toS){mTensionUniaxialStress = toS;}
    void SetCompressionStress(const double toS){mCompressionUniaxialStress = toS;}

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

    // Non Converged values
    double mNonConvTensionDamage = 0.0;
    double mNonConvTensionThreshold = 0.0;

    double mCompressionDamage = 0.0;
    double mCompressionThreshold = 0.0;

    // Non Converged values
    double mNonConvCompressionDamage = 0.0;
    double mNonConvCompressionThreshold = 0.0;

    double mTensionUniaxialStress = 0.0;
    double mCompressionUniaxialStress = 0.0;
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

    /**
     * @brief This method computes the secant tensor
     * @param rValues The constitutive law parameters and flags
     */
    void CalculateSecantTensor(ConstitutiveLaw::Parameters& rValues, Matrix& rSecantTensor);


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

	/**
	 * @brief This method computes the equivalent stress in Tension
	 * @param rValues The constitutive law parameters and flags
	 * 		  rPredictiveStressVector Predictive or effective Stress Vector
	 *        rEquivalentStress The equivalent Stress to be filled by method
	 */
	void CalculateEquivalentStressTension(
	array_1d<double, 3>& rPredictiveStressVector,
	double& rEquivalentStress,
	ConstitutiveLaw::Parameters& rValues);

	/**
	 * @brief This method computes the equivalent stress in Compression
	 * @param rValues The constitutive law parameters and flags
	 * 		  rPredictiveStressVector Predictive or effective Stress Vector
	 *        rEquivalentStress The equivalent Stress to be filled by method
	 */
	void CalculateEquivalentStressCompression(
	array_1d<double, 3>& rPredictiveStressVector,
	double& rEquivalentStress,
	ConstitutiveLaw::Parameters& rValues);

	/**
	 * @brief This method computes the final stress vector in Tension
	 * @param rValues The constitutive law parameters and flags
	 * 		  rPredictiveStressVector Tension Part of the predictive or effective stress vector
	 *        UniaxialStress The equivalent uniaxial stress in Tension
	 *        rDamage The damage variable in Tension
	 *		  rThreshold The Damage Threshold in Tension
	 *		  CharacteristicLength The finite element charecteristic length
	 */
	void IntegrateStressVectorTension(
	array_1d<double,3>& rPredictiveStressVector,
    const double UniaxialStress,
    double& rDamage,
    double& rThreshold,
    ConstitutiveLaw::Parameters& rValues,
    const double CharacteristicLength);

	/**
	 *  @brief This method computes the damage parameter for the exponential softening behavior in tension
	 *  @params rValues The constitutive law parameters and flags
	 *          rAParameter The damage parameter filled by method
	 *          CharacteristicLength The finite element charecteristic length
	 */
	void CalculateDamageParameterTension(
	ConstitutiveLaw::Parameters& rValues,
    double& rAParameter,
    const double CharacteristicLength);

	/**
	 * @brief This method computes the tension damage variable for the exponential softening law in tension
	 * @params rValues The constitutive law parameters and flags
	 * 		   UniaxialStress The equivalent uniaxial stress in Tension
	 *         Threshold The damage threshold in Tension
	 *         DamageParameter The damage parameter for the exponential softening law
	 *         CharacteristicLength The finite element charecteristic length
	 *         rDamage The tension damage variable filled by the method
	 */

	void CalculateExponentialDamageTension(
	const double UniaxialStress,
	const double Threshold,
	const double DamageParameter,
	const double CharacteristicLength,
	ConstitutiveLaw::Parameters& rValues,
	double& rDamage);

	/**
	 * @brief This method computes the final stress vector in Tension
	 * @param rValues The constitutive law parameters and flags
	 * 		  rPredictiveStressVector Compression Part of the predictive or effective stress vector
	 *        UniaxialStress The equivalent uniaxial stress in Compression
	 *        rDamage The damage variable in Compression
	 *		  rThreshold The Damage Threshold in Compression
	 *		  CharacteristicLength The finite element charecteristic length
	 */
	void IntegrateStressVectorCompression(
	array_1d<double,3>& rPredictiveStressVector,
	const double UniaxialStress,
    double& rDamage,
    double& rThreshold,
    ConstitutiveLaw::Parameters& rValues,
    const double CharacteristicLength);


    /**
     *  BRIEF DOCUMENTATION OF THE USED UNIAXIAL SOFTENING BEHAVIOR IN COMPRESSION
     *  Entire documentation can be found in the the Phd Thesis of Massimo Petracca
     *  << Computational Multiscale Analysis of Masonry Structures>>
     *
     *  UNIAXIAL BEZIER COMPRESSION DAMAGE
     *  {I}   Linear Elastic
     *  {II}  Hardening Quadratic Bezier Curve
     *          Control nodes:  0=(e_0,s_0); I=(e_i,s_p); P=(e_p,s_p)
     *  {III} Softening Quadratic Bezier Curve
     *          Control nodes:  P=(e_p,s_p); J=(e_j,s_j); K=(e_k,s_k)
     *  {IV}  Softening Quadratic Bezier Curve
     *          Control nodes:  K=(e_k,s_k); R=(e_r,s_r); U=(e_u,s_u)
     *  {V}   Residual Strength
     *
     *    STRESS
     *       ^
     *      /|\
     *       |                     (P)
     * s_p = |------------(I)+----#####--+(J)
     * s_i = |               ' ###  ' ####
     * s_j   |              ###     '    ####
     *       |            ###'      '    ' ###
     * s_k   |-----------##--+------+----+--## (K)
     * s_0   |---------##(0) '      '    '   ###
     *       |        ## '   '      '    '    '##
     *       |       ##  '   '      '    '    '   ####
     *       |      ##   '   '      '    '    '      #####
     *       |     ##    '   '      '    '    '          #####
     *       |    ##     '   '      '    '    '    (R)       ######## (U)
     * s_r = |---##------+---+------'----+----+-----+-----------------######################
     * s_u   |  ##       '   '      '    '    '     '                 '
     *       |_##________+___+______+____+____+_____+_________________+______________________________\
     *                  e_0 e_i    e_p  e_j  e_k   e_r               e_u                             / STRAIN
     *        '          '          '         '                       '
     *        '   {I}    '   {II}   '  {III}  '        {IV}           '          {V}
     *        '          '          '         '                       '
     *
     */

    /**
	 * @brief This method computes the Damage Variable in Compression by considering three Bezier curves (hardening + softening + softening + residual)
	 * @param rValues The constitutive law parameters and flags
	 *        UniaxialStress The equivalent uniaxial stress in Compression
	 *        rDamage The damage variable in Compression
	 *		  rThreshold The Damage Threshold in Compression
	 *        CharacteristicLength The finite element charecteristic length
	 */

	void CalculateBezier3DamageCompression(
	const double UniaxialStress,
	double& rDamage,
	double& rThreshold,
	const double CharacteristicLength,
	ConstitutiveLaw::Parameters& rValues);

	/**
	 * @brief This method regulates the four bezier control strains to avoid a constitutive snap-back (fracture energy considerations)
	 * @param specific_dissipated_fracture_energy FRACTURE_ENERGY_CMOPRESSION devided by CharacteristicLength
	 *        sp, sk, sr stress Values to control the bezier curves
	 *        ep strain pproperty to control the bezier curve
     *        ej, ek, er, eu strain properties to be regulated in method
	 */
	void RegulateBezierDeterminators(
	const double specific_dissipated_fracture_energy,
	const double sp, const double sk, const double sr, const double ep,
	double& ej, double& ek, double& er, double& eu);

	/**
	 * @brief This method computes the area beneath the parts of the bezier curves, respectively
     * @param BezierG Area beneath the curve, to be filled by method
	 *        x1, x2, x3, y1, y2, y3 coordinates of the control points of the bezier
	 */
	void ComputeBezierEnergy(
	double& BezierG,
	const double x1, const double x2, const double x3,
	const double y1, const double y2, const double y3);

    /**
     * @brief This method returns the bezier damage parameter
     * @param Xi Strain-like counterpart of the uniaxial compression stress
     *        x1, x2, x3 Necesarry Stress values to define the uniaxial compression damage bezier curve
     *        y1, y2, y3 Necesarry Strain vlaues to define the uniaxial compression damage bezier curve
     */
    double EvaluateBezierCurve(
    const double Xi,
    const double x1, double x2, const double x3,
    const double y1, const double y2, const double y3);


    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{



    ///@}

}; // Class DamageDPlusDMinusMasonry2DLaw
}// namespace Kratos
#endif

