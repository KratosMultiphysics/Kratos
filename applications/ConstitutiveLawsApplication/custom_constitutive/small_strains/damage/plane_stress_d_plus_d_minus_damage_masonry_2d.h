// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Philip Kalkbrenner
//                   Massimo Petracca
//                   Alejandro Cornejo
//
//
#pragma once

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos
{
/**
* @class DamageDPlusDMinusMasonry2DLaw
* @ingroup StructuralMechanicsApplication
*/
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) DamageDPlusDMinusMasonry2DLaw
	: public ConstitutiveLaw
{
public:

	KRATOS_CLASS_POINTER_DEFINITION(DamageDPlusDMinusMasonry2DLaw);

	///@name Type Definitions
	///@{

	/// We define the working dimension size, already defined in the integrator
	static constexpr SizeType Dimension = 2;

	/// We define the Voigt size, already defined in the  integrator
	static constexpr SizeType VoigtSize = 3;

	/// Definition of the machine precision tolerance
	static constexpr double tolerance = std::numeric_limits<double>::epsilon();

	///@}

	/**
	* Default Constructor.
	*/
	DamageDPlusDMinusMasonry2DLaw();

	/**
	* Destructor.
	*/
	~DamageDPlusDMinusMasonry2DLaw() override
	{
	}

	/**
	* Clone.
	*/
		ConstitutiveLaw::Pointer Clone() const override;

	/**
	* @brief returns the working space dimension of the current constitutive law
	*/
	SizeType WorkingSpaceDimension() override
	{
		return Dimension;
	};

	/**
	* @brief returns the size of the strain vector of the current constitutive law
	*/
	SizeType GetStrainSize() const override
	{
		return VoigtSize;
	};

	struct CalculationData{

		// Elastic Properties
		double YoungModulus;								double PoissonRatio;
		Matrix ElasticityMatrix;

		// Tension Damage Properties
		double YieldStressTension;							double FractureEnergyTension;

		// Compression Damage Properties
		double DamageOnsetStressCompression;				double YieldStressCompression;
		double ResidualStressCompression;					double YieldStrainCompression;
		double BezierControllerC1;							double BezierControllerC2;
		double BezierControllerC3;							double FractureEnergyCompression;
		double BiaxialCompressionMultiplier;				double ShearCompressionReductor;

		// Effective Stress Data
		array_1d<double,3> EffectiveStressVector;			array_1d<double,2> PrincipalStressVector;
		array_1d<double,3> EffectiveTensionStressVector;	array_1d<double,3> EffectiveCompressionStressVector;
		Matrix ProjectionTensorTension;						Matrix ProjectionTensorCompression;

		// Misc
		double CharacteristicLength;						double DeltaTime;
		int TensionYieldModel;

	};

	/**
	* returns whether this constitutive Law has specified variable
	* @param rThisVariable the variable to be checked for
	* @return true if the variable is defined in the constitutive law
	*/
	bool Has(const Variable<double>& rThisVariable) override;

	/**
	* returns whether this constitutive Law has specified variable
	* @param rThisVariable the variable to be checked for
	* @return true if the variable is defined in the constitutive law
	*/
	bool Has(const Variable<Vector>& rThisVariable) override;

	/**
	* returns whether this constitutive Law has specified variable
	* @param rThisVariable the variable to be checked for
	* @return true if the variable is defined in the constitutive law
	*/
	bool Has(const Variable<Matrix>& rThisVariable) override;

	/**
	* returns whether this constitutive Law has specified variable
	* @param rThisVariable the variable to be checked for
	* @return true if the variable is defined in the constitutive law
	* NOTE: fixed size array of 3 doubles (e.g. for 2D stresses, plastic strains, ...)
	*/
	bool Has(const Variable<array_1d<double, 3 > >& rThisVariable) override;

	/**
	* returns whether this constitutive Law has specified variable
	* @param rThisVariable the variable to be checked for
	* @return true if the variable is defined in the constitutive law
	* NOTE: fixed size array of 6 doubles (e.g. for stresses, plastic strains, ...)
	*/
	bool Has(const Variable<array_1d<double, 6 > >& rThisVariable) override;

	/**
	* returns the value of a specified variable
	* @param rThisVariable the variable to be returned
	* @param rValue a reference to the returned value
	* @param rValue output: the value of the specified variable
	*/
	double& GetValue(
		const Variable<double>& rThisVariable,
		double& rValue) override;

	/**
	* returns the value of a specified variable
	* @param rThisVariable the variable to be returned
	* @param rValue a reference to the returned value
	* @return the value of the specified variable
	*/
	Vector& GetValue(
		const Variable<Vector>& rThisVariable,
		Vector& rValue) override;

	/**
	* returns the value of a specified variable
	* @param rThisVariable the variable to be returned
	* @return the value of the specified variable
	*/
	Matrix& GetValue(
		const Variable<Matrix>& rThisVariable,
		Matrix& rValue) override;

	/**
	* returns the value of a specified variable
	* @param rThisVariable the variable to be returned
	* @param rValue a reference to the returned value
	* @return the value of the specified variable
	*/
	array_1d<double, 3 > & GetValue(
		const Variable<array_1d<double, 3 > >& rVariable,
		array_1d<double, 3 > & rValue) override;

	/**
	* returns the value of a specified variable
	* @param rThisVariable the variable to be returned
	* @param rValue a reference to the returned value
	* @return the value of the specified variable
	*/
	array_1d<double, 6 > & GetValue(
		const Variable<array_1d<double, 6 > >& rVariable,
		array_1d<double, 6 > & rValue) override;

	/**
	* sets the value of a specified variable
	* @param rVariable the variable to be returned
	* @param rValue new value of the specified variable
	* @param rCurrentProcessInfo the process info
	*/
	void SetValue(
		const Variable<double>& rVariable,
		const double& rValue,
		const ProcessInfo& rCurrentProcessInfo) override;

	/**
	* sets the value of a specified variable
	* @param rVariable the variable to be returned
	* @param rValue new value of the specified variable
	* @param rCurrentProcessInfo the process info
	*/
	void SetValue(
		const Variable<Vector >& rVariable,
		const Vector& rValue,
		const ProcessInfo& rCurrentProcessInfo) override;

	/**
	* sets the value of a specified variable
	* @param rVariable the variable to be returned
	* @param rValue new value of the specified variable
	* @param rCurrentProcessInfo the process info
	*/
	void SetValue(
		const Variable<Matrix >& rVariable,
		const Matrix& rValue,
		const ProcessInfo& rCurrentProcessInfo) override;

	/**
	* sets the value of a specified variable
	* @param rVariable the variable to be returned
	* @param rValue new value of the specified variable
	* @param rCurrentProcessInfo the process info
	*/
	void SetValue(
		const Variable<array_1d<double, 3 > >& rVariable,
		const array_1d<double, 3 > & rValue,
		const ProcessInfo& rCurrentProcessInfo) override;

	/**
	* sets the value of a specified variable
	* @param rVariable the variable to be returned
	* @param rValue new value of the specified variable
	* @param rCurrentProcessInfo the process info
	*/
	void SetValue(
		const Variable<array_1d<double, 6 > >& rVariable,
		const array_1d<double, 6 > & rValue,
		const ProcessInfo& rCurrentProcessInfo) override;

	/**
	* Is called to check whether the provided material parameters in the Properties
	* match the requirements of current constitutive model.
	* @param rMaterialProperties the current Properties to be validated against.
	* @return true, if parameters are correct; false, if parameters are insufficient / faulty
	* NOTE: this has to be implemented by each constitutive model. Returns false in base class since
	* no valid implementation is contained here.
	*/
	bool ValidateInput(const Properties& rMaterialProperties) override;

	/**
	* returns the expected strain measure of this constitutive law (by default linear strains)
	* @return the expected strain measure
	*/
	StrainMeasure GetStrainMeasure() override;

	/**
	* returns the stress measure of this constitutive law (by default 1st Piola-Kirchhoff stress in voigt notation)
	* @return the expected stress measure
	*/
	StressMeasure GetStressMeasure() override;

	/**
	* returns whether this constitutive model is formulated in incremental strains/stresses
	* NOTE: by default, all constitutive models should be formulated in total strains
	* @return true, if formulated in incremental strains/stresses, false otherwise
	*/
	bool IsIncremental() override;

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
		const Vector& rShapeFunctionsValues) override;

	/**
	* Initialize the material response in terms of 2nd Piola-Kirchhoff stresses
	* @see Parameters
	*/
		void InitializeMaterialResponsePK2 (
			ConstitutiveLaw::Parameters& rValues) override;

	/**
	* to be called at the beginning of each solution step
	* (e.g. from Element::InitializeSolutionStep)
	* @param rMaterialProperties the Properties instance of the current element
	* @param rElementGeometry the geometry of the current element
	* @param rShapeFunctionsValues the shape functions values in the current integration point
	* @param the current ProcessInfo instance
	*/
	void InitializeSolutionStep(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo) override;

	/**
	* to be called at the end of each solution step
	* (e.g. from Element::FinalizeSolutionStep)
	* @param rMaterialProperties the Properties instance of the current element
	* @param rElementGeometry the geometry of the current element
	* @param rShapeFunctionsValues the shape functions values in the current integration point
	* @param the current ProcessInfo instance
	*/
	void FinalizeSolutionStep(
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		const ProcessInfo& rCurrentProcessInfo) override;

	/**
	* Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
	* @see Parameters
	*/
	void CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues) override;

	/**
	* Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
	* @see Parameters
	*/
	void CalculateMaterialResponsePK2 (ConstitutiveLaw::Parameters& rValues) override;

	/**
	* Computes the material response in terms of Kirchhoff stresses and constitutive tensor
	* @see Parameters
	*/
	void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

	/**
	* Computes the material response in terms of Cauchy stresses and constitutive tensor
	* @see Parameters
	*/
	void CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues) override;

	/**
	* Updates the material response in terms of 1st Piola-Kirchhoff stresses
	* @see Parameters
	*/
	void FinalizeMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues) override;

	/**
	* Updates the material response in terms of 2nd Piola-Kirchhoff stresses
	* @see Parameters
	*/
	void FinalizeMaterialResponsePK2 (ConstitutiveLaw::Parameters& rValues) override;

	/**
	* Updates the material response in terms of Kirchhoff stresses
	* @see Parameters
	*/
	void FinalizeMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues) override;

	/**
	* Updates the material response in terms of Cauchy stresses
	* @see Parameters
	*/
	void FinalizeMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues) override;

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
		const Vector& rShapeFunctionsValues) override;

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
		const ProcessInfo& rCurrentProcessInfo) const override;


	void CalculateMaterialResponse(const Vector& StrainVector,
		const Matrix& DeformationGradient,
		Vector& StressVector,
		Matrix& AlgorithmicTangent,
		const ProcessInfo& rCurrentProcessInfo,
		const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const Vector& rShapeFunctionsValues,
		bool CalculateStresses = true,
		int CalculateTangent = true,
		bool SaveInternalVariables = true);

protected:

	///@name Protected member Variables
	///@{

	// Initialization
	bool InitializeDamageLaw = false;

	// Tension & Compression Thresholds

	// for IMPLEX_Integration:
	double PreviousThresholdTension = 0.0;			double PreviousThresholdCompression = 0.0;		    // at time step n - 1:
	// end for IMPLEX_Integration

	double CurrentThresholdTension = 0.0;     		double CurrentThresholdCompression = 0.0;		    // at time step n:
	double ThresholdTension        = 0.0;			double ThresholdCompression        = 0.0;			// at time step n + 1:

	// Damage Parameters & Uniaxial Stresses
	double DamageParameterTension = 0.0;			double DamageParameterCompression = 0.0;
	double UniaxialStressTension  = 0.0;			double UniaxialStressCompression  = 0.0;

	// Misc
	double InitialCharacteristicLength = 0.0;

	// for IMPLEX Integration
	double CurrentDeltaTime  = 0.0;					// at time step n + 1
	double PreviousDeltaTime = 0.0;					// at time step n

	// temporary internal variables saved for the implicit step in the finalize solution step of the IMPLEX
	double TemporaryImplicitThresholdTension = 0.0;
	double TemporaryImplicitThresholdTCompression = 0.0;
    // end for IMPLEX Integration

	///@}


	///@name Protected Operators
	///@{

	/**
	* @brief Initializes the CalculationData at the beginning of each SolutionStep
	* @param Properties& props 		The Material Properties of the constitutive law from rValues
	* 		 GeometryType& geom 	The Element Geometry from rValues
	*		 ProcessInfo& pinfo		The ProcessInfo from rValues
	*		 CalculationData
	*/
	void InitializeCalculationData(
		const Properties& props,
		const GeometryType& geom,
		const ProcessInfo& pinfo,
		CalculationData& data);

	/**
	* @brief Constructs the Linear Elasticity Matrix and stores it in the CalculationData
	* @param CalculationData
	*/
	void CalculateElasticityMatrix(
		CalculationData& data);

	/**
	* @brief Splits the Effective Stress Vector into a positive (tension) and negative (compression) part
	* @param CalculationData
	*/
	void TensionCompressionSplit(
		CalculationData& data);

	/**
	* @brief This method computes the positive (tension) and negative (compression) parts of the ProjectionMatrix
	* @param CalculationData
	*/
	void ConstructProjectionTensors(
		CalculationData& data);

	/**
	* @brief This method computes the equivalent stress in Tension
	* @param CalculationData
	* 		 UniaxialStressTension The variable to be filled with the resulting value
	*/
	void CalculateEquivalentStressTension(
		CalculationData& data,
		double& UniaxialStressTension);

	/**
	* @brief This method computes the equivalent stress in Compression
	* @param CalculationData
	* 		 UniaxialStressCompression The variable to be filled with the resulting value
	*/
	void CalculateEquivalentStressCompression(
		CalculationData& data,
		double& UniaxialStressCompression);

	/**
	* @brief This method computes the damage variable d+ of the tension law
	*        by considering an exponential softening behavior
	* @param CalculationData
	* 		 internal_variable 	The internal variable that considers the considers the state of damage
	*        rDamage 			The final damage variable to be filled by this method
	*/
	void CalculateDamageTension(
		CalculationData& data,
		double internal_variable,
		double& rDamageTension);

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
	* s_j   |              ###     '    ###
	*       |            ###'      '    ' ##
	* s_k   |-----------##--+------+----+--## (K)
	* s_0   |---------##(0) '      '    '   ###
	*       |        ## '   '      '    '    '###
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
	* @brief This method computes the damage variable d- of the compression law
	*        by considering the above explained bezier curved uniaxial behavior
	* @param CalculationData 	Calculation Data for the CL
	* 		 internal_variable	The internal variable that considers the considers the state of damage
	*        rDamage 			The final damage variable to be filled by this method
	*/
	void CalculateDamageCompression(
		CalculationData& data,
		double internal_variable,
		double& rDamage);

	/**
	* @brief This method computes the energy of the uniaxial damage law before regularization
	* @param rBezierEnergy
	*		 rBezierEnergy1
	* 		 double s_p, double s_k, double s_r, double e_p, double e_j, double e_k, double e_r, double e_u
	*						As inputs for the energy calculation
	*/
	void ComputeBezierEnergy(
		double& rBezierEnergy,
		double& rBezierEnergy1,
		double s_p, double s_k, double s_r,
		double e_p, double e_j, double e_k, double e_r, double e_u);

	/**
	* @brief This method evaluates the area below the bezier curves
	* @param x1, x2, x3 x-coordinates of the Bezier control points
	* 		 y1, y2, y3 y-coordinates of the Bezier control points
	*/
	double EvaluateBezierArea(
		double x1, double x2, double x3,
		double y1, double y2, double y3);

	/**
	* @brief This method applies the stretcher to the strains, to regularize the fracture energy
	* @param stretcher			The stretch factor
	* 		 e_p				The reference strain from which on the stretching is applied
	*		 e_j, e_k, e_r, e_u	The strains that have to be modified by the stretcher
	*/
	void ApplyBezierStretcherToStrains(
		double stretcher, double e_p, double& e_j,
		double& e_k, double& e_r, double& e_u);

	/**
	* @brief This method evaluates the damage parameter by considering the Bezier law explained above
	* @param rDamageParameter	The parameter to obtain the damage d-
	* 		 xi					The strain like counterpart
	*		 x1, x2, x3 		x-coordinates of the Bezier control points
	* 		 y1, y2, y3 		y-coordinates of the Bezier control points
	*/
	void EvaluateBezierCurve(
		double& rDamageParameter, double xi,
		double x1, double x2, double x3,
		double y1, double y2, double y3);

	/**
	* @brief This method computes the Charcteristic element length
	* @param GeometryType& geom		The element geometry data from rValues
	*	  	 rCharacteristicLength	The characteristic Length
	*/
	void ComputeCharacteristicLength(
		const GeometryType& geom,
		double& rCharacteristicLength);

	/**
	* @brief This method computes the internal material response strain to stress by applying cl
	* @param strain_vector		The strain vector
	*	  	 stress_vector		The stress vector
	*		 CalculationData	Calculation Data for the CL
	*/
	void CalculateMaterialResponseInternal(
		const Vector& strain_vector,
		Vector& stress_vector,
		CalculationData& data,
		const Properties props);

	/**
     * @brief This method checks whether we are in loading/unloading/damage state
     * @param is_damaging_tension	The damage tension bool
	 *		  is_damaging_tension	The damage compression bool
     */
	void CheckDamageLoadingUnloading(
		bool& is_damaging_tension,
		bool& is_damaging_compression);

	 /**
     * @brief This method computes the tangent tensor
     * @param rValues 			The constitutive law parameters and flags
	 *		  strain_vector		The strain vector
	 *	  	  stress_vector		The stress vector
	 *		  CalculationData	Calculation Data for the CL
     */
	void CalculateTangentTensor(
		Parameters& rValues,
		Vector strain_vector,
		Vector stress_vector,
		CalculationData& data,
		const Properties& props);

	/**
     * @brief This method computes the secant tensor
     * @param rValues 			The constitutive law parameters and flags
	 *		  CalculationData	Calculation Data for the CL
     */
	void CalculateSecantTensor(
		Parameters& rValues,
		CalculationData& data);


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
	///@}

	///@name Private Operators
	///@{
	///@}

	///@name Private Operations
	///@{
	///@}

	///@name Private  Access
	///@{
	///@}

	///@name Serialization
	///@{

	friend class Serializer;

	void save(Serializer& rSerializer) const override
	{
		KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw );
	}

	void load(Serializer& rSerializer) override
	{
		KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw );
	}

	///@}

}; // Class DamageDPlusDMinusMasonry2DLaw

} // namespace Kratos