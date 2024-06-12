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
//					 Massimo Petracca
//                   Alejandro Cornejo
//
//

// System includes

// External includes

// Project includes
#include "plane_stress_d_plus_d_minus_damage_masonry_2d.h"
#include "includes/model_part.h"
#include "constitutive_laws_application_variables.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"


#define OPTIMIZE_CHARACTERISTIC_LENGTH
#define HEAVISIDE(X) ( X >= 0.0 ? 1.0 : 0.0)
#define MACAULAY(X)  ( X >= 0.0 ? X : 0.0)
#define PROJECTION_OPERATOR_CERVERA_2003
//#define PROJECTION_OPERATOR_CERVERA_2017

namespace Kratos
{
/***********************************************************************************/
/***********************************************************************************/
DamageDPlusDMinusMasonry2DLaw::DamageDPlusDMinusMasonry2DLaw()
	: ConstitutiveLaw()
{
}
/***********************************************************************************/
/***********************************************************************************/
ConstitutiveLaw::Pointer DamageDPlusDMinusMasonry2DLaw::Clone() const
{
	return ConstitutiveLaw::Pointer( new DamageDPlusDMinusMasonry2DLaw() );
}
/***********************************************************************************/
/***********************************************************************************/
bool DamageDPlusDMinusMasonry2DLaw::Has(
	const Variable<double>& rThisVariable)
{
	if(rThisVariable == DAMAGE_TENSION)
		return true;
	if(rThisVariable == UNIAXIAL_STRESS_TENSION)
		return true;
	if(rThisVariable == THRESHOLD_TENSION)
		return true;
	if(rThisVariable == DAMAGE_COMPRESSION)
		return true;
	if(rThisVariable == UNIAXIAL_STRESS_COMPRESSION)
		return true;
	if(rThisVariable == THRESHOLD_COMPRESSION)
		return true;
	return false;
}
/***********************************************************************************/
/***********************************************************************************/

bool DamageDPlusDMinusMasonry2DLaw::Has(
	const Variable<Vector>& rThisVariable)
{
	if(rThisVariable == INTERNAL_VARIABLES)
		return true;
	return false;
}
/***********************************************************************************/
/***********************************************************************************/
bool DamageDPlusDMinusMasonry2DLaw::Has(
	const Variable<Matrix>& rThisVariable)
{
	return false;
}
/***********************************************************************************/
/***********************************************************************************/
bool DamageDPlusDMinusMasonry2DLaw::Has(
	const Variable<array_1d<double, 3 > >& rThisVariable)
{
	return false;
}
/***********************************************************************************/
/***********************************************************************************/
bool DamageDPlusDMinusMasonry2DLaw::Has(
	const Variable<array_1d<double, 6 > >& rThisVariable)
{
	return false;
}
/***********************************************************************************/
/***********************************************************************************/
double& DamageDPlusDMinusMasonry2DLaw::GetValue(
	const Variable<double>& rThisVariable,
	double& rValue)
{
	rValue = 0.0;
	if(rThisVariable == DAMAGE_TENSION)
		rValue = DamageParameterTension;
	else if(rThisVariable == DAMAGE_COMPRESSION)
		rValue = DamageParameterCompression;
	else if(rThisVariable == UNIAXIAL_STRESS_TENSION)
		rValue = UniaxialStressTension;
	else if(rThisVariable == UNIAXIAL_STRESS_COMPRESSION)
		rValue = UniaxialStressCompression;
	else if(rThisVariable == THRESHOLD_TENSION)
		rValue = ThresholdTension;
	else if(rThisVariable == THRESHOLD_COMPRESSION)
		rValue = ThresholdCompression;
	return rValue;
}
/***********************************************************************************/
/***********************************************************************************/
Vector& DamageDPlusDMinusMasonry2DLaw::GetValue(
	const Variable<Vector>& rThisVariable,
	Vector& rValue)
{
	return rValue;
}
/***********************************************************************************/
/***********************************************************************************/
Matrix& DamageDPlusDMinusMasonry2DLaw::GetValue(
	const Variable<Matrix>& rThisVariable,
	Matrix& rValue)
{
	return rValue;
}
/***********************************************************************************/
/***********************************************************************************/
array_1d<double, 3 > & DamageDPlusDMinusMasonry2DLaw::GetValue(
	const Variable<array_1d<double, 3 > >& rVariable,
	array_1d<double, 3 > & rValue)
{
	return rValue;
}
/***********************************************************************************/
/***********************************************************************************/
array_1d<double, 6 > & DamageDPlusDMinusMasonry2DLaw::GetValue(
	const Variable<array_1d<double, 6 > >& rVariable,
	array_1d<double, 6 > & rValue)
{
	return rValue;
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::SetValue(
	const Variable<double>& rVariable,
	const double& rValue,
	const ProcessInfo& rCurrentProcessInfo)
{
	if(rVariable == DAMAGE_TENSION)
		DamageParameterTension = rValue;
	else if(rVariable == DAMAGE_COMPRESSION)
		DamageParameterCompression = rValue;
	else if(rVariable == UNIAXIAL_STRESS_TENSION)
		UniaxialStressTension = rValue;
	else if(rVariable == UNIAXIAL_STRESS_COMPRESSION)
		UniaxialStressCompression = rValue;
	else if(rVariable == THRESHOLD_TENSION)
		ThresholdTension = rValue;
	else if(rVariable == THRESHOLD_COMPRESSION)
		ThresholdCompression = rValue;
}

/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::SetValue(
	const Variable<Vector >& rVariable,
	const Vector& rValue,
	const ProcessInfo& rCurrentProcessInfo)
{
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::SetValue(
	const Variable<Matrix >& rVariable,
	const Matrix& rValue,
	const ProcessInfo& rCurrentProcessInfo)
{
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::SetValue(
	const Variable<array_1d<double, 3 > >& rVariable,
	const array_1d<double, 3 > & rValue,
	const ProcessInfo& rCurrentProcessInfo)
{
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::SetValue(
	const Variable<array_1d<double, 6 > >& rVariable,
	const array_1d<double, 6 > & rValue,
	const ProcessInfo& rCurrentProcessInfo)
{
}
/***********************************************************************************/
/***********************************************************************************/
bool DamageDPlusDMinusMasonry2DLaw::ValidateInput(
	const Properties& rMaterialProperties)
{
	if( !rMaterialProperties.Has(YOUNG_MODULUS) ) 					return false;
	if( !rMaterialProperties.Has(POISSON_RATIO) ) 					return false;
	if( !rMaterialProperties.Has(YIELD_STRESS_TENSION) ) 			return false;
	if( !rMaterialProperties.Has(FRACTURE_ENERGY_TENSION) ) 		return false;
	if( !rMaterialProperties.Has(DAMAGE_ONSET_STRESS_COMPRESSION) ) return false;
	if( !rMaterialProperties.Has(YIELD_STRESS_COMPRESSION) ) 		return false;
	if( !rMaterialProperties.Has(RESIDUAL_STRESS_COMPRESSION) )	 	return false;
	if( !rMaterialProperties.Has(YIELD_STRAIN_COMPRESSION) ) 		return false;
	if( !rMaterialProperties.Has(FRACTURE_ENERGY_COMPRESSION) ) 	return false;
	if( !rMaterialProperties.Has(BIAXIAL_COMPRESSION_MULTIPLIER) ) 	return false;
	return true;
}
/***********************************************************************************/
/***********************************************************************************/
DamageDPlusDMinusMasonry2DLaw::StrainMeasure DamageDPlusDMinusMasonry2DLaw::GetStrainMeasure()
{
	return ConstitutiveLaw::StrainMeasure_Infinitesimal;
}
/***********************************************************************************/
/***********************************************************************************/
DamageDPlusDMinusMasonry2DLaw::StressMeasure DamageDPlusDMinusMasonry2DLaw::GetStressMeasure()
{
	return ConstitutiveLaw::StressMeasure_Cauchy;
}
/***********************************************************************************/
/***********************************************************************************/
bool DamageDPlusDMinusMasonry2DLaw::IsIncremental()
{
	return false;
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::InitializeMaterial(
	const Properties& rMaterialProperties,
	const GeometryType& rElementGeometry,
	const Vector& rShapeFunctionsValues)
{
	if (InitializeDamageLaw == false){
		ThresholdTension            	= rMaterialProperties[YIELD_STRESS_TENSION];
		CurrentThresholdTension   		= ThresholdTension;
		ThresholdCompression        	= rMaterialProperties[DAMAGE_ONSET_STRESS_COMPRESSION];
		CurrentThresholdCompression 	= ThresholdCompression;
		DamageParameterTension      	= 0.0;
		DamageParameterCompression  	= 0.0;
		UniaxialStressTension			= 0.0;
		UniaxialStressCompression		= 0.0;

		this->ComputeCharacteristicLength(rElementGeometry, InitialCharacteristicLength);

		// Begin IMPLEX Integration - Only if switched on
		if (rMaterialProperties[INTEGRATION_IMPLEX] != 0){
			PreviousThresholdTension 		= ThresholdTension;
			PreviousThresholdCompression 	= ThresholdCompression;
			CurrentDeltaTime 				= 0.0;
			PreviousDeltaTime 				= 0.0;
	    }
		// End IMPLEX Integration

		InitializeDamageLaw    			= true;
	}
}

/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::InitializeSolutionStep(
	const Properties& rMaterialProperties,
	const GeometryType& rElementGeometry,
	const Vector& rShapeFunctionsValues,
	const ProcessInfo& rCurrentProcessInfo)
{
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::FinalizeSolutionStep(
	const Properties& rMaterialProperties,
	const GeometryType& rElementGeometry,
	const Vector& rShapeFunctionsValues,
	const ProcessInfo& rCurrentProcessInfo)
{

}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CalculateMaterialResponsePK1 (
	Parameters& rValues)
{
	CalculateMaterialResponseCauchy(rValues);
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CalculateMaterialResponsePK2 (
	Parameters& rValues)
{
	CalculateMaterialResponseCauchy(rValues);
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CalculateMaterialResponseKirchhoff (
	Parameters& rValues)
{
	CalculateMaterialResponseCauchy(rValues);
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CalculateMaterialResponseCauchy (
	Parameters& rValues)
{
	const ProcessInfo&  pinfo = rValues.GetProcessInfo();
	const GeometryType& geom  = rValues.GetElementGeometry();
	const Properties&   props = rValues.GetMaterialProperties();
	Matrix r_tangent_tensor = rValues.GetConstitutiveMatrix();

	const Vector& StrainVector   	= rValues.GetStrainVector();
	Vector& PredictiveStressVector	= rValues.GetStressVector();

	Matrix T(3, 3);
	noalias(T) = GetTransformationMatrix(rValues.GetMaterialProperties());
    const Vector strain_iso = prod(T, StrainVector);
	Vector stress_iso(3);

	CalculationData data;
	this->InitializeCalculationData(props, geom, pinfo, data);

	this->CalculateMaterialResponseInternal(strain_iso, stress_iso, data, props);
	PredictiveStressVector = prod(trans(T), stress_iso);

	bool is_damaging_tension = false;
	bool is_damaging_compression = false;
	this->CheckDamageLoadingUnloading(is_damaging_tension, is_damaging_compression);

	// Computation of the Constitutive Tensor
	if (rValues.GetOptions().Is(COMPUTE_CONSTITUTIVE_TENSOR)) {
		if (is_damaging_tension || is_damaging_compression) {
			this->CalculateTangentTensor(rValues, strain_iso, stress_iso, data, props);
		} else {
			this->CalculateSecantTensor(rValues, data);
		}
		r_tangent_tensor = prod(trans(T), Matrix(prod(r_tangent_tensor, T)));
	}
}

/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::FinalizeMaterialResponsePK1 (
	Parameters& rValues)
{
	FinalizeMaterialResponseCauchy(rValues);
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::FinalizeMaterialResponsePK2 (
	Parameters& rValues)
{
	FinalizeMaterialResponseCauchy(rValues);
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::FinalizeMaterialResponseKirchhoff (
	Parameters& rValues)
{
	FinalizeMaterialResponseCauchy(rValues);
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::FinalizeMaterialResponseCauchy (
	Parameters& rValues)
{
	// Begin IMPLEX Integration - Only if switched on
	if (rValues.GetMaterialProperties()[INTEGRATION_IMPLEX] != 0){
		ThresholdTension 		= TemporaryImplicitThresholdTension;
		ThresholdCompression 	= TemporaryImplicitThresholdTCompression;

		// move from n to n-1
		PreviousThresholdTension  		= CurrentThresholdTension;
		PreviousThresholdCompression  	= CurrentThresholdCompression;
		PreviousDeltaTime 				= CurrentDeltaTime;
	}
	// End IMPLEX Integration

	// save converged values
	CurrentThresholdTension 		= ThresholdTension;
	CurrentThresholdCompression 	= ThresholdCompression;
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::ResetMaterial(
	const Properties& rMaterialProperties,
	const GeometryType& rElementGeometry,
	const Vector& rShapeFunctionsValues)
{
	ThresholdTension 			= 0.0;
	CurrentThresholdTension 	= 0.0;
	ThresholdCompression 		= 0.0;
	CurrentThresholdCompression = 0.0;
	DamageParameterTension 		= 0.0;
	DamageParameterCompression 	= 0.0;
	InitialCharacteristicLength = 0.0;
	InitializeDamageLaw 		= false;
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::GetLawFeatures(
	Features& rFeatures)
{
	//Set the type of law
	rFeatures.mOptions.Set( PLANE_STRESS_LAW );
	rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);

	//Set the strain size
	rFeatures.mStrainSize = GetStrainSize();

	//Set the space dimension
	rFeatures.mSpaceDimension = WorkingSpaceDimension();
}
/***********************************************************************************/
/***********************************************************************************/
int DamageDPlusDMinusMasonry2DLaw::Check(
	const Properties& rMaterialProperties,
	const GeometryType& rElementGeometry,
	const ProcessInfo& rCurrentProcessInfo) const
{
	KRATOS_TRY

	if( !rMaterialProperties.Has(YOUNG_MODULUS) )
		KRATOS_THROW_ERROR(std::logic_error, "Missing variable: YOUNG_MODULUS", "");

	if( !rMaterialProperties.Has(POISSON_RATIO) )
		KRATOS_THROW_ERROR(std::logic_error, "Missing variable: POISSON_RATIO", "");

	if( !rMaterialProperties.Has(YIELD_STRESS_TENSION) )
		KRATOS_THROW_ERROR(std::logic_error, "Missing variable: YIELD_STRESS_TENSION", "");

	if( !rMaterialProperties.Has(FRACTURE_ENERGY_TENSION) )
		KRATOS_THROW_ERROR(std::logic_error, "Missing variable: FRACTURE_ENERGY_TENSION", "");

	if( !rMaterialProperties.Has(DAMAGE_ONSET_STRESS_COMPRESSION) )
		KRATOS_THROW_ERROR(std::logic_error, "Missing variable: DAMAGE_ONSET_STRESS_COMPRESSION", "");

	if( !rMaterialProperties.Has(YIELD_STRESS_COMPRESSION) )
		KRATOS_THROW_ERROR(std::logic_error, "Missing variable: YIELD_STRESS_COMPRESSION", "");

	if( !rMaterialProperties.Has(RESIDUAL_STRESS_COMPRESSION) )
		KRATOS_THROW_ERROR(std::logic_error, "Missing variable: RESIDUAL_STRESS_COMPRESSION", "");

	if( !rMaterialProperties.Has(YIELD_STRAIN_COMPRESSION) )
		KRATOS_THROW_ERROR(std::logic_error, "Missing variable: YIELD_STRAIN_COMPRESSION", "");

	if( !rMaterialProperties.Has(FRACTURE_ENERGY_COMPRESSION) )
		KRATOS_THROW_ERROR(std::logic_error, "Missing variable: FRACTURE_ENERGY_COMPRESSION", "");

	if( !rMaterialProperties.Has(BIAXIAL_COMPRESSION_MULTIPLIER) )
		KRATOS_THROW_ERROR(std::logic_error, "Missing variable: BIAXIAL_COMPRESSION_MULTIPLIER", "");

	return 0;

	KRATOS_CATCH("");

}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CalculateMaterialResponse(const Vector& StrainVector,
	const Matrix& DeformationGradient,
	Vector& StressVector,
	Matrix& AlgorithmicTangent,
	const ProcessInfo& rCurrentProcessInfo,
	const Properties& rMaterialProperties,
	const GeometryType& rElementGeometry,
	const Vector& rShapeFunctionsValues,
	bool CalculateStresses,
	int CalculateTangent,
	bool SaveInternalVariables)
{
	ConstitutiveLaw::Parameters parameters(rElementGeometry, rMaterialProperties, rCurrentProcessInfo);
	Vector E(StrainVector);
	parameters.SetStrainVector( E );
	parameters.SetStressVector( StressVector );
	parameters.SetConstitutiveMatrix( AlgorithmicTangent );
	Flags& options = parameters.GetOptions();
	options.Set(ConstitutiveLaw::COMPUTE_STRESS, CalculateStresses);
	options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, CalculateTangent);
	double detF = 1.0;
	Matrix F(IdentityMatrix(2,2));
	parameters.SetDeterminantF(detF);
	parameters.SetDeformationGradientF(F);
	parameters.SetShapeFunctionsValues(rShapeFunctionsValues);
	this->CalculateMaterialResponseCauchy(parameters);
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::InitializeCalculationData(
	const Properties& props,
	const GeometryType& geom,
	const ProcessInfo& pinfo,
	CalculationData& data)
{
	// elasticity
	data.YoungModulus   			= props[YOUNG_MODULUS];
	data.PoissonRatio  				= props[POISSON_RATIO];
	this->CalculateElasticityMatrix(data);

	// Tension Damage Properties
	data.YieldStressTension 		= props[YIELD_STRESS_TENSION];
	data.FractureEnergyTension 		= props[FRACTURE_ENERGY_TENSION];

	// Compression Damage Properties
	data.DamageOnsetStressCompression 	= props[DAMAGE_ONSET_STRESS_COMPRESSION];
	data.YieldStressCompression 		= props[YIELD_STRESS_COMPRESSION];
	data.ResidualStressCompression 		= props[RESIDUAL_STRESS_COMPRESSION];
	data.YieldStrainCompression  		= props[YIELD_STRAIN_COMPRESSION];
	data.BezierControllerC1  			= props.Has(BEZIER_CONTROLLER_C1) ? props[BEZIER_CONTROLLER_C1] : 0.65;
	data.BezierControllerC2  			= props.Has(BEZIER_CONTROLLER_C2) ? props[BEZIER_CONTROLLER_C2] : 0.50;
	data.BezierControllerC3  			= props.Has(BEZIER_CONTROLLER_C3) ? props[BEZIER_CONTROLLER_C3] : 1.50;
	data.FractureEnergyCompression  	= props[FRACTURE_ENERGY_COMPRESSION];
	data.BiaxialCompressionMultiplier  	= props[BIAXIAL_COMPRESSION_MULTIPLIER];
	data.ShearCompressionReductor  		= props.Has(SHEAR_COMPRESSION_REDUCTOR) ? props[SHEAR_COMPRESSION_REDUCTOR] : 0.5;
	data.ShearCompressionReductor  		= std::min(std::max(data.ShearCompressionReductor,0.0),1.0);

	// Effective Stress Data
	data.EffectiveStressVector.resize(3,false);
	data.PrincipalStressVector.resize(2,false);
	data.EffectiveTensionStressVector.resize(3,false);
	data.EffectiveCompressionStressVector.resize(3,false);
	data.ProjectionTensorTension.resize(3,3,false);
	data.ProjectionTensorCompression.resize(3,3,false);

	// Misc
	data.CharacteristicLength = InitialCharacteristicLength;
	data.DeltaTime = pinfo[DELTA_TIME];
	data.TensionYieldModel = props.Has(TENSION_YIELD_MODEL) ? props[TENSION_YIELD_MODEL] : 0;
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CalculateElasticityMatrix(
	CalculationData& data)
{
	if(data.ElasticityMatrix.size1() != 3 || data.ElasticityMatrix.size2() != 3)
		data.ElasticityMatrix.resize(3,3,false);

	double c1 = data.YoungModulus / (1.0 - data.PoissonRatio * data.PoissonRatio);
	double c2 = c1 * data.PoissonRatio;
	double c3 = c1 * (1.0 - data.PoissonRatio) / 2.0;

	data.ElasticityMatrix(0,0) = c1;	data.ElasticityMatrix(0,1) = c2;	data.ElasticityMatrix(0,2) = 0.0;
	data.ElasticityMatrix(1,0) = c2;	data.ElasticityMatrix(1,1) = c1;	data.ElasticityMatrix(1,2) = 0.0;
	data.ElasticityMatrix(2,0) = 0.0;	data.ElasticityMatrix(2,1) = 0.0;	data.ElasticityMatrix(2,2) = c3;
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::TensionCompressionSplit(
	CalculationData& data)
{
	const array_1d<double,3>& effective_stress_vector 		= data.EffectiveStressVector;
	array_1d<double,2>& principal_stress_vector 			= data.PrincipalStressVector;
	array_1d<double,3>& effective_tension_stress_vector 	= data.EffectiveTensionStressVector;
	array_1d<double,3>& effective_compression_stress_vector = data.EffectiveCompressionStressVector;

	AdvancedConstitutiveLawUtilities<3>::CalculatePrincipalStresses(
		principal_stress_vector, effective_stress_vector);
	AdvancedConstitutiveLawUtilities<3>::SpectralDecomposition(
		effective_stress_vector, effective_tension_stress_vector, effective_compression_stress_vector);
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::ConstructProjectionTensors(
	CalculationData& data)
{
	Matrix& projection_tensor_tension 	  = data.ProjectionTensorTension;
	Matrix& projection_tensor_compression = data.ProjectionTensorCompression;

	const array_1d<double,3>& effective_stress_vector = data.EffectiveStressVector;

	Matrix effective_stress_tensor = MathUtils<double>::StressVectorToTensor(effective_stress_vector);
    BoundedMatrix<double, Dimension, Dimension> eigen_vectors_matrix;
    BoundedMatrix<double, Dimension, Dimension> eigen_values_matrix;

	MathUtils<double>::GaussSeidelEigenSystem(effective_stress_tensor, eigen_vectors_matrix, eigen_values_matrix, 1.0e-16, 20);

	array_1d<double,2> eigen_vector_1;
	array_1d<double,2> eigen_vector_2;

	for (IndexType i = 0; i < Dimension; ++i) {
		eigen_vector_1[i] = eigen_vectors_matrix(0, i);
	}
	for (IndexType i = 0; i < Dimension; ++i) {
		eigen_vector_2[i] = eigen_vectors_matrix(1, i);
	}

	array_1d<double,3> projection_vector_11;
	Matrix projection_tensor_11;
	projection_tensor_11 = outer_prod(eigen_vector_1, eigen_vector_1);
	projection_vector_11 = MathUtils<double>::StressTensorToVector(projection_tensor_11);

	array_1d<double,3> projection_vector_22;
	Matrix projection_tensor_22;
	projection_tensor_22 = outer_prod(eigen_vector_2, eigen_vector_2);
	projection_vector_22 = MathUtils<double>::StressTensorToVector(projection_tensor_22);

	projection_tensor_tension.clear();
	noalias(projection_tensor_tension) += HEAVISIDE(eigen_values_matrix(0, 0)) *
										  outer_prod(projection_vector_11, projection_vector_11);
	noalias(projection_tensor_tension) += HEAVISIDE(eigen_values_matrix(1, 1)) *
										  outer_prod(projection_vector_22, projection_vector_22);

#ifdef PROJECTION_OPERATOR_CERVERA_2003
/*
Theory from: "Viscoelasticity and rate-dependent continuum damage models"
			  Miguel Cervera
			  2003 (page 58)
*/
	array_1d<double,3> projection_vector_12;
	Matrix projection_tensor_12;
	array_1d<double,3> projection_vector_21;
	Matrix projection_tensor_21;
	projection_tensor_12 = outer_prod(eigen_vector_1, eigen_vector_2);
	projection_tensor_21 = outer_prod(eigen_vector_2, eigen_vector_1);

	array_1d<double,3> projection_vector_cross;
	Matrix projection_tensor_cross;
	projection_tensor_cross = 0.5 * (projection_tensor_12 + projection_tensor_21);
	projection_vector_cross = MathUtils<double>::StressTensorToVector(projection_tensor_cross);

	double factor_12;
	factor_12 = MACAULAY(eigen_values_matrix(0, 0)) - MACAULAY(eigen_values_matrix(1, 1));

	if (std::abs(eigen_values_matrix(0, 0) - eigen_values_matrix(1, 1)) > 0.0) {
		factor_12 *= 2.0;
		factor_12 /= eigen_values_matrix(0, 0) - eigen_values_matrix(1, 1);
	} else {
		factor_12 = 1.0;
	}
	noalias(projection_tensor_tension) += factor_12 * outer_prod(projection_vector_cross, projection_vector_cross);
#endif //PROJECTION_OPERATOR_CERVERA_2003

#ifdef PROJECTION_OPERATOR_CERVERA_2017
/*
Theory from: 	"An Energy-Equivalent d+/d- Damage Model with
				Enhanced Microcrack Closure-Reopening
				Capabilities for Cohesive-Frictional Materials"
				Miguel Cervera & Claudia Tesei
				2017 (page 7/30)
*/
	array_1d<double,3> projection_vector_12;
	Matrix projection_tensor_12;
	array_1d<double,3> projection_vector_21;
	Matrix projection_tensor_21;
	projection_tensor_12 = outer_prod(eigen_vector_1, eigen_vector_2);
	projection_tensor_21 = outer_prod(eigen_vector_2, eigen_vector_1);

	array_1d<double,3> projection_vector_cross;
	Matrix projection_tensor_cross;
	projection_tensor_cross = 0.5 * (projection_tensor_12 + projection_tensor_21);
	projection_vector_cross = MathUtils<double>::StressTensorToVector(projection_tensor_cross);

	noalias(projection_tensor_tension) += (HEAVISIDE(eigen_values_matrix(0, 0)) + HEAVISIDE(eigen_values_matrix(1, 1))) *
										  outer_prod(projection_vector_cross, projection_vector_cross);
#endif //PROJECTION_OPERATOR_CERVERA_2017

	noalias(projection_tensor_compression) = IdentityMatrix(3,3) - projection_tensor_tension;
}

/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CalculateEquivalentStressTension(CalculationData& data, double& UniaxialStressTension)
{
	UniaxialStressTension = 0.0;
	if(data.PrincipalStressVector(0) > 0.0){
		if (data.TensionYieldModel == 0) {
			// Lubliner Yield Criteria
			const double yield_compression 		 	= 	data.YieldStressCompression;
			const double yield_tension 				= 	data.YieldStressTension;
			const double alpha 						= 	(data.BiaxialCompressionMultiplier - 1.0) /
														(2.0* data.BiaxialCompressionMultiplier - 1.0);
			double I1, J2;
			array_1d<double, 3> deviator = ZeroVector(3);

			AdvancedConstitutiveLawUtilities<3>::CalculateI1Invariant(data.EffectiveStressVector, I1);
			AdvancedConstitutiveLawUtilities<3>::CalculateJ2Invariant(data.EffectiveStressVector, I1, deviator, J2);

			const double beta 	= yield_compression / yield_tension * (1.0 - alpha) - (1.0 + alpha);
			const double smax 	= std::max(std::max(data.PrincipalStressVector(0), data.PrincipalStressVector(1)),0.0);

			UniaxialStressTension = 1.0 / (1.0-alpha) * (alpha * I1 + std::sqrt(3.0 * J2) + beta * smax) /
									yield_compression * yield_tension;
		}
		else if (data.TensionYieldModel == 1) {
			// Rankine Yield Criteria
			UniaxialStressTension = std::max(std::max(data.PrincipalStressVector(0), data.PrincipalStressVector(1)), 0.0);
		}
	}
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CalculateEquivalentStressCompression(CalculationData& data, double& UniaxialStressCompression)
{
	UniaxialStressCompression = 0.0;
	if(data.PrincipalStressVector(1) < 0.0){
		const double yield_compression 			= 	data.DamageOnsetStressCompression;
		const double yield_tension 				= 	data.YieldStressTension;
		const double alpha 						= 	(data.BiaxialCompressionMultiplier - 1.0) /
													(2.0* data.BiaxialCompressionMultiplier - 1.0);
		double I1, J2;
		array_1d<double, 3> deviator = ZeroVector(3);

		AdvancedConstitutiveLawUtilities<3>::CalculateI1Invariant(data.EffectiveStressVector, I1);
		AdvancedConstitutiveLawUtilities<3>::CalculateJ2Invariant(data.EffectiveStressVector, I1, deviator, J2);

		const double beta = (yield_compression / yield_tension) * (1.0 - alpha) - (1.0 + alpha);
		const double smax = std::max(std::max(data.PrincipalStressVector(0), data.PrincipalStressVector(1)),0.0);

		UniaxialStressCompression = 1.0 / (1.0-alpha) * (alpha * I1 + std::sqrt(3.0 * J2) +
															data.ShearCompressionReductor * beta * smax);
	}
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CalculateDamageTension(
	CalculationData& data,
	double internal_variable,
	double& rDamage)
{
	if(internal_variable <= data.YieldStressTension) {
		rDamage = 0.0;
	}
	else {
		const double characteristic_length 		= 	data.CharacteristicLength;
		const double young_modulus   			= 	data.YoungModulus;
		const double yield_tension  			= 	data.YieldStressTension;
		const double fracture_energy_tension  	=  	data.FractureEnergyTension;
		const double initial_internal_variable  =  	yield_tension;
		const double material_length = 2.0 * young_modulus * fracture_energy_tension /
									   (yield_tension * yield_tension);

		if(characteristic_length >= material_length){
			std::stringstream ss;
			ss << "FRACTURE_ENERGY_TENSION is too low:  2*E*Gt/(ft*ft) = " << material_length
				<< ",   Characteristic Length = " << characteristic_length
				<< ",   FRACTURE_ENERGY_TENSION should be at least = " << (characteristic_length * yield_tension * yield_tension) / (2.0 * young_modulus) <<std::endl;
			std::cout << ss.str();
			exit(-1);
		}

		const double damage_parameter  = 2.0 * characteristic_length /
											(material_length - characteristic_length);

		rDamage = 	1.0 - initial_internal_variable / internal_variable *
							std::exp(damage_parameter *
							(1.0 - internal_variable / initial_internal_variable));

		const double internal_variable_min = 1.0e-2 * yield_tension;
		if ((1.0 - rDamage) * internal_variable < internal_variable_min) {
			rDamage = 1.0- internal_variable_min / internal_variable;
		}
	}
}

/***********************************************************************************/
/***********************************************************************************/

void DamageDPlusDMinusMasonry2DLaw::CalculateDamageCompression(
	CalculationData& data,
	double internal_variable,
	double& rDamage
	)
{
	if (internal_variable <= data.DamageOnsetStressCompression) {
		rDamage = 0.0;
	}
	else {
		// extract material parameters
		const double c1 = data.BezierControllerC1;
		const double c2 = data.BezierControllerC2;
		const double c3 = data.BezierControllerC3;

		const double specific_fracture_energy = data.FractureEnergyCompression / data.CharacteristicLength;
		const double young_modulus = data.YoungModulus;
		const double s_0 = data.DamageOnsetStressCompression;
		const double s_p = data.YieldStressCompression;
		const double s_i = s_p;
		const double s_j = s_p;
		const double s_u = std::abs(data.ResidualStressCompression);
		const double s_r = s_u;
		const double s_k = s_r + (s_p - s_r) * c1;

		const double e_0 = s_0 / young_modulus;
		double e_p = std::abs(data.YieldStrainCompression);
		const double alpha  = 2.0 * (e_p - s_p / young_modulus);
		double e_j = e_p + alpha * c2;
		double e_k = e_j + alpha * (1.0 - c2);
		double e_r = (e_k - e_j) / (s_p - s_k) * (s_p - s_r) + e_j; // 0.035;
		double e_u = e_r * c3;
		double e_i = s_p / young_modulus;

		const double Gc1 = 0.5 * s_p * e_p;
		const double Gc2 = (EvaluateBezierArea(e_p, e_j, e_k, s_p, s_j, s_k));
		const double Gc3 = (EvaluateBezierArea(e_k, e_r, e_u, s_k, s_r, s_u));
		const double bezier_fracture_energy = Gc1 + Gc2 + Gc3;

		const double stretcher = (specific_fracture_energy - Gc1) / (bezier_fracture_energy - Gc1) - 1.0;

		if (stretcher <= -1.0){
			std::stringstream ss;
			ss << "\n FRACTURE_ENERGY_COMPRESSION is too low" << std::endl;
			ss << "C1 = " << data.BezierControllerC1 <<std::endl;
			ss << "C2 = " << data.BezierControllerC2 <<std::endl;
			ss << "C3 = " << data.BezierControllerC3 <<std::endl;
			ss << "Characteristic Length = " << data.CharacteristicLength <<std::endl;
			ss << "Gc = " << data.FractureEnergyCompression <<std::endl;
			ss << "E = " << young_modulus << std::endl;
			ss << "yield compr = " << data.YieldStressCompression << std::endl;
			ss << "strain yield compr = " << data.YieldStrainCompression << std::endl;
			ss << "stretcher = " << stretcher << std::endl;
			ss << "e_0 = " << e_0 << std::endl;
			ss << "e_i = " << e_i << std::endl;
			ss << "e_p = " << e_p << std::endl;
			ss << "e_j = " << e_j << std::endl;
			ss << "e_k = " << e_k << std::endl;
			ss << "e_r = " << e_r << std::endl;
			ss << "e_u = " << e_u << std::endl;
			ss << "Gc1 = " << Gc1 << std::endl;
			ss << "Gc2 = " << Gc2 << std::endl;
			ss << "Gc3 = " << Gc3 << std::endl;
			ss << "Total energy = " << bezier_fracture_energy << std::endl;
			std::cout << ss.str();
			exit(-1);
		}

		this->ApplyBezierStretcherToStrains(stretcher, e_p, e_j, e_k, e_r, e_u);

		// current abscissa
		const double strain_like_counterpart = internal_variable / young_modulus;

		// compute damage
		double damage_variable = internal_variable;
		if (strain_like_counterpart <= e_p) {
			this->EvaluateBezierCurve(damage_variable, strain_like_counterpart, e_0, e_i, e_p, s_0, s_i, s_p);
		} else if (strain_like_counterpart <= e_k) {
			this->EvaluateBezierCurve(damage_variable, strain_like_counterpart, e_p, e_j, e_k, s_p, s_j, s_k);
		} else if (strain_like_counterpart <= e_u) {
			this->EvaluateBezierCurve(damage_variable, strain_like_counterpart, e_k, e_r, e_u, s_k, s_r, s_u);
		} else {
			damage_variable = s_r;
		}
		rDamage = 1.0 - damage_variable / internal_variable;
		rDamage = (rDamage > 0.99999) ? 0.99999 : rDamage;
	}
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::ComputeBezierEnergy(
	double& rBezierEnergy,
	double& rBezierEnergy1,
	double s_p,
	double s_k,
	double s_r,
	double e_p,
	double e_j,
	double e_k,
	double e_r,
	double e_u
	)
{
	rBezierEnergy1 = e_p * s_p / 2.0;
	double bezier_energy_2 = this->EvaluateBezierArea(e_p, e_j, e_k, s_p, s_p, s_k);
	double bezier_energy_3 = this->EvaluateBezierArea(e_k, e_r, e_u, s_k, s_r, s_r);
	rBezierEnergy = rBezierEnergy1 + bezier_energy_2 + bezier_energy_3;
}
/***********************************************************************************/
/***********************************************************************************/
double DamageDPlusDMinusMasonry2DLaw::EvaluateBezierArea(
	double x1,
	double x2,
	double x3,
	double y1,
	double y2,
	double y3
	)
{
	return x2 * y1 / 3.0 +
		   x3 * y1 / 6.0 -
		   x2 * y3 / 3.0 +
		   x3 * y2 / 3.0 +
		   x3 * y3 / 2.0 -
		   x1 * (y1 / 2.0 + y2 / 3.0 + y3 / 6.0);
}

/***********************************************************************************/
/***********************************************************************************/

void DamageDPlusDMinusMasonry2DLaw::ApplyBezierStretcherToStrains(
	double stretcher,
	double e_p,
	double& e_j,
	double& e_k,
	double& e_r,
	double& e_u
	)
{
	e_j += (e_j - e_p) * stretcher;
	e_k += (e_k - e_p) * stretcher;
	e_r += (e_r - e_p) * stretcher;
	e_u += (e_u - e_p) * stretcher;
}

/***********************************************************************************/
/***********************************************************************************/

void DamageDPlusDMinusMasonry2DLaw::EvaluateBezierCurve(
	double& rDamageParameter,
	double xi,
	double x1,
	double x2,
	double x3,
	double y1,
	double y2,
	double y3
	)
{
	double bezier_law_param_A = x1 - 2.0 * x2 + x3;
	double bezier_law_param_B = 2.0 * (x2 - x1);
	double bezier_law_param_C = x1 - xi;

	if(std::abs(bezier_law_param_A) < 1.0E-12){
		x2 					= x2 + 1.0E-6 * (x3 - x1);
		bezier_law_param_A 	= x1 - 2.0 * x2 + x3;
		bezier_law_param_B 	= 2.0 * (x2 - x1);
		bezier_law_param_C 	= x1 - xi;
	}

	const double bezier_law_param_D = bezier_law_param_B * bezier_law_param_B -
								4.0 * bezier_law_param_A * bezier_law_param_C;
	const double bezier_law_param_t = (-bezier_law_param_B + std::sqrt(bezier_law_param_D)) /
								(2.0 * bezier_law_param_A);
	rDamageParameter = (y1 - 2.0 * y2 + y3) * bezier_law_param_t * bezier_law_param_t +
					   2.0 * (y2 - y1) * bezier_law_param_t + y1;
}

/***********************************************************************************/
/***********************************************************************************/

void DamageDPlusDMinusMasonry2DLaw::ComputeCharacteristicLength(
    const GeometryType& geom,
    double& rCharacteristicLength)
{
    rCharacteristicLength = geom.Length();

#ifdef OPTIMIZE_CHARACTERISTIC_LENGTH
    if(geom. WorkingSpaceDimension() == 2 && geom.PointsNumber() == 4){
        //2D Element with 4 Nodes
        double aX = (geom[0].X0() + geom[3].X0())/2.0;  //center_coord_X_Node1Node4
        double aY = (geom[0].Y0() + geom[3].Y0())/2.0;  //center_coord_Y_Node1Node4
        double bX = (geom[1].X0() + geom[2].X0())/2.0;  //center_coord_X_Node2Node3
        double bY = (geom[1].Y0() + geom[2].Y0())/2.0;  //center_coord_Y_Node2Node3
        double cX = (geom[0].X0() + geom[1].X0())/2.0;  //center_coord_X_Node1Node2
        double cY = (geom[0].Y0() + geom[1].Y0())/2.0;  //center_coord_Y_Node1Node2
        double dX = (geom[2].X0() + geom[3].X0())/2.0;  //center_coord_X_Node3Node4
        double dY = (geom[2].Y0() + geom[3].Y0())/2.0;  //center_coord_Y_Node3Node4

        double SabX = aX - bX;
        double SabY = aY - bY;
        double ScdX = cX - dX;
        double ScdY = cY - dY;

        double length_ab = std::sqrt(std::pow(SabX,2) + std::pow(SabY,2));
        double length_cd = std::sqrt(std::pow(ScdX,2) + std::pow(ScdY,2));

        rCharacteristicLength = std::min(length_ab, length_cd);
    }
#endif
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CalculateMaterialResponseInternal(
	const Vector& StrainVector,
	Vector& PredictiveStressVector,
	CalculationData& data,
	Properties props)
{
	if(PredictiveStressVector.size() != VoigtSize)
		PredictiveStressVector.resize(VoigtSize,false);

	ThresholdTension     = CurrentThresholdTension;
	ThresholdCompression = CurrentThresholdCompression;

	noalias(data.EffectiveStressVector) = prod(data.ElasticityMatrix, StrainVector);

	if(std::abs(data.EffectiveStressVector(0)) < tolerance) {data.EffectiveStressVector(0) = 0.0;}
	if(std::abs(data.EffectiveStressVector(1)) < tolerance) {data.EffectiveStressVector(1) = 0.0;}
	if(std::abs(data.EffectiveStressVector(2)) < tolerance) {data.EffectiveStressVector(2) = 0.0;}

	this->TensionCompressionSplit(data);
	this->ConstructProjectionTensors(data);

	// compute the equivalent stress measures
	this->CalculateEquivalentStressTension(data, UniaxialStressTension);

	this->CalculateEquivalentStressCompression(data, UniaxialStressCompression);

	// damage update
	if (props[INTEGRATION_IMPLEX] != 0){ //IMPLEX Integration
		// time factor
		double time_factor = 0.0;
		if(PreviousDeltaTime > 0.0) time_factor = data.DeltaTime / PreviousDeltaTime;
		CurrentDeltaTime = data.DeltaTime;

		// explicit evaluation
		ThresholdTension = CurrentThresholdTension + time_factor * (CurrentThresholdTension - PreviousThresholdTension);
		ThresholdCompression = CurrentThresholdCompression + time_factor * (CurrentThresholdCompression - PreviousThresholdCompression);

		// save implicit variables for the finalize_solution_step
		double implicit_threshold_tension 		= CurrentThresholdTension;
		double implicit_threshold_compression 	= CurrentThresholdCompression;

		if(UniaxialStressTension > implicit_threshold_tension)
			implicit_threshold_tension = UniaxialStressTension;

		if(UniaxialStressCompression > implicit_threshold_compression)
			implicit_threshold_compression = UniaxialStressCompression;

		TemporaryImplicitThresholdTension 		= implicit_threshold_tension;
		TemporaryImplicitThresholdTCompression 	= implicit_threshold_compression;

		// new damage variables (explicit)
		this->CalculateDamageTension(data, ThresholdTension, DamageParameterTension);
		this->CalculateDamageCompression(data, ThresholdCompression, DamageParameterCompression);
	}
	else { // IMPLICIT Integration
		if(UniaxialStressTension > ThresholdTension)
			ThresholdTension = UniaxialStressTension;
		this->CalculateDamageTension(data, ThresholdTension, DamageParameterTension);

		if(UniaxialStressCompression > ThresholdCompression)
			ThresholdCompression = UniaxialStressCompression;
		this->CalculateDamageCompression(data, ThresholdCompression, DamageParameterCompression);

		TemporaryImplicitThresholdTension = ThresholdTension;
		TemporaryImplicitThresholdTCompression = ThresholdCompression;
	}

	// calculation of stress tensor
	noalias(PredictiveStressVector)  = (1.0 - DamageParameterTension)     * data.EffectiveTensionStressVector;
	noalias(PredictiveStressVector) += (1.0 - DamageParameterCompression) * data.EffectiveCompressionStressVector;
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CheckDamageLoadingUnloading(
	bool& is_damaging_tension,
	bool& is_damaging_compression)
{
	const double F_tension 		= UniaxialStressTension - CurrentThresholdTension;
	const double F_compression 	= UniaxialStressCompression - CurrentThresholdCompression;

	is_damaging_tension = false;
	if (F_tension > 0.0){
		is_damaging_tension = true;
	}
	is_damaging_compression = false;
	if (F_compression > 0.0){
		is_damaging_compression = true;
	}
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CalculateTangentTensor(
	Parameters& rValues,
	Vector StrainVector,
	Vector PredictiveStressVector,
	CalculationData& data,
	const Properties& props)
{
	// prepare constitutive matrix
	Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
	if(constitutive_matrix.size1() != VoigtSize || constitutive_matrix.size2() != VoigtSize)
		constitutive_matrix.resize(VoigtSize, VoigtSize);

	// save internal variables
	double save_threshold_tension			= ThresholdTension;
	double save_threshold_compression		= ThresholdCompression;
	double save_damage_tension				= DamageParameterTension;
	double save_damage_compression			= DamageParameterCompression;
	double save_uniaxial_stress_tension		= UniaxialStressTension;
	double save_uniaxial_stress_compression	= UniaxialStressCompression;

	// perturbation parameter
	double perturbation_factor = 1.0E-8;

	// perturbed vectors
	Vector perturbated_strain_vector(VoigtSize);
	Vector perturbated_stress_vector(VoigtSize);

	// apply perturbation to each strain component...
	for(size_t j = 0; j < VoigtSize; j++)
	{
		noalias(perturbated_strain_vector) = StrainVector;

		perturbated_strain_vector(j) = StrainVector(j) + perturbation_factor;
		this->CalculateMaterialResponseInternal(perturbated_strain_vector, perturbated_stress_vector, data, props);

		for(size_t i = 0; i < VoigtSize; i++)
			constitutive_matrix(i,j) = ( perturbated_stress_vector(i) - PredictiveStressVector(i) ) /
										perturbation_factor;
	}

	// restore internal variables
	ThresholdTension	   		= save_threshold_tension;
	ThresholdCompression	   	= save_threshold_compression;
	DamageParameterTension 		= save_damage_tension;
	DamageParameterCompression 	= save_damage_compression;
	UniaxialStressTension 		= save_uniaxial_stress_tension;
	UniaxialStressCompression 	= save_uniaxial_stress_compression;
}
/***********************************************************************************/
/***********************************************************************************/
void DamageDPlusDMinusMasonry2DLaw::CalculateSecantTensor(
	Parameters& rValues,
	CalculationData& data)
{
	Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
	if(constitutive_matrix.size1() != VoigtSize || constitutive_matrix.size2() != VoigtSize)
		constitutive_matrix.resize(VoigtSize, VoigtSize);

	Matrix DamageMatrix( IdentityMatrix(3,3) );
	noalias(DamageMatrix) -= DamageParameterTension     * data.ProjectionTensorTension;
	noalias(DamageMatrix) -= DamageParameterCompression * data.ProjectionTensorCompression;

	noalias(constitutive_matrix) = prod(DamageMatrix, data.ElasticityMatrix);
}

} // namespace Kratos
