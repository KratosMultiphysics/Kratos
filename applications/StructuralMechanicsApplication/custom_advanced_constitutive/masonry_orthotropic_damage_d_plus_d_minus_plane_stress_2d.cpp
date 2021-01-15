// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Philip Kalkbrenner
//                   Massimo Petracca
//                   Alejandro Cornejo 
//  
//

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "masonry_orthotropic_damage_d_plus_d_minus_plane_stress_2d.h"
#include "includes/model_part.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/constitutive_law_utilities.h"

#define OPTIMIZE_CHARACTERISTIC_LENGTH
#define HEAVISIDE(X) ( X >= 0.0 ? 1.0 : 0.0)
#define MACAULAY(X)  ( X >= 0.0 ? X : 0.0)
#define PROJECTION_OPERATOR_CERVERA_2003
//#define PROJECTION_OPERATOR_CERVERA_2017

namespace Kratos
{
    /// Has variable <double>
    bool MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::Has(
        const Variable<double>& rThisVariable)
    {
        // TENSION
        if (rThisVariable == DAMAGE_TENSION)
            return true;
        if (rThisVariable == UNIAXIAL_STRESS_TENSION)
            return true;
        if (rThisVariable == THRESHOLD_TENSION)
            return true;
        // COMPRESSION
        if (rThisVariable == DAMAGE_COMPRESSION)
            return true;
        if (rThisVariable == UNIAXIAL_STRESS_COMPRESSION)
            return true;
        if (rThisVariable == THRESHOLD_COMPRESSION)
            return true;
        return false;
    }

    /// GetValue variable <double>
    double& MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::GetValue(
        const Variable<double>& rThisVariable,
        double& rValue)
    {
        rValue = 0.0;
        // TENSION
        if (rThisVariable == DAMAGE_TENSION)
            rValue = mDamageTension;
        else if (rThisVariable == UNIAXIAL_STRESS_TENSION)
            rValue = UniaxialStressTension;
        else if (rThisVariable == THRESHOLD_TENSION)
            rValue = mThresholdTension;
        // COMPRESSION
        else if (rThisVariable == DAMAGE_COMPRESSION)
            rValue = mDamageCompression;
        else if (rThisVariable == UNIAXIAL_STRESS_COMPRESSION)
            rValue = UniaxialStressCompression;
        else if (rThisVariable == THRESHOLD_COMPRESSION)
            rValue = mThresholdCompression;
        return rValue;
    }

    /// SetValue variable <double>
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::SetValue(
        const Variable<double>& rVariable,
        const double& rValue,
        const ProcessInfo& rCurrentProcessInfo)
    {
        // TENSION
        if (rVariable == DAMAGE_TENSION)
            mDamageTension = rValue;
        else if (rVariable == UNIAXIAL_STRESS_TENSION)
            UniaxialStressTension = rValue;
        else if (rVariable == THRESHOLD_TENSION)
            mThresholdTension = rValue;
        // COMPRESSION
        else if (rVariable == DAMAGE_COMPRESSION)
            mDamageCompression = rValue;
        else if (rVariable == UNIAXIAL_STRESS_COMPRESSION)
            UniaxialStressCompression = rValue;
        else if (rVariable == THRESHOLD_COMPRESSION)
            mThresholdCompression = rValue;
    }

    /***********************************************************************************/
    /***********************************************************************************/
    bool MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::ValidateInput(
        const Properties& rMaterialProperties)
    {
        // ELASTIC PARAMETERS
        if (!rMaterialProperties.Has(YOUNG_MODULUS_1))                        return false;
        if (!rMaterialProperties.Has(YOUNG_MODULUS_2))                        return false;
        if (!rMaterialProperties.Has(POISSON_RATIO_12))                       return false;
        if (!rMaterialProperties.Has(POISSON_RATIO_21))                       return false;
        // TENSION PARAMETERS
        if (!rMaterialProperties.Has(YIELD_STRESS_TENSION))                   return false;
        if( !rMaterialProperties.Has(FRACTURE_ENERGY_TENSION) )               return false;
        // COMPRESSION PARAMETERS
        if (!rMaterialProperties.Has(DAMAGE_ONSET_STRESS_COMPRESSION))        return false;
        if (!rMaterialProperties.Has(YIELD_STRESS_COMPRESSION))               return false;
        if (!rMaterialProperties.Has(YIELD_STRAIN_COMPRESSION))               return false;
        if (!rMaterialProperties.Has(RESIDUAL_STRESS_COMPRESSION))            return false;
        if( !rMaterialProperties.Has(FRACTURE_ENERGY_COMPRESSION) )           return false;
        // BIAXIAL PARAMETERS
        if (!rMaterialProperties.Has(BIAXIAL_COMPRESSION_MULTIPLIER))         return false;
        return true;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::StrainMeasure MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::GetStrainMeasure()
    {
        return ConstitutiveLaw::StrainMeasure_Infinitesimal;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::StressMeasure MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::GetStressMeasure()
    {
        return ConstitutiveLaw::StressMeasure_Cauchy;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    bool MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::IsIncremental()
    {
        return false;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::InitializeMaterial(
        const Properties& rProperties,
        const GeometryType& rGeometry,
        const Vector& rShapeFunctionsValues)
    {
        if (!InitializeDamageLaw) {
            mThresholdTension = rProperties[YIELD_STRESS_TENSION];
            mCurrentThresholdTension = mThresholdTension;

            mThresholdCompression = rProperties[DAMAGE_ONSET_STRESS_COMPRESSION];
            mCurrentThresholdCompression = mThresholdCompression;

            mDamageTension = 0.0;
            mDamageCompression = 0.0;
            UniaxialStressTension = 0.0;
            UniaxialStressCompression = 0.0;

            this->ComputeCharacteristicLength(
                rGeometry,
                mInitialCharacteristicLength);

            // Begin IMPLEX Integration - Only if switched on
            if (rProperties.Has(INTEGRATION_IMPLEX))
            {
                if (rProperties[INTEGRATION_IMPLEX] != 0) {
                    PreviousThresholdTension = mThresholdTension;
                    PreviousThresholdCompression = mThresholdCompression;
                    CurrentDeltaTime = 0.0;
                    PreviousDeltaTime = 0.0;
                }
            }
            // End IMPLEX Integration

            InitializeDamageLaw = true;
        }
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::FinalizeSolutionStep(
        const Properties& rProperties,
        const GeometryType& rGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo)
    {
        // Begin IMPLEX Integration - Only if switched on
        if (rProperties.Has(INTEGRATION_IMPLEX))
        {
            if (rProperties[INTEGRATION_IMPLEX] != 0) {
                mThresholdTension = TemporaryImplicitThresholdTension;
                mThresholdCompression = TemporaryImplicitThresholdTCompression;

                // move from n to n-1
                PreviousThresholdTension = mCurrentThresholdTension;
                PreviousThresholdCompression = mCurrentThresholdCompression;
                PreviousDeltaTime = CurrentDeltaTime;
            }
        }
        // End IMPLEX Integration


        // save converged values
        mCurrentThresholdTension = mThresholdTension;
        mCurrentThresholdCompression = mThresholdCompression;

    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::CalculateMaterialResponsePK1(
        Parameters& rValues)
    {
        CalculateMaterialResponseCauchy(rValues);
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::CalculateMaterialResponsePK2(
        Parameters& rValues)
    {
        CalculateMaterialResponseCauchy(rValues);
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::CalculateMaterialResponseKirchhoff(
        Parameters& rValues)
    {
        CalculateMaterialResponseCauchy(rValues);
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::CalculateMaterialResponseCauchy(
        Parameters& rParameters)
    {
        const ProcessInfo& r_process_info = rParameters.GetProcessInfo();
        const GeometryType& r_geometry = rParameters.GetElementGeometry();
        const Properties& r_properties = rParameters.GetMaterialProperties();

        const Vector& r_strain_vector = rParameters.GetStrainVector();
        Vector& PredictiveStressVector = rParameters.GetStressVector();

        CalculationData data;
        this->InitializeCalculationData(r_properties, r_geometry, r_process_info, data);

        // Transform the FEM StrainVector to the isotropic StrainVector 
        Vector StrainVectorIso = prod(data.TransformationMatrix, r_strain_vector);

        this->CalculateMaterialResponseInternal(StrainVectorIso, PredictiveStressVector, data, 0);

        bool is_damaging_tension = false;
        bool is_damaging_compression = false;
        this->CheckDamageLoadingUnloading(is_damaging_tension, is_damaging_compression);


        // Computation of the Constitutive Tensor
        if (rParameters.GetOptions().Is(COMPUTE_CONSTITUTIVE_TENSOR)) {
            if (is_damaging_tension || is_damaging_compression) {
                this->CalculateTangentTensor(rParameters, StrainVectorIso, PredictiveStressVector, data, r_properties);
            }
            else {
                this->CalculateSecantTensor(rParameters, data);
            }
        }
        // Back rotation of final Predictive Stress
        noalias(PredictiveStressVector) = prod(data.TransformationMatrixTranspose, PredictiveStressVector);
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::ResetMaterial(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues)
    {
        mThresholdTension = 0.0;
        mCurrentThresholdTension = 0.0;

        mThresholdCompression = 0.0;
        mCurrentThresholdCompression = 0.0;

        mDamageTension = 0.0;
        mDamageCompression = 0.0;

        mInitialCharacteristicLength = 0.0;
        InitializeDamageLaw = false;

    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::GetLawFeatures(
        Features& rFeatures)
    {
        //Set the type of law
        rFeatures.mOptions.Set(PLANE_STRESS_LAW);
        rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
        rFeatures.mOptions.Set(ISOTROPIC);

        //Set strain measure required by the consitutive law
        rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);

        //Set the strain size
        rFeatures.mStrainSize = GetStrainSize();

        //Set the space dimension
        rFeatures.mSpaceDimension = WorkingSpaceDimension();
    }
    /***********************************************************************************/
    /***********************************************************************************/
    int MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

            if (!rMaterialProperties.Has(YOUNG_MODULUS))
                KRATOS_THROW_ERROR(std::logic_error, "Missing variable: YOUNG_MODULUS", "");

        if (!rMaterialProperties.Has(POISSON_RATIO))
            KRATOS_THROW_ERROR(std::logic_error, "Missing variable: POISSON_RATIO", "");

        ////if (!rMaterialProperties.Has(TRANSFORMATION_MATRIX_COMP_11))
        ////    KRATOS_THROW_ERROR(std::logic_error, "Missing variable: TRANSFORMATION_MATRIX_COMP_11", "");

        ////if (!rMaterialProperties.Has(TRANSFORMATION_MATRIX_COMP_12))
        ////    KRATOS_THROW_ERROR(std::logic_error, "Missing variable: TRANSFORMATION_MATRIX_COMP_12", "");

        ////if (!rMaterialProperties.Has(TRANSFORMATION_MATRIX_COMP_13))
        ////    KRATOS_THROW_ERROR(std::logic_error, "Missing variable: TRANSFORMATION_MATRIX_COMP_13", "");

        ////if (!rMaterialProperties.Has(TRANSFORMATION_MATRIX_COMP_21))
        ////    KRATOS_THROW_ERROR(std::logic_error, "Missing variable: TRANSFORMATION_MATRIX_COMP_21", "");

        ////if (!rMaterialProperties.Has(TRANSFORMATION_MATRIX_COMP_22))
        ////    KRATOS_THROW_ERROR(std::logic_error, "Missing variable: TRANSFORMATION_MATRIX_COMP_22", "");

        ////if (!rMaterialProperties.Has(TRANSFORMATION_MATRIX_COMP_23))
        ////    KRATOS_THROW_ERROR(std::logic_error, "Missing variable: TRANSFORMATION_MATRIX_COMP_23", "");

        ////if (!rMaterialProperties.Has(TRANSFORMATION_MATRIX_COMP_31))
        ////    KRATOS_THROW_ERROR(std::logic_error, "Missing variable: TRANSFORMATION_MATRIX_COMP_31", "");

        ////if (!rMaterialProperties.Has(TRANSFORMATION_MATRIX_COMP_32))
        ////    KRATOS_THROW_ERROR(std::logic_error, "Missing variable: TRANSFORMATION_MATRIX_COMP_32", "");

        ////if (!rMaterialProperties.Has(TRANSFORMATION_MATRIX_COMP_33))
        ////    KRATOS_THROW_ERROR(std::logic_error, "Missing variable: TRANSFORMATION_MATRIX_COMP_33", "");

        if (!rMaterialProperties.Has(YIELD_STRESS_TENSION))
            KRATOS_THROW_ERROR(std::logic_error, "Missing variable: YIELD_STRESS_TENSION", "");

        ////if (!rMaterialProperties.Has(SPECIFIC_FRACTURE_ENERGY_TENSION))
        ////    KRATOS_THROW_ERROR(std::logic_error, "Missing variable: SPECIFIC_FRACTURE_ENERGY_TENSION", "");

        if( !rMaterialProperties.Has(FRACTURE_ENERGY_TENSION) )
            KRATOS_THROW_ERROR(std::logic_error, "Missing variable: FRACTURE_ENERGY_TENSION", "");

        if (!rMaterialProperties.Has(DAMAGE_ONSET_STRESS_COMPRESSION))
            KRATOS_THROW_ERROR(std::logic_error, "Missing variable: DAMAGE_ONSET_STRESS_COMPRESSION", "");

        if (!rMaterialProperties.Has(YIELD_STRESS_COMPRESSION))
            KRATOS_THROW_ERROR(std::logic_error, "Missing variable: YIELD_STRESS_COMPRESSION", "");

        if (!rMaterialProperties.Has(RESIDUAL_STRESS_COMPRESSION))
            KRATOS_THROW_ERROR(std::logic_error, "Missing variable: RESIDUAL_STRESS_COMPRESSION", "");

        if (!rMaterialProperties.Has(YIELD_STRAIN_COMPRESSION))
            KRATOS_THROW_ERROR(std::logic_error, "Missing variable: YIELD_STRAIN_COMPRESSION", "");

        //if( !rMaterialProperties.Has(FRACTURE_ENERGY_COMPRESSION) )
        //    KRATOS_THROW_ERROR(std::logic_error, "Missing variable: FRACTURE_ENERGY_COMPRESSION", "");

        if (!rMaterialProperties.Has(BIAXIAL_COMPRESSION_MULTIPLIER))
            KRATOS_THROW_ERROR(std::logic_error, "Missing variable: BIAXIAL_COMPRESSION_MULTIPLIER", "");

        return 0;

        KRATOS_CATCH("");
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::CalculateMaterialResponse(const Vector& StrainVector,
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
        parameters.SetStrainVector(E);
        parameters.SetStressVector(StressVector);
        parameters.SetConstitutiveMatrix(AlgorithmicTangent);
        Flags& options = parameters.GetOptions();
        options.Set(ConstitutiveLaw::COMPUTE_STRESS, CalculateStresses);
        options.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, CalculateTangent);
        double detF = 1.0;
        Matrix F(IdentityMatrix(2, 2));
        parameters.SetDeterminantF(detF);
        parameters.SetDeformationGradientF(F);
        parameters.SetShapeFunctionsValues(rShapeFunctionsValues);
        this->CalculateMaterialResponseCauchy(parameters);
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::InitializeCalculationData(
        const Properties& rProperties,
        const GeometryType& rGeometry,
        const ProcessInfo& pinfo,
        CalculationData& data)
    {
        this->CalculateOrthotropicElasticityMatrix(
            data.ElasticityMatrix,
            rProperties[YOUNG_MODULUS_1],
            rProperties[YOUNG_MODULUS_2],
            rProperties[POISSON_RATIO_12],
            rProperties[POISSON_RATIO_21],
            rProperties[SHEAR_MODULUS]);

        // transformation
        ////data.TransformationMatrixComp11 = rProperties[TRANSFORMATION_MATRIX_COMP_11];
        ////data.TransformationMatrixComp12 = rProperties[TRANSFORMATION_MATRIX_COMP_12];
        ////data.TransformationMatrixComp13 = rProperties[TRANSFORMATION_MATRIX_COMP_13];
        ////data.TransformationMatrixComp21 = rProperties[TRANSFORMATION_MATRIX_COMP_21];
        ////data.TransformationMatrixComp22 = rProperties[TRANSFORMATION_MATRIX_COMP_22];
        ////data.TransformationMatrixComp23 = rProperties[TRANSFORMATION_MATRIX_COMP_23];
        ////data.TransformationMatrixComp31 = rProperties[TRANSFORMATION_MATRIX_COMP_31];
        ////data.TransformationMatrixComp32 = rProperties[TRANSFORMATION_MATRIX_COMP_32];
        ////data.TransformationMatrixComp33 = rProperties[TRANSFORMATION_MATRIX_COMP_33];

        this->AssambleTransformationMatrix(data);

        // Tension Damage Properties
        data.YieldStressTension = rProperties[YIELD_STRESS_TENSION];
        ////data.SpecificFractureEnergyTension = rProperties[SPECIFIC_FRACTURE_ENERGY_TENSION];
        data.FractureEnergyTension = rProperties[FRACTURE_ENERGY_TENSION];

        // Compression Damage Properties
        data.DamageOnsetStressCompression = rProperties[DAMAGE_ONSET_STRESS_COMPRESSION];
        data.YieldStressCompression = rProperties[YIELD_STRESS_COMPRESSION];
        data.ResidualStressCompression = rProperties[RESIDUAL_STRESS_COMPRESSION];
        data.YieldStrainCompression = rProperties[YIELD_STRAIN_COMPRESSION];
        data.FractureEnergyCompression      = rProperties[FRACTURE_ENERGY_COMPRESSION];

        data.BezierControllerC1 = rProperties.Has(BEZIER_CONTROLLER_C1) ? rProperties[BEZIER_CONTROLLER_C1] : 0.65;
        data.BezierControllerC2 = rProperties.Has(BEZIER_CONTROLLER_C2) ? rProperties[BEZIER_CONTROLLER_C2] : 0.50;
        data.BezierControllerC3 = rProperties.Has(BEZIER_CONTROLLER_C3) ? rProperties[BEZIER_CONTROLLER_C3] : 1.50;
        data.BiaxialCompressionMultiplier = rProperties[BIAXIAL_COMPRESSION_MULTIPLIER];
        data.ShearCompressionReductor = rProperties.Has(SHEAR_COMPRESSION_REDUCTOR) ? rProperties[SHEAR_COMPRESSION_REDUCTOR] : 0.5;
        data.ShearCompressionReductor = std::min(std::max(data.ShearCompressionReductor, 0.0), 1.0);

        // Effective Stress Data
        data.EffectiveStressVector.resize(3, false);
        data.PrincipalStressVector.resize(2, false);
        data.EffectiveTensionStressVector.resize(3, false);
        data.EffectiveCompressionStressVector.resize(3, false);
        data.ProjectionTensorTension.resize(3, 3, false);
        data.ProjectionTensorCompression.resize(3, 3, false);

        // Misc
        this->ComputeCharacteristicLength(
            rGeometry,
            data.CharacteristicLength);
        data.DeltaTime = pinfo[DELTA_TIME];
        data.TensionYieldModel = rProperties.Has(TENSION_YIELD_MODEL) ? rProperties[TENSION_YIELD_MODEL] : 0;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::CalculateOrthotropicElasticityMatrix(
        Matrix& rC,
        double E_1,
        double E_2,
        double nu_12,
        double nu_21,
        double G_12)
    {
        if (rC.size1() != 3 || rC.size2() != 3)
            rC.resize(3, 3, false);
        rC = ZeroMatrix(3, 3);

        const double nu_12_nu_21 = 1.0 - nu_12 * nu_21;

        rC(0, 0) = E_1 / nu_12_nu_21; rC(0, 1) = nu_12 * rC(0, 0);
        rC(1, 1) = E_2 / nu_12_nu_21; rC(1, 0) = nu_21 * rC(1, 1);
        rC(2, 2) = G_12;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::AssambleTransformationMatrix(
        CalculationData& data)
    {
        if (data.TransformationMatrix.size1() != 3 || data.TransformationMatrix.size2() != 3)
            data.TransformationMatrix.resize(3, 3, false);

        data.TransformationMatrix(0, 0) = 1; data.TransformationMatrix(0, 1) = 0; data.TransformationMatrix(0, 2) = 0;
        data.TransformationMatrix(1, 0) = 0; data.TransformationMatrix(1, 1) = 1; data.TransformationMatrix(1, 2) = 0;
        data.TransformationMatrix(2, 0) = 0; data.TransformationMatrix(2, 1) = 0; data.TransformationMatrix(2, 2) = 1;

        if (data.TransformationMatrixTranspose.size1() != 3 || data.TransformationMatrixTranspose.size2() != 3)
            data.TransformationMatrixTranspose.resize(3, 3, false);

        data.TransformationMatrixTranspose(0, 0) = 1;        data.TransformationMatrixTranspose(0, 1) = 0;        data.TransformationMatrixTranspose(0, 2) = 0;
        data.TransformationMatrixTranspose(1, 0) = 0;        data.TransformationMatrixTranspose(1, 1) = 1;        data.TransformationMatrixTranspose(1, 2) = 0;
        data.TransformationMatrixTranspose(2, 0) = 0;        data.TransformationMatrixTranspose(2, 1) = 0;        data.TransformationMatrixTranspose(2, 2) = 1;
    }
    /***********************************************************************************/
    /***********************************************************************************/

    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::TensionCompressionSplit(
        CalculationData& data)
    {
        const array_1d<double, 3>& effective_stress_vector = data.EffectiveStressVector;
        array_1d<double, 2>& principal_stress_vector = data.PrincipalStressVector;
        array_1d<double, 3>& effective_tension_stress_vector = data.EffectiveTensionStressVector;
        array_1d<double, 3>& effective_compression_stress_vector = data.EffectiveCompressionStressVector;

        ConstitutiveLawUtilities<3>::CalculatePrincipalStresses(
            principal_stress_vector, effective_stress_vector);
        ConstitutiveLawUtilities<3>::SpectralDecomposition(
            effective_stress_vector, effective_tension_stress_vector, effective_compression_stress_vector);
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::ConstructProjectionTensors(
        CalculationData& data)
    {
        Matrix& projection_tensor_tension = data.ProjectionTensorTension;
        Matrix& projection_tensor_compression = data.ProjectionTensorCompression;

        const array_1d<double, 3>& effective_stress_vector = data.EffectiveStressVector;

        Matrix effective_stress_tensor = MathUtils<double>::StressVectorToTensor(effective_stress_vector);
        BoundedMatrix<double, Dimension, Dimension> eigen_vectors_matrix;
        BoundedMatrix<double, Dimension, Dimension> eigen_values_matrix;

        MathUtils<double>::GaussSeidelEigenSystem(effective_stress_tensor, eigen_vectors_matrix, eigen_values_matrix, 1.0e-16, 20);

        array_1d<double, 2> eigen_vector_1;
        array_1d<double, 2> eigen_vector_2;

        for (IndexType i = 0; i < Dimension; ++i) {
            eigen_vector_1[i] = eigen_vectors_matrix(0, i);
        }
        for (IndexType i = 0; i < Dimension; ++i) {
            eigen_vector_2[i] = eigen_vectors_matrix(1, i);
        }

        array_1d<double, 3> projection_vector_11;
        Matrix projection_tensor_11;
        projection_tensor_11 = outer_prod(eigen_vector_1, eigen_vector_1);
        projection_vector_11 = MathUtils<double>::StressTensorToVector(projection_tensor_11);

        array_1d<double, 3> projection_vector_22;
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
        array_1d<double, 3> projection_vector_12;
        Matrix projection_tensor_12;
        array_1d<double, 3> projection_vector_21;
        Matrix projection_tensor_21;
        projection_tensor_12 = outer_prod(eigen_vector_1, eigen_vector_2);
        projection_tensor_21 = outer_prod(eigen_vector_2, eigen_vector_1);

        array_1d<double, 3> projection_vector_cross;
        Matrix projection_tensor_cross;
        projection_tensor_cross = 0.5 * (projection_tensor_12 + projection_tensor_21);
        projection_vector_cross = MathUtils<double>::StressTensorToVector(projection_tensor_cross);

        double factor_12;
        factor_12 = MACAULAY(eigen_values_matrix(0, 0)) - MACAULAY(eigen_values_matrix(1, 1));

        if (std::abs(eigen_values_matrix(0, 0) - eigen_values_matrix(1, 1)) > 0.0) {
            factor_12 *= 2.0;
            factor_12 /= eigen_values_matrix(0, 0) - eigen_values_matrix(1, 1);
        }
        else {
            factor_12 = 1.0;
        }
        noalias(projection_tensor_tension) += factor_12 * outer_prod(projection_vector_cross, projection_vector_cross);
#endif //PROJECTION_OPERATOR_CERVERA_2003

#ifdef PROJECTION_OPERATOR_CERVERA_2017
        /*
        Theory from:     "An Energy-Equivalent d+/d- Damage Model with
                        Enhanced Microcrack Closure-Reopening
                        Capabilities for Cohesive-Frictional Materials"
                        Miguel Cervera & Claudia Tesei
                        2017 (page 7/30)
        */
        array_1d<double, 3> projection_vector_12;
        Matrix projection_tensor_12;
        array_1d<double, 3> projection_vector_21;
        Matrix projection_tensor_21;
        projection_tensor_12 = outer_prod(eigen_vector_1, eigen_vector_2);
        projection_tensor_21 = outer_prod(eigen_vector_2, eigen_vector_1);

        array_1d<double, 3> projection_vector_cross;
        Matrix projection_tensor_cross;
        projection_tensor_cross = 0.5 * (projection_tensor_12 + projection_tensor_21);
        projection_vector_cross = MathUtils<double>::StressTensorToVector(projection_tensor_cross);

        noalias(projection_tensor_tension) += (HEAVISIDE(eigen_values_matrix(0, 0)) + HEAVISIDE(eigen_values_matrix(1, 1))) *
            outer_prod(projection_vector_cross, projection_vector_cross);
#endif //PROJECTION_OPERATOR_CERVERA_2017

        noalias(projection_tensor_compression) = IdentityMatrix(3, 3) - projection_tensor_tension;
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::CalculateEquivalentStressTension(CalculationData& data, double& UniaxialStressTension)
    {
        UniaxialStressTension = 0.0;
        if (data.PrincipalStressVector(0) > 0.0) {
            if (data.TensionYieldModel == 0) {
                // Lubliner Yield Criteria 
                const double yield_compression = data.YieldStressCompression;
                const double yield_tension = data.YieldStressTension;
                const double alpha = (data.BiaxialCompressionMultiplier - 1.0) /
                    (2.0 * data.BiaxialCompressionMultiplier - 1.0);
                double I1, J2;
                array_1d<double, 3> deviator = ZeroVector(3);

                ConstitutiveLawUtilities<3>::CalculateI1Invariant(data.EffectiveStressVector, I1);
                ConstitutiveLawUtilities<3>::CalculateJ2Invariant(data.EffectiveStressVector, I1, deviator, J2);

                const double beta = yield_compression / yield_tension * (1.0 - alpha) - (1.0 + alpha);
                const double smax = std::max(std::max(data.PrincipalStressVector(0), data.PrincipalStressVector(1)), 0.0);

                UniaxialStressTension = 1.0 / (1.0 - alpha) * (alpha * I1 + std::sqrt(3.0 * J2) + beta * smax) /
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
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::CalculateEquivalentStressCompression(CalculationData& data, double& UniaxialStressCompression)
    {
        UniaxialStressCompression = 0.0;
        if (data.PrincipalStressVector(1) < 0.0) {
            const double yield_compression = data.DamageOnsetStressCompression;
            const double yield_tension = data.YieldStressTension;
            const double alpha = (data.BiaxialCompressionMultiplier - 1.0) /
                (2.0 * data.BiaxialCompressionMultiplier - 1.0);
            double I1, J2;
            array_1d<double, 3> deviator = ZeroVector(3);

            ConstitutiveLawUtilities<3>::CalculateI1Invariant(data.EffectiveStressVector, I1);
            ConstitutiveLawUtilities<3>::CalculateJ2Invariant(data.EffectiveStressVector, I1, deviator, J2);

            const double beta = (yield_compression / yield_tension) * (1.0 - alpha) - (1.0 + alpha);
            const double smax = std::max(std::max(data.PrincipalStressVector(0), data.PrincipalStressVector(1)), 0.0);

            UniaxialStressCompression = 1.0 / (1.0 - alpha) * (alpha * I1 + std::sqrt(3.0 * J2) +
                data.ShearCompressionReductor * beta * smax);
        }
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::CalculateDamageTension(
        CalculationData& data,
        double internal_variable,
        double& rDamage)
    {
        if (internal_variable <= data.YieldStressTension) {
            rDamage = 0.0;
        }
        else {
            const double characteristic_length = data.CharacteristicLength;
            const double young_modulus = data.YoungModulus;
            const double yield_tension = data.YieldStressTension;
            const double initial_internal_variable = yield_tension;
            const double material_length = 2.0 * young_modulus * data.FractureEnergyTension /
                (yield_tension * yield_tension);

            if (characteristic_length >= material_length) {
                std::stringstream ss;
                ss << "FRACTURE_ENERGY_TENSION is to low:  2*E*Gt/(ft*ft) = " << material_length
                    << ",   Characteristic Length = " << characteristic_length << std::endl;
                std::cout << ss.str();
                exit(-1);
            }

            const double damage_parameter = 2.0 * characteristic_length /
                (material_length - characteristic_length);

            rDamage = 1.0 - initial_internal_variable / internal_variable *
                std::exp(damage_parameter *
                (1.0 - internal_variable / initial_internal_variable));

            const double internal_variable_min = 1.0e-2 * yield_tension;
            if ((1.0 - rDamage) * internal_variable < internal_variable_min) {
                rDamage = 1.0 - internal_variable_min / internal_variable;
            }
        }
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::CalculateDamageCompression(
        CalculationData& data,
        double internal_variable,
        double& rDamage)
    {
        if (internal_variable <= data.DamageOnsetStressCompression) {
            rDamage = 0.0;
        }
        else {
            // extract material parameters
            const double young_modulus = data.YoungModulus;
            const double s_0 = data.DamageOnsetStressCompression;
            const double s_p = data.YieldStressCompression;
            const double s_r = data.ResidualStressCompression;
            const double e_p = data.YieldStrainCompression;
            const double c1 = data.BezierControllerC1;
            const double c2 = data.BezierControllerC2;
            const double c3 = data.BezierControllerC3;
            const double specific_fracture_energy = data.FractureEnergyCompression / 
                                                    data.CharacteristicLength;

            // Auto-computation of remaining constitutive law parameters
            const double s_k = s_r + (s_p - s_r) * c1;
            const double e_0 = s_0 / young_modulus;
            const double e_i = s_p / young_modulus;
            const double alpha = 2.0 * (e_p - s_p / young_modulus);
            double e_j = e_p + alpha * c2;
            double e_k = e_j + alpha * (1 - c2);
            double e_r = (e_k - e_j) / (s_p - s_k) * (s_p - s_r) + e_j;
            double e_u = e_r * c3;

            // regularization
            double bezier_fracture_energy, bezier_energy_1;
            this->ComputeBezierEnergy(bezier_fracture_energy, bezier_energy_1,
                                        s_p, s_k, s_r, e_p, e_j, e_k, e_r, e_u);

            const double stretcher = (specific_fracture_energy - bezier_energy_1) /
                                     (bezier_fracture_energy - bezier_energy_1) - 1.0;

            if(stretcher <= -1.0){
                std::stringstream ss;
                ss << "Damage TC Error: Compressive fracture energy is too low" << std::endl;
                ss << "Input Gc/lch = " << specific_fracture_energy << std::endl;
                ss << "Minimum Gc to avoid constitutive snap-back = " << bezier_energy_1 << std::endl;
                std::cout << ss.str();
                exit(-1);
                }

            this->ApplyBezierStretcherToStrains(stretcher, e_p, e_j, e_k, e_r, e_u);

            // current abscissa
            const double strain_like_counterpart = internal_variable / young_modulus;

            // compute damage
            double damage_variable = internal_variable;
            if (strain_like_counterpart <= e_p) {
                this->EvaluateBezierCurve(damage_variable, strain_like_counterpart, e_0, e_i, e_p, s_0, s_p, s_p);
            }
            else if (strain_like_counterpart <= e_k) {
                this->EvaluateBezierCurve(damage_variable, strain_like_counterpart, e_p, e_j, e_k, s_p, s_p, s_k);
            }
            else if (strain_like_counterpart <= e_u) {
                this->EvaluateBezierCurve(damage_variable, strain_like_counterpart, e_k, e_r, e_u, s_k, s_r, s_r);
            }
            else {
                damage_variable = s_r;
            }
            rDamage = 1.0 - damage_variable / internal_variable;
        }
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::ComputeBezierEnergy(double& rBezierEnergy, double& rBezierEnergy1,
        double s_p, double s_k, double s_r,
        double e_p, double e_j, double e_k, double e_r, double e_u)
    {
        rBezierEnergy1 = e_p * s_p / 2.0;
        double bezier_energy_2 = this->EvaluateBezierArea(e_p, e_j, e_k, s_p, s_p, s_k);
        double bezier_energy_3 = this->EvaluateBezierArea(e_k, e_r, e_u, s_k, s_r, s_r);
        rBezierEnergy = rBezierEnergy1 + bezier_energy_2 + bezier_energy_3;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    double MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::EvaluateBezierArea(
        double x1, double x2, double x3, double y1, double y2, double y3)
    {
        double bezier_area = x2 * y1 / 3.0 +
            x3 * y1 / 6.0 -
            x2 * y3 / 3.0 +
            x3 * y2 / 3.0 +
            x3 * y3 / 2.0 -
            x1 * (y1 / 2.0 + y2 / 3.0 + y3 / 6.0);
        return bezier_area;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::ApplyBezierStretcherToStrains(
        double stretcher, double e_p, double& e_j, double& e_k, double& e_r, double& e_u)
    {
        e_j += (e_j - e_p) * stretcher;
        e_k += (e_k - e_p) * stretcher;
        e_r += (e_r - e_p) * stretcher;
        e_u += (e_u - e_p) * stretcher;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::EvaluateBezierCurve(
        double& rDamageParameter, double xi,
        double x1, double x2, double x3,
        double y1, double y2, double y3)
    {
        double bezier_law_param_A = x1 - 2.0 * x2 + x3;
        double bezier_law_param_B = 2.0 * (x2 - x1);
        double bezier_law_param_C = x1 - xi;

        if (std::abs(bezier_law_param_A) < 1.0E-12) {
            x2 = x2 + 1.0E-6 * (x3 - x1);
            bezier_law_param_A = x1 - 2.0 * x2 + x3;
            bezier_law_param_B = 2.0 * (x2 - x1);
            bezier_law_param_C = x1 - xi;
        }

        double bezier_law_param_D = bezier_law_param_B * bezier_law_param_B -
            4.0 * bezier_law_param_A * bezier_law_param_C;
        double bezier_law_param_t = (-bezier_law_param_B + std::sqrt(bezier_law_param_D)) /
            (2.0 * bezier_law_param_A);
        rDamageParameter = (y1 - 2.0 * y2 + y3) * bezier_law_param_t * bezier_law_param_t +
            2.0 * (y2 - y1) * bezier_law_param_t + y1;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::ComputeCharacteristicLength(
        const GeometryType& rGeometry,
        double& rCharacteristicLength)
    {
        //rCharacteristicLength = geom.Length();
//
//#ifdef OPTIMIZE_CHARACTERISTIC_LENGTH
        //if (geom.WorkingSpaceDimension() == 2 && geom.PointsNumber() == 4) {
            //2D Element with 4 Nodes
            //double aX = (geom[0].X0() + geom[2].X0()) / 2.0;  //center_coord_X_Node1Node4
            //double aY = (geom[0].Y0() + geom[2].Y0()) / 2.0;  //center_coord_Y_Node1Node4
            //double bX = (geom[1].X0() + geom[3].X0()) / 2.0;  //center_coord_X_Node2Node3
            //double bY = (geom[1].Y0() + geom[3].Y0()) / 2.0;  //center_coord_Y_Node2Node3
            //double cX = (geom[0].X0() + geom[1].X0()) / 2.0;  //center_coord_X_Node1Node2
            //double cY = (geom[0].Y0() + geom[1].Y0()) / 2.0;  //center_coord_Y_Node1Node2
            //double dX = (geom[3].X0() + geom[2].X0()) / 2.0;  //center_coord_X_Node3Node4
            //double dY = (geom[3].Y0() + geom[2].Y0()) / 2.0;  //center_coord_Y_Node3Node4

            //double SabX = aX - bX;
            //double SabY = aY - bY;
            //double ScdX = cX - dX;
            //double ScdY = cY - dY;

            //double length_ab = std::sqrt(std::pow(SabX, 2) + std::pow(SabY, 2));
            //double length_cd = std::sqrt(std::pow(ScdX, 2) + std::pow(ScdY, 2));

            //KRATOS_WATCH(geom.size())

            //KRATOS_WATCH(length_ab)
            //    KRATOS_WATCH(length_cd)

            //rCharacteristicLength = std::min(length_ab, length_cd);

            double radius = 0.0;

            const SizeType points_number = rGeometry.size();
            array_1d<double, 3> center = rGeometry[0].Coordinates();
            for (IndexType i_node = 1; i_node < points_number; ++i_node) {
                center += rGeometry[i_node].Coordinates();
            }
            center /= static_cast<double>(points_number);

            for (IndexType i_node = 0; i_node < points_number; ++i_node) {
                const array_1d<double, 3>& aux_vector = center - rGeometry[i_node].Coordinates();
                double aux_value = inner_prod(aux_vector, aux_vector);
                if (aux_value > radius) radius = aux_value;
            }

            rCharacteristicLength = std::sqrt(radius);
        //}
//#endif
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::CalculateMaterialResponseInternal(
        const Vector StrainVectorIso,
        Vector& PredictiveStressVector,
        CalculationData& data,
        int IntegrationImplex)
    {
        if (PredictiveStressVector.size() != VoigtSize)
            PredictiveStressVector.resize(VoigtSize, false);

        ///////////////////////////////////////////////////////////////////////////
        /// WHY
        mThresholdTension = mCurrentThresholdTension;
        mThresholdCompression = mCurrentThresholdCompression;

        //Vector StrainVectorIso = prod(data.TransformationMatrix, StrainVector);

        noalias(data.EffectiveStressVector) = prod(data.ElasticityMatrix, StrainVectorIso);

        if (std::abs(data.EffectiveStressVector(0)) < tolerance) { data.EffectiveStressVector(0) = 0.0; }
        if (std::abs(data.EffectiveStressVector(1)) < tolerance) { data.EffectiveStressVector(1) = 0.0; }
        if (std::abs(data.EffectiveStressVector(2)) < tolerance) { data.EffectiveStressVector(2) = 0.0; }

        this->TensionCompressionSplit(data);
        this->ConstructProjectionTensors(data);

        // compute the equivalent stress measures
        this->CalculateEquivalentStressTension(data, UniaxialStressTension);

        this->CalculateEquivalentStressCompression(data, UniaxialStressCompression);

        // damage update
        if (IntegrationImplex != 0) { //IMPLEX Integration
            // time factor
            double time_factor = 0.0;
            if (PreviousDeltaTime > 0.0) time_factor = data.DeltaTime / PreviousDeltaTime;
            CurrentDeltaTime = data.DeltaTime;

            // explicit evaluation
            mThresholdTension = mCurrentThresholdTension + time_factor * (mCurrentThresholdTension - PreviousThresholdTension);
            mThresholdCompression = mCurrentThresholdCompression + time_factor * (mCurrentThresholdCompression - PreviousThresholdCompression);

            // save implicit variables for the finalize_solution_step
            double implicit_threshold_tension = mCurrentThresholdTension;
            double implicit_threshold_compression = mCurrentThresholdCompression;

            if (UniaxialStressTension > implicit_threshold_tension)
                implicit_threshold_tension = UniaxialStressTension;

            if (UniaxialStressCompression > implicit_threshold_compression)
                implicit_threshold_compression = UniaxialStressCompression;

            TemporaryImplicitThresholdTension = implicit_threshold_tension;
            TemporaryImplicitThresholdTCompression = implicit_threshold_compression;

            // new damage variables (explicit)
            this->CalculateDamageTension(data, mThresholdTension, mDamageTension);
            this->CalculateDamageCompression(data, mThresholdCompression, mDamageCompression);
        }
        else { // IMPLICIT Integration
            if (UniaxialStressTension > mThresholdTension)
                mThresholdTension = UniaxialStressTension;
            this->CalculateDamageTension(data, mThresholdTension, mDamageTension);

            if (UniaxialStressCompression > mThresholdCompression)
                mThresholdCompression = UniaxialStressCompression;
            this->CalculateDamageCompression(data, mThresholdCompression, mDamageCompression);

            mCurrentThresholdTension = mThresholdTension;
            mCurrentThresholdCompression = mThresholdCompression;
        }

        // calculation of isotropic stress tensor
        noalias(PredictiveStressVector) = (1.0 - mDamageTension) * data.EffectiveTensionStressVector;
        noalias(PredictiveStressVector) += (1.0 - mDamageCompression) * data.EffectiveCompressionStressVector;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::CheckDamageLoadingUnloading(
        bool& is_damaging_tension,
        bool& is_damaging_compression)
    {
        const double F_tension = UniaxialStressTension - mCurrentThresholdTension;
        const double F_compression = UniaxialStressCompression - mCurrentThresholdCompression;

        is_damaging_tension = (F_tension > 0.0)
            ? true : false;
        is_damaging_compression = (F_compression > 0.0)
            ? true : false;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::CalculateTangentTensor(
        Parameters& rValues,
        Vector StrainVector,
        Vector PredictiveStressVector,
        CalculationData& data,
        Properties props)
    {
        // prepare constitutive matrix
        Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
        if (constitutive_matrix.size1() != VoigtSize || constitutive_matrix.size2() != VoigtSize)
            constitutive_matrix.resize(VoigtSize, VoigtSize);

        // save internal variables
        double save_threshold_tension = mThresholdTension;
        double save_threshold_compression = mThresholdCompression;
        double save_damage_tension = mDamageTension;
        double save_damage_compression = mDamageCompression;
        double save_uniaxial_stress_tension = UniaxialStressTension;
        double save_uniaxial_stress_compression = UniaxialStressCompression;

        // perturbation parameter
        double perturbation_factor = 1.0E-8;

        // perturbed vectors
        Vector perturbated_strain_vector(VoigtSize);
        Vector perturbated_stress_vector(VoigtSize);

        // apply perturbation to each strain component...
        for (size_t j = 0; j < VoigtSize; j++)
        {
            noalias(perturbated_strain_vector) = StrainVector;

            perturbated_strain_vector(j) = StrainVector(j) + perturbation_factor;
            int integration_implex = props[INTEGRATION_IMPLEX];
            this->CalculateMaterialResponseInternal(
                perturbated_strain_vector,
                perturbated_stress_vector,
                data,
                integration_implex);

            for (size_t i = 0; i < VoigtSize; i++)
                constitutive_matrix(i, j) = (perturbated_stress_vector(i) - PredictiveStressVector(i)) /
                perturbation_factor;
        }

        // restore internal variables
        mThresholdTension = save_threshold_tension;
        mThresholdCompression = save_threshold_compression;
        mDamageTension = save_damage_tension;
        mDamageCompression = save_damage_compression;
        UniaxialStressTension = save_uniaxial_stress_tension;
        UniaxialStressCompression = save_uniaxial_stress_compression;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2D::CalculateSecantTensor(
        Parameters& rParameters,
        CalculationData& data)
    {
        Matrix& constitutive_matrix = rParameters.GetConstitutiveMatrix();
        if (constitutive_matrix.size1() != VoigtSize || constitutive_matrix.size2() != VoigtSize)
            constitutive_matrix.resize(VoigtSize, VoigtSize);

        Matrix DamageMatrix(IdentityMatrix(3, 3));
        noalias(DamageMatrix) -= mDamageTension * data.ProjectionTensorTension;
        noalias(DamageMatrix) -= mDamageCompression * data.ProjectionTensorCompression;

        noalias(constitutive_matrix) = prod(DamageMatrix, data.ElasticityMatrix);
    }

} // namespace Kratos
