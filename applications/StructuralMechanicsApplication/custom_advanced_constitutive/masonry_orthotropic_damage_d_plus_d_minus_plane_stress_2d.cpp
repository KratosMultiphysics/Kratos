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
    bool MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::Has(
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
    double& MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::GetValue(
        const Variable<double>& rThisVariable,
        double& rValue)
    {
        rValue = 0.0;
        // TENSION
        if (rThisVariable == DAMAGE_TENSION)
            rValue = mDamageTension;
        else if (rThisVariable == UNIAXIAL_STRESS_TENSION)
            rValue = mUniaxialStressTension;
        else if (rThisVariable == THRESHOLD_TENSION)
            rValue = mThresholdTension;
        // COMPRESSION
        else if (rThisVariable == DAMAGE_COMPRESSION)
            rValue = mDamageCompression;
        else if (rThisVariable == UNIAXIAL_STRESS_COMPRESSION)
            rValue = mUniaxialStressCompression;
        else if (rThisVariable == THRESHOLD_COMPRESSION)
            rValue = mThresholdCompression;
        return rValue;
    }

    /// SetValue variable <double>
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::SetValue(
        const Variable<double>& rVariable,
        const double& rValue,
        const ProcessInfo& rCurrentProcessInfo)
    {
        // TENSION
        if (rVariable == DAMAGE_TENSION)
            mDamageTension = rValue;
        else if (rVariable == UNIAXIAL_STRESS_TENSION)
            mUniaxialStressTension = rValue;
        else if (rVariable == THRESHOLD_TENSION)
            mThresholdTension = rValue;
        // COMPRESSION
        else if (rVariable == DAMAGE_COMPRESSION)
            mDamageCompression = rValue;
        else if (rVariable == UNIAXIAL_STRESS_COMPRESSION)
            mUniaxialStressCompression = rValue;
        else if (rVariable == THRESHOLD_COMPRESSION)
            mThresholdCompression = rValue;
    }

    /***********************************************************************************/
    /***********************************************************************************/
    bool MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::ValidateInput(
        const Properties& rProperties)
    {
        KRATOS_ERROR << "just because I can";
        // ELASTIC PARAMETERS
        if (!rProperties.Has(YOUNG_MODULUS_1))                        return false;
        if (!rProperties.Has(YOUNG_MODULUS_2))                        return false;
        if (!rProperties.Has(POISSON_RATIO_12))                       return false;
        if (!rProperties.Has(POISSON_RATIO_21))                       return false;
        // TENSION PARAMETERS
        if (!rProperties.Has(YIELD_STRESS_TENSION))                   return false;
        if( !rProperties.Has(FRACTURE_ENERGY_TENSION) )               return false;
        if (!rProperties.Has(YIELD_STRESS_TENSION_2))                 return false;
        if (!rProperties.Has(FRACTURE_ENERGY_TENSION_2))              return false;
        // COMPRESSION PARAMETERS
        if (!rProperties.Has(ELASTIC_LIMIT_STRESS_COMPRESSION))       return false;
        if (!rProperties.Has(YIELD_STRESS_COMPRESSION))               return false;
        if (!rProperties.Has(YIELD_STRAIN_COMPRESSION))               return false;
        if (!rProperties.Has(RESIDUAL_STRESS_COMPRESSION))            return false;
        if( !rProperties.Has(FRACTURE_ENERGY_COMPRESSION) )           return false;
        if (!rProperties.Has(ELASTIC_LIMIT_STRESS_COMPRESSION_2))     return false;
        if (!rProperties.Has(YIELD_STRESS_COMPRESSION_2))             return false;
        if (!rProperties.Has(YIELD_STRAIN_COMPRESSION_2))             return false;
        if (!rProperties.Has(RESIDUAL_STRESS_COMPRESSION_2))          return false;
        if (!rProperties.Has(FRACTURE_ENERGY_COMPRESSION_2))          return false;
        // BIAXIAL PARAMETERS
        if (!rProperties.Has(BIAXIAL_COMPRESSION_MULTIPLIER))         return false;
        return true;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::InitializeMaterial(
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
            mUniaxialStressTension = 0.0;
            mUniaxialStressCompression = 0.0;

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
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::FinalizeSolutionStep(
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
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::CalculateMaterialResponsePK1(
        Parameters& rValues)
    {
        CalculateMaterialResponseCauchy(rValues);
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::CalculateMaterialResponsePK2(
        Parameters& rValues)
    {
        CalculateMaterialResponseCauchy(rValues);
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::CalculateMaterialResponseKirchhoff(
        Parameters& rValues)
    {
        CalculateMaterialResponseCauchy(rValues);
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::CalculateMaterialResponseCauchy(
        Parameters& rParameters)
    {
        const ProcessInfo& r_process_info = rParameters.GetProcessInfo();
        const GeometryType& r_geometry = rParameters.GetElementGeometry();
        const Properties& r_properties = rParameters.GetMaterialProperties();

        const Vector& r_strain_vector = rParameters.GetStrainVector();
        Vector& r_stress_vector = rParameters.GetStressVector();

        CalculationData data = GetCalculationData(r_properties, r_geometry, r_process_info);

        Matrix ElasticityMatrix = ZeroMatrix(3, 3);
        CalculateOrthotropicElasticityMatrix(
            ElasticityMatrix, data.MaterialProperties1, data.MaterialProperties2);

        noalias(data.EffectiveStressVector) = prod(ElasticityMatrix, r_strain_vector);

        this->TensionCompressionSplit(
            data.EffectiveStressVector,
            data.PrincipalStressVector,
            data.EffectiveTensionStressVector,
            data.EffectiveCompressionStressVector);
        this->ConstructProjectionTensors(
            data.EffectiveStressVector,
            data.ProjectionTensorTension,
            data.ProjectionTensorCompression);

        // Compute transformation from real anisotropic stress space to mapped isotropic stress space
        TransformationMatrices transformation_matrices;
        AssembleTransformationMatrix(
            data.MaterialProperties1, data.MaterialProperties2, transformation_matrices);

        // Transform real effective stresses to the mapped space
        array_1d<double, 3> effective_tension_stress_mapped_isotropic = prod(transformation_matrices.TransformationMatrixTension, data.EffectiveTensionStressVector);
        array_1d<double, 3> effective_compression_stress_mapped_isotropic = prod(transformation_matrices.TransformationMatrixCompression, data.EffectiveCompressionStressVector);

        // compute the equivalent stress measures
        this->CalculateEquivalentStressTension(data, data.MaterialProperties1, effective_tension_stress_mapped_isotropic, mUniaxialStressTension);

        this->CalculateEquivalentStressCompression(data, data.MaterialProperties1, effective_compression_stress_mapped_isotropic, mUniaxialStressCompression);

        // damage update
        if (r_properties[INTEGRATION_IMPLEX] != 0) { //IMPLEX Integration
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

            if (mUniaxialStressTension > implicit_threshold_tension)
                implicit_threshold_tension = mUniaxialStressTension;

            if (mUniaxialStressCompression > implicit_threshold_compression)
                implicit_threshold_compression = mUniaxialStressCompression;

            TemporaryImplicitThresholdTension = implicit_threshold_tension;
            TemporaryImplicitThresholdTCompression = implicit_threshold_compression;

            // new damage variables (explicit)
            this->CalculateDamageTension(data.MaterialProperties1, mThresholdTension, mDamageTension);
            this->CalculateDamageCompression(data.MaterialProperties1, mThresholdCompression, mDamageCompression);
        } // IMPLICIT Integration
        else {
            if (mUniaxialStressTension > mThresholdTension)
                mThresholdTension = mUniaxialStressTension;
            this->CalculateDamageTension(data.MaterialProperties1, mThresholdTension, mDamageTension);

            if (mUniaxialStressCompression > mThresholdCompression)
                mThresholdCompression = mUniaxialStressCompression;
            this->CalculateDamageCompression(data.MaterialProperties1, mThresholdCompression, mDamageCompression);

            mCurrentThresholdTension = mThresholdTension;
            mCurrentThresholdCompression = mThresholdCompression;
        }

        // Compute the stresses to isotropic space
        array_1d<double, 3> total_tension_stress_mapped_isotropic = (1.0 - mDamageTension) * effective_tension_stress_mapped_isotropic;
        array_1d<double, 3> total_compression_stress_mapped_isotropic = (1.0 - mDamageCompression) * effective_compression_stress_mapped_isotropic;

        // Back rotation of final Predictive Stress
        noalias(r_stress_vector) = prod(transformation_matrices.InverseTransformationMatrixTension, total_tension_stress_mapped_isotropic);
        noalias(r_stress_vector) += prod(transformation_matrices.InverseTransformationMatrixCompression, total_compression_stress_mapped_isotropic);

        bool is_damaging_tension = false;
        bool is_damaging_compression = false;
        this->CheckDamageLoadingUnloading(is_damaging_tension, is_damaging_compression);

        // Computation of the Constitutive Tensor
        if (rParameters.GetOptions().Is(COMPUTE_CONSTITUTIVE_TENSOR)) {
            if (is_damaging_tension || is_damaging_compression) {
                this->CalculateSecantTensor(rParameters, ElasticityMatrix, data);
                //this->CalculateTangentTensor(rParameters, strain_mapped_isotropic, r_predictive_stress_vector, data, r_properties[INTEGRATION_IMPLEX]);
            }
            else {
                this->CalculateSecantTensor(rParameters, ElasticityMatrix, data);
            }
        }
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::ResetMaterial(
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
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::GetLawFeatures(
        Features& rFeatures)
    {
        //Set the type of law
        rFeatures.mOptions.Set(PLANE_STRESS_LAW);
        rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
        rFeatures.mOptions.Set(ANISOTROPIC);

        //Set strain measure required by the consitutive law
        rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);

        //Set the strain size
        rFeatures.mStrainSize = GetStrainSize();

        //Set the space dimension
        rFeatures.mSpaceDimension = WorkingSpaceDimension();
    }
    /***********************************************************************************/
    /***********************************************************************************/
    int MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::Check(
        const Properties& rProperties,
        const GeometryType& rGeometry,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        // ELASTIC PARAMETERS
        if (!rProperties.Has(YOUNG_MODULUS_1))
            KRATOS_THROW_ERROR(std::logic_error, "Missing variable: YOUNG_MODULUS_1", "");
        if (!rProperties.Has(YOUNG_MODULUS_2))
            KRATOS_THROW_ERROR(std::logic_error, "Missing variable: YOUNG_MODULUS_2", "");
        if (!rProperties.Has(POISSON_RATIO_12))
            KRATOS_THROW_ERROR(std::logic_error, "Missing variable: POISSON_RATIO_12", "");
        if (!rProperties.Has(POISSON_RATIO_12))
            KRATOS_THROW_ERROR(std::logic_error, "Missing variable: POISSON_RATIO_21", "");
        if (!rProperties.Has(SHEAR_MODULUS))
            KRATOS_THROW_ERROR(std::logic_error, "Missing variable: SHEAR_MODULUS", "");
        // TENSION
        if (!rProperties.Has(YIELD_STRESS_TENSION))
            KRATOS_THROW_ERROR(std::logic_error, "Missing variable: YIELD_STRESS_TENSION", "");
        if( !rProperties.Has(FRACTURE_ENERGY_TENSION) )
            KRATOS_THROW_ERROR(std::logic_error, "Missing variable: FRACTURE_ENERGY_TENSION", "");
        // COMPRESSION
        if (!rProperties.Has(DAMAGE_ONSET_STRESS_COMPRESSION))
            KRATOS_THROW_ERROR(std::logic_error, "Missing variable: DAMAGE_ONSET_STRESS_COMPRESSION", "");
        if (!rProperties.Has(YIELD_STRESS_COMPRESSION))
            KRATOS_THROW_ERROR(std::logic_error, "Missing variable: YIELD_STRESS_COMPRESSION", "");
        if (!rProperties.Has(YIELD_STRAIN_COMPRESSION))
            KRATOS_THROW_ERROR(std::logic_error, "Missing variable: YIELD_STRAIN_COMPRESSION", "");
        if (!rProperties.Has(RESIDUAL_STRESS_COMPRESSION))
            KRATOS_THROW_ERROR(std::logic_error, "Missing variable: RESIDUAL_STRESS_COMPRESSION", "");
        if( !rProperties.Has(FRACTURE_ENERGY_COMPRESSION) )
            KRATOS_THROW_ERROR(std::logic_error, "Missing variable: FRACTURE_ENERGY_COMPRESSION", "");
        // BIAXIAL
        if (!rProperties.Has(BIAXIAL_COMPRESSION_MULTIPLIER))
            KRATOS_THROW_ERROR(std::logic_error, "Missing variable: BIAXIAL_COMPRESSION_MULTIPLIER", "");

        return 0;

        KRATOS_CATCH("");
    }
    /***********************************************************************************/
    /***********************************************************************************/
    MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::CalculationData MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::GetCalculationData(
        const Properties& rProperties,
        const GeometryType& rGeometry,
        const ProcessInfo& rProcessInfo)
    {
        auto MaterialProperties1 = DirectionalMaterialProperties(
            this->ComputeCharacteristicLength(
                rGeometry, 0),
            rProperties[YOUNG_MODULUS_1],
            rProperties[POISSON_RATIO_12],
            rProperties[SHEAR_MODULUS],
            rProperties[YIELD_STRESS_TENSION],
            rProperties[FRACTURE_ENERGY_TENSION],
            rProperties[ELASTIC_LIMIT_STRESS_COMPRESSION],
            rProperties[YIELD_STRESS_COMPRESSION],
            rProperties[YIELD_STRAIN_COMPRESSION],
            rProperties[RESIDUAL_STRESS_COMPRESSION],
            rProperties[FRACTURE_ENERGY_COMPRESSION],
            rProperties[BEZIER_CONTROLLER_C1],
            rProperties[BEZIER_CONTROLLER_C2],
            rProperties[BEZIER_CONTROLLER_C3]);

        auto MaterialProperties2 = DirectionalMaterialProperties(
            this->ComputeCharacteristicLength(
                rGeometry, 1),
            rProperties[YOUNG_MODULUS_2],
            rProperties[POISSON_RATIO_21],
            rProperties[SHEAR_MODULUS],
            rProperties[YIELD_STRESS_TENSION_2],
            rProperties[FRACTURE_ENERGY_TENSION_2],
            rProperties[ELASTIC_LIMIT_STRESS_COMPRESSION_2],
            rProperties[YIELD_STRESS_COMPRESSION_2],
            rProperties[YIELD_STRAIN_COMPRESSION_2],
            rProperties[RESIDUAL_STRESS_COMPRESSION_2],
            rProperties[FRACTURE_ENERGY_COMPRESSION_2],
            rProperties[BEZIER_CONTROLLER_C1_2],
            rProperties[BEZIER_CONTROLLER_C2_2],
            rProperties[BEZIER_CONTROLLER_C3_2]);

        CalculationData data(MaterialProperties1, MaterialProperties2);

        // Biaxial
        data.BiaxialCompressionMultiplier = rProperties[BIAXIAL_COMPRESSION_MULTIPLIER];
        data.ShearCompressionReductor = rProperties.Has(SHEAR_COMPRESSION_REDUCTOR) ? rProperties[SHEAR_COMPRESSION_REDUCTOR] : 0.5;
        data.ShearCompressionReductor = std::min(std::max(data.ShearCompressionReductor, 0.0), 1.0);

        // Misc
        data.DeltaTime = rProcessInfo[DELTA_TIME];
        data.TensionYieldModel = rProperties.Has(TENSION_YIELD_MODEL) ? rProperties[TENSION_YIELD_MODEL] : 0;

        return data;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::CalculateOrthotropicElasticityMatrix(
        Matrix& rC,
        const DirectionalMaterialProperties& rMaterialProperties1,
        const DirectionalMaterialProperties& rMaterialProperties2)
    {
        if (rC.size1() != 3 || rC.size2() != 3)
            rC.resize(3, 3, false);
        rC = ZeroMatrix(3, 3);

        const double nu_12_nu_21 = 1.0 - rMaterialProperties1.nu * rMaterialProperties2.nu;

        rC(0, 0) = rMaterialProperties1.E / nu_12_nu_21; rC(0, 1) = rMaterialProperties1.nu * rC(0, 0);
        rC(1, 1) = rMaterialProperties2.E / nu_12_nu_21; rC(1, 0) = rMaterialProperties2.nu * rC(1, 1);
        rC(2, 2) = rMaterialProperties1.G;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::AssembleTransformationMatrix(
        const DirectionalMaterialProperties& rMaterialProperties1,
        const DirectionalMaterialProperties& rMaterialProperties2,
        TransformationMatrices& rTransformationMatrices)
    {
        double shear_resistance = rMaterialProperties1.YieldStressTension;
        rTransformationMatrices.TransformationMatrixTension(0, 0) = 1;
        rTransformationMatrices.TransformationMatrixTension(1, 1) = rMaterialProperties1.YieldStressTension / rMaterialProperties2.YieldStressTension;
        rTransformationMatrices.TransformationMatrixTension(2, 2) = rMaterialProperties2.YieldStressTension / shear_resistance;

        rTransformationMatrices.InverseTransformationMatrixTension(0, 0) = 1;
        rTransformationMatrices.InverseTransformationMatrixTension(1, 1) = rMaterialProperties2.YieldStressTension / rMaterialProperties1.YieldStressTension;
        rTransformationMatrices.InverseTransformationMatrixTension(2, 2) = shear_resistance / rMaterialProperties1.YieldStressTension;

        double shear_resistance_compression = rMaterialProperties1.YieldStressCompression;
        rTransformationMatrices.TransformationMatrixCompression(0, 0) = 1;
        rTransformationMatrices.TransformationMatrixCompression(1, 1) = rMaterialProperties1.YieldStressCompression / rMaterialProperties2.YieldStressCompression;
        rTransformationMatrices.TransformationMatrixCompression(2, 2) = rMaterialProperties2.YieldStressCompression / shear_resistance_compression;

        rTransformationMatrices.InverseTransformationMatrixCompression(0, 0) = 1;
        rTransformationMatrices.InverseTransformationMatrixCompression(1, 1) = rMaterialProperties2.YieldStressCompression / rMaterialProperties1.YieldStressCompression;
        rTransformationMatrices.InverseTransformationMatrixCompression(2, 2) = shear_resistance_compression / rMaterialProperties1.YieldStressCompression;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::AssembleTransformationMatrixEnergyEquivalent(
        const DirectionalMaterialProperties& rMaterialProperties1,
        const DirectionalMaterialProperties& rMaterialProperties2,
        TransformationMatrices& rTransformationMatrices)
    {

    }
    /***********************************************************************************/
    /***********************************************************************************/

    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::TensionCompressionSplit(
        const array_1d<double, 3>& rEffectiveStressVector,
        array_1d<double, 2>& rPrincipalStressVector,
        array_1d<double, 3>& rEffectiveTensionStressVector,
        array_1d<double, 3>& rEffectiveCompressionStressVector)
    {
        ConstitutiveLawUtilities<3>::CalculatePrincipalStresses(
            rPrincipalStressVector, rEffectiveStressVector);
        ConstitutiveLawUtilities<3>::SpectralDecomposition(
            rEffectiveStressVector, rEffectiveTensionStressVector, rEffectiveCompressionStressVector);
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::ConstructProjectionTensors(
        const array_1d<double, 3>& rEffectiveStressVector,
        Matrix& rProjectionTensorTension,
        Matrix& rProjectionTensorCompression)
    {
        Matrix effective_stress_tensor = MathUtils<double>::StressVectorToTensor(rEffectiveStressVector);
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

        rProjectionTensorTension = ZeroMatrix(3, 3);
        noalias(rProjectionTensorTension) += HEAVISIDE(eigen_values_matrix(0, 0)) *
            outer_prod(projection_vector_11, projection_vector_11);
        noalias(rProjectionTensorTension) += HEAVISIDE(eigen_values_matrix(1, 1)) *
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
        noalias(rProjectionTensorTension) += factor_12 * outer_prod(projection_vector_cross, projection_vector_cross);
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

        noalias(rProjectionTensorCompression) = IdentityMatrix(3, 3) - rProjectionTensorTension;
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::CalculateEquivalentStressTension(
        const CalculationData& data,
        const DirectionalMaterialProperties& rMaterialProperties,
        const array_1d<double, 3> rEffectiveStressVector,
        double& rUniaxialStressTension) const
    {
        rUniaxialStressTension = 0.0;
        if (data.PrincipalStressVector(0) > 0.0) {
            if (data.TensionYieldModel == 0) {
                // Lubliner Yield Criteria 
                const double yield_compression = rMaterialProperties.YieldStressCompression;
                const double yield_tension = rMaterialProperties.YieldStressTension;
                const double alpha = (data.BiaxialCompressionMultiplier - 1.0) /
                    (2.0 * data.BiaxialCompressionMultiplier - 1.0);
                double I1, J2;
                array_1d<double, 3> deviator = ZeroVector(3);

                ConstitutiveLawUtilities<3>::CalculateI1Invariant(rEffectiveStressVector, I1);
                ConstitutiveLawUtilities<3>::CalculateJ2Invariant(rEffectiveStressVector, I1, deviator, J2);

                const double beta = yield_compression / yield_tension * (1.0 - alpha) - (1.0 + alpha);
                const double smax = std::max(std::max(data.PrincipalStressVector(0), data.PrincipalStressVector(1)), 0.0);

                rUniaxialStressTension = 1.0 / (1.0 - alpha) * (alpha * I1 + std::sqrt(3.0 * J2) + beta * smax) /
                    yield_compression * yield_tension;
            }
            else if (data.TensionYieldModel == 1) {
                // Rankine Yield Criteria
                rUniaxialStressTension = std::max(std::max(data.PrincipalStressVector(0), data.PrincipalStressVector(1)), 0.0);
            }
        }
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::CalculateEquivalentStressCompression(
        const CalculationData& data,
        const DirectionalMaterialProperties& rMaterialProperties,
        const array_1d<double, 3> rEffectiveStressVector,
        double& UniaxialStressCompression) const
    {
        UniaxialStressCompression = 0.0;
        if (data.PrincipalStressVector(1) < 0.0) {
            const double yield_compression = rMaterialProperties.ElasticLimitStressCompression;
            const double yield_tension = rMaterialProperties.YieldStressTension;
            const double alpha = (data.BiaxialCompressionMultiplier - 1.0) /
                (2.0 * data.BiaxialCompressionMultiplier - 1.0);
            double I1, J2;
            array_1d<double, 3> deviator = ZeroVector(3);

            ConstitutiveLawUtilities<3>::CalculateI1Invariant(rEffectiveStressVector, I1);
            ConstitutiveLawUtilities<3>::CalculateJ2Invariant(rEffectiveStressVector, I1, deviator, J2);

            const double beta = (yield_compression / yield_tension) * (1.0 - alpha) - (1.0 + alpha);
            const double smax = std::max(std::max(data.PrincipalStressVector(0), data.PrincipalStressVector(1)), 0.0);

            UniaxialStressCompression = 1.0 / (1.0 - alpha) * (alpha * I1 + std::sqrt(3.0 * J2) +
                data.ShearCompressionReductor * beta * smax);
        }
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::CalculateDamageTension(
        const DirectionalMaterialProperties& rMaterialProperties,
        double TresholdTension,
        double& rDamageTension) const
    {
        if (TresholdTension <= rMaterialProperties.YieldStressTension) {
            rDamageTension = 0.0;
        }
        else {
            const double yield_tension = rMaterialProperties.YieldStressTension;
            const double initial_treshold_tension = yield_tension;
            const double material_length = 2.0 * rMaterialProperties.E * rMaterialProperties.FractureEnergyTension /
                (yield_tension * yield_tension);

            KRATOS_WATCH(rMaterialProperties.CharacteristicLength)

            KRATOS_ERROR_IF(rMaterialProperties.CharacteristicLength >= material_length)
                << "FRACTURE_ENERGY_TENSION is to low:  2*E*Gt/(ft*ft) = " << material_length
                << ",   Characteristic Length = " << rMaterialProperties.CharacteristicLength << std::endl;

            const double damage_parameter = 2.0 * rMaterialProperties.CharacteristicLength /
                (material_length - rMaterialProperties.CharacteristicLength);

            rDamageTension = 1.0 - initial_treshold_tension / TresholdTension *
                std::exp(damage_parameter *
                (1.0 - TresholdTension / initial_treshold_tension));

            const double min_treshold_tension = 1.0e-2 * yield_tension;
            if ((1.0 - rDamageTension) * TresholdTension < min_treshold_tension) {
                rDamageTension = 1.0 - min_treshold_tension / TresholdTension;
            }
        }
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::CalculateDamageCompression(
        const DirectionalMaterialProperties& rMaterialProperties,
        const double TresholdCompression,
        double& rDamageCompression) const
    {
        if (TresholdCompression <= rMaterialProperties.ElasticLimitStressCompression) {
            rDamageCompression = 0.0;
        }
        else {
            // extract material parameters
            const double young_modulus = rMaterialProperties.E;
            const double s_0 = rMaterialProperties.ElasticLimitStressCompression;
            const double s_p = rMaterialProperties.YieldStressCompression;
            const double s_r = rMaterialProperties.ResidualStressCompression;
            const double e_p = rMaterialProperties.YieldStrainCompression;
            const double c1 = rMaterialProperties.BezierControllerC1;
            const double c2 = rMaterialProperties.BezierControllerC2;
            const double c3 = rMaterialProperties.BezierControllerC3;
            const double specific_fracture_energy = rMaterialProperties.FractureEnergyCompression /
                rMaterialProperties.CharacteristicLength;

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

            KRATOS_ERROR_IF(stretcher <= -1.0)
                << "Damage TC Error: Compressive fracture energy is too low" << std::endl
                << "Input Gc/lch = " << specific_fracture_energy << std::endl
                << "Minimum Gc to avoid constitutive snap-back = " << bezier_energy_1 << std::endl;

            this->ApplyBezierStretcherToStrains(stretcher, e_p, e_j, e_k, e_r, e_u);

            // current abscissa
            const double strain_like_counterpart = TresholdCompression / young_modulus;

            // compute damage
            double damage_variable = TresholdCompression;
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
            rDamageCompression = 1.0 - damage_variable / TresholdCompression;
        }
    }
    /***********************************************************************************/
    /***********************************************************************************/
    inline void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::ComputeBezierEnergy(
        double& rBezierEnergy, double& rBezierEnergy1,
        double s_p, double s_k, double s_r,
        double e_p, double e_j, double e_k, double e_r, double e_u) const
    {
        rBezierEnergy1 = e_p * s_p / 2.0;
        double bezier_energy_2 = this->EvaluateBezierArea(e_p, e_j, e_k, s_p, s_p, s_k);
        double bezier_energy_3 = this->EvaluateBezierArea(e_k, e_r, e_u, s_k, s_r, s_r);
        rBezierEnergy = rBezierEnergy1 + bezier_energy_2 + bezier_energy_3;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    inline double MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::EvaluateBezierArea(
        double x1, double x2, double x3, double y1, double y2, double y3) const
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
    inline void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::ApplyBezierStretcherToStrains(
        double stretcher, double e_p, double& e_j, double& e_k, double& e_r, double& e_u) const
    {
        e_j += (e_j - e_p) * stretcher;
        e_k += (e_k - e_p) * stretcher;
        e_r += (e_r - e_p) * stretcher;
        e_u += (e_u - e_p) * stretcher;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    inline void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::EvaluateBezierCurve(
        double& rDamageParameter, double xi,
        double x1, double x2, double x3,
        double y1, double y2, double y3) const
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
    double MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::ComputeCharacteristicLength(
        const GeometryType& rGeometry,
        int DirectionIndex)
    {
        array_1d<double, 3> characteristic_lengthness;
        rGeometry.Calculate(CHARACTERISTIC_GEOMETRY_LENGTH, characteristic_lengthness);
        SizeType polynomial_degree = rGeometry.PolynomialDegree(0);
        return characteristic_lengthness[DirectionIndex] / std::sqrt(polynomial_degree);
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::CalculateMaterialResponseInternal(
        const Vector& rStrainVectorIsotropic,
        Vector& rPredictiveStressVector,
        CalculationData& data,
        int IntegrationImplex)
    {
        //if (rPredictiveStressVector.size() != VoigtSize)
        //    rPredictiveStressVector.resize(VoigtSize, false);

        //mThresholdTension = mCurrentThresholdTension;
        //mThresholdCompression = mCurrentThresholdCompression;

        //noalias(data.EffectiveStressVector) = prod(data.ElasticityMatrix, rStrainVectorIsotropic);

        //if (std::abs(data.EffectiveStressVector(0)) < tolerance) { data.EffectiveStressVector(0) = 0.0; }
        //if (std::abs(data.EffectiveStressVector(1)) < tolerance) { data.EffectiveStressVector(1) = 0.0; }
        //if (std::abs(data.EffectiveStressVector(2)) < tolerance) { data.EffectiveStressVector(2) = 0.0; }

        //this->TensionCompressionSplit(data);
        //this->ConstructProjectionTensors(data);

        //// compute the equivalent stress measures
        //this->CalculateEquivalentStressTension(data, data.MaterialProperties1, mUniaxialStressTension);

        //this->CalculateEquivalentStressCompression(data, data.MaterialProperties1, mUniaxialStressCompression);

        //// damage update
        //if (IntegrationImplex != 0) { //IMPLEX Integration
        //    // time factor
        //    double time_factor = 0.0;
        //    if (PreviousDeltaTime > 0.0) time_factor = data.DeltaTime / PreviousDeltaTime;
        //    CurrentDeltaTime = data.DeltaTime;

        //    // explicit evaluation
        //    mThresholdTension = mCurrentThresholdTension + time_factor * (mCurrentThresholdTension - PreviousThresholdTension);
        //    mThresholdCompression = mCurrentThresholdCompression + time_factor * (mCurrentThresholdCompression - PreviousThresholdCompression);

        //    // save implicit variables for the finalize_solution_step
        //    double implicit_threshold_tension = mCurrentThresholdTension;
        //    double implicit_threshold_compression = mCurrentThresholdCompression;

        //    if (mUniaxialStressTension > implicit_threshold_tension)
        //        implicit_threshold_tension = mUniaxialStressTension;

        //    if (mUniaxialStressCompression > implicit_threshold_compression)
        //        implicit_threshold_compression = mUniaxialStressCompression;

        //    TemporaryImplicitThresholdTension = implicit_threshold_tension;
        //    TemporaryImplicitThresholdTCompression = implicit_threshold_compression;

        //    // new damage variables (explicit)
        //    this->CalculateDamageTension(data.MaterialProperties1, mThresholdTension, mDamageTension);
        //    this->CalculateDamageCompression(data.MaterialProperties1, mThresholdCompression, mDamageCompression);
        //} // IMPLICIT Integration
        //else {
        //    if (mUniaxialStressTension > mThresholdTension)
        //        mThresholdTension = mUniaxialStressTension;
        //    this->CalculateDamageTension(data.MaterialProperties1, mThresholdTension, mDamageTension);

        //    if (mUniaxialStressCompression > mThresholdCompression)
        //        mThresholdCompression = mUniaxialStressCompression;
        //    this->CalculateDamageCompression(data.MaterialProperties1, mThresholdCompression, mDamageCompression);

        //    mCurrentThresholdTension = mThresholdTension;
        //    mCurrentThresholdCompression = mThresholdCompression;
        //}

        //// calculation of isotropic stress tensor
        //noalias(rPredictiveStressVector) = (1.0 - mDamageTension) * data.EffectiveTensionStressVector;
        //noalias(rPredictiveStressVector) += (1.0 - mDamageCompression) * data.EffectiveCompressionStressVector;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::CheckDamageLoadingUnloading(
        bool& is_damaging_tension,
        bool& is_damaging_compression)
    {
        const double F_tension = mUniaxialStressTension - mCurrentThresholdTension;
        const double F_compression = mUniaxialStressCompression - mCurrentThresholdCompression;

        is_damaging_tension = (F_tension > 0.0)
            ? true : false;
        is_damaging_compression = (F_compression > 0.0)
            ? true : false;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::CalculateTangentTensor(
        Parameters& rValues,
        const array_1d<double, 3>& rStrainVector,
        const array_1d<double, 3>& rPredictiveStressVector,
        CalculationData& data,
        int IntegrationImplex)
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
        double save_uniaxial_stress_tension = mUniaxialStressTension;
        double save_uniaxial_stress_compression = mUniaxialStressCompression;

        // perturbation parameter
        double perturbation_factor = 1.0E-8;

        // perturbed vectors
        Vector perturbated_strain_vector(VoigtSize);
        Vector perturbated_stress_vector(VoigtSize);

        // apply perturbation to each strain component...
        for (size_t j = 0; j < VoigtSize; j++)
        {
            noalias(perturbated_strain_vector) = rStrainVector;

            perturbated_strain_vector(j) = rStrainVector(j) + perturbation_factor;
            int integration_implex = IntegrationImplex;
            //this->CalculateMaterialResponseInternal(
            //    perturbated_strain_vector,
            //    perturbated_stress_vector,
            //    data,
            //    integration_implex);

            for (size_t i = 0; i < VoigtSize; i++)
                constitutive_matrix(i, j) = (perturbated_stress_vector(i) - rPredictiveStressVector(i)) /
                perturbation_factor;
        }

        // restore internal variables
        mThresholdTension = save_threshold_tension;
        mThresholdCompression = save_threshold_compression;
        mDamageTension = save_damage_tension;
        mDamageCompression = save_damage_compression;
        mUniaxialStressTension = save_uniaxial_stress_tension;
        mUniaxialStressCompression = save_uniaxial_stress_compression;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamageDPlusDMinusPlaneStress2DLaw::CalculateSecantTensor(
        Parameters& rParameters,
        const Matrix& rElasticityMatrix,
        const CalculationData& data)
    {
        Matrix& constitutive_matrix = rParameters.GetConstitutiveMatrix();
        if (constitutive_matrix.size1() != VoigtSize || constitutive_matrix.size2() != VoigtSize)
            constitutive_matrix.resize(VoigtSize, VoigtSize);

        Matrix DamageMatrix(IdentityMatrix(3, 3));
        noalias(DamageMatrix) -= mDamageTension * data.ProjectionTensorTension;
        noalias(DamageMatrix) -= mDamageCompression * data.ProjectionTensorCompression;

        noalias(constitutive_matrix) = prod(DamageMatrix, rElasticityMatrix);
    }

} // namespace Kratos
