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

// System includes

// External includes

// Project includes
//#include "utilities/math_utils.h"
#include "masonry_orthotropic_damage_plane_stress_2d.h"

#include <cmath>

#define OPTIMIZE_CHARACTERISTIC_LENGTH
#define HEAVISIDE(X) ( X >= 0.0 ? 1.0 : 0.0)
#define MACAULAY(X)  ( X >= 0.0 ? X : 0.0)
//#define PROJECTION_OPERATOR_CERVERA_2003_ORTHORTOPIC_CL
//#define PROJECTION_OPERATOR_CERVERA_2017_ORTHORTOPIC_CL

namespace Kratos
{
    /// Has variable <double>
    bool MasonryOrthotropicDamagePlaneStress2DLaw::Has(
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
    double& MasonryOrthotropicDamagePlaneStress2DLaw::GetValue(
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
    void MasonryOrthotropicDamagePlaneStress2DLaw::SetValue(
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

    void MasonryOrthotropicDamagePlaneStress2DLaw::InitializeMaterial(
        const Properties& rProperties,
        const GeometryType& rGeometry,
        const Vector& rShapeFunctionsValues)
    {
        if (!mInitializeDamageLaw) {
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

            mInitializeDamageLaw = true;
        }
    }

    void MasonryOrthotropicDamagePlaneStress2DLaw::FinalizeMaterialResponseCauchy(
        Parameters& rParameters)
    {
        // save converged values
        mCurrentThresholdTension = mThresholdTension;
        mCurrentThresholdCompression = mThresholdCompression;
    }

    void MasonryOrthotropicDamagePlaneStress2DLaw::CalculateMaterialResponseCauchy(
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

        this->CalculateMaterialResponseInternal(
            r_stress_vector,
            data);

        bool is_damaging_tension = false;
        bool is_damaging_compression = false;
        this->CheckDamageLoadingUnloading(is_damaging_tension, is_damaging_compression);

        // Computation of the Constitutive Tensor
        if (rParameters.GetOptions().Is(COMPUTE_CONSTITUTIVE_TENSOR)) {
            if (is_damaging_tension || is_damaging_compression) {
                this->CalculateTangentTensor(rParameters);
            }
            else {
                this->CalculateSecantTensor(rParameters, ElasticityMatrix, data);
            }
        }
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamagePlaneStress2DLaw::ResetMaterial(
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
        mInitializeDamageLaw = false;

    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamagePlaneStress2DLaw::GetLawFeatures(
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
    int MasonryOrthotropicDamagePlaneStress2DLaw::Check(
        const Properties& rProperties,
        const GeometryType& rGeometry,
        const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        KRATOS_ERROR_IF(rProperties.NumberOfSubproperties() != SizeType(2))
            << "MasonryOrthotropicDamagePlaneStress2DLaw requires 2 sub properties.";

        // ELASTIC PARAMETERS
        CheckOrthotropicParameter(*rProperties.GetSubProperties().begin());
        CheckOrthotropicParameter(*(rProperties.GetSubProperties().begin() + 1));

        return 0;

        KRATOS_CATCH("");
    }

    int MasonryOrthotropicDamagePlaneStress2DLaw::CheckOrthotropicParameter(
        const Properties& rProperties)
    {
        KRATOS_ERROR_IF_NOT(rProperties.Has(YOUNG_MODULUS)) << "Missing variable: YOUNG_MODULUS";
        KRATOS_ERROR_IF_NOT(rProperties.Has(POISSON_RATIO)) << "Missing variable: POISSON_RATIO";
        KRATOS_ERROR_IF_NOT(rProperties.Has(SHEAR_MODULUS)) << "Missing variable: SHEAR_MODULUS";
        // TENSION
        KRATOS_ERROR_IF_NOT(rProperties.Has(YIELD_STRESS_TENSION)) << "Missing variable: YIELD_STRESS_TENSION";
        KRATOS_ERROR_IF_NOT(rProperties.Has(FRACTURE_ENERGY_TENSION)) << "Missing variable: FRACTURE_ENERGY_TENSION";
        // COMPRESSION
        KRATOS_ERROR_IF_NOT(rProperties.Has(DAMAGE_ONSET_STRESS_COMPRESSION)) << "Missing variable: DAMAGE_ONSET_STRESS_COMPRESSION";
        KRATOS_ERROR_IF_NOT(rProperties.Has(YIELD_STRESS_COMPRESSION)) << "Missing variable: YIELD_STRESS_COMPRESSION";
        KRATOS_ERROR_IF_NOT(rProperties.Has(YIELD_STRAIN_COMPRESSION)) << "Missing variable: YIELD_STRAIN_COMPRESSION";
        KRATOS_ERROR_IF_NOT(rProperties.Has(RESIDUAL_STRESS_COMPRESSION)) << "Missing variable: RESIDUAL_STRESS_COMPRESSION";
        KRATOS_ERROR_IF_NOT(rProperties.Has(FRACTURE_ENERGY_COMPRESSION)) << "Missing variable: FRACTURE_ENERGY_COMPRESSION";
        KRATOS_ERROR_IF_NOT(rProperties.Has(BEZIER_CONTROLLER_C1)) << "Missing variable: BEZIER_CONTROLLER_C1";
        KRATOS_ERROR_IF_NOT(rProperties.Has(BEZIER_CONTROLLER_C2)) << "Missing variable: BEZIER_CONTROLLER_C2";
        KRATOS_ERROR_IF_NOT(rProperties.Has(BEZIER_CONTROLLER_C3)) << "Missing variable: BEZIER_CONTROLLER_C3";
        // SHEAR
        KRATOS_ERROR_IF_NOT(rProperties.Has(YIELD_STRESS_SHEAR_TENSION)) << "Missing variable: YIELD_STRESS_SHEAR_TENSION";
        KRATOS_ERROR_IF_NOT(rProperties.Has(YIELD_STRESS_SHEAR_COMPRESSION)) << "Missing variable: YIELD_STRESS_SHEAR_COMPRESSION";
        KRATOS_ERROR_IF_NOT(rProperties.Has(SHEAR_COMPRESSION_REDUCTOR)) << "Missing variable: SHEAR_COMPRESSION_REDUCTOR";
        // BIAXIAL
        KRATOS_ERROR_IF_NOT(rProperties.Has(BIAXIAL_COMPRESSION_MULTIPLIER)) << "Missing variable: BIAXIAL_COMPRESSION_MULTIPLIER";

        return 0;
    }

    /***********************************************************************************/
    /***********************************************************************************/
    MasonryOrthotropicDamagePlaneStress2DLaw::CalculationData MasonryOrthotropicDamagePlaneStress2DLaw::GetCalculationData(
        const Properties& rProperties,
        const GeometryType& rGeometry,
        const ProcessInfo& rProcessInfo)
    {
        const DirectionalMaterialProperties MaterialProperties1(
            this->ComputeCharacteristicLength(
                rGeometry, 0),
            *rProperties.GetSubProperties().begin());

        const DirectionalMaterialProperties MaterialProperties2(
            this->ComputeCharacteristicLength(
                rGeometry, 1),
            *(rProperties.GetSubProperties().begin() + 1));

        CalculationData data(MaterialProperties1, MaterialProperties2);

        // Misc
        data.DeltaTime = rProcessInfo[DELTA_TIME];
        data.TensionYieldModel = rProperties.Has(TENSION_YIELD_MODEL)
            ? rProperties[TENSION_YIELD_MODEL]
            : 0;

        return data;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamagePlaneStress2DLaw::CalculateOrthotropicElasticityMatrix(
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

    double MasonryOrthotropicDamagePlaneStress2DLaw::CalculateDamageAngle(
        const array_1d<double, 3>& rStressVector)
    {
        array_1d<double, 2> local_stress_1, local_coordinate_1;
        local_stress_1[0] = rStressVector[0];
        local_stress_1[1] = rStressVector[2]/2;
        local_coordinate_1[0] = 0;
        local_coordinate_1[1] = 1;
        if (norm_2(local_stress_1) < 1e-2 && sqrt(pow(rStressVector[1], 2) + pow(rStressVector[2] / 2, 2)) > 1e-2) {
            local_stress_1[0] = rStressVector[2]/2;
            local_stress_1[1] = rStressVector[1];
            local_coordinate_1[0] = 1;
            local_coordinate_1[1] = 0;
        }
        else return 0;

        return acos(inner_prod(local_stress_1, local_coordinate_1)/(norm_2(local_stress_1)));

    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamagePlaneStress2DLaw::AssembleTransformationMatrix(
        const DirectionalMaterialProperties& rMaterialProperties1,
        const DirectionalMaterialProperties& rMaterialProperties2,
        TransformationMatrices& rTransformationMatrices)
    {
        rTransformationMatrices.TransformationMatrixTension(0, 0) = 1;
        rTransformationMatrices.TransformationMatrixTension(1, 1) = rMaterialProperties1.YieldStressTension / rMaterialProperties2.YieldStressTension;
        rTransformationMatrices.TransformationMatrixTension(2, 2) = rMaterialProperties1.YieldStressShearTension / rMaterialProperties1.YieldStressShearTension;

        rTransformationMatrices.InverseTransformationMatrixTension(0, 0) = 1;
        rTransformationMatrices.InverseTransformationMatrixTension(1, 1) = rMaterialProperties2.YieldStressTension / rMaterialProperties1.YieldStressTension;
        rTransformationMatrices.InverseTransformationMatrixTension(2, 2) = rMaterialProperties1.YieldStressShearTension / rMaterialProperties1.YieldStressShearTension;

        rTransformationMatrices.TransformationMatrixCompression(0, 0) = 1;
        rTransformationMatrices.TransformationMatrixCompression(1, 1) = rMaterialProperties1.YieldStressCompression / rMaterialProperties2.YieldStressCompression;
        rTransformationMatrices.TransformationMatrixCompression(2, 2) = rMaterialProperties1.YieldStressShearCompression / rMaterialProperties2.YieldStressShearCompression;

        rTransformationMatrices.InverseTransformationMatrixCompression(0, 0) = 1;
        rTransformationMatrices.InverseTransformationMatrixCompression(1, 1) = rMaterialProperties2.YieldStressCompression / rMaterialProperties1.YieldStressCompression;
        rTransformationMatrices.InverseTransformationMatrixCompression(2, 2) =  rMaterialProperties2.YieldStressShearCompression / rMaterialProperties1.YieldStressShearCompression;

    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamagePlaneStress2DLaw::CalculateProjectedIsotropicMaterial(
        const DirectionalMaterialProperties& rMaterialProperties1,
        const DirectionalMaterialProperties& rMaterialProperties2,
        double AngleToDamage,
        DirectionalMaterialProperties& rProjectedProperties)
    {
        const double material_length_1_tension = 2.0 * rMaterialProperties1.E * rMaterialProperties1.FractureEnergyTension /
            pow(rMaterialProperties1.YieldStressTension, 2);
        const double material_length_2_tension = 2.0 * rMaterialProperties2.E * rMaterialProperties2.FractureEnergyTension /
            pow(rMaterialProperties2.YieldStressTension, 2);

        const double material_length_projected_tension = sqrt(1 / (
            (1 / pow(material_length_1_tension, 2)) * pow(cos(AngleToDamage), 2)
            + (1 / pow(material_length_2_tension, 2)) * pow(sin(AngleToDamage), 2)));

        rProjectedProperties.FractureEnergyTension = (pow(rProjectedProperties.YieldStressTension, 2) / (2 * rProjectedProperties.E)) * material_length_projected_tension;

        const double material_length_1_compression = 2.0 * rMaterialProperties1.E * rMaterialProperties1.FractureEnergyCompression /
            pow(rMaterialProperties1.YieldStressCompression, 2);
        const double material_length_2_compression = 2.0 * rMaterialProperties2.E * rMaterialProperties2.FractureEnergyCompression /
            pow(rMaterialProperties2.YieldStressCompression, 2);

        const double material_length_projected_compression = sqrt(1 / (
            (1 / pow(material_length_1_compression, 2)) * pow(cos(AngleToDamage), 2)
            + (1 / pow(material_length_2_compression, 2)) * pow(sin(AngleToDamage), 2)));

        // Fracture Energy
        rProjectedProperties.FractureEnergyCompression = (pow(rProjectedProperties.YieldStressCompression, 2) / (2 * rProjectedProperties.E)) * material_length_projected_compression;

        // Elastic Limit Stress Compression
        rProjectedProperties.ElasticLimitStressCompression = (rProjectedProperties.YieldStressCompression / rMaterialProperties1.YieldStressCompression)
            * rMaterialProperties1.ElasticLimitStressCompression * pow(cos(AngleToDamage), 2) +
            (rProjectedProperties.YieldStressCompression / rMaterialProperties2.YieldStressCompression)
            * rMaterialProperties2.ElasticLimitStressCompression * pow(sin(AngleToDamage), 2);

        // Residual Stress Compression
        rProjectedProperties.ResidualStressCompression = (rProjectedProperties.YieldStressCompression / rMaterialProperties1.YieldStressCompression)
            * rMaterialProperties1.ResidualStressCompression * pow(cos(AngleToDamage), 2) +
            (rProjectedProperties.YieldStressCompression / rMaterialProperties2.YieldStressCompression)
            * rMaterialProperties2.ResidualStressCompression * pow(sin(AngleToDamage), 2);

        // Yield Strain
        rProjectedProperties.YieldStrainCompression = rMaterialProperties1.YieldStrainCompression * pow(cos(AngleToDamage), 2) +
            ((rProjectedProperties.ResidualStressCompression / rProjectedProperties.E)/ (rMaterialProperties2.ResidualStressCompression / rMaterialProperties2.E))
            * rMaterialProperties2.YieldStrainCompression * pow(sin(AngleToDamage), 2);

        // C1
        rProjectedProperties.BezierControllerC1 = rMaterialProperties1.BezierControllerC1 * pow(cos(AngleToDamage), 2) +
            rMaterialProperties2.BezierControllerC1 * pow(sin(AngleToDamage), 2);

        // C2
        rProjectedProperties.BezierControllerC2 = rMaterialProperties1.BezierControllerC2 * pow(cos(AngleToDamage), 2) +
            rMaterialProperties2.BezierControllerC2 * pow(sin(AngleToDamage), 2);

        // C3
        rProjectedProperties.BezierControllerC3 = rMaterialProperties1.BezierControllerC3 * pow(cos(AngleToDamage), 2) +
            rMaterialProperties2.BezierControllerC3 * pow(sin(AngleToDamage), 2);

        // ShearCompressionReductor
        rProjectedProperties.ShearCompressionReductor = rMaterialProperties1.ShearCompressionReductor * pow(cos(AngleToDamage), 2) +
            rMaterialProperties2.ShearCompressionReductor * pow(sin(AngleToDamage), 2);

        // ShearCompressionReductor
        rProjectedProperties.BiaxialCompressionMultiplier = rMaterialProperties1.BiaxialCompressionMultiplier * pow(cos(AngleToDamage), 2) +
            rMaterialProperties2.BiaxialCompressionMultiplier * pow(sin(AngleToDamage), 2);
    }
    /***********************************************************************************/
    /***********************************************************************************/

    void MasonryOrthotropicDamagePlaneStress2DLaw::TensionCompressionSplit(
        const array_1d<double, 3>& rEffectiveStressVector,
        array_1d<double, 2>& rPrincipalStressVector,
        array_1d<double, 3>& rEffectiveTensionStressVector,
        array_1d<double, 3>& rEffectiveCompressionStressVector)
    {
        AdvancedConstitutiveLawUtilities<3>::CalculatePrincipalStresses(
            rPrincipalStressVector, rEffectiveStressVector);
        AdvancedConstitutiveLawUtilities<3>::SpectralDecomposition(
            rEffectiveStressVector, rEffectiveTensionStressVector, rEffectiveCompressionStressVector);
    } 

    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamagePlaneStress2DLaw::ConstructProjectionTensors(
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

#ifdef PROJECTION_OPERATOR_CERVERA_2003_ORTHORTOPIC_CL
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
#endif //PROJECTION_OPERATOR_CERVERA_2003_ORTHORTOPIC_CL

#ifdef PROJECTION_OPERATOR_CERVERA_2017_ORTHORTOPIC_CL
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

        noalias(rProjectionTensorTension) += (HEAVISIDE(eigen_values_matrix(0, 0)) + HEAVISIDE(eigen_values_matrix(1, 1))) *
            outer_prod(projection_vector_cross, projection_vector_cross);
#endif //PROJECTION_OPERATOR_CERVERA_2017_ORTHORTOPIC_CL

        noalias(rProjectionTensorCompression) = IdentityMatrix(3, 3) - rProjectionTensorTension;
    }

    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamagePlaneStress2DLaw::CalculateEquivalentStressTension(
        const CalculationData& data,
        const DirectionalMaterialProperties& rMaterialProperties,
        const array_1d<double, 3> rEffectiveStressVector,
        const array_1d<double, 2> rPrincipalStressVector,
        double& rUniaxialStressTension) const
    {
        rUniaxialStressTension = 0.0;
        if (rPrincipalStressVector(0) > 0.0) {
            if (data.TensionYieldModel == 0) {
                // Lubliner Yield Criteria 
                const double yield_compression = rMaterialProperties.YieldStressCompression;
                const double yield_tension = rMaterialProperties.YieldStressTension;
                const double alpha = (rMaterialProperties.BiaxialCompressionMultiplier - 1.0) /
                    (2.0 * rMaterialProperties.BiaxialCompressionMultiplier - 1.0);
                double I1, J2;
                array_1d<double, 3> deviator = ZeroVector(3);

                ConstitutiveLawUtilities<3>::CalculateI1Invariant(rEffectiveStressVector, I1);
                ConstitutiveLawUtilities<3>::CalculateJ2Invariant(rEffectiveStressVector, I1, deviator, J2);

                const double beta = yield_compression / yield_tension * (1.0 - alpha) - (1.0 + alpha);
                const double smax = std::max(std::max(rPrincipalStressVector(0), rPrincipalStressVector(1)), 0.0);

                rUniaxialStressTension = 1.0 / (1.0 - alpha) * (alpha * I1 + std::sqrt(3.0 * J2) + beta * smax) /
                    yield_compression * yield_tension;
            }
            else if (data.TensionYieldModel == 1) {
                // Rankine Yield Criteria
                rUniaxialStressTension = std::max(std::max(rPrincipalStressVector(0), rPrincipalStressVector(1)), 0.0);
            }
        }
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamagePlaneStress2DLaw::CalculateEquivalentStressCompression(
        const CalculationData& data,
        const DirectionalMaterialProperties& rMaterialProperties,
        const array_1d<double, 3> rEffectiveStressVector,
        const array_1d<double, 2> rPrincipalStressVector,
        double& UniaxialStressCompression) const
    {
        UniaxialStressCompression = 0.0;
        if (rPrincipalStressVector(1) < 0.0) {
            const double yield_compression = rMaterialProperties.YieldStressCompression;
            const double yield_tension = rMaterialProperties.YieldStressTension;
            const double alpha = (rMaterialProperties.BiaxialCompressionMultiplier - 1.0) /
                (2.0 * rMaterialProperties.BiaxialCompressionMultiplier - 1.0);
            double I1, J2;
            array_1d<double, 3> deviator = ZeroVector(3);

            ConstitutiveLawUtilities<3>::CalculateI1Invariant(rEffectiveStressVector, I1);
            ConstitutiveLawUtilities<3>::CalculateJ2Invariant(rEffectiveStressVector, I1, deviator, J2);

            const double beta = (yield_compression / yield_tension) * (1.0 - alpha) - (1.0 + alpha);
            const double smax = std::max(std::max(rPrincipalStressVector(0), rPrincipalStressVector(1)), 0.0);

            UniaxialStressCompression = 1.0 / (1.0 - alpha) * (alpha * I1 + std::sqrt(3.0 * J2) +
                rMaterialProperties.ShearCompressionReductor * beta * smax);
        }
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamagePlaneStress2DLaw::CalculateDamageTension(
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

            KRATOS_WARNING_IF("::[MasonryOrthotropicDamagePlaneStress2DLaw]::CalculateDamageCompression", rMaterialProperties.CharacteristicLength >= material_length)
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
    void MasonryOrthotropicDamagePlaneStress2DLaw::CalculateDamageCompression(
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
            const double e_0 = s_0 / young_modulus;
            const double e_i = s_p / young_modulus;

            // current abscissa
            const double strain_like_counterpart = TresholdCompression / young_modulus;

            // compute damage
            double damage_variable = TresholdCompression;
            if (strain_like_counterpart <= e_p) {
                KRATOS_WARNING_IF("::[MasonryOrthotropicDamagePlaneStress2DLaw]::CalculateDamageCompression", e_0 > e_p)
                    << "Residual strain e_p: " << e_p << " is smaller than the elastic extension e_0: " << e_0 << "." << std::endl;

                this->EvaluateBezierCurve(damage_variable, strain_like_counterpart, e_0, e_i, e_p, s_0, s_p, s_p);
            } else {
                const double s_k = s_r + (s_p - s_r) * c1;
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

                KRATOS_WARNING_IF("::[MasonryOrthotropicDamagePlaneStress2DLaw]::CalculateDamageCompression", stretcher <= -1.0)
                    << "Compressive fracture energy is too low" << std::endl
                    << "Input Gc/lch = " << specific_fracture_energy << std::endl
                    << "Minimum Gc to avoid constitutive snap-back = " << bezier_energy_1 << std::endl;

                this->ApplyBezierStretcherToStrains(stretcher, e_p, e_j, e_k, e_r, e_u);


                if (strain_like_counterpart <= e_k) {
                    this->EvaluateBezierCurve(damage_variable, strain_like_counterpart, e_p, e_j, e_k, s_p, s_p, s_k);
                }
                else if (strain_like_counterpart <= e_u) {
                    this->EvaluateBezierCurve(damage_variable, strain_like_counterpart, e_k, e_r, e_u, s_k, s_r, s_r);
                }
                else {
                    damage_variable = s_r;
                }
            }
            rDamageCompression = 1.0 - damage_variable / TresholdCompression;
        }
    }
    /***********************************************************************************/
    /***********************************************************************************/
    inline void MasonryOrthotropicDamagePlaneStress2DLaw::ComputeBezierEnergy(
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
    inline double MasonryOrthotropicDamagePlaneStress2DLaw::EvaluateBezierArea(
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
    inline void MasonryOrthotropicDamagePlaneStress2DLaw::ApplyBezierStretcherToStrains(
        double stretcher, double e_p, double& e_j, double& e_k, double& e_r, double& e_u) const
    {
        e_j += (e_j - e_p) * stretcher;
        e_k += (e_k - e_p) * stretcher;
        e_r += (e_r - e_p) * stretcher;
        e_u += (e_u - e_p) * stretcher;
    }
    /***********************************************************************************/
    /***********************************************************************************/
    inline void MasonryOrthotropicDamagePlaneStress2DLaw::EvaluateBezierCurve(
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
    double MasonryOrthotropicDamagePlaneStress2DLaw::ComputeCharacteristicLength(
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
    void MasonryOrthotropicDamagePlaneStress2DLaw::CalculateMaterialResponseInternal(
        Vector& rPredictiveStressVector,
        CalculationData& data)
    {
        AdvancedConstitutiveLawUtilities<3>::SpectralDecomposition(
            data.EffectiveStressVector, data.EffectiveTensionStressVector, data.EffectiveCompressionStressVector);
        AdvancedConstitutiveLawUtilities<3>::CalculatePrincipalStresses(
            data.PrincipalStressVector, data.EffectiveStressVector);
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
        array_1d<double, 2> principal_stresses_mapped_isotropic;
        array_1d<double, 3> effective_stress_mapped_isotropic = effective_tension_stress_mapped_isotropic + effective_compression_stress_mapped_isotropic;
        AdvancedConstitutiveLawUtilities<3>::CalculatePrincipalStresses(
            principal_stresses_mapped_isotropic, effective_stress_mapped_isotropic);

        double angle_to_damage = CalculateDamageAngle(effective_stress_mapped_isotropic);

        DirectionalMaterialProperties projected_material(data.MaterialProperties1);
        CalculateProjectedIsotropicMaterial(
            data.MaterialProperties1, data.MaterialProperties2,
            angle_to_damage, projected_material);

        // compute the equivalent stress measures
        this->CalculateEquivalentStressTension(data,
            projected_material, effective_stress_mapped_isotropic, principal_stresses_mapped_isotropic, mUniaxialStressTension);
        this->CalculateEquivalentStressCompression(data,
            projected_material, effective_stress_mapped_isotropic, principal_stresses_mapped_isotropic, mUniaxialStressCompression);

        if (mUniaxialStressTension > mThresholdTension)
            mThresholdTension = mUniaxialStressTension;
        this->CalculateDamageTension(projected_material, mThresholdTension, mDamageTension);

        if (mUniaxialStressCompression > mThresholdCompression)
            mThresholdCompression = mUniaxialStressCompression;
        this->CalculateDamageCompression(projected_material, mThresholdCompression, mDamageCompression);

        mCurrentThresholdTension = mThresholdTension;
        mCurrentThresholdCompression = mThresholdCompression;

        //mDamageTension = std::max(mDamageTension, mDamageCompression);

        // Compute the stresses to isotropic space
        array_1d<double, 3> total_tension_stress_mapped_isotropic = (1.0 - mDamageTension) * effective_tension_stress_mapped_isotropic;
        array_1d<double, 3> total_compression_stress_mapped_isotropic = (1.0 - mDamageCompression) * effective_compression_stress_mapped_isotropic;

        // Back rotation of final Predictive Stress
        noalias(rPredictiveStressVector) = prod(transformation_matrices.InverseTransformationMatrixTension, total_tension_stress_mapped_isotropic);
        noalias(rPredictiveStressVector) += prod(transformation_matrices.InverseTransformationMatrixCompression, total_compression_stress_mapped_isotropic);

    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamagePlaneStress2DLaw::CheckDamageLoadingUnloading(
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
    void MasonryOrthotropicDamagePlaneStress2DLaw::CalculateTangentTensor(
        Parameters& rValues)
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();

        const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD)
            ? r_material_properties[CONSIDER_PERTURBATION_THRESHOLD]
            : true;
        const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION)
            ? static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION])
            : TangentOperatorEstimation::FirstOrderPerturbation;

        if (tangent_operator_estimation == TangentOperatorEstimation::Analytic) {
            // Already stored in rValues.GetConstitutiveMatrix()...
        }
        else if (tangent_operator_estimation == TangentOperatorEstimation::FirstOrderPerturbation) {
            // Calculates the Tangent Constitutive Tensor by perturbation (first order)
            TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 1);
        }
        else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbation) {
            // Calculates the Tangent Constitutive Tensor by perturbation (second order)
            TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 2);
        }
        else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbationV2) {
            // Calculates the Tangent Constitutive Tensor by perturbation (second order)
            TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 4);
        }
    }
    /***********************************************************************************/
    /***********************************************************************************/
    void MasonryOrthotropicDamagePlaneStress2DLaw::CalculateSecantTensor(
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
