// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//

#include "custom_constitutive/reference_3d_mohr_coulomb_plane_strain_law.h"

#include "constitutive_laws_application_variables.h"
#include "geo_mechanics_application_variables.h"
#include "includes/kratos_components.h"

namespace Kratos
{
namespace
{

constexpr const char* ReferenceLawName =
    "SmallStrainIsotropicPlasticity3DMohrCoulombMohrCoulomb";

Properties CreateReferenceProperties(const Properties& rGeoProperties)
{
    auto reference_properties = Properties{rGeoProperties};

    if (rGeoProperties.Has(GEO_COHESION)) {
        reference_properties.SetValue(COHESION, rGeoProperties[GEO_COHESION]);
    }
    if (rGeoProperties.Has(GEO_FRICTION_ANGLE)) {
        reference_properties.SetValue(FRICTION_ANGLE, rGeoProperties[GEO_FRICTION_ANGLE]);
    }
    if (rGeoProperties.Has(GEO_DILATANCY_ANGLE)) {
        reference_properties.SetValue(DILATANCY_ANGLE, rGeoProperties[GEO_DILATANCY_ANGLE]);
    }
    if (rGeoProperties.Has(GEO_TENSILE_STRENGTH)) {
        reference_properties.SetValue(YIELD_STRESS, rGeoProperties[GEO_TENSILE_STRENGTH]);
    }
    if (!rGeoProperties.Has(DENSITY) && rGeoProperties.Has(DENSITY_SOLID)) {
        reference_properties.SetValue(DENSITY, rGeoProperties[DENSITY_SOLID]);
    }

    return reference_properties;
}

void CopyCommonParameters(ConstitutiveLaw::Parameters& rDestination,
                          ConstitutiveLaw::Parameters& rSource)
{
    rDestination.SetOptions(rSource.GetOptions());

    if (rSource.IsSetShapeFunctionsValues()) {
        rDestination.SetShapeFunctionsValues(rSource.GetShapeFunctionsValues());
    }
    if (rSource.IsSetShapeFunctionsDerivatives()) {
        rDestination.SetShapeFunctionsDerivatives(rSource.GetShapeFunctionsDerivatives());
    }
    if (rSource.IsSetDeformationGradientF()) {
        rDestination.SetDeformationGradientF(rSource.GetDeformationGradientF());
    }
    if (rSource.IsSetDeterminantF()) {
        rDestination.SetDeterminantF(rSource.GetDeterminantF());
    }
}

} // namespace

ConstitutiveLaw::Pointer Reference3DMohrCoulombPlaneStrainLaw::Clone() const
{
    auto p_result = Kratos::make_shared<Reference3DMohrCoulombPlaneStrainLaw>();
    if (mpReferenceLaw) {
        p_result->mpReferenceLaw = mpReferenceLaw->Clone();
    }
    p_result->mHasReferenceStrengthParameters = mHasReferenceStrengthParameters;
    p_result->mLastReferenceCohesion          = mLastReferenceCohesion;
    p_result->mLastReferenceFrictionAngle     = mLastReferenceFrictionAngle;
    return p_result;
}

ConstitutiveLaw::SizeType Reference3DMohrCoulombPlaneStrainLaw::WorkingSpaceDimension()
{
    return 2;
}

ConstitutiveLaw::SizeType Reference3DMohrCoulombPlaneStrainLaw::GetStrainSize() const
{
    return GeoPlaneStrainSize;
}

ConstitutiveLaw::StrainMeasure Reference3DMohrCoulombPlaneStrainLaw::GetStrainMeasure()
{
    return ConstitutiveLaw::StrainMeasure_Infinitesimal;
}

ConstitutiveLaw::StressMeasure Reference3DMohrCoulombPlaneStrainLaw::GetStressMeasure()
{
    return ConstitutiveLaw::StressMeasure_Cauchy;
}

void Reference3DMohrCoulombPlaneStrainLaw::GetLawFeatures(Features& rFeatures)
{
    auto options = Flags{};
    options.Set(ConstitutiveLaw::PLANE_STRAIN_LAW);
    options.Set(ConstitutiveLaw::INFINITESIMAL_STRAINS);
    options.Set(ConstitutiveLaw::ISOTROPIC);

    rFeatures.SetOptions(options);
    rFeatures.SetStrainMeasure(GetStrainMeasure());
    rFeatures.SetStrainSize(GetStrainSize());
    rFeatures.SetSpaceDimension(WorkingSpaceDimension());
}

void Reference3DMohrCoulombPlaneStrainLaw::InitializeMaterial(
    const Properties& rMaterialProperties,
    const Geometry<Node>& rElementGeometry,
    const Vector& rShapeFunctionsValues)
{
    EnsureReferenceLawIsCreated();

    const auto reference_properties = CreateReferenceProperties(rMaterialProperties);
    SynchronizeReferenceStrengthParameters(reference_properties, rElementGeometry,
                                           rShapeFunctionsValues);
}

bool Reference3DMohrCoulombPlaneStrainLaw::RequiresInitializeMaterialResponse()
{
    EnsureReferenceLawIsCreated();
    return mpReferenceLaw->RequiresInitializeMaterialResponse();
}

bool Reference3DMohrCoulombPlaneStrainLaw::RequiresFinalizeMaterialResponse()
{
    EnsureReferenceLawIsCreated();
    return mpReferenceLaw->RequiresFinalizeMaterialResponse();
}

void Reference3DMohrCoulombPlaneStrainLaw::InitializeMaterialResponseCauchy(
    Parameters& rValues)
{
    EnsureReferenceLawIsCreated();

    auto strain_vector_3d = MapPlaneStrainVectorTo3D(rValues.GetStrainVector());
    auto stress_vector_3d = rValues.IsSetStressVector()
                                ? MapPlaneStrainVectorTo3D(rValues.GetStressVector())
                                : Vector{ZeroVector(Reference3DSize)};
    auto constitutive_matrix_3d = Matrix{ZeroMatrix(Reference3DSize, Reference3DSize)};

    const auto reference_properties = CreateReferenceProperties(rValues.GetMaterialProperties());
    const auto shape_functions_values = rValues.IsSetShapeFunctionsValues()
                                            ? rValues.GetShapeFunctionsValues()
                                            : Vector{};
    SynchronizeReferenceStrengthParameters(reference_properties, rValues.GetElementGeometry(),
                                           shape_functions_values);
    Parameters reference_values(rValues.GetElementGeometry(), reference_properties,
                                rValues.GetProcessInfo());
    CopyCommonParameters(reference_values, rValues);
    reference_values.SetStrainVector(strain_vector_3d);
    reference_values.SetStressVector(stress_vector_3d);
    reference_values.SetConstitutiveMatrix(constitutive_matrix_3d);

    mpReferenceLaw->InitializeMaterialResponseCauchy(reference_values);
}

void Reference3DMohrCoulombPlaneStrainLaw::CalculateMaterialResponseCauchy(
    Parameters& rValues)
{
    EnsureReferenceLawIsCreated();

    auto strain_vector_3d = MapPlaneStrainVectorTo3D(rValues.GetStrainVector());
    auto stress_vector_3d = rValues.IsSetStressVector()
                                ? MapPlaneStrainVectorTo3D(rValues.GetStressVector())
                                : Vector{ZeroVector(Reference3DSize)};
    auto constitutive_matrix_3d = Matrix{ZeroMatrix(Reference3DSize, Reference3DSize)};

    const auto reference_properties = CreateReferenceProperties(rValues.GetMaterialProperties());
    const auto shape_functions_values = rValues.IsSetShapeFunctionsValues()
                                            ? rValues.GetShapeFunctionsValues()
                                            : Vector{};
    SynchronizeReferenceStrengthParameters(reference_properties, rValues.GetElementGeometry(),
                                           shape_functions_values);
    Parameters reference_values(rValues.GetElementGeometry(), reference_properties,
                                rValues.GetProcessInfo());
    CopyCommonParameters(reference_values, rValues);
    reference_values.SetStrainVector(strain_vector_3d);
    reference_values.SetStressVector(stress_vector_3d);
    reference_values.SetConstitutiveMatrix(constitutive_matrix_3d);

    mpReferenceLaw->CalculateMaterialResponseCauchy(reference_values);

    if (rValues.IsSetStressVector()) {
        rValues.GetStressVector() = Map3DVectorToPlaneStrain(stress_vector_3d);
    }
    if (rValues.IsSetConstitutiveMatrix()) {
        rValues.GetConstitutiveMatrix() = Map3DMatrixToPlaneStrain(constitutive_matrix_3d);
    }
}

void Reference3DMohrCoulombPlaneStrainLaw::FinalizeMaterialResponseCauchy(
    Parameters& rValues)
{
    EnsureReferenceLawIsCreated();

    auto strain_vector_3d = MapPlaneStrainVectorTo3D(rValues.GetStrainVector());
    auto stress_vector_3d = rValues.IsSetStressVector()
                                ? MapPlaneStrainVectorTo3D(rValues.GetStressVector())
                                : Vector{ZeroVector(Reference3DSize)};
    auto constitutive_matrix_3d = Matrix{ZeroMatrix(Reference3DSize, Reference3DSize)};

    const auto reference_properties = CreateReferenceProperties(rValues.GetMaterialProperties());
    const auto shape_functions_values = rValues.IsSetShapeFunctionsValues()
                                            ? rValues.GetShapeFunctionsValues()
                                            : Vector{};
    SynchronizeReferenceStrengthParameters(reference_properties, rValues.GetElementGeometry(),
                                           shape_functions_values);
    Parameters reference_values(rValues.GetElementGeometry(), reference_properties,
                                rValues.GetProcessInfo());
    CopyCommonParameters(reference_values, rValues);
    reference_values.SetStrainVector(strain_vector_3d);
    reference_values.SetStressVector(stress_vector_3d);
    reference_values.SetConstitutiveMatrix(constitutive_matrix_3d);

    mpReferenceLaw->FinalizeMaterialResponseCauchy(reference_values);

    if (rValues.IsSetStressVector()) {
        rValues.GetStressVector() = Map3DVectorToPlaneStrain(stress_vector_3d);
    }
    if (rValues.IsSetConstitutiveMatrix()) {
        rValues.GetConstitutiveMatrix() = Map3DMatrixToPlaneStrain(constitutive_matrix_3d);
    }
}

bool Reference3DMohrCoulombPlaneStrainLaw::Has(const Variable<double>& rVariable)
{
    EnsureReferenceLawIsCreated();
    return mpReferenceLaw->Has(rVariable);
}

bool Reference3DMohrCoulombPlaneStrainLaw::Has(const Variable<Vector>& rVariable)
{
    EnsureReferenceLawIsCreated();
    return mpReferenceLaw->Has(rVariable);
}

double& Reference3DMohrCoulombPlaneStrainLaw::GetValue(
    const Variable<double>& rVariable,
    double& rValue)
{
    EnsureReferenceLawIsCreated();
    return mpReferenceLaw->GetValue(rVariable, rValue);
}

Vector& Reference3DMohrCoulombPlaneStrainLaw::GetValue(
    const Variable<Vector>& rVariable,
    Vector& rValue)
{
    EnsureReferenceLawIsCreated();

    Vector reference_value;
    mpReferenceLaw->GetValue(rVariable, reference_value);
    rValue = reference_value.size() == Reference3DSize
                 ? Map3DVectorToPlaneStrain(reference_value)
                 : reference_value;
    return rValue;
}

void Reference3DMohrCoulombPlaneStrainLaw::EnsureReferenceLawIsCreated()
{
    if (mpReferenceLaw) {
        return;
    }

    KRATOS_ERROR_IF_NOT(KratosComponents<ConstitutiveLaw>::Has(ReferenceLawName))
        << "The reference constitutive law \"" << ReferenceLawName
        << "\" is not registered." << std::endl;

    mpReferenceLaw = KratosComponents<ConstitutiveLaw>::Get(ReferenceLawName).Clone();
}

void Reference3DMohrCoulombPlaneStrainLaw::SynchronizeReferenceStrengthParameters(
    const Properties& rReferenceProperties,
    const Geometry<Node>& rElementGeometry,
    const Vector& rShapeFunctionsValues)
{
    KRATOS_ERROR_IF_NOT(rReferenceProperties.Has(COHESION))
        << "The translated reference properties do not contain COHESION." << std::endl;
    KRATOS_ERROR_IF_NOT(rReferenceProperties.Has(FRICTION_ANGLE))
        << "The translated reference properties do not contain FRICTION_ANGLE." << std::endl;

    const double cohesion       = rReferenceProperties[COHESION];
    const double friction_angle = rReferenceProperties[FRICTION_ANGLE];
    const bool strength_has_changed =
        !mHasReferenceStrengthParameters || cohesion != mLastReferenceCohesion ||
        friction_angle != mLastReferenceFrictionAngle;

    if (!strength_has_changed) {
        return;
    }

    mpReferenceLaw->InitializeMaterial(rReferenceProperties, rElementGeometry,
                                       rShapeFunctionsValues);
    mLastReferenceCohesion          = cohesion;
    mLastReferenceFrictionAngle     = friction_angle;
    mHasReferenceStrengthParameters = true;
}

Vector Reference3DMohrCoulombPlaneStrainLaw::MapPlaneStrainVectorTo3D(
    const Vector& rPlaneStrainVector)
{
    KRATOS_ERROR_IF(rPlaneStrainVector.size() != GeoPlaneStrainSize)
        << "Expected a plane-strain vector of size " << GeoPlaneStrainSize
        << ", got " << rPlaneStrainVector.size() << "." << std::endl;

    auto result = Vector{ZeroVector(Reference3DSize)};
    for (IndexType i = 0; i < GeoPlaneStrainSize; ++i) {
        result[i] = rPlaneStrainVector[i];
    }
    return result;
}

Vector Reference3DMohrCoulombPlaneStrainLaw::Map3DVectorToPlaneStrain(
    const Vector& rThreeDimensionalVector)
{
    KRATOS_ERROR_IF(rThreeDimensionalVector.size() != Reference3DSize)
        << "Expected a 3D vector of size " << Reference3DSize
        << ", got " << rThreeDimensionalVector.size() << "." << std::endl;

    auto result = Vector(GeoPlaneStrainSize);
    for (IndexType i = 0; i < GeoPlaneStrainSize; ++i) {
        result[i] = rThreeDimensionalVector[i];
    }
    return result;
}

Matrix Reference3DMohrCoulombPlaneStrainLaw::Map3DMatrixToPlaneStrain(
    const Matrix& rThreeDimensionalMatrix)
{
    KRATOS_ERROR_IF(rThreeDimensionalMatrix.size1() != Reference3DSize ||
                    rThreeDimensionalMatrix.size2() != Reference3DSize)
        << "Expected a 3D matrix of size " << Reference3DSize << "x"
        << Reference3DSize << ", got " << rThreeDimensionalMatrix.size1() << "x"
        << rThreeDimensionalMatrix.size2() << "." << std::endl;

    auto result = Matrix(GeoPlaneStrainSize, GeoPlaneStrainSize);
    for (IndexType i = 0; i < GeoPlaneStrainSize; ++i) {
        for (IndexType j = 0; j < GeoPlaneStrainSize; ++j) {
            result(i, j) = rThreeDimensionalMatrix(i, j);
        }
    }
    return result;
}

std::string Reference3DMohrCoulombPlaneStrainLaw::Info() const
{
    return "Reference3DMohrCoulombPlaneStrainLaw";
}

} // namespace Kratos
