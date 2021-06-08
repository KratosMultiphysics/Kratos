// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

// System includes

// External includes

// Project includes
#include "custom_elements/geo_truss_element_3D2N.hpp"
#include "includes/define.h"
#include "geo_mechanics_application_variables.h"
// This will be used as soon as all structural elements are derived
// #include "../StructuralMechanicsApplication/custom_utilities/structural_mechanics_element_utilities.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos {
GeoTrussElement3D2N::
    GeoTrussElement3D2N(IndexType NewId,
                        GeometryType::Pointer pGeometry)
    : TrussElement3D2N(NewId, pGeometry) {}

GeoTrussElement3D2N::
    GeoTrussElement3D2N(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties)
    : TrussElement3D2N(NewId, pGeometry, pProperties) {}

Element::Pointer
    GeoTrussElement3D2N::Create(IndexType NewId, NodesArrayType const& rThisNodes,
                         PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = GetGeometry();
    return Kratos::make_intrusive<GeoTrussElement3D2N>(NewId, rGeom.Create(rThisNodes), pProperties);
}

Element::Pointer
    GeoTrussElement3D2N::Create(IndexType NewId, GeometryType::Pointer pGeom,
                         PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GeoTrussElement3D2N>(NewId, pGeom, pProperties);
}

GeoTrussElement3D2N::~GeoTrussElement3D2N() {}


void GeoTrussElement3D2N::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    TrussElement3D2N::Initialize(rCurrentProcessInfo);

    mIsInitialization = true;

    KRATOS_CATCH("")
}

void GeoTrussElement3D2N::
    CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                 std::vector<Vector>& rOutput,
                                 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints();
    if (rOutput.size() != integration_points.size()) {
        rOutput.resize(integration_points.size());
    }
    if (rVariable == GREEN_LAGRANGE_STRAIN_VECTOR) {
        Vector strain = ZeroVector(msDimension);
        strain[0] = CalculateGreenLagrangeStrain();
        strain[1] = 0.00;
        strain[2] = 0.00;
        rOutput[0] = strain;
    }
    if (rVariable == PK2_STRESS_VECTOR) {

        array_1d<double, msDimension> temp_internal_stresses = ZeroVector(msDimension);
        ProcessInfo temp_process_information;

        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),temp_process_information);
        Vector temp_strain = ZeroVector(1);
        temp_strain[0] = CalculateGreenLagrangeStrain();
        Values.SetStrainVector(temp_strain);
        mpConstitutiveLaw->CalculateValue(Values,FORCE,temp_internal_stresses);

        temp_internal_stresses += mInternalStressesFinalizedPrevious;

        rOutput[0] = temp_internal_stresses;
    }
    if (rVariable == CAUCHY_STRESS_VECTOR) {

        array_1d<double, msDimension> temp_internal_stresses = ZeroVector(msDimension);
        ProcessInfo temp_process_information;

        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),temp_process_information);
        Vector temp_strain = ZeroVector(1);
        temp_strain[0] = CalculateGreenLagrangeStrain();
        Values.SetStrainVector(temp_strain);
        mpConstitutiveLaw->CalculateValue(Values,FORCE,temp_internal_stresses);

        temp_internal_stresses += mInternalStressesFinalizedPrevious;

        const double l = GeoStructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
        const double L0 = GeoStructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);

        rOutput[0] = temp_internal_stresses*l/L0;
    }

    KRATOS_CATCH("")
}

void GeoTrussElement3D2N::
    CalculateOnIntegrationPoints(const Variable<array_1d<double, 3>>& rVariable,
                                 std::vector<array_1d<double, 3>>& rOutput,
                                 const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints();
    if (rOutput.size() != integration_points.size()) {
        rOutput.resize(integration_points.size());
    }

    if (rVariable == FORCE) {
        BoundedVector<double, msDimension> truss_forces = ZeroVector(msDimension);
        truss_forces[2] = 0.00;
        truss_forces[1] = 0.00;
        const double A = GetProperties()[CROSS_AREA];

        double prestress = 0.00;
        if (GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
            prestress = GetProperties()[TRUSS_PRESTRESS_PK2];
        }

        const double L0 = GeoStructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
        const double l = GeoStructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);


        array_1d<double, msDimension> temp_internal_stresses = ZeroVector(msDimension);
        ProcessInfo temp_process_information;
        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),temp_process_information);

        Vector temp_strain = ZeroVector(1);
        temp_strain[0] = CalculateGreenLagrangeStrain();
        Values.SetStrainVector(temp_strain);
        mpConstitutiveLaw->CalculateValue(Values,FORCE,temp_internal_stresses);

        temp_internal_stresses += mInternalStressesFinalizedPrevious;

        truss_forces[0] = ((temp_internal_stresses[0] + prestress) * l * A) / L0;

        rOutput[0] = truss_forces;
    }

    KRATOS_CATCH("")

}


void GeoTrussElement3D2N::UpdateInternalForces(
    BoundedVector<double, TrussElement3D2N::msLocalSize>& rInternalForces,
    const ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY
    BoundedMatrix<double, msLocalSize, msLocalSize> transformation_matrix =
        ZeroMatrix(msLocalSize, msLocalSize);

    CreateTransformationMatrix(transformation_matrix);

    const double l = GeoStructuralMechanicsElementUtilities::CalculateCurrentLength3D2N(*this);
    const double L0 = GeoStructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double A = GetProperties()[CROSS_AREA];

    double prestress = 0.00;
    if (GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
        prestress = GetProperties()[TRUSS_PRESTRESS_PK2];
    }

    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);
    Vector temp_strain = ZeroVector(1);
    Vector temp_stress = ZeroVector(1);
    temp_strain[0] = CalculateGreenLagrangeStrain();
    Values.SetStrainVector(temp_strain);
    Values.SetStressVector(temp_stress);
    mpConstitutiveLaw->CalculateMaterialResponse(Values,ConstitutiveLaw::StressMeasure_PK2);

    mInternalStresses = temp_stress;

    temp_stress += mInternalStressesFinalizedPrevious;

    const double normal_force =
        ((temp_stress[0] + prestress) * l * A) / L0;

    // internal force vectors
    BoundedVector<double, msLocalSize> f_local = ZeroVector(msLocalSize);
    f_local[0] = -1.00 * normal_force;
    f_local[3] = 1.00 * normal_force;
    rInternalForces = ZeroVector(msLocalSize);
    noalias(rInternalForces) = prod(transformation_matrix, f_local);
    KRATOS_CATCH("");
}


void GeoTrussElement3D2N::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    TrussElement3D2N::InitializeSolutionStep(rCurrentProcessInfo);

    if (mIsInitialization)
    {
        if (rCurrentProcessInfo.Has(RESET_DISPLACEMENTS))
        {
            bool ResetDisplacement = rCurrentProcessInfo[RESET_DISPLACEMENTS];
            if (ResetDisplacement)
            {
                mInternalStressesFinalizedPrevious = mInternalStressesFinalized;
            }
            else
            {
                mInternalStressesFinalized = mInternalStressesFinalizedPrevious;
            }
        }
    }
    mIsInitialization = false;

    KRATOS_CATCH("")

}

void GeoTrussElement3D2N::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    TrussElement3D2N::FinalizeSolutionStep(rCurrentProcessInfo);
    mInternalStressesFinalized = mInternalStresses + mInternalStressesFinalizedPrevious;

    KRATOS_CATCH("");
}

void GeoTrussElement3D2N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TrussElement3D2N);
    rSerializer.save("mpConstitutiveLaw", mpConstitutiveLaw);
    rSerializer.save("InternalStresses", mInternalStresses);
    rSerializer.save("InternalStressesFinalized", mInternalStressesFinalized);
    rSerializer.save("InternalStressesFinalizedPrevious", mInternalStressesFinalizedPrevious);
}

void GeoTrussElement3D2N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TrussElement3D2N);
    rSerializer.load("mpConstitutiveLaw", mpConstitutiveLaw);
    rSerializer.load("InternalStresses", mInternalStresses);
    rSerializer.load("InternalStressesFinalized", mInternalStressesFinalized);
    rSerializer.load("InternalStressesFinalizedPrevious", mInternalStressesFinalizedPrevious);
}

} // namespace Kratos.
