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
#include "custom_elements/geo_truss_element_linear_3D2N.hpp"
#include "includes/define.h"
#include "geo_mechanics_application_variables.h"
#include "../StructuralMechanicsApplication/custom_utilities/structural_mechanics_element_utilities.h"
#include "includes/checks.h"

namespace Kratos {
GeoTrussElementLinear3D2N::
    GeoTrussElementLinear3D2N(IndexType NewId,
                              GeometryType::Pointer pGeometry)
    : TrussElementLinear3D2N(NewId, pGeometry) {}

GeoTrussElementLinear3D2N::
    GeoTrussElementLinear3D2N(IndexType NewId,
                              GeometryType::Pointer pGeometry,
                              PropertiesType::Pointer pProperties)
    : TrussElementLinear3D2N(NewId, pGeometry, pProperties) {}

Element::Pointer
    GeoTrussElementLinear3D2N::
        Create(IndexType NewId,
               NodesArrayType const& rThisNodes,
               PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = GetGeometry();
    return Kratos::make_intrusive<GeoTrussElementLinear3D2N>(NewId, rGeom.Create(rThisNodes), pProperties);
}

Element::Pointer
    GeoTrussElementLinear3D2N::
        Create(IndexType NewId,
               GeometryType::Pointer pGeom,
               PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GeoTrussElementLinear3D2N>(NewId, pGeom, pProperties);
}

GeoTrussElementLinear3D2N::~GeoTrussElementLinear3D2N() {}


void GeoTrussElementLinear3D2N::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    TrussElementLinear3D2N::Initialize(rCurrentProcessInfo);

    mIsInitialization = true;

    KRATOS_CATCH("")
}

void GeoTrussElementLinear3D2N::
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

        array_1d<double, msDimension> temp_internal_stresses = ZeroVector(msDimension);

        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);
        Vector temp_strain = ZeroVector(1);
        Vector temp_stress = ZeroVector(1);
        temp_strain[0] = CalculateLinearStrain();
        Values.SetStrainVector(temp_strain);
        Values.SetStressVector(temp_stress);
        mpConstitutiveLaw->CalculateMaterialResponse(Values,ConstitutiveLaw::StressMeasure_PK2);

        temp_stress += mInternalStressesFinalizedPrevious;

        truss_forces[0] = (temp_stress[0] + prestress) * A;

        rOutput[0] = truss_forces;
    }

    KRATOS_CATCH("")

}


void GeoTrussElementLinear3D2N::
    UpdateInternalForces(BoundedVector<double,msLocalSize>& rInternalForces,
                         const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    Vector temp_internal_stresses = ZeroVector(msLocalSize);
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    Vector temp_strain = ZeroVector(1);
    Vector temp_stress = ZeroVector(1);
    temp_strain[0] = CalculateLinearStrain();
    Values.SetStrainVector(temp_strain);
    Values.SetStressVector(temp_stress);
    mpConstitutiveLaw->CalculateMaterialResponse(Values,ConstitutiveLaw::StressMeasure_PK2);

    mInternalStresses = temp_stress;

    temp_stress += mInternalStressesFinalizedPrevious;

    temp_internal_stresses[0] = -1.0*temp_stress[0];
    temp_internal_stresses[3] = temp_stress[0];

    rInternalForces = temp_internal_stresses*GetProperties()[CROSS_AREA];


    BoundedMatrix<double, msLocalSize, msLocalSize> transformation_matrix =
        ZeroMatrix(msLocalSize, msLocalSize);
    CreateTransformationMatrix(transformation_matrix);

    rInternalForces = prod(transformation_matrix, rInternalForces);

    KRATOS_CATCH("");
}


void GeoTrussElementLinear3D2N::InitializeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    TrussElementLinear3D2N::InitializeSolutionStep(rCurrentProcessInfo);

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

void GeoTrussElementLinear3D2N::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    TrussElementLinear3D2N::FinalizeSolutionStep(rCurrentProcessInfo);
    mInternalStressesFinalized = mInternalStresses + mInternalStressesFinalizedPrevious;

    KRATOS_CATCH("");
}

void GeoTrussElementLinear3D2N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TrussElementLinear3D2N);
    rSerializer.save("mpConstitutiveLaw", mpConstitutiveLaw);
    rSerializer.save("InternalStresses", mInternalStresses);
    rSerializer.save("InternalStressesFinalized", mInternalStressesFinalized);
    rSerializer.save("InternalStressesFinalizedPrevious", mInternalStressesFinalizedPrevious);
}

void GeoTrussElementLinear3D2N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TrussElementLinear3D2N);
    rSerializer.load("mpConstitutiveLaw", mpConstitutiveLaw);
    rSerializer.load("InternalStresses", mInternalStresses);
    rSerializer.load("InternalStressesFinalized", mInternalStressesFinalized);
    rSerializer.load("InternalStressesFinalizedPrevious", mInternalStressesFinalizedPrevious);
}

} // namespace Kratos.
