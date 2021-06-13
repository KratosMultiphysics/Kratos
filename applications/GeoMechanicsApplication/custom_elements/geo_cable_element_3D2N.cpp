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
#include "custom_elements/geo_cable_element_3D2N.hpp"
#include "includes/define.h"
#include "geo_mechanics_application_variables.h"
// This will be used as soon as all structural elements are derived
// #include "../StructuralMechanicsApplication/custom_utilities/structural_mechanics_element_utilities.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

#include "includes/checks.h"

namespace Kratos
{
GeoCableElement3D2N::
    GeoCableElement3D2N(IndexType NewId,
                        GeometryType::Pointer pGeometry)
    : GeoTrussElement3D2N(NewId, pGeometry) {}

GeoCableElement3D2N::
    GeoCableElement3D2N(IndexType NewId,
                        GeometryType::Pointer pGeometry,
                        PropertiesType::Pointer pProperties)
    : GeoTrussElement3D2N(NewId, pGeometry, pProperties) {}

Element::Pointer
    GeoCableElement3D2N::Create(IndexType NewId,
                                NodesArrayType const& rThisNodes,
                                PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = GetGeometry();
    return Kratos::make_intrusive<GeoCableElement3D2N>(NewId, rGeom.Create(rThisNodes), pProperties);
}

Element::Pointer
    GeoCableElement3D2N::Create(IndexType NewId,
                                GeometryType::Pointer pGeom,
                                PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<GeoCableElement3D2N>(NewId, pGeom, pProperties);
}

GeoCableElement3D2N::~GeoCableElement3D2N() {}

BoundedMatrix<double, TrussElement3D2N::msLocalSize, TrussElement3D2N::msLocalSize>
GeoCableElement3D2N::CreateElementStiffnessMatrix(const ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY
    BoundedMatrix<double, msLocalSize, msLocalSize> local_stiffness_matrix =
        ZeroMatrix(msLocalSize, msLocalSize);

    if (mIsCompressed) {
        local_stiffness_matrix = ZeroMatrix(msLocalSize, msLocalSize);
    }

    else {
        CalculateElasticStiffnessMatrix(local_stiffness_matrix,
                                        rCurrentProcessInfo);
        BoundedMatrix<double, msLocalSize, msLocalSize> K_geo =
            ZeroMatrix(msLocalSize, msLocalSize);
        CalculateGeometricStiffnessMatrix(K_geo, rCurrentProcessInfo);

        local_stiffness_matrix += K_geo;
    }

    return local_stiffness_matrix;
    KRATOS_CATCH("")
}

void GeoCableElement3D2N::
    CalculateRightHandSide(VectorType& rRightHandSideVector,
                           const ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY
    rRightHandSideVector = ZeroVector(msLocalSize);

    BoundedVector<double, msLocalSize> internal_forces = ZeroVector(msLocalSize);
    UpdateInternalForces(internal_forces, rCurrentProcessInfo);

    if (!mIsCompressed) {
        noalias(rRightHandSideVector) -= internal_forces;
    }
    // add bodyforces
    if (HasSelfWeight()) {
        noalias(rRightHandSideVector) += CalculateBodyForces();
    }
    KRATOS_CATCH("")
}

void GeoCableElement3D2N::
    UpdateInternalForces(BoundedVector<double, TrussElement3D2N::msLocalSize>& rInternalForces,
                         const ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY;
    const double numerical_limit = std::numeric_limits<double>::epsilon();
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

    Vector temp_internal_stresses = ZeroVector(msLocalSize);
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

    mIsCompressed = false;
    if ((normal_force < 0.00)&&(std::abs(l-L0)>numerical_limit)) {
        mIsCompressed = true;
    }

    // internal force vectors
    BoundedVector<double, msLocalSize> f_local = ZeroVector(msLocalSize);
    f_local[0] = -1.00 * normal_force;
    f_local[3] = 1.00 * normal_force;
    rInternalForces = ZeroVector(msLocalSize);
    noalias(rInternalForces) = prod(transformation_matrix, f_local);
    KRATOS_CATCH("");
}

void GeoCableElement3D2N::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == FORCE){
        GeoTrussElement3D2N::CalculateOnIntegrationPoints(rVariable,rOutput,rCurrentProcessInfo);
        if (rOutput[0][0]<0.0){
            rOutput[0]=ZeroVector(msDimension);
        }
    }
}

void GeoCableElement3D2N::
    CalculateOnIntegrationPoints(const Variable<Vector>& rVariable,
                                 std::vector<Vector>& rOutput,
                                 const ProcessInfo& rCurrentProcessInfo)
{
    if ( rVariable == GREEN_LAGRANGE_STRAIN_VECTOR ||
         rVariable == CAUCHY_STRESS_VECTOR         ||
         rVariable == PK2_STRESS_VECTOR ) {
        GeoTrussElement3D2N::CalculateOnIntegrationPoints(rVariable,rOutput,rCurrentProcessInfo);
        if (rOutput[0][0]<0.0) {
            rOutput[0]=ZeroVector(msDimension);
        }
    }
}

void GeoCableElement3D2N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, GeoTrussElement3D2N);
    rSerializer.save("mIscompressed", mIsCompressed);
}
void GeoCableElement3D2N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, GeoTrussElement3D2N);
    rSerializer.load("mIscompressed", mIsCompressed);
}
} // namespace Kratos.
