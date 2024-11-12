// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Klaus B. Sautter
//

// System includes

// External includes

// Project includes
#include "custom_elements/truss_element_linear_3D2N.hpp"
#include "includes/define.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/structural_mechanics_element_utilities.h"

namespace Kratos {
TrussElementLinear3D2N::TrussElementLinear3D2N(IndexType NewId,
        GeometryType::Pointer pGeometry)
    : TrussElement3D2N(NewId, pGeometry) {}

TrussElementLinear3D2N::TrussElementLinear3D2N(
    IndexType NewId, GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties)
    : TrussElement3D2N(NewId, pGeometry, pProperties) {}

Element::Pointer
TrussElementLinear3D2N::Create(IndexType NewId,
                               NodesArrayType const& rThisNodes,
                               PropertiesType::Pointer pProperties) const
{
    const GeometryType& rGeom = GetGeometry();
    return Kratos::make_intrusive<TrussElementLinear3D2N>(
               NewId, rGeom.Create(rThisNodes), pProperties);
}

Element::Pointer
TrussElementLinear3D2N::Create(IndexType NewId,
                               GeometryType::Pointer pGeom,
                               PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<TrussElementLinear3D2N>(
               NewId, pGeom, pProperties);
}

TrussElementLinear3D2N::~TrussElementLinear3D2N() {}

BoundedMatrix<double, TrussElement3D2N::msLocalSize,
TrussElement3D2N::msLocalSize>
TrussElementLinear3D2N::CreateElementStiffnessMatrix(
    const ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY
    BoundedMatrix<double, msLocalSize, msLocalSize> LocalStiffnessMatrix =
        ZeroMatrix(msLocalSize, msLocalSize);
    CalculateElasticStiffnessMatrix(LocalStiffnessMatrix,
                                    rCurrentProcessInfo);

    return LocalStiffnessMatrix;
    KRATOS_CATCH("")
}

void TrussElementLinear3D2N::AddPrestressLinear(
    VectorType& rRightHandSideVector)
{
    KRATOS_TRY
    BoundedMatrix<double, msLocalSize, msLocalSize> transformation_matrix =
        ZeroMatrix(msLocalSize, msLocalSize);
    CreateTransformationMatrix(transformation_matrix);
    double prestress = 0.00;
    if (GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
        prestress = GetProperties()[TRUSS_PRESTRESS_PK2];
    }
    const double A = GetProperties()[CROSS_AREA];
    const double N = prestress * A;

    // internal force vectors
    BoundedVector<double, msLocalSize> f_local = ZeroVector(msLocalSize);
    f_local[0] = -1.00 * N;
    f_local[3] = 1.00 * N;
    rRightHandSideVector -= prod(transformation_matrix, f_local);
    KRATOS_CATCH("")
}

void TrussElementLinear3D2N::CalculateRightHandSide(
    VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY
    rRightHandSideVector = ZeroVector(msLocalSize);

    BoundedVector<double, msLocalSize> internal_forces = ZeroVector(msLocalSize);
    UpdateInternalForces(internal_forces, rCurrentProcessInfo);

    noalias(rRightHandSideVector) -= internal_forces;

    AddPrestressLinear(rRightHandSideVector);

    // add bodyforces
    noalias(rRightHandSideVector) += CalculateBodyForces();
    KRATOS_CATCH("")
}

void TrussElementLinear3D2N::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY

    // resizing the matrices + create memory for LHS
    rLeftHandSideMatrix = ZeroMatrix(msLocalSize, msLocalSize);
    // creating LHS
    rLeftHandSideMatrix = CreateElementStiffnessMatrix(rCurrentProcessInfo);
    KRATOS_CATCH("")
}

void TrussElementLinear3D2N::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{

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

        const auto stress = CalculateStressFromLinearStrain(rCurrentProcessInfo);
        truss_forces[0] = (stress + prestress) * A;

        rOutput[0] = truss_forces;
    }
}


void TrussElementLinear3D2N::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable, std::vector<Vector>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const GeometryType::IntegrationPointsArrayType& integration_points =
        GetGeometry().IntegrationPoints();
    if (rOutput.size() != integration_points.size()) {
        rOutput.resize(integration_points.size());
    }
    if (rVariable == STRAIN) {
        Vector Strain = ZeroVector(msDimension);
        Strain[0] = CalculateLinearStrain();
        Strain[1] = 0.00;
        Strain[2] = 0.00;
        rOutput[0] = Strain;
    }
    else if (rVariable == PK2_STRESS_VECTOR) {
        auto stress = CalculateStressFromLinearStrain(rCurrentProcessInfo);
        if (GetProperties().Has(TRUSS_PRESTRESS_PK2)) {
            stress += GetProperties()[TRUSS_PRESTRESS_PK2];
        }

        rOutput[0] = ScalarVector(1, stress);
    }

    TrussElement3D2N::CalculateOnIntegrationPoints(rVariable, rOutput, rCurrentProcessInfo);
    KRATOS_CATCH("")
}




void TrussElementLinear3D2N::WriteTransformationCoordinates(
    BoundedVector<double, TrussElement3D2N::msLocalSize>
    & rReferenceCoordinates)
{
    KRATOS_TRY
    const auto& r_initial_position_node_1 =  GetGeometry()[0].GetInitialPosition();
    rReferenceCoordinates[0] = r_initial_position_node_1.X();
    rReferenceCoordinates[1] = r_initial_position_node_1.Y();
    rReferenceCoordinates[2] = r_initial_position_node_1.Z();

    const auto& r_initial_position_node_2 = GetGeometry()[1].GetInitialPosition();
    rReferenceCoordinates[3] = r_initial_position_node_2.X();
    rReferenceCoordinates[4] = r_initial_position_node_2.Y();
    rReferenceCoordinates[5] = r_initial_position_node_2.Z();
    KRATOS_CATCH("")
}

double TrussElementLinear3D2N::CalculateLinearStrain()
{
    KRATOS_TRY

    Vector current_disp = ZeroVector(msLocalSize);
    GetValuesVector(current_disp);
    BoundedMatrix<double, msLocalSize, msLocalSize> transformation_matrix =
        ZeroMatrix(msLocalSize, msLocalSize);
    CreateTransformationMatrix(transformation_matrix);

    current_disp = prod(Matrix(trans(transformation_matrix)),current_disp);
    const double length_0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double e = (current_disp[3]-current_disp[0])/length_0;

    return e;
    KRATOS_CATCH("")
}


void TrussElementLinear3D2N::UpdateInternalForces(BoundedVector<double,msLocalSize>& rInternalForces, const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto stress = CalculateStressFromLinearStrain(rCurrentProcessInfo);

    Vector temp_internal_stresses = ZeroVector(msLocalSize);
    temp_internal_stresses[0] = -1.0 * stress;
    temp_internal_stresses[3] = stress;

    rInternalForces = temp_internal_stresses*GetProperties()[CROSS_AREA];


    BoundedMatrix<double, msLocalSize, msLocalSize> transformation_matrix =
        ZeroMatrix(msLocalSize, msLocalSize);
    CreateTransformationMatrix(transformation_matrix);

    rInternalForces = prod(transformation_matrix, rInternalForces);

    KRATOS_CATCH("")
}

void TrussElementLinear3D2N::FinalizeSolutionStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);
    Vector temp_strain = ZeroVector(1);
    Vector temp_stress = ZeroVector(1);
    temp_strain[0] = CalculateLinearStrain();
    Values.SetStrainVector(temp_strain);
    Values.SetStressVector(temp_stress);
    mpConstitutiveLaw->FinalizeMaterialResponse(Values,ConstitutiveLaw::StressMeasure_PK2);
    KRATOS_CATCH("")
}

double TrussElementLinear3D2N::ReturnTangentModulus1D(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    return ReturnTangentModulus1D(CalculateLinearStrain(), rCurrentProcessInfo);
    KRATOS_CATCH("")
}

double TrussElementLinear3D2N::CalculateStressFromLinearStrain(const ProcessInfo &rCurrentProcessInfo)
{
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),rCurrentProcessInfo);

    Vector temp_strain = ScalarVector(1, CalculateLinearStrain());
    Vector temp_stress = ZeroVector(1);

    Values.SetStrainVector(temp_strain);
    Values.SetStressVector(temp_stress);
    mpConstitutiveLaw->CalculateMaterialResponse(Values,ConstitutiveLaw::StressMeasure_PK2);

    return temp_stress[0];
}

void TrussElementLinear3D2N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TrussElement3D2N)
    rSerializer.save("mConstitutiveLaw", mpConstitutiveLaw);
}
void TrussElementLinear3D2N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TrussElement3D2N)
    rSerializer.load("mConstitutiveLaw", mpConstitutiveLaw);
}

} // namespace Kratos.
