// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:     BSD License
//           license: structural_mechanics_application/license.txt
//
//  Main authors: Klaus B. Sautter
//
//
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
    ProcessInfo& rCurrentProcessInfo)
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
    KRATOS_TRY;
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
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY
    rRightHandSideVector = ZeroVector(msLocalSize);

    BoundedVector<double, msLocalSize> internal_forces = ZeroVector(msLocalSize);
    UpdateInternalForces(internal_forces);

    noalias(rRightHandSideVector) -= internal_forces;

    AddPrestressLinear(rRightHandSideVector);

    // add bodyforces
    noalias(rRightHandSideVector) += CalculateBodyForces();
    KRATOS_CATCH("")
}

void TrussElementLinear3D2N::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
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
        ProcessInfo temp_process_information;

        ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),temp_process_information);
        Vector temp_strain = ZeroVector(1);
        Vector temp_stress = ZeroVector(1);
        temp_strain[0] = CalculateLinearStrain();
        Values.SetStrainVector(temp_strain);
        Values.SetStressVector(temp_stress);
        mpConstitutiveLaw->CalculateMaterialResponse(Values,ConstitutiveLaw::StressMeasure_PK2);



        truss_forces[0] = (temp_stress[0] + prestress) * A;

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
    KRATOS_CATCH("")
}




void TrussElementLinear3D2N::WriteTransformationCoordinates(
    BoundedVector<double, TrussElement3D2N::msLocalSize>
    & rReferenceCoordinates)
{
    KRATOS_TRY;
    rReferenceCoordinates = ZeroVector(msLocalSize);
    rReferenceCoordinates[0] = GetGeometry()[0].X0();
    rReferenceCoordinates[1] = GetGeometry()[0].Y0();
    rReferenceCoordinates[2] = GetGeometry()[0].Z0();
    rReferenceCoordinates[3] = GetGeometry()[1].X0();
    rReferenceCoordinates[4] = GetGeometry()[1].Y0();
    rReferenceCoordinates[5] = GetGeometry()[1].Z0();
    KRATOS_CATCH("");
}

double TrussElementLinear3D2N::CalculateLinearStrain()
{
    KRATOS_TRY;

    Vector current_disp = ZeroVector(msLocalSize);
    GetValuesVector(current_disp);
    BoundedMatrix<double, msLocalSize, msLocalSize> transformation_matrix =
        ZeroMatrix(msLocalSize, msLocalSize);
    CreateTransformationMatrix(transformation_matrix);

    current_disp = prod(Matrix(trans(transformation_matrix)),current_disp);
    const double length_0 = StructuralMechanicsElementUtilities::CalculateReferenceLength3D2N(*this);
    const double e = (current_disp[3]-current_disp[0])/length_0;

    return e;
    KRATOS_CATCH("");
}


void TrussElementLinear3D2N::UpdateInternalForces(BoundedVector<double,msLocalSize>& rInternalForces)
{
    KRATOS_TRY;

    Vector temp_internal_stresses = ZeroVector(msLocalSize);
    ProcessInfo temp_process_information;
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),temp_process_information);

    Vector temp_strain = ZeroVector(1);
    Vector temp_stress = ZeroVector(1);
    temp_strain[0] = CalculateLinearStrain();
    Values.SetStrainVector(temp_strain);
    Values.SetStressVector(temp_stress);
    mpConstitutiveLaw->CalculateMaterialResponse(Values,ConstitutiveLaw::StressMeasure_PK2);

    temp_internal_stresses[0] = -1.0*temp_stress[0];
    temp_internal_stresses[3] = temp_stress[0];

    rInternalForces = temp_internal_stresses*GetProperties()[CROSS_AREA];


    BoundedMatrix<double, msLocalSize, msLocalSize> transformation_matrix =
        ZeroMatrix(msLocalSize, msLocalSize);
    CreateTransformationMatrix(transformation_matrix);

    rInternalForces = prod(transformation_matrix, rInternalForces);

    KRATOS_CATCH("");
}

void TrussElementLinear3D2N::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;
    ProcessInfo temp_process_information;
    ConstitutiveLaw::Parameters Values(GetGeometry(),GetProperties(),temp_process_information);
    Vector temp_strain = ZeroVector(1);
    Vector temp_stress = ZeroVector(1);
    temp_strain[0] = CalculateLinearStrain();
    Values.SetStrainVector(temp_strain);
    Values.SetStressVector(temp_stress);
    mpConstitutiveLaw->FinalizeMaterialResponse(Values,ConstitutiveLaw::StressMeasure_PK2);
    KRATOS_CATCH("");
}


void TrussElementLinear3D2N::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, TrussElement3D2N);
    rSerializer.save("mConstitutiveLaw", mpConstitutiveLaw);
}
void TrussElementLinear3D2N::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, TrussElement3D2N);
    rSerializer.load("mConstitutiveLaw", mpConstitutiveLaw);
}

} // namespace Kratos.
