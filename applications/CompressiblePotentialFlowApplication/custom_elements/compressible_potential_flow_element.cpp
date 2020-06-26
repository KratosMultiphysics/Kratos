//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Inigo Lopez and Riccardo Rossi
//

#include "compressible_potential_flow_element.h"
#include "compressible_potential_flow_application_variables.h"
#include "fluid_dynamics_application_variables.h"
#include "includes/cfd_variables.h"
#include "fluid_dynamics_application_variables.h"
#include "custom_utilities/potential_flow_utilities.h"

namespace Kratos
{
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <int Dim, int NumNodes>
Element::Pointer CompressiblePotentialFlowElement<Dim, NumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<CompressiblePotentialFlowElement>(
        NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer CompressiblePotentialFlowElement<Dim, NumNodes>::Create(
    IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<CompressiblePotentialFlowElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer CompressiblePotentialFlowElement<Dim, NumNodes>::Clone(
    IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<CompressiblePotentialFlowElement>(
        NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    const CompressiblePotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake == 0) // Normal element (non-wake) - eventually an embedded
        CalculateLocalSystemNormalElement(
            rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
    else // Wake element
        CalculateLocalSystemWakeElement(
            rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    // TODO: improve speed
    Matrix tmp;
    CalculateLocalSystem(tmp, rRightHandSideVector, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    // TODO: improve speed
    VectorType tmp;
    CalculateLocalSystem(rLeftHandSideMatrix, tmp, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::EquationIdVector(
    EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    const CompressiblePotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake == 0) // Normal element
    {
        if (rResult.size() != NumNodes)
            rResult.resize(NumNodes, false);

        const int kutta = r_this.GetValue(KUTTA);

        if (kutta == 0)
            GetEquationIdVectorNormalElement(rResult);
        else
            GetEquationIdVectorKuttaElement(rResult);
    }
    else // Wake element
    {
        if (rResult.size() != 2 * NumNodes)
            rResult.resize(2 * NumNodes, false);

        GetEquationIdVectorWakeElement(rResult);
    }
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::GetDofList(DofsVectorType& rElementalDofList,
                                                                 ProcessInfo& CurrentProcessInfo)
{
    const CompressiblePotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake == 0) // Normal element
    {
        if (rElementalDofList.size() != NumNodes)
            rElementalDofList.resize(NumNodes);

        const int kutta = r_this.GetValue(KUTTA);

        if (kutta == 0)
            GetDofListNormalElement(rElementalDofList);
        else
            GetDofListKuttaElement(rElementalDofList);
    }
    else // wake element
    {
        if (rElementalDofList.size() != 2 * NumNodes)
            rElementalDofList.resize(2 * NumNodes);

        GetDofListWakeElement(rElementalDofList);
    }
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    bool active = true;
    if ((this)->IsDefined(ACTIVE))
        active = (this)->Is(ACTIVE);

    const CompressiblePotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    if (wake != 0 && active == true)
    {
        ComputePotentialJump(rCurrentProcessInfo);
    }
    ComputeElementInternalEnergy();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template <int Dim, int NumNodes>
int CompressiblePotentialFlowElement<Dim, NumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Generic geometry check
    int out = Element::Check(rCurrentProcessInfo);
    if (out != 0)
    {
        return out;
    }

    KRATOS_ERROR_IF(GetGeometry().Area() <= 0.0)
        << this->Id() << "Area cannot be less than or equal to 0" << std::endl;

    for (unsigned int i = 0; i < this->GetGeometry().size(); i++)
    {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY_POTENTIAL, this->GetGeometry()[i]);
    }

    return out;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::GetValueOnIntegrationPoints(
    const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);

    if (rVariable == PRESSURE_COEFFICIENT)
    {
        rValues[0] = PotentialFlowUtilities::ComputeCompressiblePressureCoefficient<Dim, NumNodes>(*this, rCurrentProcessInfo);
    }
    else if (rVariable == DENSITY)
    {
        rValues[0] = ComputeDensity(rCurrentProcessInfo);
    }
    else if (rVariable == MACH)
    {
        rValues[0] = PotentialFlowUtilities::ComputeLocalMachNumber<Dim, NumNodes>(*this, rCurrentProcessInfo);
    }
    else if (rVariable == SOUND_VELOCITY)
    {
        rValues[0] = PotentialFlowUtilities::ComputeLocalSpeedOfSound<Dim, NumNodes>(*this, rCurrentProcessInfo);
    }
    else if (rVariable == WAKE)
    {
        const CompressiblePotentialFlowElement& r_this = *this;
        rValues[0] = r_this.GetValue(WAKE);
    }
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::GetValueOnIntegrationPoints(
    const Variable<int>& rVariable, std::vector<int>& rValues, const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);
    if (rVariable == TRAILING_EDGE)
        rValues[0] = this->GetValue(TRAILING_EDGE);
    else if (rVariable == KUTTA)
        rValues[0] = this->GetValue(KUTTA);
    else if (rVariable == WAKE)
        rValues[0] = this->GetValue(WAKE);
    else if (rVariable == ZERO_VELOCITY_CONDITION)
        rValues[0] = this->GetValue(ZERO_VELOCITY_CONDITION);
    else if (rVariable == TRAILING_EDGE_ELEMENT)
        rValues[0] = this->GetValue(TRAILING_EDGE_ELEMENT);
    else if (rVariable == DECOUPLED_TRAILING_EDGE_ELEMENT)
        rValues[0] = this->GetValue(DECOUPLED_TRAILING_EDGE_ELEMENT);
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rValues.size() != 1)
        rValues.resize(1);
    if (rVariable == VELOCITY)
    {
        array_1d<double, 3> v(3, 0.0);
        array_1d<double, Dim> vaux = PotentialFlowUtilities::ComputeVelocity<Dim, NumNodes>(*this);
        for (unsigned int k = 0; k < Dim; k++)
            v[k] = vaux[k];
        rValues[0] = v;
    }
    else if (rVariable == PERTURBATION_VELOCITY)
    {
        const array_1d<double, Dim>& free_stream_velocity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
        array_1d<double, 3> v(3, 0.0);
        array_1d<double, Dim> vaux = PotentialFlowUtilities::ComputeVelocity<Dim, NumNodes>(*this);
        for (unsigned int k = 0; k < Dim; k++)
            v[k] = vaux[k] - free_stream_velocity[k];
        rValues[0] = v;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <int Dim, int NumNodes>
std::string CompressiblePotentialFlowElement<Dim, NumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "CompressiblePotentialFlowElement #" << Id();
    return buffer.str();
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "CompressiblePotentialFlowElement #" << Id();
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::PrintData(std::ostream& rOStream) const
{
    pGetGeometry()->PrintData(rOStream);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::GetWakeDistances(
    array_1d<double, NumNodes>& distances) const
{
    noalias(distances) = GetValue(WAKE_ELEMENTAL_DISTANCES);
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::GetEquationIdVectorNormalElement(
    EquationIdVectorType& rResult) const
{
    for (unsigned int i = 0; i < NumNodes; i++)
        rResult[i] = GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::GetEquationIdVectorKuttaElement(
    EquationIdVectorType& rResult) const
{
    // Kutta elements have only negative part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (!GetGeometry()[i].GetValue(TRAILING_EDGE))
            rResult[i] = GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
        else
            rResult[i] = GetGeometry()[i].GetDof(AUXILIARY_VELOCITY_POTENTIAL).EquationId();
    }
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::GetEquationIdVectorWakeElement(
    EquationIdVectorType& rResult) const
{
    array_1d<double, NumNodes> distances;
    GetWakeDistances(distances);

    // Positive part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (distances[i] > 0.0)
            rResult[i] = GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
        else
            rResult[i] =
                GetGeometry()[i].GetDof(AUXILIARY_VELOCITY_POTENTIAL, 0).EquationId();
    }

    // Negative part - sign is opposite to the previous case
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (distances[i] < 0.0)
            rResult[NumNodes + i] =
                GetGeometry()[i].GetDof(VELOCITY_POTENTIAL).EquationId();
        else
            rResult[NumNodes + i] =
                GetGeometry()[i].GetDof(AUXILIARY_VELOCITY_POTENTIAL).EquationId();
    }
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::GetDofListNormalElement(DofsVectorType& rElementalDofList) const
{
    for (unsigned int i = 0; i < NumNodes; i++)
        rElementalDofList[i] = GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::GetDofListKuttaElement(DofsVectorType& rElementalDofList) const
{
    // Kutta elements have only negative part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (!GetGeometry()[i].GetValue(TRAILING_EDGE))
            rElementalDofList[i] = GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
        else
            rElementalDofList[i] = GetGeometry()[i].pGetDof(AUXILIARY_VELOCITY_POTENTIAL);
    }
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::GetDofListWakeElement(DofsVectorType& rElementalDofList) const
{
    array_1d<double, NumNodes> distances;
    GetWakeDistances(distances);

    // Positive part
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (distances[i] > 0)
            rElementalDofList[i] = GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
        else
            rElementalDofList[i] = GetGeometry()[i].pGetDof(AUXILIARY_VELOCITY_POTENTIAL);
    }

    // Negative part - sign is opposite to the previous case
    for (unsigned int i = 0; i < NumNodes; i++)
    {
        if (distances[i] < 0)
            rElementalDofList[NumNodes + i] = GetGeometry()[i].pGetDof(VELOCITY_POTENTIAL);
        else
            rElementalDofList[NumNodes + i] =
                GetGeometry()[i].pGetDof(AUXILIARY_VELOCITY_POTENTIAL);
    }
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::CalculateLocalSystemNormalElement(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != NumNodes || rLeftHandSideMatrix.size2() != NumNodes)
        rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);
    if (rRightHandSideVector.size() != NumNodes)
        rRightHandSideVector.resize(NumNodes, false);
    rLeftHandSideMatrix.clear();

    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    const double density = ComputeDensity(rCurrentProcessInfo);
    const double DrhoDu2 = ComputeDensityDerivative(density, rCurrentProcessInfo);

    // Computing local velocity
    array_1d<double, Dim> v =  PotentialFlowUtilities::ComputeVelocity<Dim, NumNodes> (*this);

    const BoundedVector<double, NumNodes> DNV = prod(data.DN_DX, v);

    noalias(rLeftHandSideMatrix) +=
        data.vol * density * prod(data.DN_DX, trans(data.DN_DX));
    noalias(rLeftHandSideMatrix) += data.vol * 2 * DrhoDu2 * outer_prod(DNV, trans(DNV));

    const BoundedMatrix<double, NumNodes, NumNodes> rLaplacianMatrix =
        data.vol * density * prod(data.DN_DX, trans(data.DN_DX));

    data.potentials= PotentialFlowUtilities::GetPotentialOnNormalElement<Dim, NumNodes>(*this);
    noalias(rRightHandSideVector) = -prod(rLaplacianMatrix, data.potentials);
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::CalculateLocalSystemWakeElement(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    // Note that the lhs and rhs have double the size
    if (rLeftHandSideMatrix.size1() != 2 * NumNodes ||
        rLeftHandSideMatrix.size2() != 2 * NumNodes)
        rLeftHandSideMatrix.resize(2 * NumNodes, 2 * NumNodes, false);
    if (rRightHandSideVector.size() != 2 * NumNodes)
        rRightHandSideVector.resize(2 * NumNodes, false);
    rLeftHandSideMatrix.clear();
    rRightHandSideVector.clear();

    MatrixType rLaplacianMatrix = ZeroMatrix(2 * NumNodes, 2 * NumNodes);

    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);
    GetWakeDistances(data.distances);

    const double density = ComputeDensity(rCurrentProcessInfo);
    const double DrhoDu2 = ComputeDensityDerivative(density, rCurrentProcessInfo);

    // Computing local velocity
    array_1d<double, Dim> v = PotentialFlowUtilities::ComputeVelocityUpperWakeElement<Dim, NumNodes>(*this);

    const BoundedVector<double, NumNodes> DNV = prod(data.DN_DX, v);

    const BoundedMatrix<double, NumNodes, NumNodes> laplacian_total =
        data.vol * density * prod(data.DN_DX, trans(data.DN_DX));

    const BoundedMatrix<double, NumNodes, NumNodes> lhs_total =
        data.vol * density * prod(data.DN_DX, trans(data.DN_DX)) +
        data.vol * 2 * DrhoDu2 * outer_prod(DNV, trans(DNV));

    if (this->Is(STRUCTURE))
    {
        Matrix lhs_positive = ZeroMatrix(NumNodes, NumNodes);
        Matrix lhs_negative = ZeroMatrix(NumNodes, NumNodes);

        Matrix laplacian_positive = ZeroMatrix(NumNodes, NumNodes);
        Matrix laplacian_negative = ZeroMatrix(NumNodes, NumNodes);

        CalculateLocalSystemSubdividedElement(lhs_positive, lhs_negative, laplacian_positive,
                                              laplacian_negative, rCurrentProcessInfo);
        AssignLocalSystemSubdividedElement(
            rLeftHandSideMatrix, lhs_positive, lhs_negative, lhs_total, rLaplacianMatrix,
            laplacian_positive, laplacian_negative, laplacian_total, data);
    }
    else
    {
        AssignLocalSystemWakeElement(rLeftHandSideMatrix, lhs_total, data);
        AssignLocalSystemWakeElement(rLaplacianMatrix, laplacian_total, data);
    }

    Vector split_element_values(2 * NumNodes);
    split_element_values = PotentialFlowUtilities::GetPotentialOnWakeElement<Dim, NumNodes>(*this, data.distances);
    noalias(rRightHandSideVector) = -prod(rLaplacianMatrix, split_element_values);
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::CalculateLocalSystemSubdividedElement(
    Matrix& lhs_positive,
    Matrix& lhs_negative,
    Matrix& laplacian_positive,
    Matrix& laplacian_negative,
    const ProcessInfo& rCurrentProcessInfo)
{
    ElementalData<NumNodes, Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(GetGeometry(), data.DN_DX, data.N, data.vol);

    GetWakeDistances(data.distances);

    // Subdivide the element
    constexpr unsigned int nvolumes = 3 * (Dim - 1);
    BoundedMatrix<double, NumNodes, Dim> Points;
    array_1d<double, nvolumes> PartitionsSign;
    BoundedMatrix<double, nvolumes, NumNodes> GPShapeFunctionValues;
    array_1d<double, nvolumes> Volumes;
    std::vector<Matrix> GradientsValue(nvolumes);
    BoundedMatrix<double, nvolumes, 2> NEnriched;
    for (unsigned int i = 0; i < GradientsValue.size(); ++i)
        GradientsValue[i].resize(2, Dim, false);
    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        const array_1d<double, 3>& coords = GetGeometry()[i].Coordinates();
        for (unsigned int k = 0; k < Dim; ++k)
        {
            Points(i, k) = coords[k];
        }
    }

    const unsigned int nsubdivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(
        Points, data.DN_DX, data.distances, Volumes, GPShapeFunctionValues,
        PartitionsSign, GradientsValue, NEnriched);

    const double density = ComputeDensity(rCurrentProcessInfo);
    const double DrhoDu2 = ComputeDensityDerivative(density, rCurrentProcessInfo);

    // Computing local velocity
    array_1d<double, Dim> v = PotentialFlowUtilities::ComputeVelocityUpperWakeElement<Dim, NumNodes>(*this);

    const BoundedVector<double, NumNodes> DNV = prod(data.DN_DX, v);

    // Compute the lhs and rhs that would correspond to it being divided
    for (unsigned int i = 0; i < nsubdivisions; ++i)
    {
        if (PartitionsSign[i] > 0)
        {
            noalias(lhs_positive) +=
                Volumes[i] * density * prod(data.DN_DX, trans(data.DN_DX));
            noalias(lhs_positive) +=
                Volumes[i] * 2 * DrhoDu2 * outer_prod(DNV, trans(DNV));

            noalias(laplacian_positive) +=
                Volumes[i] * density * prod(data.DN_DX, trans(data.DN_DX));
        }
        else
        {
            noalias(lhs_negative) +=
                Volumes[i] * density * prod(data.DN_DX, trans(data.DN_DX));
            noalias(lhs_negative) +=
                Volumes[i] * 2 * DrhoDu2 * outer_prod(DNV, trans(DNV));

            noalias(laplacian_negative) +=
                Volumes[i] * density * prod(data.DN_DX, trans(data.DN_DX));
        }
    }
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::ComputeLHSGaussPointContribution(
    const double weight, Matrix& lhs, const ElementalData<NumNodes, Dim>& data) const
{
    noalias(lhs) += weight * prod(data.DN_DX, trans(data.DN_DX));
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::AssignLocalSystemSubdividedElement(
    Matrix& rLeftHandSideMatrix,
    Matrix& lhs_positive,
    Matrix& lhs_negative,
    const BoundedMatrix<double, NumNodes, NumNodes>& lhs_total,
    MatrixType& rLaplacianMatrix,
    Matrix& laplacian_positive,
    Matrix& laplacian_negative,
    const BoundedMatrix<double, NumNodes, NumNodes>& laplacian_total,
    const ElementalData<NumNodes, Dim>& data) const
{
    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        // The TE node takes the contribution of the subdivided element and
        // we do not apply the wake condition on the TE node
        if (GetGeometry()[i].GetValue(TRAILING_EDGE))
        {
            for (unsigned int j = 0; j < NumNodes; ++j)
            {
                rLeftHandSideMatrix(i, j) = lhs_positive(i, j);
                rLeftHandSideMatrix(i + NumNodes, j + NumNodes) = lhs_negative(i, j);

                rLaplacianMatrix(i, j) = laplacian_positive(i, j);
                rLaplacianMatrix(i + NumNodes, j + NumNodes) = laplacian_negative(i, j);
            }
        }
        else
        {
            AssignLocalSystemWakeNode(rLeftHandSideMatrix, lhs_total, data, i);
            AssignLocalSystemWakeNode(rLaplacianMatrix, laplacian_total, data, i);
        }
    }
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::AssignLocalSystemWakeElement(
    MatrixType& rLeftHandSideMatrix,
    const BoundedMatrix<double, NumNodes, NumNodes>& lhs_total,
    const ElementalData<NumNodes, Dim>& data) const
{
    for (unsigned int row = 0; row < NumNodes; ++row)
        AssignLocalSystemWakeNode(rLeftHandSideMatrix, lhs_total, data, row);
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::AssignLocalSystemWakeNode(
    MatrixType& rLeftHandSideMatrix,
    const BoundedMatrix<double, NumNodes, NumNodes>& lhs_total,
    const ElementalData<NumNodes, Dim>& data,
    unsigned int& row) const
{
    // Filling the diagonal blocks (i.e. decoupling upper and lower dofs)
    for (unsigned int column = 0; column < NumNodes; ++column)
    {
        rLeftHandSideMatrix(row, column) = lhs_total(row, column);
        rLeftHandSideMatrix(row + NumNodes, column + NumNodes) = lhs_total(row, column);
    }

    // Applying wake condition on the AUXILIARY_VELOCITY_POTENTIAL dofs
    if (data.distances[row] < 0.0)
        for (unsigned int column = 0; column < NumNodes; ++column)
            rLeftHandSideMatrix(row, column + NumNodes) = -lhs_total(row, column); // Side 1
    else if (data.distances[row] > 0.0)
        for (unsigned int column = 0; column < NumNodes; ++column)
            rLeftHandSideMatrix(row + NumNodes, column) = -lhs_total(row, column); // Side 2
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::ComputePotentialJump(const ProcessInfo& rCurrentProcessInfo)
{
    const array_1d<double, 3>& vinfinity = rCurrentProcessInfo[FREE_STREAM_VELOCITY];
    const double vinfinity_norm = sqrt(inner_prod(vinfinity, vinfinity));

    array_1d<double, NumNodes> distances;
    GetWakeDistances(distances);

    for (unsigned int i = 0; i < NumNodes; i++)
    {
        double aux_potential =
            GetGeometry()[i].FastGetSolutionStepValue(AUXILIARY_VELOCITY_POTENTIAL);
        double potential = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_POTENTIAL);
        double potential_jump = aux_potential - potential;

        if (distances[i] > 0)
        {
            GetGeometry()[i].SetValue(POTENTIAL_JUMP,
                                      -2.0 / vinfinity_norm * (potential_jump));
        }
        else
        {
            GetGeometry()[i].SetValue(POTENTIAL_JUMP, 2.0 / vinfinity_norm * (potential_jump));
        }
    }
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::ComputeElementInternalEnergy()
{
    array_1d<double, Dim> velocity = PotentialFlowUtilities::ComputeVelocity<Dim, NumNodes>(*this);

    double internal_energy = 0.5 * inner_prod(velocity, velocity);
    this->SetValue(INTERNAL_ENERGY, std::abs(internal_energy));
}

template <int Dim, int NumNodes>
double CompressiblePotentialFlowElement<Dim, NumNodes>::ComputeDensity(const ProcessInfo& rCurrentProcessInfo) const
{
    // Reading free stream conditions
    const double rho_inf = rCurrentProcessInfo[FREE_STREAM_DENSITY];
    const double M_inf = rCurrentProcessInfo[FREE_STREAM_MACH];
    const double heat_capacity_ratio = rCurrentProcessInfo[HEAT_CAPACITY_RATIO];
    const double mach_number_limit = rCurrentProcessInfo[MACH_LIMIT];

    // Computing local mach number
    double local_mach_number = PotentialFlowUtilities::ComputeLocalMachNumber<Dim, NumNodes>(*this, rCurrentProcessInfo);

    if (local_mach_number > mach_number_limit)
    { // Clamping the mach number to mach_number_limit
        KRATOS_WARNING("ComputeDensity") << "Clamping the local mach number to " << mach_number_limit << std::endl;
        local_mach_number = mach_number_limit;
    }

    // Computing squares
    const double M_inf_2 = M_inf * M_inf;
    const double M_2 = local_mach_number * local_mach_number;

    // Computing density according to Equation 8.9 of Drela, M. (2014) Flight Vehicle
    // Aerodynamics, The MIT Press, London
    const double numerator = 1 + (heat_capacity_ratio - 1) * M_inf_2 / 2;
    const double denominator = 1 + (heat_capacity_ratio - 1) * M_2 / 2;
    const double base = numerator / denominator;

    if (base > 0.0)
    {
        return rho_inf * pow(base, 1 / (heat_capacity_ratio - 1));
    }
    else
    {
        KRATOS_WARNING("ComputeDensity") << "Using density correction" << std::endl;
        return rho_inf * 0.00001;
    }
}

template <int Dim, int NumNodes>
double CompressiblePotentialFlowElement<Dim, NumNodes>::ComputeDensityDerivative(
    const double rho, const ProcessInfo& rCurrentProcessInfo) const
{
    // Reading free stream conditions
    const double rho_inf = rCurrentProcessInfo[FREE_STREAM_DENSITY];
    const double heat_capacity_ratio = rCurrentProcessInfo[HEAT_CAPACITY_RATIO];
    const double a_inf = rCurrentProcessInfo[SOUND_VELOCITY];

    return -pow(rho_inf, heat_capacity_ratio - 1) *
           pow(rho, 2 - heat_capacity_ratio) / (2 * a_inf * a_inf);
}

// serializer

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}

template <int Dim, int NumNodes>
void CompressiblePotentialFlowElement<Dim, NumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class CompressiblePotentialFlowElement<2, 3>;
template class CompressiblePotentialFlowElement<3, 4>;

} // namespace Kratos
