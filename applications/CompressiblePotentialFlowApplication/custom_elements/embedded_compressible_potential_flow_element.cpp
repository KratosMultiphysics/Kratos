//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Nunez, based on Inigo Lopez and Riccardo Rossi work
//
#include "embedded_compressible_potential_flow_element.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_utilities/potential_flow_utilities.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <int Dim, int NumNodes>
Element::Pointer EmbeddedCompressiblePotentialFlowElement<Dim, NumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, typename PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<EmbeddedCompressiblePotentialFlowElement>(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer EmbeddedCompressiblePotentialFlowElement<Dim, NumNodes>::Create(
    IndexType NewId, typename GeometryType::Pointer pGeom, typename PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<EmbeddedCompressiblePotentialFlowElement>(
        NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer EmbeddedCompressiblePotentialFlowElement<Dim, NumNodes>::Clone(
    IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<EmbeddedCompressiblePotentialFlowElement>(
        NewId, this->GetGeometry().Create(ThisNodes), this->pGetProperties());
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
void EmbeddedCompressiblePotentialFlowElement<Dim, NumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    const EmbeddedCompressiblePotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    BoundedVector<double,NumNodes> distances;
    for(unsigned int i_node = 0; i_node<NumNodes; i_node++){
        distances[i_node] = this->GetGeometry()[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);
    }
    const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<Dim,NumNodes>(distances);

    if (is_embedded && wake == 0) {
        CalculateEmbeddedLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
        if (std::abs(rCurrentProcessInfo[STABILIZATION_FACTOR]) > std::numeric_limits<double>::epsilon()) {
            PotentialFlowUtilities::AddPotentialGradientStabilizationTerm<Dim, NumNodes>(*this,rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
        }
    }
    else {
        if (this->Is(STRUCTURE)) {
            CalculateKuttaWakeLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
        } else {
            BaseType::CalculateLeftHandSide(rLeftHandSideMatrix, rCurrentProcessInfo);
            BaseType::CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
        }
    }

    if (std::abs(rCurrentProcessInfo[PENALTY_COEFFICIENT]) > std::numeric_limits<double>::epsilon()) {
        PotentialFlowUtilities::AddKuttaConditionPenaltyTerm<Dim,NumNodes>(r_this,rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
    }

}

template <int Dim, int NumNodes>
void EmbeddedCompressiblePotentialFlowElement<Dim, NumNodes>::CalculateRightHandSide(
    VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType tmp;
    CalculateLocalSystem(tmp, rRightHandSideVector, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void EmbeddedCompressiblePotentialFlowElement<Dim, NumNodes>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    VectorType tmp;
    CalculateLocalSystem(rLeftHandSideMatrix, tmp, rCurrentProcessInfo);
}

template <int Dim, int NumNodes>
void EmbeddedCompressiblePotentialFlowElement<Dim, NumNodes>::CalculateEmbeddedLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != NumNodes || rLeftHandSideMatrix.size2() != NumNodes)
        rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);
    if (rRightHandSideVector.size() != NumNodes)
        rRightHandSideVector.resize(NumNodes, false);
    rLeftHandSideMatrix.clear();

    array_1d<double, NumNodes> potential;
    Vector distances(NumNodes);
    for(unsigned int i_node = 0; i_node<NumNodes; i_node++)
        distances(i_node) = this->GetGeometry()[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);

    const double density = BaseType::ComputeDensity(rCurrentProcessInfo);
    const double DrhoDu2 = BaseType::ComputeDensityDerivative(density, rCurrentProcessInfo);

    potential = PotentialFlowUtilities::GetPotentialOnNormalElement<Dim, NumNodes>(*this);

    ModifiedShapeFunctions::Pointer pModifiedShFunc = this->pGetModifiedShapeFunctions(distances);
    Matrix positive_side_sh_func;
    ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients;
    Vector positive_side_weights;
    pModifiedShFunc -> ComputePositiveSideShapeFunctionsAndGradientsValues(
        positive_side_sh_func,
        positive_side_sh_func_gradients,
        positive_side_weights,
        GeometryData::GI_GAUSS_2);

    // Computing local velocity
    const array_1d<double, Dim>& local_velocity = PotentialFlowUtilities::ComputeVelocityNormalElement<Dim,NumNodes>(*this);
    const double local_velocity_squared = inner_prod(local_velocity, local_velocity);
    const double max_velocity_squared = PotentialFlowUtilities::ComputeMaximumVelocitySquared<Dim, NumNodes>(rCurrentProcessInfo);

    BoundedMatrix<double,NumNodes,Dim> DN_DX;
    BoundedVector<double, NumNodes> DNV;
    BoundedMatrix<double, NumNodes, NumNodes> laplacian;
    laplacian.clear();

    for (unsigned int i_gauss=0;i_gauss<positive_side_sh_func_gradients.size();i_gauss++){
        DN_DX = positive_side_sh_func_gradients(i_gauss);
        DNV = prod(DN_DX, local_velocity);
        BoundedMatrix<double, NumNodes, NumNodes> aux_matrix = positive_side_weights(i_gauss) * density * prod(DN_DX, trans(DN_DX));
        noalias(laplacian) += aux_matrix;
        noalias(rLeftHandSideMatrix) += aux_matrix;
        if (local_velocity_squared < max_velocity_squared){
            noalias(rLeftHandSideMatrix) += positive_side_weights(i_gauss) * 2 * DrhoDu2 * outer_prod(DNV, trans(DNV));
        }
    }

    noalias(rRightHandSideVector) = -prod(laplacian, potential);
}

template <int Dim, int NumNodes>
void EmbeddedCompressiblePotentialFlowElement<Dim, NumNodes>::CalculateKuttaWakeLocalSystem(
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

    MatrixType laplacian_matrix = ZeroMatrix(2 * NumNodes, 2 * NumNodes);

    PotentialFlowUtilities::ElementalData<NumNodes,Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), data.DN_DX, data.N, data.vol);
    data.distances = PotentialFlowUtilities::GetWakeDistances<Dim, NumNodes>(*this);

    const double density = BaseType::ComputeDensity(rCurrentProcessInfo);
    const double DrhoDu2 = BaseType::ComputeDensityDerivative(density, rCurrentProcessInfo);

    // Computing local velocity
    array_1d<double, Dim> v_upper = PotentialFlowUtilities::ComputeVelocityUpperWakeElement<Dim, NumNodes>(*this);
    array_1d<double, Dim> v_lower = PotentialFlowUtilities::ComputeVelocityLowerWakeElement<Dim, NumNodes>(*this);

    const BoundedVector<double, NumNodes> DNV_upper = prod(data.DN_DX, v_upper);
    const BoundedVector<double, NumNodes> DNV_lower = prod(data.DN_DX, v_lower);

    const BoundedMatrix<double, NumNodes, NumNodes> laplacian_total =
        data.vol * density * prod(data.DN_DX, trans(data.DN_DX));

    const BoundedMatrix<double, NumNodes, NumNodes> lhs_total_upper =
        data.vol * density * prod(data.DN_DX, trans(data.DN_DX)) +
        data.vol * 2 * DrhoDu2 * outer_prod(DNV_upper, trans(DNV_upper));

    const BoundedMatrix<double, NumNodes, NumNodes> lhs_total_lower =
        data.vol * density * prod(data.DN_DX, trans(data.DN_DX)) +
        data.vol * 2 * DrhoDu2 * outer_prod(DNV_lower, trans(DNV_lower));

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        for (unsigned int j = 0; j < NumNodes; ++j)
        {
            rLeftHandSideMatrix(i, j) = lhs_total_upper(i, j);
            rLeftHandSideMatrix(i + NumNodes, j + NumNodes) = lhs_total_lower(i, j);

            laplacian_matrix(i, j) = laplacian_total(i, j);
            laplacian_matrix(i + NumNodes, j + NumNodes) = laplacian_total(i, j);
        }
    }

    BoundedVector<double, 2*NumNodes> split_element_values;
    split_element_values = PotentialFlowUtilities::GetPotentialOnWakeElement<Dim, NumNodes>(*this, data.distances);
    noalias(rRightHandSideVector) = -prod(laplacian_matrix, split_element_values);
}

template <>
ModifiedShapeFunctions::Pointer EmbeddedCompressiblePotentialFlowElement<2,3>::pGetModifiedShapeFunctions(Vector& rDistances) {
    return Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(this->pGetGeometry(), rDistances);
}

template <>
ModifiedShapeFunctions::Pointer EmbeddedCompressiblePotentialFlowElement<3,4>::pGetModifiedShapeFunctions(Vector& rDistances) {
    return Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(this->pGetGeometry(), rDistances);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template <int Dim, int NumNodes>
int EmbeddedCompressiblePotentialFlowElement<Dim, NumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Generic geometry check
    int out = BaseType::Check(rCurrentProcessInfo);
    if (out != 0)
    {
        return out;
    }

    for (unsigned int i = 0; i < this->GetGeometry().size(); i++)
    {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(GEOMETRY_DISTANCE,this->GetGeometry()[i]);
    }

    return out;

    KRATOS_CATCH("");
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <int Dim, int NumNodes>
std::string EmbeddedCompressiblePotentialFlowElement<Dim, NumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "EmbeddedCompressiblePotentialFlowElement #" << this->Id();
    return buffer.str();
}

template <int Dim, int NumNodes>
void EmbeddedCompressiblePotentialFlowElement<Dim, NumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "EmbeddedCompressiblePotentialFlowElement #" << this->Id();
}

template <int Dim, int NumNodes>
void EmbeddedCompressiblePotentialFlowElement<Dim, NumNodes>::PrintData(std::ostream& rOStream) const
{
    this->pGetGeometry()->PrintData(rOStream);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template <int Dim, int NumNodes>
void EmbeddedCompressiblePotentialFlowElement<Dim, NumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}

template <int Dim, int NumNodes>
void EmbeddedCompressiblePotentialFlowElement<Dim, NumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class EmbeddedCompressiblePotentialFlowElement<2, 3>;
template class EmbeddedCompressiblePotentialFlowElement<3, 4>;

} // namespace Kratos
