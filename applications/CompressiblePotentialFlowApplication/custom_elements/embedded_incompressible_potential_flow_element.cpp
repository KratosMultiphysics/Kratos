//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Núñez, based on Iñigo Lopez and Riccardo Rossi work
//
#include "embedded_incompressible_potential_flow_element.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_utilities/potential_flow_utilities.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <int Dim, int NumNodes>
Element::Pointer EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, typename PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<EmbeddedIncompressiblePotentialFlowElement>(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::Create(
    IndexType NewId, typename GeometryType::Pointer pGeom, typename PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<EmbeddedIncompressiblePotentialFlowElement>(
        NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::Clone(
    IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<EmbeddedIncompressiblePotentialFlowElement>(
        NewId, this->GetGeometry().Create(ThisNodes), this->pGetProperties());
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    const EmbeddedIncompressiblePotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);
    const int kutta = r_this.GetValue(KUTTA);

    BoundedVector<double,NumNodes> distances;
    for(unsigned int i_node = 0; i_node<NumNodes; i_node++){
        distances[i_node] = this->GetGeometry()[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);
    }
    const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<Dim,NumNodes>(distances);

    if (is_embedded && wake == 0 && kutta == 0) {
        CalculateEmbeddedLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
        if (std::abs(rCurrentProcessInfo[PENALTY_COEFFICIENT]) > std::numeric_limits<double>::epsilon()) {
            AddPotentialGradientStabilizationTerm(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
        }
    }
    else {
        BaseType::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
    }

}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::CalculateEmbeddedLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
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

    potential = PotentialFlowUtilities::GetPotentialOnNormalElement<Dim, NumNodes>(*this);

    ModifiedShapeFunctions::Pointer pModifiedShFunc = this->pGetModifiedShapeFunctions(distances);
    Matrix positive_side_sh_func;
    ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients;
    Vector positive_side_weights;
    pModifiedShFunc -> ComputePositiveSideShapeFunctionsAndGradientsValues(
        positive_side_sh_func,
        positive_side_sh_func_gradients,
        positive_side_weights,
        GeometryData::GI_GAUSS_1);

    const double free_stream_density = rCurrentProcessInfo[FREE_STREAM_DENSITY];

    BoundedMatrix<double,NumNodes,Dim> DN_DX;
    for (unsigned int i_gauss=0;i_gauss<positive_side_sh_func_gradients.size();i_gauss++){
        DN_DX=positive_side_sh_func_gradients(i_gauss);
        noalias(rLeftHandSideMatrix) += free_stream_density*prod(DN_DX,trans(DN_DX))*positive_side_weights(i_gauss);
    }

    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, potential);
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::AddPotentialGradientStabilizationTerm(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    array_1d<double, NumNodes> potential;
    potential = PotentialFlowUtilities::GetPotentialOnNormalElement<Dim, NumNodes>(*this);

    std::vector<array_1d<double, Dim>> nodal_gradient_vector(NumNodes);
    for(std::size_t i_node=0; i_node<NumNodes; ++i_node) {
        auto& nodal_gradient = nodal_gradient_vector[i_node];
        nodal_gradient.clear();
        double neighbour_elements_total_area = 0.0;
        auto neighbour_elem = this->GetGeometry()[i_node].GetValue(NEIGHBOUR_ELEMENTS);
		KRATOS_ERROR_IF(neighbour_elem.size() == 0) << this->Info() << " neighbour elements were not computed\n";

        for (const auto r_elem : neighbour_elem){

            BoundedVector<double,NumNodes> neighbour_distances;
            for(unsigned int i = 0; i<NumNodes; i++){
                neighbour_distances[i] = r_elem.GetGeometry()[i].GetSolutionStepValue(GEOMETRY_DISTANCE);
            }
            const bool is_neighbour_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<Dim,NumNodes>(neighbour_distances);
            if(!is_neighbour_embedded && r_elem.Is(ACTIVE)) {
                auto r_geometry = r_elem.GetGeometry();
                const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod();
                const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method);
                Vector detJ0;
                PotentialFlowUtilities::ElementalData<NumNodes,Dim> neighbour_data;

                GeometryUtils::CalculateGeometryData(r_geometry, neighbour_data.DN_DX, neighbour_data.N, neighbour_data.vol);
                neighbour_data.potentials = PotentialFlowUtilities::GetPotentialOnNormalElement<Dim, NumNodes>(r_elem);
                r_geometry.DeterminantOfJacobian(detJ0, r_integration_method);

                const int is_neighbour_wake = r_elem.GetValue(WAKE);
                Vector neighbour_elemental_gradient;
                if (is_neighbour_wake == 0) {
                    neighbour_elemental_gradient = PotentialFlowUtilities::ComputeVelocityNormalElement<Dim,NumNodes>(r_elem);
                }
                else {
                    neighbour_elemental_gradient = PotentialFlowUtilities::ComputeVelocityUpperWakeElement<Dim,NumNodes>(r_elem);
                }

                for (IndexType i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss){
                    const double gauss_point_volume = r_integration_points[i_gauss].Weight() * detJ0[i_gauss];
                    IndexType neighbour_node_id = -1;
                    for(std::size_t j=0; j<NumNodes; ++j) {
                        if (this->GetGeometry()[i_node].Id() == r_elem.GetGeometry()[j].Id()){
                            neighbour_node_id = j;
                            break;
                        }
                    }

                    KRATOS_ERROR_IF(neighbour_node_id<0)<<"No neighbour node was found for neighbour element " << r_elem.Id() << " and element " << this-> Id() <<std::endl;

                    for(std::size_t k=0; k<Dim; ++k) {
                        nodal_gradient[k] += neighbour_data.N[neighbour_node_id] * gauss_point_volume * neighbour_elemental_gradient[k];
                    }
                    neighbour_elements_total_area += neighbour_data.N[neighbour_node_id] * gauss_point_volume;
                }
            }
        }
        nodal_gradient = nodal_gradient/neighbour_elements_total_area;
    }

    array_1d<double,Dim> averaged_nodal_gradient;
    averaged_nodal_gradient.clear();
    int number_of_positive_nodes = 0;

    for (IndexType i_node=0; i_node<NumNodes; i_node++){
        if (this->GetGeometry()[i_node].FastGetSolutionStepValue(GEOMETRY_DISTANCE)>0.0){
            number_of_positive_nodes += 1;
            averaged_nodal_gradient += nodal_gradient_vector[i_node];
        }
    }
    averaged_nodal_gradient = averaged_nodal_gradient/number_of_positive_nodes;

    PotentialFlowUtilities::ElementalData<NumNodes,Dim> data;
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), data.DN_DX, data.N, data.vol);

    auto penalty_term_nodal_gradient = data.vol*prod(data.DN_DX, averaged_nodal_gradient);
    auto penalty_term_potential = data.vol*prod(data.DN_DX,trans(data.DN_DX));
    auto penalty_coefficient = rCurrentProcessInfo[PENALTY_COEFFICIENT];

    noalias(rLeftHandSideMatrix) +=  penalty_coefficient*penalty_term_potential;
    noalias(rRightHandSideVector) += penalty_coefficient*(penalty_term_nodal_gradient-prod(penalty_term_potential, potential));
}

template <>
ModifiedShapeFunctions::Pointer EmbeddedIncompressiblePotentialFlowElement<2,3>::pGetModifiedShapeFunctions(Vector& rDistances) {
    return Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(this->pGetGeometry(), rDistances);
}

template <>
ModifiedShapeFunctions::Pointer EmbeddedIncompressiblePotentialFlowElement<3,4>::pGetModifiedShapeFunctions(Vector& rDistances) {
    return Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(this->pGetGeometry(), rDistances);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template <int Dim, int NumNodes>
int EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::Check(const ProcessInfo& rCurrentProcessInfo)
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
std::string EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "EmbeddedIncompressiblePotentialFlowElement #" << this->Id();
    return buffer.str();
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "EmbeddedIncompressiblePotentialFlowElement #" << this->Id();
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::PrintData(std::ostream& rOStream) const
{
    this->pGetGeometry()->PrintData(rOStream);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class EmbeddedIncompressiblePotentialFlowElement<2, 3>;
template class EmbeddedIncompressiblePotentialFlowElement<3, 4>;


} // namespace Kratos
