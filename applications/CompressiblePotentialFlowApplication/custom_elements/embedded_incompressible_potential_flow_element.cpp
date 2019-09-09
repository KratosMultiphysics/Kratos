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

    if (is_embedded && wake == 0 && kutta == 0)
        CalculateEmbeddedLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
    else
        BaseType::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);

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

    potential = PotentialFlowUtilities::GetPotentialOnNormalElement<2,3>(*this);


    ModifiedShapeFunctions::Pointer pModifiedShFunc = this->pGetModifiedShapeFunctions(distances);
    Matrix positive_side_sh_func;
    ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients;
    Vector positive_side_weights;
    pModifiedShFunc -> ComputePositiveSideShapeFunctionsAndGradientsValues(
        positive_side_sh_func,
        positive_side_sh_func_gradients,
        positive_side_weights,
        GeometryData::GI_GAUSS_1);

    BoundedMatrix<double,NumNodes,Dim> DN_DX;
    for (unsigned int i_gauss=0;i_gauss<positive_side_sh_func_gradients.size();i_gauss++){
        DN_DX=positive_side_sh_func_gradients(i_gauss);
        noalias(rLeftHandSideMatrix) += prod(DN_DX,trans(DN_DX))*positive_side_weights(i_gauss);;
    }

    std::vector<array_1d<double, Dim>> xi_vector(NumNodes);
    for(std::size_t i_node=0; i_node<NumNodes; ++i_node) {
        auto& xi = xi_vector[i_node]    ;
        xi.clear();
        double this_area = 0.0;
        auto neighbour_elem = this->GetGeometry()[i_node].GetValue(NEIGHBOUR_ELEMENTS);

        for (auto r_elem : neighbour_elem){

            BoundedVector<double,NumNodes> distances;
            for(unsigned int i = 0; i<NumNodes; i++){
                distances[i] = r_elem.GetGeometry()[i].GetSolutionStepValue(GEOMETRY_DISTANCE);
            }
            const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<Dim,NumNodes>(distances);
            if(!is_embedded) {
                auto r_geometry = r_elem.GetGeometry();
                const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod();
                const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method);
                Vector detJ0;
                GeometryData::ShapeFunctionsGradientsType DN_DX;

                PotentialFlowUtilities::ElementalData<NumNodes,Dim> this_data;
                GeometryUtils::CalculateGeometryData(r_geometry, this_data.DN_DX, this_data.N, this_data.vol);
                this_data.potentials = PotentialFlowUtilities::GetPotentialOnNormalElement<2,3>(r_elem);
                r_geometry.ShapeFunctionsIntegrationPointsGradients(DN_DX, detJ0, r_integration_method);


                const Vector elemental_gradient = prod(trans(this_data.DN_DX), this_data.potentials);

                for (IndexType i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss){
                    const double gauss_point_volume = r_integration_points[i_gauss].Weight() * detJ0[i_gauss];
                    IndexType this_node = -1;
                    for(std::size_t j=0; j<NumNodes; ++j) {
                        if (this->GetGeometry()[i_node].Id() == r_elem.GetGeometry()[j].Id()){
                            this_node = j;
                            break;
                        }
                    }
                    for(std::size_t k=0; k<Dim; ++k) {
                        xi[k] += this_data.N[this_node] * gauss_point_volume * elemental_gradient[k];
                    }
                    this_area += this_data.N[this_node] * gauss_point_volume;
                }
            }
        }
        xi = xi/this_area;
        xi[2] = 0.0;
        this->GetGeometry()[i_node].GetValue(VELOCITY_LOWER) = xi;
    }

    array_1d<double,Dim> averaged_xi;
    averaged_xi.clear();
    int number_of_positive_nodes = 0;

    for (IndexType i_node=0; i_node<NumNodes; i_node++){
        if (this->GetGeometry()[i_node].FastGetSolutionStepValue(GEOMETRY_DISTANCE)>0.0){
            number_of_positive_nodes += 1;
            averaged_xi += xi_vector[i_node];
        }
    }
    averaged_xi = averaged_xi/number_of_positive_nodes;

    PotentialFlowUtilities::ElementalData<NumNodes,Dim> data;
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), data.DN_DX, data.N, data.vol);
    auto penaltyTerm = data.vol*prod(DN_DX, averaged_xi);
    const auto penalty = rCurrentProcessInfo[INITIAL_PENALTY];

    noalias(rLeftHandSideMatrix) = (1+penalty) * rLeftHandSideMatrix;
    noalias(rRightHandSideVector) = penalty*penaltyTerm-prod(rLeftHandSideMatrix, potential);

}

template <>
ModifiedShapeFunctions::Pointer EmbeddedIncompressiblePotentialFlowElement<2,3>::pGetModifiedShapeFunctions(Vector& rDistances) {
    return Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(this->pGetGeometry(), rDistances);
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

} // namespace Kratos
