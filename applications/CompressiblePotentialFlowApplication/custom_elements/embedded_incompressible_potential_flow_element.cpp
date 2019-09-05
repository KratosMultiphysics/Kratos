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
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::InitializeNonLinearIteration(ProcessInfo& rCurrentProcessInfo){

        KRATOS_INFO_ONCE("")   << "ENTERING INITIALIZE NON LINEAR ITERATION"<< std::endl;
        for(std::size_t i_node=0; i_node<NumNodes; ++i_node) {
                array_1d<double, 3>& r_gradient = this->GetGeometry()[i_node].GetValue(VELOCITY);
                r_gradient.clear();
        }
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


    // Compute gradient
    // Calculate shape functions
    // PotentialFlowUtilities::ElementalData<NumNodes,Dim> data;

    // GeometryUtils::CalculateGeometryData(this->GetGeometry(), data.DN_DX, data.N, data.vol);

    // data.potentials = PotentialFlowUtilities::GetPotentialOnNormalElement<2,3>(*this);

    // const Vector grad = prod(trans(data.DN_DX), data.potentials);

//     auto r_geometry = this->GetGeometry();
//     const auto& r_integration_method = r_geometry.GetDefaultIntegrationMethod();
//     const auto& r_integration_points = r_geometry.IntegrationPoints(r_integration_method);
//     Vector detJ0;
//     GeometryData::ShapeFunctionsGradientsType DN_DX;

//     r_geometry.ShapeFunctionsIntegrationPointsGradients(DN_DX, detJ0, r_integration_method);
//     for (IndexType i_gauss = 0; i_gauss < r_integration_points.size(); ++i_gauss){

//         const double gauss_point_volume = r_integration_points[i_gauss].Weight() * detJ0[i_gauss];
//         for(std::size_t i_node=0; i_node<NumNodes; ++i_node) {
//             array_1d<double, 3>& r_gradient = r_geometry[i_node].GetValue(VELOCITY);
//             const double nodal_area = r_geometry[i_node].GetValue(NODAL_AREA);
//             KRATOS_ERROR_IF(nodal_area<1e-12) << "please compute nodal area" << std::endl;
//             for(std::size_t k=0; k<Dim; ++k) {
//                 // DIVIDING BY PRE-COMPUTED NODAL AREA
//                 r_gradient[k] += data.N[i_node] * gauss_point_volume * grad[k] / nodal_area;
//             }
//             // if (r_geometry[i_node].Id() == 3256) {
//             // KRATOS_WATCH(this->Id())
//             // KRATOS_WATCH(grad)
//             // KRATOS_WATCH(gauss_point_volume)
//             // KRATOS_WATCH(nodal_area)
//             // KRATOS_WATCH(r_gradient)
//             // }
//             // double& vol = r_geometry[i_node].GetValue(NODAL_AREA);
//             // vol += data.N[i_node] * gauss_point_volume;
//         }
//     }

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
    // KRATOS_WATCH(averaged_xi)
    // KRATOS_WATCH(number_of_positive_nodes)


    // array_1d<double,Dim> gauss_point_nodal_gradient;
    // gauss_point_nodal_gradient.clear();

    // for (unsigned int i_gauss=0;i_gauss<positive_side_sh_func_gradients.size();i_gauss++){
    //     for(std::size_t i_node=0; i_node<NumNodes; ++i_node) {
    //         array_1d<double, 3> nodal_gradient = this->GetGeometry()[i_node].GetValue(VELOCITY);
    //         for(std::size_t i_dim=0; i_dim<Dim; ++i_dim) {
    //             //Computing gradient at the gauss pints and averaging.
    //             gauss_point_nodal_gradient[i_dim] += positive_side_sh_func(i_gauss,i_node)*nodal_gradient[i_dim]/positive_side_weights.size(); //*positive_side_weights(i_gauss)
    //         }
    //     }
    // }

    // array_1d<double,3> elemental_gradient;
    // elemental_gradient.clear();
    // DN_DX=positive_side_sh_func_gradients(0);
    // noalias(elemental_gradient) = prod(trans(DN_DX), potential);

    PotentialFlowUtilities::ElementalData<NumNodes,Dim> data;
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), data.DN_DX, data.N, data.vol);
    auto penaltyTerm = data.vol*prod(DN_DX, averaged_xi);
    double penalty = 0.0;

    auto convergence_check = penaltyTerm-prod(rLeftHandSideMatrix, potential);
    if (this->Id() == 3144) {
        KRATOS_WATCH(prod(rLeftHandSideMatrix, potential))
        KRATOS_WATCH(penaltyTerm)
        KRATOS_WATCH(convergence_check)
    }
    noalias(rRightHandSideVector) = penalty*penaltyTerm-(1+penalty)*prod(rLeftHandSideMatrix, potential);
    rLeftHandSideMatrix = (1+penalty) * rLeftHandSideMatrix;

    // array_1d<double,Dim> gauss_point_nodal_gradient;
    // gauss_point_nodal_gradient.clear();

    // for (unsigned int i_gauss=0;i_gauss<positive_side_sh_func_gradients.size();i_gauss++){
    //     for(std::size_t i_node=0; i_node<NumNodes; ++i_node) {
    //         array_1d<double, 3> nodal_gradient = this->GetGeometry()[i_node].GetValue(VELOCITY);
    //         // if (this->GetGeometry()[i_node].Id() == 3256) {
    //         //     KRATOS_WATCH(nodal_gradient)
    //         //     KRATOS_WATCH(positive_side_sh_func_gradients)
    //         //     KRATOS_WATCH(positive_side_sh_func)
    //         // }
    //         for(std::size_t i_dim=0; i_dim<Dim; ++i_dim) {
    //             //Computing gradient at the gauss pints and averaging.
    //             gauss_point_nodal_gradient[i_dim] += positive_side_sh_func(i_gauss,i_node)*nodal_gradient[i_dim]/positive_side_weights.size(); //*positive_side_weights(i_gauss)
    //         }
    //     }
    // }

    // array_1d<double,3> elemental_gradient;
    // elemental_gradient.clear();
    // DN_DX=positive_side_sh_func_gradients(0);
    // noalias(elemental_gradient) = prod(trans(DN_DX), potential);

    // // const Vector grad =
    // array_1d<double,3> difference = gauss_point_nodal_gradient - elemental_gradient;

    // // Matrix penaltyTerm(NumNodes, NumNodes);

    // // double positive_weight = 0.0
    // // for (unsigned int i_gauss=0;i_gauss<positive_side_weight.size();i_gauss++){
    // //     positive_weight += positive_side_weight(i_gauss)
    // // }
    // PotentialFlowUtilities::ElementalData<NumNodes,Dim> data;
    // GeometryUtils::CalculateGeometryData(this->GetGeometry(), data.DN_DX, data.N, data.vol);
    // auto penaltyTerm = data.vol*prod(DN_DX, gauss_point_nodal_gradient);
    // double penalty = 0.0;

    // auto convergence_check = penaltyTerm-prod(rLeftHandSideMatrix, potential);
    // if (this->Id() == 3144) {
    //     KRATOS_WATCH(prod(rLeftHandSideMatrix, potential))
    //     KRATOS_WATCH(penaltyTerm)
    //     KRATOS_WATCH(convergence_check)
    // }
    // noalias(rRightHandSideVector) = penalty*penaltyTerm-(1+penalty)*prod(rLeftHandSideMatrix, potential);
    // rLeftHandSideMatrix = (1+penalty) * rLeftHandSideMatrix;
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
