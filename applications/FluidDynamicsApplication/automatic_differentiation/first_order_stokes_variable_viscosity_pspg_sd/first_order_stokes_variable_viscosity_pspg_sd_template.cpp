//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolas Sibuet
//

// System includes


// External includes


// Project includes
#include "utilities/element_size_calculator.h"

// Application includes
#include "first_order_stokes_variable_viscosity_pspg_sd.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< unsigned int TDim >
FirstOrderStokesVariableViscosityPspgSd<TDim>::FirstOrderStokesVariableViscosityPspgSd(IndexType NewId)
    : Element(NewId)
{}

template< unsigned int TDim >
FirstOrderStokesVariableViscosityPspgSd<TDim>::FirstOrderStokesVariableViscosityPspgSd(
    IndexType NewId,
    const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes)
{}

template< unsigned int TDim >
FirstOrderStokesVariableViscosityPspgSd<TDim>::FirstOrderStokesVariableViscosityPspgSd(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{}

template< unsigned int TDim >
FirstOrderStokesVariableViscosityPspgSd<TDim>::FirstOrderStokesVariableViscosityPspgSd(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{}

template< unsigned int TDim >
FirstOrderStokesVariableViscosityPspgSd<TDim>::~FirstOrderStokesVariableViscosityPspgSd()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< unsigned int TDim >
Element::Pointer FirstOrderStokesVariableViscosityPspgSd<TDim>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<FirstOrderStokesVariableViscosityPspgSd<TDim>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template< unsigned int TDim >
Element::Pointer FirstOrderStokesVariableViscosityPspgSd<TDim>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<FirstOrderStokesVariableViscosityPspgSd<TDim>>(NewId, pGeom, pProperties);
}

template< unsigned int TDim >
void FirstOrderStokesVariableViscosityPspgSd<TDim>::EquationIdVector(
    EquationIdVectorType &rResult,
    const ProcessInfo &rCurrentProcessInfo) const
{
    if (rResult.size() != LocalSize) {
        rResult.resize(LocalSize, false);
    }

    IndexType local_index = 0;
    const auto& r_geometry = this->GetGeometry();
    for (IndexType i = 0; i < NumNodes; ++i) {
        rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_X).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_Y).EquationId();
        if constexpr (TDim == 3) {
            rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_Z).EquationId();
        }
        rResult[local_index++] = r_geometry[i].GetDof(PRESSURE).EquationId();
    }
}

template< unsigned int TDim >
void FirstOrderStokesVariableViscosityPspgSd<TDim>::GetDofList(
    DofsVectorType &rElementalDofList,
    const ProcessInfo &rCurrentProcessInfo) const
{
    if (rElementalDofList.size() != LocalSize) {
        rElementalDofList.resize(LocalSize);
    }

    IndexType local_index = 0;
    const auto& r_geometry = this->GetGeometry();
    for (IndexType i = 0; i < NumNodes; ++i) {
        rElementalDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_X);
        rElementalDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_Y);
        if constexpr (TDim == 3) {
            rElementalDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_Z);
        }
    rElementalDofList[local_index++] = r_geometry[i].pGetDof(PRESSURE);
    }
}

template< unsigned int TDim >
void FirstOrderStokesVariableViscosityPspgSd<TDim>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Resize and initialize output
    if (rLeftHandSideMatrix.size1() != LocalSize || rLeftHandSideMatrix.size2() != LocalSize) {
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
    }

    if (rRightHandSideVector.size() != LocalSize) {
        rRightHandSideVector.resize(LocalSize, false);
    }

    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    // Initialize element data
    ElementDataContainer aux_data;
    SetElementData(rCurrentProcessInfo, aux_data);

    // Initialize constitutive law parameters
    const auto& r_geom = this->GetGeometry();
    const auto p_prop = this->GetProperties();

    // Calculate kinematics
    Vector weights;
    Matrix N;
    GeometryType::ShapeFunctionsGradientsType DN;
    CalculateKinematics(weights, N, DN);

    // Loop Gauss points
    const SizeType n_gauss = r_geom.IntegrationPointsNumber(IntegrationMethod);
    for (IndexType g = 0; g < n_gauss; ++g) {
        // Set current Gauss point kinematics
        noalias(aux_data.N) = row(N, g);
        noalias(aux_data.DN) = DN[g];
        aux_data.Weight = weights[g];

        // Assemble standard Galerkin contribution
        AddGaussPointLeftHandSideContribution(aux_data, rLeftHandSideMatrix);
        AddGaussPointRightHandSideContribution(aux_data, rRightHandSideVector);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template< unsigned int TDim >
int FirstOrderStokesVariableViscosityPspgSd<TDim>::Check(const ProcessInfo &rCurrentProcessInfo) const
{
    KRATOS_TRY;
    int out = Element::Check(rCurrentProcessInfo);
    KRATOS_ERROR_IF_NOT(out == 0)
        << "Error in base class Check for Element " << this->Info() << std::endl
        << "Error code is " << out << std::endl;

    return 0;

    KRATOS_CATCH("");
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public I/O

template< unsigned int TDim >
const Parameters FirstOrderStokesVariableViscosityPspgSd<TDim>::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["static"],
        "framework"                  : "Eulerian",
        "output"                     : {
            "gauss_point"            : [""],
            "nodal_historical"       : ["VELOCITY","PRESSURE"],
            "nodal_non_historical"   : [],
            "entity"                 : []
        },
        "required_variables"         : ["VELOCITY","PRESSURE","BODY_FORCE","DYNAMIC_VISCOSITY"]
        "required_dofs"              : [],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle2D3","Tetrahedra3D4"],
        "required_polynomial_degree_of_geometry" : 1,
        "documentation"   :
            "This implements the element from paper: Fully consistent lowest-order finite element methods for generalised Stokes flows with variable viscosity. Authors: Felipe Galarce et al."
    })");

    if (TDim == 2) {
        std::vector<std::string> dofs_2d({"VELOCITY_X","VELOCITY_Y","PRESSURE"});
        specifications["required_dofs"].SetStringArray(dofs_2d);
    } else {
        std::vector<std::string> dofs_3d({"VELOCITY_X","VELOCITY_Y","VELOCITY_Z","PRESSURE"});
        specifications["required_dofs"].SetStringArray(dofs_3d);
    }

    return specifications;
}

template< unsigned int TDim >
std::string FirstOrderStokesVariableViscosityPspgSd<TDim>::Info() const
{
    std::stringstream buffer;
    buffer << "FirstOrderStokesVariableViscosityPspgSd" << TDim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}

template< unsigned int TDim >
void FirstOrderStokesVariableViscosityPspgSd<TDim>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info() << std::endl;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private operations

template <unsigned int TDim>
void FirstOrderStokesVariableViscosityPspgSd<TDim>::SetElementData(
    const ProcessInfo& rProcessInfo,
    ElementDataContainer &rData)
{
    // Set nodal data
    const auto& r_geom = this->GetGeometry();
    for (IndexType i = 0; i < NumNodes; ++i) {
        const auto& r_v = r_geom[i].FastGetSolutionStepValue(VELOCITY);
        const auto& r_body_force = r_geom[i].FastGetSolutionStepValue(BODY_FORCE);

        for (IndexType d = 0; d < TDim; ++d) {
            rData.Velocity(i, d) = r_v[d];
            rData.BodyForce(i, d) = r_body_force[d];
        }
    }

    for (IndexType i = 0; i < NumNodes; ++i) {
        rData.Pressure[i] = r_geom[i].FastGetSolutionStepValue(PRESSURE);
        rData.DynamicViscosity[i] = r_geom[i].FastGetSolutionStepValue(DYNAMIC_VISCOSITY);
    }


    // Set stabilization values
    rData.Delta = rProcessInfo.GetValue(TAUONE); // It would be better to define a unique new variable for this parameter
    rData.Sigma = rProcessInfo.GetValue(TAUTWO); // It would be better to define a unique new variable for this parameter
}

template< unsigned int TDim >
void FirstOrderStokesVariableViscosityPspgSd<TDim>::CalculateKinematics(
    Vector& rGaussWeights,
    Matrix& rN,
    GeometryType::ShapeFunctionsGradientsType& rDNDX)
{
    // Get element geometry
    const auto& r_geom = this->GetGeometry();

    // Integration rule data
    // Note that we use the same for both velocity and pressure interpolations
    const SizeType n_gauss = r_geom.IntegrationPointsNumber(IntegrationMethod);
    const auto integration_points = r_geom.IntegrationPoints(IntegrationMethod);

    // Calculate Jacobians at integration points
    Matrix J;
    Matrix inv_J;
    double det_J;
    Vector det_J_vect(n_gauss);
    std::vector<BoundedMatrix<double, TDim, TDim>> inv_J_vect(n_gauss);
    for (IndexType g = 0; g < n_gauss; ++g) {
        r_geom.Jacobian(J, g, IntegrationMethod);
        MathUtils<double>::InvertMatrix(J, inv_J, det_J);
        det_J_vect[g] = det_J;
        noalias(inv_J_vect[g]) = inv_J;
    }

    // Calculate velocity and pressure kinematics from the geometry (P1 interpolation)
    GeometryType::UniquePointer aux_geom = nullptr;
    if constexpr (TDim == 2) {
        aux_geom = Kratos::make_unique<Triangle2D3<NodeType>>(r_geom(0), r_geom(1), r_geom(2));
    } else {
        aux_geom = Kratos::make_unique<Tetrahedra3D4<NodeType>>(r_geom(0), r_geom(1), r_geom(2), r_geom(3));
    }
    rN = aux_geom->ShapeFunctionsValues(IntegrationMethod);
    if (rDNDX.size() != n_gauss) {
        rDNDX.resize(n_gauss, false);
    }
    const auto& r_DN_De = aux_geom->ShapeFunctionsLocalGradients(IntegrationMethod);
    for (IndexType g = 0; g < n_gauss; ++g) {
        rDNDX[g] = prod(r_DN_De[g], inv_J_vect[g]);
    }

    // Calculate integration points weight
    if (rGaussWeights.size() != n_gauss) {
        rGaussWeights.resize(n_gauss, false);
    }
    for (IndexType g = 0; g < n_gauss; ++g) {
        rGaussWeights[g] = det_J_vect[g] * integration_points[g].Weight();
    }
}

template <>
void FirstOrderStokesVariableViscosityPspgSd<2>::AddGaussPointLeftHandSideContribution(
    const ElementDataContainer& rData,
    MatrixType& rLeftHandSideMatrix)
{
    //Get nodal data
    const auto& u_nodes = rData.Velocity;
    const auto& p_nodes = rData.Pressure; 
    const auto& f_nodes = rData.BodyForce;

    // Get material data
    const auto& nu_nodes = rData.DynamicViscosity;

    // Get parameters data
    const double sigma_gauss = rData.Sigma;

    // Get stabilization data
    const double delta_gauss = rData.Delta;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN;

    // Assemble LHS contribution
    const double w_g = rData.Weight;

    //substitute_element_lhs_2D3N
    
}

template <>
void FirstOrderStokesVariableViscosityPspgSd<3>::AddGaussPointLeftHandSideContribution(
    const ElementDataContainer& rData,
    MatrixType& rLeftHandSideMatrix)
{
    //Get nodal data
    const auto& u_nodes = rData.Velocity;
    const auto& p_nodes = rData.Pressure; 
    const auto& f_nodes = rData.BodyForce;

    // Get material data
    const auto& nu_nodes = rData.DynamicViscosity;

    // Get parameters data
    const double sigma_gauss = rData.Sigma;

    // Get stabilization data
    const double delta_gauss = rData.Delta;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN;

    // Assemble LHS contribution
    const double w_g = rData.Weight;

    //substitute_element_lhs_3D4N

}


template <>
void FirstOrderStokesVariableViscosityPspgSd<2>::AddGaussPointRightHandSideContribution(
    const ElementDataContainer& rData,
    VectorType& rRightHandSideVector)
{

    //Get nodal data
    const auto& u_nodes = rData.Velocity;
    const auto& p_nodes = rData.Pressure; 
    const auto& f_nodes = rData.BodyForce;

    // Get material data
    const auto& nu_nodes = rData.DynamicViscosity;

    // Get parameters data
    const double sigma_gauss = rData.Sigma;

    // Get stabilization data
    const double delta_gauss = rData.Delta;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN;

    // Assemble LHS contribution
    const double w_g = rData.Weight;

    //substitute_element_rhs_2D3N

}

template <>
void FirstOrderStokesVariableViscosityPspgSd<3>::AddGaussPointRightHandSideContribution(
    const ElementDataContainer& rData,
    VectorType& rRightHandSideVector)
{
    //Get nodal data
    const auto& u_nodes = rData.Velocity;
    const auto& p_nodes = rData.Pressure; 
    const auto& f_nodes = rData.BodyForce;

    // Get material data
    const auto& nu_nodes = rData.DynamicViscosity;

    // Get parameters data
    const double sigma_gauss = rData.Sigma;

    // Get stabilization data
    const double delta_gauss = rData.Delta;

    // Get shape function values
    const auto& N = rData.N;
    const auto& DN = rData.DN;

    // Assemble LHS contribution
    const double w_g = rData.Weight;

    //substitute_element_rhs_3D4N
    
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private serialization

template< unsigned int TDim >
void FirstOrderStokesVariableViscosityPspgSd<TDim>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}


template< unsigned int TDim >
void FirstOrderStokesVariableViscosityPspgSd<TDim>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class FirstOrderStokesVariableViscosityPspgSd<2>;
template class FirstOrderStokesVariableViscosityPspgSd<3>;

}