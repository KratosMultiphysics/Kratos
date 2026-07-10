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

        const double crLeftHandSideMatrix0 = sigma_gauss*(N[0]*N[0]);
    const double crLeftHandSideMatrix1 = DN(0,1)*DN(0,1);
    const double crLeftHandSideMatrix2 = N[0]*nu_nodes[0] + N[1]*nu_nodes[1] + N[2]*nu_nodes[2];
    const double crLeftHandSideMatrix3 = crLeftHandSideMatrix1*crLeftHandSideMatrix2;
    const double crLeftHandSideMatrix4 = DN(0,0)*DN(0,0);
    const double crLeftHandSideMatrix5 = crLeftHandSideMatrix2*crLeftHandSideMatrix4;
    const double crLeftHandSideMatrix6 = crLeftHandSideMatrix2*w_g;
    const double crLeftHandSideMatrix7 = DN(0,1)*crLeftHandSideMatrix6;
    const double crLeftHandSideMatrix8 = -DN(0,0)*crLeftHandSideMatrix7;
    const double crLeftHandSideMatrix9 = N[0]*sigma_gauss;
    const double crLeftHandSideMatrix10 = N[1]*crLeftHandSideMatrix9;
    const double crLeftHandSideMatrix11 = DN(0,1)*DN(1,1);
    const double crLeftHandSideMatrix12 = crLeftHandSideMatrix11*crLeftHandSideMatrix2;
    const double crLeftHandSideMatrix13 = DN(0,0)*DN(1,0);
    const double crLeftHandSideMatrix14 = crLeftHandSideMatrix13*crLeftHandSideMatrix2;
    const double crLeftHandSideMatrix15 = -w_g*(crLeftHandSideMatrix10 + crLeftHandSideMatrix12 + 2*crLeftHandSideMatrix14);
    const double crLeftHandSideMatrix16 = -DN(1,0)*crLeftHandSideMatrix7;
    const double crLeftHandSideMatrix17 = N[2]*crLeftHandSideMatrix9;
    const double crLeftHandSideMatrix18 = DN(0,1)*DN(2,1);
    const double crLeftHandSideMatrix19 = crLeftHandSideMatrix18*crLeftHandSideMatrix2;
    const double crLeftHandSideMatrix20 = DN(0,0)*DN(2,0);
    const double crLeftHandSideMatrix21 = crLeftHandSideMatrix2*crLeftHandSideMatrix20;
    const double crLeftHandSideMatrix22 = -w_g*(crLeftHandSideMatrix17 + crLeftHandSideMatrix19 + 2*crLeftHandSideMatrix21);
    const double crLeftHandSideMatrix23 = -DN(2,0)*crLeftHandSideMatrix7;
    const double crLeftHandSideMatrix24 = DN(0,0)*crLeftHandSideMatrix6;
    const double crLeftHandSideMatrix25 = -DN(1,1)*crLeftHandSideMatrix24;
    const double crLeftHandSideMatrix26 = -w_g*(crLeftHandSideMatrix10 + 2*crLeftHandSideMatrix12 + crLeftHandSideMatrix14);
    const double crLeftHandSideMatrix27 = -DN(2,1)*crLeftHandSideMatrix24;
    const double crLeftHandSideMatrix28 = -w_g*(crLeftHandSideMatrix17 + 2*crLeftHandSideMatrix19 + crLeftHandSideMatrix21);
    const double crLeftHandSideMatrix29 = DN(0,0)*nu_nodes[0] + DN(1,0)*nu_nodes[1] + DN(2,0)*nu_nodes[2];
    const double crLeftHandSideMatrix30 = crLeftHandSideMatrix29*delta_gauss;
    const double crLeftHandSideMatrix31 = -crLeftHandSideMatrix9;
    const double crLeftHandSideMatrix32 = DN(0,1)*nu_nodes[0] + DN(1,1)*nu_nodes[1] + DN(2,1)*nu_nodes[2];
    const double crLeftHandSideMatrix33 = DN(0,1)*crLeftHandSideMatrix32;
    const double crLeftHandSideMatrix34 = DN(0,0)*crLeftHandSideMatrix29;
    const double crLeftHandSideMatrix35 = crLeftHandSideMatrix31 + crLeftHandSideMatrix33 + 2*crLeftHandSideMatrix34;
    const double crLeftHandSideMatrix36 = DN(0,0)*delta_gauss;
    const double crLeftHandSideMatrix37 = crLeftHandSideMatrix32*delta_gauss;
    const double crLeftHandSideMatrix38 = crLeftHandSideMatrix31 + 2*crLeftHandSideMatrix33 + crLeftHandSideMatrix34;
    const double crLeftHandSideMatrix39 = DN(0,1)*delta_gauss;
    const double crLeftHandSideMatrix40 = delta_gauss*w_g;
    const double crLeftHandSideMatrix41 = crLeftHandSideMatrix11*crLeftHandSideMatrix30;
    const double crLeftHandSideMatrix42 = N[1]*sigma_gauss;
    const double crLeftHandSideMatrix43 = -crLeftHandSideMatrix42;
    const double crLeftHandSideMatrix44 = DN(1,1)*crLeftHandSideMatrix32;
    const double crLeftHandSideMatrix45 = DN(1,0)*crLeftHandSideMatrix29;
    const double crLeftHandSideMatrix46 = crLeftHandSideMatrix43 + crLeftHandSideMatrix44 + 2*crLeftHandSideMatrix45;
    const double crLeftHandSideMatrix47 = crLeftHandSideMatrix13*crLeftHandSideMatrix37;
    const double crLeftHandSideMatrix48 = crLeftHandSideMatrix43 + 2*crLeftHandSideMatrix44 + crLeftHandSideMatrix45;
    const double crLeftHandSideMatrix49 = -crLeftHandSideMatrix40*(crLeftHandSideMatrix11 + crLeftHandSideMatrix13);
    const double crLeftHandSideMatrix50 = crLeftHandSideMatrix18*crLeftHandSideMatrix30;
    const double crLeftHandSideMatrix51 = -N[2]*sigma_gauss;
    const double crLeftHandSideMatrix52 = DN(2,1)*crLeftHandSideMatrix32;
    const double crLeftHandSideMatrix53 = DN(2,0)*crLeftHandSideMatrix29;
    const double crLeftHandSideMatrix54 = crLeftHandSideMatrix51 + crLeftHandSideMatrix52 + 2*crLeftHandSideMatrix53;
    const double crLeftHandSideMatrix55 = crLeftHandSideMatrix20*crLeftHandSideMatrix37;
    const double crLeftHandSideMatrix56 = crLeftHandSideMatrix51 + 2*crLeftHandSideMatrix52 + crLeftHandSideMatrix53;
    const double crLeftHandSideMatrix57 = -crLeftHandSideMatrix40*(crLeftHandSideMatrix18 + crLeftHandSideMatrix20);
    const double crLeftHandSideMatrix58 = sigma_gauss*(N[1]*N[1]);
    const double crLeftHandSideMatrix59 = DN(1,1)*DN(1,1);
    const double crLeftHandSideMatrix60 = crLeftHandSideMatrix2*crLeftHandSideMatrix59;
    const double crLeftHandSideMatrix61 = DN(1,0)*DN(1,0);
    const double crLeftHandSideMatrix62 = crLeftHandSideMatrix2*crLeftHandSideMatrix61;
    const double crLeftHandSideMatrix63 = DN(1,1)*crLeftHandSideMatrix6;
    const double crLeftHandSideMatrix64 = -DN(1,0)*crLeftHandSideMatrix63;
    const double crLeftHandSideMatrix65 = N[2]*crLeftHandSideMatrix42;
    const double crLeftHandSideMatrix66 = DN(1,1)*DN(2,1);
    const double crLeftHandSideMatrix67 = crLeftHandSideMatrix2*crLeftHandSideMatrix66;
    const double crLeftHandSideMatrix68 = DN(1,0)*DN(2,0);
    const double crLeftHandSideMatrix69 = crLeftHandSideMatrix2*crLeftHandSideMatrix68;
    const double crLeftHandSideMatrix70 = -w_g*(crLeftHandSideMatrix65 + crLeftHandSideMatrix67 + 2*crLeftHandSideMatrix69);
    const double crLeftHandSideMatrix71 = -DN(2,0)*crLeftHandSideMatrix63;
    const double crLeftHandSideMatrix72 = DN(2,1)*crLeftHandSideMatrix6;
    const double crLeftHandSideMatrix73 = -DN(1,0)*crLeftHandSideMatrix72;
    const double crLeftHandSideMatrix74 = -w_g*(crLeftHandSideMatrix65 + 2*crLeftHandSideMatrix67 + crLeftHandSideMatrix69);
    const double crLeftHandSideMatrix75 = DN(1,0)*delta_gauss;
    const double crLeftHandSideMatrix76 = DN(1,1)*delta_gauss;
    const double crLeftHandSideMatrix77 = crLeftHandSideMatrix30*crLeftHandSideMatrix66;
    const double crLeftHandSideMatrix78 = crLeftHandSideMatrix37*crLeftHandSideMatrix68;
    const double crLeftHandSideMatrix79 = -crLeftHandSideMatrix40*(crLeftHandSideMatrix66 + crLeftHandSideMatrix68);
    const double crLeftHandSideMatrix80 = sigma_gauss*(N[2]*N[2]);
    const double crLeftHandSideMatrix81 = DN(2,1)*DN(2,1);
    const double crLeftHandSideMatrix82 = crLeftHandSideMatrix2*crLeftHandSideMatrix81;
    const double crLeftHandSideMatrix83 = DN(2,0)*DN(2,0);
    const double crLeftHandSideMatrix84 = crLeftHandSideMatrix2*crLeftHandSideMatrix83;
    const double crLeftHandSideMatrix85 = -DN(2,0)*crLeftHandSideMatrix72;
    const double crLeftHandSideMatrix86 = DN(2,0)*delta_gauss;
    const double crLeftHandSideMatrix87 = DN(2,1)*delta_gauss;
    rLeftHandSideMatrix(0,0)+=-w_g*(crLeftHandSideMatrix0 + crLeftHandSideMatrix3 + 2*crLeftHandSideMatrix5);
    rLeftHandSideMatrix(0,1)+=crLeftHandSideMatrix8;
    rLeftHandSideMatrix(0,2)+=DN(0,0)*N[0]*w_g;
    rLeftHandSideMatrix(0,3)+=crLeftHandSideMatrix15;
    rLeftHandSideMatrix(0,4)+=crLeftHandSideMatrix16;
    rLeftHandSideMatrix(0,5)+=DN(0,0)*N[1]*w_g;
    rLeftHandSideMatrix(0,6)+=crLeftHandSideMatrix22;
    rLeftHandSideMatrix(0,7)+=crLeftHandSideMatrix23;
    rLeftHandSideMatrix(0,8)+=DN(0,0)*N[2]*w_g;
    rLeftHandSideMatrix(1,0)+=crLeftHandSideMatrix8;
    rLeftHandSideMatrix(1,1)+=-w_g*(crLeftHandSideMatrix0 + 2*crLeftHandSideMatrix3 + crLeftHandSideMatrix5);
    rLeftHandSideMatrix(1,2)+=DN(0,1)*N[0]*w_g;
    rLeftHandSideMatrix(1,3)+=crLeftHandSideMatrix25;
    rLeftHandSideMatrix(1,4)+=crLeftHandSideMatrix26;
    rLeftHandSideMatrix(1,5)+=DN(0,1)*N[1]*w_g;
    rLeftHandSideMatrix(1,6)+=crLeftHandSideMatrix27;
    rLeftHandSideMatrix(1,7)+=crLeftHandSideMatrix28;
    rLeftHandSideMatrix(1,8)+=DN(0,1)*N[2]*w_g;
    rLeftHandSideMatrix(2,0)+=-w_g*(DN(0,0)*N[0] - crLeftHandSideMatrix1*crLeftHandSideMatrix30 - crLeftHandSideMatrix35*crLeftHandSideMatrix36);
    rLeftHandSideMatrix(2,1)+=-w_g*(DN(0,1)*N[0] - crLeftHandSideMatrix37*crLeftHandSideMatrix4 - crLeftHandSideMatrix38*crLeftHandSideMatrix39);
    rLeftHandSideMatrix(2,2)+=-crLeftHandSideMatrix40*(crLeftHandSideMatrix1 + crLeftHandSideMatrix4);
    rLeftHandSideMatrix(2,3)+=-w_g*(DN(1,0)*N[0] - crLeftHandSideMatrix36*crLeftHandSideMatrix46 - crLeftHandSideMatrix41);
    rLeftHandSideMatrix(2,4)+=-w_g*(DN(1,1)*N[0] - crLeftHandSideMatrix39*crLeftHandSideMatrix48 - crLeftHandSideMatrix47);
    rLeftHandSideMatrix(2,5)+=crLeftHandSideMatrix49;
    rLeftHandSideMatrix(2,6)+=-w_g*(DN(2,0)*N[0] - crLeftHandSideMatrix36*crLeftHandSideMatrix54 - crLeftHandSideMatrix50);
    rLeftHandSideMatrix(2,7)+=-w_g*(DN(2,1)*N[0] - crLeftHandSideMatrix39*crLeftHandSideMatrix56 - crLeftHandSideMatrix55);
    rLeftHandSideMatrix(2,8)+=crLeftHandSideMatrix57;
    rLeftHandSideMatrix(3,0)+=crLeftHandSideMatrix15;
    rLeftHandSideMatrix(3,1)+=crLeftHandSideMatrix25;
    rLeftHandSideMatrix(3,2)+=DN(1,0)*N[0]*w_g;
    rLeftHandSideMatrix(3,3)+=-w_g*(crLeftHandSideMatrix58 + crLeftHandSideMatrix60 + 2*crLeftHandSideMatrix62);
    rLeftHandSideMatrix(3,4)+=crLeftHandSideMatrix64;
    rLeftHandSideMatrix(3,5)+=DN(1,0)*N[1]*w_g;
    rLeftHandSideMatrix(3,6)+=crLeftHandSideMatrix70;
    rLeftHandSideMatrix(3,7)+=crLeftHandSideMatrix71;
    rLeftHandSideMatrix(3,8)+=DN(1,0)*N[2]*w_g;
    rLeftHandSideMatrix(4,0)+=crLeftHandSideMatrix16;
    rLeftHandSideMatrix(4,1)+=crLeftHandSideMatrix26;
    rLeftHandSideMatrix(4,2)+=DN(1,1)*N[0]*w_g;
    rLeftHandSideMatrix(4,3)+=crLeftHandSideMatrix64;
    rLeftHandSideMatrix(4,4)+=-w_g*(crLeftHandSideMatrix58 + 2*crLeftHandSideMatrix60 + crLeftHandSideMatrix62);
    rLeftHandSideMatrix(4,5)+=DN(1,1)*N[1]*w_g;
    rLeftHandSideMatrix(4,6)+=crLeftHandSideMatrix73;
    rLeftHandSideMatrix(4,7)+=crLeftHandSideMatrix74;
    rLeftHandSideMatrix(4,8)+=DN(1,1)*N[2]*w_g;
    rLeftHandSideMatrix(5,0)+=-w_g*(DN(0,0)*N[1] - crLeftHandSideMatrix35*crLeftHandSideMatrix75 - crLeftHandSideMatrix41);
    rLeftHandSideMatrix(5,1)+=-w_g*(DN(0,1)*N[1] - crLeftHandSideMatrix38*crLeftHandSideMatrix76 - crLeftHandSideMatrix47);
    rLeftHandSideMatrix(5,2)+=crLeftHandSideMatrix49;
    rLeftHandSideMatrix(5,3)+=-w_g*(DN(1,0)*N[1] - crLeftHandSideMatrix30*crLeftHandSideMatrix59 - crLeftHandSideMatrix46*crLeftHandSideMatrix75);
    rLeftHandSideMatrix(5,4)+=-w_g*(DN(1,1)*N[1] - crLeftHandSideMatrix37*crLeftHandSideMatrix61 - crLeftHandSideMatrix48*crLeftHandSideMatrix76);
    rLeftHandSideMatrix(5,5)+=-crLeftHandSideMatrix40*(crLeftHandSideMatrix59 + crLeftHandSideMatrix61);
    rLeftHandSideMatrix(5,6)+=-w_g*(DN(2,0)*N[1] - crLeftHandSideMatrix54*crLeftHandSideMatrix75 - crLeftHandSideMatrix77);
    rLeftHandSideMatrix(5,7)+=-w_g*(DN(2,1)*N[1] - crLeftHandSideMatrix56*crLeftHandSideMatrix76 - crLeftHandSideMatrix78);
    rLeftHandSideMatrix(5,8)+=crLeftHandSideMatrix79;
    rLeftHandSideMatrix(6,0)+=crLeftHandSideMatrix22;
    rLeftHandSideMatrix(6,1)+=crLeftHandSideMatrix27;
    rLeftHandSideMatrix(6,2)+=DN(2,0)*N[0]*w_g;
    rLeftHandSideMatrix(6,3)+=crLeftHandSideMatrix70;
    rLeftHandSideMatrix(6,4)+=crLeftHandSideMatrix73;
    rLeftHandSideMatrix(6,5)+=DN(2,0)*N[1]*w_g;
    rLeftHandSideMatrix(6,6)+=-w_g*(crLeftHandSideMatrix80 + crLeftHandSideMatrix82 + 2*crLeftHandSideMatrix84);
    rLeftHandSideMatrix(6,7)+=crLeftHandSideMatrix85;
    rLeftHandSideMatrix(6,8)+=DN(2,0)*N[2]*w_g;
    rLeftHandSideMatrix(7,0)+=crLeftHandSideMatrix23;
    rLeftHandSideMatrix(7,1)+=crLeftHandSideMatrix28;
    rLeftHandSideMatrix(7,2)+=DN(2,1)*N[0]*w_g;
    rLeftHandSideMatrix(7,3)+=crLeftHandSideMatrix71;
    rLeftHandSideMatrix(7,4)+=crLeftHandSideMatrix74;
    rLeftHandSideMatrix(7,5)+=DN(2,1)*N[1]*w_g;
    rLeftHandSideMatrix(7,6)+=crLeftHandSideMatrix85;
    rLeftHandSideMatrix(7,7)+=-w_g*(crLeftHandSideMatrix80 + 2*crLeftHandSideMatrix82 + crLeftHandSideMatrix84);
    rLeftHandSideMatrix(7,8)+=DN(2,1)*N[2]*w_g;
    rLeftHandSideMatrix(8,0)+=-w_g*(DN(0,0)*N[2] - crLeftHandSideMatrix35*crLeftHandSideMatrix86 - crLeftHandSideMatrix50);
    rLeftHandSideMatrix(8,1)+=-w_g*(DN(0,1)*N[2] - crLeftHandSideMatrix38*crLeftHandSideMatrix87 - crLeftHandSideMatrix55);
    rLeftHandSideMatrix(8,2)+=crLeftHandSideMatrix57;
    rLeftHandSideMatrix(8,3)+=-w_g*(DN(1,0)*N[2] - crLeftHandSideMatrix46*crLeftHandSideMatrix86 - crLeftHandSideMatrix77);
    rLeftHandSideMatrix(8,4)+=-w_g*(DN(1,1)*N[2] - crLeftHandSideMatrix48*crLeftHandSideMatrix87 - crLeftHandSideMatrix78);
    rLeftHandSideMatrix(8,5)+=crLeftHandSideMatrix79;
    rLeftHandSideMatrix(8,6)+=-w_g*(DN(2,0)*N[2] - crLeftHandSideMatrix30*crLeftHandSideMatrix81 - crLeftHandSideMatrix54*crLeftHandSideMatrix86);
    rLeftHandSideMatrix(8,7)+=-w_g*(DN(2,1)*N[2] - crLeftHandSideMatrix37*crLeftHandSideMatrix83 - crLeftHandSideMatrix56*crLeftHandSideMatrix87);
    rLeftHandSideMatrix(8,8)+=-crLeftHandSideMatrix40*(crLeftHandSideMatrix81 + crLeftHandSideMatrix83);
    
    
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

        const double crLeftHandSideMatrix0 = DN(0,1)*DN(0,1);
    const double crLeftHandSideMatrix1 = N[0]*nu_nodes[0] + N[1]*nu_nodes[1] + N[2]*nu_nodes[2] + N[3]*nu_nodes[3];
    const double crLeftHandSideMatrix2 = crLeftHandSideMatrix0*crLeftHandSideMatrix1;
    const double crLeftHandSideMatrix3 = DN(0,0)*DN(0,0);
    const double crLeftHandSideMatrix4 = crLeftHandSideMatrix1*crLeftHandSideMatrix3;
    const double crLeftHandSideMatrix5 = sigma_gauss*(N[0]*N[0]);
    const double crLeftHandSideMatrix6 = DN(0,2)*DN(0,2);
    const double crLeftHandSideMatrix7 = crLeftHandSideMatrix1*crLeftHandSideMatrix6;
    const double crLeftHandSideMatrix8 = crLeftHandSideMatrix5 + crLeftHandSideMatrix7;
    const double crLeftHandSideMatrix9 = crLeftHandSideMatrix1*w_g;
    const double crLeftHandSideMatrix10 = DN(0,0)*crLeftHandSideMatrix9;
    const double crLeftHandSideMatrix11 = -DN(0,1)*crLeftHandSideMatrix10;
    const double crLeftHandSideMatrix12 = -DN(0,2)*crLeftHandSideMatrix10;
    const double crLeftHandSideMatrix13 = DN(0,1)*DN(1,1);
    const double crLeftHandSideMatrix14 = crLeftHandSideMatrix1*crLeftHandSideMatrix13;
    const double crLeftHandSideMatrix15 = DN(0,0)*DN(1,0);
    const double crLeftHandSideMatrix16 = crLeftHandSideMatrix1*crLeftHandSideMatrix15;
    const double crLeftHandSideMatrix17 = DN(0,2)*DN(1,2);
    const double crLeftHandSideMatrix18 = crLeftHandSideMatrix1*crLeftHandSideMatrix17;
    const double crLeftHandSideMatrix19 = N[0]*sigma_gauss;
    const double crLeftHandSideMatrix20 = N[1]*crLeftHandSideMatrix19;
    const double crLeftHandSideMatrix21 = crLeftHandSideMatrix18 + crLeftHandSideMatrix20;
    const double crLeftHandSideMatrix22 = -w_g*(crLeftHandSideMatrix14 + 2*crLeftHandSideMatrix16 + crLeftHandSideMatrix21);
    const double crLeftHandSideMatrix23 = DN(1,0)*crLeftHandSideMatrix9;
    const double crLeftHandSideMatrix24 = -DN(0,1)*crLeftHandSideMatrix23;
    const double crLeftHandSideMatrix25 = -DN(0,2)*crLeftHandSideMatrix23;
    const double crLeftHandSideMatrix26 = DN(0,1)*DN(2,1);
    const double crLeftHandSideMatrix27 = crLeftHandSideMatrix1*crLeftHandSideMatrix26;
    const double crLeftHandSideMatrix28 = DN(0,0)*DN(2,0);
    const double crLeftHandSideMatrix29 = crLeftHandSideMatrix1*crLeftHandSideMatrix28;
    const double crLeftHandSideMatrix30 = DN(0,2)*DN(2,2);
    const double crLeftHandSideMatrix31 = crLeftHandSideMatrix1*crLeftHandSideMatrix30;
    const double crLeftHandSideMatrix32 = N[2]*crLeftHandSideMatrix19;
    const double crLeftHandSideMatrix33 = crLeftHandSideMatrix31 + crLeftHandSideMatrix32;
    const double crLeftHandSideMatrix34 = -w_g*(crLeftHandSideMatrix27 + 2*crLeftHandSideMatrix29 + crLeftHandSideMatrix33);
    const double crLeftHandSideMatrix35 = DN(2,0)*crLeftHandSideMatrix9;
    const double crLeftHandSideMatrix36 = -DN(0,1)*crLeftHandSideMatrix35;
    const double crLeftHandSideMatrix37 = -DN(0,2)*crLeftHandSideMatrix35;
    const double crLeftHandSideMatrix38 = DN(0,1)*DN(3,1);
    const double crLeftHandSideMatrix39 = crLeftHandSideMatrix1*crLeftHandSideMatrix38;
    const double crLeftHandSideMatrix40 = DN(0,0)*DN(3,0);
    const double crLeftHandSideMatrix41 = crLeftHandSideMatrix1*crLeftHandSideMatrix40;
    const double crLeftHandSideMatrix42 = DN(0,2)*DN(3,2);
    const double crLeftHandSideMatrix43 = crLeftHandSideMatrix1*crLeftHandSideMatrix42;
    const double crLeftHandSideMatrix44 = N[3]*crLeftHandSideMatrix19;
    const double crLeftHandSideMatrix45 = crLeftHandSideMatrix43 + crLeftHandSideMatrix44;
    const double crLeftHandSideMatrix46 = -w_g*(crLeftHandSideMatrix39 + 2*crLeftHandSideMatrix41 + crLeftHandSideMatrix45);
    const double crLeftHandSideMatrix47 = DN(3,0)*crLeftHandSideMatrix9;
    const double crLeftHandSideMatrix48 = -DN(0,1)*crLeftHandSideMatrix47;
    const double crLeftHandSideMatrix49 = -DN(0,2)*crLeftHandSideMatrix47;
    const double crLeftHandSideMatrix50 = DN(0,2)*crLeftHandSideMatrix9;
    const double crLeftHandSideMatrix51 = -DN(0,1)*crLeftHandSideMatrix50;
    const double crLeftHandSideMatrix52 = DN(0,1)*N[0];
    const double crLeftHandSideMatrix53 = -DN(1,1)*crLeftHandSideMatrix10;
    const double crLeftHandSideMatrix54 = -w_g*(2*crLeftHandSideMatrix14 + crLeftHandSideMatrix16 + crLeftHandSideMatrix21);
    const double crLeftHandSideMatrix55 = -DN(1,1)*crLeftHandSideMatrix50;
    const double crLeftHandSideMatrix56 = DN(0,1)*N[1];
    const double crLeftHandSideMatrix57 = -DN(2,1)*crLeftHandSideMatrix10;
    const double crLeftHandSideMatrix58 = -w_g*(2*crLeftHandSideMatrix27 + crLeftHandSideMatrix29 + crLeftHandSideMatrix33);
    const double crLeftHandSideMatrix59 = -DN(2,1)*crLeftHandSideMatrix50;
    const double crLeftHandSideMatrix60 = DN(0,1)*N[2];
    const double crLeftHandSideMatrix61 = -DN(3,1)*crLeftHandSideMatrix10;
    const double crLeftHandSideMatrix62 = -w_g*(2*crLeftHandSideMatrix39 + crLeftHandSideMatrix41 + crLeftHandSideMatrix45);
    const double crLeftHandSideMatrix63 = -DN(3,1)*crLeftHandSideMatrix50;
    const double crLeftHandSideMatrix64 = DN(0,1)*N[3];
    const double crLeftHandSideMatrix65 = DN(0,2)*N[0];
    const double crLeftHandSideMatrix66 = -DN(1,2)*crLeftHandSideMatrix10;
    const double crLeftHandSideMatrix67 = DN(0,1)*crLeftHandSideMatrix9;
    const double crLeftHandSideMatrix68 = -DN(1,2)*crLeftHandSideMatrix67;
    const double crLeftHandSideMatrix69 = -w_g*(crLeftHandSideMatrix14 + crLeftHandSideMatrix16 + 2*crLeftHandSideMatrix18 + crLeftHandSideMatrix20);
    const double crLeftHandSideMatrix70 = DN(0,2)*N[1];
    const double crLeftHandSideMatrix71 = -DN(2,2)*crLeftHandSideMatrix10;
    const double crLeftHandSideMatrix72 = -DN(2,2)*crLeftHandSideMatrix67;
    const double crLeftHandSideMatrix73 = -w_g*(crLeftHandSideMatrix27 + crLeftHandSideMatrix29 + 2*crLeftHandSideMatrix31 + crLeftHandSideMatrix32);
    const double crLeftHandSideMatrix74 = DN(0,2)*N[2];
    const double crLeftHandSideMatrix75 = -DN(3,2)*crLeftHandSideMatrix10;
    const double crLeftHandSideMatrix76 = -DN(3,2)*crLeftHandSideMatrix67;
    const double crLeftHandSideMatrix77 = -w_g*(crLeftHandSideMatrix39 + crLeftHandSideMatrix41 + 2*crLeftHandSideMatrix43 + crLeftHandSideMatrix44);
    const double crLeftHandSideMatrix78 = DN(0,2)*N[3];
    const double crLeftHandSideMatrix79 = DN(0,0)*nu_nodes[0] + DN(1,0)*nu_nodes[1] + DN(2,0)*nu_nodes[2] + DN(3,0)*nu_nodes[3];
    const double crLeftHandSideMatrix80 = crLeftHandSideMatrix79*delta_gauss;
    const double crLeftHandSideMatrix81 = DN(0,1)*nu_nodes[0] + DN(1,1)*nu_nodes[1] + DN(2,1)*nu_nodes[2] + DN(3,1)*nu_nodes[3];
    const double crLeftHandSideMatrix82 = DN(0,1)*crLeftHandSideMatrix81;
    const double crLeftHandSideMatrix83 = DN(0,0)*crLeftHandSideMatrix79;
    const double crLeftHandSideMatrix84 = DN(0,2)*nu_nodes[0] + DN(1,2)*nu_nodes[1] + DN(2,2)*nu_nodes[2] + DN(3,2)*nu_nodes[3];
    const double crLeftHandSideMatrix85 = DN(0,2)*crLeftHandSideMatrix84;
    const double crLeftHandSideMatrix86 = -crLeftHandSideMatrix19;
    const double crLeftHandSideMatrix87 = crLeftHandSideMatrix85 + crLeftHandSideMatrix86;
    const double crLeftHandSideMatrix88 = crLeftHandSideMatrix82 + 2*crLeftHandSideMatrix83 + crLeftHandSideMatrix87;
    const double crLeftHandSideMatrix89 = DN(0,0)*delta_gauss;
    const double crLeftHandSideMatrix90 = crLeftHandSideMatrix81*delta_gauss;
    const double crLeftHandSideMatrix91 = 2*crLeftHandSideMatrix82 + crLeftHandSideMatrix83 + crLeftHandSideMatrix87;
    const double crLeftHandSideMatrix92 = DN(0,1)*delta_gauss;
    const double crLeftHandSideMatrix93 = crLeftHandSideMatrix84*delta_gauss;
    const double crLeftHandSideMatrix94 = crLeftHandSideMatrix82 + crLeftHandSideMatrix83 + 2*crLeftHandSideMatrix85 + crLeftHandSideMatrix86;
    const double crLeftHandSideMatrix95 = DN(0,2)*delta_gauss;
    const double crLeftHandSideMatrix96 = delta_gauss*w_g;
    const double crLeftHandSideMatrix97 = DN(1,1)*crLeftHandSideMatrix81;
    const double crLeftHandSideMatrix98 = DN(1,0)*crLeftHandSideMatrix79;
    const double crLeftHandSideMatrix99 = DN(1,2)*crLeftHandSideMatrix84;
    const double crLeftHandSideMatrix100 = N[1]*sigma_gauss;
    const double crLeftHandSideMatrix101 = -crLeftHandSideMatrix100;
    const double crLeftHandSideMatrix102 = crLeftHandSideMatrix101 + crLeftHandSideMatrix99;
    const double crLeftHandSideMatrix103 = crLeftHandSideMatrix102 + crLeftHandSideMatrix97 + 2*crLeftHandSideMatrix98;
    const double crLeftHandSideMatrix104 = crLeftHandSideMatrix13*crLeftHandSideMatrix80 + crLeftHandSideMatrix17*crLeftHandSideMatrix80;
    const double crLeftHandSideMatrix105 = DN(1,1)*N[0];
    const double crLeftHandSideMatrix106 = crLeftHandSideMatrix102 + 2*crLeftHandSideMatrix97 + crLeftHandSideMatrix98;
    const double crLeftHandSideMatrix107 = crLeftHandSideMatrix15*crLeftHandSideMatrix90 + crLeftHandSideMatrix17*crLeftHandSideMatrix90;
    const double crLeftHandSideMatrix108 = DN(1,2)*N[0];
    const double crLeftHandSideMatrix109 = crLeftHandSideMatrix101 + crLeftHandSideMatrix97 + crLeftHandSideMatrix98 + 2*crLeftHandSideMatrix99;
    const double crLeftHandSideMatrix110 = crLeftHandSideMatrix13*crLeftHandSideMatrix93 + crLeftHandSideMatrix15*crLeftHandSideMatrix93;
    const double crLeftHandSideMatrix111 = -crLeftHandSideMatrix96*(crLeftHandSideMatrix13 + crLeftHandSideMatrix15 + crLeftHandSideMatrix17);
    const double crLeftHandSideMatrix112 = DN(2,1)*crLeftHandSideMatrix81;
    const double crLeftHandSideMatrix113 = DN(2,0)*crLeftHandSideMatrix79;
    const double crLeftHandSideMatrix114 = DN(2,2)*crLeftHandSideMatrix84;
    const double crLeftHandSideMatrix115 = N[2]*sigma_gauss;
    const double crLeftHandSideMatrix116 = -crLeftHandSideMatrix115;
    const double crLeftHandSideMatrix117 = crLeftHandSideMatrix114 + crLeftHandSideMatrix116;
    const double crLeftHandSideMatrix118 = crLeftHandSideMatrix112 + 2*crLeftHandSideMatrix113 + crLeftHandSideMatrix117;
    const double crLeftHandSideMatrix119 = crLeftHandSideMatrix26*crLeftHandSideMatrix80 + crLeftHandSideMatrix30*crLeftHandSideMatrix80;
    const double crLeftHandSideMatrix120 = DN(2,1)*N[0];
    const double crLeftHandSideMatrix121 = 2*crLeftHandSideMatrix112 + crLeftHandSideMatrix113 + crLeftHandSideMatrix117;
    const double crLeftHandSideMatrix122 = crLeftHandSideMatrix28*crLeftHandSideMatrix90 + crLeftHandSideMatrix30*crLeftHandSideMatrix90;
    const double crLeftHandSideMatrix123 = DN(2,2)*N[0];
    const double crLeftHandSideMatrix124 = crLeftHandSideMatrix112 + crLeftHandSideMatrix113 + 2*crLeftHandSideMatrix114 + crLeftHandSideMatrix116;
    const double crLeftHandSideMatrix125 = crLeftHandSideMatrix26*crLeftHandSideMatrix93 + crLeftHandSideMatrix28*crLeftHandSideMatrix93;
    const double crLeftHandSideMatrix126 = -crLeftHandSideMatrix96*(crLeftHandSideMatrix26 + crLeftHandSideMatrix28 + crLeftHandSideMatrix30);
    const double crLeftHandSideMatrix127 = DN(3,1)*crLeftHandSideMatrix81;
    const double crLeftHandSideMatrix128 = DN(3,0)*crLeftHandSideMatrix79;
    const double crLeftHandSideMatrix129 = DN(3,2)*crLeftHandSideMatrix84;
    const double crLeftHandSideMatrix130 = -N[3]*sigma_gauss;
    const double crLeftHandSideMatrix131 = crLeftHandSideMatrix129 + crLeftHandSideMatrix130;
    const double crLeftHandSideMatrix132 = crLeftHandSideMatrix127 + 2*crLeftHandSideMatrix128 + crLeftHandSideMatrix131;
    const double crLeftHandSideMatrix133 = crLeftHandSideMatrix38*crLeftHandSideMatrix80 + crLeftHandSideMatrix42*crLeftHandSideMatrix80;
    const double crLeftHandSideMatrix134 = DN(3,1)*N[0];
    const double crLeftHandSideMatrix135 = 2*crLeftHandSideMatrix127 + crLeftHandSideMatrix128 + crLeftHandSideMatrix131;
    const double crLeftHandSideMatrix136 = crLeftHandSideMatrix40*crLeftHandSideMatrix90 + crLeftHandSideMatrix42*crLeftHandSideMatrix90;
    const double crLeftHandSideMatrix137 = DN(3,2)*N[0];
    const double crLeftHandSideMatrix138 = crLeftHandSideMatrix127 + crLeftHandSideMatrix128 + 2*crLeftHandSideMatrix129 + crLeftHandSideMatrix130;
    const double crLeftHandSideMatrix139 = crLeftHandSideMatrix38*crLeftHandSideMatrix93 + crLeftHandSideMatrix40*crLeftHandSideMatrix93;
    const double crLeftHandSideMatrix140 = -crLeftHandSideMatrix96*(crLeftHandSideMatrix38 + crLeftHandSideMatrix40 + crLeftHandSideMatrix42);
    const double crLeftHandSideMatrix141 = DN(1,1)*DN(1,1);
    const double crLeftHandSideMatrix142 = crLeftHandSideMatrix1*crLeftHandSideMatrix141;
    const double crLeftHandSideMatrix143 = DN(1,0)*DN(1,0);
    const double crLeftHandSideMatrix144 = crLeftHandSideMatrix1*crLeftHandSideMatrix143;
    const double crLeftHandSideMatrix145 = sigma_gauss*(N[1]*N[1]);
    const double crLeftHandSideMatrix146 = DN(1,2)*DN(1,2);
    const double crLeftHandSideMatrix147 = crLeftHandSideMatrix1*crLeftHandSideMatrix146;
    const double crLeftHandSideMatrix148 = crLeftHandSideMatrix145 + crLeftHandSideMatrix147;
    const double crLeftHandSideMatrix149 = -DN(1,1)*crLeftHandSideMatrix23;
    const double crLeftHandSideMatrix150 = -DN(1,2)*crLeftHandSideMatrix23;
    const double crLeftHandSideMatrix151 = DN(1,1)*DN(2,1);
    const double crLeftHandSideMatrix152 = crLeftHandSideMatrix1*crLeftHandSideMatrix151;
    const double crLeftHandSideMatrix153 = DN(1,0)*DN(2,0);
    const double crLeftHandSideMatrix154 = crLeftHandSideMatrix1*crLeftHandSideMatrix153;
    const double crLeftHandSideMatrix155 = DN(1,2)*DN(2,2);
    const double crLeftHandSideMatrix156 = crLeftHandSideMatrix1*crLeftHandSideMatrix155;
    const double crLeftHandSideMatrix157 = N[2]*crLeftHandSideMatrix100;
    const double crLeftHandSideMatrix158 = crLeftHandSideMatrix156 + crLeftHandSideMatrix157;
    const double crLeftHandSideMatrix159 = -w_g*(crLeftHandSideMatrix152 + 2*crLeftHandSideMatrix154 + crLeftHandSideMatrix158);
    const double crLeftHandSideMatrix160 = -DN(1,1)*crLeftHandSideMatrix35;
    const double crLeftHandSideMatrix161 = -DN(1,2)*crLeftHandSideMatrix35;
    const double crLeftHandSideMatrix162 = DN(1,1)*DN(3,1);
    const double crLeftHandSideMatrix163 = crLeftHandSideMatrix1*crLeftHandSideMatrix162;
    const double crLeftHandSideMatrix164 = DN(1,0)*DN(3,0);
    const double crLeftHandSideMatrix165 = crLeftHandSideMatrix1*crLeftHandSideMatrix164;
    const double crLeftHandSideMatrix166 = DN(1,2)*DN(3,2);
    const double crLeftHandSideMatrix167 = crLeftHandSideMatrix1*crLeftHandSideMatrix166;
    const double crLeftHandSideMatrix168 = N[3]*crLeftHandSideMatrix100;
    const double crLeftHandSideMatrix169 = crLeftHandSideMatrix167 + crLeftHandSideMatrix168;
    const double crLeftHandSideMatrix170 = -w_g*(crLeftHandSideMatrix163 + 2*crLeftHandSideMatrix165 + crLeftHandSideMatrix169);
    const double crLeftHandSideMatrix171 = -DN(1,1)*crLeftHandSideMatrix47;
    const double crLeftHandSideMatrix172 = -DN(1,2)*crLeftHandSideMatrix47;
    const double crLeftHandSideMatrix173 = DN(1,2)*crLeftHandSideMatrix9;
    const double crLeftHandSideMatrix174 = -DN(1,1)*crLeftHandSideMatrix173;
    const double crLeftHandSideMatrix175 = DN(1,1)*N[1];
    const double crLeftHandSideMatrix176 = -DN(2,1)*crLeftHandSideMatrix23;
    const double crLeftHandSideMatrix177 = -w_g*(2*crLeftHandSideMatrix152 + crLeftHandSideMatrix154 + crLeftHandSideMatrix158);
    const double crLeftHandSideMatrix178 = -DN(2,1)*crLeftHandSideMatrix173;
    const double crLeftHandSideMatrix179 = DN(1,1)*N[2];
    const double crLeftHandSideMatrix180 = -DN(3,1)*crLeftHandSideMatrix23;
    const double crLeftHandSideMatrix181 = -w_g*(2*crLeftHandSideMatrix163 + crLeftHandSideMatrix165 + crLeftHandSideMatrix169);
    const double crLeftHandSideMatrix182 = -DN(3,1)*crLeftHandSideMatrix173;
    const double crLeftHandSideMatrix183 = DN(1,1)*N[3];
    const double crLeftHandSideMatrix184 = DN(1,2)*N[1];
    const double crLeftHandSideMatrix185 = -DN(2,2)*crLeftHandSideMatrix23;
    const double crLeftHandSideMatrix186 = DN(1,1)*crLeftHandSideMatrix9;
    const double crLeftHandSideMatrix187 = -DN(2,2)*crLeftHandSideMatrix186;
    const double crLeftHandSideMatrix188 = -w_g*(crLeftHandSideMatrix152 + crLeftHandSideMatrix154 + 2*crLeftHandSideMatrix156 + crLeftHandSideMatrix157);
    const double crLeftHandSideMatrix189 = DN(1,2)*N[2];
    const double crLeftHandSideMatrix190 = -DN(3,2)*crLeftHandSideMatrix23;
    const double crLeftHandSideMatrix191 = -DN(3,2)*crLeftHandSideMatrix186;
    const double crLeftHandSideMatrix192 = -w_g*(crLeftHandSideMatrix163 + crLeftHandSideMatrix165 + 2*crLeftHandSideMatrix167 + crLeftHandSideMatrix168);
    const double crLeftHandSideMatrix193 = DN(1,2)*N[3];
    const double crLeftHandSideMatrix194 = DN(1,0)*delta_gauss;
    const double crLeftHandSideMatrix195 = DN(1,1)*delta_gauss;
    const double crLeftHandSideMatrix196 = DN(1,2)*delta_gauss;
    const double crLeftHandSideMatrix197 = crLeftHandSideMatrix151*crLeftHandSideMatrix80 + crLeftHandSideMatrix155*crLeftHandSideMatrix80;
    const double crLeftHandSideMatrix198 = DN(2,1)*N[1];
    const double crLeftHandSideMatrix199 = crLeftHandSideMatrix153*crLeftHandSideMatrix90 + crLeftHandSideMatrix155*crLeftHandSideMatrix90;
    const double crLeftHandSideMatrix200 = DN(2,2)*N[1];
    const double crLeftHandSideMatrix201 = crLeftHandSideMatrix151*crLeftHandSideMatrix93 + crLeftHandSideMatrix153*crLeftHandSideMatrix93;
    const double crLeftHandSideMatrix202 = -crLeftHandSideMatrix96*(crLeftHandSideMatrix151 + crLeftHandSideMatrix153 + crLeftHandSideMatrix155);
    const double crLeftHandSideMatrix203 = crLeftHandSideMatrix162*crLeftHandSideMatrix80 + crLeftHandSideMatrix166*crLeftHandSideMatrix80;
    const double crLeftHandSideMatrix204 = DN(3,1)*N[1];
    const double crLeftHandSideMatrix205 = crLeftHandSideMatrix164*crLeftHandSideMatrix90 + crLeftHandSideMatrix166*crLeftHandSideMatrix90;
    const double crLeftHandSideMatrix206 = DN(3,2)*N[1];
    const double crLeftHandSideMatrix207 = crLeftHandSideMatrix162*crLeftHandSideMatrix93 + crLeftHandSideMatrix164*crLeftHandSideMatrix93;
    const double crLeftHandSideMatrix208 = -crLeftHandSideMatrix96*(crLeftHandSideMatrix162 + crLeftHandSideMatrix164 + crLeftHandSideMatrix166);
    const double crLeftHandSideMatrix209 = DN(2,1)*DN(2,1);
    const double crLeftHandSideMatrix210 = crLeftHandSideMatrix1*crLeftHandSideMatrix209;
    const double crLeftHandSideMatrix211 = DN(2,0)*DN(2,0);
    const double crLeftHandSideMatrix212 = crLeftHandSideMatrix1*crLeftHandSideMatrix211;
    const double crLeftHandSideMatrix213 = sigma_gauss*(N[2]*N[2]);
    const double crLeftHandSideMatrix214 = DN(2,2)*DN(2,2);
    const double crLeftHandSideMatrix215 = crLeftHandSideMatrix1*crLeftHandSideMatrix214;
    const double crLeftHandSideMatrix216 = crLeftHandSideMatrix213 + crLeftHandSideMatrix215;
    const double crLeftHandSideMatrix217 = -DN(2,1)*crLeftHandSideMatrix35;
    const double crLeftHandSideMatrix218 = -DN(2,2)*crLeftHandSideMatrix35;
    const double crLeftHandSideMatrix219 = DN(2,1)*DN(3,1);
    const double crLeftHandSideMatrix220 = crLeftHandSideMatrix1*crLeftHandSideMatrix219;
    const double crLeftHandSideMatrix221 = DN(2,0)*DN(3,0);
    const double crLeftHandSideMatrix222 = crLeftHandSideMatrix1*crLeftHandSideMatrix221;
    const double crLeftHandSideMatrix223 = DN(2,2)*DN(3,2);
    const double crLeftHandSideMatrix224 = crLeftHandSideMatrix1*crLeftHandSideMatrix223;
    const double crLeftHandSideMatrix225 = N[3]*crLeftHandSideMatrix115;
    const double crLeftHandSideMatrix226 = crLeftHandSideMatrix224 + crLeftHandSideMatrix225;
    const double crLeftHandSideMatrix227 = -w_g*(crLeftHandSideMatrix220 + 2*crLeftHandSideMatrix222 + crLeftHandSideMatrix226);
    const double crLeftHandSideMatrix228 = -DN(2,1)*crLeftHandSideMatrix47;
    const double crLeftHandSideMatrix229 = -DN(2,2)*crLeftHandSideMatrix47;
    const double crLeftHandSideMatrix230 = DN(2,2)*crLeftHandSideMatrix9;
    const double crLeftHandSideMatrix231 = -DN(2,1)*crLeftHandSideMatrix230;
    const double crLeftHandSideMatrix232 = DN(2,1)*N[2];
    const double crLeftHandSideMatrix233 = -DN(3,1)*crLeftHandSideMatrix35;
    const double crLeftHandSideMatrix234 = -w_g*(2*crLeftHandSideMatrix220 + crLeftHandSideMatrix222 + crLeftHandSideMatrix226);
    const double crLeftHandSideMatrix235 = -DN(3,1)*crLeftHandSideMatrix230;
    const double crLeftHandSideMatrix236 = DN(2,1)*N[3];
    const double crLeftHandSideMatrix237 = DN(2,2)*N[2];
    const double crLeftHandSideMatrix238 = -DN(3,2)*crLeftHandSideMatrix35;
    const double crLeftHandSideMatrix239 = DN(3,2)*crLeftHandSideMatrix9;
    const double crLeftHandSideMatrix240 = -DN(2,1)*crLeftHandSideMatrix239;
    const double crLeftHandSideMatrix241 = -w_g*(crLeftHandSideMatrix220 + crLeftHandSideMatrix222 + 2*crLeftHandSideMatrix224 + crLeftHandSideMatrix225);
    const double crLeftHandSideMatrix242 = DN(2,2)*N[3];
    const double crLeftHandSideMatrix243 = DN(2,0)*delta_gauss;
    const double crLeftHandSideMatrix244 = DN(2,1)*delta_gauss;
    const double crLeftHandSideMatrix245 = DN(2,2)*delta_gauss;
    const double crLeftHandSideMatrix246 = crLeftHandSideMatrix219*crLeftHandSideMatrix80 + crLeftHandSideMatrix223*crLeftHandSideMatrix80;
    const double crLeftHandSideMatrix247 = DN(3,1)*N[2];
    const double crLeftHandSideMatrix248 = crLeftHandSideMatrix221*crLeftHandSideMatrix90 + crLeftHandSideMatrix223*crLeftHandSideMatrix90;
    const double crLeftHandSideMatrix249 = DN(3,2)*N[2];
    const double crLeftHandSideMatrix250 = crLeftHandSideMatrix219*crLeftHandSideMatrix93 + crLeftHandSideMatrix221*crLeftHandSideMatrix93;
    const double crLeftHandSideMatrix251 = -crLeftHandSideMatrix96*(crLeftHandSideMatrix219 + crLeftHandSideMatrix221 + crLeftHandSideMatrix223);
    const double crLeftHandSideMatrix252 = DN(3,1)*DN(3,1);
    const double crLeftHandSideMatrix253 = crLeftHandSideMatrix1*crLeftHandSideMatrix252;
    const double crLeftHandSideMatrix254 = DN(3,0)*DN(3,0);
    const double crLeftHandSideMatrix255 = crLeftHandSideMatrix1*crLeftHandSideMatrix254;
    const double crLeftHandSideMatrix256 = sigma_gauss*(N[3]*N[3]);
    const double crLeftHandSideMatrix257 = DN(3,2)*DN(3,2);
    const double crLeftHandSideMatrix258 = crLeftHandSideMatrix1*crLeftHandSideMatrix257;
    const double crLeftHandSideMatrix259 = crLeftHandSideMatrix256 + crLeftHandSideMatrix258;
    const double crLeftHandSideMatrix260 = -DN(3,1)*crLeftHandSideMatrix47;
    const double crLeftHandSideMatrix261 = -DN(3,2)*crLeftHandSideMatrix47;
    const double crLeftHandSideMatrix262 = -DN(3,1)*crLeftHandSideMatrix239;
    const double crLeftHandSideMatrix263 = DN(3,1)*N[3];
    const double crLeftHandSideMatrix264 = DN(3,2)*N[3];
    const double crLeftHandSideMatrix265 = DN(3,0)*delta_gauss;
    const double crLeftHandSideMatrix266 = DN(3,1)*delta_gauss;
    const double crLeftHandSideMatrix267 = DN(3,2)*delta_gauss;
    rLeftHandSideMatrix(0,0)+=-w_g*(crLeftHandSideMatrix2 + 2*crLeftHandSideMatrix4 + crLeftHandSideMatrix8);
    rLeftHandSideMatrix(0,1)+=crLeftHandSideMatrix11;
    rLeftHandSideMatrix(0,2)+=crLeftHandSideMatrix12;
    rLeftHandSideMatrix(0,3)+=DN(0,0)*N[0]*w_g;
    rLeftHandSideMatrix(0,4)+=crLeftHandSideMatrix22;
    rLeftHandSideMatrix(0,5)+=crLeftHandSideMatrix24;
    rLeftHandSideMatrix(0,6)+=crLeftHandSideMatrix25;
    rLeftHandSideMatrix(0,7)+=DN(0,0)*N[1]*w_g;
    rLeftHandSideMatrix(0,8)+=crLeftHandSideMatrix34;
    rLeftHandSideMatrix(0,9)+=crLeftHandSideMatrix36;
    rLeftHandSideMatrix(0,10)+=crLeftHandSideMatrix37;
    rLeftHandSideMatrix(0,11)+=DN(0,0)*N[2]*w_g;
    rLeftHandSideMatrix(0,12)+=crLeftHandSideMatrix46;
    rLeftHandSideMatrix(0,13)+=crLeftHandSideMatrix48;
    rLeftHandSideMatrix(0,14)+=crLeftHandSideMatrix49;
    rLeftHandSideMatrix(0,15)+=DN(0,0)*N[3]*w_g;
    rLeftHandSideMatrix(1,0)+=crLeftHandSideMatrix11;
    rLeftHandSideMatrix(1,1)+=-w_g*(2*crLeftHandSideMatrix2 + crLeftHandSideMatrix4 + crLeftHandSideMatrix8);
    rLeftHandSideMatrix(1,2)+=crLeftHandSideMatrix51;
    rLeftHandSideMatrix(1,3)+=crLeftHandSideMatrix52*w_g;
    rLeftHandSideMatrix(1,4)+=crLeftHandSideMatrix53;
    rLeftHandSideMatrix(1,5)+=crLeftHandSideMatrix54;
    rLeftHandSideMatrix(1,6)+=crLeftHandSideMatrix55;
    rLeftHandSideMatrix(1,7)+=crLeftHandSideMatrix56*w_g;
    rLeftHandSideMatrix(1,8)+=crLeftHandSideMatrix57;
    rLeftHandSideMatrix(1,9)+=crLeftHandSideMatrix58;
    rLeftHandSideMatrix(1,10)+=crLeftHandSideMatrix59;
    rLeftHandSideMatrix(1,11)+=crLeftHandSideMatrix60*w_g;
    rLeftHandSideMatrix(1,12)+=crLeftHandSideMatrix61;
    rLeftHandSideMatrix(1,13)+=crLeftHandSideMatrix62;
    rLeftHandSideMatrix(1,14)+=crLeftHandSideMatrix63;
    rLeftHandSideMatrix(1,15)+=crLeftHandSideMatrix64*w_g;
    rLeftHandSideMatrix(2,0)+=crLeftHandSideMatrix12;
    rLeftHandSideMatrix(2,1)+=crLeftHandSideMatrix51;
    rLeftHandSideMatrix(2,2)+=-w_g*(crLeftHandSideMatrix2 + crLeftHandSideMatrix4 + crLeftHandSideMatrix5 + 2*crLeftHandSideMatrix7);
    rLeftHandSideMatrix(2,3)+=crLeftHandSideMatrix65*w_g;
    rLeftHandSideMatrix(2,4)+=crLeftHandSideMatrix66;
    rLeftHandSideMatrix(2,5)+=crLeftHandSideMatrix68;
    rLeftHandSideMatrix(2,6)+=crLeftHandSideMatrix69;
    rLeftHandSideMatrix(2,7)+=crLeftHandSideMatrix70*w_g;
    rLeftHandSideMatrix(2,8)+=crLeftHandSideMatrix71;
    rLeftHandSideMatrix(2,9)+=crLeftHandSideMatrix72;
    rLeftHandSideMatrix(2,10)+=crLeftHandSideMatrix73;
    rLeftHandSideMatrix(2,11)+=crLeftHandSideMatrix74*w_g;
    rLeftHandSideMatrix(2,12)+=crLeftHandSideMatrix75;
    rLeftHandSideMatrix(2,13)+=crLeftHandSideMatrix76;
    rLeftHandSideMatrix(2,14)+=crLeftHandSideMatrix77;
    rLeftHandSideMatrix(2,15)+=crLeftHandSideMatrix78*w_g;
    rLeftHandSideMatrix(3,0)+=-w_g*(DN(0,0)*N[0] - crLeftHandSideMatrix0*crLeftHandSideMatrix80 - crLeftHandSideMatrix6*crLeftHandSideMatrix80 - crLeftHandSideMatrix88*crLeftHandSideMatrix89);
    rLeftHandSideMatrix(3,1)+=w_g*(crLeftHandSideMatrix3*crLeftHandSideMatrix90 - crLeftHandSideMatrix52 + crLeftHandSideMatrix6*crLeftHandSideMatrix90 + crLeftHandSideMatrix91*crLeftHandSideMatrix92);
    rLeftHandSideMatrix(3,2)+=w_g*(crLeftHandSideMatrix0*crLeftHandSideMatrix93 + crLeftHandSideMatrix3*crLeftHandSideMatrix93 - crLeftHandSideMatrix65 + crLeftHandSideMatrix94*crLeftHandSideMatrix95);
    rLeftHandSideMatrix(3,3)+=-crLeftHandSideMatrix96*(crLeftHandSideMatrix0 + crLeftHandSideMatrix3 + crLeftHandSideMatrix6);
    rLeftHandSideMatrix(3,4)+=-w_g*(DN(1,0)*N[0] - crLeftHandSideMatrix103*crLeftHandSideMatrix89 - crLeftHandSideMatrix104);
    rLeftHandSideMatrix(3,5)+=w_g*(-crLeftHandSideMatrix105 + crLeftHandSideMatrix106*crLeftHandSideMatrix92 + crLeftHandSideMatrix107);
    rLeftHandSideMatrix(3,6)+=w_g*(-crLeftHandSideMatrix108 + crLeftHandSideMatrix109*crLeftHandSideMatrix95 + crLeftHandSideMatrix110);
    rLeftHandSideMatrix(3,7)+=crLeftHandSideMatrix111;
    rLeftHandSideMatrix(3,8)+=-w_g*(DN(2,0)*N[0] - crLeftHandSideMatrix118*crLeftHandSideMatrix89 - crLeftHandSideMatrix119);
    rLeftHandSideMatrix(3,9)+=w_g*(-crLeftHandSideMatrix120 + crLeftHandSideMatrix121*crLeftHandSideMatrix92 + crLeftHandSideMatrix122);
    rLeftHandSideMatrix(3,10)+=w_g*(-crLeftHandSideMatrix123 + crLeftHandSideMatrix124*crLeftHandSideMatrix95 + crLeftHandSideMatrix125);
    rLeftHandSideMatrix(3,11)+=crLeftHandSideMatrix126;
    rLeftHandSideMatrix(3,12)+=-w_g*(DN(3,0)*N[0] - crLeftHandSideMatrix132*crLeftHandSideMatrix89 - crLeftHandSideMatrix133);
    rLeftHandSideMatrix(3,13)+=w_g*(-crLeftHandSideMatrix134 + crLeftHandSideMatrix135*crLeftHandSideMatrix92 + crLeftHandSideMatrix136);
    rLeftHandSideMatrix(3,14)+=w_g*(-crLeftHandSideMatrix137 + crLeftHandSideMatrix138*crLeftHandSideMatrix95 + crLeftHandSideMatrix139);
    rLeftHandSideMatrix(3,15)+=crLeftHandSideMatrix140;
    rLeftHandSideMatrix(4,0)+=crLeftHandSideMatrix22;
    rLeftHandSideMatrix(4,1)+=crLeftHandSideMatrix53;
    rLeftHandSideMatrix(4,2)+=crLeftHandSideMatrix66;
    rLeftHandSideMatrix(4,3)+=DN(1,0)*N[0]*w_g;
    rLeftHandSideMatrix(4,4)+=-w_g*(crLeftHandSideMatrix142 + 2*crLeftHandSideMatrix144 + crLeftHandSideMatrix148);
    rLeftHandSideMatrix(4,5)+=crLeftHandSideMatrix149;
    rLeftHandSideMatrix(4,6)+=crLeftHandSideMatrix150;
    rLeftHandSideMatrix(4,7)+=DN(1,0)*N[1]*w_g;
    rLeftHandSideMatrix(4,8)+=crLeftHandSideMatrix159;
    rLeftHandSideMatrix(4,9)+=crLeftHandSideMatrix160;
    rLeftHandSideMatrix(4,10)+=crLeftHandSideMatrix161;
    rLeftHandSideMatrix(4,11)+=DN(1,0)*N[2]*w_g;
    rLeftHandSideMatrix(4,12)+=crLeftHandSideMatrix170;
    rLeftHandSideMatrix(4,13)+=crLeftHandSideMatrix171;
    rLeftHandSideMatrix(4,14)+=crLeftHandSideMatrix172;
    rLeftHandSideMatrix(4,15)+=DN(1,0)*N[3]*w_g;
    rLeftHandSideMatrix(5,0)+=crLeftHandSideMatrix24;
    rLeftHandSideMatrix(5,1)+=crLeftHandSideMatrix54;
    rLeftHandSideMatrix(5,2)+=crLeftHandSideMatrix68;
    rLeftHandSideMatrix(5,3)+=crLeftHandSideMatrix105*w_g;
    rLeftHandSideMatrix(5,4)+=crLeftHandSideMatrix149;
    rLeftHandSideMatrix(5,5)+=-w_g*(2*crLeftHandSideMatrix142 + crLeftHandSideMatrix144 + crLeftHandSideMatrix148);
    rLeftHandSideMatrix(5,6)+=crLeftHandSideMatrix174;
    rLeftHandSideMatrix(5,7)+=crLeftHandSideMatrix175*w_g;
    rLeftHandSideMatrix(5,8)+=crLeftHandSideMatrix176;
    rLeftHandSideMatrix(5,9)+=crLeftHandSideMatrix177;
    rLeftHandSideMatrix(5,10)+=crLeftHandSideMatrix178;
    rLeftHandSideMatrix(5,11)+=crLeftHandSideMatrix179*w_g;
    rLeftHandSideMatrix(5,12)+=crLeftHandSideMatrix180;
    rLeftHandSideMatrix(5,13)+=crLeftHandSideMatrix181;
    rLeftHandSideMatrix(5,14)+=crLeftHandSideMatrix182;
    rLeftHandSideMatrix(5,15)+=crLeftHandSideMatrix183*w_g;
    rLeftHandSideMatrix(6,0)+=crLeftHandSideMatrix25;
    rLeftHandSideMatrix(6,1)+=crLeftHandSideMatrix55;
    rLeftHandSideMatrix(6,2)+=crLeftHandSideMatrix69;
    rLeftHandSideMatrix(6,3)+=crLeftHandSideMatrix108*w_g;
    rLeftHandSideMatrix(6,4)+=crLeftHandSideMatrix150;
    rLeftHandSideMatrix(6,5)+=crLeftHandSideMatrix174;
    rLeftHandSideMatrix(6,6)+=-w_g*(crLeftHandSideMatrix142 + crLeftHandSideMatrix144 + crLeftHandSideMatrix145 + 2*crLeftHandSideMatrix147);
    rLeftHandSideMatrix(6,7)+=crLeftHandSideMatrix184*w_g;
    rLeftHandSideMatrix(6,8)+=crLeftHandSideMatrix185;
    rLeftHandSideMatrix(6,9)+=crLeftHandSideMatrix187;
    rLeftHandSideMatrix(6,10)+=crLeftHandSideMatrix188;
    rLeftHandSideMatrix(6,11)+=crLeftHandSideMatrix189*w_g;
    rLeftHandSideMatrix(6,12)+=crLeftHandSideMatrix190;
    rLeftHandSideMatrix(6,13)+=crLeftHandSideMatrix191;
    rLeftHandSideMatrix(6,14)+=crLeftHandSideMatrix192;
    rLeftHandSideMatrix(6,15)+=crLeftHandSideMatrix193*w_g;
    rLeftHandSideMatrix(7,0)+=-w_g*(DN(0,0)*N[1] - crLeftHandSideMatrix104 - crLeftHandSideMatrix194*crLeftHandSideMatrix88);
    rLeftHandSideMatrix(7,1)+=w_g*(crLeftHandSideMatrix107 + crLeftHandSideMatrix195*crLeftHandSideMatrix91 - crLeftHandSideMatrix56);
    rLeftHandSideMatrix(7,2)+=w_g*(crLeftHandSideMatrix110 + crLeftHandSideMatrix196*crLeftHandSideMatrix94 - crLeftHandSideMatrix70);
    rLeftHandSideMatrix(7,3)+=crLeftHandSideMatrix111;
    rLeftHandSideMatrix(7,4)+=-w_g*(DN(1,0)*N[1] - crLeftHandSideMatrix103*crLeftHandSideMatrix194 - crLeftHandSideMatrix141*crLeftHandSideMatrix80 - crLeftHandSideMatrix146*crLeftHandSideMatrix80);
    rLeftHandSideMatrix(7,5)+=w_g*(crLeftHandSideMatrix106*crLeftHandSideMatrix195 + crLeftHandSideMatrix143*crLeftHandSideMatrix90 + crLeftHandSideMatrix146*crLeftHandSideMatrix90 - crLeftHandSideMatrix175);
    rLeftHandSideMatrix(7,6)+=w_g*(crLeftHandSideMatrix109*crLeftHandSideMatrix196 + crLeftHandSideMatrix141*crLeftHandSideMatrix93 + crLeftHandSideMatrix143*crLeftHandSideMatrix93 - crLeftHandSideMatrix184);
    rLeftHandSideMatrix(7,7)+=-crLeftHandSideMatrix96*(crLeftHandSideMatrix141 + crLeftHandSideMatrix143 + crLeftHandSideMatrix146);
    rLeftHandSideMatrix(7,8)+=-w_g*(DN(2,0)*N[1] - crLeftHandSideMatrix118*crLeftHandSideMatrix194 - crLeftHandSideMatrix197);
    rLeftHandSideMatrix(7,9)+=w_g*(crLeftHandSideMatrix121*crLeftHandSideMatrix195 - crLeftHandSideMatrix198 + crLeftHandSideMatrix199);
    rLeftHandSideMatrix(7,10)+=w_g*(crLeftHandSideMatrix124*crLeftHandSideMatrix196 - crLeftHandSideMatrix200 + crLeftHandSideMatrix201);
    rLeftHandSideMatrix(7,11)+=crLeftHandSideMatrix202;
    rLeftHandSideMatrix(7,12)+=-w_g*(DN(3,0)*N[1] - crLeftHandSideMatrix132*crLeftHandSideMatrix194 - crLeftHandSideMatrix203);
    rLeftHandSideMatrix(7,13)+=w_g*(crLeftHandSideMatrix135*crLeftHandSideMatrix195 - crLeftHandSideMatrix204 + crLeftHandSideMatrix205);
    rLeftHandSideMatrix(7,14)+=w_g*(crLeftHandSideMatrix138*crLeftHandSideMatrix196 - crLeftHandSideMatrix206 + crLeftHandSideMatrix207);
    rLeftHandSideMatrix(7,15)+=crLeftHandSideMatrix208;
    rLeftHandSideMatrix(8,0)+=crLeftHandSideMatrix34;
    rLeftHandSideMatrix(8,1)+=crLeftHandSideMatrix57;
    rLeftHandSideMatrix(8,2)+=crLeftHandSideMatrix71;
    rLeftHandSideMatrix(8,3)+=DN(2,0)*N[0]*w_g;
    rLeftHandSideMatrix(8,4)+=crLeftHandSideMatrix159;
    rLeftHandSideMatrix(8,5)+=crLeftHandSideMatrix176;
    rLeftHandSideMatrix(8,6)+=crLeftHandSideMatrix185;
    rLeftHandSideMatrix(8,7)+=DN(2,0)*N[1]*w_g;
    rLeftHandSideMatrix(8,8)+=-w_g*(crLeftHandSideMatrix210 + 2*crLeftHandSideMatrix212 + crLeftHandSideMatrix216);
    rLeftHandSideMatrix(8,9)+=crLeftHandSideMatrix217;
    rLeftHandSideMatrix(8,10)+=crLeftHandSideMatrix218;
    rLeftHandSideMatrix(8,11)+=DN(2,0)*N[2]*w_g;
    rLeftHandSideMatrix(8,12)+=crLeftHandSideMatrix227;
    rLeftHandSideMatrix(8,13)+=crLeftHandSideMatrix228;
    rLeftHandSideMatrix(8,14)+=crLeftHandSideMatrix229;
    rLeftHandSideMatrix(8,15)+=DN(2,0)*N[3]*w_g;
    rLeftHandSideMatrix(9,0)+=crLeftHandSideMatrix36;
    rLeftHandSideMatrix(9,1)+=crLeftHandSideMatrix58;
    rLeftHandSideMatrix(9,2)+=crLeftHandSideMatrix72;
    rLeftHandSideMatrix(9,3)+=crLeftHandSideMatrix120*w_g;
    rLeftHandSideMatrix(9,4)+=crLeftHandSideMatrix160;
    rLeftHandSideMatrix(9,5)+=crLeftHandSideMatrix177;
    rLeftHandSideMatrix(9,6)+=crLeftHandSideMatrix187;
    rLeftHandSideMatrix(9,7)+=crLeftHandSideMatrix198*w_g;
    rLeftHandSideMatrix(9,8)+=crLeftHandSideMatrix217;
    rLeftHandSideMatrix(9,9)+=-w_g*(2*crLeftHandSideMatrix210 + crLeftHandSideMatrix212 + crLeftHandSideMatrix216);
    rLeftHandSideMatrix(9,10)+=crLeftHandSideMatrix231;
    rLeftHandSideMatrix(9,11)+=crLeftHandSideMatrix232*w_g;
    rLeftHandSideMatrix(9,12)+=crLeftHandSideMatrix233;
    rLeftHandSideMatrix(9,13)+=crLeftHandSideMatrix234;
    rLeftHandSideMatrix(9,14)+=crLeftHandSideMatrix235;
    rLeftHandSideMatrix(9,15)+=crLeftHandSideMatrix236*w_g;
    rLeftHandSideMatrix(10,0)+=crLeftHandSideMatrix37;
    rLeftHandSideMatrix(10,1)+=crLeftHandSideMatrix59;
    rLeftHandSideMatrix(10,2)+=crLeftHandSideMatrix73;
    rLeftHandSideMatrix(10,3)+=crLeftHandSideMatrix123*w_g;
    rLeftHandSideMatrix(10,4)+=crLeftHandSideMatrix161;
    rLeftHandSideMatrix(10,5)+=crLeftHandSideMatrix178;
    rLeftHandSideMatrix(10,6)+=crLeftHandSideMatrix188;
    rLeftHandSideMatrix(10,7)+=crLeftHandSideMatrix200*w_g;
    rLeftHandSideMatrix(10,8)+=crLeftHandSideMatrix218;
    rLeftHandSideMatrix(10,9)+=crLeftHandSideMatrix231;
    rLeftHandSideMatrix(10,10)+=-w_g*(crLeftHandSideMatrix210 + crLeftHandSideMatrix212 + crLeftHandSideMatrix213 + 2*crLeftHandSideMatrix215);
    rLeftHandSideMatrix(10,11)+=crLeftHandSideMatrix237*w_g;
    rLeftHandSideMatrix(10,12)+=crLeftHandSideMatrix238;
    rLeftHandSideMatrix(10,13)+=crLeftHandSideMatrix240;
    rLeftHandSideMatrix(10,14)+=crLeftHandSideMatrix241;
    rLeftHandSideMatrix(10,15)+=crLeftHandSideMatrix242*w_g;
    rLeftHandSideMatrix(11,0)+=-w_g*(DN(0,0)*N[2] - crLeftHandSideMatrix119 - crLeftHandSideMatrix243*crLeftHandSideMatrix88);
    rLeftHandSideMatrix(11,1)+=w_g*(crLeftHandSideMatrix122 + crLeftHandSideMatrix244*crLeftHandSideMatrix91 - crLeftHandSideMatrix60);
    rLeftHandSideMatrix(11,2)+=w_g*(crLeftHandSideMatrix125 + crLeftHandSideMatrix245*crLeftHandSideMatrix94 - crLeftHandSideMatrix74);
    rLeftHandSideMatrix(11,3)+=crLeftHandSideMatrix126;
    rLeftHandSideMatrix(11,4)+=-w_g*(DN(1,0)*N[2] - crLeftHandSideMatrix103*crLeftHandSideMatrix243 - crLeftHandSideMatrix197);
    rLeftHandSideMatrix(11,5)+=w_g*(crLeftHandSideMatrix106*crLeftHandSideMatrix244 - crLeftHandSideMatrix179 + crLeftHandSideMatrix199);
    rLeftHandSideMatrix(11,6)+=w_g*(crLeftHandSideMatrix109*crLeftHandSideMatrix245 - crLeftHandSideMatrix189 + crLeftHandSideMatrix201);
    rLeftHandSideMatrix(11,7)+=crLeftHandSideMatrix202;
    rLeftHandSideMatrix(11,8)+=-w_g*(DN(2,0)*N[2] - crLeftHandSideMatrix118*crLeftHandSideMatrix243 - crLeftHandSideMatrix209*crLeftHandSideMatrix80 - crLeftHandSideMatrix214*crLeftHandSideMatrix80);
    rLeftHandSideMatrix(11,9)+=w_g*(crLeftHandSideMatrix121*crLeftHandSideMatrix244 + crLeftHandSideMatrix211*crLeftHandSideMatrix90 + crLeftHandSideMatrix214*crLeftHandSideMatrix90 - crLeftHandSideMatrix232);
    rLeftHandSideMatrix(11,10)+=w_g*(crLeftHandSideMatrix124*crLeftHandSideMatrix245 + crLeftHandSideMatrix209*crLeftHandSideMatrix93 + crLeftHandSideMatrix211*crLeftHandSideMatrix93 - crLeftHandSideMatrix237);
    rLeftHandSideMatrix(11,11)+=-crLeftHandSideMatrix96*(crLeftHandSideMatrix209 + crLeftHandSideMatrix211 + crLeftHandSideMatrix214);
    rLeftHandSideMatrix(11,12)+=-w_g*(DN(3,0)*N[2] - crLeftHandSideMatrix132*crLeftHandSideMatrix243 - crLeftHandSideMatrix246);
    rLeftHandSideMatrix(11,13)+=w_g*(crLeftHandSideMatrix135*crLeftHandSideMatrix244 - crLeftHandSideMatrix247 + crLeftHandSideMatrix248);
    rLeftHandSideMatrix(11,14)+=w_g*(crLeftHandSideMatrix138*crLeftHandSideMatrix245 - crLeftHandSideMatrix249 + crLeftHandSideMatrix250);
    rLeftHandSideMatrix(11,15)+=crLeftHandSideMatrix251;
    rLeftHandSideMatrix(12,0)+=crLeftHandSideMatrix46;
    rLeftHandSideMatrix(12,1)+=crLeftHandSideMatrix61;
    rLeftHandSideMatrix(12,2)+=crLeftHandSideMatrix75;
    rLeftHandSideMatrix(12,3)+=DN(3,0)*N[0]*w_g;
    rLeftHandSideMatrix(12,4)+=crLeftHandSideMatrix170;
    rLeftHandSideMatrix(12,5)+=crLeftHandSideMatrix180;
    rLeftHandSideMatrix(12,6)+=crLeftHandSideMatrix190;
    rLeftHandSideMatrix(12,7)+=DN(3,0)*N[1]*w_g;
    rLeftHandSideMatrix(12,8)+=crLeftHandSideMatrix227;
    rLeftHandSideMatrix(12,9)+=crLeftHandSideMatrix233;
    rLeftHandSideMatrix(12,10)+=crLeftHandSideMatrix238;
    rLeftHandSideMatrix(12,11)+=DN(3,0)*N[2]*w_g;
    rLeftHandSideMatrix(12,12)+=-w_g*(crLeftHandSideMatrix253 + 2*crLeftHandSideMatrix255 + crLeftHandSideMatrix259);
    rLeftHandSideMatrix(12,13)+=crLeftHandSideMatrix260;
    rLeftHandSideMatrix(12,14)+=crLeftHandSideMatrix261;
    rLeftHandSideMatrix(12,15)+=DN(3,0)*N[3]*w_g;
    rLeftHandSideMatrix(13,0)+=crLeftHandSideMatrix48;
    rLeftHandSideMatrix(13,1)+=crLeftHandSideMatrix62;
    rLeftHandSideMatrix(13,2)+=crLeftHandSideMatrix76;
    rLeftHandSideMatrix(13,3)+=crLeftHandSideMatrix134*w_g;
    rLeftHandSideMatrix(13,4)+=crLeftHandSideMatrix171;
    rLeftHandSideMatrix(13,5)+=crLeftHandSideMatrix181;
    rLeftHandSideMatrix(13,6)+=crLeftHandSideMatrix191;
    rLeftHandSideMatrix(13,7)+=crLeftHandSideMatrix204*w_g;
    rLeftHandSideMatrix(13,8)+=crLeftHandSideMatrix228;
    rLeftHandSideMatrix(13,9)+=crLeftHandSideMatrix234;
    rLeftHandSideMatrix(13,10)+=crLeftHandSideMatrix240;
    rLeftHandSideMatrix(13,11)+=crLeftHandSideMatrix247*w_g;
    rLeftHandSideMatrix(13,12)+=crLeftHandSideMatrix260;
    rLeftHandSideMatrix(13,13)+=-w_g*(2*crLeftHandSideMatrix253 + crLeftHandSideMatrix255 + crLeftHandSideMatrix259);
    rLeftHandSideMatrix(13,14)+=crLeftHandSideMatrix262;
    rLeftHandSideMatrix(13,15)+=crLeftHandSideMatrix263*w_g;
    rLeftHandSideMatrix(14,0)+=crLeftHandSideMatrix49;
    rLeftHandSideMatrix(14,1)+=crLeftHandSideMatrix63;
    rLeftHandSideMatrix(14,2)+=crLeftHandSideMatrix77;
    rLeftHandSideMatrix(14,3)+=crLeftHandSideMatrix137*w_g;
    rLeftHandSideMatrix(14,4)+=crLeftHandSideMatrix172;
    rLeftHandSideMatrix(14,5)+=crLeftHandSideMatrix182;
    rLeftHandSideMatrix(14,6)+=crLeftHandSideMatrix192;
    rLeftHandSideMatrix(14,7)+=crLeftHandSideMatrix206*w_g;
    rLeftHandSideMatrix(14,8)+=crLeftHandSideMatrix229;
    rLeftHandSideMatrix(14,9)+=crLeftHandSideMatrix235;
    rLeftHandSideMatrix(14,10)+=crLeftHandSideMatrix241;
    rLeftHandSideMatrix(14,11)+=crLeftHandSideMatrix249*w_g;
    rLeftHandSideMatrix(14,12)+=crLeftHandSideMatrix261;
    rLeftHandSideMatrix(14,13)+=crLeftHandSideMatrix262;
    rLeftHandSideMatrix(14,14)+=-w_g*(crLeftHandSideMatrix253 + crLeftHandSideMatrix255 + crLeftHandSideMatrix256 + 2*crLeftHandSideMatrix258);
    rLeftHandSideMatrix(14,15)+=crLeftHandSideMatrix264*w_g;
    rLeftHandSideMatrix(15,0)+=-w_g*(DN(0,0)*N[3] - crLeftHandSideMatrix133 - crLeftHandSideMatrix265*crLeftHandSideMatrix88);
    rLeftHandSideMatrix(15,1)+=w_g*(crLeftHandSideMatrix136 + crLeftHandSideMatrix266*crLeftHandSideMatrix91 - crLeftHandSideMatrix64);
    rLeftHandSideMatrix(15,2)+=w_g*(crLeftHandSideMatrix139 + crLeftHandSideMatrix267*crLeftHandSideMatrix94 - crLeftHandSideMatrix78);
    rLeftHandSideMatrix(15,3)+=crLeftHandSideMatrix140;
    rLeftHandSideMatrix(15,4)+=-w_g*(DN(1,0)*N[3] - crLeftHandSideMatrix103*crLeftHandSideMatrix265 - crLeftHandSideMatrix203);
    rLeftHandSideMatrix(15,5)+=w_g*(crLeftHandSideMatrix106*crLeftHandSideMatrix266 - crLeftHandSideMatrix183 + crLeftHandSideMatrix205);
    rLeftHandSideMatrix(15,6)+=w_g*(crLeftHandSideMatrix109*crLeftHandSideMatrix267 - crLeftHandSideMatrix193 + crLeftHandSideMatrix207);
    rLeftHandSideMatrix(15,7)+=crLeftHandSideMatrix208;
    rLeftHandSideMatrix(15,8)+=-w_g*(DN(2,0)*N[3] - crLeftHandSideMatrix118*crLeftHandSideMatrix265 - crLeftHandSideMatrix246);
    rLeftHandSideMatrix(15,9)+=w_g*(crLeftHandSideMatrix121*crLeftHandSideMatrix266 - crLeftHandSideMatrix236 + crLeftHandSideMatrix248);
    rLeftHandSideMatrix(15,10)+=w_g*(crLeftHandSideMatrix124*crLeftHandSideMatrix267 - crLeftHandSideMatrix242 + crLeftHandSideMatrix250);
    rLeftHandSideMatrix(15,11)+=crLeftHandSideMatrix251;
    rLeftHandSideMatrix(15,12)+=-w_g*(DN(3,0)*N[3] - crLeftHandSideMatrix132*crLeftHandSideMatrix265 - crLeftHandSideMatrix252*crLeftHandSideMatrix80 - crLeftHandSideMatrix257*crLeftHandSideMatrix80);
    rLeftHandSideMatrix(15,13)+=w_g*(crLeftHandSideMatrix135*crLeftHandSideMatrix266 + crLeftHandSideMatrix254*crLeftHandSideMatrix90 + crLeftHandSideMatrix257*crLeftHandSideMatrix90 - crLeftHandSideMatrix263);
    rLeftHandSideMatrix(15,14)+=w_g*(crLeftHandSideMatrix138*crLeftHandSideMatrix267 + crLeftHandSideMatrix252*crLeftHandSideMatrix93 + crLeftHandSideMatrix254*crLeftHandSideMatrix93 - crLeftHandSideMatrix264);
    rLeftHandSideMatrix(15,15)+=-crLeftHandSideMatrix96*(crLeftHandSideMatrix252 + crLeftHandSideMatrix254 + crLeftHandSideMatrix257);
    

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

        const double crRightHandSideVector0 = N[0]*p_nodes[0] + N[1]*p_nodes[1] + N[2]*p_nodes[2];
    const double crRightHandSideVector1 = N[0]*f_nodes(0,0) + N[1]*f_nodes(1,0) + N[2]*f_nodes(2,0);
    const double crRightHandSideVector2 = sigma_gauss*(N[0]*u_nodes(0,0) + N[1]*u_nodes(1,0) + N[2]*u_nodes(2,0));
    const double crRightHandSideVector3 = N[0]*nu_nodes[0] + N[1]*nu_nodes[1] + N[2]*nu_nodes[2];
    const double crRightHandSideVector4 = DN(0,0)*u_nodes(0,0) + DN(1,0)*u_nodes(1,0) + DN(2,0)*u_nodes(2,0);
    const double crRightHandSideVector5 = 2*crRightHandSideVector4;
    const double crRightHandSideVector6 = crRightHandSideVector3*crRightHandSideVector5;
    const double crRightHandSideVector7 = DN(0,0)*u_nodes(0,1) + DN(0,1)*u_nodes(0,0) + DN(1,0)*u_nodes(1,1) + DN(1,1)*u_nodes(1,0) + DN(2,0)*u_nodes(2,1) + DN(2,1)*u_nodes(2,0);
    const double crRightHandSideVector8 = crRightHandSideVector3*crRightHandSideVector7;
    const double crRightHandSideVector9 = N[0]*f_nodes(0,1) + N[1]*f_nodes(1,1) + N[2]*f_nodes(2,1);
    const double crRightHandSideVector10 = sigma_gauss*(N[0]*u_nodes(0,1) + N[1]*u_nodes(1,1) + N[2]*u_nodes(2,1));
    const double crRightHandSideVector11 = DN(0,1)*u_nodes(0,1) + DN(1,1)*u_nodes(1,1) + DN(2,1)*u_nodes(2,1);
    const double crRightHandSideVector12 = 2*crRightHandSideVector11;
    const double crRightHandSideVector13 = crRightHandSideVector12*crRightHandSideVector3;
    const double crRightHandSideVector14 = crRightHandSideVector11 + crRightHandSideVector4;
    const double crRightHandSideVector15 = DN(0,1)*nu_nodes[0] + DN(1,1)*nu_nodes[1] + DN(2,1)*nu_nodes[2];
    const double crRightHandSideVector16 = DN(0,0)*nu_nodes[0] + DN(1,0)*nu_nodes[1] + DN(2,0)*nu_nodes[2];
    const double crRightHandSideVector17 = delta_gauss*(-DN(0,0)*p_nodes[0] - DN(1,0)*p_nodes[1] - DN(2,0)*p_nodes[2] + crRightHandSideVector1 + crRightHandSideVector15*crRightHandSideVector7 + crRightHandSideVector16*crRightHandSideVector5 - crRightHandSideVector2);
    const double crRightHandSideVector18 = delta_gauss*(-DN(0,1)*p_nodes[0] - DN(1,1)*p_nodes[1] - DN(2,1)*p_nodes[2] - crRightHandSideVector10 + crRightHandSideVector12*crRightHandSideVector15 + crRightHandSideVector16*crRightHandSideVector7 + crRightHandSideVector9);
    rRightHandSideVector[0]+=w_g*(-DN(0,0)*crRightHandSideVector0 + DN(0,0)*crRightHandSideVector6 + DN(0,1)*crRightHandSideVector8 - N[0]*crRightHandSideVector1 + N[0]*crRightHandSideVector2);
    rRightHandSideVector[1]+=w_g*(DN(0,0)*crRightHandSideVector8 - DN(0,1)*crRightHandSideVector0 + DN(0,1)*crRightHandSideVector13 + N[0]*crRightHandSideVector10 - N[0]*crRightHandSideVector9);
    rRightHandSideVector[2]+=-w_g*(DN(0,0)*crRightHandSideVector17 + DN(0,1)*crRightHandSideVector18 - N[0]*crRightHandSideVector14);
    rRightHandSideVector[3]+=w_g*(-DN(1,0)*crRightHandSideVector0 + DN(1,0)*crRightHandSideVector6 + DN(1,1)*crRightHandSideVector8 - N[1]*crRightHandSideVector1 + N[1]*crRightHandSideVector2);
    rRightHandSideVector[4]+=w_g*(DN(1,0)*crRightHandSideVector8 - DN(1,1)*crRightHandSideVector0 + DN(1,1)*crRightHandSideVector13 + N[1]*crRightHandSideVector10 - N[1]*crRightHandSideVector9);
    rRightHandSideVector[5]+=-w_g*(DN(1,0)*crRightHandSideVector17 + DN(1,1)*crRightHandSideVector18 - N[1]*crRightHandSideVector14);
    rRightHandSideVector[6]+=w_g*(-DN(2,0)*crRightHandSideVector0 + DN(2,0)*crRightHandSideVector6 + DN(2,1)*crRightHandSideVector8 - N[2]*crRightHandSideVector1 + N[2]*crRightHandSideVector2);
    rRightHandSideVector[7]+=w_g*(DN(2,0)*crRightHandSideVector8 - DN(2,1)*crRightHandSideVector0 + DN(2,1)*crRightHandSideVector13 + N[2]*crRightHandSideVector10 - N[2]*crRightHandSideVector9);
    rRightHandSideVector[8]+=-w_g*(DN(2,0)*crRightHandSideVector17 + DN(2,1)*crRightHandSideVector18 - N[2]*crRightHandSideVector14);
    

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

        const double crRightHandSideVector0 = N[0]*p_nodes[0] + N[1]*p_nodes[1] + N[2]*p_nodes[2] + N[3]*p_nodes[3];
    const double crRightHandSideVector1 = N[0]*f_nodes(0,0) + N[1]*f_nodes(1,0) + N[2]*f_nodes(2,0) + N[3]*f_nodes(3,0);
    const double crRightHandSideVector2 = sigma_gauss*(N[0]*u_nodes(0,0) + N[1]*u_nodes(1,0) + N[2]*u_nodes(2,0) + N[3]*u_nodes(3,0));
    const double crRightHandSideVector3 = N[0]*nu_nodes[0] + N[1]*nu_nodes[1] + N[2]*nu_nodes[2] + N[3]*nu_nodes[3];
    const double crRightHandSideVector4 = DN(0,0)*u_nodes(0,0) + DN(1,0)*u_nodes(1,0) + DN(2,0)*u_nodes(2,0) + DN(3,0)*u_nodes(3,0);
    const double crRightHandSideVector5 = 2*crRightHandSideVector4;
    const double crRightHandSideVector6 = crRightHandSideVector3*crRightHandSideVector5;
    const double crRightHandSideVector7 = DN(0,0)*u_nodes(0,1) + DN(0,1)*u_nodes(0,0) + DN(1,0)*u_nodes(1,1) + DN(1,1)*u_nodes(1,0) + DN(2,0)*u_nodes(2,1) + DN(2,1)*u_nodes(2,0) + DN(3,0)*u_nodes(3,1) + DN(3,1)*u_nodes(3,0);
    const double crRightHandSideVector8 = crRightHandSideVector3*crRightHandSideVector7;
    const double crRightHandSideVector9 = DN(0,0)*u_nodes(0,2) + DN(0,2)*u_nodes(0,0) + DN(1,0)*u_nodes(1,2) + DN(1,2)*u_nodes(1,0) + DN(2,0)*u_nodes(2,2) + DN(2,2)*u_nodes(2,0) + DN(3,0)*u_nodes(3,2) + DN(3,2)*u_nodes(3,0);
    const double crRightHandSideVector10 = DN(0,2)*crRightHandSideVector3;
    const double crRightHandSideVector11 = N[0]*f_nodes(0,1) + N[1]*f_nodes(1,1) + N[2]*f_nodes(2,1) + N[3]*f_nodes(3,1);
    const double crRightHandSideVector12 = sigma_gauss*(N[0]*u_nodes(0,1) + N[1]*u_nodes(1,1) + N[2]*u_nodes(2,1) + N[3]*u_nodes(3,1));
    const double crRightHandSideVector13 = DN(0,1)*u_nodes(0,1) + DN(1,1)*u_nodes(1,1) + DN(2,1)*u_nodes(2,1) + DN(3,1)*u_nodes(3,1);
    const double crRightHandSideVector14 = 2*crRightHandSideVector13;
    const double crRightHandSideVector15 = crRightHandSideVector14*crRightHandSideVector3;
    const double crRightHandSideVector16 = DN(0,1)*u_nodes(0,2) + DN(0,2)*u_nodes(0,1) + DN(1,1)*u_nodes(1,2) + DN(1,2)*u_nodes(1,1) + DN(2,1)*u_nodes(2,2) + DN(2,2)*u_nodes(2,1) + DN(3,1)*u_nodes(3,2) + DN(3,2)*u_nodes(3,1);
    const double crRightHandSideVector17 = N[0]*f_nodes(0,2) + N[1]*f_nodes(1,2) + N[2]*f_nodes(2,2) + N[3]*f_nodes(3,2);
    const double crRightHandSideVector18 = sigma_gauss*(N[0]*u_nodes(0,2) + N[1]*u_nodes(1,2) + N[2]*u_nodes(2,2) + N[3]*u_nodes(3,2));
    const double crRightHandSideVector19 = DN(0,2)*u_nodes(0,2) + DN(1,2)*u_nodes(1,2) + DN(2,2)*u_nodes(2,2) + DN(3,2)*u_nodes(3,2);
    const double crRightHandSideVector20 = 2*crRightHandSideVector19;
    const double crRightHandSideVector21 = crRightHandSideVector3*crRightHandSideVector9;
    const double crRightHandSideVector22 = crRightHandSideVector16*crRightHandSideVector3;
    const double crRightHandSideVector23 = crRightHandSideVector13 + crRightHandSideVector19 + crRightHandSideVector4;
    const double crRightHandSideVector24 = DN(0,1)*nu_nodes[0] + DN(1,1)*nu_nodes[1] + DN(2,1)*nu_nodes[2] + DN(3,1)*nu_nodes[3];
    const double crRightHandSideVector25 = DN(0,2)*nu_nodes[0] + DN(1,2)*nu_nodes[1] + DN(2,2)*nu_nodes[2] + DN(3,2)*nu_nodes[3];
    const double crRightHandSideVector26 = DN(0,0)*nu_nodes[0] + DN(1,0)*nu_nodes[1] + DN(2,0)*nu_nodes[2] + DN(3,0)*nu_nodes[3];
    const double crRightHandSideVector27 = delta_gauss*(-DN(0,0)*p_nodes[0] - DN(1,0)*p_nodes[1] - DN(2,0)*p_nodes[2] - DN(3,0)*p_nodes[3] + crRightHandSideVector1 - crRightHandSideVector2 + crRightHandSideVector24*crRightHandSideVector7 + crRightHandSideVector25*crRightHandSideVector9 + crRightHandSideVector26*crRightHandSideVector5);
    const double crRightHandSideVector28 = delta_gauss*(-DN(0,1)*p_nodes[0] - DN(1,1)*p_nodes[1] - DN(2,1)*p_nodes[2] - DN(3,1)*p_nodes[3] + crRightHandSideVector11 - crRightHandSideVector12 + crRightHandSideVector14*crRightHandSideVector24 + crRightHandSideVector16*crRightHandSideVector25 + crRightHandSideVector26*crRightHandSideVector7);
    const double crRightHandSideVector29 = delta_gauss*(-DN(0,2)*p_nodes[0] - DN(1,2)*p_nodes[1] - DN(2,2)*p_nodes[2] - DN(3,2)*p_nodes[3] + crRightHandSideVector16*crRightHandSideVector24 + crRightHandSideVector17 - crRightHandSideVector18 + crRightHandSideVector20*crRightHandSideVector25 + crRightHandSideVector26*crRightHandSideVector9);
    const double crRightHandSideVector30 = crRightHandSideVector20*crRightHandSideVector3;
    rRightHandSideVector[0]+=w_g*(-DN(0,0)*crRightHandSideVector0 + DN(0,0)*crRightHandSideVector6 + DN(0,1)*crRightHandSideVector8 - N[0]*crRightHandSideVector1 + N[0]*crRightHandSideVector2 + crRightHandSideVector10*crRightHandSideVector9);
    rRightHandSideVector[1]+=w_g*(DN(0,0)*crRightHandSideVector8 - DN(0,1)*crRightHandSideVector0 + DN(0,1)*crRightHandSideVector15 - N[0]*crRightHandSideVector11 + N[0]*crRightHandSideVector12 + crRightHandSideVector10*crRightHandSideVector16);
    rRightHandSideVector[2]+=w_g*(DN(0,0)*crRightHandSideVector21 + DN(0,1)*crRightHandSideVector22 - DN(0,2)*crRightHandSideVector0 - N[0]*crRightHandSideVector17 + N[0]*crRightHandSideVector18 + crRightHandSideVector10*crRightHandSideVector20);
    rRightHandSideVector[3]+=-w_g*(DN(0,0)*crRightHandSideVector27 + DN(0,1)*crRightHandSideVector28 + DN(0,2)*crRightHandSideVector29 - N[0]*crRightHandSideVector23);
    rRightHandSideVector[4]+=w_g*(-DN(1,0)*crRightHandSideVector0 + DN(1,0)*crRightHandSideVector6 + DN(1,1)*crRightHandSideVector8 + DN(1,2)*crRightHandSideVector21 - N[1]*crRightHandSideVector1 + N[1]*crRightHandSideVector2);
    rRightHandSideVector[5]+=w_g*(DN(1,0)*crRightHandSideVector8 - DN(1,1)*crRightHandSideVector0 + DN(1,1)*crRightHandSideVector15 + DN(1,2)*crRightHandSideVector22 - N[1]*crRightHandSideVector11 + N[1]*crRightHandSideVector12);
    rRightHandSideVector[6]+=w_g*(DN(1,0)*crRightHandSideVector21 + DN(1,1)*crRightHandSideVector22 - DN(1,2)*crRightHandSideVector0 + DN(1,2)*crRightHandSideVector30 - N[1]*crRightHandSideVector17 + N[1]*crRightHandSideVector18);
    rRightHandSideVector[7]+=-w_g*(DN(1,0)*crRightHandSideVector27 + DN(1,1)*crRightHandSideVector28 + DN(1,2)*crRightHandSideVector29 - N[1]*crRightHandSideVector23);
    rRightHandSideVector[8]+=w_g*(-DN(2,0)*crRightHandSideVector0 + DN(2,0)*crRightHandSideVector6 + DN(2,1)*crRightHandSideVector8 + DN(2,2)*crRightHandSideVector21 - N[2]*crRightHandSideVector1 + N[2]*crRightHandSideVector2);
    rRightHandSideVector[9]+=w_g*(DN(2,0)*crRightHandSideVector8 - DN(2,1)*crRightHandSideVector0 + DN(2,1)*crRightHandSideVector15 + DN(2,2)*crRightHandSideVector22 - N[2]*crRightHandSideVector11 + N[2]*crRightHandSideVector12);
    rRightHandSideVector[10]+=w_g*(DN(2,0)*crRightHandSideVector21 + DN(2,1)*crRightHandSideVector22 - DN(2,2)*crRightHandSideVector0 + DN(2,2)*crRightHandSideVector30 - N[2]*crRightHandSideVector17 + N[2]*crRightHandSideVector18);
    rRightHandSideVector[11]+=-w_g*(DN(2,0)*crRightHandSideVector27 + DN(2,1)*crRightHandSideVector28 + DN(2,2)*crRightHandSideVector29 - N[2]*crRightHandSideVector23);
    rRightHandSideVector[12]+=w_g*(-DN(3,0)*crRightHandSideVector0 + DN(3,0)*crRightHandSideVector6 + DN(3,1)*crRightHandSideVector8 + DN(3,2)*crRightHandSideVector21 - N[3]*crRightHandSideVector1 + N[3]*crRightHandSideVector2);
    rRightHandSideVector[13]+=w_g*(DN(3,0)*crRightHandSideVector8 - DN(3,1)*crRightHandSideVector0 + DN(3,1)*crRightHandSideVector15 + DN(3,2)*crRightHandSideVector22 - N[3]*crRightHandSideVector11 + N[3]*crRightHandSideVector12);
    rRightHandSideVector[14]+=w_g*(DN(3,0)*crRightHandSideVector21 + DN(3,1)*crRightHandSideVector22 - DN(3,2)*crRightHandSideVector0 + DN(3,2)*crRightHandSideVector30 - N[3]*crRightHandSideVector17 + N[3]*crRightHandSideVector18);
    rRightHandSideVector[15]+=-w_g*(DN(3,0)*crRightHandSideVector27 + DN(3,1)*crRightHandSideVector28 + DN(3,2)*crRightHandSideVector29 - N[3]*crRightHandSideVector23);
    
    
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