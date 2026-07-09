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
#include "first_order_stokes_variable_viscosity_bvs_gl.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< unsigned int TDim >
FirstOrderStokesVariableViscosityBvsGl<TDim>::FirstOrderStokesVariableViscosityBvsGl(IndexType NewId)
    : Element(NewId)
{}

template< unsigned int TDim >
FirstOrderStokesVariableViscosityBvsGl<TDim>::FirstOrderStokesVariableViscosityBvsGl(
    IndexType NewId,
    const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes)
{}

template< unsigned int TDim >
FirstOrderStokesVariableViscosityBvsGl<TDim>::FirstOrderStokesVariableViscosityBvsGl(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{}

template< unsigned int TDim >
FirstOrderStokesVariableViscosityBvsGl<TDim>::FirstOrderStokesVariableViscosityBvsGl(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{}

template< unsigned int TDim >
FirstOrderStokesVariableViscosityBvsGl<TDim>::~FirstOrderStokesVariableViscosityBvsGl()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< unsigned int TDim >
Element::Pointer FirstOrderStokesVariableViscosityBvsGl<TDim>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<FirstOrderStokesVariableViscosityBvsGl<TDim>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template< unsigned int TDim >
Element::Pointer FirstOrderStokesVariableViscosityBvsGl<TDim>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<FirstOrderStokesVariableViscosityBvsGl<TDim>>(NewId, pGeom, pProperties);
}

template< unsigned int TDim >
void FirstOrderStokesVariableViscosityBvsGl<TDim>::EquationIdVector(
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
void FirstOrderStokesVariableViscosityBvsGl<TDim>::GetDofList(
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

template< unsigned int TDim>
void FirstOrderStokesVariableViscosityBvsGl<TDim>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(TDim!=2) << "Method Initialize in Element of type FirstOrderStokesVariableViscosityBvsGl is only implemented for 2D yet" << std::endl;

    const auto& rGeom = this->GetGeometry();
    // Get corresponding list of child conditions
    // Start by getting the neighbor conditions of the three nodes of the element
    const GlobalPointersVector<Condition>& rNodeConditionCandidates0 = rGeom[0].GetValue(NEIGHBOUR_CONDITIONS);
    std::vector<IndexType> neighbour_ids_0(rNodeConditionCandidates0.size());
    for (IndexType i = 0; i < rNodeConditionCandidates0.size(); ++i) {
        neighbour_ids_0[i] = rNodeConditionCandidates0[i].Id();
    }
    const GlobalPointersVector<Condition>& rNodeConditionCandidates1 = rGeom[1].GetValue(NEIGHBOUR_CONDITIONS);
    std::vector<IndexType> neighbour_ids_1(rNodeConditionCandidates1.size());
    for (IndexType i = 0; i < rNodeConditionCandidates1.size(); ++i) {
        neighbour_ids_1[i] = rNodeConditionCandidates1[i].Id();
    }
    const GlobalPointersVector<Condition>& rNodeConditionCandidates2 = rGeom[2].GetValue(NEIGHBOUR_CONDITIONS);
    std::vector<IndexType> neighbour_ids_2(rNodeConditionCandidates2.size());
    for (IndexType i = 0; i < rNodeConditionCandidates2.size(); ++i) {
        neighbour_ids_2[i] = rNodeConditionCandidates2[i].Id();
    }
    // Get the intersection each pair of sets of condition IDs
    std::sort(neighbour_ids_0.begin(), neighbour_ids_0.end());
    std::sort(neighbour_ids_1.begin(), neighbour_ids_1.end());
    std::sort(neighbour_ids_2.begin(), neighbour_ids_2.end());
    std::vector<IndexType> common_conditions_0;
    std::set_intersection(neighbour_ids_0.begin(), neighbour_ids_0.end(),
                          neighbour_ids_1.begin(), neighbour_ids_1.end(),
                          std::back_inserter(common_conditions_0));
    std::vector<IndexType> common_conditions_1;
    std::set_intersection(neighbour_ids_1.begin(), neighbour_ids_1.end(),
                          neighbour_ids_2.begin(), neighbour_ids_2.end(),
                          std::back_inserter(common_conditions_1));
    std::vector<IndexType> common_conditions_2;
    std::set_intersection(neighbour_ids_2.begin(), neighbour_ids_2.end(),
                          neighbour_ids_0.begin(), neighbour_ids_0.end(),
                          std::back_inserter(common_conditions_2));

    KRATOS_ERROR_IF(common_conditions_0.size() > 1 || common_conditions_1.size() > 1 || common_conditions_2.size() > 1) << "Multiple conditions found within two single nodes of Element " << this->Id() << std::endl;
    
    // Make sure mChildConditionsList is empty and then fill it with the child condition pointers
    mChildConditionsList.resize(0);
    if (common_conditions_0.size() == 1) {
        for (IndexType i=0; i < rNodeConditionCandidates0.size(); i++) {
            if (rNodeConditionCandidates0(i)->Id() == common_conditions_0[0]) {
                mChildConditionsList.push_back(rNodeConditionCandidates0(i));
                break;
            }
        }
    }
    if (common_conditions_1.size() == 1) {
        for (IndexType i=0; i < rNodeConditionCandidates1.size(); i++) {
            if (rNodeConditionCandidates1(i)->Id() == common_conditions_1[0]) {
                mChildConditionsList.push_back(rNodeConditionCandidates1(i));
                break;
            }
        }
    }
    if (common_conditions_2.size() == 1) {
        for (IndexType i=0; i < rNodeConditionCandidates2.size(); i++) {
            if (rNodeConditionCandidates2(i)->Id() == common_conditions_2[0]) {
                mChildConditionsList.push_back(rNodeConditionCandidates2(i));
                break;
            }
        }
    }
}

template< unsigned int TDim >
void FirstOrderStokesVariableViscosityBvsGl<TDim>::CalculateLocalSystem(
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

    for (IndexType i=0; i < mChildConditionsList.size(); i++) {
        AddConditionGaussPointLeftHandSideContribution(mChildConditionsList(i), aux_data, rLeftHandSideMatrix);
        AddConditionGaussPointRightHandSideContribution(mChildConditionsList(i), aux_data, rRightHandSideVector);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template< unsigned int TDim >
int FirstOrderStokesVariableViscosityBvsGl<TDim>::Check(const ProcessInfo &rCurrentProcessInfo) const
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
const Parameters FirstOrderStokesVariableViscosityBvsGl<TDim>::GetSpecifications() const
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
std::string FirstOrderStokesVariableViscosityBvsGl<TDim>::Info() const
{
    std::stringstream buffer;
    buffer << "FirstOrderStokesVariableViscosityBvsGl" << TDim << "D" << NumNodes << "N #" << this->Id();
    return buffer.str();
}

template< unsigned int TDim >
void FirstOrderStokesVariableViscosityBvsGl<TDim>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info() << std::endl;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private operations

template <unsigned int TDim>
void FirstOrderStokesVariableViscosityBvsGl<TDim>::SetElementData(
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
void FirstOrderStokesVariableViscosityBvsGl<TDim>::CalculateKinematics(
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
void FirstOrderStokesVariableViscosityBvsGl<2>::AddGaussPointLeftHandSideContribution(
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

        const double crLeftHandSideMatrix0 = DN(0,0)*nu_nodes[0] + DN(1,0)*nu_nodes[1] + DN(2,0)*nu_nodes[2];
    const double crLeftHandSideMatrix1 = DN(0,0)*N[0];
    const double crLeftHandSideMatrix2 = DN(0,0)*DN(0,0);
    const double crLeftHandSideMatrix3 = N[0]*nu_nodes[0] + N[1]*nu_nodes[1] + N[2]*nu_nodes[2];
    const double crLeftHandSideMatrix4 = DN(0,1)*DN(0,1);
    const double crLeftHandSideMatrix5 = crLeftHandSideMatrix2*crLeftHandSideMatrix3 + crLeftHandSideMatrix3*crLeftHandSideMatrix4 + sigma_gauss*(N[0]*N[0]);
    const double crLeftHandSideMatrix6 = DN(0,1)*nu_nodes[0] + DN(1,1)*nu_nodes[1] + DN(2,1)*nu_nodes[2];
    const double crLeftHandSideMatrix7 = crLeftHandSideMatrix1*w_g;
    const double crLeftHandSideMatrix8 = DN(1,0)*N[0];
    const double crLeftHandSideMatrix9 = DN(0,0)*DN(1,0);
    const double crLeftHandSideMatrix10 = DN(0,1)*DN(1,1);
    const double crLeftHandSideMatrix11 = N[0]*sigma_gauss;
    const double crLeftHandSideMatrix12 = N[1]*crLeftHandSideMatrix11 + crLeftHandSideMatrix10*crLeftHandSideMatrix3 + crLeftHandSideMatrix3*crLeftHandSideMatrix9;
    const double crLeftHandSideMatrix13 = crLeftHandSideMatrix8*w_g;
    const double crLeftHandSideMatrix14 = DN(0,0)*N[1];
    const double crLeftHandSideMatrix15 = crLeftHandSideMatrix14*w_g;
    const double crLeftHandSideMatrix16 = DN(2,0)*N[0];
    const double crLeftHandSideMatrix17 = DN(0,0)*DN(2,0);
    const double crLeftHandSideMatrix18 = DN(0,1)*DN(2,1);
    const double crLeftHandSideMatrix19 = N[2]*crLeftHandSideMatrix11 + crLeftHandSideMatrix17*crLeftHandSideMatrix3 + crLeftHandSideMatrix18*crLeftHandSideMatrix3;
    const double crLeftHandSideMatrix20 = crLeftHandSideMatrix16*w_g;
    const double crLeftHandSideMatrix21 = DN(0,0)*N[2];
    const double crLeftHandSideMatrix22 = crLeftHandSideMatrix21*w_g;
    const double crLeftHandSideMatrix23 = DN(0,1)*N[0];
    const double crLeftHandSideMatrix24 = crLeftHandSideMatrix23*w_g;
    const double crLeftHandSideMatrix25 = DN(1,1)*N[0];
    const double crLeftHandSideMatrix26 = crLeftHandSideMatrix25*w_g;
    const double crLeftHandSideMatrix27 = DN(0,1)*N[1];
    const double crLeftHandSideMatrix28 = crLeftHandSideMatrix27*w_g;
    const double crLeftHandSideMatrix29 = DN(2,1)*N[0];
    const double crLeftHandSideMatrix30 = crLeftHandSideMatrix29*w_g;
    const double crLeftHandSideMatrix31 = DN(0,1)*N[2];
    const double crLeftHandSideMatrix32 = crLeftHandSideMatrix31*w_g;
    const double crLeftHandSideMatrix33 = 2*crLeftHandSideMatrix0;
    const double crLeftHandSideMatrix34 = -crLeftHandSideMatrix11;
    const double crLeftHandSideMatrix35 = DN(0,0)*crLeftHandSideMatrix33 + crLeftHandSideMatrix34;
    const double crLeftHandSideMatrix36 = 2*crLeftHandSideMatrix6;
    const double crLeftHandSideMatrix37 = -DN(0,1)*crLeftHandSideMatrix36 - crLeftHandSideMatrix34;
    const double crLeftHandSideMatrix38 = delta_gauss*w_g;
    const double crLeftHandSideMatrix39 = crLeftHandSideMatrix10*crLeftHandSideMatrix33;
    const double crLeftHandSideMatrix40 = N[1]*sigma_gauss;
    const double crLeftHandSideMatrix41 = -crLeftHandSideMatrix40;
    const double crLeftHandSideMatrix42 = DN(1,0)*crLeftHandSideMatrix33 + crLeftHandSideMatrix41;
    const double crLeftHandSideMatrix43 = crLeftHandSideMatrix36*crLeftHandSideMatrix9;
    const double crLeftHandSideMatrix44 = -DN(1,1)*crLeftHandSideMatrix36 - crLeftHandSideMatrix41;
    const double crLeftHandSideMatrix45 = -crLeftHandSideMatrix38*(crLeftHandSideMatrix10 + crLeftHandSideMatrix9);
    const double crLeftHandSideMatrix46 = crLeftHandSideMatrix18*crLeftHandSideMatrix33;
    const double crLeftHandSideMatrix47 = -N[2]*sigma_gauss;
    const double crLeftHandSideMatrix48 = DN(2,0)*crLeftHandSideMatrix33 + crLeftHandSideMatrix47;
    const double crLeftHandSideMatrix49 = crLeftHandSideMatrix17*crLeftHandSideMatrix36;
    const double crLeftHandSideMatrix50 = -DN(2,1)*crLeftHandSideMatrix36 - crLeftHandSideMatrix47;
    const double crLeftHandSideMatrix51 = -crLeftHandSideMatrix38*(crLeftHandSideMatrix17 + crLeftHandSideMatrix18);
    const double crLeftHandSideMatrix52 = DN(1,0)*N[1];
    const double crLeftHandSideMatrix53 = DN(1,0)*DN(1,0);
    const double crLeftHandSideMatrix54 = DN(1,1)*DN(1,1);
    const double crLeftHandSideMatrix55 = crLeftHandSideMatrix3*crLeftHandSideMatrix53 + crLeftHandSideMatrix3*crLeftHandSideMatrix54 + sigma_gauss*(N[1]*N[1]);
    const double crLeftHandSideMatrix56 = crLeftHandSideMatrix52*w_g;
    const double crLeftHandSideMatrix57 = DN(2,0)*N[1];
    const double crLeftHandSideMatrix58 = DN(1,0)*DN(2,0);
    const double crLeftHandSideMatrix59 = DN(1,1)*DN(2,1);
    const double crLeftHandSideMatrix60 = N[2]*crLeftHandSideMatrix40 + crLeftHandSideMatrix3*crLeftHandSideMatrix58 + crLeftHandSideMatrix3*crLeftHandSideMatrix59;
    const double crLeftHandSideMatrix61 = crLeftHandSideMatrix57*w_g;
    const double crLeftHandSideMatrix62 = DN(1,0)*N[2];
    const double crLeftHandSideMatrix63 = crLeftHandSideMatrix62*w_g;
    const double crLeftHandSideMatrix64 = DN(1,1)*N[1];
    const double crLeftHandSideMatrix65 = crLeftHandSideMatrix64*w_g;
    const double crLeftHandSideMatrix66 = DN(2,1)*N[1];
    const double crLeftHandSideMatrix67 = crLeftHandSideMatrix66*w_g;
    const double crLeftHandSideMatrix68 = DN(1,1)*N[2];
    const double crLeftHandSideMatrix69 = crLeftHandSideMatrix68*w_g;
    const double crLeftHandSideMatrix70 = -crLeftHandSideMatrix35;
    const double crLeftHandSideMatrix71 = crLeftHandSideMatrix33*crLeftHandSideMatrix59;
    const double crLeftHandSideMatrix72 = crLeftHandSideMatrix36*crLeftHandSideMatrix58;
    const double crLeftHandSideMatrix73 = -crLeftHandSideMatrix38*(crLeftHandSideMatrix58 + crLeftHandSideMatrix59);
    const double crLeftHandSideMatrix74 = DN(2,0)*N[2];
    const double crLeftHandSideMatrix75 = DN(2,0)*DN(2,0);
    const double crLeftHandSideMatrix76 = DN(2,1)*DN(2,1);
    const double crLeftHandSideMatrix77 = crLeftHandSideMatrix3*crLeftHandSideMatrix75 + crLeftHandSideMatrix3*crLeftHandSideMatrix76 + sigma_gauss*(N[2]*N[2]);
    const double crLeftHandSideMatrix78 = crLeftHandSideMatrix74*w_g;
    const double crLeftHandSideMatrix79 = DN(2,1)*N[2];
    const double crLeftHandSideMatrix80 = crLeftHandSideMatrix79*w_g;
    rLeftHandSideMatrix(0,0)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix1 + crLeftHandSideMatrix5);
    rLeftHandSideMatrix(0,1)+=crLeftHandSideMatrix6*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(0,2)+=crLeftHandSideMatrix7;
    rLeftHandSideMatrix(0,3)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix8 + crLeftHandSideMatrix12);
    rLeftHandSideMatrix(0,4)+=crLeftHandSideMatrix13*crLeftHandSideMatrix6;
    rLeftHandSideMatrix(0,5)+=crLeftHandSideMatrix15;
    rLeftHandSideMatrix(0,6)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix16 + crLeftHandSideMatrix19);
    rLeftHandSideMatrix(0,7)+=crLeftHandSideMatrix20*crLeftHandSideMatrix6;
    rLeftHandSideMatrix(0,8)+=crLeftHandSideMatrix22;
    rLeftHandSideMatrix(1,0)+=crLeftHandSideMatrix0*crLeftHandSideMatrix24;
    rLeftHandSideMatrix(1,1)+=-w_g*(-crLeftHandSideMatrix23*crLeftHandSideMatrix6 + crLeftHandSideMatrix5);
    rLeftHandSideMatrix(1,2)+=crLeftHandSideMatrix24;
    rLeftHandSideMatrix(1,3)+=crLeftHandSideMatrix0*crLeftHandSideMatrix26;
    rLeftHandSideMatrix(1,4)+=-w_g*(crLeftHandSideMatrix12 - crLeftHandSideMatrix25*crLeftHandSideMatrix6);
    rLeftHandSideMatrix(1,5)+=crLeftHandSideMatrix28;
    rLeftHandSideMatrix(1,6)+=crLeftHandSideMatrix0*crLeftHandSideMatrix30;
    rLeftHandSideMatrix(1,7)+=-w_g*(crLeftHandSideMatrix19 - crLeftHandSideMatrix29*crLeftHandSideMatrix6);
    rLeftHandSideMatrix(1,8)+=crLeftHandSideMatrix32;
    rLeftHandSideMatrix(2,0)+=-w_g*(crLeftHandSideMatrix1 - delta_gauss*(DN(0,0)*crLeftHandSideMatrix35 + crLeftHandSideMatrix33*crLeftHandSideMatrix4));
    rLeftHandSideMatrix(2,1)+=-w_g*(crLeftHandSideMatrix23 - delta_gauss*(-DN(0,1)*crLeftHandSideMatrix37 + crLeftHandSideMatrix2*crLeftHandSideMatrix36));
    rLeftHandSideMatrix(2,2)+=-crLeftHandSideMatrix38*(crLeftHandSideMatrix2 + crLeftHandSideMatrix4);
    rLeftHandSideMatrix(2,3)+=-w_g*(crLeftHandSideMatrix8 - delta_gauss*(DN(0,0)*crLeftHandSideMatrix42 + crLeftHandSideMatrix39));
    rLeftHandSideMatrix(2,4)+=-w_g*(crLeftHandSideMatrix25 - delta_gauss*(-DN(0,1)*crLeftHandSideMatrix44 + crLeftHandSideMatrix43));
    rLeftHandSideMatrix(2,5)+=crLeftHandSideMatrix45;
    rLeftHandSideMatrix(2,6)+=-w_g*(crLeftHandSideMatrix16 - delta_gauss*(DN(0,0)*crLeftHandSideMatrix48 + crLeftHandSideMatrix46));
    rLeftHandSideMatrix(2,7)+=-w_g*(crLeftHandSideMatrix29 - delta_gauss*(-DN(0,1)*crLeftHandSideMatrix50 + crLeftHandSideMatrix49));
    rLeftHandSideMatrix(2,8)+=crLeftHandSideMatrix51;
    rLeftHandSideMatrix(3,0)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix14 + crLeftHandSideMatrix12);
    rLeftHandSideMatrix(3,1)+=crLeftHandSideMatrix15*crLeftHandSideMatrix6;
    rLeftHandSideMatrix(3,2)+=crLeftHandSideMatrix13;
    rLeftHandSideMatrix(3,3)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix52 + crLeftHandSideMatrix55);
    rLeftHandSideMatrix(3,4)+=crLeftHandSideMatrix56*crLeftHandSideMatrix6;
    rLeftHandSideMatrix(3,5)+=crLeftHandSideMatrix56;
    rLeftHandSideMatrix(3,6)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix57 + crLeftHandSideMatrix60);
    rLeftHandSideMatrix(3,7)+=crLeftHandSideMatrix6*crLeftHandSideMatrix61;
    rLeftHandSideMatrix(3,8)+=crLeftHandSideMatrix63;
    rLeftHandSideMatrix(4,0)+=crLeftHandSideMatrix0*crLeftHandSideMatrix28;
    rLeftHandSideMatrix(4,1)+=-w_g*(crLeftHandSideMatrix12 - crLeftHandSideMatrix27*crLeftHandSideMatrix6);
    rLeftHandSideMatrix(4,2)+=crLeftHandSideMatrix26;
    rLeftHandSideMatrix(4,3)+=crLeftHandSideMatrix0*crLeftHandSideMatrix65;
    rLeftHandSideMatrix(4,4)+=-w_g*(crLeftHandSideMatrix55 - crLeftHandSideMatrix6*crLeftHandSideMatrix64);
    rLeftHandSideMatrix(4,5)+=crLeftHandSideMatrix65;
    rLeftHandSideMatrix(4,6)+=crLeftHandSideMatrix0*crLeftHandSideMatrix67;
    rLeftHandSideMatrix(4,7)+=-w_g*(-crLeftHandSideMatrix6*crLeftHandSideMatrix66 + crLeftHandSideMatrix60);
    rLeftHandSideMatrix(4,8)+=crLeftHandSideMatrix69;
    rLeftHandSideMatrix(5,0)+=-w_g*(crLeftHandSideMatrix14 - delta_gauss*(-DN(1,0)*crLeftHandSideMatrix70 + crLeftHandSideMatrix39));
    rLeftHandSideMatrix(5,1)+=-w_g*(crLeftHandSideMatrix27 - delta_gauss*(-DN(1,1)*crLeftHandSideMatrix37 + crLeftHandSideMatrix43));
    rLeftHandSideMatrix(5,2)+=crLeftHandSideMatrix45;
    rLeftHandSideMatrix(5,3)+=-w_g*(crLeftHandSideMatrix52 - delta_gauss*(DN(1,0)*crLeftHandSideMatrix42 + crLeftHandSideMatrix33*crLeftHandSideMatrix54));
    rLeftHandSideMatrix(5,4)+=-w_g*(crLeftHandSideMatrix64 - delta_gauss*(-DN(1,1)*crLeftHandSideMatrix44 + crLeftHandSideMatrix36*crLeftHandSideMatrix53));
    rLeftHandSideMatrix(5,5)+=-crLeftHandSideMatrix38*(crLeftHandSideMatrix53 + crLeftHandSideMatrix54);
    rLeftHandSideMatrix(5,6)+=-w_g*(crLeftHandSideMatrix57 - delta_gauss*(DN(1,0)*crLeftHandSideMatrix48 + crLeftHandSideMatrix71));
    rLeftHandSideMatrix(5,7)+=-w_g*(crLeftHandSideMatrix66 - delta_gauss*(-DN(1,1)*crLeftHandSideMatrix50 + crLeftHandSideMatrix72));
    rLeftHandSideMatrix(5,8)+=crLeftHandSideMatrix73;
    rLeftHandSideMatrix(6,0)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix21 + crLeftHandSideMatrix19);
    rLeftHandSideMatrix(6,1)+=crLeftHandSideMatrix22*crLeftHandSideMatrix6;
    rLeftHandSideMatrix(6,2)+=crLeftHandSideMatrix20;
    rLeftHandSideMatrix(6,3)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix62 + crLeftHandSideMatrix60);
    rLeftHandSideMatrix(6,4)+=crLeftHandSideMatrix6*crLeftHandSideMatrix63;
    rLeftHandSideMatrix(6,5)+=crLeftHandSideMatrix61;
    rLeftHandSideMatrix(6,6)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix74 + crLeftHandSideMatrix77);
    rLeftHandSideMatrix(6,7)+=crLeftHandSideMatrix6*crLeftHandSideMatrix78;
    rLeftHandSideMatrix(6,8)+=crLeftHandSideMatrix78;
    rLeftHandSideMatrix(7,0)+=crLeftHandSideMatrix0*crLeftHandSideMatrix32;
    rLeftHandSideMatrix(7,1)+=-w_g*(crLeftHandSideMatrix19 - crLeftHandSideMatrix31*crLeftHandSideMatrix6);
    rLeftHandSideMatrix(7,2)+=crLeftHandSideMatrix30;
    rLeftHandSideMatrix(7,3)+=crLeftHandSideMatrix0*crLeftHandSideMatrix69;
    rLeftHandSideMatrix(7,4)+=-w_g*(-crLeftHandSideMatrix6*crLeftHandSideMatrix68 + crLeftHandSideMatrix60);
    rLeftHandSideMatrix(7,5)+=crLeftHandSideMatrix67;
    rLeftHandSideMatrix(7,6)+=crLeftHandSideMatrix0*crLeftHandSideMatrix80;
    rLeftHandSideMatrix(7,7)+=-w_g*(-crLeftHandSideMatrix6*crLeftHandSideMatrix79 + crLeftHandSideMatrix77);
    rLeftHandSideMatrix(7,8)+=crLeftHandSideMatrix80;
    rLeftHandSideMatrix(8,0)+=-w_g*(crLeftHandSideMatrix21 - delta_gauss*(-DN(2,0)*crLeftHandSideMatrix70 + crLeftHandSideMatrix46));
    rLeftHandSideMatrix(8,1)+=-w_g*(crLeftHandSideMatrix31 - delta_gauss*(-DN(2,1)*crLeftHandSideMatrix37 + crLeftHandSideMatrix49));
    rLeftHandSideMatrix(8,2)+=crLeftHandSideMatrix51;
    rLeftHandSideMatrix(8,3)+=-w_g*(crLeftHandSideMatrix62 - delta_gauss*(DN(2,0)*crLeftHandSideMatrix42 + crLeftHandSideMatrix71));
    rLeftHandSideMatrix(8,4)+=-w_g*(crLeftHandSideMatrix68 - delta_gauss*(-DN(2,1)*crLeftHandSideMatrix44 + crLeftHandSideMatrix72));
    rLeftHandSideMatrix(8,5)+=crLeftHandSideMatrix73;
    rLeftHandSideMatrix(8,6)+=-w_g*(crLeftHandSideMatrix74 - delta_gauss*(DN(2,0)*crLeftHandSideMatrix48 + crLeftHandSideMatrix33*crLeftHandSideMatrix76));
    rLeftHandSideMatrix(8,7)+=-w_g*(crLeftHandSideMatrix79 - delta_gauss*(-DN(2,1)*crLeftHandSideMatrix50 + crLeftHandSideMatrix36*crLeftHandSideMatrix75));
    rLeftHandSideMatrix(8,8)+=-crLeftHandSideMatrix38*(crLeftHandSideMatrix75 + crLeftHandSideMatrix76);
    
    
}

template <>
void FirstOrderStokesVariableViscosityBvsGl<2>::AddConditionGaussPointLeftHandSideContribution(
    const GlobalPointer<Condition>& rConditionPointer,
    const ElementDataContainer& rData,
    MatrixType& rLeftHandSideMatrix)
{
    const unsigned int TDim = 2;

    // Get element geometry
    const auto& r_geom = this->GetGeometry();
    const auto& r_condition_geom = rConditionPointer->GetGeometry();

    // Integration rule data: WE GET THE GAUSS POINTS FROM THE CONDITION, NOT THE ELEMENT
    // Note that we use the same for both velocity and pressure interpolations
    const auto integration_points = r_condition_geom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);
    const SizeType n_gauss = integration_points.size();

    std::vector<Vector> NContainer(n_gauss);
    std::vector<Matrix> DNContainer(n_gauss);
    Vector det_J_vect(n_gauss);
    for (IndexType g=0; g < n_gauss; g++) {
        const auto condition_local_coords = integration_points[g].Coordinates();
        Point::CoordinatesArrayType global_coords;
        r_condition_geom.GlobalCoordinates(global_coords, condition_local_coords);
        // Now we have the global coordinates of the gauss point in the condition
        // We can use these coordinates to find the corresponding local coordinates in the parent element
        Point::CoordinatesArrayType local_coords;
        r_geom.PointLocalCoordinates(local_coords, global_coords);
        // Now we have the local coordinates of the gauss point in the parent element
        // We can use these coordinates to evaluate shape functions and their gradients in the parent element
        
        // Get Shape functions evaluation and their derivatives at these local coords
        NContainer[g] = r_geom.ShapeFunctionsValues(NContainer[g], local_coords);
        Matrix DN_De;
        DN_De = r_geom.ShapeFunctionsLocalGradients(DN_De, local_coords);

        // Calculate Jacobians at integration points (with respect to the parent element)
        Matrix J;
        J = r_geom.Jacobian(J, local_coords);

        Matrix inv_J;
        double det_J;
        MathUtils<double>::InvertMatrix(J, inv_J, det_J);
        det_J_vect[g] = det_J;

        DNContainer[g].resize(NumNodes, TDim);
        DNContainer[g] = prod(DN_De, inv_J);
    }

    for (IndexType g=0; g < n_gauss; g++) {
        Vector N = NContainer[g];
        BoundedMatrix<double, NumNodes, TDim> DN = DNContainer[g];
        double w_g = det_J_vect[g] * integration_points[g].Weight();

        // Compute Unit normal
        array_1d<double, 3> n_gauss;
        n_gauss = r_condition_geom.Normal(integration_points[g].Coordinates());
        double A = norm_2(n_gauss);
        n_gauss /= A;

        array_1d<double, NumNodes> nu_nodes;
        BoundedMatrix<double, NumNodes, TDim> u_nodes;
        for (IndexType i = 0; i < NumNodes; ++i) {
            nu_nodes[i] = r_geom[i].FastGetSolutionStepValue(DYNAMIC_VISCOSITY);
            const auto& r_v = r_geom[i].FastGetSolutionStepValue(VELOCITY);
            for (IndexType d = 0; d < TDim; ++d) {
                u_nodes(i, d) = r_v[d];
            }
        }

        const double delta_gauss = rData.Delta;
        
                const double crLeftHandSideMatrix0 = delta_gauss*w_g*(N[0]*nu_nodes[0] + N[1]*nu_nodes[1] + N[2]*nu_nodes[2]);
        const double crLeftHandSideMatrix1 = crLeftHandSideMatrix0*(DN(0,0)*n_gauss[1] - DN(0,1)*n_gauss[0]);
        const double crLeftHandSideMatrix2 = crLeftHandSideMatrix0*(DN(1,0)*n_gauss[1] - DN(1,1)*n_gauss[0]);
        const double crLeftHandSideMatrix3 = crLeftHandSideMatrix0*(DN(2,0)*n_gauss[1] - DN(2,1)*n_gauss[0]);
        rLeftHandSideMatrix(0,0)+=0;
        rLeftHandSideMatrix(0,1)+=0;
        rLeftHandSideMatrix(0,2)+=0;
        rLeftHandSideMatrix(0,3)+=0;
        rLeftHandSideMatrix(0,4)+=0;
        rLeftHandSideMatrix(0,5)+=0;
        rLeftHandSideMatrix(0,6)+=0;
        rLeftHandSideMatrix(0,7)+=0;
        rLeftHandSideMatrix(0,8)+=0;
        rLeftHandSideMatrix(1,0)+=0;
        rLeftHandSideMatrix(1,1)+=0;
        rLeftHandSideMatrix(1,2)+=0;
        rLeftHandSideMatrix(1,3)+=0;
        rLeftHandSideMatrix(1,4)+=0;
        rLeftHandSideMatrix(1,5)+=0;
        rLeftHandSideMatrix(1,6)+=0;
        rLeftHandSideMatrix(1,7)+=0;
        rLeftHandSideMatrix(1,8)+=0;
        rLeftHandSideMatrix(2,0)+=DN(0,1)*crLeftHandSideMatrix1;
        rLeftHandSideMatrix(2,1)+=-DN(0,0)*crLeftHandSideMatrix1;
        rLeftHandSideMatrix(2,2)+=0;
        rLeftHandSideMatrix(2,3)+=DN(1,1)*crLeftHandSideMatrix1;
        rLeftHandSideMatrix(2,4)+=-DN(1,0)*crLeftHandSideMatrix1;
        rLeftHandSideMatrix(2,5)+=0;
        rLeftHandSideMatrix(2,6)+=DN(2,1)*crLeftHandSideMatrix1;
        rLeftHandSideMatrix(2,7)+=-DN(2,0)*crLeftHandSideMatrix1;
        rLeftHandSideMatrix(2,8)+=0;
        rLeftHandSideMatrix(3,0)+=0;
        rLeftHandSideMatrix(3,1)+=0;
        rLeftHandSideMatrix(3,2)+=0;
        rLeftHandSideMatrix(3,3)+=0;
        rLeftHandSideMatrix(3,4)+=0;
        rLeftHandSideMatrix(3,5)+=0;
        rLeftHandSideMatrix(3,6)+=0;
        rLeftHandSideMatrix(3,7)+=0;
        rLeftHandSideMatrix(3,8)+=0;
        rLeftHandSideMatrix(4,0)+=0;
        rLeftHandSideMatrix(4,1)+=0;
        rLeftHandSideMatrix(4,2)+=0;
        rLeftHandSideMatrix(4,3)+=0;
        rLeftHandSideMatrix(4,4)+=0;
        rLeftHandSideMatrix(4,5)+=0;
        rLeftHandSideMatrix(4,6)+=0;
        rLeftHandSideMatrix(4,7)+=0;
        rLeftHandSideMatrix(4,8)+=0;
        rLeftHandSideMatrix(5,0)+=DN(0,1)*crLeftHandSideMatrix2;
        rLeftHandSideMatrix(5,1)+=-DN(0,0)*crLeftHandSideMatrix2;
        rLeftHandSideMatrix(5,2)+=0;
        rLeftHandSideMatrix(5,3)+=DN(1,1)*crLeftHandSideMatrix2;
        rLeftHandSideMatrix(5,4)+=-DN(1,0)*crLeftHandSideMatrix2;
        rLeftHandSideMatrix(5,5)+=0;
        rLeftHandSideMatrix(5,6)+=DN(2,1)*crLeftHandSideMatrix2;
        rLeftHandSideMatrix(5,7)+=-DN(2,0)*crLeftHandSideMatrix2;
        rLeftHandSideMatrix(5,8)+=0;
        rLeftHandSideMatrix(6,0)+=0;
        rLeftHandSideMatrix(6,1)+=0;
        rLeftHandSideMatrix(6,2)+=0;
        rLeftHandSideMatrix(6,3)+=0;
        rLeftHandSideMatrix(6,4)+=0;
        rLeftHandSideMatrix(6,5)+=0;
        rLeftHandSideMatrix(6,6)+=0;
        rLeftHandSideMatrix(6,7)+=0;
        rLeftHandSideMatrix(6,8)+=0;
        rLeftHandSideMatrix(7,0)+=0;
        rLeftHandSideMatrix(7,1)+=0;
        rLeftHandSideMatrix(7,2)+=0;
        rLeftHandSideMatrix(7,3)+=0;
        rLeftHandSideMatrix(7,4)+=0;
        rLeftHandSideMatrix(7,5)+=0;
        rLeftHandSideMatrix(7,6)+=0;
        rLeftHandSideMatrix(7,7)+=0;
        rLeftHandSideMatrix(7,8)+=0;
        rLeftHandSideMatrix(8,0)+=DN(0,1)*crLeftHandSideMatrix3;
        rLeftHandSideMatrix(8,1)+=-DN(0,0)*crLeftHandSideMatrix3;
        rLeftHandSideMatrix(8,2)+=0;
        rLeftHandSideMatrix(8,3)+=DN(1,1)*crLeftHandSideMatrix3;
        rLeftHandSideMatrix(8,4)+=-DN(1,0)*crLeftHandSideMatrix3;
        rLeftHandSideMatrix(8,5)+=0;
        rLeftHandSideMatrix(8,6)+=DN(2,1)*crLeftHandSideMatrix3;
        rLeftHandSideMatrix(8,7)+=-DN(2,0)*crLeftHandSideMatrix3;
        rLeftHandSideMatrix(8,8)+=0;
        

    }
}

template <>
void FirstOrderStokesVariableViscosityBvsGl<3>::AddGaussPointLeftHandSideContribution(
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

        const double crLeftHandSideMatrix0 = DN(0,0)*nu_nodes[0] + DN(1,0)*nu_nodes[1] + DN(2,0)*nu_nodes[2] + DN(3,0)*nu_nodes[3];
    const double crLeftHandSideMatrix1 = DN(0,0)*N[0];
    const double crLeftHandSideMatrix2 = DN(0,0)*DN(0,0);
    const double crLeftHandSideMatrix3 = N[0]*nu_nodes[0] + N[1]*nu_nodes[1] + N[2]*nu_nodes[2] + N[3]*nu_nodes[3];
    const double crLeftHandSideMatrix4 = DN(0,1)*DN(0,1);
    const double crLeftHandSideMatrix5 = DN(0,2)*DN(0,2);
    const double crLeftHandSideMatrix6 = crLeftHandSideMatrix2*crLeftHandSideMatrix3 + crLeftHandSideMatrix3*crLeftHandSideMatrix4 + crLeftHandSideMatrix3*crLeftHandSideMatrix5 + sigma_gauss*(N[0]*N[0]);
    const double crLeftHandSideMatrix7 = DN(0,1)*nu_nodes[0] + DN(1,1)*nu_nodes[1] + DN(2,1)*nu_nodes[2] + DN(3,1)*nu_nodes[3];
    const double crLeftHandSideMatrix8 = crLeftHandSideMatrix1*w_g;
    const double crLeftHandSideMatrix9 = DN(0,2)*nu_nodes[0] + DN(1,2)*nu_nodes[1] + DN(2,2)*nu_nodes[2] + DN(3,2)*nu_nodes[3];
    const double crLeftHandSideMatrix10 = DN(1,0)*N[0];
    const double crLeftHandSideMatrix11 = DN(0,0)*DN(1,0);
    const double crLeftHandSideMatrix12 = DN(0,1)*DN(1,1);
    const double crLeftHandSideMatrix13 = DN(0,2)*DN(1,2);
    const double crLeftHandSideMatrix14 = N[0]*sigma_gauss;
    const double crLeftHandSideMatrix15 = N[1]*crLeftHandSideMatrix14 + crLeftHandSideMatrix11*crLeftHandSideMatrix3 + crLeftHandSideMatrix12*crLeftHandSideMatrix3 + crLeftHandSideMatrix13*crLeftHandSideMatrix3;
    const double crLeftHandSideMatrix16 = crLeftHandSideMatrix10*w_g;
    const double crLeftHandSideMatrix17 = DN(0,0)*N[1];
    const double crLeftHandSideMatrix18 = crLeftHandSideMatrix17*w_g;
    const double crLeftHandSideMatrix19 = DN(2,0)*N[0];
    const double crLeftHandSideMatrix20 = DN(0,0)*DN(2,0);
    const double crLeftHandSideMatrix21 = DN(0,1)*DN(2,1);
    const double crLeftHandSideMatrix22 = DN(0,2)*DN(2,2);
    const double crLeftHandSideMatrix23 = N[2]*crLeftHandSideMatrix14 + crLeftHandSideMatrix20*crLeftHandSideMatrix3 + crLeftHandSideMatrix21*crLeftHandSideMatrix3 + crLeftHandSideMatrix22*crLeftHandSideMatrix3;
    const double crLeftHandSideMatrix24 = crLeftHandSideMatrix19*w_g;
    const double crLeftHandSideMatrix25 = DN(0,0)*N[2];
    const double crLeftHandSideMatrix26 = crLeftHandSideMatrix25*w_g;
    const double crLeftHandSideMatrix27 = DN(3,0)*N[0];
    const double crLeftHandSideMatrix28 = DN(0,0)*DN(3,0);
    const double crLeftHandSideMatrix29 = DN(0,1)*DN(3,1);
    const double crLeftHandSideMatrix30 = DN(0,2)*DN(3,2);
    const double crLeftHandSideMatrix31 = N[3]*crLeftHandSideMatrix14 + crLeftHandSideMatrix28*crLeftHandSideMatrix3 + crLeftHandSideMatrix29*crLeftHandSideMatrix3 + crLeftHandSideMatrix3*crLeftHandSideMatrix30;
    const double crLeftHandSideMatrix32 = crLeftHandSideMatrix27*w_g;
    const double crLeftHandSideMatrix33 = DN(0,0)*N[3];
    const double crLeftHandSideMatrix34 = crLeftHandSideMatrix33*w_g;
    const double crLeftHandSideMatrix35 = DN(0,1)*N[0];
    const double crLeftHandSideMatrix36 = crLeftHandSideMatrix35*w_g;
    const double crLeftHandSideMatrix37 = DN(1,1)*N[0];
    const double crLeftHandSideMatrix38 = crLeftHandSideMatrix37*w_g;
    const double crLeftHandSideMatrix39 = DN(0,1)*N[1];
    const double crLeftHandSideMatrix40 = crLeftHandSideMatrix39*w_g;
    const double crLeftHandSideMatrix41 = DN(2,1)*N[0];
    const double crLeftHandSideMatrix42 = crLeftHandSideMatrix41*w_g;
    const double crLeftHandSideMatrix43 = DN(0,1)*N[2];
    const double crLeftHandSideMatrix44 = crLeftHandSideMatrix43*w_g;
    const double crLeftHandSideMatrix45 = DN(3,1)*N[0];
    const double crLeftHandSideMatrix46 = crLeftHandSideMatrix45*w_g;
    const double crLeftHandSideMatrix47 = DN(0,1)*N[3];
    const double crLeftHandSideMatrix48 = crLeftHandSideMatrix47*w_g;
    const double crLeftHandSideMatrix49 = DN(0,2)*N[0];
    const double crLeftHandSideMatrix50 = crLeftHandSideMatrix49*w_g;
    const double crLeftHandSideMatrix51 = DN(1,2)*N[0];
    const double crLeftHandSideMatrix52 = crLeftHandSideMatrix51*w_g;
    const double crLeftHandSideMatrix53 = DN(0,2)*N[1];
    const double crLeftHandSideMatrix54 = crLeftHandSideMatrix53*w_g;
    const double crLeftHandSideMatrix55 = DN(2,2)*N[0];
    const double crLeftHandSideMatrix56 = crLeftHandSideMatrix55*w_g;
    const double crLeftHandSideMatrix57 = DN(0,2)*N[2];
    const double crLeftHandSideMatrix58 = crLeftHandSideMatrix57*w_g;
    const double crLeftHandSideMatrix59 = DN(3,2)*N[0];
    const double crLeftHandSideMatrix60 = crLeftHandSideMatrix59*w_g;
    const double crLeftHandSideMatrix61 = DN(0,2)*N[3];
    const double crLeftHandSideMatrix62 = crLeftHandSideMatrix61*w_g;
    const double crLeftHandSideMatrix63 = 2*crLeftHandSideMatrix0;
    const double crLeftHandSideMatrix64 = -N[0]*sigma_gauss;
    const double crLeftHandSideMatrix65 = -DN(0,0)*crLeftHandSideMatrix63 - crLeftHandSideMatrix64;
    const double crLeftHandSideMatrix66 = 2*crLeftHandSideMatrix7;
    const double crLeftHandSideMatrix67 = -DN(0,1)*crLeftHandSideMatrix66 - crLeftHandSideMatrix64;
    const double crLeftHandSideMatrix68 = 2*crLeftHandSideMatrix9;
    const double crLeftHandSideMatrix69 = -DN(0,2)*crLeftHandSideMatrix68 - crLeftHandSideMatrix64;
    const double crLeftHandSideMatrix70 = delta_gauss*w_g;
    const double crLeftHandSideMatrix71 = -N[1]*sigma_gauss;
    const double crLeftHandSideMatrix72 = -DN(1,0)*crLeftHandSideMatrix63 - crLeftHandSideMatrix71;
    const double crLeftHandSideMatrix73 = crLeftHandSideMatrix12*crLeftHandSideMatrix63 + crLeftHandSideMatrix13*crLeftHandSideMatrix63;
    const double crLeftHandSideMatrix74 = -DN(1,1)*crLeftHandSideMatrix66 - crLeftHandSideMatrix71;
    const double crLeftHandSideMatrix75 = crLeftHandSideMatrix11*crLeftHandSideMatrix66 + crLeftHandSideMatrix13*crLeftHandSideMatrix66;
    const double crLeftHandSideMatrix76 = -DN(1,2)*crLeftHandSideMatrix68 - crLeftHandSideMatrix71;
    const double crLeftHandSideMatrix77 = crLeftHandSideMatrix11*crLeftHandSideMatrix68 + crLeftHandSideMatrix12*crLeftHandSideMatrix68;
    const double crLeftHandSideMatrix78 = -crLeftHandSideMatrix70*(crLeftHandSideMatrix11 + crLeftHandSideMatrix12 + crLeftHandSideMatrix13);
    const double crLeftHandSideMatrix79 = -N[2]*sigma_gauss;
    const double crLeftHandSideMatrix80 = -DN(2,0)*crLeftHandSideMatrix63 - crLeftHandSideMatrix79;
    const double crLeftHandSideMatrix81 = crLeftHandSideMatrix21*crLeftHandSideMatrix63 + crLeftHandSideMatrix22*crLeftHandSideMatrix63;
    const double crLeftHandSideMatrix82 = -DN(2,1)*crLeftHandSideMatrix66 - crLeftHandSideMatrix79;
    const double crLeftHandSideMatrix83 = crLeftHandSideMatrix20*crLeftHandSideMatrix66 + crLeftHandSideMatrix22*crLeftHandSideMatrix66;
    const double crLeftHandSideMatrix84 = -DN(2,2)*crLeftHandSideMatrix68 - crLeftHandSideMatrix79;
    const double crLeftHandSideMatrix85 = crLeftHandSideMatrix20*crLeftHandSideMatrix68 + crLeftHandSideMatrix21*crLeftHandSideMatrix68;
    const double crLeftHandSideMatrix86 = -crLeftHandSideMatrix70*(crLeftHandSideMatrix20 + crLeftHandSideMatrix21 + crLeftHandSideMatrix22);
    const double crLeftHandSideMatrix87 = -N[3]*sigma_gauss;
    const double crLeftHandSideMatrix88 = -DN(3,0)*crLeftHandSideMatrix63 - crLeftHandSideMatrix87;
    const double crLeftHandSideMatrix89 = crLeftHandSideMatrix29*crLeftHandSideMatrix63 + crLeftHandSideMatrix30*crLeftHandSideMatrix63;
    const double crLeftHandSideMatrix90 = -DN(3,1)*crLeftHandSideMatrix66 - crLeftHandSideMatrix87;
    const double crLeftHandSideMatrix91 = crLeftHandSideMatrix28*crLeftHandSideMatrix66 + crLeftHandSideMatrix30*crLeftHandSideMatrix66;
    const double crLeftHandSideMatrix92 = -DN(3,2)*crLeftHandSideMatrix68 - crLeftHandSideMatrix87;
    const double crLeftHandSideMatrix93 = crLeftHandSideMatrix28*crLeftHandSideMatrix68 + crLeftHandSideMatrix29*crLeftHandSideMatrix68;
    const double crLeftHandSideMatrix94 = -crLeftHandSideMatrix70*(crLeftHandSideMatrix28 + crLeftHandSideMatrix29 + crLeftHandSideMatrix30);
    const double crLeftHandSideMatrix95 = DN(1,0)*N[1];
    const double crLeftHandSideMatrix96 = DN(1,0)*DN(1,0);
    const double crLeftHandSideMatrix97 = DN(1,1)*DN(1,1);
    const double crLeftHandSideMatrix98 = DN(1,2)*DN(1,2);
    const double crLeftHandSideMatrix99 = crLeftHandSideMatrix3*crLeftHandSideMatrix96 + crLeftHandSideMatrix3*crLeftHandSideMatrix97 + crLeftHandSideMatrix3*crLeftHandSideMatrix98 + sigma_gauss*(N[1]*N[1]);
    const double crLeftHandSideMatrix100 = crLeftHandSideMatrix95*w_g;
    const double crLeftHandSideMatrix101 = DN(2,0)*N[1];
    const double crLeftHandSideMatrix102 = DN(1,0)*DN(2,0);
    const double crLeftHandSideMatrix103 = DN(1,1)*DN(2,1);
    const double crLeftHandSideMatrix104 = DN(1,2)*DN(2,2);
    const double crLeftHandSideMatrix105 = N[1]*sigma_gauss;
    const double crLeftHandSideMatrix106 = N[2]*crLeftHandSideMatrix105 + crLeftHandSideMatrix102*crLeftHandSideMatrix3 + crLeftHandSideMatrix103*crLeftHandSideMatrix3 + crLeftHandSideMatrix104*crLeftHandSideMatrix3;
    const double crLeftHandSideMatrix107 = crLeftHandSideMatrix101*w_g;
    const double crLeftHandSideMatrix108 = DN(1,0)*N[2];
    const double crLeftHandSideMatrix109 = crLeftHandSideMatrix108*w_g;
    const double crLeftHandSideMatrix110 = DN(3,0)*N[1];
    const double crLeftHandSideMatrix111 = DN(1,0)*DN(3,0);
    const double crLeftHandSideMatrix112 = DN(1,1)*DN(3,1);
    const double crLeftHandSideMatrix113 = DN(1,2)*DN(3,2);
    const double crLeftHandSideMatrix114 = N[3]*crLeftHandSideMatrix105 + crLeftHandSideMatrix111*crLeftHandSideMatrix3 + crLeftHandSideMatrix112*crLeftHandSideMatrix3 + crLeftHandSideMatrix113*crLeftHandSideMatrix3;
    const double crLeftHandSideMatrix115 = crLeftHandSideMatrix110*w_g;
    const double crLeftHandSideMatrix116 = DN(1,0)*N[3];
    const double crLeftHandSideMatrix117 = crLeftHandSideMatrix116*w_g;
    const double crLeftHandSideMatrix118 = DN(1,1)*N[1];
    const double crLeftHandSideMatrix119 = crLeftHandSideMatrix118*w_g;
    const double crLeftHandSideMatrix120 = DN(2,1)*N[1];
    const double crLeftHandSideMatrix121 = crLeftHandSideMatrix120*w_g;
    const double crLeftHandSideMatrix122 = DN(1,1)*N[2];
    const double crLeftHandSideMatrix123 = crLeftHandSideMatrix122*w_g;
    const double crLeftHandSideMatrix124 = DN(3,1)*N[1];
    const double crLeftHandSideMatrix125 = crLeftHandSideMatrix124*w_g;
    const double crLeftHandSideMatrix126 = DN(1,1)*N[3];
    const double crLeftHandSideMatrix127 = crLeftHandSideMatrix126*w_g;
    const double crLeftHandSideMatrix128 = DN(1,2)*N[1];
    const double crLeftHandSideMatrix129 = crLeftHandSideMatrix128*w_g;
    const double crLeftHandSideMatrix130 = DN(2,2)*N[1];
    const double crLeftHandSideMatrix131 = crLeftHandSideMatrix130*w_g;
    const double crLeftHandSideMatrix132 = DN(1,2)*N[2];
    const double crLeftHandSideMatrix133 = crLeftHandSideMatrix132*w_g;
    const double crLeftHandSideMatrix134 = DN(3,2)*N[1];
    const double crLeftHandSideMatrix135 = crLeftHandSideMatrix134*w_g;
    const double crLeftHandSideMatrix136 = DN(1,2)*N[3];
    const double crLeftHandSideMatrix137 = crLeftHandSideMatrix136*w_g;
    const double crLeftHandSideMatrix138 = crLeftHandSideMatrix103*crLeftHandSideMatrix63 + crLeftHandSideMatrix104*crLeftHandSideMatrix63;
    const double crLeftHandSideMatrix139 = crLeftHandSideMatrix102*crLeftHandSideMatrix66 + crLeftHandSideMatrix104*crLeftHandSideMatrix66;
    const double crLeftHandSideMatrix140 = crLeftHandSideMatrix102*crLeftHandSideMatrix68 + crLeftHandSideMatrix103*crLeftHandSideMatrix68;
    const double crLeftHandSideMatrix141 = -crLeftHandSideMatrix70*(crLeftHandSideMatrix102 + crLeftHandSideMatrix103 + crLeftHandSideMatrix104);
    const double crLeftHandSideMatrix142 = crLeftHandSideMatrix112*crLeftHandSideMatrix63 + crLeftHandSideMatrix113*crLeftHandSideMatrix63;
    const double crLeftHandSideMatrix143 = crLeftHandSideMatrix111*crLeftHandSideMatrix66 + crLeftHandSideMatrix113*crLeftHandSideMatrix66;
    const double crLeftHandSideMatrix144 = crLeftHandSideMatrix111*crLeftHandSideMatrix68 + crLeftHandSideMatrix112*crLeftHandSideMatrix68;
    const double crLeftHandSideMatrix145 = -crLeftHandSideMatrix70*(crLeftHandSideMatrix111 + crLeftHandSideMatrix112 + crLeftHandSideMatrix113);
    const double crLeftHandSideMatrix146 = DN(2,0)*N[2];
    const double crLeftHandSideMatrix147 = DN(2,0)*DN(2,0);
    const double crLeftHandSideMatrix148 = DN(2,1)*DN(2,1);
    const double crLeftHandSideMatrix149 = DN(2,2)*DN(2,2);
    const double crLeftHandSideMatrix150 = crLeftHandSideMatrix147*crLeftHandSideMatrix3 + crLeftHandSideMatrix148*crLeftHandSideMatrix3 + crLeftHandSideMatrix149*crLeftHandSideMatrix3 + sigma_gauss*(N[2]*N[2]);
    const double crLeftHandSideMatrix151 = crLeftHandSideMatrix146*w_g;
    const double crLeftHandSideMatrix152 = DN(3,0)*N[2];
    const double crLeftHandSideMatrix153 = DN(2,0)*DN(3,0);
    const double crLeftHandSideMatrix154 = DN(2,1)*DN(3,1);
    const double crLeftHandSideMatrix155 = DN(2,2)*DN(3,2);
    const double crLeftHandSideMatrix156 = N[2]*N[3]*sigma_gauss + crLeftHandSideMatrix153*crLeftHandSideMatrix3 + crLeftHandSideMatrix154*crLeftHandSideMatrix3 + crLeftHandSideMatrix155*crLeftHandSideMatrix3;
    const double crLeftHandSideMatrix157 = crLeftHandSideMatrix152*w_g;
    const double crLeftHandSideMatrix158 = DN(2,0)*N[3];
    const double crLeftHandSideMatrix159 = crLeftHandSideMatrix158*w_g;
    const double crLeftHandSideMatrix160 = DN(2,1)*N[2];
    const double crLeftHandSideMatrix161 = crLeftHandSideMatrix160*w_g;
    const double crLeftHandSideMatrix162 = DN(3,1)*N[2];
    const double crLeftHandSideMatrix163 = crLeftHandSideMatrix162*w_g;
    const double crLeftHandSideMatrix164 = DN(2,1)*N[3];
    const double crLeftHandSideMatrix165 = crLeftHandSideMatrix164*w_g;
    const double crLeftHandSideMatrix166 = DN(2,2)*N[2];
    const double crLeftHandSideMatrix167 = crLeftHandSideMatrix166*w_g;
    const double crLeftHandSideMatrix168 = DN(3,2)*N[2];
    const double crLeftHandSideMatrix169 = crLeftHandSideMatrix168*w_g;
    const double crLeftHandSideMatrix170 = DN(2,2)*N[3];
    const double crLeftHandSideMatrix171 = crLeftHandSideMatrix170*w_g;
    const double crLeftHandSideMatrix172 = crLeftHandSideMatrix154*crLeftHandSideMatrix63 + crLeftHandSideMatrix155*crLeftHandSideMatrix63;
    const double crLeftHandSideMatrix173 = crLeftHandSideMatrix153*crLeftHandSideMatrix66 + crLeftHandSideMatrix155*crLeftHandSideMatrix66;
    const double crLeftHandSideMatrix174 = crLeftHandSideMatrix153*crLeftHandSideMatrix68 + crLeftHandSideMatrix154*crLeftHandSideMatrix68;
    const double crLeftHandSideMatrix175 = -crLeftHandSideMatrix70*(crLeftHandSideMatrix153 + crLeftHandSideMatrix154 + crLeftHandSideMatrix155);
    const double crLeftHandSideMatrix176 = DN(3,0)*N[3];
    const double crLeftHandSideMatrix177 = DN(3,0)*DN(3,0);
    const double crLeftHandSideMatrix178 = DN(3,1)*DN(3,1);
    const double crLeftHandSideMatrix179 = DN(3,2)*DN(3,2);
    const double crLeftHandSideMatrix180 = crLeftHandSideMatrix177*crLeftHandSideMatrix3 + crLeftHandSideMatrix178*crLeftHandSideMatrix3 + crLeftHandSideMatrix179*crLeftHandSideMatrix3 + sigma_gauss*(N[3]*N[3]);
    const double crLeftHandSideMatrix181 = crLeftHandSideMatrix176*w_g;
    const double crLeftHandSideMatrix182 = DN(3,1)*N[3];
    const double crLeftHandSideMatrix183 = crLeftHandSideMatrix182*w_g;
    const double crLeftHandSideMatrix184 = DN(3,2)*N[3];
    const double crLeftHandSideMatrix185 = crLeftHandSideMatrix184*w_g;
    rLeftHandSideMatrix(0,0)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix1 + crLeftHandSideMatrix6);
    rLeftHandSideMatrix(0,1)+=crLeftHandSideMatrix7*crLeftHandSideMatrix8;
    rLeftHandSideMatrix(0,2)+=crLeftHandSideMatrix8*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(0,3)+=crLeftHandSideMatrix8;
    rLeftHandSideMatrix(0,4)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix10 + crLeftHandSideMatrix15);
    rLeftHandSideMatrix(0,5)+=crLeftHandSideMatrix16*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(0,6)+=crLeftHandSideMatrix16*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(0,7)+=crLeftHandSideMatrix18;
    rLeftHandSideMatrix(0,8)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix19 + crLeftHandSideMatrix23);
    rLeftHandSideMatrix(0,9)+=crLeftHandSideMatrix24*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(0,10)+=crLeftHandSideMatrix24*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(0,11)+=crLeftHandSideMatrix26;
    rLeftHandSideMatrix(0,12)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix27 + crLeftHandSideMatrix31);
    rLeftHandSideMatrix(0,13)+=crLeftHandSideMatrix32*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(0,14)+=crLeftHandSideMatrix32*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(0,15)+=crLeftHandSideMatrix34;
    rLeftHandSideMatrix(1,0)+=crLeftHandSideMatrix0*crLeftHandSideMatrix36;
    rLeftHandSideMatrix(1,1)+=-w_g*(-crLeftHandSideMatrix35*crLeftHandSideMatrix7 + crLeftHandSideMatrix6);
    rLeftHandSideMatrix(1,2)+=crLeftHandSideMatrix36*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(1,3)+=crLeftHandSideMatrix36;
    rLeftHandSideMatrix(1,4)+=crLeftHandSideMatrix0*crLeftHandSideMatrix38;
    rLeftHandSideMatrix(1,5)+=-w_g*(crLeftHandSideMatrix15 - crLeftHandSideMatrix37*crLeftHandSideMatrix7);
    rLeftHandSideMatrix(1,6)+=crLeftHandSideMatrix38*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(1,7)+=crLeftHandSideMatrix40;
    rLeftHandSideMatrix(1,8)+=crLeftHandSideMatrix0*crLeftHandSideMatrix42;
    rLeftHandSideMatrix(1,9)+=-w_g*(crLeftHandSideMatrix23 - crLeftHandSideMatrix41*crLeftHandSideMatrix7);
    rLeftHandSideMatrix(1,10)+=crLeftHandSideMatrix42*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(1,11)+=crLeftHandSideMatrix44;
    rLeftHandSideMatrix(1,12)+=crLeftHandSideMatrix0*crLeftHandSideMatrix46;
    rLeftHandSideMatrix(1,13)+=-w_g*(crLeftHandSideMatrix31 - crLeftHandSideMatrix45*crLeftHandSideMatrix7);
    rLeftHandSideMatrix(1,14)+=crLeftHandSideMatrix46*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(1,15)+=crLeftHandSideMatrix48;
    rLeftHandSideMatrix(2,0)+=crLeftHandSideMatrix0*crLeftHandSideMatrix50;
    rLeftHandSideMatrix(2,1)+=crLeftHandSideMatrix50*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(2,2)+=-w_g*(-crLeftHandSideMatrix49*crLeftHandSideMatrix9 + crLeftHandSideMatrix6);
    rLeftHandSideMatrix(2,3)+=crLeftHandSideMatrix50;
    rLeftHandSideMatrix(2,4)+=crLeftHandSideMatrix0*crLeftHandSideMatrix52;
    rLeftHandSideMatrix(2,5)+=crLeftHandSideMatrix52*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(2,6)+=-w_g*(crLeftHandSideMatrix15 - crLeftHandSideMatrix51*crLeftHandSideMatrix9);
    rLeftHandSideMatrix(2,7)+=crLeftHandSideMatrix54;
    rLeftHandSideMatrix(2,8)+=crLeftHandSideMatrix0*crLeftHandSideMatrix56;
    rLeftHandSideMatrix(2,9)+=crLeftHandSideMatrix56*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(2,10)+=-w_g*(crLeftHandSideMatrix23 - crLeftHandSideMatrix55*crLeftHandSideMatrix9);
    rLeftHandSideMatrix(2,11)+=crLeftHandSideMatrix58;
    rLeftHandSideMatrix(2,12)+=crLeftHandSideMatrix0*crLeftHandSideMatrix60;
    rLeftHandSideMatrix(2,13)+=crLeftHandSideMatrix60*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(2,14)+=-w_g*(crLeftHandSideMatrix31 - crLeftHandSideMatrix59*crLeftHandSideMatrix9);
    rLeftHandSideMatrix(2,15)+=crLeftHandSideMatrix62;
    rLeftHandSideMatrix(3,0)+=-w_g*(crLeftHandSideMatrix1 - delta_gauss*(-DN(0,0)*crLeftHandSideMatrix65 + crLeftHandSideMatrix4*crLeftHandSideMatrix63 + crLeftHandSideMatrix5*crLeftHandSideMatrix63));
    rLeftHandSideMatrix(3,1)+=-w_g*(crLeftHandSideMatrix35 - delta_gauss*(-DN(0,1)*crLeftHandSideMatrix67 + crLeftHandSideMatrix2*crLeftHandSideMatrix66 + crLeftHandSideMatrix5*crLeftHandSideMatrix66));
    rLeftHandSideMatrix(3,2)+=-w_g*(crLeftHandSideMatrix49 - delta_gauss*(-DN(0,2)*crLeftHandSideMatrix69 + crLeftHandSideMatrix2*crLeftHandSideMatrix68 + crLeftHandSideMatrix4*crLeftHandSideMatrix68));
    rLeftHandSideMatrix(3,3)+=-crLeftHandSideMatrix70*(crLeftHandSideMatrix2 + crLeftHandSideMatrix4 + crLeftHandSideMatrix5);
    rLeftHandSideMatrix(3,4)+=-w_g*(crLeftHandSideMatrix10 - delta_gauss*(-DN(0,0)*crLeftHandSideMatrix72 + crLeftHandSideMatrix73));
    rLeftHandSideMatrix(3,5)+=-w_g*(crLeftHandSideMatrix37 - delta_gauss*(-DN(0,1)*crLeftHandSideMatrix74 + crLeftHandSideMatrix75));
    rLeftHandSideMatrix(3,6)+=-w_g*(crLeftHandSideMatrix51 - delta_gauss*(-DN(0,2)*crLeftHandSideMatrix76 + crLeftHandSideMatrix77));
    rLeftHandSideMatrix(3,7)+=crLeftHandSideMatrix78;
    rLeftHandSideMatrix(3,8)+=-w_g*(crLeftHandSideMatrix19 - delta_gauss*(-DN(0,0)*crLeftHandSideMatrix80 + crLeftHandSideMatrix81));
    rLeftHandSideMatrix(3,9)+=-w_g*(crLeftHandSideMatrix41 - delta_gauss*(-DN(0,1)*crLeftHandSideMatrix82 + crLeftHandSideMatrix83));
    rLeftHandSideMatrix(3,10)+=-w_g*(crLeftHandSideMatrix55 - delta_gauss*(-DN(0,2)*crLeftHandSideMatrix84 + crLeftHandSideMatrix85));
    rLeftHandSideMatrix(3,11)+=crLeftHandSideMatrix86;
    rLeftHandSideMatrix(3,12)+=-w_g*(crLeftHandSideMatrix27 - delta_gauss*(-DN(0,0)*crLeftHandSideMatrix88 + crLeftHandSideMatrix89));
    rLeftHandSideMatrix(3,13)+=-w_g*(crLeftHandSideMatrix45 - delta_gauss*(-DN(0,1)*crLeftHandSideMatrix90 + crLeftHandSideMatrix91));
    rLeftHandSideMatrix(3,14)+=-w_g*(crLeftHandSideMatrix59 - delta_gauss*(-DN(0,2)*crLeftHandSideMatrix92 + crLeftHandSideMatrix93));
    rLeftHandSideMatrix(3,15)+=crLeftHandSideMatrix94;
    rLeftHandSideMatrix(4,0)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix17 + crLeftHandSideMatrix15);
    rLeftHandSideMatrix(4,1)+=crLeftHandSideMatrix18*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(4,2)+=crLeftHandSideMatrix18*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(4,3)+=crLeftHandSideMatrix16;
    rLeftHandSideMatrix(4,4)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix95 + crLeftHandSideMatrix99);
    rLeftHandSideMatrix(4,5)+=crLeftHandSideMatrix100*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(4,6)+=crLeftHandSideMatrix100*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(4,7)+=crLeftHandSideMatrix100;
    rLeftHandSideMatrix(4,8)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix101 + crLeftHandSideMatrix106);
    rLeftHandSideMatrix(4,9)+=crLeftHandSideMatrix107*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(4,10)+=crLeftHandSideMatrix107*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(4,11)+=crLeftHandSideMatrix109;
    rLeftHandSideMatrix(4,12)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix110 + crLeftHandSideMatrix114);
    rLeftHandSideMatrix(4,13)+=crLeftHandSideMatrix115*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(4,14)+=crLeftHandSideMatrix115*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(4,15)+=crLeftHandSideMatrix117;
    rLeftHandSideMatrix(5,0)+=crLeftHandSideMatrix0*crLeftHandSideMatrix40;
    rLeftHandSideMatrix(5,1)+=-w_g*(crLeftHandSideMatrix15 - crLeftHandSideMatrix39*crLeftHandSideMatrix7);
    rLeftHandSideMatrix(5,2)+=crLeftHandSideMatrix40*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(5,3)+=crLeftHandSideMatrix38;
    rLeftHandSideMatrix(5,4)+=crLeftHandSideMatrix0*crLeftHandSideMatrix119;
    rLeftHandSideMatrix(5,5)+=-w_g*(-crLeftHandSideMatrix118*crLeftHandSideMatrix7 + crLeftHandSideMatrix99);
    rLeftHandSideMatrix(5,6)+=crLeftHandSideMatrix119*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(5,7)+=crLeftHandSideMatrix119;
    rLeftHandSideMatrix(5,8)+=crLeftHandSideMatrix0*crLeftHandSideMatrix121;
    rLeftHandSideMatrix(5,9)+=-w_g*(crLeftHandSideMatrix106 - crLeftHandSideMatrix120*crLeftHandSideMatrix7);
    rLeftHandSideMatrix(5,10)+=crLeftHandSideMatrix121*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(5,11)+=crLeftHandSideMatrix123;
    rLeftHandSideMatrix(5,12)+=crLeftHandSideMatrix0*crLeftHandSideMatrix125;
    rLeftHandSideMatrix(5,13)+=-w_g*(crLeftHandSideMatrix114 - crLeftHandSideMatrix124*crLeftHandSideMatrix7);
    rLeftHandSideMatrix(5,14)+=crLeftHandSideMatrix125*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(5,15)+=crLeftHandSideMatrix127;
    rLeftHandSideMatrix(6,0)+=crLeftHandSideMatrix0*crLeftHandSideMatrix54;
    rLeftHandSideMatrix(6,1)+=crLeftHandSideMatrix54*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(6,2)+=-w_g*(crLeftHandSideMatrix15 - crLeftHandSideMatrix53*crLeftHandSideMatrix9);
    rLeftHandSideMatrix(6,3)+=crLeftHandSideMatrix52;
    rLeftHandSideMatrix(6,4)+=crLeftHandSideMatrix0*crLeftHandSideMatrix129;
    rLeftHandSideMatrix(6,5)+=crLeftHandSideMatrix129*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(6,6)+=-w_g*(-crLeftHandSideMatrix128*crLeftHandSideMatrix9 + crLeftHandSideMatrix99);
    rLeftHandSideMatrix(6,7)+=crLeftHandSideMatrix129;
    rLeftHandSideMatrix(6,8)+=crLeftHandSideMatrix0*crLeftHandSideMatrix131;
    rLeftHandSideMatrix(6,9)+=crLeftHandSideMatrix131*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(6,10)+=-w_g*(crLeftHandSideMatrix106 - crLeftHandSideMatrix130*crLeftHandSideMatrix9);
    rLeftHandSideMatrix(6,11)+=crLeftHandSideMatrix133;
    rLeftHandSideMatrix(6,12)+=crLeftHandSideMatrix0*crLeftHandSideMatrix135;
    rLeftHandSideMatrix(6,13)+=crLeftHandSideMatrix135*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(6,14)+=-w_g*(crLeftHandSideMatrix114 - crLeftHandSideMatrix134*crLeftHandSideMatrix9);
    rLeftHandSideMatrix(6,15)+=crLeftHandSideMatrix137;
    rLeftHandSideMatrix(7,0)+=-w_g*(crLeftHandSideMatrix17 - delta_gauss*(-DN(1,0)*crLeftHandSideMatrix65 + crLeftHandSideMatrix73));
    rLeftHandSideMatrix(7,1)+=-w_g*(crLeftHandSideMatrix39 - delta_gauss*(-DN(1,1)*crLeftHandSideMatrix67 + crLeftHandSideMatrix75));
    rLeftHandSideMatrix(7,2)+=-w_g*(crLeftHandSideMatrix53 - delta_gauss*(-DN(1,2)*crLeftHandSideMatrix69 + crLeftHandSideMatrix77));
    rLeftHandSideMatrix(7,3)+=crLeftHandSideMatrix78;
    rLeftHandSideMatrix(7,4)+=-w_g*(crLeftHandSideMatrix95 - delta_gauss*(-DN(1,0)*crLeftHandSideMatrix72 + crLeftHandSideMatrix63*crLeftHandSideMatrix97 + crLeftHandSideMatrix63*crLeftHandSideMatrix98));
    rLeftHandSideMatrix(7,5)+=-w_g*(crLeftHandSideMatrix118 - delta_gauss*(-DN(1,1)*crLeftHandSideMatrix74 + crLeftHandSideMatrix66*crLeftHandSideMatrix96 + crLeftHandSideMatrix66*crLeftHandSideMatrix98));
    rLeftHandSideMatrix(7,6)+=-w_g*(crLeftHandSideMatrix128 - delta_gauss*(-DN(1,2)*crLeftHandSideMatrix76 + crLeftHandSideMatrix68*crLeftHandSideMatrix96 + crLeftHandSideMatrix68*crLeftHandSideMatrix97));
    rLeftHandSideMatrix(7,7)+=-crLeftHandSideMatrix70*(crLeftHandSideMatrix96 + crLeftHandSideMatrix97 + crLeftHandSideMatrix98);
    rLeftHandSideMatrix(7,8)+=-w_g*(crLeftHandSideMatrix101 - delta_gauss*(-DN(1,0)*crLeftHandSideMatrix80 + crLeftHandSideMatrix138));
    rLeftHandSideMatrix(7,9)+=-w_g*(crLeftHandSideMatrix120 - delta_gauss*(-DN(1,1)*crLeftHandSideMatrix82 + crLeftHandSideMatrix139));
    rLeftHandSideMatrix(7,10)+=-w_g*(crLeftHandSideMatrix130 - delta_gauss*(-DN(1,2)*crLeftHandSideMatrix84 + crLeftHandSideMatrix140));
    rLeftHandSideMatrix(7,11)+=crLeftHandSideMatrix141;
    rLeftHandSideMatrix(7,12)+=-w_g*(crLeftHandSideMatrix110 - delta_gauss*(-DN(1,0)*crLeftHandSideMatrix88 + crLeftHandSideMatrix142));
    rLeftHandSideMatrix(7,13)+=-w_g*(crLeftHandSideMatrix124 - delta_gauss*(-DN(1,1)*crLeftHandSideMatrix90 + crLeftHandSideMatrix143));
    rLeftHandSideMatrix(7,14)+=-w_g*(crLeftHandSideMatrix134 - delta_gauss*(-DN(1,2)*crLeftHandSideMatrix92 + crLeftHandSideMatrix144));
    rLeftHandSideMatrix(7,15)+=crLeftHandSideMatrix145;
    rLeftHandSideMatrix(8,0)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix25 + crLeftHandSideMatrix23);
    rLeftHandSideMatrix(8,1)+=crLeftHandSideMatrix26*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(8,2)+=crLeftHandSideMatrix26*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(8,3)+=crLeftHandSideMatrix24;
    rLeftHandSideMatrix(8,4)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix108 + crLeftHandSideMatrix106);
    rLeftHandSideMatrix(8,5)+=crLeftHandSideMatrix109*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(8,6)+=crLeftHandSideMatrix109*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(8,7)+=crLeftHandSideMatrix107;
    rLeftHandSideMatrix(8,8)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix146 + crLeftHandSideMatrix150);
    rLeftHandSideMatrix(8,9)+=crLeftHandSideMatrix151*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(8,10)+=crLeftHandSideMatrix151*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(8,11)+=crLeftHandSideMatrix151;
    rLeftHandSideMatrix(8,12)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix152 + crLeftHandSideMatrix156);
    rLeftHandSideMatrix(8,13)+=crLeftHandSideMatrix157*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(8,14)+=crLeftHandSideMatrix157*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(8,15)+=crLeftHandSideMatrix159;
    rLeftHandSideMatrix(9,0)+=crLeftHandSideMatrix0*crLeftHandSideMatrix44;
    rLeftHandSideMatrix(9,1)+=-w_g*(crLeftHandSideMatrix23 - crLeftHandSideMatrix43*crLeftHandSideMatrix7);
    rLeftHandSideMatrix(9,2)+=crLeftHandSideMatrix44*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(9,3)+=crLeftHandSideMatrix42;
    rLeftHandSideMatrix(9,4)+=crLeftHandSideMatrix0*crLeftHandSideMatrix123;
    rLeftHandSideMatrix(9,5)+=-w_g*(crLeftHandSideMatrix106 - crLeftHandSideMatrix122*crLeftHandSideMatrix7);
    rLeftHandSideMatrix(9,6)+=crLeftHandSideMatrix123*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(9,7)+=crLeftHandSideMatrix121;
    rLeftHandSideMatrix(9,8)+=crLeftHandSideMatrix0*crLeftHandSideMatrix161;
    rLeftHandSideMatrix(9,9)+=-w_g*(crLeftHandSideMatrix150 - crLeftHandSideMatrix160*crLeftHandSideMatrix7);
    rLeftHandSideMatrix(9,10)+=crLeftHandSideMatrix161*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(9,11)+=crLeftHandSideMatrix161;
    rLeftHandSideMatrix(9,12)+=crLeftHandSideMatrix0*crLeftHandSideMatrix163;
    rLeftHandSideMatrix(9,13)+=-w_g*(crLeftHandSideMatrix156 - crLeftHandSideMatrix162*crLeftHandSideMatrix7);
    rLeftHandSideMatrix(9,14)+=crLeftHandSideMatrix163*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(9,15)+=crLeftHandSideMatrix165;
    rLeftHandSideMatrix(10,0)+=crLeftHandSideMatrix0*crLeftHandSideMatrix58;
    rLeftHandSideMatrix(10,1)+=crLeftHandSideMatrix58*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(10,2)+=-w_g*(crLeftHandSideMatrix23 - crLeftHandSideMatrix57*crLeftHandSideMatrix9);
    rLeftHandSideMatrix(10,3)+=crLeftHandSideMatrix56;
    rLeftHandSideMatrix(10,4)+=crLeftHandSideMatrix0*crLeftHandSideMatrix133;
    rLeftHandSideMatrix(10,5)+=crLeftHandSideMatrix133*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(10,6)+=-w_g*(crLeftHandSideMatrix106 - crLeftHandSideMatrix132*crLeftHandSideMatrix9);
    rLeftHandSideMatrix(10,7)+=crLeftHandSideMatrix131;
    rLeftHandSideMatrix(10,8)+=crLeftHandSideMatrix0*crLeftHandSideMatrix167;
    rLeftHandSideMatrix(10,9)+=crLeftHandSideMatrix167*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(10,10)+=-w_g*(crLeftHandSideMatrix150 - crLeftHandSideMatrix166*crLeftHandSideMatrix9);
    rLeftHandSideMatrix(10,11)+=crLeftHandSideMatrix167;
    rLeftHandSideMatrix(10,12)+=crLeftHandSideMatrix0*crLeftHandSideMatrix169;
    rLeftHandSideMatrix(10,13)+=crLeftHandSideMatrix169*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(10,14)+=-w_g*(crLeftHandSideMatrix156 - crLeftHandSideMatrix168*crLeftHandSideMatrix9);
    rLeftHandSideMatrix(10,15)+=crLeftHandSideMatrix171;
    rLeftHandSideMatrix(11,0)+=-w_g*(crLeftHandSideMatrix25 - delta_gauss*(-DN(2,0)*crLeftHandSideMatrix65 + crLeftHandSideMatrix81));
    rLeftHandSideMatrix(11,1)+=-w_g*(crLeftHandSideMatrix43 - delta_gauss*(-DN(2,1)*crLeftHandSideMatrix67 + crLeftHandSideMatrix83));
    rLeftHandSideMatrix(11,2)+=-w_g*(crLeftHandSideMatrix57 - delta_gauss*(-DN(2,2)*crLeftHandSideMatrix69 + crLeftHandSideMatrix85));
    rLeftHandSideMatrix(11,3)+=crLeftHandSideMatrix86;
    rLeftHandSideMatrix(11,4)+=-w_g*(crLeftHandSideMatrix108 - delta_gauss*(-DN(2,0)*crLeftHandSideMatrix72 + crLeftHandSideMatrix138));
    rLeftHandSideMatrix(11,5)+=-w_g*(crLeftHandSideMatrix122 - delta_gauss*(-DN(2,1)*crLeftHandSideMatrix74 + crLeftHandSideMatrix139));
    rLeftHandSideMatrix(11,6)+=-w_g*(crLeftHandSideMatrix132 - delta_gauss*(-DN(2,2)*crLeftHandSideMatrix76 + crLeftHandSideMatrix140));
    rLeftHandSideMatrix(11,7)+=crLeftHandSideMatrix141;
    rLeftHandSideMatrix(11,8)+=-w_g*(crLeftHandSideMatrix146 - delta_gauss*(-DN(2,0)*crLeftHandSideMatrix80 + crLeftHandSideMatrix148*crLeftHandSideMatrix63 + crLeftHandSideMatrix149*crLeftHandSideMatrix63));
    rLeftHandSideMatrix(11,9)+=-w_g*(crLeftHandSideMatrix160 - delta_gauss*(-DN(2,1)*crLeftHandSideMatrix82 + crLeftHandSideMatrix147*crLeftHandSideMatrix66 + crLeftHandSideMatrix149*crLeftHandSideMatrix66));
    rLeftHandSideMatrix(11,10)+=-w_g*(crLeftHandSideMatrix166 - delta_gauss*(-DN(2,2)*crLeftHandSideMatrix84 + crLeftHandSideMatrix147*crLeftHandSideMatrix68 + crLeftHandSideMatrix148*crLeftHandSideMatrix68));
    rLeftHandSideMatrix(11,11)+=-crLeftHandSideMatrix70*(crLeftHandSideMatrix147 + crLeftHandSideMatrix148 + crLeftHandSideMatrix149);
    rLeftHandSideMatrix(11,12)+=-w_g*(crLeftHandSideMatrix152 - delta_gauss*(-DN(2,0)*crLeftHandSideMatrix88 + crLeftHandSideMatrix172));
    rLeftHandSideMatrix(11,13)+=-w_g*(crLeftHandSideMatrix162 - delta_gauss*(-DN(2,1)*crLeftHandSideMatrix90 + crLeftHandSideMatrix173));
    rLeftHandSideMatrix(11,14)+=-w_g*(crLeftHandSideMatrix168 - delta_gauss*(-DN(2,2)*crLeftHandSideMatrix92 + crLeftHandSideMatrix174));
    rLeftHandSideMatrix(11,15)+=crLeftHandSideMatrix175;
    rLeftHandSideMatrix(12,0)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix33 + crLeftHandSideMatrix31);
    rLeftHandSideMatrix(12,1)+=crLeftHandSideMatrix34*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(12,2)+=crLeftHandSideMatrix34*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(12,3)+=crLeftHandSideMatrix32;
    rLeftHandSideMatrix(12,4)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix116 + crLeftHandSideMatrix114);
    rLeftHandSideMatrix(12,5)+=crLeftHandSideMatrix117*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(12,6)+=crLeftHandSideMatrix117*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(12,7)+=crLeftHandSideMatrix115;
    rLeftHandSideMatrix(12,8)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix158 + crLeftHandSideMatrix156);
    rLeftHandSideMatrix(12,9)+=crLeftHandSideMatrix159*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(12,10)+=crLeftHandSideMatrix159*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(12,11)+=crLeftHandSideMatrix157;
    rLeftHandSideMatrix(12,12)+=-w_g*(-crLeftHandSideMatrix0*crLeftHandSideMatrix176 + crLeftHandSideMatrix180);
    rLeftHandSideMatrix(12,13)+=crLeftHandSideMatrix181*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(12,14)+=crLeftHandSideMatrix181*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(12,15)+=crLeftHandSideMatrix181;
    rLeftHandSideMatrix(13,0)+=crLeftHandSideMatrix0*crLeftHandSideMatrix48;
    rLeftHandSideMatrix(13,1)+=-w_g*(crLeftHandSideMatrix31 - crLeftHandSideMatrix47*crLeftHandSideMatrix7);
    rLeftHandSideMatrix(13,2)+=crLeftHandSideMatrix48*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(13,3)+=crLeftHandSideMatrix46;
    rLeftHandSideMatrix(13,4)+=crLeftHandSideMatrix0*crLeftHandSideMatrix127;
    rLeftHandSideMatrix(13,5)+=-w_g*(crLeftHandSideMatrix114 - crLeftHandSideMatrix126*crLeftHandSideMatrix7);
    rLeftHandSideMatrix(13,6)+=crLeftHandSideMatrix127*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(13,7)+=crLeftHandSideMatrix125;
    rLeftHandSideMatrix(13,8)+=crLeftHandSideMatrix0*crLeftHandSideMatrix165;
    rLeftHandSideMatrix(13,9)+=-w_g*(crLeftHandSideMatrix156 - crLeftHandSideMatrix164*crLeftHandSideMatrix7);
    rLeftHandSideMatrix(13,10)+=crLeftHandSideMatrix165*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(13,11)+=crLeftHandSideMatrix163;
    rLeftHandSideMatrix(13,12)+=crLeftHandSideMatrix0*crLeftHandSideMatrix183;
    rLeftHandSideMatrix(13,13)+=-w_g*(crLeftHandSideMatrix180 - crLeftHandSideMatrix182*crLeftHandSideMatrix7);
    rLeftHandSideMatrix(13,14)+=crLeftHandSideMatrix183*crLeftHandSideMatrix9;
    rLeftHandSideMatrix(13,15)+=crLeftHandSideMatrix183;
    rLeftHandSideMatrix(14,0)+=crLeftHandSideMatrix0*crLeftHandSideMatrix62;
    rLeftHandSideMatrix(14,1)+=crLeftHandSideMatrix62*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(14,2)+=-w_g*(crLeftHandSideMatrix31 - crLeftHandSideMatrix61*crLeftHandSideMatrix9);
    rLeftHandSideMatrix(14,3)+=crLeftHandSideMatrix60;
    rLeftHandSideMatrix(14,4)+=crLeftHandSideMatrix0*crLeftHandSideMatrix137;
    rLeftHandSideMatrix(14,5)+=crLeftHandSideMatrix137*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(14,6)+=-w_g*(crLeftHandSideMatrix114 - crLeftHandSideMatrix136*crLeftHandSideMatrix9);
    rLeftHandSideMatrix(14,7)+=crLeftHandSideMatrix135;
    rLeftHandSideMatrix(14,8)+=crLeftHandSideMatrix0*crLeftHandSideMatrix171;
    rLeftHandSideMatrix(14,9)+=crLeftHandSideMatrix171*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(14,10)+=-w_g*(crLeftHandSideMatrix156 - crLeftHandSideMatrix170*crLeftHandSideMatrix9);
    rLeftHandSideMatrix(14,11)+=crLeftHandSideMatrix169;
    rLeftHandSideMatrix(14,12)+=crLeftHandSideMatrix0*crLeftHandSideMatrix185;
    rLeftHandSideMatrix(14,13)+=crLeftHandSideMatrix185*crLeftHandSideMatrix7;
    rLeftHandSideMatrix(14,14)+=-w_g*(crLeftHandSideMatrix180 - crLeftHandSideMatrix184*crLeftHandSideMatrix9);
    rLeftHandSideMatrix(14,15)+=crLeftHandSideMatrix185;
    rLeftHandSideMatrix(15,0)+=-w_g*(crLeftHandSideMatrix33 - delta_gauss*(-DN(3,0)*crLeftHandSideMatrix65 + crLeftHandSideMatrix89));
    rLeftHandSideMatrix(15,1)+=-w_g*(crLeftHandSideMatrix47 - delta_gauss*(-DN(3,1)*crLeftHandSideMatrix67 + crLeftHandSideMatrix91));
    rLeftHandSideMatrix(15,2)+=-w_g*(crLeftHandSideMatrix61 - delta_gauss*(-DN(3,2)*crLeftHandSideMatrix69 + crLeftHandSideMatrix93));
    rLeftHandSideMatrix(15,3)+=crLeftHandSideMatrix94;
    rLeftHandSideMatrix(15,4)+=-w_g*(crLeftHandSideMatrix116 - delta_gauss*(-DN(3,0)*crLeftHandSideMatrix72 + crLeftHandSideMatrix142));
    rLeftHandSideMatrix(15,5)+=-w_g*(crLeftHandSideMatrix126 - delta_gauss*(-DN(3,1)*crLeftHandSideMatrix74 + crLeftHandSideMatrix143));
    rLeftHandSideMatrix(15,6)+=-w_g*(crLeftHandSideMatrix136 - delta_gauss*(-DN(3,2)*crLeftHandSideMatrix76 + crLeftHandSideMatrix144));
    rLeftHandSideMatrix(15,7)+=crLeftHandSideMatrix145;
    rLeftHandSideMatrix(15,8)+=-w_g*(crLeftHandSideMatrix158 - delta_gauss*(-DN(3,0)*crLeftHandSideMatrix80 + crLeftHandSideMatrix172));
    rLeftHandSideMatrix(15,9)+=-w_g*(crLeftHandSideMatrix164 - delta_gauss*(-DN(3,1)*crLeftHandSideMatrix82 + crLeftHandSideMatrix173));
    rLeftHandSideMatrix(15,10)+=-w_g*(crLeftHandSideMatrix170 - delta_gauss*(-DN(3,2)*crLeftHandSideMatrix84 + crLeftHandSideMatrix174));
    rLeftHandSideMatrix(15,11)+=crLeftHandSideMatrix175;
    rLeftHandSideMatrix(15,12)+=-w_g*(crLeftHandSideMatrix176 - delta_gauss*(-DN(3,0)*crLeftHandSideMatrix88 + crLeftHandSideMatrix178*crLeftHandSideMatrix63 + crLeftHandSideMatrix179*crLeftHandSideMatrix63));
    rLeftHandSideMatrix(15,13)+=-w_g*(crLeftHandSideMatrix182 - delta_gauss*(-DN(3,1)*crLeftHandSideMatrix90 + crLeftHandSideMatrix177*crLeftHandSideMatrix66 + crLeftHandSideMatrix179*crLeftHandSideMatrix66));
    rLeftHandSideMatrix(15,14)+=-w_g*(crLeftHandSideMatrix184 - delta_gauss*(-DN(3,2)*crLeftHandSideMatrix92 + crLeftHandSideMatrix177*crLeftHandSideMatrix68 + crLeftHandSideMatrix178*crLeftHandSideMatrix68));
    rLeftHandSideMatrix(15,15)+=-crLeftHandSideMatrix70*(crLeftHandSideMatrix177 + crLeftHandSideMatrix178 + crLeftHandSideMatrix179);
    

}

template <>
void FirstOrderStokesVariableViscosityBvsGl<3>::AddConditionGaussPointLeftHandSideContribution(
    const GlobalPointer<Condition>& rConditionPointer,
    const ElementDataContainer& rData,
    MatrixType& rLHS)
{
    KRATOS_ERROR << "Method AddConditionGaussPointLeftHandSideContribution for element type FirstOrderStokesVariableViscosityBvsGl is not implemented for 3D yet" << std::endl;


    const auto& u_nodes = rData.Velocity;
    const auto& nu_nodes = rData.DynamicViscosity;

    const auto& N = rData.N;
    const auto& DN = rData.DN;

    // Assemble LHS contribution
    const double w_g = rData.Weight;

    // Get stabilization data
    const double delta_gauss = rData.Delta;

    const Vector n_gauss = ZeroVector(3);

    Matrix rLeftHandSideMatrix = ZeroMatrix(LocalSize);

        const double crLeftHandSideMatrix0 = DN(0,0)*n_gauss[1] - DN(0,1)*n_gauss[0];
    const double crLeftHandSideMatrix1 = DN(0,0)*n_gauss[2] - DN(0,2)*n_gauss[0];
    const double crLeftHandSideMatrix2 = delta_gauss*w_g*(N[0]*nu_nodes[0] + N[1]*nu_nodes[1] + N[2]*nu_nodes[2] + N[3]*nu_nodes[3]);
    const double crLeftHandSideMatrix3 = DN(0,1)*n_gauss[2] - DN(0,2)*n_gauss[1];
    const double crLeftHandSideMatrix4 = DN(1,0)*n_gauss[1] - DN(1,1)*n_gauss[0];
    const double crLeftHandSideMatrix5 = DN(1,0)*n_gauss[2] - DN(1,2)*n_gauss[0];
    const double crLeftHandSideMatrix6 = DN(1,1)*n_gauss[2] - DN(1,2)*n_gauss[1];
    const double crLeftHandSideMatrix7 = DN(2,0)*n_gauss[1] - DN(2,1)*n_gauss[0];
    const double crLeftHandSideMatrix8 = DN(2,0)*n_gauss[2] - DN(2,2)*n_gauss[0];
    const double crLeftHandSideMatrix9 = DN(2,1)*n_gauss[2] - DN(2,2)*n_gauss[1];
    const double crLeftHandSideMatrix10 = DN(3,0)*n_gauss[1] - DN(3,1)*n_gauss[0];
    const double crLeftHandSideMatrix11 = DN(3,0)*n_gauss[2] - DN(3,2)*n_gauss[0];
    const double crLeftHandSideMatrix12 = DN(3,1)*n_gauss[2] - DN(3,2)*n_gauss[1];
    rLeftHandSideMatrix(0,0)+=0;
    rLeftHandSideMatrix(0,1)+=0;
    rLeftHandSideMatrix(0,2)+=0;
    rLeftHandSideMatrix(0,3)+=0;
    rLeftHandSideMatrix(0,4)+=0;
    rLeftHandSideMatrix(0,5)+=0;
    rLeftHandSideMatrix(0,6)+=0;
    rLeftHandSideMatrix(0,7)+=0;
    rLeftHandSideMatrix(0,8)+=0;
    rLeftHandSideMatrix(0,9)+=0;
    rLeftHandSideMatrix(0,10)+=0;
    rLeftHandSideMatrix(0,11)+=0;
    rLeftHandSideMatrix(0,12)+=0;
    rLeftHandSideMatrix(0,13)+=0;
    rLeftHandSideMatrix(0,14)+=0;
    rLeftHandSideMatrix(0,15)+=0;
    rLeftHandSideMatrix(1,0)+=0;
    rLeftHandSideMatrix(1,1)+=0;
    rLeftHandSideMatrix(1,2)+=0;
    rLeftHandSideMatrix(1,3)+=0;
    rLeftHandSideMatrix(1,4)+=0;
    rLeftHandSideMatrix(1,5)+=0;
    rLeftHandSideMatrix(1,6)+=0;
    rLeftHandSideMatrix(1,7)+=0;
    rLeftHandSideMatrix(1,8)+=0;
    rLeftHandSideMatrix(1,9)+=0;
    rLeftHandSideMatrix(1,10)+=0;
    rLeftHandSideMatrix(1,11)+=0;
    rLeftHandSideMatrix(1,12)+=0;
    rLeftHandSideMatrix(1,13)+=0;
    rLeftHandSideMatrix(1,14)+=0;
    rLeftHandSideMatrix(1,15)+=0;
    rLeftHandSideMatrix(2,0)+=0;
    rLeftHandSideMatrix(2,1)+=0;
    rLeftHandSideMatrix(2,2)+=0;
    rLeftHandSideMatrix(2,3)+=0;
    rLeftHandSideMatrix(2,4)+=0;
    rLeftHandSideMatrix(2,5)+=0;
    rLeftHandSideMatrix(2,6)+=0;
    rLeftHandSideMatrix(2,7)+=0;
    rLeftHandSideMatrix(2,8)+=0;
    rLeftHandSideMatrix(2,9)+=0;
    rLeftHandSideMatrix(2,10)+=0;
    rLeftHandSideMatrix(2,11)+=0;
    rLeftHandSideMatrix(2,12)+=0;
    rLeftHandSideMatrix(2,13)+=0;
    rLeftHandSideMatrix(2,14)+=0;
    rLeftHandSideMatrix(2,15)+=0;
    rLeftHandSideMatrix(3,0)+=crLeftHandSideMatrix2*(DN(0,1)*crLeftHandSideMatrix0 + DN(0,2)*crLeftHandSideMatrix1);
    rLeftHandSideMatrix(3,1)+=-crLeftHandSideMatrix2*(DN(0,0)*crLeftHandSideMatrix0 - DN(0,2)*crLeftHandSideMatrix3);
    rLeftHandSideMatrix(3,2)+=-crLeftHandSideMatrix2*(DN(0,0)*crLeftHandSideMatrix1 + DN(0,1)*crLeftHandSideMatrix3);
    rLeftHandSideMatrix(3,3)+=0;
    rLeftHandSideMatrix(3,4)+=crLeftHandSideMatrix2*(DN(1,1)*crLeftHandSideMatrix0 + DN(1,2)*crLeftHandSideMatrix1);
    rLeftHandSideMatrix(3,5)+=-crLeftHandSideMatrix2*(DN(1,0)*crLeftHandSideMatrix0 - DN(1,2)*crLeftHandSideMatrix3);
    rLeftHandSideMatrix(3,6)+=-crLeftHandSideMatrix2*(DN(1,0)*crLeftHandSideMatrix1 + DN(1,1)*crLeftHandSideMatrix3);
    rLeftHandSideMatrix(3,7)+=0;
    rLeftHandSideMatrix(3,8)+=crLeftHandSideMatrix2*(DN(2,1)*crLeftHandSideMatrix0 + DN(2,2)*crLeftHandSideMatrix1);
    rLeftHandSideMatrix(3,9)+=-crLeftHandSideMatrix2*(DN(2,0)*crLeftHandSideMatrix0 - DN(2,2)*crLeftHandSideMatrix3);
    rLeftHandSideMatrix(3,10)+=-crLeftHandSideMatrix2*(DN(2,0)*crLeftHandSideMatrix1 + DN(2,1)*crLeftHandSideMatrix3);
    rLeftHandSideMatrix(3,11)+=0;
    rLeftHandSideMatrix(3,12)+=crLeftHandSideMatrix2*(DN(3,1)*crLeftHandSideMatrix0 + DN(3,2)*crLeftHandSideMatrix1);
    rLeftHandSideMatrix(3,13)+=-crLeftHandSideMatrix2*(DN(3,0)*crLeftHandSideMatrix0 - DN(3,2)*crLeftHandSideMatrix3);
    rLeftHandSideMatrix(3,14)+=-crLeftHandSideMatrix2*(DN(3,0)*crLeftHandSideMatrix1 + DN(3,1)*crLeftHandSideMatrix3);
    rLeftHandSideMatrix(3,15)+=0;
    rLeftHandSideMatrix(4,0)+=0;
    rLeftHandSideMatrix(4,1)+=0;
    rLeftHandSideMatrix(4,2)+=0;
    rLeftHandSideMatrix(4,3)+=0;
    rLeftHandSideMatrix(4,4)+=0;
    rLeftHandSideMatrix(4,5)+=0;
    rLeftHandSideMatrix(4,6)+=0;
    rLeftHandSideMatrix(4,7)+=0;
    rLeftHandSideMatrix(4,8)+=0;
    rLeftHandSideMatrix(4,9)+=0;
    rLeftHandSideMatrix(4,10)+=0;
    rLeftHandSideMatrix(4,11)+=0;
    rLeftHandSideMatrix(4,12)+=0;
    rLeftHandSideMatrix(4,13)+=0;
    rLeftHandSideMatrix(4,14)+=0;
    rLeftHandSideMatrix(4,15)+=0;
    rLeftHandSideMatrix(5,0)+=0;
    rLeftHandSideMatrix(5,1)+=0;
    rLeftHandSideMatrix(5,2)+=0;
    rLeftHandSideMatrix(5,3)+=0;
    rLeftHandSideMatrix(5,4)+=0;
    rLeftHandSideMatrix(5,5)+=0;
    rLeftHandSideMatrix(5,6)+=0;
    rLeftHandSideMatrix(5,7)+=0;
    rLeftHandSideMatrix(5,8)+=0;
    rLeftHandSideMatrix(5,9)+=0;
    rLeftHandSideMatrix(5,10)+=0;
    rLeftHandSideMatrix(5,11)+=0;
    rLeftHandSideMatrix(5,12)+=0;
    rLeftHandSideMatrix(5,13)+=0;
    rLeftHandSideMatrix(5,14)+=0;
    rLeftHandSideMatrix(5,15)+=0;
    rLeftHandSideMatrix(6,0)+=0;
    rLeftHandSideMatrix(6,1)+=0;
    rLeftHandSideMatrix(6,2)+=0;
    rLeftHandSideMatrix(6,3)+=0;
    rLeftHandSideMatrix(6,4)+=0;
    rLeftHandSideMatrix(6,5)+=0;
    rLeftHandSideMatrix(6,6)+=0;
    rLeftHandSideMatrix(6,7)+=0;
    rLeftHandSideMatrix(6,8)+=0;
    rLeftHandSideMatrix(6,9)+=0;
    rLeftHandSideMatrix(6,10)+=0;
    rLeftHandSideMatrix(6,11)+=0;
    rLeftHandSideMatrix(6,12)+=0;
    rLeftHandSideMatrix(6,13)+=0;
    rLeftHandSideMatrix(6,14)+=0;
    rLeftHandSideMatrix(6,15)+=0;
    rLeftHandSideMatrix(7,0)+=crLeftHandSideMatrix2*(DN(0,1)*crLeftHandSideMatrix4 + DN(0,2)*crLeftHandSideMatrix5);
    rLeftHandSideMatrix(7,1)+=-crLeftHandSideMatrix2*(DN(0,0)*crLeftHandSideMatrix4 - DN(0,2)*crLeftHandSideMatrix6);
    rLeftHandSideMatrix(7,2)+=-crLeftHandSideMatrix2*(DN(0,0)*crLeftHandSideMatrix5 + DN(0,1)*crLeftHandSideMatrix6);
    rLeftHandSideMatrix(7,3)+=0;
    rLeftHandSideMatrix(7,4)+=crLeftHandSideMatrix2*(DN(1,1)*crLeftHandSideMatrix4 + DN(1,2)*crLeftHandSideMatrix5);
    rLeftHandSideMatrix(7,5)+=-crLeftHandSideMatrix2*(DN(1,0)*crLeftHandSideMatrix4 - DN(1,2)*crLeftHandSideMatrix6);
    rLeftHandSideMatrix(7,6)+=-crLeftHandSideMatrix2*(DN(1,0)*crLeftHandSideMatrix5 + DN(1,1)*crLeftHandSideMatrix6);
    rLeftHandSideMatrix(7,7)+=0;
    rLeftHandSideMatrix(7,8)+=crLeftHandSideMatrix2*(DN(2,1)*crLeftHandSideMatrix4 + DN(2,2)*crLeftHandSideMatrix5);
    rLeftHandSideMatrix(7,9)+=-crLeftHandSideMatrix2*(DN(2,0)*crLeftHandSideMatrix4 - DN(2,2)*crLeftHandSideMatrix6);
    rLeftHandSideMatrix(7,10)+=-crLeftHandSideMatrix2*(DN(2,0)*crLeftHandSideMatrix5 + DN(2,1)*crLeftHandSideMatrix6);
    rLeftHandSideMatrix(7,11)+=0;
    rLeftHandSideMatrix(7,12)+=crLeftHandSideMatrix2*(DN(3,1)*crLeftHandSideMatrix4 + DN(3,2)*crLeftHandSideMatrix5);
    rLeftHandSideMatrix(7,13)+=-crLeftHandSideMatrix2*(DN(3,0)*crLeftHandSideMatrix4 - DN(3,2)*crLeftHandSideMatrix6);
    rLeftHandSideMatrix(7,14)+=-crLeftHandSideMatrix2*(DN(3,0)*crLeftHandSideMatrix5 + DN(3,1)*crLeftHandSideMatrix6);
    rLeftHandSideMatrix(7,15)+=0;
    rLeftHandSideMatrix(8,0)+=0;
    rLeftHandSideMatrix(8,1)+=0;
    rLeftHandSideMatrix(8,2)+=0;
    rLeftHandSideMatrix(8,3)+=0;
    rLeftHandSideMatrix(8,4)+=0;
    rLeftHandSideMatrix(8,5)+=0;
    rLeftHandSideMatrix(8,6)+=0;
    rLeftHandSideMatrix(8,7)+=0;
    rLeftHandSideMatrix(8,8)+=0;
    rLeftHandSideMatrix(8,9)+=0;
    rLeftHandSideMatrix(8,10)+=0;
    rLeftHandSideMatrix(8,11)+=0;
    rLeftHandSideMatrix(8,12)+=0;
    rLeftHandSideMatrix(8,13)+=0;
    rLeftHandSideMatrix(8,14)+=0;
    rLeftHandSideMatrix(8,15)+=0;
    rLeftHandSideMatrix(9,0)+=0;
    rLeftHandSideMatrix(9,1)+=0;
    rLeftHandSideMatrix(9,2)+=0;
    rLeftHandSideMatrix(9,3)+=0;
    rLeftHandSideMatrix(9,4)+=0;
    rLeftHandSideMatrix(9,5)+=0;
    rLeftHandSideMatrix(9,6)+=0;
    rLeftHandSideMatrix(9,7)+=0;
    rLeftHandSideMatrix(9,8)+=0;
    rLeftHandSideMatrix(9,9)+=0;
    rLeftHandSideMatrix(9,10)+=0;
    rLeftHandSideMatrix(9,11)+=0;
    rLeftHandSideMatrix(9,12)+=0;
    rLeftHandSideMatrix(9,13)+=0;
    rLeftHandSideMatrix(9,14)+=0;
    rLeftHandSideMatrix(9,15)+=0;
    rLeftHandSideMatrix(10,0)+=0;
    rLeftHandSideMatrix(10,1)+=0;
    rLeftHandSideMatrix(10,2)+=0;
    rLeftHandSideMatrix(10,3)+=0;
    rLeftHandSideMatrix(10,4)+=0;
    rLeftHandSideMatrix(10,5)+=0;
    rLeftHandSideMatrix(10,6)+=0;
    rLeftHandSideMatrix(10,7)+=0;
    rLeftHandSideMatrix(10,8)+=0;
    rLeftHandSideMatrix(10,9)+=0;
    rLeftHandSideMatrix(10,10)+=0;
    rLeftHandSideMatrix(10,11)+=0;
    rLeftHandSideMatrix(10,12)+=0;
    rLeftHandSideMatrix(10,13)+=0;
    rLeftHandSideMatrix(10,14)+=0;
    rLeftHandSideMatrix(10,15)+=0;
    rLeftHandSideMatrix(11,0)+=crLeftHandSideMatrix2*(DN(0,1)*crLeftHandSideMatrix7 + DN(0,2)*crLeftHandSideMatrix8);
    rLeftHandSideMatrix(11,1)+=-crLeftHandSideMatrix2*(DN(0,0)*crLeftHandSideMatrix7 - DN(0,2)*crLeftHandSideMatrix9);
    rLeftHandSideMatrix(11,2)+=-crLeftHandSideMatrix2*(DN(0,0)*crLeftHandSideMatrix8 + DN(0,1)*crLeftHandSideMatrix9);
    rLeftHandSideMatrix(11,3)+=0;
    rLeftHandSideMatrix(11,4)+=crLeftHandSideMatrix2*(DN(1,1)*crLeftHandSideMatrix7 + DN(1,2)*crLeftHandSideMatrix8);
    rLeftHandSideMatrix(11,5)+=-crLeftHandSideMatrix2*(DN(1,0)*crLeftHandSideMatrix7 - DN(1,2)*crLeftHandSideMatrix9);
    rLeftHandSideMatrix(11,6)+=-crLeftHandSideMatrix2*(DN(1,0)*crLeftHandSideMatrix8 + DN(1,1)*crLeftHandSideMatrix9);
    rLeftHandSideMatrix(11,7)+=0;
    rLeftHandSideMatrix(11,8)+=crLeftHandSideMatrix2*(DN(2,1)*crLeftHandSideMatrix7 + DN(2,2)*crLeftHandSideMatrix8);
    rLeftHandSideMatrix(11,9)+=-crLeftHandSideMatrix2*(DN(2,0)*crLeftHandSideMatrix7 - DN(2,2)*crLeftHandSideMatrix9);
    rLeftHandSideMatrix(11,10)+=-crLeftHandSideMatrix2*(DN(2,0)*crLeftHandSideMatrix8 + DN(2,1)*crLeftHandSideMatrix9);
    rLeftHandSideMatrix(11,11)+=0;
    rLeftHandSideMatrix(11,12)+=crLeftHandSideMatrix2*(DN(3,1)*crLeftHandSideMatrix7 + DN(3,2)*crLeftHandSideMatrix8);
    rLeftHandSideMatrix(11,13)+=-crLeftHandSideMatrix2*(DN(3,0)*crLeftHandSideMatrix7 - DN(3,2)*crLeftHandSideMatrix9);
    rLeftHandSideMatrix(11,14)+=-crLeftHandSideMatrix2*(DN(3,0)*crLeftHandSideMatrix8 + DN(3,1)*crLeftHandSideMatrix9);
    rLeftHandSideMatrix(11,15)+=0;
    rLeftHandSideMatrix(12,0)+=0;
    rLeftHandSideMatrix(12,1)+=0;
    rLeftHandSideMatrix(12,2)+=0;
    rLeftHandSideMatrix(12,3)+=0;
    rLeftHandSideMatrix(12,4)+=0;
    rLeftHandSideMatrix(12,5)+=0;
    rLeftHandSideMatrix(12,6)+=0;
    rLeftHandSideMatrix(12,7)+=0;
    rLeftHandSideMatrix(12,8)+=0;
    rLeftHandSideMatrix(12,9)+=0;
    rLeftHandSideMatrix(12,10)+=0;
    rLeftHandSideMatrix(12,11)+=0;
    rLeftHandSideMatrix(12,12)+=0;
    rLeftHandSideMatrix(12,13)+=0;
    rLeftHandSideMatrix(12,14)+=0;
    rLeftHandSideMatrix(12,15)+=0;
    rLeftHandSideMatrix(13,0)+=0;
    rLeftHandSideMatrix(13,1)+=0;
    rLeftHandSideMatrix(13,2)+=0;
    rLeftHandSideMatrix(13,3)+=0;
    rLeftHandSideMatrix(13,4)+=0;
    rLeftHandSideMatrix(13,5)+=0;
    rLeftHandSideMatrix(13,6)+=0;
    rLeftHandSideMatrix(13,7)+=0;
    rLeftHandSideMatrix(13,8)+=0;
    rLeftHandSideMatrix(13,9)+=0;
    rLeftHandSideMatrix(13,10)+=0;
    rLeftHandSideMatrix(13,11)+=0;
    rLeftHandSideMatrix(13,12)+=0;
    rLeftHandSideMatrix(13,13)+=0;
    rLeftHandSideMatrix(13,14)+=0;
    rLeftHandSideMatrix(13,15)+=0;
    rLeftHandSideMatrix(14,0)+=0;
    rLeftHandSideMatrix(14,1)+=0;
    rLeftHandSideMatrix(14,2)+=0;
    rLeftHandSideMatrix(14,3)+=0;
    rLeftHandSideMatrix(14,4)+=0;
    rLeftHandSideMatrix(14,5)+=0;
    rLeftHandSideMatrix(14,6)+=0;
    rLeftHandSideMatrix(14,7)+=0;
    rLeftHandSideMatrix(14,8)+=0;
    rLeftHandSideMatrix(14,9)+=0;
    rLeftHandSideMatrix(14,10)+=0;
    rLeftHandSideMatrix(14,11)+=0;
    rLeftHandSideMatrix(14,12)+=0;
    rLeftHandSideMatrix(14,13)+=0;
    rLeftHandSideMatrix(14,14)+=0;
    rLeftHandSideMatrix(14,15)+=0;
    rLeftHandSideMatrix(15,0)+=crLeftHandSideMatrix2*(DN(0,1)*crLeftHandSideMatrix10 + DN(0,2)*crLeftHandSideMatrix11);
    rLeftHandSideMatrix(15,1)+=-crLeftHandSideMatrix2*(DN(0,0)*crLeftHandSideMatrix10 - DN(0,2)*crLeftHandSideMatrix12);
    rLeftHandSideMatrix(15,2)+=-crLeftHandSideMatrix2*(DN(0,0)*crLeftHandSideMatrix11 + DN(0,1)*crLeftHandSideMatrix12);
    rLeftHandSideMatrix(15,3)+=0;
    rLeftHandSideMatrix(15,4)+=crLeftHandSideMatrix2*(DN(1,1)*crLeftHandSideMatrix10 + DN(1,2)*crLeftHandSideMatrix11);
    rLeftHandSideMatrix(15,5)+=-crLeftHandSideMatrix2*(DN(1,0)*crLeftHandSideMatrix10 - DN(1,2)*crLeftHandSideMatrix12);
    rLeftHandSideMatrix(15,6)+=-crLeftHandSideMatrix2*(DN(1,0)*crLeftHandSideMatrix11 + DN(1,1)*crLeftHandSideMatrix12);
    rLeftHandSideMatrix(15,7)+=0;
    rLeftHandSideMatrix(15,8)+=crLeftHandSideMatrix2*(DN(2,1)*crLeftHandSideMatrix10 + DN(2,2)*crLeftHandSideMatrix11);
    rLeftHandSideMatrix(15,9)+=-crLeftHandSideMatrix2*(DN(2,0)*crLeftHandSideMatrix10 - DN(2,2)*crLeftHandSideMatrix12);
    rLeftHandSideMatrix(15,10)+=-crLeftHandSideMatrix2*(DN(2,0)*crLeftHandSideMatrix11 + DN(2,1)*crLeftHandSideMatrix12);
    rLeftHandSideMatrix(15,11)+=0;
    rLeftHandSideMatrix(15,12)+=crLeftHandSideMatrix2*(DN(3,1)*crLeftHandSideMatrix10 + DN(3,2)*crLeftHandSideMatrix11);
    rLeftHandSideMatrix(15,13)+=-crLeftHandSideMatrix2*(DN(3,0)*crLeftHandSideMatrix10 - DN(3,2)*crLeftHandSideMatrix12);
    rLeftHandSideMatrix(15,14)+=-crLeftHandSideMatrix2*(DN(3,0)*crLeftHandSideMatrix11 + DN(3,1)*crLeftHandSideMatrix12);
    rLeftHandSideMatrix(15,15)+=0;
    

}

template <>
void FirstOrderStokesVariableViscosityBvsGl<2>::AddGaussPointRightHandSideContribution(
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
    const double crRightHandSideVector3 = DN(0,0)*u_nodes(0,0) + DN(1,0)*u_nodes(1,0) + DN(2,0)*u_nodes(2,0);
    const double crRightHandSideVector4 = N[0]*nu_nodes[0] + N[1]*nu_nodes[1] + N[2]*nu_nodes[2];
    const double crRightHandSideVector5 = DN(0,0)*crRightHandSideVector4;
    const double crRightHandSideVector6 = DN(0,1)*u_nodes(0,0) + DN(1,1)*u_nodes(1,0) + DN(2,1)*u_nodes(2,0);
    const double crRightHandSideVector7 = DN(0,1)*crRightHandSideVector4;
    const double crRightHandSideVector8 = DN(0,0)*nu_nodes[0] + DN(1,0)*nu_nodes[1] + DN(2,0)*nu_nodes[2];
    const double crRightHandSideVector9 = crRightHandSideVector3*crRightHandSideVector8;
    const double crRightHandSideVector10 = DN(0,0)*u_nodes(0,1) + DN(1,0)*u_nodes(1,1) + DN(2,0)*u_nodes(2,1);
    const double crRightHandSideVector11 = DN(0,1)*nu_nodes[0] + DN(1,1)*nu_nodes[1] + DN(2,1)*nu_nodes[2];
    const double crRightHandSideVector12 = crRightHandSideVector10*crRightHandSideVector11;
    const double crRightHandSideVector13 = crRightHandSideVector12 + crRightHandSideVector9;
    const double crRightHandSideVector14 = N[0]*f_nodes(0,1) + N[1]*f_nodes(1,1) + N[2]*f_nodes(2,1);
    const double crRightHandSideVector15 = sigma_gauss*(N[0]*u_nodes(0,1) + N[1]*u_nodes(1,1) + N[2]*u_nodes(2,1));
    const double crRightHandSideVector16 = DN(0,1)*u_nodes(0,1) + DN(1,1)*u_nodes(1,1) + DN(2,1)*u_nodes(2,1);
    const double crRightHandSideVector17 = crRightHandSideVector6*crRightHandSideVector8;
    const double crRightHandSideVector18 = crRightHandSideVector11*crRightHandSideVector16;
    const double crRightHandSideVector19 = crRightHandSideVector17 + crRightHandSideVector18;
    const double crRightHandSideVector20 = crRightHandSideVector16 + crRightHandSideVector3;
    const double crRightHandSideVector21 = DN(0,0)*p_nodes[0] + DN(1,0)*p_nodes[1] + DN(2,0)*p_nodes[2] - 2*crRightHandSideVector12 + crRightHandSideVector2 - 2*crRightHandSideVector9;
    const double crRightHandSideVector22 = DN(0,1)*p_nodes[0] + DN(1,1)*p_nodes[1] + DN(2,1)*p_nodes[2] + crRightHandSideVector15 - 2*crRightHandSideVector17 - 2*crRightHandSideVector18;
    const double crRightHandSideVector23 = DN(1,0)*crRightHandSideVector4;
    const double crRightHandSideVector24 = DN(1,1)*crRightHandSideVector4;
    const double crRightHandSideVector25 = DN(2,0)*crRightHandSideVector4;
    const double crRightHandSideVector26 = DN(2,1)*crRightHandSideVector4;
    rRightHandSideVector[0]+=w_g*(-DN(0,0)*crRightHandSideVector0 - N[0]*crRightHandSideVector1 - N[0]*crRightHandSideVector13 + N[0]*crRightHandSideVector2 + crRightHandSideVector3*crRightHandSideVector5 + crRightHandSideVector6*crRightHandSideVector7);
    rRightHandSideVector[1]+=w_g*(-DN(0,1)*crRightHandSideVector0 - N[0]*crRightHandSideVector14 + N[0]*crRightHandSideVector15 - N[0]*crRightHandSideVector19 + crRightHandSideVector10*crRightHandSideVector5 + crRightHandSideVector16*crRightHandSideVector7);
    rRightHandSideVector[2]+=w_g*(N[0]*crRightHandSideVector20 - delta_gauss*(DN(0,0)*crRightHandSideVector1 + DN(0,1)*crRightHandSideVector14) + delta_gauss*(DN(0,0)*crRightHandSideVector21 + DN(0,1)*crRightHandSideVector22));
    rRightHandSideVector[3]+=w_g*(-DN(1,0)*crRightHandSideVector0 - N[1]*crRightHandSideVector1 - N[1]*crRightHandSideVector13 + N[1]*crRightHandSideVector2 + crRightHandSideVector23*crRightHandSideVector3 + crRightHandSideVector24*crRightHandSideVector6);
    rRightHandSideVector[4]+=w_g*(-DN(1,1)*crRightHandSideVector0 - N[1]*crRightHandSideVector14 + N[1]*crRightHandSideVector15 - N[1]*crRightHandSideVector19 + crRightHandSideVector10*crRightHandSideVector23 + crRightHandSideVector16*crRightHandSideVector24);
    rRightHandSideVector[5]+=w_g*(N[1]*crRightHandSideVector20 - delta_gauss*(DN(1,0)*crRightHandSideVector1 + DN(1,1)*crRightHandSideVector14) + delta_gauss*(DN(1,0)*crRightHandSideVector21 + DN(1,1)*crRightHandSideVector22));
    rRightHandSideVector[6]+=w_g*(-DN(2,0)*crRightHandSideVector0 - N[2]*crRightHandSideVector1 - N[2]*crRightHandSideVector13 + N[2]*crRightHandSideVector2 + crRightHandSideVector25*crRightHandSideVector3 + crRightHandSideVector26*crRightHandSideVector6);
    rRightHandSideVector[7]+=w_g*(-DN(2,1)*crRightHandSideVector0 - N[2]*crRightHandSideVector14 + N[2]*crRightHandSideVector15 - N[2]*crRightHandSideVector19 + crRightHandSideVector10*crRightHandSideVector25 + crRightHandSideVector16*crRightHandSideVector26);
    rRightHandSideVector[8]+=w_g*(N[2]*crRightHandSideVector20 - delta_gauss*(DN(2,0)*crRightHandSideVector1 + DN(2,1)*crRightHandSideVector14) + delta_gauss*(DN(2,0)*crRightHandSideVector21 + DN(2,1)*crRightHandSideVector22));
    

}

template <>
void FirstOrderStokesVariableViscosityBvsGl<2>::AddConditionGaussPointRightHandSideContribution(
    const GlobalPointer<Condition>& rConditionPointer,
    const ElementDataContainer& rData,
    VectorType& rRightHandSideVector)
{
    const unsigned int TDim = 2;

    // Get element geometry
    const auto& r_geom = this->GetGeometry();
    const auto& r_condition_geom = rConditionPointer->GetGeometry();

    // Integration rule data: WE GET THE GAUSS POINTS FROM THE CONDITION, NOT THE ELEMENT
    // Note that we use the same for both velocity and pressure interpolations
    const auto integration_points = r_condition_geom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);
    const SizeType n_gauss = integration_points.size();

    std::vector<Vector> NContainer(n_gauss);
    std::vector<Matrix> DNContainer(n_gauss);
    Vector det_J_vect(n_gauss);
    for (IndexType g=0; g < n_gauss; g++) {
        const auto condition_local_coords = integration_points[g].Coordinates();
        Point::CoordinatesArrayType global_coords;
        r_condition_geom.GlobalCoordinates(global_coords, condition_local_coords);
        // Now we have the global coordinates of the gauss point in the condition
        // We can use these coordinates to find the corresponding local coordinates in the parent element
        Point::CoordinatesArrayType local_coords;
        r_geom.PointLocalCoordinates(local_coords, global_coords);
        // Now we have the local coordinates of the gauss point in the parent element
        // We can use these coordinates to evaluate shape functions and their gradients in the parent element
        
        // Get Shape functions evaluation and their derivatives at these local coords
        NContainer[g] = r_geom.ShapeFunctionsValues(NContainer[g], local_coords);
        Matrix DN_De;
        DN_De = r_geom.ShapeFunctionsLocalGradients(DN_De, local_coords);

        // Calculate Jacobians at integration points (with respect to the parent element)
        Matrix J;
        J = r_geom.Jacobian(J, local_coords);

        Matrix inv_J;
        double det_J;
        MathUtils<double>::InvertMatrix(J, inv_J, det_J);
        det_J_vect[g] = det_J;

        DNContainer[g].resize(NumNodes, TDim);
        DNContainer[g] = prod(DN_De, inv_J);
    }

    for (IndexType g=0; g < n_gauss; g++) {
        Vector N = NContainer[g];
        BoundedMatrix<double, NumNodes, TDim> DN = DNContainer[g];
        double w_g = det_J_vect[g] * integration_points[g].Weight();

        // Compute Unit normal
        array_1d<double, 3> n_gauss;
        n_gauss = r_condition_geom.Normal(integration_points[g].Coordinates());
        double A = norm_2(n_gauss);
        n_gauss /= A;

        array_1d<double, NumNodes> nu_nodes;
        BoundedMatrix<double, NumNodes, TDim> u_nodes;
        for (IndexType i = 0; i < NumNodes; ++i) {
            nu_nodes[i] = r_geom[i].FastGetSolutionStepValue(DYNAMIC_VISCOSITY);
            const auto& r_v = r_geom[i].FastGetSolutionStepValue(VELOCITY);
            for (IndexType d = 0; d < TDim; ++d) {
                u_nodes(i, d) = r_v[d];
            }
        }

        const double delta_gauss = rData.Delta;

                const double crRightHandSideVector0 = delta_gauss*w_g*(N[0]*nu_nodes[0] + N[1]*nu_nodes[1] + N[2]*nu_nodes[2])*(DN(0,0)*u_nodes(0,1) - DN(0,1)*u_nodes(0,0) + DN(1,0)*u_nodes(1,1) - DN(1,1)*u_nodes(1,0) + DN(2,0)*u_nodes(2,1) - DN(2,1)*u_nodes(2,0));
        rRightHandSideVector[0]+=0;
        rRightHandSideVector[1]+=0;
        rRightHandSideVector[2]+=crRightHandSideVector0*(DN(0,0)*n_gauss[1] - DN(0,1)*n_gauss[0]);
        rRightHandSideVector[3]+=0;
        rRightHandSideVector[4]+=0;
        rRightHandSideVector[5]+=crRightHandSideVector0*(DN(1,0)*n_gauss[1] - DN(1,1)*n_gauss[0]);
        rRightHandSideVector[6]+=0;
        rRightHandSideVector[7]+=0;
        rRightHandSideVector[8]+=crRightHandSideVector0*(DN(2,0)*n_gauss[1] - DN(2,1)*n_gauss[0]);
        

    }
}

template <>
void FirstOrderStokesVariableViscosityBvsGl<3>::AddGaussPointRightHandSideContribution(
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
    const double crRightHandSideVector3 = DN(0,0)*u_nodes(0,0) + DN(1,0)*u_nodes(1,0) + DN(2,0)*u_nodes(2,0) + DN(3,0)*u_nodes(3,0);
    const double crRightHandSideVector4 = N[0]*nu_nodes[0] + N[1]*nu_nodes[1] + N[2]*nu_nodes[2] + N[3]*nu_nodes[3];
    const double crRightHandSideVector5 = DN(0,0)*crRightHandSideVector4;
    const double crRightHandSideVector6 = DN(0,1)*u_nodes(0,0) + DN(1,1)*u_nodes(1,0) + DN(2,1)*u_nodes(2,0) + DN(3,1)*u_nodes(3,0);
    const double crRightHandSideVector7 = DN(0,1)*crRightHandSideVector4;
    const double crRightHandSideVector8 = DN(0,2)*u_nodes(0,0) + DN(1,2)*u_nodes(1,0) + DN(2,2)*u_nodes(2,0) + DN(3,2)*u_nodes(3,0);
    const double crRightHandSideVector9 = DN(0,2)*crRightHandSideVector4;
    const double crRightHandSideVector10 = DN(0,0)*nu_nodes[0] + DN(1,0)*nu_nodes[1] + DN(2,0)*nu_nodes[2] + DN(3,0)*nu_nodes[3];
    const double crRightHandSideVector11 = crRightHandSideVector10*crRightHandSideVector3;
    const double crRightHandSideVector12 = DN(0,0)*u_nodes(0,1) + DN(1,0)*u_nodes(1,1) + DN(2,0)*u_nodes(2,1) + DN(3,0)*u_nodes(3,1);
    const double crRightHandSideVector13 = DN(0,1)*nu_nodes[0] + DN(1,1)*nu_nodes[1] + DN(2,1)*nu_nodes[2] + DN(3,1)*nu_nodes[3];
    const double crRightHandSideVector14 = crRightHandSideVector12*crRightHandSideVector13;
    const double crRightHandSideVector15 = DN(0,0)*u_nodes(0,2) + DN(1,0)*u_nodes(1,2) + DN(2,0)*u_nodes(2,2) + DN(3,0)*u_nodes(3,2);
    const double crRightHandSideVector16 = DN(0,2)*nu_nodes[0] + DN(1,2)*nu_nodes[1] + DN(2,2)*nu_nodes[2] + DN(3,2)*nu_nodes[3];
    const double crRightHandSideVector17 = crRightHandSideVector15*crRightHandSideVector16;
    const double crRightHandSideVector18 = crRightHandSideVector11 + crRightHandSideVector14 + crRightHandSideVector17;
    const double crRightHandSideVector19 = N[0]*f_nodes(0,1) + N[1]*f_nodes(1,1) + N[2]*f_nodes(2,1) + N[3]*f_nodes(3,1);
    const double crRightHandSideVector20 = sigma_gauss*(N[0]*u_nodes(0,1) + N[1]*u_nodes(1,1) + N[2]*u_nodes(2,1) + N[3]*u_nodes(3,1));
    const double crRightHandSideVector21 = DN(0,1)*u_nodes(0,1) + DN(1,1)*u_nodes(1,1) + DN(2,1)*u_nodes(2,1) + DN(3,1)*u_nodes(3,1);
    const double crRightHandSideVector22 = DN(0,2)*u_nodes(0,1) + DN(1,2)*u_nodes(1,1) + DN(2,2)*u_nodes(2,1) + DN(3,2)*u_nodes(3,1);
    const double crRightHandSideVector23 = crRightHandSideVector10*crRightHandSideVector6;
    const double crRightHandSideVector24 = crRightHandSideVector13*crRightHandSideVector21;
    const double crRightHandSideVector25 = DN(0,1)*u_nodes(0,2) + DN(1,1)*u_nodes(1,2) + DN(2,1)*u_nodes(2,2) + DN(3,1)*u_nodes(3,2);
    const double crRightHandSideVector26 = crRightHandSideVector16*crRightHandSideVector25;
    const double crRightHandSideVector27 = crRightHandSideVector23 + crRightHandSideVector24 + crRightHandSideVector26;
    const double crRightHandSideVector28 = N[0]*f_nodes(0,2) + N[1]*f_nodes(1,2) + N[2]*f_nodes(2,2) + N[3]*f_nodes(3,2);
    const double crRightHandSideVector29 = sigma_gauss*(N[0]*u_nodes(0,2) + N[1]*u_nodes(1,2) + N[2]*u_nodes(2,2) + N[3]*u_nodes(3,2));
    const double crRightHandSideVector30 = DN(0,2)*u_nodes(0,2) + DN(1,2)*u_nodes(1,2) + DN(2,2)*u_nodes(2,2) + DN(3,2)*u_nodes(3,2);
    const double crRightHandSideVector31 = crRightHandSideVector10*crRightHandSideVector8;
    const double crRightHandSideVector32 = crRightHandSideVector13*crRightHandSideVector22;
    const double crRightHandSideVector33 = crRightHandSideVector16*crRightHandSideVector30;
    const double crRightHandSideVector34 = crRightHandSideVector31 + crRightHandSideVector32 + crRightHandSideVector33;
    const double crRightHandSideVector35 = crRightHandSideVector21 + crRightHandSideVector3 + crRightHandSideVector30;
    const double crRightHandSideVector36 = DN(0,0)*p_nodes[0] + DN(1,0)*p_nodes[1] + DN(2,0)*p_nodes[2] + DN(3,0)*p_nodes[3] - 2*crRightHandSideVector11 - 2*crRightHandSideVector14 - 2*crRightHandSideVector17 + crRightHandSideVector2;
    const double crRightHandSideVector37 = DN(0,1)*p_nodes[0] + DN(1,1)*p_nodes[1] + DN(2,1)*p_nodes[2] + DN(3,1)*p_nodes[3] + crRightHandSideVector20 - 2*crRightHandSideVector23 - 2*crRightHandSideVector24 - 2*crRightHandSideVector26;
    const double crRightHandSideVector38 = DN(0,2)*p_nodes[0] + DN(1,2)*p_nodes[1] + DN(2,2)*p_nodes[2] + DN(3,2)*p_nodes[3] + crRightHandSideVector29 - 2*crRightHandSideVector31 - 2*crRightHandSideVector32 - 2*crRightHandSideVector33;
    const double crRightHandSideVector39 = DN(1,0)*crRightHandSideVector4;
    const double crRightHandSideVector40 = DN(1,1)*crRightHandSideVector4;
    const double crRightHandSideVector41 = DN(1,2)*crRightHandSideVector4;
    const double crRightHandSideVector42 = DN(2,0)*crRightHandSideVector4;
    const double crRightHandSideVector43 = DN(2,1)*crRightHandSideVector4;
    const double crRightHandSideVector44 = DN(2,2)*crRightHandSideVector4;
    const double crRightHandSideVector45 = DN(3,0)*crRightHandSideVector4;
    const double crRightHandSideVector46 = DN(3,1)*crRightHandSideVector4;
    const double crRightHandSideVector47 = DN(3,2)*crRightHandSideVector4;
    rRightHandSideVector[0]+=w_g*(-DN(0,0)*crRightHandSideVector0 - N[0]*crRightHandSideVector1 - N[0]*crRightHandSideVector18 + N[0]*crRightHandSideVector2 + crRightHandSideVector3*crRightHandSideVector5 + crRightHandSideVector6*crRightHandSideVector7 + crRightHandSideVector8*crRightHandSideVector9);
    rRightHandSideVector[1]+=w_g*(-DN(0,1)*crRightHandSideVector0 - N[0]*crRightHandSideVector19 + N[0]*crRightHandSideVector20 - N[0]*crRightHandSideVector27 + crRightHandSideVector12*crRightHandSideVector5 + crRightHandSideVector21*crRightHandSideVector7 + crRightHandSideVector22*crRightHandSideVector9);
    rRightHandSideVector[2]+=w_g*(-DN(0,2)*crRightHandSideVector0 - N[0]*crRightHandSideVector28 + N[0]*crRightHandSideVector29 - N[0]*crRightHandSideVector34 + crRightHandSideVector15*crRightHandSideVector5 + crRightHandSideVector25*crRightHandSideVector7 + crRightHandSideVector30*crRightHandSideVector9);
    rRightHandSideVector[3]+=w_g*(N[0]*crRightHandSideVector35 - delta_gauss*(DN(0,0)*crRightHandSideVector1 + DN(0,1)*crRightHandSideVector19 + DN(0,2)*crRightHandSideVector28) + delta_gauss*(DN(0,0)*crRightHandSideVector36 + DN(0,1)*crRightHandSideVector37 + DN(0,2)*crRightHandSideVector38));
    rRightHandSideVector[4]+=w_g*(-DN(1,0)*crRightHandSideVector0 - N[1]*crRightHandSideVector1 - N[1]*crRightHandSideVector18 + N[1]*crRightHandSideVector2 + crRightHandSideVector3*crRightHandSideVector39 + crRightHandSideVector40*crRightHandSideVector6 + crRightHandSideVector41*crRightHandSideVector8);
    rRightHandSideVector[5]+=w_g*(-DN(1,1)*crRightHandSideVector0 - N[1]*crRightHandSideVector19 + N[1]*crRightHandSideVector20 - N[1]*crRightHandSideVector27 + crRightHandSideVector12*crRightHandSideVector39 + crRightHandSideVector21*crRightHandSideVector40 + crRightHandSideVector22*crRightHandSideVector41);
    rRightHandSideVector[6]+=w_g*(-DN(1,2)*crRightHandSideVector0 - N[1]*crRightHandSideVector28 + N[1]*crRightHandSideVector29 - N[1]*crRightHandSideVector34 + crRightHandSideVector15*crRightHandSideVector39 + crRightHandSideVector25*crRightHandSideVector40 + crRightHandSideVector30*crRightHandSideVector41);
    rRightHandSideVector[7]+=w_g*(N[1]*crRightHandSideVector35 - delta_gauss*(DN(1,0)*crRightHandSideVector1 + DN(1,1)*crRightHandSideVector19 + DN(1,2)*crRightHandSideVector28) + delta_gauss*(DN(1,0)*crRightHandSideVector36 + DN(1,1)*crRightHandSideVector37 + DN(1,2)*crRightHandSideVector38));
    rRightHandSideVector[8]+=w_g*(-DN(2,0)*crRightHandSideVector0 - N[2]*crRightHandSideVector1 - N[2]*crRightHandSideVector18 + N[2]*crRightHandSideVector2 + crRightHandSideVector3*crRightHandSideVector42 + crRightHandSideVector43*crRightHandSideVector6 + crRightHandSideVector44*crRightHandSideVector8);
    rRightHandSideVector[9]+=w_g*(-DN(2,1)*crRightHandSideVector0 - N[2]*crRightHandSideVector19 + N[2]*crRightHandSideVector20 - N[2]*crRightHandSideVector27 + crRightHandSideVector12*crRightHandSideVector42 + crRightHandSideVector21*crRightHandSideVector43 + crRightHandSideVector22*crRightHandSideVector44);
    rRightHandSideVector[10]+=w_g*(-DN(2,2)*crRightHandSideVector0 - N[2]*crRightHandSideVector28 + N[2]*crRightHandSideVector29 - N[2]*crRightHandSideVector34 + crRightHandSideVector15*crRightHandSideVector42 + crRightHandSideVector25*crRightHandSideVector43 + crRightHandSideVector30*crRightHandSideVector44);
    rRightHandSideVector[11]+=w_g*(N[2]*crRightHandSideVector35 - delta_gauss*(DN(2,0)*crRightHandSideVector1 + DN(2,1)*crRightHandSideVector19 + DN(2,2)*crRightHandSideVector28) + delta_gauss*(DN(2,0)*crRightHandSideVector36 + DN(2,1)*crRightHandSideVector37 + DN(2,2)*crRightHandSideVector38));
    rRightHandSideVector[12]+=w_g*(-DN(3,0)*crRightHandSideVector0 - N[3]*crRightHandSideVector1 - N[3]*crRightHandSideVector18 + N[3]*crRightHandSideVector2 + crRightHandSideVector3*crRightHandSideVector45 + crRightHandSideVector46*crRightHandSideVector6 + crRightHandSideVector47*crRightHandSideVector8);
    rRightHandSideVector[13]+=w_g*(-DN(3,1)*crRightHandSideVector0 - N[3]*crRightHandSideVector19 + N[3]*crRightHandSideVector20 - N[3]*crRightHandSideVector27 + crRightHandSideVector12*crRightHandSideVector45 + crRightHandSideVector21*crRightHandSideVector46 + crRightHandSideVector22*crRightHandSideVector47);
    rRightHandSideVector[14]+=w_g*(-DN(3,2)*crRightHandSideVector0 - N[3]*crRightHandSideVector28 + N[3]*crRightHandSideVector29 - N[3]*crRightHandSideVector34 + crRightHandSideVector15*crRightHandSideVector45 + crRightHandSideVector25*crRightHandSideVector46 + crRightHandSideVector30*crRightHandSideVector47);
    rRightHandSideVector[15]+=w_g*(N[3]*crRightHandSideVector35 - delta_gauss*(DN(3,0)*crRightHandSideVector1 + DN(3,1)*crRightHandSideVector19 + DN(3,2)*crRightHandSideVector28) + delta_gauss*(DN(3,0)*crRightHandSideVector36 + DN(3,1)*crRightHandSideVector37 + DN(3,2)*crRightHandSideVector38));
    
    
}

template <>
void FirstOrderStokesVariableViscosityBvsGl<3>::AddConditionGaussPointRightHandSideContribution(
    const GlobalPointer<Condition>& rConditionPointer,
    const ElementDataContainer& rData,
    VectorType& rRHS)
{
    KRATOS_ERROR << "Method AddConditionGaussPointRightHandSideContribution for element type FirstOrderStokesVariableViscosityBvsGl is not implemented for 3D yet" << std::endl;

    const auto& u_nodes = rData.Velocity;
    const auto& nu_nodes = rData.DynamicViscosity;

    const auto& N = rData.N;
    const auto& DN = rData.DN;

    // Assemble LHS contribution
    const double w_g = rData.Weight;

    // Get stabilization data
    const double delta_gauss = rData.Delta;

    const Vector n_gauss = ZeroVector(3);

    Vector rRightHandSideVector = ZeroVector(LocalSize);

        const double crRightHandSideVector0 = DN(0,0)*u_nodes(0,1) - DN(0,1)*u_nodes(0,0) + DN(1,0)*u_nodes(1,1) - DN(1,1)*u_nodes(1,0) + DN(2,0)*u_nodes(2,1) - DN(2,1)*u_nodes(2,0) + DN(3,0)*u_nodes(3,1) - DN(3,1)*u_nodes(3,0);
    const double crRightHandSideVector1 = DN(0,0)*u_nodes(0,2) - DN(0,2)*u_nodes(0,0) + DN(1,0)*u_nodes(1,2) - DN(1,2)*u_nodes(1,0) + DN(2,0)*u_nodes(2,2) - DN(2,2)*u_nodes(2,0) + DN(3,0)*u_nodes(3,2) - DN(3,2)*u_nodes(3,0);
    const double crRightHandSideVector2 = DN(0,1)*u_nodes(0,2) - DN(0,2)*u_nodes(0,1) + DN(1,1)*u_nodes(1,2) - DN(1,2)*u_nodes(1,1) + DN(2,1)*u_nodes(2,2) - DN(2,2)*u_nodes(2,1) + DN(3,1)*u_nodes(3,2) - DN(3,2)*u_nodes(3,1);
    const double crRightHandSideVector3 = delta_gauss*w_g*(N[0]*nu_nodes[0] + N[1]*nu_nodes[1] + N[2]*nu_nodes[2] + N[3]*nu_nodes[3]);
    rRightHandSideVector[0]+=0;
    rRightHandSideVector[1]+=0;
    rRightHandSideVector[2]+=0;
    rRightHandSideVector[3]+=crRightHandSideVector3*(crRightHandSideVector0*(DN(0,0)*n_gauss[1] - DN(0,1)*n_gauss[0]) + crRightHandSideVector1*(DN(0,0)*n_gauss[2] - DN(0,2)*n_gauss[0]) + crRightHandSideVector2*(DN(0,1)*n_gauss[2] - DN(0,2)*n_gauss[1]));
    rRightHandSideVector[4]+=0;
    rRightHandSideVector[5]+=0;
    rRightHandSideVector[6]+=0;
    rRightHandSideVector[7]+=crRightHandSideVector3*(crRightHandSideVector0*(DN(1,0)*n_gauss[1] - DN(1,1)*n_gauss[0]) + crRightHandSideVector1*(DN(1,0)*n_gauss[2] - DN(1,2)*n_gauss[0]) + crRightHandSideVector2*(DN(1,1)*n_gauss[2] - DN(1,2)*n_gauss[1]));
    rRightHandSideVector[8]+=0;
    rRightHandSideVector[9]+=0;
    rRightHandSideVector[10]+=0;
    rRightHandSideVector[11]+=crRightHandSideVector3*(crRightHandSideVector0*(DN(2,0)*n_gauss[1] - DN(2,1)*n_gauss[0]) + crRightHandSideVector1*(DN(2,0)*n_gauss[2] - DN(2,2)*n_gauss[0]) + crRightHandSideVector2*(DN(2,1)*n_gauss[2] - DN(2,2)*n_gauss[1]));
    rRightHandSideVector[12]+=0;
    rRightHandSideVector[13]+=0;
    rRightHandSideVector[14]+=0;
    rRightHandSideVector[15]+=crRightHandSideVector3*(crRightHandSideVector0*(DN(3,0)*n_gauss[1] - DN(3,1)*n_gauss[0]) + crRightHandSideVector1*(DN(3,0)*n_gauss[2] - DN(3,2)*n_gauss[0]) + crRightHandSideVector2*(DN(3,1)*n_gauss[2] - DN(3,2)*n_gauss[1]));
    

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private serialization

template< unsigned int TDim >
void FirstOrderStokesVariableViscosityBvsGl<TDim>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}


template< unsigned int TDim >
void FirstOrderStokesVariableViscosityBvsGl<TDim>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class FirstOrderStokesVariableViscosityBvsGl<2>;
template class FirstOrderStokesVariableViscosityBvsGl<3>;

}