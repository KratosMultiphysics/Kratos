//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Nicolas Sibuet Ruiz
//

// System includes


// External includes


// Project includes
#include "utilities/element_size_calculator.h"

// Application includes
#include "thermally_driven_active_fluid.h"

namespace Kratos
{

///////////////////////////////////////////////////////////////////////////////////////////////////
// Life cycle

template< unsigned int TDim >
ThermallyDrivenActiveFluid<TDim>::ThermallyDrivenActiveFluid(IndexType NewId)
    : Element(NewId)
{}

template< unsigned int TDim >
ThermallyDrivenActiveFluid<TDim>::ThermallyDrivenActiveFluid(
    IndexType NewId,
    const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes)
{}

template< unsigned int TDim >
ThermallyDrivenActiveFluid<TDim>::ThermallyDrivenActiveFluid(
    IndexType NewId,
    GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{}

template< unsigned int TDim >
ThermallyDrivenActiveFluid<TDim>::ThermallyDrivenActiveFluid(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    Properties::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{}

template< unsigned int TDim >
ThermallyDrivenActiveFluid<TDim>::~ThermallyDrivenActiveFluid()
{}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template< unsigned int TDim >
Element::Pointer ThermallyDrivenActiveFluid<TDim>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<ThermallyDrivenActiveFluid<TDim>>(NewId, this->GetGeometry().Create(ThisNodes), pProperties);
}

template< unsigned int TDim >
Element::Pointer ThermallyDrivenActiveFluid<TDim>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    Properties::Pointer pProperties) const
{
    return Kratos::make_intrusive<ThermallyDrivenActiveFluid<TDim>>(NewId, pGeom, pProperties);
}

template< unsigned int TDim >
void ThermallyDrivenActiveFluid<TDim>::EquationIdVector(
    EquationIdVectorType &rResult,
    const ProcessInfo &rCurrentProcessInfo) const
{
    if (rResult.size() != LocalSize) {
        rResult.resize(LocalSize, false);
    }

    IndexType local_index = 0;
    const auto& r_geometry = this->GetGeometry();
    for (IndexType i = 0; i < VelocityNumNodes; ++i) {
        rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_X).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_Y).EquationId();
        if constexpr (TDim == 3) {
            rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_Z).EquationId();
        }
        rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_LAPLACIAN_X).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_LAPLACIAN_Y).EquationId();
        if constexpr (TDim == 3) {
            rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_LAPLACIAN_Z).EquationId();
        }
        rResult[local_index++] = r_geometry[i].GetDof(TEMPERATURE).EquationId();
    }
    
    for (IndexType i = 0; i < PressureNumNodes; ++i) {
        rResult[local_index++] = r_geometry[i].GetDof(PRESSURE).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(PRESSUREAUX).EquationId();
    }
}

template< unsigned int TDim >
void ThermallyDrivenActiveFluid<TDim>::GetDofList(
    DofsVectorType &rElementalDofList,
    const ProcessInfo &rCurrentProcessInfo) const
{
    if (rElementalDofList.size() != LocalSize) {
        rElementalDofList.resize(LocalSize);
    }

    IndexType local_index = 0;
    const auto& r_geometry = this->GetGeometry();
    for (IndexType i = 0; i < VelocityNumNodes; ++i) {
        rElementalDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_X);
        rElementalDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_Y);
        if constexpr (TDim == 3) {
            rElementalDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_Z);
        }
        rElementalDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_LAPLACIAN_X);
        rElementalDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_LAPLACIAN_Y);
        if constexpr (TDim == 3) {
            rElementalDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_LAPLACIAN_Z);
        }
        rElementalDofList[local_index++] = r_geometry[i].pGetDof(TEMPERATURE);
    }
    
    for (IndexType i = 0; i < PressureNumNodes; ++i) {
        rElementalDofList[local_index++] = r_geometry[i].pGetDof(PRESSURE);
        rElementalDofList[local_index++] = r_geometry[i].pGetDof(PRESSUREAUX);
    }
}

template< unsigned int TDim >
void ThermallyDrivenActiveFluid<TDim>::CalculateLocalSystem(
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

    // Calculate kinematics
    Vector weights;
    Matrix velocity_N;
    Matrix pressure_N;
    GeometryType::ShapeFunctionsGradientsType velocity_DN;
    GeometryType::ShapeFunctionsGradientsType pressure_DN;
    CalculateKinematics(weights, velocity_N, pressure_N, velocity_DN, pressure_DN);

    // Loop Gauss points
    const auto& r_geom = this->GetGeometry();
    const SizeType n_gauss = r_geom.IntegrationPointsNumber(IntegrationMethod);
    for (IndexType g = 0; g < n_gauss; ++g) {
        // Set current Gauss point kinematics
        noalias(aux_data.N_v) = row(velocity_N, g);
        noalias(aux_data.N_p) = row(pressure_N, g);
        noalias(aux_data.DN_v) = velocity_DN[g];
        noalias(aux_data.DN_p) = pressure_DN[g];
        aux_data.Weight = weights[g];

        // Assemble standard Galerkin contribution
        AddGaussPointLeftHandSideContribution(aux_data, rLeftHandSideMatrix);
        AddGaussPointRightHandSideContribution(aux_data, rRightHandSideVector);

        AddGaussPointRightHandSideMomentumForcingContribution(aux_data, rRightHandSideVector);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Inquiry

template< unsigned int TDim >
int ThermallyDrivenActiveFluid<TDim>::Check(const ProcessInfo &rCurrentProcessInfo) const
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
const Parameters ThermallyDrivenActiveFluid<TDim>::GetSpecifications() const
{
    const Parameters specifications = Parameters(R"({
        "time_integration"           : ["implicit"],
        "framework"                  : "ale",
        "symmetric_lhs"              : false,
        "positive_definite_lhs"      : false,
        "output"                     : {
            "gauss_point"            : [""],
            "nodal_historical"       : ["VELOCITY","VELOCITY_LAPLACIAN","TEMPERATURE","PRESSURE","BODY_FORCE","HEAT_FLUX"],
            "nodal_non_historical"   : ["PRESSUREAUX"],
            "entity"                 : []
        },
        "required_variables"         : ["VELOCITY","VELOCITY_LAPLACIAN","TEMPERATURE","PRESSURE","PRESSUREAUX","BODY_FORCE","HEAT_FLUX","REACTION","REACTION_WATER_PRESSURE"]
        "required_dofs"              : [],
        "flags_used"                 : [],
        "compatible_geometries"      : ["Triangle2D6","Tetrahedra3D10"],
        "element_integrates_in_time" : true,
        "required_polynomial_degree_of_geometry" : 2,
        "documentation"   :
            "This implements a divergence constrained mixed element for thermally driven active fluid models. Taken from paper with DOI: https://doi.org/10.1016/j.apnum.2026.04.008"
    })");

    if (TDim == 2) {
        std::vector<std::string> dofs_2d({"VELOCITY_X","VELOCITY_Y","VELOCITY_LAPLACIAN_X","VELOCITY_LAPLACIAN_Y","TEMPERATURE","PRESSURE","PRESSUREAUX"});
        specifications["required_dofs"].SetStringArray(dofs_2d);
    } else {
        std::vector<std::string> dofs_3d({"VELOCITY_X","VELOCITY_Y","VELOCITY_Z","VELOCITY_LAPLACIAN_X","VELOCITY_LAPLACIAN_Y","VELOCITY_LAPLACIAN_Z","TEMPERATURE","PRESSURE","PRESSUREAUX"});
        specifications["required_dofs"].SetStringArray(dofs_3d);
    }

    return specifications;
}

template< unsigned int TDim >
std::string ThermallyDrivenActiveFluid<TDim>::Info() const
{
    std::stringstream buffer;
    buffer << "ThermallyDrivenActiveFluid" << TDim << "D" << VelocityNumNodes << "N #" << this->Id();
    return buffer.str();
}

template< unsigned int TDim >
void ThermallyDrivenActiveFluid<TDim>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info() << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private operations

template <unsigned int TDim>
void ThermallyDrivenActiveFluid<TDim>::SetElementData(
    const ProcessInfo& rProcessInfo,
    ElementDataContainer &rData)
{
    // Set nodal data. First, for P2 variables taking the values in all nodes
    const auto& r_geom = this->GetGeometry();
    for (IndexType i = 0; i < VelocityNumNodes; ++i) {
        const auto& r_v = r_geom[i].FastGetSolutionStepValue(VELOCITY);
        const auto& r_v_n = r_geom[i].FastGetSolutionStepValue(VELOCITY, 1);
        const auto& r_v_nn = r_geom[i].FastGetSolutionStepValue(VELOCITY, 2);
        const auto& r_vlap = r_geom[i].FastGetSolutionStepValue(VELOCITY_LAPLACIAN);
        const auto& r_vlap_n = r_geom[i].FastGetSolutionStepValue(VELOCITY_LAPLACIAN, 1);
        const auto& r_vlap_nn = r_geom[i].FastGetSolutionStepValue(VELOCITY_LAPLACIAN, 2);
        const auto& r_F_u = r_geom[i].FastGetSolutionStepValue(BODY_FORCE);
        const auto& r_F_u_n = r_geom[i].FastGetSolutionStepValue(BODY_FORCE, 1);
        const auto& r_F_u_nn = r_geom[i].FastGetSolutionStepValue(BODY_FORCE, 2);

        for (IndexType d = 0; d < TDim; ++d) {

            // Nodal values
            rData.Velocity(i, d) = r_v[d];
            rData.VelocityOld1(i, d) = r_v_n[d];
            rData.VelocityOld2(i, d) = r_v_nn[d];
            rData.LaplaceVel(i, d) = r_vlap[d];
            rData.LaplaceVelOld1(i, d) = r_vlap_n[d];
            rData.LaplaceVelOld2(i, d) = r_vlap_nn[d];
            rData.MomentumForcing(i, d) = r_F_u[d];
            rData.MomentumForcingOld1(i, d) = r_F_u_n[d];
            rData.MomentumForcingOld2(i, d) = r_F_u_nn[d];

        }

        rData.Temperature[i]  = r_geom[i].FastGetSolutionStepValue(TEMPERATURE);
        rData.TemperatureOld1[i]  = r_geom[i].FastGetSolutionStepValue(TEMPERATURE, 1);
        rData.TemperatureOld2[i]  = r_geom[i].FastGetSolutionStepValue(TEMPERATURE, 2);
        rData.BodySource[i] = r_geom[i].FastGetSolutionStepValue(HEAT_FLUX);
        rData.BodySourceOld1[i] = r_geom[i].FastGetSolutionStepValue(HEAT_FLUX, 1);
        rData.BodySourceOld2[i] = r_geom[i].FastGetSolutionStepValue(HEAT_FLUX, 2);
        
    }

    // Set nodal data. Second, for P1 variables taking only the values from the first 3 nodes (corners)
    for (IndexType i = 0; i < PressureNumNodes; ++i) {
        rData.Pressure[i] = r_geom[i].FastGetSolutionStepValue(PRESSURE);
        rData.PressureOld1[i] = r_geom[i].FastGetSolutionStepValue(PRESSURE, 1);
        rData.PressureOld2[i] = r_geom[i].FastGetSolutionStepValue(PRESSURE, 2);
        rData.Phi[i] =  r_geom[i].FastGetSolutionStepValue(PRESSUREAUX);
    }

    // Set time integration (DML) values

    const double delta_time = rProcessInfo[DELTA_TIME];
    double delta_time_old;

    delta_time_old = rProcessInfo.GetPreviousSolutionStepInfo()[DELTA_TIME];

    const double dln_theta = rProcessInfo[TIME_INTEGRATION_THETA];
    const double dln_epsilon = (delta_time-delta_time_old)/(delta_time+delta_time_old);

    const double alpha_0 = 1/2*(dln_theta-1);
    const double alpha_2 = 1/2*(dln_theta+1);
    rData.alpha_0 = alpha_0;
    rData.alpha_1 = -dln_theta;
    rData.alpha_2 = alpha_2;
    const double beta_aux = (1-std::pow(dln_theta,2.0))/std::pow(1+dln_epsilon*dln_theta,2);
    rData.beta_0 = 1/4*(1+beta_aux-std::pow(dln_epsilon,2.0)*dln_theta*beta_aux-dln_theta);
    rData.beta_1 = 1/2*(1-beta_aux);
    rData.beta_2 = 1/4*(1+beta_aux+std::pow(dln_epsilon,2.0)*dln_theta*beta_aux+dln_theta);

    rData.k_hat_n = alpha_2*delta_time-alpha_0*delta_time_old;

    // Set material parameters
    rData.Viscosity = this->GetProperties().GetValue(VISCOSITY);
    rData.Density = this->GetProperties().GetValue(DENSITY);
    rData.ThermalConductivity = this->GetProperties().GetValue(CONDUCTIVITY);

    // Set formulation parameters
    rData.GammaStabilityCoefficient = rProcessInfo[STABILIZATION_FACTOR];
    rData.RhoParameter = rProcessInfo[Y1];
    rData.LambdaParameter = rProcessInfo[Y2];
    rData.SigmaBuoyancyParameter = rProcessInfo[YF];
    const auto& gravity_unary_direction = rProcessInfo[GRAVITY];
    for (IndexType d = 0; d < TDim; ++d) {
        rData.GravityUnitaryDirection[d] = gravity_unary_direction[d];
    }
}

template< unsigned int TDim >
void ThermallyDrivenActiveFluid<TDim>::CalculateKinematics(
    Vector& rGaussWeights,
    Matrix& rVelocityN,
    Matrix& rPressureN,
    GeometryType::ShapeFunctionsGradientsType& rVelocityDNDX,
    GeometryType::ShapeFunctionsGradientsType& rPressureDNDX)
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

    // Calculate velocity kinematics from the geometry (P2 interpolation)
    rVelocityN = r_geom.ShapeFunctionsValues(IntegrationMethod);
    const auto& r_DN_De_v = r_geom.ShapeFunctionsLocalGradients(IntegrationMethod);
    if (rVelocityDNDX.size() != n_gauss) {
        rVelocityDNDX.resize(n_gauss, false);
    }
    for (IndexType g = 0; g < n_gauss; ++g) {
        rVelocityDNDX[g] = prod(r_DN_De_v[g], inv_J_vect[g]);
    }

    // Calculate pressure kinematics from an auxiliary geometry (P1 interpolation)
    GeometryType::UniquePointer p_aux_geom = nullptr;
    if constexpr (TDim == 2) {
        p_aux_geom = Kratos::make_unique<Triangle2D3<NodeType>>(r_geom(0), r_geom(1), r_geom(2));
    } else {
        p_aux_geom = Kratos::make_unique<Tetrahedra3D4<NodeType>>(r_geom(0), r_geom(1), r_geom(2), r_geom(3));
    }
    rPressureN = p_aux_geom->ShapeFunctionsValues(IntegrationMethod);
    if (rPressureDNDX.size() != n_gauss) {
        rPressureDNDX.resize(n_gauss, false);
    }
    const auto& r_DN_De_p = p_aux_geom->ShapeFunctionsLocalGradients(IntegrationMethod);
    for (IndexType g = 0; g < n_gauss; ++g) {
        rPressureDNDX[g] = prod(r_DN_De_p[g], inv_J_vect[g]);
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
void ThermallyDrivenActiveFluid<2>::AddGaussPointLeftHandSideContribution(
    const ElementDataContainer& rData,
    MatrixType& rLeftHandSideMatrix)
{
    const SizeType TDim = 2;

    // Get material data
    const double mu_gauss = rData.Viscosity;
    const double nu_gauss = rData.Density;
    const double kappa_gauss = rData.ThermalConductivity;

    // Get formulation parameters
    const double gamma_gauss = rData.GammaStabilityCoefficient;
    const double rgo_gauss = rData.RhoParameter;
    const double lambda_gauss = rData.LambdaParameter;
    const double sigma_gauss = rData.SigmaBuoyancyParameter;
    const array_1d<double, TDim> xi_gauss = rData.GravityUnitaryDirection;

    // Get time integration scheme data
    const double alpha0_gauss = rData.alpha_0;
    const double alpha1_gauss = rData.alpha_1;
    const double alpha2_gauss = rData.alpha_2;
    const double beta0_gauss = rData.beta_0;
    const double beta1_gauss = rData.beta_1;
    const double beta2_gauss = rData.beta_2;

    const double k_hat_n_gauss = rData.k_hat_n;

    // Get shape function values
    const auto& N_p1 = rData.N_p;
    const auto& N_p2 = rData.N_v;
    const auto& DN_p1 = rData.DN_p;
    const auto& DN_p2 = rData.DN_v;

    // Get nodal values
    const auto& u_nodes = rData.Velocity;
    const auto& u_n_nodes = rData.VelocityOld1;
    const auto& u_nn_nodes = rData.VelocityOld2;
    const auto& w_nodes = rData.LaplaceVel;
    const auto& w_n_nodes = rData.LaplaceVelOld1;
    const auto& w_nn_nodes = rData.LaplaceVelOld2;
    const auto& T_nodes = rData.Temperature;
    const auto& T_n_nodes = rData.TemperatureOld1;
    const auto& T_nn_nodes = rData.TemperatureOld2;
    const auto& p_nodes = rData.Pressure;
    const auto& p_n_nodes = rData.PressureOld1;
    const auto& p_nn_nodes = rData.PressureOld2;
    const auto& Phi_nodes = rData.Phi;
    const auto& f_nodes = rData.BodySource;
    const auto& f_n_nodes = rData.BodySourceOld1;
    const auto& f_nn_nodes = rData.BodySourceOld2;

    // Assemble LHS contribution
    const double w_g = rData.Weight;

    //substitute_element_lhs_2D_p2-6N_p1-3N

}

template <>
void ThermallyDrivenActiveFluid<3>::AddGaussPointLeftHandSideContribution(
    const ElementDataContainer& rData,
    MatrixType& rLeftHandSideMatrix)
{
    const SizeType TDim = 3;

    // Get material data
    const double mu_gauss = rData.Viscosity;
    const double nu_gauss = rData.Density;
    const double kappa_gauss = rData.ThermalConductivity;

    // Get formulation parameters
    const double gamma_gauss = rData.GammaStabilityCoefficient;
    const double rgo_gauss = rData.RhoParameter;
    const double lambda_gauss = rData.LambdaParameter;
    const double sigma_gauss = rData.SigmaBuoyancyParameter;
    const array_1d<double, TDim> xi_gauss = rData.GravityUnitaryDirection;

    // Get time integration scheme data
    const double alpha0_gauss = rData.alpha_0;
    const double alpha1_gauss = rData.alpha_1;
    const double alpha2_gauss = rData.alpha_2;
    const double beta0_gauss = rData.beta_0;
    const double beta1_gauss = rData.beta_1;
    const double beta2_gauss = rData.beta_2;

    const double k_hat_n_gauss = rData.k_hat_n;

    // Get shape function values
    const auto& N_p1 = rData.N_p;
    const auto& N_p2 = rData.N_v;
    const auto& DN_p1 = rData.DN_p;
    const auto& DN_p2 = rData.DN_v;

    // Get nodal values
    const auto& u_nodes = rData.Velocity;
    const auto& u_n_nodes = rData.VelocityOld1;
    const auto& u_nn_nodes = rData.VelocityOld2;
    const auto& w_nodes = rData.LaplaceVel;
    const auto& w_n_nodes = rData.LaplaceVelOld1;
    const auto& w_nn_nodes = rData.LaplaceVelOld2;
    const auto& T_nodes = rData.Temperature;
    const auto& T_n_nodes = rData.TemperatureOld1;
    const auto& T_nn_nodes = rData.TemperatureOld2;
    const auto& p_nodes = rData.Pressure;
    const auto& p_n_nodes = rData.PressureOld1;
    const auto& p_nn_nodes = rData.PressureOld2;
    const auto& Phi_nodes = rData.Phi;
    const auto& f_nodes = rData.BodySource;
    const auto& f_n_nodes = rData.BodySourceOld1;
    const auto& f_nn_nodes = rData.BodySourceOld2;

    // Assemble LHS contribution
    const double w_g = rData.Weight;

    //substitute_element_lhs_3D_p2-10N_p1-4N

}

template <>
void ThermallyDrivenActiveFluid<2>::AddGaussPointRightHandSideContribution(
    const ElementDataContainer& rData,
    VectorType& rRightHandSideVector)
{
    const SizeType TDim = 2;

    // Get material data
    const double mu_gauss = rData.Viscosity;
    const double nu_gauss = rData.Density;
    const double kappa_gauss = rData.ThermalConductivity;

    // Get formulation parameters
    const double gamma_gauss = rData.GammaStabilityCoefficient;
    const double rgo_gauss = rData.RhoParameter;
    const double lambda_gauss = rData.LambdaParameter;
    const double sigma_gauss = rData.SigmaBuoyancyParameter;
    const array_1d<double, TDim> xi_gauss = rData.GravityUnitaryDirection;

    // Get time integration scheme data
    const double alpha0_gauss = rData.alpha_0;
    const double alpha1_gauss = rData.alpha_1;
    const double alpha2_gauss = rData.alpha_2;
    const double beta0_gauss = rData.beta_0;
    const double beta1_gauss = rData.beta_1;
    const double beta2_gauss = rData.beta_2;

    const double k_hat_n_gauss = rData.k_hat_n;

    // Get shape function values
    const auto& N_p1 = rData.N_p;
    const auto& N_p2 = rData.N_v;
    const auto& DN_p1 = rData.DN_p;
    const auto& DN_p2 = rData.DN_v;

    // Get nodal values
    const auto& u_nodes = rData.Velocity;
    const auto& u_n_nodes = rData.VelocityOld1;
    const auto& u_nn_nodes = rData.VelocityOld2;
    const auto& w_nodes = rData.LaplaceVel;
    const auto& w_n_nodes = rData.LaplaceVelOld1;
    const auto& w_nn_nodes = rData.LaplaceVelOld2;
    const auto& T_nodes = rData.Temperature;
    const auto& T_n_nodes = rData.TemperatureOld1;
    const auto& T_nn_nodes = rData.TemperatureOld2;
    const auto& p_nodes = rData.Pressure;
    const auto& p_n_nodes = rData.PressureOld1;
    const auto& p_nn_nodes = rData.PressureOld2;
    const auto& Phi_nodes = rData.Phi;
    const auto& f_nodes = rData.BodySource;
    const auto& f_n_nodes = rData.BodySourceOld1;
    const auto& f_nn_nodes = rData.BodySourceOld2;

    // Assemble RHS contribution
    const double w_g = rData.Weight;

    //substitute_element_rhs_2D_p2-6N_p1-3N

}

template <>
void ThermallyDrivenActiveFluid<3>::AddGaussPointRightHandSideContribution(
    const ElementDataContainer& rData,
    VectorType& rRightHandSideVector)
{
    const SizeType TDim = 3;

    // Get material data
    const double mu_gauss = rData.Viscosity;
    const double nu_gauss = rData.Density;
    const double kappa_gauss = rData.ThermalConductivity;

    // Get formulation parameters
    const double gamma_gauss = rData.GammaStabilityCoefficient;
    const double rgo_gauss = rData.RhoParameter;
    const double lambda_gauss = rData.LambdaParameter;
    const double sigma_gauss = rData.SigmaBuoyancyParameter;
    const array_1d<double, TDim> xi_gauss = rData.GravityUnitaryDirection;

    // Get time integration scheme data
    const double alpha0_gauss = rData.alpha_0;
    const double alpha1_gauss = rData.alpha_1;
    const double alpha2_gauss = rData.alpha_2;
    const double beta0_gauss = rData.beta_0;
    const double beta1_gauss = rData.beta_1;
    const double beta2_gauss = rData.beta_2;

    const double k_hat_n_gauss = rData.k_hat_n;

    // Get shape function values
    const auto& N_p1 = rData.N_p;
    const auto& N_p2 = rData.N_v;
    const auto& DN_p1 = rData.DN_p;
    const auto& DN_p2 = rData.DN_v;

    // Get nodal values
    const auto& u_nodes = rData.Velocity;
    const auto& u_n_nodes = rData.VelocityOld1;
    const auto& u_nn_nodes = rData.VelocityOld2;
    const auto& w_nodes = rData.LaplaceVel;
    const auto& w_n_nodes = rData.LaplaceVelOld1;
    const auto& w_nn_nodes = rData.LaplaceVelOld2;
    const auto& T_nodes = rData.Temperature;
    const auto& T_n_nodes = rData.TemperatureOld1;
    const auto& T_nn_nodes = rData.TemperatureOld2;
    const auto& p_nodes = rData.Pressure;
    const auto& p_n_nodes = rData.PressureOld1;
    const auto& p_nn_nodes = rData.PressureOld2;
    const auto& Phi_nodes = rData.Phi;
    const auto& f_nodes = rData.BodySource;
    const auto& f_n_nodes = rData.BodySourceOld1;
    const auto& f_nn_nodes = rData.BodySourceOld2;

    // Assemble RHS contribution
    const double w_g = rData.Weight;

    //substitute_element_rhs_3D_p2-10N_p1-4N

}

template <>
void ThermallyDrivenActiveFluid<2>::AddGaussPointRightHandSideMomentumForcingContribution(
    const ElementDataContainer& rData,
    VectorType& rRightHandSideVector)
{
    const SizeType TDim = 2;

    // Get shape function values
    const auto& N_p2 = rData.N_v;
    const auto& DN_p2 = rData.DN_v;

    // Get time integration scheme data
    const double beta0_gauss = rData.beta_0;
    const double beta1_gauss = rData.beta_1;
    const double beta2_gauss = rData.beta_2;

    // Get nodal values
    const auto& F_u_nodes = rData.MomentumForcing;
    const auto& F_u_n_nodes = rData.MomentumForcingOld1;
    const auto& F_u_nn_nodes = rData.MomentumForcingOld2;

    // Assemble RHS contribution
    const double w_g = rData.Weight;

    const double crRightHandSideVector0 = beta0_gauss*(F_u_nn_nodes(0,0)*N_p2[0] + F_u_nn_nodes(1,0)*N_p2[1] + F_u_nn_nodes(2,0)*N_p2[2] + F_u_nn_nodes(3,0)*N_p2[3] + F_u_nn_nodes(4,0)*N_p2[4] + F_u_nn_nodes(5,0)*N_p2[5]) + beta1_gauss*(F_u_n_nodes(0,0)*N_p2[0] + F_u_n_nodes(1,0)*N_p2[1] + F_u_n_nodes(2,0)*N_p2[2] + F_u_n_nodes(3,0)*N_p2[3] + F_u_n_nodes(4,0)*N_p2[4] + F_u_n_nodes(5,0)*N_p2[5]) + beta2_gauss*(F_u_nodes(0,0)*N_p2[0] + F_u_nodes(1,0)*N_p2[1] + F_u_nodes(2,0)*N_p2[2] + F_u_nodes(3,0)*N_p2[3] + F_u_nodes(4,0)*N_p2[4] + F_u_nodes(5,0)*N_p2[5]);
    const double crRightHandSideVector1 = N_p2[0]*w_g;
    const double crRightHandSideVector2 = beta0_gauss*(F_u_nn_nodes(0,1)*N_p2[0] + F_u_nn_nodes(1,1)*N_p2[1] + F_u_nn_nodes(2,1)*N_p2[2] + F_u_nn_nodes(3,1)*N_p2[3] + F_u_nn_nodes(4,1)*N_p2[4] + F_u_nn_nodes(5,1)*N_p2[5]) + beta1_gauss*(F_u_n_nodes(0,1)*N_p2[0] + F_u_n_nodes(1,1)*N_p2[1] + F_u_n_nodes(2,1)*N_p2[2] + F_u_n_nodes(3,1)*N_p2[3] + F_u_n_nodes(4,1)*N_p2[4] + F_u_n_nodes(5,1)*N_p2[5]) + beta2_gauss*(F_u_nodes(0,1)*N_p2[0] + F_u_nodes(1,1)*N_p2[1] + F_u_nodes(2,1)*N_p2[2] + F_u_nodes(3,1)*N_p2[3] + F_u_nodes(4,1)*N_p2[4] + F_u_nodes(5,1)*N_p2[5]);
    const double crRightHandSideVector3 = N_p2[1]*w_g;
    const double crRightHandSideVector4 = N_p2[2]*w_g;
    const double crRightHandSideVector5 = N_p2[3]*w_g;
    const double crRightHandSideVector6 = N_p2[4]*w_g;
    const double crRightHandSideVector7 = N_p2[5]*w_g;
    rRightHandSideVector[0]+=-crRightHandSideVector0*crRightHandSideVector1;
    rRightHandSideVector[1]+=-crRightHandSideVector1*crRightHandSideVector2;
    rRightHandSideVector[2]+=0;
    rRightHandSideVector[3]+=0;
    rRightHandSideVector[4]+=0;
    rRightHandSideVector[5]+=-crRightHandSideVector0*crRightHandSideVector3;
    rRightHandSideVector[6]+=-crRightHandSideVector2*crRightHandSideVector3;
    rRightHandSideVector[7]+=0;
    rRightHandSideVector[8]+=0;
    rRightHandSideVector[9]+=0;
    rRightHandSideVector[10]+=-crRightHandSideVector0*crRightHandSideVector4;
    rRightHandSideVector[11]+=-crRightHandSideVector2*crRightHandSideVector4;
    rRightHandSideVector[12]+=0;
    rRightHandSideVector[13]+=0;
    rRightHandSideVector[14]+=0;
    rRightHandSideVector[15]+=-crRightHandSideVector0*crRightHandSideVector5;
    rRightHandSideVector[16]+=-crRightHandSideVector2*crRightHandSideVector5;
    rRightHandSideVector[17]+=0;
    rRightHandSideVector[18]+=0;
    rRightHandSideVector[19]+=0;
    rRightHandSideVector[20]+=-crRightHandSideVector0*crRightHandSideVector6;
    rRightHandSideVector[21]+=-crRightHandSideVector2*crRightHandSideVector6;
    rRightHandSideVector[22]+=0;
    rRightHandSideVector[23]+=0;
    rRightHandSideVector[24]+=0;
    rRightHandSideVector[25]+=-crRightHandSideVector0*crRightHandSideVector7;
    rRightHandSideVector[26]+=-crRightHandSideVector2*crRightHandSideVector7;
    rRightHandSideVector[27]+=0;
    rRightHandSideVector[28]+=0;
    rRightHandSideVector[29]+=0;
    rRightHandSideVector[30]+=0;
    rRightHandSideVector[31]+=0;
    rRightHandSideVector[32]+=0;
    rRightHandSideVector[33]+=0;
    rRightHandSideVector[34]+=0;
    rRightHandSideVector[35]+=0;

}

template <>
void ThermallyDrivenActiveFluid<3>::AddGaussPointRightHandSideMomentumForcingContribution(
    const ElementDataContainer& rData,
    VectorType& rRightHandSideVector)
{
    const SizeType TDim = 3;

    // Get shape function values
    const auto& N_p2 = rData.N_v;
    const auto& DN_p2 = rData.DN_v;

    // Get time integration scheme data
    const double beta0_gauss = rData.beta_0;
    const double beta1_gauss = rData.beta_1;
    const double beta2_gauss = rData.beta_2;

    // Get nodal values
    const auto& F_u_nodes = rData.MomentumForcing;
    const auto& F_u_n_nodes = rData.MomentumForcingOld1;
    const auto& F_u_nn_nodes = rData.MomentumForcingOld2;

    // Assemble RHS contribution
    const double w_g = rData.Weight;

    const double crRightHandSideVector0 = beta0_gauss*(F_u_nn_nodes(0,0)*N_p2[0] + F_u_nn_nodes(1,0)*N_p2[1] + F_u_nn_nodes(2,0)*N_p2[2] + F_u_nn_nodes(3,0)*N_p2[3] + F_u_nn_nodes(4,0)*N_p2[4] + F_u_nn_nodes(5,0)*N_p2[5] + F_u_nn_nodes(6,0)*N_p2[6] + F_u_nn_nodes(7,0)*N_p2[7] + F_u_nn_nodes(8,0)*N_p2[8] + F_u_nn_nodes(9,0)*N_p2[9]) + beta1_gauss*(F_u_n_nodes(0,0)*N_p2[0] + F_u_n_nodes(1,0)*N_p2[1] + F_u_n_nodes(2,0)*N_p2[2] + F_u_n_nodes(3,0)*N_p2[3] + F_u_n_nodes(4,0)*N_p2[4] + F_u_n_nodes(5,0)*N_p2[5] + F_u_n_nodes(6,0)*N_p2[6] + F_u_n_nodes(7,0)*N_p2[7] + F_u_n_nodes(8,0)*N_p2[8] + F_u_n_nodes(9,0)*N_p2[9]) + beta2_gauss*(F_u_nodes(0,0)*N_p2[0] + F_u_nodes(1,0)*N_p2[1] + F_u_nodes(2,0)*N_p2[2] + F_u_nodes(3,0)*N_p2[3] + F_u_nodes(4,0)*N_p2[4] + F_u_nodes(5,0)*N_p2[5] + F_u_nodes(6,0)*N_p2[6] + F_u_nodes(7,0)*N_p2[7] + F_u_nodes(8,0)*N_p2[8] + F_u_nodes(9,0)*N_p2[9]);
    const double crRightHandSideVector1 = N_p2[0]*w_g;
    const double crRightHandSideVector2 = beta0_gauss*(F_u_nn_nodes(0,1)*N_p2[0] + F_u_nn_nodes(1,1)*N_p2[1] + F_u_nn_nodes(2,1)*N_p2[2] + F_u_nn_nodes(3,1)*N_p2[3] + F_u_nn_nodes(4,1)*N_p2[4] + F_u_nn_nodes(5,1)*N_p2[5] + F_u_nn_nodes(6,1)*N_p2[6] + F_u_nn_nodes(7,1)*N_p2[7] + F_u_nn_nodes(8,1)*N_p2[8] + F_u_nn_nodes(9,1)*N_p2[9]) + beta1_gauss*(F_u_n_nodes(0,1)*N_p2[0] + F_u_n_nodes(1,1)*N_p2[1] + F_u_n_nodes(2,1)*N_p2[2] + F_u_n_nodes(3,1)*N_p2[3] + F_u_n_nodes(4,1)*N_p2[4] + F_u_n_nodes(5,1)*N_p2[5] + F_u_n_nodes(6,1)*N_p2[6] + F_u_n_nodes(7,1)*N_p2[7] + F_u_n_nodes(8,1)*N_p2[8] + F_u_n_nodes(9,1)*N_p2[9]) + beta2_gauss*(F_u_nodes(0,1)*N_p2[0] + F_u_nodes(1,1)*N_p2[1] + F_u_nodes(2,1)*N_p2[2] + F_u_nodes(3,1)*N_p2[3] + F_u_nodes(4,1)*N_p2[4] + F_u_nodes(5,1)*N_p2[5] + F_u_nodes(6,1)*N_p2[6] + F_u_nodes(7,1)*N_p2[7] + F_u_nodes(8,1)*N_p2[8] + F_u_nodes(9,1)*N_p2[9]);
    const double crRightHandSideVector3 = beta0_gauss*(F_u_nn_nodes(0,2)*N_p2[0] + F_u_nn_nodes(1,2)*N_p2[1] + F_u_nn_nodes(2,2)*N_p2[2] + F_u_nn_nodes(3,2)*N_p2[3] + F_u_nn_nodes(4,2)*N_p2[4] + F_u_nn_nodes(5,2)*N_p2[5] + F_u_nn_nodes(6,2)*N_p2[6] + F_u_nn_nodes(7,2)*N_p2[7] + F_u_nn_nodes(8,2)*N_p2[8] + F_u_nn_nodes(9,2)*N_p2[9]) + beta1_gauss*(F_u_n_nodes(0,2)*N_p2[0] + F_u_n_nodes(1,2)*N_p2[1] + F_u_n_nodes(2,2)*N_p2[2] + F_u_n_nodes(3,2)*N_p2[3] + F_u_n_nodes(4,2)*N_p2[4] + F_u_n_nodes(5,2)*N_p2[5] + F_u_n_nodes(6,2)*N_p2[6] + F_u_n_nodes(7,2)*N_p2[7] + F_u_n_nodes(8,2)*N_p2[8] + F_u_n_nodes(9,2)*N_p2[9]) + beta2_gauss*(F_u_nodes(0,2)*N_p2[0] + F_u_nodes(1,2)*N_p2[1] + F_u_nodes(2,2)*N_p2[2] + F_u_nodes(3,2)*N_p2[3] + F_u_nodes(4,2)*N_p2[4] + F_u_nodes(5,2)*N_p2[5] + F_u_nodes(6,2)*N_p2[6] + F_u_nodes(7,2)*N_p2[7] + F_u_nodes(8,2)*N_p2[8] + F_u_nodes(9,2)*N_p2[9]);
    const double crRightHandSideVector4 = N_p2[1]*w_g;
    const double crRightHandSideVector5 = N_p2[2]*w_g;
    const double crRightHandSideVector6 = N_p2[3]*w_g;
    const double crRightHandSideVector7 = N_p2[4]*w_g;
    const double crRightHandSideVector8 = N_p2[5]*w_g;
    const double crRightHandSideVector9 = N_p2[6]*w_g;
    const double crRightHandSideVector10 = N_p2[7]*w_g;
    const double crRightHandSideVector11 = N_p2[8]*w_g;
    const double crRightHandSideVector12 = N_p2[9]*w_g;
    rRightHandSideVector[0]+=-crRightHandSideVector0*crRightHandSideVector1;
    rRightHandSideVector[1]+=-crRightHandSideVector1*crRightHandSideVector2;
    rRightHandSideVector[2]+=-crRightHandSideVector1*crRightHandSideVector3;
    rRightHandSideVector[3]+=0;
    rRightHandSideVector[4]+=0;
    rRightHandSideVector[5]+=0;
    rRightHandSideVector[6]+=0;
    rRightHandSideVector[7]+=-crRightHandSideVector0*crRightHandSideVector4;
    rRightHandSideVector[8]+=-crRightHandSideVector2*crRightHandSideVector4;
    rRightHandSideVector[9]+=-crRightHandSideVector3*crRightHandSideVector4;
    rRightHandSideVector[10]+=0;
    rRightHandSideVector[11]+=0;
    rRightHandSideVector[12]+=0;
    rRightHandSideVector[13]+=0;
    rRightHandSideVector[14]+=-crRightHandSideVector0*crRightHandSideVector5;
    rRightHandSideVector[15]+=-crRightHandSideVector2*crRightHandSideVector5;
    rRightHandSideVector[16]+=-crRightHandSideVector3*crRightHandSideVector5;
    rRightHandSideVector[17]+=0;
    rRightHandSideVector[18]+=0;
    rRightHandSideVector[19]+=0;
    rRightHandSideVector[20]+=0;
    rRightHandSideVector[21]+=-crRightHandSideVector0*crRightHandSideVector6;
    rRightHandSideVector[22]+=-crRightHandSideVector2*crRightHandSideVector6;
    rRightHandSideVector[23]+=-crRightHandSideVector3*crRightHandSideVector6;
    rRightHandSideVector[24]+=0;
    rRightHandSideVector[25]+=0;
    rRightHandSideVector[26]+=0;
    rRightHandSideVector[27]+=0;
    rRightHandSideVector[28]+=-crRightHandSideVector0*crRightHandSideVector7;
    rRightHandSideVector[29]+=-crRightHandSideVector2*crRightHandSideVector7;
    rRightHandSideVector[30]+=-crRightHandSideVector3*crRightHandSideVector7;
    rRightHandSideVector[31]+=0;
    rRightHandSideVector[32]+=0;
    rRightHandSideVector[33]+=0;
    rRightHandSideVector[34]+=0;
    rRightHandSideVector[35]+=-crRightHandSideVector0*crRightHandSideVector8;
    rRightHandSideVector[36]+=-crRightHandSideVector2*crRightHandSideVector8;
    rRightHandSideVector[37]+=-crRightHandSideVector3*crRightHandSideVector8;
    rRightHandSideVector[38]+=0;
    rRightHandSideVector[39]+=0;
    rRightHandSideVector[40]+=0;
    rRightHandSideVector[41]+=0;
    rRightHandSideVector[42]+=-crRightHandSideVector0*crRightHandSideVector9;
    rRightHandSideVector[43]+=-crRightHandSideVector2*crRightHandSideVector9;
    rRightHandSideVector[44]+=-crRightHandSideVector3*crRightHandSideVector9;
    rRightHandSideVector[45]+=0;
    rRightHandSideVector[46]+=0;
    rRightHandSideVector[47]+=0;
    rRightHandSideVector[48]+=0;
    rRightHandSideVector[49]+=-crRightHandSideVector0*crRightHandSideVector10;
    rRightHandSideVector[50]+=-crRightHandSideVector10*crRightHandSideVector2;
    rRightHandSideVector[51]+=-crRightHandSideVector10*crRightHandSideVector3;
    rRightHandSideVector[52]+=0;
    rRightHandSideVector[53]+=0;
    rRightHandSideVector[54]+=0;
    rRightHandSideVector[55]+=0;
    rRightHandSideVector[56]+=-crRightHandSideVector0*crRightHandSideVector11;
    rRightHandSideVector[57]+=-crRightHandSideVector11*crRightHandSideVector2;
    rRightHandSideVector[58]+=-crRightHandSideVector11*crRightHandSideVector3;
    rRightHandSideVector[59]+=0;
    rRightHandSideVector[60]+=0;
    rRightHandSideVector[61]+=0;
    rRightHandSideVector[62]+=0;
    rRightHandSideVector[63]+=-crRightHandSideVector0*crRightHandSideVector12;
    rRightHandSideVector[64]+=-crRightHandSideVector12*crRightHandSideVector2;
    rRightHandSideVector[65]+=-crRightHandSideVector12*crRightHandSideVector3;
    rRightHandSideVector[66]+=0;
    rRightHandSideVector[67]+=0;
    rRightHandSideVector[68]+=0;
    rRightHandSideVector[69]+=0;
    rRightHandSideVector[70]+=0;
    rRightHandSideVector[71]+=0;
    rRightHandSideVector[72]+=0;
    rRightHandSideVector[73]+=0;
    rRightHandSideVector[74]+=0;
    rRightHandSideVector[75]+=0;
    rRightHandSideVector[76]+=0;
    rRightHandSideVector[77]+=0;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private serialization

template< unsigned int TDim >
void ThermallyDrivenActiveFluid<TDim>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}


template< unsigned int TDim >
void ThermallyDrivenActiveFluid<TDim>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation

template class ThermallyDrivenActiveFluid<2>;
template class ThermallyDrivenActiveFluid<3>;

}