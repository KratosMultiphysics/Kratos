//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

// System includes


// External includes


// Project includes
#include "includes/checks.h"

// Application includes
#include "navier_stokes_p2_p1_continuous_wall_condition.h"
#include "wall_laws/linear_log_wall_law.h"
#include "wall_laws/navier_slip_wall_law.h"

namespace Kratos
{

template<unsigned int TDim, class... TWallModel>
void NavierStokesP2P1ContinuousWallCondition<TDim, TWallModel...>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rResult.size() != LocalSize) {
        rResult.resize(LocalSize, false);
    }

    IndexType local_index = 0;
    const auto& r_geometry = this->GetGeometry();
    const IndexType x_pos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    for (IndexType i = 0; i < VelocityNumNodes; ++i) {
        rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_X, x_pos).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_Y, x_pos+1).EquationId();
        if constexpr (TDim == 3) {
            rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_Z, x_pos+2).EquationId();
        }
    }

    const IndexType p_pos = this->GetGeometry()[0].GetDofPosition(PRESSURE);
    for (IndexType i = 0; i < PressureNumNodes; ++i) {
        rResult[local_index++] = r_geometry[i].GetDof(PRESSURE, p_pos).EquationId();
    }
}

template<unsigned int TDim, class... TWallModel>
void NavierStokesP2P1ContinuousWallCondition<TDim, TWallModel...>::GetDofList(
    DofsVectorType& rConditionDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rConditionDofList.size() != LocalSize) {
        rConditionDofList.resize(LocalSize);
    }

    IndexType local_index = 0;
    const auto& r_geometry = this->GetGeometry();
    const IndexType x_pos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    for (IndexType i = 0; i < VelocityNumNodes; ++i) {
        rConditionDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_X, x_pos);
        rConditionDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_Y, x_pos+1);
        if constexpr (TDim == 3) {
            rConditionDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_Z, x_pos+2);
        }
    }

    const IndexType p_pos = this->GetGeometry()[0].GetDofPosition(PRESSURE);
    for (IndexType i = 0; i < PressureNumNodes; ++i) {
        rConditionDofList[local_index++] = r_geometry[i].pGetDof(PRESSURE, p_pos);
    }
}

template<unsigned int TDim, class... TWallModel>
void NavierStokesP2P1ContinuousWallCondition<TDim, TWallModel...>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check (and resize) RHS and LHS arrays
    if (rLeftHandSideMatrix.size1() != LocalSize) {
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false); //false says not to preserve existing storage!!
    }
    if (rRightHandSideVector.size() != LocalSize) {
        rRightHandSideVector.resize(LocalSize, false); //false says not to preserve existing storage!!
    }

    // Struct to pass around the data
    ConditionDataStruct data;

    // LHS and RHS contributions initialization
    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    // Compute condition unit normal vector
    this->CalculateUnitNormal(data.UnitNormal);

    // Gauss point information
    auto& r_geom = this->GetGeometry();
    const auto& r_integration_points = r_geom.IntegrationPoints(IntegrationMethod);
    const SizeType num_gauss = r_integration_points.size();
    Vector gauss_pts_det = ZeroVector(num_gauss);
    r_geom.DeterminantOfJacobian(gauss_pts_det, IntegrationMethod);
    const MatrixType N_v_container = r_geom.ShapeFunctionsValues(IntegrationMethod);

    // Calculate viscous stress for the slip tangential correction
    if (rCurrentProcessInfo.Has(SLIP_TANGENTIAL_CORRECTION_SWITCH)) {
        if (this->Is(SLIP) && rCurrentProcessInfo.GetValue(SLIP_TANGENTIAL_CORRECTION_SWITCH)) {
            KRATOS_WARNING("NavierStokesP2P1ContinuousWallCondition")
                << "Slip tangential correction is not implemented. Please switch off 'slip_tangential_correction' in 'ApplySlipProcess'." << std::endl;
        }
    }

    // Loop on gauss points
    for(IndexType i_gauss = 0; i_gauss< num_gauss; ++i_gauss) {
        data.N_v = row(N_v_container, i_gauss);
        data.Weight = gauss_pts_det[i_gauss] * r_integration_points[i_gauss].Weight();
        AddGaussPointRHSContribution(rRightHandSideVector, data, rCurrentProcessInfo);
    }

    // Add the wall law contribution
    constexpr SizeType n_wall_models = sizeof...(TWallModel);
    static_assert(n_wall_models < 2, "More than one template wall model argument in 'NavierStokesWallCondition'.");
    if (this->Is(WALL) && n_wall_models != 0) {
        (AddWallModelLocalSystemCall<TWallModel>(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo), ...);
    }

    KRATOS_CATCH("")
}

template<unsigned int TDim, class... TWallModel>
void NavierStokesP2P1ContinuousWallCondition<TDim, TWallModel...>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check (and resize) LHS matrix
    if (rLeftHandSideMatrix.size1() != LocalSize) {
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false); //false says not to preserve existing storage!!
    }

    // LHS contribution initialization
    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);

    // Add the wall law contribution
    constexpr SizeType n_wall_models = sizeof...(TWallModel);
    static_assert(n_wall_models < 2, "More than one template wall model argument in 'NavierStokesWallCondition'.");
    if (this->Is(WALL) && n_wall_models != 0) {
        (AddWallModelLeftHandSideCall<TWallModel>(rLeftHandSideMatrix, rCurrentProcessInfo), ...);
    }

    KRATOS_CATCH("")
}

template<unsigned int TDim, class... TWallModel>
void NavierStokesP2P1ContinuousWallCondition<TDim, TWallModel...>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Check (and resize) RHS vector
    if (rRightHandSideVector.size() != LocalSize) {
        rRightHandSideVector.resize(LocalSize, false); //false says not to preserve existing storage!!
    }

    // Struct to pass around the data
    ConditionDataStruct data;

    // RHS contribution initialization
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    // Compute condition unit normal vector
    this->CalculateUnitNormal(data.UnitNormal);

    // Gauss point information
    auto& r_geom = this->GetGeometry();
    const auto& r_integration_points = r_geom.IntegrationPoints(IntegrationMethod);
    const SizeType num_gauss = r_integration_points.size();
    Vector gauss_pts_det = ZeroVector(num_gauss);
    r_geom.DeterminantOfJacobian(gauss_pts_det, IntegrationMethod);
    const MatrixType N_v_container = r_geom.ShapeFunctionsValues(IntegrationMethod);

    // Calculate viscous stress for the slip tangential correction
    if (rCurrentProcessInfo.Has(SLIP_TANGENTIAL_CORRECTION_SWITCH)) {
        if (this->Is(SLIP) && rCurrentProcessInfo.GetValue(SLIP_TANGENTIAL_CORRECTION_SWITCH)) {
            KRATOS_WARNING("NavierStokesP2P1ContinuousWallCondition")
                << "Slip tangential correction is not implemented. Please switch off 'slip_tangential_correction' in 'ApplySlipProcess'." << std::endl;
        }
    }

    // Loop on gauss points
    for(IndexType i_gauss = 0; i_gauss< num_gauss; ++i_gauss) {
        data.N_v = row(N_v_container, i_gauss);
        data.Weight = gauss_pts_det[i_gauss] * r_integration_points[i_gauss].Weight();
        AddGaussPointRHSContribution(rRightHandSideVector, data, rCurrentProcessInfo);
    }

    // Add the wall law contribution
    constexpr SizeType n_wall_models = sizeof...(TWallModel);
    static_assert(n_wall_models < 2, "More than one template wall model argument in 'NavierStokesWallCondition'.");
    if (this->Is(WALL) && n_wall_models != 0) {
        (AddWallModelRightHandSideCall<TWallModel>(rRightHandSideVector, rCurrentProcessInfo), ...);
    }

    KRATOS_CATCH("")
}

template<unsigned int TDim, class... TWallModel>
int NavierStokesP2P1ContinuousWallCondition<TDim, TWallModel...>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;
    int check = BaseType::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0
    if (check != 0) {
        return check;
    } else {
        // Check that geometry is coplanar (i.e. edges midpoint nodes are aligned)
        // Note that though Kratos geometry supports non coplanar, this is assumed througout current implementation
        const auto& r_geometry = this->GetGeometry();
        if constexpr (TDim == 2) {
            array_1d<double, 3> vect_01 = r_geometry[1].Coordinates() - r_geometry[0].Coordinates();
            vect_01 /= norm_2(vect_01);
            array_1d<double, 3> vect_02 = r_geometry[1].Coordinates() - r_geometry[0].Coordinates();
            vect_02 /= norm_2(vect_02);
            KRATOS_CHECK_VECTOR_NEAR(vect_01, vect_02, 1.0e-12)
        } else {
            array_1d<double, 3> vect_01 = r_geometry[1].Coordinates() - r_geometry[0].Coordinates();
            vect_01 /= norm_2(vect_01);
            array_1d<double, 3> vect_03 = r_geometry[3].Coordinates() - r_geometry[0].Coordinates();
            vect_03 /= norm_2(vect_03);
            KRATOS_CHECK_VECTOR_NEAR(vect_01, vect_03, 1.0e-12)
            array_1d<double, 3> vect_02 = r_geometry[2].Coordinates() - r_geometry[0].Coordinates();
            vect_02 /= norm_2(vect_02);
            array_1d<double, 3> vect_05 = r_geometry[5].Coordinates() - r_geometry[0].Coordinates();
            vect_05 /= norm_2(vect_05);
            KRATOS_CHECK_VECTOR_NEAR(vect_02, vect_05, 1.0e-12)
            array_1d<double, 3> vect_12 = r_geometry[2].Coordinates() - r_geometry[1].Coordinates();
            vect_12 /= norm_2(vect_12);
            array_1d<double, 3> vect_14 = r_geometry[4].Coordinates() - r_geometry[1].Coordinates();
            vect_14 /= norm_2(vect_14);
            KRATOS_CHECK_VECTOR_NEAR(vect_12, vect_14, 1.0e-12)
        }

        // Checks on nodes
        // Check that the element's nodes contain all required SolutionStepData and Degrees Of Freedom variables
        for (const auto& r_node : r_geometry) {
            // Check variables
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY, r_node)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(PRESSURE, r_node)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(MESH_VELOCITY, r_node)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(EXTERNAL_PRESSURE, r_node)
            // Check DOFs
            KRATOS_CHECK_DOF_IN_NODE(VELOCITY_X, r_node)
            KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Y, r_node)
            KRATOS_CHECK_DOF_IN_NODE(VELOCITY_Z, r_node)
            KRATOS_CHECK_DOF_IN_NODE(PRESSURE, r_node)
        }

        // Check that parents have been computed
        // These are required to retrieve the material properties and the viscous stress
        auto& parent_elements = this->GetValue(NEIGHBOUR_ELEMENTS);
        KRATOS_ERROR_IF(parent_elements.size() > 1) << "Condition " << this->Id() << " was assigned more than one parent element." << std::endl;
        KRATOS_ERROR_IF(parent_elements.size() == 0) << "Condition " << this->Id() << " has no parent element. Please execute 'check_and_prepare_model_process_fluid' process." << std::endl;

        // If provided, check wall law
        constexpr SizeType n_wall_models = sizeof...(TWallModel);
        static_assert(n_wall_models < 2, "More than one template wall model argument in 'NavierStokesWallCondition'.");
        if constexpr (n_wall_models != 0) {
            ((check = WallModelCheckCall<TWallModel>(rCurrentProcessInfo)), ...);
        }

        return check;
    }

    KRATOS_CATCH("");
}

template<unsigned int TDim, class... TWallModel>
void NavierStokesP2P1ContinuousWallCondition<TDim, TWallModel...>::Calculate(
    const Variable< array_1d<double,3> >& rVariable,
    array_1d<double,3>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    rOutput = ZeroVector(3);

    if (rVariable == DRAG_FORCE) {
        KRATOS_ERROR << "'DRAG_FORCE' variable is not implemented for NavierStokesP2P1ContinuousWallCondition" << TDim << "D." << std::endl;
    } else {
        BaseType::Calculate(rVariable, rOutput, rCurrentProcessInfo);
    }
}

template<unsigned int TDim, class... TWallModel>
void NavierStokesP2P1ContinuousWallCondition<TDim, TWallModel...>::CalculateUnitNormal(array_1d<double,TDim>& rUnitNormal)
{
    const auto& r_geom = GetGeometry();
    if constexpr (TDim == 2) {
        rUnitNormal[0] = r_geom[1].Y() - r_geom[0].Y();
        rUnitNormal[1] = - (r_geom[1].X() - r_geom[0].X());
    } else if constexpr (TDim == 3) {
        array_1d<double,3> v1,v2;
        v1[0] = r_geom[1].X() - r_geom[0].X();
        v1[1] = r_geom[1].Y() - r_geom[0].Y();
        v1[2] = r_geom[1].Z() - r_geom[0].Z();

        v2[0] = r_geom[2].X() - r_geom[0].X();
        v2[1] = r_geom[2].Y() - r_geom[0].Y();
        v2[2] = r_geom[2].Z() - r_geom[0].Z();

        MathUtils<double>::CrossProduct(rUnitNormal,v1,v2);
        rUnitNormal *= 0.5;
    } else {
        KRATOS_ERROR << "'CalculateUnitNormal' is not implemented for current geometry." << std::endl;
    }
    rUnitNormal /= norm_2(rUnitNormal);
}

template<unsigned int TDim, class... TWallModel>
void NavierStokesP2P1ContinuousWallCondition<TDim, TWallModel...>::AddGaussPointRHSContribution(
    VectorType& rRHS,
    const ConditionDataStruct& rData,
    const ProcessInfo& rProcessInfo)
{
    // Gauss pt. Neumann BC contribution
    this->ComputeRHSNeumannContribution(rRHS, rData);

    // Gauss pt. outlet inflow prevention contribution
    if (rProcessInfo.Has(OUTLET_INFLOW_CONTRIBUTION_SWITCH)) {
        if (this->Is(OUTLET) && rProcessInfo[OUTLET_INFLOW_CONTRIBUTION_SWITCH]){
            this->ComputeRHSOutletInflowContribution(rRHS, rData, rProcessInfo);
        }
    }
}

template<unsigned int TDim, class... TWallModel>
void NavierStokesP2P1ContinuousWallCondition<TDim, TWallModel...>::ComputeRHSNeumannContribution(
    VectorType& rRHS,
    const ConditionDataStruct& rData)
{
    const auto& r_geom = this->GetGeometry();
    for (IndexType i = 0; i < VelocityNumNodes; ++i) {
        const double p_ext = r_geom[i].FastGetSolutionStepValue(EXTERNAL_PRESSURE);
        for (IndexType j = 0; j < VelocityNumNodes; ++j) {
            for (IndexType d = 0; d < TDim; ++d) {
                rRHS[j*TDim + d] -= rData.Weight * rData.N_v[j] * rData.N_v[i] * p_ext * rData.UnitNormal[d];
            }
        }
    }
}

template<unsigned int TDim, class... TWallModel>
void NavierStokesP2P1ContinuousWallCondition<TDim, TWallModel...>::ComputeRHSOutletInflowContribution(
    VectorType& rRHS,
    const ConditionDataStruct& rData,
    const ProcessInfo& rProcessInfo)
{
    // Get DENSITY from parent element properties
    const auto& r_neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);
    const double rho = r_neighbours[0].GetProperties().GetValue(DENSITY);

    // Compute Gauss velocity norm and velocity projection
    const auto& r_geom = this->GetGeometry();
    array_1d<double, 3> v_gauss = ZeroVector(3);
    for (IndexType i = 0; i < VelocityNumNodes; ++i) {
        const auto& r_v = r_geom[i].FastGetSolutionStepValue(VELOCITY);
        v_gauss += rData.N_v[i] * r_v;
    }
    const double v_gauss_proj = std::inner_product(rData.UnitNormal.begin(), rData.UnitNormal.end(), v_gauss.begin(), 0.0);
    const double v_gauss_squared_norm = std::pow(v_gauss[0],2) + std::pow(v_gauss[1],2) + std::pow(v_gauss[2],2);

    // Add outlet inflow prevention contribution
    const double delta = 1.0e-2;
    const double U_0 = rProcessInfo[CHARACTERISTIC_VELOCITY];
    const double S_0 = 0.5*(1.0 - std::tanh(v_gauss_proj/(U_0*delta)));
    const double aux = rData.Weight * 0.5 * rho * v_gauss_squared_norm * S_0;
    for (IndexType i = 0; i < VelocityNumNodes; ++i) {
        for (IndexType d = 0; d < TDim; ++d) {
            rRHS[i*TDim + d] += aux * rData.N_v[i] * rData.UnitNormal[d];
        }
    }
}

template class NavierStokesP2P1ContinuousWallCondition<2>;
template class NavierStokesP2P1ContinuousWallCondition<3>;

} // namespace Kratos
