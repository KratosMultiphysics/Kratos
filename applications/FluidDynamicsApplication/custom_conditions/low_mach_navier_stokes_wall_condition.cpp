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
#include "low_mach_navier_stokes_wall_condition.h"

namespace Kratos
{

///@name Specialized implementation of VMS for functions that depend on TDim
///@{

template<unsigned int TDim, unsigned int TNumNodes, class... TWallModel>
void LowMachNavierStokesWallCondition<TDim,TNumNodes,TWallModel...>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_geometry = this->GetGeometry();

    if (rResult.size() != LocalSize) {
        rResult.resize(LocalSize, false);
    }

    const unsigned int p_pos = this->GetGeometry()[0].GetDofPosition(PRESSURE);
    const unsigned int x_pos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    const unsigned int t_pos = this->GetGeometry()[0].GetDofPosition(TEMPERATURE);

    unsigned int local_index = 0;
    for (unsigned int i = 0; i < NumNodes; ++i) {
        rResult[local_index++] = r_geometry[i].GetDof(PRESSURE, p_pos).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_X, x_pos).EquationId();
        rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_Y, x_pos+1).EquationId();
        if constexpr (Dim == 3) {
            rResult[local_index++] = r_geometry[i].GetDof(VELOCITY_Z, x_pos+2).EquationId();
        }
        rResult[local_index++] = r_geometry[i].GetDof(TEMPERATURE, t_pos).EquationId();
    }
}

template<unsigned int TDim, unsigned int TNumNodes, class... TWallModel>
void LowMachNavierStokesWallCondition<TDim,TNumNodes,TWallModel...>::GetDofList(
    DofsVectorType& rConditionDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const auto& r_geometry = this->GetGeometry();

    if (rConditionDofList.size() != LocalSize) {
         rConditionDofList.resize(LocalSize);
    }

    const unsigned int p_pos = this->GetGeometry()[0].GetDofPosition(PRESSURE);
    const unsigned int x_pos = this->GetGeometry()[0].GetDofPosition(VELOCITY_X);
    const unsigned int t_pos = this->GetGeometry()[0].GetDofPosition(TEMPERATURE);

    unsigned int local_index = 0;
    for (unsigned int i = 0; i < NumNodes; ++i) {
        rConditionDofList[local_index++] = r_geometry[i].pGetDof(PRESSURE, p_pos);
        rConditionDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_X, x_pos);
        rConditionDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_Y, x_pos+1);
        if constexpr (Dim == 3) {
            rConditionDofList[local_index++] = r_geometry[i].pGetDof(VELOCITY_Z,x_pos+2);
        }
        rConditionDofList[local_index++] = r_geometry[i].pGetDof(TEMPERATURE, t_pos);
    }
}

template<unsigned int TDim, unsigned int TNumNodes, class... TWallModel>
void LowMachNavierStokesWallCondition<TDim,TNumNodes,TWallModel...>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Resize if needed
    if (rLeftHandSideMatrix.size1() != LocalSize || rLeftHandSideMatrix.size2() != LocalSize) {
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
    }
    if (rRightHandSideVector.size() != LocalSize) {
        rRightHandSideVector.resize(LocalSize, false);
    }

    // LHS and RHS contributions initialization
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);

    // // Struct to pass around the data
    // ConditionDataStruct data;

    // // Allocate memory needed
    // array_1d<double,MatrixSize> rhs_gauss;
    // BoundedMatrix<double,MatrixSize, MatrixSize> lhs_gauss;

    // // Compute condition unit normal vector
    // this->CalculateNormal(data.Normal); //this already contains the area
    // const double A = norm_2(data.Normal);
    // data.Normal /= A;

    // // Gauss point information
    // GeometryType& rGeom = this->GetGeometry();
    // const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);
    // const unsigned int NumGauss = IntegrationPoints.size();
    // Vector GaussPtsJDet = ZeroVector(NumGauss);
    // rGeom.DeterminantOfJacobian(GaussPtsJDet, GeometryData::IntegrationMethod::GI_GAUSS_2);
    // const MatrixType Ncontainer = rGeom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);

    // // Calculate viscous stress for the slip tangential correction
    // if (rCurrentProcessInfo.Has(SLIP_TANGENTIAL_CORRECTION_SWITCH)) {
    //     if (this->Is(SLIP) && rCurrentProcessInfo.GetValue(SLIP_TANGENTIAL_CORRECTION_SWITCH)) {
    //         // Finding parent element to retrieve viscous stresses which are later stored in "data"
    //         auto& r_parent = this->GetValue(NEIGHBOUR_ELEMENTS)[0];
    //         data.ViscousStress = ZeroVector(VoigtSize);
    //         r_parent.Calculate(FLUID_STRESS, data.ViscousStress, rCurrentProcessInfo);
    //     }
    // }

    // // Loop on gauss points
    // for(unsigned int igauss = 0; igauss<NumGauss; igauss++)
    // {
    //     data.N = row(Ncontainer, igauss);
    //     const double J = GaussPtsJDet[igauss];
    //     data.wGauss = J * IntegrationPoints[igauss].Weight();
    //     ComputeGaussPointRHSContribution(rhs_gauss, data, rCurrentProcessInfo);
    //     ComputeGaussPointLHSContribution(lhs_gauss, data, rCurrentProcessInfo);
    //     //TODO: Implement a ComputeGausPointLocalSystemContribution
    //     noalias(rLeftHandSideMatrix) += lhs_gauss;
    //     noalias(rRightHandSideVector) += rhs_gauss;
    // }

    // // Add the wall law contribution
    // constexpr SizeType n_wall_models = sizeof...(TWallModel);
    // static_assert(n_wall_models < 2, "More than one template wall model argument in 'LowMachNavierStokesWallCondition'.");
    // if (this->Is(WALL) && n_wall_models != 0) {
    //     (AddWallModelLocalSystemCall<TWallModel>(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo), ...);
    // }

    KRATOS_CATCH("")
}


template<unsigned int TDim, unsigned int TNumNodes, class... TWallModel>
void LowMachNavierStokesWallCondition<TDim,TNumNodes,TWallModel...>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Resize if needed
    if (rLeftHandSideMatrix.size1() != LocalSize || rLeftHandSideMatrix.size2() != LocalSize) {
        rLeftHandSideMatrix.resize(LocalSize, LocalSize, false);
    }

    // LHS contributions initialization
    noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize, LocalSize);

    // // Add the wall law contribution
    // constexpr SizeType n_wall_models = sizeof...(TWallModel);
    // static_assert(n_wall_models < 2, "More than one template wall model argument in 'LowMachNavierStokesWallCondition'.");
    // if (this->Is(WALL) && n_wall_models != 0) {
    //     (AddWallModelLeftHandSideCall<TWallModel>(rLeftHandSideMatrix, rCurrentProcessInfo), ...);
    // }

    KRATOS_CATCH("")
}



template<unsigned int TDim, unsigned int TNumNodes, class... TWallModel>
void LowMachNavierStokesWallCondition<TDim,TNumNodes,TWallModel...>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Resize if needed
    if (rRightHandSideVector.size() != LocalSize) {
        rRightHandSideVector.resize(LocalSize, false);
    }

    // Loop on gauss points
    noalias(rRightHandSideVector) = ZeroVector(LocalSize);

    // // Struct to pass around the data
    // ConditionDataStruct data;
    // // Allocate memory needed
    // array_1d<double,MatrixSize> rhs_gauss;

    // // Compute condition normal
    // this->CalculateNormal(data.Normal); //this already contains the area
    // const double A = norm_2(data.Normal);
    // data.Normal /= A;

    // // Gauss point information
    // GeometryType& rGeom = this->GetGeometry();
    // const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);
    // const unsigned int NumGauss = IntegrationPoints.size();
    // Vector GaussPtsJDet = ZeroVector(NumGauss);
    // rGeom.DeterminantOfJacobian(GaussPtsJDet, GeometryData::IntegrationMethod::GI_GAUSS_2);
    // const MatrixType Ncontainer = rGeom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);

    // // Calculate viscous stress for the slip tangential correction
    // if (rCurrentProcessInfo.Has(SLIP_TANGENTIAL_CORRECTION_SWITCH)) {
    //     if (this->Is(SLIP) && rCurrentProcessInfo.GetValue(SLIP_TANGENTIAL_CORRECTION_SWITCH)) {
    //         // Finding parent element to retrieve viscous stresses which are later stored in "data"
    //         auto& r_parent = this->GetValue(NEIGHBOUR_ELEMENTS)[0];
    //         data.ViscousStress = ZeroVector(VoigtSize);
    //         r_parent.Calculate(FLUID_STRESS, data.ViscousStress, rCurrentProcessInfo);
    //     }
    // }

    // for(unsigned int igauss = 0; igauss<NumGauss; igauss++)
    // {
    //     data.N = row(Ncontainer, igauss);
    //     const double J = GaussPtsJDet[igauss];
    //     data.wGauss = J * IntegrationPoints[igauss].Weight();
    //     ComputeGaussPointRHSContribution(rhs_gauss, data, rCurrentProcessInfo);
    //     noalias(rRightHandSideVector) += rhs_gauss;
    // }

    // // Add the wall law contribution
    // constexpr SizeType n_wall_models = sizeof...(TWallModel);
    // static_assert(n_wall_models < 2, "More than one template wall model argument in 'LowMachNavierStokesWallCondition'.");
    // if (this->Is(WALL) && n_wall_models != 0) {
    //     (AddWallModelRightHandSideCall<TWallModel>(rRightHandSideVector, rCurrentProcessInfo), ...);
    // }

    KRATOS_CATCH("")
}


/// Condition check
/**
 * @param rCurrentProcessInfo reference to the ProcessInfo
 */
template<unsigned int TDim, unsigned int TNumNodes, class... TWallModel>
int LowMachNavierStokesWallCondition<TDim,TNumNodes,TWallModel...>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;
    int check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0
    if (check != 0) {
        return check;
    } else {
        // Checks on nodes
        // Check that the element's nodes contain all required SolutionStepData and Degrees Of Freedom variables
        for (const auto& r_node : this->GetGeometry()) {
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
        static_assert(n_wall_models < 2, "More than one template wall model argument in 'LowMachNavierStokesWallCondition'.");
        if constexpr (n_wall_models != 0) {
            ((check = WallModelCheckCall<TWallModel>(rCurrentProcessInfo)), ...);
        }

        return check;
    }

    KRATOS_CATCH("");
}

template<unsigned int TDim, unsigned int TNumNodes, class... TWallModel>
void LowMachNavierStokesWallCondition<TDim, TNumNodes,TWallModel...>::Calculate(
    const Variable< array_1d<double,3> >& rVariable,
    array_1d<double,3>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    rOutput = ZeroVector(3);

    if (rVariable == DRAG_FORCE) {
        // const auto& r_geom = GetGeometry();
        // const auto& r_integration_points = r_geom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);
        // unsigned int n_gauss = r_integration_points.size();
        // Vector det_J_vect = ZeroVector(n_gauss);
        // r_geom.DeterminantOfJacobian(det_J_vect, GeometryData::IntegrationMethod::GI_GAUSS_2);
        // const auto N_container = r_geom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);

        // // Calculate normal
        // array_1d<double,3> normal;
        // CalculateNormal(normal);
        // normal /= norm_2(normal);

        // // Finding parent element to retrieve viscous stresses
        // // Note that we assume in here that the shear stress is constant inside the element
        // auto& r_neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);
        // KRATOS_ERROR_IF(r_neighbours.size() > 1) << "A condition was assigned more than one parent element." << std::endl;
        // KRATOS_ERROR_IF(r_neighbours.size() == 0) << "A condition was NOT assigned a parent element. "
        // << "This leads to errors for the slip condition [BEHR2004] "
        // << "Please execute the check_and_prepare_model_process_fluid process." << std::endl;

        // auto& r_parent = r_neighbours[0];
        // Vector shear_stress;
        // r_parent.Calculate(FLUID_STRESS, shear_stress, rCurrentProcessInfo);
        // array_1d<double,3> shear_stress_n;
        // ProjectViscousStress(shear_stress, normal, shear_stress_n);

        // // Loop the Gauss pts
        // for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
        //     const double w = det_J_vect[i_gauss] * r_integration_points[i_gauss].Weight();
        //     double p = 0.0;
        //     const auto& r_N = row(N_container, i_gauss);
        //     for (unsigned int i_node = 0; i_node < r_geom.PointsNumber(); ++i_node) {
        //         p += r_N[i_node] * r_geom[i_node].FastGetSolutionStepValue(PRESSURE);
        //     }
        //     rOutput += w * (p * normal - shear_stress_n);
        // }
    } else {
        Condition::Calculate(rVariable, rOutput, rCurrentProcessInfo);
    }
}

template<unsigned int TDim, unsigned int TNumNodes, class... TWallModel>
void LowMachNavierStokesWallCondition<TDim,TNumNodes,TWallModel...>::CalculateNormal(array_1d<double,3>& rAreaNormal)
{
    const auto& r_geom = GetGeometry();
    if constexpr (TDim == 2) {
        rAreaNormal[0] = r_geom[1].Y() - r_geom[0].Y();
        rAreaNormal[1] = - (r_geom[1].X() - r_geom[0].X());
        rAreaNormal[2] = 0.0;
    } else if constexpr (TDim == 3 && TNumNodes == 3) {
        array_1d<double,3> v1,v2;
        v1[0] = r_geom[1].X() - r_geom[0].X();
        v1[1] = r_geom[1].Y() - r_geom[0].Y();
        v1[2] = r_geom[1].Z() - r_geom[0].Z();

        v2[0] = r_geom[2].X() - r_geom[0].X();
        v2[1] = r_geom[2].Y() - r_geom[0].Y();
        v2[2] = r_geom[2].Z() - r_geom[0].Z();

        MathUtils<double>::CrossProduct(rAreaNormal,v1,v2);
        rAreaNormal *= 0.5;
    } else {
        KRATOS_ERROR << "'CalculateNormal' is not implemented for current geometry." << std::endl;
    }
}

// template<unsigned int TDim, unsigned int TNumNodes, class... TWallModel>
// void LowMachNavierStokesWallCondition<TDim,TNumNodes,TWallModel...>::ComputeGaussPointLHSContribution(
//     BoundedMatrix<double, LocalSize, LocalSize>& lhs_gauss,
//     const ConditionDataStruct& data,
//     const ProcessInfo& rProcessInfo)
// {
//     const unsigned int LocalSize = TDim+1;
//     lhs_gauss = ZeroMatrix(TNumNodes*LocalSize, TNumNodes*LocalSize);

//     // Contribution to avoid spurious tangential components in the pure-slip residual
//     if (rProcessInfo.Has(SLIP_TANGENTIAL_CORRECTION_SWITCH)) {
//         if (this->Is(SLIP) && rProcessInfo.GetValue(SLIP_TANGENTIAL_CORRECTION_SWITCH)) {
//             CalculateGaussPointSlipTangentialCorrectionLHSContribution(lhs_gauss, data);
//         }
//     }
// }



// template<unsigned int TDim, unsigned int TNumNodes, class... TWallModel>
// void LowMachNavierStokesWallCondition<TDim,TNumNodes,TWallModel...>::ComputeGaussPointRHSContribution(
//     array_1d<double, LocalSize>& rhs_gauss,
//     const ConditionDataStruct& data,
//     const ProcessInfo& rProcessInfo)
// {
//     // Initialize the local RHS
//     const unsigned int LocalSize = TDim+1;
//     noalias(rhs_gauss) = ZeroVector(TNumNodes*LocalSize);

//     // Gauss pt. Neumann BC contribution
//     this->ComputeRHSNeumannContribution(rhs_gauss, data);

//     // Gauss pt. outlet inflow prevention contribution
//     if (rProcessInfo.Has(OUTLET_INFLOW_CONTRIBUTION_SWITCH)) {
//         if (this->Is(OUTLET) && rProcessInfo[OUTLET_INFLOW_CONTRIBUTION_SWITCH]){
//             this->ComputeRHSOutletInflowContribution(rhs_gauss, data, rProcessInfo);
//         }
//     }

//     // Contribution to avoid spurious tangential components in the pure-slip residual
//     if (rProcessInfo.Has(SLIP_TANGENTIAL_CORRECTION_SWITCH)) {
//         if (this->Is(SLIP) && rProcessInfo[SLIP_TANGENTIAL_CORRECTION_SWITCH]) {
//             CalculateGaussPointSlipTangentialCorrectionRHSContribution(rhs_gauss, data);
//         }
//     }
// }



// template<unsigned int TDim, unsigned int TNumNodes, class... TWallModel>
// void LowMachNavierStokesWallCondition<TDim,TNumNodes,TWallModel...>::ComputeRHSNeumannContribution(
//     array_1d<double,LocalSize>& rhs_gauss,
//     const ConditionDataStruct& data)
// {
//     const unsigned int LocalSize = TDim+1;
//     const GeometryType& rGeom = this->GetGeometry();

//     // Add Neumann BC contribution
//     for (unsigned int i=0; i<TNumNodes; ++i)
//     {
//         const double pext = rGeom[i].FastGetSolutionStepValue(EXTERNAL_PRESSURE);

//         for (unsigned int j=0; j<TNumNodes; ++j)
//         {
//             unsigned int row = j*LocalSize;
//             for (unsigned int d=0; d<TDim; ++d)
//             {
//                 rhs_gauss[row+d] -= data.wGauss*data.N[j]*data.N[i]*pext*data.Normal[d];
//             }
//         }
//     }
// }


// template<unsigned int TDim, unsigned int TNumNodes, class... TWallModel>
// void LowMachNavierStokesWallCondition<TDim,TNumNodes,TWallModel...>::ComputeRHSOutletInflowContribution(
//     array_1d<double,LocalSize>& rhs_gauss,
//     const ConditionDataStruct& data,
//     const ProcessInfo& rProcessInfo)
// {
//     constexpr SizeType LocalSize = TDim+1;
//     const GeometryType& rGeom = this->GetGeometry();

//     // Get DENSITY from parent element properties
//     auto & r_neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);
//     const double rho = r_neighbours[0].GetProperties().GetValue(DENSITY);

//     // Compute Gauss pt. density, velocity norm and velocity projection
//     array_1d<double, 3> vGauss = ZeroVector(3);
//     for (unsigned int i=0; i<TNumNodes; ++i)
//     {
//         const array_1d<double, 3>& rVelNode = rGeom[i].FastGetSolutionStepValue(VELOCITY);
//         vGauss += data.N[i]*rVelNode;
//     }

//     const double vGaussProj = inner_prod(vGauss, data.Normal);
//     const double vGaussSquaredNorm = std::pow(vGauss[0],2) + std::pow(vGauss[1],2) + std::pow(vGauss[2],2);

//     // Add outlet inflow prevention contribution
//     const double delta = 1.0e-2;
//     const double U_0 = rProcessInfo[CHARACTERISTIC_VELOCITY];
//     const double S_0 = 0.5*(1-tanh(vGaussProj/(U_0*delta)));

//     for (unsigned int i=0; i<TNumNodes; ++i)
//     {
//         unsigned int row = i*LocalSize;
//         for (unsigned int d=0; d<TDim; ++d)
//         {
//             rhs_gauss[row+d] += data.wGauss*data.N[i]*0.5*rho*vGaussSquaredNorm*S_0*data.Normal[d];
//         }
//     }
// }

// template<unsigned int TDim, unsigned int TNumNodes, class... TWallModel>
// void LowMachNavierStokesWallCondition<TDim,TNumNodes,TWallModel...>::CalculateGaussPointSlipTangentialCorrectionLHSContribution(
//     BoundedMatrix<double,LocalSize,LocalSize>& rLeftHandSideMatrix,
//     const ConditionDataStruct& rDataStruct)
// {
//     KRATOS_TRY

//     // Get element data
//     const auto& r_geom = this->GetGeometry();
//     const auto& r_N = rDataStruct.N;
//     const auto& r_cond_normal = rDataStruct.Normal;

//     // Set auxiliary condition normal to match array sizes
//     array_1d<double, TDim> aux_cond_normal;
//     if constexpr (TDim == 2) {
//         aux_cond_normal[0] = r_cond_normal[0];
//         aux_cond_normal[1] = r_cond_normal[1];
//     } else {
//         noalias(aux_cond_normal) = r_cond_normal;
//     }

//     // Allocate auxiliary arrays
//     array_1d<double, 3> i_node_unit_normal;
//     BoundedMatrix<double,TDim,TDim> tang_proj_mat;
//     array_1d<double, TDim> cauchy_traction_tang_proj;

//     for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node) {
//         // Set the nodal tangential projection matrix
//         noalias(i_node_unit_normal) = r_geom[i_node].FastGetSolutionStepValue(NORMAL);
//         i_node_unit_normal /= norm_2(i_node_unit_normal);
//         this->SetTangentialProjectionMatrix(i_node_unit_normal, tang_proj_mat);

//         // Get the spurious tangential component of the traction vector
//         // Note that in here we are projecting with the nodal tangential operator
//         noalias(cauchy_traction_tang_proj) = prod(tang_proj_mat, aux_cond_normal);

//         // Assemble the LHS contribution
//         // Note that only the pressure stress contribution is included in the linearisation
//         // The viscous stress contribution is dropped as it comes from the parent element
//         for (std::size_t j_node = 0; j_node < TNumNodes; ++j_node) {
//             for (std::size_t d = 0; d < TDim; ++d) {
//                 rLeftHandSideMatrix(i_node*BlockSize + d, j_node*BlockSize + TDim) += rDataStruct.wGauss * r_N[i_node] * cauchy_traction_tang_proj[d] *r_N[j_node];
//             }
//         }
//     }

//     KRATOS_CATCH("");
// }

// template<unsigned int TDim, unsigned int TNumNodes, class... TWallModel>
// void LowMachNavierStokesWallCondition<TDim,TNumNodes,TWallModel...>::CalculateGaussPointSlipTangentialCorrectionRHSContribution(
//     array_1d<double,LocalSize>& rRightHandSideVector,
//     const ConditionDataStruct& rDataStruct)
// {
//     KRATOS_TRY

//     // Get element data
//     const auto& r_geom = this->GetGeometry();
//     const auto& r_N = rDataStruct.N;
//     const auto& r_cond_normal = rDataStruct.Normal;
//     const auto& r_viscous_stress = rDataStruct.ViscousStress;

//     // Allocate auxiliary arrays
//     array_1d<double,3> i_node_unit_normal;
//     array_1d<double,VoigtSize> voigt_stress;
//     array_1d<double,TDim> cauchy_traction_vect;
//     BoundedMatrix<double,TDim,TDim> tang_proj_mat;
//     array_1d<double,TDim> cauchy_traction_tang_proj;

//     for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node) {
//         // Set the nodal tangential projection matrix
//         noalias(i_node_unit_normal) = r_geom[i_node].FastGetSolutionStepValue(NORMAL);
//         i_node_unit_normal /= norm_2(i_node_unit_normal);
//         this->SetTangentialProjectionMatrix(i_node_unit_normal, tang_proj_mat);

//         // Set the current Gauss point Cauchy traction vector with the condition normal
//         // Note that we add the corresponding nodal pressure to the constant viscous traction
//         cauchy_traction_vect = ZeroVector(TDim);
//         for (std::size_t j_node = 0; j_node < TNumNodes; ++j_node) {
//             if constexpr (VoigtSize == 3) {
//                 // Voigt stress
//                 voigt_stress[0] = r_viscous_stress[0] - r_geom[j_node].FastGetSolutionStepValue(PRESSURE);
//                 voigt_stress[1] = r_viscous_stress[1] - r_geom[j_node].FastGetSolutionStepValue(PRESSURE);
//                 voigt_stress[2] = r_viscous_stress[2]; // no pressure in shear component
//                 // Projection along the condition normal
//                 cauchy_traction_vect[0] += r_N[j_node]*(voigt_stress[0]*r_cond_normal[0] + voigt_stress[2]*r_cond_normal[1]);
//                 cauchy_traction_vect[1] += r_N[j_node]*(voigt_stress[2]*r_cond_normal[0] + voigt_stress[1]*r_cond_normal[1]);
//             } else if constexpr (VoigtSize == 6) {
//                 // Voigt stress
//                 voigt_stress[0] = r_viscous_stress[0] - r_geom[j_node].FastGetSolutionStepValue(PRESSURE);
//                 voigt_stress[1] = r_viscous_stress[1] - r_geom[j_node].FastGetSolutionStepValue(PRESSURE);
//                 voigt_stress[2] = r_viscous_stress[2] - r_geom[j_node].FastGetSolutionStepValue(PRESSURE);
//                 voigt_stress[3] = r_viscous_stress[3]; // no pressure in shear component
//                 voigt_stress[4] = r_viscous_stress[4]; // no pressure in shear component
//                 voigt_stress[5] = r_viscous_stress[5]; // no pressure in shear component
//                 // Projection along the condition normal
//                 cauchy_traction_vect[0] += r_N[j_node]*(voigt_stress[0]*r_cond_normal[0] + voigt_stress[3]*r_cond_normal[1] + voigt_stress[5]*r_cond_normal[2]);
//                 cauchy_traction_vect[1] += r_N[j_node]*(voigt_stress[3]*r_cond_normal[0] + voigt_stress[1]*r_cond_normal[1] + voigt_stress[4]*r_cond_normal[2]);
//                 cauchy_traction_vect[2] += r_N[j_node]*(voigt_stress[5]*r_cond_normal[0] + voigt_stress[4]*r_cond_normal[1] + voigt_stress[2]*r_cond_normal[2]);
//             }
//         }

//         // Get the spurious tangential component of the traction vector
//         // Note that in here we are projecting with the nodal tangential operator
//         noalias(cauchy_traction_tang_proj) = prod(tang_proj_mat, cauchy_traction_vect);

//         // Assemble the RHS contribution
//         for (std::size_t d = 0; d < TDim; ++d) {
//             rRightHandSideVector[i_node*BlockSize + d] += rDataStruct.wGauss * r_N[i_node] * cauchy_traction_tang_proj[d];
//         }
//     }

//     KRATOS_CATCH("")
// }

// template<unsigned int TDim, unsigned int TNumNodes, class... TWallModel>
// void LowMachNavierStokesWallCondition<TDim,TNumNodes,TWallModel...>::ProjectViscousStress(
//     const Vector& rViscousStress,
//     const array_1d<double,3> rNormal,
//     array_1d<double,3>& rProjectedViscousStress)
// {
//     if constexpr (TDim == 2) {
//         rProjectedViscousStress[0] = rViscousStress[0] * rNormal[0] + rViscousStress[2] * rNormal[1];
//         rProjectedViscousStress[1] = rViscousStress[2] * rNormal[0] + rViscousStress[1] * rNormal[1];
//         rProjectedViscousStress[2] = 0.0;
//     } else {
//         rProjectedViscousStress[0] = rViscousStress[0] * rNormal[0] + rViscousStress[3] * rNormal[1] + rViscousStress[5] * rNormal[2];
//         rProjectedViscousStress[1] = rViscousStress[3] * rNormal[0] + rViscousStress[1] * rNormal[1] + rViscousStress[4] * rNormal[2];
//         rProjectedViscousStress[2] = rViscousStress[5] * rNormal[0] + rViscousStress[4] * rNormal[1] + rViscousStress[2] * rNormal[2];
//     }
// }

template class LowMachNavierStokesWallCondition<2,2>;
// template class LowMachNavierStokesWallCondition<3,3>;
// template class LowMachNavierStokesWallCondition<2,2,LinearLogWallLaw<2,2>>;
// template class LowMachNavierStokesWallCondition<3,3,LinearLogWallLaw<3,3>>;
// template class LowMachNavierStokesWallCondition<2,2,NavierSlipWallLaw<2,2>>;
// template class LowMachNavierStokesWallCondition<3,3,NavierSlipWallLaw<3,3>>;

} // namespace Kratos
