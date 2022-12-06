//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#include "navier_stokes_wall_condition.h"
#include "includes/checks.h"

namespace Kratos
{

///@name Specialized implementation of VMS for functions that depend on TDim
///@{


/**
 * @see NavierStokesWallCondition::EquationIdVector
 */
template <>
void NavierStokesWallCondition<2,2>::EquationIdVector(EquationIdVectorType& rResult,
                                                      const ProcessInfo& rCurrentProcessInfo) const
{
    const unsigned int NumNodes = 2;
    const unsigned int LocalSize = 6;
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
    }
}


/**
 * @see NavierStokesWallCondition::EquationIdVector
 */
template <>
void NavierStokesWallCondition<3,3>::EquationIdVector(EquationIdVectorType& rResult,
                                                      const ProcessInfo& rCurrentProcessInfo) const
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 12;
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z).EquationId();
        rResult[LocalIndex++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
    }
}



template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                      VectorType& rRightHandSideVector,
                                      const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);
    if (rLeftHandSideMatrix.size1() != MatrixSize)
        rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); //false says not to preserve existing storage!!
    if (rRightHandSideVector.size() != MatrixSize)
        rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

    // Check that parents have been computed
    // These are required to retrieve the material properties and the viscous stress
    auto& parentElement = this->GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_ERROR_IF(parentElement.size() > 1) << "A condition was assigned more than one parent element." << std::endl;
    KRATOS_ERROR_IF(parentElement.size() == 0) << "A condition was NOT assigned a parent element. Please execute the check_and_prepare_model_process_fluid process." << std::endl;

    // Struct to pass around the data
    ConditionDataStruct data;
    // Allocate memory needed
    array_1d<double,MatrixSize> rhs_gauss;
    BoundedMatrix<double,MatrixSize, MatrixSize> lhs_gauss;

    // LHS and RHS contributions initialization
    noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);
    noalias(rRightHandSideVector) = ZeroVector(MatrixSize);

    // Compute condition unit normal vector
    this->CalculateNormal(data.Normal); //this already contains the area
    const double A = norm_2(data.Normal);
    data.Normal /= A;

    // Gauss point information
    GeometryType& rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);
    const unsigned int NumGauss = IntegrationPoints.size();
    Vector GaussPtsJDet = ZeroVector(NumGauss);
    rGeom.DeterminantOfJacobian(GaussPtsJDet, GeometryData::IntegrationMethod::GI_GAUSS_2);
    const MatrixType Ncontainer = rGeom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);

    if ( this->Is(SLIP) ){
        // finding parent element to retrieve viscous stresses which are later stored in "data"
        Element& parent = parentElement[0];
        data.ViscousStress = ZeroVector( 3*(TDim-1) );
        parent.Calculate(FLUID_STRESS, data.ViscousStress, rCurrentProcessInfo);
    }

    // Loop on gauss points
    for(unsigned int igauss = 0; igauss<NumGauss; igauss++)
    {
        data.N = row(Ncontainer, igauss);
        const double J = GaussPtsJDet[igauss];
        data.wGauss = J * IntegrationPoints[igauss].Weight();
        ComputeGaussPointRHSContribution(rhs_gauss, data, rCurrentProcessInfo);
        ComputeGaussPointLHSContribution(lhs_gauss, data, rCurrentProcessInfo);
        noalias(rLeftHandSideMatrix) += lhs_gauss;
        noalias(rRightHandSideVector) += rhs_gauss;
    }

    KRATOS_CATCH("")
}



template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
                                   const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

    if (rLeftHandSideMatrix.size1() != MatrixSize)
        rLeftHandSideMatrix.resize(MatrixSize, MatrixSize, false); //false says not to preserve existing storage!!

    // LHS contributions initialization
    noalias(rLeftHandSideMatrix) = ZeroMatrix(MatrixSize,MatrixSize);

    KRATOS_CATCH("")
}



template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int MatrixSize = TNumNodes*(TDim+1);

    if (rRightHandSideVector.size() != MatrixSize)
        rRightHandSideVector.resize(MatrixSize, false); //false says not to preserve existing storage!!

    // Struct to pass around the data
    ConditionDataStruct data;
    // Allocate memory needed
    array_1d<double,MatrixSize> rhs_gauss;
    // Loop on gauss points
    noalias(rRightHandSideVector) = ZeroVector(MatrixSize);

    // Compute condition normal
    this->CalculateNormal(data.Normal); //this already contains the area
    const double A = norm_2(data.Normal);
    data.Normal /= A;

    // Gauss point information
    GeometryType& rGeom = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);
    const unsigned int NumGauss = IntegrationPoints.size();
    Vector GaussPtsJDet = ZeroVector(NumGauss);
    rGeom.DeterminantOfJacobian(GaussPtsJDet, GeometryData::IntegrationMethod::GI_GAUSS_2);
    const MatrixType Ncontainer = rGeom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);

    if ( this->Is(SLIP) ){
        // finding parent element to retrieve viscous stresses which are later stored in "data"
        GlobalPointersVector<Element> parentElement = this->GetValue( NEIGHBOUR_ELEMENTS );
        KRATOS_ERROR_IF( parentElement.size() > 1 ) << "A condition was assigned more than one parent element." << std::endl;
        KRATOS_ERROR_IF( parentElement.size() == 0 ) << "A condition was NOT assigned a parent element. "
        << "This leads to errors for the slip condition [BEHR2004] "
        << "Please execute the check_and_prepare_model_process_fluid process." << std::endl;

        Element& parent = parentElement[0];
        data.ViscousStress = ZeroVector( 3*(TDim-1) );
        parent.Calculate(FLUID_STRESS, data.ViscousStress, rCurrentProcessInfo);
    }

    for(unsigned int igauss = 0; igauss<NumGauss; igauss++)
    {
        data.N = row(Ncontainer, igauss);
        const double J = GaussPtsJDet[igauss];
        data.wGauss = J * IntegrationPoints[igauss].Weight();
        ComputeGaussPointRHSContribution(rhs_gauss, data, rCurrentProcessInfo);
        noalias(rRightHandSideVector) += rhs_gauss;
    }

    KRATOS_CATCH("")
}


/// Condition check
/**
 * @param rCurrentProcessInfo reference to the ProcessInfo
 */
template<unsigned int TDim, unsigned int TNumNodes>
int NavierStokesWallCondition<TDim,TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;
    int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0
    if (Check != 0) {
        return Check;
    }
    else {
        // Checks on nodes
        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {
            if(this->GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_ERROR << "missing VELOCITY variable on solution step data for node " << this->GetGeometry()[i].Id();
            if(this->GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
                KRATOS_ERROR << "missing PRESSURE variable on solution step data for node " << this->GetGeometry()[i].Id();
            if(this->GetGeometry()[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
                KRATOS_ERROR << "missing MESH_VELOCITY variable on solution step data for node " << this->GetGeometry()[i].Id();
            if(this->GetGeometry()[i].SolutionStepsDataHas(ACCELERATION) == false)
                KRATOS_ERROR << "missing ACCELERATION variable on solution step data for node " << this->GetGeometry()[i].Id();
            if(this->GetGeometry()[i].SolutionStepsDataHas(EXTERNAL_PRESSURE) == false)
                KRATOS_ERROR << "missing EXTERNAL_PRESSURE variable on solution step data for node " << this->GetGeometry()[i].Id();
            if(this->GetGeometry()[i].HasDofFor(VELOCITY_X) == false ||
               this->GetGeometry()[i].HasDofFor(VELOCITY_Y) == false ||
               this->GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
                KRATOS_ERROR << "missing VELOCITY component degree of freedom on node " << this->GetGeometry()[i].Id();
            if(this->GetGeometry()[i].HasDofFor(PRESSURE) == false)
                KRATOS_ERROR << "missing PRESSURE component degree of freedom on node " << this->GetGeometry()[i].Id();
        }

        return Check;
    }

    KRATOS_CATCH("");
}


/**
 * @see NavierStokesWallCondition::GetDofList
 */
template <>
void NavierStokesWallCondition<2,2>::GetDofList(DofsVectorType& rElementalDofList,
                                                const ProcessInfo& rCurrentProcessInfo) const
{
    const SizeType NumNodes = 2;
    const SizeType LocalSize = 6;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
    }
}



template <>
void NavierStokesWallCondition<3,3>::GetDofList(DofsVectorType& rElementalDofList,
                                                const ProcessInfo& rCurrentProcessInfo) const
{
    const SizeType NumNodes = 3;
    const SizeType LocalSize = 12;

    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for (unsigned int iNode = 0; iNode < NumNodes; ++iNode)
    {
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Z);
        rElementalDofList[LocalIndex++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
    }
}

template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim, TNumNodes>::Calculate(
    const Variable< array_1d<double,3> >& rVariable,
    array_1d<double,3>& rOutput,
    const ProcessInfo& rCurrentProcessInfo)
{
    rOutput = ZeroVector(3);

    if (rVariable == DRAG_FORCE) {
        const auto& r_geom = GetGeometry();
        const auto& r_integration_points = r_geom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);
        unsigned int n_gauss = r_integration_points.size();
        Vector det_J_vect = ZeroVector(n_gauss);
        r_geom.DeterminantOfJacobian(det_J_vect, GeometryData::IntegrationMethod::GI_GAUSS_2);
        const auto N_container = r_geom.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);

        // Calculate normal
        array_1d<double,3> normal;
        CalculateNormal(normal);
        normal /= norm_2(normal);

        // Finding parent element to retrieve viscous stresses
        // Note that we assume in here that the shear stress is constant inside the element
        auto& r_neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);
        KRATOS_ERROR_IF(r_neighbours.size() > 1) << "A condition was assigned more than one parent element." << std::endl;
        KRATOS_ERROR_IF(r_neighbours.size() == 0) << "A condition was NOT assigned a parent element. "
        << "This leads to errors for the slip condition [BEHR2004] "
        << "Please execute the check_and_prepare_model_process_fluid process." << std::endl;

        auto& r_parent = r_neighbours[0];
        Vector shear_stress;
        r_parent.Calculate(FLUID_STRESS, shear_stress, rCurrentProcessInfo);
        array_1d<double,3> shear_stress_n;
        ProjectViscousStress(shear_stress, normal, shear_stress_n);

        // Loop the Gauss pts
        for (unsigned int i_gauss = 0; i_gauss < n_gauss; ++i_gauss) {
            const double w = det_J_vect[i_gauss] * r_integration_points[i_gauss].Weight();
            double p = 0.0;
            const auto& r_N = row(N_container, i_gauss);
            for (unsigned int i_node = 0; i_node < r_geom.PointsNumber(); ++i_node) {
                p += r_N[i_node] * r_geom[i_node].FastGetSolutionStepValue(PRESSURE);
            }
            rOutput += w * (p * normal - shear_stress_n);
        }
    } else {
        Condition::Calculate(rVariable, rOutput, rCurrentProcessInfo);
    }
}


template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::ComputeGaussPointLHSContribution(
    BoundedMatrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)>& lhs_gauss,
    const ConditionDataStruct& data,
    const ProcessInfo& rProcessInfo)
{
    const unsigned int LocalSize = TDim+1;
    lhs_gauss = ZeroMatrix(TNumNodes*LocalSize, TNumNodes*LocalSize);

    //TODO: Add a proper switch to activate this
    // if (this->Is(SLIP)){
    //     ComputeGaussPointNavierSlipLHSContribution( lhs_gauss, data );
    // }

    // Contribution to avoid spurious tangential components in the pure-slip residual
    if (rProcessInfo.Has(SLIP_TANGENTIAL_CORRECTION_SWITCH)) {
        if (this->Is(SLIP) && rProcessInfo.GetValue(SLIP_TANGENTIAL_CORRECTION_SWITCH)) {
            CalculateGaussPointSlipTangentialCorrectionLHSContribution(lhs_gauss, data);
        }
    }
}



template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::ComputeGaussPointRHSContribution(
    array_1d<double,TNumNodes*(TDim+1)>& rhs_gauss,
    const ConditionDataStruct& data,
    const ProcessInfo& rProcessInfo)
{
    // Initialize the local RHS
    const unsigned int LocalSize = TDim+1;
    noalias(rhs_gauss) = ZeroVector(TNumNodes*LocalSize);

    // Gauss pt. Neumann BC contribution
    this->ComputeRHSNeumannContribution(rhs_gauss, data);

    // Gauss pt. outlet inflow prevention contribution
    if (rProcessInfo.Has(OUTLET_INFLOW_CONTRIBUTION_SWITCH)) {
        if (this->Is(OUTLET) && rProcessInfo[OUTLET_INFLOW_CONTRIBUTION_SWITCH]){
            this->ComputeRHSOutletInflowContribution(rhs_gauss, data, rProcessInfo);
        }
    }

    //TODO: Add a proper switch to activate this
    // if (this->Is(SLIP)){
    //     ComputeGaussPointNavierSlipRHSContribution( rhs_gauss, data );
    // }

    // Contribution to avoid spurious tangential components in the pure-slip residual
    if (rProcessInfo.Has(SLIP_TANGENTIAL_CORRECTION_SWITCH)) {
        if (this->Is(SLIP) && rProcessInfo[SLIP_TANGENTIAL_CORRECTION_SWITCH]) {
            CalculateGaussPointSlipTangentialCorrectionRHSContribution(rhs_gauss, data);
        }
    }

}



template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::ComputeRHSNeumannContribution(array_1d<double,TNumNodes*(TDim+1)>& rhs_gauss,
                                                                              const ConditionDataStruct& data)
{
    const unsigned int LocalSize = TDim+1;
    const GeometryType& rGeom = this->GetGeometry();

    // Add Neumann BC contribution
    for (unsigned int i=0; i<TNumNodes; ++i)
    {
        const double pext = rGeom[i].FastGetSolutionStepValue(EXTERNAL_PRESSURE);

        for (unsigned int j=0; j<TNumNodes; ++j)
        {
            unsigned int row = j*LocalSize;
            for (unsigned int d=0; d<TDim; ++d)
            {
                rhs_gauss[row+d] -= data.wGauss*data.N[j]*data.N[i]*pext*data.Normal[d];
            }
        }
    }
}


template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::ComputeRHSOutletInflowContribution(
    array_1d<double,TNumNodes*(TDim+1)>& rhs_gauss,
    const ConditionDataStruct& data,
    const ProcessInfo& rProcessInfo)
{
    constexpr SizeType LocalSize = TDim+1;
    const GeometryType& rGeom = this->GetGeometry();

    // Get DENSITY from parent element properties
    auto & r_neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);
    const double rho = r_neighbours[0].GetProperties().GetValue(DENSITY);

    // Compute Gauss pt. density, velocity norm and velocity projection
    array_1d<double, 3> vGauss = ZeroVector(3);
    for (unsigned int i=0; i<TNumNodes; ++i)
    {
        const array_1d<double, 3>& rVelNode = rGeom[i].FastGetSolutionStepValue(VELOCITY);
        vGauss += data.N[i]*rVelNode;
    }

    const double vGaussProj = inner_prod(vGauss, data.Normal);
    const double vGaussSquaredNorm = std::pow(vGauss[0],2) + std::pow(vGauss[1],2) + std::pow(vGauss[2],2);

    // Add outlet inflow prevention contribution
    const double delta = 1.0e-2;
    const double U_0 = rProcessInfo[CHARACTERISTIC_VELOCITY];
    const double S_0 = 0.5*(1-tanh(vGaussProj/(U_0*delta)));

    for (unsigned int i=0; i<TNumNodes; ++i)
    {
        unsigned int row = i*LocalSize;
        for (unsigned int d=0; d<TDim; ++d)
        {
            rhs_gauss[row+d] += data.wGauss*data.N[i]*0.5*rho*vGaussSquaredNorm*S_0*data.Normal[d];
        }
    }
}


/// Computes the 2D condition normal
/**
* @param An reference to condition normal vector
*/
template <>
void NavierStokesWallCondition<2,2>::CalculateNormal(array_1d<double,3>& An)
{
    Geometry<Node<3> >& pGeometry = this->GetGeometry();

    An[0] =   pGeometry[1].Y() - pGeometry[0].Y();
    An[1] = - (pGeometry[1].X() - pGeometry[0].X());
    An[2] =    0.0;

}


/// Computes the 3D condition normal
/**
* @param An reference to condition normal vector
*/
template <>
void NavierStokesWallCondition<3,3>::CalculateNormal(array_1d<double,3>& An )
{
    Geometry<Node<3> >& pGeometry = this->GetGeometry();

    array_1d<double,3> v1,v2;
    v1[0] = pGeometry[1].X() - pGeometry[0].X();
    v1[1] = pGeometry[1].Y() - pGeometry[0].Y();
    v1[2] = pGeometry[1].Z() - pGeometry[0].Z();

    v2[0] = pGeometry[2].X() - pGeometry[0].X();
    v2[1] = pGeometry[2].Y() - pGeometry[0].Y();
    v2[2] = pGeometry[2].Z() - pGeometry[0].Z();

    MathUtils<double>::CrossProduct(An,v1,v2);
    An *= 0.5;
}

template<>
void NavierStokesWallCondition<2,2>::ProjectViscousStress(
    const Vector& rViscousStress,
    const array_1d<double,3> rNormal,
    array_1d<double,3>& rProjectedViscousStress)
{
    rProjectedViscousStress[0] = rViscousStress[0] * rNormal[0] + rViscousStress[2] * rNormal[1];
    rProjectedViscousStress[1] = rViscousStress[2] * rNormal[0] + rViscousStress[1] * rNormal[1];
    rProjectedViscousStress[2] = 0.0;
}

template<>
void NavierStokesWallCondition<3,3>::ProjectViscousStress(
    const Vector& rViscousStress,
    const array_1d<double,3> rNormal,
    array_1d<double,3>& rProjectedViscousStress)
{
    rProjectedViscousStress[0] = rViscousStress[0] * rNormal[0] + rViscousStress[3] * rNormal[1] + rViscousStress[5] * rNormal[2];
    rProjectedViscousStress[1] = rViscousStress[3] * rNormal[0] + rViscousStress[1] * rNormal[1] + rViscousStress[4] * rNormal[2];
    rProjectedViscousStress[2] = rViscousStress[5] * rNormal[0] + rViscousStress[4] * rNormal[1] + rViscousStress[2] * rNormal[2];
}

template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::CalculateGaussPointSlipTangentialCorrectionLHSContribution(
    BoundedMatrix<double,LocalSize,LocalSize>& rLeftHandSideMatrix,
    const ConditionDataStruct& rDataStruct)
{
    KRATOS_TRY

    // Get element data
    const auto& r_geom = this->GetGeometry();
    const auto& r_N = rDataStruct.N;
    const auto& r_cond_normal = rDataStruct.Normal;

    // Set auxiliary condition normal to match array sizes
    array_1d<double, TDim> aux_cond_normal;
    if constexpr (TDim == 2) {
        aux_cond_normal[0] = r_cond_normal[0];
        aux_cond_normal[1] = r_cond_normal[1];
    } else {
        noalias(aux_cond_normal) = r_cond_normal;
    }

    // Allocate auxiliary arrays
    array_1d<double, 3> i_node_unit_normal;
    BoundedMatrix<double,TDim,TDim> tang_proj_mat;
    array_1d<double, TDim> cauchy_traction_tang_proj;

    for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node) {
        // Set the nodal tangential projection matrix
        noalias(i_node_unit_normal) = r_geom[i_node].FastGetSolutionStepValue(NORMAL);
        i_node_unit_normal /= norm_2(i_node_unit_normal);
        this->SetTangentialProjectionMatrix(i_node_unit_normal, tang_proj_mat);

        // Get the spurious tangential component of the traction vector
        // Note that in here we are projecting with the nodal tangential operator
        noalias(cauchy_traction_tang_proj) = prod(tang_proj_mat, aux_cond_normal);

        // Assemble the LHS contribution
        // Note that only the pressure stress contribution is included in the linearisation
        // The viscous stress contribution is dropped as it comes from the parent element
        for (std::size_t j_node = 0; j_node < TNumNodes; ++j_node) {
            for (std::size_t d = 0; d < TDim; ++d) {
                rLeftHandSideMatrix(i_node*BlockSize + d, j_node*BlockSize + TDim) += rDataStruct.wGauss * r_N[i_node] * cauchy_traction_tang_proj[d] *r_N[j_node];
            }
        }
    }

    KRATOS_CATCH("");
}




template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::CalculateGaussPointSlipTangentialCorrectionRHSContribution(
    array_1d<double,LocalSize>& rRightHandSideVector,
    const ConditionDataStruct& rDataStruct)
{
    KRATOS_TRY

    // Get element data
    const auto& r_geom = this->GetGeometry();
    const auto& r_N = rDataStruct.N;
    const auto& r_cond_normal = rDataStruct.Normal;
    const auto& r_viscous_stress = rDataStruct.ViscousStress;

    // Allocate auxiliary arrays
    array_1d<double,3> i_node_unit_normal;
    array_1d<double,VoigtSize> voigt_stress;
    array_1d<double,TDim> cauchy_traction_vect;
    BoundedMatrix<double,TDim,TDim> tang_proj_mat;
    array_1d<double,TDim> cauchy_traction_tang_proj;

    for (std::size_t i_node = 0; i_node < TNumNodes; ++i_node) {
        // Set the nodal tangential projection matrix
        noalias(i_node_unit_normal) = r_geom[i_node].FastGetSolutionStepValue(NORMAL);
        i_node_unit_normal /= norm_2(i_node_unit_normal);
        this->SetTangentialProjectionMatrix(i_node_unit_normal, tang_proj_mat);

        // Set the current Gauss point Cauchy traction vector with the condition normal
        // Note that we add the corresponding nodal pressure to the constant viscous traction
        cauchy_traction_vect = ZeroVector(TDim);
        for (std::size_t j_node = 0; j_node < TNumNodes; ++j_node) {
            if constexpr (VoigtSize == 3) {
                // Voigt stress
                voigt_stress[0] = r_viscous_stress[0] - r_geom[j_node].FastGetSolutionStepValue(PRESSURE);
                voigt_stress[1] = r_viscous_stress[1] - r_geom[j_node].FastGetSolutionStepValue(PRESSURE);
                voigt_stress[2] = r_viscous_stress[2]; // no pressure in shear component
                // Projection along the condition normal
                cauchy_traction_vect[0] += r_N[j_node]*(voigt_stress[0]*r_cond_normal[0] + voigt_stress[2]*r_cond_normal[1]);
                cauchy_traction_vect[1] += r_N[j_node]*(voigt_stress[2]*r_cond_normal[0] + voigt_stress[1]*r_cond_normal[1]);
            } else if constexpr (VoigtSize == 6) {
                // Voigt stress
                voigt_stress[0] = r_viscous_stress[0] - r_geom[j_node].FastGetSolutionStepValue(PRESSURE);
                voigt_stress[1] = r_viscous_stress[1] - r_geom[j_node].FastGetSolutionStepValue(PRESSURE);
                voigt_stress[2] = r_viscous_stress[2] - r_geom[j_node].FastGetSolutionStepValue(PRESSURE);
                voigt_stress[3] = r_viscous_stress[3]; // no pressure in shear component
                voigt_stress[4] = r_viscous_stress[4]; // no pressure in shear component
                voigt_stress[5] = r_viscous_stress[5]; // no pressure in shear component
                // Projection along the condition normal
                cauchy_traction_vect[0] += r_N[j_node]*(voigt_stress[0]*r_cond_normal[0] + voigt_stress[3]*r_cond_normal[1] + voigt_stress[5]*r_cond_normal[2]);
                cauchy_traction_vect[1] += r_N[j_node]*(voigt_stress[3]*r_cond_normal[0] + voigt_stress[1]*r_cond_normal[1] + voigt_stress[4]*r_cond_normal[2]);
                cauchy_traction_vect[2] += r_N[j_node]*(voigt_stress[5]*r_cond_normal[0] + voigt_stress[4]*r_cond_normal[1] + voigt_stress[2]*r_cond_normal[2]);
            }
        }

        // Get the spurious tangential component of the traction vector
        // Note that in here we are projecting with the nodal tangential operator
        noalias(cauchy_traction_tang_proj) = prod(tang_proj_mat, cauchy_traction_vect);

        // Assemble the RHS contribution
        for (std::size_t d = 0; d < TDim; ++d) {
            rRightHandSideVector[i_node*BlockSize + d] += rDataStruct.wGauss * r_N[i_node] * cauchy_traction_tang_proj[d];
        }
    }

    KRATOS_CATCH("")
}



template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::ComputeGaussPointNavierSlipRHSContribution(
    array_1d<double,TNumNodes*(TDim+1)>& rRightHandSideVector,
    const ConditionDataStruct& rDataStruct)
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();
    GlobalPointersVector<Element> parentElement = this->GetValue(NEIGHBOUR_ELEMENTS);
    const double viscosity = parentElement[0].GetProperties().GetValue(DYNAMIC_VISCOSITY);

    const array_1d<double, TNumNodes> N = rDataStruct.N;
    const double wGauss = rDataStruct.wGauss;

    for (unsigned int nnode = 0; nnode < TNumNodes; nnode++){

        // finding the nodal projection matrix nodal_projection_matrix = ( [I] - (na)(na) )
        BoundedMatrix<double, TNumNodes, TNumNodes> nodal_projection_matrix;
        array_1d<double,3> nodal_normal = rGeom[nnode].FastGetSolutionStepValue(NORMAL);
        double sum_of_squares = 0.0;
        for (unsigned int j = 0; j < 3; j++){
            sum_of_squares += nodal_normal[j] * nodal_normal[j];
        }
        nodal_normal /= sqrt(sum_of_squares);
        FluidElementUtilities<3>::SetTangentialProjectionMatrix( nodal_normal, nodal_projection_matrix );

        // finding the coefficent to relate velocity to drag
        const double navier_slip_length = rGeom[nnode].GetValue(SLIP_LENGTH);
        KRATOS_ERROR_IF_NOT( navier_slip_length > 0.0 ) << "Negative or zero slip length was defined" << std::endl;
        const double nodal_beta = viscosity / navier_slip_length;


        Vector interpolated_velocity = ZeroVector(TNumNodes);
        for( unsigned int comp = 0; comp < TNumNodes; comp++){
            for (unsigned int i = 0; i < TNumNodes; i++){
                // necessary because VELOCITY with 3 entries even in 2D case
                interpolated_velocity[i] -= N[comp] * rGeom[comp].FastGetSolutionStepValue(VELOCITY)[i];
            }
        }
        // application of the nodal projection matrix
        const array_1d<double,TNumNodes> nodal_entry_rhs = prod( nodal_projection_matrix, (wGauss * N[nnode] * nodal_beta * interpolated_velocity) );
        for (unsigned int entry = 0; entry < TNumNodes; entry++){
            rRightHandSideVector( nnode*(TNumNodes+1) + entry ) += nodal_entry_rhs[entry];
        }
    }

    KRATOS_CATCH("")
}


template<unsigned int TDim, unsigned int TNumNodes>
void NavierStokesWallCondition<TDim,TNumNodes>::ComputeGaussPointNavierSlipLHSContribution(
    BoundedMatrix<double,TNumNodes*(TDim+1),TNumNodes*(TDim+1)>& rLeftHandSideMatrix,
    const ConditionDataStruct& rDataStruct)
{
    KRATOS_TRY

    const GeometryType& rGeom = this->GetGeometry();
    GlobalPointersVector<Element> parentElement = this->GetValue(NEIGHBOUR_ELEMENTS);
    const double viscosity = parentElement[0].GetProperties().GetValue(DYNAMIC_VISCOSITY);

    array_1d<double, TNumNodes> N = rDataStruct.N;
    const double wGauss = rDataStruct.wGauss;

    for(unsigned int inode = 0; inode < TNumNodes; inode++){

        // finding the nodal projection matrix nodal_projection_matrix = ( [I] - (na)(na) )
        BoundedMatrix<double, TNumNodes, TNumNodes> nodal_projection_matrix;
        array_1d<double,3> nodal_normal = rGeom[inode].FastGetSolutionStepValue(NORMAL);
        double sum_of_squares = 0.0;
        for (unsigned int j = 0; j < 3; j++){
            sum_of_squares += nodal_normal[j] * nodal_normal[j];
        }
        nodal_normal /= sqrt(sum_of_squares);
        FluidElementUtilities<3>::SetTangentialProjectionMatrix( nodal_normal, nodal_projection_matrix );

        // finding the coefficent to relate velocity to drag
        const double navier_slip_length = rGeom[inode].GetValue(SLIP_LENGTH);
        KRATOS_ERROR_IF_NOT( navier_slip_length > 0.0 ) << "Negative or zero slip length was defined" << std::endl;
        const double nodal_beta = viscosity / navier_slip_length;

        for(unsigned int jnode = 0; jnode < TNumNodes; jnode++){

            const BoundedMatrix<double, TNumNodes, TNumNodes> nodal_lhs_contribution = wGauss * nodal_beta * N[inode] * N[jnode] * nodal_projection_matrix;

            for( unsigned int i = 0; i < TNumNodes; i++){
                for( unsigned int j = 0; j < TNumNodes; j++){

                    const unsigned int istart = inode * (TNumNodes+1);
                    const unsigned int jstart = jnode * (TNumNodes+1);
                    rLeftHandSideMatrix(istart + i, jstart + j) += nodal_lhs_contribution(i,j);
                }
            }
        }
    }

    KRATOS_CATCH("")
}


template class NavierStokesWallCondition<2,2>;
template class NavierStokesWallCondition<3,3>;

} // namespace Kratos
