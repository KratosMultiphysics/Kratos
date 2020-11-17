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

#include "compressible_navier_stokes_wall_condition.h"
#include "includes/checks.h"

namespace Kratos
{

template <>
void CompressibleNavierStokesWallCondition<2,2>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const unsigned int n_nodes = 2;
    const unsigned int local_size = 8;

    if (rResult.size() != local_size) {
        rResult.resize(local_size, false);
    }

    unsigned int local_index = 0;
    const auto& r_geom = this->GetGeometry();
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_node = r_geom[i_node];
        rResult[local_index++] = r_node.GetDof(DENSITY).EquationId();
        rResult[local_index++] = r_node.GetDof(MOMENTUM_X).EquationId();
        rResult[local_index++] = r_node.GetDof(MOMENTUM_Y).EquationId();
        rResult[local_index++] = r_node.GetDof(TOTAL_ENERGY).EquationId();
    }
}

template <>
void CompressibleNavierStokesWallCondition<3,3>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const unsigned int n_nodes = 3;
    const unsigned int local_size = 15;

    if (rResult.size() != local_size) {
        rResult.resize(local_size, false);
    }

    unsigned int local_index = 0;
    const auto& r_geom = this->GetGeometry();
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_node = r_geom[i_node];
        rResult[local_index++] = r_node.GetDof(DENSITY).EquationId();
        rResult[local_index++] = r_node.GetDof(MOMENTUM_X).EquationId();
        rResult[local_index++] = r_node.GetDof(MOMENTUM_Y).EquationId();
        rResult[local_index++] = r_node.GetDof(MOMENTUM_Z).EquationId();
        rResult[local_index++] = r_node.GetDof(TOTAL_ENERGY).EquationId();
    }
}

template <>
void CompressibleNavierStokesWallCondition<2,2>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const SizeType n_nodes = 2;
    const SizeType local_size = 8;

    if (rElementalDofList.size() != local_size) {
        rElementalDofList.resize(local_size);
    }

    unsigned int local_index = 0;
    const auto& r_geom = this->GetGeometry();
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_node = r_geom[i_node];
        rElementalDofList[local_index++] = r_node.pGetDof(DENSITY);
        rElementalDofList[local_index++] = r_node.pGetDof(MOMENTUM_X);
        rElementalDofList[local_index++] = r_node.pGetDof(MOMENTUM_Y);
        rElementalDofList[local_index++] = r_node.pGetDof(TOTAL_ENERGY);
    }
}

template <>
void CompressibleNavierStokesWallCondition<3,3>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const SizeType n_nodes = 3;
    const SizeType local_size = 15;

    if (rElementalDofList.size() != local_size) {
        rElementalDofList.resize(local_size);
    }

    unsigned int local_index = 0;
    const auto& r_geom = this->GetGeometry();
    for (unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_node = r_geom[i_node];
        rElementalDofList[local_index++] = r_node.pGetDof(DENSITY);
        rElementalDofList[local_index++] = r_node.pGetDof(MOMENTUM_X);
        rElementalDofList[local_index++] = r_node.pGetDof(MOMENTUM_Y);
        rElementalDofList[local_index++] = r_node.pGetDof(MOMENTUM_Z);
        rElementalDofList[local_index++] = r_node.pGetDof(TOTAL_ENERGY);
    }
}

template<unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void CompressibleNavierStokesWallCondition<TDim,TNumNodes,TBlockSize>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "Implicit \'CalculateLocalSystem\' has not been implemented yet." << std::endl;

    KRATOS_CATCH("")
}

template<unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void CompressibleNavierStokesWallCondition<TDim,TNumNodes,TBlockSize>::CalculateLeftHandSide(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "Implicit \'CalculateLeftHandSide\' has not been implemented yet." << std::endl;

    KRATOS_CATCH("")
}

template<unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void CompressibleNavierStokesWallCondition<TDim,TNumNodes,TBlockSize>::CalculateRightHandSide(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR << "Implicit \'CalculateRightHandSide\' has not been implemented yet." << std::endl;

    KRATOS_CATCH("")
}



template<unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
int CompressibleNavierStokesWallCondition<TDim,TNumNodes,TBlockSize>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    // Base condition check (id and negative domain size)
    int check = Condition::Check(rCurrentProcessInfo);

    if (check != 0) {
        return check;
    } else {
        // Check that the element's nodes contain all required SolutionStepData and DOFs
        const auto& r_geom = this->GetGeometry();
        const unsigned int n_nodes  = r_geom.PointsNumber();
        for(unsigned int i_node = 0; i_node < n_nodes; ++i_node) {
            const auto& r_node = r_geom[i_node];
            // TODO: Add missing variables
            if(!r_node.SolutionStepsDataHas(DENSITY)) {
                KRATOS_ERROR << "missing DENSITY variable on solution step data for node " << r_node.Id();
            }
            if(!r_node.SolutionStepsDataHas(MOMENTUM)) {
                KRATOS_ERROR << "missing MOMENTUM variable on solution step data for node " << r_node.Id();
            }
            if(!r_node.SolutionStepsDataHas(TOTAL_ENERGY)) {
                KRATOS_ERROR << "missing TOTAL_ENERGY variable on solution step data for node " << r_node.Id();
            }
            if(!r_node.HasDofFor(DENSITY)) {
                KRATOS_ERROR << "missing DENSITY component degree of freedom on node " << r_node.Id();
            }
            if(!r_node.HasDofFor(MOMENTUM_X) || !r_node.HasDofFor(MOMENTUM_Y) || (TDim == 3 ? !r_node.HasDofFor(MOMENTUM_Z) : false)) {
                KRATOS_ERROR << "missing MOMENTUM_X component degree of freedom on node " << r_node.Id();
            }
            if(!r_node.HasDofFor(TOTAL_ENERGY)) {
                KRATOS_ERROR << "missing TOTAL_ENERGY component degree of freedom on node " << r_node.Id();
            }
        }

        return check;
    }

    KRATOS_CATCH("");
}

template<unsigned int TDim, unsigned int TNumNodes, unsigned int TBlockSize>
void CompressibleNavierStokesWallCondition<TDim,TNumNodes,TBlockSize>::FillConditionData(
    ConditionDataStruct &rData,
    const ProcessInfo &rCurrentProcessInfo)
{
    // Getting data for the given geometry
    const auto& r_geometry = GetGeometry();
    GeometryUtils::CalculateGeometryData(r_geometry, rData.DN_DX, rData.N, rData.Area);
    CalculateNormal(rData.Normal);

    // Check that parents have been computed
    // These are required to retrieve the material properties and the viscous stress
    auto& r_neighbours = this->GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_ERROR_IF(r_neighbours.size() > 1) << "A condition was assigned more than one parent element." << std::endl;
    KRATOS_ERROR_IF(r_neighbours.size() == 0) << "A condition was not assigned a parent element. Please execute the check_and_prepare_model_process_fluid process." << std::endl;

    // Database access to all of the variables needed
    Properties& r_parent_properties = r_neighbours[0].GetProperties();
    rData.mu = r_parent_properties.GetValue(DYNAMIC_VISCOSITY);
    rData.c_v = r_parent_properties.GetValue(SPECIFIC_HEAT); // TODO: WE SHOULD SPECIFY WHICH ONE --> CREATE SPECIFIC_HEAT_CONSTANT_VOLUME
    rData.lambda = r_parent_properties.GetValue(CONDUCTIVITY);
    rData.gamma = r_parent_properties.GetValue(HEAT_CAPACITY_RATIO);
    // rData.ShockCapturing = rCurrentProcessInfo[SHOCK_CAPTURING_SWITCH];

    for (unsigned int i = 0; i < TNumNodes; ++i) {
        const auto& r_node = r_geometry[i];
        rData.U(i, 0) = r_node.FastGetSolutionStepValue(DENSITY);
        const auto& r_momentum = r_node.FastGetSolutionStepValue(MOMENTUM);
        for (unsigned int k = 0; k < TDim; ++k) {
            rData.U(i, k + 1) = r_momentum[k];
        }
        rData.U(i, TDim + 1) = r_node.FastGetSolutionStepValue(TOTAL_ENERGY);
        rData.beta_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_BULK_VISCOSITY);
        rData.lamb_sc_nodes(i) = r_node.GetValue(ARTIFICIAL_CONDUCTIVITY);
    }
}

template<>
void CompressibleNavierStokesWallCondition<2,2>::CalculateRightHandSideInternal(
    BoundedVector<double, 8>& rRightHandSideBoundedVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    constexpr unsigned int dim = 2;
    constexpr unsigned int n_nodes = 2;
    constexpr unsigned int block_size = 4;

    // Struct to pass around the data
    ConditionDataStruct data;
    this->FillConditionData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const double mu = data.mu;
    const double c_v = data.c_v;
    const double gamma = data.gamma;
    const double lambda = data.lambda;
    const array_1d<double, 3> normal = data.Normal;
    const BoundedMatrix<double, n_nodes, block_size>& U = data.U;
    const array_1d<double, n_nodes > N = data.N;
    const BoundedMatrix<double, n_nodes, dim> DN = data.DN_DX;
    const array_1d<double, n_nodes> &beta_sc_nodes = data.beta_sc_nodes;
    const array_1d<double, n_nodes> &lamb_sc_nodes = data.lamb_sc_nodes;

    const double crRightHandSideBoundedVector0 =             0.26794940881858814*U(0,0) + U(1,0);
const double crRightHandSideBoundedVector1 =             pow(crRightHandSideBoundedVector0, -2);
const double crRightHandSideBoundedVector2 =             0.21132500000000001*U(0,0) + 0.78867500000000001*U(1,0);
const double crRightHandSideBoundedVector3 =             DN(0,0)*U(0,2) + DN(1,0)*U(1,2);
const double crRightHandSideBoundedVector4 =             DN(0,1)*U(0,1) + DN(1,1)*U(1,1);
const double crRightHandSideBoundedVector5 =             0.21132500000000001*U(0,1) + 0.78867500000000001*U(1,1);
const double crRightHandSideBoundedVector6 =             DN(0,1)*U(0,0) + DN(1,1)*U(1,0);
const double crRightHandSideBoundedVector7 =             0.21132500000000001*U(0,2) + 0.78867500000000001*U(1,2);
const double crRightHandSideBoundedVector8 =             DN(0,0)*U(0,0) + DN(1,0)*U(1,0);
const double crRightHandSideBoundedVector9 =             crRightHandSideBoundedVector2*crRightHandSideBoundedVector3 + crRightHandSideBoundedVector2*crRightHandSideBoundedVector4 - crRightHandSideBoundedVector5*crRightHandSideBoundedVector6 - crRightHandSideBoundedVector7*crRightHandSideBoundedVector8;
const double crRightHandSideBoundedVector10 =             crRightHandSideBoundedVector1*crRightHandSideBoundedVector9;
const double crRightHandSideBoundedVector11 =             mu*normal[1];
const double crRightHandSideBoundedVector12 =             0.33974629450481902*crRightHandSideBoundedVector11;
const double crRightHandSideBoundedVector13 =             1.2679494088185881*crRightHandSideBoundedVector11;
const double crRightHandSideBoundedVector14 =             U(0,0) + 0.26794940881858814*U(1,0);
const double crRightHandSideBoundedVector15 =             pow(crRightHandSideBoundedVector14, -2);
const double crRightHandSideBoundedVector16 =             0.78867500000000001*U(0,0) + 0.21132500000000001*U(1,0);
const double crRightHandSideBoundedVector17 =             0.78867500000000001*U(0,1) + 0.21132500000000001*U(1,1);
const double crRightHandSideBoundedVector18 =             0.78867500000000001*U(0,2) + 0.21132500000000001*U(1,2);
const double crRightHandSideBoundedVector19 =             crRightHandSideBoundedVector16*crRightHandSideBoundedVector3 + crRightHandSideBoundedVector16*crRightHandSideBoundedVector4 - crRightHandSideBoundedVector17*crRightHandSideBoundedVector6 - crRightHandSideBoundedVector18*crRightHandSideBoundedVector8;
const double crRightHandSideBoundedVector20 =             crRightHandSideBoundedVector15*crRightHandSideBoundedVector19;
const double crRightHandSideBoundedVector21 =             DN(0,0)*U(0,1) + DN(1,0)*U(1,1);
const double crRightHandSideBoundedVector22 =             crRightHandSideBoundedVector2*crRightHandSideBoundedVector21 - crRightHandSideBoundedVector5*crRightHandSideBoundedVector8;
const double crRightHandSideBoundedVector23 =             3.2153914066468143*mu;
const double crRightHandSideBoundedVector24 =             -0.66666666666666663*mu;
const double crRightHandSideBoundedVector25 =             DN(0,1)*U(0,2) + DN(1,1)*U(1,2);
const double crRightHandSideBoundedVector26 =             crRightHandSideBoundedVector2*crRightHandSideBoundedVector25 - crRightHandSideBoundedVector6*crRightHandSideBoundedVector7;
const double crRightHandSideBoundedVector27 =             1.6076957033234072*(crRightHandSideBoundedVector22 + crRightHandSideBoundedVector26)*(0.21132500000000001*beta_sc_nodes[0] + 0.78867500000000001*beta_sc_nodes[1] + crRightHandSideBoundedVector24);
const double crRightHandSideBoundedVector28 =             crRightHandSideBoundedVector22*crRightHandSideBoundedVector23 + crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector29 =             0.21132500000000001*crRightHandSideBoundedVector1;
const double crRightHandSideBoundedVector30 =             crRightHandSideBoundedVector29*normal[0];
const double crRightHandSideBoundedVector31 =             crRightHandSideBoundedVector16*crRightHandSideBoundedVector21 - crRightHandSideBoundedVector17*crRightHandSideBoundedVector8;
const double crRightHandSideBoundedVector32 =             crRightHandSideBoundedVector16*crRightHandSideBoundedVector25 - crRightHandSideBoundedVector18*crRightHandSideBoundedVector6;
const double crRightHandSideBoundedVector33 =             1.6076957033234072*(crRightHandSideBoundedVector31 + crRightHandSideBoundedVector32)*(0.78867500000000001*beta_sc_nodes[0] + 0.21132500000000001*beta_sc_nodes[1] + crRightHandSideBoundedVector24);
const double crRightHandSideBoundedVector34 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector31 + crRightHandSideBoundedVector33;
const double crRightHandSideBoundedVector35 =             0.78867500000000001*crRightHandSideBoundedVector15;
const double crRightHandSideBoundedVector36 =             crRightHandSideBoundedVector35*normal[0];
const double crRightHandSideBoundedVector37 =             mu*normal[0];
const double crRightHandSideBoundedVector38 =             0.33974629450481902*crRightHandSideBoundedVector37;
const double crRightHandSideBoundedVector39 =             1.2679494088185881*crRightHandSideBoundedVector37;
const double crRightHandSideBoundedVector40 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector26 + crRightHandSideBoundedVector27;
const double crRightHandSideBoundedVector41 =             crRightHandSideBoundedVector29*normal[1];
const double crRightHandSideBoundedVector42 =             crRightHandSideBoundedVector23*crRightHandSideBoundedVector32 + crRightHandSideBoundedVector33;
const double crRightHandSideBoundedVector43 =             crRightHandSideBoundedVector35*normal[1];
const double crRightHandSideBoundedVector44 =             1.6076957033234072*crRightHandSideBoundedVector7;
const double crRightHandSideBoundedVector45 =             1.0/crRightHandSideBoundedVector2;
const double crRightHandSideBoundedVector46 =             crRightHandSideBoundedVector45*crRightHandSideBoundedVector9*mu;
const double crRightHandSideBoundedVector47 =             DN(0,0)*U(0,3) + DN(1,0)*U(1,3);
const double crRightHandSideBoundedVector48 =             1.6076957033234072*crRightHandSideBoundedVector2;
const double crRightHandSideBoundedVector49 =             1.6076957033234072*crRightHandSideBoundedVector5;
const double crRightHandSideBoundedVector50 =             0.33974629450481902*U(0,3) + 1.2679494088185881*U(1,3);
const double crRightHandSideBoundedVector51 =             pow(0.26794940881858814*U(0,1) + U(1,1), 2);
const double crRightHandSideBoundedVector52 =             1.2679494088185883/crRightHandSideBoundedVector0;
const double crRightHandSideBoundedVector53 =             crRightHandSideBoundedVector52*crRightHandSideBoundedVector8;
const double crRightHandSideBoundedVector54 =             pow(0.26794940881858814*U(0,2) + U(1,2), 2);
const double crRightHandSideBoundedVector55 =             1.0/c_v;
const double crRightHandSideBoundedVector56 =             crRightHandSideBoundedVector55*(0.21132500000000001*lamb_sc_nodes[0] + 0.78867500000000001*lamb_sc_nodes[1] + lambda);
const double crRightHandSideBoundedVector57 =             crRightHandSideBoundedVector28*crRightHandSideBoundedVector45*crRightHandSideBoundedVector5 + crRightHandSideBoundedVector44*crRightHandSideBoundedVector46 + crRightHandSideBoundedVector56*(-crRightHandSideBoundedVector21*crRightHandSideBoundedVector49 - crRightHandSideBoundedVector3*crRightHandSideBoundedVector44 + crRightHandSideBoundedVector47*crRightHandSideBoundedVector48 - crRightHandSideBoundedVector50*crRightHandSideBoundedVector8 + crRightHandSideBoundedVector51*crRightHandSideBoundedVector53 + crRightHandSideBoundedVector53*crRightHandSideBoundedVector54);
const double crRightHandSideBoundedVector58 =             1.6076957033234072*crRightHandSideBoundedVector18;
const double crRightHandSideBoundedVector59 =             1.0/crRightHandSideBoundedVector16;
const double crRightHandSideBoundedVector60 =             crRightHandSideBoundedVector19*crRightHandSideBoundedVector59*mu;
const double crRightHandSideBoundedVector61 =             1.6076957033234072*crRightHandSideBoundedVector16;
const double crRightHandSideBoundedVector62 =             1.6076957033234072*crRightHandSideBoundedVector17;
const double crRightHandSideBoundedVector63 =             1.2679494088185881*U(0,3) + 0.33974629450481902*U(1,3);
const double crRightHandSideBoundedVector64 =             pow(U(0,1) + 0.26794940881858814*U(1,1), 2);
const double crRightHandSideBoundedVector65 =             1.2679494088185883/crRightHandSideBoundedVector14;
const double crRightHandSideBoundedVector66 =             crRightHandSideBoundedVector65*crRightHandSideBoundedVector8;
const double crRightHandSideBoundedVector67 =             pow(U(0,2) + 0.26794940881858814*U(1,2), 2);
const double crRightHandSideBoundedVector68 =             crRightHandSideBoundedVector55*(0.78867500000000001*lamb_sc_nodes[0] + 0.21132500000000001*lamb_sc_nodes[1] + lambda);
const double crRightHandSideBoundedVector69 =             crRightHandSideBoundedVector17*crRightHandSideBoundedVector34*crRightHandSideBoundedVector59 + crRightHandSideBoundedVector58*crRightHandSideBoundedVector60 + crRightHandSideBoundedVector68*(-crRightHandSideBoundedVector21*crRightHandSideBoundedVector62 - crRightHandSideBoundedVector3*crRightHandSideBoundedVector58 + crRightHandSideBoundedVector47*crRightHandSideBoundedVector61 - crRightHandSideBoundedVector63*crRightHandSideBoundedVector8 + crRightHandSideBoundedVector64*crRightHandSideBoundedVector66 + crRightHandSideBoundedVector66*crRightHandSideBoundedVector67);
const double crRightHandSideBoundedVector70 =             DN(0,1)*U(0,3) + DN(1,1)*U(1,3);
const double crRightHandSideBoundedVector71 =             crRightHandSideBoundedVector52*crRightHandSideBoundedVector6;
const double crRightHandSideBoundedVector72 =             crRightHandSideBoundedVector40*crRightHandSideBoundedVector45*crRightHandSideBoundedVector7 + crRightHandSideBoundedVector46*crRightHandSideBoundedVector49 + crRightHandSideBoundedVector56*(-crRightHandSideBoundedVector25*crRightHandSideBoundedVector44 - crRightHandSideBoundedVector4*crRightHandSideBoundedVector49 + crRightHandSideBoundedVector48*crRightHandSideBoundedVector70 - crRightHandSideBoundedVector50*crRightHandSideBoundedVector6 + crRightHandSideBoundedVector51*crRightHandSideBoundedVector71 + crRightHandSideBoundedVector54*crRightHandSideBoundedVector71);
const double crRightHandSideBoundedVector73 =             crRightHandSideBoundedVector6*crRightHandSideBoundedVector65;
const double crRightHandSideBoundedVector74 =             crRightHandSideBoundedVector18*crRightHandSideBoundedVector42*crRightHandSideBoundedVector59 + crRightHandSideBoundedVector60*crRightHandSideBoundedVector62 + crRightHandSideBoundedVector68*(-crRightHandSideBoundedVector25*crRightHandSideBoundedVector58 - crRightHandSideBoundedVector4*crRightHandSideBoundedVector62 - crRightHandSideBoundedVector6*crRightHandSideBoundedVector63 + crRightHandSideBoundedVector61*crRightHandSideBoundedVector70 + crRightHandSideBoundedVector64*crRightHandSideBoundedVector73 + crRightHandSideBoundedVector67*crRightHandSideBoundedVector73);
const double crRightHandSideBoundedVector75 =             0.78867500000000001*crRightHandSideBoundedVector1;
const double crRightHandSideBoundedVector76 =             crRightHandSideBoundedVector75*normal[0];
const double crRightHandSideBoundedVector77 =             0.21132500000000001*crRightHandSideBoundedVector15;
const double crRightHandSideBoundedVector78 =             crRightHandSideBoundedVector77*normal[0];
const double crRightHandSideBoundedVector79 =             crRightHandSideBoundedVector75*normal[1];
const double crRightHandSideBoundedVector80 =             crRightHandSideBoundedVector77*normal[1];
            rRightHandSideBoundedVector[0]=0;
            rRightHandSideBoundedVector[1]=crRightHandSideBoundedVector10*crRightHandSideBoundedVector12 + crRightHandSideBoundedVector13*crRightHandSideBoundedVector20 + crRightHandSideBoundedVector28*crRightHandSideBoundedVector30 + crRightHandSideBoundedVector34*crRightHandSideBoundedVector36;
            rRightHandSideBoundedVector[2]=crRightHandSideBoundedVector10*crRightHandSideBoundedVector38 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector39 + crRightHandSideBoundedVector40*crRightHandSideBoundedVector41 + crRightHandSideBoundedVector42*crRightHandSideBoundedVector43;
            rRightHandSideBoundedVector[3]=crRightHandSideBoundedVector30*crRightHandSideBoundedVector57 + crRightHandSideBoundedVector36*crRightHandSideBoundedVector69 + crRightHandSideBoundedVector41*crRightHandSideBoundedVector72 + crRightHandSideBoundedVector43*crRightHandSideBoundedVector74;
            rRightHandSideBoundedVector[4]=0;
            rRightHandSideBoundedVector[5]=crRightHandSideBoundedVector10*crRightHandSideBoundedVector13 + crRightHandSideBoundedVector12*crRightHandSideBoundedVector20 + crRightHandSideBoundedVector28*crRightHandSideBoundedVector76 + crRightHandSideBoundedVector34*crRightHandSideBoundedVector78;
            rRightHandSideBoundedVector[6]=crRightHandSideBoundedVector10*crRightHandSideBoundedVector39 + crRightHandSideBoundedVector20*crRightHandSideBoundedVector38 + crRightHandSideBoundedVector40*crRightHandSideBoundedVector79 + crRightHandSideBoundedVector42*crRightHandSideBoundedVector80;
            rRightHandSideBoundedVector[7]=crRightHandSideBoundedVector57*crRightHandSideBoundedVector76 + crRightHandSideBoundedVector69*crRightHandSideBoundedVector78 + crRightHandSideBoundedVector72*crRightHandSideBoundedVector79 + crRightHandSideBoundedVector74*crRightHandSideBoundedVector80;


    data.Area / static_cast<double>(n_nodes);
}

template<>
void CompressibleNavierStokesWallCondition<3,3>::CalculateRightHandSideInternal(
    BoundedVector<double, 15>& rRightHandSideBoundedVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    constexpr unsigned int dim = 3;
    constexpr unsigned int n_nodes = 3;
    constexpr unsigned int block_size = 5;

    // Struct to pass around the data
    ConditionDataStruct data;
    this->FillConditionData(data, rCurrentProcessInfo);

    // Substitute the formulation symbols by the data structure values
    const double mu = data.mu;
    const double c_v = data.c_v;
    const double gamma = data.gamma;
    const double lambda = data.lambda;
    const array_1d<double, 3> normal = data.Normal;
    const BoundedMatrix<double, n_nodes, block_size>& U = data.U;
    const array_1d<double, n_nodes > N = data.N;
    const BoundedMatrix<double, n_nodes, dim> DN = data.DN_DX;
    const array_1d<double, n_nodes> &beta_sc_nodes = data.beta_sc_nodes;
    const array_1d<double, n_nodes> &lamb_sc_nodes = data.lamb_sc_nodes;

    //substitute_rhs_3D

    data.Area / static_cast<double>(n_nodes);
}

template <>
void CompressibleNavierStokesWallCondition<2,2>::AddExplicitContribution(const ProcessInfo &rCurrentProcessInfo)
{
    constexpr IndexType dim = 2;
    constexpr IndexType n_nodes = 2;
    constexpr IndexType block_size = 4;

    // Calculate the explicit residual vector
    BoundedVector<double, n_nodes * block_size> rhs;
    CalculateRightHandSideInternal(rhs, rCurrentProcessInfo);

    // Add the residual contribution
    // Note that the reaction is indeed the formulation residual
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const IndexType aux = i_node * block_size;
#pragma omp atomic
        r_geometry[i_node].FastGetSolutionStepValue(REACTION_DENSITY) += rhs[aux];
        auto& r_mom = r_geometry[i_node].FastGetSolutionStepValue(REACTION);
        for (IndexType d = 0; d < dim; ++d) {
#pragma omp atomic
            r_mom[d] += rhs[aux + (d + 1)];
        }
#pragma omp atomic
        r_geometry[i_node].FastGetSolutionStepValue(REACTION_ENERGY) += rhs[aux + 3];
    }
}

template <>
void CompressibleNavierStokesWallCondition<3,3>::AddExplicitContribution(const ProcessInfo &rCurrentProcessInfo)
{
    constexpr IndexType dim = 3;
    constexpr IndexType n_nodes = 3;
    constexpr IndexType block_size = 5;

    // Calculate the explicit residual vector
    BoundedVector<double, n_nodes * block_size> rhs;
    CalculateRightHandSideInternal(rhs, rCurrentProcessInfo);

    // Add the residual contribution
    // Note that the reaction is indeed the formulation residual
    auto& r_geometry = GetGeometry();
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const IndexType aux = i_node * block_size;
#pragma omp atomic
        r_geometry[i_node].FastGetSolutionStepValue(REACTION_DENSITY) += rhs[aux];
        auto& r_mom = r_geometry[i_node].FastGetSolutionStepValue(REACTION);
        for (IndexType d = 0; d < dim; ++d) {
#pragma omp atomic
            r_mom[d] += rhs[aux + (d + 1)];
        }
#pragma omp atomic
        r_geometry[i_node].FastGetSolutionStepValue(REACTION_ENERGY) += rhs[aux + 4];
    }
}

template <>
void CompressibleNavierStokesWallCondition<2,2>::CalculateNormal(array_1d<double,3>& rNormal) const
{
    const auto& r_geometry = this->GetGeometry();
    rNormal[0] = r_geometry[1].Y() - r_geometry[0].Y();
    rNormal[1] = -(r_geometry[1].X() - r_geometry[0].X());
    rNormal[2] = 0.0;
}

template <>
void CompressibleNavierStokesWallCondition<3,3>::CalculateNormal(array_1d<double,3>& rNormal) const
{
    array_1d<double,3> v1,v2;
    const auto& r_geometry = this->GetGeometry();

    v1[0] = r_geometry[1].X() - r_geometry[0].X();
    v1[1] = r_geometry[1].Y() - r_geometry[0].Y();
    v1[2] = r_geometry[1].Z() - r_geometry[0].Z();

    v2[0] = r_geometry[2].X() - r_geometry[0].X();
    v2[1] = r_geometry[2].Y() - r_geometry[0].Y();
    v2[2] = r_geometry[2].Z() - r_geometry[0].Z();

    MathUtils<double>::CrossProduct(rNormal, v1, v2);
    rNormal *= 0.5;
}

template class CompressibleNavierStokesWallCondition<2,2>;
template class CompressibleNavierStokesWallCondition<3,3>;

} // namespace Kratos
