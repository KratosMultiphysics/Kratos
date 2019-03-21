//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla, Mengjie Zhao
//

#include "custom_elements/time_averaged_navier_stokes.h"

namespace Kratos {

template<>
void TimeAveragedNavierStokes<3>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int dim = 3;
    constexpr unsigned int num_nodes = 4;
    constexpr unsigned int dof_size = num_nodes*(dim+1);

    if (rResult.size() != dof_size){
        rResult.resize(dof_size, false);
    }

    for(unsigned int i=0; i<num_nodes; i++){
        rResult[i*(dim+1)  ]  =  this->GetGeometry()[i].GetDof(TIME_AVERAGED_VELOCITY_X).EquationId();
        rResult[i*(dim+1)+1]  =  this->GetGeometry()[i].GetDof(TIME_AVERAGED_VELOCITY_Y).EquationId();
        rResult[i*(dim+1)+2]  =  this->GetGeometry()[i].GetDof(TIME_AVERAGED_VELOCITY_Z).EquationId();
        rResult[i*(dim+1)+3]  =  this->GetGeometry()[i].GetDof(TIME_AVERAGED_PRESSURE).EquationId();
    }

    KRATOS_CATCH("ERROR in constrcuting EquationIdVector")
}


template<>
void TimeAveragedNavierStokes<2>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int dim = 2;
    constexpr unsigned int num_nodes = 3;
    constexpr unsigned int dof_size  = num_nodes*(dim+1);

    if (rResult.size() != dof_size){
        rResult.resize(dof_size, false);
    }

    for(unsigned int i=0; i<num_nodes; i++){
        rResult[i*(dim+1)  ]  =  this->GetGeometry()[i].GetDof(TIME_AVERAGED_VELOCITY_X).EquationId();
        rResult[i*(dim+1)+1]  =  this->GetGeometry()[i].GetDof(TIME_AVERAGED_VELOCITY_Y).EquationId();
        rResult[i*(dim+1)+2]  =  this->GetGeometry()[i].GetDof(TIME_AVERAGED_PRESSURE).EquationId();
    }

    KRATOS_CATCH("ERROR in constrcuting EquationIdVector")
}


template<>
void TimeAveragedNavierStokes<3>::GetDofList(
    DofsVectorType& ElementalDofList,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int dim = 3;
    constexpr unsigned int num_nodes = 4;
    constexpr unsigned int dof_size  = num_nodes*(dim+1);

    if (ElementalDofList.size() != dof_size){
        ElementalDofList.resize(dof_size);
    }

    for(unsigned int i=0; i<num_nodes; i++){
        ElementalDofList[i*(dim+1)  ]  =  this->GetGeometry()[i].pGetDof(TIME_AVERAGED_VELOCITY_X);
        ElementalDofList[i*(dim+1)+1]  =  this->GetGeometry()[i].pGetDof(TIME_AVERAGED_VELOCITY_Y);
        ElementalDofList[i*(dim+1)+2]  =  this->GetGeometry()[i].pGetDof(TIME_AVERAGED_VELOCITY_Z);
        ElementalDofList[i*(dim+1)+3]  =  this->GetGeometry()[i].pGetDof(TIME_AVERAGED_PRESSURE);
    }
    KRATOS_CATCH("ERROR in constrcuting GetDofList")
}


template<>
void TimeAveragedNavierStokes<2>::GetDofList(
    DofsVectorType& ElementalDofList,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int dim = 2;
    constexpr unsigned int num_nodes = 3;
    constexpr unsigned int dof_size  = num_nodes*(dim+1);

    if (ElementalDofList.size() != dof_size){
        ElementalDofList.resize(dof_size);
    }

    for(unsigned int i=0; i<num_nodes; i++){
        ElementalDofList[i*(dim+1)  ]  =  this->GetGeometry()[i].pGetDof(TIME_AVERAGED_VELOCITY_X);
        ElementalDofList[i*(dim+1)+1]  =  this->GetGeometry()[i].pGetDof(TIME_AVERAGED_VELOCITY_Y);
        ElementalDofList[i*(dim+1)+2]  =  this->GetGeometry()[i].pGetDof(TIME_AVERAGED_PRESSURE);
    }

    KRATOS_CATCH("ERROR in constrcuting GetDofList")
}


template<>
void TimeAveragedNavierStokes<3>::ComputeGaussPointLHSContribution(
    BoundedMatrix<double,16,16>& lhs,
    const ElementDataStruct& data)
{
    constexpr int dim = 3;
    constexpr int nnodes = 4;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size
    const double c = data.c;                                // Wave velocity

    //const double& dts = data.dts;                           // The averaging time period
    const double& dtn = data.dtn;                           // Time increment: notice t = tn + dtn

    const double& bdf0 = data.bdf0;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v_ave = data.v_ave;
    const BoundedMatrix<double,nnodes,dim>& vn_ave = data.vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;

    // time averaging parameters
    double ave_c1 = data.ave_c1;
    double ave_c2 = data.ave_c2;

    // get time accurate expression
    const BoundedMatrix<double,nnodes,dim>& v = ave_c1 * v_ave - ave_c2 * vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vconv = v - vmesh;

    // Get constitutive matrix
    const Matrix& C = data.C;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    const double clhs0 =             C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs1 =             C(0,3)*DN(0,0);
const double clhs2 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs1;
const double clhs3 =             C(0,5)*DN(0,0);
const double clhs4 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs3;
const double clhs5 =             pow(DN(0,0), 2);
const double clhs6 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double clhs7 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double clhs8 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double clhs9 =             rho*stab_c2*sqrt(pow(clhs6, 2) + pow(clhs7, 2) + pow(clhs8, 2));
const double clhs10 =             clhs9*h/stab_c1 + mu;
const double clhs11 =             DN(0,0)*clhs6 + DN(0,1)*clhs7 + DN(0,2)*clhs8;
const double clhs12 =             N[0]*rho;
const double clhs13 =             pow(N[0], 2);
const double clhs14 =             bdf0*rho;
const double clhs15 =             N[0]*bdf0;
const double clhs16 =             clhs11 + clhs15;
const double clhs17 =             1.0*h/(clhs9 + mu*stab_c1/h);
const double clhs18 =             clhs17*pow(rho, 2);
const double clhs19 =             clhs11*clhs18;
const double clhs20 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double clhs21 =             clhs18*clhs20;
const double clhs22 =             N[0]*clhs21;
const double clhs23 =             clhs11*clhs12 + clhs13*clhs14 + clhs16*clhs19 + clhs16*clhs22;
const double clhs24 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs1;
const double clhs25 =             C(1,3)*DN(0,1);
const double clhs26 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs25;
const double clhs27 =             C(3,5)*DN(0,0);
const double clhs28 =             C(4,5)*DN(0,2);
const double clhs29 =             C(1,5)*DN(0,1) + clhs27 + clhs28;
const double clhs30 =             DN(0,0)*clhs10;
const double clhs31 =             DN(0,1)*clhs30;
const double clhs32 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs3;
const double clhs33 =             C(3,4)*DN(0,1);
const double clhs34 =             C(2,3)*DN(0,2) + clhs27 + clhs33;
const double clhs35 =             C(2,5)*DN(0,2);
const double clhs36 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs35;
const double clhs37 =             DN(0,2)*clhs30;
const double clhs38 =             1/(pow(c, 2)*rho);
const double clhs39 =             clhs15*clhs38;
const double clhs40 =             clhs10*clhs39;
const double clhs41 =             clhs17*clhs20;
const double clhs42 =             clhs17*rho;
const double clhs43 =             clhs11*clhs42;
const double clhs44 =             ave_c1*(-N[0] + clhs12*clhs41 + clhs40 + clhs43);
const double clhs45 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs46 =             C(0,3)*DN(1,0);
const double clhs47 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs46;
const double clhs48 =             C(0,5)*DN(1,0);
const double clhs49 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs48;
const double clhs50 =             N[1]*rho;
const double clhs51 =             clhs15*clhs50;
const double clhs52 =             DN(1,0)*clhs30 + clhs51;
const double clhs53 =             DN(1,0)*clhs6 + DN(1,1)*clhs7 + DN(1,2)*clhs8;
const double clhs54 =             N[1]*bdf0;
const double clhs55 =             clhs53 + clhs54;
const double clhs56 =             clhs12*clhs53 + clhs19*clhs55 + clhs22*clhs55;
const double clhs57 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs46;
const double clhs58 =             C(1,3)*DN(1,1);
const double clhs59 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs58;
const double clhs60 =             C(3,5)*DN(1,0);
const double clhs61 =             C(4,5)*DN(1,2);
const double clhs62 =             C(1,5)*DN(1,1) + clhs60 + clhs61;
const double clhs63 =             DN(1,1)*clhs30;
const double clhs64 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs48;
const double clhs65 =             C(3,4)*DN(1,1);
const double clhs66 =             C(2,3)*DN(1,2) + clhs60 + clhs65;
const double clhs67 =             C(2,5)*DN(1,2);
const double clhs68 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs67;
const double clhs69 =             DN(1,2)*clhs30;
const double clhs70 =             DN(0,0)*N[1];
const double clhs71 =             clhs38*clhs54;
const double clhs72 =             DN(1,0)*N[0];
const double clhs73 =             clhs20*clhs42;
const double clhs74 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs75 =             C(0,3)*DN(2,0);
const double clhs76 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs75;
const double clhs77 =             C(0,5)*DN(2,0);
const double clhs78 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs77;
const double clhs79 =             N[2]*rho;
const double clhs80 =             clhs15*clhs79;
const double clhs81 =             DN(2,0)*clhs30 + clhs80;
const double clhs82 =             DN(2,0)*clhs6 + DN(2,1)*clhs7 + DN(2,2)*clhs8;
const double clhs83 =             N[2]*bdf0;
const double clhs84 =             clhs82 + clhs83;
const double clhs85 =             clhs12*clhs82 + clhs19*clhs84 + clhs22*clhs84;
const double clhs86 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs75;
const double clhs87 =             C(1,3)*DN(2,1);
const double clhs88 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs87;
const double clhs89 =             C(3,5)*DN(2,0);
const double clhs90 =             C(4,5)*DN(2,2);
const double clhs91 =             C(1,5)*DN(2,1) + clhs89 + clhs90;
const double clhs92 =             DN(2,1)*clhs30;
const double clhs93 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs77;
const double clhs94 =             C(3,4)*DN(2,1);
const double clhs95 =             C(2,3)*DN(2,2) + clhs89 + clhs94;
const double clhs96 =             C(2,5)*DN(2,2);
const double clhs97 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs96;
const double clhs98 =             DN(2,2)*clhs30;
const double clhs99 =             DN(0,0)*N[2];
const double clhs100 =             clhs38*clhs83;
const double clhs101 =             DN(2,0)*N[0];
const double clhs102 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs103 =             C(0,3)*DN(3,0);
const double clhs104 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs103;
const double clhs105 =             C(0,5)*DN(3,0);
const double clhs106 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs105;
const double clhs107 =             N[3]*rho;
const double clhs108 =             clhs107*clhs15;
const double clhs109 =             DN(3,0)*clhs30 + clhs108;
const double clhs110 =             DN(3,0)*clhs6 + DN(3,1)*clhs7 + DN(3,2)*clhs8;
const double clhs111 =             N[3]*bdf0;
const double clhs112 =             clhs110 + clhs111;
const double clhs113 =             clhs110*clhs12 + clhs112*clhs19 + clhs112*clhs22;
const double clhs114 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs103;
const double clhs115 =             C(1,3)*DN(3,1);
const double clhs116 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs115;
const double clhs117 =             C(3,5)*DN(3,0);
const double clhs118 =             C(4,5)*DN(3,2);
const double clhs119 =             C(1,5)*DN(3,1) + clhs117 + clhs118;
const double clhs120 =             DN(3,1)*clhs30;
const double clhs121 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs105;
const double clhs122 =             C(3,4)*DN(3,1);
const double clhs123 =             C(2,3)*DN(3,2) + clhs117 + clhs122;
const double clhs124 =             C(2,5)*DN(3,2);
const double clhs125 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs124;
const double clhs126 =             DN(3,2)*clhs30;
const double clhs127 =             DN(0,0)*N[3];
const double clhs128 =             clhs111*clhs38;
const double clhs129 =             DN(3,0)*N[0];
const double clhs130 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs25;
const double clhs131 =             C(0,4)*DN(0,0) + clhs28 + clhs33;
const double clhs132 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs133 =             C(1,4)*DN(0,1);
const double clhs134 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs133;
const double clhs135 =             pow(DN(0,1), 2);
const double clhs136 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs133;
const double clhs137 =             C(2,4)*DN(0,2);
const double clhs138 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs137;
const double clhs139 =             DN(0,1)*clhs10;
const double clhs140 =             DN(0,2)*clhs139;
const double clhs141 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs58;
const double clhs142 =             C(0,4)*DN(1,0) + clhs61 + clhs65;
const double clhs143 =             DN(1,0)*clhs139;
const double clhs144 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs145 =             C(1,4)*DN(1,1);
const double clhs146 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs145;
const double clhs147 =             DN(1,1)*clhs139;
const double clhs148 =             clhs51 + clhs56;
const double clhs149 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs145;
const double clhs150 =             C(2,4)*DN(1,2);
const double clhs151 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs150;
const double clhs152 =             DN(1,2)*clhs139;
const double clhs153 =             DN(0,1)*N[1];
const double clhs154 =             DN(1,1)*N[0];
const double clhs155 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs87;
const double clhs156 =             C(0,4)*DN(2,0) + clhs90 + clhs94;
const double clhs157 =             DN(2,0)*clhs139;
const double clhs158 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs159 =             C(1,4)*DN(2,1);
const double clhs160 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs159;
const double clhs161 =             DN(2,1)*clhs139;
const double clhs162 =             clhs80 + clhs85;
const double clhs163 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs159;
const double clhs164 =             C(2,4)*DN(2,2);
const double clhs165 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs164;
const double clhs166 =             DN(2,2)*clhs139;
const double clhs167 =             DN(0,1)*N[2];
const double clhs168 =             DN(2,1)*N[0];
const double clhs169 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs115;
const double clhs170 =             C(0,4)*DN(3,0) + clhs118 + clhs122;
const double clhs171 =             DN(3,0)*clhs139;
const double clhs172 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs173 =             C(1,4)*DN(3,1);
const double clhs174 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs173;
const double clhs175 =             DN(3,1)*clhs139;
const double clhs176 =             clhs108 + clhs113;
const double clhs177 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs173;
const double clhs178 =             C(2,4)*DN(3,2);
const double clhs179 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs178;
const double clhs180 =             DN(3,2)*clhs139;
const double clhs181 =             DN(0,1)*N[3];
const double clhs182 =             DN(3,1)*N[0];
const double clhs183 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs35;
const double clhs184 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs137;
const double clhs185 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs186 =             pow(DN(0,2), 2);
const double clhs187 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs67;
const double clhs188 =             DN(0,2)*clhs10;
const double clhs189 =             DN(1,0)*clhs188;
const double clhs190 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs150;
const double clhs191 =             DN(1,1)*clhs188;
const double clhs192 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs193 =             DN(1,2)*clhs188;
const double clhs194 =             DN(0,2)*N[1];
const double clhs195 =             DN(1,2)*N[0];
const double clhs196 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs96;
const double clhs197 =             DN(2,0)*clhs188;
const double clhs198 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs164;
const double clhs199 =             DN(2,1)*clhs188;
const double clhs200 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs201 =             DN(2,2)*clhs188;
const double clhs202 =             DN(0,2)*N[2];
const double clhs203 =             DN(2,2)*N[0];
const double clhs204 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs124;
const double clhs205 =             DN(3,0)*clhs188;
const double clhs206 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs178;
const double clhs207 =             DN(3,1)*clhs188;
const double clhs208 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs209 =             DN(3,2)*clhs188;
const double clhs210 =             DN(0,2)*N[3];
const double clhs211 =             DN(3,2)*N[0];
const double clhs212 =             clhs16*clhs42;
const double clhs213 =             ave_c1*(N[0] + clhs212);
const double clhs214 =             bdf0*clhs38;
const double clhs215 =             clhs42*clhs55;
const double clhs216 =             DN(0,0)*clhs17;
const double clhs217 =             DN(0,1)*clhs17;
const double clhs218 =             DN(0,2)*clhs17;
const double clhs219 =             ave_c1*(DN(1,0)*clhs216 + DN(1,1)*clhs217 + DN(1,2)*clhs218 + N[1]*clhs39);
const double clhs220 =             clhs42*clhs84;
const double clhs221 =             ave_c1*(DN(2,0)*clhs216 + DN(2,1)*clhs217 + DN(2,2)*clhs218 + N[2]*clhs39);
const double clhs222 =             clhs112*clhs42;
const double clhs223 =             ave_c1*(DN(3,0)*clhs216 + DN(3,1)*clhs217 + DN(3,2)*clhs218 + N[3]*clhs39);
const double clhs224 =             clhs18*clhs53;
const double clhs225 =             N[1]*clhs21;
const double clhs226 =             clhs11*clhs50 + clhs16*clhs224 + clhs16*clhs225;
const double clhs227 =             DN(1,0)*clhs10;
const double clhs228 =             clhs42*clhs53;
const double clhs229 =             pow(DN(1,0), 2);
const double clhs230 =             pow(N[1], 2);
const double clhs231 =             clhs14*clhs230 + clhs224*clhs55 + clhs225*clhs55 + clhs50*clhs53;
const double clhs232 =             DN(1,1)*clhs227;
const double clhs233 =             DN(1,2)*clhs227;
const double clhs234 =             clhs10*clhs71;
const double clhs235 =             ave_c1*(-N[1] + clhs228 + clhs234 + clhs41*clhs50);
const double clhs236 =             clhs54*clhs79;
const double clhs237 =             DN(2,0)*clhs227 + clhs236;
const double clhs238 =             clhs224*clhs84 + clhs225*clhs84 + clhs50*clhs82;
const double clhs239 =             DN(2,1)*clhs227;
const double clhs240 =             DN(2,2)*clhs227;
const double clhs241 =             DN(1,0)*N[2];
const double clhs242 =             DN(2,0)*N[1];
const double clhs243 =             clhs107*clhs54;
const double clhs244 =             DN(3,0)*clhs227 + clhs243;
const double clhs245 =             clhs110*clhs50 + clhs112*clhs224 + clhs112*clhs225;
const double clhs246 =             DN(3,1)*clhs227;
const double clhs247 =             DN(3,2)*clhs227;
const double clhs248 =             DN(1,0)*N[3];
const double clhs249 =             DN(3,0)*N[1];
const double clhs250 =             clhs226 + clhs51;
const double clhs251 =             DN(1,1)*clhs10;
const double clhs252 =             pow(DN(1,1), 2);
const double clhs253 =             DN(1,2)*clhs251;
const double clhs254 =             DN(2,0)*clhs251;
const double clhs255 =             DN(2,1)*clhs251;
const double clhs256 =             clhs236 + clhs238;
const double clhs257 =             DN(2,2)*clhs251;
const double clhs258 =             DN(1,1)*N[2];
const double clhs259 =             DN(2,1)*N[1];
const double clhs260 =             DN(3,0)*clhs251;
const double clhs261 =             DN(3,1)*clhs251;
const double clhs262 =             clhs243 + clhs245;
const double clhs263 =             DN(3,2)*clhs251;
const double clhs264 =             DN(1,1)*N[3];
const double clhs265 =             DN(3,1)*N[1];
const double clhs266 =             DN(1,2)*clhs10;
const double clhs267 =             pow(DN(1,2), 2);
const double clhs268 =             DN(2,0)*clhs266;
const double clhs269 =             DN(2,1)*clhs266;
const double clhs270 =             DN(2,2)*clhs266;
const double clhs271 =             DN(1,2)*N[2];
const double clhs272 =             DN(2,2)*N[1];
const double clhs273 =             DN(3,0)*clhs266;
const double clhs274 =             DN(3,1)*clhs266;
const double clhs275 =             DN(3,2)*clhs266;
const double clhs276 =             DN(1,2)*N[3];
const double clhs277 =             DN(3,2)*N[1];
const double clhs278 =             ave_c1*(N[1] + clhs215);
const double clhs279 =             DN(1,0)*clhs17;
const double clhs280 =             DN(1,1)*clhs17;
const double clhs281 =             DN(1,2)*clhs17;
const double clhs282 =             ave_c1*(DN(2,0)*clhs279 + DN(2,1)*clhs280 + DN(2,2)*clhs281 + N[2]*clhs71);
const double clhs283 =             ave_c1*(DN(3,0)*clhs279 + DN(3,1)*clhs280 + DN(3,2)*clhs281 + N[3]*clhs71);
const double clhs284 =             clhs18*clhs82;
const double clhs285 =             N[2]*clhs21;
const double clhs286 =             clhs11*clhs79 + clhs16*clhs284 + clhs16*clhs285;
const double clhs287 =             DN(2,0)*clhs10;
const double clhs288 =             clhs42*clhs82;
const double clhs289 =             clhs284*clhs55 + clhs285*clhs55 + clhs53*clhs79;
const double clhs290 =             pow(DN(2,0), 2);
const double clhs291 =             pow(N[2], 2);
const double clhs292 =             clhs14*clhs291 + clhs284*clhs84 + clhs285*clhs84 + clhs79*clhs82;
const double clhs293 =             DN(2,1)*clhs287;
const double clhs294 =             DN(2,2)*clhs287;
const double clhs295 =             clhs10*clhs100;
const double clhs296 =             ave_c1*(-N[2] + clhs288 + clhs295 + clhs41*clhs79);
const double clhs297 =             clhs107*clhs83;
const double clhs298 =             DN(3,0)*clhs287 + clhs297;
const double clhs299 =             clhs110*clhs79 + clhs112*clhs284 + clhs112*clhs285;
const double clhs300 =             DN(3,1)*clhs287;
const double clhs301 =             DN(3,2)*clhs287;
const double clhs302 =             DN(2,0)*N[3];
const double clhs303 =             DN(3,0)*N[2];
const double clhs304 =             clhs286 + clhs80;
const double clhs305 =             DN(2,1)*clhs10;
const double clhs306 =             clhs236 + clhs289;
const double clhs307 =             pow(DN(2,1), 2);
const double clhs308 =             DN(2,2)*clhs305;
const double clhs309 =             DN(3,0)*clhs305;
const double clhs310 =             DN(3,1)*clhs305;
const double clhs311 =             clhs297 + clhs299;
const double clhs312 =             DN(3,2)*clhs305;
const double clhs313 =             DN(2,1)*N[3];
const double clhs314 =             DN(3,1)*N[2];
const double clhs315 =             DN(2,2)*clhs10;
const double clhs316 =             pow(DN(2,2), 2);
const double clhs317 =             DN(3,0)*clhs315;
const double clhs318 =             DN(3,1)*clhs315;
const double clhs319 =             DN(3,2)*clhs315;
const double clhs320 =             DN(2,2)*N[3];
const double clhs321 =             DN(3,2)*N[2];
const double clhs322 =             ave_c1*(N[2] + clhs220);
const double clhs323 =             ave_c1*(DN(2,0)*DN(3,0)*clhs17 + DN(2,1)*DN(3,1)*clhs17 + DN(2,2)*DN(3,2)*clhs17 + N[3]*clhs100);
const double clhs324 =             clhs110*clhs18;
const double clhs325 =             N[3]*clhs21;
const double clhs326 =             clhs107*clhs11 + clhs16*clhs324 + clhs16*clhs325;
const double clhs327 =             DN(3,0)*clhs10;
const double clhs328 =             clhs110*clhs42;
const double clhs329 =             clhs107*clhs53 + clhs324*clhs55 + clhs325*clhs55;
const double clhs330 =             clhs107*clhs82 + clhs324*clhs84 + clhs325*clhs84;
const double clhs331 =             pow(DN(3,0), 2);
const double clhs332 =             pow(N[3], 2);
const double clhs333 =             clhs107*clhs110 + clhs112*clhs324 + clhs112*clhs325 + clhs14*clhs332;
const double clhs334 =             DN(3,1)*clhs327;
const double clhs335 =             DN(3,2)*clhs327;
const double clhs336 =             ave_c1*(-N[3] + clhs10*clhs128 + clhs107*clhs41 + clhs328);
const double clhs337 =             clhs108 + clhs326;
const double clhs338 =             DN(3,1)*clhs10;
const double clhs339 =             clhs243 + clhs329;
const double clhs340 =             clhs297 + clhs330;
const double clhs341 =             pow(DN(3,1), 2);
const double clhs342 =             DN(3,2)*clhs338;
const double clhs343 =             pow(DN(3,2), 2);
const double clhs344 =             ave_c1*(N[3] + clhs222);
            lhs(0,0)=ave_c1*(DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs10*clhs5 + clhs23);
            lhs(0,1)=ave_c1*(DN(0,0)*clhs24 + DN(0,1)*clhs26 + DN(0,2)*clhs29 + clhs31);
            lhs(0,2)=ave_c1*(DN(0,0)*clhs32 + DN(0,1)*clhs34 + DN(0,2)*clhs36 + clhs37);
            lhs(0,3)=DN(0,0)*clhs44;
            lhs(0,4)=ave_c1*(DN(0,0)*clhs45 + DN(0,1)*clhs47 + DN(0,2)*clhs49 + clhs52 + clhs56);
            lhs(0,5)=ave_c1*(DN(0,0)*clhs57 + DN(0,1)*clhs59 + DN(0,2)*clhs62 + clhs63);
            lhs(0,6)=ave_c1*(DN(0,0)*clhs64 + DN(0,1)*clhs66 + DN(0,2)*clhs68 + clhs69);
            lhs(0,7)=ave_c1*(DN(1,0)*clhs43 + clhs30*clhs71 - clhs70 + clhs72*clhs73);
            lhs(0,8)=ave_c1*(DN(0,0)*clhs74 + DN(0,1)*clhs76 + DN(0,2)*clhs78 + clhs81 + clhs85);
            lhs(0,9)=ave_c1*(DN(0,0)*clhs86 + DN(0,1)*clhs88 + DN(0,2)*clhs91 + clhs92);
            lhs(0,10)=ave_c1*(DN(0,0)*clhs93 + DN(0,1)*clhs95 + DN(0,2)*clhs97 + clhs98);
            lhs(0,11)=ave_c1*(DN(2,0)*clhs43 + clhs100*clhs30 + clhs101*clhs73 - clhs99);
            lhs(0,12)=ave_c1*(DN(0,0)*clhs102 + DN(0,1)*clhs104 + DN(0,2)*clhs106 + clhs109 + clhs113);
            lhs(0,13)=ave_c1*(DN(0,0)*clhs114 + DN(0,1)*clhs116 + DN(0,2)*clhs119 + clhs120);
            lhs(0,14)=ave_c1*(DN(0,0)*clhs121 + DN(0,1)*clhs123 + DN(0,2)*clhs125 + clhs126);
            lhs(0,15)=ave_c1*(DN(3,0)*clhs43 - clhs127 + clhs128*clhs30 + clhs129*clhs73);
            lhs(1,0)=ave_c1*(DN(0,0)*clhs2 + DN(0,1)*clhs130 + DN(0,2)*clhs131 + clhs31);
            lhs(1,1)=ave_c1*(DN(0,0)*clhs26 + DN(0,1)*clhs132 + DN(0,2)*clhs134 + clhs10*clhs135 + clhs23);
            lhs(1,2)=ave_c1*(DN(0,0)*clhs34 + DN(0,1)*clhs136 + DN(0,2)*clhs138 + clhs140);
            lhs(1,3)=DN(0,1)*clhs44;
            lhs(1,4)=ave_c1*(DN(0,0)*clhs47 + DN(0,1)*clhs141 + DN(0,2)*clhs142 + clhs143);
            lhs(1,5)=ave_c1*(DN(0,0)*clhs59 + DN(0,1)*clhs144 + DN(0,2)*clhs146 + clhs147 + clhs148);
            lhs(1,6)=ave_c1*(DN(0,0)*clhs66 + DN(0,1)*clhs149 + DN(0,2)*clhs151 + clhs152);
            lhs(1,7)=ave_c1*(DN(1,1)*clhs43 + clhs139*clhs71 - clhs153 + clhs154*clhs73);
            lhs(1,8)=ave_c1*(DN(0,0)*clhs76 + DN(0,1)*clhs155 + DN(0,2)*clhs156 + clhs157);
            lhs(1,9)=ave_c1*(DN(0,0)*clhs88 + DN(0,1)*clhs158 + DN(0,2)*clhs160 + clhs161 + clhs162);
            lhs(1,10)=ave_c1*(DN(0,0)*clhs95 + DN(0,1)*clhs163 + DN(0,2)*clhs165 + clhs166);
            lhs(1,11)=ave_c1*(DN(2,1)*clhs43 + clhs100*clhs139 - clhs167 + clhs168*clhs73);
            lhs(1,12)=ave_c1*(DN(0,0)*clhs104 + DN(0,1)*clhs169 + DN(0,2)*clhs170 + clhs171);
            lhs(1,13)=ave_c1*(DN(0,0)*clhs116 + DN(0,1)*clhs172 + DN(0,2)*clhs174 + clhs175 + clhs176);
            lhs(1,14)=ave_c1*(DN(0,0)*clhs123 + DN(0,1)*clhs177 + DN(0,2)*clhs179 + clhs180);
            lhs(1,15)=ave_c1*(DN(3,1)*clhs43 + clhs128*clhs139 - clhs181 + clhs182*clhs73);
            lhs(2,0)=ave_c1*(DN(0,0)*clhs4 + DN(0,1)*clhs131 + DN(0,2)*clhs183 + clhs37);
            lhs(2,1)=ave_c1*(DN(0,0)*clhs29 + DN(0,1)*clhs134 + DN(0,2)*clhs184 + clhs140);
            lhs(2,2)=ave_c1*(DN(0,0)*clhs36 + DN(0,1)*clhs138 + DN(0,2)*clhs185 + clhs10*clhs186 + clhs23);
            lhs(2,3)=DN(0,2)*clhs44;
            lhs(2,4)=ave_c1*(DN(0,0)*clhs49 + DN(0,1)*clhs142 + DN(0,2)*clhs187 + clhs189);
            lhs(2,5)=ave_c1*(DN(0,0)*clhs62 + DN(0,1)*clhs146 + DN(0,2)*clhs190 + clhs191);
            lhs(2,6)=ave_c1*(DN(0,0)*clhs68 + DN(0,1)*clhs151 + DN(0,2)*clhs192 + clhs148 + clhs193);
            lhs(2,7)=ave_c1*(DN(1,2)*clhs43 + clhs188*clhs71 - clhs194 + clhs195*clhs73);
            lhs(2,8)=ave_c1*(DN(0,0)*clhs78 + DN(0,1)*clhs156 + DN(0,2)*clhs196 + clhs197);
            lhs(2,9)=ave_c1*(DN(0,0)*clhs91 + DN(0,1)*clhs160 + DN(0,2)*clhs198 + clhs199);
            lhs(2,10)=ave_c1*(DN(0,0)*clhs97 + DN(0,1)*clhs165 + DN(0,2)*clhs200 + clhs162 + clhs201);
            lhs(2,11)=ave_c1*(DN(2,2)*clhs43 + clhs100*clhs188 - clhs202 + clhs203*clhs73);
            lhs(2,12)=ave_c1*(DN(0,0)*clhs106 + DN(0,1)*clhs170 + DN(0,2)*clhs204 + clhs205);
            lhs(2,13)=ave_c1*(DN(0,0)*clhs119 + DN(0,1)*clhs174 + DN(0,2)*clhs206 + clhs207);
            lhs(2,14)=ave_c1*(DN(0,0)*clhs125 + DN(0,1)*clhs179 + DN(0,2)*clhs208 + clhs176 + clhs209);
            lhs(2,15)=ave_c1*(DN(3,2)*clhs43 + clhs128*clhs188 - clhs210 + clhs211*clhs73);
            lhs(3,0)=DN(0,0)*clhs213;
            lhs(3,1)=DN(0,1)*clhs213;
            lhs(3,2)=DN(0,2)*clhs213;
            lhs(3,3)=ave_c1*(clhs13*clhs214 + clhs135*clhs17 + clhs17*clhs186 + clhs17*clhs5);
            lhs(3,4)=ave_c1*(DN(0,0)*clhs215 + clhs72);
            lhs(3,5)=ave_c1*(DN(0,1)*clhs215 + clhs154);
            lhs(3,6)=ave_c1*(DN(0,2)*clhs215 + clhs195);
            lhs(3,7)=clhs219;
            lhs(3,8)=ave_c1*(DN(0,0)*clhs220 + clhs101);
            lhs(3,9)=ave_c1*(DN(0,1)*clhs220 + clhs168);
            lhs(3,10)=ave_c1*(DN(0,2)*clhs220 + clhs203);
            lhs(3,11)=clhs221;
            lhs(3,12)=ave_c1*(DN(0,0)*clhs222 + clhs129);
            lhs(3,13)=ave_c1*(DN(0,1)*clhs222 + clhs182);
            lhs(3,14)=ave_c1*(DN(0,2)*clhs222 + clhs211);
            lhs(3,15)=clhs223;
            lhs(4,0)=ave_c1*(DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs226 + clhs52);
            lhs(4,1)=ave_c1*(DN(1,0)*clhs24 + DN(1,1)*clhs26 + DN(1,2)*clhs29 + clhs143);
            lhs(4,2)=ave_c1*(DN(1,0)*clhs32 + DN(1,1)*clhs34 + DN(1,2)*clhs36 + clhs189);
            lhs(4,3)=ave_c1*(DN(0,0)*clhs228 + clhs227*clhs39 + clhs70*clhs73 - clhs72);
            lhs(4,4)=ave_c1*(DN(1,0)*clhs45 + DN(1,1)*clhs47 + DN(1,2)*clhs49 + clhs10*clhs229 + clhs231);
            lhs(4,5)=ave_c1*(DN(1,0)*clhs57 + DN(1,1)*clhs59 + DN(1,2)*clhs62 + clhs232);
            lhs(4,6)=ave_c1*(DN(1,0)*clhs64 + DN(1,1)*clhs66 + DN(1,2)*clhs68 + clhs233);
            lhs(4,7)=DN(1,0)*clhs235;
            lhs(4,8)=ave_c1*(DN(1,0)*clhs74 + DN(1,1)*clhs76 + DN(1,2)*clhs78 + clhs237 + clhs238);
            lhs(4,9)=ave_c1*(DN(1,0)*clhs86 + DN(1,1)*clhs88 + DN(1,2)*clhs91 + clhs239);
            lhs(4,10)=ave_c1*(DN(1,0)*clhs93 + DN(1,1)*clhs95 + DN(1,2)*clhs97 + clhs240);
            lhs(4,11)=ave_c1*(DN(2,0)*clhs228 + clhs100*clhs227 - clhs241 + clhs242*clhs73);
            lhs(4,12)=ave_c1*(DN(1,0)*clhs102 + DN(1,1)*clhs104 + DN(1,2)*clhs106 + clhs244 + clhs245);
            lhs(4,13)=ave_c1*(DN(1,0)*clhs114 + DN(1,1)*clhs116 + DN(1,2)*clhs119 + clhs246);
            lhs(4,14)=ave_c1*(DN(1,0)*clhs121 + DN(1,1)*clhs123 + DN(1,2)*clhs125 + clhs247);
            lhs(4,15)=ave_c1*(DN(3,0)*clhs228 + clhs128*clhs227 - clhs248 + clhs249*clhs73);
            lhs(5,0)=ave_c1*(DN(1,0)*clhs2 + DN(1,1)*clhs130 + DN(1,2)*clhs131 + clhs63);
            lhs(5,1)=ave_c1*(DN(1,0)*clhs26 + DN(1,1)*clhs132 + DN(1,2)*clhs134 + clhs147 + clhs250);
            lhs(5,2)=ave_c1*(DN(1,0)*clhs34 + DN(1,1)*clhs136 + DN(1,2)*clhs138 + clhs191);
            lhs(5,3)=ave_c1*(DN(0,1)*clhs228 + clhs153*clhs73 - clhs154 + clhs251*clhs39);
            lhs(5,4)=ave_c1*(DN(1,0)*clhs47 + DN(1,1)*clhs141 + DN(1,2)*clhs142 + clhs232);
            lhs(5,5)=ave_c1*(DN(1,0)*clhs59 + DN(1,1)*clhs144 + DN(1,2)*clhs146 + clhs10*clhs252 + clhs231);
            lhs(5,6)=ave_c1*(DN(1,0)*clhs66 + DN(1,1)*clhs149 + DN(1,2)*clhs151 + clhs253);
            lhs(5,7)=DN(1,1)*clhs235;
            lhs(5,8)=ave_c1*(DN(1,0)*clhs76 + DN(1,1)*clhs155 + DN(1,2)*clhs156 + clhs254);
            lhs(5,9)=ave_c1*(DN(1,0)*clhs88 + DN(1,1)*clhs158 + DN(1,2)*clhs160 + clhs255 + clhs256);
            lhs(5,10)=ave_c1*(DN(1,0)*clhs95 + DN(1,1)*clhs163 + DN(1,2)*clhs165 + clhs257);
            lhs(5,11)=ave_c1*(DN(2,1)*clhs228 + clhs100*clhs251 - clhs258 + clhs259*clhs73);
            lhs(5,12)=ave_c1*(DN(1,0)*clhs104 + DN(1,1)*clhs169 + DN(1,2)*clhs170 + clhs260);
            lhs(5,13)=ave_c1*(DN(1,0)*clhs116 + DN(1,1)*clhs172 + DN(1,2)*clhs174 + clhs261 + clhs262);
            lhs(5,14)=ave_c1*(DN(1,0)*clhs123 + DN(1,1)*clhs177 + DN(1,2)*clhs179 + clhs263);
            lhs(5,15)=ave_c1*(DN(3,1)*clhs228 + clhs128*clhs251 - clhs264 + clhs265*clhs73);
            lhs(6,0)=ave_c1*(DN(1,0)*clhs4 + DN(1,1)*clhs131 + DN(1,2)*clhs183 + clhs69);
            lhs(6,1)=ave_c1*(DN(1,0)*clhs29 + DN(1,1)*clhs134 + DN(1,2)*clhs184 + clhs152);
            lhs(6,2)=ave_c1*(DN(1,0)*clhs36 + DN(1,1)*clhs138 + DN(1,2)*clhs185 + clhs193 + clhs250);
            lhs(6,3)=ave_c1*(DN(0,2)*clhs228 + clhs194*clhs73 - clhs195 + clhs266*clhs39);
            lhs(6,4)=ave_c1*(DN(1,0)*clhs49 + DN(1,1)*clhs142 + DN(1,2)*clhs187 + clhs233);
            lhs(6,5)=ave_c1*(DN(1,0)*clhs62 + DN(1,1)*clhs146 + DN(1,2)*clhs190 + clhs253);
            lhs(6,6)=ave_c1*(DN(1,0)*clhs68 + DN(1,1)*clhs151 + DN(1,2)*clhs192 + clhs10*clhs267 + clhs231);
            lhs(6,7)=DN(1,2)*clhs235;
            lhs(6,8)=ave_c1*(DN(1,0)*clhs78 + DN(1,1)*clhs156 + DN(1,2)*clhs196 + clhs268);
            lhs(6,9)=ave_c1*(DN(1,0)*clhs91 + DN(1,1)*clhs160 + DN(1,2)*clhs198 + clhs269);
            lhs(6,10)=ave_c1*(DN(1,0)*clhs97 + DN(1,1)*clhs165 + DN(1,2)*clhs200 + clhs256 + clhs270);
            lhs(6,11)=ave_c1*(DN(2,2)*clhs228 + clhs100*clhs266 - clhs271 + clhs272*clhs73);
            lhs(6,12)=ave_c1*(DN(1,0)*clhs106 + DN(1,1)*clhs170 + DN(1,2)*clhs204 + clhs273);
            lhs(6,13)=ave_c1*(DN(1,0)*clhs119 + DN(1,1)*clhs174 + DN(1,2)*clhs206 + clhs274);
            lhs(6,14)=ave_c1*(DN(1,0)*clhs125 + DN(1,1)*clhs179 + DN(1,2)*clhs208 + clhs262 + clhs275);
            lhs(6,15)=ave_c1*(DN(3,2)*clhs228 + clhs128*clhs266 - clhs276 + clhs277*clhs73);
            lhs(7,0)=ave_c1*(DN(1,0)*clhs212 + clhs70);
            lhs(7,1)=ave_c1*(DN(1,1)*clhs212 + clhs153);
            lhs(7,2)=ave_c1*(DN(1,2)*clhs212 + clhs194);
            lhs(7,3)=clhs219;
            lhs(7,4)=DN(1,0)*clhs278;
            lhs(7,5)=DN(1,1)*clhs278;
            lhs(7,6)=DN(1,2)*clhs278;
            lhs(7,7)=ave_c1*(clhs17*clhs229 + clhs17*clhs252 + clhs17*clhs267 + clhs214*clhs230);
            lhs(7,8)=ave_c1*(DN(1,0)*clhs220 + clhs242);
            lhs(7,9)=ave_c1*(DN(1,1)*clhs220 + clhs259);
            lhs(7,10)=ave_c1*(DN(1,2)*clhs220 + clhs272);
            lhs(7,11)=clhs282;
            lhs(7,12)=ave_c1*(DN(1,0)*clhs222 + clhs249);
            lhs(7,13)=ave_c1*(DN(1,1)*clhs222 + clhs265);
            lhs(7,14)=ave_c1*(DN(1,2)*clhs222 + clhs277);
            lhs(7,15)=clhs283;
            lhs(8,0)=ave_c1*(DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs286 + clhs81);
            lhs(8,1)=ave_c1*(DN(2,0)*clhs24 + DN(2,1)*clhs26 + DN(2,2)*clhs29 + clhs157);
            lhs(8,2)=ave_c1*(DN(2,0)*clhs32 + DN(2,1)*clhs34 + DN(2,2)*clhs36 + clhs197);
            lhs(8,3)=ave_c1*(DN(0,0)*clhs288 - clhs101 + clhs287*clhs39 + clhs73*clhs99);
            lhs(8,4)=ave_c1*(DN(2,0)*clhs45 + DN(2,1)*clhs47 + DN(2,2)*clhs49 + clhs237 + clhs289);
            lhs(8,5)=ave_c1*(DN(2,0)*clhs57 + DN(2,1)*clhs59 + DN(2,2)*clhs62 + clhs254);
            lhs(8,6)=ave_c1*(DN(2,0)*clhs64 + DN(2,1)*clhs66 + DN(2,2)*clhs68 + clhs268);
            lhs(8,7)=ave_c1*(DN(1,0)*clhs288 + clhs241*clhs73 - clhs242 + clhs287*clhs71);
            lhs(8,8)=ave_c1*(DN(2,0)*clhs74 + DN(2,1)*clhs76 + DN(2,2)*clhs78 + clhs10*clhs290 + clhs292);
            lhs(8,9)=ave_c1*(DN(2,0)*clhs86 + DN(2,1)*clhs88 + DN(2,2)*clhs91 + clhs293);
            lhs(8,10)=ave_c1*(DN(2,0)*clhs93 + DN(2,1)*clhs95 + DN(2,2)*clhs97 + clhs294);
            lhs(8,11)=DN(2,0)*clhs296;
            lhs(8,12)=ave_c1*(DN(2,0)*clhs102 + DN(2,1)*clhs104 + DN(2,2)*clhs106 + clhs298 + clhs299);
            lhs(8,13)=ave_c1*(DN(2,0)*clhs114 + DN(2,1)*clhs116 + DN(2,2)*clhs119 + clhs300);
            lhs(8,14)=ave_c1*(DN(2,0)*clhs121 + DN(2,1)*clhs123 + DN(2,2)*clhs125 + clhs301);
            lhs(8,15)=ave_c1*(DN(3,0)*clhs288 + clhs128*clhs287 - clhs302 + clhs303*clhs73);
            lhs(9,0)=ave_c1*(DN(2,0)*clhs2 + DN(2,1)*clhs130 + DN(2,2)*clhs131 + clhs92);
            lhs(9,1)=ave_c1*(DN(2,0)*clhs26 + DN(2,1)*clhs132 + DN(2,2)*clhs134 + clhs161 + clhs304);
            lhs(9,2)=ave_c1*(DN(2,0)*clhs34 + DN(2,1)*clhs136 + DN(2,2)*clhs138 + clhs199);
            lhs(9,3)=ave_c1*(DN(0,1)*clhs288 + clhs167*clhs73 - clhs168 + clhs305*clhs39);
            lhs(9,4)=ave_c1*(DN(2,0)*clhs47 + DN(2,1)*clhs141 + DN(2,2)*clhs142 + clhs239);
            lhs(9,5)=ave_c1*(DN(2,0)*clhs59 + DN(2,1)*clhs144 + DN(2,2)*clhs146 + clhs255 + clhs306);
            lhs(9,6)=ave_c1*(DN(2,0)*clhs66 + DN(2,1)*clhs149 + DN(2,2)*clhs151 + clhs269);
            lhs(9,7)=ave_c1*(DN(1,1)*clhs288 + clhs258*clhs73 - clhs259 + clhs305*clhs71);
            lhs(9,8)=ave_c1*(DN(2,0)*clhs76 + DN(2,1)*clhs155 + DN(2,2)*clhs156 + clhs293);
            lhs(9,9)=ave_c1*(DN(2,0)*clhs88 + DN(2,1)*clhs158 + DN(2,2)*clhs160 + clhs10*clhs307 + clhs292);
            lhs(9,10)=ave_c1*(DN(2,0)*clhs95 + DN(2,1)*clhs163 + DN(2,2)*clhs165 + clhs308);
            lhs(9,11)=DN(2,1)*clhs296;
            lhs(9,12)=ave_c1*(DN(2,0)*clhs104 + DN(2,1)*clhs169 + DN(2,2)*clhs170 + clhs309);
            lhs(9,13)=ave_c1*(DN(2,0)*clhs116 + DN(2,1)*clhs172 + DN(2,2)*clhs174 + clhs310 + clhs311);
            lhs(9,14)=ave_c1*(DN(2,0)*clhs123 + DN(2,1)*clhs177 + DN(2,2)*clhs179 + clhs312);
            lhs(9,15)=ave_c1*(DN(3,1)*clhs288 + clhs128*clhs305 - clhs313 + clhs314*clhs73);
            lhs(10,0)=ave_c1*(DN(2,0)*clhs4 + DN(2,1)*clhs131 + DN(2,2)*clhs183 + clhs98);
            lhs(10,1)=ave_c1*(DN(2,0)*clhs29 + DN(2,1)*clhs134 + DN(2,2)*clhs184 + clhs166);
            lhs(10,2)=ave_c1*(DN(2,0)*clhs36 + DN(2,1)*clhs138 + DN(2,2)*clhs185 + clhs201 + clhs304);
            lhs(10,3)=ave_c1*(DN(0,2)*clhs288 + clhs202*clhs73 - clhs203 + clhs315*clhs39);
            lhs(10,4)=ave_c1*(DN(2,0)*clhs49 + DN(2,1)*clhs142 + DN(2,2)*clhs187 + clhs240);
            lhs(10,5)=ave_c1*(DN(2,0)*clhs62 + DN(2,1)*clhs146 + DN(2,2)*clhs190 + clhs257);
            lhs(10,6)=ave_c1*(DN(2,0)*clhs68 + DN(2,1)*clhs151 + DN(2,2)*clhs192 + clhs270 + clhs306);
            lhs(10,7)=ave_c1*(DN(1,2)*clhs288 + clhs271*clhs73 - clhs272 + clhs315*clhs71);
            lhs(10,8)=ave_c1*(DN(2,0)*clhs78 + DN(2,1)*clhs156 + DN(2,2)*clhs196 + clhs294);
            lhs(10,9)=ave_c1*(DN(2,0)*clhs91 + DN(2,1)*clhs160 + DN(2,2)*clhs198 + clhs308);
            lhs(10,10)=ave_c1*(DN(2,0)*clhs97 + DN(2,1)*clhs165 + DN(2,2)*clhs200 + clhs10*clhs316 + clhs292);
            lhs(10,11)=DN(2,2)*clhs296;
            lhs(10,12)=ave_c1*(DN(2,0)*clhs106 + DN(2,1)*clhs170 + DN(2,2)*clhs204 + clhs317);
            lhs(10,13)=ave_c1*(DN(2,0)*clhs119 + DN(2,1)*clhs174 + DN(2,2)*clhs206 + clhs318);
            lhs(10,14)=ave_c1*(DN(2,0)*clhs125 + DN(2,1)*clhs179 + DN(2,2)*clhs208 + clhs311 + clhs319);
            lhs(10,15)=ave_c1*(DN(3,2)*clhs288 + clhs128*clhs315 - clhs320 + clhs321*clhs73);
            lhs(11,0)=ave_c1*(DN(2,0)*clhs212 + clhs99);
            lhs(11,1)=ave_c1*(DN(2,1)*clhs212 + clhs167);
            lhs(11,2)=ave_c1*(DN(2,2)*clhs212 + clhs202);
            lhs(11,3)=clhs221;
            lhs(11,4)=ave_c1*(DN(2,0)*clhs215 + clhs241);
            lhs(11,5)=ave_c1*(DN(2,1)*clhs215 + clhs258);
            lhs(11,6)=ave_c1*(DN(2,2)*clhs215 + clhs271);
            lhs(11,7)=clhs282;
            lhs(11,8)=DN(2,0)*clhs322;
            lhs(11,9)=DN(2,1)*clhs322;
            lhs(11,10)=DN(2,2)*clhs322;
            lhs(11,11)=ave_c1*(clhs17*clhs290 + clhs17*clhs307 + clhs17*clhs316 + clhs214*clhs291);
            lhs(11,12)=ave_c1*(DN(2,0)*clhs222 + clhs303);
            lhs(11,13)=ave_c1*(DN(2,1)*clhs222 + clhs314);
            lhs(11,14)=ave_c1*(DN(2,2)*clhs222 + clhs321);
            lhs(11,15)=clhs323;
            lhs(12,0)=ave_c1*(DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs109 + clhs326);
            lhs(12,1)=ave_c1*(DN(3,0)*clhs24 + DN(3,1)*clhs26 + DN(3,2)*clhs29 + clhs171);
            lhs(12,2)=ave_c1*(DN(3,0)*clhs32 + DN(3,1)*clhs34 + DN(3,2)*clhs36 + clhs205);
            lhs(12,3)=ave_c1*(DN(0,0)*clhs328 + clhs127*clhs73 - clhs129 + clhs327*clhs39);
            lhs(12,4)=ave_c1*(DN(3,0)*clhs45 + DN(3,1)*clhs47 + DN(3,2)*clhs49 + clhs244 + clhs329);
            lhs(12,5)=ave_c1*(DN(3,0)*clhs57 + DN(3,1)*clhs59 + DN(3,2)*clhs62 + clhs260);
            lhs(12,6)=ave_c1*(DN(3,0)*clhs64 + DN(3,1)*clhs66 + DN(3,2)*clhs68 + clhs273);
            lhs(12,7)=ave_c1*(DN(1,0)*clhs328 + clhs248*clhs73 - clhs249 + clhs327*clhs71);
            lhs(12,8)=ave_c1*(DN(3,0)*clhs74 + DN(3,1)*clhs76 + DN(3,2)*clhs78 + clhs298 + clhs330);
            lhs(12,9)=ave_c1*(DN(3,0)*clhs86 + DN(3,1)*clhs88 + DN(3,2)*clhs91 + clhs309);
            lhs(12,10)=ave_c1*(DN(3,0)*clhs93 + DN(3,1)*clhs95 + DN(3,2)*clhs97 + clhs317);
            lhs(12,11)=ave_c1*(DN(2,0)*clhs328 + clhs100*clhs327 + clhs302*clhs73 - clhs303);
            lhs(12,12)=ave_c1*(DN(3,0)*clhs102 + DN(3,1)*clhs104 + DN(3,2)*clhs106 + clhs10*clhs331 + clhs333);
            lhs(12,13)=ave_c1*(DN(3,0)*clhs114 + DN(3,1)*clhs116 + DN(3,2)*clhs119 + clhs334);
            lhs(12,14)=ave_c1*(DN(3,0)*clhs121 + DN(3,1)*clhs123 + DN(3,2)*clhs125 + clhs335);
            lhs(12,15)=DN(3,0)*clhs336;
            lhs(13,0)=ave_c1*(DN(3,0)*clhs2 + DN(3,1)*clhs130 + DN(3,2)*clhs131 + clhs120);
            lhs(13,1)=ave_c1*(DN(3,0)*clhs26 + DN(3,1)*clhs132 + DN(3,2)*clhs134 + clhs175 + clhs337);
            lhs(13,2)=ave_c1*(DN(3,0)*clhs34 + DN(3,1)*clhs136 + DN(3,2)*clhs138 + clhs207);
            lhs(13,3)=ave_c1*(DN(0,1)*clhs328 + clhs181*clhs73 - clhs182 + clhs338*clhs39);
            lhs(13,4)=ave_c1*(DN(3,0)*clhs47 + DN(3,1)*clhs141 + DN(3,2)*clhs142 + clhs246);
            lhs(13,5)=ave_c1*(DN(3,0)*clhs59 + DN(3,1)*clhs144 + DN(3,2)*clhs146 + clhs261 + clhs339);
            lhs(13,6)=ave_c1*(DN(3,0)*clhs66 + DN(3,1)*clhs149 + DN(3,2)*clhs151 + clhs274);
            lhs(13,7)=ave_c1*(DN(1,1)*clhs328 + clhs264*clhs73 - clhs265 + clhs338*clhs71);
            lhs(13,8)=ave_c1*(DN(3,0)*clhs76 + DN(3,1)*clhs155 + DN(3,2)*clhs156 + clhs300);
            lhs(13,9)=ave_c1*(DN(3,0)*clhs88 + DN(3,1)*clhs158 + DN(3,2)*clhs160 + clhs310 + clhs340);
            lhs(13,10)=ave_c1*(DN(3,0)*clhs95 + DN(3,1)*clhs163 + DN(3,2)*clhs165 + clhs318);
            lhs(13,11)=ave_c1*(DN(2,1)*clhs328 + clhs100*clhs338 + clhs313*clhs73 - clhs314);
            lhs(13,12)=ave_c1*(DN(3,0)*clhs104 + DN(3,1)*clhs169 + DN(3,2)*clhs170 + clhs334);
            lhs(13,13)=ave_c1*(DN(3,0)*clhs116 + DN(3,1)*clhs172 + DN(3,2)*clhs174 + clhs10*clhs341 + clhs333);
            lhs(13,14)=ave_c1*(DN(3,0)*clhs123 + DN(3,1)*clhs177 + DN(3,2)*clhs179 + clhs342);
            lhs(13,15)=DN(3,1)*clhs336;
            lhs(14,0)=ave_c1*(DN(3,0)*clhs4 + DN(3,1)*clhs131 + DN(3,2)*clhs183 + clhs126);
            lhs(14,1)=ave_c1*(DN(3,0)*clhs29 + DN(3,1)*clhs134 + DN(3,2)*clhs184 + clhs180);
            lhs(14,2)=ave_c1*(DN(3,0)*clhs36 + DN(3,1)*clhs138 + DN(3,2)*clhs185 + clhs209 + clhs337);
            lhs(14,3)=ave_c1*(DN(0,2)*clhs328 + DN(3,2)*clhs40 + clhs210*clhs73 - clhs211);
            lhs(14,4)=ave_c1*(DN(3,0)*clhs49 + DN(3,1)*clhs142 + DN(3,2)*clhs187 + clhs247);
            lhs(14,5)=ave_c1*(DN(3,0)*clhs62 + DN(3,1)*clhs146 + DN(3,2)*clhs190 + clhs263);
            lhs(14,6)=ave_c1*(DN(3,0)*clhs68 + DN(3,1)*clhs151 + DN(3,2)*clhs192 + clhs275 + clhs339);
            lhs(14,7)=ave_c1*(DN(1,2)*clhs328 + DN(3,2)*clhs234 + clhs276*clhs73 - clhs277);
            lhs(14,8)=ave_c1*(DN(3,0)*clhs78 + DN(3,1)*clhs156 + DN(3,2)*clhs196 + clhs301);
            lhs(14,9)=ave_c1*(DN(3,0)*clhs91 + DN(3,1)*clhs160 + DN(3,2)*clhs198 + clhs312);
            lhs(14,10)=ave_c1*(DN(3,0)*clhs97 + DN(3,1)*clhs165 + DN(3,2)*clhs200 + clhs319 + clhs340);
            lhs(14,11)=ave_c1*(DN(2,2)*clhs328 + DN(3,2)*clhs295 + clhs320*clhs73 - clhs321);
            lhs(14,12)=ave_c1*(DN(3,0)*clhs106 + DN(3,1)*clhs170 + DN(3,2)*clhs204 + clhs335);
            lhs(14,13)=ave_c1*(DN(3,0)*clhs119 + DN(3,1)*clhs174 + DN(3,2)*clhs206 + clhs342);
            lhs(14,14)=ave_c1*(DN(3,0)*clhs125 + DN(3,1)*clhs179 + DN(3,2)*clhs208 + clhs10*clhs343 + clhs333);
            lhs(14,15)=DN(3,2)*clhs336;
            lhs(15,0)=ave_c1*(DN(3,0)*clhs212 + clhs127);
            lhs(15,1)=ave_c1*(DN(3,1)*clhs212 + clhs181);
            lhs(15,2)=ave_c1*(DN(3,2)*clhs212 + clhs210);
            lhs(15,3)=clhs223;
            lhs(15,4)=ave_c1*(DN(3,0)*clhs215 + clhs248);
            lhs(15,5)=ave_c1*(DN(3,1)*clhs215 + clhs264);
            lhs(15,6)=ave_c1*(DN(3,2)*clhs215 + clhs276);
            lhs(15,7)=clhs283;
            lhs(15,8)=ave_c1*(DN(3,0)*clhs220 + clhs302);
            lhs(15,9)=ave_c1*(DN(3,1)*clhs220 + clhs313);
            lhs(15,10)=ave_c1*(DN(3,2)*clhs220 + clhs320);
            lhs(15,11)=clhs323;
            lhs(15,12)=DN(3,0)*clhs344;
            lhs(15,13)=DN(3,1)*clhs344;
            lhs(15,14)=DN(3,2)*clhs344;
            lhs(15,15)=ave_c1*(clhs17*clhs331 + clhs17*clhs341 + clhs17*clhs343 + clhs214*clhs332);

}


template<>
void TimeAveragedNavierStokes<2>::ComputeGaussPointLHSContribution(
    BoundedMatrix<double,9,9>& lhs,
    const ElementDataStruct& data)
{
    constexpr int dim = 2;
    constexpr int nnodes = 3;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size
    const double c = data.c;                                // Wave velocity

    //const double& dts = data.dts;                           // The averaging time period
    const double& dtn = data.dtn;                           // Time increment: notice t = tn + dtn

    const double& bdf0 = data.bdf0;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v_ave = data.v_ave;
    const BoundedMatrix<double,nnodes,dim>& vn_ave = data.vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;

    // time averaging parameters
    double ave_c1 = data.ave_c1;
    double ave_c2 = data.ave_c2;

    // get time accurate expression
    const BoundedMatrix<double,nnodes,dim>& v = ave_c1 * v_ave - ave_c2 * vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vconv = v - vmesh;

    // get constitutive matrix
    const Matrix& C = data.C;

    // get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    const double clhs0 =             C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double clhs1 =             C(0,2)*DN(0,0);
const double clhs2 =             C(2,2)*DN(0,1) + clhs1;
const double clhs3 =             pow(DN(0,0), 2);
const double clhs4 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double clhs5 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double clhs6 =             rho*stab_c2*sqrt(pow(clhs4, 2) + pow(clhs5, 2));
const double clhs7 =             clhs6*h/stab_c1 + mu;
const double clhs8 =             DN(0,0)*clhs4 + DN(0,1)*clhs5;
const double clhs9 =             N[0]*rho;
const double clhs10 =             pow(N[0], 2);
const double clhs11 =             bdf0*rho;
const double clhs12 =             N[0]*bdf0;
const double clhs13 =             clhs12 + clhs8;
const double clhs14 =             1.0*h/(clhs6 + mu*stab_c1/h);
const double clhs15 =             clhs14*pow(rho, 2);
const double clhs16 =             clhs15*clhs8;
const double clhs17 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double clhs18 =             clhs15*clhs17;
const double clhs19 =             N[0]*clhs18;
const double clhs20 =             clhs10*clhs11 + clhs13*clhs16 + clhs13*clhs19 + clhs8*clhs9;
const double clhs21 =             C(0,1)*DN(0,1) + clhs1;
const double clhs22 =             C(1,2)*DN(0,1);
const double clhs23 =             C(2,2)*DN(0,0) + clhs22;
const double clhs24 =             DN(0,0)*clhs7;
const double clhs25 =             DN(0,1)*clhs24;
const double clhs26 =             1/(pow(c, 2)*rho);
const double clhs27 =             clhs12*clhs26;
const double clhs28 =             clhs27*clhs7;
const double clhs29 =             clhs14*clhs17;
const double clhs30 =             clhs14*rho;
const double clhs31 =             clhs30*clhs8;
const double clhs32 =             ave_c1*(-N[0] + clhs28 + clhs29*clhs9 + clhs31);
const double clhs33 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs34 =             C(0,2)*DN(1,0);
const double clhs35 =             C(2,2)*DN(1,1) + clhs34;
const double clhs36 =             N[1]*rho;
const double clhs37 =             clhs12*clhs36;
const double clhs38 =             DN(1,0)*clhs24 + clhs37;
const double clhs39 =             DN(1,0)*clhs4 + DN(1,1)*clhs5;
const double clhs40 =             N[1]*bdf0;
const double clhs41 =             clhs39 + clhs40;
const double clhs42 =             clhs16*clhs41 + clhs19*clhs41 + clhs39*clhs9;
const double clhs43 =             C(0,1)*DN(1,1) + clhs34;
const double clhs44 =             C(1,2)*DN(1,1);
const double clhs45 =             C(2,2)*DN(1,0) + clhs44;
const double clhs46 =             DN(1,1)*clhs24;
const double clhs47 =             DN(0,0)*N[1];
const double clhs48 =             clhs26*clhs40;
const double clhs49 =             DN(1,0)*N[0];
const double clhs50 =             clhs17*clhs30;
const double clhs51 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs52 =             C(0,2)*DN(2,0);
const double clhs53 =             C(2,2)*DN(2,1) + clhs52;
const double clhs54 =             N[2]*rho;
const double clhs55 =             clhs12*clhs54;
const double clhs56 =             DN(2,0)*clhs24 + clhs55;
const double clhs57 =             DN(2,0)*clhs4 + DN(2,1)*clhs5;
const double clhs58 =             N[2]*bdf0;
const double clhs59 =             clhs57 + clhs58;
const double clhs60 =             clhs16*clhs59 + clhs19*clhs59 + clhs57*clhs9;
const double clhs61 =             C(0,1)*DN(2,1) + clhs52;
const double clhs62 =             C(1,2)*DN(2,1);
const double clhs63 =             C(2,2)*DN(2,0) + clhs62;
const double clhs64 =             DN(2,1)*clhs24;
const double clhs65 =             DN(0,0)*N[2];
const double clhs66 =             clhs26*clhs58;
const double clhs67 =             DN(2,0)*N[0];
const double clhs68 =             C(0,1)*DN(0,0) + clhs22;
const double clhs69 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs70 =             pow(DN(0,1), 2);
const double clhs71 =             C(0,1)*DN(1,0) + clhs44;
const double clhs72 =             DN(0,1)*clhs7;
const double clhs73 =             DN(1,0)*clhs72;
const double clhs74 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs75 =             DN(1,1)*clhs72 + clhs37;
const double clhs76 =             DN(0,1)*N[1];
const double clhs77 =             DN(1,1)*N[0];
const double clhs78 =             C(0,1)*DN(2,0) + clhs62;
const double clhs79 =             DN(2,0)*clhs72;
const double clhs80 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs81 =             DN(2,1)*clhs72 + clhs55;
const double clhs82 =             DN(0,1)*N[2];
const double clhs83 =             DN(2,1)*N[0];
const double clhs84 =             clhs13*clhs30;
const double clhs85 =             ave_c1*(N[0] + clhs84);
const double clhs86 =             bdf0*clhs26;
const double clhs87 =             clhs30*clhs41;
const double clhs88 =             DN(0,0)*clhs14;
const double clhs89 =             DN(0,1)*clhs14;
const double clhs90 =             ave_c1*(DN(1,0)*clhs88 + DN(1,1)*clhs89 + N[1]*clhs27);
const double clhs91 =             clhs30*clhs59;
const double clhs92 =             ave_c1*(DN(2,0)*clhs88 + DN(2,1)*clhs89 + N[2]*clhs27);
const double clhs93 =             clhs15*clhs39;
const double clhs94 =             N[1]*clhs18;
const double clhs95 =             clhs13*clhs93 + clhs13*clhs94 + clhs36*clhs8;
const double clhs96 =             DN(1,0)*clhs7;
const double clhs97 =             clhs30*clhs39;
const double clhs98 =             pow(DN(1,0), 2);
const double clhs99 =             pow(N[1], 2);
const double clhs100 =             clhs11*clhs99 + clhs36*clhs39 + clhs41*clhs93 + clhs41*clhs94;
const double clhs101 =             DN(1,1)*clhs96;
const double clhs102 =             clhs48*clhs7;
const double clhs103 =             ave_c1*(-N[1] + clhs102 + clhs29*clhs36 + clhs97);
const double clhs104 =             clhs40*clhs54;
const double clhs105 =             DN(2,0)*clhs96 + clhs104;
const double clhs106 =             clhs36*clhs57 + clhs59*clhs93 + clhs59*clhs94;
const double clhs107 =             DN(2,1)*clhs96;
const double clhs108 =             DN(1,0)*N[2];
const double clhs109 =             DN(2,0)*N[1];
const double clhs110 =             DN(1,1)*clhs7;
const double clhs111 =             pow(DN(1,1), 2);
const double clhs112 =             DN(2,0)*clhs110;
const double clhs113 =             DN(2,1)*clhs110 + clhs104;
const double clhs114 =             DN(1,1)*N[2];
const double clhs115 =             DN(2,1)*N[1];
const double clhs116 =             ave_c1*(N[1] + clhs87);
const double clhs117 =             ave_c1*(DN(1,0)*DN(2,0)*clhs14 + DN(1,1)*DN(2,1)*clhs14 + N[2]*clhs48);
const double clhs118 =             clhs15*clhs57;
const double clhs119 =             N[2]*clhs18;
const double clhs120 =             clhs118*clhs13 + clhs119*clhs13 + clhs54*clhs8;
const double clhs121 =             DN(2,0)*clhs7;
const double clhs122 =             clhs30*clhs57;
const double clhs123 =             clhs118*clhs41 + clhs119*clhs41 + clhs39*clhs54;
const double clhs124 =             pow(DN(2,0), 2);
const double clhs125 =             pow(N[2], 2);
const double clhs126 =             clhs11*clhs125 + clhs118*clhs59 + clhs119*clhs59 + clhs54*clhs57;
const double clhs127 =             DN(2,1)*clhs121;
const double clhs128 =             ave_c1*(-N[2] + clhs122 + clhs29*clhs54 + clhs66*clhs7);
const double clhs129 =             pow(DN(2,1), 2);
const double clhs130 =             ave_c1*(N[2] + clhs91);
            lhs(0,0)=ave_c1*(DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs20 + clhs3*clhs7);
            lhs(0,1)=ave_c1*(DN(0,0)*clhs21 + DN(0,1)*clhs23 + clhs25);
            lhs(0,2)=DN(0,0)*clhs32;
            lhs(0,3)=ave_c1*(DN(0,0)*clhs33 + DN(0,1)*clhs35 + clhs38 + clhs42);
            lhs(0,4)=ave_c1*(DN(0,0)*clhs43 + DN(0,1)*clhs45 + clhs46);
            lhs(0,5)=ave_c1*(DN(1,0)*clhs31 + clhs24*clhs48 - clhs47 + clhs49*clhs50);
            lhs(0,6)=ave_c1*(DN(0,0)*clhs51 + DN(0,1)*clhs53 + clhs56 + clhs60);
            lhs(0,7)=ave_c1*(DN(0,0)*clhs61 + DN(0,1)*clhs63 + clhs64);
            lhs(0,8)=ave_c1*(DN(2,0)*clhs31 + clhs24*clhs66 + clhs50*clhs67 - clhs65);
            lhs(1,0)=ave_c1*(DN(0,0)*clhs2 + DN(0,1)*clhs68 + clhs25);
            lhs(1,1)=ave_c1*(DN(0,0)*clhs23 + DN(0,1)*clhs69 + clhs20 + clhs7*clhs70);
            lhs(1,2)=DN(0,1)*clhs32;
            lhs(1,3)=ave_c1*(DN(0,0)*clhs35 + DN(0,1)*clhs71 + clhs73);
            lhs(1,4)=ave_c1*(DN(0,0)*clhs45 + DN(0,1)*clhs74 + clhs42 + clhs75);
            lhs(1,5)=ave_c1*(DN(1,1)*clhs31 + clhs48*clhs72 + clhs50*clhs77 - clhs76);
            lhs(1,6)=ave_c1*(DN(0,0)*clhs53 + DN(0,1)*clhs78 + clhs79);
            lhs(1,7)=ave_c1*(DN(0,0)*clhs63 + DN(0,1)*clhs80 + clhs60 + clhs81);
            lhs(1,8)=ave_c1*(DN(2,1)*clhs31 + clhs50*clhs83 + clhs66*clhs72 - clhs82);
            lhs(2,0)=DN(0,0)*clhs85;
            lhs(2,1)=DN(0,1)*clhs85;
            lhs(2,2)=ave_c1*(clhs10*clhs86 + clhs14*clhs3 + clhs14*clhs70);
            lhs(2,3)=ave_c1*(DN(0,0)*clhs87 + clhs49);
            lhs(2,4)=ave_c1*(DN(0,1)*clhs87 + clhs77);
            lhs(2,5)=clhs90;
            lhs(2,6)=ave_c1*(DN(0,0)*clhs91 + clhs67);
            lhs(2,7)=ave_c1*(DN(0,1)*clhs91 + clhs83);
            lhs(2,8)=clhs92;
            lhs(3,0)=ave_c1*(DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs38 + clhs95);
            lhs(3,1)=ave_c1*(DN(1,0)*clhs21 + DN(1,1)*clhs23 + clhs73);
            lhs(3,2)=ave_c1*(DN(0,0)*clhs97 + clhs27*clhs96 + clhs47*clhs50 - clhs49);
            lhs(3,3)=ave_c1*(DN(1,0)*clhs33 + DN(1,1)*clhs35 + clhs100 + clhs7*clhs98);
            lhs(3,4)=ave_c1*(DN(1,0)*clhs43 + DN(1,1)*clhs45 + clhs101);
            lhs(3,5)=DN(1,0)*clhs103;
            lhs(3,6)=ave_c1*(DN(1,0)*clhs51 + DN(1,1)*clhs53 + clhs105 + clhs106);
            lhs(3,7)=ave_c1*(DN(1,0)*clhs61 + DN(1,1)*clhs63 + clhs107);
            lhs(3,8)=ave_c1*(DN(2,0)*clhs97 - clhs108 + clhs109*clhs50 + clhs66*clhs96);
            lhs(4,0)=ave_c1*(DN(1,0)*clhs2 + DN(1,1)*clhs68 + clhs46);
            lhs(4,1)=ave_c1*(DN(1,0)*clhs23 + DN(1,1)*clhs69 + clhs75 + clhs95);
            lhs(4,2)=ave_c1*(DN(0,1)*clhs97 + clhs110*clhs27 + clhs50*clhs76 - clhs77);
            lhs(4,3)=ave_c1*(DN(1,0)*clhs35 + DN(1,1)*clhs71 + clhs101);
            lhs(4,4)=ave_c1*(DN(1,0)*clhs45 + DN(1,1)*clhs74 + clhs100 + clhs111*clhs7);
            lhs(4,5)=DN(1,1)*clhs103;
            lhs(4,6)=ave_c1*(DN(1,0)*clhs53 + DN(1,1)*clhs78 + clhs112);
            lhs(4,7)=ave_c1*(DN(1,0)*clhs63 + DN(1,1)*clhs80 + clhs106 + clhs113);
            lhs(4,8)=ave_c1*(DN(2,1)*clhs97 + clhs110*clhs66 - clhs114 + clhs115*clhs50);
            lhs(5,0)=ave_c1*(DN(1,0)*clhs84 + clhs47);
            lhs(5,1)=ave_c1*(DN(1,1)*clhs84 + clhs76);
            lhs(5,2)=clhs90;
            lhs(5,3)=DN(1,0)*clhs116;
            lhs(5,4)=DN(1,1)*clhs116;
            lhs(5,5)=ave_c1*(clhs111*clhs14 + clhs14*clhs98 + clhs86*clhs99);
            lhs(5,6)=ave_c1*(DN(1,0)*clhs91 + clhs109);
            lhs(5,7)=ave_c1*(DN(1,1)*clhs91 + clhs115);
            lhs(5,8)=clhs117;
            lhs(6,0)=ave_c1*(DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs120 + clhs56);
            lhs(6,1)=ave_c1*(DN(2,0)*clhs21 + DN(2,1)*clhs23 + clhs79);
            lhs(6,2)=ave_c1*(DN(0,0)*clhs122 + clhs121*clhs27 + clhs50*clhs65 - clhs67);
            lhs(6,3)=ave_c1*(DN(2,0)*clhs33 + DN(2,1)*clhs35 + clhs105 + clhs123);
            lhs(6,4)=ave_c1*(DN(2,0)*clhs43 + DN(2,1)*clhs45 + clhs112);
            lhs(6,5)=ave_c1*(DN(1,0)*clhs122 + clhs108*clhs50 - clhs109 + clhs121*clhs48);
            lhs(6,6)=ave_c1*(DN(2,0)*clhs51 + DN(2,1)*clhs53 + clhs124*clhs7 + clhs126);
            lhs(6,7)=ave_c1*(DN(2,0)*clhs61 + DN(2,1)*clhs63 + clhs127);
            lhs(6,8)=DN(2,0)*clhs128;
            lhs(7,0)=ave_c1*(DN(2,0)*clhs2 + DN(2,1)*clhs68 + clhs64);
            lhs(7,1)=ave_c1*(DN(2,0)*clhs23 + DN(2,1)*clhs69 + clhs120 + clhs81);
            lhs(7,2)=ave_c1*(DN(0,1)*clhs122 + DN(2,1)*clhs28 + clhs50*clhs82 - clhs83);
            lhs(7,3)=ave_c1*(DN(2,0)*clhs35 + DN(2,1)*clhs71 + clhs107);
            lhs(7,4)=ave_c1*(DN(2,0)*clhs45 + DN(2,1)*clhs74 + clhs113 + clhs123);
            lhs(7,5)=ave_c1*(DN(1,1)*clhs122 + DN(2,1)*clhs102 + clhs114*clhs50 - clhs115);
            lhs(7,6)=ave_c1*(DN(2,0)*clhs53 + DN(2,1)*clhs78 + clhs127);
            lhs(7,7)=ave_c1*(DN(2,0)*clhs63 + DN(2,1)*clhs80 + clhs126 + clhs129*clhs7);
            lhs(7,8)=DN(2,1)*clhs128;
            lhs(8,0)=ave_c1*(DN(2,0)*clhs84 + clhs65);
            lhs(8,1)=ave_c1*(DN(2,1)*clhs84 + clhs82);
            lhs(8,2)=clhs92;
            lhs(8,3)=ave_c1*(DN(2,0)*clhs87 + clhs108);
            lhs(8,4)=ave_c1*(DN(2,1)*clhs87 + clhs114);
            lhs(8,5)=clhs117;
            lhs(8,6)=DN(2,0)*clhs130;
            lhs(8,7)=DN(2,1)*clhs130;
            lhs(8,8)=ave_c1*(clhs124*clhs14 + clhs125*clhs86 + clhs129*clhs14);
            for (int i=0; i<9; i++){    
                for (int j=0; j<9; j++){    
                    if(isnan(lhs(i,j)))
                    {
                        std::cout<<"##### NaN FOUND ("<<this->GetId()<<") ####"<<std::endl;
                    }
                }
            }
}


template<>
void TimeAveragedNavierStokes<3>::ComputeGaussPointRHSContribution(
    array_1d<double,16>& rhs,
    const ElementDataStruct& data)
{
    constexpr int dim = 3;
    constexpr int nnodes = 4;
    constexpr int strain_size = 6;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size
    const double c = data.c;                                // Wave velocity

    //const double& dts = data.dts;                           // The averaging time period
    const double& dtn = data.dtn;                           // Time increment: notice t = tn + dtn

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v_ave = data.v_ave;
    const BoundedMatrix<double,nnodes,dim>& vn_ave = data.vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnn_ave = data.vnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnnn_ave = data.vnnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& p_ave = data.p_ave;
    const array_1d<double,nnodes>& pn_ave = data.pn_ave;
    const array_1d<double,nnodes>& pnn_ave = data.pnn_ave;
    const array_1d<double,nnodes>& pnnn_ave = data.pnnn_ave;
    const array_1d<double,strain_size>& stress = data.stress;

    // time averaging parameters
    double ave_c1 = data.ave_c1;
    double ave_c2 = data.ave_c2;
    double ave_n_c1 = data.ave_n_c1;
    double ave_n_c2 = data.ave_n_c2;
    double ave_nn_c1 = data.ave_nn_c1;
    double ave_nn_c2 = data.ave_nn_c2;

    // get time accurate expression
    const array_1d<double,nnodes>& p = ave_c1 * p_ave - ave_c2 * pn_ave;
    const BoundedMatrix<double,nnodes,dim>& v = ave_c1 * v_ave - ave_c2 * vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vconv = v - vmesh;

    // get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // auxiliary variables used in the calculation of the RHS
    const array_1d<double,dim> f_gauss = prod(trans(f), N);
    const array_1d<double,dim> grad_p = prod(trans(DN), p);

    // stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    const double crhs0 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crhs1 =             ave_c1*p_ave[0] - ave_c2*pn_ave[0];
const double crhs2 =             ave_c1*p_ave[1] - ave_c2*pn_ave[1];
const double crhs3 =             ave_c1*p_ave[2] - ave_c2*pn_ave[2];
const double crhs4 =             ave_c1*p_ave[3] - ave_c2*pn_ave[3];
const double crhs5 =             N[0]*crhs1 + N[1]*crhs2 + N[2]*crhs3 + N[3]*crhs4;
const double crhs6 =             ave_c1*v_ave(0,0) - ave_c2*vn_ave(0,0);
const double crhs7 =             ave_c1*v_ave(1,0) - ave_c2*vn_ave(1,0);
const double crhs8 =             ave_c1*v_ave(2,0) - ave_c2*vn_ave(2,0);
const double crhs9 =             ave_c1*v_ave(3,0) - ave_c2*vn_ave(3,0);
const double crhs10 =             rho*(N[0]*(bdf0*crhs6 + bdf1*(ave_n_c1*vn_ave(0,0) - ave_n_c2*vnn_ave(0,0)) + bdf2*(ave_nn_c1*vnn_ave(0,0) - ave_nn_c2*vnnn_ave(0,0))) + N[1]*(bdf0*crhs7 + bdf1*(ave_n_c1*vn_ave(1,0) - ave_n_c2*vnn_ave(1,0)) + bdf2*(ave_nn_c1*vnn_ave(1,0) - ave_nn_c2*vnnn_ave(1,0))) + N[2]*(bdf0*crhs8 + bdf1*(ave_n_c1*vn_ave(2,0) - ave_n_c2*vnn_ave(2,0)) + bdf2*(ave_nn_c1*vnn_ave(2,0) - ave_nn_c2*vnnn_ave(2,0))) + N[3]*(bdf0*crhs9 + bdf1*(ave_n_c1*vn_ave(3,0) - ave_n_c2*vnn_ave(3,0)) + bdf2*(ave_nn_c1*vnn_ave(3,0) - ave_nn_c2*vnnn_ave(3,0))));
const double crhs11 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs12 =             DN(0,0)*crhs6 + DN(1,0)*crhs7 + DN(2,0)*crhs8 + DN(3,0)*crhs9;
const double crhs13 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs14 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs15 =             rho*(crhs11*crhs12 + crhs13*(DN(0,1)*crhs6 + DN(1,1)*crhs7 + DN(2,1)*crhs8 + DN(3,1)*crhs9) + crhs14*(DN(0,2)*crhs6 + DN(1,2)*crhs7 + DN(2,2)*crhs8 + DN(3,2)*crhs9));
const double crhs16 =             rho*stab_c2*sqrt(pow(crhs11, 2) + pow(crhs13, 2) + pow(crhs14, 2));
const double crhs17 =             (N[0]*(bdf0*crhs1 + bdf1*(ave_n_c1*pn_ave[0] - ave_n_c2*pnn_ave[0]) + bdf2*(ave_nn_c1*pnn_ave[0] - ave_nn_c2*pnnn_ave[0])) + N[1]*(bdf0*crhs2 + bdf1*(ave_n_c1*pn_ave[1] - ave_n_c2*pnn_ave[1]) + bdf2*(ave_nn_c1*pnn_ave[1] - ave_nn_c2*pnnn_ave[1])) + N[2]*(bdf0*crhs3 + bdf1*(ave_n_c1*pn_ave[2] - ave_n_c2*pnn_ave[2]) + bdf2*(ave_nn_c1*pnn_ave[2] - ave_nn_c2*pnnn_ave[2])) + N[3]*(bdf0*crhs4 + bdf1*(ave_n_c1*pn_ave[3] - ave_n_c2*pnn_ave[3]) + bdf2*(ave_nn_c1*pnn_ave[3] - ave_nn_c2*pnnn_ave[3])))/(pow(c, 2)*rho);
const double crhs18 =             ave_c1*v_ave(0,1) - ave_c2*vn_ave(0,1);
const double crhs19 =             ave_c1*v_ave(1,1) - ave_c2*vn_ave(1,1);
const double crhs20 =             ave_c1*v_ave(2,1) - ave_c2*vn_ave(2,1);
const double crhs21 =             ave_c1*v_ave(3,1) - ave_c2*vn_ave(3,1);
const double crhs22 =             DN(0,1)*crhs18 + DN(1,1)*crhs19 + DN(2,1)*crhs20 + DN(3,1)*crhs21;
const double crhs23 =             ave_c1*v_ave(0,2) - ave_c2*vn_ave(0,2);
const double crhs24 =             ave_c1*v_ave(1,2) - ave_c2*vn_ave(1,2);
const double crhs25 =             ave_c1*v_ave(2,2) - ave_c2*vn_ave(2,2);
const double crhs26 =             ave_c1*v_ave(3,2) - ave_c2*vn_ave(3,2);
const double crhs27 =             DN(0,2)*crhs23 + DN(1,2)*crhs24 + DN(2,2)*crhs25 + DN(3,2)*crhs26;
const double crhs28 =             crhs12 + crhs22 + crhs27;
const double crhs29 =             (crhs17 + crhs28)*(crhs16*h/stab_c1 + mu);
const double crhs30 =             1.0*h/(crhs16 + mu*stab_c1/h);
const double crhs31 =             crhs30*(DN(0,0)*crhs1 + DN(1,0)*crhs2 + DN(2,0)*crhs3 + DN(3,0)*crhs4 - crhs0 + crhs10 + crhs15);
const double crhs32 =             rho*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double crhs33 =             N[0]*crhs32;
const double crhs34 =             rho*(DN(0,0)*crhs11 + DN(0,1)*crhs13 + DN(0,2)*crhs14);
const double crhs35 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs36 =             rho*(N[0]*(bdf0*crhs18 + bdf1*(ave_n_c1*vn_ave(0,1) - ave_n_c2*vnn_ave(0,1)) + bdf2*(ave_nn_c1*vnn_ave(0,1) - ave_nn_c2*vnnn_ave(0,1))) + N[1]*(bdf0*crhs19 + bdf1*(ave_n_c1*vn_ave(1,1) - ave_n_c2*vnn_ave(1,1)) + bdf2*(ave_nn_c1*vnn_ave(1,1) - ave_nn_c2*vnnn_ave(1,1))) + N[2]*(bdf0*crhs20 + bdf1*(ave_n_c1*vn_ave(2,1) - ave_n_c2*vnn_ave(2,1)) + bdf2*(ave_nn_c1*vnn_ave(2,1) - ave_nn_c2*vnnn_ave(2,1))) + N[3]*(bdf0*crhs21 + bdf1*(ave_n_c1*vn_ave(3,1) - ave_n_c2*vnn_ave(3,1)) + bdf2*(ave_nn_c1*vnn_ave(3,1) - ave_nn_c2*vnnn_ave(3,1))));
const double crhs37 =             rho*(crhs11*(DN(0,0)*crhs18 + DN(1,0)*crhs19 + DN(2,0)*crhs20 + DN(3,0)*crhs21) + crhs13*crhs22 + crhs14*(DN(0,2)*crhs18 + DN(1,2)*crhs19 + DN(2,2)*crhs20 + DN(3,2)*crhs21));
const double crhs38 =             crhs30*(DN(0,1)*crhs1 + DN(1,1)*crhs2 + DN(2,1)*crhs3 + DN(3,1)*crhs4 - crhs35 + crhs36 + crhs37);
const double crhs39 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs40 =             rho*(N[0]*(bdf0*crhs23 + bdf1*(ave_n_c1*vn_ave(0,2) - ave_n_c2*vnn_ave(0,2)) + bdf2*(ave_nn_c1*vnn_ave(0,2) - ave_nn_c2*vnnn_ave(0,2))) + N[1]*(bdf0*crhs24 + bdf1*(ave_n_c1*vn_ave(1,2) - ave_n_c2*vnn_ave(1,2)) + bdf2*(ave_nn_c1*vnn_ave(1,2) - ave_nn_c2*vnnn_ave(1,2))) + N[2]*(bdf0*crhs25 + bdf1*(ave_n_c1*vn_ave(2,2) - ave_n_c2*vnn_ave(2,2)) + bdf2*(ave_nn_c1*vnn_ave(2,2) - ave_nn_c2*vnnn_ave(2,2))) + N[3]*(bdf0*crhs26 + bdf1*(ave_n_c1*vn_ave(3,2) - ave_n_c2*vnn_ave(3,2)) + bdf2*(ave_nn_c1*vnn_ave(3,2) - ave_nn_c2*vnnn_ave(3,2))));
const double crhs41 =             rho*(crhs11*(DN(0,0)*crhs23 + DN(1,0)*crhs24 + DN(2,0)*crhs25 + DN(3,0)*crhs26) + crhs13*(DN(0,1)*crhs23 + DN(1,1)*crhs24 + DN(2,1)*crhs25 + DN(3,1)*crhs26) + crhs14*crhs27);
const double crhs42 =             crhs30*(DN(0,2)*crhs1 + DN(1,2)*crhs2 + DN(2,2)*crhs3 + DN(3,2)*crhs4 - crhs39 + crhs40 + crhs41);
const double crhs43 =             N[1]*crhs32;
const double crhs44 =             rho*(DN(1,0)*crhs11 + DN(1,1)*crhs13 + DN(1,2)*crhs14);
const double crhs45 =             N[2]*crhs32;
const double crhs46 =             rho*(DN(2,0)*crhs11 + DN(2,1)*crhs13 + DN(2,2)*crhs14);
const double crhs47 =             N[3]*crhs32;
const double crhs48 =             rho*(DN(3,0)*crhs11 + DN(3,1)*crhs13 + DN(3,2)*crhs14);
            rhs[0]=-DN(0,0)*crhs29 + DN(0,0)*crhs5 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs0 - N[0]*crhs10 - N[0]*crhs15 - crhs31*crhs33 - crhs31*crhs34;
            rhs[1]=-DN(0,0)*stress[3] - DN(0,1)*crhs29 + DN(0,1)*crhs5 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs35 - N[0]*crhs36 - N[0]*crhs37 - crhs33*crhs38 - crhs34*crhs38;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] - DN(0,2)*crhs29 + DN(0,2)*crhs5 - DN(0,2)*stress[2] + N[0]*crhs39 - N[0]*crhs40 - N[0]*crhs41 - crhs33*crhs42 - crhs34*crhs42;
            rhs[3]=-DN(0,0)*crhs31 - DN(0,1)*crhs38 - DN(0,2)*crhs42 - N[0]*crhs17 - N[0]*crhs28;
            rhs[4]=-DN(1,0)*crhs29 + DN(1,0)*crhs5 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs0 - N[1]*crhs10 - N[1]*crhs15 - crhs31*crhs43 - crhs31*crhs44;
            rhs[5]=-DN(1,0)*stress[3] - DN(1,1)*crhs29 + DN(1,1)*crhs5 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs35 - N[1]*crhs36 - N[1]*crhs37 - crhs38*crhs43 - crhs38*crhs44;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] - DN(1,2)*crhs29 + DN(1,2)*crhs5 - DN(1,2)*stress[2] + N[1]*crhs39 - N[1]*crhs40 - N[1]*crhs41 - crhs42*crhs43 - crhs42*crhs44;
            rhs[7]=-DN(1,0)*crhs31 - DN(1,1)*crhs38 - DN(1,2)*crhs42 - N[1]*crhs17 - N[1]*crhs28;
            rhs[8]=-DN(2,0)*crhs29 + DN(2,0)*crhs5 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs0 - N[2]*crhs10 - N[2]*crhs15 - crhs31*crhs45 - crhs31*crhs46;
            rhs[9]=-DN(2,0)*stress[3] - DN(2,1)*crhs29 + DN(2,1)*crhs5 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs35 - N[2]*crhs36 - N[2]*crhs37 - crhs38*crhs45 - crhs38*crhs46;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] - DN(2,2)*crhs29 + DN(2,2)*crhs5 - DN(2,2)*stress[2] + N[2]*crhs39 - N[2]*crhs40 - N[2]*crhs41 - crhs42*crhs45 - crhs42*crhs46;
            rhs[11]=-DN(2,0)*crhs31 - DN(2,1)*crhs38 - DN(2,2)*crhs42 - N[2]*crhs17 - N[2]*crhs28;
            rhs[12]=-DN(3,0)*crhs29 + DN(3,0)*crhs5 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs0 - N[3]*crhs10 - N[3]*crhs15 - crhs31*crhs47 - crhs31*crhs48;
            rhs[13]=-DN(3,0)*stress[3] - DN(3,1)*crhs29 + DN(3,1)*crhs5 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs35 - N[3]*crhs36 - N[3]*crhs37 - crhs38*crhs47 - crhs38*crhs48;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] - DN(3,2)*crhs29 + DN(3,2)*crhs5 - DN(3,2)*stress[2] + N[3]*crhs39 - N[3]*crhs40 - N[3]*crhs41 - crhs42*crhs47 - crhs42*crhs48;
            rhs[15]=-DN(3,0)*crhs31 - DN(3,1)*crhs38 - DN(3,2)*crhs42 - N[3]*crhs17 - N[3]*crhs28;


}


template<>
void TimeAveragedNavierStokes<2>::ComputeGaussPointRHSContribution(
    array_1d<double,9>& rhs,
    const ElementDataStruct& data)
{
    constexpr int dim = 2;
    constexpr int nnodes = 3;
    constexpr int strain_size = 3;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size
    const double c = data.c;                                // Wave velocity

    //const double& dts = data.dts;                           // The averaging time period

    const double& dtn = data.dtn;                           // Time increment: notice t = tn + dtn

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v_ave = data.v_ave;
    const BoundedMatrix<double,nnodes,dim>& vn_ave = data.vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnn_ave = data.vnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnnn_ave = data.vnnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& p_ave = data.p_ave;
    const array_1d<double,nnodes>& pn_ave = data.pn_ave;
    const array_1d<double,nnodes>& pnn_ave = data.pnn_ave;
    const array_1d<double,nnodes>& pnnn_ave = data.pnnn_ave;
    const array_1d<double,strain_size>& stress = data.stress;

    // time averaging parameters
    double ave_c1 = data.ave_c1;
    double ave_c2 = data.ave_c2;
    double ave_n_c1 = data.ave_n_c1;
    double ave_n_c2 = data.ave_n_c2;
    double ave_nn_c1 = data.ave_nn_c1;
    double ave_nn_c2 = data.ave_nn_c2;

    // get time accurate expression
    const array_1d<double,nnodes>& p = ave_c1 * p_ave - ave_c2 * pn_ave;
    const BoundedMatrix<double,nnodes,dim>& v = ave_c1 * v_ave - ave_c2 * vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vconv = v - vmesh;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,dim> f_gauss = prod(trans(f), N);
    const array_1d<double,dim> grad_p = prod(trans(DN), p);

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    const double crhs0 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crhs1 =             ave_c1*p_ave[0] - ave_c2*pn_ave[0];
const double crhs2 =             ave_c1*p_ave[1] - ave_c2*pn_ave[1];
const double crhs3 =             ave_c1*p_ave[2] - ave_c2*pn_ave[2];
const double crhs4 =             N[0]*crhs1 + N[1]*crhs2 + N[2]*crhs3;
const double crhs5 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs6 =             ave_c1*v_ave(0,0) - ave_c2*vn_ave(0,0);
const double crhs7 =             ave_c1*v_ave(1,0) - ave_c2*vn_ave(1,0);
const double crhs8 =             ave_c1*v_ave(2,0) - ave_c2*vn_ave(2,0);
const double crhs9 =             DN(0,0)*crhs6 + DN(1,0)*crhs7 + DN(2,0)*crhs8;
const double crhs10 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs11 =             rho*(crhs10*(DN(0,1)*crhs6 + DN(1,1)*crhs7 + DN(2,1)*crhs8) + crhs5*crhs9);
const double crhs12 =             rho*(N[0]*(bdf0*crhs6 + bdf1*(ave_n_c1*vn_ave(0,0) - ave_n_c2*vnn_ave(0,0)) + bdf2*(ave_nn_c1*vnn_ave(0,0) - ave_nn_c2*vnnn_ave(0,0))) + N[1]*(bdf0*crhs7 + bdf1*(ave_n_c1*vn_ave(1,0) - ave_n_c2*vnn_ave(1,0)) + bdf2*(ave_nn_c1*vnn_ave(1,0) - ave_nn_c2*vnnn_ave(1,0))) + N[2]*(bdf0*crhs8 + bdf1*(ave_n_c1*vn_ave(2,0) - ave_n_c2*vnn_ave(2,0)) + bdf2*(ave_nn_c1*vnn_ave(2,0) - ave_nn_c2*vnnn_ave(2,0))));
const double crhs13 =             rho*stab_c2*sqrt(pow(crhs10, 2) + pow(crhs5, 2));
const double crhs14 =             (N[0]*(bdf0*crhs1 + bdf1*(ave_n_c1*pn_ave[0] - ave_n_c2*pnn_ave[0]) + bdf2*(ave_nn_c1*pnn_ave[0] - ave_nn_c2*pnnn_ave[0])) + N[1]*(bdf0*crhs2 + bdf1*(ave_n_c1*pn_ave[1] - ave_n_c2*pnn_ave[1]) + bdf2*(ave_nn_c1*pnn_ave[1] - ave_nn_c2*pnnn_ave[1])) + N[2]*(bdf0*crhs3 + bdf1*(ave_n_c1*pn_ave[2] - ave_n_c2*pnn_ave[2]) + bdf2*(ave_nn_c1*pnn_ave[2] - ave_nn_c2*pnnn_ave[2])))/(pow(c, 2)*rho);
const double crhs15 =             ave_c1*v_ave(0,1) - ave_c2*vn_ave(0,1);
const double crhs16 =             ave_c1*v_ave(1,1) - ave_c2*vn_ave(1,1);
const double crhs17 =             ave_c1*v_ave(2,1) - ave_c2*vn_ave(2,1);
const double crhs18 =             DN(0,1)*crhs15 + DN(1,1)*crhs16 + DN(2,1)*crhs17;
const double crhs19 =             crhs18 + crhs9;
const double crhs20 =             (crhs14 + crhs19)*(crhs13*h/stab_c1 + mu);
const double crhs21 =             1.0*h/(crhs13 + mu*stab_c1/h);
const double crhs22 =             crhs21*(DN(0,0)*crhs1 + DN(1,0)*crhs2 + DN(2,0)*crhs3 - crhs0 + crhs11 + crhs12);
const double crhs23 =             rho*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double crhs24 =             N[0]*crhs23;
const double crhs25 =             rho*(DN(0,0)*crhs5 + DN(0,1)*crhs10);
const double crhs26 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs27 =             rho*(crhs10*crhs18 + crhs5*(DN(0,0)*crhs15 + DN(1,0)*crhs16 + DN(2,0)*crhs17));
const double crhs28 =             rho*(N[0]*(bdf0*crhs15 + bdf1*(ave_n_c1*vn_ave(0,1) - ave_n_c2*vnn_ave(0,1)) + bdf2*(ave_nn_c1*vnn_ave(0,1) - ave_nn_c2*vnnn_ave(0,1))) + N[1]*(bdf0*crhs16 + bdf1*(ave_n_c1*vn_ave(1,1) - ave_n_c2*vnn_ave(1,1)) + bdf2*(ave_nn_c1*vnn_ave(1,1) - ave_nn_c2*vnnn_ave(1,1))) + N[2]*(bdf0*crhs17 + bdf1*(ave_n_c1*vn_ave(2,1) - ave_n_c2*vnn_ave(2,1)) + bdf2*(ave_nn_c1*vnn_ave(2,1) - ave_nn_c2*vnnn_ave(2,1))));
const double crhs29 =             crhs21*(DN(0,1)*crhs1 + DN(1,1)*crhs2 + DN(2,1)*crhs3 - crhs26 + crhs27 + crhs28);
const double crhs30 =             N[1]*crhs23;
const double crhs31 =             rho*(DN(1,0)*crhs5 + DN(1,1)*crhs10);
const double crhs32 =             N[2]*crhs23;
const double crhs33 =             rho*(DN(2,0)*crhs5 + DN(2,1)*crhs10);
            rhs[0]=-DN(0,0)*crhs20 + DN(0,0)*crhs4 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs0 - N[0]*crhs11 - N[0]*crhs12 - crhs22*crhs24 - crhs22*crhs25;
            rhs[1]=-DN(0,0)*stress[2] - DN(0,1)*crhs20 + DN(0,1)*crhs4 - DN(0,1)*stress[1] + N[0]*crhs26 - N[0]*crhs27 - N[0]*crhs28 - crhs24*crhs29 - crhs25*crhs29;
            rhs[2]=-DN(0,0)*crhs22 - DN(0,1)*crhs29 - N[0]*crhs14 - N[0]*crhs19;
            rhs[3]=-DN(1,0)*crhs20 + DN(1,0)*crhs4 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs0 - N[1]*crhs11 - N[1]*crhs12 - crhs22*crhs30 - crhs22*crhs31;
            rhs[4]=-DN(1,0)*stress[2] - DN(1,1)*crhs20 + DN(1,1)*crhs4 - DN(1,1)*stress[1] + N[1]*crhs26 - N[1]*crhs27 - N[1]*crhs28 - crhs29*crhs30 - crhs29*crhs31;
            rhs[5]=-DN(1,0)*crhs22 - DN(1,1)*crhs29 - N[1]*crhs14 - N[1]*crhs19;
            rhs[6]=-DN(2,0)*crhs20 + DN(2,0)*crhs4 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs0 - N[2]*crhs11 - N[2]*crhs12 - crhs22*crhs32 - crhs22*crhs33;
            rhs[7]=-DN(2,0)*stress[2] - DN(2,1)*crhs20 + DN(2,1)*crhs4 - DN(2,1)*stress[1] + N[2]*crhs26 - N[2]*crhs27 - N[2]*crhs28 - crhs29*crhs32 - crhs29*crhs33;
            rhs[8]=-DN(2,0)*crhs22 - DN(2,1)*crhs29 - N[2]*crhs14 - N[2]*crhs19;


}


template<>
double TimeAveragedNavierStokes<3>::SubscaleErrorEstimate(const ElementDataStruct& data)
{
    constexpr int dim = 3;
    constexpr int nnodes = 4;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size
    // const double c = data.c;                                // Wave velocity

    //const double& dts = data.dts;                           // The averaging time period
    const double& dtn = data.dtn;                           // Time increment: notice t = tn + dtn

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v_ave = data.v_ave;
    const BoundedMatrix<double,nnodes,dim>& vn_ave = data.vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnn_ave = data.vnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnnn_ave = data.vnnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& p_ave = data.p_ave;
    const array_1d<double,nnodes>& pn_ave = data.pn_ave;

    // time averaging parameters
    double ave_c1 = data.ave_c1;
    double ave_c2 = data.ave_c2;
    double ave_n_c1 = data.ave_n_c1;
    double ave_n_c2 = data.ave_n_c2;
    double ave_nn_c1 = data.ave_nn_c1;
    double ave_nn_c2 = data.ave_nn_c2;

    // get time accurate expression
    const array_1d<double,nnodes>& p = data.p;
    const BoundedMatrix<double,nnodes,dim>& v = data.v;
    const BoundedMatrix<double,nnodes,dim>& vconv = v - vmesh;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the error estimator
    array_1d<double,dim> v_s_gauss;
    const array_1d<double,dim> v_gauss = prod(trans(v), N);
    const array_1d<double,dim> f_gauss = prod(trans(f), N);
    const array_1d<double,dim> grad_p = prod(trans(DN), p);

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Gauss point velocity subscale value computation
    const double cv_s_gauss0 =             ave_c1*p_ave[0] - ave_c2*pn_ave[0];
const double cv_s_gauss1 =             ave_c1*p_ave[1] - ave_c2*pn_ave[1];
const double cv_s_gauss2 =             ave_c1*p_ave[2] - ave_c2*pn_ave[2];
const double cv_s_gauss3 =             ave_c1*p_ave[3] - ave_c2*pn_ave[3];
const double cv_s_gauss4 =             ave_c1*v_ave(0,0) - ave_c2*vn_ave(0,0);
const double cv_s_gauss5 =             ave_c1*v_ave(1,0) - ave_c2*vn_ave(1,0);
const double cv_s_gauss6 =             ave_c1*v_ave(2,0) - ave_c2*vn_ave(2,0);
const double cv_s_gauss7 =             ave_c1*v_ave(3,0) - ave_c2*vn_ave(3,0);
const double cv_s_gauss8 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cv_s_gauss9 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cv_s_gauss10 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cv_s_gauss11 =             1.0*h/(rho*stab_c2*sqrt(pow(cv_s_gauss10, 2) + pow(cv_s_gauss8, 2) + pow(cv_s_gauss9, 2)) + mu*stab_c1/h);
const double cv_s_gauss12 =             ave_c1*v_ave(0,1) - ave_c2*vn_ave(0,1);
const double cv_s_gauss13 =             ave_c1*v_ave(1,1) - ave_c2*vn_ave(1,1);
const double cv_s_gauss14 =             ave_c1*v_ave(2,1) - ave_c2*vn_ave(2,1);
const double cv_s_gauss15 =             ave_c1*v_ave(3,1) - ave_c2*vn_ave(3,1);
const double cv_s_gauss16 =             ave_c1*v_ave(0,2) - ave_c2*vn_ave(0,2);
const double cv_s_gauss17 =             ave_c1*v_ave(1,2) - ave_c2*vn_ave(1,2);
const double cv_s_gauss18 =             ave_c1*v_ave(2,2) - ave_c2*vn_ave(2,2);
const double cv_s_gauss19 =             ave_c1*v_ave(3,2) - ave_c2*vn_ave(3,2);
            v_s_gauss[0]=-cv_s_gauss11*(DN(0,0)*cv_s_gauss0 + DN(1,0)*cv_s_gauss1 + DN(2,0)*cv_s_gauss2 + DN(3,0)*cv_s_gauss3 + rho*(-N[0]*f(0,0) + N[0]*(bdf0*cv_s_gauss4 + bdf1*(ave_n_c1*vn_ave(0,0) - ave_n_c2*vnn_ave(0,0)) + bdf2*(ave_nn_c1*vnn_ave(0,0) - ave_nn_c2*vnnn_ave(0,0))) - N[1]*f(1,0) + N[1]*(bdf0*cv_s_gauss5 + bdf1*(ave_n_c1*vn_ave(1,0) - ave_n_c2*vnn_ave(1,0)) + bdf2*(ave_nn_c1*vnn_ave(1,0) - ave_nn_c2*vnnn_ave(1,0))) - N[2]*f(2,0) + N[2]*(bdf0*cv_s_gauss6 + bdf1*(ave_n_c1*vn_ave(2,0) - ave_n_c2*vnn_ave(2,0)) + bdf2*(ave_nn_c1*vnn_ave(2,0) - ave_nn_c2*vnnn_ave(2,0))) - N[3]*f(3,0) + N[3]*(bdf0*cv_s_gauss7 + bdf1*(ave_n_c1*vn_ave(3,0) - ave_n_c2*vnn_ave(3,0)) + bdf2*(ave_nn_c1*vnn_ave(3,0) - ave_nn_c2*vnnn_ave(3,0))) + cv_s_gauss10*(DN(0,2)*cv_s_gauss4 + DN(1,2)*cv_s_gauss5 + DN(2,2)*cv_s_gauss6 + DN(3,2)*cv_s_gauss7) + cv_s_gauss8*(DN(0,0)*cv_s_gauss4 + DN(1,0)*cv_s_gauss5 + DN(2,0)*cv_s_gauss6 + DN(3,0)*cv_s_gauss7) + cv_s_gauss9*(DN(0,1)*cv_s_gauss4 + DN(1,1)*cv_s_gauss5 + DN(2,1)*cv_s_gauss6 + DN(3,1)*cv_s_gauss7)));
            v_s_gauss[1]=-cv_s_gauss11*(DN(0,1)*cv_s_gauss0 + DN(1,1)*cv_s_gauss1 + DN(2,1)*cv_s_gauss2 + DN(3,1)*cv_s_gauss3 + rho*(-N[0]*f(0,1) + N[0]*(bdf0*cv_s_gauss12 + bdf1*(ave_n_c1*vn_ave(0,1) - ave_n_c2*vnn_ave(0,1)) + bdf2*(ave_nn_c1*vnn_ave(0,1) - ave_nn_c2*vnnn_ave(0,1))) - N[1]*f(1,1) + N[1]*(bdf0*cv_s_gauss13 + bdf1*(ave_n_c1*vn_ave(1,1) - ave_n_c2*vnn_ave(1,1)) + bdf2*(ave_nn_c1*vnn_ave(1,1) - ave_nn_c2*vnnn_ave(1,1))) - N[2]*f(2,1) + N[2]*(bdf0*cv_s_gauss14 + bdf1*(ave_n_c1*vn_ave(2,1) - ave_n_c2*vnn_ave(2,1)) + bdf2*(ave_nn_c1*vnn_ave(2,1) - ave_nn_c2*vnnn_ave(2,1))) - N[3]*f(3,1) + N[3]*(bdf0*cv_s_gauss15 + bdf1*(ave_n_c1*vn_ave(3,1) - ave_n_c2*vnn_ave(3,1)) + bdf2*(ave_nn_c1*vnn_ave(3,1) - ave_nn_c2*vnnn_ave(3,1))) + cv_s_gauss10*(DN(0,2)*cv_s_gauss12 + DN(1,2)*cv_s_gauss13 + DN(2,2)*cv_s_gauss14 + DN(3,2)*cv_s_gauss15) + cv_s_gauss8*(DN(0,0)*cv_s_gauss12 + DN(1,0)*cv_s_gauss13 + DN(2,0)*cv_s_gauss14 + DN(3,0)*cv_s_gauss15) + cv_s_gauss9*(DN(0,1)*cv_s_gauss12 + DN(1,1)*cv_s_gauss13 + DN(2,1)*cv_s_gauss14 + DN(3,1)*cv_s_gauss15)));
            v_s_gauss[2]=-cv_s_gauss11*(DN(0,2)*cv_s_gauss0 + DN(1,2)*cv_s_gauss1 + DN(2,2)*cv_s_gauss2 + DN(3,2)*cv_s_gauss3 + rho*(-N[0]*f(0,2) + N[0]*(bdf0*cv_s_gauss16 + bdf1*(ave_n_c1*vn_ave(0,2) - ave_n_c2*vnn_ave(0,2)) + bdf2*(ave_nn_c1*vnn_ave(0,2) - ave_nn_c2*vnnn_ave(0,2))) - N[1]*f(1,2) + N[1]*(bdf0*cv_s_gauss17 + bdf1*(ave_n_c1*vn_ave(1,2) - ave_n_c2*vnn_ave(1,2)) + bdf2*(ave_nn_c1*vnn_ave(1,2) - ave_nn_c2*vnnn_ave(1,2))) - N[2]*f(2,2) + N[2]*(bdf0*cv_s_gauss18 + bdf1*(ave_n_c1*vn_ave(2,2) - ave_n_c2*vnn_ave(2,2)) + bdf2*(ave_nn_c1*vnn_ave(2,2) - ave_nn_c2*vnnn_ave(2,2))) - N[3]*f(3,2) + N[3]*(bdf0*cv_s_gauss19 + bdf1*(ave_n_c1*vn_ave(3,2) - ave_n_c2*vnn_ave(3,2)) + bdf2*(ave_nn_c1*vnn_ave(3,2) - ave_nn_c2*vnnn_ave(3,2))) + cv_s_gauss10*(DN(0,2)*cv_s_gauss16 + DN(1,2)*cv_s_gauss17 + DN(2,2)*cv_s_gauss18 + DN(3,2)*cv_s_gauss19) + cv_s_gauss8*(DN(0,0)*cv_s_gauss16 + DN(1,0)*cv_s_gauss17 + DN(2,0)*cv_s_gauss18 + DN(3,0)*cv_s_gauss19) + cv_s_gauss9*(DN(0,1)*cv_s_gauss16 + DN(1,1)*cv_s_gauss17 + DN(2,1)*cv_s_gauss18 + DN(3,1)*cv_s_gauss19)));


    const double v_gauss_norm = norm_2(v_gauss);
    const double v_s_gauss_norm = norm_2(v_s_gauss);

    return v_s_gauss_norm/v_gauss_norm;
}


template<>
double TimeAveragedNavierStokes<2>::SubscaleErrorEstimate(const ElementDataStruct& data)
{
    constexpr int dim = 2;
    constexpr int nnodes = 3;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size

    //const double& dts = data.dts;                           // The averaging time period
    const double& dtn = data.dtn;                           // Time increment: notice t = tn + dtn

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v_ave    = data.v_ave;
    const BoundedMatrix<double,nnodes,dim>& vn_ave   = data.vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnn_ave  = data.vnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnnn_ave = data.vnnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& p_ave = data.p_ave;
    const array_1d<double,nnodes>& pn_ave = data.pn_ave;

    // time averaging parameters
    double ave_c1 = data.ave_c1;
    double ave_c2 = data.ave_c2;
    double ave_n_c1 = data.ave_n_c1;
    double ave_n_c2 = data.ave_n_c2;
    double ave_nn_c1 = data.ave_nn_c1;
    double ave_nn_c2 = data.ave_nn_c2;

    // get time accurate expression
    const array_1d<double,nnodes>& p = data.p;
    const BoundedMatrix<double,nnodes,dim>& v = data.v;
    const BoundedMatrix<double,nnodes,dim>& vconv = v - vmesh;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the error estimator
    array_1d<double,dim> v_s_gauss;
    const array_1d<double,dim> v_gauss = prod(trans(v), N);
    const array_1d<double,dim> f_gauss = prod(trans(f), N);
    const array_1d<double,dim> grad_p = prod(trans(DN), p);

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Gauss point velocity subscale value computation
    const double cv_s_gauss0 =             ave_c1*p_ave[0] - ave_c2*pn_ave[0];
const double cv_s_gauss1 =             ave_c1*p_ave[1] - ave_c2*pn_ave[1];
const double cv_s_gauss2 =             ave_c1*p_ave[2] - ave_c2*pn_ave[2];
const double cv_s_gauss3 =             ave_c1*v_ave(0,0) - ave_c2*vn_ave(0,0);
const double cv_s_gauss4 =             ave_c1*v_ave(1,0) - ave_c2*vn_ave(1,0);
const double cv_s_gauss5 =             ave_c1*v_ave(2,0) - ave_c2*vn_ave(2,0);
const double cv_s_gauss6 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cv_s_gauss7 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cv_s_gauss8 =             1.0*h/(rho*stab_c2*sqrt(pow(cv_s_gauss6, 2) + pow(cv_s_gauss7, 2)) + mu*stab_c1/h);
const double cv_s_gauss9 =             ave_c1*v_ave(0,1) - ave_c2*vn_ave(0,1);
const double cv_s_gauss10 =             ave_c1*v_ave(1,1) - ave_c2*vn_ave(1,1);
const double cv_s_gauss11 =             ave_c1*v_ave(2,1) - ave_c2*vn_ave(2,1);
            v_s_gauss[0]=-cv_s_gauss8*(DN(0,0)*cv_s_gauss0 + DN(1,0)*cv_s_gauss1 + DN(2,0)*cv_s_gauss2 + rho*(-N[0]*f(0,0) + N[0]*(bdf0*cv_s_gauss3 + bdf1*(ave_n_c1*vn_ave(0,0) - ave_n_c2*vnn_ave(0,0)) + bdf2*(ave_nn_c1*vnn_ave(0,0) - ave_nn_c2*vnnn_ave(0,0))) - N[1]*f(1,0) + N[1]*(bdf0*cv_s_gauss4 + bdf1*(ave_n_c1*vn_ave(1,0) - ave_n_c2*vnn_ave(1,0)) + bdf2*(ave_nn_c1*vnn_ave(1,0) - ave_nn_c2*vnnn_ave(1,0))) - N[2]*f(2,0) + N[2]*(bdf0*cv_s_gauss5 + bdf1*(ave_n_c1*vn_ave(2,0) - ave_n_c2*vnn_ave(2,0)) + bdf2*(ave_nn_c1*vnn_ave(2,0) - ave_nn_c2*vnnn_ave(2,0))) + cv_s_gauss6*(DN(0,0)*cv_s_gauss3 + DN(1,0)*cv_s_gauss4 + DN(2,0)*cv_s_gauss5) + cv_s_gauss7*(DN(0,1)*cv_s_gauss3 + DN(1,1)*cv_s_gauss4 + DN(2,1)*cv_s_gauss5)));
            v_s_gauss[1]=-cv_s_gauss8*(DN(0,1)*cv_s_gauss0 + DN(1,1)*cv_s_gauss1 + DN(2,1)*cv_s_gauss2 + rho*(-N[0]*f(0,1) + N[0]*(bdf0*cv_s_gauss9 + bdf1*(ave_n_c1*vn_ave(0,1) - ave_n_c2*vnn_ave(0,1)) + bdf2*(ave_nn_c1*vnn_ave(0,1) - ave_nn_c2*vnnn_ave(0,1))) - N[1]*f(1,1) + N[1]*(bdf0*cv_s_gauss10 + bdf1*(ave_n_c1*vn_ave(1,1) - ave_n_c2*vnn_ave(1,1)) + bdf2*(ave_nn_c1*vnn_ave(1,1) - ave_nn_c2*vnnn_ave(1,1))) - N[2]*f(2,1) + N[2]*(bdf0*cv_s_gauss11 + bdf1*(ave_n_c1*vn_ave(2,1) - ave_n_c2*vnn_ave(2,1)) + bdf2*(ave_nn_c1*vnn_ave(2,1) - ave_nn_c2*vnnn_ave(2,1))) + cv_s_gauss6*(DN(0,0)*cv_s_gauss9 + DN(1,0)*cv_s_gauss10 + DN(2,0)*cv_s_gauss11) + cv_s_gauss7*(DN(0,1)*cv_s_gauss9 + DN(1,1)*cv_s_gauss10 + DN(2,1)*cv_s_gauss11)));


    const double v_gauss_norm = norm_2(v_gauss);
    const double v_s_gauss_norm = norm_2(v_s_gauss);

    return v_s_gauss_norm/v_gauss_norm;
}

}
