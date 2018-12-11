//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
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

    KRATOS_CATCH("")
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

    KRATOS_CATCH("")
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

    KRATOS_CATCH("");
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

    KRATOS_CATCH("");
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

    const double& dt = data.dt;
    const double& bdf0 = data.bdf0;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v_ave = data.v_ave;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& vconv = v_ave - vmesh;

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
const double clhs11 =             pow(N[0], 2);
const double clhs12 =             bdf0*n*rho;
const double clhs13 =             N[0]*rho;
const double clhs14 =             DN(0,0)*clhs6 + DN(0,1)*clhs7 + DN(0,2)*clhs8;
const double clhs15 =             N[0]*n;
const double clhs16 =             bdf0*clhs15 + clhs14;
const double clhs17 =             pow(rho, 2);
const double clhs18 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double clhs19 =             1.0/(clhs9/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs20 =             1.0*N[0]*clhs17*clhs18*clhs19;
const double clhs21 =             1.0*clhs14*clhs17*clhs19;
const double clhs22 =             clhs11*clhs12 + clhs13*clhs14 + clhs16*clhs20 + clhs16*clhs21;
const double clhs23 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs1;
const double clhs24 =             C(1,3)*DN(0,1);
const double clhs25 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs24;
const double clhs26 =             C(3,5)*DN(0,0);
const double clhs27 =             C(4,5)*DN(0,2);
const double clhs28 =             C(1,5)*DN(0,1) + clhs26 + clhs27;
const double clhs29 =             DN(0,0)*clhs10;
const double clhs30 =             DN(0,1)*clhs29;
const double clhs31 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs3;
const double clhs32 =             C(3,4)*DN(0,1);
const double clhs33 =             C(2,3)*DN(0,2) + clhs26 + clhs32;
const double clhs34 =             C(2,5)*DN(0,2);
const double clhs35 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs34;
const double clhs36 =             DN(0,2)*clhs29;
const double clhs37 =             pow(c, -2);
const double clhs38 =             1.0/rho;
const double clhs39 =             N[0]*bdf0*clhs37*clhs38*n;
const double clhs40 =             1.0*clhs18*clhs19;
const double clhs41 =             1.0*clhs19*rho;
const double clhs42 =             clhs14*clhs41;
const double clhs43 =             clhs10*clhs39 + clhs13*clhs40 - clhs15 + clhs42;
const double clhs44 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs45 =             C(0,3)*DN(1,0);
const double clhs46 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs45;
const double clhs47 =             C(0,5)*DN(1,0);
const double clhs48 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs47;
const double clhs49 =             DN(1,0)*clhs29;
const double clhs50 =             N[0]*bdf0*n*rho;
const double clhs51 =             N[1]*clhs50;
const double clhs52 =             DN(1,0)*clhs6 + DN(1,1)*clhs7 + DN(1,2)*clhs8;
const double clhs53 =             N[1]*n;
const double clhs54 =             bdf0*clhs53 + clhs52;
const double clhs55 =             clhs13*clhs52 + clhs20*clhs54 + clhs21*clhs54 + clhs51;
const double clhs56 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs45;
const double clhs57 =             C(1,3)*DN(1,1);
const double clhs58 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs57;
const double clhs59 =             C(3,5)*DN(1,0);
const double clhs60 =             C(4,5)*DN(1,2);
const double clhs61 =             C(1,5)*DN(1,1) + clhs59 + clhs60;
const double clhs62 =             DN(1,1)*clhs29;
const double clhs63 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs47;
const double clhs64 =             C(3,4)*DN(1,1);
const double clhs65 =             C(2,3)*DN(1,2) + clhs59 + clhs64;
const double clhs66 =             C(2,5)*DN(1,2);
const double clhs67 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs66;
const double clhs68 =             DN(1,2)*clhs29;
const double clhs69 =             DN(0,0)*N[1];
const double clhs70 =             clhs69*n;
const double clhs71 =             bdf0*clhs10*clhs37*clhs38;
const double clhs72 =             DN(1,0)*N[0];
const double clhs73 =             1.0*clhs18*clhs19*rho;
const double clhs74 =             1.0*DN(1,0)*clhs19*rho;
const double clhs75 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs76 =             C(0,3)*DN(2,0);
const double clhs77 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs76;
const double clhs78 =             C(0,5)*DN(2,0);
const double clhs79 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs78;
const double clhs80 =             DN(2,0)*clhs29;
const double clhs81 =             N[2]*clhs50;
const double clhs82 =             DN(2,0)*clhs6 + DN(2,1)*clhs7 + DN(2,2)*clhs8;
const double clhs83 =             N[2]*n;
const double clhs84 =             bdf0*clhs83;
const double clhs85 =             clhs82 + clhs84;
const double clhs86 =             clhs13*clhs82 + clhs20*clhs85 + clhs21*clhs85 + clhs81;
const double clhs87 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs76;
const double clhs88 =             C(1,3)*DN(2,1);
const double clhs89 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs88;
const double clhs90 =             C(3,5)*DN(2,0);
const double clhs91 =             C(4,5)*DN(2,2);
const double clhs92 =             C(1,5)*DN(2,1) + clhs90 + clhs91;
const double clhs93 =             DN(2,1)*clhs29;
const double clhs94 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs78;
const double clhs95 =             C(3,4)*DN(2,1);
const double clhs96 =             C(2,3)*DN(2,2) + clhs90 + clhs95;
const double clhs97 =             C(2,5)*DN(2,2);
const double clhs98 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs97;
const double clhs99 =             DN(2,2)*clhs29;
const double clhs100 =             DN(0,0)*N[2];
const double clhs101 =             clhs100*n;
const double clhs102 =             DN(2,0)*N[0];
const double clhs103 =             1.0*DN(2,0)*clhs19*rho;
const double clhs104 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs105 =             C(0,3)*DN(3,0);
const double clhs106 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs105;
const double clhs107 =             C(0,5)*DN(3,0);
const double clhs108 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs107;
const double clhs109 =             DN(3,0)*clhs29;
const double clhs110 =             N[3]*clhs50;
const double clhs111 =             DN(3,0)*clhs6 + DN(3,1)*clhs7 + DN(3,2)*clhs8;
const double clhs112 =             N[3]*n;
const double clhs113 =             bdf0*clhs112;
const double clhs114 =             clhs111 + clhs113;
const double clhs115 =             clhs110 + clhs111*clhs13 + clhs114*clhs20 + clhs114*clhs21;
const double clhs116 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs105;
const double clhs117 =             C(1,3)*DN(3,1);
const double clhs118 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs117;
const double clhs119 =             C(3,5)*DN(3,0);
const double clhs120 =             C(4,5)*DN(3,2);
const double clhs121 =             C(1,5)*DN(3,1) + clhs119 + clhs120;
const double clhs122 =             DN(3,1)*clhs29;
const double clhs123 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs107;
const double clhs124 =             C(3,4)*DN(3,1);
const double clhs125 =             C(2,3)*DN(3,2) + clhs119 + clhs124;
const double clhs126 =             C(2,5)*DN(3,2);
const double clhs127 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs126;
const double clhs128 =             DN(3,2)*clhs29;
const double clhs129 =             DN(0,0)*N[3];
const double clhs130 =             clhs129*n;
const double clhs131 =             DN(3,0)*N[0];
const double clhs132 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs24;
const double clhs133 =             C(0,4)*DN(0,0) + clhs27 + clhs32;
const double clhs134 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs135 =             C(1,4)*DN(0,1);
const double clhs136 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs135;
const double clhs137 =             pow(DN(0,1), 2);
const double clhs138 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs135;
const double clhs139 =             C(2,4)*DN(0,2);
const double clhs140 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs139;
const double clhs141 =             DN(0,1)*clhs10;
const double clhs142 =             DN(0,2)*clhs141;
const double clhs143 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs57;
const double clhs144 =             C(0,4)*DN(1,0) + clhs60 + clhs64;
const double clhs145 =             DN(1,0)*clhs141;
const double clhs146 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs147 =             C(1,4)*DN(1,1);
const double clhs148 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs147;
const double clhs149 =             DN(1,1)*clhs141;
const double clhs150 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs147;
const double clhs151 =             C(2,4)*DN(1,2);
const double clhs152 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs151;
const double clhs153 =             DN(1,2)*clhs141;
const double clhs154 =             DN(0,1)*N[1];
const double clhs155 =             clhs154*n;
const double clhs156 =             DN(1,1)*N[0];
const double clhs157 =             1.0*DN(1,1)*clhs19*rho;
const double clhs158 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs88;
const double clhs159 =             C(0,4)*DN(2,0) + clhs91 + clhs95;
const double clhs160 =             DN(2,0)*clhs141;
const double clhs161 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs162 =             C(1,4)*DN(2,1);
const double clhs163 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs162;
const double clhs164 =             DN(2,1)*clhs141;
const double clhs165 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs162;
const double clhs166 =             C(2,4)*DN(2,2);
const double clhs167 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs166;
const double clhs168 =             DN(2,2)*clhs141;
const double clhs169 =             DN(0,1)*N[2];
const double clhs170 =             clhs169*n;
const double clhs171 =             DN(2,1)*N[0];
const double clhs172 =             1.0*DN(2,1)*clhs19*rho;
const double clhs173 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs117;
const double clhs174 =             C(0,4)*DN(3,0) + clhs120 + clhs124;
const double clhs175 =             DN(3,0)*clhs141;
const double clhs176 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs177 =             C(1,4)*DN(3,1);
const double clhs178 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs177;
const double clhs179 =             DN(3,1)*clhs141;
const double clhs180 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs177;
const double clhs181 =             C(2,4)*DN(3,2);
const double clhs182 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs181;
const double clhs183 =             DN(3,2)*clhs141;
const double clhs184 =             DN(0,1)*N[3];
const double clhs185 =             clhs184*n;
const double clhs186 =             DN(3,1)*N[0];
const double clhs187 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs34;
const double clhs188 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs139;
const double clhs189 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs190 =             pow(DN(0,2), 2);
const double clhs191 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs66;
const double clhs192 =             DN(0,2)*clhs10;
const double clhs193 =             DN(1,0)*clhs192;
const double clhs194 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs151;
const double clhs195 =             DN(1,1)*clhs192;
const double clhs196 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs197 =             DN(1,2)*clhs192;
const double clhs198 =             DN(0,2)*N[1];
const double clhs199 =             clhs198*n;
const double clhs200 =             DN(1,2)*N[0];
const double clhs201 =             1.0*DN(1,2)*clhs19*rho;
const double clhs202 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs97;
const double clhs203 =             DN(2,0)*clhs192;
const double clhs204 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs166;
const double clhs205 =             DN(2,1)*clhs192;
const double clhs206 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs207 =             DN(2,2)*clhs192;
const double clhs208 =             DN(0,2)*N[2];
const double clhs209 =             clhs208*n;
const double clhs210 =             DN(2,2)*N[0];
const double clhs211 =             1.0*DN(2,2)*clhs19*rho;
const double clhs212 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs126;
const double clhs213 =             DN(3,0)*clhs192;
const double clhs214 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs181;
const double clhs215 =             DN(3,1)*clhs192;
const double clhs216 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs217 =             DN(3,2)*clhs192;
const double clhs218 =             DN(0,2)*N[3];
const double clhs219 =             clhs218*n;
const double clhs220 =             DN(3,2)*N[0];
const double clhs221 =             clhs16*clhs41;
const double clhs222 =             N[0] + clhs221;
const double clhs223 =             bdf0*clhs37*clhs38*n;
const double clhs224 =             1.0*clhs19;
const double clhs225 =             1.0*DN(0,0)*clhs19*rho;
const double clhs226 =             1.0*DN(0,1)*clhs19*rho;
const double clhs227 =             1.0*DN(0,2)*clhs19*rho;
const double clhs228 =             1.0*DN(0,0)*clhs19;
const double clhs229 =             1.0*DN(0,1)*clhs19;
const double clhs230 =             1.0*DN(0,2)*clhs19;
const double clhs231 =             DN(1,0)*clhs228 + DN(1,1)*clhs229 + DN(1,2)*clhs230 + N[1]*clhs39;
const double clhs232 =             DN(2,0)*clhs228 + DN(2,1)*clhs229 + DN(2,2)*clhs230 + N[2]*clhs39;
const double clhs233 =             DN(3,0)*clhs228 + DN(3,1)*clhs229 + DN(3,2)*clhs230 + N[3]*clhs39;
const double clhs234 =             N[1]*rho;
const double clhs235 =             1.0*N[1]*clhs17*clhs18*clhs19;
const double clhs236 =             1.0*clhs17*clhs19*clhs52;
const double clhs237 =             clhs14*clhs234 + clhs16*clhs235 + clhs16*clhs236 + clhs51;
const double clhs238 =             clhs72*n;
const double clhs239 =             pow(DN(1,0), 2);
const double clhs240 =             pow(N[1], 2);
const double clhs241 =             clhs12*clhs240 + clhs234*clhs52 + clhs235*clhs54 + clhs236*clhs54;
const double clhs242 =             DN(1,0)*clhs10;
const double clhs243 =             DN(1,1)*clhs242;
const double clhs244 =             DN(1,2)*clhs242;
const double clhs245 =             N[1]*bdf0*clhs37*clhs38*n;
const double clhs246 =             clhs41*clhs52;
const double clhs247 =             clhs10*clhs245 + clhs234*clhs40 + clhs246 - clhs53;
const double clhs248 =             DN(2,0)*clhs242;
const double clhs249 =             N[1]*bdf0*n*rho;
const double clhs250 =             N[2]*clhs249;
const double clhs251 =             clhs234*clhs82 + clhs235*clhs85 + clhs236*clhs85 + clhs250;
const double clhs252 =             DN(2,1)*clhs242;
const double clhs253 =             DN(2,2)*clhs242;
const double clhs254 =             DN(1,0)*N[2];
const double clhs255 =             clhs254*n;
const double clhs256 =             DN(2,0)*N[1];
const double clhs257 =             DN(3,0)*clhs242;
const double clhs258 =             N[3]*clhs249;
const double clhs259 =             clhs111*clhs234 + clhs114*clhs235 + clhs114*clhs236 + clhs258;
const double clhs260 =             DN(3,1)*clhs242;
const double clhs261 =             DN(3,2)*clhs242;
const double clhs262 =             DN(1,0)*N[3];
const double clhs263 =             clhs262*n;
const double clhs264 =             DN(3,0)*N[1];
const double clhs265 =             clhs156*n;
const double clhs266 =             pow(DN(1,1), 2);
const double clhs267 =             DN(1,1)*clhs10;
const double clhs268 =             DN(1,2)*clhs267;
const double clhs269 =             DN(2,0)*clhs267;
const double clhs270 =             DN(2,1)*clhs267;
const double clhs271 =             DN(2,2)*clhs267;
const double clhs272 =             DN(1,1)*N[2];
const double clhs273 =             clhs272*n;
const double clhs274 =             DN(2,1)*N[1];
const double clhs275 =             DN(3,0)*clhs267;
const double clhs276 =             DN(3,1)*clhs267;
const double clhs277 =             DN(3,2)*clhs267;
const double clhs278 =             DN(1,1)*N[3];
const double clhs279 =             clhs278*n;
const double clhs280 =             DN(3,1)*N[1];
const double clhs281 =             clhs200*n;
const double clhs282 =             pow(DN(1,2), 2);
const double clhs283 =             DN(1,2)*clhs10;
const double clhs284 =             DN(2,0)*clhs283;
const double clhs285 =             DN(2,1)*clhs283;
const double clhs286 =             DN(2,2)*clhs283;
const double clhs287 =             DN(1,2)*N[2];
const double clhs288 =             clhs287*n;
const double clhs289 =             DN(2,2)*N[1];
const double clhs290 =             DN(3,0)*clhs283;
const double clhs291 =             DN(3,1)*clhs283;
const double clhs292 =             DN(3,2)*clhs283;
const double clhs293 =             DN(1,2)*N[3];
const double clhs294 =             clhs293*n;
const double clhs295 =             DN(3,2)*N[1];
const double clhs296 =             clhs41*clhs54;
const double clhs297 =             N[1] + clhs296;
const double clhs298 =             1.0*DN(1,0)*clhs19;
const double clhs299 =             1.0*DN(1,1)*clhs19;
const double clhs300 =             1.0*DN(1,2)*clhs19;
const double clhs301 =             DN(2,0)*clhs298 + DN(2,1)*clhs299 + DN(2,2)*clhs300 + N[2]*clhs245;
const double clhs302 =             DN(3,0)*clhs298 + DN(3,1)*clhs299 + DN(3,2)*clhs300 + N[3]*clhs245;
const double clhs303 =             N[2]*rho;
const double clhs304 =             1.0*N[2]*clhs17*clhs18*clhs19;
const double clhs305 =             1.0*clhs17*clhs19*clhs82;
const double clhs306 =             clhs14*clhs303 + clhs16*clhs304 + clhs16*clhs305 + clhs81;
const double clhs307 =             clhs102*n;
const double clhs308 =             clhs250 + clhs303*clhs52 + clhs304*clhs54 + clhs305*clhs54;
const double clhs309 =             clhs256*n;
const double clhs310 =             pow(DN(2,0), 2);
const double clhs311 =             pow(N[2], 2);
const double clhs312 =             clhs12*clhs311 + clhs303*clhs82 + clhs304*clhs85 + clhs305*clhs85;
const double clhs313 =             DN(2,0)*clhs10;
const double clhs314 =             DN(2,1)*clhs313;
const double clhs315 =             DN(2,2)*clhs313;
const double clhs316 =             clhs10*clhs37*clhs38;
const double clhs317 =             clhs41*clhs82;
const double clhs318 =             clhs303*clhs40 + clhs316*clhs84 + clhs317 - clhs83;
const double clhs319 =             DN(3,0)*clhs313;
const double clhs320 =             N[2]*N[3]*bdf0*n;
const double clhs321 =             clhs320*rho;
const double clhs322 =             clhs111*clhs303 + clhs114*clhs304 + clhs114*clhs305 + clhs321;
const double clhs323 =             DN(3,1)*clhs313;
const double clhs324 =             DN(3,2)*clhs313;
const double clhs325 =             DN(2,0)*N[3];
const double clhs326 =             clhs325*n;
const double clhs327 =             DN(3,0)*N[2];
const double clhs328 =             clhs171*n;
const double clhs329 =             clhs274*n;
const double clhs330 =             pow(DN(2,1), 2);
const double clhs331 =             DN(2,1)*clhs10;
const double clhs332 =             DN(2,2)*clhs331;
const double clhs333 =             DN(3,0)*clhs331;
const double clhs334 =             DN(3,1)*clhs331;
const double clhs335 =             DN(3,2)*clhs331;
const double clhs336 =             DN(2,1)*N[3];
const double clhs337 =             clhs336*n;
const double clhs338 =             DN(3,1)*N[2];
const double clhs339 =             clhs210*n;
const double clhs340 =             clhs289*n;
const double clhs341 =             pow(DN(2,2), 2);
const double clhs342 =             DN(2,2)*clhs10;
const double clhs343 =             DN(3,0)*clhs342;
const double clhs344 =             DN(3,1)*clhs342;
const double clhs345 =             DN(3,2)*clhs342;
const double clhs346 =             DN(2,2)*N[3];
const double clhs347 =             clhs346*n;
const double clhs348 =             DN(3,2)*N[2];
const double clhs349 =             clhs41*clhs85;
const double clhs350 =             N[2] + clhs349;
const double clhs351 =             1.0*DN(2,0)*DN(3,0)*clhs19 + 1.0*DN(2,1)*DN(3,1)*clhs19 + 1.0*DN(2,2)*DN(3,2)*clhs19 + clhs320*clhs37*clhs38;
const double clhs352 =             N[3]*rho;
const double clhs353 =             1.0*N[3]*clhs17*clhs18*clhs19;
const double clhs354 =             1.0*clhs111*clhs17*clhs19;
const double clhs355 =             clhs110 + clhs14*clhs352 + clhs16*clhs353 + clhs16*clhs354;
const double clhs356 =             clhs131*n;
const double clhs357 =             clhs258 + clhs352*clhs52 + clhs353*clhs54 + clhs354*clhs54;
const double clhs358 =             clhs264*n;
const double clhs359 =             clhs321 + clhs352*clhs82 + clhs353*clhs85 + clhs354*clhs85;
const double clhs360 =             clhs327*n;
const double clhs361 =             pow(DN(3,0), 2);
const double clhs362 =             pow(N[3], 2);
const double clhs363 =             clhs111*clhs352 + clhs114*clhs353 + clhs114*clhs354 + clhs12*clhs362;
const double clhs364 =             DN(3,0)*clhs10;
const double clhs365 =             DN(3,1)*clhs364;
const double clhs366 =             DN(3,2)*clhs364;
const double clhs367 =             clhs111*clhs41 - clhs112 + clhs113*clhs316 + clhs352*clhs40;
const double clhs368 =             clhs186*n;
const double clhs369 =             clhs280*n;
const double clhs370 =             clhs338*n;
const double clhs371 =             pow(DN(3,1), 2);
const double clhs372 =             DN(3,1)*DN(3,2)*clhs10;
const double clhs373 =             clhs220*n;
const double clhs374 =             clhs295*n;
const double clhs375 =             clhs348*n;
const double clhs376 =             pow(DN(3,2), 2);
const double clhs377 =             N[3] + clhs114*clhs41;
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs10*clhs5 + clhs22;
            lhs(0,1)=DN(0,0)*clhs23 + DN(0,1)*clhs25 + DN(0,2)*clhs28 + clhs30;
            lhs(0,2)=DN(0,0)*clhs31 + DN(0,1)*clhs33 + DN(0,2)*clhs35 + clhs36;
            lhs(0,3)=DN(0,0)*clhs43;
            lhs(0,4)=DN(0,0)*clhs44 + DN(0,1)*clhs46 + DN(0,2)*clhs48 + clhs49 + clhs55;
            lhs(0,5)=DN(0,0)*clhs56 + DN(0,1)*clhs58 + DN(0,2)*clhs61 + clhs62;
            lhs(0,6)=DN(0,0)*clhs63 + DN(0,1)*clhs65 + DN(0,2)*clhs67 + clhs68;
            lhs(0,7)=clhs14*clhs74 + clhs70*clhs71 - clhs70 + clhs72*clhs73;
            lhs(0,8)=DN(0,0)*clhs75 + DN(0,1)*clhs77 + DN(0,2)*clhs79 + clhs80 + clhs86;
            lhs(0,9)=DN(0,0)*clhs87 + DN(0,1)*clhs89 + DN(0,2)*clhs92 + clhs93;
            lhs(0,10)=DN(0,0)*clhs94 + DN(0,1)*clhs96 + DN(0,2)*clhs98 + clhs99;
            lhs(0,11)=clhs101*clhs71 - clhs101 + clhs102*clhs73 + clhs103*clhs14;
            lhs(0,12)=DN(0,0)*clhs104 + DN(0,1)*clhs106 + DN(0,2)*clhs108 + clhs109 + clhs115;
            lhs(0,13)=DN(0,0)*clhs116 + DN(0,1)*clhs118 + DN(0,2)*clhs121 + clhs122;
            lhs(0,14)=DN(0,0)*clhs123 + DN(0,1)*clhs125 + DN(0,2)*clhs127 + clhs128;
            lhs(0,15)=DN(3,0)*clhs42 + clhs130*clhs71 - clhs130 + clhs131*clhs73;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs132 + DN(0,2)*clhs133 + clhs30;
            lhs(1,1)=DN(0,0)*clhs25 + DN(0,1)*clhs134 + DN(0,2)*clhs136 + clhs10*clhs137 + clhs22;
            lhs(1,2)=DN(0,0)*clhs33 + DN(0,1)*clhs138 + DN(0,2)*clhs140 + clhs142;
            lhs(1,3)=DN(0,1)*clhs43;
            lhs(1,4)=DN(0,0)*clhs46 + DN(0,1)*clhs143 + DN(0,2)*clhs144 + clhs145;
            lhs(1,5)=DN(0,0)*clhs58 + DN(0,1)*clhs146 + DN(0,2)*clhs148 + clhs149 + clhs55;
            lhs(1,6)=DN(0,0)*clhs65 + DN(0,1)*clhs150 + DN(0,2)*clhs152 + clhs153;
            lhs(1,7)=clhs14*clhs157 + clhs155*clhs71 - clhs155 + clhs156*clhs73;
            lhs(1,8)=DN(0,0)*clhs77 + DN(0,1)*clhs158 + DN(0,2)*clhs159 + clhs160;
            lhs(1,9)=DN(0,0)*clhs89 + DN(0,1)*clhs161 + DN(0,2)*clhs163 + clhs164 + clhs86;
            lhs(1,10)=DN(0,0)*clhs96 + DN(0,1)*clhs165 + DN(0,2)*clhs167 + clhs168;
            lhs(1,11)=clhs14*clhs172 + clhs170*clhs71 - clhs170 + clhs171*clhs73;
            lhs(1,12)=DN(0,0)*clhs106 + DN(0,1)*clhs173 + DN(0,2)*clhs174 + clhs175;
            lhs(1,13)=DN(0,0)*clhs118 + DN(0,1)*clhs176 + DN(0,2)*clhs178 + clhs115 + clhs179;
            lhs(1,14)=DN(0,0)*clhs125 + DN(0,1)*clhs180 + DN(0,2)*clhs182 + clhs183;
            lhs(1,15)=DN(3,1)*clhs42 + clhs185*clhs71 - clhs185 + clhs186*clhs73;
            lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs133 + DN(0,2)*clhs187 + clhs36;
            lhs(2,1)=DN(0,0)*clhs28 + DN(0,1)*clhs136 + DN(0,2)*clhs188 + clhs142;
            lhs(2,2)=DN(0,0)*clhs35 + DN(0,1)*clhs140 + DN(0,2)*clhs189 + clhs10*clhs190 + clhs22;
            lhs(2,3)=DN(0,2)*clhs43;
            lhs(2,4)=DN(0,0)*clhs48 + DN(0,1)*clhs144 + DN(0,2)*clhs191 + clhs193;
            lhs(2,5)=DN(0,0)*clhs61 + DN(0,1)*clhs148 + DN(0,2)*clhs194 + clhs195;
            lhs(2,6)=DN(0,0)*clhs67 + DN(0,1)*clhs152 + DN(0,2)*clhs196 + clhs197 + clhs55;
            lhs(2,7)=clhs14*clhs201 + clhs199*clhs71 - clhs199 + clhs200*clhs73;
            lhs(2,8)=DN(0,0)*clhs79 + DN(0,1)*clhs159 + DN(0,2)*clhs202 + clhs203;
            lhs(2,9)=DN(0,0)*clhs92 + DN(0,1)*clhs163 + DN(0,2)*clhs204 + clhs205;
            lhs(2,10)=DN(0,0)*clhs98 + DN(0,1)*clhs167 + DN(0,2)*clhs206 + clhs207 + clhs86;
            lhs(2,11)=clhs14*clhs211 + clhs209*clhs71 - clhs209 + clhs210*clhs73;
            lhs(2,12)=DN(0,0)*clhs108 + DN(0,1)*clhs174 + DN(0,2)*clhs212 + clhs213;
            lhs(2,13)=DN(0,0)*clhs121 + DN(0,1)*clhs178 + DN(0,2)*clhs214 + clhs215;
            lhs(2,14)=DN(0,0)*clhs127 + DN(0,1)*clhs182 + DN(0,2)*clhs216 + clhs115 + clhs217;
            lhs(2,15)=DN(3,2)*clhs42 + clhs219*clhs71 - clhs219 + clhs220*clhs73;
            lhs(3,0)=DN(0,0)*clhs222;
            lhs(3,1)=DN(0,1)*clhs222;
            lhs(3,2)=DN(0,2)*clhs222;
            lhs(3,3)=clhs11*clhs223 + clhs137*clhs224 + clhs190*clhs224 + clhs224*clhs5;
            lhs(3,4)=clhs225*clhs54 + clhs72;
            lhs(3,5)=clhs156 + clhs226*clhs54;
            lhs(3,6)=clhs200 + clhs227*clhs54;
            lhs(3,7)=clhs231;
            lhs(3,8)=clhs102 + clhs225*clhs85;
            lhs(3,9)=clhs171 + clhs226*clhs85;
            lhs(3,10)=clhs210 + clhs227*clhs85;
            lhs(3,11)=clhs232;
            lhs(3,12)=clhs114*clhs225 + clhs131;
            lhs(3,13)=clhs114*clhs226 + clhs186;
            lhs(3,14)=clhs114*clhs227 + clhs220;
            lhs(3,15)=clhs233;
            lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs237 + clhs49;
            lhs(4,1)=DN(1,0)*clhs23 + DN(1,1)*clhs25 + DN(1,2)*clhs28 + clhs145;
            lhs(4,2)=DN(1,0)*clhs31 + DN(1,1)*clhs33 + DN(1,2)*clhs35 + clhs193;
            lhs(4,3)=clhs225*clhs52 + clhs238*clhs71 - clhs238 + clhs69*clhs73;
            lhs(4,4)=DN(1,0)*clhs44 + DN(1,1)*clhs46 + DN(1,2)*clhs48 + clhs10*clhs239 + clhs241;
            lhs(4,5)=DN(1,0)*clhs56 + DN(1,1)*clhs58 + DN(1,2)*clhs61 + clhs243;
            lhs(4,6)=DN(1,0)*clhs63 + DN(1,1)*clhs65 + DN(1,2)*clhs67 + clhs244;
            lhs(4,7)=DN(1,0)*clhs247;
            lhs(4,8)=DN(1,0)*clhs75 + DN(1,1)*clhs77 + DN(1,2)*clhs79 + clhs248 + clhs251;
            lhs(4,9)=DN(1,0)*clhs87 + DN(1,1)*clhs89 + DN(1,2)*clhs92 + clhs252;
            lhs(4,10)=DN(1,0)*clhs94 + DN(1,1)*clhs96 + DN(1,2)*clhs98 + clhs253;
            lhs(4,11)=clhs103*clhs52 + clhs255*clhs71 - clhs255 + clhs256*clhs73;
            lhs(4,12)=DN(1,0)*clhs104 + DN(1,1)*clhs106 + DN(1,2)*clhs108 + clhs257 + clhs259;
            lhs(4,13)=DN(1,0)*clhs116 + DN(1,1)*clhs118 + DN(1,2)*clhs121 + clhs260;
            lhs(4,14)=DN(1,0)*clhs123 + DN(1,1)*clhs125 + DN(1,2)*clhs127 + clhs261;
            lhs(4,15)=DN(3,0)*clhs246 + clhs263*clhs71 - clhs263 + clhs264*clhs73;
            lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs132 + DN(1,2)*clhs133 + clhs62;
            lhs(5,1)=DN(1,0)*clhs25 + DN(1,1)*clhs134 + DN(1,2)*clhs136 + clhs149 + clhs237;
            lhs(5,2)=DN(1,0)*clhs33 + DN(1,1)*clhs138 + DN(1,2)*clhs140 + clhs195;
            lhs(5,3)=clhs154*clhs73 + clhs226*clhs52 + clhs265*clhs71 - clhs265;
            lhs(5,4)=DN(1,0)*clhs46 + DN(1,1)*clhs143 + DN(1,2)*clhs144 + clhs243;
            lhs(5,5)=DN(1,0)*clhs58 + DN(1,1)*clhs146 + DN(1,2)*clhs148 + clhs10*clhs266 + clhs241;
            lhs(5,6)=DN(1,0)*clhs65 + DN(1,1)*clhs150 + DN(1,2)*clhs152 + clhs268;
            lhs(5,7)=DN(1,1)*clhs247;
            lhs(5,8)=DN(1,0)*clhs77 + DN(1,1)*clhs158 + DN(1,2)*clhs159 + clhs269;
            lhs(5,9)=DN(1,0)*clhs89 + DN(1,1)*clhs161 + DN(1,2)*clhs163 + clhs251 + clhs270;
            lhs(5,10)=DN(1,0)*clhs96 + DN(1,1)*clhs165 + DN(1,2)*clhs167 + clhs271;
            lhs(5,11)=clhs172*clhs52 + clhs273*clhs71 - clhs273 + clhs274*clhs73;
            lhs(5,12)=DN(1,0)*clhs106 + DN(1,1)*clhs173 + DN(1,2)*clhs174 + clhs275;
            lhs(5,13)=DN(1,0)*clhs118 + DN(1,1)*clhs176 + DN(1,2)*clhs178 + clhs259 + clhs276;
            lhs(5,14)=DN(1,0)*clhs125 + DN(1,1)*clhs180 + DN(1,2)*clhs182 + clhs277;
            lhs(5,15)=DN(3,1)*clhs246 + clhs279*clhs71 - clhs279 + clhs280*clhs73;
            lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs133 + DN(1,2)*clhs187 + clhs68;
            lhs(6,1)=DN(1,0)*clhs28 + DN(1,1)*clhs136 + DN(1,2)*clhs188 + clhs153;
            lhs(6,2)=DN(1,0)*clhs35 + DN(1,1)*clhs140 + DN(1,2)*clhs189 + clhs197 + clhs237;
            lhs(6,3)=clhs198*clhs73 + clhs227*clhs52 + clhs281*clhs71 - clhs281;
            lhs(6,4)=DN(1,0)*clhs48 + DN(1,1)*clhs144 + DN(1,2)*clhs191 + clhs244;
            lhs(6,5)=DN(1,0)*clhs61 + DN(1,1)*clhs148 + DN(1,2)*clhs194 + clhs268;
            lhs(6,6)=DN(1,0)*clhs67 + DN(1,1)*clhs152 + DN(1,2)*clhs196 + clhs10*clhs282 + clhs241;
            lhs(6,7)=DN(1,2)*clhs247;
            lhs(6,8)=DN(1,0)*clhs79 + DN(1,1)*clhs159 + DN(1,2)*clhs202 + clhs284;
            lhs(6,9)=DN(1,0)*clhs92 + DN(1,1)*clhs163 + DN(1,2)*clhs204 + clhs285;
            lhs(6,10)=DN(1,0)*clhs98 + DN(1,1)*clhs167 + DN(1,2)*clhs206 + clhs251 + clhs286;
            lhs(6,11)=clhs211*clhs52 + clhs288*clhs71 - clhs288 + clhs289*clhs73;
            lhs(6,12)=DN(1,0)*clhs108 + DN(1,1)*clhs174 + DN(1,2)*clhs212 + clhs290;
            lhs(6,13)=DN(1,0)*clhs121 + DN(1,1)*clhs178 + DN(1,2)*clhs214 + clhs291;
            lhs(6,14)=DN(1,0)*clhs127 + DN(1,1)*clhs182 + DN(1,2)*clhs216 + clhs259 + clhs292;
            lhs(6,15)=DN(3,2)*clhs246 + clhs294*clhs71 - clhs294 + clhs295*clhs73;
            lhs(7,0)=clhs16*clhs74 + clhs69;
            lhs(7,1)=clhs154 + clhs157*clhs16;
            lhs(7,2)=clhs16*clhs201 + clhs198;
            lhs(7,3)=clhs231;
            lhs(7,4)=DN(1,0)*clhs297;
            lhs(7,5)=DN(1,1)*clhs297;
            lhs(7,6)=DN(1,2)*clhs297;
            lhs(7,7)=clhs223*clhs240 + clhs224*clhs239 + clhs224*clhs266 + clhs224*clhs282;
            lhs(7,8)=clhs256 + clhs74*clhs85;
            lhs(7,9)=clhs157*clhs85 + clhs274;
            lhs(7,10)=clhs201*clhs85 + clhs289;
            lhs(7,11)=clhs301;
            lhs(7,12)=clhs114*clhs74 + clhs264;
            lhs(7,13)=clhs114*clhs157 + clhs280;
            lhs(7,14)=clhs114*clhs201 + clhs295;
            lhs(7,15)=clhs302;
            lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs306 + clhs80;
            lhs(8,1)=DN(2,0)*clhs23 + DN(2,1)*clhs25 + DN(2,2)*clhs28 + clhs160;
            lhs(8,2)=DN(2,0)*clhs31 + DN(2,1)*clhs33 + DN(2,2)*clhs35 + clhs203;
            lhs(8,3)=clhs100*clhs73 + clhs225*clhs82 + clhs307*clhs71 - clhs307;
            lhs(8,4)=DN(2,0)*clhs44 + DN(2,1)*clhs46 + DN(2,2)*clhs48 + clhs248 + clhs308;
            lhs(8,5)=DN(2,0)*clhs56 + DN(2,1)*clhs58 + DN(2,2)*clhs61 + clhs269;
            lhs(8,6)=DN(2,0)*clhs63 + DN(2,1)*clhs65 + DN(2,2)*clhs67 + clhs284;
            lhs(8,7)=clhs254*clhs73 + clhs309*clhs71 - clhs309 + clhs74*clhs82;
            lhs(8,8)=DN(2,0)*clhs75 + DN(2,1)*clhs77 + DN(2,2)*clhs79 + clhs10*clhs310 + clhs312;
            lhs(8,9)=DN(2,0)*clhs87 + DN(2,1)*clhs89 + DN(2,2)*clhs92 + clhs314;
            lhs(8,10)=DN(2,0)*clhs94 + DN(2,1)*clhs96 + DN(2,2)*clhs98 + clhs315;
            lhs(8,11)=DN(2,0)*clhs318;
            lhs(8,12)=DN(2,0)*clhs104 + DN(2,1)*clhs106 + DN(2,2)*clhs108 + clhs319 + clhs322;
            lhs(8,13)=DN(2,0)*clhs116 + DN(2,1)*clhs118 + DN(2,2)*clhs121 + clhs323;
            lhs(8,14)=DN(2,0)*clhs123 + DN(2,1)*clhs125 + DN(2,2)*clhs127 + clhs324;
            lhs(8,15)=DN(3,0)*clhs317 + clhs326*clhs71 - clhs326 + clhs327*clhs73;
            lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs132 + DN(2,2)*clhs133 + clhs93;
            lhs(9,1)=DN(2,0)*clhs25 + DN(2,1)*clhs134 + DN(2,2)*clhs136 + clhs164 + clhs306;
            lhs(9,2)=DN(2,0)*clhs33 + DN(2,1)*clhs138 + DN(2,2)*clhs140 + clhs205;
            lhs(9,3)=clhs169*clhs73 + clhs226*clhs82 + clhs328*clhs71 - clhs328;
            lhs(9,4)=DN(2,0)*clhs46 + DN(2,1)*clhs143 + DN(2,2)*clhs144 + clhs252;
            lhs(9,5)=DN(2,0)*clhs58 + DN(2,1)*clhs146 + DN(2,2)*clhs148 + clhs270 + clhs308;
            lhs(9,6)=DN(2,0)*clhs65 + DN(2,1)*clhs150 + DN(2,2)*clhs152 + clhs285;
            lhs(9,7)=clhs157*clhs82 + clhs272*clhs73 + clhs329*clhs71 - clhs329;
            lhs(9,8)=DN(2,0)*clhs77 + DN(2,1)*clhs158 + DN(2,2)*clhs159 + clhs314;
            lhs(9,9)=DN(2,0)*clhs89 + DN(2,1)*clhs161 + DN(2,2)*clhs163 + clhs10*clhs330 + clhs312;
            lhs(9,10)=DN(2,0)*clhs96 + DN(2,1)*clhs165 + DN(2,2)*clhs167 + clhs332;
            lhs(9,11)=DN(2,1)*clhs318;
            lhs(9,12)=DN(2,0)*clhs106 + DN(2,1)*clhs173 + DN(2,2)*clhs174 + clhs333;
            lhs(9,13)=DN(2,0)*clhs118 + DN(2,1)*clhs176 + DN(2,2)*clhs178 + clhs322 + clhs334;
            lhs(9,14)=DN(2,0)*clhs125 + DN(2,1)*clhs180 + DN(2,2)*clhs182 + clhs335;
            lhs(9,15)=DN(3,1)*clhs317 + clhs337*clhs71 - clhs337 + clhs338*clhs73;
            lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs133 + DN(2,2)*clhs187 + clhs99;
            lhs(10,1)=DN(2,0)*clhs28 + DN(2,1)*clhs136 + DN(2,2)*clhs188 + clhs168;
            lhs(10,2)=DN(2,0)*clhs35 + DN(2,1)*clhs140 + DN(2,2)*clhs189 + clhs207 + clhs306;
            lhs(10,3)=clhs208*clhs73 + clhs227*clhs82 + clhs339*clhs71 - clhs339;
            lhs(10,4)=DN(2,0)*clhs48 + DN(2,1)*clhs144 + DN(2,2)*clhs191 + clhs253;
            lhs(10,5)=DN(2,0)*clhs61 + DN(2,1)*clhs148 + DN(2,2)*clhs194 + clhs271;
            lhs(10,6)=DN(2,0)*clhs67 + DN(2,1)*clhs152 + DN(2,2)*clhs196 + clhs286 + clhs308;
            lhs(10,7)=clhs201*clhs82 + clhs287*clhs73 + clhs340*clhs71 - clhs340;
            lhs(10,8)=DN(2,0)*clhs79 + DN(2,1)*clhs159 + DN(2,2)*clhs202 + clhs315;
            lhs(10,9)=DN(2,0)*clhs92 + DN(2,1)*clhs163 + DN(2,2)*clhs204 + clhs332;
            lhs(10,10)=DN(2,0)*clhs98 + DN(2,1)*clhs167 + DN(2,2)*clhs206 + clhs10*clhs341 + clhs312;
            lhs(10,11)=DN(2,2)*clhs318;
            lhs(10,12)=DN(2,0)*clhs108 + DN(2,1)*clhs174 + DN(2,2)*clhs212 + clhs343;
            lhs(10,13)=DN(2,0)*clhs121 + DN(2,1)*clhs178 + DN(2,2)*clhs214 + clhs344;
            lhs(10,14)=DN(2,0)*clhs127 + DN(2,1)*clhs182 + DN(2,2)*clhs216 + clhs322 + clhs345;
            lhs(10,15)=DN(3,2)*clhs317 + clhs347*clhs71 - clhs347 + clhs348*clhs73;
            lhs(11,0)=clhs100 + clhs103*clhs16;
            lhs(11,1)=clhs16*clhs172 + clhs169;
            lhs(11,2)=clhs16*clhs211 + clhs208;
            lhs(11,3)=clhs232;
            lhs(11,4)=clhs103*clhs54 + clhs254;
            lhs(11,5)=clhs172*clhs54 + clhs272;
            lhs(11,6)=clhs211*clhs54 + clhs287;
            lhs(11,7)=clhs301;
            lhs(11,8)=DN(2,0)*clhs350;
            lhs(11,9)=DN(2,1)*clhs350;
            lhs(11,10)=DN(2,2)*clhs350;
            lhs(11,11)=clhs223*clhs311 + clhs224*clhs310 + clhs224*clhs330 + clhs224*clhs341;
            lhs(11,12)=clhs103*clhs114 + clhs327;
            lhs(11,13)=clhs114*clhs172 + clhs338;
            lhs(11,14)=clhs114*clhs211 + clhs348;
            lhs(11,15)=clhs351;
            lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs109 + clhs355;
            lhs(12,1)=DN(3,0)*clhs23 + DN(3,1)*clhs25 + DN(3,2)*clhs28 + clhs175;
            lhs(12,2)=DN(3,0)*clhs31 + DN(3,1)*clhs33 + DN(3,2)*clhs35 + clhs213;
            lhs(12,3)=clhs111*clhs225 + clhs129*clhs73 + clhs356*clhs71 - clhs356;
            lhs(12,4)=DN(3,0)*clhs44 + DN(3,1)*clhs46 + DN(3,2)*clhs48 + clhs257 + clhs357;
            lhs(12,5)=DN(3,0)*clhs56 + DN(3,1)*clhs58 + DN(3,2)*clhs61 + clhs275;
            lhs(12,6)=DN(3,0)*clhs63 + DN(3,1)*clhs65 + DN(3,2)*clhs67 + clhs290;
            lhs(12,7)=clhs111*clhs74 + clhs262*clhs73 + clhs358*clhs71 - clhs358;
            lhs(12,8)=DN(3,0)*clhs75 + DN(3,1)*clhs77 + DN(3,2)*clhs79 + clhs319 + clhs359;
            lhs(12,9)=DN(3,0)*clhs87 + DN(3,1)*clhs89 + DN(3,2)*clhs92 + clhs333;
            lhs(12,10)=DN(3,0)*clhs94 + DN(3,1)*clhs96 + DN(3,2)*clhs98 + clhs343;
            lhs(12,11)=clhs103*clhs111 + clhs325*clhs73 + clhs360*clhs71 - clhs360;
            lhs(12,12)=DN(3,0)*clhs104 + DN(3,1)*clhs106 + DN(3,2)*clhs108 + clhs10*clhs361 + clhs363;
            lhs(12,13)=DN(3,0)*clhs116 + DN(3,1)*clhs118 + DN(3,2)*clhs121 + clhs365;
            lhs(12,14)=DN(3,0)*clhs123 + DN(3,1)*clhs125 + DN(3,2)*clhs127 + clhs366;
            lhs(12,15)=DN(3,0)*clhs367;
            lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs132 + DN(3,2)*clhs133 + clhs122;
            lhs(13,1)=DN(3,0)*clhs25 + DN(3,1)*clhs134 + DN(3,2)*clhs136 + clhs179 + clhs355;
            lhs(13,2)=DN(3,0)*clhs33 + DN(3,1)*clhs138 + DN(3,2)*clhs140 + clhs215;
            lhs(13,3)=clhs111*clhs226 + clhs184*clhs73 + clhs368*clhs71 - clhs368;
            lhs(13,4)=DN(3,0)*clhs46 + DN(3,1)*clhs143 + DN(3,2)*clhs144 + clhs260;
            lhs(13,5)=DN(3,0)*clhs58 + DN(3,1)*clhs146 + DN(3,2)*clhs148 + clhs276 + clhs357;
            lhs(13,6)=DN(3,0)*clhs65 + DN(3,1)*clhs150 + DN(3,2)*clhs152 + clhs291;
            lhs(13,7)=clhs111*clhs157 + clhs278*clhs73 + clhs369*clhs71 - clhs369;
            lhs(13,8)=DN(3,0)*clhs77 + DN(3,1)*clhs158 + DN(3,2)*clhs159 + clhs323;
            lhs(13,9)=DN(3,0)*clhs89 + DN(3,1)*clhs161 + DN(3,2)*clhs163 + clhs334 + clhs359;
            lhs(13,10)=DN(3,0)*clhs96 + DN(3,1)*clhs165 + DN(3,2)*clhs167 + clhs344;
            lhs(13,11)=clhs111*clhs172 + clhs336*clhs73 + clhs370*clhs71 - clhs370;
            lhs(13,12)=DN(3,0)*clhs106 + DN(3,1)*clhs173 + DN(3,2)*clhs174 + clhs365;
            lhs(13,13)=DN(3,0)*clhs118 + DN(3,1)*clhs176 + DN(3,2)*clhs178 + clhs10*clhs371 + clhs363;
            lhs(13,14)=DN(3,0)*clhs125 + DN(3,1)*clhs180 + DN(3,2)*clhs182 + clhs372;
            lhs(13,15)=DN(3,1)*clhs367;
            lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs133 + DN(3,2)*clhs187 + clhs128;
            lhs(14,1)=DN(3,0)*clhs28 + DN(3,1)*clhs136 + DN(3,2)*clhs188 + clhs183;
            lhs(14,2)=DN(3,0)*clhs35 + DN(3,1)*clhs140 + DN(3,2)*clhs189 + clhs217 + clhs355;
            lhs(14,3)=clhs111*clhs227 + clhs218*clhs73 + clhs373*clhs71 - clhs373;
            lhs(14,4)=DN(3,0)*clhs48 + DN(3,1)*clhs144 + DN(3,2)*clhs191 + clhs261;
            lhs(14,5)=DN(3,0)*clhs61 + DN(3,1)*clhs148 + DN(3,2)*clhs194 + clhs277;
            lhs(14,6)=DN(3,0)*clhs67 + DN(3,1)*clhs152 + DN(3,2)*clhs196 + clhs292 + clhs357;
            lhs(14,7)=clhs111*clhs201 + clhs293*clhs73 + clhs374*clhs71 - clhs374;
            lhs(14,8)=DN(3,0)*clhs79 + DN(3,1)*clhs159 + DN(3,2)*clhs202 + clhs324;
            lhs(14,9)=DN(3,0)*clhs92 + DN(3,1)*clhs163 + DN(3,2)*clhs204 + clhs335;
            lhs(14,10)=DN(3,0)*clhs98 + DN(3,1)*clhs167 + DN(3,2)*clhs206 + clhs345 + clhs359;
            lhs(14,11)=clhs111*clhs211 + clhs346*clhs73 + clhs375*clhs71 - clhs375;
            lhs(14,12)=DN(3,0)*clhs108 + DN(3,1)*clhs174 + DN(3,2)*clhs212 + clhs366;
            lhs(14,13)=DN(3,0)*clhs121 + DN(3,1)*clhs178 + DN(3,2)*clhs214 + clhs372;
            lhs(14,14)=DN(3,0)*clhs127 + DN(3,1)*clhs182 + DN(3,2)*clhs216 + clhs10*clhs376 + clhs363;
            lhs(14,15)=DN(3,2)*clhs367;
            lhs(15,0)=DN(3,0)*clhs221 + clhs129;
            lhs(15,1)=DN(3,1)*clhs221 + clhs184;
            lhs(15,2)=DN(3,2)*clhs221 + clhs218;
            lhs(15,3)=clhs233;
            lhs(15,4)=DN(3,0)*clhs296 + clhs262;
            lhs(15,5)=DN(3,1)*clhs296 + clhs278;
            lhs(15,6)=DN(3,2)*clhs296 + clhs293;
            lhs(15,7)=clhs302;
            lhs(15,8)=DN(3,0)*clhs349 + clhs325;
            lhs(15,9)=DN(3,1)*clhs349 + clhs336;
            lhs(15,10)=DN(3,2)*clhs349 + clhs346;
            lhs(15,11)=clhs351;
            lhs(15,12)=DN(3,0)*clhs377;
            lhs(15,13)=DN(3,1)*clhs377;
            lhs(15,14)=DN(3,2)*clhs377;
            lhs(15,15)=clhs223*clhs362 + clhs224*clhs361 + clhs224*clhs371 + clhs224*clhs376;


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

    const double& dt = data.dt;
    const double& n = data.n;
    const double& bdf0 = data.bdf0;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v_ave = data.v_ave;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& vconv = v_ave - vmesh;

    // Get constitutive matrix
    const Matrix& C = data.C;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Stabilization parameters
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
const double clhs8 =             pow(N[0], 2);
const double clhs9 =             bdf0*n*rho;
const double clhs10 =             N[0]*rho;
const double clhs11 =             DN(0,0)*clhs4 + DN(0,1)*clhs5;
const double clhs12 =             N[0]*n;
const double clhs13 =             bdf0*clhs12 + clhs11;
const double clhs14 =             pow(rho, 2);
const double clhs15 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double clhs16 =             1.0/(clhs6/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs17 =             1.0*N[0]*clhs14*clhs15*clhs16;
const double clhs18 =             1.0*clhs11*clhs14*clhs16;
const double clhs19 =             clhs10*clhs11 + clhs13*clhs17 + clhs13*clhs18 + clhs8*clhs9;
const double clhs20 =             C(0,1)*DN(0,1) + clhs1;
const double clhs21 =             C(1,2)*DN(0,1);
const double clhs22 =             C(2,2)*DN(0,0) + clhs21;
const double clhs23 =             DN(0,0)*clhs7;
const double clhs24 =             DN(0,1)*clhs23;
const double clhs25 =             pow(c, -2);
const double clhs26 =             1.0/rho;
const double clhs27 =             N[0]*bdf0*clhs25*clhs26*n;
const double clhs28 =             1.0*clhs15*clhs16;
const double clhs29 =             1.0*clhs16*rho;
const double clhs30 =             clhs11*clhs29;
const double clhs31 =             clhs10*clhs28 - clhs12 + clhs27*clhs7 + clhs30;
const double clhs32 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs33 =             C(0,2)*DN(1,0);
const double clhs34 =             C(2,2)*DN(1,1) + clhs33;
const double clhs35 =             DN(1,0)*clhs23;
const double clhs36 =             N[0]*bdf0*n*rho;
const double clhs37 =             N[1]*clhs36;
const double clhs38 =             DN(1,0)*clhs4 + DN(1,1)*clhs5;
const double clhs39 =             N[1]*n;
const double clhs40 =             bdf0*clhs39;
const double clhs41 =             clhs38 + clhs40;
const double clhs42 =             clhs10*clhs38 + clhs17*clhs41 + clhs18*clhs41 + clhs37;
const double clhs43 =             C(0,1)*DN(1,1) + clhs33;
const double clhs44 =             C(1,2)*DN(1,1);
const double clhs45 =             C(2,2)*DN(1,0) + clhs44;
const double clhs46 =             DN(1,1)*clhs23;
const double clhs47 =             DN(0,0)*N[1];
const double clhs48 =             clhs47*n;
const double clhs49 =             bdf0*clhs25*clhs26*clhs7;
const double clhs50 =             DN(1,0)*N[0];
const double clhs51 =             1.0*clhs15*clhs16*rho;
const double clhs52 =             1.0*DN(1,0)*clhs16*rho;
const double clhs53 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs54 =             C(0,2)*DN(2,0);
const double clhs55 =             C(2,2)*DN(2,1) + clhs54;
const double clhs56 =             DN(2,0)*clhs23;
const double clhs57 =             N[2]*clhs36;
const double clhs58 =             DN(2,0)*clhs4 + DN(2,1)*clhs5;
const double clhs59 =             N[2]*n;
const double clhs60 =             bdf0*clhs59;
const double clhs61 =             clhs58 + clhs60;
const double clhs62 =             clhs10*clhs58 + clhs17*clhs61 + clhs18*clhs61 + clhs57;
const double clhs63 =             C(0,1)*DN(2,1) + clhs54;
const double clhs64 =             C(1,2)*DN(2,1);
const double clhs65 =             C(2,2)*DN(2,0) + clhs64;
const double clhs66 =             DN(2,1)*clhs23;
const double clhs67 =             DN(0,0)*N[2];
const double clhs68 =             clhs67*n;
const double clhs69 =             DN(2,0)*N[0];
const double clhs70 =             C(0,1)*DN(0,0) + clhs21;
const double clhs71 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs72 =             pow(DN(0,1), 2);
const double clhs73 =             C(0,1)*DN(1,0) + clhs44;
const double clhs74 =             DN(0,1)*clhs7;
const double clhs75 =             DN(1,0)*clhs74;
const double clhs76 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs77 =             DN(1,1)*clhs74;
const double clhs78 =             DN(0,1)*N[1];
const double clhs79 =             clhs78*n;
const double clhs80 =             DN(1,1)*N[0];
const double clhs81 =             1.0*DN(1,1)*clhs16*rho;
const double clhs82 =             C(0,1)*DN(2,0) + clhs64;
const double clhs83 =             DN(2,0)*clhs74;
const double clhs84 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs85 =             DN(2,1)*clhs74;
const double clhs86 =             DN(0,1)*N[2];
const double clhs87 =             clhs86*n;
const double clhs88 =             DN(2,1)*N[0];
const double clhs89 =             clhs13*clhs29;
const double clhs90 =             N[0] + clhs89;
const double clhs91 =             bdf0*clhs25*clhs26*n;
const double clhs92 =             1.0*clhs16;
const double clhs93 =             1.0*DN(0,0)*clhs16*rho;
const double clhs94 =             1.0*DN(0,1)*clhs16*rho;
const double clhs95 =             1.0*DN(0,0)*clhs16;
const double clhs96 =             1.0*DN(0,1)*clhs16;
const double clhs97 =             DN(1,0)*clhs95 + DN(1,1)*clhs96 + N[1]*clhs27;
const double clhs98 =             DN(2,0)*clhs95 + DN(2,1)*clhs96 + N[2]*clhs27;
const double clhs99 =             N[1]*rho;
const double clhs100 =             1.0*N[1]*clhs14*clhs15*clhs16;
const double clhs101 =             1.0*clhs14*clhs16*clhs38;
const double clhs102 =             clhs100*clhs13 + clhs101*clhs13 + clhs11*clhs99 + clhs37;
const double clhs103 =             clhs50*n;
const double clhs104 =             pow(DN(1,0), 2);
const double clhs105 =             pow(N[1], 2);
const double clhs106 =             clhs100*clhs41 + clhs101*clhs41 + clhs105*clhs9 + clhs38*clhs99;
const double clhs107 =             DN(1,0)*clhs7;
const double clhs108 =             DN(1,1)*clhs107;
const double clhs109 =             clhs25*clhs26*clhs7;
const double clhs110 =             clhs29*clhs38;
const double clhs111 =             clhs109*clhs40 + clhs110 + clhs28*clhs99 - clhs39;
const double clhs112 =             DN(2,0)*clhs107;
const double clhs113 =             N[1]*N[2]*bdf0*n;
const double clhs114 =             clhs113*rho;
const double clhs115 =             clhs100*clhs61 + clhs101*clhs61 + clhs114 + clhs58*clhs99;
const double clhs116 =             DN(2,1)*clhs107;
const double clhs117 =             DN(1,0)*N[2];
const double clhs118 =             clhs117*n;
const double clhs119 =             DN(2,0)*N[1];
const double clhs120 =             clhs80*n;
const double clhs121 =             pow(DN(1,1), 2);
const double clhs122 =             DN(1,1)*clhs7;
const double clhs123 =             DN(2,0)*clhs122;
const double clhs124 =             DN(2,1)*clhs122;
const double clhs125 =             DN(1,1)*N[2];
const double clhs126 =             clhs125*n;
const double clhs127 =             DN(2,1)*N[1];
const double clhs128 =             clhs29*clhs41;
const double clhs129 =             N[1] + clhs128;
const double clhs130 =             1.0*DN(1,0)*DN(2,0)*clhs16 + 1.0*DN(1,1)*DN(2,1)*clhs16 + clhs113*clhs25*clhs26;
const double clhs131 =             N[2]*rho;
const double clhs132 =             1.0*N[2]*clhs14*clhs15*clhs16;
const double clhs133 =             1.0*clhs14*clhs16*clhs58;
const double clhs134 =             clhs11*clhs131 + clhs13*clhs132 + clhs13*clhs133 + clhs57;
const double clhs135 =             clhs69*n;
const double clhs136 =             clhs114 + clhs131*clhs38 + clhs132*clhs41 + clhs133*clhs41;
const double clhs137 =             clhs119*n;
const double clhs138 =             pow(DN(2,0), 2);
const double clhs139 =             pow(N[2], 2);
const double clhs140 =             clhs131*clhs58 + clhs132*clhs61 + clhs133*clhs61 + clhs139*clhs9;
const double clhs141 =             DN(2,0)*DN(2,1)*clhs7;
const double clhs142 =             clhs109*clhs60 + clhs131*clhs28 + clhs29*clhs58 - clhs59;
const double clhs143 =             clhs88*n;
const double clhs144 =             clhs127*n;
const double clhs145 =             pow(DN(2,1), 2);
const double clhs146 =             N[2] + clhs29*clhs61;
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs19 + clhs3*clhs7;
            lhs(0,1)=DN(0,0)*clhs20 + DN(0,1)*clhs22 + clhs24;
            lhs(0,2)=DN(0,0)*clhs31;
            lhs(0,3)=DN(0,0)*clhs32 + DN(0,1)*clhs34 + clhs35 + clhs42;
            lhs(0,4)=DN(0,0)*clhs43 + DN(0,1)*clhs45 + clhs46;
            lhs(0,5)=clhs11*clhs52 + clhs48*clhs49 - clhs48 + clhs50*clhs51;
            lhs(0,6)=DN(0,0)*clhs53 + DN(0,1)*clhs55 + clhs56 + clhs62;
            lhs(0,7)=DN(0,0)*clhs63 + DN(0,1)*clhs65 + clhs66;
            lhs(0,8)=DN(2,0)*clhs30 + clhs49*clhs68 + clhs51*clhs69 - clhs68;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs70 + clhs24;
            lhs(1,1)=DN(0,0)*clhs22 + DN(0,1)*clhs71 + clhs19 + clhs7*clhs72;
            lhs(1,2)=DN(0,1)*clhs31;
            lhs(1,3)=DN(0,0)*clhs34 + DN(0,1)*clhs73 + clhs75;
            lhs(1,4)=DN(0,0)*clhs45 + DN(0,1)*clhs76 + clhs42 + clhs77;
            lhs(1,5)=clhs11*clhs81 + clhs49*clhs79 + clhs51*clhs80 - clhs79;
            lhs(1,6)=DN(0,0)*clhs55 + DN(0,1)*clhs82 + clhs83;
            lhs(1,7)=DN(0,0)*clhs65 + DN(0,1)*clhs84 + clhs62 + clhs85;
            lhs(1,8)=DN(2,1)*clhs30 + clhs49*clhs87 + clhs51*clhs88 - clhs87;
            lhs(2,0)=DN(0,0)*clhs90;
            lhs(2,1)=DN(0,1)*clhs90;
            lhs(2,2)=clhs3*clhs92 + clhs72*clhs92 + clhs8*clhs91;
            lhs(2,3)=clhs41*clhs93 + clhs50;
            lhs(2,4)=clhs41*clhs94 + clhs80;
            lhs(2,5)=clhs97;
            lhs(2,6)=clhs61*clhs93 + clhs69;
            lhs(2,7)=clhs61*clhs94 + clhs88;
            lhs(2,8)=clhs98;
            lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs102 + clhs35;
            lhs(3,1)=DN(1,0)*clhs20 + DN(1,1)*clhs22 + clhs75;
            lhs(3,2)=clhs103*clhs49 - clhs103 + clhs38*clhs93 + clhs47*clhs51;
            lhs(3,3)=DN(1,0)*clhs32 + DN(1,1)*clhs34 + clhs104*clhs7 + clhs106;
            lhs(3,4)=DN(1,0)*clhs43 + DN(1,1)*clhs45 + clhs108;
            lhs(3,5)=DN(1,0)*clhs111;
            lhs(3,6)=DN(1,0)*clhs53 + DN(1,1)*clhs55 + clhs112 + clhs115;
            lhs(3,7)=DN(1,0)*clhs63 + DN(1,1)*clhs65 + clhs116;
            lhs(3,8)=DN(2,0)*clhs110 + clhs118*clhs49 - clhs118 + clhs119*clhs51;
            lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs70 + clhs46;
            lhs(4,1)=DN(1,0)*clhs22 + DN(1,1)*clhs71 + clhs102 + clhs77;
            lhs(4,2)=clhs120*clhs49 - clhs120 + clhs38*clhs94 + clhs51*clhs78;
            lhs(4,3)=DN(1,0)*clhs34 + DN(1,1)*clhs73 + clhs108;
            lhs(4,4)=DN(1,0)*clhs45 + DN(1,1)*clhs76 + clhs106 + clhs121*clhs7;
            lhs(4,5)=DN(1,1)*clhs111;
            lhs(4,6)=DN(1,0)*clhs55 + DN(1,1)*clhs82 + clhs123;
            lhs(4,7)=DN(1,0)*clhs65 + DN(1,1)*clhs84 + clhs115 + clhs124;
            lhs(4,8)=DN(2,1)*clhs110 + clhs126*clhs49 - clhs126 + clhs127*clhs51;
            lhs(5,0)=clhs13*clhs52 + clhs47;
            lhs(5,1)=clhs13*clhs81 + clhs78;
            lhs(5,2)=clhs97;
            lhs(5,3)=DN(1,0)*clhs129;
            lhs(5,4)=DN(1,1)*clhs129;
            lhs(5,5)=clhs104*clhs92 + clhs105*clhs91 + clhs121*clhs92;
            lhs(5,6)=clhs119 + clhs52*clhs61;
            lhs(5,7)=clhs127 + clhs61*clhs81;
            lhs(5,8)=clhs130;
            lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs134 + clhs56;
            lhs(6,1)=DN(2,0)*clhs20 + DN(2,1)*clhs22 + clhs83;
            lhs(6,2)=clhs135*clhs49 - clhs135 + clhs51*clhs67 + clhs58*clhs93;
            lhs(6,3)=DN(2,0)*clhs32 + DN(2,1)*clhs34 + clhs112 + clhs136;
            lhs(6,4)=DN(2,0)*clhs43 + DN(2,1)*clhs45 + clhs123;
            lhs(6,5)=clhs117*clhs51 + clhs137*clhs49 - clhs137 + clhs52*clhs58;
            lhs(6,6)=DN(2,0)*clhs53 + DN(2,1)*clhs55 + clhs138*clhs7 + clhs140;
            lhs(6,7)=DN(2,0)*clhs63 + DN(2,1)*clhs65 + clhs141;
            lhs(6,8)=DN(2,0)*clhs142;
            lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs70 + clhs66;
            lhs(7,1)=DN(2,0)*clhs22 + DN(2,1)*clhs71 + clhs134 + clhs85;
            lhs(7,2)=clhs143*clhs49 - clhs143 + clhs51*clhs86 + clhs58*clhs94;
            lhs(7,3)=DN(2,0)*clhs34 + DN(2,1)*clhs73 + clhs116;
            lhs(7,4)=DN(2,0)*clhs45 + DN(2,1)*clhs76 + clhs124 + clhs136;
            lhs(7,5)=clhs125*clhs51 + clhs144*clhs49 - clhs144 + clhs58*clhs81;
            lhs(7,6)=DN(2,0)*clhs55 + DN(2,1)*clhs82 + clhs141;
            lhs(7,7)=DN(2,0)*clhs65 + DN(2,1)*clhs84 + clhs140 + clhs145*clhs7;
            lhs(7,8)=DN(2,1)*clhs142;
            lhs(8,0)=DN(2,0)*clhs89 + clhs67;
            lhs(8,1)=DN(2,1)*clhs89 + clhs86;
            lhs(8,2)=clhs98;
            lhs(8,3)=DN(2,0)*clhs128 + clhs117;
            lhs(8,4)=DN(2,1)*clhs128 + clhs125;
            lhs(8,5)=clhs130;
            lhs(8,6)=DN(2,0)*clhs146;
            lhs(8,7)=DN(2,1)*clhs146;
            lhs(8,8)=clhs138*clhs92 + clhs139*clhs91 + clhs145*clhs92;


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

    const double& dt = data.dt;
    const double& n = data.n;
    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v_ave = data.v_ave;
    const BoundedMatrix<double,nnodes,dim>& vn_ave = data.vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnn_ave = data.vnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnnn_ave = data.vnnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& vconv = v_ave - vmesh;
    const BoundedMatrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& p_ave = data.p_ave;
    const array_1d<double,nnodes>& pn_ave = data.pn_ave;
    const array_1d<double,nnodes>& pnn_ave = data.pnn_ave;
    const array_1d<double,nnodes>& pnnn_ave = data.pnnn_ave;
    const array_1d<double,strain_size>& stress = data.stress;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,dim> f_gauss = prod(trans(f), N);
    const array_1d<double,dim> grad_p = prod(trans(DN), p_ave);

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    const double crhs0 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crhs1 =             n - 1;
const double crhs2 =             crhs1/n;
const double crhs3 =             -crhs2*pn_ave[0] + p_ave[0];
const double crhs4 =             -crhs2*pn_ave[1] + p_ave[1];
const double crhs5 =             -crhs2*pn_ave[2] + p_ave[2];
const double crhs6 =             -crhs2*pn_ave[3] + p_ave[3];
const double crhs7 =             n*(N[0]*crhs3 + N[1]*crhs4 + N[2]*crhs5 + N[3]*crhs6);
const double crhs8 =             DN(0,0)*v_ave(0,0);
const double crhs9 =             DN(1,0)*v_ave(1,0);
const double crhs10 =             DN(2,0)*v_ave(2,0);
const double crhs11 =             DN(3,0)*v_ave(3,0);
const double crhs12 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs13 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs14 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs15 =             rho*(crhs12*(crhs10 + crhs11 + crhs8 + crhs9) + crhs13*(DN(0,1)*v_ave(0,0) + DN(1,1)*v_ave(1,0) + DN(2,1)*v_ave(2,0) + DN(3,1)*v_ave(3,0)) + crhs14*(DN(0,2)*v_ave(0,0) + DN(1,2)*v_ave(1,0) + DN(2,2)*v_ave(2,0) + DN(3,2)*v_ave(3,0)));
const double crhs16 =             bdf0*n;
const double crhs17 =             bdf1*crhs1;
const double crhs18 =             n - 2;
const double crhs19 =             crhs18/crhs1;
const double crhs20 =             bdf2*crhs18;
const double crhs21 =             (n - 3)/crhs18;
const double crhs22 =             rho*(N[0]*(crhs16*(-crhs2*vn_ave(0,0) + v_ave(0,0)) + crhs17*(-crhs19*vnn_ave(0,0) + vn_ave(0,0)) + crhs20*(-crhs21*vnnn_ave(0,0) + vnn_ave(0,0))) + N[1]*(crhs16*(-crhs2*vn_ave(1,0) + v_ave(1,0)) + crhs17*(-crhs19*vnn_ave(1,0) + vn_ave(1,0)) + crhs20*(-crhs21*vnnn_ave(1,0) + vnn_ave(1,0))) + N[2]*(crhs16*(-crhs2*vn_ave(2,0) + v_ave(2,0)) + crhs17*(-crhs19*vnn_ave(2,0) + vn_ave(2,0)) + crhs20*(-crhs21*vnnn_ave(2,0) + vnn_ave(2,0))) + N[3]*(crhs16*(-crhs2*vn_ave(3,0) + v_ave(3,0)) + crhs17*(-crhs19*vnn_ave(3,0) + vn_ave(3,0)) + crhs20*(-crhs21*vnnn_ave(3,0) + vnn_ave(3,0))));
const double crhs23 =             rho*stab_c2*sqrt(pow(crhs12, 2) + pow(crhs13, 2) + pow(crhs14, 2));
const double crhs24 =             DN(0,2)*v_ave(0,2) + DN(1,2)*v_ave(1,2) + DN(2,2)*v_ave(2,2) + DN(3,2)*v_ave(3,2);
const double crhs25 =             DN(0,1)*v_ave(0,1);
const double crhs26 =             DN(1,1)*v_ave(1,1);
const double crhs27 =             DN(2,1)*v_ave(2,1);
const double crhs28 =             DN(3,1)*v_ave(3,1);
const double crhs29 =             crhs10 + crhs11 + crhs24 + crhs25 + crhs26 + crhs27 + crhs28 + crhs8 + crhs9;
const double crhs30 =             (N[0]*(crhs16*crhs3 + crhs17*(-crhs19*pnn_ave[0] + pn_ave[0]) + crhs20*(-crhs21*pnnn_ave[0] + pnn_ave[0])) + N[1]*(crhs16*crhs4 + crhs17*(-crhs19*pnn_ave[1] + pn_ave[1]) + crhs20*(-crhs21*pnnn_ave[1] + pnn_ave[1])) + N[2]*(crhs16*crhs5 + crhs17*(-crhs19*pnn_ave[2] + pn_ave[2]) + crhs20*(-crhs21*pnnn_ave[2] + pnn_ave[2])) + N[3]*(crhs16*crhs6 + crhs17*(-crhs19*pnn_ave[3] + pn_ave[3]) + crhs20*(-crhs21*pnnn_ave[3] + pnn_ave[3])))/(pow(c, 2)*rho);
const double crhs31 =             (crhs29 + crhs30)*(crhs23*h/stab_c1 + mu);
const double crhs32 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crhs33 =             N[0]*crhs32*rho;
const double crhs34 =             1.0/(crhs23/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs35 =             1.0*crhs34*(DN(0,0)*p_ave[0] + DN(1,0)*p_ave[1] + DN(2,0)*p_ave[2] + DN(3,0)*p_ave[3] - crhs0 + crhs15 + crhs22);
const double crhs36 =             rho*(DN(0,0)*crhs12 + DN(0,1)*crhs13 + DN(0,2)*crhs14);
const double crhs37 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs38 =             rho*(crhs12*(DN(0,0)*v_ave(0,1) + DN(1,0)*v_ave(1,1) + DN(2,0)*v_ave(2,1) + DN(3,0)*v_ave(3,1)) + crhs13*(crhs25 + crhs26 + crhs27 + crhs28) + crhs14*(DN(0,2)*v_ave(0,1) + DN(1,2)*v_ave(1,1) + DN(2,2)*v_ave(2,1) + DN(3,2)*v_ave(3,1)));
const double crhs39 =             rho*(N[0]*(crhs16*(-crhs2*vn_ave(0,1) + v_ave(0,1)) + crhs17*(-crhs19*vnn_ave(0,1) + vn_ave(0,1)) + crhs20*(-crhs21*vnnn_ave(0,1) + vnn_ave(0,1))) + N[1]*(crhs16*(-crhs2*vn_ave(1,1) + v_ave(1,1)) + crhs17*(-crhs19*vnn_ave(1,1) + vn_ave(1,1)) + crhs20*(-crhs21*vnnn_ave(1,1) + vnn_ave(1,1))) + N[2]*(crhs16*(-crhs2*vn_ave(2,1) + v_ave(2,1)) + crhs17*(-crhs19*vnn_ave(2,1) + vn_ave(2,1)) + crhs20*(-crhs21*vnnn_ave(2,1) + vnn_ave(2,1))) + N[3]*(crhs16*(-crhs2*vn_ave(3,1) + v_ave(3,1)) + crhs17*(-crhs19*vnn_ave(3,1) + vn_ave(3,1)) + crhs20*(-crhs21*vnnn_ave(3,1) + vnn_ave(3,1))));
const double crhs40 =             1.0*crhs34*(DN(0,1)*p_ave[0] + DN(1,1)*p_ave[1] + DN(2,1)*p_ave[2] + DN(3,1)*p_ave[3] - crhs37 + crhs38 + crhs39);
const double crhs41 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs42 =             rho*(crhs12*(DN(0,0)*v_ave(0,2) + DN(1,0)*v_ave(1,2) + DN(2,0)*v_ave(2,2) + DN(3,0)*v_ave(3,2)) + crhs13*(DN(0,1)*v_ave(0,2) + DN(1,1)*v_ave(1,2) + DN(2,1)*v_ave(2,2) + DN(3,1)*v_ave(3,2)) + crhs14*crhs24);
const double crhs43 =             rho*(N[0]*(crhs16*(-crhs2*vn_ave(0,2) + v_ave(0,2)) + crhs17*(-crhs19*vnn_ave(0,2) + vn_ave(0,2)) + crhs20*(-crhs21*vnnn_ave(0,2) + vnn_ave(0,2))) + N[1]*(crhs16*(-crhs2*vn_ave(1,2) + v_ave(1,2)) + crhs17*(-crhs19*vnn_ave(1,2) + vn_ave(1,2)) + crhs20*(-crhs21*vnnn_ave(1,2) + vnn_ave(1,2))) + N[2]*(crhs16*(-crhs2*vn_ave(2,2) + v_ave(2,2)) + crhs17*(-crhs19*vnn_ave(2,2) + vn_ave(2,2)) + crhs20*(-crhs21*vnnn_ave(2,2) + vnn_ave(2,2))) + N[3]*(crhs16*(-crhs2*vn_ave(3,2) + v_ave(3,2)) + crhs17*(-crhs19*vnn_ave(3,2) + vn_ave(3,2)) + crhs20*(-crhs21*vnnn_ave(3,2) + vnn_ave(3,2))));
const double crhs44 =             1.0*crhs34*(DN(0,2)*p_ave[0] + DN(1,2)*p_ave[1] + DN(2,2)*p_ave[2] + DN(3,2)*p_ave[3] - crhs41 + crhs42 + crhs43);
const double crhs45 =             N[1]*crhs32*rho;
const double crhs46 =             rho*(DN(1,0)*crhs12 + DN(1,1)*crhs13 + DN(1,2)*crhs14);
const double crhs47 =             N[2]*crhs32*rho;
const double crhs48 =             rho*(DN(2,0)*crhs12 + DN(2,1)*crhs13 + DN(2,2)*crhs14);
const double crhs49 =             N[3]*crhs32*rho;
const double crhs50 =             rho*(DN(3,0)*crhs12 + DN(3,1)*crhs13 + DN(3,2)*crhs14);
            rhs[0]=-DN(0,0)*crhs31 + DN(0,0)*crhs7 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs0 - N[0]*crhs15 - N[0]*crhs22 - crhs33*crhs35 - crhs35*crhs36;
            rhs[1]=-DN(0,0)*stress[3] - DN(0,1)*crhs31 + DN(0,1)*crhs7 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs37 - N[0]*crhs38 - N[0]*crhs39 - crhs33*crhs40 - crhs36*crhs40;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] - DN(0,2)*crhs31 + DN(0,2)*crhs7 - DN(0,2)*stress[2] + N[0]*crhs41 - N[0]*crhs42 - N[0]*crhs43 - crhs33*crhs44 - crhs36*crhs44;
            rhs[3]=-DN(0,0)*crhs35 - DN(0,1)*crhs40 - DN(0,2)*crhs44 - N[0]*crhs29 - N[0]*crhs30;
            rhs[4]=-DN(1,0)*crhs31 + DN(1,0)*crhs7 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs0 - N[1]*crhs15 - N[1]*crhs22 - crhs35*crhs45 - crhs35*crhs46;
            rhs[5]=-DN(1,0)*stress[3] - DN(1,1)*crhs31 + DN(1,1)*crhs7 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs37 - N[1]*crhs38 - N[1]*crhs39 - crhs40*crhs45 - crhs40*crhs46;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] - DN(1,2)*crhs31 + DN(1,2)*crhs7 - DN(1,2)*stress[2] + N[1]*crhs41 - N[1]*crhs42 - N[1]*crhs43 - crhs44*crhs45 - crhs44*crhs46;
            rhs[7]=-DN(1,0)*crhs35 - DN(1,1)*crhs40 - DN(1,2)*crhs44 - N[1]*crhs29 - N[1]*crhs30;
            rhs[8]=-DN(2,0)*crhs31 + DN(2,0)*crhs7 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs0 - N[2]*crhs15 - N[2]*crhs22 - crhs35*crhs47 - crhs35*crhs48;
            rhs[9]=-DN(2,0)*stress[3] - DN(2,1)*crhs31 + DN(2,1)*crhs7 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs37 - N[2]*crhs38 - N[2]*crhs39 - crhs40*crhs47 - crhs40*crhs48;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] - DN(2,2)*crhs31 + DN(2,2)*crhs7 - DN(2,2)*stress[2] + N[2]*crhs41 - N[2]*crhs42 - N[2]*crhs43 - crhs44*crhs47 - crhs44*crhs48;
            rhs[11]=-DN(2,0)*crhs35 - DN(2,1)*crhs40 - DN(2,2)*crhs44 - N[2]*crhs29 - N[2]*crhs30;
            rhs[12]=-DN(3,0)*crhs31 + DN(3,0)*crhs7 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs0 - N[3]*crhs15 - N[3]*crhs22 - crhs35*crhs49 - crhs35*crhs50;
            rhs[13]=-DN(3,0)*stress[3] - DN(3,1)*crhs31 + DN(3,1)*crhs7 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs37 - N[3]*crhs38 - N[3]*crhs39 - crhs40*crhs49 - crhs40*crhs50;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] - DN(3,2)*crhs31 + DN(3,2)*crhs7 - DN(3,2)*stress[2] + N[3]*crhs41 - N[3]*crhs42 - N[3]*crhs43 - crhs44*crhs49 - crhs44*crhs50;
            rhs[15]=-DN(3,0)*crhs35 - DN(3,1)*crhs40 - DN(3,2)*crhs44 - N[3]*crhs29 - N[3]*crhs30;


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

    const double& dt = data.dt;
    const double& n = data.n;
    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v_ave = data.v_ave;
    const BoundedMatrix<double,nnodes,dim>& vn_ave = data.vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnn_ave = data.vnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnnn_ave = data.vnnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& vconv = v_ave - vmesh;
    const BoundedMatrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& p_ave = data.p_ave;
    const array_1d<double,nnodes>& pn_ave = data.pn_ave;
    const array_1d<double,nnodes>& pnn_ave = data.pnn_ave;
    const array_1d<double,nnodes>& pnnn_ave = data.pnnn_ave;
    const array_1d<double,strain_size>& stress = data.stress;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,dim> f_gauss = prod(trans(f), N);
    const array_1d<double,dim> grad_p = prod(trans(DN), p_ave);

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    const double crhs0 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crhs1 =             n - 1;
const double crhs2 =             crhs1/n;
const double crhs3 =             -crhs2*pn_ave[0] + p_ave[0];
const double crhs4 =             -crhs2*pn_ave[1] + p_ave[1];
const double crhs5 =             -crhs2*pn_ave[2] + p_ave[2];
const double crhs6 =             n*(N[0]*crhs3 + N[1]*crhs4 + N[2]*crhs5);
const double crhs7 =             DN(0,0)*v_ave(0,0) + DN(1,0)*v_ave(1,0) + DN(2,0)*v_ave(2,0);
const double crhs8 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs9 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs10 =             rho*(crhs7*crhs8 + crhs9*(DN(0,1)*v_ave(0,0) + DN(1,1)*v_ave(1,0) + DN(2,1)*v_ave(2,0)));
const double crhs11 =             bdf0*n;
const double crhs12 =             bdf1*crhs1;
const double crhs13 =             n - 2;
const double crhs14 =             crhs13/crhs1;
const double crhs15 =             bdf2*crhs13;
const double crhs16 =             (n - 3)/crhs13;
const double crhs17 =             rho*(N[0]*(crhs11*(-crhs2*vn_ave(0,0) + v_ave(0,0)) + crhs12*(-crhs14*vnn_ave(0,0) + vn_ave(0,0)) + crhs15*(-crhs16*vnnn_ave(0,0) + vnn_ave(0,0))) + N[1]*(crhs11*(-crhs2*vn_ave(1,0) + v_ave(1,0)) + crhs12*(-crhs14*vnn_ave(1,0) + vn_ave(1,0)) + crhs15*(-crhs16*vnnn_ave(1,0) + vnn_ave(1,0))) + N[2]*(crhs11*(-crhs2*vn_ave(2,0) + v_ave(2,0)) + crhs12*(-crhs14*vnn_ave(2,0) + vn_ave(2,0)) + crhs15*(-crhs16*vnnn_ave(2,0) + vnn_ave(2,0))));
const double crhs18 =             rho*stab_c2*sqrt(pow(crhs8, 2) + pow(crhs9, 2));
const double crhs19 =             DN(0,1)*v_ave(0,1) + DN(1,1)*v_ave(1,1) + DN(2,1)*v_ave(2,1);
const double crhs20 =             crhs19 + crhs7;
const double crhs21 =             (N[0]*(crhs11*crhs3 + crhs12*(-crhs14*pnn_ave[0] + pn_ave[0]) + crhs15*(-crhs16*pnnn_ave[0] + pnn_ave[0])) + N[1]*(crhs11*crhs4 + crhs12*(-crhs14*pnn_ave[1] + pn_ave[1]) + crhs15*(-crhs16*pnnn_ave[1] + pnn_ave[1])) + N[2]*(crhs11*crhs5 + crhs12*(-crhs14*pnn_ave[2] + pn_ave[2]) + crhs15*(-crhs16*pnnn_ave[2] + pnn_ave[2])))/(pow(c, 2)*rho);
const double crhs22 =             (crhs20 + crhs21)*(crhs18*h/stab_c1 + mu);
const double crhs23 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crhs24 =             N[0]*crhs23*rho;
const double crhs25 =             1.0/(crhs18/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs26 =             1.0*crhs25*(DN(0,0)*p_ave[0] + DN(1,0)*p_ave[1] + DN(2,0)*p_ave[2] - crhs0 + crhs10 + crhs17);
const double crhs27 =             rho*(DN(0,0)*crhs8 + DN(0,1)*crhs9);
const double crhs28 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs29 =             rho*(crhs19*crhs9 + crhs8*(DN(0,0)*v_ave(0,1) + DN(1,0)*v_ave(1,1) + DN(2,0)*v_ave(2,1)));
const double crhs30 =             rho*(N[0]*(crhs11*(-crhs2*vn_ave(0,1) + v_ave(0,1)) + crhs12*(-crhs14*vnn_ave(0,1) + vn_ave(0,1)) + crhs15*(-crhs16*vnnn_ave(0,1) + vnn_ave(0,1))) + N[1]*(crhs11*(-crhs2*vn_ave(1,1) + v_ave(1,1)) + crhs12*(-crhs14*vnn_ave(1,1) + vn_ave(1,1)) + crhs15*(-crhs16*vnnn_ave(1,1) + vnn_ave(1,1))) + N[2]*(crhs11*(-crhs2*vn_ave(2,1) + v_ave(2,1)) + crhs12*(-crhs14*vnn_ave(2,1) + vn_ave(2,1)) + crhs15*(-crhs16*vnnn_ave(2,1) + vnn_ave(2,1))));
const double crhs31 =             1.0*crhs25*(DN(0,1)*p_ave[0] + DN(1,1)*p_ave[1] + DN(2,1)*p_ave[2] - crhs28 + crhs29 + crhs30);
const double crhs32 =             N[1]*crhs23*rho;
const double crhs33 =             rho*(DN(1,0)*crhs8 + DN(1,1)*crhs9);
const double crhs34 =             N[2]*crhs23*rho;
const double crhs35 =             rho*(DN(2,0)*crhs8 + DN(2,1)*crhs9);
            rhs[0]=-DN(0,0)*crhs22 + DN(0,0)*crhs6 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs0 - N[0]*crhs10 - N[0]*crhs17 - crhs24*crhs26 - crhs26*crhs27;
            rhs[1]=-DN(0,0)*stress[2] - DN(0,1)*crhs22 + DN(0,1)*crhs6 - DN(0,1)*stress[1] + N[0]*crhs28 - N[0]*crhs29 - N[0]*crhs30 - crhs24*crhs31 - crhs27*crhs31;
            rhs[2]=-DN(0,0)*crhs26 - DN(0,1)*crhs31 - N[0]*crhs20 - N[0]*crhs21;
            rhs[3]=-DN(1,0)*crhs22 + DN(1,0)*crhs6 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs0 - N[1]*crhs10 - N[1]*crhs17 - crhs26*crhs32 - crhs26*crhs33;
            rhs[4]=-DN(1,0)*stress[2] - DN(1,1)*crhs22 + DN(1,1)*crhs6 - DN(1,1)*stress[1] + N[1]*crhs28 - N[1]*crhs29 - N[1]*crhs30 - crhs31*crhs32 - crhs31*crhs33;
            rhs[5]=-DN(1,0)*crhs26 - DN(1,1)*crhs31 - N[1]*crhs20 - N[1]*crhs21;
            rhs[6]=-DN(2,0)*crhs22 + DN(2,0)*crhs6 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs0 - N[2]*crhs10 - N[2]*crhs17 - crhs26*crhs34 - crhs26*crhs35;
            rhs[7]=-DN(2,0)*stress[2] - DN(2,1)*crhs22 + DN(2,1)*crhs6 - DN(2,1)*stress[1] + N[2]*crhs28 - N[2]*crhs29 - N[2]*crhs30 - crhs31*crhs34 - crhs31*crhs35;
            rhs[8]=-DN(2,0)*crhs26 - DN(2,1)*crhs31 - N[2]*crhs20 - N[2]*crhs21;


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

    const double& dt = data.dt;
    const double& n = data.n;
    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v_ave = data.v_ave;
    const BoundedMatrix<double,nnodes,dim>& vn_ave = data.vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnn_ave = data.vnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnnn_ave = data.vnnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& vconv = v_ave - vmesh;
    const BoundedMatrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& p_ave = data.p_ave;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the error estimator
    array_1d<double,dim> v_s_gauss;
    const array_1d<double,dim> v_gauss = prod(trans(v_ave), N);
    const array_1d<double,dim> f_gauss = prod(trans(f), N);
    const array_1d<double,dim> grad_p = prod(trans(DN), p_ave);

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Gauss point velocity subscale value computation
    const double cv_s_gauss0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cv_s_gauss1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cv_s_gauss2 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cv_s_gauss3 =             1.0/(rho*stab_c2*sqrt(pow(cv_s_gauss0, 2) + pow(cv_s_gauss1, 2) + pow(cv_s_gauss2, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cv_s_gauss4 =             bdf0*n;
const double cv_s_gauss5 =             n - 1;
const double cv_s_gauss6 =             cv_s_gauss5/n;
const double cv_s_gauss7 =             bdf1*cv_s_gauss5;
const double cv_s_gauss8 =             n - 2;
const double cv_s_gauss9 =             cv_s_gauss8/cv_s_gauss5;
const double cv_s_gauss10 =             bdf2*cv_s_gauss8;
const double cv_s_gauss11 =             (n - 3)/cv_s_gauss8;
            v_s_gauss[0]=-cv_s_gauss3*(DN(0,0)*p_ave[0] + DN(1,0)*p_ave[1] + DN(2,0)*p_ave[2] + DN(3,0)*p_ave[3] + rho*(-N[0]*f(0,0) + N[0]*(cv_s_gauss10*(-cv_s_gauss11*vnnn_ave(0,0) + vnn_ave(0,0)) + cv_s_gauss4*(-cv_s_gauss6*vn_ave(0,0) + v_ave(0,0)) + cv_s_gauss7*(-cv_s_gauss9*vnn_ave(0,0) + vn_ave(0,0))) - N[1]*f(1,0) + N[1]*(cv_s_gauss10*(-cv_s_gauss11*vnnn_ave(1,0) + vnn_ave(1,0)) + cv_s_gauss4*(-cv_s_gauss6*vn_ave(1,0) + v_ave(1,0)) + cv_s_gauss7*(-cv_s_gauss9*vnn_ave(1,0) + vn_ave(1,0))) - N[2]*f(2,0) + N[2]*(cv_s_gauss10*(-cv_s_gauss11*vnnn_ave(2,0) + vnn_ave(2,0)) + cv_s_gauss4*(-cv_s_gauss6*vn_ave(2,0) + v_ave(2,0)) + cv_s_gauss7*(-cv_s_gauss9*vnn_ave(2,0) + vn_ave(2,0))) - N[3]*f(3,0) + N[3]*(cv_s_gauss10*(-cv_s_gauss11*vnnn_ave(3,0) + vnn_ave(3,0)) + cv_s_gauss4*(-cv_s_gauss6*vn_ave(3,0) + v_ave(3,0)) + cv_s_gauss7*(-cv_s_gauss9*vnn_ave(3,0) + vn_ave(3,0))) + cv_s_gauss0*(DN(0,0)*v_ave(0,0) + DN(1,0)*v_ave(1,0) + DN(2,0)*v_ave(2,0) + DN(3,0)*v_ave(3,0)) + cv_s_gauss1*(DN(0,1)*v_ave(0,0) + DN(1,1)*v_ave(1,0) + DN(2,1)*v_ave(2,0) + DN(3,1)*v_ave(3,0)) + cv_s_gauss2*(DN(0,2)*v_ave(0,0) + DN(1,2)*v_ave(1,0) + DN(2,2)*v_ave(2,0) + DN(3,2)*v_ave(3,0))));
            v_s_gauss[1]=-cv_s_gauss3*(DN(0,1)*p_ave[0] + DN(1,1)*p_ave[1] + DN(2,1)*p_ave[2] + DN(3,1)*p_ave[3] + rho*(-N[0]*f(0,1) + N[0]*(cv_s_gauss10*(-cv_s_gauss11*vnnn_ave(0,1) + vnn_ave(0,1)) + cv_s_gauss4*(-cv_s_gauss6*vn_ave(0,1) + v_ave(0,1)) + cv_s_gauss7*(-cv_s_gauss9*vnn_ave(0,1) + vn_ave(0,1))) - N[1]*f(1,1) + N[1]*(cv_s_gauss10*(-cv_s_gauss11*vnnn_ave(1,1) + vnn_ave(1,1)) + cv_s_gauss4*(-cv_s_gauss6*vn_ave(1,1) + v_ave(1,1)) + cv_s_gauss7*(-cv_s_gauss9*vnn_ave(1,1) + vn_ave(1,1))) - N[2]*f(2,1) + N[2]*(cv_s_gauss10*(-cv_s_gauss11*vnnn_ave(2,1) + vnn_ave(2,1)) + cv_s_gauss4*(-cv_s_gauss6*vn_ave(2,1) + v_ave(2,1)) + cv_s_gauss7*(-cv_s_gauss9*vnn_ave(2,1) + vn_ave(2,1))) - N[3]*f(3,1) + N[3]*(cv_s_gauss10*(-cv_s_gauss11*vnnn_ave(3,1) + vnn_ave(3,1)) + cv_s_gauss4*(-cv_s_gauss6*vn_ave(3,1) + v_ave(3,1)) + cv_s_gauss7*(-cv_s_gauss9*vnn_ave(3,1) + vn_ave(3,1))) + cv_s_gauss0*(DN(0,0)*v_ave(0,1) + DN(1,0)*v_ave(1,1) + DN(2,0)*v_ave(2,1) + DN(3,0)*v_ave(3,1)) + cv_s_gauss1*(DN(0,1)*v_ave(0,1) + DN(1,1)*v_ave(1,1) + DN(2,1)*v_ave(2,1) + DN(3,1)*v_ave(3,1)) + cv_s_gauss2*(DN(0,2)*v_ave(0,1) + DN(1,2)*v_ave(1,1) + DN(2,2)*v_ave(2,1) + DN(3,2)*v_ave(3,1))));
            v_s_gauss[2]=-cv_s_gauss3*(DN(0,2)*p_ave[0] + DN(1,2)*p_ave[1] + DN(2,2)*p_ave[2] + DN(3,2)*p_ave[3] + rho*(-N[0]*f(0,2) + N[0]*(cv_s_gauss10*(-cv_s_gauss11*vnnn_ave(0,2) + vnn_ave(0,2)) + cv_s_gauss4*(-cv_s_gauss6*vn_ave(0,2) + v_ave(0,2)) + cv_s_gauss7*(-cv_s_gauss9*vnn_ave(0,2) + vn_ave(0,2))) - N[1]*f(1,2) + N[1]*(cv_s_gauss10*(-cv_s_gauss11*vnnn_ave(1,2) + vnn_ave(1,2)) + cv_s_gauss4*(-cv_s_gauss6*vn_ave(1,2) + v_ave(1,2)) + cv_s_gauss7*(-cv_s_gauss9*vnn_ave(1,2) + vn_ave(1,2))) - N[2]*f(2,2) + N[2]*(cv_s_gauss10*(-cv_s_gauss11*vnnn_ave(2,2) + vnn_ave(2,2)) + cv_s_gauss4*(-cv_s_gauss6*vn_ave(2,2) + v_ave(2,2)) + cv_s_gauss7*(-cv_s_gauss9*vnn_ave(2,2) + vn_ave(2,2))) - N[3]*f(3,2) + N[3]*(cv_s_gauss10*(-cv_s_gauss11*vnnn_ave(3,2) + vnn_ave(3,2)) + cv_s_gauss4*(-cv_s_gauss6*vn_ave(3,2) + v_ave(3,2)) + cv_s_gauss7*(-cv_s_gauss9*vnn_ave(3,2) + vn_ave(3,2))) + cv_s_gauss0*(DN(0,0)*v_ave(0,2) + DN(1,0)*v_ave(1,2) + DN(2,0)*v_ave(2,2) + DN(3,0)*v_ave(3,2)) + cv_s_gauss1*(DN(0,1)*v_ave(0,2) + DN(1,1)*v_ave(1,2) + DN(2,1)*v_ave(2,2) + DN(3,1)*v_ave(3,2)) + cv_s_gauss2*(DN(0,2)*v_ave(0,2) + DN(1,2)*v_ave(1,2) + DN(2,2)*v_ave(2,2) + DN(3,2)*v_ave(3,2))));


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

    const double& dt = data.dt;
    const double& n = data.n;
    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v_ave = data.v_ave;
    const BoundedMatrix<double,nnodes,dim>& vn_ave = data.vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnn_ave = data.vnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnnn_ave = data.vnnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& vconv = v_ave - vmesh;
    const BoundedMatrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& p_ave = data.p_ave;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the error estimator
    array_1d<double,dim> v_s_gauss;
    const array_1d<double,dim> v_gauss = prod(trans(v_ave), N);
    const array_1d<double,dim> f_gauss = prod(trans(f), N);
    const array_1d<double,dim> grad_p = prod(trans(DN), p_ave);

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Gauss point velocity subscale value computation
    const double cv_s_gauss0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cv_s_gauss1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cv_s_gauss2 =             1.0/(rho*stab_c2*sqrt(pow(cv_s_gauss0, 2) + pow(cv_s_gauss1, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cv_s_gauss3 =             bdf0*n;
const double cv_s_gauss4 =             n - 1;
const double cv_s_gauss5 =             cv_s_gauss4/n;
const double cv_s_gauss6 =             bdf1*cv_s_gauss4;
const double cv_s_gauss7 =             n - 2;
const double cv_s_gauss8 =             cv_s_gauss7/cv_s_gauss4;
const double cv_s_gauss9 =             bdf2*cv_s_gauss7;
const double cv_s_gauss10 =             (n - 3)/cv_s_gauss7;
            v_s_gauss[0]=-cv_s_gauss2*(DN(0,0)*p_ave[0] + DN(1,0)*p_ave[1] + DN(2,0)*p_ave[2] + rho*(-N[0]*f(0,0) + N[0]*(cv_s_gauss3*(-cv_s_gauss5*vn_ave(0,0) + v_ave(0,0)) + cv_s_gauss6*(-cv_s_gauss8*vnn_ave(0,0) + vn_ave(0,0)) + cv_s_gauss9*(-cv_s_gauss10*vnnn_ave(0,0) + vnn_ave(0,0))) - N[1]*f(1,0) + N[1]*(cv_s_gauss3*(-cv_s_gauss5*vn_ave(1,0) + v_ave(1,0)) + cv_s_gauss6*(-cv_s_gauss8*vnn_ave(1,0) + vn_ave(1,0)) + cv_s_gauss9*(-cv_s_gauss10*vnnn_ave(1,0) + vnn_ave(1,0))) - N[2]*f(2,0) + N[2]*(cv_s_gauss3*(-cv_s_gauss5*vn_ave(2,0) + v_ave(2,0)) + cv_s_gauss6*(-cv_s_gauss8*vnn_ave(2,0) + vn_ave(2,0)) + cv_s_gauss9*(-cv_s_gauss10*vnnn_ave(2,0) + vnn_ave(2,0))) + cv_s_gauss0*(DN(0,0)*v_ave(0,0) + DN(1,0)*v_ave(1,0) + DN(2,0)*v_ave(2,0)) + cv_s_gauss1*(DN(0,1)*v_ave(0,0) + DN(1,1)*v_ave(1,0) + DN(2,1)*v_ave(2,0))));
            v_s_gauss[1]=-cv_s_gauss2*(DN(0,1)*p_ave[0] + DN(1,1)*p_ave[1] + DN(2,1)*p_ave[2] + rho*(-N[0]*f(0,1) + N[0]*(cv_s_gauss3*(-cv_s_gauss5*vn_ave(0,1) + v_ave(0,1)) + cv_s_gauss6*(-cv_s_gauss8*vnn_ave(0,1) + vn_ave(0,1)) + cv_s_gauss9*(-cv_s_gauss10*vnnn_ave(0,1) + vnn_ave(0,1))) - N[1]*f(1,1) + N[1]*(cv_s_gauss3*(-cv_s_gauss5*vn_ave(1,1) + v_ave(1,1)) + cv_s_gauss6*(-cv_s_gauss8*vnn_ave(1,1) + vn_ave(1,1)) + cv_s_gauss9*(-cv_s_gauss10*vnnn_ave(1,1) + vnn_ave(1,1))) - N[2]*f(2,1) + N[2]*(cv_s_gauss3*(-cv_s_gauss5*vn_ave(2,1) + v_ave(2,1)) + cv_s_gauss6*(-cv_s_gauss8*vnn_ave(2,1) + vn_ave(2,1)) + cv_s_gauss9*(-cv_s_gauss10*vnnn_ave(2,1) + vnn_ave(2,1))) + cv_s_gauss0*(DN(0,0)*v_ave(0,1) + DN(1,0)*v_ave(1,1) + DN(2,0)*v_ave(2,1)) + cv_s_gauss1*(DN(0,1)*v_ave(0,1) + DN(1,1)*v_ave(1,1) + DN(2,1)*v_ave(2,1))));


    const double v_gauss_norm = norm_2(v_gauss);
    const double v_s_gauss_norm = norm_2(v_s_gauss);

    return v_s_gauss_norm/v_gauss_norm;
}

}
