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
    const double& dtn = data.dtn;
    const double& dtnn = data.dtnn;
    const double& dtnnn = data.dtnnn;

    const double& t = data.t;
    const double& tn = data.tn;
    const double& tnn = data.tnn;
    const double& tnnn = data.tnnn;

    const double& bdf0 = data.bdf0;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v_ave = data.v_ave;
    const BoundedMatrix<double,nnodes,dim>& vn_ave = data.vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;

    // Get time accurate expression
    const BoundedMatrix<double,nnodes,dim>& v = (t * v_ave - tn * vn_ave) * ( 1.0 / dt );
    const BoundedMatrix<double, nnodes, dim>& vconv = v - vmesh;

    // Get constitutive matrix
    const Matrix& C = data.C;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    const double clhs0 =             1.0/dt;
const double clhs1 =             clhs0*t;
const double clhs2 =             C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs3 =             C(0,3)*DN(0,0);
const double clhs4 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs3;
const double clhs5 =             C(0,5)*DN(0,0);
const double clhs6 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs5;
const double clhs7 =             pow(DN(0,0), 2);
const double clhs8 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double clhs9 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double clhs10 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double clhs11 =             rho*stab_c2*sqrt(pow(clhs10, 2) + pow(clhs8, 2) + pow(clhs9, 2));
const double clhs12 =             clhs11*h/stab_c1 + mu;
const double clhs13 =             pow(N[0], 2);
const double clhs14 =             bdf0*rho;
const double clhs15 =             N[0]*rho;
const double clhs16 =             DN(0,0)*clhs8 + DN(0,1)*clhs9 + DN(0,2)*clhs10;
const double clhs17 =             N[0]*bdf0 + clhs16;
const double clhs18 =             pow(rho, 2);
const double clhs19 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double clhs20 =             1.0/(clhs11/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs21 =             1.0*N[0]*clhs18*clhs19*clhs20;
const double clhs22 =             1.0*clhs16*clhs18*clhs20;
const double clhs23 =             clhs13*clhs14 + clhs15*clhs16 + clhs17*clhs21 + clhs17*clhs22;
const double clhs24 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs3;
const double clhs25 =             C(1,3)*DN(0,1);
const double clhs26 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs25;
const double clhs27 =             C(3,5)*DN(0,0);
const double clhs28 =             C(4,5)*DN(0,2);
const double clhs29 =             C(1,5)*DN(0,1) + clhs27 + clhs28;
const double clhs30 =             DN(0,0)*clhs12;
const double clhs31 =             DN(0,1)*clhs30;
const double clhs32 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs5;
const double clhs33 =             C(3,4)*DN(0,1);
const double clhs34 =             C(2,3)*DN(0,2) + clhs27 + clhs33;
const double clhs35 =             C(2,5)*DN(0,2);
const double clhs36 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs35;
const double clhs37 =             DN(0,2)*clhs30;
const double clhs38 =             pow(c, -2);
const double clhs39 =             1.0/rho;
const double clhs40 =             N[0]*bdf0*clhs38*clhs39;
const double clhs41 =             1.0*clhs19*clhs20;
const double clhs42 =             1.0*clhs20*rho;
const double clhs43 =             clhs16*clhs42;
const double clhs44 =             clhs0*t*(-N[0] + clhs12*clhs40 + clhs15*clhs41 + clhs43);
const double clhs45 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs46 =             C(0,3)*DN(1,0);
const double clhs47 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs46;
const double clhs48 =             C(0,5)*DN(1,0);
const double clhs49 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs48;
const double clhs50 =             DN(1,0)*clhs30;
const double clhs51 =             N[0]*bdf0*rho;
const double clhs52 =             N[1]*clhs51;
const double clhs53 =             DN(1,0)*clhs8 + DN(1,1)*clhs9 + DN(1,2)*clhs10;
const double clhs54 =             N[1]*bdf0 + clhs53;
const double clhs55 =             clhs15*clhs53 + clhs21*clhs54 + clhs22*clhs54 + clhs52;
const double clhs56 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs46;
const double clhs57 =             C(1,3)*DN(1,1);
const double clhs58 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs57;
const double clhs59 =             C(3,5)*DN(1,0);
const double clhs60 =             C(4,5)*DN(1,2);
const double clhs61 =             C(1,5)*DN(1,1) + clhs59 + clhs60;
const double clhs62 =             DN(1,1)*clhs30;
const double clhs63 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs48;
const double clhs64 =             C(3,4)*DN(1,1);
const double clhs65 =             C(2,3)*DN(1,2) + clhs59 + clhs64;
const double clhs66 =             C(2,5)*DN(1,2);
const double clhs67 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs66;
const double clhs68 =             DN(1,2)*clhs30;
const double clhs69 =             DN(0,0)*N[1];
const double clhs70 =             bdf0*clhs12*clhs38*clhs39;
const double clhs71 =             DN(1,0)*N[0];
const double clhs72 =             1.0*clhs19*clhs20*rho;
const double clhs73 =             1.0*DN(1,0)*clhs20*rho;
const double clhs74 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs75 =             C(0,3)*DN(2,0);
const double clhs76 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs75;
const double clhs77 =             C(0,5)*DN(2,0);
const double clhs78 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs77;
const double clhs79 =             DN(2,0)*clhs30;
const double clhs80 =             N[2]*clhs51;
const double clhs81 =             DN(2,0)*clhs8 + DN(2,1)*clhs9 + DN(2,2)*clhs10;
const double clhs82 =             N[2]*bdf0;
const double clhs83 =             clhs81 + clhs82;
const double clhs84 =             clhs15*clhs81 + clhs21*clhs83 + clhs22*clhs83 + clhs80;
const double clhs85 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs75;
const double clhs86 =             C(1,3)*DN(2,1);
const double clhs87 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs86;
const double clhs88 =             C(3,5)*DN(2,0);
const double clhs89 =             C(4,5)*DN(2,2);
const double clhs90 =             C(1,5)*DN(2,1) + clhs88 + clhs89;
const double clhs91 =             DN(2,1)*clhs30;
const double clhs92 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs77;
const double clhs93 =             C(3,4)*DN(2,1);
const double clhs94 =             C(2,3)*DN(2,2) + clhs88 + clhs93;
const double clhs95 =             C(2,5)*DN(2,2);
const double clhs96 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs95;
const double clhs97 =             DN(2,2)*clhs30;
const double clhs98 =             DN(0,0)*N[2];
const double clhs99 =             DN(2,0)*N[0];
const double clhs100 =             1.0*DN(2,0)*clhs20*rho;
const double clhs101 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs102 =             C(0,3)*DN(3,0);
const double clhs103 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs102;
const double clhs104 =             C(0,5)*DN(3,0);
const double clhs105 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs104;
const double clhs106 =             DN(3,0)*clhs30;
const double clhs107 =             N[3]*clhs51;
const double clhs108 =             DN(3,0)*clhs8 + DN(3,1)*clhs9 + DN(3,2)*clhs10;
const double clhs109 =             N[3]*bdf0;
const double clhs110 =             clhs108 + clhs109;
const double clhs111 =             clhs107 + clhs108*clhs15 + clhs110*clhs21 + clhs110*clhs22;
const double clhs112 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs102;
const double clhs113 =             C(1,3)*DN(3,1);
const double clhs114 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs113;
const double clhs115 =             C(3,5)*DN(3,0);
const double clhs116 =             C(4,5)*DN(3,2);
const double clhs117 =             C(1,5)*DN(3,1) + clhs115 + clhs116;
const double clhs118 =             DN(3,1)*clhs30;
const double clhs119 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs104;
const double clhs120 =             C(3,4)*DN(3,1);
const double clhs121 =             C(2,3)*DN(3,2) + clhs115 + clhs120;
const double clhs122 =             C(2,5)*DN(3,2);
const double clhs123 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs122;
const double clhs124 =             DN(3,2)*clhs30;
const double clhs125 =             DN(0,0)*N[3];
const double clhs126 =             DN(3,0)*N[0];
const double clhs127 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs25;
const double clhs128 =             C(0,4)*DN(0,0) + clhs28 + clhs33;
const double clhs129 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs130 =             C(1,4)*DN(0,1);
const double clhs131 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs130;
const double clhs132 =             pow(DN(0,1), 2);
const double clhs133 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs130;
const double clhs134 =             C(2,4)*DN(0,2);
const double clhs135 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs134;
const double clhs136 =             DN(0,1)*clhs12;
const double clhs137 =             DN(0,2)*clhs136;
const double clhs138 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs57;
const double clhs139 =             C(0,4)*DN(1,0) + clhs60 + clhs64;
const double clhs140 =             DN(1,0)*clhs136;
const double clhs141 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs142 =             C(1,4)*DN(1,1);
const double clhs143 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs142;
const double clhs144 =             DN(1,1)*clhs136;
const double clhs145 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs142;
const double clhs146 =             C(2,4)*DN(1,2);
const double clhs147 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs146;
const double clhs148 =             DN(1,2)*clhs136;
const double clhs149 =             DN(0,1)*N[1];
const double clhs150 =             DN(1,1)*N[0];
const double clhs151 =             1.0*DN(1,1)*clhs20*rho;
const double clhs152 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs86;
const double clhs153 =             C(0,4)*DN(2,0) + clhs89 + clhs93;
const double clhs154 =             DN(2,0)*clhs136;
const double clhs155 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs156 =             C(1,4)*DN(2,1);
const double clhs157 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs156;
const double clhs158 =             DN(2,1)*clhs136;
const double clhs159 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs156;
const double clhs160 =             C(2,4)*DN(2,2);
const double clhs161 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs160;
const double clhs162 =             DN(2,2)*clhs136;
const double clhs163 =             DN(0,1)*N[2];
const double clhs164 =             DN(2,1)*N[0];
const double clhs165 =             1.0*DN(2,1)*clhs20*rho;
const double clhs166 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs113;
const double clhs167 =             C(0,4)*DN(3,0) + clhs116 + clhs120;
const double clhs168 =             DN(3,0)*clhs136;
const double clhs169 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs170 =             C(1,4)*DN(3,1);
const double clhs171 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs170;
const double clhs172 =             DN(3,1)*clhs136;
const double clhs173 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs170;
const double clhs174 =             C(2,4)*DN(3,2);
const double clhs175 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs174;
const double clhs176 =             DN(3,2)*clhs136;
const double clhs177 =             DN(0,1)*N[3];
const double clhs178 =             DN(3,1)*N[0];
const double clhs179 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs35;
const double clhs180 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs134;
const double clhs181 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs182 =             pow(DN(0,2), 2);
const double clhs183 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs66;
const double clhs184 =             DN(0,2)*clhs12;
const double clhs185 =             DN(1,0)*clhs184;
const double clhs186 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs146;
const double clhs187 =             DN(1,1)*clhs184;
const double clhs188 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs189 =             DN(1,2)*clhs184;
const double clhs190 =             DN(0,2)*N[1];
const double clhs191 =             DN(1,2)*N[0];
const double clhs192 =             1.0*DN(1,2)*clhs20*rho;
const double clhs193 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs95;
const double clhs194 =             DN(2,0)*clhs184;
const double clhs195 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs160;
const double clhs196 =             DN(2,1)*clhs184;
const double clhs197 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs198 =             DN(2,2)*clhs184;
const double clhs199 =             DN(0,2)*N[2];
const double clhs200 =             DN(2,2)*N[0];
const double clhs201 =             1.0*DN(2,2)*clhs20*rho;
const double clhs202 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs122;
const double clhs203 =             DN(3,0)*clhs184;
const double clhs204 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs174;
const double clhs205 =             DN(3,1)*clhs184;
const double clhs206 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs207 =             DN(3,2)*clhs184;
const double clhs208 =             DN(0,2)*N[3];
const double clhs209 =             DN(3,2)*N[0];
const double clhs210 =             clhs17*clhs42;
const double clhs211 =             clhs0*t*(N[0] + clhs210);
const double clhs212 =             bdf0*clhs38*clhs39;
const double clhs213 =             1.0*clhs20;
const double clhs214 =             1.0*DN(0,0)*clhs20*rho;
const double clhs215 =             1.0*DN(0,1)*clhs20*rho;
const double clhs216 =             1.0*DN(0,2)*clhs20*rho;
const double clhs217 =             1.0*DN(0,0)*clhs20;
const double clhs218 =             1.0*DN(0,1)*clhs20;
const double clhs219 =             1.0*DN(0,2)*clhs20;
const double clhs220 =             clhs1*(DN(1,0)*clhs217 + DN(1,1)*clhs218 + DN(1,2)*clhs219 + N[1]*clhs40);
const double clhs221 =             clhs1*(DN(2,0)*clhs217 + DN(2,1)*clhs218 + DN(2,2)*clhs219 + N[2]*clhs40);
const double clhs222 =             clhs1*(DN(3,0)*clhs217 + DN(3,1)*clhs218 + DN(3,2)*clhs219 + N[3]*clhs40);
const double clhs223 =             N[1]*rho;
const double clhs224 =             1.0*N[1]*clhs18*clhs19*clhs20;
const double clhs225 =             1.0*clhs18*clhs20*clhs53;
const double clhs226 =             clhs16*clhs223 + clhs17*clhs224 + clhs17*clhs225 + clhs52;
const double clhs227 =             pow(DN(1,0), 2);
const double clhs228 =             pow(N[1], 2);
const double clhs229 =             clhs14*clhs228 + clhs223*clhs53 + clhs224*clhs54 + clhs225*clhs54;
const double clhs230 =             DN(1,0)*clhs12;
const double clhs231 =             DN(1,1)*clhs230;
const double clhs232 =             DN(1,2)*clhs230;
const double clhs233 =             N[1]*bdf0*clhs38*clhs39;
const double clhs234 =             clhs42*clhs53;
const double clhs235 =             clhs0*t*(-N[1] + clhs12*clhs233 + clhs223*clhs41 + clhs234);
const double clhs236 =             DN(2,0)*clhs230;
const double clhs237 =             N[1]*bdf0*rho;
const double clhs238 =             N[2]*clhs237;
const double clhs239 =             clhs223*clhs81 + clhs224*clhs83 + clhs225*clhs83 + clhs238;
const double clhs240 =             DN(2,1)*clhs230;
const double clhs241 =             DN(2,2)*clhs230;
const double clhs242 =             DN(1,0)*N[2];
const double clhs243 =             DN(2,0)*N[1];
const double clhs244 =             DN(3,0)*clhs230;
const double clhs245 =             N[3]*clhs237;
const double clhs246 =             clhs108*clhs223 + clhs110*clhs224 + clhs110*clhs225 + clhs245;
const double clhs247 =             DN(3,1)*clhs230;
const double clhs248 =             DN(3,2)*clhs230;
const double clhs249 =             DN(1,0)*N[3];
const double clhs250 =             DN(3,0)*N[1];
const double clhs251 =             pow(DN(1,1), 2);
const double clhs252 =             DN(1,1)*clhs12;
const double clhs253 =             DN(1,2)*clhs252;
const double clhs254 =             DN(2,0)*clhs252;
const double clhs255 =             DN(2,1)*clhs252;
const double clhs256 =             DN(2,2)*clhs252;
const double clhs257 =             DN(1,1)*N[2];
const double clhs258 =             DN(2,1)*N[1];
const double clhs259 =             DN(3,0)*clhs252;
const double clhs260 =             DN(3,1)*clhs252;
const double clhs261 =             DN(3,2)*clhs252;
const double clhs262 =             DN(1,1)*N[3];
const double clhs263 =             DN(3,1)*N[1];
const double clhs264 =             pow(DN(1,2), 2);
const double clhs265 =             DN(1,2)*clhs12;
const double clhs266 =             DN(2,0)*clhs265;
const double clhs267 =             DN(2,1)*clhs265;
const double clhs268 =             DN(2,2)*clhs265;
const double clhs269 =             DN(1,2)*N[2];
const double clhs270 =             DN(2,2)*N[1];
const double clhs271 =             DN(3,0)*clhs265;
const double clhs272 =             DN(3,1)*clhs265;
const double clhs273 =             DN(3,2)*clhs265;
const double clhs274 =             DN(1,2)*N[3];
const double clhs275 =             DN(3,2)*N[1];
const double clhs276 =             clhs42*clhs54;
const double clhs277 =             clhs0*t*(N[1] + clhs276);
const double clhs278 =             1.0*DN(1,0)*clhs20;
const double clhs279 =             1.0*DN(1,1)*clhs20;
const double clhs280 =             1.0*DN(1,2)*clhs20;
const double clhs281 =             clhs1*(DN(2,0)*clhs278 + DN(2,1)*clhs279 + DN(2,2)*clhs280 + N[2]*clhs233);
const double clhs282 =             clhs1*(DN(3,0)*clhs278 + DN(3,1)*clhs279 + DN(3,2)*clhs280 + N[3]*clhs233);
const double clhs283 =             N[2]*rho;
const double clhs284 =             1.0*N[2]*clhs18*clhs19*clhs20;
const double clhs285 =             1.0*clhs18*clhs20*clhs81;
const double clhs286 =             clhs16*clhs283 + clhs17*clhs284 + clhs17*clhs285 + clhs80;
const double clhs287 =             clhs238 + clhs283*clhs53 + clhs284*clhs54 + clhs285*clhs54;
const double clhs288 =             pow(DN(2,0), 2);
const double clhs289 =             pow(N[2], 2);
const double clhs290 =             clhs14*clhs289 + clhs283*clhs81 + clhs284*clhs83 + clhs285*clhs83;
const double clhs291 =             DN(2,0)*clhs12;
const double clhs292 =             DN(2,1)*clhs291;
const double clhs293 =             DN(2,2)*clhs291;
const double clhs294 =             clhs12*clhs38*clhs39;
const double clhs295 =             clhs42*clhs81;
const double clhs296 =             clhs0*t*(-N[2] + clhs283*clhs41 + clhs294*clhs82 + clhs295);
const double clhs297 =             DN(3,0)*clhs291;
const double clhs298 =             N[2]*N[3]*bdf0;
const double clhs299 =             clhs298*rho;
const double clhs300 =             clhs108*clhs283 + clhs110*clhs284 + clhs110*clhs285 + clhs299;
const double clhs301 =             DN(3,1)*clhs291;
const double clhs302 =             DN(3,2)*clhs291;
const double clhs303 =             DN(2,0)*N[3];
const double clhs304 =             DN(3,0)*N[2];
const double clhs305 =             pow(DN(2,1), 2);
const double clhs306 =             DN(2,1)*clhs12;
const double clhs307 =             DN(2,2)*clhs306;
const double clhs308 =             DN(3,0)*clhs306;
const double clhs309 =             DN(3,1)*clhs306;
const double clhs310 =             DN(3,2)*clhs306;
const double clhs311 =             DN(2,1)*N[3];
const double clhs312 =             DN(3,1)*N[2];
const double clhs313 =             pow(DN(2,2), 2);
const double clhs314 =             DN(2,2)*clhs12;
const double clhs315 =             DN(3,0)*clhs314;
const double clhs316 =             DN(3,1)*clhs314;
const double clhs317 =             DN(3,2)*clhs314;
const double clhs318 =             DN(2,2)*N[3];
const double clhs319 =             DN(3,2)*N[2];
const double clhs320 =             clhs42*clhs83;
const double clhs321 =             clhs0*t*(N[2] + clhs320);
const double clhs322 =             clhs1*(1.0*DN(2,0)*DN(3,0)*clhs20 + 1.0*DN(2,1)*DN(3,1)*clhs20 + 1.0*DN(2,2)*DN(3,2)*clhs20 + clhs298*clhs38*clhs39);
const double clhs323 =             N[3]*rho;
const double clhs324 =             1.0*N[3]*clhs18*clhs19*clhs20;
const double clhs325 =             1.0*clhs108*clhs18*clhs20;
const double clhs326 =             clhs107 + clhs16*clhs323 + clhs17*clhs324 + clhs17*clhs325;
const double clhs327 =             clhs245 + clhs323*clhs53 + clhs324*clhs54 + clhs325*clhs54;
const double clhs328 =             clhs299 + clhs323*clhs81 + clhs324*clhs83 + clhs325*clhs83;
const double clhs329 =             pow(DN(3,0), 2);
const double clhs330 =             pow(N[3], 2);
const double clhs331 =             clhs108*clhs323 + clhs110*clhs324 + clhs110*clhs325 + clhs14*clhs330;
const double clhs332 =             DN(3,0)*clhs12;
const double clhs333 =             DN(3,1)*clhs332;
const double clhs334 =             DN(3,2)*clhs332;
const double clhs335 =             clhs0*t*(-N[3] + clhs108*clhs42 + clhs109*clhs294 + clhs323*clhs41);
const double clhs336 =             pow(DN(3,1), 2);
const double clhs337 =             DN(3,1)*DN(3,2)*clhs12;
const double clhs338 =             pow(DN(3,2), 2);
const double clhs339 =             clhs0*t*(N[3] + clhs110*clhs42);
            lhs(0,0)=clhs1*(DN(0,0)*clhs2 + DN(0,1)*clhs4 + DN(0,2)*clhs6 + clhs12*clhs7 + clhs23);
            lhs(0,1)=clhs1*(DN(0,0)*clhs24 + DN(0,1)*clhs26 + DN(0,2)*clhs29 + clhs31);
            lhs(0,2)=clhs1*(DN(0,0)*clhs32 + DN(0,1)*clhs34 + DN(0,2)*clhs36 + clhs37);
            lhs(0,3)=DN(0,0)*clhs44;
            lhs(0,4)=clhs1*(DN(0,0)*clhs45 + DN(0,1)*clhs47 + DN(0,2)*clhs49 + clhs50 + clhs55);
            lhs(0,5)=clhs1*(DN(0,0)*clhs56 + DN(0,1)*clhs58 + DN(0,2)*clhs61 + clhs62);
            lhs(0,6)=clhs1*(DN(0,0)*clhs63 + DN(0,1)*clhs65 + DN(0,2)*clhs67 + clhs68);
            lhs(0,7)=clhs1*(clhs16*clhs73 + clhs69*clhs70 - clhs69 + clhs71*clhs72);
            lhs(0,8)=clhs1*(DN(0,0)*clhs74 + DN(0,1)*clhs76 + DN(0,2)*clhs78 + clhs79 + clhs84);
            lhs(0,9)=clhs1*(DN(0,0)*clhs85 + DN(0,1)*clhs87 + DN(0,2)*clhs90 + clhs91);
            lhs(0,10)=clhs1*(DN(0,0)*clhs92 + DN(0,1)*clhs94 + DN(0,2)*clhs96 + clhs97);
            lhs(0,11)=clhs1*(clhs100*clhs16 + clhs70*clhs98 + clhs72*clhs99 - clhs98);
            lhs(0,12)=clhs1*(DN(0,0)*clhs101 + DN(0,1)*clhs103 + DN(0,2)*clhs105 + clhs106 + clhs111);
            lhs(0,13)=clhs1*(DN(0,0)*clhs112 + DN(0,1)*clhs114 + DN(0,2)*clhs117 + clhs118);
            lhs(0,14)=clhs1*(DN(0,0)*clhs119 + DN(0,1)*clhs121 + DN(0,2)*clhs123 + clhs124);
            lhs(0,15)=clhs1*(DN(3,0)*clhs43 + clhs125*clhs70 - clhs125 + clhs126*clhs72);
            lhs(1,0)=clhs1*(DN(0,0)*clhs4 + DN(0,1)*clhs127 + DN(0,2)*clhs128 + clhs31);
            lhs(1,1)=clhs1*(DN(0,0)*clhs26 + DN(0,1)*clhs129 + DN(0,2)*clhs131 + clhs12*clhs132 + clhs23);
            lhs(1,2)=clhs1*(DN(0,0)*clhs34 + DN(0,1)*clhs133 + DN(0,2)*clhs135 + clhs137);
            lhs(1,3)=DN(0,1)*clhs44;
            lhs(1,4)=clhs1*(DN(0,0)*clhs47 + DN(0,1)*clhs138 + DN(0,2)*clhs139 + clhs140);
            lhs(1,5)=clhs1*(DN(0,0)*clhs58 + DN(0,1)*clhs141 + DN(0,2)*clhs143 + clhs144 + clhs55);
            lhs(1,6)=clhs1*(DN(0,0)*clhs65 + DN(0,1)*clhs145 + DN(0,2)*clhs147 + clhs148);
            lhs(1,7)=clhs1*(clhs149*clhs70 - clhs149 + clhs150*clhs72 + clhs151*clhs16);
            lhs(1,8)=clhs1*(DN(0,0)*clhs76 + DN(0,1)*clhs152 + DN(0,2)*clhs153 + clhs154);
            lhs(1,9)=clhs1*(DN(0,0)*clhs87 + DN(0,1)*clhs155 + DN(0,2)*clhs157 + clhs158 + clhs84);
            lhs(1,10)=clhs1*(DN(0,0)*clhs94 + DN(0,1)*clhs159 + DN(0,2)*clhs161 + clhs162);
            lhs(1,11)=clhs1*(clhs16*clhs165 + clhs163*clhs70 - clhs163 + clhs164*clhs72);
            lhs(1,12)=clhs1*(DN(0,0)*clhs103 + DN(0,1)*clhs166 + DN(0,2)*clhs167 + clhs168);
            lhs(1,13)=clhs1*(DN(0,0)*clhs114 + DN(0,1)*clhs169 + DN(0,2)*clhs171 + clhs111 + clhs172);
            lhs(1,14)=clhs1*(DN(0,0)*clhs121 + DN(0,1)*clhs173 + DN(0,2)*clhs175 + clhs176);
            lhs(1,15)=clhs1*(DN(3,1)*clhs43 + clhs177*clhs70 - clhs177 + clhs178*clhs72);
            lhs(2,0)=clhs1*(DN(0,0)*clhs6 + DN(0,1)*clhs128 + DN(0,2)*clhs179 + clhs37);
            lhs(2,1)=clhs1*(DN(0,0)*clhs29 + DN(0,1)*clhs131 + DN(0,2)*clhs180 + clhs137);
            lhs(2,2)=clhs1*(DN(0,0)*clhs36 + DN(0,1)*clhs135 + DN(0,2)*clhs181 + clhs12*clhs182 + clhs23);
            lhs(2,3)=DN(0,2)*clhs44;
            lhs(2,4)=clhs1*(DN(0,0)*clhs49 + DN(0,1)*clhs139 + DN(0,2)*clhs183 + clhs185);
            lhs(2,5)=clhs1*(DN(0,0)*clhs61 + DN(0,1)*clhs143 + DN(0,2)*clhs186 + clhs187);
            lhs(2,6)=clhs1*(DN(0,0)*clhs67 + DN(0,1)*clhs147 + DN(0,2)*clhs188 + clhs189 + clhs55);
            lhs(2,7)=clhs1*(clhs16*clhs192 + clhs190*clhs70 - clhs190 + clhs191*clhs72);
            lhs(2,8)=clhs1*(DN(0,0)*clhs78 + DN(0,1)*clhs153 + DN(0,2)*clhs193 + clhs194);
            lhs(2,9)=clhs1*(DN(0,0)*clhs90 + DN(0,1)*clhs157 + DN(0,2)*clhs195 + clhs196);
            lhs(2,10)=clhs1*(DN(0,0)*clhs96 + DN(0,1)*clhs161 + DN(0,2)*clhs197 + clhs198 + clhs84);
            lhs(2,11)=clhs1*(clhs16*clhs201 + clhs199*clhs70 - clhs199 + clhs200*clhs72);
            lhs(2,12)=clhs1*(DN(0,0)*clhs105 + DN(0,1)*clhs167 + DN(0,2)*clhs202 + clhs203);
            lhs(2,13)=clhs1*(DN(0,0)*clhs117 + DN(0,1)*clhs171 + DN(0,2)*clhs204 + clhs205);
            lhs(2,14)=clhs1*(DN(0,0)*clhs123 + DN(0,1)*clhs175 + DN(0,2)*clhs206 + clhs111 + clhs207);
            lhs(2,15)=clhs1*(DN(3,2)*clhs43 + clhs208*clhs70 - clhs208 + clhs209*clhs72);
            lhs(3,0)=DN(0,0)*clhs211;
            lhs(3,1)=DN(0,1)*clhs211;
            lhs(3,2)=DN(0,2)*clhs211;
            lhs(3,3)=clhs1*(clhs13*clhs212 + clhs132*clhs213 + clhs182*clhs213 + clhs213*clhs7);
            lhs(3,4)=clhs1*(clhs214*clhs54 + clhs71);
            lhs(3,5)=clhs1*(clhs150 + clhs215*clhs54);
            lhs(3,6)=clhs1*(clhs191 + clhs216*clhs54);
            lhs(3,7)=clhs220;
            lhs(3,8)=clhs1*(clhs214*clhs83 + clhs99);
            lhs(3,9)=clhs1*(clhs164 + clhs215*clhs83);
            lhs(3,10)=clhs1*(clhs200 + clhs216*clhs83);
            lhs(3,11)=clhs221;
            lhs(3,12)=clhs1*(clhs110*clhs214 + clhs126);
            lhs(3,13)=clhs1*(clhs110*clhs215 + clhs178);
            lhs(3,14)=clhs1*(clhs110*clhs216 + clhs209);
            lhs(3,15)=clhs222;
            lhs(4,0)=clhs1*(DN(1,0)*clhs2 + DN(1,1)*clhs4 + DN(1,2)*clhs6 + clhs226 + clhs50);
            lhs(4,1)=clhs1*(DN(1,0)*clhs24 + DN(1,1)*clhs26 + DN(1,2)*clhs29 + clhs140);
            lhs(4,2)=clhs1*(DN(1,0)*clhs32 + DN(1,1)*clhs34 + DN(1,2)*clhs36 + clhs185);
            lhs(4,3)=clhs1*(clhs214*clhs53 + clhs69*clhs72 + clhs70*clhs71 - clhs71);
            lhs(4,4)=clhs1*(DN(1,0)*clhs45 + DN(1,1)*clhs47 + DN(1,2)*clhs49 + clhs12*clhs227 + clhs229);
            lhs(4,5)=clhs1*(DN(1,0)*clhs56 + DN(1,1)*clhs58 + DN(1,2)*clhs61 + clhs231);
            lhs(4,6)=clhs1*(DN(1,0)*clhs63 + DN(1,1)*clhs65 + DN(1,2)*clhs67 + clhs232);
            lhs(4,7)=DN(1,0)*clhs235;
            lhs(4,8)=clhs1*(DN(1,0)*clhs74 + DN(1,1)*clhs76 + DN(1,2)*clhs78 + clhs236 + clhs239);
            lhs(4,9)=clhs1*(DN(1,0)*clhs85 + DN(1,1)*clhs87 + DN(1,2)*clhs90 + clhs240);
            lhs(4,10)=clhs1*(DN(1,0)*clhs92 + DN(1,1)*clhs94 + DN(1,2)*clhs96 + clhs241);
            lhs(4,11)=clhs1*(clhs100*clhs53 + clhs242*clhs70 - clhs242 + clhs243*clhs72);
            lhs(4,12)=clhs1*(DN(1,0)*clhs101 + DN(1,1)*clhs103 + DN(1,2)*clhs105 + clhs244 + clhs246);
            lhs(4,13)=clhs1*(DN(1,0)*clhs112 + DN(1,1)*clhs114 + DN(1,2)*clhs117 + clhs247);
            lhs(4,14)=clhs1*(DN(1,0)*clhs119 + DN(1,1)*clhs121 + DN(1,2)*clhs123 + clhs248);
            lhs(4,15)=clhs1*(DN(3,0)*clhs234 + clhs249*clhs70 - clhs249 + clhs250*clhs72);
            lhs(5,0)=clhs1*(DN(1,0)*clhs4 + DN(1,1)*clhs127 + DN(1,2)*clhs128 + clhs62);
            lhs(5,1)=clhs1*(DN(1,0)*clhs26 + DN(1,1)*clhs129 + DN(1,2)*clhs131 + clhs144 + clhs226);
            lhs(5,2)=clhs1*(DN(1,0)*clhs34 + DN(1,1)*clhs133 + DN(1,2)*clhs135 + clhs187);
            lhs(5,3)=clhs1*(clhs149*clhs72 + clhs150*clhs70 - clhs150 + clhs215*clhs53);
            lhs(5,4)=clhs1*(DN(1,0)*clhs47 + DN(1,1)*clhs138 + DN(1,2)*clhs139 + clhs231);
            lhs(5,5)=clhs1*(DN(1,0)*clhs58 + DN(1,1)*clhs141 + DN(1,2)*clhs143 + clhs12*clhs251 + clhs229);
            lhs(5,6)=clhs1*(DN(1,0)*clhs65 + DN(1,1)*clhs145 + DN(1,2)*clhs147 + clhs253);
            lhs(5,7)=DN(1,1)*clhs235;
            lhs(5,8)=clhs1*(DN(1,0)*clhs76 + DN(1,1)*clhs152 + DN(1,2)*clhs153 + clhs254);
            lhs(5,9)=clhs1*(DN(1,0)*clhs87 + DN(1,1)*clhs155 + DN(1,2)*clhs157 + clhs239 + clhs255);
            lhs(5,10)=clhs1*(DN(1,0)*clhs94 + DN(1,1)*clhs159 + DN(1,2)*clhs161 + clhs256);
            lhs(5,11)=clhs1*(clhs165*clhs53 + clhs257*clhs70 - clhs257 + clhs258*clhs72);
            lhs(5,12)=clhs1*(DN(1,0)*clhs103 + DN(1,1)*clhs166 + DN(1,2)*clhs167 + clhs259);
            lhs(5,13)=clhs1*(DN(1,0)*clhs114 + DN(1,1)*clhs169 + DN(1,2)*clhs171 + clhs246 + clhs260);
            lhs(5,14)=clhs1*(DN(1,0)*clhs121 + DN(1,1)*clhs173 + DN(1,2)*clhs175 + clhs261);
            lhs(5,15)=clhs1*(DN(3,1)*clhs234 + clhs262*clhs70 - clhs262 + clhs263*clhs72);
            lhs(6,0)=clhs1*(DN(1,0)*clhs6 + DN(1,1)*clhs128 + DN(1,2)*clhs179 + clhs68);
            lhs(6,1)=clhs1*(DN(1,0)*clhs29 + DN(1,1)*clhs131 + DN(1,2)*clhs180 + clhs148);
            lhs(6,2)=clhs1*(DN(1,0)*clhs36 + DN(1,1)*clhs135 + DN(1,2)*clhs181 + clhs189 + clhs226);
            lhs(6,3)=clhs1*(clhs190*clhs72 + clhs191*clhs70 - clhs191 + clhs216*clhs53);
            lhs(6,4)=clhs1*(DN(1,0)*clhs49 + DN(1,1)*clhs139 + DN(1,2)*clhs183 + clhs232);
            lhs(6,5)=clhs1*(DN(1,0)*clhs61 + DN(1,1)*clhs143 + DN(1,2)*clhs186 + clhs253);
            lhs(6,6)=clhs1*(DN(1,0)*clhs67 + DN(1,1)*clhs147 + DN(1,2)*clhs188 + clhs12*clhs264 + clhs229);
            lhs(6,7)=DN(1,2)*clhs235;
            lhs(6,8)=clhs1*(DN(1,0)*clhs78 + DN(1,1)*clhs153 + DN(1,2)*clhs193 + clhs266);
            lhs(6,9)=clhs1*(DN(1,0)*clhs90 + DN(1,1)*clhs157 + DN(1,2)*clhs195 + clhs267);
            lhs(6,10)=clhs1*(DN(1,0)*clhs96 + DN(1,1)*clhs161 + DN(1,2)*clhs197 + clhs239 + clhs268);
            lhs(6,11)=clhs1*(clhs201*clhs53 + clhs269*clhs70 - clhs269 + clhs270*clhs72);
            lhs(6,12)=clhs1*(DN(1,0)*clhs105 + DN(1,1)*clhs167 + DN(1,2)*clhs202 + clhs271);
            lhs(6,13)=clhs1*(DN(1,0)*clhs117 + DN(1,1)*clhs171 + DN(1,2)*clhs204 + clhs272);
            lhs(6,14)=clhs1*(DN(1,0)*clhs123 + DN(1,1)*clhs175 + DN(1,2)*clhs206 + clhs246 + clhs273);
            lhs(6,15)=clhs1*(DN(3,2)*clhs234 + clhs274*clhs70 - clhs274 + clhs275*clhs72);
            lhs(7,0)=clhs1*(clhs17*clhs73 + clhs69);
            lhs(7,1)=clhs1*(clhs149 + clhs151*clhs17);
            lhs(7,2)=clhs1*(clhs17*clhs192 + clhs190);
            lhs(7,3)=clhs220;
            lhs(7,4)=DN(1,0)*clhs277;
            lhs(7,5)=DN(1,1)*clhs277;
            lhs(7,6)=DN(1,2)*clhs277;
            lhs(7,7)=clhs1*(clhs212*clhs228 + clhs213*clhs227 + clhs213*clhs251 + clhs213*clhs264);
            lhs(7,8)=clhs1*(clhs243 + clhs73*clhs83);
            lhs(7,9)=clhs1*(clhs151*clhs83 + clhs258);
            lhs(7,10)=clhs1*(clhs192*clhs83 + clhs270);
            lhs(7,11)=clhs281;
            lhs(7,12)=clhs1*(clhs110*clhs73 + clhs250);
            lhs(7,13)=clhs1*(clhs110*clhs151 + clhs263);
            lhs(7,14)=clhs1*(clhs110*clhs192 + clhs275);
            lhs(7,15)=clhs282;
            lhs(8,0)=clhs1*(DN(2,0)*clhs2 + DN(2,1)*clhs4 + DN(2,2)*clhs6 + clhs286 + clhs79);
            lhs(8,1)=clhs1*(DN(2,0)*clhs24 + DN(2,1)*clhs26 + DN(2,2)*clhs29 + clhs154);
            lhs(8,2)=clhs1*(DN(2,0)*clhs32 + DN(2,1)*clhs34 + DN(2,2)*clhs36 + clhs194);
            lhs(8,3)=clhs1*(clhs214*clhs81 + clhs70*clhs99 + clhs72*clhs98 - clhs99);
            lhs(8,4)=clhs1*(DN(2,0)*clhs45 + DN(2,1)*clhs47 + DN(2,2)*clhs49 + clhs236 + clhs287);
            lhs(8,5)=clhs1*(DN(2,0)*clhs56 + DN(2,1)*clhs58 + DN(2,2)*clhs61 + clhs254);
            lhs(8,6)=clhs1*(DN(2,0)*clhs63 + DN(2,1)*clhs65 + DN(2,2)*clhs67 + clhs266);
            lhs(8,7)=clhs1*(clhs242*clhs72 + clhs243*clhs70 - clhs243 + clhs73*clhs81);
            lhs(8,8)=clhs1*(DN(2,0)*clhs74 + DN(2,1)*clhs76 + DN(2,2)*clhs78 + clhs12*clhs288 + clhs290);
            lhs(8,9)=clhs1*(DN(2,0)*clhs85 + DN(2,1)*clhs87 + DN(2,2)*clhs90 + clhs292);
            lhs(8,10)=clhs1*(DN(2,0)*clhs92 + DN(2,1)*clhs94 + DN(2,2)*clhs96 + clhs293);
            lhs(8,11)=DN(2,0)*clhs296;
            lhs(8,12)=clhs1*(DN(2,0)*clhs101 + DN(2,1)*clhs103 + DN(2,2)*clhs105 + clhs297 + clhs300);
            lhs(8,13)=clhs1*(DN(2,0)*clhs112 + DN(2,1)*clhs114 + DN(2,2)*clhs117 + clhs301);
            lhs(8,14)=clhs1*(DN(2,0)*clhs119 + DN(2,1)*clhs121 + DN(2,2)*clhs123 + clhs302);
            lhs(8,15)=clhs1*(DN(3,0)*clhs295 + clhs303*clhs70 - clhs303 + clhs304*clhs72);
            lhs(9,0)=clhs1*(DN(2,0)*clhs4 + DN(2,1)*clhs127 + DN(2,2)*clhs128 + clhs91);
            lhs(9,1)=clhs1*(DN(2,0)*clhs26 + DN(2,1)*clhs129 + DN(2,2)*clhs131 + clhs158 + clhs286);
            lhs(9,2)=clhs1*(DN(2,0)*clhs34 + DN(2,1)*clhs133 + DN(2,2)*clhs135 + clhs196);
            lhs(9,3)=clhs1*(clhs163*clhs72 + clhs164*clhs70 - clhs164 + clhs215*clhs81);
            lhs(9,4)=clhs1*(DN(2,0)*clhs47 + DN(2,1)*clhs138 + DN(2,2)*clhs139 + clhs240);
            lhs(9,5)=clhs1*(DN(2,0)*clhs58 + DN(2,1)*clhs141 + DN(2,2)*clhs143 + clhs255 + clhs287);
            lhs(9,6)=clhs1*(DN(2,0)*clhs65 + DN(2,1)*clhs145 + DN(2,2)*clhs147 + clhs267);
            lhs(9,7)=clhs1*(clhs151*clhs81 + clhs257*clhs72 + clhs258*clhs70 - clhs258);
            lhs(9,8)=clhs1*(DN(2,0)*clhs76 + DN(2,1)*clhs152 + DN(2,2)*clhs153 + clhs292);
            lhs(9,9)=clhs1*(DN(2,0)*clhs87 + DN(2,1)*clhs155 + DN(2,2)*clhs157 + clhs12*clhs305 + clhs290);
            lhs(9,10)=clhs1*(DN(2,0)*clhs94 + DN(2,1)*clhs159 + DN(2,2)*clhs161 + clhs307);
            lhs(9,11)=DN(2,1)*clhs296;
            lhs(9,12)=clhs1*(DN(2,0)*clhs103 + DN(2,1)*clhs166 + DN(2,2)*clhs167 + clhs308);
            lhs(9,13)=clhs1*(DN(2,0)*clhs114 + DN(2,1)*clhs169 + DN(2,2)*clhs171 + clhs300 + clhs309);
            lhs(9,14)=clhs1*(DN(2,0)*clhs121 + DN(2,1)*clhs173 + DN(2,2)*clhs175 + clhs310);
            lhs(9,15)=clhs1*(DN(3,1)*clhs295 + clhs311*clhs70 - clhs311 + clhs312*clhs72);
            lhs(10,0)=clhs1*(DN(2,0)*clhs6 + DN(2,1)*clhs128 + DN(2,2)*clhs179 + clhs97);
            lhs(10,1)=clhs1*(DN(2,0)*clhs29 + DN(2,1)*clhs131 + DN(2,2)*clhs180 + clhs162);
            lhs(10,2)=clhs1*(DN(2,0)*clhs36 + DN(2,1)*clhs135 + DN(2,2)*clhs181 + clhs198 + clhs286);
            lhs(10,3)=clhs1*(clhs199*clhs72 + clhs200*clhs70 - clhs200 + clhs216*clhs81);
            lhs(10,4)=clhs1*(DN(2,0)*clhs49 + DN(2,1)*clhs139 + DN(2,2)*clhs183 + clhs241);
            lhs(10,5)=clhs1*(DN(2,0)*clhs61 + DN(2,1)*clhs143 + DN(2,2)*clhs186 + clhs256);
            lhs(10,6)=clhs1*(DN(2,0)*clhs67 + DN(2,1)*clhs147 + DN(2,2)*clhs188 + clhs268 + clhs287);
            lhs(10,7)=clhs1*(clhs192*clhs81 + clhs269*clhs72 + clhs270*clhs70 - clhs270);
            lhs(10,8)=clhs1*(DN(2,0)*clhs78 + DN(2,1)*clhs153 + DN(2,2)*clhs193 + clhs293);
            lhs(10,9)=clhs1*(DN(2,0)*clhs90 + DN(2,1)*clhs157 + DN(2,2)*clhs195 + clhs307);
            lhs(10,10)=clhs1*(DN(2,0)*clhs96 + DN(2,1)*clhs161 + DN(2,2)*clhs197 + clhs12*clhs313 + clhs290);
            lhs(10,11)=DN(2,2)*clhs296;
            lhs(10,12)=clhs1*(DN(2,0)*clhs105 + DN(2,1)*clhs167 + DN(2,2)*clhs202 + clhs315);
            lhs(10,13)=clhs1*(DN(2,0)*clhs117 + DN(2,1)*clhs171 + DN(2,2)*clhs204 + clhs316);
            lhs(10,14)=clhs1*(DN(2,0)*clhs123 + DN(2,1)*clhs175 + DN(2,2)*clhs206 + clhs300 + clhs317);
            lhs(10,15)=clhs1*(DN(3,2)*clhs295 + clhs318*clhs70 - clhs318 + clhs319*clhs72);
            lhs(11,0)=clhs1*(clhs100*clhs17 + clhs98);
            lhs(11,1)=clhs1*(clhs163 + clhs165*clhs17);
            lhs(11,2)=clhs1*(clhs17*clhs201 + clhs199);
            lhs(11,3)=clhs221;
            lhs(11,4)=clhs1*(clhs100*clhs54 + clhs242);
            lhs(11,5)=clhs1*(clhs165*clhs54 + clhs257);
            lhs(11,6)=clhs1*(clhs201*clhs54 + clhs269);
            lhs(11,7)=clhs281;
            lhs(11,8)=DN(2,0)*clhs321;
            lhs(11,9)=DN(2,1)*clhs321;
            lhs(11,10)=DN(2,2)*clhs321;
            lhs(11,11)=clhs1*(clhs212*clhs289 + clhs213*clhs288 + clhs213*clhs305 + clhs213*clhs313);
            lhs(11,12)=clhs1*(clhs100*clhs110 + clhs304);
            lhs(11,13)=clhs1*(clhs110*clhs165 + clhs312);
            lhs(11,14)=clhs1*(clhs110*clhs201 + clhs319);
            lhs(11,15)=clhs322;
            lhs(12,0)=clhs1*(DN(3,0)*clhs2 + DN(3,1)*clhs4 + DN(3,2)*clhs6 + clhs106 + clhs326);
            lhs(12,1)=clhs1*(DN(3,0)*clhs24 + DN(3,1)*clhs26 + DN(3,2)*clhs29 + clhs168);
            lhs(12,2)=clhs1*(DN(3,0)*clhs32 + DN(3,1)*clhs34 + DN(3,2)*clhs36 + clhs203);
            lhs(12,3)=clhs1*(clhs108*clhs214 + clhs125*clhs72 + clhs126*clhs70 - clhs126);
            lhs(12,4)=clhs1*(DN(3,0)*clhs45 + DN(3,1)*clhs47 + DN(3,2)*clhs49 + clhs244 + clhs327);
            lhs(12,5)=clhs1*(DN(3,0)*clhs56 + DN(3,1)*clhs58 + DN(3,2)*clhs61 + clhs259);
            lhs(12,6)=clhs1*(DN(3,0)*clhs63 + DN(3,1)*clhs65 + DN(3,2)*clhs67 + clhs271);
            lhs(12,7)=clhs1*(clhs108*clhs73 + clhs249*clhs72 + clhs250*clhs70 - clhs250);
            lhs(12,8)=clhs1*(DN(3,0)*clhs74 + DN(3,1)*clhs76 + DN(3,2)*clhs78 + clhs297 + clhs328);
            lhs(12,9)=clhs1*(DN(3,0)*clhs85 + DN(3,1)*clhs87 + DN(3,2)*clhs90 + clhs308);
            lhs(12,10)=clhs1*(DN(3,0)*clhs92 + DN(3,1)*clhs94 + DN(3,2)*clhs96 + clhs315);
            lhs(12,11)=clhs1*(clhs100*clhs108 + clhs303*clhs72 + clhs304*clhs70 - clhs304);
            lhs(12,12)=clhs1*(DN(3,0)*clhs101 + DN(3,1)*clhs103 + DN(3,2)*clhs105 + clhs12*clhs329 + clhs331);
            lhs(12,13)=clhs1*(DN(3,0)*clhs112 + DN(3,1)*clhs114 + DN(3,2)*clhs117 + clhs333);
            lhs(12,14)=clhs1*(DN(3,0)*clhs119 + DN(3,1)*clhs121 + DN(3,2)*clhs123 + clhs334);
            lhs(12,15)=DN(3,0)*clhs335;
            lhs(13,0)=clhs1*(DN(3,0)*clhs4 + DN(3,1)*clhs127 + DN(3,2)*clhs128 + clhs118);
            lhs(13,1)=clhs1*(DN(3,0)*clhs26 + DN(3,1)*clhs129 + DN(3,2)*clhs131 + clhs172 + clhs326);
            lhs(13,2)=clhs1*(DN(3,0)*clhs34 + DN(3,1)*clhs133 + DN(3,2)*clhs135 + clhs205);
            lhs(13,3)=clhs1*(clhs108*clhs215 + clhs177*clhs72 + clhs178*clhs70 - clhs178);
            lhs(13,4)=clhs1*(DN(3,0)*clhs47 + DN(3,1)*clhs138 + DN(3,2)*clhs139 + clhs247);
            lhs(13,5)=clhs1*(DN(3,0)*clhs58 + DN(3,1)*clhs141 + DN(3,2)*clhs143 + clhs260 + clhs327);
            lhs(13,6)=clhs1*(DN(3,0)*clhs65 + DN(3,1)*clhs145 + DN(3,2)*clhs147 + clhs272);
            lhs(13,7)=clhs1*(clhs108*clhs151 + clhs262*clhs72 + clhs263*clhs70 - clhs263);
            lhs(13,8)=clhs1*(DN(3,0)*clhs76 + DN(3,1)*clhs152 + DN(3,2)*clhs153 + clhs301);
            lhs(13,9)=clhs1*(DN(3,0)*clhs87 + DN(3,1)*clhs155 + DN(3,2)*clhs157 + clhs309 + clhs328);
            lhs(13,10)=clhs1*(DN(3,0)*clhs94 + DN(3,1)*clhs159 + DN(3,2)*clhs161 + clhs316);
            lhs(13,11)=clhs1*(clhs108*clhs165 + clhs311*clhs72 + clhs312*clhs70 - clhs312);
            lhs(13,12)=clhs1*(DN(3,0)*clhs103 + DN(3,1)*clhs166 + DN(3,2)*clhs167 + clhs333);
            lhs(13,13)=clhs1*(DN(3,0)*clhs114 + DN(3,1)*clhs169 + DN(3,2)*clhs171 + clhs12*clhs336 + clhs331);
            lhs(13,14)=clhs1*(DN(3,0)*clhs121 + DN(3,1)*clhs173 + DN(3,2)*clhs175 + clhs337);
            lhs(13,15)=DN(3,1)*clhs335;
            lhs(14,0)=clhs1*(DN(3,0)*clhs6 + DN(3,1)*clhs128 + DN(3,2)*clhs179 + clhs124);
            lhs(14,1)=clhs1*(DN(3,0)*clhs29 + DN(3,1)*clhs131 + DN(3,2)*clhs180 + clhs176);
            lhs(14,2)=clhs1*(DN(3,0)*clhs36 + DN(3,1)*clhs135 + DN(3,2)*clhs181 + clhs207 + clhs326);
            lhs(14,3)=clhs1*(clhs108*clhs216 + clhs208*clhs72 + clhs209*clhs70 - clhs209);
            lhs(14,4)=clhs1*(DN(3,0)*clhs49 + DN(3,1)*clhs139 + DN(3,2)*clhs183 + clhs248);
            lhs(14,5)=clhs1*(DN(3,0)*clhs61 + DN(3,1)*clhs143 + DN(3,2)*clhs186 + clhs261);
            lhs(14,6)=clhs1*(DN(3,0)*clhs67 + DN(3,1)*clhs147 + DN(3,2)*clhs188 + clhs273 + clhs327);
            lhs(14,7)=clhs1*(clhs108*clhs192 + clhs274*clhs72 + clhs275*clhs70 - clhs275);
            lhs(14,8)=clhs1*(DN(3,0)*clhs78 + DN(3,1)*clhs153 + DN(3,2)*clhs193 + clhs302);
            lhs(14,9)=clhs1*(DN(3,0)*clhs90 + DN(3,1)*clhs157 + DN(3,2)*clhs195 + clhs310);
            lhs(14,10)=clhs1*(DN(3,0)*clhs96 + DN(3,1)*clhs161 + DN(3,2)*clhs197 + clhs317 + clhs328);
            lhs(14,11)=clhs1*(clhs108*clhs201 + clhs318*clhs72 + clhs319*clhs70 - clhs319);
            lhs(14,12)=clhs1*(DN(3,0)*clhs105 + DN(3,1)*clhs167 + DN(3,2)*clhs202 + clhs334);
            lhs(14,13)=clhs1*(DN(3,0)*clhs117 + DN(3,1)*clhs171 + DN(3,2)*clhs204 + clhs337);
            lhs(14,14)=clhs1*(DN(3,0)*clhs123 + DN(3,1)*clhs175 + DN(3,2)*clhs206 + clhs12*clhs338 + clhs331);
            lhs(14,15)=DN(3,2)*clhs335;
            lhs(15,0)=clhs1*(DN(3,0)*clhs210 + clhs125);
            lhs(15,1)=clhs1*(DN(3,1)*clhs210 + clhs177);
            lhs(15,2)=clhs1*(DN(3,2)*clhs210 + clhs208);
            lhs(15,3)=clhs222;
            lhs(15,4)=clhs1*(DN(3,0)*clhs276 + clhs249);
            lhs(15,5)=clhs1*(DN(3,1)*clhs276 + clhs262);
            lhs(15,6)=clhs1*(DN(3,2)*clhs276 + clhs274);
            lhs(15,7)=clhs282;
            lhs(15,8)=clhs1*(DN(3,0)*clhs320 + clhs303);
            lhs(15,9)=clhs1*(DN(3,1)*clhs320 + clhs311);
            lhs(15,10)=clhs1*(DN(3,2)*clhs320 + clhs318);
            lhs(15,11)=clhs322;
            lhs(15,12)=DN(3,0)*clhs339;
            lhs(15,13)=DN(3,1)*clhs339;
            lhs(15,14)=DN(3,2)*clhs339;
            lhs(15,15)=clhs1*(clhs212*clhs330 + clhs213*clhs329 + clhs213*clhs336 + clhs213*clhs338);


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
    const double& dtn = data.dtn;
    const double& dtnn = data.dtnn;
    const double& dtnnn = data.dtnnn;

    const double& t = data.t;
    const double& tn = data.tn;
    const double& tnn = data.tnn;
    const double& tnnn = data.tnnn;
    const double& bdf0 = data.bdf0;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v_ave = data.v_ave;
    const BoundedMatrix<double,nnodes,dim>& vn_ave = data.vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;

    // Get time accurate expression
    const BoundedMatrix<double,nnodes,dim>& v = (t * v_ave - tn * vn_ave) * ( 1.0 / dt );
    const BoundedMatrix<double,nnodes,dim>& vconv = v - vmesh;

    // Get constitutive matrix
    const Matrix& C = data.C;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    const double clhs0 =             1.0/dt;
const double clhs1 =             clhs0*t;
const double clhs2 =             C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double clhs3 =             C(0,2)*DN(0,0);
const double clhs4 =             C(2,2)*DN(0,1) + clhs3;
const double clhs5 =             pow(DN(0,0), 2);
const double clhs6 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double clhs7 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double clhs8 =             rho*stab_c2*sqrt(pow(clhs6, 2) + pow(clhs7, 2));
const double clhs9 =             clhs8*h/stab_c1 + mu;
const double clhs10 =             pow(N[0], 2);
const double clhs11 =             bdf0*rho;
const double clhs12 =             N[0]*rho;
const double clhs13 =             DN(0,0)*clhs6 + DN(0,1)*clhs7;
const double clhs14 =             N[0]*bdf0 + clhs13;
const double clhs15 =             pow(rho, 2);
const double clhs16 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double clhs17 =             1.0/(clhs8/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs18 =             1.0*N[0]*clhs15*clhs16*clhs17;
const double clhs19 =             1.0*clhs13*clhs15*clhs17;
const double clhs20 =             clhs10*clhs11 + clhs12*clhs13 + clhs14*clhs18 + clhs14*clhs19;
const double clhs21 =             C(0,1)*DN(0,1) + clhs3;
const double clhs22 =             C(1,2)*DN(0,1);
const double clhs23 =             C(2,2)*DN(0,0) + clhs22;
const double clhs24 =             DN(0,0)*clhs9;
const double clhs25 =             DN(0,1)*clhs24;
const double clhs26 =             pow(c, -2);
const double clhs27 =             1.0/rho;
const double clhs28 =             N[0]*bdf0*clhs26*clhs27;
const double clhs29 =             1.0*clhs16*clhs17;
const double clhs30 =             1.0*clhs17*rho;
const double clhs31 =             clhs13*clhs30;
const double clhs32 =             clhs0*t*(-N[0] + clhs12*clhs29 + clhs28*clhs9 + clhs31);
const double clhs33 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs34 =             C(0,2)*DN(1,0);
const double clhs35 =             C(2,2)*DN(1,1) + clhs34;
const double clhs36 =             DN(1,0)*clhs24;
const double clhs37 =             N[0]*bdf0*rho;
const double clhs38 =             N[1]*clhs37;
const double clhs39 =             DN(1,0)*clhs6 + DN(1,1)*clhs7;
const double clhs40 =             N[1]*bdf0;
const double clhs41 =             clhs39 + clhs40;
const double clhs42 =             clhs12*clhs39 + clhs18*clhs41 + clhs19*clhs41 + clhs38;
const double clhs43 =             C(0,1)*DN(1,1) + clhs34;
const double clhs44 =             C(1,2)*DN(1,1);
const double clhs45 =             C(2,2)*DN(1,0) + clhs44;
const double clhs46 =             DN(1,1)*clhs24;
const double clhs47 =             DN(0,0)*N[1];
const double clhs48 =             bdf0*clhs26*clhs27*clhs9;
const double clhs49 =             DN(1,0)*N[0];
const double clhs50 =             1.0*clhs16*clhs17*rho;
const double clhs51 =             1.0*DN(1,0)*clhs17*rho;
const double clhs52 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs53 =             C(0,2)*DN(2,0);
const double clhs54 =             C(2,2)*DN(2,1) + clhs53;
const double clhs55 =             DN(2,0)*clhs24;
const double clhs56 =             N[2]*clhs37;
const double clhs57 =             DN(2,0)*clhs6 + DN(2,1)*clhs7;
const double clhs58 =             N[2]*bdf0;
const double clhs59 =             clhs57 + clhs58;
const double clhs60 =             clhs12*clhs57 + clhs18*clhs59 + clhs19*clhs59 + clhs56;
const double clhs61 =             C(0,1)*DN(2,1) + clhs53;
const double clhs62 =             C(1,2)*DN(2,1);
const double clhs63 =             C(2,2)*DN(2,0) + clhs62;
const double clhs64 =             DN(2,1)*clhs24;
const double clhs65 =             DN(0,0)*N[2];
const double clhs66 =             DN(2,0)*N[0];
const double clhs67 =             C(0,1)*DN(0,0) + clhs22;
const double clhs68 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs69 =             pow(DN(0,1), 2);
const double clhs70 =             C(0,1)*DN(1,0) + clhs44;
const double clhs71 =             DN(0,1)*clhs9;
const double clhs72 =             DN(1,0)*clhs71;
const double clhs73 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs74 =             DN(1,1)*clhs71;
const double clhs75 =             DN(0,1)*N[1];
const double clhs76 =             DN(1,1)*N[0];
const double clhs77 =             1.0*DN(1,1)*clhs17*rho;
const double clhs78 =             C(0,1)*DN(2,0) + clhs62;
const double clhs79 =             DN(2,0)*clhs71;
const double clhs80 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs81 =             DN(2,1)*clhs71;
const double clhs82 =             DN(0,1)*N[2];
const double clhs83 =             DN(2,1)*N[0];
const double clhs84 =             clhs14*clhs30;
const double clhs85 =             clhs0*t*(N[0] + clhs84);
const double clhs86 =             bdf0*clhs26*clhs27;
const double clhs87 =             1.0*clhs17;
const double clhs88 =             1.0*DN(0,0)*clhs17*rho;
const double clhs89 =             1.0*DN(0,1)*clhs17*rho;
const double clhs90 =             1.0*DN(0,0)*clhs17;
const double clhs91 =             1.0*DN(0,1)*clhs17;
const double clhs92 =             clhs1*(DN(1,0)*clhs90 + DN(1,1)*clhs91 + N[1]*clhs28);
const double clhs93 =             clhs1*(DN(2,0)*clhs90 + DN(2,1)*clhs91 + N[2]*clhs28);
const double clhs94 =             N[1]*rho;
const double clhs95 =             1.0*N[1]*clhs15*clhs16*clhs17;
const double clhs96 =             1.0*clhs15*clhs17*clhs39;
const double clhs97 =             clhs13*clhs94 + clhs14*clhs95 + clhs14*clhs96 + clhs38;
const double clhs98 =             pow(DN(1,0), 2);
const double clhs99 =             pow(N[1], 2);
const double clhs100 =             clhs11*clhs99 + clhs39*clhs94 + clhs41*clhs95 + clhs41*clhs96;
const double clhs101 =             DN(1,0)*clhs9;
const double clhs102 =             DN(1,1)*clhs101;
const double clhs103 =             clhs26*clhs27*clhs9;
const double clhs104 =             clhs30*clhs39;
const double clhs105 =             clhs0*t*(-N[1] + clhs103*clhs40 + clhs104 + clhs29*clhs94);
const double clhs106 =             DN(2,0)*clhs101;
const double clhs107 =             N[1]*N[2]*bdf0;
const double clhs108 =             clhs107*rho;
const double clhs109 =             clhs108 + clhs57*clhs94 + clhs59*clhs95 + clhs59*clhs96;
const double clhs110 =             DN(2,1)*clhs101;
const double clhs111 =             DN(1,0)*N[2];
const double clhs112 =             DN(2,0)*N[1];
const double clhs113 =             pow(DN(1,1), 2);
const double clhs114 =             DN(1,1)*clhs9;
const double clhs115 =             DN(2,0)*clhs114;
const double clhs116 =             DN(2,1)*clhs114;
const double clhs117 =             DN(1,1)*N[2];
const double clhs118 =             DN(2,1)*N[1];
const double clhs119 =             clhs30*clhs41;
const double clhs120 =             clhs0*t*(N[1] + clhs119);
const double clhs121 =             clhs1*(1.0*DN(1,0)*DN(2,0)*clhs17 + 1.0*DN(1,1)*DN(2,1)*clhs17 + clhs107*clhs26*clhs27);
const double clhs122 =             N[2]*rho;
const double clhs123 =             1.0*N[2]*clhs15*clhs16*clhs17;
const double clhs124 =             1.0*clhs15*clhs17*clhs57;
const double clhs125 =             clhs122*clhs13 + clhs123*clhs14 + clhs124*clhs14 + clhs56;
const double clhs126 =             clhs108 + clhs122*clhs39 + clhs123*clhs41 + clhs124*clhs41;
const double clhs127 =             pow(DN(2,0), 2);
const double clhs128 =             pow(N[2], 2);
const double clhs129 =             clhs11*clhs128 + clhs122*clhs57 + clhs123*clhs59 + clhs124*clhs59;
const double clhs130 =             DN(2,0)*DN(2,1)*clhs9;
const double clhs131 =             clhs0*t*(-N[2] + clhs103*clhs58 + clhs122*clhs29 + clhs30*clhs57);
const double clhs132 =             pow(DN(2,1), 2);
const double clhs133 =             clhs0*t*(N[2] + clhs30*clhs59);
            lhs(0,0)=clhs1*(DN(0,0)*clhs2 + DN(0,1)*clhs4 + clhs20 + clhs5*clhs9);
            lhs(0,1)=clhs1*(DN(0,0)*clhs21 + DN(0,1)*clhs23 + clhs25);
            lhs(0,2)=DN(0,0)*clhs32;
            lhs(0,3)=clhs1*(DN(0,0)*clhs33 + DN(0,1)*clhs35 + clhs36 + clhs42);
            lhs(0,4)=clhs1*(DN(0,0)*clhs43 + DN(0,1)*clhs45 + clhs46);
            lhs(0,5)=clhs1*(clhs13*clhs51 + clhs47*clhs48 - clhs47 + clhs49*clhs50);
            lhs(0,6)=clhs1*(DN(0,0)*clhs52 + DN(0,1)*clhs54 + clhs55 + clhs60);
            lhs(0,7)=clhs1*(DN(0,0)*clhs61 + DN(0,1)*clhs63 + clhs64);
            lhs(0,8)=clhs1*(DN(2,0)*clhs31 + clhs48*clhs65 + clhs50*clhs66 - clhs65);
            lhs(1,0)=clhs1*(DN(0,0)*clhs4 + DN(0,1)*clhs67 + clhs25);
            lhs(1,1)=clhs1*(DN(0,0)*clhs23 + DN(0,1)*clhs68 + clhs20 + clhs69*clhs9);
            lhs(1,2)=DN(0,1)*clhs32;
            lhs(1,3)=clhs1*(DN(0,0)*clhs35 + DN(0,1)*clhs70 + clhs72);
            lhs(1,4)=clhs1*(DN(0,0)*clhs45 + DN(0,1)*clhs73 + clhs42 + clhs74);
            lhs(1,5)=clhs1*(clhs13*clhs77 + clhs48*clhs75 + clhs50*clhs76 - clhs75);
            lhs(1,6)=clhs1*(DN(0,0)*clhs54 + DN(0,1)*clhs78 + clhs79);
            lhs(1,7)=clhs1*(DN(0,0)*clhs63 + DN(0,1)*clhs80 + clhs60 + clhs81);
            lhs(1,8)=clhs1*(DN(2,1)*clhs31 + clhs48*clhs82 + clhs50*clhs83 - clhs82);
            lhs(2,0)=DN(0,0)*clhs85;
            lhs(2,1)=DN(0,1)*clhs85;
            lhs(2,2)=clhs1*(clhs10*clhs86 + clhs5*clhs87 + clhs69*clhs87);
            lhs(2,3)=clhs1*(clhs41*clhs88 + clhs49);
            lhs(2,4)=clhs1*(clhs41*clhs89 + clhs76);
            lhs(2,5)=clhs92;
            lhs(2,6)=clhs1*(clhs59*clhs88 + clhs66);
            lhs(2,7)=clhs1*(clhs59*clhs89 + clhs83);
            lhs(2,8)=clhs93;
            lhs(3,0)=clhs1*(DN(1,0)*clhs2 + DN(1,1)*clhs4 + clhs36 + clhs97);
            lhs(3,1)=clhs1*(DN(1,0)*clhs21 + DN(1,1)*clhs23 + clhs72);
            lhs(3,2)=clhs1*(clhs39*clhs88 + clhs47*clhs50 + clhs48*clhs49 - clhs49);
            lhs(3,3)=clhs1*(DN(1,0)*clhs33 + DN(1,1)*clhs35 + clhs100 + clhs9*clhs98);
            lhs(3,4)=clhs1*(DN(1,0)*clhs43 + DN(1,1)*clhs45 + clhs102);
            lhs(3,5)=DN(1,0)*clhs105;
            lhs(3,6)=clhs1*(DN(1,0)*clhs52 + DN(1,1)*clhs54 + clhs106 + clhs109);
            lhs(3,7)=clhs1*(DN(1,0)*clhs61 + DN(1,1)*clhs63 + clhs110);
            lhs(3,8)=clhs1*(DN(2,0)*clhs104 + clhs111*clhs48 - clhs111 + clhs112*clhs50);
            lhs(4,0)=clhs1*(DN(1,0)*clhs4 + DN(1,1)*clhs67 + clhs46);
            lhs(4,1)=clhs1*(DN(1,0)*clhs23 + DN(1,1)*clhs68 + clhs74 + clhs97);
            lhs(4,2)=clhs1*(clhs39*clhs89 + clhs48*clhs76 + clhs50*clhs75 - clhs76);
            lhs(4,3)=clhs1*(DN(1,0)*clhs35 + DN(1,1)*clhs70 + clhs102);
            lhs(4,4)=clhs1*(DN(1,0)*clhs45 + DN(1,1)*clhs73 + clhs100 + clhs113*clhs9);
            lhs(4,5)=DN(1,1)*clhs105;
            lhs(4,6)=clhs1*(DN(1,0)*clhs54 + DN(1,1)*clhs78 + clhs115);
            lhs(4,7)=clhs1*(DN(1,0)*clhs63 + DN(1,1)*clhs80 + clhs109 + clhs116);
            lhs(4,8)=clhs1*(DN(2,1)*clhs104 + clhs117*clhs48 - clhs117 + clhs118*clhs50);
            lhs(5,0)=clhs1*(clhs14*clhs51 + clhs47);
            lhs(5,1)=clhs1*(clhs14*clhs77 + clhs75);
            lhs(5,2)=clhs92;
            lhs(5,3)=DN(1,0)*clhs120;
            lhs(5,4)=DN(1,1)*clhs120;
            lhs(5,5)=clhs1*(clhs113*clhs87 + clhs86*clhs99 + clhs87*clhs98);
            lhs(5,6)=clhs1*(clhs112 + clhs51*clhs59);
            lhs(5,7)=clhs1*(clhs118 + clhs59*clhs77);
            lhs(5,8)=clhs121;
            lhs(6,0)=clhs1*(DN(2,0)*clhs2 + DN(2,1)*clhs4 + clhs125 + clhs55);
            lhs(6,1)=clhs1*(DN(2,0)*clhs21 + DN(2,1)*clhs23 + clhs79);
            lhs(6,2)=clhs1*(clhs48*clhs66 + clhs50*clhs65 + clhs57*clhs88 - clhs66);
            lhs(6,3)=clhs1*(DN(2,0)*clhs33 + DN(2,1)*clhs35 + clhs106 + clhs126);
            lhs(6,4)=clhs1*(DN(2,0)*clhs43 + DN(2,1)*clhs45 + clhs115);
            lhs(6,5)=clhs1*(clhs111*clhs50 + clhs112*clhs48 - clhs112 + clhs51*clhs57);
            lhs(6,6)=clhs1*(DN(2,0)*clhs52 + DN(2,1)*clhs54 + clhs127*clhs9 + clhs129);
            lhs(6,7)=clhs1*(DN(2,0)*clhs61 + DN(2,1)*clhs63 + clhs130);
            lhs(6,8)=DN(2,0)*clhs131;
            lhs(7,0)=clhs1*(DN(2,0)*clhs4 + DN(2,1)*clhs67 + clhs64);
            lhs(7,1)=clhs1*(DN(2,0)*clhs23 + DN(2,1)*clhs68 + clhs125 + clhs81);
            lhs(7,2)=clhs1*(clhs48*clhs83 + clhs50*clhs82 + clhs57*clhs89 - clhs83);
            lhs(7,3)=clhs1*(DN(2,0)*clhs35 + DN(2,1)*clhs70 + clhs110);
            lhs(7,4)=clhs1*(DN(2,0)*clhs45 + DN(2,1)*clhs73 + clhs116 + clhs126);
            lhs(7,5)=clhs1*(clhs117*clhs50 + clhs118*clhs48 - clhs118 + clhs57*clhs77);
            lhs(7,6)=clhs1*(DN(2,0)*clhs54 + DN(2,1)*clhs78 + clhs130);
            lhs(7,7)=clhs1*(DN(2,0)*clhs63 + DN(2,1)*clhs80 + clhs129 + clhs132*clhs9);
            lhs(7,8)=DN(2,1)*clhs131;
            lhs(8,0)=clhs1*(DN(2,0)*clhs84 + clhs65);
            lhs(8,1)=clhs1*(DN(2,1)*clhs84 + clhs82);
            lhs(8,2)=clhs93;
            lhs(8,3)=clhs1*(DN(2,0)*clhs119 + clhs111);
            lhs(8,4)=clhs1*(DN(2,1)*clhs119 + clhs117);
            lhs(8,5)=clhs121;
            lhs(8,6)=DN(2,0)*clhs133;
            lhs(8,7)=DN(2,1)*clhs133;
            lhs(8,8)=clhs1*(clhs127*clhs87 + clhs128*clhs86 + clhs132*clhs87);


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
    const double& dtn = data.dtn;
    const double& dtnn = data.dtnn;
    const double& dtnnn = data.dtnnn;

    const double& t = data.t;
    const double& tn = data.tn;
    const double& tnn = data.tnn;
    const double& tnnn = data.tnnn;

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

    // Get time accurate expression
    const array_1d<double,nnodes>& p = (t * p_ave - tn * pn_ave) * ( 1.0 / dt );
    const BoundedMatrix<double,nnodes,dim>& v = (t * v_ave - tn * vn_ave) * ( 1.0 / dt );
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

    const double crhs0 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crhs1 =             1.0/dt;
const double crhs2 =             pn_ave[0]*tn;
const double crhs3 =             -crhs2 + p_ave[0]*t;
const double crhs4 =             pn_ave[1]*tn;
const double crhs5 =             -crhs4 + p_ave[1]*t;
const double crhs6 =             pn_ave[2]*tn;
const double crhs7 =             -crhs6 + p_ave[2]*t;
const double crhs8 =             pn_ave[3]*tn;
const double crhs9 =             -crhs8 + p_ave[3]*t;
const double crhs10 =             crhs1*(N[0]*crhs3 + N[1]*crhs5 + N[2]*crhs7 + N[3]*crhs9);
const double crhs11 =             bdf0*crhs1;
const double crhs12 =             tn*vn_ave(0,0);
const double crhs13 =             -crhs12 + t*v_ave(0,0);
const double crhs14 =             bdf1/dtn;
const double crhs15 =             tnn*vnn_ave(0,0);
const double crhs16 =             bdf2/dtnn;
const double crhs17 =             tn*vn_ave(1,0);
const double crhs18 =             -crhs17 + t*v_ave(1,0);
const double crhs19 =             tnn*vnn_ave(1,0);
const double crhs20 =             tn*vn_ave(2,0);
const double crhs21 =             -crhs20 + t*v_ave(2,0);
const double crhs22 =             tnn*vnn_ave(2,0);
const double crhs23 =             tn*vn_ave(3,0);
const double crhs24 =             -crhs23 + t*v_ave(3,0);
const double crhs25 =             tnn*vnn_ave(3,0);
const double crhs26 =             rho*(N[0]*(crhs11*crhs13 + crhs14*(crhs12 - crhs15) + crhs16*(crhs15 - tnnn*vnnn_ave(0,0))) + N[1]*(crhs11*crhs18 + crhs14*(crhs17 - crhs19) + crhs16*(crhs19 - tnnn*vnnn_ave(1,0))) + N[2]*(crhs11*crhs21 + crhs14*(crhs20 - crhs22) + crhs16*(crhs22 - tnnn*vnnn_ave(2,0))) + N[3]*(crhs11*crhs24 + crhs14*(crhs23 - crhs25) + crhs16*(crhs25 - tnnn*vnnn_ave(3,0))));
const double crhs27 =             crhs1*rho;
const double crhs28 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs29 =             DN(0,0)*crhs13;
const double crhs30 =             DN(1,0)*crhs18;
const double crhs31 =             DN(2,0)*crhs21;
const double crhs32 =             DN(3,0)*crhs24;
const double crhs33 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs34 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs35 =             crhs27*(crhs28*(crhs29 + crhs30 + crhs31 + crhs32) + crhs33*(DN(0,1)*crhs13 + DN(1,1)*crhs18 + DN(2,1)*crhs21 + DN(3,1)*crhs24) + crhs34*(DN(0,2)*crhs13 + DN(1,2)*crhs18 + DN(2,2)*crhs21 + DN(3,2)*crhs24));
const double crhs36 =             rho*stab_c2*sqrt(pow(crhs28, 2) + pow(crhs33, 2) + pow(crhs34, 2));
const double crhs37 =             tn*vn_ave(0,2);
const double crhs38 =             -crhs37 + t*v_ave(0,2);
const double crhs39 =             tn*vn_ave(1,2);
const double crhs40 =             -crhs39 + t*v_ave(1,2);
const double crhs41 =             tn*vn_ave(2,2);
const double crhs42 =             -crhs41 + t*v_ave(2,2);
const double crhs43 =             tn*vn_ave(3,2);
const double crhs44 =             -crhs43 + t*v_ave(3,2);
const double crhs45 =             DN(0,2)*crhs38 + DN(1,2)*crhs40 + DN(2,2)*crhs42 + DN(3,2)*crhs44;
const double crhs46 =             tn*vn_ave(0,1);
const double crhs47 =             -crhs46 + t*v_ave(0,1);
const double crhs48 =             DN(0,1)*crhs47;
const double crhs49 =             tn*vn_ave(1,1);
const double crhs50 =             -crhs49 + t*v_ave(1,1);
const double crhs51 =             DN(1,1)*crhs50;
const double crhs52 =             tn*vn_ave(2,1);
const double crhs53 =             -crhs52 + t*v_ave(2,1);
const double crhs54 =             DN(2,1)*crhs53;
const double crhs55 =             tn*vn_ave(3,1);
const double crhs56 =             -crhs55 + t*v_ave(3,1);
const double crhs57 =             DN(3,1)*crhs56;
const double crhs58 =             crhs1*(crhs29 + crhs30 + crhs31 + crhs32 + crhs45 + crhs48 + crhs51 + crhs54 + crhs57);
const double crhs59 =             pnn_ave[0]*tnn;
const double crhs60 =             pnn_ave[1]*tnn;
const double crhs61 =             pnn_ave[2]*tnn;
const double crhs62 =             pnn_ave[3]*tnn;
const double crhs63 =             (N[0]*(crhs11*crhs3 + crhs14*(crhs2 - crhs59) + crhs16*(crhs59 - pnnn_ave[0]*tnnn)) + N[1]*(crhs11*crhs5 + crhs14*(crhs4 - crhs60) + crhs16*(crhs60 - pnnn_ave[1]*tnnn)) + N[2]*(crhs11*crhs7 + crhs14*(crhs6 - crhs61) + crhs16*(crhs61 - pnnn_ave[2]*tnnn)) + N[3]*(crhs11*crhs9 + crhs14*(-crhs62 + crhs8) + crhs16*(crhs62 - pnnn_ave[3]*tnnn)))/(pow(c, 2)*rho);
const double crhs64 =             (crhs58 + crhs63)*(crhs36*h/stab_c1 + mu);
const double crhs65 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crhs66 =             N[0]*crhs65*rho;
const double crhs67 =             1.0/(crhs36/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs68 =             1.0*crhs67*(-crhs0 + crhs1*(DN(0,0)*crhs3 + DN(1,0)*crhs5 + DN(2,0)*crhs7 + DN(3,0)*crhs9) + crhs26 + crhs35);
const double crhs69 =             rho*(DN(0,0)*crhs28 + DN(0,1)*crhs33 + DN(0,2)*crhs34);
const double crhs70 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs71 =             tnn*vnn_ave(0,1);
const double crhs72 =             tnn*vnn_ave(1,1);
const double crhs73 =             tnn*vnn_ave(2,1);
const double crhs74 =             tnn*vnn_ave(3,1);
const double crhs75 =             rho*(N[0]*(crhs11*crhs47 + crhs14*(crhs46 - crhs71) + crhs16*(crhs71 - tnnn*vnnn_ave(0,1))) + N[1]*(crhs11*crhs50 + crhs14*(crhs49 - crhs72) + crhs16*(crhs72 - tnnn*vnnn_ave(1,1))) + N[2]*(crhs11*crhs53 + crhs14*(crhs52 - crhs73) + crhs16*(crhs73 - tnnn*vnnn_ave(2,1))) + N[3]*(crhs11*crhs56 + crhs14*(crhs55 - crhs74) + crhs16*(crhs74 - tnnn*vnnn_ave(3,1))));
const double crhs76 =             crhs27*(crhs28*(DN(0,0)*crhs47 + DN(1,0)*crhs50 + DN(2,0)*crhs53 + DN(3,0)*crhs56) + crhs33*(crhs48 + crhs51 + crhs54 + crhs57) + crhs34*(DN(0,2)*crhs47 + DN(1,2)*crhs50 + DN(2,2)*crhs53 + DN(3,2)*crhs56));
const double crhs77 =             1.0*crhs67*(crhs1*(DN(0,1)*crhs3 + DN(1,1)*crhs5 + DN(2,1)*crhs7 + DN(3,1)*crhs9) - crhs70 + crhs75 + crhs76);
const double crhs78 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs79 =             tnn*vnn_ave(0,2);
const double crhs80 =             tnn*vnn_ave(1,2);
const double crhs81 =             tnn*vnn_ave(2,2);
const double crhs82 =             tnn*vnn_ave(3,2);
const double crhs83 =             rho*(N[0]*(crhs11*crhs38 + crhs14*(crhs37 - crhs79) + crhs16*(crhs79 - tnnn*vnnn_ave(0,2))) + N[1]*(crhs11*crhs40 + crhs14*(crhs39 - crhs80) + crhs16*(crhs80 - tnnn*vnnn_ave(1,2))) + N[2]*(crhs11*crhs42 + crhs14*(crhs41 - crhs81) + crhs16*(crhs81 - tnnn*vnnn_ave(2,2))) + N[3]*(crhs11*crhs44 + crhs14*(crhs43 - crhs82) + crhs16*(crhs82 - tnnn*vnnn_ave(3,2))));
const double crhs84 =             crhs27*(crhs28*(DN(0,0)*crhs38 + DN(1,0)*crhs40 + DN(2,0)*crhs42 + DN(3,0)*crhs44) + crhs33*(DN(0,1)*crhs38 + DN(1,1)*crhs40 + DN(2,1)*crhs42 + DN(3,1)*crhs44) + crhs34*crhs45);
const double crhs85 =             1.0*crhs67*(crhs1*(DN(0,2)*crhs3 + DN(1,2)*crhs5 + DN(2,2)*crhs7 + DN(3,2)*crhs9) - crhs78 + crhs83 + crhs84);
const double crhs86 =             N[1]*crhs65*rho;
const double crhs87 =             rho*(DN(1,0)*crhs28 + DN(1,1)*crhs33 + DN(1,2)*crhs34);
const double crhs88 =             N[2]*crhs65*rho;
const double crhs89 =             rho*(DN(2,0)*crhs28 + DN(2,1)*crhs33 + DN(2,2)*crhs34);
const double crhs90 =             N[3]*crhs65*rho;
const double crhs91 =             rho*(DN(3,0)*crhs28 + DN(3,1)*crhs33 + DN(3,2)*crhs34);
            rhs[0]=DN(0,0)*crhs10 - DN(0,0)*crhs64 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs0 - N[0]*crhs26 - N[0]*crhs35 - crhs66*crhs68 - crhs68*crhs69;
            rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs10 - DN(0,1)*crhs64 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs70 - N[0]*crhs75 - N[0]*crhs76 - crhs66*crhs77 - crhs69*crhs77;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs10 - DN(0,2)*crhs64 - DN(0,2)*stress[2] + N[0]*crhs78 - N[0]*crhs83 - N[0]*crhs84 - crhs66*crhs85 - crhs69*crhs85;
            rhs[3]=-DN(0,0)*crhs68 - DN(0,1)*crhs77 - DN(0,2)*crhs85 - N[0]*crhs58 - N[0]*crhs63;
            rhs[4]=DN(1,0)*crhs10 - DN(1,0)*crhs64 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs0 - N[1]*crhs26 - N[1]*crhs35 - crhs68*crhs86 - crhs68*crhs87;
            rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs10 - DN(1,1)*crhs64 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs70 - N[1]*crhs75 - N[1]*crhs76 - crhs77*crhs86 - crhs77*crhs87;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs10 - DN(1,2)*crhs64 - DN(1,2)*stress[2] + N[1]*crhs78 - N[1]*crhs83 - N[1]*crhs84 - crhs85*crhs86 - crhs85*crhs87;
            rhs[7]=-DN(1,0)*crhs68 - DN(1,1)*crhs77 - DN(1,2)*crhs85 - N[1]*crhs58 - N[1]*crhs63;
            rhs[8]=DN(2,0)*crhs10 - DN(2,0)*crhs64 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs0 - N[2]*crhs26 - N[2]*crhs35 - crhs68*crhs88 - crhs68*crhs89;
            rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs10 - DN(2,1)*crhs64 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs70 - N[2]*crhs75 - N[2]*crhs76 - crhs77*crhs88 - crhs77*crhs89;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs10 - DN(2,2)*crhs64 - DN(2,2)*stress[2] + N[2]*crhs78 - N[2]*crhs83 - N[2]*crhs84 - crhs85*crhs88 - crhs85*crhs89;
            rhs[11]=-DN(2,0)*crhs68 - DN(2,1)*crhs77 - DN(2,2)*crhs85 - N[2]*crhs58 - N[2]*crhs63;
            rhs[12]=DN(3,0)*crhs10 - DN(3,0)*crhs64 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs0 - N[3]*crhs26 - N[3]*crhs35 - crhs68*crhs90 - crhs68*crhs91;
            rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs10 - DN(3,1)*crhs64 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs70 - N[3]*crhs75 - N[3]*crhs76 - crhs77*crhs90 - crhs77*crhs91;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs10 - DN(3,2)*crhs64 - DN(3,2)*stress[2] + N[3]*crhs78 - N[3]*crhs83 - N[3]*crhs84 - crhs85*crhs90 - crhs85*crhs91;
            rhs[15]=-DN(3,0)*crhs68 - DN(3,1)*crhs77 - DN(3,2)*crhs85 - N[3]*crhs58 - N[3]*crhs63;


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
    const double& dtn = data.dtn;
    const double& dtnn = data.dtnn;
    const double& dtnnn = data.dtnnn;

    const double& t = data.t;
    const double& tn = data.tn;
    const double& tnn = data.tnn;
    const double& tnnn = data.tnnn;

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

    // Get time accurate expression
    const array_1d<double,nnodes>& p = (t * p_ave - tn * pn_ave) * ( 1.0 / dt );
    const BoundedMatrix<double,nnodes,dim>& v = (t * v_ave - tn * vn_ave) * ( 1.0 / dt );
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
const double crhs1 =             1.0/dt;
const double crhs2 =             pn_ave[0]*tn;
const double crhs3 =             -crhs2 + p_ave[0]*t;
const double crhs4 =             pn_ave[1]*tn;
const double crhs5 =             -crhs4 + p_ave[1]*t;
const double crhs6 =             pn_ave[2]*tn;
const double crhs7 =             -crhs6 + p_ave[2]*t;
const double crhs8 =             crhs1*(N[0]*crhs3 + N[1]*crhs5 + N[2]*crhs7);
const double crhs9 =             crhs1*rho;
const double crhs10 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs11 =             tn*vn_ave(0,0);
const double crhs12 =             -crhs11 + t*v_ave(0,0);
const double crhs13 =             tn*vn_ave(1,0);
const double crhs14 =             -crhs13 + t*v_ave(1,0);
const double crhs15 =             tn*vn_ave(2,0);
const double crhs16 =             -crhs15 + t*v_ave(2,0);
const double crhs17 =             DN(0,0)*crhs12 + DN(1,0)*crhs14 + DN(2,0)*crhs16;
const double crhs18 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs19 =             crhs9*(crhs10*crhs17 + crhs18*(DN(0,1)*crhs12 + DN(1,1)*crhs14 + DN(2,1)*crhs16));
const double crhs20 =             bdf0*crhs1;
const double crhs21 =             bdf1/dtn;
const double crhs22 =             tnn*vnn_ave(0,0);
const double crhs23 =             bdf2/dtnn;
const double crhs24 =             tnn*vnn_ave(1,0);
const double crhs25 =             tnn*vnn_ave(2,0);
const double crhs26 =             rho*(N[0]*(crhs12*crhs20 + crhs21*(crhs11 - crhs22) + crhs23*(crhs22 - tnnn*vnnn_ave(0,0))) + N[1]*(crhs14*crhs20 + crhs21*(crhs13 - crhs24) + crhs23*(crhs24 - tnnn*vnnn_ave(1,0))) + N[2]*(crhs16*crhs20 + crhs21*(crhs15 - crhs25) + crhs23*(crhs25 - tnnn*vnnn_ave(2,0))));
const double crhs27 =             rho*stab_c2*sqrt(pow(crhs10, 2) + pow(crhs18, 2));
const double crhs28 =             tn*vn_ave(0,1);
const double crhs29 =             -crhs28 + t*v_ave(0,1);
const double crhs30 =             tn*vn_ave(1,1);
const double crhs31 =             -crhs30 + t*v_ave(1,1);
const double crhs32 =             tn*vn_ave(2,1);
const double crhs33 =             -crhs32 + t*v_ave(2,1);
const double crhs34 =             DN(0,1)*crhs29 + DN(1,1)*crhs31 + DN(2,1)*crhs33;
const double crhs35 =             crhs1*(crhs17 + crhs34);
const double crhs36 =             pnn_ave[0]*tnn;
const double crhs37 =             pnn_ave[1]*tnn;
const double crhs38 =             pnn_ave[2]*tnn;
const double crhs39 =             (N[0]*(crhs20*crhs3 + crhs21*(crhs2 - crhs36) + crhs23*(crhs36 - pnnn_ave[0]*tnnn)) + N[1]*(crhs20*crhs5 + crhs21*(-crhs37 + crhs4) + crhs23*(crhs37 - pnnn_ave[1]*tnnn)) + N[2]*(crhs20*crhs7 + crhs21*(-crhs38 + crhs6) + crhs23*(crhs38 - pnnn_ave[2]*tnnn)))/(pow(c, 2)*rho);
const double crhs40 =             (crhs35 + crhs39)*(crhs27*h/stab_c1 + mu);
const double crhs41 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crhs42 =             N[0]*crhs41*rho;
const double crhs43 =             1.0/(crhs27/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs44 =             1.0*crhs43*(-crhs0 + crhs1*(DN(0,0)*crhs3 + DN(1,0)*crhs5 + DN(2,0)*crhs7) + crhs19 + crhs26);
const double crhs45 =             rho*(DN(0,0)*crhs10 + DN(0,1)*crhs18);
const double crhs46 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs47 =             crhs9*(crhs10*(DN(0,0)*crhs29 + DN(1,0)*crhs31 + DN(2,0)*crhs33) + crhs18*crhs34);
const double crhs48 =             tnn*vnn_ave(0,1);
const double crhs49 =             tnn*vnn_ave(1,1);
const double crhs50 =             tnn*vnn_ave(2,1);
const double crhs51 =             rho*(N[0]*(crhs20*crhs29 + crhs21*(crhs28 - crhs48) + crhs23*(crhs48 - tnnn*vnnn_ave(0,1))) + N[1]*(crhs20*crhs31 + crhs21*(crhs30 - crhs49) + crhs23*(crhs49 - tnnn*vnnn_ave(1,1))) + N[2]*(crhs20*crhs33 + crhs21*(crhs32 - crhs50) + crhs23*(crhs50 - tnnn*vnnn_ave(2,1))));
const double crhs52 =             1.0*crhs43*(crhs1*(DN(0,1)*crhs3 + DN(1,1)*crhs5 + DN(2,1)*crhs7) - crhs46 + crhs47 + crhs51);
const double crhs53 =             N[1]*crhs41*rho;
const double crhs54 =             rho*(DN(1,0)*crhs10 + DN(1,1)*crhs18);
const double crhs55 =             N[2]*crhs41*rho;
const double crhs56 =             rho*(DN(2,0)*crhs10 + DN(2,1)*crhs18);
            rhs[0]=-DN(0,0)*crhs40 + DN(0,0)*crhs8 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs0 - N[0]*crhs19 - N[0]*crhs26 - crhs42*crhs44 - crhs44*crhs45;
            rhs[1]=-DN(0,0)*stress[2] - DN(0,1)*crhs40 + DN(0,1)*crhs8 - DN(0,1)*stress[1] + N[0]*crhs46 - N[0]*crhs47 - N[0]*crhs51 - crhs42*crhs52 - crhs45*crhs52;
            rhs[2]=-DN(0,0)*crhs44 - DN(0,1)*crhs52 - N[0]*crhs35 - N[0]*crhs39;
            rhs[3]=-DN(1,0)*crhs40 + DN(1,0)*crhs8 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs0 - N[1]*crhs19 - N[1]*crhs26 - crhs44*crhs53 - crhs44*crhs54;
            rhs[4]=-DN(1,0)*stress[2] - DN(1,1)*crhs40 + DN(1,1)*crhs8 - DN(1,1)*stress[1] + N[1]*crhs46 - N[1]*crhs47 - N[1]*crhs51 - crhs52*crhs53 - crhs52*crhs54;
            rhs[5]=-DN(1,0)*crhs44 - DN(1,1)*crhs52 - N[1]*crhs35 - N[1]*crhs39;
            rhs[6]=-DN(2,0)*crhs40 + DN(2,0)*crhs8 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs0 - N[2]*crhs19 - N[2]*crhs26 - crhs44*crhs55 - crhs44*crhs56;
            rhs[7]=-DN(2,0)*stress[2] - DN(2,1)*crhs40 + DN(2,1)*crhs8 - DN(2,1)*stress[1] + N[2]*crhs46 - N[2]*crhs47 - N[2]*crhs51 - crhs52*crhs55 - crhs52*crhs56;
            rhs[8]=-DN(2,0)*crhs44 - DN(2,1)*crhs52 - N[2]*crhs35 - N[2]*crhs39;


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
    const double& dtn = data.dtn;
    const double& dtnn = data.dtnn;
    const double& dtnnn = data.dtnnn;

    const double& t = data.t;
    const double& tn = data.tn;
    const double& tnn = data.tnn;
    const double& tnnn = data.tnnn;

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

    // Get time accurate expression
    const array_1d<double,nnodes>& p = (t * p_ave - tn * pn_ave) * ( 1.0 / dt );
    const BoundedMatrix<double,nnodes,dim>& v = (t * v_ave - tn * vn_ave) * ( 1.0 / dt );
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
    const double cv_s_gauss0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double cv_s_gauss1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double cv_s_gauss2 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double cv_s_gauss3 =             1.0/(rho*stab_c2*sqrt(pow(cv_s_gauss0, 2) + pow(cv_s_gauss1, 2) + pow(cv_s_gauss2, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cv_s_gauss4 =             1.0/dt;
const double cv_s_gauss5 =             p_ave[0]*t - pn_ave[0]*tn;
const double cv_s_gauss6 =             p_ave[1]*t - pn_ave[1]*tn;
const double cv_s_gauss7 =             p_ave[2]*t - pn_ave[2]*tn;
const double cv_s_gauss8 =             p_ave[3]*t - pn_ave[3]*tn;
const double cv_s_gauss9 =             bdf0*cv_s_gauss4;
const double cv_s_gauss10 =             tn*vn_ave(0,0);
const double cv_s_gauss11 =             -cv_s_gauss10 + t*v_ave(0,0);
const double cv_s_gauss12 =             bdf1/dtn;
const double cv_s_gauss13 =             tnn*vnn_ave(0,0);
const double cv_s_gauss14 =             bdf2/dtnn;
const double cv_s_gauss15 =             tn*vn_ave(1,0);
const double cv_s_gauss16 =             -cv_s_gauss15 + t*v_ave(1,0);
const double cv_s_gauss17 =             tnn*vnn_ave(1,0);
const double cv_s_gauss18 =             tn*vn_ave(2,0);
const double cv_s_gauss19 =             -cv_s_gauss18 + t*v_ave(2,0);
const double cv_s_gauss20 =             tnn*vnn_ave(2,0);
const double cv_s_gauss21 =             tn*vn_ave(3,0);
const double cv_s_gauss22 =             -cv_s_gauss21 + t*v_ave(3,0);
const double cv_s_gauss23 =             tnn*vnn_ave(3,0);
const double cv_s_gauss24 =             cv_s_gauss0*cv_s_gauss4;
const double cv_s_gauss25 =             cv_s_gauss1*cv_s_gauss4;
const double cv_s_gauss26 =             cv_s_gauss2*cv_s_gauss4;
const double cv_s_gauss27 =             tn*vn_ave(0,1);
const double cv_s_gauss28 =             -cv_s_gauss27 + t*v_ave(0,1);
const double cv_s_gauss29 =             tnn*vnn_ave(0,1);
const double cv_s_gauss30 =             tn*vn_ave(1,1);
const double cv_s_gauss31 =             -cv_s_gauss30 + t*v_ave(1,1);
const double cv_s_gauss32 =             tnn*vnn_ave(1,1);
const double cv_s_gauss33 =             tn*vn_ave(2,1);
const double cv_s_gauss34 =             -cv_s_gauss33 + t*v_ave(2,1);
const double cv_s_gauss35 =             tnn*vnn_ave(2,1);
const double cv_s_gauss36 =             tn*vn_ave(3,1);
const double cv_s_gauss37 =             -cv_s_gauss36 + t*v_ave(3,1);
const double cv_s_gauss38 =             tnn*vnn_ave(3,1);
const double cv_s_gauss39 =             tn*vn_ave(0,2);
const double cv_s_gauss40 =             -cv_s_gauss39 + t*v_ave(0,2);
const double cv_s_gauss41 =             tnn*vnn_ave(0,2);
const double cv_s_gauss42 =             tn*vn_ave(1,2);
const double cv_s_gauss43 =             -cv_s_gauss42 + t*v_ave(1,2);
const double cv_s_gauss44 =             tnn*vnn_ave(1,2);
const double cv_s_gauss45 =             tn*vn_ave(2,2);
const double cv_s_gauss46 =             -cv_s_gauss45 + t*v_ave(2,2);
const double cv_s_gauss47 =             tnn*vnn_ave(2,2);
const double cv_s_gauss48 =             tn*vn_ave(3,2);
const double cv_s_gauss49 =             -cv_s_gauss48 + t*v_ave(3,2);
const double cv_s_gauss50 =             tnn*vnn_ave(3,2);
            v_s_gauss[0]=-cv_s_gauss3*(cv_s_gauss4*(DN(0,0)*cv_s_gauss5 + DN(1,0)*cv_s_gauss6 + DN(2,0)*cv_s_gauss7 + DN(3,0)*cv_s_gauss8) + rho*(-N[0]*f(0,0) + N[0]*(cv_s_gauss11*cv_s_gauss9 + cv_s_gauss12*(cv_s_gauss10 - cv_s_gauss13) + cv_s_gauss14*(cv_s_gauss13 - tnnn*vnnn_ave(0,0))) - N[1]*f(1,0) + N[1]*(cv_s_gauss12*(cv_s_gauss15 - cv_s_gauss17) + cv_s_gauss14*(cv_s_gauss17 - tnnn*vnnn_ave(1,0)) + cv_s_gauss16*cv_s_gauss9) - N[2]*f(2,0) + N[2]*(cv_s_gauss12*(cv_s_gauss18 - cv_s_gauss20) + cv_s_gauss14*(cv_s_gauss20 - tnnn*vnnn_ave(2,0)) + cv_s_gauss19*cv_s_gauss9) - N[3]*f(3,0) + N[3]*(cv_s_gauss12*(cv_s_gauss21 - cv_s_gauss23) + cv_s_gauss14*(cv_s_gauss23 - tnnn*vnnn_ave(3,0)) + cv_s_gauss22*cv_s_gauss9) + cv_s_gauss24*(DN(0,0)*cv_s_gauss11 + DN(1,0)*cv_s_gauss16 + DN(2,0)*cv_s_gauss19 + DN(3,0)*cv_s_gauss22) + cv_s_gauss25*(DN(0,1)*cv_s_gauss11 + DN(1,1)*cv_s_gauss16 + DN(2,1)*cv_s_gauss19 + DN(3,1)*cv_s_gauss22) + cv_s_gauss26*(DN(0,2)*cv_s_gauss11 + DN(1,2)*cv_s_gauss16 + DN(2,2)*cv_s_gauss19 + DN(3,2)*cv_s_gauss22)));
            v_s_gauss[1]=-cv_s_gauss3*(cv_s_gauss4*(DN(0,1)*cv_s_gauss5 + DN(1,1)*cv_s_gauss6 + DN(2,1)*cv_s_gauss7 + DN(3,1)*cv_s_gauss8) + rho*(-N[0]*f(0,1) + N[0]*(cv_s_gauss12*(cv_s_gauss27 - cv_s_gauss29) + cv_s_gauss14*(cv_s_gauss29 - tnnn*vnnn_ave(0,1)) + cv_s_gauss28*cv_s_gauss9) - N[1]*f(1,1) + N[1]*(cv_s_gauss12*(cv_s_gauss30 - cv_s_gauss32) + cv_s_gauss14*(cv_s_gauss32 - tnnn*vnnn_ave(1,1)) + cv_s_gauss31*cv_s_gauss9) - N[2]*f(2,1) + N[2]*(cv_s_gauss12*(cv_s_gauss33 - cv_s_gauss35) + cv_s_gauss14*(cv_s_gauss35 - tnnn*vnnn_ave(2,1)) + cv_s_gauss34*cv_s_gauss9) - N[3]*f(3,1) + N[3]*(cv_s_gauss12*(cv_s_gauss36 - cv_s_gauss38) + cv_s_gauss14*(cv_s_gauss38 - tnnn*vnnn_ave(3,1)) + cv_s_gauss37*cv_s_gauss9) + cv_s_gauss24*(DN(0,0)*cv_s_gauss28 + DN(1,0)*cv_s_gauss31 + DN(2,0)*cv_s_gauss34 + DN(3,0)*cv_s_gauss37) + cv_s_gauss25*(DN(0,1)*cv_s_gauss28 + DN(1,1)*cv_s_gauss31 + DN(2,1)*cv_s_gauss34 + DN(3,1)*cv_s_gauss37) + cv_s_gauss26*(DN(0,2)*cv_s_gauss28 + DN(1,2)*cv_s_gauss31 + DN(2,2)*cv_s_gauss34 + DN(3,2)*cv_s_gauss37)));
            v_s_gauss[2]=-cv_s_gauss3*(cv_s_gauss4*(DN(0,2)*cv_s_gauss5 + DN(1,2)*cv_s_gauss6 + DN(2,2)*cv_s_gauss7 + DN(3,2)*cv_s_gauss8) + rho*(-N[0]*f(0,2) + N[0]*(cv_s_gauss12*(cv_s_gauss39 - cv_s_gauss41) + cv_s_gauss14*(cv_s_gauss41 - tnnn*vnnn_ave(0,2)) + cv_s_gauss40*cv_s_gauss9) - N[1]*f(1,2) + N[1]*(cv_s_gauss12*(cv_s_gauss42 - cv_s_gauss44) + cv_s_gauss14*(cv_s_gauss44 - tnnn*vnnn_ave(1,2)) + cv_s_gauss43*cv_s_gauss9) - N[2]*f(2,2) + N[2]*(cv_s_gauss12*(cv_s_gauss45 - cv_s_gauss47) + cv_s_gauss14*(cv_s_gauss47 - tnnn*vnnn_ave(2,2)) + cv_s_gauss46*cv_s_gauss9) - N[3]*f(3,2) + N[3]*(cv_s_gauss12*(cv_s_gauss48 - cv_s_gauss50) + cv_s_gauss14*(cv_s_gauss50 - tnnn*vnnn_ave(3,2)) + cv_s_gauss49*cv_s_gauss9) + cv_s_gauss24*(DN(0,0)*cv_s_gauss40 + DN(1,0)*cv_s_gauss43 + DN(2,0)*cv_s_gauss46 + DN(3,0)*cv_s_gauss49) + cv_s_gauss25*(DN(0,1)*cv_s_gauss40 + DN(1,1)*cv_s_gauss43 + DN(2,1)*cv_s_gauss46 + DN(3,1)*cv_s_gauss49) + cv_s_gauss26*(DN(0,2)*cv_s_gauss40 + DN(1,2)*cv_s_gauss43 + DN(2,2)*cv_s_gauss46 + DN(3,2)*cv_s_gauss49)));


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
    const double& dtn = data.dtn;
    const double& dtnn = data.dtnn;
    const double& dtnnn = data.dtnnn;

    const double& t = data.t;
    const double& tn = data.tn;
    const double& tnn = data.tnn;
    const double& tnnn = data.tnnn;

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

    // Get time accurate expression
    const array_1d<double,nnodes> p = (t * p_ave - tn * pn_ave) * (1.0 / dt );
    const BoundedMatrix<double,nnodes,dim> v = (t * v_ave - tn * vn_ave) * (1.0 / dt );
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
    const double cv_s_gauss0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double cv_s_gauss1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double cv_s_gauss2 =             1.0/(rho*stab_c2*sqrt(pow(cv_s_gauss0, 2) + pow(cv_s_gauss1, 2))/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double cv_s_gauss3 =             1.0/dt;
const double cv_s_gauss4 =             p_ave[0]*t - pn_ave[0]*tn;
const double cv_s_gauss5 =             p_ave[1]*t - pn_ave[1]*tn;
const double cv_s_gauss6 =             p_ave[2]*t - pn_ave[2]*tn;
const double cv_s_gauss7 =             bdf0*cv_s_gauss3;
const double cv_s_gauss8 =             tn*vn_ave(0,0);
const double cv_s_gauss9 =             -cv_s_gauss8 + t*v_ave(0,0);
const double cv_s_gauss10 =             bdf1/dtn;
const double cv_s_gauss11 =             tnn*vnn_ave(0,0);
const double cv_s_gauss12 =             bdf2/dtnn;
const double cv_s_gauss13 =             tn*vn_ave(1,0);
const double cv_s_gauss14 =             -cv_s_gauss13 + t*v_ave(1,0);
const double cv_s_gauss15 =             tnn*vnn_ave(1,0);
const double cv_s_gauss16 =             tn*vn_ave(2,0);
const double cv_s_gauss17 =             -cv_s_gauss16 + t*v_ave(2,0);
const double cv_s_gauss18 =             tnn*vnn_ave(2,0);
const double cv_s_gauss19 =             cv_s_gauss0*cv_s_gauss3;
const double cv_s_gauss20 =             cv_s_gauss1*cv_s_gauss3;
const double cv_s_gauss21 =             tn*vn_ave(0,1);
const double cv_s_gauss22 =             -cv_s_gauss21 + t*v_ave(0,1);
const double cv_s_gauss23 =             tnn*vnn_ave(0,1);
const double cv_s_gauss24 =             tn*vn_ave(1,1);
const double cv_s_gauss25 =             -cv_s_gauss24 + t*v_ave(1,1);
const double cv_s_gauss26 =             tnn*vnn_ave(1,1);
const double cv_s_gauss27 =             tn*vn_ave(2,1);
const double cv_s_gauss28 =             -cv_s_gauss27 + t*v_ave(2,1);
const double cv_s_gauss29 =             tnn*vnn_ave(2,1);
            v_s_gauss[0]=-cv_s_gauss2*(cv_s_gauss3*(DN(0,0)*cv_s_gauss4 + DN(1,0)*cv_s_gauss5 + DN(2,0)*cv_s_gauss6) + rho*(-N[0]*f(0,0) + N[0]*(cv_s_gauss10*(-cv_s_gauss11 + cv_s_gauss8) + cv_s_gauss12*(cv_s_gauss11 - tnnn*vnnn_ave(0,0)) + cv_s_gauss7*cv_s_gauss9) - N[1]*f(1,0) + N[1]*(cv_s_gauss10*(cv_s_gauss13 - cv_s_gauss15) + cv_s_gauss12*(cv_s_gauss15 - tnnn*vnnn_ave(1,0)) + cv_s_gauss14*cv_s_gauss7) - N[2]*f(2,0) + N[2]*(cv_s_gauss10*(cv_s_gauss16 - cv_s_gauss18) + cv_s_gauss12*(cv_s_gauss18 - tnnn*vnnn_ave(2,0)) + cv_s_gauss17*cv_s_gauss7) + cv_s_gauss19*(DN(0,0)*cv_s_gauss9 + DN(1,0)*cv_s_gauss14 + DN(2,0)*cv_s_gauss17) + cv_s_gauss20*(DN(0,1)*cv_s_gauss9 + DN(1,1)*cv_s_gauss14 + DN(2,1)*cv_s_gauss17)));
            v_s_gauss[1]=-cv_s_gauss2*(cv_s_gauss3*(DN(0,1)*cv_s_gauss4 + DN(1,1)*cv_s_gauss5 + DN(2,1)*cv_s_gauss6) + rho*(-N[0]*f(0,1) + N[0]*(cv_s_gauss10*(cv_s_gauss21 - cv_s_gauss23) + cv_s_gauss12*(cv_s_gauss23 - tnnn*vnnn_ave(0,1)) + cv_s_gauss22*cv_s_gauss7) - N[1]*f(1,1) + N[1]*(cv_s_gauss10*(cv_s_gauss24 - cv_s_gauss26) + cv_s_gauss12*(cv_s_gauss26 - tnnn*vnnn_ave(1,1)) + cv_s_gauss25*cv_s_gauss7) - N[2]*f(2,1) + N[2]*(cv_s_gauss10*(cv_s_gauss27 - cv_s_gauss29) + cv_s_gauss12*(cv_s_gauss29 - tnnn*vnnn_ave(2,1)) + cv_s_gauss28*cv_s_gauss7) + cv_s_gauss19*(DN(0,0)*cv_s_gauss22 + DN(1,0)*cv_s_gauss25 + DN(2,0)*cv_s_gauss28) + cv_s_gauss20*(DN(0,1)*cv_s_gauss22 + DN(1,1)*cv_s_gauss25 + DN(2,1)*cv_s_gauss28)));


    const double v_gauss_norm = norm_2(v_gauss);
    const double v_s_gauss_norm = norm_2(v_s_gauss);

    return v_s_gauss_norm/v_gauss_norm;
}

}
