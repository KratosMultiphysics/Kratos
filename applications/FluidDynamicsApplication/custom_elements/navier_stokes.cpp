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

#include "custom_elements/navier_stokes.h"

namespace Kratos {

template<>
void NavierStokes<3>::EquationIdVector(
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
        rResult[i*(dim+1)  ]  =  this->GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
        rResult[i*(dim+1)+1]  =  this->GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
        rResult[i*(dim+1)+2]  =  this->GetGeometry()[i].GetDof(VELOCITY_Z).EquationId();
        rResult[i*(dim+1)+3]  =  this->GetGeometry()[i].GetDof(PRESSURE).EquationId();
    }

    KRATOS_CATCH("")
}


template<>
void NavierStokes<2>::EquationIdVector(
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
        rResult[i*(dim+1)  ]  =  this->GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
        rResult[i*(dim+1)+1]  =  this->GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
        rResult[i*(dim+1)+2]  =  this->GetGeometry()[i].GetDof(PRESSURE).EquationId();
    }

    KRATOS_CATCH("")
}


template<>
void NavierStokes<3>::GetDofList(
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
        ElementalDofList[i*(dim+1)  ]  =  this->GetGeometry()[i].pGetDof(VELOCITY_X);
        ElementalDofList[i*(dim+1)+1]  =  this->GetGeometry()[i].pGetDof(VELOCITY_Y);
        ElementalDofList[i*(dim+1)+2]  =  this->GetGeometry()[i].pGetDof(VELOCITY_Z);
        ElementalDofList[i*(dim+1)+3]  =  this->GetGeometry()[i].pGetDof(PRESSURE);
    }

    KRATOS_CATCH("");
}


template<>
void NavierStokes<2>::GetDofList(
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
        ElementalDofList[i*(dim+1)  ]  =  this->GetGeometry()[i].pGetDof(VELOCITY_X);
        ElementalDofList[i*(dim+1)+1]  =  this->GetGeometry()[i].pGetDof(VELOCITY_Y);
        ElementalDofList[i*(dim+1)+2]  =  this->GetGeometry()[i].pGetDof(PRESSURE);
    }

    KRATOS_CATCH("");
}


template<>
void NavierStokes<3>::ComputeGaussPointLHSContribution(
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

    const BoundedMatrix<double,nnodes,dim>& v = data.v;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
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
const double clhs11 =             pow(N[0], 2);
const double clhs12 =             bdf0*rho;
const double clhs13 =             N[0]*rho;
const double clhs14 =             DN(0,0)*clhs6 + DN(0,1)*clhs7 + DN(0,2)*clhs8;
const double clhs15 =             N[0]*bdf0 + clhs14;
const double clhs16 =             pow(rho, 2);
const double clhs17 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double clhs18 =             1.0/(clhs9/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs19 =             1.0*N[0]*clhs16*clhs17*clhs18;
const double clhs20 =             1.0*clhs14*clhs16*clhs18;
const double clhs21 =             clhs11*clhs12 + clhs13*clhs14 + clhs15*clhs19 + clhs15*clhs20;
const double clhs22 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs1;
const double clhs23 =             C(1,3)*DN(0,1);
const double clhs24 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs23;
const double clhs25 =             C(3,5)*DN(0,0);
const double clhs26 =             C(4,5)*DN(0,2);
const double clhs27 =             C(1,5)*DN(0,1) + clhs25 + clhs26;
const double clhs28 =             DN(0,0)*clhs10;
const double clhs29 =             DN(0,1)*clhs28;
const double clhs30 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs3;
const double clhs31 =             C(3,4)*DN(0,1);
const double clhs32 =             C(2,3)*DN(0,2) + clhs25 + clhs31;
const double clhs33 =             C(2,5)*DN(0,2);
const double clhs34 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs33;
const double clhs35 =             DN(0,2)*clhs28;
const double clhs36 =             pow(c, -2);
const double clhs37 =             1.0/rho;
const double clhs38 =             N[0]*bdf0*clhs36*clhs37;
const double clhs39 =             1.0*clhs17*clhs18;
const double clhs40 =             1.0*clhs18*rho;
const double clhs41 =             clhs14*clhs40;
const double clhs42 =             -N[0] + clhs10*clhs38 + clhs13*clhs39 + clhs41;
const double clhs43 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs44 =             C(0,3)*DN(1,0);
const double clhs45 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs44;
const double clhs46 =             C(0,5)*DN(1,0);
const double clhs47 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs46;
const double clhs48 =             DN(1,0)*clhs28;
const double clhs49 =             N[0]*bdf0*rho;
const double clhs50 =             N[1]*clhs49;
const double clhs51 =             DN(1,0)*clhs6 + DN(1,1)*clhs7 + DN(1,2)*clhs8;
const double clhs52 =             N[1]*bdf0 + clhs51;
const double clhs53 =             clhs13*clhs51 + clhs19*clhs52 + clhs20*clhs52 + clhs50;
const double clhs54 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs44;
const double clhs55 =             C(1,3)*DN(1,1);
const double clhs56 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs55;
const double clhs57 =             C(3,5)*DN(1,0);
const double clhs58 =             C(4,5)*DN(1,2);
const double clhs59 =             C(1,5)*DN(1,1) + clhs57 + clhs58;
const double clhs60 =             DN(1,1)*clhs28;
const double clhs61 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs46;
const double clhs62 =             C(3,4)*DN(1,1);
const double clhs63 =             C(2,3)*DN(1,2) + clhs57 + clhs62;
const double clhs64 =             C(2,5)*DN(1,2);
const double clhs65 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs64;
const double clhs66 =             DN(1,2)*clhs28;
const double clhs67 =             DN(0,0)*N[1];
const double clhs68 =             bdf0*clhs10*clhs36*clhs37;
const double clhs69 =             DN(1,0)*N[0];
const double clhs70 =             1.0*clhs17*clhs18*rho;
const double clhs71 =             1.0*DN(1,0)*clhs18*rho;
const double clhs72 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs73 =             C(0,3)*DN(2,0);
const double clhs74 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs73;
const double clhs75 =             C(0,5)*DN(2,0);
const double clhs76 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs75;
const double clhs77 =             DN(2,0)*clhs28;
const double clhs78 =             N[2]*clhs49;
const double clhs79 =             DN(2,0)*clhs6 + DN(2,1)*clhs7 + DN(2,2)*clhs8;
const double clhs80 =             N[2]*bdf0;
const double clhs81 =             clhs79 + clhs80;
const double clhs82 =             clhs13*clhs79 + clhs19*clhs81 + clhs20*clhs81 + clhs78;
const double clhs83 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs73;
const double clhs84 =             C(1,3)*DN(2,1);
const double clhs85 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs84;
const double clhs86 =             C(3,5)*DN(2,0);
const double clhs87 =             C(4,5)*DN(2,2);
const double clhs88 =             C(1,5)*DN(2,1) + clhs86 + clhs87;
const double clhs89 =             DN(2,1)*clhs28;
const double clhs90 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs75;
const double clhs91 =             C(3,4)*DN(2,1);
const double clhs92 =             C(2,3)*DN(2,2) + clhs86 + clhs91;
const double clhs93 =             C(2,5)*DN(2,2);
const double clhs94 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs93;
const double clhs95 =             DN(2,2)*clhs28;
const double clhs96 =             DN(0,0)*N[2];
const double clhs97 =             DN(2,0)*N[0];
const double clhs98 =             1.0*DN(2,0)*clhs18*rho;
const double clhs99 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs100 =             C(0,3)*DN(3,0);
const double clhs101 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs100;
const double clhs102 =             C(0,5)*DN(3,0);
const double clhs103 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs102;
const double clhs104 =             DN(3,0)*clhs28;
const double clhs105 =             N[3]*clhs49;
const double clhs106 =             DN(3,0)*clhs6 + DN(3,1)*clhs7 + DN(3,2)*clhs8;
const double clhs107 =             N[3]*bdf0;
const double clhs108 =             clhs106 + clhs107;
const double clhs109 =             clhs105 + clhs106*clhs13 + clhs108*clhs19 + clhs108*clhs20;
const double clhs110 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs100;
const double clhs111 =             C(1,3)*DN(3,1);
const double clhs112 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs111;
const double clhs113 =             C(3,5)*DN(3,0);
const double clhs114 =             C(4,5)*DN(3,2);
const double clhs115 =             C(1,5)*DN(3,1) + clhs113 + clhs114;
const double clhs116 =             DN(3,1)*clhs28;
const double clhs117 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs102;
const double clhs118 =             C(3,4)*DN(3,1);
const double clhs119 =             C(2,3)*DN(3,2) + clhs113 + clhs118;
const double clhs120 =             C(2,5)*DN(3,2);
const double clhs121 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs120;
const double clhs122 =             DN(3,2)*clhs28;
const double clhs123 =             DN(0,0)*N[3];
const double clhs124 =             DN(3,0)*N[0];
const double clhs125 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs23;
const double clhs126 =             C(0,4)*DN(0,0) + clhs26 + clhs31;
const double clhs127 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs128 =             C(1,4)*DN(0,1);
const double clhs129 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs128;
const double clhs130 =             pow(DN(0,1), 2);
const double clhs131 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs128;
const double clhs132 =             C(2,4)*DN(0,2);
const double clhs133 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs132;
const double clhs134 =             DN(0,1)*clhs10;
const double clhs135 =             DN(0,2)*clhs134;
const double clhs136 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs55;
const double clhs137 =             C(0,4)*DN(1,0) + clhs58 + clhs62;
const double clhs138 =             DN(1,0)*clhs134;
const double clhs139 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs140 =             C(1,4)*DN(1,1);
const double clhs141 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs140;
const double clhs142 =             DN(1,1)*clhs134;
const double clhs143 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs140;
const double clhs144 =             C(2,4)*DN(1,2);
const double clhs145 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs144;
const double clhs146 =             DN(1,2)*clhs134;
const double clhs147 =             DN(0,1)*N[1];
const double clhs148 =             DN(1,1)*N[0];
const double clhs149 =             1.0*DN(1,1)*clhs18*rho;
const double clhs150 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs84;
const double clhs151 =             C(0,4)*DN(2,0) + clhs87 + clhs91;
const double clhs152 =             DN(2,0)*clhs134;
const double clhs153 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs154 =             C(1,4)*DN(2,1);
const double clhs155 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs154;
const double clhs156 =             DN(2,1)*clhs134;
const double clhs157 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs154;
const double clhs158 =             C(2,4)*DN(2,2);
const double clhs159 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs158;
const double clhs160 =             DN(2,2)*clhs134;
const double clhs161 =             DN(0,1)*N[2];
const double clhs162 =             DN(2,1)*N[0];
const double clhs163 =             1.0*DN(2,1)*clhs18*rho;
const double clhs164 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs111;
const double clhs165 =             C(0,4)*DN(3,0) + clhs114 + clhs118;
const double clhs166 =             DN(3,0)*clhs134;
const double clhs167 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs168 =             C(1,4)*DN(3,1);
const double clhs169 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs168;
const double clhs170 =             DN(3,1)*clhs134;
const double clhs171 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs168;
const double clhs172 =             C(2,4)*DN(3,2);
const double clhs173 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs172;
const double clhs174 =             DN(3,2)*clhs134;
const double clhs175 =             DN(0,1)*N[3];
const double clhs176 =             DN(3,1)*N[0];
const double clhs177 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs33;
const double clhs178 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs132;
const double clhs179 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs180 =             pow(DN(0,2), 2);
const double clhs181 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs64;
const double clhs182 =             DN(0,2)*clhs10;
const double clhs183 =             DN(1,0)*clhs182;
const double clhs184 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs144;
const double clhs185 =             DN(1,1)*clhs182;
const double clhs186 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs187 =             DN(1,2)*clhs182;
const double clhs188 =             DN(0,2)*N[1];
const double clhs189 =             DN(1,2)*N[0];
const double clhs190 =             1.0*DN(1,2)*clhs18*rho;
const double clhs191 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs93;
const double clhs192 =             DN(2,0)*clhs182;
const double clhs193 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs158;
const double clhs194 =             DN(2,1)*clhs182;
const double clhs195 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs196 =             DN(2,2)*clhs182;
const double clhs197 =             DN(0,2)*N[2];
const double clhs198 =             DN(2,2)*N[0];
const double clhs199 =             1.0*DN(2,2)*clhs18*rho;
const double clhs200 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs120;
const double clhs201 =             DN(3,0)*clhs182;
const double clhs202 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs172;
const double clhs203 =             DN(3,1)*clhs182;
const double clhs204 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs205 =             DN(3,2)*clhs182;
const double clhs206 =             DN(0,2)*N[3];
const double clhs207 =             DN(3,2)*N[0];
const double clhs208 =             clhs15*clhs40;
const double clhs209 =             N[0] + clhs208;
const double clhs210 =             bdf0*clhs36*clhs37;
const double clhs211 =             1.0*clhs18;
const double clhs212 =             1.0*DN(0,0)*clhs18*rho;
const double clhs213 =             1.0*DN(0,1)*clhs18*rho;
const double clhs214 =             1.0*DN(0,2)*clhs18*rho;
const double clhs215 =             1.0*DN(0,0)*clhs18;
const double clhs216 =             1.0*DN(0,1)*clhs18;
const double clhs217 =             1.0*DN(0,2)*clhs18;
const double clhs218 =             DN(1,0)*clhs215 + DN(1,1)*clhs216 + DN(1,2)*clhs217 + N[1]*clhs38;
const double clhs219 =             DN(2,0)*clhs215 + DN(2,1)*clhs216 + DN(2,2)*clhs217 + N[2]*clhs38;
const double clhs220 =             DN(3,0)*clhs215 + DN(3,1)*clhs216 + DN(3,2)*clhs217 + N[3]*clhs38;
const double clhs221 =             N[1]*rho;
const double clhs222 =             1.0*N[1]*clhs16*clhs17*clhs18;
const double clhs223 =             1.0*clhs16*clhs18*clhs51;
const double clhs224 =             clhs14*clhs221 + clhs15*clhs222 + clhs15*clhs223 + clhs50;
const double clhs225 =             pow(DN(1,0), 2);
const double clhs226 =             pow(N[1], 2);
const double clhs227 =             clhs12*clhs226 + clhs221*clhs51 + clhs222*clhs52 + clhs223*clhs52;
const double clhs228 =             DN(1,0)*clhs10;
const double clhs229 =             DN(1,1)*clhs228;
const double clhs230 =             DN(1,2)*clhs228;
const double clhs231 =             N[1]*bdf0*clhs36*clhs37;
const double clhs232 =             clhs40*clhs51;
const double clhs233 =             -N[1] + clhs10*clhs231 + clhs221*clhs39 + clhs232;
const double clhs234 =             DN(2,0)*clhs228;
const double clhs235 =             N[1]*bdf0*rho;
const double clhs236 =             N[2]*clhs235;
const double clhs237 =             clhs221*clhs79 + clhs222*clhs81 + clhs223*clhs81 + clhs236;
const double clhs238 =             DN(2,1)*clhs228;
const double clhs239 =             DN(2,2)*clhs228;
const double clhs240 =             DN(1,0)*N[2];
const double clhs241 =             DN(2,0)*N[1];
const double clhs242 =             DN(3,0)*clhs228;
const double clhs243 =             N[3]*clhs235;
const double clhs244 =             clhs106*clhs221 + clhs108*clhs222 + clhs108*clhs223 + clhs243;
const double clhs245 =             DN(3,1)*clhs228;
const double clhs246 =             DN(3,2)*clhs228;
const double clhs247 =             DN(1,0)*N[3];
const double clhs248 =             DN(3,0)*N[1];
const double clhs249 =             pow(DN(1,1), 2);
const double clhs250 =             DN(1,1)*clhs10;
const double clhs251 =             DN(1,2)*clhs250;
const double clhs252 =             DN(2,0)*clhs250;
const double clhs253 =             DN(2,1)*clhs250;
const double clhs254 =             DN(2,2)*clhs250;
const double clhs255 =             DN(1,1)*N[2];
const double clhs256 =             DN(2,1)*N[1];
const double clhs257 =             DN(3,0)*clhs250;
const double clhs258 =             DN(3,1)*clhs250;
const double clhs259 =             DN(3,2)*clhs250;
const double clhs260 =             DN(1,1)*N[3];
const double clhs261 =             DN(3,1)*N[1];
const double clhs262 =             pow(DN(1,2), 2);
const double clhs263 =             DN(1,2)*clhs10;
const double clhs264 =             DN(2,0)*clhs263;
const double clhs265 =             DN(2,1)*clhs263;
const double clhs266 =             DN(2,2)*clhs263;
const double clhs267 =             DN(1,2)*N[2];
const double clhs268 =             DN(2,2)*N[1];
const double clhs269 =             DN(3,0)*clhs263;
const double clhs270 =             DN(3,1)*clhs263;
const double clhs271 =             DN(3,2)*clhs263;
const double clhs272 =             DN(1,2)*N[3];
const double clhs273 =             DN(3,2)*N[1];
const double clhs274 =             clhs40*clhs52;
const double clhs275 =             N[1] + clhs274;
const double clhs276 =             1.0*DN(1,0)*clhs18;
const double clhs277 =             1.0*DN(1,1)*clhs18;
const double clhs278 =             1.0*DN(1,2)*clhs18;
const double clhs279 =             DN(2,0)*clhs276 + DN(2,1)*clhs277 + DN(2,2)*clhs278 + N[2]*clhs231;
const double clhs280 =             DN(3,0)*clhs276 + DN(3,1)*clhs277 + DN(3,2)*clhs278 + N[3]*clhs231;
const double clhs281 =             N[2]*rho;
const double clhs282 =             1.0*N[2]*clhs16*clhs17*clhs18;
const double clhs283 =             1.0*clhs16*clhs18*clhs79;
const double clhs284 =             clhs14*clhs281 + clhs15*clhs282 + clhs15*clhs283 + clhs78;
const double clhs285 =             clhs236 + clhs281*clhs51 + clhs282*clhs52 + clhs283*clhs52;
const double clhs286 =             pow(DN(2,0), 2);
const double clhs287 =             pow(N[2], 2);
const double clhs288 =             clhs12*clhs287 + clhs281*clhs79 + clhs282*clhs81 + clhs283*clhs81;
const double clhs289 =             DN(2,0)*clhs10;
const double clhs290 =             DN(2,1)*clhs289;
const double clhs291 =             DN(2,2)*clhs289;
const double clhs292 =             clhs10*clhs36*clhs37;
const double clhs293 =             clhs40*clhs79;
const double clhs294 =             -N[2] + clhs281*clhs39 + clhs292*clhs80 + clhs293;
const double clhs295 =             DN(3,0)*clhs289;
const double clhs296 =             N[2]*N[3]*bdf0;
const double clhs297 =             clhs296*rho;
const double clhs298 =             clhs106*clhs281 + clhs108*clhs282 + clhs108*clhs283 + clhs297;
const double clhs299 =             DN(3,1)*clhs289;
const double clhs300 =             DN(3,2)*clhs289;
const double clhs301 =             DN(2,0)*N[3];
const double clhs302 =             DN(3,0)*N[2];
const double clhs303 =             pow(DN(2,1), 2);
const double clhs304 =             DN(2,1)*clhs10;
const double clhs305 =             DN(2,2)*clhs304;
const double clhs306 =             DN(3,0)*clhs304;
const double clhs307 =             DN(3,1)*clhs304;
const double clhs308 =             DN(3,2)*clhs304;
const double clhs309 =             DN(2,1)*N[3];
const double clhs310 =             DN(3,1)*N[2];
const double clhs311 =             pow(DN(2,2), 2);
const double clhs312 =             DN(2,2)*clhs10;
const double clhs313 =             DN(3,0)*clhs312;
const double clhs314 =             DN(3,1)*clhs312;
const double clhs315 =             DN(3,2)*clhs312;
const double clhs316 =             DN(2,2)*N[3];
const double clhs317 =             DN(3,2)*N[2];
const double clhs318 =             clhs40*clhs81;
const double clhs319 =             N[2] + clhs318;
const double clhs320 =             1.0*DN(2,0)*DN(3,0)*clhs18 + 1.0*DN(2,1)*DN(3,1)*clhs18 + 1.0*DN(2,2)*DN(3,2)*clhs18 + clhs296*clhs36*clhs37;
const double clhs321 =             N[3]*rho;
const double clhs322 =             1.0*N[3]*clhs16*clhs17*clhs18;
const double clhs323 =             1.0*clhs106*clhs16*clhs18;
const double clhs324 =             clhs105 + clhs14*clhs321 + clhs15*clhs322 + clhs15*clhs323;
const double clhs325 =             clhs243 + clhs321*clhs51 + clhs322*clhs52 + clhs323*clhs52;
const double clhs326 =             clhs297 + clhs321*clhs79 + clhs322*clhs81 + clhs323*clhs81;
const double clhs327 =             pow(DN(3,0), 2);
const double clhs328 =             pow(N[3], 2);
const double clhs329 =             clhs106*clhs321 + clhs108*clhs322 + clhs108*clhs323 + clhs12*clhs328;
const double clhs330 =             DN(3,0)*clhs10;
const double clhs331 =             DN(3,1)*clhs330;
const double clhs332 =             DN(3,2)*clhs330;
const double clhs333 =             -N[3] + clhs106*clhs40 + clhs107*clhs292 + clhs321*clhs39;
const double clhs334 =             pow(DN(3,1), 2);
const double clhs335 =             DN(3,1)*DN(3,2)*clhs10;
const double clhs336 =             pow(DN(3,2), 2);
const double clhs337 =             N[3] + clhs108*clhs40;
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs10*clhs5 + clhs21;
            lhs(0,1)=DN(0,0)*clhs22 + DN(0,1)*clhs24 + DN(0,2)*clhs27 + clhs29;
            lhs(0,2)=DN(0,0)*clhs30 + DN(0,1)*clhs32 + DN(0,2)*clhs34 + clhs35;
            lhs(0,3)=DN(0,0)*clhs42;
            lhs(0,4)=DN(0,0)*clhs43 + DN(0,1)*clhs45 + DN(0,2)*clhs47 + clhs48 + clhs53;
            lhs(0,5)=DN(0,0)*clhs54 + DN(0,1)*clhs56 + DN(0,2)*clhs59 + clhs60;
            lhs(0,6)=DN(0,0)*clhs61 + DN(0,1)*clhs63 + DN(0,2)*clhs65 + clhs66;
            lhs(0,7)=clhs14*clhs71 + clhs67*clhs68 - clhs67 + clhs69*clhs70;
            lhs(0,8)=DN(0,0)*clhs72 + DN(0,1)*clhs74 + DN(0,2)*clhs76 + clhs77 + clhs82;
            lhs(0,9)=DN(0,0)*clhs83 + DN(0,1)*clhs85 + DN(0,2)*clhs88 + clhs89;
            lhs(0,10)=DN(0,0)*clhs90 + DN(0,1)*clhs92 + DN(0,2)*clhs94 + clhs95;
            lhs(0,11)=clhs14*clhs98 + clhs68*clhs96 + clhs70*clhs97 - clhs96;
            lhs(0,12)=DN(0,0)*clhs99 + DN(0,1)*clhs101 + DN(0,2)*clhs103 + clhs104 + clhs109;
            lhs(0,13)=DN(0,0)*clhs110 + DN(0,1)*clhs112 + DN(0,2)*clhs115 + clhs116;
            lhs(0,14)=DN(0,0)*clhs117 + DN(0,1)*clhs119 + DN(0,2)*clhs121 + clhs122;
            lhs(0,15)=DN(3,0)*clhs41 + clhs123*clhs68 - clhs123 + clhs124*clhs70;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs125 + DN(0,2)*clhs126 + clhs29;
            lhs(1,1)=DN(0,0)*clhs24 + DN(0,1)*clhs127 + DN(0,2)*clhs129 + clhs10*clhs130 + clhs21;
            lhs(1,2)=DN(0,0)*clhs32 + DN(0,1)*clhs131 + DN(0,2)*clhs133 + clhs135;
            lhs(1,3)=DN(0,1)*clhs42;
            lhs(1,4)=DN(0,0)*clhs45 + DN(0,1)*clhs136 + DN(0,2)*clhs137 + clhs138;
            lhs(1,5)=DN(0,0)*clhs56 + DN(0,1)*clhs139 + DN(0,2)*clhs141 + clhs142 + clhs53;
            lhs(1,6)=DN(0,0)*clhs63 + DN(0,1)*clhs143 + DN(0,2)*clhs145 + clhs146;
            lhs(1,7)=clhs14*clhs149 + clhs147*clhs68 - clhs147 + clhs148*clhs70;
            lhs(1,8)=DN(0,0)*clhs74 + DN(0,1)*clhs150 + DN(0,2)*clhs151 + clhs152;
            lhs(1,9)=DN(0,0)*clhs85 + DN(0,1)*clhs153 + DN(0,2)*clhs155 + clhs156 + clhs82;
            lhs(1,10)=DN(0,0)*clhs92 + DN(0,1)*clhs157 + DN(0,2)*clhs159 + clhs160;
            lhs(1,11)=clhs14*clhs163 + clhs161*clhs68 - clhs161 + clhs162*clhs70;
            lhs(1,12)=DN(0,0)*clhs101 + DN(0,1)*clhs164 + DN(0,2)*clhs165 + clhs166;
            lhs(1,13)=DN(0,0)*clhs112 + DN(0,1)*clhs167 + DN(0,2)*clhs169 + clhs109 + clhs170;
            lhs(1,14)=DN(0,0)*clhs119 + DN(0,1)*clhs171 + DN(0,2)*clhs173 + clhs174;
            lhs(1,15)=DN(3,1)*clhs41 + clhs175*clhs68 - clhs175 + clhs176*clhs70;
            lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs126 + DN(0,2)*clhs177 + clhs35;
            lhs(2,1)=DN(0,0)*clhs27 + DN(0,1)*clhs129 + DN(0,2)*clhs178 + clhs135;
            lhs(2,2)=DN(0,0)*clhs34 + DN(0,1)*clhs133 + DN(0,2)*clhs179 + clhs10*clhs180 + clhs21;
            lhs(2,3)=DN(0,2)*clhs42;
            lhs(2,4)=DN(0,0)*clhs47 + DN(0,1)*clhs137 + DN(0,2)*clhs181 + clhs183;
            lhs(2,5)=DN(0,0)*clhs59 + DN(0,1)*clhs141 + DN(0,2)*clhs184 + clhs185;
            lhs(2,6)=DN(0,0)*clhs65 + DN(0,1)*clhs145 + DN(0,2)*clhs186 + clhs187 + clhs53;
            lhs(2,7)=clhs14*clhs190 + clhs188*clhs68 - clhs188 + clhs189*clhs70;
            lhs(2,8)=DN(0,0)*clhs76 + DN(0,1)*clhs151 + DN(0,2)*clhs191 + clhs192;
            lhs(2,9)=DN(0,0)*clhs88 + DN(0,1)*clhs155 + DN(0,2)*clhs193 + clhs194;
            lhs(2,10)=DN(0,0)*clhs94 + DN(0,1)*clhs159 + DN(0,2)*clhs195 + clhs196 + clhs82;
            lhs(2,11)=clhs14*clhs199 + clhs197*clhs68 - clhs197 + clhs198*clhs70;
            lhs(2,12)=DN(0,0)*clhs103 + DN(0,1)*clhs165 + DN(0,2)*clhs200 + clhs201;
            lhs(2,13)=DN(0,0)*clhs115 + DN(0,1)*clhs169 + DN(0,2)*clhs202 + clhs203;
            lhs(2,14)=DN(0,0)*clhs121 + DN(0,1)*clhs173 + DN(0,2)*clhs204 + clhs109 + clhs205;
            lhs(2,15)=DN(3,2)*clhs41 + clhs206*clhs68 - clhs206 + clhs207*clhs70;
            lhs(3,0)=DN(0,0)*clhs209;
            lhs(3,1)=DN(0,1)*clhs209;
            lhs(3,2)=DN(0,2)*clhs209;
            lhs(3,3)=clhs11*clhs210 + clhs130*clhs211 + clhs180*clhs211 + clhs211*clhs5;
            lhs(3,4)=clhs212*clhs52 + clhs69;
            lhs(3,5)=clhs148 + clhs213*clhs52;
            lhs(3,6)=clhs189 + clhs214*clhs52;
            lhs(3,7)=clhs218;
            lhs(3,8)=clhs212*clhs81 + clhs97;
            lhs(3,9)=clhs162 + clhs213*clhs81;
            lhs(3,10)=clhs198 + clhs214*clhs81;
            lhs(3,11)=clhs219;
            lhs(3,12)=clhs108*clhs212 + clhs124;
            lhs(3,13)=clhs108*clhs213 + clhs176;
            lhs(3,14)=clhs108*clhs214 + clhs207;
            lhs(3,15)=clhs220;
            lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs224 + clhs48;
            lhs(4,1)=DN(1,0)*clhs22 + DN(1,1)*clhs24 + DN(1,2)*clhs27 + clhs138;
            lhs(4,2)=DN(1,0)*clhs30 + DN(1,1)*clhs32 + DN(1,2)*clhs34 + clhs183;
            lhs(4,3)=clhs212*clhs51 + clhs67*clhs70 + clhs68*clhs69 - clhs69;
            lhs(4,4)=DN(1,0)*clhs43 + DN(1,1)*clhs45 + DN(1,2)*clhs47 + clhs10*clhs225 + clhs227;
            lhs(4,5)=DN(1,0)*clhs54 + DN(1,1)*clhs56 + DN(1,2)*clhs59 + clhs229;
            lhs(4,6)=DN(1,0)*clhs61 + DN(1,1)*clhs63 + DN(1,2)*clhs65 + clhs230;
            lhs(4,7)=DN(1,0)*clhs233;
            lhs(4,8)=DN(1,0)*clhs72 + DN(1,1)*clhs74 + DN(1,2)*clhs76 + clhs234 + clhs237;
            lhs(4,9)=DN(1,0)*clhs83 + DN(1,1)*clhs85 + DN(1,2)*clhs88 + clhs238;
            lhs(4,10)=DN(1,0)*clhs90 + DN(1,1)*clhs92 + DN(1,2)*clhs94 + clhs239;
            lhs(4,11)=clhs240*clhs68 - clhs240 + clhs241*clhs70 + clhs51*clhs98;
            lhs(4,12)=DN(1,0)*clhs99 + DN(1,1)*clhs101 + DN(1,2)*clhs103 + clhs242 + clhs244;
            lhs(4,13)=DN(1,0)*clhs110 + DN(1,1)*clhs112 + DN(1,2)*clhs115 + clhs245;
            lhs(4,14)=DN(1,0)*clhs117 + DN(1,1)*clhs119 + DN(1,2)*clhs121 + clhs246;
            lhs(4,15)=DN(3,0)*clhs232 + clhs247*clhs68 - clhs247 + clhs248*clhs70;
            lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs125 + DN(1,2)*clhs126 + clhs60;
            lhs(5,1)=DN(1,0)*clhs24 + DN(1,1)*clhs127 + DN(1,2)*clhs129 + clhs142 + clhs224;
            lhs(5,2)=DN(1,0)*clhs32 + DN(1,1)*clhs131 + DN(1,2)*clhs133 + clhs185;
            lhs(5,3)=clhs147*clhs70 + clhs148*clhs68 - clhs148 + clhs213*clhs51;
            lhs(5,4)=DN(1,0)*clhs45 + DN(1,1)*clhs136 + DN(1,2)*clhs137 + clhs229;
            lhs(5,5)=DN(1,0)*clhs56 + DN(1,1)*clhs139 + DN(1,2)*clhs141 + clhs10*clhs249 + clhs227;
            lhs(5,6)=DN(1,0)*clhs63 + DN(1,1)*clhs143 + DN(1,2)*clhs145 + clhs251;
            lhs(5,7)=DN(1,1)*clhs233;
            lhs(5,8)=DN(1,0)*clhs74 + DN(1,1)*clhs150 + DN(1,2)*clhs151 + clhs252;
            lhs(5,9)=DN(1,0)*clhs85 + DN(1,1)*clhs153 + DN(1,2)*clhs155 + clhs237 + clhs253;
            lhs(5,10)=DN(1,0)*clhs92 + DN(1,1)*clhs157 + DN(1,2)*clhs159 + clhs254;
            lhs(5,11)=clhs163*clhs51 + clhs255*clhs68 - clhs255 + clhs256*clhs70;
            lhs(5,12)=DN(1,0)*clhs101 + DN(1,1)*clhs164 + DN(1,2)*clhs165 + clhs257;
            lhs(5,13)=DN(1,0)*clhs112 + DN(1,1)*clhs167 + DN(1,2)*clhs169 + clhs244 + clhs258;
            lhs(5,14)=DN(1,0)*clhs119 + DN(1,1)*clhs171 + DN(1,2)*clhs173 + clhs259;
            lhs(5,15)=DN(3,1)*clhs232 + clhs260*clhs68 - clhs260 + clhs261*clhs70;
            lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs126 + DN(1,2)*clhs177 + clhs66;
            lhs(6,1)=DN(1,0)*clhs27 + DN(1,1)*clhs129 + DN(1,2)*clhs178 + clhs146;
            lhs(6,2)=DN(1,0)*clhs34 + DN(1,1)*clhs133 + DN(1,2)*clhs179 + clhs187 + clhs224;
            lhs(6,3)=clhs188*clhs70 + clhs189*clhs68 - clhs189 + clhs214*clhs51;
            lhs(6,4)=DN(1,0)*clhs47 + DN(1,1)*clhs137 + DN(1,2)*clhs181 + clhs230;
            lhs(6,5)=DN(1,0)*clhs59 + DN(1,1)*clhs141 + DN(1,2)*clhs184 + clhs251;
            lhs(6,6)=DN(1,0)*clhs65 + DN(1,1)*clhs145 + DN(1,2)*clhs186 + clhs10*clhs262 + clhs227;
            lhs(6,7)=DN(1,2)*clhs233;
            lhs(6,8)=DN(1,0)*clhs76 + DN(1,1)*clhs151 + DN(1,2)*clhs191 + clhs264;
            lhs(6,9)=DN(1,0)*clhs88 + DN(1,1)*clhs155 + DN(1,2)*clhs193 + clhs265;
            lhs(6,10)=DN(1,0)*clhs94 + DN(1,1)*clhs159 + DN(1,2)*clhs195 + clhs237 + clhs266;
            lhs(6,11)=clhs199*clhs51 + clhs267*clhs68 - clhs267 + clhs268*clhs70;
            lhs(6,12)=DN(1,0)*clhs103 + DN(1,1)*clhs165 + DN(1,2)*clhs200 + clhs269;
            lhs(6,13)=DN(1,0)*clhs115 + DN(1,1)*clhs169 + DN(1,2)*clhs202 + clhs270;
            lhs(6,14)=DN(1,0)*clhs121 + DN(1,1)*clhs173 + DN(1,2)*clhs204 + clhs244 + clhs271;
            lhs(6,15)=DN(3,2)*clhs232 + clhs272*clhs68 - clhs272 + clhs273*clhs70;
            lhs(7,0)=clhs15*clhs71 + clhs67;
            lhs(7,1)=clhs147 + clhs149*clhs15;
            lhs(7,2)=clhs15*clhs190 + clhs188;
            lhs(7,3)=clhs218;
            lhs(7,4)=DN(1,0)*clhs275;
            lhs(7,5)=DN(1,1)*clhs275;
            lhs(7,6)=DN(1,2)*clhs275;
            lhs(7,7)=clhs210*clhs226 + clhs211*clhs225 + clhs211*clhs249 + clhs211*clhs262;
            lhs(7,8)=clhs241 + clhs71*clhs81;
            lhs(7,9)=clhs149*clhs81 + clhs256;
            lhs(7,10)=clhs190*clhs81 + clhs268;
            lhs(7,11)=clhs279;
            lhs(7,12)=clhs108*clhs71 + clhs248;
            lhs(7,13)=clhs108*clhs149 + clhs261;
            lhs(7,14)=clhs108*clhs190 + clhs273;
            lhs(7,15)=clhs280;
            lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs284 + clhs77;
            lhs(8,1)=DN(2,0)*clhs22 + DN(2,1)*clhs24 + DN(2,2)*clhs27 + clhs152;
            lhs(8,2)=DN(2,0)*clhs30 + DN(2,1)*clhs32 + DN(2,2)*clhs34 + clhs192;
            lhs(8,3)=clhs212*clhs79 + clhs68*clhs97 + clhs70*clhs96 - clhs97;
            lhs(8,4)=DN(2,0)*clhs43 + DN(2,1)*clhs45 + DN(2,2)*clhs47 + clhs234 + clhs285;
            lhs(8,5)=DN(2,0)*clhs54 + DN(2,1)*clhs56 + DN(2,2)*clhs59 + clhs252;
            lhs(8,6)=DN(2,0)*clhs61 + DN(2,1)*clhs63 + DN(2,2)*clhs65 + clhs264;
            lhs(8,7)=clhs240*clhs70 + clhs241*clhs68 - clhs241 + clhs71*clhs79;
            lhs(8,8)=DN(2,0)*clhs72 + DN(2,1)*clhs74 + DN(2,2)*clhs76 + clhs10*clhs286 + clhs288;
            lhs(8,9)=DN(2,0)*clhs83 + DN(2,1)*clhs85 + DN(2,2)*clhs88 + clhs290;
            lhs(8,10)=DN(2,0)*clhs90 + DN(2,1)*clhs92 + DN(2,2)*clhs94 + clhs291;
            lhs(8,11)=DN(2,0)*clhs294;
            lhs(8,12)=DN(2,0)*clhs99 + DN(2,1)*clhs101 + DN(2,2)*clhs103 + clhs295 + clhs298;
            lhs(8,13)=DN(2,0)*clhs110 + DN(2,1)*clhs112 + DN(2,2)*clhs115 + clhs299;
            lhs(8,14)=DN(2,0)*clhs117 + DN(2,1)*clhs119 + DN(2,2)*clhs121 + clhs300;
            lhs(8,15)=DN(3,0)*clhs293 + clhs301*clhs68 - clhs301 + clhs302*clhs70;
            lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs125 + DN(2,2)*clhs126 + clhs89;
            lhs(9,1)=DN(2,0)*clhs24 + DN(2,1)*clhs127 + DN(2,2)*clhs129 + clhs156 + clhs284;
            lhs(9,2)=DN(2,0)*clhs32 + DN(2,1)*clhs131 + DN(2,2)*clhs133 + clhs194;
            lhs(9,3)=clhs161*clhs70 + clhs162*clhs68 - clhs162 + clhs213*clhs79;
            lhs(9,4)=DN(2,0)*clhs45 + DN(2,1)*clhs136 + DN(2,2)*clhs137 + clhs238;
            lhs(9,5)=DN(2,0)*clhs56 + DN(2,1)*clhs139 + DN(2,2)*clhs141 + clhs253 + clhs285;
            lhs(9,6)=DN(2,0)*clhs63 + DN(2,1)*clhs143 + DN(2,2)*clhs145 + clhs265;
            lhs(9,7)=clhs149*clhs79 + clhs255*clhs70 + clhs256*clhs68 - clhs256;
            lhs(9,8)=DN(2,0)*clhs74 + DN(2,1)*clhs150 + DN(2,2)*clhs151 + clhs290;
            lhs(9,9)=DN(2,0)*clhs85 + DN(2,1)*clhs153 + DN(2,2)*clhs155 + clhs10*clhs303 + clhs288;
            lhs(9,10)=DN(2,0)*clhs92 + DN(2,1)*clhs157 + DN(2,2)*clhs159 + clhs305;
            lhs(9,11)=DN(2,1)*clhs294;
            lhs(9,12)=DN(2,0)*clhs101 + DN(2,1)*clhs164 + DN(2,2)*clhs165 + clhs306;
            lhs(9,13)=DN(2,0)*clhs112 + DN(2,1)*clhs167 + DN(2,2)*clhs169 + clhs298 + clhs307;
            lhs(9,14)=DN(2,0)*clhs119 + DN(2,1)*clhs171 + DN(2,2)*clhs173 + clhs308;
            lhs(9,15)=DN(3,1)*clhs293 + clhs309*clhs68 - clhs309 + clhs310*clhs70;
            lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs126 + DN(2,2)*clhs177 + clhs95;
            lhs(10,1)=DN(2,0)*clhs27 + DN(2,1)*clhs129 + DN(2,2)*clhs178 + clhs160;
            lhs(10,2)=DN(2,0)*clhs34 + DN(2,1)*clhs133 + DN(2,2)*clhs179 + clhs196 + clhs284;
            lhs(10,3)=clhs197*clhs70 + clhs198*clhs68 - clhs198 + clhs214*clhs79;
            lhs(10,4)=DN(2,0)*clhs47 + DN(2,1)*clhs137 + DN(2,2)*clhs181 + clhs239;
            lhs(10,5)=DN(2,0)*clhs59 + DN(2,1)*clhs141 + DN(2,2)*clhs184 + clhs254;
            lhs(10,6)=DN(2,0)*clhs65 + DN(2,1)*clhs145 + DN(2,2)*clhs186 + clhs266 + clhs285;
            lhs(10,7)=clhs190*clhs79 + clhs267*clhs70 + clhs268*clhs68 - clhs268;
            lhs(10,8)=DN(2,0)*clhs76 + DN(2,1)*clhs151 + DN(2,2)*clhs191 + clhs291;
            lhs(10,9)=DN(2,0)*clhs88 + DN(2,1)*clhs155 + DN(2,2)*clhs193 + clhs305;
            lhs(10,10)=DN(2,0)*clhs94 + DN(2,1)*clhs159 + DN(2,2)*clhs195 + clhs10*clhs311 + clhs288;
            lhs(10,11)=DN(2,2)*clhs294;
            lhs(10,12)=DN(2,0)*clhs103 + DN(2,1)*clhs165 + DN(2,2)*clhs200 + clhs313;
            lhs(10,13)=DN(2,0)*clhs115 + DN(2,1)*clhs169 + DN(2,2)*clhs202 + clhs314;
            lhs(10,14)=DN(2,0)*clhs121 + DN(2,1)*clhs173 + DN(2,2)*clhs204 + clhs298 + clhs315;
            lhs(10,15)=DN(3,2)*clhs293 + clhs316*clhs68 - clhs316 + clhs317*clhs70;
            lhs(11,0)=clhs15*clhs98 + clhs96;
            lhs(11,1)=clhs15*clhs163 + clhs161;
            lhs(11,2)=clhs15*clhs199 + clhs197;
            lhs(11,3)=clhs219;
            lhs(11,4)=clhs240 + clhs52*clhs98;
            lhs(11,5)=clhs163*clhs52 + clhs255;
            lhs(11,6)=clhs199*clhs52 + clhs267;
            lhs(11,7)=clhs279;
            lhs(11,8)=DN(2,0)*clhs319;
            lhs(11,9)=DN(2,1)*clhs319;
            lhs(11,10)=DN(2,2)*clhs319;
            lhs(11,11)=clhs210*clhs287 + clhs211*clhs286 + clhs211*clhs303 + clhs211*clhs311;
            lhs(11,12)=clhs108*clhs98 + clhs302;
            lhs(11,13)=clhs108*clhs163 + clhs310;
            lhs(11,14)=clhs108*clhs199 + clhs317;
            lhs(11,15)=clhs320;
            lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs104 + clhs324;
            lhs(12,1)=DN(3,0)*clhs22 + DN(3,1)*clhs24 + DN(3,2)*clhs27 + clhs166;
            lhs(12,2)=DN(3,0)*clhs30 + DN(3,1)*clhs32 + DN(3,2)*clhs34 + clhs201;
            lhs(12,3)=clhs106*clhs212 + clhs123*clhs70 + clhs124*clhs68 - clhs124;
            lhs(12,4)=DN(3,0)*clhs43 + DN(3,1)*clhs45 + DN(3,2)*clhs47 + clhs242 + clhs325;
            lhs(12,5)=DN(3,0)*clhs54 + DN(3,1)*clhs56 + DN(3,2)*clhs59 + clhs257;
            lhs(12,6)=DN(3,0)*clhs61 + DN(3,1)*clhs63 + DN(3,2)*clhs65 + clhs269;
            lhs(12,7)=clhs106*clhs71 + clhs247*clhs70 + clhs248*clhs68 - clhs248;
            lhs(12,8)=DN(3,0)*clhs72 + DN(3,1)*clhs74 + DN(3,2)*clhs76 + clhs295 + clhs326;
            lhs(12,9)=DN(3,0)*clhs83 + DN(3,1)*clhs85 + DN(3,2)*clhs88 + clhs306;
            lhs(12,10)=DN(3,0)*clhs90 + DN(3,1)*clhs92 + DN(3,2)*clhs94 + clhs313;
            lhs(12,11)=clhs106*clhs98 + clhs301*clhs70 + clhs302*clhs68 - clhs302;
            lhs(12,12)=DN(3,0)*clhs99 + DN(3,1)*clhs101 + DN(3,2)*clhs103 + clhs10*clhs327 + clhs329;
            lhs(12,13)=DN(3,0)*clhs110 + DN(3,1)*clhs112 + DN(3,2)*clhs115 + clhs331;
            lhs(12,14)=DN(3,0)*clhs117 + DN(3,1)*clhs119 + DN(3,2)*clhs121 + clhs332;
            lhs(12,15)=DN(3,0)*clhs333;
            lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs125 + DN(3,2)*clhs126 + clhs116;
            lhs(13,1)=DN(3,0)*clhs24 + DN(3,1)*clhs127 + DN(3,2)*clhs129 + clhs170 + clhs324;
            lhs(13,2)=DN(3,0)*clhs32 + DN(3,1)*clhs131 + DN(3,2)*clhs133 + clhs203;
            lhs(13,3)=clhs106*clhs213 + clhs175*clhs70 + clhs176*clhs68 - clhs176;
            lhs(13,4)=DN(3,0)*clhs45 + DN(3,1)*clhs136 + DN(3,2)*clhs137 + clhs245;
            lhs(13,5)=DN(3,0)*clhs56 + DN(3,1)*clhs139 + DN(3,2)*clhs141 + clhs258 + clhs325;
            lhs(13,6)=DN(3,0)*clhs63 + DN(3,1)*clhs143 + DN(3,2)*clhs145 + clhs270;
            lhs(13,7)=clhs106*clhs149 + clhs260*clhs70 + clhs261*clhs68 - clhs261;
            lhs(13,8)=DN(3,0)*clhs74 + DN(3,1)*clhs150 + DN(3,2)*clhs151 + clhs299;
            lhs(13,9)=DN(3,0)*clhs85 + DN(3,1)*clhs153 + DN(3,2)*clhs155 + clhs307 + clhs326;
            lhs(13,10)=DN(3,0)*clhs92 + DN(3,1)*clhs157 + DN(3,2)*clhs159 + clhs314;
            lhs(13,11)=clhs106*clhs163 + clhs309*clhs70 + clhs310*clhs68 - clhs310;
            lhs(13,12)=DN(3,0)*clhs101 + DN(3,1)*clhs164 + DN(3,2)*clhs165 + clhs331;
            lhs(13,13)=DN(3,0)*clhs112 + DN(3,1)*clhs167 + DN(3,2)*clhs169 + clhs10*clhs334 + clhs329;
            lhs(13,14)=DN(3,0)*clhs119 + DN(3,1)*clhs171 + DN(3,2)*clhs173 + clhs335;
            lhs(13,15)=DN(3,1)*clhs333;
            lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs126 + DN(3,2)*clhs177 + clhs122;
            lhs(14,1)=DN(3,0)*clhs27 + DN(3,1)*clhs129 + DN(3,2)*clhs178 + clhs174;
            lhs(14,2)=DN(3,0)*clhs34 + DN(3,1)*clhs133 + DN(3,2)*clhs179 + clhs205 + clhs324;
            lhs(14,3)=clhs106*clhs214 + clhs206*clhs70 + clhs207*clhs68 - clhs207;
            lhs(14,4)=DN(3,0)*clhs47 + DN(3,1)*clhs137 + DN(3,2)*clhs181 + clhs246;
            lhs(14,5)=DN(3,0)*clhs59 + DN(3,1)*clhs141 + DN(3,2)*clhs184 + clhs259;
            lhs(14,6)=DN(3,0)*clhs65 + DN(3,1)*clhs145 + DN(3,2)*clhs186 + clhs271 + clhs325;
            lhs(14,7)=clhs106*clhs190 + clhs272*clhs70 + clhs273*clhs68 - clhs273;
            lhs(14,8)=DN(3,0)*clhs76 + DN(3,1)*clhs151 + DN(3,2)*clhs191 + clhs300;
            lhs(14,9)=DN(3,0)*clhs88 + DN(3,1)*clhs155 + DN(3,2)*clhs193 + clhs308;
            lhs(14,10)=DN(3,0)*clhs94 + DN(3,1)*clhs159 + DN(3,2)*clhs195 + clhs315 + clhs326;
            lhs(14,11)=clhs106*clhs199 + clhs316*clhs70 + clhs317*clhs68 - clhs317;
            lhs(14,12)=DN(3,0)*clhs103 + DN(3,1)*clhs165 + DN(3,2)*clhs200 + clhs332;
            lhs(14,13)=DN(3,0)*clhs115 + DN(3,1)*clhs169 + DN(3,2)*clhs202 + clhs335;
            lhs(14,14)=DN(3,0)*clhs121 + DN(3,1)*clhs173 + DN(3,2)*clhs204 + clhs10*clhs336 + clhs329;
            lhs(14,15)=DN(3,2)*clhs333;
            lhs(15,0)=DN(3,0)*clhs208 + clhs123;
            lhs(15,1)=DN(3,1)*clhs208 + clhs175;
            lhs(15,2)=DN(3,2)*clhs208 + clhs206;
            lhs(15,3)=clhs220;
            lhs(15,4)=DN(3,0)*clhs274 + clhs247;
            lhs(15,5)=DN(3,1)*clhs274 + clhs260;
            lhs(15,6)=DN(3,2)*clhs274 + clhs272;
            lhs(15,7)=clhs280;
            lhs(15,8)=DN(3,0)*clhs318 + clhs301;
            lhs(15,9)=DN(3,1)*clhs318 + clhs309;
            lhs(15,10)=DN(3,2)*clhs318 + clhs316;
            lhs(15,11)=clhs320;
            lhs(15,12)=DN(3,0)*clhs337;
            lhs(15,13)=DN(3,1)*clhs337;
            lhs(15,14)=DN(3,2)*clhs337;
            lhs(15,15)=clhs210*clhs328 + clhs211*clhs327 + clhs211*clhs334 + clhs211*clhs336;


}


template<>
void NavierStokes<2>::ComputeGaussPointLHSContribution(
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
    const double& bdf0 = data.bdf0;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v = data.v;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& vconv = v - vmesh;

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
const double clhs9 =             bdf0*rho;
const double clhs10 =             N[0]*rho;
const double clhs11 =             DN(0,0)*clhs4 + DN(0,1)*clhs5;
const double clhs12 =             N[0]*bdf0 + clhs11;
const double clhs13 =             pow(rho, 2);
const double clhs14 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double clhs15 =             1.0/(clhs6/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs16 =             1.0*N[0]*clhs13*clhs14*clhs15;
const double clhs17 =             1.0*clhs11*clhs13*clhs15;
const double clhs18 =             clhs10*clhs11 + clhs12*clhs16 + clhs12*clhs17 + clhs8*clhs9;
const double clhs19 =             C(0,1)*DN(0,1) + clhs1;
const double clhs20 =             C(1,2)*DN(0,1);
const double clhs21 =             C(2,2)*DN(0,0) + clhs20;
const double clhs22 =             DN(0,0)*clhs7;
const double clhs23 =             DN(0,1)*clhs22;
const double clhs24 =             pow(c, -2);
const double clhs25 =             1.0/rho;
const double clhs26 =             N[0]*bdf0*clhs24*clhs25;
const double clhs27 =             1.0*clhs14*clhs15;
const double clhs28 =             1.0*clhs15*rho;
const double clhs29 =             clhs11*clhs28;
const double clhs30 =             -N[0] + clhs10*clhs27 + clhs26*clhs7 + clhs29;
const double clhs31 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs32 =             C(0,2)*DN(1,0);
const double clhs33 =             C(2,2)*DN(1,1) + clhs32;
const double clhs34 =             DN(1,0)*clhs22;
const double clhs35 =             N[0]*bdf0*rho;
const double clhs36 =             N[1]*clhs35;
const double clhs37 =             DN(1,0)*clhs4 + DN(1,1)*clhs5;
const double clhs38 =             N[1]*bdf0;
const double clhs39 =             clhs37 + clhs38;
const double clhs40 =             clhs10*clhs37 + clhs16*clhs39 + clhs17*clhs39 + clhs36;
const double clhs41 =             C(0,1)*DN(1,1) + clhs32;
const double clhs42 =             C(1,2)*DN(1,1);
const double clhs43 =             C(2,2)*DN(1,0) + clhs42;
const double clhs44 =             DN(1,1)*clhs22;
const double clhs45 =             DN(0,0)*N[1];
const double clhs46 =             bdf0*clhs24*clhs25*clhs7;
const double clhs47 =             DN(1,0)*N[0];
const double clhs48 =             1.0*clhs14*clhs15*rho;
const double clhs49 =             1.0*DN(1,0)*clhs15*rho;
const double clhs50 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs51 =             C(0,2)*DN(2,0);
const double clhs52 =             C(2,2)*DN(2,1) + clhs51;
const double clhs53 =             DN(2,0)*clhs22;
const double clhs54 =             N[2]*clhs35;
const double clhs55 =             DN(2,0)*clhs4 + DN(2,1)*clhs5;
const double clhs56 =             N[2]*bdf0;
const double clhs57 =             clhs55 + clhs56;
const double clhs58 =             clhs10*clhs55 + clhs16*clhs57 + clhs17*clhs57 + clhs54;
const double clhs59 =             C(0,1)*DN(2,1) + clhs51;
const double clhs60 =             C(1,2)*DN(2,1);
const double clhs61 =             C(2,2)*DN(2,0) + clhs60;
const double clhs62 =             DN(2,1)*clhs22;
const double clhs63 =             DN(0,0)*N[2];
const double clhs64 =             DN(2,0)*N[0];
const double clhs65 =             C(0,1)*DN(0,0) + clhs20;
const double clhs66 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs67 =             pow(DN(0,1), 2);
const double clhs68 =             C(0,1)*DN(1,0) + clhs42;
const double clhs69 =             DN(0,1)*clhs7;
const double clhs70 =             DN(1,0)*clhs69;
const double clhs71 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs72 =             DN(1,1)*clhs69;
const double clhs73 =             DN(0,1)*N[1];
const double clhs74 =             DN(1,1)*N[0];
const double clhs75 =             1.0*DN(1,1)*clhs15*rho;
const double clhs76 =             C(0,1)*DN(2,0) + clhs60;
const double clhs77 =             DN(2,0)*clhs69;
const double clhs78 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs79 =             DN(2,1)*clhs69;
const double clhs80 =             DN(0,1)*N[2];
const double clhs81 =             DN(2,1)*N[0];
const double clhs82 =             clhs12*clhs28;
const double clhs83 =             N[0] + clhs82;
const double clhs84 =             bdf0*clhs24*clhs25;
const double clhs85 =             1.0*clhs15;
const double clhs86 =             1.0*DN(0,0)*clhs15*rho;
const double clhs87 =             1.0*DN(0,1)*clhs15*rho;
const double clhs88 =             1.0*DN(0,0)*clhs15;
const double clhs89 =             1.0*DN(0,1)*clhs15;
const double clhs90 =             DN(1,0)*clhs88 + DN(1,1)*clhs89 + N[1]*clhs26;
const double clhs91 =             DN(2,0)*clhs88 + DN(2,1)*clhs89 + N[2]*clhs26;
const double clhs92 =             N[1]*rho;
const double clhs93 =             1.0*N[1]*clhs13*clhs14*clhs15;
const double clhs94 =             1.0*clhs13*clhs15*clhs37;
const double clhs95 =             clhs11*clhs92 + clhs12*clhs93 + clhs12*clhs94 + clhs36;
const double clhs96 =             pow(DN(1,0), 2);
const double clhs97 =             pow(N[1], 2);
const double clhs98 =             clhs37*clhs92 + clhs39*clhs93 + clhs39*clhs94 + clhs9*clhs97;
const double clhs99 =             DN(1,0)*clhs7;
const double clhs100 =             DN(1,1)*clhs99;
const double clhs101 =             clhs24*clhs25*clhs7;
const double clhs102 =             clhs28*clhs37;
const double clhs103 =             -N[1] + clhs101*clhs38 + clhs102 + clhs27*clhs92;
const double clhs104 =             DN(2,0)*clhs99;
const double clhs105 =             N[1]*N[2]*bdf0;
const double clhs106 =             clhs105*rho;
const double clhs107 =             clhs106 + clhs55*clhs92 + clhs57*clhs93 + clhs57*clhs94;
const double clhs108 =             DN(2,1)*clhs99;
const double clhs109 =             DN(1,0)*N[2];
const double clhs110 =             DN(2,0)*N[1];
const double clhs111 =             pow(DN(1,1), 2);
const double clhs112 =             DN(1,1)*clhs7;
const double clhs113 =             DN(2,0)*clhs112;
const double clhs114 =             DN(2,1)*clhs112;
const double clhs115 =             DN(1,1)*N[2];
const double clhs116 =             DN(2,1)*N[1];
const double clhs117 =             clhs28*clhs39;
const double clhs118 =             N[1] + clhs117;
const double clhs119 =             1.0*DN(1,0)*DN(2,0)*clhs15 + 1.0*DN(1,1)*DN(2,1)*clhs15 + clhs105*clhs24*clhs25;
const double clhs120 =             N[2]*rho;
const double clhs121 =             1.0*N[2]*clhs13*clhs14*clhs15;
const double clhs122 =             1.0*clhs13*clhs15*clhs55;
const double clhs123 =             clhs11*clhs120 + clhs12*clhs121 + clhs12*clhs122 + clhs54;
const double clhs124 =             clhs106 + clhs120*clhs37 + clhs121*clhs39 + clhs122*clhs39;
const double clhs125 =             pow(DN(2,0), 2);
const double clhs126 =             pow(N[2], 2);
const double clhs127 =             clhs120*clhs55 + clhs121*clhs57 + clhs122*clhs57 + clhs126*clhs9;
const double clhs128 =             DN(2,0)*DN(2,1)*clhs7;
const double clhs129 =             -N[2] + clhs101*clhs56 + clhs120*clhs27 + clhs28*clhs55;
const double clhs130 =             pow(DN(2,1), 2);
const double clhs131 =             N[2] + clhs28*clhs57;
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs18 + clhs3*clhs7;
            lhs(0,1)=DN(0,0)*clhs19 + DN(0,1)*clhs21 + clhs23;
            lhs(0,2)=DN(0,0)*clhs30;
            lhs(0,3)=DN(0,0)*clhs31 + DN(0,1)*clhs33 + clhs34 + clhs40;
            lhs(0,4)=DN(0,0)*clhs41 + DN(0,1)*clhs43 + clhs44;
            lhs(0,5)=clhs11*clhs49 + clhs45*clhs46 - clhs45 + clhs47*clhs48;
            lhs(0,6)=DN(0,0)*clhs50 + DN(0,1)*clhs52 + clhs53 + clhs58;
            lhs(0,7)=DN(0,0)*clhs59 + DN(0,1)*clhs61 + clhs62;
            lhs(0,8)=DN(2,0)*clhs29 + clhs46*clhs63 + clhs48*clhs64 - clhs63;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs65 + clhs23;
            lhs(1,1)=DN(0,0)*clhs21 + DN(0,1)*clhs66 + clhs18 + clhs67*clhs7;
            lhs(1,2)=DN(0,1)*clhs30;
            lhs(1,3)=DN(0,0)*clhs33 + DN(0,1)*clhs68 + clhs70;
            lhs(1,4)=DN(0,0)*clhs43 + DN(0,1)*clhs71 + clhs40 + clhs72;
            lhs(1,5)=clhs11*clhs75 + clhs46*clhs73 + clhs48*clhs74 - clhs73;
            lhs(1,6)=DN(0,0)*clhs52 + DN(0,1)*clhs76 + clhs77;
            lhs(1,7)=DN(0,0)*clhs61 + DN(0,1)*clhs78 + clhs58 + clhs79;
            lhs(1,8)=DN(2,1)*clhs29 + clhs46*clhs80 + clhs48*clhs81 - clhs80;
            lhs(2,0)=DN(0,0)*clhs83;
            lhs(2,1)=DN(0,1)*clhs83;
            lhs(2,2)=clhs3*clhs85 + clhs67*clhs85 + clhs8*clhs84;
            lhs(2,3)=clhs39*clhs86 + clhs47;
            lhs(2,4)=clhs39*clhs87 + clhs74;
            lhs(2,5)=clhs90;
            lhs(2,6)=clhs57*clhs86 + clhs64;
            lhs(2,7)=clhs57*clhs87 + clhs81;
            lhs(2,8)=clhs91;
            lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs34 + clhs95;
            lhs(3,1)=DN(1,0)*clhs19 + DN(1,1)*clhs21 + clhs70;
            lhs(3,2)=clhs37*clhs86 + clhs45*clhs48 + clhs46*clhs47 - clhs47;
            lhs(3,3)=DN(1,0)*clhs31 + DN(1,1)*clhs33 + clhs7*clhs96 + clhs98;
            lhs(3,4)=DN(1,0)*clhs41 + DN(1,1)*clhs43 + clhs100;
            lhs(3,5)=DN(1,0)*clhs103;
            lhs(3,6)=DN(1,0)*clhs50 + DN(1,1)*clhs52 + clhs104 + clhs107;
            lhs(3,7)=DN(1,0)*clhs59 + DN(1,1)*clhs61 + clhs108;
            lhs(3,8)=DN(2,0)*clhs102 + clhs109*clhs46 - clhs109 + clhs110*clhs48;
            lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs65 + clhs44;
            lhs(4,1)=DN(1,0)*clhs21 + DN(1,1)*clhs66 + clhs72 + clhs95;
            lhs(4,2)=clhs37*clhs87 + clhs46*clhs74 + clhs48*clhs73 - clhs74;
            lhs(4,3)=DN(1,0)*clhs33 + DN(1,1)*clhs68 + clhs100;
            lhs(4,4)=DN(1,0)*clhs43 + DN(1,1)*clhs71 + clhs111*clhs7 + clhs98;
            lhs(4,5)=DN(1,1)*clhs103;
            lhs(4,6)=DN(1,0)*clhs52 + DN(1,1)*clhs76 + clhs113;
            lhs(4,7)=DN(1,0)*clhs61 + DN(1,1)*clhs78 + clhs107 + clhs114;
            lhs(4,8)=DN(2,1)*clhs102 + clhs115*clhs46 - clhs115 + clhs116*clhs48;
            lhs(5,0)=clhs12*clhs49 + clhs45;
            lhs(5,1)=clhs12*clhs75 + clhs73;
            lhs(5,2)=clhs90;
            lhs(5,3)=DN(1,0)*clhs118;
            lhs(5,4)=DN(1,1)*clhs118;
            lhs(5,5)=clhs111*clhs85 + clhs84*clhs97 + clhs85*clhs96;
            lhs(5,6)=clhs110 + clhs49*clhs57;
            lhs(5,7)=clhs116 + clhs57*clhs75;
            lhs(5,8)=clhs119;
            lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs123 + clhs53;
            lhs(6,1)=DN(2,0)*clhs19 + DN(2,1)*clhs21 + clhs77;
            lhs(6,2)=clhs46*clhs64 + clhs48*clhs63 + clhs55*clhs86 - clhs64;
            lhs(6,3)=DN(2,0)*clhs31 + DN(2,1)*clhs33 + clhs104 + clhs124;
            lhs(6,4)=DN(2,0)*clhs41 + DN(2,1)*clhs43 + clhs113;
            lhs(6,5)=clhs109*clhs48 + clhs110*clhs46 - clhs110 + clhs49*clhs55;
            lhs(6,6)=DN(2,0)*clhs50 + DN(2,1)*clhs52 + clhs125*clhs7 + clhs127;
            lhs(6,7)=DN(2,0)*clhs59 + DN(2,1)*clhs61 + clhs128;
            lhs(6,8)=DN(2,0)*clhs129;
            lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs65 + clhs62;
            lhs(7,1)=DN(2,0)*clhs21 + DN(2,1)*clhs66 + clhs123 + clhs79;
            lhs(7,2)=clhs46*clhs81 + clhs48*clhs80 + clhs55*clhs87 - clhs81;
            lhs(7,3)=DN(2,0)*clhs33 + DN(2,1)*clhs68 + clhs108;
            lhs(7,4)=DN(2,0)*clhs43 + DN(2,1)*clhs71 + clhs114 + clhs124;
            lhs(7,5)=clhs115*clhs48 + clhs116*clhs46 - clhs116 + clhs55*clhs75;
            lhs(7,6)=DN(2,0)*clhs52 + DN(2,1)*clhs76 + clhs128;
            lhs(7,7)=DN(2,0)*clhs61 + DN(2,1)*clhs78 + clhs127 + clhs130*clhs7;
            lhs(7,8)=DN(2,1)*clhs129;
            lhs(8,0)=DN(2,0)*clhs82 + clhs63;
            lhs(8,1)=DN(2,1)*clhs82 + clhs80;
            lhs(8,2)=clhs91;
            lhs(8,3)=DN(2,0)*clhs117 + clhs109;
            lhs(8,4)=DN(2,1)*clhs117 + clhs115;
            lhs(8,5)=clhs119;
            lhs(8,6)=DN(2,0)*clhs131;
            lhs(8,7)=DN(2,1)*clhs131;
            lhs(8,8)=clhs125*clhs85 + clhs126*clhs84 + clhs130*clhs85;


}


template<>
void NavierStokes<3>::ComputeGaussPointRHSContribution(
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
    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v = data.v;
    const BoundedMatrix<double,nnodes,dim>& vn = data.vn;
    const BoundedMatrix<double,nnodes,dim>& vnn = data.vnn;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& vconv = v - vmesh;
    const BoundedMatrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& p = data.p;
    const array_1d<double,nnodes>& pn = data.pn;
    const array_1d<double,nnodes>& pnn = data.pnn;
    const array_1d<double,strain_size>& stress = data.stress;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,dim> f_gauss = prod(trans(f), N);
    const array_1d<double,dim> grad_p = prod(trans(DN), p);

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crhs2 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crhs3 =             DN(0,0)*v(0,0);
const double crhs4 =             DN(1,0)*v(1,0);
const double crhs5 =             DN(2,0)*v(2,0);
const double crhs6 =             DN(3,0)*v(3,0);
const double crhs7 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs8 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs9 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs10 =             rho*(crhs7*(crhs3 + crhs4 + crhs5 + crhs6) + crhs8*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs9*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)));
const double crhs11 =             rho*stab_c2*sqrt(pow(crhs7, 2) + pow(crhs8, 2) + pow(crhs9, 2));
const double crhs12 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs13 =             DN(0,1)*v(0,1);
const double crhs14 =             DN(1,1)*v(1,1);
const double crhs15 =             DN(2,1)*v(2,1);
const double crhs16 =             DN(3,1)*v(3,1);
const double crhs17 =             crhs12 + crhs13 + crhs14 + crhs15 + crhs16 + crhs3 + crhs4 + crhs5 + crhs6;
const double crhs18 =             (N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]) + N[3]*(bdf0*p[3] + bdf1*pn[3] + bdf2*pnn[3]))/(pow(c, 2)*rho);
const double crhs19 =             (crhs17 + crhs18)*(crhs11*h/stab_c1 + mu);
const double crhs20 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crhs21 =             N[0]*crhs20*rho;
const double crhs22 =             1.0/(crhs11/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs23 =             1.0*crhs22*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + crhs10 + crhs2);
const double crhs24 =             rho*(DN(0,0)*crhs7 + DN(0,1)*crhs8 + DN(0,2)*crhs9);
const double crhs25 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs26 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crhs27 =             rho*(crhs7*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs8*(crhs13 + crhs14 + crhs15 + crhs16) + crhs9*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
const double crhs28 =             1.0*crhs22*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs25 + crhs26 + crhs27);
const double crhs29 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs30 =             rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crhs31 =             rho*(crhs12*crhs9 + crhs7*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs8*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)));
const double crhs32 =             1.0*crhs22*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs29 + crhs30 + crhs31);
const double crhs33 =             N[1]*crhs20*rho;
const double crhs34 =             rho*(DN(1,0)*crhs7 + DN(1,1)*crhs8 + DN(1,2)*crhs9);
const double crhs35 =             N[2]*crhs20*rho;
const double crhs36 =             rho*(DN(2,0)*crhs7 + DN(2,1)*crhs8 + DN(2,2)*crhs9);
const double crhs37 =             N[3]*crhs20*rho;
const double crhs38 =             rho*(DN(3,0)*crhs7 + DN(3,1)*crhs8 + DN(3,2)*crhs9);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs19 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - N[0]*crhs10 - N[0]*crhs2 - crhs21*crhs23 - crhs23*crhs24;
            rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs19 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs25 - N[0]*crhs26 - N[0]*crhs27 - crhs21*crhs28 - crhs24*crhs28;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs19 - DN(0,2)*stress[2] + N[0]*crhs29 - N[0]*crhs30 - N[0]*crhs31 - crhs21*crhs32 - crhs24*crhs32;
            rhs[3]=-DN(0,0)*crhs23 - DN(0,1)*crhs28 - DN(0,2)*crhs32 - N[0]*crhs17 - N[0]*crhs18;
            rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs19 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - N[1]*crhs10 - N[1]*crhs2 - crhs23*crhs33 - crhs23*crhs34;
            rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs19 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs25 - N[1]*crhs26 - N[1]*crhs27 - crhs28*crhs33 - crhs28*crhs34;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs19 - DN(1,2)*stress[2] + N[1]*crhs29 - N[1]*crhs30 - N[1]*crhs31 - crhs32*crhs33 - crhs32*crhs34;
            rhs[7]=-DN(1,0)*crhs23 - DN(1,1)*crhs28 - DN(1,2)*crhs32 - N[1]*crhs17 - N[1]*crhs18;
            rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs19 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - N[2]*crhs10 - N[2]*crhs2 - crhs23*crhs35 - crhs23*crhs36;
            rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs19 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs25 - N[2]*crhs26 - N[2]*crhs27 - crhs28*crhs35 - crhs28*crhs36;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs19 - DN(2,2)*stress[2] + N[2]*crhs29 - N[2]*crhs30 - N[2]*crhs31 - crhs32*crhs35 - crhs32*crhs36;
            rhs[11]=-DN(2,0)*crhs23 - DN(2,1)*crhs28 - DN(2,2)*crhs32 - N[2]*crhs17 - N[2]*crhs18;
            rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs19 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - N[3]*crhs10 - N[3]*crhs2 - crhs23*crhs37 - crhs23*crhs38;
            rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs19 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs25 - N[3]*crhs26 - N[3]*crhs27 - crhs28*crhs37 - crhs28*crhs38;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs19 - DN(3,2)*stress[2] + N[3]*crhs29 - N[3]*crhs30 - N[3]*crhs31 - crhs32*crhs37 - crhs32*crhs38;
            rhs[15]=-DN(3,0)*crhs23 - DN(3,1)*crhs28 - DN(3,2)*crhs32 - N[3]*crhs17 - N[3]*crhs18;


}


template<>
void NavierStokes<2>::ComputeGaussPointRHSContribution(
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
    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v = data.v;
    const BoundedMatrix<double,nnodes,dim>& vn = data.vn;
    const BoundedMatrix<double,nnodes,dim>& vnn = data.vnn;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& vconv = v - vmesh;
    const BoundedMatrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& p = data.p;
    const array_1d<double,nnodes>& pn = data.pn;
    const array_1d<double,nnodes>& pnn = data.pnn;
    const array_1d<double,strain_size>& stress = data.stress;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,dim> f_gauss = prod(trans(f), N);
    const array_1d<double,dim> grad_p = prod(trans(DN), p);

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crhs2 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crhs3 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs4 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs5 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs6 =             rho*(crhs3*crhs4 + crhs5*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crhs7 =             rho*stab_c2*sqrt(pow(crhs4, 2) + pow(crhs5, 2));
const double crhs8 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs9 =             crhs3 + crhs8;
const double crhs10 =             (N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]))/(pow(c, 2)*rho);
const double crhs11 =             (crhs10 + crhs9)*(crhs7*h/stab_c1 + mu);
const double crhs12 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crhs13 =             N[0]*crhs12*rho;
const double crhs14 =             1.0/(crhs7/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs15 =             1.0*crhs14*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + crhs2 + crhs6);
const double crhs16 =             rho*(DN(0,0)*crhs4 + DN(0,1)*crhs5);
const double crhs17 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs18 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crhs19 =             rho*(crhs4*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + crhs5*crhs8);
const double crhs20 =             1.0*crhs14*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs17 + crhs18 + crhs19);
const double crhs21 =             N[1]*crhs12*rho;
const double crhs22 =             rho*(DN(1,0)*crhs4 + DN(1,1)*crhs5);
const double crhs23 =             N[2]*crhs12*rho;
const double crhs24 =             rho*(DN(2,0)*crhs4 + DN(2,1)*crhs5);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs11 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - N[0]*crhs2 - N[0]*crhs6 - crhs13*crhs15 - crhs15*crhs16;
            rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs11 - DN(0,1)*stress[1] + N[0]*crhs17 - N[0]*crhs18 - N[0]*crhs19 - crhs13*crhs20 - crhs16*crhs20;
            rhs[2]=-DN(0,0)*crhs15 - DN(0,1)*crhs20 - N[0]*crhs10 - N[0]*crhs9;
            rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs11 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - N[1]*crhs2 - N[1]*crhs6 - crhs15*crhs21 - crhs15*crhs22;
            rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs11 - DN(1,1)*stress[1] + N[1]*crhs17 - N[1]*crhs18 - N[1]*crhs19 - crhs20*crhs21 - crhs20*crhs22;
            rhs[5]=-DN(1,0)*crhs15 - DN(1,1)*crhs20 - N[1]*crhs10 - N[1]*crhs9;
            rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs11 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - N[2]*crhs2 - N[2]*crhs6 - crhs15*crhs23 - crhs15*crhs24;
            rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs11 - DN(2,1)*stress[1] + N[2]*crhs17 - N[2]*crhs18 - N[2]*crhs19 - crhs20*crhs23 - crhs20*crhs24;
            rhs[8]=-DN(2,0)*crhs15 - DN(2,1)*crhs20 - N[2]*crhs10 - N[2]*crhs9;


}


template<>
double NavierStokes<3>::SubscaleErrorEstimate(const ElementDataStruct& data)
{
    constexpr int dim = 3;
    constexpr int nnodes = 4;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size
    // const double c = data.c;                                // Wave velocity

    const double& dt = data.dt;
    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v = data.v;
    const BoundedMatrix<double,nnodes,dim>& vn = data.vn;
    const BoundedMatrix<double,nnodes,dim>& vnn = data.vnn;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& vconv = v - vmesh;
    const BoundedMatrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& p = data.p;

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
            v_s_gauss[0]=-cv_s_gauss3*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + rho*(-N[0]*f(0,0) + N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) - N[1]*f(1,0) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) - N[2]*f(2,0) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) - N[3]*f(3,0) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)) + cv_s_gauss0*(DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0)) + cv_s_gauss1*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + cv_s_gauss2*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0))));
            v_s_gauss[1]=-cv_s_gauss3*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + rho*(-N[0]*f(0,1) + N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) - N[1]*f(1,1) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) - N[2]*f(2,1) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) - N[3]*f(3,1) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)) + cv_s_gauss0*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + cv_s_gauss1*(DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1)) + cv_s_gauss2*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1))));
            v_s_gauss[2]=-cv_s_gauss3*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + rho*(-N[0]*f(0,2) + N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) - N[1]*f(1,2) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) - N[2]*f(2,2) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) - N[3]*f(3,2) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)) + cv_s_gauss0*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + cv_s_gauss1*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + cv_s_gauss2*(DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2))));


    const double v_gauss_norm = norm_2(v_gauss);
    const double v_s_gauss_norm = norm_2(v_s_gauss);

    return v_s_gauss_norm/v_gauss_norm;
}


template<>
double NavierStokes<2>::SubscaleErrorEstimate(const ElementDataStruct& data)
{
    constexpr int dim = 2;
    constexpr int nnodes = 3;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size

    const double& dt = data.dt;
    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v = data.v;
    const BoundedMatrix<double,nnodes,dim>& vn = data.vn;
    const BoundedMatrix<double,nnodes,dim>& vnn = data.vnn;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& vconv = v - vmesh;
    const BoundedMatrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& p = data.p;

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
            v_s_gauss[0]=-cv_s_gauss2*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + rho*(-N[0]*f(0,0) + N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) - N[1]*f(1,0) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) - N[2]*f(2,0) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + cv_s_gauss0*(DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0)) + cv_s_gauss1*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0))));
            v_s_gauss[1]=-cv_s_gauss2*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + rho*(-N[0]*f(0,1) + N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) - N[1]*f(1,1) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) - N[2]*f(2,1) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + cv_s_gauss0*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + cv_s_gauss1*(DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1))));


    const double v_gauss_norm = norm_2(v_gauss);
    const double v_s_gauss_norm = norm_2(v_s_gauss);

    return v_s_gauss_norm/v_gauss_norm;
}

}
