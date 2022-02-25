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

#include "custom_elements/embedded_ausas_navier_stokes.h"

namespace Kratos {

template<>
void EmbeddedAusasNavierStokes<3>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const 
{
    KRATOS_TRY

    constexpr unsigned int dim = 3;
    constexpr unsigned int num_nodes = 4;
    constexpr unsigned int dof_size  = num_nodes*(dim+1);

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
void EmbeddedAusasNavierStokes<2>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const 
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
void EmbeddedAusasNavierStokes<3>::GetDofList(
    DofsVectorType& ElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const 
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
void EmbeddedAusasNavierStokes<2>::GetDofList(
    DofsVectorType& ElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const 
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
void EmbeddedAusasNavierStokes<3>::ComputeGaussPointLHSContribution(
    BoundedMatrix<double,16,16>& lhs,
    const EmbeddedAusasElementDataStruct& data)
{
    constexpr unsigned int dim = 3;
    constexpr unsigned int nnodes = 4;

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

    const double clhs0 = C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs1 = C(0,3)*DN(0,0);
const double clhs2 = C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs1;
const double clhs3 = C(0,5)*DN(0,0);
const double clhs4 = C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs3;
const double clhs5 = pow(DN(0,0), 2);
const double clhs6 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double clhs7 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double clhs8 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double clhs9 = rho*stab_c2*sqrt(pow(clhs6, 2) + pow(clhs7, 2) + pow(clhs8, 2));
const double clhs10 = clhs9*h/stab_c1 + mu;
const double clhs11 = DN(0,0)*clhs6 + DN(0,1)*clhs7 + DN(0,2)*clhs8;
const double clhs12 = N[0]*rho;
const double clhs13 = pow(N[0], 2);
const double clhs14 = bdf0*rho;
const double clhs15 = N[0]*bdf0;
const double clhs16 = clhs11 + clhs15;
const double clhs17 = 1.0/(clhs9/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs18 = clhs17*pow(rho, 2);
const double clhs19 = clhs11*clhs18;
const double clhs20 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double clhs21 = clhs18*clhs20;
const double clhs22 = N[0]*clhs21;
const double clhs23 = clhs11*clhs12 + clhs13*clhs14 + clhs16*clhs19 + clhs16*clhs22;
const double clhs24 = C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs1;
const double clhs25 = C(1,3)*DN(0,1);
const double clhs26 = C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs25;
const double clhs27 = C(3,5)*DN(0,0);
const double clhs28 = C(4,5)*DN(0,2);
const double clhs29 = C(1,5)*DN(0,1) + clhs27 + clhs28;
const double clhs30 = DN(0,0)*clhs10;
const double clhs31 = DN(0,1)*clhs30;
const double clhs32 = C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs3;
const double clhs33 = C(3,4)*DN(0,1);
const double clhs34 = C(2,3)*DN(0,2) + clhs27 + clhs33;
const double clhs35 = C(2,5)*DN(0,2);
const double clhs36 = C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs35;
const double clhs37 = DN(0,2)*clhs30;
const double clhs38 = -N[0];
const double clhs39 = 1/(pow(c, 2)*rho);
const double clhs40 = clhs15*clhs39;
const double clhs41 = clhs10*clhs40;
const double clhs42 = clhs17*clhs20;
const double clhs43 = clhs17*rho;
const double clhs44 = clhs11*clhs43;
const double clhs45 = clhs12*clhs42 + clhs38 + clhs41 + clhs44;
const double clhs46 = C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs47 = C(0,3)*DN(1,0);
const double clhs48 = C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs47;
const double clhs49 = C(0,5)*DN(1,0);
const double clhs50 = C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs49;
const double clhs51 = N[1]*rho;
const double clhs52 = clhs15*clhs51;
const double clhs53 = DN(1,0)*clhs30 + clhs52;
const double clhs54 = DN(1,0)*clhs6 + DN(1,1)*clhs7 + DN(1,2)*clhs8;
const double clhs55 = N[1]*bdf0;
const double clhs56 = clhs54 + clhs55;
const double clhs57 = clhs12*clhs54 + clhs19*clhs56 + clhs22*clhs56;
const double clhs58 = C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs47;
const double clhs59 = C(1,3)*DN(1,1);
const double clhs60 = C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs59;
const double clhs61 = C(3,5)*DN(1,0);
const double clhs62 = C(4,5)*DN(1,2);
const double clhs63 = C(1,5)*DN(1,1) + clhs61 + clhs62;
const double clhs64 = DN(1,1)*clhs30;
const double clhs65 = C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs49;
const double clhs66 = C(3,4)*DN(1,1);
const double clhs67 = C(2,3)*DN(1,2) + clhs61 + clhs66;
const double clhs68 = C(2,5)*DN(1,2);
const double clhs69 = C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs68;
const double clhs70 = DN(1,2)*clhs30;
const double clhs71 = DN(0,0)*N[1];
const double clhs72 = clhs39*clhs55;
const double clhs73 = DN(1,0)*N[0];
const double clhs74 = clhs20*clhs43;
const double clhs75 = C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs76 = C(0,3)*DN(2,0);
const double clhs77 = C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs76;
const double clhs78 = C(0,5)*DN(2,0);
const double clhs79 = C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs78;
const double clhs80 = N[2]*rho;
const double clhs81 = clhs15*clhs80;
const double clhs82 = DN(2,0)*clhs30 + clhs81;
const double clhs83 = DN(2,0)*clhs6 + DN(2,1)*clhs7 + DN(2,2)*clhs8;
const double clhs84 = N[2]*bdf0;
const double clhs85 = clhs83 + clhs84;
const double clhs86 = clhs12*clhs83 + clhs19*clhs85 + clhs22*clhs85;
const double clhs87 = C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs76;
const double clhs88 = C(1,3)*DN(2,1);
const double clhs89 = C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs88;
const double clhs90 = C(3,5)*DN(2,0);
const double clhs91 = C(4,5)*DN(2,2);
const double clhs92 = C(1,5)*DN(2,1) + clhs90 + clhs91;
const double clhs93 = DN(2,1)*clhs30;
const double clhs94 = C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs78;
const double clhs95 = C(3,4)*DN(2,1);
const double clhs96 = C(2,3)*DN(2,2) + clhs90 + clhs95;
const double clhs97 = C(2,5)*DN(2,2);
const double clhs98 = C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs97;
const double clhs99 = DN(2,2)*clhs30;
const double clhs100 = DN(0,0)*N[2];
const double clhs101 = clhs39*clhs84;
const double clhs102 = DN(2,0)*N[0];
const double clhs103 = C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs104 = C(0,3)*DN(3,0);
const double clhs105 = C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs104;
const double clhs106 = C(0,5)*DN(3,0);
const double clhs107 = C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs106;
const double clhs108 = N[3]*rho;
const double clhs109 = clhs108*clhs15;
const double clhs110 = DN(3,0)*clhs30 + clhs109;
const double clhs111 = DN(3,0)*clhs6 + DN(3,1)*clhs7 + DN(3,2)*clhs8;
const double clhs112 = N[3]*bdf0;
const double clhs113 = clhs111 + clhs112;
const double clhs114 = clhs111*clhs12 + clhs113*clhs19 + clhs113*clhs22;
const double clhs115 = C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs104;
const double clhs116 = C(1,3)*DN(3,1);
const double clhs117 = C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs116;
const double clhs118 = C(3,5)*DN(3,0);
const double clhs119 = C(4,5)*DN(3,2);
const double clhs120 = C(1,5)*DN(3,1) + clhs118 + clhs119;
const double clhs121 = DN(3,1)*clhs30;
const double clhs122 = C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs106;
const double clhs123 = C(3,4)*DN(3,1);
const double clhs124 = C(2,3)*DN(3,2) + clhs118 + clhs123;
const double clhs125 = C(2,5)*DN(3,2);
const double clhs126 = C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs125;
const double clhs127 = DN(3,2)*clhs30;
const double clhs128 = DN(0,0)*N[3];
const double clhs129 = clhs112*clhs39;
const double clhs130 = DN(3,0)*N[0];
const double clhs131 = C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs25;
const double clhs132 = C(0,4)*DN(0,0) + clhs28 + clhs33;
const double clhs133 = C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs134 = C(1,4)*DN(0,1);
const double clhs135 = C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs134;
const double clhs136 = pow(DN(0,1), 2);
const double clhs137 = C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs134;
const double clhs138 = C(2,4)*DN(0,2);
const double clhs139 = C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs138;
const double clhs140 = DN(0,1)*clhs10;
const double clhs141 = DN(0,2)*clhs140;
const double clhs142 = C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs59;
const double clhs143 = C(0,4)*DN(1,0) + clhs62 + clhs66;
const double clhs144 = DN(1,0)*clhs140;
const double clhs145 = C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs146 = C(1,4)*DN(1,1);
const double clhs147 = C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs146;
const double clhs148 = DN(1,1)*clhs140;
const double clhs149 = clhs52 + clhs57;
const double clhs150 = C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs146;
const double clhs151 = C(2,4)*DN(1,2);
const double clhs152 = C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs151;
const double clhs153 = DN(1,2)*clhs140;
const double clhs154 = DN(0,1)*N[1];
const double clhs155 = DN(1,1)*N[0];
const double clhs156 = C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs88;
const double clhs157 = C(0,4)*DN(2,0) + clhs91 + clhs95;
const double clhs158 = DN(2,0)*clhs140;
const double clhs159 = C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs160 = C(1,4)*DN(2,1);
const double clhs161 = C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs160;
const double clhs162 = DN(2,1)*clhs140;
const double clhs163 = clhs81 + clhs86;
const double clhs164 = C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs160;
const double clhs165 = C(2,4)*DN(2,2);
const double clhs166 = C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs165;
const double clhs167 = DN(2,2)*clhs140;
const double clhs168 = DN(0,1)*N[2];
const double clhs169 = DN(2,1)*N[0];
const double clhs170 = C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs116;
const double clhs171 = C(0,4)*DN(3,0) + clhs119 + clhs123;
const double clhs172 = DN(3,0)*clhs140;
const double clhs173 = C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs174 = C(1,4)*DN(3,1);
const double clhs175 = C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs174;
const double clhs176 = DN(3,1)*clhs140;
const double clhs177 = clhs109 + clhs114;
const double clhs178 = C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs174;
const double clhs179 = C(2,4)*DN(3,2);
const double clhs180 = C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs179;
const double clhs181 = DN(3,2)*clhs140;
const double clhs182 = DN(0,1)*N[3];
const double clhs183 = DN(3,1)*N[0];
const double clhs184 = C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs35;
const double clhs185 = C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs138;
const double clhs186 = C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs187 = pow(DN(0,2), 2);
const double clhs188 = C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs68;
const double clhs189 = DN(0,2)*clhs10;
const double clhs190 = DN(1,0)*clhs189;
const double clhs191 = C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs151;
const double clhs192 = DN(1,1)*clhs189;
const double clhs193 = C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs194 = DN(1,2)*clhs189;
const double clhs195 = DN(0,2)*N[1];
const double clhs196 = DN(1,2)*N[0];
const double clhs197 = C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs97;
const double clhs198 = DN(2,0)*clhs189;
const double clhs199 = C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs165;
const double clhs200 = DN(2,1)*clhs189;
const double clhs201 = C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs202 = DN(2,2)*clhs189;
const double clhs203 = DN(0,2)*N[2];
const double clhs204 = DN(2,2)*N[0];
const double clhs205 = C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs125;
const double clhs206 = DN(3,0)*clhs189;
const double clhs207 = C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs179;
const double clhs208 = DN(3,1)*clhs189;
const double clhs209 = C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs210 = DN(3,2)*clhs189;
const double clhs211 = DN(0,2)*N[3];
const double clhs212 = DN(3,2)*N[0];
const double clhs213 = clhs16*clhs43 + clhs38;
const double clhs214 = bdf0*clhs39;
const double clhs215 = -N[1];
const double clhs216 = clhs215 + clhs43*clhs56;
const double clhs217 = DN(0,0)*clhs17;
const double clhs218 = DN(0,1)*clhs17;
const double clhs219 = DN(0,2)*clhs17;
const double clhs220 = DN(1,0)*clhs217 + DN(1,1)*clhs218 + DN(1,2)*clhs219 + N[1]*clhs40;
const double clhs221 = -N[2];
const double clhs222 = clhs221 + clhs43*clhs85;
const double clhs223 = DN(2,0)*clhs217 + DN(2,1)*clhs218 + DN(2,2)*clhs219 + N[2]*clhs40;
const double clhs224 = -N[3];
const double clhs225 = clhs113*clhs43 + clhs224;
const double clhs226 = DN(3,0)*clhs217 + DN(3,1)*clhs218 + DN(3,2)*clhs219 + N[3]*clhs40;
const double clhs227 = clhs18*clhs54;
const double clhs228 = N[1]*clhs21;
const double clhs229 = clhs11*clhs51 + clhs16*clhs227 + clhs16*clhs228;
const double clhs230 = DN(1,0)*clhs10;
const double clhs231 = clhs43*clhs54;
const double clhs232 = pow(DN(1,0), 2);
const double clhs233 = pow(N[1], 2);
const double clhs234 = clhs14*clhs233 + clhs227*clhs56 + clhs228*clhs56 + clhs51*clhs54;
const double clhs235 = DN(1,1)*clhs230;
const double clhs236 = DN(1,2)*clhs230;
const double clhs237 = clhs10*clhs72;
const double clhs238 = clhs215 + clhs231 + clhs237 + clhs42*clhs51;
const double clhs239 = clhs55*clhs80;
const double clhs240 = DN(2,0)*clhs230 + clhs239;
const double clhs241 = clhs227*clhs85 + clhs228*clhs85 + clhs51*clhs83;
const double clhs242 = DN(2,1)*clhs230;
const double clhs243 = DN(2,2)*clhs230;
const double clhs244 = DN(1,0)*N[2];
const double clhs245 = DN(2,0)*N[1];
const double clhs246 = clhs108*clhs55;
const double clhs247 = DN(3,0)*clhs230 + clhs246;
const double clhs248 = clhs111*clhs51 + clhs113*clhs227 + clhs113*clhs228;
const double clhs249 = DN(3,1)*clhs230;
const double clhs250 = DN(3,2)*clhs230;
const double clhs251 = DN(1,0)*N[3];
const double clhs252 = DN(3,0)*N[1];
const double clhs253 = clhs229 + clhs52;
const double clhs254 = DN(1,1)*clhs10;
const double clhs255 = pow(DN(1,1), 2);
const double clhs256 = DN(1,2)*clhs254;
const double clhs257 = DN(2,0)*clhs254;
const double clhs258 = DN(2,1)*clhs254;
const double clhs259 = clhs239 + clhs241;
const double clhs260 = DN(2,2)*clhs254;
const double clhs261 = DN(1,1)*N[2];
const double clhs262 = DN(2,1)*N[1];
const double clhs263 = DN(3,0)*clhs254;
const double clhs264 = DN(3,1)*clhs254;
const double clhs265 = clhs246 + clhs248;
const double clhs266 = DN(3,2)*clhs254;
const double clhs267 = DN(1,1)*N[3];
const double clhs268 = DN(3,1)*N[1];
const double clhs269 = DN(1,2)*clhs10;
const double clhs270 = pow(DN(1,2), 2);
const double clhs271 = DN(2,0)*clhs269;
const double clhs272 = DN(2,1)*clhs269;
const double clhs273 = DN(2,2)*clhs269;
const double clhs274 = DN(1,2)*N[2];
const double clhs275 = DN(2,2)*N[1];
const double clhs276 = DN(3,0)*clhs269;
const double clhs277 = DN(3,1)*clhs269;
const double clhs278 = DN(3,2)*clhs269;
const double clhs279 = DN(1,2)*N[3];
const double clhs280 = DN(3,2)*N[1];
const double clhs281 = DN(1,0)*clhs17;
const double clhs282 = DN(1,1)*clhs17;
const double clhs283 = DN(1,2)*clhs17;
const double clhs284 = DN(2,0)*clhs281 + DN(2,1)*clhs282 + DN(2,2)*clhs283 + N[2]*clhs72;
const double clhs285 = DN(3,0)*clhs281 + DN(3,1)*clhs282 + DN(3,2)*clhs283 + N[3]*clhs72;
const double clhs286 = clhs18*clhs83;
const double clhs287 = N[2]*clhs21;
const double clhs288 = clhs11*clhs80 + clhs16*clhs286 + clhs16*clhs287;
const double clhs289 = DN(2,0)*clhs10;
const double clhs290 = clhs43*clhs83;
const double clhs291 = clhs286*clhs56 + clhs287*clhs56 + clhs54*clhs80;
const double clhs292 = pow(DN(2,0), 2);
const double clhs293 = pow(N[2], 2);
const double clhs294 = clhs14*clhs293 + clhs286*clhs85 + clhs287*clhs85 + clhs80*clhs83;
const double clhs295 = DN(2,1)*clhs289;
const double clhs296 = DN(2,2)*clhs289;
const double clhs297 = clhs10*clhs101;
const double clhs298 = clhs221 + clhs290 + clhs297 + clhs42*clhs80;
const double clhs299 = clhs108*clhs84;
const double clhs300 = DN(3,0)*clhs289 + clhs299;
const double clhs301 = clhs111*clhs80 + clhs113*clhs286 + clhs113*clhs287;
const double clhs302 = DN(3,1)*clhs289;
const double clhs303 = DN(3,2)*clhs289;
const double clhs304 = DN(2,0)*N[3];
const double clhs305 = DN(3,0)*N[2];
const double clhs306 = clhs288 + clhs81;
const double clhs307 = DN(2,1)*clhs10;
const double clhs308 = clhs239 + clhs291;
const double clhs309 = pow(DN(2,1), 2);
const double clhs310 = DN(2,2)*clhs307;
const double clhs311 = DN(3,0)*clhs307;
const double clhs312 = DN(3,1)*clhs307;
const double clhs313 = clhs299 + clhs301;
const double clhs314 = DN(3,2)*clhs307;
const double clhs315 = DN(2,1)*N[3];
const double clhs316 = DN(3,1)*N[2];
const double clhs317 = DN(2,2)*clhs10;
const double clhs318 = pow(DN(2,2), 2);
const double clhs319 = DN(3,0)*clhs317;
const double clhs320 = DN(3,1)*clhs317;
const double clhs321 = DN(3,2)*clhs317;
const double clhs322 = DN(2,2)*N[3];
const double clhs323 = DN(3,2)*N[2];
const double clhs324 = DN(2,0)*DN(3,0)*clhs17 + DN(2,1)*DN(3,1)*clhs17 + DN(2,2)*DN(3,2)*clhs17 + N[3]*clhs101;
const double clhs325 = clhs111*clhs18;
const double clhs326 = N[3]*clhs21;
const double clhs327 = clhs108*clhs11 + clhs16*clhs325 + clhs16*clhs326;
const double clhs328 = DN(3,0)*clhs10;
const double clhs329 = clhs111*clhs43;
const double clhs330 = clhs108*clhs54 + clhs325*clhs56 + clhs326*clhs56;
const double clhs331 = clhs108*clhs83 + clhs325*clhs85 + clhs326*clhs85;
const double clhs332 = pow(DN(3,0), 2);
const double clhs333 = pow(N[3], 2);
const double clhs334 = clhs108*clhs111 + clhs113*clhs325 + clhs113*clhs326 + clhs14*clhs333;
const double clhs335 = DN(3,1)*clhs328;
const double clhs336 = DN(3,2)*clhs328;
const double clhs337 = clhs10*clhs129 + clhs108*clhs42 + clhs224 + clhs329;
const double clhs338 = clhs109 + clhs327;
const double clhs339 = DN(3,1)*clhs10;
const double clhs340 = clhs246 + clhs330;
const double clhs341 = clhs299 + clhs331;
const double clhs342 = pow(DN(3,1), 2);
const double clhs343 = DN(3,2)*clhs339;
const double clhs344 = pow(DN(3,2), 2);
lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs10*clhs5 + clhs23;
lhs(0,1)=DN(0,0)*clhs24 + DN(0,1)*clhs26 + DN(0,2)*clhs29 + clhs31;
lhs(0,2)=DN(0,0)*clhs32 + DN(0,1)*clhs34 + DN(0,2)*clhs36 + clhs37;
lhs(0,3)=DN(0,0)*clhs45;
lhs(0,4)=DN(0,0)*clhs46 + DN(0,1)*clhs48 + DN(0,2)*clhs50 + clhs53 + clhs57;
lhs(0,5)=DN(0,0)*clhs58 + DN(0,1)*clhs60 + DN(0,2)*clhs63 + clhs64;
lhs(0,6)=DN(0,0)*clhs65 + DN(0,1)*clhs67 + DN(0,2)*clhs69 + clhs70;
lhs(0,7)=DN(1,0)*clhs44 + clhs30*clhs72 - clhs71 + clhs73*clhs74;
lhs(0,8)=DN(0,0)*clhs75 + DN(0,1)*clhs77 + DN(0,2)*clhs79 + clhs82 + clhs86;
lhs(0,9)=DN(0,0)*clhs87 + DN(0,1)*clhs89 + DN(0,2)*clhs92 + clhs93;
lhs(0,10)=DN(0,0)*clhs94 + DN(0,1)*clhs96 + DN(0,2)*clhs98 + clhs99;
lhs(0,11)=DN(2,0)*clhs44 - clhs100 + clhs101*clhs30 + clhs102*clhs74;
lhs(0,12)=DN(0,0)*clhs103 + DN(0,1)*clhs105 + DN(0,2)*clhs107 + clhs110 + clhs114;
lhs(0,13)=DN(0,0)*clhs115 + DN(0,1)*clhs117 + DN(0,2)*clhs120 + clhs121;
lhs(0,14)=DN(0,0)*clhs122 + DN(0,1)*clhs124 + DN(0,2)*clhs126 + clhs127;
lhs(0,15)=DN(3,0)*clhs44 - clhs128 + clhs129*clhs30 + clhs130*clhs74;
lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs131 + DN(0,2)*clhs132 + clhs31;
lhs(1,1)=DN(0,0)*clhs26 + DN(0,1)*clhs133 + DN(0,2)*clhs135 + clhs10*clhs136 + clhs23;
lhs(1,2)=DN(0,0)*clhs34 + DN(0,1)*clhs137 + DN(0,2)*clhs139 + clhs141;
lhs(1,3)=DN(0,1)*clhs45;
lhs(1,4)=DN(0,0)*clhs48 + DN(0,1)*clhs142 + DN(0,2)*clhs143 + clhs144;
lhs(1,5)=DN(0,0)*clhs60 + DN(0,1)*clhs145 + DN(0,2)*clhs147 + clhs148 + clhs149;
lhs(1,6)=DN(0,0)*clhs67 + DN(0,1)*clhs150 + DN(0,2)*clhs152 + clhs153;
lhs(1,7)=DN(1,1)*clhs44 + clhs140*clhs72 - clhs154 + clhs155*clhs74;
lhs(1,8)=DN(0,0)*clhs77 + DN(0,1)*clhs156 + DN(0,2)*clhs157 + clhs158;
lhs(1,9)=DN(0,0)*clhs89 + DN(0,1)*clhs159 + DN(0,2)*clhs161 + clhs162 + clhs163;
lhs(1,10)=DN(0,0)*clhs96 + DN(0,1)*clhs164 + DN(0,2)*clhs166 + clhs167;
lhs(1,11)=DN(2,1)*clhs44 + clhs101*clhs140 - clhs168 + clhs169*clhs74;
lhs(1,12)=DN(0,0)*clhs105 + DN(0,1)*clhs170 + DN(0,2)*clhs171 + clhs172;
lhs(1,13)=DN(0,0)*clhs117 + DN(0,1)*clhs173 + DN(0,2)*clhs175 + clhs176 + clhs177;
lhs(1,14)=DN(0,0)*clhs124 + DN(0,1)*clhs178 + DN(0,2)*clhs180 + clhs181;
lhs(1,15)=DN(3,1)*clhs44 + clhs129*clhs140 - clhs182 + clhs183*clhs74;
lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs132 + DN(0,2)*clhs184 + clhs37;
lhs(2,1)=DN(0,0)*clhs29 + DN(0,1)*clhs135 + DN(0,2)*clhs185 + clhs141;
lhs(2,2)=DN(0,0)*clhs36 + DN(0,1)*clhs139 + DN(0,2)*clhs186 + clhs10*clhs187 + clhs23;
lhs(2,3)=DN(0,2)*clhs45;
lhs(2,4)=DN(0,0)*clhs50 + DN(0,1)*clhs143 + DN(0,2)*clhs188 + clhs190;
lhs(2,5)=DN(0,0)*clhs63 + DN(0,1)*clhs147 + DN(0,2)*clhs191 + clhs192;
lhs(2,6)=DN(0,0)*clhs69 + DN(0,1)*clhs152 + DN(0,2)*clhs193 + clhs149 + clhs194;
lhs(2,7)=DN(1,2)*clhs44 + clhs189*clhs72 - clhs195 + clhs196*clhs74;
lhs(2,8)=DN(0,0)*clhs79 + DN(0,1)*clhs157 + DN(0,2)*clhs197 + clhs198;
lhs(2,9)=DN(0,0)*clhs92 + DN(0,1)*clhs161 + DN(0,2)*clhs199 + clhs200;
lhs(2,10)=DN(0,0)*clhs98 + DN(0,1)*clhs166 + DN(0,2)*clhs201 + clhs163 + clhs202;
lhs(2,11)=DN(2,2)*clhs44 + clhs101*clhs189 - clhs203 + clhs204*clhs74;
lhs(2,12)=DN(0,0)*clhs107 + DN(0,1)*clhs171 + DN(0,2)*clhs205 + clhs206;
lhs(2,13)=DN(0,0)*clhs120 + DN(0,1)*clhs175 + DN(0,2)*clhs207 + clhs208;
lhs(2,14)=DN(0,0)*clhs126 + DN(0,1)*clhs180 + DN(0,2)*clhs209 + clhs177 + clhs210;
lhs(2,15)=DN(3,2)*clhs44 + clhs129*clhs189 - clhs211 + clhs212*clhs74;
lhs(3,0)=DN(0,0)*clhs213;
lhs(3,1)=DN(0,1)*clhs213;
lhs(3,2)=DN(0,2)*clhs213;
lhs(3,3)=clhs13*clhs214 + clhs136*clhs17 + clhs17*clhs187 + clhs17*clhs5;
lhs(3,4)=DN(0,0)*clhs216;
lhs(3,5)=DN(0,1)*clhs216;
lhs(3,6)=DN(0,2)*clhs216;
lhs(3,7)=clhs220;
lhs(3,8)=DN(0,0)*clhs222;
lhs(3,9)=DN(0,1)*clhs222;
lhs(3,10)=DN(0,2)*clhs222;
lhs(3,11)=clhs223;
lhs(3,12)=DN(0,0)*clhs225;
lhs(3,13)=DN(0,1)*clhs225;
lhs(3,14)=DN(0,2)*clhs225;
lhs(3,15)=clhs226;
lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs229 + clhs53;
lhs(4,1)=DN(1,0)*clhs24 + DN(1,1)*clhs26 + DN(1,2)*clhs29 + clhs144;
lhs(4,2)=DN(1,0)*clhs32 + DN(1,1)*clhs34 + DN(1,2)*clhs36 + clhs190;
lhs(4,3)=DN(0,0)*clhs231 + clhs230*clhs40 + clhs71*clhs74 - clhs73;
lhs(4,4)=DN(1,0)*clhs46 + DN(1,1)*clhs48 + DN(1,2)*clhs50 + clhs10*clhs232 + clhs234;
lhs(4,5)=DN(1,0)*clhs58 + DN(1,1)*clhs60 + DN(1,2)*clhs63 + clhs235;
lhs(4,6)=DN(1,0)*clhs65 + DN(1,1)*clhs67 + DN(1,2)*clhs69 + clhs236;
lhs(4,7)=DN(1,0)*clhs238;
lhs(4,8)=DN(1,0)*clhs75 + DN(1,1)*clhs77 + DN(1,2)*clhs79 + clhs240 + clhs241;
lhs(4,9)=DN(1,0)*clhs87 + DN(1,1)*clhs89 + DN(1,2)*clhs92 + clhs242;
lhs(4,10)=DN(1,0)*clhs94 + DN(1,1)*clhs96 + DN(1,2)*clhs98 + clhs243;
lhs(4,11)=DN(2,0)*clhs231 + clhs101*clhs230 - clhs244 + clhs245*clhs74;
lhs(4,12)=DN(1,0)*clhs103 + DN(1,1)*clhs105 + DN(1,2)*clhs107 + clhs247 + clhs248;
lhs(4,13)=DN(1,0)*clhs115 + DN(1,1)*clhs117 + DN(1,2)*clhs120 + clhs249;
lhs(4,14)=DN(1,0)*clhs122 + DN(1,1)*clhs124 + DN(1,2)*clhs126 + clhs250;
lhs(4,15)=DN(3,0)*clhs231 + clhs129*clhs230 - clhs251 + clhs252*clhs74;
lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs131 + DN(1,2)*clhs132 + clhs64;
lhs(5,1)=DN(1,0)*clhs26 + DN(1,1)*clhs133 + DN(1,2)*clhs135 + clhs148 + clhs253;
lhs(5,2)=DN(1,0)*clhs34 + DN(1,1)*clhs137 + DN(1,2)*clhs139 + clhs192;
lhs(5,3)=DN(0,1)*clhs231 + clhs154*clhs74 - clhs155 + clhs254*clhs40;
lhs(5,4)=DN(1,0)*clhs48 + DN(1,1)*clhs142 + DN(1,2)*clhs143 + clhs235;
lhs(5,5)=DN(1,0)*clhs60 + DN(1,1)*clhs145 + DN(1,2)*clhs147 + clhs10*clhs255 + clhs234;
lhs(5,6)=DN(1,0)*clhs67 + DN(1,1)*clhs150 + DN(1,2)*clhs152 + clhs256;
lhs(5,7)=DN(1,1)*clhs238;
lhs(5,8)=DN(1,0)*clhs77 + DN(1,1)*clhs156 + DN(1,2)*clhs157 + clhs257;
lhs(5,9)=DN(1,0)*clhs89 + DN(1,1)*clhs159 + DN(1,2)*clhs161 + clhs258 + clhs259;
lhs(5,10)=DN(1,0)*clhs96 + DN(1,1)*clhs164 + DN(1,2)*clhs166 + clhs260;
lhs(5,11)=DN(2,1)*clhs231 + clhs101*clhs254 - clhs261 + clhs262*clhs74;
lhs(5,12)=DN(1,0)*clhs105 + DN(1,1)*clhs170 + DN(1,2)*clhs171 + clhs263;
lhs(5,13)=DN(1,0)*clhs117 + DN(1,1)*clhs173 + DN(1,2)*clhs175 + clhs264 + clhs265;
lhs(5,14)=DN(1,0)*clhs124 + DN(1,1)*clhs178 + DN(1,2)*clhs180 + clhs266;
lhs(5,15)=DN(3,1)*clhs231 + clhs129*clhs254 - clhs267 + clhs268*clhs74;
lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs132 + DN(1,2)*clhs184 + clhs70;
lhs(6,1)=DN(1,0)*clhs29 + DN(1,1)*clhs135 + DN(1,2)*clhs185 + clhs153;
lhs(6,2)=DN(1,0)*clhs36 + DN(1,1)*clhs139 + DN(1,2)*clhs186 + clhs194 + clhs253;
lhs(6,3)=DN(0,2)*clhs231 + clhs195*clhs74 - clhs196 + clhs269*clhs40;
lhs(6,4)=DN(1,0)*clhs50 + DN(1,1)*clhs143 + DN(1,2)*clhs188 + clhs236;
lhs(6,5)=DN(1,0)*clhs63 + DN(1,1)*clhs147 + DN(1,2)*clhs191 + clhs256;
lhs(6,6)=DN(1,0)*clhs69 + DN(1,1)*clhs152 + DN(1,2)*clhs193 + clhs10*clhs270 + clhs234;
lhs(6,7)=DN(1,2)*clhs238;
lhs(6,8)=DN(1,0)*clhs79 + DN(1,1)*clhs157 + DN(1,2)*clhs197 + clhs271;
lhs(6,9)=DN(1,0)*clhs92 + DN(1,1)*clhs161 + DN(1,2)*clhs199 + clhs272;
lhs(6,10)=DN(1,0)*clhs98 + DN(1,1)*clhs166 + DN(1,2)*clhs201 + clhs259 + clhs273;
lhs(6,11)=DN(2,2)*clhs231 + clhs101*clhs269 - clhs274 + clhs275*clhs74;
lhs(6,12)=DN(1,0)*clhs107 + DN(1,1)*clhs171 + DN(1,2)*clhs205 + clhs276;
lhs(6,13)=DN(1,0)*clhs120 + DN(1,1)*clhs175 + DN(1,2)*clhs207 + clhs277;
lhs(6,14)=DN(1,0)*clhs126 + DN(1,1)*clhs180 + DN(1,2)*clhs209 + clhs265 + clhs278;
lhs(6,15)=DN(3,2)*clhs231 + clhs129*clhs269 - clhs279 + clhs280*clhs74;
lhs(7,0)=DN(1,0)*clhs213;
lhs(7,1)=DN(1,1)*clhs213;
lhs(7,2)=DN(1,2)*clhs213;
lhs(7,3)=clhs220;
lhs(7,4)=DN(1,0)*clhs216;
lhs(7,5)=DN(1,1)*clhs216;
lhs(7,6)=DN(1,2)*clhs216;
lhs(7,7)=clhs17*clhs232 + clhs17*clhs255 + clhs17*clhs270 + clhs214*clhs233;
lhs(7,8)=DN(1,0)*clhs222;
lhs(7,9)=DN(1,1)*clhs222;
lhs(7,10)=DN(1,2)*clhs222;
lhs(7,11)=clhs284;
lhs(7,12)=DN(1,0)*clhs225;
lhs(7,13)=DN(1,1)*clhs225;
lhs(7,14)=DN(1,2)*clhs225;
lhs(7,15)=clhs285;
lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs288 + clhs82;
lhs(8,1)=DN(2,0)*clhs24 + DN(2,1)*clhs26 + DN(2,2)*clhs29 + clhs158;
lhs(8,2)=DN(2,0)*clhs32 + DN(2,1)*clhs34 + DN(2,2)*clhs36 + clhs198;
lhs(8,3)=DN(0,0)*clhs290 + clhs100*clhs74 - clhs102 + clhs289*clhs40;
lhs(8,4)=DN(2,0)*clhs46 + DN(2,1)*clhs48 + DN(2,2)*clhs50 + clhs240 + clhs291;
lhs(8,5)=DN(2,0)*clhs58 + DN(2,1)*clhs60 + DN(2,2)*clhs63 + clhs257;
lhs(8,6)=DN(2,0)*clhs65 + DN(2,1)*clhs67 + DN(2,2)*clhs69 + clhs271;
lhs(8,7)=DN(1,0)*clhs290 + clhs244*clhs74 - clhs245 + clhs289*clhs72;
lhs(8,8)=DN(2,0)*clhs75 + DN(2,1)*clhs77 + DN(2,2)*clhs79 + clhs10*clhs292 + clhs294;
lhs(8,9)=DN(2,0)*clhs87 + DN(2,1)*clhs89 + DN(2,2)*clhs92 + clhs295;
lhs(8,10)=DN(2,0)*clhs94 + DN(2,1)*clhs96 + DN(2,2)*clhs98 + clhs296;
lhs(8,11)=DN(2,0)*clhs298;
lhs(8,12)=DN(2,0)*clhs103 + DN(2,1)*clhs105 + DN(2,2)*clhs107 + clhs300 + clhs301;
lhs(8,13)=DN(2,0)*clhs115 + DN(2,1)*clhs117 + DN(2,2)*clhs120 + clhs302;
lhs(8,14)=DN(2,0)*clhs122 + DN(2,1)*clhs124 + DN(2,2)*clhs126 + clhs303;
lhs(8,15)=DN(3,0)*clhs290 + clhs129*clhs289 - clhs304 + clhs305*clhs74;
lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs131 + DN(2,2)*clhs132 + clhs93;
lhs(9,1)=DN(2,0)*clhs26 + DN(2,1)*clhs133 + DN(2,2)*clhs135 + clhs162 + clhs306;
lhs(9,2)=DN(2,0)*clhs34 + DN(2,1)*clhs137 + DN(2,2)*clhs139 + clhs200;
lhs(9,3)=DN(0,1)*clhs290 + clhs168*clhs74 - clhs169 + clhs307*clhs40;
lhs(9,4)=DN(2,0)*clhs48 + DN(2,1)*clhs142 + DN(2,2)*clhs143 + clhs242;
lhs(9,5)=DN(2,0)*clhs60 + DN(2,1)*clhs145 + DN(2,2)*clhs147 + clhs258 + clhs308;
lhs(9,6)=DN(2,0)*clhs67 + DN(2,1)*clhs150 + DN(2,2)*clhs152 + clhs272;
lhs(9,7)=DN(1,1)*clhs290 + clhs261*clhs74 - clhs262 + clhs307*clhs72;
lhs(9,8)=DN(2,0)*clhs77 + DN(2,1)*clhs156 + DN(2,2)*clhs157 + clhs295;
lhs(9,9)=DN(2,0)*clhs89 + DN(2,1)*clhs159 + DN(2,2)*clhs161 + clhs10*clhs309 + clhs294;
lhs(9,10)=DN(2,0)*clhs96 + DN(2,1)*clhs164 + DN(2,2)*clhs166 + clhs310;
lhs(9,11)=DN(2,1)*clhs298;
lhs(9,12)=DN(2,0)*clhs105 + DN(2,1)*clhs170 + DN(2,2)*clhs171 + clhs311;
lhs(9,13)=DN(2,0)*clhs117 + DN(2,1)*clhs173 + DN(2,2)*clhs175 + clhs312 + clhs313;
lhs(9,14)=DN(2,0)*clhs124 + DN(2,1)*clhs178 + DN(2,2)*clhs180 + clhs314;
lhs(9,15)=DN(3,1)*clhs290 + clhs129*clhs307 - clhs315 + clhs316*clhs74;
lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs132 + DN(2,2)*clhs184 + clhs99;
lhs(10,1)=DN(2,0)*clhs29 + DN(2,1)*clhs135 + DN(2,2)*clhs185 + clhs167;
lhs(10,2)=DN(2,0)*clhs36 + DN(2,1)*clhs139 + DN(2,2)*clhs186 + clhs202 + clhs306;
lhs(10,3)=DN(0,2)*clhs290 + clhs203*clhs74 - clhs204 + clhs317*clhs40;
lhs(10,4)=DN(2,0)*clhs50 + DN(2,1)*clhs143 + DN(2,2)*clhs188 + clhs243;
lhs(10,5)=DN(2,0)*clhs63 + DN(2,1)*clhs147 + DN(2,2)*clhs191 + clhs260;
lhs(10,6)=DN(2,0)*clhs69 + DN(2,1)*clhs152 + DN(2,2)*clhs193 + clhs273 + clhs308;
lhs(10,7)=DN(1,2)*clhs290 + clhs274*clhs74 - clhs275 + clhs317*clhs72;
lhs(10,8)=DN(2,0)*clhs79 + DN(2,1)*clhs157 + DN(2,2)*clhs197 + clhs296;
lhs(10,9)=DN(2,0)*clhs92 + DN(2,1)*clhs161 + DN(2,2)*clhs199 + clhs310;
lhs(10,10)=DN(2,0)*clhs98 + DN(2,1)*clhs166 + DN(2,2)*clhs201 + clhs10*clhs318 + clhs294;
lhs(10,11)=DN(2,2)*clhs298;
lhs(10,12)=DN(2,0)*clhs107 + DN(2,1)*clhs171 + DN(2,2)*clhs205 + clhs319;
lhs(10,13)=DN(2,0)*clhs120 + DN(2,1)*clhs175 + DN(2,2)*clhs207 + clhs320;
lhs(10,14)=DN(2,0)*clhs126 + DN(2,1)*clhs180 + DN(2,2)*clhs209 + clhs313 + clhs321;
lhs(10,15)=DN(3,2)*clhs290 + clhs129*clhs317 - clhs322 + clhs323*clhs74;
lhs(11,0)=DN(2,0)*clhs213;
lhs(11,1)=DN(2,1)*clhs213;
lhs(11,2)=DN(2,2)*clhs213;
lhs(11,3)=clhs223;
lhs(11,4)=DN(2,0)*clhs216;
lhs(11,5)=DN(2,1)*clhs216;
lhs(11,6)=DN(2,2)*clhs216;
lhs(11,7)=clhs284;
lhs(11,8)=DN(2,0)*clhs222;
lhs(11,9)=DN(2,1)*clhs222;
lhs(11,10)=DN(2,2)*clhs222;
lhs(11,11)=clhs17*clhs292 + clhs17*clhs309 + clhs17*clhs318 + clhs214*clhs293;
lhs(11,12)=DN(2,0)*clhs225;
lhs(11,13)=DN(2,1)*clhs225;
lhs(11,14)=DN(2,2)*clhs225;
lhs(11,15)=clhs324;
lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs110 + clhs327;
lhs(12,1)=DN(3,0)*clhs24 + DN(3,1)*clhs26 + DN(3,2)*clhs29 + clhs172;
lhs(12,2)=DN(3,0)*clhs32 + DN(3,1)*clhs34 + DN(3,2)*clhs36 + clhs206;
lhs(12,3)=DN(0,0)*clhs329 + clhs128*clhs74 - clhs130 + clhs328*clhs40;
lhs(12,4)=DN(3,0)*clhs46 + DN(3,1)*clhs48 + DN(3,2)*clhs50 + clhs247 + clhs330;
lhs(12,5)=DN(3,0)*clhs58 + DN(3,1)*clhs60 + DN(3,2)*clhs63 + clhs263;
lhs(12,6)=DN(3,0)*clhs65 + DN(3,1)*clhs67 + DN(3,2)*clhs69 + clhs276;
lhs(12,7)=DN(1,0)*clhs329 + clhs251*clhs74 - clhs252 + clhs328*clhs72;
lhs(12,8)=DN(3,0)*clhs75 + DN(3,1)*clhs77 + DN(3,2)*clhs79 + clhs300 + clhs331;
lhs(12,9)=DN(3,0)*clhs87 + DN(3,1)*clhs89 + DN(3,2)*clhs92 + clhs311;
lhs(12,10)=DN(3,0)*clhs94 + DN(3,1)*clhs96 + DN(3,2)*clhs98 + clhs319;
lhs(12,11)=DN(2,0)*clhs329 + clhs101*clhs328 + clhs304*clhs74 - clhs305;
lhs(12,12)=DN(3,0)*clhs103 + DN(3,1)*clhs105 + DN(3,2)*clhs107 + clhs10*clhs332 + clhs334;
lhs(12,13)=DN(3,0)*clhs115 + DN(3,1)*clhs117 + DN(3,2)*clhs120 + clhs335;
lhs(12,14)=DN(3,0)*clhs122 + DN(3,1)*clhs124 + DN(3,2)*clhs126 + clhs336;
lhs(12,15)=DN(3,0)*clhs337;
lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs131 + DN(3,2)*clhs132 + clhs121;
lhs(13,1)=DN(3,0)*clhs26 + DN(3,1)*clhs133 + DN(3,2)*clhs135 + clhs176 + clhs338;
lhs(13,2)=DN(3,0)*clhs34 + DN(3,1)*clhs137 + DN(3,2)*clhs139 + clhs208;
lhs(13,3)=DN(0,1)*clhs329 + clhs182*clhs74 - clhs183 + clhs339*clhs40;
lhs(13,4)=DN(3,0)*clhs48 + DN(3,1)*clhs142 + DN(3,2)*clhs143 + clhs249;
lhs(13,5)=DN(3,0)*clhs60 + DN(3,1)*clhs145 + DN(3,2)*clhs147 + clhs264 + clhs340;
lhs(13,6)=DN(3,0)*clhs67 + DN(3,1)*clhs150 + DN(3,2)*clhs152 + clhs277;
lhs(13,7)=DN(1,1)*clhs329 + clhs267*clhs74 - clhs268 + clhs339*clhs72;
lhs(13,8)=DN(3,0)*clhs77 + DN(3,1)*clhs156 + DN(3,2)*clhs157 + clhs302;
lhs(13,9)=DN(3,0)*clhs89 + DN(3,1)*clhs159 + DN(3,2)*clhs161 + clhs312 + clhs341;
lhs(13,10)=DN(3,0)*clhs96 + DN(3,1)*clhs164 + DN(3,2)*clhs166 + clhs320;
lhs(13,11)=DN(2,1)*clhs329 + clhs101*clhs339 + clhs315*clhs74 - clhs316;
lhs(13,12)=DN(3,0)*clhs105 + DN(3,1)*clhs170 + DN(3,2)*clhs171 + clhs335;
lhs(13,13)=DN(3,0)*clhs117 + DN(3,1)*clhs173 + DN(3,2)*clhs175 + clhs10*clhs342 + clhs334;
lhs(13,14)=DN(3,0)*clhs124 + DN(3,1)*clhs178 + DN(3,2)*clhs180 + clhs343;
lhs(13,15)=DN(3,1)*clhs337;
lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs132 + DN(3,2)*clhs184 + clhs127;
lhs(14,1)=DN(3,0)*clhs29 + DN(3,1)*clhs135 + DN(3,2)*clhs185 + clhs181;
lhs(14,2)=DN(3,0)*clhs36 + DN(3,1)*clhs139 + DN(3,2)*clhs186 + clhs210 + clhs338;
lhs(14,3)=DN(0,2)*clhs329 + DN(3,2)*clhs41 + clhs211*clhs74 - clhs212;
lhs(14,4)=DN(3,0)*clhs50 + DN(3,1)*clhs143 + DN(3,2)*clhs188 + clhs250;
lhs(14,5)=DN(3,0)*clhs63 + DN(3,1)*clhs147 + DN(3,2)*clhs191 + clhs266;
lhs(14,6)=DN(3,0)*clhs69 + DN(3,1)*clhs152 + DN(3,2)*clhs193 + clhs278 + clhs340;
lhs(14,7)=DN(1,2)*clhs329 + DN(3,2)*clhs237 + clhs279*clhs74 - clhs280;
lhs(14,8)=DN(3,0)*clhs79 + DN(3,1)*clhs157 + DN(3,2)*clhs197 + clhs303;
lhs(14,9)=DN(3,0)*clhs92 + DN(3,1)*clhs161 + DN(3,2)*clhs199 + clhs314;
lhs(14,10)=DN(3,0)*clhs98 + DN(3,1)*clhs166 + DN(3,2)*clhs201 + clhs321 + clhs341;
lhs(14,11)=DN(2,2)*clhs329 + DN(3,2)*clhs297 + clhs322*clhs74 - clhs323;
lhs(14,12)=DN(3,0)*clhs107 + DN(3,1)*clhs171 + DN(3,2)*clhs205 + clhs336;
lhs(14,13)=DN(3,0)*clhs120 + DN(3,1)*clhs175 + DN(3,2)*clhs207 + clhs343;
lhs(14,14)=DN(3,0)*clhs126 + DN(3,1)*clhs180 + DN(3,2)*clhs209 + clhs10*clhs344 + clhs334;
lhs(14,15)=DN(3,2)*clhs337;
lhs(15,0)=DN(3,0)*clhs213;
lhs(15,1)=DN(3,1)*clhs213;
lhs(15,2)=DN(3,2)*clhs213;
lhs(15,3)=clhs226;
lhs(15,4)=DN(3,0)*clhs216;
lhs(15,5)=DN(3,1)*clhs216;
lhs(15,6)=DN(3,2)*clhs216;
lhs(15,7)=clhs285;
lhs(15,8)=DN(3,0)*clhs222;
lhs(15,9)=DN(3,1)*clhs222;
lhs(15,10)=DN(3,2)*clhs222;
lhs(15,11)=clhs324;
lhs(15,12)=DN(3,0)*clhs225;
lhs(15,13)=DN(3,1)*clhs225;
lhs(15,14)=DN(3,2)*clhs225;
lhs(15,15)=clhs17*clhs332 + clhs17*clhs342 + clhs17*clhs344 + clhs214*clhs333;


}


template<>
void EmbeddedAusasNavierStokes<2>::ComputeGaussPointLHSContribution(
    BoundedMatrix<double,9,9>& lhs,
    const EmbeddedAusasElementDataStruct& data)
{
    constexpr unsigned int dim = 2;
    constexpr unsigned int nnodes = 3;

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

    const double clhs0 = C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double clhs1 = C(0,2)*DN(0,0);
const double clhs2 = C(2,2)*DN(0,1) + clhs1;
const double clhs3 = pow(DN(0,0), 2);
const double clhs4 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double clhs5 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double clhs6 = rho*stab_c2*sqrt(pow(clhs4, 2) + pow(clhs5, 2));
const double clhs7 = clhs6*h/stab_c1 + mu;
const double clhs8 = DN(0,0)*clhs4 + DN(0,1)*clhs5;
const double clhs9 = N[0]*rho;
const double clhs10 = pow(N[0], 2);
const double clhs11 = bdf0*rho;
const double clhs12 = N[0]*bdf0;
const double clhs13 = clhs12 + clhs8;
const double clhs14 = 1.0/(clhs6/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double clhs15 = clhs14*pow(rho, 2);
const double clhs16 = clhs15*clhs8;
const double clhs17 = DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double clhs18 = clhs15*clhs17;
const double clhs19 = N[0]*clhs18;
const double clhs20 = clhs10*clhs11 + clhs13*clhs16 + clhs13*clhs19 + clhs8*clhs9;
const double clhs21 = C(0,1)*DN(0,1) + clhs1;
const double clhs22 = C(1,2)*DN(0,1);
const double clhs23 = C(2,2)*DN(0,0) + clhs22;
const double clhs24 = DN(0,0)*clhs7;
const double clhs25 = DN(0,1)*clhs24;
const double clhs26 = -N[0];
const double clhs27 = 1/(pow(c, 2)*rho);
const double clhs28 = clhs12*clhs27;
const double clhs29 = clhs28*clhs7;
const double clhs30 = clhs14*clhs17;
const double clhs31 = clhs14*rho;
const double clhs32 = clhs31*clhs8;
const double clhs33 = clhs26 + clhs29 + clhs30*clhs9 + clhs32;
const double clhs34 = C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs35 = C(0,2)*DN(1,0);
const double clhs36 = C(2,2)*DN(1,1) + clhs35;
const double clhs37 = N[1]*rho;
const double clhs38 = clhs12*clhs37;
const double clhs39 = DN(1,0)*clhs24 + clhs38;
const double clhs40 = DN(1,0)*clhs4 + DN(1,1)*clhs5;
const double clhs41 = N[1]*bdf0;
const double clhs42 = clhs40 + clhs41;
const double clhs43 = clhs16*clhs42 + clhs19*clhs42 + clhs40*clhs9;
const double clhs44 = C(0,1)*DN(1,1) + clhs35;
const double clhs45 = C(1,2)*DN(1,1);
const double clhs46 = C(2,2)*DN(1,0) + clhs45;
const double clhs47 = DN(1,1)*clhs24;
const double clhs48 = DN(0,0)*N[1];
const double clhs49 = clhs27*clhs41;
const double clhs50 = DN(1,0)*N[0];
const double clhs51 = clhs17*clhs31;
const double clhs52 = C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs53 = C(0,2)*DN(2,0);
const double clhs54 = C(2,2)*DN(2,1) + clhs53;
const double clhs55 = N[2]*rho;
const double clhs56 = clhs12*clhs55;
const double clhs57 = DN(2,0)*clhs24 + clhs56;
const double clhs58 = DN(2,0)*clhs4 + DN(2,1)*clhs5;
const double clhs59 = N[2]*bdf0;
const double clhs60 = clhs58 + clhs59;
const double clhs61 = clhs16*clhs60 + clhs19*clhs60 + clhs58*clhs9;
const double clhs62 = C(0,1)*DN(2,1) + clhs53;
const double clhs63 = C(1,2)*DN(2,1);
const double clhs64 = C(2,2)*DN(2,0) + clhs63;
const double clhs65 = DN(2,1)*clhs24;
const double clhs66 = DN(0,0)*N[2];
const double clhs67 = clhs27*clhs59;
const double clhs68 = DN(2,0)*N[0];
const double clhs69 = C(0,1)*DN(0,0) + clhs22;
const double clhs70 = C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs71 = pow(DN(0,1), 2);
const double clhs72 = C(0,1)*DN(1,0) + clhs45;
const double clhs73 = DN(0,1)*clhs7;
const double clhs74 = DN(1,0)*clhs73;
const double clhs75 = C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs76 = DN(1,1)*clhs73 + clhs38;
const double clhs77 = DN(0,1)*N[1];
const double clhs78 = DN(1,1)*N[0];
const double clhs79 = C(0,1)*DN(2,0) + clhs63;
const double clhs80 = DN(2,0)*clhs73;
const double clhs81 = C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs82 = DN(2,1)*clhs73 + clhs56;
const double clhs83 = DN(0,1)*N[2];
const double clhs84 = DN(2,1)*N[0];
const double clhs85 = clhs13*clhs31 + clhs26;
const double clhs86 = bdf0*clhs27;
const double clhs87 = -N[1];
const double clhs88 = clhs31*clhs42 + clhs87;
const double clhs89 = DN(0,0)*clhs14;
const double clhs90 = DN(0,1)*clhs14;
const double clhs91 = DN(1,0)*clhs89 + DN(1,1)*clhs90 + N[1]*clhs28;
const double clhs92 = -N[2];
const double clhs93 = clhs31*clhs60 + clhs92;
const double clhs94 = DN(2,0)*clhs89 + DN(2,1)*clhs90 + N[2]*clhs28;
const double clhs95 = clhs15*clhs40;
const double clhs96 = N[1]*clhs18;
const double clhs97 = clhs13*clhs95 + clhs13*clhs96 + clhs37*clhs8;
const double clhs98 = DN(1,0)*clhs7;
const double clhs99 = clhs31*clhs40;
const double clhs100 = pow(DN(1,0), 2);
const double clhs101 = pow(N[1], 2);
const double clhs102 = clhs101*clhs11 + clhs37*clhs40 + clhs42*clhs95 + clhs42*clhs96;
const double clhs103 = DN(1,1)*clhs98;
const double clhs104 = clhs49*clhs7;
const double clhs105 = clhs104 + clhs30*clhs37 + clhs87 + clhs99;
const double clhs106 = clhs41*clhs55;
const double clhs107 = DN(2,0)*clhs98 + clhs106;
const double clhs108 = clhs37*clhs58 + clhs60*clhs95 + clhs60*clhs96;
const double clhs109 = DN(2,1)*clhs98;
const double clhs110 = DN(1,0)*N[2];
const double clhs111 = DN(2,0)*N[1];
const double clhs112 = DN(1,1)*clhs7;
const double clhs113 = pow(DN(1,1), 2);
const double clhs114 = DN(2,0)*clhs112;
const double clhs115 = DN(2,1)*clhs112 + clhs106;
const double clhs116 = DN(1,1)*N[2];
const double clhs117 = DN(2,1)*N[1];
const double clhs118 = DN(1,0)*DN(2,0)*clhs14 + DN(1,1)*DN(2,1)*clhs14 + N[2]*clhs49;
const double clhs119 = clhs15*clhs58;
const double clhs120 = N[2]*clhs18;
const double clhs121 = clhs119*clhs13 + clhs120*clhs13 + clhs55*clhs8;
const double clhs122 = DN(2,0)*clhs7;
const double clhs123 = clhs31*clhs58;
const double clhs124 = clhs119*clhs42 + clhs120*clhs42 + clhs40*clhs55;
const double clhs125 = pow(DN(2,0), 2);
const double clhs126 = pow(N[2], 2);
const double clhs127 = clhs11*clhs126 + clhs119*clhs60 + clhs120*clhs60 + clhs55*clhs58;
const double clhs128 = DN(2,1)*clhs122;
const double clhs129 = clhs123 + clhs30*clhs55 + clhs67*clhs7 + clhs92;
const double clhs130 = pow(DN(2,1), 2);
lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs20 + clhs3*clhs7;
lhs(0,1)=DN(0,0)*clhs21 + DN(0,1)*clhs23 + clhs25;
lhs(0,2)=DN(0,0)*clhs33;
lhs(0,3)=DN(0,0)*clhs34 + DN(0,1)*clhs36 + clhs39 + clhs43;
lhs(0,4)=DN(0,0)*clhs44 + DN(0,1)*clhs46 + clhs47;
lhs(0,5)=DN(1,0)*clhs32 + clhs24*clhs49 - clhs48 + clhs50*clhs51;
lhs(0,6)=DN(0,0)*clhs52 + DN(0,1)*clhs54 + clhs57 + clhs61;
lhs(0,7)=DN(0,0)*clhs62 + DN(0,1)*clhs64 + clhs65;
lhs(0,8)=DN(2,0)*clhs32 + clhs24*clhs67 + clhs51*clhs68 - clhs66;
lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs69 + clhs25;
lhs(1,1)=DN(0,0)*clhs23 + DN(0,1)*clhs70 + clhs20 + clhs7*clhs71;
lhs(1,2)=DN(0,1)*clhs33;
lhs(1,3)=DN(0,0)*clhs36 + DN(0,1)*clhs72 + clhs74;
lhs(1,4)=DN(0,0)*clhs46 + DN(0,1)*clhs75 + clhs43 + clhs76;
lhs(1,5)=DN(1,1)*clhs32 + clhs49*clhs73 + clhs51*clhs78 - clhs77;
lhs(1,6)=DN(0,0)*clhs54 + DN(0,1)*clhs79 + clhs80;
lhs(1,7)=DN(0,0)*clhs64 + DN(0,1)*clhs81 + clhs61 + clhs82;
lhs(1,8)=DN(2,1)*clhs32 + clhs51*clhs84 + clhs67*clhs73 - clhs83;
lhs(2,0)=DN(0,0)*clhs85;
lhs(2,1)=DN(0,1)*clhs85;
lhs(2,2)=clhs10*clhs86 + clhs14*clhs3 + clhs14*clhs71;
lhs(2,3)=DN(0,0)*clhs88;
lhs(2,4)=DN(0,1)*clhs88;
lhs(2,5)=clhs91;
lhs(2,6)=DN(0,0)*clhs93;
lhs(2,7)=DN(0,1)*clhs93;
lhs(2,8)=clhs94;
lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs39 + clhs97;
lhs(3,1)=DN(1,0)*clhs21 + DN(1,1)*clhs23 + clhs74;
lhs(3,2)=DN(0,0)*clhs99 + clhs28*clhs98 + clhs48*clhs51 - clhs50;
lhs(3,3)=DN(1,0)*clhs34 + DN(1,1)*clhs36 + clhs100*clhs7 + clhs102;
lhs(3,4)=DN(1,0)*clhs44 + DN(1,1)*clhs46 + clhs103;
lhs(3,5)=DN(1,0)*clhs105;
lhs(3,6)=DN(1,0)*clhs52 + DN(1,1)*clhs54 + clhs107 + clhs108;
lhs(3,7)=DN(1,0)*clhs62 + DN(1,1)*clhs64 + clhs109;
lhs(3,8)=DN(2,0)*clhs99 - clhs110 + clhs111*clhs51 + clhs67*clhs98;
lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs69 + clhs47;
lhs(4,1)=DN(1,0)*clhs23 + DN(1,1)*clhs70 + clhs76 + clhs97;
lhs(4,2)=DN(0,1)*clhs99 + clhs112*clhs28 + clhs51*clhs77 - clhs78;
lhs(4,3)=DN(1,0)*clhs36 + DN(1,1)*clhs72 + clhs103;
lhs(4,4)=DN(1,0)*clhs46 + DN(1,1)*clhs75 + clhs102 + clhs113*clhs7;
lhs(4,5)=DN(1,1)*clhs105;
lhs(4,6)=DN(1,0)*clhs54 + DN(1,1)*clhs79 + clhs114;
lhs(4,7)=DN(1,0)*clhs64 + DN(1,1)*clhs81 + clhs108 + clhs115;
lhs(4,8)=DN(2,1)*clhs99 + clhs112*clhs67 - clhs116 + clhs117*clhs51;
lhs(5,0)=DN(1,0)*clhs85;
lhs(5,1)=DN(1,1)*clhs85;
lhs(5,2)=clhs91;
lhs(5,3)=DN(1,0)*clhs88;
lhs(5,4)=DN(1,1)*clhs88;
lhs(5,5)=clhs100*clhs14 + clhs101*clhs86 + clhs113*clhs14;
lhs(5,6)=DN(1,0)*clhs93;
lhs(5,7)=DN(1,1)*clhs93;
lhs(5,8)=clhs118;
lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs121 + clhs57;
lhs(6,1)=DN(2,0)*clhs21 + DN(2,1)*clhs23 + clhs80;
lhs(6,2)=DN(0,0)*clhs123 + clhs122*clhs28 + clhs51*clhs66 - clhs68;
lhs(6,3)=DN(2,0)*clhs34 + DN(2,1)*clhs36 + clhs107 + clhs124;
lhs(6,4)=DN(2,0)*clhs44 + DN(2,1)*clhs46 + clhs114;
lhs(6,5)=DN(1,0)*clhs123 + clhs110*clhs51 - clhs111 + clhs122*clhs49;
lhs(6,6)=DN(2,0)*clhs52 + DN(2,1)*clhs54 + clhs125*clhs7 + clhs127;
lhs(6,7)=DN(2,0)*clhs62 + DN(2,1)*clhs64 + clhs128;
lhs(6,8)=DN(2,0)*clhs129;
lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs69 + clhs65;
lhs(7,1)=DN(2,0)*clhs23 + DN(2,1)*clhs70 + clhs121 + clhs82;
lhs(7,2)=DN(0,1)*clhs123 + DN(2,1)*clhs29 + clhs51*clhs83 - clhs84;
lhs(7,3)=DN(2,0)*clhs36 + DN(2,1)*clhs72 + clhs109;
lhs(7,4)=DN(2,0)*clhs46 + DN(2,1)*clhs75 + clhs115 + clhs124;
lhs(7,5)=DN(1,1)*clhs123 + DN(2,1)*clhs104 + clhs116*clhs51 - clhs117;
lhs(7,6)=DN(2,0)*clhs54 + DN(2,1)*clhs79 + clhs128;
lhs(7,7)=DN(2,0)*clhs64 + DN(2,1)*clhs81 + clhs127 + clhs130*clhs7;
lhs(7,8)=DN(2,1)*clhs129;
lhs(8,0)=DN(2,0)*clhs85;
lhs(8,1)=DN(2,1)*clhs85;
lhs(8,2)=clhs94;
lhs(8,3)=DN(2,0)*clhs88;
lhs(8,4)=DN(2,1)*clhs88;
lhs(8,5)=clhs118;
lhs(8,6)=DN(2,0)*clhs93;
lhs(8,7)=DN(2,1)*clhs93;
lhs(8,8)=clhs125*clhs14 + clhs126*clhs86 + clhs130*clhs14;


}


template<>
void EmbeddedAusasNavierStokes<3>::ComputeGaussPointRHSContribution(
    array_1d<double,16>& rhs,
    const EmbeddedAusasElementDataStruct& data)
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

    const double crhs0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 = rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crhs2 = rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)));
const double crhs3 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double crhs4 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
const double crhs5 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
const double crhs6 = N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
const double crhs7 = rho*(crhs3*crhs4 + crhs5*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs6*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)));
const double crhs8 = rho*stab_c2*sqrt(pow(crhs4, 2) + pow(crhs5, 2) + pow(crhs6, 2));
const double crhs9 = (N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]) + N[3]*(bdf0*p[3] + bdf1*pn[3] + bdf2*pnn[3]))/(pow(c, 2)*rho);
const double crhs10 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double crhs11 = DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs12 = (crhs8*h/stab_c1 + mu)*(crhs10 + crhs11 + crhs3 + crhs9);
const double crhs13 = 1.0/(crhs8/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs14 = crhs13*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + crhs2 + crhs7);
const double crhs15 = rho*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2));
const double crhs16 = N[0]*crhs15;
const double crhs17 = rho*(DN(0,0)*crhs4 + DN(0,1)*crhs5 + DN(0,2)*crhs6);
const double crhs18 = rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs19 = rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crhs20 = rho*(crhs10*crhs5 + crhs4*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs6*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
const double crhs21 = crhs13*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs18 + crhs19 + crhs20);
const double crhs22 = rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs23 = rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crhs24 = rho*(crhs11*crhs6 + crhs4*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs5*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)));
const double crhs25 = crhs13*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs22 + crhs23 + crhs24);
const double crhs26 = N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0);
const double crhs27 = N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
const double crhs28 = N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2);
const double crhs29 = N[1]*crhs15;
const double crhs30 = rho*(DN(1,0)*crhs4 + DN(1,1)*crhs5 + DN(1,2)*crhs6);
const double crhs31 = N[2]*crhs15;
const double crhs32 = rho*(DN(2,0)*crhs4 + DN(2,1)*crhs5 + DN(2,2)*crhs6);
const double crhs33 = N[3]*crhs15;
const double crhs34 = rho*(DN(3,0)*crhs4 + DN(3,1)*crhs5 + DN(3,2)*crhs6);
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs12 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - N[0]*crhs2 - N[0]*crhs7 - crhs14*crhs16 - crhs14*crhs17;
rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs12 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs18 - N[0]*crhs19 - N[0]*crhs20 - crhs16*crhs21 - crhs17*crhs21;
rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs12 - DN(0,2)*stress[2] + N[0]*crhs22 - N[0]*crhs23 - N[0]*crhs24 - crhs16*crhs25 - crhs17*crhs25;
rhs[3]=-DN(0,0)*crhs14 + DN(0,0)*crhs26 - DN(0,1)*crhs21 + DN(0,1)*crhs27 - DN(0,2)*crhs25 + DN(0,2)*crhs28 - N[0]*crhs9;
rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs12 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - N[1]*crhs2 - N[1]*crhs7 - crhs14*crhs29 - crhs14*crhs30;
rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs12 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs18 - N[1]*crhs19 - N[1]*crhs20 - crhs21*crhs29 - crhs21*crhs30;
rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs12 - DN(1,2)*stress[2] + N[1]*crhs22 - N[1]*crhs23 - N[1]*crhs24 - crhs25*crhs29 - crhs25*crhs30;
rhs[7]=-DN(1,0)*crhs14 + DN(1,0)*crhs26 - DN(1,1)*crhs21 + DN(1,1)*crhs27 - DN(1,2)*crhs25 + DN(1,2)*crhs28 - N[1]*crhs9;
rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs12 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - N[2]*crhs2 - N[2]*crhs7 - crhs14*crhs31 - crhs14*crhs32;
rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs12 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs18 - N[2]*crhs19 - N[2]*crhs20 - crhs21*crhs31 - crhs21*crhs32;
rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs12 - DN(2,2)*stress[2] + N[2]*crhs22 - N[2]*crhs23 - N[2]*crhs24 - crhs25*crhs31 - crhs25*crhs32;
rhs[11]=-DN(2,0)*crhs14 + DN(2,0)*crhs26 - DN(2,1)*crhs21 + DN(2,1)*crhs27 - DN(2,2)*crhs25 + DN(2,2)*crhs28 - N[2]*crhs9;
rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs12 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - N[3]*crhs2 - N[3]*crhs7 - crhs14*crhs33 - crhs14*crhs34;
rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs12 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs18 - N[3]*crhs19 - N[3]*crhs20 - crhs21*crhs33 - crhs21*crhs34;
rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs12 - DN(3,2)*stress[2] + N[3]*crhs22 - N[3]*crhs23 - N[3]*crhs24 - crhs25*crhs33 - crhs25*crhs34;
rhs[15]=-DN(3,0)*crhs14 + DN(3,0)*crhs26 - DN(3,1)*crhs21 + DN(3,1)*crhs27 - DN(3,2)*crhs25 + DN(3,2)*crhs28 - N[3]*crhs9;


}


template<>
void EmbeddedAusasNavierStokes<2>::ComputeGaussPointRHSContribution(
    array_1d<double,9>& rhs,
    const EmbeddedAusasElementDataStruct& data)
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

    const double crhs0 = N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 = rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crhs2 = rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crhs3 = DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs4 = N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs5 = N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs6 = rho*(crhs3*crhs4 + crhs5*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crhs7 = rho*stab_c2*sqrt(pow(crhs4, 2) + pow(crhs5, 2));
const double crhs8 = (N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]))/(pow(c, 2)*rho);
const double crhs9 = DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs10 = (crhs7*h/stab_c1 + mu)*(crhs3 + crhs8 + crhs9);
const double crhs11 = 1.0/(crhs7/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs12 = crhs11*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + crhs2 + crhs6);
const double crhs13 = rho*(DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1));
const double crhs14 = N[0]*crhs13;
const double crhs15 = rho*(DN(0,0)*crhs4 + DN(0,1)*crhs5);
const double crhs16 = rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs17 = rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crhs18 = rho*(crhs4*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + crhs5*crhs9);
const double crhs19 = crhs11*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs16 + crhs17 + crhs18);
const double crhs20 = N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0);
const double crhs21 = N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
const double crhs22 = N[1]*crhs13;
const double crhs23 = rho*(DN(1,0)*crhs4 + DN(1,1)*crhs5);
const double crhs24 = N[2]*crhs13;
const double crhs25 = rho*(DN(2,0)*crhs4 + DN(2,1)*crhs5);
rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs10 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - N[0]*crhs2 - N[0]*crhs6 - crhs12*crhs14 - crhs12*crhs15;
rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs10 - DN(0,1)*stress[1] + N[0]*crhs16 - N[0]*crhs17 - N[0]*crhs18 - crhs14*crhs19 - crhs15*crhs19;
rhs[2]=-DN(0,0)*crhs12 + DN(0,0)*crhs20 - DN(0,1)*crhs19 + DN(0,1)*crhs21 - N[0]*crhs8;
rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs10 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - N[1]*crhs2 - N[1]*crhs6 - crhs12*crhs22 - crhs12*crhs23;
rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs10 - DN(1,1)*stress[1] + N[1]*crhs16 - N[1]*crhs17 - N[1]*crhs18 - crhs19*crhs22 - crhs19*crhs23;
rhs[5]=-DN(1,0)*crhs12 + DN(1,0)*crhs20 - DN(1,1)*crhs19 + DN(1,1)*crhs21 - N[1]*crhs8;
rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs10 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - N[2]*crhs2 - N[2]*crhs6 - crhs12*crhs24 - crhs12*crhs25;
rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs10 - DN(2,1)*stress[1] + N[2]*crhs16 - N[2]*crhs17 - N[2]*crhs18 - crhs19*crhs24 - crhs19*crhs25;
rhs[8]=-DN(2,0)*crhs12 + DN(2,0)*crhs20 - DN(2,1)*crhs19 + DN(2,1)*crhs21 - N[2]*crhs8;


}

}
