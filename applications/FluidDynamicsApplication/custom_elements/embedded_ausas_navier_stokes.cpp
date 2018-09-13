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
    ProcessInfo& rCurrentProcessInfo)
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
void EmbeddedAusasNavierStokes<3>::GetDofList(
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
void EmbeddedAusasNavierStokes<2>::GetDofList(
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
const double clhs36 =             -N[0];
const double clhs37 =             pow(c, -2);
const double clhs38 =             1.0/rho;
const double clhs39 =             N[0]*bdf0*clhs37*clhs38;
const double clhs40 =             1.0*clhs17*clhs18;
const double clhs41 =             1.0*clhs18*rho;
const double clhs42 =             clhs14*clhs41;
const double clhs43 =             clhs10*clhs39 + clhs13*clhs40 + clhs36 + clhs42;
const double clhs44 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs45 =             C(0,3)*DN(1,0);
const double clhs46 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs45;
const double clhs47 =             C(0,5)*DN(1,0);
const double clhs48 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs47;
const double clhs49 =             DN(1,0)*clhs28;
const double clhs50 =             N[0]*bdf0*rho;
const double clhs51 =             N[1]*clhs50;
const double clhs52 =             DN(1,0)*clhs6 + DN(1,1)*clhs7 + DN(1,2)*clhs8;
const double clhs53 =             N[1]*bdf0 + clhs52;
const double clhs54 =             clhs13*clhs52 + clhs19*clhs53 + clhs20*clhs53 + clhs51;
const double clhs55 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs45;
const double clhs56 =             C(1,3)*DN(1,1);
const double clhs57 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs56;
const double clhs58 =             C(3,5)*DN(1,0);
const double clhs59 =             C(4,5)*DN(1,2);
const double clhs60 =             C(1,5)*DN(1,1) + clhs58 + clhs59;
const double clhs61 =             DN(1,1)*clhs28;
const double clhs62 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs47;
const double clhs63 =             C(3,4)*DN(1,1);
const double clhs64 =             C(2,3)*DN(1,2) + clhs58 + clhs63;
const double clhs65 =             C(2,5)*DN(1,2);
const double clhs66 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs65;
const double clhs67 =             DN(1,2)*clhs28;
const double clhs68 =             DN(0,0)*N[1];
const double clhs69 =             bdf0*clhs10*clhs37*clhs38;
const double clhs70 =             DN(1,0)*N[0];
const double clhs71 =             1.0*clhs17*clhs18*rho;
const double clhs72 =             1.0*DN(1,0)*clhs18*rho;
const double clhs73 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs74 =             C(0,3)*DN(2,0);
const double clhs75 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs74;
const double clhs76 =             C(0,5)*DN(2,0);
const double clhs77 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs76;
const double clhs78 =             DN(2,0)*clhs28;
const double clhs79 =             N[2]*clhs50;
const double clhs80 =             DN(2,0)*clhs6 + DN(2,1)*clhs7 + DN(2,2)*clhs8;
const double clhs81 =             N[2]*bdf0;
const double clhs82 =             clhs80 + clhs81;
const double clhs83 =             clhs13*clhs80 + clhs19*clhs82 + clhs20*clhs82 + clhs79;
const double clhs84 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs74;
const double clhs85 =             C(1,3)*DN(2,1);
const double clhs86 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs85;
const double clhs87 =             C(3,5)*DN(2,0);
const double clhs88 =             C(4,5)*DN(2,2);
const double clhs89 =             C(1,5)*DN(2,1) + clhs87 + clhs88;
const double clhs90 =             DN(2,1)*clhs28;
const double clhs91 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs76;
const double clhs92 =             C(3,4)*DN(2,1);
const double clhs93 =             C(2,3)*DN(2,2) + clhs87 + clhs92;
const double clhs94 =             C(2,5)*DN(2,2);
const double clhs95 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs94;
const double clhs96 =             DN(2,2)*clhs28;
const double clhs97 =             DN(0,0)*N[2];
const double clhs98 =             DN(2,0)*N[0];
const double clhs99 =             1.0*DN(2,0)*clhs18*rho;
const double clhs100 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs101 =             C(0,3)*DN(3,0);
const double clhs102 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs101;
const double clhs103 =             C(0,5)*DN(3,0);
const double clhs104 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs103;
const double clhs105 =             DN(3,0)*clhs28;
const double clhs106 =             N[3]*clhs50;
const double clhs107 =             DN(3,0)*clhs6 + DN(3,1)*clhs7 + DN(3,2)*clhs8;
const double clhs108 =             N[3]*bdf0;
const double clhs109 =             clhs107 + clhs108;
const double clhs110 =             clhs106 + clhs107*clhs13 + clhs109*clhs19 + clhs109*clhs20;
const double clhs111 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs101;
const double clhs112 =             C(1,3)*DN(3,1);
const double clhs113 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs112;
const double clhs114 =             C(3,5)*DN(3,0);
const double clhs115 =             C(4,5)*DN(3,2);
const double clhs116 =             C(1,5)*DN(3,1) + clhs114 + clhs115;
const double clhs117 =             DN(3,1)*clhs28;
const double clhs118 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs103;
const double clhs119 =             C(3,4)*DN(3,1);
const double clhs120 =             C(2,3)*DN(3,2) + clhs114 + clhs119;
const double clhs121 =             C(2,5)*DN(3,2);
const double clhs122 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs121;
const double clhs123 =             DN(3,2)*clhs28;
const double clhs124 =             DN(0,0)*N[3];
const double clhs125 =             DN(3,0)*N[0];
const double clhs126 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs23;
const double clhs127 =             C(0,4)*DN(0,0) + clhs26 + clhs31;
const double clhs128 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs129 =             C(1,4)*DN(0,1);
const double clhs130 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs129;
const double clhs131 =             pow(DN(0,1), 2);
const double clhs132 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs129;
const double clhs133 =             C(2,4)*DN(0,2);
const double clhs134 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs133;
const double clhs135 =             DN(0,1)*clhs10;
const double clhs136 =             DN(0,2)*clhs135;
const double clhs137 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs56;
const double clhs138 =             C(0,4)*DN(1,0) + clhs59 + clhs63;
const double clhs139 =             DN(1,0)*clhs135;
const double clhs140 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs141 =             C(1,4)*DN(1,1);
const double clhs142 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs141;
const double clhs143 =             DN(1,1)*clhs135;
const double clhs144 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs141;
const double clhs145 =             C(2,4)*DN(1,2);
const double clhs146 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs145;
const double clhs147 =             DN(1,2)*clhs135;
const double clhs148 =             DN(0,1)*N[1];
const double clhs149 =             DN(1,1)*N[0];
const double clhs150 =             1.0*DN(1,1)*clhs18*rho;
const double clhs151 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs85;
const double clhs152 =             C(0,4)*DN(2,0) + clhs88 + clhs92;
const double clhs153 =             DN(2,0)*clhs135;
const double clhs154 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs155 =             C(1,4)*DN(2,1);
const double clhs156 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs155;
const double clhs157 =             DN(2,1)*clhs135;
const double clhs158 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs155;
const double clhs159 =             C(2,4)*DN(2,2);
const double clhs160 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs159;
const double clhs161 =             DN(2,2)*clhs135;
const double clhs162 =             DN(0,1)*N[2];
const double clhs163 =             DN(2,1)*N[0];
const double clhs164 =             1.0*DN(2,1)*clhs18*rho;
const double clhs165 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs112;
const double clhs166 =             C(0,4)*DN(3,0) + clhs115 + clhs119;
const double clhs167 =             DN(3,0)*clhs135;
const double clhs168 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs169 =             C(1,4)*DN(3,1);
const double clhs170 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs169;
const double clhs171 =             DN(3,1)*clhs135;
const double clhs172 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs169;
const double clhs173 =             C(2,4)*DN(3,2);
const double clhs174 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs173;
const double clhs175 =             DN(3,2)*clhs135;
const double clhs176 =             DN(0,1)*N[3];
const double clhs177 =             DN(3,1)*N[0];
const double clhs178 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs33;
const double clhs179 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs133;
const double clhs180 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs181 =             pow(DN(0,2), 2);
const double clhs182 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs65;
const double clhs183 =             DN(0,2)*clhs10;
const double clhs184 =             DN(1,0)*clhs183;
const double clhs185 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs145;
const double clhs186 =             DN(1,1)*clhs183;
const double clhs187 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs188 =             DN(1,2)*clhs183;
const double clhs189 =             DN(0,2)*N[1];
const double clhs190 =             DN(1,2)*N[0];
const double clhs191 =             1.0*DN(1,2)*clhs18*rho;
const double clhs192 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs94;
const double clhs193 =             DN(2,0)*clhs183;
const double clhs194 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs159;
const double clhs195 =             DN(2,1)*clhs183;
const double clhs196 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs197 =             DN(2,2)*clhs183;
const double clhs198 =             DN(0,2)*N[2];
const double clhs199 =             DN(2,2)*N[0];
const double clhs200 =             1.0*DN(2,2)*clhs18*rho;
const double clhs201 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs121;
const double clhs202 =             DN(3,0)*clhs183;
const double clhs203 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs173;
const double clhs204 =             DN(3,1)*clhs183;
const double clhs205 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs206 =             DN(3,2)*clhs183;
const double clhs207 =             DN(0,2)*N[3];
const double clhs208 =             DN(3,2)*N[0];
const double clhs209 =             clhs15*clhs41 + clhs36;
const double clhs210 =             bdf0*clhs37*clhs38;
const double clhs211 =             1.0*clhs18;
const double clhs212 =             -N[1];
const double clhs213 =             clhs212 + clhs41*clhs53;
const double clhs214 =             1.0*DN(0,0)*clhs18;
const double clhs215 =             1.0*DN(0,1)*clhs18;
const double clhs216 =             1.0*DN(0,2)*clhs18;
const double clhs217 =             DN(1,0)*clhs214 + DN(1,1)*clhs215 + DN(1,2)*clhs216 + N[1]*clhs39;
const double clhs218 =             -N[2];
const double clhs219 =             clhs218 + clhs41*clhs82;
const double clhs220 =             DN(2,0)*clhs214 + DN(2,1)*clhs215 + DN(2,2)*clhs216 + N[2]*clhs39;
const double clhs221 =             -N[3];
const double clhs222 =             clhs109*clhs41 + clhs221;
const double clhs223 =             DN(3,0)*clhs214 + DN(3,1)*clhs215 + DN(3,2)*clhs216 + N[3]*clhs39;
const double clhs224 =             N[1]*rho;
const double clhs225 =             1.0*N[1]*clhs16*clhs17*clhs18;
const double clhs226 =             1.0*clhs16*clhs18*clhs52;
const double clhs227 =             clhs14*clhs224 + clhs15*clhs225 + clhs15*clhs226 + clhs51;
const double clhs228 =             1.0*DN(0,0)*clhs18*rho;
const double clhs229 =             pow(DN(1,0), 2);
const double clhs230 =             pow(N[1], 2);
const double clhs231 =             clhs12*clhs230 + clhs224*clhs52 + clhs225*clhs53 + clhs226*clhs53;
const double clhs232 =             DN(1,0)*clhs10;
const double clhs233 =             DN(1,1)*clhs232;
const double clhs234 =             DN(1,2)*clhs232;
const double clhs235 =             N[1]*bdf0*clhs37*clhs38;
const double clhs236 =             clhs41*clhs52;
const double clhs237 =             clhs10*clhs235 + clhs212 + clhs224*clhs40 + clhs236;
const double clhs238 =             DN(2,0)*clhs232;
const double clhs239 =             N[1]*bdf0*rho;
const double clhs240 =             N[2]*clhs239;
const double clhs241 =             clhs224*clhs80 + clhs225*clhs82 + clhs226*clhs82 + clhs240;
const double clhs242 =             DN(2,1)*clhs232;
const double clhs243 =             DN(2,2)*clhs232;
const double clhs244 =             DN(1,0)*N[2];
const double clhs245 =             DN(2,0)*N[1];
const double clhs246 =             DN(3,0)*clhs232;
const double clhs247 =             N[3]*clhs239;
const double clhs248 =             clhs107*clhs224 + clhs109*clhs225 + clhs109*clhs226 + clhs247;
const double clhs249 =             DN(3,1)*clhs232;
const double clhs250 =             DN(3,2)*clhs232;
const double clhs251 =             DN(1,0)*N[3];
const double clhs252 =             DN(3,0)*N[1];
const double clhs253 =             1.0*DN(0,1)*clhs18*rho;
const double clhs254 =             pow(DN(1,1), 2);
const double clhs255 =             DN(1,1)*clhs10;
const double clhs256 =             DN(1,2)*clhs255;
const double clhs257 =             DN(2,0)*clhs255;
const double clhs258 =             DN(2,1)*clhs255;
const double clhs259 =             DN(2,2)*clhs255;
const double clhs260 =             DN(1,1)*N[2];
const double clhs261 =             DN(2,1)*N[1];
const double clhs262 =             DN(3,0)*clhs255;
const double clhs263 =             DN(3,1)*clhs255;
const double clhs264 =             DN(3,2)*clhs255;
const double clhs265 =             DN(1,1)*N[3];
const double clhs266 =             DN(3,1)*N[1];
const double clhs267 =             1.0*DN(0,2)*clhs18*rho;
const double clhs268 =             pow(DN(1,2), 2);
const double clhs269 =             DN(1,2)*clhs10;
const double clhs270 =             DN(2,0)*clhs269;
const double clhs271 =             DN(2,1)*clhs269;
const double clhs272 =             DN(2,2)*clhs269;
const double clhs273 =             DN(1,2)*N[2];
const double clhs274 =             DN(2,2)*N[1];
const double clhs275 =             DN(3,0)*clhs269;
const double clhs276 =             DN(3,1)*clhs269;
const double clhs277 =             DN(3,2)*clhs269;
const double clhs278 =             DN(1,2)*N[3];
const double clhs279 =             DN(3,2)*N[1];
const double clhs280 =             1.0*DN(1,0)*clhs18;
const double clhs281 =             1.0*DN(1,1)*clhs18;
const double clhs282 =             1.0*DN(1,2)*clhs18;
const double clhs283 =             DN(2,0)*clhs280 + DN(2,1)*clhs281 + DN(2,2)*clhs282 + N[2]*clhs235;
const double clhs284 =             DN(3,0)*clhs280 + DN(3,1)*clhs281 + DN(3,2)*clhs282 + N[3]*clhs235;
const double clhs285 =             N[2]*rho;
const double clhs286 =             1.0*N[2]*clhs16*clhs17*clhs18;
const double clhs287 =             1.0*clhs16*clhs18*clhs80;
const double clhs288 =             clhs14*clhs285 + clhs15*clhs286 + clhs15*clhs287 + clhs79;
const double clhs289 =             clhs240 + clhs285*clhs52 + clhs286*clhs53 + clhs287*clhs53;
const double clhs290 =             pow(DN(2,0), 2);
const double clhs291 =             pow(N[2], 2);
const double clhs292 =             clhs12*clhs291 + clhs285*clhs80 + clhs286*clhs82 + clhs287*clhs82;
const double clhs293 =             DN(2,0)*clhs10;
const double clhs294 =             DN(2,1)*clhs293;
const double clhs295 =             DN(2,2)*clhs293;
const double clhs296 =             clhs10*clhs37*clhs38;
const double clhs297 =             clhs41*clhs80;
const double clhs298 =             clhs218 + clhs285*clhs40 + clhs296*clhs81 + clhs297;
const double clhs299 =             DN(3,0)*clhs293;
const double clhs300 =             N[2]*N[3]*bdf0;
const double clhs301 =             clhs300*rho;
const double clhs302 =             clhs107*clhs285 + clhs109*clhs286 + clhs109*clhs287 + clhs301;
const double clhs303 =             DN(3,1)*clhs293;
const double clhs304 =             DN(3,2)*clhs293;
const double clhs305 =             DN(2,0)*N[3];
const double clhs306 =             DN(3,0)*N[2];
const double clhs307 =             pow(DN(2,1), 2);
const double clhs308 =             DN(2,1)*clhs10;
const double clhs309 =             DN(2,2)*clhs308;
const double clhs310 =             DN(3,0)*clhs308;
const double clhs311 =             DN(3,1)*clhs308;
const double clhs312 =             DN(3,2)*clhs308;
const double clhs313 =             DN(2,1)*N[3];
const double clhs314 =             DN(3,1)*N[2];
const double clhs315 =             pow(DN(2,2), 2);
const double clhs316 =             DN(2,2)*clhs10;
const double clhs317 =             DN(3,0)*clhs316;
const double clhs318 =             DN(3,1)*clhs316;
const double clhs319 =             DN(3,2)*clhs316;
const double clhs320 =             DN(2,2)*N[3];
const double clhs321 =             DN(3,2)*N[2];
const double clhs322 =             1.0*DN(2,0)*DN(3,0)*clhs18 + 1.0*DN(2,1)*DN(3,1)*clhs18 + 1.0*DN(2,2)*DN(3,2)*clhs18 + clhs300*clhs37*clhs38;
const double clhs323 =             N[3]*rho;
const double clhs324 =             1.0*N[3]*clhs16*clhs17*clhs18;
const double clhs325 =             1.0*clhs107*clhs16*clhs18;
const double clhs326 =             clhs106 + clhs14*clhs323 + clhs15*clhs324 + clhs15*clhs325;
const double clhs327 =             clhs247 + clhs323*clhs52 + clhs324*clhs53 + clhs325*clhs53;
const double clhs328 =             clhs301 + clhs323*clhs80 + clhs324*clhs82 + clhs325*clhs82;
const double clhs329 =             pow(DN(3,0), 2);
const double clhs330 =             pow(N[3], 2);
const double clhs331 =             clhs107*clhs323 + clhs109*clhs324 + clhs109*clhs325 + clhs12*clhs330;
const double clhs332 =             DN(3,0)*clhs10;
const double clhs333 =             DN(3,1)*clhs332;
const double clhs334 =             DN(3,2)*clhs332;
const double clhs335 =             clhs107*clhs41 + clhs108*clhs296 + clhs221 + clhs323*clhs40;
const double clhs336 =             pow(DN(3,1), 2);
const double clhs337 =             DN(3,1)*DN(3,2)*clhs10;
const double clhs338 =             pow(DN(3,2), 2);
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + DN(0,2)*clhs4 + clhs10*clhs5 + clhs21;
            lhs(0,1)=DN(0,0)*clhs22 + DN(0,1)*clhs24 + DN(0,2)*clhs27 + clhs29;
            lhs(0,2)=DN(0,0)*clhs30 + DN(0,1)*clhs32 + DN(0,2)*clhs34 + clhs35;
            lhs(0,3)=DN(0,0)*clhs43;
            lhs(0,4)=DN(0,0)*clhs44 + DN(0,1)*clhs46 + DN(0,2)*clhs48 + clhs49 + clhs54;
            lhs(0,5)=DN(0,0)*clhs55 + DN(0,1)*clhs57 + DN(0,2)*clhs60 + clhs61;
            lhs(0,6)=DN(0,0)*clhs62 + DN(0,1)*clhs64 + DN(0,2)*clhs66 + clhs67;
            lhs(0,7)=clhs14*clhs72 + clhs68*clhs69 - clhs68 + clhs70*clhs71;
            lhs(0,8)=DN(0,0)*clhs73 + DN(0,1)*clhs75 + DN(0,2)*clhs77 + clhs78 + clhs83;
            lhs(0,9)=DN(0,0)*clhs84 + DN(0,1)*clhs86 + DN(0,2)*clhs89 + clhs90;
            lhs(0,10)=DN(0,0)*clhs91 + DN(0,1)*clhs93 + DN(0,2)*clhs95 + clhs96;
            lhs(0,11)=clhs14*clhs99 + clhs69*clhs97 + clhs71*clhs98 - clhs97;
            lhs(0,12)=DN(0,0)*clhs100 + DN(0,1)*clhs102 + DN(0,2)*clhs104 + clhs105 + clhs110;
            lhs(0,13)=DN(0,0)*clhs111 + DN(0,1)*clhs113 + DN(0,2)*clhs116 + clhs117;
            lhs(0,14)=DN(0,0)*clhs118 + DN(0,1)*clhs120 + DN(0,2)*clhs122 + clhs123;
            lhs(0,15)=DN(3,0)*clhs42 + clhs124*clhs69 - clhs124 + clhs125*clhs71;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs126 + DN(0,2)*clhs127 + clhs29;
            lhs(1,1)=DN(0,0)*clhs24 + DN(0,1)*clhs128 + DN(0,2)*clhs130 + clhs10*clhs131 + clhs21;
            lhs(1,2)=DN(0,0)*clhs32 + DN(0,1)*clhs132 + DN(0,2)*clhs134 + clhs136;
            lhs(1,3)=DN(0,1)*clhs43;
            lhs(1,4)=DN(0,0)*clhs46 + DN(0,1)*clhs137 + DN(0,2)*clhs138 + clhs139;
            lhs(1,5)=DN(0,0)*clhs57 + DN(0,1)*clhs140 + DN(0,2)*clhs142 + clhs143 + clhs54;
            lhs(1,6)=DN(0,0)*clhs64 + DN(0,1)*clhs144 + DN(0,2)*clhs146 + clhs147;
            lhs(1,7)=clhs14*clhs150 + clhs148*clhs69 - clhs148 + clhs149*clhs71;
            lhs(1,8)=DN(0,0)*clhs75 + DN(0,1)*clhs151 + DN(0,2)*clhs152 + clhs153;
            lhs(1,9)=DN(0,0)*clhs86 + DN(0,1)*clhs154 + DN(0,2)*clhs156 + clhs157 + clhs83;
            lhs(1,10)=DN(0,0)*clhs93 + DN(0,1)*clhs158 + DN(0,2)*clhs160 + clhs161;
            lhs(1,11)=clhs14*clhs164 + clhs162*clhs69 - clhs162 + clhs163*clhs71;
            lhs(1,12)=DN(0,0)*clhs102 + DN(0,1)*clhs165 + DN(0,2)*clhs166 + clhs167;
            lhs(1,13)=DN(0,0)*clhs113 + DN(0,1)*clhs168 + DN(0,2)*clhs170 + clhs110 + clhs171;
            lhs(1,14)=DN(0,0)*clhs120 + DN(0,1)*clhs172 + DN(0,2)*clhs174 + clhs175;
            lhs(1,15)=DN(3,1)*clhs42 + clhs176*clhs69 - clhs176 + clhs177*clhs71;
            lhs(2,0)=DN(0,0)*clhs4 + DN(0,1)*clhs127 + DN(0,2)*clhs178 + clhs35;
            lhs(2,1)=DN(0,0)*clhs27 + DN(0,1)*clhs130 + DN(0,2)*clhs179 + clhs136;
            lhs(2,2)=DN(0,0)*clhs34 + DN(0,1)*clhs134 + DN(0,2)*clhs180 + clhs10*clhs181 + clhs21;
            lhs(2,3)=DN(0,2)*clhs43;
            lhs(2,4)=DN(0,0)*clhs48 + DN(0,1)*clhs138 + DN(0,2)*clhs182 + clhs184;
            lhs(2,5)=DN(0,0)*clhs60 + DN(0,1)*clhs142 + DN(0,2)*clhs185 + clhs186;
            lhs(2,6)=DN(0,0)*clhs66 + DN(0,1)*clhs146 + DN(0,2)*clhs187 + clhs188 + clhs54;
            lhs(2,7)=clhs14*clhs191 + clhs189*clhs69 - clhs189 + clhs190*clhs71;
            lhs(2,8)=DN(0,0)*clhs77 + DN(0,1)*clhs152 + DN(0,2)*clhs192 + clhs193;
            lhs(2,9)=DN(0,0)*clhs89 + DN(0,1)*clhs156 + DN(0,2)*clhs194 + clhs195;
            lhs(2,10)=DN(0,0)*clhs95 + DN(0,1)*clhs160 + DN(0,2)*clhs196 + clhs197 + clhs83;
            lhs(2,11)=clhs14*clhs200 + clhs198*clhs69 - clhs198 + clhs199*clhs71;
            lhs(2,12)=DN(0,0)*clhs104 + DN(0,1)*clhs166 + DN(0,2)*clhs201 + clhs202;
            lhs(2,13)=DN(0,0)*clhs116 + DN(0,1)*clhs170 + DN(0,2)*clhs203 + clhs204;
            lhs(2,14)=DN(0,0)*clhs122 + DN(0,1)*clhs174 + DN(0,2)*clhs205 + clhs110 + clhs206;
            lhs(2,15)=DN(3,2)*clhs42 + clhs207*clhs69 - clhs207 + clhs208*clhs71;
            lhs(3,0)=DN(0,0)*clhs209;
            lhs(3,1)=DN(0,1)*clhs209;
            lhs(3,2)=DN(0,2)*clhs209;
            lhs(3,3)=clhs11*clhs210 + clhs131*clhs211 + clhs181*clhs211 + clhs211*clhs5;
            lhs(3,4)=DN(0,0)*clhs213;
            lhs(3,5)=DN(0,1)*clhs213;
            lhs(3,6)=DN(0,2)*clhs213;
            lhs(3,7)=clhs217;
            lhs(3,8)=DN(0,0)*clhs219;
            lhs(3,9)=DN(0,1)*clhs219;
            lhs(3,10)=DN(0,2)*clhs219;
            lhs(3,11)=clhs220;
            lhs(3,12)=DN(0,0)*clhs222;
            lhs(3,13)=DN(0,1)*clhs222;
            lhs(3,14)=DN(0,2)*clhs222;
            lhs(3,15)=clhs223;
            lhs(4,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + DN(1,2)*clhs4 + clhs227 + clhs49;
            lhs(4,1)=DN(1,0)*clhs22 + DN(1,1)*clhs24 + DN(1,2)*clhs27 + clhs139;
            lhs(4,2)=DN(1,0)*clhs30 + DN(1,1)*clhs32 + DN(1,2)*clhs34 + clhs184;
            lhs(4,3)=clhs228*clhs52 + clhs68*clhs71 + clhs69*clhs70 - clhs70;
            lhs(4,4)=DN(1,0)*clhs44 + DN(1,1)*clhs46 + DN(1,2)*clhs48 + clhs10*clhs229 + clhs231;
            lhs(4,5)=DN(1,0)*clhs55 + DN(1,1)*clhs57 + DN(1,2)*clhs60 + clhs233;
            lhs(4,6)=DN(1,0)*clhs62 + DN(1,1)*clhs64 + DN(1,2)*clhs66 + clhs234;
            lhs(4,7)=DN(1,0)*clhs237;
            lhs(4,8)=DN(1,0)*clhs73 + DN(1,1)*clhs75 + DN(1,2)*clhs77 + clhs238 + clhs241;
            lhs(4,9)=DN(1,0)*clhs84 + DN(1,1)*clhs86 + DN(1,2)*clhs89 + clhs242;
            lhs(4,10)=DN(1,0)*clhs91 + DN(1,1)*clhs93 + DN(1,2)*clhs95 + clhs243;
            lhs(4,11)=clhs244*clhs69 - clhs244 + clhs245*clhs71 + clhs52*clhs99;
            lhs(4,12)=DN(1,0)*clhs100 + DN(1,1)*clhs102 + DN(1,2)*clhs104 + clhs246 + clhs248;
            lhs(4,13)=DN(1,0)*clhs111 + DN(1,1)*clhs113 + DN(1,2)*clhs116 + clhs249;
            lhs(4,14)=DN(1,0)*clhs118 + DN(1,1)*clhs120 + DN(1,2)*clhs122 + clhs250;
            lhs(4,15)=DN(3,0)*clhs236 + clhs251*clhs69 - clhs251 + clhs252*clhs71;
            lhs(5,0)=DN(1,0)*clhs2 + DN(1,1)*clhs126 + DN(1,2)*clhs127 + clhs61;
            lhs(5,1)=DN(1,0)*clhs24 + DN(1,1)*clhs128 + DN(1,2)*clhs130 + clhs143 + clhs227;
            lhs(5,2)=DN(1,0)*clhs32 + DN(1,1)*clhs132 + DN(1,2)*clhs134 + clhs186;
            lhs(5,3)=clhs148*clhs71 + clhs149*clhs69 - clhs149 + clhs253*clhs52;
            lhs(5,4)=DN(1,0)*clhs46 + DN(1,1)*clhs137 + DN(1,2)*clhs138 + clhs233;
            lhs(5,5)=DN(1,0)*clhs57 + DN(1,1)*clhs140 + DN(1,2)*clhs142 + clhs10*clhs254 + clhs231;
            lhs(5,6)=DN(1,0)*clhs64 + DN(1,1)*clhs144 + DN(1,2)*clhs146 + clhs256;
            lhs(5,7)=DN(1,1)*clhs237;
            lhs(5,8)=DN(1,0)*clhs75 + DN(1,1)*clhs151 + DN(1,2)*clhs152 + clhs257;
            lhs(5,9)=DN(1,0)*clhs86 + DN(1,1)*clhs154 + DN(1,2)*clhs156 + clhs241 + clhs258;
            lhs(5,10)=DN(1,0)*clhs93 + DN(1,1)*clhs158 + DN(1,2)*clhs160 + clhs259;
            lhs(5,11)=clhs164*clhs52 + clhs260*clhs69 - clhs260 + clhs261*clhs71;
            lhs(5,12)=DN(1,0)*clhs102 + DN(1,1)*clhs165 + DN(1,2)*clhs166 + clhs262;
            lhs(5,13)=DN(1,0)*clhs113 + DN(1,1)*clhs168 + DN(1,2)*clhs170 + clhs248 + clhs263;
            lhs(5,14)=DN(1,0)*clhs120 + DN(1,1)*clhs172 + DN(1,2)*clhs174 + clhs264;
            lhs(5,15)=DN(3,1)*clhs236 + clhs265*clhs69 - clhs265 + clhs266*clhs71;
            lhs(6,0)=DN(1,0)*clhs4 + DN(1,1)*clhs127 + DN(1,2)*clhs178 + clhs67;
            lhs(6,1)=DN(1,0)*clhs27 + DN(1,1)*clhs130 + DN(1,2)*clhs179 + clhs147;
            lhs(6,2)=DN(1,0)*clhs34 + DN(1,1)*clhs134 + DN(1,2)*clhs180 + clhs188 + clhs227;
            lhs(6,3)=clhs189*clhs71 + clhs190*clhs69 - clhs190 + clhs267*clhs52;
            lhs(6,4)=DN(1,0)*clhs48 + DN(1,1)*clhs138 + DN(1,2)*clhs182 + clhs234;
            lhs(6,5)=DN(1,0)*clhs60 + DN(1,1)*clhs142 + DN(1,2)*clhs185 + clhs256;
            lhs(6,6)=DN(1,0)*clhs66 + DN(1,1)*clhs146 + DN(1,2)*clhs187 + clhs10*clhs268 + clhs231;
            lhs(6,7)=DN(1,2)*clhs237;
            lhs(6,8)=DN(1,0)*clhs77 + DN(1,1)*clhs152 + DN(1,2)*clhs192 + clhs270;
            lhs(6,9)=DN(1,0)*clhs89 + DN(1,1)*clhs156 + DN(1,2)*clhs194 + clhs271;
            lhs(6,10)=DN(1,0)*clhs95 + DN(1,1)*clhs160 + DN(1,2)*clhs196 + clhs241 + clhs272;
            lhs(6,11)=clhs200*clhs52 + clhs273*clhs69 - clhs273 + clhs274*clhs71;
            lhs(6,12)=DN(1,0)*clhs104 + DN(1,1)*clhs166 + DN(1,2)*clhs201 + clhs275;
            lhs(6,13)=DN(1,0)*clhs116 + DN(1,1)*clhs170 + DN(1,2)*clhs203 + clhs276;
            lhs(6,14)=DN(1,0)*clhs122 + DN(1,1)*clhs174 + DN(1,2)*clhs205 + clhs248 + clhs277;
            lhs(6,15)=DN(3,2)*clhs236 + clhs278*clhs69 - clhs278 + clhs279*clhs71;
            lhs(7,0)=DN(1,0)*clhs209;
            lhs(7,1)=DN(1,1)*clhs209;
            lhs(7,2)=DN(1,2)*clhs209;
            lhs(7,3)=clhs217;
            lhs(7,4)=DN(1,0)*clhs213;
            lhs(7,5)=DN(1,1)*clhs213;
            lhs(7,6)=DN(1,2)*clhs213;
            lhs(7,7)=clhs210*clhs230 + clhs211*clhs229 + clhs211*clhs254 + clhs211*clhs268;
            lhs(7,8)=DN(1,0)*clhs219;
            lhs(7,9)=DN(1,1)*clhs219;
            lhs(7,10)=DN(1,2)*clhs219;
            lhs(7,11)=clhs283;
            lhs(7,12)=DN(1,0)*clhs222;
            lhs(7,13)=DN(1,1)*clhs222;
            lhs(7,14)=DN(1,2)*clhs222;
            lhs(7,15)=clhs284;
            lhs(8,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + DN(2,2)*clhs4 + clhs288 + clhs78;
            lhs(8,1)=DN(2,0)*clhs22 + DN(2,1)*clhs24 + DN(2,2)*clhs27 + clhs153;
            lhs(8,2)=DN(2,0)*clhs30 + DN(2,1)*clhs32 + DN(2,2)*clhs34 + clhs193;
            lhs(8,3)=clhs228*clhs80 + clhs69*clhs98 + clhs71*clhs97 - clhs98;
            lhs(8,4)=DN(2,0)*clhs44 + DN(2,1)*clhs46 + DN(2,2)*clhs48 + clhs238 + clhs289;
            lhs(8,5)=DN(2,0)*clhs55 + DN(2,1)*clhs57 + DN(2,2)*clhs60 + clhs257;
            lhs(8,6)=DN(2,0)*clhs62 + DN(2,1)*clhs64 + DN(2,2)*clhs66 + clhs270;
            lhs(8,7)=clhs244*clhs71 + clhs245*clhs69 - clhs245 + clhs72*clhs80;
            lhs(8,8)=DN(2,0)*clhs73 + DN(2,1)*clhs75 + DN(2,2)*clhs77 + clhs10*clhs290 + clhs292;
            lhs(8,9)=DN(2,0)*clhs84 + DN(2,1)*clhs86 + DN(2,2)*clhs89 + clhs294;
            lhs(8,10)=DN(2,0)*clhs91 + DN(2,1)*clhs93 + DN(2,2)*clhs95 + clhs295;
            lhs(8,11)=DN(2,0)*clhs298;
            lhs(8,12)=DN(2,0)*clhs100 + DN(2,1)*clhs102 + DN(2,2)*clhs104 + clhs299 + clhs302;
            lhs(8,13)=DN(2,0)*clhs111 + DN(2,1)*clhs113 + DN(2,2)*clhs116 + clhs303;
            lhs(8,14)=DN(2,0)*clhs118 + DN(2,1)*clhs120 + DN(2,2)*clhs122 + clhs304;
            lhs(8,15)=DN(3,0)*clhs297 + clhs305*clhs69 - clhs305 + clhs306*clhs71;
            lhs(9,0)=DN(2,0)*clhs2 + DN(2,1)*clhs126 + DN(2,2)*clhs127 + clhs90;
            lhs(9,1)=DN(2,0)*clhs24 + DN(2,1)*clhs128 + DN(2,2)*clhs130 + clhs157 + clhs288;
            lhs(9,2)=DN(2,0)*clhs32 + DN(2,1)*clhs132 + DN(2,2)*clhs134 + clhs195;
            lhs(9,3)=clhs162*clhs71 + clhs163*clhs69 - clhs163 + clhs253*clhs80;
            lhs(9,4)=DN(2,0)*clhs46 + DN(2,1)*clhs137 + DN(2,2)*clhs138 + clhs242;
            lhs(9,5)=DN(2,0)*clhs57 + DN(2,1)*clhs140 + DN(2,2)*clhs142 + clhs258 + clhs289;
            lhs(9,6)=DN(2,0)*clhs64 + DN(2,1)*clhs144 + DN(2,2)*clhs146 + clhs271;
            lhs(9,7)=clhs150*clhs80 + clhs260*clhs71 + clhs261*clhs69 - clhs261;
            lhs(9,8)=DN(2,0)*clhs75 + DN(2,1)*clhs151 + DN(2,2)*clhs152 + clhs294;
            lhs(9,9)=DN(2,0)*clhs86 + DN(2,1)*clhs154 + DN(2,2)*clhs156 + clhs10*clhs307 + clhs292;
            lhs(9,10)=DN(2,0)*clhs93 + DN(2,1)*clhs158 + DN(2,2)*clhs160 + clhs309;
            lhs(9,11)=DN(2,1)*clhs298;
            lhs(9,12)=DN(2,0)*clhs102 + DN(2,1)*clhs165 + DN(2,2)*clhs166 + clhs310;
            lhs(9,13)=DN(2,0)*clhs113 + DN(2,1)*clhs168 + DN(2,2)*clhs170 + clhs302 + clhs311;
            lhs(9,14)=DN(2,0)*clhs120 + DN(2,1)*clhs172 + DN(2,2)*clhs174 + clhs312;
            lhs(9,15)=DN(3,1)*clhs297 + clhs313*clhs69 - clhs313 + clhs314*clhs71;
            lhs(10,0)=DN(2,0)*clhs4 + DN(2,1)*clhs127 + DN(2,2)*clhs178 + clhs96;
            lhs(10,1)=DN(2,0)*clhs27 + DN(2,1)*clhs130 + DN(2,2)*clhs179 + clhs161;
            lhs(10,2)=DN(2,0)*clhs34 + DN(2,1)*clhs134 + DN(2,2)*clhs180 + clhs197 + clhs288;
            lhs(10,3)=clhs198*clhs71 + clhs199*clhs69 - clhs199 + clhs267*clhs80;
            lhs(10,4)=DN(2,0)*clhs48 + DN(2,1)*clhs138 + DN(2,2)*clhs182 + clhs243;
            lhs(10,5)=DN(2,0)*clhs60 + DN(2,1)*clhs142 + DN(2,2)*clhs185 + clhs259;
            lhs(10,6)=DN(2,0)*clhs66 + DN(2,1)*clhs146 + DN(2,2)*clhs187 + clhs272 + clhs289;
            lhs(10,7)=clhs191*clhs80 + clhs273*clhs71 + clhs274*clhs69 - clhs274;
            lhs(10,8)=DN(2,0)*clhs77 + DN(2,1)*clhs152 + DN(2,2)*clhs192 + clhs295;
            lhs(10,9)=DN(2,0)*clhs89 + DN(2,1)*clhs156 + DN(2,2)*clhs194 + clhs309;
            lhs(10,10)=DN(2,0)*clhs95 + DN(2,1)*clhs160 + DN(2,2)*clhs196 + clhs10*clhs315 + clhs292;
            lhs(10,11)=DN(2,2)*clhs298;
            lhs(10,12)=DN(2,0)*clhs104 + DN(2,1)*clhs166 + DN(2,2)*clhs201 + clhs317;
            lhs(10,13)=DN(2,0)*clhs116 + DN(2,1)*clhs170 + DN(2,2)*clhs203 + clhs318;
            lhs(10,14)=DN(2,0)*clhs122 + DN(2,1)*clhs174 + DN(2,2)*clhs205 + clhs302 + clhs319;
            lhs(10,15)=DN(3,2)*clhs297 + clhs320*clhs69 - clhs320 + clhs321*clhs71;
            lhs(11,0)=DN(2,0)*clhs209;
            lhs(11,1)=DN(2,1)*clhs209;
            lhs(11,2)=DN(2,2)*clhs209;
            lhs(11,3)=clhs220;
            lhs(11,4)=DN(2,0)*clhs213;
            lhs(11,5)=DN(2,1)*clhs213;
            lhs(11,6)=DN(2,2)*clhs213;
            lhs(11,7)=clhs283;
            lhs(11,8)=DN(2,0)*clhs219;
            lhs(11,9)=DN(2,1)*clhs219;
            lhs(11,10)=DN(2,2)*clhs219;
            lhs(11,11)=clhs210*clhs291 + clhs211*clhs290 + clhs211*clhs307 + clhs211*clhs315;
            lhs(11,12)=DN(2,0)*clhs222;
            lhs(11,13)=DN(2,1)*clhs222;
            lhs(11,14)=DN(2,2)*clhs222;
            lhs(11,15)=clhs322;
            lhs(12,0)=DN(3,0)*clhs0 + DN(3,1)*clhs2 + DN(3,2)*clhs4 + clhs105 + clhs326;
            lhs(12,1)=DN(3,0)*clhs22 + DN(3,1)*clhs24 + DN(3,2)*clhs27 + clhs167;
            lhs(12,2)=DN(3,0)*clhs30 + DN(3,1)*clhs32 + DN(3,2)*clhs34 + clhs202;
            lhs(12,3)=clhs107*clhs228 + clhs124*clhs71 + clhs125*clhs69 - clhs125;
            lhs(12,4)=DN(3,0)*clhs44 + DN(3,1)*clhs46 + DN(3,2)*clhs48 + clhs246 + clhs327;
            lhs(12,5)=DN(3,0)*clhs55 + DN(3,1)*clhs57 + DN(3,2)*clhs60 + clhs262;
            lhs(12,6)=DN(3,0)*clhs62 + DN(3,1)*clhs64 + DN(3,2)*clhs66 + clhs275;
            lhs(12,7)=clhs107*clhs72 + clhs251*clhs71 + clhs252*clhs69 - clhs252;
            lhs(12,8)=DN(3,0)*clhs73 + DN(3,1)*clhs75 + DN(3,2)*clhs77 + clhs299 + clhs328;
            lhs(12,9)=DN(3,0)*clhs84 + DN(3,1)*clhs86 + DN(3,2)*clhs89 + clhs310;
            lhs(12,10)=DN(3,0)*clhs91 + DN(3,1)*clhs93 + DN(3,2)*clhs95 + clhs317;
            lhs(12,11)=clhs107*clhs99 + clhs305*clhs71 + clhs306*clhs69 - clhs306;
            lhs(12,12)=DN(3,0)*clhs100 + DN(3,1)*clhs102 + DN(3,2)*clhs104 + clhs10*clhs329 + clhs331;
            lhs(12,13)=DN(3,0)*clhs111 + DN(3,1)*clhs113 + DN(3,2)*clhs116 + clhs333;
            lhs(12,14)=DN(3,0)*clhs118 + DN(3,1)*clhs120 + DN(3,2)*clhs122 + clhs334;
            lhs(12,15)=DN(3,0)*clhs335;
            lhs(13,0)=DN(3,0)*clhs2 + DN(3,1)*clhs126 + DN(3,2)*clhs127 + clhs117;
            lhs(13,1)=DN(3,0)*clhs24 + DN(3,1)*clhs128 + DN(3,2)*clhs130 + clhs171 + clhs326;
            lhs(13,2)=DN(3,0)*clhs32 + DN(3,1)*clhs132 + DN(3,2)*clhs134 + clhs204;
            lhs(13,3)=clhs107*clhs253 + clhs176*clhs71 + clhs177*clhs69 - clhs177;
            lhs(13,4)=DN(3,0)*clhs46 + DN(3,1)*clhs137 + DN(3,2)*clhs138 + clhs249;
            lhs(13,5)=DN(3,0)*clhs57 + DN(3,1)*clhs140 + DN(3,2)*clhs142 + clhs263 + clhs327;
            lhs(13,6)=DN(3,0)*clhs64 + DN(3,1)*clhs144 + DN(3,2)*clhs146 + clhs276;
            lhs(13,7)=clhs107*clhs150 + clhs265*clhs71 + clhs266*clhs69 - clhs266;
            lhs(13,8)=DN(3,0)*clhs75 + DN(3,1)*clhs151 + DN(3,2)*clhs152 + clhs303;
            lhs(13,9)=DN(3,0)*clhs86 + DN(3,1)*clhs154 + DN(3,2)*clhs156 + clhs311 + clhs328;
            lhs(13,10)=DN(3,0)*clhs93 + DN(3,1)*clhs158 + DN(3,2)*clhs160 + clhs318;
            lhs(13,11)=clhs107*clhs164 + clhs313*clhs71 + clhs314*clhs69 - clhs314;
            lhs(13,12)=DN(3,0)*clhs102 + DN(3,1)*clhs165 + DN(3,2)*clhs166 + clhs333;
            lhs(13,13)=DN(3,0)*clhs113 + DN(3,1)*clhs168 + DN(3,2)*clhs170 + clhs10*clhs336 + clhs331;
            lhs(13,14)=DN(3,0)*clhs120 + DN(3,1)*clhs172 + DN(3,2)*clhs174 + clhs337;
            lhs(13,15)=DN(3,1)*clhs335;
            lhs(14,0)=DN(3,0)*clhs4 + DN(3,1)*clhs127 + DN(3,2)*clhs178 + clhs123;
            lhs(14,1)=DN(3,0)*clhs27 + DN(3,1)*clhs130 + DN(3,2)*clhs179 + clhs175;
            lhs(14,2)=DN(3,0)*clhs34 + DN(3,1)*clhs134 + DN(3,2)*clhs180 + clhs206 + clhs326;
            lhs(14,3)=clhs107*clhs267 + clhs207*clhs71 + clhs208*clhs69 - clhs208;
            lhs(14,4)=DN(3,0)*clhs48 + DN(3,1)*clhs138 + DN(3,2)*clhs182 + clhs250;
            lhs(14,5)=DN(3,0)*clhs60 + DN(3,1)*clhs142 + DN(3,2)*clhs185 + clhs264;
            lhs(14,6)=DN(3,0)*clhs66 + DN(3,1)*clhs146 + DN(3,2)*clhs187 + clhs277 + clhs327;
            lhs(14,7)=clhs107*clhs191 + clhs278*clhs71 + clhs279*clhs69 - clhs279;
            lhs(14,8)=DN(3,0)*clhs77 + DN(3,1)*clhs152 + DN(3,2)*clhs192 + clhs304;
            lhs(14,9)=DN(3,0)*clhs89 + DN(3,1)*clhs156 + DN(3,2)*clhs194 + clhs312;
            lhs(14,10)=DN(3,0)*clhs95 + DN(3,1)*clhs160 + DN(3,2)*clhs196 + clhs319 + clhs328;
            lhs(14,11)=clhs107*clhs200 + clhs320*clhs71 + clhs321*clhs69 - clhs321;
            lhs(14,12)=DN(3,0)*clhs104 + DN(3,1)*clhs166 + DN(3,2)*clhs201 + clhs334;
            lhs(14,13)=DN(3,0)*clhs116 + DN(3,1)*clhs170 + DN(3,2)*clhs203 + clhs337;
            lhs(14,14)=DN(3,0)*clhs122 + DN(3,1)*clhs174 + DN(3,2)*clhs205 + clhs10*clhs338 + clhs331;
            lhs(14,15)=DN(3,2)*clhs335;
            lhs(15,0)=DN(3,0)*clhs209;
            lhs(15,1)=DN(3,1)*clhs209;
            lhs(15,2)=DN(3,2)*clhs209;
            lhs(15,3)=clhs223;
            lhs(15,4)=DN(3,0)*clhs213;
            lhs(15,5)=DN(3,1)*clhs213;
            lhs(15,6)=DN(3,2)*clhs213;
            lhs(15,7)=clhs284;
            lhs(15,8)=DN(3,0)*clhs219;
            lhs(15,9)=DN(3,1)*clhs219;
            lhs(15,10)=DN(3,2)*clhs219;
            lhs(15,11)=clhs322;
            lhs(15,12)=DN(3,0)*clhs222;
            lhs(15,13)=DN(3,1)*clhs222;
            lhs(15,14)=DN(3,2)*clhs222;
            lhs(15,15)=clhs210*clhs330 + clhs211*clhs329 + clhs211*clhs336 + clhs211*clhs338;


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
const double clhs24 =             -N[0];
const double clhs25 =             pow(c, -2);
const double clhs26 =             1.0/rho;
const double clhs27 =             N[0]*bdf0*clhs25*clhs26;
const double clhs28 =             1.0*clhs14*clhs15;
const double clhs29 =             1.0*clhs15*rho;
const double clhs30 =             clhs11*clhs29;
const double clhs31 =             clhs10*clhs28 + clhs24 + clhs27*clhs7 + clhs30;
const double clhs32 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs33 =             C(0,2)*DN(1,0);
const double clhs34 =             C(2,2)*DN(1,1) + clhs33;
const double clhs35 =             DN(1,0)*clhs22;
const double clhs36 =             N[0]*bdf0*rho;
const double clhs37 =             N[1]*clhs36;
const double clhs38 =             DN(1,0)*clhs4 + DN(1,1)*clhs5;
const double clhs39 =             N[1]*bdf0;
const double clhs40 =             clhs38 + clhs39;
const double clhs41 =             clhs10*clhs38 + clhs16*clhs40 + clhs17*clhs40 + clhs37;
const double clhs42 =             C(0,1)*DN(1,1) + clhs33;
const double clhs43 =             C(1,2)*DN(1,1);
const double clhs44 =             C(2,2)*DN(1,0) + clhs43;
const double clhs45 =             DN(1,1)*clhs22;
const double clhs46 =             DN(0,0)*N[1];
const double clhs47 =             bdf0*clhs25*clhs26*clhs7;
const double clhs48 =             DN(1,0)*N[0];
const double clhs49 =             1.0*clhs14*clhs15*rho;
const double clhs50 =             1.0*DN(1,0)*clhs15*rho;
const double clhs51 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs52 =             C(0,2)*DN(2,0);
const double clhs53 =             C(2,2)*DN(2,1) + clhs52;
const double clhs54 =             DN(2,0)*clhs22;
const double clhs55 =             N[2]*clhs36;
const double clhs56 =             DN(2,0)*clhs4 + DN(2,1)*clhs5;
const double clhs57 =             N[2]*bdf0;
const double clhs58 =             clhs56 + clhs57;
const double clhs59 =             clhs10*clhs56 + clhs16*clhs58 + clhs17*clhs58 + clhs55;
const double clhs60 =             C(0,1)*DN(2,1) + clhs52;
const double clhs61 =             C(1,2)*DN(2,1);
const double clhs62 =             C(2,2)*DN(2,0) + clhs61;
const double clhs63 =             DN(2,1)*clhs22;
const double clhs64 =             DN(0,0)*N[2];
const double clhs65 =             DN(2,0)*N[0];
const double clhs66 =             C(0,1)*DN(0,0) + clhs20;
const double clhs67 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs68 =             pow(DN(0,1), 2);
const double clhs69 =             C(0,1)*DN(1,0) + clhs43;
const double clhs70 =             DN(0,1)*clhs7;
const double clhs71 =             DN(1,0)*clhs70;
const double clhs72 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs73 =             DN(1,1)*clhs70;
const double clhs74 =             DN(0,1)*N[1];
const double clhs75 =             DN(1,1)*N[0];
const double clhs76 =             1.0*DN(1,1)*clhs15*rho;
const double clhs77 =             C(0,1)*DN(2,0) + clhs61;
const double clhs78 =             DN(2,0)*clhs70;
const double clhs79 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs80 =             DN(2,1)*clhs70;
const double clhs81 =             DN(0,1)*N[2];
const double clhs82 =             DN(2,1)*N[0];
const double clhs83 =             clhs12*clhs29 + clhs24;
const double clhs84 =             bdf0*clhs25*clhs26;
const double clhs85 =             1.0*clhs15;
const double clhs86 =             -N[1];
const double clhs87 =             clhs29*clhs40 + clhs86;
const double clhs88 =             1.0*DN(0,0)*clhs15;
const double clhs89 =             1.0*DN(0,1)*clhs15;
const double clhs90 =             DN(1,0)*clhs88 + DN(1,1)*clhs89 + N[1]*clhs27;
const double clhs91 =             -N[2];
const double clhs92 =             clhs29*clhs58 + clhs91;
const double clhs93 =             DN(2,0)*clhs88 + DN(2,1)*clhs89 + N[2]*clhs27;
const double clhs94 =             N[1]*rho;
const double clhs95 =             1.0*N[1]*clhs13*clhs14*clhs15;
const double clhs96 =             1.0*clhs13*clhs15*clhs38;
const double clhs97 =             clhs11*clhs94 + clhs12*clhs95 + clhs12*clhs96 + clhs37;
const double clhs98 =             1.0*DN(0,0)*clhs15*rho;
const double clhs99 =             pow(DN(1,0), 2);
const double clhs100 =             pow(N[1], 2);
const double clhs101 =             clhs100*clhs9 + clhs38*clhs94 + clhs40*clhs95 + clhs40*clhs96;
const double clhs102 =             DN(1,0)*clhs7;
const double clhs103 =             DN(1,1)*clhs102;
const double clhs104 =             clhs25*clhs26*clhs7;
const double clhs105 =             clhs29*clhs38;
const double clhs106 =             clhs104*clhs39 + clhs105 + clhs28*clhs94 + clhs86;
const double clhs107 =             DN(2,0)*clhs102;
const double clhs108 =             N[1]*N[2]*bdf0;
const double clhs109 =             clhs108*rho;
const double clhs110 =             clhs109 + clhs56*clhs94 + clhs58*clhs95 + clhs58*clhs96;
const double clhs111 =             DN(2,1)*clhs102;
const double clhs112 =             DN(1,0)*N[2];
const double clhs113 =             DN(2,0)*N[1];
const double clhs114 =             1.0*DN(0,1)*clhs15*rho;
const double clhs115 =             pow(DN(1,1), 2);
const double clhs116 =             DN(1,1)*clhs7;
const double clhs117 =             DN(2,0)*clhs116;
const double clhs118 =             DN(2,1)*clhs116;
const double clhs119 =             DN(1,1)*N[2];
const double clhs120 =             DN(2,1)*N[1];
const double clhs121 =             1.0*DN(1,0)*DN(2,0)*clhs15 + 1.0*DN(1,1)*DN(2,1)*clhs15 + clhs108*clhs25*clhs26;
const double clhs122 =             N[2]*rho;
const double clhs123 =             1.0*N[2]*clhs13*clhs14*clhs15;
const double clhs124 =             1.0*clhs13*clhs15*clhs56;
const double clhs125 =             clhs11*clhs122 + clhs12*clhs123 + clhs12*clhs124 + clhs55;
const double clhs126 =             clhs109 + clhs122*clhs38 + clhs123*clhs40 + clhs124*clhs40;
const double clhs127 =             pow(DN(2,0), 2);
const double clhs128 =             pow(N[2], 2);
const double clhs129 =             clhs122*clhs56 + clhs123*clhs58 + clhs124*clhs58 + clhs128*clhs9;
const double clhs130 =             DN(2,0)*DN(2,1)*clhs7;
const double clhs131 =             clhs104*clhs57 + clhs122*clhs28 + clhs29*clhs56 + clhs91;
const double clhs132 =             pow(DN(2,1), 2);
            lhs(0,0)=DN(0,0)*clhs0 + DN(0,1)*clhs2 + clhs18 + clhs3*clhs7;
            lhs(0,1)=DN(0,0)*clhs19 + DN(0,1)*clhs21 + clhs23;
            lhs(0,2)=DN(0,0)*clhs31;
            lhs(0,3)=DN(0,0)*clhs32 + DN(0,1)*clhs34 + clhs35 + clhs41;
            lhs(0,4)=DN(0,0)*clhs42 + DN(0,1)*clhs44 + clhs45;
            lhs(0,5)=clhs11*clhs50 + clhs46*clhs47 - clhs46 + clhs48*clhs49;
            lhs(0,6)=DN(0,0)*clhs51 + DN(0,1)*clhs53 + clhs54 + clhs59;
            lhs(0,7)=DN(0,0)*clhs60 + DN(0,1)*clhs62 + clhs63;
            lhs(0,8)=DN(2,0)*clhs30 + clhs47*clhs64 + clhs49*clhs65 - clhs64;
            lhs(1,0)=DN(0,0)*clhs2 + DN(0,1)*clhs66 + clhs23;
            lhs(1,1)=DN(0,0)*clhs21 + DN(0,1)*clhs67 + clhs18 + clhs68*clhs7;
            lhs(1,2)=DN(0,1)*clhs31;
            lhs(1,3)=DN(0,0)*clhs34 + DN(0,1)*clhs69 + clhs71;
            lhs(1,4)=DN(0,0)*clhs44 + DN(0,1)*clhs72 + clhs41 + clhs73;
            lhs(1,5)=clhs11*clhs76 + clhs47*clhs74 + clhs49*clhs75 - clhs74;
            lhs(1,6)=DN(0,0)*clhs53 + DN(0,1)*clhs77 + clhs78;
            lhs(1,7)=DN(0,0)*clhs62 + DN(0,1)*clhs79 + clhs59 + clhs80;
            lhs(1,8)=DN(2,1)*clhs30 + clhs47*clhs81 + clhs49*clhs82 - clhs81;
            lhs(2,0)=DN(0,0)*clhs83;
            lhs(2,1)=DN(0,1)*clhs83;
            lhs(2,2)=clhs3*clhs85 + clhs68*clhs85 + clhs8*clhs84;
            lhs(2,3)=DN(0,0)*clhs87;
            lhs(2,4)=DN(0,1)*clhs87;
            lhs(2,5)=clhs90;
            lhs(2,6)=DN(0,0)*clhs92;
            lhs(2,7)=DN(0,1)*clhs92;
            lhs(2,8)=clhs93;
            lhs(3,0)=DN(1,0)*clhs0 + DN(1,1)*clhs2 + clhs35 + clhs97;
            lhs(3,1)=DN(1,0)*clhs19 + DN(1,1)*clhs21 + clhs71;
            lhs(3,2)=clhs38*clhs98 + clhs46*clhs49 + clhs47*clhs48 - clhs48;
            lhs(3,3)=DN(1,0)*clhs32 + DN(1,1)*clhs34 + clhs101 + clhs7*clhs99;
            lhs(3,4)=DN(1,0)*clhs42 + DN(1,1)*clhs44 + clhs103;
            lhs(3,5)=DN(1,0)*clhs106;
            lhs(3,6)=DN(1,0)*clhs51 + DN(1,1)*clhs53 + clhs107 + clhs110;
            lhs(3,7)=DN(1,0)*clhs60 + DN(1,1)*clhs62 + clhs111;
            lhs(3,8)=DN(2,0)*clhs105 + clhs112*clhs47 - clhs112 + clhs113*clhs49;
            lhs(4,0)=DN(1,0)*clhs2 + DN(1,1)*clhs66 + clhs45;
            lhs(4,1)=DN(1,0)*clhs21 + DN(1,1)*clhs67 + clhs73 + clhs97;
            lhs(4,2)=clhs114*clhs38 + clhs47*clhs75 + clhs49*clhs74 - clhs75;
            lhs(4,3)=DN(1,0)*clhs34 + DN(1,1)*clhs69 + clhs103;
            lhs(4,4)=DN(1,0)*clhs44 + DN(1,1)*clhs72 + clhs101 + clhs115*clhs7;
            lhs(4,5)=DN(1,1)*clhs106;
            lhs(4,6)=DN(1,0)*clhs53 + DN(1,1)*clhs77 + clhs117;
            lhs(4,7)=DN(1,0)*clhs62 + DN(1,1)*clhs79 + clhs110 + clhs118;
            lhs(4,8)=DN(2,1)*clhs105 + clhs119*clhs47 - clhs119 + clhs120*clhs49;
            lhs(5,0)=DN(1,0)*clhs83;
            lhs(5,1)=DN(1,1)*clhs83;
            lhs(5,2)=clhs90;
            lhs(5,3)=DN(1,0)*clhs87;
            lhs(5,4)=DN(1,1)*clhs87;
            lhs(5,5)=clhs100*clhs84 + clhs115*clhs85 + clhs85*clhs99;
            lhs(5,6)=DN(1,0)*clhs92;
            lhs(5,7)=DN(1,1)*clhs92;
            lhs(5,8)=clhs121;
            lhs(6,0)=DN(2,0)*clhs0 + DN(2,1)*clhs2 + clhs125 + clhs54;
            lhs(6,1)=DN(2,0)*clhs19 + DN(2,1)*clhs21 + clhs78;
            lhs(6,2)=clhs47*clhs65 + clhs49*clhs64 + clhs56*clhs98 - clhs65;
            lhs(6,3)=DN(2,0)*clhs32 + DN(2,1)*clhs34 + clhs107 + clhs126;
            lhs(6,4)=DN(2,0)*clhs42 + DN(2,1)*clhs44 + clhs117;
            lhs(6,5)=clhs112*clhs49 + clhs113*clhs47 - clhs113 + clhs50*clhs56;
            lhs(6,6)=DN(2,0)*clhs51 + DN(2,1)*clhs53 + clhs127*clhs7 + clhs129;
            lhs(6,7)=DN(2,0)*clhs60 + DN(2,1)*clhs62 + clhs130;
            lhs(6,8)=DN(2,0)*clhs131;
            lhs(7,0)=DN(2,0)*clhs2 + DN(2,1)*clhs66 + clhs63;
            lhs(7,1)=DN(2,0)*clhs21 + DN(2,1)*clhs67 + clhs125 + clhs80;
            lhs(7,2)=clhs114*clhs56 + clhs47*clhs82 + clhs49*clhs81 - clhs82;
            lhs(7,3)=DN(2,0)*clhs34 + DN(2,1)*clhs69 + clhs111;
            lhs(7,4)=DN(2,0)*clhs44 + DN(2,1)*clhs72 + clhs118 + clhs126;
            lhs(7,5)=clhs119*clhs49 + clhs120*clhs47 - clhs120 + clhs56*clhs76;
            lhs(7,6)=DN(2,0)*clhs53 + DN(2,1)*clhs77 + clhs130;
            lhs(7,7)=DN(2,0)*clhs62 + DN(2,1)*clhs79 + clhs129 + clhs132*clhs7;
            lhs(7,8)=DN(2,1)*clhs131;
            lhs(8,0)=DN(2,0)*clhs83;
            lhs(8,1)=DN(2,1)*clhs83;
            lhs(8,2)=clhs93;
            lhs(8,3)=DN(2,0)*clhs87;
            lhs(8,4)=DN(2,1)*clhs87;
            lhs(8,5)=clhs121;
            lhs(8,6)=DN(2,0)*clhs92;
            lhs(8,7)=DN(2,1)*clhs92;
            lhs(8,8)=clhs127*clhs85 + clhs128*clhs84 + clhs132*clhs85;


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
const double crhs17 =             (N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]) + N[3]*(bdf0*p[3] + bdf1*pn[3] + bdf2*pnn[3]))/(pow(c, 2)*rho);
const double crhs18 =             (crhs11*h/stab_c1 + mu)*(crhs12 + crhs13 + crhs14 + crhs15 + crhs16 + crhs17 + crhs3 + crhs4 + crhs5 + crhs6);
const double crhs19 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
const double crhs20 =             N[0]*crhs19*rho;
const double crhs21 =             1.0/(crhs11/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs22 =             1.0*crhs21*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + crhs10 + crhs2);
const double crhs23 =             rho*(DN(0,0)*crhs7 + DN(0,1)*crhs8 + DN(0,2)*crhs9);
const double crhs24 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs25 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)));
const double crhs26 =             rho*(crhs7*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs8*(crhs13 + crhs14 + crhs15 + crhs16) + crhs9*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
const double crhs27 =             1.0*crhs21*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs24 + crhs25 + crhs26);
const double crhs28 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs29 =             rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)));
const double crhs30 =             rho*(crhs12*crhs9 + crhs7*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs8*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)));
const double crhs31 =             1.0*crhs21*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs28 + crhs29 + crhs30);
const double crhs32 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0);
const double crhs33 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
const double crhs34 =             N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2);
const double crhs35 =             N[1]*crhs19*rho;
const double crhs36 =             rho*(DN(1,0)*crhs7 + DN(1,1)*crhs8 + DN(1,2)*crhs9);
const double crhs37 =             N[2]*crhs19*rho;
const double crhs38 =             rho*(DN(2,0)*crhs7 + DN(2,1)*crhs8 + DN(2,2)*crhs9);
const double crhs39 =             N[3]*crhs19*rho;
const double crhs40 =             rho*(DN(3,0)*crhs7 + DN(3,1)*crhs8 + DN(3,2)*crhs9);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs18 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - N[0]*crhs10 - N[0]*crhs2 - crhs20*crhs22 - crhs22*crhs23;
            rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs18 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs24 - N[0]*crhs25 - N[0]*crhs26 - crhs20*crhs27 - crhs23*crhs27;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs18 - DN(0,2)*stress[2] + N[0]*crhs28 - N[0]*crhs29 - N[0]*crhs30 - crhs20*crhs31 - crhs23*crhs31;
            rhs[3]=-DN(0,0)*crhs22 + DN(0,0)*crhs32 - DN(0,1)*crhs27 + DN(0,1)*crhs33 - DN(0,2)*crhs31 + DN(0,2)*crhs34 - N[0]*crhs17;
            rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs18 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - N[1]*crhs10 - N[1]*crhs2 - crhs22*crhs35 - crhs22*crhs36;
            rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs18 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs24 - N[1]*crhs25 - N[1]*crhs26 - crhs27*crhs35 - crhs27*crhs36;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs18 - DN(1,2)*stress[2] + N[1]*crhs28 - N[1]*crhs29 - N[1]*crhs30 - crhs31*crhs35 - crhs31*crhs36;
            rhs[7]=-DN(1,0)*crhs22 + DN(1,0)*crhs32 - DN(1,1)*crhs27 + DN(1,1)*crhs33 - DN(1,2)*crhs31 + DN(1,2)*crhs34 - N[1]*crhs17;
            rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs18 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - N[2]*crhs10 - N[2]*crhs2 - crhs22*crhs37 - crhs22*crhs38;
            rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs18 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs24 - N[2]*crhs25 - N[2]*crhs26 - crhs27*crhs37 - crhs27*crhs38;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs18 - DN(2,2)*stress[2] + N[2]*crhs28 - N[2]*crhs29 - N[2]*crhs30 - crhs31*crhs37 - crhs31*crhs38;
            rhs[11]=-DN(2,0)*crhs22 + DN(2,0)*crhs32 - DN(2,1)*crhs27 + DN(2,1)*crhs33 - DN(2,2)*crhs31 + DN(2,2)*crhs34 - N[2]*crhs17;
            rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs18 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - N[3]*crhs10 - N[3]*crhs2 - crhs22*crhs39 - crhs22*crhs40;
            rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs18 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs24 - N[3]*crhs25 - N[3]*crhs26 - crhs27*crhs39 - crhs27*crhs40;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs18 - DN(3,2)*stress[2] + N[3]*crhs28 - N[3]*crhs29 - N[3]*crhs30 - crhs31*crhs39 - crhs31*crhs40;
            rhs[15]=-DN(3,0)*crhs22 + DN(3,0)*crhs32 - DN(3,1)*crhs27 + DN(3,1)*crhs33 - DN(3,2)*crhs31 + DN(3,2)*crhs34 - N[3]*crhs17;


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

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crhs2 =             rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)));
const double crhs3 =             DN(0,0)*v(0,0);
const double crhs4 =             DN(1,0)*v(1,0);
const double crhs5 =             DN(2,0)*v(2,0);
const double crhs6 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
const double crhs7 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
const double crhs8 =             rho*(crhs6*(crhs3 + crhs4 + crhs5) + crhs7*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
const double crhs9 =             rho*stab_c2*sqrt(pow(crhs6, 2) + pow(crhs7, 2));
const double crhs10 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs11 =             (N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]))/(pow(c, 2)*rho);
const double crhs12 =             (crhs9*h/stab_c1 + mu)*(crhs10 + crhs11 + crhs3 + crhs4 + crhs5);
const double crhs13 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
const double crhs14 =             N[0]*crhs13*rho;
const double crhs15 =             1.0/(crhs9/h + mu*stab_c1/pow(h, 2) + dyn_tau*rho/dt);
const double crhs16 =             1.0*crhs15*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + crhs2 + crhs8);
const double crhs17 =             rho*(DN(0,0)*crhs6 + DN(0,1)*crhs7);
const double crhs18 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs19 =             rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)));
const double crhs20 =             rho*(crhs10*crhs7 + crhs6*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)));
const double crhs21 =             1.0*crhs15*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs18 + crhs19 + crhs20);
const double crhs22 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0);
const double crhs23 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
const double crhs24 =             N[1]*crhs13*rho;
const double crhs25 =             rho*(DN(1,0)*crhs6 + DN(1,1)*crhs7);
const double crhs26 =             N[2]*crhs13*rho;
const double crhs27 =             rho*(DN(2,0)*crhs6 + DN(2,1)*crhs7);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs12 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - N[0]*crhs2 - N[0]*crhs8 - crhs14*crhs16 - crhs16*crhs17;
            rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs12 - DN(0,1)*stress[1] + N[0]*crhs18 - N[0]*crhs19 - N[0]*crhs20 - crhs14*crhs21 - crhs17*crhs21;
            rhs[2]=-DN(0,0)*crhs16 + DN(0,0)*crhs22 - DN(0,1)*crhs21 + DN(0,1)*crhs23 - N[0]*crhs11;
            rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs12 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - N[1]*crhs2 - N[1]*crhs8 - crhs16*crhs24 - crhs16*crhs25;
            rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs12 - DN(1,1)*stress[1] + N[1]*crhs18 - N[1]*crhs19 - N[1]*crhs20 - crhs21*crhs24 - crhs21*crhs25;
            rhs[5]=-DN(1,0)*crhs16 + DN(1,0)*crhs22 - DN(1,1)*crhs21 + DN(1,1)*crhs23 - N[1]*crhs11;
            rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs12 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - N[2]*crhs2 - N[2]*crhs8 - crhs16*crhs26 - crhs16*crhs27;
            rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs12 - DN(2,1)*stress[1] + N[2]*crhs18 - N[2]*crhs19 - N[2]*crhs20 - crhs21*crhs26 - crhs21*crhs27;
            rhs[8]=-DN(2,0)*crhs16 + DN(2,0)*crhs22 - DN(2,1)*crhs21 + DN(2,1)*crhs23 - N[2]*crhs11;


}

}
