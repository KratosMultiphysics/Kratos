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
void NavierStokes<3>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int Dim = 3;
    unsigned int NumNodes = 4;
    unsigned int DofSize  = NumNodes*(Dim+1);

    if (rResult.size() != DofSize)
        rResult.resize(DofSize, false);

    for(unsigned int i=0; i<NumNodes; i++)
    {
        rResult[i*(Dim+1)  ]  =  this->GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
        rResult[i*(Dim+1)+1]  =  this->GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
        rResult[i*(Dim+1)+2]  =  this->GetGeometry()[i].GetDof(VELOCITY_Z).EquationId();
        rResult[i*(Dim+1)+3]  =  this->GetGeometry()[i].GetDof(PRESSURE).EquationId();
    }

    KRATOS_CATCH("")
}


template<>
void NavierStokes<2>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int Dim = 2;
    unsigned int NumNodes = 3;
    unsigned int DofSize  = NumNodes*(Dim+1);

    if (rResult.size() != DofSize)
        rResult.resize(DofSize, false);

    for(unsigned int i=0; i<NumNodes; i++)
    {
        rResult[i*(Dim+1)  ]  =  this->GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
        rResult[i*(Dim+1)+1]  =  this->GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
        rResult[i*(Dim+1)+2]  =  this->GetGeometry()[i].GetDof(PRESSURE).EquationId();
    }

    KRATOS_CATCH("")
}


template<>
void NavierStokes<3>::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int Dim = 3;
    unsigned int NumNodes = 4;
    unsigned int DofSize  = NumNodes*(Dim+1);

    if (ElementalDofList.size() != DofSize)
        ElementalDofList.resize(DofSize);

    for(unsigned int i=0; i<NumNodes; i++)
    {
        ElementalDofList[i*(Dim+1)  ]  =  this->GetGeometry()[i].pGetDof(VELOCITY_X);
        ElementalDofList[i*(Dim+1)+1]  =  this->GetGeometry()[i].pGetDof(VELOCITY_Y);
        ElementalDofList[i*(Dim+1)+2]  =  this->GetGeometry()[i].pGetDof(VELOCITY_Z);
        ElementalDofList[i*(Dim+1)+3]  =  this->GetGeometry()[i].pGetDof(PRESSURE);
    }

    KRATOS_CATCH("");
}


template<>
void NavierStokes<2>::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int Dim = 2;
    unsigned int NumNodes = 3;
    unsigned int DofSize  = NumNodes*(Dim+1);

    if (ElementalDofList.size() != DofSize)
        ElementalDofList.resize(DofSize);

    for(unsigned int i=0; i<NumNodes; i++)
    {
        ElementalDofList[i*(Dim+1)  ]  =  this->GetGeometry()[i].pGetDof(VELOCITY_X);
        ElementalDofList[i*(Dim+1)+1]  =  this->GetGeometry()[i].pGetDof(VELOCITY_Y);
        ElementalDofList[i*(Dim+1)+2]  =  this->GetGeometry()[i].pGetDof(PRESSURE);
    }

    KRATOS_CATCH("");
}


template<>
void NavierStokes<3>::ComputeGaussPointLHSContribution(bounded_matrix<double,16,16>& lhs, const ElementDataStruct& data)
{
    const int nnodes = 4;
    const int dim = 3;
    //~ const int strain_size = 6;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size
    const double c = data.c;                                // Wave velocity

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& delta_t = data.delta_t;
    const double& dyn_tau_coeff = data.dyn_tau_coeff;

    const bounded_matrix<double,nnodes,dim>& v = data.v;
    const bounded_matrix<double,nnodes,dim>& vn = data.vn;
    const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
    const bounded_matrix<double,nnodes,dim>& vmesh = data.vmesh;
    const bounded_matrix<double,nnodes,dim>& vconv = v - vmesh;
    const bounded_matrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& p = data.p;
    const array_1d<double,nnodes>& pn = data.pn;
    const array_1d<double,nnodes>& pnn = data.pnn;
    //~ const array_1d<double,strain_size>& stress = data.stress;

    // Get constitutive matrix
    const Matrix& C = data.C;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    const array_1d<double,dim> vconv_gauss = prod(trans(vconv), N);

    const double vconv_norm = norm_2(vconv_gauss);

    // Stabilization parameters
    const double c1 = 4.0;
    const double c2 = 2.0;
    const double tau1 = 1.0/((rho*dyn_tau_coeff)/delta_t + (c2*rho*vconv_norm)/h + (c1*mu)/(h*h));
    const double tau2 = (h*h)/(c1*tau1);

    const double clhs0 =             pow(DN(0,0), 2);
const double clhs1 =             pow(N[0], 2);
const double clhs2 =             bdf0*rho;
const double clhs3 =             clhs1*clhs2;
const double clhs4 =             C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs5 =             C(0,3)*DN(0,0);
const double clhs6 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs5;
const double clhs7 =             C(0,5)*DN(0,0);
const double clhs8 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs7;
const double clhs9 =             N[0]*rho;
const double clhs10 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0);
const double clhs11 =             N[0]*(v(0,0) - vmesh(0,0)) + N[1]*(v(1,0) - vmesh(1,0)) + N[2]*(v(2,0) - vmesh(2,0)) + N[3]*(v(3,0) - vmesh(3,0));
const double clhs12 =             N[0]*(v(0,1) - vmesh(0,1)) + N[1]*(v(1,1) - vmesh(1,1)) + N[2]*(v(2,1) - vmesh(2,1)) + N[3]*(v(3,1) - vmesh(3,1));
const double clhs13 =             N[0]*(v(0,2) - vmesh(0,2)) + N[1]*(v(1,2) - vmesh(1,2)) + N[2]*(v(2,2) - vmesh(2,2)) + N[3]*(v(3,2) - vmesh(3,2));
const double clhs14 =             DN(0,0)*clhs11 + DN(0,1)*clhs12 + DN(0,2)*clhs13;
const double clhs15 =             N[0]*clhs10 + clhs14;
const double clhs16 =             pow(rho, 2);
const double clhs17 =             clhs14*clhs16*tau1;
const double clhs18 =             N[0]*bdf0;
const double clhs19 =             clhs15 + clhs18;
const double clhs20 =             DN(0,0)*N[0]*rho*tau1;
const double clhs21 =             DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0);
const double clhs22 =             DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0);
const double clhs23 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)) + clhs10*clhs11 + clhs12*clhs21 + clhs13*clhs22);
const double clhs24 =             DN(0,0)*tau2;
const double clhs25 =             DN(0,1)*clhs24;
const double clhs26 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs5;
const double clhs27 =             C(1,3)*DN(0,1);
const double clhs28 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs27;
const double clhs29 =             C(3,5)*DN(0,0);
const double clhs30 =             C(4,5)*DN(0,2);
const double clhs31 =             C(1,5)*DN(0,1) + clhs29 + clhs30;
const double clhs32 =             clhs1*rho;
const double clhs33 =             N[0]*clhs14*clhs16*tau1;
const double clhs34 =             DN(0,1)*N[0]*rho*tau1;
const double clhs35 =             DN(0,2)*clhs24;
const double clhs36 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs7;
const double clhs37 =             C(3,4)*DN(0,1);
const double clhs38 =             C(2,3)*DN(0,2) + clhs29 + clhs37;
const double clhs39 =             C(2,5)*DN(0,2);
const double clhs40 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs39;
const double clhs41 =             DN(0,2)*N[0]*rho*tau1;
const double clhs42 =             rho*tau1;
const double clhs43 =             clhs14*clhs42;
const double clhs44 =             -N[0] + clhs43;
const double clhs45 =             DN(0,0)*DN(1,0);
const double clhs46 =             N[0]*bdf0*rho;
const double clhs47 =             N[1]*clhs46;
const double clhs48 =             clhs45*tau2 + clhs47;
const double clhs49 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs50 =             C(0,3)*DN(1,0);
const double clhs51 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs50;
const double clhs52 =             C(0,5)*DN(1,0);
const double clhs53 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs52;
const double clhs54 =             DN(1,0)*clhs11 + DN(1,1)*clhs12 + DN(1,2)*clhs13;
const double clhs55 =             N[1]*clhs10 + clhs54;
const double clhs56 =             N[1]*bdf0;
const double clhs57 =             clhs55 + clhs56;
const double clhs58 =             DN(0,0)*N[1]*rho*tau1;
const double clhs59 =             DN(1,1)*clhs24;
const double clhs60 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs50;
const double clhs61 =             C(1,3)*DN(1,1);
const double clhs62 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs61;
const double clhs63 =             C(3,5)*DN(1,0);
const double clhs64 =             C(4,5)*DN(1,2);
const double clhs65 =             C(1,5)*DN(1,1) + clhs63 + clhs64;
const double clhs66 =             N[0]*N[1]*rho;
const double clhs67 =             clhs21*clhs66;
const double clhs68 =             N[1]*clhs14*clhs16*tau1;
const double clhs69 =             DN(0,1)*N[1]*rho*tau1;
const double clhs70 =             DN(1,2)*clhs24;
const double clhs71 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs52;
const double clhs72 =             C(3,4)*DN(1,1);
const double clhs73 =             C(2,3)*DN(1,2) + clhs63 + clhs72;
const double clhs74 =             C(2,5)*DN(1,2);
const double clhs75 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs74;
const double clhs76 =             clhs22*clhs66;
const double clhs77 =             DN(0,2)*N[1]*rho*tau1;
const double clhs78 =             DN(0,0)*N[1];
const double clhs79 =             DN(0,0)*DN(2,0);
const double clhs80 =             N[2]*clhs46;
const double clhs81 =             clhs79*tau2 + clhs80;
const double clhs82 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs83 =             C(0,3)*DN(2,0);
const double clhs84 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs83;
const double clhs85 =             C(0,5)*DN(2,0);
const double clhs86 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs85;
const double clhs87 =             DN(2,0)*clhs11 + DN(2,1)*clhs12 + DN(2,2)*clhs13;
const double clhs88 =             N[2]*clhs10 + clhs87;
const double clhs89 =             N[2]*bdf0;
const double clhs90 =             clhs88 + clhs89;
const double clhs91 =             DN(0,0)*N[2]*rho*tau1;
const double clhs92 =             DN(2,1)*clhs24;
const double clhs93 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs83;
const double clhs94 =             C(1,3)*DN(2,1);
const double clhs95 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs94;
const double clhs96 =             C(3,5)*DN(2,0);
const double clhs97 =             C(4,5)*DN(2,2);
const double clhs98 =             C(1,5)*DN(2,1) + clhs96 + clhs97;
const double clhs99 =             N[0]*N[2]*rho;
const double clhs100 =             clhs21*clhs99;
const double clhs101 =             N[2]*clhs14*clhs16*tau1;
const double clhs102 =             DN(0,1)*N[2]*rho*tau1;
const double clhs103 =             DN(2,2)*clhs24;
const double clhs104 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs85;
const double clhs105 =             C(3,4)*DN(2,1);
const double clhs106 =             C(2,3)*DN(2,2) + clhs105 + clhs96;
const double clhs107 =             C(2,5)*DN(2,2);
const double clhs108 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs107;
const double clhs109 =             clhs22*clhs99;
const double clhs110 =             DN(0,2)*N[2]*rho*tau1;
const double clhs111 =             DN(0,0)*N[2];
const double clhs112 =             DN(0,0)*DN(3,0);
const double clhs113 =             N[3]*clhs46;
const double clhs114 =             clhs112*tau2 + clhs113;
const double clhs115 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs116 =             C(0,3)*DN(3,0);
const double clhs117 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs116;
const double clhs118 =             C(0,5)*DN(3,0);
const double clhs119 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs118;
const double clhs120 =             DN(3,0)*clhs11 + DN(3,1)*clhs12 + DN(3,2)*clhs13;
const double clhs121 =             N[3]*clhs10 + clhs120;
const double clhs122 =             N[3]*bdf0;
const double clhs123 =             clhs121 + clhs122;
const double clhs124 =             DN(0,0)*N[3]*rho*tau1;
const double clhs125 =             DN(3,1)*clhs24;
const double clhs126 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs116;
const double clhs127 =             C(1,3)*DN(3,1);
const double clhs128 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs127;
const double clhs129 =             C(3,5)*DN(3,0);
const double clhs130 =             C(4,5)*DN(3,2);
const double clhs131 =             C(1,5)*DN(3,1) + clhs129 + clhs130;
const double clhs132 =             N[0]*N[3]*rho;
const double clhs133 =             clhs132*clhs21;
const double clhs134 =             N[3]*clhs14*clhs16*tau1;
const double clhs135 =             DN(0,1)*N[3]*rho*tau1;
const double clhs136 =             DN(3,2)*clhs24;
const double clhs137 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs118;
const double clhs138 =             C(3,4)*DN(3,1);
const double clhs139 =             C(2,3)*DN(3,2) + clhs129 + clhs138;
const double clhs140 =             C(2,5)*DN(3,2);
const double clhs141 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs140;
const double clhs142 =             clhs132*clhs22;
const double clhs143 =             DN(0,2)*N[3]*rho*tau1;
const double clhs144 =             DN(0,0)*N[3];
const double clhs145 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs27;
const double clhs146 =             C(0,4)*DN(0,0) + clhs30 + clhs37;
const double clhs147 =             DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1);
const double clhs148 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
const double clhs149 =             DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1);
const double clhs150 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)) + clhs11*clhs147 + clhs12*clhs148 + clhs13*clhs149);
const double clhs151 =             pow(DN(0,1), 2);
const double clhs152 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs153 =             C(1,4)*DN(0,1);
const double clhs154 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs153;
const double clhs155 =             N[0]*clhs148 + clhs14;
const double clhs156 =             clhs155 + clhs18;
const double clhs157 =             DN(0,1)*tau2;
const double clhs158 =             DN(0,2)*clhs157;
const double clhs159 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs153;
const double clhs160 =             C(2,4)*DN(0,2);
const double clhs161 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs160;
const double clhs162 =             DN(1,0)*clhs157;
const double clhs163 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs61;
const double clhs164 =             C(0,4)*DN(1,0) + clhs64 + clhs72;
const double clhs165 =             clhs147*clhs66;
const double clhs166 =             DN(0,1)*DN(1,1);
const double clhs167 =             clhs166*tau2 + clhs47;
const double clhs168 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs169 =             C(1,4)*DN(1,1);
const double clhs170 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs169;
const double clhs171 =             N[1]*clhs148 + clhs54;
const double clhs172 =             clhs171 + clhs56;
const double clhs173 =             DN(1,2)*clhs157;
const double clhs174 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs169;
const double clhs175 =             C(2,4)*DN(1,2);
const double clhs176 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs175;
const double clhs177 =             clhs149*clhs66;
const double clhs178 =             DN(0,1)*N[1];
const double clhs179 =             DN(2,0)*clhs157;
const double clhs180 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs94;
const double clhs181 =             C(0,4)*DN(2,0) + clhs105 + clhs97;
const double clhs182 =             clhs147*clhs99;
const double clhs183 =             DN(0,1)*DN(2,1);
const double clhs184 =             clhs183*tau2 + clhs80;
const double clhs185 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs186 =             C(1,4)*DN(2,1);
const double clhs187 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs186;
const double clhs188 =             N[2]*clhs148 + clhs87;
const double clhs189 =             clhs188 + clhs89;
const double clhs190 =             DN(2,2)*clhs157;
const double clhs191 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs186;
const double clhs192 =             C(2,4)*DN(2,2);
const double clhs193 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs192;
const double clhs194 =             clhs149*clhs99;
const double clhs195 =             DN(0,1)*N[2];
const double clhs196 =             DN(3,0)*clhs157;
const double clhs197 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs127;
const double clhs198 =             C(0,4)*DN(3,0) + clhs130 + clhs138;
const double clhs199 =             clhs132*clhs147;
const double clhs200 =             DN(0,1)*DN(3,1);
const double clhs201 =             clhs113 + clhs200*tau2;
const double clhs202 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs203 =             C(1,4)*DN(3,1);
const double clhs204 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs203;
const double clhs205 =             N[3]*clhs148 + clhs120;
const double clhs206 =             clhs122 + clhs205;
const double clhs207 =             DN(3,2)*clhs157;
const double clhs208 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs203;
const double clhs209 =             C(2,4)*DN(3,2);
const double clhs210 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs209;
const double clhs211 =             clhs132*clhs149;
const double clhs212 =             DN(0,1)*N[3];
const double clhs213 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs39;
const double clhs214 =             DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2);
const double clhs215 =             DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2);
const double clhs216 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double clhs217 =             DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)) + clhs11*clhs214 + clhs12*clhs215 + clhs13*clhs216);
const double clhs218 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs160;
const double clhs219 =             pow(DN(0,2), 2);
const double clhs220 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs221 =             N[0]*clhs216 + clhs14;
const double clhs222 =             clhs18 + clhs221;
const double clhs223 =             DN(0,2)*tau2;
const double clhs224 =             DN(1,0)*clhs223;
const double clhs225 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs74;
const double clhs226 =             clhs214*clhs66;
const double clhs227 =             DN(1,1)*clhs223;
const double clhs228 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs175;
const double clhs229 =             clhs215*clhs66;
const double clhs230 =             DN(0,2)*DN(1,2);
const double clhs231 =             clhs230*tau2 + clhs47;
const double clhs232 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs233 =             N[1]*clhs216 + clhs54;
const double clhs234 =             clhs233 + clhs56;
const double clhs235 =             DN(0,2)*N[1];
const double clhs236 =             DN(2,0)*clhs223;
const double clhs237 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs107;
const double clhs238 =             clhs214*clhs99;
const double clhs239 =             DN(2,1)*clhs223;
const double clhs240 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs192;
const double clhs241 =             clhs215*clhs99;
const double clhs242 =             DN(0,2)*DN(2,2);
const double clhs243 =             clhs242*tau2 + clhs80;
const double clhs244 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs245 =             N[2]*clhs216 + clhs87;
const double clhs246 =             clhs245 + clhs89;
const double clhs247 =             DN(0,2)*N[2];
const double clhs248 =             DN(3,0)*clhs223;
const double clhs249 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs140;
const double clhs250 =             clhs132*clhs214;
const double clhs251 =             DN(3,1)*clhs223;
const double clhs252 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs209;
const double clhs253 =             clhs132*clhs215;
const double clhs254 =             DN(0,2)*DN(3,2);
const double clhs255 =             clhs113 + clhs254*tau2;
const double clhs256 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs257 =             N[3]*clhs216 + clhs120;
const double clhs258 =             clhs122 + clhs257;
const double clhs259 =             DN(0,2)*N[3];
const double clhs260 =             DN(0,0)*rho*tau1;
const double clhs261 =             DN(0,1)*rho*tau1;
const double clhs262 =             DN(0,2)*rho*tau1;
const double clhs263 =             DN(1,0)*N[0];
const double clhs264 =             DN(1,1)*N[0];
const double clhs265 =             DN(1,2)*N[0];
const double clhs266 =             tau1*(clhs166 + clhs230 + clhs45);
const double clhs267 =             DN(2,0)*N[0];
const double clhs268 =             DN(2,1)*N[0];
const double clhs269 =             DN(2,2)*N[0];
const double clhs270 =             tau1*(clhs183 + clhs242 + clhs79);
const double clhs271 =             DN(3,0)*N[0];
const double clhs272 =             DN(3,1)*N[0];
const double clhs273 =             DN(3,2)*N[0];
const double clhs274 =             tau1*(clhs112 + clhs200 + clhs254);
const double clhs275 =             N[1]*rho;
const double clhs276 =             clhs16*clhs54*tau1;
const double clhs277 =             DN(1,0)*N[0]*rho*tau1;
const double clhs278 =             N[0]*clhs16*clhs54*tau1;
const double clhs279 =             DN(1,1)*N[0]*rho*tau1;
const double clhs280 =             DN(1,2)*N[0]*rho*tau1;
const double clhs281 =             clhs42*clhs54;
const double clhs282 =             pow(DN(1,0), 2);
const double clhs283 =             pow(N[1], 2);
const double clhs284 =             clhs2*clhs283;
const double clhs285 =             DN(1,0)*N[1]*rho*tau1;
const double clhs286 =             DN(1,0)*tau2;
const double clhs287 =             DN(1,1)*clhs286;
const double clhs288 =             clhs283*rho;
const double clhs289 =             N[1]*clhs16*clhs54*tau1;
const double clhs290 =             DN(1,1)*N[1]*rho*tau1;
const double clhs291 =             DN(1,2)*clhs286;
const double clhs292 =             DN(1,2)*N[1]*rho*tau1;
const double clhs293 =             -N[1] + clhs281;
const double clhs294 =             DN(1,0)*DN(2,0);
const double clhs295 =             N[1]*bdf0*rho;
const double clhs296 =             N[2]*clhs295;
const double clhs297 =             clhs294*tau2 + clhs296;
const double clhs298 =             DN(1,0)*N[2]*rho*tau1;
const double clhs299 =             DN(2,1)*clhs286;
const double clhs300 =             N[1]*N[2]*rho;
const double clhs301 =             clhs21*clhs300;
const double clhs302 =             N[2]*clhs16*clhs54*tau1;
const double clhs303 =             DN(1,1)*N[2]*rho*tau1;
const double clhs304 =             DN(2,2)*clhs286;
const double clhs305 =             clhs22*clhs300;
const double clhs306 =             DN(1,2)*N[2]*rho*tau1;
const double clhs307 =             DN(1,0)*N[2];
const double clhs308 =             DN(1,0)*DN(3,0);
const double clhs309 =             N[3]*clhs295;
const double clhs310 =             clhs308*tau2 + clhs309;
const double clhs311 =             DN(1,0)*N[3]*rho*tau1;
const double clhs312 =             DN(3,1)*clhs286;
const double clhs313 =             N[1]*N[3]*rho;
const double clhs314 =             clhs21*clhs313;
const double clhs315 =             N[3]*clhs16*clhs54*tau1;
const double clhs316 =             DN(1,1)*N[3]*rho*tau1;
const double clhs317 =             DN(3,2)*clhs286;
const double clhs318 =             clhs22*clhs313;
const double clhs319 =             DN(1,2)*N[3]*rho*tau1;
const double clhs320 =             DN(1,0)*N[3];
const double clhs321 =             pow(DN(1,1), 2);
const double clhs322 =             DN(1,1)*tau2;
const double clhs323 =             DN(1,2)*clhs322;
const double clhs324 =             DN(2,0)*clhs322;
const double clhs325 =             clhs147*clhs300;
const double clhs326 =             DN(1,1)*DN(2,1);
const double clhs327 =             clhs296 + clhs326*tau2;
const double clhs328 =             DN(2,2)*clhs322;
const double clhs329 =             clhs149*clhs300;
const double clhs330 =             DN(1,1)*N[2];
const double clhs331 =             DN(3,0)*clhs322;
const double clhs332 =             clhs147*clhs313;
const double clhs333 =             DN(1,1)*DN(3,1);
const double clhs334 =             clhs309 + clhs333*tau2;
const double clhs335 =             DN(3,2)*clhs322;
const double clhs336 =             clhs149*clhs313;
const double clhs337 =             DN(1,1)*N[3];
const double clhs338 =             pow(DN(1,2), 2);
const double clhs339 =             DN(1,2)*tau2;
const double clhs340 =             DN(2,0)*clhs339;
const double clhs341 =             clhs214*clhs300;
const double clhs342 =             DN(2,1)*clhs339;
const double clhs343 =             clhs215*clhs300;
const double clhs344 =             DN(1,2)*DN(2,2);
const double clhs345 =             clhs296 + clhs344*tau2;
const double clhs346 =             DN(1,2)*N[2];
const double clhs347 =             DN(3,0)*clhs339;
const double clhs348 =             clhs214*clhs313;
const double clhs349 =             DN(3,1)*clhs339;
const double clhs350 =             clhs215*clhs313;
const double clhs351 =             DN(1,2)*DN(3,2);
const double clhs352 =             clhs309 + clhs351*tau2;
const double clhs353 =             DN(1,2)*N[3];
const double clhs354 =             DN(1,0)*rho*tau1;
const double clhs355 =             DN(1,1)*rho*tau1;
const double clhs356 =             DN(1,2)*rho*tau1;
const double clhs357 =             DN(2,0)*N[1];
const double clhs358 =             DN(2,1)*N[1];
const double clhs359 =             DN(2,2)*N[1];
const double clhs360 =             tau1*(clhs294 + clhs326 + clhs344);
const double clhs361 =             DN(3,0)*N[1];
const double clhs362 =             DN(3,1)*N[1];
const double clhs363 =             DN(3,2)*N[1];
const double clhs364 =             tau1*(clhs308 + clhs333 + clhs351);
const double clhs365 =             N[2]*rho;
const double clhs366 =             clhs16*clhs87*tau1;
const double clhs367 =             DN(2,0)*N[0]*rho*tau1;
const double clhs368 =             N[0]*clhs16*clhs87*tau1;
const double clhs369 =             DN(2,1)*N[0]*rho*tau1;
const double clhs370 =             DN(2,2)*N[0]*rho*tau1;
const double clhs371 =             clhs42*clhs87;
const double clhs372 =             DN(2,0)*N[1]*rho*tau1;
const double clhs373 =             N[1]*clhs16*clhs87*tau1;
const double clhs374 =             DN(2,1)*N[1]*rho*tau1;
const double clhs375 =             DN(2,2)*N[1]*rho*tau1;
const double clhs376 =             pow(DN(2,0), 2);
const double clhs377 =             pow(N[2], 2);
const double clhs378 =             clhs2*clhs377;
const double clhs379 =             DN(2,0)*N[2]*rho*tau1;
const double clhs380 =             DN(2,0)*tau2;
const double clhs381 =             DN(2,1)*clhs380;
const double clhs382 =             clhs377*rho;
const double clhs383 =             N[2]*clhs16*clhs87*tau1;
const double clhs384 =             DN(2,1)*N[2]*rho*tau1;
const double clhs385 =             DN(2,2)*clhs380;
const double clhs386 =             DN(2,2)*N[2]*rho*tau1;
const double clhs387 =             -N[2] + clhs371;
const double clhs388 =             DN(2,0)*DN(3,0);
const double clhs389 =             N[3]*rho;
const double clhs390 =             clhs389*clhs89;
const double clhs391 =             clhs388*tau2 + clhs390;
const double clhs392 =             DN(2,0)*N[3]*rho*tau1;
const double clhs393 =             DN(3,1)*clhs380;
const double clhs394 =             N[2]*N[3]*rho;
const double clhs395 =             clhs21*clhs394;
const double clhs396 =             N[3]*clhs16*clhs87*tau1;
const double clhs397 =             DN(2,1)*N[3]*rho*tau1;
const double clhs398 =             DN(3,2)*clhs380;
const double clhs399 =             clhs22*clhs394;
const double clhs400 =             DN(2,2)*N[3]*rho*tau1;
const double clhs401 =             DN(2,0)*N[3];
const double clhs402 =             pow(DN(2,1), 2);
const double clhs403 =             DN(2,1)*tau2;
const double clhs404 =             DN(2,2)*clhs403;
const double clhs405 =             DN(3,0)*clhs403;
const double clhs406 =             clhs147*clhs394;
const double clhs407 =             DN(2,1)*DN(3,1);
const double clhs408 =             clhs390 + clhs407*tau2;
const double clhs409 =             DN(3,2)*clhs403;
const double clhs410 =             clhs149*clhs394;
const double clhs411 =             DN(2,1)*N[3];
const double clhs412 =             pow(DN(2,2), 2);
const double clhs413 =             DN(2,2)*tau2;
const double clhs414 =             DN(3,0)*clhs413;
const double clhs415 =             clhs214*clhs394;
const double clhs416 =             DN(3,1)*clhs413;
const double clhs417 =             clhs215*clhs394;
const double clhs418 =             DN(2,2)*DN(3,2);
const double clhs419 =             clhs390 + clhs418*tau2;
const double clhs420 =             DN(2,2)*N[3];
const double clhs421 =             DN(2,0)*rho*tau1;
const double clhs422 =             DN(2,1)*rho*tau1;
const double clhs423 =             DN(2,2)*rho*tau1;
const double clhs424 =             DN(3,0)*N[2];
const double clhs425 =             DN(3,1)*N[2];
const double clhs426 =             DN(3,2)*N[2];
const double clhs427 =             tau1*(clhs388 + clhs407 + clhs418);
const double clhs428 =             clhs120*clhs16*tau1;
const double clhs429 =             DN(3,0)*N[0]*rho*tau1;
const double clhs430 =             N[0]*clhs120*clhs16*tau1;
const double clhs431 =             DN(3,1)*N[0]*rho*tau1;
const double clhs432 =             DN(3,2)*N[0]*rho*tau1;
const double clhs433 =             clhs120*clhs42;
const double clhs434 =             DN(3,0)*N[1]*rho*tau1;
const double clhs435 =             N[1]*clhs120*clhs16*tau1;
const double clhs436 =             DN(3,1)*N[1]*rho*tau1;
const double clhs437 =             DN(3,2)*N[1]*rho*tau1;
const double clhs438 =             DN(3,0)*N[2]*rho*tau1;
const double clhs439 =             N[2]*clhs120*clhs16*tau1;
const double clhs440 =             DN(3,1)*N[2]*rho*tau1;
const double clhs441 =             DN(3,2)*N[2]*rho*tau1;
const double clhs442 =             pow(DN(3,0), 2);
const double clhs443 =             pow(N[3], 2);
const double clhs444 =             clhs2*clhs443;
const double clhs445 =             DN(3,0)*N[3]*rho*tau1;
const double clhs446 =             DN(3,0)*tau2;
const double clhs447 =             DN(3,1)*clhs446;
const double clhs448 =             clhs443*rho;
const double clhs449 =             N[3]*clhs120*clhs16*tau1;
const double clhs450 =             DN(3,1)*N[3]*rho*tau1;
const double clhs451 =             DN(3,2)*clhs446;
const double clhs452 =             DN(3,2)*N[3]*rho*tau1;
const double clhs453 =             -N[3] + clhs433;
const double clhs454 =             pow(DN(3,1), 2);
const double clhs455 =             DN(3,1)*DN(3,2)*tau2;
const double clhs456 =             pow(DN(3,2), 2);
const double clhs457 =             DN(3,0)*rho*tau1;
const double clhs458 =             DN(3,1)*rho*tau1;
const double clhs459 =             DN(3,2)*rho*tau1;
            lhs(0,0)=DN(0,0)*clhs4 + DN(0,1)*clhs6 + DN(0,2)*clhs8 + clhs0*tau2 + clhs15*clhs9 + clhs17*clhs19 + clhs20*clhs23 + clhs3;
            lhs(0,1)=DN(0,0)*clhs26 + DN(0,1)*clhs28 + DN(0,2)*clhs31 + clhs21*clhs32 + clhs21*clhs33 + clhs23*clhs34 + clhs25;
            lhs(0,2)=DN(0,0)*clhs36 + DN(0,1)*clhs38 + DN(0,2)*clhs40 + clhs22*clhs32 + clhs22*clhs33 + clhs23*clhs41 + clhs35;
            lhs(0,3)=DN(0,0)*clhs44;
            lhs(0,4)=DN(0,0)*clhs49 + DN(0,1)*clhs51 + DN(0,2)*clhs53 + clhs17*clhs57 + clhs23*clhs58 + clhs48 + clhs55*clhs9;
            lhs(0,5)=DN(0,0)*clhs60 + DN(0,1)*clhs62 + DN(0,2)*clhs65 + clhs21*clhs68 + clhs23*clhs69 + clhs59 + clhs67;
            lhs(0,6)=DN(0,0)*clhs71 + DN(0,1)*clhs73 + DN(0,2)*clhs75 + clhs22*clhs68 + clhs23*clhs77 + clhs70 + clhs76;
            lhs(0,7)=DN(1,0)*clhs43 - clhs78;
            lhs(0,8)=DN(0,0)*clhs82 + DN(0,1)*clhs84 + DN(0,2)*clhs86 + clhs17*clhs90 + clhs23*clhs91 + clhs81 + clhs88*clhs9;
            lhs(0,9)=DN(0,0)*clhs93 + DN(0,1)*clhs95 + DN(0,2)*clhs98 + clhs100 + clhs101*clhs21 + clhs102*clhs23 + clhs92;
            lhs(0,10)=DN(0,0)*clhs104 + DN(0,1)*clhs106 + DN(0,2)*clhs108 + clhs101*clhs22 + clhs103 + clhs109 + clhs110*clhs23;
            lhs(0,11)=DN(2,0)*clhs43 - clhs111;
            lhs(0,12)=DN(0,0)*clhs115 + DN(0,1)*clhs117 + DN(0,2)*clhs119 + clhs114 + clhs121*clhs9 + clhs123*clhs17 + clhs124*clhs23;
            lhs(0,13)=DN(0,0)*clhs126 + DN(0,1)*clhs128 + DN(0,2)*clhs131 + clhs125 + clhs133 + clhs134*clhs21 + clhs135*clhs23;
            lhs(0,14)=DN(0,0)*clhs137 + DN(0,1)*clhs139 + DN(0,2)*clhs141 + clhs134*clhs22 + clhs136 + clhs142 + clhs143*clhs23;
            lhs(0,15)=DN(3,0)*clhs43 - clhs144;
            lhs(1,0)=DN(0,0)*clhs6 + DN(0,1)*clhs145 + DN(0,2)*clhs146 + clhs147*clhs32 + clhs147*clhs33 + clhs150*clhs20 + clhs25;
            lhs(1,1)=DN(0,0)*clhs28 + DN(0,1)*clhs152 + DN(0,2)*clhs154 + clhs150*clhs34 + clhs151*tau2 + clhs155*clhs9 + clhs156*clhs17 + clhs3;
            lhs(1,2)=DN(0,0)*clhs38 + DN(0,1)*clhs159 + DN(0,2)*clhs161 + clhs149*clhs32 + clhs149*clhs33 + clhs150*clhs41 + clhs158;
            lhs(1,3)=DN(0,1)*clhs44;
            lhs(1,4)=DN(0,0)*clhs51 + DN(0,1)*clhs163 + DN(0,2)*clhs164 + clhs147*clhs68 + clhs150*clhs58 + clhs162 + clhs165;
            lhs(1,5)=DN(0,0)*clhs62 + DN(0,1)*clhs168 + DN(0,2)*clhs170 + clhs150*clhs69 + clhs167 + clhs17*clhs172 + clhs171*clhs9;
            lhs(1,6)=DN(0,0)*clhs73 + DN(0,1)*clhs174 + DN(0,2)*clhs176 + clhs149*clhs68 + clhs150*clhs77 + clhs173 + clhs177;
            lhs(1,7)=DN(1,1)*clhs43 - clhs178;
            lhs(1,8)=DN(0,0)*clhs84 + DN(0,1)*clhs180 + DN(0,2)*clhs181 + clhs101*clhs147 + clhs150*clhs91 + clhs179 + clhs182;
            lhs(1,9)=DN(0,0)*clhs95 + DN(0,1)*clhs185 + DN(0,2)*clhs187 + clhs102*clhs150 + clhs17*clhs189 + clhs184 + clhs188*clhs9;
            lhs(1,10)=DN(0,0)*clhs106 + DN(0,1)*clhs191 + DN(0,2)*clhs193 + clhs101*clhs149 + clhs110*clhs150 + clhs190 + clhs194;
            lhs(1,11)=DN(2,1)*clhs43 - clhs195;
            lhs(1,12)=DN(0,0)*clhs117 + DN(0,1)*clhs197 + DN(0,2)*clhs198 + clhs124*clhs150 + clhs134*clhs147 + clhs196 + clhs199;
            lhs(1,13)=DN(0,0)*clhs128 + DN(0,1)*clhs202 + DN(0,2)*clhs204 + clhs135*clhs150 + clhs17*clhs206 + clhs201 + clhs205*clhs9;
            lhs(1,14)=DN(0,0)*clhs139 + DN(0,1)*clhs208 + DN(0,2)*clhs210 + clhs134*clhs149 + clhs143*clhs150 + clhs207 + clhs211;
            lhs(1,15)=DN(3,1)*clhs43 - clhs212;
            lhs(2,0)=DN(0,0)*clhs8 + DN(0,1)*clhs146 + DN(0,2)*clhs213 + clhs20*clhs217 + clhs214*clhs32 + clhs214*clhs33 + clhs35;
            lhs(2,1)=DN(0,0)*clhs31 + DN(0,1)*clhs154 + DN(0,2)*clhs218 + clhs158 + clhs215*clhs32 + clhs215*clhs33 + clhs217*clhs34;
            lhs(2,2)=DN(0,0)*clhs40 + DN(0,1)*clhs161 + DN(0,2)*clhs220 + clhs17*clhs222 + clhs217*clhs41 + clhs219*tau2 + clhs221*clhs9 + clhs3;
            lhs(2,3)=DN(0,2)*clhs44;
            lhs(2,4)=DN(0,0)*clhs53 + DN(0,1)*clhs164 + DN(0,2)*clhs225 + clhs214*clhs68 + clhs217*clhs58 + clhs224 + clhs226;
            lhs(2,5)=DN(0,0)*clhs65 + DN(0,1)*clhs170 + DN(0,2)*clhs228 + clhs215*clhs68 + clhs217*clhs69 + clhs227 + clhs229;
            lhs(2,6)=DN(0,0)*clhs75 + DN(0,1)*clhs176 + DN(0,2)*clhs232 + clhs17*clhs234 + clhs217*clhs77 + clhs231 + clhs233*clhs9;
            lhs(2,7)=DN(1,2)*clhs43 - clhs235;
            lhs(2,8)=DN(0,0)*clhs86 + DN(0,1)*clhs181 + DN(0,2)*clhs237 + clhs101*clhs214 + clhs217*clhs91 + clhs236 + clhs238;
            lhs(2,9)=DN(0,0)*clhs98 + DN(0,1)*clhs187 + DN(0,2)*clhs240 + clhs101*clhs215 + clhs102*clhs217 + clhs239 + clhs241;
            lhs(2,10)=DN(0,0)*clhs108 + DN(0,1)*clhs193 + DN(0,2)*clhs244 + clhs110*clhs217 + clhs17*clhs246 + clhs243 + clhs245*clhs9;
            lhs(2,11)=DN(2,2)*clhs43 - clhs247;
            lhs(2,12)=DN(0,0)*clhs119 + DN(0,1)*clhs198 + DN(0,2)*clhs249 + clhs124*clhs217 + clhs134*clhs214 + clhs248 + clhs250;
            lhs(2,13)=DN(0,0)*clhs131 + DN(0,1)*clhs204 + DN(0,2)*clhs252 + clhs134*clhs215 + clhs135*clhs217 + clhs251 + clhs253;
            lhs(2,14)=DN(0,0)*clhs141 + DN(0,1)*clhs210 + DN(0,2)*clhs256 + clhs143*clhs217 + clhs17*clhs258 + clhs255 + clhs257*clhs9;
            lhs(2,15)=DN(3,2)*clhs43 - clhs259;
            lhs(3,0)=DN(0,0)*N[0] + clhs147*clhs34 + clhs19*clhs260 + clhs214*clhs41;
            lhs(3,1)=DN(0,1)*N[0] + clhs156*clhs261 + clhs20*clhs21 + clhs215*clhs41;
            lhs(3,2)=DN(0,2)*N[0] + clhs149*clhs34 + clhs20*clhs22 + clhs222*clhs262;
            lhs(3,3)=tau1*(clhs0 + clhs151 + clhs219);
            lhs(3,4)=clhs147*clhs69 + clhs214*clhs77 + clhs260*clhs57 + clhs263;
            lhs(3,5)=clhs172*clhs261 + clhs21*clhs58 + clhs215*clhs77 + clhs264;
            lhs(3,6)=clhs149*clhs69 + clhs22*clhs58 + clhs234*clhs262 + clhs265;
            lhs(3,7)=clhs266;
            lhs(3,8)=clhs102*clhs147 + clhs110*clhs214 + clhs260*clhs90 + clhs267;
            lhs(3,9)=clhs110*clhs215 + clhs189*clhs261 + clhs21*clhs91 + clhs268;
            lhs(3,10)=clhs102*clhs149 + clhs22*clhs91 + clhs246*clhs262 + clhs269;
            lhs(3,11)=clhs270;
            lhs(3,12)=clhs123*clhs260 + clhs135*clhs147 + clhs143*clhs214 + clhs271;
            lhs(3,13)=clhs124*clhs21 + clhs143*clhs215 + clhs206*clhs261 + clhs272;
            lhs(3,14)=clhs124*clhs22 + clhs135*clhs149 + clhs258*clhs262 + clhs273;
            lhs(3,15)=clhs274;
            lhs(4,0)=DN(1,0)*clhs4 + DN(1,1)*clhs6 + DN(1,2)*clhs8 + clhs15*clhs275 + clhs19*clhs276 + clhs23*clhs277 + clhs48;
            lhs(4,1)=DN(1,0)*clhs26 + DN(1,1)*clhs28 + DN(1,2)*clhs31 + clhs162 + clhs21*clhs278 + clhs23*clhs279 + clhs67;
            lhs(4,2)=DN(1,0)*clhs36 + DN(1,1)*clhs38 + DN(1,2)*clhs40 + clhs22*clhs278 + clhs224 + clhs23*clhs280 + clhs76;
            lhs(4,3)=DN(0,0)*clhs281 - clhs263;
            lhs(4,4)=DN(1,0)*clhs49 + DN(1,1)*clhs51 + DN(1,2)*clhs53 + clhs23*clhs285 + clhs275*clhs55 + clhs276*clhs57 + clhs282*tau2 + clhs284;
            lhs(4,5)=DN(1,0)*clhs60 + DN(1,1)*clhs62 + DN(1,2)*clhs65 + clhs21*clhs288 + clhs21*clhs289 + clhs23*clhs290 + clhs287;
            lhs(4,6)=DN(1,0)*clhs71 + DN(1,1)*clhs73 + DN(1,2)*clhs75 + clhs22*clhs288 + clhs22*clhs289 + clhs23*clhs292 + clhs291;
            lhs(4,7)=DN(1,0)*clhs293;
            lhs(4,8)=DN(1,0)*clhs82 + DN(1,1)*clhs84 + DN(1,2)*clhs86 + clhs23*clhs298 + clhs275*clhs88 + clhs276*clhs90 + clhs297;
            lhs(4,9)=DN(1,0)*clhs93 + DN(1,1)*clhs95 + DN(1,2)*clhs98 + clhs21*clhs302 + clhs23*clhs303 + clhs299 + clhs301;
            lhs(4,10)=DN(1,0)*clhs104 + DN(1,1)*clhs106 + DN(1,2)*clhs108 + clhs22*clhs302 + clhs23*clhs306 + clhs304 + clhs305;
            lhs(4,11)=DN(2,0)*clhs281 - clhs307;
            lhs(4,12)=DN(1,0)*clhs115 + DN(1,1)*clhs117 + DN(1,2)*clhs119 + clhs121*clhs275 + clhs123*clhs276 + clhs23*clhs311 + clhs310;
            lhs(4,13)=DN(1,0)*clhs126 + DN(1,1)*clhs128 + DN(1,2)*clhs131 + clhs21*clhs315 + clhs23*clhs316 + clhs312 + clhs314;
            lhs(4,14)=DN(1,0)*clhs137 + DN(1,1)*clhs139 + DN(1,2)*clhs141 + clhs22*clhs315 + clhs23*clhs319 + clhs317 + clhs318;
            lhs(4,15)=DN(3,0)*clhs281 - clhs320;
            lhs(5,0)=DN(1,0)*clhs6 + DN(1,1)*clhs145 + DN(1,2)*clhs146 + clhs147*clhs278 + clhs150*clhs277 + clhs165 + clhs59;
            lhs(5,1)=DN(1,0)*clhs28 + DN(1,1)*clhs152 + DN(1,2)*clhs154 + clhs150*clhs279 + clhs155*clhs275 + clhs156*clhs276 + clhs167;
            lhs(5,2)=DN(1,0)*clhs38 + DN(1,1)*clhs159 + DN(1,2)*clhs161 + clhs149*clhs278 + clhs150*clhs280 + clhs177 + clhs227;
            lhs(5,3)=DN(0,1)*clhs281 - clhs264;
            lhs(5,4)=DN(1,0)*clhs51 + DN(1,1)*clhs163 + DN(1,2)*clhs164 + clhs147*clhs288 + clhs147*clhs289 + clhs150*clhs285 + clhs287;
            lhs(5,5)=DN(1,0)*clhs62 + DN(1,1)*clhs168 + DN(1,2)*clhs170 + clhs150*clhs290 + clhs171*clhs275 + clhs172*clhs276 + clhs284 + clhs321*tau2;
            lhs(5,6)=DN(1,0)*clhs73 + DN(1,1)*clhs174 + DN(1,2)*clhs176 + clhs149*clhs288 + clhs149*clhs289 + clhs150*clhs292 + clhs323;
            lhs(5,7)=DN(1,1)*clhs293;
            lhs(5,8)=DN(1,0)*clhs84 + DN(1,1)*clhs180 + DN(1,2)*clhs181 + clhs147*clhs302 + clhs150*clhs298 + clhs324 + clhs325;
            lhs(5,9)=DN(1,0)*clhs95 + DN(1,1)*clhs185 + DN(1,2)*clhs187 + clhs150*clhs303 + clhs188*clhs275 + clhs189*clhs276 + clhs327;
            lhs(5,10)=DN(1,0)*clhs106 + DN(1,1)*clhs191 + DN(1,2)*clhs193 + clhs149*clhs302 + clhs150*clhs306 + clhs328 + clhs329;
            lhs(5,11)=DN(2,1)*clhs281 - clhs330;
            lhs(5,12)=DN(1,0)*clhs117 + DN(1,1)*clhs197 + DN(1,2)*clhs198 + clhs147*clhs315 + clhs150*clhs311 + clhs331 + clhs332;
            lhs(5,13)=DN(1,0)*clhs128 + DN(1,1)*clhs202 + DN(1,2)*clhs204 + clhs150*clhs316 + clhs205*clhs275 + clhs206*clhs276 + clhs334;
            lhs(5,14)=DN(1,0)*clhs139 + DN(1,1)*clhs208 + DN(1,2)*clhs210 + clhs149*clhs315 + clhs150*clhs319 + clhs335 + clhs336;
            lhs(5,15)=DN(3,1)*clhs281 - clhs337;
            lhs(6,0)=DN(1,0)*clhs8 + DN(1,1)*clhs146 + DN(1,2)*clhs213 + clhs214*clhs278 + clhs217*clhs277 + clhs226 + clhs70;
            lhs(6,1)=DN(1,0)*clhs31 + DN(1,1)*clhs154 + DN(1,2)*clhs218 + clhs173 + clhs215*clhs278 + clhs217*clhs279 + clhs229;
            lhs(6,2)=DN(1,0)*clhs40 + DN(1,1)*clhs161 + DN(1,2)*clhs220 + clhs217*clhs280 + clhs221*clhs275 + clhs222*clhs276 + clhs231;
            lhs(6,3)=DN(0,2)*clhs281 - clhs265;
            lhs(6,4)=DN(1,0)*clhs53 + DN(1,1)*clhs164 + DN(1,2)*clhs225 + clhs214*clhs288 + clhs214*clhs289 + clhs217*clhs285 + clhs291;
            lhs(6,5)=DN(1,0)*clhs65 + DN(1,1)*clhs170 + DN(1,2)*clhs228 + clhs215*clhs288 + clhs215*clhs289 + clhs217*clhs290 + clhs323;
            lhs(6,6)=DN(1,0)*clhs75 + DN(1,1)*clhs176 + DN(1,2)*clhs232 + clhs217*clhs292 + clhs233*clhs275 + clhs234*clhs276 + clhs284 + clhs338*tau2;
            lhs(6,7)=DN(1,2)*clhs293;
            lhs(6,8)=DN(1,0)*clhs86 + DN(1,1)*clhs181 + DN(1,2)*clhs237 + clhs214*clhs302 + clhs217*clhs298 + clhs340 + clhs341;
            lhs(6,9)=DN(1,0)*clhs98 + DN(1,1)*clhs187 + DN(1,2)*clhs240 + clhs215*clhs302 + clhs217*clhs303 + clhs342 + clhs343;
            lhs(6,10)=DN(1,0)*clhs108 + DN(1,1)*clhs193 + DN(1,2)*clhs244 + clhs217*clhs306 + clhs245*clhs275 + clhs246*clhs276 + clhs345;
            lhs(6,11)=DN(2,2)*clhs281 - clhs346;
            lhs(6,12)=DN(1,0)*clhs119 + DN(1,1)*clhs198 + DN(1,2)*clhs249 + clhs214*clhs315 + clhs217*clhs311 + clhs347 + clhs348;
            lhs(6,13)=DN(1,0)*clhs131 + DN(1,1)*clhs204 + DN(1,2)*clhs252 + clhs215*clhs315 + clhs217*clhs316 + clhs349 + clhs350;
            lhs(6,14)=DN(1,0)*clhs141 + DN(1,1)*clhs210 + DN(1,2)*clhs256 + clhs217*clhs319 + clhs257*clhs275 + clhs258*clhs276 + clhs352;
            lhs(6,15)=DN(3,2)*clhs281 - clhs353;
            lhs(7,0)=clhs147*clhs279 + clhs19*clhs354 + clhs214*clhs280 + clhs78;
            lhs(7,1)=clhs156*clhs355 + clhs178 + clhs21*clhs277 + clhs215*clhs280;
            lhs(7,2)=clhs149*clhs279 + clhs22*clhs277 + clhs222*clhs356 + clhs235;
            lhs(7,3)=clhs266;
            lhs(7,4)=DN(1,0)*N[1] + clhs147*clhs290 + clhs214*clhs292 + clhs354*clhs57;
            lhs(7,5)=DN(1,1)*N[1] + clhs172*clhs355 + clhs21*clhs285 + clhs215*clhs292;
            lhs(7,6)=DN(1,2)*N[1] + clhs149*clhs290 + clhs22*clhs285 + clhs234*clhs356;
            lhs(7,7)=tau1*(clhs282 + clhs321 + clhs338);
            lhs(7,8)=clhs147*clhs303 + clhs214*clhs306 + clhs354*clhs90 + clhs357;
            lhs(7,9)=clhs189*clhs355 + clhs21*clhs298 + clhs215*clhs306 + clhs358;
            lhs(7,10)=clhs149*clhs303 + clhs22*clhs298 + clhs246*clhs356 + clhs359;
            lhs(7,11)=clhs360;
            lhs(7,12)=clhs123*clhs354 + clhs147*clhs316 + clhs214*clhs319 + clhs361;
            lhs(7,13)=clhs206*clhs355 + clhs21*clhs311 + clhs215*clhs319 + clhs362;
            lhs(7,14)=clhs149*clhs316 + clhs22*clhs311 + clhs258*clhs356 + clhs363;
            lhs(7,15)=clhs364;
            lhs(8,0)=DN(2,0)*clhs4 + DN(2,1)*clhs6 + DN(2,2)*clhs8 + clhs15*clhs365 + clhs19*clhs366 + clhs23*clhs367 + clhs81;
            lhs(8,1)=DN(2,0)*clhs26 + DN(2,1)*clhs28 + DN(2,2)*clhs31 + clhs100 + clhs179 + clhs21*clhs368 + clhs23*clhs369;
            lhs(8,2)=DN(2,0)*clhs36 + DN(2,1)*clhs38 + DN(2,2)*clhs40 + clhs109 + clhs22*clhs368 + clhs23*clhs370 + clhs236;
            lhs(8,3)=DN(0,0)*clhs371 - clhs267;
            lhs(8,4)=DN(2,0)*clhs49 + DN(2,1)*clhs51 + DN(2,2)*clhs53 + clhs23*clhs372 + clhs297 + clhs365*clhs55 + clhs366*clhs57;
            lhs(8,5)=DN(2,0)*clhs60 + DN(2,1)*clhs62 + DN(2,2)*clhs65 + clhs21*clhs373 + clhs23*clhs374 + clhs301 + clhs324;
            lhs(8,6)=DN(2,0)*clhs71 + DN(2,1)*clhs73 + DN(2,2)*clhs75 + clhs22*clhs373 + clhs23*clhs375 + clhs305 + clhs340;
            lhs(8,7)=DN(1,0)*clhs371 - clhs357;
            lhs(8,8)=DN(2,0)*clhs82 + DN(2,1)*clhs84 + DN(2,2)*clhs86 + clhs23*clhs379 + clhs365*clhs88 + clhs366*clhs90 + clhs376*tau2 + clhs378;
            lhs(8,9)=DN(2,0)*clhs93 + DN(2,1)*clhs95 + DN(2,2)*clhs98 + clhs21*clhs382 + clhs21*clhs383 + clhs23*clhs384 + clhs381;
            lhs(8,10)=DN(2,0)*clhs104 + DN(2,1)*clhs106 + DN(2,2)*clhs108 + clhs22*clhs382 + clhs22*clhs383 + clhs23*clhs386 + clhs385;
            lhs(8,11)=DN(2,0)*clhs387;
            lhs(8,12)=DN(2,0)*clhs115 + DN(2,1)*clhs117 + DN(2,2)*clhs119 + clhs121*clhs365 + clhs123*clhs366 + clhs23*clhs392 + clhs391;
            lhs(8,13)=DN(2,0)*clhs126 + DN(2,1)*clhs128 + DN(2,2)*clhs131 + clhs21*clhs396 + clhs23*clhs397 + clhs393 + clhs395;
            lhs(8,14)=DN(2,0)*clhs137 + DN(2,1)*clhs139 + DN(2,2)*clhs141 + clhs22*clhs396 + clhs23*clhs400 + clhs398 + clhs399;
            lhs(8,15)=DN(3,0)*clhs371 - clhs401;
            lhs(9,0)=DN(2,0)*clhs6 + DN(2,1)*clhs145 + DN(2,2)*clhs146 + clhs147*clhs368 + clhs150*clhs367 + clhs182 + clhs92;
            lhs(9,1)=DN(2,0)*clhs28 + DN(2,1)*clhs152 + DN(2,2)*clhs154 + clhs150*clhs369 + clhs155*clhs365 + clhs156*clhs366 + clhs184;
            lhs(9,2)=DN(2,0)*clhs38 + DN(2,1)*clhs159 + DN(2,2)*clhs161 + clhs149*clhs368 + clhs150*clhs370 + clhs194 + clhs239;
            lhs(9,3)=DN(0,1)*clhs371 - clhs268;
            lhs(9,4)=DN(2,0)*clhs51 + DN(2,1)*clhs163 + DN(2,2)*clhs164 + clhs147*clhs373 + clhs150*clhs372 + clhs299 + clhs325;
            lhs(9,5)=DN(2,0)*clhs62 + DN(2,1)*clhs168 + DN(2,2)*clhs170 + clhs150*clhs374 + clhs171*clhs365 + clhs172*clhs366 + clhs327;
            lhs(9,6)=DN(2,0)*clhs73 + DN(2,1)*clhs174 + DN(2,2)*clhs176 + clhs149*clhs373 + clhs150*clhs375 + clhs329 + clhs342;
            lhs(9,7)=DN(1,1)*clhs371 - clhs358;
            lhs(9,8)=DN(2,0)*clhs84 + DN(2,1)*clhs180 + DN(2,2)*clhs181 + clhs147*clhs382 + clhs147*clhs383 + clhs150*clhs379 + clhs381;
            lhs(9,9)=DN(2,0)*clhs95 + DN(2,1)*clhs185 + DN(2,2)*clhs187 + clhs150*clhs384 + clhs188*clhs365 + clhs189*clhs366 + clhs378 + clhs402*tau2;
            lhs(9,10)=DN(2,0)*clhs106 + DN(2,1)*clhs191 + DN(2,2)*clhs193 + clhs149*clhs382 + clhs149*clhs383 + clhs150*clhs386 + clhs404;
            lhs(9,11)=DN(2,1)*clhs387;
            lhs(9,12)=DN(2,0)*clhs117 + DN(2,1)*clhs197 + DN(2,2)*clhs198 + clhs147*clhs396 + clhs150*clhs392 + clhs405 + clhs406;
            lhs(9,13)=DN(2,0)*clhs128 + DN(2,1)*clhs202 + DN(2,2)*clhs204 + clhs150*clhs397 + clhs205*clhs365 + clhs206*clhs366 + clhs408;
            lhs(9,14)=DN(2,0)*clhs139 + DN(2,1)*clhs208 + DN(2,2)*clhs210 + clhs149*clhs396 + clhs150*clhs400 + clhs409 + clhs410;
            lhs(9,15)=DN(3,1)*clhs371 - clhs411;
            lhs(10,0)=DN(2,0)*clhs8 + DN(2,1)*clhs146 + DN(2,2)*clhs213 + clhs103 + clhs214*clhs368 + clhs217*clhs367 + clhs238;
            lhs(10,1)=DN(2,0)*clhs31 + DN(2,1)*clhs154 + DN(2,2)*clhs218 + clhs190 + clhs215*clhs368 + clhs217*clhs369 + clhs241;
            lhs(10,2)=DN(2,0)*clhs40 + DN(2,1)*clhs161 + DN(2,2)*clhs220 + clhs217*clhs370 + clhs221*clhs365 + clhs222*clhs366 + clhs243;
            lhs(10,3)=DN(0,2)*clhs371 - clhs269;
            lhs(10,4)=DN(2,0)*clhs53 + DN(2,1)*clhs164 + DN(2,2)*clhs225 + clhs214*clhs373 + clhs217*clhs372 + clhs304 + clhs341;
            lhs(10,5)=DN(2,0)*clhs65 + DN(2,1)*clhs170 + DN(2,2)*clhs228 + clhs215*clhs373 + clhs217*clhs374 + clhs328 + clhs343;
            lhs(10,6)=DN(2,0)*clhs75 + DN(2,1)*clhs176 + DN(2,2)*clhs232 + clhs217*clhs375 + clhs233*clhs365 + clhs234*clhs366 + clhs345;
            lhs(10,7)=DN(1,2)*clhs371 - clhs359;
            lhs(10,8)=DN(2,0)*clhs86 + DN(2,1)*clhs181 + DN(2,2)*clhs237 + clhs214*clhs382 + clhs214*clhs383 + clhs217*clhs379 + clhs385;
            lhs(10,9)=DN(2,0)*clhs98 + DN(2,1)*clhs187 + DN(2,2)*clhs240 + clhs215*clhs382 + clhs215*clhs383 + clhs217*clhs384 + clhs404;
            lhs(10,10)=DN(2,0)*clhs108 + DN(2,1)*clhs193 + DN(2,2)*clhs244 + clhs217*clhs386 + clhs245*clhs365 + clhs246*clhs366 + clhs378 + clhs412*tau2;
            lhs(10,11)=DN(2,2)*clhs387;
            lhs(10,12)=DN(2,0)*clhs119 + DN(2,1)*clhs198 + DN(2,2)*clhs249 + clhs214*clhs396 + clhs217*clhs392 + clhs414 + clhs415;
            lhs(10,13)=DN(2,0)*clhs131 + DN(2,1)*clhs204 + DN(2,2)*clhs252 + clhs215*clhs396 + clhs217*clhs397 + clhs416 + clhs417;
            lhs(10,14)=DN(2,0)*clhs141 + DN(2,1)*clhs210 + DN(2,2)*clhs256 + clhs217*clhs400 + clhs257*clhs365 + clhs258*clhs366 + clhs419;
            lhs(10,15)=DN(3,2)*clhs371 - clhs420;
            lhs(11,0)=clhs111 + clhs147*clhs369 + clhs19*clhs421 + clhs214*clhs370;
            lhs(11,1)=clhs156*clhs422 + clhs195 + clhs21*clhs367 + clhs215*clhs370;
            lhs(11,2)=clhs149*clhs369 + clhs22*clhs367 + clhs222*clhs423 + clhs247;
            lhs(11,3)=clhs270;
            lhs(11,4)=clhs147*clhs374 + clhs214*clhs375 + clhs307 + clhs421*clhs57;
            lhs(11,5)=clhs172*clhs422 + clhs21*clhs372 + clhs215*clhs375 + clhs330;
            lhs(11,6)=clhs149*clhs374 + clhs22*clhs372 + clhs234*clhs423 + clhs346;
            lhs(11,7)=clhs360;
            lhs(11,8)=DN(2,0)*N[2] + clhs147*clhs384 + clhs214*clhs386 + clhs421*clhs90;
            lhs(11,9)=DN(2,1)*N[2] + clhs189*clhs422 + clhs21*clhs379 + clhs215*clhs386;
            lhs(11,10)=DN(2,2)*N[2] + clhs149*clhs384 + clhs22*clhs379 + clhs246*clhs423;
            lhs(11,11)=tau1*(clhs376 + clhs402 + clhs412);
            lhs(11,12)=clhs123*clhs421 + clhs147*clhs397 + clhs214*clhs400 + clhs424;
            lhs(11,13)=clhs206*clhs422 + clhs21*clhs392 + clhs215*clhs400 + clhs425;
            lhs(11,14)=clhs149*clhs397 + clhs22*clhs392 + clhs258*clhs423 + clhs426;
            lhs(11,15)=clhs427;
            lhs(12,0)=DN(3,0)*clhs4 + DN(3,1)*clhs6 + DN(3,2)*clhs8 + clhs114 + clhs15*clhs389 + clhs19*clhs428 + clhs23*clhs429;
            lhs(12,1)=DN(3,0)*clhs26 + DN(3,1)*clhs28 + DN(3,2)*clhs31 + clhs133 + clhs196 + clhs21*clhs430 + clhs23*clhs431;
            lhs(12,2)=DN(3,0)*clhs36 + DN(3,1)*clhs38 + DN(3,2)*clhs40 + clhs142 + clhs22*clhs430 + clhs23*clhs432 + clhs248;
            lhs(12,3)=DN(0,0)*clhs433 - clhs271;
            lhs(12,4)=DN(3,0)*clhs49 + DN(3,1)*clhs51 + DN(3,2)*clhs53 + clhs23*clhs434 + clhs310 + clhs389*clhs55 + clhs428*clhs57;
            lhs(12,5)=DN(3,0)*clhs60 + DN(3,1)*clhs62 + DN(3,2)*clhs65 + clhs21*clhs435 + clhs23*clhs436 + clhs314 + clhs331;
            lhs(12,6)=DN(3,0)*clhs71 + DN(3,1)*clhs73 + DN(3,2)*clhs75 + clhs22*clhs435 + clhs23*clhs437 + clhs318 + clhs347;
            lhs(12,7)=DN(1,0)*clhs433 - clhs361;
            lhs(12,8)=DN(3,0)*clhs82 + DN(3,1)*clhs84 + DN(3,2)*clhs86 + clhs23*clhs438 + clhs389*clhs88 + clhs391 + clhs428*clhs90;
            lhs(12,9)=DN(3,0)*clhs93 + DN(3,1)*clhs95 + DN(3,2)*clhs98 + clhs21*clhs439 + clhs23*clhs440 + clhs395 + clhs405;
            lhs(12,10)=DN(3,0)*clhs104 + DN(3,1)*clhs106 + DN(3,2)*clhs108 + clhs22*clhs439 + clhs23*clhs441 + clhs399 + clhs414;
            lhs(12,11)=DN(2,0)*clhs433 - clhs424;
            lhs(12,12)=DN(3,0)*clhs115 + DN(3,1)*clhs117 + DN(3,2)*clhs119 + clhs121*clhs389 + clhs123*clhs428 + clhs23*clhs445 + clhs442*tau2 + clhs444;
            lhs(12,13)=DN(3,0)*clhs126 + DN(3,1)*clhs128 + DN(3,2)*clhs131 + clhs21*clhs448 + clhs21*clhs449 + clhs23*clhs450 + clhs447;
            lhs(12,14)=DN(3,0)*clhs137 + DN(3,1)*clhs139 + DN(3,2)*clhs141 + clhs22*clhs448 + clhs22*clhs449 + clhs23*clhs452 + clhs451;
            lhs(12,15)=DN(3,0)*clhs453;
            lhs(13,0)=DN(3,0)*clhs6 + DN(3,1)*clhs145 + DN(3,2)*clhs146 + clhs125 + clhs147*clhs430 + clhs150*clhs429 + clhs199;
            lhs(13,1)=DN(3,0)*clhs28 + DN(3,1)*clhs152 + DN(3,2)*clhs154 + clhs150*clhs431 + clhs155*clhs389 + clhs156*clhs428 + clhs201;
            lhs(13,2)=DN(3,0)*clhs38 + DN(3,1)*clhs159 + DN(3,2)*clhs161 + clhs149*clhs430 + clhs150*clhs432 + clhs211 + clhs251;
            lhs(13,3)=DN(0,1)*clhs433 - clhs272;
            lhs(13,4)=DN(3,0)*clhs51 + DN(3,1)*clhs163 + DN(3,2)*clhs164 + clhs147*clhs435 + clhs150*clhs434 + clhs312 + clhs332;
            lhs(13,5)=DN(3,0)*clhs62 + DN(3,1)*clhs168 + DN(3,2)*clhs170 + clhs150*clhs436 + clhs171*clhs389 + clhs172*clhs428 + clhs334;
            lhs(13,6)=DN(3,0)*clhs73 + DN(3,1)*clhs174 + DN(3,2)*clhs176 + clhs149*clhs435 + clhs150*clhs437 + clhs336 + clhs349;
            lhs(13,7)=DN(1,1)*clhs433 - clhs362;
            lhs(13,8)=DN(3,0)*clhs84 + DN(3,1)*clhs180 + DN(3,2)*clhs181 + clhs147*clhs439 + clhs150*clhs438 + clhs393 + clhs406;
            lhs(13,9)=DN(3,0)*clhs95 + DN(3,1)*clhs185 + DN(3,2)*clhs187 + clhs150*clhs440 + clhs188*clhs389 + clhs189*clhs428 + clhs408;
            lhs(13,10)=DN(3,0)*clhs106 + DN(3,1)*clhs191 + DN(3,2)*clhs193 + clhs149*clhs439 + clhs150*clhs441 + clhs410 + clhs416;
            lhs(13,11)=DN(2,1)*clhs433 - clhs425;
            lhs(13,12)=DN(3,0)*clhs117 + DN(3,1)*clhs197 + DN(3,2)*clhs198 + clhs147*clhs448 + clhs147*clhs449 + clhs150*clhs445 + clhs447;
            lhs(13,13)=DN(3,0)*clhs128 + DN(3,1)*clhs202 + DN(3,2)*clhs204 + clhs150*clhs450 + clhs205*clhs389 + clhs206*clhs428 + clhs444 + clhs454*tau2;
            lhs(13,14)=DN(3,0)*clhs139 + DN(3,1)*clhs208 + DN(3,2)*clhs210 + clhs149*clhs448 + clhs149*clhs449 + clhs150*clhs452 + clhs455;
            lhs(13,15)=DN(3,1)*clhs453;
            lhs(14,0)=DN(3,0)*clhs8 + DN(3,1)*clhs146 + DN(3,2)*clhs213 + clhs136 + clhs214*clhs430 + clhs217*clhs429 + clhs250;
            lhs(14,1)=DN(3,0)*clhs31 + DN(3,1)*clhs154 + DN(3,2)*clhs218 + clhs207 + clhs215*clhs430 + clhs217*clhs431 + clhs253;
            lhs(14,2)=DN(3,0)*clhs40 + DN(3,1)*clhs161 + DN(3,2)*clhs220 + clhs217*clhs432 + clhs221*clhs389 + clhs222*clhs428 + clhs255;
            lhs(14,3)=DN(0,2)*clhs433 - clhs273;
            lhs(14,4)=DN(3,0)*clhs53 + DN(3,1)*clhs164 + DN(3,2)*clhs225 + clhs214*clhs435 + clhs217*clhs434 + clhs317 + clhs348;
            lhs(14,5)=DN(3,0)*clhs65 + DN(3,1)*clhs170 + DN(3,2)*clhs228 + clhs215*clhs435 + clhs217*clhs436 + clhs335 + clhs350;
            lhs(14,6)=DN(3,0)*clhs75 + DN(3,1)*clhs176 + DN(3,2)*clhs232 + clhs217*clhs437 + clhs233*clhs389 + clhs234*clhs428 + clhs352;
            lhs(14,7)=DN(1,2)*clhs433 - clhs363;
            lhs(14,8)=DN(3,0)*clhs86 + DN(3,1)*clhs181 + DN(3,2)*clhs237 + clhs214*clhs439 + clhs217*clhs438 + clhs398 + clhs415;
            lhs(14,9)=DN(3,0)*clhs98 + DN(3,1)*clhs187 + DN(3,2)*clhs240 + clhs215*clhs439 + clhs217*clhs440 + clhs409 + clhs417;
            lhs(14,10)=DN(3,0)*clhs108 + DN(3,1)*clhs193 + DN(3,2)*clhs244 + clhs217*clhs441 + clhs245*clhs389 + clhs246*clhs428 + clhs419;
            lhs(14,11)=DN(2,2)*clhs433 - clhs426;
            lhs(14,12)=DN(3,0)*clhs119 + DN(3,1)*clhs198 + DN(3,2)*clhs249 + clhs214*clhs448 + clhs214*clhs449 + clhs217*clhs445 + clhs451;
            lhs(14,13)=DN(3,0)*clhs131 + DN(3,1)*clhs204 + DN(3,2)*clhs252 + clhs215*clhs448 + clhs215*clhs449 + clhs217*clhs450 + clhs455;
            lhs(14,14)=DN(3,0)*clhs141 + DN(3,1)*clhs210 + DN(3,2)*clhs256 + clhs217*clhs452 + clhs257*clhs389 + clhs258*clhs428 + clhs444 + clhs456*tau2;
            lhs(14,15)=DN(3,2)*clhs453;
            lhs(15,0)=clhs144 + clhs147*clhs431 + clhs19*clhs457 + clhs214*clhs432;
            lhs(15,1)=clhs156*clhs458 + clhs21*clhs429 + clhs212 + clhs215*clhs432;
            lhs(15,2)=clhs149*clhs431 + clhs22*clhs429 + clhs222*clhs459 + clhs259;
            lhs(15,3)=clhs274;
            lhs(15,4)=clhs147*clhs436 + clhs214*clhs437 + clhs320 + clhs457*clhs57;
            lhs(15,5)=clhs172*clhs458 + clhs21*clhs434 + clhs215*clhs437 + clhs337;
            lhs(15,6)=clhs149*clhs436 + clhs22*clhs434 + clhs234*clhs459 + clhs353;
            lhs(15,7)=clhs364;
            lhs(15,8)=clhs147*clhs440 + clhs214*clhs441 + clhs401 + clhs457*clhs90;
            lhs(15,9)=clhs189*clhs458 + clhs21*clhs438 + clhs215*clhs441 + clhs411;
            lhs(15,10)=clhs149*clhs440 + clhs22*clhs438 + clhs246*clhs459 + clhs420;
            lhs(15,11)=clhs427;
            lhs(15,12)=DN(3,0)*N[3] + clhs123*clhs457 + clhs147*clhs450 + clhs214*clhs452;
            lhs(15,13)=DN(3,1)*N[3] + clhs206*clhs458 + clhs21*clhs445 + clhs215*clhs452;
            lhs(15,14)=DN(3,2)*N[3] + clhs149*clhs450 + clhs22*clhs445 + clhs258*clhs459;
            lhs(15,15)=tau1*(clhs442 + clhs454 + clhs456);


}


template<>
void NavierStokes<2>::ComputeGaussPointLHSContribution(bounded_matrix<double,9,9>& lhs, const ElementDataStruct& data)
{
    const int nnodes = 3;
    const int dim = 2;
    //~ const int strain_size = 3;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size
    const double c = data.c;                                // Wave velocity

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& delta_t = data.delta_t;
    const double& dyn_tau_coeff = data.dyn_tau_coeff;

    const bounded_matrix<double,nnodes,dim>& v = data.v;
    const bounded_matrix<double,nnodes,dim>& vn = data.vn;
    const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
    const bounded_matrix<double,nnodes,dim>& vmesh = data.vmesh;
    const bounded_matrix<double,nnodes,dim>& vconv = v - vmesh;
    const bounded_matrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& p = data.p;
    const array_1d<double,nnodes>& pn = data.pn;
    const array_1d<double,nnodes>& pnn = data.pnn;
    //~ const array_1d<double,strain_size>& stress = data.stress;

    // Get constitutive matrix
    const Matrix& C = data.C;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    const array_1d<double,dim> vconv_gauss = prod(trans(vconv), N);

    const double vconv_norm = norm_2(vconv_gauss);

    // Stabilization parameters
    const double c1 = 4.0;
    const double c2 = 2.0;
    const double tau1 = 1.0/((rho*dyn_tau_coeff)/delta_t + (c2*rho*vconv_norm)/h + (c1*mu)/(h*h));
    const double tau2 = (h*h)/(c1*tau1);

    const double clhs0 =             pow(DN(0,0), 2);
const double clhs1 =             pow(N[0], 2);
const double clhs2 =             bdf0*rho;
const double clhs3 =             clhs1*clhs2;
const double clhs4 =             C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
const double clhs5 =             C(0,2)*DN(0,0);
const double clhs6 =             C(2,2)*DN(0,1) + clhs5;
const double clhs7 =             N[0]*rho;
const double clhs8 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double clhs9 =             N[0]*(v(0,0) - vmesh(0,0)) + N[1]*(v(1,0) - vmesh(1,0)) + N[2]*(v(2,0) - vmesh(2,0));
const double clhs10 =             N[0]*(v(0,1) - vmesh(0,1)) + N[1]*(v(1,1) - vmesh(1,1)) + N[2]*(v(2,1) - vmesh(2,1));
const double clhs11 =             DN(0,0)*clhs9 + DN(0,1)*clhs10;
const double clhs12 =             N[0]*clhs8 + clhs11;
const double clhs13 =             pow(rho, 2);
const double clhs14 =             clhs11*clhs13*tau1;
const double clhs15 =             N[0]*bdf0;
const double clhs16 =             clhs12 + clhs15;
const double clhs17 =             DN(0,0)*N[0]*rho*tau1;
const double clhs18 =             DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0);
const double clhs19 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + clhs10*clhs18 + clhs8*clhs9);
const double clhs20 =             DN(0,0)*tau2;
const double clhs21 =             DN(0,1)*clhs20;
const double clhs22 =             C(0,1)*DN(0,1) + clhs5;
const double clhs23 =             C(1,2)*DN(0,1);
const double clhs24 =             C(2,2)*DN(0,0) + clhs23;
const double clhs25 =             clhs1*rho;
const double clhs26 =             N[0]*clhs11*clhs13*tau1;
const double clhs27 =             DN(0,1)*N[0]*rho*tau1;
const double clhs28 =             rho*tau1;
const double clhs29 =             clhs11*clhs28;
const double clhs30 =             -N[0] + clhs29;
const double clhs31 =             DN(0,0)*DN(1,0);
const double clhs32 =             N[0]*bdf0*rho;
const double clhs33 =             N[1]*clhs32;
const double clhs34 =             clhs31*tau2 + clhs33;
const double clhs35 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
const double clhs36 =             C(0,2)*DN(1,0);
const double clhs37 =             C(2,2)*DN(1,1) + clhs36;
const double clhs38 =             DN(1,0)*clhs9 + DN(1,1)*clhs10;
const double clhs39 =             N[1]*clhs8 + clhs38;
const double clhs40 =             N[1]*bdf0;
const double clhs41 =             clhs39 + clhs40;
const double clhs42 =             DN(0,0)*N[1]*rho*tau1;
const double clhs43 =             DN(1,1)*clhs20;
const double clhs44 =             C(0,1)*DN(1,1) + clhs36;
const double clhs45 =             C(1,2)*DN(1,1);
const double clhs46 =             C(2,2)*DN(1,0) + clhs45;
const double clhs47 =             N[0]*N[1]*rho;
const double clhs48 =             clhs18*clhs47;
const double clhs49 =             N[1]*clhs11*clhs13*tau1;
const double clhs50 =             DN(0,1)*N[1]*rho*tau1;
const double clhs51 =             DN(0,0)*N[1];
const double clhs52 =             DN(0,0)*DN(2,0);
const double clhs53 =             N[2]*clhs32;
const double clhs54 =             clhs52*tau2 + clhs53;
const double clhs55 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
const double clhs56 =             C(0,2)*DN(2,0);
const double clhs57 =             C(2,2)*DN(2,1) + clhs56;
const double clhs58 =             DN(2,0)*clhs9 + DN(2,1)*clhs10;
const double clhs59 =             N[2]*clhs8 + clhs58;
const double clhs60 =             N[2]*bdf0;
const double clhs61 =             clhs59 + clhs60;
const double clhs62 =             DN(0,0)*N[2]*rho*tau1;
const double clhs63 =             DN(2,1)*clhs20;
const double clhs64 =             C(0,1)*DN(2,1) + clhs56;
const double clhs65 =             C(1,2)*DN(2,1);
const double clhs66 =             C(2,2)*DN(2,0) + clhs65;
const double clhs67 =             N[0]*N[2]*rho;
const double clhs68 =             clhs18*clhs67;
const double clhs69 =             N[2]*clhs11*clhs13*tau1;
const double clhs70 =             DN(0,1)*N[2]*rho*tau1;
const double clhs71 =             DN(0,0)*N[2];
const double clhs72 =             C(0,1)*DN(0,0) + clhs23;
const double clhs73 =             DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1);
const double clhs74 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double clhs75 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + clhs10*clhs74 + clhs73*clhs9);
const double clhs76 =             pow(DN(0,1), 2);
const double clhs77 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
const double clhs78 =             N[0]*clhs74 + clhs11;
const double clhs79 =             clhs15 + clhs78;
const double clhs80 =             DN(0,1)*tau2;
const double clhs81 =             DN(1,0)*clhs80;
const double clhs82 =             C(0,1)*DN(1,0) + clhs45;
const double clhs83 =             clhs47*clhs73;
const double clhs84 =             DN(0,1)*DN(1,1);
const double clhs85 =             clhs33 + clhs84*tau2;
const double clhs86 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
const double clhs87 =             N[1]*clhs74 + clhs38;
const double clhs88 =             clhs40 + clhs87;
const double clhs89 =             DN(0,1)*N[1];
const double clhs90 =             DN(2,0)*clhs80;
const double clhs91 =             C(0,1)*DN(2,0) + clhs65;
const double clhs92 =             clhs67*clhs73;
const double clhs93 =             DN(0,1)*DN(2,1);
const double clhs94 =             clhs53 + clhs93*tau2;
const double clhs95 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
const double clhs96 =             N[2]*clhs74 + clhs58;
const double clhs97 =             clhs60 + clhs96;
const double clhs98 =             DN(0,1)*N[2];
const double clhs99 =             DN(0,0)*N[0];
const double clhs100 =             DN(0,1)*N[0];
const double clhs101 =             clhs73*rho*tau1;
const double clhs102 =             DN(0,0)*rho*tau1;
const double clhs103 =             clhs18*rho*tau1;
const double clhs104 =             DN(0,1)*rho*tau1;
const double clhs105 =             DN(1,0)*N[0];
const double clhs106 =             DN(1,1)*N[0];
const double clhs107 =             tau1*(clhs31 + clhs84);
const double clhs108 =             DN(2,0)*N[0];
const double clhs109 =             DN(2,1)*N[0];
const double clhs110 =             tau1*(clhs52 + clhs93);
const double clhs111 =             N[1]*rho;
const double clhs112 =             clhs13*clhs38*tau1;
const double clhs113 =             DN(1,0)*N[0]*rho*tau1;
const double clhs114 =             N[0]*clhs13*clhs38*tau1;
const double clhs115 =             DN(1,1)*N[0]*rho*tau1;
const double clhs116 =             clhs28*clhs38;
const double clhs117 =             pow(DN(1,0), 2);
const double clhs118 =             pow(N[1], 2);
const double clhs119 =             clhs118*clhs2;
const double clhs120 =             DN(1,0)*N[1]*rho*tau1;
const double clhs121 =             DN(1,0)*tau2;
const double clhs122 =             DN(1,1)*clhs121;
const double clhs123 =             clhs118*rho;
const double clhs124 =             N[1]*clhs13*clhs38*tau1;
const double clhs125 =             DN(1,1)*N[1]*rho*tau1;
const double clhs126 =             -N[1] + clhs116;
const double clhs127 =             DN(1,0)*DN(2,0);
const double clhs128 =             N[2]*rho;
const double clhs129 =             clhs128*clhs40;
const double clhs130 =             clhs127*tau2 + clhs129;
const double clhs131 =             DN(1,0)*N[2]*rho*tau1;
const double clhs132 =             DN(2,1)*clhs121;
const double clhs133 =             N[1]*N[2]*rho;
const double clhs134 =             clhs133*clhs18;
const double clhs135 =             N[2]*clhs13*clhs38*tau1;
const double clhs136 =             DN(1,1)*N[2]*rho*tau1;
const double clhs137 =             DN(1,0)*N[2];
const double clhs138 =             pow(DN(1,1), 2);
const double clhs139 =             DN(2,0)*tau2;
const double clhs140 =             DN(1,1)*clhs139;
const double clhs141 =             clhs133*clhs73;
const double clhs142 =             DN(1,1)*DN(2,1);
const double clhs143 =             clhs129 + clhs142*tau2;
const double clhs144 =             DN(1,1)*N[2];
const double clhs145 =             DN(1,0)*rho*tau1;
const double clhs146 =             DN(1,1)*rho*tau1;
const double clhs147 =             DN(1,0)*N[1];
const double clhs148 =             DN(1,1)*N[1];
const double clhs149 =             DN(2,0)*N[1];
const double clhs150 =             DN(2,1)*N[1];
const double clhs151 =             tau1*(clhs127 + clhs142);
const double clhs152 =             clhs13*clhs58*tau1;
const double clhs153 =             DN(2,0)*N[0]*rho*tau1;
const double clhs154 =             N[0]*clhs13*clhs58*tau1;
const double clhs155 =             DN(2,1)*N[0]*rho*tau1;
const double clhs156 =             clhs28*clhs58;
const double clhs157 =             DN(2,0)*N[1]*rho*tau1;
const double clhs158 =             N[1]*clhs13*clhs58*tau1;
const double clhs159 =             DN(2,1)*N[1]*rho*tau1;
const double clhs160 =             pow(DN(2,0), 2);
const double clhs161 =             pow(N[2], 2);
const double clhs162 =             clhs161*clhs2;
const double clhs163 =             DN(2,0)*N[2]*rho*tau1;
const double clhs164 =             DN(2,1)*clhs139;
const double clhs165 =             clhs161*rho;
const double clhs166 =             N[2]*clhs13*clhs58*tau1;
const double clhs167 =             DN(2,1)*N[2]*rho*tau1;
const double clhs168 =             -N[2] + clhs156;
const double clhs169 =             pow(DN(2,1), 2);
const double clhs170 =             DN(2,0)*rho*tau1;
const double clhs171 =             DN(2,1)*rho*tau1;
const double clhs172 =             DN(2,0)*N[2];
const double clhs173 =             DN(2,1)*N[2];
            lhs(0,0)=DN(0,0)*clhs4 + DN(0,1)*clhs6 + clhs0*tau2 + clhs12*clhs7 + clhs14*clhs16 + clhs17*clhs19 + clhs3;
            lhs(0,1)=DN(0,0)*clhs22 + DN(0,1)*clhs24 + clhs18*clhs25 + clhs18*clhs26 + clhs19*clhs27 + clhs21;
            lhs(0,2)=DN(0,0)*clhs30;
            lhs(0,3)=DN(0,0)*clhs35 + DN(0,1)*clhs37 + clhs14*clhs41 + clhs19*clhs42 + clhs34 + clhs39*clhs7;
            lhs(0,4)=DN(0,0)*clhs44 + DN(0,1)*clhs46 + clhs18*clhs49 + clhs19*clhs50 + clhs43 + clhs48;
            lhs(0,5)=DN(1,0)*clhs29 - clhs51;
            lhs(0,6)=DN(0,0)*clhs55 + DN(0,1)*clhs57 + clhs14*clhs61 + clhs19*clhs62 + clhs54 + clhs59*clhs7;
            lhs(0,7)=DN(0,0)*clhs64 + DN(0,1)*clhs66 + clhs18*clhs69 + clhs19*clhs70 + clhs63 + clhs68;
            lhs(0,8)=DN(2,0)*clhs29 - clhs71;
            lhs(1,0)=DN(0,0)*clhs6 + DN(0,1)*clhs72 + clhs17*clhs75 + clhs21 + clhs25*clhs73 + clhs26*clhs73;
            lhs(1,1)=DN(0,0)*clhs24 + DN(0,1)*clhs77 + clhs14*clhs79 + clhs27*clhs75 + clhs3 + clhs7*clhs78 + clhs76*tau2;
            lhs(1,2)=DN(0,1)*clhs30;
            lhs(1,3)=DN(0,0)*clhs37 + DN(0,1)*clhs82 + clhs42*clhs75 + clhs49*clhs73 + clhs81 + clhs83;
            lhs(1,4)=DN(0,0)*clhs46 + DN(0,1)*clhs86 + clhs14*clhs88 + clhs50*clhs75 + clhs7*clhs87 + clhs85;
            lhs(1,5)=DN(1,1)*clhs29 - clhs89;
            lhs(1,6)=DN(0,0)*clhs57 + DN(0,1)*clhs91 + clhs62*clhs75 + clhs69*clhs73 + clhs90 + clhs92;
            lhs(1,7)=DN(0,0)*clhs66 + DN(0,1)*clhs95 + clhs14*clhs97 + clhs7*clhs96 + clhs70*clhs75 + clhs94;
            lhs(1,8)=DN(2,1)*clhs29 - clhs98;
            lhs(2,0)=clhs100*clhs101 + clhs102*clhs16 + clhs99;
            lhs(2,1)=clhs100 + clhs103*clhs99 + clhs104*clhs79;
            lhs(2,2)=tau1*(clhs0 + clhs76);
            lhs(2,3)=clhs101*clhs89 + clhs102*clhs41 + clhs105;
            lhs(2,4)=clhs103*clhs51 + clhs104*clhs88 + clhs106;
            lhs(2,5)=clhs107;
            lhs(2,6)=clhs101*clhs98 + clhs102*clhs61 + clhs108;
            lhs(2,7)=clhs103*clhs71 + clhs104*clhs97 + clhs109;
            lhs(2,8)=clhs110;
            lhs(3,0)=DN(1,0)*clhs4 + DN(1,1)*clhs6 + clhs111*clhs12 + clhs112*clhs16 + clhs113*clhs19 + clhs34;
            lhs(3,1)=DN(1,0)*clhs22 + DN(1,1)*clhs24 + clhs114*clhs18 + clhs115*clhs19 + clhs48 + clhs81;
            lhs(3,2)=DN(0,0)*clhs116 - clhs105;
            lhs(3,3)=DN(1,0)*clhs35 + DN(1,1)*clhs37 + clhs111*clhs39 + clhs112*clhs41 + clhs117*tau2 + clhs119 + clhs120*clhs19;
            lhs(3,4)=DN(1,0)*clhs44 + DN(1,1)*clhs46 + clhs122 + clhs123*clhs18 + clhs124*clhs18 + clhs125*clhs19;
            lhs(3,5)=DN(1,0)*clhs126;
            lhs(3,6)=DN(1,0)*clhs55 + DN(1,1)*clhs57 + clhs111*clhs59 + clhs112*clhs61 + clhs130 + clhs131*clhs19;
            lhs(3,7)=DN(1,0)*clhs64 + DN(1,1)*clhs66 + clhs132 + clhs134 + clhs135*clhs18 + clhs136*clhs19;
            lhs(3,8)=DN(2,0)*clhs116 - clhs137;
            lhs(4,0)=DN(1,0)*clhs6 + DN(1,1)*clhs72 + clhs113*clhs75 + clhs114*clhs73 + clhs43 + clhs83;
            lhs(4,1)=DN(1,0)*clhs24 + DN(1,1)*clhs77 + clhs111*clhs78 + clhs112*clhs79 + clhs115*clhs75 + clhs85;
            lhs(4,2)=DN(0,1)*clhs116 - clhs106;
            lhs(4,3)=DN(1,0)*clhs37 + DN(1,1)*clhs82 + clhs120*clhs75 + clhs122 + clhs123*clhs73 + clhs124*clhs73;
            lhs(4,4)=DN(1,0)*clhs46 + DN(1,1)*clhs86 + clhs111*clhs87 + clhs112*clhs88 + clhs119 + clhs125*clhs75 + clhs138*tau2;
            lhs(4,5)=DN(1,1)*clhs126;
            lhs(4,6)=DN(1,0)*clhs57 + DN(1,1)*clhs91 + clhs131*clhs75 + clhs135*clhs73 + clhs140 + clhs141;
            lhs(4,7)=DN(1,0)*clhs66 + DN(1,1)*clhs95 + clhs111*clhs96 + clhs112*clhs97 + clhs136*clhs75 + clhs143;
            lhs(4,8)=DN(2,1)*clhs116 - clhs144;
            lhs(5,0)=clhs101*clhs106 + clhs145*clhs16 + clhs51;
            lhs(5,1)=clhs103*clhs105 + clhs146*clhs79 + clhs89;
            lhs(5,2)=clhs107;
            lhs(5,3)=clhs101*clhs148 + clhs145*clhs41 + clhs147;
            lhs(5,4)=clhs103*clhs147 + clhs146*clhs88 + clhs148;
            lhs(5,5)=tau1*(clhs117 + clhs138);
            lhs(5,6)=clhs101*clhs144 + clhs145*clhs61 + clhs149;
            lhs(5,7)=clhs103*clhs137 + clhs146*clhs97 + clhs150;
            lhs(5,8)=clhs151;
            lhs(6,0)=DN(2,0)*clhs4 + DN(2,1)*clhs6 + clhs12*clhs128 + clhs152*clhs16 + clhs153*clhs19 + clhs54;
            lhs(6,1)=DN(2,0)*clhs22 + DN(2,1)*clhs24 + clhs154*clhs18 + clhs155*clhs19 + clhs68 + clhs90;
            lhs(6,2)=DN(0,0)*clhs156 - clhs108;
            lhs(6,3)=DN(2,0)*clhs35 + DN(2,1)*clhs37 + clhs128*clhs39 + clhs130 + clhs152*clhs41 + clhs157*clhs19;
            lhs(6,4)=DN(2,0)*clhs44 + DN(2,1)*clhs46 + clhs134 + clhs140 + clhs158*clhs18 + clhs159*clhs19;
            lhs(6,5)=DN(1,0)*clhs156 - clhs149;
            lhs(6,6)=DN(2,0)*clhs55 + DN(2,1)*clhs57 + clhs128*clhs59 + clhs152*clhs61 + clhs160*tau2 + clhs162 + clhs163*clhs19;
            lhs(6,7)=DN(2,0)*clhs64 + DN(2,1)*clhs66 + clhs164 + clhs165*clhs18 + clhs166*clhs18 + clhs167*clhs19;
            lhs(6,8)=DN(2,0)*clhs168;
            lhs(7,0)=DN(2,0)*clhs6 + DN(2,1)*clhs72 + clhs153*clhs75 + clhs154*clhs73 + clhs63 + clhs92;
            lhs(7,1)=DN(2,0)*clhs24 + DN(2,1)*clhs77 + clhs128*clhs78 + clhs152*clhs79 + clhs155*clhs75 + clhs94;
            lhs(7,2)=DN(0,1)*clhs156 - clhs109;
            lhs(7,3)=DN(2,0)*clhs37 + DN(2,1)*clhs82 + clhs132 + clhs141 + clhs157*clhs75 + clhs158*clhs73;
            lhs(7,4)=DN(2,0)*clhs46 + DN(2,1)*clhs86 + clhs128*clhs87 + clhs143 + clhs152*clhs88 + clhs159*clhs75;
            lhs(7,5)=DN(1,1)*clhs156 - clhs150;
            lhs(7,6)=DN(2,0)*clhs57 + DN(2,1)*clhs91 + clhs163*clhs75 + clhs164 + clhs165*clhs73 + clhs166*clhs73;
            lhs(7,7)=DN(2,0)*clhs66 + DN(2,1)*clhs95 + clhs128*clhs96 + clhs152*clhs97 + clhs162 + clhs167*clhs75 + clhs169*tau2;
            lhs(7,8)=DN(2,1)*clhs168;
            lhs(8,0)=clhs101*clhs109 + clhs16*clhs170 + clhs71;
            lhs(8,1)=clhs103*clhs108 + clhs171*clhs79 + clhs98;
            lhs(8,2)=clhs110;
            lhs(8,3)=clhs101*clhs150 + clhs137 + clhs170*clhs41;
            lhs(8,4)=clhs103*clhs149 + clhs144 + clhs171*clhs88;
            lhs(8,5)=clhs151;
            lhs(8,6)=clhs101*clhs173 + clhs170*clhs61 + clhs172;
            lhs(8,7)=clhs103*clhs172 + clhs171*clhs97 + clhs173;
            lhs(8,8)=tau1*(clhs160 + clhs169);


}


template<>
void NavierStokes<3>::ComputeGaussPointRHSContribution(array_1d<double,16>& rhs, const ElementDataStruct& data)
{
    const int nnodes = 4;
    const int dim = 3;
    const int strain_size = 6;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size
    const double c = data.c;                                // Wave velocity

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& delta_t = data.delta_t;
    const double& dyn_tau_coeff = data.dyn_tau_coeff;

    const bounded_matrix<double,nnodes,dim>& v = data.v;
    const bounded_matrix<double,nnodes,dim>& vn = data.vn;
    const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
    const bounded_matrix<double,nnodes,dim>& vmesh = data.vmesh;
    const bounded_matrix<double,nnodes,dim>& vconv = v - vmesh;
    const bounded_matrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& p = data.p;
    const array_1d<double,nnodes>& pn = data.pn;
    const array_1d<double,nnodes>& pnn = data.pnn;
    const array_1d<double,strain_size>& stress = data.stress;

    // Get constitutive matrix
    // const Matrix& C = data.C;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,dim> f_gauss = prod(trans(f), N);
    const array_1d<double,dim> grad_p = prod(trans(DN), p);
    const array_1d<double,dim> vconv_gauss = prod(trans(vconv), N);
    //~ const double p_gauss = inner_prod(N,p);

    const double vconv_norm = norm_2(vconv_gauss);

    //~ array_1d<double,dim> accel_gauss = bdf0*v_gauss;
    //~ noalias(accel_gauss) += bdf1*prod(trans(vn), N);
    //~ noalias(accel_gauss) += bdf2*prod(trans(vnn), N);

    // Stabilization parameters
    const double c1 = 4.0;
    const double c2 = 2.0;
    const double tau1 = 1.0/((rho*dyn_tau_coeff)/delta_t + (c2*rho*vconv_norm)/h + (c1*mu)/(h*h));
    const double tau2 = (h*h)/(c1*tau1);

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0));
const double crhs2 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
const double crhs3 =             DN(0,0)*v(0,0);
const double crhs4 =             DN(0,1)*v(0,1);
const double crhs5 =             DN(1,0)*v(1,0);
const double crhs6 =             DN(1,1)*v(1,1);
const double crhs7 =             DN(2,0)*v(2,0);
const double crhs8 =             DN(2,1)*v(2,1);
const double crhs9 =             DN(3,0)*v(3,0);
const double crhs10 =             DN(3,1)*v(3,1);
const double crhs11 =             crhs10 + crhs2 + crhs3 + crhs4 + crhs5 + crhs6 + crhs7 + crhs8 + crhs9;
const double crhs12 =             crhs11*tau2;
const double crhs13 =             N[0]*rho;
const double crhs14 =             N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0));
const double crhs15 =             N[0]*(v(0,0) - vmesh(0,0)) + N[1]*(v(1,0) - vmesh(1,0)) + N[2]*(v(2,0) - vmesh(2,0)) + N[3]*(v(3,0) - vmesh(3,0));
const double crhs16 =             N[0]*(v(0,1) - vmesh(0,1)) + N[1]*(v(1,1) - vmesh(1,1)) + N[2]*(v(2,1) - vmesh(2,1)) + N[3]*(v(3,1) - vmesh(3,1));
const double crhs17 =             N[0]*(v(0,2) - vmesh(0,2)) + N[1]*(v(1,2) - vmesh(1,2)) + N[2]*(v(2,2) - vmesh(2,2)) + N[3]*(v(3,2) - vmesh(3,2));
const double crhs18 =             crhs15*(crhs3 + crhs5 + crhs7 + crhs9) + crhs16*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs17*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0));
const double crhs19 =             rho*(DN(0,0)*crhs15 + DN(0,1)*crhs16 + DN(0,2)*crhs17);
const double crhs20 =             tau1*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + rho*(crhs14 + crhs18));
const double crhs21 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
const double crhs22 =             N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1));
const double crhs23 =             crhs15*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs16*(crhs10 + crhs4 + crhs6 + crhs8) + crhs17*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1));
const double crhs24 =             tau1*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs21 + rho*(crhs22 + crhs23));
const double crhs25 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
const double crhs26 =             N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2));
const double crhs27 =             crhs15*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs16*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + crhs17*crhs2;
const double crhs28 =             tau1*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs25 + rho*(crhs26 + crhs27));
const double crhs29 =             N[1]*rho;
const double crhs30 =             rho*(DN(1,0)*crhs15 + DN(1,1)*crhs16 + DN(1,2)*crhs17);
const double crhs31 =             N[2]*rho;
const double crhs32 =             rho*(DN(2,0)*crhs15 + DN(2,1)*crhs16 + DN(2,2)*crhs17);
const double crhs33 =             N[3]*rho;
const double crhs34 =             rho*(DN(3,0)*crhs15 + DN(3,1)*crhs16 + DN(3,2)*crhs17);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs12 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - crhs13*crhs14 - crhs13*crhs18 - crhs19*crhs20;
            rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs12 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs21 - crhs13*crhs22 - crhs13*crhs23 - crhs19*crhs24;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs12 - DN(0,2)*stress[2] + N[0]*crhs25 - crhs13*crhs26 - crhs13*crhs27 - crhs19*crhs28;
            rhs[3]=-DN(0,0)*crhs20 - DN(0,1)*crhs24 - DN(0,2)*crhs28 - N[0]*crhs11;
            rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs12 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - crhs14*crhs29 - crhs18*crhs29 - crhs20*crhs30;
            rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs12 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs21 - crhs22*crhs29 - crhs23*crhs29 - crhs24*crhs30;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs12 - DN(1,2)*stress[2] + N[1]*crhs25 - crhs26*crhs29 - crhs27*crhs29 - crhs28*crhs30;
            rhs[7]=-DN(1,0)*crhs20 - DN(1,1)*crhs24 - DN(1,2)*crhs28 - N[1]*crhs11;
            rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs12 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - crhs14*crhs31 - crhs18*crhs31 - crhs20*crhs32;
            rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs12 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs21 - crhs22*crhs31 - crhs23*crhs31 - crhs24*crhs32;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs12 - DN(2,2)*stress[2] + N[2]*crhs25 - crhs26*crhs31 - crhs27*crhs31 - crhs28*crhs32;
            rhs[11]=-DN(2,0)*crhs20 - DN(2,1)*crhs24 - DN(2,2)*crhs28 - N[2]*crhs11;
            rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs12 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - crhs14*crhs33 - crhs18*crhs33 - crhs20*crhs34;
            rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs12 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs21 - crhs22*crhs33 - crhs23*crhs33 - crhs24*crhs34;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs12 - DN(3,2)*stress[2] + N[3]*crhs25 - crhs26*crhs33 - crhs27*crhs33 - crhs28*crhs34;
            rhs[15]=-DN(3,0)*crhs20 - DN(3,1)*crhs24 - DN(3,2)*crhs28 - N[3]*crhs11;

}


template<>
void NavierStokes<2>::ComputeGaussPointRHSContribution(array_1d<double,9>& rhs, const ElementDataStruct& data)
{
    const int nnodes = 3;
    const int dim = 2;
    const int strain_size = 3;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size
    const double c = data.c;                                // Wave velocity

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& delta_t = data.delta_t;
    const double& dyn_tau_coeff = data.dyn_tau_coeff;

    const bounded_matrix<double,nnodes,dim>& v = data.v;
    const bounded_matrix<double,nnodes,dim>& vn = data.vn;
    const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
    const bounded_matrix<double,nnodes,dim>& vmesh = data.vmesh;
    const bounded_matrix<double,nnodes,dim>& vconv = v - vmesh;
    const bounded_matrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& p = data.p;
    const array_1d<double,nnodes>& pn = data.pn;
    const array_1d<double,nnodes>& pnn = data.pnn;
    const array_1d<double,strain_size>& stress = data.stress;

    // Get constitutive matrix
    // const Matrix& C = data.C;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,dim> f_gauss = prod(trans(f), N);
    const array_1d<double,dim> grad_p = prod(trans(DN), p);
    const array_1d<double,dim> vconv_gauss = prod(trans(vconv), N);
    //~ const double p_gauss = inner_prod(N,p);

    const double vconv_norm = norm_2(vconv_gauss);

    //~ array_1d<double,dim> accel_gauss = bdf0*v_gauss;
    //~ noalias(accel_gauss) += bdf1*prod(trans(vn), N);
    //~ noalias(accel_gauss) += bdf2*prod(trans(vnn), N);

    // Stabilization parameters
    const double c1 = 4.0;
    const double c2 = 2.0;
    const double tau1 = 1.0/((rho*dyn_tau_coeff)/delta_t + (c2*rho*vconv_norm)/h + (c1*mu)/(h*h));
    const double tau2 = (h*h)/(c1*tau1);

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
const double crhs1 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
const double crhs2 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
const double crhs3 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
const double crhs4 =             crhs2 + crhs3;
const double crhs5 =             crhs4*tau2;
const double crhs6 =             N[0]*rho;
const double crhs7 =             N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0));
const double crhs8 =             N[0]*(v(0,0) - vmesh(0,0)) + N[1]*(v(1,0) - vmesh(1,0)) + N[2]*(v(2,0) - vmesh(2,0));
const double crhs9 =             N[0]*(v(0,1) - vmesh(0,1)) + N[1]*(v(1,1) - vmesh(1,1)) + N[2]*(v(2,1) - vmesh(2,1));
const double crhs10 =             crhs2*crhs8 + crhs9*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0));
const double crhs11 =             rho*(DN(0,0)*crhs8 + DN(0,1)*crhs9);
const double crhs12 =             tau1*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + rho*(crhs10 + crhs7));
const double crhs13 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
const double crhs14 =             N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1));
const double crhs15 =             crhs3*crhs9 + crhs8*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1));
const double crhs16 =             tau1*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs13 + rho*(crhs14 + crhs15));
const double crhs17 =             N[1]*rho;
const double crhs18 =             rho*(DN(1,0)*crhs8 + DN(1,1)*crhs9);
const double crhs19 =             N[2]*rho;
const double crhs20 =             rho*(DN(2,0)*crhs8 + DN(2,1)*crhs9);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs5 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - crhs10*crhs6 - crhs11*crhs12 - crhs6*crhs7;
            rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs5 - DN(0,1)*stress[1] + N[0]*crhs13 - crhs11*crhs16 - crhs14*crhs6 - crhs15*crhs6;
            rhs[2]=-DN(0,0)*crhs12 - DN(0,1)*crhs16 - N[0]*crhs4;
            rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs5 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - crhs10*crhs17 - crhs12*crhs18 - crhs17*crhs7;
            rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs5 - DN(1,1)*stress[1] + N[1]*crhs13 - crhs14*crhs17 - crhs15*crhs17 - crhs16*crhs18;
            rhs[5]=-DN(1,0)*crhs12 - DN(1,1)*crhs16 - N[1]*crhs4;
            rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs5 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - crhs10*crhs19 - crhs12*crhs20 - crhs19*crhs7;
            rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs5 - DN(2,1)*stress[1] + N[2]*crhs13 - crhs14*crhs19 - crhs15*crhs19 - crhs16*crhs20;
            rhs[8]=-DN(2,0)*crhs12 - DN(2,1)*crhs16 - N[2]*crhs4;

}


template<>
double NavierStokes<3>::SubscaleErrorEstimate(const ElementDataStruct& data)
{
    const int nnodes = 4;
    const int dim = 3;
    // const int strain_size = 3;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size
    const double c = data.c;                                // Wave velocity

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& delta_t = data.delta_t;
    const double& dyn_tau_coeff = data.dyn_tau_coeff;

    const bounded_matrix<double,nnodes,dim>& v = data.v;
    const bounded_matrix<double,nnodes,dim>& vn = data.vn;
    const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
    const bounded_matrix<double,nnodes,dim>& vmesh = data.vmesh;
    const bounded_matrix<double,nnodes,dim>& vconv = v - vmesh;
    const bounded_matrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& p = data.p;
    const array_1d<double,nnodes>& pn = data.pn;
    const array_1d<double,nnodes>& pnn = data.pnn;
    // const array_1d<double,strain_size>& stress = data.stress;

    // // Get constitutive matrix
    // const Matrix& C = data.C;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the error estimator
    array_1d<double,dim> v_s_gauss;
    const array_1d<double,dim> v_gauss = prod(trans(v), N);
    const array_1d<double,dim> f_gauss = prod(trans(f), N);
    const array_1d<double,dim> grad_p = prod(trans(DN), p);
    const array_1d<double,dim> vconv_gauss = prod(trans(vconv), N);

    const double vconv_norm = norm_2(vconv_gauss);

    // Stabilization parameters
    const double c1 = 4.0;
    const double c2 = 2.0;
    const double tau1 = 1.0/((rho*dyn_tau_coeff)/delta_t + (c2*rho*vconv_norm)/h + (c1*mu)/(h*h));
    // const double tau2 = (h*h)/(c1*tau1);

    // Gauss point velocity subscale value computation
    const double cv_s_gauss0 =             N[0]*(v(0,0) - vmesh(0,0)) + N[1]*(v(1,0) - vmesh(1,0)) + N[2]*(v(2,0) - vmesh(2,0)) + N[3]*(v(3,0) - vmesh(3,0));
const double cv_s_gauss1 =             N[0]*(v(0,1) - vmesh(0,1)) + N[1]*(v(1,1) - vmesh(1,1)) + N[2]*(v(2,1) - vmesh(2,1)) + N[3]*(v(3,1) - vmesh(3,1));
const double cv_s_gauss2 =             N[0]*(v(0,2) - vmesh(0,2)) + N[1]*(v(1,2) - vmesh(1,2)) + N[2]*(v(2,2) - vmesh(2,2)) + N[3]*(v(3,2) - vmesh(3,2));
            v_s_gauss[0]=-tau1*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + rho*(-N[0]*f(0,0) + N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) - N[1]*f(1,0) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) - N[2]*f(2,0) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) - N[3]*f(3,0) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)) + cv_s_gauss0*(DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0)) + cv_s_gauss1*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + cv_s_gauss2*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0))));
            v_s_gauss[1]=-tau1*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + rho*(-N[0]*f(0,1) + N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) - N[1]*f(1,1) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) - N[2]*f(2,1) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) - N[3]*f(3,1) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)) + cv_s_gauss0*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + cv_s_gauss1*(DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1)) + cv_s_gauss2*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1))));
            v_s_gauss[2]=-tau1*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + rho*(-N[0]*f(0,2) + N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) - N[1]*f(1,2) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) - N[2]*f(2,2) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) - N[3]*f(3,2) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)) + cv_s_gauss0*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + cv_s_gauss1*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + cv_s_gauss2*(DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2))));


    const double v_gauss_norm = norm_2(v_gauss);
    const double v_s_gauss_norm = norm_2(v_s_gauss);

    return v_s_gauss_norm/v_gauss_norm;
}


template<>
double NavierStokes<2>::SubscaleErrorEstimate(const ElementDataStruct& data)
{
    const int nnodes = 3;
    const int dim = 2;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size
    const double c = data.c;                                // Wave velocity

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& delta_t = data.delta_t;
    const double& dyn_tau_coeff = data.dyn_tau_coeff;

    const bounded_matrix<double,nnodes,dim>& v = data.v;
    const bounded_matrix<double,nnodes,dim>& vn = data.vn;
    const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
    const bounded_matrix<double,nnodes,dim>& vmesh = data.vmesh;
    const bounded_matrix<double,nnodes,dim>& vconv = v - vmesh;
    const bounded_matrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& p = data.p;
    const array_1d<double,nnodes>& pn = data.pn;
    const array_1d<double,nnodes>& pnn = data.pnn;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the error estimator
    array_1d<double,dim> v_s_gauss;
    const array_1d<double,dim> v_gauss = prod(trans(v), N);
    const array_1d<double,dim> f_gauss = prod(trans(f), N);
    const array_1d<double,dim> grad_p = prod(trans(DN), p);
    const array_1d<double,dim> vconv_gauss = prod(trans(vconv), N);

    const double vconv_norm = norm_2(vconv_gauss);

    // Stabilization parameters
    const double c1 = 4.0;
    const double c2 = 2.0;
    const double tau1 = 1.0/((rho*dyn_tau_coeff)/delta_t + (c2*rho*vconv_norm)/h + (c1*mu)/(h*h));
    // const double tau2 = (h*h)/(c1*tau1);

    // Gauss point velocity subscale value computation
    const double cv_s_gauss0 =             N[0]*(v(0,0) - vmesh(0,0)) + N[1]*(v(1,0) - vmesh(1,0)) + N[2]*(v(2,0) - vmesh(2,0));
const double cv_s_gauss1 =             N[0]*(v(0,1) - vmesh(0,1)) + N[1]*(v(1,1) - vmesh(1,1)) + N[2]*(v(2,1) - vmesh(2,1));
            v_s_gauss[0]=-tau1*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + rho*(-N[0]*f(0,0) + N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) - N[1]*f(1,0) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) - N[2]*f(2,0) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + cv_s_gauss0*(DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0)) + cv_s_gauss1*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0))));
            v_s_gauss[1]=-tau1*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + rho*(-N[0]*f(0,1) + N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) - N[1]*f(1,1) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) - N[2]*f(2,1) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + cv_s_gauss0*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + cv_s_gauss1*(DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1))));


    const double v_gauss_norm = norm_2(v_gauss);
    const double v_s_gauss_norm = norm_2(v_s_gauss);

    return v_s_gauss_norm/v_gauss_norm;
}

}
