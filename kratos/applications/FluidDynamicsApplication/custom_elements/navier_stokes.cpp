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
    // const array_1d<double,nnodes>& pn = data.pn;
    // const array_1d<double,nnodes>& pnn = data.pnn;
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
    const double clhs11 =             v(0,0) - vmesh(0,0);
    const double clhs12 =             v(1,0) - vmesh(1,0);
    const double clhs13 =             v(2,0) - vmesh(2,0);
    const double clhs14 =             v(3,0) - vmesh(3,0);
    const double clhs15 =             N[0]*clhs11 + N[1]*clhs12 + N[2]*clhs13 + N[3]*clhs14;
    const double clhs16 =             v(0,1) - vmesh(0,1);
    const double clhs17 =             v(1,1) - vmesh(1,1);
    const double clhs18 =             v(2,1) - vmesh(2,1);
    const double clhs19 =             v(3,1) - vmesh(3,1);
    const double clhs20 =             N[0]*clhs16 + N[1]*clhs17 + N[2]*clhs18 + N[3]*clhs19;
    const double clhs21 =             v(0,2) - vmesh(0,2);
    const double clhs22 =             v(1,2) - vmesh(1,2);
    const double clhs23 =             v(2,2) - vmesh(2,2);
    const double clhs24 =             v(3,2) - vmesh(3,2);
    const double clhs25 =             N[0]*clhs21 + N[1]*clhs22 + N[2]*clhs23 + N[3]*clhs24;
    const double clhs26 =             DN(0,0)*clhs15 + DN(0,1)*clhs20 + DN(0,2)*clhs25;
    const double clhs27 =             N[0]*clhs10 + clhs26;
    const double clhs28 =             2*DN(0,0)*N[0]*rho*tau1;
    const double clhs29 =             DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0);
    const double clhs30 =             DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0);
    const double clhs31 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)) + clhs10*clhs15 + clhs20*clhs29 + clhs25*clhs30);
    const double clhs32 =             N[0]*bdf0;
    const double clhs33 =             clhs27 + clhs32;
    const double clhs34 =             pow(rho, 2);
    const double clhs35 =             DN(0,0)*clhs11 + DN(0,1)*clhs16 + DN(0,2)*clhs21 + DN(1,0)*clhs12 + DN(1,1)*clhs17 + DN(1,2)*clhs22 + DN(2,0)*clhs13 + DN(2,1)*clhs18 + DN(2,2)*clhs23 + DN(3,0)*clhs14 + DN(3,1)*clhs19 + DN(3,2)*clhs24;
    const double clhs36 =             N[0]*clhs35 + clhs26;
    const double clhs37 =             clhs34*clhs36*tau1;
    const double clhs38 =             DN(0,0)*tau2;
    const double clhs39 =             DN(0,1)*clhs38;
    const double clhs40 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs5;
    const double clhs41 =             C(1,3)*DN(0,1);
    const double clhs42 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs41;
    const double clhs43 =             C(3,5)*DN(0,0);
    const double clhs44 =             C(4,5)*DN(0,2);
    const double clhs45 =             C(1,5)*DN(0,1) + clhs43 + clhs44;
    const double clhs46 =             clhs1*rho;
    const double clhs47 =             N[0]*clhs34*clhs36*tau1;
    const double clhs48 =             2*DN(0,1)*N[0]*rho*tau1;
    const double clhs49 =             DN(0,2)*clhs38;
    const double clhs50 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs7;
    const double clhs51 =             C(3,4)*DN(0,1);
    const double clhs52 =             C(2,3)*DN(0,2) + clhs43 + clhs51;
    const double clhs53 =             C(2,5)*DN(0,2);
    const double clhs54 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs53;
    const double clhs55 =             2*DN(0,2)*N[0]*rho*tau1;
    const double clhs56 =             pow(c, -2);
    const double clhs57 =             1.0/rho;
    const double clhs58 =             N[0]*bdf0*clhs56*clhs57;
    const double clhs59 =             rho*tau1;
    const double clhs60 =             -N[0] + clhs36*clhs59 + clhs58*tau2;
    const double clhs61 =             N[0]*bdf0*rho;
    const double clhs62 =             N[1]*clhs61;
    const double clhs63 =             DN(0,0)*N[1];
    const double clhs64 =             DN(1,0)*N[0];
    const double clhs65 =             rho*tau1*(clhs63 + clhs64);
    const double clhs66 =             DN(1,0)*clhs38 + clhs31*clhs65 + clhs62;
    const double clhs67 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
    const double clhs68 =             C(0,3)*DN(1,0);
    const double clhs69 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs68;
    const double clhs70 =             C(0,5)*DN(1,0);
    const double clhs71 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs70;
    const double clhs72 =             DN(1,0)*clhs15 + DN(1,1)*clhs20 + DN(1,2)*clhs25;
    const double clhs73 =             N[1]*clhs10 + clhs72;
    const double clhs74 =             N[1]*bdf0;
    const double clhs75 =             clhs73 + clhs74;
    const double clhs76 =             N[0]*N[1]*rho;
    const double clhs77 =             DN(0,1)*N[1];
    const double clhs78 =             DN(1,1)*N[0];
    const double clhs79 =             rho*tau1*(clhs77 + clhs78);
    const double clhs80 =             clhs29*clhs76 + clhs31*clhs79;
    const double clhs81 =             DN(1,1)*clhs38;
    const double clhs82 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs68;
    const double clhs83 =             C(1,3)*DN(1,1);
    const double clhs84 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs83;
    const double clhs85 =             C(3,5)*DN(1,0);
    const double clhs86 =             C(4,5)*DN(1,2);
    const double clhs87 =             C(1,5)*DN(1,1) + clhs85 + clhs86;
    const double clhs88 =             N[1]*clhs34*clhs36*tau1;
    const double clhs89 =             DN(0,2)*N[1];
    const double clhs90 =             DN(1,2)*N[0];
    const double clhs91 =             rho*tau1*(clhs89 + clhs90);
    const double clhs92 =             clhs30*clhs76 + clhs31*clhs91;
    const double clhs93 =             DN(1,2)*clhs38;
    const double clhs94 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs70;
    const double clhs95 =             C(3,4)*DN(1,1);
    const double clhs96 =             C(2,3)*DN(1,2) + clhs85 + clhs95;
    const double clhs97 =             C(2,5)*DN(1,2);
    const double clhs98 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs97;
    const double clhs99 =             bdf0*clhs56*clhs57*tau2;
    const double clhs100 =             DN(1,0)*rho*tau1;
    const double clhs101 =             N[2]*clhs61;
    const double clhs102 =             DN(0,0)*N[2];
    const double clhs103 =             DN(2,0)*N[0];
    const double clhs104 =             rho*tau1*(clhs102 + clhs103);
    const double clhs105 =             DN(2,0)*clhs38 + clhs101 + clhs104*clhs31;
    const double clhs106 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
    const double clhs107 =             C(0,3)*DN(2,0);
    const double clhs108 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs107;
    const double clhs109 =             C(0,5)*DN(2,0);
    const double clhs110 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs109;
    const double clhs111 =             DN(2,0)*clhs15 + DN(2,1)*clhs20 + DN(2,2)*clhs25;
    const double clhs112 =             N[2]*clhs10 + clhs111;
    const double clhs113 =             N[2]*bdf0;
    const double clhs114 =             clhs112 + clhs113;
    const double clhs115 =             N[0]*N[2]*rho;
    const double clhs116 =             DN(0,1)*N[2];
    const double clhs117 =             DN(2,1)*N[0];
    const double clhs118 =             rho*tau1*(clhs116 + clhs117);
    const double clhs119 =             clhs115*clhs29 + clhs118*clhs31;
    const double clhs120 =             DN(2,1)*clhs38;
    const double clhs121 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs107;
    const double clhs122 =             C(1,3)*DN(2,1);
    const double clhs123 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs122;
    const double clhs124 =             C(3,5)*DN(2,0);
    const double clhs125 =             C(4,5)*DN(2,2);
    const double clhs126 =             C(1,5)*DN(2,1) + clhs124 + clhs125;
    const double clhs127 =             N[2]*clhs34*clhs36*tau1;
    const double clhs128 =             DN(0,2)*N[2];
    const double clhs129 =             DN(2,2)*N[0];
    const double clhs130 =             rho*tau1*(clhs128 + clhs129);
    const double clhs131 =             clhs115*clhs30 + clhs130*clhs31;
    const double clhs132 =             DN(2,2)*clhs38;
    const double clhs133 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs109;
    const double clhs134 =             C(3,4)*DN(2,1);
    const double clhs135 =             C(2,3)*DN(2,2) + clhs124 + clhs134;
    const double clhs136 =             C(2,5)*DN(2,2);
    const double clhs137 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs136;
    const double clhs138 =             DN(2,0)*rho*tau1;
    const double clhs139 =             N[3]*clhs61;
    const double clhs140 =             DN(0,0)*N[3];
    const double clhs141 =             DN(3,0)*N[0];
    const double clhs142 =             rho*tau1*(clhs140 + clhs141);
    const double clhs143 =             DN(3,0)*clhs38 + clhs139 + clhs142*clhs31;
    const double clhs144 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
    const double clhs145 =             C(0,3)*DN(3,0);
    const double clhs146 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs145;
    const double clhs147 =             C(0,5)*DN(3,0);
    const double clhs148 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs147;
    const double clhs149 =             DN(3,0)*clhs15 + DN(3,1)*clhs20 + DN(3,2)*clhs25;
    const double clhs150 =             N[3]*clhs10 + clhs149;
    const double clhs151 =             N[3]*bdf0;
    const double clhs152 =             clhs150 + clhs151;
    const double clhs153 =             N[0]*N[3]*rho;
    const double clhs154 =             DN(0,1)*N[3];
    const double clhs155 =             DN(3,1)*N[0];
    const double clhs156 =             rho*tau1*(clhs154 + clhs155);
    const double clhs157 =             clhs153*clhs29 + clhs156*clhs31;
    const double clhs158 =             DN(3,1)*clhs38;
    const double clhs159 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs145;
    const double clhs160 =             C(1,3)*DN(3,1);
    const double clhs161 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs160;
    const double clhs162 =             C(3,5)*DN(3,0);
    const double clhs163 =             C(4,5)*DN(3,2);
    const double clhs164 =             C(1,5)*DN(3,1) + clhs162 + clhs163;
    const double clhs165 =             N[3]*clhs34*clhs36*tau1;
    const double clhs166 =             DN(0,2)*N[3];
    const double clhs167 =             DN(3,2)*N[0];
    const double clhs168 =             rho*tau1*(clhs166 + clhs167);
    const double clhs169 =             clhs153*clhs30 + clhs168*clhs31;
    const double clhs170 =             DN(3,2)*clhs38;
    const double clhs171 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs147;
    const double clhs172 =             C(3,4)*DN(3,1);
    const double clhs173 =             C(2,3)*DN(3,2) + clhs162 + clhs172;
    const double clhs174 =             C(2,5)*DN(3,2);
    const double clhs175 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs174;
    const double clhs176 =             DN(3,0)*rho*tau1;
    const double clhs177 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs41;
    const double clhs178 =             C(0,4)*DN(0,0) + clhs44 + clhs51;
    const double clhs179 =             DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1);
    const double clhs180 =             N[0]*clhs179*clhs34*tau1;
    const double clhs181 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
    const double clhs182 =             DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1);
    const double clhs183 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)) + clhs15*clhs179 + clhs181*clhs20 + clhs182*clhs25);
    const double clhs184 =             pow(DN(0,1), 2);
    const double clhs185 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
    const double clhs186 =             C(1,4)*DN(0,1);
    const double clhs187 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs186;
    const double clhs188 =             N[0]*clhs181 + clhs26;
    const double clhs189 =             clhs188 + clhs32;
    const double clhs190 =             DN(0,1)*tau2;
    const double clhs191 =             DN(0,2)*clhs190;
    const double clhs192 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs186;
    const double clhs193 =             C(2,4)*DN(0,2);
    const double clhs194 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs193;
    const double clhs195 =             clhs179*clhs76 + clhs183*clhs65;
    const double clhs196 =             DN(1,0)*clhs190;
    const double clhs197 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs83;
    const double clhs198 =             C(0,4)*DN(1,0) + clhs86 + clhs95;
    const double clhs199 =             N[1]*clhs179*clhs34*tau1;
    const double clhs200 =             DN(1,1)*clhs190 + clhs183*clhs79 + clhs62;
    const double clhs201 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
    const double clhs202 =             C(1,4)*DN(1,1);
    const double clhs203 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs202;
    const double clhs204 =             N[1]*clhs181 + clhs72;
    const double clhs205 =             clhs204 + clhs74;
    const double clhs206 =             clhs182*clhs76 + clhs183*clhs91;
    const double clhs207 =             DN(1,2)*clhs190;
    const double clhs208 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs202;
    const double clhs209 =             C(2,4)*DN(1,2);
    const double clhs210 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs209;
    const double clhs211 =             DN(1,1)*rho*tau1;
    const double clhs212 =             clhs104*clhs183 + clhs115*clhs179;
    const double clhs213 =             DN(2,0)*clhs190;
    const double clhs214 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs122;
    const double clhs215 =             C(0,4)*DN(2,0) + clhs125 + clhs134;
    const double clhs216 =             N[2]*clhs179*clhs34*tau1;
    const double clhs217 =             DN(2,1)*clhs190 + clhs101 + clhs118*clhs183;
    const double clhs218 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
    const double clhs219 =             C(1,4)*DN(2,1);
    const double clhs220 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs219;
    const double clhs221 =             N[2]*clhs181 + clhs111;
    const double clhs222 =             clhs113 + clhs221;
    const double clhs223 =             clhs115*clhs182 + clhs130*clhs183;
    const double clhs224 =             DN(2,2)*clhs190;
    const double clhs225 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs219;
    const double clhs226 =             C(2,4)*DN(2,2);
    const double clhs227 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs226;
    const double clhs228 =             DN(2,1)*rho*tau1;
    const double clhs229 =             clhs142*clhs183 + clhs153*clhs179;
    const double clhs230 =             DN(3,0)*clhs190;
    const double clhs231 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs160;
    const double clhs232 =             C(0,4)*DN(3,0) + clhs163 + clhs172;
    const double clhs233 =             N[3]*clhs179*clhs34*tau1;
    const double clhs234 =             DN(3,1)*clhs190 + clhs139 + clhs156*clhs183;
    const double clhs235 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
    const double clhs236 =             C(1,4)*DN(3,1);
    const double clhs237 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs236;
    const double clhs238 =             N[3]*clhs181 + clhs149;
    const double clhs239 =             clhs151 + clhs238;
    const double clhs240 =             clhs153*clhs182 + clhs168*clhs183;
    const double clhs241 =             DN(3,2)*clhs190;
    const double clhs242 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs236;
    const double clhs243 =             C(2,4)*DN(3,2);
    const double clhs244 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs243;
    const double clhs245 =             DN(3,1)*rho*tau1;
    const double clhs246 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs53;
    const double clhs247 =             DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2);
    const double clhs248 =             N[0]*clhs247*clhs34*tau1;
    const double clhs249 =             DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2);
    const double clhs250 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
    const double clhs251 =             DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)) + clhs15*clhs247 + clhs20*clhs249 + clhs25*clhs250);
    const double clhs252 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs193;
    const double clhs253 =             pow(DN(0,2), 2);
    const double clhs254 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
    const double clhs255 =             N[0]*clhs250 + clhs26;
    const double clhs256 =             clhs255 + clhs32;
    const double clhs257 =             clhs247*clhs76 + clhs251*clhs65;
    const double clhs258 =             DN(0,2)*tau2;
    const double clhs259 =             DN(1,0)*clhs258;
    const double clhs260 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs97;
    const double clhs261 =             N[1]*clhs247*clhs34*tau1;
    const double clhs262 =             clhs249*clhs76 + clhs251*clhs79;
    const double clhs263 =             DN(1,1)*clhs258;
    const double clhs264 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs209;
    const double clhs265 =             DN(1,2)*clhs258 + clhs251*clhs91 + clhs62;
    const double clhs266 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
    const double clhs267 =             N[1]*clhs250 + clhs72;
    const double clhs268 =             clhs267 + clhs74;
    const double clhs269 =             DN(1,2)*rho*tau1;
    const double clhs270 =             clhs104*clhs251 + clhs115*clhs247;
    const double clhs271 =             DN(2,0)*clhs258;
    const double clhs272 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs136;
    const double clhs273 =             N[2]*clhs247*clhs34*tau1;
    const double clhs274 =             clhs115*clhs249 + clhs118*clhs251;
    const double clhs275 =             DN(2,1)*clhs258;
    const double clhs276 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs226;
    const double clhs277 =             DN(2,2)*clhs258 + clhs101 + clhs130*clhs251;
    const double clhs278 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
    const double clhs279 =             N[2]*clhs250 + clhs111;
    const double clhs280 =             clhs113 + clhs279;
    const double clhs281 =             DN(2,2)*rho*tau1;
    const double clhs282 =             clhs142*clhs251 + clhs153*clhs247;
    const double clhs283 =             DN(3,0)*clhs258;
    const double clhs284 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs174;
    const double clhs285 =             N[3]*clhs247*clhs34*tau1;
    const double clhs286 =             clhs153*clhs249 + clhs156*clhs251;
    const double clhs287 =             DN(3,1)*clhs258;
    const double clhs288 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs243;
    const double clhs289 =             DN(3,2)*clhs258 + clhs139 + clhs168*clhs251;
    const double clhs290 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
    const double clhs291 =             N[3]*clhs250 + clhs149;
    const double clhs292 =             clhs151 + clhs291;
    const double clhs293 =             DN(3,2)*rho*tau1;
    const double clhs294 =             DN(0,1)*N[0]*rho*tau1;
    const double clhs295 =             DN(0,2)*N[0]*rho*tau1;
    const double clhs296 =             DN(0,0)*rho*tau1;
    const double clhs297 =             DN(0,0)*N[0]*rho*tau1;
    const double clhs298 =             DN(0,1)*rho*tau1;
    const double clhs299 =             DN(0,2)*rho*tau1;
    const double clhs300 =             bdf0*clhs56*clhs57;
    const double clhs301 =             DN(0,1)*N[1]*rho*tau1;
    const double clhs302 =             DN(0,2)*N[1]*rho*tau1;
    const double clhs303 =             DN(0,0)*N[1]*rho*tau1;
    const double clhs304 =             DN(0,0)*tau1;
    const double clhs305 =             DN(0,1)*tau1;
    const double clhs306 =             DN(0,2)*tau1;
    const double clhs307 =             DN(1,0)*clhs304 + DN(1,1)*clhs305 + DN(1,2)*clhs306 + N[1]*clhs58;
    const double clhs308 =             DN(0,1)*N[2]*rho*tau1;
    const double clhs309 =             DN(0,2)*N[2]*rho*tau1;
    const double clhs310 =             DN(0,0)*N[2]*rho*tau1;
    const double clhs311 =             DN(2,0)*clhs304 + DN(2,1)*clhs305 + DN(2,2)*clhs306 + N[2]*clhs58;
    const double clhs312 =             DN(0,1)*N[3]*rho*tau1;
    const double clhs313 =             DN(0,2)*N[3]*rho*tau1;
    const double clhs314 =             DN(0,0)*N[3]*rho*tau1;
    const double clhs315 =             DN(3,0)*clhs304 + DN(3,1)*clhs305 + DN(3,2)*clhs306 + N[3]*clhs58;
    const double clhs316 =             N[1]*rho;
    const double clhs317 =             N[1]*clhs35 + clhs72;
    const double clhs318 =             clhs317*clhs34*tau1;
    const double clhs319 =             N[0]*clhs29*clhs34*tau1;
    const double clhs320 =             N[0]*clhs30*clhs34*tau1;
    const double clhs321 =             pow(DN(1,0), 2);
    const double clhs322 =             pow(N[1], 2);
    const double clhs323 =             clhs2*clhs322;
    const double clhs324 =             2*DN(1,0)*N[1]*rho*tau1;
    const double clhs325 =             DN(1,0)*tau2;
    const double clhs326 =             DN(1,1)*clhs325;
    const double clhs327 =             clhs322*rho;
    const double clhs328 =             N[1]*clhs29*clhs34*tau1;
    const double clhs329 =             2*DN(1,1)*N[1]*rho*tau1;
    const double clhs330 =             DN(1,2)*clhs325;
    const double clhs331 =             N[1]*clhs30*clhs34*tau1;
    const double clhs332 =             2*DN(1,2)*N[1]*rho*tau1;
    const double clhs333 =             N[1]*bdf0*clhs56*clhs57;
    const double clhs334 =             -N[1] + clhs317*clhs59 + clhs333*tau2;
    const double clhs335 =             N[1]*bdf0*rho;
    const double clhs336 =             N[2]*clhs335;
    const double clhs337 =             DN(1,0)*N[2];
    const double clhs338 =             DN(2,0)*N[1];
    const double clhs339 =             rho*tau1*(clhs337 + clhs338);
    const double clhs340 =             DN(2,0)*clhs325 + clhs31*clhs339 + clhs336;
    const double clhs341 =             N[1]*N[2]*rho;
    const double clhs342 =             DN(1,1)*N[2];
    const double clhs343 =             DN(2,1)*N[1];
    const double clhs344 =             rho*tau1*(clhs342 + clhs343);
    const double clhs345 =             clhs29*clhs341 + clhs31*clhs344;
    const double clhs346 =             DN(2,1)*clhs325;
    const double clhs347 =             N[2]*clhs29*clhs34*tau1;
    const double clhs348 =             DN(1,2)*N[2];
    const double clhs349 =             DN(2,2)*N[1];
    const double clhs350 =             rho*tau1*(clhs348 + clhs349);
    const double clhs351 =             clhs30*clhs341 + clhs31*clhs350;
    const double clhs352 =             DN(2,2)*clhs325;
    const double clhs353 =             N[2]*clhs30*clhs34*tau1;
    const double clhs354 =             N[3]*clhs335;
    const double clhs355 =             DN(1,0)*N[3];
    const double clhs356 =             DN(3,0)*N[1];
    const double clhs357 =             rho*tau1*(clhs355 + clhs356);
    const double clhs358 =             DN(3,0)*clhs325 + clhs31*clhs357 + clhs354;
    const double clhs359 =             N[1]*N[3]*rho;
    const double clhs360 =             DN(1,1)*N[3];
    const double clhs361 =             DN(3,1)*N[1];
    const double clhs362 =             rho*tau1*(clhs360 + clhs361);
    const double clhs363 =             clhs29*clhs359 + clhs31*clhs362;
    const double clhs364 =             DN(3,1)*clhs325;
    const double clhs365 =             N[3]*clhs29*clhs34*tau1;
    const double clhs366 =             DN(1,2)*N[3];
    const double clhs367 =             DN(3,2)*N[1];
    const double clhs368 =             rho*tau1*(clhs366 + clhs367);
    const double clhs369 =             clhs30*clhs359 + clhs31*clhs368;
    const double clhs370 =             DN(3,2)*clhs325;
    const double clhs371 =             N[3]*clhs30*clhs34*tau1;
    const double clhs372 =             N[0]*clhs182*clhs34*tau1;
    const double clhs373 =             pow(DN(1,1), 2);
    const double clhs374 =             DN(1,1)*tau2;
    const double clhs375 =             DN(1,2)*clhs374;
    const double clhs376 =             N[1]*clhs182*clhs34*tau1;
    const double clhs377 =             clhs179*clhs341 + clhs183*clhs339;
    const double clhs378 =             DN(2,0)*clhs374;
    const double clhs379 =             DN(2,1)*clhs374 + clhs183*clhs344 + clhs336;
    const double clhs380 =             clhs182*clhs341 + clhs183*clhs350;
    const double clhs381 =             DN(2,2)*clhs374;
    const double clhs382 =             N[2]*clhs182*clhs34*tau1;
    const double clhs383 =             clhs179*clhs359 + clhs183*clhs357;
    const double clhs384 =             DN(3,0)*clhs374;
    const double clhs385 =             DN(3,1)*clhs374 + clhs183*clhs362 + clhs354;
    const double clhs386 =             clhs182*clhs359 + clhs183*clhs368;
    const double clhs387 =             DN(3,2)*clhs374;
    const double clhs388 =             N[3]*clhs182*clhs34*tau1;
    const double clhs389 =             N[0]*clhs249*clhs34*tau1;
    const double clhs390 =             N[1]*clhs249*clhs34*tau1;
    const double clhs391 =             pow(DN(1,2), 2);
    const double clhs392 =             clhs247*clhs341 + clhs251*clhs339;
    const double clhs393 =             DN(1,2)*tau2;
    const double clhs394 =             DN(2,0)*clhs393;
    const double clhs395 =             clhs249*clhs341 + clhs251*clhs344;
    const double clhs396 =             DN(2,1)*clhs393;
    const double clhs397 =             N[2]*clhs249*clhs34*tau1;
    const double clhs398 =             DN(2,2)*clhs393 + clhs251*clhs350 + clhs336;
    const double clhs399 =             clhs247*clhs359 + clhs251*clhs357;
    const double clhs400 =             DN(3,0)*clhs393;
    const double clhs401 =             clhs249*clhs359 + clhs251*clhs362;
    const double clhs402 =             DN(3,1)*clhs393;
    const double clhs403 =             N[3]*clhs249*clhs34*tau1;
    const double clhs404 =             DN(3,2)*clhs393 + clhs251*clhs368 + clhs354;
    const double clhs405 =             DN(1,1)*N[0]*rho*tau1;
    const double clhs406 =             DN(1,2)*N[0]*rho*tau1;
    const double clhs407 =             DN(1,0)*N[0]*rho*tau1;
    const double clhs408 =             DN(1,1)*N[1]*rho*tau1;
    const double clhs409 =             DN(1,2)*N[1]*rho*tau1;
    const double clhs410 =             DN(1,0)*N[1]*rho*tau1;
    const double clhs411 =             DN(1,1)*N[2]*rho*tau1;
    const double clhs412 =             DN(1,2)*N[2]*rho*tau1;
    const double clhs413 =             DN(1,0)*N[2]*rho*tau1;
    const double clhs414 =             DN(1,0)*tau1;
    const double clhs415 =             DN(1,1)*tau1;
    const double clhs416 =             DN(1,2)*tau1;
    const double clhs417 =             DN(2,0)*clhs414 + DN(2,1)*clhs415 + DN(2,2)*clhs416 + N[2]*clhs333;
    const double clhs418 =             DN(1,1)*N[3]*rho*tau1;
    const double clhs419 =             DN(1,2)*N[3]*rho*tau1;
    const double clhs420 =             DN(1,0)*N[3]*rho*tau1;
    const double clhs421 =             DN(3,0)*clhs414 + DN(3,1)*clhs415 + DN(3,2)*clhs416 + N[3]*clhs333;
    const double clhs422 =             N[2]*rho;
    const double clhs423 =             N[2]*clhs35 + clhs111;
    const double clhs424 =             clhs34*clhs423*tau1;
    const double clhs425 =             pow(DN(2,0), 2);
    const double clhs426 =             pow(N[2], 2);
    const double clhs427 =             clhs2*clhs426;
    const double clhs428 =             2*DN(2,0)*N[2]*rho*tau1;
    const double clhs429 =             DN(2,0)*tau2;
    const double clhs430 =             DN(2,1)*clhs429;
    const double clhs431 =             clhs426*rho;
    const double clhs432 =             2*DN(2,1)*N[2]*rho*tau1;
    const double clhs433 =             DN(2,2)*clhs429;
    const double clhs434 =             2*DN(2,2)*N[2]*rho*tau1;
    const double clhs435 =             clhs56*clhs57*tau2;
    const double clhs436 =             -N[2] + clhs113*clhs435 + clhs423*clhs59;
    const double clhs437 =             N[2]*N[3]*bdf0;
    const double clhs438 =             clhs437*rho;
    const double clhs439 =             DN(2,0)*N[3];
    const double clhs440 =             DN(3,0)*N[2];
    const double clhs441 =             rho*tau1*(clhs439 + clhs440);
    const double clhs442 =             DN(3,0)*clhs429 + clhs31*clhs441 + clhs438;
    const double clhs443 =             N[2]*N[3]*rho;
    const double clhs444 =             DN(2,1)*N[3];
    const double clhs445 =             DN(3,1)*N[2];
    const double clhs446 =             rho*tau1*(clhs444 + clhs445);
    const double clhs447 =             clhs29*clhs443 + clhs31*clhs446;
    const double clhs448 =             DN(3,1)*clhs429;
    const double clhs449 =             DN(2,2)*N[3];
    const double clhs450 =             DN(3,2)*N[2];
    const double clhs451 =             rho*tau1*(clhs449 + clhs450);
    const double clhs452 =             clhs30*clhs443 + clhs31*clhs451;
    const double clhs453 =             DN(3,2)*clhs429;
    const double clhs454 =             pow(DN(2,1), 2);
    const double clhs455 =             DN(2,1)*tau2;
    const double clhs456 =             DN(2,2)*clhs455;
    const double clhs457 =             clhs179*clhs443 + clhs183*clhs441;
    const double clhs458 =             DN(3,0)*clhs455;
    const double clhs459 =             DN(3,1)*clhs455 + clhs183*clhs446 + clhs438;
    const double clhs460 =             clhs182*clhs443 + clhs183*clhs451;
    const double clhs461 =             DN(3,2)*clhs455;
    const double clhs462 =             pow(DN(2,2), 2);
    const double clhs463 =             clhs247*clhs443 + clhs251*clhs441;
    const double clhs464 =             DN(2,2)*tau2;
    const double clhs465 =             DN(3,0)*clhs464;
    const double clhs466 =             clhs249*clhs443 + clhs251*clhs446;
    const double clhs467 =             DN(3,1)*clhs464;
    const double clhs468 =             DN(3,2)*clhs464 + clhs251*clhs451 + clhs438;
    const double clhs469 =             DN(2,1)*N[0]*rho*tau1;
    const double clhs470 =             DN(2,2)*N[0]*rho*tau1;
    const double clhs471 =             DN(2,0)*N[0]*rho*tau1;
    const double clhs472 =             DN(2,1)*N[1]*rho*tau1;
    const double clhs473 =             DN(2,2)*N[1]*rho*tau1;
    const double clhs474 =             DN(2,0)*N[1]*rho*tau1;
    const double clhs475 =             DN(2,1)*N[2]*rho*tau1;
    const double clhs476 =             DN(2,2)*N[2]*rho*tau1;
    const double clhs477 =             DN(2,0)*N[2]*rho*tau1;
    const double clhs478 =             DN(2,1)*N[3]*rho*tau1;
    const double clhs479 =             DN(2,2)*N[3]*rho*tau1;
    const double clhs480 =             DN(2,0)*N[3]*rho*tau1;
    const double clhs481 =             DN(3,0)*(DN(2,0)*tau1) + DN(3,1)*(DN(2,1)*tau1) + DN(3,2)*(DN(2,2)*tau1) + clhs437*clhs56*clhs57;
    const double clhs482 =             N[3]*rho;
    const double clhs483 =             N[3]*clhs35 + clhs149;
    const double clhs484 =             clhs34*clhs483*tau1;
    const double clhs485 =             pow(DN(3,0), 2);
    const double clhs486 =             pow(N[3], 2);
    const double clhs487 =             clhs2*clhs486;
    const double clhs488 =             2*DN(3,0)*N[3]*rho*tau1;
    const double clhs489 =             DN(3,0)*tau2;
    const double clhs490 =             DN(3,1)*clhs489;
    const double clhs491 =             clhs486*rho;
    const double clhs492 =             2*DN(3,1)*N[3]*rho*tau1;
    const double clhs493 =             DN(3,2)*clhs489;
    const double clhs494 =             2*DN(3,2)*N[3]*rho*tau1;
    const double clhs495 =             -N[3] + clhs151*clhs435 + clhs483*clhs59;
    const double clhs496 =             pow(DN(3,1), 2);
    const double clhs497 =             DN(3,1)*DN(3,2)*tau2;
    const double clhs498 =             pow(DN(3,2), 2);
    const double clhs499 =             DN(3,1)*N[0]*rho*tau1;
    const double clhs500 =             DN(3,2)*N[0]*rho*tau1;
    const double clhs501 =             DN(3,0)*N[0]*rho*tau1;
    const double clhs502 =             DN(3,1)*N[1]*rho*tau1;
    const double clhs503 =             DN(3,2)*N[1]*rho*tau1;
    const double clhs504 =             DN(3,0)*N[1]*rho*tau1;
    const double clhs505 =             DN(3,1)*N[2]*rho*tau1;
    const double clhs506 =             DN(3,2)*N[2]*rho*tau1;
    const double clhs507 =             DN(3,0)*N[2]*rho*tau1;
    const double clhs508 =             DN(3,1)*N[3]*rho*tau1;
    const double clhs509 =             DN(3,2)*N[3]*rho*tau1;
    const double clhs510 =             DN(3,0)*N[3]*rho*tau1;

    lhs(0,0)=DN(0,0)*clhs4 + DN(0,1)*clhs6 + DN(0,2)*clhs8 + clhs0*tau2 + clhs27*clhs9 + clhs28*clhs31 + clhs3 + clhs33*clhs37;
    lhs(0,1)=DN(0,0)*clhs40 + DN(0,1)*clhs42 + DN(0,2)*clhs45 + clhs29*clhs46 + clhs29*clhs47 + clhs31*clhs48 + clhs39;
    lhs(0,2)=DN(0,0)*clhs50 + DN(0,1)*clhs52 + DN(0,2)*clhs54 + clhs30*clhs46 + clhs30*clhs47 + clhs31*clhs55 + clhs49;
    lhs(0,3)=DN(0,0)*clhs60;
    lhs(0,4)=DN(0,0)*clhs67 + DN(0,1)*clhs69 + DN(0,2)*clhs71 + clhs37*clhs75 + clhs66 + clhs73*clhs9;
    lhs(0,5)=DN(0,0)*clhs82 + DN(0,1)*clhs84 + DN(0,2)*clhs87 + clhs29*clhs88 + clhs80 + clhs81;
    lhs(0,6)=DN(0,0)*clhs94 + DN(0,1)*clhs96 + DN(0,2)*clhs98 + clhs30*clhs88 + clhs92 + clhs93;
    lhs(0,7)=clhs100*clhs36 + clhs63*clhs99 - clhs63;
    lhs(0,8)=DN(0,0)*clhs106 + DN(0,1)*clhs108 + DN(0,2)*clhs110 + clhs105 + clhs112*clhs9 + clhs114*clhs37;
    lhs(0,9)=DN(0,0)*clhs121 + DN(0,1)*clhs123 + DN(0,2)*clhs126 + clhs119 + clhs120 + clhs127*clhs29;
    lhs(0,10)=DN(0,0)*clhs133 + DN(0,1)*clhs135 + DN(0,2)*clhs137 + clhs127*clhs30 + clhs131 + clhs132;
    lhs(0,11)=clhs102*clhs99 - clhs102 + clhs138*clhs36;
    lhs(0,12)=DN(0,0)*clhs144 + DN(0,1)*clhs146 + DN(0,2)*clhs148 + clhs143 + clhs150*clhs9 + clhs152*clhs37;
    lhs(0,13)=DN(0,0)*clhs159 + DN(0,1)*clhs161 + DN(0,2)*clhs164 + clhs157 + clhs158 + clhs165*clhs29;
    lhs(0,14)=DN(0,0)*clhs171 + DN(0,1)*clhs173 + DN(0,2)*clhs175 + clhs165*clhs30 + clhs169 + clhs170;
    lhs(0,15)=clhs140*clhs99 - clhs140 + clhs176*clhs36;
    lhs(1,0)=DN(0,0)*clhs6 + DN(0,1)*clhs177 + DN(0,2)*clhs178 + clhs179*clhs46 + clhs180*clhs36 + clhs183*clhs28 + clhs39;
    lhs(1,1)=DN(0,0)*clhs42 + DN(0,1)*clhs185 + DN(0,2)*clhs187 + clhs183*clhs48 + clhs184*tau2 + clhs188*clhs9 + clhs189*clhs37 + clhs3;
    lhs(1,2)=DN(0,0)*clhs52 + DN(0,1)*clhs192 + DN(0,2)*clhs194 + clhs182*clhs46 + clhs182*clhs47 + clhs183*clhs55 + clhs191;
    lhs(1,3)=DN(0,1)*clhs60;
    lhs(1,4)=DN(0,0)*clhs69 + DN(0,1)*clhs197 + DN(0,2)*clhs198 + clhs195 + clhs196 + clhs199*clhs36;
    lhs(1,5)=DN(0,0)*clhs84 + DN(0,1)*clhs201 + DN(0,2)*clhs203 + clhs200 + clhs204*clhs9 + clhs205*clhs37;
    lhs(1,6)=DN(0,0)*clhs96 + DN(0,1)*clhs208 + DN(0,2)*clhs210 + clhs182*clhs88 + clhs206 + clhs207;
    lhs(1,7)=clhs211*clhs36 + clhs77*clhs99 - clhs77;
    lhs(1,8)=DN(0,0)*clhs108 + DN(0,1)*clhs214 + DN(0,2)*clhs215 + clhs212 + clhs213 + clhs216*clhs36;
    lhs(1,9)=DN(0,0)*clhs123 + DN(0,1)*clhs218 + DN(0,2)*clhs220 + clhs217 + clhs221*clhs9 + clhs222*clhs37;
    lhs(1,10)=DN(0,0)*clhs135 + DN(0,1)*clhs225 + DN(0,2)*clhs227 + clhs127*clhs182 + clhs223 + clhs224;
    lhs(1,11)=clhs116*clhs99 - clhs116 + clhs228*clhs36;
    lhs(1,12)=DN(0,0)*clhs146 + DN(0,1)*clhs231 + DN(0,2)*clhs232 + clhs229 + clhs230 + clhs233*clhs36;
    lhs(1,13)=DN(0,0)*clhs161 + DN(0,1)*clhs235 + DN(0,2)*clhs237 + clhs234 + clhs238*clhs9 + clhs239*clhs37;
    lhs(1,14)=DN(0,0)*clhs173 + DN(0,1)*clhs242 + DN(0,2)*clhs244 + clhs165*clhs182 + clhs240 + clhs241;
    lhs(1,15)=clhs154*clhs99 - clhs154 + clhs245*clhs36;
    lhs(2,0)=DN(0,0)*clhs8 + DN(0,1)*clhs178 + DN(0,2)*clhs246 + clhs247*clhs46 + clhs248*clhs36 + clhs251*clhs28 + clhs49;
    lhs(2,1)=DN(0,0)*clhs45 + DN(0,1)*clhs187 + DN(0,2)*clhs252 + clhs191 + clhs249*clhs46 + clhs249*clhs47 + clhs251*clhs48;
    lhs(2,2)=DN(0,0)*clhs54 + DN(0,1)*clhs194 + DN(0,2)*clhs254 + clhs251*clhs55 + clhs253*tau2 + clhs255*clhs9 + clhs256*clhs37 + clhs3;
    lhs(2,3)=DN(0,2)*clhs60;
    lhs(2,4)=DN(0,0)*clhs71 + DN(0,1)*clhs198 + DN(0,2)*clhs260 + clhs257 + clhs259 + clhs261*clhs36;
    lhs(2,5)=DN(0,0)*clhs87 + DN(0,1)*clhs203 + DN(0,2)*clhs264 + clhs249*clhs88 + clhs262 + clhs263;
    lhs(2,6)=DN(0,0)*clhs98 + DN(0,1)*clhs210 + DN(0,2)*clhs266 + clhs265 + clhs267*clhs9 + clhs268*clhs37;
    lhs(2,7)=clhs269*clhs36 + clhs89*clhs99 - clhs89;
    lhs(2,8)=DN(0,0)*clhs110 + DN(0,1)*clhs215 + DN(0,2)*clhs272 + clhs270 + clhs271 + clhs273*clhs36;
    lhs(2,9)=DN(0,0)*clhs126 + DN(0,1)*clhs220 + DN(0,2)*clhs276 + clhs127*clhs249 + clhs274 + clhs275;
    lhs(2,10)=DN(0,0)*clhs137 + DN(0,1)*clhs227 + DN(0,2)*clhs278 + clhs277 + clhs279*clhs9 + clhs280*clhs37;
    lhs(2,11)=clhs128*clhs99 - clhs128 + clhs281*clhs36;
    lhs(2,12)=DN(0,0)*clhs148 + DN(0,1)*clhs232 + DN(0,2)*clhs284 + clhs282 + clhs283 + clhs285*clhs36;
    lhs(2,13)=DN(0,0)*clhs164 + DN(0,1)*clhs237 + DN(0,2)*clhs288 + clhs165*clhs249 + clhs286 + clhs287;
    lhs(2,14)=DN(0,0)*clhs175 + DN(0,1)*clhs244 + DN(0,2)*clhs290 + clhs289 + clhs291*clhs9 + clhs292*clhs37;
    lhs(2,15)=clhs166*clhs99 - clhs166 + clhs293*clhs36;
    lhs(3,0)=DN(0,0)*N[0] + clhs179*clhs294 + clhs247*clhs295 + clhs296*clhs33;
    lhs(3,1)=DN(0,1)*N[0] + clhs189*clhs298 + clhs249*clhs295 + clhs29*clhs297;
    lhs(3,2)=DN(0,2)*N[0] + clhs182*clhs294 + clhs256*clhs299 + clhs297*clhs30;
    lhs(3,3)=clhs0*tau1 + clhs1*clhs300 + clhs184*tau1 + clhs253*tau1;
    lhs(3,4)=clhs179*clhs301 + clhs247*clhs302 + clhs296*clhs75 + clhs64;
    lhs(3,5)=clhs205*clhs298 + clhs249*clhs302 + clhs29*clhs303 + clhs78;
    lhs(3,6)=clhs182*clhs301 + clhs268*clhs299 + clhs30*clhs303 + clhs90;
    lhs(3,7)=clhs307;
    lhs(3,8)=clhs103 + clhs114*clhs296 + clhs179*clhs308 + clhs247*clhs309;
    lhs(3,9)=clhs117 + clhs222*clhs298 + clhs249*clhs309 + clhs29*clhs310;
    lhs(3,10)=clhs129 + clhs182*clhs308 + clhs280*clhs299 + clhs30*clhs310;
    lhs(3,11)=clhs311;
    lhs(3,12)=clhs141 + clhs152*clhs296 + clhs179*clhs312 + clhs247*clhs313;
    lhs(3,13)=clhs155 + clhs239*clhs298 + clhs249*clhs313 + clhs29*clhs314;
    lhs(3,14)=clhs167 + clhs182*clhs312 + clhs292*clhs299 + clhs30*clhs314;
    lhs(3,15)=clhs315;
    lhs(4,0)=DN(1,0)*clhs4 + DN(1,1)*clhs6 + DN(1,2)*clhs8 + clhs27*clhs316 + clhs318*clhs33 + clhs66;
    lhs(4,1)=DN(1,0)*clhs40 + DN(1,1)*clhs42 + DN(1,2)*clhs45 + clhs196 + clhs317*clhs319 + clhs80;
    lhs(4,2)=DN(1,0)*clhs50 + DN(1,1)*clhs52 + DN(1,2)*clhs54 + clhs259 + clhs317*clhs320 + clhs92;
    lhs(4,3)=clhs296*clhs317 + clhs64*clhs99 - clhs64;
    lhs(4,4)=DN(1,0)*clhs67 + DN(1,1)*clhs69 + DN(1,2)*clhs71 + clhs31*clhs324 + clhs316*clhs73 + clhs318*clhs75 + clhs321*tau2 + clhs323;
    lhs(4,5)=DN(1,0)*clhs82 + DN(1,1)*clhs84 + DN(1,2)*clhs87 + clhs29*clhs327 + clhs31*clhs329 + clhs317*clhs328 + clhs326;
    lhs(4,6)=DN(1,0)*clhs94 + DN(1,1)*clhs96 + DN(1,2)*clhs98 + clhs30*clhs327 + clhs31*clhs332 + clhs317*clhs331 + clhs330;
    lhs(4,7)=DN(1,0)*clhs334;
    lhs(4,8)=DN(1,0)*clhs106 + DN(1,1)*clhs108 + DN(1,2)*clhs110 + clhs112*clhs316 + clhs114*clhs318 + clhs340;
    lhs(4,9)=DN(1,0)*clhs121 + DN(1,1)*clhs123 + DN(1,2)*clhs126 + clhs317*clhs347 + clhs345 + clhs346;
    lhs(4,10)=DN(1,0)*clhs133 + DN(1,1)*clhs135 + DN(1,2)*clhs137 + clhs317*clhs353 + clhs351 + clhs352;
    lhs(4,11)=clhs138*clhs317 + clhs337*clhs99 - clhs337;
    lhs(4,12)=DN(1,0)*clhs144 + DN(1,1)*clhs146 + DN(1,2)*clhs148 + clhs150*clhs316 + clhs152*clhs318 + clhs358;
    lhs(4,13)=DN(1,0)*clhs159 + DN(1,1)*clhs161 + DN(1,2)*clhs164 + clhs317*clhs365 + clhs363 + clhs364;
    lhs(4,14)=DN(1,0)*clhs171 + DN(1,1)*clhs173 + DN(1,2)*clhs175 + clhs317*clhs371 + clhs369 + clhs370;
    lhs(4,15)=clhs176*clhs317 + clhs355*clhs99 - clhs355;
    lhs(5,0)=DN(1,0)*clhs6 + DN(1,1)*clhs177 + DN(1,2)*clhs178 + clhs180*clhs317 + clhs195 + clhs81;
    lhs(5,1)=DN(1,0)*clhs42 + DN(1,1)*clhs185 + DN(1,2)*clhs187 + clhs188*clhs316 + clhs189*clhs318 + clhs200;
    lhs(5,2)=DN(1,0)*clhs52 + DN(1,1)*clhs192 + DN(1,2)*clhs194 + clhs206 + clhs263 + clhs317*clhs372;
    lhs(5,3)=clhs298*clhs317 + clhs78*clhs99 - clhs78;
    lhs(5,4)=DN(1,0)*clhs69 + DN(1,1)*clhs197 + DN(1,2)*clhs198 + clhs179*clhs327 + clhs183*clhs324 + clhs199*clhs317 + clhs326;
    lhs(5,5)=DN(1,0)*clhs84 + DN(1,1)*clhs201 + DN(1,2)*clhs203 + clhs183*clhs329 + clhs204*clhs316 + clhs205*clhs318 + clhs323 + clhs373*tau2;
    lhs(5,6)=DN(1,0)*clhs96 + DN(1,1)*clhs208 + DN(1,2)*clhs210 + clhs182*clhs327 + clhs183*clhs332 + clhs317*clhs376 + clhs375;
    lhs(5,7)=DN(1,1)*clhs334;
    lhs(5,8)=DN(1,0)*clhs108 + DN(1,1)*clhs214 + DN(1,2)*clhs215 + clhs216*clhs317 + clhs377 + clhs378;
    lhs(5,9)=DN(1,0)*clhs123 + DN(1,1)*clhs218 + DN(1,2)*clhs220 + clhs221*clhs316 + clhs222*clhs318 + clhs379;
    lhs(5,10)=DN(1,0)*clhs135 + DN(1,1)*clhs225 + DN(1,2)*clhs227 + clhs317*clhs382 + clhs380 + clhs381;
    lhs(5,11)=clhs228*clhs317 + clhs342*clhs99 - clhs342;
    lhs(5,12)=DN(1,0)*clhs146 + DN(1,1)*clhs231 + DN(1,2)*clhs232 + clhs233*clhs317 + clhs383 + clhs384;
    lhs(5,13)=DN(1,0)*clhs161 + DN(1,1)*clhs235 + DN(1,2)*clhs237 + clhs238*clhs316 + clhs239*clhs318 + clhs385;
    lhs(5,14)=DN(1,0)*clhs173 + DN(1,1)*clhs242 + DN(1,2)*clhs244 + clhs317*clhs388 + clhs386 + clhs387;
    lhs(5,15)=clhs245*clhs317 + clhs360*clhs99 - clhs360;
    lhs(6,0)=DN(1,0)*clhs8 + DN(1,1)*clhs178 + DN(1,2)*clhs246 + clhs248*clhs317 + clhs257 + clhs93;
    lhs(6,1)=DN(1,0)*clhs45 + DN(1,1)*clhs187 + DN(1,2)*clhs252 + clhs207 + clhs262 + clhs317*clhs389;
    lhs(6,2)=DN(1,0)*clhs54 + DN(1,1)*clhs194 + DN(1,2)*clhs254 + clhs255*clhs316 + clhs256*clhs318 + clhs265;
    lhs(6,3)=clhs299*clhs317 + clhs90*clhs99 - clhs90;
    lhs(6,4)=DN(1,0)*clhs71 + DN(1,1)*clhs198 + DN(1,2)*clhs260 + clhs247*clhs327 + clhs251*clhs324 + clhs261*clhs317 + clhs330;
    lhs(6,5)=DN(1,0)*clhs87 + DN(1,1)*clhs203 + DN(1,2)*clhs264 + clhs249*clhs327 + clhs251*clhs329 + clhs317*clhs390 + clhs375;
    lhs(6,6)=DN(1,0)*clhs98 + DN(1,1)*clhs210 + DN(1,2)*clhs266 + clhs251*clhs332 + clhs267*clhs316 + clhs268*clhs318 + clhs323 + clhs391*tau2;
    lhs(6,7)=DN(1,2)*clhs334;
    lhs(6,8)=DN(1,0)*clhs110 + DN(1,1)*clhs215 + DN(1,2)*clhs272 + clhs273*clhs317 + clhs392 + clhs394;
    lhs(6,9)=DN(1,0)*clhs126 + DN(1,1)*clhs220 + DN(1,2)*clhs276 + clhs317*clhs397 + clhs395 + clhs396;
    lhs(6,10)=DN(1,0)*clhs137 + DN(1,1)*clhs227 + DN(1,2)*clhs278 + clhs279*clhs316 + clhs280*clhs318 + clhs398;
    lhs(6,11)=clhs281*clhs317 + clhs348*clhs99 - clhs348;
    lhs(6,12)=DN(1,0)*clhs148 + DN(1,1)*clhs232 + DN(1,2)*clhs284 + clhs285*clhs317 + clhs399 + clhs400;
    lhs(6,13)=DN(1,0)*clhs164 + DN(1,1)*clhs237 + DN(1,2)*clhs288 + clhs317*clhs403 + clhs401 + clhs402;
    lhs(6,14)=DN(1,0)*clhs175 + DN(1,1)*clhs244 + DN(1,2)*clhs290 + clhs291*clhs316 + clhs292*clhs318 + clhs404;
    lhs(6,15)=clhs293*clhs317 + clhs366*clhs99 - clhs366;
    lhs(7,0)=clhs100*clhs33 + clhs179*clhs405 + clhs247*clhs406 + clhs63;
    lhs(7,1)=clhs189*clhs211 + clhs249*clhs406 + clhs29*clhs407 + clhs77;
    lhs(7,2)=clhs182*clhs405 + clhs256*clhs269 + clhs30*clhs407 + clhs89;
    lhs(7,3)=clhs307;
    lhs(7,4)=DN(1,0)*N[1] + clhs100*clhs75 + clhs179*clhs408 + clhs247*clhs409;
    lhs(7,5)=DN(1,1)*N[1] + clhs205*clhs211 + clhs249*clhs409 + clhs29*clhs410;
    lhs(7,6)=DN(1,2)*N[1] + clhs182*clhs408 + clhs268*clhs269 + clhs30*clhs410;
    lhs(7,7)=clhs300*clhs322 + clhs321*tau1 + clhs373*tau1 + clhs391*tau1;
    lhs(7,8)=clhs100*clhs114 + clhs179*clhs411 + clhs247*clhs412 + clhs338;
    lhs(7,9)=clhs211*clhs222 + clhs249*clhs412 + clhs29*clhs413 + clhs343;
    lhs(7,10)=clhs182*clhs411 + clhs269*clhs280 + clhs30*clhs413 + clhs349;
    lhs(7,11)=clhs417;
    lhs(7,12)=clhs100*clhs152 + clhs179*clhs418 + clhs247*clhs419 + clhs356;
    lhs(7,13)=clhs211*clhs239 + clhs249*clhs419 + clhs29*clhs420 + clhs361;
    lhs(7,14)=clhs182*clhs418 + clhs269*clhs292 + clhs30*clhs420 + clhs367;
    lhs(7,15)=clhs421;
    lhs(8,0)=DN(2,0)*clhs4 + DN(2,1)*clhs6 + DN(2,2)*clhs8 + clhs105 + clhs27*clhs422 + clhs33*clhs424;
    lhs(8,1)=DN(2,0)*clhs40 + DN(2,1)*clhs42 + DN(2,2)*clhs45 + clhs119 + clhs213 + clhs319*clhs423;
    lhs(8,2)=DN(2,0)*clhs50 + DN(2,1)*clhs52 + DN(2,2)*clhs54 + clhs131 + clhs271 + clhs320*clhs423;
    lhs(8,3)=clhs103*clhs99 - clhs103 + clhs296*clhs423;
    lhs(8,4)=DN(2,0)*clhs67 + DN(2,1)*clhs69 + DN(2,2)*clhs71 + clhs340 + clhs422*clhs73 + clhs424*clhs75;
    lhs(8,5)=DN(2,0)*clhs82 + DN(2,1)*clhs84 + DN(2,2)*clhs87 + clhs328*clhs423 + clhs345 + clhs378;
    lhs(8,6)=DN(2,0)*clhs94 + DN(2,1)*clhs96 + DN(2,2)*clhs98 + clhs331*clhs423 + clhs351 + clhs394;
    lhs(8,7)=clhs100*clhs423 + clhs338*clhs99 - clhs338;
    lhs(8,8)=DN(2,0)*clhs106 + DN(2,1)*clhs108 + DN(2,2)*clhs110 + clhs112*clhs422 + clhs114*clhs424 + clhs31*clhs428 + clhs425*tau2 + clhs427;
    lhs(8,9)=DN(2,0)*clhs121 + DN(2,1)*clhs123 + DN(2,2)*clhs126 + clhs29*clhs431 + clhs31*clhs432 + clhs347*clhs423 + clhs430;
    lhs(8,10)=DN(2,0)*clhs133 + DN(2,1)*clhs135 + DN(2,2)*clhs137 + clhs30*clhs431 + clhs31*clhs434 + clhs353*clhs423 + clhs433;
    lhs(8,11)=DN(2,0)*clhs436;
    lhs(8,12)=DN(2,0)*clhs144 + DN(2,1)*clhs146 + DN(2,2)*clhs148 + clhs150*clhs422 + clhs152*clhs424 + clhs442;
    lhs(8,13)=DN(2,0)*clhs159 + DN(2,1)*clhs161 + DN(2,2)*clhs164 + clhs365*clhs423 + clhs447 + clhs448;
    lhs(8,14)=DN(2,0)*clhs171 + DN(2,1)*clhs173 + DN(2,2)*clhs175 + clhs371*clhs423 + clhs452 + clhs453;
    lhs(8,15)=clhs176*clhs423 + clhs439*clhs99 - clhs439;
    lhs(9,0)=DN(2,0)*clhs6 + DN(2,1)*clhs177 + DN(2,2)*clhs178 + clhs120 + clhs180*clhs423 + clhs212;
    lhs(9,1)=DN(2,0)*clhs42 + DN(2,1)*clhs185 + DN(2,2)*clhs187 + clhs188*clhs422 + clhs189*clhs424 + clhs217;
    lhs(9,2)=DN(2,0)*clhs52 + DN(2,1)*clhs192 + DN(2,2)*clhs194 + clhs223 + clhs275 + clhs372*clhs423;
    lhs(9,3)=clhs117*clhs99 - clhs117 + clhs298*clhs423;
    lhs(9,4)=DN(2,0)*clhs69 + DN(2,1)*clhs197 + DN(2,2)*clhs198 + clhs199*clhs423 + clhs346 + clhs377;
    lhs(9,5)=DN(2,0)*clhs84 + DN(2,1)*clhs201 + DN(2,2)*clhs203 + clhs204*clhs422 + clhs205*clhs424 + clhs379;
    lhs(9,6)=DN(2,0)*clhs96 + DN(2,1)*clhs208 + DN(2,2)*clhs210 + clhs376*clhs423 + clhs380 + clhs396;
    lhs(9,7)=clhs211*clhs423 + clhs343*clhs99 - clhs343;
    lhs(9,8)=DN(2,0)*clhs108 + DN(2,1)*clhs214 + DN(2,2)*clhs215 + clhs179*clhs431 + clhs183*clhs428 + clhs216*clhs423 + clhs430;
    lhs(9,9)=DN(2,0)*clhs123 + DN(2,1)*clhs218 + DN(2,2)*clhs220 + clhs183*clhs432 + clhs221*clhs422 + clhs222*clhs424 + clhs427 + clhs454*tau2;
    lhs(9,10)=DN(2,0)*clhs135 + DN(2,1)*clhs225 + DN(2,2)*clhs227 + clhs182*clhs431 + clhs183*clhs434 + clhs382*clhs423 + clhs456;
    lhs(9,11)=DN(2,1)*clhs436;
    lhs(9,12)=DN(2,0)*clhs146 + DN(2,1)*clhs231 + DN(2,2)*clhs232 + clhs233*clhs423 + clhs457 + clhs458;
    lhs(9,13)=DN(2,0)*clhs161 + DN(2,1)*clhs235 + DN(2,2)*clhs237 + clhs238*clhs422 + clhs239*clhs424 + clhs459;
    lhs(9,14)=DN(2,0)*clhs173 + DN(2,1)*clhs242 + DN(2,2)*clhs244 + clhs388*clhs423 + clhs460 + clhs461;
    lhs(9,15)=clhs245*clhs423 + clhs444*clhs99 - clhs444;
    lhs(10,0)=DN(2,0)*clhs8 + DN(2,1)*clhs178 + DN(2,2)*clhs246 + clhs132 + clhs248*clhs423 + clhs270;
    lhs(10,1)=DN(2,0)*clhs45 + DN(2,1)*clhs187 + DN(2,2)*clhs252 + clhs224 + clhs274 + clhs389*clhs423;
    lhs(10,2)=DN(2,0)*clhs54 + DN(2,1)*clhs194 + DN(2,2)*clhs254 + clhs255*clhs422 + clhs256*clhs424 + clhs277;
    lhs(10,3)=clhs129*clhs99 - clhs129 + clhs299*clhs423;
    lhs(10,4)=DN(2,0)*clhs71 + DN(2,1)*clhs198 + DN(2,2)*clhs260 + clhs261*clhs423 + clhs352 + clhs392;
    lhs(10,5)=DN(2,0)*clhs87 + DN(2,1)*clhs203 + DN(2,2)*clhs264 + clhs381 + clhs390*clhs423 + clhs395;
    lhs(10,6)=DN(2,0)*clhs98 + DN(2,1)*clhs210 + DN(2,2)*clhs266 + clhs267*clhs422 + clhs268*clhs424 + clhs398;
    lhs(10,7)=clhs269*clhs423 + clhs349*clhs99 - clhs349;
    lhs(10,8)=DN(2,0)*clhs110 + DN(2,1)*clhs215 + DN(2,2)*clhs272 + clhs247*clhs431 + clhs251*clhs428 + clhs273*clhs423 + clhs433;
    lhs(10,9)=DN(2,0)*clhs126 + DN(2,1)*clhs220 + DN(2,2)*clhs276 + clhs249*clhs431 + clhs251*clhs432 + clhs397*clhs423 + clhs456;
    lhs(10,10)=DN(2,0)*clhs137 + DN(2,1)*clhs227 + DN(2,2)*clhs278 + clhs251*clhs434 + clhs279*clhs422 + clhs280*clhs424 + clhs427 + clhs462*tau2;
    lhs(10,11)=DN(2,2)*clhs436;
    lhs(10,12)=DN(2,0)*clhs148 + DN(2,1)*clhs232 + DN(2,2)*clhs284 + clhs285*clhs423 + clhs463 + clhs465;
    lhs(10,13)=DN(2,0)*clhs164 + DN(2,1)*clhs237 + DN(2,2)*clhs288 + clhs403*clhs423 + clhs466 + clhs467;
    lhs(10,14)=DN(2,0)*clhs175 + DN(2,1)*clhs244 + DN(2,2)*clhs290 + clhs291*clhs422 + clhs292*clhs424 + clhs468;
    lhs(10,15)=clhs293*clhs423 + clhs449*clhs99 - clhs449;
    lhs(11,0)=clhs102 + clhs138*clhs33 + clhs179*clhs469 + clhs247*clhs470;
    lhs(11,1)=clhs116 + clhs189*clhs228 + clhs249*clhs470 + clhs29*clhs471;
    lhs(11,2)=clhs128 + clhs182*clhs469 + clhs256*clhs281 + clhs30*clhs471;
    lhs(11,3)=clhs311;
    lhs(11,4)=clhs138*clhs75 + clhs179*clhs472 + clhs247*clhs473 + clhs337;
    lhs(11,5)=clhs205*clhs228 + clhs249*clhs473 + clhs29*clhs474 + clhs342;
    lhs(11,6)=clhs182*clhs472 + clhs268*clhs281 + clhs30*clhs474 + clhs348;
    lhs(11,7)=clhs417;
    lhs(11,8)=DN(2,0)*N[2] + clhs114*clhs138 + clhs179*clhs475 + clhs247*clhs476;
    lhs(11,9)=DN(2,1)*N[2] + clhs222*clhs228 + clhs249*clhs476 + clhs29*clhs477;
    lhs(11,10)=DN(2,2)*N[2] + clhs182*clhs475 + clhs280*clhs281 + clhs30*clhs477;
    lhs(11,11)=clhs300*clhs426 + clhs425*tau1 + clhs454*tau1 + clhs462*tau1;
    lhs(11,12)=clhs138*clhs152 + clhs179*clhs478 + clhs247*clhs479 + clhs440;
    lhs(11,13)=clhs228*clhs239 + clhs249*clhs479 + clhs29*clhs480 + clhs445;
    lhs(11,14)=clhs182*clhs478 + clhs281*clhs292 + clhs30*clhs480 + clhs450;
    lhs(11,15)=clhs481;
    lhs(12,0)=DN(3,0)*clhs4 + DN(3,1)*clhs6 + DN(3,2)*clhs8 + clhs143 + clhs27*clhs482 + clhs33*clhs484;
    lhs(12,1)=DN(3,0)*clhs40 + DN(3,1)*clhs42 + DN(3,2)*clhs45 + clhs157 + clhs230 + clhs319*clhs483;
    lhs(12,2)=DN(3,0)*clhs50 + DN(3,1)*clhs52 + DN(3,2)*clhs54 + clhs169 + clhs283 + clhs320*clhs483;
    lhs(12,3)=clhs141*clhs99 - clhs141 + clhs296*clhs483;
    lhs(12,4)=DN(3,0)*clhs67 + DN(3,1)*clhs69 + DN(3,2)*clhs71 + clhs358 + clhs482*clhs73 + clhs484*clhs75;
    lhs(12,5)=DN(3,0)*clhs82 + DN(3,1)*clhs84 + DN(3,2)*clhs87 + clhs328*clhs483 + clhs363 + clhs384;
    lhs(12,6)=DN(3,0)*clhs94 + DN(3,1)*clhs96 + DN(3,2)*clhs98 + clhs331*clhs483 + clhs369 + clhs400;
    lhs(12,7)=clhs100*clhs483 + clhs356*clhs99 - clhs356;
    lhs(12,8)=DN(3,0)*clhs106 + DN(3,1)*clhs108 + DN(3,2)*clhs110 + clhs112*clhs482 + clhs114*clhs484 + clhs442;
    lhs(12,9)=DN(3,0)*clhs121 + DN(3,1)*clhs123 + DN(3,2)*clhs126 + clhs347*clhs483 + clhs447 + clhs458;
    lhs(12,10)=DN(3,0)*clhs133 + DN(3,1)*clhs135 + DN(3,2)*clhs137 + clhs353*clhs483 + clhs452 + clhs465;
    lhs(12,11)=clhs138*clhs483 + clhs440*clhs99 - clhs440;
    lhs(12,12)=DN(3,0)*clhs144 + DN(3,1)*clhs146 + DN(3,2)*clhs148 + clhs150*clhs482 + clhs152*clhs484 + clhs31*clhs488 + clhs485*tau2 + clhs487;
    lhs(12,13)=DN(3,0)*clhs159 + DN(3,1)*clhs161 + DN(3,2)*clhs164 + clhs29*clhs491 + clhs31*clhs492 + clhs365*clhs483 + clhs490;
    lhs(12,14)=DN(3,0)*clhs171 + DN(3,1)*clhs173 + DN(3,2)*clhs175 + clhs30*clhs491 + clhs31*clhs494 + clhs371*clhs483 + clhs493;
    lhs(12,15)=DN(3,0)*clhs495;
    lhs(13,0)=DN(3,0)*clhs6 + DN(3,1)*clhs177 + DN(3,2)*clhs178 + clhs158 + clhs180*clhs483 + clhs229;
    lhs(13,1)=DN(3,0)*clhs42 + DN(3,1)*clhs185 + DN(3,2)*clhs187 + clhs188*clhs482 + clhs189*clhs484 + clhs234;
    lhs(13,2)=DN(3,0)*clhs52 + DN(3,1)*clhs192 + DN(3,2)*clhs194 + clhs240 + clhs287 + clhs372*clhs483;
    lhs(13,3)=clhs155*clhs99 - clhs155 + clhs298*clhs483;
    lhs(13,4)=DN(3,0)*clhs69 + DN(3,1)*clhs197 + DN(3,2)*clhs198 + clhs199*clhs483 + clhs364 + clhs383;
    lhs(13,5)=DN(3,0)*clhs84 + DN(3,1)*clhs201 + DN(3,2)*clhs203 + clhs204*clhs482 + clhs205*clhs484 + clhs385;
    lhs(13,6)=DN(3,0)*clhs96 + DN(3,1)*clhs208 + DN(3,2)*clhs210 + clhs376*clhs483 + clhs386 + clhs402;
    lhs(13,7)=clhs211*clhs483 + clhs361*clhs99 - clhs361;
    lhs(13,8)=DN(3,0)*clhs108 + DN(3,1)*clhs214 + DN(3,2)*clhs215 + clhs216*clhs483 + clhs448 + clhs457;
    lhs(13,9)=DN(3,0)*clhs123 + DN(3,1)*clhs218 + DN(3,2)*clhs220 + clhs221*clhs482 + clhs222*clhs484 + clhs459;
    lhs(13,10)=DN(3,0)*clhs135 + DN(3,1)*clhs225 + DN(3,2)*clhs227 + clhs382*clhs483 + clhs460 + clhs467;
    lhs(13,11)=clhs228*clhs483 + clhs445*clhs99 - clhs445;
    lhs(13,12)=DN(3,0)*clhs146 + DN(3,1)*clhs231 + DN(3,2)*clhs232 + clhs179*clhs491 + clhs183*clhs488 + clhs233*clhs483 + clhs490;
    lhs(13,13)=DN(3,0)*clhs161 + DN(3,1)*clhs235 + DN(3,2)*clhs237 + clhs183*clhs492 + clhs238*clhs482 + clhs239*clhs484 + clhs487 + clhs496*tau2;
    lhs(13,14)=DN(3,0)*clhs173 + DN(3,1)*clhs242 + DN(3,2)*clhs244 + clhs182*clhs491 + clhs183*clhs494 + clhs388*clhs483 + clhs497;
    lhs(13,15)=DN(3,1)*clhs495;
    lhs(14,0)=DN(3,0)*clhs8 + DN(3,1)*clhs178 + DN(3,2)*clhs246 + clhs170 + clhs248*clhs483 + clhs282;
    lhs(14,1)=DN(3,0)*clhs45 + DN(3,1)*clhs187 + DN(3,2)*clhs252 + clhs241 + clhs286 + clhs389*clhs483;
    lhs(14,2)=DN(3,0)*clhs54 + DN(3,1)*clhs194 + DN(3,2)*clhs254 + clhs255*clhs482 + clhs256*clhs484 + clhs289;
    lhs(14,3)=clhs167*clhs99 - clhs167 + clhs299*clhs483;
    lhs(14,4)=DN(3,0)*clhs71 + DN(3,1)*clhs198 + DN(3,2)*clhs260 + clhs261*clhs483 + clhs370 + clhs399;
    lhs(14,5)=DN(3,0)*clhs87 + DN(3,1)*clhs203 + DN(3,2)*clhs264 + clhs387 + clhs390*clhs483 + clhs401;
    lhs(14,6)=DN(3,0)*clhs98 + DN(3,1)*clhs210 + DN(3,2)*clhs266 + clhs267*clhs482 + clhs268*clhs484 + clhs404;
    lhs(14,7)=clhs269*clhs483 + clhs367*clhs99 - clhs367;
    lhs(14,8)=DN(3,0)*clhs110 + DN(3,1)*clhs215 + DN(3,2)*clhs272 + clhs273*clhs483 + clhs453 + clhs463;
    lhs(14,9)=DN(3,0)*clhs126 + DN(3,1)*clhs220 + DN(3,2)*clhs276 + clhs397*clhs483 + clhs461 + clhs466;
    lhs(14,10)=DN(3,0)*clhs137 + DN(3,1)*clhs227 + DN(3,2)*clhs278 + clhs279*clhs482 + clhs280*clhs484 + clhs468;
    lhs(14,11)=clhs281*clhs483 + clhs450*clhs99 - clhs450;
    lhs(14,12)=DN(3,0)*clhs148 + DN(3,1)*clhs232 + DN(3,2)*clhs284 + clhs247*clhs491 + clhs251*clhs488 + clhs285*clhs483 + clhs493;
    lhs(14,13)=DN(3,0)*clhs164 + DN(3,1)*clhs237 + DN(3,2)*clhs288 + clhs249*clhs491 + clhs251*clhs492 + clhs403*clhs483 + clhs497;
    lhs(14,14)=DN(3,0)*clhs175 + DN(3,1)*clhs244 + DN(3,2)*clhs290 + clhs251*clhs494 + clhs291*clhs482 + clhs292*clhs484 + clhs487 + clhs498*tau2;
    lhs(14,15)=DN(3,2)*clhs495;
    lhs(15,0)=clhs140 + clhs176*clhs33 + clhs179*clhs499 + clhs247*clhs500;
    lhs(15,1)=clhs154 + clhs189*clhs245 + clhs249*clhs500 + clhs29*clhs501;
    lhs(15,2)=clhs166 + clhs182*clhs499 + clhs256*clhs293 + clhs30*clhs501;
    lhs(15,3)=clhs315;
    lhs(15,4)=clhs176*clhs75 + clhs179*clhs502 + clhs247*clhs503 + clhs355;
    lhs(15,5)=clhs205*clhs245 + clhs249*clhs503 + clhs29*clhs504 + clhs360;
    lhs(15,6)=clhs182*clhs502 + clhs268*clhs293 + clhs30*clhs504 + clhs366;
    lhs(15,7)=clhs421;
    lhs(15,8)=clhs114*clhs176 + clhs179*clhs505 + clhs247*clhs506 + clhs439;
    lhs(15,9)=clhs222*clhs245 + clhs249*clhs506 + clhs29*clhs507 + clhs444;
    lhs(15,10)=clhs182*clhs505 + clhs280*clhs293 + clhs30*clhs507 + clhs449;
    lhs(15,11)=clhs481;
    lhs(15,12)=DN(3,0)*N[3] + clhs152*clhs176 + clhs179*clhs508 + clhs247*clhs509;
    lhs(15,13)=DN(3,1)*N[3] + clhs239*clhs245 + clhs249*clhs509 + clhs29*clhs510;
    lhs(15,14)=DN(3,2)*N[3] + clhs182*clhs508 + clhs292*clhs293 + clhs30*clhs510;
    lhs(15,15)=clhs300*clhs486 + clhs485*tau1 + clhs496*tau1 + clhs498*tau1;

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
    // const array_1d<double,nnodes>& pn = data.pn;
    // const array_1d<double,nnodes>& pnn = data.pnn;
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
    const double clhs9 =             v(0,0) - vmesh(0,0);
    const double clhs10 =             v(1,0) - vmesh(1,0);
    const double clhs11 =             v(2,0) - vmesh(2,0);
    const double clhs12 =             N[0]*clhs9 + N[1]*clhs10 + N[2]*clhs11;
    const double clhs13 =             v(0,1) - vmesh(0,1);
    const double clhs14 =             v(1,1) - vmesh(1,1);
    const double clhs15 =             v(2,1) - vmesh(2,1);
    const double clhs16 =             N[0]*clhs13 + N[1]*clhs14 + N[2]*clhs15;
    const double clhs17 =             DN(0,0)*clhs12 + DN(0,1)*clhs16;
    const double clhs18 =             N[0]*clhs8 + clhs17;
    const double clhs19 =             2*DN(0,0)*N[0]*rho*tau1;
    const double clhs20 =             DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0);
    const double clhs21 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + clhs12*clhs8 + clhs16*clhs20);
    const double clhs22 =             N[0]*bdf0;
    const double clhs23 =             clhs18 + clhs22;
    const double clhs24 =             pow(rho, 2);
    const double clhs25 =             DN(0,0)*clhs9 + DN(0,1)*clhs13 + DN(1,0)*clhs10 + DN(1,1)*clhs14 + DN(2,0)*clhs11 + DN(2,1)*clhs15;
    const double clhs26 =             N[0]*clhs25 + clhs17;
    const double clhs27 =             clhs24*clhs26*tau1;
    const double clhs28 =             DN(0,0)*tau2;
    const double clhs29 =             DN(0,1)*clhs28;
    const double clhs30 =             C(0,1)*DN(0,1) + clhs5;
    const double clhs31 =             C(1,2)*DN(0,1);
    const double clhs32 =             C(2,2)*DN(0,0) + clhs31;
    const double clhs33 =             clhs1*rho;
    const double clhs34 =             N[0]*clhs20*clhs24*tau1;
    const double clhs35 =             2*DN(0,1)*N[0]*rho*tau1;
    const double clhs36 =             pow(c, -2);
    const double clhs37 =             1.0/rho;
    const double clhs38 =             N[0]*bdf0*clhs36*clhs37;
    const double clhs39 =             rho*tau1;
    const double clhs40 =             -N[0] + clhs26*clhs39 + clhs38*tau2;
    const double clhs41 =             N[0]*bdf0*rho;
    const double clhs42 =             N[1]*clhs41;
    const double clhs43 =             DN(0,0)*N[1];
    const double clhs44 =             DN(1,0)*N[0];
    const double clhs45 =             rho*tau1*(clhs43 + clhs44);
    const double clhs46 =             DN(1,0)*clhs28 + clhs21*clhs45 + clhs42;
    const double clhs47 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
    const double clhs48 =             C(0,2)*DN(1,0);
    const double clhs49 =             C(2,2)*DN(1,1) + clhs48;
    const double clhs50 =             DN(1,0)*clhs12 + DN(1,1)*clhs16;
    const double clhs51 =             N[1]*clhs8 + clhs50;
    const double clhs52 =             N[1]*bdf0;
    const double clhs53 =             clhs51 + clhs52;
    const double clhs54 =             DN(1,1)*clhs28;
    const double clhs55 =             C(0,1)*DN(1,1) + clhs48;
    const double clhs56 =             C(1,2)*DN(1,1);
    const double clhs57 =             C(2,2)*DN(1,0) + clhs56;
    const double clhs58 =             N[1]*clhs20*clhs24*tau1;
    const double clhs59 =             N[0]*N[1]*rho;
    const double clhs60 =             DN(0,1)*N[1];
    const double clhs61 =             DN(1,1)*N[0];
    const double clhs62 =             rho*tau1*(clhs60 + clhs61);
    const double clhs63 =             clhs20*clhs59 + clhs21*clhs62;
    const double clhs64 =             bdf0*clhs36*clhs37*tau2;
    const double clhs65 =             DN(1,0)*rho*tau1;
    const double clhs66 =             N[2]*clhs41;
    const double clhs67 =             DN(0,0)*N[2];
    const double clhs68 =             DN(2,0)*N[0];
    const double clhs69 =             rho*tau1*(clhs67 + clhs68);
    const double clhs70 =             DN(2,0)*clhs28 + clhs21*clhs69 + clhs66;
    const double clhs71 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
    const double clhs72 =             C(0,2)*DN(2,0);
    const double clhs73 =             C(2,2)*DN(2,1) + clhs72;
    const double clhs74 =             DN(2,0)*clhs12 + DN(2,1)*clhs16;
    const double clhs75 =             N[2]*clhs8 + clhs74;
    const double clhs76 =             N[2]*bdf0;
    const double clhs77 =             clhs75 + clhs76;
    const double clhs78 =             DN(2,1)*clhs28;
    const double clhs79 =             C(0,1)*DN(2,1) + clhs72;
    const double clhs80 =             C(1,2)*DN(2,1);
    const double clhs81 =             C(2,2)*DN(2,0) + clhs80;
    const double clhs82 =             N[2]*clhs20*clhs24*tau1;
    const double clhs83 =             N[0]*N[2]*rho;
    const double clhs84 =             DN(0,1)*N[2];
    const double clhs85 =             DN(2,1)*N[0];
    const double clhs86 =             rho*tau1*(clhs84 + clhs85);
    const double clhs87 =             clhs20*clhs83 + clhs21*clhs86;
    const double clhs88 =             DN(2,0)*rho*tau1;
    const double clhs89 =             C(0,1)*DN(0,0) + clhs31;
    const double clhs90 =             DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1);
    const double clhs91 =             N[0]*clhs24*clhs90*tau1;
    const double clhs92 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
    const double clhs93 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + clhs12*clhs90 + clhs16*clhs92);
    const double clhs94 =             pow(DN(0,1), 2);
    const double clhs95 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
    const double clhs96 =             N[0]*clhs92 + clhs17;
    const double clhs97 =             clhs22 + clhs96;
    const double clhs98 =             DN(0,1)*tau2;
    const double clhs99 =             DN(1,0)*clhs98;
    const double clhs100 =             C(0,1)*DN(1,0) + clhs56;
    const double clhs101 =             N[1]*clhs24*clhs90*tau1;
    const double clhs102 =             clhs45*clhs93 + clhs59*clhs90;
    const double clhs103 =             DN(1,1)*clhs98 + clhs42 + clhs62*clhs93;
    const double clhs104 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
    const double clhs105 =             N[1]*clhs92 + clhs50;
    const double clhs106 =             clhs105 + clhs52;
    const double clhs107 =             DN(1,1)*rho*tau1;
    const double clhs108 =             DN(2,0)*clhs98;
    const double clhs109 =             C(0,1)*DN(2,0) + clhs80;
    const double clhs110 =             N[2]*clhs24*clhs90*tau1;
    const double clhs111 =             clhs69*clhs93 + clhs83*clhs90;
    const double clhs112 =             DN(2,1)*clhs98 + clhs66 + clhs86*clhs93;
    const double clhs113 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
    const double clhs114 =             N[2]*clhs92 + clhs74;
    const double clhs115 =             clhs114 + clhs76;
    const double clhs116 =             DN(2,1)*rho*tau1;
    const double clhs117 =             DN(0,0)*N[0];
    const double clhs118 =             DN(0,1)*N[0];
    const double clhs119 =             clhs90*rho*tau1;
    const double clhs120 =             DN(0,0)*rho*tau1;
    const double clhs121 =             clhs20*rho*tau1;
    const double clhs122 =             DN(0,1)*rho*tau1;
    const double clhs123 =             bdf0*clhs36*clhs37;
    const double clhs124 =             DN(0,0)*tau1;
    const double clhs125 =             DN(0,1)*tau1;
    const double clhs126 =             DN(1,0)*clhs124 + DN(1,1)*clhs125 + N[1]*clhs38;
    const double clhs127 =             DN(2,0)*clhs124 + DN(2,1)*clhs125 + N[2]*clhs38;
    const double clhs128 =             N[1]*rho;
    const double clhs129 =             N[1]*clhs25 + clhs50;
    const double clhs130 =             clhs129*clhs24*tau1;
    const double clhs131 =             pow(DN(1,0), 2);
    const double clhs132 =             pow(N[1], 2);
    const double clhs133 =             clhs132*clhs2;
    const double clhs134 =             2*DN(1,0)*N[1]*rho*tau1;
    const double clhs135 =             DN(1,0)*tau2;
    const double clhs136 =             DN(1,1)*clhs135;
    const double clhs137 =             clhs132*rho;
    const double clhs138 =             2*DN(1,1)*N[1]*rho*tau1;
    const double clhs139 =             clhs36*clhs37*tau2;
    const double clhs140 =             -N[1] + clhs129*clhs39 + clhs139*clhs52;
    const double clhs141 =             N[1]*N[2]*bdf0;
    const double clhs142 =             clhs141*rho;
    const double clhs143 =             DN(1,0)*N[2];
    const double clhs144 =             DN(2,0)*N[1];
    const double clhs145 =             rho*tau1*(clhs143 + clhs144);
    const double clhs146 =             DN(2,0)*clhs135 + clhs142 + clhs145*clhs21;
    const double clhs147 =             DN(2,1)*clhs135;
    const double clhs148 =             N[1]*N[2]*rho;
    const double clhs149 =             DN(1,1)*N[2];
    const double clhs150 =             DN(2,1)*N[1];
    const double clhs151 =             rho*tau1*(clhs149 + clhs150);
    const double clhs152 =             clhs148*clhs20 + clhs151*clhs21;
    const double clhs153 =             pow(DN(1,1), 2);
    const double clhs154 =             DN(1,1)*tau2;
    const double clhs155 =             DN(2,0)*clhs154;
    const double clhs156 =             clhs145*clhs93 + clhs148*clhs90;
    const double clhs157 =             DN(2,1)*clhs154 + clhs142 + clhs151*clhs93;
    const double clhs158 =             DN(1,0)*N[1];
    const double clhs159 =             DN(1,1)*N[1];
    const double clhs160 =             DN(2,0)*(DN(1,0)*tau1) + DN(2,1)*(DN(1,1)*tau1) + clhs141*clhs36*clhs37;
    const double clhs161 =             N[2]*rho;
    const double clhs162 =             N[2]*clhs25 + clhs74;
    const double clhs163 =             clhs162*clhs24*tau1;
    const double clhs164 =             pow(DN(2,0), 2);
    const double clhs165 =             pow(N[2], 2);
    const double clhs166 =             clhs165*clhs2;
    const double clhs167 =             2*DN(2,0)*N[2]*rho*tau1;
    const double clhs168 =             DN(2,0)*DN(2,1)*tau2;
    const double clhs169 =             clhs165*rho;
    const double clhs170 =             2*DN(2,1)*N[2]*rho*tau1;
    const double clhs171 =             -N[2] + clhs139*clhs76 + clhs162*clhs39;
    const double clhs172 =             pow(DN(2,1), 2);
    const double clhs173 =             DN(2,0)*N[2];
    const double clhs174 =             DN(2,1)*N[2];

    lhs(0,0)=DN(0,0)*clhs4 + DN(0,1)*clhs6 + clhs0*tau2 + clhs18*clhs7 + clhs19*clhs21 + clhs23*clhs27 + clhs3;
    lhs(0,1)=DN(0,0)*clhs30 + DN(0,1)*clhs32 + clhs20*clhs33 + clhs21*clhs35 + clhs26*clhs34 + clhs29;
    lhs(0,2)=DN(0,0)*clhs40;
    lhs(0,3)=DN(0,0)*clhs47 + DN(0,1)*clhs49 + clhs27*clhs53 + clhs46 + clhs51*clhs7;
    lhs(0,4)=DN(0,0)*clhs55 + DN(0,1)*clhs57 + clhs26*clhs58 + clhs54 + clhs63;
    lhs(0,5)=clhs26*clhs65 + clhs43*clhs64 - clhs43;
    lhs(0,6)=DN(0,0)*clhs71 + DN(0,1)*clhs73 + clhs27*clhs77 + clhs7*clhs75 + clhs70;
    lhs(0,7)=DN(0,0)*clhs79 + DN(0,1)*clhs81 + clhs26*clhs82 + clhs78 + clhs87;
    lhs(0,8)=clhs26*clhs88 + clhs64*clhs67 - clhs67;
    lhs(1,0)=DN(0,0)*clhs6 + DN(0,1)*clhs89 + clhs19*clhs93 + clhs26*clhs91 + clhs29 + clhs33*clhs90;
    lhs(1,1)=DN(0,0)*clhs32 + DN(0,1)*clhs95 + clhs27*clhs97 + clhs3 + clhs35*clhs93 + clhs7*clhs96 + clhs94*tau2;
    lhs(1,2)=DN(0,1)*clhs40;
    lhs(1,3)=DN(0,0)*clhs49 + DN(0,1)*clhs100 + clhs101*clhs26 + clhs102 + clhs99;
    lhs(1,4)=DN(0,0)*clhs57 + DN(0,1)*clhs104 + clhs103 + clhs105*clhs7 + clhs106*clhs27;
    lhs(1,5)=clhs107*clhs26 + clhs60*clhs64 - clhs60;
    lhs(1,6)=DN(0,0)*clhs73 + DN(0,1)*clhs109 + clhs108 + clhs110*clhs26 + clhs111;
    lhs(1,7)=DN(0,0)*clhs81 + DN(0,1)*clhs113 + clhs112 + clhs114*clhs7 + clhs115*clhs27;
    lhs(1,8)=clhs116*clhs26 + clhs64*clhs84 - clhs84;
    lhs(2,0)=clhs117 + clhs118*clhs119 + clhs120*clhs23;
    lhs(2,1)=clhs117*clhs121 + clhs118 + clhs122*clhs97;
    lhs(2,2)=clhs0*tau1 + clhs1*clhs123 + clhs94*tau1;
    lhs(2,3)=clhs119*clhs60 + clhs120*clhs53 + clhs44;
    lhs(2,4)=clhs106*clhs122 + clhs121*clhs43 + clhs61;
    lhs(2,5)=clhs126;
    lhs(2,6)=clhs119*clhs84 + clhs120*clhs77 + clhs68;
    lhs(2,7)=clhs115*clhs122 + clhs121*clhs67 + clhs85;
    lhs(2,8)=clhs127;
    lhs(3,0)=DN(1,0)*clhs4 + DN(1,1)*clhs6 + clhs128*clhs18 + clhs130*clhs23 + clhs46;
    lhs(3,1)=DN(1,0)*clhs30 + DN(1,1)*clhs32 + clhs129*clhs34 + clhs63 + clhs99;
    lhs(3,2)=clhs120*clhs129 + clhs44*clhs64 - clhs44;
    lhs(3,3)=DN(1,0)*clhs47 + DN(1,1)*clhs49 + clhs128*clhs51 + clhs130*clhs53 + clhs131*tau2 + clhs133 + clhs134*clhs21;
    lhs(3,4)=DN(1,0)*clhs55 + DN(1,1)*clhs57 + clhs129*clhs58 + clhs136 + clhs137*clhs20 + clhs138*clhs21;
    lhs(3,5)=DN(1,0)*clhs140;
    lhs(3,6)=DN(1,0)*clhs71 + DN(1,1)*clhs73 + clhs128*clhs75 + clhs130*clhs77 + clhs146;
    lhs(3,7)=DN(1,0)*clhs79 + DN(1,1)*clhs81 + clhs129*clhs82 + clhs147 + clhs152;
    lhs(3,8)=clhs129*clhs88 + clhs143*clhs64 - clhs143;
    lhs(4,0)=DN(1,0)*clhs6 + DN(1,1)*clhs89 + clhs102 + clhs129*clhs91 + clhs54;
    lhs(4,1)=DN(1,0)*clhs32 + DN(1,1)*clhs95 + clhs103 + clhs128*clhs96 + clhs130*clhs97;
    lhs(4,2)=clhs122*clhs129 + clhs61*clhs64 - clhs61;
    lhs(4,3)=DN(1,0)*clhs49 + DN(1,1)*clhs100 + clhs101*clhs129 + clhs134*clhs93 + clhs136 + clhs137*clhs90;
    lhs(4,4)=DN(1,0)*clhs57 + DN(1,1)*clhs104 + clhs105*clhs128 + clhs106*clhs130 + clhs133 + clhs138*clhs93 + clhs153*tau2;
    lhs(4,5)=DN(1,1)*clhs140;
    lhs(4,6)=DN(1,0)*clhs73 + DN(1,1)*clhs109 + clhs110*clhs129 + clhs155 + clhs156;
    lhs(4,7)=DN(1,0)*clhs81 + DN(1,1)*clhs113 + clhs114*clhs128 + clhs115*clhs130 + clhs157;
    lhs(4,8)=clhs116*clhs129 + clhs149*clhs64 - clhs149;
    lhs(5,0)=clhs119*clhs61 + clhs23*clhs65 + clhs43;
    lhs(5,1)=clhs107*clhs97 + clhs121*clhs44 + clhs60;
    lhs(5,2)=clhs126;
    lhs(5,3)=clhs119*clhs159 + clhs158 + clhs53*clhs65;
    lhs(5,4)=clhs106*clhs107 + clhs121*clhs158 + clhs159;
    lhs(5,5)=clhs123*clhs132 + clhs131*tau1 + clhs153*tau1;
    lhs(5,6)=clhs119*clhs149 + clhs144 + clhs65*clhs77;
    lhs(5,7)=clhs107*clhs115 + clhs121*clhs143 + clhs150;
    lhs(5,8)=clhs160;
    lhs(6,0)=DN(2,0)*clhs4 + DN(2,1)*clhs6 + clhs161*clhs18 + clhs163*clhs23 + clhs70;
    lhs(6,1)=DN(2,0)*clhs30 + DN(2,1)*clhs32 + clhs108 + clhs162*clhs34 + clhs87;
    lhs(6,2)=clhs120*clhs162 + clhs64*clhs68 - clhs68;
    lhs(6,3)=DN(2,0)*clhs47 + DN(2,1)*clhs49 + clhs146 + clhs161*clhs51 + clhs163*clhs53;
    lhs(6,4)=DN(2,0)*clhs55 + DN(2,1)*clhs57 + clhs152 + clhs155 + clhs162*clhs58;
    lhs(6,5)=clhs144*clhs64 - clhs144 + clhs162*clhs65;
    lhs(6,6)=DN(2,0)*clhs71 + DN(2,1)*clhs73 + clhs161*clhs75 + clhs163*clhs77 + clhs164*tau2 + clhs166 + clhs167*clhs21;
    lhs(6,7)=DN(2,0)*clhs79 + DN(2,1)*clhs81 + clhs162*clhs82 + clhs168 + clhs169*clhs20 + clhs170*clhs21;
    lhs(6,8)=DN(2,0)*clhs171;
    lhs(7,0)=DN(2,0)*clhs6 + DN(2,1)*clhs89 + clhs111 + clhs162*clhs91 + clhs78;
    lhs(7,1)=DN(2,0)*clhs32 + DN(2,1)*clhs95 + clhs112 + clhs161*clhs96 + clhs163*clhs97;
    lhs(7,2)=clhs122*clhs162 + clhs64*clhs85 - clhs85;
    lhs(7,3)=DN(2,0)*clhs49 + DN(2,1)*clhs100 + clhs101*clhs162 + clhs147 + clhs156;
    lhs(7,4)=DN(2,0)*clhs57 + DN(2,1)*clhs104 + clhs105*clhs161 + clhs106*clhs163 + clhs157;
    lhs(7,5)=clhs107*clhs162 + clhs150*clhs64 - clhs150;
    lhs(7,6)=DN(2,0)*clhs73 + DN(2,1)*clhs109 + clhs110*clhs162 + clhs167*clhs93 + clhs168 + clhs169*clhs90;
    lhs(7,7)=DN(2,0)*clhs81 + DN(2,1)*clhs113 + clhs114*clhs161 + clhs115*clhs163 + clhs166 + clhs170*clhs93 + clhs172*tau2;
    lhs(7,8)=DN(2,1)*clhs171;
    lhs(8,0)=clhs119*clhs85 + clhs23*clhs88 + clhs67;
    lhs(8,1)=clhs116*clhs97 + clhs121*clhs68 + clhs84;
    lhs(8,2)=clhs127;
    lhs(8,3)=clhs119*clhs150 + clhs143 + clhs53*clhs88;
    lhs(8,4)=clhs106*clhs116 + clhs121*clhs144 + clhs149;
    lhs(8,5)=clhs160;
    lhs(8,6)=clhs119*clhs174 + clhs173 + clhs77*clhs88;
    lhs(8,7)=clhs115*clhs116 + clhs121*clhs173 + clhs174;
    lhs(8,8)=clhs123*clhs165 + clhs164*tau1 + clhs172*tau1;

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
    const double crhs2 =             N[0]*rho;
    const double crhs3 =             N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0));
    const double crhs4 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
    const double crhs5 =             DN(0,0)*v(0,0);
    const double crhs6 =             DN(0,1)*v(0,1);
    const double crhs7 =             DN(1,0)*v(1,0);
    const double crhs8 =             DN(1,1)*v(1,1);
    const double crhs9 =             DN(2,0)*v(2,0);
    const double crhs10 =             DN(2,1)*v(2,1);
    const double crhs11 =             DN(3,0)*v(3,0);
    const double crhs12 =             DN(3,1)*v(3,1);
    const double crhs13 =             crhs10 + crhs11 + crhs12 + crhs4 + crhs5 + crhs6 + crhs7 + crhs8 + crhs9;
    const double crhs14 =             (N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]) + N[3]*(bdf0*p[3] + bdf1*pn[3] + bdf2*pnn[3]))/(pow(c, 2)*rho);
    const double crhs15 =             tau2*(crhs13 + crhs14);
    const double crhs16 =             v(0,0) - vmesh(0,0);
    const double crhs17 =             v(1,0) - vmesh(1,0);
    const double crhs18 =             v(2,0) - vmesh(2,0);
    const double crhs19 =             v(3,0) - vmesh(3,0);
    const double crhs20 =             N[0]*crhs16 + N[1]*crhs17 + N[2]*crhs18 + N[3]*crhs19;
    const double crhs21 =             v(0,1) - vmesh(0,1);
    const double crhs22 =             v(1,1) - vmesh(1,1);
    const double crhs23 =             v(2,1) - vmesh(2,1);
    const double crhs24 =             v(3,1) - vmesh(3,1);
    const double crhs25 =             N[0]*crhs21 + N[1]*crhs22 + N[2]*crhs23 + N[3]*crhs24;
    const double crhs26 =             v(0,2) - vmesh(0,2);
    const double crhs27 =             v(1,2) - vmesh(1,2);
    const double crhs28 =             v(2,2) - vmesh(2,2);
    const double crhs29 =             v(3,2) - vmesh(3,2);
    const double crhs30 =             N[0]*crhs26 + N[1]*crhs27 + N[2]*crhs28 + N[3]*crhs29;
    const double crhs31 =             crhs20*(crhs11 + crhs5 + crhs7 + crhs9) + crhs25*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs30*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0));
    const double crhs32 =             DN(0,0)*crhs16 + DN(0,1)*crhs21 + DN(0,2)*crhs26 + DN(1,0)*crhs17 + DN(1,1)*crhs22 + DN(1,2)*crhs27 + DN(2,0)*crhs18 + DN(2,1)*crhs23 + DN(2,2)*crhs28 + DN(3,0)*crhs19 + DN(3,1)*crhs24 + DN(3,2)*crhs29;
    const double crhs33 =             rho*(DN(0,0)*crhs20 + DN(0,1)*crhs25 + DN(0,2)*crhs30 + N[0]*crhs32);
    const double crhs34 =             tau1*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + rho*(crhs3 + crhs31));
    const double crhs35 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
    const double crhs36 =             N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1));
    const double crhs37 =             crhs20*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs25*(crhs10 + crhs12 + crhs6 + crhs8) + crhs30*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1));
    const double crhs38 =             tau1*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs35 + rho*(crhs36 + crhs37));
    const double crhs39 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
    const double crhs40 =             N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2));
    const double crhs41 =             crhs20*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs25*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + crhs30*crhs4;
    const double crhs42 =             tau1*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs39 + rho*(crhs40 + crhs41));
    const double crhs43 =             N[1]*rho;
    const double crhs44 =             rho*(DN(1,0)*crhs20 + DN(1,1)*crhs25 + DN(1,2)*crhs30 + N[1]*crhs32);
    const double crhs45 =             N[2]*rho;
    const double crhs46 =             rho*(DN(2,0)*crhs20 + DN(2,1)*crhs25 + DN(2,2)*crhs30 + N[2]*crhs32);
    const double crhs47 =             N[3]*rho;
    const double crhs48 =             rho*(DN(3,0)*crhs20 + DN(3,1)*crhs25 + DN(3,2)*crhs30 + N[3]*crhs32);

    rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs15 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - crhs2*crhs3 - crhs2*crhs31 - crhs33*crhs34;
    rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs15 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs35 - crhs2*crhs36 - crhs2*crhs37 - crhs33*crhs38;
    rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs15 - DN(0,2)*stress[2] + N[0]*crhs39 - crhs2*crhs40 - crhs2*crhs41 - crhs33*crhs42;
    rhs[3]=-DN(0,0)*crhs34 - DN(0,1)*crhs38 - DN(0,2)*crhs42 - N[0]*crhs13 - N[0]*crhs14;
    rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs15 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - crhs3*crhs43 - crhs31*crhs43 - crhs34*crhs44;
    rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs15 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs35 - crhs36*crhs43 - crhs37*crhs43 - crhs38*crhs44;
    rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs15 - DN(1,2)*stress[2] + N[1]*crhs39 - crhs40*crhs43 - crhs41*crhs43 - crhs42*crhs44;
    rhs[7]=-DN(1,0)*crhs34 - DN(1,1)*crhs38 - DN(1,2)*crhs42 - N[1]*crhs13 - N[1]*crhs14;
    rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs15 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - crhs3*crhs45 - crhs31*crhs45 - crhs34*crhs46;
    rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs15 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs35 - crhs36*crhs45 - crhs37*crhs45 - crhs38*crhs46;
    rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs15 - DN(2,2)*stress[2] + N[2]*crhs39 - crhs40*crhs45 - crhs41*crhs45 - crhs42*crhs46;
    rhs[11]=-DN(2,0)*crhs34 - DN(2,1)*crhs38 - DN(2,2)*crhs42 - N[2]*crhs13 - N[2]*crhs14;
    rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs15 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - crhs3*crhs47 - crhs31*crhs47 - crhs34*crhs48;
    rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs15 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs35 - crhs36*crhs47 - crhs37*crhs47 - crhs38*crhs48;
    rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs15 - DN(3,2)*stress[2] + N[3]*crhs39 - crhs40*crhs47 - crhs41*crhs47 - crhs42*crhs48;
    rhs[15]=-DN(3,0)*crhs34 - DN(3,1)*crhs38 - DN(3,2)*crhs42 - N[3]*crhs13 - N[3]*crhs14;

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
    const double crhs2 =             N[0]*rho;
    const double crhs3 =             N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0));
    const double crhs4 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
    const double crhs5 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
    const double crhs6 =             crhs4 + crhs5;
    const double crhs7 =             (N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]))/(pow(c, 2)*rho);
    const double crhs8 =             tau2*(crhs6 + crhs7);
    const double crhs9 =             v(0,0) - vmesh(0,0);
    const double crhs10 =             v(1,0) - vmesh(1,0);
    const double crhs11 =             v(2,0) - vmesh(2,0);
    const double crhs12 =             N[0]*crhs9 + N[1]*crhs10 + N[2]*crhs11;
    const double crhs13 =             v(0,1) - vmesh(0,1);
    const double crhs14 =             v(1,1) - vmesh(1,1);
    const double crhs15 =             v(2,1) - vmesh(2,1);
    const double crhs16 =             N[0]*crhs13 + N[1]*crhs14 + N[2]*crhs15;
    const double crhs17 =             crhs12*crhs4 + crhs16*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0));
    const double crhs18 =             DN(0,0)*crhs9 + DN(0,1)*crhs13 + DN(1,0)*crhs10 + DN(1,1)*crhs14 + DN(2,0)*crhs11 + DN(2,1)*crhs15;
    const double crhs19 =             rho*(DN(0,0)*crhs12 + DN(0,1)*crhs16 + N[0]*crhs18);
    const double crhs20 =             tau1*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + rho*(crhs17 + crhs3));
    const double crhs21 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
    const double crhs22 =             N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1));
    const double crhs23 =             crhs12*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + crhs16*crhs5;
    const double crhs24 =             tau1*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs21 + rho*(crhs22 + crhs23));
    const double crhs25 =             N[1]*rho;
    const double crhs26 =             rho*(DN(1,0)*crhs12 + DN(1,1)*crhs16 + N[1]*crhs18);
    const double crhs27 =             N[2]*rho;
    const double crhs28 =             rho*(DN(2,0)*crhs12 + DN(2,1)*crhs16 + N[2]*crhs18);

    rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs8 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - crhs17*crhs2 - crhs19*crhs20 - crhs2*crhs3;
    rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs8 - DN(0,1)*stress[1] + N[0]*crhs21 - crhs19*crhs24 - crhs2*crhs22 - crhs2*crhs23;
    rhs[2]=-DN(0,0)*crhs20 - DN(0,1)*crhs24 - N[0]*crhs6 - N[0]*crhs7;
    rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs8 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - crhs17*crhs25 - crhs20*crhs26 - crhs25*crhs3;
    rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs8 - DN(1,1)*stress[1] + N[1]*crhs21 - crhs22*crhs25 - crhs23*crhs25 - crhs24*crhs26;
    rhs[5]=-DN(1,0)*crhs20 - DN(1,1)*crhs24 - N[1]*crhs6 - N[1]*crhs7;
    rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs8 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - crhs17*crhs27 - crhs20*crhs28 - crhs27*crhs3;
    rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs8 - DN(2,1)*stress[1] + N[2]*crhs21 - crhs22*crhs27 - crhs23*crhs27 - crhs24*crhs28;
    rhs[8]=-DN(2,0)*crhs20 - DN(2,1)*crhs24 - N[2]*crhs6 - N[2]*crhs7;

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
    // const double c = data.c;                                // Wave velocity

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
    // const array_1d<double,nnodes>& pn = data.pn;
    // const array_1d<double,nnodes>& pnn = data.pnn;
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
    // const double c = data.c;                                // Wave velocity

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
    // const array_1d<double,nnodes>& pn = data.pn;
    // const array_1d<double,nnodes>& pnn = data.pnn;

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
