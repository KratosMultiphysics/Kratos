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
void NavierStokes<3>::ComputeGaussPointLHSContribution(bounded_matrix<double,16,16>& lhs, const element_data& data)
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
    //~ const Matrix& C = data.C;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    const array_1d<double,dim> vconv_gauss = prod(trans(vconv), N);

    const double vconv_norm = norm_2(vconv_gauss);

    // Stabilization parameters
    const double tau1 = 1.0/((rho*dyn_tau_coeff)/delta_t + (2*rho*vconv_norm)/h + (4*mu)/(h*h));
    const double tau2 = mu + 0.5*h*vconv_norm;

    const double clhs0 =             pow(DN(0,0), 2);
    const double clhs1 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
    const double clhs2 =             DN(0,0)*clhs1*rho;
    const double clhs3 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
    const double clhs4 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
    const double clhs5 =             DN(0,0)*clhs1 + DN(0,1)*clhs3 + DN(0,2)*clhs4;
    const double clhs6 =             pow(c, -2);
    const double clhs7 =             N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]) + N[3]*(bdf0*p[3] + bdf1*pn[3] + bdf2*pnn[3]);
    const double clhs8 =             clhs6*clhs7;
    const double clhs9 =             N[0]*clhs8;
    const double clhs10 =             clhs9 + rho*(N[0]*bdf0 + clhs5);
    const double clhs11 =             clhs10*tau1;
    const double clhs12 =             pow(N[0], 2);
    const double clhs13 =             bdf0*rho;
    const double clhs14 =             N[0]*rho;
    const double clhs15 =             N[0]*clhs6*clhs7*tau1;
    const double clhs16 =             -clhs10*clhs15 + clhs12*clhs13 + clhs12*clhs8 + clhs14*clhs5;
    const double clhs17 =             DN(0,0)*rho;
    const double clhs18 =             DN(0,0)*tau2;
    const double clhs19 =             clhs10*rho*tau1;
    const double clhs20 =             clhs1*clhs19 + clhs17 - clhs18;
    const double clhs21 =             DN(0,0)*N[0];
    const double clhs22 =             1.0/rho;
    const double clhs23 =             bdf0*clhs22*clhs6*tau2;
    const double clhs24 =             bdf0*clhs12*clhs6;
    const double clhs25 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0);
    const double clhs26 =             N[0]*bdf0*clhs6;
    const double clhs27 =             clhs25*clhs26;
    const double clhs28 =             DN(0,0) + clhs27;
    const double clhs29 =             DN(0,0)*clhs28*rho*tau1;
    const double clhs30 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
    const double clhs31 =             clhs26*clhs30;
    const double clhs32 =             DN(0,1) + clhs31;
    const double clhs33 =             DN(0,1)*clhs32*rho*tau1;
    const double clhs34 =             N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2);
    const double clhs35 =             clhs26*clhs34;
    const double clhs36 =             DN(0,2) + clhs35;
    const double clhs37 =             DN(0,2)*clhs36*rho*tau1;
    const double clhs38 =             bdf0*clhs12*clhs6*tau1;
    const double clhs39 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + clhs25*clhs8 - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)) + clhs1*(DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0)) + clhs3*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + clhs4*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)));
    const double clhs40 =             DN(1,0)*clhs17;
    const double clhs41 =             -DN(1,0)*clhs18;
    const double clhs42 =             N[0]*bdf0*rho;
    const double clhs43 =             N[1]*clhs42;
    const double clhs44 =             N[1]*clhs9;
    const double clhs45 =             DN(1,0)*clhs1 + DN(1,1)*clhs3 + DN(1,2)*clhs4;
    const double clhs46 =             N[1]*clhs8;
    const double clhs47 =             clhs46 + rho*(N[1]*bdf0 + clhs45);
    const double clhs48 =             clhs47*tau1;
    const double clhs49 =             clhs14*clhs45;
    const double clhs50 =             -clhs15*clhs47;
    const double clhs51 =             DN(1,1)*clhs17 - DN(1,1)*clhs18;
    const double clhs52 =             DN(0,1)*clhs47*rho*tau1;
    const double clhs53 =             DN(1,2)*clhs17 - DN(1,2)*clhs18;
    const double clhs54 =             DN(0,2)*clhs47*rho*tau1;
    const double clhs55 =             DN(0,0)*N[1];
    const double clhs56 =             N[1]*bdf0*clhs6;
    const double clhs57 =             clhs25*clhs56;
    const double clhs58 =             DN(1,0) + clhs57;
    const double clhs59 =             DN(0,0)*clhs58*rho*tau1;
    const double clhs60 =             clhs30*clhs56;
    const double clhs61 =             DN(1,1) + clhs60;
    const double clhs62 =             DN(0,1)*clhs61*rho*tau1;
    const double clhs63 =             clhs34*clhs56;
    const double clhs64 =             DN(1,2) + clhs63;
    const double clhs65 =             DN(0,2)*clhs64*rho*tau1;
    const double clhs66 =             N[0]*N[1]*bdf0*clhs6*tau1;
    const double clhs67 =             N[1]*clhs27 - clhs39*clhs66;
    const double clhs68 =             DN(2,0)*clhs17;
    const double clhs69 =             -DN(2,0)*clhs18;
    const double clhs70 =             N[2]*clhs42;
    const double clhs71 =             N[2]*clhs9;
    const double clhs72 =             DN(2,0)*clhs1 + DN(2,1)*clhs3 + DN(2,2)*clhs4;
    const double clhs73 =             N[2]*clhs8;
    const double clhs74 =             clhs73 + rho*(N[2]*bdf0 + clhs72);
    const double clhs75 =             clhs74*tau1;
    const double clhs76 =             clhs14*clhs72;
    const double clhs77 =             -clhs15*clhs74;
    const double clhs78 =             DN(2,1)*clhs17 - DN(2,1)*clhs18;
    const double clhs79 =             DN(0,1)*clhs74*rho*tau1;
    const double clhs80 =             DN(2,2)*clhs17 - DN(2,2)*clhs18;
    const double clhs81 =             DN(0,2)*clhs74*rho*tau1;
    const double clhs82 =             DN(0,0)*N[2];
    const double clhs83 =             N[2]*bdf0*clhs6;
    const double clhs84 =             DN(2,0) + clhs25*clhs83;
    const double clhs85 =             DN(0,0)*clhs84*rho*tau1;
    const double clhs86 =             DN(2,1) + clhs30*clhs83;
    const double clhs87 =             DN(0,1)*clhs86*rho*tau1;
    const double clhs88 =             DN(2,2) + clhs34*clhs83;
    const double clhs89 =             DN(0,2)*clhs88*rho*tau1;
    const double clhs90 =             N[0]*N[2]*bdf0*clhs6*tau1;
    const double clhs91 =             N[2]*clhs27 - clhs39*clhs90;
    const double clhs92 =             DN(3,0)*clhs17;
    const double clhs93 =             -DN(3,0)*clhs18;
    const double clhs94 =             N[3]*clhs42;
    const double clhs95 =             N[3]*clhs9;
    const double clhs96 =             DN(3,0)*clhs1 + DN(3,1)*clhs3 + DN(3,2)*clhs4;
    const double clhs97 =             N[3]*clhs8 + rho*(N[3]*bdf0 + clhs96);
    const double clhs98 =             clhs97*tau1;
    const double clhs99 =             clhs14*clhs96;
    const double clhs100 =             -clhs15*clhs97;
    const double clhs101 =             DN(3,1)*clhs17 - DN(3,1)*clhs18;
    const double clhs102 =             DN(0,1)*clhs97*rho*tau1;
    const double clhs103 =             DN(3,2)*clhs17 - DN(3,2)*clhs18;
    const double clhs104 =             DN(0,2)*clhs97*rho*tau1;
    const double clhs105 =             DN(0,0)*N[3];
    const double clhs106 =             N[3]*bdf0*clhs6;
    const double clhs107 =             DN(3,0) + clhs106*clhs25;
    const double clhs108 =             DN(0,0)*clhs107*rho*tau1;
    const double clhs109 =             DN(3,1) + clhs106*clhs30;
    const double clhs110 =             DN(0,1)*clhs109*rho*tau1;
    const double clhs111 =             DN(3,2) + clhs106*clhs34;
    const double clhs112 =             DN(0,2)*clhs111*rho*tau1;
    const double clhs113 =             N[0]*N[3]*bdf0*clhs6*tau1;
    const double clhs114 =             N[3]*clhs27 - clhs113*clhs39;
    const double clhs115 =             DN(0,1)*rho;
    const double clhs116 =             DN(0,1)*tau2;
    const double clhs117 =             clhs115 - clhs116 + clhs19*clhs3;
    const double clhs118 =             pow(DN(0,1), 2);
    const double clhs119 =             DN(0,1)*clhs3*rho;
    const double clhs120 =             DN(0,1)*N[0];
    const double clhs121 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + clhs30*clhs8 - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)) + clhs1*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + clhs3*(DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1)) + clhs4*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
    const double clhs122 =             DN(1,0)*clhs115 - DN(1,0)*clhs116;
    const double clhs123 =             DN(0,0)*clhs47*rho*tau1;
    const double clhs124 =             DN(1,1)*clhs115 - DN(1,1)*clhs116 + clhs43 + clhs44;
    const double clhs125 =             DN(1,2)*clhs115 - DN(1,2)*clhs116;
    const double clhs126 =             DN(0,1)*N[1];
    const double clhs127 =             N[1]*clhs31 - clhs121*clhs66;
    const double clhs128 =             DN(2,0)*clhs115 - DN(2,0)*clhs116;
    const double clhs129 =             DN(0,0)*clhs74*rho*tau1;
    const double clhs130 =             DN(2,1)*clhs115 - DN(2,1)*clhs116 + clhs70 + clhs71;
    const double clhs131 =             DN(2,2)*clhs115 - DN(2,2)*clhs116;
    const double clhs132 =             DN(0,1)*N[2];
    const double clhs133 =             N[2]*clhs31 - clhs121*clhs90;
    const double clhs134 =             DN(3,0)*clhs115 - DN(3,0)*clhs116;
    const double clhs135 =             DN(0,0)*clhs97*rho*tau1;
    const double clhs136 =             DN(3,1)*clhs115 - DN(3,1)*clhs116 + clhs94 + clhs95;
    const double clhs137 =             DN(3,2)*clhs115 - DN(3,2)*clhs116;
    const double clhs138 =             DN(0,1)*N[3];
    const double clhs139 =             N[3]*clhs31 - clhs113*clhs121;
    const double clhs140 =             DN(0,2)*rho;
    const double clhs141 =             DN(0,2)*tau2;
    const double clhs142 =             clhs140 - clhs141 + clhs19*clhs4;
    const double clhs143 =             pow(DN(0,2), 2);
    const double clhs144 =             DN(0,2)*clhs4*rho;
    const double clhs145 =             DN(0,2)*N[0];
    const double clhs146 =             DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + clhs34*clhs8 - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)) + clhs1*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + clhs3*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + clhs4*(DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2)));
    const double clhs147 =             DN(1,0)*clhs140 - DN(1,0)*clhs141;
    const double clhs148 =             DN(1,1)*clhs140 - DN(1,1)*clhs141;
    const double clhs149 =             DN(1,2)*clhs140 - DN(1,2)*clhs141 + clhs43 + clhs44;
    const double clhs150 =             DN(0,2)*N[1];
    const double clhs151 =             N[1]*clhs35 - clhs146*clhs66;
    const double clhs152 =             DN(2,0)*clhs140 - DN(2,0)*clhs141;
    const double clhs153 =             DN(2,1)*clhs140 - DN(2,1)*clhs141;
    const double clhs154 =             DN(2,2)*clhs140 - DN(2,2)*clhs141 + clhs70 + clhs71;
    const double clhs155 =             DN(0,2)*N[2];
    const double clhs156 =             N[2]*clhs35 - clhs146*clhs90;
    const double clhs157 =             DN(3,0)*clhs140 - DN(3,0)*clhs141;
    const double clhs158 =             DN(3,1)*clhs140 - DN(3,1)*clhs141;
    const double clhs159 =             DN(3,2)*clhs140 - DN(3,2)*clhs141 + clhs94 + clhs95;
    const double clhs160 =             DN(0,2)*N[3];
    const double clhs161 =             N[3]*clhs35 - clhs113*clhs146;
    const double clhs162 =             2*N[0];
    const double clhs163 =             clhs11 + clhs162;
    const double clhs164 =             bdf0*clhs22*clhs6;
    const double clhs165 =             DN(0,0)*tau1;
    const double clhs166 =             DN(0,1)*tau1;
    const double clhs167 =             DN(0,2)*tau1;
    const double clhs168 =             N[0]*bdf0*clhs22*clhs6;
    const double clhs169 =             N[1]*clhs168;
    const double clhs170 =             N[2]*clhs168;
    const double clhs171 =             N[3]*clhs168;
    const double clhs172 =             DN(1,0)*clhs10*rho*tau1;
    const double clhs173 =             N[1]*rho;
    const double clhs174 =             clhs173*clhs5;
    const double clhs175 =             N[1]*clhs6*clhs7*tau1;
    const double clhs176 =             -clhs10*clhs175;
    const double clhs177 =             DN(1,1)*clhs10*rho*tau1;
    const double clhs178 =             DN(1,2)*clhs10*rho*tau1;
    const double clhs179 =             DN(1,0)*N[0];
    const double clhs180 =             DN(1,0)*clhs28*rho*tau1;
    const double clhs181 =             DN(1,1)*clhs32*rho*tau1;
    const double clhs182 =             DN(1,2)*clhs36*rho*tau1;
    const double clhs183 =             pow(DN(1,0), 2);
    const double clhs184 =             DN(1,0)*clhs1*rho;
    const double clhs185 =             pow(N[1], 2);
    const double clhs186 =             clhs13*clhs185 + clhs173*clhs45 - clhs175*clhs47 + clhs185*clhs8;
    const double clhs187 =             DN(1,0)*rho;
    const double clhs188 =             DN(1,0)*tau2;
    const double clhs189 =             clhs47*rho*tau1;
    const double clhs190 =             clhs1*clhs189 + clhs187 - clhs188;
    const double clhs191 =             DN(1,0)*N[1];
    const double clhs192 =             bdf0*clhs185*clhs6;
    const double clhs193 =             DN(1,0)*clhs58*rho*tau1;
    const double clhs194 =             DN(1,1)*clhs61*rho*tau1;
    const double clhs195 =             DN(1,2)*clhs64*rho*tau1;
    const double clhs196 =             bdf0*clhs185*clhs6*tau1;
    const double clhs197 =             DN(2,0)*clhs187;
    const double clhs198 =             -DN(2,0)*clhs188;
    const double clhs199 =             N[1]*bdf0*rho;
    const double clhs200 =             N[2]*clhs199;
    const double clhs201 =             N[2]*clhs46;
    const double clhs202 =             clhs173*clhs72;
    const double clhs203 =             -clhs175*clhs74;
    const double clhs204 =             DN(2,1)*clhs187 - DN(2,1)*clhs188;
    const double clhs205 =             DN(1,1)*clhs74*rho*tau1;
    const double clhs206 =             DN(2,2)*clhs187 - DN(2,2)*clhs188;
    const double clhs207 =             DN(1,2)*clhs74*rho*tau1;
    const double clhs208 =             DN(1,0)*N[2];
    const double clhs209 =             DN(1,0)*clhs84*rho*tau1;
    const double clhs210 =             DN(1,1)*clhs86*rho*tau1;
    const double clhs211 =             DN(1,2)*clhs88*rho*tau1;
    const double clhs212 =             N[1]*N[2]*bdf0*clhs6*tau1;
    const double clhs213 =             N[2]*clhs57 - clhs212*clhs39;
    const double clhs214 =             DN(3,0)*clhs187;
    const double clhs215 =             -DN(3,0)*clhs188;
    const double clhs216 =             N[3]*clhs199;
    const double clhs217 =             N[3]*clhs46;
    const double clhs218 =             clhs173*clhs96;
    const double clhs219 =             -clhs175*clhs97;
    const double clhs220 =             DN(3,1)*clhs187 - DN(3,1)*clhs188;
    const double clhs221 =             DN(1,1)*clhs97*rho*tau1;
    const double clhs222 =             DN(3,2)*clhs187 - DN(3,2)*clhs188;
    const double clhs223 =             DN(1,2)*clhs97*rho*tau1;
    const double clhs224 =             DN(1,0)*N[3];
    const double clhs225 =             DN(1,0)*clhs107*rho*tau1;
    const double clhs226 =             DN(1,1)*clhs109*rho*tau1;
    const double clhs227 =             DN(1,2)*clhs111*rho*tau1;
    const double clhs228 =             N[1]*N[3]*bdf0*clhs6*tau1;
    const double clhs229 =             N[3]*clhs57 - clhs228*clhs39;
    const double clhs230 =             DN(1,1)*N[0];
    const double clhs231 =             DN(1,1)*rho;
    const double clhs232 =             DN(1,1)*tau2;
    const double clhs233 =             clhs189*clhs3 + clhs231 - clhs232;
    const double clhs234 =             pow(DN(1,1), 2);
    const double clhs235 =             DN(1,1)*clhs3*rho;
    const double clhs236 =             DN(1,1)*N[1];
    const double clhs237 =             DN(2,0)*clhs231 - DN(2,0)*clhs232;
    const double clhs238 =             DN(1,0)*clhs74*rho*tau1;
    const double clhs239 =             DN(2,1)*clhs231 - DN(2,1)*clhs232 + clhs200 + clhs201;
    const double clhs240 =             DN(2,2)*clhs231 - DN(2,2)*clhs232;
    const double clhs241 =             DN(1,1)*N[2];
    const double clhs242 =             N[2]*clhs60 - clhs121*clhs212;
    const double clhs243 =             DN(3,0)*clhs231 - DN(3,0)*clhs232;
    const double clhs244 =             DN(1,0)*clhs97*rho*tau1;
    const double clhs245 =             DN(3,1)*clhs231 - DN(3,1)*clhs232 + clhs216 + clhs217;
    const double clhs246 =             DN(3,2)*clhs231 - DN(3,2)*clhs232;
    const double clhs247 =             DN(1,1)*N[3];
    const double clhs248 =             N[3]*clhs60 - clhs121*clhs228;
    const double clhs249 =             DN(1,2)*N[0];
    const double clhs250 =             DN(1,2)*rho;
    const double clhs251 =             DN(1,2)*tau2;
    const double clhs252 =             clhs189*clhs4 + clhs250 - clhs251;
    const double clhs253 =             pow(DN(1,2), 2);
    const double clhs254 =             DN(1,2)*clhs4*rho;
    const double clhs255 =             DN(1,2)*N[1];
    const double clhs256 =             DN(2,0)*clhs250 - DN(2,0)*clhs251;
    const double clhs257 =             DN(2,1)*clhs250 - DN(2,1)*clhs251;
    const double clhs258 =             DN(2,2)*clhs250 - DN(2,2)*clhs251 + clhs200 + clhs201;
    const double clhs259 =             DN(1,2)*N[2];
    const double clhs260 =             N[2]*clhs63 - clhs146*clhs212;
    const double clhs261 =             DN(3,0)*clhs250 - DN(3,0)*clhs251;
    const double clhs262 =             DN(3,1)*clhs250 - DN(3,1)*clhs251;
    const double clhs263 =             DN(3,2)*clhs250 - DN(3,2)*clhs251 + clhs216 + clhs217;
    const double clhs264 =             DN(1,2)*N[3];
    const double clhs265 =             N[3]*clhs63 - clhs146*clhs228;
    const double clhs266 =             2*N[1];
    const double clhs267 =             DN(1,0)*tau1;
    const double clhs268 =             DN(1,1)*tau1;
    const double clhs269 =             DN(1,2)*tau1;
    const double clhs270 =             clhs266 + clhs48;
    const double clhs271 =             N[1]*bdf0*clhs22*clhs6;
    const double clhs272 =             N[2]*clhs271;
    const double clhs273 =             N[3]*clhs271;
    const double clhs274 =             DN(2,0)*clhs10*rho*tau1;
    const double clhs275 =             N[2]*rho;
    const double clhs276 =             clhs275*clhs5;
    const double clhs277 =             N[2]*clhs6*clhs7*tau1;
    const double clhs278 =             -clhs10*clhs277;
    const double clhs279 =             DN(2,1)*clhs10*rho*tau1;
    const double clhs280 =             DN(2,2)*clhs10*rho*tau1;
    const double clhs281 =             DN(2,0)*N[0];
    const double clhs282 =             DN(2,0)*clhs28*rho*tau1;
    const double clhs283 =             DN(2,1)*clhs32*rho*tau1;
    const double clhs284 =             DN(2,2)*clhs36*rho*tau1;
    const double clhs285 =             DN(2,0)*clhs47*rho*tau1;
    const double clhs286 =             clhs275*clhs45;
    const double clhs287 =             -clhs277*clhs47;
    const double clhs288 =             DN(2,1)*clhs47*rho*tau1;
    const double clhs289 =             DN(2,2)*clhs47*rho*tau1;
    const double clhs290 =             DN(2,0)*N[1];
    const double clhs291 =             DN(2,0)*clhs58*rho*tau1;
    const double clhs292 =             DN(2,1)*clhs61*rho*tau1;
    const double clhs293 =             DN(2,2)*clhs64*rho*tau1;
    const double clhs294 =             pow(DN(2,0), 2);
    const double clhs295 =             DN(2,0)*clhs1*rho;
    const double clhs296 =             pow(N[2], 2);
    const double clhs297 =             clhs13*clhs296 + clhs275*clhs72 - clhs277*clhs74 + clhs296*clhs8;
    const double clhs298 =             DN(2,0)*rho;
    const double clhs299 =             DN(2,0)*tau2;
    const double clhs300 =             clhs74*rho*tau1;
    const double clhs301 =             clhs1*clhs300 + clhs298 - clhs299;
    const double clhs302 =             DN(2,0)*N[2];
    const double clhs303 =             bdf0*clhs296*clhs6;
    const double clhs304 =             DN(2,0)*clhs84*rho*tau1;
    const double clhs305 =             DN(2,1)*clhs86*rho*tau1;
    const double clhs306 =             DN(2,2)*clhs88*rho*tau1;
    const double clhs307 =             bdf0*clhs296*clhs6*tau1;
    const double clhs308 =             DN(3,0)*clhs298;
    const double clhs309 =             -DN(3,0)*clhs299;
    const double clhs310 =             N[2]*N[3]*bdf0;
    const double clhs311 =             clhs310*rho;
    const double clhs312 =             N[3]*clhs73;
    const double clhs313 =             clhs275*clhs96;
    const double clhs314 =             -clhs277*clhs97;
    const double clhs315 =             DN(3,1)*clhs298 - DN(3,1)*clhs299;
    const double clhs316 =             DN(2,1)*clhs97*rho*tau1;
    const double clhs317 =             DN(3,2)*clhs298 - DN(3,2)*clhs299;
    const double clhs318 =             DN(2,2)*clhs97*rho*tau1;
    const double clhs319 =             DN(2,0)*N[3];
    const double clhs320 =             DN(2,0)*clhs107*rho*tau1;
    const double clhs321 =             DN(2,1)*clhs109*rho*tau1;
    const double clhs322 =             DN(2,2)*clhs111*rho*tau1;
    const double clhs323 =             N[2]*N[3]*bdf0*clhs6;
    const double clhs324 =             N[2]*N[3]*bdf0*clhs6*tau1;
    const double clhs325 =             clhs25*clhs323 - clhs324*clhs39;
    const double clhs326 =             DN(2,1)*N[0];
    const double clhs327 =             DN(2,1)*N[1];
    const double clhs328 =             DN(2,1)*rho;
    const double clhs329 =             DN(2,1)*tau2;
    const double clhs330 =             clhs3*clhs300 + clhs328 - clhs329;
    const double clhs331 =             pow(DN(2,1), 2);
    const double clhs332 =             DN(2,1)*clhs3*rho;
    const double clhs333 =             DN(2,1)*N[2];
    const double clhs334 =             DN(3,0)*clhs328 - DN(3,0)*clhs329;
    const double clhs335 =             DN(2,0)*clhs97*rho*tau1;
    const double clhs336 =             DN(3,1)*clhs328 - DN(3,1)*clhs329 + clhs311 + clhs312;
    const double clhs337 =             DN(3,2)*clhs328 - DN(3,2)*clhs329;
    const double clhs338 =             DN(2,1)*N[3];
    const double clhs339 =             -clhs121*clhs324 + clhs30*clhs323;
    const double clhs340 =             DN(2,2)*N[0];
    const double clhs341 =             DN(2,2)*N[1];
    const double clhs342 =             DN(2,2)*rho;
    const double clhs343 =             DN(2,2)*tau2;
    const double clhs344 =             clhs300*clhs4 + clhs342 - clhs343;
    const double clhs345 =             pow(DN(2,2), 2);
    const double clhs346 =             DN(2,2)*clhs4*rho;
    const double clhs347 =             DN(2,2)*N[2];
    const double clhs348 =             DN(3,0)*clhs342 - DN(3,0)*clhs343;
    const double clhs349 =             DN(3,1)*clhs342 - DN(3,1)*clhs343;
    const double clhs350 =             DN(3,2)*clhs342 - DN(3,2)*clhs343 + clhs311 + clhs312;
    const double clhs351 =             DN(2,2)*N[3];
    const double clhs352 =             -clhs146*clhs324 + clhs323*clhs34;
    const double clhs353 =             2*N[2];
    const double clhs354 =             DN(2,0)*tau1;
    const double clhs355 =             DN(2,1)*tau1;
    const double clhs356 =             DN(2,2)*tau1;
    const double clhs357 =             clhs353 + clhs75;
    const double clhs358 =             clhs22*clhs310*clhs6;
    const double clhs359 =             DN(3,0)*clhs10*rho*tau1;
    const double clhs360 =             N[3]*rho;
    const double clhs361 =             clhs360*clhs5;
    const double clhs362 =             N[3]*clhs6*clhs7*tau1;
    const double clhs363 =             -clhs10*clhs362;
    const double clhs364 =             DN(3,1)*clhs10*rho*tau1;
    const double clhs365 =             DN(3,2)*clhs10*rho*tau1;
    const double clhs366 =             DN(3,0)*N[0];
    const double clhs367 =             DN(3,0)*clhs28*rho*tau1;
    const double clhs368 =             DN(3,1)*clhs32*rho*tau1;
    const double clhs369 =             DN(3,2)*clhs36*rho*tau1;
    const double clhs370 =             DN(3,0)*clhs47*rho*tau1;
    const double clhs371 =             clhs360*clhs45;
    const double clhs372 =             -clhs362*clhs47;
    const double clhs373 =             DN(3,1)*clhs47*rho*tau1;
    const double clhs374 =             DN(3,2)*clhs47*rho*tau1;
    const double clhs375 =             DN(3,0)*N[1];
    const double clhs376 =             DN(3,0)*clhs58*rho*tau1;
    const double clhs377 =             DN(3,1)*clhs61*rho*tau1;
    const double clhs378 =             DN(3,2)*clhs64*rho*tau1;
    const double clhs379 =             DN(3,0)*clhs74*rho*tau1;
    const double clhs380 =             clhs360*clhs72;
    const double clhs381 =             -clhs362*clhs74;
    const double clhs382 =             DN(3,1)*clhs74*rho*tau1;
    const double clhs383 =             DN(3,2)*clhs74*rho*tau1;
    const double clhs384 =             DN(3,0)*N[2];
    const double clhs385 =             DN(3,0)*clhs84*rho*tau1;
    const double clhs386 =             DN(3,1)*clhs86*rho*tau1;
    const double clhs387 =             DN(3,2)*clhs88*rho*tau1;
    const double clhs388 =             pow(DN(3,0), 2);
    const double clhs389 =             pow(N[3], 2);
    const double clhs390 =             clhs13*clhs389 + clhs360*clhs96 - clhs362*clhs97 + clhs389*clhs8;
    const double clhs391 =             clhs97*rho*tau1;
    const double clhs392 =             DN(3,0)*rho - DN(3,0)*tau2 + clhs1*clhs391;
    const double clhs393 =             DN(3,0)*N[3];
    const double clhs394 =             bdf0*clhs389*clhs6;
    const double clhs395 =             DN(3,0)*clhs107*rho*tau1;
    const double clhs396 =             DN(3,1)*clhs109*rho*tau1;
    const double clhs397 =             DN(3,2)*clhs111*rho*tau1;
    const double clhs398 =             bdf0*clhs389*clhs6*tau1;
    const double clhs399 =             DN(3,1)*N[0];
    const double clhs400 =             DN(3,1)*N[1];
    const double clhs401 =             DN(3,1)*N[2];
    const double clhs402 =             DN(3,1)*rho - DN(3,1)*tau2 + clhs3*clhs391;
    const double clhs403 =             pow(DN(3,1), 2);
    const double clhs404 =             DN(3,1)*N[3];
    const double clhs405 =             DN(3,2)*N[0];
    const double clhs406 =             DN(3,2)*N[1];
    const double clhs407 =             DN(3,2)*N[2];
    const double clhs408 =             DN(3,2)*rho - DN(3,2)*tau2 + clhs391*clhs4;
    const double clhs409 =             pow(DN(3,2), 2);
    const double clhs410 =             DN(3,2)*N[3];
    const double clhs411 =             2*N[3];
    const double clhs412 =             DN(3,0)*tau1;
    const double clhs413 =             DN(3,1)*tau1;
    const double clhs414 =             DN(3,2)*tau1;
    const double clhs415 =             clhs411 + clhs98;

    lhs(0,0)=clhs0*rho - clhs0*tau2 + clhs11*clhs2 + clhs16;
    lhs(0,1)=DN(0,1)*clhs20;
    lhs(0,2)=DN(0,2)*clhs20;
    lhs(0,3)=clhs1*clhs29 + clhs1*clhs33 + clhs1*clhs37 - clhs15*clhs28 - clhs21*clhs23 - clhs21 + clhs24*clhs25 - clhs38*clhs39;
    lhs(0,4)=clhs2*clhs48 + clhs40 + clhs41 + clhs43 + clhs44 + clhs49 + clhs50;
    lhs(0,5)=clhs1*clhs52 + clhs51;
    lhs(0,6)=clhs1*clhs54 + clhs53;
    lhs(0,7)=clhs1*clhs59 + clhs1*clhs62 + clhs1*clhs65 - clhs15*clhs58 - clhs23*clhs55 - clhs55 + clhs67;
    lhs(0,8)=clhs2*clhs75 + clhs68 + clhs69 + clhs70 + clhs71 + clhs76 + clhs77;
    lhs(0,9)=clhs1*clhs79 + clhs78;
    lhs(0,10)=clhs1*clhs81 + clhs80;
    lhs(0,11)=clhs1*clhs85 + clhs1*clhs87 + clhs1*clhs89 - clhs15*clhs84 - clhs23*clhs82 - clhs82 + clhs91;
    lhs(0,12)=clhs100 + clhs2*clhs98 + clhs92 + clhs93 + clhs94 + clhs95 + clhs99;
    lhs(0,13)=clhs1*clhs102 + clhs101;
    lhs(0,14)=clhs1*clhs104 + clhs103;
    lhs(0,15)=clhs1*clhs108 + clhs1*clhs110 + clhs1*clhs112 - clhs105*clhs23 - clhs105 - clhs107*clhs15 + clhs114;
    lhs(1,0)=DN(0,0)*clhs117;
    lhs(1,1)=clhs11*clhs119 + clhs118*rho - clhs118*tau2 + clhs16;
    lhs(1,2)=DN(0,2)*clhs117;
    lhs(1,3)=-clhs120*clhs23 - clhs120 - clhs121*clhs38 - clhs15*clhs32 + clhs24*clhs30 + clhs29*clhs3 + clhs3*clhs33 + clhs3*clhs37;
    lhs(1,4)=clhs122 + clhs123*clhs3;
    lhs(1,5)=clhs119*clhs48 + clhs124 + clhs49 + clhs50;
    lhs(1,6)=clhs125 + clhs3*clhs54;
    lhs(1,7)=-clhs126*clhs23 - clhs126 + clhs127 - clhs15*clhs61 + clhs3*clhs59 + clhs3*clhs62 + clhs3*clhs65;
    lhs(1,8)=clhs128 + clhs129*clhs3;
    lhs(1,9)=clhs119*clhs75 + clhs130 + clhs76 + clhs77;
    lhs(1,10)=clhs131 + clhs3*clhs81;
    lhs(1,11)=-clhs132*clhs23 - clhs132 + clhs133 - clhs15*clhs86 + clhs3*clhs85 + clhs3*clhs87 + clhs3*clhs89;
    lhs(1,12)=clhs134 + clhs135*clhs3;
    lhs(1,13)=clhs100 + clhs119*clhs98 + clhs136 + clhs99;
    lhs(1,14)=clhs104*clhs3 + clhs137;
    lhs(1,15)=clhs108*clhs3 - clhs109*clhs15 + clhs110*clhs3 + clhs112*clhs3 - clhs138*clhs23 - clhs138 + clhs139;
    lhs(2,0)=DN(0,0)*clhs142;
    lhs(2,1)=DN(0,1)*clhs142;
    lhs(2,2)=clhs11*clhs144 + clhs143*rho - clhs143*tau2 + clhs16;
    lhs(2,3)=-clhs145*clhs23 - clhs145 - clhs146*clhs38 - clhs15*clhs36 + clhs24*clhs34 + clhs29*clhs4 + clhs33*clhs4 + clhs37*clhs4;
    lhs(2,4)=clhs123*clhs4 + clhs147;
    lhs(2,5)=clhs148 + clhs4*clhs52;
    lhs(2,6)=clhs144*clhs48 + clhs149 + clhs49 + clhs50;
    lhs(2,7)=-clhs15*clhs64 - clhs150*clhs23 - clhs150 + clhs151 + clhs4*clhs59 + clhs4*clhs62 + clhs4*clhs65;
    lhs(2,8)=clhs129*clhs4 + clhs152;
    lhs(2,9)=clhs153 + clhs4*clhs79;
    lhs(2,10)=clhs144*clhs75 + clhs154 + clhs76 + clhs77;
    lhs(2,11)=-clhs15*clhs88 - clhs155*clhs23 - clhs155 + clhs156 + clhs4*clhs85 + clhs4*clhs87 + clhs4*clhs89;
    lhs(2,12)=clhs135*clhs4 + clhs157;
    lhs(2,13)=clhs102*clhs4 + clhs158;
    lhs(2,14)=clhs100 + clhs144*clhs98 + clhs159 + clhs99;
    lhs(2,15)=clhs108*clhs4 + clhs110*clhs4 - clhs111*clhs15 + clhs112*clhs4 - clhs160*clhs23 - clhs160 + clhs161;
    lhs(3,0)=DN(0,0)*clhs163;
    lhs(3,1)=DN(0,1)*clhs163;
    lhs(3,2)=DN(0,2)*clhs163;
    lhs(3,3)=clhs12*clhs164 + clhs165*clhs28 + clhs166*clhs32 + clhs167*clhs36;
    lhs(3,4)=DN(1,0)*clhs162 + clhs165*clhs47;
    lhs(3,5)=DN(1,1)*clhs162 + clhs166*clhs47;
    lhs(3,6)=DN(1,2)*clhs162 + clhs167*clhs47;
    lhs(3,7)=clhs165*clhs58 + clhs166*clhs61 + clhs167*clhs64 + clhs169;
    lhs(3,8)=DN(2,0)*clhs162 + clhs165*clhs74;
    lhs(3,9)=DN(2,1)*clhs162 + clhs166*clhs74;
    lhs(3,10)=DN(2,2)*clhs162 + clhs167*clhs74;
    lhs(3,11)=clhs165*clhs84 + clhs166*clhs86 + clhs167*clhs88 + clhs170;
    lhs(3,12)=DN(3,0)*clhs162 + clhs165*clhs97;
    lhs(3,13)=DN(3,1)*clhs162 + clhs166*clhs97;
    lhs(3,14)=DN(3,2)*clhs162 + clhs167*clhs97;
    lhs(3,15)=clhs107*clhs165 + clhs109*clhs166 + clhs111*clhs167 + clhs171;
    lhs(4,0)=clhs1*clhs172 + clhs174 + clhs176 + clhs40 + clhs41 + clhs43 + clhs44;
    lhs(4,1)=clhs1*clhs177 + clhs122;
    lhs(4,2)=clhs1*clhs178 + clhs147;
    lhs(4,3)=clhs1*clhs180 + clhs1*clhs181 + clhs1*clhs182 - clhs175*clhs28 - clhs179*clhs23 - clhs179 + clhs67;
    lhs(4,4)=clhs183*rho - clhs183*tau2 + clhs184*clhs48 + clhs186;
    lhs(4,5)=DN(1,1)*clhs190;
    lhs(4,6)=DN(1,2)*clhs190;
    lhs(4,7)=clhs1*clhs193 + clhs1*clhs194 + clhs1*clhs195 - clhs175*clhs58 - clhs191*clhs23 - clhs191 + clhs192*clhs25 - clhs196*clhs39;
    lhs(4,8)=clhs184*clhs75 + clhs197 + clhs198 + clhs200 + clhs201 + clhs202 + clhs203;
    lhs(4,9)=clhs1*clhs205 + clhs204;
    lhs(4,10)=clhs1*clhs207 + clhs206;
    lhs(4,11)=clhs1*clhs209 + clhs1*clhs210 + clhs1*clhs211 - clhs175*clhs84 - clhs208*clhs23 - clhs208 + clhs213;
    lhs(4,12)=clhs184*clhs98 + clhs214 + clhs215 + clhs216 + clhs217 + clhs218 + clhs219;
    lhs(4,13)=clhs1*clhs221 + clhs220;
    lhs(4,14)=clhs1*clhs223 + clhs222;
    lhs(4,15)=clhs1*clhs225 + clhs1*clhs226 + clhs1*clhs227 - clhs107*clhs175 - clhs224*clhs23 - clhs224 + clhs229;
    lhs(5,0)=clhs172*clhs3 + clhs51;
    lhs(5,1)=clhs124 + clhs174 + clhs176 + clhs177*clhs3;
    lhs(5,2)=clhs148 + clhs178*clhs3;
    lhs(5,3)=clhs127 - clhs175*clhs32 + clhs180*clhs3 + clhs181*clhs3 + clhs182*clhs3 - clhs23*clhs230 - clhs230;
    lhs(5,4)=DN(1,0)*clhs233;
    lhs(5,5)=clhs186 + clhs234*rho - clhs234*tau2 + clhs235*clhs48;
    lhs(5,6)=DN(1,2)*clhs233;
    lhs(5,7)=-clhs121*clhs196 - clhs175*clhs61 + clhs192*clhs30 + clhs193*clhs3 + clhs194*clhs3 + clhs195*clhs3 - clhs23*clhs236 - clhs236;
    lhs(5,8)=clhs237 + clhs238*clhs3;
    lhs(5,9)=clhs202 + clhs203 + clhs235*clhs75 + clhs239;
    lhs(5,10)=clhs207*clhs3 + clhs240;
    lhs(5,11)=-clhs175*clhs86 + clhs209*clhs3 + clhs210*clhs3 + clhs211*clhs3 - clhs23*clhs241 - clhs241 + clhs242;
    lhs(5,12)=clhs243 + clhs244*clhs3;
    lhs(5,13)=clhs218 + clhs219 + clhs235*clhs98 + clhs245;
    lhs(5,14)=clhs223*clhs3 + clhs246;
    lhs(5,15)=-clhs109*clhs175 + clhs225*clhs3 + clhs226*clhs3 + clhs227*clhs3 - clhs23*clhs247 - clhs247 + clhs248;
    lhs(6,0)=clhs172*clhs4 + clhs53;
    lhs(6,1)=clhs125 + clhs177*clhs4;
    lhs(6,2)=clhs149 + clhs174 + clhs176 + clhs178*clhs4;
    lhs(6,3)=clhs151 - clhs175*clhs36 + clhs180*clhs4 + clhs181*clhs4 + clhs182*clhs4 - clhs23*clhs249 - clhs249;
    lhs(6,4)=DN(1,0)*clhs252;
    lhs(6,5)=DN(1,1)*clhs252;
    lhs(6,6)=clhs186 + clhs253*rho - clhs253*tau2 + clhs254*clhs48;
    lhs(6,7)=-clhs146*clhs196 - clhs175*clhs64 + clhs192*clhs34 + clhs193*clhs4 + clhs194*clhs4 + clhs195*clhs4 - clhs23*clhs255 - clhs255;
    lhs(6,8)=clhs238*clhs4 + clhs256;
    lhs(6,9)=clhs205*clhs4 + clhs257;
    lhs(6,10)=clhs202 + clhs203 + clhs254*clhs75 + clhs258;
    lhs(6,11)=-clhs175*clhs88 + clhs209*clhs4 + clhs210*clhs4 + clhs211*clhs4 - clhs23*clhs259 - clhs259 + clhs260;
    lhs(6,12)=clhs244*clhs4 + clhs261;
    lhs(6,13)=clhs221*clhs4 + clhs262;
    lhs(6,14)=clhs218 + clhs219 + clhs254*clhs98 + clhs263;
    lhs(6,15)=-clhs111*clhs175 + clhs225*clhs4 + clhs226*clhs4 + clhs227*clhs4 - clhs23*clhs264 - clhs264 + clhs265;
    lhs(7,0)=DN(0,0)*clhs266 + clhs10*clhs267;
    lhs(7,1)=DN(0,1)*clhs266 + clhs10*clhs268;
    lhs(7,2)=DN(0,2)*clhs266 + clhs10*clhs269;
    lhs(7,3)=clhs169 + clhs267*clhs28 + clhs268*clhs32 + clhs269*clhs36;
    lhs(7,4)=DN(1,0)*clhs270;
    lhs(7,5)=DN(1,1)*clhs270;
    lhs(7,6)=DN(1,2)*clhs270;
    lhs(7,7)=clhs164*clhs185 + clhs267*clhs58 + clhs268*clhs61 + clhs269*clhs64;
    lhs(7,8)=DN(2,0)*clhs266 + clhs267*clhs74;
    lhs(7,9)=DN(2,1)*clhs266 + clhs268*clhs74;
    lhs(7,10)=DN(2,2)*clhs266 + clhs269*clhs74;
    lhs(7,11)=clhs267*clhs84 + clhs268*clhs86 + clhs269*clhs88 + clhs272;
    lhs(7,12)=DN(3,0)*clhs266 + clhs267*clhs97;
    lhs(7,13)=DN(3,1)*clhs266 + clhs268*clhs97;
    lhs(7,14)=DN(3,2)*clhs266 + clhs269*clhs97;
    lhs(7,15)=clhs107*clhs267 + clhs109*clhs268 + clhs111*clhs269 + clhs273;
    lhs(8,0)=clhs1*clhs274 + clhs276 + clhs278 + clhs68 + clhs69 + clhs70 + clhs71;
    lhs(8,1)=clhs1*clhs279 + clhs128;
    lhs(8,2)=clhs1*clhs280 + clhs152;
    lhs(8,3)=clhs1*clhs282 + clhs1*clhs283 + clhs1*clhs284 - clhs23*clhs281 - clhs277*clhs28 - clhs281 + clhs91;
    lhs(8,4)=clhs1*clhs285 + clhs197 + clhs198 + clhs200 + clhs201 + clhs286 + clhs287;
    lhs(8,5)=clhs1*clhs288 + clhs237;
    lhs(8,6)=clhs1*clhs289 + clhs256;
    lhs(8,7)=clhs1*clhs291 + clhs1*clhs292 + clhs1*clhs293 + clhs213 - clhs23*clhs290 - clhs277*clhs58 - clhs290;
    lhs(8,8)=clhs294*rho - clhs294*tau2 + clhs295*clhs75 + clhs297;
    lhs(8,9)=DN(2,1)*clhs301;
    lhs(8,10)=DN(2,2)*clhs301;
    lhs(8,11)=clhs1*clhs304 + clhs1*clhs305 + clhs1*clhs306 - clhs23*clhs302 + clhs25*clhs303 - clhs277*clhs84 - clhs302 - clhs307*clhs39;
    lhs(8,12)=clhs295*clhs98 + clhs308 + clhs309 + clhs311 + clhs312 + clhs313 + clhs314;
    lhs(8,13)=clhs1*clhs316 + clhs315;
    lhs(8,14)=clhs1*clhs318 + clhs317;
    lhs(8,15)=clhs1*clhs320 + clhs1*clhs321 + clhs1*clhs322 - clhs107*clhs277 - clhs23*clhs319 - clhs319 + clhs325;
    lhs(9,0)=clhs274*clhs3 + clhs78;
    lhs(9,1)=clhs130 + clhs276 + clhs278 + clhs279*clhs3;
    lhs(9,2)=clhs153 + clhs280*clhs3;
    lhs(9,3)=clhs133 - clhs23*clhs326 - clhs277*clhs32 + clhs282*clhs3 + clhs283*clhs3 + clhs284*clhs3 - clhs326;
    lhs(9,4)=clhs204 + clhs285*clhs3;
    lhs(9,5)=clhs239 + clhs286 + clhs287 + clhs288*clhs3;
    lhs(9,6)=clhs257 + clhs289*clhs3;
    lhs(9,7)=-clhs23*clhs327 + clhs242 - clhs277*clhs61 + clhs291*clhs3 + clhs292*clhs3 + clhs293*clhs3 - clhs327;
    lhs(9,8)=DN(2,0)*clhs330;
    lhs(9,9)=clhs297 + clhs331*rho - clhs331*tau2 + clhs332*clhs75;
    lhs(9,10)=DN(2,2)*clhs330;
    lhs(9,11)=-clhs121*clhs307 - clhs23*clhs333 - clhs277*clhs86 + clhs3*clhs304 + clhs3*clhs305 + clhs3*clhs306 + clhs30*clhs303 - clhs333;
    lhs(9,12)=clhs3*clhs335 + clhs334;
    lhs(9,13)=clhs313 + clhs314 + clhs332*clhs98 + clhs336;
    lhs(9,14)=clhs3*clhs318 + clhs337;
    lhs(9,15)=-clhs109*clhs277 - clhs23*clhs338 + clhs3*clhs320 + clhs3*clhs321 + clhs3*clhs322 - clhs338 + clhs339;
    lhs(10,0)=clhs274*clhs4 + clhs80;
    lhs(10,1)=clhs131 + clhs279*clhs4;
    lhs(10,2)=clhs154 + clhs276 + clhs278 + clhs280*clhs4;
    lhs(10,3)=clhs156 - clhs23*clhs340 - clhs277*clhs36 + clhs282*clhs4 + clhs283*clhs4 + clhs284*clhs4 - clhs340;
    lhs(10,4)=clhs206 + clhs285*clhs4;
    lhs(10,5)=clhs240 + clhs288*clhs4;
    lhs(10,6)=clhs258 + clhs286 + clhs287 + clhs289*clhs4;
    lhs(10,7)=-clhs23*clhs341 + clhs260 - clhs277*clhs64 + clhs291*clhs4 + clhs292*clhs4 + clhs293*clhs4 - clhs341;
    lhs(10,8)=DN(2,0)*clhs344;
    lhs(10,9)=DN(2,1)*clhs344;
    lhs(10,10)=clhs297 + clhs345*rho - clhs345*tau2 + clhs346*clhs75;
    lhs(10,11)=-clhs146*clhs307 - clhs23*clhs347 - clhs277*clhs88 + clhs303*clhs34 + clhs304*clhs4 + clhs305*clhs4 + clhs306*clhs4 - clhs347;
    lhs(10,12)=clhs335*clhs4 + clhs348;
    lhs(10,13)=clhs316*clhs4 + clhs349;
    lhs(10,14)=clhs313 + clhs314 + clhs346*clhs98 + clhs350;
    lhs(10,15)=-clhs111*clhs277 - clhs23*clhs351 + clhs320*clhs4 + clhs321*clhs4 + clhs322*clhs4 - clhs351 + clhs352;
    lhs(11,0)=DN(0,0)*clhs353 + clhs10*clhs354;
    lhs(11,1)=DN(0,1)*clhs353 + clhs10*clhs355;
    lhs(11,2)=DN(0,2)*clhs353 + clhs10*clhs356;
    lhs(11,3)=clhs170 + clhs28*clhs354 + clhs32*clhs355 + clhs356*clhs36;
    lhs(11,4)=DN(1,0)*clhs353 + clhs354*clhs47;
    lhs(11,5)=DN(1,1)*clhs353 + clhs355*clhs47;
    lhs(11,6)=DN(1,2)*clhs353 + clhs356*clhs47;
    lhs(11,7)=clhs272 + clhs354*clhs58 + clhs355*clhs61 + clhs356*clhs64;
    lhs(11,8)=DN(2,0)*clhs357;
    lhs(11,9)=DN(2,1)*clhs357;
    lhs(11,10)=DN(2,2)*clhs357;
    lhs(11,11)=clhs164*clhs296 + clhs354*clhs84 + clhs355*clhs86 + clhs356*clhs88;
    lhs(11,12)=DN(3,0)*clhs353 + clhs354*clhs97;
    lhs(11,13)=DN(3,1)*clhs353 + clhs355*clhs97;
    lhs(11,14)=DN(3,2)*clhs353 + clhs356*clhs97;
    lhs(11,15)=clhs107*clhs354 + clhs109*clhs355 + clhs111*clhs356 + clhs358;
    lhs(12,0)=clhs1*clhs359 + clhs361 + clhs363 + clhs92 + clhs93 + clhs94 + clhs95;
    lhs(12,1)=clhs1*clhs364 + clhs134;
    lhs(12,2)=clhs1*clhs365 + clhs157;
    lhs(12,3)=clhs1*clhs367 + clhs1*clhs368 + clhs1*clhs369 + clhs114 - clhs23*clhs366 - clhs28*clhs362 - clhs366;
    lhs(12,4)=clhs1*clhs370 + clhs214 + clhs215 + clhs216 + clhs217 + clhs371 + clhs372;
    lhs(12,5)=clhs1*clhs373 + clhs243;
    lhs(12,6)=clhs1*clhs374 + clhs261;
    lhs(12,7)=clhs1*clhs376 + clhs1*clhs377 + clhs1*clhs378 + clhs229 - clhs23*clhs375 - clhs362*clhs58 - clhs375;
    lhs(12,8)=clhs1*clhs379 + clhs308 + clhs309 + clhs311 + clhs312 + clhs380 + clhs381;
    lhs(12,9)=clhs1*clhs382 + clhs334;
    lhs(12,10)=clhs1*clhs383 + clhs348;
    lhs(12,11)=clhs1*clhs385 + clhs1*clhs386 + clhs1*clhs387 - clhs23*clhs384 + clhs325 - clhs362*clhs84 - clhs384;
    lhs(12,12)=DN(3,0)*clhs1*clhs98*rho + clhs388*rho - clhs388*tau2 + clhs390;
    lhs(12,13)=DN(3,1)*clhs392;
    lhs(12,14)=DN(3,2)*clhs392;
    lhs(12,15)=clhs1*clhs395 + clhs1*clhs396 + clhs1*clhs397 - clhs107*clhs362 - clhs23*clhs393 + clhs25*clhs394 - clhs39*clhs398 - clhs393;
    lhs(13,0)=clhs101 + clhs3*clhs359;
    lhs(13,1)=clhs136 + clhs3*clhs364 + clhs361 + clhs363;
    lhs(13,2)=clhs158 + clhs3*clhs365;
    lhs(13,3)=clhs139 - clhs23*clhs399 + clhs3*clhs367 + clhs3*clhs368 + clhs3*clhs369 - clhs32*clhs362 - clhs399;
    lhs(13,4)=clhs220 + clhs3*clhs370;
    lhs(13,5)=clhs245 + clhs3*clhs373 + clhs371 + clhs372;
    lhs(13,6)=clhs262 + clhs3*clhs374;
    lhs(13,7)=-clhs23*clhs400 + clhs248 + clhs3*clhs376 + clhs3*clhs377 + clhs3*clhs378 - clhs362*clhs61 - clhs400;
    lhs(13,8)=clhs3*clhs379 + clhs315;
    lhs(13,9)=clhs3*clhs382 + clhs336 + clhs380 + clhs381;
    lhs(13,10)=clhs3*clhs383 + clhs349;
    lhs(13,11)=-clhs23*clhs401 + clhs3*clhs385 + clhs3*clhs386 + clhs3*clhs387 + clhs339 - clhs362*clhs86 - clhs401;
    lhs(13,12)=DN(3,0)*clhs402;
    lhs(13,13)=DN(3,1)*clhs3*clhs98*rho + clhs390 + clhs403*rho - clhs403*tau2;
    lhs(13,14)=DN(3,2)*clhs402;
    lhs(13,15)=-clhs109*clhs362 - clhs121*clhs398 - clhs23*clhs404 + clhs3*clhs395 + clhs3*clhs396 + clhs3*clhs397 + clhs30*clhs394 - clhs404;
    lhs(14,0)=clhs103 + clhs359*clhs4;
    lhs(14,1)=clhs137 + clhs364*clhs4;
    lhs(14,2)=clhs159 + clhs361 + clhs363 + clhs365*clhs4;
    lhs(14,3)=clhs161 - clhs23*clhs405 - clhs36*clhs362 + clhs367*clhs4 + clhs368*clhs4 + clhs369*clhs4 - clhs405;
    lhs(14,4)=clhs222 + clhs370*clhs4;
    lhs(14,5)=clhs246 + clhs373*clhs4;
    lhs(14,6)=clhs263 + clhs371 + clhs372 + clhs374*clhs4;
    lhs(14,7)=-clhs23*clhs406 + clhs265 - clhs362*clhs64 + clhs376*clhs4 + clhs377*clhs4 + clhs378*clhs4 - clhs406;
    lhs(14,8)=clhs317 + clhs379*clhs4;
    lhs(14,9)=clhs337 + clhs382*clhs4;
    lhs(14,10)=clhs350 + clhs380 + clhs381 + clhs383*clhs4;
    lhs(14,11)=-clhs23*clhs407 + clhs352 - clhs362*clhs88 + clhs385*clhs4 + clhs386*clhs4 + clhs387*clhs4 - clhs407;
    lhs(14,12)=DN(3,0)*clhs408;
    lhs(14,13)=DN(3,1)*clhs408;
    lhs(14,14)=DN(3,2)*clhs4*clhs98*rho + clhs390 + clhs409*rho - clhs409*tau2;
    lhs(14,15)=-clhs111*clhs362 - clhs146*clhs398 - clhs23*clhs410 + clhs34*clhs394 + clhs395*clhs4 + clhs396*clhs4 + clhs397*clhs4 - clhs410;
    lhs(15,0)=DN(0,0)*clhs411 + clhs10*clhs412;
    lhs(15,1)=DN(0,1)*clhs411 + clhs10*clhs413;
    lhs(15,2)=DN(0,2)*clhs411 + clhs10*clhs414;
    lhs(15,3)=clhs171 + clhs28*clhs412 + clhs32*clhs413 + clhs36*clhs414;
    lhs(15,4)=DN(1,0)*clhs411 + clhs412*clhs47;
    lhs(15,5)=DN(1,1)*clhs411 + clhs413*clhs47;
    lhs(15,6)=DN(1,2)*clhs411 + clhs414*clhs47;
    lhs(15,7)=clhs273 + clhs412*clhs58 + clhs413*clhs61 + clhs414*clhs64;
    lhs(15,8)=DN(2,0)*clhs411 + clhs412*clhs74;
    lhs(15,9)=DN(2,1)*clhs411 + clhs413*clhs74;
    lhs(15,10)=DN(2,2)*clhs411 + clhs414*clhs74;
    lhs(15,11)=clhs358 + clhs412*clhs84 + clhs413*clhs86 + clhs414*clhs88;
    lhs(15,12)=DN(3,0)*clhs415;
    lhs(15,13)=DN(3,1)*clhs415;
    lhs(15,14)=DN(3,2)*clhs415;
    lhs(15,15)=clhs107*clhs412 + clhs109*clhs413 + clhs111*clhs414 + clhs164*clhs389;

}


template<>
void NavierStokes<2>::ComputeGaussPointLHSContribution(bounded_matrix<double,9,9>& lhs, const element_data& data)
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
    //~ const Matrix& C = data.C;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    const array_1d<double,dim> vconv_gauss = prod(trans(vconv), N);

    const double vconv_norm = norm_2(vconv_gauss);

    // Stabilization parameters
    const double tau1 = 1.0/((rho*dyn_tau_coeff)/delta_t + (2*rho*vconv_norm)/h + (4*mu)/(h*h));
    const double tau2 = mu + 0.5*h*vconv_norm;

    const double clhs0 =             pow(DN(0,0), 2);
    const double clhs1 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
    const double clhs2 =             DN(0,0)*clhs1*rho;
    const double clhs3 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
    const double clhs4 =             DN(0,0)*clhs1 + DN(0,1)*clhs3;
    const double clhs5 =             pow(c, -2);
    const double clhs6 =             N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]);
    const double clhs7 =             clhs5*clhs6;
    const double clhs8 =             N[0]*clhs7;
    const double clhs9 =             clhs8 + rho*(N[0]*bdf0 + clhs4);
    const double clhs10 =             clhs9*tau1;
    const double clhs11 =             pow(N[0], 2);
    const double clhs12 =             bdf0*rho;
    const double clhs13 =             N[0]*rho;
    const double clhs14 =             N[0]*clhs5*clhs6*tau1;
    const double clhs15 =             clhs11*clhs12 + clhs11*clhs7 + clhs13*clhs4 - clhs14*clhs9;
    const double clhs16 =             DN(0,0)*rho;
    const double clhs17 =             DN(0,0)*tau2;
    const double clhs18 =             clhs9*rho*tau1;
    const double clhs19 =             DN(0,0)*N[0];
    const double clhs20 =             1.0/rho;
    const double clhs21 =             bdf0*clhs20*clhs5*tau2;
    const double clhs22 =             bdf0*clhs11*clhs5;
    const double clhs23 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0);
    const double clhs24 =             N[0]*bdf0*clhs5;
    const double clhs25 =             clhs23*clhs24;
    const double clhs26 =             DN(0,0) + clhs25;
    const double clhs27 =             DN(0,0)*clhs26*rho*tau1;
    const double clhs28 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
    const double clhs29 =             clhs24*clhs28;
    const double clhs30 =             DN(0,1) + clhs29;
    const double clhs31 =             DN(0,1)*clhs30*rho*tau1;
    const double clhs32 =             bdf0*clhs11*clhs5*tau1;
    const double clhs33 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + clhs23*clhs7 - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + clhs1*(DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0)) + clhs3*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
    const double clhs34 =             DN(1,0)*clhs16;
    const double clhs35 =             -DN(1,0)*clhs17;
    const double clhs36 =             N[0]*bdf0*rho;
    const double clhs37 =             N[1]*clhs36;
    const double clhs38 =             N[1]*clhs8;
    const double clhs39 =             DN(1,0)*clhs1 + DN(1,1)*clhs3;
    const double clhs40 =             N[1]*clhs7;
    const double clhs41 =             clhs40 + rho*(N[1]*bdf0 + clhs39);
    const double clhs42 =             clhs41*tau1;
    const double clhs43 =             clhs13*clhs39;
    const double clhs44 =             -clhs14*clhs41;
    const double clhs45 =             DN(1,1)*clhs16 - DN(1,1)*clhs17;
    const double clhs46 =             DN(0,1)*clhs1*rho;
    const double clhs47 =             DN(0,0)*N[1];
    const double clhs48 =             N[1]*bdf0*clhs5;
    const double clhs49 =             DN(1,0) + clhs23*clhs48;
    const double clhs50 =             DN(0,0)*clhs49*rho*tau1;
    const double clhs51 =             DN(1,1) + clhs28*clhs48;
    const double clhs52 =             DN(0,1)*clhs51*rho*tau1;
    const double clhs53 =             N[0]*N[1]*bdf0*clhs5*tau1;
    const double clhs54 =             N[1]*clhs25 - clhs33*clhs53;
    const double clhs55 =             DN(2,0)*clhs16;
    const double clhs56 =             -DN(2,0)*clhs17;
    const double clhs57 =             N[2]*clhs36;
    const double clhs58 =             N[2]*clhs8;
    const double clhs59 =             DN(2,0)*clhs1 + DN(2,1)*clhs3;
    const double clhs60 =             N[2]*clhs7 + rho*(N[2]*bdf0 + clhs59);
    const double clhs61 =             clhs60*tau1;
    const double clhs62 =             clhs13*clhs59;
    const double clhs63 =             -clhs14*clhs60;
    const double clhs64 =             DN(2,1)*clhs16 - DN(2,1)*clhs17;
    const double clhs65 =             DN(0,0)*N[2];
    const double clhs66 =             N[2]*bdf0*clhs5;
    const double clhs67 =             DN(2,0) + clhs23*clhs66;
    const double clhs68 =             DN(0,0)*clhs67*rho*tau1;
    const double clhs69 =             DN(2,1) + clhs28*clhs66;
    const double clhs70 =             DN(0,1)*clhs69*rho*tau1;
    const double clhs71 =             N[0]*N[2]*bdf0*clhs5*tau1;
    const double clhs72 =             N[2]*clhs25 - clhs33*clhs71;
    const double clhs73 =             DN(0,1)*rho;
    const double clhs74 =             DN(0,1)*tau2;
    const double clhs75 =             pow(DN(0,1), 2);
    const double clhs76 =             DN(0,1)*clhs3*rho;
    const double clhs77 =             DN(0,1)*N[0];
    const double clhs78 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + clhs28*clhs7 - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + clhs1*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + clhs3*(DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1)));
    const double clhs79 =             DN(1,0)*clhs73 - DN(1,0)*clhs74;
    const double clhs80 =             DN(0,0)*clhs3*rho;
    const double clhs81 =             DN(1,1)*clhs73 - DN(1,1)*clhs74 + clhs37 + clhs38;
    const double clhs82 =             DN(0,1)*N[1];
    const double clhs83 =             N[1]*clhs29 - clhs53*clhs78;
    const double clhs84 =             DN(2,0)*clhs73 - DN(2,0)*clhs74;
    const double clhs85 =             DN(2,1)*clhs73 - DN(2,1)*clhs74 + clhs57 + clhs58;
    const double clhs86 =             DN(0,1)*N[2];
    const double clhs87 =             N[2]*clhs29 - clhs71*clhs78;
    const double clhs88 =             2*N[0];
    const double clhs89 =             clhs10 + clhs88;
    const double clhs90 =             bdf0*clhs20*clhs5;
    const double clhs91 =             DN(0,0)*tau1;
    const double clhs92 =             DN(0,1)*tau1;
    const double clhs93 =             N[0]*bdf0*clhs20*clhs5;
    const double clhs94 =             N[1]*clhs93;
    const double clhs95 =             N[2]*clhs93;
    const double clhs96 =             DN(1,0)*clhs9*rho*tau1;
    const double clhs97 =             N[1]*rho;
    const double clhs98 =             clhs4*clhs97;
    const double clhs99 =             N[1]*clhs5*clhs6*tau1;
    const double clhs100 =             -clhs9*clhs99;
    const double clhs101 =             DN(1,1)*clhs9*rho*tau1;
    const double clhs102 =             DN(1,0)*N[0];
    const double clhs103 =             DN(1,0)*clhs26*rho*tau1;
    const double clhs104 =             DN(1,1)*clhs30*rho*tau1;
    const double clhs105 =             pow(DN(1,0), 2);
    const double clhs106 =             DN(1,0)*clhs1*rho;
    const double clhs107 =             pow(N[1], 2);
    const double clhs108 =             clhs107*clhs12 + clhs107*clhs7 + clhs39*clhs97 - clhs41*clhs99;
    const double clhs109 =             DN(1,0)*rho;
    const double clhs110 =             DN(1,0)*tau2;
    const double clhs111 =             clhs41*rho*tau1;
    const double clhs112 =             DN(1,0)*N[1];
    const double clhs113 =             bdf0*clhs107*clhs5;
    const double clhs114 =             DN(1,0)*clhs49*rho*tau1;
    const double clhs115 =             DN(1,1)*clhs51*rho*tau1;
    const double clhs116 =             bdf0*clhs107*clhs5*tau1;
    const double clhs117 =             DN(2,0)*clhs109;
    const double clhs118 =             -DN(2,0)*clhs110;
    const double clhs119 =             N[1]*N[2]*bdf0;
    const double clhs120 =             clhs119*rho;
    const double clhs121 =             N[2]*clhs40;
    const double clhs122 =             clhs59*clhs97;
    const double clhs123 =             -clhs60*clhs99;
    const double clhs124 =             DN(2,1)*clhs109 - DN(2,1)*clhs110;
    const double clhs125 =             DN(1,1)*rho;
    const double clhs126 =             clhs1*clhs60*tau1;
    const double clhs127 =             DN(1,0)*N[2];
    const double clhs128 =             DN(1,0)*clhs67*rho*tau1;
    const double clhs129 =             DN(1,1)*clhs69*rho*tau1;
    const double clhs130 =             N[1]*N[2]*bdf0*clhs5;
    const double clhs131 =             N[1]*N[2]*bdf0*clhs5*tau1;
    const double clhs132 =             clhs130*clhs23 - clhs131*clhs33;
    const double clhs133 =             DN(1,1)*N[0];
    const double clhs134 =             DN(1,1)*tau2;
    const double clhs135 =             pow(DN(1,1), 2);
    const double clhs136 =             clhs3*clhs41*tau1;
    const double clhs137 =             DN(1,1)*N[1];
    const double clhs138 =             DN(2,0)*clhs125 - DN(2,0)*clhs134;
    const double clhs139 =             clhs3*clhs60*tau1;
    const double clhs140 =             DN(2,1)*clhs125 - DN(2,1)*clhs134 + clhs120 + clhs121;
    const double clhs141 =             DN(1,1)*N[2];
    const double clhs142 =             clhs130*clhs28 - clhs131*clhs78;
    const double clhs143 =             2*N[1];
    const double clhs144 =             DN(1,0)*tau1;
    const double clhs145 =             DN(1,1)*tau1;
    const double clhs146 =             clhs143 + clhs42;
    const double clhs147 =             clhs119*clhs20*clhs5;
    const double clhs148 =             DN(2,0)*clhs9*rho*tau1;
    const double clhs149 =             N[2]*rho;
    const double clhs150 =             clhs149*clhs4;
    const double clhs151 =             N[2]*clhs5*clhs6*tau1;
    const double clhs152 =             -clhs151*clhs9;
    const double clhs153 =             DN(2,1)*clhs9*rho*tau1;
    const double clhs154 =             DN(2,0)*N[0];
    const double clhs155 =             DN(2,0)*clhs26*rho*tau1;
    const double clhs156 =             DN(2,1)*clhs30*rho*tau1;
    const double clhs157 =             DN(2,0)*rho;
    const double clhs158 =             clhs1*clhs41*tau1;
    const double clhs159 =             clhs149*clhs39;
    const double clhs160 =             -clhs151*clhs41;
    const double clhs161 =             DN(2,1)*rho;
    const double clhs162 =             DN(2,0)*N[1];
    const double clhs163 =             DN(2,0)*clhs49*rho*tau1;
    const double clhs164 =             DN(2,1)*clhs51*rho*tau1;
    const double clhs165 =             pow(DN(2,0), 2);
    const double clhs166 =             pow(N[2], 2);
    const double clhs167 =             clhs12*clhs166 + clhs149*clhs59 - clhs151*clhs60 + clhs166*clhs7;
    const double clhs168 =             clhs60*rho*tau1;
    const double clhs169 =             DN(2,0)*N[2];
    const double clhs170 =             bdf0*clhs166*clhs5;
    const double clhs171 =             DN(2,0)*clhs67*rho*tau1;
    const double clhs172 =             DN(2,1)*clhs69*rho*tau1;
    const double clhs173 =             bdf0*clhs166*clhs5*tau1;
    const double clhs174 =             DN(2,1)*N[0];
    const double clhs175 =             DN(2,1)*N[1];
    const double clhs176 =             pow(DN(2,1), 2);
    const double clhs177 =             DN(2,1)*N[2];
    const double clhs178 =             2*N[2];
    const double clhs179 =             DN(2,0)*tau1;
    const double clhs180 =             DN(2,1)*tau1;
    const double clhs181 =             clhs178 + clhs61;

    lhs(0,0)=clhs0*rho - clhs0*tau2 + clhs10*clhs2 + clhs15;
    lhs(0,1)=DN(0,1)*(clhs1*clhs18 + clhs16 - clhs17);
    lhs(0,2)=clhs1*clhs27 + clhs1*clhs31 - clhs14*clhs26 - clhs19*clhs21 - clhs19 + clhs22*clhs23 - clhs32*clhs33;
    lhs(0,3)=clhs2*clhs42 + clhs34 + clhs35 + clhs37 + clhs38 + clhs43 + clhs44;
    lhs(0,4)=clhs42*clhs46 + clhs45;
    lhs(0,5)=clhs1*clhs50 + clhs1*clhs52 - clhs14*clhs49 - clhs21*clhs47 - clhs47 + clhs54;
    lhs(0,6)=clhs2*clhs61 + clhs55 + clhs56 + clhs57 + clhs58 + clhs62 + clhs63;
    lhs(0,7)=clhs46*clhs61 + clhs64;
    lhs(0,8)=clhs1*clhs68 + clhs1*clhs70 - clhs14*clhs67 - clhs21*clhs65 - clhs65 + clhs72;
    lhs(1,0)=DN(0,0)*(clhs18*clhs3 + clhs73 - clhs74);
    lhs(1,1)=clhs10*clhs76 + clhs15 + clhs75*rho - clhs75*tau2;
    lhs(1,2)=-clhs14*clhs30 - clhs21*clhs77 + clhs22*clhs28 + clhs27*clhs3 + clhs3*clhs31 - clhs32*clhs78 - clhs77;
    lhs(1,3)=clhs42*clhs80 + clhs79;
    lhs(1,4)=clhs42*clhs76 + clhs43 + clhs44 + clhs81;
    lhs(1,5)=-clhs14*clhs51 - clhs21*clhs82 + clhs3*clhs50 + clhs3*clhs52 - clhs82 + clhs83;
    lhs(1,6)=clhs61*clhs80 + clhs84;
    lhs(1,7)=clhs61*clhs76 + clhs62 + clhs63 + clhs85;
    lhs(1,8)=-clhs14*clhs69 - clhs21*clhs86 + clhs3*clhs68 + clhs3*clhs70 - clhs86 + clhs87;
    lhs(2,0)=DN(0,0)*clhs89;
    lhs(2,1)=DN(0,1)*clhs89;
    lhs(2,2)=clhs11*clhs90 + clhs26*clhs91 + clhs30*clhs92;
    lhs(2,3)=DN(1,0)*clhs88 + clhs41*clhs91;
    lhs(2,4)=DN(1,1)*clhs88 + clhs41*clhs92;
    lhs(2,5)=clhs49*clhs91 + clhs51*clhs92 + clhs94;
    lhs(2,6)=DN(2,0)*clhs88 + clhs60*clhs91;
    lhs(2,7)=DN(2,1)*clhs88 + clhs60*clhs92;
    lhs(2,8)=clhs67*clhs91 + clhs69*clhs92 + clhs95;
    lhs(3,0)=clhs1*clhs96 + clhs100 + clhs34 + clhs35 + clhs37 + clhs38 + clhs98;
    lhs(3,1)=clhs1*clhs101 + clhs79;
    lhs(3,2)=clhs1*clhs103 + clhs1*clhs104 - clhs102*clhs21 - clhs102 - clhs26*clhs99 + clhs54;
    lhs(3,3)=clhs105*rho - clhs105*tau2 + clhs106*clhs42 + clhs108;
    lhs(3,4)=DN(1,1)*(clhs1*clhs111 + clhs109 - clhs110);
    lhs(3,5)=clhs1*clhs114 + clhs1*clhs115 - clhs112*clhs21 - clhs112 + clhs113*clhs23 - clhs116*clhs33 - clhs49*clhs99;
    lhs(3,6)=clhs106*clhs61 + clhs117 + clhs118 + clhs120 + clhs121 + clhs122 + clhs123;
    lhs(3,7)=clhs124 + clhs125*clhs126;
    lhs(3,8)=clhs1*clhs128 + clhs1*clhs129 - clhs127*clhs21 - clhs127 + clhs132 - clhs67*clhs99;
    lhs(4,0)=clhs3*clhs96 + clhs45;
    lhs(4,1)=clhs100 + clhs101*clhs3 + clhs81 + clhs98;
    lhs(4,2)=clhs103*clhs3 + clhs104*clhs3 - clhs133*clhs21 - clhs133 - clhs30*clhs99 + clhs83;
    lhs(4,3)=DN(1,0)*(clhs111*clhs3 + clhs125 - clhs134);
    lhs(4,4)=clhs108 + clhs125*clhs136 + clhs135*rho - clhs135*tau2;
    lhs(4,5)=clhs113*clhs28 + clhs114*clhs3 + clhs115*clhs3 - clhs116*clhs78 - clhs137*clhs21 - clhs137 - clhs51*clhs99;
    lhs(4,6)=clhs109*clhs139 + clhs138;
    lhs(4,7)=clhs122 + clhs123 + clhs125*clhs139 + clhs140;
    lhs(4,8)=clhs128*clhs3 + clhs129*clhs3 - clhs141*clhs21 - clhs141 + clhs142 - clhs69*clhs99;
    lhs(5,0)=DN(0,0)*clhs143 + clhs144*clhs9;
    lhs(5,1)=DN(0,1)*clhs143 + clhs145*clhs9;
    lhs(5,2)=clhs144*clhs26 + clhs145*clhs30 + clhs94;
    lhs(5,3)=DN(1,0)*clhs146;
    lhs(5,4)=DN(1,1)*clhs146;
    lhs(5,5)=clhs107*clhs90 + clhs144*clhs49 + clhs145*clhs51;
    lhs(5,6)=DN(2,0)*clhs143 + clhs144*clhs60;
    lhs(5,7)=DN(2,1)*clhs143 + clhs145*clhs60;
    lhs(5,8)=clhs144*clhs67 + clhs145*clhs69 + clhs147;
    lhs(6,0)=clhs1*clhs148 + clhs150 + clhs152 + clhs55 + clhs56 + clhs57 + clhs58;
    lhs(6,1)=clhs1*clhs153 + clhs84;
    lhs(6,2)=clhs1*clhs155 + clhs1*clhs156 - clhs151*clhs26 - clhs154*clhs21 - clhs154 + clhs72;
    lhs(6,3)=clhs117 + clhs118 + clhs120 + clhs121 + clhs157*clhs158 + clhs159 + clhs160;
    lhs(6,4)=clhs138 + clhs158*clhs161;
    lhs(6,5)=clhs1*clhs163 + clhs1*clhs164 + clhs132 - clhs151*clhs49 - clhs162*clhs21 - clhs162;
    lhs(6,6)=clhs126*clhs157 + clhs165*rho - clhs165*tau2 + clhs167;
    lhs(6,7)=DN(2,1)*(-DN(2,0)*tau2 + clhs1*clhs168 + clhs157);
    lhs(6,8)=clhs1*clhs171 + clhs1*clhs172 - clhs151*clhs67 - clhs169*clhs21 - clhs169 + clhs170*clhs23 - clhs173*clhs33;
    lhs(7,0)=clhs148*clhs3 + clhs64;
    lhs(7,1)=clhs150 + clhs152 + clhs153*clhs3 + clhs85;
    lhs(7,2)=-clhs151*clhs30 + clhs155*clhs3 + clhs156*clhs3 - clhs174*clhs21 - clhs174 + clhs87;
    lhs(7,3)=clhs124 + clhs136*clhs157;
    lhs(7,4)=clhs136*clhs161 + clhs140 + clhs159 + clhs160;
    lhs(7,5)=clhs142 - clhs151*clhs51 + clhs163*clhs3 + clhs164*clhs3 - clhs175*clhs21 - clhs175;
    lhs(7,6)=DN(2,0)*(-DN(2,1)*tau2 + clhs161 + clhs168*clhs3);
    lhs(7,7)=clhs139*clhs161 + clhs167 + clhs176*rho - clhs176*tau2;
    lhs(7,8)=-clhs151*clhs69 + clhs170*clhs28 + clhs171*clhs3 + clhs172*clhs3 - clhs173*clhs78 - clhs177*clhs21 - clhs177;
    lhs(8,0)=DN(0,0)*clhs178 + clhs179*clhs9;
    lhs(8,1)=DN(0,1)*clhs178 + clhs180*clhs9;
    lhs(8,2)=clhs179*clhs26 + clhs180*clhs30 + clhs95;
    lhs(8,3)=DN(1,0)*clhs178 + clhs179*clhs41;
    lhs(8,4)=DN(1,1)*clhs178 + clhs180*clhs41;
    lhs(8,5)=clhs147 + clhs179*clhs49 + clhs180*clhs51;
    lhs(8,6)=DN(2,0)*clhs181;
    lhs(8,7)=DN(2,1)*clhs181;
    lhs(8,8)=clhs166*clhs90 + clhs179*clhs67 + clhs180*clhs69;

}


template<>
void NavierStokes<3>::ComputeGaussPointRHSContribution(array_1d<double,16>& rhs, const element_data& data)
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
    const double tau1 = 1.0/((rho*dyn_tau_coeff)/delta_t + (2*rho*vconv_norm)/h + (4*mu)/(h*h));
    const double tau2 = mu + 0.5*h*vconv_norm;

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
    const double crhs12 =             crhs11*rho;
    const double crhs13 =             N[0]*rho;
    const double crhs14 =             N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0));
    const double crhs15 =             pow(c, -2);
    const double crhs16 =             N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]) + N[3]*(bdf0*p[3] + bdf1*pn[3] + bdf2*pnn[3]);
    const double crhs17 =             crhs15*crhs16;
    const double crhs18 =             crhs17*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0));
    const double crhs19 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
    const double crhs20 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
    const double crhs21 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
    const double crhs22 =             crhs19*(crhs3 + crhs5 + crhs7 + crhs9) + crhs20*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs21*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0));
    const double crhs23 =             crhs15*crhs16/rho;
    const double crhs24 =             tau2*(crhs11 + crhs23);
    const double crhs25 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + crhs18 + rho*(crhs14 + crhs22);
    const double crhs26 =             DN(0,0)*crhs25*rho*tau1;
    const double crhs27 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
    const double crhs28 =             crhs17*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1));
    const double crhs29 =             N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1));
    const double crhs30 =             crhs19*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs20*(crhs10 + crhs4 + crhs6 + crhs8) + crhs21*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1));
    const double crhs31 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs27 + crhs28 + rho*(crhs29 + crhs30);
    const double crhs32 =             DN(0,1)*crhs31*rho*tau1;
    const double crhs33 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
    const double crhs34 =             crhs17*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2));
    const double crhs35 =             N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2));
    const double crhs36 =             crhs19*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs2*crhs21 + crhs20*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2));
    const double crhs37 =             DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs33 + crhs34 + rho*(crhs35 + crhs36);
    const double crhs38 =             DN(0,2)*crhs37*rho*tau1;
    const double crhs39 =             N[0]*crhs15*crhs16;
    const double crhs40 =             crhs25*tau1;
    const double crhs41 =             crhs31*tau1;
    const double crhs42 =             crhs37*tau1;
    const double crhs43 =             2*DN(0,0)*v(0,0) + 2*DN(0,1)*v(0,1) + 2*DN(0,2)*v(0,2) + 2*DN(1,0)*v(1,0) + 2*DN(1,1)*v(1,1) + 2*DN(1,2)*v(1,2) + 2*DN(2,0)*v(2,0) + 2*DN(2,1)*v(2,1) + 2*DN(2,2)*v(2,2) + 2*DN(3,0)*v(3,0) + 2*DN(3,1)*v(3,1) + 2*DN(3,2)*v(3,2);
    const double crhs44 =             N[1]*rho;
    const double crhs45 =             DN(1,0)*crhs25*rho*tau1;
    const double crhs46 =             DN(1,1)*crhs31*rho*tau1;
    const double crhs47 =             DN(1,2)*crhs37*rho*tau1;
    const double crhs48 =             N[1]*crhs15*crhs16;
    const double crhs49 =             N[2]*rho;
    const double crhs50 =             DN(2,0)*crhs25*rho*tau1;
    const double crhs51 =             DN(2,1)*crhs31*rho*tau1;
    const double crhs52 =             DN(2,2)*crhs37*rho*tau1;
    const double crhs53 =             N[2]*crhs15*crhs16;
    const double crhs54 =             N[3]*rho;
    const double crhs55 =             DN(3,0)*crhs25*rho*tau1;
    const double crhs56 =             DN(3,1)*crhs31*rho*tau1;
    const double crhs57 =             DN(3,2)*crhs37*rho*tau1;
    const double crhs58 =             N[3]*crhs15*crhs16;

    rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs12 + DN(0,0)*crhs24 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - N[0]*crhs18 - crhs13*crhs14 - crhs13*crhs22 - crhs19*crhs26 - crhs19*crhs32 - crhs19*crhs38 + crhs39*crhs40;
    rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs12 + DN(0,1)*crhs24 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs27 - N[0]*crhs28 - crhs13*crhs29 - crhs13*crhs30 - crhs20*crhs26 - crhs20*crhs32 - crhs20*crhs38 + crhs39*crhs41;
    rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs12 + DN(0,2)*crhs24 - DN(0,2)*stress[2] + N[0]*crhs33 - N[0]*crhs34 - crhs13*crhs35 - crhs13*crhs36 - crhs21*crhs26 - crhs21*crhs32 - crhs21*crhs38 + crhs39*crhs42;
    rhs[3]=-DN(0,0)*crhs40 - DN(0,1)*crhs41 - DN(0,2)*crhs42 - N[0]*crhs23 - N[0]*crhs43;
    rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs12 + DN(1,0)*crhs24 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - N[1]*crhs18 - crhs14*crhs44 - crhs19*crhs45 - crhs19*crhs46 - crhs19*crhs47 - crhs22*crhs44 + crhs40*crhs48;
    rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs12 + DN(1,1)*crhs24 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs27 - N[1]*crhs28 - crhs20*crhs45 - crhs20*crhs46 - crhs20*crhs47 - crhs29*crhs44 - crhs30*crhs44 + crhs41*crhs48;
    rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs12 + DN(1,2)*crhs24 - DN(1,2)*stress[2] + N[1]*crhs33 - N[1]*crhs34 - crhs21*crhs45 - crhs21*crhs46 - crhs21*crhs47 - crhs35*crhs44 - crhs36*crhs44 + crhs42*crhs48;
    rhs[7]=-DN(1,0)*crhs40 - DN(1,1)*crhs41 - DN(1,2)*crhs42 - N[1]*crhs23 - N[1]*crhs43;
    rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs12 + DN(2,0)*crhs24 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - N[2]*crhs18 - crhs14*crhs49 - crhs19*crhs50 - crhs19*crhs51 - crhs19*crhs52 - crhs22*crhs49 + crhs40*crhs53;
    rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs12 + DN(2,1)*crhs24 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs27 - N[2]*crhs28 - crhs20*crhs50 - crhs20*crhs51 - crhs20*crhs52 - crhs29*crhs49 - crhs30*crhs49 + crhs41*crhs53;
    rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs12 + DN(2,2)*crhs24 - DN(2,2)*stress[2] + N[2]*crhs33 - N[2]*crhs34 - crhs21*crhs50 - crhs21*crhs51 - crhs21*crhs52 - crhs35*crhs49 - crhs36*crhs49 + crhs42*crhs53;
    rhs[11]=-DN(2,0)*crhs40 - DN(2,1)*crhs41 - DN(2,2)*crhs42 - N[2]*crhs23 - N[2]*crhs43;
    rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs12 + DN(3,0)*crhs24 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - N[3]*crhs18 - crhs14*crhs54 - crhs19*crhs55 - crhs19*crhs56 - crhs19*crhs57 - crhs22*crhs54 + crhs40*crhs58;
    rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs12 + DN(3,1)*crhs24 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs27 - N[3]*crhs28 - crhs20*crhs55 - crhs20*crhs56 - crhs20*crhs57 - crhs29*crhs54 - crhs30*crhs54 + crhs41*crhs58;
    rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs12 + DN(3,2)*crhs24 - DN(3,2)*stress[2] + N[3]*crhs33 - N[3]*crhs34 - crhs21*crhs55 - crhs21*crhs56 - crhs21*crhs57 - crhs35*crhs54 - crhs36*crhs54 + crhs42*crhs58;
    rhs[15]=-DN(3,0)*crhs40 - DN(3,1)*crhs41 - DN(3,2)*crhs42 - N[3]*crhs23 - N[3]*crhs43;

}


template<>
void NavierStokes<2>::ComputeGaussPointRHSContribution(array_1d<double,9>& rhs, const element_data& data)
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
    const double tau1 = 1.0/((rho*dyn_tau_coeff)/delta_t + (2*rho*vconv_norm)/h + (4*mu)/(h*h));
    const double tau2 = mu + 0.5*h*vconv_norm;

    const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2];
    const double crhs1 =             rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0));
    const double crhs2 =             DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0);
    const double crhs3 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
    const double crhs4 =             crhs2 + crhs3;
    const double crhs5 =             crhs4*rho;
    const double crhs6 =             N[0]*rho;
    const double crhs7 =             N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0));
    const double crhs8 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
    const double crhs9 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
    const double crhs10 =             crhs2*crhs8 + crhs9*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0));
    const double crhs11 =             pow(c, -2);
    const double crhs12 =             N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]);
    const double crhs13 =             crhs11*crhs12;
    const double crhs14 =             crhs13*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0));
    const double crhs15 =             crhs11*crhs12/rho;
    const double crhs16 =             tau2*(crhs15 + crhs4);
    const double crhs17 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + crhs14 + rho*(crhs10 + crhs7);
    const double crhs18 =             DN(0,0)*crhs17*rho*tau1;
    const double crhs19 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
    const double crhs20 =             crhs13*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1));
    const double crhs21 =             N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1));
    const double crhs22 =             crhs3*crhs9 + crhs8*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1));
    const double crhs23 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs19 + crhs20 + rho*(crhs21 + crhs22);
    const double crhs24 =             DN(0,1)*crhs23*rho*tau1;
    const double crhs25 =             N[0]*crhs11*crhs12;
    const double crhs26 =             crhs17*tau1;
    const double crhs27 =             crhs23*tau1;
    const double crhs28 =             2*DN(0,0)*v(0,0) + 2*DN(0,1)*v(0,1) + 2*DN(1,0)*v(1,0) + 2*DN(1,1)*v(1,1) + 2*DN(2,0)*v(2,0) + 2*DN(2,1)*v(2,1);
    const double crhs29 =             N[1]*rho;
    const double crhs30 =             DN(1,0)*crhs17*rho*tau1;
    const double crhs31 =             DN(1,1)*crhs23*rho*tau1;
    const double crhs32 =             N[1]*crhs11*crhs12;
    const double crhs33 =             N[2]*rho;
    const double crhs34 =             DN(2,0)*crhs17*rho*tau1;
    const double crhs35 =             DN(2,1)*crhs23*rho*tau1;
    const double crhs36 =             N[2]*crhs11*crhs12;

    rhs[0]=DN(0,0)*crhs0 + DN(0,0)*crhs16 - DN(0,0)*crhs5 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - N[0]*crhs14 - crhs10*crhs6 - crhs18*crhs8 - crhs24*crhs8 + crhs25*crhs26 - crhs6*crhs7;
    rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 + DN(0,1)*crhs16 - DN(0,1)*crhs5 - DN(0,1)*stress[1] + N[0]*crhs19 - N[0]*crhs20 - crhs18*crhs9 - crhs21*crhs6 - crhs22*crhs6 - crhs24*crhs9 + crhs25*crhs27;
    rhs[2]=-DN(0,0)*crhs26 - DN(0,1)*crhs27 - N[0]*crhs15 - N[0]*crhs28;
    rhs[3]=DN(1,0)*crhs0 + DN(1,0)*crhs16 - DN(1,0)*crhs5 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - N[1]*crhs14 - crhs10*crhs29 + crhs26*crhs32 - crhs29*crhs7 - crhs30*crhs8 - crhs31*crhs8;
    rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 + DN(1,1)*crhs16 - DN(1,1)*crhs5 - DN(1,1)*stress[1] + N[1]*crhs19 - N[1]*crhs20 - crhs21*crhs29 - crhs22*crhs29 + crhs27*crhs32 - crhs30*crhs9 - crhs31*crhs9;
    rhs[5]=-DN(1,0)*crhs26 - DN(1,1)*crhs27 - N[1]*crhs15 - N[1]*crhs28;
    rhs[6]=DN(2,0)*crhs0 + DN(2,0)*crhs16 - DN(2,0)*crhs5 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - N[2]*crhs14 - crhs10*crhs33 + crhs26*crhs36 - crhs33*crhs7 - crhs34*crhs8 - crhs35*crhs8;
    rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 + DN(2,1)*crhs16 - DN(2,1)*crhs5 - DN(2,1)*stress[1] + N[2]*crhs19 - N[2]*crhs20 - crhs21*crhs33 - crhs22*crhs33 + crhs27*crhs36 - crhs34*crhs9 - crhs35*crhs9;
    rhs[8]=-DN(2,0)*crhs26 - DN(2,1)*crhs27 - N[2]*crhs15 - N[2]*crhs28;

}


template<>
double NavierStokes<3>::SubscaleErrorEstimate(const element_data& data)
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

    // Get constitutive matrix
    //~ const Matrix& C = data.C;

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
    const double tau1 = 1.0/((rho*dyn_tau_coeff)/delta_t + (2*rho*vconv_norm)/h + (4*mu)/(h*h));
    // const double tau2 = mu + 0.5*h*vconv_norm;

    // Gauss point velocity subscale value computation
    const double cv_s_gauss0 =             (N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]) + N[3]*(bdf0*p[3] + bdf1*pn[3] + bdf2*pnn[3]))/pow(c, 2);
    const double cv_s_gauss1 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
    const double cv_s_gauss2 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
    const double cv_s_gauss3 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);

    v_s_gauss[0]=-tau1*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + cv_s_gauss0*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0)) + rho*(-N[0]*f(0,0) + N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) - N[1]*f(1,0) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) - N[2]*f(2,0) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) - N[3]*f(3,0) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)) + cv_s_gauss1*(DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0)) + cv_s_gauss2*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + cv_s_gauss3*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0))));
    v_s_gauss[1]=-tau1*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + cv_s_gauss0*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1)) + rho*(-N[0]*f(0,1) + N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) - N[1]*f(1,1) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) - N[2]*f(2,1) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) - N[3]*f(3,1) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)) + cv_s_gauss1*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + cv_s_gauss2*(DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1)) + cv_s_gauss3*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1))));
    v_s_gauss[2]=-tau1*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + cv_s_gauss0*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2)) + rho*(-N[0]*f(0,2) + N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) - N[1]*f(1,2) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) - N[2]*f(2,2) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) - N[3]*f(3,2) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)) + cv_s_gauss1*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + cv_s_gauss2*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + cv_s_gauss3*(DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2))));

    const double v_gauss_norm = norm_2(v_gauss);
    const double v_s_gauss_norm = norm_2(v_s_gauss);

    return v_s_gauss_norm/v_gauss_norm;
}


template<>
double NavierStokes<2>::SubscaleErrorEstimate(const element_data& data)
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
    const double tau1 = 1.0/((rho*dyn_tau_coeff)/delta_t + (2*rho*vconv_norm)/h + (4*mu)/(h*h));
    // const double tau2 = mu + 0.5*h*vconv_norm;

    // Gauss point velocity subscale value computation
    const double cv_s_gauss0 =             (N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]))/pow(c, 2);
    const double cv_s_gauss1 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
    const double cv_s_gauss2 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);

    v_s_gauss[0]=-tau1*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + cv_s_gauss0*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0)) + rho*(-N[0]*f(0,0) + N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) - N[1]*f(1,0) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) - N[2]*f(2,0) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + cv_s_gauss1*(DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0)) + cv_s_gauss2*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0))));
    v_s_gauss[1]=-tau1*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + cv_s_gauss0*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1)) + rho*(-N[0]*f(0,1) + N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) - N[1]*f(1,1) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) - N[2]*f(2,1) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + cv_s_gauss1*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + cv_s_gauss2*(DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1))));

    const double v_gauss_norm = norm_2(v_gauss);
    const double v_s_gauss_norm = norm_2(v_s_gauss);

    return v_s_gauss_norm/v_gauss_norm;
}

}
