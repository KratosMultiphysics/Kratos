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
        const double clhs1 =             C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
        const double clhs2 =             C(0,3)*DN(0,0);
        const double clhs3 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs2;
        const double clhs4 =             C(0,5)*DN(0,0);
        const double clhs5 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs4;
        const double clhs6 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
        const double clhs7 =             DN(0,0)*clhs6*rho*tau1;
        const double clhs8 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
        const double clhs9 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
        const double clhs10 =             DN(0,0)*clhs6 + DN(0,1)*clhs8 + DN(0,2)*clhs9;
        const double clhs11 =             pow(c, -2);
        const double clhs12 =             N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]) + N[3]*(bdf0*p[3] + bdf1*pn[3] + bdf2*pnn[3]);
        const double clhs13 =             clhs11*clhs12;
        const double clhs14 =             N[0]*clhs13;
        const double clhs15 =             clhs14 + rho*(N[0]*bdf0 + clhs10);
        const double clhs16 =             pow(N[0], 2);
        const double clhs17 =             bdf0*rho;
        const double clhs18 =             N[0]*rho;
        const double clhs19 =             N[0]*clhs11*clhs12*tau1;
        const double clhs20 =             clhs10*clhs18 + clhs13*clhs16 - clhs15*clhs19 + clhs16*clhs17;
        const double clhs21 =             DN(0,0)*DN(0,1);
        const double clhs22 =             clhs21*rho + clhs21*tau2;
        const double clhs23 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs2;
        const double clhs24 =             C(1,3)*DN(0,1);
        const double clhs25 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs24;
        const double clhs26 =             C(3,5)*DN(0,0);
        const double clhs27 =             C(4,5)*DN(0,2);
        const double clhs28 =             C(1,5)*DN(0,1) + clhs26 + clhs27;
        const double clhs29 =             DN(0,1)*clhs15*rho*tau1;
        const double clhs30 =             DN(0,0)*DN(0,2);
        const double clhs31 =             clhs30*rho + clhs30*tau2;
        const double clhs32 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs4;
        const double clhs33 =             C(3,4)*DN(0,1);
        const double clhs34 =             C(2,3)*DN(0,2) + clhs26 + clhs33;
        const double clhs35 =             C(2,5)*DN(0,2);
        const double clhs36 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs35;
        const double clhs37 =             DN(0,2)*clhs15*rho*tau1;
        const double clhs38 =             DN(0,0)*N[0];
        const double clhs39 =             1.0/rho;
        const double clhs40 =             bdf0*clhs11*clhs39*tau2;
        const double clhs41 =             bdf0*clhs11*clhs16;
        const double clhs42 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0);
        const double clhs43 =             N[0]*bdf0*clhs11;
        const double clhs44 =             clhs42*clhs43;
        const double clhs45 =             DN(0,0) + clhs44;
        const double clhs46 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
        const double clhs47 =             clhs43*clhs46;
        const double clhs48 =             DN(0,1) + clhs47;
        const double clhs49 =             DN(0,1)*clhs48*rho*tau1;
        const double clhs50 =             N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2);
        const double clhs51 =             clhs43*clhs50;
        const double clhs52 =             DN(0,2) + clhs51;
        const double clhs53 =             DN(0,2)*clhs52*rho*tau1;
        const double clhs54 =             bdf0*clhs11*clhs16*tau1;
        const double clhs55 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + clhs13*clhs42 - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)) + clhs6*(DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0)) + clhs8*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + clhs9*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)));
        const double clhs56 =             DN(0,0)*DN(1,0);
        const double clhs57 =             clhs56*rho;
        const double clhs58 =             clhs56*tau2;
        const double clhs59 =             N[0]*bdf0*rho;
        const double clhs60 =             N[1]*clhs59;
        const double clhs61 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
        const double clhs62 =             C(0,3)*DN(1,0);
        const double clhs63 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs62;
        const double clhs64 =             C(0,5)*DN(1,0);
        const double clhs65 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs64;
        const double clhs66 =             DN(1,0)*clhs6 + DN(1,1)*clhs8 + DN(1,2)*clhs9;
        const double clhs67 =             clhs18*clhs66;
        const double clhs68 =             N[1]*clhs14;
        const double clhs69 =             N[1]*clhs13;
        const double clhs70 =             clhs69 + rho*(N[1]*bdf0 + clhs66);
        const double clhs71 =             -clhs19*clhs70;
        const double clhs72 =             DN(0,0)*DN(1,1);
        const double clhs73 =             clhs72*rho + clhs72*tau2;
        const double clhs74 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs62;
        const double clhs75 =             C(1,3)*DN(1,1);
        const double clhs76 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs75;
        const double clhs77 =             C(3,5)*DN(1,0);
        const double clhs78 =             C(4,5)*DN(1,2);
        const double clhs79 =             C(1,5)*DN(1,1) + clhs77 + clhs78;
        const double clhs80 =             DN(0,1)*clhs70*rho*tau1;
        const double clhs81 =             DN(0,0)*DN(1,2);
        const double clhs82 =             clhs81*rho + clhs81*tau2;
        const double clhs83 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs64;
        const double clhs84 =             C(3,4)*DN(1,1);
        const double clhs85 =             C(2,3)*DN(1,2) + clhs77 + clhs84;
        const double clhs86 =             C(2,5)*DN(1,2);
        const double clhs87 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs86;
        const double clhs88 =             DN(0,2)*clhs70*rho*tau1;
        const double clhs89 =             DN(0,0)*N[1];
        const double clhs90 =             N[1]*bdf0*clhs11;
        const double clhs91 =             clhs42*clhs90;
        const double clhs92 =             DN(1,0) + clhs91;
        const double clhs93 =             clhs46*clhs90;
        const double clhs94 =             DN(1,1) + clhs93;
        const double clhs95 =             DN(0,1)*clhs94*rho*tau1;
        const double clhs96 =             clhs50*clhs90;
        const double clhs97 =             DN(1,2) + clhs96;
        const double clhs98 =             DN(0,2)*clhs97*rho*tau1;
        const double clhs99 =             N[0]*N[1]*bdf0*clhs11*tau1;
        const double clhs100 =             N[1]*clhs44 - clhs55*clhs99;
        const double clhs101 =             DN(0,0)*DN(2,0);
        const double clhs102 =             clhs101*rho;
        const double clhs103 =             clhs101*tau2;
        const double clhs104 =             N[2]*clhs59;
        const double clhs105 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
        const double clhs106 =             C(0,3)*DN(2,0);
        const double clhs107 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs106;
        const double clhs108 =             C(0,5)*DN(2,0);
        const double clhs109 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs108;
        const double clhs110 =             DN(2,0)*clhs6 + DN(2,1)*clhs8 + DN(2,2)*clhs9;
        const double clhs111 =             clhs110*clhs18;
        const double clhs112 =             N[2]*clhs14;
        const double clhs113 =             N[2]*clhs13;
        const double clhs114 =             clhs113 + rho*(N[2]*bdf0 + clhs110);
        const double clhs115 =             -clhs114*clhs19;
        const double clhs116 =             DN(0,0)*DN(2,1);
        const double clhs117 =             clhs116*rho + clhs116*tau2;
        const double clhs118 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs106;
        const double clhs119 =             C(1,3)*DN(2,1);
        const double clhs120 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs119;
        const double clhs121 =             C(3,5)*DN(2,0);
        const double clhs122 =             C(4,5)*DN(2,2);
        const double clhs123 =             C(1,5)*DN(2,1) + clhs121 + clhs122;
        const double clhs124 =             DN(0,1)*clhs114*rho*tau1;
        const double clhs125 =             DN(0,0)*DN(2,2);
        const double clhs126 =             clhs125*rho + clhs125*tau2;
        const double clhs127 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs108;
        const double clhs128 =             C(3,4)*DN(2,1);
        const double clhs129 =             C(2,3)*DN(2,2) + clhs121 + clhs128;
        const double clhs130 =             C(2,5)*DN(2,2);
        const double clhs131 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs130;
        const double clhs132 =             DN(0,2)*clhs114*rho*tau1;
        const double clhs133 =             DN(0,0)*N[2];
        const double clhs134 =             N[2]*bdf0*clhs11;
        const double clhs135 =             DN(2,0) + clhs134*clhs42;
        const double clhs136 =             DN(2,1) + clhs134*clhs46;
        const double clhs137 =             DN(0,1)*clhs136*rho*tau1;
        const double clhs138 =             DN(2,2) + clhs134*clhs50;
        const double clhs139 =             DN(0,2)*clhs138*rho*tau1;
        const double clhs140 =             N[0]*N[2]*bdf0*clhs11*tau1;
        const double clhs141 =             N[2]*clhs44 - clhs140*clhs55;
        const double clhs142 =             DN(0,0)*DN(3,0);
        const double clhs143 =             clhs142*rho;
        const double clhs144 =             clhs142*tau2;
        const double clhs145 =             N[3]*clhs59;
        const double clhs146 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
        const double clhs147 =             C(0,3)*DN(3,0);
        const double clhs148 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs147;
        const double clhs149 =             C(0,5)*DN(3,0);
        const double clhs150 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs149;
        const double clhs151 =             DN(3,0)*clhs6 + DN(3,1)*clhs8 + DN(3,2)*clhs9;
        const double clhs152 =             clhs151*clhs18;
        const double clhs153 =             N[3]*clhs14;
        const double clhs154 =             N[3]*clhs13 + rho*(N[3]*bdf0 + clhs151);
        const double clhs155 =             -clhs154*clhs19;
        const double clhs156 =             DN(0,0)*DN(3,1);
        const double clhs157 =             clhs156*rho + clhs156*tau2;
        const double clhs158 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs147;
        const double clhs159 =             C(1,3)*DN(3,1);
        const double clhs160 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs159;
        const double clhs161 =             C(3,5)*DN(3,0);
        const double clhs162 =             C(4,5)*DN(3,2);
        const double clhs163 =             C(1,5)*DN(3,1) + clhs161 + clhs162;
        const double clhs164 =             DN(0,1)*clhs154*rho*tau1;
        const double clhs165 =             DN(0,0)*DN(3,2);
        const double clhs166 =             clhs165*rho + clhs165*tau2;
        const double clhs167 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs149;
        const double clhs168 =             C(3,4)*DN(3,1);
        const double clhs169 =             C(2,3)*DN(3,2) + clhs161 + clhs168;
        const double clhs170 =             C(2,5)*DN(3,2);
        const double clhs171 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs170;
        const double clhs172 =             DN(0,2)*clhs154*rho*tau1;
        const double clhs173 =             DN(0,0)*N[3];
        const double clhs174 =             N[3]*bdf0*clhs11;
        const double clhs175 =             DN(3,0) + clhs174*clhs42;
        const double clhs176 =             DN(3,1) + clhs174*clhs46;
        const double clhs177 =             DN(0,1)*clhs176*rho*tau1;
        const double clhs178 =             DN(3,2) + clhs174*clhs50;
        const double clhs179 =             DN(0,2)*clhs178*rho*tau1;
        const double clhs180 =             N[0]*N[3]*bdf0*clhs11*tau1;
        const double clhs181 =             N[3]*clhs44 - clhs180*clhs55;
        const double clhs182 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs24;
        const double clhs183 =             C(0,4)*DN(0,0) + clhs27 + clhs33;
        const double clhs184 =             DN(0,0)*clhs15*rho*tau1;
        const double clhs185 =             pow(DN(0,1), 2);
        const double clhs186 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
        const double clhs187 =             C(1,4)*DN(0,1);
        const double clhs188 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs187;
        const double clhs189 =             DN(0,1)*clhs8*rho*tau1;
        const double clhs190 =             DN(0,1)*DN(0,2);
        const double clhs191 =             clhs190*rho + clhs190*tau2;
        const double clhs192 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs187;
        const double clhs193 =             C(2,4)*DN(0,2);
        const double clhs194 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs193;
        const double clhs195 =             DN(0,1)*N[0];
        const double clhs196 =             DN(0,0)*clhs45*rho*tau1;
        const double clhs197 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + clhs13*clhs46 - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)) + clhs6*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + clhs8*(DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1)) + clhs9*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
        const double clhs198 =             DN(0,1)*DN(1,0);
        const double clhs199 =             clhs198*rho + clhs198*tau2;
        const double clhs200 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs75;
        const double clhs201 =             C(0,4)*DN(1,0) + clhs78 + clhs84;
        const double clhs202 =             DN(0,0)*clhs70*rho*tau1;
        const double clhs203 =             DN(0,1)*DN(1,1);
        const double clhs204 =             clhs203*rho;
        const double clhs205 =             clhs203*tau2;
        const double clhs206 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
        const double clhs207 =             C(1,4)*DN(1,1);
        const double clhs208 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs207;
        const double clhs209 =             DN(0,1)*DN(1,2);
        const double clhs210 =             clhs209*rho + clhs209*tau2;
        const double clhs211 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs207;
        const double clhs212 =             C(2,4)*DN(1,2);
        const double clhs213 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs212;
        const double clhs214 =             DN(0,1)*N[1];
        const double clhs215 =             DN(0,0)*clhs92*rho*tau1;
        const double clhs216 =             N[1]*clhs47 - clhs197*clhs99;
        const double clhs217 =             DN(0,1)*DN(2,0);
        const double clhs218 =             clhs217*rho + clhs217*tau2;
        const double clhs219 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs119;
        const double clhs220 =             C(0,4)*DN(2,0) + clhs122 + clhs128;
        const double clhs221 =             DN(0,0)*clhs114*rho*tau1;
        const double clhs222 =             DN(0,1)*DN(2,1);
        const double clhs223 =             clhs222*rho;
        const double clhs224 =             clhs222*tau2;
        const double clhs225 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
        const double clhs226 =             C(1,4)*DN(2,1);
        const double clhs227 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs226;
        const double clhs228 =             DN(0,1)*DN(2,2);
        const double clhs229 =             clhs228*rho + clhs228*tau2;
        const double clhs230 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs226;
        const double clhs231 =             C(2,4)*DN(2,2);
        const double clhs232 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs231;
        const double clhs233 =             DN(0,1)*N[2];
        const double clhs234 =             DN(0,0)*clhs135*rho*tau1;
        const double clhs235 =             N[2]*clhs47 - clhs140*clhs197;
        const double clhs236 =             DN(0,1)*DN(3,0);
        const double clhs237 =             clhs236*rho + clhs236*tau2;
        const double clhs238 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs159;
        const double clhs239 =             C(0,4)*DN(3,0) + clhs162 + clhs168;
        const double clhs240 =             DN(0,0)*clhs154*rho*tau1;
        const double clhs241 =             DN(0,1)*DN(3,1);
        const double clhs242 =             clhs241*rho;
        const double clhs243 =             clhs241*tau2;
        const double clhs244 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
        const double clhs245 =             C(1,4)*DN(3,1);
        const double clhs246 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs245;
        const double clhs247 =             DN(0,1)*DN(3,2);
        const double clhs248 =             clhs247*rho + clhs247*tau2;
        const double clhs249 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs245;
        const double clhs250 =             C(2,4)*DN(3,2);
        const double clhs251 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs250;
        const double clhs252 =             DN(0,1)*N[3];
        const double clhs253 =             DN(0,0)*clhs175*rho*tau1;
        const double clhs254 =             N[3]*clhs47 - clhs180*clhs197;
        const double clhs255 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs35;
        const double clhs256 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs193;
        const double clhs257 =             pow(DN(0,2), 2);
        const double clhs258 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
        const double clhs259 =             DN(0,2)*clhs9*rho*tau1;
        const double clhs260 =             DN(0,2)*N[0];
        const double clhs261 =             DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + clhs13*clhs50 - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)) + clhs6*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + clhs8*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + clhs9*(DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2)));
        const double clhs262 =             DN(0,2)*DN(1,0);
        const double clhs263 =             clhs262*rho + clhs262*tau2;
        const double clhs264 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs86;
        const double clhs265 =             DN(0,2)*DN(1,1);
        const double clhs266 =             clhs265*rho + clhs265*tau2;
        const double clhs267 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs212;
        const double clhs268 =             DN(0,2)*DN(1,2);
        const double clhs269 =             clhs268*rho;
        const double clhs270 =             clhs268*tau2;
        const double clhs271 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
        const double clhs272 =             DN(0,2)*N[1];
        const double clhs273 =             N[1]*clhs51 - clhs261*clhs99;
        const double clhs274 =             DN(0,2)*DN(2,0);
        const double clhs275 =             clhs274*rho + clhs274*tau2;
        const double clhs276 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs130;
        const double clhs277 =             DN(0,2)*DN(2,1);
        const double clhs278 =             clhs277*rho + clhs277*tau2;
        const double clhs279 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs231;
        const double clhs280 =             DN(0,2)*DN(2,2);
        const double clhs281 =             clhs280*rho;
        const double clhs282 =             clhs280*tau2;
        const double clhs283 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
        const double clhs284 =             DN(0,2)*N[2];
        const double clhs285 =             N[2]*clhs51 - clhs140*clhs261;
        const double clhs286 =             DN(0,2)*DN(3,0);
        const double clhs287 =             clhs286*rho + clhs286*tau2;
        const double clhs288 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs170;
        const double clhs289 =             DN(0,2)*DN(3,1);
        const double clhs290 =             clhs289*rho + clhs289*tau2;
        const double clhs291 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs250;
        const double clhs292 =             DN(0,2)*DN(3,2);
        const double clhs293 =             clhs292*rho;
        const double clhs294 =             clhs292*tau2;
        const double clhs295 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
        const double clhs296 =             DN(0,2)*N[3];
        const double clhs297 =             N[3]*clhs51 - clhs180*clhs261;
        const double clhs298 =             2*N[0];
        const double clhs299 =             clhs15*tau1 + clhs298;
        const double clhs300 =             bdf0*clhs11*clhs39;
        const double clhs301 =             DN(0,0)*tau1;
        const double clhs302 =             DN(0,1)*tau1;
        const double clhs303 =             DN(0,2)*tau1;
        const double clhs304 =             N[0]*bdf0*clhs11*clhs39;
        const double clhs305 =             N[1]*clhs304;
        const double clhs306 =             N[2]*clhs304;
        const double clhs307 =             N[3]*clhs304;
        const double clhs308 =             DN(1,0)*clhs6*rho*tau1;
        const double clhs309 =             N[1]*rho;
        const double clhs310 =             N[1]*clhs11*clhs12*tau1;
        const double clhs311 =             clhs10*clhs309 - clhs15*clhs310 + clhs60 + clhs68;
        const double clhs312 =             DN(1,1)*clhs15*rho*tau1;
        const double clhs313 =             DN(1,2)*clhs15*rho*tau1;
        const double clhs314 =             DN(1,0)*N[0];
        const double clhs315 =             DN(1,1)*clhs48*rho*tau1;
        const double clhs316 =             DN(1,2)*clhs52*rho*tau1;
        const double clhs317 =             pow(DN(1,0), 2);
        const double clhs318 =             pow(N[1], 2);
        const double clhs319 =             clhs13*clhs318 + clhs17*clhs318 + clhs309*clhs66 - clhs310*clhs70;
        const double clhs320 =             DN(1,0)*DN(1,1);
        const double clhs321 =             clhs320*rho + clhs320*tau2;
        const double clhs322 =             DN(1,1)*clhs70*rho*tau1;
        const double clhs323 =             DN(1,0)*DN(1,2);
        const double clhs324 =             clhs323*rho + clhs323*tau2;
        const double clhs325 =             DN(1,2)*clhs70*rho*tau1;
        const double clhs326 =             DN(1,0)*N[1];
        const double clhs327 =             bdf0*clhs11*clhs318;
        const double clhs328 =             DN(1,1)*clhs94*rho*tau1;
        const double clhs329 =             DN(1,2)*clhs97*rho*tau1;
        const double clhs330 =             bdf0*clhs11*clhs318*tau1;
        const double clhs331 =             DN(1,0)*DN(2,0);
        const double clhs332 =             clhs331*rho;
        const double clhs333 =             clhs331*tau2;
        const double clhs334 =             N[1]*bdf0*rho;
        const double clhs335 =             N[2]*clhs334;
        const double clhs336 =             clhs110*clhs309;
        const double clhs337 =             N[2]*clhs69;
        const double clhs338 =             -clhs114*clhs310;
        const double clhs339 =             DN(1,0)*DN(2,1);
        const double clhs340 =             clhs339*rho + clhs339*tau2;
        const double clhs341 =             DN(1,1)*clhs114*rho*tau1;
        const double clhs342 =             DN(1,0)*DN(2,2);
        const double clhs343 =             clhs342*rho + clhs342*tau2;
        const double clhs344 =             DN(1,2)*clhs114*rho*tau1;
        const double clhs345 =             DN(1,0)*N[2];
        const double clhs346 =             DN(1,1)*clhs136*rho*tau1;
        const double clhs347 =             DN(1,2)*clhs138*rho*tau1;
        const double clhs348 =             N[1]*N[2]*bdf0*clhs11*tau1;
        const double clhs349 =             N[2]*clhs91 - clhs348*clhs55;
        const double clhs350 =             DN(1,0)*DN(3,0);
        const double clhs351 =             clhs350*rho;
        const double clhs352 =             clhs350*tau2;
        const double clhs353 =             N[3]*clhs334;
        const double clhs354 =             clhs151*clhs309;
        const double clhs355 =             N[3]*clhs69;
        const double clhs356 =             -clhs154*clhs310;
        const double clhs357 =             DN(1,0)*DN(3,1);
        const double clhs358 =             clhs357*rho + clhs357*tau2;
        const double clhs359 =             DN(1,1)*clhs154*rho*tau1;
        const double clhs360 =             DN(1,0)*DN(3,2);
        const double clhs361 =             clhs360*rho + clhs360*tau2;
        const double clhs362 =             DN(1,2)*clhs154*rho*tau1;
        const double clhs363 =             DN(1,0)*N[3];
        const double clhs364 =             DN(1,1)*clhs176*rho*tau1;
        const double clhs365 =             DN(1,2)*clhs178*rho*tau1;
        const double clhs366 =             N[1]*N[3]*bdf0*clhs11*tau1;
        const double clhs367 =             N[3]*clhs91 - clhs366*clhs55;
        const double clhs368 =             DN(1,0)*clhs15*rho*tau1;
        const double clhs369 =             DN(1,1)*clhs8*rho*tau1;
        const double clhs370 =             DN(1,1)*N[0];
        const double clhs371 =             DN(1,0)*clhs45*rho*tau1;
        const double clhs372 =             DN(1,0)*clhs70*rho*tau1;
        const double clhs373 =             pow(DN(1,1), 2);
        const double clhs374 =             DN(1,1)*DN(1,2);
        const double clhs375 =             clhs374*rho + clhs374*tau2;
        const double clhs376 =             DN(1,1)*N[1];
        const double clhs377 =             DN(1,0)*clhs92*rho*tau1;
        const double clhs378 =             DN(1,1)*DN(2,0);
        const double clhs379 =             clhs378*rho + clhs378*tau2;
        const double clhs380 =             DN(1,0)*clhs114*rho*tau1;
        const double clhs381 =             DN(1,1)*DN(2,1);
        const double clhs382 =             clhs381*rho;
        const double clhs383 =             clhs381*tau2;
        const double clhs384 =             DN(1,1)*DN(2,2);
        const double clhs385 =             clhs384*rho + clhs384*tau2;
        const double clhs386 =             DN(1,1)*N[2];
        const double clhs387 =             DN(1,0)*clhs135*rho*tau1;
        const double clhs388 =             N[2]*clhs93 - clhs197*clhs348;
        const double clhs389 =             DN(1,1)*DN(3,0);
        const double clhs390 =             clhs389*rho + clhs389*tau2;
        const double clhs391 =             DN(1,0)*clhs154*rho*tau1;
        const double clhs392 =             DN(1,1)*DN(3,1);
        const double clhs393 =             clhs392*rho;
        const double clhs394 =             clhs392*tau2;
        const double clhs395 =             DN(1,1)*DN(3,2);
        const double clhs396 =             clhs395*rho + clhs395*tau2;
        const double clhs397 =             DN(1,1)*N[3];
        const double clhs398 =             DN(1,0)*clhs175*rho*tau1;
        const double clhs399 =             N[3]*clhs93 - clhs197*clhs366;
        const double clhs400 =             DN(1,2)*clhs9*rho*tau1;
        const double clhs401 =             DN(1,2)*N[0];
        const double clhs402 =             pow(DN(1,2), 2);
        const double clhs403 =             DN(1,2)*N[1];
        const double clhs404 =             DN(1,2)*DN(2,0);
        const double clhs405 =             clhs404*rho + clhs404*tau2;
        const double clhs406 =             DN(1,2)*DN(2,1);
        const double clhs407 =             clhs406*rho + clhs406*tau2;
        const double clhs408 =             DN(1,2)*DN(2,2);
        const double clhs409 =             clhs408*rho;
        const double clhs410 =             clhs408*tau2;
        const double clhs411 =             DN(1,2)*N[2];
        const double clhs412 =             N[2]*clhs96 - clhs261*clhs348;
        const double clhs413 =             DN(1,2)*DN(3,0);
        const double clhs414 =             clhs413*rho + clhs413*tau2;
        const double clhs415 =             DN(1,2)*DN(3,1);
        const double clhs416 =             clhs415*rho + clhs415*tau2;
        const double clhs417 =             DN(1,2)*DN(3,2);
        const double clhs418 =             clhs417*rho;
        const double clhs419 =             clhs417*tau2;
        const double clhs420 =             DN(1,2)*N[3];
        const double clhs421 =             N[3]*clhs96 - clhs261*clhs366;
        const double clhs422 =             2*N[1];
        const double clhs423 =             DN(1,0)*tau1;
        const double clhs424 =             DN(1,1)*tau1;
        const double clhs425 =             DN(1,2)*tau1;
        const double clhs426 =             clhs422 + clhs70*tau1;
        const double clhs427 =             N[1]*bdf0*clhs11*clhs39;
        const double clhs428 =             N[2]*clhs427;
        const double clhs429 =             N[3]*clhs427;
        const double clhs430 =             DN(2,0)*clhs6*rho*tau1;
        const double clhs431 =             N[2]*rho;
        const double clhs432 =             N[2]*clhs11*clhs12*tau1;
        const double clhs433 =             clhs10*clhs431 + clhs104 + clhs112 - clhs15*clhs432;
        const double clhs434 =             DN(2,1)*clhs15*rho*tau1;
        const double clhs435 =             DN(2,2)*clhs15*rho*tau1;
        const double clhs436 =             DN(2,0)*N[0];
        const double clhs437 =             DN(2,1)*clhs48*rho*tau1;
        const double clhs438 =             DN(2,2)*clhs52*rho*tau1;
        const double clhs439 =             clhs335 + clhs337 + clhs431*clhs66 - clhs432*clhs70;
        const double clhs440 =             DN(2,1)*clhs70*rho*tau1;
        const double clhs441 =             DN(2,2)*clhs70*rho*tau1;
        const double clhs442 =             DN(2,0)*N[1];
        const double clhs443 =             DN(2,1)*clhs94*rho*tau1;
        const double clhs444 =             DN(2,2)*clhs97*rho*tau1;
        const double clhs445 =             pow(DN(2,0), 2);
        const double clhs446 =             pow(N[2], 2);
        const double clhs447 =             clhs110*clhs431 - clhs114*clhs432 + clhs13*clhs446 + clhs17*clhs446;
        const double clhs448 =             DN(2,0)*DN(2,1);
        const double clhs449 =             clhs448*rho + clhs448*tau2;
        const double clhs450 =             DN(2,1)*clhs114*rho*tau1;
        const double clhs451 =             DN(2,0)*DN(2,2);
        const double clhs452 =             clhs451*rho + clhs451*tau2;
        const double clhs453 =             DN(2,2)*clhs114*rho*tau1;
        const double clhs454 =             DN(2,0)*N[2];
        const double clhs455 =             bdf0*clhs11*clhs446;
        const double clhs456 =             DN(2,1)*clhs136*rho*tau1;
        const double clhs457 =             DN(2,2)*clhs138*rho*tau1;
        const double clhs458 =             bdf0*clhs11*clhs446*tau1;
        const double clhs459 =             DN(2,0)*DN(3,0);
        const double clhs460 =             clhs459*rho;
        const double clhs461 =             clhs459*tau2;
        const double clhs462 =             N[2]*N[3]*bdf0;
        const double clhs463 =             clhs462*rho;
        const double clhs464 =             clhs151*clhs431;
        const double clhs465 =             N[3]*clhs113;
        const double clhs466 =             -clhs154*clhs432;
        const double clhs467 =             DN(2,0)*DN(3,1);
        const double clhs468 =             clhs467*rho + clhs467*tau2;
        const double clhs469 =             DN(2,1)*clhs154*rho*tau1;
        const double clhs470 =             DN(2,0)*DN(3,2);
        const double clhs471 =             clhs470*rho + clhs470*tau2;
        const double clhs472 =             DN(2,2)*clhs154*rho*tau1;
        const double clhs473 =             DN(2,0)*N[3];
        const double clhs474 =             DN(2,1)*clhs176*rho*tau1;
        const double clhs475 =             DN(2,2)*clhs178*rho*tau1;
        const double clhs476 =             N[2]*N[3]*bdf0*clhs11;
        const double clhs477 =             N[2]*N[3]*bdf0*clhs11*tau1;
        const double clhs478 =             clhs42*clhs476 - clhs477*clhs55;
        const double clhs479 =             DN(2,0)*clhs15*rho*tau1;
        const double clhs480 =             DN(2,1)*clhs8*rho*tau1;
        const double clhs481 =             DN(2,1)*N[0];
        const double clhs482 =             DN(2,0)*clhs45*rho*tau1;
        const double clhs483 =             DN(2,0)*clhs70*rho*tau1;
        const double clhs484 =             DN(2,1)*N[1];
        const double clhs485 =             DN(2,0)*clhs92*rho*tau1;
        const double clhs486 =             DN(2,0)*clhs114*rho*tau1;
        const double clhs487 =             pow(DN(2,1), 2);
        const double clhs488 =             DN(2,1)*DN(2,2);
        const double clhs489 =             clhs488*rho + clhs488*tau2;
        const double clhs490 =             DN(2,1)*N[2];
        const double clhs491 =             DN(2,0)*clhs135*rho*tau1;
        const double clhs492 =             DN(2,1)*DN(3,0);
        const double clhs493 =             clhs492*rho + clhs492*tau2;
        const double clhs494 =             DN(2,0)*clhs154*rho*tau1;
        const double clhs495 =             DN(2,1)*DN(3,1);
        const double clhs496 =             clhs495*rho;
        const double clhs497 =             clhs495*tau2;
        const double clhs498 =             DN(2,1)*DN(3,2);
        const double clhs499 =             clhs498*rho + clhs498*tau2;
        const double clhs500 =             DN(2,1)*N[3];
        const double clhs501 =             DN(2,0)*clhs175*rho*tau1;
        const double clhs502 =             -clhs197*clhs477 + clhs46*clhs476;
        const double clhs503 =             DN(2,2)*clhs9*rho*tau1;
        const double clhs504 =             DN(2,2)*N[0];
        const double clhs505 =             DN(2,2)*N[1];
        const double clhs506 =             pow(DN(2,2), 2);
        const double clhs507 =             DN(2,2)*N[2];
        const double clhs508 =             DN(2,2)*DN(3,0);
        const double clhs509 =             clhs508*rho + clhs508*tau2;
        const double clhs510 =             DN(2,2)*DN(3,1);
        const double clhs511 =             clhs510*rho + clhs510*tau2;
        const double clhs512 =             DN(2,2)*DN(3,2);
        const double clhs513 =             clhs512*rho;
        const double clhs514 =             clhs512*tau2;
        const double clhs515 =             DN(2,2)*N[3];
        const double clhs516 =             -clhs261*clhs477 + clhs476*clhs50;
        const double clhs517 =             2*N[2];
        const double clhs518 =             DN(2,0)*tau1;
        const double clhs519 =             DN(2,1)*tau1;
        const double clhs520 =             DN(2,2)*tau1;
        const double clhs521 =             clhs114*tau1 + clhs517;
        const double clhs522 =             clhs11*clhs39*clhs462;
        const double clhs523 =             DN(3,0)*clhs6*rho*tau1;
        const double clhs524 =             N[3]*rho;
        const double clhs525 =             N[3]*clhs11*clhs12*tau1;
        const double clhs526 =             clhs10*clhs524 + clhs145 - clhs15*clhs525 + clhs153;
        const double clhs527 =             DN(3,1)*clhs15*rho*tau1;
        const double clhs528 =             DN(3,2)*clhs15*rho*tau1;
        const double clhs529 =             DN(3,0)*N[0];
        const double clhs530 =             DN(3,1)*clhs48*rho*tau1;
        const double clhs531 =             DN(3,2)*clhs52*rho*tau1;
        const double clhs532 =             clhs353 + clhs355 + clhs524*clhs66 - clhs525*clhs70;
        const double clhs533 =             DN(3,1)*clhs70*rho*tau1;
        const double clhs534 =             DN(3,2)*clhs70*rho*tau1;
        const double clhs535 =             DN(3,0)*N[1];
        const double clhs536 =             DN(3,1)*clhs94*rho*tau1;
        const double clhs537 =             DN(3,2)*clhs97*rho*tau1;
        const double clhs538 =             clhs110*clhs524 - clhs114*clhs525 + clhs463 + clhs465;
        const double clhs539 =             DN(3,1)*clhs114*rho*tau1;
        const double clhs540 =             DN(3,2)*clhs114*rho*tau1;
        const double clhs541 =             DN(3,0)*N[2];
        const double clhs542 =             DN(3,1)*clhs136*rho*tau1;
        const double clhs543 =             DN(3,2)*clhs138*rho*tau1;
        const double clhs544 =             pow(DN(3,0), 2);
        const double clhs545 =             pow(N[3], 2);
        const double clhs546 =             clhs13*clhs545 + clhs151*clhs524 - clhs154*clhs525 + clhs17*clhs545;
        const double clhs547 =             DN(3,0)*DN(3,1);
        const double clhs548 =             clhs547*rho + clhs547*tau2;
        const double clhs549 =             DN(3,1)*clhs154*rho*tau1;
        const double clhs550 =             DN(3,0)*DN(3,2);
        const double clhs551 =             clhs550*rho + clhs550*tau2;
        const double clhs552 =             DN(3,2)*clhs154*rho*tau1;
        const double clhs553 =             DN(3,0)*N[3];
        const double clhs554 =             bdf0*clhs11*clhs545;
        const double clhs555 =             DN(3,1)*clhs176*rho*tau1;
        const double clhs556 =             DN(3,2)*clhs178*rho*tau1;
        const double clhs557 =             bdf0*clhs11*clhs545*tau1;
        const double clhs558 =             DN(3,0)*clhs15*rho*tau1;
        const double clhs559 =             DN(3,1)*clhs8*rho*tau1;
        const double clhs560 =             DN(3,1)*N[0];
        const double clhs561 =             DN(3,0)*clhs45*rho*tau1;
        const double clhs562 =             DN(3,0)*clhs70*rho*tau1;
        const double clhs563 =             DN(3,1)*N[1];
        const double clhs564 =             DN(3,0)*clhs92*rho*tau1;
        const double clhs565 =             DN(3,0)*clhs114*rho*tau1;
        const double clhs566 =             DN(3,1)*N[2];
        const double clhs567 =             DN(3,0)*clhs135*rho*tau1;
        const double clhs568 =             DN(3,0)*clhs154*rho*tau1;
        const double clhs569 =             pow(DN(3,1), 2);
        const double clhs570 =             DN(3,1)*DN(3,2);
        const double clhs571 =             clhs570*rho + clhs570*tau2;
        const double clhs572 =             DN(3,1)*N[3];
        const double clhs573 =             DN(3,0)*clhs175*rho*tau1;
        const double clhs574 =             DN(3,2)*clhs9*rho*tau1;
        const double clhs575 =             DN(3,2)*N[0];
        const double clhs576 =             DN(3,2)*N[1];
        const double clhs577 =             DN(3,2)*N[2];
        const double clhs578 =             pow(DN(3,2), 2);
        const double clhs579 =             DN(3,2)*N[3];
        const double clhs580 =             2*N[3];
        const double clhs581 =             DN(3,0)*tau1;
        const double clhs582 =             DN(3,1)*tau1;
        const double clhs583 =             DN(3,2)*tau1;
        const double clhs584 =             clhs154*tau1 + clhs580;

        lhs(0,0)=DN(0,0)*clhs1 + DN(0,1)*clhs3 + DN(0,2)*clhs5 + clhs0*rho + clhs0*tau2 + clhs15*clhs7 + clhs20;
        lhs(0,1)=DN(0,0)*clhs23 + DN(0,1)*clhs25 + DN(0,2)*clhs28 + clhs22 + clhs29*clhs6;
        lhs(0,2)=DN(0,0)*clhs32 + DN(0,1)*clhs34 + DN(0,2)*clhs36 + clhs31 + clhs37*clhs6;
        lhs(0,3)=-clhs19*clhs45 + clhs38*clhs40 - clhs38 + clhs41*clhs42 + clhs45*clhs7 + clhs49*clhs6 + clhs53*clhs6 - clhs54*clhs55;
        lhs(0,4)=DN(0,0)*clhs61 + DN(0,1)*clhs63 + DN(0,2)*clhs65 + clhs57 + clhs58 + clhs60 + clhs67 + clhs68 + clhs7*clhs70 + clhs71;
        lhs(0,5)=DN(0,0)*clhs74 + DN(0,1)*clhs76 + DN(0,2)*clhs79 + clhs6*clhs80 + clhs73;
        lhs(0,6)=DN(0,0)*clhs83 + DN(0,1)*clhs85 + DN(0,2)*clhs87 + clhs6*clhs88 + clhs82;
        lhs(0,7)=clhs100 - clhs19*clhs92 + clhs40*clhs89 + clhs6*clhs95 + clhs6*clhs98 + clhs7*clhs92 - clhs89;
        lhs(0,8)=DN(0,0)*clhs105 + DN(0,1)*clhs107 + DN(0,2)*clhs109 + clhs102 + clhs103 + clhs104 + clhs111 + clhs112 + clhs114*clhs7 + clhs115;
        lhs(0,9)=DN(0,0)*clhs118 + DN(0,1)*clhs120 + DN(0,2)*clhs123 + clhs117 + clhs124*clhs6;
        lhs(0,10)=DN(0,0)*clhs127 + DN(0,1)*clhs129 + DN(0,2)*clhs131 + clhs126 + clhs132*clhs6;
        lhs(0,11)=clhs133*clhs40 - clhs133 - clhs135*clhs19 + clhs135*clhs7 + clhs137*clhs6 + clhs139*clhs6 + clhs141;
        lhs(0,12)=DN(0,0)*clhs146 + DN(0,1)*clhs148 + DN(0,2)*clhs150 + clhs143 + clhs144 + clhs145 + clhs152 + clhs153 + clhs154*clhs7 + clhs155;
        lhs(0,13)=DN(0,0)*clhs158 + DN(0,1)*clhs160 + DN(0,2)*clhs163 + clhs157 + clhs164*clhs6;
        lhs(0,14)=DN(0,0)*clhs167 + DN(0,1)*clhs169 + DN(0,2)*clhs171 + clhs166 + clhs172*clhs6;
        lhs(0,15)=clhs173*clhs40 - clhs173 - clhs175*clhs19 + clhs175*clhs7 + clhs177*clhs6 + clhs179*clhs6 + clhs181;
        lhs(1,0)=DN(0,0)*clhs3 + DN(0,1)*clhs182 + DN(0,2)*clhs183 + clhs184*clhs8 + clhs22;
        lhs(1,1)=DN(0,0)*clhs25 + DN(0,1)*clhs186 + DN(0,2)*clhs188 + clhs15*clhs189 + clhs185*rho + clhs185*tau2 + clhs20;
        lhs(1,2)=DN(0,0)*clhs34 + DN(0,1)*clhs192 + DN(0,2)*clhs194 + clhs191 + clhs37*clhs8;
        lhs(1,3)=clhs189*clhs48 - clhs19*clhs48 + clhs195*clhs40 - clhs195 + clhs196*clhs8 - clhs197*clhs54 + clhs41*clhs46 + clhs53*clhs8;
        lhs(1,4)=DN(0,0)*clhs63 + DN(0,1)*clhs200 + DN(0,2)*clhs201 + clhs199 + clhs202*clhs8;
        lhs(1,5)=DN(0,0)*clhs76 + DN(0,1)*clhs206 + DN(0,2)*clhs208 + clhs189*clhs70 + clhs204 + clhs205 + clhs60 + clhs67 + clhs68 + clhs71;
        lhs(1,6)=DN(0,0)*clhs85 + DN(0,1)*clhs211 + DN(0,2)*clhs213 + clhs210 + clhs8*clhs88;
        lhs(1,7)=clhs189*clhs94 - clhs19*clhs94 + clhs214*clhs40 - clhs214 + clhs215*clhs8 + clhs216 + clhs8*clhs98;
        lhs(1,8)=DN(0,0)*clhs107 + DN(0,1)*clhs219 + DN(0,2)*clhs220 + clhs218 + clhs221*clhs8;
        lhs(1,9)=DN(0,0)*clhs120 + DN(0,1)*clhs225 + DN(0,2)*clhs227 + clhs104 + clhs111 + clhs112 + clhs114*clhs189 + clhs115 + clhs223 + clhs224;
        lhs(1,10)=DN(0,0)*clhs129 + DN(0,1)*clhs230 + DN(0,2)*clhs232 + clhs132*clhs8 + clhs229;
        lhs(1,11)=clhs136*clhs189 - clhs136*clhs19 + clhs139*clhs8 + clhs233*clhs40 - clhs233 + clhs234*clhs8 + clhs235;
        lhs(1,12)=DN(0,0)*clhs148 + DN(0,1)*clhs238 + DN(0,2)*clhs239 + clhs237 + clhs240*clhs8;
        lhs(1,13)=DN(0,0)*clhs160 + DN(0,1)*clhs244 + DN(0,2)*clhs246 + clhs145 + clhs152 + clhs153 + clhs154*clhs189 + clhs155 + clhs242 + clhs243;
        lhs(1,14)=DN(0,0)*clhs169 + DN(0,1)*clhs249 + DN(0,2)*clhs251 + clhs172*clhs8 + clhs248;
        lhs(1,15)=clhs176*clhs189 - clhs176*clhs19 + clhs179*clhs8 + clhs252*clhs40 - clhs252 + clhs253*clhs8 + clhs254;
        lhs(2,0)=DN(0,0)*clhs5 + DN(0,1)*clhs183 + DN(0,2)*clhs255 + clhs184*clhs9 + clhs31;
        lhs(2,1)=DN(0,0)*clhs28 + DN(0,1)*clhs188 + DN(0,2)*clhs256 + clhs191 + clhs29*clhs9;
        lhs(2,2)=DN(0,0)*clhs36 + DN(0,1)*clhs194 + DN(0,2)*clhs258 + clhs15*clhs259 + clhs20 + clhs257*rho + clhs257*tau2;
        lhs(2,3)=-clhs19*clhs52 + clhs196*clhs9 + clhs259*clhs52 + clhs260*clhs40 - clhs260 - clhs261*clhs54 + clhs41*clhs50 + clhs49*clhs9;
        lhs(2,4)=DN(0,0)*clhs65 + DN(0,1)*clhs201 + DN(0,2)*clhs264 + clhs202*clhs9 + clhs263;
        lhs(2,5)=DN(0,0)*clhs79 + DN(0,1)*clhs208 + DN(0,2)*clhs267 + clhs266 + clhs80*clhs9;
        lhs(2,6)=DN(0,0)*clhs87 + DN(0,1)*clhs213 + DN(0,2)*clhs271 + clhs259*clhs70 + clhs269 + clhs270 + clhs60 + clhs67 + clhs68 + clhs71;
        lhs(2,7)=-clhs19*clhs97 + clhs215*clhs9 + clhs259*clhs97 + clhs272*clhs40 - clhs272 + clhs273 + clhs9*clhs95;
        lhs(2,8)=DN(0,0)*clhs109 + DN(0,1)*clhs220 + DN(0,2)*clhs276 + clhs221*clhs9 + clhs275;
        lhs(2,9)=DN(0,0)*clhs123 + DN(0,1)*clhs227 + DN(0,2)*clhs279 + clhs124*clhs9 + clhs278;
        lhs(2,10)=DN(0,0)*clhs131 + DN(0,1)*clhs232 + DN(0,2)*clhs283 + clhs104 + clhs111 + clhs112 + clhs114*clhs259 + clhs115 + clhs281 + clhs282;
        lhs(2,11)=clhs137*clhs9 - clhs138*clhs19 + clhs138*clhs259 + clhs234*clhs9 + clhs284*clhs40 - clhs284 + clhs285;
        lhs(2,12)=DN(0,0)*clhs150 + DN(0,1)*clhs239 + DN(0,2)*clhs288 + clhs240*clhs9 + clhs287;
        lhs(2,13)=DN(0,0)*clhs163 + DN(0,1)*clhs246 + DN(0,2)*clhs291 + clhs164*clhs9 + clhs290;
        lhs(2,14)=DN(0,0)*clhs171 + DN(0,1)*clhs251 + DN(0,2)*clhs295 + clhs145 + clhs152 + clhs153 + clhs154*clhs259 + clhs155 + clhs293 + clhs294;
        lhs(2,15)=clhs177*clhs9 - clhs178*clhs19 + clhs178*clhs259 + clhs253*clhs9 + clhs296*clhs40 - clhs296 + clhs297;
        lhs(3,0)=DN(0,0)*clhs299;
        lhs(3,1)=DN(0,1)*clhs299;
        lhs(3,2)=DN(0,2)*clhs299;
        lhs(3,3)=clhs16*clhs300 + clhs301*clhs45 + clhs302*clhs48 + clhs303*clhs52;
        lhs(3,4)=DN(1,0)*clhs298 + clhs301*clhs70;
        lhs(3,5)=DN(1,1)*clhs298 + clhs302*clhs70;
        lhs(3,6)=DN(1,2)*clhs298 + clhs303*clhs70;
        lhs(3,7)=clhs301*clhs92 + clhs302*clhs94 + clhs303*clhs97 + clhs305;
        lhs(3,8)=DN(2,0)*clhs298 + clhs114*clhs301;
        lhs(3,9)=DN(2,1)*clhs298 + clhs114*clhs302;
        lhs(3,10)=DN(2,2)*clhs298 + clhs114*clhs303;
        lhs(3,11)=clhs135*clhs301 + clhs136*clhs302 + clhs138*clhs303 + clhs306;
        lhs(3,12)=DN(3,0)*clhs298 + clhs154*clhs301;
        lhs(3,13)=DN(3,1)*clhs298 + clhs154*clhs302;
        lhs(3,14)=DN(3,2)*clhs298 + clhs154*clhs303;
        lhs(3,15)=clhs175*clhs301 + clhs176*clhs302 + clhs178*clhs303 + clhs307;
        lhs(4,0)=DN(1,0)*clhs1 + DN(1,1)*clhs3 + DN(1,2)*clhs5 + clhs15*clhs308 + clhs311 + clhs57 + clhs58;
        lhs(4,1)=DN(1,0)*clhs23 + DN(1,1)*clhs25 + DN(1,2)*clhs28 + clhs199 + clhs312*clhs6;
        lhs(4,2)=DN(1,0)*clhs32 + DN(1,1)*clhs34 + DN(1,2)*clhs36 + clhs263 + clhs313*clhs6;
        lhs(4,3)=clhs100 + clhs308*clhs45 - clhs310*clhs45 + clhs314*clhs40 - clhs314 + clhs315*clhs6 + clhs316*clhs6;
        lhs(4,4)=DN(1,0)*clhs61 + DN(1,1)*clhs63 + DN(1,2)*clhs65 + clhs308*clhs70 + clhs317*rho + clhs317*tau2 + clhs319;
        lhs(4,5)=DN(1,0)*clhs74 + DN(1,1)*clhs76 + DN(1,2)*clhs79 + clhs321 + clhs322*clhs6;
        lhs(4,6)=DN(1,0)*clhs83 + DN(1,1)*clhs85 + DN(1,2)*clhs87 + clhs324 + clhs325*clhs6;
        lhs(4,7)=clhs308*clhs92 - clhs310*clhs92 + clhs326*clhs40 - clhs326 + clhs327*clhs42 + clhs328*clhs6 + clhs329*clhs6 - clhs330*clhs55;
        lhs(4,8)=DN(1,0)*clhs105 + DN(1,1)*clhs107 + DN(1,2)*clhs109 + clhs114*clhs308 + clhs332 + clhs333 + clhs335 + clhs336 + clhs337 + clhs338;
        lhs(4,9)=DN(1,0)*clhs118 + DN(1,1)*clhs120 + DN(1,2)*clhs123 + clhs340 + clhs341*clhs6;
        lhs(4,10)=DN(1,0)*clhs127 + DN(1,1)*clhs129 + DN(1,2)*clhs131 + clhs343 + clhs344*clhs6;
        lhs(4,11)=clhs135*clhs308 - clhs135*clhs310 + clhs345*clhs40 - clhs345 + clhs346*clhs6 + clhs347*clhs6 + clhs349;
        lhs(4,12)=DN(1,0)*clhs146 + DN(1,1)*clhs148 + DN(1,2)*clhs150 + clhs154*clhs308 + clhs351 + clhs352 + clhs353 + clhs354 + clhs355 + clhs356;
        lhs(4,13)=DN(1,0)*clhs158 + DN(1,1)*clhs160 + DN(1,2)*clhs163 + clhs358 + clhs359*clhs6;
        lhs(4,14)=DN(1,0)*clhs167 + DN(1,1)*clhs169 + DN(1,2)*clhs171 + clhs361 + clhs362*clhs6;
        lhs(4,15)=clhs175*clhs308 - clhs175*clhs310 + clhs363*clhs40 - clhs363 + clhs364*clhs6 + clhs365*clhs6 + clhs367;
        lhs(5,0)=DN(1,0)*clhs3 + DN(1,1)*clhs182 + DN(1,2)*clhs183 + clhs368*clhs8 + clhs73;
        lhs(5,1)=DN(1,0)*clhs25 + DN(1,1)*clhs186 + DN(1,2)*clhs188 + clhs15*clhs369 + clhs204 + clhs205 + clhs311;
        lhs(5,2)=DN(1,0)*clhs34 + DN(1,1)*clhs192 + DN(1,2)*clhs194 + clhs266 + clhs313*clhs8;
        lhs(5,3)=clhs216 - clhs310*clhs48 + clhs316*clhs8 + clhs369*clhs48 + clhs370*clhs40 - clhs370 + clhs371*clhs8;
        lhs(5,4)=DN(1,0)*clhs63 + DN(1,1)*clhs200 + DN(1,2)*clhs201 + clhs321 + clhs372*clhs8;
        lhs(5,5)=DN(1,0)*clhs76 + DN(1,1)*clhs206 + DN(1,2)*clhs208 + clhs319 + clhs369*clhs70 + clhs373*rho + clhs373*tau2;
        lhs(5,6)=DN(1,0)*clhs85 + DN(1,1)*clhs211 + DN(1,2)*clhs213 + clhs325*clhs8 + clhs375;
        lhs(5,7)=-clhs197*clhs330 - clhs310*clhs94 + clhs327*clhs46 + clhs329*clhs8 + clhs369*clhs94 + clhs376*clhs40 - clhs376 + clhs377*clhs8;
        lhs(5,8)=DN(1,0)*clhs107 + DN(1,1)*clhs219 + DN(1,2)*clhs220 + clhs379 + clhs380*clhs8;
        lhs(5,9)=DN(1,0)*clhs120 + DN(1,1)*clhs225 + DN(1,2)*clhs227 + clhs114*clhs369 + clhs335 + clhs336 + clhs337 + clhs338 + clhs382 + clhs383;
        lhs(5,10)=DN(1,0)*clhs129 + DN(1,1)*clhs230 + DN(1,2)*clhs232 + clhs344*clhs8 + clhs385;
        lhs(5,11)=-clhs136*clhs310 + clhs136*clhs369 + clhs347*clhs8 + clhs386*clhs40 - clhs386 + clhs387*clhs8 + clhs388;
        lhs(5,12)=DN(1,0)*clhs148 + DN(1,1)*clhs238 + DN(1,2)*clhs239 + clhs390 + clhs391*clhs8;
        lhs(5,13)=DN(1,0)*clhs160 + DN(1,1)*clhs244 + DN(1,2)*clhs246 + clhs154*clhs369 + clhs353 + clhs354 + clhs355 + clhs356 + clhs393 + clhs394;
        lhs(5,14)=DN(1,0)*clhs169 + DN(1,1)*clhs249 + DN(1,2)*clhs251 + clhs362*clhs8 + clhs396;
        lhs(5,15)=-clhs176*clhs310 + clhs176*clhs369 + clhs365*clhs8 + clhs397*clhs40 - clhs397 + clhs398*clhs8 + clhs399;
        lhs(6,0)=DN(1,0)*clhs5 + DN(1,1)*clhs183 + DN(1,2)*clhs255 + clhs368*clhs9 + clhs82;
        lhs(6,1)=DN(1,0)*clhs28 + DN(1,1)*clhs188 + DN(1,2)*clhs256 + clhs210 + clhs312*clhs9;
        lhs(6,2)=DN(1,0)*clhs36 + DN(1,1)*clhs194 + DN(1,2)*clhs258 + clhs15*clhs400 + clhs269 + clhs270 + clhs311;
        lhs(6,3)=clhs273 - clhs310*clhs52 + clhs315*clhs9 + clhs371*clhs9 + clhs40*clhs401 + clhs400*clhs52 - clhs401;
        lhs(6,4)=DN(1,0)*clhs65 + DN(1,1)*clhs201 + DN(1,2)*clhs264 + clhs324 + clhs372*clhs9;
        lhs(6,5)=DN(1,0)*clhs79 + DN(1,1)*clhs208 + DN(1,2)*clhs267 + clhs322*clhs9 + clhs375;
        lhs(6,6)=DN(1,0)*clhs87 + DN(1,1)*clhs213 + DN(1,2)*clhs271 + clhs319 + clhs400*clhs70 + clhs402*rho + clhs402*tau2;
        lhs(6,7)=-clhs261*clhs330 - clhs310*clhs97 + clhs327*clhs50 + clhs328*clhs9 + clhs377*clhs9 + clhs40*clhs403 + clhs400*clhs97 - clhs403;
        lhs(6,8)=DN(1,0)*clhs109 + DN(1,1)*clhs220 + DN(1,2)*clhs276 + clhs380*clhs9 + clhs405;
        lhs(6,9)=DN(1,0)*clhs123 + DN(1,1)*clhs227 + DN(1,2)*clhs279 + clhs341*clhs9 + clhs407;
        lhs(6,10)=DN(1,0)*clhs131 + DN(1,1)*clhs232 + DN(1,2)*clhs283 + clhs114*clhs400 + clhs335 + clhs336 + clhs337 + clhs338 + clhs409 + clhs410;
        lhs(6,11)=-clhs138*clhs310 + clhs138*clhs400 + clhs346*clhs9 + clhs387*clhs9 + clhs40*clhs411 - clhs411 + clhs412;
        lhs(6,12)=DN(1,0)*clhs150 + DN(1,1)*clhs239 + DN(1,2)*clhs288 + clhs391*clhs9 + clhs414;
        lhs(6,13)=DN(1,0)*clhs163 + DN(1,1)*clhs246 + DN(1,2)*clhs291 + clhs359*clhs9 + clhs416;
        lhs(6,14)=DN(1,0)*clhs171 + DN(1,1)*clhs251 + DN(1,2)*clhs295 + clhs154*clhs400 + clhs353 + clhs354 + clhs355 + clhs356 + clhs418 + clhs419;
        lhs(6,15)=-clhs178*clhs310 + clhs178*clhs400 + clhs364*clhs9 + clhs398*clhs9 + clhs40*clhs420 - clhs420 + clhs421;
        lhs(7,0)=DN(0,0)*clhs422 + clhs15*clhs423;
        lhs(7,1)=DN(0,1)*clhs422 + clhs15*clhs424;
        lhs(7,2)=DN(0,2)*clhs422 + clhs15*clhs425;
        lhs(7,3)=clhs305 + clhs423*clhs45 + clhs424*clhs48 + clhs425*clhs52;
        lhs(7,4)=DN(1,0)*clhs426;
        lhs(7,5)=DN(1,1)*clhs426;
        lhs(7,6)=DN(1,2)*clhs426;
        lhs(7,7)=clhs300*clhs318 + clhs423*clhs92 + clhs424*clhs94 + clhs425*clhs97;
        lhs(7,8)=DN(2,0)*clhs422 + clhs114*clhs423;
        lhs(7,9)=DN(2,1)*clhs422 + clhs114*clhs424;
        lhs(7,10)=DN(2,2)*clhs422 + clhs114*clhs425;
        lhs(7,11)=clhs135*clhs423 + clhs136*clhs424 + clhs138*clhs425 + clhs428;
        lhs(7,12)=DN(3,0)*clhs422 + clhs154*clhs423;
        lhs(7,13)=DN(3,1)*clhs422 + clhs154*clhs424;
        lhs(7,14)=DN(3,2)*clhs422 + clhs154*clhs425;
        lhs(7,15)=clhs175*clhs423 + clhs176*clhs424 + clhs178*clhs425 + clhs429;
        lhs(8,0)=DN(2,0)*clhs1 + DN(2,1)*clhs3 + DN(2,2)*clhs5 + clhs102 + clhs103 + clhs15*clhs430 + clhs433;
        lhs(8,1)=DN(2,0)*clhs23 + DN(2,1)*clhs25 + DN(2,2)*clhs28 + clhs218 + clhs434*clhs6;
        lhs(8,2)=DN(2,0)*clhs32 + DN(2,1)*clhs34 + DN(2,2)*clhs36 + clhs275 + clhs435*clhs6;
        lhs(8,3)=clhs141 + clhs40*clhs436 + clhs430*clhs45 - clhs432*clhs45 - clhs436 + clhs437*clhs6 + clhs438*clhs6;
        lhs(8,4)=DN(2,0)*clhs61 + DN(2,1)*clhs63 + DN(2,2)*clhs65 + clhs332 + clhs333 + clhs430*clhs70 + clhs439;
        lhs(8,5)=DN(2,0)*clhs74 + DN(2,1)*clhs76 + DN(2,2)*clhs79 + clhs379 + clhs440*clhs6;
        lhs(8,6)=DN(2,0)*clhs83 + DN(2,1)*clhs85 + DN(2,2)*clhs87 + clhs405 + clhs441*clhs6;
        lhs(8,7)=clhs349 + clhs40*clhs442 + clhs430*clhs92 - clhs432*clhs92 - clhs442 + clhs443*clhs6 + clhs444*clhs6;
        lhs(8,8)=DN(2,0)*clhs105 + DN(2,1)*clhs107 + DN(2,2)*clhs109 + clhs114*clhs430 + clhs445*rho + clhs445*tau2 + clhs447;
        lhs(8,9)=DN(2,0)*clhs118 + DN(2,1)*clhs120 + DN(2,2)*clhs123 + clhs449 + clhs450*clhs6;
        lhs(8,10)=DN(2,0)*clhs127 + DN(2,1)*clhs129 + DN(2,2)*clhs131 + clhs452 + clhs453*clhs6;
        lhs(8,11)=clhs135*clhs430 - clhs135*clhs432 + clhs40*clhs454 + clhs42*clhs455 - clhs454 + clhs456*clhs6 + clhs457*clhs6 - clhs458*clhs55;
        lhs(8,12)=DN(2,0)*clhs146 + DN(2,1)*clhs148 + DN(2,2)*clhs150 + clhs154*clhs430 + clhs460 + clhs461 + clhs463 + clhs464 + clhs465 + clhs466;
        lhs(8,13)=DN(2,0)*clhs158 + DN(2,1)*clhs160 + DN(2,2)*clhs163 + clhs468 + clhs469*clhs6;
        lhs(8,14)=DN(2,0)*clhs167 + DN(2,1)*clhs169 + DN(2,2)*clhs171 + clhs471 + clhs472*clhs6;
        lhs(8,15)=clhs175*clhs430 - clhs175*clhs432 + clhs40*clhs473 - clhs473 + clhs474*clhs6 + clhs475*clhs6 + clhs478;
        lhs(9,0)=DN(2,0)*clhs3 + DN(2,1)*clhs182 + DN(2,2)*clhs183 + clhs117 + clhs479*clhs8;
        lhs(9,1)=DN(2,0)*clhs25 + DN(2,1)*clhs186 + DN(2,2)*clhs188 + clhs15*clhs480 + clhs223 + clhs224 + clhs433;
        lhs(9,2)=DN(2,0)*clhs34 + DN(2,1)*clhs192 + DN(2,2)*clhs194 + clhs278 + clhs435*clhs8;
        lhs(9,3)=clhs235 + clhs40*clhs481 - clhs432*clhs48 + clhs438*clhs8 + clhs48*clhs480 - clhs481 + clhs482*clhs8;
        lhs(9,4)=DN(2,0)*clhs63 + DN(2,1)*clhs200 + DN(2,2)*clhs201 + clhs340 + clhs483*clhs8;
        lhs(9,5)=DN(2,0)*clhs76 + DN(2,1)*clhs206 + DN(2,2)*clhs208 + clhs382 + clhs383 + clhs439 + clhs480*clhs70;
        lhs(9,6)=DN(2,0)*clhs85 + DN(2,1)*clhs211 + DN(2,2)*clhs213 + clhs407 + clhs441*clhs8;
        lhs(9,7)=clhs388 + clhs40*clhs484 - clhs432*clhs94 + clhs444*clhs8 + clhs480*clhs94 - clhs484 + clhs485*clhs8;
        lhs(9,8)=DN(2,0)*clhs107 + DN(2,1)*clhs219 + DN(2,2)*clhs220 + clhs449 + clhs486*clhs8;
        lhs(9,9)=DN(2,0)*clhs120 + DN(2,1)*clhs225 + DN(2,2)*clhs227 + clhs114*clhs480 + clhs447 + clhs487*rho + clhs487*tau2;
        lhs(9,10)=DN(2,0)*clhs129 + DN(2,1)*clhs230 + DN(2,2)*clhs232 + clhs453*clhs8 + clhs489;
        lhs(9,11)=-clhs136*clhs432 + clhs136*clhs480 - clhs197*clhs458 + clhs40*clhs490 + clhs455*clhs46 + clhs457*clhs8 - clhs490 + clhs491*clhs8;
        lhs(9,12)=DN(2,0)*clhs148 + DN(2,1)*clhs238 + DN(2,2)*clhs239 + clhs493 + clhs494*clhs8;
        lhs(9,13)=DN(2,0)*clhs160 + DN(2,1)*clhs244 + DN(2,2)*clhs246 + clhs154*clhs480 + clhs463 + clhs464 + clhs465 + clhs466 + clhs496 + clhs497;
        lhs(9,14)=DN(2,0)*clhs169 + DN(2,1)*clhs249 + DN(2,2)*clhs251 + clhs472*clhs8 + clhs499;
        lhs(9,15)=-clhs176*clhs432 + clhs176*clhs480 + clhs40*clhs500 + clhs475*clhs8 - clhs500 + clhs501*clhs8 + clhs502;
        lhs(10,0)=DN(2,0)*clhs5 + DN(2,1)*clhs183 + DN(2,2)*clhs255 + clhs126 + clhs479*clhs9;
        lhs(10,1)=DN(2,0)*clhs28 + DN(2,1)*clhs188 + DN(2,2)*clhs256 + clhs229 + clhs434*clhs9;
        lhs(10,2)=DN(2,0)*clhs36 + DN(2,1)*clhs194 + DN(2,2)*clhs258 + clhs15*clhs503 + clhs281 + clhs282 + clhs433;
        lhs(10,3)=clhs285 + clhs40*clhs504 - clhs432*clhs52 + clhs437*clhs9 + clhs482*clhs9 + clhs503*clhs52 - clhs504;
        lhs(10,4)=DN(2,0)*clhs65 + DN(2,1)*clhs201 + DN(2,2)*clhs264 + clhs343 + clhs483*clhs9;
        lhs(10,5)=DN(2,0)*clhs79 + DN(2,1)*clhs208 + DN(2,2)*clhs267 + clhs385 + clhs440*clhs9;
        lhs(10,6)=DN(2,0)*clhs87 + DN(2,1)*clhs213 + DN(2,2)*clhs271 + clhs409 + clhs410 + clhs439 + clhs503*clhs70;
        lhs(10,7)=clhs40*clhs505 + clhs412 - clhs432*clhs97 + clhs443*clhs9 + clhs485*clhs9 + clhs503*clhs97 - clhs505;
        lhs(10,8)=DN(2,0)*clhs109 + DN(2,1)*clhs220 + DN(2,2)*clhs276 + clhs452 + clhs486*clhs9;
        lhs(10,9)=DN(2,0)*clhs123 + DN(2,1)*clhs227 + DN(2,2)*clhs279 + clhs450*clhs9 + clhs489;
        lhs(10,10)=DN(2,0)*clhs131 + DN(2,1)*clhs232 + DN(2,2)*clhs283 + clhs114*clhs503 + clhs447 + clhs506*rho + clhs506*tau2;
        lhs(10,11)=-clhs138*clhs432 + clhs138*clhs503 - clhs261*clhs458 + clhs40*clhs507 + clhs455*clhs50 + clhs456*clhs9 + clhs491*clhs9 - clhs507;
        lhs(10,12)=DN(2,0)*clhs150 + DN(2,1)*clhs239 + DN(2,2)*clhs288 + clhs494*clhs9 + clhs509;
        lhs(10,13)=DN(2,0)*clhs163 + DN(2,1)*clhs246 + DN(2,2)*clhs291 + clhs469*clhs9 + clhs511;
        lhs(10,14)=DN(2,0)*clhs171 + DN(2,1)*clhs251 + DN(2,2)*clhs295 + clhs154*clhs503 + clhs463 + clhs464 + clhs465 + clhs466 + clhs513 + clhs514;
        lhs(10,15)=-clhs178*clhs432 + clhs178*clhs503 + clhs40*clhs515 + clhs474*clhs9 + clhs501*clhs9 - clhs515 + clhs516;
        lhs(11,0)=DN(0,0)*clhs517 + clhs15*clhs518;
        lhs(11,1)=DN(0,1)*clhs517 + clhs15*clhs519;
        lhs(11,2)=DN(0,2)*clhs517 + clhs15*clhs520;
        lhs(11,3)=clhs306 + clhs45*clhs518 + clhs48*clhs519 + clhs52*clhs520;
        lhs(11,4)=DN(1,0)*clhs517 + clhs518*clhs70;
        lhs(11,5)=DN(1,1)*clhs517 + clhs519*clhs70;
        lhs(11,6)=DN(1,2)*clhs517 + clhs520*clhs70;
        lhs(11,7)=clhs428 + clhs518*clhs92 + clhs519*clhs94 + clhs520*clhs97;
        lhs(11,8)=DN(2,0)*clhs521;
        lhs(11,9)=DN(2,1)*clhs521;
        lhs(11,10)=DN(2,2)*clhs521;
        lhs(11,11)=clhs135*clhs518 + clhs136*clhs519 + clhs138*clhs520 + clhs300*clhs446;
        lhs(11,12)=DN(3,0)*clhs517 + clhs154*clhs518;
        lhs(11,13)=DN(3,1)*clhs517 + clhs154*clhs519;
        lhs(11,14)=DN(3,2)*clhs517 + clhs154*clhs520;
        lhs(11,15)=clhs175*clhs518 + clhs176*clhs519 + clhs178*clhs520 + clhs522;
        lhs(12,0)=DN(3,0)*clhs1 + DN(3,1)*clhs3 + DN(3,2)*clhs5 + clhs143 + clhs144 + clhs15*clhs523 + clhs526;
        lhs(12,1)=DN(3,0)*clhs23 + DN(3,1)*clhs25 + DN(3,2)*clhs28 + clhs237 + clhs527*clhs6;
        lhs(12,2)=DN(3,0)*clhs32 + DN(3,1)*clhs34 + DN(3,2)*clhs36 + clhs287 + clhs528*clhs6;
        lhs(12,3)=clhs181 + clhs40*clhs529 + clhs45*clhs523 - clhs45*clhs525 - clhs529 + clhs530*clhs6 + clhs531*clhs6;
        lhs(12,4)=DN(3,0)*clhs61 + DN(3,1)*clhs63 + DN(3,2)*clhs65 + clhs351 + clhs352 + clhs523*clhs70 + clhs532;
        lhs(12,5)=DN(3,0)*clhs74 + DN(3,1)*clhs76 + DN(3,2)*clhs79 + clhs390 + clhs533*clhs6;
        lhs(12,6)=DN(3,0)*clhs83 + DN(3,1)*clhs85 + DN(3,2)*clhs87 + clhs414 + clhs534*clhs6;
        lhs(12,7)=clhs367 + clhs40*clhs535 + clhs523*clhs92 - clhs525*clhs92 - clhs535 + clhs536*clhs6 + clhs537*clhs6;
        lhs(12,8)=DN(3,0)*clhs105 + DN(3,1)*clhs107 + DN(3,2)*clhs109 + clhs114*clhs523 + clhs460 + clhs461 + clhs538;
        lhs(12,9)=DN(3,0)*clhs118 + DN(3,1)*clhs120 + DN(3,2)*clhs123 + clhs493 + clhs539*clhs6;
        lhs(12,10)=DN(3,0)*clhs127 + DN(3,1)*clhs129 + DN(3,2)*clhs131 + clhs509 + clhs540*clhs6;
        lhs(12,11)=clhs135*clhs523 - clhs135*clhs525 + clhs40*clhs541 + clhs478 - clhs541 + clhs542*clhs6 + clhs543*clhs6;
        lhs(12,12)=DN(3,0)*clhs146 + DN(3,1)*clhs148 + DN(3,2)*clhs150 + clhs154*clhs523 + clhs544*rho + clhs544*tau2 + clhs546;
        lhs(12,13)=DN(3,0)*clhs158 + DN(3,1)*clhs160 + DN(3,2)*clhs163 + clhs548 + clhs549*clhs6;
        lhs(12,14)=DN(3,0)*clhs167 + DN(3,1)*clhs169 + DN(3,2)*clhs171 + clhs551 + clhs552*clhs6;
        lhs(12,15)=clhs175*clhs523 - clhs175*clhs525 + clhs40*clhs553 + clhs42*clhs554 - clhs55*clhs557 - clhs553 + clhs555*clhs6 + clhs556*clhs6;
        lhs(13,0)=DN(3,0)*clhs3 + DN(3,1)*clhs182 + DN(3,2)*clhs183 + clhs157 + clhs558*clhs8;
        lhs(13,1)=DN(3,0)*clhs25 + DN(3,1)*clhs186 + DN(3,2)*clhs188 + clhs15*clhs559 + clhs242 + clhs243 + clhs526;
        lhs(13,2)=DN(3,0)*clhs34 + DN(3,1)*clhs192 + DN(3,2)*clhs194 + clhs290 + clhs528*clhs8;
        lhs(13,3)=clhs254 + clhs40*clhs560 - clhs48*clhs525 + clhs48*clhs559 + clhs531*clhs8 - clhs560 + clhs561*clhs8;
        lhs(13,4)=DN(3,0)*clhs63 + DN(3,1)*clhs200 + DN(3,2)*clhs201 + clhs358 + clhs562*clhs8;
        lhs(13,5)=DN(3,0)*clhs76 + DN(3,1)*clhs206 + DN(3,2)*clhs208 + clhs393 + clhs394 + clhs532 + clhs559*clhs70;
        lhs(13,6)=DN(3,0)*clhs85 + DN(3,1)*clhs211 + DN(3,2)*clhs213 + clhs416 + clhs534*clhs8;
        lhs(13,7)=clhs399 + clhs40*clhs563 - clhs525*clhs94 + clhs537*clhs8 + clhs559*clhs94 - clhs563 + clhs564*clhs8;
        lhs(13,8)=DN(3,0)*clhs107 + DN(3,1)*clhs219 + DN(3,2)*clhs220 + clhs468 + clhs565*clhs8;
        lhs(13,9)=DN(3,0)*clhs120 + DN(3,1)*clhs225 + DN(3,2)*clhs227 + clhs114*clhs559 + clhs496 + clhs497 + clhs538;
        lhs(13,10)=DN(3,0)*clhs129 + DN(3,1)*clhs230 + DN(3,2)*clhs232 + clhs511 + clhs540*clhs8;
        lhs(13,11)=-clhs136*clhs525 + clhs136*clhs559 + clhs40*clhs566 + clhs502 + clhs543*clhs8 - clhs566 + clhs567*clhs8;
        lhs(13,12)=DN(3,0)*clhs148 + DN(3,1)*clhs238 + DN(3,2)*clhs239 + clhs548 + clhs568*clhs8;
        lhs(13,13)=DN(3,0)*clhs160 + DN(3,1)*clhs244 + DN(3,2)*clhs246 + clhs154*clhs559 + clhs546 + clhs569*rho + clhs569*tau2;
        lhs(13,14)=DN(3,0)*clhs169 + DN(3,1)*clhs249 + DN(3,2)*clhs251 + clhs552*clhs8 + clhs571;
        lhs(13,15)=-clhs176*clhs525 + clhs176*clhs559 - clhs197*clhs557 + clhs40*clhs572 + clhs46*clhs554 + clhs556*clhs8 - clhs572 + clhs573*clhs8;
        lhs(14,0)=DN(3,0)*clhs5 + DN(3,1)*clhs183 + DN(3,2)*clhs255 + clhs166 + clhs558*clhs9;
        lhs(14,1)=DN(3,0)*clhs28 + DN(3,1)*clhs188 + DN(3,2)*clhs256 + clhs248 + clhs527*clhs9;
        lhs(14,2)=DN(3,0)*clhs36 + DN(3,1)*clhs194 + DN(3,2)*clhs258 + clhs15*clhs574 + clhs293 + clhs294 + clhs526;
        lhs(14,3)=clhs297 + clhs40*clhs575 - clhs52*clhs525 + clhs52*clhs574 + clhs530*clhs9 + clhs561*clhs9 - clhs575;
        lhs(14,4)=DN(3,0)*clhs65 + DN(3,1)*clhs201 + DN(3,2)*clhs264 + clhs361 + clhs562*clhs9;
        lhs(14,5)=DN(3,0)*clhs79 + DN(3,1)*clhs208 + DN(3,2)*clhs267 + clhs396 + clhs533*clhs9;
        lhs(14,6)=DN(3,0)*clhs87 + DN(3,1)*clhs213 + DN(3,2)*clhs271 + clhs418 + clhs419 + clhs532 + clhs574*clhs70;
        lhs(14,7)=clhs40*clhs576 + clhs421 - clhs525*clhs97 + clhs536*clhs9 + clhs564*clhs9 + clhs574*clhs97 - clhs576;
        lhs(14,8)=DN(3,0)*clhs109 + DN(3,1)*clhs220 + DN(3,2)*clhs276 + clhs471 + clhs565*clhs9;
        lhs(14,9)=DN(3,0)*clhs123 + DN(3,1)*clhs227 + DN(3,2)*clhs279 + clhs499 + clhs539*clhs9;
        lhs(14,10)=DN(3,0)*clhs131 + DN(3,1)*clhs232 + DN(3,2)*clhs283 + clhs114*clhs574 + clhs513 + clhs514 + clhs538;
        lhs(14,11)=-clhs138*clhs525 + clhs138*clhs574 + clhs40*clhs577 + clhs516 + clhs542*clhs9 + clhs567*clhs9 - clhs577;
        lhs(14,12)=DN(3,0)*clhs150 + DN(3,1)*clhs239 + DN(3,2)*clhs288 + clhs551 + clhs568*clhs9;
        lhs(14,13)=DN(3,0)*clhs163 + DN(3,1)*clhs246 + DN(3,2)*clhs291 + clhs549*clhs9 + clhs571;
        lhs(14,14)=DN(3,0)*clhs171 + DN(3,1)*clhs251 + DN(3,2)*clhs295 + clhs154*clhs574 + clhs546 + clhs578*rho + clhs578*tau2;
        lhs(14,15)=-clhs178*clhs525 + clhs178*clhs574 - clhs261*clhs557 + clhs40*clhs579 + clhs50*clhs554 + clhs555*clhs9 + clhs573*clhs9 - clhs579;
        lhs(15,0)=DN(0,0)*clhs580 + clhs15*clhs581;
        lhs(15,1)=DN(0,1)*clhs580 + clhs15*clhs582;
        lhs(15,2)=DN(0,2)*clhs580 + clhs15*clhs583;
        lhs(15,3)=clhs307 + clhs45*clhs581 + clhs48*clhs582 + clhs52*clhs583;
        lhs(15,4)=DN(1,0)*clhs580 + clhs581*clhs70;
        lhs(15,5)=DN(1,1)*clhs580 + clhs582*clhs70;
        lhs(15,6)=DN(1,2)*clhs580 + clhs583*clhs70;
        lhs(15,7)=clhs429 + clhs581*clhs92 + clhs582*clhs94 + clhs583*clhs97;
        lhs(15,8)=DN(2,0)*clhs580 + clhs114*clhs581;
        lhs(15,9)=DN(2,1)*clhs580 + clhs114*clhs582;
        lhs(15,10)=DN(2,2)*clhs580 + clhs114*clhs583;
        lhs(15,11)=clhs135*clhs581 + clhs136*clhs582 + clhs138*clhs583 + clhs522;
        lhs(15,12)=DN(3,0)*clhs584;
        lhs(15,13)=DN(3,1)*clhs584;
        lhs(15,14)=DN(3,2)*clhs584;
        lhs(15,15)=clhs175*clhs581 + clhs176*clhs582 + clhs178*clhs583 + clhs300*clhs545;

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
        const double clhs1 =             C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
        const double clhs2 =             C(0,2)*DN(0,0);
        const double clhs3 =             C(2,2)*DN(0,1) + clhs2;
        const double clhs4 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
        const double clhs5 =             DN(0,0)*clhs4*rho*tau1;
        const double clhs6 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
        const double clhs7 =             DN(0,0)*clhs4 + DN(0,1)*clhs6;
        const double clhs8 =             pow(c, -2);
        const double clhs9 =             N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]);
        const double clhs10 =             clhs8*clhs9;
        const double clhs11 =             N[0]*clhs10;
        const double clhs12 =             clhs11 + rho*(N[0]*bdf0 + clhs7);
        const double clhs13 =             pow(N[0], 2);
        const double clhs14 =             bdf0*rho;
        const double clhs15 =             N[0]*rho;
        const double clhs16 =             N[0]*clhs8*clhs9*tau1;
        const double clhs17 =             clhs10*clhs13 - clhs12*clhs16 + clhs13*clhs14 + clhs15*clhs7;
        const double clhs18 =             DN(0,0)*DN(0,1);
        const double clhs19 =             clhs18*rho + clhs18*tau2;
        const double clhs20 =             C(0,1)*DN(0,1) + clhs2;
        const double clhs21 =             C(1,2)*DN(0,1);
        const double clhs22 =             C(2,2)*DN(0,0) + clhs21;
        const double clhs23 =             clhs4*rho;
        const double clhs24 =             DN(0,0)*N[0];
        const double clhs25 =             1.0/rho;
        const double clhs26 =             bdf0*clhs25*clhs8*tau2;
        const double clhs27 =             bdf0*clhs13*clhs8;
        const double clhs28 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0);
        const double clhs29 =             N[0]*bdf0*clhs8;
        const double clhs30 =             clhs28*clhs29;
        const double clhs31 =             DN(0,0) + clhs30;
        const double clhs32 =             DN(0,1)*tau1;
        const double clhs33 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
        const double clhs34 =             clhs29*clhs33;
        const double clhs35 =             DN(0,1) + clhs34;
        const double clhs36 =             clhs32*clhs35;
        const double clhs37 =             bdf0*clhs13*clhs8*tau1;
        const double clhs38 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + clhs10*clhs28 - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + clhs4*(DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0)) + clhs6*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
        const double clhs39 =             DN(0,0)*DN(1,0);
        const double clhs40 =             clhs39*rho;
        const double clhs41 =             clhs39*tau2;
        const double clhs42 =             N[0]*bdf0*rho;
        const double clhs43 =             N[1]*clhs42;
        const double clhs44 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
        const double clhs45 =             C(0,2)*DN(1,0);
        const double clhs46 =             C(2,2)*DN(1,1) + clhs45;
        const double clhs47 =             DN(1,0)*clhs4 + DN(1,1)*clhs6;
        const double clhs48 =             clhs15*clhs47;
        const double clhs49 =             N[1]*clhs11;
        const double clhs50 =             N[1]*clhs10;
        const double clhs51 =             clhs50 + rho*(N[1]*bdf0 + clhs47);
        const double clhs52 =             -clhs16*clhs51;
        const double clhs53 =             DN(0,0)*DN(1,1);
        const double clhs54 =             clhs53*rho + clhs53*tau2;
        const double clhs55 =             C(0,1)*DN(1,1) + clhs45;
        const double clhs56 =             C(1,2)*DN(1,1);
        const double clhs57 =             C(2,2)*DN(1,0) + clhs56;
        const double clhs58 =             clhs32*clhs51;
        const double clhs59 =             DN(0,0)*N[1];
        const double clhs60 =             N[1]*bdf0*clhs8;
        const double clhs61 =             DN(1,0) + clhs28*clhs60;
        const double clhs62 =             DN(1,1) + clhs33*clhs60;
        const double clhs63 =             clhs32*clhs62;
        const double clhs64 =             N[0]*N[1]*bdf0*clhs8*tau1;
        const double clhs65 =             N[1]*clhs30 - clhs38*clhs64;
        const double clhs66 =             DN(0,0)*DN(2,0);
        const double clhs67 =             clhs66*rho;
        const double clhs68 =             clhs66*tau2;
        const double clhs69 =             N[2]*clhs42;
        const double clhs70 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
        const double clhs71 =             C(0,2)*DN(2,0);
        const double clhs72 =             C(2,2)*DN(2,1) + clhs71;
        const double clhs73 =             DN(2,0)*clhs4 + DN(2,1)*clhs6;
        const double clhs74 =             clhs15*clhs73;
        const double clhs75 =             N[2]*clhs11;
        const double clhs76 =             N[2]*clhs10 + rho*(N[2]*bdf0 + clhs73);
        const double clhs77 =             -clhs16*clhs76;
        const double clhs78 =             DN(0,0)*DN(2,1);
        const double clhs79 =             clhs78*rho + clhs78*tau2;
        const double clhs80 =             C(0,1)*DN(2,1) + clhs71;
        const double clhs81 =             C(1,2)*DN(2,1);
        const double clhs82 =             C(2,2)*DN(2,0) + clhs81;
        const double clhs83 =             clhs32*clhs76;
        const double clhs84 =             DN(0,0)*N[2];
        const double clhs85 =             N[2]*bdf0*clhs8;
        const double clhs86 =             DN(2,0) + clhs28*clhs85;
        const double clhs87 =             DN(2,1) + clhs33*clhs85;
        const double clhs88 =             clhs32*clhs87;
        const double clhs89 =             N[0]*N[2]*bdf0*clhs8*tau1;
        const double clhs90 =             N[2]*clhs30 - clhs38*clhs89;
        const double clhs91 =             C(0,1)*DN(0,0) + clhs21;
        const double clhs92 =             clhs6*rho;
        const double clhs93 =             pow(DN(0,1), 2);
        const double clhs94 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
        const double clhs95 =             DN(0,1)*clhs6*rho*tau1;
        const double clhs96 =             DN(0,1)*N[0];
        const double clhs97 =             DN(0,0)*tau1;
        const double clhs98 =             clhs31*clhs97;
        const double clhs99 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + clhs10*clhs33 - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + clhs4*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + clhs6*(DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1)));
        const double clhs100 =             DN(0,1)*DN(1,0);
        const double clhs101 =             clhs100*rho + clhs100*tau2;
        const double clhs102 =             C(0,1)*DN(1,0) + clhs56;
        const double clhs103 =             clhs51*clhs97;
        const double clhs104 =             DN(0,1)*DN(1,1);
        const double clhs105 =             clhs104*rho;
        const double clhs106 =             clhs104*tau2;
        const double clhs107 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
        const double clhs108 =             DN(0,1)*N[1];
        const double clhs109 =             clhs61*clhs97;
        const double clhs110 =             N[1]*clhs34 - clhs64*clhs99;
        const double clhs111 =             DN(0,1)*DN(2,0);
        const double clhs112 =             clhs111*rho + clhs111*tau2;
        const double clhs113 =             C(0,1)*DN(2,0) + clhs81;
        const double clhs114 =             clhs76*clhs97;
        const double clhs115 =             DN(0,1)*DN(2,1);
        const double clhs116 =             clhs115*rho;
        const double clhs117 =             clhs115*tau2;
        const double clhs118 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
        const double clhs119 =             DN(0,1)*N[2];
        const double clhs120 =             clhs86*clhs97;
        const double clhs121 =             N[2]*clhs34 - clhs89*clhs99;
        const double clhs122 =             2*N[0];
        const double clhs123 =             clhs12*tau1 + clhs122;
        const double clhs124 =             bdf0*clhs25*clhs8;
        const double clhs125 =             N[0]*bdf0*clhs25*clhs8;
        const double clhs126 =             N[1]*clhs125;
        const double clhs127 =             N[2]*clhs125;
        const double clhs128 =             DN(1,0)*clhs4*rho*tau1;
        const double clhs129 =             N[1]*rho;
        const double clhs130 =             N[1]*clhs8*clhs9*tau1;
        const double clhs131 =             -clhs12*clhs130 + clhs129*clhs7 + clhs43 + clhs49;
        const double clhs132 =             DN(1,1)*tau1;
        const double clhs133 =             clhs12*clhs132;
        const double clhs134 =             DN(1,0)*N[0];
        const double clhs135 =             clhs132*clhs35;
        const double clhs136 =             pow(DN(1,0), 2);
        const double clhs137 =             pow(N[1], 2);
        const double clhs138 =             clhs10*clhs137 + clhs129*clhs47 - clhs130*clhs51 + clhs137*clhs14;
        const double clhs139 =             DN(1,0)*DN(1,1);
        const double clhs140 =             clhs139*rho + clhs139*tau2;
        const double clhs141 =             DN(1,0)*N[1];
        const double clhs142 =             bdf0*clhs137*clhs8;
        const double clhs143 =             clhs132*clhs62;
        const double clhs144 =             bdf0*clhs137*clhs8*tau1;
        const double clhs145 =             DN(1,0)*DN(2,0);
        const double clhs146 =             clhs145*rho;
        const double clhs147 =             clhs145*tau2;
        const double clhs148 =             N[1]*N[2]*bdf0;
        const double clhs149 =             clhs148*rho;
        const double clhs150 =             clhs129*clhs73;
        const double clhs151 =             N[2]*clhs50;
        const double clhs152 =             -clhs130*clhs76;
        const double clhs153 =             DN(1,0)*DN(2,1);
        const double clhs154 =             clhs153*rho + clhs153*tau2;
        const double clhs155 =             clhs132*clhs76;
        const double clhs156 =             DN(1,0)*N[2];
        const double clhs157 =             clhs132*clhs87;
        const double clhs158 =             N[1]*N[2]*bdf0*clhs8;
        const double clhs159 =             N[1]*N[2]*bdf0*clhs8*tau1;
        const double clhs160 =             clhs158*clhs28 - clhs159*clhs38;
        const double clhs161 =             DN(1,0)*tau1;
        const double clhs162 =             clhs12*clhs161;
        const double clhs163 =             DN(1,1)*clhs6*rho*tau1;
        const double clhs164 =             DN(1,1)*N[0];
        const double clhs165 =             clhs161*clhs31;
        const double clhs166 =             pow(DN(1,1), 2);
        const double clhs167 =             DN(1,1)*N[1];
        const double clhs168 =             clhs161*clhs61;
        const double clhs169 =             DN(1,1)*DN(2,0);
        const double clhs170 =             clhs169*rho + clhs169*tau2;
        const double clhs171 =             clhs161*clhs76;
        const double clhs172 =             DN(1,1)*DN(2,1);
        const double clhs173 =             clhs172*rho;
        const double clhs174 =             clhs172*tau2;
        const double clhs175 =             DN(1,1)*N[2];
        const double clhs176 =             clhs161*clhs86;
        const double clhs177 =             clhs158*clhs33 - clhs159*clhs99;
        const double clhs178 =             2*N[1];
        const double clhs179 =             clhs178 + clhs51*tau1;
        const double clhs180 =             clhs148*clhs25*clhs8;
        const double clhs181 =             DN(2,0)*clhs4*rho*tau1;
        const double clhs182 =             N[2]*rho;
        const double clhs183 =             N[2]*clhs8*clhs9*tau1;
        const double clhs184 =             -clhs12*clhs183 + clhs182*clhs7 + clhs69 + clhs75;
        const double clhs185 =             DN(2,1)*tau1;
        const double clhs186 =             clhs12*clhs185;
        const double clhs187 =             DN(2,0)*N[0];
        const double clhs188 =             clhs185*clhs35;
        const double clhs189 =             clhs149 + clhs151 + clhs182*clhs47 - clhs183*clhs51;
        const double clhs190 =             clhs185*clhs51;
        const double clhs191 =             DN(2,0)*N[1];
        const double clhs192 =             clhs185*clhs62;
        const double clhs193 =             pow(DN(2,0), 2);
        const double clhs194 =             pow(N[2], 2);
        const double clhs195 =             clhs10*clhs194 + clhs14*clhs194 + clhs182*clhs73 - clhs183*clhs76;
        const double clhs196 =             DN(2,0)*DN(2,1);
        const double clhs197 =             clhs196*rho + clhs196*tau2;
        const double clhs198 =             DN(2,0)*N[2];
        const double clhs199 =             bdf0*clhs194*clhs8;
        const double clhs200 =             clhs185*clhs87;
        const double clhs201 =             bdf0*clhs194*clhs8*tau1;
        const double clhs202 =             DN(2,0)*tau1;
        const double clhs203 =             clhs12*clhs202;
        const double clhs204 =             DN(2,1)*clhs6*rho*tau1;
        const double clhs205 =             DN(2,1)*N[0];
        const double clhs206 =             clhs202*clhs31;
        const double clhs207 =             clhs202*clhs51;
        const double clhs208 =             DN(2,1)*N[1];
        const double clhs209 =             clhs202*clhs61;
        const double clhs210 =             pow(DN(2,1), 2);
        const double clhs211 =             DN(2,1)*N[2];
        const double clhs212 =             clhs202*clhs86;
        const double clhs213 =             2*N[2];
        const double clhs214 =             clhs213 + clhs76*tau1;

        lhs(0,0)=DN(0,0)*clhs1 + DN(0,1)*clhs3 + clhs0*rho + clhs0*tau2 + clhs12*clhs5 + clhs17;
        lhs(0,1)=DN(0,0)*clhs20 + DN(0,1)*clhs12*clhs23*tau1 + DN(0,1)*clhs22 + clhs19;
        lhs(0,2)=-clhs16*clhs31 + clhs23*clhs36 + clhs24*clhs26 - clhs24 + clhs27*clhs28 + clhs31*clhs5 - clhs37*clhs38;
        lhs(0,3)=DN(0,0)*clhs44 + DN(0,1)*clhs46 + clhs40 + clhs41 + clhs43 + clhs48 + clhs49 + clhs5*clhs51 + clhs52;
        lhs(0,4)=DN(0,0)*clhs55 + DN(0,1)*clhs57 + clhs23*clhs58 + clhs54;
        lhs(0,5)=-clhs16*clhs61 + clhs23*clhs63 + clhs26*clhs59 + clhs5*clhs61 - clhs59 + clhs65;
        lhs(0,6)=DN(0,0)*clhs70 + DN(0,1)*clhs72 + clhs5*clhs76 + clhs67 + clhs68 + clhs69 + clhs74 + clhs75 + clhs77;
        lhs(0,7)=DN(0,0)*clhs80 + DN(0,1)*clhs82 + clhs23*clhs83 + clhs79;
        lhs(0,8)=-clhs16*clhs86 + clhs23*clhs88 + clhs26*clhs84 + clhs5*clhs86 - clhs84 + clhs90;
        lhs(1,0)=DN(0,0)*clhs12*clhs92*tau1 + DN(0,0)*clhs3 + DN(0,1)*clhs91 + clhs19;
        lhs(1,1)=DN(0,0)*clhs22 + DN(0,1)*clhs94 + clhs12*clhs95 + clhs17 + clhs93*rho + clhs93*tau2;
        lhs(1,2)=-clhs16*clhs35 + clhs26*clhs96 + clhs27*clhs33 + clhs35*clhs95 - clhs37*clhs99 + clhs92*clhs98 - clhs96;
        lhs(1,3)=DN(0,0)*clhs46 + DN(0,1)*clhs102 + clhs101 + clhs103*clhs92;
        lhs(1,4)=DN(0,0)*clhs57 + DN(0,1)*clhs107 + clhs105 + clhs106 + clhs43 + clhs48 + clhs49 + clhs51*clhs95 + clhs52;
        lhs(1,5)=clhs108*clhs26 - clhs108 + clhs109*clhs92 + clhs110 - clhs16*clhs62 + clhs62*clhs95;
        lhs(1,6)=DN(0,0)*clhs72 + DN(0,1)*clhs113 + clhs112 + clhs114*clhs92;
        lhs(1,7)=DN(0,0)*clhs82 + DN(0,1)*clhs118 + clhs116 + clhs117 + clhs69 + clhs74 + clhs75 + clhs76*clhs95 + clhs77;
        lhs(1,8)=clhs119*clhs26 - clhs119 + clhs120*clhs92 + clhs121 - clhs16*clhs87 + clhs87*clhs95;
        lhs(2,0)=DN(0,0)*clhs123;
        lhs(2,1)=DN(0,1)*clhs123;
        lhs(2,2)=clhs124*clhs13 + clhs36 + clhs98;
        lhs(2,3)=DN(1,0)*clhs122 + clhs103;
        lhs(2,4)=DN(1,1)*clhs122 + clhs58;
        lhs(2,5)=clhs109 + clhs126 + clhs63;
        lhs(2,6)=DN(2,0)*clhs122 + clhs114;
        lhs(2,7)=DN(2,1)*clhs122 + clhs83;
        lhs(2,8)=clhs120 + clhs127 + clhs88;
        lhs(3,0)=DN(1,0)*clhs1 + DN(1,1)*clhs3 + clhs12*clhs128 + clhs131 + clhs40 + clhs41;
        lhs(3,1)=DN(1,0)*clhs20 + DN(1,1)*clhs22 + clhs101 + clhs133*clhs23;
        lhs(3,2)=clhs128*clhs31 - clhs130*clhs31 + clhs134*clhs26 - clhs134 + clhs135*clhs23 + clhs65;
        lhs(3,3)=DN(1,0)*clhs44 + DN(1,1)*clhs46 + clhs128*clhs51 + clhs136*rho + clhs136*tau2 + clhs138;
        lhs(3,4)=DN(1,0)*clhs55 + DN(1,1)*clhs23*clhs51*tau1 + DN(1,1)*clhs57 + clhs140;
        lhs(3,5)=clhs128*clhs61 - clhs130*clhs61 + clhs141*clhs26 - clhs141 + clhs142*clhs28 + clhs143*clhs23 - clhs144*clhs38;
        lhs(3,6)=DN(1,0)*clhs70 + DN(1,1)*clhs72 + clhs128*clhs76 + clhs146 + clhs147 + clhs149 + clhs150 + clhs151 + clhs152;
        lhs(3,7)=DN(1,0)*clhs80 + DN(1,1)*clhs82 + clhs154 + clhs155*clhs23;
        lhs(3,8)=clhs128*clhs86 - clhs130*clhs86 + clhs156*clhs26 - clhs156 + clhs157*clhs23 + clhs160;
        lhs(4,0)=DN(1,0)*clhs3 + DN(1,1)*clhs91 + clhs162*clhs92 + clhs54;
        lhs(4,1)=DN(1,0)*clhs22 + DN(1,1)*clhs94 + clhs105 + clhs106 + clhs12*clhs163 + clhs131;
        lhs(4,2)=clhs110 - clhs130*clhs35 + clhs163*clhs35 + clhs164*clhs26 - clhs164 + clhs165*clhs92;
        lhs(4,3)=DN(1,0)*clhs46 + DN(1,0)*clhs51*clhs92*tau1 + DN(1,1)*clhs102 + clhs140;
        lhs(4,4)=DN(1,0)*clhs57 + DN(1,1)*clhs107 + clhs138 + clhs163*clhs51 + clhs166*rho + clhs166*tau2;
        lhs(4,5)=-clhs130*clhs62 + clhs142*clhs33 - clhs144*clhs99 + clhs163*clhs62 + clhs167*clhs26 - clhs167 + clhs168*clhs92;
        lhs(4,6)=DN(1,0)*clhs72 + DN(1,1)*clhs113 + clhs170 + clhs171*clhs92;
        lhs(4,7)=DN(1,0)*clhs82 + DN(1,1)*clhs118 + clhs149 + clhs150 + clhs151 + clhs152 + clhs163*clhs76 + clhs173 + clhs174;
        lhs(4,8)=-clhs130*clhs87 + clhs163*clhs87 + clhs175*clhs26 - clhs175 + clhs176*clhs92 + clhs177;
        lhs(5,0)=DN(0,0)*clhs178 + clhs162;
        lhs(5,1)=DN(0,1)*clhs178 + clhs133;
        lhs(5,2)=clhs126 + clhs135 + clhs165;
        lhs(5,3)=DN(1,0)*clhs179;
        lhs(5,4)=DN(1,1)*clhs179;
        lhs(5,5)=clhs124*clhs137 + clhs143 + clhs168;
        lhs(5,6)=DN(2,0)*clhs178 + clhs171;
        lhs(5,7)=DN(2,1)*clhs178 + clhs155;
        lhs(5,8)=clhs157 + clhs176 + clhs180;
        lhs(6,0)=DN(2,0)*clhs1 + DN(2,1)*clhs3 + clhs12*clhs181 + clhs184 + clhs67 + clhs68;
        lhs(6,1)=DN(2,0)*clhs20 + DN(2,1)*clhs22 + clhs112 + clhs186*clhs23;
        lhs(6,2)=clhs181*clhs31 - clhs183*clhs31 + clhs187*clhs26 - clhs187 + clhs188*clhs23 + clhs90;
        lhs(6,3)=DN(2,0)*clhs44 + DN(2,1)*clhs46 + clhs146 + clhs147 + clhs181*clhs51 + clhs189;
        lhs(6,4)=DN(2,0)*clhs55 + DN(2,1)*clhs57 + clhs170 + clhs190*clhs23;
        lhs(6,5)=clhs160 + clhs181*clhs61 - clhs183*clhs61 + clhs191*clhs26 - clhs191 + clhs192*clhs23;
        lhs(6,6)=DN(2,0)*clhs70 + DN(2,1)*clhs72 + clhs181*clhs76 + clhs193*rho + clhs193*tau2 + clhs195;
        lhs(6,7)=DN(2,0)*clhs80 + DN(2,1)*clhs23*clhs76*tau1 + DN(2,1)*clhs82 + clhs197;
        lhs(6,8)=clhs181*clhs86 - clhs183*clhs86 + clhs198*clhs26 - clhs198 + clhs199*clhs28 + clhs200*clhs23 - clhs201*clhs38;
        lhs(7,0)=DN(2,0)*clhs3 + DN(2,1)*clhs91 + clhs203*clhs92 + clhs79;
        lhs(7,1)=DN(2,0)*clhs22 + DN(2,1)*clhs94 + clhs116 + clhs117 + clhs12*clhs204 + clhs184;
        lhs(7,2)=clhs121 - clhs183*clhs35 + clhs204*clhs35 + clhs205*clhs26 - clhs205 + clhs206*clhs92;
        lhs(7,3)=DN(2,0)*clhs46 + DN(2,1)*clhs102 + clhs154 + clhs207*clhs92;
        lhs(7,4)=DN(2,0)*clhs57 + DN(2,1)*clhs107 + clhs173 + clhs174 + clhs189 + clhs204*clhs51;
        lhs(7,5)=clhs177 - clhs183*clhs62 + clhs204*clhs62 + clhs208*clhs26 - clhs208 + clhs209*clhs92;
        lhs(7,6)=DN(2,0)*clhs72 + DN(2,0)*clhs76*clhs92*tau1 + DN(2,1)*clhs113 + clhs197;
        lhs(7,7)=DN(2,0)*clhs82 + DN(2,1)*clhs118 + clhs195 + clhs204*clhs76 + clhs210*rho + clhs210*tau2;
        lhs(7,8)=-clhs183*clhs87 + clhs199*clhs33 - clhs201*clhs99 + clhs204*clhs87 + clhs211*clhs26 - clhs211 + clhs212*clhs92;
        lhs(8,0)=DN(0,0)*clhs213 + clhs203;
        lhs(8,1)=DN(0,1)*clhs213 + clhs186;
        lhs(8,2)=clhs127 + clhs188 + clhs206;
        lhs(8,3)=DN(1,0)*clhs213 + clhs207;
        lhs(8,4)=DN(1,1)*clhs213 + clhs190;
        lhs(8,5)=clhs180 + clhs192 + clhs209;
        lhs(8,6)=DN(2,0)*clhs214;
        lhs(8,7)=DN(2,1)*clhs214;
        lhs(8,8)=clhs124*clhs194 + clhs200 + clhs212;

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
    // const array_1d<double,strain_size>& stress = data.stress;

    // Get constitutive matrix
    const Matrix& C = data.C;

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
        const double crhs12 =             crhs11*rho;
        const double crhs13 =             N[0]*rho;
        const double crhs14 =             N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0));
        const double crhs15 =             pow(c, -2);
        const double crhs16 =             N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]) + N[3]*(bdf0*p[3] + bdf1*pn[3] + bdf2*pnn[3]);
        const double crhs17 =             crhs15*crhs16;
        const double crhs18 =             crhs17*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0));
        const double crhs19 =             crhs3 + crhs5 + crhs7 + crhs9;
        const double crhs20 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
        const double crhs21 =             DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0);
        const double crhs22 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
        const double crhs23 =             DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0);
        const double crhs24 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
        const double crhs25 =             crhs19*crhs20 + crhs21*crhs22 + crhs23*crhs24;
        const double crhs26 =             crhs15*crhs16/rho;
        const double crhs27 =             tau2*(crhs11 + crhs26);
        const double crhs28 =             crhs10 + crhs4 + crhs6 + crhs8;
        const double crhs29 =             DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1);
        const double crhs30 =             crhs21 + crhs29;
        const double crhs31 =             DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2);
        const double crhs32 =             DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1);
        const double crhs33 =             crhs31 + crhs32;
        const double crhs34 =             DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2);
        const double crhs35 =             crhs23 + crhs34;
        const double crhs36 =             C(0,0)*crhs19 + C(0,1)*crhs28 + C(0,2)*crhs2 + C(0,3)*crhs30 + C(0,4)*crhs33 + C(0,5)*crhs35;
        const double crhs37 =             C(0,3)*crhs19 + C(1,3)*crhs28 + C(2,3)*crhs2 + C(3,3)*crhs30 + C(3,4)*crhs33 + C(3,5)*crhs35;
        const double crhs38 =             C(0,5)*crhs19 + C(1,5)*crhs28 + C(2,5)*crhs2 + C(3,5)*crhs30 + C(4,5)*crhs33 + C(5,5)*crhs35;
        const double crhs39 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + crhs18 + rho*(crhs14 + crhs25);
        const double crhs40 =             DN(0,0)*crhs39*rho*tau1;
        const double crhs41 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
        const double crhs42 =             crhs17*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1));
        const double crhs43 =             N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1));
        const double crhs44 =             crhs20*crhs29 + crhs22*crhs28 + crhs24*crhs32;
        const double crhs45 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs41 + crhs42 + rho*(crhs43 + crhs44);
        const double crhs46 =             DN(0,1)*crhs45*rho*tau1;
        const double crhs47 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
        const double crhs48 =             crhs17*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2));
        const double crhs49 =             N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2));
        const double crhs50 =             crhs2*crhs24 + crhs20*crhs34 + crhs22*crhs31;
        const double crhs51 =             DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs47 + crhs48 + rho*(crhs49 + crhs50);
        const double crhs52 =             DN(0,2)*crhs51*rho*tau1;
        const double crhs53 =             N[0]*crhs15*crhs16;
        const double crhs54 =             crhs39*tau1;
        const double crhs55 =             C(0,1)*crhs19 + C(1,1)*crhs28 + C(1,2)*crhs2 + C(1,3)*crhs30 + C(1,4)*crhs33 + C(1,5)*crhs35;
        const double crhs56 =             C(0,4)*crhs19 + C(1,4)*crhs28 + C(2,4)*crhs2 + C(3,4)*crhs30 + C(4,4)*crhs33 + C(4,5)*crhs35;
        const double crhs57 =             crhs45*tau1;
        const double crhs58 =             C(0,2)*crhs19 + C(1,2)*crhs28 + C(2,2)*crhs2 + C(2,3)*crhs30 + C(2,4)*crhs33 + C(2,5)*crhs35;
        const double crhs59 =             crhs51*tau1;
        const double crhs60 =             2*DN(0,0)*v(0,0) + 2*DN(0,1)*v(0,1) + 2*DN(0,2)*v(0,2) + 2*DN(1,0)*v(1,0) + 2*DN(1,1)*v(1,1) + 2*DN(1,2)*v(1,2) + 2*DN(2,0)*v(2,0) + 2*DN(2,1)*v(2,1) + 2*DN(2,2)*v(2,2) + 2*DN(3,0)*v(3,0) + 2*DN(3,1)*v(3,1) + 2*DN(3,2)*v(3,2);
        const double crhs61 =             N[1]*rho;
        const double crhs62 =             DN(1,0)*crhs39*rho*tau1;
        const double crhs63 =             DN(1,1)*crhs45*rho*tau1;
        const double crhs64 =             DN(1,2)*crhs51*rho*tau1;
        const double crhs65 =             N[1]*crhs15*crhs16;
        const double crhs66 =             N[2]*rho;
        const double crhs67 =             DN(2,0)*crhs39*rho*tau1;
        const double crhs68 =             DN(2,1)*crhs45*rho*tau1;
        const double crhs69 =             DN(2,2)*crhs51*rho*tau1;
        const double crhs70 =             N[2]*crhs15*crhs16;
        const double crhs71 =             N[3]*rho;
        const double crhs72 =             DN(3,0)*crhs39*rho*tau1;
        const double crhs73 =             DN(3,1)*crhs45*rho*tau1;
        const double crhs74 =             DN(3,2)*crhs51*rho*tau1;
        const double crhs75 =             N[3]*crhs15*crhs16;

        rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs12 - DN(0,0)*crhs27 - DN(0,0)*crhs36 - DN(0,1)*crhs37 - DN(0,2)*crhs38 + N[0]*crhs1 - N[0]*crhs18 - crhs13*crhs14 - crhs13*crhs25 - crhs20*crhs40 - crhs20*crhs46 - crhs20*crhs52 + crhs53*crhs54;
        rhs[1]=-DN(0,0)*crhs37 + DN(0,1)*crhs0 - DN(0,1)*crhs12 - DN(0,1)*crhs27 - DN(0,1)*crhs55 - DN(0,2)*crhs56 + N[0]*crhs41 - N[0]*crhs42 - crhs13*crhs43 - crhs13*crhs44 - crhs22*crhs40 - crhs22*crhs46 - crhs22*crhs52 + crhs53*crhs57;
        rhs[2]=-DN(0,0)*crhs38 - DN(0,1)*crhs56 + DN(0,2)*crhs0 - DN(0,2)*crhs12 - DN(0,2)*crhs27 - DN(0,2)*crhs58 + N[0]*crhs47 - N[0]*crhs48 - crhs13*crhs49 - crhs13*crhs50 - crhs24*crhs40 - crhs24*crhs46 - crhs24*crhs52 + crhs53*crhs59;
        rhs[3]=-DN(0,0)*crhs54 - DN(0,1)*crhs57 - DN(0,2)*crhs59 - N[0]*crhs26 - N[0]*crhs60;
        rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs12 - DN(1,0)*crhs27 - DN(1,0)*crhs36 - DN(1,1)*crhs37 - DN(1,2)*crhs38 + N[1]*crhs1 - N[1]*crhs18 - crhs14*crhs61 - crhs20*crhs62 - crhs20*crhs63 - crhs20*crhs64 - crhs25*crhs61 + crhs54*crhs65;
        rhs[5]=-DN(1,0)*crhs37 + DN(1,1)*crhs0 - DN(1,1)*crhs12 - DN(1,1)*crhs27 - DN(1,1)*crhs55 - DN(1,2)*crhs56 + N[1]*crhs41 - N[1]*crhs42 - crhs22*crhs62 - crhs22*crhs63 - crhs22*crhs64 - crhs43*crhs61 - crhs44*crhs61 + crhs57*crhs65;
        rhs[6]=-DN(1,0)*crhs38 - DN(1,1)*crhs56 + DN(1,2)*crhs0 - DN(1,2)*crhs12 - DN(1,2)*crhs27 - DN(1,2)*crhs58 + N[1]*crhs47 - N[1]*crhs48 - crhs24*crhs62 - crhs24*crhs63 - crhs24*crhs64 - crhs49*crhs61 - crhs50*crhs61 + crhs59*crhs65;
        rhs[7]=-DN(1,0)*crhs54 - DN(1,1)*crhs57 - DN(1,2)*crhs59 - N[1]*crhs26 - N[1]*crhs60;
        rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs12 - DN(2,0)*crhs27 - DN(2,0)*crhs36 - DN(2,1)*crhs37 - DN(2,2)*crhs38 + N[2]*crhs1 - N[2]*crhs18 - crhs14*crhs66 - crhs20*crhs67 - crhs20*crhs68 - crhs20*crhs69 - crhs25*crhs66 + crhs54*crhs70;
        rhs[9]=-DN(2,0)*crhs37 + DN(2,1)*crhs0 - DN(2,1)*crhs12 - DN(2,1)*crhs27 - DN(2,1)*crhs55 - DN(2,2)*crhs56 + N[2]*crhs41 - N[2]*crhs42 - crhs22*crhs67 - crhs22*crhs68 - crhs22*crhs69 - crhs43*crhs66 - crhs44*crhs66 + crhs57*crhs70;
        rhs[10]=-DN(2,0)*crhs38 - DN(2,1)*crhs56 + DN(2,2)*crhs0 - DN(2,2)*crhs12 - DN(2,2)*crhs27 - DN(2,2)*crhs58 + N[2]*crhs47 - N[2]*crhs48 - crhs24*crhs67 - crhs24*crhs68 - crhs24*crhs69 - crhs49*crhs66 - crhs50*crhs66 + crhs59*crhs70;
        rhs[11]=-DN(2,0)*crhs54 - DN(2,1)*crhs57 - DN(2,2)*crhs59 - N[2]*crhs26 - N[2]*crhs60;
        rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs12 - DN(3,0)*crhs27 - DN(3,0)*crhs36 - DN(3,1)*crhs37 - DN(3,2)*crhs38 + N[3]*crhs1 - N[3]*crhs18 - crhs14*crhs71 - crhs20*crhs72 - crhs20*crhs73 - crhs20*crhs74 - crhs25*crhs71 + crhs54*crhs75;
        rhs[13]=-DN(3,0)*crhs37 + DN(3,1)*crhs0 - DN(3,1)*crhs12 - DN(3,1)*crhs27 - DN(3,1)*crhs55 - DN(3,2)*crhs56 + N[3]*crhs41 - N[3]*crhs42 - crhs22*crhs72 - crhs22*crhs73 - crhs22*crhs74 - crhs43*crhs71 - crhs44*crhs71 + crhs57*crhs75;
        rhs[14]=-DN(3,0)*crhs38 - DN(3,1)*crhs56 + DN(3,2)*crhs0 - DN(3,2)*crhs12 - DN(3,2)*crhs27 - DN(3,2)*crhs58 + N[3]*crhs47 - N[3]*crhs48 - crhs24*crhs72 - crhs24*crhs73 - crhs24*crhs74 - crhs49*crhs71 - crhs50*crhs71 + crhs59*crhs75;
        rhs[15]=-DN(3,0)*crhs54 - DN(3,1)*crhs57 - DN(3,2)*crhs59 - N[3]*crhs26 - N[3]*crhs60;

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
    // const array_1d<double,strain_size>& stress = data.stress;

    // Get constitutive matrix
    const Matrix& C = data.C;

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
        const double crhs5 =             crhs4*rho;
        const double crhs6 =             N[0]*rho;
        const double crhs7 =             N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0));
        const double crhs8 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
        const double crhs9 =             DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0);
        const double crhs10 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
        const double crhs11 =             crhs10*crhs9 + crhs2*crhs8;
        const double crhs12 =             DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1);
        const double crhs13 =             crhs12 + crhs9;
        const double crhs14 =             C(0,0)*crhs2 + C(0,1)*crhs3 + C(0,2)*crhs13;
        const double crhs15 =             C(0,2)*crhs2 + C(1,2)*crhs3 + C(2,2)*crhs13;
        const double crhs16 =             pow(c, -2);
        const double crhs17 =             N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]);
        const double crhs18 =             crhs16*crhs17;
        const double crhs19 =             crhs18*(N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0));
        const double crhs20 =             crhs16*crhs17/rho;
        const double crhs21 =             tau2*(crhs20 + crhs4);
        const double crhs22 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + crhs19 + rho*(crhs11 + crhs7);
        const double crhs23 =             DN(0,0)*crhs22*rho*tau1;
        const double crhs24 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
        const double crhs25 =             crhs18*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1));
        const double crhs26 =             N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1));
        const double crhs27 =             crhs10*crhs3 + crhs12*crhs8;
        const double crhs28 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs24 + crhs25 + rho*(crhs26 + crhs27);
        const double crhs29 =             DN(0,1)*crhs28*rho*tau1;
        const double crhs30 =             N[0]*crhs16*crhs17;
        const double crhs31 =             crhs22*tau1;
        const double crhs32 =             C(0,1)*crhs2 + C(1,1)*crhs3 + C(1,2)*crhs13;
        const double crhs33 =             crhs28*tau1;
        const double crhs34 =             2*DN(0,0)*v(0,0) + 2*DN(0,1)*v(0,1) + 2*DN(1,0)*v(1,0) + 2*DN(1,1)*v(1,1) + 2*DN(2,0)*v(2,0) + 2*DN(2,1)*v(2,1);
        const double crhs35 =             N[1]*rho;
        const double crhs36 =             DN(1,0)*crhs22*rho*tau1;
        const double crhs37 =             DN(1,1)*crhs28*rho*tau1;
        const double crhs38 =             N[1]*crhs16*crhs17;
        const double crhs39 =             N[2]*rho;
        const double crhs40 =             DN(2,0)*crhs22*rho*tau1;
        const double crhs41 =             DN(2,1)*crhs28*rho*tau1;
        const double crhs42 =             N[2]*crhs16*crhs17;

        rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs14 - DN(0,0)*crhs21 - DN(0,0)*crhs5 - DN(0,1)*crhs15 + N[0]*crhs1 - N[0]*crhs19 - crhs11*crhs6 - crhs23*crhs8 - crhs29*crhs8 + crhs30*crhs31 - crhs6*crhs7;
        rhs[1]=-DN(0,0)*crhs15 + DN(0,1)*crhs0 - DN(0,1)*crhs21 - DN(0,1)*crhs32 - DN(0,1)*crhs5 + N[0]*crhs24 - N[0]*crhs25 - crhs10*crhs23 - crhs10*crhs29 - crhs26*crhs6 - crhs27*crhs6 + crhs30*crhs33;
        rhs[2]=-DN(0,0)*crhs31 - DN(0,1)*crhs33 - N[0]*crhs20 - N[0]*crhs34;
        rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs14 - DN(1,0)*crhs21 - DN(1,0)*crhs5 - DN(1,1)*crhs15 + N[1]*crhs1 - N[1]*crhs19 - crhs11*crhs35 + crhs31*crhs38 - crhs35*crhs7 - crhs36*crhs8 - crhs37*crhs8;
        rhs[4]=-DN(1,0)*crhs15 + DN(1,1)*crhs0 - DN(1,1)*crhs21 - DN(1,1)*crhs32 - DN(1,1)*crhs5 + N[1]*crhs24 - N[1]*crhs25 - crhs10*crhs36 - crhs10*crhs37 - crhs26*crhs35 - crhs27*crhs35 + crhs33*crhs38;
        rhs[5]=-DN(1,0)*crhs31 - DN(1,1)*crhs33 - N[1]*crhs20 - N[1]*crhs34;
        rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs14 - DN(2,0)*crhs21 - DN(2,0)*crhs5 - DN(2,1)*crhs15 + N[2]*crhs1 - N[2]*crhs19 - crhs11*crhs39 + crhs31*crhs42 - crhs39*crhs7 - crhs40*crhs8 - crhs41*crhs8;
        rhs[7]=-DN(2,0)*crhs15 + DN(2,1)*crhs0 - DN(2,1)*crhs21 - DN(2,1)*crhs32 - DN(2,1)*crhs5 + N[2]*crhs24 - N[2]*crhs25 - crhs10*crhs40 - crhs10*crhs41 - crhs26*crhs39 - crhs27*crhs39 + crhs33*crhs42;
        rhs[8]=-DN(2,0)*crhs31 - DN(2,1)*crhs33 - N[2]*crhs20 - N[2]*crhs34;

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

    // Get constitutive matrix
    const Matrix& C = data.C;

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
