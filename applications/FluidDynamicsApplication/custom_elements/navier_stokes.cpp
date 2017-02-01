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
        const double clhs42 =             pow(c, -2);
        const double clhs43 =             1.0/rho;
        const double clhs44 =             N[0]*bdf0*clhs42*clhs43;
        const double clhs45 =             rho*tau1;
        const double clhs46 =             clhs14*clhs45;
        const double clhs47 =             -N[0] + clhs44*tau2 + clhs46;
        const double clhs48 =             N[0]*bdf0*rho;
        const double clhs49 =             N[1]*clhs48;
        const double clhs50 =             DN(1,0)*clhs24 + clhs49;
        const double clhs51 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
        const double clhs52 =             C(0,3)*DN(1,0);
        const double clhs53 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs52;
        const double clhs54 =             C(0,5)*DN(1,0);
        const double clhs55 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs54;
        const double clhs56 =             DN(1,0)*clhs11 + DN(1,1)*clhs12 + DN(1,2)*clhs13;
        const double clhs57 =             N[1]*clhs10 + clhs56;
        const double clhs58 =             N[1]*bdf0;
        const double clhs59 =             clhs57 + clhs58;
        const double clhs60 =             DN(0,0)*N[1]*rho*tau1;
        const double clhs61 =             DN(1,1)*clhs24;
        const double clhs62 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs52;
        const double clhs63 =             C(1,3)*DN(1,1);
        const double clhs64 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs63;
        const double clhs65 =             C(3,5)*DN(1,0);
        const double clhs66 =             C(4,5)*DN(1,2);
        const double clhs67 =             C(1,5)*DN(1,1) + clhs65 + clhs66;
        const double clhs68 =             N[0]*N[1]*rho;
        const double clhs69 =             clhs21*clhs68;
        const double clhs70 =             N[1]*clhs14*clhs16*tau1;
        const double clhs71 =             DN(0,1)*N[1]*rho*tau1;
        const double clhs72 =             DN(1,2)*clhs24;
        const double clhs73 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs54;
        const double clhs74 =             C(3,4)*DN(1,1);
        const double clhs75 =             C(2,3)*DN(1,2) + clhs65 + clhs74;
        const double clhs76 =             C(2,5)*DN(1,2);
        const double clhs77 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs76;
        const double clhs78 =             clhs22*clhs68;
        const double clhs79 =             DN(0,2)*N[1]*rho*tau1;
        const double clhs80 =             DN(0,0)*N[1];
        const double clhs81 =             bdf0*clhs42*clhs43*tau2;
        const double clhs82 =             DN(1,0)*rho*tau1;
        const double clhs83 =             N[2]*clhs48;
        const double clhs84 =             DN(2,0)*clhs24 + clhs83;
        const double clhs85 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
        const double clhs86 =             C(0,3)*DN(2,0);
        const double clhs87 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs86;
        const double clhs88 =             C(0,5)*DN(2,0);
        const double clhs89 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs88;
        const double clhs90 =             DN(2,0)*clhs11 + DN(2,1)*clhs12 + DN(2,2)*clhs13;
        const double clhs91 =             N[2]*clhs10 + clhs90;
        const double clhs92 =             N[2]*bdf0;
        const double clhs93 =             clhs91 + clhs92;
        const double clhs94 =             DN(0,0)*N[2]*rho*tau1;
        const double clhs95 =             DN(2,1)*clhs24;
        const double clhs96 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs86;
        const double clhs97 =             C(1,3)*DN(2,1);
        const double clhs98 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs97;
        const double clhs99 =             C(3,5)*DN(2,0);
        const double clhs100 =             C(4,5)*DN(2,2);
        const double clhs101 =             C(1,5)*DN(2,1) + clhs100 + clhs99;
        const double clhs102 =             N[0]*N[2]*rho;
        const double clhs103 =             clhs102*clhs21;
        const double clhs104 =             N[2]*clhs14*clhs16*tau1;
        const double clhs105 =             DN(0,1)*N[2]*rho*tau1;
        const double clhs106 =             DN(2,2)*clhs24;
        const double clhs107 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs88;
        const double clhs108 =             C(3,4)*DN(2,1);
        const double clhs109 =             C(2,3)*DN(2,2) + clhs108 + clhs99;
        const double clhs110 =             C(2,5)*DN(2,2);
        const double clhs111 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs110;
        const double clhs112 =             clhs102*clhs22;
        const double clhs113 =             DN(0,2)*N[2]*rho*tau1;
        const double clhs114 =             DN(0,0)*N[2];
        const double clhs115 =             DN(2,0)*rho*tau1;
        const double clhs116 =             N[3]*clhs48;
        const double clhs117 =             DN(3,0)*clhs24 + clhs116;
        const double clhs118 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
        const double clhs119 =             C(0,3)*DN(3,0);
        const double clhs120 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs119;
        const double clhs121 =             C(0,5)*DN(3,0);
        const double clhs122 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs121;
        const double clhs123 =             DN(3,0)*clhs11 + DN(3,1)*clhs12 + DN(3,2)*clhs13;
        const double clhs124 =             N[3]*clhs10 + clhs123;
        const double clhs125 =             N[3]*bdf0;
        const double clhs126 =             clhs124 + clhs125;
        const double clhs127 =             DN(0,0)*N[3]*rho*tau1;
        const double clhs128 =             DN(3,1)*clhs24;
        const double clhs129 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs119;
        const double clhs130 =             C(1,3)*DN(3,1);
        const double clhs131 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs130;
        const double clhs132 =             C(3,5)*DN(3,0);
        const double clhs133 =             C(4,5)*DN(3,2);
        const double clhs134 =             C(1,5)*DN(3,1) + clhs132 + clhs133;
        const double clhs135 =             N[0]*N[3]*rho;
        const double clhs136 =             clhs135*clhs21;
        const double clhs137 =             N[3]*clhs14*clhs16*tau1;
        const double clhs138 =             DN(0,1)*N[3]*rho*tau1;
        const double clhs139 =             DN(3,2)*clhs24;
        const double clhs140 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs121;
        const double clhs141 =             C(3,4)*DN(3,1);
        const double clhs142 =             C(2,3)*DN(3,2) + clhs132 + clhs141;
        const double clhs143 =             C(2,5)*DN(3,2);
        const double clhs144 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs143;
        const double clhs145 =             clhs135*clhs22;
        const double clhs146 =             DN(0,2)*N[3]*rho*tau1;
        const double clhs147 =             DN(0,0)*N[3];
        const double clhs148 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs27;
        const double clhs149 =             C(0,4)*DN(0,0) + clhs30 + clhs37;
        const double clhs150 =             DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1);
        const double clhs151 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1);
        const double clhs152 =             DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1);
        const double clhs153 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)) + clhs11*clhs150 + clhs12*clhs151 + clhs13*clhs152);
        const double clhs154 =             pow(DN(0,1), 2);
        const double clhs155 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
        const double clhs156 =             C(1,4)*DN(0,1);
        const double clhs157 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs156;
        const double clhs158 =             N[0]*clhs151 + clhs14;
        const double clhs159 =             clhs158 + clhs18;
        const double clhs160 =             DN(0,1)*tau2;
        const double clhs161 =             DN(0,2)*clhs160;
        const double clhs162 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs156;
        const double clhs163 =             C(2,4)*DN(0,2);
        const double clhs164 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs163;
        const double clhs165 =             DN(1,0)*clhs160;
        const double clhs166 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs63;
        const double clhs167 =             C(0,4)*DN(1,0) + clhs66 + clhs74;
        const double clhs168 =             clhs150*clhs68;
        const double clhs169 =             DN(1,1)*clhs160 + clhs49;
        const double clhs170 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
        const double clhs171 =             C(1,4)*DN(1,1);
        const double clhs172 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs171;
        const double clhs173 =             N[1]*clhs151 + clhs56;
        const double clhs174 =             clhs173 + clhs58;
        const double clhs175 =             DN(1,2)*clhs160;
        const double clhs176 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs171;
        const double clhs177 =             C(2,4)*DN(1,2);
        const double clhs178 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs177;
        const double clhs179 =             clhs152*clhs68;
        const double clhs180 =             DN(0,1)*N[1];
        const double clhs181 =             DN(1,1)*rho*tau1;
        const double clhs182 =             DN(2,0)*clhs160;
        const double clhs183 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs97;
        const double clhs184 =             C(0,4)*DN(2,0) + clhs100 + clhs108;
        const double clhs185 =             clhs102*clhs150;
        const double clhs186 =             DN(2,1)*clhs160 + clhs83;
        const double clhs187 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
        const double clhs188 =             C(1,4)*DN(2,1);
        const double clhs189 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs188;
        const double clhs190 =             N[2]*clhs151 + clhs90;
        const double clhs191 =             clhs190 + clhs92;
        const double clhs192 =             DN(2,2)*clhs160;
        const double clhs193 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs188;
        const double clhs194 =             C(2,4)*DN(2,2);
        const double clhs195 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs194;
        const double clhs196 =             clhs102*clhs152;
        const double clhs197 =             DN(0,1)*N[2];
        const double clhs198 =             DN(2,1)*rho*tau1;
        const double clhs199 =             DN(3,0)*clhs160;
        const double clhs200 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs130;
        const double clhs201 =             C(0,4)*DN(3,0) + clhs133 + clhs141;
        const double clhs202 =             clhs135*clhs150;
        const double clhs203 =             DN(3,1)*clhs160 + clhs116;
        const double clhs204 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
        const double clhs205 =             C(1,4)*DN(3,1);
        const double clhs206 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs205;
        const double clhs207 =             N[3]*clhs151 + clhs123;
        const double clhs208 =             clhs125 + clhs207;
        const double clhs209 =             DN(3,2)*clhs160;
        const double clhs210 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs205;
        const double clhs211 =             C(2,4)*DN(3,2);
        const double clhs212 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs211;
        const double clhs213 =             clhs135*clhs152;
        const double clhs214 =             DN(0,1)*N[3];
        const double clhs215 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs39;
        const double clhs216 =             DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2);
        const double clhs217 =             DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2);
        const double clhs218 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
        const double clhs219 =             DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)) + clhs11*clhs216 + clhs12*clhs217 + clhs13*clhs218);
        const double clhs220 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs163;
        const double clhs221 =             pow(DN(0,2), 2);
        const double clhs222 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
        const double clhs223 =             N[0]*clhs218 + clhs14;
        const double clhs224 =             clhs18 + clhs223;
        const double clhs225 =             DN(0,2)*tau2;
        const double clhs226 =             DN(1,0)*clhs225;
        const double clhs227 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs76;
        const double clhs228 =             clhs216*clhs68;
        const double clhs229 =             DN(1,1)*clhs225;
        const double clhs230 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs177;
        const double clhs231 =             clhs217*clhs68;
        const double clhs232 =             DN(1,2)*clhs225 + clhs49;
        const double clhs233 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
        const double clhs234 =             N[1]*clhs218 + clhs56;
        const double clhs235 =             clhs234 + clhs58;
        const double clhs236 =             DN(0,2)*N[1];
        const double clhs237 =             DN(1,2)*rho*tau1;
        const double clhs238 =             DN(2,0)*clhs225;
        const double clhs239 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs110;
        const double clhs240 =             clhs102*clhs216;
        const double clhs241 =             DN(2,1)*clhs225;
        const double clhs242 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs194;
        const double clhs243 =             clhs102*clhs217;
        const double clhs244 =             DN(2,2)*clhs225 + clhs83;
        const double clhs245 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
        const double clhs246 =             N[2]*clhs218 + clhs90;
        const double clhs247 =             clhs246 + clhs92;
        const double clhs248 =             DN(0,2)*N[2];
        const double clhs249 =             DN(2,2)*rho*tau1;
        const double clhs250 =             DN(3,0)*clhs225;
        const double clhs251 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs143;
        const double clhs252 =             clhs135*clhs216;
        const double clhs253 =             DN(3,1)*clhs225;
        const double clhs254 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs211;
        const double clhs255 =             clhs135*clhs217;
        const double clhs256 =             DN(3,2)*clhs225 + clhs116;
        const double clhs257 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
        const double clhs258 =             N[3]*clhs218 + clhs123;
        const double clhs259 =             clhs125 + clhs258;
        const double clhs260 =             DN(0,2)*N[3];
        const double clhs261 =             DN(0,0)*rho*tau1;
        const double clhs262 =             DN(0,1)*rho*tau1;
        const double clhs263 =             DN(0,2)*rho*tau1;
        const double clhs264 =             bdf0*clhs42*clhs43;
        const double clhs265 =             DN(1,0)*N[0];
        const double clhs266 =             DN(1,1)*N[0];
        const double clhs267 =             DN(1,2)*N[0];
        const double clhs268 =             DN(0,0)*tau1;
        const double clhs269 =             DN(0,1)*tau1;
        const double clhs270 =             DN(0,2)*tau1;
        const double clhs271 =             DN(1,0)*clhs268 + DN(1,1)*clhs269 + DN(1,2)*clhs270 + N[1]*clhs44;
        const double clhs272 =             DN(2,0)*N[0];
        const double clhs273 =             DN(2,1)*N[0];
        const double clhs274 =             DN(2,2)*N[0];
        const double clhs275 =             DN(2,0)*clhs268 + DN(2,1)*clhs269 + DN(2,2)*clhs270 + N[2]*clhs44;
        const double clhs276 =             DN(3,0)*N[0];
        const double clhs277 =             DN(3,1)*N[0];
        const double clhs278 =             DN(3,2)*N[0];
        const double clhs279 =             DN(3,0)*clhs268 + DN(3,1)*clhs269 + DN(3,2)*clhs270 + N[3]*clhs44;
        const double clhs280 =             N[1]*rho;
        const double clhs281 =             clhs16*clhs56*tau1;
        const double clhs282 =             DN(1,0)*N[0]*rho*tau1;
        const double clhs283 =             N[0]*clhs16*clhs56*tau1;
        const double clhs284 =             DN(1,1)*N[0]*rho*tau1;
        const double clhs285 =             DN(1,2)*N[0]*rho*tau1;
        const double clhs286 =             pow(DN(1,0), 2);
        const double clhs287 =             pow(N[1], 2);
        const double clhs288 =             clhs2*clhs287;
        const double clhs289 =             DN(1,0)*N[1]*rho*tau1;
        const double clhs290 =             DN(1,0)*tau2;
        const double clhs291 =             DN(1,1)*clhs290;
        const double clhs292 =             clhs287*rho;
        const double clhs293 =             N[1]*clhs16*clhs56*tau1;
        const double clhs294 =             DN(1,1)*N[1]*rho*tau1;
        const double clhs295 =             DN(1,2)*clhs290;
        const double clhs296 =             DN(1,2)*N[1]*rho*tau1;
        const double clhs297 =             N[1]*bdf0*clhs42*clhs43;
        const double clhs298 =             clhs45*clhs56;
        const double clhs299 =             -N[1] + clhs297*tau2 + clhs298;
        const double clhs300 =             N[1]*bdf0*rho;
        const double clhs301 =             N[2]*clhs300;
        const double clhs302 =             DN(2,0)*clhs290 + clhs301;
        const double clhs303 =             DN(1,0)*N[2]*rho*tau1;
        const double clhs304 =             DN(2,1)*clhs290;
        const double clhs305 =             N[1]*N[2]*rho;
        const double clhs306 =             clhs21*clhs305;
        const double clhs307 =             N[2]*clhs16*clhs56*tau1;
        const double clhs308 =             DN(1,1)*N[2]*rho*tau1;
        const double clhs309 =             DN(2,2)*clhs290;
        const double clhs310 =             clhs22*clhs305;
        const double clhs311 =             DN(1,2)*N[2]*rho*tau1;
        const double clhs312 =             DN(1,0)*N[2];
        const double clhs313 =             N[3]*clhs300;
        const double clhs314 =             DN(3,0)*clhs290 + clhs313;
        const double clhs315 =             DN(1,0)*N[3]*rho*tau1;
        const double clhs316 =             DN(3,1)*clhs290;
        const double clhs317 =             N[1]*N[3]*rho;
        const double clhs318 =             clhs21*clhs317;
        const double clhs319 =             N[3]*clhs16*clhs56*tau1;
        const double clhs320 =             DN(1,1)*N[3]*rho*tau1;
        const double clhs321 =             DN(3,2)*clhs290;
        const double clhs322 =             clhs22*clhs317;
        const double clhs323 =             DN(1,2)*N[3]*rho*tau1;
        const double clhs324 =             DN(1,0)*N[3];
        const double clhs325 =             pow(DN(1,1), 2);
        const double clhs326 =             DN(1,1)*tau2;
        const double clhs327 =             DN(1,2)*clhs326;
        const double clhs328 =             DN(2,0)*clhs326;
        const double clhs329 =             clhs150*clhs305;
        const double clhs330 =             DN(2,1)*clhs326 + clhs301;
        const double clhs331 =             DN(2,2)*clhs326;
        const double clhs332 =             clhs152*clhs305;
        const double clhs333 =             DN(1,1)*N[2];
        const double clhs334 =             DN(3,0)*clhs326;
        const double clhs335 =             clhs150*clhs317;
        const double clhs336 =             DN(3,1)*clhs326 + clhs313;
        const double clhs337 =             DN(3,2)*clhs326;
        const double clhs338 =             clhs152*clhs317;
        const double clhs339 =             DN(1,1)*N[3];
        const double clhs340 =             pow(DN(1,2), 2);
        const double clhs341 =             DN(1,2)*tau2;
        const double clhs342 =             DN(2,0)*clhs341;
        const double clhs343 =             clhs216*clhs305;
        const double clhs344 =             DN(2,1)*clhs341;
        const double clhs345 =             clhs217*clhs305;
        const double clhs346 =             DN(2,2)*clhs341 + clhs301;
        const double clhs347 =             DN(1,2)*N[2];
        const double clhs348 =             DN(3,0)*clhs341;
        const double clhs349 =             clhs216*clhs317;
        const double clhs350 =             DN(3,1)*clhs341;
        const double clhs351 =             clhs217*clhs317;
        const double clhs352 =             DN(3,2)*clhs341 + clhs313;
        const double clhs353 =             DN(1,2)*N[3];
        const double clhs354 =             DN(2,0)*N[1];
        const double clhs355 =             DN(2,1)*N[1];
        const double clhs356 =             DN(2,2)*N[1];
        const double clhs357 =             DN(1,0)*tau1;
        const double clhs358 =             DN(1,1)*tau1;
        const double clhs359 =             DN(1,2)*tau1;
        const double clhs360 =             DN(2,0)*clhs357 + DN(2,1)*clhs358 + DN(2,2)*clhs359 + N[2]*clhs297;
        const double clhs361 =             DN(3,0)*N[1];
        const double clhs362 =             DN(3,1)*N[1];
        const double clhs363 =             DN(3,2)*N[1];
        const double clhs364 =             DN(3,0)*clhs357 + DN(3,1)*clhs358 + DN(3,2)*clhs359 + N[3]*clhs297;
        const double clhs365 =             N[2]*rho;
        const double clhs366 =             clhs16*clhs90*tau1;
        const double clhs367 =             DN(2,0)*N[0]*rho*tau1;
        const double clhs368 =             N[0]*clhs16*clhs90*tau1;
        const double clhs369 =             DN(2,1)*N[0]*rho*tau1;
        const double clhs370 =             DN(2,2)*N[0]*rho*tau1;
        const double clhs371 =             DN(2,0)*N[1]*rho*tau1;
        const double clhs372 =             N[1]*clhs16*clhs90*tau1;
        const double clhs373 =             DN(2,1)*N[1]*rho*tau1;
        const double clhs374 =             DN(2,2)*N[1]*rho*tau1;
        const double clhs375 =             pow(DN(2,0), 2);
        const double clhs376 =             pow(N[2], 2);
        const double clhs377 =             clhs2*clhs376;
        const double clhs378 =             DN(2,0)*N[2]*rho*tau1;
        const double clhs379 =             DN(2,0)*tau2;
        const double clhs380 =             DN(2,1)*clhs379;
        const double clhs381 =             clhs376*rho;
        const double clhs382 =             N[2]*clhs16*clhs90*tau1;
        const double clhs383 =             DN(2,1)*N[2]*rho*tau1;
        const double clhs384 =             DN(2,2)*clhs379;
        const double clhs385 =             DN(2,2)*N[2]*rho*tau1;
        const double clhs386 =             clhs42*clhs43*tau2;
        const double clhs387 =             clhs45*clhs90;
        const double clhs388 =             -N[2] + clhs386*clhs92 + clhs387;
        const double clhs389 =             N[2]*N[3]*bdf0;
        const double clhs390 =             clhs389*rho;
        const double clhs391 =             DN(3,0)*clhs379 + clhs390;
        const double clhs392 =             DN(2,0)*N[3]*rho*tau1;
        const double clhs393 =             DN(3,1)*clhs379;
        const double clhs394 =             N[2]*N[3]*rho;
        const double clhs395 =             clhs21*clhs394;
        const double clhs396 =             N[3]*clhs16*clhs90*tau1;
        const double clhs397 =             DN(2,1)*N[3]*rho*tau1;
        const double clhs398 =             DN(3,2)*clhs379;
        const double clhs399 =             clhs22*clhs394;
        const double clhs400 =             DN(2,2)*N[3]*rho*tau1;
        const double clhs401 =             DN(2,0)*N[3];
        const double clhs402 =             pow(DN(2,1), 2);
        const double clhs403 =             DN(2,1)*tau2;
        const double clhs404 =             DN(2,2)*clhs403;
        const double clhs405 =             DN(3,0)*clhs403;
        const double clhs406 =             clhs150*clhs394;
        const double clhs407 =             DN(3,1)*clhs403 + clhs390;
        const double clhs408 =             DN(3,2)*clhs403;
        const double clhs409 =             clhs152*clhs394;
        const double clhs410 =             DN(2,1)*N[3];
        const double clhs411 =             pow(DN(2,2), 2);
        const double clhs412 =             DN(2,2)*tau2;
        const double clhs413 =             DN(3,0)*clhs412;
        const double clhs414 =             clhs216*clhs394;
        const double clhs415 =             DN(3,1)*clhs412;
        const double clhs416 =             clhs217*clhs394;
        const double clhs417 =             DN(3,2)*clhs412 + clhs390;
        const double clhs418 =             DN(2,2)*N[3];
        const double clhs419 =             DN(3,0)*N[2];
        const double clhs420 =             DN(3,1)*N[2];
        const double clhs421 =             DN(3,2)*N[2];
        const double clhs422 =             DN(3,0)*(DN(2,0)*tau1) + DN(3,1)*(DN(2,1)*tau1) + DN(3,2)*(DN(2,2)*tau1) + clhs389*clhs42*clhs43;
        const double clhs423 =             N[3]*rho;
        const double clhs424 =             clhs123*clhs16*tau1;
        const double clhs425 =             DN(3,0)*N[0]*rho*tau1;
        const double clhs426 =             N[0]*clhs123*clhs16*tau1;
        const double clhs427 =             DN(3,1)*N[0]*rho*tau1;
        const double clhs428 =             DN(3,2)*N[0]*rho*tau1;
        const double clhs429 =             DN(3,0)*N[1]*rho*tau1;
        const double clhs430 =             N[1]*clhs123*clhs16*tau1;
        const double clhs431 =             DN(3,1)*N[1]*rho*tau1;
        const double clhs432 =             DN(3,2)*N[1]*rho*tau1;
        const double clhs433 =             DN(3,0)*N[2]*rho*tau1;
        const double clhs434 =             N[2]*clhs123*clhs16*tau1;
        const double clhs435 =             DN(3,1)*N[2]*rho*tau1;
        const double clhs436 =             DN(3,2)*N[2]*rho*tau1;
        const double clhs437 =             pow(DN(3,0), 2);
        const double clhs438 =             pow(N[3], 2);
        const double clhs439 =             clhs2*clhs438;
        const double clhs440 =             DN(3,0)*N[3]*rho*tau1;
        const double clhs441 =             DN(3,0)*tau2;
        const double clhs442 =             DN(3,1)*clhs441;
        const double clhs443 =             clhs438*rho;
        const double clhs444 =             N[3]*clhs123*clhs16*tau1;
        const double clhs445 =             DN(3,1)*N[3]*rho*tau1;
        const double clhs446 =             DN(3,2)*clhs441;
        const double clhs447 =             DN(3,2)*N[3]*rho*tau1;
        const double clhs448 =             -N[3] + clhs123*clhs45 + clhs125*clhs386;
        const double clhs449 =             pow(DN(3,1), 2);
        const double clhs450 =             DN(3,1)*DN(3,2)*tau2;
        const double clhs451 =             pow(DN(3,2), 2);
        const double clhs452 =             DN(3,0)*rho*tau1;
        const double clhs453 =             DN(3,1)*rho*tau1;
        const double clhs454 =             DN(3,2)*rho*tau1;

        lhs(0,0)=DN(0,0)*clhs4 + DN(0,1)*clhs6 + DN(0,2)*clhs8 + clhs0*tau2 + clhs15*clhs9 + clhs17*clhs19 + clhs20*clhs23 + clhs3;
        lhs(0,1)=DN(0,0)*clhs26 + DN(0,1)*clhs28 + DN(0,2)*clhs31 + clhs21*clhs32 + clhs21*clhs33 + clhs23*clhs34 + clhs25;
        lhs(0,2)=DN(0,0)*clhs36 + DN(0,1)*clhs38 + DN(0,2)*clhs40 + clhs22*clhs32 + clhs22*clhs33 + clhs23*clhs41 + clhs35;
        lhs(0,3)=DN(0,0)*clhs47;
        lhs(0,4)=DN(0,0)*clhs51 + DN(0,1)*clhs53 + DN(0,2)*clhs55 + clhs17*clhs59 + clhs23*clhs60 + clhs50 + clhs57*clhs9;
        lhs(0,5)=DN(0,0)*clhs62 + DN(0,1)*clhs64 + DN(0,2)*clhs67 + clhs21*clhs70 + clhs23*clhs71 + clhs61 + clhs69;
        lhs(0,6)=DN(0,0)*clhs73 + DN(0,1)*clhs75 + DN(0,2)*clhs77 + clhs22*clhs70 + clhs23*clhs79 + clhs72 + clhs78;
        lhs(0,7)=clhs14*clhs82 + clhs80*clhs81 - clhs80;
        lhs(0,8)=DN(0,0)*clhs85 + DN(0,1)*clhs87 + DN(0,2)*clhs89 + clhs17*clhs93 + clhs23*clhs94 + clhs84 + clhs9*clhs91;
        lhs(0,9)=DN(0,0)*clhs96 + DN(0,1)*clhs98 + DN(0,2)*clhs101 + clhs103 + clhs104*clhs21 + clhs105*clhs23 + clhs95;
        lhs(0,10)=DN(0,0)*clhs107 + DN(0,1)*clhs109 + DN(0,2)*clhs111 + clhs104*clhs22 + clhs106 + clhs112 + clhs113*clhs23;
        lhs(0,11)=clhs114*clhs81 - clhs114 + clhs115*clhs14;
        lhs(0,12)=DN(0,0)*clhs118 + DN(0,1)*clhs120 + DN(0,2)*clhs122 + clhs117 + clhs124*clhs9 + clhs126*clhs17 + clhs127*clhs23;
        lhs(0,13)=DN(0,0)*clhs129 + DN(0,1)*clhs131 + DN(0,2)*clhs134 + clhs128 + clhs136 + clhs137*clhs21 + clhs138*clhs23;
        lhs(0,14)=DN(0,0)*clhs140 + DN(0,1)*clhs142 + DN(0,2)*clhs144 + clhs137*clhs22 + clhs139 + clhs145 + clhs146*clhs23;
        lhs(0,15)=DN(3,0)*clhs46 + clhs147*clhs81 - clhs147;
        lhs(1,0)=DN(0,0)*clhs6 + DN(0,1)*clhs148 + DN(0,2)*clhs149 + clhs150*clhs32 + clhs150*clhs33 + clhs153*clhs20 + clhs25;
        lhs(1,1)=DN(0,0)*clhs28 + DN(0,1)*clhs155 + DN(0,2)*clhs157 + clhs153*clhs34 + clhs154*tau2 + clhs158*clhs9 + clhs159*clhs17 + clhs3;
        lhs(1,2)=DN(0,0)*clhs38 + DN(0,1)*clhs162 + DN(0,2)*clhs164 + clhs152*clhs32 + clhs152*clhs33 + clhs153*clhs41 + clhs161;
        lhs(1,3)=DN(0,1)*clhs47;
        lhs(1,4)=DN(0,0)*clhs53 + DN(0,1)*clhs166 + DN(0,2)*clhs167 + clhs150*clhs70 + clhs153*clhs60 + clhs165 + clhs168;
        lhs(1,5)=DN(0,0)*clhs64 + DN(0,1)*clhs170 + DN(0,2)*clhs172 + clhs153*clhs71 + clhs169 + clhs17*clhs174 + clhs173*clhs9;
        lhs(1,6)=DN(0,0)*clhs75 + DN(0,1)*clhs176 + DN(0,2)*clhs178 + clhs152*clhs70 + clhs153*clhs79 + clhs175 + clhs179;
        lhs(1,7)=clhs14*clhs181 + clhs180*clhs81 - clhs180;
        lhs(1,8)=DN(0,0)*clhs87 + DN(0,1)*clhs183 + DN(0,2)*clhs184 + clhs104*clhs150 + clhs153*clhs94 + clhs182 + clhs185;
        lhs(1,9)=DN(0,0)*clhs98 + DN(0,1)*clhs187 + DN(0,2)*clhs189 + clhs105*clhs153 + clhs17*clhs191 + clhs186 + clhs190*clhs9;
        lhs(1,10)=DN(0,0)*clhs109 + DN(0,1)*clhs193 + DN(0,2)*clhs195 + clhs104*clhs152 + clhs113*clhs153 + clhs192 + clhs196;
        lhs(1,11)=clhs14*clhs198 + clhs197*clhs81 - clhs197;
        lhs(1,12)=DN(0,0)*clhs120 + DN(0,1)*clhs200 + DN(0,2)*clhs201 + clhs127*clhs153 + clhs137*clhs150 + clhs199 + clhs202;
        lhs(1,13)=DN(0,0)*clhs131 + DN(0,1)*clhs204 + DN(0,2)*clhs206 + clhs138*clhs153 + clhs17*clhs208 + clhs203 + clhs207*clhs9;
        lhs(1,14)=DN(0,0)*clhs142 + DN(0,1)*clhs210 + DN(0,2)*clhs212 + clhs137*clhs152 + clhs146*clhs153 + clhs209 + clhs213;
        lhs(1,15)=DN(3,1)*clhs46 + clhs214*clhs81 - clhs214;
        lhs(2,0)=DN(0,0)*clhs8 + DN(0,1)*clhs149 + DN(0,2)*clhs215 + clhs20*clhs219 + clhs216*clhs32 + clhs216*clhs33 + clhs35;
        lhs(2,1)=DN(0,0)*clhs31 + DN(0,1)*clhs157 + DN(0,2)*clhs220 + clhs161 + clhs217*clhs32 + clhs217*clhs33 + clhs219*clhs34;
        lhs(2,2)=DN(0,0)*clhs40 + DN(0,1)*clhs164 + DN(0,2)*clhs222 + clhs17*clhs224 + clhs219*clhs41 + clhs221*tau2 + clhs223*clhs9 + clhs3;
        lhs(2,3)=DN(0,2)*clhs47;
        lhs(2,4)=DN(0,0)*clhs55 + DN(0,1)*clhs167 + DN(0,2)*clhs227 + clhs216*clhs70 + clhs219*clhs60 + clhs226 + clhs228;
        lhs(2,5)=DN(0,0)*clhs67 + DN(0,1)*clhs172 + DN(0,2)*clhs230 + clhs217*clhs70 + clhs219*clhs71 + clhs229 + clhs231;
        lhs(2,6)=DN(0,0)*clhs77 + DN(0,1)*clhs178 + DN(0,2)*clhs233 + clhs17*clhs235 + clhs219*clhs79 + clhs232 + clhs234*clhs9;
        lhs(2,7)=clhs14*clhs237 + clhs236*clhs81 - clhs236;
        lhs(2,8)=DN(0,0)*clhs89 + DN(0,1)*clhs184 + DN(0,2)*clhs239 + clhs104*clhs216 + clhs219*clhs94 + clhs238 + clhs240;
        lhs(2,9)=DN(0,0)*clhs101 + DN(0,1)*clhs189 + DN(0,2)*clhs242 + clhs104*clhs217 + clhs105*clhs219 + clhs241 + clhs243;
        lhs(2,10)=DN(0,0)*clhs111 + DN(0,1)*clhs195 + DN(0,2)*clhs245 + clhs113*clhs219 + clhs17*clhs247 + clhs244 + clhs246*clhs9;
        lhs(2,11)=clhs14*clhs249 + clhs248*clhs81 - clhs248;
        lhs(2,12)=DN(0,0)*clhs122 + DN(0,1)*clhs201 + DN(0,2)*clhs251 + clhs127*clhs219 + clhs137*clhs216 + clhs250 + clhs252;
        lhs(2,13)=DN(0,0)*clhs134 + DN(0,1)*clhs206 + DN(0,2)*clhs254 + clhs137*clhs217 + clhs138*clhs219 + clhs253 + clhs255;
        lhs(2,14)=DN(0,0)*clhs144 + DN(0,1)*clhs212 + DN(0,2)*clhs257 + clhs146*clhs219 + clhs17*clhs259 + clhs256 + clhs258*clhs9;
        lhs(2,15)=DN(3,2)*clhs46 + clhs260*clhs81 - clhs260;
        lhs(3,0)=DN(0,0)*N[0] + clhs150*clhs34 + clhs19*clhs261 + clhs216*clhs41;
        lhs(3,1)=DN(0,1)*N[0] + clhs159*clhs262 + clhs20*clhs21 + clhs217*clhs41;
        lhs(3,2)=DN(0,2)*N[0] + clhs152*clhs34 + clhs20*clhs22 + clhs224*clhs263;
        lhs(3,3)=clhs0*tau1 + clhs1*clhs264 + clhs154*tau1 + clhs221*tau1;
        lhs(3,4)=clhs150*clhs71 + clhs216*clhs79 + clhs261*clhs59 + clhs265;
        lhs(3,5)=clhs174*clhs262 + clhs21*clhs60 + clhs217*clhs79 + clhs266;
        lhs(3,6)=clhs152*clhs71 + clhs22*clhs60 + clhs235*clhs263 + clhs267;
        lhs(3,7)=clhs271;
        lhs(3,8)=clhs105*clhs150 + clhs113*clhs216 + clhs261*clhs93 + clhs272;
        lhs(3,9)=clhs113*clhs217 + clhs191*clhs262 + clhs21*clhs94 + clhs273;
        lhs(3,10)=clhs105*clhs152 + clhs22*clhs94 + clhs247*clhs263 + clhs274;
        lhs(3,11)=clhs275;
        lhs(3,12)=clhs126*clhs261 + clhs138*clhs150 + clhs146*clhs216 + clhs276;
        lhs(3,13)=clhs127*clhs21 + clhs146*clhs217 + clhs208*clhs262 + clhs277;
        lhs(3,14)=clhs127*clhs22 + clhs138*clhs152 + clhs259*clhs263 + clhs278;
        lhs(3,15)=clhs279;
        lhs(4,0)=DN(1,0)*clhs4 + DN(1,1)*clhs6 + DN(1,2)*clhs8 + clhs15*clhs280 + clhs19*clhs281 + clhs23*clhs282 + clhs50;
        lhs(4,1)=DN(1,0)*clhs26 + DN(1,1)*clhs28 + DN(1,2)*clhs31 + clhs165 + clhs21*clhs283 + clhs23*clhs284 + clhs69;
        lhs(4,2)=DN(1,0)*clhs36 + DN(1,1)*clhs38 + DN(1,2)*clhs40 + clhs22*clhs283 + clhs226 + clhs23*clhs285 + clhs78;
        lhs(4,3)=clhs261*clhs56 + clhs265*clhs81 - clhs265;
        lhs(4,4)=DN(1,0)*clhs51 + DN(1,1)*clhs53 + DN(1,2)*clhs55 + clhs23*clhs289 + clhs280*clhs57 + clhs281*clhs59 + clhs286*tau2 + clhs288;
        lhs(4,5)=DN(1,0)*clhs62 + DN(1,1)*clhs64 + DN(1,2)*clhs67 + clhs21*clhs292 + clhs21*clhs293 + clhs23*clhs294 + clhs291;
        lhs(4,6)=DN(1,0)*clhs73 + DN(1,1)*clhs75 + DN(1,2)*clhs77 + clhs22*clhs292 + clhs22*clhs293 + clhs23*clhs296 + clhs295;
        lhs(4,7)=DN(1,0)*clhs299;
        lhs(4,8)=DN(1,0)*clhs85 + DN(1,1)*clhs87 + DN(1,2)*clhs89 + clhs23*clhs303 + clhs280*clhs91 + clhs281*clhs93 + clhs302;
        lhs(4,9)=DN(1,0)*clhs96 + DN(1,1)*clhs98 + DN(1,2)*clhs101 + clhs21*clhs307 + clhs23*clhs308 + clhs304 + clhs306;
        lhs(4,10)=DN(1,0)*clhs107 + DN(1,1)*clhs109 + DN(1,2)*clhs111 + clhs22*clhs307 + clhs23*clhs311 + clhs309 + clhs310;
        lhs(4,11)=clhs115*clhs56 + clhs312*clhs81 - clhs312;
        lhs(4,12)=DN(1,0)*clhs118 + DN(1,1)*clhs120 + DN(1,2)*clhs122 + clhs124*clhs280 + clhs126*clhs281 + clhs23*clhs315 + clhs314;
        lhs(4,13)=DN(1,0)*clhs129 + DN(1,1)*clhs131 + DN(1,2)*clhs134 + clhs21*clhs319 + clhs23*clhs320 + clhs316 + clhs318;
        lhs(4,14)=DN(1,0)*clhs140 + DN(1,1)*clhs142 + DN(1,2)*clhs144 + clhs22*clhs319 + clhs23*clhs323 + clhs321 + clhs322;
        lhs(4,15)=DN(3,0)*clhs298 + clhs324*clhs81 - clhs324;
        lhs(5,0)=DN(1,0)*clhs6 + DN(1,1)*clhs148 + DN(1,2)*clhs149 + clhs150*clhs283 + clhs153*clhs282 + clhs168 + clhs61;
        lhs(5,1)=DN(1,0)*clhs28 + DN(1,1)*clhs155 + DN(1,2)*clhs157 + clhs153*clhs284 + clhs158*clhs280 + clhs159*clhs281 + clhs169;
        lhs(5,2)=DN(1,0)*clhs38 + DN(1,1)*clhs162 + DN(1,2)*clhs164 + clhs152*clhs283 + clhs153*clhs285 + clhs179 + clhs229;
        lhs(5,3)=clhs262*clhs56 + clhs266*clhs81 - clhs266;
        lhs(5,4)=DN(1,0)*clhs53 + DN(1,1)*clhs166 + DN(1,2)*clhs167 + clhs150*clhs292 + clhs150*clhs293 + clhs153*clhs289 + clhs291;
        lhs(5,5)=DN(1,0)*clhs64 + DN(1,1)*clhs170 + DN(1,2)*clhs172 + clhs153*clhs294 + clhs173*clhs280 + clhs174*clhs281 + clhs288 + clhs325*tau2;
        lhs(5,6)=DN(1,0)*clhs75 + DN(1,1)*clhs176 + DN(1,2)*clhs178 + clhs152*clhs292 + clhs152*clhs293 + clhs153*clhs296 + clhs327;
        lhs(5,7)=DN(1,1)*clhs299;
        lhs(5,8)=DN(1,0)*clhs87 + DN(1,1)*clhs183 + DN(1,2)*clhs184 + clhs150*clhs307 + clhs153*clhs303 + clhs328 + clhs329;
        lhs(5,9)=DN(1,0)*clhs98 + DN(1,1)*clhs187 + DN(1,2)*clhs189 + clhs153*clhs308 + clhs190*clhs280 + clhs191*clhs281 + clhs330;
        lhs(5,10)=DN(1,0)*clhs109 + DN(1,1)*clhs193 + DN(1,2)*clhs195 + clhs152*clhs307 + clhs153*clhs311 + clhs331 + clhs332;
        lhs(5,11)=clhs198*clhs56 + clhs333*clhs81 - clhs333;
        lhs(5,12)=DN(1,0)*clhs120 + DN(1,1)*clhs200 + DN(1,2)*clhs201 + clhs150*clhs319 + clhs153*clhs315 + clhs334 + clhs335;
        lhs(5,13)=DN(1,0)*clhs131 + DN(1,1)*clhs204 + DN(1,2)*clhs206 + clhs153*clhs320 + clhs207*clhs280 + clhs208*clhs281 + clhs336;
        lhs(5,14)=DN(1,0)*clhs142 + DN(1,1)*clhs210 + DN(1,2)*clhs212 + clhs152*clhs319 + clhs153*clhs323 + clhs337 + clhs338;
        lhs(5,15)=DN(3,1)*clhs298 + clhs339*clhs81 - clhs339;
        lhs(6,0)=DN(1,0)*clhs8 + DN(1,1)*clhs149 + DN(1,2)*clhs215 + clhs216*clhs283 + clhs219*clhs282 + clhs228 + clhs72;
        lhs(6,1)=DN(1,0)*clhs31 + DN(1,1)*clhs157 + DN(1,2)*clhs220 + clhs175 + clhs217*clhs283 + clhs219*clhs284 + clhs231;
        lhs(6,2)=DN(1,0)*clhs40 + DN(1,1)*clhs164 + DN(1,2)*clhs222 + clhs219*clhs285 + clhs223*clhs280 + clhs224*clhs281 + clhs232;
        lhs(6,3)=clhs263*clhs56 + clhs267*clhs81 - clhs267;
        lhs(6,4)=DN(1,0)*clhs55 + DN(1,1)*clhs167 + DN(1,2)*clhs227 + clhs216*clhs292 + clhs216*clhs293 + clhs219*clhs289 + clhs295;
        lhs(6,5)=DN(1,0)*clhs67 + DN(1,1)*clhs172 + DN(1,2)*clhs230 + clhs217*clhs292 + clhs217*clhs293 + clhs219*clhs294 + clhs327;
        lhs(6,6)=DN(1,0)*clhs77 + DN(1,1)*clhs178 + DN(1,2)*clhs233 + clhs219*clhs296 + clhs234*clhs280 + clhs235*clhs281 + clhs288 + clhs340*tau2;
        lhs(6,7)=DN(1,2)*clhs299;
        lhs(6,8)=DN(1,0)*clhs89 + DN(1,1)*clhs184 + DN(1,2)*clhs239 + clhs216*clhs307 + clhs219*clhs303 + clhs342 + clhs343;
        lhs(6,9)=DN(1,0)*clhs101 + DN(1,1)*clhs189 + DN(1,2)*clhs242 + clhs217*clhs307 + clhs219*clhs308 + clhs344 + clhs345;
        lhs(6,10)=DN(1,0)*clhs111 + DN(1,1)*clhs195 + DN(1,2)*clhs245 + clhs219*clhs311 + clhs246*clhs280 + clhs247*clhs281 + clhs346;
        lhs(6,11)=clhs249*clhs56 + clhs347*clhs81 - clhs347;
        lhs(6,12)=DN(1,0)*clhs122 + DN(1,1)*clhs201 + DN(1,2)*clhs251 + clhs216*clhs319 + clhs219*clhs315 + clhs348 + clhs349;
        lhs(6,13)=DN(1,0)*clhs134 + DN(1,1)*clhs206 + DN(1,2)*clhs254 + clhs217*clhs319 + clhs219*clhs320 + clhs350 + clhs351;
        lhs(6,14)=DN(1,0)*clhs144 + DN(1,1)*clhs212 + DN(1,2)*clhs257 + clhs219*clhs323 + clhs258*clhs280 + clhs259*clhs281 + clhs352;
        lhs(6,15)=DN(3,2)*clhs298 + clhs353*clhs81 - clhs353;
        lhs(7,0)=clhs150*clhs284 + clhs19*clhs82 + clhs216*clhs285 + clhs80;
        lhs(7,1)=clhs159*clhs181 + clhs180 + clhs21*clhs282 + clhs217*clhs285;
        lhs(7,2)=clhs152*clhs284 + clhs22*clhs282 + clhs224*clhs237 + clhs236;
        lhs(7,3)=clhs271;
        lhs(7,4)=DN(1,0)*N[1] + clhs150*clhs294 + clhs216*clhs296 + clhs59*clhs82;
        lhs(7,5)=DN(1,1)*N[1] + clhs174*clhs181 + clhs21*clhs289 + clhs217*clhs296;
        lhs(7,6)=DN(1,2)*N[1] + clhs152*clhs294 + clhs22*clhs289 + clhs235*clhs237;
        lhs(7,7)=clhs264*clhs287 + clhs286*tau1 + clhs325*tau1 + clhs340*tau1;
        lhs(7,8)=clhs150*clhs308 + clhs216*clhs311 + clhs354 + clhs82*clhs93;
        lhs(7,9)=clhs181*clhs191 + clhs21*clhs303 + clhs217*clhs311 + clhs355;
        lhs(7,10)=clhs152*clhs308 + clhs22*clhs303 + clhs237*clhs247 + clhs356;
        lhs(7,11)=clhs360;
        lhs(7,12)=clhs126*clhs82 + clhs150*clhs320 + clhs216*clhs323 + clhs361;
        lhs(7,13)=clhs181*clhs208 + clhs21*clhs315 + clhs217*clhs323 + clhs362;
        lhs(7,14)=clhs152*clhs320 + clhs22*clhs315 + clhs237*clhs259 + clhs363;
        lhs(7,15)=clhs364;
        lhs(8,0)=DN(2,0)*clhs4 + DN(2,1)*clhs6 + DN(2,2)*clhs8 + clhs15*clhs365 + clhs19*clhs366 + clhs23*clhs367 + clhs84;
        lhs(8,1)=DN(2,0)*clhs26 + DN(2,1)*clhs28 + DN(2,2)*clhs31 + clhs103 + clhs182 + clhs21*clhs368 + clhs23*clhs369;
        lhs(8,2)=DN(2,0)*clhs36 + DN(2,1)*clhs38 + DN(2,2)*clhs40 + clhs112 + clhs22*clhs368 + clhs23*clhs370 + clhs238;
        lhs(8,3)=clhs261*clhs90 + clhs272*clhs81 - clhs272;
        lhs(8,4)=DN(2,0)*clhs51 + DN(2,1)*clhs53 + DN(2,2)*clhs55 + clhs23*clhs371 + clhs302 + clhs365*clhs57 + clhs366*clhs59;
        lhs(8,5)=DN(2,0)*clhs62 + DN(2,1)*clhs64 + DN(2,2)*clhs67 + clhs21*clhs372 + clhs23*clhs373 + clhs306 + clhs328;
        lhs(8,6)=DN(2,0)*clhs73 + DN(2,1)*clhs75 + DN(2,2)*clhs77 + clhs22*clhs372 + clhs23*clhs374 + clhs310 + clhs342;
        lhs(8,7)=clhs354*clhs81 - clhs354 + clhs82*clhs90;
        lhs(8,8)=DN(2,0)*clhs85 + DN(2,1)*clhs87 + DN(2,2)*clhs89 + clhs23*clhs378 + clhs365*clhs91 + clhs366*clhs93 + clhs375*tau2 + clhs377;
        lhs(8,9)=DN(2,0)*clhs96 + DN(2,1)*clhs98 + DN(2,2)*clhs101 + clhs21*clhs381 + clhs21*clhs382 + clhs23*clhs383 + clhs380;
        lhs(8,10)=DN(2,0)*clhs107 + DN(2,1)*clhs109 + DN(2,2)*clhs111 + clhs22*clhs381 + clhs22*clhs382 + clhs23*clhs385 + clhs384;
        lhs(8,11)=DN(2,0)*clhs388;
        lhs(8,12)=DN(2,0)*clhs118 + DN(2,1)*clhs120 + DN(2,2)*clhs122 + clhs124*clhs365 + clhs126*clhs366 + clhs23*clhs392 + clhs391;
        lhs(8,13)=DN(2,0)*clhs129 + DN(2,1)*clhs131 + DN(2,2)*clhs134 + clhs21*clhs396 + clhs23*clhs397 + clhs393 + clhs395;
        lhs(8,14)=DN(2,0)*clhs140 + DN(2,1)*clhs142 + DN(2,2)*clhs144 + clhs22*clhs396 + clhs23*clhs400 + clhs398 + clhs399;
        lhs(8,15)=DN(3,0)*clhs387 + clhs401*clhs81 - clhs401;
        lhs(9,0)=DN(2,0)*clhs6 + DN(2,1)*clhs148 + DN(2,2)*clhs149 + clhs150*clhs368 + clhs153*clhs367 + clhs185 + clhs95;
        lhs(9,1)=DN(2,0)*clhs28 + DN(2,1)*clhs155 + DN(2,2)*clhs157 + clhs153*clhs369 + clhs158*clhs365 + clhs159*clhs366 + clhs186;
        lhs(9,2)=DN(2,0)*clhs38 + DN(2,1)*clhs162 + DN(2,2)*clhs164 + clhs152*clhs368 + clhs153*clhs370 + clhs196 + clhs241;
        lhs(9,3)=clhs262*clhs90 + clhs273*clhs81 - clhs273;
        lhs(9,4)=DN(2,0)*clhs53 + DN(2,1)*clhs166 + DN(2,2)*clhs167 + clhs150*clhs372 + clhs153*clhs371 + clhs304 + clhs329;
        lhs(9,5)=DN(2,0)*clhs64 + DN(2,1)*clhs170 + DN(2,2)*clhs172 + clhs153*clhs373 + clhs173*clhs365 + clhs174*clhs366 + clhs330;
        lhs(9,6)=DN(2,0)*clhs75 + DN(2,1)*clhs176 + DN(2,2)*clhs178 + clhs152*clhs372 + clhs153*clhs374 + clhs332 + clhs344;
        lhs(9,7)=clhs181*clhs90 + clhs355*clhs81 - clhs355;
        lhs(9,8)=DN(2,0)*clhs87 + DN(2,1)*clhs183 + DN(2,2)*clhs184 + clhs150*clhs381 + clhs150*clhs382 + clhs153*clhs378 + clhs380;
        lhs(9,9)=DN(2,0)*clhs98 + DN(2,1)*clhs187 + DN(2,2)*clhs189 + clhs153*clhs383 + clhs190*clhs365 + clhs191*clhs366 + clhs377 + clhs402*tau2;
        lhs(9,10)=DN(2,0)*clhs109 + DN(2,1)*clhs193 + DN(2,2)*clhs195 + clhs152*clhs381 + clhs152*clhs382 + clhs153*clhs385 + clhs404;
        lhs(9,11)=DN(2,1)*clhs388;
        lhs(9,12)=DN(2,0)*clhs120 + DN(2,1)*clhs200 + DN(2,2)*clhs201 + clhs150*clhs396 + clhs153*clhs392 + clhs405 + clhs406;
        lhs(9,13)=DN(2,0)*clhs131 + DN(2,1)*clhs204 + DN(2,2)*clhs206 + clhs153*clhs397 + clhs207*clhs365 + clhs208*clhs366 + clhs407;
        lhs(9,14)=DN(2,0)*clhs142 + DN(2,1)*clhs210 + DN(2,2)*clhs212 + clhs152*clhs396 + clhs153*clhs400 + clhs408 + clhs409;
        lhs(9,15)=DN(3,1)*clhs387 + clhs410*clhs81 - clhs410;
        lhs(10,0)=DN(2,0)*clhs8 + DN(2,1)*clhs149 + DN(2,2)*clhs215 + clhs106 + clhs216*clhs368 + clhs219*clhs367 + clhs240;
        lhs(10,1)=DN(2,0)*clhs31 + DN(2,1)*clhs157 + DN(2,2)*clhs220 + clhs192 + clhs217*clhs368 + clhs219*clhs369 + clhs243;
        lhs(10,2)=DN(2,0)*clhs40 + DN(2,1)*clhs164 + DN(2,2)*clhs222 + clhs219*clhs370 + clhs223*clhs365 + clhs224*clhs366 + clhs244;
        lhs(10,3)=clhs263*clhs90 + clhs274*clhs81 - clhs274;
        lhs(10,4)=DN(2,0)*clhs55 + DN(2,1)*clhs167 + DN(2,2)*clhs227 + clhs216*clhs372 + clhs219*clhs371 + clhs309 + clhs343;
        lhs(10,5)=DN(2,0)*clhs67 + DN(2,1)*clhs172 + DN(2,2)*clhs230 + clhs217*clhs372 + clhs219*clhs373 + clhs331 + clhs345;
        lhs(10,6)=DN(2,0)*clhs77 + DN(2,1)*clhs178 + DN(2,2)*clhs233 + clhs219*clhs374 + clhs234*clhs365 + clhs235*clhs366 + clhs346;
        lhs(10,7)=clhs237*clhs90 + clhs356*clhs81 - clhs356;
        lhs(10,8)=DN(2,0)*clhs89 + DN(2,1)*clhs184 + DN(2,2)*clhs239 + clhs216*clhs381 + clhs216*clhs382 + clhs219*clhs378 + clhs384;
        lhs(10,9)=DN(2,0)*clhs101 + DN(2,1)*clhs189 + DN(2,2)*clhs242 + clhs217*clhs381 + clhs217*clhs382 + clhs219*clhs383 + clhs404;
        lhs(10,10)=DN(2,0)*clhs111 + DN(2,1)*clhs195 + DN(2,2)*clhs245 + clhs219*clhs385 + clhs246*clhs365 + clhs247*clhs366 + clhs377 + clhs411*tau2;
        lhs(10,11)=DN(2,2)*clhs388;
        lhs(10,12)=DN(2,0)*clhs122 + DN(2,1)*clhs201 + DN(2,2)*clhs251 + clhs216*clhs396 + clhs219*clhs392 + clhs413 + clhs414;
        lhs(10,13)=DN(2,0)*clhs134 + DN(2,1)*clhs206 + DN(2,2)*clhs254 + clhs217*clhs396 + clhs219*clhs397 + clhs415 + clhs416;
        lhs(10,14)=DN(2,0)*clhs144 + DN(2,1)*clhs212 + DN(2,2)*clhs257 + clhs219*clhs400 + clhs258*clhs365 + clhs259*clhs366 + clhs417;
        lhs(10,15)=DN(3,2)*clhs387 + clhs418*clhs81 - clhs418;
        lhs(11,0)=clhs114 + clhs115*clhs19 + clhs150*clhs369 + clhs216*clhs370;
        lhs(11,1)=clhs159*clhs198 + clhs197 + clhs21*clhs367 + clhs217*clhs370;
        lhs(11,2)=clhs152*clhs369 + clhs22*clhs367 + clhs224*clhs249 + clhs248;
        lhs(11,3)=clhs275;
        lhs(11,4)=clhs115*clhs59 + clhs150*clhs373 + clhs216*clhs374 + clhs312;
        lhs(11,5)=clhs174*clhs198 + clhs21*clhs371 + clhs217*clhs374 + clhs333;
        lhs(11,6)=clhs152*clhs373 + clhs22*clhs371 + clhs235*clhs249 + clhs347;
        lhs(11,7)=clhs360;
        lhs(11,8)=DN(2,0)*N[2] + clhs115*clhs93 + clhs150*clhs383 + clhs216*clhs385;
        lhs(11,9)=DN(2,1)*N[2] + clhs191*clhs198 + clhs21*clhs378 + clhs217*clhs385;
        lhs(11,10)=DN(2,2)*N[2] + clhs152*clhs383 + clhs22*clhs378 + clhs247*clhs249;
        lhs(11,11)=clhs264*clhs376 + clhs375*tau1 + clhs402*tau1 + clhs411*tau1;
        lhs(11,12)=clhs115*clhs126 + clhs150*clhs397 + clhs216*clhs400 + clhs419;
        lhs(11,13)=clhs198*clhs208 + clhs21*clhs392 + clhs217*clhs400 + clhs420;
        lhs(11,14)=clhs152*clhs397 + clhs22*clhs392 + clhs249*clhs259 + clhs421;
        lhs(11,15)=clhs422;
        lhs(12,0)=DN(3,0)*clhs4 + DN(3,1)*clhs6 + DN(3,2)*clhs8 + clhs117 + clhs15*clhs423 + clhs19*clhs424 + clhs23*clhs425;
        lhs(12,1)=DN(3,0)*clhs26 + DN(3,1)*clhs28 + DN(3,2)*clhs31 + clhs136 + clhs199 + clhs21*clhs426 + clhs23*clhs427;
        lhs(12,2)=DN(3,0)*clhs36 + DN(3,1)*clhs38 + DN(3,2)*clhs40 + clhs145 + clhs22*clhs426 + clhs23*clhs428 + clhs250;
        lhs(12,3)=clhs123*clhs261 + clhs276*clhs81 - clhs276;
        lhs(12,4)=DN(3,0)*clhs51 + DN(3,1)*clhs53 + DN(3,2)*clhs55 + clhs23*clhs429 + clhs314 + clhs423*clhs57 + clhs424*clhs59;
        lhs(12,5)=DN(3,0)*clhs62 + DN(3,1)*clhs64 + DN(3,2)*clhs67 + clhs21*clhs430 + clhs23*clhs431 + clhs318 + clhs334;
        lhs(12,6)=DN(3,0)*clhs73 + DN(3,1)*clhs75 + DN(3,2)*clhs77 + clhs22*clhs430 + clhs23*clhs432 + clhs322 + clhs348;
        lhs(12,7)=clhs123*clhs82 + clhs361*clhs81 - clhs361;
        lhs(12,8)=DN(3,0)*clhs85 + DN(3,1)*clhs87 + DN(3,2)*clhs89 + clhs23*clhs433 + clhs391 + clhs423*clhs91 + clhs424*clhs93;
        lhs(12,9)=DN(3,0)*clhs96 + DN(3,1)*clhs98 + DN(3,2)*clhs101 + clhs21*clhs434 + clhs23*clhs435 + clhs395 + clhs405;
        lhs(12,10)=DN(3,0)*clhs107 + DN(3,1)*clhs109 + DN(3,2)*clhs111 + clhs22*clhs434 + clhs23*clhs436 + clhs399 + clhs413;
        lhs(12,11)=clhs115*clhs123 + clhs419*clhs81 - clhs419;
        lhs(12,12)=DN(3,0)*clhs118 + DN(3,1)*clhs120 + DN(3,2)*clhs122 + clhs124*clhs423 + clhs126*clhs424 + clhs23*clhs440 + clhs437*tau2 + clhs439;
        lhs(12,13)=DN(3,0)*clhs129 + DN(3,1)*clhs131 + DN(3,2)*clhs134 + clhs21*clhs443 + clhs21*clhs444 + clhs23*clhs445 + clhs442;
        lhs(12,14)=DN(3,0)*clhs140 + DN(3,1)*clhs142 + DN(3,2)*clhs144 + clhs22*clhs443 + clhs22*clhs444 + clhs23*clhs447 + clhs446;
        lhs(12,15)=DN(3,0)*clhs448;
        lhs(13,0)=DN(3,0)*clhs6 + DN(3,1)*clhs148 + DN(3,2)*clhs149 + clhs128 + clhs150*clhs426 + clhs153*clhs425 + clhs202;
        lhs(13,1)=DN(3,0)*clhs28 + DN(3,1)*clhs155 + DN(3,2)*clhs157 + clhs153*clhs427 + clhs158*clhs423 + clhs159*clhs424 + clhs203;
        lhs(13,2)=DN(3,0)*clhs38 + DN(3,1)*clhs162 + DN(3,2)*clhs164 + clhs152*clhs426 + clhs153*clhs428 + clhs213 + clhs253;
        lhs(13,3)=clhs123*clhs262 + clhs277*clhs81 - clhs277;
        lhs(13,4)=DN(3,0)*clhs53 + DN(3,1)*clhs166 + DN(3,2)*clhs167 + clhs150*clhs430 + clhs153*clhs429 + clhs316 + clhs335;
        lhs(13,5)=DN(3,0)*clhs64 + DN(3,1)*clhs170 + DN(3,2)*clhs172 + clhs153*clhs431 + clhs173*clhs423 + clhs174*clhs424 + clhs336;
        lhs(13,6)=DN(3,0)*clhs75 + DN(3,1)*clhs176 + DN(3,2)*clhs178 + clhs152*clhs430 + clhs153*clhs432 + clhs338 + clhs350;
        lhs(13,7)=clhs123*clhs181 + clhs362*clhs81 - clhs362;
        lhs(13,8)=DN(3,0)*clhs87 + DN(3,1)*clhs183 + DN(3,2)*clhs184 + clhs150*clhs434 + clhs153*clhs433 + clhs393 + clhs406;
        lhs(13,9)=DN(3,0)*clhs98 + DN(3,1)*clhs187 + DN(3,2)*clhs189 + clhs153*clhs435 + clhs190*clhs423 + clhs191*clhs424 + clhs407;
        lhs(13,10)=DN(3,0)*clhs109 + DN(3,1)*clhs193 + DN(3,2)*clhs195 + clhs152*clhs434 + clhs153*clhs436 + clhs409 + clhs415;
        lhs(13,11)=clhs123*clhs198 + clhs420*clhs81 - clhs420;
        lhs(13,12)=DN(3,0)*clhs120 + DN(3,1)*clhs200 + DN(3,2)*clhs201 + clhs150*clhs443 + clhs150*clhs444 + clhs153*clhs440 + clhs442;
        lhs(13,13)=DN(3,0)*clhs131 + DN(3,1)*clhs204 + DN(3,2)*clhs206 + clhs153*clhs445 + clhs207*clhs423 + clhs208*clhs424 + clhs439 + clhs449*tau2;
        lhs(13,14)=DN(3,0)*clhs142 + DN(3,1)*clhs210 + DN(3,2)*clhs212 + clhs152*clhs443 + clhs152*clhs444 + clhs153*clhs447 + clhs450;
        lhs(13,15)=DN(3,1)*clhs448;
        lhs(14,0)=DN(3,0)*clhs8 + DN(3,1)*clhs149 + DN(3,2)*clhs215 + clhs139 + clhs216*clhs426 + clhs219*clhs425 + clhs252;
        lhs(14,1)=DN(3,0)*clhs31 + DN(3,1)*clhs157 + DN(3,2)*clhs220 + clhs209 + clhs217*clhs426 + clhs219*clhs427 + clhs255;
        lhs(14,2)=DN(3,0)*clhs40 + DN(3,1)*clhs164 + DN(3,2)*clhs222 + clhs219*clhs428 + clhs223*clhs423 + clhs224*clhs424 + clhs256;
        lhs(14,3)=clhs123*clhs263 + clhs278*clhs81 - clhs278;
        lhs(14,4)=DN(3,0)*clhs55 + DN(3,1)*clhs167 + DN(3,2)*clhs227 + clhs216*clhs430 + clhs219*clhs429 + clhs321 + clhs349;
        lhs(14,5)=DN(3,0)*clhs67 + DN(3,1)*clhs172 + DN(3,2)*clhs230 + clhs217*clhs430 + clhs219*clhs431 + clhs337 + clhs351;
        lhs(14,6)=DN(3,0)*clhs77 + DN(3,1)*clhs178 + DN(3,2)*clhs233 + clhs219*clhs432 + clhs234*clhs423 + clhs235*clhs424 + clhs352;
        lhs(14,7)=clhs123*clhs237 + clhs363*clhs81 - clhs363;
        lhs(14,8)=DN(3,0)*clhs89 + DN(3,1)*clhs184 + DN(3,2)*clhs239 + clhs216*clhs434 + clhs219*clhs433 + clhs398 + clhs414;
        lhs(14,9)=DN(3,0)*clhs101 + DN(3,1)*clhs189 + DN(3,2)*clhs242 + clhs217*clhs434 + clhs219*clhs435 + clhs408 + clhs416;
        lhs(14,10)=DN(3,0)*clhs111 + DN(3,1)*clhs195 + DN(3,2)*clhs245 + clhs219*clhs436 + clhs246*clhs423 + clhs247*clhs424 + clhs417;
        lhs(14,11)=clhs123*clhs249 + clhs421*clhs81 - clhs421;
        lhs(14,12)=DN(3,0)*clhs122 + DN(3,1)*clhs201 + DN(3,2)*clhs251 + clhs216*clhs443 + clhs216*clhs444 + clhs219*clhs440 + clhs446;
        lhs(14,13)=DN(3,0)*clhs134 + DN(3,1)*clhs206 + DN(3,2)*clhs254 + clhs217*clhs443 + clhs217*clhs444 + clhs219*clhs445 + clhs450;
        lhs(14,14)=DN(3,0)*clhs144 + DN(3,1)*clhs212 + DN(3,2)*clhs257 + clhs219*clhs447 + clhs258*clhs423 + clhs259*clhs424 + clhs439 + clhs451*tau2;
        lhs(14,15)=DN(3,2)*clhs448;
        lhs(15,0)=clhs147 + clhs150*clhs427 + clhs19*clhs452 + clhs216*clhs428;
        lhs(15,1)=clhs159*clhs453 + clhs21*clhs425 + clhs214 + clhs217*clhs428;
        lhs(15,2)=clhs152*clhs427 + clhs22*clhs425 + clhs224*clhs454 + clhs260;
        lhs(15,3)=clhs279;
        lhs(15,4)=clhs150*clhs431 + clhs216*clhs432 + clhs324 + clhs452*clhs59;
        lhs(15,5)=clhs174*clhs453 + clhs21*clhs429 + clhs217*clhs432 + clhs339;
        lhs(15,6)=clhs152*clhs431 + clhs22*clhs429 + clhs235*clhs454 + clhs353;
        lhs(15,7)=clhs364;
        lhs(15,8)=clhs150*clhs435 + clhs216*clhs436 + clhs401 + clhs452*clhs93;
        lhs(15,9)=clhs191*clhs453 + clhs21*clhs433 + clhs217*clhs436 + clhs410;
        lhs(15,10)=clhs152*clhs435 + clhs22*clhs433 + clhs247*clhs454 + clhs418;
        lhs(15,11)=clhs422;
        lhs(15,12)=DN(3,0)*N[3] + clhs126*clhs452 + clhs150*clhs445 + clhs216*clhs447;
        lhs(15,13)=DN(3,1)*N[3] + clhs208*clhs453 + clhs21*clhs440 + clhs217*clhs447;
        lhs(15,14)=DN(3,2)*N[3] + clhs152*clhs445 + clhs22*clhs440 + clhs259*clhs454;
        lhs(15,15)=clhs264*clhs438 + clhs437*tau1 + clhs449*tau1 + clhs451*tau1;

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
        const double clhs28 =             pow(c, -2);
        const double clhs29 =             1.0/rho;
        const double clhs30 =             N[0]*bdf0*clhs28*clhs29;
        const double clhs31 =             rho*tau1;
        const double clhs32 =             clhs11*clhs31;
        const double clhs33 =             -N[0] + clhs30*tau2 + clhs32;
        const double clhs34 =             N[0]*bdf0*rho;
        const double clhs35 =             N[1]*clhs34;
        const double clhs36 =             DN(1,0)*clhs20 + clhs35;
        const double clhs37 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
        const double clhs38 =             C(0,2)*DN(1,0);
        const double clhs39 =             C(2,2)*DN(1,1) + clhs38;
        const double clhs40 =             DN(1,0)*clhs9 + DN(1,1)*clhs10;
        const double clhs41 =             N[1]*clhs8 + clhs40;
        const double clhs42 =             N[1]*bdf0;
        const double clhs43 =             clhs41 + clhs42;
        const double clhs44 =             DN(0,0)*N[1]*rho*tau1;
        const double clhs45 =             DN(1,1)*clhs20;
        const double clhs46 =             C(0,1)*DN(1,1) + clhs38;
        const double clhs47 =             C(1,2)*DN(1,1);
        const double clhs48 =             C(2,2)*DN(1,0) + clhs47;
        const double clhs49 =             N[0]*N[1]*rho;
        const double clhs50 =             clhs18*clhs49;
        const double clhs51 =             N[1]*clhs11*clhs13*tau1;
        const double clhs52 =             DN(0,1)*N[1]*rho*tau1;
        const double clhs53 =             DN(0,0)*N[1];
        const double clhs54 =             bdf0*clhs28*clhs29*tau2;
        const double clhs55 =             DN(1,0)*rho*tau1;
        const double clhs56 =             N[2]*clhs34;
        const double clhs57 =             DN(2,0)*clhs20 + clhs56;
        const double clhs58 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
        const double clhs59 =             C(0,2)*DN(2,0);
        const double clhs60 =             C(2,2)*DN(2,1) + clhs59;
        const double clhs61 =             DN(2,0)*clhs9 + DN(2,1)*clhs10;
        const double clhs62 =             N[2]*clhs8 + clhs61;
        const double clhs63 =             N[2]*bdf0;
        const double clhs64 =             clhs62 + clhs63;
        const double clhs65 =             DN(0,0)*N[2]*rho*tau1;
        const double clhs66 =             DN(2,1)*clhs20;
        const double clhs67 =             C(0,1)*DN(2,1) + clhs59;
        const double clhs68 =             C(1,2)*DN(2,1);
        const double clhs69 =             C(2,2)*DN(2,0) + clhs68;
        const double clhs70 =             N[0]*N[2]*rho;
        const double clhs71 =             clhs18*clhs70;
        const double clhs72 =             N[2]*clhs11*clhs13*tau1;
        const double clhs73 =             DN(0,1)*N[2]*rho*tau1;
        const double clhs74 =             DN(0,0)*N[2];
        const double clhs75 =             C(0,1)*DN(0,0) + clhs23;
        const double clhs76 =             DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1);
        const double clhs77 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
        const double clhs78 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + clhs10*clhs77 + clhs76*clhs9);
        const double clhs79 =             pow(DN(0,1), 2);
        const double clhs80 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
        const double clhs81 =             N[0]*clhs77 + clhs11;
        const double clhs82 =             clhs15 + clhs81;
        const double clhs83 =             DN(0,1)*tau2;
        const double clhs84 =             DN(1,0)*clhs83;
        const double clhs85 =             C(0,1)*DN(1,0) + clhs47;
        const double clhs86 =             clhs49*clhs76;
        const double clhs87 =             DN(1,1)*clhs83 + clhs35;
        const double clhs88 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
        const double clhs89 =             N[1]*clhs77 + clhs40;
        const double clhs90 =             clhs42 + clhs89;
        const double clhs91 =             DN(0,1)*N[1];
        const double clhs92 =             DN(1,1)*rho*tau1;
        const double clhs93 =             DN(2,0)*clhs83;
        const double clhs94 =             C(0,1)*DN(2,0) + clhs68;
        const double clhs95 =             clhs70*clhs76;
        const double clhs96 =             DN(2,1)*clhs83 + clhs56;
        const double clhs97 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
        const double clhs98 =             N[2]*clhs77 + clhs61;
        const double clhs99 =             clhs63 + clhs98;
        const double clhs100 =             DN(0,1)*N[2];
        const double clhs101 =             DN(0,0)*N[0];
        const double clhs102 =             DN(0,1)*N[0];
        const double clhs103 =             clhs76*rho*tau1;
        const double clhs104 =             DN(0,0)*rho*tau1;
        const double clhs105 =             clhs18*rho*tau1;
        const double clhs106 =             DN(0,1)*rho*tau1;
        const double clhs107 =             bdf0*clhs28*clhs29;
        const double clhs108 =             DN(1,0)*N[0];
        const double clhs109 =             DN(1,1)*N[0];
        const double clhs110 =             DN(0,0)*tau1;
        const double clhs111 =             DN(0,1)*tau1;
        const double clhs112 =             DN(1,0)*clhs110 + DN(1,1)*clhs111 + N[1]*clhs30;
        const double clhs113 =             DN(2,0)*N[0];
        const double clhs114 =             DN(2,1)*N[0];
        const double clhs115 =             DN(2,0)*clhs110 + DN(2,1)*clhs111 + N[2]*clhs30;
        const double clhs116 =             N[1]*rho;
        const double clhs117 =             clhs13*clhs40*tau1;
        const double clhs118 =             DN(1,0)*N[0]*rho*tau1;
        const double clhs119 =             N[0]*clhs13*clhs40*tau1;
        const double clhs120 =             DN(1,1)*N[0]*rho*tau1;
        const double clhs121 =             pow(DN(1,0), 2);
        const double clhs122 =             pow(N[1], 2);
        const double clhs123 =             clhs122*clhs2;
        const double clhs124 =             DN(1,0)*N[1]*rho*tau1;
        const double clhs125 =             DN(1,0)*tau2;
        const double clhs126 =             DN(1,1)*clhs125;
        const double clhs127 =             clhs122*rho;
        const double clhs128 =             N[1]*clhs13*clhs40*tau1;
        const double clhs129 =             DN(1,1)*N[1]*rho*tau1;
        const double clhs130 =             clhs28*clhs29*tau2;
        const double clhs131 =             clhs31*clhs40;
        const double clhs132 =             -N[1] + clhs130*clhs42 + clhs131;
        const double clhs133 =             N[1]*N[2]*bdf0;
        const double clhs134 =             clhs133*rho;
        const double clhs135 =             DN(2,0)*clhs125 + clhs134;
        const double clhs136 =             DN(1,0)*N[2]*rho*tau1;
        const double clhs137 =             DN(2,1)*clhs125;
        const double clhs138 =             N[1]*N[2]*rho;
        const double clhs139 =             clhs138*clhs18;
        const double clhs140 =             N[2]*clhs13*clhs40*tau1;
        const double clhs141 =             DN(1,1)*N[2]*rho*tau1;
        const double clhs142 =             DN(1,0)*N[2];
        const double clhs143 =             pow(DN(1,1), 2);
        const double clhs144 =             DN(1,1)*tau2;
        const double clhs145 =             DN(2,0)*clhs144;
        const double clhs146 =             clhs138*clhs76;
        const double clhs147 =             DN(2,1)*clhs144 + clhs134;
        const double clhs148 =             DN(1,1)*N[2];
        const double clhs149 =             DN(1,0)*N[1];
        const double clhs150 =             DN(1,1)*N[1];
        const double clhs151 =             DN(2,0)*N[1];
        const double clhs152 =             DN(2,1)*N[1];
        const double clhs153 =             DN(2,0)*(DN(1,0)*tau1) + DN(2,1)*(DN(1,1)*tau1) + clhs133*clhs28*clhs29;
        const double clhs154 =             N[2]*rho;
        const double clhs155 =             clhs13*clhs61*tau1;
        const double clhs156 =             DN(2,0)*N[0]*rho*tau1;
        const double clhs157 =             N[0]*clhs13*clhs61*tau1;
        const double clhs158 =             DN(2,1)*N[0]*rho*tau1;
        const double clhs159 =             DN(2,0)*N[1]*rho*tau1;
        const double clhs160 =             N[1]*clhs13*clhs61*tau1;
        const double clhs161 =             DN(2,1)*N[1]*rho*tau1;
        const double clhs162 =             pow(DN(2,0), 2);
        const double clhs163 =             pow(N[2], 2);
        const double clhs164 =             clhs163*clhs2;
        const double clhs165 =             DN(2,0)*N[2]*rho*tau1;
        const double clhs166 =             DN(2,0)*DN(2,1)*tau2;
        const double clhs167 =             clhs163*rho;
        const double clhs168 =             N[2]*clhs13*clhs61*tau1;
        const double clhs169 =             DN(2,1)*N[2]*rho*tau1;
        const double clhs170 =             -N[2] + clhs130*clhs63 + clhs31*clhs61;
        const double clhs171 =             pow(DN(2,1), 2);
        const double clhs172 =             DN(2,0)*rho*tau1;
        const double clhs173 =             DN(2,1)*rho*tau1;
        const double clhs174 =             DN(2,0)*N[2];
        const double clhs175 =             DN(2,1)*N[2];

        lhs(0,0)=DN(0,0)*clhs4 + DN(0,1)*clhs6 + clhs0*tau2 + clhs12*clhs7 + clhs14*clhs16 + clhs17*clhs19 + clhs3;
        lhs(0,1)=DN(0,0)*clhs22 + DN(0,1)*clhs24 + clhs18*clhs25 + clhs18*clhs26 + clhs19*clhs27 + clhs21;
        lhs(0,2)=DN(0,0)*clhs33;
        lhs(0,3)=DN(0,0)*clhs37 + DN(0,1)*clhs39 + clhs14*clhs43 + clhs19*clhs44 + clhs36 + clhs41*clhs7;
        lhs(0,4)=DN(0,0)*clhs46 + DN(0,1)*clhs48 + clhs18*clhs51 + clhs19*clhs52 + clhs45 + clhs50;
        lhs(0,5)=clhs11*clhs55 + clhs53*clhs54 - clhs53;
        lhs(0,6)=DN(0,0)*clhs58 + DN(0,1)*clhs60 + clhs14*clhs64 + clhs19*clhs65 + clhs57 + clhs62*clhs7;
        lhs(0,7)=DN(0,0)*clhs67 + DN(0,1)*clhs69 + clhs18*clhs72 + clhs19*clhs73 + clhs66 + clhs71;
        lhs(0,8)=DN(2,0)*clhs32 + clhs54*clhs74 - clhs74;
        lhs(1,0)=DN(0,0)*clhs6 + DN(0,1)*clhs75 + clhs17*clhs78 + clhs21 + clhs25*clhs76 + clhs26*clhs76;
        lhs(1,1)=DN(0,0)*clhs24 + DN(0,1)*clhs80 + clhs14*clhs82 + clhs27*clhs78 + clhs3 + clhs7*clhs81 + clhs79*tau2;
        lhs(1,2)=DN(0,1)*clhs33;
        lhs(1,3)=DN(0,0)*clhs39 + DN(0,1)*clhs85 + clhs44*clhs78 + clhs51*clhs76 + clhs84 + clhs86;
        lhs(1,4)=DN(0,0)*clhs48 + DN(0,1)*clhs88 + clhs14*clhs90 + clhs52*clhs78 + clhs7*clhs89 + clhs87;
        lhs(1,5)=clhs11*clhs92 + clhs54*clhs91 - clhs91;
        lhs(1,6)=DN(0,0)*clhs60 + DN(0,1)*clhs94 + clhs65*clhs78 + clhs72*clhs76 + clhs93 + clhs95;
        lhs(1,7)=DN(0,0)*clhs69 + DN(0,1)*clhs97 + clhs14*clhs99 + clhs7*clhs98 + clhs73*clhs78 + clhs96;
        lhs(1,8)=DN(2,1)*clhs32 + clhs100*clhs54 - clhs100;
        lhs(2,0)=clhs101 + clhs102*clhs103 + clhs104*clhs16;
        lhs(2,1)=clhs101*clhs105 + clhs102 + clhs106*clhs82;
        lhs(2,2)=clhs0*tau1 + clhs1*clhs107 + clhs79*tau1;
        lhs(2,3)=clhs103*clhs91 + clhs104*clhs43 + clhs108;
        lhs(2,4)=clhs105*clhs53 + clhs106*clhs90 + clhs109;
        lhs(2,5)=clhs112;
        lhs(2,6)=clhs100*clhs103 + clhs104*clhs64 + clhs113;
        lhs(2,7)=clhs105*clhs74 + clhs106*clhs99 + clhs114;
        lhs(2,8)=clhs115;
        lhs(3,0)=DN(1,0)*clhs4 + DN(1,1)*clhs6 + clhs116*clhs12 + clhs117*clhs16 + clhs118*clhs19 + clhs36;
        lhs(3,1)=DN(1,0)*clhs22 + DN(1,1)*clhs24 + clhs119*clhs18 + clhs120*clhs19 + clhs50 + clhs84;
        lhs(3,2)=clhs104*clhs40 + clhs108*clhs54 - clhs108;
        lhs(3,3)=DN(1,0)*clhs37 + DN(1,1)*clhs39 + clhs116*clhs41 + clhs117*clhs43 + clhs121*tau2 + clhs123 + clhs124*clhs19;
        lhs(3,4)=DN(1,0)*clhs46 + DN(1,1)*clhs48 + clhs126 + clhs127*clhs18 + clhs128*clhs18 + clhs129*clhs19;
        lhs(3,5)=DN(1,0)*clhs132;
        lhs(3,6)=DN(1,0)*clhs58 + DN(1,1)*clhs60 + clhs116*clhs62 + clhs117*clhs64 + clhs135 + clhs136*clhs19;
        lhs(3,7)=DN(1,0)*clhs67 + DN(1,1)*clhs69 + clhs137 + clhs139 + clhs140*clhs18 + clhs141*clhs19;
        lhs(3,8)=DN(2,0)*clhs131 + clhs142*clhs54 - clhs142;
        lhs(4,0)=DN(1,0)*clhs6 + DN(1,1)*clhs75 + clhs118*clhs78 + clhs119*clhs76 + clhs45 + clhs86;
        lhs(4,1)=DN(1,0)*clhs24 + DN(1,1)*clhs80 + clhs116*clhs81 + clhs117*clhs82 + clhs120*clhs78 + clhs87;
        lhs(4,2)=clhs106*clhs40 + clhs109*clhs54 - clhs109;
        lhs(4,3)=DN(1,0)*clhs39 + DN(1,1)*clhs85 + clhs124*clhs78 + clhs126 + clhs127*clhs76 + clhs128*clhs76;
        lhs(4,4)=DN(1,0)*clhs48 + DN(1,1)*clhs88 + clhs116*clhs89 + clhs117*clhs90 + clhs123 + clhs129*clhs78 + clhs143*tau2;
        lhs(4,5)=DN(1,1)*clhs132;
        lhs(4,6)=DN(1,0)*clhs60 + DN(1,1)*clhs94 + clhs136*clhs78 + clhs140*clhs76 + clhs145 + clhs146;
        lhs(4,7)=DN(1,0)*clhs69 + DN(1,1)*clhs97 + clhs116*clhs98 + clhs117*clhs99 + clhs141*clhs78 + clhs147;
        lhs(4,8)=DN(2,1)*clhs131 + clhs148*clhs54 - clhs148;
        lhs(5,0)=clhs103*clhs109 + clhs16*clhs55 + clhs53;
        lhs(5,1)=clhs105*clhs108 + clhs82*clhs92 + clhs91;
        lhs(5,2)=clhs112;
        lhs(5,3)=clhs103*clhs150 + clhs149 + clhs43*clhs55;
        lhs(5,4)=clhs105*clhs149 + clhs150 + clhs90*clhs92;
        lhs(5,5)=clhs107*clhs122 + clhs121*tau1 + clhs143*tau1;
        lhs(5,6)=clhs103*clhs148 + clhs151 + clhs55*clhs64;
        lhs(5,7)=clhs105*clhs142 + clhs152 + clhs92*clhs99;
        lhs(5,8)=clhs153;
        lhs(6,0)=DN(2,0)*clhs4 + DN(2,1)*clhs6 + clhs12*clhs154 + clhs155*clhs16 + clhs156*clhs19 + clhs57;
        lhs(6,1)=DN(2,0)*clhs22 + DN(2,1)*clhs24 + clhs157*clhs18 + clhs158*clhs19 + clhs71 + clhs93;
        lhs(6,2)=clhs104*clhs61 + clhs113*clhs54 - clhs113;
        lhs(6,3)=DN(2,0)*clhs37 + DN(2,1)*clhs39 + clhs135 + clhs154*clhs41 + clhs155*clhs43 + clhs159*clhs19;
        lhs(6,4)=DN(2,0)*clhs46 + DN(2,1)*clhs48 + clhs139 + clhs145 + clhs160*clhs18 + clhs161*clhs19;
        lhs(6,5)=clhs151*clhs54 - clhs151 + clhs55*clhs61;
        lhs(6,6)=DN(2,0)*clhs58 + DN(2,1)*clhs60 + clhs154*clhs62 + clhs155*clhs64 + clhs162*tau2 + clhs164 + clhs165*clhs19;
        lhs(6,7)=DN(2,0)*clhs67 + DN(2,1)*clhs69 + clhs166 + clhs167*clhs18 + clhs168*clhs18 + clhs169*clhs19;
        lhs(6,8)=DN(2,0)*clhs170;
        lhs(7,0)=DN(2,0)*clhs6 + DN(2,1)*clhs75 + clhs156*clhs78 + clhs157*clhs76 + clhs66 + clhs95;
        lhs(7,1)=DN(2,0)*clhs24 + DN(2,1)*clhs80 + clhs154*clhs81 + clhs155*clhs82 + clhs158*clhs78 + clhs96;
        lhs(7,2)=clhs106*clhs61 + clhs114*clhs54 - clhs114;
        lhs(7,3)=DN(2,0)*clhs39 + DN(2,1)*clhs85 + clhs137 + clhs146 + clhs159*clhs78 + clhs160*clhs76;
        lhs(7,4)=DN(2,0)*clhs48 + DN(2,1)*clhs88 + clhs147 + clhs154*clhs89 + clhs155*clhs90 + clhs161*clhs78;
        lhs(7,5)=clhs152*clhs54 - clhs152 + clhs61*clhs92;
        lhs(7,6)=DN(2,0)*clhs60 + DN(2,1)*clhs94 + clhs165*clhs78 + clhs166 + clhs167*clhs76 + clhs168*clhs76;
        lhs(7,7)=DN(2,0)*clhs69 + DN(2,1)*clhs97 + clhs154*clhs98 + clhs155*clhs99 + clhs164 + clhs169*clhs78 + clhs171*tau2;
        lhs(7,8)=DN(2,1)*clhs170;
        lhs(8,0)=clhs103*clhs114 + clhs16*clhs172 + clhs74;
        lhs(8,1)=clhs100 + clhs105*clhs113 + clhs173*clhs82;
        lhs(8,2)=clhs115;
        lhs(8,3)=clhs103*clhs152 + clhs142 + clhs172*clhs43;
        lhs(8,4)=clhs105*clhs151 + clhs148 + clhs173*clhs90;
        lhs(8,5)=clhs153;
        lhs(8,6)=clhs103*clhs175 + clhs172*clhs64 + clhs174;
        lhs(8,7)=clhs105*clhs174 + clhs173*clhs99 + clhs175;
        lhs(8,8)=clhs107*clhs163 + clhs162*tau1 + clhs171*tau1;

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
        const double crhs16 =             N[0]*(v(0,0) - vmesh(0,0)) + N[1]*(v(1,0) - vmesh(1,0)) + N[2]*(v(2,0) - vmesh(2,0)) + N[3]*(v(3,0) - vmesh(3,0));
        const double crhs17 =             N[0]*(v(0,1) - vmesh(0,1)) + N[1]*(v(1,1) - vmesh(1,1)) + N[2]*(v(2,1) - vmesh(2,1)) + N[3]*(v(3,1) - vmesh(3,1));
        const double crhs18 =             N[0]*(v(0,2) - vmesh(0,2)) + N[1]*(v(1,2) - vmesh(1,2)) + N[2]*(v(2,2) - vmesh(2,2)) + N[3]*(v(3,2) - vmesh(3,2));
        const double crhs19 =             crhs16*(crhs11 + crhs5 + crhs7 + crhs9) + crhs17*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs18*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0));
        const double crhs20 =             rho*(DN(0,0)*crhs16 + DN(0,1)*crhs17 + DN(0,2)*crhs18);
        const double crhs21 =             tau1*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + rho*(crhs19 + crhs3));
        const double crhs22 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
        const double crhs23 =             N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1));
        const double crhs24 =             crhs16*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs17*(crhs10 + crhs12 + crhs6 + crhs8) + crhs18*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1));
        const double crhs25 =             tau1*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs22 + rho*(crhs23 + crhs24));
        const double crhs26 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
        const double crhs27 =             N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2));
        const double crhs28 =             crhs16*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs17*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + crhs18*crhs4;
        const double crhs29 =             tau1*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs26 + rho*(crhs27 + crhs28));
        const double crhs30 =             N[1]*rho;
        const double crhs31 =             rho*(DN(1,0)*crhs16 + DN(1,1)*crhs17 + DN(1,2)*crhs18);
        const double crhs32 =             N[2]*rho;
        const double crhs33 =             rho*(DN(2,0)*crhs16 + DN(2,1)*crhs17 + DN(2,2)*crhs18);
        const double crhs34 =             N[3]*rho;
        const double crhs35 =             rho*(DN(3,0)*crhs16 + DN(3,1)*crhs17 + DN(3,2)*crhs18);

        rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs15 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - crhs19*crhs2 - crhs2*crhs3 - crhs20*crhs21;
        rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs15 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs22 - crhs2*crhs23 - crhs2*crhs24 - crhs20*crhs25;
        rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs15 - DN(0,2)*stress[2] + N[0]*crhs26 - crhs2*crhs27 - crhs2*crhs28 - crhs20*crhs29;
        rhs[3]=-DN(0,0)*crhs21 - DN(0,1)*crhs25 - DN(0,2)*crhs29 - N[0]*crhs13 - N[0]*crhs14;
        rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs15 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - crhs19*crhs30 - crhs21*crhs31 - crhs3*crhs30;
        rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs15 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs22 - crhs23*crhs30 - crhs24*crhs30 - crhs25*crhs31;
        rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs15 - DN(1,2)*stress[2] + N[1]*crhs26 - crhs27*crhs30 - crhs28*crhs30 - crhs29*crhs31;
        rhs[7]=-DN(1,0)*crhs21 - DN(1,1)*crhs25 - DN(1,2)*crhs29 - N[1]*crhs13 - N[1]*crhs14;
        rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs15 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - crhs19*crhs32 - crhs21*crhs33 - crhs3*crhs32;
        rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs15 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs22 - crhs23*crhs32 - crhs24*crhs32 - crhs25*crhs33;
        rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs15 - DN(2,2)*stress[2] + N[2]*crhs26 - crhs27*crhs32 - crhs28*crhs32 - crhs29*crhs33;
        rhs[11]=-DN(2,0)*crhs21 - DN(2,1)*crhs25 - DN(2,2)*crhs29 - N[2]*crhs13 - N[2]*crhs14;
        rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs15 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - crhs19*crhs34 - crhs21*crhs35 - crhs3*crhs34;
        rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs15 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs22 - crhs23*crhs34 - crhs24*crhs34 - crhs25*crhs35;
        rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs15 - DN(3,2)*stress[2] + N[3]*crhs26 - crhs27*crhs34 - crhs28*crhs34 - crhs29*crhs35;
        rhs[15]=-DN(3,0)*crhs21 - DN(3,1)*crhs25 - DN(3,2)*crhs29 - N[3]*crhs13 - N[3]*crhs14;

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
        const double crhs9 =             N[0]*(v(0,0) - vmesh(0,0)) + N[1]*(v(1,0) - vmesh(1,0)) + N[2]*(v(2,0) - vmesh(2,0));
        const double crhs10 =             N[0]*(v(0,1) - vmesh(0,1)) + N[1]*(v(1,1) - vmesh(1,1)) + N[2]*(v(2,1) - vmesh(2,1));
        const double crhs11 =             crhs10*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)) + crhs4*crhs9;
        const double crhs12 =             rho*(DN(0,0)*crhs9 + DN(0,1)*crhs10);
        const double crhs13 =             tau1*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + rho*(crhs11 + crhs3));
        const double crhs14 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
        const double crhs15 =             N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1));
        const double crhs16 =             crhs10*crhs5 + crhs9*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1));
        const double crhs17 =             tau1*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs14 + rho*(crhs15 + crhs16));
        const double crhs18 =             N[1]*rho;
        const double crhs19 =             rho*(DN(1,0)*crhs9 + DN(1,1)*crhs10);
        const double crhs20 =             N[2]*rho;
        const double crhs21 =             rho*(DN(2,0)*crhs9 + DN(2,1)*crhs10);

        rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs8 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - crhs11*crhs2 - crhs12*crhs13 - crhs2*crhs3;
        rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs8 - DN(0,1)*stress[1] + N[0]*crhs14 - crhs12*crhs17 - crhs15*crhs2 - crhs16*crhs2;
        rhs[2]=-DN(0,0)*crhs13 - DN(0,1)*crhs17 - N[0]*crhs6 - N[0]*crhs7;
        rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs8 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - crhs11*crhs18 - crhs13*crhs19 - crhs18*crhs3;
        rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs8 - DN(1,1)*stress[1] + N[1]*crhs14 - crhs15*crhs18 - crhs16*crhs18 - crhs17*crhs19;
        rhs[5]=-DN(1,0)*crhs13 - DN(1,1)*crhs17 - N[1]*crhs6 - N[1]*crhs7;
        rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs8 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - crhs11*crhs20 - crhs13*crhs21 - crhs20*crhs3;
        rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs8 - DN(2,1)*stress[1] + N[2]*crhs14 - crhs15*crhs20 - crhs16*crhs20 - crhs17*crhs21;
        rhs[8]=-DN(2,0)*crhs13 - DN(2,1)*crhs17 - N[2]*crhs6 - N[2]*crhs7;

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
