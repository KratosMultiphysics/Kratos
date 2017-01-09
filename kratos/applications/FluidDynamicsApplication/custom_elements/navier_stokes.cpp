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
        const double clhs1 =             DN(0,0)*rho*tau1;
        const double clhs2 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
        const double clhs3 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
        const double clhs4 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
        const double clhs5 =             DN(0,0)*clhs2 + DN(0,1)*clhs3 + DN(0,2)*clhs4;
        const double clhs6 =             pow(c, -2);
        const double clhs7 =             N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]) + N[3]*(bdf0*p[3] + bdf1*pn[3] + bdf2*pnn[3]);
        const double clhs8 =             clhs6*clhs7;
        const double clhs9 =             N[0]*clhs8;
        const double clhs10 =             clhs9 + rho*(N[0]*bdf0 + clhs5);
        const double clhs11 =             pow(N[0], 2);
        const double clhs12 =             bdf0*rho;
        const double clhs13 =             N[0]*rho;
        const double clhs14 =             N[0]*clhs6*clhs7*tau1;
        const double clhs15 =             -clhs10*clhs14 + clhs11*clhs12 + clhs11*clhs8 + clhs13*clhs5;
        const double clhs16 =             DN(0,0)*rho;
        const double clhs17 =             DN(0,0)*tau2;
        const double clhs18 =             rho*tau1;
        const double clhs19 =             clhs10*clhs18;
        const double clhs20 =             clhs16 - clhs17 + clhs19*clhs2;
        const double clhs21 =             DN(0,0)*N[0];
        const double clhs22 =             1.0/rho;
        const double clhs23 =             bdf0*clhs22*clhs6*tau2;
        const double clhs24 =             bdf0*clhs11*clhs6;
        const double clhs25 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0);
        const double clhs26 =             N[0]*bdf0*clhs6;
        const double clhs27 =             clhs25*clhs26;
        const double clhs28 =             DN(0,0) + clhs27;
        const double clhs29 =             clhs1*clhs28;
        const double clhs30 =             DN(0,1)*rho*tau1;
        const double clhs31 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
        const double clhs32 =             clhs26*clhs31;
        const double clhs33 =             DN(0,1) + clhs32;
        const double clhs34 =             clhs30*clhs33;
        const double clhs35 =             DN(0,2)*rho*tau1;
        const double clhs36 =             N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2);
        const double clhs37 =             clhs26*clhs36;
        const double clhs38 =             DN(0,2) + clhs37;
        const double clhs39 =             clhs35*clhs38;
        const double clhs40 =             bdf0*clhs11*clhs6*tau1;
        const double clhs41 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + clhs25*clhs8 - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)) + clhs2*(DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0)) + clhs3*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + clhs4*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)));
        const double clhs42 =             DN(1,0)*clhs16;
        const double clhs43 =             -DN(1,0)*clhs17;
        const double clhs44 =             N[0]*bdf0*rho;
        const double clhs45 =             N[1]*clhs44;
        const double clhs46 =             N[1]*clhs9;
        const double clhs47 =             DN(1,0)*clhs2 + DN(1,1)*clhs3 + DN(1,2)*clhs4;
        const double clhs48 =             N[1]*clhs8;
        const double clhs49 =             clhs48 + rho*(N[1]*bdf0 + clhs47);
        const double clhs50 =             clhs1*clhs49;
        const double clhs51 =             clhs13*clhs47;
        const double clhs52 =             -clhs14*clhs49;
        const double clhs53 =             DN(1,1)*clhs16 - DN(1,1)*clhs17;
        const double clhs54 =             clhs30*clhs49;
        const double clhs55 =             DN(1,2)*clhs16 - DN(1,2)*clhs17;
        const double clhs56 =             clhs35*clhs49;
        const double clhs57 =             DN(0,0)*N[1];
        const double clhs58 =             N[1]*bdf0*clhs6;
        const double clhs59 =             clhs25*clhs58;
        const double clhs60 =             DN(1,0) + clhs59;
        const double clhs61 =             clhs1*clhs60;
        const double clhs62 =             clhs31*clhs58;
        const double clhs63 =             DN(1,1) + clhs62;
        const double clhs64 =             clhs30*clhs63;
        const double clhs65 =             clhs36*clhs58;
        const double clhs66 =             DN(1,2) + clhs65;
        const double clhs67 =             clhs35*clhs66;
        const double clhs68 =             N[0]*N[1]*bdf0*clhs6*tau1;
        const double clhs69 =             N[1]*clhs27 - clhs41*clhs68;
        const double clhs70 =             DN(2,0)*clhs16;
        const double clhs71 =             -DN(2,0)*clhs17;
        const double clhs72 =             N[2]*clhs44;
        const double clhs73 =             N[2]*clhs9;
        const double clhs74 =             DN(2,0)*clhs2 + DN(2,1)*clhs3 + DN(2,2)*clhs4;
        const double clhs75 =             N[2]*clhs8;
        const double clhs76 =             clhs75 + rho*(N[2]*bdf0 + clhs74);
        const double clhs77 =             clhs1*clhs76;
        const double clhs78 =             clhs13*clhs74;
        const double clhs79 =             -clhs14*clhs76;
        const double clhs80 =             DN(2,1)*clhs16 - DN(2,1)*clhs17;
        const double clhs81 =             clhs30*clhs76;
        const double clhs82 =             DN(2,2)*clhs16 - DN(2,2)*clhs17;
        const double clhs83 =             clhs35*clhs76;
        const double clhs84 =             DN(0,0)*N[2];
        const double clhs85 =             N[2]*bdf0*clhs6;
        const double clhs86 =             DN(2,0) + clhs25*clhs85;
        const double clhs87 =             clhs1*clhs86;
        const double clhs88 =             DN(2,1) + clhs31*clhs85;
        const double clhs89 =             clhs30*clhs88;
        const double clhs90 =             DN(2,2) + clhs36*clhs85;
        const double clhs91 =             clhs35*clhs90;
        const double clhs92 =             N[0]*N[2]*bdf0*clhs6*tau1;
        const double clhs93 =             N[2]*clhs27 - clhs41*clhs92;
        const double clhs94 =             DN(3,0)*clhs16;
        const double clhs95 =             -DN(3,0)*clhs17;
        const double clhs96 =             N[3]*clhs44;
        const double clhs97 =             N[3]*clhs9;
        const double clhs98 =             DN(3,0)*clhs2 + DN(3,1)*clhs3 + DN(3,2)*clhs4;
        const double clhs99 =             N[3]*clhs8 + rho*(N[3]*bdf0 + clhs98);
        const double clhs100 =             clhs1*clhs99;
        const double clhs101 =             clhs13*clhs98;
        const double clhs102 =             -clhs14*clhs99;
        const double clhs103 =             DN(3,1)*clhs16 - DN(3,1)*clhs17;
        const double clhs104 =             clhs30*clhs99;
        const double clhs105 =             DN(3,2)*clhs16 - DN(3,2)*clhs17;
        const double clhs106 =             clhs35*clhs99;
        const double clhs107 =             DN(0,0)*N[3];
        const double clhs108 =             N[3]*bdf0*clhs6;
        const double clhs109 =             DN(3,0) + clhs108*clhs25;
        const double clhs110 =             clhs1*clhs109;
        const double clhs111 =             DN(3,1) + clhs108*clhs31;
        const double clhs112 =             clhs111*clhs30;
        const double clhs113 =             DN(3,2) + clhs108*clhs36;
        const double clhs114 =             clhs113*clhs35;
        const double clhs115 =             N[0]*N[3]*bdf0*clhs6*tau1;
        const double clhs116 =             N[3]*clhs27 - clhs115*clhs41;
        const double clhs117 =             DN(0,1)*rho;
        const double clhs118 =             DN(0,1)*tau2;
        const double clhs119 =             clhs117 - clhs118 + clhs19*clhs3;
        const double clhs120 =             pow(DN(0,1), 2);
        const double clhs121 =             DN(0,1)*N[0];
        const double clhs122 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + clhs31*clhs8 - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)) + clhs2*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + clhs3*(DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1)) + clhs4*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
        const double clhs123 =             DN(1,0)*clhs117 - DN(1,0)*clhs118;
        const double clhs124 =             DN(1,1)*clhs117 - DN(1,1)*clhs118 + clhs45 + clhs46;
        const double clhs125 =             DN(1,2)*clhs117 - DN(1,2)*clhs118;
        const double clhs126 =             DN(0,1)*N[1];
        const double clhs127 =             N[1]*clhs32 - clhs122*clhs68;
        const double clhs128 =             DN(2,0)*clhs117 - DN(2,0)*clhs118;
        const double clhs129 =             DN(2,1)*clhs117 - DN(2,1)*clhs118 + clhs72 + clhs73;
        const double clhs130 =             DN(2,2)*clhs117 - DN(2,2)*clhs118;
        const double clhs131 =             DN(0,1)*N[2];
        const double clhs132 =             N[2]*clhs32 - clhs122*clhs92;
        const double clhs133 =             DN(3,0)*clhs117 - DN(3,0)*clhs118;
        const double clhs134 =             DN(3,1)*clhs117 - DN(3,1)*clhs118 + clhs96 + clhs97;
        const double clhs135 =             DN(3,2)*clhs117 - DN(3,2)*clhs118;
        const double clhs136 =             DN(0,1)*N[3];
        const double clhs137 =             N[3]*clhs32 - clhs115*clhs122;
        const double clhs138 =             DN(0,2)*rho;
        const double clhs139 =             DN(0,2)*tau2;
        const double clhs140 =             clhs138 - clhs139 + clhs19*clhs4;
        const double clhs141 =             pow(DN(0,2), 2);
        const double clhs142 =             DN(0,2)*N[0];
        const double clhs143 =             DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + clhs36*clhs8 - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)) + clhs2*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + clhs3*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + clhs4*(DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2)));
        const double clhs144 =             DN(1,0)*clhs138 - DN(1,0)*clhs139;
        const double clhs145 =             DN(1,1)*clhs138 - DN(1,1)*clhs139;
        const double clhs146 =             DN(1,2)*clhs138 - DN(1,2)*clhs139 + clhs45 + clhs46;
        const double clhs147 =             DN(0,2)*N[1];
        const double clhs148 =             N[1]*clhs37 - clhs143*clhs68;
        const double clhs149 =             DN(2,0)*clhs138 - DN(2,0)*clhs139;
        const double clhs150 =             DN(2,1)*clhs138 - DN(2,1)*clhs139;
        const double clhs151 =             DN(2,2)*clhs138 - DN(2,2)*clhs139 + clhs72 + clhs73;
        const double clhs152 =             DN(0,2)*N[2];
        const double clhs153 =             N[2]*clhs37 - clhs143*clhs92;
        const double clhs154 =             DN(3,0)*clhs138 - DN(3,0)*clhs139;
        const double clhs155 =             DN(3,1)*clhs138 - DN(3,1)*clhs139;
        const double clhs156 =             DN(3,2)*clhs138 - DN(3,2)*clhs139 + clhs96 + clhs97;
        const double clhs157 =             DN(0,2)*N[3];
        const double clhs158 =             N[3]*clhs37 - clhs115*clhs143;
        const double clhs159 =             2*N[0];
        const double clhs160 =             clhs159 + clhs19;
        const double clhs161 =             bdf0*clhs22*clhs6;
        const double clhs162 =             N[0]*bdf0*clhs22*clhs6;
        const double clhs163 =             N[1]*clhs162;
        const double clhs164 =             N[2]*clhs162;
        const double clhs165 =             N[3]*clhs162;
        const double clhs166 =             DN(1,0)*rho*tau1;
        const double clhs167 =             clhs10*clhs166;
        const double clhs168 =             N[1]*rho;
        const double clhs169 =             clhs168*clhs5;
        const double clhs170 =             N[1]*clhs6*clhs7*tau1;
        const double clhs171 =             -clhs10*clhs170;
        const double clhs172 =             DN(1,1)*rho*tau1;
        const double clhs173 =             clhs10*clhs172;
        const double clhs174 =             DN(1,2)*rho*tau1;
        const double clhs175 =             clhs10*clhs174;
        const double clhs176 =             DN(1,0)*N[0];
        const double clhs177 =             clhs166*clhs28;
        const double clhs178 =             clhs172*clhs33;
        const double clhs179 =             clhs174*clhs38;
        const double clhs180 =             pow(DN(1,0), 2);
        const double clhs181 =             pow(N[1], 2);
        const double clhs182 =             clhs12*clhs181 + clhs168*clhs47 - clhs170*clhs49 + clhs181*clhs8;
        const double clhs183 =             DN(1,0)*rho;
        const double clhs184 =             DN(1,0)*tau2;
        const double clhs185 =             clhs18*clhs49;
        const double clhs186 =             clhs183 - clhs184 + clhs185*clhs2;
        const double clhs187 =             DN(1,0)*N[1];
        const double clhs188 =             bdf0*clhs181*clhs6;
        const double clhs189 =             clhs166*clhs60;
        const double clhs190 =             clhs172*clhs63;
        const double clhs191 =             clhs174*clhs66;
        const double clhs192 =             bdf0*clhs181*clhs6*tau1;
        const double clhs193 =             DN(2,0)*clhs183;
        const double clhs194 =             -DN(2,0)*clhs184;
        const double clhs195 =             N[1]*bdf0*rho;
        const double clhs196 =             N[2]*clhs195;
        const double clhs197 =             N[2]*clhs48;
        const double clhs198 =             clhs166*clhs76;
        const double clhs199 =             clhs168*clhs74;
        const double clhs200 =             -clhs170*clhs76;
        const double clhs201 =             DN(2,1)*clhs183 - DN(2,1)*clhs184;
        const double clhs202 =             clhs172*clhs76;
        const double clhs203 =             DN(2,2)*clhs183 - DN(2,2)*clhs184;
        const double clhs204 =             clhs174*clhs76;
        const double clhs205 =             DN(1,0)*N[2];
        const double clhs206 =             clhs166*clhs86;
        const double clhs207 =             clhs172*clhs88;
        const double clhs208 =             clhs174*clhs90;
        const double clhs209 =             N[1]*N[2]*bdf0*clhs6*tau1;
        const double clhs210 =             N[2]*clhs59 - clhs209*clhs41;
        const double clhs211 =             DN(3,0)*clhs183;
        const double clhs212 =             -DN(3,0)*clhs184;
        const double clhs213 =             N[3]*clhs195;
        const double clhs214 =             N[3]*clhs48;
        const double clhs215 =             clhs166*clhs99;
        const double clhs216 =             clhs168*clhs98;
        const double clhs217 =             -clhs170*clhs99;
        const double clhs218 =             DN(3,1)*clhs183 - DN(3,1)*clhs184;
        const double clhs219 =             clhs172*clhs99;
        const double clhs220 =             DN(3,2)*clhs183 - DN(3,2)*clhs184;
        const double clhs221 =             clhs174*clhs99;
        const double clhs222 =             DN(1,0)*N[3];
        const double clhs223 =             clhs109*clhs166;
        const double clhs224 =             clhs111*clhs172;
        const double clhs225 =             clhs113*clhs174;
        const double clhs226 =             N[1]*N[3]*bdf0*clhs6*tau1;
        const double clhs227 =             N[3]*clhs59 - clhs226*clhs41;
        const double clhs228 =             DN(1,1)*N[0];
        const double clhs229 =             DN(1,1)*rho;
        const double clhs230 =             DN(1,1)*tau2;
        const double clhs231 =             clhs185*clhs3 + clhs229 - clhs230;
        const double clhs232 =             pow(DN(1,1), 2);
        const double clhs233 =             DN(1,1)*N[1];
        const double clhs234 =             DN(2,0)*clhs229 - DN(2,0)*clhs230;
        const double clhs235 =             DN(2,1)*clhs229 - DN(2,1)*clhs230 + clhs196 + clhs197;
        const double clhs236 =             DN(2,2)*clhs229 - DN(2,2)*clhs230;
        const double clhs237 =             DN(1,1)*N[2];
        const double clhs238 =             N[2]*clhs62 - clhs122*clhs209;
        const double clhs239 =             DN(3,0)*clhs229 - DN(3,0)*clhs230;
        const double clhs240 =             DN(3,1)*clhs229 - DN(3,1)*clhs230 + clhs213 + clhs214;
        const double clhs241 =             DN(3,2)*clhs229 - DN(3,2)*clhs230;
        const double clhs242 =             DN(1,1)*N[3];
        const double clhs243 =             N[3]*clhs62 - clhs122*clhs226;
        const double clhs244 =             DN(1,2)*N[0];
        const double clhs245 =             DN(1,2)*rho;
        const double clhs246 =             DN(1,2)*tau2;
        const double clhs247 =             clhs185*clhs4 + clhs245 - clhs246;
        const double clhs248 =             pow(DN(1,2), 2);
        const double clhs249 =             DN(1,2)*N[1];
        const double clhs250 =             DN(2,0)*clhs245 - DN(2,0)*clhs246;
        const double clhs251 =             DN(2,1)*clhs245 - DN(2,1)*clhs246;
        const double clhs252 =             DN(2,2)*clhs245 - DN(2,2)*clhs246 + clhs196 + clhs197;
        const double clhs253 =             DN(1,2)*N[2];
        const double clhs254 =             N[2]*clhs65 - clhs143*clhs209;
        const double clhs255 =             DN(3,0)*clhs245 - DN(3,0)*clhs246;
        const double clhs256 =             DN(3,1)*clhs245 - DN(3,1)*clhs246;
        const double clhs257 =             DN(3,2)*clhs245 - DN(3,2)*clhs246 + clhs213 + clhs214;
        const double clhs258 =             DN(1,2)*N[3];
        const double clhs259 =             N[3]*clhs65 - clhs143*clhs226;
        const double clhs260 =             2*N[1];
        const double clhs261 =             clhs185 + clhs260;
        const double clhs262 =             N[1]*bdf0*clhs22*clhs6;
        const double clhs263 =             N[2]*clhs262;
        const double clhs264 =             N[3]*clhs262;
        const double clhs265 =             DN(2,0)*rho*tau1;
        const double clhs266 =             clhs10*clhs265;
        const double clhs267 =             N[2]*rho;
        const double clhs268 =             clhs267*clhs5;
        const double clhs269 =             N[2]*clhs6*clhs7*tau1;
        const double clhs270 =             -clhs10*clhs269;
        const double clhs271 =             DN(2,1)*rho*tau1;
        const double clhs272 =             clhs10*clhs271;
        const double clhs273 =             DN(2,2)*rho*tau1;
        const double clhs274 =             clhs10*clhs273;
        const double clhs275 =             DN(2,0)*N[0];
        const double clhs276 =             clhs265*clhs28;
        const double clhs277 =             clhs271*clhs33;
        const double clhs278 =             clhs273*clhs38;
        const double clhs279 =             clhs265*clhs49;
        const double clhs280 =             clhs267*clhs47;
        const double clhs281 =             -clhs269*clhs49;
        const double clhs282 =             clhs271*clhs49;
        const double clhs283 =             clhs273*clhs49;
        const double clhs284 =             DN(2,0)*N[1];
        const double clhs285 =             clhs265*clhs60;
        const double clhs286 =             clhs271*clhs63;
        const double clhs287 =             clhs273*clhs66;
        const double clhs288 =             pow(DN(2,0), 2);
        const double clhs289 =             pow(N[2], 2);
        const double clhs290 =             clhs12*clhs289 + clhs267*clhs74 - clhs269*clhs76 + clhs289*clhs8;
        const double clhs291 =             DN(2,0)*rho;
        const double clhs292 =             DN(2,0)*tau2;
        const double clhs293 =             clhs18*clhs76;
        const double clhs294 =             clhs2*clhs293 + clhs291 - clhs292;
        const double clhs295 =             DN(2,0)*N[2];
        const double clhs296 =             bdf0*clhs289*clhs6;
        const double clhs297 =             clhs265*clhs86;
        const double clhs298 =             clhs271*clhs88;
        const double clhs299 =             clhs273*clhs90;
        const double clhs300 =             bdf0*clhs289*clhs6*tau1;
        const double clhs301 =             DN(3,0)*clhs291;
        const double clhs302 =             -DN(3,0)*clhs292;
        const double clhs303 =             N[2]*N[3]*bdf0;
        const double clhs304 =             clhs303*rho;
        const double clhs305 =             N[3]*clhs75;
        const double clhs306 =             clhs265*clhs99;
        const double clhs307 =             clhs267*clhs98;
        const double clhs308 =             -clhs269*clhs99;
        const double clhs309 =             DN(3,1)*clhs291 - DN(3,1)*clhs292;
        const double clhs310 =             clhs271*clhs99;
        const double clhs311 =             DN(3,2)*clhs291 - DN(3,2)*clhs292;
        const double clhs312 =             clhs273*clhs99;
        const double clhs313 =             DN(2,0)*N[3];
        const double clhs314 =             clhs109*clhs265;
        const double clhs315 =             clhs111*clhs271;
        const double clhs316 =             clhs113*clhs273;
        const double clhs317 =             N[2]*N[3]*bdf0*clhs6;
        const double clhs318 =             N[2]*N[3]*bdf0*clhs6*tau1;
        const double clhs319 =             clhs25*clhs317 - clhs318*clhs41;
        const double clhs320 =             DN(2,1)*N[0];
        const double clhs321 =             DN(2,1)*N[1];
        const double clhs322 =             DN(2,1)*rho;
        const double clhs323 =             DN(2,1)*tau2;
        const double clhs324 =             clhs293*clhs3 + clhs322 - clhs323;
        const double clhs325 =             pow(DN(2,1), 2);
        const double clhs326 =             DN(2,1)*N[2];
        const double clhs327 =             DN(3,0)*clhs322 - DN(3,0)*clhs323;
        const double clhs328 =             DN(3,1)*clhs322 - DN(3,1)*clhs323 + clhs304 + clhs305;
        const double clhs329 =             DN(3,2)*clhs322 - DN(3,2)*clhs323;
        const double clhs330 =             DN(2,1)*N[3];
        const double clhs331 =             -clhs122*clhs318 + clhs31*clhs317;
        const double clhs332 =             DN(2,2)*N[0];
        const double clhs333 =             DN(2,2)*N[1];
        const double clhs334 =             DN(2,2)*rho;
        const double clhs335 =             DN(2,2)*tau2;
        const double clhs336 =             clhs293*clhs4 + clhs334 - clhs335;
        const double clhs337 =             pow(DN(2,2), 2);
        const double clhs338 =             DN(2,2)*N[2];
        const double clhs339 =             DN(3,0)*clhs334 - DN(3,0)*clhs335;
        const double clhs340 =             DN(3,1)*clhs334 - DN(3,1)*clhs335;
        const double clhs341 =             DN(3,2)*clhs334 - DN(3,2)*clhs335 + clhs304 + clhs305;
        const double clhs342 =             DN(2,2)*N[3];
        const double clhs343 =             -clhs143*clhs318 + clhs317*clhs36;
        const double clhs344 =             2*N[2];
        const double clhs345 =             clhs293 + clhs344;
        const double clhs346 =             clhs22*clhs303*clhs6;
        const double clhs347 =             DN(3,0)*rho*tau1;
        const double clhs348 =             clhs10*clhs347;
        const double clhs349 =             N[3]*rho;
        const double clhs350 =             clhs349*clhs5;
        const double clhs351 =             N[3]*clhs6*clhs7*tau1;
        const double clhs352 =             -clhs10*clhs351;
        const double clhs353 =             DN(3,1)*rho*tau1;
        const double clhs354 =             clhs10*clhs353;
        const double clhs355 =             DN(3,2)*rho*tau1;
        const double clhs356 =             clhs10*clhs355;
        const double clhs357 =             DN(3,0)*N[0];
        const double clhs358 =             clhs28*clhs347;
        const double clhs359 =             clhs33*clhs353;
        const double clhs360 =             clhs355*clhs38;
        const double clhs361 =             clhs347*clhs49;
        const double clhs362 =             clhs349*clhs47;
        const double clhs363 =             -clhs351*clhs49;
        const double clhs364 =             clhs353*clhs49;
        const double clhs365 =             clhs355*clhs49;
        const double clhs366 =             DN(3,0)*N[1];
        const double clhs367 =             clhs347*clhs60;
        const double clhs368 =             clhs353*clhs63;
        const double clhs369 =             clhs355*clhs66;
        const double clhs370 =             clhs347*clhs76;
        const double clhs371 =             clhs349*clhs74;
        const double clhs372 =             -clhs351*clhs76;
        const double clhs373 =             clhs353*clhs76;
        const double clhs374 =             clhs355*clhs76;
        const double clhs375 =             DN(3,0)*N[2];
        const double clhs376 =             clhs347*clhs86;
        const double clhs377 =             clhs353*clhs88;
        const double clhs378 =             clhs355*clhs90;
        const double clhs379 =             pow(DN(3,0), 2);
        const double clhs380 =             pow(N[3], 2);
        const double clhs381 =             clhs12*clhs380 + clhs349*clhs98 - clhs351*clhs99 + clhs380*clhs8;
        const double clhs382 =             clhs18*clhs99;
        const double clhs383 =             DN(3,0)*rho - DN(3,0)*tau2 + clhs2*clhs382;
        const double clhs384 =             DN(3,0)*N[3];
        const double clhs385 =             bdf0*clhs380*clhs6;
        const double clhs386 =             clhs109*clhs347;
        const double clhs387 =             clhs111*clhs353;
        const double clhs388 =             clhs113*clhs355;
        const double clhs389 =             bdf0*clhs380*clhs6*tau1;
        const double clhs390 =             DN(3,1)*N[0];
        const double clhs391 =             DN(3,1)*N[1];
        const double clhs392 =             DN(3,1)*N[2];
        const double clhs393 =             DN(3,1)*rho - DN(3,1)*tau2 + clhs3*clhs382;
        const double clhs394 =             pow(DN(3,1), 2);
        const double clhs395 =             DN(3,1)*N[3];
        const double clhs396 =             DN(3,2)*N[0];
        const double clhs397 =             DN(3,2)*N[1];
        const double clhs398 =             DN(3,2)*N[2];
        const double clhs399 =             DN(3,2)*rho - DN(3,2)*tau2 + clhs382*clhs4;
        const double clhs400 =             pow(DN(3,2), 2);
        const double clhs401 =             DN(3,2)*N[3];
        const double clhs402 =             2*N[3];
        const double clhs403 =             clhs382 + clhs402;

        lhs(0,0)=clhs0*rho - clhs0*tau2 + clhs1*clhs10*clhs2 + clhs15;
        lhs(0,1)=DN(0,1)*clhs20;
        lhs(0,2)=DN(0,2)*clhs20;
        lhs(0,3)=-clhs14*clhs28 + clhs2*clhs29 + clhs2*clhs34 + clhs2*clhs39 - clhs21*clhs23 - clhs21 + clhs24*clhs25 - clhs40*clhs41;
        lhs(0,4)=clhs2*clhs50 + clhs42 + clhs43 + clhs45 + clhs46 + clhs51 + clhs52;
        lhs(0,5)=clhs2*clhs54 + clhs53;
        lhs(0,6)=clhs2*clhs56 + clhs55;
        lhs(0,7)=-clhs14*clhs60 + clhs2*clhs61 + clhs2*clhs64 + clhs2*clhs67 - clhs23*clhs57 - clhs57 + clhs69;
        lhs(0,8)=clhs2*clhs77 + clhs70 + clhs71 + clhs72 + clhs73 + clhs78 + clhs79;
        lhs(0,9)=clhs2*clhs81 + clhs80;
        lhs(0,10)=clhs2*clhs83 + clhs82;
        lhs(0,11)=-clhs14*clhs86 + clhs2*clhs87 + clhs2*clhs89 + clhs2*clhs91 - clhs23*clhs84 - clhs84 + clhs93;
        lhs(0,12)=clhs100*clhs2 + clhs101 + clhs102 + clhs94 + clhs95 + clhs96 + clhs97;
        lhs(0,13)=clhs103 + clhs104*clhs2;
        lhs(0,14)=clhs105 + clhs106*clhs2;
        lhs(0,15)=-clhs107*clhs23 - clhs107 - clhs109*clhs14 + clhs110*clhs2 + clhs112*clhs2 + clhs114*clhs2 + clhs116;
        lhs(1,0)=DN(0,0)*clhs119;
        lhs(1,1)=clhs10*clhs3*clhs30 + clhs120*rho - clhs120*tau2 + clhs15;
        lhs(1,2)=DN(0,2)*clhs119;
        lhs(1,3)=-clhs121*clhs23 - clhs121 - clhs122*clhs40 - clhs14*clhs33 + clhs24*clhs31 + clhs29*clhs3 + clhs3*clhs34 + clhs3*clhs39;
        lhs(1,4)=clhs123 + clhs3*clhs50;
        lhs(1,5)=clhs124 + clhs3*clhs54 + clhs51 + clhs52;
        lhs(1,6)=clhs125 + clhs3*clhs56;
        lhs(1,7)=-clhs126*clhs23 - clhs126 + clhs127 - clhs14*clhs63 + clhs3*clhs61 + clhs3*clhs64 + clhs3*clhs67;
        lhs(1,8)=clhs128 + clhs3*clhs77;
        lhs(1,9)=clhs129 + clhs3*clhs81 + clhs78 + clhs79;
        lhs(1,10)=clhs130 + clhs3*clhs83;
        lhs(1,11)=-clhs131*clhs23 - clhs131 + clhs132 - clhs14*clhs88 + clhs3*clhs87 + clhs3*clhs89 + clhs3*clhs91;
        lhs(1,12)=clhs100*clhs3 + clhs133;
        lhs(1,13)=clhs101 + clhs102 + clhs104*clhs3 + clhs134;
        lhs(1,14)=clhs106*clhs3 + clhs135;
        lhs(1,15)=clhs110*clhs3 - clhs111*clhs14 + clhs112*clhs3 + clhs114*clhs3 - clhs136*clhs23 - clhs136 + clhs137;
        lhs(2,0)=DN(0,0)*clhs140;
        lhs(2,1)=DN(0,1)*clhs140;
        lhs(2,2)=clhs10*clhs35*clhs4 + clhs141*rho - clhs141*tau2 + clhs15;
        lhs(2,3)=-clhs14*clhs38 - clhs142*clhs23 - clhs142 - clhs143*clhs40 + clhs24*clhs36 + clhs29*clhs4 + clhs34*clhs4 + clhs39*clhs4;
        lhs(2,4)=clhs144 + clhs4*clhs50;
        lhs(2,5)=clhs145 + clhs4*clhs54;
        lhs(2,6)=clhs146 + clhs4*clhs56 + clhs51 + clhs52;
        lhs(2,7)=-clhs14*clhs66 - clhs147*clhs23 - clhs147 + clhs148 + clhs4*clhs61 + clhs4*clhs64 + clhs4*clhs67;
        lhs(2,8)=clhs149 + clhs4*clhs77;
        lhs(2,9)=clhs150 + clhs4*clhs81;
        lhs(2,10)=clhs151 + clhs4*clhs83 + clhs78 + clhs79;
        lhs(2,11)=-clhs14*clhs90 - clhs152*clhs23 - clhs152 + clhs153 + clhs4*clhs87 + clhs4*clhs89 + clhs4*clhs91;
        lhs(2,12)=clhs100*clhs4 + clhs154;
        lhs(2,13)=clhs104*clhs4 + clhs155;
        lhs(2,14)=clhs101 + clhs102 + clhs106*clhs4 + clhs156;
        lhs(2,15)=clhs110*clhs4 + clhs112*clhs4 - clhs113*clhs14 + clhs114*clhs4 - clhs157*clhs23 - clhs157 + clhs158;
        lhs(3,0)=DN(0,0)*clhs160;
        lhs(3,1)=DN(0,1)*clhs160;
        lhs(3,2)=DN(0,2)*clhs160;
        lhs(3,3)=clhs11*clhs161 + clhs29 + clhs34 + clhs39;
        lhs(3,4)=DN(1,0)*clhs159 + clhs50;
        lhs(3,5)=DN(1,1)*clhs159 + clhs54;
        lhs(3,6)=DN(1,2)*clhs159 + clhs56;
        lhs(3,7)=clhs163 + clhs61 + clhs64 + clhs67;
        lhs(3,8)=DN(2,0)*clhs159 + clhs77;
        lhs(3,9)=DN(2,1)*clhs159 + clhs81;
        lhs(3,10)=DN(2,2)*clhs159 + clhs83;
        lhs(3,11)=clhs164 + clhs87 + clhs89 + clhs91;
        lhs(3,12)=DN(3,0)*clhs159 + clhs100;
        lhs(3,13)=DN(3,1)*clhs159 + clhs104;
        lhs(3,14)=DN(3,2)*clhs159 + clhs106;
        lhs(3,15)=clhs110 + clhs112 + clhs114 + clhs165;
        lhs(4,0)=clhs167*clhs2 + clhs169 + clhs171 + clhs42 + clhs43 + clhs45 + clhs46;
        lhs(4,1)=clhs123 + clhs173*clhs2;
        lhs(4,2)=clhs144 + clhs175*clhs2;
        lhs(4,3)=-clhs170*clhs28 - clhs176*clhs23 - clhs176 + clhs177*clhs2 + clhs178*clhs2 + clhs179*clhs2 + clhs69;
        lhs(4,4)=clhs166*clhs2*clhs49 + clhs180*rho - clhs180*tau2 + clhs182;
        lhs(4,5)=DN(1,1)*clhs186;
        lhs(4,6)=DN(1,2)*clhs186;
        lhs(4,7)=-clhs170*clhs60 - clhs187*clhs23 - clhs187 + clhs188*clhs25 + clhs189*clhs2 + clhs190*clhs2 + clhs191*clhs2 - clhs192*clhs41;
        lhs(4,8)=clhs193 + clhs194 + clhs196 + clhs197 + clhs198*clhs2 + clhs199 + clhs200;
        lhs(4,9)=clhs2*clhs202 + clhs201;
        lhs(4,10)=clhs2*clhs204 + clhs203;
        lhs(4,11)=-clhs170*clhs86 + clhs2*clhs206 + clhs2*clhs207 + clhs2*clhs208 - clhs205*clhs23 - clhs205 + clhs210;
        lhs(4,12)=clhs2*clhs215 + clhs211 + clhs212 + clhs213 + clhs214 + clhs216 + clhs217;
        lhs(4,13)=clhs2*clhs219 + clhs218;
        lhs(4,14)=clhs2*clhs221 + clhs220;
        lhs(4,15)=-clhs109*clhs170 + clhs2*clhs223 + clhs2*clhs224 + clhs2*clhs225 - clhs222*clhs23 - clhs222 + clhs227;
        lhs(5,0)=clhs167*clhs3 + clhs53;
        lhs(5,1)=clhs124 + clhs169 + clhs171 + clhs173*clhs3;
        lhs(5,2)=clhs145 + clhs175*clhs3;
        lhs(5,3)=clhs127 - clhs170*clhs33 + clhs177*clhs3 + clhs178*clhs3 + clhs179*clhs3 - clhs228*clhs23 - clhs228;
        lhs(5,4)=DN(1,0)*clhs231;
        lhs(5,5)=clhs172*clhs3*clhs49 + clhs182 + clhs232*rho - clhs232*tau2;
        lhs(5,6)=DN(1,2)*clhs231;
        lhs(5,7)=-clhs122*clhs192 - clhs170*clhs63 + clhs188*clhs31 + clhs189*clhs3 + clhs190*clhs3 + clhs191*clhs3 - clhs23*clhs233 - clhs233;
        lhs(5,8)=clhs198*clhs3 + clhs234;
        lhs(5,9)=clhs199 + clhs200 + clhs202*clhs3 + clhs235;
        lhs(5,10)=clhs204*clhs3 + clhs236;
        lhs(5,11)=-clhs170*clhs88 + clhs206*clhs3 + clhs207*clhs3 + clhs208*clhs3 - clhs23*clhs237 - clhs237 + clhs238;
        lhs(5,12)=clhs215*clhs3 + clhs239;
        lhs(5,13)=clhs216 + clhs217 + clhs219*clhs3 + clhs240;
        lhs(5,14)=clhs221*clhs3 + clhs241;
        lhs(5,15)=-clhs111*clhs170 + clhs223*clhs3 + clhs224*clhs3 + clhs225*clhs3 - clhs23*clhs242 - clhs242 + clhs243;
        lhs(6,0)=clhs167*clhs4 + clhs55;
        lhs(6,1)=clhs125 + clhs173*clhs4;
        lhs(6,2)=clhs146 + clhs169 + clhs171 + clhs175*clhs4;
        lhs(6,3)=clhs148 - clhs170*clhs38 + clhs177*clhs4 + clhs178*clhs4 + clhs179*clhs4 - clhs23*clhs244 - clhs244;
        lhs(6,4)=DN(1,0)*clhs247;
        lhs(6,5)=DN(1,1)*clhs247;
        lhs(6,6)=clhs174*clhs4*clhs49 + clhs182 + clhs248*rho - clhs248*tau2;
        lhs(6,7)=-clhs143*clhs192 - clhs170*clhs66 + clhs188*clhs36 + clhs189*clhs4 + clhs190*clhs4 + clhs191*clhs4 - clhs23*clhs249 - clhs249;
        lhs(6,8)=clhs198*clhs4 + clhs250;
        lhs(6,9)=clhs202*clhs4 + clhs251;
        lhs(6,10)=clhs199 + clhs200 + clhs204*clhs4 + clhs252;
        lhs(6,11)=-clhs170*clhs90 + clhs206*clhs4 + clhs207*clhs4 + clhs208*clhs4 - clhs23*clhs253 - clhs253 + clhs254;
        lhs(6,12)=clhs215*clhs4 + clhs255;
        lhs(6,13)=clhs219*clhs4 + clhs256;
        lhs(6,14)=clhs216 + clhs217 + clhs221*clhs4 + clhs257;
        lhs(6,15)=-clhs113*clhs170 + clhs223*clhs4 + clhs224*clhs4 + clhs225*clhs4 - clhs23*clhs258 - clhs258 + clhs259;
        lhs(7,0)=DN(0,0)*clhs260 + clhs167;
        lhs(7,1)=DN(0,1)*clhs260 + clhs173;
        lhs(7,2)=DN(0,2)*clhs260 + clhs175;
        lhs(7,3)=clhs163 + clhs177 + clhs178 + clhs179;
        lhs(7,4)=DN(1,0)*clhs261;
        lhs(7,5)=DN(1,1)*clhs261;
        lhs(7,6)=DN(1,2)*clhs261;
        lhs(7,7)=clhs161*clhs181 + clhs189 + clhs190 + clhs191;
        lhs(7,8)=DN(2,0)*clhs260 + clhs198;
        lhs(7,9)=DN(2,1)*clhs260 + clhs202;
        lhs(7,10)=DN(2,2)*clhs260 + clhs204;
        lhs(7,11)=clhs206 + clhs207 + clhs208 + clhs263;
        lhs(7,12)=DN(3,0)*clhs260 + clhs215;
        lhs(7,13)=DN(3,1)*clhs260 + clhs219;
        lhs(7,14)=DN(3,2)*clhs260 + clhs221;
        lhs(7,15)=clhs223 + clhs224 + clhs225 + clhs264;
        lhs(8,0)=clhs2*clhs266 + clhs268 + clhs270 + clhs70 + clhs71 + clhs72 + clhs73;
        lhs(8,1)=clhs128 + clhs2*clhs272;
        lhs(8,2)=clhs149 + clhs2*clhs274;
        lhs(8,3)=clhs2*clhs276 + clhs2*clhs277 + clhs2*clhs278 - clhs23*clhs275 - clhs269*clhs28 - clhs275 + clhs93;
        lhs(8,4)=clhs193 + clhs194 + clhs196 + clhs197 + clhs2*clhs279 + clhs280 + clhs281;
        lhs(8,5)=clhs2*clhs282 + clhs234;
        lhs(8,6)=clhs2*clhs283 + clhs250;
        lhs(8,7)=clhs2*clhs285 + clhs2*clhs286 + clhs2*clhs287 + clhs210 - clhs23*clhs284 - clhs269*clhs60 - clhs284;
        lhs(8,8)=clhs2*clhs265*clhs76 + clhs288*rho - clhs288*tau2 + clhs290;
        lhs(8,9)=DN(2,1)*clhs294;
        lhs(8,10)=DN(2,2)*clhs294;
        lhs(8,11)=clhs2*clhs297 + clhs2*clhs298 + clhs2*clhs299 - clhs23*clhs295 + clhs25*clhs296 - clhs269*clhs86 - clhs295 - clhs300*clhs41;
        lhs(8,12)=clhs2*clhs306 + clhs301 + clhs302 + clhs304 + clhs305 + clhs307 + clhs308;
        lhs(8,13)=clhs2*clhs310 + clhs309;
        lhs(8,14)=clhs2*clhs312 + clhs311;
        lhs(8,15)=-clhs109*clhs269 + clhs2*clhs314 + clhs2*clhs315 + clhs2*clhs316 - clhs23*clhs313 - clhs313 + clhs319;
        lhs(9,0)=clhs266*clhs3 + clhs80;
        lhs(9,1)=clhs129 + clhs268 + clhs270 + clhs272*clhs3;
        lhs(9,2)=clhs150 + clhs274*clhs3;
        lhs(9,3)=clhs132 - clhs23*clhs320 - clhs269*clhs33 + clhs276*clhs3 + clhs277*clhs3 + clhs278*clhs3 - clhs320;
        lhs(9,4)=clhs201 + clhs279*clhs3;
        lhs(9,5)=clhs235 + clhs280 + clhs281 + clhs282*clhs3;
        lhs(9,6)=clhs251 + clhs283*clhs3;
        lhs(9,7)=-clhs23*clhs321 + clhs238 - clhs269*clhs63 + clhs285*clhs3 + clhs286*clhs3 + clhs287*clhs3 - clhs321;
        lhs(9,8)=DN(2,0)*clhs324;
        lhs(9,9)=clhs271*clhs3*clhs76 + clhs290 + clhs325*rho - clhs325*tau2;
        lhs(9,10)=DN(2,2)*clhs324;
        lhs(9,11)=-clhs122*clhs300 - clhs23*clhs326 - clhs269*clhs88 + clhs296*clhs31 + clhs297*clhs3 + clhs298*clhs3 + clhs299*clhs3 - clhs326;
        lhs(9,12)=clhs3*clhs306 + clhs327;
        lhs(9,13)=clhs3*clhs310 + clhs307 + clhs308 + clhs328;
        lhs(9,14)=clhs3*clhs312 + clhs329;
        lhs(9,15)=-clhs111*clhs269 - clhs23*clhs330 + clhs3*clhs314 + clhs3*clhs315 + clhs3*clhs316 - clhs330 + clhs331;
        lhs(10,0)=clhs266*clhs4 + clhs82;
        lhs(10,1)=clhs130 + clhs272*clhs4;
        lhs(10,2)=clhs151 + clhs268 + clhs270 + clhs274*clhs4;
        lhs(10,3)=clhs153 - clhs23*clhs332 - clhs269*clhs38 + clhs276*clhs4 + clhs277*clhs4 + clhs278*clhs4 - clhs332;
        lhs(10,4)=clhs203 + clhs279*clhs4;
        lhs(10,5)=clhs236 + clhs282*clhs4;
        lhs(10,6)=clhs252 + clhs280 + clhs281 + clhs283*clhs4;
        lhs(10,7)=-clhs23*clhs333 + clhs254 - clhs269*clhs66 + clhs285*clhs4 + clhs286*clhs4 + clhs287*clhs4 - clhs333;
        lhs(10,8)=DN(2,0)*clhs336;
        lhs(10,9)=DN(2,1)*clhs336;
        lhs(10,10)=clhs273*clhs4*clhs76 + clhs290 + clhs337*rho - clhs337*tau2;
        lhs(10,11)=-clhs143*clhs300 - clhs23*clhs338 - clhs269*clhs90 + clhs296*clhs36 + clhs297*clhs4 + clhs298*clhs4 + clhs299*clhs4 - clhs338;
        lhs(10,12)=clhs306*clhs4 + clhs339;
        lhs(10,13)=clhs310*clhs4 + clhs340;
        lhs(10,14)=clhs307 + clhs308 + clhs312*clhs4 + clhs341;
        lhs(10,15)=-clhs113*clhs269 - clhs23*clhs342 + clhs314*clhs4 + clhs315*clhs4 + clhs316*clhs4 - clhs342 + clhs343;
        lhs(11,0)=DN(0,0)*clhs344 + clhs266;
        lhs(11,1)=DN(0,1)*clhs344 + clhs272;
        lhs(11,2)=DN(0,2)*clhs344 + clhs274;
        lhs(11,3)=clhs164 + clhs276 + clhs277 + clhs278;
        lhs(11,4)=DN(1,0)*clhs344 + clhs279;
        lhs(11,5)=DN(1,1)*clhs344 + clhs282;
        lhs(11,6)=DN(1,2)*clhs344 + clhs283;
        lhs(11,7)=clhs263 + clhs285 + clhs286 + clhs287;
        lhs(11,8)=DN(2,0)*clhs345;
        lhs(11,9)=DN(2,1)*clhs345;
        lhs(11,10)=DN(2,2)*clhs345;
        lhs(11,11)=clhs161*clhs289 + clhs297 + clhs298 + clhs299;
        lhs(11,12)=DN(3,0)*clhs344 + clhs306;
        lhs(11,13)=DN(3,1)*clhs344 + clhs310;
        lhs(11,14)=DN(3,2)*clhs344 + clhs312;
        lhs(11,15)=clhs314 + clhs315 + clhs316 + clhs346;
        lhs(12,0)=clhs2*clhs348 + clhs350 + clhs352 + clhs94 + clhs95 + clhs96 + clhs97;
        lhs(12,1)=clhs133 + clhs2*clhs354;
        lhs(12,2)=clhs154 + clhs2*clhs356;
        lhs(12,3)=clhs116 + clhs2*clhs358 + clhs2*clhs359 + clhs2*clhs360 - clhs23*clhs357 - clhs28*clhs351 - clhs357;
        lhs(12,4)=clhs2*clhs361 + clhs211 + clhs212 + clhs213 + clhs214 + clhs362 + clhs363;
        lhs(12,5)=clhs2*clhs364 + clhs239;
        lhs(12,6)=clhs2*clhs365 + clhs255;
        lhs(12,7)=clhs2*clhs367 + clhs2*clhs368 + clhs2*clhs369 + clhs227 - clhs23*clhs366 - clhs351*clhs60 - clhs366;
        lhs(12,8)=clhs2*clhs370 + clhs301 + clhs302 + clhs304 + clhs305 + clhs371 + clhs372;
        lhs(12,9)=clhs2*clhs373 + clhs327;
        lhs(12,10)=clhs2*clhs374 + clhs339;
        lhs(12,11)=clhs2*clhs376 + clhs2*clhs377 + clhs2*clhs378 - clhs23*clhs375 + clhs319 - clhs351*clhs86 - clhs375;
        lhs(12,12)=clhs2*clhs347*clhs99 + clhs379*rho - clhs379*tau2 + clhs381;
        lhs(12,13)=DN(3,1)*clhs383;
        lhs(12,14)=DN(3,2)*clhs383;
        lhs(12,15)=-clhs109*clhs351 + clhs2*clhs386 + clhs2*clhs387 + clhs2*clhs388 - clhs23*clhs384 + clhs25*clhs385 - clhs384 - clhs389*clhs41;
        lhs(13,0)=clhs103 + clhs3*clhs348;
        lhs(13,1)=clhs134 + clhs3*clhs354 + clhs350 + clhs352;
        lhs(13,2)=clhs155 + clhs3*clhs356;
        lhs(13,3)=clhs137 - clhs23*clhs390 + clhs3*clhs358 + clhs3*clhs359 + clhs3*clhs360 - clhs33*clhs351 - clhs390;
        lhs(13,4)=clhs218 + clhs3*clhs361;
        lhs(13,5)=clhs240 + clhs3*clhs364 + clhs362 + clhs363;
        lhs(13,6)=clhs256 + clhs3*clhs365;
        lhs(13,7)=-clhs23*clhs391 + clhs243 + clhs3*clhs367 + clhs3*clhs368 + clhs3*clhs369 - clhs351*clhs63 - clhs391;
        lhs(13,8)=clhs3*clhs370 + clhs309;
        lhs(13,9)=clhs3*clhs373 + clhs328 + clhs371 + clhs372;
        lhs(13,10)=clhs3*clhs374 + clhs340;
        lhs(13,11)=-clhs23*clhs392 + clhs3*clhs376 + clhs3*clhs377 + clhs3*clhs378 + clhs331 - clhs351*clhs88 - clhs392;
        lhs(13,12)=DN(3,0)*clhs393;
        lhs(13,13)=clhs3*clhs353*clhs99 + clhs381 + clhs394*rho - clhs394*tau2;
        lhs(13,14)=DN(3,2)*clhs393;
        lhs(13,15)=-clhs111*clhs351 - clhs122*clhs389 - clhs23*clhs395 + clhs3*clhs386 + clhs3*clhs387 + clhs3*clhs388 + clhs31*clhs385 - clhs395;
        lhs(14,0)=clhs105 + clhs348*clhs4;
        lhs(14,1)=clhs135 + clhs354*clhs4;
        lhs(14,2)=clhs156 + clhs350 + clhs352 + clhs356*clhs4;
        lhs(14,3)=clhs158 - clhs23*clhs396 - clhs351*clhs38 + clhs358*clhs4 + clhs359*clhs4 + clhs360*clhs4 - clhs396;
        lhs(14,4)=clhs220 + clhs361*clhs4;
        lhs(14,5)=clhs241 + clhs364*clhs4;
        lhs(14,6)=clhs257 + clhs362 + clhs363 + clhs365*clhs4;
        lhs(14,7)=-clhs23*clhs397 + clhs259 - clhs351*clhs66 + clhs367*clhs4 + clhs368*clhs4 + clhs369*clhs4 - clhs397;
        lhs(14,8)=clhs311 + clhs370*clhs4;
        lhs(14,9)=clhs329 + clhs373*clhs4;
        lhs(14,10)=clhs341 + clhs371 + clhs372 + clhs374*clhs4;
        lhs(14,11)=-clhs23*clhs398 + clhs343 - clhs351*clhs90 + clhs376*clhs4 + clhs377*clhs4 + clhs378*clhs4 - clhs398;
        lhs(14,12)=DN(3,0)*clhs399;
        lhs(14,13)=DN(3,1)*clhs399;
        lhs(14,14)=clhs355*clhs4*clhs99 + clhs381 + clhs400*rho - clhs400*tau2;
        lhs(14,15)=-clhs113*clhs351 - clhs143*clhs389 - clhs23*clhs401 + clhs36*clhs385 + clhs386*clhs4 + clhs387*clhs4 + clhs388*clhs4 - clhs401;
        lhs(15,0)=DN(0,0)*clhs402 + clhs348;
        lhs(15,1)=DN(0,1)*clhs402 + clhs354;
        lhs(15,2)=DN(0,2)*clhs402 + clhs356;
        lhs(15,3)=clhs165 + clhs358 + clhs359 + clhs360;
        lhs(15,4)=DN(1,0)*clhs402 + clhs361;
        lhs(15,5)=DN(1,1)*clhs402 + clhs364;
        lhs(15,6)=DN(1,2)*clhs402 + clhs365;
        lhs(15,7)=clhs264 + clhs367 + clhs368 + clhs369;
        lhs(15,8)=DN(2,0)*clhs402 + clhs370;
        lhs(15,9)=DN(2,1)*clhs402 + clhs373;
        lhs(15,10)=DN(2,2)*clhs402 + clhs374;
        lhs(15,11)=clhs346 + clhs376 + clhs377 + clhs378;
        lhs(15,12)=DN(3,0)*clhs403;
        lhs(15,13)=DN(3,1)*clhs403;
        lhs(15,14)=DN(3,2)*clhs403;
        lhs(15,15)=clhs161*clhs380 + clhs386 + clhs387 + clhs388;

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
        const double clhs1 =             DN(0,0)*rho*tau1;
        const double clhs2 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
        const double clhs3 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
        const double clhs4 =             DN(0,0)*clhs2 + DN(0,1)*clhs3;
        const double clhs5 =             pow(c, -2);
        const double clhs6 =             N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]);
        const double clhs7 =             clhs5*clhs6;
        const double clhs8 =             N[0]*clhs7;
        const double clhs9 =             clhs8 + rho*(N[0]*bdf0 + clhs4);
        const double clhs10 =             pow(N[0], 2);
        const double clhs11 =             bdf0*rho;
        const double clhs12 =             N[0]*rho;
        const double clhs13 =             N[0]*clhs5*clhs6*tau1;
        const double clhs14 =             clhs10*clhs11 + clhs10*clhs7 + clhs12*clhs4 - clhs13*clhs9;
        const double clhs15 =             DN(0,0)*rho;
        const double clhs16 =             DN(0,0)*tau2;
        const double clhs17 =             rho*tau1;
        const double clhs18 =             clhs17*clhs9;
        const double clhs19 =             DN(0,0)*N[0];
        const double clhs20 =             1.0/rho;
        const double clhs21 =             bdf0*clhs20*clhs5*tau2;
        const double clhs22 =             bdf0*clhs10*clhs5;
        const double clhs23 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0);
        const double clhs24 =             N[0]*bdf0*clhs5;
        const double clhs25 =             clhs23*clhs24;
        const double clhs26 =             DN(0,0) + clhs25;
        const double clhs27 =             clhs1*clhs26;
        const double clhs28 =             DN(0,1)*rho*tau1;
        const double clhs29 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
        const double clhs30 =             clhs24*clhs29;
        const double clhs31 =             DN(0,1) + clhs30;
        const double clhs32 =             clhs28*clhs31;
        const double clhs33 =             bdf0*clhs10*clhs5*tau1;
        const double clhs34 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + clhs23*clhs7 - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + clhs2*(DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0)) + clhs3*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
        const double clhs35 =             DN(1,0)*clhs15;
        const double clhs36 =             -DN(1,0)*clhs16;
        const double clhs37 =             N[0]*bdf0*rho;
        const double clhs38 =             N[1]*clhs37;
        const double clhs39 =             N[1]*clhs8;
        const double clhs40 =             DN(1,0)*clhs2 + DN(1,1)*clhs3;
        const double clhs41 =             N[1]*clhs7;
        const double clhs42 =             clhs41 + rho*(N[1]*bdf0 + clhs40);
        const double clhs43 =             clhs1*clhs42;
        const double clhs44 =             clhs12*clhs40;
        const double clhs45 =             -clhs13*clhs42;
        const double clhs46 =             DN(1,1)*clhs15 - DN(1,1)*clhs16;
        const double clhs47 =             clhs28*clhs42;
        const double clhs48 =             DN(0,0)*N[1];
        const double clhs49 =             N[1]*bdf0*clhs5;
        const double clhs50 =             DN(1,0) + clhs23*clhs49;
        const double clhs51 =             clhs1*clhs50;
        const double clhs52 =             DN(1,1) + clhs29*clhs49;
        const double clhs53 =             clhs28*clhs52;
        const double clhs54 =             N[0]*N[1]*bdf0*clhs5*tau1;
        const double clhs55 =             N[1]*clhs25 - clhs34*clhs54;
        const double clhs56 =             DN(2,0)*clhs15;
        const double clhs57 =             -DN(2,0)*clhs16;
        const double clhs58 =             N[2]*clhs37;
        const double clhs59 =             N[2]*clhs8;
        const double clhs60 =             DN(2,0)*clhs2 + DN(2,1)*clhs3;
        const double clhs61 =             N[2]*clhs7 + rho*(N[2]*bdf0 + clhs60);
        const double clhs62 =             clhs1*clhs61;
        const double clhs63 =             clhs12*clhs60;
        const double clhs64 =             -clhs13*clhs61;
        const double clhs65 =             DN(2,1)*clhs15 - DN(2,1)*clhs16;
        const double clhs66 =             clhs28*clhs61;
        const double clhs67 =             DN(0,0)*N[2];
        const double clhs68 =             N[2]*bdf0*clhs5;
        const double clhs69 =             DN(2,0) + clhs23*clhs68;
        const double clhs70 =             clhs1*clhs69;
        const double clhs71 =             DN(2,1) + clhs29*clhs68;
        const double clhs72 =             clhs28*clhs71;
        const double clhs73 =             N[0]*N[2]*bdf0*clhs5*tau1;
        const double clhs74 =             N[2]*clhs25 - clhs34*clhs73;
        const double clhs75 =             DN(0,1)*rho;
        const double clhs76 =             DN(0,1)*tau2;
        const double clhs77 =             pow(DN(0,1), 2);
        const double clhs78 =             DN(0,1)*N[0];
        const double clhs79 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + clhs29*clhs7 - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + clhs2*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + clhs3*(DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1)));
        const double clhs80 =             DN(1,0)*clhs75 - DN(1,0)*clhs76;
        const double clhs81 =             DN(1,1)*clhs75 - DN(1,1)*clhs76 + clhs38 + clhs39;
        const double clhs82 =             DN(0,1)*N[1];
        const double clhs83 =             N[1]*clhs30 - clhs54*clhs79;
        const double clhs84 =             DN(2,0)*clhs75 - DN(2,0)*clhs76;
        const double clhs85 =             DN(2,1)*clhs75 - DN(2,1)*clhs76 + clhs58 + clhs59;
        const double clhs86 =             DN(0,1)*N[2];
        const double clhs87 =             N[2]*clhs30 - clhs73*clhs79;
        const double clhs88 =             2*N[0];
        const double clhs89 =             clhs18 + clhs88;
        const double clhs90 =             bdf0*clhs20*clhs5;
        const double clhs91 =             N[0]*bdf0*clhs20*clhs5;
        const double clhs92 =             N[1]*clhs91;
        const double clhs93 =             N[2]*clhs91;
        const double clhs94 =             DN(1,0)*rho*tau1;
        const double clhs95 =             clhs9*clhs94;
        const double clhs96 =             N[1]*rho;
        const double clhs97 =             clhs4*clhs96;
        const double clhs98 =             N[1]*clhs5*clhs6*tau1;
        const double clhs99 =             -clhs9*clhs98;
        const double clhs100 =             DN(1,1)*rho*tau1;
        const double clhs101 =             clhs100*clhs9;
        const double clhs102 =             DN(1,0)*N[0];
        const double clhs103 =             clhs26*clhs94;
        const double clhs104 =             clhs100*clhs31;
        const double clhs105 =             pow(DN(1,0), 2);
        const double clhs106 =             pow(N[1], 2);
        const double clhs107 =             clhs106*clhs11 + clhs106*clhs7 + clhs40*clhs96 - clhs42*clhs98;
        const double clhs108 =             DN(1,0)*rho;
        const double clhs109 =             DN(1,0)*tau2;
        const double clhs110 =             clhs17*clhs42;
        const double clhs111 =             DN(1,0)*N[1];
        const double clhs112 =             bdf0*clhs106*clhs5;
        const double clhs113 =             clhs50*clhs94;
        const double clhs114 =             clhs100*clhs52;
        const double clhs115 =             bdf0*clhs106*clhs5*tau1;
        const double clhs116 =             DN(2,0)*clhs108;
        const double clhs117 =             -DN(2,0)*clhs109;
        const double clhs118 =             N[1]*N[2]*bdf0;
        const double clhs119 =             clhs118*rho;
        const double clhs120 =             N[2]*clhs41;
        const double clhs121 =             clhs61*clhs94;
        const double clhs122 =             clhs60*clhs96;
        const double clhs123 =             -clhs61*clhs98;
        const double clhs124 =             DN(2,1)*clhs108 - DN(2,1)*clhs109;
        const double clhs125 =             clhs100*clhs61;
        const double clhs126 =             DN(1,0)*N[2];
        const double clhs127 =             clhs69*clhs94;
        const double clhs128 =             clhs100*clhs71;
        const double clhs129 =             N[1]*N[2]*bdf0*clhs5;
        const double clhs130 =             N[1]*N[2]*bdf0*clhs5*tau1;
        const double clhs131 =             clhs129*clhs23 - clhs130*clhs34;
        const double clhs132 =             DN(1,1)*N[0];
        const double clhs133 =             DN(1,1)*rho;
        const double clhs134 =             DN(1,1)*tau2;
        const double clhs135 =             pow(DN(1,1), 2);
        const double clhs136 =             DN(1,1)*N[1];
        const double clhs137 =             DN(2,0)*clhs133 - DN(2,0)*clhs134;
        const double clhs138 =             DN(2,1)*clhs133 - DN(2,1)*clhs134 + clhs119 + clhs120;
        const double clhs139 =             DN(1,1)*N[2];
        const double clhs140 =             clhs129*clhs29 - clhs130*clhs79;
        const double clhs141 =             2*N[1];
        const double clhs142 =             clhs110 + clhs141;
        const double clhs143 =             clhs118*clhs20*clhs5;
        const double clhs144 =             DN(2,0)*rho*tau1;
        const double clhs145 =             clhs144*clhs9;
        const double clhs146 =             N[2]*rho;
        const double clhs147 =             clhs146*clhs4;
        const double clhs148 =             N[2]*clhs5*clhs6*tau1;
        const double clhs149 =             -clhs148*clhs9;
        const double clhs150 =             DN(2,1)*rho*tau1;
        const double clhs151 =             clhs150*clhs9;
        const double clhs152 =             DN(2,0)*N[0];
        const double clhs153 =             clhs144*clhs26;
        const double clhs154 =             clhs150*clhs31;
        const double clhs155 =             clhs144*clhs42;
        const double clhs156 =             clhs146*clhs40;
        const double clhs157 =             -clhs148*clhs42;
        const double clhs158 =             clhs150*clhs42;
        const double clhs159 =             DN(2,0)*N[1];
        const double clhs160 =             clhs144*clhs50;
        const double clhs161 =             clhs150*clhs52;
        const double clhs162 =             pow(DN(2,0), 2);
        const double clhs163 =             pow(N[2], 2);
        const double clhs164 =             clhs11*clhs163 + clhs146*clhs60 - clhs148*clhs61 + clhs163*clhs7;
        const double clhs165 =             clhs17*clhs61;
        const double clhs166 =             DN(2,0)*N[2];
        const double clhs167 =             bdf0*clhs163*clhs5;
        const double clhs168 =             clhs144*clhs69;
        const double clhs169 =             clhs150*clhs71;
        const double clhs170 =             bdf0*clhs163*clhs5*tau1;
        const double clhs171 =             DN(2,1)*N[0];
        const double clhs172 =             DN(2,1)*N[1];
        const double clhs173 =             pow(DN(2,1), 2);
        const double clhs174 =             DN(2,1)*N[2];
        const double clhs175 =             2*N[2];
        const double clhs176 =             clhs165 + clhs175;

        lhs(0,0)=clhs0*rho - clhs0*tau2 + clhs1*clhs2*clhs9 + clhs14;
        lhs(0,1)=DN(0,1)*(clhs15 - clhs16 + clhs18*clhs2);
        lhs(0,2)=-clhs13*clhs26 - clhs19*clhs21 - clhs19 + clhs2*clhs27 + clhs2*clhs32 + clhs22*clhs23 - clhs33*clhs34;
        lhs(0,3)=clhs2*clhs43 + clhs35 + clhs36 + clhs38 + clhs39 + clhs44 + clhs45;
        lhs(0,4)=clhs2*clhs47 + clhs46;
        lhs(0,5)=-clhs13*clhs50 + clhs2*clhs51 + clhs2*clhs53 - clhs21*clhs48 - clhs48 + clhs55;
        lhs(0,6)=clhs2*clhs62 + clhs56 + clhs57 + clhs58 + clhs59 + clhs63 + clhs64;
        lhs(0,7)=clhs2*clhs66 + clhs65;
        lhs(0,8)=-clhs13*clhs69 + clhs2*clhs70 + clhs2*clhs72 - clhs21*clhs67 - clhs67 + clhs74;
        lhs(1,0)=DN(0,0)*(clhs18*clhs3 + clhs75 - clhs76);
        lhs(1,1)=clhs14 + clhs28*clhs3*clhs9 + clhs77*rho - clhs77*tau2;
        lhs(1,2)=-clhs13*clhs31 - clhs21*clhs78 + clhs22*clhs29 + clhs27*clhs3 + clhs3*clhs32 - clhs33*clhs79 - clhs78;
        lhs(1,3)=clhs3*clhs43 + clhs80;
        lhs(1,4)=clhs3*clhs47 + clhs44 + clhs45 + clhs81;
        lhs(1,5)=-clhs13*clhs52 - clhs21*clhs82 + clhs3*clhs51 + clhs3*clhs53 - clhs82 + clhs83;
        lhs(1,6)=clhs3*clhs62 + clhs84;
        lhs(1,7)=clhs3*clhs66 + clhs63 + clhs64 + clhs85;
        lhs(1,8)=-clhs13*clhs71 - clhs21*clhs86 + clhs3*clhs70 + clhs3*clhs72 - clhs86 + clhs87;
        lhs(2,0)=DN(0,0)*clhs89;
        lhs(2,1)=DN(0,1)*clhs89;
        lhs(2,2)=clhs10*clhs90 + clhs27 + clhs32;
        lhs(2,3)=DN(1,0)*clhs88 + clhs43;
        lhs(2,4)=DN(1,1)*clhs88 + clhs47;
        lhs(2,5)=clhs51 + clhs53 + clhs92;
        lhs(2,6)=DN(2,0)*clhs88 + clhs62;
        lhs(2,7)=DN(2,1)*clhs88 + clhs66;
        lhs(2,8)=clhs70 + clhs72 + clhs93;
        lhs(3,0)=clhs2*clhs95 + clhs35 + clhs36 + clhs38 + clhs39 + clhs97 + clhs99;
        lhs(3,1)=clhs101*clhs2 + clhs80;
        lhs(3,2)=-clhs102*clhs21 - clhs102 + clhs103*clhs2 + clhs104*clhs2 - clhs26*clhs98 + clhs55;
        lhs(3,3)=clhs105*rho - clhs105*tau2 + clhs107 + clhs2*clhs42*clhs94;
        lhs(3,4)=DN(1,1)*(clhs108 - clhs109 + clhs110*clhs2);
        lhs(3,5)=-clhs111*clhs21 - clhs111 + clhs112*clhs23 + clhs113*clhs2 + clhs114*clhs2 - clhs115*clhs34 - clhs50*clhs98;
        lhs(3,6)=clhs116 + clhs117 + clhs119 + clhs120 + clhs121*clhs2 + clhs122 + clhs123;
        lhs(3,7)=clhs124 + clhs125*clhs2;
        lhs(3,8)=-clhs126*clhs21 - clhs126 + clhs127*clhs2 + clhs128*clhs2 + clhs131 - clhs69*clhs98;
        lhs(4,0)=clhs3*clhs95 + clhs46;
        lhs(4,1)=clhs101*clhs3 + clhs81 + clhs97 + clhs99;
        lhs(4,2)=clhs103*clhs3 + clhs104*clhs3 - clhs132*clhs21 - clhs132 - clhs31*clhs98 + clhs83;
        lhs(4,3)=DN(1,0)*(clhs110*clhs3 + clhs133 - clhs134);
        lhs(4,4)=clhs100*clhs3*clhs42 + clhs107 + clhs135*rho - clhs135*tau2;
        lhs(4,5)=clhs112*clhs29 + clhs113*clhs3 + clhs114*clhs3 - clhs115*clhs79 - clhs136*clhs21 - clhs136 - clhs52*clhs98;
        lhs(4,6)=clhs121*clhs3 + clhs137;
        lhs(4,7)=clhs122 + clhs123 + clhs125*clhs3 + clhs138;
        lhs(4,8)=clhs127*clhs3 + clhs128*clhs3 - clhs139*clhs21 - clhs139 + clhs140 - clhs71*clhs98;
        lhs(5,0)=DN(0,0)*clhs141 + clhs95;
        lhs(5,1)=DN(0,1)*clhs141 + clhs101;
        lhs(5,2)=clhs103 + clhs104 + clhs92;
        lhs(5,3)=DN(1,0)*clhs142;
        lhs(5,4)=DN(1,1)*clhs142;
        lhs(5,5)=clhs106*clhs90 + clhs113 + clhs114;
        lhs(5,6)=DN(2,0)*clhs141 + clhs121;
        lhs(5,7)=DN(2,1)*clhs141 + clhs125;
        lhs(5,8)=clhs127 + clhs128 + clhs143;
        lhs(6,0)=clhs145*clhs2 + clhs147 + clhs149 + clhs56 + clhs57 + clhs58 + clhs59;
        lhs(6,1)=clhs151*clhs2 + clhs84;
        lhs(6,2)=-clhs148*clhs26 - clhs152*clhs21 - clhs152 + clhs153*clhs2 + clhs154*clhs2 + clhs74;
        lhs(6,3)=clhs116 + clhs117 + clhs119 + clhs120 + clhs155*clhs2 + clhs156 + clhs157;
        lhs(6,4)=clhs137 + clhs158*clhs2;
        lhs(6,5)=clhs131 - clhs148*clhs50 - clhs159*clhs21 - clhs159 + clhs160*clhs2 + clhs161*clhs2;
        lhs(6,6)=clhs144*clhs2*clhs61 + clhs162*rho - clhs162*tau2 + clhs164;
        lhs(6,7)=DN(2,1)*(DN(2,0)*rho - DN(2,0)*tau2 + clhs165*clhs2);
        lhs(6,8)=-clhs148*clhs69 - clhs166*clhs21 - clhs166 + clhs167*clhs23 + clhs168*clhs2 + clhs169*clhs2 - clhs170*clhs34;
        lhs(7,0)=clhs145*clhs3 + clhs65;
        lhs(7,1)=clhs147 + clhs149 + clhs151*clhs3 + clhs85;
        lhs(7,2)=-clhs148*clhs31 + clhs153*clhs3 + clhs154*clhs3 - clhs171*clhs21 - clhs171 + clhs87;
        lhs(7,3)=clhs124 + clhs155*clhs3;
        lhs(7,4)=clhs138 + clhs156 + clhs157 + clhs158*clhs3;
        lhs(7,5)=clhs140 - clhs148*clhs52 + clhs160*clhs3 + clhs161*clhs3 - clhs172*clhs21 - clhs172;
        lhs(7,6)=DN(2,0)*(DN(2,1)*rho - DN(2,1)*tau2 + clhs165*clhs3);
        lhs(7,7)=clhs150*clhs3*clhs61 + clhs164 + clhs173*rho - clhs173*tau2;
        lhs(7,8)=-clhs148*clhs71 + clhs167*clhs29 + clhs168*clhs3 + clhs169*clhs3 - clhs170*clhs79 - clhs174*clhs21 - clhs174;
        lhs(8,0)=DN(0,0)*clhs175 + clhs145;
        lhs(8,1)=DN(0,1)*clhs175 + clhs151;
        lhs(8,2)=clhs153 + clhs154 + clhs93;
        lhs(8,3)=DN(1,0)*clhs175 + clhs155;
        lhs(8,4)=DN(1,1)*clhs175 + clhs158;
        lhs(8,5)=clhs143 + clhs160 + clhs161;
        lhs(8,6)=DN(2,0)*clhs176;
        lhs(8,7)=DN(2,1)*clhs176;
        lhs(8,8)=clhs163*clhs90 + clhs168 + clhs169;

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

        // Auxiliary variables used in the calculation of the RHS
        const array_1d<double,dim> f_gauss = prod(trans(f), N);
        const array_1d<double,dim> grad_p = prod(trans(DN), p);
        //~ const double p_gauss = inner_prod(N,p);

        //~ array_1d<double,dim> accel_gauss = bdf0*v_gauss;
        //~ noalias(accel_gauss) += bdf1*prod(trans(vn), N);
        //~ noalias(accel_gauss) += bdf2*prod(trans(vnn), N);

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
        const double crhs26 =             crhs25*rho*tau1;
        const double crhs27 =             DN(0,0)*crhs26;
        const double crhs28 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
        const double crhs29 =             crhs17*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1));
        const double crhs30 =             N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1));
        const double crhs31 =             crhs19*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs20*(crhs10 + crhs4 + crhs6 + crhs8) + crhs21*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1));
        const double crhs32 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs28 + crhs29 + rho*(crhs30 + crhs31);
        const double crhs33 =             crhs32*rho*tau1;
        const double crhs34 =             DN(0,1)*crhs33;
        const double crhs35 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
        const double crhs36 =             crhs17*(N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2));
        const double crhs37 =             N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2));
        const double crhs38 =             crhs19*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs2*crhs21 + crhs20*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2));
        const double crhs39 =             DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs35 + crhs36 + rho*(crhs37 + crhs38);
        const double crhs40 =             crhs39*rho*tau1;
        const double crhs41 =             DN(0,2)*crhs40;
        const double crhs42 =             N[0]*crhs15*crhs16*tau1;
        const double crhs43 =             2*DN(0,0)*v(0,0) + 2*DN(0,1)*v(0,1) + 2*DN(0,2)*v(0,2) + 2*DN(1,0)*v(1,0) + 2*DN(1,1)*v(1,1) + 2*DN(1,2)*v(1,2) + 2*DN(2,0)*v(2,0) + 2*DN(2,1)*v(2,1) + 2*DN(2,2)*v(2,2) + 2*DN(3,0)*v(3,0) + 2*DN(3,1)*v(3,1) + 2*DN(3,2)*v(3,2);
        const double crhs44 =             N[1]*rho;
        const double crhs45 =             DN(1,0)*crhs26;
        const double crhs46 =             DN(1,1)*crhs33;
        const double crhs47 =             DN(1,2)*crhs40;
        const double crhs48 =             N[1]*crhs15*crhs16*tau1;
        const double crhs49 =             N[2]*rho;
        const double crhs50 =             DN(2,0)*crhs26;
        const double crhs51 =             DN(2,1)*crhs33;
        const double crhs52 =             DN(2,2)*crhs40;
        const double crhs53 =             N[2]*crhs15*crhs16*tau1;
        const double crhs54 =             N[3]*rho;
        const double crhs55 =             DN(3,0)*crhs26;
        const double crhs56 =             DN(3,1)*crhs33;
        const double crhs57 =             DN(3,2)*crhs40;
        const double crhs58 =             N[3]*crhs15*crhs16*tau1;

        rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs12 + DN(0,0)*crhs24 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - N[0]*crhs18 - crhs13*crhs14 - crhs13*crhs22 - crhs19*crhs27 - crhs19*crhs34 - crhs19*crhs41 + crhs25*crhs42;
        rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs12 + DN(0,1)*crhs24 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs28 - N[0]*crhs29 - crhs13*crhs30 - crhs13*crhs31 - crhs20*crhs27 - crhs20*crhs34 - crhs20*crhs41 + crhs32*crhs42;
        rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs12 + DN(0,2)*crhs24 - DN(0,2)*stress[2] + N[0]*crhs35 - N[0]*crhs36 - crhs13*crhs37 - crhs13*crhs38 - crhs21*crhs27 - crhs21*crhs34 - crhs21*crhs41 + crhs39*crhs42;
        rhs[3]=-N[0]*crhs23 - N[0]*crhs43 - crhs27 - crhs34 - crhs41;
        rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs12 + DN(1,0)*crhs24 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - N[1]*crhs18 - crhs14*crhs44 - crhs19*crhs45 - crhs19*crhs46 - crhs19*crhs47 - crhs22*crhs44 + crhs25*crhs48;
        rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs12 + DN(1,1)*crhs24 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs28 - N[1]*crhs29 - crhs20*crhs45 - crhs20*crhs46 - crhs20*crhs47 - crhs30*crhs44 - crhs31*crhs44 + crhs32*crhs48;
        rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs12 + DN(1,2)*crhs24 - DN(1,2)*stress[2] + N[1]*crhs35 - N[1]*crhs36 - crhs21*crhs45 - crhs21*crhs46 - crhs21*crhs47 - crhs37*crhs44 - crhs38*crhs44 + crhs39*crhs48;
        rhs[7]=-N[1]*crhs23 - N[1]*crhs43 - crhs45 - crhs46 - crhs47;
        rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs12 + DN(2,0)*crhs24 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - N[2]*crhs18 - crhs14*crhs49 - crhs19*crhs50 - crhs19*crhs51 - crhs19*crhs52 - crhs22*crhs49 + crhs25*crhs53;
        rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs12 + DN(2,1)*crhs24 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs28 - N[2]*crhs29 - crhs20*crhs50 - crhs20*crhs51 - crhs20*crhs52 - crhs30*crhs49 - crhs31*crhs49 + crhs32*crhs53;
        rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs12 + DN(2,2)*crhs24 - DN(2,2)*stress[2] + N[2]*crhs35 - N[2]*crhs36 - crhs21*crhs50 - crhs21*crhs51 - crhs21*crhs52 - crhs37*crhs49 - crhs38*crhs49 + crhs39*crhs53;
        rhs[11]=-N[2]*crhs23 - N[2]*crhs43 - crhs50 - crhs51 - crhs52;
        rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs12 + DN(3,0)*crhs24 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - N[3]*crhs18 - crhs14*crhs54 - crhs19*crhs55 - crhs19*crhs56 - crhs19*crhs57 - crhs22*crhs54 + crhs25*crhs58;
        rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs12 + DN(3,1)*crhs24 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs28 - N[3]*crhs29 - crhs20*crhs55 - crhs20*crhs56 - crhs20*crhs57 - crhs30*crhs54 - crhs31*crhs54 + crhs32*crhs58;
        rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs12 + DN(3,2)*crhs24 - DN(3,2)*stress[2] + N[3]*crhs35 - N[3]*crhs36 - crhs21*crhs55 - crhs21*crhs56 - crhs21*crhs57 - crhs37*crhs54 - crhs38*crhs54 + crhs39*crhs58;
        rhs[15]=-N[3]*crhs23 - N[3]*crhs43 - crhs55 - crhs56 - crhs57;

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

        // Auxiliary variables used in the calculation of the RHS
        const array_1d<double,dim> f_gauss = prod(trans(f), N);
        const array_1d<double,dim> grad_p = prod(trans(DN), p);
        //~ const double p_gauss = inner_prod(N,p);

        //~ array_1d<double,dim> accel_gauss = bdf0*v_gauss;
        //~ noalias(accel_gauss) += bdf1*prod(trans(vn), N);
        //~ noalias(accel_gauss) += bdf2*prod(trans(vnn), N);

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
        const double crhs18 =             crhs17*rho*tau1;
        const double crhs19 =             DN(0,0)*crhs18;
        const double crhs20 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
        const double crhs21 =             crhs13*(N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1));
        const double crhs22 =             N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1));
        const double crhs23 =             crhs3*crhs9 + crhs8*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1));
        const double crhs24 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs20 + crhs21 + rho*(crhs22 + crhs23);
        const double crhs25 =             crhs24*rho*tau1;
        const double crhs26 =             DN(0,1)*crhs25;
        const double crhs27 =             N[0]*crhs11*crhs12*tau1;
        const double crhs28 =             2*DN(0,0)*v(0,0) + 2*DN(0,1)*v(0,1) + 2*DN(1,0)*v(1,0) + 2*DN(1,1)*v(1,1) + 2*DN(2,0)*v(2,0) + 2*DN(2,1)*v(2,1);
        const double crhs29 =             N[1]*rho;
        const double crhs30 =             DN(1,0)*crhs18;
        const double crhs31 =             DN(1,1)*crhs25;
        const double crhs32 =             N[1]*crhs11*crhs12*tau1;
        const double crhs33 =             N[2]*rho;
        const double crhs34 =             DN(2,0)*crhs18;
        const double crhs35 =             DN(2,1)*crhs25;
        const double crhs36 =             N[2]*crhs11*crhs12*tau1;

        rhs[0]=DN(0,0)*crhs0 + DN(0,0)*crhs16 - DN(0,0)*crhs5 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - N[0]*crhs14 - crhs10*crhs6 + crhs17*crhs27 - crhs19*crhs8 - crhs26*crhs8 - crhs6*crhs7;
        rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 + DN(0,1)*crhs16 - DN(0,1)*crhs5 - DN(0,1)*stress[1] + N[0]*crhs20 - N[0]*crhs21 - crhs19*crhs9 - crhs22*crhs6 - crhs23*crhs6 + crhs24*crhs27 - crhs26*crhs9;
        rhs[2]=-N[0]*crhs15 - N[0]*crhs28 - crhs19 - crhs26;
        rhs[3]=DN(1,0)*crhs0 + DN(1,0)*crhs16 - DN(1,0)*crhs5 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - N[1]*crhs14 - crhs10*crhs29 + crhs17*crhs32 - crhs29*crhs7 - crhs30*crhs8 - crhs31*crhs8;
        rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 + DN(1,1)*crhs16 - DN(1,1)*crhs5 - DN(1,1)*stress[1] + N[1]*crhs20 - N[1]*crhs21 - crhs22*crhs29 - crhs23*crhs29 + crhs24*crhs32 - crhs30*crhs9 - crhs31*crhs9;
        rhs[5]=-N[1]*crhs15 - N[1]*crhs28 - crhs30 - crhs31;
        rhs[6]=DN(2,0)*crhs0 + DN(2,0)*crhs16 - DN(2,0)*crhs5 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - N[2]*crhs14 - crhs10*crhs33 + crhs17*crhs36 - crhs33*crhs7 - crhs34*crhs8 - crhs35*crhs8;
        rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 + DN(2,1)*crhs16 - DN(2,1)*crhs5 - DN(2,1)*stress[1] + N[2]*crhs20 - N[2]*crhs21 - crhs22*crhs33 - crhs23*crhs33 + crhs24*crhs36 - crhs34*crhs9 - crhs35*crhs9;
        rhs[8]=-N[2]*crhs15 - N[2]*crhs28 - crhs34 - crhs35;

    }
}
