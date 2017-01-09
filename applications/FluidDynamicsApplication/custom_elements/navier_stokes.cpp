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

        const double clhs0 =             pow(DN(0,0), 2)*rho;
        const double clhs1 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
        const double clhs2 =             DN(0,0)*clhs1*rho*tau1;
        const double clhs3 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
        const double clhs4 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
        const double clhs5 =             DN(0,0)*clhs1 + DN(0,1)*clhs3 + DN(0,2)*clhs4;
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
        const double clhs16 =             DN(0,1)*rho;
        const double clhs17 =             DN(0,0)*tau2;
        const double clhs18 =             clhs10*tau1;
        const double clhs19 =             DN(0,0) + clhs1*clhs18 - clhs17;
        const double clhs20 =             DN(0,2)*rho;
        const double clhs21 =             DN(0,0)*N[0];
        const double clhs22 =             bdf0*clhs6*tau2;
        const double clhs23 =             bdf0*clhs6;
        const double clhs24 =             clhs11*clhs23;
        const double clhs25 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0);
        const double clhs26 =             N[0]*bdf0*clhs6;
        const double clhs27 =             DN(0,0) + clhs25*clhs26;
        const double clhs28 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
        const double clhs29 =             DN(0,1) + clhs26*clhs28;
        const double clhs30 =             N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2);
        const double clhs31 =             DN(0,2) + clhs26*clhs30;
        const double clhs32 =             rho*tau1*(DN(0,0)*clhs27 + DN(0,1)*clhs29 + DN(0,2)*clhs31);
        const double clhs33 =             bdf0*clhs11*clhs6*tau1;
        const double clhs34 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] + clhs25*clhs8 - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0) + N[3]*f(3,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0)) + clhs1*(DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0) + DN(3,0)*v(3,0)) + clhs3*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + clhs4*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)));
        const double clhs35 =             DN(0,0)*rho;
        const double clhs36 =             DN(1,0)*clhs35;
        const double clhs37 =             N[0]*bdf0*rho;
        const double clhs38 =             N[1]*clhs37;
        const double clhs39 =             DN(1,0)*rho;
        const double clhs40 =             -clhs17*clhs39;
        const double clhs41 =             N[1]*clhs9;
        const double clhs42 =             DN(1,0)*clhs1 + DN(1,1)*clhs3 + DN(1,2)*clhs4;
        const double clhs43 =             N[1]*clhs8;
        const double clhs44 =             clhs43 + rho*(N[1]*bdf0 + clhs42);
        const double clhs45 =             clhs13*clhs42;
        const double clhs46 =             -clhs14*clhs44;
        const double clhs47 =             DN(0,0)*DN(1,1);
        const double clhs48 =             -clhs47*tau2 + clhs47;
        const double clhs49 =             clhs44*tau1;
        const double clhs50 =             DN(0,1)*clhs49;
        const double clhs51 =             DN(0,0)*DN(1,2);
        const double clhs52 =             -clhs51*tau2 + clhs51;
        const double clhs53 =             DN(0,2)*clhs49;
        const double clhs54 =             DN(0,0)*N[1];
        const double clhs55 =             N[1]*bdf0*clhs6;
        const double clhs56 =             DN(1,0) + clhs25*clhs55;
        const double clhs57 =             DN(1,1) + clhs28*clhs55;
        const double clhs58 =             DN(1,2) + clhs30*clhs55;
        const double clhs59 =             rho*tau1*(DN(0,0)*clhs56 + DN(0,1)*clhs57 + DN(0,2)*clhs58);
        const double clhs60 =             N[1]*clhs26;
        const double clhs61 =             N[0]*N[1]*bdf0*clhs6*tau1;
        const double clhs62 =             clhs25*clhs60 - clhs34*clhs61;
        const double clhs63 =             DN(2,0)*clhs35;
        const double clhs64 =             N[2]*clhs37;
        const double clhs65 =             DN(2,0)*rho;
        const double clhs66 =             -clhs17*clhs65;
        const double clhs67 =             N[2]*clhs9;
        const double clhs68 =             DN(2,0)*clhs1 + DN(2,1)*clhs3 + DN(2,2)*clhs4;
        const double clhs69 =             N[2]*clhs8;
        const double clhs70 =             clhs69 + rho*(N[2]*bdf0 + clhs68);
        const double clhs71 =             clhs13*clhs68;
        const double clhs72 =             -clhs14*clhs70;
        const double clhs73 =             DN(0,0)*DN(2,1);
        const double clhs74 =             -clhs73*tau2 + clhs73;
        const double clhs75 =             clhs70*tau1;
        const double clhs76 =             DN(0,1)*clhs75;
        const double clhs77 =             DN(0,0)*DN(2,2);
        const double clhs78 =             -clhs77*tau2 + clhs77;
        const double clhs79 =             DN(0,2)*clhs75;
        const double clhs80 =             DN(0,0)*N[2];
        const double clhs81 =             N[2]*bdf0*clhs6;
        const double clhs82 =             DN(2,0) + clhs25*clhs81;
        const double clhs83 =             DN(2,1) + clhs28*clhs81;
        const double clhs84 =             DN(2,2) + clhs30*clhs81;
        const double clhs85 =             rho*tau1*(DN(0,0)*clhs82 + DN(0,1)*clhs83 + DN(0,2)*clhs84);
        const double clhs86 =             N[2]*clhs26;
        const double clhs87 =             N[0]*N[2]*bdf0*clhs6*tau1;
        const double clhs88 =             clhs25*clhs86 - clhs34*clhs87;
        const double clhs89 =             DN(3,0)*clhs35;
        const double clhs90 =             N[3]*clhs37;
        const double clhs91 =             DN(3,0)*rho;
        const double clhs92 =             -clhs17*clhs91;
        const double clhs93 =             N[3]*clhs9;
        const double clhs94 =             DN(3,0)*clhs1 + DN(3,1)*clhs3 + DN(3,2)*clhs4;
        const double clhs95 =             N[3]*clhs8 + rho*(N[3]*bdf0 + clhs94);
        const double clhs96 =             clhs13*clhs94;
        const double clhs97 =             -clhs14*clhs95;
        const double clhs98 =             DN(0,0)*DN(3,1);
        const double clhs99 =             -clhs98*tau2 + clhs98;
        const double clhs100 =             clhs95*tau1;
        const double clhs101 =             DN(0,1)*clhs100;
        const double clhs102 =             DN(0,0)*DN(3,2);
        const double clhs103 =             -clhs102*tau2 + clhs102;
        const double clhs104 =             DN(0,2)*clhs100;
        const double clhs105 =             DN(0,0)*N[3];
        const double clhs106 =             N[3]*bdf0*clhs6;
        const double clhs107 =             DN(3,0) + clhs106*clhs25;
        const double clhs108 =             DN(3,1) + clhs106*clhs28;
        const double clhs109 =             DN(3,2) + clhs106*clhs30;
        const double clhs110 =             rho*tau1*(DN(0,0)*clhs107 + DN(0,1)*clhs108 + DN(0,2)*clhs109);
        const double clhs111 =             N[3]*clhs26;
        const double clhs112 =             N[0]*N[3]*bdf0*clhs6*tau1;
        const double clhs113 =             clhs111*clhs25 - clhs112*clhs34;
        const double clhs114 =             DN(0,1)*tau2;
        const double clhs115 =             DN(0,1) - clhs114 + clhs18*clhs3;
        const double clhs116 =             pow(DN(0,1), 2)*rho;
        const double clhs117 =             DN(0,1)*clhs3*rho*tau1;
        const double clhs118 =             DN(0,1)*N[0];
        const double clhs119 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + clhs28*clhs8 - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1)) + clhs1*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + clhs3*(DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1) + DN(3,1)*v(3,1)) + clhs4*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)));
        const double clhs120 =             DN(0,1)*DN(1,0);
        const double clhs121 =             -clhs120*tau2 + clhs120;
        const double clhs122 =             DN(0,0)*clhs49;
        const double clhs123 =             DN(1,1)*rho;
        const double clhs124 =             DN(1,1)*clhs16 - clhs114*clhs123 + clhs38 + clhs41;
        const double clhs125 =             DN(0,1)*DN(1,2);
        const double clhs126 =             -clhs125*tau2 + clhs125;
        const double clhs127 =             DN(0,1)*N[1];
        const double clhs128 =             -clhs119*clhs61 + clhs28*clhs60;
        const double clhs129 =             DN(0,1)*DN(2,0);
        const double clhs130 =             -clhs129*tau2 + clhs129;
        const double clhs131 =             DN(0,0)*clhs75;
        const double clhs132 =             DN(2,1)*rho;
        const double clhs133 =             DN(2,1)*clhs16 - clhs114*clhs132 + clhs64 + clhs67;
        const double clhs134 =             DN(0,1)*DN(2,2);
        const double clhs135 =             -clhs134*tau2 + clhs134;
        const double clhs136 =             DN(0,1)*N[2];
        const double clhs137 =             -clhs119*clhs87 + clhs28*clhs86;
        const double clhs138 =             DN(0,1)*DN(3,0);
        const double clhs139 =             -clhs138*tau2 + clhs138;
        const double clhs140 =             DN(0,0)*clhs100;
        const double clhs141 =             DN(3,1)*rho;
        const double clhs142 =             DN(3,1)*clhs16 - clhs114*clhs141 + clhs90 + clhs93;
        const double clhs143 =             DN(0,1)*DN(3,2);
        const double clhs144 =             -clhs143*tau2 + clhs143;
        const double clhs145 =             DN(0,1)*N[3];
        const double clhs146 =             clhs111*clhs28 - clhs112*clhs119;
        const double clhs147 =             DN(0,2)*tau2;
        const double clhs148 =             DN(0,2) - clhs147 + clhs18*clhs4;
        const double clhs149 =             pow(DN(0,2), 2)*rho;
        const double clhs150 =             DN(0,2)*clhs4*rho*tau1;
        const double clhs151 =             DN(0,2)*N[0];
        const double clhs152 =             DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + clhs30*clhs8 - rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2)) + rho*(N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2)) + clhs1*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + clhs3*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + clhs4*(DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2)));
        const double clhs153 =             DN(0,2)*DN(1,0);
        const double clhs154 =             -clhs153*tau2 + clhs153;
        const double clhs155 =             DN(0,2)*DN(1,1);
        const double clhs156 =             -clhs155*tau2 + clhs155;
        const double clhs157 =             DN(1,2)*rho;
        const double clhs158 =             DN(1,2)*clhs20 - clhs147*clhs157 + clhs38 + clhs41;
        const double clhs159 =             DN(0,2)*N[1];
        const double clhs160 =             -clhs152*clhs61 + clhs30*clhs60;
        const double clhs161 =             DN(0,2)*DN(2,0);
        const double clhs162 =             -clhs161*tau2 + clhs161;
        const double clhs163 =             DN(0,2)*DN(2,1);
        const double clhs164 =             -clhs163*tau2 + clhs163;
        const double clhs165 =             DN(2,2)*rho;
        const double clhs166 =             DN(2,2)*clhs20 - clhs147*clhs165 + clhs64 + clhs67;
        const double clhs167 =             DN(0,2)*N[2];
        const double clhs168 =             -clhs152*clhs87 + clhs30*clhs86;
        const double clhs169 =             DN(0,2)*DN(3,0);
        const double clhs170 =             -clhs169*tau2 + clhs169;
        const double clhs171 =             DN(0,2)*DN(3,1);
        const double clhs172 =             -clhs171*tau2 + clhs171;
        const double clhs173 =             DN(3,2)*rho;
        const double clhs174 =             DN(3,2)*clhs20 - clhs147*clhs173 + clhs90 + clhs93;
        const double clhs175 =             DN(0,2)*N[3];
        const double clhs176 =             clhs111*clhs30 - clhs112*clhs152;
        const double clhs177 =             2*N[0];
        const double clhs178 =             clhs177 + clhs18;
        const double clhs179 =             DN(0,0)*rho*tau1;
        const double clhs180 =             DN(0,1)*rho*tau1;
        const double clhs181 =             DN(0,2)*rho*tau1;
        const double clhs182 =             DN(1,0)*clhs1*rho*tau1;
        const double clhs183 =             N[1]*rho;
        const double clhs184 =             clhs183*clhs5;
        const double clhs185 =             N[1]*clhs6*clhs7*tau1;
        const double clhs186 =             -clhs10*clhs185;
        const double clhs187 =             DN(1,1)*clhs18;
        const double clhs188 =             DN(1,2)*clhs18;
        const double clhs189 =             DN(1,0)*N[0];
        const double clhs190 =             rho*tau1*(DN(1,0)*clhs27 + DN(1,1)*clhs29 + DN(1,2)*clhs31);
        const double clhs191 =             pow(DN(1,0), 2)*rho;
        const double clhs192 =             pow(N[1], 2);
        const double clhs193 =             clhs12*clhs192 + clhs183*clhs42 - clhs185*clhs44 + clhs192*clhs8;
        const double clhs194 =             DN(1,0)*tau2;
        const double clhs195 =             DN(1,0) + clhs1*clhs49 - clhs194;
        const double clhs196 =             DN(1,0)*N[1];
        const double clhs197 =             clhs192*clhs23;
        const double clhs198 =             rho*tau1*(DN(1,0)*clhs56 + DN(1,1)*clhs57 + DN(1,2)*clhs58);
        const double clhs199 =             bdf0*clhs192*clhs6*tau1;
        const double clhs200 =             DN(2,0)*clhs39;
        const double clhs201 =             N[1]*bdf0*rho;
        const double clhs202 =             N[2]*clhs201;
        const double clhs203 =             -clhs194*clhs65;
        const double clhs204 =             N[2]*clhs43;
        const double clhs205 =             clhs183*clhs68;
        const double clhs206 =             -clhs185*clhs70;
        const double clhs207 =             DN(1,0)*DN(2,1);
        const double clhs208 =             -clhs207*tau2 + clhs207;
        const double clhs209 =             DN(1,1)*clhs75;
        const double clhs210 =             DN(1,0)*DN(2,2);
        const double clhs211 =             -clhs210*tau2 + clhs210;
        const double clhs212 =             DN(1,2)*clhs75;
        const double clhs213 =             DN(1,0)*N[2];
        const double clhs214 =             rho*tau1*(DN(1,0)*clhs82 + DN(1,1)*clhs83 + DN(1,2)*clhs84);
        const double clhs215 =             N[2]*clhs55;
        const double clhs216 =             N[1]*N[2]*bdf0*clhs6*tau1;
        const double clhs217 =             clhs215*clhs25 - clhs216*clhs34;
        const double clhs218 =             DN(3,0)*clhs39;
        const double clhs219 =             N[3]*clhs201;
        const double clhs220 =             -clhs194*clhs91;
        const double clhs221 =             N[3]*clhs43;
        const double clhs222 =             clhs183*clhs94;
        const double clhs223 =             -clhs185*clhs95;
        const double clhs224 =             DN(1,0)*DN(3,1);
        const double clhs225 =             -clhs224*tau2 + clhs224;
        const double clhs226 =             DN(1,1)*clhs100;
        const double clhs227 =             DN(1,0)*DN(3,2);
        const double clhs228 =             -clhs227*tau2 + clhs227;
        const double clhs229 =             DN(1,2)*clhs100;
        const double clhs230 =             DN(1,0)*N[3];
        const double clhs231 =             rho*tau1*(DN(1,0)*clhs107 + DN(1,1)*clhs108 + DN(1,2)*clhs109);
        const double clhs232 =             N[3]*clhs55;
        const double clhs233 =             N[1]*N[3]*bdf0*clhs6*tau1;
        const double clhs234 =             clhs232*clhs25 - clhs233*clhs34;
        const double clhs235 =             DN(1,0)*clhs18;
        const double clhs236 =             DN(1,1)*clhs3*rho*tau1;
        const double clhs237 =             DN(1,1)*N[0];
        const double clhs238 =             DN(1,1)*tau2;
        const double clhs239 =             DN(1,1) - clhs238 + clhs3*clhs49;
        const double clhs240 =             pow(DN(1,1), 2)*rho;
        const double clhs241 =             DN(1,1)*N[1];
        const double clhs242 =             DN(1,1)*DN(2,0);
        const double clhs243 =             -clhs242*tau2 + clhs242;
        const double clhs244 =             DN(1,0)*clhs75;
        const double clhs245 =             DN(2,1)*clhs123 - clhs132*clhs238 + clhs202 + clhs204;
        const double clhs246 =             DN(1,1)*DN(2,2);
        const double clhs247 =             -clhs246*tau2 + clhs246;
        const double clhs248 =             DN(1,1)*N[2];
        const double clhs249 =             -clhs119*clhs216 + clhs215*clhs28;
        const double clhs250 =             DN(1,1)*DN(3,0);
        const double clhs251 =             -clhs250*tau2 + clhs250;
        const double clhs252 =             DN(1,0)*clhs100;
        const double clhs253 =             DN(3,1)*clhs123 - clhs141*clhs238 + clhs219 + clhs221;
        const double clhs254 =             DN(1,1)*DN(3,2);
        const double clhs255 =             -clhs254*tau2 + clhs254;
        const double clhs256 =             DN(1,1)*N[3];
        const double clhs257 =             -clhs119*clhs233 + clhs232*clhs28;
        const double clhs258 =             DN(1,2)*clhs4*rho*tau1;
        const double clhs259 =             DN(1,2)*N[0];
        const double clhs260 =             DN(1,2)*tau2;
        const double clhs261 =             DN(1,2) - clhs260 + clhs4*clhs49;
        const double clhs262 =             pow(DN(1,2), 2)*rho;
        const double clhs263 =             DN(1,2)*N[1];
        const double clhs264 =             DN(1,2)*DN(2,0);
        const double clhs265 =             -clhs264*tau2 + clhs264;
        const double clhs266 =             DN(1,2)*DN(2,1);
        const double clhs267 =             -clhs266*tau2 + clhs266;
        const double clhs268 =             DN(2,2)*clhs157 - clhs165*clhs260 + clhs202 + clhs204;
        const double clhs269 =             DN(1,2)*N[2];
        const double clhs270 =             -clhs152*clhs216 + clhs215*clhs30;
        const double clhs271 =             DN(1,2)*DN(3,0);
        const double clhs272 =             -clhs271*tau2 + clhs271;
        const double clhs273 =             DN(1,2)*DN(3,1);
        const double clhs274 =             -clhs273*tau2 + clhs273;
        const double clhs275 =             DN(3,2)*clhs157 - clhs173*clhs260 + clhs219 + clhs221;
        const double clhs276 =             DN(1,2)*N[3];
        const double clhs277 =             -clhs152*clhs233 + clhs232*clhs30;
        const double clhs278 =             2*N[1];
        const double clhs279 =             DN(1,0)*rho*tau1;
        const double clhs280 =             DN(1,1)*rho*tau1;
        const double clhs281 =             DN(1,2)*rho*tau1;
        const double clhs282 =             clhs278 + clhs49;
        const double clhs283 =             DN(2,0)*clhs1*rho*tau1;
        const double clhs284 =             N[2]*rho;
        const double clhs285 =             clhs284*clhs5;
        const double clhs286 =             N[2]*clhs6*clhs7*tau1;
        const double clhs287 =             -clhs10*clhs286;
        const double clhs288 =             DN(2,1)*clhs18;
        const double clhs289 =             DN(2,2)*clhs18;
        const double clhs290 =             DN(2,0)*N[0];
        const double clhs291 =             rho*tau1*(DN(2,0)*clhs27 + DN(2,1)*clhs29 + DN(2,2)*clhs31);
        const double clhs292 =             clhs284*clhs42;
        const double clhs293 =             -clhs286*clhs44;
        const double clhs294 =             DN(2,1)*clhs49;
        const double clhs295 =             DN(2,2)*clhs49;
        const double clhs296 =             DN(2,0)*N[1];
        const double clhs297 =             rho*tau1*(DN(2,0)*clhs56 + DN(2,1)*clhs57 + DN(2,2)*clhs58);
        const double clhs298 =             pow(DN(2,0), 2)*rho;
        const double clhs299 =             pow(N[2], 2);
        const double clhs300 =             clhs12*clhs299 + clhs284*clhs68 - clhs286*clhs70 + clhs299*clhs8;
        const double clhs301 =             DN(2,0)*tau2;
        const double clhs302 =             DN(2,0) + clhs1*clhs75 - clhs301;
        const double clhs303 =             DN(2,0)*N[2];
        const double clhs304 =             clhs23*clhs299;
        const double clhs305 =             rho*tau1*(DN(2,0)*clhs82 + DN(2,1)*clhs83 + DN(2,2)*clhs84);
        const double clhs306 =             bdf0*clhs299*clhs6*tau1;
        const double clhs307 =             DN(3,0)*clhs65;
        const double clhs308 =             N[2]*N[3]*bdf0;
        const double clhs309 =             clhs308*rho;
        const double clhs310 =             -clhs301*clhs91;
        const double clhs311 =             N[3]*clhs69;
        const double clhs312 =             clhs284*clhs94;
        const double clhs313 =             -clhs286*clhs95;
        const double clhs314 =             DN(2,0)*DN(3,1);
        const double clhs315 =             -clhs314*tau2 + clhs314;
        const double clhs316 =             DN(2,1)*clhs100;
        const double clhs317 =             DN(2,0)*DN(3,2);
        const double clhs318 =             -clhs317*tau2 + clhs317;
        const double clhs319 =             DN(2,2)*clhs100;
        const double clhs320 =             DN(2,0)*N[3];
        const double clhs321 =             rho*tau1*(DN(2,0)*clhs107 + DN(2,1)*clhs108 + DN(2,2)*clhs109);
        const double clhs322 =             clhs308*clhs6;
        const double clhs323 =             N[2]*N[3]*bdf0*clhs6*tau1;
        const double clhs324 =             clhs25*clhs322 - clhs323*clhs34;
        const double clhs325 =             DN(2,0)*clhs18;
        const double clhs326 =             DN(2,1)*clhs3*rho*tau1;
        const double clhs327 =             DN(2,1)*N[0];
        const double clhs328 =             DN(2,0)*clhs49;
        const double clhs329 =             DN(2,1)*N[1];
        const double clhs330 =             DN(2,1)*tau2;
        const double clhs331 =             DN(2,1) + clhs3*clhs75 - clhs330;
        const double clhs332 =             pow(DN(2,1), 2)*rho;
        const double clhs333 =             DN(2,1)*N[2];
        const double clhs334 =             DN(2,1)*DN(3,0);
        const double clhs335 =             -clhs334*tau2 + clhs334;
        const double clhs336 =             DN(2,0)*clhs100;
        const double clhs337 =             DN(3,1)*clhs132 - clhs141*clhs330 + clhs309 + clhs311;
        const double clhs338 =             DN(2,1)*DN(3,2);
        const double clhs339 =             -clhs338*tau2 + clhs338;
        const double clhs340 =             DN(2,1)*N[3];
        const double clhs341 =             -clhs119*clhs323 + clhs28*clhs322;
        const double clhs342 =             DN(2,2)*clhs4*rho*tau1;
        const double clhs343 =             DN(2,2)*N[0];
        const double clhs344 =             DN(2,2)*N[1];
        const double clhs345 =             DN(2,2)*tau2;
        const double clhs346 =             DN(2,2) - clhs345 + clhs4*clhs75;
        const double clhs347 =             pow(DN(2,2), 2)*rho;
        const double clhs348 =             DN(2,2)*N[2];
        const double clhs349 =             DN(2,2)*DN(3,0);
        const double clhs350 =             -clhs349*tau2 + clhs349;
        const double clhs351 =             DN(2,2)*DN(3,1);
        const double clhs352 =             -clhs351*tau2 + clhs351;
        const double clhs353 =             DN(3,2)*clhs165 - clhs173*clhs345 + clhs309 + clhs311;
        const double clhs354 =             DN(2,2)*N[3];
        const double clhs355 =             -clhs152*clhs323 + clhs30*clhs322;
        const double clhs356 =             2*N[2];
        const double clhs357 =             DN(2,0)*rho*tau1;
        const double clhs358 =             DN(2,1)*rho*tau1;
        const double clhs359 =             DN(2,2)*rho*tau1;
        const double clhs360 =             clhs356 + clhs75;
        const double clhs361 =             DN(3,0)*clhs1*rho*tau1;
        const double clhs362 =             N[3]*rho;
        const double clhs363 =             clhs362*clhs5;
        const double clhs364 =             N[3]*clhs6*clhs7*tau1;
        const double clhs365 =             -clhs10*clhs364;
        const double clhs366 =             DN(3,1)*clhs18;
        const double clhs367 =             DN(3,2)*clhs18;
        const double clhs368 =             DN(3,0)*N[0];
        const double clhs369 =             rho*tau1*(DN(3,0)*clhs27 + DN(3,1)*clhs29 + DN(3,2)*clhs31);
        const double clhs370 =             clhs362*clhs42;
        const double clhs371 =             -clhs364*clhs44;
        const double clhs372 =             DN(3,1)*clhs49;
        const double clhs373 =             DN(3,2)*clhs49;
        const double clhs374 =             DN(3,0)*N[1];
        const double clhs375 =             rho*tau1*(DN(3,0)*clhs56 + DN(3,1)*clhs57 + DN(3,2)*clhs58);
        const double clhs376 =             clhs362*clhs68;
        const double clhs377 =             -clhs364*clhs70;
        const double clhs378 =             DN(3,1)*clhs75;
        const double clhs379 =             DN(3,2)*clhs75;
        const double clhs380 =             DN(3,0)*N[2];
        const double clhs381 =             rho*tau1*(DN(3,0)*clhs82 + DN(3,1)*clhs83 + DN(3,2)*clhs84);
        const double clhs382 =             pow(DN(3,0), 2)*rho;
        const double clhs383 =             pow(N[3], 2);
        const double clhs384 =             clhs12*clhs383 + clhs362*clhs94 - clhs364*clhs95 + clhs383*clhs8;
        const double clhs385 =             -DN(3,0)*tau2 + DN(3,0) + clhs1*clhs100;
        const double clhs386 =             DN(3,0)*N[3];
        const double clhs387 =             clhs23*clhs383;
        const double clhs388 =             rho*tau1*(DN(3,0)*clhs107 + DN(3,1)*clhs108 + DN(3,2)*clhs109);
        const double clhs389 =             bdf0*clhs383*clhs6*tau1;
        const double clhs390 =             DN(3,0)*clhs18;
        const double clhs391 =             DN(3,1)*clhs3*rho*tau1;
        const double clhs392 =             DN(3,1)*N[0];
        const double clhs393 =             DN(3,0)*clhs49;
        const double clhs394 =             DN(3,1)*N[1];
        const double clhs395 =             DN(3,0)*clhs75;
        const double clhs396 =             DN(3,1)*N[2];
        const double clhs397 =             -DN(3,1)*tau2 + DN(3,1) + clhs100*clhs3;
        const double clhs398 =             pow(DN(3,1), 2)*rho;
        const double clhs399 =             DN(3,1)*N[3];
        const double clhs400 =             DN(3,2)*clhs4*rho*tau1;
        const double clhs401 =             DN(3,2)*N[0];
        const double clhs402 =             DN(3,2)*N[1];
        const double clhs403 =             DN(3,2)*N[2];
        const double clhs404 =             -DN(3,2)*tau2 + DN(3,2) + clhs100*clhs4;
        const double clhs405 =             pow(DN(3,2), 2)*rho;
        const double clhs406 =             DN(3,2)*N[3];
        const double clhs407 =             2*N[3];
        const double clhs408 =             DN(3,0)*rho*tau1;
        const double clhs409 =             DN(3,1)*rho*tau1;
        const double clhs410 =             DN(3,2)*rho*tau1;
        const double clhs411 =             clhs100 + clhs407;

        lhs(0,0)=-clhs0*tau2 + clhs0 + clhs10*clhs2 + clhs15;
        lhs(0,1)=clhs16*clhs19;
        lhs(0,2)=clhs19*clhs20;
        lhs(0,3)=clhs1*clhs32 - clhs14*clhs27 - clhs21*clhs22 - clhs21 + clhs24*clhs25 - clhs33*clhs34;
        lhs(0,4)=clhs2*clhs44 + clhs36 + clhs38 + clhs40 + clhs41 + clhs45 + clhs46;
        lhs(0,5)=rho*(clhs1*clhs50 + clhs48);
        lhs(0,6)=rho*(clhs1*clhs53 + clhs52);
        lhs(0,7)=clhs1*clhs59 - clhs14*clhs56 - clhs22*clhs54 - clhs54 + clhs62;
        lhs(0,8)=clhs2*clhs70 + clhs63 + clhs64 + clhs66 + clhs67 + clhs71 + clhs72;
        lhs(0,9)=rho*(clhs1*clhs76 + clhs74);
        lhs(0,10)=rho*(clhs1*clhs79 + clhs78);
        lhs(0,11)=clhs1*clhs85 - clhs14*clhs82 - clhs22*clhs80 - clhs80 + clhs88;
        lhs(0,12)=clhs2*clhs95 + clhs89 + clhs90 + clhs92 + clhs93 + clhs96 + clhs97;
        lhs(0,13)=rho*(clhs1*clhs101 + clhs99);
        lhs(0,14)=rho*(clhs1*clhs104 + clhs103);
        lhs(0,15)=clhs1*clhs110 - clhs105*clhs22 - clhs105 - clhs107*clhs14 + clhs113;
        lhs(1,0)=clhs115*clhs35;
        lhs(1,1)=clhs10*clhs117 - clhs116*tau2 + clhs116 + clhs15;
        lhs(1,2)=clhs115*clhs20;
        lhs(1,3)=-clhs118*clhs22 - clhs118 - clhs119*clhs33 - clhs14*clhs29 + clhs24*clhs28 + clhs3*clhs32;
        lhs(1,4)=rho*(clhs121 + clhs122*clhs3);
        lhs(1,5)=clhs117*clhs44 + clhs124 + clhs45 + clhs46;
        lhs(1,6)=rho*(clhs126 + clhs3*clhs53);
        lhs(1,7)=-clhs127*clhs22 - clhs127 + clhs128 - clhs14*clhs57 + clhs3*clhs59;
        lhs(1,8)=rho*(clhs130 + clhs131*clhs3);
        lhs(1,9)=clhs117*clhs70 + clhs133 + clhs71 + clhs72;
        lhs(1,10)=rho*(clhs135 + clhs3*clhs79);
        lhs(1,11)=-clhs136*clhs22 - clhs136 + clhs137 - clhs14*clhs83 + clhs3*clhs85;
        lhs(1,12)=rho*(clhs139 + clhs140*clhs3);
        lhs(1,13)=clhs117*clhs95 + clhs142 + clhs96 + clhs97;
        lhs(1,14)=rho*(clhs104*clhs3 + clhs144);
        lhs(1,15)=-clhs108*clhs14 + clhs110*clhs3 - clhs145*clhs22 - clhs145 + clhs146;
        lhs(2,0)=clhs148*clhs35;
        lhs(2,1)=clhs148*clhs16;
        lhs(2,2)=clhs10*clhs150 - clhs149*tau2 + clhs149 + clhs15;
        lhs(2,3)=-clhs14*clhs31 - clhs151*clhs22 - clhs151 - clhs152*clhs33 + clhs24*clhs30 + clhs32*clhs4;
        lhs(2,4)=rho*(clhs122*clhs4 + clhs154);
        lhs(2,5)=rho*(clhs156 + clhs4*clhs50);
        lhs(2,6)=clhs150*clhs44 + clhs158 + clhs45 + clhs46;
        lhs(2,7)=-clhs14*clhs58 - clhs159*clhs22 - clhs159 + clhs160 + clhs4*clhs59;
        lhs(2,8)=rho*(clhs131*clhs4 + clhs162);
        lhs(2,9)=rho*(clhs164 + clhs4*clhs76);
        lhs(2,10)=clhs150*clhs70 + clhs166 + clhs71 + clhs72;
        lhs(2,11)=-clhs14*clhs84 - clhs167*clhs22 - clhs167 + clhs168 + clhs4*clhs85;
        lhs(2,12)=rho*(clhs140*clhs4 + clhs170);
        lhs(2,13)=rho*(clhs101*clhs4 + clhs172);
        lhs(2,14)=clhs150*clhs95 + clhs174 + clhs96 + clhs97;
        lhs(2,15)=-clhs109*clhs14 + clhs110*clhs4 - clhs175*clhs22 - clhs175 + clhs176;
        lhs(3,0)=clhs178*clhs35;
        lhs(3,1)=clhs16*clhs178;
        lhs(3,2)=clhs178*clhs20;
        lhs(3,3)=clhs179*clhs27 + clhs180*clhs29 + clhs181*clhs31 + clhs24;
        lhs(3,4)=rho*(DN(1,0)*clhs177 + clhs122);
        lhs(3,5)=rho*(DN(1,1)*clhs177 + clhs50);
        lhs(3,6)=rho*(DN(1,2)*clhs177 + clhs53);
        lhs(3,7)=clhs179*clhs56 + clhs180*clhs57 + clhs181*clhs58 + clhs60;
        lhs(3,8)=rho*(DN(2,0)*clhs177 + clhs131);
        lhs(3,9)=rho*(DN(2,1)*clhs177 + clhs76);
        lhs(3,10)=rho*(DN(2,2)*clhs177 + clhs79);
        lhs(3,11)=clhs179*clhs82 + clhs180*clhs83 + clhs181*clhs84 + clhs86;
        lhs(3,12)=rho*(DN(3,0)*clhs177 + clhs140);
        lhs(3,13)=rho*(DN(3,1)*clhs177 + clhs101);
        lhs(3,14)=rho*(DN(3,2)*clhs177 + clhs104);
        lhs(3,15)=clhs107*clhs179 + clhs108*clhs180 + clhs109*clhs181 + clhs111;
        lhs(4,0)=clhs10*clhs182 + clhs184 + clhs186 + clhs36 + clhs38 + clhs40 + clhs41;
        lhs(4,1)=rho*(clhs1*clhs187 + clhs121);
        lhs(4,2)=rho*(clhs1*clhs188 + clhs154);
        lhs(4,3)=clhs1*clhs190 - clhs185*clhs27 - clhs189*clhs22 - clhs189 + clhs62;
        lhs(4,4)=clhs182*clhs44 - clhs191*tau2 + clhs191 + clhs193;
        lhs(4,5)=clhs123*clhs195;
        lhs(4,6)=clhs157*clhs195;
        lhs(4,7)=clhs1*clhs198 - clhs185*clhs56 - clhs196*clhs22 - clhs196 + clhs197*clhs25 - clhs199*clhs34;
        lhs(4,8)=clhs182*clhs70 + clhs200 + clhs202 + clhs203 + clhs204 + clhs205 + clhs206;
        lhs(4,9)=rho*(clhs1*clhs209 + clhs208);
        lhs(4,10)=rho*(clhs1*clhs212 + clhs211);
        lhs(4,11)=clhs1*clhs214 - clhs185*clhs82 - clhs213*clhs22 - clhs213 + clhs217;
        lhs(4,12)=clhs182*clhs95 + clhs218 + clhs219 + clhs220 + clhs221 + clhs222 + clhs223;
        lhs(4,13)=rho*(clhs1*clhs226 + clhs225);
        lhs(4,14)=rho*(clhs1*clhs229 + clhs228);
        lhs(4,15)=clhs1*clhs231 - clhs107*clhs185 - clhs22*clhs230 - clhs230 + clhs234;
        lhs(5,0)=rho*(clhs235*clhs3 + clhs48);
        lhs(5,1)=clhs10*clhs236 + clhs124 + clhs184 + clhs186;
        lhs(5,2)=rho*(clhs156 + clhs188*clhs3);
        lhs(5,3)=clhs128 - clhs185*clhs29 + clhs190*clhs3 - clhs22*clhs237 - clhs237;
        lhs(5,4)=clhs239*clhs39;
        lhs(5,5)=clhs193 + clhs236*clhs44 - clhs240*tau2 + clhs240;
        lhs(5,6)=clhs157*clhs239;
        lhs(5,7)=-clhs119*clhs199 - clhs185*clhs57 + clhs197*clhs28 + clhs198*clhs3 - clhs22*clhs241 - clhs241;
        lhs(5,8)=rho*(clhs243 + clhs244*clhs3);
        lhs(5,9)=clhs205 + clhs206 + clhs236*clhs70 + clhs245;
        lhs(5,10)=rho*(clhs212*clhs3 + clhs247);
        lhs(5,11)=-clhs185*clhs83 + clhs214*clhs3 - clhs22*clhs248 - clhs248 + clhs249;
        lhs(5,12)=rho*(clhs251 + clhs252*clhs3);
        lhs(5,13)=clhs222 + clhs223 + clhs236*clhs95 + clhs253;
        lhs(5,14)=rho*(clhs229*clhs3 + clhs255);
        lhs(5,15)=-clhs108*clhs185 - clhs22*clhs256 + clhs231*clhs3 - clhs256 + clhs257;
        lhs(6,0)=rho*(clhs235*clhs4 + clhs52);
        lhs(6,1)=rho*(clhs126 + clhs187*clhs4);
        lhs(6,2)=clhs10*clhs258 + clhs158 + clhs184 + clhs186;
        lhs(6,3)=clhs160 - clhs185*clhs31 + clhs190*clhs4 - clhs22*clhs259 - clhs259;
        lhs(6,4)=clhs261*clhs39;
        lhs(6,5)=clhs123*clhs261;
        lhs(6,6)=clhs193 + clhs258*clhs44 - clhs262*tau2 + clhs262;
        lhs(6,7)=-clhs152*clhs199 - clhs185*clhs58 + clhs197*clhs30 + clhs198*clhs4 - clhs22*clhs263 - clhs263;
        lhs(6,8)=rho*(clhs244*clhs4 + clhs265);
        lhs(6,9)=rho*(clhs209*clhs4 + clhs267);
        lhs(6,10)=clhs205 + clhs206 + clhs258*clhs70 + clhs268;
        lhs(6,11)=-clhs185*clhs84 + clhs214*clhs4 - clhs22*clhs269 - clhs269 + clhs270;
        lhs(6,12)=rho*(clhs252*clhs4 + clhs272);
        lhs(6,13)=rho*(clhs226*clhs4 + clhs274);
        lhs(6,14)=clhs222 + clhs223 + clhs258*clhs95 + clhs275;
        lhs(6,15)=-clhs109*clhs185 - clhs22*clhs276 + clhs231*clhs4 - clhs276 + clhs277;
        lhs(7,0)=rho*(DN(0,0)*clhs278 + clhs235);
        lhs(7,1)=rho*(DN(0,1)*clhs278 + clhs187);
        lhs(7,2)=rho*(DN(0,2)*clhs278 + clhs188);
        lhs(7,3)=clhs27*clhs279 + clhs280*clhs29 + clhs281*clhs31 + clhs60;
        lhs(7,4)=clhs282*clhs39;
        lhs(7,5)=clhs123*clhs282;
        lhs(7,6)=clhs157*clhs282;
        lhs(7,7)=clhs197 + clhs279*clhs56 + clhs280*clhs57 + clhs281*clhs58;
        lhs(7,8)=rho*(DN(2,0)*clhs278 + clhs244);
        lhs(7,9)=rho*(DN(2,1)*clhs278 + clhs209);
        lhs(7,10)=rho*(DN(2,2)*clhs278 + clhs212);
        lhs(7,11)=clhs215 + clhs279*clhs82 + clhs280*clhs83 + clhs281*clhs84;
        lhs(7,12)=rho*(DN(3,0)*clhs278 + clhs252);
        lhs(7,13)=rho*(DN(3,1)*clhs278 + clhs226);
        lhs(7,14)=rho*(DN(3,2)*clhs278 + clhs229);
        lhs(7,15)=clhs107*clhs279 + clhs108*clhs280 + clhs109*clhs281 + clhs232;
        lhs(8,0)=clhs10*clhs283 + clhs285 + clhs287 + clhs63 + clhs64 + clhs66 + clhs67;
        lhs(8,1)=rho*(clhs1*clhs288 + clhs130);
        lhs(8,2)=rho*(clhs1*clhs289 + clhs162);
        lhs(8,3)=clhs1*clhs291 - clhs22*clhs290 - clhs27*clhs286 - clhs290 + clhs88;
        lhs(8,4)=clhs200 + clhs202 + clhs203 + clhs204 + clhs283*clhs44 + clhs292 + clhs293;
        lhs(8,5)=rho*(clhs1*clhs294 + clhs243);
        lhs(8,6)=rho*(clhs1*clhs295 + clhs265);
        lhs(8,7)=clhs1*clhs297 + clhs217 - clhs22*clhs296 - clhs286*clhs56 - clhs296;
        lhs(8,8)=clhs283*clhs70 - clhs298*tau2 + clhs298 + clhs300;
        lhs(8,9)=clhs132*clhs302;
        lhs(8,10)=clhs165*clhs302;
        lhs(8,11)=clhs1*clhs305 - clhs22*clhs303 + clhs25*clhs304 - clhs286*clhs82 - clhs303 - clhs306*clhs34;
        lhs(8,12)=clhs283*clhs95 + clhs307 + clhs309 + clhs310 + clhs311 + clhs312 + clhs313;
        lhs(8,13)=rho*(clhs1*clhs316 + clhs315);
        lhs(8,14)=rho*(clhs1*clhs319 + clhs318);
        lhs(8,15)=clhs1*clhs321 - clhs107*clhs286 - clhs22*clhs320 - clhs320 + clhs324;
        lhs(9,0)=rho*(clhs3*clhs325 + clhs74);
        lhs(9,1)=clhs10*clhs326 + clhs133 + clhs285 + clhs287;
        lhs(9,2)=rho*(clhs164 + clhs289*clhs3);
        lhs(9,3)=clhs137 - clhs22*clhs327 - clhs286*clhs29 + clhs291*clhs3 - clhs327;
        lhs(9,4)=rho*(clhs208 + clhs3*clhs328);
        lhs(9,5)=clhs245 + clhs292 + clhs293 + clhs326*clhs44;
        lhs(9,6)=rho*(clhs267 + clhs295*clhs3);
        lhs(9,7)=-clhs22*clhs329 + clhs249 - clhs286*clhs57 + clhs297*clhs3 - clhs329;
        lhs(9,8)=clhs331*clhs65;
        lhs(9,9)=clhs300 + clhs326*clhs70 - clhs332*tau2 + clhs332;
        lhs(9,10)=clhs165*clhs331;
        lhs(9,11)=-clhs119*clhs306 - clhs22*clhs333 + clhs28*clhs304 - clhs286*clhs83 + clhs3*clhs305 - clhs333;
        lhs(9,12)=rho*(clhs3*clhs336 + clhs335);
        lhs(9,13)=clhs312 + clhs313 + clhs326*clhs95 + clhs337;
        lhs(9,14)=rho*(clhs3*clhs319 + clhs339);
        lhs(9,15)=-clhs108*clhs286 - clhs22*clhs340 + clhs3*clhs321 - clhs340 + clhs341;
        lhs(10,0)=rho*(clhs325*clhs4 + clhs78);
        lhs(10,1)=rho*(clhs135 + clhs288*clhs4);
        lhs(10,2)=clhs10*clhs342 + clhs166 + clhs285 + clhs287;
        lhs(10,3)=clhs168 - clhs22*clhs343 - clhs286*clhs31 + clhs291*clhs4 - clhs343;
        lhs(10,4)=rho*(clhs211 + clhs328*clhs4);
        lhs(10,5)=rho*(clhs247 + clhs294*clhs4);
        lhs(10,6)=clhs268 + clhs292 + clhs293 + clhs342*clhs44;
        lhs(10,7)=-clhs22*clhs344 + clhs270 - clhs286*clhs58 + clhs297*clhs4 - clhs344;
        lhs(10,8)=clhs346*clhs65;
        lhs(10,9)=clhs132*clhs346;
        lhs(10,10)=clhs300 + clhs342*clhs70 - clhs347*tau2 + clhs347;
        lhs(10,11)=-clhs152*clhs306 - clhs22*clhs348 - clhs286*clhs84 + clhs30*clhs304 + clhs305*clhs4 - clhs348;
        lhs(10,12)=rho*(clhs336*clhs4 + clhs350);
        lhs(10,13)=rho*(clhs316*clhs4 + clhs352);
        lhs(10,14)=clhs312 + clhs313 + clhs342*clhs95 + clhs353;
        lhs(10,15)=-clhs109*clhs286 - clhs22*clhs354 + clhs321*clhs4 - clhs354 + clhs355;
        lhs(11,0)=rho*(DN(0,0)*clhs356 + clhs325);
        lhs(11,1)=rho*(DN(0,1)*clhs356 + clhs288);
        lhs(11,2)=rho*(DN(0,2)*clhs356 + clhs289);
        lhs(11,3)=clhs27*clhs357 + clhs29*clhs358 + clhs31*clhs359 + clhs86;
        lhs(11,4)=rho*(DN(1,0)*clhs356 + clhs328);
        lhs(11,5)=rho*(DN(1,1)*clhs356 + clhs294);
        lhs(11,6)=rho*(DN(1,2)*clhs356 + clhs295);
        lhs(11,7)=clhs215 + clhs357*clhs56 + clhs358*clhs57 + clhs359*clhs58;
        lhs(11,8)=clhs360*clhs65;
        lhs(11,9)=clhs132*clhs360;
        lhs(11,10)=clhs165*clhs360;
        lhs(11,11)=clhs304 + clhs357*clhs82 + clhs358*clhs83 + clhs359*clhs84;
        lhs(11,12)=rho*(DN(3,0)*clhs356 + clhs336);
        lhs(11,13)=rho*(DN(3,1)*clhs356 + clhs316);
        lhs(11,14)=rho*(DN(3,2)*clhs356 + clhs319);
        lhs(11,15)=clhs107*clhs357 + clhs108*clhs358 + clhs109*clhs359 + clhs322;
        lhs(12,0)=clhs10*clhs361 + clhs363 + clhs365 + clhs89 + clhs90 + clhs92 + clhs93;
        lhs(12,1)=rho*(clhs1*clhs366 + clhs139);
        lhs(12,2)=rho*(clhs1*clhs367 + clhs170);
        lhs(12,3)=clhs1*clhs369 + clhs113 - clhs22*clhs368 - clhs27*clhs364 - clhs368;
        lhs(12,4)=clhs218 + clhs219 + clhs220 + clhs221 + clhs361*clhs44 + clhs370 + clhs371;
        lhs(12,5)=rho*(clhs1*clhs372 + clhs251);
        lhs(12,6)=rho*(clhs1*clhs373 + clhs272);
        lhs(12,7)=clhs1*clhs375 - clhs22*clhs374 + clhs234 - clhs364*clhs56 - clhs374;
        lhs(12,8)=clhs307 + clhs309 + clhs310 + clhs311 + clhs361*clhs70 + clhs376 + clhs377;
        lhs(12,9)=rho*(clhs1*clhs378 + clhs335);
        lhs(12,10)=rho*(clhs1*clhs379 + clhs350);
        lhs(12,11)=clhs1*clhs381 - clhs22*clhs380 + clhs324 - clhs364*clhs82 - clhs380;
        lhs(12,12)=clhs361*clhs95 - clhs382*tau2 + clhs382 + clhs384;
        lhs(12,13)=clhs141*clhs385;
        lhs(12,14)=clhs173*clhs385;
        lhs(12,15)=clhs1*clhs388 - clhs107*clhs364 - clhs22*clhs386 + clhs25*clhs387 - clhs34*clhs389 - clhs386;
        lhs(13,0)=rho*(clhs3*clhs390 + clhs99);
        lhs(13,1)=clhs10*clhs391 + clhs142 + clhs363 + clhs365;
        lhs(13,2)=rho*(clhs172 + clhs3*clhs367);
        lhs(13,3)=clhs146 - clhs22*clhs392 - clhs29*clhs364 + clhs3*clhs369 - clhs392;
        lhs(13,4)=rho*(clhs225 + clhs3*clhs393);
        lhs(13,5)=clhs253 + clhs370 + clhs371 + clhs391*clhs44;
        lhs(13,6)=rho*(clhs274 + clhs3*clhs373);
        lhs(13,7)=-clhs22*clhs394 + clhs257 + clhs3*clhs375 - clhs364*clhs57 - clhs394;
        lhs(13,8)=rho*(clhs3*clhs395 + clhs315);
        lhs(13,9)=clhs337 + clhs376 + clhs377 + clhs391*clhs70;
        lhs(13,10)=rho*(clhs3*clhs379 + clhs352);
        lhs(13,11)=-clhs22*clhs396 + clhs3*clhs381 + clhs341 - clhs364*clhs83 - clhs396;
        lhs(13,12)=clhs397*clhs91;
        lhs(13,13)=clhs384 + clhs391*clhs95 - clhs398*tau2 + clhs398;
        lhs(13,14)=clhs173*clhs397;
        lhs(13,15)=-clhs108*clhs364 - clhs119*clhs389 - clhs22*clhs399 + clhs28*clhs387 + clhs3*clhs388 - clhs399;
        lhs(14,0)=rho*(clhs103 + clhs390*clhs4);
        lhs(14,1)=rho*(clhs144 + clhs366*clhs4);
        lhs(14,2)=clhs10*clhs400 + clhs174 + clhs363 + clhs365;
        lhs(14,3)=clhs176 - clhs22*clhs401 - clhs31*clhs364 + clhs369*clhs4 - clhs401;
        lhs(14,4)=rho*(clhs228 + clhs393*clhs4);
        lhs(14,5)=rho*(clhs255 + clhs372*clhs4);
        lhs(14,6)=clhs275 + clhs370 + clhs371 + clhs400*clhs44;
        lhs(14,7)=-clhs22*clhs402 + clhs277 - clhs364*clhs58 + clhs375*clhs4 - clhs402;
        lhs(14,8)=rho*(clhs318 + clhs395*clhs4);
        lhs(14,9)=rho*(clhs339 + clhs378*clhs4);
        lhs(14,10)=clhs353 + clhs376 + clhs377 + clhs400*clhs70;
        lhs(14,11)=-clhs22*clhs403 + clhs355 - clhs364*clhs84 + clhs381*clhs4 - clhs403;
        lhs(14,12)=clhs404*clhs91;
        lhs(14,13)=clhs141*clhs404;
        lhs(14,14)=clhs384 + clhs400*clhs95 - clhs405*tau2 + clhs405;
        lhs(14,15)=-clhs109*clhs364 - clhs152*clhs389 - clhs22*clhs406 + clhs30*clhs387 + clhs388*clhs4 - clhs406;
        lhs(15,0)=rho*(DN(0,0)*clhs407 + clhs390);
        lhs(15,1)=rho*(DN(0,1)*clhs407 + clhs366);
        lhs(15,2)=rho*(DN(0,2)*clhs407 + clhs367);
        lhs(15,3)=clhs111 + clhs27*clhs408 + clhs29*clhs409 + clhs31*clhs410;
        lhs(15,4)=rho*(DN(1,0)*clhs407 + clhs393);
        lhs(15,5)=rho*(DN(1,1)*clhs407 + clhs372);
        lhs(15,6)=rho*(DN(1,2)*clhs407 + clhs373);
        lhs(15,7)=clhs232 + clhs408*clhs56 + clhs409*clhs57 + clhs410*clhs58;
        lhs(15,8)=rho*(DN(2,0)*clhs407 + clhs395);
        lhs(15,9)=rho*(DN(2,1)*clhs407 + clhs378);
        lhs(15,10)=rho*(DN(2,2)*clhs407 + clhs379);
        lhs(15,11)=clhs322 + clhs408*clhs82 + clhs409*clhs83 + clhs410*clhs84;
        lhs(15,12)=clhs411*clhs91;
        lhs(15,13)=clhs141*clhs411;
        lhs(15,14)=clhs173*clhs411;
        lhs(15,15)=clhs107*clhs408 + clhs108*clhs409 + clhs109*clhs410 + clhs387;

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

        const double clhs0 =             pow(DN(0,0), 2)*rho;
        const double clhs1 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
        const double clhs2 =             DN(0,0)*clhs1*rho*tau1;
        const double clhs3 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
        const double clhs4 =             DN(0,0)*clhs1 + DN(0,1)*clhs3;
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
        const double clhs15 =             DN(0,1)*rho;
        const double clhs16 =             DN(0,0)*tau2;
        const double clhs17 =             clhs9*tau1;
        const double clhs18 =             DN(0,0)*N[0];
        const double clhs19 =             bdf0*clhs5*tau2;
        const double clhs20 =             bdf0*clhs5;
        const double clhs21 =             clhs10*clhs20;
        const double clhs22 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0);
        const double clhs23 =             N[0]*bdf0*clhs5;
        const double clhs24 =             DN(0,0) + clhs22*clhs23;
        const double clhs25 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
        const double clhs26 =             DN(0,1) + clhs23*clhs25;
        const double clhs27 =             rho*tau1*(DN(0,0)*clhs24 + DN(0,1)*clhs26);
        const double clhs28 =             bdf0*clhs10*clhs5*tau1;
        const double clhs29 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + clhs22*clhs7 - rho*(N[0]*f(0,0) + N[1]*f(1,0) + N[2]*f(2,0)) + rho*(N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + clhs1*(DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0)) + clhs3*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0)));
        const double clhs30 =             DN(1,0)*clhs1 + DN(1,1)*clhs3;
        const double clhs31 =             clhs12*clhs30;
        const double clhs32 =             N[1]*clhs7;
        const double clhs33 =             clhs32 + rho*(N[1]*bdf0 + clhs30);
        const double clhs34 =             -clhs13*clhs33;
        const double clhs35 =             DN(0,0)*rho;
        const double clhs36 =             DN(1,0)*clhs35;
        const double clhs37 =             N[0]*bdf0*rho;
        const double clhs38 =             N[1]*clhs37;
        const double clhs39 =             DN(1,0)*rho;
        const double clhs40 =             -clhs16*clhs39;
        const double clhs41 =             N[1]*clhs8;
        const double clhs42 =             DN(0,0)*DN(1,1);
        const double clhs43 =             -clhs42*tau2 + clhs42;
        const double clhs44 =             clhs33*tau1;
        const double clhs45 =             DN(0,1)*clhs44;
        const double clhs46 =             DN(0,0)*N[1];
        const double clhs47 =             N[1]*bdf0*clhs5;
        const double clhs48 =             DN(1,0) + clhs22*clhs47;
        const double clhs49 =             DN(1,1) + clhs25*clhs47;
        const double clhs50 =             rho*tau1*(DN(0,0)*clhs48 + DN(0,1)*clhs49);
        const double clhs51 =             N[1]*clhs23;
        const double clhs52 =             N[0]*N[1]*bdf0*clhs5*tau1;
        const double clhs53 =             clhs22*clhs51 - clhs29*clhs52;
        const double clhs54 =             DN(2,0)*clhs1 + DN(2,1)*clhs3;
        const double clhs55 =             clhs12*clhs54;
        const double clhs56 =             N[2]*clhs7 + rho*(N[2]*bdf0 + clhs54);
        const double clhs57 =             -clhs13*clhs56;
        const double clhs58 =             DN(2,0)*clhs35;
        const double clhs59 =             N[2]*clhs37;
        const double clhs60 =             DN(2,0)*rho;
        const double clhs61 =             -clhs16*clhs60;
        const double clhs62 =             N[2]*clhs8;
        const double clhs63 =             DN(0,0)*DN(2,1);
        const double clhs64 =             -clhs63*tau2 + clhs63;
        const double clhs65 =             clhs56*tau1;
        const double clhs66 =             DN(0,1)*clhs65;
        const double clhs67 =             DN(0,0)*N[2];
        const double clhs68 =             N[2]*bdf0*clhs5;
        const double clhs69 =             DN(2,0) + clhs22*clhs68;
        const double clhs70 =             DN(2,1) + clhs25*clhs68;
        const double clhs71 =             rho*tau1*(DN(0,0)*clhs69 + DN(0,1)*clhs70);
        const double clhs72 =             N[2]*clhs23;
        const double clhs73 =             N[0]*N[2]*bdf0*clhs5*tau1;
        const double clhs74 =             clhs22*clhs72 - clhs29*clhs73;
        const double clhs75 =             DN(0,1)*tau2;
        const double clhs76 =             pow(DN(0,1), 2)*rho;
        const double clhs77 =             DN(0,1)*clhs3*rho*tau1;
        const double clhs78 =             DN(0,1)*N[0];
        const double clhs79 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + clhs25*clhs7 - rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1)) + rho*(N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + clhs1*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + clhs3*(DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1)));
        const double clhs80 =             DN(0,1)*DN(1,0);
        const double clhs81 =             -clhs80*tau2 + clhs80;
        const double clhs82 =             DN(0,0)*clhs44;
        const double clhs83 =             DN(1,1)*rho;
        const double clhs84 =             DN(1,1)*clhs15 + clhs38 + clhs41 - clhs75*clhs83;
        const double clhs85 =             DN(0,1)*N[1];
        const double clhs86 =             clhs25*clhs51 - clhs52*clhs79;
        const double clhs87 =             DN(0,1)*DN(2,0);
        const double clhs88 =             -clhs87*tau2 + clhs87;
        const double clhs89 =             DN(0,0)*clhs65;
        const double clhs90 =             DN(2,1)*rho;
        const double clhs91 =             DN(2,1)*clhs15 + clhs59 + clhs62 - clhs75*clhs90;
        const double clhs92 =             DN(0,1)*N[2];
        const double clhs93 =             clhs25*clhs72 - clhs73*clhs79;
        const double clhs94 =             2*N[0];
        const double clhs95 =             clhs17 + clhs94;
        const double clhs96 =             DN(0,0)*rho*tau1;
        const double clhs97 =             DN(0,1)*rho*tau1;
        const double clhs98 =             N[1]*rho;
        const double clhs99 =             clhs4*clhs98;
        const double clhs100 =             N[1]*clhs5*clhs6*tau1;
        const double clhs101 =             -clhs100*clhs9;
        const double clhs102 =             DN(1,0)*clhs1*rho*tau1;
        const double clhs103 =             DN(1,1)*clhs17;
        const double clhs104 =             DN(1,0)*N[0];
        const double clhs105 =             rho*tau1*(DN(1,0)*clhs24 + DN(1,1)*clhs26);
        const double clhs106 =             pow(DN(1,0), 2)*rho;
        const double clhs107 =             pow(N[1], 2);
        const double clhs108 =             -clhs100*clhs33 + clhs107*clhs11 + clhs107*clhs7 + clhs30*clhs98;
        const double clhs109 =             DN(1,0)*tau2;
        const double clhs110 =             DN(1,0)*N[1];
        const double clhs111 =             clhs107*clhs20;
        const double clhs112 =             rho*tau1*(DN(1,0)*clhs48 + DN(1,1)*clhs49);
        const double clhs113 =             bdf0*clhs107*clhs5*tau1;
        const double clhs114 =             clhs54*clhs98;
        const double clhs115 =             -clhs100*clhs56;
        const double clhs116 =             DN(2,0)*clhs39;
        const double clhs117 =             N[1]*N[2]*bdf0;
        const double clhs118 =             clhs117*rho;
        const double clhs119 =             -clhs109*clhs60;
        const double clhs120 =             N[2]*clhs32;
        const double clhs121 =             DN(1,0)*DN(2,1);
        const double clhs122 =             -clhs121*tau2 + clhs121;
        const double clhs123 =             DN(1,1)*clhs65;
        const double clhs124 =             DN(1,0)*N[2];
        const double clhs125 =             rho*tau1*(DN(1,0)*clhs69 + DN(1,1)*clhs70);
        const double clhs126 =             clhs117*clhs5;
        const double clhs127 =             N[1]*N[2]*bdf0*clhs5*tau1;
        const double clhs128 =             clhs126*clhs22 - clhs127*clhs29;
        const double clhs129 =             DN(1,0)*clhs17;
        const double clhs130 =             DN(1,1)*clhs3*rho*tau1;
        const double clhs131 =             DN(1,1)*N[0];
        const double clhs132 =             DN(1,1)*tau2;
        const double clhs133 =             pow(DN(1,1), 2)*rho;
        const double clhs134 =             DN(1,1)*N[1];
        const double clhs135 =             DN(1,1)*DN(2,0);
        const double clhs136 =             -clhs135*tau2 + clhs135;
        const double clhs137 =             DN(1,0)*clhs65;
        const double clhs138 =             DN(2,1)*clhs83 + clhs118 + clhs120 - clhs132*clhs90;
        const double clhs139 =             DN(1,1)*N[2];
        const double clhs140 =             clhs126*clhs25 - clhs127*clhs79;
        const double clhs141 =             2*N[1];
        const double clhs142 =             DN(1,0)*rho*tau1;
        const double clhs143 =             DN(1,1)*rho*tau1;
        const double clhs144 =             clhs141 + clhs44;
        const double clhs145 =             N[2]*rho;
        const double clhs146 =             clhs145*clhs4;
        const double clhs147 =             N[2]*clhs5*clhs6*tau1;
        const double clhs148 =             -clhs147*clhs9;
        const double clhs149 =             DN(2,0)*clhs1*rho*tau1;
        const double clhs150 =             DN(2,1)*clhs17;
        const double clhs151 =             DN(2,0)*N[0];
        const double clhs152 =             rho*tau1*(DN(2,0)*clhs24 + DN(2,1)*clhs26);
        const double clhs153 =             clhs145*clhs30;
        const double clhs154 =             -clhs147*clhs33;
        const double clhs155 =             DN(2,1)*clhs44;
        const double clhs156 =             DN(2,0)*N[1];
        const double clhs157 =             rho*tau1*(DN(2,0)*clhs48 + DN(2,1)*clhs49);
        const double clhs158 =             pow(DN(2,0), 2)*rho;
        const double clhs159 =             pow(N[2], 2);
        const double clhs160 =             clhs11*clhs159 + clhs145*clhs54 - clhs147*clhs56 + clhs159*clhs7;
        const double clhs161 =             DN(2,0)*N[2];
        const double clhs162 =             clhs159*clhs20;
        const double clhs163 =             rho*tau1*(DN(2,0)*clhs69 + DN(2,1)*clhs70);
        const double clhs164 =             bdf0*clhs159*clhs5*tau1;
        const double clhs165 =             DN(2,0)*clhs17;
        const double clhs166 =             DN(2,1)*clhs3*rho*tau1;
        const double clhs167 =             DN(2,1)*N[0];
        const double clhs168 =             DN(2,0)*clhs44;
        const double clhs169 =             DN(2,1)*N[1];
        const double clhs170 =             pow(DN(2,1), 2)*rho;
        const double clhs171 =             DN(2,1)*N[2];
        const double clhs172 =             2*N[2];
        const double clhs173 =             DN(2,0)*rho*tau1;
        const double clhs174 =             DN(2,1)*rho*tau1;
        const double clhs175 =             clhs172 + clhs65;

        lhs(0,0)=-clhs0*tau2 + clhs0 + clhs14 + clhs2*clhs9;
        lhs(0,1)=clhs15*(DN(0,0) + clhs1*clhs17 - clhs16);
        lhs(0,2)=clhs1*clhs27 - clhs13*clhs24 - clhs18*clhs19 - clhs18 + clhs21*clhs22 - clhs28*clhs29;
        lhs(0,3)=clhs2*clhs33 + clhs31 + clhs34 + clhs36 + clhs38 + clhs40 + clhs41;
        lhs(0,4)=rho*(clhs1*clhs45 + clhs43);
        lhs(0,5)=clhs1*clhs50 - clhs13*clhs48 - clhs19*clhs46 - clhs46 + clhs53;
        lhs(0,6)=clhs2*clhs56 + clhs55 + clhs57 + clhs58 + clhs59 + clhs61 + clhs62;
        lhs(0,7)=rho*(clhs1*clhs66 + clhs64);
        lhs(0,8)=clhs1*clhs71 - clhs13*clhs69 - clhs19*clhs67 - clhs67 + clhs74;
        lhs(1,0)=clhs35*(DN(0,1) + clhs17*clhs3 - clhs75);
        lhs(1,1)=clhs14 - clhs76*tau2 + clhs76 + clhs77*clhs9;
        lhs(1,2)=-clhs13*clhs26 - clhs19*clhs78 + clhs21*clhs25 + clhs27*clhs3 - clhs28*clhs79 - clhs78;
        lhs(1,3)=rho*(clhs3*clhs82 + clhs81);
        lhs(1,4)=clhs31 + clhs33*clhs77 + clhs34 + clhs84;
        lhs(1,5)=-clhs13*clhs49 - clhs19*clhs85 + clhs3*clhs50 - clhs85 + clhs86;
        lhs(1,6)=rho*(clhs3*clhs89 + clhs88);
        lhs(1,7)=clhs55 + clhs56*clhs77 + clhs57 + clhs91;
        lhs(1,8)=-clhs13*clhs70 - clhs19*clhs92 + clhs3*clhs71 - clhs92 + clhs93;
        lhs(2,0)=clhs35*clhs95;
        lhs(2,1)=clhs15*clhs95;
        lhs(2,2)=clhs21 + clhs24*clhs96 + clhs26*clhs97;
        lhs(2,3)=rho*(DN(1,0)*clhs94 + clhs82);
        lhs(2,4)=rho*(DN(1,1)*clhs94 + clhs45);
        lhs(2,5)=clhs48*clhs96 + clhs49*clhs97 + clhs51;
        lhs(2,6)=rho*(DN(2,0)*clhs94 + clhs89);
        lhs(2,7)=rho*(DN(2,1)*clhs94 + clhs66);
        lhs(2,8)=clhs69*clhs96 + clhs70*clhs97 + clhs72;
        lhs(3,0)=clhs101 + clhs102*clhs9 + clhs36 + clhs38 + clhs40 + clhs41 + clhs99;
        lhs(3,1)=rho*(clhs1*clhs103 + clhs81);
        lhs(3,2)=clhs1*clhs105 - clhs100*clhs24 - clhs104*clhs19 - clhs104 + clhs53;
        lhs(3,3)=clhs102*clhs33 - clhs106*tau2 + clhs106 + clhs108;
        lhs(3,4)=clhs83*(DN(1,0) + clhs1*clhs44 - clhs109);
        lhs(3,5)=clhs1*clhs112 - clhs100*clhs48 - clhs110*clhs19 - clhs110 + clhs111*clhs22 - clhs113*clhs29;
        lhs(3,6)=clhs102*clhs56 + clhs114 + clhs115 + clhs116 + clhs118 + clhs119 + clhs120;
        lhs(3,7)=rho*(clhs1*clhs123 + clhs122);
        lhs(3,8)=clhs1*clhs125 - clhs100*clhs69 - clhs124*clhs19 - clhs124 + clhs128;
        lhs(4,0)=rho*(clhs129*clhs3 + clhs43);
        lhs(4,1)=clhs101 + clhs130*clhs9 + clhs84 + clhs99;
        lhs(4,2)=-clhs100*clhs26 + clhs105*clhs3 - clhs131*clhs19 - clhs131 + clhs86;
        lhs(4,3)=clhs39*(DN(1,1) - clhs132 + clhs3*clhs44);
        lhs(4,4)=clhs108 + clhs130*clhs33 - clhs133*tau2 + clhs133;
        lhs(4,5)=-clhs100*clhs49 + clhs111*clhs25 + clhs112*clhs3 - clhs113*clhs79 - clhs134*clhs19 - clhs134;
        lhs(4,6)=rho*(clhs136 + clhs137*clhs3);
        lhs(4,7)=clhs114 + clhs115 + clhs130*clhs56 + clhs138;
        lhs(4,8)=-clhs100*clhs70 + clhs125*clhs3 - clhs139*clhs19 - clhs139 + clhs140;
        lhs(5,0)=rho*(DN(0,0)*clhs141 + clhs129);
        lhs(5,1)=rho*(DN(0,1)*clhs141 + clhs103);
        lhs(5,2)=clhs142*clhs24 + clhs143*clhs26 + clhs51;
        lhs(5,3)=clhs144*clhs39;
        lhs(5,4)=clhs144*clhs83;
        lhs(5,5)=clhs111 + clhs142*clhs48 + clhs143*clhs49;
        lhs(5,6)=rho*(DN(2,0)*clhs141 + clhs137);
        lhs(5,7)=rho*(DN(2,1)*clhs141 + clhs123);
        lhs(5,8)=clhs126 + clhs142*clhs69 + clhs143*clhs70;
        lhs(6,0)=clhs146 + clhs148 + clhs149*clhs9 + clhs58 + clhs59 + clhs61 + clhs62;
        lhs(6,1)=rho*(clhs1*clhs150 + clhs88);
        lhs(6,2)=clhs1*clhs152 - clhs147*clhs24 - clhs151*clhs19 - clhs151 + clhs74;
        lhs(6,3)=clhs116 + clhs118 + clhs119 + clhs120 + clhs149*clhs33 + clhs153 + clhs154;
        lhs(6,4)=rho*(clhs1*clhs155 + clhs136);
        lhs(6,5)=clhs1*clhs157 + clhs128 - clhs147*clhs48 - clhs156*clhs19 - clhs156;
        lhs(6,6)=clhs149*clhs56 - clhs158*tau2 + clhs158 + clhs160;
        lhs(6,7)=clhs90*(-DN(2,0)*tau2 + DN(2,0) + clhs1*clhs65);
        lhs(6,8)=clhs1*clhs163 - clhs147*clhs69 - clhs161*clhs19 - clhs161 + clhs162*clhs22 - clhs164*clhs29;
        lhs(7,0)=rho*(clhs165*clhs3 + clhs64);
        lhs(7,1)=clhs146 + clhs148 + clhs166*clhs9 + clhs91;
        lhs(7,2)=-clhs147*clhs26 + clhs152*clhs3 - clhs167*clhs19 - clhs167 + clhs93;
        lhs(7,3)=rho*(clhs122 + clhs168*clhs3);
        lhs(7,4)=clhs138 + clhs153 + clhs154 + clhs166*clhs33;
        lhs(7,5)=clhs140 - clhs147*clhs49 + clhs157*clhs3 - clhs169*clhs19 - clhs169;
        lhs(7,6)=clhs60*(-DN(2,1)*tau2 + DN(2,1) + clhs3*clhs65);
        lhs(7,7)=clhs160 + clhs166*clhs56 - clhs170*tau2 + clhs170;
        lhs(7,8)=-clhs147*clhs70 + clhs162*clhs25 + clhs163*clhs3 - clhs164*clhs79 - clhs171*clhs19 - clhs171;
        lhs(8,0)=rho*(DN(0,0)*clhs172 + clhs165);
        lhs(8,1)=rho*(DN(0,1)*clhs172 + clhs150);
        lhs(8,2)=clhs173*clhs24 + clhs174*clhs26 + clhs72;
        lhs(8,3)=rho*(DN(1,0)*clhs172 + clhs168);
        lhs(8,4)=rho*(DN(1,1)*clhs172 + clhs155);
        lhs(8,5)=clhs126 + clhs173*clhs48 + clhs174*clhs49;
        lhs(8,6)=clhs175*clhs60;
        lhs(8,7)=clhs175*clhs90;
        lhs(8,8)=clhs162 + clhs173*clhs69 + clhs174*clhs70;

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
        const double crhs15 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0) + N[3]*v(3,0);
        const double crhs16 =             pow(c, -2);
        const double crhs17 =             N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]) + N[3]*(bdf0*p[3] + bdf1*pn[3] + bdf2*pnn[3]);
        const double crhs18 =             crhs16*crhs17;
        const double crhs19 =             N[0]*crhs18;
        const double crhs20 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
        const double crhs21 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
        const double crhs22 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
        const double crhs23 =             crhs20*(crhs3 + crhs5 + crhs7 + crhs9) + crhs21*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs22*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0));
        const double crhs24 =             tau2*(crhs12 + crhs18);
        const double crhs25 =             N[0]*crhs16*crhs17*tau1;
        const double crhs26 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + crhs15*crhs18 + rho*(crhs14 + crhs23);
        const double crhs27 =             DN(0,0)*crhs26;
        const double crhs28 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
        const double crhs29 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1) + N[3]*v(3,1);
        const double crhs30 =             N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1));
        const double crhs31 =             crhs20*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs21*(crhs10 + crhs4 + crhs6 + crhs8) + crhs22*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1));
        const double crhs32 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] + crhs18*crhs29 - crhs28 + rho*(crhs30 + crhs31);
        const double crhs33 =             DN(0,1)*crhs32;
        const double crhs34 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
        const double crhs35 =             N[0]*v(0,2) + N[1]*v(1,2) + N[2]*v(2,2) + N[3]*v(3,2);
        const double crhs36 =             N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2));
        const double crhs37 =             crhs2*crhs22 + crhs20*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs21*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2));
        const double crhs38 =             DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] + crhs18*crhs35 - crhs34 + rho*(crhs36 + crhs37);
        const double crhs39 =             DN(0,2)*crhs38;
        const double crhs40 =             rho*tau1*(crhs27 + crhs33 + crhs39);
        const double crhs41 =             2*crhs11*rho;
        const double crhs42 =             rho*tau1;
        const double crhs43 =             N[1]*rho;
        const double crhs44 =             N[1]*crhs18;
        const double crhs45 =             N[1]*crhs16*crhs17*tau1;
        const double crhs46 =             DN(1,0)*crhs26;
        const double crhs47 =             DN(1,1)*crhs32;
        const double crhs48 =             DN(1,2)*crhs38;
        const double crhs49 =             rho*tau1*(crhs46 + crhs47 + crhs48);
        const double crhs50 =             N[2]*rho;
        const double crhs51 =             N[2]*crhs18;
        const double crhs52 =             N[2]*crhs16*crhs17*tau1;
        const double crhs53 =             DN(2,0)*crhs26;
        const double crhs54 =             DN(2,1)*crhs32;
        const double crhs55 =             DN(2,2)*crhs38;
        const double crhs56 =             rho*tau1*(crhs53 + crhs54 + crhs55);
        const double crhs57 =             N[3]*rho;
        const double crhs58 =             N[3]*crhs18;
        const double crhs59 =             N[3]*crhs16*crhs17*tau1;
        const double crhs60 =             DN(3,0)*crhs26;
        const double crhs61 =             DN(3,1)*crhs32;
        const double crhs62 =             DN(3,2)*crhs38;
        const double crhs63 =             rho*tau1*(crhs60 + crhs61 + crhs62);

        rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs12 + DN(0,0)*crhs24 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - crhs13*crhs14 - crhs13*crhs23 - crhs15*crhs19 - crhs20*crhs40 + crhs25*crhs26;
        rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs12 + DN(0,1)*crhs24 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs28 - crhs13*crhs30 - crhs13*crhs31 - crhs19*crhs29 - crhs21*crhs40 + crhs25*crhs32;
        rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs12 + DN(0,2)*crhs24 - DN(0,2)*stress[2] + N[0]*crhs34 - crhs13*crhs36 - crhs13*crhs37 - crhs19*crhs35 - crhs22*crhs40 + crhs25*crhs38;
        rhs[3]=-N[0]*crhs41 - crhs19 - crhs27*crhs42 - crhs33*crhs42 - crhs39*crhs42;
        rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs12 + DN(1,0)*crhs24 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - crhs14*crhs43 - crhs15*crhs44 - crhs20*crhs49 - crhs23*crhs43 + crhs26*crhs45;
        rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs12 + DN(1,1)*crhs24 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs28 - crhs21*crhs49 - crhs29*crhs44 - crhs30*crhs43 - crhs31*crhs43 + crhs32*crhs45;
        rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs12 + DN(1,2)*crhs24 - DN(1,2)*stress[2] + N[1]*crhs34 - crhs22*crhs49 - crhs35*crhs44 - crhs36*crhs43 - crhs37*crhs43 + crhs38*crhs45;
        rhs[7]=-N[1]*crhs41 - crhs42*crhs46 - crhs42*crhs47 - crhs42*crhs48 - crhs44;
        rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs12 + DN(2,0)*crhs24 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - crhs14*crhs50 - crhs15*crhs51 - crhs20*crhs56 - crhs23*crhs50 + crhs26*crhs52;
        rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs12 + DN(2,1)*crhs24 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs28 - crhs21*crhs56 - crhs29*crhs51 - crhs30*crhs50 - crhs31*crhs50 + crhs32*crhs52;
        rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs12 + DN(2,2)*crhs24 - DN(2,2)*stress[2] + N[2]*crhs34 - crhs22*crhs56 - crhs35*crhs51 - crhs36*crhs50 - crhs37*crhs50 + crhs38*crhs52;
        rhs[11]=-N[2]*crhs41 - crhs42*crhs53 - crhs42*crhs54 - crhs42*crhs55 - crhs51;
        rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs12 + DN(3,0)*crhs24 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - crhs14*crhs57 - crhs15*crhs58 - crhs20*crhs63 - crhs23*crhs57 + crhs26*crhs59;
        rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs12 + DN(3,1)*crhs24 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs28 - crhs21*crhs63 - crhs29*crhs58 - crhs30*crhs57 - crhs31*crhs57 + crhs32*crhs59;
        rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs12 + DN(3,2)*crhs24 - DN(3,2)*stress[2] + N[3]*crhs34 - crhs22*crhs63 - crhs35*crhs58 - crhs36*crhs57 - crhs37*crhs57 + crhs38*crhs59;
        rhs[15]=-N[3]*crhs41 - crhs42*crhs60 - crhs42*crhs61 - crhs42*crhs62 - crhs58;

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
        const double crhs11 =             N[0]*v(0,0) + N[1]*v(1,0) + N[2]*v(2,0);
        const double crhs12 =             pow(c, -2);
        const double crhs13 =             N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]);
        const double crhs14 =             crhs12*crhs13;
        const double crhs15 =             N[0]*crhs14;
        const double crhs16 =             tau2*(crhs14 + crhs5);
        const double crhs17 =             N[0]*crhs12*crhs13*tau1;
        const double crhs18 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + crhs11*crhs14 + rho*(crhs10 + crhs7);
        const double crhs19 =             DN(0,0)*crhs18;
        const double crhs20 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
        const double crhs21 =             N[0]*v(0,1) + N[1]*v(1,1) + N[2]*v(2,1);
        const double crhs22 =             N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1));
        const double crhs23 =             crhs3*crhs9 + crhs8*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1));
        const double crhs24 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + crhs14*crhs21 - crhs20 + rho*(crhs22 + crhs23);
        const double crhs25 =             DN(0,1)*crhs24;
        const double crhs26 =             rho*tau1*(crhs19 + crhs25);
        const double crhs27 =             2*crhs4*rho;
        const double crhs28 =             rho*tau1;
        const double crhs29 =             N[1]*rho;
        const double crhs30 =             N[1]*crhs14;
        const double crhs31 =             N[1]*crhs12*crhs13*tau1;
        const double crhs32 =             DN(1,0)*crhs18;
        const double crhs33 =             DN(1,1)*crhs24;
        const double crhs34 =             rho*tau1*(crhs32 + crhs33);
        const double crhs35 =             N[2]*rho;
        const double crhs36 =             N[2]*crhs14;
        const double crhs37 =             N[2]*crhs12*crhs13*tau1;
        const double crhs38 =             DN(2,0)*crhs18;
        const double crhs39 =             DN(2,1)*crhs24;
        const double crhs40 =             rho*tau1*(crhs38 + crhs39);

        rhs[0]=DN(0,0)*crhs0 + DN(0,0)*crhs16 - DN(0,0)*crhs5 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - crhs10*crhs6 - crhs11*crhs15 + crhs17*crhs18 - crhs26*crhs8 - crhs6*crhs7;
        rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 + DN(0,1)*crhs16 - DN(0,1)*crhs5 - DN(0,1)*stress[1] + N[0]*crhs20 - crhs15*crhs21 + crhs17*crhs24 - crhs22*crhs6 - crhs23*crhs6 - crhs26*crhs9;
        rhs[2]=-N[0]*crhs27 - crhs15 - crhs19*crhs28 - crhs25*crhs28;
        rhs[3]=DN(1,0)*crhs0 + DN(1,0)*crhs16 - DN(1,0)*crhs5 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - crhs10*crhs29 - crhs11*crhs30 + crhs18*crhs31 - crhs29*crhs7 - crhs34*crhs8;
        rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 + DN(1,1)*crhs16 - DN(1,1)*crhs5 - DN(1,1)*stress[1] + N[1]*crhs20 - crhs21*crhs30 - crhs22*crhs29 - crhs23*crhs29 + crhs24*crhs31 - crhs34*crhs9;
        rhs[5]=-N[1]*crhs27 - crhs28*crhs32 - crhs28*crhs33 - crhs30;
        rhs[6]=DN(2,0)*crhs0 + DN(2,0)*crhs16 - DN(2,0)*crhs5 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - crhs10*crhs35 - crhs11*crhs36 + crhs18*crhs37 - crhs35*crhs7 - crhs40*crhs8;
        rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 + DN(2,1)*crhs16 - DN(2,1)*crhs5 - DN(2,1)*stress[1] + N[2]*crhs20 - crhs21*crhs36 - crhs22*crhs35 - crhs23*crhs35 + crhs24*crhs37 - crhs40*crhs9;
        rhs[8]=-N[2]*crhs27 - crhs28*crhs38 - crhs28*crhs39 - crhs36;

    }

}
