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
    // const double& bdf1 = data.bdf1;
    // const double& bdf2 = data.bdf2;
    const double& delta_t = data.delta_t;
    const double& dyn_tau_coeff = data.dyn_tau_coeff;

    const bounded_matrix<double,nnodes,dim>& v = data.v;
    // const bounded_matrix<double,nnodes,dim>& vn = data.vn;
    // const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
    const bounded_matrix<double,nnodes,dim>& vmesh = data.vmesh;
    const bounded_matrix<double,nnodes,dim>& vconv = v - vmesh;
    // const bounded_matrix<double,nnodes,dim>& f = data.f;
    // const array_1d<double,nnodes>& p = data.p;
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
    const double clhs1 =             C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
    const double clhs2 =             C(0,3)*DN(0,0);
    const double clhs3 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs2;
    const double clhs4 =             C(0,5)*DN(0,0);
    const double clhs5 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs4;
    const double clhs6 =             pow(N[0], 2);
    const double clhs7 =             bdf0*rho;
    const double clhs8 =             N[0]*rho;
    const double clhs9 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
    const double clhs10 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
    const double clhs11 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
    const double clhs12 =             DN(0,0)*clhs9 + DN(0,1)*clhs10 + DN(0,2)*clhs11;
    const double clhs13 =             pow(rho, 2);
    const double clhs14 =             N[0]*bdf0 + clhs12;
    const double clhs15 =             clhs13*clhs14*tau1;
    const double clhs16 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
    const double clhs17 =             N[0]*clhs16 + clhs12;
    const double clhs18 =             clhs12*clhs8 + clhs15*clhs17 + clhs6*clhs7;
    const double clhs19 =             DN(0,0)*tau2;
    const double clhs20 =             DN(0,1)*clhs19;
    const double clhs21 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs2;
    const double clhs22 =             C(1,3)*DN(0,1);
    const double clhs23 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs22;
    const double clhs24 =             C(3,5)*DN(0,0);
    const double clhs25 =             C(4,5)*DN(0,2);
    const double clhs26 =             C(1,5)*DN(0,1) + clhs24 + clhs25;
    const double clhs27 =             DN(0,2)*clhs19;
    const double clhs28 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs4;
    const double clhs29 =             C(3,4)*DN(0,1);
    const double clhs30 =             C(2,3)*DN(0,2) + clhs24 + clhs29;
    const double clhs31 =             C(2,5)*DN(0,2);
    const double clhs32 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs31;
    const double clhs33 =             pow(c, -2);
    const double clhs34 =             1.0/rho;
    const double clhs35 =             N[0]*bdf0*clhs33*clhs34;
    const double clhs36 =             rho*tau1;
    const double clhs37 =             clhs17*clhs36;
    const double clhs38 =             -N[0] + clhs35*tau2 + clhs37;
    const double clhs39 =             DN(1,0)*clhs19;
    const double clhs40 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
    const double clhs41 =             C(0,3)*DN(1,0);
    const double clhs42 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs41;
    const double clhs43 =             C(0,5)*DN(1,0);
    const double clhs44 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs43;
    const double clhs45 =             N[0]*bdf0*rho;
    const double clhs46 =             N[1]*clhs45;
    const double clhs47 =             DN(1,0)*clhs9 + DN(1,1)*clhs10 + DN(1,2)*clhs11;
    const double clhs48 =             N[1]*bdf0 + clhs47;
    const double clhs49 =             clhs13*clhs17*tau1;
    const double clhs50 =             clhs46 + clhs47*clhs8 + clhs48*clhs49;
    const double clhs51 =             DN(1,1)*clhs19;
    const double clhs52 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs41;
    const double clhs53 =             C(1,3)*DN(1,1);
    const double clhs54 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs53;
    const double clhs55 =             C(3,5)*DN(1,0);
    const double clhs56 =             C(4,5)*DN(1,2);
    const double clhs57 =             C(1,5)*DN(1,1) + clhs55 + clhs56;
    const double clhs58 =             DN(1,2)*clhs19;
    const double clhs59 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs43;
    const double clhs60 =             C(3,4)*DN(1,1);
    const double clhs61 =             C(2,3)*DN(1,2) + clhs55 + clhs60;
    const double clhs62 =             C(2,5)*DN(1,2);
    const double clhs63 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs62;
    const double clhs64 =             DN(0,0)*N[1];
    const double clhs65 =             bdf0*clhs33*clhs34*tau2;
    const double clhs66 =             DN(1,0)*rho*tau1;
    const double clhs67 =             DN(2,0)*clhs19;
    const double clhs68 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
    const double clhs69 =             C(0,3)*DN(2,0);
    const double clhs70 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs69;
    const double clhs71 =             C(0,5)*DN(2,0);
    const double clhs72 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs71;
    const double clhs73 =             N[2]*clhs45;
    const double clhs74 =             DN(2,0)*clhs9 + DN(2,1)*clhs10 + DN(2,2)*clhs11;
    const double clhs75 =             N[2]*bdf0;
    const double clhs76 =             clhs74 + clhs75;
    const double clhs77 =             clhs49*clhs76 + clhs73 + clhs74*clhs8;
    const double clhs78 =             DN(2,1)*clhs19;
    const double clhs79 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs69;
    const double clhs80 =             C(1,3)*DN(2,1);
    const double clhs81 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs80;
    const double clhs82 =             C(3,5)*DN(2,0);
    const double clhs83 =             C(4,5)*DN(2,2);
    const double clhs84 =             C(1,5)*DN(2,1) + clhs82 + clhs83;
    const double clhs85 =             DN(2,2)*clhs19;
    const double clhs86 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs71;
    const double clhs87 =             C(3,4)*DN(2,1);
    const double clhs88 =             C(2,3)*DN(2,2) + clhs82 + clhs87;
    const double clhs89 =             C(2,5)*DN(2,2);
    const double clhs90 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs89;
    const double clhs91 =             DN(0,0)*N[2];
    const double clhs92 =             DN(2,0)*rho*tau1;
    const double clhs93 =             DN(3,0)*clhs19;
    const double clhs94 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
    const double clhs95 =             C(0,3)*DN(3,0);
    const double clhs96 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs95;
    const double clhs97 =             C(0,5)*DN(3,0);
    const double clhs98 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs97;
    const double clhs99 =             N[3]*clhs45;
    const double clhs100 =             DN(3,0)*clhs9 + DN(3,1)*clhs10 + DN(3,2)*clhs11;
    const double clhs101 =             N[3]*bdf0;
    const double clhs102 =             clhs100 + clhs101;
    const double clhs103 =             clhs100*clhs8 + clhs102*clhs49 + clhs99;
    const double clhs104 =             DN(3,1)*clhs19;
    const double clhs105 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs95;
    const double clhs106 =             C(1,3)*DN(3,1);
    const double clhs107 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs106;
    const double clhs108 =             C(3,5)*DN(3,0);
    const double clhs109 =             C(4,5)*DN(3,2);
    const double clhs110 =             C(1,5)*DN(3,1) + clhs108 + clhs109;
    const double clhs111 =             DN(3,2)*clhs19;
    const double clhs112 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs97;
    const double clhs113 =             C(3,4)*DN(3,1);
    const double clhs114 =             C(2,3)*DN(3,2) + clhs108 + clhs113;
    const double clhs115 =             C(2,5)*DN(3,2);
    const double clhs116 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs115;
    const double clhs117 =             DN(0,0)*N[3];
    const double clhs118 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs22;
    const double clhs119 =             C(0,4)*DN(0,0) + clhs25 + clhs29;
    const double clhs120 =             pow(DN(0,1), 2);
    const double clhs121 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
    const double clhs122 =             C(1,4)*DN(0,1);
    const double clhs123 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs122;
    const double clhs124 =             DN(0,1)*tau2;
    const double clhs125 =             DN(0,2)*clhs124;
    const double clhs126 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs122;
    const double clhs127 =             C(2,4)*DN(0,2);
    const double clhs128 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs127;
    const double clhs129 =             DN(1,0)*clhs124;
    const double clhs130 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs53;
    const double clhs131 =             C(0,4)*DN(1,0) + clhs56 + clhs60;
    const double clhs132 =             DN(1,1)*clhs124;
    const double clhs133 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
    const double clhs134 =             C(1,4)*DN(1,1);
    const double clhs135 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs134;
    const double clhs136 =             DN(1,2)*clhs124;
    const double clhs137 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs134;
    const double clhs138 =             C(2,4)*DN(1,2);
    const double clhs139 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs138;
    const double clhs140 =             DN(0,1)*N[1];
    const double clhs141 =             DN(1,1)*rho*tau1;
    const double clhs142 =             DN(2,0)*clhs124;
    const double clhs143 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs80;
    const double clhs144 =             C(0,4)*DN(2,0) + clhs83 + clhs87;
    const double clhs145 =             DN(2,1)*clhs124;
    const double clhs146 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
    const double clhs147 =             C(1,4)*DN(2,1);
    const double clhs148 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs147;
    const double clhs149 =             DN(2,2)*clhs124;
    const double clhs150 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs147;
    const double clhs151 =             C(2,4)*DN(2,2);
    const double clhs152 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs151;
    const double clhs153 =             DN(0,1)*N[2];
    const double clhs154 =             DN(2,1)*rho*tau1;
    const double clhs155 =             DN(3,0)*clhs124;
    const double clhs156 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs106;
    const double clhs157 =             C(0,4)*DN(3,0) + clhs109 + clhs113;
    const double clhs158 =             DN(3,1)*clhs124;
    const double clhs159 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
    const double clhs160 =             C(1,4)*DN(3,1);
    const double clhs161 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs160;
    const double clhs162 =             DN(3,2)*clhs124;
    const double clhs163 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs160;
    const double clhs164 =             C(2,4)*DN(3,2);
    const double clhs165 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs164;
    const double clhs166 =             DN(0,1)*N[3];
    const double clhs167 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs31;
    const double clhs168 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs127;
    const double clhs169 =             pow(DN(0,2), 2);
    const double clhs170 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
    const double clhs171 =             DN(0,2)*tau2;
    const double clhs172 =             DN(1,0)*clhs171;
    const double clhs173 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs62;
    const double clhs174 =             DN(1,1)*clhs171;
    const double clhs175 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs138;
    const double clhs176 =             DN(1,2)*clhs171;
    const double clhs177 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
    const double clhs178 =             DN(0,2)*N[1];
    const double clhs179 =             DN(1,2)*rho*tau1;
    const double clhs180 =             DN(2,0)*clhs171;
    const double clhs181 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs89;
    const double clhs182 =             DN(2,1)*clhs171;
    const double clhs183 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs151;
    const double clhs184 =             DN(2,2)*clhs171;
    const double clhs185 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
    const double clhs186 =             DN(0,2)*N[2];
    const double clhs187 =             DN(2,2)*rho*tau1;
    const double clhs188 =             DN(3,0)*clhs171;
    const double clhs189 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs115;
    const double clhs190 =             DN(3,1)*clhs171;
    const double clhs191 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs164;
    const double clhs192 =             DN(3,2)*clhs171;
    const double clhs193 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
    const double clhs194 =             DN(0,2)*N[3];
    const double clhs195 =             clhs14*clhs36;
    const double clhs196 =             N[0] + clhs195;
    const double clhs197 =             bdf0*clhs33*clhs34;
    const double clhs198 =             DN(1,0)*N[0];
    const double clhs199 =             DN(0,0)*rho*tau1;
    const double clhs200 =             DN(1,1)*N[0];
    const double clhs201 =             DN(0,1)*rho*tau1;
    const double clhs202 =             DN(1,2)*N[0];
    const double clhs203 =             DN(0,2)*rho*tau1;
    const double clhs204 =             DN(0,0)*tau1;
    const double clhs205 =             DN(0,1)*tau1;
    const double clhs206 =             DN(0,2)*tau1;
    const double clhs207 =             DN(1,0)*clhs204 + DN(1,1)*clhs205 + DN(1,2)*clhs206 + N[1]*clhs35;
    const double clhs208 =             DN(2,0)*N[0];
    const double clhs209 =             DN(2,1)*N[0];
    const double clhs210 =             DN(2,2)*N[0];
    const double clhs211 =             DN(2,0)*clhs204 + DN(2,1)*clhs205 + DN(2,2)*clhs206 + N[2]*clhs35;
    const double clhs212 =             DN(3,0)*N[0];
    const double clhs213 =             DN(3,1)*N[0];
    const double clhs214 =             DN(3,2)*N[0];
    const double clhs215 =             DN(3,0)*clhs204 + DN(3,1)*clhs205 + DN(3,2)*clhs206 + N[3]*clhs35;
    const double clhs216 =             N[1]*rho;
    const double clhs217 =             N[1]*clhs16 + clhs47;
    const double clhs218 =             clhs12*clhs216 + clhs15*clhs217 + clhs46;
    const double clhs219 =             pow(DN(1,0), 2);
    const double clhs220 =             pow(N[1], 2);
    const double clhs221 =             clhs13*clhs48*tau1;
    const double clhs222 =             clhs216*clhs47 + clhs217*clhs221 + clhs220*clhs7;
    const double clhs223 =             DN(1,0)*tau2;
    const double clhs224 =             DN(1,1)*clhs223;
    const double clhs225 =             DN(1,2)*clhs223;
    const double clhs226 =             N[1]*bdf0*clhs33*clhs34;
    const double clhs227 =             clhs217*clhs36;
    const double clhs228 =             -N[1] + clhs226*tau2 + clhs227;
    const double clhs229 =             DN(2,0)*clhs223;
    const double clhs230 =             N[1]*bdf0*rho;
    const double clhs231 =             N[2]*clhs230;
    const double clhs232 =             clhs13*clhs217*tau1;
    const double clhs233 =             clhs216*clhs74 + clhs231 + clhs232*clhs76;
    const double clhs234 =             DN(2,1)*clhs223;
    const double clhs235 =             DN(2,2)*clhs223;
    const double clhs236 =             DN(1,0)*N[2];
    const double clhs237 =             DN(3,0)*clhs223;
    const double clhs238 =             N[3]*clhs230;
    const double clhs239 =             clhs100*clhs216 + clhs102*clhs232 + clhs238;
    const double clhs240 =             DN(3,1)*clhs223;
    const double clhs241 =             DN(3,2)*clhs223;
    const double clhs242 =             DN(1,0)*N[3];
    const double clhs243 =             pow(DN(1,1), 2);
    const double clhs244 =             DN(1,1)*tau2;
    const double clhs245 =             DN(1,2)*clhs244;
    const double clhs246 =             DN(2,0)*clhs244;
    const double clhs247 =             DN(2,1)*clhs244;
    const double clhs248 =             DN(2,2)*clhs244;
    const double clhs249 =             DN(1,1)*N[2];
    const double clhs250 =             DN(3,0)*clhs244;
    const double clhs251 =             DN(3,1)*clhs244;
    const double clhs252 =             DN(3,2)*clhs244;
    const double clhs253 =             DN(1,1)*N[3];
    const double clhs254 =             pow(DN(1,2), 2);
    const double clhs255 =             DN(1,2)*tau2;
    const double clhs256 =             DN(2,0)*clhs255;
    const double clhs257 =             DN(2,1)*clhs255;
    const double clhs258 =             DN(2,2)*clhs255;
    const double clhs259 =             DN(1,2)*N[2];
    const double clhs260 =             DN(3,0)*clhs255;
    const double clhs261 =             DN(3,1)*clhs255;
    const double clhs262 =             DN(3,2)*clhs255;
    const double clhs263 =             DN(1,2)*N[3];
    const double clhs264 =             clhs36*clhs48;
    const double clhs265 =             N[1] + clhs264;
    const double clhs266 =             DN(2,0)*N[1];
    const double clhs267 =             DN(2,1)*N[1];
    const double clhs268 =             DN(2,2)*N[1];
    const double clhs269 =             DN(1,0)*tau1;
    const double clhs270 =             DN(1,1)*tau1;
    const double clhs271 =             DN(1,2)*tau1;
    const double clhs272 =             DN(2,0)*clhs269 + DN(2,1)*clhs270 + DN(2,2)*clhs271 + N[2]*clhs226;
    const double clhs273 =             DN(3,0)*N[1];
    const double clhs274 =             DN(3,1)*N[1];
    const double clhs275 =             DN(3,2)*N[1];
    const double clhs276 =             DN(3,0)*clhs269 + DN(3,1)*clhs270 + DN(3,2)*clhs271 + N[3]*clhs226;
    const double clhs277 =             N[2]*rho;
    const double clhs278 =             N[2]*clhs16 + clhs74;
    const double clhs279 =             clhs12*clhs277 + clhs15*clhs278 + clhs73;
    const double clhs280 =             clhs221*clhs278 + clhs231 + clhs277*clhs47;
    const double clhs281 =             pow(DN(2,0), 2);
    const double clhs282 =             pow(N[2], 2);
    const double clhs283 =             clhs13*clhs76*tau1;
    const double clhs284 =             clhs277*clhs74 + clhs278*clhs283 + clhs282*clhs7;
    const double clhs285 =             DN(2,0)*tau2;
    const double clhs286 =             DN(2,1)*clhs285;
    const double clhs287 =             DN(2,2)*clhs285;
    const double clhs288 =             clhs33*clhs34*tau2;
    const double clhs289 =             clhs278*clhs36;
    const double clhs290 =             -N[2] + clhs288*clhs75 + clhs289;
    const double clhs291 =             DN(3,0)*clhs285;
    const double clhs292 =             N[2]*N[3]*bdf0;
    const double clhs293 =             clhs292*rho;
    const double clhs294 =             clhs102*clhs13*tau1;
    const double clhs295 =             clhs100*clhs277 + clhs278*clhs294 + clhs293;
    const double clhs296 =             DN(3,1)*clhs285;
    const double clhs297 =             DN(3,2)*clhs285;
    const double clhs298 =             DN(2,0)*N[3];
    const double clhs299 =             pow(DN(2,1), 2);
    const double clhs300 =             DN(2,1)*tau2;
    const double clhs301 =             DN(2,2)*clhs300;
    const double clhs302 =             DN(3,0)*clhs300;
    const double clhs303 =             DN(3,1)*clhs300;
    const double clhs304 =             DN(3,2)*clhs300;
    const double clhs305 =             DN(2,1)*N[3];
    const double clhs306 =             pow(DN(2,2), 2);
    const double clhs307 =             DN(2,2)*tau2;
    const double clhs308 =             DN(3,0)*clhs307;
    const double clhs309 =             DN(3,1)*clhs307;
    const double clhs310 =             DN(3,2)*clhs307;
    const double clhs311 =             DN(2,2)*N[3];
    const double clhs312 =             clhs36*clhs76;
    const double clhs313 =             N[2] + clhs312;
    const double clhs314 =             DN(3,0)*N[2];
    const double clhs315 =             DN(3,1)*N[2];
    const double clhs316 =             DN(3,2)*N[2];
    const double clhs317 =             DN(3,0)*(DN(2,0)*tau1) + DN(3,1)*(DN(2,1)*tau1) + DN(3,2)*(DN(2,2)*tau1) + clhs292*clhs33*clhs34;
    const double clhs318 =             N[3]*rho;
    const double clhs319 =             N[3]*clhs16 + clhs100;
    const double clhs320 =             clhs12*clhs318 + clhs15*clhs319 + clhs99;
    const double clhs321 =             clhs221*clhs319 + clhs238 + clhs318*clhs47;
    const double clhs322 =             clhs283*clhs319 + clhs293 + clhs318*clhs74;
    const double clhs323 =             pow(DN(3,0), 2);
    const double clhs324 =             pow(N[3], 2);
    const double clhs325 =             clhs100*clhs318 + clhs294*clhs319 + clhs324*clhs7;
    const double clhs326 =             DN(3,0)*tau2;
    const double clhs327 =             DN(3,1)*clhs326;
    const double clhs328 =             DN(3,2)*clhs326;
    const double clhs329 =             -N[3] + clhs101*clhs288 + clhs319*clhs36;
    const double clhs330 =             pow(DN(3,1), 2);
    const double clhs331 =             DN(3,1)*DN(3,2)*tau2;
    const double clhs332 =             pow(DN(3,2), 2);
    const double clhs333 =             N[3] + clhs102*clhs36;

    lhs(0,0)=DN(0,0)*clhs1 + DN(0,1)*clhs3 + DN(0,2)*clhs5 + clhs0*tau2 + clhs18;
    lhs(0,1)=DN(0,0)*clhs21 + DN(0,1)*clhs23 + DN(0,2)*clhs26 + clhs20;
    lhs(0,2)=DN(0,0)*clhs28 + DN(0,1)*clhs30 + DN(0,2)*clhs32 + clhs27;
    lhs(0,3)=DN(0,0)*clhs38;
    lhs(0,4)=DN(0,0)*clhs40 + DN(0,1)*clhs42 + DN(0,2)*clhs44 + clhs39 + clhs50;
    lhs(0,5)=DN(0,0)*clhs52 + DN(0,1)*clhs54 + DN(0,2)*clhs57 + clhs51;
    lhs(0,6)=DN(0,0)*clhs59 + DN(0,1)*clhs61 + DN(0,2)*clhs63 + clhs58;
    lhs(0,7)=clhs17*clhs66 + clhs64*clhs65 - clhs64;
    lhs(0,8)=DN(0,0)*clhs68 + DN(0,1)*clhs70 + DN(0,2)*clhs72 + clhs67 + clhs77;
    lhs(0,9)=DN(0,0)*clhs79 + DN(0,1)*clhs81 + DN(0,2)*clhs84 + clhs78;
    lhs(0,10)=DN(0,0)*clhs86 + DN(0,1)*clhs88 + DN(0,2)*clhs90 + clhs85;
    lhs(0,11)=clhs17*clhs92 + clhs65*clhs91 - clhs91;
    lhs(0,12)=DN(0,0)*clhs94 + DN(0,1)*clhs96 + DN(0,2)*clhs98 + clhs103 + clhs93;
    lhs(0,13)=DN(0,0)*clhs105 + DN(0,1)*clhs107 + DN(0,2)*clhs110 + clhs104;
    lhs(0,14)=DN(0,0)*clhs112 + DN(0,1)*clhs114 + DN(0,2)*clhs116 + clhs111;
    lhs(0,15)=DN(3,0)*clhs37 + clhs117*clhs65 - clhs117;
    lhs(1,0)=DN(0,0)*clhs3 + DN(0,1)*clhs118 + DN(0,2)*clhs119 + clhs20;
    lhs(1,1)=DN(0,0)*clhs23 + DN(0,1)*clhs121 + DN(0,2)*clhs123 + clhs120*tau2 + clhs18;
    lhs(1,2)=DN(0,0)*clhs30 + DN(0,1)*clhs126 + DN(0,2)*clhs128 + clhs125;
    lhs(1,3)=DN(0,1)*clhs38;
    lhs(1,4)=DN(0,0)*clhs42 + DN(0,1)*clhs130 + DN(0,2)*clhs131 + clhs129;
    lhs(1,5)=DN(0,0)*clhs54 + DN(0,1)*clhs133 + DN(0,2)*clhs135 + clhs132 + clhs50;
    lhs(1,6)=DN(0,0)*clhs61 + DN(0,1)*clhs137 + DN(0,2)*clhs139 + clhs136;
    lhs(1,7)=clhs140*clhs65 - clhs140 + clhs141*clhs17;
    lhs(1,8)=DN(0,0)*clhs70 + DN(0,1)*clhs143 + DN(0,2)*clhs144 + clhs142;
    lhs(1,9)=DN(0,0)*clhs81 + DN(0,1)*clhs146 + DN(0,2)*clhs148 + clhs145 + clhs77;
    lhs(1,10)=DN(0,0)*clhs88 + DN(0,1)*clhs150 + DN(0,2)*clhs152 + clhs149;
    lhs(1,11)=clhs153*clhs65 - clhs153 + clhs154*clhs17;
    lhs(1,12)=DN(0,0)*clhs96 + DN(0,1)*clhs156 + DN(0,2)*clhs157 + clhs155;
    lhs(1,13)=DN(0,0)*clhs107 + DN(0,1)*clhs159 + DN(0,2)*clhs161 + clhs103 + clhs158;
    lhs(1,14)=DN(0,0)*clhs114 + DN(0,1)*clhs163 + DN(0,2)*clhs165 + clhs162;
    lhs(1,15)=DN(3,1)*clhs37 + clhs166*clhs65 - clhs166;
    lhs(2,0)=DN(0,0)*clhs5 + DN(0,1)*clhs119 + DN(0,2)*clhs167 + clhs27;
    lhs(2,1)=DN(0,0)*clhs26 + DN(0,1)*clhs123 + DN(0,2)*clhs168 + clhs125;
    lhs(2,2)=DN(0,0)*clhs32 + DN(0,1)*clhs128 + DN(0,2)*clhs170 + clhs169*tau2 + clhs18;
    lhs(2,3)=DN(0,2)*clhs38;
    lhs(2,4)=DN(0,0)*clhs44 + DN(0,1)*clhs131 + DN(0,2)*clhs173 + clhs172;
    lhs(2,5)=DN(0,0)*clhs57 + DN(0,1)*clhs135 + DN(0,2)*clhs175 + clhs174;
    lhs(2,6)=DN(0,0)*clhs63 + DN(0,1)*clhs139 + DN(0,2)*clhs177 + clhs176 + clhs50;
    lhs(2,7)=clhs17*clhs179 + clhs178*clhs65 - clhs178;
    lhs(2,8)=DN(0,0)*clhs72 + DN(0,1)*clhs144 + DN(0,2)*clhs181 + clhs180;
    lhs(2,9)=DN(0,0)*clhs84 + DN(0,1)*clhs148 + DN(0,2)*clhs183 + clhs182;
    lhs(2,10)=DN(0,0)*clhs90 + DN(0,1)*clhs152 + DN(0,2)*clhs185 + clhs184 + clhs77;
    lhs(2,11)=clhs17*clhs187 + clhs186*clhs65 - clhs186;
    lhs(2,12)=DN(0,0)*clhs98 + DN(0,1)*clhs157 + DN(0,2)*clhs189 + clhs188;
    lhs(2,13)=DN(0,0)*clhs110 + DN(0,1)*clhs161 + DN(0,2)*clhs191 + clhs190;
    lhs(2,14)=DN(0,0)*clhs116 + DN(0,1)*clhs165 + DN(0,2)*clhs193 + clhs103 + clhs192;
    lhs(2,15)=DN(3,2)*clhs37 + clhs194*clhs65 - clhs194;
    lhs(3,0)=DN(0,0)*clhs196;
    lhs(3,1)=DN(0,1)*clhs196;
    lhs(3,2)=DN(0,2)*clhs196;
    lhs(3,3)=clhs0*tau1 + clhs120*tau1 + clhs169*tau1 + clhs197*clhs6;
    lhs(3,4)=clhs198 + clhs199*clhs48;
    lhs(3,5)=clhs200 + clhs201*clhs48;
    lhs(3,6)=clhs202 + clhs203*clhs48;
    lhs(3,7)=clhs207;
    lhs(3,8)=clhs199*clhs76 + clhs208;
    lhs(3,9)=clhs201*clhs76 + clhs209;
    lhs(3,10)=clhs203*clhs76 + clhs210;
    lhs(3,11)=clhs211;
    lhs(3,12)=clhs102*clhs199 + clhs212;
    lhs(3,13)=clhs102*clhs201 + clhs213;
    lhs(3,14)=clhs102*clhs203 + clhs214;
    lhs(3,15)=clhs215;
    lhs(4,0)=DN(1,0)*clhs1 + DN(1,1)*clhs3 + DN(1,2)*clhs5 + clhs218 + clhs39;
    lhs(4,1)=DN(1,0)*clhs21 + DN(1,1)*clhs23 + DN(1,2)*clhs26 + clhs129;
    lhs(4,2)=DN(1,0)*clhs28 + DN(1,1)*clhs30 + DN(1,2)*clhs32 + clhs172;
    lhs(4,3)=clhs198*clhs65 - clhs198 + clhs199*clhs217;
    lhs(4,4)=DN(1,0)*clhs40 + DN(1,1)*clhs42 + DN(1,2)*clhs44 + clhs219*tau2 + clhs222;
    lhs(4,5)=DN(1,0)*clhs52 + DN(1,1)*clhs54 + DN(1,2)*clhs57 + clhs224;
    lhs(4,6)=DN(1,0)*clhs59 + DN(1,1)*clhs61 + DN(1,2)*clhs63 + clhs225;
    lhs(4,7)=DN(1,0)*clhs228;
    lhs(4,8)=DN(1,0)*clhs68 + DN(1,1)*clhs70 + DN(1,2)*clhs72 + clhs229 + clhs233;
    lhs(4,9)=DN(1,0)*clhs79 + DN(1,1)*clhs81 + DN(1,2)*clhs84 + clhs234;
    lhs(4,10)=DN(1,0)*clhs86 + DN(1,1)*clhs88 + DN(1,2)*clhs90 + clhs235;
    lhs(4,11)=clhs217*clhs92 + clhs236*clhs65 - clhs236;
    lhs(4,12)=DN(1,0)*clhs94 + DN(1,1)*clhs96 + DN(1,2)*clhs98 + clhs237 + clhs239;
    lhs(4,13)=DN(1,0)*clhs105 + DN(1,1)*clhs107 + DN(1,2)*clhs110 + clhs240;
    lhs(4,14)=DN(1,0)*clhs112 + DN(1,1)*clhs114 + DN(1,2)*clhs116 + clhs241;
    lhs(4,15)=DN(3,0)*clhs227 + clhs242*clhs65 - clhs242;
    lhs(5,0)=DN(1,0)*clhs3 + DN(1,1)*clhs118 + DN(1,2)*clhs119 + clhs51;
    lhs(5,1)=DN(1,0)*clhs23 + DN(1,1)*clhs121 + DN(1,2)*clhs123 + clhs132 + clhs218;
    lhs(5,2)=DN(1,0)*clhs30 + DN(1,1)*clhs126 + DN(1,2)*clhs128 + clhs174;
    lhs(5,3)=clhs200*clhs65 - clhs200 + clhs201*clhs217;
    lhs(5,4)=DN(1,0)*clhs42 + DN(1,1)*clhs130 + DN(1,2)*clhs131 + clhs224;
    lhs(5,5)=DN(1,0)*clhs54 + DN(1,1)*clhs133 + DN(1,2)*clhs135 + clhs222 + clhs243*tau2;
    lhs(5,6)=DN(1,0)*clhs61 + DN(1,1)*clhs137 + DN(1,2)*clhs139 + clhs245;
    lhs(5,7)=DN(1,1)*clhs228;
    lhs(5,8)=DN(1,0)*clhs70 + DN(1,1)*clhs143 + DN(1,2)*clhs144 + clhs246;
    lhs(5,9)=DN(1,0)*clhs81 + DN(1,1)*clhs146 + DN(1,2)*clhs148 + clhs233 + clhs247;
    lhs(5,10)=DN(1,0)*clhs88 + DN(1,1)*clhs150 + DN(1,2)*clhs152 + clhs248;
    lhs(5,11)=clhs154*clhs217 + clhs249*clhs65 - clhs249;
    lhs(5,12)=DN(1,0)*clhs96 + DN(1,1)*clhs156 + DN(1,2)*clhs157 + clhs250;
    lhs(5,13)=DN(1,0)*clhs107 + DN(1,1)*clhs159 + DN(1,2)*clhs161 + clhs239 + clhs251;
    lhs(5,14)=DN(1,0)*clhs114 + DN(1,1)*clhs163 + DN(1,2)*clhs165 + clhs252;
    lhs(5,15)=DN(3,1)*clhs227 + clhs253*clhs65 - clhs253;
    lhs(6,0)=DN(1,0)*clhs5 + DN(1,1)*clhs119 + DN(1,2)*clhs167 + clhs58;
    lhs(6,1)=DN(1,0)*clhs26 + DN(1,1)*clhs123 + DN(1,2)*clhs168 + clhs136;
    lhs(6,2)=DN(1,0)*clhs32 + DN(1,1)*clhs128 + DN(1,2)*clhs170 + clhs176 + clhs218;
    lhs(6,3)=clhs202*clhs65 - clhs202 + clhs203*clhs217;
    lhs(6,4)=DN(1,0)*clhs44 + DN(1,1)*clhs131 + DN(1,2)*clhs173 + clhs225;
    lhs(6,5)=DN(1,0)*clhs57 + DN(1,1)*clhs135 + DN(1,2)*clhs175 + clhs245;
    lhs(6,6)=DN(1,0)*clhs63 + DN(1,1)*clhs139 + DN(1,2)*clhs177 + clhs222 + clhs254*tau2;
    lhs(6,7)=DN(1,2)*clhs228;
    lhs(6,8)=DN(1,0)*clhs72 + DN(1,1)*clhs144 + DN(1,2)*clhs181 + clhs256;
    lhs(6,9)=DN(1,0)*clhs84 + DN(1,1)*clhs148 + DN(1,2)*clhs183 + clhs257;
    lhs(6,10)=DN(1,0)*clhs90 + DN(1,1)*clhs152 + DN(1,2)*clhs185 + clhs233 + clhs258;
    lhs(6,11)=clhs187*clhs217 + clhs259*clhs65 - clhs259;
    lhs(6,12)=DN(1,0)*clhs98 + DN(1,1)*clhs157 + DN(1,2)*clhs189 + clhs260;
    lhs(6,13)=DN(1,0)*clhs110 + DN(1,1)*clhs161 + DN(1,2)*clhs191 + clhs261;
    lhs(6,14)=DN(1,0)*clhs116 + DN(1,1)*clhs165 + DN(1,2)*clhs193 + clhs239 + clhs262;
    lhs(6,15)=DN(3,2)*clhs227 + clhs263*clhs65 - clhs263;
    lhs(7,0)=clhs14*clhs66 + clhs64;
    lhs(7,1)=clhs14*clhs141 + clhs140;
    lhs(7,2)=clhs14*clhs179 + clhs178;
    lhs(7,3)=clhs207;
    lhs(7,4)=DN(1,0)*clhs265;
    lhs(7,5)=DN(1,1)*clhs265;
    lhs(7,6)=DN(1,2)*clhs265;
    lhs(7,7)=clhs197*clhs220 + clhs219*tau1 + clhs243*tau1 + clhs254*tau1;
    lhs(7,8)=clhs266 + clhs66*clhs76;
    lhs(7,9)=clhs141*clhs76 + clhs267;
    lhs(7,10)=clhs179*clhs76 + clhs268;
    lhs(7,11)=clhs272;
    lhs(7,12)=clhs102*clhs66 + clhs273;
    lhs(7,13)=clhs102*clhs141 + clhs274;
    lhs(7,14)=clhs102*clhs179 + clhs275;
    lhs(7,15)=clhs276;
    lhs(8,0)=DN(2,0)*clhs1 + DN(2,1)*clhs3 + DN(2,2)*clhs5 + clhs279 + clhs67;
    lhs(8,1)=DN(2,0)*clhs21 + DN(2,1)*clhs23 + DN(2,2)*clhs26 + clhs142;
    lhs(8,2)=DN(2,0)*clhs28 + DN(2,1)*clhs30 + DN(2,2)*clhs32 + clhs180;
    lhs(8,3)=clhs199*clhs278 + clhs208*clhs65 - clhs208;
    lhs(8,4)=DN(2,0)*clhs40 + DN(2,1)*clhs42 + DN(2,2)*clhs44 + clhs229 + clhs280;
    lhs(8,5)=DN(2,0)*clhs52 + DN(2,1)*clhs54 + DN(2,2)*clhs57 + clhs246;
    lhs(8,6)=DN(2,0)*clhs59 + DN(2,1)*clhs61 + DN(2,2)*clhs63 + clhs256;
    lhs(8,7)=clhs266*clhs65 - clhs266 + clhs278*clhs66;
    lhs(8,8)=DN(2,0)*clhs68 + DN(2,1)*clhs70 + DN(2,2)*clhs72 + clhs281*tau2 + clhs284;
    lhs(8,9)=DN(2,0)*clhs79 + DN(2,1)*clhs81 + DN(2,2)*clhs84 + clhs286;
    lhs(8,10)=DN(2,0)*clhs86 + DN(2,1)*clhs88 + DN(2,2)*clhs90 + clhs287;
    lhs(8,11)=DN(2,0)*clhs290;
    lhs(8,12)=DN(2,0)*clhs94 + DN(2,1)*clhs96 + DN(2,2)*clhs98 + clhs291 + clhs295;
    lhs(8,13)=DN(2,0)*clhs105 + DN(2,1)*clhs107 + DN(2,2)*clhs110 + clhs296;
    lhs(8,14)=DN(2,0)*clhs112 + DN(2,1)*clhs114 + DN(2,2)*clhs116 + clhs297;
    lhs(8,15)=DN(3,0)*clhs289 + clhs298*clhs65 - clhs298;
    lhs(9,0)=DN(2,0)*clhs3 + DN(2,1)*clhs118 + DN(2,2)*clhs119 + clhs78;
    lhs(9,1)=DN(2,0)*clhs23 + DN(2,1)*clhs121 + DN(2,2)*clhs123 + clhs145 + clhs279;
    lhs(9,2)=DN(2,0)*clhs30 + DN(2,1)*clhs126 + DN(2,2)*clhs128 + clhs182;
    lhs(9,3)=clhs201*clhs278 + clhs209*clhs65 - clhs209;
    lhs(9,4)=DN(2,0)*clhs42 + DN(2,1)*clhs130 + DN(2,2)*clhs131 + clhs234;
    lhs(9,5)=DN(2,0)*clhs54 + DN(2,1)*clhs133 + DN(2,2)*clhs135 + clhs247 + clhs280;
    lhs(9,6)=DN(2,0)*clhs61 + DN(2,1)*clhs137 + DN(2,2)*clhs139 + clhs257;
    lhs(9,7)=clhs141*clhs278 + clhs267*clhs65 - clhs267;
    lhs(9,8)=DN(2,0)*clhs70 + DN(2,1)*clhs143 + DN(2,2)*clhs144 + clhs286;
    lhs(9,9)=DN(2,0)*clhs81 + DN(2,1)*clhs146 + DN(2,2)*clhs148 + clhs284 + clhs299*tau2;
    lhs(9,10)=DN(2,0)*clhs88 + DN(2,1)*clhs150 + DN(2,2)*clhs152 + clhs301;
    lhs(9,11)=DN(2,1)*clhs290;
    lhs(9,12)=DN(2,0)*clhs96 + DN(2,1)*clhs156 + DN(2,2)*clhs157 + clhs302;
    lhs(9,13)=DN(2,0)*clhs107 + DN(2,1)*clhs159 + DN(2,2)*clhs161 + clhs295 + clhs303;
    lhs(9,14)=DN(2,0)*clhs114 + DN(2,1)*clhs163 + DN(2,2)*clhs165 + clhs304;
    lhs(9,15)=DN(3,1)*clhs289 + clhs305*clhs65 - clhs305;
    lhs(10,0)=DN(2,0)*clhs5 + DN(2,1)*clhs119 + DN(2,2)*clhs167 + clhs85;
    lhs(10,1)=DN(2,0)*clhs26 + DN(2,1)*clhs123 + DN(2,2)*clhs168 + clhs149;
    lhs(10,2)=DN(2,0)*clhs32 + DN(2,1)*clhs128 + DN(2,2)*clhs170 + clhs184 + clhs279;
    lhs(10,3)=clhs203*clhs278 + clhs210*clhs65 - clhs210;
    lhs(10,4)=DN(2,0)*clhs44 + DN(2,1)*clhs131 + DN(2,2)*clhs173 + clhs235;
    lhs(10,5)=DN(2,0)*clhs57 + DN(2,1)*clhs135 + DN(2,2)*clhs175 + clhs248;
    lhs(10,6)=DN(2,0)*clhs63 + DN(2,1)*clhs139 + DN(2,2)*clhs177 + clhs258 + clhs280;
    lhs(10,7)=clhs179*clhs278 + clhs268*clhs65 - clhs268;
    lhs(10,8)=DN(2,0)*clhs72 + DN(2,1)*clhs144 + DN(2,2)*clhs181 + clhs287;
    lhs(10,9)=DN(2,0)*clhs84 + DN(2,1)*clhs148 + DN(2,2)*clhs183 + clhs301;
    lhs(10,10)=DN(2,0)*clhs90 + DN(2,1)*clhs152 + DN(2,2)*clhs185 + clhs284 + clhs306*tau2;
    lhs(10,11)=DN(2,2)*clhs290;
    lhs(10,12)=DN(2,0)*clhs98 + DN(2,1)*clhs157 + DN(2,2)*clhs189 + clhs308;
    lhs(10,13)=DN(2,0)*clhs110 + DN(2,1)*clhs161 + DN(2,2)*clhs191 + clhs309;
    lhs(10,14)=DN(2,0)*clhs116 + DN(2,1)*clhs165 + DN(2,2)*clhs193 + clhs295 + clhs310;
    lhs(10,15)=DN(3,2)*clhs289 + clhs311*clhs65 - clhs311;
    lhs(11,0)=clhs14*clhs92 + clhs91;
    lhs(11,1)=clhs14*clhs154 + clhs153;
    lhs(11,2)=clhs14*clhs187 + clhs186;
    lhs(11,3)=clhs211;
    lhs(11,4)=clhs236 + clhs48*clhs92;
    lhs(11,5)=clhs154*clhs48 + clhs249;
    lhs(11,6)=clhs187*clhs48 + clhs259;
    lhs(11,7)=clhs272;
    lhs(11,8)=DN(2,0)*clhs313;
    lhs(11,9)=DN(2,1)*clhs313;
    lhs(11,10)=DN(2,2)*clhs313;
    lhs(11,11)=clhs197*clhs282 + clhs281*tau1 + clhs299*tau1 + clhs306*tau1;
    lhs(11,12)=clhs102*clhs92 + clhs314;
    lhs(11,13)=clhs102*clhs154 + clhs315;
    lhs(11,14)=clhs102*clhs187 + clhs316;
    lhs(11,15)=clhs317;
    lhs(12,0)=DN(3,0)*clhs1 + DN(3,1)*clhs3 + DN(3,2)*clhs5 + clhs320 + clhs93;
    lhs(12,1)=DN(3,0)*clhs21 + DN(3,1)*clhs23 + DN(3,2)*clhs26 + clhs155;
    lhs(12,2)=DN(3,0)*clhs28 + DN(3,1)*clhs30 + DN(3,2)*clhs32 + clhs188;
    lhs(12,3)=clhs199*clhs319 + clhs212*clhs65 - clhs212;
    lhs(12,4)=DN(3,0)*clhs40 + DN(3,1)*clhs42 + DN(3,2)*clhs44 + clhs237 + clhs321;
    lhs(12,5)=DN(3,0)*clhs52 + DN(3,1)*clhs54 + DN(3,2)*clhs57 + clhs250;
    lhs(12,6)=DN(3,0)*clhs59 + DN(3,1)*clhs61 + DN(3,2)*clhs63 + clhs260;
    lhs(12,7)=clhs273*clhs65 - clhs273 + clhs319*clhs66;
    lhs(12,8)=DN(3,0)*clhs68 + DN(3,1)*clhs70 + DN(3,2)*clhs72 + clhs291 + clhs322;
    lhs(12,9)=DN(3,0)*clhs79 + DN(3,1)*clhs81 + DN(3,2)*clhs84 + clhs302;
    lhs(12,10)=DN(3,0)*clhs86 + DN(3,1)*clhs88 + DN(3,2)*clhs90 + clhs308;
    lhs(12,11)=clhs314*clhs65 - clhs314 + clhs319*clhs92;
    lhs(12,12)=DN(3,0)*clhs94 + DN(3,1)*clhs96 + DN(3,2)*clhs98 + clhs323*tau2 + clhs325;
    lhs(12,13)=DN(3,0)*clhs105 + DN(3,1)*clhs107 + DN(3,2)*clhs110 + clhs327;
    lhs(12,14)=DN(3,0)*clhs112 + DN(3,1)*clhs114 + DN(3,2)*clhs116 + clhs328;
    lhs(12,15)=DN(3,0)*clhs329;
    lhs(13,0)=DN(3,0)*clhs3 + DN(3,1)*clhs118 + DN(3,2)*clhs119 + clhs104;
    lhs(13,1)=DN(3,0)*clhs23 + DN(3,1)*clhs121 + DN(3,2)*clhs123 + clhs158 + clhs320;
    lhs(13,2)=DN(3,0)*clhs30 + DN(3,1)*clhs126 + DN(3,2)*clhs128 + clhs190;
    lhs(13,3)=clhs201*clhs319 + clhs213*clhs65 - clhs213;
    lhs(13,4)=DN(3,0)*clhs42 + DN(3,1)*clhs130 + DN(3,2)*clhs131 + clhs240;
    lhs(13,5)=DN(3,0)*clhs54 + DN(3,1)*clhs133 + DN(3,2)*clhs135 + clhs251 + clhs321;
    lhs(13,6)=DN(3,0)*clhs61 + DN(3,1)*clhs137 + DN(3,2)*clhs139 + clhs261;
    lhs(13,7)=clhs141*clhs319 + clhs274*clhs65 - clhs274;
    lhs(13,8)=DN(3,0)*clhs70 + DN(3,1)*clhs143 + DN(3,2)*clhs144 + clhs296;
    lhs(13,9)=DN(3,0)*clhs81 + DN(3,1)*clhs146 + DN(3,2)*clhs148 + clhs303 + clhs322;
    lhs(13,10)=DN(3,0)*clhs88 + DN(3,1)*clhs150 + DN(3,2)*clhs152 + clhs309;
    lhs(13,11)=clhs154*clhs319 + clhs315*clhs65 - clhs315;
    lhs(13,12)=DN(3,0)*clhs96 + DN(3,1)*clhs156 + DN(3,2)*clhs157 + clhs327;
    lhs(13,13)=DN(3,0)*clhs107 + DN(3,1)*clhs159 + DN(3,2)*clhs161 + clhs325 + clhs330*tau2;
    lhs(13,14)=DN(3,0)*clhs114 + DN(3,1)*clhs163 + DN(3,2)*clhs165 + clhs331;
    lhs(13,15)=DN(3,1)*clhs329;
    lhs(14,0)=DN(3,0)*clhs5 + DN(3,1)*clhs119 + DN(3,2)*clhs167 + clhs111;
    lhs(14,1)=DN(3,0)*clhs26 + DN(3,1)*clhs123 + DN(3,2)*clhs168 + clhs162;
    lhs(14,2)=DN(3,0)*clhs32 + DN(3,1)*clhs128 + DN(3,2)*clhs170 + clhs192 + clhs320;
    lhs(14,3)=clhs203*clhs319 + clhs214*clhs65 - clhs214;
    lhs(14,4)=DN(3,0)*clhs44 + DN(3,1)*clhs131 + DN(3,2)*clhs173 + clhs241;
    lhs(14,5)=DN(3,0)*clhs57 + DN(3,1)*clhs135 + DN(3,2)*clhs175 + clhs252;
    lhs(14,6)=DN(3,0)*clhs63 + DN(3,1)*clhs139 + DN(3,2)*clhs177 + clhs262 + clhs321;
    lhs(14,7)=clhs179*clhs319 + clhs275*clhs65 - clhs275;
    lhs(14,8)=DN(3,0)*clhs72 + DN(3,1)*clhs144 + DN(3,2)*clhs181 + clhs297;
    lhs(14,9)=DN(3,0)*clhs84 + DN(3,1)*clhs148 + DN(3,2)*clhs183 + clhs304;
    lhs(14,10)=DN(3,0)*clhs90 + DN(3,1)*clhs152 + DN(3,2)*clhs185 + clhs310 + clhs322;
    lhs(14,11)=clhs187*clhs319 + clhs316*clhs65 - clhs316;
    lhs(14,12)=DN(3,0)*clhs98 + DN(3,1)*clhs157 + DN(3,2)*clhs189 + clhs328;
    lhs(14,13)=DN(3,0)*clhs110 + DN(3,1)*clhs161 + DN(3,2)*clhs191 + clhs331;
    lhs(14,14)=DN(3,0)*clhs116 + DN(3,1)*clhs165 + DN(3,2)*clhs193 + clhs325 + clhs332*tau2;
    lhs(14,15)=DN(3,2)*clhs329;
    lhs(15,0)=DN(3,0)*clhs195 + clhs117;
    lhs(15,1)=DN(3,1)*clhs195 + clhs166;
    lhs(15,2)=DN(3,2)*clhs195 + clhs194;
    lhs(15,3)=clhs215;
    lhs(15,4)=DN(3,0)*clhs264 + clhs242;
    lhs(15,5)=DN(3,1)*clhs264 + clhs253;
    lhs(15,6)=DN(3,2)*clhs264 + clhs263;
    lhs(15,7)=clhs276;
    lhs(15,8)=DN(3,0)*clhs312 + clhs298;
    lhs(15,9)=DN(3,1)*clhs312 + clhs305;
    lhs(15,10)=DN(3,2)*clhs312 + clhs311;
    lhs(15,11)=clhs317;
    lhs(15,12)=DN(3,0)*clhs333;
    lhs(15,13)=DN(3,1)*clhs333;
    lhs(15,14)=DN(3,2)*clhs333;
    lhs(15,15)=clhs197*clhs324 + clhs323*tau1 + clhs330*tau1 + clhs332*tau1;

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
    // const double& bdf1 = data.bdf1;
    // const double& bdf2 = data.bdf2;
    const double& delta_t = data.delta_t;
    const double& dyn_tau_coeff = data.dyn_tau_coeff;

    const bounded_matrix<double,nnodes,dim>& v = data.v;
    // const bounded_matrix<double,nnodes,dim>& vn = data.vn;
    // const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
    const bounded_matrix<double,nnodes,dim>& vmesh = data.vmesh;
    const bounded_matrix<double,nnodes,dim>& vconv = v - vmesh;
    // const bounded_matrix<double,nnodes,dim>& f = data.f;
    // const array_1d<double,nnodes>& p = data.p;
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
    const double clhs1 =             C(0,0)*DN(0,0) + C(0,2)*DN(0,1);
    const double clhs2 =             C(0,2)*DN(0,0);
    const double clhs3 =             C(2,2)*DN(0,1) + clhs2;
    const double clhs4 =             pow(N[0], 2);
    const double clhs5 =             bdf0*rho;
    const double clhs6 =             N[0]*rho;
    const double clhs7 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
    const double clhs8 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
    const double clhs9 =             DN(0,0)*clhs7 + DN(0,1)*clhs8;
    const double clhs10 =             pow(rho, 2);
    const double clhs11 =             N[0]*bdf0 + clhs9;
    const double clhs12 =             clhs10*clhs11*tau1;
    const double clhs13 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
    const double clhs14 =             N[0]*clhs13 + clhs9;
    const double clhs15 =             clhs12*clhs14 + clhs4*clhs5 + clhs6*clhs9;
    const double clhs16 =             DN(0,0)*tau2;
    const double clhs17 =             DN(0,1)*clhs16;
    const double clhs18 =             C(0,1)*DN(0,1) + clhs2;
    const double clhs19 =             C(1,2)*DN(0,1);
    const double clhs20 =             C(2,2)*DN(0,0) + clhs19;
    const double clhs21 =             pow(c, -2);
    const double clhs22 =             1.0/rho;
    const double clhs23 =             N[0]*bdf0*clhs21*clhs22;
    const double clhs24 =             rho*tau1;
    const double clhs25 =             clhs14*clhs24;
    const double clhs26 =             -N[0] + clhs23*tau2 + clhs25;
    const double clhs27 =             DN(1,0)*clhs16;
    const double clhs28 =             C(0,0)*DN(1,0) + C(0,2)*DN(1,1);
    const double clhs29 =             C(0,2)*DN(1,0);
    const double clhs30 =             C(2,2)*DN(1,1) + clhs29;
    const double clhs31 =             N[0]*bdf0*rho;
    const double clhs32 =             N[1]*clhs31;
    const double clhs33 =             DN(1,0)*clhs7 + DN(1,1)*clhs8;
    const double clhs34 =             N[1]*bdf0;
    const double clhs35 =             clhs33 + clhs34;
    const double clhs36 =             clhs10*clhs14*tau1;
    const double clhs37 =             clhs32 + clhs33*clhs6 + clhs35*clhs36;
    const double clhs38 =             DN(1,1)*clhs16;
    const double clhs39 =             C(0,1)*DN(1,1) + clhs29;
    const double clhs40 =             C(1,2)*DN(1,1);
    const double clhs41 =             C(2,2)*DN(1,0) + clhs40;
    const double clhs42 =             DN(0,0)*N[1];
    const double clhs43 =             bdf0*clhs21*clhs22*tau2;
    const double clhs44 =             DN(1,0)*rho*tau1;
    const double clhs45 =             DN(2,0)*clhs16;
    const double clhs46 =             C(0,0)*DN(2,0) + C(0,2)*DN(2,1);
    const double clhs47 =             C(0,2)*DN(2,0);
    const double clhs48 =             C(2,2)*DN(2,1) + clhs47;
    const double clhs49 =             N[2]*clhs31;
    const double clhs50 =             DN(2,0)*clhs7 + DN(2,1)*clhs8;
    const double clhs51 =             N[2]*bdf0;
    const double clhs52 =             clhs50 + clhs51;
    const double clhs53 =             clhs36*clhs52 + clhs49 + clhs50*clhs6;
    const double clhs54 =             DN(2,1)*clhs16;
    const double clhs55 =             C(0,1)*DN(2,1) + clhs47;
    const double clhs56 =             C(1,2)*DN(2,1);
    const double clhs57 =             C(2,2)*DN(2,0) + clhs56;
    const double clhs58 =             DN(0,0)*N[2];
    const double clhs59 =             C(0,1)*DN(0,0) + clhs19;
    const double clhs60 =             pow(DN(0,1), 2);
    const double clhs61 =             C(1,1)*DN(0,1) + C(1,2)*DN(0,0);
    const double clhs62 =             DN(0,1)*tau2;
    const double clhs63 =             DN(1,0)*clhs62;
    const double clhs64 =             C(0,1)*DN(1,0) + clhs40;
    const double clhs65 =             DN(1,1)*clhs62;
    const double clhs66 =             C(1,1)*DN(1,1) + C(1,2)*DN(1,0);
    const double clhs67 =             DN(0,1)*N[1];
    const double clhs68 =             DN(1,1)*rho*tau1;
    const double clhs69 =             DN(2,0)*clhs62;
    const double clhs70 =             C(0,1)*DN(2,0) + clhs56;
    const double clhs71 =             DN(2,1)*clhs62;
    const double clhs72 =             C(1,1)*DN(2,1) + C(1,2)*DN(2,0);
    const double clhs73 =             DN(0,1)*N[2];
    const double clhs74 =             clhs11*clhs24;
    const double clhs75 =             N[0] + clhs74;
    const double clhs76 =             bdf0*clhs21*clhs22;
    const double clhs77 =             DN(1,0)*N[0];
    const double clhs78 =             DN(0,0)*rho*tau1;
    const double clhs79 =             DN(1,1)*N[0];
    const double clhs80 =             DN(0,1)*rho*tau1;
    const double clhs81 =             DN(0,0)*tau1;
    const double clhs82 =             DN(0,1)*tau1;
    const double clhs83 =             DN(1,0)*clhs81 + DN(1,1)*clhs82 + N[1]*clhs23;
    const double clhs84 =             DN(2,0)*N[0];
    const double clhs85 =             DN(2,1)*N[0];
    const double clhs86 =             DN(2,0)*clhs81 + DN(2,1)*clhs82 + N[2]*clhs23;
    const double clhs87 =             N[1]*rho;
    const double clhs88 =             N[1]*clhs13 + clhs33;
    const double clhs89 =             clhs12*clhs88 + clhs32 + clhs87*clhs9;
    const double clhs90 =             pow(DN(1,0), 2);
    const double clhs91 =             pow(N[1], 2);
    const double clhs92 =             clhs10*clhs35*tau1;
    const double clhs93 =             clhs33*clhs87 + clhs5*clhs91 + clhs88*clhs92;
    const double clhs94 =             DN(1,0)*tau2;
    const double clhs95 =             DN(1,1)*clhs94;
    const double clhs96 =             clhs21*clhs22*tau2;
    const double clhs97 =             clhs24*clhs88;
    const double clhs98 =             -N[1] + clhs34*clhs96 + clhs97;
    const double clhs99 =             DN(2,0)*clhs94;
    const double clhs100 =             N[1]*N[2]*bdf0;
    const double clhs101 =             clhs100*rho;
    const double clhs102 =             clhs10*clhs52*tau1;
    const double clhs103 =             clhs101 + clhs102*clhs88 + clhs50*clhs87;
    const double clhs104 =             DN(2,1)*clhs94;
    const double clhs105 =             DN(1,0)*N[2];
    const double clhs106 =             pow(DN(1,1), 2);
    const double clhs107 =             DN(1,1)*tau2;
    const double clhs108 =             DN(2,0)*clhs107;
    const double clhs109 =             DN(2,1)*clhs107;
    const double clhs110 =             DN(1,1)*N[2];
    const double clhs111 =             clhs24*clhs35;
    const double clhs112 =             N[1] + clhs111;
    const double clhs113 =             DN(2,0)*N[1];
    const double clhs114 =             DN(2,1)*N[1];
    const double clhs115 =             DN(2,0)*(DN(1,0)*tau1) + DN(2,1)*(DN(1,1)*tau1) + clhs100*clhs21*clhs22;
    const double clhs116 =             N[2]*rho;
    const double clhs117 =             N[2]*clhs13 + clhs50;
    const double clhs118 =             clhs116*clhs9 + clhs117*clhs12 + clhs49;
    const double clhs119 =             clhs101 + clhs116*clhs33 + clhs117*clhs92;
    const double clhs120 =             pow(DN(2,0), 2);
    const double clhs121 =             pow(N[2], 2);
    const double clhs122 =             clhs102*clhs117 + clhs116*clhs50 + clhs121*clhs5;
    const double clhs123 =             DN(2,0)*DN(2,1)*tau2;
    const double clhs124 =             -N[2] + clhs117*clhs24 + clhs51*clhs96;
    const double clhs125 =             pow(DN(2,1), 2);
    const double clhs126 =             N[2] + clhs24*clhs52;

    lhs(0,0)=DN(0,0)*clhs1 + DN(0,1)*clhs3 + clhs0*tau2 + clhs15;
    lhs(0,1)=DN(0,0)*clhs18 + DN(0,1)*clhs20 + clhs17;
    lhs(0,2)=DN(0,0)*clhs26;
    lhs(0,3)=DN(0,0)*clhs28 + DN(0,1)*clhs30 + clhs27 + clhs37;
    lhs(0,4)=DN(0,0)*clhs39 + DN(0,1)*clhs41 + clhs38;
    lhs(0,5)=clhs14*clhs44 + clhs42*clhs43 - clhs42;
    lhs(0,6)=DN(0,0)*clhs46 + DN(0,1)*clhs48 + clhs45 + clhs53;
    lhs(0,7)=DN(0,0)*clhs55 + DN(0,1)*clhs57 + clhs54;
    lhs(0,8)=DN(2,0)*clhs25 + clhs43*clhs58 - clhs58;
    lhs(1,0)=DN(0,0)*clhs3 + DN(0,1)*clhs59 + clhs17;
    lhs(1,1)=DN(0,0)*clhs20 + DN(0,1)*clhs61 + clhs15 + clhs60*tau2;
    lhs(1,2)=DN(0,1)*clhs26;
    lhs(1,3)=DN(0,0)*clhs30 + DN(0,1)*clhs64 + clhs63;
    lhs(1,4)=DN(0,0)*clhs41 + DN(0,1)*clhs66 + clhs37 + clhs65;
    lhs(1,5)=clhs14*clhs68 + clhs43*clhs67 - clhs67;
    lhs(1,6)=DN(0,0)*clhs48 + DN(0,1)*clhs70 + clhs69;
    lhs(1,7)=DN(0,0)*clhs57 + DN(0,1)*clhs72 + clhs53 + clhs71;
    lhs(1,8)=DN(2,1)*clhs25 + clhs43*clhs73 - clhs73;
    lhs(2,0)=DN(0,0)*clhs75;
    lhs(2,1)=DN(0,1)*clhs75;
    lhs(2,2)=clhs0*tau1 + clhs4*clhs76 + clhs60*tau1;
    lhs(2,3)=clhs35*clhs78 + clhs77;
    lhs(2,4)=clhs35*clhs80 + clhs79;
    lhs(2,5)=clhs83;
    lhs(2,6)=clhs52*clhs78 + clhs84;
    lhs(2,7)=clhs52*clhs80 + clhs85;
    lhs(2,8)=clhs86;
    lhs(3,0)=DN(1,0)*clhs1 + DN(1,1)*clhs3 + clhs27 + clhs89;
    lhs(3,1)=DN(1,0)*clhs18 + DN(1,1)*clhs20 + clhs63;
    lhs(3,2)=clhs43*clhs77 - clhs77 + clhs78*clhs88;
    lhs(3,3)=DN(1,0)*clhs28 + DN(1,1)*clhs30 + clhs90*tau2 + clhs93;
    lhs(3,4)=DN(1,0)*clhs39 + DN(1,1)*clhs41 + clhs95;
    lhs(3,5)=DN(1,0)*clhs98;
    lhs(3,6)=DN(1,0)*clhs46 + DN(1,1)*clhs48 + clhs103 + clhs99;
    lhs(3,7)=DN(1,0)*clhs55 + DN(1,1)*clhs57 + clhs104;
    lhs(3,8)=DN(2,0)*clhs97 + clhs105*clhs43 - clhs105;
    lhs(4,0)=DN(1,0)*clhs3 + DN(1,1)*clhs59 + clhs38;
    lhs(4,1)=DN(1,0)*clhs20 + DN(1,1)*clhs61 + clhs65 + clhs89;
    lhs(4,2)=clhs43*clhs79 - clhs79 + clhs80*clhs88;
    lhs(4,3)=DN(1,0)*clhs30 + DN(1,1)*clhs64 + clhs95;
    lhs(4,4)=DN(1,0)*clhs41 + DN(1,1)*clhs66 + clhs106*tau2 + clhs93;
    lhs(4,5)=DN(1,1)*clhs98;
    lhs(4,6)=DN(1,0)*clhs48 + DN(1,1)*clhs70 + clhs108;
    lhs(4,7)=DN(1,0)*clhs57 + DN(1,1)*clhs72 + clhs103 + clhs109;
    lhs(4,8)=DN(2,1)*clhs97 + clhs110*clhs43 - clhs110;
    lhs(5,0)=clhs11*clhs44 + clhs42;
    lhs(5,1)=clhs11*clhs68 + clhs67;
    lhs(5,2)=clhs83;
    lhs(5,3)=DN(1,0)*clhs112;
    lhs(5,4)=DN(1,1)*clhs112;
    lhs(5,5)=clhs106*tau1 + clhs76*clhs91 + clhs90*tau1;
    lhs(5,6)=clhs113 + clhs44*clhs52;
    lhs(5,7)=clhs114 + clhs52*clhs68;
    lhs(5,8)=clhs115;
    lhs(6,0)=DN(2,0)*clhs1 + DN(2,1)*clhs3 + clhs118 + clhs45;
    lhs(6,1)=DN(2,0)*clhs18 + DN(2,1)*clhs20 + clhs69;
    lhs(6,2)=clhs117*clhs78 + clhs43*clhs84 - clhs84;
    lhs(6,3)=DN(2,0)*clhs28 + DN(2,1)*clhs30 + clhs119 + clhs99;
    lhs(6,4)=DN(2,0)*clhs39 + DN(2,1)*clhs41 + clhs108;
    lhs(6,5)=clhs113*clhs43 - clhs113 + clhs117*clhs44;
    lhs(6,6)=DN(2,0)*clhs46 + DN(2,1)*clhs48 + clhs120*tau2 + clhs122;
    lhs(6,7)=DN(2,0)*clhs55 + DN(2,1)*clhs57 + clhs123;
    lhs(6,8)=DN(2,0)*clhs124;
    lhs(7,0)=DN(2,0)*clhs3 + DN(2,1)*clhs59 + clhs54;
    lhs(7,1)=DN(2,0)*clhs20 + DN(2,1)*clhs61 + clhs118 + clhs71;
    lhs(7,2)=clhs117*clhs80 + clhs43*clhs85 - clhs85;
    lhs(7,3)=DN(2,0)*clhs30 + DN(2,1)*clhs64 + clhs104;
    lhs(7,4)=DN(2,0)*clhs41 + DN(2,1)*clhs66 + clhs109 + clhs119;
    lhs(7,5)=clhs114*clhs43 - clhs114 + clhs117*clhs68;
    lhs(7,6)=DN(2,0)*clhs48 + DN(2,1)*clhs70 + clhs123;
    lhs(7,7)=DN(2,0)*clhs57 + DN(2,1)*clhs72 + clhs122 + clhs125*tau2;
    lhs(7,8)=DN(2,1)*clhs124;
    lhs(8,0)=DN(2,0)*clhs74 + clhs58;
    lhs(8,1)=DN(2,1)*clhs74 + clhs73;
    lhs(8,2)=clhs86;
    lhs(8,3)=DN(2,0)*clhs111 + clhs105;
    lhs(8,4)=DN(2,1)*clhs111 + clhs110;
    lhs(8,5)=clhs115;
    lhs(8,6)=DN(2,0)*clhs126;
    lhs(8,7)=DN(2,1)*clhs126;
    lhs(8,8)=clhs120*tau1 + clhs121*clhs76 + clhs125*tau1;

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
    const double crhs4 =             DN(0,0)*v(0,0);
    const double crhs5 =             DN(1,0)*v(1,0);
    const double crhs6 =             DN(2,0)*v(2,0);
    const double crhs7 =             DN(3,0)*v(3,0);
    const double crhs8 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
    const double crhs9 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
    const double crhs10 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
    const double crhs11 =             crhs10*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0)) + crhs8*(crhs4 + crhs5 + crhs6 + crhs7) + crhs9*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0));
    const double crhs12 =             DN(0,2)*v(0,2) + DN(1,2)*v(1,2) + DN(2,2)*v(2,2) + DN(3,2)*v(3,2);
    const double crhs13 =             DN(0,1)*v(0,1);
    const double crhs14 =             DN(1,1)*v(1,1);
    const double crhs15 =             DN(2,1)*v(2,1);
    const double crhs16 =             DN(3,1)*v(3,1);
    const double crhs17 =             crhs12 + crhs13 + crhs14 + crhs15 + crhs16 + crhs4 + crhs5 + crhs6 + crhs7;
    const double crhs18 =             (N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]) + N[3]*(bdf0*p[3] + bdf1*pn[3] + bdf2*pnn[3]))/(pow(c, 2)*rho);
    const double crhs19 =             tau2*(crhs17 + crhs18);
    const double crhs20 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(0,2)*vconv(0,2) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(1,2)*vconv(1,2) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1) + DN(2,2)*vconv(2,2) + DN(3,0)*vconv(3,0) + DN(3,1)*vconv(3,1) + DN(3,2)*vconv(3,2);
    const double crhs21 =             rho*(DN(0,0)*crhs8 + DN(0,1)*crhs9 + DN(0,2)*crhs10 + N[0]*crhs20);
    const double crhs22 =             tau1*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + rho*(crhs11 + crhs3));
    const double crhs23 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
    const double crhs24 =             N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1));
    const double crhs25 =             crhs10*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1)) + crhs8*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs9*(crhs13 + crhs14 + crhs15 + crhs16);
    const double crhs26 =             tau1*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs23 + rho*(crhs24 + crhs25));
    const double crhs27 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
    const double crhs28 =             N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2));
    const double crhs29 =             crhs10*crhs12 + crhs8*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs9*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2));
    const double crhs30 =             tau1*(DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs27 + rho*(crhs28 + crhs29));
    const double crhs31 =             N[1]*rho;
    const double crhs32 =             rho*(DN(1,0)*crhs8 + DN(1,1)*crhs9 + DN(1,2)*crhs10 + N[1]*crhs20);
    const double crhs33 =             N[2]*rho;
    const double crhs34 =             rho*(DN(2,0)*crhs8 + DN(2,1)*crhs9 + DN(2,2)*crhs10 + N[2]*crhs20);
    const double crhs35 =             N[3]*rho;
    const double crhs36 =             rho*(DN(3,0)*crhs8 + DN(3,1)*crhs9 + DN(3,2)*crhs10 + N[3]*crhs20);

    rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs19 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - crhs11*crhs2 - crhs2*crhs3 - crhs21*crhs22;
    rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs19 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs23 - crhs2*crhs24 - crhs2*crhs25 - crhs21*crhs26;
    rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs19 - DN(0,2)*stress[2] + N[0]*crhs27 - crhs2*crhs28 - crhs2*crhs29 - crhs21*crhs30;
    rhs[3]=-DN(0,0)*crhs22 - DN(0,1)*crhs26 - DN(0,2)*crhs30 - N[0]*crhs17 - N[0]*crhs18;
    rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs19 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - crhs11*crhs31 - crhs22*crhs32 - crhs3*crhs31;
    rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs19 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs23 - crhs24*crhs31 - crhs25*crhs31 - crhs26*crhs32;
    rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs19 - DN(1,2)*stress[2] + N[1]*crhs27 - crhs28*crhs31 - crhs29*crhs31 - crhs30*crhs32;
    rhs[7]=-DN(1,0)*crhs22 - DN(1,1)*crhs26 - DN(1,2)*crhs30 - N[1]*crhs17 - N[1]*crhs18;
    rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs19 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - crhs11*crhs33 - crhs22*crhs34 - crhs3*crhs33;
    rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs19 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs23 - crhs24*crhs33 - crhs25*crhs33 - crhs26*crhs34;
    rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs19 - DN(2,2)*stress[2] + N[2]*crhs27 - crhs28*crhs33 - crhs29*crhs33 - crhs30*crhs34;
    rhs[11]=-DN(2,0)*crhs22 - DN(2,1)*crhs26 - DN(2,2)*crhs30 - N[2]*crhs17 - N[2]*crhs18;
    rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs19 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - crhs11*crhs35 - crhs22*crhs36 - crhs3*crhs35;
    rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs19 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs23 - crhs24*crhs35 - crhs25*crhs35 - crhs26*crhs36;
    rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs19 - DN(3,2)*stress[2] + N[3]*crhs27 - crhs28*crhs35 - crhs29*crhs35 - crhs30*crhs36;
    rhs[15]=-DN(3,0)*crhs22 - DN(3,1)*crhs26 - DN(3,2)*crhs30 - N[3]*crhs17 - N[3]*crhs18;

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
    const double crhs5 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
    const double crhs6 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
    const double crhs7 =             crhs4*crhs5 + crhs6*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0));
    const double crhs8 =             DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1);
    const double crhs9 =             crhs4 + crhs8;
    const double crhs10 =             (N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]))/(pow(c, 2)*rho);
    const double crhs11 =             tau2*(crhs10 + crhs9);
    const double crhs12 =             DN(0,0)*vconv(0,0) + DN(0,1)*vconv(0,1) + DN(1,0)*vconv(1,0) + DN(1,1)*vconv(1,1) + DN(2,0)*vconv(2,0) + DN(2,1)*vconv(2,1);
    const double crhs13 =             rho*(DN(0,0)*crhs5 + DN(0,1)*crhs6 + N[0]*crhs12);
    const double crhs14 =             tau1*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + rho*(crhs3 + crhs7));
    const double crhs15 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
    const double crhs16 =             N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1));
    const double crhs17 =             crhs5*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + crhs6*crhs8;
    const double crhs18 =             tau1*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs15 + rho*(crhs16 + crhs17));
    const double crhs19 =             N[1]*rho;
    const double crhs20 =             rho*(DN(1,0)*crhs5 + DN(1,1)*crhs6 + N[1]*crhs12);
    const double crhs21 =             N[2]*rho;
    const double crhs22 =             rho*(DN(2,0)*crhs5 + DN(2,1)*crhs6 + N[2]*crhs12);

    rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs11 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - crhs13*crhs14 - crhs2*crhs3 - crhs2*crhs7;
    rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs11 - DN(0,1)*stress[1] + N[0]*crhs15 - crhs13*crhs18 - crhs16*crhs2 - crhs17*crhs2;
    rhs[2]=-DN(0,0)*crhs14 - DN(0,1)*crhs18 - N[0]*crhs10 - N[0]*crhs9;
    rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs11 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - crhs14*crhs20 - crhs19*crhs3 - crhs19*crhs7;
    rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs11 - DN(1,1)*stress[1] + N[1]*crhs15 - crhs16*crhs19 - crhs17*crhs19 - crhs18*crhs20;
    rhs[5]=-DN(1,0)*crhs14 - DN(1,1)*crhs18 - N[1]*crhs10 - N[1]*crhs9;
    rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs11 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - crhs14*crhs22 - crhs21*crhs3 - crhs21*crhs7;
    rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs11 - DN(2,1)*stress[1] + N[2]*crhs15 - crhs16*crhs21 - crhs17*crhs21 - crhs18*crhs22;
    rhs[8]=-DN(2,0)*crhs14 - DN(2,1)*crhs18 - N[2]*crhs10 - N[2]*crhs9;

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
    const double cv_s_gauss0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
    const double cv_s_gauss1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
    const double cv_s_gauss2 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);

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
    const double cv_s_gauss0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
    const double cv_s_gauss1 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);

    v_s_gauss[0]=-tau1*(DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + rho*(-N[0]*f(0,0) + N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) - N[1]*f(1,0) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) - N[2]*f(2,0) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + cv_s_gauss0*(DN(0,0)*v(0,0) + DN(1,0)*v(1,0) + DN(2,0)*v(2,0)) + cv_s_gauss1*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0))));
    v_s_gauss[1]=-tau1*(DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + rho*(-N[0]*f(0,1) + N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) - N[1]*f(1,1) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) - N[2]*f(2,1) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + cv_s_gauss0*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1)) + cv_s_gauss1*(DN(0,1)*v(0,1) + DN(1,1)*v(1,1) + DN(2,1)*v(2,1))));

    const double v_gauss_norm = norm_2(v_gauss);
    const double v_s_gauss_norm = norm_2(v_s_gauss);

    return v_s_gauss_norm/v_gauss_norm;
}

}
