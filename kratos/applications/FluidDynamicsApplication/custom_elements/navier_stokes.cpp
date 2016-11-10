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
        const int strain_size = 6;
        
        const double rho = inner_prod(data.N, data.rho);        // Density
        const double nu = inner_prod(data.N, data.nu);          // Kinematic viscosity
        const double h = data.h;                                // Characteristic element size
        
        const double& bdf0 = data.bdf0;
        //~ const double& bdf1 = data.bdf1;
        //~ const double& bdf2 = data.bdf2;
        const double& delta_t = data.delta_t;
        const double& dyn_tau_coeff = data.dyn_tau_coeff;
        
        const bounded_matrix<double,nnodes,dim>& v = data.v;
        //~ const bounded_matrix<double,nnodes,dim>& vn = data.vn;
        //~ const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
        const bounded_matrix<double,nnodes,dim>& vmesh = data.vmesh;
        const bounded_matrix<double,nnodes,dim>& vconv = v - vmesh;
        //~ const bounded_matrix<double,nnodes,dim>& f = data.f;
        //~ const array_1d<double,nnodes>& p = data.p;
        const array_1d<double,strain_size>& stress = data.stress;
        
        // Get constitutive matrix 
        //~ const Matrix& C = data.C;
        
        // Get shape function values
        const array_1d<double,nnodes>& N = data.N;
        const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;
        
        const array_1d<double,dim> vconv_gauss = prod(trans(vconv), N);
        
        const double vconv_norm = norm_2(vconv_gauss);
                
        // Stabilization parameters
        const double tau1 = 1.0/((rho*dyn_tau_coeff)/delta_t + (2*rho*vconv_norm)/h + (4*rho*nu)/(h*h));
        const double tau2 = (rho*nu) + 0.5*h*vconv_norm;
        
        const double clhs0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
        const double clhs1 =             DN(0,0)*clhs0;
        const double clhs2 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
        const double clhs3 =             DN(0,1)*clhs2;
        const double clhs4 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
        const double clhs5 =             DN(0,2)*clhs4;
        const double clhs6 =             clhs1 + clhs3 + clhs5;
        const double clhs7 =             pow(N[0], 2)*bdf0 + N[0]*clhs6;
        const double clhs8 =             pow(DN(0,0), 2);
        const double clhs9 =             rho*tau1;
        const double clhs10 =             N[0]*bdf0;
        const double clhs11 =             clhs9*(clhs10 + clhs6);
        const double clhs12 =             DN(0,1)*rho;
        const double clhs13 =             DN(0,0)*tau2;
        const double clhs14 =             clhs0*clhs11 + clhs13;
        const double clhs15 =             DN(0,2)*rho;
        const double clhs16 =             pow(DN(0,1), 2);
        const double clhs17 =             pow(DN(0,2), 2);
        const double clhs18 =             clhs9*(clhs16 + clhs17 + clhs8);
        const double clhs19 =             DN(0,0)*DN(1,0);
        const double clhs20 =             N[1]*clhs10;
        const double clhs21 =             clhs19*tau2 + clhs20;
        const double clhs22 =             DN(1,0)*clhs0;
        const double clhs23 =             DN(1,1)*clhs2;
        const double clhs24 =             DN(1,2)*clhs4;
        const double clhs25 =             clhs22 + clhs23 + clhs24;
        const double clhs26 =             N[0]*clhs25;
        const double clhs27 =             N[1]*bdf0;
        const double clhs28 =             clhs9*(clhs25 + clhs27);
        const double clhs29 =             DN(1,1)*clhs13;
        const double clhs30 =             DN(0,1)*clhs28;
        const double clhs31 =             DN(1,2)*clhs13;
        const double clhs32 =             DN(0,2)*clhs28;
        const double clhs33 =             DN(0,0)*N[1];
        const double clhs34 =             DN(0,1)*DN(1,1);
        const double clhs35 =             DN(0,2)*DN(1,2);
        const double clhs36 =             clhs9*(clhs19 + clhs34 + clhs35);
        const double clhs37 =             clhs0*clhs36;
        const double clhs38 =             DN(0,0)*DN(2,0);
        const double clhs39 =             N[2]*clhs10;
        const double clhs40 =             clhs38*tau2 + clhs39;
        const double clhs41 =             DN(2,0)*clhs0;
        const double clhs42 =             DN(2,1)*clhs2;
        const double clhs43 =             DN(2,2)*clhs4;
        const double clhs44 =             clhs41 + clhs42 + clhs43;
        const double clhs45 =             N[0]*clhs44;
        const double clhs46 =             N[2]*bdf0;
        const double clhs47 =             clhs9*(clhs44 + clhs46);
        const double clhs48 =             DN(2,1)*clhs13;
        const double clhs49 =             DN(0,1)*clhs47;
        const double clhs50 =             DN(2,2)*clhs13;
        const double clhs51 =             DN(0,2)*clhs47;
        const double clhs52 =             DN(0,0)*N[2];
        const double clhs53 =             DN(0,1)*DN(2,1);
        const double clhs54 =             DN(0,2)*DN(2,2);
        const double clhs55 =             clhs9*(clhs38 + clhs53 + clhs54);
        const double clhs56 =             clhs0*clhs55;
        const double clhs57 =             DN(0,0)*DN(3,0);
        const double clhs58 =             N[3]*clhs10;
        const double clhs59 =             clhs57*tau2 + clhs58;
        const double clhs60 =             DN(3,0)*clhs0;
        const double clhs61 =             DN(3,1)*clhs2;
        const double clhs62 =             DN(3,2)*clhs4;
        const double clhs63 =             clhs60 + clhs61 + clhs62;
        const double clhs64 =             N[0]*clhs63;
        const double clhs65 =             clhs9*(N[3]*bdf0 + clhs63);
        const double clhs66 =             DN(3,1)*clhs13;
        const double clhs67 =             DN(0,1)*clhs65;
        const double clhs68 =             DN(3,2)*clhs13;
        const double clhs69 =             DN(0,2)*clhs65;
        const double clhs70 =             DN(0,0)*N[3];
        const double clhs71 =             DN(0,1)*DN(3,1);
        const double clhs72 =             DN(0,2)*DN(3,2);
        const double clhs73 =             clhs9*(clhs57 + clhs71 + clhs72);
        const double clhs74 =             clhs0*clhs73;
        const double clhs75 =             DN(0,0)*rho;
        const double clhs76 =             DN(0,1)*tau2;
        const double clhs77 =             clhs11*clhs2 + clhs76;
        const double clhs78 =             DN(1,0)*clhs76;
        const double clhs79 =             DN(0,0)*clhs28;
        const double clhs80 =             clhs20 + clhs34*tau2;
        const double clhs81 =             DN(1,2)*clhs76;
        const double clhs82 =             DN(0,1)*N[1];
        const double clhs83 =             clhs2*clhs36;
        const double clhs84 =             DN(2,0)*clhs76;
        const double clhs85 =             DN(0,0)*clhs47;
        const double clhs86 =             clhs39 + clhs53*tau2;
        const double clhs87 =             DN(2,2)*clhs76;
        const double clhs88 =             DN(0,1)*N[2];
        const double clhs89 =             clhs2*clhs55;
        const double clhs90 =             DN(3,0)*clhs76;
        const double clhs91 =             DN(0,0)*clhs65;
        const double clhs92 =             clhs58 + clhs71*tau2;
        const double clhs93 =             DN(3,2)*clhs76;
        const double clhs94 =             DN(0,1)*N[3];
        const double clhs95 =             clhs2*clhs73;
        const double clhs96 =             DN(0,2)*tau2;
        const double clhs97 =             clhs11*clhs4 + clhs96;
        const double clhs98 =             DN(1,0)*clhs96;
        const double clhs99 =             DN(1,1)*clhs96;
        const double clhs100 =             clhs20 + clhs35*tau2;
        const double clhs101 =             DN(0,2)*N[1];
        const double clhs102 =             clhs36*clhs4;
        const double clhs103 =             DN(2,0)*clhs96;
        const double clhs104 =             DN(2,1)*clhs96;
        const double clhs105 =             clhs39 + clhs54*tau2;
        const double clhs106 =             DN(0,2)*N[2];
        const double clhs107 =             clhs4*clhs55;
        const double clhs108 =             DN(3,0)*clhs96;
        const double clhs109 =             DN(3,1)*clhs96;
        const double clhs110 =             clhs58 + clhs72*tau2;
        const double clhs111 =             DN(0,2)*N[3];
        const double clhs112 =             clhs4*clhs73;
        const double clhs113 =             rho*(N[0] + clhs11);
        const double clhs114 =             DN(1,0)*N[0];
        const double clhs115 =             DN(1,1)*N[0];
        const double clhs116 =             DN(1,2)*N[0];
        const double clhs117 =             DN(2,0)*N[0];
        const double clhs118 =             DN(2,1)*N[0];
        const double clhs119 =             DN(2,2)*N[0];
        const double clhs120 =             DN(3,0)*N[0];
        const double clhs121 =             DN(3,1)*N[0];
        const double clhs122 =             DN(3,2)*N[0];
        const double clhs123 =             N[1]*clhs6;
        const double clhs124 =             DN(1,1)*clhs11;
        const double clhs125 =             DN(1,2)*clhs11;
        const double clhs126 =             pow(N[1], 2)*bdf0 + N[1]*clhs25;
        const double clhs127 =             pow(DN(1,0), 2);
        const double clhs128 =             DN(1,1)*rho;
        const double clhs129 =             DN(1,0)*tau2;
        const double clhs130 =             clhs0*clhs28 + clhs129;
        const double clhs131 =             DN(1,2)*rho;
        const double clhs132 =             pow(DN(1,1), 2);
        const double clhs133 =             pow(DN(1,2), 2);
        const double clhs134 =             clhs9*(clhs127 + clhs132 + clhs133);
        const double clhs135 =             DN(1,0)*DN(2,0);
        const double clhs136 =             N[2]*clhs27;
        const double clhs137 =             clhs135*tau2 + clhs136;
        const double clhs138 =             N[1]*clhs44;
        const double clhs139 =             DN(2,1)*clhs129;
        const double clhs140 =             DN(1,1)*clhs47;
        const double clhs141 =             DN(2,2)*clhs129;
        const double clhs142 =             DN(1,2)*clhs47;
        const double clhs143 =             DN(1,0)*N[2];
        const double clhs144 =             DN(1,1)*DN(2,1);
        const double clhs145 =             DN(1,2)*DN(2,2);
        const double clhs146 =             clhs9*(clhs135 + clhs144 + clhs145);
        const double clhs147 =             clhs0*clhs146;
        const double clhs148 =             DN(1,0)*DN(3,0);
        const double clhs149 =             N[3]*clhs27;
        const double clhs150 =             clhs148*tau2 + clhs149;
        const double clhs151 =             N[1]*clhs63;
        const double clhs152 =             DN(3,1)*clhs129;
        const double clhs153 =             DN(1,1)*clhs65;
        const double clhs154 =             DN(3,2)*clhs129;
        const double clhs155 =             DN(1,2)*clhs65;
        const double clhs156 =             DN(1,0)*N[3];
        const double clhs157 =             DN(1,1)*DN(3,1);
        const double clhs158 =             DN(1,2)*DN(3,2);
        const double clhs159 =             clhs9*(clhs148 + clhs157 + clhs158);
        const double clhs160 =             clhs0*clhs159;
        const double clhs161 =             DN(1,0)*clhs11;
        const double clhs162 =             DN(1,0)*rho;
        const double clhs163 =             DN(1,1)*tau2;
        const double clhs164 =             clhs163 + clhs2*clhs28;
        const double clhs165 =             DN(2,0)*clhs163;
        const double clhs166 =             DN(1,0)*clhs47;
        const double clhs167 =             clhs136 + clhs144*tau2;
        const double clhs168 =             DN(2,2)*clhs163;
        const double clhs169 =             DN(1,1)*N[2];
        const double clhs170 =             clhs146*clhs2;
        const double clhs171 =             DN(3,0)*clhs163;
        const double clhs172 =             DN(1,0)*clhs65;
        const double clhs173 =             clhs149 + clhs157*tau2;
        const double clhs174 =             DN(3,2)*clhs163;
        const double clhs175 =             DN(1,1)*N[3];
        const double clhs176 =             clhs159*clhs2;
        const double clhs177 =             DN(1,2)*tau2;
        const double clhs178 =             clhs177 + clhs28*clhs4;
        const double clhs179 =             DN(2,0)*clhs177;
        const double clhs180 =             DN(2,1)*clhs177;
        const double clhs181 =             clhs136 + clhs145*tau2;
        const double clhs182 =             DN(1,2)*N[2];
        const double clhs183 =             clhs146*clhs4;
        const double clhs184 =             DN(3,0)*clhs177;
        const double clhs185 =             DN(3,1)*clhs177;
        const double clhs186 =             clhs149 + clhs158*tau2;
        const double clhs187 =             DN(1,2)*N[3];
        const double clhs188 =             clhs159*clhs4;
        const double clhs189 =             rho*(N[1] + clhs28);
        const double clhs190 =             DN(2,0)*N[1];
        const double clhs191 =             DN(2,1)*N[1];
        const double clhs192 =             DN(2,2)*N[1];
        const double clhs193 =             DN(3,0)*N[1];
        const double clhs194 =             DN(3,1)*N[1];
        const double clhs195 =             DN(3,2)*N[1];
        const double clhs196 =             N[2]*clhs6;
        const double clhs197 =             DN(2,1)*clhs11;
        const double clhs198 =             DN(2,2)*clhs11;
        const double clhs199 =             N[2]*clhs25;
        const double clhs200 =             DN(2,1)*clhs28;
        const double clhs201 =             DN(2,2)*clhs28;
        const double clhs202 =             pow(N[2], 2)*bdf0 + N[2]*clhs44;
        const double clhs203 =             pow(DN(2,0), 2);
        const double clhs204 =             DN(2,1)*rho;
        const double clhs205 =             DN(2,0)*tau2;
        const double clhs206 =             clhs0*clhs47 + clhs205;
        const double clhs207 =             DN(2,2)*rho;
        const double clhs208 =             pow(DN(2,1), 2);
        const double clhs209 =             pow(DN(2,2), 2);
        const double clhs210 =             clhs9*(clhs203 + clhs208 + clhs209);
        const double clhs211 =             DN(2,0)*DN(3,0);
        const double clhs212 =             N[3]*clhs46;
        const double clhs213 =             clhs211*tau2 + clhs212;
        const double clhs214 =             N[2]*clhs63;
        const double clhs215 =             DN(3,1)*clhs205;
        const double clhs216 =             DN(2,1)*clhs65;
        const double clhs217 =             DN(3,2)*clhs205;
        const double clhs218 =             DN(2,2)*clhs65;
        const double clhs219 =             DN(2,0)*N[3];
        const double clhs220 =             DN(2,1)*DN(3,1);
        const double clhs221 =             DN(2,2)*DN(3,2);
        const double clhs222 =             clhs9*(clhs211 + clhs220 + clhs221);
        const double clhs223 =             clhs0*clhs222;
        const double clhs224 =             DN(2,0)*clhs11;
        const double clhs225 =             DN(2,0)*clhs28;
        const double clhs226 =             DN(2,0)*rho;
        const double clhs227 =             DN(2,1)*tau2;
        const double clhs228 =             clhs2*clhs47 + clhs227;
        const double clhs229 =             DN(3,0)*clhs227;
        const double clhs230 =             DN(2,0)*clhs65;
        const double clhs231 =             clhs212 + clhs220*tau2;
        const double clhs232 =             DN(3,2)*clhs227;
        const double clhs233 =             DN(2,1)*N[3];
        const double clhs234 =             clhs2*clhs222;
        const double clhs235 =             DN(2,2)*tau2;
        const double clhs236 =             clhs235 + clhs4*clhs47;
        const double clhs237 =             DN(3,0)*clhs235;
        const double clhs238 =             DN(3,1)*clhs235;
        const double clhs239 =             clhs212 + clhs221*tau2;
        const double clhs240 =             DN(2,2)*N[3];
        const double clhs241 =             clhs222*clhs4;
        const double clhs242 =             rho*(N[2] + clhs47);
        const double clhs243 =             DN(3,0)*N[2];
        const double clhs244 =             DN(3,1)*N[2];
        const double clhs245 =             DN(3,2)*N[2];
        const double clhs246 =             N[3]*clhs6;
        const double clhs247 =             DN(3,1)*clhs11;
        const double clhs248 =             DN(3,2)*clhs11;
        const double clhs249 =             N[3]*clhs25;
        const double clhs250 =             DN(3,1)*clhs28;
        const double clhs251 =             DN(3,2)*clhs28;
        const double clhs252 =             N[3]*clhs44;
        const double clhs253 =             DN(3,1)*clhs47;
        const double clhs254 =             DN(3,2)*clhs47;
        const double clhs255 =             pow(N[3], 2)*bdf0 + N[3]*clhs63;
        const double clhs256 =             pow(DN(3,0), 2);
        const double clhs257 =             DN(3,1)*rho;
        const double clhs258 =             DN(3,0)*tau2 + clhs0*clhs65;
        const double clhs259 =             DN(3,2)*rho;
        const double clhs260 =             pow(DN(3,1), 2);
        const double clhs261 =             pow(DN(3,2), 2);
        const double clhs262 =             clhs9*(clhs256 + clhs260 + clhs261);
        const double clhs263 =             DN(3,0)*clhs11;
        const double clhs264 =             DN(3,0)*clhs28;
        const double clhs265 =             DN(3,0)*clhs47;
        const double clhs266 =             DN(3,0)*rho;
        const double clhs267 =             DN(3,1)*tau2 + clhs2*clhs65;
        const double clhs268 =             DN(3,2)*tau2 + clhs4*clhs65;
        const double clhs269 =             rho*(N[3] + clhs65);

        lhs(0,0)=rho*(clhs1*clhs11 + clhs7 + clhs8*tau2);
        lhs(0,1)=clhs12*clhs14;
        lhs(0,2)=clhs14*clhs15;
        lhs(0,3)=-DN(0,0)*N[0] + clhs0*clhs18;
        lhs(0,4)=rho*(clhs1*clhs28 + clhs21 + clhs26);
        lhs(0,5)=rho*(clhs0*clhs30 + clhs29);
        lhs(0,6)=rho*(clhs0*clhs32 + clhs31);
        lhs(0,7)=-clhs33 + clhs37;
        lhs(0,8)=rho*(clhs1*clhs47 + clhs40 + clhs45);
        lhs(0,9)=rho*(clhs0*clhs49 + clhs48);
        lhs(0,10)=rho*(clhs0*clhs51 + clhs50);
        lhs(0,11)=-clhs52 + clhs56;
        lhs(0,12)=rho*(clhs1*clhs65 + clhs59 + clhs64);
        lhs(0,13)=rho*(clhs0*clhs67 + clhs66);
        lhs(0,14)=rho*(clhs0*clhs69 + clhs68);
        lhs(0,15)=-clhs70 + clhs74;
        lhs(1,0)=clhs75*clhs77;
        lhs(1,1)=rho*(clhs11*clhs3 + clhs16*tau2 + clhs7);
        lhs(1,2)=clhs15*clhs77;
        lhs(1,3)=-DN(0,1)*N[0] + clhs18*clhs2;
        lhs(1,4)=rho*(clhs2*clhs79 + clhs78);
        lhs(1,5)=rho*(clhs26 + clhs28*clhs3 + clhs80);
        lhs(1,6)=rho*(clhs2*clhs32 + clhs81);
        lhs(1,7)=-clhs82 + clhs83;
        lhs(1,8)=rho*(clhs2*clhs85 + clhs84);
        lhs(1,9)=rho*(clhs3*clhs47 + clhs45 + clhs86);
        lhs(1,10)=rho*(clhs2*clhs51 + clhs87);
        lhs(1,11)=-clhs88 + clhs89;
        lhs(1,12)=rho*(clhs2*clhs91 + clhs90);
        lhs(1,13)=rho*(clhs3*clhs65 + clhs64 + clhs92);
        lhs(1,14)=rho*(clhs2*clhs69 + clhs93);
        lhs(1,15)=-clhs94 + clhs95;
        lhs(2,0)=clhs75*clhs97;
        lhs(2,1)=clhs12*clhs97;
        lhs(2,2)=rho*(clhs11*clhs5 + clhs17*tau2 + clhs7);
        lhs(2,3)=-DN(0,2)*N[0] + clhs18*clhs4;
        lhs(2,4)=rho*(clhs4*clhs79 + clhs98);
        lhs(2,5)=rho*(clhs30*clhs4 + clhs99);
        lhs(2,6)=rho*(clhs100 + clhs26 + clhs28*clhs5);
        lhs(2,7)=-clhs101 + clhs102;
        lhs(2,8)=rho*(clhs103 + clhs4*clhs85);
        lhs(2,9)=rho*(clhs104 + clhs4*clhs49);
        lhs(2,10)=rho*(clhs105 + clhs45 + clhs47*clhs5);
        lhs(2,11)=-clhs106 + clhs107;
        lhs(2,12)=rho*(clhs108 + clhs4*clhs91);
        lhs(2,13)=rho*(clhs109 + clhs4*clhs67);
        lhs(2,14)=rho*(clhs110 + clhs5*clhs65 + clhs64);
        lhs(2,15)=-clhs111 + clhs112;
        lhs(3,0)=DN(0,0)*clhs113;
        lhs(3,1)=DN(0,1)*clhs113;
        lhs(3,2)=DN(0,2)*clhs113;
        lhs(3,3)=clhs18;
        lhs(3,4)=rho*(clhs114 + clhs79);
        lhs(3,5)=rho*(clhs115 + clhs30);
        lhs(3,6)=rho*(clhs116 + clhs32);
        lhs(3,7)=clhs36;
        lhs(3,8)=rho*(clhs117 + clhs85);
        lhs(3,9)=rho*(clhs118 + clhs49);
        lhs(3,10)=rho*(clhs119 + clhs51);
        lhs(3,11)=clhs55;
        lhs(3,12)=rho*(clhs120 + clhs91);
        lhs(3,13)=rho*(clhs121 + clhs67);
        lhs(3,14)=rho*(clhs122 + clhs69);
        lhs(3,15)=clhs73;
        lhs(4,0)=rho*(clhs11*clhs22 + clhs123 + clhs21);
        lhs(4,1)=rho*(clhs0*clhs124 + clhs78);
        lhs(4,2)=rho*(clhs0*clhs125 + clhs98);
        lhs(4,3)=-clhs114 + clhs37;
        lhs(4,4)=rho*(clhs126 + clhs127*tau2 + clhs22*clhs28);
        lhs(4,5)=clhs128*clhs130;
        lhs(4,6)=clhs130*clhs131;
        lhs(4,7)=-DN(1,0)*N[1] + clhs0*clhs134;
        lhs(4,8)=rho*(clhs137 + clhs138 + clhs22*clhs47);
        lhs(4,9)=rho*(clhs0*clhs140 + clhs139);
        lhs(4,10)=rho*(clhs0*clhs142 + clhs141);
        lhs(4,11)=-clhs143 + clhs147;
        lhs(4,12)=rho*(clhs150 + clhs151 + clhs22*clhs65);
        lhs(4,13)=rho*(clhs0*clhs153 + clhs152);
        lhs(4,14)=rho*(clhs0*clhs155 + clhs154);
        lhs(4,15)=-clhs156 + clhs160;
        lhs(5,0)=rho*(clhs161*clhs2 + clhs29);
        lhs(5,1)=rho*(clhs11*clhs23 + clhs123 + clhs80);
        lhs(5,2)=rho*(clhs125*clhs2 + clhs99);
        lhs(5,3)=-clhs115 + clhs83;
        lhs(5,4)=clhs162*clhs164;
        lhs(5,5)=rho*(clhs126 + clhs132*tau2 + clhs23*clhs28);
        lhs(5,6)=clhs131*clhs164;
        lhs(5,7)=-DN(1,1)*N[1] + clhs134*clhs2;
        lhs(5,8)=rho*(clhs165 + clhs166*clhs2);
        lhs(5,9)=rho*(clhs138 + clhs167 + clhs23*clhs47);
        lhs(5,10)=rho*(clhs142*clhs2 + clhs168);
        lhs(5,11)=-clhs169 + clhs170;
        lhs(5,12)=rho*(clhs171 + clhs172*clhs2);
        lhs(5,13)=rho*(clhs151 + clhs173 + clhs23*clhs65);
        lhs(5,14)=rho*(clhs155*clhs2 + clhs174);
        lhs(5,15)=-clhs175 + clhs176;
        lhs(6,0)=rho*(clhs161*clhs4 + clhs31);
        lhs(6,1)=rho*(clhs124*clhs4 + clhs81);
        lhs(6,2)=rho*(clhs100 + clhs11*clhs24 + clhs123);
        lhs(6,3)=clhs102 - clhs116;
        lhs(6,4)=clhs162*clhs178;
        lhs(6,5)=clhs128*clhs178;
        lhs(6,6)=rho*(clhs126 + clhs133*tau2 + clhs24*clhs28);
        lhs(6,7)=-DN(1,2)*N[1] + clhs134*clhs4;
        lhs(6,8)=rho*(clhs166*clhs4 + clhs179);
        lhs(6,9)=rho*(clhs140*clhs4 + clhs180);
        lhs(6,10)=rho*(clhs138 + clhs181 + clhs24*clhs47);
        lhs(6,11)=-clhs182 + clhs183;
        lhs(6,12)=rho*(clhs172*clhs4 + clhs184);
        lhs(6,13)=rho*(clhs153*clhs4 + clhs185);
        lhs(6,14)=rho*(clhs151 + clhs186 + clhs24*clhs65);
        lhs(6,15)=-clhs187 + clhs188;
        lhs(7,0)=rho*(clhs161 + clhs33);
        lhs(7,1)=rho*(clhs124 + clhs82);
        lhs(7,2)=rho*(clhs101 + clhs125);
        lhs(7,3)=clhs36;
        lhs(7,4)=DN(1,0)*clhs189;
        lhs(7,5)=DN(1,1)*clhs189;
        lhs(7,6)=DN(1,2)*clhs189;
        lhs(7,7)=clhs134;
        lhs(7,8)=rho*(clhs166 + clhs190);
        lhs(7,9)=rho*(clhs140 + clhs191);
        lhs(7,10)=rho*(clhs142 + clhs192);
        lhs(7,11)=clhs146;
        lhs(7,12)=rho*(clhs172 + clhs193);
        lhs(7,13)=rho*(clhs153 + clhs194);
        lhs(7,14)=rho*(clhs155 + clhs195);
        lhs(7,15)=clhs159;
        lhs(8,0)=rho*(clhs11*clhs41 + clhs196 + clhs40);
        lhs(8,1)=rho*(clhs0*clhs197 + clhs84);
        lhs(8,2)=rho*(clhs0*clhs198 + clhs103);
        lhs(8,3)=-clhs117 + clhs56;
        lhs(8,4)=rho*(clhs137 + clhs199 + clhs28*clhs41);
        lhs(8,5)=rho*(clhs0*clhs200 + clhs165);
        lhs(8,6)=rho*(clhs0*clhs201 + clhs179);
        lhs(8,7)=clhs147 - clhs190;
        lhs(8,8)=rho*(clhs202 + clhs203*tau2 + clhs41*clhs47);
        lhs(8,9)=clhs204*clhs206;
        lhs(8,10)=clhs206*clhs207;
        lhs(8,11)=-DN(2,0)*N[2] + clhs0*clhs210;
        lhs(8,12)=rho*(clhs213 + clhs214 + clhs41*clhs65);
        lhs(8,13)=rho*(clhs0*clhs216 + clhs215);
        lhs(8,14)=rho*(clhs0*clhs218 + clhs217);
        lhs(8,15)=-clhs219 + clhs223;
        lhs(9,0)=rho*(clhs2*clhs224 + clhs48);
        lhs(9,1)=rho*(clhs11*clhs42 + clhs196 + clhs86);
        lhs(9,2)=rho*(clhs104 + clhs198*clhs2);
        lhs(9,3)=-clhs118 + clhs89;
        lhs(9,4)=rho*(clhs139 + clhs2*clhs225);
        lhs(9,5)=rho*(clhs167 + clhs199 + clhs28*clhs42);
        lhs(9,6)=rho*(clhs180 + clhs2*clhs201);
        lhs(9,7)=clhs170 - clhs191;
        lhs(9,8)=clhs226*clhs228;
        lhs(9,9)=rho*(clhs202 + clhs208*tau2 + clhs42*clhs47);
        lhs(9,10)=clhs207*clhs228;
        lhs(9,11)=-DN(2,1)*N[2] + clhs2*clhs210;
        lhs(9,12)=rho*(clhs2*clhs230 + clhs229);
        lhs(9,13)=rho*(clhs214 + clhs231 + clhs42*clhs65);
        lhs(9,14)=rho*(clhs2*clhs218 + clhs232);
        lhs(9,15)=-clhs233 + clhs234;
        lhs(10,0)=rho*(clhs224*clhs4 + clhs50);
        lhs(10,1)=rho*(clhs197*clhs4 + clhs87);
        lhs(10,2)=rho*(clhs105 + clhs11*clhs43 + clhs196);
        lhs(10,3)=clhs107 - clhs119;
        lhs(10,4)=rho*(clhs141 + clhs225*clhs4);
        lhs(10,5)=rho*(clhs168 + clhs200*clhs4);
        lhs(10,6)=rho*(clhs181 + clhs199 + clhs28*clhs43);
        lhs(10,7)=clhs183 - clhs192;
        lhs(10,8)=clhs226*clhs236;
        lhs(10,9)=clhs204*clhs236;
        lhs(10,10)=rho*(clhs202 + clhs209*tau2 + clhs43*clhs47);
        lhs(10,11)=-DN(2,2)*N[2] + clhs210*clhs4;
        lhs(10,12)=rho*(clhs230*clhs4 + clhs237);
        lhs(10,13)=rho*(clhs216*clhs4 + clhs238);
        lhs(10,14)=rho*(clhs214 + clhs239 + clhs43*clhs65);
        lhs(10,15)=-clhs240 + clhs241;
        lhs(11,0)=rho*(clhs224 + clhs52);
        lhs(11,1)=rho*(clhs197 + clhs88);
        lhs(11,2)=rho*(clhs106 + clhs198);
        lhs(11,3)=clhs55;
        lhs(11,4)=rho*(clhs143 + clhs225);
        lhs(11,5)=rho*(clhs169 + clhs200);
        lhs(11,6)=rho*(clhs182 + clhs201);
        lhs(11,7)=clhs146;
        lhs(11,8)=DN(2,0)*clhs242;
        lhs(11,9)=DN(2,1)*clhs242;
        lhs(11,10)=DN(2,2)*clhs242;
        lhs(11,11)=clhs210;
        lhs(11,12)=rho*(clhs230 + clhs243);
        lhs(11,13)=rho*(clhs216 + clhs244);
        lhs(11,14)=rho*(clhs218 + clhs245);
        lhs(11,15)=clhs222;
        lhs(12,0)=rho*(clhs11*clhs60 + clhs246 + clhs59);
        lhs(12,1)=rho*(clhs0*clhs247 + clhs90);
        lhs(12,2)=rho*(clhs0*clhs248 + clhs108);
        lhs(12,3)=-clhs120 + clhs74;
        lhs(12,4)=rho*(clhs150 + clhs249 + clhs28*clhs60);
        lhs(12,5)=rho*(clhs0*clhs250 + clhs171);
        lhs(12,6)=rho*(clhs0*clhs251 + clhs184);
        lhs(12,7)=clhs160 - clhs193;
        lhs(12,8)=rho*(clhs213 + clhs252 + clhs47*clhs60);
        lhs(12,9)=rho*(clhs0*clhs253 + clhs229);
        lhs(12,10)=rho*(clhs0*clhs254 + clhs237);
        lhs(12,11)=clhs223 - clhs243;
        lhs(12,12)=rho*(clhs255 + clhs256*tau2 + clhs60*clhs65);
        lhs(12,13)=clhs257*clhs258;
        lhs(12,14)=clhs258*clhs259;
        lhs(12,15)=-DN(3,0)*N[3] + clhs0*clhs262;
        lhs(13,0)=rho*(clhs2*clhs263 + clhs66);
        lhs(13,1)=rho*(clhs11*clhs61 + clhs246 + clhs92);
        lhs(13,2)=rho*(clhs109 + clhs2*clhs248);
        lhs(13,3)=-clhs121 + clhs95;
        lhs(13,4)=rho*(clhs152 + clhs2*clhs264);
        lhs(13,5)=rho*(clhs173 + clhs249 + clhs28*clhs61);
        lhs(13,6)=rho*(clhs185 + clhs2*clhs251);
        lhs(13,7)=clhs176 - clhs194;
        lhs(13,8)=rho*(clhs2*clhs265 + clhs215);
        lhs(13,9)=rho*(clhs231 + clhs252 + clhs47*clhs61);
        lhs(13,10)=rho*(clhs2*clhs254 + clhs238);
        lhs(13,11)=clhs234 - clhs244;
        lhs(13,12)=clhs266*clhs267;
        lhs(13,13)=rho*(clhs255 + clhs260*tau2 + clhs61*clhs65);
        lhs(13,14)=clhs259*clhs267;
        lhs(13,15)=-DN(3,1)*N[3] + clhs2*clhs262;
        lhs(14,0)=rho*(clhs263*clhs4 + clhs68);
        lhs(14,1)=rho*(clhs247*clhs4 + clhs93);
        lhs(14,2)=rho*(clhs11*clhs62 + clhs110 + clhs246);
        lhs(14,3)=clhs112 - clhs122;
        lhs(14,4)=rho*(clhs154 + clhs264*clhs4);
        lhs(14,5)=rho*(clhs174 + clhs250*clhs4);
        lhs(14,6)=rho*(clhs186 + clhs249 + clhs28*clhs62);
        lhs(14,7)=clhs188 - clhs195;
        lhs(14,8)=rho*(clhs217 + clhs265*clhs4);
        lhs(14,9)=rho*(clhs232 + clhs253*clhs4);
        lhs(14,10)=rho*(clhs239 + clhs252 + clhs47*clhs62);
        lhs(14,11)=clhs241 - clhs245;
        lhs(14,12)=clhs266*clhs268;
        lhs(14,13)=clhs257*clhs268;
        lhs(14,14)=rho*(clhs255 + clhs261*tau2 + clhs62*clhs65);
        lhs(14,15)=-DN(3,2)*N[3] + clhs262*clhs4;
        lhs(15,0)=rho*(clhs263 + clhs70);
        lhs(15,1)=rho*(clhs247 + clhs94);
        lhs(15,2)=rho*(clhs111 + clhs248);
        lhs(15,3)=clhs73;
        lhs(15,4)=rho*(clhs156 + clhs264);
        lhs(15,5)=rho*(clhs175 + clhs250);
        lhs(15,6)=rho*(clhs187 + clhs251);
        lhs(15,7)=clhs159;
        lhs(15,8)=rho*(clhs219 + clhs265);
        lhs(15,9)=rho*(clhs233 + clhs253);
        lhs(15,10)=rho*(clhs240 + clhs254);
        lhs(15,11)=clhs222;
        lhs(15,12)=DN(3,0)*clhs269;
        lhs(15,13)=DN(3,1)*clhs269;
        lhs(15,14)=DN(3,2)*clhs269;
        lhs(15,15)=clhs262;

    }


template<>
void NavierStokes<2>::ComputeGaussPointLHSContribution(bounded_matrix<double,9,9>& lhs, const element_data& data)
    {
        const int nnodes = 3;
        const int dim = 2;
        const int strain_size = 3;
        
        const double rho = inner_prod(data.N, data.rho);        // Density
        const double nu = inner_prod(data.N, data.nu);          // Kinematic viscosity
        const double h = data.h;                                // Characteristic element size
        
        const double& bdf0 = data.bdf0;
        //~ const double& bdf1 = data.bdf1;
        //~ const double& bdf2 = data.bdf2;
        const double& delta_t = data.delta_t;
        const double& dyn_tau_coeff = data.dyn_tau_coeff;
        
        const bounded_matrix<double,nnodes,dim>& v = data.v;
        //~ const bounded_matrix<double,nnodes,dim>& vn = data.vn;
        //~ const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
        const bounded_matrix<double,nnodes,dim>& vmesh = data.vmesh;
        const bounded_matrix<double,nnodes,dim>& vconv = v - vmesh;
        //~ const bounded_matrix<double,nnodes,dim>& f = data.f;
        //~ const array_1d<double,nnodes>& p = data.p;
        const array_1d<double,strain_size>& stress = data.stress;
        
        // Get constitutive matrix 
        //~ const Matrix& C = data.C;
        
        // Get shape function values
        const array_1d<double,nnodes>& N = data.N;
        const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;
        
        const array_1d<double,dim> vconv_gauss = prod(trans(vconv), N);
        
        const double vconv_norm = norm_2(vconv_gauss);
                
        // Stabilization parameters
        const double tau1 = 1.0/((rho*dyn_tau_coeff)/delta_t + (2*rho*vconv_norm)/h + (4*rho*nu)/(h*h));
        const double tau2 = (rho*nu) + 0.5*h*vconv_norm;
        
        const double clhs0 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
        const double clhs1 =             DN(0,0)*clhs0;
        const double clhs2 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
        const double clhs3 =             DN(0,1)*clhs2;
        const double clhs4 =             clhs1 + clhs3;
        const double clhs5 =             pow(N[0], 2)*bdf0 + N[0]*clhs4;
        const double clhs6 =             pow(DN(0,0), 2);
        const double clhs7 =             rho*tau1;
        const double clhs8 =             N[0]*bdf0;
        const double clhs9 =             clhs7*(clhs4 + clhs8);
        const double clhs10 =             DN(0,0)*tau2;
        const double clhs11 =             pow(DN(0,1), 2);
        const double clhs12 =             clhs7*(clhs11 + clhs6);
        const double clhs13 =             DN(0,0)*DN(1,0);
        const double clhs14 =             N[1]*clhs8;
        const double clhs15 =             clhs13*tau2 + clhs14;
        const double clhs16 =             DN(1,0)*clhs0;
        const double clhs17 =             DN(1,1)*clhs2;
        const double clhs18 =             clhs16 + clhs17;
        const double clhs19 =             N[0]*clhs18;
        const double clhs20 =             N[1]*bdf0;
        const double clhs21 =             clhs7*(clhs18 + clhs20);
        const double clhs22 =             DN(1,1)*clhs10;
        const double clhs23 =             DN(0,1)*clhs21;
        const double clhs24 =             DN(0,0)*N[1];
        const double clhs25 =             DN(0,1)*DN(1,1);
        const double clhs26 =             clhs7*(clhs13 + clhs25);
        const double clhs27 =             clhs0*clhs26;
        const double clhs28 =             DN(0,0)*DN(2,0);
        const double clhs29 =             N[2]*clhs8;
        const double clhs30 =             clhs28*tau2 + clhs29;
        const double clhs31 =             DN(2,0)*clhs0;
        const double clhs32 =             DN(2,1)*clhs2;
        const double clhs33 =             clhs31 + clhs32;
        const double clhs34 =             N[0]*clhs33;
        const double clhs35 =             clhs7*(N[2]*bdf0 + clhs33);
        const double clhs36 =             DN(2,1)*clhs10;
        const double clhs37 =             DN(0,1)*clhs35;
        const double clhs38 =             DN(0,0)*N[2];
        const double clhs39 =             DN(0,1)*DN(2,1);
        const double clhs40 =             clhs7*(clhs28 + clhs39);
        const double clhs41 =             clhs0*clhs40;
        const double clhs42 =             DN(0,1)*tau2;
        const double clhs43 =             DN(1,0)*clhs42;
        const double clhs44 =             DN(0,0)*clhs21;
        const double clhs45 =             clhs14 + clhs25*tau2;
        const double clhs46 =             DN(0,1)*N[1];
        const double clhs47 =             clhs2*clhs26;
        const double clhs48 =             DN(2,0)*clhs42;
        const double clhs49 =             DN(0,0)*clhs35;
        const double clhs50 =             clhs29 + clhs39*tau2;
        const double clhs51 =             DN(0,1)*N[2];
        const double clhs52 =             clhs2*clhs40;
        const double clhs53 =             rho*(N[0] + clhs9);
        const double clhs54 =             DN(1,0)*N[0];
        const double clhs55 =             DN(1,1)*N[0];
        const double clhs56 =             DN(2,0)*N[0];
        const double clhs57 =             DN(2,1)*N[0];
        const double clhs58 =             N[1]*clhs4;
        const double clhs59 =             DN(1,1)*clhs9;
        const double clhs60 =             pow(N[1], 2)*bdf0 + N[1]*clhs18;
        const double clhs61 =             pow(DN(1,0), 2);
        const double clhs62 =             DN(1,0)*tau2;
        const double clhs63 =             pow(DN(1,1), 2);
        const double clhs64 =             clhs7*(clhs61 + clhs63);
        const double clhs65 =             DN(1,0)*DN(2,0);
        const double clhs66 =             N[2]*clhs20;
        const double clhs67 =             clhs65*tau2 + clhs66;
        const double clhs68 =             N[1]*clhs33;
        const double clhs69 =             DN(2,1)*clhs62;
        const double clhs70 =             DN(1,1)*clhs35;
        const double clhs71 =             DN(1,0)*N[2];
        const double clhs72 =             DN(1,1)*DN(2,1);
        const double clhs73 =             clhs7*(clhs65 + clhs72);
        const double clhs74 =             clhs0*clhs73;
        const double clhs75 =             DN(1,0)*clhs9;
        const double clhs76 =             DN(1,1)*tau2;
        const double clhs77 =             DN(2,0)*clhs76;
        const double clhs78 =             DN(1,0)*clhs35;
        const double clhs79 =             clhs66 + clhs72*tau2;
        const double clhs80 =             DN(1,1)*N[2];
        const double clhs81 =             clhs2*clhs73;
        const double clhs82 =             rho*(N[1] + clhs21);
        const double clhs83 =             DN(2,0)*N[1];
        const double clhs84 =             DN(2,1)*N[1];
        const double clhs85 =             N[2]*clhs4;
        const double clhs86 =             DN(2,1)*clhs9;
        const double clhs87 =             N[2]*clhs18;
        const double clhs88 =             DN(2,1)*clhs21;
        const double clhs89 =             pow(N[2], 2)*bdf0 + N[2]*clhs33;
        const double clhs90 =             pow(DN(2,0), 2);
        const double clhs91 =             pow(DN(2,1), 2);
        const double clhs92 =             clhs7*(clhs90 + clhs91);
        const double clhs93 =             DN(2,0)*clhs9;
        const double clhs94 =             DN(2,0)*clhs21;
        const double clhs95 =             rho*(N[2] + clhs35);

        lhs(0,0)=rho*(clhs1*clhs9 + clhs5 + clhs6*tau2);
        lhs(0,1)=DN(0,1)*rho*(clhs0*clhs9 + clhs10);
        lhs(0,2)=-DN(0,0)*N[0] + clhs0*clhs12;
        lhs(0,3)=rho*(clhs1*clhs21 + clhs15 + clhs19);
        lhs(0,4)=rho*(clhs0*clhs23 + clhs22);
        lhs(0,5)=-clhs24 + clhs27;
        lhs(0,6)=rho*(clhs1*clhs35 + clhs30 + clhs34);
        lhs(0,7)=rho*(clhs0*clhs37 + clhs36);
        lhs(0,8)=-clhs38 + clhs41;
        lhs(1,0)=DN(0,0)*rho*(clhs2*clhs9 + clhs42);
        lhs(1,1)=rho*(clhs11*tau2 + clhs3*clhs9 + clhs5);
        lhs(1,2)=-DN(0,1)*N[0] + clhs12*clhs2;
        lhs(1,3)=rho*(clhs2*clhs44 + clhs43);
        lhs(1,4)=rho*(clhs19 + clhs21*clhs3 + clhs45);
        lhs(1,5)=-clhs46 + clhs47;
        lhs(1,6)=rho*(clhs2*clhs49 + clhs48);
        lhs(1,7)=rho*(clhs3*clhs35 + clhs34 + clhs50);
        lhs(1,8)=-clhs51 + clhs52;
        lhs(2,0)=DN(0,0)*clhs53;
        lhs(2,1)=DN(0,1)*clhs53;
        lhs(2,2)=clhs12;
        lhs(2,3)=rho*(clhs44 + clhs54);
        lhs(2,4)=rho*(clhs23 + clhs55);
        lhs(2,5)=clhs26;
        lhs(2,6)=rho*(clhs49 + clhs56);
        lhs(2,7)=rho*(clhs37 + clhs57);
        lhs(2,8)=clhs40;
        lhs(3,0)=rho*(clhs15 + clhs16*clhs9 + clhs58);
        lhs(3,1)=rho*(clhs0*clhs59 + clhs43);
        lhs(3,2)=clhs27 - clhs54;
        lhs(3,3)=rho*(clhs16*clhs21 + clhs60 + clhs61*tau2);
        lhs(3,4)=DN(1,1)*rho*(clhs0*clhs21 + clhs62);
        lhs(3,5)=-DN(1,0)*N[1] + clhs0*clhs64;
        lhs(3,6)=rho*(clhs16*clhs35 + clhs67 + clhs68);
        lhs(3,7)=rho*(clhs0*clhs70 + clhs69);
        lhs(3,8)=-clhs71 + clhs74;
        lhs(4,0)=rho*(clhs2*clhs75 + clhs22);
        lhs(4,1)=rho*(clhs17*clhs9 + clhs45 + clhs58);
        lhs(4,2)=clhs47 - clhs55;
        lhs(4,3)=DN(1,0)*rho*(clhs2*clhs21 + clhs76);
        lhs(4,4)=rho*(clhs17*clhs21 + clhs60 + clhs63*tau2);
        lhs(4,5)=-DN(1,1)*N[1] + clhs2*clhs64;
        lhs(4,6)=rho*(clhs2*clhs78 + clhs77);
        lhs(4,7)=rho*(clhs17*clhs35 + clhs68 + clhs79);
        lhs(4,8)=-clhs80 + clhs81;
        lhs(5,0)=rho*(clhs24 + clhs75);
        lhs(5,1)=rho*(clhs46 + clhs59);
        lhs(5,2)=clhs26;
        lhs(5,3)=DN(1,0)*clhs82;
        lhs(5,4)=DN(1,1)*clhs82;
        lhs(5,5)=clhs64;
        lhs(5,6)=rho*(clhs78 + clhs83);
        lhs(5,7)=rho*(clhs70 + clhs84);
        lhs(5,8)=clhs73;
        lhs(6,0)=rho*(clhs30 + clhs31*clhs9 + clhs85);
        lhs(6,1)=rho*(clhs0*clhs86 + clhs48);
        lhs(6,2)=clhs41 - clhs56;
        lhs(6,3)=rho*(clhs21*clhs31 + clhs67 + clhs87);
        lhs(6,4)=rho*(clhs0*clhs88 + clhs77);
        lhs(6,5)=clhs74 - clhs83;
        lhs(6,6)=rho*(clhs31*clhs35 + clhs89 + clhs90*tau2);
        lhs(6,7)=DN(2,1)*rho*(DN(2,0)*tau2 + clhs0*clhs35);
        lhs(6,8)=-DN(2,0)*N[2] + clhs0*clhs92;
        lhs(7,0)=rho*(clhs2*clhs93 + clhs36);
        lhs(7,1)=rho*(clhs32*clhs9 + clhs50 + clhs85);
        lhs(7,2)=clhs52 - clhs57;
        lhs(7,3)=rho*(clhs2*clhs94 + clhs69);
        lhs(7,4)=rho*(clhs21*clhs32 + clhs79 + clhs87);
        lhs(7,5)=clhs81 - clhs84;
        lhs(7,6)=DN(2,0)*rho*(DN(2,1)*tau2 + clhs2*clhs35);
        lhs(7,7)=rho*(clhs32*clhs35 + clhs89 + clhs91*tau2);
        lhs(7,8)=-DN(2,1)*N[2] + clhs2*clhs92;
        lhs(8,0)=rho*(clhs38 + clhs93);
        lhs(8,1)=rho*(clhs51 + clhs86);
        lhs(8,2)=clhs40;
        lhs(8,3)=rho*(clhs71 + clhs94);
        lhs(8,4)=rho*(clhs80 + clhs88);
        lhs(8,5)=clhs73;
        lhs(8,6)=DN(2,0)*clhs95;
        lhs(8,7)=DN(2,1)*clhs95;
        lhs(8,8)=clhs92;

    }


template<>
void NavierStokes<3>::ComputeGaussPointRHSContribution(array_1d<double,16>& rhs, const element_data& data)
    {
        const int nnodes = 4;
        const int dim = 3;
        const int strain_size = 6;
        
        const double rho = inner_prod(data.N, data.rho);        // Density
        const double nu = inner_prod(data.N, data.nu);          // Kinematic viscosity
        const double h = data.h;                                // Characteristic element size
        
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
        const array_1d<double,strain_size>& stress = data.stress;
        
        // Get constitutive matrix 
        //~ const Matrix& C = data.C;
        
        // Get shape function values
        const array_1d<double,nnodes>& N = data.N;
        const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;
        
        const array_1d<double,dim> vconv_gauss = prod(trans(vconv), N);
        
        const double vconv_norm = norm_2(vconv_gauss);
                
        // Stabilization parameters
        const double tau1 = 1.0/((rho*dyn_tau_coeff)/delta_t + (2*rho*vconv_norm)/h + (4*rho*nu)/(h*h));
        const double tau2 = (rho*nu) + 0.5*h*vconv_norm;
        
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
        const double crhs12 =             crhs11*rho*tau2;
        const double crhs13 =             N[0]*rho;
        const double crhs14 =             N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0)) + N[3]*(bdf0*v(3,0) + bdf1*vn(3,0) + bdf2*vnn(3,0));
        const double crhs15 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0) + N[3]*vconv(3,0);
        const double crhs16 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1) + N[3]*vconv(3,1);
        const double crhs17 =             N[0]*vconv(0,2) + N[1]*vconv(1,2) + N[2]*vconv(2,2) + N[3]*vconv(3,2);
        const double crhs18 =             crhs15*(crhs3 + crhs5 + crhs7 + crhs9) + crhs16*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0) + DN(3,1)*v(3,0)) + crhs17*(DN(0,2)*v(0,0) + DN(1,2)*v(1,0) + DN(2,2)*v(2,0) + DN(3,2)*v(3,0));
        const double crhs19 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] + DN(3,0)*p[3] - crhs1 + rho*(crhs14 + crhs18);
        const double crhs20 =             DN(0,0)*crhs19;
        const double crhs21 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1) + N[3]*f(3,1));
        const double crhs22 =             N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1)) + N[3]*(bdf0*v(3,1) + bdf1*vn(3,1) + bdf2*vnn(3,1));
        const double crhs23 =             crhs15*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1) + DN(3,0)*v(3,1)) + crhs16*(crhs10 + crhs4 + crhs6 + crhs8) + crhs17*(DN(0,2)*v(0,1) + DN(1,2)*v(1,1) + DN(2,2)*v(2,1) + DN(3,2)*v(3,1));
        const double crhs24 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] + DN(3,1)*p[3] - crhs21 + rho*(crhs22 + crhs23);
        const double crhs25 =             DN(0,1)*crhs24;
        const double crhs26 =             rho*(N[0]*f(0,2) + N[1]*f(1,2) + N[2]*f(2,2) + N[3]*f(3,2));
        const double crhs27 =             N[0]*(bdf0*v(0,2) + bdf1*vn(0,2) + bdf2*vnn(0,2)) + N[1]*(bdf0*v(1,2) + bdf1*vn(1,2) + bdf2*vnn(1,2)) + N[2]*(bdf0*v(2,2) + bdf1*vn(2,2) + bdf2*vnn(2,2)) + N[3]*(bdf0*v(3,2) + bdf1*vn(3,2) + bdf2*vnn(3,2));
        const double crhs28 =             crhs15*(DN(0,0)*v(0,2) + DN(1,0)*v(1,2) + DN(2,0)*v(2,2) + DN(3,0)*v(3,2)) + crhs16*(DN(0,1)*v(0,2) + DN(1,1)*v(1,2) + DN(2,1)*v(2,2) + DN(3,1)*v(3,2)) + crhs17*crhs2;
        const double crhs29 =             DN(0,2)*p[0] + DN(1,2)*p[1] + DN(2,2)*p[2] + DN(3,2)*p[3] - crhs26 + rho*(crhs27 + crhs28);
        const double crhs30 =             DN(0,2)*crhs29;
        const double crhs31 =             rho*tau1*(crhs20 + crhs25 + crhs30);
        const double crhs32 =             N[1]*rho;
        const double crhs33 =             DN(1,0)*crhs19;
        const double crhs34 =             DN(1,1)*crhs24;
        const double crhs35 =             DN(1,2)*crhs29;
        const double crhs36 =             rho*tau1*(crhs33 + crhs34 + crhs35);
        const double crhs37 =             N[2]*rho;
        const double crhs38 =             DN(2,0)*crhs19;
        const double crhs39 =             DN(2,1)*crhs24;
        const double crhs40 =             DN(2,2)*crhs29;
        const double crhs41 =             rho*tau1*(crhs38 + crhs39 + crhs40);
        const double crhs42 =             N[3]*rho;
        const double crhs43 =             DN(3,0)*crhs19;
        const double crhs44 =             DN(3,1)*crhs24;
        const double crhs45 =             DN(3,2)*crhs29;
        const double crhs46 =             rho*tau1*(crhs43 + crhs44 + crhs45);

        rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs12 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs1 - crhs13*crhs14 - crhs13*crhs18 - crhs15*crhs31;
        rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs12 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs21 - crhs13*crhs22 - crhs13*crhs23 - crhs16*crhs31;
        rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs12 - DN(0,2)*stress[2] + N[0]*crhs26 - crhs13*crhs27 - crhs13*crhs28 - crhs17*crhs31;
        rhs[3]=-rho*(N[0]*crhs11 + crhs20*tau1 + crhs25*tau1 + crhs30*tau1);
        rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs12 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs1 - crhs14*crhs32 - crhs15*crhs36 - crhs18*crhs32;
        rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs12 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs21 - crhs16*crhs36 - crhs22*crhs32 - crhs23*crhs32;
        rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs12 - DN(1,2)*stress[2] + N[1]*crhs26 - crhs17*crhs36 - crhs27*crhs32 - crhs28*crhs32;
        rhs[7]=-rho*(N[1]*crhs11 + crhs33*tau1 + crhs34*tau1 + crhs35*tau1);
        rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs12 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs1 - crhs14*crhs37 - crhs15*crhs41 - crhs18*crhs37;
        rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs12 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs21 - crhs16*crhs41 - crhs22*crhs37 - crhs23*crhs37;
        rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs12 - DN(2,2)*stress[2] + N[2]*crhs26 - crhs17*crhs41 - crhs27*crhs37 - crhs28*crhs37;
        rhs[11]=-rho*(N[2]*crhs11 + crhs38*tau1 + crhs39*tau1 + crhs40*tau1);
        rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs12 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs1 - crhs14*crhs42 - crhs15*crhs46 - crhs18*crhs42;
        rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs12 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs21 - crhs16*crhs46 - crhs22*crhs42 - crhs23*crhs42;
        rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs12 - DN(3,2)*stress[2] + N[3]*crhs26 - crhs17*crhs46 - crhs27*crhs42 - crhs28*crhs42;
        rhs[15]=-rho*(N[3]*crhs11 + crhs43*tau1 + crhs44*tau1 + crhs45*tau1);

    }


template<>
void NavierStokes<2>::ComputeGaussPointRHSContribution(array_1d<double,9>& rhs, const element_data& data)
    {
        const int nnodes = 3;
        const int dim = 2;
        const int strain_size = 3;
        
        const double rho = inner_prod(data.N, data.rho);        // Density
        const double nu = inner_prod(data.N, data.nu);          // Kinematic viscosity
        const double h = data.h;                                // Characteristic element size
        
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
        const array_1d<double,strain_size>& stress = data.stress;
        
        // Get constitutive matrix 
        //~ const Matrix& C = data.C;
        
        // Get shape function values
        const array_1d<double,nnodes>& N = data.N;
        const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;
        
        const array_1d<double,dim> vconv_gauss = prod(trans(vconv), N);
        
        const double vconv_norm = norm_2(vconv_gauss);
                
        // Stabilization parameters
        const double tau1 = 1.0/((rho*dyn_tau_coeff)/delta_t + (2*rho*vconv_norm)/h + (4*rho*nu)/(h*h));
        const double tau2 = (rho*nu) + 0.5*h*vconv_norm;
        
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
        const double crhs5 =             crhs4*rho*tau2;
        const double crhs6 =             N[0]*rho;
        const double crhs7 =             N[0]*(bdf0*v(0,0) + bdf1*vn(0,0) + bdf2*vnn(0,0)) + N[1]*(bdf0*v(1,0) + bdf1*vn(1,0) + bdf2*vnn(1,0)) + N[2]*(bdf0*v(2,0) + bdf1*vn(2,0) + bdf2*vnn(2,0));
        const double crhs8 =             N[0]*vconv(0,0) + N[1]*vconv(1,0) + N[2]*vconv(2,0);
        const double crhs9 =             N[0]*vconv(0,1) + N[1]*vconv(1,1) + N[2]*vconv(2,1);
        const double crhs10 =             crhs2*crhs8 + crhs9*(DN(0,1)*v(0,0) + DN(1,1)*v(1,0) + DN(2,1)*v(2,0));
        const double crhs11 =             DN(0,0)*p[0] + DN(1,0)*p[1] + DN(2,0)*p[2] - crhs1 + rho*(crhs10 + crhs7);
        const double crhs12 =             DN(0,0)*crhs11;
        const double crhs13 =             rho*(N[0]*f(0,1) + N[1]*f(1,1) + N[2]*f(2,1));
        const double crhs14 =             N[0]*(bdf0*v(0,1) + bdf1*vn(0,1) + bdf2*vnn(0,1)) + N[1]*(bdf0*v(1,1) + bdf1*vn(1,1) + bdf2*vnn(1,1)) + N[2]*(bdf0*v(2,1) + bdf1*vn(2,1) + bdf2*vnn(2,1));
        const double crhs15 =             crhs3*crhs9 + crhs8*(DN(0,0)*v(0,1) + DN(1,0)*v(1,1) + DN(2,0)*v(2,1));
        const double crhs16 =             DN(0,1)*p[0] + DN(1,1)*p[1] + DN(2,1)*p[2] - crhs13 + rho*(crhs14 + crhs15);
        const double crhs17 =             DN(0,1)*crhs16;
        const double crhs18 =             rho*tau1*(crhs12 + crhs17);
        const double crhs19 =             N[1]*rho;
        const double crhs20 =             DN(1,0)*crhs11;
        const double crhs21 =             DN(1,1)*crhs16;
        const double crhs22 =             rho*tau1*(crhs20 + crhs21);
        const double crhs23 =             N[2]*rho;
        const double crhs24 =             DN(2,0)*crhs11;
        const double crhs25 =             DN(2,1)*crhs16;
        const double crhs26 =             rho*tau1*(crhs24 + crhs25);
        
        rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs5 - DN(0,0)*stress[0] - DN(0,1)*stress[2] + N[0]*crhs1 - crhs10*crhs6 - crhs18*crhs8 - crhs6*crhs7;
        rhs[1]=-DN(0,0)*stress[2] + DN(0,1)*crhs0 - DN(0,1)*crhs5 - DN(0,1)*stress[1] + N[0]*crhs13 - crhs14*crhs6 - crhs15*crhs6 - crhs18*crhs9;
        rhs[2]=-rho*(N[0]*crhs4 + crhs12*tau1 + crhs17*tau1);
        rhs[3]=DN(1,0)*crhs0 - DN(1,0)*crhs5 - DN(1,0)*stress[0] - DN(1,1)*stress[2] + N[1]*crhs1 - crhs10*crhs19 - crhs19*crhs7 - crhs22*crhs8;
        rhs[4]=-DN(1,0)*stress[2] + DN(1,1)*crhs0 - DN(1,1)*crhs5 - DN(1,1)*stress[1] + N[1]*crhs13 - crhs14*crhs19 - crhs15*crhs19 - crhs22*crhs9;
        rhs[5]=-rho*(N[1]*crhs4 + crhs20*tau1 + crhs21*tau1);
        rhs[6]=DN(2,0)*crhs0 - DN(2,0)*crhs5 - DN(2,0)*stress[0] - DN(2,1)*stress[2] + N[2]*crhs1 - crhs10*crhs23 - crhs23*crhs7 - crhs26*crhs8;
        rhs[7]=-DN(2,0)*stress[2] + DN(2,1)*crhs0 - DN(2,1)*crhs5 - DN(2,1)*stress[1] + N[2]*crhs13 - crhs14*crhs23 - crhs15*crhs23 - crhs26*crhs9;
        rhs[8]=-rho*(N[2]*crhs4 + crhs24*tau1 + crhs25*tau1);

    }
    
}

