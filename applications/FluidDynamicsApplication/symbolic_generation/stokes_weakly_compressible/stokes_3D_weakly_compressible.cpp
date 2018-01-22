#include "custom_elements/stokes_3D_twofluid.h"

namespace Kratos {


void Stokes3DTwoFluid::ComputeGaussPointLHSContribution(bounded_matrix<double,16,16>& lhs, const element_data<4,3>& data)
    {
        const int nnodes = 4;
        const int dim = 3;
        
        const double rho = inner_prod(data.N, data.rho);
        const double& bdf0 = data.bdf0;
        
        //get constitutive matrix 
        const Matrix& C = data.C;
        
        //get shape function values
        const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;
        const array_1d<double,nnodes>& N = data.N;
        
        
        //compute an equivalent tau by Bitrans*c*Bi
        const double tau_denom = replace_tau_denom
        const double tau1 = 1.0/(tau_denom*rho);
        const double tau2 = (C(3,3) + C(4,4) + C(5,5))/(6.0*rho);
        
        const double clhs0 =             pow(N[0], 2)*bdf0;
const double clhs1 =             clhs0*r[0];
const double clhs2 =             k*r[0]*tau2;
const double clhs3 =             C(0,0)*DN(0,0) + C(0,3)*DN(0,1) + C(0,5)*DN(0,2);
const double clhs4 =             C(0,3)*DN(0,0);
const double clhs5 =             C(3,3)*DN(0,1) + C(3,5)*DN(0,2) + clhs4;
const double clhs6 =             C(0,5)*DN(0,0);
const double clhs7 =             C(3,5)*DN(0,1) + C(5,5)*DN(0,2) + clhs6;
const double clhs8 =             DN(0,0)*k*r[0]*tau2;
const double clhs9 =             DN(0,1)*clhs8;
const double clhs10 =             C(0,1)*DN(0,1) + C(0,4)*DN(0,2) + clhs4;
const double clhs11 =             C(1,3)*DN(0,1);
const double clhs12 =             C(3,3)*DN(0,0) + C(3,4)*DN(0,2) + clhs11;
const double clhs13 =             C(3,5)*DN(0,0);
const double clhs14 =             C(4,5)*DN(0,2);
const double clhs15 =             C(1,5)*DN(0,1) + clhs13 + clhs14;
const double clhs16 =             DN(0,2)*clhs8;
const double clhs17 =             C(0,2)*DN(0,2) + C(0,4)*DN(0,1) + clhs6;
const double clhs18 =             C(3,4)*DN(0,1);
const double clhs19 =             C(2,3)*DN(0,2) + clhs13 + clhs18;
const double clhs20 =             C(2,5)*DN(0,2);
const double clhs21 =             C(4,5)*DN(0,1) + C(5,5)*DN(0,0) + clhs20;
const double clhs22 =             bdf0*tau2 - 1;
const double clhs23 =             DN(0,0)*clhs22;
const double clhs24 =             N[0]*bdf0;
const double clhs25 =             N[1]*clhs24;
const double clhs26 =             clhs25*r[1];
const double clhs27 =             DN(0,0)*k*r[1]*tau2;
const double clhs28 =             C(0,0)*DN(1,0) + C(0,3)*DN(1,1) + C(0,5)*DN(1,2);
const double clhs29 =             C(0,3)*DN(1,0);
const double clhs30 =             C(3,3)*DN(1,1) + C(3,5)*DN(1,2) + clhs29;
const double clhs31 =             C(0,5)*DN(1,0);
const double clhs32 =             C(3,5)*DN(1,1) + C(5,5)*DN(1,2) + clhs31;
const double clhs33 =             C(0,1)*DN(1,1) + C(0,4)*DN(1,2) + clhs29;
const double clhs34 =             C(1,3)*DN(1,1);
const double clhs35 =             C(3,3)*DN(1,0) + C(3,4)*DN(1,2) + clhs34;
const double clhs36 =             C(3,5)*DN(1,0);
const double clhs37 =             C(4,5)*DN(1,2);
const double clhs38 =             C(1,5)*DN(1,1) + clhs36 + clhs37;
const double clhs39 =             C(0,2)*DN(1,2) + C(0,4)*DN(1,1) + clhs31;
const double clhs40 =             C(3,4)*DN(1,1);
const double clhs41 =             C(2,3)*DN(1,2) + clhs36 + clhs40;
const double clhs42 =             C(2,5)*DN(1,2);
const double clhs43 =             C(4,5)*DN(1,1) + C(5,5)*DN(1,0) + clhs42;
const double clhs44 =             N[2]*clhs24;
const double clhs45 =             clhs44*r[2];
const double clhs46 =             DN(0,0)*k*r[2]*tau2;
const double clhs47 =             C(0,0)*DN(2,0) + C(0,3)*DN(2,1) + C(0,5)*DN(2,2);
const double clhs48 =             C(0,3)*DN(2,0);
const double clhs49 =             C(3,3)*DN(2,1) + C(3,5)*DN(2,2) + clhs48;
const double clhs50 =             C(0,5)*DN(2,0);
const double clhs51 =             C(3,5)*DN(2,1) + C(5,5)*DN(2,2) + clhs50;
const double clhs52 =             C(0,1)*DN(2,1) + C(0,4)*DN(2,2) + clhs48;
const double clhs53 =             C(1,3)*DN(2,1);
const double clhs54 =             C(3,3)*DN(2,0) + C(3,4)*DN(2,2) + clhs53;
const double clhs55 =             C(3,5)*DN(2,0);
const double clhs56 =             C(4,5)*DN(2,2);
const double clhs57 =             C(1,5)*DN(2,1) + clhs55 + clhs56;
const double clhs58 =             C(0,2)*DN(2,2) + C(0,4)*DN(2,1) + clhs50;
const double clhs59 =             C(3,4)*DN(2,1);
const double clhs60 =             C(2,3)*DN(2,2) + clhs55 + clhs59;
const double clhs61 =             C(2,5)*DN(2,2);
const double clhs62 =             C(4,5)*DN(2,1) + C(5,5)*DN(2,0) + clhs61;
const double clhs63 =             N[3]*clhs24;
const double clhs64 =             clhs63*r[3];
const double clhs65 =             DN(0,0)*k*r[3]*tau2;
const double clhs66 =             C(0,0)*DN(3,0) + C(0,3)*DN(3,1) + C(0,5)*DN(3,2);
const double clhs67 =             C(0,3)*DN(3,0);
const double clhs68 =             C(3,3)*DN(3,1) + C(3,5)*DN(3,2) + clhs67;
const double clhs69 =             C(0,5)*DN(3,0);
const double clhs70 =             C(3,5)*DN(3,1) + C(5,5)*DN(3,2) + clhs69;
const double clhs71 =             C(0,1)*DN(3,1) + C(0,4)*DN(3,2) + clhs67;
const double clhs72 =             C(1,3)*DN(3,1);
const double clhs73 =             C(3,3)*DN(3,0) + C(3,4)*DN(3,2) + clhs72;
const double clhs74 =             C(3,5)*DN(3,0);
const double clhs75 =             C(4,5)*DN(3,2);
const double clhs76 =             C(1,5)*DN(3,1) + clhs74 + clhs75;
const double clhs77 =             C(0,2)*DN(3,2) + C(0,4)*DN(3,1) + clhs69;
const double clhs78 =             C(3,4)*DN(3,1);
const double clhs79 =             C(2,3)*DN(3,2) + clhs74 + clhs78;
const double clhs80 =             C(2,5)*DN(3,2);
const double clhs81 =             C(4,5)*DN(3,1) + C(5,5)*DN(3,0) + clhs80;
const double clhs82 =             C(0,1)*DN(0,0) + C(1,5)*DN(0,2) + clhs11;
const double clhs83 =             C(0,4)*DN(0,0) + clhs14 + clhs18;
const double clhs84 =             C(1,1)*DN(0,1) + C(1,3)*DN(0,0) + C(1,4)*DN(0,2);
const double clhs85 =             C(1,4)*DN(0,1);
const double clhs86 =             C(3,4)*DN(0,0) + C(4,4)*DN(0,2) + clhs85;
const double clhs87 =             DN(0,1)*k*r[0]*tau2;
const double clhs88 =             DN(0,2)*clhs87;
const double clhs89 =             C(1,2)*DN(0,2) + C(1,5)*DN(0,0) + clhs85;
const double clhs90 =             C(2,4)*DN(0,2);
const double clhs91 =             C(4,4)*DN(0,1) + C(4,5)*DN(0,0) + clhs90;
const double clhs92 =             DN(0,1)*clhs22;
const double clhs93 =             DN(0,1)*k*r[1]*tau2;
const double clhs94 =             C(0,1)*DN(1,0) + C(1,5)*DN(1,2) + clhs34;
const double clhs95 =             C(0,4)*DN(1,0) + clhs37 + clhs40;
const double clhs96 =             C(1,1)*DN(1,1) + C(1,3)*DN(1,0) + C(1,4)*DN(1,2);
const double clhs97 =             C(1,4)*DN(1,1);
const double clhs98 =             C(3,4)*DN(1,0) + C(4,4)*DN(1,2) + clhs97;
const double clhs99 =             C(1,2)*DN(1,2) + C(1,5)*DN(1,0) + clhs97;
const double clhs100 =             C(2,4)*DN(1,2);
const double clhs101 =             C(4,4)*DN(1,1) + C(4,5)*DN(1,0) + clhs100;
const double clhs102 =             DN(0,1)*k*r[2]*tau2;
const double clhs103 =             C(0,1)*DN(2,0) + C(1,5)*DN(2,2) + clhs53;
const double clhs104 =             C(0,4)*DN(2,0) + clhs56 + clhs59;
const double clhs105 =             C(1,1)*DN(2,1) + C(1,3)*DN(2,0) + C(1,4)*DN(2,2);
const double clhs106 =             C(1,4)*DN(2,1);
const double clhs107 =             C(3,4)*DN(2,0) + C(4,4)*DN(2,2) + clhs106;
const double clhs108 =             C(1,2)*DN(2,2) + C(1,5)*DN(2,0) + clhs106;
const double clhs109 =             C(2,4)*DN(2,2);
const double clhs110 =             C(4,4)*DN(2,1) + C(4,5)*DN(2,0) + clhs109;
const double clhs111 =             DN(0,1)*k*r[3]*tau2;
const double clhs112 =             C(0,1)*DN(3,0) + C(1,5)*DN(3,2) + clhs72;
const double clhs113 =             C(0,4)*DN(3,0) + clhs75 + clhs78;
const double clhs114 =             C(1,1)*DN(3,1) + C(1,3)*DN(3,0) + C(1,4)*DN(3,2);
const double clhs115 =             C(1,4)*DN(3,1);
const double clhs116 =             C(3,4)*DN(3,0) + C(4,4)*DN(3,2) + clhs115;
const double clhs117 =             C(1,2)*DN(3,2) + C(1,5)*DN(3,0) + clhs115;
const double clhs118 =             C(2,4)*DN(3,2);
const double clhs119 =             C(4,4)*DN(3,1) + C(4,5)*DN(3,0) + clhs118;
const double clhs120 =             C(0,2)*DN(0,0) + C(2,3)*DN(0,1) + clhs20;
const double clhs121 =             C(1,2)*DN(0,1) + C(2,3)*DN(0,0) + clhs90;
const double clhs122 =             C(2,2)*DN(0,2) + C(2,4)*DN(0,1) + C(2,5)*DN(0,0);
const double clhs123 =             DN(0,2)*clhs22;
const double clhs124 =             DN(0,2)*DN(1,0)*k*tau2;
const double clhs125 =             C(0,2)*DN(1,0) + C(2,3)*DN(1,1) + clhs42;
const double clhs126 =             DN(0,2)*DN(1,1)*k*tau2;
const double clhs127 =             C(1,2)*DN(1,1) + C(2,3)*DN(1,0) + clhs100;
const double clhs128 =             DN(0,2)*DN(1,2)*k*tau2;
const double clhs129 =             C(2,2)*DN(1,2) + C(2,4)*DN(1,1) + C(2,5)*DN(1,0);
const double clhs130 =             DN(0,2)*DN(2,0)*k*tau2;
const double clhs131 =             C(0,2)*DN(2,0) + C(2,3)*DN(2,1) + clhs61;
const double clhs132 =             DN(0,2)*DN(2,1)*k*tau2;
const double clhs133 =             C(1,2)*DN(2,1) + C(2,3)*DN(2,0) + clhs109;
const double clhs134 =             DN(0,2)*DN(2,2)*k*tau2;
const double clhs135 =             C(2,2)*DN(2,2) + C(2,4)*DN(2,1) + C(2,5)*DN(2,0);
const double clhs136 =             DN(0,2)*DN(3,0)*k*tau2;
const double clhs137 =             C(0,2)*DN(3,0) + C(2,3)*DN(3,1) + clhs80;
const double clhs138 =             DN(0,2)*DN(3,1)*k*tau2;
const double clhs139 =             C(1,2)*DN(3,1) + C(2,3)*DN(3,0) + clhs118;
const double clhs140 =             DN(0,2)*DN(3,2)*k*tau2;
const double clhs141 =             C(2,2)*DN(3,2) + C(2,4)*DN(3,1) + C(2,5)*DN(3,0);
const double clhs142 =             bdf0*tau1 + k;
const double clhs143 =             N[0]*clhs142*r[0];
const double clhs144 =             DN(1,0)*k;
const double clhs145 =             DN(0,0)*bdf0*tau1;
const double clhs146 =             DN(1,1)*k;
const double clhs147 =             DN(0,1)*bdf0*tau1;
const double clhs148 =             DN(1,2)*k;
const double clhs149 =             DN(0,2)*bdf0*tau1;
const double clhs150 =             DN(2,0)*k;
const double clhs151 =             DN(2,1)*k;
const double clhs152 =             DN(2,2)*k;
const double clhs153 =             DN(3,0)*k;
const double clhs154 =             DN(3,1)*k;
const double clhs155 =             DN(3,2)*k;
const double clhs156 =             clhs25*r[0];
const double clhs157 =             DN(1,0)*clhs22;
const double clhs158 =             pow(N[1], 2)*bdf0;
const double clhs159 =             clhs158*r[1];
const double clhs160 =             k*r[1]*tau2;
const double clhs161 =             DN(1,0)*k*r[1]*tau2;
const double clhs162 =             DN(1,1)*clhs161;
const double clhs163 =             DN(1,2)*clhs161;
const double clhs164 =             N[1]*bdf0;
const double clhs165 =             N[2]*clhs164;
const double clhs166 =             clhs165*r[2];
const double clhs167 =             DN(1,0)*k*r[2]*tau2;
const double clhs168 =             N[3]*clhs164;
const double clhs169 =             clhs168*r[3];
const double clhs170 =             DN(1,0)*k*r[3]*tau2;
const double clhs171 =             DN(1,1)*clhs22;
const double clhs172 =             DN(1,1)*k*r[1]*tau2;
const double clhs173 =             DN(1,2)*clhs172;
const double clhs174 =             DN(1,1)*k*r[2]*tau2;
const double clhs175 =             DN(1,1)*k*r[3]*tau2;
const double clhs176 =             DN(1,2)*clhs22;
const double clhs177 =             DN(1,2)*DN(2,0)*k*tau2;
const double clhs178 =             DN(1,2)*DN(2,1)*k*tau2;
const double clhs179 =             DN(1,2)*DN(2,2)*k*tau2;
const double clhs180 =             DN(1,2)*DN(3,0)*k*tau2;
const double clhs181 =             DN(1,2)*DN(3,1)*k*tau2;
const double clhs182 =             DN(1,2)*DN(3,2)*k*tau2;
const double clhs183 =             DN(0,0)*k;
const double clhs184 =             DN(1,0)*bdf0*tau1;
const double clhs185 =             DN(0,1)*k;
const double clhs186 =             DN(1,1)*bdf0*tau1;
const double clhs187 =             DN(0,2)*k;
const double clhs188 =             DN(1,2)*bdf0*tau1;
const double clhs189 =             N[1]*clhs142*r[1];
const double clhs190 =             clhs44*r[0];
const double clhs191 =             DN(2,0)*clhs22;
const double clhs192 =             clhs165*r[1];
const double clhs193 =             pow(N[2], 2)*bdf0;
const double clhs194 =             clhs193*r[2];
const double clhs195 =             k*r[2]*tau2;
const double clhs196 =             DN(2,0)*k*r[2]*tau2;
const double clhs197 =             DN(2,1)*clhs196;
const double clhs198 =             DN(2,2)*clhs196;
const double clhs199 =             N[2]*N[3]*bdf0;
const double clhs200 =             clhs199*r[3];
const double clhs201 =             DN(2,0)*k*r[3]*tau2;
const double clhs202 =             DN(2,1)*clhs22;
const double clhs203 =             DN(2,1)*k*r[2]*tau2;
const double clhs204 =             DN(2,2)*clhs203;
const double clhs205 =             DN(2,1)*k*r[3]*tau2;
const double clhs206 =             DN(2,2)*clhs22;
const double clhs207 =             DN(2,2)*DN(3,0)*k*tau2;
const double clhs208 =             DN(2,2)*DN(3,1)*k*tau2;
const double clhs209 =             DN(2,2)*DN(3,2)*k*tau2;
const double clhs210 =             DN(2,0)*bdf0*tau1;
const double clhs211 =             DN(2,1)*bdf0*tau1;
const double clhs212 =             DN(2,2)*bdf0*tau1;
const double clhs213 =             N[2]*clhs142*r[2];
const double clhs214 =             clhs63*r[0];
const double clhs215 =             DN(3,0)*clhs22;
const double clhs216 =             clhs168*r[1];
const double clhs217 =             clhs199*r[2];
const double clhs218 =             pow(N[3], 2)*bdf0;
const double clhs219 =             clhs218*r[3];
const double clhs220 =             k*r[3]*tau2;
const double clhs221 =             DN(3,0)*k*r[3]*tau2;
const double clhs222 =             DN(3,1)*clhs221;
const double clhs223 =             DN(3,2)*clhs221;
const double clhs224 =             DN(3,1)*clhs22;
const double clhs225 =             (r[3]*tau2)*(DN(3,1)*DN(3,2)*k);
const double clhs226 =             DN(3,2)*clhs22;
const double clhs227 =             DN(3,0)*bdf0*tau1;
const double clhs228 =             DN(3,1)*bdf0*tau1;
const double clhs229 =             DN(3,2)*bdf0*tau1;
const double clhs230 =             N[3]*clhs142*r[3];
            lhs(0,0)=pow(DN(0,0), 2)*clhs2 + DN(0,0)*clhs3 + DN(0,1)*clhs5 + DN(0,2)*clhs7 + clhs1;
            lhs(0,1)=DN(0,0)*clhs10 + DN(0,1)*clhs12 + DN(0,2)*clhs15 + clhs9;
            lhs(0,2)=DN(0,0)*clhs17 + DN(0,1)*clhs19 + DN(0,2)*clhs21 + clhs16;
            lhs(0,3)=N[0]*clhs23;
            lhs(0,4)=DN(0,0)*clhs28 + DN(0,1)*clhs30 + DN(0,2)*clhs32 + DN(1,0)*clhs27 + clhs26;
            lhs(0,5)=DN(0,0)*clhs33 + DN(0,1)*clhs35 + DN(0,2)*clhs38 + DN(1,1)*clhs27;
            lhs(0,6)=DN(0,0)*clhs39 + DN(0,1)*clhs41 + DN(0,2)*clhs43 + DN(1,2)*clhs27;
            lhs(0,7)=N[1]*clhs23;
            lhs(0,8)=DN(0,0)*clhs47 + DN(0,1)*clhs49 + DN(0,2)*clhs51 + DN(2,0)*clhs46 + clhs45;
            lhs(0,9)=DN(0,0)*clhs52 + DN(0,1)*clhs54 + DN(0,2)*clhs57 + DN(2,1)*clhs46;
            lhs(0,10)=DN(0,0)*clhs58 + DN(0,1)*clhs60 + DN(0,2)*clhs62 + DN(2,2)*clhs46;
            lhs(0,11)=N[2]*clhs23;
            lhs(0,12)=DN(0,0)*clhs66 + DN(0,1)*clhs68 + DN(0,2)*clhs70 + DN(3,0)*clhs65 + clhs64;
            lhs(0,13)=DN(0,0)*clhs71 + DN(0,1)*clhs73 + DN(0,2)*clhs76 + DN(3,1)*clhs65;
            lhs(0,14)=DN(0,0)*clhs77 + DN(0,1)*clhs79 + DN(0,2)*clhs81 + DN(3,2)*clhs65;
            lhs(0,15)=N[3]*clhs23;
            lhs(1,0)=DN(0,0)*clhs5 + DN(0,1)*clhs82 + DN(0,2)*clhs83 + clhs9;
            lhs(1,1)=DN(0,0)*clhs12 + pow(DN(0,1), 2)*clhs2 + DN(0,1)*clhs84 + DN(0,2)*clhs86 + clhs1;
            lhs(1,2)=DN(0,0)*clhs19 + DN(0,1)*clhs89 + DN(0,2)*clhs91 + clhs88;
            lhs(1,3)=N[0]*clhs92;
            lhs(1,4)=DN(0,0)*clhs30 + DN(0,1)*clhs94 + DN(0,2)*clhs95 + DN(1,0)*clhs93;
            lhs(1,5)=DN(0,0)*clhs35 + DN(0,1)*clhs96 + DN(0,2)*clhs98 + DN(1,1)*clhs93 + clhs26;
            lhs(1,6)=DN(0,0)*clhs41 + DN(0,1)*clhs99 + DN(0,2)*clhs101 + DN(1,2)*clhs93;
            lhs(1,7)=N[1]*clhs92;
            lhs(1,8)=DN(0,0)*clhs49 + DN(0,1)*clhs103 + DN(0,2)*clhs104 + DN(2,0)*clhs102;
            lhs(1,9)=DN(0,0)*clhs54 + DN(0,1)*clhs105 + DN(0,2)*clhs107 + DN(2,1)*clhs102 + clhs45;
            lhs(1,10)=DN(0,0)*clhs60 + DN(0,1)*clhs108 + DN(0,2)*clhs110 + DN(2,2)*clhs102;
            lhs(1,11)=N[2]*clhs92;
            lhs(1,12)=DN(0,0)*clhs68 + DN(0,1)*clhs112 + DN(0,2)*clhs113 + DN(3,0)*clhs111;
            lhs(1,13)=DN(0,0)*clhs73 + DN(0,1)*clhs114 + DN(0,2)*clhs116 + DN(3,1)*clhs111 + clhs64;
            lhs(1,14)=DN(0,0)*clhs79 + DN(0,1)*clhs117 + DN(0,2)*clhs119 + DN(3,2)*clhs111;
            lhs(1,15)=N[3]*clhs92;
            lhs(2,0)=DN(0,0)*clhs7 + DN(0,1)*clhs83 + DN(0,2)*clhs120 + clhs16;
            lhs(2,1)=DN(0,0)*clhs15 + DN(0,1)*clhs86 + DN(0,2)*clhs121 + clhs88;
            lhs(2,2)=DN(0,0)*clhs21 + DN(0,1)*clhs91 + pow(DN(0,2), 2)*clhs2 + DN(0,2)*clhs122 + clhs1;
            lhs(2,3)=N[0]*clhs123;
            lhs(2,4)=DN(0,0)*clhs32 + DN(0,1)*clhs95 + DN(0,2)*clhs125 + clhs124*r[1];
            lhs(2,5)=DN(0,0)*clhs38 + DN(0,1)*clhs98 + DN(0,2)*clhs127 + clhs126*r[1];
            lhs(2,6)=DN(0,0)*clhs43 + DN(0,1)*clhs101 + DN(0,2)*clhs129 + clhs128*r[1] + clhs26;
            lhs(2,7)=N[1]*clhs123;
            lhs(2,8)=DN(0,0)*clhs51 + DN(0,1)*clhs104 + DN(0,2)*clhs131 + clhs130*r[2];
            lhs(2,9)=DN(0,0)*clhs57 + DN(0,1)*clhs107 + DN(0,2)*clhs133 + clhs132*r[2];
            lhs(2,10)=DN(0,0)*clhs62 + DN(0,1)*clhs110 + DN(0,2)*clhs135 + clhs134*r[2] + clhs45;
            lhs(2,11)=N[2]*clhs123;
            lhs(2,12)=DN(0,0)*clhs70 + DN(0,1)*clhs113 + DN(0,2)*clhs137 + clhs136*r[3];
            lhs(2,13)=DN(0,0)*clhs76 + DN(0,1)*clhs116 + DN(0,2)*clhs139 + clhs138*r[3];
            lhs(2,14)=DN(0,0)*clhs81 + DN(0,1)*clhs119 + DN(0,2)*clhs141 + clhs140*r[3] + clhs64;
            lhs(2,15)=N[3]*clhs123;
            lhs(3,0)=DN(0,0)*clhs143;
            lhs(3,1)=DN(0,1)*clhs143;
            lhs(3,2)=DN(0,2)*clhs143;
            lhs(3,3)=clhs0;
            lhs(3,4)=r[1]*(N[0]*clhs144 + N[1]*clhs145);
            lhs(3,5)=r[1]*(N[0]*clhs146 + N[1]*clhs147);
            lhs(3,6)=r[1]*(N[0]*clhs148 + N[1]*clhs149);
            lhs(3,7)=clhs25;
            lhs(3,8)=r[2]*(N[0]*clhs150 + N[2]*clhs145);
            lhs(3,9)=r[2]*(N[0]*clhs151 + N[2]*clhs147);
            lhs(3,10)=r[2]*(N[0]*clhs152 + N[2]*clhs149);
            lhs(3,11)=clhs44;
            lhs(3,12)=r[3]*(N[0]*clhs153 + N[3]*clhs145);
            lhs(3,13)=r[3]*(N[0]*clhs154 + N[3]*clhs147);
            lhs(3,14)=r[3]*(N[0]*clhs155 + N[3]*clhs149);
            lhs(3,15)=clhs63;
            lhs(4,0)=DN(1,0)*clhs3 + DN(1,0)*clhs8 + DN(1,1)*clhs5 + DN(1,2)*clhs7 + clhs156;
            lhs(4,1)=DN(1,0)*clhs10 + DN(1,0)*clhs87 + DN(1,1)*clhs12 + DN(1,2)*clhs15;
            lhs(4,2)=DN(1,0)*clhs17 + DN(1,1)*clhs19 + DN(1,2)*clhs21 + clhs124*r[0];
            lhs(4,3)=N[0]*clhs157;
            lhs(4,4)=pow(DN(1,0), 2)*clhs160 + DN(1,0)*clhs28 + DN(1,1)*clhs30 + DN(1,2)*clhs32 + clhs159;
            lhs(4,5)=DN(1,0)*clhs33 + DN(1,1)*clhs35 + DN(1,2)*clhs38 + clhs162;
            lhs(4,6)=DN(1,0)*clhs39 + DN(1,1)*clhs41 + DN(1,2)*clhs43 + clhs163;
            lhs(4,7)=N[1]*clhs157;
            lhs(4,8)=DN(1,0)*clhs47 + DN(1,1)*clhs49 + DN(1,2)*clhs51 + DN(2,0)*clhs167 + clhs166;
            lhs(4,9)=DN(1,0)*clhs52 + DN(1,1)*clhs54 + DN(1,2)*clhs57 + DN(2,1)*clhs167;
            lhs(4,10)=DN(1,0)*clhs58 + DN(1,1)*clhs60 + DN(1,2)*clhs62 + DN(2,2)*clhs167;
            lhs(4,11)=N[2]*clhs157;
            lhs(4,12)=DN(1,0)*clhs66 + DN(1,1)*clhs68 + DN(1,2)*clhs70 + DN(3,0)*clhs170 + clhs169;
            lhs(4,13)=DN(1,0)*clhs71 + DN(1,1)*clhs73 + DN(1,2)*clhs76 + DN(3,1)*clhs170;
            lhs(4,14)=DN(1,0)*clhs77 + DN(1,1)*clhs79 + DN(1,2)*clhs81 + DN(3,2)*clhs170;
            lhs(4,15)=N[3]*clhs157;
            lhs(5,0)=DN(1,0)*clhs5 + DN(1,1)*clhs8 + DN(1,1)*clhs82 + DN(1,2)*clhs83;
            lhs(5,1)=DN(1,0)*clhs12 + DN(1,1)*clhs84 + DN(1,1)*clhs87 + DN(1,2)*clhs86 + clhs156;
            lhs(5,2)=DN(1,0)*clhs19 + DN(1,1)*clhs89 + DN(1,2)*clhs91 + clhs126*r[0];
            lhs(5,3)=N[0]*clhs171;
            lhs(5,4)=DN(1,0)*clhs30 + DN(1,1)*clhs94 + DN(1,2)*clhs95 + clhs162;
            lhs(5,5)=DN(1,0)*clhs35 + pow(DN(1,1), 2)*clhs160 + DN(1,1)*clhs96 + DN(1,2)*clhs98 + clhs159;
            lhs(5,6)=DN(1,0)*clhs41 + DN(1,1)*clhs99 + DN(1,2)*clhs101 + clhs173;
            lhs(5,7)=N[1]*clhs171;
            lhs(5,8)=DN(1,0)*clhs49 + DN(1,1)*clhs103 + DN(1,2)*clhs104 + DN(2,0)*clhs174;
            lhs(5,9)=DN(1,0)*clhs54 + DN(1,1)*clhs105 + DN(1,2)*clhs107 + DN(2,1)*clhs174 + clhs166;
            lhs(5,10)=DN(1,0)*clhs60 + DN(1,1)*clhs108 + DN(1,2)*clhs110 + DN(2,2)*clhs174;
            lhs(5,11)=N[2]*clhs171;
            lhs(5,12)=DN(1,0)*clhs68 + DN(1,1)*clhs112 + DN(1,2)*clhs113 + DN(3,0)*clhs175;
            lhs(5,13)=DN(1,0)*clhs73 + DN(1,1)*clhs114 + DN(1,2)*clhs116 + DN(3,1)*clhs175 + clhs169;
            lhs(5,14)=DN(1,0)*clhs79 + DN(1,1)*clhs117 + DN(1,2)*clhs119 + DN(3,2)*clhs175;
            lhs(5,15)=N[3]*clhs171;
            lhs(6,0)=DN(1,0)*clhs7 + DN(1,1)*clhs83 + DN(1,2)*clhs120 + DN(1,2)*clhs8;
            lhs(6,1)=DN(1,0)*clhs15 + DN(1,1)*clhs86 + DN(1,2)*clhs121 + DN(1,2)*clhs87;
            lhs(6,2)=DN(1,0)*clhs21 + DN(1,1)*clhs91 + DN(1,2)*clhs122 + clhs128*r[0] + clhs156;
            lhs(6,3)=N[0]*clhs176;
            lhs(6,4)=DN(1,0)*clhs32 + DN(1,1)*clhs95 + DN(1,2)*clhs125 + clhs163;
            lhs(6,5)=DN(1,0)*clhs38 + DN(1,1)*clhs98 + DN(1,2)*clhs127 + clhs173;
            lhs(6,6)=DN(1,0)*clhs43 + DN(1,1)*clhs101 + pow(DN(1,2), 2)*clhs160 + DN(1,2)*clhs129 + clhs159;
            lhs(6,7)=N[1]*clhs176;
            lhs(6,8)=DN(1,0)*clhs51 + DN(1,1)*clhs104 + DN(1,2)*clhs131 + clhs177*r[2];
            lhs(6,9)=DN(1,0)*clhs57 + DN(1,1)*clhs107 + DN(1,2)*clhs133 + clhs178*r[2];
            lhs(6,10)=DN(1,0)*clhs62 + DN(1,1)*clhs110 + DN(1,2)*clhs135 + clhs166 + clhs179*r[2];
            lhs(6,11)=N[2]*clhs176;
            lhs(6,12)=DN(1,0)*clhs70 + DN(1,1)*clhs113 + DN(1,2)*clhs137 + clhs180*r[3];
            lhs(6,13)=DN(1,0)*clhs76 + DN(1,1)*clhs116 + DN(1,2)*clhs139 + clhs181*r[3];
            lhs(6,14)=DN(1,0)*clhs81 + DN(1,1)*clhs119 + DN(1,2)*clhs141 + clhs169 + clhs182*r[3];
            lhs(6,15)=N[3]*clhs176;
            lhs(7,0)=r[0]*(N[0]*clhs184 + N[1]*clhs183);
            lhs(7,1)=r[0]*(N[0]*clhs186 + N[1]*clhs185);
            lhs(7,2)=r[0]*(N[0]*clhs188 + N[1]*clhs187);
            lhs(7,3)=clhs25;
            lhs(7,4)=DN(1,0)*clhs189;
            lhs(7,5)=DN(1,1)*clhs189;
            lhs(7,6)=DN(1,2)*clhs189;
            lhs(7,7)=clhs158;
            lhs(7,8)=r[2]*(N[1]*clhs150 + N[2]*clhs184);
            lhs(7,9)=r[2]*(N[1]*clhs151 + N[2]*clhs186);
            lhs(7,10)=r[2]*(N[1]*clhs152 + N[2]*clhs188);
            lhs(7,11)=clhs165;
            lhs(7,12)=r[3]*(N[1]*clhs153 + N[3]*clhs184);
            lhs(7,13)=r[3]*(N[1]*clhs154 + N[3]*clhs186);
            lhs(7,14)=r[3]*(N[1]*clhs155 + N[3]*clhs188);
            lhs(7,15)=clhs168;
            lhs(8,0)=DN(2,0)*clhs3 + DN(2,0)*clhs8 + DN(2,1)*clhs5 + DN(2,2)*clhs7 + clhs190;
            lhs(8,1)=DN(2,0)*clhs10 + DN(2,0)*clhs87 + DN(2,1)*clhs12 + DN(2,2)*clhs15;
            lhs(8,2)=DN(2,0)*clhs17 + DN(2,1)*clhs19 + DN(2,2)*clhs21 + clhs130*r[0];
            lhs(8,3)=N[0]*clhs191;
            lhs(8,4)=DN(2,0)*clhs161 + DN(2,0)*clhs28 + DN(2,1)*clhs30 + DN(2,2)*clhs32 + clhs192;
            lhs(8,5)=DN(2,0)*clhs172 + DN(2,0)*clhs33 + DN(2,1)*clhs35 + DN(2,2)*clhs38;
            lhs(8,6)=DN(2,0)*clhs39 + DN(2,1)*clhs41 + DN(2,2)*clhs43 + clhs177*r[1];
            lhs(8,7)=N[1]*clhs191;
            lhs(8,8)=pow(DN(2,0), 2)*clhs195 + DN(2,0)*clhs47 + DN(2,1)*clhs49 + DN(2,2)*clhs51 + clhs194;
            lhs(8,9)=DN(2,0)*clhs52 + DN(2,1)*clhs54 + DN(2,2)*clhs57 + clhs197;
            lhs(8,10)=DN(2,0)*clhs58 + DN(2,1)*clhs60 + DN(2,2)*clhs62 + clhs198;
            lhs(8,11)=N[2]*clhs191;
            lhs(8,12)=DN(2,0)*clhs66 + DN(2,1)*clhs68 + DN(2,2)*clhs70 + DN(3,0)*clhs201 + clhs200;
            lhs(8,13)=DN(2,0)*clhs71 + DN(2,1)*clhs73 + DN(2,2)*clhs76 + DN(3,1)*clhs201;
            lhs(8,14)=DN(2,0)*clhs77 + DN(2,1)*clhs79 + DN(2,2)*clhs81 + DN(3,2)*clhs201;
            lhs(8,15)=N[3]*clhs191;
            lhs(9,0)=DN(2,0)*clhs5 + DN(2,1)*clhs8 + DN(2,1)*clhs82 + DN(2,2)*clhs83;
            lhs(9,1)=DN(2,0)*clhs12 + DN(2,1)*clhs84 + DN(2,1)*clhs87 + DN(2,2)*clhs86 + clhs190;
            lhs(9,2)=DN(2,0)*clhs19 + DN(2,1)*clhs89 + DN(2,2)*clhs91 + clhs132*r[0];
            lhs(9,3)=N[0]*clhs202;
            lhs(9,4)=DN(2,0)*clhs30 + DN(2,1)*clhs161 + DN(2,1)*clhs94 + DN(2,2)*clhs95;
            lhs(9,5)=DN(2,0)*clhs35 + DN(2,1)*clhs172 + DN(2,1)*clhs96 + DN(2,2)*clhs98 + clhs192;
            lhs(9,6)=DN(2,0)*clhs41 + DN(2,1)*clhs99 + DN(2,2)*clhs101 + clhs178*r[1];
            lhs(9,7)=N[1]*clhs202;
            lhs(9,8)=DN(2,0)*clhs49 + DN(2,1)*clhs103 + DN(2,2)*clhs104 + clhs197;
            lhs(9,9)=DN(2,0)*clhs54 + pow(DN(2,1), 2)*clhs195 + DN(2,1)*clhs105 + DN(2,2)*clhs107 + clhs194;
            lhs(9,10)=DN(2,0)*clhs60 + DN(2,1)*clhs108 + DN(2,2)*clhs110 + clhs204;
            lhs(9,11)=N[2]*clhs202;
            lhs(9,12)=DN(2,0)*clhs68 + DN(2,1)*clhs112 + DN(2,2)*clhs113 + DN(3,0)*clhs205;
            lhs(9,13)=DN(2,0)*clhs73 + DN(2,1)*clhs114 + DN(2,2)*clhs116 + DN(3,1)*clhs205 + clhs200;
            lhs(9,14)=DN(2,0)*clhs79 + DN(2,1)*clhs117 + DN(2,2)*clhs119 + DN(3,2)*clhs205;
            lhs(9,15)=N[3]*clhs202;
            lhs(10,0)=DN(2,0)*clhs7 + DN(2,1)*clhs83 + DN(2,2)*clhs120 + DN(2,2)*clhs8;
            lhs(10,1)=DN(2,0)*clhs15 + DN(2,1)*clhs86 + DN(2,2)*clhs121 + DN(2,2)*clhs87;
            lhs(10,2)=DN(2,0)*clhs21 + DN(2,1)*clhs91 + DN(2,2)*clhs122 + clhs134*r[0] + clhs190;
            lhs(10,3)=N[0]*clhs206;
            lhs(10,4)=DN(2,0)*clhs32 + DN(2,1)*clhs95 + DN(2,2)*clhs125 + DN(2,2)*clhs161;
            lhs(10,5)=DN(2,0)*clhs38 + DN(2,1)*clhs98 + DN(2,2)*clhs127 + DN(2,2)*clhs172;
            lhs(10,6)=DN(2,0)*clhs43 + DN(2,1)*clhs101 + DN(2,2)*clhs129 + clhs179*r[1] + clhs192;
            lhs(10,7)=N[1]*clhs206;
            lhs(10,8)=DN(2,0)*clhs51 + DN(2,1)*clhs104 + DN(2,2)*clhs131 + clhs198;
            lhs(10,9)=DN(2,0)*clhs57 + DN(2,1)*clhs107 + DN(2,2)*clhs133 + clhs204;
            lhs(10,10)=DN(2,0)*clhs62 + DN(2,1)*clhs110 + pow(DN(2,2), 2)*clhs195 + DN(2,2)*clhs135 + clhs194;
            lhs(10,11)=N[2]*clhs206;
            lhs(10,12)=DN(2,0)*clhs70 + DN(2,1)*clhs113 + DN(2,2)*clhs137 + clhs207*r[3];
            lhs(10,13)=DN(2,0)*clhs76 + DN(2,1)*clhs116 + DN(2,2)*clhs139 + clhs208*r[3];
            lhs(10,14)=DN(2,0)*clhs81 + DN(2,1)*clhs119 + DN(2,2)*clhs141 + clhs200 + clhs209*r[3];
            lhs(10,15)=N[3]*clhs206;
            lhs(11,0)=r[0]*(N[0]*clhs210 + N[2]*clhs183);
            lhs(11,1)=r[0]*(N[0]*clhs211 + N[2]*clhs185);
            lhs(11,2)=r[0]*(N[0]*clhs212 + N[2]*clhs187);
            lhs(11,3)=clhs44;
            lhs(11,4)=r[1]*(N[1]*clhs210 + N[2]*clhs144);
            lhs(11,5)=r[1]*(N[1]*clhs211 + N[2]*clhs146);
            lhs(11,6)=r[1]*(N[1]*clhs212 + N[2]*clhs148);
            lhs(11,7)=clhs165;
            lhs(11,8)=DN(2,0)*clhs213;
            lhs(11,9)=DN(2,1)*clhs213;
            lhs(11,10)=DN(2,2)*clhs213;
            lhs(11,11)=clhs193;
            lhs(11,12)=r[3]*(N[2]*clhs153 + N[3]*clhs210);
            lhs(11,13)=r[3]*(N[2]*clhs154 + N[3]*clhs211);
            lhs(11,14)=r[3]*(N[2]*clhs155 + N[3]*clhs212);
            lhs(11,15)=clhs199;
            lhs(12,0)=DN(3,0)*clhs3 + DN(3,0)*clhs8 + DN(3,1)*clhs5 + DN(3,2)*clhs7 + clhs214;
            lhs(12,1)=DN(3,0)*clhs10 + DN(3,0)*clhs87 + DN(3,1)*clhs12 + DN(3,2)*clhs15;
            lhs(12,2)=DN(3,0)*clhs17 + DN(3,1)*clhs19 + DN(3,2)*clhs21 + clhs136*r[0];
            lhs(12,3)=N[0]*clhs215;
            lhs(12,4)=DN(3,0)*clhs161 + DN(3,0)*clhs28 + DN(3,1)*clhs30 + DN(3,2)*clhs32 + clhs216;
            lhs(12,5)=DN(3,0)*clhs172 + DN(3,0)*clhs33 + DN(3,1)*clhs35 + DN(3,2)*clhs38;
            lhs(12,6)=DN(3,0)*clhs39 + DN(3,1)*clhs41 + DN(3,2)*clhs43 + clhs180*r[1];
            lhs(12,7)=N[1]*clhs215;
            lhs(12,8)=DN(3,0)*clhs196 + DN(3,0)*clhs47 + DN(3,1)*clhs49 + DN(3,2)*clhs51 + clhs217;
            lhs(12,9)=DN(3,0)*clhs203 + DN(3,0)*clhs52 + DN(3,1)*clhs54 + DN(3,2)*clhs57;
            lhs(12,10)=DN(3,0)*clhs58 + DN(3,1)*clhs60 + DN(3,2)*clhs62 + clhs207*r[2];
            lhs(12,11)=N[2]*clhs215;
            lhs(12,12)=pow(DN(3,0), 2)*clhs220 + DN(3,0)*clhs66 + DN(3,1)*clhs68 + DN(3,2)*clhs70 + clhs219;
            lhs(12,13)=DN(3,0)*clhs71 + DN(3,1)*clhs73 + DN(3,2)*clhs76 + clhs222;
            lhs(12,14)=DN(3,0)*clhs77 + DN(3,1)*clhs79 + DN(3,2)*clhs81 + clhs223;
            lhs(12,15)=N[3]*clhs215;
            lhs(13,0)=DN(3,0)*clhs5 + DN(3,1)*clhs8 + DN(3,1)*clhs82 + DN(3,2)*clhs83;
            lhs(13,1)=DN(3,0)*clhs12 + DN(3,1)*clhs84 + DN(3,1)*clhs87 + DN(3,2)*clhs86 + clhs214;
            lhs(13,2)=DN(3,0)*clhs19 + DN(3,1)*clhs89 + DN(3,2)*clhs91 + clhs138*r[0];
            lhs(13,3)=N[0]*clhs224;
            lhs(13,4)=DN(3,0)*clhs30 + DN(3,1)*clhs161 + DN(3,1)*clhs94 + DN(3,2)*clhs95;
            lhs(13,5)=DN(3,0)*clhs35 + DN(3,1)*clhs172 + DN(3,1)*clhs96 + DN(3,2)*clhs98 + clhs216;
            lhs(13,6)=DN(3,0)*clhs41 + DN(3,1)*clhs99 + DN(3,2)*clhs101 + clhs181*r[1];
            lhs(13,7)=N[1]*clhs224;
            lhs(13,8)=DN(3,0)*clhs49 + DN(3,1)*clhs103 + DN(3,1)*clhs196 + DN(3,2)*clhs104;
            lhs(13,9)=DN(3,0)*clhs54 + DN(3,1)*clhs105 + DN(3,1)*clhs203 + DN(3,2)*clhs107 + clhs217;
            lhs(13,10)=DN(3,0)*clhs60 + DN(3,1)*clhs108 + DN(3,2)*clhs110 + clhs208*r[2];
            lhs(13,11)=N[2]*clhs224;
            lhs(13,12)=DN(3,0)*clhs68 + DN(3,1)*clhs112 + DN(3,2)*clhs113 + clhs222;
            lhs(13,13)=DN(3,0)*clhs73 + pow(DN(3,1), 2)*clhs220 + DN(3,1)*clhs114 + DN(3,2)*clhs116 + clhs219;
            lhs(13,14)=DN(3,0)*clhs79 + DN(3,1)*clhs117 + DN(3,2)*clhs119 + clhs225;
            lhs(13,15)=N[3]*clhs224;
            lhs(14,0)=DN(3,0)*clhs7 + DN(3,1)*clhs83 + DN(3,2)*clhs120 + DN(3,2)*clhs8;
            lhs(14,1)=DN(3,0)*clhs15 + DN(3,1)*clhs86 + DN(3,2)*clhs121 + DN(3,2)*clhs87;
            lhs(14,2)=DN(3,0)*clhs21 + DN(3,1)*clhs91 + DN(3,2)*clhs122 + clhs140*r[0] + clhs214;
            lhs(14,3)=N[0]*clhs226;
            lhs(14,4)=DN(3,0)*clhs32 + DN(3,1)*clhs95 + DN(3,2)*clhs125 + DN(3,2)*clhs161;
            lhs(14,5)=DN(3,0)*clhs38 + DN(3,1)*clhs98 + DN(3,2)*clhs127 + DN(3,2)*clhs172;
            lhs(14,6)=DN(3,0)*clhs43 + DN(3,1)*clhs101 + DN(3,2)*clhs129 + clhs182*r[1] + clhs216;
            lhs(14,7)=N[1]*clhs226;
            lhs(14,8)=DN(3,0)*clhs51 + DN(3,1)*clhs104 + DN(3,2)*clhs131 + DN(3,2)*clhs196;
            lhs(14,9)=DN(3,0)*clhs57 + DN(3,1)*clhs107 + DN(3,2)*clhs133 + DN(3,2)*clhs203;
            lhs(14,10)=DN(3,0)*clhs62 + DN(3,1)*clhs110 + DN(3,2)*clhs135 + clhs209*r[2] + clhs217;
            lhs(14,11)=N[2]*clhs226;
            lhs(14,12)=DN(3,0)*clhs70 + DN(3,1)*clhs113 + DN(3,2)*clhs137 + clhs223;
            lhs(14,13)=DN(3,0)*clhs76 + DN(3,1)*clhs116 + DN(3,2)*clhs139 + clhs225;
            lhs(14,14)=DN(3,0)*clhs81 + DN(3,1)*clhs119 + pow(DN(3,2), 2)*clhs220 + DN(3,2)*clhs141 + clhs219;
            lhs(14,15)=N[3]*clhs226;
            lhs(15,0)=r[0]*(N[0]*clhs227 + N[3]*clhs183);
            lhs(15,1)=r[0]*(N[0]*clhs228 + N[3]*clhs185);
            lhs(15,2)=r[0]*(N[0]*clhs229 + N[3]*clhs187);
            lhs(15,3)=clhs63;
            lhs(15,4)=r[1]*(N[1]*clhs227 + N[3]*clhs144);
            lhs(15,5)=r[1]*(N[1]*clhs228 + N[3]*clhs146);
            lhs(15,6)=r[1]*(N[1]*clhs229 + N[3]*clhs148);
            lhs(15,7)=clhs168;
            lhs(15,8)=r[2]*(N[2]*clhs227 + N[3]*clhs150);
            lhs(15,9)=r[2]*(N[2]*clhs228 + N[3]*clhs151);
            lhs(15,10)=r[2]*(N[2]*clhs229 + N[3]*clhs152);
            lhs(15,11)=clhs199;
            lhs(15,12)=DN(3,0)*clhs230;
            lhs(15,13)=DN(3,1)*clhs230;
            lhs(15,14)=DN(3,2)*clhs230;
            lhs(15,15)=clhs218;


    }

void Stokes3DTwoFluid::ComputeGaussPointRHSContribution(array_1d<double,16>& rhs, const element_data<4,3>& data)
    {
        const int nnodes = 4;
        const int dim = 3;
        
        const double rho = inner_prod(data.N, data.rho);
        const double& bdf0 = data.bdf0;
        const double& bdf1 = data.bdf1;
        const double& bdf2 = data.bdf2;
        
        const bounded_matrix<double,nnodes,dim>& v = data.v;
        const bounded_matrix<double,nnodes,dim>& vn = data.vn;
        const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
        const bounded_matrix<double,nnodes,dim>& f = data.f;
        const array_1d<double,nnodes>& p = data.p;
        
        //get constitutive matrix 
        const Matrix& C = data.C;
        const Vector& stress = data.stress;
        
        //get shape function values
        const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;
        const array_1d<double,nnodes>& N = data.N;
        
        //compute an equivalent tau by Bitrans*c*Bi
        const double tau_denom = replace_tau_denom
        const double tau1 = 1.0/(tau_denom*rho);
        const double tau2 = (C(3,3) + C(4,4) + C(5,5))/(6.0*rho);
        
        //auxiliary variables used in the calculation of the RHS
        const array_1d<double,dim> fgauss = prod(trans(f), N);
        const array_1d<double,dim> vgauss = prod(trans(v), N);
        const array_1d<double,dim> grad_p = prod(trans(DN), p);
        const double pgauss = inner_prod(N,p);
        
        array_1d<double,dim> acch = bdf0*vgauss;
        noalias(acch) += bdf1*prod(trans(vn), N);
        noalias(acch) += bdf2*prod(trans(vnn), N);

        const double crhs0 =             N[0]*p[0] + N[1]*p[1] + N[2]*p[2] + N[3]*p[3];
const double crhs1 =             N[0]*f(0,0);
const double crhs2 =             N[1]*f(1,0);
const double crhs3 =             N[2]*f(2,0);
const double crhs4 =             N[3]*f(3,0);
const double crhs5 =             r[0]*v(0,0);
const double crhs6 =             bdf1*rn[0];
const double crhs7 =             bdf2*rnn[0];
const double crhs8 =             N[0]*(bdf0*crhs5 + crhs6*vn(0,0) + crhs7*vnn(0,0));
const double crhs9 =             r[1]*v(1,0);
const double crhs10 =             bdf1*rn[1];
const double crhs11 =             bdf2*rnn[1];
const double crhs12 =             N[1]*(bdf0*crhs9 + crhs10*vn(1,0) + crhs11*vnn(1,0));
const double crhs13 =             r[2]*v(2,0);
const double crhs14 =             bdf1*rn[2];
const double crhs15 =             bdf2*rnn[2];
const double crhs16 =             N[2]*(bdf0*crhs13 + crhs14*vn(2,0) + crhs15*vnn(2,0));
const double crhs17 =             r[3]*v(3,0);
const double crhs18 =             bdf1*rn[3];
const double crhs19 =             bdf2*rnn[3];
const double crhs20 =             N[3]*(bdf0*crhs17 + crhs18*vn(3,0) + crhs19*vnn(3,0));
const double crhs21 =             crhs1 - crhs12 - crhs16 + crhs2 - crhs20 + crhs3 + crhs4 - crhs8;
const double crhs22 =             r[0]*v(0,1);
const double crhs23 =             r[0]*v(0,2);
const double crhs24 =             r[1]*v(1,1);
const double crhs25 =             r[1]*v(1,2);
const double crhs26 =             r[2]*v(2,1);
const double crhs27 =             r[2]*v(2,2);
const double crhs28 =             r[3]*v(3,1);
const double crhs29 =             r[3]*v(3,2);
const double crhs30 =             N[0]*(bdf0*p[0] + bdf1*pn[0] + bdf2*pnn[0]) + N[1]*(bdf0*p[1] + bdf1*pn[1] + bdf2*pnn[1]) + N[2]*(bdf0*p[2] + bdf1*pn[2] + bdf2*pnn[2]) + N[3]*(bdf0*p[3] + bdf1*pn[3] + bdf2*pnn[3]) + k*(DN(0,0)*crhs5 + DN(0,1)*crhs22 + DN(0,2)*crhs23 + DN(1,0)*crhs9 + DN(1,1)*crhs24 + DN(1,2)*crhs25 + DN(2,0)*crhs13 + DN(2,1)*crhs26 + DN(2,2)*crhs27 + DN(3,0)*crhs17 + DN(3,1)*crhs28 + DN(3,2)*crhs29);
const double crhs31 =             crhs30*tau2;
const double crhs32 =             N[0]*f(0,1);
const double crhs33 =             N[1]*f(1,1);
const double crhs34 =             N[2]*f(2,1);
const double crhs35 =             N[3]*f(3,1);
const double crhs36 =             N[0]*(bdf0*crhs22 + crhs6*vn(0,1) + crhs7*vnn(0,1));
const double crhs37 =             N[1]*(bdf0*crhs24 + crhs10*vn(1,1) + crhs11*vnn(1,1));
const double crhs38 =             N[2]*(bdf0*crhs26 + crhs14*vn(2,1) + crhs15*vnn(2,1));
const double crhs39 =             N[3]*(bdf0*crhs28 + crhs18*vn(3,1) + crhs19*vnn(3,1));
const double crhs40 =             crhs32 + crhs33 + crhs34 + crhs35 - crhs36 - crhs37 - crhs38 - crhs39;
const double crhs41 =             N[0]*f(0,2);
const double crhs42 =             N[1]*f(1,2);
const double crhs43 =             N[2]*f(2,2);
const double crhs44 =             N[3]*f(3,2);
const double crhs45 =             N[0]*(bdf0*crhs23 + crhs6*vn(0,2) + crhs7*vnn(0,2));
const double crhs46 =             N[1]*(bdf0*crhs25 + crhs10*vn(1,2) + crhs11*vnn(1,2));
const double crhs47 =             N[2]*(bdf0*crhs27 + crhs14*vn(2,2) + crhs15*vnn(2,2));
const double crhs48 =             N[3]*(bdf0*crhs29 + crhs18*vn(3,2) + crhs19*vnn(3,2));
const double crhs49 =             crhs41 + crhs42 + crhs43 + crhs44 - crhs45 - crhs46 - crhs47 - crhs48;
const double crhs50 =             tau1*(-crhs1 + crhs12 + crhs16 - crhs2 + crhs20 - crhs3 - crhs4 + crhs8 + grad_p[0]);
const double crhs51 =             tau1*(-crhs32 - crhs33 - crhs34 - crhs35 + crhs36 + crhs37 + crhs38 + crhs39 + grad_p[1]);
const double crhs52 =             tau1*(-crhs41 - crhs42 - crhs43 - crhs44 + crhs45 + crhs46 + crhs47 + crhs48 + grad_p[2]);
            rhs[0]=DN(0,0)*crhs0 - DN(0,0)*crhs31 - DN(0,0)*stress[0] - DN(0,1)*stress[3] - DN(0,2)*stress[5] + N[0]*crhs21;
            rhs[1]=-DN(0,0)*stress[3] + DN(0,1)*crhs0 - DN(0,1)*crhs31 - DN(0,1)*stress[1] - DN(0,2)*stress[4] + N[0]*crhs40;
            rhs[2]=-DN(0,0)*stress[5] - DN(0,1)*stress[4] + DN(0,2)*crhs0 - DN(0,2)*crhs31 - DN(0,2)*stress[2] + N[0]*crhs49;
            rhs[3]=-DN(0,0)*crhs50 - DN(0,1)*crhs51 - DN(0,2)*crhs52 - N[0]*crhs30;
            rhs[4]=DN(1,0)*crhs0 - DN(1,0)*crhs31 - DN(1,0)*stress[0] - DN(1,1)*stress[3] - DN(1,2)*stress[5] + N[1]*crhs21;
            rhs[5]=-DN(1,0)*stress[3] + DN(1,1)*crhs0 - DN(1,1)*crhs31 - DN(1,1)*stress[1] - DN(1,2)*stress[4] + N[1]*crhs40;
            rhs[6]=-DN(1,0)*stress[5] - DN(1,1)*stress[4] + DN(1,2)*crhs0 - DN(1,2)*crhs31 - DN(1,2)*stress[2] + N[1]*crhs49;
            rhs[7]=-DN(1,0)*crhs50 - DN(1,1)*crhs51 - DN(1,2)*crhs52 - N[1]*crhs30;
            rhs[8]=DN(2,0)*crhs0 - DN(2,0)*crhs31 - DN(2,0)*stress[0] - DN(2,1)*stress[3] - DN(2,2)*stress[5] + N[2]*crhs21;
            rhs[9]=-DN(2,0)*stress[3] + DN(2,1)*crhs0 - DN(2,1)*crhs31 - DN(2,1)*stress[1] - DN(2,2)*stress[4] + N[2]*crhs40;
            rhs[10]=-DN(2,0)*stress[5] - DN(2,1)*stress[4] + DN(2,2)*crhs0 - DN(2,2)*crhs31 - DN(2,2)*stress[2] + N[2]*crhs49;
            rhs[11]=-DN(2,0)*crhs50 - DN(2,1)*crhs51 - DN(2,2)*crhs52 - N[2]*crhs30;
            rhs[12]=DN(3,0)*crhs0 - DN(3,0)*crhs31 - DN(3,0)*stress[0] - DN(3,1)*stress[3] - DN(3,2)*stress[5] + N[3]*crhs21;
            rhs[13]=-DN(3,0)*stress[3] + DN(3,1)*crhs0 - DN(3,1)*crhs31 - DN(3,1)*stress[1] - DN(3,2)*stress[4] + N[3]*crhs40;
            rhs[14]=-DN(3,0)*stress[5] - DN(3,1)*stress[4] + DN(3,2)*crhs0 - DN(3,2)*crhs31 - DN(3,2)*stress[2] + N[3]*crhs49;
            rhs[15]=-DN(3,0)*crhs50 - DN(3,1)*crhs51 - DN(3,2)*crhs52 - N[3]*crhs30;

    }
    
                            
void Stokes3DTwoFluid::ComputeGaussPointEnrichmentContributions(
    boost::numeric::ublas::bounded_matrix<double,4,16>& H,
    boost::numeric::ublas::bounded_matrix<double,16,4>& V,
    boost::numeric::ublas::bounded_matrix<double,4,4>&  Kee,
    array_1d<double,4>& rhs_ee,
    const element_data<4,3>& data,
    const array_1d<double,4>& distances,
    const array_1d<double,4>& Nenr,
    const boost::numeric::ublas::bounded_matrix<double,4,4>& DNenr
    )
    {
        const int nnodes = 4;
        const int dim = 3;
        
        const double rho = inner_prod(data.N, data.rho);
        const double& bdf0 = data.bdf0;
        const double& bdf1 = data.bdf1;
        const double& bdf2 = data.bdf2;
        
        const bounded_matrix<double,nnodes,dim>& v = data.v;
        const bounded_matrix<double,nnodes,dim>& vn = data.vn;
        const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
        const bounded_matrix<double,nnodes,dim>& f = data.f;
        const array_1d<double,nnodes>& p = data.p;
        
        //get constitutive matrix 
        const Matrix& C = data.C;
//         const Vector& stress = data.stress;
        
        //get shape function values
        const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;
        const array_1d<double,nnodes>& N = data.N;
        
        //compute an equivalent tau by Bitrans*c*Bi
        const double tau_denom = replace_tau_denom
        const double tau1 = 1.0/(tau_denom*rho);
        
        //auxiliary variables used in the calculation of the RHS
        const array_1d<double,dim> fgauss = prod(trans(f), N);
        const array_1d<double,dim> vgauss = prod(trans(v), N);
        const array_1d<double,dim> grad_p = prod(trans(DN), p);
//         const double pgauss = inner_prod(N,p);
        
        array_1d<double,dim> acch = bdf0*vgauss;
        noalias(acch) += bdf1*prod(trans(vn), N);
        noalias(acch) += bdf2*prod(trans(vnn), N);
        
        array_1d<double,4> penr = ZeroVector(4); //penriched is considered to be zero as we do not want to store it

        
        
    }

}

