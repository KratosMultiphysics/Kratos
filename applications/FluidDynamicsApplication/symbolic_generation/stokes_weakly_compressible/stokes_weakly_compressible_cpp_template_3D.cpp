#include "custom_elements/stokes_3D_weakly_compressible.h"

namespace Kratos {


void Stokes3DWeaklyCompressible::ComputeGaussPointLHSContribution(bounded_matrix<double,16,16>& lhs, const element_data<4,3>& data)
    {
        const int nnodes = 4;
        const int dim = 3;
        
        const double& bdf0 = data.bdf0;
//         const double& bdf1 = data.bdf1;
//         const double& bdf2 = data.bdf2;
        
        const auto& v = data.v;
//         const auto& vn = data.vn;
//         const auto& vnn = data.vnn;
        
//         const auto& p = data.p;
//         const auto& pn = data.pn;
//         const auto& pnn = data.pnn;
        
         const auto& r = data.r;
//         const auto& rn = data.rn;
//         const auto& rnn = data.rnn;   
        
        const double rho_gauss = inner_prod(data.N, data.r);
         
        const double k = data.k;
        
        //get constitutive matrix 
        const Matrix& C = data.C;
        
        //get shape function values
        const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;
        const array_1d<double,nnodes>& N = data.N;
        
        
        //compute an equivalent tau by Bitrans*c*Bi
        const double tau_denom = //replace_tau_denom
        const double tau1 = 1.0/(tau_denom*rho_gauss);
//         const double tau2 = (C(3,3) + C(4,4) + C(5,5))/(6.0*rho_gauss);
        
        //substitute_lhs

    }

void Stokes3DWeaklyCompressible::ComputeGaussPointRHSContribution(array_1d<double,16>& rhs, const element_data<4,3>& data)
    {
        const int nnodes = 4;
        const int dim = 3;
        
        
        const double& bdf0 = data.bdf0;
        const double& bdf1 = data.bdf1;
        const double& bdf2 = data.bdf2;
        
        const bounded_matrix<double,nnodes,dim>& v = data.v;
        const bounded_matrix<double,nnodes,dim>& vn = data.vn;
        const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
        const bounded_matrix<double,nnodes,dim>& f = data.f;
        const auto& p = data.p;
        const auto& pn = data.pn;
        const auto& pnn = data.pnn;
        const auto& r = data.r;
        const auto& rn = data.rn;
        const auto& rnn = data.rnn; 
        
        const double rho_gauss = inner_prod(data.N, data.r);
        
        const double k = data.k;
        
        //get constitutive matrix 
        const Matrix& C = data.C;
        const Vector& stress = data.stress;
        
        //get shape function values
        const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;
        const array_1d<double,nnodes>& N = data.N;
        
        //compute an equivalent tau by Bitrans*c*Bi
        const double tau_denom = //replace_tau_denom
        const double tau1 = 1.0/(tau_denom*rho_gauss);
//         const double tau2 = (C(3,3) + C(4,4) + C(5,5))/(6.0*rho_gauss);
        
        

        //substitute_rhs
        
    }
   
}
