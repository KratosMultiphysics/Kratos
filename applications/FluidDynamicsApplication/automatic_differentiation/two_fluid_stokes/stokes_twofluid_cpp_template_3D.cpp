#include "custom_elements/stokes_3D_twofluid.h"

namespace Kratos {


void Stokes3DTwoFluid::ComputeGaussPointLHSContribution(BoundedMatrix<double,16,16>& lhs, const element_data<4,3>& data)
    {
        const int nnodes = 4;
        const int dim = 3;

        const double rho = inner_prod(data.N, data.rho);
        const double& bdf0 = data.bdf0;

        //get constitutive matrix
        const Matrix& C = data.C;

        //get shape function values
        const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;
        const array_1d<double,nnodes>& N = data.N;


        //compute an equivalent tau by Bitrans*c*Bi
        const double tau_denom = replace_tau_denom
        const double tau1 = 1.0/(tau_denom*rho);
        const double tau2 = (C(3,3) + C(4,4) + C(5,5))/(6.0*rho);

        const double DN_0_0 = DN(0,0); const double DN_0_1 = DN(0,1); const double DN_0_2 = DN(0,2);
        const double DN_1_0 = DN(1,0); const double DN_1_1 = DN(1,1); const double DN_1_2 = DN(1,2);
        const double DN_2_0 = DN(2,0); const double DN_2_1 = DN(2,1); const double DN_2_2 = DN(2,2);
        const double DN_3_0 = DN(3,0); const double DN_3_1 = DN(3,1); const double DN_3_2 = DN(3,2);

        const double N_0 = N[0];
        const double N_1 = N[1];
        const double N_2 = N[2];
        const double N_3 = N[3];

        const double C_0_0 = C(0,0); const double C_0_1 = C(0,1); const double C_0_2 = C(0,2); const double C_0_3 = C(0,3); const double C_0_4 = C(0,4); const double C_0_5 = C(0,5);
        const double C_1_1 = C(1,1); const double C_1_2 = C(1,2); const double C_1_3 = C(1,3); const double C_1_4 = C(1,4); const double C_1_5 = C(1,5);
        const double C_2_2 = C(2,2); const double C_2_3 = C(2,3); const double C_2_4 = C(2,4); const double C_2_5 = C(2,5);
        const double C_3_3 = C(3,3); const double C_3_4 = C(3,4); const double C_3_5 = C(3,5);
        const double C_4_4 = C(4,4); const double C_4_5 = C(4,5); const double C_5_5 = C(5,5);

        //substitute_lhs

    }

void Stokes3DTwoFluid::ComputeGaussPointRHSContribution(array_1d<double,16>& rhs, const element_data<4,3>& data)
    {
        const int nnodes = 4;
        const int dim = 3;

        const double rho = inner_prod(data.N, data.rho);
        const double& bdf0 = data.bdf0;
        const double& bdf1 = data.bdf1;
        const double& bdf2 = data.bdf2;

        const BoundedMatrix<double,nnodes,dim>& v = data.v;
        const BoundedMatrix<double,nnodes,dim>& vn = data.vn;
        const BoundedMatrix<double,nnodes,dim>& vnn = data.vnn;
        const BoundedMatrix<double,nnodes,dim>& f = data.f;
        const array_1d<double,nnodes>& p = data.p;

        //get constitutive matrix
        const Matrix& C = data.C;
        const Vector& stress = data.stress;

        //get shape function values
        const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;
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

        const double fgauss_0 = fgauss[0];
        const double fgauss_1 = fgauss[1];
        const double fgauss_2 = fgauss[2];
        
        const double grad_p_0 = grad_p[0];
        const double grad_p_1 = grad_p[1];
        const double grad_p_2 = grad_p[2];

        const double stress_0 = stress[0];
        const double stress_1 = stress[1];
        const double stress_2 = stress[2];
        const double stress_3 = stress[3];
        const double stress_4 = stress[4];
        const double stress_5 = stress[5];

        const double N_0 = N[0];
        const double N_1 = N[1];
        const double N_2 = N[2];
        const double N_3 = N[3];

        const double f_0_0 = f(0,0); const double f_0_1 = f(0,1); const double f_0_2 = f(0,2);
        const double f_1_0 = f(0,0); const double f_1_1 = f(0,1); const double f_1_2 = f(0,2);
        const double f_2_0 = f(0,0); const double f_2_1 = f(0,1); const double f_2_2 = f(0,2);
        const double f_3_0 = f(0,0); const double f_3_1 = f(0,1); const double f_3_2 = f(0,2);

        const double vn_0_0 = vn(0,0); const double vn_0_1 = vn(0,1); const double vn_0_2 = vn(0,2);
        const double vn_1_0 = vn(1,0); const double vn_1_1 = vn(1,1); const double vn_1_2 = vn(1,2);
        const double vn_2_0 = vn(2,0); const double vn_2_1 = vn(2,1); const double vn_2_2 = vn(2,2);
        const double vn_3_0 = vn(3,0); const double vn_3_1 = vn(3,1); const double vn_3_2 = vn(3,2);

        const double vnn_0_0 = vnn(0,0); const double vnn_0_1 = vnn(0,1); const double vnn_0_2 = vnn(0,2);
        const double vnn_1_0 = vnn(1,0); const double vnn_1_1 = vnn(1,1); const double vnn_1_2 = vnn(1,2);
        const double vnn_2_0 = vnn(2,0); const double vnn_2_1 = vnn(2,1); const double vnn_2_2 = vnn(2,2);
        const double vnn_3_0 = vnn(3,0); const double vnn_3_1 = vnn(3,1); const double vnn_3_2 = vnn(3,2);

        const double acch_0 = acch[0];
        const double acch_1 = acch[1];
        const double acch_2 = acch[2];

        const double C_0_0 = C(0,0); const double C_0_1 = C(0,1); const double C_0_2 = C(0,2); const double C_0_3 = C(0,3); const double C_0_4 = C(0,4); const double C_0_5 = C(0,5);
        const double C_1_1 = C(1,1); const double C_1_2 = C(1,2); const double C_1_3 = C(1,3); const double C_1_4 = C(1,4); const double C_1_5 = C(1,5);
        const double C_2_2 = C(2,2); const double C_2_3 = C(2,3); const double C_2_4 = C(2,4); const double C_2_5 = C(2,5);
        const double C_3_3 = C(3,3); const double C_3_4 = C(3,4); const double C_3_5 = C(3,5);
        const double C_4_4 = C(4,4); const double C_4_5 = C(4,5); const double C_5_5 = C(5,5);
        
        const double DN_0_0 = DN(0,0); const double DN_0_1 = DN(0,1); const double DN_0_2 = DN(0,2);
        const double DN_1_0 = DN(1,0); const double DN_1_1 = DN(1,1); const double DN_1_2 = DN(1,2);
        const double DN_2_0 = DN(2,0); const double DN_2_1 = DN(2,1); const double DN_2_2 = DN(2,2);
        const double DN_3_0 = DN(3,0); const double DN_3_1 = DN(3,1); const double DN_3_2 = DN(3,2);
        
        const double v_0_0 = v(0,0); const double v_0_1 = v(0,1); const double v_0_2 = v(0,2);
        const double v_1_0 = v(1,0); const double v_1_1 = v(1,1); const double v_1_2 = v(1,2);
        const double v_2_0 = v(2,0); const double v_2_1 = v(2,1); const double v_2_2 = v(2,2);
        const double v_3_0 = v(3,0); const double v_3_1 = v(3,1); const double v_3_2 = v(3,2);

        //substitute_rhs
    }


void Stokes3DTwoFluid::ComputeGaussPointEnrichmentContributions(
    BoundedMatrix<double,4,16>& H,
    BoundedMatrix<double,16,4>& V,
    BoundedMatrix<double,4,4>&  Kee,
    array_1d<double,4>& rhs_ee,
    const element_data<4,3>& data,
    const array_1d<double,4>& distances,
    const array_1d<double,4>& Nenr,
    const BoundedMatrix<double,4,4>& DNenr
    )
    {
        const int nnodes = 4;
        const int dim = 3;

        const double rho = inner_prod(data.N, data.rho);
        const double& bdf0 = data.bdf0;
        const double& bdf1 = data.bdf1;
        const double& bdf2 = data.bdf2;

        const BoundedMatrix<double,nnodes,dim>& v = data.v;
        const BoundedMatrix<double,nnodes,dim>& vn = data.vn;
        const BoundedMatrix<double,nnodes,dim>& vnn = data.vnn;
        const BoundedMatrix<double,nnodes,dim>& f = data.f;
        const array_1d<double,nnodes>& p = data.p;

        //get constitutive matrix
        const Matrix& C = data.C;
//         const Vector& stress = data.stress;

        //get shape function values
        const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;
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

        const double DN_0_0 = DN(0,0); const double DN_0_1 = DN(0,1); const double DN_0_2 = DN(0,2);
        const double DN_1_0 = DN(1,0); const double DN_1_1 = DN(1,1); const double DN_1_2 = DN(1,2);
        const double DN_2_0 = DN(2,0); const double DN_2_1 = DN(2,1); const double DN_2_2 = DN(2,2);
        const double DN_3_0 = DN(3,0); const double DN_3_1 = DN(3,1); const double DN_3_2 = DN(3,2);

        const double Nenr_0 = Nenr[0];
        const double Nenr_1 = Nenr[1];
        const double Nenr_2 = Nenr[2];
        const double Nenr_3 = Nenr[3];

        const double N_0 = N[0];
        const double N_1 = N[1];
        const double N_2 = N[2];
        const double N_3 = N[3];

        const double p_0 = p[0];
        const double p_1 = p[1];
        const double p_2 = p[2];
        const double p_3 = p[3];
        
        const double penr_0 = penr[0];
        const double penr_1 = penr[1];
        const double penr_2 = penr[2];
        const double penr_3 = penr[3];

        const double DNenr_0_0 = DNenr(0,0); const double DNenr_0_1 = DNenr(0,1); const double DNenr_0_2 = DNenr(0,2);
        const double DNenr_1_0 = DNenr(1,0); const double DNenr_1_1 = DNenr(1,1); const double DNenr_1_2 = DNenr(1,2);
        const double DNenr_2_0 = DNenr(2,0); const double DNenr_2_1 = DNenr(2,1); const double DNenr_2_2 = DNenr(2,2);
        const double DNenr_3_0 = DNenr(3,0); const double DNenr_3_1 = DNenr(3,1); const double DNenr_3_2 = DNenr(3,2);

        const double f_0_0 = f(0,0); const double f_0_1 = f(0,1); const double f_0_2 = f(0,2);
        const double f_1_0 = f(1,0); const double f_1_1 = f(1,1); const double f_1_2 = f(1,2);
        const double f_2_0 = f(2,0); const double f_2_1 = f(2,1); const double f_2_2 = f(2,2);
        const double f_3_0 = f(3,0); const double f_3_1 = f(3,1); const double f_3_2 = f(3,2);

        const double v_0_0 = v(0,0); const double v_0_1 = v(0,1); const double v_0_2 = v(0,2); 
        const double v_1_0 = v(1,0); const double v_1_1 = v(1,1); const double v_1_2 = v(1,2);
        const double v_2_0 = v(2,0); const double v_2_1 = v(2,1); const double v_2_2 = v(2,2);
        const double v_3_0 = v(3,0); const double v_3_1 = v(3,1); const double v_3_2 = v(3,2);

        const double vn_0_0 = vn(0,0); const double vn_0_1 = vn(0,1); const double vn_0_2 = vn(0,2); 
        const double vn_1_0 = vn(1,0); const double vn_1_1 = vn(1,1); const double vn_1_2 = vn(1,2);
        const double vn_2_0 = vn(2,0); const double vn_2_1 = vn(2,1); const double vn_2_2 = vn(2,2);
        const double vn_3_0 = vn(3,0); const double vn_3_1 = vn(3,1); const double vn_3_2 = vn(3,2);

        const double vnn_0_0 = vnn(0,0); const double vnn_0_1 = vnn(0,1); const double vnn_0_2 = vnn(0,2); 
        const double vnn_1_0 = vnn(1,0); const double vnn_1_1 = vnn(1,1); const double vnn_1_2 = vnn(1,2);
        const double vnn_2_0 = vnn(2,0); const double vnn_2_1 = vnn(2,1); const double vnn_2_2 = vnn(2,2);
        const double vnn_3_0 = vnn(3,0); const double vnn_3_1 = vnn(3,1); const double vnn_3_2 = vnn(3,2);

        //substitute_enrichment_V

        //substitute_enrichment_H

        //substitute_enrichment_Kee

        //substitute_enrichment_rhs_ee


    }

}