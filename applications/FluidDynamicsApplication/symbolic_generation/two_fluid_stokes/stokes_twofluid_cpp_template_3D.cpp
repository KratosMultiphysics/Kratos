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

        //substitute_enrichment_V

        //substitute_enrichment_H

        //substitute_enrichment_Kee

        //substitute_enrichment_rhs_ee


    }

}

