//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#include "custom_elements/helmholtz.h"  // See in helmoltz.h: FillElementaData function

namespace Kratos {

template<>
void Helmholtz<3>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int Dim = 3;
    unsigned int NumVar = 1;
    unsigned int NumNodes = Dim+1;
    unsigned int DofSize  = NumNodes*NumVar;

    if (rResult.size() != DofSize)
        rResult.resize(DofSize, false);

    for(unsigned int i=0; i<NumNodes; i++)
    {
        rResult[i*(NumVar)  ]  =  this->GetGeometry()[i].GetDof(ETA).EquationId();
        //~ rResult[i*(NumVar)+1]  =  this->GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
        //~ rResult[i*(NumVar)+2]  =  this->GetGeometry()[i].GetDof(VELOCITY_Z).EquationId();
        //~ rResult[i*(NumVar)+3]  =  this->GetGeometry()[i].GetDof(PRESSURE).EquationId();
    }

    KRATOS_CATCH("")
}


template<>
void Helmholtz<2>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int Dim = 2;
    unsigned int NumVar = 1;
    unsigned int NumNodes = Dim+1;
    unsigned int DofSize  = NumNodes*NumVar;

    if (rResult.size() != DofSize)
        rResult.resize(DofSize, false);

    for(unsigned int i=0; i<NumNodes; i++)
    {
        rResult[i*(NumVar)  ]  =  this->GetGeometry()[i].GetDof(ETA).EquationId();
        //~ rResult[i*(NumVar)+1]  =  this->GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
        //~ rResult[i*(NumVar)+2]  =  this->GetGeometry()[i].GetDof(PRESSURE).EquationId();
    }

    KRATOS_CATCH("")
}


template<>
void Helmholtz<3>::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int Dim = 3;
    unsigned int NumVar = 1;
    unsigned int NumNodes = Dim+1;
    unsigned int DofSize  = NumNodes*NumVar;

    if (ElementalDofList.size() != DofSize)
        ElementalDofList.resize(DofSize);

    for(unsigned int i=0; i<NumNodes; i++)
    {
        ElementalDofList[i*(NumVar)  ]  =  this->GetGeometry()[i].pGetDof(ETA);
        //~ ElementalDofList[i*(NumVar)+1]  =  this->GetGeometry()[i].pGetDof(VELOCITY_Y);
        //~ ElementalDofList[i*(NumVar)+2]  =  this->GetGeometry()[i].pGetDof(VELOCITY_Z);
        //~ ElementalDofList[i*(NumVar)+3]  =  this->GetGeometry()[i].pGetDof(PRESSURE);
    }

    KRATOS_CATCH("");
}


template<>
void Helmholtz<2>::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int Dim = 2;
    unsigned int NumVar = 1;
    unsigned int NumNodes = Dim+1;
    unsigned int DofSize  = NumNodes*NumVar;

    if (ElementalDofList.size() != DofSize)
        ElementalDofList.resize(DofSize);

    for(unsigned int i=0; i<NumNodes; i++)
    {
        ElementalDofList[i*(NumVar)  ]  =  this->GetGeometry()[i].pGetDof(ETA);
        //~ ElementalDofList[i*(NumVar)+1]  =  this->GetGeometry()[i].pGetDof(VELOCITY_Y);
        //~ ElementalDofList[i*(NumVar)+2]  =  this->GetGeometry()[i].pGetDof(PRESSURE);
    }

    KRATOS_CATCH("");
}


//~ template<>
//~ void Helmholtz<3>::ComputeGaussPointLHSContribution(bounded_matrix<double,16,16>& lhs, const ElementDataStruct& data) // NOT USED, SINCE SHALLOW WATER ARE 2D
//~ {
    //~ const int nnodes = 4;
    //~ const int dim = 3;
    //~ const int strain_size = 6;
//~ 
    //~ const double rho = inner_prod(data.N, data.rho);        // Density
    //~ const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    //~ const double h = data.h;                                // Characteristic element size
    //~ const double c = data.c;                                // Wave velocity
//~ 
    //~ const double& dt = data.dt;
    //~ const double& bdf0 = data.bdf0;
    //~ // const double& bdf1 = data.bdf1;
    //~ // const double& bdf2 = data.bdf2;
    //~ const double& dyn_tau = data.dyn_tau;
//~ 
    //~ const bounded_matrix<double,nnodes,dim>& v = data.v;
    //~ // const bounded_matrix<double,nnodes,dim>& vn = data.vn;
    //~ // const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
    //~ const bounded_matrix<double,nnodes,dim>& vmesh = data.vmesh;
    //~ const bounded_matrix<double,nnodes,dim>& vconv = v - vmesh;
    //~ // const bounded_matrix<double,nnodes,dim>& f = data.f;
    //~ // const array_1d<double,nnodes>& p = data.p;
    //~ // const array_1d<double,nnodes>& pn = data.pn;
    //~ // const array_1d<double,nnodes>& pnn = data.pnn;
    //~ const array_1d<double,strain_size>& stress = data.stress;
//~ 
    //~ // Get constitutive matrix
    //~ const Matrix& C = data.C;
//~ 
    //~ // Get shape function values
    //~ const array_1d<double,nnodes>& N = data.N;
    //~ const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;
//~ 
    //~ // const array_1d<double,dim> vconv_gauss = prod(trans(vconv), N);
//~ 
    //~ // const double vconv_norm = norm_2(vconv_gauss);
//~ 
    //~ // Stabilization parameters
    //~ const double stab_c1 = 4.0;
    //~ const double stab_c2 = 2.0;
    //~ // const double tau1 = 1.0/((rho*dyn_tau_coeff)/delta_t + (c2*rho*vconv_norm)/h + (c1*mu)/(h*h));
    //~ // const double tau2 = (h*h)/(c1*tau1);
//~ 
    //~ //substitute_lhs_3D
//~ 
//~ }


template<>
void Helmholtz<2>::ComputeGaussPointLHSContribution(bounded_matrix<double,3,3>& lhs, const ElementDataStruct& data)
{
    const int nnodes = 3;
    const int dim = 2;
    //~ const int strain_size = 3;

    const double H_depth = inner_prod(data.N, data.H_depth);      // Density
    const double h = data.h;                                      // Characteristic element size

    const double& dt = data.dt;
    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;

    //~ const bounded_matrix<double,nnodes,dim>& v = data.v;
    //~ const bounded_matrix<double,nnodes,dim>& vn = data.vn;
    //~ const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
    //~ const bounded_matrix<double,nnodes,dim>& vmesh = data.vmesh;
    //~ const bounded_matrix<double,nnodes,dim>& vconv = v - vmesh;
    const bounded_matrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& eta = data.eta;
    const array_1d<double,nnodes>& etan = data.etan;
    const array_1d<double,nnodes>& etann = data.etann;
    //~ const array_1d<double,strain_size>& stress = data.stress;

    // Get constitutive matrix
    const Matrix& C = data.C;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    // const array_1d<double,dim> vconv_gauss = prod(trans(vconv), N);

    // const double vconv_norm = norm_2(vconv_gauss);

    //~ // Stabilization parameters
    //~ const double stab_c1 = 4.0;
    //~ const double stab_c2 = 2.0;
    //~ // const double tau1 = 1.0/((rho*dyn_tau_coeff)/delta_t + (c2*rho*vconv_norm)/h + (c1*mu)/(h*h));
    //~ // const double tau2 = (h*h)/(c1*tau1);

    const double clhs0 =             H_depth*h;
const double clhs1 =             -clhs0*(DN(0,0)*DN(1,0) + DN(0,1)*DN(1,1));
const double clhs2 =             -clhs0*(DN(0,0)*DN(2,0) + DN(0,1)*DN(2,1));
const double clhs3 =             -clhs0*(DN(1,0)*DN(2,0) + DN(1,1)*DN(2,1));
            lhs(0,0)=-clhs0*(pow(DN(0,0), 2) + pow(DN(0,1), 2));
            lhs(0,1)=clhs1;
            lhs(0,2)=clhs2;
            lhs(1,0)=clhs1;
            lhs(1,1)=-clhs0*(pow(DN(1,0), 2) + pow(DN(1,1), 2));
            lhs(1,2)=clhs3;
            lhs(2,0)=clhs2;
            lhs(2,1)=clhs3;
            lhs(2,2)=-clhs0*(pow(DN(2,0), 2) + pow(DN(2,1), 2));


}


//~ template<>
//~ void NavierStokes<3>::ComputeGaussPointRHSContribution(array_1d<double,16>& rhs, const ElementDataStruct& data) // Not used, since shallow water are 2D
//~ {
    //~ const int nnodes = 4;
    //~ const int dim = 3;
    //~ const int strain_size = 6;
//~ 
    //~ const double rho = inner_prod(data.N, data.rho);        // Density
    //~ const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    //~ const double h = data.h;                                // Characteristic element size
    //~ const double c = data.c;                                // Wave velocity
//~ 
    //~ const double& dt = data.dt;
    //~ const double& bdf0 = data.bdf0;
    //~ const double& bdf1 = data.bdf1;
    //~ const double& bdf2 = data.bdf2;
    //~ const double& dyn_tau = data.dyn_tau;
//~ 
    //~ const bounded_matrix<double,nnodes,dim>& v = data.v;
    //~ const bounded_matrix<double,nnodes,dim>& vn = data.vn;
    //~ const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
    //~ const bounded_matrix<double,nnodes,dim>& vmesh = data.vmesh;
    //~ const bounded_matrix<double,nnodes,dim>& vconv = v - vmesh;
    //~ const bounded_matrix<double,nnodes,dim>& f = data.f;
    //~ const array_1d<double,nnodes>& p = data.p;
    //~ const array_1d<double,nnodes>& pn = data.pn;
    //~ const array_1d<double,nnodes>& pnn = data.pnn;
    //~ const array_1d<double,strain_size>& stress = data.stress;
//~ 
    //~ // Get constitutive matrix
    //~ // const Matrix& C = data.C;
//~ 
    //~ // Get shape function values
    //~ const array_1d<double,nnodes>& N = data.N;
    //~ const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;
//~ 
    //~ // Auxiliary variables used in the calculation of the RHS
    //~ const array_1d<double,dim> f_gauss = prod(trans(f), N);
    //~ const array_1d<double,dim> grad_p = prod(trans(DN), p);
    //~ // const array_1d<double,dim> vconv_gauss = prod(trans(vconv), N);
    //~ const double p_gauss = inner_prod(N,p);
//~ 
    //~ // const double vconv_norm = norm_2(vconv_gauss);
//~ 
    //~ array_1d<double,dim> accel_gauss = bdf0*v_gauss;
    //~ noalias(accel_gauss) += bdf1*prod(trans(vn), N);
    //~ noalias(accel_gauss) += bdf2*prod(trans(vnn), N);
//~ 
    //~ // Stabilization parameters
    //~ const double stab_c1 = 4.0;
    //~ const double stab_c2 = 2.0;
    //~ // const double tau1 = 1.0/((rho*dyn_tau_coeff)/delta_t + (c2*rho*vconv_norm)/h + (c1*mu)/(h*h));
    //~ // const double tau2 = (h*h)/(c1*tau1);
//~ 
    //~ //substitute_rhs_3D
//~ }


template<>
void Helmholtz<2>::ComputeGaussPointRHSContribution(array_1d<double,3>& rhs, const ElementDataStruct& data)
{
    const int nnodes = 3;
    const int dim = 2;
    // const int strain_size = 3;

    const double H_depth = inner_prod(data.N, data.H_depth);      // Density
    const double h = data.h;                                      // Characteristic element size

    const double& dt = data.dt;
    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;

    //~ const bounded_matrix<double,nnodes,dim>& v = data.v;
    //~ const bounded_matrix<double,nnodes,dim>& vn = data.vn;
    //~ const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
    //~ const bounded_matrix<double,nnodes,dim>& vmesh = data.vmesh;
    //~ const bounded_matrix<double,nnodes,dim>& vconv = v - vmesh;
    const bounded_matrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& eta = data.eta;
    const array_1d<double,nnodes>& etan = data.etan;
    const array_1d<double,nnodes>& etann = data.etann;
    //~ const array_1d<double,strain_size>& stress = data.stress;

    // Get constitutive matrix
    // const Matrix& C = data.C;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,dim> f_gauss = prod(trans(f), N);
    const array_1d<double,dim> grad_eta = prod(trans(DN), eta);
    // const array_1d<double,dim> vconv_gauss = prod(trans(vconv), N);
    //~ const double p_gauss = inner_prod(N,p);

    // const double vconv_norm = norm_2(vconv_gauss);

    //~ array_1d<double,dim> accel_gauss = bdf0*v_gauss;
    //~ noalias(accel_gauss) += bdf1*prod(trans(vn), N);
    //~ noalias(accel_gauss) += bdf2*prod(trans(vnn), N);

    //~ // Stabilization parameters
    //~ const double stab_c1 = 4.0;
    //~ const double stab_c2 = 2.0;
    //~ // const double tau1 = 1.0/((rho*dyn_tau_coeff)/delta_t + (c2*rho*vconv_norm)/h + (c1*mu)/(h*h));
    //~ // const double tau2 = (h*h)/(c1*tau1);

    const double crhs0 =             H_depth*h;
const double crhs1 =             DN(0,0)*eta[0] + DN(1,0)*eta[1] + DN(2,0)*eta[2];
const double crhs2 =             DN(0,1)*eta[0] + DN(1,1)*eta[1] + DN(2,1)*eta[2];
            rhs[0]=crhs0*(DN(0,0)*crhs1 + DN(0,1)*crhs2);
            rhs[1]=crhs0*(DN(1,0)*crhs1 + DN(1,1)*crhs2);
            rhs[2]=crhs0*(DN(2,0)*crhs1 + DN(2,1)*crhs2);

}


//~ template<>
//~ double NavierStokes<3>::SubscaleErrorEstimate(const ElementDataStruct& data) // Not used, since shallow water are 2D
//~ {
    //~ const int nnodes = 4;
    //~ const int dim = 3;
    //~ // const int strain_size = 3;
//~ 
    //~ const double rho = inner_prod(data.N, data.rho);        // Density
    //~ const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    //~ const double h = data.h;                                // Characteristic element size
    //~ // const double c = data.c;                                // Wave velocity
//~ 
    //~ const double& dt = data.dt;
    //~ const double& bdf0 = data.bdf0;
    //~ const double& bdf1 = data.bdf1;
    //~ const double& bdf2 = data.bdf2;
    //~ const double& dyn_tau = data.dyn_tau;
//~ 
    //~ const bounded_matrix<double,nnodes,dim>& v = data.v;
    //~ const bounded_matrix<double,nnodes,dim>& vn = data.vn;
    //~ const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
    //~ const bounded_matrix<double,nnodes,dim>& vmesh = data.vmesh;
    //~ const bounded_matrix<double,nnodes,dim>& vconv = v - vmesh;
    //~ const bounded_matrix<double,nnodes,dim>& f = data.f;
    //~ const array_1d<double,nnodes>& p = data.p;
    //~ // const array_1d<double,nnodes>& pn = data.pn;
    //~ // const array_1d<double,nnodes>& pnn = data.pnn;
    //~ // const array_1d<double,strain_size>& stress = data.stress;
//~ 
    //~ // // Get constitutive matrix
    //~ // const Matrix& C = data.C;
//~ 
    //~ // Get shape function values
    //~ const array_1d<double,nnodes>& N = data.N;
    //~ const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;
//~ 
    //~ // Auxiliary variables used in the calculation of the error estimator
    //~ array_1d<double,dim> v_s_gauss;
    //~ const array_1d<double,dim> v_gauss = prod(trans(v), N);
    //~ const array_1d<double,dim> f_gauss = prod(trans(f), N);
    //~ const array_1d<double,dim> grad_p = prod(trans(DN), p);
    //~ // const array_1d<double,dim> vconv_gauss = prod(trans(vconv), N);
//~ 
    //~ // const double vconv_norm = norm_2(vconv_gauss);
//~ 
    //~ // Stabilization parameters
    //~ const double stab_c1 = 4.0;
    //~ const double stab_c2 = 2.0;
    //~ // const double tau1 = 1.0/((rho*dyn_tau_coeff)/delta_t + (c2*rho*vconv_norm)/h + (c1*mu)/(h*h));
    //~ // // const double tau2 = (h*h)/(c1*tau1);
//~ 
    //~ // Gauss point velocity subscale value computation
    //~ //substitute_gausspt_subscale_3D
//~ 
    //~ const double v_gauss_norm = norm_2(v_gauss);
    //~ const double v_s_gauss_norm = norm_2(v_s_gauss);
//~ 
    //~ return v_s_gauss_norm/v_gauss_norm;
//~ }


template<>
double Helmholtz<2>::SubscaleErrorEstimate(const ElementDataStruct& data)  // Not finished
{
    const int nnodes = 3;
    const int dim = 2;

    const double h_depth = inner_prod(data.N, data.H_depth);     // Bathymetry
    const double h = data.h;                                     // Characteristic element size

    const double& dt = data.dt;
    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& dyn_tau = data.dyn_tau;

    //~ const bounded_matrix<double,nnodes,dim>& v = data.v;
    //~ const bounded_matrix<double,nnodes,dim>& vn = data.vn;
    //~ const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
    //~ const bounded_matrix<double,nnodes,dim>& vmesh = data.vmesh;
    //~ const bounded_matrix<double,nnodes,dim>& vconv = v - vmesh;
    const bounded_matrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& eta = data.eta;
    const array_1d<double,nnodes>& etan = data.etan;
    const array_1d<double,nnodes>& etann = data.etann;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the error estimator
    // array_1d<double,dim> v_s_gauss;
    // const array_1d<double,dim> v_gauss = prod(trans(v), N);
    const array_1d<double,dim> f_gauss = prod(trans(f), N);
    const array_1d<double,dim> grad_eta = prod(trans(DN), eta);
    // const array_1d<double,dim> vconv_gauss = prod(trans(vconv), N);

    // const double vconv_norm = norm_2(vconv_gauss);

    //~ // Stabilization parameters
    //~ const double stab_c1 = 4.0;
    //~ const double stab_c2 = 2.0;
    //~ // const double tau1 = 1.0/((rho*dyn_tau_coeff)/delta_t + (c2*rho*vconv_norm)/h + (c1*mu)/(h*h));
    //~ // // const double tau2 = (h*h)/(c1*tau1);

    // Gauss point velocity subscale value computation
    //substitute_gausspt_subscale_2D

    // const double v_gauss_norm = norm_2(v_gauss);
    // const double v_s_gauss_norm = norm_2(v_s_gauss);

    return v_s_gauss_norm/v_gauss_norm;
}

}
