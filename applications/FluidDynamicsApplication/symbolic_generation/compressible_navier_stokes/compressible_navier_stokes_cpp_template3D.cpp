//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Elisa Magliozzi
//

#include "custom_elements/compressible_navier_stokes.h"

namespace Kratos {

template<>
void CompressibleNavierStokes<3>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int Dim = 3;
    unsigned int BlockSize = Dim+2;
    unsigned int NumNodes = 4;
    unsigned int DofSize  = NumNodes*(BlockSize);

    if (rResult.size() != DofSize)
        rResult.resize(DofSize, false);

    for(unsigned int i=0; i<NumNodes; i++)
    {
        rResult[i*(BlockSize)  ]  =  this->GetGeometry()[i].GetDof(DENSITY).EquationId();
        rResult[i*(BlockSize)+1]  =  this->GetGeometry()[i].GetDof(MOMENTUM_X).EquationId();
        rResult[i*(BlockSize)+2]  =  this->GetGeometry()[i].GetDof(MOMENTUM_Y).EquationId();
        rResult[i*(BlockSize)+3]  =  this->GetGeometry()[i].GetDof(MOMENTUM_Z).EquationId(); 
        rResult[i*(BlockSize)+4]  =  this->GetGeometry()[i].GetDof(TOTAL_ENERGY).EquationId();
    }

    KRATOS_CATCH("")
}

template<>
void CompressibleNavierStokes<3>::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    unsigned int Dim = 3;
    unsigned int BlockSize = Dim+2;
    unsigned int NumNodes = 4;
    unsigned int DofSize  = NumNodes*(BlockSize);

    if (ElementalDofList.size() != DofSize)
        ElementalDofList.resize(DofSize);

    for(unsigned int i=0; i<NumNodes; i++)
    {
        ElementalDofList[i*(BlockSize)  ]  =  this->GetGeometry()[i].pGetDof(DENSITY);
        ElementalDofList[i*(BlockSize)+1]  =  this->GetGeometry()[i].pGetDof(MOMENTUM_X);
        ElementalDofList[i*(BlockSize)+2]  =  this->GetGeometry()[i].pGetDof(MOMENTUM_Y);
        ElementalDofList[i*(BlockSize)+3]  =  this->GetGeometry()[i].pGetDof(MOMENTUM_Z);
        ElementalDofList[i*(BlockSize)+4]  =  this->GetGeometry()[i].pGetDof(TOTAL_ENERGY);
    }

    KRATOS_CATCH("");
}

template<>
void CompressibleNavierStokes<3>::ComputeGaussPointLHSContribution(BoundedMatrix<double,20,20>& lhs, const ElementDataStruct& data,double data_v_sc, double data_k_sc)
{
    const int nnodes = 4;
    const int dim = 3;
    const int BlockSize = dim+2;
    const double h = data.h; 
 
    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
 
    const BoundedMatrix<double,nnodes,BlockSize>& U = data.U;
    const BoundedMatrix<double,nnodes,BlockSize>& Un = data.Un;
    const BoundedMatrix<double,nnodes,BlockSize>& Unn = data.Unn;
    const BoundedMatrix<double,nnodes,dim>& f_ext = data.f_ext;
    const array_1d<double,nnodes>& r = data.r;
    const double mu = data.mu;
    const double nu = data.nu;
    const double lambda = data.lambda;
    const double c_v = data.c_v;
    const double gamma = data.gamma;
    const double cp = c_v*gamma;
    const double v_sc = data_v_sc;
    const double k_sc = data_k_sc;
 
    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;
    
    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;
    
    const array_1d<double,BlockSize> U_gauss= prod(trans(U),N);
    
    double tmp = U_gauss(dim+1)/U_gauss(0);
    for(unsigned int ll=0; ll<dim; ll++)
        tmp -=(U_gauss(ll+1)*U_gauss(ll+1))/(2*U_gauss(0)*U_gauss(0));
    double c = sqrt(gamma*(gamma-1)*tmp);

    double tau1inv = 0.0;
    for(unsigned int ll=0; ll<dim; ll++)
        tau1inv += (U_gauss(ll+1)/U_gauss(0))*(U_gauss(ll+1)/U_gauss(0));
    tau1inv = (sqrt(tau1inv)+c)*stab_c2/h;
    double tau2inv = stab_c1*nu/(h*h)+tau1inv;
    double tau3inv = stab_c1*lambda/(U_gauss(0)*cp*h*h)+tau1inv;
        
    const double tau1 = 1/tau1inv;
    const double tau2 = 1/tau2inv;
    const double tau3 = 1/tau3inv;
    
    //substitute_lhs_3D

}

template<>
void CompressibleNavierStokes<3>::ComputeGaussPointRHSContribution(array_1d<double,20>& rhs, const ElementDataStruct& data,double data_v_sc, double data_k_sc)
{
    const int nnodes = 4;
    const int dim = 3;
    const int BlockSize = dim+2;
    const double h = data.h; 

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    
    const BoundedMatrix<double,nnodes,BlockSize>& U = data.U;
    const BoundedMatrix<double,nnodes,BlockSize>& Un = data.Un;
    const BoundedMatrix<double,nnodes,BlockSize>& Unn = data.Unn;
    const BoundedMatrix<double,nnodes,dim>& f_ext = data.f_ext;
    const array_1d<double,nnodes>& r = data.r;
    const double mu = data.mu;
    const double nu = data.nu;
    const double lambda = data.lambda;
    const double c_v = data.c_v;
    const double gamma = data.gamma;
    const double cp = c_v*gamma;
    const double v_sc = data_v_sc;
    const double k_sc = data_k_sc;
    
    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,BlockSize> U_gauss = prod(trans(U), N);
    const array_1d<double,dim> f_gauss = prod(trans(f_ext), N);
    const BoundedMatrix<double,dim,BlockSize> grad_U = prod(trans(U), DN); 
    const array_1d<double,BlockSize> accel_gauss = bdf0*U_gauss+bdf1*prod(trans(Un), N)+bdf2*prod(trans(Unn), N);
    
    // Stabilization parameters
    const double stab_c1 = 4.0;
    const double stab_c2 = 2.0;
    double tmp = U_gauss(dim+1)/U_gauss(0);
    for(unsigned int ll=0; ll<dim; ll++)
        tmp -=(U_gauss(ll+1)*U_gauss(ll+1))/(2*U_gauss(0)*U_gauss(0));
    double c = sqrt(gamma*(gamma-1)*tmp);

    double tau1inv = 0.0;
    for(unsigned int ll=0; ll<dim; ll++)
        tau1inv += (U_gauss(ll+1)/U_gauss(0))*(U_gauss(ll+1)/U_gauss(0));
    tau1inv = (sqrt(tau1inv)+c)*stab_c2/h;
    double tau2inv = stab_c1*nu/(h*h)+tau1inv;
    double tau3inv = stab_c1*lambda/(U_gauss(0)*cp*h*h)+tau1inv;
        
    const double tau1 = 1/tau1inv;
    const double tau2 = 1/tau2inv;
    const double tau3 = 1/tau3inv;

    //substitute_rhs_3D
}

template<>
double CompressibleNavierStokes<3>::ShockCapturingViscosity(const ElementDataStruct& data)
{
    const int nnodes = 4;
    const int dim = 3;
    const int BlockSize = dim+2;

   const double h = data.h;                                // Characteristic element size
   const double alpha = 0.8;                               // Algorithm constant
   const double tol = 0.001;                               

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;

    const BoundedMatrix<double,nnodes,BlockSize>& U = data.U;
    const BoundedMatrix<double,nnodes,BlockSize>& Un = data.Un;
    const BoundedMatrix<double,nnodes,BlockSize>& Unn = data.Unn;
    const BoundedMatrix<double,nnodes,dim>& f_ext = data.f_ext;
    const double gamma = data.gamma;
    double v_sc = 0.0;                                      //Shock capturing viscosity
    BoundedMatrix<double,dim,1> res_m; 
    res_m(0,0)= 0; res_m(1,0)= 0; res_m(2,0)= 0;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,BlockSize> U_gauss = prod(trans(U), N);
    const array_1d<double,dim> f_gauss = prod(trans(f_ext), N);
    const BoundedMatrix<double,BlockSize,dim> grad_U = prod(trans(U), DN);     // Dfi/Dxj
    const array_1d<double,BlockSize> accel_gauss = bdf0*U_gauss+bdf1*prod(trans(Un), N)+bdf2*prod(trans(Unn), N);
    
    //substitute_res_m_3D

    double norm_res_m;
    norm_res_m = sqrt(res_m(0,0)*res_m(0,0)+res_m(1,0)*res_m(1,0)+res_m(2,0)*res_m(2,0));

    double norm_gradm = 0.0;                                    // Frobenius norm of momentum gradient
    for (unsigned int i=1; i<dim+1; i++){
        for (unsigned int j=0; j<dim; j++)
            norm_gradm += grad_U(i,j)*grad_U(i,j);
    }
    norm_gradm = sqrt(norm_gradm);
    
    if (norm_gradm>tol)
        v_sc = 0.5*h*alpha*(norm_res_m/norm_gradm);
    
    return v_sc;
}

template<>
double CompressibleNavierStokes<3>::ShockCapturingConductivity(const ElementDataStruct& data)
{
    const int nnodes = 4;
    const int dim = 3;
    const int BlockSize = dim+2;

   const double h = data.h;                                // Characteristic element size
   const double alpha = 0.8;                               // Algorithm constant
   const double tol = 0.001;                               

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;

    const BoundedMatrix<double,nnodes,BlockSize>& U = data.U;
    const BoundedMatrix<double,nnodes,BlockSize>& Un = data.Un;
    const BoundedMatrix<double,nnodes,BlockSize>& Unn = data.Unn;
    const BoundedMatrix<double,nnodes,dim>& f_ext = data.f_ext;
    const array_1d<double,nnodes>& r = data.r;
    const double gamma = data.gamma;
    double k_sc = 0.0;          // Shock Capturing Conductivity
    BoundedMatrix<double,dim,1> res_e;
    res_e(0,0) = 0;
    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,BlockSize> U_gauss = prod(trans(U), N);
    const array_1d<double,dim> f_gauss = prod(trans(f_ext), N);
    const BoundedMatrix<double,BlockSize,dim> grad_U = prod(trans(U), DN);     // Dfi/Dxj
    const array_1d<double,BlockSize> accel_gauss = bdf0*U_gauss+bdf1*prod(trans(Un), N)+bdf2*prod(trans(Unn), N);
    
    //substitute_res_e_3D

    double norm_res_e;
    norm_res_e = sqrt(res_e(0,0)*res_e(0,0));

    double norm_grade = 0.0;              // Frobenius norm of total energy gradient
    for (unsigned int i=0; i<dim; i++)
        norm_grade += grad_U(dim+1,i)*grad_U(dim+1,i);
    norm_grade = sqrt(norm_grade);
    
 
    if (norm_grade>tol)
        k_sc = 0.5*h*alpha*(norm_res_e/norm_grade);

    return k_sc;
}

}
