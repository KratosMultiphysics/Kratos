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
        
        //substitute_lhs_3D

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
        
        //substitute_lhs_2D

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
        
        //substitute_rhs_3D
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
        
        //substitute_rhs_2D
    }
    
}

