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
        //~ rResult[i*(NumVar)  ]  =  this->GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
        //~ rResult[i*(NumVar)+1]  =  this->GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
        rResult[i*NumVar  ]  =  this->GetGeometry()[i].GetDof(ELEVATION).EquationId();
    }

    KRATOS_CATCH("")
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
        ElementalDofList[i*(NumVar)  ]  =  this->GetGeometry()[i].pGetDof(ELEVATION);
        //~ ElementalDofList[i*(NumVar)+1]  =  this->GetGeometry()[i].pGetDof(VELOCITY_Y);
        //~ ElementalDofList[i*(NumVar)+2]  =  this->GetGeometry()[i].pGetDof(PRESSURE);
    }

    KRATOS_CATCH("");
}


template<>
void Helmholtz<2>::ComputeGaussPointLHSContribution(bounded_matrix<double,3,3>& lhs, const ElementDataStruct& data)
{
    const int nnodes = 3;
    const int dim = 2;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    // Get problem data
    const double H = inner_prod(N, data.H);                // Bathymetry
    const double g = data.g;                               // Gravity
    //~ const double l = data.l;                               // Characteristic element size  // Unused variable

    //~ const double& dt = data.dt;              // Unused variable
    const double& bdf0 = data.bdf0;          // Unused variable
    //~ const double& bdf1 = data.bdf1;          // Unused variable
    //~ const double& bdf2 = data.bdf2;          // Unused variable

    //~ const bounded_matrix<double,nnodes,dim>& v = data.v;
    //~ const bounded_matrix<double,nnodes,dim>& vn = data.vn;
    //~ const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
    //~ const bounded_matrix<double,nnodes,dim>& vmesh = data.vmesh;
    //~ const bounded_matrix<double,nnodes,dim>& vconv = v - vmesh;
    //~ const array_1d<double,nnodes>& eta = data.eta;          // Unused variable
    //~ const array_1d<double,nnodes>& etan = data.etan;        // Unused variable
    //~ const array_1d<double,nnodes>& etann = data.etann;      // Unused variable

    // Get constitutive matrix
    //~ const Matrix& C = data.C;                               // Unused variable

    // const array_1d<double,dim> vconv_gauss = prod(trans(vconv), N);

    // const double vconv_norm = norm_2(vconv_gauss);

    const double clhs0 =             H*g;
const double clhs1 =             N[0]*bdf0;
const double clhs2 =             DN(0,0)*H*g;
const double clhs3 =             DN(0,1)*H*g;
const double clhs4 =             -DN(1,0)*clhs2 - DN(1,1)*clhs3 - N[1]*clhs1;
const double clhs5 =             -DN(2,0)*clhs2 - DN(2,1)*clhs3 - N[2]*clhs1;
const double clhs6 =             -DN(1,0)*DN(2,0)*clhs0 - DN(1,1)*DN(2,1)*clhs0 - N[1]*N[2]*bdf0;
            lhs(0,0)=-pow(DN(0,0), 2)*clhs0 - pow(DN(0,1), 2)*clhs0 - pow(N[0], 2)*bdf0;
            lhs(0,1)=clhs4;
            lhs(0,2)=clhs5;
            lhs(1,0)=clhs4;
            lhs(1,1)=-pow(DN(1,0), 2)*clhs0 - pow(DN(1,1), 2)*clhs0 - pow(N[1], 2)*bdf0;
            lhs(1,2)=clhs6;
            lhs(2,0)=clhs5;
            lhs(2,1)=clhs6;
            lhs(2,2)=-pow(DN(2,0), 2)*clhs0 - pow(DN(2,1), 2)*clhs0 - pow(N[2], 2)*bdf0;


}


template<>
void Helmholtz<2>::ComputeGaussPointRHSContribution(array_1d<double,3>& rhs, const ElementDataStruct& data)
{
    const int nnodes = 3;
    const int dim = 2;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;
    
    // Get problem data
    const double H = inner_prod(N, data.H);                // Bathymetry
    const double g = data.g;                               // Gravity
    //~ const double l = data.l;                               // Characteristic element size    // Unused variable

    //~ const double& dt = data.dt;             // Unused variable
    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;

    //~ const bounded_matrix<double,nnodes,dim>& v = data.v;
    //~ const bounded_matrix<double,nnodes,dim>& vn = data.vn;
    //~ const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
    //~ const bounded_matrix<double,nnodes,dim>& vmesh = data.vmesh;
    //~ const bounded_matrix<double,nnodes,dim>& vconv = v - vmesh;
    const array_1d<double,nnodes>& eta = data.eta;
    const array_1d<double,nnodes>& etan = data.etan;
    const array_1d<double,nnodes>& etann = data.etann;

    // Get constitutive matrix
    // const Matrix& C = data.C;

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,dim> grad_eta = prod(trans(DN), eta);
    // const array_1d<double,dim> vconv_gauss = prod(trans(vconv), N);
    //~ const double eta_gauss = inner_prod(N,eta);             // Unused variable

    // const double vconv_norm = norm_2(vconv_gauss);

    //~ array_1d<double,dim> accel_gauss = bdf0*v_gauss;
    //~ noalias(accel_gauss) += bdf1*prod(trans(vn), N);
    //~ noalias(accel_gauss) += bdf2*prod(trans(vnn), N);

    const double crhs0 =             H*g*(DN(0,0)*eta[0] + DN(1,0)*eta[1] + DN(2,0)*eta[2]);
const double crhs1 =             H*g*(DN(0,1)*eta[0] + DN(1,1)*eta[1] + DN(2,1)*eta[2]);
const double crhs2 =             N[0]*(bdf0*eta[0] + bdf1*etan[0] + bdf2*etann[0]) + N[1]*(bdf0*eta[1] + bdf1*etan[1] + bdf2*etann[1]) + N[2]*(bdf0*eta[2] + bdf1*etan[2] + bdf2*etann[2]);
            rhs[0]=DN(0,0)*crhs0 + DN(0,1)*crhs1 + N[0]*crhs2;
            rhs[1]=DN(1,0)*crhs0 + DN(1,1)*crhs1 + N[1]*crhs2;
            rhs[2]=DN(2,0)*crhs0 + DN(2,1)*crhs1 + N[2]*crhs2;

}


template<>
double Helmholtz<2>::SubscaleErrorEstimate(const ElementDataStruct& data)  // Not finished
{
    const int nnodes = 3;
    const int dim = 2;

    // Get shape function values
    //~ const array_1d<double,nnodes>& N = data.N;                          // Unused variable
    const bounded_matrix<double,nnodes,dim>& DN = data.DN_DX;

    // Get problem data
    //~ const double H = inner_prod(data.N, data.H);           // Bathymetry                        // Unused variable
    //~ const double g = data.g;                               // Gravity                           // Unused variable
    //~ const double l = data.l;                               // Characteristic element size       // Unused variable

    //~ const double& dt = data.dt;             // Unused variable
    //~ const double& bdf0 = data.bdf0;         // Unused variable
    //~ const double& bdf1 = data.bdf1;         // Unused variable
    //~ const double& bdf2 = data.bdf2;         // Unused variable

    //~ const bounded_matrix<double,nnodes,dim>& v = data.v;
    //~ const bounded_matrix<double,nnodes,dim>& vn = data.vn;
    //~ const bounded_matrix<double,nnodes,dim>& vnn = data.vnn;
    //~ const bounded_matrix<double,nnodes,dim>& vmesh = data.vmesh;
    //~ const bounded_matrix<double,nnodes,dim>& vconv = v - vmesh;
    const array_1d<double,nnodes>& eta = data.eta;
    //~ const array_1d<double,nnodes>& etan = data.etan;                   // Unused variable
    //~ const array_1d<double,nnodes>& etann = data.etann;                 // Unused variable

    // Auxiliary variables used in the calculation of the error estimator
    // array_1d<double,dim> v_s_gauss;
    // const array_1d<double,dim> v_gauss = prod(trans(v), N);
    const array_1d<double,dim> grad_eta = prod(trans(DN), eta);
    // const array_1d<double,dim> vconv_gauss = prod(trans(vconv), N);

    // const double vconv_norm = norm_2(vconv_gauss);

    // Gauss point velocity subscale value computation
    //substitute_gausspt_subscale_2D

    // const double v_gauss_norm = norm_2(v_gauss);
    // const double v_s_gauss_norm = norm_2(v_s_gauss);

    return 1;
}

}
