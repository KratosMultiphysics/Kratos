//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla, Mengjie Zhao
//

#include "custom_elements/time_averaged_navier_stokes.h"

namespace Kratos {

template<>
void TimeAveragedNavierStokes<3>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int dim = 3;
    constexpr unsigned int num_nodes = 4;
    constexpr unsigned int dof_size = num_nodes*(dim+1);

    if (rResult.size() != dof_size){
        rResult.resize(dof_size, false);
    }

    for(unsigned int i=0; i<num_nodes; i++){
        rResult[i*(dim+1)  ]  =  this->GetGeometry()[i].GetDof(TIME_AVERAGED_VELOCITY_X).EquationId();
        rResult[i*(dim+1)+1]  =  this->GetGeometry()[i].GetDof(TIME_AVERAGED_VELOCITY_Y).EquationId();
        rResult[i*(dim+1)+2]  =  this->GetGeometry()[i].GetDof(TIME_AVERAGED_VELOCITY_Z).EquationId();
        rResult[i*(dim+1)+3]  =  this->GetGeometry()[i].GetDof(TIME_AVERAGED_PRESSURE).EquationId();
    }

    KRATOS_CATCH("")
}


template<>
void TimeAveragedNavierStokes<2>::EquationIdVector(
    EquationIdVectorType& rResult,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int dim = 2;
    constexpr unsigned int num_nodes = 3;
    constexpr unsigned int dof_size  = num_nodes*(dim+1);

    if (rResult.size() != dof_size){
        rResult.resize(dof_size, false);
    }

    for(unsigned int i=0; i<num_nodes; i++){
        rResult[i*(dim+1)  ]  =  this->GetGeometry()[i].GetDof(TIME_AVERAGED_VELOCITY_X).EquationId();
        rResult[i*(dim+1)+1]  =  this->GetGeometry()[i].GetDof(TIME_AVERAGED_VELOCITY_Y).EquationId();
        rResult[i*(dim+1)+2]  =  this->GetGeometry()[i].GetDof(TIME_AVERAGED_PRESSURE).EquationId();
    }

    KRATOS_CATCH("")
}


template<>
void TimeAveragedNavierStokes<3>::GetDofList(
    DofsVectorType& ElementalDofList,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int dim = 3;
    constexpr unsigned int num_nodes = 4;
    constexpr unsigned int dof_size  = num_nodes*(dim+1);

    if (ElementalDofList.size() != dof_size){
        ElementalDofList.resize(dof_size);
    }

    for(unsigned int i=0; i<num_nodes; i++){
        ElementalDofList[i*(dim+1)  ]  =  this->GetGeometry()[i].pGetDof(TIME_AVERAGED_VELOCITY_X);
        ElementalDofList[i*(dim+1)+1]  =  this->GetGeometry()[i].pGetDof(TIME_AVERAGED_VELOCITY_Y);
        ElementalDofList[i*(dim+1)+2]  =  this->GetGeometry()[i].pGetDof(TIME_AVERAGED_VELOCITY_Z);
        ElementalDofList[i*(dim+1)+3]  =  this->GetGeometry()[i].pGetDof(TIME_AVERAGED_PRESSURE);
    }

    KRATOS_CATCH("");
}


template<>
void TimeAveragedNavierStokes<2>::GetDofList(
    DofsVectorType& ElementalDofList,
    ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    constexpr unsigned int dim = 2;
    constexpr unsigned int num_nodes = 3;
    constexpr unsigned int dof_size  = num_nodes*(dim+1);

    if (ElementalDofList.size() != dof_size){
        ElementalDofList.resize(dof_size);
    }

    for(unsigned int i=0; i<num_nodes; i++){
        ElementalDofList[i*(dim+1)  ]  =  this->GetGeometry()[i].pGetDof(TIME_AVERAGED_VELOCITY_X);
        ElementalDofList[i*(dim+1)+1]  =  this->GetGeometry()[i].pGetDof(TIME_AVERAGED_VELOCITY_Y);
        ElementalDofList[i*(dim+1)+2]  =  this->GetGeometry()[i].pGetDof(TIME_AVERAGED_PRESSURE);
    }

    KRATOS_CATCH("");
}


template<>
void TimeAveragedNavierStokes<3>::ComputeGaussPointLHSContribution(
    BoundedMatrix<double,16,16>& lhs,
    const ElementDataStruct& data)
{
    constexpr int dim = 3;
    constexpr int nnodes = 4;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size
    const double c = data.c;                                // Wave velocity

    const double& dts = data.dts;                           // The averaging time period

    // time step history info
    const double& dt = data.dt;
    const double& dtn = data.dtn;
    const double& dtnn = data.dtnn;
    const double& dtnnn = data.dtnnn;

    // time history info
    const double& t = data.t;
    const double& tn = data.tn;
    const double& tnn = data.tnn;
    const double& tnnn = data.tnnn;

    const double& bdf0 = data.bdf0;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v_ave = data.v_ave;
    const BoundedMatrix<double,nnodes,dim>& vn_ave = data.vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;

    // time averaging parameters
    double ave_c1 = data.ave_c1;
    double ave_c2 = data.ave_c2;
    double ave_n_c1 = data.ave_n_c1;
    double ave_n_c2 = data.ave_n_c2;
    double ave_nn_c1 = data.ave_nn_c1;
    double ave_nn_c2 = data.ave_nn_c2;

    // get time accurate expression
    const BoundedMatrix<double,nnodes,dim>& v = ave_c1 * v_ave - ave_c2 * vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vconv = v - vmesh;

    // Get constitutive matrix
    const Matrix& C = data.C;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    //substitute_lhs_3D

}


template<>
void TimeAveragedNavierStokes<2>::ComputeGaussPointLHSContribution(
    BoundedMatrix<double,9,9>& lhs,
    const ElementDataStruct& data)
{
    constexpr int dim = 2;
    constexpr int nnodes = 3;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size
    const double c = data.c;                                // Wave velocity

    const double& dts = data.dts;                           // The averaging time period

    // time step history info
    const double& dt = data.dt;
    const double& dtn = data.dtn;
    const double& dtnn = data.dtnn;
    const double& dtnnn = data.dtnnn;

    // time history info
    const double& t = data.t;
    const double& tn = data.tn;
    const double& tnn = data.tnn;
    const double& tnnn = data.tnnn;

    const double& bdf0 = data.bdf0;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v_ave = data.v_ave;
    const BoundedMatrix<double,nnodes,dim>& vn_ave = data.vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;

    // time averaging parameters
    double ave_c1 = data.ave_c1;
    double ave_c2 = data.ave_c2;
    double ave_n_c1 = data.ave_n_c1;
    double ave_n_c2 = data.ave_n_c2;
    double ave_nn_c1 = data.ave_nn_c1;
    double ave_nn_c2 = data.ave_nn_c2;

    // get time accurate expression
    const BoundedMatrix<double,nnodes,dim>& v = ave_c1 * v_ave - ave_c2 * vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vconv = v - vmesh;

    // get constitutive matrix
    const Matrix& C = data.C;

    // get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    //substitute_lhs_2D

}


template<>
void TimeAveragedNavierStokes<3>::ComputeGaussPointRHSContribution(
    array_1d<double,16>& rhs,
    const ElementDataStruct& data)
{
    constexpr int dim = 3;
    constexpr int nnodes = 4;
    constexpr int strain_size = 6;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size
    const double c = data.c;                                // Wave velocity

    const double& dts = data.dts;                           // The averaging time period

    // time step history info
    const double& dt = data.dt;
    const double& dtn = data.dtn;
    const double& dtnn = data.dtnn;
    const double& dtnnn = data.dtnnn;

    // time history info
    const double& t = data.t;
    const double& tn = data.tn;
    const double& tnn = data.tnn;
    const double& tnnn = data.tnnn;

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v_ave = data.v_ave;
    const BoundedMatrix<double,nnodes,dim>& vn_ave = data.vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnn_ave = data.vnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnnn_ave = data.vnnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& p_ave = data.p_ave;
    const array_1d<double,nnodes>& pn_ave = data.pn_ave;
    const array_1d<double,nnodes>& pnn_ave = data.pnn_ave;
    const array_1d<double,nnodes>& pnnn_ave = data.pnnn_ave;
    const array_1d<double,strain_size>& stress = data.stress;

    // time averaging parameters
    double ave_c1 = data.ave_c1;
    double ave_c2 = data.ave_c2;
    double ave_n_c1 = data.ave_n_c1;
    double ave_n_c2 = data.ave_n_c2;
    double ave_nn_c1 = data.ave_nn_c1;
    double ave_nn_c2 = data.ave_nn_c2;

    // get time accurate expression
    const array_1d<double,nnodes>& p = ave_c1 * p_ave - ave_c2 * pn_ave;
    const BoundedMatrix<double,nnodes,dim>& v = ave_c1 * v_ave - ave_c2 * vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vconv = v - vmesh;

    // get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // auxiliary variables used in the calculation of the RHS
    const array_1d<double,dim> f_gauss = prod(trans(f), N);
    const array_1d<double,dim> grad_p = prod(trans(DN), p);

    // stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    //substitute_rhs_3D

}


template<>
void TimeAveragedNavierStokes<2>::ComputeGaussPointRHSContribution(
    array_1d<double,9>& rhs,
    const ElementDataStruct& data)
{
    constexpr int dim = 2;
    constexpr int nnodes = 3;
    constexpr int strain_size = 3;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size
    const double c = data.c;                                // Wave velocity

    const double& dts = data.dts;                           // The averaging time period

    // time step history info
    const double& dt = data.dt;
    const double& dtn = data.dtn;
    const double& dtnn = data.dtnn;
    const double& dtnnn = data.dtnnn;

    // time history info
    const double& t = data.t;
    const double& tn = data.tn;
    const double& tnn = data.tnn;
    const double& tnnn = data.tnnn;

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v_ave = data.v_ave;
    const BoundedMatrix<double,nnodes,dim>& vn_ave = data.vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnn_ave = data.vnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnnn_ave = data.vnnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& p_ave = data.p_ave;
    const array_1d<double,nnodes>& pn_ave = data.pn_ave;
    const array_1d<double,nnodes>& pnn_ave = data.pnn_ave;
    const array_1d<double,nnodes>& pnnn_ave = data.pnnn_ave;
    const array_1d<double,strain_size>& stress = data.stress;

    // time averaging parameters
    double ave_c1 = data.ave_c1;
    double ave_c2 = data.ave_c2;
    double ave_n_c1 = data.ave_n_c1;
    double ave_n_c2 = data.ave_n_c2;
    double ave_nn_c1 = data.ave_nn_c1;
    double ave_nn_c2 = data.ave_nn_c2;

    // get time accurate expression
    const array_1d<double,nnodes>& p = ave_c1 * p_ave - ave_c2 * pn_ave;
    const BoundedMatrix<double,nnodes,dim>& v = ave_c1 * v_ave - ave_c2 * vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vconv = v - vmesh;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the RHS
    const array_1d<double,dim> f_gauss = prod(trans(f), N);
    const array_1d<double,dim> grad_p = prod(trans(DN), p);

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    //substitute_rhs_2D

}


template<>
double TimeAveragedNavierStokes<3>::SubscaleErrorEstimate(const ElementDataStruct& data)
{
    constexpr int dim = 3;
    constexpr int nnodes = 4;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size
    // const double c = data.c;                                // Wave velocity

    const double& dts = data.dts;                           // The averaging time period

    // time step history info
    const double& dt = data.dt;
    const double& dtn = data.dtn;
    const double& dtnn = data.dtnn;
    const double& dtnnn = data.dtnnn;

    // time history info
    const double& t = data.t;
    const double& tn = data.tn;
    const double& tnn = data.tnn;
    const double& tnnn = data.tnnn;

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v_ave = data.v_ave;
    const BoundedMatrix<double,nnodes,dim>& vn_ave = data.vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnn_ave = data.vnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnnn_ave = data.vnnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& p_ave = data.p_ave;
    const array_1d<double,nnodes>& pn_ave = data.pn_ave;

    // time averaging parameters
    double ave_c1 = data.ave_c1;
    double ave_c2 = data.ave_c2;
    double ave_n_c1 = data.ave_n_c1;
    double ave_n_c2 = data.ave_n_c2;
    double ave_nn_c1 = data.ave_nn_c1;
    double ave_nn_c2 = data.ave_nn_c2;

    // get time accurate expression
    const array_1d<double,nnodes>& p = ave_c1 * p_ave - ave_c2 * pn_ave;
    const BoundedMatrix<double,nnodes,dim>& v = ave_c1 * v_ave - ave_c2 * vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vconv = v - vmesh;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the error estimator
    array_1d<double,dim> v_s_gauss;
    const array_1d<double,dim> v_gauss = prod(trans(v), N);
    const array_1d<double,dim> f_gauss = prod(trans(f), N);
    const array_1d<double,dim> grad_p = prod(trans(DN), p);

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Gauss point velocity subscale value computation
    //substitute_gausspt_subscale_3D

    const double v_gauss_norm = norm_2(v_gauss);
    const double v_s_gauss_norm = norm_2(v_s_gauss);

    return v_s_gauss_norm/v_gauss_norm;
}


template<>
double TimeAveragedNavierStokes<2>::SubscaleErrorEstimate(const ElementDataStruct& data)
{
    constexpr int dim = 2;
    constexpr int nnodes = 3;

    const double rho = inner_prod(data.N, data.rho);        // Density
    const double mu = inner_prod(data.N, data.mu);          // Dynamic viscosity
    const double h = data.h;                                // Characteristic element size

    const double& dts = data.dts;                           // The averaging time period

    // time step history info
    const double& dt = data.dt;
    const double& dtn = data.dtn;
    const double& dtnn = data.dtnn;
    const double& dtnnn = data.dtnnn;

    // time history info
    const double& t = data.t;
    const double& tn = data.tn;
    const double& tnn = data.tnn;
    const double& tnnn = data.tnnn;

    const double& bdf0 = data.bdf0;
    const double& bdf1 = data.bdf1;
    const double& bdf2 = data.bdf2;
    const double& dyn_tau = data.dyn_tau;

    const BoundedMatrix<double,nnodes,dim>& v_ave    = data.v_ave;
    const BoundedMatrix<double,nnodes,dim>& vn_ave   = data.vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnn_ave  = data.vnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vnnn_ave = data.vnnn_ave;
    const BoundedMatrix<double,nnodes,dim>& vmesh = data.vmesh;
    const BoundedMatrix<double,nnodes,dim>& f = data.f;
    const array_1d<double,nnodes>& p_ave = data.p_ave;
    const array_1d<double,nnodes>& pn_ave = data.pn_ave;

    // time averaging parameters
    double ave_c1 = data.ave_c1;
    double ave_c2 = data.ave_c2;
    double ave_n_c1 = data.ave_n_c1;
    double ave_n_c2 = data.ave_n_c2;
    double ave_nn_c1 = data.ave_nn_c1;
    double ave_nn_c2 = data.ave_nn_c2;

    // get time accurate expression
    const array_1d<double,nnodes>& p = ave_c1 * p_ave - ave_c2 * pn_ave;
    const BoundedMatrix<double,nnodes,dim>& v = ave_c1 * v_ave - ave_c2 * vn_ave;
    const BoundedMatrix<double,nnodes,dim>& vconv = v - vmesh;

    // Get shape function values
    const array_1d<double,nnodes>& N = data.N;
    const BoundedMatrix<double,nnodes,dim>& DN = data.DN_DX;

    // Auxiliary variables used in the calculation of the error estimator
    array_1d<double,dim> v_s_gauss;
    const array_1d<double,dim> v_gauss = prod(trans(v), N);
    const array_1d<double,dim> f_gauss = prod(trans(f), N);
    const array_1d<double,dim> grad_p = prod(trans(DN), p);

    // Stabilization parameters
    constexpr double stab_c1 = 4.0;
    constexpr double stab_c2 = 2.0;

    // Gauss point velocity subscale value computation
    //substitute_gausspt_subscale_2D

    const double v_gauss_norm = norm_2(v_gauss);
    const double v_s_gauss_norm = norm_2(v_s_gauss);

    return v_s_gauss_norm/v_gauss_norm;
}

}
