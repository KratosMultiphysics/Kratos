/*
==============================================================================
KratosR1IncompressibleFluidApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: pryzhakov $
//   Date:                $Date: 2008-11-26 08:17:41 $
//   Revision:            $Revision: 1.0 $
//
//

//#define GRADPN_FORM

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/fluid_3dGLS_expl_comp.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{


//THIS IS A COMPRESSIBLE FLUID ELEMENT, WITH GLS STABILIZATION, RUNGE-KUTTA Momentum Time integration, FRACTIONAL STEP
//************************************************************************************
//************************************************************************************
Fluid3DGLS_expl_comp::Fluid3DGLS_expl_comp(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
Fluid3DGLS_expl_comp::Fluid3DGLS_expl_comp(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{


}
void Fluid3DGLS_expl_comp::CalculateLumpedMass()
{
    //note that for the compressible case, rho will be also a variable

    const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
    const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
    const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
    const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);

    //double Volume;
    //GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Volume);
    double Volume = GeometryUtils::CalculateVolume3D(GetGeometry());
    double lumped_mass_fac = Volume * 0.25;
    //filling in the diagonal of the lumped mass matrix,  (later I can change it to vector...)
    GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS)+=lumped_mass_fac*rho0;
    GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS)+=lumped_mass_fac*rho1;
    GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS)+=lumped_mass_fac*rho2;
    GetGeometry()[3].FastGetSolutionStepValue(NODAL_MASS)+=lumped_mass_fac*rho3;
}

Element::Pointer Fluid3DGLS_expl_comp::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY

    return Element::Pointer(new Fluid3DGLS_expl_comp(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
}

Fluid3DGLS_expl_comp::~Fluid3DGLS_expl_comp()
{
}

//************************************************************************************
//************************************************************************************


//************************************************************************************
//************************************************************************************
void Fluid3DGLS_expl_comp::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
}

//************************************************************************************
//************************************************************************************

void Fluid3DGLS_expl_comp::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
{
    //KRATOS_WATCH("Empty function for this element")
    //KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
}

void Fluid3DGLS_expl_comp::CalculateGalerkinMomentumResidual(VectorType& GalerkinRHS)
{
    KRATOS_TRY

    boost::numeric::ublas::bounded_matrix<double,12,12> msAuxMat;
    boost::numeric::ublas::bounded_matrix<double,12,12> msAuxMat1;
    boost::numeric::ublas::bounded_matrix<double,12,3> msAuxMat2;
    boost::numeric::ublas::bounded_matrix<double,3,3> msGrad_ug;
    array_1d<double,12> msAuxVec;
    array_1d<double,12> msStabMomRes;
    boost::numeric::ublas::bounded_matrix<double,4,4> msWorkMatrix;
    boost::numeric::ublas::bounded_matrix<double,12,3> msShapeFunc;
    boost::numeric::ublas::bounded_matrix<double,3,12> msConvOp;
    boost::numeric::ublas::bounded_matrix<double,12,4> msGradOp;
    boost::numeric::ublas::bounded_matrix<double,4,3> msDN_DX;
    array_1d<double,3> ms_adv_vel;
    array_1d<double,4> msN;
    array_1d<double,3> ms_vel_gauss;
    array_1d<double,4> ms_temp_vec_np;
    array_1d<double,4> ms_aux0;
    array_1d<double,4> ms_aux1;
    array_1d<double,3> ms_aux2;

    //first we compute  the force term and pressure gradient terms:
    //getting data for the given geometry
    double Volume;
    GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Volume);

    //getting the velocity on the nodes and other necessary variabless
    const array_1d<double,3> vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    double p_n0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
    const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
    const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);

    const array_1d<double,3> vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
    double p_n1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
    const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
    const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);

    const array_1d<double,3>& vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);
    double p_n2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
    const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
    const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);

    const array_1d<double,3>& vel3 = GetGeometry()[3].FastGetSolutionStepValue(VELOCITY);
    double p_n3 = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE);
    const double nu3 = GetGeometry()[3].FastGetSolutionStepValue(VISCOSITY);
    const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);

    //====================================================================
    //calculating viscosity and density
    double nu = 0.25*(nu0 + nu1 + nu2 +nu3);
    double density = 0.25*(rho0 + rho1 + rho2 + rho3);

    //VISCOUS CONTRIBUTION
    // += Laplacian * nu; --> ONE GAUSS POINT
    //msWorkMatrix is used now to store the element laplacian 4x4
    noalias(msWorkMatrix) = Volume*density*nu * prod(msDN_DX,trans(msDN_DX));

    //x comp
    GalerkinRHS[0]=-1.0*(msWorkMatrix(0,0)*vel0[0]+msWorkMatrix(0,1)*vel1[0]+msWorkMatrix(0,2)*vel2[0] + msWorkMatrix(0,3)*vel3[0]);
    //y comp
    GalerkinRHS[1]=-1.0*(msWorkMatrix(0,0)*vel0[1]+msWorkMatrix(0,1)*vel1[1]+msWorkMatrix(0,2)*vel2[1] + msWorkMatrix(0,3)*vel3[1]);
    //z comp
    GalerkinRHS[2]=-1.0*(msWorkMatrix(0,0)*vel0[2]+msWorkMatrix(0,1)*vel1[2]+msWorkMatrix(0,2)*vel2[2] + msWorkMatrix(0,3)*vel3[2]);

    //x comp
    GalerkinRHS[3]=-1.0*(msWorkMatrix(1,0)*vel0[0]+msWorkMatrix(1,1)*vel1[0]+msWorkMatrix(1,2)*vel2[0] + msWorkMatrix(1,3)*vel3[0]);
    //y comp
    GalerkinRHS[4]=-1.0*(msWorkMatrix(1,0)*vel0[1]+msWorkMatrix(1,1)*vel1[1]+msWorkMatrix(1,2)*vel2[1] + msWorkMatrix(1,3)*vel3[1]);
    //z comp
    GalerkinRHS[5]=-1.0*(msWorkMatrix(1,0)*vel0[2]+msWorkMatrix(1,1)*vel1[2]+msWorkMatrix(1,2)*vel2[2] + msWorkMatrix(1,3)*vel3[2]);

    //x comp
    GalerkinRHS[6]=-1.0*(msWorkMatrix(2,0)*vel0[0]+msWorkMatrix(2,1)*vel1[0]+msWorkMatrix(2,2)*vel2[0] + msWorkMatrix(2,3)*vel3[0]);
    //y comp
    GalerkinRHS[7]=-1.0*(msWorkMatrix(2,0)*vel0[1]+msWorkMatrix(2,1)*vel1[1]+msWorkMatrix(2,2)*vel2[1] + msWorkMatrix(2,3)*vel3[1]);
    //z comp
    GalerkinRHS[8]=-1.0*(msWorkMatrix(2,0)*vel0[2]+msWorkMatrix(2,1)*vel1[2]+msWorkMatrix(2,2)*vel2[2] + msWorkMatrix(2,3)*vel3[2]);

    //x comp
    GalerkinRHS[9]=-1.0*(msWorkMatrix(3,0)*vel0[0]+msWorkMatrix(3,1)*vel1[0]+msWorkMatrix(3,2)*vel2[0] + msWorkMatrix(3,3)*vel3[0]);
    //y comp
    GalerkinRHS[10]=-1.0*(msWorkMatrix(3,0)*vel0[1]+msWorkMatrix(3,1)*vel1[1]+msWorkMatrix(3,2)*vel2[1] + msWorkMatrix(3,3)*vel3[1]);
    //z comp
    GalerkinRHS[11]=-1.0*(msWorkMatrix(3,0)*vel0[2]+msWorkMatrix(3,1)*vel1[2]+msWorkMatrix(3,2)*vel2[2] + msWorkMatrix(3,3)*vel3[2]);



    //convective contribution. Note that N[0]=N[1]=N[2]=0.33333 for our 1-Gauss Point integration
    //WATCH OUT THAT I AM NOT PLUGGING THE NEW ADVECTIVE VELOCITY WHEN EXECUTING RUNGE KUTTA - I use u_n (and not the intermediate u_aux)
    //KRATOS_WATCH("Now lets see N - they should be all equal 0.3333")
    //KRATOS_WATCH(msN)
    //KRATOS_WATCH(msDN_DX)
    ms_adv_vel[0] = msN[0]*(vel0[0])+msN[1]*(vel1[0])+msN[2]*(vel2[0])+msN[3]*(vel3[0]);
    ms_adv_vel[1] = msN[0]*(vel0[1])+msN[1]*(vel1[1])+msN[2]*(vel2[1])+msN[3]*(vel3[1]);
    ms_adv_vel[2] = msN[0]*(vel0[2])+msN[1]*(vel1[2])+msN[2]*(vel2[2])+msN[3]*(vel3[2]);

    //calculate convective term
    int nodes_number = 4;
    int dof = 3;


    for (int ii = 0; ii< nodes_number; ii++)
    {
        int column = ii*dof;
        msConvOp(0,column) = msDN_DX(ii,0)*ms_adv_vel[0] + msDN_DX(ii,1)*ms_adv_vel[1] + msDN_DX(ii, 2) * ms_adv_vel[2];
        msConvOp(1,column + 1) = msConvOp(0,column);
        msConvOp(2,column + 2) = msConvOp(0,column);

        msShapeFunc(column,0) = msN[ii];
        msShapeFunc(column + 1, 1) = msShapeFunc(column,0);
        msShapeFunc(column + 2, 2) = msShapeFunc(column,0);
    }
    //and now we store the convective term to the msAuxMat
    noalias(msAuxMat) = prod(msShapeFunc, msConvOp);

    //here we multiply the msAuxMat (convective term) with the KNOWN primary variable velocity (or Momentum in case we switch to the rho*u formulation)at time n and get RHS

    //in this vector msAuxVec now we store the velocities
    msAuxVec[0]=vel0[0];
    msAuxVec[1]=vel0[1];
    msAuxVec[2]=vel0[2];

    msAuxVec[3]=vel1[0];
    msAuxVec[4]=vel1[1];
    msAuxVec[5]=vel1[2];

    msAuxVec[6]=vel2[0];
    msAuxVec[7]=vel2[1];
    msAuxVec[8]=vel2[2];

    msAuxVec[9]=vel3[0];
    msAuxVec[10]=vel3[1];
    msAuxVec[11]=vel3[2];

    noalias(GalerkinRHS)-=(Volume * density)*prod(msAuxMat, msAuxVec);//0.3333333333333*prod(AUX, u_n);
    //GalerkinRHS=prod(AUX, Mom_n)

    //and now we add the pressure gradient and the force term
    //external forces (component)
    const array_1d<double,3> body_force = 0.25*(GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE)+
                                          GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE) +
                                          GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE)+
                                          GetGeometry()[3].FastGetSolutionStepValue(BODY_FORCE));
    unsigned int number_of_nodes=4;
    for(unsigned int i = 0; i<number_of_nodes; i++)
    {
        //f=A*N_I*b, N_I=0.33333333 for 1 Gauss point
        GalerkinRHS[i*3] += body_force[0]* density * Volume * 0.25;
        GalerkinRHS[i*3+1] += body_force[1] * density * Volume * 0.25;
        GalerkinRHS[i*3+2] += body_force[2] * density * Volume * 0.25;
    }

    //pressure gradient (sign is positive), ITS NOT INTEGRATED BY PARTS!!!!0.33333 stands for N at the integration point
    //double a = 0.33333333*Volume*(msDN_DX(0,0)*p_n0 + msDN_DX(1,0)*p_n1 + msDN_DX(2,0)*p_n2);
    //double b = 0.33333333*Volume*(msDN_DX(0,1)*p_n0 + msDN_DX(1,1)*p_n1 + msDN_DX(2,1)*p_n2);

    //Now we shall add the Gp term

    double p_avg=0.25*(p_n0+p_n1+p_n2+p_n3)*Volume;
    GalerkinRHS[0]+=msDN_DX(0,0)*p_avg;
    GalerkinRHS[1]+=msDN_DX(0,1)*p_avg;
    GalerkinRHS[2]+=msDN_DX(0,2)*p_avg;

    GalerkinRHS[3]+=msDN_DX(1,0)*p_avg;
    GalerkinRHS[4]+=msDN_DX(1,1)*p_avg;
    GalerkinRHS[5]+=msDN_DX(1,2)*p_avg;

    GalerkinRHS[6]+=msDN_DX(2,0)*p_avg;
    GalerkinRHS[7]+=msDN_DX(2,1)*p_avg;
    GalerkinRHS[8]+=msDN_DX(2,2)*p_avg;

    GalerkinRHS[9]+=msDN_DX(3,0)*p_avg;
    GalerkinRHS[10]+=msDN_DX(3,1)*p_avg;
    GalerkinRHS[11]+=msDN_DX(3,2)*p_avg;

    /*
    //the below is Grad form (not integrated by parts)
    //G=DN_DX*N
    noalias(msGradOp)=prod(msShapeFunc, trans(msDN_DX));
    msGradOp*=Volume;

    //we use now ms_aux0 for storing the 3 nodal pressure in a vector
    ms_aux0[0]=p_n0;
    ms_aux0[1]=p_n1;
    ms_aux0[2]=p_n2;

    //and now we reuse msAuxVec for storing the product of the gradient operraor with the pressure vector - this will give us RHS contr of Gp
    noalias(msAuxVec) = prod(msGradOp,ms_aux0);
    noalias(GalerkinRHS)-=msAuxVec;
    */

    GetGeometry()[0].FastGetSolutionStepValue(FORCE_X)=GalerkinRHS[0];
    GetGeometry()[0].FastGetSolutionStepValue(FORCE_Y)=GalerkinRHS[1];
    GetGeometry()[0].FastGetSolutionStepValue(FORCE_Z)=GalerkinRHS[2];

    GetGeometry()[1].FastGetSolutionStepValue(FORCE_X)=GalerkinRHS[3];
    GetGeometry()[1].FastGetSolutionStepValue(FORCE_Y)=GalerkinRHS[4];
    GetGeometry()[1].FastGetSolutionStepValue(FORCE_Z)=GalerkinRHS[5];

    GetGeometry()[2].FastGetSolutionStepValue(FORCE_X)=GalerkinRHS[6];
    GetGeometry()[2].FastGetSolutionStepValue(FORCE_Y)=GalerkinRHS[7];
    GetGeometry()[2].FastGetSolutionStepValue(FORCE_Z)=GalerkinRHS[8];

    GetGeometry()[3].FastGetSolutionStepValue(FORCE_X)=GalerkinRHS[9];
    GetGeometry()[3].FastGetSolutionStepValue(FORCE_Y)=GalerkinRHS[10];
    GetGeometry()[3].FastGetSolutionStepValue(FORCE_Z)=GalerkinRHS[11];

    KRATOS_CATCH("")

}

void Fluid3DGLS_expl_comp::CalculateRHSVector(VectorType& Galerkin_RHS, double& dt)
{
    KRATOS_TRY

    boost::numeric::ublas::bounded_matrix<double,12,12> msAuxMat;
    boost::numeric::ublas::bounded_matrix<double,12,12> msAuxMat1;
    boost::numeric::ublas::bounded_matrix<double,12,3> msAuxMat2;
    boost::numeric::ublas::bounded_matrix<double,3,3> msGrad_ug;
    array_1d<double,12> msAuxVec;
    array_1d<double,12> msStabMomRes;
    boost::numeric::ublas::bounded_matrix<double,4,4> msWorkMatrix;
    boost::numeric::ublas::bounded_matrix<double,12,3> msShapeFunc;
    boost::numeric::ublas::bounded_matrix<double,3,12> msConvOp;
    boost::numeric::ublas::bounded_matrix<double,12,4> msGradOp;
    boost::numeric::ublas::bounded_matrix<double,4,3> msDN_DX;
    array_1d<double,3> ms_adv_vel;
    array_1d<double,4> msN;
    array_1d<double,3> ms_vel_gauss;
    array_1d<double,4> ms_temp_vec_np;
    array_1d<double,4> ms_aux0;
    array_1d<double,4> ms_aux1;
    array_1d<double,3> ms_aux2;

    //first we compute  the force term and pressure gradient terms:
    //getting data for the given geometry
    double Volume;
    GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Volume);

    //getting the velocity vector on the nodes

    //getting the velocity on the nodes
    const array_1d<double,3>& adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& vel_old0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,1);
    //const array_1d<double,3>& vel_old0_n1 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,2);
    const array_1d<double,3>& ff0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
    const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
    const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
    const double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);

    const array_1d<double,3>& adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& vel_old1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,1);
    //const array_1d<double,3>& vel_old1_n1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,2);
    const array_1d<double,3>& ff1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
    const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
    const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
    const double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);

    const array_1d<double,3>& adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& vel_old2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,1);
    //const array_1d<double,3>& vel_old2_n1 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,2);
    const array_1d<double,3>& ff2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);
    const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
    const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
    const double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);

    const array_1d<double,3>& adv_vel3 = GetGeometry()[3].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& vel_old3 = GetGeometry()[3].FastGetSolutionStepValue(VELOCITY,1);
    //const array_1d<double,3>& vel_old3_n1 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,2);
    const array_1d<double,3>& ff3 = GetGeometry()[3].FastGetSolutionStepValue(BODY_FORCE);
    const double nu3 = GetGeometry()[3].FastGetSolutionStepValue(VISCOSITY);
    const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);
    const double p3 = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE);
    //====================================================================
    //calculating viscosity
    double nu = 0.25*(nu0 + nu1 + nu2 + nu3 );
    double density = 0.25*(rho0 + rho1 + rho2 + rho3);

    //advective velocity at the gauss point
    ms_adv_vel[0] = msN[0]*(adv_vel0[0])+msN[1]*(adv_vel1[0])+msN[2]*(adv_vel2[0])+msN[3]*(adv_vel3[0]);
    ms_adv_vel[1] = msN[0]*(adv_vel0[1])+msN[1]*(adv_vel1[1])+msN[2]*(adv_vel2[1])+msN[3]*(adv_vel3[1]);
    ms_adv_vel[2] = msN[0]*(adv_vel0[2])+msN[1]*(adv_vel1[2])+msN[2]*(adv_vel2[2])+msN[3]*(adv_vel3[2]);

    //now we use msAuxVec to store the nodal velocity vector 6x1
    msAuxVec[0]=adv_vel0[0];
    msAuxVec[1]=adv_vel0[1];
    msAuxVec[2]=adv_vel0[2];

    msAuxVec[3]=adv_vel1[0];
    msAuxVec[4]=adv_vel1[1];
    msAuxVec[5]=adv_vel1[2];

    msAuxVec[6]=adv_vel2[0];
    msAuxVec[7]=adv_vel2[1];
    msAuxVec[8]=adv_vel2[2];

    msAuxVec[9]=adv_vel3[0];
    msAuxVec[10]=adv_vel3[1];
    msAuxVec[11]=adv_vel3[2];

    //ms_temp_vec_np is the nodal pressure vector
    ms_temp_vec_np[0]=p0;
    ms_temp_vec_np[1]=p1;
    ms_temp_vec_np[2]=p2;
    ms_temp_vec_np[3]=p3;




    //KRATOS_WATCH(pressure)
    //KRATOS_WATCH(Visc_and_Conv)
    //convective contribution


    //initialize the RHS with the Galerkin residual
    noalias(msStabMomRes)=Galerkin_RHS;

    //now compute all the stabilization terms: tau(a*nabla v) (-rho*du/dt + mLapl u - nabla p - a nabla u + f)

    //calculating parameter tau (saved internally to each element). to decrease the stabilization, you can itroduce 1/dt in the denominator
    //double h = pow(6.00*Volume,0.3333333);
    double h = pow(12.00*Volume,0.3333333);
    h*=(2.0/sqrt(3.0));
    //we are doing it explicitly, so adv vel is equal to the v_n
    double norm_u = ms_adv_vel[0]*ms_adv_vel[0] + ms_adv_vel[1]*ms_adv_vel[1] + ms_adv_vel[2]*ms_adv_vel[2];
    norm_u = sqrt(norm_u);
    double tau = 1.00 / ( 4.00*nu/(h*h) + 2.00*norm_u/h+1.0/dt);
    //calculate convective term
    int nodes_number = 4;
    int dof = 3;
//		int matsize = dof*nodes_number;

    for (int ii = 0; ii< nodes_number; ii++)
    {
        int column = ii*dof;
        msConvOp(0,column) = msDN_DX(ii,0)*ms_adv_vel[0] + msDN_DX(ii,1)*ms_adv_vel[1] + msDN_DX(ii, 2) * ms_adv_vel[2];
        msConvOp(1,column + 1) = msConvOp(0,column);
        msConvOp(2,column + 2) = msConvOp(0,column);

        msShapeFunc(column,0) = msN[ii];
        msShapeFunc(column + 1, 1) = msShapeFunc(column,0);
        msShapeFunc(column + 2, 2) = msShapeFunc(column,0);
    }

    //stabilization that comes from product of two convections - conv_conv_stab -> (a*nabla v, a*nabla u): 6x6
    noalias(msAuxMat) = prod(trans(msConvOp), msConvOp);
    msAuxMat*= (tau*density*Volume);

    //add the convective term(msAuxMat), multiplied by the nodal velocity vector (stored in msAuxVec) to the RHS
    noalias(msStabMomRes)-= prod(msAuxMat, msAuxVec);

    //now the conv_grad_stab
    noalias(msAuxMat2) = prod(trans(msConvOp), trans(msDN_DX));
    msAuxMat2*= (tau*Volume);

    //add the conv_grad term(msAuxMat2), multiplied by the nodal pressure (stored in ms_temp_vec_np) to the RHS
    noalias(msStabMomRes) -= prod(msAuxMat2, ms_temp_vec_np);

    //stabilization that comes from product of convection and body force
    noalias(msAuxMat1) = prod(trans(msConvOp), trans(msShapeFunc));
    msAuxMat1*= tau*density*Volume;
    //we reuse msAuxVec
    msAuxVec[0]=ff0[0];
    msAuxVec[1]=ff0[1];
    msAuxVec[2]=ff0[2];

    msAuxVec[3]=ff1[0];
    msAuxVec[4]=ff1[1];
    msAuxVec[5]=ff1[2];

    msAuxVec[6]=ff2[0];
    msAuxVec[7]=ff2[1];
    msAuxVec[8]=ff2[2];

    msAuxVec[9]=ff3[0];
    msAuxVec[10]=ff3[1];
    msAuxVec[11]=ff3[2];

    noalias(msStabMomRes)+=prod(msAuxMat1, msAuxVec);
    //and again reuse it now to store accelerations
    //WE SHOULD ADD THE INERTIA CONTRIB in a smart way(coz its intrinscially implicit)

    msAuxVec[0]=adv_vel0[0]-vel_old0[0];
    msAuxVec[1]=adv_vel0[1]-vel_old0[1];
    msAuxVec[2]=adv_vel0[2]-vel_old0[2];

    msAuxVec[3]=adv_vel1[0]-vel_old1[0];
    msAuxVec[4]=adv_vel1[1]-vel_old1[1];
    msAuxVec[5]=adv_vel1[2]-vel_old1[2];

    msAuxVec[6]=adv_vel2[0]-vel_old2[0];
    msAuxVec[7]=adv_vel2[1]-vel_old2[1];
    msAuxVec[8]=adv_vel2[2]-vel_old2[2];

    msAuxVec[9]=adv_vel3[0]-vel_old3[0];
    msAuxVec[10]=adv_vel3[1]-vel_old3[1];
    msAuxVec[11]=adv_vel3[2]-vel_old3[2];

    //always using the a_n in the stabilization
    /*
    msAuxVec[0]=vel_old0[0]-vel_old0_n1[0];
    msAuxVec[1]=vel_old0[1]-vel_old0_n1[1];
    msAuxVec[2]=vel_old1[0]-vel_old1_n1[0];
    msAuxVec[3]=vel_old1[1]-vel_old1_n1[1];
    msAuxVec[4]=vel_old2[0]-vel_old2_n1[0];
    msAuxVec[5]=vel_old2[1]-vel_old2_n1[1];
    */
    msAuxVec*=(1.00/dt);
    //WE SHOULD ADD THE INERTIA CONTRIB in a smart way(coz its intrinscially implicit)
    //we again reuse msAuMat to store cinv_inret_stab, and msAuxVec - to store the accelerations
    //stabilization that comes from product of convection and inertia term
    noalias(msAuxMat) = prod(trans(msConvOp), trans(msShapeFunc));
    msAuxMat*= (tau*density*Volume);

    noalias(msStabMomRes)-=prod(msAuxMat, msAuxVec);

    //here we multiply the AUX (sum of viscous and convective terms) with the Galerkin Residual.. and add to Galerkin residual, to finally get
    //the stabilized RHS
    array_1d<double,3>& rhs0 = GetGeometry()[0].FastGetSolutionStepValue(RHS_VECTOR);
    rhs0[0] += msStabMomRes[0];
    rhs0[1] += msStabMomRes[1];
    rhs0[2] += msStabMomRes[2];

    array_1d<double,3>& rhs1 = GetGeometry()[1].FastGetSolutionStepValue(RHS_VECTOR);
    rhs1[0] += msStabMomRes[3];
    rhs1[1] += msStabMomRes[4];
    rhs1[2] += msStabMomRes[5];

    array_1d<double,3>& rhs2 = GetGeometry()[2].FastGetSolutionStepValue(RHS_VECTOR);
    rhs2[0] += msStabMomRes[6];
    rhs2[1] += msStabMomRes[7];
    rhs2[2] += msStabMomRes[8];

    array_1d<double,3>& rhs3 = GetGeometry()[3].FastGetSolutionStepValue(RHS_VECTOR);
    rhs3[0] += msStabMomRes[9];
    rhs3[1] += msStabMomRes[10];
    rhs3[2] += msStabMomRes[11];

    KRATOS_CATCH("")

}


void Fluid3DGLS_expl_comp::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    boost::numeric::ublas::bounded_matrix<double,12,12> msAuxMat;
    boost::numeric::ublas::bounded_matrix<double,12,12> msAuxMat1;
    boost::numeric::ublas::bounded_matrix<double,12,3> msAuxMat2;
    boost::numeric::ublas::bounded_matrix<double,3,3> msGrad_ug;
    array_1d<double,12> msAuxVec;
    array_1d<double,12> msStabMomRes;
    boost::numeric::ublas::bounded_matrix<double,4,4> msWorkMatrix;
    boost::numeric::ublas::bounded_matrix<double,12,3> msShapeFunc;
    boost::numeric::ublas::bounded_matrix<double,3,12> msConvOp;
    boost::numeric::ublas::bounded_matrix<double,12,4> msGradOp;
    boost::numeric::ublas::bounded_matrix<double,4,3> msDN_DX;
    array_1d<double,3> ms_adv_vel;
    array_1d<double,4> msN;
    array_1d<double,3> ms_vel_gauss;
    array_1d<double,4> ms_temp_vec_np;
    array_1d<double,4> ms_aux0;
    array_1d<double,4> ms_aux1;
    array_1d<double,3> ms_aux2;

    if(rRightHandSideVector.size() != 4)
    {
        rLeftHandSideMatrix.resize(4,4,false);
        rRightHandSideVector.resize(4,false);
    }

    double dt = rCurrentProcessInfo[DELTA_TIME];

    //fract. vel, that is calculated in the first Fractional Step.. but is saved inside the "VELOCITY" VARIABLE
    //so, u_n os VELOCITY, 1 and u_n-1 VELOCITY,2
    const array_1d<double,3>& fv0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& fv0_old = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,1);
    const array_1d<double,3>& fv0_n_1 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,2);
    const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
    const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
    double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
    //double p0_old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT);
    double p0_old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1);
    const array_1d<double,3>& ff0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);

    const array_1d<double,3>& fv1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& fv1_old = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,1);
    const array_1d<double,3>& fv1_n_1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,2);
    const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
    const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
    double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
    //double p1_old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT);
    double p1_old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1);
    const array_1d<double,3>& ff1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);

    const array_1d<double,3>& fv2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& fv2_old = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,1);
    const array_1d<double,3>& fv2_n_1 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,2);
    const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
    const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
    double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
    double p2_old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1);
    //old iteration can be used if we want to iterate between 1st and 2nd fractional steps
    //double p2_old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT);
    const array_1d<double,3>& ff2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);

    const array_1d<double,3>& fv3 = GetGeometry()[3].FastGetSolutionStepValue(VELOCITY);
    const array_1d<double,3>& fv3_old = GetGeometry()[3].FastGetSolutionStepValue(VELOCITY,1);
    const array_1d<double,3>& fv3_n_1 = GetGeometry()[3].FastGetSolutionStepValue(VELOCITY,2);
    const double nu3 = GetGeometry()[3].FastGetSolutionStepValue(VISCOSITY);
    const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);
    double p3 = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE);
    double p3_old = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE,1);
    //old iteration can be used if we want to iterate between 1st and 2nd fractional steps
    //double p2_old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT);
    const array_1d<double,3>& ff3 = GetGeometry()[3].FastGetSolutionStepValue(BODY_FORCE);


    //in msAuxVec we store the velocity, (not the frac step vel, but u_n, the one that enters the stabilization)
    msAuxVec[0]=fv0[0];
    msAuxVec[1]=fv0[1];
    msAuxVec[2]=fv0[2];

    msAuxVec[3]=fv1[0];
    msAuxVec[4]=fv1[1];
    msAuxVec[5]=fv1[2];

    msAuxVec[6]=fv2[0];
    msAuxVec[7]=fv2[1];
    msAuxVec[8]=fv2[2];

    msAuxVec[9]=fv3[0];
    msAuxVec[10]=fv3[1];
    msAuxVec[11]=fv3[2];

    //getting data for the given geometry
    double Volume;
    GeometryUtils::CalculateGeometryData(GetGeometry(),msDN_DX,msN,Volume);

    //calculating average density and viscosity
    double nu = 0.25*(nu0 + nu1 + nu2 + nu3);
    double density = 0.25*(rho0 + rho1 + rho2 + rho3);

    //ms_vel_gauss[i] =  msN[0]*(fv0[i]) + msN[1]*(fv1[i]) +  msN[2]*(fv2[i]);
    //but with one integration N=0.333333333
    ms_vel_gauss[0] =  0.25*(fv0[0]+fv1[0]+fv2[0]+fv3[0]);
    ms_vel_gauss[1] =  0.25*(fv0[1]+fv1[1]+fv2[1]+fv3[1]);
    ms_vel_gauss[2] =  0.25*(fv0[2]+fv1[2]+fv2[2]+fv3[2]);

    //calculating parameter tau (saved internally to each element)
    //h in 3D is calculated like this!
    //double h = pow(6.00*Volume,0.3333333);
    double h = pow(12.00*Volume,0.3333333);
    h*=(2.0/sqrt(3.0));

    double norm_u = ms_vel_gauss[0]*ms_vel_gauss[0] + ms_vel_gauss[1]*ms_vel_gauss[1] + ms_vel_gauss[2]*ms_vel_gauss[2];
    norm_u = sqrt(norm_u);
    double tau = 1.00 / ( 4.00*nu/(h*h) + 2.00*norm_u/h+ 1.0/dt);


    //AND NOW WE ADD THE RESPECTIVE CONTRIBUTIONS TO THE RHS AND LHS of THE SECOND FRAC STEP
    //we use Backward Euler for this step, therefore stab. contribution no RHS +=Tau1*(gradQ, residual)
    //								   and LHS +=Tau1*(gradQ, gradP)
    //laplacian term	       L = Dt * gradN * trans(gradN);
    //stabilization term       Spp = tau * gradN * trans(gradN);
    //WATCH OUT for DIVISION with RHO - check if it changes or not in case of Momentum being the primary Variable
    //
    //	msWorkMatrix stores the element laplacian
    //
    noalias(msWorkMatrix)=prod(msDN_DX,trans(msDN_DX));
    noalias(rLeftHandSideMatrix) = (dt + tau) * Volume*msWorkMatrix;

    //****************************************************************************
    //***************************************************************************
    //
    //		IDEAL GAS CASE (add only in this case)
    //
    //**************************************************************************
    //this is for the compressible case

    Matrix Mass(4,4);
    Mass(0,0)=2.0;
    Mass(0,1)=1.0;
    Mass(0,2)=1.0;
    Mass(0,3)=1.0;
    Mass(1,0)=1.0;
    Mass(1,1)=2.0;
    Mass(1,2)=1.0;
    Mass(1,3)=1.0;
    Mass(2,0)=1.0;
    Mass(2,1)=1.0;
    Mass(2,2)=2.0;
    Mass(2,3)=1.0;
    Mass(3,0)=1.0;
    Mass(3,1)=1.0;
    Mass(3,2)=1.0;
    Mass(3,3)=2.0;

    Mass/=20.0;

    //double p_old_avg=0.33333333*(p0_old+p1_old+p2_old);

    //adiabatic constant gamma
    //double gamma=1.40;
    //double A = 77510.000;
    //const double A=77510.00;
    //gas constant for dry air
    const double R=287.06;
    //*************************************************************************///////////////////////////////////////
    double t0 = GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE);
    double t1 = GetGeometry()[1].FastGetSolutionStepValue(TEMPERATURE);						//
    double t2 = GetGeometry()[2].FastGetSolutionStepValue(TEMPERATURE);
    double t3 = GetGeometry()[3].FastGetSolutionStepValue(TEMPERATURE);
    //
    double t0_old = GetGeometry()[0].FastGetSolutionStepValue(TEMPERATURE,1);
    double t1_old = GetGeometry()[1].FastGetSolutionStepValue(TEMPERATURE,1);					//
    double t2_old = GetGeometry()[2].FastGetSolutionStepValue(TEMPERATURE,1);
    double t3_old = GetGeometry()[3].FastGetSolutionStepValue(TEMPERATURE,1);

    double t_g=0.25*(t0+t1+t2+t3);										//
    double t_n=0.25*(t0_old+t1_old+t2_old+t3_old);
    Matrix TempCheck(4,4);												//
    TempCheck=1.0/(dt*R*t_g)*Volume*Mass;
    noalias(rLeftHandSideMatrix)+=TempCheck;									//
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    ////////////		AND NOW RHS	//////////////////
    //////////////////////////////////////////////////////////

    //Dirichlet contribution  (that is: LHS*p_new)
    ms_temp_vec_np[0] = p0;
    ms_temp_vec_np[1] = p1;
    ms_temp_vec_np[2] = p2;
    ms_temp_vec_np[3] = p3;
    //LHS is already multiplied by AREA
    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix,ms_temp_vec_np);

    //NOW RHS-=dt L p_old
    //changing the meaning of temp_vec_np
    ms_temp_vec_np[0] = p0_old;
    ms_temp_vec_np[1] = p1_old;
    ms_temp_vec_np[2] = p2_old;
    ms_temp_vec_np[3] = p3_old;

    noalias(rRightHandSideVector) += Volume*dt* (prod(msWorkMatrix,ms_temp_vec_np)) ;

    //***************************************************************************
    //***************************************************************************
    //****************************************************************************
    //***************************************************************************
    //
    //		IDEAL GAS CASE (add only in this case)
    //
    //**************************************************************************
    //*************************************************************************
    //p_atm=100000.0
    Vector TempVec(4);
    Vector LumpedM(4);
    LumpedM[0]=0.25;
    LumpedM[1]=0.25;
    LumpedM[2]=0.25;
    LumpedM[4]=0.25;
    //rhs term
    TempVec=(1.0/(dt*R*t_n) )*Volume*(prod(Mass, ms_temp_vec_np));



    //directly using the equation of state (adding 100000 - we have to use absolute pressure)
    //double p_avg=0.33333333333*(p0+p1+p2);
    //TempVec=1.0*(100000.0+p_avg-p_old_avg)/(R*t_n*dt)*Area*LumpedM;
    noalias(rRightHandSideVector)+= TempVec;
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //rhs consists of D*u_tilda (divergence of the Fractional velocity) and the term: Tau1*(nabla_q, residual_gausspoint)
    //fv is u_tilda


    //***************************************************************************

    //here we have the Du_tila term
    double Gaux;
    Gaux =  msDN_DX(0,0)*fv0[0] + msDN_DX(0,1)*fv0[1] + msDN_DX(0,2)*fv0[2];
    Gaux += msDN_DX(1,0)*fv1[0] + msDN_DX(1,1)*fv1[1] + msDN_DX(1,2)*fv1[2];
    Gaux += msDN_DX(2,0)*fv2[0] + msDN_DX(2,1)*fv2[1] + msDN_DX(2,2)*fv2[2];
    Gaux += msDN_DX(3,0)*fv3[0] + msDN_DX(3,1)*fv3[1] + msDN_DX(3,2)*fv3[2];
    rRightHandSideVector[0] -= density*Volume*Gaux * msN[0];
    rRightHandSideVector[1] -= density*Volume*Gaux * msN[1];
    rRightHandSideVector[2] -= density*Volume*Gaux * msN[2];
    rRightHandSideVector[3] -= density*Volume*Gaux * msN[3];



    //RHS+=-Dv
    /*
    boost::numeric::ublas::bounded_matrix<double,3,6> D = ZeroMatrix(3,6);
    boost::numeric::ublas::bounded_matrix<double,6,2> shape_func = ZeroMatrix(6, 2);
    for (int ii = 0; ii< 3; ii++)
        {
    	int column = ii*2;
    	shape_func(column,0) = msN[ii];
    	shape_func(column + 1, 1) = shape_func(column,0);
        }
    noalias(D)=prod(msDN_DX, trans(shape_func));
    temp = prod(D, u_n);
    rRightHandSideVector -= tau*density*Volume*temp;
    */

    //RHS = +tau*nablaN*f, we reuse aux
    //ms_aux0 stores ff_gauss;

    ms_aux2=0.25*(ff0+ff1+ff2+ff3);


    //ms_aux1 - is the product of: (nabla q, f)
    ms_aux1[0]=msDN_DX(0,0)*ms_aux2[0]+msDN_DX(0,1)*ms_aux2[1]+msDN_DX(0,2)*ms_aux2[2];
    ms_aux1[1]=msDN_DX(1,0)*ms_aux2[0]+msDN_DX(1,1)*ms_aux2[1]+msDN_DX(1,2)*ms_aux2[2];
    ms_aux1[2]=msDN_DX(2,0)*ms_aux2[0]+msDN_DX(2,1)*ms_aux2[1]+msDN_DX(2,2)*ms_aux2[2];
    ms_aux1[3]=msDN_DX(3,0)*ms_aux2[0]+msDN_DX(3,1)*ms_aux2[1]+msDN_DX(3,2)*ms_aux2[2];
    //KRATOS_WATCH(temp)

    rRightHandSideVector += tau*density*Volume*ms_aux1;


    //RHS += -tau*nablaN*du_gausspoint/dt

    //we reuse ms_vel_gauss to store the accelerations( (u_n - u_n-1)/dt)

    ms_vel_gauss[0]=0.25*(fv0_old[0]+fv1_old[0]+fv2_old[0]+fv3_old[0]-fv0_n_1[0]-fv1_n_1[0]-fv2_n_1[0]-fv3_n_1[0])/dt;
    ms_vel_gauss[1]=0.25*(fv0_old[1]+fv1_old[1]+fv2_old[1]+fv3_old[1]-fv0_n_1[1]-fv1_n_1[1]-fv2_n_1[1]-fv3_n_1[1])/dt;
    ms_vel_gauss[2]=0.25*(fv0_old[2]+fv1_old[2]+fv2_old[2]+fv3_old[2]-fv0_n_1[2]-fv1_n_1[2]-fv2_n_1[2]-fv3_n_1[2])/dt;

    //and now we reuse ms_aux1

    ms_aux1=prod(msDN_DX,ms_vel_gauss);

    noalias(rRightHandSideVector) -= tau*density*Volume*ms_aux1;


    //and now the stabilization term referring to the convective operator
    //RHS+=nablaq*convop (convetcion of the Gauss point vel)

    //contains d_ug/dx
    //x-component of u derived with resp to x, remember msAuxVec stores velocity
    msGrad_ug(0,0)=msDN_DX(0,0)*msAuxVec[0]+msDN_DX(1,0)*msAuxVec[3]+msDN_DX(2,0)*msAuxVec[6]+msDN_DX(3,0)*msAuxVec[9];
    //x-component of u derived with resp to y
    msGrad_ug(0,1)=msDN_DX(0,1)*msAuxVec[0]+msDN_DX(1,1)*msAuxVec[3]+msDN_DX(2,1)*msAuxVec[6]+msDN_DX(3,1)*msAuxVec[9];
    //x-component of u derived with resp to z
    msGrad_ug(0,2)=msDN_DX(0,2)*msAuxVec[0]+msDN_DX(1,2)*msAuxVec[3]+msDN_DX(2,2)*msAuxVec[6]+msDN_DX(3,2)*msAuxVec[9];


    //y-component of u derived with resp to x
    msGrad_ug(1,0)=msDN_DX(0,0)*msAuxVec[1]+msDN_DX(1,0)*msAuxVec[4]+msDN_DX(2,0)*msAuxVec[7]+msDN_DX(3,0)*msAuxVec[10];
    //y-component of u derived with resp to y
    msGrad_ug(1,1)=msDN_DX(0,1)*msAuxVec[1]+msDN_DX(1,1)*msAuxVec[4]+msDN_DX(2,1)*msAuxVec[7]+msDN_DX(3,1)*msAuxVec[10];
    //y-component of u derived with resp to z
    msGrad_ug(1,2)=msDN_DX(0,2)*msAuxVec[1]+msDN_DX(1,2)*msAuxVec[4]+msDN_DX(2,2)*msAuxVec[7]+msDN_DX(3,2)*msAuxVec[10];


    msGrad_ug(2,0)=msDN_DX(0,0)*msAuxVec[2]+msDN_DX(1,0)*msAuxVec[5]+msDN_DX(2,0)*msAuxVec[8]+msDN_DX(3,0)*msAuxVec[11];
    //y-component of u derived with resp to y
    msGrad_ug(2,1)=msDN_DX(0,1)*msAuxVec[2]+msDN_DX(1,1)*msAuxVec[5]+msDN_DX(2,1)*msAuxVec[8]+msDN_DX(3,1)*msAuxVec[11];
    //y-component of u derived with resp to z
    msGrad_ug(2,2)=msDN_DX(0,2)*msAuxVec[2]+msDN_DX(1,2)*msAuxVec[5]+msDN_DX(2,2)*msAuxVec[8]+msDN_DX(3,2)*msAuxVec[11];


    array_1d<double,3> a;
    a[0]=0.25*(msAuxVec[0]+msAuxVec[3]+msAuxVec[6]+msAuxVec[9])*msGrad_ug(0,0) +
         0.25*(msAuxVec[1]+msAuxVec[4]+msAuxVec[7]+msAuxVec[10])*msGrad_ug(0,1) +
         0.25*(msAuxVec[2]+msAuxVec[5]+msAuxVec[8]+msAuxVec[11])*msGrad_ug(0,2);

    a[1]=0.25*(msAuxVec[0]+msAuxVec[3]+msAuxVec[6]+msAuxVec[9])*msGrad_ug(1,0) +
         0.25*(msAuxVec[1]+msAuxVec[4]+msAuxVec[7]+msAuxVec[10])*msGrad_ug(1,1)+
         0.25*(msAuxVec[2]+msAuxVec[5]+msAuxVec[8]+msAuxVec[11])*msGrad_ug(1,2);

    a[2]=0.25*(msAuxVec[0]+msAuxVec[3]+msAuxVec[6]+msAuxVec[9])*msGrad_ug(2,0)+
         0.25*(msAuxVec[1]+msAuxVec[4]+msAuxVec[7]+msAuxVec[10])*msGrad_ug(2,1)+
         0.25*(msAuxVec[2]+msAuxVec[5]+msAuxVec[8]+msAuxVec[11])*msGrad_ug(2,2);


    //we again reuse ms_aux0
    noalias(ms_aux0) = prod(msDN_DX,a);

    noalias(rRightHandSideVector) -= tau*density*Volume*ms_aux0;

    KRATOS_CATCH("")
}
//************************************************************************************
//************************************************************************************


void Fluid3DGLS_expl_comp::FinalFractionalStep(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    KRATOS_THROW_ERROR(std::logic_error,  "METHOD NOT IMPL inside the element Final Fractional Step is done within the low_mach strategy.. " , "");

    KRATOS_CATCH("")
}
//************************************************************************************
//************************************************************************************






//************************************************************************************
//************************************************************************************
void Fluid3DGLS_expl_comp::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
{
    //if the VAr is NODAL_MASS, we calculate the lumped mass
    if(rVariable == NODAL_MASS)
    {
        CalculateLumpedMass();
    }
    else
        KRATOS_THROW_ERROR(std::logic_error,  "You are doing something wrong  FCT calculate... of nodal_mass with wring parameters.. " , "");
}
//************************************************************************************
//************************************************************************************
void Fluid3DGLS_expl_comp::Calculate(const Variable<array_1d<double,3> >& rVariable, array_1d<double,3>& Output, const ProcessInfo& rCurrentProcessInfo)
{

    if(rVariable == VELOCITY && Output[0]==1.0)
        //we use "Output" as a switch between 1st Frac Step and last Frac Step(Correction Step)
    {
        //here the residual will be temporarily written
        Vector TmpRhs(12);
        // first we write the Galerkin contributions to the momentum residual
        CalculateGalerkinMomentumResidual(TmpRhs);
        //and now the stabilization terms added
        double dt = rCurrentProcessInfo[DELTA_TIME];
        CalculateRHSVector(TmpRhs,  dt);

    }
    else if(rVariable == VELOCITY && Output[0]==2.0)
    {
        FinalFractionalStep(rCurrentProcessInfo);
    }
    else
    {
        KRATOS_THROW_ERROR(std::logic_error,  "You are doing something wrong in ur fractional step.... " , "");
    }


}
//************************************************************************************
//************************************************************************************
//************************************************************************************
//************************************************************************************
void Fluid3DGLS_expl_comp::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().PointsNumber();


    if(rResult.size() != number_of_nodes)
        rResult.resize(number_of_nodes,false);

    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        rResult[i] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
    }
    KRATOS_CATCH("")

}

//************************************************************************************
//************************************************************************************
void Fluid3DGLS_expl_comp::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().PointsNumber();

    if(ElementalDofList.size() != number_of_nodes)
        ElementalDofList.resize(number_of_nodes);

    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        ElementalDofList[i] = GetGeometry()[i].pGetDof(PRESSURE);
    }
    KRATOS_CATCH("");
}



} // Namespace Kratos


