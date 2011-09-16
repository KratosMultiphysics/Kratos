/* b
==============================================================================
KratosIncompressibleFluidApplication 
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
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2009-01-13 15:39:56 $
//   Revision:            $Revision: 1.14 $
//
//

//#define GRADPN_FORM

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/fluid_2d.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{

    //************************************************************************************
    //************************************************************************************

    Fluid2D::Fluid2D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    Fluid2D::Fluid2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {
    }

    Element::Pointer Fluid2D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        KRATOS_TRY

        return Element::Pointer(new Fluid2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        KRATOS_CATCH("");
    }

    Fluid2D::~Fluid2D()
    {
    }

    //************************************************************************************
    //************************************************************************************

    void Fluid2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

                int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];

        if (FractionalStepNumber < 3) //first step of the fractional step solution
        {
            int ComponentIndex = FractionalStepNumber - 1;
            Stage1(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, ComponentIndex);
        } else if (FractionalStepNumber == 4)//second step of the fractional step solution
        {
            Stage2(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
        }

        KRATOS_CATCH("")
    }

    //************************************************************************************
    //************************************************************************************

    void Fluid2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR(std::logic_error, "method not implemented", "");
    }

    //************************************************************************************
    //************************************************************************************
    //calculation by component of the fractional step velocity corresponding to the first stage

    void Fluid2D::Stage1(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo, unsigned int ComponentIndex)
    {
        KRATOS_TRY;

        const unsigned int number_of_points = 3;

        if (rLeftHandSideMatrix.size1() != number_of_points)
            rLeftHandSideMatrix.resize(number_of_points, number_of_points, false); //false says not to preserve existing storage!!

        if (rRightHandSideVector.size() != number_of_points)
            rRightHandSideVector.resize(number_of_points, false); //false says not to preserve existing storage!!

        boost::numeric::ublas::bounded_matrix<double, 3, 3 > MassFactors = 1.0 / 3.0 * IdentityMatrix(3, 3);
        boost::numeric::ublas::bounded_matrix<double, 3, 3 > WorkMatrix;
        boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
        array_1d<double, 3 > N;
        array_1d<double, 2 > vel_gauss;
        array_1d<double, 3 > temp_vec_np;
        array_1d<double, 3 > u_DN;
        array_1d<double, 3 > aux0;
        array_1d<double, 3 > aux1;
        array_1d<double, 3 > aux2;

        //getting data for the given geometry
        double Area;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);

        //getting the velocity vector on the nodes

        //getting the velocity on the nodes
        const array_1d<double, 3 > & fv0 = GetGeometry()[0].FastGetSolutionStepValue(FRACT_VEL, 0);
        const array_1d<double, 3 > & w0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
        const array_1d<double, 3 > & proj0 = GetGeometry()[0].FastGetSolutionStepValue(CONV_PROJ);
        const double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT);
        const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
        const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);
        const double fcomp0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE)[ComponentIndex];

        const array_1d<double, 3 > & fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
        const array_1d<double, 3 > & w1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
        const array_1d<double, 3 > & proj1 = GetGeometry()[1].FastGetSolutionStepValue(CONV_PROJ);
        const double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT);
        const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
        const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);
        const double fcomp1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE)[ComponentIndex];

        const array_1d<double, 3 > & fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
        const array_1d<double, 3 > & w2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);
        const array_1d<double, 3 > & proj2 = GetGeometry()[2].FastGetSolutionStepValue(CONV_PROJ);
        const double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT);
        const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
        const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);
        const double fcomp2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE)[ComponentIndex];

        //====================================================================
        //calculating viscosity
        const double nu = 0.333333333333333333333333 * (nu0 + nu1 + nu2);
        const double density = 0.3333333333333333333333 * (rho0 + rho1 + rho2);
        const double h  = CalculateH(DN_DX,Area);

        //getting the BDF2 coefficients (not fixed to allow variable time step)
        //the coefficients INCLUDE the time step
        const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

        //compute turbulent viscosity
        const double Cs = this->GetValue(C_SMAGORINSKY);
	

        double nu_turbulent = 0.0;
        if (Cs != 0.0){
            nu_turbulent = ComputeSmagorinskyViscosity(DN_DX, h, Cs, nu);
	
}


        //VISCOUS CONTRIBUTION TO THE STIFFNESS MATRIX
        // rLeftHandSideMatrix += Laplacian * nu; --> ONE GAUSS POINT
        noalias(rLeftHandSideMatrix) = (nu + nu_turbulent) * prod(DN_DX, trans(DN_DX));

        //INERTIA CONTRIBUTION
        //filling the mass factors

        //  rLeftHandSideMatrix += M/Dt
        noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * MassFactors;

        // *****************************************
        //CALCULATION OF THE RHS

        //external forces (component)
        double force_component = 0.3333333333333333 * (fcomp0 + fcomp1 + fcomp2);

#if defined(GRADPN_FORM)
        //std::cout << "grad pn form " << std::endl;
        //calculating pressure grad (component)
        double p_grad = DN_DX(0, ComponentIndex) * p0old
                + DN_DX(1, ComponentIndex) * p1old
                + DN_DX(2, ComponentIndex) * p2old;
        p_grad /= density;
        // RHS = Fext - grad*pn
        noalias(rRightHandSideVector) = (force_component - p_grad) * N;
#else 
        //adding pressure gradient (integrated by parts)
        noalias(rRightHandSideVector) = (force_component) * N;
        double p_avg = N[0] * p0old + N[1] * p1old + N[2] * p2old;
        p_avg /= density;
        rRightHandSideVector[0] += DN_DX(0, ComponentIndex) * p_avg;
        rRightHandSideVector[1] += DN_DX(1, ComponentIndex) * p_avg;
        rRightHandSideVector[2] += DN_DX(2, ComponentIndex) * p_avg;
#endif

        //adding the inertia terms
        // RHS += M*vhistory
        //calculating the historical velocity
        noalias(temp_vec_np) = ZeroVector(3);
        for (unsigned int step = 1; step < BDFcoeffs.size(); step++)
        {
            for (int iii = 0; iii < 3; iii++)
            {
                const array_1d<double, 3 > & v = (GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY, step));
                temp_vec_np[iii] -= BDFcoeffs[step] * v[ComponentIndex];
            }

        }
        noalias(rRightHandSideVector) += prod(MassFactors, temp_vec_np);

        //multiplying by area
        rLeftHandSideMatrix *= (Area * density);
        rRightHandSideVector *= (Area * density);



        //====================================================================
        //calculation of convective and stabilizing terms ... using 3 gauss points (on the sides)
        //		double c1 = 4.00; double c2 = 2.00;

        double norm_u = 0.0;
        double tau = 0.0;
        double area_density_third = Area * density * 0.33333333333333333333;

        // ******************* GAUSS1 ***************************
        N[0] = 0.5;
        N[1] = 0.5;
        N[2] = 0.0;

        // vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points)
        //(note that the fractional step vel is used)
        vel_gauss[0] = N[0]*(fv0[0] - w0[0]) + N[1]*(fv1[0] - w1[0]) + N[2]*(fv2[0] - w2[0]);
        vel_gauss[1] = N[0]*(fv0[1] - w0[1]) + N[1]*(fv1[1] - w1[1]) + N[2]*(fv2[1] - w2[1]);
        // KRATOS_WATCH(vel_gauss);

        //calculating parameter tau
        norm_u = vel_gauss[0] * vel_gauss[0] + vel_gauss[1] * vel_gauss[1];
        norm_u = sqrt(norm_u);

//         tau = CalculateTau(h, nu, norm_u, rCurrentProcessInfo);
        tau = CalculateTau(DN_DX, vel_gauss, h, nu, norm_u, rCurrentProcessInfo);

        //CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
        noalias(u_DN) = prod(DN_DX, vel_gauss);
        noalias(WorkMatrix) = outer_prod(N, u_DN);

        //CONVECTION STABILIZING CONTRIBUTION (Suu)
        noalias(WorkMatrix) += tau * outer_prod(u_DN, u_DN);

        //adding gauss point contribution
        noalias(rLeftHandSideMatrix) += area_density_third * WorkMatrix;

        //RHS += Suy * proj[component]
        double proj_component = N[0] * proj0[ComponentIndex]
                + N[1] * proj1[ComponentIndex]
                + N[2] * proj2[ComponentIndex];
        noalias(rRightHandSideVector) += (area_density_third * tau * proj_component) * u_DN;


        // ******************* GAUSS2 ***************************
        N[0] = 0.0;
        N[1] = 0.5;
        N[2] = 0.5;

        // vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points)
        //(note that the fractional step vel is used)
        vel_gauss[0] = N[0]*(fv0[0] - w0[0]) + N[1]*(fv1[0] - w1[0]) + N[2]*(fv2[0] - w2[0]);
        vel_gauss[1] = N[0]*(fv0[1] - w0[1]) + N[1]*(fv1[1] - w1[1]) + N[2]*(fv2[1] - w2[1]);
        // KRATOS_WATCH(vel_gauss);

        //calculating parameter tau
        norm_u = vel_gauss[0] * vel_gauss[0] + vel_gauss[1] * vel_gauss[1];
        norm_u = sqrt(norm_u);

//         tau = CalculateTau(h, nu, norm_u, rCurrentProcessInfo);
        tau = CalculateTau(DN_DX, vel_gauss, h, nu, norm_u, rCurrentProcessInfo);

        //CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
        noalias(u_DN) = prod(DN_DX, vel_gauss);
        noalias(WorkMatrix) = outer_prod(N, u_DN);

        //CONVECTION STABILIZING CONTRIBUTION (Suu)
        noalias(WorkMatrix) += tau * outer_prod(u_DN, u_DN);

        //adding gauss point contribution
        noalias(rLeftHandSideMatrix) += area_density_third * WorkMatrix;

        //RHS += Suy * proj[component]
        proj_component = N[0] * proj0[ComponentIndex]
                + N[1] * proj1[ComponentIndex]
                + N[2] * proj2[ComponentIndex];
        noalias(rRightHandSideVector) += (area_density_third * tau * proj_component) * u_DN;

        // ******************* GAUSS3 ***************************
        N[0] = 0.5;
        N[1] = 0.0;
        N[2] = 0.5;

        // vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points)
        //(note that the fractional step vel is used)
        vel_gauss[0] = N[0]*(fv0[0] - w0[0]) + N[1]*(fv1[0] - w1[0]) + N[2]*(fv2[0] - w2[0]);
        vel_gauss[1] = N[0]*(fv0[1] - w0[1]) + N[1]*(fv1[1] - w1[1]) + N[2]*(fv2[1] - w2[1]);
        // KRATOS_WATCH(vel_gauss);

        //calculating parameter tau
        norm_u = vel_gauss[0] * vel_gauss[0] + vel_gauss[1] * vel_gauss[1];
        norm_u = sqrt(norm_u);


//         tau = CalculateTau(h, nu, norm_u, rCurrentProcessInfo);
        tau = CalculateTau(DN_DX, vel_gauss, h, nu, norm_u, rCurrentProcessInfo);

        // KRATOS_WATCH(tau);

        //CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
        noalias(u_DN) = prod(DN_DX, vel_gauss);
        noalias(WorkMatrix) = outer_prod(N, u_DN);

        //CONVECTION STABILIZING CONTRIBUTION (Suu)
        noalias(WorkMatrix) += tau * outer_prod(u_DN, u_DN);

        //adding gauss point contribution
        noalias(rLeftHandSideMatrix) += area_density_third * WorkMatrix;

        //RHS += Suy * proj[component]
        proj_component = N[0] * proj0[ComponentIndex]
                + N[1] * proj1[ComponentIndex]
                + N[2] * proj2[ComponentIndex];
        noalias(rRightHandSideVector) += (area_density_third * tau * proj_component) * u_DN;



        //suubtracting the dirichlet term
        // RHS -= LHS*FracVel
        temp_vec_np[0] = fv0[ComponentIndex];
        temp_vec_np[1] = fv1[ComponentIndex];
        temp_vec_np[2] = fv2[ComponentIndex];
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, temp_vec_np);
        // rLeftHandSideMatrix *= density;
        // rRightHandSideVector *= density;

        // 		KRATOS_WATCH(rLeftHandSideMatrix);
        // 		KRATOS_WATCH(rRightHandSideVector);

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************
    //calculation by component of the fractional step velocity corresponding to the first stage

    void Fluid2D::Stage2(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        unsigned int number_of_points = 3;

        if (rLeftHandSideMatrix.size1() != number_of_points)
            rLeftHandSideMatrix.resize(number_of_points, number_of_points, false);

        if (rRightHandSideVector.size() != number_of_points)
            rRightHandSideVector.resize(number_of_points, false);


        //boost::numeric::ublas::bounded_matrix<double,3,3> MassFactors = 1.0/3.0*IdentityMatrix(3,3);
        boost::numeric::ublas::bounded_matrix<double, 3, 3 > WorkMatrix;
        boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
        array_1d<double, 3 > N;
        array_1d<double, 2 > vel_gauss;
        array_1d<double, 3 > temp_vec_np;
        array_1d<double, 3 > u_DN;
        double Area;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);


        const array_1d<double, 3 > & fv0 = GetGeometry()[0].FastGetSolutionStepValue(FRACT_VEL);
        const array_1d<double, 3 > & w0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
        const array_1d<double, 3 > & proj0 = GetGeometry()[0].FastGetSolutionStepValue(PRESS_PROJ);
        double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
        double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT);
        const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
        const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);

        const array_1d<double, 3 > & fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
        const array_1d<double, 3 > & w1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
        const array_1d<double, 3 > & proj1 = GetGeometry()[1].FastGetSolutionStepValue(PRESS_PROJ);
        double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
        double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT);
        const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
        const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);

        const array_1d<double, 3 > & fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
        const array_1d<double, 3 > & w2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);
        const array_1d<double, 3 > & proj2 = GetGeometry()[2].FastGetSolutionStepValue(PRESS_PROJ);
        double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
        double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT);
        const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
        const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);

        // vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) (note that the fractional step vel is used)
        vel_gauss[0] = N[0]*(fv0[0] - w0[0]) + N[1]*(fv1[0] - w1[0]) + N[2]*(fv2[0] - w2[0]);
        vel_gauss[1] = N[0]*(fv0[1] - w0[1]) + N[1]*(fv1[1] - w1[1]) + N[2]*(fv2[1] - w2[1]);

        //calculating convective auxiliary vector

        noalias(u_DN) = prod(DN_DX, vel_gauss);

        //calculating average density and viscosity
        double nu = 0.33333333333333 * (nu0 + nu1 + nu2);
        double density = 0.33333333333333 * (rho0 + rho1 + rho2);

        //calculating parameter tau (saved internally to each element)
        double h  = CalculateH(DN_DX,Area);
        double norm_u = vel_gauss[0] * vel_gauss[0] + vel_gauss[1] * vel_gauss[1];
        norm_u = sqrt(norm_u);
//         double tau = CalculateTau(h, nu, norm_u, rCurrentProcessInfo);
        double tau = CalculateTau(DN_DX, vel_gauss, h, nu, norm_u, rCurrentProcessInfo);

        //getting the BDF2 coefficients (not fixed to allow variable time step)
        //the coefficients INCLUDE the time step
        const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

        //CALCULATION OF THE LEFT HAND SIDE
        //laplacian term	       L = Dt * gradN * trans(gradN);
        //stabilization term       Spp = tau * gradN * trans(gradN);
        noalias(rLeftHandSideMatrix) = ((1.00 / BDFcoeffs[0] + tau) / density) * prod(DN_DX, trans(DN_DX));

        //calculation of the RHS
        // RHS = -G*vfrac
        double Gaux;
        Gaux = DN_DX(0, 0) * fv0[0] + DN_DX(0, 1) * fv0[1];
        Gaux += DN_DX(1, 0) * fv1[0] + DN_DX(1, 1) * fv1[1];
        Gaux += DN_DX(2, 0) * fv2[0] + DN_DX(2, 1) * fv2[1];
        rRightHandSideVector[0] = -Gaux * N[0];
        rRightHandSideVector[1] = -Gaux * N[1];
        rRightHandSideVector[2] = -Gaux * N[2];

        //attention!! changing the meaning of vel_gauss
        // RHS += Sz * proj
        //contrib of Spy*proj
        vel_gauss[0] = N[0] * proj0[0] + N[1] * proj1[0] + N[2] * proj2[0];
        vel_gauss[1] = N[0] * proj0[1] + N[1] * proj1[1] + N[2] * proj2[1];
        vel_gauss *= tau;
        noalias(rRightHandSideVector) += prod(DN_DX, vel_gauss);

        //dirichlet contribution
        temp_vec_np[0] = p0;
        temp_vec_np[1] = p1;
        temp_vec_np[2] = p2;
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, temp_vec_np);

        // RHS += dt * L * pold
        temp_vec_np[0] = p0old;
        temp_vec_np[1] = p1old;
        temp_vec_np[2] = p2old;
        noalias(vel_gauss) = prod(trans(DN_DX), temp_vec_np);
        noalias(rRightHandSideVector) += (1.00 / BDFcoeffs[0] / density) * prod(DN_DX, vel_gauss);

        //multiplicating by the area
        rLeftHandSideMatrix *= Area;
        rRightHandSideVector *= Area;

        //adding contributions to nodal areas following the corresponding lumping term
        double nodal_contrib = 0.333333333333333333333333333 * Area*density;
        //		GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
        //		GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
        //		GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;

        GetGeometry()[0].SetLock();
        GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
        GetGeometry()[0].UnSetLock();

        GetGeometry()[1].SetLock();
        GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
        GetGeometry()[1].UnSetLock();

        GetGeometry()[2].SetLock();
        GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
        GetGeometry()[2].UnSetLock();

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************
    // this subroutine calculates the nodal contributions for the explicit steps of the
    // fractional step procedure

    void Fluid2D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY
                int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];
        boost::numeric::ublas::bounded_matrix<double, 3, 3 > MassFactors = 1.0 / 3.0 * IdentityMatrix(3, 3);
        boost::numeric::ublas::bounded_matrix<double, 3, 3 > WorkMatrix;
        boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
        array_1d<double, 3 > N;
        array_1d<double, 2 > vel_gauss;
        array_1d<double, 3 > temp_vec_np;
        array_1d<double, 3 > u_DN;
        array_1d<double, 3 > aux0;
        array_1d<double, 3 > aux1;
        array_1d<double, 3 > aux2;

        //getting data for the given geometry
        double Area;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);

        if (FractionalStepNumber == 5) //calculation of stabilization terms
        {

            array_1d<double, 3 > & fv0 = GetGeometry()[0].FastGetSolutionStepValue(FRACT_VEL);
            array_1d<double, 3 > & w0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
            array_1d<double, 3 > & press_proj0 = GetGeometry()[0].FastGetSolutionStepValue(PRESS_PROJ);
            array_1d<double, 3 > & conv_proj0 = GetGeometry()[0].FastGetSolutionStepValue(CONV_PROJ);
            double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
            const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);

            array_1d<double, 3 > & fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
            array_1d<double, 3 > & w1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
            array_1d<double, 3 > & press_proj1 = GetGeometry()[1].FastGetSolutionStepValue(PRESS_PROJ);
            array_1d<double, 3 > & conv_proj1 = GetGeometry()[1].FastGetSolutionStepValue(CONV_PROJ);
            double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
            const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);

            array_1d<double, 3 > & fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
            array_1d<double, 3 > & w2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);
            array_1d<double, 3 > & press_proj2 = GetGeometry()[2].FastGetSolutionStepValue(PRESS_PROJ);
            array_1d<double, 3 > & conv_proj2 = GetGeometry()[2].FastGetSolutionStepValue(CONV_PROJ);
            double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
            const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);

            double density = 0.3333333333333333333333 * (rho0 + rho1 + rho2);

            //calculation of the pressure gradient (saved in vel_gauss)
            //note that here we calculate it "strong"
            vel_gauss[0] = DN_DX(0, 0)*(p0) + DN_DX(1, 0)*(p1) + DN_DX(2, 0)*(p2);
            vel_gauss[1] = DN_DX(0, 1)*(p0) + DN_DX(1, 1)*(p1) + DN_DX(2, 1)*(p2);
            vel_gauss *= Area;

            //press_proj += G*p
            GetGeometry()[0].SetLock();
            press_proj0[0] += N[0] * vel_gauss[0];
            press_proj0[1] += N[0] * vel_gauss[1];
            GetGeometry()[0].UnSetLock();

            GetGeometry()[1].SetLock();
            press_proj1[0] += N[1] * vel_gauss[0];
            press_proj1[1] += N[1] * vel_gauss[1];
            GetGeometry()[1].UnSetLock();

            GetGeometry()[2].SetLock();
            press_proj2[0] += N[2] * vel_gauss[0];
            press_proj2[1] += N[2] * vel_gauss[1];
            GetGeometry()[2].UnSetLock();

            //CONVECTIVE PROJECTION

            // ******************* GAUSS1 ***************************
            N[0] = 0.5;
            N[1] = 0.5;
            N[2] = 0.0;

            // vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points)
            //(note that the fractional step vel is used)
            vel_gauss[0] = N[0]*(fv0[0] - w0[0]) + N[1]*(fv1[0] - w1[0]) + N[2]*(fv2[0] - w2[0]);
            vel_gauss[1] = N[0]*(fv0[1] - w0[1]) + N[1]*(fv1[1] - w1[1]) + N[2]*(fv2[1] - w2[1]);

            //calculating convective auxiliary vector
            noalias(u_DN) = prod(DN_DX, vel_gauss);

            //attention changing the meaning of vel_gauss!!
            vel_gauss[0] = u_DN[0] * fv0[0] + u_DN[1] * fv1[0] + u_DN[2] * fv2[0];
            vel_gauss[1] = u_DN[0] * fv0[1] + u_DN[1] * fv1[1] + u_DN[2] * fv2[1];

            // conv_proj += C*u
            aux0[0] = N[0] * vel_gauss[0];
            aux0[1] = N[0] * vel_gauss[1];

            aux1[0] = N[1] * vel_gauss[0];
            aux1[1] = N[1] * vel_gauss[1];

            aux2[0] = N[2] * vel_gauss[0];
            aux2[1] = N[2] * vel_gauss[1];


            // ******************* GAUSS2 ***************************
            N[0] = 0.0;
            N[1] = 0.5;
            N[2] = 0.5;

            // vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points)
            //(note that the fractional step vel is used)
            vel_gauss[0] = N[0]*(fv0[0] - w0[0]) + N[1]*(fv1[0] - w1[0]) + N[2]*(fv2[0] - w2[0]);
            vel_gauss[1] = N[0]*(fv0[1] - w0[1]) + N[1]*(fv1[1] - w1[1]) + N[2]*(fv2[1] - w2[1]);

            //calculating convective auxiliary vector
            noalias(u_DN) = prod(DN_DX, vel_gauss);

            //attention changing the meaning of vel_gauss!!
            vel_gauss[0] = u_DN[0] * fv0[0] + u_DN[1] * fv1[0] + u_DN[2] * fv2[0];
            vel_gauss[1] = u_DN[0] * fv0[1] + u_DN[1] * fv1[1] + u_DN[2] * fv2[1];

            // conv_proj += C*u
            aux0[0] += N[0] * vel_gauss[0];
            aux0[1] += N[0] * vel_gauss[1];

            aux1[0] += N[1] * vel_gauss[0];
            aux1[1] += N[1] * vel_gauss[1];

            aux2[0] += N[2] * vel_gauss[0];
            aux2[1] += N[2] * vel_gauss[1];

            // ******************* GAUSS3 ***************************
            N[0] = 0.5;
            N[1] = 0.0;
            N[2] = 0.5;

            // vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points)
            //(note that the fractional step vel is used)
            vel_gauss[0] = N[0]*(fv0[0] - w0[0]) + N[1]*(fv1[0] - w1[0]) + N[2]*(fv2[0] - w2[0]);
            vel_gauss[1] = N[0]*(fv0[1] - w0[1]) + N[1]*(fv1[1] - w1[1]) + N[2]*(fv2[1] - w2[1]);

            //calculating convective auxiliary vector
            noalias(u_DN) = prod(DN_DX, vel_gauss);

            //attention changing the meaning of vel_gauss!!
            vel_gauss[0] = u_DN[0] * fv0[0] + u_DN[1] * fv1[0] + u_DN[2] * fv2[0];
            vel_gauss[1] = u_DN[0] * fv0[1] + u_DN[1] * fv1[1] + u_DN[2] * fv2[1];

            // conv_proj += C*u
            aux0[0] += N[0] * vel_gauss[0];
            aux0[1] += N[0] * vel_gauss[1];

            aux1[0] += N[1] * vel_gauss[0];
            aux1[1] += N[1] * vel_gauss[1];

            aux2[0] += N[2] * vel_gauss[0];
            aux2[1] += N[2] * vel_gauss[1];

            //access to database

            // conv_proj += C*u
            double temp = Area * density * 0.33333333333333333333333333;

            GetGeometry()[0].SetLock();
            conv_proj0[0] += temp * aux0[0];
            conv_proj0[1] += temp * aux0[1];
            GetGeometry()[0].UnSetLock();

            GetGeometry()[1].SetLock();
            conv_proj1[0] += temp * aux1[0];
            conv_proj1[1] += temp * aux1[1];
            GetGeometry()[1].UnSetLock();

            GetGeometry()[2].SetLock();
            conv_proj2[0] += temp * aux2[0];
            conv_proj2[1] += temp * aux2[1];
            GetGeometry()[2].UnSetLock();
        }
        else if (FractionalStepNumber == 6) //calculation of velocities
        {
            //outside of the element it is performed a loop on the elements in which it is considered nodally
            //v = vfrac - Dt/Mnodal * G*(p-pold)
            // the term G*(p-pold) needs to be assembled by element. the contribution is thus directly added to the nodal contrib

            double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
            double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT);
            array_1d<double, 3 > & fv0 = GetGeometry()[0].FastGetSolutionStepValue(FRACT_VEL);
            //const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);

            double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
            double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT);
            array_1d<double, 3 > & fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
            //const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);

            double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
            double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT);
            array_1d<double, 3 > & fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
            //const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);

            //double density = 0.33333333333333333333333*(rho0 + rho1 + rho2);


#if defined(GRADPN_FORM)
            //adding pressure gradient integrated by parts
            // fv += G*p(n+1) - grad p(n)
            double p_avg = N[0] * p0 + N[1] * p1 + N[2] * p2;
            p_avg *= Area / density;
            fv0[0] += DN_DX(0, 0) * p_avg;
            fv0[1] += DN_DX(0, 1) * p_avg;
            fv1[0] += DN_DX(1, 0) * p_avg;
            fv1[1] += DN_DX(1, 1) * p_avg;
            fv2[0] += DN_DX(2, 0) * p_avg;
            fv2[1] += DN_DX(2, 1) * p_avg;
            //fv -= grad_pn
            temp_vec_np[0] = p0old;
            temp_vec_np[1] = p1old;
            temp_vec_np[2] = p2old;
            noalias(vel_gauss) = prod(trans(DN_DX), temp_vec_np);
            vel_gauss *= (Area);
            fv0[0] += (N[0]) * vel_gauss[0];
            fv0[1] += (N[0]) * vel_gauss[1];
            fv1[0] += (N[1]) * vel_gauss[0];
            fv1[1] += (N[1]) * vel_gauss[1];
            fv2[0] += (N[2]) * vel_gauss[0];
            fv2[1] += (N[2]) * vel_gauss[1];
#else //G pn form
            double p_avg = N[0]*(p0 - p0old)
                    + N[1]*(p1 - p1old)
                    + N[2]*(p2 - p2old);
            p_avg *= Area;
            fv0[0] += DN_DX(0, 0) * p_avg;
            fv0[1] += DN_DX(0, 1) * p_avg;
            fv1[0] += DN_DX(1, 0) * p_avg;
            fv1[1] += DN_DX(1, 1) * p_avg;
            fv2[0] += DN_DX(2, 0) * p_avg;
            fv2[1] += DN_DX(2, 1) * p_avg;
#endif
        }

        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************

    void Fluid2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
    {
        unsigned int number_of_nodes = GetGeometry().PointsNumber();
        if (rResult.size() != number_of_nodes)
            rResult.resize(number_of_nodes, false);

        unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

        if (FractionalStepNumber == 1) //step 1
            for (unsigned int i = 0; i < number_of_nodes; i++)
                rResult[i] = GetGeometry()[i].GetDof(FRACT_VEL_X).EquationId();
        else if (FractionalStepNumber == 2) //step 2
            for (unsigned int i = 0; i < number_of_nodes; i++)
                rResult[i] = GetGeometry()[i].GetDof(FRACT_VEL_Y).EquationId();

        else if (FractionalStepNumber == 4) // pressure correction step
            for (unsigned int i = 0; i < number_of_nodes; i++)
                rResult[i] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
    }

    //************************************************************************************
    //************************************************************************************

    void Fluid2D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
    {
        unsigned int number_of_nodes = GetGeometry().PointsNumber();
        if (ElementalDofList.size() != number_of_nodes)
            ElementalDofList.resize(number_of_nodes);

        unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

        if (FractionalStepNumber == 1) //step 1
            for (unsigned int i = 0; i < number_of_nodes; i++)
                ElementalDofList[i] = GetGeometry()[i].pGetDof(FRACT_VEL_X);
        else if (FractionalStepNumber == 2) //step 2
            for (unsigned int i = 0; i < number_of_nodes; i++)
                ElementalDofList[i] = GetGeometry()[i].pGetDof(FRACT_VEL_Y);
        else if (FractionalStepNumber == 4) // pressure correction step
            for (unsigned int i = 0; i < number_of_nodes; i++)
                ElementalDofList[i] = GetGeometry()[i].pGetDof(PRESSURE);

    }

    //************************************************************************************
    //************************************************************************************
    double Fluid2D::CalculateTau(boost::numeric::ublas::bounded_matrix<double, 3,2> & DN_DX, array_1d<double, 2 > & vel_gauss, const double h, const double nu, const double norm_u, const ProcessInfo& CurrentProcessInfo)
    {
	
	double temp=0.0;
	
	//viscous parts
	double viscous_part=0.0;
	for(unsigned int i=0; i<3; i++)
	  for(unsigned int k=0; k<2; k++)
	    viscous_part += DN_DX(i,k)*DN_DX(i,k);
	viscous_part *= nu;
	  
	//convective parts
	array_1d<double,2> aux = ZeroVector(3);
	for(unsigned int i=0; i<3; i++)
	{
	  const array_1d<double,3>& v = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
	  double tmp = 0.0;
	  for(unsigned int k=0; k<2; k++)
	    tmp += DN_DX(i,k)*vel_gauss[k];
	  
	  for(unsigned int l=0; l<2; l++)
	    aux[l] += tmp*v[l];
	}
	aux /= (norm_u + 1e-9);
	double conv_part=norm_2(aux);
	
	const double dyn_st_beta = CurrentProcessInfo[DYNAMIC_TAU];
        const double inv_dt_coeff = CurrentProcessInfo[BDF_COEFFICIENTS][0];
        double tau = 1.00 / (dyn_st_beta * inv_dt_coeff + viscous_part + conv_part);

        return tau;
    }

    //************************************************************************************
    //************************************************************************************

    double Fluid2D::ComputeSmagorinskyViscosity(const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX,
            const double& h,
            const double& C,
            const double nu
            )
    {
        boost::numeric::ublas::bounded_matrix<double, 2, 2 > dv_dx = ZeroMatrix(2, 2);

        // Compute Symmetric Grad(u). Note that only the lower half of the matrix is filled
        for (unsigned int k = 0; k < 3; ++k)
        {
            const array_1d< double, 3 > & rNodeVel = this->GetGeometry()[k].FastGetSolutionStepValue(FRACT_VEL);
            for (unsigned int i = 0; i < 2; ++i)
            {
                for (unsigned int j = 0; j < i; ++j) // Off-diagonal
                    dv_dx(i, j) += 0.5 * (DN_DX(k, j) * rNodeVel[i] + DN_DX(k, i) * rNodeVel[j]);
                dv_dx(i, i) += DN_DX(k, i) * rNodeVel[i]; // Diagonal
            }
        }

        // Norm[ Grad(u) ]
        double NormS(0.0);
        for (unsigned int i = 0; i < 2; ++i)
        {
            for (unsigned int j = 0; j < i; ++j)
                NormS += 2.0 * dv_dx(i, j) * dv_dx(i, j); // Using symmetry, lower half terms of the matrix are added twice
            NormS += dv_dx(i, i) * dv_dx(i, i); // Diagonal terms
        }

        NormS = sqrt(NormS);

        // Total Viscosity
        return 2.0 * C * C * h * h * NormS;
    }

    //************************************************************************************
    //************************************************************************************
    void Fluid2D::Calculate(const Variable<double >& rVariable,
            double& Output,
            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        //the variable error_ratio is here the norm of the subscale velocity as computed at the level of the gauss point
        if (rVariable == ERROR_RATIO)
        {
            double Volume;
            array_1d<double, 3 > N;
            boost::numeric::ublas::bounded_matrix<double, 3,2 > DN_DX;
            GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Volume);

            //getting the velocity vector on the nodes

            //getting the velocity on the nodes
            const array_1d<double, 3 > & fv0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3 > & w0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
            const array_1d<double, 3 > & proj0 = GetGeometry()[0].FastGetSolutionStepValue(CONV_PROJ);
            const double& p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
            const array_1d<double, 3 > & press_proj0 = GetGeometry()[0].FastGetSolutionStepValue(PRESS_PROJ);
            const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
            const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);

            const array_1d<double, 3 > & fv1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3 > & w1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
            const array_1d<double, 3 > & proj1 = GetGeometry()[1].FastGetSolutionStepValue(CONV_PROJ);
            const double& p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
            const array_1d<double, 3 > & press_proj1 = GetGeometry()[1].FastGetSolutionStepValue(PRESS_PROJ);
            const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
            const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);

            const array_1d<double, 3 > & fv2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3 > & w2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);
            const array_1d<double, 3 > & proj2 = GetGeometry()[2].FastGetSolutionStepValue(CONV_PROJ);
            const double& p2 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
            const array_1d<double, 3 > & press_proj2 = GetGeometry()[2].FastGetSolutionStepValue(PRESS_PROJ);
            const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
            const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);

            // vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) (note that the fractional step vel is used)
            array_1d<double, 2 > vel_gauss;
 	    
            vel_gauss[0] = N[0] * ( fv0[0] - w0[0]) + N[1] * ( fv1[0] - w1[0]) + N[2] * ( fv2[0] - w2[0]);
            vel_gauss[1] = N[1] * ( fv0[1] - w0[1]) + N[1] * ( fv1[1] - w1[1]) + N[2] * ( fv2[1] - w2[1]);


            //calculating viscosity
            const double one_third = 1.0/3.0;
            double nu = one_third * (nu0 + nu1 + nu2 );
            double density = one_third * (rho0 + rho1 + rho2 );

            //getting the BDF2 coefficients (not fixed to allow variable time step)
            //the coefficients INCLUDE the time step
            //const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

            //calculating parameter tau (saved internally to each element)
            double h  = CalculateH(DN_DX,Volume);
            double norm_u = vel_gauss[0] * vel_gauss[0] + vel_gauss[1] * vel_gauss[1];
            norm_u = sqrt(norm_u);
//             double tau = CalculateTau(h, nu, norm_u, rCurrentProcessInfo);
            double tau = CalculateTau(DN_DX, vel_gauss, h, nu, norm_u, rCurrentProcessInfo);

            //calculating the subscale velocity assuming that it is governed by the convection part
            array_1d<double, 3 > u_DN;
            noalias(u_DN) = prod(DN_DX, vel_gauss);

            array_1d<double, 3 > vel_subscale = u_DN[0] * fv0;
            noalias(vel_subscale) += u_DN[1] * fv1;
            noalias(vel_subscale) += u_DN[2] * fv2;
            vel_subscale *= density;

            vel_subscale[0] -= one_third * (proj0[0] + proj1[0] + proj2[0] );
            vel_subscale[1] -= one_third * (proj0[1] + proj1[1] + proj2[1] );

            //pressure gradient contributions
            vel_subscale[0] += (DN_DX(0, 0) * p0 + DN_DX(1, 0) * p1 + DN_DX(2, 0) * p2 );
            vel_subscale[1] += (DN_DX(0, 1) * p0 + DN_DX(1, 1) * p1 + DN_DX(2, 1) * p2 );

            vel_subscale[0] -= 0.25 * (press_proj0[0] + press_proj1[0] + press_proj2[0] );
            vel_subscale[1] -= 0.25 * (press_proj0[1] + press_proj1[1] + press_proj2[1] );

            vel_subscale *= (tau / density);

            vel_subscale[2] = 0.0;

            Output = norm_2(vel_subscale);

        }
        KRATOS_CATCH("");

    }

    int Fluid2D::Check(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        //check the area

        //check if if is in the XY plane

        //check that no variable has zero key


        return 0;


        KRATOS_CATCH("");
    }

    inline double Fluid2D::CalculateH(boost::numeric::ublas::bounded_matrix<double, 3,2 > & DN_DX, double Volume)
    {
	double inv_h_max = 0.0;
	for(unsigned int i=0; i<3; i++)
	{
	  double inv_h = 0.0;
	  for(unsigned int k=0; k<2; k++)
	    inv_h += DN_DX(i,k)*DN_DX(i,k);
	  
	  if(inv_h > inv_h_max) inv_h_max = inv_h;
	}
	inv_h_max = sqrt(inv_h_max);
	double h = 1.0/inv_h_max;
	
        return h ;
    }

    //#undef GRADPN_FORM
} // Namespace Kratos



