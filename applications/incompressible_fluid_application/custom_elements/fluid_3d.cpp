/*
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
//   Date:                $Date: 2009-01-14 08:32:06 $
//   Revision:            $Revision: 1.11 $
//
//


// System includes 
//#define GRADPN_FORM //the grad(pn) is used instead of the G(pn) in doing the splitting

// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/fluid_3d.h"
#include "incompressible_fluid_application.h"
#include "utilities/math_utils.h" 
#include "utilities/geometry_utilities.h" 

namespace Kratos
{
    //************************************************************************************
    //************************************************************************************

    Fluid3D::Fluid3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    Fluid3D::Fluid3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {
    }

    Element::Pointer Fluid3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        KRATOS_TRY
        return Element::Pointer(new Fluid3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        KRATOS_CATCH("");
    }

    Fluid3D::~Fluid3D()
    {
    }

    //************************************************************************************
    //************************************************************************************

    void Fluid3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

                int FractionalStepNumber = rCurrentProcessInfo[FRACTIONAL_STEP];

        if (FractionalStepNumber <= 3) //first step of the fractional step solution
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

    void Fluid3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR(std::logic_error, "method not implemented", "");
    }

    //************************************************************************************
    //************************************************************************************
    //calculation by component of the fractional step velocity corresponding to the first stage

    void Fluid3D::Stage1(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo, unsigned int ComponentIndex)
    {
        KRATOS_TRY;

        unsigned int number_of_points = 4;

        if (rLeftHandSideMatrix.size1() != number_of_points)
            rLeftHandSideMatrix.resize(number_of_points, number_of_points, false);

        if (rRightHandSideVector.size() != number_of_points)
            rRightHandSideVector.resize(number_of_points, false);

        //getting data for the given geometry
        double Volume;
        boost::numeric::ublas::bounded_matrix<double, 4, 3 > DN_DX;
        array_1d<double, 4 > N;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Volume);
        //CalculateGeometryData(DN_DX,N,Volume);

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

        const array_1d<double, 3 > & fv3 = GetGeometry()[3].FastGetSolutionStepValue(FRACT_VEL);
        const array_1d<double, 3 > & w3 = GetGeometry()[3].FastGetSolutionStepValue(MESH_VELOCITY);
        const array_1d<double, 3 > & proj3 = GetGeometry()[3].FastGetSolutionStepValue(CONV_PROJ);
        const double p3old = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE_OLD_IT);
        const double nu3 = GetGeometry()[3].FastGetSolutionStepValue(VISCOSITY);
        const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);
        const double fcomp3 = GetGeometry()[3].FastGetSolutionStepValue(BODY_FORCE)[ComponentIndex];

        // vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) (note that the fractional step vel is used)
        array_1d<double, 3 > aux;
        array_1d<double, 3 > vel_gauss;
        noalias(aux) = fv0;
        noalias(aux) -= w0;
        noalias(vel_gauss) = N[0] * aux;

        noalias(aux) = fv1;
        noalias(aux) -= w1;
        noalias(vel_gauss) += N[1] * aux;

        noalias(aux) = fv2;
        noalias(aux) -= w2;
        noalias(vel_gauss) += N[2] * aux;

        noalias(aux) = fv3;
        noalias(aux) -= w3;
        noalias(vel_gauss) += N[3] * aux;

        //calculating viscosity
        double nu = 0.25 * (nu0 + nu1 + nu2 + nu3);
        double density = 0.25 * (rho0 + rho1 + rho2 + rho3);

        //getting the BDF2 coefficients (not fixed to allow variable time step)
        //the coefficients INCLUDE the time step
        const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

        //calculating parameter tau (saved internally to each element)
        //		double c1 = 4.00;
        //		double c2 = 2.00;
        double h = CalculateH(DN_DX,Volume);
        //double h = pow(6.00*Volume,0.3333333);
        double norm_u = vel_gauss[0] * vel_gauss[0] + vel_gauss[1] * vel_gauss[1] + vel_gauss[2] * vel_gauss[2];
        norm_u = sqrt(norm_u);
        double tau = CalculateTau(DN_DX, vel_gauss, h, nu, norm_u, rCurrentProcessInfo);

        //compute turbulent viscosity
        const double Cs = this->GetValue(C_SMAGORINSKY);
        double nu_turbulent = 0.0;
        if (Cs != 0.0)
            nu_turbulent = ComputeSmagorinskyViscosity(DN_DX, h, Cs, nu);

        //CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
        array_1d<double, 4 > u_DN;
        noalias(u_DN) = prod(DN_DX, vel_gauss);
        noalias(rLeftHandSideMatrix) = outer_prod(N, u_DN);

        //CONVECTION STABILIZING CONTRIBUTION (Suu)
        noalias(rLeftHandSideMatrix) += tau * outer_prod(u_DN, u_DN);

        //VISCOUS CONTRIBUTION TO THE STIFFNESS MATRIX
        // rLeftHandSideMatrix += Laplacian * nu;
        noalias(rLeftHandSideMatrix) += (nu + nu_turbulent) * prod(DN_DX, trans(DN_DX));

        //INERTIA CONTRIBUTION
        //  rLeftHandSideMatrix += M*BDFcoeffs[0]
        boost::numeric::ublas::bounded_matrix<double, 4, 4 > MassFactors = 0.25 * IdentityMatrix(4, 4);
        noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * MassFactors;

        //multiplication by the Volume
        rLeftHandSideMatrix *= (Volume * density);

        // *****************************************
        //CALCULATION OF THE RHS

        //external forces (component)
        double force_component = 0.25 * (fcomp0 + fcomp1 + fcomp2 + fcomp3);


#if defined(GRADPN_FORM)
        //std::cout << "grad pn form " << std::endl;
        //calculating pressure grad (component)
        double p_grad = DN_DX(0, ComponentIndex) * p0old
                + DN_DX(1, ComponentIndex) * p1old
                + DN_DX(2, ComponentIndex) * p2old
                + DN_DX(3, ComponentIndex) * p3old;
        p_grad /= density;
        // RHS = Fext - grad*pn
        noalias(rRightHandSideVector) = (force_component - p_grad) * N;
#else 
        //adding pressure gradient (integrated by parts)
        noalias(rRightHandSideVector) = (force_component) * N;
        double p_avg = p0old + p1old + p2old + p3old;
        p_avg *= 0.25 / density;
        rRightHandSideVector[0] += DN_DX(0, ComponentIndex) * p_avg;
        rRightHandSideVector[1] += DN_DX(1, ComponentIndex) * p_avg;
        rRightHandSideVector[2] += DN_DX(2, ComponentIndex) * p_avg;
        rRightHandSideVector[3] += DN_DX(3, ComponentIndex) * p_avg;
#endif

        //adding the inertia terms
        // RHS += M*vhistory
        //calculating the historical velocity
        array_1d<double, 4 > temp_vec_np = ZeroVector(4);
        for (unsigned int step = 1; step < BDFcoeffs.size(); step++)
        {
            for (int iii = 0; iii < 4; iii++)
            {
                const array_1d<double, 3 > & v = (GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY, step));
                temp_vec_np[iii] -= BDFcoeffs[step] * v[ComponentIndex];
            }

        }
        noalias(rRightHandSideVector) += prod(MassFactors, temp_vec_np);

        //RHS += Suy * proj[component]
        double proj_component = N[0] * proj0[ComponentIndex]
                + N[1] * proj1[ComponentIndex]
                + N[2] * proj2[ComponentIndex]
                + N[3] * proj3[ComponentIndex];
        noalias(rRightHandSideVector) += (tau * proj_component) * u_DN;

        //multiplying by Volume
        rRightHandSideVector *= (Volume * density);

        //suubtracting the dirichlet term
        // RHS -= LHS*FracVel
        temp_vec_np[0] = fv0[ComponentIndex];
        temp_vec_np[1] = fv1[ComponentIndex];
        temp_vec_np[2] = fv2[ComponentIndex];
        temp_vec_np[3] = fv3[ComponentIndex];
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, temp_vec_np);
	
	//tau2 stabilization
	const bool activate_tau2 = rCurrentProcessInfo.GetValue(ACTIVATE_TAU2);
	if(activate_tau2 == true)
	{
	    double tau2 = nu + 0.5*h*norm_u;
	    double div_v =  DN_DX(0,0)*fv0[0] + DN_DX(0,1)*fv0[1] + DN_DX(0,2)*fv0[2] +
			    DN_DX(1,0)*fv1[0] + DN_DX(1,1)*fv1[1] + DN_DX(0,2)*fv1[2] +
			    DN_DX(2,0)*fv2[0] + DN_DX(2,1)*fv2[1] + DN_DX(0,2)*fv2[2] +
			    DN_DX(3,0)*fv3[0] + DN_DX(3,1)*fv3[1] + DN_DX(0,3)*fv3[2] 			    
			    ;
	    div_v *= tau2;
	    rRightHandSideVector[0] += Volume*density*DN_DX(0,ComponentIndex)*div_v;
	    rRightHandSideVector[1] -= Volume*density*DN_DX(1,ComponentIndex)*div_v;
	    rRightHandSideVector[2] -= Volume*density*DN_DX(2,ComponentIndex)*div_v;
	    rRightHandSideVector[3] -= Volume*density*DN_DX(3,ComponentIndex)*div_v;
	}


        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************
    //calculation by component of the fractional step velocity corresponding to the first stage

    void Fluid3D::Stage2(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
            ProcessInfo& rCurrentProcessInfo)
    {
        //KRATOS_TRY;

        unsigned int number_of_points = 4;

        if (rLeftHandSideMatrix.size1() != number_of_points)
            rLeftHandSideMatrix.resize(number_of_points, number_of_points);

        if (rRightHandSideVector.size() != number_of_points)
            rRightHandSideVector.resize(number_of_points);

        array_1d<double, 3 > vel_gauss;
        array_1d<double, 3 > aux;
        array_1d<double, 4 > temp_vec_np;

        //getting data for the given geometry
        double Volume;
        array_1d<double, 4 > N;
        boost::numeric::ublas::bounded_matrix<double, 4, 3 > DN_DX;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Volume);

        const array_1d<double, 3 > & fv0 = GetGeometry()[0].FastGetSolutionStepValue(FRACT_VEL);
        const array_1d<double, 3 > & w0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
        const array_1d<double, 3 > & proj0 = GetGeometry()[0].FastGetSolutionStepValue(PRESS_PROJ);
        const double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
        const double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT);
        const double nu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
        const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);

        const array_1d<double, 3 > & fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
        const array_1d<double, 3 > & w1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
        const array_1d<double, 3 > & proj1 = GetGeometry()[1].FastGetSolutionStepValue(PRESS_PROJ);
        const double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
        const double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT);
        const double nu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
        const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);

        const array_1d<double, 3 > & fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
        const array_1d<double, 3 > & w2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);
        const array_1d<double, 3 > & proj2 = GetGeometry()[2].FastGetSolutionStepValue(PRESS_PROJ);
        const double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
        const double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT);
        const double nu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
        const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);

        const array_1d<double, 3 > & fv3 = GetGeometry()[3].FastGetSolutionStepValue(FRACT_VEL);
        const array_1d<double, 3 > & w3 = GetGeometry()[3].FastGetSolutionStepValue(MESH_VELOCITY);
        const array_1d<double, 3 > & proj3 = GetGeometry()[3].FastGetSolutionStepValue(PRESS_PROJ);
        const double p3 = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE);
        const double p3old = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE_OLD_IT);
        const double nu3 = GetGeometry()[3].FastGetSolutionStepValue(VISCOSITY);
        const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);

        // vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) (note that the fractional step vel is used)

        noalias(aux) = fv0;
        aux[0] -= w0[0];
        aux[1] -= w0[1];
        aux[2] -= w0[2];
        noalias(vel_gauss) = N[0] * aux;

        noalias(aux) = fv1;
        aux[0] -= w1[0];
        aux[1] -= w1[1];
        aux[2] -= w1[2];
        noalias(vel_gauss) += N[1] * aux;

        noalias(aux) = fv2;
        aux[0] -= w2[0];
        aux[1] -= w2[1];
        aux[2] -= w2[2];
        noalias(vel_gauss) += N[2] * aux;

        noalias(aux) = fv3;
        aux[0] -= w3[0];
        aux[1] -= w3[1];
        aux[2] -= w3[2];
        noalias(vel_gauss) += N[3] * aux;

        //calculating avergage density and viscosity
        double nu = 0.25 * (nu0 + nu1 + nu2 + nu3);
        double density = 0.25 * (rho0 + rho1 + rho2 + rho3);

        //calculating parameter tau
        double h = CalculateH(DN_DX,Volume);
        double norm_u = vel_gauss[0] * vel_gauss[0] + vel_gauss[1] * vel_gauss[1] + vel_gauss[2] * vel_gauss[2];
        norm_u = sqrt(norm_u);
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
        Gaux = DN_DX(0, 0) * fv0[0] + DN_DX(0, 1) * fv0[1] + DN_DX(0, 2) * fv0[2];
        Gaux += DN_DX(1, 0) * fv1[0] + DN_DX(1, 1) * fv1[1] + DN_DX(1, 2) * fv1[2];
        Gaux += DN_DX(2, 0) * fv2[0] + DN_DX(2, 1) * fv2[1] + DN_DX(2, 2) * fv2[2];
        Gaux += DN_DX(3, 0) * fv3[0] + DN_DX(3, 1) * fv3[1] + DN_DX(3, 2) * fv3[2];

        //**************************** EXPERIMENTAL *****************************
        //**************************** EXPERIMENTAL *****************************
        //try ... should conserve the mass much better!
        // 		const array_1d<double,3>& v0old = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY,1);
        // 		const array_1d<double,3>& v1old = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY,1);
        // 		const array_1d<double,3>& v2old = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY,1);
        // 		const array_1d<double,3>& v3old = GetGeometry()[3].FastGetSolutionStepValue(VELOCITY,1);
        // 		Gaux +=  DN_DX(0,0)*v0old[0] + DN_DX(0,1)*v0old[1] + DN_DX(0,2)*v0old[2];
        // 		Gaux += DN_DX(1,0)*v1old[0] + DN_DX(1,1)*v1old[1] + DN_DX(1,2)*v1old[2];
        // 		Gaux += DN_DX(2,0)*v2old[0] + DN_DX(2,1)*v2old[1] + DN_DX(2,2)*v2old[2];
        // 		Gaux += DN_DX(3,0)*v3old[0] + DN_DX(3,1)*v3old[1] + DN_DX(3,2)*v3old[2];
        //**************************** EXPERIMENTAL *****************************
        //**************************** EXPERIMENTAL *****************************

        rRightHandSideVector[0] = -Gaux * N[0];
        rRightHandSideVector[1] = -Gaux * N[1];
        rRightHandSideVector[2] = -Gaux * N[2];
        rRightHandSideVector[3] = -Gaux * N[3];
        //std::cout << Id() << " Gtrans fv " << rRightHandSideVector << std::endl;

        //attention!! changing the meaning of vel_gauss
        // RHS += Sz * proj
        //contrib of Spy*proj
        noalias(vel_gauss) = N[0] * proj0;
        noalias(vel_gauss) += N[1] * proj1;
        noalias(vel_gauss) += N[2] * proj2;
        noalias(vel_gauss) += N[3] * proj3;
        vel_gauss *= tau;
        noalias(rRightHandSideVector) += prod(DN_DX, vel_gauss);

        //dirichlet contribution
        temp_vec_np[0] = p0;
        temp_vec_np[1] = p1;
        temp_vec_np[2] = p2;
        temp_vec_np[3] = p3;
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, temp_vec_np);

        // RHS += dt * L * pold
        temp_vec_np[0] = p0old;
        temp_vec_np[1] = p1old;
        temp_vec_np[2] = p2old;
        temp_vec_np[3] = p3old;
        noalias(vel_gauss) = prod(trans(DN_DX), temp_vec_np);
        noalias(rRightHandSideVector) += (1.00 / (BDFcoeffs[0] * density)) * prod(DN_DX, vel_gauss);

        //multiplicating by the Volume
        rLeftHandSideMatrix *= Volume;
        rRightHandSideVector *= Volume;

        //adding contributions to nodal Volumes following the corresponding lumping term
        double nodal_contrib = 0.25 * Volume * density;

        //                GetGeometry()[0].SetLock();
        //                GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
        //                GetGeometry()[0].UnSetLock();
        //
        //                GetGeometry()[1].SetLock();
        //                GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
        //                GetGeometry()[1].UnSetLock();
        //
        //                GetGeometry()[2].SetLock();
        //                GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
        //                GetGeometry()[2].UnSetLock();
        //
        //                GetGeometry()[3].SetLock();
        //                GetGeometry()[3].FastGetSolutionStepValue(NODAL_MASS) += nodal_contrib;
        //                GetGeometry()[3].UnSetLock();

        double& m0 = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS);
#pragma omp atomic
        m0 += nodal_contrib;

        double& m1 = GetGeometry()[1].FastGetSolutionStepValue(NODAL_MASS);
#pragma omp atomic
        m1 += nodal_contrib;

        double& m2 = GetGeometry()[2].FastGetSolutionStepValue(NODAL_MASS);
#pragma omp atomic
        m2 += nodal_contrib;

        double& m3 = GetGeometry()[3].FastGetSolutionStepValue(NODAL_MASS);
#pragma omp atomic
        m3 += nodal_contrib;


        //KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************
    // this subroutine calculates the nodal contributions for the explicit steps of the
    // fractional step procedure

    void Fluid3D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY
                int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

        //getting data for the given geometry
        double Volume;
        array_1d<double, 4 > N;
        boost::numeric::ublas::bounded_matrix<double, 4, 3 > DN_DX;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Volume);
        //CalculateGeometryData(DN_DX,N,Volume);

        if (FractionalStepNumber == 5) //calculation of stabilization terms
        {

            const array_1d<double, 3 > & fv0 = GetGeometry()[0].FastGetSolutionStepValue(FRACT_VEL);
            const array_1d<double, 3 > & w0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
            array_1d<double, 3 > & press_proj0 = GetGeometry()[0].FastGetSolutionStepValue(PRESS_PROJ);
            array_1d<double, 3 > & conv_proj0 = GetGeometry()[0].FastGetSolutionStepValue(CONV_PROJ);
            double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
            const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);

            const array_1d<double, 3 > & fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
            const array_1d<double, 3 > & w1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
            array_1d<double, 3 > & press_proj1 = GetGeometry()[1].FastGetSolutionStepValue(PRESS_PROJ);
            array_1d<double, 3 > & conv_proj1 = GetGeometry()[1].FastGetSolutionStepValue(CONV_PROJ);
            double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
            const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);

            const array_1d<double, 3 > & fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
            const array_1d<double, 3 > & w2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);
            array_1d<double, 3 > & press_proj2 = GetGeometry()[2].FastGetSolutionStepValue(PRESS_PROJ);
            array_1d<double, 3 > & conv_proj2 = GetGeometry()[2].FastGetSolutionStepValue(CONV_PROJ);
            double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
            const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);

            const array_1d<double, 3 > & fv3 = GetGeometry()[3].FastGetSolutionStepValue(FRACT_VEL);
            const array_1d<double, 3 > & w3 = GetGeometry()[3].FastGetSolutionStepValue(MESH_VELOCITY);
            array_1d<double, 3 > & press_proj3 = GetGeometry()[3].FastGetSolutionStepValue(PRESS_PROJ);
            array_1d<double, 3 > & conv_proj3 = GetGeometry()[3].FastGetSolutionStepValue(CONV_PROJ);
            double p3 = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE);
            const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);

            double density = 0.25 * (rho0 + rho1 + rho2 + rho3);

            //calculation of the pressure gradient (saved in vel_gauss)
            array_1d<double, 3 > temp;
            array_1d<double, 4 > temp_vec_np;
            temp_vec_np[0] = p0;
            temp_vec_np[1] = p1;
            temp_vec_np[2] = p2;
            temp_vec_np[3] = p3;
            noalias(temp) = prod(trans(DN_DX), temp_vec_np);
            // 			vel_gauss *= Volume/density;
            temp *= Volume * 0.25;

            //press_proj += G*p
            //			noalias(press_proj0) += N[0]*vel_gauss;
            //			noalias(press_proj1) += N[1]*vel_gauss;
            //			noalias(press_proj2) += N[2]*vel_gauss;
            //			noalias(press_proj3) += N[3]*vel_gauss;

            ////                        GetGeometry()[0].SetLock();
            ////                        noalias(press_proj0) += N[0]*vel_gauss;
            ////                        GetGeometry()[0].UnSetLock();
            ////
            ////                        GetGeometry()[1].SetLock();
            ////                        noalias(press_proj1) += N[1]*vel_gauss;
            ////                        GetGeometry()[1].UnSetLock();
            ////
            ////                        GetGeometry()[2].SetLock();
            ////                        noalias(press_proj2) += N[2]*vel_gauss;
            ////                        GetGeometry()[2].UnSetLock();
            ////
            ////                        GetGeometry()[3].SetLock();
            ////                        noalias(press_proj3) += N[3]*vel_gauss;
            ////                        GetGeometry()[3].UnSetLock();

            //                        #pragma omp atomic
            //                        press_proj0[0] += N[0]*vel_gauss[0];
            //                        #pragma omp atomic
            //                        press_proj0[1] += N[0]*vel_gauss[1];
            //                        #pragma omp atomic
            //                        press_proj0[2] += N[0]*vel_gauss[2];
            //
            //                        #pragma omp atomic
            //                        press_proj1[0] += N[1]*vel_gauss[0];
            //                        #pragma omp atomic
            //                        press_proj1[1] += N[1]*vel_gauss[1];
            //                        #pragma omp atomic
            //                        press_proj1[2] += N[1]*vel_gauss[2];
            //
            //                        #pragma omp atomic
            //                        press_proj2[0] += N[2]*vel_gauss[0];
            //                        #pragma omp atomic
            //                        press_proj2[1] += N[2]*vel_gauss[1];
            //                        #pragma omp atomic
            //                        press_proj2[2] += N[2]*vel_gauss[2];
            //
            //                        #pragma omp atomic
            //                        press_proj3[0] += N[3]*vel_gauss[0];
            //                        #pragma omp atomic
            //                        press_proj3[1] += N[3]*vel_gauss[1];
            //                        #pragma omp atomic
            //                        press_proj3[2] += N[3]*vel_gauss[2];

            // vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) (note that the fractional step vel is used)
            array_1d<double, 3 > vel_gauss;
            array_1d<double, 3 > aux;
            noalias(aux) = fv0;
            aux[0] -= w0[0];
            aux[1] -= w0[1];
            aux[2] -= w0[2];
            noalias(vel_gauss) = N[0] * aux;

            noalias(aux) = fv1;
            aux[0] -= w1[0];
            aux[1] -= w1[1];
            aux[2] -= w1[2];
            noalias(vel_gauss) += N[1] * aux;

            noalias(aux) = fv2;
            aux[0] -= w2[0];
            aux[1] -= w2[1];
            aux[2] -= w2[2];
            noalias(vel_gauss) += N[2] * aux;

            noalias(aux) = fv3;
            aux[0] -= w3[0];
            aux[1] -= w3[1];
            aux[2] -= w3[2];
            noalias(vel_gauss) += N[3] * aux;


            //calculating convective auxiliary vector
            array_1d<double, 4 > u_DN;
            noalias(u_DN) = prod(DN_DX, vel_gauss);

            //attention changing the meaning of vel_gauss!!
            noalias(vel_gauss) = u_DN[0] * fv0;
            noalias(vel_gauss) += u_DN[1] * fv1;
            noalias(vel_gauss) += u_DN[2] * fv2;
            noalias(vel_gauss) += u_DN[3] * fv3;
            // 			vel_gauss *= Volume;
            vel_gauss *= Volume * density;

            // conv_proj += C*u
            GetGeometry()[0].SetLock();
            noalias(conv_proj0) += N[0] * vel_gauss;
            noalias(press_proj0) += temp;
            GetGeometry()[0].UnSetLock();

            GetGeometry()[1].SetLock();
            noalias(conv_proj1) += N[1] * vel_gauss;
            noalias(press_proj1) += temp;
            GetGeometry()[1].UnSetLock();

            GetGeometry()[2].SetLock();
            noalias(conv_proj2) += N[2] * vel_gauss;
            noalias(press_proj2) += temp;
            GetGeometry()[2].UnSetLock();

            GetGeometry()[3].SetLock();
            noalias(conv_proj3) += N[3] * vel_gauss;
            noalias(press_proj3) += temp;
            GetGeometry()[3].UnSetLock();




        }
        else if (FractionalStepNumber == 6) //calculation of velocities
        {
            //outside of the element it is performed a loop on the elements in which it is considered nodally
            //v = vfrac - Dt/Mnodal * G*(p-pold)
            // the term G*(p-pold) needs to be assembled by element. the contribution is thus directly added to the nodal contrib

            double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
            double p0old = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE_OLD_IT);
            array_1d<double, 3 > & fv0 = GetGeometry()[0].FastGetSolutionStepValue(FRACT_VEL);

            double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
            double p1old = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE_OLD_IT);
            array_1d<double, 3 > & fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);

            double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
            double p2old = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE_OLD_IT);
            array_1d<double, 3 > & fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);

            double p3 = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE);
            double p3old = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE_OLD_IT);
            array_1d<double, 3 > & fv3 = GetGeometry()[3].FastGetSolutionStepValue(FRACT_VEL);

#if defined(GRADPN_FORM)
            //adding pressure gradient integrated by parts
            // fv += G*p(n+1) - grad p(n)
            double p_avg = N[0] * p0 + N[1] * p1 + N[2] * p2 + N[3] * p3;
            // 			p_avg *= Volume/density;
            p_avg *= Volume;
            noalias(fv0) += DN_DX.row(0) * p_avg;
            noalias(fv1) += DN_DX.row(1) * p_avg;
            noalias(fv2) += DN_DX.row(2) * p_avg;
            noalias(fv3) += DN_DX.row(3) * p_avg;
            //fv -= grad_pn
            temp_vec_np[0] = p0old;
            temp_vec_np[1] = p1old;
            temp_vec_np[2] = p2old;
            temp_vec_np[3] = p3old;
            noalias(vel_gauss) = prod(trans(DN_DX), temp_vec_np);
            // 			vel_gauss *= (Volume/density);
            vel_gauss *= (Volume);
            noalias(fv0) += (N[0]) * vel_gauss;
            noalias(fv1) += (N[1]) * vel_gauss;
            noalias(fv2) += (N[2]) * vel_gauss;
            noalias(fv3) += (N[3]) * vel_gauss;
#else //G pn form

            double p_avg = N[0]*(p0 - p0old)
                    + N[1]*(p1 - p1old)
                    + N[2]*(p2 - p2old)
                    + N[3]*(p3 - p3old);
            // 			p_avg *= Volume/density;
            p_avg *= Volume;
            noalias(fv0) += p_avg * row(DN_DX, 0);
            noalias(fv1) += p_avg * row(DN_DX, 1);
            noalias(fv2) += p_avg * row(DN_DX, 2);
            noalias(fv3) += p_avg * row(DN_DX, 3);
#endif
        }


        KRATOS_CATCH("");
    }

    //************************************************************************************
    //************************************************************************************

    void Fluid3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
    {
        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        if (rResult.size() != number_of_nodes)
            rResult.resize(number_of_nodes);

        int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

        if (FractionalStepNumber == 1) //step 1
            for (unsigned int i = 0; i < number_of_nodes; i++)
                rResult[i] = GetGeometry()[i].GetDof(FRACT_VEL_X).EquationId();
        else if (FractionalStepNumber == 2) //step 2
            for (unsigned int i = 0; i < number_of_nodes; i++)
                rResult[i] = GetGeometry()[i].GetDof(FRACT_VEL_Y).EquationId();
        else if (FractionalStepNumber == 3) //step 3
            for (unsigned int i = 0; i < number_of_nodes; i++)
                rResult[i] = GetGeometry()[i].GetDof(FRACT_VEL_Z).EquationId();

        else if (FractionalStepNumber == 4) // pressure correction step
            for (unsigned int i = 0; i < number_of_nodes; i++)
                rResult[i] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
    }

    //************************************************************************************
    //************************************************************************************

    void Fluid3D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
    {
        const unsigned int number_of_nodes = GetGeometry().PointsNumber();
        if (ElementalDofList.size() != number_of_nodes)
            ElementalDofList.resize(number_of_nodes);

        unsigned int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

        if (FractionalStepNumber == 1) //step 1
            for (unsigned int i = 0; i < number_of_nodes; i++)
                ElementalDofList[i] = GetGeometry()[i].pGetDof(FRACT_VEL_X);
        else if (FractionalStepNumber == 2) //step 2
            for (unsigned int i = 0; i < number_of_nodes; i++)
                ElementalDofList[i] = GetGeometry()[i].pGetDof(FRACT_VEL_Y);
        else if (FractionalStepNumber == 3) //step 2
            for (unsigned int i = 0; i < number_of_nodes; i++)
                ElementalDofList[i] = GetGeometry()[i].pGetDof(FRACT_VEL_Z);
        else if (FractionalStepNumber == 4) // pressure correction step
            for (unsigned int i = 0; i < number_of_nodes; i++)
                ElementalDofList[i] = GetGeometry()[i].pGetDof(PRESSURE);

    }

    inline double Fluid3D::CalculateH(boost::numeric::ublas::bounded_matrix<double, 4, 3 > & DN_DX, double Volume)
    {
//         double h = pow(6.00 * Volume, 0.3333333);

	double inv_h_max = 0.0;
	for(unsigned int i=0; i<4; i++)
	{
	  double inv_h = 0.0;
	  for(unsigned int k=0; k<3; k++)
	    inv_h += DN_DX(i,k)*DN_DX(i,k);
	  
	  if(inv_h > inv_h_max) inv_h_max = inv_h;
	}
	inv_h_max = sqrt(inv_h_max);
	double h = 1.0/inv_h_max;
	
        return h ;
    }

    //************************************************************************************
    //************************************************************************************

    inline double Fluid3D::CalculateTau(boost::numeric::ublas::bounded_matrix<double, 4, 3 > & DN_DX, array_1d<double, 3 > & vel_gauss, const double h, const double nu, const double norm_u, const ProcessInfo& CurrentProcessInfo)
    {
/*        const double c1 = 4.00;
        const double c2 = 2.00;

        const double dyn_st_beta = CurrentProcessInfo[DYNAMIC_TAU];
        const double inv_dt_coeff = CurrentProcessInfo[BDF_COEFFICIENTS][0];
        double tau = 1.00 / (dyn_st_beta * inv_dt_coeff + c1 * nu / (h * h) + c2 * norm_u / h);*/
	
	
	
	
	
	double temp=0.0;
	
	//viscous parts
	double viscous_part=0.0;
	for(unsigned int i=0; i<4; i++)
	  for(unsigned int k=0; k<3; k++)
	    viscous_part += DN_DX(i,k)*DN_DX(i,k);
	viscous_part *= nu;
	  
	//convective part
	array_1d<double,3> aux = ZeroVector(3);
	for(unsigned int i=0; i<4; i++)
	{
	  const array_1d<double,3>& v = GetGeometry()[i].FastGetSolutionStepValue(FRACT_VEL);
	  double tmp = 0.0;
	  for(unsigned int k=0; k<3; k++)
	    tmp += DN_DX(i,k)*vel_gauss[k];
	  
	  for(unsigned int l=0; l<3; l++)
	    aux[l] += tmp*v[l];
	}
	aux /= (norm_u + 1e-9);
	double conv_part=norm_2(aux);
	
	const double dyn_st_beta = CurrentProcessInfo[DYNAMIC_TAU];
        const double inv_dt_coeff = CurrentProcessInfo[BDF_COEFFICIENTS][0];
        double tau = 1.00 / (dyn_st_beta * inv_dt_coeff + viscous_part + conv_part);

        // 	      boost::numeric::ublas::bounded_matrix<double,4,3> aux;
        // 	      for(unsigned int i=0; i<4; i++)
        // 	      {
        // 		array_1d<double,3>& vv = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
        // 		aux(i,0) = vv[0];
        // 		aux(i,1) = vv[1];
        //  		aux(i,2) = vv[2];
        // 	      }
        // 	      boost::numeric::ublas::bounded_matrix<double,3,3> grad_u = prod(trans(DN_DX),aux);
        // 	      array_1d<double,3> uDu = prod(grad_u, vel_gauss);
        // 	      double uDu_norm = norm_2(uDu);
        //
        // 	      double nn = 0;
        // 	      for(unsigned int i=0; i<3; i++)
        // 		for(unsigned int j=0; j<3; j++)
        // 		  nn += grad_u(i,j)*grad_u(i,j);
        //
        // 	      double tau = norm_u*norm_u/(dyn_st_beta*inv_dt_coeff*norm_u*norm_u + c1*nu*(nn) + c2*norm_u*uDu_norm + 1e-10);

        //	      const double dyn_st_beta = CurrentProcessInfo[DYNAMIC_TAU];
        //              const double inv_dt_coeff = CurrentProcessInfo[BDF_COEFFICIENTS][0];
        //	      double viscous_part = 0.0;
        //	      double convective_part=0.0;
        // 	      for(unsigned int i=0; i<4; i++)
        //	      {
        //		  double aaa = 0.0;
        //		  for(unsigned int j=0; j<3; j++)
        //		  {
        //		      viscous_part    += DN_DX(i,j)*DN_DX(i,j);
        //		      aaa += /*fabs*/(vel_gauss[j]*DN_DX(i,j));
        //		  }
        //		  convective_part += fabs(aaa);
        //	      }
        ///*	      viscous_part *= 4.0;
        //	      convective_part *= 0.5;*/
        //	      viscous_part *= 16.0;
        //// 	      convective_part = fabs(convective_part);
        //	      double tau = 1.00 / (dyn_st_beta*inv_dt_coeff +  nu*viscous_part +  convective_part);




        return tau;

    }

    //************************************************************************************
    //************************************************************************************
    void Fluid3D::Calculate(const Variable<double >& rVariable,
            double& Output,
            const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        //the variable error_ratio is here the norm of the subscale velocity as computed at the level of the gauss point
        if (rVariable == ERROR_RATIO)
        {
            double Volume;
            array_1d<double, 4 > N;
            boost::numeric::ublas::bounded_matrix<double, 4, 3 > DN_DX;
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

            const array_1d<double, 3 > & fv3 = GetGeometry()[3].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3 > & w3 = GetGeometry()[3].FastGetSolutionStepValue(MESH_VELOCITY);
            const array_1d<double, 3 > & proj3 = GetGeometry()[3].FastGetSolutionStepValue(CONV_PROJ);
            const double& p3 = GetGeometry()[3].FastGetSolutionStepValue(PRESSURE);
            const array_1d<double, 3 > & press_proj3 = GetGeometry()[3].FastGetSolutionStepValue(PRESS_PROJ);
            const double nu3 = GetGeometry()[3].FastGetSolutionStepValue(VISCOSITY);
            const double rho3 = GetGeometry()[3].FastGetSolutionStepValue(DENSITY);

            // vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) (note that the fractional step vel is used)
            array_1d<double, 3 > vel_gauss;
            array_1d<double, 3 > aux;
            noalias(aux) = fv0;
            noalias(aux) -= w0;
            noalias(vel_gauss) = N[0] * aux;

            noalias(aux) = fv1;
            noalias(aux) -= w1;
            noalias(vel_gauss) += N[1] * aux;

            noalias(aux) = fv2;
            noalias(aux) -= w2;
            noalias(vel_gauss) += N[2] * aux;

            noalias(aux) = fv3;
            noalias(aux) -= w3;
            noalias(vel_gauss) += N[3] * aux;

            //calculating viscosity
            double nu = 0.25 * (nu0 + nu1 + nu2 + nu3);
            double density = 0.25 * (rho0 + rho1 + rho2 + rho3);

            //getting the BDF2 coefficients (not fixed to allow variable time step)
            //the coefficients INCLUDE the time step
            //const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

            //calculating parameter tau (saved internally to each element)
            double h = CalculateH(DN_DX,Volume);
            double norm_u = vel_gauss[0] * vel_gauss[0] + vel_gauss[1] * vel_gauss[1] + vel_gauss[2] * vel_gauss[2];
            norm_u = sqrt(norm_u);
            double tau = CalculateTau(DN_DX, vel_gauss, h, nu, norm_u, rCurrentProcessInfo);

            //calculating the subscale velocity assuming that it is governed by the convection part
            array_1d<double, 4 > u_DN;
            noalias(u_DN) = prod(DN_DX, vel_gauss);

            array_1d<double, 3 > vel_subscale = u_DN[0] * fv0;
            noalias(vel_subscale) += u_DN[1] * fv1;
            noalias(vel_subscale) += u_DN[2] * fv2;
            noalias(vel_subscale) += u_DN[3] * fv3;
            vel_subscale *= density;

            vel_subscale[0] -= 0.25 * (proj0[0] + proj1[0] + proj2[0] + proj3[0]);
            vel_subscale[1] -= 0.25 * (proj0[1] + proj1[1] + proj2[1] + proj3[1]);
            vel_subscale[2] -= 0.25 * (proj0[2] + proj1[2] + proj2[2] + proj3[2]);

            //pressure gradient contributions
            vel_subscale[0] += (DN_DX(0, 0) * p0 + DN_DX(1, 0) * p1 + DN_DX(2, 0) * p2 + DN_DX(3, 0) * p3);
            vel_subscale[1] += (DN_DX(0, 1) * p0 + DN_DX(1, 1) * p1 + DN_DX(2, 1) * p2 + DN_DX(3, 1) * p3);
            vel_subscale[2] += (DN_DX(0, 2) * p0 + DN_DX(1, 2) * p1 + DN_DX(2, 2) * p2 + DN_DX(3, 2) * p3);
            vel_subscale[0] -= 0.25 * (press_proj0[0] + press_proj1[0] + press_proj2[0] + press_proj3[0]);
            vel_subscale[1] -= 0.25 * (press_proj0[1] + press_proj1[1] + press_proj2[1] + press_proj3[1]);
            vel_subscale[2] -= 0.25 * (press_proj0[2] + press_proj1[2] + press_proj2[2] + press_proj3[2]);

            vel_subscale *= (tau / density);

            Output = norm_2(vel_subscale);

        }
        KRATOS_CATCH("");

    }

    //************************************************************************************
    //************************************************************************************

    double Fluid3D::ComputeSmagorinskyViscosity(const boost::numeric::ublas::bounded_matrix<double, 4, 3 > & DN_DX,
            const double& h,
            const double& C,
            const double nu
            )
    {
        boost::numeric::ublas::bounded_matrix<double, 3, 3 > dv_dx = ZeroMatrix(3, 3);

        const unsigned int nnodes = 4;
        // Compute Symmetric Grad(u). Note that only the lower half of the matrix is filled
        for (unsigned int k = 0; k < nnodes; ++k)
        {
            const array_1d< double, 3 > & rNodeVel = this->GetGeometry()[k].FastGetSolutionStepValue(FRACT_VEL);
            for (unsigned int i = 0; i < 3; ++i)
            {
                for (unsigned int j = 0; j < i; ++j) // Off-diagonal
                    dv_dx(i, j) += 0.5 * (DN_DX(k, j) * rNodeVel[i] + DN_DX(k, i) * rNodeVel[j]);
                dv_dx(i, i) += DN_DX(k, i) * rNodeVel[i]; // Diagonal
            }
        }

        // Norm[ Grad(u) ]
        double NormS(0.0);
        for (unsigned int i = 0; i < 3; ++i)
        {
            for (unsigned int j = 0; j < i; ++j)
                NormS += 2.0 * dv_dx(i, j) * dv_dx(i, j); // Using symmetry, lower half terms of the matrix are added twice
            NormS += dv_dx(i, i) * dv_dx(i, i); // Diagonal terms
        }

        NormS = sqrt(NormS);

        // Total Viscosity
        return 2.0 * C * C * h * h * NormS;
    }

    int Fluid3D::Check(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

        int check_result = Element::Check(rCurrentProcessInfo);
        if(check_result != 0)
            return 1;

        if(VELOCITY.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"VELOCITY has Key zero! (check if the application is correctly registered","");
        if(FRACT_VEL.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"FRACT_VEL has Key zero! (check if the application is correctly registered","");
        if(MESH_VELOCITY.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"MESH_VELOCITY has Key zero! (check if the application is correctly registered","");
        if(PRESSURE.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"PRESSURE has Key zero! (check if the application is correctly registered","");
        if(DENSITY.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"DENSITY has Key zero! (check if the application is correctly registered","");
        if(VISCOSITY.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"VISCOSITY has Key zero! (check if the application is correctly registered","");
        if(CONV_PROJ.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"CONV_PROJ has Key zero! (check if the application is correctly registered","");
        if(PRESS_PROJ.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"PRESS_PROJ has Key zero! (check if the application is correctly registered","");
        if(BODY_FORCE.Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"BODY_FORCE has Key zero! (check if the application is correctly registered","");




        for(unsigned int i=0; i<this->GetGeometry().size(); i++)
        {
            
            if(GetGeometry()[i].SolutionStepsDataHas(VELOCITY) == false)
                KRATOS_ERROR(std::invalid_argument,"missing variable VELOCITY","");
            if(GetGeometry()[i].HasDofFor(VELOCITY_X) == false || GetGeometry()[i].HasDofFor(VELOCITY_Y) == false || GetGeometry()[i].HasDofFor(VELOCITY_Z) == false)
                KRATOS_ERROR(std::invalid_argument,"missing dof VELOCITY","");

            if(GetGeometry()[i].SolutionStepsDataHas(FRACT_VEL) == false)
                KRATOS_ERROR(std::invalid_argument,"missing variable FRACT_VEL","");
            if(GetGeometry()[i].HasDofFor(FRACT_VEL_X) == false || GetGeometry()[i].HasDofFor(FRACT_VEL_Y) == false || GetGeometry()[i].HasDofFor(FRACT_VEL_Z) == false)
                KRATOS_ERROR(std::invalid_argument,"missing dof FRACT_VEL","");
            
            if(GetGeometry()[i].SolutionStepsDataHas(PRESSURE) == false)
                KRATOS_ERROR(std::invalid_argument,"missing variable PRESSURE","");
            if(GetGeometry()[i].HasDofFor(PRESSURE) == false)
                KRATOS_ERROR(std::invalid_argument,"missing dof PRESSURE","");

            if(GetGeometry()[i].SolutionStepsDataHas(DENSITY) == false)
                KRATOS_ERROR(std::invalid_argument,"missing variable DENSITY","");
            if(GetGeometry()[i].SolutionStepsDataHas(VISCOSITY) == false)
                KRATOS_ERROR(std::invalid_argument,"missing variable VISCOSITY","");
            if(GetGeometry()[i].SolutionStepsDataHas(MESH_VELOCITY) == false)
                KRATOS_ERROR(std::invalid_argument,"missing variable FRACT_VEL","");
            if(GetGeometry()[i].SolutionStepsDataHas(FRACT_VEL) == false)
                KRATOS_ERROR(std::invalid_argument,"missing variable FRACT_VEL","");
            if(GetGeometry()[i].SolutionStepsDataHas(PRESS_PROJ) == false)
                KRATOS_ERROR(std::invalid_argument,"missing variable PRESS_PROJ","")
            if(GetGeometry()[i].SolutionStepsDataHas(CONV_PROJ) == false)
                KRATOS_ERROR(std::invalid_argument,"missing variable CONV_PROJ","")
            if(GetGeometry()[i].SolutionStepsDataHas(BODY_FORCE) == false)
                KRATOS_ERROR(std::invalid_argument,"missing variable CONV_PROJ","")
        }

        return 0;
        KRATOS_CATCH("");
    }

} // Namespace Kratos


