/*
==============================================================================
KratosConvectionDiffusionApplication 
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
//   Date:                $Date: 2007-03-06 10:30:32 $
//   Revision:            $Revision: 1.4 $
//
//


// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/spalart_allmaras_element.h"
#include "fluid_dynamics_application_variables.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "spalart_allmaras_element.h"
#include "includes/element.h"
#include "includes/node.h"

namespace Kratos
{

    //************************************************************************************
    //************************************************************************************

    SpalartAllmaras2D::SpalartAllmaras2D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    SpalartAllmaras2D::SpalartAllmaras2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {

    }

    Element::Pointer SpalartAllmaras2D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new SpalartAllmaras2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    SpalartAllmaras2D::~SpalartAllmaras2D()
    {
    }

    //************************************************************************************
    //************************************************************************************

    void SpalartAllmaras2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const double sigma = 2.0 / 3.0;

        const unsigned int number_of_points = GetGeometry().size();
        const double lumping_factor = 1.00 / double(number_of_points);
        const unsigned int TDim = 2;

        if (rLeftHandSideMatrix.size1() != number_of_points)
            rLeftHandSideMatrix.resize(number_of_points, number_of_points, false);

        if (rRightHandSideVector.size() != number_of_points)
            rRightHandSideVector.resize(number_of_points, false);

        //mThisIntegrationMethod= GeometryData::GI_GAUSS_1;

        boost::numeric::ublas::bounded_matrix<double, 3, 3 > MassFactors = 1.0 / 3.0 * IdentityMatrix(3, 3);
        boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
        array_1d<double, 3 > N;
        array_1d<double, 2 > vel_gauss;
        array_1d<double, 3 > temp_vec_np;
        array_1d<double, 3 > u_DN;
        array_1d<double, 2 > grad_g;
        boost::numeric::ublas::bounded_matrix<double, 2, 2 > Identity = IdentityMatrix(2, 2);
        boost::numeric::ublas::bounded_matrix<double, 2, 2 > First;
        boost::numeric::ublas::bounded_matrix<double, 2, 2 > Second;
        boost::numeric::ublas::bounded_matrix<double, 2, 3 > Third;


        //getting data for the given geometry
        double Area;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);


        //calculating viscosity
        double molecular_viscosity = GetGeometry()[0].FastGetSolutionStepValue(MOLECULAR_VISCOSITY);
        double turbulent_viscosity = GetGeometry()[0].FastGetSolutionStepValue(TURBULENT_VISCOSITY);
        double proj = GetGeometry()[0].FastGetSolutionStepValue(TEMP_CONV_PROJ);
        const array_1d<double, 3 > & v = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
        const array_1d<double, 3 > & w = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);


        for (unsigned int j = 0; j < TDim; j++)
            vel_gauss[j] = v[j] - w[j];

        for (unsigned int i = 1; i < number_of_points; i++)
        {
            molecular_viscosity += GetGeometry()[i].FastGetSolutionStepValue(MOLECULAR_VISCOSITY);
            turbulent_viscosity += GetGeometry()[i].FastGetSolutionStepValue(TURBULENT_VISCOSITY);
            proj += GetGeometry()[i].FastGetSolutionStepValue(TEMP_CONV_PROJ);

            const array_1d<double, 3 > & v = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3 > & w = GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY);
            for (unsigned int j = 0; j < TDim; j++)
                vel_gauss[j] += v[j] - w[j];

        }
        molecular_viscosity *= lumping_factor;
        turbulent_viscosity *= lumping_factor;
        const double conductivity = (molecular_viscosity + turbulent_viscosity) / sigma;
        proj *= lumping_factor;
        vel_gauss *= lumping_factor;
        const double c1 = 4.00;
        const double c2 = 2.00;
        const double h = sqrt(2.00 * Area);
        const double norm_u = norm_2(vel_gauss);
        double tau1 = (h * h) / (c1 * conductivity + c2 * norm_u * h);

        // 		double g=0.0;
        double p1 = DN_DX(0, 0) * GetGeometry()[0].FastGetSolutionStepValue(TURBULENT_VISCOSITY) + DN_DX(1, 0) * GetGeometry()[1].FastGetSolutionStepValue(TURBULENT_VISCOSITY) + DN_DX(2, 0) * GetGeometry()[2].FastGetSolutionStepValue(TURBULENT_VISCOSITY);
        double p2 = DN_DX(0, 1) * GetGeometry()[0].FastGetSolutionStepValue(TURBULENT_VISCOSITY) + DN_DX(1, 1) * GetGeometry()[1].FastGetSolutionStepValue(TURBULENT_VISCOSITY) + DN_DX(2, 1) * GetGeometry()[2].FastGetSolutionStepValue(TURBULENT_VISCOSITY);
        grad_g[0] = p1;
        grad_g[1] = p2;
        //     		double norm_g =norm_2(grad_g);

        const double source_term = this->CalculateSourceTerm(DN_DX, N, molecular_viscosity, turbulent_viscosity);

        double res = (inner_prod(vel_gauss, grad_g));
        double norm_grad = norm_2(grad_g);
        double k_aux = fabs(res) / (norm_grad + 0.000000000001);


        noalias(First) = outer_prod(vel_gauss, trans(vel_gauss));
        First /= ((norm_u + 0.0000000001)*(norm_u + 0.0000000001));
        noalias(Second) = Identity - First;
        noalias(Third) = prod(Second, trans(DN_DX));


        //calculating parameter tau


        //getting the BDF2 coefficients (not fixed to allow variable time step)
        //the coefficients INCLUDE the time step
        const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

        //CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
        noalias(u_DN) = prod(DN_DX, vel_gauss);
        noalias(rLeftHandSideMatrix) = outer_prod(N, u_DN);

        //CONVECTION STABILIZING CONTRIBUTION (Suu)
        noalias(rLeftHandSideMatrix) += tau1 * outer_prod(u_DN, u_DN);

        //VISCOUS CONTRIBUTION TO THE STIFFNESS MATRIX
        noalias(rLeftHandSideMatrix) += (conductivity * prod(DN_DX, trans(DN_DX)) + k_aux * h * prod(DN_DX, Third));

        //filling the mass factors
        MassFactors(0, 0) = 1.00 / 3.00;
        MassFactors(0, 1) = 0.00;
        MassFactors(0, 2) = 0.00;
        MassFactors(1, 0) = 0.00;
        MassFactors(1, 1) = 1.00 / 3.00;
        MassFactors(1, 2) = 0.00;
        MassFactors(2, 0) = 0.00;
        MassFactors(2, 1) = 0.00;
        MassFactors(2, 2) = 1.00 / 3.00;

        //INERTIA CONTRIBUTION
        noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * MassFactors;

        // RHS = Fext
        noalias(rRightHandSideVector) = source_term*N;

        //RHS += Suy * proj[component]
        noalias(rRightHandSideVector) += (tau1 * proj) * u_DN;

        //adding the inertia ter
        // RHS += M*vhistory
        //calculating the historical velocity
        for (unsigned int iii = 0; iii < number_of_points; iii++)
            temp_vec_np[iii] = BDFcoeffs[1] * GetGeometry()[iii].FastGetSolutionStepValue(TURBULENT_VISCOSITY, 1);
        for (unsigned int step = 2; step < BDFcoeffs.size(); step++)
        {
            for (unsigned int iii = 0; iii < number_of_points; iii++)
                temp_vec_np[iii] += BDFcoeffs[step] * GetGeometry()[iii].FastGetSolutionStepValue(TURBULENT_VISCOSITY, step);
        }
        noalias(rRightHandSideVector) -= prod(MassFactors, temp_vec_np);

        //subtracting the dirichlet term
        // RHS -= LHS*temperatures
        for (unsigned int iii = 0; iii < number_of_points; iii++)
            temp_vec_np[iii] = GetGeometry()[iii].FastGetSolutionStepValue(TURBULENT_VISCOSITY);

        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, temp_vec_np);

        rRightHandSideVector *= Area;
        rLeftHandSideMatrix *= Area;

        KRATOS_CATCH("");
    }



    //************************************************************************************
    //************************************************************************************

    void SpalartAllmaras2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR(std::logic_error, "method not implemented", "");
    }

    //************************************************************************************
    //************************************************************************************
    // this subroutine calculates the nodal contributions for the explicit steps of the
    // fractional step procedure

    void SpalartAllmaras2D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY
        int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

        boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
        array_1d<double, 3 > N;
        array_1d<double, 2 > vel_gauss;
        array_1d<double, 3 > temp_vec_np;
        array_1d<double, 3 > u_DN;

        //getting data for the given geometry
        double Area;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);


        if (FractionalStepNumber == 2) //calculation of temperature convective projection
        {
            const unsigned int number_of_points = GetGeometry().size();
            const double lumping_factor = 1.00 / double(number_of_points);
            unsigned int TDim = 2;

            //calculating viscosity
            temp_vec_np[0] = GetGeometry()[0].FastGetSolutionStepValue(TURBULENT_VISCOSITY);
            const array_1d<double, 3 > & v = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3 > & w = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
            for (unsigned int j = 0; j < TDim; j++)
                vel_gauss[j] = v[j] - w[j];

            for (unsigned int i = 1; i < number_of_points; i++)
            {
                temp_vec_np[i] = GetGeometry()[i].FastGetSolutionStepValue(TURBULENT_VISCOSITY);
                const array_1d<double, 3 > & v = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
                const array_1d<double, 3 > & w = GetGeometry()[i].FastGetSolutionStepValue(MESH_VELOCITY);
                for (unsigned int j = 0; j < TDim; j++)
                    vel_gauss[j] += v[j] - w[j];

            }
            vel_gauss *= lumping_factor;

            //calculating convective auxiliary vector
            noalias(u_DN) = prod(DN_DX, vel_gauss);
            double temp_conv = inner_prod(u_DN, temp_vec_np);
            temp_conv *= Area;

            for (unsigned int i = 0; i < number_of_points; i++)
            {
                GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += lumping_factor*Area;
                GetGeometry()[i].FastGetSolutionStepValue(TEMP_CONV_PROJ) += lumping_factor*temp_conv;
                ;
            }
        }
        KRATOS_CATCH("");
    }


    //************************************************************************************
    //************************************************************************************

    void SpalartAllmaras2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
    {
        unsigned int number_of_nodes = GetGeometry().PointsNumber();
        if (rResult.size() != number_of_nodes)
            rResult.resize(number_of_nodes, false);

        for (unsigned int i = 0; i < number_of_nodes; i++)
            rResult[i] = GetGeometry()[i].GetDof(TURBULENT_VISCOSITY).EquationId();
    }

    //************************************************************************************
    //************************************************************************************

    void SpalartAllmaras2D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
    {
        unsigned int number_of_nodes = GetGeometry().PointsNumber();

        if (ElementalDofList.size() != number_of_nodes)
            ElementalDofList.resize(number_of_nodes);

        for (unsigned int i = 0; i < number_of_nodes; i++)
            ElementalDofList[i] = GetGeometry()[i].pGetDof(TURBULENT_VISCOSITY);

    }

    double SpalartAllmaras2D::CalculateSourceTerm(const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX,
            const array_1d<double,3>& N,
            const double molecular_viscosity,
            const double turbulent_viscosity)
    {
        const double sigma = 2.0 / 3.0;
        const double cb1 = 0.1355;
        const double cb2 = 0.622;
        const double kappa = 0.41;
        const double cw1 = cb1 / (kappa * kappa) + (1.0 + cb2) / sigma;
        const double cw2 = 0.3;
        const double cw3 = 2.0;
        const double cv1 = 7.1;
//        const double ct1 = 1.0; // For ft1 (trip term, not implemented)
//        const double ct2 = 2.0;
        const double ct3 = 1.2; // Original value is 1.1, correction by Spalart
        const double ct4 = 0.5; // Original value is 2.0, correction by Spalart

        const double xi = turbulent_viscosity/molecular_viscosity;
        const double fv1 = xi * xi * xi / (xi * xi * xi + cv1 * cv1 * cv1);
        const double fv2 = 1.0 - xi / ( 1.0 + xi * fv1);
        const double S = this->AntimetricGradientNorm(DN_DX);

        double distance = N[0]*this->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE) +
                          N[1]*this->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE) +
                          N[2]*this->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE);

        const double S_hat = S + fv2 * turbulent_viscosity / (distance * distance * kappa * kappa);

        const double ft2 = ct3 * exp(-ct4*xi*xi);

        const double r = turbulent_viscosity / ( S_hat * kappa * kappa * distance * distance);
        const double g = r + cw2 * (pow(r,6) - r);
        const double fw = g * pow( (1.0+pow(cw3,6)) / ( pow(g,6) + pow(cw3,6) ) , 1.0/6.0);

        double p1 = DN_DX(0, 0) * GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY) + DN_DX(1, 0) * GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY) + DN_DX(2, 0) * GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
        double p2 = DN_DX(0, 1) * GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY) + DN_DX(1, 1) * GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY) + DN_DX(2, 1) * GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
        const double norm2_gradnu = p1*p1 + p2*p2;

        double source_term = cb1 * (1.0 - ft2) * S_hat * turbulent_viscosity;
        source_term += cb2 * norm2_gradnu / sigma;

        source_term -= (cw1 * fw - ft2 * cb1/(kappa*kappa) ) * pow(molecular_viscosity/distance,2);

        return source_term;
    }

    double SpalartAllmaras2D::AntimetricGradientNorm(const boost::numeric::ublas::bounded_matrix<double,3,2>& rShapeDeriv)
        {
            const unsigned int GradientSize = 3; // Number of different terms in the symmetric gradient matrix
            array_1d<double,GradientSize> GradientVector( GradientSize, 0.0 );
            unsigned int Index;

            // Compute Symmetric Grad(u). Note that only the lower half of the matrix is calculated
            for (unsigned int k = 0; k < 3; ++k)
            {
                const array_1d<double, 3> & rNodeVel = this->GetGeometry()[k].FastGetSolutionStepValue(VELOCITY);
                Index = 0;
                for (unsigned int i = 0; i < 2; ++i)
                {
                    for (unsigned int j = 0; j < i; ++j) // Off-diagonal
                        GradientVector[Index++] += 0.5 * (rShapeDeriv(k, j) * rNodeVel[i] - rShapeDeriv(k, i) * rNodeVel[j]);
                    GradientVector[Index++] = 0.0; // Diagonal
                }
            }

            // Norm[ Antimetric Grad(u) ] = ( 2 * Omegaij * Omegaij )^(1/2)
            Index = 0;
            double NormS(0.0);
            for (unsigned int i = 0; i < 2; ++i)
            {
                for (unsigned int j = 0; j < i; ++j)
                {
                    NormS += 2.0 * GradientVector[Index] * GradientVector[Index]; // Using symmetry, lower half terms of the matrix are added twice
                    ++Index;
                }
                NormS += GradientVector[Index] * GradientVector[Index]; // Diagonal terms
                ++Index; // Diagonal terms
            }

            NormS = sqrt( 2.0 * NormS );
            return NormS;
        }

} // Namespace Kratos


