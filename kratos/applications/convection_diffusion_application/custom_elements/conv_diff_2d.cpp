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
#include "custom_elements/conv_diff_2d.h"
#include "convection_diffusion_application.h"
#include "includes/convection_diffusion_settings.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{


    //************************************************************************************
    //************************************************************************************

    ConvDiff2D::ConvDiff2D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    ConvDiff2D::ConvDiff2D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {

    }

    Element::Pointer ConvDiff2D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new ConvDiff2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    ConvDiff2D::~ConvDiff2D()
    {
    }

    //************************************************************************************
    //************************************************************************************

    void ConvDiff2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY

                //getting the BDF2 coefficients (not fixed to allow variable time step)
                //the coefficients INCLUDE the time step
                const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];

        const unsigned int number_of_points = GetGeometry().size();
        const double lumping_factor = 1.00 / double(number_of_points);
        unsigned int TDim = 2;

        if (rLeftHandSideMatrix.size1() != number_of_points)
            rLeftHandSideMatrix.resize(number_of_points, number_of_points, false);

        if (rRightHandSideVector.size() != number_of_points)
            rRightHandSideVector.resize(number_of_points, false);

        //mThisIntegrationMethod= GeometryData::GI_GAUSS_1;

        boost::numeric::ublas::bounded_matrix<double, 3, 3 > msMassFactors = 1.0 / 3.0 * IdentityMatrix(3, 3);
        boost::numeric::ublas::bounded_matrix<double, 3, 2 > msDN_DX;
        array_1d<double, 3 > msN;
        array_1d<double, 2 > ms_vel_gauss;
        array_1d<double, 3 > ms_temp_vec_np;
        array_1d<double, 3 > ms_u_DN;
        array_1d<double, 2 > grad_g;
        boost::numeric::ublas::bounded_matrix<double, 2, 2 > Identity = IdentityMatrix(2, 2);
        boost::numeric::ublas::bounded_matrix<double, 2, 2 > First;
        boost::numeric::ublas::bounded_matrix<double, 2, 2 > Second;
        boost::numeric::ublas::bounded_matrix<double, 2, 3 > Third;


        //getting data for the given geometry
        double Area;
        GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);

        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);


        //calculating viscosity

        const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
        //const Variable<array_1d<double,3> >& rConvectionVar = my_settings->GetConvectionVariable();
        const Variable<double>& rDiffusionVar = my_settings->GetDiffusionVariable();
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
        const Variable<double>& rSourceVar = my_settings->GetVolumeSourceVariable();
        //const Variable<double>& rSurfaceSourceVar =my_settings->GetSurfaceSourceVariable();
        const Variable<array_1d<double, 3 > >& rMeshVelocityVar = my_settings->GetMeshVelocityVariable();
	const Variable<double>& rProjectionVariable = my_settings->GetProjectionVariable();

        double conductivity = GetGeometry()[0].FastGetSolutionStepValue(rDiffusionVar);
        double specific_heat = GetGeometry()[0].FastGetSolutionStepValue(SPECIFIC_HEAT);
        double density = GetGeometry()[0].FastGetSolutionStepValue(rDensityVar);
        double heat_flux = GetGeometry()[0].FastGetSolutionStepValue(rSourceVar);
        double proj = GetGeometry()[0].FastGetSolutionStepValue(rProjectionVariable);
        const array_1d<double, 3 > & v = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY); //VELOCITY
        const array_1d<double, 3 > & w = GetGeometry()[0].FastGetSolutionStepValue(rMeshVelocityVar); //
	double nu = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);	

        for (unsigned int j = 0; j < TDim; j++)
            ms_vel_gauss[j] = v[j] - w[j];

        for (unsigned int i = 1; i < number_of_points; i++)
        {
            conductivity += GetGeometry()[i].FastGetSolutionStepValue(rDiffusionVar);
            density += GetGeometry()[i].FastGetSolutionStepValue(rDensityVar);
            specific_heat += GetGeometry()[i].FastGetSolutionStepValue(SPECIFIC_HEAT);
            heat_flux += GetGeometry()[i].FastGetSolutionStepValue(rSourceVar);
            //proj += GetGeometry()[i].FastGetSolutionStepValue(TEMP_CONV_PROJ);
	    proj += GetGeometry()[i].FastGetSolutionStepValue(rProjectionVariable);
	    nu += GetGeometry()[i].FastGetSolutionStepValue(VISCOSITY);
            const array_1d<double, 3 > & v = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3 > & w = GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar);
            for (unsigned int j = 0; j < TDim; j++)
                ms_vel_gauss[j] += v[j] - w[j];

        }
        conductivity *= lumping_factor;
        density *= lumping_factor;
        specific_heat *= lumping_factor;
        heat_flux *= lumping_factor;
        proj *= lumping_factor;
        ms_vel_gauss *= lumping_factor;
	nu *= lumping_factor;



        double c1 = 4.00;
        double c2 = 2.00;
        double h = sqrt(2.00 * Area);
        double norm_u = norm_2(ms_vel_gauss);
        double tau1 = (h * h) / (density * specific_heat * BDFcoeffs[0] * h * h + c1 * conductivity + c2 * density * specific_heat * (norm_u + 1e-6) * h);
        //        tau1 *= 0.1;
        // 		double alpha = conductivity/(density*specific_heat);

        // 		double g=0.0;
        double p1 = msDN_DX(0, 0) * GetGeometry()[0].FastGetSolutionStepValue(rUnknownVar) + msDN_DX(1, 0) * GetGeometry()[1].FastGetSolutionStepValue(rUnknownVar) + msDN_DX(2, 0) * GetGeometry()[2].FastGetSolutionStepValue(rUnknownVar);
        double p2 = msDN_DX(0, 1) * GetGeometry()[0].FastGetSolutionStepValue(rUnknownVar) + msDN_DX(1, 1) * GetGeometry()[1].FastGetSolutionStepValue(rUnknownVar) + msDN_DX(2, 1) * GetGeometry()[2].FastGetSolutionStepValue(rUnknownVar);
        grad_g[0] = p1;
        grad_g[1] = p2;
        //     		double norm_g =norm_2(grad_g);

//        double res = (inner_prod(ms_vel_gauss, grad_g)); //+ 0.333333333333333 * (t0media+t1media+t2media)*(1/dt)*density*conductivity;
//        double aux_res = fabs(res - proj);
//        if (fabs(res) > aux_res)
//            res = aux_res;
        //        res -= proj;
//        res *= density*specific_heat;
//        double norm_grad = norm_2(grad_g);
//        double k_aux = fabs(res) / (norm_grad + 1e-6);
//        k_aux *= 0.707;

	double Schmidt_Prandl=1.0;
	double nu_turbulent = 0.0;
	const double Cs = this->GetValue(C_SMAGORINSKY);
	if (Cs != 0.0){
        	nu_turbulent = ComputeSmagorinskyViscosity(msDN_DX, h, Cs, nu);
	}


	double norm_g =norm_2(grad_g);

	double res = density * specific_heat*(inner_prod(ms_vel_gauss,grad_g)) ;//+ 0.333333333333333 * (t0media+t1media+t2media)*(1/dt)*density*conductivity;
	double norm_grad=norm_2(grad_g);
	double k_aux=fabs(res) /(norm_grad + 0.000000000001);


        noalias(First) = outer_prod(ms_vel_gauss, trans(ms_vel_gauss));
        First /= ((norm_u + 1e-6)*(norm_u + 1e-6));
        noalias(Second) = Identity - First;
        noalias(Third) = prod(Second, trans(msDN_DX));


        //calculating parameter tau



        //CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
        noalias(ms_u_DN) = prod(msDN_DX, ms_vel_gauss);
        noalias(rLeftHandSideMatrix) = (density * specific_heat) * outer_prod(msN, ms_u_DN);

        //CONVECTION STABILIZING CONTRIBUTION (Suu)
        noalias(rLeftHandSideMatrix) += density * specific_heat * density * specific_heat * tau1 * outer_prod(ms_u_DN, ms_u_DN);

        //VISCOUS CONTRIBUTION TO THE STIFFNESS MATRIX
        noalias(rLeftHandSideMatrix) += (conductivity * prod(msDN_DX, trans(msDN_DX)) + k_aux * h * prod(msDN_DX, Third));

	//noalias(rLeftHandSideMatrix) += density * (nu_turbulent /Schmidt_Prandl ) * prod(msDN_DX,trans(msDN_DX));
        //filling the mass factors
        //        msMassFactors(0, 0) = 1/6.0; 1.00 / 3.00;
        //        msMassFactors(0, 1) = 0.00;
        //        msMassFactors(0, 2) = 0.00;
        //        msMassFactors(1, 0) = 0.00;
        //        msMassFactors(1, 1) = 1.00 / 3.00;
        //        msMassFactors(1, 2) = 0.00;
        //        msMassFactors(2, 0) = 0.00;
        //        msMassFactors(2, 1) = 0.00;
        //        msMassFactors(2, 2) = 1.00 / 3.00;

        //INERTIA CONTRIBUTION
        noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * (density * specific_heat) * msMassFactors;

        // RHS = Fext (heat per unit mass)
        noalias(rRightHandSideVector) = (heat_flux * density) * msN;

        //RHS += Suy * proj[component]
        noalias(rRightHandSideVector) += density * specific_heat * density * specific_heat * (tau1 * proj) * ms_u_DN;

        //adding the inertia terms
        // RHS += M*vhistory
        //calculating the historical velocity
        for (unsigned int iii = 0; iii < number_of_points; iii++)
            ms_temp_vec_np[iii] = BDFcoeffs[1] * GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar, 1);
        for (unsigned int step = 2; step < BDFcoeffs.size(); step++)
        {
            for (unsigned int iii = 0; iii < number_of_points; iii++)
                ms_temp_vec_np[iii] += BDFcoeffs[step] * GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar, step);
        }
        noalias(rRightHandSideVector) -= prod(msMassFactors, ms_temp_vec_np * density * specific_heat);

        //subtracting the dirichlet term
        // RHS -= LHS*temperatures
        for (unsigned int iii = 0; iii < number_of_points; iii++)
            ms_temp_vec_np[iii] = GetGeometry()[iii].FastGetSolutionStepValue(rUnknownVar);
        noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, ms_temp_vec_np);


        rRightHandSideVector *= Area;
        rLeftHandSideMatrix *= Area;

        KRATOS_CATCH("");
    }



    //************************************************************************************
    //************************************************************************************

    void ConvDiff2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_ERROR(std::logic_error, "method not implemented", "");
    }

    //************************************************************************************
    //************************************************************************************
    // this subroutine calculates the nodal contributions for the explicit steps of the
    // fractional step procedure

    void ConvDiff2D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
        KRATOS_TRY
                int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

        boost::numeric::ublas::bounded_matrix<double, 3, 2 > msDN_DX;
        array_1d<double, 3 > msN;
        array_1d<double, 2 > ms_vel_gauss;
        array_1d<double, 3 > ms_temp_vec_np;
        array_1d<double, 3 > ms_u_DN;

        //getting data for the given geometry
        double Area;
        GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);
        ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
        //const Variable<array_1d<double,3> >& rConvectionVar = my_settings->GetConvectionVariable();
        const Variable<array_1d<double, 3 > >& rMeshVelocityVar = my_settings->GetMeshVelocityVariable();
	const Variable<double>& rProjectionVariable = my_settings->GetProjectionVariable();

        if (FractionalStepNumber == 2) //calculation of temperature convective projection
        {
            const unsigned int number_of_points = GetGeometry().size();
            const double lumping_factor = 1.00 / double(number_of_points);
            unsigned int TDim = 2;

            //calculating viscosity
            ms_temp_vec_np[0] = GetGeometry()[0].FastGetSolutionStepValue(rUnknownVar);
            const array_1d<double, 3 > & v = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3 > & w = GetGeometry()[0].FastGetSolutionStepValue(rMeshVelocityVar);
            for (unsigned int j = 0; j < TDim; j++)
                ms_vel_gauss[j] = v[j] - w[j];

            for (unsigned int i = 1; i < number_of_points; i++)
            {
                ms_temp_vec_np[i] = GetGeometry()[i].FastGetSolutionStepValue(rUnknownVar);
                const array_1d<double, 3 > & v = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
                const array_1d<double, 3 > & w = GetGeometry()[i].FastGetSolutionStepValue(rMeshVelocityVar);
                for (unsigned int j = 0; j < TDim; j++)
                    ms_vel_gauss[j] += v[j] - w[j];

            }
            ms_vel_gauss *= lumping_factor;

            //calculating convective auxiliary vector
            noalias(ms_u_DN) = prod(msDN_DX, ms_vel_gauss);
            double temp_conv = inner_prod(ms_u_DN, ms_temp_vec_np);
            temp_conv *= Area;

            for (unsigned int i = 0; i < number_of_points; i++)
            {
                GetGeometry()[i].FastGetSolutionStepValue(NODAL_AREA) += lumping_factor*Area;
                //GetGeometry()[i].FastGetSolutionStepValue(TEMP_CONV_PROJ) += lumping_factor*temp_conv;
		GetGeometry()[i].FastGetSolutionStepValue(rProjectionVariable) += lumping_factor*temp_conv;
            }
        }
        KRATOS_CATCH("");
    }


    //************************************************************************************
    //************************************************************************************

    void ConvDiff2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
    {
        ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();
        unsigned int number_of_nodes = GetGeometry().PointsNumber();
        if (rResult.size() != number_of_nodes)
            rResult.resize(number_of_nodes, false);

        for (unsigned int i = 0; i < number_of_nodes; i++)
            rResult[i] = GetGeometry()[i].GetDof(rUnknownVar).EquationId();
    }

    //************************************************************************************
    //************************************************************************************

    void ConvDiff2D::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
    {
        unsigned int number_of_nodes = GetGeometry().PointsNumber();
        ConvectionDiffusionSettings::Pointer my_settings = CurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);
        const Variable<double>& rUnknownVar = my_settings->GetUnknownVariable();

        if (ElementalDofList.size() != number_of_nodes)
            ElementalDofList.resize(number_of_nodes);

        for (unsigned int i = 0; i < number_of_nodes; i++)
            ElementalDofList[i] = GetGeometry()[i].pGetDof(rUnknownVar);

    }
 
    //************************************************************************************
    //************************************************************************************
	double ConvDiff2D::ComputeSmagorinskyViscosity(const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const double& h, const double& C, const double nu )
  {
    boost::numeric::ublas::bounded_matrix<double, 2, 2 > dv_dx = ZeroMatrix(2, 2);
    
    // Compute Symmetric Grad(u). Note that only the lower half of the matrix is filled
    for (unsigned int k = 0; k < 3; ++k)
        {
	  const array_1d< double, 3 > & rNodeVel = this->GetGeometry()[k].FastGetSolutionStepValue(VELOCITY);
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
    
    return 2.0 * C * C * h * h * NormS;
  }


} // Namespace Kratos


