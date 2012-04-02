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
//   Last modified by:    $Author: kazem $
//   Date:                $Date: 2009-01-21 14:15:02 $
//   Revision:            $Revision: 1.6 $
//
//

//#define GRADPN_FORM
//#define STOKES

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/asgs_compressible_2d.h"
#include "utilities/math_utils.h"
#include "incompressible_fluid_application.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
ASGSCompressible2D::ASGSCompressible2D(IndexType NewId, GeometryType::Pointer pGeometry)
        : ASGS2D(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!

}

//************************************************************************************
//************************************************************************************
ASGSCompressible2D::ASGSCompressible2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
        : ASGS2D(NewId, pGeometry, pProperties)
{

}

Element::Pointer ASGSCompressible2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{

    KRATOS_TRY
    return Element::Pointer(new ASGSCompressible2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
    KRATOS_CATCH("");
}

ASGSCompressible2D::~ASGSCompressible2D()
{
}

//************************************************************************************
//************************************************************************************

void ASGSCompressible2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    int nodes_number = 3;
    int dim = 2;
    unsigned int matsize = nodes_number * (dim + 1);

    if (rLeftHandSideMatrix.size1() != matsize)
        rLeftHandSideMatrix.resize(matsize, matsize,false); //false says not to preserve existing storage!!

    if (rRightHandSideVector.size() != matsize)
        rRightHandSideVector.resize(matsize,false); //false says not to preserve existing storage!!


    noalias(rLeftHandSideMatrix) = ZeroMatrix(matsize, matsize);
    noalias(rRightHandSideVector) = ZeroVector(matsize);

    double delta_t = rCurrentProcessInfo[DELTA_TIME];



    //getting data for the given geometry
    double Area;
    boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX = ZeroMatrix(3,2);
    array_1d<double,3> N = ZeroVector(3);
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);

    double tauone;
    double tautwo;
    CalculateTau(N,tauone, tautwo, delta_t, Area, rCurrentProcessInfo);


    //add body force and momentum
    AddBodyForceAndMomentum(rRightHandSideVector, DN_DX, N, delta_t, Area, tauone, tautwo);


    KRATOS_CATCH("")
}
//************************************************************************************
//************************************************************************************
void ASGSCompressible2D::CalculateLocalVelocityContribution(MatrixType& rDampMatrix,VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    int nodes_number = 3;
    int dim = 2;
    unsigned int matsize = nodes_number*(dim+1);

    if (rDampMatrix.size1() != matsize)
        rDampMatrix.resize(matsize,matsize,false); //false says not to preserve existing storage!!


    noalias(rDampMatrix) = ZeroMatrix(matsize,matsize);

    double delta_t= rCurrentProcessInfo[DELTA_TIME];



    //getting data for the given geometry
    double Area;
    boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX = ZeroMatrix(3,2);
    array_1d<double,3> N = ZeroVector(3);
    GeometryUtils::CalculateGeometryData(GetGeometry(),DN_DX,N,Area);


    //viscous term
    CalculateViscousTerm(rDampMatrix, DN_DX, Area);

    // DCoperator
    CalcualteDCOperator(rDampMatrix, DN_DX, Area);

    //Advective term
    double tauone;
    double tautwo;
    CalculateTau(N,tauone, tautwo, delta_t, Area, rCurrentProcessInfo);

    double VC2;
    CalculateSoundVelocity(GetGeometry(), VC2);

    CalculateAdvectiveTerm(rDampMatrix, DN_DX,N, tauone, tautwo, delta_t, Area);
    //calculate pressure term
    CalculatePressureTerm(rDampMatrix, DN_DX, N, delta_t,Area);

    //compute projections
    //stabilization terms
    //KRATOS_WATCH("ñññññññññ calculate stabilizing terms ñññññññññ");
    double tau_div = tautwo*VC2;
    CalculateDivStblTerm(rDampMatrix, DN_DX, tau_div, Area);
    //CalculateAdvStblAllTerms(rDampMatrix,rRightHandSideVector, DN_DX, N, tauone,delta_t, Area);
    CalculateGradStblAllTerms(rDampMatrix,rRightHandSideVector,DN_DX,N, delta_t, tauone, Area);

    CalculateDivPdotStblTerms(rDampMatrix,rRightHandSideVector, DN_DX,N, delta_t, tautwo, Area);
    CalculateNonlinearStblTerm(rRightHandSideVector, DN_DX,N,  delta_t, tautwo, Area);

    //KRATOS_WATCH(rRightHandSideVector);


    //add body force and momentum
//         AddBodyForceAndMomentum(rRightHandSideVector, N, delta_t, Area, tauone, tautwo);

    CalculateResidual(rDampMatrix, rRightHandSideVector);

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void ASGSCompressible2D::CalculateMassContribution(MatrixType& K,const double time,const double area)
{
    KRATOS_TRY
    double lump_mass_fac = area * 0.333333333333333333333333;
    double density;
    double mu;
    calculatedensity(GetGeometry(), density, mu);

    //update density and calculate sound velocity
    //double VC2;
    //CalculateSoundVelocity(GetGeometry(), VC2);

    int nodes_number = 3;
    int dof = 2;
    for ( int nd = 0; nd< nodes_number; nd++)
    {
        int row = nd*(dof + 1);
        for ( int jj=0; jj< dof; jj++)
            K(row + jj, row + jj) += density*lump_mass_fac;

        //add pressure mass
        K(row + dof , row + dof ) += lump_mass_fac;
    }


    KRATOS_CATCH("")

}
//************************************************************************************
//************************************************************************************
void ASGSCompressible2D::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //lumped
    unsigned int dimension = GetGeometry().WorkingSpaceDimension();
    unsigned int NumberOfNodes = GetGeometry().size();
    unsigned int MatSize = (dimension + 1) * NumberOfNodes;
    if (rMassMatrix.size1() != MatSize)
        rMassMatrix.resize(MatSize,MatSize,false);

    rMassMatrix = ZeroMatrix(MatSize,MatSize);
    double delta_t= rCurrentProcessInfo[DELTA_TIME];


    //getting data for the given geometry
    double Area;
    boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX = ZeroMatrix(3,2);
    array_1d<double,3> N = ZeroVector(3);
    GeometryUtils::CalculateGeometryData(GetGeometry(),DN_DX,N,Area);

    //Calculate tau
    double tauone;
    double tautwo;
    CalculateTau(N,tauone, tautwo, delta_t, Area, rCurrentProcessInfo);

    CalculateMassContribution(rMassMatrix,delta_t,Area);
    //add stablilization terms due to advective term (a)grad(V) * ro*Acce
    //CalculateAdvMassStblTerms(rMassMatrix, DN_DX, N,tauone,Area);
    //add stablilization terms due to grad term grad(q) * ro*Acce
    //CalculateGradMassStblTerms(rMassMatrix, DN_DX, tauone,Area);
    //add compressible stabilization terms
    //CalculateCompressibleStblTerms(rMassMatrix, DN_DX,N, tautwo,Area);

    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void ASGSCompressible2D::CalculateGradMassStblTerms(MatrixType& M, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const array_1d<double,3>& N, const double tauone, const double area)
{
    KRATOS_TRY
    int nodes_number = 3;
    int dof = 2;

    //   double density;
    // double mu;
    // calculatedensity(GetGeometry(), density, mu);

    double VC2;
    CalculateSoundVelocity(GetGeometry(), VC2);

    //build 1*tau1*ro Nacc grad q)
    double fac = area*tauone*VC2;
    for (int ii = 0; ii < nodes_number; ii++)
    {
        int row = ii * (dof + 1);
        for (int jj = 0; jj < nodes_number; jj++)
        {
            int column = jj * (dof + 1) + dof;

            //K(row,column) += -1*area * fac* N(ii) * DN_DX(jj,0);
            M(column, row) +=  fac * N[ii] * DN_DX(jj, 0);

            //K(row + 1,column) += -1*area * fac* N(ii) * DN_DX(jj,1);
            M(column, row + 1) +=  fac * N[ii] * DN_DX(jj, 1);
        }
    }

    KRATOS_CATCH("")

}

//*************************************************************************************
//*************************************************************************************


void ASGSCompressible2D::CalculateGradStblAllTerms(MatrixType& K, VectorType& F, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX,const array_1d<double,3>& N, const double time, const double tauone, const double area)
{
    KRATOS_TRY
    int nodes_number = 3;
    int dof = 2;


    double VC2;
    CalculateSoundVelocity(GetGeometry(), VC2);

    double density;
    double mu;
    calculatedensity(GetGeometry(), density, mu);

    //build 1*(grad q . grad p) stabilization term & assemble
    boost::numeric::ublas::bounded_matrix<double, 3, 3 > gard_opr = ZeroMatrix(nodes_number, nodes_number);
    gard_opr =  prod(DN_DX, trans(DN_DX));

    for (int ii = 0; ii < nodes_number; ii++)
    {
        int row = ii * (dof + 1) + dof;

        for (int jj = 0; jj < nodes_number; jj++)
        {
            int column = jj * (dof + 1) + dof;

            K(row, column) += area * tauone * gard_opr(ii, jj);

        }
    }

    //build 1*(grad q) (Fbody ) stabilization term & assemble
    array_1d<double, 2 > bdf = ZeroVector(2);
    const array_1d<double, 3 > bdf0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);
    const array_1d<double, 3 > bdf1 = GetGeometry()[1].FastGetSolutionStepValue(BODY_FORCE);
    const array_1d<double, 3 > bdf2 = GetGeometry()[2].FastGetSolutionStepValue(BODY_FORCE);

    //add acceleration
    const array_1d<double,3>& acce0 = GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION);
    const array_1d<double,3>& acce1 = GetGeometry()[1].FastGetSolutionStepValue(ACCELERATION);
    const array_1d<double,3>& acce2 = GetGeometry()[2].FastGetSolutionStepValue(ACCELERATION);


    bdf[0] = N[0]*(density * (bdf0[0] - acce0[0])) + N[1]*(density * (bdf1[0] - acce1[0])) + N[2]*(density * (bdf2[0] - acce2[0]));
    bdf[1] = N[0]*(density * (bdf0[1] - acce0[1])) + N[1]*(density * (bdf1[1] - acce1[1])) + N[2]*(density * (bdf2[1] - acce2[1]));


    array_1d<double, 3 > fbd_stblterm = ZeroVector(nodes_number);
    fbd_stblterm = tauone * prod(DN_DX, bdf);


    for (int ii = 0; ii < nodes_number; ++ii)
    {
        int index = ii * (dof + 1) + dof;
        F[index] +=  area * fbd_stblterm[ii];
    }

    KRATOS_CATCH("")
}
//*************************************************************************************
//*************************************************************************************
void ASGSCompressible2D::CalculateDivPdotStblTerms(MatrixType& K,VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,const array_1d<double,3>& N, const double time,const double tautwo,const double area)
{
    KRATOS_TRY
    //tau*div(V).P_dot
    double stbl_fac = tautwo * area;
    int nodes_number = 3;
    int dof = 2;

    double mean_pressure_rate = N(0) * (GetGeometry()[0].FastGetSolutionStepValue(AIR_PRESSURE_DT));
    for ( int ii=1; ii < nodes_number; ++ii)
        mean_pressure_rate +=  N(ii) * (GetGeometry()[ii].FastGetSolutionStepValue(AIR_PRESSURE_DT));

    mean_pressure_rate *= stbl_fac;

    for (int ii = 0; ii < nodes_number; ++ii)
    {
        int index = ii * (dof + 1);
        F[index] -=  DN_DX(ii,0) * mean_pressure_rate;
        F[index+1] -=  DN_DX(ii,1) * mean_pressure_rate;
    }
// 	  double VC2;
// 	  CalculateSoundVelocity(GetGeometry(), VC2);
// 	  double stbl_fac = tautwo/VC2 * volume;
// 	int nodes_number = 4;
// 	int dof = 3;
// 	for ( int ii = 0; ii < nodes_number; ii++)
// 	    {
// 		//LHS contribution
// 		int row = ii*(dof+1);
// 		for( int jj=0; jj < nodes_number; jj++)
// 		   {
// 			int column = jj*(dof+1) + dof;
//
// 			K(row,column) +=  stbl_fac* N(jj) * DN_DX(ii,0);
// 			K(row + 1,column) +=  stbl_fac* N(jj) * DN_DX(ii,1);
// 			K(row + 2,column) +=  stbl_fac* N(jj) * DN_DX(ii,2);
//
// 		   }
// 	    }

    KRATOS_CATCH("")
}
//*************************************************************************************
//*************************************************************************************
void ASGSCompressible2D::GetFirstDerivativesVector(Vector& values, int Step)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * (dim + 1);
    if (values.size() != MatSize)   values.resize(MatSize,false);
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        unsigned int index = i * (dim + 1);
        values[index] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_X,Step);
        values[index + 1] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_Y,Step);
        values[index + 2] = GetGeometry()[i].GetSolutionStepValue(AIR_PRESSURE,Step);

    }
}
//************************************************************************************
//************************************************************************************
void ASGSCompressible2D::GetSecondDerivativesVector(Vector& values, int Step)
{
    const unsigned int number_of_nodes = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();
    unsigned int MatSize = number_of_nodes * (dim + 1);
    if (values.size() != MatSize) values.resize(MatSize,false);

    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        unsigned int index = i * (dim + 1);
        values[index] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_X,Step);
        values[index + 1] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Y,Step);
        values[index + 2] = GetGeometry()[i].GetSolutionStepValue(AIR_PRESSURE_DT,Step);
    }

}
//************************************************************************************
//************************************************************************************
void ASGSCompressible2D::CalculateSoundVelocity(Geometry< Node<3> > geom, double& vc2)
{

    /*	vc2 = 1441.000 * 1441.000;
      double K1 = 2070000000;
          double K2 = 7.15;

    double mean_rho_w = 0.0;
    double mean_old_rho_w = 0.0;
    double mean_old_pr_w = 0.0;

    for(int ii = 0; ii<3; ii++)
       {
    	mean_rho_w += GetGeometry()[ii].FastGetSolutionStepValue(DENSITY );
    	mean_old_rho_w += GetGeometry()[ii].FastGetSolutionStepValue(DENSITY,1 );
    	mean_old_pr_w += GetGeometry()[ii].FastGetSolutionStepValue(PRESSURE,1);
       }

    mean_rho_w *= 0.333333333333333333333333333333333333;
    mean_old_rho_w *=0.333333333333333333333333333333333333;
    mean_old_pr_w *=0.333333333333333333333333333333333333;


    double alpha = (mean_old_pr_w * K2 + K1)/mean_old_rho_w;

    vc2 = alpha*pow(mean_rho_w/mean_old_rho_w, K2-1.0);*/


    double mean_vc2=0.0;
    double mean_vc=0.0;
    double max_vc = 0.0;
    for (int ii = 0; ii<3; ii++)
    {
        mean_vc = GetGeometry()[ii].FastGetSolutionStepValue(AIR_SOUND_VELOCITY ) ;
        if (mean_vc > max_vc)
            max_vc = mean_vc;

        mean_vc2 += (mean_vc) * (mean_vc);
    }

    vc2 = max_vc*max_vc;
// 	 vc2 =mean_vc2*0.333333333333333333333333;

// 	  vc2 = 5.0/9.0;

}
//*************************************************************************************
//*************************************************************************************
void ASGSCompressible2D::CalculateCompressibleStblTerms(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,const array_1d<double,3>& N,const double tautwo,const double area)
{

    //double VC2;
    //CalculateSoundVelocity(GetGeometry(), VC2);
    double stbl_fac = tautwo * area;
    int nodes_number = 3;
    int dof = 2;
    for ( int ii = 0; ii < nodes_number; ii++)
    {
        //LHS contribution
        int row = ii*(dof+1);
        for ( int jj=0; jj < nodes_number; jj++)
        {
            int column = jj*(dof+1) + dof;

            K(row,column) +=  stbl_fac* N(jj) * DN_DX(ii,0);
            K(row + 1,column) +=  stbl_fac* N(jj) * DN_DX(ii,1);

        }
    }

}

//************************************************************************************
//************************************************************************************
void ASGSCompressible2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dim = 2;
    unsigned int node_size = dim+1;


    if (rResult.size() != number_of_nodes*node_size)
        rResult.resize(number_of_nodes*node_size,false);

    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        rResult[i*node_size] = GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
        rResult[i*node_size+1] = GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
        rResult[i*node_size+2] = GetGeometry()[i].GetDof(AIR_PRESSURE).EquationId();
    }
    KRATOS_CATCH("")

}

//************************************************************************************
//************************************************************************************
void ASGSCompressible2D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dim = 2;
    unsigned int node_size = dim+1;


    if (ElementalDofList.size() != number_of_nodes*node_size)
        ElementalDofList.resize(number_of_nodes*node_size);

    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        ElementalDofList[i*node_size] = GetGeometry()[i].pGetDof(VELOCITY_X);
        ElementalDofList[i*node_size+1] = GetGeometry()[i].pGetDof(VELOCITY_Y);
        ElementalDofList[i*node_size+2] = GetGeometry()[i].pGetDof(AIR_PRESSURE);
    }
    KRATOS_CATCH("");

}
//*************************************************************************************
//*************************************************************************************

void ASGSCompressible2D::CalculateTau(const array_1d<double,3>& N, double& tauone, double& tautwo, const double time, const double area, const ProcessInfo& rCurrentProcessInfo)
{

    KRATOS_TRY
    const array_1d<double, 3 > & adv_vel0 = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY, 0);
    const array_1d<double, 3 > & mesh_vel0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
    const array_1d<double, 3 > & adv_vel1 = GetGeometry()[1].FastGetSolutionStepValue(VELOCITY, 0);
    const array_1d<double, 3 > & mesh_vel1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
    const array_1d<double, 3 > & adv_vel2 = GetGeometry()[2].FastGetSolutionStepValue(VELOCITY, 0);
    const array_1d<double, 3 > & mesh_vel2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);

    array_1d<double, 2 > ms_adv_vel;

    ms_adv_vel[0] = N[0]*(adv_vel0[0] - mesh_vel0[0]) + N[1]*(adv_vel1[0] - mesh_vel1[0]) + N[2]*(adv_vel2[0] - mesh_vel2[0]);
    ms_adv_vel[1] = N[0]*(adv_vel0[1] - mesh_vel0[1]) + N[1]*(adv_vel1[1] - mesh_vel1[1]) + N[2]*(adv_vel2[1] - mesh_vel2[1]);

    //ms_adv_vel[0] = 0.0;
    //ms_adv_vel[1] = 0.0;


    double advvel_norm = ms_adv_vel[0] * ms_adv_vel[0] + ms_adv_vel[1] * ms_adv_vel[1];
    advvel_norm = sqrt(advvel_norm);

    double ele_length = 2.0 * sqrt(area / 3.00);

    double mu;
    //const double mu0 = GetGeometry()[0].FastGetSolutionStepValue(VISCOSITY);
    //const double mu1 = GetGeometry()[1].FastGetSolutionStepValue(VISCOSITY);
    //const double mu2 = GetGeometry()[2].FastGetSolutionStepValue(VISCOSITY);
    //mu = 0.333333333333333333333333*(mu0 + mu1 + mu2);

    double density;
    calculatedensity(GetGeometry(), density, mu);

    double VC2;
    CalculateSoundVelocity(GetGeometry(), VC2);
    double vc = sqrt(VC2);

// 	double int_time = ele_length/vc;
//
//         const double dyn_st_beta = rCurrentProcessInfo[DYNAMIC_TAU];
//         //tauone = 1.0 / (dyn_st_beta / time + 4.0 * mu / (ele_length * ele_length * density) + 2.0 * advvel_norm  / ele_length);
//         tauone = 1.0 / (dyn_st_beta / int_time + 4.0 * mu / (ele_length * ele_length * density) + 2.0 * advvel_norm  / ele_length);

    tauone = time*VC2;
    tautwo = time;

    KRATOS_CATCH("")

}

//*************************************************************************************
//*************************************************************************************
void ASGSCompressible2D::calculatedensity(Geometry< Node<3> > geom, double& density, double& viscosity)
{

    /*double kk = 0.0;
    for(int ii=0;ii<3;++ii)
    	if(geom[ii].GetSolutionStepValue(IS_STRUCTURE) != 1.0)
    		{
    			kk++;
    			density +=geom[ii].FastGetSolutionStepValue(DENSITY);
    		}

    density/=kk;*/
    /*
    	density = ZeroVector(3);
    for(int ii=0;ii<3;++ii)
    	density[ii] = geom[ii].FastGetSolutionStepValue(DENSITY);*/

    const double rho0 = geom[0].FastGetSolutionStepValue(DENSITY_AIR);
    const double rho1 = geom[1].FastGetSolutionStepValue(DENSITY_AIR);
    const double rho2 = geom[2].FastGetSolutionStepValue(DENSITY_AIR);
    density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );

    const double visc0 = geom[0].FastGetSolutionStepValue(VISCOSITY_AIR);
    const double visc1 = geom[1].FastGetSolutionStepValue(VISCOSITY_AIR);
    const double visc2 = geom[2].FastGetSolutionStepValue(VISCOSITY_AIR);
    viscosity = 0.3333333333333333333333*(visc0 + visc1 + visc2 );

    //density = 10.0;

    /*
    	double first = geom[0].FastGetSolutionStepValue(IS_WATER);
    	double second = geom[1].FastGetSolutionStepValue(IS_WATER);
    	double third = geom[2].FastGetSolutionStepValue(IS_WATER);

    	density = 0.0;

    	if(first == second && second==third)
    	  {
    		//for inside the domain totally inside one fluid
    		density = geom[0].FastGetSolutionStepValue(DENSITY);
    		viscosity = geom[0].FastGetSolutionStepValue(VISCOSITY);
    	  }
    	else
    	  {
    		//for element having common node between two fluids or boundary element with IS_WATER==1 inside the domain
    		for(int ii=0;ii<3;++ii)
    			{
    			  if(geom[ii].GetSolutionStepValue(IS_WATER) == 1.0 && geom[ii].GetSolutionStepValue(IS_STRUCTURE) != 1.0)
    				{
    			  	 density = geom[ii].FastGetSolutionStepValue(DENSITY_WATER);
    				 viscosity = geom[ii].FastGetSolutionStepValue(VISCOSITY_WATER);

    				}
    			}
    		//for boundary element with IS_WATER==1 on the boundary
    		if(density == 0.0)
    			for(int ii=0;ii<3;++ii)
    			 {
    			  if(geom[ii].GetSolutionStepValue(IS_WATER) == 0.0)
    				{
    				 density = geom[ii].FastGetSolutionStepValue(DENSITY_AIR);
    				 viscosity = geom[ii].FastGetSolutionStepValue(VISCOSITY_AIR);


    				}
    			 }

    	  }
    */

    //Here we calculate Dynamic viscosity from Kinemeatic viscosity
    viscosity *= density;
}
//************************************************************************************
//************************************************************************************

void ASGSCompressible2D::AddBodyForceAndMomentum(VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& msDN_DX, const array_1d<double, 3 > & N, const double time, const double area, const double tauone, const double tautwo)
{
    KRATOS_TRY
    int nodes_number = 3;
    int dof = 2;


    //double lump_mass_fac = area * 0.333333333333333333333333;

    double density;
    double mu;
    calculatedensity(GetGeometry(), density, mu);
// KRATOS_WATCH(density);
    //body  & momentum term force
    for (int ii = 0; ii < nodes_number; ii++)
    {
        int index = ii * (dof + 1);
        //int loc_index = ii * dof;
        const array_1d<double, 2 > bdf = GetGeometry()[ii].FastGetSolutionStepValue(BODY_FORCE);

        F[index] += area * N[ii] * density * bdf[0];
        F[index + 1] += area * N[ii] * density * bdf[1];

    }


    KRATOS_CATCH("")
}
//************************************************************************************
//************************************************************************************

void ASGSCompressible2D::CalculateResidual(const MatrixType& K, VectorType& F)
{
    KRATOS_TRY

    int nodes_number = 3;
    int dof = 2;


    array_1d<double, 9 > UP = ZeroVector(9);
    for (int ii = 0; ii < nodes_number; ++ii)
    {
        int index = ii * (dof + 1);
        UP[index] = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY, 0)[0];
        UP[index + 1] = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY, 0)[1];
        UP[index + 2] = GetGeometry()[ii].FastGetSolutionStepValue(AIR_PRESSURE, 0);
    }

    F -= prod(K, UP);

    KRATOS_CATCH("")
}
//************************************************************************************
//************************************************************************************

void ASGSCompressible2D::Calculate(const Variable<array_1d<double, 3 > >& rVariable,
                                   array_1d<double, 3 > & Output,
                                   const ProcessInfo& rCurrentProcessInfo)
{
    Output = ZeroVector(3);

    double Area;
    boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX = ZeroMatrix(3,2);
    array_1d<double,3> N = ZeroVector(3);
    GeometryUtils::CalculateGeometryData(GetGeometry(),DN_DX,N,Area);


    double lump_mass_fac = Area * 0.333333333333333333333333;
    double density;
    double mu;

    calculatedensity(GetGeometry(), density, mu);

    //fill velocity and pressure mass

    Output[0] = density*lump_mass_fac;
    Output[1] = lump_mass_fac;
//KRATOS_WATCH(Output[0]);

}
//************************************************************************************
//************************************************************************************
void ASGSCompressible2D::CalculatePressureTerm(MatrixType& K, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const array_1d<double, 3 > & N, const double time, const double area)
{
    KRATOS_TRY
    int nodes_number = 3;
    int dof = 2;

    double density;
    double mu;
    calculatedensity(GetGeometry(), density, mu);

    double VC2;
    CalculateSoundVelocity(GetGeometry(), VC2);
// KRATOS_WATCH(VC2);
    for (int ii = 0; ii < nodes_number; ii++)
    {
        int row = ii * (dof + 1);
        for (int jj = 0; jj < nodes_number; jj++)
        {
            int column = jj * (dof + 1) + dof;

            K(row, column) -=  area * N(jj) * DN_DX(ii, 0);
            K(column, row) += VC2 * area * density * N(jj) * DN_DX(ii, 0);

            K(row + 1, column) -=  area * N(jj) * DN_DX(ii, 1);
            K(column, row + 1) += VC2 * area * density * N(jj) * DN_DX(ii, 1);
        }
    }


    KRATOS_CATCH("")
}
//************************************************************************************
//************************************************************************************
void ASGSCompressible2D::Calculate( const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
{

    Output = 100.0;
    double Area = 0.0;

    boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX = ZeroMatrix(3, 2);
    array_1d<double, 3 > N = ZeroVector(3); //dimension = number of nodes
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);

    double calc_t = 0.0;
    for (unsigned int ii = 0; ii < 3; ii++)
    {
        double inv_max_h = DN_DX(ii,0)*DN_DX(ii,0) + DN_DX(ii,1)*DN_DX(ii,1);
        double VC = GetGeometry()[ii].FastGetSolutionStepValue(AIR_SOUND_VELOCITY ) ;
        inv_max_h = sqrt(inv_max_h);
        if (VC!=0.0)
            calc_t = 1.0/(inv_max_h * VC );
        else
            calc_t = 2*Output;

        if ( calc_t < Output)
            Output = calc_t;

    }
//KRATOS_WATCH(Output);

}
//************************************************************************************
//************************************************************************************
void ASGSCompressible2D::CalcualteDCOperator(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double area)
{
    KRATOS_TRY
    int nodes_number = 3;
    int dof = 2;

    double Vel_art_visc = 0.0;
    double Pr_art_visc = 0.0;
    CalculateArtifitialViscosity(Vel_art_visc,Pr_art_visc,DN_DX);

    double mu;
    double density;
    calculatedensity(GetGeometry(), density, mu);

    double Vel_fac = Vel_art_visc * density * area;
    double Pr_fac = Pr_art_visc*area;

    boost::numeric::ublas::bounded_matrix<double, 3, 3 > gard_opr = ZeroMatrix(nodes_number, nodes_number);
    gard_opr =  prod(DN_DX, trans(DN_DX));

    //nu = nu/density;


    for (int ii = 0; ii < nodes_number; ii++)
    {
        int row = ii * (dof + 1);
        int grad_row = row+dof;
        for (int jj = 0; jj < nodes_number; jj++)
        {
	  Output = ZeroVector(3);

	double Area;
        boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX = ZeroMatrix(3,2);
        array_1d<double,3> N = ZeroVector(3);
	GeometryUtils::CalculateGeometryData(GetGeometry(),DN_DX,N,Area);


	  double lump_mass_fac = Area * 0.333333333333333333333333;
	  double density;
	  double mu;

	  calculatedensity(GetGeometry(), density, mu);

	  //fill velocity and pressure mass

	  Output[0] = density*lump_mass_fac;
	  Output[1] = lump_mass_fac;
//KRATOS_WATCH(Output[0]);

       }
	//************************************************************************************
	//************************************************************************************
    void ASGSCompressible2D::CalculatePressureTerm(MatrixType& K, const boost::numeric::ublas::bounded_matrix<double, 3, 2 > & DN_DX, const array_1d<double, 3 > & N, const double time, const double area) {
        KRATOS_TRY
                int nodes_number = 3;
        int dof = 2;

        double density;
        double mu;
        calculatedensity(GetGeometry(), density, mu);

	double VC2;
	CalculateSoundVelocity(GetGeometry(), VC2);
// KRATOS_WATCH(VC2);
        for (int ii = 0; ii < nodes_number; ii++) {
            int row = ii * (dof + 1);
            for (int jj = 0; jj < nodes_number; jj++) {
                int column = jj * (dof + 1) + dof;

                K(row, column) -=  area * N(jj) * DN_DX(ii, 0);
                K(column, row) += VC2 * area * density * N(jj) * DN_DX(ii, 0);

                K(row + 1, column) -=  area * N(jj) * DN_DX(ii, 1);
                K(column, row + 1) += VC2 * area * density * N(jj) * DN_DX(ii, 1);
            }
        }


        KRATOS_CATCH("")
    }
	//************************************************************************************
	//************************************************************************************
        void ASGSCompressible2D::Calculate( const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
       {

	   Output = 100.0;
	   double Area = 0.0;

	  boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX = ZeroMatrix(3, 2);
	  array_1d<double, 3 > N = ZeroVector(3); //dimension = number of nodes
	  GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);

	double calc_t = 0.0;
        for (unsigned int ii = 0; ii < 3; ii++)
        {
            double inv_max_h = DN_DX(ii,0)*DN_DX(ii,0) + DN_DX(ii,1)*DN_DX(ii,1);
            double VC = GetGeometry()[ii].FastGetSolutionStepValue(AIR_SOUND_VELOCITY ) ;
	    inv_max_h = sqrt(inv_max_h);
	    if(VC!=0.0)
	     calc_t = 1.0/(inv_max_h * VC );
	    else
	      calc_t = 2*Output;	
	    
	    if( calc_t < Output)
		  Output = calc_t;
              
        }
//KRATOS_WATCH(Output);

        }
        //************************************************************************************
    //************************************************************************************
     void ASGSCompressible2D::CalcualteDCOperator(MatrixType& K,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX, const double area)
      {
        KRATOS_TRY 
        int nodes_number = 3;
        int dof = 2;

	double Vel_art_visc = 0.0;
        double Pr_art_visc = 0.0;
	CalculateArtifitialViscosity(Vel_art_visc,Pr_art_visc,DN_DX);
	
        double mu;
        double density;
        calculatedensity(GetGeometry(), density, mu);

       double Vel_fac = Vel_art_visc * density * area;
       double Pr_fac = Pr_art_visc*area;
       
        boost::numeric::ublas::bounded_matrix<double, 3, 3 > gard_opr = ZeroMatrix(nodes_number, nodes_number);
        gard_opr =  prod(DN_DX, trans(DN_DX));

        //nu = nu/density;


        for (int ii = 0; ii < nodes_number; ii++) {
            int row = ii * (dof + 1);
	    int grad_row = row+dof;
            for (int jj = 0; jj < nodes_number; jj++) {
                int column = jj * (dof + 1);
		int gard_column = column + dof;
                K(row, column) +=  Vel_fac * (DN_DX(ii, 0) * DN_DX(jj, 0) + 0.5 * DN_DX(ii, 1) * DN_DX(jj, 1));
                K(row, column + 1) +=  Vel_fac * 0.5 * DN_DX(ii, 1) * DN_DX(jj, 0);		
		K(row + 1, column ) +=  Vel_fac * 0.5 * DN_DX(ii, 0) * DN_DX(jj, 1);	
                K(row + 1, column + 1) +=  Vel_fac * (DN_DX(ii, 1) * DN_DX(jj, 1) + 0.5*DN_DX(ii, 0) * DN_DX(jj, 0));
		
		K(grad_row, gard_column) += Pr_fac * gard_opr(ii, jj);
		
            }
        }
                    
	    KRATOS_CATCH("")
      }
    //************************************************************************************
    //************************************************************************************
     void ASGSCompressible2D::CalculateArtifitialViscosity(double& Vel_art_visc ,double& Pr_art_visc ,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX)
 	{
	    KRATOS_TRY    
	  
	 Vel_art_visc = 0.0; 
	 Pr_art_visc = 0.0;
         int nodes_number = 3;
	 double div_vel = 0.0;
         for( int ii=0; ii < nodes_number; ++ii)
	 {
	    const array_1d<double,3>& vel = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY);
	    div_vel += (vel[0]*DN_DX(ii,0) + vel[1]*DN_DX(ii,1));
	    
	 }  
	 
	 double H=0.0;
	 double norm_grad_p = 0.0;
	 
	double mu;
        double density;
        calculatedensity(GetGeometry(), density, mu);
	 
	 if( div_vel < 0.0)
	    {
             CalculateCharectristicLength(H,DN_DX,norm_grad_p);	
             Vel_art_visc =  0.8* abs(div_vel) * pow(H,2);  //Noh: 1d 0.4 2d 0.65 Sedov_2d 0.5  last:20.0 * 0.5

	     Pr_art_visc = 0.6*sqrt(norm_grad_p/density) * pow(H,1.5);  //Noh: 1d 0.3  2d 0.3 sedov_2d 0.3 last 10.0 * 0.3
	    } 
	   
	   this->GetValue(VEL_ART_VISC)=Vel_art_visc;
	   this->GetValue(PR_ART_VISC)=Pr_art_visc;
	    
	    
	    KRATOS_CATCH("")
	} 
	
    //************************************************************************************
    //************************************************************************************
     void ASGSCompressible2D::CalculateCharectristicLength(double& ch_length, const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,double& norm_grad )
 	{
	    KRATOS_TRY 
	    
  	  GeometryType::JacobiansType J;
	  GetGeometry().Jacobian(J); 

	  Matrix CC;
	  CC.resize(2,2,false);	  	
	  
	  Matrix DD;
	  DD.resize(2,2,false);	  	  

	  noalias(CC) = prod(J[0],trans(J[0]));	 
	  DD=CC;
	  
	  double det = CC(1,1)*CC(0,0) - CC(0,1)*CC(1,0);
	  
	  if(det == 0.0)
	      	KRATOS_ERROR(std::logic_error,"ZERO DETERMINANT IN ARTIFICIAL VISCOSITY ",det)
	  else 
	     det = 1.0/det;
	  
          double zarf = CC(0,0);	  
	  CC(0,0) = det*CC(1,1);
	  CC(1,1) = det*zarf;
	  CC(0,1) = -CC(0,1)*det;
	  CC(1,0) = -CC(1,0)*det;
	  //KRATOS_WATCH( prod(DD,CC));	  

	  int nodes_number = 3;
	  array_1d<double,3> mean_acc =  GetGeometry()[0].FastGetSolutionStepValue(ACCELERATION);
	  double rho = GetGeometry()[0].FastGetSolutionStepValue(DENSITY_AIR);
	  double pr = GetGeometry()[0].FastGetSolutionStepValue(AIR_PRESSURE);
	  array_1d<double,3> grad_rho,grad_pr;
	  grad_rho[0] = DN_DX(0,0)*rho;
	  grad_rho[1] = DN_DX(0,1)*rho;
	  
	  grad_pr[0] = DN_DX(0,0)*pr;
	  grad_pr[1] = DN_DX(0,1)*pr;	  
	  
		  
          for (int ii = 1; ii < nodes_number; ii++) {	  
	        mean_acc += GetGeometry()[ii].FastGetSolutionStepValue(ACCELERATION);  
		rho = GetGeometry()[ii].FastGetSolutionStepValue(DENSITY_AIR);
		pr = GetGeometry()[ii].FastGetSolutionStepValue(AIR_PRESSURE);		
		grad_rho[0] += DN_DX(ii,0)*rho;
		grad_rho[1] += DN_DX(ii,1)*rho;	
		grad_pr[0] += DN_DX(ii,0)*pr;
		grad_pr[1] += DN_DX(ii,1)*pr;			
	  }
	  mean_acc *= 0.33333333333333333333333333333333333;
	  grad_rho *= 0.33333333333333333333333333333333333;
	  grad_pr  *= 0.33333333333333333333333333333333333;	  
	  grad_rho[2] = 0.0;
	  grad_pr[2] = 0.0;
	  
	  double norm_acc =  MathUtils<double>::Norm3(mean_acc);
	  double norm_grad_rho =  MathUtils<double>::Norm3(grad_rho);
	  norm_grad =  MathUtils<double>::Norm3(grad_pr);

	  array_1d<double,3>  n_dir= ZeroVector(3);
	  if(norm_acc != 0.0)
	             n_dir = 0.75/norm_acc*mean_acc;
	  if(norm_grad_rho != 0.0)
	             n_dir +=  0.25/norm_grad_rho*grad_rho;
	  
	  double norm_n_dir =  MathUtils<double>::Norm3(n_dir);	
 	  
	  if(norm_n_dir != 0.0)
	    {
	      n_dir/=norm_n_dir;	  
	      array_1d<double,3>  CC_n;	  
	      CC_n = prod(CC,n_dir);	  
	      double denom = n_dir[0]*CC_n[0] + n_dir[1]*CC_n[1] + n_dir[2]*CC_n[2];

	      if(denom <= 0.0)
		    KRATOS_ERROR(std::logic_error,"CalculateCharectristicLength zero or negative denominator ",denom)	
	      else
		    ch_length = 2.0/sqrt(denom);
	    }
	   else{
	         ch_length = 0.0;
/*		 KRATOS_WATCH("CalculateCharectristicLength zero norm_n_dir ");
		 KRATOS_WATCH(this->Id());
		 KRATOS_WATCH(norm_grad_rho);
		 KRATOS_WATCH(norm_acc);*/		 
	   }
	   
	   this->GetValue(NORMAL) = n_dir;
	   
	    KRATOS_CATCH("")
	}
    //************************************************************************************
    //************************************************************************************

    void ASGSCompressible2D::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo) {

/*        double delta_t = rCurrentProcessInfo[DELTA_TIME];*/
	boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
        array_1d<double, 3 > N;
        //getting data for the given geometry
        double Area;
        GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);

        if (rVariable == VEL_ART_VISC) {
            for (unsigned int PointNumber = 0;
                    PointNumber < 1; PointNumber++) {

                rValues[PointNumber] = this->GetValue(VEL_ART_VISC);
            }
        }
        if (rVariable == PR_ART_VISC) {
            for (unsigned int PointNumber = 0;
                    PointNumber < 1; PointNumber++) {
	//	KRATOS_WATCH(this->GetValue(IS_WATER));
	//	KRATOS_WATCH(this->Info());
                rValues[PointNumber] = this->GetValue(PR_ART_VISC);
            }
        }
         if (rVariable == IS_WATER_ELEMENT) {
            for (unsigned int PointNumber = 0;
                    PointNumber < 1; PointNumber++) {
	//	KRATOS_WATCH(this->GetValue(IS_WATER));
// 		KRATOS_WATCH(this->Info());
                rValues[PointNumber] = this->GetValue(IS_WATER_ELEMENT);
            }
        }       

    }
     	//*************************************************************************************
	//*************************************************************************************	   
    void ASGSCompressible2D::GetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
                std::vector<array_1d<double, 3 > >& rValues,
                const ProcessInfo& rCurrentProcessInfo)
        {
	  if(rVariable == NORMAL)
	    rValues[0] = this->GetValue(NORMAL);
	  
	  
        }
    	//*************************************************************************************
	//*************************************************************************************	  	  
      void ASGSCompressible2D::CalculateNonlinearStblTerm(VectorType& F,const boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,const array_1d<double,3>& N, const double time,const double tautwo,const double area)
	{
	double lump_mass_fac = area * 0.333333333333333333333333;
	int nodes_number = 3;
        int dof = 2;	
	
        double density;
        double mu;
        calculatedensity(GetGeometry(), density, mu);
	
	double VC2;
	CalculateSoundVelocity(GetGeometry(), VC2);
	
         double mean_pressure_rate = N(0) * (GetGeometry()[0].FastGetSolutionStepValue(AIR_PRESSURE_DT));
	 for( int ii=1; ii < nodes_number; ++ii)
            mean_pressure_rate +=  N(ii) * (GetGeometry()[ii].FastGetSolutionStepValue(AIR_PRESSURE_DT));
	 
	 

	 double div_vel = 0.0;
         for( int ii=0; ii < nodes_number; ++ii)
	 {
	    const array_1d<double,3>& vel = GetGeometry()[ii].FastGetSolutionStepValue(VELOCITY);
	    div_vel += (vel[0]*DN_DX(ii,0) + vel[1]*DN_DX(ii,1));
	    
	 }
	 
	 double nonlinear_term = div_vel *(mean_pressure_rate + density*VC2*div_vel);
						
        for (int ii = 0; ii < nodes_number; ++ii) {
            int index = ii * (dof + 1) + dof;	    
            F[index] -= 0.4*tautwo * lump_mass_fac * nonlinear_term; // gamma = 5/3 ---> 0.66666666666667 gamma = 7/5 ----> 0.4
        }	 
	 
	 
	  
	}     
        //************************************************************************************
    //************************************************************************************




} // Namespace Kratos


