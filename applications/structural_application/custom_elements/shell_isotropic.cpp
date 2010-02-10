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
//   Date:                $Date: 2008-02-21 09:36:17 $
//   Revision:            $Revision: 1.7 $
//
//


// System includes 
 
// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/shell_isotropic.h"
#include "structural_application.h"
#include "utilities/math_utils.h" 

                                                                                                                              namespace Kratos
{
	namespace ShellIsotropicAuxiliaries
	{
		array_1d<int,9> local_indices;
		#pragma omp threadprivate(local_indices)
		boost::numeric::ublas::bounded_matrix<double,9,3> mBm = ZeroMatrix(9,3); //membrane displacement-strain matrix
		#pragma omp threadprivate(mBm)
		boost::numeric::ublas::bounded_matrix<double,9,3> mBb = ZeroMatrix(9,3); //bending displacement-strain matrix
		#pragma omp threadprivate(mBb)

		boost::numeric::ublas::bounded_matrix<double,3,3> Q = ZeroMatrix(3,3); 
		#pragma omp threadprivate(Q)
		boost::numeric::ublas::bounded_matrix<double,3,3> Q1 = ZeroMatrix(3,3); 
		#pragma omp threadprivate(Q1)
		boost::numeric::ublas::bounded_matrix<double,3,3> Q2 = ZeroMatrix(3,3); 
		#pragma omp threadprivate(Q2)
		boost::numeric::ublas::bounded_matrix<double,3,3> Q3 = ZeroMatrix(3,3); 
		#pragma omp threadprivate(Q3)
		boost::numeric::ublas::bounded_matrix<double,3,3> aux33 = ZeroMatrix(3,3);
		#pragma omp threadprivate(aux33)
		boost::numeric::ublas::bounded_matrix<double,3,3> Te = ZeroMatrix(3,3);
		#pragma omp threadprivate(Te)

		boost::numeric::ublas::bounded_matrix<double,3,3> mEm = ZeroMatrix(3,3);
		#pragma omp threadprivate(mEm)
		boost::numeric::ublas::bounded_matrix<double,3,3> mEb = ZeroMatrix(3,3);
		#pragma omp threadprivate(mEb)

		array_1d<double,9> H1 = ZeroVector(9); 
		#pragma omp threadprivate(H1)
		array_1d<double,9> H2 = ZeroVector(9); 
		#pragma omp threadprivate(H2)
		array_1d<double,9> H3 = ZeroVector(9);
		#pragma omp threadprivate(H3)
		array_1d<double,9> H4 = ZeroVector(9);
		#pragma omp threadprivate(H4)

		boost::numeric::ublas::bounded_matrix<double,9,3> TTu  = ZeroMatrix(9,3); //attention 9*3 and not 3*9
		#pragma omp threadprivate(TTu)
		boost::numeric::ublas::bounded_matrix<double,3,9> aux39  = ZeroMatrix(3,9);
		#pragma omp threadprivate(aux39)
		
		boost::numeric::ublas::bounded_matrix<double,18,18> mKloc_system  = ZeroMatrix(18,18); //stiffness matrix in the local reference system
		#pragma omp threadprivate(mKloc_system)
		boost::numeric::ublas::bounded_matrix<double,18,18> rot18  = ZeroMatrix(18,18); //TAKE CARE!! this is VERY inefficient
		#pragma omp threadprivate(rot18)

		array_1d<double,3> temp = ZeroVector(3);
		#pragma omp threadprivate(temp)
		boost::numeric::ublas::bounded_matrix<double,9,9> mKloc99 = ZeroMatrix(9,9);
		#pragma omp threadprivate(mKloc99)
		Vector values = ZeroVector(18); //help vector to get the current values of displacements (and rotations)
		#pragma omp threadprivate(values)
	}
	
	using namespace ShellIsotropicAuxiliaries;
	
	//************************************************************************************
	//************************************************************************************
	ShellIsotropic::ShellIsotropic(IndexType NewId, GeometryType::Pointer pGeometry)
	: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	ShellIsotropic::ShellIsotropic(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
	}

	Element::Pointer ShellIsotropic::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new ShellIsotropic(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	ShellIsotropic::~ShellIsotropic()
	{
	}
	
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		CalculateAllMatrices(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
	}	
	
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		Matrix lhs(18,18);
		CalculateAllMatrices(lhs,rRightHandSideVector,rCurrentProcessInfo);
	}
	

	//************************************************************************************
	//************************************************************************************
	//getting the local coordinates and calculating the base vectors
	void ShellIsotropic::CalculateLocalGlobalTransformation(
		double& x12,
		double& x23,
		double& x31,
		double& y12,
		double& y23,
		double& y31,
		array_1d<double,3>& v1,
		array_1d<double,3>& v2,
		array_1d<double,3>& v3,
		double& area
		)
	{
		KRATOS_TRY
		//forming the base vectors - align with the side 1-2
		v1[0] = GetGeometry()[1].X()-GetGeometry()[0].X(); 
		v1[1] = GetGeometry()[1].Y()-GetGeometry()[0].Y();
		v1[2] = GetGeometry()[1].Z()-GetGeometry()[0].Z();
		double x21=norm_2(v1); //here we assume the basis aligned with 1-2
		double y21 = 0.0; 
		x12=-x21;
		y12=-y21;

		//forming the (area) normal
		temp[0]=GetGeometry()[2].X()-GetGeometry()[0].X();
		temp[1]=GetGeometry()[2].Y()-GetGeometry()[0].Y();
		temp[2]=GetGeometry()[2].Z()-GetGeometry()[0].Z();
		MathUtils<double>::CrossProduct(v3,v1,temp);
		area = 0.5 * norm_2(v3);

		//normalizing base vectors
		v1 /= x21;
		v3 /= (2.0*area);

		//forming the "second" base vector - it is already normalized
		MathUtils<double>::CrossProduct(v2,v3,v1);
		x31 = inner_prod(temp,v1);
		y31 = inner_prod(temp,v2);

		x23 = x21-x31;
		y23 = y21-y31;

		//calculating 
		KRATOS_CATCH("");
	}


	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::CalculateMembraneB( boost::numeric::ublas::bounded_matrix<double,9,3>& Bm,
				const double&  beta0,
				const double& loc1,
				const double& loc2,
				const double& loc3,
				const double& x12,
				const double& x23,
				const double& x31,
				const double& y12,
				const double& y23,
				const double& y31
				)
	{
		KRATOS_TRY
	
		//template parameters
		const double alpha = 1.5;
		const double b1 = 1;
		const double b2 = 2;
		const double b3 = 1;
		const double b4 = 0;
		const double b5 = 1;
		const double b6 = -1;
		const double b7 = -1;
		const double b8 = -1;
		const double b9 = -2;
	
		const double x21 = -x12;
		const double x32 = -x23;
		const double x13 = -x31;
		const double y21 = -y12;
		const double y32 = -y23;
		const double y13 = -y31;
	
		const double A = 0.5*(y21*x13 - x21*y13);
		const double A2 = 2.00*A;
		const double A4 = 4.00*A;
		
		//calculating L and writing it directly in Bm
		double temp = alpha/6.00;
		Bm(0,0) = y23;                       Bm(0,1) = 0.00;                  Bm(0,2) = x32;
		Bm(1,0) = 0.00;                      Bm(1,1) = x32;                   Bm(1,2) = y23;
		Bm(2,0) = y23*(y13-y21)*temp;        Bm(2,1) = x32*(x31-x12)*temp;    Bm(2,2) = 2.00*(x31*y13-x12*y21)*temp;
		Bm(3,0) = y31;                       Bm(3,1) = 0.00;                  Bm(3,2) = x13;
		Bm(4,0) = 0.00;                      Bm(4,1) = x13;                   Bm(4,2) = y31;
		Bm(5,0) = y31*(y21-y32)*temp;        Bm(5,1) = x13*(x12-x23)*temp;    Bm(5,2) = 2.00*(x12*y21-x23*y32)*temp;
		Bm(6,0) = y12;                       Bm(6,1) = 0.00;                  Bm(6,2) = x21;
		Bm(7,0) = 0.00;                      Bm(7,1) = x21;                   Bm(7,2) = y12;
		Bm(8,0) = y12*(y32-y13)*temp;        Bm(8,1) = x21*(x23-x31)*temp;    Bm(8,2) = 2.00*(x23*y32-x31*y13)*temp;
		Bm*=(0.5/(A)); //multiply by h and divide by volume (the h erases)

		//calculation of higher order  contributions............
		//Calculate matrix Te
		double LL21 = x21*x21 + y21*y21;
		double LL32 = x32*x32 + y32*y32;
		double LL13 = x13*x13 + y13*y13;
		Te(0,0) = y23*y13*LL21;             Te(0,1) = y31*y21*LL32;             Te(0,2) = y12*y32*LL13;
		Te(1,0) = x23*x13*LL21;             Te(1,1) = x31*x21*LL32;             Te(1,2) = x12*x32*LL13;
		Te(2,0) = (y23*x31+x32*y13)*LL21;   Te(2,1) = (y31*x12+x13*y21)*LL32;   Te(2,2) = (y12*x23+x21*y32)*LL13;
		Te /= (A*A4);
	
		//matrices Q
		Q1(0,0) = b1*A2/(LL21*3.00); Q1(0,1) = b2*A2/(LL21*3.00); Q1(0,2) = b3*A2/(LL21*3.00);
		Q1(1,0) = b4*A2/(LL32*3.00); Q1(1,1) = b5*A2/(LL32*3.00); Q1(1,2) = b6*A2/(LL32*3.00);
		Q1(2,0) = b7*A2/(LL13*3.00); Q1(2,1) = b8*A2/(LL13*3.00); Q1(2,2) = b9*A2/(LL13*3.00);
	
		Q2(0,0) = b9*A2/(LL21*3.00); Q2(0,1) = b7*A2/(LL21*3.00); Q2(0,2) = b8*A2/(LL21*3.00);
		Q2(1,0) = b3*A2/(LL32*3.00); Q2(1,1) = b1*A2/(LL32*3.00); Q2(1,2) = b2*A2/(LL32*3.00);
		Q2(2,0) = b6*A2/(LL13*3.00); Q2(2,1) = b4*A2/(LL13*3.00); Q2(2,2) = b5*A2/(LL13*3.00);
	
		Q3(0,0) = b5*A2/(LL21*3.00); Q3(0,1) = b6*A2/(LL21*3.00); Q3(0,2) = b4*A2/(LL21*3.00);
		Q3(1,0) = b8*A2/(LL32*3.00); Q3(1,1) = b9*A2/(LL32*3.00); Q3(1,2) = b7*A2/(LL32*3.00);
		Q3(2,0) = b2*A2/(LL13*3.00); Q3(2,1) = b3*A2/(LL13*3.00); Q3(2,2) = b1*A2/(LL13*3.00);
	
		noalias(Q) = loc1 * Q1;
		noalias(Q) += loc2 * Q2;
		noalias(Q) += loc3 * Q3;
	
		//constructing ttilda
		for(unsigned int i=0; i<3; i++)
		{
			TTu(0,i) = x32;
			TTu(1,i) = y32;
			TTu(2,i) = 0.0;
		
			TTu(3,i) = x13;
			TTu(4,i) = y13;
			TTu(5,i) = 0.0;
		
			TTu(6,i) = x21;
			TTu(7,i) = y21;
			TTu(8,i) = 0.0;
		}
		TTu(2,0) = A4;
		TTu(5,1) = A4;
		TTu(8,2) = A4;
		TTu *= (1.0/A4);

		//adding the higher order stiffness to Bm
		//Bm = L + 3/2*sqrt(beta0) * TTu * Q * Te;
		noalias(aux33) = (1.5*sqrt(beta0)) * prod( trans(Q),trans(Te) );
		noalias(Bm) += prod(TTu,aux33);

		KRATOS_CATCH("")
	}


	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::CalculateBendingB( boost::numeric::ublas::bounded_matrix<double,9,3>& Bb,
				const double& loc2,
				const double& loc3,
				const double& x12,
				const double& x23,
				const double& x31,
				const double& y12,
				const double& y23,
				const double& y31
				)
	{
		KRATOS_TRY

		//calculate auxiliary parameters
		double LL12 = x12*x12 + y12*y12;
		double LL23 = x23*x23 + y23*y23;
		double LL31 = x31*x31 + y31*y31;
	
		const double p4 = -6.0*x23/(LL23);
		const double p5 = -6.0*x31/(LL31);
		const double p6 = -6.0*x12/(LL12);
	
		const double t4 = -6.0*y23/(LL23);
		const double t5 = -6.0*y31/(LL31);
		const double t6 = -6.0*y12/(LL12);
	
		const double q4 = 3.0*x23*y23/(LL23);
		const double q5 = 3.0*x31*y31/(LL31);
		const double q6 = 3.0*x12*y12/(LL12);
	
		const double r4 = 3.0*y23*y23/(LL23);
		const double r5 = 3.0*y31*y31/(LL31);
		const double r6 = 3.0*y12*y12/(LL12);
		
		const double area = 0.5*(y12*x31 - x12*y31);
		
		//calculating auxiliary vectors H1
		H1[0] = p6*(1.0-2.0*loc2) + (p5-p6)*loc3;
		H1[1] = q6*(1.0-2.0*loc2) - (q5+q6)*loc3;
		H1[2]= -4.0+6.0*(loc2+loc3) + r6*(1.0-2*loc2) - loc3*(r5+r6);
		H1[3] = -p6*(1.0-2.0*loc2) + (p4+p6)*loc3;
		H1[4] = q6*(1.0-2.0*loc2) - (q6-q4)*loc3;
		H1[5]= -2.0+6.0*loc2  + r6*(1-2.0*loc2) + loc3*(r4-r6);
		H1[6] = -(p4+p5)*loc3;
		H1[7] = (q4-q5)*loc3;
		H1[8]= -loc3*(r5-r4);
	
		//calculating auxiliary vectors H2
		H2[0] = t6*(1.0-2.0*loc2) + (t5-t6)*loc3;
		H2[1] = 1.0 + r6*(1.0-2.0*loc2) - loc3*(r5+r6);
		H2[2]= -q6*(1.0-2.0*loc2)+loc3*(q5+q6);
		H2[3] = -t6*(1.0-2.0*loc2) + (t4+t6)*loc3;
		H2[4] = -1.0 + r6*(1-2*loc2) + loc3*(r4-r6);
		H2[5]= -q6*(1.0-2.0*loc2)-loc3*(q4-q6);
		H2[6] = -(t5+t4)*loc3;
		H2[7] = (r4-r5)*loc3;
		H2[8]= -loc3*(q4-q5);
	
		//calculating auxiliary vectors H3
		H3[0] = -p5*(1.0-2.0*loc3) - (p6-p5)*loc2;
		H3[1] = q5*(1.0-2.0*loc3) - (q5+q6)*loc2;
		H3[2]= -4.0 + 6.0*(loc2+loc3) + r5*(1.0-2.0*loc3) - loc2*(r5+r6);
		H3[3] = loc2*(p4+p6);
		H3[4] = loc2*(q4-q6);
		H3[5]= -loc2*(r6-r4);
		H3[6] = p5*(1.0-2.0*loc3) - (p4+p5)*loc2;
		H3[7] = q5*(1.0-2.0*loc3)+loc2*(q4-q5);
		H3[8]= -2.0+6.0*loc3+r5*(1.0-2.0*loc3)+loc2*(r4-r5);
	
		//calculating auxiliary vectors H1
		H4[0] = -t5*(1.0-2.0*loc3) - (t6-t5)*loc2;
		H4[1] = 1.0+r5*(1.0-2.0*loc3)-loc2*(r5+r6);
		H4[2]= -q5*(1.0-2.0*loc3)+loc2*(q5+q6);
		H4[3] = (t4+t6)*loc2;
		H4[4] = (r4-r6)*loc2;
		H4[5]= -loc2*(q4-q6);
		H4[6] = t5*(1.0-2.0*loc3) - (t5+t4)*loc2;
		H4[7] = -1.0+r5*(1.0-2.0*loc3)+loc2*(r4-r5);
		H4[8]= -q5*(1.0-2.0*loc3)-loc2*(q4-q5);
	
		//forming Bb_trans
		double temp = 0.5 / area;
		for(unsigned int i =0;i<9;i++)
		{
			Bb(i,0) = temp * (y31*H1[i] + y12*H3[i]);
			Bb(i,1) = temp * (-x31*H2[i] - x12*H4[i]);
			Bb(i,2) = temp * (-x31*H1[i] - x12*H3[i] + y31*H2[i] + y12*H4[i] );
		}
		
		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::CalculateMembraneContribution( 
					const boost::numeric::ublas::bounded_matrix<double,9,3>& Bm,
					const boost::numeric::ublas::bounded_matrix<double,3,3>& Em,
					boost::numeric::ublas::bounded_matrix<double,9,9>& Km)
	{
		KRATOS_TRY
		noalias(aux39) = prod(Em,trans(Bm) );
		noalias(Km) = prod(Bm, aux39);
		
		KRATOS_CATCH("")
	}
	
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::AssembleMembraneContribution( 
					const boost::numeric::ublas::bounded_matrix<double,9,9>& Km,
					const double& coeff, 
					boost::numeric::ublas::bounded_matrix<double,18,18>& Kloc_system )
	{
		KRATOS_TRY
	
		local_indices[0] = 0;
		local_indices[1] = 1;
		local_indices[2] = 5;
		local_indices[3] = 6;
		local_indices[4] = 7;
		local_indices[5] = 11;
		local_indices[6] = 12;
		local_indices[7] = 13;
		local_indices[8] = 17;
	
		for(unsigned int i = 0; i<9; i++)
			{
			for(unsigned int j = 0; j<9; j++)
			{
				Kloc_system( local_indices[i],local_indices[j] ) += coeff*Km(i,j);
			}
		}
		KRATOS_CATCH("")
	}
	
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::CalculateBendingContribution(
					const boost::numeric::ublas::bounded_matrix<double,9,3>& Bb, 
					const boost::numeric::ublas::bounded_matrix<double,3,3>& Eb, 
       					boost::numeric::ublas::bounded_matrix<double,9,9>& Kb
					 )
	{
		KRATOS_TRY
		noalias(aux39) = prod(Eb,trans(Bb) );
		noalias(Kb) = prod(Bb, aux39);
		KRATOS_CATCH("")
	}
	
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::AssembleBendingContribution( 
					const boost::numeric::ublas::bounded_matrix<double,9,9>& Kb, 
					const double& coeff, 
					boost::numeric::ublas::bounded_matrix<double,18,18>& Kloc_system )
	{
		KRATOS_TRY
	
		local_indices[0] = 2;
		local_indices[1] = 3;
		local_indices[2] = 4;
		local_indices[3] = 8;
		local_indices[4] = 9;
		local_indices[5] = 10;
		local_indices[6] = 14;
		local_indices[7] = 15;
		local_indices[8] = 16;
	
		for(unsigned int i = 0; i<9; i++)
		{
			for(unsigned int j = 0; j<9; j++)
			{
				Kloc_system( local_indices[i],local_indices[j] ) += coeff*Kb(i,j);
			}
		}
		KRATOS_CATCH("")
	}
	
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::CalculateGaussPointContribution(
				boost::numeric::ublas::bounded_matrix<double,18,18>& Kloc_system ,
    				const boost::numeric::ublas::bounded_matrix<double,3,3>& Em,
				const boost::numeric::ublas::bounded_matrix<double,3,3>& Eb,				
				const double& weight,
				const double& h, /*thickness*/
				const double& loc1, /*local coords*/
				const double& loc2,
				const double& loc3,
				const double& x12,
				const double& x23,
				const double& x31,
				const double& y12,
				const double& y23,
				const double& y31
				)
	{
		//membrane stiffness
		double beta0 = CalculateBeta( Em );
		CalculateMembraneB(mBm, beta0,  loc1, loc2, loc3, x12, x23, x31, y12, y23, y31 );	
		CalculateMembraneContribution( mBm, Em, mKloc99);
		AssembleMembraneContribution(   mKloc99, weight, Kloc_system  );
		

		//bending stiffness 
		CalculateBendingB(mBb, loc2, loc3, x12, x23, x31, y12, y23, y31 );
		CalculateBendingContribution( mBb, Eb, mKloc99);

		double bending_weight = weight; //attention!!
		AssembleBendingContribution(   mKloc99, bending_weight, Kloc_system  );
		//membrane+bending stiffness
	
	}
	
	//************************************************************************************
	//************************************************************************************
	double ShellIsotropic::CalculateBeta( const boost::numeric::ublas::bounded_matrix<double,3,3>& Em )
	{
			double nu = GetProperties()[POISSON_RATIO];
			double beta0 = 0.5*(1-4*nu*nu);
			return beta0;
	}
	
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::CalculateMembraneElasticityTensor( 
						boost::numeric::ublas::bounded_matrix<double,3,3>& Em, 
						const double& h 
						)
	{
		double NU = GetProperties()[POISSON_RATIO];
		double E = GetProperties()[YOUNG_MODULUS];
		double coeff = E  * h/ (1 - NU*NU); //we take in account here the thickness
		Em(0,0) = coeff; 
		Em(0,1) = NU*coeff;
		Em(0,2) = 0.0;
	
		Em(1,0) = NU*coeff;
		Em(1,1) = coeff;
		Em(1,2) = 0.0;
	
		Em(2,0) = 0.0;
		Em(2,1) = 0.0;
		Em(2,2) = 0.5*(1-NU)*coeff;
	}
	
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::CalculateBendingElasticityTensor( boost::numeric::ublas::bounded_matrix<double,3,3>& Eb, const double& h )
	{
		double NU = GetProperties()[POISSON_RATIO];
		double E = GetProperties()[YOUNG_MODULUS];
		double coeff = E * h*h*h / (12.0*(1 - NU*NU) );
		Eb(0,0) = coeff; 
		Eb(0,1) = NU*coeff;
		Eb(0,2) = 0.0;
	
		Eb(1,0) = NU*coeff;
		Eb(1,1) = coeff;
		Eb(1,2) = 0.0;
	
		Eb(2,0) = 0.0;
		Eb(2,1) = 0.0;
		Eb(2,2) = 0.5*(1-NU)*coeff;
	}
	
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::CalculateAllMatrices(	
			MatrixType& rLeftHandSideMatrix,
			VectorType& rRightHandSideVector,
			ProcessInfo& rCurrentProcessInfo)
	{
		//calculate local coordinates and rotation matrix
		array_1d<double,3> v1;
		array_1d<double,3> v2;
		array_1d<double,3> v3;
		double x12, x23, x31, y12, y23, y31;
		double A;
	
		CalculateLocalGlobalTransformation( x12, x23, x31, y12, y23, y31,v1,v2,v3,A);
	
		double weight,loc1,loc2,loc3;
	
		//resizing and resetting the system matrix
		if(rLeftHandSideMatrix.size1() != 18)
			rLeftHandSideMatrix.resize(18,18,false);
		if(rRightHandSideVector.size() != 18)
			rRightHandSideVector.resize(18,false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(18,18);
		noalias(rRightHandSideVector) = ZeroVector(18);
		
		noalias(mKloc_system) = ZeroMatrix(18,18);
		
		double h = GetProperties()[THICKNESS];

		//calculating material matrices
		CalculateMembraneElasticityTensor(mEm,h);
		CalculateBendingElasticityTensor(mEb,h);
		
		//here we calculate everything in the local system of coordinates
		//calculate integration point 1
		weight = A/3;
		loc1= 0.5;
		loc2 = 0.5;
		loc3 = 0.0;
		CalculateGaussPointContribution( mKloc_system, mEm, mEb, weight, h, loc1, loc2, loc3, x12, x23, x31, y12, y23, y31);

		//calculate integration point 2
		weight = A/3;
		loc1= 0.0;
		loc2 = 0.5;
		loc3 = 0.5;
		CalculateGaussPointContribution( mKloc_system, mEm, mEb, weight, h, loc1, loc2, loc3, x12, x23, x31, y12, y23, y31);
	
		//calculate integration point 3
		weight = A/3;
		loc1= 0.5;
		loc2 = 0.0;
		loc3 = 0.5;
		CalculateGaussPointContribution( mKloc_system, mEm, mEb, weight, h, loc1, loc2, loc3, x12, x23, x31, y12, y23, y31);
		
		//calculating the contribution of the internal forces (in the local system)
		// to the RHS
		CalculatePureDisplacement(values,v1,v2,v3);
		noalias(rRightHandSideVector) = - prod(mKloc_system,values);		
		
		boost::numeric::ublas::bounded_matrix<double,18,18> ProjOperator, WorkMatrix;
		array_1d<double,18> WorkArray;
		
		//add geometric stiffness to the local system
		Matrix Kg = ZeroMatrix(18,18);
		CalculateAndAddKg(Kg,WorkMatrix,x12,x23,x31,y12,y23,y31,v1,v2,v3,A); //take care.. this should be optimized
		noalias(mKloc_system) += Kg;
		
		CalculateProjectionOperator(ProjOperator, x12, x23, x31, y12, y23, y31);
		ApplyProjection(mKloc_system, rRightHandSideVector, WorkMatrix, WorkArray,ProjOperator );

		//rotate from local coordinates to global coordinates
		RotateToGlobal(v1,v2,v3,mKloc_system,rLeftHandSideMatrix,rRightHandSideVector);
		
		//adding the body force (which is already given in global coordinates)		
		AddBodyForce(h,A,rRightHandSideVector);

		//adding pressure force
		const double one_third = 1.00 / 3.00;
		double element_pressure = (GetGeometry()[0].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE)
		    + GetGeometry()[1].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE)
		    + GetGeometry()[2].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE)) * one_third;

	        array_1d<double,3> pressure_force = element_pressure * A * one_third * v3;

		rRightHandSideVector[0] += pressure_force[0];
		rRightHandSideVector[1] += pressure_force[1];
		rRightHandSideVector[2] += pressure_force[2];
		
		rRightHandSideVector[6] += pressure_force[0];
		rRightHandSideVector[7] += pressure_force[1];
		rRightHandSideVector[8] += pressure_force[2];
		
		rRightHandSideVector[12] += pressure_force[0];
		rRightHandSideVector[13] += pressure_force[1];
		rRightHandSideVector[14] += pressure_force[2];
		
		
	} 
	
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		int number_of_nodes = 3;
		if(rResult.size() != 18)
			rResult.resize(18,false);

		for (int i=0;i<number_of_nodes;i++)
		{
			int index = i*6;
			rResult[index]	   = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
			rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
			rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();

			rResult[index + 3] = GetGeometry()[i].GetDof(ROTATION_X).EquationId();
			rResult[index + 4] = GetGeometry()[i].GetDof(ROTATION_Y).EquationId();
			rResult[index + 5] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
		}

	}

	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		ElementalDofList.resize(0);

		for (unsigned int i=0;i<GetGeometry().size();i++)
		{
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));

			ElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Y));
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Z));
		}
	}

	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::GetValuesVector(Vector& values, int Step)
	{
		const unsigned int number_of_nodes = 3;
		//const unsigned int dim = 3;
		unsigned int MatSize = 18;
		if(values.size() != MatSize)	values.resize(MatSize,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			const array_1d<double,3>& disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,Step);
			const array_1d<double,3>& rot = GetGeometry()[i].FastGetSolutionStepValue(ROTATION,Step);
				
			unsigned int index = i*6;
			values[index] 	  = disp[0];
			values[index + 1] = disp[1];
			values[index + 2] = disp[2];

			values[index + 3] = rot[0];
			values[index + 4] = rot[1];
			values[index + 5] = rot[2];

		}
	}


	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::RotateToGlobal(
			const array_1d<double,3>& v1,
			const array_1d<double,3>& v2,
			const array_1d<double,3>& v3,
			const boost::numeric::ublas::bounded_matrix<double,18,18>& Kloc_system,
			Matrix& rLeftHandSideMatrix
			)
	{
		KRATOS_TRY
				
		//calculate local rotation matrix
		aux33(0,0) = v1[0];  aux33(0,1) = v1[1];  aux33(0,2) = v1[2]; 
		aux33(1,0) = v2[0];  aux33(1,1) = v2[1];  aux33(1,2) = v2[2]; 
		aux33(2,0) = v3[0];  aux33(2,1) = v3[1];  aux33(2,2) = v3[2]; 
		
		//calculate the global rotation matrix ... VERY INEFFICIENT!!!!
		//this should be done block by block
		
		for(unsigned int I = 0; I<6; I++)
		{
			unsigned int index = I*3;
			for(unsigned int i = 0; i<3; i++)
			{
				for(unsigned int j=0; j<3; j++)
				{
					rot18(index+i,index+j) = aux33(i,j);
				}
				
			}
		}
		
		boost::numeric::ublas::bounded_matrix<double,18,18> temp = prod(Kloc_system,rot18);
		noalias(rLeftHandSideMatrix) = prod(trans(rot18),temp);
		
		KRATOS_CATCH("");
	}
	
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::RotateToGlobal(
			const array_1d<double,3>& v1,
			const array_1d<double,3>& v2,
			const array_1d<double,3>& v3,
			const boost::numeric::ublas::bounded_matrix<double,18,18>& Kloc_system,
			Matrix& rLeftHandSideMatrix,
			VectorType& rRightHandSideVector
		)
		{
			KRATOS_TRY
	
			//calculate local rotation matrix
			aux33(0,0) = v1[0];  aux33(0,1) = v1[1];  aux33(0,2) = v1[2]; 
			aux33(1,0) = v2[0];  aux33(1,1) = v2[1];  aux33(1,2) = v2[2]; 
			aux33(2,0) = v3[0];  aux33(2,1) = v3[1];  aux33(2,2) = v3[2]; 
	
			//calculate the global rotation matrix ... VERY INEFFICIENT!!!!
			//this should be done block by block
	
			for(unsigned int I = 0; I<6; I++)
			{
				unsigned int index = I*3;
				for(unsigned int i = 0; i<3; i++)
				{
					for(unsigned int j=0; j<3; j++)
					{
						rot18(index+i,index+j) = aux33(i,j);
					}
	
				}
			}
	
			boost::numeric::ublas::bounded_matrix<double,18,18> temp = prod(Kloc_system,rot18);
			noalias(rLeftHandSideMatrix) = prod(trans(rot18),temp);
	
			array_1d<double,18> aaa;
			noalias(aaa) = prod(trans(rot18),rRightHandSideVector);
			noalias(rRightHandSideVector) = aaa;
	
			KRATOS_CATCH("");
		}
					    
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::AddBodyForce(	
			    	const double& h,
				const double& Area,
				VectorType& rRightHandSideVector
			   )
	{
		KRATOS_TRY
				
		array_1d<double,3> bf;
		noalias(bf) = GetProperties()[BODY_FORCE];
		const double& density = GetProperties()[DENSITY];
		
		bf *= (density * h * 0.333333333333333333) * Area;
		
		rRightHandSideVector[0] += bf[0];
		rRightHandSideVector[1] += bf[1];
		rRightHandSideVector[2] += bf[2];
		
		rRightHandSideVector[6] += bf[0];
		rRightHandSideVector[7] += bf[1];
		rRightHandSideVector[8] += bf[2];		
		
		rRightHandSideVector[12] += bf[0];
		rRightHandSideVector[13] += bf[1];
		rRightHandSideVector[14] += bf[2];		
		
		KRATOS_CATCH("");
		
	}
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::NicePrint(const Matrix& A)
	{
		for(unsigned int i = 0; i<A.size1(); i++)
		{
			for(unsigned int j = 0; j<A.size2(); j++)
			{
				std::cout << A(i,j) << " " ;
			}
			std::cout << std::endl;
		}
	}
	
	
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& Output, const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY


		if(Output.size() != 1)
			Output.resize(1);

		if(rVariable==PK2_STRESS_TENSOR)
		{
			//calculate local coordinates and rotation matrix
			array_1d<double,3> v1;
			array_1d<double,3> v2;
			array_1d<double,3> v3;
			double x12, x23, x31, y12, y23, y31;
			double A;
	
			CalculateLocalGlobalTransformation( x12, x23, x31, y12, y23, y31,v1,v2,v3,A);
	
			double weight,loc1,loc2,loc3;		
			//double h = GetProperties()[THICKNESS];
			double h = 1.0; //note that we want the stress NOT the stress multiplied by the thickness
			
			//calculating material matrices
			CalculateMembraneElasticityTensor(mEm,h);

			//here we calculate everything in the local system of coordinates
			//calculate integration point 1
			weight = A;
			loc1 = 0.33333333333333;
			loc2 = 0.33333333333333;
			loc3 = 0.33333333333333;
			
			
			
			double beta0 = 1.5; //note that this is just for stress evaluation!
			CalculateMembraneB(mBm, beta0,  loc1, loc2, loc3, x12, x23, x31, y12, y23, y31 );	

			array_1d<double,9> local_disp;
			array_1d<double,3> local_stress, local_strain;
			array_1d<double,6> rotated_stress = ZeroVector(6);
			
			CalculatePureMembraneDisplacement(local_disp,v1,v2,v3);
			
			noalias(local_strain) = prod(trans(mBm),local_disp);
			noalias(local_stress) = prod(mEm,local_strain);
			
			//rotating the stress to global coordinates
			AddVoigtTensorComponents(local_stress[0],rotated_stress,v1,v1);
			AddVoigtTensorComponents(local_stress[1],rotated_stress,v2,v2);
			
			AddVoigtTensorComponents(local_stress[2],rotated_stress,v1,v2);
			AddVoigtTensorComponents(local_stress[2],rotated_stress,v2,v1);

			if(Output[0].size2() != 6)
				Output[0].resize(1,6,false);
			
			for(unsigned int ii = 0; ii<rotated_stress.size(); ii++)
				Output[0](0,ii) = rotated_stress[ii];

		}

		KRATOS_CATCH("")

	}
	
	
	void ShellIsotropic::CalculateOnIntegrationPoints(const Variable<double >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY


		if(Output.size() != 1)
			Output.resize(1);

		if(rVariable==TEMPERATURE)
		{
			//calculate local coordinates and rotation matrix
			array_1d<double,3> v1;
			array_1d<double,3> v2;
			array_1d<double,3> v3;
			double x12, x23, x31, y12, y23, y31;
			double A;
	
			CalculateLocalGlobalTransformation( x12, x23, x31, y12, y23, y31,v1,v2,v3,A);
	
			double weight,loc1,loc2,loc3;		
			//double h = GetProperties()[THICKNESS];
			double h = 1.0; //note that we want the stress NOT the stress multiplied by the thickness
			
			//calculating material matrices
			CalculateMembraneElasticityTensor(mEm,h);

			//here we calculate everything in the local system of coordinates
			//calculate integration point 1
			weight = A;
			loc1 = 0.33333333333333;
			loc2 = 0.33333333333333;
			loc3 = 0.33333333333333;
			
			
			
			double beta0 = 1.5; //note that this is just for stress evaluation!
			CalculateMembraneB(mBm, beta0,  loc1, loc2, loc3, x12, x23, x31, y12, y23, y31 );	

			array_1d<double,9> local_disp;
			array_1d<double,3> local_stress, local_strain, membrane_stress;
			array_1d<double,6> rotated_stress = ZeroVector(6);
			
			CalculatePureMembraneDisplacement(local_disp,v1,v2,v3);
			
			noalias(local_strain) = prod(trans(mBm),local_disp);
			noalias(local_stress) = prod(mEm,local_strain);

			noalias(membrane_stress) = local_stress;

			//bending part
//AHHHHHHHHHHHHH here i am changing the meaning
			h = GetProperties()[THICKNESS];
			array_1d<double,3> bending_stress;
			CalculateBendingElasticityTensor(mEb,h);
			CalculateBendingB(mBb, loc2, loc3, x12, x23, x31, y12, y23, y31 );
			CalculatePureBendingDisplacement(local_disp,v1,v2,v3);
			noalias(local_strain) = prod(trans(mBb),local_disp);
			noalias(bending_stress) = prod(mEb,local_strain); //TAKE CARE !! this should be changed ... added or removed
			bending_stress *= 6.00 / (h*h);



			noalias(local_stress) = membrane_stress + bending_stress;

			double von_mises_face1 = sqrt((local_stress[0] * local_stress[0] + local_stress[1] * local_stress[1] + pow(local_stress[0] - local_stress[1], 2.00) + 6 * local_stress[2] * local_stress[2])/2.00);


			noalias(local_stress) = membrane_stress - bending_stress;

			double von_mises_face2 = sqrt((local_stress[0] * local_stress[0] + local_stress[1] * local_stress[1] + pow(local_stress[0] - local_stress[1], 2.00) + 6 * local_stress[2] * local_stress[2])/2.00);
			
// 			//rotating the stress to global coordinates
// 			AddVoigtTensorComponents(local_stress[0],rotated_stress,v1,v1);
// 			AddVoigtTensorComponents(local_stress[1],rotated_stress,v2,v2);
			
// 			AddVoigtTensorComponents(local_stress[2],rotated_stress,v1,v2);
// 			AddVoigtTensorComponents(local_stress[2],rotated_stress,v2,v1);

			Output[0] = std::max(von_mises_face1,von_mises_face2);
			//Output[0] = von_mises;
			
//  			KRATOS_WATCH(von_mises);
			

		}

		KRATOS_CATCH("")

	}
	
	//***********************************************************************************
	//***********************************************************************************

	//auxiliary function needed in the calculation of output stresses
	inline void ShellIsotropic::AddVoigtTensorComponents(
			const double local_component,
			array_1d<double,6>& v,
			const array_1d<double,3>& a,
   			const array_1d<double,3>& b)
	{
				v[0] += local_component*a[0]*b[0];
				v[1] += local_component*a[1]*b[1];
				v[2] += local_component*a[2]*b[2];
				v[3] += local_component*a[0]*b[1];
				v[4] += local_component*a[1]*b[2];
				v[5] += local_component*a[0]*b[2];
				
	}
	
	
	
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::CalculateKg_GaussPointContribution(
			boost::numeric::ublas::bounded_matrix<double,18,18>& Kloc_system ,
			const boost::numeric::ublas::bounded_matrix<double,3,3>& Em,
			const double& weight,
			const double& h, /*thickness*/
			const double& loc1, /*local coords*/
			const double& loc2,
			const double& loc3,
			const double& x12,
			const double& x23,
			const double& x31,
			const double& y12,
			const double& y23,
			const double& y31,
			const array_1d<double,9>& membrane_disp
			)	
		{
		KRATOS_TRY
		//definition of auxiliaries ... should be put global to optimize
		array_1d<double,3> local_stress, local_strain;
		boost::numeric::ublas::bounded_matrix<double,2,2> stress_tensor, Jinv;
		boost::numeric::ublas::bounded_matrix<double,9,2> DN_Dx;
		boost::numeric::ublas::bounded_matrix<double,2,9> temp29;

		//membrane stresses and strains
 		double beta0 = CalculateBeta( Em );
 		CalculateMembraneB(mBm, beta0,  loc1, loc2, loc3, x12, x23, x31, y12, y23, y31 );	
		noalias(local_strain) = prod(trans(mBm),membrane_disp);
		noalias(local_stress) = prod(mEm,local_strain); //this stress is STRESS*thickness

		//write the equivalent stress tensor
		stress_tensor(0,0) = local_stress[0]; stress_tensor(0,1) = local_stress[2];
		stress_tensor(1,0) = local_stress[2]; stress_tensor(1,1) = local_stress[1];

		boost::numeric::ublas::bounded_matrix<double,2,9> DNu_loc,DNv_loc,DNw_loc,DN;

		//calculation of the local derivatives
		double alpha = 1.5;
		CalculateLocalShapeDerivatives(alpha,DNu_loc,DNv_loc,DNw_loc, loc1, loc2, loc3, x12, x23, x31, y12, y23, y31 );
	
		//calculate Jinv
		const double x21 = -x12;
		const double y21 = -y12;
		double aaa = (x21 * y31 - y21 * x31);
		Jinv(0,0) = y31 / aaa;
		Jinv(0,1) = -y21 / aaa;
		Jinv(1,0) = -x31 / aaa;
		Jinv(1,1) = x21 / aaa;

		//******************************************************************
		//assembling membrane part of the geometric stiffness

		//u-part ---> Kg_membrane = DNu * stress_tensor * trans(DNu);
		noalias(DN) = prod(Jinv,DNu_loc);
		noalias(temp29) = prod(stress_tensor,DN );
		noalias(mKloc99) = prod(trans(DN),temp29);

		//v-part ---> Kg_membrane += DNv * stress_tensor * trans(DNv);
		noalias(DN) = prod(Jinv,DNv_loc);
		noalias(temp29) = prod(stress_tensor,DN );
		noalias(mKloc99) += prod(trans(DN),temp29);

		AssembleMembraneContribution(   mKloc99, weight, Kloc_system  );

		//******************************************************************
		//w-part ---> Kg_bending += DNw * stress_tensor * trans(DNw);
		noalias(DN) = prod(Jinv,DNw_loc);
		noalias(temp29) = prod(stress_tensor,DN );
		noalias(mKloc99) = prod(trans(DN),temp29);

		AssembleBendingContribution(   mKloc99, weight, Kloc_system  );	

		KRATOS_CATCH("")
					
	}
				
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::CalculateAndAddKg(
				MatrixType& LHS, 
				boost::numeric::ublas::bounded_matrix<double,18,18>& rWorkMatrix,
				const double& x12,
				const double& x23,
				const double& x31,
				const double& y12,
				const double& y23,
				const double& y31,
				const array_1d<double,3>& v1,
				const array_1d<double,3>& v2,
				const array_1d<double,3>& v3,
				const double& A
			      )
	{
		KRATOS_TRY

		double weight,loc1,loc2,loc3;
		
		//calculate displacements in local coordinates
		array_1d<double,9> membrane_disp; //, bending_disp;
		
		CalculatePureMembraneDisplacement(membrane_disp,v1,v2,v3);
		noalias(rWorkMatrix) = ZeroMatrix(18,18);

		double h = GetProperties()[THICKNESS];

		//calculating material matrices
		CalculateMembraneElasticityTensor(mEm,h);

		//here we calculate everything in the local system of coordinates
// 		//calculate integration point 1
// 		weight = A/3;
// 		loc1= 0.5;
// 		loc2 = 0.5;
// 		loc3 = 0.0;
// 		CalculateKg_GaussPointContribution( Kg, mEm, mEb, weight, h, loc1, loc2, loc3, x12, x23, x31, y12, y23, y31,membrane_disp, bending_disp);
// 
// 		//calculate integration point 2
// 		weight = A/3;
// 		loc1= 0.0;
// 		loc2 = 0.5;
// 		loc3 = 0.5;
// 		CalculateKg_GaussPointContribution( Kg, mEm, mEb, weight, h, loc1, loc2, loc3, x12, x23, x31, y12, y23, y31,membrane_disp, bending_disp);
// 
// 		//calculate integration point 3
// 		weight = A/3;
// 		loc1= 0.5;
// 		loc2 = 0.0;
// 		loc3 = 0.5;
// 		CalculateKg_GaussPointContribution( Kg, mEm, mEb, weight, h, loc1, loc2, loc3, x12, x23, x31, y12, y23, y31,membrane_disp, bending_disp);
		//calculate integration point 1
		    double alpha1 =0.8168475730;
		    double beta1 = 0.0915762135;
		    double alpha2 = 0.1081030182;
		    double beta2 = 0.4459484909;
		    double gamma3 = 0.1099517437; // / 2.00;
		    double gamma4 = 0.2233815897; // / 2.00;
		weight = A*gamma3;
		loc1= alpha1;
		loc2 = beta1;
		loc3 = beta1;
		CalculateKg_GaussPointContribution( rWorkMatrix, mEm, weight, h, loc1, loc2, loc3, x12, x23, x31, y12, y23, y31,membrane_disp);

		//calculate integration point 2
		weight = A*gamma3;
		loc1= beta1;
		loc2 = alpha1;
		loc3 = beta1;
		CalculateKg_GaussPointContribution( rWorkMatrix, mEm, weight, h, loc1, loc2, loc3, x12, x23, x31, y12, y23, y31,membrane_disp);

		//calculate integration point 3
		weight = A*gamma3;
		loc1= beta1;
		loc2 = beta1;
		loc3 = alpha1;
		CalculateKg_GaussPointContribution( rWorkMatrix, mEm, weight, h, loc1, loc2, loc3, x12, x23, x31, y12, y23, y31,membrane_disp);

		//calculate integration point 4
		weight = A*gamma4;
		loc1= alpha2;
		loc2 = beta2;
		loc3 = beta2;
		CalculateKg_GaussPointContribution( rWorkMatrix, mEm, weight, h, loc1, loc2, loc3, x12, x23, x31, y12, y23, y31,membrane_disp);

		//calculate integration point 5
		weight = A*gamma4;
		loc1= beta2;
		loc2 = alpha2;
		loc3 = beta2;
		CalculateKg_GaussPointContribution( rWorkMatrix, mEm, weight, h, loc1, loc2, loc3, x12, x23, x31, y12, y23, y31,membrane_disp);

		//calculate integration point 4
		weight = A*gamma4;
		loc1= beta2;
		loc2 = beta2;
		loc3 = alpha2;
		CalculateKg_GaussPointContribution( rWorkMatrix, mEm, weight, h, loc1, loc2, loc3, x12, x23, x31, y12, y23, y31,membrane_disp);
		
		noalias(LHS) += rWorkMatrix;
		
		KRATOS_CATCH("")

	} 

	
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo)
	{
		if(rVariable==GEOMETRIC_STIFFNESS)
		{
			//resize Kg if needed
			if(Output.size1() != 18) Output.resize(18,18,false);
			boost::numeric::ublas::bounded_matrix<double,18,18> WorkMatrix;
			
			//calculate local coordinates and rotation matrix
			array_1d<double,3> v1;
			array_1d<double,3> v2;
			array_1d<double,3> v3;
			double x12, x23, x31, y12, y23, y31;
			double A;
			CalculateLocalGlobalTransformation( x12, x23, x31, y12, y23, y31,v1,v2,v3,A);
			
			//calculating the Kg in local coordinates
			noalias(Output) = ZeroMatrix(18,18);
			CalculateAndAddKg(Output,WorkMatrix,x12,x23,x31,y12,y23,y31,v1,v2,v3,A); 
			
			//performing rotation as needed
			noalias(WorkMatrix) = Output;
			RotateToGlobal(v1,v2,v3,WorkMatrix,Output);			
		}
	}




	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::CalculateLocalShapeDerivatives(
			double alpha,
			boost::numeric::ublas::bounded_matrix<double,2,9>& DNu_loc ,
			boost::numeric::ublas::bounded_matrix<double,2,9>& DNv_loc ,
			boost::numeric::ublas::bounded_matrix<double,2,9>& DNw_loc ,
			const double& a, /*local coords*/ //loc1
			const double& b, //loc2
			const double& c, //loc3
			const double& x12,
			const double& x23,
			const double& x31,
			const double& y12,
			const double& y23,
			const double& y31
			)	
			{
				KRATOS_TRY

				const double x21 = -x12;
				const double x32 = -x23;
				const double x13 = -x31;
				const double y21 = -y12;
				const double y32 = -y23;
				const double y13 = -y31;
				
				//***************************************************************************
				//DNu_Dloc
				DNu_loc(0,0) = -1;
				DNu_loc(0,1) = 0;
				DNu_loc(0,2) = -alpha * y12 * b + alpha * y31 * c / 0.2e1 + alpha * y12 / 0.2e1 - alpha * y12 * c / 0.2e1;
				DNu_loc(0,3) = 1;
				DNu_loc(0,4) = 0;
				DNu_loc(0,5) = alpha * y23 * c / 0.2e1 - alpha * y12 / 0.2e1 + alpha * y12 * b + alpha * y12 * c / 0.2e1;
				DNu_loc(0,6) = 0;
				DNu_loc(0,7) = 0;
				DNu_loc(0,8) = -alpha * c * (y31 + y23) / 0.2e1;;

				DNu_loc(1,0) = -1;
				DNu_loc(1,1) = 0;
				DNu_loc(1,2) = -alpha * y12 * b / 0.2e1 + alpha * y31 * c - alpha * y31 / 0.2e1 + alpha * y31 * b / 0.2e1;
				DNu_loc(1,3) = 0;
				DNu_loc(1,4) = 0;
				DNu_loc(1,5) = alpha * b * (y23 + y12) / 0.2e1;
				DNu_loc(1,6) = 1;
				DNu_loc(1,7) = 0;
				DNu_loc(1,8) = alpha * y31 / 0.2e1 - alpha * y31 * b / 0.2e1 - alpha * y31 * c - alpha * y23 * b / 0.2e1;

				//***************************************************************************
				//DNv_Dloc
				DNv_loc(0,0) = 0;
				DNv_loc(0,1) = -1;
				DNv_loc(0,2) = -alpha * x12 * b + alpha * x31 * c / 0.2e1 + alpha * x12 / 0.2e1 - alpha * x12 * c / 0.2e1;
				DNv_loc(0,3) = 0;
				DNv_loc(0,4) = 1;
				DNv_loc(0,5) = alpha * x23 * c / 0.2e1 - alpha * x12 / 0.2e1 + alpha * x12 * b + alpha * x12 * c / 0.2e1;
				DNv_loc(0,6) = 0;
				DNv_loc(0,7) = 0;
				DNv_loc(0,8) = -alpha * c * (x31 + x23) / 0.2e1;

				DNv_loc(1,0) = 0;
				DNv_loc(1,1) = -1;
				DNv_loc(1,2) = -alpha * x12 * b / 0.2e1 + alpha * x31 * c - alpha * x31 / 0.2e1 + alpha * x31 * b / 0.2e1;
				DNv_loc(1,3) = 0;
				DNv_loc(1,4) = 0;
				DNv_loc(1,5) = alpha * b * (x23 + x12) / 0.2e1;
				DNv_loc(1,6) = 0;
				DNv_loc(1,7) = 1;
				DNv_loc(1,8) = alpha * x31 / 0.2e1 - alpha * x31 * b / 0.2e1 - alpha * x31 * c - alpha * x23 * b / 0.2e1;




				//***************************************************************************
				//Dw 
				DNw_loc(0,0) = -0.6e1 * b - 0.4e1 * c + 0.6e1 * b * b + 0.8e1 * b * c + 0.4e1 * c * c;
				DNw_loc(0,1) = 0.4e1 * y12 * b + 0.1500000000e1 * y13 * c - 0.3e1 * y12 * b * b - 0.1e1 * b * y13 * c - 0.3e1 * c * y12 * b - 0.1500000000e1 * y13 * c * c - 0.1e1 * y12 + 0.1500000000e1 * y12 * c - 0.5000000000e0 * y12 * c * c;
				DNw_loc(0,2) = -0.4e1 * x12 * b - 0.1500000000e1 * x13 * c + 0.3e1 * x12 * b * b + b * x13 * c + 0.3e1 * c * x12 * b + 0.1500000000e1 * x13 * c * c + x12 - 0.1500000000e1 * x12 * c + 0.5000000000e0 * x12 * c * c;
				DNw_loc(0,3) = 0.6e1 * b - 0.6e1 * b * b - 0.4e1 * b * c + 0.2e1 * c - 0.2e1 * c * c;
				DNw_loc(0,4) = -0.1e1 * b * c * y23 - 0.2e1 * y21 * b + 0.3e1 * b * b * y21 + 0.3e1 * b * c * y21 - 0.5000000000e0 * y23 * c + 0.5000000000e0 * y23 * c * c - 0.5000000000e0 * y21 * c + 0.5000000000e0 * y21 * c * c;
				DNw_loc(0,5) = b * c * x23 + 0.2e1 * x21 * b - 0.3e1 * b * b * x21 - 0.3e1 * b * c * x21 + 0.5000000000e0 * x23 * c - 0.5000000000e0 * x23 * c * c + 0.5000000000e0 * x21 * c - 0.5000000000e0 * x21 * c * c;
				DNw_loc(0,6) = -0.4e1 * b * c + 0.2e1 * c - 0.2e1 * c * c;
				DNw_loc(0,7) = 0.1500000000e1 * c * c * y31 - 0.5000000000e0 * c * c * y32 + b * c * y31 + b * c * y32 - 0.5000000000e0 * y31 * c - 0.5000000000e0 * y32 * c;
				DNw_loc(0,8) = -0.1500000000e1 * c * c * x31 + 0.5000000000e0 * c * c * x32 - 0.1e1 * b * c * x31 - 0.1e1 * b * c * x32 + 0.5000000000e0 * x31 * c + 0.5000000000e0 * x32 * c;

				DNw_loc(1,0) = -0.4e1 * b - 0.6e1 * c + 0.4e1 * b * b + 0.8e1 * b * c + 0.6e1 * c * c;
				DNw_loc(1,1) = 0.1500000000e1 * y12 * b + 0.4e1 * y13 * c - 0.1500000000e1 * y12 * b * b - 0.3e1 * b * y13 * c - 0.1e1 * c * y12 * b - 0.3e1 * y13 * c * c - 0.1e1 * y13 + 0.1500000000e1 * y13 * b - 0.5000000000e0 * y13 * b * b;
				DNw_loc(1,2) = -0.1500000000e1 * x12 * b - 0.4e1 * x13 * c + 0.1500000000e1 * x12 * b * b + 0.3e1 * b * x13 * c + c * x12 * b + 0.3e1 * x13 * c * c + x13 - 0.1500000000e1 * x13 * b + 0.5000000000e0 * x13 * b * b;
				DNw_loc(1,3) = -0.4e1 * b * c + 0.2e1 * b - 0.2e1 * b * b;
				DNw_loc(1,4) = -0.5000000000e0 * b * b * y23 + 0.1500000000e1 * b * b * y21 + b * c * y23 + b * c * y21 - 0.5000000000e0 * y23 * b - 0.5000000000e0 * y21 * b;
				DNw_loc(1,5) = 0.5000000000e0 * b * b * x23 - 0.1500000000e1 * b * b * x21 - 0.1e1 * b * c * x23 - 0.1e1 * b * c * x21 + 0.5000000000e0 * x23 * b + 0.5000000000e0 * x21 * b;
				DNw_loc(1,6) = 0.6e1 * c - 0.6e1 * c * c - 0.4e1 * b * c + 0.2e1 * b - 0.2e1 * b * b;
				DNw_loc(1,7) = -0.2e1 * y31 * c + 0.3e1 * b * c * y31 + 0.3e1 * c * c * y31 - 0.1e1 * b * c * y32 - 0.5000000000e0 * y31 * b + 0.5000000000e0 * y31 * b * b - 0.5000000000e0 * y32 * b + 0.5000000000e0 * y32 * b * b;
				DNw_loc(1,8) = 0.2e1 * x31 * c - 0.3e1 * b * c * x31 - 0.3e1 * c * c * x31 + b * c * x32 + 0.5000000000e0 * x31 * b - 0.5000000000e0 * x31 * b * b + 0.5000000000e0 * x32 * b - 0.5000000000e0 * x32 * b * b;

				KRATOS_CATCH("")
			}
			
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::CalculateProjectionOperator(
			boost::numeric::ublas::bounded_matrix<double,18,18>& rProjOperator, 
			const double& x12,
			const double& x23,
			const double& x31,
			const double& y12,
			const double& y23,
			const double& y31
			)
			{
				KRATOS_TRY
						
				boost::numeric::ublas::bounded_matrix<double,18,3> psi = ZeroMatrix(18,3);
				boost::numeric::ublas::bounded_matrix<double,18,3> rho = ZeroMatrix(18,3);
				
				double x21 = -x12;
								
				//ATTENTION THIS MAY GIVE PROBLEMS FOR THE ANISOTROPIC CASE AS 
				//THE FIRST BASE VECTOR IS NOT PARALLEL TO THE SIDE 21
				//filling psi
				psi(3,0) = 1.0; psi(4,1) = 1.0; psi(5,2) = 1.0; //second block
				psi(6,2) = 0.0; psi(7,2) = x21; //third block
				psi(8,0) = 0.0; psi(8,1) = -x21; 	
				psi(9,0) = 1.0; psi(10,1) = 1.0; psi(11,2) = 1.0; //fourth block
				psi(12,2) = -y31; psi(13,2) = x31; //fifth block
				psi(14,0) = y31; psi(14,1) = -x31; 	
				psi(15,0) = 1.0; psi(16,1) = 1.0; psi(17,2) = 1.0; //sixt block
			
				//filling rho
				rho(1,2) = -1.0/x21; rho(2,0) = (x31-x21)/(y31*x21); rho(2,1) = (1.0)/(x21);
				rho(7,2) =  1.0/x21; rho(8,0) = (-x31)/(y31*x21);    rho(8,1) =-(1.0)/(x21);
				rho(14,0) = 1.0/y31;
				
				//completing the calculation of the projections
				noalias(rProjOperator) = IdentityMatrix(18,18);
				//noalias(rProjOperator) -= prod(trans(psi),rho);
				noalias(rProjOperator) -= prod(psi,trans(rho));
				
				
				KRATOS_CATCH("")
			}
			
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::ApplyProjection(
		boost::numeric::ublas::bounded_matrix<double,18,18>& rLeftHandSideMatrix, 
  		VectorType& rRightHandSideVector,
    		boost::numeric::ublas::bounded_matrix<double,18,18>& rWorkMatrix, 
		array_1d<double,18>& rWorkArray,
  		const boost::numeric::ublas::bounded_matrix<double,18,18>& rProjOperator
		)
		{
			KRATOS_TRY
					
			//LHS = Ptrans * LHS * P
			noalias(rWorkMatrix) = prod(rLeftHandSideMatrix,rProjOperator);
			noalias(rLeftHandSideMatrix) = prod(trans(rProjOperator),rWorkMatrix);
			
			//RHS = Ptrans * RHS
			noalias(rWorkArray) = prod(trans(rProjOperator),rRightHandSideVector);
			noalias(rRightHandSideVector) = rWorkArray;
						
			KRATOS_CATCH("")
		}
			
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::UpdateNodalReferenceSystem(
		const double& x12,
		const double& x23,
		const double& x31,
		const double& y12,
		const double& y23,
		const double& y31
		)
		{
			KRATOS_TRY
					
			boost::numeric::ublas::bounded_matrix<double,3,3> Ttilde;
			boost::numeric::ublas::bounded_matrix<double,3,3> Omega = ZeroMatrix(3,3);
			array_1d<double,3> rot; //this will contain the incremental rotation in the last iteration
					
			//loop over nodes
			for(unsigned int I=0; I<3; I++)
			{
				//calculating the incremental rotation during the last iteration
				noalias(rot) = GetGeometry()[I].FastGetSolutionStepValue(ROTATION);
				noalias(rot) -= rot_oldit[I];
				
				//saving the value of the nodal rotation at the last iteration
				noalias(rot_oldit[I]) = GetGeometry()[I].FastGetSolutionStepValue(ROTATION);
				
				double omega_scalar_2 = rot[0]*rot[0] + rot[1]*rot[1] + rot[2]*rot[2];
				
				Omega(0,0) = 0.0;     Omega(0,1) = -rot[2]; 	Omega(0,2)=rot[1];
				Omega(1,0) =  rot[2]; Omega(1,1) = 0.0;		Omega(1,2)=-rot[0];
				Omega(2,0) = -rot[1]; Omega(2,1)=rot[0];	Omega(2,2) = 0.0;	
				
				double temp = 1.0/(1.0 + 0.25*omega_scalar_2);
				
				noalias(Ttilde) = IdentityMatrix(3,3);
				noalias(Ttilde) += temp * Omega;
				noalias(Ttilde) += 0.5*temp * prod( Omega, Omega);
				
				//attention here i overwrite Omega with something different
				//to do efficiently Tsnew = Ttilde*Ts
				noalias(Omega) = prod(Ttilde,mTs[I]);
				noalias(mTs[I]) = Omega;
			}
				
				
			KRATOS_CATCH("")
		}
		
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::SaveOriginalReference(
		const array_1d<double,3>& v1,
		const array_1d<double,3>& v2,
		const array_1d<double,3>& v3
		)
		{
			KRATOS_TRY
					
			//calculating original reference matrix
			// if vE is expressed in the original system of local coords
			// mTE0*vE gives its expression in global coordinates
			for(unsigned int i=0; i<3;i++)
			{
				mTE0(i,0) = v1[i];
				mTE0(i,1) = v2[i];
				mTE0(i,2) = v3[i];
			}
			
			//initializing nodal triad matrices
			noalias(mTs[0]) = IdentityMatrix(3,3);
			noalias(mTs[1]) = IdentityMatrix(3,3);
			noalias(mTs[2]) = IdentityMatrix(3,3);
			
		
			KRATOS_CATCH("")
		}
			
	//************************************************************************************
	//************************************************************************************
	//this returns "pure" displacements and rotations
		void ShellIsotropic::CalculatePureDisplacement(
				Vector& values,
				const array_1d<double,3>& v1,
				const array_1d<double,3>& v2,
				const array_1d<double,3>& v3
				)
		{
			const unsigned int number_of_nodes = 3;
			unsigned int MatSize = 18;
			if(values.size() != MatSize)	values.resize(MatSize,false);
			
			array_1d<double,3> temp_vec, temp_vec1;
			boost::numeric::ublas::bounded_matrix<double,3,3> Ttilde, temp, Omega, TE, aux;
			
			//calculating original reference matrix
			for(unsigned int i=0; i<3;i++)
			{
				TE(i,0) = v1[i];
				TE(i,1) = v2[i];
				TE(i,2) = v3[i];
			}
			
			double aaa;
			
			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				//calculate matrix T = TEtrans * Ts * TE0
				noalias(temp) = prod(mTs[i],mTE0);
				noalias(Ttilde) = prod(trans(TE),temp);
				
				//calculate Omega = 2.0*(T-I)*(T+I)^-1
				noalias(aux) = IdentityMatrix(3,3);
				noalias(aux) += Ttilde;
				InvertMatrix(aux,temp,aaa); //now temp contains the inverse
				noalias(aux) = Ttilde;
				noalias(aux) -= IdentityMatrix(3,3);
				noalias(Omega) = 2.0 * prod(aux,temp);
				
				//extract pure rotations from Omega
				temp_vec[0] = GetGeometry()[i].X() - GetGeometry()[0].X();
				temp_vec[1] = GetGeometry()[i].Y() - GetGeometry()[0].Y();
				temp_vec[2] = GetGeometry()[i].Z() - GetGeometry()[0].Z();
				noalias(temp_vec1) = prod(trans(TE),temp_vec);
				
				//this are the node position in the original reference system
				temp_vec[0] = GetGeometry()[i].X0() - GetGeometry()[0].X0();
				temp_vec[1] = GetGeometry()[i].Y0() - GetGeometry()[0].Y0();
				temp_vec[2] = GetGeometry()[i].Z0() - GetGeometry()[0].Z0();
				//noalias(temp_vec1) -= prod(mTE0,temp_vec); //ojo ... sin transpuesta!
				noalias(temp_vec1) -= prod(trans(mTE0),temp_vec); 
								
				unsigned int index = i*6;
				values[index] 	  = temp_vec1[0];
				values[index + 1] = temp_vec1[1];
				values[index + 2] = 0.0; //temp_vec1[2]; //the v3 disp has to be 0!

				//extract pure rotations from Omega
				values[index + 3] = -Omega(1,2);
				values[index + 4] = Omega(0,2);
				values[index + 5] = -Omega(0,1);
			}
		}
		
	//************************************************************************************
	//************************************************************************************
	//this returns "pure" displacements and rotations
	void ShellIsotropic::CalculatePureMembraneDisplacement(
		array_1d<double,9>& values,
		const array_1d<double,3>& v1,
		const array_1d<double,3>& v2,
		const array_1d<double,3>& v3
		)
		{
		const unsigned int number_of_nodes = 3;
		//const unsigned int dim = 3;
		unsigned int MatSize = 18;
		if(values.size() != MatSize)	values.resize(MatSize,false);

		array_1d<double,3> temp_vec, temp_vec1;
		boost::numeric::ublas::bounded_matrix<double,3,3> Ttilde, temp, Omega, TE, aux;

		//calculating original reference matrix
		for(unsigned int i=0; i<3;i++)
		{
			TE(i,0) = v1[i];
			TE(i,1) = v2[i];
			TE(i,2) = v3[i];
		}

		double aaa;

		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			//calculate matrix T = TEtrans * Ts * TE0
			noalias(temp) = prod(mTs[i],mTE0);
			noalias(Ttilde) = prod(trans(TE),temp);
	
			//calculate Omega = 2.0*(T-I)*(T+I)^-1
			noalias(aux) = IdentityMatrix(3,3);
			noalias(aux) += Ttilde;
			InvertMatrix(aux,temp,aaa); //now temp contains the inverse
			noalias(aux) = Ttilde;
			noalias(aux) -= IdentityMatrix(3,3);
			noalias(Omega) = 2.0 * prod(aux,temp);
	
			//node pos in the current config
			temp_vec[0] = GetGeometry()[i].X() - GetGeometry()[0].X();
			temp_vec[1] = GetGeometry()[i].Y() - GetGeometry()[0].Y();
			temp_vec[2] = GetGeometry()[i].Z() - GetGeometry()[0].Z();
			noalias(temp_vec1) = prod(trans(TE),temp_vec);
	
			//this are the node position in the original reference system
			temp_vec[0] = GetGeometry()[i].X0() - GetGeometry()[0].X0();
			temp_vec[1] = GetGeometry()[i].Y0() - GetGeometry()[0].Y0();
			temp_vec[2] = GetGeometry()[i].Z0() - GetGeometry()[0].Z0();
			
			//finalizing calculation of pure nodal displacements
			noalias(temp_vec1) -= prod(trans(mTE0),temp_vec); 
					
			unsigned int index = i*3;
			values[index] 	  = temp_vec1[0];
			values[index + 1] = temp_vec1[1];
			values[index + 2] = -Omega(0,1); //this contains the rot in dir V3
		}
	}
	
	//************************************************************************************
	//************************************************************************************
	//this returns "pure" displacements and rotations
	void ShellIsotropic::CalculatePureBendingDisplacement(
		array_1d<double,9>& values,
		const array_1d<double,3>& v1,
		const array_1d<double,3>& v2,
		const array_1d<double,3>& v3
		)
		{
		const unsigned int number_of_nodes = 3;
		//const unsigned int dim = 3;
		unsigned int MatSize = 18;
		if(values.size() != MatSize)	values.resize(MatSize,false);

		array_1d<double,3> temp_vec, temp_vec1;
		boost::numeric::ublas::bounded_matrix<double,3,3> Ttilde, temp, Omega, TE, aux;

		//calculating original reference matrix
		for(unsigned int i=0; i<3;i++)
		{
			TE(i,0) = v1[i];
			TE(i,1) = v2[i];
			TE(i,2) = v3[i];
		}

		double aaa;

		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			//calculate matrix T = TEtrans * Ts * TE0
			noalias(temp) = prod(mTs[i],mTE0);
			noalias(Ttilde) = prod(trans(TE),temp);
	
			//calculate Omega = 2.0*(T-I)*(T+I)^-1
			noalias(aux) = IdentityMatrix(3,3);
			noalias(aux) += Ttilde;
			InvertMatrix(aux,temp,aaa); //now temp contains the inverse
			noalias(aux) = Ttilde;
			noalias(aux) -= IdentityMatrix(3,3);
			noalias(Omega) = 2.0 * prod(aux,temp);
	
			//node pos in the current config
			temp_vec[0] = GetGeometry()[i].X() - GetGeometry()[0].X();
			temp_vec[1] = GetGeometry()[i].Y() - GetGeometry()[0].Y();
			temp_vec[2] = GetGeometry()[i].Z() - GetGeometry()[0].Z();
			noalias(temp_vec1) = prod(trans(TE),temp_vec);
	
			//this are the node position in the original reference system
			temp_vec[0] = GetGeometry()[i].X0() - GetGeometry()[0].X0();
			temp_vec[1] = GetGeometry()[i].Y0() - GetGeometry()[0].Y0();
			temp_vec[2] = GetGeometry()[i].Z0() - GetGeometry()[0].Z0();
			
			//finalizing calculation of pure nodal displacements
			noalias(temp_vec1) -= prod(trans(mTE0),temp_vec); 
					
			unsigned int index = i*3;

			values[index] 	  = temp_vec1[2];
			values[index + 1 ] = -Omega(1,2);
			values[index + 2] = Omega(0,2);
		}
	}
			
	//************************************************************************************
	//************************************************************************************
		void ShellIsotropic::InvertMatrix(const boost::numeric::ublas::bounded_matrix<double,3,3>& InputMatrix, 
					  boost::numeric::ublas::bounded_matrix<double,3,3>& InvertedMatrix,
				         double& InputMatrixDet) 
		{
			KRATOS_TRY
					if(InvertedMatrix.size1() != 3 || InvertedMatrix.size2() != 3)
					InvertedMatrix.resize(3,3);
			
			//filling the inverted matrix with the algebraic complements
			//first column
			InvertedMatrix(0,0) = InputMatrix(1,1)*InputMatrix(2,2) - InputMatrix(1,2)*InputMatrix(2,1);
			InvertedMatrix(1,0) = -InputMatrix(1,0)*InputMatrix(2,2) + InputMatrix(1,2)*InputMatrix(2,0);
			InvertedMatrix(2,0) = InputMatrix(1,0)*InputMatrix(2,1) - InputMatrix(1,1)*InputMatrix(2,0);
			
			//second column
			InvertedMatrix(0,1) = -InputMatrix(0,1)*InputMatrix(2,2) + InputMatrix(0,2)*InputMatrix(2,1);
			InvertedMatrix(1,1) = InputMatrix(0,0)*InputMatrix(2,2) - InputMatrix(0,2)*InputMatrix(2,0);
			InvertedMatrix(2,1) = -InputMatrix(0,0)*InputMatrix(2,1) + InputMatrix(0,1)*InputMatrix(2,0);
			
			//third column
			InvertedMatrix(0,2) = InputMatrix(0,1)*InputMatrix(1,2) - InputMatrix(0,2)*InputMatrix(1,1);
			InvertedMatrix(1,2) = -InputMatrix(0,0)*InputMatrix(1,2) + InputMatrix(0,2)*InputMatrix(1,0);
			InvertedMatrix(2,2) = InputMatrix(0,0)*InputMatrix(1,1) - InputMatrix(0,1)*InputMatrix(1,0);
			
			//calculation of determinant (of the input matrix)
			InputMatrixDet = InputMatrix(0,0)*InvertedMatrix(0,0) + InputMatrix(0,1)*InvertedMatrix(1,0) + InputMatrix(0,2)*InvertedMatrix(2,0);
			
			//finalizing the calculation of the inverted matrix
			InvertedMatrix /= InputMatrixDet;
			KRATOS_CATCH("")
		}		
		
		
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::Initialize()
		{
			KRATOS_TRY
			//calculate local coordinates and rotation matrix
			array_1d<double,3> v1;
			array_1d<double,3> v2;
			array_1d<double,3> v3;
			double x12, x23, x31, y12, y23, y31;
			double A;
	
			CalculateLocalGlobalTransformation( x12, x23, x31, y12, y23, y31,v1,v2,v3,A);
		
			SaveOriginalReference(v1,v2,v3);
			
			noalias(rot_oldit[0]) = GetGeometry()[0].FastGetSolutionStepValue(ROTATION);
			noalias(rot_oldit[1]) = GetGeometry()[1].FastGetSolutionStepValue(ROTATION);
			noalias(rot_oldit[2]) = GetGeometry()[2].FastGetSolutionStepValue(ROTATION);
			
			KRATOS_CATCH("")
		}
		
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
			
		//calculate local coordinates and rotation matrix
		array_1d<double,3> v1;
		array_1d<double,3> v2;
		array_1d<double,3> v3;
		double x12, x23, x31, y12, y23, y31;
		double A;
	
		CalculateLocalGlobalTransformation( x12, x23, x31, y12, y23, y31,v1,v2,v3,A);
				
		UpdateNodalReferenceSystem(x12,x23,x31,y12,y23,y31);

		KRATOS_CATCH("")

	}


	

	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		//calculate local coordinates and rotation matrix
		array_1d<double,3> v1;
		array_1d<double,3> v2;
		array_1d<double,3> v3;
		double x12, x23, x31, y12, y23, y31;
		double area;
	
		CalculateLocalGlobalTransformation( x12, x23, x31, y12, y23, y31,v1,v2,v3,area);

		double h = GetProperties()[THICKNESS];
		double density = GetProperties()[DENSITY];
		double node_mass = area * density * h / 3.00;
	
		//lumped
		unsigned int MatSize = 18;
		if(rMassMatrix.size1() != MatSize)
			rMassMatrix.resize(MatSize,MatSize,false);

		rMassMatrix = ZeroMatrix(MatSize,MatSize);

		rMassMatrix(0,0) = node_mass;
		rMassMatrix(1,1) = node_mass;
		rMassMatrix(2,2) = node_mass;

		rMassMatrix(6,6) = node_mass;
		rMassMatrix(7,7) = node_mass;
		rMassMatrix(8,8) = node_mass;

		rMassMatrix(12,12) = node_mass;
		rMassMatrix(13,13) = node_mass;
		rMassMatrix(14,14) = node_mass;

		KRATOS_CATCH("")
	}

		//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::GetFirstDerivativesVector(Vector& values, int Step)
	{
		unsigned int MatSize = 18;
		if(values.size() != MatSize)   values.resize(MatSize,false);
		for (unsigned int i=0;i<3;i++)
		{
			unsigned int index = i*6;
			values[index] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_X,Step);
			values[index + 1] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_Y,Step);
			values[index + 2] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_Z,Step);
			values[index + 3] = 0.0;
			values[index + 4] = 0.0;
			values[index + 5] = 0.0;
		}
	}
	//************************************************************************************
	//************************************************************************************
	void ShellIsotropic::GetSecondDerivativesVector(Vector& values, int Step)
	{
		unsigned int MatSize = 18;
		if(values.size() != MatSize) values.resize(MatSize,false);
		for (unsigned int i=0;i<3;i++)
		{
			unsigned int index = i*6;
			values[index] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_X,Step);
			values[index + 1] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Y,Step);
			values[index + 2] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Z,Step);
			values[index + 3] = 0.0;
			values[index + 4] = 0.0;
			values[index + 5] = 0.0;
		}
	}

}



