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
//   Last modified by:    $Author: anonymous $
//   Date:                $Date: 2007-12-12 08:24:59 $
//   Revision:            $Revision: 1.3 $
//
//


// System includes 
 
// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/shell_anisotropic_linear.h"
#include "structural_application.h"
#include "utilities/math_utils.h" 

                                                                                                                              
namespace Kratos
{
	namespace ShellAnisotropicLinearAuxiliaries
	{
		array_1d<int,9> local_indices;
		array_1d<int,9> local_j;
		boost::numeric::ublas::bounded_matrix<double,9,3> mBm = ZeroMatrix(9,3); //Membrane displacement-strain matrix
		boost::numeric::ublas::bounded_matrix<double,9,3> mBb = ZeroMatrix(9,3); //bending displacement-strain matrix

		boost::numeric::ublas::bounded_matrix<double,3,3> Q = ZeroMatrix(3,3); 
		boost::numeric::ublas::bounded_matrix<double,3,3> Q1 = ZeroMatrix(3,3); 
		boost::numeric::ublas::bounded_matrix<double,3,3> Q2 = ZeroMatrix(3,3); 
		boost::numeric::ublas::bounded_matrix<double,3,3> Q3 = ZeroMatrix(3,3); 
		boost::numeric::ublas::bounded_matrix<double,3,3> aux33 = ZeroMatrix(3,3);
		boost::numeric::ublas::bounded_matrix<double,3,3> Te = ZeroMatrix(3,3);

		boost::numeric::ublas::bounded_matrix<double,3,3> msA = ZeroMatrix(3,3);
		boost::numeric::ublas::bounded_matrix<double,3,3> msB = ZeroMatrix(3,3);
		boost::numeric::ublas::bounded_matrix<double,3,3> msD = ZeroMatrix(3,3);

		array_1d<double,9> H1 = ZeroVector(9); 
		array_1d<double,9> H2 = ZeroVector(9); 
		array_1d<double,9> H3 = ZeroVector(9);
		array_1d<double,9> H4 = ZeroVector(9);

		boost::numeric::ublas::bounded_matrix<double,9,3> TTu  = ZeroMatrix(9,3); //attention 9*3 and not 3*9
		boost::numeric::ublas::bounded_matrix<double,3,9> aux39  = ZeroMatrix(3,9);
		
		boost::numeric::ublas::bounded_matrix<double,18,18> mKloc_system  = ZeroMatrix(18,18); //stiffness matrix in the local reference system
		boost::numeric::ublas::bounded_matrix<double,18,18> rot18  = ZeroMatrix(18,18); //TAKE CARE!! this is VERY inefficient

		boost::numeric::ublas::bounded_matrix<double,9,9> mKloc99 = ZeroMatrix(9,9);
		Vector values = ZeroVector(18); //help vector to get the current values of displacements (and rotations)
	}
	
	using namespace ShellAnisotropicLinearAuxiliaries;
	
	//************************************************************************************
	//************************************************************************************
	ShellAnisotropicLinear::ShellAnisotropicLinear(IndexType NewId, GeometryType::Pointer pGeometry)
	: Element(NewId, pGeometry) 
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	ShellAnisotropicLinear::ShellAnisotropicLinear(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
	}

	Element::Pointer ShellAnisotropicLinear::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new ShellAnisotropicLinear(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	ShellAnisotropicLinear::~ShellAnisotropicLinear()
	{
	}
	
	//************************************************************************************
	//************************************************************************************
	void ShellAnisotropicLinear::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
	
	CalculateAllMatrices(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
	}	
	
	
	

	//************************************************************************************
	//************************************************************************************
	//getting the local coordinates and calculating the base vectors
	void ShellAnisotropicLinear::CalculateLocalGlobalTransformation(
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
		//getting the composite direction and aligning the base with that
		noalias(v1) = GetValue(COMPOSITE_DIRECTION);
				
		//vector aligned with the side 1-2
		array_1d<double,3> v21;
		v21[0] = GetGeometry()[1].X()-GetGeometry()[0].X(); 
		v21[1] = GetGeometry()[1].Y()-GetGeometry()[0].Y();
		v21[2] = GetGeometry()[1].Z()-GetGeometry()[0].Z();

		//vector aligned with the side 1-3
		array_1d<double,3> v31;
		v31[0]=GetGeometry()[2].X()-GetGeometry()[0].X();
		v31[1]=GetGeometry()[2].Y()-GetGeometry()[0].Y();
		v31[2]=GetGeometry()[2].Z()-GetGeometry()[0].Z();
		
		//forming the (area) normal
		MathUtils<double>::CrossProduct(v3,v21,v31);
		area = 0.5 * norm_2(v3);

		//normalizing normal vector
		v3 /= (2.0*area);
		
		//projecting the composite direction to the plane of the element and re-normalizing
		double aaa = inner_prod(v3,v1);
		noalias(v1) -= aaa*v3;
		double norm_v1 = norm_2(v1);
		v1 /= norm_v1;

		//forming the "second" base vector - it is already normalized
		//MathUtils<double>::CrossProduct(v2,v3,v21);
		MathUtils<double>::CrossProduct(v2,v3,v1);
		
		//calculating the local coordinates
		double x21 	= inner_prod(v21,v1);
		double y21 	= inner_prod(v21,v2);
		x31 		= inner_prod(v31,v1);
		y31 		= inner_prod(v31,v2);

		x23 = x21-x31;
		y23 = y21-y31;
		
		x12=-x21;
		y12=-y21;

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	//getting the local coordinates and calculating the base vectors
	/*void ShellAnisotropicLinear::CalculateLocalGlobalTransformation(
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
						v1[0] = GetGeometry()[1].X()-
GetGeometry()[0].X(); 
				v1[1] = GetGeometry()[1].Y()-GetGeometry()[0].Y();
				v1[2] = GetGeometry()[1].Z()-GetGeometry()[0].Z();
				double x21=norm_2(v1); //here we assume the basis 
aligned with 1-2
				double y21 = 0.0; 
				x12=-x21;
				y12=-y21;

		//forming the (area) normal
array_1d<double,3> temp;
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
	*/
			
	//************************************************************************************
	//************************************************************************************
	void ShellAnisotropicLinear::CalculateMembraneB( 
boost::numeric::ublas::bounded_matrix<double,9,3>& Bm,
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
	void ShellAnisotropicLinear::CalculateBendingB( 
boost::numeric::ublas::bounded_matrix<double,9,3>& Bb,
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
	void ShellAnisotropicLinear::CalculateMembraneContribution( 
					const boost::numeric::ublas::bounded_matrix<double,9,3>& Bm,
					const boost::numeric::ublas::bounded_matrix<double,3,3>& A,
					boost::numeric::ublas::bounded_matrix<double,9,9>& Km)
	{
		KRATOS_TRY
		noalias(aux39) = prod(A,trans(Bm) );
		noalias(Km) = prod(Bm, aux39);
		
		KRATOS_CATCH("")
	}
	
	//************************************************************************************
	//************************************************************************************
	void ShellAnisotropicLinear::AssembleMembraneContribution( 
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
	void ShellAnisotropicLinear::CalculateBendingContribution(
					const boost::numeric::ublas::bounded_matrix<double,9,3>& Bb, 
					const boost::numeric::ublas::bounded_matrix<double,3,3>& D, 
       				boost::numeric::ublas::bounded_matrix<double,9,9>& Kb
					 )
	{
		KRATOS_TRY
		noalias(aux39) = prod(D,trans(Bb) );
		noalias(Kb) = prod(Bb, aux39);
		KRATOS_CATCH("")
	}
	
	//************************************************************************************
	//************************************************************************************
	void ShellAnisotropicLinear::AssembleBendingContribution( 
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
	void ShellAnisotropicLinear::CalculateMixedContribution(
			const boost::numeric::ublas::bounded_matrix<double,9,3>& Bm, 
			const boost::numeric::ublas::bounded_matrix<double,9,3>& Bb, 
   			const boost::numeric::ublas::bounded_matrix<double,3,3>& B, 
   			boost::numeric::ublas::bounded_matrix<double,9,9>& Kmix_top
			)
	{
		KRATOS_TRY
		noalias(aux39) = prod(B,trans(Bb) );
		noalias(Kmix_top) = prod(Bm, aux39);
		KRATOS_CATCH("")
	}
			
	//************************************************************************************
	//************************************************************************************
	void ShellAnisotropicLinear::AssembleMixedContribution( 
			const boost::numeric::ublas::bounded_matrix<double,9,9>& Kmix_top, 
   			const double& coeff, 
 			boost::numeric::ublas::bounded_matrix<double,18,18>& Kloc_system 
			)
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
					
				local_j[0] = 2;
				local_j[1] = 3;
				local_j[2] = 4;
				local_j[3] = 8;
				local_j[4] = 9;
				local_j[5] = 10;
				local_j[6] = 14;
				local_j[7] = 15;
				local_j[8] = 16;
	
				for(unsigned int i = 0; i<9; i++)
				{
					for(unsigned int j = 0; j<9; j++)
					{
						Kloc_system( local_indices[i],local_j[j] ) += coeff*Kmix_top(i,j);
						Kloc_system( local_j[j],local_indices[i] ) += coeff*Kmix_top(i,j);
					}
				}
				KRATOS_CATCH("")
			}
			
			
	//************************************************************************************
	//************************************************************************************
	void ShellAnisotropicLinear::CalculateGaussPointContribution(
				boost::numeric::ublas::bounded_matrix<double,18,18>& Kloc_system ,
    				const boost::numeric::ublas::bounded_matrix<double,3,3>& A,
				const boost::numeric::ublas::bounded_matrix<double,3,3>& B,				
 				const boost::numeric::ublas::bounded_matrix<double,3,3>& D,				
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
		//Membrane stiffness
		double beta0 = CalculateBeta( A );
		CalculateMembraneB(mBm, beta0,  loc1, loc2, loc3, x12, x23, x31, y12, y23, y31 );	
		CalculateMembraneContribution( mBm, A, mKloc99);
		AssembleMembraneContribution(   mKloc99, weight, Kloc_system  );
		
		//bending stiffness 
		CalculateBendingB(mBb, loc2, loc3, x12, x23, x31, y12, y23, y31 );
		CalculateBendingContribution( mBb, D, mKloc99);
		double bending_weight = weight; //attention!!
		AssembleBendingContribution(   mKloc99, bending_weight, Kloc_system  );
	
		//Membrane+bending stiffness
		CalculateMixedContribution( mBm, mBb, B, mKloc99);
		AssembleMixedContribution(   mKloc99, weight, Kloc_system  );
	
	}
	
	//************************************************************************************
	//************************************************************************************
	double ShellAnisotropicLinear::CalculateBeta( const 
boost::numeric::ublas::bounded_matrix<double,3,3>& A )
	{
		double W = -6.0*A(0,1)*A(0,1)*A(0,1) + 5.0*A(0,0)*A(0,0)*A(1,1) - 5.0*A(0,1)*A(0,1)*A(1,1)
					-A(1,1)*( 75.0*A(0,2)*A(0,2) + 14.0*A(0,2)*A(1,2) + 3.0*A(1,2)*A(1,2)  )
					+2.0*A(0,1)*( 7.0*A(0,2)*A(0,2)  + 46.0*A(0,2)*A(1,2) + 7.0*A(1,2)*A(1,2)  )	
					-A(0,0)*( 5.0*A(0,1)*A(0,1)+3.0*A(0,2)*A(0,2)-6.0*A(0,1)*A(1,1)-5.0*A(1,1)*A(1,1)+14.0*A(0,2)*A(1,2)+75.0*A(1,2)*A(1,2))	
						+( 3.0*A(0,0)+82.0*A(0,0)*A(1,1) + 3.0*A(1,1)*A(1,1)  
							-4.0*(6.0*A(0,1)*A(0,1) +5.0*A(0,2)*A(0,2) - 6.0*A(0,2)*A(1,2) + 5.0*A(1,2)*A(1,2) ) )*A(2,2)
					+4.0*(  5.0*A(0,0)-6.0*A(0,1)+5.0*A(1,1)  )*A(2,2)*A(2,2);
		
		double detA = MathUtils<double>::Det3(A);
		double beta0 = std::max(256.0*detA/W-1.5 , 0.01);
// 		double beta0 = std::max(fabs(256.0*detA/W-1.5) , 0.001);
// 		KRATOS_WATCH(beta0);
// KRATOS_WATCH(detA);
// KRATOS_WATCH(W);
// KRATOS_WATCH(256.0*detA/W-1.5);

		return beta0;
	}
	

	
	//************************************************************************************
	//************************************************************************************
	void ShellAnisotropicLinear::CalculateAllMatrices(	
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

		//getting material matrices - making a local copy so to allow parallelism
		noalias(msA) = GetProperties()[MATRIX_A]; //membrane
		noalias(msB) = GetProperties()[MATRIX_B]; //coupled membrane-bending
		noalias(msD) = GetProperties()[MATRIX_D]; //bending
// 		KRATOS_WATCH(msA)
// KRATOS_WATCH(msB)
// KRATOS_WATCH(msD)

		//here we calculate everything in the local system of coordinates
		//calculate integration point 1
		weight = A/3;
		loc1= 0.5;
		loc2 = 0.5;
		loc3 = 0.0;
		CalculateGaussPointContribution( mKloc_system, msA, msB, msD, weight, h, loc1, loc2, loc3, x12, x23, x31, y12, y23, y31);

		//calculate integration point 2
		weight = A/3;
		loc1= 0.0;
		loc2 = 0.5;
		loc3 = 0.5;
		CalculateGaussPointContribution( mKloc_system, msA, msB, msD, weight, h, loc1, loc2, loc3, x12, x23, x31, y12, y23, y31);
	
		//calculate integration point 3
		weight = A/3;
		loc1= 0.5;
		loc2 = 0.0;
		loc3 = 0.5;
		CalculateGaussPointContribution( mKloc_system, msA, msB, msD, weight, h, loc1, loc2, loc3, x12, x23, x31, y12, y23, y31);
	
		//rotate from local coordinates to global coordinates
		RotateToGlobal(v1,v2,v3,mKloc_system,rLeftHandSideMatrix);
		
		//calculate residua = b - Kglob*values;
		ApplyBodyForce(h,A,rRightHandSideVector);
		GetValuesVector(values,0);
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,values);
	} 
	
	//************************************************************************************
	//************************************************************************************
	void ShellAnisotropicLinear::EquationIdVector(EquationIdVectorType& 
rResult, ProcessInfo& CurrentProcessInfo)
	{
		int number_of_nodes = 3;
		if(rResult.size() != 18)
			rResult.resize(18,false);

		for (int i=0;i<number_of_nodes;i++)
		{
			int index = i*6;
			rResult[index]	   = 
GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
			rResult[index + 1] = 
GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
			rResult[index + 2] = 
GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();

			rResult[index + 3] = 
GetGeometry()[i].GetDof(ROTATION_X).EquationId();
			rResult[index + 4] = 
GetGeometry()[i].GetDof(ROTATION_Y).EquationId();
			rResult[index + 5] = 
GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
		}

	}

	//************************************************************************************
	//************************************************************************************
	void ShellAnisotropicLinear::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
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
	void ShellAnisotropicLinear::GetValuesVector(Vector& values, int Step)
	{
		const unsigned int number_of_nodes = 3;
		unsigned int MatSize = 18;
		if(values.size() != MatSize)	values.resize(MatSize,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			const array_1d<double,3>& disp = 
GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,Step);
			const array_1d<double,3>& rot = 
GetGeometry()[i].FastGetSolutionStepValue(ROTATION,Step);
				
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
	void ShellAnisotropicLinear::RotateToGlobal(
			const array_1d<double,3>& v1,
			const array_1d<double,3>& v2,
			const array_1d<double,3>& v3,
			const boost::numeric::ublas::bounded_matrix<double,18,18>& Kloc_system,
			Matrix& rLeftHandSideMatrix)
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
	void ShellAnisotropicLinear::ApplyBodyForce(	
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
		
		rRightHandSideVector[0] = bf[0];
		rRightHandSideVector[1] = bf[1];
		rRightHandSideVector[2] = bf[2];
		
		rRightHandSideVector[6] = bf[0];
		rRightHandSideVector[7] = bf[1];
		rRightHandSideVector[8] = bf[2];		
		
		rRightHandSideVector[12] = bf[0];
		rRightHandSideVector[13] = bf[1];
		rRightHandSideVector[14] = bf[2];		
		
		KRATOS_CATCH("");
		
	}
	//************************************************************************************
	//************************************************************************************
	void ShellAnisotropicLinear::NicePrint(const Matrix& A)
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
	void ShellAnisotropicLinear::CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& Output, const ProcessInfo& rCurrentProcessInfo)
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
			double h = GetProperties()[THICKNESS];
			//double h = 1.0; //note that we want the stress NOT the stress multiplied by the thickness
			

			//here we calculate everything in the local system of coordinates
			//calculate integration point 1
			weight = A;
			loc1 = 0.33333333333333;
			loc2 = 0.33333333333333;
			loc3 = 0.33333333333333;
			
			//direct membrane stress			
			double beta0 = 1.5; //note that this is just for stress evaluation!
			CalculateMembraneB(mBm, beta0,  loc1, loc2, loc3, x12, x23, x31, y12, y23, y31 );
			CalculateBendingB(mBb, loc2, loc3, x12, x23, x31, y12, y23, y31 );	

			array_1d<double,9> membrane_disp, bending_disp;
			array_1d<double,3> local_stress, local_strain;
			array_1d<double,6> rotated_stress = ZeroVector(6);
			
			//calculating stress in local coordinates
			const array_1d<double,3>& disp0 = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
			const array_1d<double,3>& disp1 = GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT);
			const array_1d<double,3>& disp2 = GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT);
			const array_1d<double,3>& rot0 = GetGeometry()[0].FastGetSolutionStepValue(ROTATION);
			const array_1d<double,3>& rot1 = GetGeometry()[1].FastGetSolutionStepValue(ROTATION);
			const array_1d<double,3>& rot2 = GetGeometry()[2].FastGetSolutionStepValue(ROTATION);
			
			membrane_disp[0] = inner_prod(disp0,v1);
			membrane_disp[1] = inner_prod(disp0,v2);
			membrane_disp[2] = inner_prod(rot0,v3);			
			membrane_disp[3] = inner_prod(disp1,v1);
			membrane_disp[4] = inner_prod(disp1,v2);
			membrane_disp[5] = inner_prod(rot1,v3);			
			membrane_disp[6] = inner_prod(disp2,v1);
			membrane_disp[7] = inner_prod(disp2,v2);
			membrane_disp[8] = inner_prod(rot2,v3);
			
			bending_disp[0] = inner_prod(disp0,v3); //node 1
			bending_disp[1] = inner_prod(rot0,v1);
			bending_disp[2] = inner_prod(rot0,v2);
			bending_disp[3] = inner_prod(disp1,v3); //node 2
			bending_disp[4] = inner_prod(rot1,v1);
			bending_disp[5] = inner_prod(rot1,v2);
			bending_disp[6] = inner_prod(disp2,v3); //node 3
			bending_disp[7] = inner_prod(rot2,v1);
			bending_disp[8] = inner_prod(rot2,v2);	
			
			//calculating material matrices
			noalias(msA) = GetProperties()[MATRIX_A]; //membrane
			noalias(msB) = GetProperties()[MATRIX_B]; //coupled membrane-bending
			
			//calculating membrane stress
			noalias(local_strain) = prod(trans(mBm),membrane_disp);
			noalias(local_stress) = prod(msA,local_strain);
			
			//adding "coupling" stress
			noalias(local_strain) = prod(trans(mBb),bending_disp);
			noalias(local_stress) += prod(msB,local_strain);
			
			local_stress /= h;
			
			//KRATOS_WATCH(local_stress);
			
			//rotating the stress to global coordinates
			AddVoigtTensorComponents(local_stress[0],rotated_stress,v1,v1);
			AddVoigtTensorComponents(local_stress[1],rotated_stress,v2,v2);			
			AddVoigtTensorComponents(local_stress[2],rotated_stress,v1,v2);
			AddVoigtTensorComponents(local_stress[2],rotated_stress,v2,v1);
			

			if(Output[0].size2() != 6)
				Output[0].resize(1,6,false);
			
			for(unsigned int ii = 0; ii<rotated_stress.size(); ii++)
				Output[0](0,ii) = rotated_stress[ii];
			
			
//PROVA!!! CANCELLARE
// 			std::cout << "testing" << std::endl;
// 			Matrix Kg;
// 			this->Calculate(GEOMETRIC_STIFFNESS,Kg,rCurrentProcessInfo);
// 			KRATOS_WATCH(Kg);

		}

		KRATOS_CATCH("")

	}

	//***********************************************************************************
	//***********************************************************************************
	//auxiliary function needed in the calculation of output stresses
	inline void ShellAnisotropicLinear::AddVoigtTensorComponents(
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
}



