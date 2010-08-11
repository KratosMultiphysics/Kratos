/*
==============================================================================
KratosR1MagnetostaticApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2010
Pooyan Dadvand, Riccardo Rossi, Javier Mora
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
mora@cimne.upc.edu
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
//   Last modified by:    $Author: rrossi  jmora $
//   Date:                $Date: 2010-02-02 $
//   Revision:            $Revision: 1.4 $
//
//
 

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/magnetostatic_3d.h"
#include "kMagnetostatic.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{
	namespace Magnetostatic3Dauxiliaries
    {
		//static variables
		boost::numeric::ublas::bounded_matrix<double,4,3> msDN_DX = ZeroMatrix(4,3);
		#pragma omp threadprivate(msDN_DX)
		boost::numeric::ublas::bounded_matrix<double,12,12> msB = ZeroMatrix(12,12);
		#pragma omp threadprivate(msB)
		boost::numeric::ublas::bounded_matrix<double,3,3> msD = ZeroMatrix(3,3);
		#pragma omp threadprivate(msD)
		array_1d<double,3> msHc = ZeroVector(3); //dimension = of the domain
		#pragma omp threadprivate(msHc)

		array_1d<double,4> msN = ZeroVector(4); //dimension = number of nodes
		#pragma omp threadprivate(msN)
		array_1d<double,12> ms_temp = ZeroVector(12); //dimension = number of nodes
		#pragma omp threadprivate(ms_temp)
		array_1d<double,4> point_sources = ZeroVector(4); //dimension = number of nodes
		#pragma omp threadprivate(point_sources)

//		array_1d<double,4> volume_sources = ZeroVector(4); //dimension = number of nodes
		array_1d<double,3> volume_sources = ZeroVector(3); //dimension = number of dof
		#pragma omp threadprivate(volume_sources)
    }
    using  namespace Magnetostatic3Dauxiliaries;


	//************************************************************************************
	//************************************************************************************
	Magnetostatic3D::Magnetostatic3D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	Magnetostatic3D::Magnetostatic3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{

	}

	Element::Pointer Magnetostatic3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new Magnetostatic3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	Magnetostatic3D::~Magnetostatic3D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void Magnetostatic3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		const unsigned int number_of_points = GetGeometry().size();

		if(rLeftHandSideMatrix.size1() != number_of_points)
			rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);

		if(rRightHandSideVector.size() != number_of_points)
			rRightHandSideVector.resize(number_of_points,false);

		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dimension = GetGeometry().WorkingSpaceDimension();

		unsigned int MatSize=number_of_nodes*dimension;

		if(rLeftHandSideMatrix.size1() != MatSize)
			rLeftHandSideMatrix.resize(MatSize,MatSize,false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize,MatSize); //resetting LHS

		//resizing as needed the RHS
		if(rRightHandSideVector.size() != MatSize)
			rRightHandSideVector.resize(MatSize,false);
		rRightHandSideVector = ZeroVector(MatSize); //resetting RHS

		//getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);

		//reading properties and conditions
		array_1d<double,3> permeability = GetProperties()[MAGNETIC_PERMEABILITY];
		msD(0,0)=1.0/permeability[0];
		msD(1,1)=1.0/permeability[1];
		msD(2,2)=1.0/permeability[2];
		msD(0,1)=0.0;
		msD(0,2)=0.0;
		msD(1,0)=0.0;
		msD(1,2)=0.0;
		msD(2,0)=0.0;
		msD(2,1)=0.0;

		//double surface_sources = (this)->GetValue(MAGNETOSTATIC_SURFACE_CURRENT);
		volume_sources = (this)->GetValue(MAGNETOSTATIC_VOLUME_CURRENT);
		//KRATOS_WATCH(volume_sources);

		//volume_sources[1] = 0.0;//(this)->GetValue(MAGNETOSTATIC_VOLUME_CURRENT);
		//volume_sources[2] = 0.0;//(this)->GetValue(MAGNETOSTATIC_VOLUME_CURRENT);
		//volume_sources[3] = 0.0;//(this)->GetValue(MAGNETOSTATIC_VOLUME_CURRENT);
		
		// main loop	
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

		// building stiffness matrix
	
		unsigned int dim2 = number_of_nodes*dimension;

		double penalty=0.0;
		double index_m,index_n;
		array_1d<double,3> DNI;
		array_1d<double,3> DNJ;
		array_1d<double,3> paux;

		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			DNI[0]=msDN_DX(i,0);
			DNI[1]=msDN_DX(i,1);
			DNI[2]=msDN_DX(i,2);

			//KRATOS_WATCH(i);

			for (unsigned int j=0;j<number_of_nodes;j++)
			{
				//KRATOS_WATCH(j);
				//KRATOS_WATCH(msDN_DX);
				DNJ[0]=msDN_DX(j,0);
				DNJ[1]=msDN_DX(j,1);
				DNJ[2]=msDN_DX(j,2);

				for (unsigned int k=0;k<dimension;k++)
				{
					//KRATOS_WATCH(k);

					for (unsigned int l=0;l<dimension;l++)
					{
						paux[0]=1.0;
						paux[1]=1.0;
						paux[2]=1.0;

						//KRATOS_WATCH(l);

						index_m=i*dimension+k;
						index_n=j*dimension+l;

						//KRATOS_WATCH(index_m);
						//KRATOS_WATCH(index_n);

						if(k==l)
						{
							paux(k)=penalty;
							//KRATOS_WATCH(k);
							//KRATOS_WATCH(paux);
							
							msB(index_m,index_n)=
								paux[0]*DNI[0]*DNJ[0]*msD(0,0)+
								paux[1]*DNI[1]*DNJ[1]*msD(1,1)+
								paux[2]*DNI[2]*DNJ[2]*msD(2,2);
/*							msB(index_m,index_n)=
								DNI[0]*DNJ[0]+
								DNI[1]*DNJ[1]+
								DNI[2]*DNJ[2];
							msB(index_m,index_n)=
								paux[0]*msDN_DX(i,0)*msDN_DX(j,0)+
								paux[1]*msDN_DX(i,1)*msDN_DX(j,1)+
								paux[2]*msDN_DX(i,2)*msDN_DX(j,2);
*/						}
						else
						{
/*							msB(index_m,index_n)=
								penalty*DNI[k]*DNJ[l]*msD(k,k)-
								DNI[l]*DNJ[k]*msD(l,l);
							msB(index_m,index_n)=
								DNI[k]*DNJ[l]-
								DNI[l]*DNJ[k];
*/							msB(index_m,index_n)=0.0;
//								penalty*DNI[k]*DNJ[l]-
//								DNI[l]*DNJ[k];
						}
					}
				}
			}
		}

		//noalias(rLeftHandSideMatrix) = prod(msDN_DX,Matrix(prod(msD,trans(msDN_DX))));
		noalias(rLeftHandSideMatrix) = msB*Area;

		array_1d<double,3> coercivity = GetProperties()[COERCIVITY];
		msHc[0]=coercivity[0];
		msHc[1]=coercivity[1];
		msHc[2]=coercivity[2];

		int dim = GetGeometry().WorkingSpaceDimension();

		for(unsigned int iii = 0; iii<number_of_points; iii++)
		{
			int index = iii*dim;

			ms_temp[index] = (msDN_DX(iii,1)*msHc[2]-msDN_DX(iii,2)*msHc[1])+volume_sources[0]/4.0;
			ms_temp[index+1] = (msDN_DX(iii,2)*msHc[0]-msDN_DX(iii,0)*msHc[2])+volume_sources[1]/4.0;;
			ms_temp[index+2] = (msDN_DX(iii,0)*msHc[1]-msDN_DX(iii,1)*msHc[0])+volume_sources[2]/4.0;;
		}

		//noalias(rRightHandSideVector) = (prod(msDN_DX,msHc)+volume_sources/4.0)*Area; 
		//noalias(rRightHandSideVector) = ms_temp+volume_sources*0.0;//(volume_sources/4.0)*Area; 
		noalias(rRightHandSideVector) = ms_temp*Area; 

		//subtracting the dirichlet term
		// RHS -= LHS*Magnetostatic_POTENTIALs

		for(unsigned int iii = 0; iii<number_of_points; iii++)
		{
			int index = iii*dim;

//						rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
			ms_temp[index] = GetGeometry()[iii].FastGetSolutionStepValue(MAGNETOSTATIC_VECTOR_POTENTIAL_X);
			ms_temp[index+1] = GetGeometry()[iii].FastGetSolutionStepValue(MAGNETOSTATIC_VECTOR_POTENTIAL_Y);
			ms_temp[index+2] = GetGeometry()[iii].FastGetSolutionStepValue(MAGNETOSTATIC_VECTOR_POTENTIAL_Z);

//			ms_temp[iii] = GetGeometry()[iii].FastGetSolutionStepValue(MAGNETOSTATIC_VECTOR_POTENTIAL);
		}
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp);
		
		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void Magnetostatic3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_ERROR(std::logic_error,  "method not implemented" , "");
	}	 

	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void Magnetostatic3D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_CATCH("");
	}


	//************************************************************************************
	//************************************************************************************
	void Magnetostatic3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		int number_of_nodes = GetGeometry().size();
		int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int dim2 = number_of_nodes*dim;
		if(rResult.size() != dim2)
			rResult.resize(dim2,false);

		for (int i=0;i<number_of_nodes;i++)
		{
			int index = i*dim;
			rResult[index] = GetGeometry()[i].GetDof(MAGNETOSTATIC_VECTOR_POTENTIAL_X).EquationId();
			rResult[index + 1] = GetGeometry()[i].GetDof(MAGNETOSTATIC_VECTOR_POTENTIAL_Y).EquationId();
			if(dim == 3)
				rResult[index + 2] = GetGeometry()[i].GetDof(MAGNETOSTATIC_VECTOR_POTENTIAL_Z).EquationId();
		}
	}

	//************************************************************************************
	//************************************************************************************
	  void Magnetostatic3D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
//		unsigned int number_of_nodes = GetGeometry().PointsNumber();
//		if(ElementalDofList.size() != number_of_nodes)
//			ElementalDofList.resize(number_of_nodes);	

//		for (unsigned int i=0;i<number_of_nodes;i++)
//			ElementalDofList[i] = GetGeometry()[i].pGetDof(MAGNETOSTATIC_VECTOR_POTENTIAL);

		ElementalDofList.resize(0);

		for (unsigned int i=0;i<GetGeometry().size();i++)
		{
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(MAGNETOSTATIC_VECTOR_POTENTIAL_X));
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(MAGNETOSTATIC_VECTOR_POTENTIAL_Y));
			if(GetGeometry().WorkingSpaceDimension() == 3)
            {
				ElementalDofList.push_back(GetGeometry()[i].pGetDof(MAGNETOSTATIC_VECTOR_POTENTIAL_Z));
            }
		}
	}


	//************************************************************************************
	//************************************************************************************
	void Magnetostatic3D::CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3> >& Output, const ProcessInfo& rCurrentProcessInfo)
	{
		IntegrationMethod mThisIntegrationMethod;

        mThisIntegrationMethod= GetGeometry().GetDefaultIntegrationMethod();
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);

		Vector vect_tmp(12);
		array_1d<double,4> Ax;
		array_1d<double,4> Ay;
		array_1d<double,4> Az;

		const unsigned int number_of_points = GetGeometry().size();
		int dim = GetGeometry().WorkingSpaceDimension();

		for(unsigned int iii = 0; iii<number_of_points; iii++)
		{
			int index = iii*dim;

			vect_tmp[index] = GetGeometry()[iii].FastGetSolutionStepValue(MAGNETOSTATIC_VECTOR_POTENTIAL_X);
			vect_tmp[index+1] = GetGeometry()[iii].FastGetSolutionStepValue(MAGNETOSTATIC_VECTOR_POTENTIAL_Y);
			vect_tmp[index+2] = GetGeometry()[iii].FastGetSolutionStepValue(MAGNETOSTATIC_VECTOR_POTENTIAL_Z);
			Ax[iii]=vect_tmp[index];
			Ay[iii]=vect_tmp[index+1];
			Az[iii]=vect_tmp[index+2];

			//			ms_temp[iii] = GetGeometry()[iii].FastGetSolutionStepValue(MAGNETOSTATIC_VECTOR_POTENTIAL);
		}

//		vect_tmp[0] = GetGeometry()[0].FastGetSolutionStepValue(MAGNETOSTATIC_VECTOR_POTENTIAL);
//		vect_tmp[1] = GetGeometry()[1].FastGetSolutionStepValue(MAGNETOSTATIC_VECTOR_POTENTIAL);
//		vect_tmp[2] = GetGeometry()[2].FastGetSolutionStepValue(MAGNETOSTATIC_VECTOR_POTENTIAL);
//		vect_tmp[3] = GetGeometry()[3].FastGetSolutionStepValue(MAGNETOSTATIC_VECTOR_POTENTIAL);

		//reading properties and conditions
		array_1d<double,3> permeability = GetProperties()[MAGNETIC_PERMEABILITY];
		msD(0,0)=1.0/permeability[0];
		msD(1,1)=1.0/permeability[1];
		msD(2,2)=1.0/permeability[2];
		msD(0,1)=0.0;
		msD(0,2)=0.0;
		msD(1,0)=0.0;
		msD(1,2)=0.0;
		msD(2,0)=0.0;
		msD(2,1)=0.0;

		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);

		if(Output.size() != GetGeometry().IntegrationPoints(mThisIntegrationMethod).size())
            Output.resize(GetGeometry().IntegrationPoints(mThisIntegrationMethod).size());

		double swap_aux;
		Vector dAx(4);
		Vector dAy(4);
		Vector dAz(4);

		array_1d<double,3> coercivity = GetProperties()[COERCIVITY];
		msHc[0]=coercivity[0];
		msHc[1]=coercivity[1];
		msHc[2]=coercivity[2];

		for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
		{
			dAx=prod(trans(msDN_DX),Ax);
			dAy=prod(trans(msDN_DX),Ay);
			dAz=prod(trans(msDN_DX),Az);

			if(rVariable==MAGNETIC_FLUX_DENSITY)
			{

				Output[PointNumber][0]=dAz[1]-dAy[2];
				Output[PointNumber][1]=dAx[2]-dAz[0];
				Output[PointNumber][2]=dAy[0]-dAx[1];

//				noalias(Output[PointNumber])=-prod(trans(msDN_DX),vect_tmp);
//				swap_aux = Output[PointNumber][0];
//				Output[PointNumber][0]=Output[PointNumber][1];
//				Output[PointNumber][1]=-swap_aux;
//				Output[PointNumber][2]=0.0;
//				KRATOS_WATCH(Output[PointNumber]);
			}
			else if(rVariable==MAGNETIC_FIELD_INTENSITY)
			{
				Output[PointNumber][0]=(dAz[1]-dAy[2])*msD(0,0)-msHc[0];
				Output[PointNumber][1]=(dAx[2]-dAz[0])*msD(1,1)-msHc[1];
				Output[PointNumber][2]=(dAy[0]-dAx[1])*msD(2,2)-msHc[2];

				//KRATOS_WATCH(rValues[PointNumber]);
//				noalias(Output[PointNumber])=-prod(prod(msD,trans(msDN_DX)),vect_tmp);
//				swap_aux = Output[PointNumber][0];
//				Output[PointNumber][0]=Output[PointNumber][1];
//				Output[PointNumber][1]=-swap_aux;
//				Output[PointNumber][2]=0.0;
//				KRATOS_WATCH(Output[PointNumber]);
			}
			else 
			{
				Output[PointNumber][0]=1/msD(0,0);
				Output[PointNumber][1]=1/msD(1,1);
				Output[PointNumber][2]=1/msD(2,2);
//				KRATOS_WATCH(Output[PointNumber]);
			}
		}



	}


 	//************************************************************************************
    //************************************************************************************

    void Magnetostatic3D::GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable,
		 std::vector<array_1d<double,3> >& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
		//KRATOS_WATCH("GiD Post Magnetostatic - GetValueOnIntegrationPoints");

        if(rVariable==MAGNETIC_FLUX_DENSITY)
        {
//			KRATOS_WATCH("GiD Post Magnetostatic - get - elec-field");
            CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
        }
        else if(rVariable==MAGNETIC_FIELD_INTENSITY)
        {
//			KRATOS_WATCH("GiD Post Magnetostatic - get - elec-disp");
            CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
        }
		else
        {
//			KRATOS_WATCH("GiD Post Magnetostatic - get - else");
            CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
        }
    }


} // Namespace Kratos


