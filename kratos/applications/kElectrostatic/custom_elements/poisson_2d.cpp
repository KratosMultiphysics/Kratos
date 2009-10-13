/*
==============================================================================
KratosR1PoissonApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2008
Pooyan Dadvand, Riccardo Rossi, Javier Mora
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
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/poisson_2d.h"
#include "kPoisson.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 

namespace Kratos
{
	//static variables
	boost::numeric::ublas::bounded_matrix<double,3,2> Poisson2D::msDN_DX;
	boost::numeric::ublas::bounded_matrix<double,2,2> Poisson2D::msD;
  	array_1d<double,3> Poisson2D::msN; //dimension = number of nodes
	array_1d<double,3> Poisson2D::ms_temp; //dimension = number of nodes
	array_1d<double,3> Poisson2D::point_sources; //dimension = number of nodes

	//************************************************************************************
	//************************************************************************************
	Poisson2D::Poisson2D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	Poisson2D::Poisson2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{

	}

	Element::Pointer Poisson2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new Poisson2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	Poisson2D::~Poisson2D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void Poisson2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		const unsigned int number_of_points = GetGeometry().size();

		if(rLeftHandSideMatrix.size1() != number_of_points)
			rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);

		if(rRightHandSideVector.size() != number_of_points)
			rRightHandSideVector.resize(number_of_points,false);

		//getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);

		//reading properties and conditions
		double permittivity = GetProperties()[DUMMY_MATERIAL];
		msD(0,0)=permittivity;
		msD(1,1)=permittivity;
		point_sources[0] = GetGeometry()[0].FastGetSolutionStepValue(DUMMY_POINT_SOURCE);
		point_sources[1] = GetGeometry()[1].FastGetSolutionStepValue(DUMMY_POINT_SOURCE);
		point_sources[2] = GetGeometry()[2].FastGetSolutionStepValue(DUMMY_POINT_SOURCE);
		
		// main loop	
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

		noalias(rLeftHandSideMatrix) = prod(msDN_DX,Matrix(prod(msD,trans(msDN_DX))));

		rLeftHandSideMatrix *= Area;

		//Point charge contribution 
		noalias(rRightHandSideVector) = point_sources;  


		//subtracting the dirichlet term
		// RHS -= LHS*DUMMY_UNKNOWNs
		for(unsigned int iii = 0; iii<number_of_points; iii++)
			ms_temp[iii] = GetGeometry()[iii].FastGetSolutionStepValue(DUMMY_UNKNOWN);
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp);
		
		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void Poisson2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_ERROR(std::logic_error,  "method not implemented" , "");
	}	 

	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void Poisson2D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void Poisson2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);	

		for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(DUMMY_UNKNOWN).EquationId();
	}

	//************************************************************************************
	//************************************************************************************
	  void Poisson2D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);	

		for (unsigned int i=0;i<number_of_nodes;i++)
			ElementalDofList[i] = GetGeometry()[i].pGetDof(DUMMY_UNKNOWN);

	}



} // Namespace Kratos