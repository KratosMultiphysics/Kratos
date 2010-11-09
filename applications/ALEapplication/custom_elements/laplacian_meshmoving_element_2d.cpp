/*
==============================================================================
KratosALEApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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
//   Date:                $Date: 2007-03-06 10:30:31 $
//   Revision:            $Revision: 1.2 $
//
//


// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/laplacian_meshmoving_element_2d.h"
#include "ale_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{
	boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX;
  	array_1d<double,3> msN; 
  	array_1d<double,3> ms_temp_vec_np; 
	//************************************************************************************
	//************************************************************************************
	LaplacianMeshMovingElem2D::LaplacianMeshMovingElem2D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	LaplacianMeshMovingElem2D::LaplacianMeshMovingElem2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
	}

	Element::Pointer LaplacianMeshMovingElem2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new LaplacianMeshMovingElem2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	LaplacianMeshMovingElem2D::~LaplacianMeshMovingElem2D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void LaplacianMeshMovingElem2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_ERROR(std::logic_error,  "method not implemented" , "");
	}	
	
  
	//************************************************************************************
	//************************************************************************************
	void LaplacianMeshMovingElem2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_points = 3;

		if(rLeftHandSideMatrix.size1() != number_of_points)
			rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);

		if(rRightHandSideVector.size() != number_of_points)
			rRightHandSideVector.resize(number_of_points,false);

		unsigned int ComponentIndex = rCurrentProcessInfo[FRACTIONAL_STEP] - 1;
		
		boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX;
		array_1d<double,3> msN; 
		array_1d<double,3> ms_temp_vec_np; 
	
		//getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);
					
		const array_1d<double,3>& disp0 = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
		const array_1d<double,3>& disp1 = GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT);
		const array_1d<double,3>& disp2 = GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT);

		noalias(rLeftHandSideMatrix) = prod(msDN_DX,trans(msDN_DX));

		//dirichlet contribution
		ms_temp_vec_np[0] = disp0[ComponentIndex]; 
		ms_temp_vec_np[1] = disp1[ComponentIndex]; 
		ms_temp_vec_np[2] = disp2[ComponentIndex]; 
		noalias(rRightHandSideVector) = ZeroVector(number_of_points);
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp_vec_np);

		//note that no multiplication by area is performed,
		//this makes smaller elements more rigid and minimizes the mesh deformation

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void LaplacianMeshMovingElem2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);	

		for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(AUX_MESH_VAR).EquationId();
	}

	//************************************************************************************
	//************************************************************************************
	  void LaplacianMeshMovingElem2D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);	

		for (unsigned int i=0;i<number_of_nodes;i++)
			ElementalDofList[i] = GetGeometry()[i].pGetDof(AUX_MESH_VAR);

	}

	

} // Namespace Kratos


