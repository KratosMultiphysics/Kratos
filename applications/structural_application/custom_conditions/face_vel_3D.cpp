/*
==============================================================================
KratosStructuralApplication 
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
/* **************************************************************************************
*          
*   Last Modified by:    $Author: rrossi $
*   Date:                $Date: 2008-02-14 09:41:09 $
*   Revision:            $Revision: 1.7 $
*
* ***************************************************************************************/


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/face_vel_3D.h"
#include "structural_application.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"

namespace Kratos
{

	//***********************************************************************************
	//***********************************************************************************
		  // -------- //
		 //  PUBLIC  //
		// -------- //

	// Constructor
	FaceVel3D::FaceVel3D()
	{
	}

	// Constructor
	FaceVel3D::FaceVel3D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Face3D(NewId, pGeometry)
	{
	}

	// Constructor
	FaceVel3D::FaceVel3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
		: Face3D(NewId, pGeometry, pProperties)

	{
		//setting up the nodal degrees of freedom
/*		for(unsigned int i = 0 ; i != GetGeometry().size() ; ++i)
		{
			(GetGeometry()[i].pAddDof(DISPLACEMENT_X,REACTION_X));
			(GetGeometry()[i].pAddDof(DISPLACEMENT_Y,REACTION_Y));
			(GetGeometry()[i].pAddDof(DISPLACEMENT_Z,REACTION_Z));
		}
*/
	}

	//***********************************************************************************
	//***********************************************************************************

	Condition::Pointer FaceVel3D::Create(
		IndexType NewId,
		NodesArrayType const& ThisNodes,
		PropertiesType::Pointer pProperties) const

	{
KRATOS_TRY

		return Condition::Pointer(new FaceVel3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
KRATOS_CATCH("")
	}

	//***********************************************************************************
	//***********************************************************************************
	// Destructor
	FaceVel3D::~FaceVel3D()
	{
	}

	//***********************************************************************************
	//***********************************************************************************

	void FaceVel3D::EquationIdVector(
		EquationIdVectorType& rResult,
		ProcessInfo& rCurrentProcessInfo)

	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().size();
		unsigned int dim = number_of_nodes*3;
		if(rResult.size() != dim)
            rResult.resize(dim);

		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			int index = i*3;
			rResult[index]   = GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
			rResult[index+1] = GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
			rResult[index+2] = GetGeometry()[i].GetDof(VELOCITY_Z).EquationId();
		}
		KRATOS_CATCH("")
	}

	//***********************************************************************************
	//***********************************************************************************

	void FaceVel3D::GetDofList(
		DofsVectorType& ElementalDofList,
		ProcessInfo& rCurrentProcessInfo)
	
	{
		ElementalDofList.resize(0);

		for (unsigned int i=0;i<GetGeometry().size();i++)
		{
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_X));
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Y));
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(VELOCITY_Z));
		}
	}

	//***********************************************************************************
	//***********************************************************************************


}	// Namespace Kratos.
