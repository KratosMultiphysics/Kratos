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
//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: kazem $
//   Date:                $Date: 2008-07-24 16:46:55 $
//   Revision:            $Revision: 1.1 $
//
//


// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_conditions/spring.h"
#include "structural_application.h"
#include "utilities/math_utils.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	Spring2D::Spring2D(IndexType NewId, GeometryType::Pointer 
        pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	Spring2D::Spring2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}

	Condition::Pointer Spring2D::Create(IndexType NewId, NodesArrayType 
        const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
          return Condition::Pointer(new Spring2D(NewId, 
          GetGeometry().Create(ThisNodes), pProperties));
	}

	Spring2D::~Spring2D()
	{
	}

        void Spring2D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
        {
	  
	  mfail = false;
	}

	//************************************************************************************
	//************************************************************************************
	void Spring2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		if(rRightHandSideVector.size() != 4)
			rRightHandSideVector.resize(4,false);
		
		double oc  = 1.5E-5; /// desplazamiento critico
		double miu = 0.60;  /// coeficiente de friccion
		
		
		noalias(rRightHandSideVector) = ZeroVector(4);
		
		if(mfail==false){
		array_1d<double,2> dcoord1;
		array_1d<double,2> dcoord2;
		array_1d<double,2> tangential;
		array_1d<double,2> normal;
	
		array_1d<double,3>& normal_1 = GetGeometry()[0].GetSolutionStepValue(NORMAL);
		array_1d<double,3>& normal_2 = GetGeometry()[1].GetSolutionStepValue(NORMAL);
		array_1d<double,3>& spring   = GetGeometry()[0].GetSolutionStepValue(SPRING_STIFFNESS);
		
		 normal_1   = ZeroVector(3);
		 normal_2   = ZeroVector(3);
		 
		 array_1d<double,2> d1_inicial;
		 array_1d<double,2> d2_inicial;
		 array_1d<double,2> dinicial;
		 d1_inicial[0]     = GetGeometry()[0].X0(); 
		 d1_inicial[1]     = GetGeometry()[0].Y0();
		 d2_inicial[0]     = GetGeometry()[1].X0(); 
		 d2_inicial[1]     = GetGeometry()[1].Y0();
		 noalias(dinicial) = d2_inicial - d1_inicial;  
 
		 
		 array_1d<double,2> d1_current;
		 array_1d<double,2> d2_current;
		 array_1d<double,2> dcurrent;
		 d1_current[0]     = GetGeometry()[0].X(); 
		 d1_current[1]     = GetGeometry()[0].Y();
		 d2_current[0]     = GetGeometry()[1].X(); 
		 d2_current[1]     = GetGeometry()[1].Y();
		 noalias(dcurrent) = d2_current - d1_current;  

		 double h    = std::sqrt(inner_prod(dinicial, dinicial));
		 noalias(tangential) = (1.00/h) * dinicial;
		 
		 normal[0] =  -tangential[1]; 
		 normal[1] =   tangential[0];
		 
		 double d_o       = inner_prod(dinicial, dinicial);  
		 double d_c       = inner_prod(dcurrent, dcurrent);  
		 double  o        = std::sqrt(d_c) - std::sqrt(d_o);    
		  
		double spring_force  = spring[0] * o; 
		double tan_force     = miu * spring_force;  
		 
		/// spring traccionado
		
		if(o>oc){ mfail = true; return;}
		else if(o<oc){
		 normal_1[0] = rRightHandSideVector[0] =  spring_force * tangential[0];
		 normal_1[1] = rRightHandSideVector[1] =  spring_force * tangential[1];
		 normal_2[0] = rRightHandSideVector[2] = -spring_force * tangential[0];
		 normal_2[1] = rRightHandSideVector[3] = -spring_force * tangential[1];
		}
		///spring comprimido 
		else if(o<0.00)
		{
		  
		 /// componente normal 
		 normal_1[0] = rRightHandSideVector[0] =  spring_force * tangential[0];
		 normal_1[1] = rRightHandSideVector[1] =  spring_force * tangential[1];
		 normal_2[0] = rRightHandSideVector[2] = -spring_force * tangential[0];
		 normal_2[1] = rRightHandSideVector[3] = -spring_force * tangential[1];
		 
		 /// componente tangencial debido a la friccion
		 rRightHandSideVector[0] =  tan_force * normal[0];
		 rRightHandSideVector[1] =  tan_force * normal[1];
		 rRightHandSideVector[2] = -tan_force * normal[0];
		 rRightHandSideVector[3] = -tan_force * normal[1];
		 }
		}
		KRATOS_CATCH("")
	      }
	

	//************************************************************************************
	//************************************************************************************
	void Spring2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
	         KRATOS_TRY
	         /*
		 if(rLeftHandSideMatrix.size1() != 4)
			rLeftHandSideMatrix.resize(4,4,false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(4,4);

		if(rRightHandSideVector.size() != 4)
			rRightHandSideVector.resize(4,false);

		
	        array_1d<double,2> coordx;
		array_1d<double,2> coordy;
		array_1d<double,3>& spring = GetGeometry()[0].GetSolutionStepValue(SPRING_STIFFNESS);
		
		coordx[0] = GetGeometry()[0].X();  
		coordx[1] = GetGeometry()[1].X();  
		coordy[0] = GetGeometry()[0].Y();  
		coordy[1] = GetGeometry()[1].Y();  
		
		double  o = coordx[1] - coordx[0];
		double  s = coordy[1] - coordy[0];
		
		rRightHandSideVector[0]   =  spring[0] * o;
		rRightHandSideVector[1]   =  spring[1] * s;
		rRightHandSideVector[2]   = -spring[0] * o;
		rRightHandSideVector[3]   = -spring[1]* s;
		*/
		KRATOS_CATCH("")
	}


	//************************************************************************************
	//************************************************************************************
	void Spring2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int index;
		unsigned int dim = 2;
		rResult.resize(number_of_nodes*dim);
		for (int i=0;i<number_of_nodes;i++)
		{
			index = i*dim;
			rResult[index] = (GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId());
			rResult[index+1] = (GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId());
			
		}
	}

	//************************************************************************************
	//************************************************************************************
	  void Spring2D::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int dim = 2;
		ConditionalDofList.resize(GetGeometry().size()*dim);
		unsigned int index;
		for (unsigned int i=0;i<GetGeometry().size();i++)
		{
			
			index = i*dim;
			ConditionalDofList[index] = (GetGeometry()[i].pGetDof(DISPLACEMENT_X));
			ConditionalDofList[index+1] = (GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
			
		}
	}
} // Namespace Kratos

 

