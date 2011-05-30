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
//   Last modified by:    $Author: virginia $
//   Date:                $Date: 2009-01-23 14:39:59 $
//   Revision:            $Revision: 1.27 $
//
//


// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/pfem_contact_element3D.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 


namespace Kratos
{

	//************************************************************************************
	//************************************************************************************
	PfemContactElement3D::PfemContactElement3D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}
	//************************************************************************************
	//************************************************************************************
    PfemContactElement3D::PfemContactElement3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {
	mcontact_is_active = false;	
	mpenetration = 100.0;
    }
	//************************************************************************************
	//************************************************************************************
	Element::Pointer PfemContactElement3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
        KRATOS_TRY
		return Element::Pointer(new PfemContactElement3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
KRATOS_CATCH("");
	}


	//constructor
	PfemContactElement3D::~PfemContactElement3D()
	{
	}


	
	//************************************************************************************
	//************************************************************************************
	void PfemContactElement3D::Initialize()
	{
		KRATOS_TRY

		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void PfemContactElement3D::CalculateAll(MatrixType& rLeftHandSideMatrix, 
                                           VectorType& rRightHandSideVector, 
                                           ProcessInfo& rCurrentProcessInfo,
                                           bool CalculateStiffnessMatrixFlag,
                                           bool CalculateResidualVectorFlag)
        {
		KRATOS_TRY
		boost::numeric::ublas::bounded_matrix<double, 4, 3 > DN_DX;
		array_1d<double,4> N;

		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize=number_of_nodes*dim;

		if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
		{
			if(rLeftHandSideMatrix.size1() != MatSize)
				rLeftHandSideMatrix.resize(MatSize,MatSize,false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize,MatSize); //resetting LHS
		}

		//resizing as needed the RHS
		if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
		{
			if(rRightHandSideVector.size() != MatSize)
				rRightHandSideVector.resize(MatSize,false);
			rRightHandSideVector = ZeroVector(MatSize); //resetting RHS
		}
		
		const double l0 = GetProperties()[THICKNESS];

		//count external faces
		unsigned int numb_ext_faces = 0;
		bool is_contact = false;
	       Geometry< Node<3> >& geom = GetGeometry();
		array_1d<double,3> n;
		double Volume;
		double h = 1.0e6;
		GeometryUtils::CalculateGeometryData(geom,DN_DX,N,Volume);

		//find reference face. If a reference face is not found, then do nothing
		int reference_face = -1;
		WeakPointerVector< Element >& neigb_el = this->GetValue(NEIGHBOUR_ELEMENTS);

	 	int ss3 = 0.0;
		
		double accepted = 1.0;
		FlagVariableCheckForNonSuitableElements(accepted);
		
	 if(Volume > 0.0 && accepted == 1.0){	
		for(unsigned int i=0;i<4; i++)
		{
		      if(neigb_el[i].Id() == this->Id()) //missing neighbour
			{
				numb_ext_faces++;
				reference_face=i;
			}
		}
		
		if(numb_ext_faces > 1) //this is most probably a sliver
		{
			h = 1e6;
			n[0] = 1.0; n[1] = 1.0; n[2] = 1.0;
			array_1d<double,3> auxn;
			is_contact = false;
			for(unsigned int i=0;i<4; i++)
			{
			  if(neigb_el[i].Id() == this->Id()) //missing neighbour
			  {
				reference_face = i;
				double aaa = norm_2(row(DN_DX,reference_face));
				h = 1.0/aaa;
				noalias(n) = row(DN_DX,reference_face);
				n /= aaa;

				//accept only if projection is inside
				bool image_inside = Check_image_inside_face(n, reference_face, GetGeometry(),true); 
				  if(image_inside)
				  {
/*				    KRATOS_WATCH(h);*/
// 				    if(h<l0){
// 				     KRATOS_WATCH("numb_ext_faces > 1");
// 				     KRATOS_WATCH(h);
// 				    }
// 				    
// 				    }
/*KRATOS_WATCH("numb_ext_faces > 1");
KRATOS_WATCH(h);
KRATOS_WATCH(n);	*/			    
					    is_contact = true;
					    break;
				  }    
			    }
			}
		}
		else if (numb_ext_faces == 1) //in this case i have a reference face
		{
			double aaa = norm_2(row(DN_DX,reference_face));
									 
			h = 1.0/aaa;
			noalias(n) = row(DN_DX,reference_face);
			n /= aaa;

			//accept only if projection is inside
// 			bool image_inside = Check_image_inside_face(n, reference_face, GetGeometry(),false); 
// KRATOS_WATCH("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<before detect");
			DetectContact(GetGeometry(),DN_DX,reference_face, n, h);
			bool image_inside = true;

			   if(image_inside){
				    is_contact = true;
/*				    				    KRATOS_WATCH(h);*/
// 				    if(h<l0){
// 				     KRATOS_WATCH("numb_ext_faces == 1");
// 				     KRATOS_WATCH(h);
// 				    }
/*KRATOS_WATCH(h);
KRATOS_WATCH(n);*/				    
			   }
			   else
				    is_contact = false;

		}
		else //case B - one edge on each side
		{	h = 1e6;
			n[0] = 1.0; n[1] = 1.0; n[2] = 1.0;
			array_1d<double,3> auxn;
			//loop on the neighbouring contact elements and try to calculate h and n basing on them
			for(unsigned int i=0;i<4; i++)
			{
				
				Geometry< Node<3> >& other_geom = neigb_el[i].GetGeometry();	
				if(neigb_el[i].Id() != this->Id())
				{
				    WeakPointerVector< Element >& other_neigb_el = neigb_el[i].GetValue(NEIGHBOUR_ELEMENTS);

				    int numb_other_ext_faces = 0;
				    int other_reference_face = -1;
				    for(unsigned int kkk=0;kkk<4; kkk++)
				    {
					    if(other_neigb_el[kkk].Id() == neigb_el[i].Id()) //missing neighbour
					    {
						    numb_other_ext_faces++;
						    other_reference_face=kkk;
					    }
				    }

				    if (numb_other_ext_faces == 1) //in this case i have a reference face
				    {
					    boost::numeric::ublas::bounded_matrix<double, 4, 3 > otherDN_DX;
					    array_1d<double,4> otherN;
					    double othervol;
					    GeometryUtils::CalculateGeometryData(other_geom,otherDN_DX,otherN,othervol);
					    double aaa = norm_2(row(otherDN_DX,other_reference_face));
					    double auxh = 1.0/aaa;
					    noalias(auxn) = row(otherDN_DX,other_reference_face);
					    auxn /= aaa;
					    bool image_inside = Check_image_inside_face(auxn, other_reference_face, other_geom,true); 
					    if(image_inside==true && auxh < h)
					    {
// std::cout << " h calc" << Id() << " " << neigb_el[i].Id() << std::endl;
// 				    KRATOS_WATCH(h);

						h = auxh;
						noalias(n) = auxn;	
/*KRATOS_WATCH("numb_ext_faces == else");
KRATOS_WATCH(h);
KRATOS_WATCH(n);*/						
//  						if(h<l0){
// 				                     KRATOS_WATCH("numb_ext_faces == else");
// 						     KRATOS_WATCH(h);
// 						}

					    }

					    is_contact = true;
					    ss3 = 500;
				    }
				 }

			}
// 	      // look for strange elements
// 	      double  primer_flag = 3.0;
// 	      double second_flag = 5.0;
// 	      int green = 0;
// 	      int red = 0;
// 
// 		boost::numeric::ublas::bounded_matrix<double, 4, 3 > ordered_points;
//               double flag;
// 		  for(int ii = 0; ii<4; ++ii){
// 			flag = geom[ii].FastGetSolutionStepValue(FLAG_VARIABLE);
// 			if(flag == primer_flag)
// 			      {
// 				    ordered_points(green,0 )= geom[ii].X();
// 				    ordered_points(green,1 )= geom[ii].Y();
// 				    ordered_points(green,2 )= geom[ii].Z();
// 				green++;
// 			      }
// 			if(flag == second_flag)
// 			      {
// 
// 				int index = 3-red;
// 				    ordered_points(index,0 )= geom[ii].X();
// 				    ordered_points(index,1 )= geom[ii].Y();
// 				    ordered_points(index,2 )= geom[ii].Z();
// 				red++;
// 			      }
// 
// 		  }
// 
// 
// 		if(green == 2 && red == 2 )
// 		  {
// 		      is_contact = true;
// 		      CalculateMinDistanceAndNormal(h,n,ordered_points);
// //KRATOS_WATCH(">>>>>>>>>>>>>>>>>>>>>>>>> IT IS AN STRANGE CONTACT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<");
// 
// 
// 		  }


		}
	
// 		if(ss3==500 && h<l0){
// 		  KRATOS_WATCH(is_contact);
// 		  KRATOS_WATCH(h);
// 		}
		
	 		
// 		if( h<l0){
// 		  KRATOS_WATCH("Good H but not contact");
// 		 KRATOS_WATCH(h); 
// 			 KRATOS_WATCH(is_contact);
// 		}
	 }
// 	    int my_flag = 0;
//           CheckIsContactMaster(my_flag);
// 	  if(my_flag ==1)
// 	    KRATOS_WATCH("BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBVVVVVVVVVVVVVVVVVVVVVVVVVVVNNNNNNNNNNNNNNNNNNNNNNN");


		if(is_contact == true) //note that if a reference face is not encountered nothing is done
		{
			//calculate h, n
			/*Geometry< Node<3> >& geom = GetGeometry();
			double Volume,norm_n;				
			GeometryUtils::CalculateGeometryData(geom,DN_DX,N,Volume);

			array_1d<double,3> n;
			noalias(n) = row(DN_DX,reference_face);
			norm_n = norm_2(n);
			double h=1.0/norm_n;
			n/=norm_n;*/
// 	KRATOS_WATCH(h);		
// 	KRATOS_WATCH(l0);
	//KRATOS_WATCH(GetProperties()[THICKNESS]);
	//KRATOS_WATCH(GetProperties()[DENSITY]);
	
	if(h<0.0)
             KRATOS_ERROR(std::logic_error,"NEGATIVE HH","")
             
			if(h>l0) //in this case it is not in contact
			{
				if(CalculateResidualVectorFlag==true) noalias(rRightHandSideVector) = ZeroVector(12);
				if(CalculateStiffnessMatrixFlag==true) noalias(rLeftHandSideMatrix) = ZeroMatrix(12,12);
			}
			else
			{
			  //bool image_inside = Check_image_inside_face(n, reference_face);  
			//   if(image_inside)
			      //{
// KRATOS_WATCH(">>>>>>>>>>>>>>>>>>>>>><<<one contact element is detected<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<");
// KRATOS_WATCH(h);
// KRATOS_WATCH(n);
// int flag = 0;
// CheckIsContactMaster(flag);
// KRATOS_WATCH(flag);
 				this->GetValue(IS_CONTACT_MASTER) = 120;
/*				KRATOS_WATCH(h);*/
				mcontact_is_active = true;
				mpenetration = l0 - h;
// if(mcontact_is_active == true)
//    KRATOS_WATCH("inside pfem_contact_element_3D");
				//calculate old iteration height and Alpha factor
				//double H_zero = l0; // It means even in teh old iteration it was in contact
				//CalculateOldIterationContactHeight(H_zero,reference_face);
				
				/*double Alpha_factor;
				if(H_zero < l0)
				    Alpha_factor = 1.0;// it is allready in contact
				else
				    Alpha_factor = (l0 - h)/(H_zero - h);*/

//KRATOS_WATCH(H_zero);
				//calculate 1d strain
				double eps= 0.5*(h*h-l0*l0)/(l0*l0);
// KRATOS_WATCH(eps);
				//calculate effective strain 
				array_1d<double,6> strain;
				noalias(strain)  = VoigtTensorComponents(n,n);
				strain *= eps;
// 	KRATOS_WATCH(h);		
// KRATOS_WATCH(strain);
				//compute elastic matrix
				//setting up material matrix
				const double E = GetProperties()[YOUNG_MODULUS];
				const double NU = GetProperties()[POISSON_RATIO];
				const double c1 = E / ((1.00+NU)*(1.0-2.0*NU));
				const double c2 = c1 * (1.0-NU);
				const double c3 = c1 * NU;
				const double c4 = c1 * 0.5 * (1.0 - 2.0*NU);
				boost::numeric::ublas::bounded_matrix<double, 6, 6 > C;
				boost::numeric::ublas::bounded_matrix<double, 6, 12 > B;
				boost::numeric::ublas::bounded_matrix<double, 6, 12 > aux;

				//filling material matrix
				C(0,0) = c2;    C(0,1) = c3;    C(0,2) = c3;    C(0,3) = 0.0;   C(0,4) = 0.0;   C(0,5) = 0.0;
				C(1,0) = c3;    C(1,1) = c2;    C(1,2) = c3;    C(1,3) = 0.0;   C(1,4) = 0.0;   C(1,5) = 0.0;
				C(2,0) = c3;    C(2,1) = c3;    C(2,2) = c2;    C(2,3) = 0.0;   C(2,4) = 0.0;   C(2,5) = 0.0;
				C(3,0) = 0.0;   C(3,1) = 0.0;   C(3,2) = 0.0;   C(3,3) = c4;    C(3,4) = 0.0;   C(3,5) = 0.0;
				C(4,0) = 0.0;   C(4,1) = 0.0;   C(4,2) = 0.0;   C(4,3) = 0.0;   C(4,4) = c4;    C(4,5) = 0.0;
				C(5,0) = 0.0;   C(5,1) = 0.0;   C(5,2) = 0.0;   C(5,3) = 0.0;   C(5,4) = 0.0;   C(5,5) = c4;

				//computer stresses
				array_1d<double,6> stress;
				noalias(stress) = prod(C,strain);
// KRATOS_WATCH(stress);				
				//calculate B
				for (unsigned int i=0;i<number_of_nodes;i++)
				{
					unsigned int index = dim*i;
					B(0,index+0)=DN_DX(i,0);	B(0,index+1)=0.0;		B(0,index+2)=0.0;
					B(1,index+0)=0.0;		B(1,index+1)=DN_DX(i,1);	B(1,index+2)=0.0;
					B(2,index+0)=0.0;		B(2,index+1)=0.0;		B(2,index+2)=DN_DX(i,2);
					B(3,index+0)=DN_DX(i,1);	B(3,index+1)=DN_DX(i,0);	B(3,index+2)=0.0;
					B(4,index+0)=0.0;		B(4,index+1)=DN_DX(i,2);	B(4,index+2)=DN_DX(i,1);
					B(5,index+0)=DN_DX(i,2);	B(5,index+1)=0.0;		B(5,index+2)=DN_DX(i,0);
				}
// KRATOS_WATCH(Volume);	
//KRATOS_WATCH(B);		
				//calculate RHS and LHS;
				if(CalculateResidualVectorFlag==true) 
					noalias(rRightHandSideVector) -= Volume*prod(trans(B),stress);
//  KRATOS_WATCH(rRightHandSideVector);					
				if(CalculateStiffnessMatrixFlag==true)
				{
					//KRATOS_WATCH(Alpha_factor);
					noalias(aux) = prod(C,B);
					noalias(rLeftHandSideMatrix) = Volume*prod(trans(B),aux);
					//noalias(rLeftHandSideMatrix) *= Volume;
				}
// KRATOS_WATCH(rLeftHandSideMatrix);
			 /* }
			  else
			      {
				if(CalculateResidualVectorFlag==true) noalias(rRightHandSideVector) = ZeroVector(12);
				if(CalculateStiffnessMatrixFlag==true) noalias(rLeftHandSideMatrix) = ZeroMatrix(12,12);
			      }*/
/*double auxL = 0.0;
double auxR = 0.0;
for(unsigned int i=0; i<12; i++)
{
    auxR += rRightHandSideVector[i];
    if(CalculateStiffnessMatrixFlag==true)	
      for(unsigned int j=0; j<12; j++)
	  auxL += rLeftHandSideMatrix(i,j);
}

std::cout << Id() << " " << auxR << " " << auxL << std::endl;*/
			}

			

			


			
		}

		

		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void PfemContactElement3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
	  KRATOS_TRY
		//calculation flags
		bool CalculateStiffnessMatrixFlag = false;
		bool CalculateResidualVectorFlag = true;
		MatrixType temp = Matrix();
			
	  
		CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
	  KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void PfemContactElement3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
	      KRATOS_TRY
		//calculation flags
		bool CalculateStiffnessMatrixFlag = true;
		bool CalculateResidualVectorFlag = true;
		
		CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
            KRATOS_CATCH("")    

	}

	

        ////************************************************************************************
	////************************************************************************************
	
        void PfemContactElement3D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)	
 	    {
	      KRATOS_TRY

          
	      KRATOS_CATCH("")
 	    }

	////************************************************************************************
	////************************************************************************************
	void PfemContactElement3D::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
	  KRATOS_TRY
	  KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void PfemContactElement3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
	    KRATOS_TRY
		int number_of_nodes = GetGeometry().size();
		int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int dim2 = number_of_nodes*dim;
		if(rResult.size() != dim2)
			rResult.resize(dim2,false);

		for (int i=0;i<number_of_nodes;i++)
		{
			int index = i*dim;
			rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
			rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
			if(dim == 3)
				rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
		}
	    KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void PfemContactElement3D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		ElementalDofList.resize(0);

		for (unsigned int i=0;i<GetGeometry().size();i++)
		{
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
			if(GetGeometry().WorkingSpaceDimension() == 3)
                        {
				ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
                        }
		}
		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void PfemContactElement3D::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
			//lumped
			unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		unsigned int NumberOfNodes = GetGeometry().size();
		unsigned int MatSize = dimension * NumberOfNodes;

		if(rMassMatrix.size1() != MatSize)
			rMassMatrix.resize(MatSize,MatSize,false);

		noalias(rMassMatrix) = ZeroMatrix(MatSize,MatSize);


		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void PfemContactElement3D::DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().size();
		unsigned int dim = GetGeometry().WorkingSpaceDimension();

		boost::numeric::ublas::bounded_matrix<double, 4, 3 > DN_DX;
		array_1d<double,4> N;

		
		//resizing as needed the LHS
		unsigned int MatSize=number_of_nodes*dim;

		if(rDampMatrix.size1() != MatSize)
			rDampMatrix.resize(MatSize,MatSize,false);

		noalias(rDampMatrix)= ZeroMatrix(MatSize,MatSize);
      


// 		if(mcontact_is_active == true)
// 		  {
// 
// 
// 		      Geometry< Node<3> >& geom = GetGeometry();
// 		      double Volume;				
// 		      GeometryUtils::CalculateGeometryData(geom,DN_DX,N,Volume);
// 
// 		      double lump_mass_fac = Volume * 0.25;
// 		      const double density=GetProperties()[DENSITY];
// 
// 		      int nodes_number = 4;
// 		      int dof = 3;
// 		  for (int nd = 0; nd < nodes_number; nd++) {
// 			  int row = nd * dof ;
// 			    for (int jj = 0; jj < dof; jj++)
// 			  rDampMatrix(row + jj, row + jj) += 100.0*density * lump_mass_fac;
// 			  }
// 
// 
// 
// 		   }

		KRATOS_CATCH("")
	}
        
        




	//************************************************************************************
	//************************************************************************************
	void PfemContactElement3D::GetValuesVector(Vector& values, int Step)
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize = number_of_nodes*dim;
		if(values.size() != MatSize)	values.resize(MatSize,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = i*dim;
			values[index] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_X,Step);
			values[index + 1] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Y,Step);
			if(dim == 3)
				values[index + 2] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Z,Step);
		}
	}
	
	
	//************************************************************************************
	//************************************************************************************
	  void PfemContactElement3D::GetFirstDerivativesVector(Vector& values, int Step)
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize = number_of_nodes*dim;
		if(values.size() != MatSize)   values.resize(MatSize,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = i*dim;
			values[index] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_X,Step);
			values[index + 1] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_Y,Step);
			if(dim == 3)
				values[index + 2] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_Z,Step);
		}
	}

	//************************************************************************************
	//************************************************************************************
	  void PfemContactElement3D::GetSecondDerivativesVector(Vector& values, int Step)
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize = number_of_nodes*dim;
//KRATOS_WATCH(dim);
//KRATOS_WATCH(number_of_nodes);
		if(values.size() != MatSize) values.resize(MatSize,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = i*dim;
			values[index] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_X,Step);
			values[index + 1] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Y,Step);
			if(dim == 3)
				values[index + 2] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Z,Step);
		}
	}
	//************************************************************************************
	//************************************************************************************
	//auxiliary function needed in the calculation of output stresses
	inline array_1d<double,6> PfemContactElement3D::VoigtTensorComponents(
		array_1d<double,3>& a,
		array_1d<double,3>& b)

	{
		array_1d<double,6> v;

		v[0] = a[0]*b[0];
		v[1] = a[1]*b[1];
		v[2] = a[2]*b[2];
		v[3] = a[0]*b[1];
		v[4] = a[1]*b[2];
		v[5] = a[0]*b[2];

		return v;
	}
	//************************************************************************************
	//************************************************************************************
	 bool PfemContactElement3D::Check_image_inside_face(const array_1d<double,3> n,
							    int reference_face,
							    const Geometry< Node<3> >& geom,
							    bool image_inside)
	{
	    
	    boost::numeric::ublas::bounded_matrix<double, 3, 3 > triangle = ZeroMatrix(3, 3);
	    array_1d<double,3> check_point, edge;
	    check_point(0) = geom[reference_face].X();
	    check_point(1) = geom[reference_face].Y();
	    check_point(2) = geom[reference_face].Z();

	    int kk = 0;
	    for(int ii=0; ii<4; ++ii)
	      {
		  if(ii != reference_face)
		      {
			//row(triangle,kk) = geom(ii);
			triangle(kk,0) = geom[ii].X();
			triangle(kk,1) = geom[ii].Y();
			triangle(kk,2) = geom[ii].Z();
			kk++;
		      }
	      }
	    edge = check_point - row(triangle,0);
// 	    check_point -= row(triangle,0);
	    double image = inner_prod(edge,n);
	    check_point -= image*n;
 
	  if(image_inside)
	    {
		if(same_side(check_point,row(triangle,0),row(triangle,1),row(triangle,2)) && 
		  same_side(check_point,row(triangle,1),row(triangle,0),row(triangle,2)) && 
		  same_side(check_point,row(triangle,2),row(triangle,0),row(triangle,1)))
		      return true;
		else
		      return false;
	    }
	  else
	    {
		//check if the image is inside circumsphere
		array_1d<double,3> center,radi_vec,dist_vec;
		for(int ii=0; ii<3;++ii)
		      center(ii) = 0.33333333333333333333*(triangle(0,ii)+triangle(1,ii)+triangle(2,ii));
		
	      double max_radi=0.0;
		for(int ii=0; ii<3;++ii)
		  {
		    radi_vec = center -  row(triangle,ii);
		    double radi = norm_2(radi_vec);	
		      if(radi>max_radi)
			  max_radi = radi;
		  }


		dist_vec = check_point - center;
		double dist = norm_2(dist_vec);

		if(dist <= 2.0*max_radi)
		  return true;
		else
		  return false;
	    }
	}
	//************************************************************************************
	//************************************************************************************
	 bool PfemContactElement3D::same_side(const array_1d<double,3> p0, const array_1d<double,3> p1,const array_1d<double,3> a,const array_1d<double,3> b)
	{
	     array_1d<double,3> check_line;
	     array_1d<double,3> diff_1;
	     array_1d<double,3> diff_2;

	      check_line = b - a;
	      diff_1 = p1 - a;
	      diff_2 = p0 - a;

	      array_1d<double,3> first_out;
	      array_1d<double,3> second_out;

	     MathUtils<double>::CrossProduct(first_out,check_line , diff_1);
	     MathUtils<double>::CrossProduct(second_out,check_line, diff_2);

	      double val = inner_prod(first_out , second_out);

	      if(val >= 0.0)
		  return true;
	      else
		  return false;
	    
	}
	//************************************************************************************
	//***********************************************************************************
      	  void PfemContactElement3D::CalculateOldIterationContactHeight(double& H_zero, const int reference_face)
	  {
	    Geometry< Node<3> >& geom = GetGeometry();
	    boost::numeric::ublas::bounded_matrix<double, 3, 3 > triangle = ZeroMatrix(3, 3);
	    array_1d<double,3> check_point;
	    check_point(0) = geom[reference_face].X0();
	    check_point(1) = geom[reference_face].Y0();
	    check_point(2) = geom[reference_face].Z0();

	    check_point(0) += geom[reference_face].FastGetSolutionStepValue(DISPLACEMENT_X,1);
	    check_point(1) += geom[reference_face].FastGetSolutionStepValue(DISPLACEMENT_Y,1);
	    check_point(2) += geom[reference_face].FastGetSolutionStepValue(DISPLACEMENT_Z,1);
				//(i)->X() = (i)->X0() + i->GetSolutionStepValue(DISPLACEMENT_X);

	    int kk = 0;
	    for(int ii=0; ii<4; ++ii)
	      {
		  if(ii != reference_face)
		      {
			//row(triangle,kk) = geom(ii);
			triangle(kk,0) = geom[ii].X0();
			triangle(kk,0) += geom[ii].FastGetSolutionStepValue(DISPLACEMENT_X,1);

			triangle(kk,1) = geom[ii].Y0();
			triangle(kk,1) += geom[ii].FastGetSolutionStepValue(DISPLACEMENT_Y,1);

			triangle(kk,2) = geom[ii].Z0();
			triangle(kk,2) += geom[ii].FastGetSolutionStepValue(DISPLACEMENT_Z,1);

			kk++;
		      }
	      }
	    
	    //calculate normal
	     array_1d<double,3> diff_1;
	     array_1d<double,3> diff_2;

	     diff_1 = row(triangle,1) - row(triangle,0);
	     diff_2 = row(triangle,2) - row(triangle,0);
	     check_point -=  row(triangle,0);

	     array_1d<double,3> normal;
	     MathUtils<double>::CrossProduct(normal,diff_1 , diff_2);

	     double area = norm_2(normal);
	     normal /= area;

	     double val = inner_prod(normal , check_point);

	     H_zero = fabs(val);

	  }
	//************************************************************************************
	//***********************************************************************************
         void PfemContactElement3D::CalculateMinDistanceAndNormal(double& h,
								  array_1d<double,3>& n,
								  const boost::numeric::ublas::bounded_matrix<double, 4, 3 > ordered_points)
	{

	    array_1d<double,3> u,v,w;
	
	      u = row(ordered_points,1) - row(ordered_points,0);
	      v = row(ordered_points,3) - row(ordered_points,2);
	      w = row(ordered_points,0) - row(ordered_points,2);

	      double a,b,c,d,e;
	      a = inner_prod(u,u);
	      b = inner_prod(u,v);
	      c = inner_prod(v,v);
	      d = inner_prod(u,w);
	      e = inner_prod(v,w);

	      double DD = a*c - b*b;
	      double sc,sN,sD =DD;
	      double tc, tN, tD =DD;

	      if(DD < 0.000000001){
		sN = 0.0;        // force using point P0 on segment S1
		sD = 1.0;        // to prevent possible division by 0.0 later
		tN = e;
		tD = c;
		  }
	      else
		  {
		sN = (b*e - c*d);
		tN = (a*e - b*d);
		if (sN < 0.0) {       // sc < 0 => the s=0 edge is visible
		    sN = 0.0;
		    tN = e;
		    tD = c;
		   }
	      else if (sN > sD) {  // sc > 1 => the s=1 edge is visible
		    sN = sD;
		    tN = e + b;
		     tD = c;
		    }
		  }

	      if (tN < 0.0) {           // tc < 0 => the t=0 edge is visible
		   tN = 0.0;
		  // recompute sc for this edge
		  if (-d < 0.0)
		    sN = 0.0;
		  else if (-d > a)
		    sN = sD;
		  else {
		    sN = -d;
		    sD = a;
		    }
		  }
	      else if (tN > tD) {      // tc > 1 => the t=1 edge is visible
		     tN = tD;
		    // recompute sc for this edge
		    if ((-d + b) < 0.0)
			  sN = 0;
		    else if ((-d + b) > a)
			  sN = sD;
		    else {
			  sN = (-d + b);
			  sD = a;
			  }
		  }
	  
	  // finally do the division to get sc and tc
	      sc = (abs(sN) < 0.000000001 ? 0.0 : sN / sD);
	      tc = (abs(tN) < 0.000000001 ? 0.0 : tN / tD);


		    n = w + (sc * u) - (tc * v);
	
		    h = norm_2(n);
		    n/=h;

 
/*	array_1d<double,3> mid_one,mid_two;
      mid_one = 0.5*(row(ordered_points,1) + row(ordered_points,0));
      mid_two = 0.5*(row(ordered_points,3) + row(ordered_points,2));

      n = mid_two - mid_one;
      h = norm_2(n);
      n/=h;*/
	      


	}




	//************************************************************************************
	//***********************************************************************************
         void PfemContactElement3D::DetectContact(Geometry< Node<3> >& geom, 
		  boost::numeric::ublas::bounded_matrix<double, 4, 3 > DN_DX, 
		  const unsigned int single_node_index,
		  array_1d<double,3>& n,
		  double& h)
		{
		    KRATOS_TRY
// 	KRATOS_WATCH("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<11111111");	    
		    array_1d<unsigned int,3> local_indices;
		    int jj=0;
		    for(unsigned int i=0; i<4; i++)
		      if(i != single_node_index)
			{
			  local_indices[jj] = i;
			  jj++;
			}
// KRATOS_WATCH(local_indices);
		    //coordinates of the nodes in the base face
		    const array_1d<double,3>& p0 = geom[local_indices[0]].Coordinates();
		    const array_1d<double,3>& p1 = geom[local_indices[1]].Coordinates();
		    const array_1d<double,3>& p2 = geom[local_indices[2]].Coordinates();
		    
		    //coordinates of the single node
		    const array_1d<double,3>& s2 = geom[single_node_index].Coordinates();
// 	KRATOS_WATCH("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<22222222");		    
		    //compute the direction and height of the tetrahedra
		    noalias(n) = row(DN_DX,single_node_index);
		    h = 1.0 / norm_2(n);
		    n *= h; //normalizing the direction    
		    
		    //check if it fills inside
		    array_1d<double,3> proj_in_face = s2 - h*n; //check the sign!!!!!
		    const array_1d<double,3> x01 = p1 - p0;
		    const array_1d<double,3> x12 = p2 - p1;
		    const array_1d<double,3> x20 = p0 - p2;
		    
		    const array_1d<double,3> x03 = proj_in_face - p0;
		    const array_1d<double,3> x13 = proj_in_face - p1;
		    const array_1d<double,3> x23 = proj_in_face - p2;
// 	KRATOS_WATCH("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<333333333");		    
		    array_1d<double,3> nel;
		    MathUtils<double>::CrossProduct(nel,x01,x20);
		    double a_el = -norm_2(nel);
		    nel /= a_el;
		
		    array_1d<double,3> aux;
		    MathUtils<double>::CrossProduct(aux,x12,x13);
		    double A0 = 0.5*inner_prod(aux,nel);
		    
		    MathUtils<double>::CrossProduct(aux,x20,x23);
		    double A1 = 0.5*inner_prod(aux,nel);

		    MathUtils<double>::CrossProduct(aux,x01,x03);
		    double A2 = 0.5*inner_prod(aux,nel);
// 	KRATOS_WATCH("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<444444444");			    
		    if(A0 < 0.0 || A1<0.0 || A2 < 0.0) 
		    {
			//the projection is outside of the base triangle ... we can not use it
			if(A0 < 0.0)
			{
			  if(A1 <= 0.0) //closest point is 2
				  noalias(n) = s2 - p2;
			  else if(A2 <= 0.0) //closest point is 1
				  noalias(n) = s2 - p1;
			  else //compute distance from segment 12
			  {
			    noalias(aux) = s2 - p1;
			    double l12 = norm_2(x12);
			    double aaa = inner_prod(x12,aux) / l12;
			    noalias(n) = aux - (aaa/l12)*x12;	    
			  }
			}
			
			else if(A1 < 0.0)
			{
			  if(A2 <= 0.0) //closest point is 0
				  noalias(n) = s2 - p0;
			  else if(A0 <= 0.0) //closest point is 2
				  noalias(n) = s2 - p2;
			  else //compute distance from segment 20
			  {
			    noalias(aux) = s2 - p2;
			    double l20 = norm_2(x20);
			    double aaa = inner_prod(x20,aux) / l20;
			    noalias(n) = aux - (aaa/l20)*x20;	    
			  }
			}

			else if(A2 < 0.0)
			{
			  if(A0 <= 0.0) //closest point is 1
				  noalias(n) = s2 - p1;
			  else if(A1 <= 0.0) //closest point is 0
				  noalias(n) = s2 - p0;
			  else //compute distance from segment 10
			  {
			    noalias(aux) = s2 - p1;
			    double l10 = norm_2(x01);
			    double aaa = inner_prod(x01,aux) / l10;
			    noalias(n) = aux - (aaa/l10)*x01;	    
			  }  
			}

			h = norm_2(n);
			n /= h;
		      }
//     	KRATOS_WATCH("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<555555555");	
    
		      KRATOS_CATCH("");
		    }
	//************************************************************************************
	//***********************************************************************************
      void PfemContactElement3D::CheckIsContactMaster(int& flag)
	  {
		    KRATOS_TRY
		boost::numeric::ublas::bounded_matrix<double, 4, 3 > DN_DX;
		array_1d<double,4> N;

// 		const unsigned int number_of_nodes = GetGeometry().size();
// 		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
// 		unsigned int MatSize=number_of_nodes*dim;
		
		const double l0 = GetProperties()[THICKNESS];

		//count external faces
		unsigned int numb_ext_faces = 0;
		bool is_contact = false;
		flag = 0;
	       Geometry< Node<3> >& geom = GetGeometry();
		array_1d<double,3> n;
		double h,Volume;
		GeometryUtils::CalculateGeometryData(geom,DN_DX,N,Volume);

		//find reference face. If a reference face is not found, then do nothing
		int reference_face = -1;
		WeakPointerVector< Element >& neigb_el = this->GetValue(NEIGHBOUR_ELEMENTS);
		for(unsigned int i=0;i<4; i++)
		{
		      if(neigb_el[i].Id() == this->Id()) //missing neighbour
			{
				numb_ext_faces++;
				reference_face=i;
			}
		}

		if(numb_ext_faces > 1) //this is most probably a sliver
		{
			h = 1e6;
			n[0] = 1.0; n[1] = 1.0; n[2] = 1.0;
			array_1d<double,3> auxn;
			is_contact = false;
			
			for(unsigned int i=0;i<4; i++)
			{
			  if(neigb_el[i].Id() == this->Id()) //missing neighbour
			  {
				reference_face = i;
				double aaa = norm_2(row(DN_DX,reference_face));
				h = 1.0/aaa;
				noalias(n) = row(DN_DX,reference_face);
				n /= aaa;

				//accept only if projection is inside
				bool image_inside = Check_image_inside_face(n, reference_face, GetGeometry(),true); 
				  if(image_inside)
				  {
					    is_contact = true;   

					    break;
				  }    
			    }
			}
		}
		else if (numb_ext_faces == 1) //in this case i have a reference face
		{
			double aaa = norm_2(row(DN_DX,reference_face));
			h = 1.0/aaa;
			noalias(n) = row(DN_DX,reference_face);
			n /= aaa;

			//accept only if projection is inside
// 			bool image_inside = Check_image_inside_face(n, reference_face, GetGeometry(),false); 
// KRATOS_WATCH("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<before detect");
			DetectContact(GetGeometry(),DN_DX,reference_face, n, h);
			bool image_inside = true;

			   if(image_inside)
				  {
				    is_contact = true;

				  }
			   else
				  {
				    is_contact = false;
				    
				  }
		}
		else //case B - one edge on each side
		{
			h = 1e6;
			n[0] = 1.0; n[1] = 1.0; n[2] = 1.0;
			array_1d<double,3> auxn;
			//loop on the neighbouring contact elements and try to calculate h and n basing on them
			for(unsigned int i=0;i<4; i++)
			{
				
				Geometry< Node<3> >& other_geom = neigb_el[i].GetGeometry();	
				if(neigb_el[i].Id() != this->Id())
				{
				    WeakPointerVector< Element >& other_neigb_el = neigb_el[i].GetValue(NEIGHBOUR_ELEMENTS);

				    int numb_other_ext_faces = 0;
				    int other_reference_face = -1;
				    for(unsigned int kkk=0;kkk<4; kkk++)
				    {
					    if(other_neigb_el[kkk].Id() == neigb_el[i].Id()) //missing neighbour
					    {
						    numb_other_ext_faces++;
						    other_reference_face=kkk;
					    }
				    }

				    if (numb_other_ext_faces == 1) //in this case i have a reference face
				    {
					    boost::numeric::ublas::bounded_matrix<double, 4, 3 > otherDN_DX;
					    array_1d<double,4> otherN;
					    double othervol;
					    GeometryUtils::CalculateGeometryData(other_geom,otherDN_DX,otherN,othervol);
					    double aaa = norm_2(row(otherDN_DX,other_reference_face));
					    double auxh = 1.0/aaa;
					    noalias(auxn) = row(otherDN_DX,other_reference_face);
					    auxn /= aaa;

					    bool image_inside = Check_image_inside_face(auxn, other_reference_face, other_geom,true); 

					    if(image_inside==true && auxh < h)
					    {
// std::cout << " h calc" << Id() << " " << neigb_el[i].Id() << std::endl;
						h = auxh;
						noalias(n) = auxn;					    
					    }
					    is_contact = true;

				    }
				 }

			}
               }

	      if(is_contact ==true && h < l0 )
	      {
                      flag = 1;
// 		      KRATOS_WATCH(h);
// 
// 		      KRATOS_WATCH(numb_ext_faces);
// 		      KRATOS_WATCH(is_contact);
	      }
		      KRATOS_CATCH("");
	  }
	//************************************************************************************
	//************************************************************************************
        void PfemContactElement3D::Calculate( const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
       {
//KRATOS_WATCH("mcontact_is_active");
          int flag = 0;
           CheckIsContactMaster(flag);


	  double sound_velocity = 4700.0;
	  //double hp = (this)->mpenetration;
	  const double hp = GetProperties()[THICKNESS];
	  if( flag == 1 )
	     {
	      Output = 1.0 * hp / sound_velocity;
 		KRATOS_WATCH(" XXXXXXXXXXXXXX CONTACT DT XXXXXXXXXXXXX");
// 		KRATOS_WATCH(Output);
	     }
	  else 
	      Output = 100.0;
	 
        }
	//************************************************************************************
	//************************************************************************************                
          void PfemContactElement3D::InitializeNonLinearIteration(ProcessInfo& CurrentProcessInfo){
	  

  
	  }
	  
	//************************************************************************************
	//************************************************************************************                
          void PfemContactElement3D::FlagVariableCheckForNonSuitableElements(double& accepted){
	    
// 	    double first = GetGeometry()[0].FastGetSolutionStepValue(FLAG_VARIABLE);
// 	    double suitable = 0.0;
// 	    for(int ii = 1; ii<4; ++ii){
// 	      double other =  GetGeometry()[ii].FastGetSolutionStepValue(FLAG_VARIABLE);
// 	      if(first != other)
// 		suitable = 1.0;
// 	      
// 	      if(other < 0.0)
// 		suitable = 0.0;
// 	    }
// 	    
// 	    if(suitable == 1.0 and first>0.0){
// 	         accepted = 1.0;
// 		this->GetValue(IS_CONTACT_MASTER) = 300;
// 	    }
// 	    else
// 	    {
// 	         accepted = 0.0;
// 
// 	    }
	    
	    for(int ii = 0; ii<4; ++ii){
	      double neg =  GetGeometry()[ii].FastGetSolutionStepValue(FLAG_VARIABLE);
	      if(neg < 0.0)
	      {
// 		this->GetValue(IS_CONTACT_MASTER) = 5000;
		accepted = 0.0;
	      }
	    }	    
	    
	    
	  }
} // Namespace Kratos