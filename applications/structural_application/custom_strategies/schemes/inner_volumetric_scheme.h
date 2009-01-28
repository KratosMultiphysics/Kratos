/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

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
 
/* *********************************************************   
*          
*   Last Modified by:    $Author: kazem $
*   Date:                $Date: 2009-01-16 10:37:04 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/


#if !defined(KRATOS_INNER_VOLUMETRIC_STATIC_SCHEME )
#define  KRATOS_INNER_VOLUMETRIC_STATIC_SCHEME


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "custom_conditions/pointforce2D.h"
#include "custom_conditions/pointforce3D.h"
#include "utilities/geometry_utilities.h" 
//#include "geometries/geometry.h"
//#include "includes/constitutive_law.h"

namespace Kratos
{

	/**@name Kratos Globals */
	/*@{ */


	/*@} */
	/**@name Type Definitions */       
	/*@{ */

	/*@} */


	/**@name  Enum's */       
	/*@{ */


	/*@} */
	/**@name  Functions */       
	/*@{ */



	/*@} */
	/**@name Kratos Classes */
	/*@{ */

	/** Short class definition.

	This class provides the implementation of the basic tasks that are needed by the solution strategy.
	It is intended to be the place for tailoring the solution strategies to problem specific tasks.

	Detail class definition.

	\URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

	\URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

	\URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

	\URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


	\URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

	\URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

	\URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

	\URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


	*/
	template<const int TDim, class TSparseSpace,
	class TDenseSpace //= DenseSpace<double>
	>
	class InnerVolumetricScheme : public Scheme<TSparseSpace,TDenseSpace>
	{

	public:
		/**@name Type Definitions */       
		/*@{ */

		//typedef boost::shared_ptr< ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace> > Pointer;
		KRATOS_CLASS_POINTER_DEFINITION( InnerVolumetricScheme);

		typedef Scheme<TSparseSpace,TDenseSpace> BaseType;

		typedef typename BaseType::TDataType TDataType;

		typedef typename BaseType::DofsArrayType DofsArrayType;

		typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

		typedef typename BaseType::TSystemVectorType TSystemVectorType;

		typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
		typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
		//typedef Geometry<TPointType> GeometryType;

		/*@} */
		/**@name Life Cycle 
		*/    
		/*@{ */

		/** Constructor.
		*/
		InnerVolumetricScheme()
			: Scheme<TSparseSpace,TDenseSpace>()
		{}

		/** Destructor.
		*/
		virtual ~InnerVolumetricScheme(){}


		/*@} */
		/**@name Operators 
		*/  
		/*@{ */

		/** 
		Performing the update of the solution.
		*/


		//void Predict(
		//	const String& ElementGroupName,
		//	DofsArrayType& rDofSet,
		//	TSystemMatrixType& A,
		//	TSystemVectorType& Dx,
		//	TSystemVectorType& b,
		//	ProcessInfo& CurrentProcessInfo
		//	) 
		//{
		//	double CurrentTime = CurrentProcessInfo.GetCurrentTime();
		//	double DeltaTime = CurrentProcessInfo.GetDeltaTime();
		//	double OldTime = CurrentTime - DeltaTime;
		//	
		//	int i;
		//	typename DofsArrayType::iterator it2;
		//	
		//	//predicting variables
		//	for (it2=rDofSet.begin();it2 != rDofSet.end(); ++it2)
		//	{
		//		// N.B. fixed values are not predicted!!
		//		if ( !(*it2)->IsFixed()  )
		//		{
		//			const Dof& X = *(*it2);		
		//			
		//			mpModel->Value(*it2) =  mpModel->Value(X.GetVariable(), X, OldTime); 
		//		}
		//	}		  
		//	
		//}		
		//***************************************************************************
		//**this function ganarates a new " point condition" called "FORCE" to save volumetric force 
		// calculated in every patch

		void Initialize( 
			ModelPart& r_model_part
			)
		{ 
			Properties::Pointer p_properties = r_model_part.pGetProperties(1);
			unsigned int condition_id = r_model_part.NumberOfConditions() + 1;
			Condition::NodesArrayType nodes_array(1);
			KRATOS_TRY
			for(typename ModelPart::NodesContainerType::iterator ind = r_model_part.NodesBegin(); ind != r_model_part.NodesEnd();ind++)
	        	 {
//			   nodes_array[0] = *(ind.base()); //verify we are not making copies here!!
			   Condition::NodesArrayType nodes_array;
			   nodes_array.push_back( *(ind.base())  );

			   Geometry< Node<3> >::Pointer cond = 
						Geometry< Node<3> >::Pointer(new Geometry< Node<3> >(nodes_array) );
			   Condition::Pointer p_cond; 
			   if(TDim == 2)
			   	p_cond = Condition::Pointer(new PointForce2D(condition_id++, cond, p_properties) );
			   else
			   	p_cond = Condition::Pointer(new PointForce3D(condition_id++, cond, p_properties) );

			   r_model_part.AddCondition(p_cond);
//			   r_model_part.AddCondition(PointForce3D(condition_id++, nodes_array, p_properties));
//			   for(typename ModelPart::NodesContainerType::iterator ind = r_model_part.NodesBegin(); ind != r_model_part.NodesEnd();ind++)
//				{
//				ind->FastGetSolutionStepValue(IS_VISITED)=100.00;
//				ind->FastGetSolutionStepValue(NODAL_AREA,2)=-1000000.00;
//				
//				}
	        	 }
			
			KRATOS_CATCH("")
		}
		//****************************************************************************
		virtual void InitializeSolutionStep(
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b
			)
		{
			for(typename ModelPart::NodesContainerType::iterator ind = r_model_part.NodesBegin(); ind != r_model_part.NodesEnd();ind++)
				{
				ind->FastGetSolutionStepValue(IS_VISITED)=100.00;
				ind-> FastGetSolutionStepValue(NODAL_H)=10000.00;
				}
			KRATOS_WATCH("*+*+*+*+*+*InitializeSolutionStep*+*+*+*+*+*");

		}
		//****************************************************************************
		// this function calculate volumetric force and at last fill "b" vector
		virtual void FinalizeNonLinIteration(
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b)
		{
			KRATOS_TRY
			bool DVnorm = true;
			int contor=0;
			int iteration = (r_model_part.GetProcessInfo())[NL_ITERATION_NUMBER];
			//DVnorm = CalculateVolumeChange(r_model_part);
			//DVnorm= false;
// 			 KRATOS_WATCH(r_model_part.Nodes()[3].FastGetSolutionStepValue(FORCE));
			if(DVnorm || iteration<6)
			 {
		           for(typename ModelPart::NodesContainerType::iterator ind = r_model_part.NodesBegin(); ind != r_model_part.NodesEnd();ind++)
	        	      {
				//if(ind->FastGetSolutionStepValue(IS_INTERFACE)==1.0)
// 				KRATOS_WATCH("before");
// 				KRATOS_WATCH(ind->Id());
// 				 KRATOS_WATCH(ind->FastGetSolutionStepValue(FORCE));

				
			         CalculatePatchParameters(ind,r_model_part,contor);
				   //conter++;
				
// 				KRATOS_WATCH("AFTER");
// 				KRATOS_WATCH(ind->Id());
// 				 KRATOS_WATCH(ind->FastGetSolutionStepValue(FORCE));
	        	      }
			   KRATOS_WATCH("**AFTER calculating patch****");


//		 for(typename ModelPart::NodesContainerType::iterator ind = r_model_part.NodesBegin(); ind != r_model_part.NodesEnd();ind++)
//	        	      {
			
//			std::cout << "node Id " << ind->Id()<<" Force X: "<< ind->FastGetSolutionStepValue(FORCE)[0]<< " "<<
//			"Force Y: "<< ind->FastGetSolutionStepValue(FORCE)[1]<<" "<< std::endl;

//			       }
//			   for(typename ModelPart::ElementsContainerType::iterator iel = r_model_part.ElementsBegin(); iel != r_model_part.ElementsEnd();iel++)
//	        	      {
				
				//CalculateCorrectionForce(iel,r_model_part);
//
//	        	      }	
			 }
// 			 KRATOS_WATCH(r_model_part.Nodes()[3].FastGetSolutionStepValue(FORCE));
			//DVnorm = CalculateVolumeChange(r_model_part);
			KRATOS_CATCH("")
		}

		//*****************************************************************************
		void CalculateCauchyStress(ModelPart& model_part)
		{
			KRATOS_TRY
			Matrix temp;
			for(ModelPart::ElementsContainerType::iterator iel = model_part.ElementsBegin(); iel != model_part.ElementsEnd();++iel)
			{
				double dim=(iel->GetGeometry()).WorkingSpaceDimension();
				temp.resize(1,(dim*dim+dim)/2,false);
				temp=ZeroMatrix(1,(dim*dim+dim)/2);
				
				iel->Calculate(CAUCHY_STRESS_TENSOR,temp,model_part.GetProcessInfo());
				iel->SetValue(PK2_STRESS_TENSOR,temp);
			}
			
			KRATOS_CATCH("")
		}




		//*****************************************************************************


		/*@} */
		/**@name Operations */
		/*@{ */


		/*@} */  
		/**@name Access */
		/*@{ */


		/*@} */
		/**@name Inquiry */
		/*@{ */


		/*@} */      
		/**@name Friends */
		/*@{ */


		/*@} */

	protected:
		/**@name Protected static Member Variables */
		/*@{ */


		/*@} */
		/**@name Protected member Variables */
		/*@{ */

		/*@} */
		/**@name Protected Operators*/
		/*@{ */

		/*@} */
		/**@name Protected Operations*/
		/*@{ */


		/*@} */
		/**@name Protected  Access */
		/*@{ */


		/*@} */     
		/**@name Protected Inquiry */
		/*@{ */


		/*@} */   
		/**@name Protected LifeCycle */  
		/*@{ */



		/*@} */    

	private:
	

	//********************************************************************
	void CalculatePatchParameters(ModelPart::NodesContainerType::iterator& inode,ModelPart& r_model_part,int& contor)

 	{

	  WeakPointerVector< Element >& neighbor_els = inode->GetValue(NEIGHBOUR_ELEMENTS);
	  int prsfac = 1;
	  
	for(WeakPointerVector< Element >::iterator elenum = neighbor_els.begin(); elenum!=neighbor_els.end(); elenum++)
	 {
	   elenum->GetValue(IS_VISITED) = 0;
	   prsfac++;
	  }
	 for(WeakPointerVector< Element >::iterator elenum = neighbor_els.begin(); elenum!=neighbor_els.end(); elenum++)
	  {
	    Geometry< Node<3> >& geom = elenum->GetGeometry();

	    //loop over element nodes 
	    for(int elend=0; elend<TDim+1; elend++)
	     {
	      if(geom[elend].Id() != inode->Id())
	       {
	        Element::WeakPointer& nghele = elenum->GetValue(NEIGHBOUR_ELEMENTS)(elend);

		 if(nghele.lock()->GetValue(IS_VISITED) != 10 && nghele.lock()->Id() != elenum->Id())  //check not to repeat	
		  {
	   Matrix auxcoord = ZeroMatrix(TDim, TDim); //save two auxillary coordinates
   
	   Geometry< Node<3> >& nghgeom = nghele.lock()->GetGeometry();
 

	   auxcoord(0,0) = 1.0/3*(geom[0].X()+geom[1].X()+geom[2].X());
	   auxcoord(0,1) = 1.0/3*(geom[0].Y()+geom[1].Y()+geom[2].Y());
	   auxcoord(1,0) = 1.0/3*(nghgeom[0].X()+nghgeom[1].X()+nghgeom[2].X());
	   auxcoord(1,1) = 1.0/3*(nghgeom[0].Y()+nghgeom[1].Y()+nghgeom[2].Y());


	   double area = 0.0;
	   double centerX = inode->X();
	   double centerY = inode->Y();

	   area = centerX*(auxcoord(0,1)-auxcoord(1,1))+auxcoord(1,0)*(centerY-auxcoord(0,1))+auxcoord(0,0)*(auxcoord(1,1)-centerY);

	   //code table
	   int turnflag = 1 ; 
	   if(area < 0.0)
	       {
		turnflag = 0;
		centerX=auxcoord(1,0);			centerY=auxcoord(1,1);
		auxcoord(1,0) = auxcoord(0,0);   	auxcoord(1,1)=auxcoord(0,1);
		auxcoord(0,0)=centerX;			auxcoord(0,1)=centerY;
	       }
   	

	   //calculate volumetric force
	   CalculatePressureForces(auxcoord, inode, elenum, nghele, turnflag, prsfac );

	  } 

	  }
	 }



	 elenum->GetValue(IS_VISITED) = 10;
	 }
//			 for(typename ModelPart::NodesContainerType::iterator ind = r_model_part.NodesBegin(); ind != r_model_part.NodesEnd();ind++)
//	        	      {
			
//			std::cout << "node Id " << ind->Id()<<" Force X: "<< ind->FastGetSolutionStepValue(FORCE)[0]<< " "<<
//			"Force Y: "<< ind->FastGetSolutionStepValue(FORCE)[1]<<" "<< std::endl;

//			       }


	  }

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	 void CalculatePressureForces(Matrix crd,ModelPart::NodesContainerType::iterator& ctrnd, WeakPointerVector< Element >::iterator& host,  Element::WeakPointer& neighbor, int flag,int pressurefac )
	 {
	 int patchsize = TDim*(TDim+1);
	 Matrix volstiffness = ZeroMatrix(patchsize,patchsize);
	 Matrix B = ZeroMatrix(TDim +1, patchsize);
	 Vector H = ZeroVector(patchsize);

	 double Xc = ctrnd->X();
	 double Yc = ctrnd->Y();
	 double area = 0.0;
 
	 area = 0.5 *( Xc*(crd(0,1)-crd(1,1))+crd(1,0)*(Yc-crd(0,1))+crd(0,0)*(crd(1,1)-Yc));
	

	 //fill B for linear element
	  B(0,0) = 0.5*(crd(0,1)-crd(1,1));		B(0,1) = 0.0;		
 	  B(1,0) = 0.0;				B(1,1) = 0.5*(crd(1,0)-crd(0,0));	
	  B(2,0) = 0.5*(crd(1,0)-crd(0,0));	B(2,1) = 0.5*(crd(0,1)-crd(1,1));	

	  B(0,2) = 0.5*(crd(1,1)-Yc);		B(0,3) = 0.0;		
	  B(1,2) = 0.0;				B(1,3) = 0.5*(Xc-crd(1,0));	
	  B(2,2) = 0.5*(Xc-crd(1,0));		B(2,3) = 0.5*(crd(1,1)-Yc);	

	  B(0,4) = 0.5*(Yc-crd(0,1));		B(0,5) = 0.0;	
	  B(1,4) = 0.0;				B(1,5) = 0.5*(crd(0,0)-Xc);
	  B(2,4) = 0.5*(crd(0,0)-Xc);		B(2,5) = 0.5*(Yc-crd(0,1));	

	 //fill H
	 H[0] = B(0,0); H[1] = B(1,1); H[2] = B(0,2); H[3] = B(1,3); H[4] = B(0,4); H[5] = B(1,5);

 
	 //get properties
	 double NU = host->GetProperties()[POISSON_RATIO];
	 double EE = host->GetProperties()[YOUNG_MODULUS];
	 double bulk = EE/(3.0*(1-2*NU));
	 double factor = bulk ;
	 double accelerator = 3.0/(2.0*(1.0 + NU));

	//KRATOS_WATCH(bulk);
	//KRATOS_WATCH(accelerator);


	 //build volstiffness
	 volstiffness = accelerator*factor/area*outer_prod(H, H);



	 Matrix AV;
	 AV = outer_prod(H,H);
	 double denom =0.0;
	 denom = inner_prod(H,H);
	
	 AV/=denom;
	 volstiffness=prod(volstiffness, AV);

	 //fill displacement vector
	 Vector disp=ZeroVector(patchsize);
	 disp[0] = ctrnd->GetSolutionStepValue(DISPLACEMENT_X);
	 disp[1] = ctrnd->GetSolutionStepValue(DISPLACEMENT_Y);

	 Geometry< Node<3> >& hostgeom = host->GetGeometry();
	 Geometry< Node<3> >& nghgeom = neighbor.lock()->GetGeometry();

	 double meanhostX = 0.0;
	 double meanhostY = 0.0;
	 double meannghX = 0.0;
	 double meannghY = 0.0;

	 for(int ii = 0; ii<TDim+1; ii++)
	 {

	   meanhostX += 1.0/3*(hostgeom[ii].GetSolutionStepValue(DISPLACEMENT_X));
	   meanhostY += 1.0/3*(hostgeom[ii].GetSolutionStepValue(DISPLACEMENT_Y));
 	   meannghX += 1.0/3*(nghgeom[ii].GetSolutionStepValue(DISPLACEMENT_X));
	   meannghY += 1.0/3*(nghgeom[ii].GetSolutionStepValue(DISPLACEMENT_Y));

 
	 }

	 if(flag ==1 )
	  {
	   disp[2] = meanhostX;
	   disp[3] = meanhostY;
	   disp[4] = meannghX;
	   disp[5] = meannghY;
	  }
	 else
	  {
	   disp[2] = meannghX ;
	   disp[3] = meannghY;
	   disp[4] = meanhostX;
	   disp[5] = meanhostY;
	  }


 	//calculate pressure force & assembling
	  Vector volforce = ZeroVector(patchsize);
	  volforce = -1*prod( volstiffness, disp);


	 for(int ii=0; ii<TDim; ++ii)
	  ctrnd->FastGetSolutionStepValue(FORCE)[ii] += volforce[ii];


  
	 for(int jj=0;jj<TDim+1; jj++)
	  {
	   int index =TDim; 
	   for(int kk=0;kk<TDim; kk++)
	    if(flag == 1)
	     {

		hostgeom[jj].FastGetSolutionStepValue(FORCE)[kk] += 1.0/3* volforce[index+kk];
	        
	        nghgeom[jj].FastGetSolutionStepValue(FORCE)[kk] += 1.0/3* volforce[index+TDim+kk];

	     }
	    else 
	     {
		nghgeom[jj].FastGetSolutionStepValue(FORCE)[kk] += 1.0/3* volforce[index+kk];
	        
	        hostgeom[jj].FastGetSolutionStepValue(FORCE)[kk] += 1.0/3* volforce[index+TDim+kk];

	     }
		
		

	  }


	 //calculate & update pressure
	 double pressure=0.0;
	 for(int ii=0; ii< patchsize; ii++)
	  if ( H[ii]*H[ii]> 10e-14)
		pressure +=volforce[ii]/(pressurefac*H[ii]); 


	  ctrnd->FastGetSolutionStepValue(PRESSURE) += pressure;

	 }

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	void localAssembleRHS(Vector& Hpatch, const Vector localH, const Vector codetable)
	{
		int size = localH.size();
		
	for ( int ii = 0; ii< size; ii++)
	   {
			
		//KRATOS_WATCH(Hpatch[3]);
		int gcounter = codetable[ii];
		Hpatch(gcounter) += localH[ii];
	   }
	}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	void factoraccelerator(WeakPointerVector< Element >::iterator ielem, double& accelerator, double& factor, const double DT,const double DetF)
	{
		
		int cnt=0;
	
		for(int ii=0;ii<TDim+1;ii++)
			if((ielem->GetGeometry())[ii].FastGetSolutionStepValue(IS_INTERFACE)!=0.0)
				cnt++;
		//KRATOS_WATCH(cnt);
		if(cnt==TDim+1)
			{
			double miu = ielem->GetProperties()[MIU];
			double bulk = ielem->GetProperties()[BULK_MODULUS];
 			//KRATOS_WATCH(miu);
 			//KRATOS_WATCH(bulk);
			double NUhypo = (3*bulk-2*miu)/(2*(3*bulk+miu));
			 factor = bulk ;
			 accelerator = 3.0/(2.0*(1.0 + NUhypo));
			}
		else
			{
			double miu = ielem->GetProperties()[MIU];
			double lambda = ielem->GetProperties()[LAMBDA];
			miu = miu*DT/DetF;
			lambda= lambda*DT/DetF;
			double Ehypo = miu*(3*lambda+2*miu)/(lambda+miu);
			double NUhypo = lambda/(2*(lambda+miu));
			 //factor = Ehypo/(3*(1-2*NUhypo));
			factor= lambda + 2.0/3.0*miu;
			 accelerator = 3.0/(2.0*(1.0 + NUhypo));
			//KRATOS_WATCH(Ehypo);
			//KRATOS_WATCH(NUhypo);
			}

		
		
	}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		void localAssembleLHS( Matrix& K_patch, const Matrix K_element, const Vector codetable)
	{
		unsigned int local_size = K_element.size1();

			for (unsigned int i_local=0; i_local<local_size; i_local++)
			{
				unsigned int i_global=codetable[i_local];
				{
					for (unsigned int j_local=0; j_local<local_size; j_local++)
					{
						unsigned int j_global=codetable[j_local];
						K_patch(i_global,j_global) += K_element(i_local,j_local);
					}
				}
			}
	}
		
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@		
		bool CalculateVolumeChange(ModelPart& r_model_part)
	{
		double current_vol=0.0;
		double old_vol=0.0;
		//ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
		
		//KRATOS_WATCH(iteration);
		
		for(typename ModelPart::NodesContainerType::iterator ind = r_model_part.NodesBegin(); ind != r_model_part.NodesEnd();ind++)
            	   {
			//if(ind->FastGetSolutionStepValue(IS_INTERFACE)==1.0)
			   {
				
				old_vol += ind->FastGetSolutionStepValue(NODAL_H,2);
				current_vol += ind->FastGetSolutionStepValue(NODAL_H);
			   }
		   }
// 		KRATOS_WATCH(1.0*old_vol);
// 		KRATOS_WATCH(1.0*current_vol);
// 		KRATOS_WATCH(1.0*(current_vol-old_vol));
		if(pow(current_vol,2) < pow(old_vol,2))
		    {
			
			
			return true;
		    }
		else
			return false;

	}	
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void CalculateCorrectionForce(ModelPart::ElementsContainerType::iterator& ielem,ModelPart& r_model_part)
	{
		Geometry< Node<3> >& geom = ielem->GetGeometry();
		boost::numeric::ublas::bounded_matrix<double, TDim+1, TDim> DN_DX;	
		array_1d<double,TDim+1> N; 
		double volume;

		GeometryUtils::CalculateGeometryData(geom, DN_DX, N, volume);
		int matsize = TDim*(TDim+1);

		Vector Prforce=ZeroVector(matsize);
		Vector pressure=ZeroVector(TDim+1);
		//Get pressure
		for(int ii=0; ii<TDim+1; ++ii)
			pressure[ii]=geom[ii].FastGetSolutionStepValue(PRESSURE);

		Matrix	BPV = ZeroMatrix(matsize,TDim+1);
		CalculateBPV(BPV,N,DN_DX);
			
		noalias(Prforce)= volume*prod(BPV,pressure);
		 KRATOS_WATCH("**Inside CalculateCorrectionForce****");
		KRATOS_WATCH(Prforce);
		KRATOS_WATCH(volume);
		KRATOS_WATCH(pressure);

// 		for(int ii=0; ii<TDim+1; ++ii)
// 		     {
// 			Vector nodeforce = ZeroVector(TDim);
// 			int index=ii*TDim;
// 
// 			for(int jj=0; jj<TDim; ++jj)
// 				nodeforce[jj]+=Prforce[index+jj];
// 
// 			(geom[ii].FastGetSolutionStepValue(FORCE)) += nodeforce;
// 		     }
	}		
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	void CalculateBPV(
		Matrix& BPV,
		const array_1d<double,TDim+1> N,
		boost::numeric::ublas::bounded_matrix<double, TDim+1, TDim> DN_DX)

	{

		KRATOS_TRY

		
		for (unsigned int i=0;i<(TDim+1);i++)
		{
			unsigned int index = TDim*i;
			if (TDim == 2)
			{
		BPV(index+0,0)=N[0]*DN_DX(i,0);	BPV(index+0,1)=N[1]*DN_DX(i,0); BPV(index+0,2)=N[2]*DN_DX(i,0);
		BPV(index+1,0)=N[0]*DN_DX(i,1);	BPV(index+1,1)=N[1]*DN_DX(i,1); BPV(index+1,2)=N[2]*DN_DX(i,1);					
			}
			else
			{
		BPV(index+0,0)=N[0]*DN_DX(i,0);	BPV(index+0,1)=N[1]*DN_DX(i,0); BPV(index+0,2)=N[2]*DN_DX(i,0);BPV(index+0,3)=N[3]*DN_DX(i,0);

		BPV(index+1,0)=N[0]*DN_DX(i,1);	BPV(index+1,1)=N[1]*DN_DX(i,1); BPV(index+1,2)=N[2]*DN_DX(i,1);BPV(index+1,3)=N[3]*DN_DX(i,1);

			}

		}

		KRATOS_CATCH("")

	}		
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// ModelPart::NodesContainerType& rPatchNodes;
//  for(ModelPart::NodesContainerType::iterator inode=rBaseNodes.NodesBegin(); inode!= rBaseNodes.NodesEnd(); inode++)
// 	{
// 	double OldNodalVolume=inode->FastGetSolutionStepValue(NODAL_AREA,2);
// 	double CurrentNodalVolume = inode->FastGetSolutionStepValue(NODAL_AREA,1);
// 	if(Abs(Current
// 
// 	}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
				/*@} */
				/**@name Member Variables */
		/*@{ */

		/*@} */
		/**@name Private Operators*/
		/*@{ */


		/*@} */
		/**@name Private Operations*/
		/*@{ */


		/*@} */
		/**@name Private  Access */
		/*@{ */


		/*@} */     
		/**@name Private Inquiry */
		/*@{ */


		/*@} */   
		/**@name Un accessible methods */
		/*@{ */


		/*@} */   

	}; /* Class Scheme */

	/*@} */

	/**@name Type Definitions */       
	/*@{ */


	/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_INNER_VOLUMETRIC_STATIC_SCHEME  defined */

