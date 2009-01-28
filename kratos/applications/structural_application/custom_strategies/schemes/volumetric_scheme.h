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
*   Date:                $Date: 2009-01-15 18:47:00 $
*   Revision:            $Revision: 1.2 $
*
* ***********************************************************/


#if !defined(KRATOS_NEW_STANDARD_STATIC_SCHEME )
#define  KRATOS_NEW_STANDARD_STATIC_SCHEME


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
	class VolumetricScheme : public Scheme<TSparseSpace,TDenseSpace>
	{

	public:
		/**@name Type Definitions */       
		/*@{ */

		//typedef boost::shared_ptr< ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace> > Pointer;
		KRATOS_CLASS_POINTER_DEFINITION( VolumetricScheme);

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
		VolumetricScheme()
			: Scheme<TSparseSpace,TDenseSpace>()
		{}

		/** Destructor.
		*/
		virtual ~VolumetricScheme(){}


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
			   for(typename ModelPart::NodesContainerType::iterator ind = r_model_part.NodesBegin(); ind != r_model_part.NodesEnd();ind++)
				{
				ind->FastGetSolutionStepValue(IS_VISITED)=100.00;
				ind->FastGetSolutionStepValue(NODAL_AREA,2)=-1000000.00;
				
				}
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
			   KRATOS_WATCH(contor);
			   for(typename ModelPart::ElementsContainerType::iterator iel = r_model_part.ElementsBegin(); iel != r_model_part.ElementsEnd();iel++)
	        	      {
				
				//CalculateCorrectionForce(iel,r_model_part);
//
	        	      }	
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
				int dim=(iel->GetGeometry()).WorkingSpaceDimension();
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
	WeakPointerVector< Node<3> >& neighbor_nds = inode->GetValue(NEIGHBOUR_NODES);
	int nghsize = (neighbor_nds).size();
	ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

	//unsigned int TDim = inode->Dimension();
	int patchsize = (nghsize+1)*TDim;
	
	//KRATOS_WATCH(neighbor_nds);

	Vector Hoperator;
	DofsArrayType patchDofSet;	
	Vector patch_disp;
	Vector volumetricforce;
	Matrix patch_volumetric_K ;
	Matrix patch_mass;

	patchDofSet.reserve(patchsize);
	volumetricforce = ZeroVector(patchsize);
	Hoperator = ZeroVector(patchsize);
	patch_volumetric_K = ZeroMatrix(patchsize, patchsize) ;
	patch_mass = ZeroMatrix(patchsize, patchsize) ;
	patch_disp = ZeroVector(patchsize);	
	double deltaT= CurrentProcessInfo[DELTA_TIME];
	

	//patch`s local code table
	 patchDofSet.push_back(inode->pGetDof(ROTATION_X));
	 patchDofSet.push_back(inode->pGetDof(ROTATION_Y));
	  if(TDim==3)
	 patchDofSet.push_back(inode->pGetDof(ROTATION_Z));

	 for( WeakPointerVector< Node<3> >::iterator ind = neighbor_nds.begin(); ind!=neighbor_nds.end(); ind++)
	    {
		patchDofSet.push_back(ind->pGetDof(ROTATION_X));
	 	patchDofSet.push_back(ind->pGetDof(ROTATION_Y));
	  	 if(TDim==3)
	 	patchDofSet.push_back(ind->pGetDof(ROTATION_Z));
            }

	 int local_index = 0;
	for( typename DofsArrayType::iterator dof_iterator = patchDofSet.begin(); dof_iterator!= patchDofSet.end();
++dof_iterator)
		dof_iterator-> SetEquationId(local_index++);
		
//  	(inode-> FastGetSolutionStepValue(NODAL_H,2)) = (inode-> FastGetSolutionStepValue(NODAL_H));
//  	inode-> FastGetSolutionStepValue(NODAL_H)=0.0;
// 	(inode-> FastGetSolutionStepValue(NODAL_AREA,2)) = (inode-> FastGetSolutionStepValue(NODAL_AREA));
//  	inode-> FastGetSolutionStepValue(NODAL_AREA)=0.0;

	//fill patch velocity vector
	array_1d<double, TDim> nddisp = inode->FastGetSolutionStepValue(VELOCITY);
	for(int ii=0; ii<TDim; ii++)
	 patch_disp[ii] = nddisp[ii];
	
	int counter = 1;	
	for(WeakPointerVector< Node<3> >::iterator ii=neighbor_nds.begin(); ii!=neighbor_nds.end(); ii++)
	    {
	      nddisp = ii->FastGetSolutionStepValue(VELOCITY);
		for(int j=0; j<TDim; j++)
		 patch_disp[counter*TDim + j] = nddisp[j];
		 counter += 1;
	    }		

	//calculate H
	for( WeakPointerVector< Element >::iterator iel = neighbor_els.begin(); iel != neighbor_els.end(); iel++)
	     {
		Geometry< Node<3> >& geom = iel->GetGeometry();
		boost::numeric::ublas::bounded_matrix<double, TDim+1, TDim> DN_DX;	
		array_1d<double,TDim+1> N; 
		double volume;

		GeometryUtils::CalculateGeometryData(geom, DN_DX, N, volume);

		//building H for element
		array_1d< double, TDim*(TDim+1)> local_H;
		array_1d< unsigned int, TDim*(TDim+1)> ielindex;
	
		   for( int ii=0; ii<TDim+1;ii++)
		     {
			int index = ii*TDim; 
			ielindex[index+0]=geom[ii].GetDof(ROTATION_X).EquationId();
			ielindex[index+1]=geom[ii].GetDof(ROTATION_Y).EquationId();
			 if ( TDim == 3)
			ielindex[ii*TDim+2]=geom[ii].GetDof(ROTATION_Z).EquationId();
		     for(int jj=0; jj<TDim; jj++)
			  local_H[index+jj]= DN_DX(ii,jj)*volume;
		      }

		
		//assembling
		localAssembleRHS(Hoperator, local_H, ielindex);
	    }
	//calculate volume
	double node_vol = inner_prod(Hoperator,patch_disp);
	(inode-> FastGetSolutionStepValue(NODAL_H,1)) = (inode-> FastGetSolutionStepValue(NODAL_H));
	inode-> FastGetSolutionStepValue(NODAL_H)=node_vol;
	inode-> FastGetSolutionStepValue(NODAL_AREA)=node_vol;

	//yes or no
// 	double OldNodalVolume=inode->FastGetSolutionStepValue(NODAL_AREA,2);
// 	double CurrentNodalVolume = inode->FastGetSolutionStepValue(NODAL_AREA,1);
	//double OldNodalVolume=inode->FastGetSolutionStepValue(NODAL_H,1);
	//double CurrentNodalVolume = inode->FastGetSolutionStepValue(NODAL_H);
// 	KRATOS_WATCH("***********THIS IS PATCH VOLUME*****************");
// 	KRATOS_WATCH(fabs(OldNodalVolume));
// 	KRATOS_WATCH(fabs(CurrentNodalVolume));
// 	KRATOS_WATCH("***********THIS IS PATCH VOLUME*****************");
	//if(OldNodalVolume == 0.0)
	//	OldNodalVolume = CurrentNodalVolume;

// 	if (fabs(OldNodalVolume)<fabs(CurrentNodalVolume))
// 		inode->FastGetSolutionStepValue(IS_VISITED)=0.0;

	//main loop
 	if(inode->FastGetSolutionStepValue(IS_INTERFACE)==10.0 && inode->FastGetSolutionStepValue(IS_VISITED)==100.0)
	{
		//patch number contor
		contor++;
	//calculate patch volumetric stiffness
	  for( WeakPointerVector< Element >::iterator iel = neighbor_els.begin(); iel != neighbor_els.end(); iel++)
	     {
		Geometry< Node<3> >& geom = iel->GetGeometry();
		boost::numeric::ublas::bounded_matrix<double, TDim+1, TDim> DN_DX;	
		array_1d<double,TDim+1> N; 
		double volume;
		
		
		GeometryUtils::CalculateGeometryData(geom, DN_DX, N, volume);

		//building H for element
		array_1d< double, TDim*(TDim+1)> local_H;
		array_1d< int, TDim*(TDim+1)> ielindex;
	
		   for( int ii=0; ii<TDim+1;ii++)
		     {
			int index = ii*TDim; 
			ielindex[index+0]=geom[ii].GetDof(ROTATION_X).EquationId();
			ielindex[index+1]=geom[ii].GetDof(ROTATION_Y).EquationId();
			 if ( TDim == 3)
			ielindex[ii*TDim+2]=geom[ii].GetDof(ROTATION_Z).EquationId();
		     for(int jj=0; jj<TDim; jj++)
			  local_H[index+jj]= DN_DX(ii,jj)*volume;
		      }
		//element volume
		array_1d< double, TDim*(TDim+1)> eledisp;
		for (unsigned int i=0;i<(TDim+1);i++)
		{
			unsigned int index = i*TDim;
			eledisp[index] = (iel->GetGeometry())[i].GetSolutionStepValue(DISPLACEMENT_X);
			eledisp[index + 1] = (iel->GetGeometry())[i].GetSolutionStepValue(DISPLACEMENT_Y);
			if(TDim == 3)
				eledisp[index + 2] = (iel->GetGeometry())[i].GetSolutionStepValue(DISPLACEMENT_Z);
	
		}
		//double ele_vol = inner_prod(local_H,eledisp);
		//inode-> FastGetSolutionStepValue(NODAL_H)+=(ele_vol);
		
		//inode-> FastGetSolutionStepValue(NODAL_H)=10000;
		//KRATOS_WATCH(ele_vol);
		//KRATOS_WATCH(inode-> FastGetSolutionStepValue(NODAL_H));
		//std::cout<<inode-> FastGetSolutionStepValue(NODAL_H)<<std::endl;
		
		Matrix FF=ZeroMatrix(2,2);
		Element::GeometryType::JacobiansType nowJ;
		geom.Jacobian(nowJ);
		std::vector< Matrix > refrenceJ;
		iel->CalculateOnIntegrationPoints(AUXILIARY_MATRIX_1, refrenceJ, CurrentProcessInfo);
		noalias(FF) = prod(nowJ[0],refrenceJ[0]);
		double detF = MathUtils<double>::Det(FF);
		//Matrix FF=ZeroMatrix(2,2);
		//KRATOS_WATCH(nowJ[0]);
		//KRATOS_WATCH(refrenceJ[0]);
		//noalias(FF) = prod(nowJ[0],mInvJ0[0]);
		//bulding volumetric stiffness for each element
		std::vector< double> B;
		B.reserve(patchsize);
//		boost::numeric::ublas::bounded_matrix<double, TDim*(TDim+1),TDim*(TDim+1)> element_volumetric_K;
		Matrix element_volumetric_K;
		
		element_volumetric_K.resize(patchsize, patchsize);
// 		double miu = iel->GetProperties()[MIU];
// 		double lambda = iel->GetProperties()[LAMBDA];
// 		miu = miu*deltaT/detF;
// 		lambda= lambda*deltaT/detF;
// 		double Ehypo = miu*(3*lambda+2*miu)/(lambda+miu);
// 		double NUhypo = lambda/(2*(lambda+miu));
// 		double factor = Ehypo/(3*(1-2*NUhypo));
// 		double accelerator = 3.0/(2.0*(1.0 + NUhypo));
		double factor=1.0/3.0;
		double accelerator=1.0;

		
		factoraccelerator(iel,accelerator,factor,deltaT,detF);
	
		//KRATOS_WATCH(accelerator);
		//KRATOS_WATCH(factor);
		element_volumetric_K = accelerator*factor/volume*outer_prod(local_H, local_H);
		//iel->DampMatrix(element_volumetric_K,  CurrentProcessInfo);
		//element_volumetric_K*=accelerator;
		//CALCULATE MASS
		Matrix element_mass;
		iel->MassMatrix(element_mass,CurrentProcessInfo);
		
		//assembling
// 		localAssembleRHS(Hoperator, local_H, ielindex);
		//KRATOS_WATCH("H is assembeled");
		localAssembleLHS(patch_volumetric_K, element_volumetric_K, ielindex);
		localAssembleLHS(patch_mass, element_mass, ielindex);

	    }

	//KRATOS_WATCH(Hoperator);
	//KRATOS_WATCH(patch_volumetric_K);
	//****finalize patch volumetric stiffness calculation
	
		Matrix AV;
		Matrix vol_dynLHS;
		vol_dynLHS=ZeroMatrix(patchsize, patchsize);


		AV = outer_prod(Hoperator,Hoperator);
	//KRATOS_WATCH(AV);
		double denom = inner_prod(Hoperator,Hoperator);
		AV/=denom;
	
		vol_dynLHS = 1/deltaT*patch_mass + patch_volumetric_K;
	
	//patch_volumetric_K = 1/denom*prod(patch_volumetric_K,AV);


	//calculate volumetric force
		double overlap_factor = 1.0/3.0;
	//KRATOS_WATCH(overlap_factor);
		Vector vol_proj;
		vol_proj = ZeroVector(patchsize);	
		vol_proj = prod( AV, patch_disp);
		double oldvol=inode->FastGetSolutionStepValue(NODAL_AREA,1);
		vol_proj += (Hoperator*oldvol)/denom;
		noalias(volumetricforce) = -1*overlap_factor*prod( vol_dynLHS, vol_proj);
	//KRATOS_WATCH(patch_disp);

	//KRATOS_WATCH(volumetricforce);

	// calculate pressure how to save? how reset?
		double pressure=0.0;
		//double nominator=0.0;
		//double denominator = 0.0;
		for(int ii=0; ii< patchsize; ii++)
		{
		//KRATOS_WATCH(Hoperator[ii]);
	 	   if ( Hoperator[ii]*Hoperator[ii]> 10e-14)
			{
			//KRATOS_WATCH( Hoperator[ii]*Hoperator[ii]);
			pressure += (1/overlap_factor)*volumetricforce[ii]/Hoperator[ii]; 
			//nominator += 3.0*volumetricforce[ii]*Hoperator[ii];
			//denominator += Hoperator[ii]*Hoperator[ii];
			}
		}	
		
	
		pressure /= Hoperator.size();
	//pressure = nominator/denominator;
	//KRATOS_WATCH(pressure);
	//KRATOS_WATCH(inode-> FastGetSolutionStepValue(PRESSURE));
		inode-> FastGetSolutionStepValue(PRESSURE) += pressure;
	//KRATOS_WATCH(inode-> FastGetSolutionStepValue(PRESSURE));
// 	double node_vol = inner_prod(Hoperator,patch_disp);
// 
// 	inode-> FastGetSolutionStepValue(NODAL_AREA)=node_vol;
	//KRATOS_WATCH(inode-> FastGetSolutionStepValue(NODAL_H));
	//KRATOS_WATCH(inode-> FastGetSolutionStepValue(PRESSURE));

	//assemble patch force to global force
	//std::vector<unsigned int> GlobalPatchIndex;
	//GlobalPatchIndex.resize(patchsize);
	//array_1d<unsigned int, (nghsize+1)*TDim> GlobalPatchIndex;
	//GlobalPatchIndex[0] = inode-> GetDof(DISPLACEMENT_X).EquationId();
	//GlobalPatchIndex[1] = inode-> GetDof(DISPLACEMENT_Y).EquationId();
	  //if (TDim == 3)
	//GlobalPatchIndex[2] = inode-> GetDof(DISPLACEMENT_Z).EquationId();
	//unsigned int index=TDim;
	//for( WeakPointerVector< Node<3> >::iterator ind = neighbor_nds.begin(); ind!=neighbor_nds.end(); ind++)
	 //   {
	//		GlobalPatchIndex[index] = ind-> GetDof(DISPLACEMENT_X).EquationId();
		//KRATOS_WATCH(ind-> GetDof(DISPLACEMENT_X).EquationId());
	//	GlobalPatchIndex[index+1] = ind-> GetDof(DISPLACEMENT_Y).EquationId();
	//		  if (TDim == 3)
	//		GlobalPatchIndex[++index+2] = ind-> GetDof(DISPLACEMENT_Z).EquationId();
	//		index += TDim;
	  //  }
	//KRATOS_WATCH(mDeltaF);
	//KRATOS_WATCH(volumetricforce);

	//**update total volumetric force***
		 for(int ii = 0; ii<TDim; ii++)
			inode->FastGetSolutionStepValue(FORCE)[ii] += volumetricforce[ii];

		int index = TDim;
		 for( WeakPointerVector< Node<3> >::iterator ind = neighbor_nds.begin(); ind!=neighbor_nds.end(); ind++)
	 	    {
			 for(int ii = 0; ii<TDim; ii++)
			   ind->FastGetSolutionStepValue(FORCE)[ii] += volumetricforce[index+ii];
			 index += TDim;
	  	   }
	//************************************






 	}

	//KRATOS_WATCH(mDeltaF);
	}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	void localAssembleRHS(Vector& Hpatch, const Vector localH, const vector< int > codetable)
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
	//		double Ehypo = miu*(3*lambda+2*miu)/(lambda+miu);
			double NUhypo = lambda/(2*(lambda+miu));
			 //factor = Ehypo/(3*(1-2*NUhypo));
			factor= lambda + 2.0/3.0*miu;
			 accelerator = 3.0/(2.0*(1.0 + NUhypo));
			//KRATOS_WATCH(Ehypo);
			//KRATOS_WATCH(NUhypo);
			}

		
		
	}

//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		void localAssembleLHS( Matrix& K_patch, const Matrix K_element, const boost::numeric::ublas::vector< int> codetable)
	{
		unsigned int local_size = K_element.size1();

			for (unsigned int i_local=0; i_local<local_size; i_local++)
			{
				 int i_global=codetable[i_local];
				{
					for (unsigned int j_local=0; j_local<local_size; j_local++)
					{
						 int j_global=codetable[j_local];
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

#endif /* KRATOS_NEW_STANDARD_STATIC_SCHEME  defined */

