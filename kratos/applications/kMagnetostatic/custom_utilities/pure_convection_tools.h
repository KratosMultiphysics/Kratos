/*
==============================================================================
KratosR1ConvectionDiffusionApplication 
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
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-11-14 15:25:09 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_PURE_CONVECTION_UTILITIES_INCLUDED )
#define  KRATOS_LARANGIAN_FLUID_UTILITIES_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// External includes 


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "utilities/geometry_utilities.h"
#include "kElectrostatic.h"

namespace Kratos
{

    template<int TDim,class TSparseSpace,class TLinearSolver>
	class PureConvectionUtilities
	{
		public:

		typedef typename TSparseSpace::MatrixType TSystemMatrixType;
		typedef typename TSparseSpace::VectorType TSystemVectorType;
        typedef PointerVectorSet< Dof<double> , IndexedObject> DofsArrayType;


		PureConvectionUtilities(){};
		~PureConvectionUtilities(){};

 

        //************************************************************************
        //************************************************************************
	void ConstructSystem(ModelPart& model_part, Variable<double>& rScalarVar, Variable< array_1d<double,3> >& rTransportVel, Variable< array_1d<double,3> >& rMeshVel)
        {
            KRATOS_TRY

            //loop over nodes with neighbours
            mDofSet = DofsArrayType();
	    //mDofSet.resize(0);
            mDofSet.reserve(model_part.Nodes().size() );

            int tot_nnz=0;
            int tot_row=0;
            for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
            {
                if( (it->GetValue(NEIGHBOUR_NODES)).size() != 0 )
                {
                    mDofSet.push_back( it->pGetDof(rScalarVar) );
                    tot_row += 1;
                    tot_nnz +=  (it->GetValue(NEIGHBOUR_NODES)).size()+1;
                }
             }
     
             //fills the DofList and give a unique progressive tag to each node 
             
			int free_id = 0;
			int fix_id = mDofSet.size();

			for (typename DofsArrayType::iterator dof_iterator = mDofSet.begin(); dof_iterator != mDofSet.end(); ++dof_iterator)
				if (dof_iterator->IsFixed())
					dof_iterator->SetEquationId(--fix_id);
				else
					dof_iterator->SetEquationId(free_id++);

			mEquationSystemSize = fix_id;
                
            std::vector< int > work_array; 
            work_array.reserve(1000);
    
            //guessing the total size of the matrix 
            mA.resize(mEquationSystemSize, mEquationSystemSize, tot_nnz); 
    
            //getting the dof position 
            //unsigned int dof_position = (model_part.NodesBegin())->GetDofPosition(rScalarVar); 

            //building up the matrix graph row by row 
            int total_size = 0; 
            for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
            {
		unsigned int index_i = it->GetDof(rScalarVar).EquationId();  

		if(index_i < mEquationSystemSize && (it->GetValue(NEIGHBOUR_NODES)).size() != 0)
                {
                    WeakPointerVector< Node<3> >& neighb_nodes = it->GetValue(NEIGHBOUR_NODES); 
        
                    //filling the first neighbours list 
                    work_array.push_back(index_i); 
                    for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin(); i != neighb_nodes.end(); i++) 
                    { 
                        int index_j = i->GetDof(rScalarVar).EquationId(); 
                        work_array.push_back(index_j); 
                    } 
                                    
                    //sorting the indices and elminating the duplicates 
                    std::sort(work_array.begin(),work_array.end()); 
                    typename std::vector<int>::iterator new_end = std::unique(work_array.begin(),work_array.end()); 
                    unsigned int number_of_entries = new_end - work_array.begin(); 
                                    
                    //filling up the matrix 
                    for(unsigned int j=0; j<number_of_entries; j++) 
                    { 
                        mA.push_back(index_i,work_array[j] , 0.00); 
                    } 
    
                    //clearing the array for the next step 
                    work_array.erase(work_array.begin(),work_array.end()); 
                    total_size += number_of_entries;   
                }
                            


             
            }
            
            mDx.resize(mA.size1(),false);
            mb.resize(mA.size1(),false);
            KRATOS_CATCH("")
        }
    
        //************************************************************************
        //************************************************************************
	void ConvectScalarVar(ModelPart& model_part, 
			        typename TLinearSolver::Pointer linear_solver, 
	 			Variable<double>& rScalarVar, 
     				const Variable< array_1d<double,3> >& rTransportVel, 
     				const Variable< array_1d<double,3> >& rMeshVel, 
	 			const Variable< double >& rProjVar, 
	 			unsigned int time_order )
        {
            KRATOS_TRY
    
                TSparseSpace::SetToZero(mA);
                TSparseSpace::SetToZero(mb);
                TSparseSpace::SetToZero(mDx);

	        ProcessInfo& CurrentProcessInfo = model_part.GetProcessInfo();
                double dt = CurrentProcessInfo[DELTA_TIME];
		
		Vector BDFcoeffs(time_order+1);
		
		if(time_order == 2)
		{
			if(model_part.GetBufferSize() < 3)
				KRATOS_ERROR(std::logic_error,"insufficient buffer size for BDF2","")

			BDFcoeffs[0] =	1.5 / dt;	//coefficient for step n+1
			BDFcoeffs[1] =	-2.0 / dt;//coefficient for step n
			BDFcoeffs[2] =	0.5 / dt;//coefficient for step n-1
		}
		else
		{
			BDFcoeffs[0] =	1.0 / dt;	//coefficient for step n+1
			BDFcoeffs[1] =	-1.0 / dt;//coefficient for step n
		}
    
                //**********************************
                //BUILD PHASE
                //**********************************
		boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim> DN_DX; 
		boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim+1> lhs_contribution; 
		array_1d<double,TDim+1> N; 
		array_1d<unsigned int ,TDim+1> local_indices; 
		array_1d<double,TDim+1> rhs_contribution; 
		array_1d<double,TDim+1> dp_vector; 
		array_1d<double,TDim+1> MassFactors; 
		array_1d<double,TDim> vel_gauss;
 		array_1d<double,TDim> proj;
		array_1d<double,TDim+1> u_DN;
		array_1d<double,TDim+1> temp_vec_np;
		

		
		for(unsigned int i=0; i<TDim+1; i++)
			MassFactors[i] = 1.0/double(TDim+1);

		//getting the dof position 
//		unsigned int dof_position = (model_part.NodesBegin())->GetDofPosition(rScalarVar); 

		double lumping_factor = 1.0/double(TDim+1.0); 
		const unsigned int number_of_points = TDim+1;
		
		array_1d<double,3> vg; 
		for(ModelPart::ElementsContainerType::iterator i = model_part.ElementsBegin();  
			i!=model_part.ElementsEnd(); i++) 
		{	 
			Geometry< Node<3> >& geom = i->GetGeometry();

			//calculating elemental values 
			double Volume; 
			GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume); 
			
			
			//finiding local indices 
			for(int ii = 0; ii<TDim+1; ii++) 
			{ 
				local_indices[ii] = geom[ii].GetDof(rScalarVar).EquationId(); 

/*				dp_vector[ii] = geom[ii].FastGetSolutionStepValue(rScalarVar) 
						- geom[ii].FastGetSolutionStepValue(rScalarVar,1);*/
			} 

			//calculating convective velocity
			noalias(vel_gauss) = ZeroVector(TDim);
			double proj = 0.0;
			for(unsigned int i = 0; i<number_of_points; i++)
			{
				const array_1d<double,3>& v = geom[i].FastGetSolutionStepValue(rTransportVel);
				const array_1d<double,3>& w = geom[i].FastGetSolutionStepValue(rMeshVel);
				proj += geom[i].FastGetSolutionStepValue(rProjVar);
				for(unsigned int j = 0; j<TDim; j++)
				{
					vel_gauss[j] += v[j] - w[j];
				}
			}
			vel_gauss *= lumping_factor;
			proj *= lumping_factor;

			
			//CONVECTIVE CONTRIBUTION TO THE STIFFNESS MATRIX
			noalias(u_DN) = prod(DN_DX , vel_gauss);
			noalias(lhs_contribution) = outer_prod(N,u_DN);

			//CONVECTION STABILIZING CONTRIBUTION (Suu)
			double tau = CalculateTAU(vel_gauss,Volume);
			noalias(lhs_contribution) += (tau) * outer_prod(u_DN,u_DN);

			//INERTIA CONTRIBUTION
			for(unsigned int iii = 0; iii<number_of_points; iii++)
				lhs_contribution(iii,iii) += BDFcoeffs[0] * MassFactors[iii];

			//adding projection contribution
			//RHS += Suy * proj[component] 
 			noalias(rhs_contribution) = (tau*proj)*u_DN;  

			
			//adding the inertia terms
			// RHS += M*vhistory 
			//calculating the historical velocity
			for(unsigned int iii = 0; iii<number_of_points; iii++)
			{
				temp_vec_np[iii] =  BDFcoeffs[1]*geom[iii].FastGetSolutionStepValue(rScalarVar,1);
			}
			for(unsigned int step = 2; step<BDFcoeffs.size(); step++)
			{
				for(unsigned int iii = 0; iii<number_of_points; iii++)
					temp_vec_np[iii] += BDFcoeffs[step]*geom[iii].FastGetSolutionStepValue(rScalarVar,step);
			}
			for(unsigned int iii = 0; iii<number_of_points; iii++)	
				rhs_contribution[iii] -= MassFactors[iii]*temp_vec_np[iii];
// 				rhs_contribution[iii] = - MassFactors[iii]*temp_vec_np[iii];				


			//subtracting the dirichlet term
			// RHS -= LHS*rScalarVars
			for(unsigned int iii = 0; iii<number_of_points; iii++)
				temp_vec_np[iii] = geom[iii].FastGetSolutionStepValue(rScalarVar);
			noalias(rhs_contribution) -= prod(lhs_contribution,temp_vec_np);
		
		//multiplying by Volume, rho and density
			rhs_contribution *= Volume;
			lhs_contribution *= Volume;

			AssembleLHS(mA,lhs_contribution,local_indices);
			AssembleRHS(mb,rhs_contribution,local_indices);



		}

                //**********************************
                //SOLVE PHASE
                //**********************************
    //             KRATOS_WATCH(mA);
  //               KRATOS_WATCH(mb);
    //             KRATOS_WATCH(mDx);
//		std::cout << "convection residual norm" << TSparseSpace::TwoNorm(mb) << std::endl;
		
                linear_solver->Solve(mA,mDx,mb);
               std::cout << *(linear_solver) << std::endl;
      //           KRATOS_WATCH(mDx);
    
                //**********************************
                //UPDATE rScalarVarS
                //**********************************
                for(typename DofsArrayType::iterator i_dof = mDofSet.begin() ; i_dof != mDofSet.end() ; ++i_dof)
                {
                    if(i_dof->IsFree())
                    {
                        i_dof->GetSolutionStepValue() += mDx[i_dof->EquationId()];
                    }
                }
   
            KRATOS_CATCH("")
        }
	
	
	
	

        //************************************************************************
        //************************************************************************
	void CalculateProjection(ModelPart& model_part, 
				const Variable<double>& rScalarVar, 
    				const Variable<double>& rNodalArea,
				const Variable< array_1d<double,3> >& rTransportVel, 
				const Variable< array_1d<double,3> >& rMeshVel, 
				Variable< double >& rProjVar )
	{
		KRATOS_TRY
    
   
                //**********************************
                //BUILD PHASE
                //**********************************
		boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim> DN_DX; 
		array_1d<double,TDim+1> N; 
		array_1d<double,TDim> vel_gauss;
		array_1d<double,TDim+1> temp_vec_np;
		array_1d<double,TDim+1> u_DN;
		
		double lumping_factor = 1.0/double(TDim+1.0); 
		const unsigned int number_of_points = TDim+1;
		
		//set to zero nodal areas
		for(ModelPart::NodesContainerType::iterator i = model_part.NodesBegin();  
				  i!=model_part.NodesEnd(); i++) 
		{
			i->FastGetSolutionStepValue(rNodalArea) = 0.0;
			i->FastGetSolutionStepValue(rProjVar) = 0.0;
		}
		
		//build up nodal areas and projections
		array_1d<double,3> vg; 
		for(ModelPart::ElementsContainerType::iterator i = model_part.ElementsBegin();  
				  i!=model_part.ElementsEnd(); i++) 
		{	 

			Geometry< Node<3> >& geom = i->GetGeometry();

			//calculating elemental values 
			double Volume; 
			GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume); 
			
			//calculating convective velocity
			double nodal_area = Volume*lumping_factor;
			noalias(vel_gauss) = ZeroVector(TDim);
			for(unsigned int i = 0; i<number_of_points; i++)
			{
				const array_1d<double,3>& v = geom[i].FastGetSolutionStepValue(rTransportVel);
				const array_1d<double,3>& w = geom[i].FastGetSolutionStepValue(rMeshVel);
				
				//adding up nodal area as needed
				geom[i].FastGetSolutionStepValue(rNodalArea) += nodal_area;
				
				for(unsigned int j = 0; j<TDim; j++)
				{
					vel_gauss[j] += v[j] - w[j];
				}
			}
			vel_gauss *= lumping_factor;

			//filling up a vector with all of the nodal values			
			for(unsigned int iii = 0; iii<number_of_points; iii++)
				temp_vec_np[iii] = geom[iii].FastGetSolutionStepValue(rScalarVar);
			
			//
			noalias(u_DN) = prod(DN_DX , vel_gauss);
			double conv_proj = inner_prod( u_DN , temp_vec_np);
			conv_proj *= Volume * lumping_factor;

			for(unsigned int iii = 0; iii<number_of_points; iii++)
				geom[iii].FastGetSolutionStepValue(rProjVar) += conv_proj;
		}
		
		//dividing by the overall mass matrix
		for(ModelPart::NodesContainerType::iterator i = model_part.NodesBegin();  
				  i!=model_part.NodesEnd(); i++) 
		{
			i->FastGetSolutionStepValue(rProjVar) /= i->FastGetSolutionStepValue(rNodalArea);
		}		

		KRATOS_CATCH("")
	}
    
        void ClearSystem()
        {
            KRATOS_TRY
	      mDofSet = DofsArrayType();
	    //		mDofSet.resize(0);
            mDx.resize(0,false);
            mb.resize(0,false);
            mA.resize(0,0,false);
    
            KRATOS_CATCH("")
        }




    private:

    unsigned int mEquationSystemSize;		
	TSystemVectorType mDx;
	TSystemVectorType mb;
	TSystemMatrixType mA;
    DofsArrayType mDofSet;





		//**************************************************************************
		void AssembleLHS(
			TSystemMatrixType& A,
			const boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim+1>& LHS_Contribution,
			const array_1d<unsigned int,TDim+1>& EquationId
			)
		{
			unsigned int local_size = LHS_Contribution.size1();

			for (unsigned int i_local=0; i_local<local_size; i_local++)
			{
				unsigned int i_global=EquationId[i_local];
				if ( i_global < mEquationSystemSize )
				{
					for (unsigned int j_local=0; j_local<local_size; j_local++)
					{
						unsigned int j_global=EquationId[j_local];
						if ( j_global < mEquationSystemSize )
							A(i_global,j_global) += LHS_Contribution(i_local,j_local);
					}
				}
			}
		}



		//**************************************************************************
		void AssembleRHS(
			TSystemVectorType& b,
			const array_1d<double,TDim+1>& RHS_Contribution,
			const array_1d<unsigned int,TDim+1>& EquationId
			)
		{
			unsigned int local_size = RHS_Contribution.size();

			for (unsigned int i_local=0; i_local<local_size; i_local++)
				{
					unsigned int i_global=EquationId[i_local];
					if ( i_global < mEquationSystemSize ) //on "free" DOFs
					{	// ASSEMBLING THE SYSTEM VECTOR
						b[i_global] += RHS_Contribution[i_local];
					}
				}

		}
		
		double CalculateTAU(const array_1d<double,2>& vel, const double& Volume)
		{
			//calculating parameter tau 
//			double c1 = 4.00;
			double c2 = 2.00;
			double h = sqrt(2.00*Volume);
			double norm_u =norm_2(vel);
			if(norm_u > 1e-15)
				return h / (  c2*norm_u );
			return 0.0;
		}
		
		double CalculateTAU(const array_1d<double,3>& vel, const double& Volume)
		{
			//calculating parameter tau 
//			double c1 = 4.00;
			double c2 = 2.00;
			double h =  pow(6.00*Volume,0.3333333);
			double norm_u =norm_2(vel);
			if(norm_u > 1e-15)
				return h / ( c2*norm_u );
			return 0.0;
		}

    };

}  // namespace Kratos.

#endif // KRATOS_LARANGIAN_FLUID_UTILITIES_INCLUDED  defined 


