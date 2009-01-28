/*
==============================================================================
KratosPFEMApplication 
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
//   Date:                $Date: 2007-10-24 08:20:31 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_LARANGIAN_FLUID_UTILITIES_INCLUDED )
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
#include "PFEM_application.h"

namespace Kratos
{

    template<int TDim,class TSparseSpace,class TLinearSolver>
	class LagrangianEulerSolver
	{
		public:

		typedef typename TSparseSpace::MatrixType TSystemMatrixType;
		typedef typename TSparseSpace::VectorType TSystemVectorType;
        typedef PointerVectorSet< Dof<double> , IndexedObject> DofsArrayType;


		LagrangianEulerSolver(){};
		~LagrangianEulerSolver(){};

        void Predict(ModelPart& model_part)
        {
            KRATOS_TRY
            KRATOS_WATCH("quii");

            //newmark coefficients
    	     ProcessInfo& CurrentProcessInfo = model_part.GetProcessInfo();
             double DeltaTime = CurrentProcessInfo[DELTA_TIME];


            double AlphaBossak = 0.0;
            double BetaNewmark = 0.25*(1.00-AlphaBossak)*(1.00-AlphaBossak);
            double GammaNewmark = 0.5-AlphaBossak;

			double a0 = 1.0/(BetaNewmark*pow(DeltaTime,2));
			double a1 = GammaNewmark / (BetaNewmark*DeltaTime);
			double a2 = 1.0/(BetaNewmark*DeltaTime);
			double a3 = 1.0/(2.0*BetaNewmark) - 1.0;
			double a4 = GammaNewmark/BetaNewmark - 1.0;	
			double a5 = DeltaTime*0.5*(GammaNewmark/BetaNewmark-2.0);
//			double am = (1.0-AlphaBossak)/(BetaNewmark*pow(DeltaTime,2));


            //predicting end of step position
            array_1d<double,3> displacement_step_increment; //this is the displacement increment from step n to n+1
			for(ModelPart::NodesContainerType::iterator in = model_part.NodesBegin(); in!=model_part.NodesEnd(); in++)
			{   

                if(in->FastGetSolutionStepValue(IS_STRUCTURE) != 1)
                {
                    array_1d<double,3>& disp = in->FastGetSolutionStepValue(DISPLACEMENT);
                    array_1d<double,3>& vel = in->FastGetSolutionStepValue(VELOCITY);
                    array_1d<double,3>& acc = in->FastGetSolutionStepValue(ACCELERATION);

                    const array_1d<double,3>& disp_old = in->FastGetSolutionStepValue(DISPLACEMENT,1);
                    const array_1d<double,3>& vel_old = in->FastGetSolutionStepValue(VELOCITY,1);
                    const array_1d<double,3>& acc_old = in->FastGetSolutionStepValue(ACCELERATION,1);

                    for(unsigned int j=0;j<TDim;j++) 
                    {
                        //predicting the displacement increment in the step
                        displacement_step_increment[j] = vel_old[j]*DeltaTime; 

                        //updating displacements
                        disp[j] = disp_old[j] + displacement_step_increment[j];

                        //updating velocity
                        vel[j] = a1*displacement_step_increment[j] - a4*vel_old[j] - a5*acc_old[j];

                        //updating acceleration
                        acc[j] = a0*displacement_step_increment[j] - a2*vel_old[j] - a3*acc_old[j];
                    }

                    //updating coordinates
                    in->X() = in->X0() + disp[0];
                    in->Y() = in->Y0() + disp[1];
                    in->Z() = in->Z0() + disp[2];
                }

            }
            KRATOS_CATCH("");
        
	}	


        void SolveConvectionStep(ModelPart& model_part)
        {
            KRATOS_TRY
            KRATOS_WATCH("quii");

            //newmark coefficients
	 ProcessInfo& CurrentProcessInfo = model_part.GetProcessInfo();
         double DeltaTime = CurrentProcessInfo[DELTA_TIME];


            double AlphaBossak = 0.0;
            double BetaNewmark = 0.25*(1.00-AlphaBossak)*(1.00-AlphaBossak);
            double GammaNewmark = 0.5-AlphaBossak;

			double a0 = 1.0/(BetaNewmark*pow(DeltaTime,2));
			double a1 = GammaNewmark / (BetaNewmark*DeltaTime);
			double a2 = 1.0/(BetaNewmark*DeltaTime);
			double a3 = 1.0/(2.0*BetaNewmark) - 1.0;
			double a4 = GammaNewmark/BetaNewmark - 1.0;	
			double a5 = DeltaTime*0.5*(GammaNewmark/BetaNewmark-2.0);
			double am = (1.0-AlphaBossak)/(BetaNewmark*pow(DeltaTime,2));

            //inverse coefficient for the acceleration term
            double inverse_inertia_coeff = 1.0/am;

            //reset nodal area and press proj
        	for(ModelPart::NodesContainerType::iterator in = model_part.NodesBegin(); in!=model_part.NodesEnd(); in++)
        	{   
                in->FastGetSolutionStepValue(NODAL_AREA) = 0.0;

                array_1d<double,3>& press_proj = in->FastGetSolutionStepValue(PRESS_PROJ);
                for(unsigned int i=0; i<TDim; i++) press_proj[i]=0.0;
            }

            //elemental calculation of pressure proj (Minverse*G*p)
			//compute the projection (first step
            double nn_inv = 1.0/(double(TDim+1));
            array_1d<double,TDim> aux;
            boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim> DN_DX; 
            array_1d<double,TDim+1> N; 

			for(ModelPart::ElementsContainerType::iterator ie = model_part.ElementsBegin(); 
				ie!=model_part.ElementsEnd(); ie++)
			{	
				//calculating shape functions values
				Geometry< Node<3> >& geom = ie->GetGeometry();
				double Volume;
				
				GeometryUtils::CalculateGeometryData(geom,DN_DX,N,Volume);

                //calculate pressure at the gauss point
                double p_avg = geom[0].FastGetSolutionStepValue(PRESSURE);
                for(unsigned int i=1;i<TDim+1;i++)
                    p_avg += geom[i].FastGetSolutionStepValue(PRESSURE);
                p_avg *= nn_inv;

                //pre multiplying by volume to minimize calculations
                p_avg *= Volume; 
                double nodal_area = Volume * nn_inv;

                //completing the calculation p_proj = A * grad(p)  
				for(int I = 0; I<TDim+1; I++)
				{
                    //adding the pressure gradient to the nodes
                    array_1d<double,3>& press_proj = geom[I].FastGetSolutionStepValue(PRESS_PROJ);
					for(unsigned int j=0;j<TDim;j++)
                        press_proj[j] += DN_DX(I,j)*p_avg; 

                    //computation of the nodal area
                    geom[I].FastGetSolutionStepValue(NODAL_AREA) += nodal_area;
				}
			}
 
            //finalizing calculation and calculating displacement increment
            array_1d<double,3> delta_disp; //this is the displacement correction
            array_1d<double,3> displacement_step_increment; //this is the displacement increment from step n to n+1
            array_1d<double,3> zero = ZeroVector(3);
			for(ModelPart::NodesContainerType::iterator in = model_part.NodesBegin(); in!=model_part.NodesEnd(); in++)
			{   
                const double nodal_area = in->FastGetSolutionStepValue(NODAL_AREA);
                const double density = in->FastGetSolutionStepValue(DENSITY);

                //finisgin calculation of Minverse*G*p
                array_1d<double,3>& press_proj = in->FastGetSolutionStepValue(PRESS_PROJ);
                if(nodal_area != 0.0)
                    press_proj *= (1.0/nodal_area);
                else
                    noalias(press_proj) = zero;

                if(in->FastGetSolutionStepValue(IS_STRUCTURE) != 1)
                {
                    const array_1d<double,3>& body_force = in->FastGetSolutionStepValue(BODY_FORCE);
                    array_1d<double,3>& disp = in->FastGetSolutionStepValue(DISPLACEMENT);
                    array_1d<double,3>& vel = in->FastGetSolutionStepValue(VELOCITY);
                    array_1d<double,3>& acc = in->FastGetSolutionStepValue(ACCELERATION);

                    const array_1d<double,3>& disp_old = in->FastGetSolutionStepValue(DISPLACEMENT,1);
                    const array_1d<double,3>& vel_old = in->FastGetSolutionStepValue(VELOCITY,1);
                    const array_1d<double,3>& acc_old = in->FastGetSolutionStepValue(ACCELERATION,1);


                    //calculating displacement increment (using newmark)
                    for(unsigned int j=0;j<TDim;j++)
                    {
                        delta_disp[j] = body_force[j] + press_proj[j]/density - nodal_area*acc[j];
                    }
                    delta_disp *= inverse_inertia_coeff;
                    
                    for(unsigned int j=0;j<TDim;j++) 
                    {
                        //updating displacements
                        disp[j] = disp[j] + delta_disp[j];

                        //getting the displacement increment in the step
                        displacement_step_increment[j] = disp[j] - disp_old[j]; 

                        //updating velocity
                        vel[j] = a1*displacement_step_increment[j] - a4*vel_old[j] - a5*acc_old[j];

                        //updating acceleration
                        acc[j] = a0*displacement_step_increment[j] - a2*vel_old[j] - a3*acc_old[j];
                    }

                    //updating coordinates
                    in->X() = in->X0() + disp[0];
                    in->Y() = in->Y0() + disp[1];
                    in->Z() = in->Z0() + disp[2];
                }

            }
            KRATOS_CATCH("");
        
	}	

        //************************************************************************
        //************************************************************************
        void ConstructLaplacianSystem(ModelPart& model_part)
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
                    mDofSet.push_back( it->pGetDof(PRESSURE) );
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
            unsigned int dof_position = (model_part.NodesBegin())->GetDofPosition(PRESSURE); 

            //building up the matrix graph row by row 
            unsigned int total_size = 0; 
            for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
            {
                unsigned int index_i = it->GetDof(PRESSURE,dof_position).EquationId(); 

		if(index_i < mEquationSystemSize && (it->GetValue(NEIGHBOUR_NODES)).size() != 0)
                {
                    WeakPointerVector< Node<3> >& neighb_nodes = it->GetValue(NEIGHBOUR_NODES); 
        
                    //filling the first neighbours list 
                    work_array.push_back(index_i); 
                    for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin(); i != neighb_nodes.end(); i++) 
                    { 
                        int index_j = i->GetDof(PRESSURE,dof_position).EquationId(); 
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
        void BuildAndSolveLaplacianSystem(ModelPart& model_part, typename TLinearSolver::Pointer linear_solver )
        {
            KRATOS_TRY
    
                TSparseSpace::SetToZero(mA);
                TSparseSpace::SetToZero(mb);
                TSparseSpace::SetToZero(mDx);

	            ProcessInfo& CurrentProcessInfo = model_part.GetProcessInfo();
                double dt = CurrentProcessInfo[DELTA_TIME];

    
                //**********************************
                //BUILD PHASE
                //**********************************
                            boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim> DN_DX; 
                            boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim+1> elemental_laplacian; 
                            array_1d<double,TDim+1> N; 
                            array_1d<unsigned int ,TDim+1> local_indices; 
                            array_1d<double,TDim+1> rhs_contribution; 
                            array_1d<double,TDim+1> dp_vector; 
    
                            //getting the dof position 
                            unsigned int dof_position = (model_part.NodesBegin())->GetDofPosition(PRESSURE); 
    
                            double aaa = 1.0/double(TDim+1.0); 
                            
                            array_1d<double,3> vg; 
                            for(ModelPart::ElementsContainerType::iterator i = model_part.ElementsBegin();  
                                    i!=model_part.ElementsEnd(); i++) 
                            {	 
    
                                    Geometry< Node<3> >& geom = i->GetGeometry();
    
                                    //calculating elemental values 
                                    double Volume; 
                                    GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume); 
    
                                    //determining elemental density and viscosity 
                                    double density = geom[0].FastGetSolutionStepValue(DENSITY); 
                                    for(int ii = 1; ii<TDim+1; ii++) 
                                    { 
                                            density += geom[ii].FastGetSolutionStepValue(DENSITY); 
                                    } 
                                    density *= aaa; 
    
                                    //finiding local indices 
                                    for(int ii = 0; ii<TDim+1; ii++) 
                                    { 
                                            local_indices[ii] = geom[ii].GetDof(PRESSURE,dof_position).EquationId(); 
    
                                            dp_vector[ii] = geom[ii].FastGetSolutionStepValue(PRESSURE) 
                                                            - geom[ii].FastGetSolutionStepValue(PRESSURE,1);
                                    } 
    
  
                                    //calculating laplacian LHS 
                                    double laplacian_coeff = dt/density; 
                                    noalias(elemental_laplacian) =  laplacian_coeff*prod(DN_DX,trans(DN_DX)); 
    
                    //calculating D*v
                                    double Gaux = 0.00; 
                                    for(int kk = 0; kk<TDim+1; kk++) 
                                    { 
                                            for(int tt = 0; tt<TDim; tt++) 
                                            { 
                                                    const array_1d<double,3>& v = geom[kk].FastGetSolutionStepValue(VELOCITY); 
                                                    Gaux += DN_DX(kk,tt)*v[tt]; 
                                            } 
                                    } 
                                    noalias(rhs_contribution) = - Gaux * N; 
    
                                    //applying dirichhlet
                                    noalias(rhs_contribution) -= prod(elemental_laplacian,dp_vector); 

                                    elemental_laplacian *= Volume;
                                    rhs_contribution *= Volume;

                                    AssembleLHS(mA,elemental_laplacian,local_indices);
                                    AssembleRHS(mb,rhs_contribution,local_indices);


	 
                            }
  
                //**********************************
                //SOLVE PHASE
                //**********************************
    //             KRATOS_WATCH(mA);
  //               KRATOS_WATCH(mb);
    //             KRATOS_WATCH(mDx);
                linear_solver->Solve(mA,mDx,mb);
               std::cout << *(linear_solver) << std::endl;
      //           KRATOS_WATCH(mDx);
    
                //**********************************
                //UPDATE PRESSURES
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
    
        void ClearLaplacianSystem(ModelPart& model_part)
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

    };

}  // namespace Kratos.

#endif // KRATOS_LARANGIAN_FLUID_UTILITIES_INCLUDED  defined 


