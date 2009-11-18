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
//   Date:                $Date: 2007-11-14 15:28:04 $
//   Revision:            $Revision: 1.1 $
//
//


#if !defined(KRATOS_BALLVERTEX_MESHMOVING3D_H_INCLUDED )
#define  KRATOS_BALLVERTEX_MESHMOVING3D_H_INCLUDED



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
#include "ale_application.h"
#include "ale_application.h"

namespace Kratos
{

    template<int TDim,class TSparseSpace,class TLinearSolver>
	class BallVertexMeshMoving3D
	{
		public:

	typedef typename TSparseSpace::MatrixType TSystemMatrixType;
	typedef typename TSparseSpace::VectorType TSystemVectorType;
        typedef PointerVectorSet< Dof<double> , IndexedObject> DofsArrayType;


		BallVertexMeshMoving3D(){};
		~BallVertexMeshMoving3D(){};

 



        //************************************************************************
        //************************************************************************
        void ConstructSystem(ModelPart& model_part)
        {
            KRATOS_TRY

            //loop over nodes with neighbours
	    mDofSet = DofsArrayType();
            //mDofSet.Resize(0);
	    //mDofSet.clear();
            mDofSet.reserve(model_part.Nodes().size() );

            int tot_nnz=0;
            int tot_row=0;
            for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
            {
                if( (it->GetValue(NEIGHBOUR_NODES)).size() != 0 )
                {
                    mDofSet.push_back( it->pGetDof(DISPLACEMENT_X) );
		    mDofSet.push_back( it->pGetDof(DISPLACEMENT_Y) );
		    mDofSet.push_back( it->pGetDof(DISPLACEMENT_Z) );
                    tot_row += TDim;
                    tot_nnz +=  ((it->GetValue(NEIGHBOUR_NODES)).size()+1)*TDim*TDim;
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
    
            //building up the matrix graph row by row 
            int total_size = 0; 
            for (typename ModelPart::NodesContainerType::iterator it=model_part.NodesBegin(); it!=model_part.NodesEnd(); ++it)
            {
                unsigned int index_i = it->GetDof(DISPLACEMENT_X).EquationId(); 
                unsigned int index_j = it->GetDof(DISPLACEMENT_Y).EquationId(); 
                unsigned int index_k = it->GetDof(DISPLACEMENT_Z).EquationId(); 

		//if(index_i < mEquationSystemSize && (it->GetValue(NEIGHBOUR_NODES)).size() != 0)
		if(index_i < mEquationSystemSize && index_j < mEquationSystemSize)
                {
                    WeakPointerVector< Node<3> >& neighb_nodes = it->GetValue(NEIGHBOUR_NODES); 
        
                    //filling the first neighbours list 
                    work_array.push_back(index_i);
		    work_array.push_back(index_j);  
		    work_array.push_back(index_k); 
                    for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin(); i != neighb_nodes.end(); i++) 
                    { 
                        unsigned int index_l = i->GetDof(DISPLACEMENT_X).EquationId();
			unsigned int index_r = i->GetDof(DISPLACEMENT_Y).EquationId();
			unsigned int index_s = i->GetDof(DISPLACEMENT_Z).EquationId();
                        if(index_l < mEquationSystemSize)
                            work_array.push_back(index_l);
                        if(index_r < mEquationSystemSize)
                            work_array.push_back(index_r);
                        if(index_s < mEquationSystemSize)
                            work_array.push_back(index_s);
                    } 
                                    
                    //sorting the indices and elminating the duplicates 
                    std::sort(work_array.begin(),work_array.end()); 
                    unsigned int number_of_entries = work_array.size(); 
                                    
                    //filling up the matrix 
                    for(unsigned int j=0; j<number_of_entries; j++) 
                    { 
                        mA.push_back(index_i,work_array[j] , 0.00);
                    } 
                    for(unsigned int j=0; j<number_of_entries; j++) 
                    { 
			mA.push_back(index_j,work_array[j] , 0.00);  
                    } 
                    for(unsigned int j=0; j<number_of_entries; j++) 
                    { 
			mA.push_back(index_k,work_array[j] , 0.00);  
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
        void BuildAndSolveSystem(ModelPart& model_part, typename TLinearSolver::Pointer linear_solver )
        {
            KRATOS_TRY
    
                TSparseSpace::SetToZero(mA);
                TSparseSpace::SetToZero(mb);
                TSparseSpace::SetToZero(mDx);

// 	            ProcessInfo& CurrentProcessInfo = model_part.GetProcessInfo();
         //       double dt = CurrentProcessInfo[DELTA_TIME];

    
                //**********************************
                //BUILD PHASE
                //**********************************
    
     
                 //           double aaa = 1.0/double(TDim+1.0); 

				boost::numeric::ublas::bounded_matrix<double,(TDim+1)*TDim,(TDim+1)*TDim>  K_matrix;
				array_1d<double ,(TDim+1)*TDim> rhs_vector;
				array_1d<double ,(TDim+1)*TDim> disps;
				array_1d<unsigned int ,(TDim+1)*TDim> local_indices; 

                            array_1d<double,3> vg; 
                            for(ModelPart::ElementsContainerType::iterator i = model_part.ElementsBegin();  
                                    i!=model_part.ElementsEnd(); i++) 
                            {	 
    
                                Geometry< Node<3> >& geom = i->GetGeometry();

				BallVertex3D( geom,  K_matrix    );
 
				const array_1d<double,3>& ddd = geom[0].FastGetSolutionStepValue(DISPLACEMENT);
				disps[0] = ddd[0];
				disps[1] = ddd[1];
				disps[2] = ddd[2];
				const array_1d<double,3>& bbb = geom[1].FastGetSolutionStepValue(DISPLACEMENT);
				disps[3] = bbb[0];
				disps[4] = bbb[1];
				disps[5] = bbb[2];
				const array_1d<double,3>& ccc = geom[2].FastGetSolutionStepValue(DISPLACEMENT);
				disps[6] = ccc[0];
				disps[7] = ccc[1];
				disps[8] = ccc[2];
				const array_1d<double,3>& aaa = geom[3].FastGetSolutionStepValue(DISPLACEMENT);
				disps[9] = aaa[0];
				disps[10] = aaa[1];
				disps[11] = aaa[2];

				for(int ii=0; ii < TDim+1; ii++){
					unsigned int base = ii*TDim;
					local_indices[base]=geom[ii].GetDof(DISPLACEMENT_X).EquationId();
					local_indices[base+1]=geom[ii].GetDof(DISPLACEMENT_Y).EquationId();
					local_indices[base+2]=geom[ii].GetDof(DISPLACEMENT_Z).EquationId();
				}


				//noalias(rhs_vector) = ZeroVector((TDim+1)*TDim);
				noalias(rhs_vector) = -prod(K_matrix,disps);
				//KRATOS_WATCH(K_matrix);
				//KRATOS_WATCH(disps);

                                    AssembleLHS(mA,K_matrix,local_indices);
                                    AssembleRHS(mb,rhs_vector,local_indices);

	 
                            }

                //**********************************
                //SOLVE PHASE
                //**********************************
                 //KRATOS_WATCH(mA.size1());
	 	 //KRATOS_WATCH(mA.size2());
  		 //KRATOS_WATCH(mA.nnz());
		 //KRATOS_WATCH(TSparseSpace::TwoNorm(mb));
           
   		//KRATOS_WATCH(mb);
                // KRATOS_WATCH(mDx);
               std::cout << *(linear_solver) << std::endl;
                linear_solver->Solve(mA,mDx,mb);
               std::cout << *(linear_solver) << std::endl;
      //           KRATOS_WATCH(mDx);
    
                //**********************************
                //UPDATE DISPLACEMENTS
                //**********************************
                for(typename DofsArrayType::iterator i_dof = mDofSet.begin() ; i_dof != mDofSet.end() ; ++i_dof)
                {
                    if(i_dof->IsFree())
                    {
                        i_dof->GetSolutionStepValue() += mDx[i_dof->EquationId()];
                    }
                }

// 		//update nodal coordinates
//                 for(ModelPart::NodesContainerType::iterator i = model_part.NodesBegin();  
//                                     i!=model_part.NodesEnd(); i++) 
//                 {
// 			const array_1d<double,3>& disp = i->FastGetSolutionStepValue(DISPLACEMENT);
// 			i->X() = i->X0() + disp[0];
// 			i->Y() = i->Y0() + disp[1];
// 			i->Z() = i->Z0() + disp[2];
// 
// //calculate mesh vel
// 		}
   
            KRATOS_CATCH("")
        }
    
        void ClearSystem()
        {
            KRATOS_TRY
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
			const boost::numeric::ublas::bounded_matrix<double,(TDim+1)*TDim,(TDim+1)*TDim> & LHS_Contribution,
			const array_1d<unsigned int,(TDim+1)*TDim>& EquationId
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
			const array_1d<double,(TDim+1)*TDim>& RHS_Contribution,
			const array_1d<unsigned int,(TDim+1)*TDim>& EquationId
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

		//*****************************************************************************	
		void BallVertex3D( const Geometry< Node<3> >& geom, boost::numeric::ublas::bounded_matrix<double,(TDim+1)*TDim,(TDim+1)*TDim>& K_matrix    ){

			array_1d<double,4> x,y,z;

			x[0] = geom[0].X();
			y[0] = geom[0].Y();
			z[0] = geom[0].Z();
			x[1] = geom[1].X();
			y[1] = geom[1].Y();
			z[1] = geom[1].Z();
			x[2] = geom[2].X();
			y[2] = geom[2].Y();
			z[2] = geom[2].Z();
			x[3] = geom[3].X();
			y[3] = geom[3].Y();
			z[3] = geom[3].Z();

			noalias(K_matrix) = ZeroMatrix(12,12);

			boost::numeric::ublas::bounded_matrix<int,6,2> index;
			index(0,0)=0;
			index(0,1)=1;
			index(1,0)=0;
			index(1,1)=2;
			index(2,0)=0;
			index(2,1)=3;
			index(3,0)=1;
			index(3,1)=2;
			index(4,0)=1;
			index(4,1)=3;
			index(5,0)=2;
			index(5,1)=3;

 			boost::numeric::ublas::bounded_matrix<double,3,3> vij_mat;
			double invlijq;
			double invlij;
			int i;
			int j;
			for (unsigned int k=0; k<6; k++)
			{
				i=index(k,0);
				j=index(k,1);

				// edge ij e ji	
				invlijq=1.0/(pow((x[j]-x[i]),2)+pow((y[j]-y[i]),2)+pow((z[j]-z[i]),2));
				invlij=1.0/sqrt(pow((x[j]-x[i]),2)+pow((y[j]-y[i]),2)+pow((z[j]-z[i]),2));
				
				vij_mat(0,0)=pow((x[j]-x[i]),2);
				vij_mat(0,1)=((x[j]-x[i])*(y[j]-y[i]));
				vij_mat(0,2)=((x[j]-x[i])*(z[j]-z[i]));
				vij_mat(1,0)=((x[j]-x[i])*(y[j]-y[i]));
				vij_mat(1,1)=pow((y[j]-y[i]),2);
				vij_mat(1,2)=((x[j]-x[i])*(z[j]-z[i]));
				vij_mat(2,0)=((x[j]-x[i])*(z[j]-z[i]));
				vij_mat(2,1)=((y[j]-y[i])*(z[j]-z[i]));
				vij_mat(2,2)=pow((z[j]-z[i]),2);
	
				vij_mat=vij_mat*invlijq;
	
				K_matrix(3*i,3*i)+=invlij*vij_mat(0,0);
				K_matrix(3*i,3*i+1)+=invlij*vij_mat(0,1);
				K_matrix(3*i,3*i+2)+=invlij*vij_mat(0,2);
				K_matrix(3*i+1,3*i)+=invlij*vij_mat(1,0);
				K_matrix(3*i+1,3*i+1)+=invlij*vij_mat(1,1);
				K_matrix(3*i+1,3*i+2)+=invlij*vij_mat(1,2);
				K_matrix(3*i+2,3*i)+=invlij*vij_mat(2,0);
				K_matrix(3*i+2,3*i+1)+=invlij*vij_mat(2,1);
				K_matrix(3*i+2,3*i+2)+=invlij*vij_mat(2,2);

				K_matrix(3*i,3*j)+=-invlij*vij_mat(0,0);
				K_matrix(3*i,3*j+1)+=-invlij*vij_mat(0,1);
				K_matrix(3*i,3*j+2)+=-invlij*vij_mat(0,2);
				K_matrix(3*i+1,3*j)+=-invlij*vij_mat(1,0);
				K_matrix(3*i+1,3*j+1)+=-invlij*vij_mat(1,1);
				K_matrix(3*i+1,3*j+2)+=-invlij*vij_mat(1,2);
				K_matrix(3*i+2,3*j)+=-invlij*vij_mat(2,0);
				K_matrix(3*i+2,3*j+1)+=-invlij*vij_mat(2,1);
				K_matrix(3*i+2,3*j+2)+=-invlij*vij_mat(2,2);

				K_matrix(3*j,3*i)+=-invlij*vij_mat(0,0);
				K_matrix(3*j,3*i+1)+=-invlij*vij_mat(0,1);
				K_matrix(3*j,3*i+2)+=-invlij*vij_mat(0,2);
				K_matrix(3*j+1,3*i)+=-invlij*vij_mat(1,0);
				K_matrix(3*j+1,3*i+1)+=-invlij*vij_mat(1,1);
				K_matrix(3*j+1,3*i+2)+=-invlij*vij_mat(1,2);
				K_matrix(3*j+2,3*i)+=-invlij*vij_mat(2,0);
				K_matrix(3*j+2,3*i+1)+=-invlij*vij_mat(2,1);
				K_matrix(3*j+2,3*i+2)+=-invlij*vij_mat(2,2);
	
				K_matrix(3*j,3*j)+=invlij*vij_mat(0,0);
				K_matrix(3*j,3*j+1)+=invlij*vij_mat(0,1);
				K_matrix(3*j,3*j+2)+=invlij*vij_mat(0,2);
				K_matrix(3*j+1,3*j)+=invlij*vij_mat(1,0);
				K_matrix(3*j+1,3*j+1)+=invlij*vij_mat(1,1);
				K_matrix(3*j+1,3*j+2)+=invlij*vij_mat(1,2);
				K_matrix(3*j+2,3*j)+=invlij*vij_mat(2,0);
				K_matrix(3*j+2,3*j+1)+=invlij*vij_mat(2,1);
				K_matrix(3*j+2,3*j+2)+=invlij*vij_mat(2,2);
 	
 			}

/////////////////////////////////////////////////////////////////////////////
			// edge 1p p1
			array_1d<double,3> v24;
			array_1d<double,3> v23;
			array_1d<double,3> n1p;
			double xp,yp,zp;
 			double area=0.0;
			triangle_area(x[1],x[2],x[3],y[1],y[2],y[3],z[1],z[2],z[3],area);

			v24[0]=x[3]-x[1];
			v24[1]=y[3]-y[1];
			v24[2]=z[3]-z[1];
			v24=v24/sqrt(pow(v24[0],2)+pow(v24[1],2)+pow(v24[2],2));
			v23[0]=x[2]-x[1];
			v23[1]=y[2]-y[1];
			v23[2]=z[2]-z[1];
			v23=v23/sqrt(pow(v23[0],2)+pow(v23[1],2)+pow(v23[2],2));

			n1p[0]=v24[1]*v23[2]-v24[2]*v23[1];
			n1p[1]=v24[2]*v23[0]-v24[0]*v23[2];
			n1p[2]=v24[0]*v23[1]-v24[1]*v23[0];
			n1p=n1p/sqrt(pow(n1p[0],2)+pow(n1p[1],2)+pow(n1p[2],2));

			boost::numeric::ublas::bounded_matrix<double,3,3> n1p_mat;
			
			n1p_mat(0,0)=pow(n1p[0],2);
			n1p_mat(0,1)=(n1p[0]*n1p[1]);
			n1p_mat(0,2)=(n1p[0]*n1p[2]);
			n1p_mat(1,0)=(n1p[1]*n1p[0]);
			n1p_mat(1,1)=pow(n1p[1],2);
			n1p_mat(1,2)=(n1p[1]*n1p[2]);
			n1p_mat(2,0)=(n1p[2]*n1p[0]);
			n1p_mat(2,1)=(n1p[2]*n1p[1]);
			n1p_mat(2,2)=pow(n1p[2],2);

			xp=x[0]-(n1p_mat(0,0)*(x[0]-x[1])+n1p_mat(0,1)*(y[0]-y[1])+n1p_mat(0,2)*(z[0]-z[1]));
			yp=y[0]-(n1p_mat(1,0)*(x[0]-x[1])+n1p_mat(1,1)*(y[0]-y[1])+n1p_mat(1,2)*(z[0]-z[1]));
			zp=z[0]-(n1p_mat(2,0)*(x[0]-x[1])+n1p_mat(2,1)*(y[0]-y[1])+n1p_mat(2,2)*(z[0]-z[1]));

			boost::numeric::ublas::bounded_matrix<double,3,3> v1p_mat;
			v1p_mat(0,0)=pow((xp-x[0]),2);
			v1p_mat(0,1)=((xp-x[0])*(yp-y[0]));
			v1p_mat(0,2)=((xp-x[0])*(zp-z[0]));
			v1p_mat(1,0)=((xp-x[0])*(yp-y[0]));
			v1p_mat(1,1)=pow((yp-y[0]),2);
			v1p_mat(1,2)=((yp-y[0])*(zp-z[0]));
			v1p_mat(2,0)=((xp-x[0])*(zp-z[0]));
			v1p_mat(2,1)=((yp-y[0])*(zp-z[0]));
			v1p_mat(2,2)=pow((zp-z[0]),2);
			
			v1p_mat=v1p_mat/(pow(xp-x[0],2)+pow(yp-y[0],2)+pow(zp-z[0],2));
			double invl1p=1.0/sqrt(pow(xp-x[0],2)+pow(yp-y[0],2)+pow(zp-z[0],2));
			
			double areap34=0.0;
			triangle_area(xp,x[2],x[3],yp,y[2],y[3],zp,z[2],z[3],areap34);
			double areap23=0.0;
			triangle_area(xp,x[1],x[2],yp,y[1],y[2],zp,z[1],z[2],areap23);
			double L2=areap34/area;
			double L4=areap23/area;
			double L3=1.0-L2-L4;

/////////////////////////////////////////////////////////////////////////////
			// edge 2r r2
			array_1d<double,3> v14;
			array_1d<double,3> v13;
			array_1d<double,3> n2r;
			double xr,yr,zr;
 			double area134=0.0;
			triangle_area(x[0],x[2],x[3],y[0],y[2],y[3],z[0],z[2],z[3],area134);

			v14[0]=x[3]-x[0];
			v14[1]=y[3]-y[0];
			v14[2]=z[3]-z[0];
			v14=v14/sqrt(pow(v14[0],2)+pow(v14[1],2)+pow(v14[2],2));
			v13[0]=x[2]-x[0];
			v13[1]=y[2]-y[0];
			v13[2]=z[2]-z[0];
			v13=v13/sqrt(pow(v13[0],2)+pow(v13[1],2)+pow(v13[2],2));

			n2r[0]=v14[1]*v13[2]-v14[2]*v13[1];
			n2r[1]=v14[2]*v13[0]-v14[0]*v13[2];
			n2r[2]=v14[0]*v13[1]-v14[1]*v13[0];
			n2r=n2r/sqrt(pow(n2r[0],2)+pow(n2r[1],2)+pow(n2r[2],2));

			boost::numeric::ublas::bounded_matrix<double,3,3> n2r_mat;
			
			n2r_mat(0,0)=pow(n2r[0],2);
			n2r_mat(0,1)=(n2r[0]*n2r[1]);
			n2r_mat(0,2)=(n2r[0]*n2r[2]);
			n2r_mat(1,0)=(n2r[1]*n2r[0]);
			n2r_mat(1,1)=pow(n2r[1],2);
			n2r_mat(1,2)=(n2r[1]*n2r[2]);
			n2r_mat(2,0)=(n2r[2]*n2r[0]);
			n2r_mat(2,1)=(n2r[2]*n2r[1]);
			n2r_mat(2,2)=pow(n2r[2],2);

			xr=x[1]-(n2r_mat(0,0)*(x[1]-x[0])+n2r_mat(0,1)*(y[1]-y[0])+n2r_mat(0,2)*(z[1]-z[0]));
			yr=y[1]-(n2r_mat(1,0)*(x[1]-x[0])+n2r_mat(1,1)*(y[1]-y[0])+n2r_mat(1,2)*(z[1]-z[0]));
			zr=z[1]-(n2r_mat(2,0)*(x[1]-x[0])+n2r_mat(2,1)*(y[1]-y[0])+n2r_mat(2,2)*(z[1]-z[0]));

			boost::numeric::ublas::bounded_matrix<double,3,3> v2r_mat;
			v2r_mat(0,0)=pow((xr-x[1]),2);
			v2r_mat(0,1)=((xr-x[1])*(yr-y[1]));
			v2r_mat(0,2)=((xr-x[1])*(zr-z[1]));
			v2r_mat(1,0)=((xr-x[1])*(yr-y[1]));
			v2r_mat(1,1)=pow((yr-y[1]),2);
			v2r_mat(1,2)=((yr-y[1])*(zr-z[1]));
			v2r_mat(2,0)=((xr-x[1])*(zr-z[1]));
			v2r_mat(2,1)=((yr-y[1])*(zr-z[1]));
			v2r_mat(2,2)=pow((zr-z[1]),2);
			
			v2r_mat=v2r_mat/(pow(xr-x[1],2)+pow(yr-y[1],2)+pow(zr-z[1],2));
			double invl2r=1.0/sqrt(pow(xr-x[1],2)+pow(yr-y[1],2)+pow(zr-z[1],2));;
			
			double arear34=0.0;
			triangle_area(xr,x[2],x[3],yr,y[2],y[3],zr,z[2],z[3],arear34);
			double arear14=0.0;
			triangle_area(xr,x[0],x[3],yr,y[0],y[3],zr,z[0],z[3],arear14);
			double M1=arear34/area134;
			double M3=arear14/area134;
			double M4=1.0-M1-M3;

/////////////////////////////////////////////////////////////////////////////
			// edge 3s 3s
			array_1d<double,3> v12;
			array_1d<double,3> n3s;
			double xs,ys,zs;
 			double area124=0.0;
			triangle_area(x[0],x[1],x[3],y[0],y[1],y[3],z[0],z[1],z[3],area124);

			v12[0]=x[1]-x[0];
			v12[1]=y[1]-y[0];
			v12[2]=z[1]-z[0];
			v12=v12/sqrt(pow(v12[0],2)+pow(v12[1],2)+pow(v12[2],2));

			n3s[0]=v12[1]*v14[2]-v12[2]*v14[1];
			n3s[1]=v12[2]*v14[0]-v12[0]*v14[2];
			n3s[2]=v12[0]*v14[1]-v12[1]*v14[0];
			n3s=n3s/sqrt(pow(n3s[0],2)+pow(n3s[1],2)+pow(n3s[2],2));

			boost::numeric::ublas::bounded_matrix<double,3,3> n3s_mat;
			
			n3s_mat(0,0)=pow(n3s[0],2);
			n3s_mat(0,1)=(n3s[0]*n3s[1]);
			n3s_mat(0,2)=(n3s[0]*n3s[2]);
			n3s_mat(1,0)=(n3s[1]*n3s[0]);
			n3s_mat(1,1)=pow(n3s[1],2);
			n3s_mat(1,2)=(n3s[1]*n3s[2]);
			n3s_mat(2,0)=(n3s[2]*n3s[0]);
			n3s_mat(2,1)=(n3s[2]*n3s[1]);
			n3s_mat(2,2)=pow(n3s[2],2);

			xs=x[2]-(n3s_mat(0,0)*(x[2]-x[0])+n3s_mat(0,1)*(y[2]-y[0])+n3s_mat(0,2)*(z[2]-z[0]));
			ys=y[2]-(n3s_mat(1,0)*(x[2]-x[0])+n3s_mat(1,1)*(y[2]-y[0])+n3s_mat(1,2)*(z[2]-z[0]));
			zs=z[2]-(n3s_mat(2,0)*(x[2]-x[0])+n3s_mat(2,1)*(y[2]-y[0])+n3s_mat(2,2)*(z[2]-z[0]));

			boost::numeric::ublas::bounded_matrix<double,3,3> v3s_mat;
			v3s_mat(0,0)=pow((xs-x[2]),2);
			v3s_mat(0,1)=((xs-x[2])*(ys-y[2]));
			v3s_mat(0,2)=((xs-x[2])*(zs-z[2]));
			v3s_mat(1,0)=((xs-x[2])*(ys-y[2]));
			v3s_mat(1,1)=pow((ys-y[2]),2);
			v3s_mat(1,2)=((ys-y[2])*(zs-z[2]));
			v3s_mat(2,0)=((xs-x[2])*(zs-z[2]));
			v3s_mat(2,1)=((ys-y[2])*(zs-z[2]));
			v3s_mat(2,2)=pow((zs-z[2]),2);
			
			v3s_mat=v3s_mat/(pow(xs-x[2],2)+pow(ys-y[2],2)+pow(zs-z[2],2));
			double invl3s=1.0/sqrt(pow(xs-x[2],2)+pow(ys-y[2],2)+pow(zs-z[2],2));
			
			double areas24=0.0;
			triangle_area(xs,x[1],x[3],ys,y[1],y[3],zs,z[1],z[3],areas24);
			double areas14=0.0;
			triangle_area(xs,x[0],x[3],ys,y[0],y[3],zs,z[0],z[3],areas14);
			double N1=areas24/area124;
			double N2=areas14/area124;
			double N4=1.0-N1-N2;


/////////////////////////////////////////////////////////////////////////////
			// edge 4t t4
			array_1d<double,3> v32;
			array_1d<double,3> v31;
			array_1d<double,3> n4t;
			double xt,yt,zt;
 			double area123=0.0;
			triangle_area(x[0],x[1],x[2],y[0],y[1],y[2],z[0],z[1],z[2],area123);

			v32[0]=x[2]-x[1];
			v32[1]=y[2]-y[1];
			v32[2]=z[2]-z[1];
			v32=v32/sqrt(pow(v32[0],2)+pow(v32[1],2)+pow(v32[2],2));
			v31[0]=x[2]-x[0];
			v31[1]=y[2]-y[0];
			v31[2]=z[2]-z[0];
			v31=v31/sqrt(pow(v31[0],2)+pow(v31[1],2)+pow(v31[2],2));

			n4t[0]=v32[1]*v31[2]-v32[2]*v31[1];
			n4t[1]=v32[2]*v31[0]-v32[0]*v31[2];
			n4t[2]=v32[0]*v31[1]-v32[1]*v31[0];
			n4t=n4t/sqrt(pow(n4t[0],2)+pow(n4t[1],2)+pow(n4t[2],2));

			boost::numeric::ublas::bounded_matrix<double,3,3> n4t_mat;
			
			n4t_mat(0,0)=pow(n4t[0],2);
			n4t_mat(0,1)=(n4t[0]*n4t[1]);
			n4t_mat(0,2)=(n4t[0]*n4t[2]);
			n4t_mat(1,0)=(n4t[1]*n4t[0]);
			n4t_mat(1,1)=pow(n4t[1],2);
			n4t_mat(1,2)=(n4t[1]*n4t[2]);
			n4t_mat(2,0)=(n4t[2]*n4t[0]);
			n4t_mat(2,1)=(n4t[2]*n4t[1]);
			n4t_mat(2,2)=pow(n4t[2],2);

			xt=x[3]-(n4t_mat(0,0)*(x[3]-x[2])+n4t_mat(0,1)*(y[3]-y[2])+n4t_mat(0,2)*(z[3]-z[2]));
			yt=y[3]-(n4t_mat(1,0)*(x[3]-x[2])+n4t_mat(1,1)*(y[3]-y[2])+n4t_mat(1,2)*(z[3]-z[2]));
			zt=z[3]-(n4t_mat(2,0)*(x[3]-x[2])+n4t_mat(2,1)*(y[3]-y[2])+n4t_mat(2,2)*(z[3]-z[2]));

			boost::numeric::ublas::bounded_matrix<double,3,3> v4t_mat;
			v4t_mat(0,0)=pow((xt-x[3]),2);
			v4t_mat(0,1)=((xt-x[3])*(yt-y[3]));
			v4t_mat(0,2)=((xt-x[3])*(zt-z[3]));
			v4t_mat(1,0)=((xt-x[3])*(yt-y[3]));
			v4t_mat(1,1)=pow((yt-y[3]),2);
			v4t_mat(1,2)=((yt-y[3])*(zt-z[3]));
			v4t_mat(2,0)=((xt-x[3])*(zt-z[3]));
			v4t_mat(2,1)=((yt-y[3])*(zt-z[3]));
			v4t_mat(2,2)=pow((zt-z[3]),2);
			
			v4t_mat=v4t_mat/(pow(xt-x[3],2)+pow(yt-y[3],2)+pow(zt-z[3],2));
			double invl4t=1.0/sqrt(pow(xt-x[3],2)+pow(yt-y[3],2)+pow(zt-z[3],2));
			
			double areat23=0.0;
			triangle_area(xt,x[1],x[2],yt,y[1],y[2],zt,z[1],z[2],areat23);
			double areat21=0.0;
			triangle_area(xt,x[1],x[0],yt,y[1],y[0],zt,z[1],z[0],areat21);
			double Q1=areat23/area123;
			double Q3=areat21/area123;
			double Q2=1.0-Q1-Q3;


			for (unsigned int ii=0; ii<3; ii++){
				for (unsigned int jj=0; jj<3; jj++){
					K_matrix(ii,jj)+=invl1p*v1p_mat(ii,jj) + pow(M1,2)*invl2r*v2r_mat(ii,jj) + pow(N1,2)*invl3s*v3s_mat(ii,jj) + pow(Q1,2)*invl4t*v4t_mat(ii,jj);
					K_matrix(ii,3+jj)+=-L2*invl1p*v1p_mat(ii,jj) - M1*invl2r*v2r_mat(ii,jj) + N1*N2*invl3s*v3s_mat(ii,jj) + Q1*Q2*invl4t*v4t_mat(ii,jj);
					K_matrix(ii,6+jj)+=-L3*invl1p*v1p_mat(ii,jj)+ M1*M3*invl2r*v2r_mat(ii,jj) - N1*invl3s*v3s_mat(ii,jj) + Q1*Q3*invl4t*v4t_mat(ii,jj);
					K_matrix(ii,9+jj)+=-L4*invl1p*v1p_mat(ii,jj)+ M1*M4*invl2r*v2r_mat(ii,jj)+ N1*N4*invl3s*v3s_mat(ii,jj) - Q1*invl4t*v4t_mat(ii,jj);

					K_matrix(3+ii,jj)+=-L2*invl1p*v1p_mat(ii,jj) - M1*invl2r*v2r_mat(ii,jj) + N1*N2*invl3s*v3s_mat(ii,jj) + Q1*Q2*invl4t*v4t_mat(ii,jj);
					K_matrix(3+ii,3+jj)+=pow(L2,2)*invl1p*v1p_mat(ii,jj) + invl2r*v2r_mat(ii,jj) + pow(N2,2)*invl3s*v3s_mat(ii,jj) + pow(Q2,2)*invl4t*v4t_mat(ii,jj);
					K_matrix(3+ii,6+jj)+=L2*L3*invl1p*v1p_mat(ii,jj) - M3*invl2r*v2r_mat(ii,jj) - N2*invl3s*v3s_mat(ii,jj) + Q3*Q2*invl4t*v4t_mat(ii,jj);
					K_matrix(3+ii,9+jj)+=L2*L4*invl1p*v1p_mat(ii,jj) - M4*invl2r*v2r_mat(ii,jj)+ N2*N4*invl3s*v3s_mat(ii,jj) - Q2*invl4t*v4t_mat(ii,jj);

					K_matrix(6+ii,jj)+=-L3*invl1p*v1p_mat(ii,jj) + M1*M3*invl2r*v2r_mat(ii,jj) - N1*invl3s*v3s_mat(ii,jj) + Q1*Q3*invl4t*v4t_mat(ii,jj);
					K_matrix(6+ii,3+jj)+=L3*L2*invl1p*v1p_mat(ii,jj) - M3*invl2r*v2r_mat(ii,jj) - N2*invl3s*v3s_mat(ii,jj) + Q3*Q2*invl4t*v4t_mat(ii,jj);
					K_matrix(6+ii,6+jj)+=pow(L3,2)*invl1p*v1p_mat(ii,jj) + pow(M3,2)*invl2r*v2r_mat(ii,jj) + invl3s*v3s_mat(ii,jj) + pow(Q3,2)*invl4t*v4t_mat(ii,jj);
					K_matrix(6+ii,9+jj)+=L3*L4*invl1p*v1p_mat(ii,jj) + M4*M3*invl2r*v2r_mat(ii,jj)- N4*invl3s*v3s_mat(ii,jj) - Q3*invl4t*v4t_mat(ii,jj);
					
					K_matrix(9+ii,jj)+=-L4*invl1p*v1p_mat(ii,jj) + M1*M4*invl2r*v2r_mat(ii,jj) + N1*N4*invl3s*v3s_mat(ii,jj) - Q1*invl4t*v4t_mat(ii,jj);
					K_matrix(9+ii,3+jj)+=L4*L2*invl1p*v1p_mat(ii,jj) - M4*invl2r*v2r_mat(ii,jj) + N2*N4*invl3s*v3s_mat(ii,jj) - Q2*invl4t*v4t_mat(ii,jj);
					K_matrix(9+ii,6+jj)+=L4*L3*invl1p*v1p_mat(ii,jj) + M4*M3*invl2r*v2r_mat(ii,jj) - N4*invl3s*v3s_mat(ii,jj) - Q3*invl4t*v4t_mat(ii,jj);
					K_matrix(9+ii,9+jj)+=pow(L4,2)*invl1p*v1p_mat(ii,jj) + pow(M4,2)*invl2r*v2r_mat(ii,jj) + pow(N4,2)*invl3s*v3s_mat(ii,jj) + invl4t*v4t_mat(ii,jj);

				}
			}

		}


		void triangle_area(double xa,double xb,double xc,double ya,double yb,double yc,double za,double zb,double zc, double& area){
		
		area=.5*sqrt(pow((xa*yc-xa*yb+xb*ya-xb*yc+xc*yb-xc*ya),2) + pow((ya*zb+yb*zc+yc*za-yc*zb-ya*zc-yb*za),2) + pow((za*xb+zb*xc+zc*xa-zc*xb-za*xc-zb*xa),2));
		}

    };

}  // namespace Kratos.

#endif // KRATOS_BALLVERTEX_MESHMOVING_H_INCLUDED  defined 

