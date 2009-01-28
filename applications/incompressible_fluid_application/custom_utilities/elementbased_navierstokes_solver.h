/*
==============================================================================
KratosIncompressibleFluidApplication 
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
//   Date:                $Date: 2009-01-13 15:39:56 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_ELEMENTBASED_NAVIER_STOKES_SOLVER_H_INCLUDED )
#define  KRATOS_ELEMENTBASED_NAVIER_STOKES_SOLVER_H_INCLUDED



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
#include "incompressible_fluid_application.h"


namespace Kratos
{
	
//variables to be used:
//VELOCITY
//NODAL_MASS
//DENSITY
//VISCOSITY
//PRESSURE
//PRESSURE_OLD_IT
//BODY_FORCE
//CONV_PROJ
//RHS_VECTOR 
//AUX_VECTOR


template<unsigned int TDim,class TSparseSpace,class TLinearSolver>
	class ElementBasedNavierStokesSolver
	{
		public:

		typedef typename TSparseSpace::MatrixType TSystemMatrixType;
		typedef typename TSparseSpace::VectorType TSystemVectorType;
	typedef PointerVectorSet< Dof<double> , IndexedObject> DofsArrayType;

		ElementBasedNavierStokesSolver(ModelPart& r_model_part)	: mr_model_part(r_model_part)
		{};
		~ElementBasedNavierStokesSolver(){};

	//************************************************************************
	//************************************************************************
	//constructs the structure for the pressure solver
	void ConstructSystemStructure()
	{
		KRATOS_TRY

		mDofSet = DofsArrayType();
		mDofSet.reserve(mr_model_part.Nodes().size()  );

		//counting number of non zeros and number of rows
		unsigned int tot_nnz=0;
		for (typename ModelPart::NodesContainerType::iterator it=mr_model_part.NodesBegin(); it!=mr_model_part.NodesEnd(); ++it)
		{
			if( (it->GetValue(NEIGHBOUR_NODES)).size() != 0 )
			{	
				mDofSet.push_back( it->pGetDof(PRESSURE) );
				tot_nnz +=  ((it->GetValue(NEIGHBOUR_NODES)).size()+1);
			}
		}

		//fills the DofList and give a unique progressive tag to each node 
		unsigned int free_id = 0;
		unsigned int fix_id = mDofSet.size();	
		//assiging Ids to pressures
		for (typename DofsArrayType::iterator dof_iterator = mDofSet.begin(); 
				   	dof_iterator != mDofSet.end(); ++dof_iterator)
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
		//int total_size = 0; 
		for (typename ModelPart::NodesContainerType::iterator it=mr_model_part.NodesBegin(); it!=mr_model_part.NodesEnd(); ++it)
		{
			unsigned int index_I = it->GetDof(PRESSURE).EquationId(); 
	
			if(index_I < mEquationSystemSize )
			{
				WeakPointerVector< Node<3> >& neighb_nodes = it->GetValue(NEIGHBOUR_NODES); 
		
				//filling the first neighbours list 
				work_array.push_back(index_I);
				for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin(); i != neighb_nodes.end(); i++) 
				{ 
					unsigned int index_J = i->GetDof(PRESSURE).EquationId(); 
					if(index_J < mEquationSystemSize )
						work_array.push_back(index_J); 
				} 
					
			//sorting 
			std::sort(work_array.begin(),work_array.end()); 
					
			//filling up the matrix 
			for(unsigned int j=0; j<work_array.size(); j++) 
				mA.push_back(index_I,work_array[j] , 0.00);
			
			//clearing the array for the next step 
			work_array.erase(work_array.begin(),work_array.end()); 
			}	
		}
	
		if(mDx.size() != mA.size1())
			mDx.resize(mA.size1(),false);
		if(mb.size() != mA.size1())
			mb.resize(mA.size1(),false);
		KRATOS_CATCH("")
	}



	void Clear()
	{
	KRATOS_TRY
	mDx.resize(0,false);
	mb.resize(0,false);
	mA.resize(0,0,false);

	KRATOS_CATCH("")
	}

	//advance in time the momentum equation using 4step runge kutta
	void SolveStep1()
	{
		KRATOS_TRY

		double lumping_factor = 1.0/double(TDim+1.0);
		array_1d<double,3> aux;	

		//getting delta time 
		double delta_t = mr_model_part.GetProcessInfo()[DELTA_TIME];
		
		//saving the old pressure
		for (typename ModelPart::NodesContainerType::iterator it=mr_model_part.NodesBegin(); it!=mr_model_part.NodesEnd(); ++it)
			it->FastGetSolutionStepValue(PRESSURE_OLD_IT) = it->FastGetSolutionStepValue(PRESSURE);
				
		//calculate nodal mass
		for (typename ModelPart::NodesContainerType::iterator it=mr_model_part.NodesBegin(); it!=mr_model_part.NodesEnd(); ++it)
			it->FastGetSolutionStepValue(NODAL_MASS) = 0.0;

		//computing the nodal mass
		for (typename ModelPart::ElementsContainerType::iterator it=mr_model_part.ElementsBegin(); it!=mr_model_part.ElementsEnd(); ++it)
		{
			//get the list of nodes of the element
			Geometry< Node<3> >& geom = it->GetGeometry();

			//compute elemental area
			double nodal_vol=0.0;
			if(TDim == 2) 
			        nodal_vol = GeometryUtils::CalculateVolume2D(geom);
			else 
				nodal_vol = GeometryUtils::CalculateVolume3D(geom);
			nodal_vol *= lumping_factor;

			for(unsigned int i=0; i<TDim+1 ; i++)
				geom[i].FastGetSolutionStepValue(NODAL_MASS) += nodal_vol * geom[i].FastGetSolutionStepValue(DENSITY);
		}
		
				
		//save velocity boundary conditions
		for (typename ModelPart::NodesContainerType::iterator it=mr_model_part.NodesBegin(); it!=mr_model_part.NodesEnd(); ++it)
		{
			if(it->pGetDof(VELOCITY_X)->IsFixed() == true)
			{
				mFixedVelocityDofSet.push_back( it->pGetDof(VELOCITY_X) );
				mFixedVelocityDofValues.push_back( it->pGetDof(VELOCITY_X)->GetSolutionStepValue() );
			}

			if(it->pGetDof(VELOCITY_Y)->IsFixed() == true)
			{
				mFixedVelocityDofSet.push_back( it->pGetDof(VELOCITY_Y) );
				mFixedVelocityDofValues.push_back( it->pGetDof(VELOCITY_Y)->GetSolutionStepValue() );
			}

			if(it->pGetDof(VELOCITY_Z)->IsFixed() == true)
			{
				mFixedVelocityDofSet.push_back( it->pGetDof(VELOCITY_Z) );
				mFixedVelocityDofValues.push_back( it->pGetDof(VELOCITY_Z)->GetSolutionStepValue() );
			}
		}
				
		//set WORK = VELOCITY of the old step
		for (typename ModelPart::NodesContainerType::iterator it=mr_model_part.NodesBegin(); it!=mr_model_part.NodesEnd(); ++it)
		{
			noalias(it->FastGetSolutionStepValue(AUX_VECTOR)) = it->FastGetSolutionStepValue(VELOCITY,1);
		}
				
		//first step of RUNGE KUTTA
		//set velocity to the velocity at the end of the old step (preserving its boundary conditions!)
		SetToZero(RHS_VECTOR,mr_model_part.Nodes());
		CalculateRHS(mr_model_part.Elements(),mr_model_part.Nodes());
		double one_sixt = 1.0/6.0; 
		for (typename ModelPart::NodesContainerType::iterator it=mr_model_part.NodesBegin(); it!=mr_model_part.NodesEnd(); ++it)
		{
			//AUX_VECTOR = AUX_VECTOR + delta_T/6 * 1/NODAL_MASS * RHS
			//VELOCITY = VELOCITY_old + delta_T/2 * 1/NODAL_MASS * RHS
			noalias(aux) = delta_t/(it->FastGetSolutionStepValue(NODAL_MASS)) * it->FastGetSolutionStepValue(RHS_VECTOR);
			noalias(it->FastGetSolutionStepValue(AUX_VECTOR)) += one_sixt * aux;
			noalias(it->FastGetSolutionStepValue(VELOCITY)) += 0.5 * aux;			
		}
		ApplyVelocityBoundaryConditions(mFixedVelocityDofSet,mFixedVelocityDofValues);
				
				
		//second step of Runge Kutta ... some of the coefficients changed
		SetToZero(RHS_VECTOR,mr_model_part.Nodes());
		CalculateRHS(mr_model_part.Elements(),mr_model_part.Nodes());
		double one_third = 1.0/3.0; 
		for (typename ModelPart::NodesContainerType::iterator it=mr_model_part.NodesBegin(); it!=mr_model_part.NodesEnd(); ++it)
		{
			//AUX_VECTOR = AUX_VECTOR + delta_T/3 * 1/NODAL_MASS * RHS
			//VELOCITY = VELOCITY_old + delta_T/2 * 1/NODAL_MASS * RHS
			noalias(aux) = delta_t/(it->FastGetSolutionStepValue(NODAL_MASS)) * it->FastGetSolutionStepValue(RHS_VECTOR);
			noalias(it->FastGetSolutionStepValue(AUX_VECTOR)) += one_third * aux;
			noalias(it->FastGetSolutionStepValue(VELOCITY)) += 0.5 * aux;		
		}
		ApplyVelocityBoundaryConditions(mFixedVelocityDofSet,mFixedVelocityDofValues);


		//third step of Runge Kutta ... some of the coefficients changed
		SetToZero(RHS_VECTOR,mr_model_part.Nodes());
		CalculateRHS(mr_model_part.Elements(),mr_model_part.Nodes());
		for (typename ModelPart::NodesContainerType::iterator it=mr_model_part.NodesBegin(); it!=mr_model_part.NodesEnd(); ++it)
		{
			//AUX_VECTOR = AUX_VECTOR + delta_T/3 * 1/NODAL_MASS * RHS
			//VELOCITY = VELOCITY_old + delta_T* 1/NODAL_MASS * RHS
			noalias(aux) = delta_t/(it->FastGetSolutionStepValue(NODAL_MASS)) * it->FastGetSolutionStepValue(RHS_VECTOR);
			noalias(it->FastGetSolutionStepValue(AUX_VECTOR)) += one_third * aux;
			noalias(it->FastGetSolutionStepValue(VELOCITY)) += aux;		
		}
		ApplyVelocityBoundaryConditions(mFixedVelocityDofSet,mFixedVelocityDofValues);
				
		//third step of Runge Kutta ... some of the coefficients changed
		SetToZero(RHS_VECTOR,mr_model_part.Nodes());
		CalculateRHS(mr_model_part.Elements(),mr_model_part.Nodes());
		for (typename ModelPart::NodesContainerType::iterator it=mr_model_part.NodesBegin(); it!=mr_model_part.NodesEnd(); ++it)
		{
			//AUX_VECTOR = AUX_VECTOR + delta_T/6 * 1/NODAL_MASS * RHS
			//VELOCITY = AUX_VSetToZero_VectorVarECTOR 
			noalias(aux) = delta_t/(it->FastGetSolutionStepValue(NODAL_MASS)) * it->FastGetSolutionStepValue(RHS_VECTOR);
			noalias(it->FastGetSolutionStepValue(AUX_VECTOR)) += one_sixt * aux;
			noalias(it->FastGetSolutionStepValue(VELOCITY)) = it->FastGetSolutionStepValue(AUX_VECTOR);			
		}
		ApplyVelocityBoundaryConditions(mFixedVelocityDofSet,mFixedVelocityDofValues);


		KRATOS_CATCH("");
	}
				
	//solve the pressure equation
	void SolveStep2(typename TLinearSolver::Pointer p_linear_solver)
	{
		KRATOS_TRY
				
		const double dt = mr_model_part.GetProcessInfo()[DELTA_TIME];
				
		boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim> DN_DX;
		boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim+1> L, lhs_contribution;
		array_1d<double,TDim+1> N;
		array_1d<unsigned int ,TDim+1> local_indices;
		array_1d<double ,TDim+1> el_pressures;
		array_1d<double,3> vel;
		array_1d<double,TDim> proj_aux;
		array_1d<double,TDim+1> rhs_contribution;
		double lumping_factor = 1.0/(TDim+1.0);
		array_1d<double,3> vg;
		
		//getting the dof position
		unsigned int dof_position = (mr_model_part.NodesBegin())->GetDofPosition(PRESSURE);
		
		for(ModelPart::ElementsContainerType::iterator i = mr_model_part.ElementsBegin(); 
				  i!=mr_model_part.ElementsEnd(); i++)
		{	

			Geometry< Node<3> >& geom = i->GetGeometry();

			//calculating elemental values
			double Volume;
			GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Volume);

			//determining elemental density and viscosity
			double density = geom[0].FastGetSolutionStepValue(DENSITY);
			double nu = geom[0].FastGetSolutionStepValue(VISCOSITY);
			for(unsigned int ii = 1; ii<TDim+1; ii++)
			{
				density += geom[ii].FastGetSolutionStepValue(DENSITY);
				nu += geom[ii].FastGetSolutionStepValue(VISCOSITY);
			}
			density *= lumping_factor;
			nu *= lumping_factor;

			//finiding local indices
			for(unsigned int ii = 0; ii<TDim+1; ii++)
			{
				local_indices[ii] = geom[ii].GetDof(PRESSURE,dof_position).EquationId();
				el_pressures[ii]  = geom[ii].FastGetSolutionStepValue(PRESSURE);
			}

			//calculating  "gauss point velocity"
			noalias(vg) = N[0]* ( geom[0].FastGetSolutionStepValue(FRACT_VEL) - 
					geom[0].FastGetSolutionStepValue(MESH_VELOCITY) );; //velocity on the gauss points
			for(unsigned int kk = 1; kk<TDim+1; kk++)
			{
				noalias(vg) +=	N[kk]* ( geom[kk].FastGetSolutionStepValue(FRACT_VEL) - 
							 geom[kk].FastGetSolutionStepValue(MESH_VELOCITY) );
			}

			//adding stabilization
			double h;

			if(TDim == 2)
				h = sqrt(2.0*Volume);
			else
				h = pow(6.00*Volume,0.33333333333);
			double c1 = 4.00;
			double c2 = 2.00;
			double norm_u = norm_2(vg);
			double tau = 1.00 / ( c1*nu/(h*h) + c2*norm_u/h );
			
			//calculating stabilization laplacian LHS
			double dt_coeff = dt; 
			double tau_dt_coeff = (tau + dt_coeff)/density;
			double tau_coeff = (tau )/density;
			noalias(L) =  prod(DN_DX,trans(DN_DX));
			

			//RHS = -Dv
			double Gaux = 0.00;
			for(unsigned int kk = 0; kk<TDim+1; kk++)
			{
				for(unsigned int tt = 0; tt<TDim; tt++)
				{
					const array_1d<double,3>& fv = geom[kk].FastGetSolutionStepValue(FRACT_VEL);
					Gaux += DN_DX(kk,tt)*fv[tt];
				}
			}
			noalias(rhs_contribution) = -Gaux * N;
			
			
			
			//RHS += D*proj 
			const array_1d<double,3>& proj_temp = geom[0].FastGetSolutionStepValue(PRESS_PROJ);
			for(unsigned int iii = 0; iii<TDim; iii++)
				proj_aux[iii] = N[0]*proj_temp[iii];
			for(unsigned int kk = 1; kk<TDim+1; kk++)
			{
				const array_1d<double,3>& proj_temp = geom[kk].FastGetSolutionStepValue(PRESS_PROJ);
				for(unsigned int iii = 0; iii<TDim; iii++)
					proj_aux[iii] += N[kk]*proj_temp[iii];
			}
			proj_aux *= tau_coeff;
			noalias(rhs_contribution) += prod(DN_DX , proj_aux);  
			
			//RHS += dt*L*pn - dt*L*pn1 - tau*L*pn1 
			// ==> RHS -=  tau*L*pn1
			noalias(rhs_contribution) -= tau_coeff * prod(DN_DX , proj_aux); 
			
			
			//completing the integration on the element
			rhs_contribution *= Volume;
			noalias(lhs_contribution) = (Volume * tau_dt_coeff)*L;

			//adding the stabilization laplacian to the system matrix and the stabilizing contribution to the LHS
			int active_size = mA.size1();
			for(unsigned int ii = 0; ii<TDim+1; ii++)
			{
				int i1 = local_indices[ii];
				if(i1 < active_size)
				{
					for(unsigned int jj = 0; jj<TDim+1; jj++)
					{
						int i2 = local_indices[jj];
						if(i2 < active_size)
							mA(i1,i2) += lhs_contribution(ii,jj);
					}
					mb[i1] += rhs_contribution[ii];
				}
			}	
		}	
			
		//solving the system
		p_linear_solver->Solve(mA,mDx,mb);
		
		//updating the pressure		
		for(typename DofsArrayType::iterator i_dof = mDofSet.begin() ; i_dof != mDofSet.end() ; ++i_dof)
				i_dof->GetSolutionStepValue() += mDx[i_dof->EquationId()];

		
		KRATOS_CATCH("");
	}
			
	//correct the velocity
	void SolveStep3()
	{
		KRATOS_TRY
				
		const double dt = mr_model_part.GetProcessInfo()[DELTA_TIME];
				
		//set to zero AUX vector
		SetToZero(AUX_VECTOR,mr_model_part.Nodes());
				
		//allocation of work space
		boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim> DN_DX;
		array_1d<double,TDim+1> N;
		array_1d<double,3> aux0, aux1, aux2, aux3; //this are sized to 3 even in 2D!!		
		//double lumping_factor = 1.0/double(TDim+1.0);
		
		
		//calculate the velocity correction and store it in AUX_VECTOR
		for (typename ModelPart::ElementsContainerType::iterator it=mr_model_part.ElementsBegin(); it!=mr_model_part.ElementsEnd(); ++it)
		{
			//get the list of nodes of the element
			Geometry< Node<3> >& geom = it->GetGeometry();

			double volume;
			GeometryUtils::CalculateGeometryData(geom, DN_DX, N, volume);
			
			//calculating average pressure variation
			double dp_avg = N[0]*(geom[0].FastGetSolutionStepValue(PRESSURE)-
					      geom[0].FastGetSolutionStepValue(PRESSURE_OLD_IT));
			for(unsigned int i=1;i<TDim+1;i++)
				dp_avg += N[i]*(geom[i].FastGetSolutionStepValue(PRESSURE)-
						geom[i].FastGetSolutionStepValue(PRESSURE_OLD_IT));
			dp_avg *= volume;
			
			//updating local variables prior to database writing access
			for(unsigned int i=0;i<TDim;i++)
			{
				aux0[i] += DN_DX(0,i)*dp_avg;
				aux1[i] += DN_DX(1,i)*dp_avg;
				aux2[i] += DN_DX(2,i)*dp_avg;
			}
			
			if(TDim == 3) for(unsigned int i=0;i<TDim;i++) aux3[i] += DN_DX(3,i)*dp_avg;
				
			//updating database
			for(unsigned int i=0;i<TDim+1;i++)
			  geom[i].FastGetSolutionStepValue(AUX_VECTOR) += aux1;
		}
		
		
		//correct the velocities
		for (typename ModelPart::NodesContainerType::iterator it=mr_model_part.NodesBegin(); it!=mr_model_part.NodesEnd(); ++it)
		{
			//VELOCITY = VELOCITY + dt * Minv * AUX_VECTOR
			double dt_Minv = dt / it->FastGetSolutionStepValue(NODAL_MASS);
			noalias(it->FastGetSolutionStepValue(VELOCITY))+=dt_Minv*it->FastGetSolutionStepValue(AUX_VECTOR);		
		}
		
				
		KRATOS_CATCH("");
	}
	
	void CalculateProjection( )
	{
		KRATOS_TRY
		if(TDim == 2)
			CalculateProjection2D(mr_model_part.Elements(),mr_model_part.Nodes());
		KRATOS_CATCH("");
	}
	
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	private:
		
	ModelPart& mr_model_part;
	
	DofsArrayType mFixedVelocityDofSet;
	std::vector<double> mFixedVelocityDofValues;
	
	//variables to be stored	
	unsigned int mEquationSystemSize;		
	TSystemVectorType mDx;
	TSystemVectorType mb;
	TSystemMatrixType mA;
	DofsArrayType mDofSet;

	void ApplyVelocityBoundaryConditions(DofsArrayType& mFixedVelocityDofSet,std::vector<double>& mFixedVelocityDofValues)
	{	
		KRATOS_TRY
				
		unsigned int i=0;
		for(typename DofsArrayType::iterator i_dof = mFixedVelocityDofSet.begin() ; i_dof != mFixedVelocityDofSet.end() ; ++i_dof)
		{
			i_dof->GetSolutionStepValue() = mFixedVelocityDofValues[i];
			i++;
		}
		
		KRATOS_CATCH("")
	}
	void SetToZero( Variable<array_1d<double,3> >& rVariable, ModelPart::NodesContainerType& rNodes)
        {
		KRATOS_TRY
		array_1d<double,3> zero = ZeroVector(3);
		for(ModelPart::NodesContainerType::iterator i = rNodes.begin(); i!=rNodes.end(); i++)
				noalias(i->FastGetSolutionStepValue(rVariable)) = zero;
		KRATOS_CATCH("")
	}
	
	void SetToZero( Variable<  double >& rVariable, ModelPart::NodesContainerType& rNodes)
        {
		KRATOS_TRY
		for(ModelPart::NodesContainerType::iterator i = rNodes.begin(); i!=rNodes.end(); i++)
				i->FastGetSolutionStepValue(rVariable) = 0.0;
		KRATOS_CATCH("")
	}
	
	

	void CalculateRHS(ModelPart::ElementsContainerType& rElements,
			  ModelPart::NodesContainerType& rNodes)
	{
		KRATOS_TRY
				
		if(TDim == 2) //2D case
			CalculateRHS2D(rElements,rNodes);
/*		else if(TDim == 3) //2D case
			CalculateRHS3D(rElements,rNodes);*/
		
		KRATOS_CATCH("");
	}
	
	void CalculateRHS2D(ModelPart::ElementsContainerType& rElements,
			  ModelPart::NodesContainerType& rNodes)
	{
		KRATOS_TRY
					
		//allocation of work space
		boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim+1> element_laplacian;
		boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim> DN_DX;
		array_1d<double,3> aux0, aux1, aux2, work; //this are sized to 3 even in 2D!!
		array_1d<double,TDim> vel_gauss;
		array_1d<double,TDim+1> N, u_DN;
		
		double lumping_factor = 1.0/double(TDim+1.0);
		
		
		//calculate RHS (and stabilization)
		for (typename ModelPart::ElementsContainerType::iterator it=mr_model_part.ElementsBegin(); it!=mr_model_part.ElementsEnd(); ++it)
		{
			//get the list of nodes of the element
			Geometry< Node<3> >& geom = it->GetGeometry();

			double volume;
			GeometryUtils::CalculateGeometryData(geom, DN_DX, N, volume);
			
			double volume_third = volume*lumping_factor;
			
			double c1 = 4.00; double c2 = 2.00;
			double h = sqrt(2.00*volume/3.0);
			double norm_u = 0.0; double tau=0.0;
			
			//getting the velocity on the nodes
			const array_1d<double,3>& 	v0 = geom[0].FastGetSolutionStepValue(VELOCITY);
			const array_1d<double,3>& 	w0 = geom[0].FastGetSolutionStepValue(MESH_VELOCITY);
			const array_1d<double,3>& 	proj0 = geom[0].FastGetSolutionStepValue(CONV_PROJ);
			const  array_1d<double,3>&	force0 = geom[0].FastGetSolutionStepValue(BODY_FORCE);
			double p0 = geom[0].FastGetSolutionStepValue(PRESSURE);
			//const double nu0 = geom[0].FastGetSolutionStepValue(VISCOSITY);
			//const double rho0 = geom[0].FastGetSolutionStepValue(DENSITY);
			
			const array_1d<double,3>& 	 v1 = geom[1].FastGetSolutionStepValue(VELOCITY);
			const array_1d<double,3>&	 w1 = geom[1].FastGetSolutionStepValue(MESH_VELOCITY);
			const array_1d<double,3>&	 proj1 = geom[1].FastGetSolutionStepValue(CONV_PROJ);
			const  array_1d<double,3>&	 force1 = geom[1].FastGetSolutionStepValue(BODY_FORCE);
			double p1 = geom[1].FastGetSolutionStepValue(PRESSURE);
			//const double nu1 = geom[1].FastGetSolutionStepValue(VISCOSITY);
			//const double rho1 = geom[1].FastGetSolutionStepValue(DENSITY);

			const array_1d<double,3>&  	v2    = geom[2].FastGetSolutionStepValue(VELOCITY);
			const array_1d<double,3>&  	w2     = geom[2].FastGetSolutionStepValue(MESH_VELOCITY);
			const array_1d<double,3>&  	proj2  = geom[2].FastGetSolutionStepValue(CONV_PROJ);
			const  array_1d<double,3>& 	force2 = geom[2].FastGetSolutionStepValue(BODY_FORCE);
			double p2 = geom[2].FastGetSolutionStepValue(PRESSURE);
			//const double nu2 = geom[2].FastGetSolutionStepValue(VISCOSITY);
			//const double rho2 = geom[2].FastGetSolutionStepValue(DENSITY);
			
			//******************************************************
			//body force  contribution to RHS (nodal integration)
			noalias(aux0) = volume_third * geom[0].FastGetSolutionStepValue(DENSITY) * force0;
			noalias(aux1) = volume_third * geom[1].FastGetSolutionStepValue(DENSITY) * force1;
			noalias(aux2) = volume_third * geom[2].FastGetSolutionStepValue(DENSITY) * force2;
			
			//pressure gradient contribution
			double p_avg = N[0]*p0 + N[1]*p1 + N[2]*p2;
			p_avg *= volume * lumping_factor; //it is divided and multiplied by density
			aux0[0] += DN_DX(0,0)*p_avg;
			aux0[1] += DN_DX(0,1)*p_avg;
			aux0[2] += DN_DX(0,2)*p_avg; 
			
			aux1[0] += DN_DX(1,0)*p_avg;
			aux1[1] += DN_DX(1,1)*p_avg;
			aux1[2] += DN_DX(1,2)*p_avg; 
			
			aux2[0] += DN_DX(2,0)*p_avg;
			aux2[1] += DN_DX(2,1)*p_avg;
			aux2[2] += DN_DX(2,2)*p_avg; 			
			
			//******************************************************
			//VISCOUS CONTRIBUTION TO THE RHS
			//calculating the elemental nu
			double nu=geom[0].FastGetSolutionStepValue(VISCOSITY);
			double density=geom[0].FastGetSolutionStepValue(DENSITY);
			for(unsigned int i=1;i<TDim+1;i++)
			{
				nu += geom[i].FastGetSolutionStepValue(VISCOSITY);
				density += geom[i].FastGetSolutionStepValue(DENSITY);
			}
			nu *= lumping_factor;
			density *= lumping_factor;
			// element_laplacian = Laplacian * nu; --> ONE GAUSS POINT
			// take care ... we are multiplying here by area and density
			noalias(element_laplacian) = (nu*volume*density) * prod(DN_DX,trans(DN_DX));
			
			//THS = - nu L * v
			noalias(aux0) -= element_laplacian(0,0)*v0; 
			noalias(aux0) -= element_laplacian(0,1)*v1;
			noalias(aux0) -= element_laplacian(0,2)*v2;
			
			noalias(aux1) -= element_laplacian(1,0)*v0; 
			noalias(aux1) -= element_laplacian(1,1)*v1;
			noalias(aux1) -= element_laplacian(1,2)*v2;
			
			noalias(aux2) -= element_laplacian(2,0)*v0; 
			noalias(aux2) -= element_laplacian(2,1)*v1;
			noalias(aux2) -= element_laplacian(2,2)*v2;
		
			
			//******************************************************
			//CONVECTIVE CONTRIBUTION (integrated using 3 gauss points)
			
			// ******************* GAUSS1 *********
			N[0]=0.5; N[1]=0.5; N[2]=0.0;
			vel_gauss[0] =  N[0]*(v0[0]-w0[0]) + N[1]*(v1[0]-w1[0]) +  N[2]*(v2[0]-w2[0]);
			vel_gauss[1] =  N[0]*(v0[1]-w0[1]) + N[1]*(v1[1]-w1[1]) +  N[2]*(v2[1]-w2[1]);

			//calculating parameter tau 		
			norm_u = vel_gauss[0]*vel_gauss[0] + vel_gauss[1]*vel_gauss[1];
			norm_u = sqrt(norm_u);
			tau = 1.00 / ( c1*nu/(h*h) + c2*norm_u/h );
			
			noalias(u_DN) = prod(DN_DX , vel_gauss);
			
			//here we sum both the convective part and the stabilization
			double temp = (N[0] + tau * u_DN[0] ) * volume * density * lumping_factor;
			noalias(aux0) -= (temp * u_DN[0]) * v0;
			noalias(aux0) -= (temp * u_DN[1]) * v1;
			noalias(aux0) -= (temp * u_DN[2]) * v2;
			
			temp = (N[1] + tau * u_DN[1] ) * volume * density * lumping_factor;
			noalias(aux1) -= (temp * u_DN[0]) * v0;
			noalias(aux1) -= (temp * u_DN[1]) * v1;
			noalias(aux1) -= (temp * u_DN[2]) * v2;
			
			temp = (N[2] + tau * u_DN[2] ) * volume * density * lumping_factor;
			noalias(aux2) -= (temp * u_DN[0]) * v0;
			noalias(aux2) -= (temp * u_DN[1]) * v1;
			noalias(aux2) -= (temp * u_DN[2]) * v2;
			
			//RHS += Suy * proj[component] 
			noalias(work) = N[0]*proj0;
			noalias(work) += N[1]*proj1;
			noalias(work) += N[2]*proj2;
			work *= volume * density * lumping_factor * tau;
			noalias(aux0) += u_DN[0] * work;
			noalias(aux1) += u_DN[1] * work;
			noalias(aux2) += u_DN[2] * work;
			
			// ******************* GAUSS2 *********
			N[0]=0.0; N[1]=0.5; N[2]=0.5;
			vel_gauss[0] =  N[0]*(v0[0]-w0[0]) + N[1]*(v1[0]-w1[0]) +  N[2]*(v2[0]-w2[0]);
			vel_gauss[1] =  N[0]*(v0[1]-w0[1]) + N[1]*(v1[1]-w1[1]) +  N[2]*(v2[1]-w2[1]);

			//calculating parameter tau 		
			norm_u = vel_gauss[0]*vel_gauss[0] + vel_gauss[1]*vel_gauss[1];
			norm_u = sqrt(norm_u);
			tau = 1.00 / ( c1*nu/(h*h) + c2*norm_u/h );
			
			noalias(u_DN) = prod(DN_DX , vel_gauss);
			
			//here we sum both the convective part and the stabilization
			temp = (N[0] + tau * u_DN[0] ) * volume * density * lumping_factor;
			noalias(aux0) -= (temp * u_DN[0]) * v0;
			noalias(aux0) -= (temp * u_DN[1]) * v1;
			noalias(aux0) -= (temp * u_DN[2]) * v2;
			
			temp = (N[1] + tau * u_DN[1] ) * volume * density * lumping_factor;
			noalias(aux1) -= (temp * u_DN[0]) * v0;
			noalias(aux1) -= (temp * u_DN[1]) * v1;
			noalias(aux1) -= (temp * u_DN[2]) * v2;
			
			temp = (N[2] + tau * u_DN[2] ) * volume * density * lumping_factor;
			noalias(aux2) -= (temp * u_DN[0]) * v0;
			noalias(aux2) -= (temp * u_DN[1]) * v1;
			noalias(aux2) -= (temp * u_DN[2]) * v2;
			
			//RHS += Suy * proj[component] 
			noalias(work) = N[0]*proj0;
			noalias(work) += N[1]*proj1;
			noalias(work) += N[2]*proj2;
			work *= volume * density * lumping_factor * tau;
			noalias(aux0) += u_DN[0] * work;
			noalias(aux1) += u_DN[1] * work;
			noalias(aux2) += u_DN[2] * work;
			
			// ******************* GAUSS3 *********
			N[0]=0.5; N[1]=0.0; N[2]=0.5;
			vel_gauss[0] =  N[0]*(v0[0]-w0[0]) + N[1]*(v1[0]-w1[0]) +  N[2]*(v2[0]-w2[0]);
			vel_gauss[1] =  N[0]*(v0[1]-w0[1]) + N[1]*(v1[1]-w1[1]) +  N[2]*(v2[1]-w2[1]);

			//calculating parameter tau 		
			norm_u = vel_gauss[0]*vel_gauss[0] + vel_gauss[1]*vel_gauss[1];
			norm_u = sqrt(norm_u);
			tau = 1.00 / ( c1*nu/(h*h) + c2*norm_u/h );
			
			noalias(u_DN) = prod(DN_DX , vel_gauss);
			
			//here we sum both the convective part and the stabilization
			temp = (N[0] + tau * u_DN[0] ) * volume * density * lumping_factor;
			noalias(aux0) -= (temp * u_DN[0]) * v0;
			noalias(aux0) -= (temp * u_DN[1]) * v1;
			noalias(aux0) -= (temp * u_DN[2]) * v2;
			
			temp = (N[1] + tau * u_DN[1] ) * volume * density * lumping_factor;
			noalias(aux1) -= (temp * u_DN[0]) * v0;
			noalias(aux1) -= (temp * u_DN[1]) * v1;
			noalias(aux1) -= (temp * u_DN[2]) * v2;
			
			temp = (N[2] + tau * u_DN[2] ) * volume * density * lumping_factor;
			noalias(aux2) -= (temp * u_DN[0]) * v0;
			noalias(aux2) -= (temp * u_DN[1]) * v1;
			noalias(aux2) -= (temp * u_DN[2]) * v2;
			
			//RHS += Suy * proj[component] 
			noalias(work) = N[0]*proj0;
			noalias(work) += N[1]*proj1;
			noalias(work) += N[2]*proj2;
			work *= volume * density * lumping_factor * tau;
			noalias(aux0) += u_DN[0] * work;
			noalias(aux1) += u_DN[1] * work;
			noalias(aux2) += u_DN[2] * work;
			
			//************************************************************	
			//************************************************************	
			// this is the only point with writing access to the database
			//************************************************************	
			//************************************************************			
			//writing on the database
			geom[0].FastGetSolutionStepValue(RHS_VECTOR) += aux0;
			geom[1].FastGetSolutionStepValue(RHS_VECTOR) += aux1;
			geom[2].FastGetSolutionStepValue(RHS_VECTOR) += aux2;
			
			
			
		}		
				
		KRATOS_CATCH("");
	}

	

	
	
	void CalculateProjection2D(ModelPart::ElementsContainerType& rElements, 				   ModelPart::NodesContainerType& rNodes )
	{
		KRATOS_TRY
				
		  //ProcessInfo& rCurrentProcessInfo = mr_model_part.GetProcessInfo();
		array_1d<double,3> zero = ZeroVector(3);

			//first of all set to zero the nodal variables to be updated nodally
		for(ModelPart::NodeIterator i = mr_model_part.NodesBegin() ; 
				  i != mr_model_part.NodesEnd() ; ++i)
		{
			noalias((i)->FastGetSolutionStepValue(PRESS_PROJ)) = zero;
			noalias((i)->FastGetSolutionStepValue(CONV_PROJ) ) = zero;
		}

		
		//allocation of work space
		boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim+1> element_laplacian;
		boost::numeric::ublas::bounded_matrix<double,TDim+1,TDim> DN_DX;
		array_1d<double,3> aux0, aux1, aux2, press_aux0, press_aux1, press_aux2; //this are sized to 3 even in 2D!!
		array_1d<double,TDim> vel_gauss;
		array_1d<double,TDim+1> N, u_DN;
		
		//double lumping_factor = 1.0/double(TDim+1.0);
		
		//calculate RHS (and stabilization)
		for (typename ModelPart::ElementsContainerType::iterator it=mr_model_part.ElementsBegin(); it!=mr_model_part.ElementsEnd(); ++it)
		{
			//get the list of nodes of the element
			Geometry< Node<3> >& geom = it->GetGeometry();

			double volume;
			GeometryUtils::CalculateGeometryData(geom, DN_DX, N, volume);
		
			array_1d<double,3>& v0 = geom[0].FastGetSolutionStepValue(VELOCITY);
			array_1d<double,3>& w0 = geom[0].FastGetSolutionStepValue(MESH_VELOCITY);
			//array_1d<double,3>& press_proj0 = geom[0].FastGetSolutionStepValue(PRESS_PROJ);
			array_1d<double,3>& conv_proj0 = geom[0].FastGetSolutionStepValue(CONV_PROJ);
			double p0 = geom[0].FastGetSolutionStepValue(PRESSURE);
			const double rho0 = geom[0].FastGetSolutionStepValue(DENSITY);
			
			array_1d<double,3>& v1 = geom[1].FastGetSolutionStepValue(VELOCITY);
			array_1d<double,3>& w1 = geom[1].FastGetSolutionStepValue(MESH_VELOCITY);
			//array_1d<double,3>& press_proj1 = geom[1].FastGetSolutionStepValue(PRESS_PROJ);
			array_1d<double,3>& conv_proj1 = geom[1].FastGetSolutionStepValue(CONV_PROJ);
			double p1 = geom[1].FastGetSolutionStepValue(PRESSURE);
			const double rho1 = geom[1].FastGetSolutionStepValue(DENSITY);

			array_1d<double,3>& v2 = geom[2].FastGetSolutionStepValue(VELOCITY);
			array_1d<double,3>& w2 = geom[2].FastGetSolutionStepValue(MESH_VELOCITY);
			//array_1d<double,3>& press_proj2 = geom[2].FastGetSolutionStepValue(PRESS_PROJ);
			array_1d<double,3>& conv_proj2 = geom[2].FastGetSolutionStepValue(CONV_PROJ);
			double p2 = geom[2].FastGetSolutionStepValue(PRESSURE);
			const double rho2 = geom[2].FastGetSolutionStepValue(DENSITY);

			double density = 0.3333333333333333333333*(rho0 + rho1 + rho2 );

			//PRESSURE PROJECTION
			//note that here we calculate it "strong"
			vel_gauss[0] = DN_DX(0,0)*(p0) + DN_DX(1,0)*(p1) + DN_DX(2,0)*(p2);
			vel_gauss[1] = DN_DX(0,1)*(p0) + DN_DX(1,1)*(p1) + DN_DX(2,1)*(p2);
			vel_gauss *= volume;

			//press_proj += G*p
			press_aux0[0] += N[0]*vel_gauss[0]; 
			press_aux0[1] += N[0]*vel_gauss[1]; 
			press_aux0[2] = 0.0;

			press_aux1[0] += N[1]*vel_gauss[0]; 
			press_aux1[1] += N[1]*vel_gauss[1]; 
			press_aux1[2] = 0.0;

			press_aux2[0] += N[2]*vel_gauss[0]; 
			press_aux2[1] += N[2]*vel_gauss[1]; 
			press_aux2[2] = 0.0;
						
			//CONVECTIVE PROJECTION			
			// ******************* GAUSS1 ***************************
			N[0]=0.5; N[1]=0.5; N[2]=0.0;
		
			// vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) 
			//(note that the fractional step vel is used)
			vel_gauss[0] =  N[0]*(v0[0]-w0[0]) + N[1]*(v1[0]-w1[0]) +  N[2]*(v2[0]-w2[0]);
			vel_gauss[1] =  N[0]*(v0[1]-w0[1]) + N[1]*(v1[1]-w1[1]) +  N[2]*(v2[1]-w2[1]);

			//calculating convective auxiliary vector
			noalias(u_DN) = prod(DN_DX , vel_gauss);

			//attention changing the meaning of vel_gauss!!
			vel_gauss[0] = u_DN[0] * v0[0] + u_DN[1] * v1[0] + u_DN[2] * v2[0];
			vel_gauss[1] = u_DN[0] * v0[1] + u_DN[1] * v1[1] + u_DN[2] * v2[1];
			
			// conv_proj += C*u
			aux0[0] = N[0]*vel_gauss[0]; 
			aux0[1] = N[0]*vel_gauss[1];

			aux1[0] = N[1]*vel_gauss[0]; 
			aux1[1] = N[1]*vel_gauss[1];

			aux2[0] = N[2]*vel_gauss[0];
			aux2[1] = N[2]*vel_gauss[1];
			
			
			// ******************* GAUSS2 ***************************
			N[0]=0.0; N[1]=0.5; N[2]=0.5;
		
			// vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) 
			//(note that the fractional step vel is used)
			vel_gauss[0] =  N[0]*(v0[0]-w0[0]) + N[1]*(v1[0]-w1[0]) +  N[2]*(v2[0]-w2[0]);
			vel_gauss[1] =  N[0]*(v0[1]-w0[1]) + N[1]*(v1[1]-w1[1]) +  N[2]*(v2[1]-w2[1]);

			//calculating convective auxiliary vector
			noalias(u_DN) = prod(DN_DX , vel_gauss);

			//attention changing the meaning of vel_gauss!!
			vel_gauss[0] = u_DN[0] * v0[0] + u_DN[1] * v1[0] + u_DN[2] * v2[0];
			vel_gauss[1] = u_DN[0] * v0[1] + u_DN[1] * v1[1] + u_DN[2] * v2[1];
			
			// conv_proj += C*u
			aux0[0] += N[0]*vel_gauss[0]; 
			aux0[1] += N[0]*vel_gauss[1];

			aux1[0] += N[1]*vel_gauss[0]; 
			aux1[1] += N[1]*vel_gauss[1];

			aux2[0] += N[2]*vel_gauss[0];
			aux2[1] += N[2]*vel_gauss[1];
						
			// ******************* GAUSS3 ***************************
			N[0]=0.5; N[1]=0.0; N[2]=0.5;
		
			// vel_gauss = sum( N[i]*(vel[i]-mesh_vel[i]), i=0, number_of_points) 
			//(note that the fractional step vel is used)
			vel_gauss[0] =  N[0]*(v0[0]-w0[0]) + N[1]*(v1[0]-w1[0]) +  N[2]*(v2[0]-w2[0]);
			vel_gauss[1] =  N[0]*(v0[1]-w0[1]) + N[1]*(v1[1]-w1[1]) +  N[2]*(v2[1]-w2[1]);

			//calculating convective auxiliary vector
			noalias(u_DN) = prod(DN_DX , vel_gauss);

			//attention changing the meaning of vel_gauss!!
			vel_gauss[0] = u_DN[0] * v0[0] + u_DN[1] * v1[0] + u_DN[2] * v2[0];
			vel_gauss[1] = u_DN[0] * v0[1] + u_DN[1] * v1[1] + u_DN[2] * v2[1];
			
			// conv_proj += C*u
			aux0[0] += N[0]*vel_gauss[0]; 
			aux0[1] += N[0]*vel_gauss[1];

			aux1[0] += N[1]*vel_gauss[0]; 
			aux1[1] += N[1]*vel_gauss[1];

			aux2[0] += N[2]*vel_gauss[0];
			aux2[1] += N[2]*vel_gauss[1];
			
			
			
			
			//********************************************
			//access to database
			
			//press_proj += grad p
			noalias(conv_proj0) += press_aux0;
			noalias(conv_proj1) += press_aux1;
			noalias(conv_proj2) += press_aux2;
			

			// conv_proj += C*u
			double temp = volume*density*0.33333333333333333333333333;
			conv_proj0[0] += temp*aux0[0]; 
			conv_proj0[1] +=  temp*aux0[1];

			conv_proj1[0] +=  temp*aux1[0]; 
			conv_proj1[1] +=  temp*aux1[1];

			conv_proj2[0] +=  temp*aux2[0];
			conv_proj2[1] +=  temp*aux2[1];			
			
			
			
		}

			//solve nodally for the velocity
		double temp;
		for(ModelPart::NodeIterator i = mr_model_part.NodesBegin() ; 
				  i != mr_model_part.NodesEnd() ; ++i)
		{
			array_1d<double,3>& press_proj = (i)->FastGetSolutionStepValue(PRESS_PROJ);
			array_1d<double,3>& conv_proj = (i)->FastGetSolutionStepValue(CONV_PROJ);
			double A = (i)->FastGetSolutionStepValue(NODAL_MASS);
				
			temp = 1.00 / A;
			press_proj *= temp; 
			conv_proj *= temp; 
				
		}

				
				
		KRATOS_CATCH("");
	}
			
			
			
			
			
			
			
	};

}  // namespace Kratos.

#endif // KRATOS_ELEMENTBASED_NAVIER_STOKES_SOLVER_H_INCLUDED  defined 

