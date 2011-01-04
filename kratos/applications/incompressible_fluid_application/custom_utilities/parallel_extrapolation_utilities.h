//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_PARALLEL_EXTRAPOLATION_UTILITIES_H_INCLUDED )
#define  KRATOS_PARALLEL_EXTRAPOLATION_UTILITIES_H_INCLUDED



// System includes
#include <string>
#include <iostream> 


// External includes 


// Project includes
#include "includes/define.h"
#include "utilities/geometry_utilities.h"


namespace Kratos
{

  ///@name Kratos Globals
  ///@{ 
  
  ///@} 
  ///@name Type Definitions
  ///@{ 
  
  ///@} 
  ///@name  Enum's
  ///@{
      
  ///@}
  ///@name  Functions 
  ///@{
      
  ///@}
  ///@name Kratos Classes
  ///@{
  
  /// Short class definition.
  /** Detail class definition.
  */
  template< unsigned int TDim>
  class ParallelExtrapolationUtilities
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of ParallelExtrapolationUtilities
      KRATOS_CLASS_POINTER_DEFINITION(ParallelExtrapolationUtilities);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      ParallelExtrapolationUtilities(){};

      /// Destructor.
      virtual ~ParallelExtrapolationUtilities(){};
      
      ///Function to extrapolate the velocity from the interior of the level set domain.
      ///the extrapolation velocity is taken from the first layer of nodes INSIDE the level set domain
      ///@param rmodel_part is the ModelPart on which we will operate
      ///@param rDistanceVar is the Variable that we will use in calculating the distance
      ///@param rVelocityVar is the Variable being extrapolated
      ///@param rAreaVar is the Variable that we will use for L2 projections
      ///@param max_levels is the number of maximum "layers" of element that will be used in the calculation of the distances
      void ExtrapolateVelocity(ModelPart& rmodel_part,
			     const Variable<double>& rDistanceVar, 
			     const Variable< array_1d<double,3> >& rVelocityVar, 
			     const Variable<double>& rAreaVar, 
			     const unsigned int max_levels)
	{
	  KRATOS_TRY
	  
	  bool is_distributed = false;
	  if(rmodel_part.GetCommunicator().TotalProcesses() > 1)
	    is_distributed = true;
	  
	   //check that variables needed are in the model part
	   if(!(rmodel_part.NodesBegin()->SolutionStepsDataHas(rDistanceVar)) )
	     KRATOS_ERROR(std::logic_error,"distance Variable is not in the model part","")
	   if(!(rmodel_part.NodesBegin()->SolutionStepsDataHas(rVelocityVar)) )
	     KRATOS_ERROR(std::logic_error,"velocity Variable is not in the model part","")
	   if(!(rmodel_part.NodesBegin()->SolutionStepsDataHas(rAreaVar)) )
	     KRATOS_ERROR(std::logic_error,"Area Variable is not in the model part","")
	     
	   if(is_distributed == true)
	    if(!(rmodel_part.NodesBegin()->SolutionStepsDataHas(PARTITION_INDEX)) )
	     KRATOS_ERROR(std::logic_error,"PARTITION_INDEX Variable is not in the model part","")

	  //set as active the internal nodes
	   int node_size = rmodel_part.Nodes().size();
	  #pragma omp parallel for
	  for(int i = 0; i<node_size; i++)
	  {
	      ModelPart::NodesContainerType::iterator it=rmodel_part.NodesBegin()+i;
	      it->FastGetSolutionStepValue(rAreaVar) = 0.0;
	      const double& dist = it->FastGetSolutionStepValue(rDistanceVar);
	      if(dist < 0.0)
		it->SetValue(IS_VISITED,1.0);
	      else
		it->SetValue(IS_VISITED,0.0);
	   }

	   //now extrapolate the velocity var from all of the nodes that have at least one internal node_size
	   array_1d<double,TDim+1> visited;
	   array_1d<double,TDim+1> N;
	   array_1d<double,TDim> avg;
	   double lumping_factor = 1.0/double(TDim+1);
	   boost::numeric::ublas::bounded_matrix <double, TDim+1,TDim> DN_DX;
	   int elem_size = rmodel_part.Elements().size();
	   for(unsigned int level=0; level<max_levels; level++)
	   {
		#pragma omp parallel for private(DN_DX,N,visited,avg) firstprivate(lumping_factor,elem_size)
		for(int i = 0; i<elem_size; i++)
		{
		    PointerVector< Element>::iterator it=rmodel_part.ElementsBegin()+i;

		    Geometry<Node<3> >&geom = it->GetGeometry();
				    
		    for(unsigned int k=0; k<TDim+1; k++)
			    visited[k] = geom[k].GetValue(IS_VISITED);
					  
		    if(IsActive(visited))
		    {	     
			double Volume;
			GeometryUtils::CalculateGeometryData(geom,DN_DX,N,Volume);
			
			//calculated average value of "velocity"
			noalias(avg) = ZeroVector(3);
			double counter = 0.0;
			for(unsigned int k=0; k<TDim+1; k++)
			{
			  if(visited[k]>0.0)
			  {
			      noalias(avg) += geom[k].FastGetSolutionStepValue(rVelocityVar);
			      counter += 1.0;
			  }
			}
			avg /= counter;
			
			for(unsigned int k=0; k<TDim+1; k++)
			{
			  if(visited[k]==0.0)
			  {
			      geom[k].SetLock();
			      geom[k].FastGetSolutionStepValue(rVelocityVar) += avg*Volume*lumping_factor;
			      geom[k].FastGetSolutionStepValue(rAreaVar) += Volume*lumping_factor;
			      geom[k].UnSetLock();
			  }
			}
		    }
		}	

		//mpi sync variables
		if(is_distributed == false)
		{
		  int my_rank = rmodel_part.GetCommunicator().MyPID();
		    #pragma omp parallel for firstprivate(node_size,my_rank)
		    for(int i = 0; i<node_size; i++)
		    {
			    ModelPart::NodesContainerType::iterator it=rmodel_part.NodesBegin()+i;
			    if(it->GetValue(IS_VISITED) == 1.0 && my_rank != it->FastGetSolutionStepValue(PARTITION_INDEX))
			    {
			      it->FastGetSolutionStepValue(rAreaVar) = 0.0;
			      noalias( it->FastGetSolutionStepValue(rVelocityVar) ) = ZeroVector(3) ;
			    }
		    }
		    rmodel_part.GetCommunicator().AssembleCurrentData(rAreaVar);
		    rmodel_part.GetCommunicator().AssembleCurrentData(rVelocityVar);
		}
		
		#pragma omp parallel for firstprivate(node_size)
		for(int i = 0; i<node_size; i++)
		{
			ModelPart::NodesContainerType::iterator it=rmodel_part.NodesBegin()+i;
			const double& area = it->FastGetSolutionStepValue(rAreaVar);
			if(area > 1e-20 && it->GetValue(IS_VISITED)==0.0) //this implies that node was computed
			{
			  it->FastGetSolutionStepValue(rVelocityVar) /= area;
			  it->SetValue(IS_VISITED,1.0); //node is marked as already visitedz
			}

		}
	   }
	         
	     KRATOS_CATCH("")
	}

      ///Function to extrapolate the pressure projection from the interior of the level set domain.
      ///the extrapolation velocity is taken from the SECOND layer of nodes INSIDE the level set domain,
      ///that is, the nodes that correspond to elements cut by the free surface ARE NOT used as a source for the extrapolation
      ///note that the body force is assigned to the pressure projection in all of the nodes which are not calculated otherwise
      ///@param rmodel_part is the ModelPart on which we will operate
      ///@param rDistanceVar is the Variable that we will use in calculating the distance
      ///@param rProjVar is the Variable being extrapolated
      ///@param rAreaVar is the Variable that we will use for L2 projections
      ///@param max_levels is the number of maximum "layers" of element that will be used in the calculation of the distances
      void ExtrapolatePressureProjection(ModelPart& rmodel_part,
			     const Variable<double>& rDistanceVar, 
			     const Variable<array_1d<double,3> >& rProjVar, 
			     const Variable<double>& rPressureVar,
			     const Variable<double>& rAreaVar, 
			     const unsigned int max_levels)
	{
	  KRATOS_TRY
	  
	  bool is_distributed = false;
	  if(rmodel_part.GetCommunicator().TotalProcesses() > 1)
	    is_distributed = true;
	  
	   //check that variables needed are in the model part
	   if(!(rmodel_part.NodesBegin()->SolutionStepsDataHas(rDistanceVar)) )
	     KRATOS_ERROR(std::logic_error,"distance Variable is not in the model part","")
	   if(!(rmodel_part.NodesBegin()->SolutionStepsDataHas(rProjVar)) )
	     KRATOS_ERROR(std::logic_error,"Pressure Projection Variable is not in the model part","")
	   if(!(rmodel_part.NodesBegin()->SolutionStepsDataHas(BODY_FORCE)) )
	     KRATOS_ERROR(std::logic_error,"BODY_FORCE Variable is not in the model part","")
	   if(!(rmodel_part.NodesBegin()->SolutionStepsDataHas(rAreaVar)) )
	     KRATOS_ERROR(std::logic_error,"Area Variable is not in the model part","")
	   if(!(rmodel_part.NodesBegin()->SolutionStepsDataHas(rPressureVar)) )
	     KRATOS_ERROR(std::logic_error,"rPressureVar Variable is not in the model part","")
	     
	   if(is_distributed == true)
	    if(!(rmodel_part.NodesBegin()->SolutionStepsDataHas(PARTITION_INDEX)) )
	     KRATOS_ERROR(std::logic_error,"PARTITION_INDEX Variable is not in the model part","")

	  //set as active the internal nodes
	  int node_size = rmodel_part.Nodes().size();
	  int elem_size = rmodel_part.Elements().size();
	  
	  #pragma omp parallel for firstprivate(node_size)
	  for(int i = 0; i<node_size; i++)
	  {
	      ModelPart::NodesContainerType::iterator it=rmodel_part.NodesBegin()+i;
	      it->FastGetSolutionStepValue(rAreaVar) = 0.0;
	      const double& dist = it->FastGetSolutionStepValue(rDistanceVar);
	      if(dist < 0.0)
		it->SetValue(IS_VISITED,1.0);
	      else
	      {
		it->SetValue(IS_VISITED,0.0);
		noalias(it->FastGetSolutionStepValue(rProjVar)) = ZeroVector(3);
		it->FastGetSolutionStepValue(rPressureVar) =  0.0;
	      }
	   }

	   array_1d<double,TDim+1> dist;
	   //now deselect all of the nodes of the divided elements
	   #pragma omp parallel for private(dist)
	   for(int i = 0; i<elem_size; i++)
	   {
	      PointerVector< Element>::iterator it=rmodel_part.ElementsBegin()+i;

	      Geometry<Node<3> >&geom = it->GetGeometry();
			      
	      for(unsigned int k=0; k<TDim+1; k++)
		dist[k] = geom[k].FastGetSolutionStepValue(rDistanceVar);
	      
	      if(IsDivided(dist))
	      {
		for(unsigned int k=0; k<TDim+1; k++)
		{
		    geom[k].SetLock();
		    geom[k].SetValue(IS_VISITED,0.0);
		    geom[k].UnSetLock();
		}
	      }
	   }

	   //now extrapolate the pressure projection var from all of the nodes that have at least one internal node_size
	   array_1d<double,TDim+1> visited;
	   array_1d<double,TDim+1> N;
	   array_1d<double,TDim> avg;
	   double lumping_factor = 1.0/double(TDim+1);
	   boost::numeric::ublas::bounded_matrix <double, TDim+1,TDim> DN_DX;
	   
	   for(unsigned int level=0; level<max_levels; level++)
	   {
		#pragma omp parallel for private(DN_DX,N,visited,avg) firstprivate(lumping_factor,elem_size)
		for(int i = 0; i<elem_size; i++)
		{
		    PointerVector< Element>::iterator it=rmodel_part.ElementsBegin()+i;

		    Geometry<Node<3> >&geom = it->GetGeometry();
				    
		    for(unsigned int k=0; k<TDim+1; k++)
			    visited[k] = geom[k].GetValue(IS_VISITED);
					  
		    if(IsActive(visited))
		    {	     
			double Volume;
			GeometryUtils::CalculateGeometryData(geom,DN_DX,N,Volume);
			
			//calculated average value of "velocity"
			noalias(avg) = ZeroVector(3);
			double counter = 0.0;
			for(unsigned int k=0; k<TDim+1; k++)
			{
			  if(visited[k]>0.0)
			  {
			      noalias(avg) += geom[k].FastGetSolutionStepValue(rProjVar);
			      counter += 1.0;
			  }
			}
			avg /= counter;
			
			for(unsigned int k=0; k<TDim+1; k++)
			{
			  if(visited[k]==0.0)
			  {
			      geom[k].SetLock();
			      geom[k].FastGetSolutionStepValue(rProjVar) += avg*Volume*lumping_factor;
			      geom[k].FastGetSolutionStepValue(rAreaVar) += Volume*lumping_factor;
			      geom[k].UnSetLock();
			  }
			}
		    }
		}	

		//mpi sync variables
		if(is_distributed == false)
		{
		  int my_rank = rmodel_part.GetCommunicator().MyPID();
		    #pragma omp parallel for firstprivate(node_size,my_rank)
		    for(int i = 0; i<node_size; i++)
		    {
			    ModelPart::NodesContainerType::iterator it=rmodel_part.NodesBegin()+i;
			    if(it->GetValue(IS_VISITED) == 1.0 && my_rank != it->FastGetSolutionStepValue(PARTITION_INDEX))
			    {
			      it->FastGetSolutionStepValue(rAreaVar) = 0.0;
			      noalias(it->FastGetSolutionStepValue(rProjVar)) = ZeroVector(3) ;
			    }
		    }
		    rmodel_part.GetCommunicator().AssembleCurrentData(rAreaVar);
		    rmodel_part.GetCommunicator().AssembleCurrentData(rProjVar);
		}
		
		//finish the computation of the nodal pressure gradient
		#pragma omp parallel for firstprivate(node_size)
		for(int i = 0; i<node_size; i++)
		{
			ModelPart::NodesContainerType::iterator it=rmodel_part.NodesBegin()+i;
			const double& area = it->FastGetSolutionStepValue(rAreaVar);
			if(area > 1e-20 && it->GetValue(IS_VISITED)==0.0) //this implies that node was computed
			{
			  it->FastGetSolutionStepValue(rProjVar) /= area;
			}
		}
		
		//now compute a nodal pressure which corresponds to the pressure gradient
		array_1d<double,3> xorigin,auxiliary_dist;
		#pragma omp parallel for private(DN_DX,N,visited,xorigin,avg,auxiliary_dist) firstprivate(lumping_factor,elem_size)
		for(int i = 0; i<elem_size; i++)
		{
		    PointerVector< Element>::iterator it=rmodel_part.ElementsBegin()+i;

		    Geometry<Node<3> >&geom = it->GetGeometry();
				    
		    for(unsigned int k=0; k<TDim+1; k++)
			    visited[k] = geom[k].GetValue(IS_VISITED);
					  
		    if(IsActive(visited))
		    {	     
			double Volume;
			GeometryUtils::CalculateGeometryData(geom,DN_DX,N,Volume);
			
			//calculated average value of "velocity"
			noalias(xorigin) = ZeroVector(3);
			noalias(avg) = ZeroVector(3);
			double porigin = 0.0;
			double counter = 0.0;
			for(unsigned int k=0; k<TDim+1; k++)
			{
			  if(visited[k]>0.0)
			  {
			      noalias(xorigin) += geom[k].Coordinates();
			      porigin += geom[k].FastGetSolutionStepValue(rPressureVar);
			      noalias(avg) += geom[k].FastGetSolutionStepValue(rProjVar);
			      counter += 1.0;
			  }
			}
			xorigin /= counter;
			avg /= counter;
			porigin /= counter;
			
			for(unsigned int k=0; k<TDim+1; k++)
			{
			  if(visited[k]==0.0)
			  {
			      noalias(auxiliary_dist) = geom[k].Coordinates();
			      noalias(auxiliary_dist) -= xorigin;
			      double deltap = inner_prod(auxiliary_dist,avg);
			      
			      geom[k].SetLock();
			      geom[k].FastGetSolutionStepValue(rPressureVar) += (porigin+deltap)*Volume*lumping_factor;
			      geom[k].UnSetLock();
			  }
			}
		    }
		}

		//mpi sync variables
		if(is_distributed == false)
		{
		  int my_rank = rmodel_part.GetCommunicator().MyPID();
		    #pragma omp parallel for firstprivate(node_size,my_rank)
		    for(int i = 0; i<node_size; i++)
		    {
			    ModelPart::NodesContainerType::iterator it=rmodel_part.NodesBegin()+i;
			    if(it->GetValue(IS_VISITED) == 1.0 && my_rank != it->FastGetSolutionStepValue(PARTITION_INDEX) )
			      it->FastGetSolutionStepValue(rPressureVar) = 0.0;
		    }
		    rmodel_part.GetCommunicator().AssembleCurrentData(rPressureVar);
		}
		
		//finish the computation of the nodal pressure
		#pragma omp parallel for firstprivate(node_size)
		for(int i = 0; i<node_size; i++)
		{
			ModelPart::NodesContainerType::iterator it=rmodel_part.NodesBegin()+i;
			const double& area = it->FastGetSolutionStepValue(rAreaVar);
			if(area > 1e-20 && it->GetValue(IS_VISITED)==0.0) //this implies that node was computed
			{
			  it->FastGetSolutionStepValue(rPressureVar) /= area;
			  it->SetValue(IS_VISITED,1.0); //node is marked as already visitedz
			}
		}
		
	   }	
	   
	   //assign to the body force the projection var for the nodes on which the PRESS_PROJ was not estimated otherwise
	  #pragma omp parallel for firstprivate(node_size)
	  for(int i = 0; i<node_size; i++)
	  {
		  ModelPart::NodesContainerType::iterator it=rmodel_part.NodesBegin()+i;
		  const double& dist = it->FastGetSolutionStepValue(rDistanceVar);
		  if(dist > 0.0 && it->GetValue(IS_VISITED)==0.0) //this implies that node was computed
		  {
		    it->FastGetSolutionStepValue(rProjVar) = it->FastGetSolutionStepValue(BODY_FORCE);
		  }
	  }
	         
	     KRATOS_CATCH("")
	}
	
      ///Function to assign a pressure to the first layer of nodes outside of the fluid.
      ///Pressure is assigned in such a way that the pressure gradient is mantained, but the pressure is 0 where the distance is 0
      ///the exact distance from the free surface is recalculated here
      ///@param rmodel_part is the ModelPart on which we will operate
      ///@param rDistanceVar is the Variable that we will use in calculating the distance
      ///@param rProjVar is the Variable being extrapolated
      ///@param rPressureVar is the Variable that we will use for L2 projections
      ///@param rAreaVar is the Variable that we will use for L2 projections
      void AssignFreeSurfacePressure(ModelPart& rmodel_part,
			     const Variable<double>& rDistanceVar, 
			     const Variable<array_1d<double,3> >& rProjVar, 
			     const Variable<double>& rPressureVar,
			     const Variable<double>& rAreaVar
			     )
	{
	  KRATOS_TRY
	  
	  bool is_distributed = false;
	  if(rmodel_part.GetCommunicator().TotalProcesses() > 1)
	    is_distributed = true;
	  
	   //check that variables needed are in the model part
	   if(!(rmodel_part.NodesBegin()->SolutionStepsDataHas(rDistanceVar)) )
	     KRATOS_ERROR(std::logic_error,"distance Variable is not in the model part","")
	   if(!(rmodel_part.NodesBegin()->SolutionStepsDataHas(rProjVar)) )
	     KRATOS_ERROR(std::logic_error,"Pressure Projection Variable is not in the model part","")
	   if(!(rmodel_part.NodesBegin()->SolutionStepsDataHas(rPressureVar)) )
	     KRATOS_ERROR(std::logic_error,"rPressureVar Variable is not in the model part","")
	   if(!(rmodel_part.NodesBegin()->SolutionStepsDataHas(rAreaVar)) )
	     KRATOS_ERROR(std::logic_error,"rAreaVar Variable is not in the model part","")

	     
	   if(is_distributed == true)
	    if(!(rmodel_part.NodesBegin()->SolutionStepsDataHas(PARTITION_INDEX)) )
	     KRATOS_ERROR(std::logic_error,"PARTITION_INDEX Variable is not in the model part","")

	  //set as active the internal nodes
	  int node_size = rmodel_part.Nodes().size();
	  int elem_size = rmodel_part.Elements().size();
	  
	  #pragma omp parallel for firstprivate(node_size)
	  for(int i = 0; i<node_size; i++)
	  {
	      ModelPart::NodesContainerType::iterator it=rmodel_part.NodesBegin()+i;
	      it->FastGetSolutionStepValue(rAreaVar) = 0.0;
	      const double& dist = it->FastGetSolutionStepValue(rDistanceVar);
	      if(dist < 0.0)
	      {
		it->SetValue(IS_VISITED,1.0);
	      }
	      else
	      {
		it->SetValue(IS_VISITED,0.0);
		it->FastGetSolutionStepValue(rPressureVar) =  0.0;
	      }
	   }
	  
	   array_1d<double,TDim+1> dist, exact_dist;
	   array_1d<double,TDim> grad_d;
	   array_1d<double,TDim+1> N;
	   double lumping_factor = 1.0/double(TDim+1);
	   boost::numeric::ublas::bounded_matrix <double, TDim+1,TDim> DN_DX;
	   
	   //now deselect all of the nodes of the divided elements
	   #pragma omp parallel for private(dist,N,DN_DX)
	   for(int i = 0; i<elem_size; i++)
	   {
	      PointerVector< Element>::iterator it=rmodel_part.ElementsBegin()+i;

	      Geometry<Node<3> >&geom = it->GetGeometry();
			      
	      for(unsigned int k=0; k<TDim+1; k++)
		dist[k] = geom[k].FastGetSolutionStepValue(rDistanceVar);
	      
	      if(IsDivided(dist))
	      {
		    //calculate exact distance
		    double Volume;
		    GeometryUtils::CalculateGeometryData(geom,DN_DX,N,Volume);
		    
		    ComputeExactDistances(DN_DX,Volume,geom,dist,exact_dist);
		    
		    //compute the direction of the gradient of the distances
		    for(unsigned int iii=0; iii<TDim; iii++)
			    grad_d[iii] = DN_DX(0,iii)*dist[0];
		    for(unsigned int k=1; k<TDim+1; k++)
			for(unsigned int iii=0; iii<TDim; iii++)
			    grad_d[iii] += DN_DX(k,iii)*dist[k];
		    double norm_grad = norm_2(grad_d);
		    grad_d/= norm_grad;
		    
		    //assign pressure contribution
		    for(unsigned int k=0; k<TDim+1; k++)
		    {
		      if(dist[k]>=0.0)
		      {
			  const array_1d<double,3>& press_proj = geom[k].FastGetSolutionStepValue(rProjVar);
			  
			  double pouter = 0.0;
			  for(unsigned int iii=0; iii<TDim; iii++)
			    pouter += grad_d[iii]*press_proj[iii];
			  pouter *= exact_dist[k];

			  geom[k].SetLock();
			  geom[k].FastGetSolutionStepValue(rPressureVar) += pouter*Volume*lumping_factor;
			  geom[k].FastGetSolutionStepValue(rAreaVar) += Volume*lumping_factor;
			  geom[k].UnSetLock();
		      }
		    }
		
	      }
	   }
	   
	   rmodel_part.GetCommunicator().AssembleCurrentData(rAreaVar);
	   rmodel_part.GetCommunicator().AssembleCurrentData(rPressureVar);

	  //finish the computation of the nodal pressure gradient
	  #pragma omp parallel for firstprivate(node_size)
	  for(int i = 0; i<node_size; i++)
	  {
		  ModelPart::NodesContainerType::iterator it=rmodel_part.NodesBegin()+i;
		  const double& area = it->FastGetSolutionStepValue(rAreaVar);
		  const double& dist = it->FastGetSolutionStepValue(rDistanceVar);
		  if(area > 1e-20 && it->GetValue(IS_VISITED)==0.0) //this implies that node was computed
		  {
		    it->SetValue(IS_VISITED,1.0);
		    it->FastGetSolutionStepValue(rPressureVar) /= area;
		  }
		  else if(dist >= 0.0)
		    it->FastGetSolutionStepValue(rPressureVar) = -0.000001;
	  }
        
	KRATOS_CATCH("")
	}	
      ///@}
      ///@name Operators 
      ///@{
      
      
      ///@}
      ///@name Operations
      ///@{
      
      
      ///@}
      ///@name Access
      ///@{ 
      
      
      ///@}
      ///@name Inquiry
      ///@{
      
      
      ///@}      
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
      {
	std::stringstream buffer;
        buffer << "ParallelExtrapolationUtilities" << TDim << "D";
        return buffer.str();
      };
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "ParallelExtrapolationUtilities" << TDim << "D";};

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const {};
      
            
      ///@}      
      ///@name Friends
      ///@{
      
            
      ///@}
      
    protected:
      ///@name Protected static Member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Protected member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Protected Operators
      ///@{ 
      
      //*******************************************************************
      bool IsDivided(array_1d<double,TDim+1>& dist)
      {
	  unsigned int positive = 0;
	  unsigned int negative = 0;
	  
	  for(unsigned int i=0; i<TDim+1; i++)
	  {
	    if(dist[i] >= 0)
	      positive++;
	    else
	      negative++;
	  }
	    
	   bool is_divided = false;
	   if(positive > 0 && negative>0)
	     is_divided = true;
	   
	   return is_divided;
      }
      
      //*******************************************************************
      bool IsActive(array_1d<double,TDim+1>& visited)
      {
	  unsigned int positive = 0;
	  
	  for(unsigned int i=0; i<TDim+1; i++)
	    if(visited[i] > 0.9999999999) //node was considered
	      positive++;
	    
	   bool is_active = false;
	   if(positive >= 1)
	     is_active = true;
	   
	   return is_active;
      }
      
      //*******************************************************************
      void ComputeExactDistances(const boost::numeric::ublas::bounded_matrix <double, TDim+1,TDim>& DN_DX,
				 const double& Area,
				 Geometry<Node<3> >& geom,
				 const array_1d<double,TDim+1>& distances,
				 array_1d<double,TDim+1>& exact_dist
				 )
      {
	array_1d<double,TDim> grad_d;
	array_1d<double,3> coord_on_0;
	coord_on_0[0] = -100000000000.0;
	coord_on_0[1] = -100000000000.0;
	coord_on_0[2] = -100000000000.0;
	array_1d<double,3> temp;
	
	//compute the gradient of the distance and normalize it
	noalias(grad_d) = prod(trans(DN_DX),distances);
	double norm = norm_2(grad_d);
	grad_d /= norm;
				    
	//find one division point on one edge
	for(unsigned int i = 1; i<TDim+1; i++)
	{
	    if(distances[0]*distances[i]<=0.0) //if the edge is divided
	      {
		    double delta_d = fabs(distances[i]) + fabs(distances[0]);
		    
		    if(delta_d>1e-20)
		    {
			double Ni = fabs(distances[0]) / delta_d;
			double N0 = fabs(distances[i]) / delta_d;

			noalias(coord_on_0) = N0 * geom[0].Coordinates();
			noalias(coord_on_0) += Ni * geom[i].Coordinates();
		    }
		    else
		        noalias(coord_on_0) = geom[0].Coordinates();

		    break;

	    }
	}
	
	//double nodal_area_factor = Area/static_cast<double>(TDim+1);

	//now calculate the distance of all the nodes from the elemental free surface
	for(unsigned int i = 0; i<TDim+1; i++)
	{
	    noalias(temp) = geom[i].Coordinates();
	    noalias(temp) -= coord_on_0 ;

	    double real_distance = 0.0;
	    for(unsigned int k=0; k<TDim; k++)
		real_distance += temp[k]*grad_d[k];
	    real_distance = fabs(real_distance);

	    exact_dist[i] = real_distance;
	}	
      }      
      ///@} 
      
      ///@name Protected Operations
      ///@{ 
        
        
      ///@} 
      ///@name Protected  Access 
      ///@{ 
        
        
      ///@}      
      ///@name Protected Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Protected LifeCycle 
      ///@{ 
      
            
      ///@}
      
    private:
      ///@name Static Member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Member Variables 
      ///@{ 
        
        
      ///@} 
      ///@name Private Operators
      ///@{ 
        
        
      ///@} 
      ///@name Private Operations
      ///@{ 
        
        
      ///@} 
      ///@name Private  Access 
      ///@{ 
        
        
      ///@}    
      ///@name Private Inquiry 
      ///@{ 
        
        
      ///@}    
      ///@name Un accessible methods 
      ///@{ 
      
      /// Assignment operator.
       ParallelExtrapolationUtilities<TDim>& operator=(ParallelExtrapolationUtilities<TDim> const& rOther){};

      /// Copy constructor.
       ParallelExtrapolationUtilities(ParallelExtrapolationUtilities<TDim> const& rOther){};

        
      ///@}    
        
    }; // Class ParallelExtrapolationUtilities 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template<unsigned int TDim>
  inline std::istream& operator >> (std::istream& rIStream, 
				    ParallelExtrapolationUtilities<TDim>& rThis){return rIStream;}

  /// output stream function
   template<unsigned int TDim>
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const ParallelExtrapolationUtilities<TDim>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_PARALLEL_EXTRAPOLATION_UTILITIES_H_INCLUDED  defined 


