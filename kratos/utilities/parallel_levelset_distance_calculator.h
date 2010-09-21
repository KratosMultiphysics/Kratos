//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_PARALLEL_DISTANCE_CALCULATOR_H_INCLUDED )
#define  KRATOS_PARALLEL_DISTANCE_CALCULATOR_H_INCLUDED



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
  class ParallelDistanceCalculator
    {
    public:
      ///@name Type Definitions
      ///@{
      
      /// Pointer definition of ParallelDistanceCalculator
      KRATOS_CLASS_POINTER_DEFINITION(ParallelDistanceCalculator);
  
      ///@}
      ///@name Life Cycle 
      ///@{ 
      
      /// Default constructor.
      ParallelDistanceCalculator(){};

      /// Destructor.
      virtual ~ParallelDistanceCalculator(){};
      
      ///Function to calculate a signed distance function suitable for calculations using the Level Set Method
      ///the function assumes given a "signed distance" distributions and recomputes the distances
      ///respecting as accurately as possible the position of the zero of the original distributions
      ///@param rmodel_part is the ModelPart on which we will operate
      ///@param rDistanceVar is the Variable that we will use in calculating the distance
      ///@param rAreaVar is the Variable that we will use for L2 projections
      ///@param max_levels is the number of maximum "layers" of element that will be used in the calculation of the distances
      ///@param max_distance distances will not be computed after reaching this limit
      void CalculateDistances(ModelPart& rmodel_part,
			     const Variable<double>& rDistanceVar, 
			     const Variable<double>& rAreaVar, 
			     const unsigned int max_levels,
			     const double max_distance)
	{
	  KRATOS_TRY
	   //check that variables needed are in the model part
	   if(!rmodel_part.Has(rDistanceVar))
	     KRATOS_ERROR(std::logic_error,"distance Variable is not in the model part","")
	   if(!rmodel_part.Has(rAreaVar))
	     KRATOS_ERROR(std::logic_error,"Area Variable is not in the model part","")
	     
	   //reset the variables needed
	   for(ModelPart::NodesContainerType::iterator it=rmodel_part.NodesBegin(); it!=rmodel_part.NodesEnd(); it++)
	   {
	      it->FastGetSolutionStepValue(rAreaVar) = 0.0;
	      double& dist = it->FastGetSolutionStepValue(rDistanceVar);
	      it->SetValue(rDistanceVar,dist); //here we copy the distance function to the fixed database
	      dist = 0.0; 
	   }
	   
	   //*****************************************************************+
	   //*****************************************************************+
	   //*****************************************************************+
	   //identify the list of elements divided by the original distance distribution and recompute an "exact" distance 
	   //attempting to mantain the original position of the free surface
	   //note that the backup value is used in calculating the position of the free surface and the divided elements
	   array_1d<double,TDim+1> dist, exact_dist;
	   array_1d<double,TDim+1> areas;
	   array_1d<double,TDim+1> N;
	   double lumping_factor = 1.0/double(TDim+1);
	   boost::numeric::ublas::bounded_matrix <double, TDim+1,TDim> DN_DX;
	   for(ModelPart::ElementsContainerType::iterator it=rmodel_part.ElementsBegin(); it!=rmodel_part.ElementsEnd(); it++)
	   {
		  Geometry<Node<3> >&geom = it->GetGeometry();
		  
		  
		  
		  for(unsigned int i=0; i<TDim+1; i++)
		    dist[i] = geom[i].GetValue(rDistanceVar);
		  
		  bool is_divided = IsDivided(dist);
		  
		  if(is_divided == true)
		  {
		    double Volume;
                    GeometryUtils::CalculateGeometryData(geom,DN_DX,N,Volume);
		    
		    ComputeExactDistances(DN_DX,Volume,geom,dist,exact_dist);
		    
		    for(unsigned int i=0; i<TDim+1; i++)
		    {
		      geom[i].SetLock();
		      geom[i].FastGetSolutionStepValue(rDistanceVar) += exact_dist[i];
		      geom[i].FastGetSolutionStepValue(rAreaVar) += Volume*lumping_factor;
		      geom[i].UnSetLock();
		    }
		  }
	   }	
	   
	   //mpi sync variables
	   rmodel_part.GetCommunicator().AssembleCurrentData(rAreaVar);
	   rmodel_part.GetCommunicator().AssembleCurrentData(rDistanceVar);
	   
	   for(ModelPart::NodesContainerType::iterator it=rmodel_part.NodesBegin(); it!=rmodel_part.NodesEnd(); it++)
	   {
		  double area = it->FastGetSolutionStepValue(rAreaVar);
		  if(area > 1e-20) //this implies that node was computed
		  {
		    double& dist = it->FastGetSolutionStepValue(rDistanceVar);
		    dist /= area;
		  }
	   }
	   
	   
	   	 
	   //*****************************************************************+
	   //*****************************************************************+
	   //*****************************************************************+
	   PointerVector< Element > ActiveElements;
	   ActiveElements.reserve( rmodel_part.Nodes().size() /4 );
	   
	   //now extend the distances layer by layer up to a maximum level of layers
	   for(unsigned int level=0; level<max_levels; level++)
	   {
		//erase the active elements array
		ActiveElements.clear();
		
		//fill the list of active elements
		for(ModelPart::ElementsContainerType::iterator it=rmodel_part.ElementsBegin(); it!=rmodel_part.ElementsEnd(); it++)
		{
		    Geometry<Node<3> >&geom = it->GetGeometry();
		    
		    for(unsigned int i=0; i<TDim+1; i++)
		      areas[i] = geom[i].FastGetSolutionStepValue(rAreaVar);
				    
		    if(IsActive(areas))
		      ActiveElements.push_back(*(it.base()));
		}
		
		//loop on active elements and advance the distance computation
		for(PointerVector< Element>::iterator it=ActiveElements.begin();it!= ActiveElements.end(); it++)
		{
		     Geometry<Node<3> >&geom = it->GetGeometry();
		     double Volume;
		     GeometryUtils::CalculateGeometryData(geom,DN_DX,N,Volume);
		     
		     AddDistanceToNodes(rDistanceVar,rAreaVar,geom,DN_DX,Volume);
		}
		
		//mpi sync variables
		rmodel_part.GetCommunicator().AssembleCurrentData(rAreaVar);
		rmodel_part.GetCommunicator().AssembleCurrentData(rDistanceVar);
		
		//finalize the computation of the distance
		for(ModelPart::NodesContainerType::iterator it=rmodel_part.NodesBegin(); it!=rmodel_part.NodesEnd(); it++)
		{
			double area = it->FastGetSolutionStepValue(rAreaVar);
			if(area > 1e-20) //this implies that node was computed
			{
			  double& dist = it->FastGetSolutionStepValue(rDistanceVar);
			  dist /= area;
			}
		}
	   }
	   
	   	 
	   //*****************************************************************+
	   //*****************************************************************+
	   //*****************************************************************+
	   //assign the sign to the distance function according to the original distribution. Set to max for nodes that were not calculated
	  for(ModelPart::NodesContainerType::iterator it=rmodel_part.NodesBegin(); it!=rmodel_part.NodesEnd(); it++)
	  {
		  const double area = it->FastGetSolutionStepValue(rAreaVar);
		  const double old_dist = it->GetValue(rDistanceVar);
		  double& dist = it->FastGetSolutionStepValue(rDistanceVar);
		  
		  if(dist > max_distance)
		    dist = max_distance;
		  
		  if(old_dist < 0)
		      dist = -max_distance;
	  }	   
	     
	     KRATOS_CATCH("")
	}
	
	
	
	
	//**********************************************************************************
	//**********************************************************************************
        double FindMaximumEdgeSize(ModelPart& r_model_part) 
        {
            KRATOS_TRY

            ModelPart::NodesContainerType& rNodes = r_model_part.Nodes();

            double h_max = 0.0;
	    
	    for(ModelPart::ElementsContainerType::iterator it=r_model_part.ElementsBegin(); it!=r_model_part.ElementsEnd(); it++)
	    {
		Geometry<Node<3> >&geom = it->GetGeometry();
		
		double h = 0.0; 
		
		for(unsigned int i=0; i<TDim+1; i++)
		{
		  
		  double xc = geom[i].X();
		  double yc = geom[i].Y();
		  double zc = geom[i].Z();
		  for(unsigned int j=i+1; j<TDim+1; j++)
		  {
		    double x = geom[j].X();
                    double y = geom[j].Y();
                    double z = geom[j].Z();
                    double l = (x - xc)*(x - xc);
                    l += (y - yc)*(y - yc);
                    l += (z - zc)*(z - zc);

                    if (l > h) h = l;
		  }
		}
		
		h = sqrt(h);

                if(h > h_max) h_max = h;
		  
	    }
	    
// 	    rmodel_part.GetCommunicator().AllReduce(h_max,KRATOS_REDUCE_MIN);
	    
            return h_max;

            KRATOS_CATCH("");
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
        buffer << "ParallelDistanceCalculator" << TDim << "D";
        return buffer.str();
      };
      
      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "ParallelDistanceCalculator" << TDim << "D";};

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
	    else;
	      negative++;
	  }
	    
	   bool is_divided = false;
	   if(positive > 0 and negative>0)
	     is_divided = true;
	   
	   return is_divided;
      }
      
      //*******************************************************************
      bool IsActive(array_1d<double,TDim+1>& areas)
      {
	  unsigned int positive = 0;
	  
	  for(unsigned int i=0; i<TDim+1; i++)
	    if(areas[i] > 1e-20) //node was considered
	      positive++;
	    
	   bool is_active = false;
	   if(positive > 0 )
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

		    double Ni = fabs(distances[0]) / delta_d;
		    double N0 = fabs(distances[i]) / delta_d;

		    noalias(coord_on_0) = N0 * geom[0].Coordinates();
		    noalias(coord_on_0) += Ni * geom[i].Coordinates();

		    break;

	    }
	}
	
	double nodal_area_factor = Area/static_cast<double>(TDim+1);

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
        
      //*******************************************************************
      void AddDistanceToNodes(const Variable<double>& rDistanceVar,
			      const Variable<double>& rAreaVar,
			      Geometry<Node<3> >& geom,
			      const boost::numeric::ublas::bounded_matrix <double, TDim+1,TDim>& DN_DX,
			      const double& Volume
			      )
      {
	      unsigned int unknown_node_index = 0;
	      array_1d<double,TDim> d;
	      double nodal_vol = Volume/static_cast<double>(TDim+1);
	      

	      //compute discriminant
	      noalias(d) = ZeroVector(TDim);
	      for (unsigned int iii = 0; iii < TDim + 1; iii++)
	      {
		  const double distance = geom[iii].FastGetSolutionStepValue(rDistanceVar);
		  double nodal_area = geom[iii].FastGetSolutionStepValue(rAreaVar);
		  
		  if (nodal_area > 1e-20) //identyfing the unknown node
		      for (unsigned int jjj = 0; jjj < TDim; jjj++)
			  d[jjj] += DN_DX(iii, jjj) * distance;
		  else
		      unknown_node_index = iii;
	      }

	      //finalizing computation of discriminant
	      double c = -1.0;
	      double a = 0.0;
	      double b = 0.0;
	      for (unsigned int jjj = 0; jjj < TDim; jjj++)
	      {
		  a += DN_DX(unknown_node_index, jjj) * DN_DX(unknown_node_index, jjj);
		  b += d[jjj] * DN_DX(unknown_node_index, jjj);
		  c += d[jjj] * d[jjj];
	      }
	      b *= 2.0;

	      //here we require (a*x^2 + b*x + c)^2 to be minimum (x represents the unknown distance)
	      //this implies setting to zero
	      //(a*x^2 + b*x + c)*(2ax+b) = 0
	      double distance = -10000.0;
	      
	      double discriminant = b * b - 4.0 * a*c;

	      if (discriminant < 0.0) //here we solve (2ax+b) = 0
	      {
		  distance = -b / (2.0*a);
	      } 
	      else //in this case we solve (a*x^2 + b*x + c)=0
	      {
		  //(accurate) computation of the distance
		  //requires the solution of a*x^2+b*x+c=0
		  double q, root1, root2;
		  if (a != 0.0)
		  {
		      if (b > 0) q = -0.5 * (b + sqrt(discriminant));
		      else q = -0.5 * (b - sqrt(discriminant));
		      root1 = q / a;
		      root2 = c / q;
		      if (root1 > root2) distance = root1;
		      else distance = root2;
		  } else //in this case we have a linear equation
		  {
		      distance = -c / b;
		  }
	      }
	      
	      geom[unknown_node_index].SetLock();
	      geom[unknown_node_index].FastGetSolutionStepValue(rDistanceVar) += distance;
	      geom[unknown_node_index].FastGetSolutionStepValue(rAreaVar) += nodal_vol;
	      geom[unknown_node_index].UnSetLock();	
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
       ParallelDistanceCalculator<TDim>& operator=(ParallelDistanceCalculator<TDim> const& rOther){};

      /// Copy constructor.
       ParallelDistanceCalculator(ParallelDistanceCalculator<TDim> const& rOther){};

        
      ///@}    
        
    }; // Class ParallelDistanceCalculator 

  ///@} 
  
  ///@name Type Definitions       
  ///@{ 
  
  
  ///@} 
  ///@name Input and output 
  ///@{ 
        
 
  /// input stream function
  template<unsigned int TDim>
  inline std::istream& operator >> (std::istream& rIStream, 
				    ParallelDistanceCalculator<TDim>& rThis){}

  /// output stream function
   template<unsigned int TDim>
  inline std::ostream& operator << (std::ostream& rOStream, 
				    const ParallelDistanceCalculator<TDim>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@} 
  
  
}  // namespace Kratos.

#endif // KRATOS_PARALLEL_DISTANCE_CALCULATOR_H_INCLUDED  defined 


