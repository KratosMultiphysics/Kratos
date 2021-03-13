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
//   Last Modified by:    $Author: dcagri $
//   Date:                $Date: 2021-03-10 15:19:00 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_INTEGRATE_UTILITY_INCLUDED)
#define  KRATOS_INTEGRATE_UTILITY_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <algorithm>

// External includes 


// Project includes
#include "includes/define.h"
#include "includes/node.h"

///
#include "includes/dof.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"
#include "includes/deprecated_variables.h"
#include "containers/array_1d.h"
#include "containers/data_value_container.h"
#include "includes/mesh.h"
#include "utilities/math_utils.h"
#include "processes/node_erase_process.h" 
///

#include "utilities/geometry_utilities.h"

#include "includes/model_part.h"




#include "spatial_containers/spatial_containers.h"
// #include "spatial_containers/bounding_box.h"
#include "spatial_containers/cell.h"
#include "spatial_containers/bins_dynamic_objects.h"

#include "utilities/spatial_containers_configure.h"

#include "geometries/line_2d_2.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/point.h"

#include "pfem_2_application.h"
#include "pfem_particle_fluidonly.h"

//#include "utilities/enrich_2d_2dofs.h"
#include "utilities/enrichment_utilities.h"
#include "utilities/divide_triangle_2d_3.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"


#include "includes/element.h"


#include "includes/serializer.h"



#include "utilities/openmp_utils.h"

#include "time.h"

//#include "processes/process.h"


namespace Kratos
{
	template< unsigned int TDim>
	class IntegrateTractionVectorUtility
	{
	public:
	
	    typedef SpatialContainersConfigure<TDim>     Configure;   
	    typedef typename Configure::PointType                      PointType; 
	    //typedef PointType::CoordinatesArrayType           CoordinatesArrayType;
        typedef typename Configure::ContainerType                  ContainerType;   
        //typedef Configure::PointerType                    PointerType;
        typedef typename Configure::IteratorType                   IteratorType; 
        typedef typename Configure::ResultContainerType            ResultContainerType;
	    //typedef Configure::ResultPointerType              ResultPointerType;
        typedef typename Configure::ResultIteratorType             ResultIteratorType; 
        typedef PointerVector< PFEM_Particle_Fluid, PFEM_Particle_Fluid*, std::vector<PFEM_Particle_Fluid*> > ParticlePointerVector;

		KRATOS_CLASS_POINTER_DEFINITION(IntegrateTractionVectorUtility);
        

        struct element_data
        {
          bounded_matrix<double,TDim+1, TDim> v, vn, vnn, f; 
          array_1d<double,TDim+1> p, rho;
          
          bounded_matrix<double, TDim+1, TDim > DN_DX;
          array_1d<double, TDim+1 > N;
          
          Matrix C;
          Vector stress;
          
          double bdf0;
          double bdf1;
          double bdf2;
          double h;
          double dyn_tau_coeff;
        };

		//template<unsigned int TDim>
		TractionVectorIntegrationUtility(ModelPart& model_part)
			: mr_model_part(model_part)
		{
			std::cout << "initializing IntegrateTractionVectorUtility" << std::endl;
			
			Check();
		
		}
		
         
		~IntegrateTractionVectorUtility()
		{}
		
		vector<double> ReturnDrag()
        {
            KRATOS_TRY
            
            const SizeType NumNodes = TDim+1;
            element_data data; 
            double drag=0.0;
            double lift=0.0;
            double Volume = 0.0;
            double Area=0.0;
            array_1d<double,TDim> traction_vector = ZeroVector(TDim); 
            array_1d<double,9> traction_vector_aux = ZeroVector(9); 
             
                  
            
            auto &rConditionsArray = mr_model_part.Conditions();
            const auto it_condition_begin = rConditionsArray.begin();
            array_1d<double, 3> dummy;
            for (int i = 0; i < static_cast<int>(rConditionsArray.size()); ++i) 
            {
             auto it_condition = it_condition_begin + i;
             boost::numeric::ublas::bounded_matrix<double,TDim*(TDim+1),TDim> Kvu_nitsche_Cmatrix = ZeroMatrix(TDim*(TDim+1),TDim);
             boost::numeric::ublas::bounded_matrix<double,TDim,TDim*(TDim+1)> Kvu_nitsche_Dmatrix = ZeroMatrix(TDim,TDim*(TDim+1)); 
             boost::numeric::ublas::bounded_matrix<double,TDim*(TDim+1),TDim*(TDim+1)> Kvu_nitsche = ZeroMatrix(TDim*(TDim+1),TDim*(TDim+1)); 
             boost::numeric::ublas::bounded_matrix<double,TDim*(TDim+1),TDim+1> Kvp_nitsche = ZeroMatrix(TDim*(TDim+1),TDim+1);
             boost::numeric::ublas::bounded_matrix<double,2,TDim+1> Kvp_nitsche_Dmatrix = ZeroMatrix(2,TDim+1);
             boost::numeric::ublas::bounded_matrix<double,TDim*(TDim+1)*3/2,TDim*(TDim+1)*3/2> Kvu_nitsche_large = ZeroMatrix(TDim*(TDim+1)*3/2,TDim*(TDim+1)*3/2); 
             boost::numeric::ublas::bounded_matrix<double,TDim*(TDim+1)*3/2,TDim*(TDim+1)*3/2> Kvp_nitsche_large = ZeroMatrix(TDim*(TDim+1)*3/2,TDim*(TDim+1)*3/2); 
             
             array_1d<double,(TDim*2)> gp = ZeroVector(TDim*2);
             array_1d<double,(TDim+1)> N = ZeroVector(TDim+1);

             array_1d<double,TDim> normal_vector = ZeroVector(TDim); 
             GlobalPointersVector< Element >& neighbor_els = it_condition->GetValue(NEIGHBOUR_ELEMENTS);
             if (it_condition->GetValue(NEIGHBOUR_ELEMENTS).size()>1)
               KRATOS_THROW_ERROR(std::invalid_argument,"NEIGHBOUR_ELEMENTS is more than one","");
             if (it_condition->GetValue(NEIGHBOUR_ELEMENTS).size()<1)
               KRATOS_THROW_ERROR(std::invalid_argument,"NEIGHBOUR_ELEMENTS does not exist","");
             for(GlobalPointersVector< Element >::iterator ielem = neighbor_els.begin(); ielem!=neighbor_els.end(); ielem++)
             {
              Geometry< Node<3> >& geom = ielem->GetGeometry();
              GeometryUtils::CalculateGeometryData(geom, data.DN_DX, data.N, Volume);
              const bounded_matrix<double,NumNodes,TDim>& DN = data.DN_DX;
              const double mu = ielem->GetProperties()[VISCOSITY];
              array_1d<double,TDim*(TDim+1)> v_aux = ZeroVector(TDim*(TDim+1));
              array_1d<double,(TDim+1)> p_aux = ZeroVector(TDim+1);
              array_1d<double,(TDim+1)*(TDim+1)> temp_vector;
    
              v_aux[0]=(ielem->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X));
              v_aux[1]=(ielem->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y));
              v_aux[2]=(ielem->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_X));
              v_aux[3]=(ielem->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY_Y));
              v_aux[4]=(ielem->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_X));
              v_aux[5]=(ielem->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY_Y));
              p_aux[0]=(ielem->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE));
              p_aux[1]=(ielem->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE));
              p_aux[2]=(ielem->GetGeometry()[2].FastGetSolutionStepValue(PRESSURE));
              int counter = 0;
              for(unsigned int iii = 0; iii<(TDim+1); iii++)
              {
               temp_vector[counter++] = v_aux[TDim*iii];
               temp_vector[counter++] = v_aux[TDim*iii+1];
               temp_vector[counter++] = p_aux[iii];             
              }

             normal_vector = it_condition->GetGeometry().Normal(dummy);
             double normal_vector_magnitude=sqrt(normal_vector[0]*normal_vector[0]+normal_vector[1]*normal_vector[1]);
             normal_vector[0]=normal_vector[0]/normal_vector_magnitude;
             normal_vector[1]=normal_vector[1]/normal_vector_magnitude;
             
             double  length =it_condition->GetGeometry().Length();
             double xdistance = it_condition->GetGeometry()[1].X() - it_condition->GetGeometry()[0].X();
             double ydistance = it_condition->GetGeometry()[1].Y() - it_condition->GetGeometry()[0].Y();
             gp[0]=it_condition->GetGeometry()[0].X()+0.211324865*xdistance;
             gp[1]=it_condition->GetGeometry()[0].Y()+0.211324865*ydistance;
             gp[2]=it_condition->GetGeometry()[0].X()+0.788675135*xdistance;
             gp[3]=it_condition->GetGeometry()[0].Y()+0.788675135*ydistance;
             

             
             double TotalArea=0.0;
             double x10 = geom[1].X() - geom[0].X();
             double y10 = geom[1].Y() - geom[0].Y();

             double x20 = geom[2].X() - geom[0].X();
             double y20 = geom[2].Y() - geom[0].Y();
             double detJ = x10 * y20-y10 * x20;
             TotalArea = 0.5*detJ;
             
             //integration over two gauss points starts here
             for(int m = 0; m<2; m++)
             {
              double Area0=0.0;
              x10 = geom[1].X() - gp[m*2];
              y10 = geom[1].Y() - gp[m*2+1];

              x20 = geom[2].X() - gp[m*2];
              y20 = geom[2].Y() - gp[m*2+1];
              detJ = x10 * y20-y10 * x20;
              Area0 = 0.5*detJ;
              N[0]=Area0/TotalArea;
              
              double Area1=0.0;
              x10 = gp[m*2] - geom[0].X();
              y10 = gp[m*2+1] - geom[0].Y();

              x20 = geom[2].X() - geom[0].X();
              y20 = geom[2].Y() - geom[0].Y();
              detJ = x10 * y20-y10 * x20;
              Area1 = 0.5*detJ;
              N[1]=Area1/TotalArea;
              
              double Area2=0.0;
              x10 = geom[1].X() - geom[0].X();
              y10 = geom[1].Y() - geom[0].Y();

              x20 = gp[m*2] - geom[0].X();
              y20 = gp[m*2+1] - geom[0].Y();
              detJ = x10 * y20-y10 * x20;
              Area2 = 0.5*detJ;
              N[2]=Area2/TotalArea;
              
              Kvu_nitsche_large = ZeroMatrix(TDim*(TDim+1)*3/2,TDim*(TDim+1)*3/2); 
              Kvp_nitsche_large = ZeroMatrix(TDim*(TDim+1)*3/2,TDim*(TDim+1)*3/2); 
              
              for (unsigned int i = 0; i < (TDim+1); i++)
              {
                Kvu_nitsche_Cmatrix(i*TDim,0)  =N[i]; 
                Kvu_nitsche_Cmatrix(i*TDim+1,1)=N[i]; 
                
                Kvu_nitsche_Dmatrix(0,i*TDim)  =2.0*mu*(DN(i,0)*normal_vector(0)+0.5*DN(i,1)*normal_vector(1));
                Kvu_nitsche_Dmatrix(0,i*TDim+1)=2.0*mu*(0.5*DN(i,0)*normal_vector(1));
                
                Kvu_nitsche_Dmatrix(1,i*TDim)  =2.0*mu*(0.5*DN(i,1)*normal_vector(0));
                Kvu_nitsche_Dmatrix(1,i*TDim+1)=2.0*mu*(DN(i,1)*normal_vector(1)+0.5*DN(i,0)*normal_vector(0));
              }
     
             Kvu_nitsche=prod(Kvu_nitsche_Cmatrix,Kvu_nitsche_Dmatrix);
             
             
             
             for (int i=0; i<(TDim+1); i++)
             {
               for (int j=0; j<(TDim+1); j++)
               {
                 for (int k=0; k<TDim; k++)
                 {
                   for (int l=0; l<TDim; l++)
                   {
                     Kvu_nitsche_large(i * (TDim+1) + k, j * (TDim+1) + l) += Kvu_nitsche(i * (TDim) + k, j * (TDim) + l);
                   }
                 }
               }
             }
             
             Kvp_nitsche_Dmatrix(0,0)=N[0]*normal_vector(0);
             Kvp_nitsche_Dmatrix(0,1)=N[1]*normal_vector(0);
             Kvp_nitsche_Dmatrix(0,2)=N[2]*normal_vector(0);
             Kvp_nitsche_Dmatrix(1,0)=N[0]*normal_vector(1);
             Kvp_nitsche_Dmatrix(1,1)=N[1]*normal_vector(1);
             Kvp_nitsche_Dmatrix(1,2)=N[2]*normal_vector(1);
             
             
             
             Kvp_nitsche=prod(Kvu_nitsche_Cmatrix,Kvp_nitsche_Dmatrix);
             
             for (int i=0; i<(TDim+1); i++)
             {
               for (int j=0; j<(TDim+1); j++)
               {
                 for (int k=0; k<TDim; k++)
                 {
                   Kvp_nitsche_large(i * (TDim+1) + k, j * (TDim+1) + TDim) -= Kvp_nitsche(i * TDim  + k, j);              
                 }
               }
             }
              
              
            
             traction_vector_aux=0.5*length*(prod(Kvu_nitsche_large,temp_vector)+prod(Kvp_nitsche_large,temp_vector));
             
               
             drag+=traction_vector_aux[0]+traction_vector_aux[3]+traction_vector_aux[6];
             lift+=traction_vector_aux[1]+traction_vector_aux[4]+traction_vector_aux[7];
               
             }
             
             

             }
            }
            
            traction_vector[0]=drag;
            traction_vector[1]=lift;
            
            return traction_vector;

            KRATOS_CATCH("")
        }
		

		

		
	protected:


	private:
	
	
	void Check()
	{
// 		if(mr_model_part.NodesBegin()->SolutionStepsDataHas(DISTANCE) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing DISTANCE variable on solution step data","");
// 		if(mr_model_part.NodesBegin()->SolutionStepsDataHas(PRESS_PROJ) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESS_PROJ variable on solution step data","");
		if(mr_model_part.NodesBegin()->SolutionStepsDataHas(VELOCITY) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data","");
		if(mr_model_part.NodesBegin()->SolutionStepsDataHas(PRESSURE) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE variable on solution step data","");
// 		if(mr_model_part.NodesBegin()->SolutionStepsDataHas(PROJECTED_VELOCITY) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing PROJECTED_VELOCITY variable on solution step data","");
// 		if(mr_model_part.NodesBegin()->SolutionStepsDataHas(DELTA_VELOCITY) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing DELTA_VELOCITY variable on solution step data","");
// 		if(mr_model_part.NodesBegin()->SolutionStepsDataHas(MESH_VELOCITY) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing MESH_VELOCITY variable on solution step data","");
// 		if(mr_model_part.NodesBegin()->SolutionStepsDataHas(YP) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing YP variable on solution step data","");
// 		if(mr_model_part.NodesBegin()->SolutionStepsDataHas(NORMAL) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing NORMAL variable on solution step data","");
// 		if(mr_model_part.NodesBegin()->SolutionStepsDataHas(NODAL_AREA) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing NODAL_AREA variable on solution step data","");
	}
	

	


		
	
		

		

		
	
	ModelPart& mr_model_part;
	ModelPart* mtopographic_model_part_pointer;
	array_1d<double, 3 > mcalculation_domain_complete_displacement;
	array_1d<double, 3 > mcalculation_domain_added_displacement;
	bool mintialized_transfer_tool;
	bool muse_mesh_velocity_to_convect;
	int m_nparticles;
	int mnelems;
	double mDENSITY_WATER;
	double mDENSITY_AIR;
	
	//vector<double> mareas_vector; UNUSED SO COMMENTED 
	int max_nsubsteps;
	double max_substep_dt;
	int mmaximum_number_of_particles;
	std::vector< PFEM_Particle_Fluid  > mparticles_vector; //Point<3>
	int mlast_elem_id;
	bool modd_timestep;
	bool mparticle_printing_tool_initialized;
	unsigned int mfilter_factor;
	unsigned int mlast_node_id;
	//ModelPart& mr_particle_model_part;
	
	vector<int> mnumber_of_particles_in_elems; 
	vector<int> mnumber_of_particles_in_elems_aux; 
	vector<ParticlePointerVector*>  mpointers_to_particle_pointers_vectors;
	
	typename BinsObjectDynamic<Configure>::Pointer  mpBinsObjectDynamic;
	typename BinsObjectDynamic<Configure>::Pointer  mpTopographicBinsObjectDynamic;


	void CalculateNormal(Geometry<Node<3> >& pGeometry, array_1d<double,3>& An );
	
	};
	

	

	
}  // namespace Kratos.

#endif // KRATOS_INTEGRATE_UTILITY_INCLUDED  defined 
