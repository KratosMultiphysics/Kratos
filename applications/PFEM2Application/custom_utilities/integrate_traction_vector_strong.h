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

#if !defined(KRATOS_INTEGRATE_UTILITY_STRONG_INCLUDED)
#define  KRATOS_INTEGRATE_UTILITY_STRONG_INCLUDED



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



namespace Kratos
{
	template< unsigned int TDim>
	class IntegrateTractionVectorStrongUtility
	{
	public:
	
	  typedef SpatialContainersConfigure<TDim>     Configure;   
	  typedef typename Configure::PointType                      PointType; 
    typedef typename Configure::ContainerType                  ContainerType;   
    typedef typename Configure::IteratorType                   IteratorType; 
    typedef typename Configure::ResultContainerType            ResultContainerType;
    typedef typename Configure::ResultIteratorType             ResultIteratorType; 

		KRATOS_CLASS_POINTER_DEFINITION(IntegrateTractionVectorStrongUtility);
        

    struct element_data
    {
      bounded_matrix<double,TDim+1, TDim> v;
      array_1d<double,TDim+1> p, rho;
          
      bounded_matrix<double, TDim+1, TDim > DN_DX;
      array_1d<double, TDim+1 > N;
          
      Matrix C;
      Vector stress;
          
      double bdf0;
      double bdf1;
      double bdf2;
      double h;
    };

		//template<unsigned int TDim>
		IntegrateTractionVectorStrongUtility(ModelPart& model_part, double mu)
			: mr_model_part(model_part), mr_mu(mu)
		{
			std::cout << "initializing IntegrateTractionVectorStrongUtility" << std::endl;
			
			Check();
		
		}
		
         
		~IntegrateTractionVectorStrongUtility()
		{}
		
		vector<double> ReturnDrag()
        {
            KRATOS_TRY
            
            const SizeType NumNodes = TDim+1;
            element_data data; 
            double drag=0.0;
            double lift=0.0;
            double Volume = 0.0;
            double pressure=0.0;
            double Area=0.0;
            array_1d<double,TDim> traction_vector = ZeroVector(TDim); 
             
                  
            
            auto &rConditionsArray = mr_model_part.Conditions();
            const auto it_condition_begin = rConditionsArray.begin();
            array_1d<double, 3> dummy;
            int conditioncounter = 0;
            for (int i = 0; i < static_cast<int>(rConditionsArray.size()); ++i) 
            {
             auto it_condition = it_condition_begin + i;
             boost::numeric::ublas::bounded_matrix<double,2, 2> CouchyStressTensor=ZeroMatrix(2,2);
             boost::numeric::ublas::bounded_matrix<double,2, 2> GradientTensor=ZeroMatrix(2,2);
             boost::numeric::ublas::bounded_matrix<double,2, 2> GradientTensorTranspose=ZeroMatrix(2,2);
             boost::numeric::ublas::bounded_matrix<double,2, 2> PressureTensor=ZeroMatrix(2,2);
             array_1d<double,TDim> normal_vector = ZeroVector(TDim); 
             GlobalPointersVector< Element >& neighbor_els = it_condition->GetValue(NEIGHBOUR_ELEMENTS);
             if (it_condition->GetValue(NEIGHBOUR_ELEMENTS).size()>1)
               KRATOS_THROW_ERROR(std::invalid_argument,"NEIGHBOUR_ELEMENTS is more than one","");
             pressure=0.0;
             for(GlobalPointersVector< Element >::iterator ielem = neighbor_els.begin(); ielem!=neighbor_els.end(); ielem++)
             {

              Geometry< Node<3> >& geom = ielem->GetGeometry();
              double gx=0.3333*(geom[0].X()+geom[1].X()+geom[2].X());
              double gy=0.3333*(geom[0].Y()+geom[1].Y()+geom[2].Y());
              if (gx<0.51 && gx>-0.51 && gy<0.51 && gy>-0.51)
              {
              conditioncounter++;
              Geometry< Node<3> >& geom = ielem->GetGeometry();
              GeometryUtils::CalculateGeometryData(geom, data.DN_DX, data.N, Volume);
            //   const double mu = ielem->GetProperties()[VISCOSITY];
              for (unsigned int i = 0; i < NumNodes; i++)
              {
               const array_1d<double,TDim>& vel = ielem->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
               for(unsigned int k=0; k<TDim; k++)
               {
                data.v(i,k)   = vel[k];
               }
//                pressure += data.N[i]*ielem->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);
              }
              pressure+=0.5*(it_condition->GetGeometry()[0].FastGetSolutionStepValue(PRESSURE)+it_condition->GetGeometry()[1].FastGetSolutionStepValue(PRESSURE));
             Area=it_condition->GetGeometry().Length();
             GradientTensor=prod(trans(data.v),data.DN_DX);
             GradientTensorTranspose=trans(GradientTensor);
             PressureTensor(0,0)=pressure;
             PressureTensor(1,1)=pressure;
             CouchyStressTensor=2.0*mr_mu*(1.0/2.0)*(GradientTensor+GradientTensorTranspose)-PressureTensor;
             normal_vector = it_condition->GetGeometry().Normal(dummy);
             double normal_vector_magnitude=sqrt(normal_vector[0]*normal_vector[0]+normal_vector[1]*normal_vector[1]);
             normal_vector[0]=normal_vector[0]/normal_vector_magnitude;
             normal_vector[1]=normal_vector[1]/normal_vector_magnitude;
             drag+=(CouchyStressTensor(0,0)*normal_vector[0]+CouchyStressTensor(0,1)*normal_vector[1])*Area;
             lift+=(CouchyStressTensor(1,0)*normal_vector[0]+CouchyStressTensor(1,1)*normal_vector[1])*Area;
             }
             }
            }
            
            traction_vector[0]=drag;
            traction_vector[1]=lift;
            KRATOS_WATCH(conditioncounter);
            return traction_vector;

            KRATOS_CATCH("")
        }
		

		

		
	protected:


	private:
	
	
	void Check()
	{
		if(mr_model_part.NodesBegin()->SolutionStepsDataHas(VELOCITY) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data","");
		if(mr_model_part.NodesBegin()->SolutionStepsDataHas(PRESSURE) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE variable on solution step data","");
  }
	

	ModelPart& mr_model_part;
  double mr_mu;
	double mDENSITY_WATER;
	double mDENSITY_AIR;

  void CalculateShapeFunctionsOnGaussPoints(Geometry< Node<3> >& geom, array_1d<double,(TDim*2)>& gp, array_1d<double,(TDim+1)>& N,int& m)
  {
    
    double x10 = geom[1].X() - geom[0].X();
    double y10 = geom[1].Y() - geom[0].Y();
    double x20 = geom[2].X() - geom[0].X();
    double y20 = geom[2].Y() - geom[0].Y();
    double detJ = x10 * y20-y10 * x20;
    double TotalArea = 0.5*detJ;

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

  }

	
	};
	

	

	
}  // namespace Kratos.

#endif // KRATOS_INTEGRATE_UTILITY_STRONG_INCLUDED  defined 
