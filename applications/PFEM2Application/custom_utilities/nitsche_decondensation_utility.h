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
//   Last Modified by:    $Author: pbecker $
//   Date:                $Date: 2011-09-21 12:30:32 $
//   Revision:            $Revision: 1.0 $
//
//


#if !defined(KRATOS_NITSCHE_DECONDENSATION_UTILITY_INCLUDED)
#define  KRATOS_NITSCHE_DECONDENSATION_UTILITY_INCLUDED



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
	//this class is to be modified by the user to customize the interpolation process
	template< unsigned int TDim>
	class NitscheDecondensationUtility
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
        //typedef Configure::ContactPairType                ContactPairType;
        //typedef Configure::ContainerContactType           ContainerContactType; 
        //typedef Configure::IteratorContactType            IteratorContactType; 
        //typedef Configure::PointerContactType             PointerContactType; 
        //typedef Configure::PointerTypeIterator            PointerTypeIterator;

		KRATOS_CLASS_POINTER_DEFINITION(NitscheDecondensationUtility);
        

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
		NitscheDecondensationUtility(ModelPart& model_part, ModelPart& model_part2, double& nodeoffset)
			: mr_model_part(model_part) , mr_model_part2(model_part2) 
		{
			std::cout << "initializing NitscheDecondensation utility" << std::endl;
			
			Check();
            
            const bool printkratoswatch=0.0;
            
         
            if (printkratoswatch)
            KRATOS_WATCH("1");
            
			ModelPart::ElementsContainerType::iterator ielembegin = model_part.ElementsBegin();

            
            
		    const unsigned int NumNodes = 3;
            const unsigned int Dim = 2;
            const int ndofs = Dim + 1;
            const unsigned int MatrixSize = NumNodes*ndofs;
            unsigned int node_id = (model_part.Nodes().end() - 1)->Id() + 1;
            unsigned int element_id = (model_part.Elements().end() - 1)->Id() + 1;
            Matrix rPositiveEdgeIntersectionsShapeFunctionsValues;
            Matrix rNegativeEdgeIntersectionsShapeFunctionsValues;
            
            if (printkratoswatch)
            KRATOS_WATCH("2");
            
            model_part2.Nodes().clear();
            model_part2.Elements().clear();
            model_part2.Nodes()=model_part.Nodes();
            model_part2.Elements()=model_part.Elements();
            model_part2.pGetProperties(0)=model_part.pGetProperties(0);
            model_part2.pGetProperties(1)=model_part.pGetProperties(1);
            
            Node < 3 > ::Pointer pnode;

 
            double node0x=0.0;
            double node0y=0.0;
            double node1x=0.0;
            double node1y=0.0;
            double node2x=0.0;
            double node2y=0.0;
            double node3x=0.0;
            double node3y=0.0;
            double node4x=0.0;
            double node4y=0.0;
            double node5x=0.0;
            double node5y=0.0;
            double node6x=0.0;
            double node6y=0.0;
            double N3 = 0.0;
            double N4 = 0.0;
            double N5 = 0.0;
            double N6 = 0.0;
            double u3 = 0.0;
            double v3 = 0.0;
            double p3 = 0.0;
            double u4 = 0.0;
            double v4 = 0.0;
            double p4 = 0.0;
            double u5 = 0.0;
            double v5 = 0.0;
            double p5 = 0.0;
            double u6 = 0.0;
            double v6 = 0.0;
            double p6 = 0.0;
            
            
            if (printkratoswatch)
            KRATOS_WATCH("3");

            
			//using loop index, DO NOT paralelize this! change lines : mparticles_in_elems_pointers((ii*mmaximum_number_of_particles)+mparticles_in_elems_integers(ii)) = pparticle; and the next one

			for(unsigned int ii=0; ii<model_part.Elements().size(); ii++)
			{

              ModelPart::ElementsContainerType::iterator ielem = ielembegin+ii; 
              Geometry< Node<3> >& geom = ielem->GetGeometry();
              array_1d<double, NumNodes> distances;
              
              if (printkratoswatch)
              KRATOS_WATCH("4");
//               
//               KRATOS_WATCH(geom);
//               KRATOS_WATCH(distances);
              
              
              double Volume;
              boost::numeric::ublas::bounded_matrix<double,(3), 2 > coords;
            element_data data;
            boost::numeric::ublas::bounded_matrix<double, 3*(2-1), (2+1) > Ngauss;   
            boost::numeric::ublas::bounded_matrix<double,3, 5> Nenriched1;
            std::vector< Matrix > gauss_gradients(3);
            array_1d<double,(3*(2-1))> volumes;
            array_1d<double,(3*(2-1))> signs; //ATTENTION: this shall be initialized of size 6
            boost::numeric::ublas::bounded_matrix<double,10, 2> rGradientpositive;
            boost::numeric::ublas::bounded_matrix<double,10, 2> rGradientnegative;
            boost::numeric::ublas::bounded_matrix<int,3,3> father_nodes;
            boost::numeric::ublas::bounded_matrix<int,3,3> subelementnodes;
            boost::numeric::ublas::bounded_matrix<double,7,2> nodecoordinates;
            array_1d<int,7> subelementnodeIds;
            array_1d<double,(12)> ui = ZeroVector(12);
            boost::numeric::ublas::bounded_matrix<double,2, 2> nodes=ZeroMatrix(2,2);
            
            const bool reduction = false;
              
            if (printkratoswatch)
            KRATOS_WATCH("5");
              
              std::vector< Matrix > gauss_gradients_discon(3);
            for (unsigned int i = 0; i < 3; i++)
        {
          gauss_gradients_discon[i].resize(7, 2, false);
        }
        for (unsigned int i = 0; i < 3; i++)
         for (unsigned int j = 0; j < 7; j++)
          for (unsigned int k = 0; k < 2; k++)
           gauss_gradients_discon[i](j,k) = 0.0; 
        boost::numeric::ublas::bounded_matrix<double,3, 7> Nenriched1_discon;
        for (unsigned int i = 0; i < 3; i++)
         for (unsigned int j = 0; j < 7; j++)
          Nenriched1_discon(i,j) = 0.0;  

        unsigned int npos=0, nneg=0;
        for (unsigned int i = 0; i < NumNodes; i++)
        {   
//           KRATOS_WATCH(ielem->GetGeometry()[i].X());
//           KRATOS_WATCH(ielem->GetGeometry()[i].Y());
          distances[i] = ielem->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
          if(distances[i] > 0.0)
            npos++;
          if (distances[i] < 0.0) 
            nneg++;
        }
        
//         KRATOS_WATCH(distances);
//         KRATOS_WATCH(npos);
//         KRATOS_WATCH(nneg);
//         KRATOS_WATCH("6");
        
        if(npos != NumNodes && nneg != NumNodes)
        {
               
          
          if (printkratoswatch)
          KRATOS_WATCH("7");
                
          GeometryUtils::CalculateGeometryData(geom, data.DN_DX, data.N, Volume);
          for (unsigned int i = 0; i < NumNodes; i++)
          {
            const array_1d<double, 3 > & xyz = geom[i].Coordinates();
            for (unsigned int j = 0; j < Dim; j++)
              coords(i, j) = xyz[j];
          }
          
          unsigned int ndivisions = CalculateEnrichedShapeFuncionsExtendedmodified(coords, data.DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched1, rGradientpositive , rGradientnegative, father_nodes);
          //                   KRATOS_WATCH("I'm gonna remove an element");
          
//           KRATOS_WATCH(geom);
//           KRATOS_WATCH(coords);
//           KRATOS_WATCH(data.DN_DX);
//           KRATOS_WATCH(distances);
//           KRATOS_WATCH(volumes);
//           KRATOS_WATCH(Ngauss);
//           KRATOS_WATCH(signs);
//           KRATOS_WATCH(gauss_gradients[0]);
//           KRATOS_WATCH(gauss_gradients[1]);
//           KRATOS_WATCH(gauss_gradients[2]);
//           KRATOS_WATCH(Nenriched1);
//           KRATOS_WATCH(rGradientpositive);
//           KRATOS_WATCH(rGradientnegative);
//           KRATOS_WATCH(father_nodes);
          if (printkratoswatch)
          KRATOS_WATCH("8");
          
//           Geometry<Node<3>>::Pointer p_geometry = ielem->pGetGeometry();
//           const Vector& r_elemental_distances = ielem->GetValue(ELEMENTAL_DISTANCES);
//           Triangle2D3ModifiedShapeFunctions triangle_shape_functions(p_geometry, r_elemental_distances);
//           triangle_shape_functions.ComputeShapeFunctionsOnPositiveEdgeIntersections(rPositiveEdgeIntersectionsShapeFunctionsValues);
// 
//           triangle_shape_functions.ComputeShapeFunctionsOnNegativeEdgeIntersections(rNegativeEdgeIntersectionsShapeFunctionsValues);
          
          
//           KRATOS_WATCH("9");
          
          
          int index = 0;
          for(int aux=0;aux<3;aux++)
          {
            if(father_nodes(aux,0)!=father_nodes(aux,1))
            {
              nodes(index,0)=father_nodes(aux,0);
              nodes(index,1)=father_nodes(aux,1);
              index++;
            }
          }
           if (printkratoswatch)
           KRATOS_WATCH("10");
//            double nodeoffset=0.001;
           
           
           
           //For Side-1
           array_1d<double,2> interface_segment=ZeroVector(2);
//            array_1d<double,2> normaledge1=ZeroVector(2);
//            array_1d<double,2> normaledge2=ZeroVector(2);
           int indexf=nodes(0,0);
           int indexs=nodes(0,1);
           double norm=0.0;
           double areaplus1=0.0;
           double areaminus1=0.0;
           double toleranceside1=0.00;
           double disttocenter = 0.0;
           double mside1=0.0;
           double distx1=0.0;
           double sinalpha1=0.0;
           double cosalpha1=0.0;
           interface_segment[0] = (ielem->GetGeometry()[indexf].X()-ielem->GetGeometry()[indexs].X());
           interface_segment[1] = (ielem->GetGeometry()[indexf].Y()-ielem->GetGeometry()[indexs].Y());
           norm=sqrt( pow((interface_segment[0]),2)+pow((interface_segment[1]),2));
           double area1=norm;
           
//            KRATOS_WATCH(indexf);
//            KRATOS_WATCH(indexs);
//            KRATOS_WATCH(area1);

           if (printkratoswatch)
           KRATOS_WATCH("11");

           boost::numeric::ublas::bounded_matrix<double, 2, 2 > InterfacePoints;
           double relative_position = fabs(ielem->GetGeometry()[indexs].FastGetSolutionStepValue(DISTANCE)/
           (ielem->GetGeometry()[indexs].FastGetSolutionStepValue(DISTANCE)-ielem->GetGeometry()[indexf].FastGetSolutionStepValue(DISTANCE)));
           InterfacePoints(0,0) = relative_position*ielem->GetGeometry()[indexf].X()+(1.0-relative_position)*ielem->GetGeometry()[indexs].X();
           InterfacePoints(0,1) = relative_position*ielem->GetGeometry()[indexf].Y()+(1.0-relative_position)*ielem->GetGeometry()[indexs].Y();

           node3x=InterfacePoints(0,0); 
           node3y=InterfacePoints(0,1);
           node5x=InterfacePoints(0,0); 
           node5y=InterfacePoints(0,1);
           
           if (printkratoswatch)
           KRATOS_WATCH("12");
           
           if(ielem->GetGeometry()[indexf].FastGetSolutionStepValue(DISTANCE)>0.0)
           {
           interface_segment[0] = (ielem->GetGeometry()[indexf].X()-InterfacePoints(0,0));
           interface_segment[1] = (ielem->GetGeometry()[indexf].Y()-InterfacePoints(0,1));
           norm=sqrt( pow((interface_segment[0]),2)+pow((interface_segment[1]),2));
           areaplus1=norm;
           areaminus1=area1-areaplus1;
           
           N3=areaminus1/area1;
           N5=areaplus1/area1;
           u3=N3*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(VELOCITY_X)+N5*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(VELOCITY_NITSCHE_X);
           v3=N3*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(VELOCITY_Y)+N5*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(VELOCITY_NITSCHE_Y);
           p3=N3*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(PRESSURE)+N5*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(PRESSURE_NITSCHE);
           
           u5=N5*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(VELOCITY_X)+N3*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(VELOCITY_NITSCHE_X);
           v5=N5*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(VELOCITY_Y)+N3*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(VELOCITY_NITSCHE_Y);
           p5=N5*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(PRESSURE)+N3*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(PRESSURE_NITSCHE);
           
           if (printkratoswatch)
           KRATOS_WATCH("13");
           
           if (reduction==true)
           {
            toleranceside1=0.01*areaminus1; 
            sinalpha1=(ielem->GetGeometry()[indexs].Y()-node5y)/areaminus1;
            cosalpha1=(ielem->GetGeometry()[indexs].X()-node5x)/areaminus1;
            
            disttocenter=sqrt((node5y-1.5)*(node5y-1.5)+(node5x-0.5)*(node5x-0.5));
            sinalpha1=(1.5-node5y)/disttocenter;
            cosalpha1=(0.5-node5x)/disttocenter;

            node5x=node5x+nodeoffset*disttocenter*cosalpha1; 
            node5y=node5y+nodeoffset*disttocenter*sinalpha1; 

            

           }
           
           if (printkratoswatch)
           KRATOS_WATCH("14");
           
           //u3=n3_indexf*GetGeometry()[indexf]velocityu+n3_indexs*GetGeometry()[indexs]velocityunitsche
//                                  KRATOS_WATCH(areaplus1);
//            KRATOS_WATCH(areaminus1);
           }
           else
           {
           interface_segment[0] = (ielem->GetGeometry()[indexf].X()-InterfacePoints(0,0));
           interface_segment[1] = (ielem->GetGeometry()[indexf].Y()-InterfacePoints(0,1));
           norm=sqrt( pow((interface_segment[0]),2)+pow((interface_segment[1]),2));
           areaminus1=norm;
           areaplus1=area1-areaminus1;
           N3=areaminus1/area1;
           N5=areaplus1/area1;
           u3=N3*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(VELOCITY_X)+N5*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(VELOCITY_NITSCHE_X);
           v3=N3*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(VELOCITY_Y)+N5*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(VELOCITY_NITSCHE_Y);
           p3=N3*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(PRESSURE)+N5*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(PRESSURE_NITSCHE);
           
           u5=N5*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(VELOCITY_X)+N3*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(VELOCITY_NITSCHE_X);
           v5=N5*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(VELOCITY_Y)+N3*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(VELOCITY_NITSCHE_Y);
           p5=N5*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(PRESSURE)+N3*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(PRESSURE_NITSCHE);
          
           if (printkratoswatch)
           KRATOS_WATCH("15");
           
           if (reduction==true)
           {
            toleranceside1=0.01*areaminus1; 
            sinalpha1=(ielem->GetGeometry()[indexf].Y()-node5y)/areaminus1;
            cosalpha1=(ielem->GetGeometry()[indexf].X()-node5x)/areaminus1;
            
            disttocenter=sqrt((node5y-1.5)*(node5y-1.5)+(node5x-0.5)*(node5x-0.5));
            sinalpha1=(1.5-node5y)/disttocenter;
            cosalpha1=(0.5-node5x)/disttocenter;

            node5x=node5x+nodeoffset*disttocenter*cosalpha1; 
            node5y=node5y+nodeoffset*disttocenter*sinalpha1; 
            
           }
           
           if (printkratoswatch)
           KRATOS_WATCH("16");
           
           
           
          }
//            KRATOS_WATCH(areaplus1);
//            KRATOS_WATCH(areaminus1);
           
           //For Side-2
           interface_segment=ZeroVector(2);
           indexf=nodes(1,0);
           indexs=nodes(1,1); 
           norm=0.0;
           double areaplus2=0.0;
           double areaminus2=0.0;
           double toleranceside2=0.00;
           double mside2=0.0;
           double distx2=0.0;
           double sinalpha2=0.0;
           double cosalpha2=0.0;
           interface_segment[0] = (ielem->GetGeometry()[indexf].X()-ielem->GetGeometry()[indexs].X());
           interface_segment[1] = (ielem->GetGeometry()[indexf].Y()-ielem->GetGeometry()[indexs].Y());
           norm=sqrt( pow((interface_segment[0]),2)+pow((interface_segment[1]),2));
           double area2=norm;

           if (printkratoswatch)
           KRATOS_WATCH("17");


           relative_position = fabs(ielem->GetGeometry()[indexs].FastGetSolutionStepValue(DISTANCE)/
           (ielem->GetGeometry()[indexs].FastGetSolutionStepValue(DISTANCE)-ielem->GetGeometry()[indexf].FastGetSolutionStepValue(DISTANCE)));
           InterfacePoints(0,0) = relative_position*ielem->GetGeometry()[indexf].X()+(1.0-relative_position)*ielem->GetGeometry()[indexs].X();
           InterfacePoints(0,1) = relative_position*ielem->GetGeometry()[indexf].Y()+(1.0-relative_position)*ielem->GetGeometry()[indexs].Y();

           node4x=InterfacePoints(0,0); 
           node4y=InterfacePoints(0,1);
           node6x=InterfacePoints(0,0); 
           node6y=InterfacePoints(0,1); 
           
           if (printkratoswatch)
           KRATOS_WATCH("18");
           
           if(ielem->GetGeometry()[indexf].FastGetSolutionStepValue(DISTANCE)>0.0)
           {
           interface_segment[0] = (ielem->GetGeometry()[indexf].X()-InterfacePoints(0,0));
           interface_segment[1] = (ielem->GetGeometry()[indexf].Y()-InterfacePoints(0,1));
           norm=sqrt( pow((interface_segment[0]),2)+pow((interface_segment[1]),2));
           areaplus2=norm;
           areaminus2=area2-areaplus2;  
           N4=areaminus2/area2;
           N6=areaplus2/area2;
           u4=N4*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(VELOCITY_X)+N6*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(VELOCITY_NITSCHE_X);
           v4=N4*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(VELOCITY_Y)+N6*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(VELOCITY_NITSCHE_Y);
           p4=N4*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(PRESSURE)+N6*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(PRESSURE_NITSCHE);
           
           u6=N6*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(VELOCITY_X)+N4*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(VELOCITY_NITSCHE_X);
           v6=N6*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(VELOCITY_Y)+N4*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(VELOCITY_NITSCHE_Y);
           p6=N6*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(PRESSURE)+N4*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(PRESSURE_NITSCHE);
//            KRATOS_WATCH(areaplus2);
//            KRATOS_WATCH(areaminus2);
           
           if (printkratoswatch)
           KRATOS_WATCH("19");
           
           if (reduction==true)
           {


            sinalpha2=(ielem->GetGeometry()[indexs].Y()-node6y)/areaminus2;
            cosalpha2=(ielem->GetGeometry()[indexs].X()-node6x)/areaminus2;
            
            disttocenter=sqrt((node6y-1.5)*(node6y-1.5)+(node6x-0.5)*(node6x-0.5));
            sinalpha2=(1.5-node6y)/disttocenter;
            cosalpha2=(0.5-node6x)/disttocenter;

            node6x=node6x+nodeoffset*disttocenter*cosalpha2; 
            node6y=node6y+nodeoffset*disttocenter*sinalpha2; 
           }
           
           
           }
           else
           {
           interface_segment[0] = (ielem->GetGeometry()[indexf].X()-InterfacePoints(0,0));
           interface_segment[1] = (ielem->GetGeometry()[indexf].Y()-InterfacePoints(0,1));
           norm=sqrt( pow((interface_segment[0]),2)+pow((interface_segment[1]),2));
           areaminus2=norm;
           areaplus2=area2-areaminus2; 
           N4=areaminus2/area2;
           N6=areaplus2/area2;
           u4=N4*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(VELOCITY_X)+N6*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(VELOCITY_NITSCHE_X);
           v4=N4*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(VELOCITY_Y)+N6*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(VELOCITY_NITSCHE_Y);
           p4=N4*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(PRESSURE)+N6*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(PRESSURE_NITSCHE);
           
           u6=N6*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(VELOCITY_X)+N4*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(VELOCITY_NITSCHE_X);
           v6=N6*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(VELOCITY_Y)+N4*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(VELOCITY_NITSCHE_Y);
           p6=N6*ielem->GetGeometry()[indexf].FastGetSolutionStepValue(PRESSURE)+N4*ielem->GetGeometry()[indexs].FastGetSolutionStepValue(PRESSURE_NITSCHE);
//            KRATOS_WATCH(areaplus2);
//            KRATOS_WATCH(areaminus2);
           
           if (printkratoswatch)
           KRATOS_WATCH("20");
             
           if (reduction==true)
           {


            sinalpha2=(ielem->GetGeometry()[indexf].Y()-node6y)/areaminus2;
            cosalpha2=(ielem->GetGeometry()[indexf].X()-node6x)/areaminus2;
            
            disttocenter=sqrt((node6y-1.5)*(node6y-1.5)+(node6x-0.5)*(node6x-0.5));
            
            sinalpha2=(1.5-node6y)/disttocenter;
            cosalpha2=(0.5-node6x)/disttocenter;

            node6x=node6x+nodeoffset*disttocenter*cosalpha2; 
            node6y=node6y+nodeoffset*disttocenter*sinalpha2; 
           }
           
           
           
          }
          
          if (printkratoswatch)
          KRATOS_WATCH("21");
             
//            KRATOS_WATCH(u3);
//            KRATOS_WATCH(v3);
//            KRATOS_WATCH(p3);

//                   KRATOS_WATCH(rPositiveEdgeIntersectionsShapeFunctionsValues);
//        KRATOS_WATCH(rNegativeEdgeIntersectionsShapeFunctionsValues);
// if (ielem->Id()==762)
//            KRATOS_THROW_ERROR(std::logic_error,"the case is not implemented yet","");
//           KRATOS_WATCH(ndivisions);
       for (unsigned int i=0;i<ndivisions;i++)
       {
         for (unsigned int j=0;j<3;j++)
         {

           gauss_gradients_discon[i](j,0)=gauss_gradients[i](j,0);

           gauss_gradients_discon[i](j,1)=gauss_gradients[i](j,1);

           Nenriched1_discon(i,j) = Nenriched1(i,j);

         }
         
         if (printkratoswatch)
         KRATOS_WATCH("22");
         
         if(signs(i)>0.0)
         {
           gauss_gradients_discon[i](3,0)=gauss_gradients[i](3,0);
           gauss_gradients_discon[i](3,1)=gauss_gradients[i](3,1);
           gauss_gradients_discon[i](4,0)=gauss_gradients[i](4,0);
           gauss_gradients_discon[i](4,1)=gauss_gradients[i](4,1);
           Nenriched1_discon(i,3) = Nenriched1(i,3);
           Nenriched1_discon(i,4) = Nenriched1(i,4);

         }
         else
         {
           gauss_gradients_discon[i](5,0)=gauss_gradients[i](3,0);
           gauss_gradients_discon[i](5,1)=gauss_gradients[i](3,1);
           gauss_gradients_discon[i](6,0)=gauss_gradients[i](4,0);
           gauss_gradients_discon[i](6,1)=gauss_gradients[i](4,1);
           Nenriched1_discon(i,5) = Nenriched1(i,3);
           Nenriched1_discon(i,6) = Nenriched1(i,4);
         }
       }
       
       if (printkratoswatch)
       KRATOS_WATCH("23");
       
       subelementnodeIds(0)=ielem->GetGeometry()[0].Id();
       //         r_model_part2.CreateNewNode(ielem->GetGeometry()[1].Id(), ielem->GetGeometry()[1].X(), ielem->GetGeometry()[1].Y(), 0.0);
       subelementnodeIds(1)=ielem->GetGeometry()[1].Id();
       //         r_model_part2.CreateNewNode(ielem->GetGeometry()[2].Id(), ielem->GetGeometry()[2].X(), ielem->GetGeometry()[2].Y(), 0.0);
       subelementnodeIds(2)=ielem->GetGeometry()[2].Id();
       //         KRATOS_WATCH(node0x);
       node0x=ielem->GetGeometry()[0].X();
       node0y=ielem->GetGeometry()[0].Y();
       node1x=ielem->GetGeometry()[1].X();
       node1y=ielem->GetGeometry()[1].Y();
       node2x=ielem->GetGeometry()[2].X();
       node2y=ielem->GetGeometry()[2].Y();
       //         KRATOS_WATCH(node2x);
       pnode = mr_model_part2.CreateNewNode(node_id, node3x, node3y, 0.0);
       pnode->FastGetSolutionStepValue(VELOCITY_X) =u3;
       pnode->FastGetSolutionStepValue(VELOCITY_Y) =v3;
       pnode->FastGetSolutionStepValue(PRESSURE) =p3;
       pnode->FastGetSolutionStepValue(DISTANCE) =100.0;
       subelementnodeIds(3)=node_id;
       node_id++;
       pnode = mr_model_part2.CreateNewNode(node_id, node4x, node4y, 0.0);
       pnode->FastGetSolutionStepValue(VELOCITY_X) =u4;
       pnode->FastGetSolutionStepValue(VELOCITY_Y) =v4;
       pnode->FastGetSolutionStepValue(PRESSURE) =p4;
       pnode->FastGetSolutionStepValue(DISTANCE) =100.0;
       //         pnode->FastGetSolutionStepValue(VELOCITY_X) =ui(0);
       subelementnodeIds(4)=node_id;
       node_id++;
       pnode = mr_model_part2.CreateNewNode(node_id, node5x, node5y, 0.0);
       pnode->FastGetSolutionStepValue(VELOCITY_X) =u5;
       pnode->FastGetSolutionStepValue(VELOCITY_Y) =v5;
       pnode->FastGetSolutionStepValue(PRESSURE) =p5;
       pnode->FastGetSolutionStepValue(DISTANCE) =-100.0;
       //         pnode->FastGetSolutionStepValue(VELOCITY_X) =ui(0);
       subelementnodeIds(5)=node_id;
       node_id++;
       pnode = mr_model_part2.CreateNewNode(node_id, node6x, node6y, 0.0);
       pnode->FastGetSolutionStepValue(VELOCITY_X) =u6;
       pnode->FastGetSolutionStepValue(VELOCITY_Y) =v6;
       pnode->FastGetSolutionStepValue(PRESSURE) =p6;
       pnode->FastGetSolutionStepValue(DISTANCE) =-100.0;
       //         pnode->FastGetSolutionStepValue(VELOCITY_X) =ui(0);
       subelementnodeIds(6)=node_id;
       node_id++;
       
       nodecoordinates(0,0)=node0x;
       nodecoordinates(0,1)=node0y;
       
       nodecoordinates(1,0)=node1x;
       nodecoordinates(1,1)=node1y;
       
       nodecoordinates(2,0)=node2x;
       nodecoordinates(2,1)=node2y;
       
       nodecoordinates(3,0)=node3x;
       nodecoordinates(3,1)=node3y;
       
       nodecoordinates(4,0)=node4x;
       nodecoordinates(4,1)=node4y;
       
       nodecoordinates(5,0)=node5x;
       nodecoordinates(5,1)=node5y;
       
       nodecoordinates(6,0)=node6x;
       nodecoordinates(6,1)=node6y;
       
       nodecoordinates(6,0)=node6x;
       nodecoordinates(6,1)=node6y;
       
//        KRATOS_WATCH("24");

//        KRATOS_THROW_ERROR(std::logic_error,"the case is not implemented yet","");
       
                  
       Properties::Pointer p_properties = model_part.pGetProperties(0);
       std::vector<ModelPart::IndexType> ElementNodeIds;
       ElementNodeIds.resize(3);
       ModelPart::IndexType default_index = 0;
//          KRATOS_WATCH(mr_model_part2.Elements());
       model_part2.RemoveElement(ielem->Id(),default_index);
       
//        KRATOS_WATCH("25");
        
//         KRATOS_WATCH(mr_model_part2.Elements());
       for (int i=0;i<3;i++)
        {
          bounded_matrix<double,3+4,2> DN_sub_aux=ZeroMatrix(3+4,2);  
          bounded_matrix<double,3,2> DN_sub=ZeroMatrix(3,2); 
          for (int j=0;j<3+4;j++)
          {
           for (int k=0;k<2;k++)
           {
            DN_sub(j,k) = gauss_gradients_discon[i](j,k);
           }
          }
          
//           KRATOS_WATCH(DN_sub);
          int counter=0;
          for (int j=0;j<7;j++)
          {
            if (gauss_gradients_discon[i](j,0)==0.0 and gauss_gradients_discon[i](j,1)==0.0)
            { 
            }
            else
            {
            subelementnodes(i,counter)=j;
            DN_sub(counter,0)=gauss_gradients_discon[i](j,0);
            DN_sub(counter,1)=gauss_gradients_discon[i](j,1);
            counter++;  
            }
          }
          

//           KRATOS_WATCH("26");
          
//           KRATOS_WATCH(i);
          
              
//          mr_model_part2.CreateNewElement("Element2D3N", element_id, [subelementnodes(i,0),subelementnodes(i,1),subelementnodes(i,2)], r_model_part.GetProperties()[1]);
          int a=subelementnodes(i,0);
          int b=subelementnodes(i,1);
          int c=subelementnodes(i,2);
          
//           KRATOS_WATCH(a);
//           KRATOS_WATCH(b);
//           KRATOS_WATCH(c);

//           ElementNodeIds[0]=subelementnodeIds(a);
//           ElementNodeIds[1]=subelementnodeIds(b);
//           ElementNodeIds[2]=subelementnodeIds(c);
          ElementNodeIds[0]=subelementnodeIds(a);
          ElementNodeIds[1]=subelementnodeIds(b);
          ElementNodeIds[2]=subelementnodeIds(c);
          
//           KRATOS_WATCH("27");

         ModelPart::ElementType::Pointer p_new_element = model_part2.CreateNewElement("MonolithicPFEM22DDenizNitsche", element_id, ElementNodeIds, p_properties, default_index);

//          

//          KRATOS_WATCH("28");
                

//          KRATOS_WATCH(ElementNodeIds);
//           mr_model_part2.CreateNewElement();
//          ElementType::Pointer p_new_element = mr_model_part2->CreateNewElement("Element2D3N", element_id, [subelementnodeIds(a),subelementnodeIds(b),subelementnodeIds(c)], r_model_part.GetProperties()[1]);
         element_id++;
         
//          KRATOS_WATCH("29");

        }

        

        //positive side
        
                  
                  
                  
                }
              
              
            }
            
            
            
           
//       KRATOS_WATCH("29");
			
			

		}
		
         
		~NitscheDecondensationUtility()
		{}

		void MountBin()
		{
			KRATOS_TRY

			//copy the elements to a new container, as the list will
			//be shuffled duringthe construction of the tree
			ContainerType& rElements           =  mr_model_part.ElementsArray();
	        IteratorType it_begin              =  rElements.begin();
	        IteratorType it_end                =  rElements.end();
	        //const int number_of_elem 		   =   rElements.size();

			typename BinsObjectDynamic<Configure>::Pointer paux = typename BinsObjectDynamic<Configure>::Pointer(new BinsObjectDynamic<Configure>(it_begin, it_end  ) );
			paux.swap(mpBinsObjectDynamic);
			//BinsObjectDynamic<Configure>  mpBinsObjectDynamic(it_begin, it_end ); 
				
			std::cout << "finished mounting Bins" << std::endl;

			KRATOS_CATCH("")
		}
		
		static int CalculateEnrichedShapeFuncionsExtendedmodified(boost::numeric::ublas::bounded_matrix<double,(2+1), 2 >& rPoints, boost::numeric::ublas::bounded_matrix<double, (2+1), 2 >& DN_DX,
                                  array_1d<double,(2+1)>& rDistances, array_1d<double,(3*(2-1))>& rVolumes, boost::numeric::ublas::bounded_matrix<double, 3*(2-1), (2+1) >& rGPShapeFunctionValues, array_1d<double,(3*(2-1))>& rPartitionsSign, std::vector<Matrix>& rGradientsValue, boost::numeric::ublas::bounded_matrix<double,3*(2-1), (5)>& NEnriched,boost::numeric::ublas::bounded_matrix<double,10, 2>& rGradientpositive,boost::numeric::ublas::bounded_matrix<double,10, 2>& rGradientnegative ,boost::numeric::ublas::bounded_matrix<int,3,3>& father_nodes) 
        {
      KRATOS_TRY
        
      const double one_third=1.0/3.0;
      boost::numeric::ublas::bounded_matrix<double,3,2> aux_points; //for auxiliary nodes 4(between 1 and 2) ,5(between 2 and 3) ,6 (between 3 and 1)
      boost::numeric::ublas::bounded_matrix<double, 3, 2 > coord_subdomain; //used to pass arguments when we must calculate areas, shape functions, etc
      boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX_subdomain; //used to retrieve derivatives




       
      double most_common_sign=0; //the side of the cut in which two nodes are found (same sign) will be the ones that remains unchanged when builing the discontinuity
      double Area;//area of the complete element
      rGPShapeFunctionValues(0,0)=one_third; rGPShapeFunctionValues(0,1)=one_third; rGPShapeFunctionValues(0,2)=one_third; //default, when no interfase has been found
      Area = CalculateVolume2D( rPoints );
      array_1d<bool,3> cut_edges;
      array_1d<double,3> aux_nodes_relative_locations;
      boost::numeric::ublas::bounded_matrix<int,3,2> aux_nodes_father_nodes;
      
      //to begin with we must check whether our element is cut or not by the interfase.
      if( (rDistances(0)*rDistances(1))>0.0 && (rDistances(0)*rDistances(2))>0.0 ) //it means that this element IS NOT cut by the interfase. we must return data of a normal, non-enriched element
        {
          rVolumes(0)=Area;
          rGPShapeFunctionValues(0,0)=one_third; rGPShapeFunctionValues(0,1)=one_third; rGPShapeFunctionValues(0,2)=one_third;
          NEnriched(0,0) = 0.0;
          //type_of_cut=1;
          for (int j = 0; j < 2; j++)
                rGradientsValue[0](0, j) = 0.0;
          if (rDistances(0) < 0.0) rPartitionsSign[0] = -1.0;
          else rPartitionsSign[0] = 1.0;
            return 1;
        }

      //else //we must create the enrichement, it can be in 2 or 3 parts. we'll start with 3 always.

      //'TRICK' TO AVOID HAVING THE INTERFASE TOO CLOSE TO THE NODES:
      //since we cannot collapse node because we have to contemplate the possibility of discontinuities, we will move a little the intefase so that it is not that close.
      
      array_1d<double,3>  rDistances_aux;
      rDistances_aux(0)=rDistances(0);
      rDistances_aux(1)=rDistances(1);
      rDistances_aux(2)=rDistances(2);
      

      
      array_1d<double,3>  sign_aux;
      double sign_side;
      sign_aux(0)=1.0;
      sign_aux(1)=1.0;
      sign_aux(2)=1.0;
      if (rDistances(0)<0.0)
        sign_aux(0)=-1.0;
      if (rDistances(1)<0.0)
        sign_aux(1)=-1.0;
      if (rDistances(2)<0.0)
        sign_aux(2)=-1.0;
       
      const double unsigned_distance0=fabs(rDistances(0));
      const double unsigned_distance1=fabs(rDistances(1));
      const double unsigned_distance2=fabs(rDistances(2));
      



      //we begin by finding the largest distance:
      double longest_distance=fabs(unsigned_distance0);
      if (unsigned_distance1>longest_distance)
        longest_distance=unsigned_distance1;
      if (unsigned_distance2>longest_distance)
        longest_distance=unsigned_distance2;
//       double h=0.15;
// //       CALCULATION OF THE H - START
    double inv_h_max = 0.0;
    for(unsigned int i=0; i<3; i++)
    {
        double inv_h = 0.0;
        for(unsigned int k=0; k<2; k++)
            inv_h += DN_DX(i,k)*DN_DX(i,k);

        if(inv_h > inv_h_max) inv_h_max = inv_h;
    }
    inv_h_max = sqrt(inv_h_max);
    double h = 1.0/inv_h_max;
    
    h=0.1;
    
    


// //       CALCULATION OF THE H - END
      const double tolerable_distance =h*0.0;
      const double ctolerable_distance =h*0.1;
      

      

      double changecounter=0.0;
      double diff=0.0;
      array_1d<double,3> changeswitch;
      changeswitch(0)=1.0;
      changeswitch(1)=1.0;
      changeswitch(2)=1.0;
      array_1d<double,3> difference=ZeroVector(3);
      

      
      //const double epsilon = 1e-15; //1.00e-9;
      //compute the gradient of the distance and normalize it
      array_1d<double, 2 > grad_d;
      noalias(grad_d) = prod(trans(DN_DX), rDistances);
      
      
      array_1d<double, 3> exact_distance = rDistances;
      array_1d<double, 3> abs_distance = ZeroVector(3);

      //KRATOS_WATCH("one element IS in the intefase")
      if ((rDistances(0)*rDistances(1))<0.0) //edge 12 is cut
        cut_edges[0]=true;
      else
        cut_edges[0]=false;
      if ((rDistances(1)*rDistances(2))<0.0) //edge 23 is cut. 
        cut_edges[1]=true;
      else
        cut_edges[1]=false;
      if ((rDistances(2)*rDistances(0))<0.0) //edge 23 is cut. 
        cut_edges[2]=true;
      else
        cut_edges[2]=false;
      
      //We have 3 edges, meaning we created 3 aux nodes. But one of them actually matches the position of a real node (the one that is not on an interface edge is displaced to one of the ends (a node)
      //the new shape functions are built by setting the values in all (real and aux) nodes to zero except in one of the interphase nodes. in the array aux_node_shape_function_index we assign this value
      array_1d<int, 3 > aux_node_enrichment_shape_function_index; //when not used, it must be -1;
      
      

      int shape_function_id=0;
      father_nodes(0,0)=-1;
      father_nodes(0,1)=-1;
      father_nodes(0,2)=-1;
      father_nodes(1,0)=-1;
      father_nodes(1,1)=-1;
      father_nodes(1,2)=-1;
      father_nodes(2,0)=-1;
      father_nodes(2,1)=-1;
      father_nodes(2,2)=-1;

      //KRATOS_WATCH(father_nodes);
      for (unsigned int i=0; i<3; i++) //we go over the 3 edges:
        {
          int edge_begin_node=i;
          int edge_end_node=i+1;
          if (edge_end_node==3) edge_end_node=0; //it's a triangle, so node 3 is actually node 0
            
          if(cut_edges(i)==true)
        {
          aux_nodes_relative_locations(i)=fabs(rDistances(edge_end_node)/(rDistances(edge_end_node)-rDistances(edge_begin_node) ) ) ; //position in 'natural' coordinates of edge 12, 1 when it passes over node 1. (it is over the edge 01)
          //KRATOS_WATCH(aux_nodes_relative_locations(i));
          
          aux_nodes_father_nodes(i,0)=edge_begin_node;
          aux_nodes_father_nodes(i,1)=edge_end_node;
          
          aux_node_enrichment_shape_function_index(i)=shape_function_id;
          father_nodes(i,0)=edge_begin_node;
          father_nodes(i,1)=edge_end_node;
          father_nodes(i,2)=shape_function_id;
          shape_function_id++;
          

        }
          else
        {//<
//           if(fabs(rDistances(edge_end_node))<fabs(rDistances(edge_begin_node)))    //old
          if(fabs(rDistances(edge_end_node))>fabs(rDistances(edge_begin_node))) //if edge is not cut, we collapse the aux node into the node which has the highest absolute value to have "nicer" (less "slivery") subelements
            {
              aux_nodes_relative_locations(i)=0.0;
              aux_nodes_father_nodes(i,0)=edge_end_node;
              aux_nodes_father_nodes(i,1)=edge_end_node;
            }
          else
            {
              aux_nodes_relative_locations(i)=1.0;
              aux_nodes_father_nodes(i,0)=edge_begin_node;
              aux_nodes_father_nodes(i,1)=edge_begin_node;
            }
          //KRATOS_WATCH("compileddd");
          aux_node_enrichment_shape_function_index(i)=-1;
        }
          
          //and we save the coordinate of the new aux nodes:
          for (unsigned int j=0;j<2;j++){   //x,y coordinates
        aux_points(i,j)= rPoints(edge_begin_node,j) * aux_nodes_relative_locations(i) + rPoints(edge_end_node,j) * (1.0- aux_nodes_relative_locations(i));
          }
        }

      array_1d<double,2> base_point;
      if (cut_edges(0)==true) // it means it is a cut edge, if it was 0.0 or 1.0 then it would be an uncut edge 
        { 
          base_point[0] = aux_points(0,0);
          base_point[1] = aux_points(0,1);
        }
      else //it means aux_point 0 is a clone of other point, so we go to the second edge.
        { 
          base_point[0] = aux_points(1,0);
          base_point[1] = aux_points(1,1);
        }
      
      for (int i_node = 0; i_node < 3; i_node++)
        {
          double d =    (rPoints(i_node,0) - base_point[0]) * grad_d[0] + (rPoints(i_node,1) - base_point[1]) * grad_d[1] ;
          abs_distance[i_node] = fabs(d);
        }
      
      //assign correct sign to exact distance
      for (int i = 0; i < 3; i++)
        {
          if (rDistances[i] < 0.0)
        {
          exact_distance[i] = -abs_distance[i];
          --most_common_sign;
        }
          else
        {
          exact_distance[i] = abs_distance[i];
          ++most_common_sign;
        }
        }
      
      //compute exact distance gradients
        array_1d<double, 2 > exact_distance_gradient;
        noalias(exact_distance_gradient) = prod(trans(DN_DX), exact_distance);
    
        array_1d<double, 2 > abs_distance_gradient;
        noalias(abs_distance_gradient) = prod(trans(DN_DX), abs_distance);
    
    
    double max_aux_dist_on_cut = -1;
    for (int edge = 0; edge < 3; edge++)
      {
        const int i = edge;
        int j = edge+1;
        if (j==3) j=0;
        if (rDistances[i] * rDistances[j] < 0.0)
          {
        const double tmp = fabs(rDistances[i]) / (fabs(rDistances[i]) + fabs(rDistances[j]));
        //compute the position of the edge node
        double abs_dist_on_cut = abs_distance[i] * tmp + abs_distance[j] * (1.00 - tmp);
                if(abs_dist_on_cut > max_aux_dist_on_cut) max_aux_dist_on_cut = abs_dist_on_cut;
             }
      }
    //we reset all data:
    rGradientsValue[0]=ZeroMatrix(5,2);
    rGradientsValue[1]=ZeroMatrix(5,2);
    rGradientsValue[2]=ZeroMatrix(5,2);

    NEnriched=ZeroMatrix(3,5);
    rGPShapeFunctionValues=ZeroMatrix(3,3);
    
    
    //now we must check the 4 created partitions of the domain. 
    //one has been collapsed, so we discard it and therefore save only one.
    unsigned int partition_number=0;        //  
    //the 3 first partitions are  created using 2 auxiliary nodes and a normal node. at least one of these will be discarded due to zero area       
    //the last one is composed by the 3 auxiliary nodes. it 'looks' wrong, but since at least one has been collapsed, it actually has a normal node.      
    bool found_empty_partition=false;
    
    //the enrichment is directly the shape functions created by using the partition.
    //we have to save for the enrichment 0 and 1 which is the node that will be active, that is, whose shape function value is zero (for some partitions all 3 shape functions are inactive)
    
//     KRATOS_WATCH("hebele");
    

    
    for (unsigned int i=0; i<4; i++) //i partition  
      {
        
        array_1d<int, 2 > active_node_in_enrichment_shape_function; 
        active_node_in_enrichment_shape_function(0)=-1;  active_node_in_enrichment_shape_function(1)=-1; //initialized as if all are inactive -> gradient=0;
        //the same but for the replacement shape functions
        array_1d<int, 3 > active_node_in_replacement_shape_function; 
        active_node_in_replacement_shape_function(0)=-1;  active_node_in_replacement_shape_function(1)=-1; active_node_in_replacement_shape_function(2)=-1; //initialized as if all are inactive -> gradient=0;
        
        int j_aux = i + 2;
        if (j_aux>2) j_aux -= 3; 
        boost::numeric::ublas::bounded_matrix<int,3,2> partition_father_nodes;
        array_1d<double,3> N; bool useful=false;
                  
          int useful_node_for_N0star=-1;
          int useful_node_for_N1star=-1;
        
        if (i<3)
          {
        partition_father_nodes(0,0)=i;
        partition_father_nodes(0,1)=i;
        partition_father_nodes(1,0)=aux_nodes_father_nodes(i,0); //we are using i aux node
        partition_father_nodes(1,1)=aux_nodes_father_nodes(i,1); //we are using i aux node
        partition_father_nodes(2,0)=aux_nodes_father_nodes(j_aux,0); //we are using j_aux node
        partition_father_nodes(2,1)=aux_nodes_father_nodes(j_aux,1); //we are using j_aux node
        
        coord_subdomain(0,0)=rPoints(i,0);
        coord_subdomain(0,1)=rPoints(i,1);
        coord_subdomain(1,0)=aux_points(i,0);
        coord_subdomain(1,1)=aux_points(i,1);
        coord_subdomain(2,0)=aux_points(j_aux,0);
        coord_subdomain(2,1)=aux_points(j_aux,1);
//         if (rPoints(0,0)<5.21527 and rPoints(0,0)>5.21525 and rPoints(1,0)<5.253 and rPoints(1,0)>5.25298 and rPoints(2,0)<5.26242 and rPoints(2,0)>5.26240)
//           KRATOS_WATCH(coord_subdomain);
//         if (rPoints(0,0)<5.26242 and rPoints(0,0)>5.2624 and rPoints(1,0)<5.253 and rPoints(1,0)>5.25298 and rPoints(2,0)<5.31215 and rPoints(2,0)>5.31213)
//           KRATOS_WATCH(coord_subdomain);
//         if (rPoints(0,0)<5.21527 and rPoints(0,0)>5.21525 and rPoints(1,0)<5.204 and rPoints(1,0)>5.2038 and rPoints(2,0)<5.253  and rPoints(2,0)>5.25298)
//           KRATOS_WATCH(coord_subdomain);
//         if (rPoints(0,0)<5.31215 and rPoints(0,0)>5.31213 and rPoints(1,0)<5.253 and rPoints(1,0)>5.25298 and rPoints(2,0)<5.29079 and rPoints(2,0)>5.29077)
//           KRATOS_WATCH(coord_subdomain);
//         if (rPoints(0,0)<5.39205 and rPoints(0,0)>5.39203 and rPoints(1,0)<5.2126 and rPoints(1,0)>5.2124 and rPoints(2,0)<5.31 and rPoints(2,0)>5.29)
//         KRATOS_WATCH(coord_subdomain);
//         if (rPoints(0,0)<5.39205 and rPoints(0,0)>5.39203 and rPoints(1,0)<5.31 and rPoints(1,0)>5.29 and rPoints(2,0)<5.451 and rPoints(2,0)>5.449)
//         KRATOS_WATCH(coord_subdomain);
//         if (rPoints(0,1)<2.0 and rPoints(0,0)>4.99 and rPoints(1,1)<2.0 and rPoints(1,0)>4.99 and rPoints(2,1)<2.0 and rPoints(2,0)>4.99)
//         KRATOS_WATCH(coord_subdomain);

//         if (rPoints(0,0)<5.39205 and rPoints(0,0)>5.39203 and rPoints(1,0)<5.2126 and rPoints(1,0)>5.2124 and rPoints(2,0)<5.31 and rPoints(2,0)>5.29)
//         KRATOS_WATCH(coord_subdomain);
//         if (rPoints(0,0)<5.39205 and rPoints(0,0)>5.39203 and rPoints(1,0)<5.31 and rPoints(1,0)>5.29 and rPoints(2,0)<5.451 and rPoints(2,0)>5.449)
//         KRATOS_WATCH(coord_subdomain);
        
//         KRATOS_WATCH(coord_subdomain);
        
//           if (unsigned_distance0<tolerable_distance or unsigned_distance1<tolerable_distance or unsigned_distance2<tolerable_distance)
//           KRATOS_WATCH(coord_subdomain);
        
        //notice that local nodes 2 and 3 and the possible candidates, with indexes i and j_aux:
        if (aux_node_enrichment_shape_function_index(i)> -1) //that is, local node 2 it is a useful node:
          active_node_in_enrichment_shape_function(  aux_node_enrichment_shape_function_index(i)  )=1;  //saving that local node 2 will be active for either enrichment 1 or 2.
        // else we do nothing, we are not saving this node and the -1 stays
        
        //now the same for the local node 3 (j_aux)
        if (aux_node_enrichment_shape_function_index(j_aux)> -1) //that is, local node 3 it is a useful node:
          active_node_in_enrichment_shape_function( aux_node_enrichment_shape_function_index(j_aux) )=2;    //saving that local node 3 will be active for either enrichment 1 or 2.
        // else we do nothing, we are not saving this node and the -1 stays
                
        active_node_in_replacement_shape_function(i)=0; //standard shape function i will be replaced by the one of local subelement node 1 
        //now local nodes 2 and 3
        if (aux_nodes_father_nodes(i,0)==aux_nodes_father_nodes(i,1))
          active_node_in_replacement_shape_function(aux_nodes_father_nodes(i,0))=1;
        if (aux_nodes_father_nodes(j_aux,0)==aux_nodes_father_nodes(j_aux,1))
          active_node_in_replacement_shape_function(aux_nodes_father_nodes(j_aux,0))=2;
        if( (aux_nodes_father_nodes(i,0)!=aux_nodes_father_nodes(i,1)) && (aux_nodes_father_nodes(j_aux,0)!=aux_nodes_father_nodes(j_aux,1)))
          {
            useful=true;
          }
        
        coord_subdomain(0,0)=rPoints(i,0);
        coord_subdomain(0,1)=rPoints(i,1);
        coord_subdomain(1,0)=aux_points(i,0);
        coord_subdomain(1,1)=aux_points(i,1);
        coord_subdomain(2,0)=aux_points(j_aux,0);
        coord_subdomain(2,1)=aux_points(j_aux,1);
//           if (rPoints(0,0)<5.21527 and rPoints(0,0)>5.21525 and rPoints(1,0)<5.253 and rPoints(1,0)>5.25298 and rPoints(2,0)<5.26242 and rPoints(2,0)>5.26240)
//           KRATOS_WATCH(coord_subdomain);
//         if (rPoints(0,0)<5.26242 and rPoints(0,0)>5.2624 and rPoints(1,0)<5.253 and rPoints(1,0)>5.25298 and rPoints(2,0)<5.31215 and rPoints(2,0)>5.31213)
//           KRATOS_WATCH(coord_subdomain);
//         if (rPoints(0,0)<5.21527 and rPoints(0,0)>5.21525 and rPoints(1,0)<5.204 and rPoints(1,0)>5.2038 and rPoints(2,0)<5.253  and rPoints(2,0)>5.25298)
//           KRATOS_WATCH(coord_subdomain);
//         if (rPoints(0,0)<5.31215 and rPoints(0,0)>5.31213 and rPoints(1,0)<5.253 and rPoints(1,0)>5.25298 and rPoints(2,0)<5.29079 and rPoints(2,0)>5.29077)
//           KRATOS_WATCH(coord_subdomain);
//         if (rPoints(0,0)<5.39205 and rPoints(0,0)>5.39203 and rPoints(1,0)<5.2126 and rPoints(1,0)>5.2124 and rPoints(2,0)<5.31 and rPoints(2,0)>5.29)
//         KRATOS_WATCH(coord_subdomain);
//         if (rPoints(0,0)<5.39205 and rPoints(0,0)>5.39203 and rPoints(1,0)<5.31 and rPoints(1,0)>5.29 and rPoints(2,0)<5.451 and rPoints(2,0)>5.449)
//         KRATOS_WATCH(coord_subdomain);
//         if (rPoints(0,1)<2.0 and rPoints(0,0)>4.99 and rPoints(1,1)<2.0 and rPoints(1,0)>4.99 and rPoints(2,1)<2.0 and rPoints(2,0)>4.99)
//         KRATOS_WATCH(coord_subdomain);
        
//       if (rPoints(0,0)<5.313 and rPoints(0,0)>5.2038 and rPoints(0,1)<1.08 and rPoints(0,1)>0.996 
//           and rPoints(1,0)<5.313 and rPoints(1,0)>5.2038 and rPoints(1,1)<1.08 and rPoints(1,1)>0.996
//           and rPoints(2,0)<5.313  and rPoints(2,0)>5.2038 and rPoints(2,1)<1.08 and rPoints(2,1)>0.996)
//       {
//          KRATOS_WATCH(coord_subdomain);
//         }
//         if (rPoints(0,0)<5.39205 and rPoints(0,0)>5.39203 and rPoints(1,0)<5.2126 and rPoints(1,0)>5.2124 and rPoints(2,0)<5.31 and rPoints(2,0)>5.29)
//         KRATOS_WATCH(coord_subdomain);
//         if (rPoints(0,0)<5.39205 and rPoints(0,0)>5.39203 and rPoints(1,0)<5.31 and rPoints(1,0)>5.29 and rPoints(2,0)<5.451 and rPoints(2,0)>5.449)
//         KRATOS_WATCH(coord_subdomain);
        
//         KRATOS_WATCH(coord_subdomain);
        

        
//         if (unsigned_distance0<tolerable_distance or unsigned_distance1<tolerable_distance or unsigned_distance2<tolerable_distance)
//           KRATOS_WATCH(coord_subdomain);
        //an aux_node is useful when one of its father nodes is the real node of the subelement. that means that the edge is part of the subelement.
        if(partition_father_nodes(1,0)==i || partition_father_nodes(1,1)==i) //if one of the father nodes of node_aux_i is equal to the real node i
        {
          if(aux_node_enrichment_shape_function_index(i)==0)
            useful_node_for_N0star=1;
          if(aux_node_enrichment_shape_function_index(i)==1)
            useful_node_for_N1star=1;
        }
        if(partition_father_nodes(2,0)==j_aux || partition_father_nodes(2,1)==j_aux) //if one of the father nodes of node_aux_i is equal to the real node i
        {
          if(aux_node_enrichment_shape_function_index(j_aux)==0)
            useful_node_for_N0star=2;
          if(aux_node_enrichment_shape_function_index(j_aux)==1)
            useful_node_for_N1star=2;
        }
          }
        else
          {
        partition_father_nodes=aux_nodes_father_nodes;
        coord_subdomain=aux_points;
//           if (rPoints(0,0)<5.21527 and rPoints(0,0)>5.21525 and rPoints(1,0)<5.253 and rPoints(1,0)>5.25298 and rPoints(2,0)<5.26242 and rPoints(2,0)>5.26240)
//           KRATOS_WATCH(coord_subdomain);
//         if (rPoints(0,0)<5.26242 and rPoints(0,0)>5.2624 and rPoints(1,0)<5.253 and rPoints(1,0)>5.25298 and rPoints(2,0)<5.31215 and rPoints(2,0)>5.31213)
//           KRATOS_WATCH(coord_subdomain);
//         if (rPoints(0,0)<5.21527 and rPoints(0,0)>5.21525 and rPoints(1,0)<5.204 and rPoints(1,0)>5.2038 and rPoints(2,0)<5.253  and rPoints(2,0)>5.25298)
//           KRATOS_WATCH(coord_subdomain);
//         if (rPoints(0,0)<5.31215 and rPoints(0,0)>5.31213 and rPoints(1,0)<5.253 and rPoints(1,0)>5.25298 and rPoints(2,0)<5.29079 and rPoints(2,0)>5.29077)
//           KRATOS_WATCH(coord_subdomain);
        
//         if (rPoints(0,0)<5.39205 and rPoints(0,0)>5.39203 and rPoints(1,0)<5.2126 and rPoints(1,0)>5.2124 and rPoints(2,0)<5.31 and rPoints(2,0)>5.29)
//         KRATOS_WATCH(coord_subdomain);
//         if (rPoints(0,0)<5.39205 and rPoints(0,0)>5.39203 and rPoints(1,0)<5.31 and rPoints(1,0)>5.29 and rPoints(2,0)<5.451 and rPoints(2,0)>5.449)
//         KRATOS_WATCH(coord_subdomain);
//         if (rPoints(0,1)<2.0 and rPoints(0,0)>4.99 and rPoints(1,1)<2.0 and rPoints(1,0)>4.99 and rPoints(2,1)<2.0 and rPoints(2,0)>4.99)
//         KRATOS_WATCH(coord_subdomain);
        

        
//       if (rPoints(0,0)<5.313 and rPoints(0,0)>5.2038 and rPoints(0,1)<1.08 and rPoints(0,1)>0.996 
//           and rPoints(1,0)<5.313 and rPoints(1,0)>5.2038 and rPoints(1,1)<1.08 and rPoints(1,1)>0.996
//           and rPoints(2,0)<5.313  and rPoints(2,0)>5.2038 and rPoints(2,1)<1.08 and rPoints(2,1)>0.996)
//       {
//          KRATOS_WATCH(coord_subdomain);
//         }
//         if (rPoints(0,0)<5.39205 and rPoints(0,0)>5.39203 and rPoints(1,0)<5.2126 and rPoints(1,0)>5.2124 and rPoints(2,0)<5.31 and rPoints(2,0)>5.29)
//         KRATOS_WATCH(coord_subdomain);
//         if (rPoints(0,0)<5.39205 and rPoints(0,0)>5.39203 and rPoints(1,0)<5.31 and rPoints(1,0)>5.29 and rPoints(2,0)<5.451 and rPoints(2,0)>5.449)
//         KRATOS_WATCH(coord_subdomain);
//        
//         if (unsigned_distance0<tolerable_distance or unsigned_distance1<tolerable_distance or unsigned_distance2<tolerable_distance)
//           KRATOS_WATCH(coord_subdomain);
        
//         KRATOS_WATCH(coord_subdomain);
        
        
        
        
        useful=true;
        int non_aux_node=-1; //one of the nodes of this partition is actually a real node, which has repeated father node 
        int non_aux_node_father_node=-1; //the real node id
        //we have to check the 3 of them: 
        for (int j = 0; j < 3; j++)
            {
              if(partition_father_nodes(j,0)==partition_father_nodes(j,1))
                {
              non_aux_node=j;
              non_aux_node_father_node=partition_father_nodes(j,0);
            }
            }
        for (int j = 0; j < 3; j++)
          {
            if (aux_node_enrichment_shape_function_index(j)> -1) //that is, local node j it is a useful node:
            {
              active_node_in_enrichment_shape_function(  aux_node_enrichment_shape_function_index(j)  ) = j;

              if(partition_father_nodes(j,0)==non_aux_node_father_node || partition_father_nodes(j,1)==non_aux_node_father_node)
            {
              if (aux_node_enrichment_shape_function_index(j)==0)
                useful_node_for_N0star=j;
              if (aux_node_enrichment_shape_function_index(j)==1)
                useful_node_for_N1star=j;
            }
            }
            //to replace the standard shape functions:
            if (aux_nodes_father_nodes(j,0)==aux_nodes_father_nodes(j,1))
              active_node_in_replacement_shape_function(aux_nodes_father_nodes(j,0))=j; 
          }
        //found_last_partition=true;
          }
        //calculate data of this partition
        double temp_area;
        CalculateGeometryData(coord_subdomain, DN_DX_subdomain, temp_area);
        if (temp_area > 1.0e-20) //ok, it does not have zero area
          {
        //KRATOS_ERROR(std::logic_error, "method not implemented", "");
        rVolumes(partition_number)=temp_area;
        //we look for the gauss point of the partition:
        double x_GP_partition =  one_third * ( coord_subdomain(0,0) + coord_subdomain(1,0) + coord_subdomain(2,0) );
        double y_GP_partition =  one_third * ( coord_subdomain(0,1) + coord_subdomain(1,1) + coord_subdomain(2,1) );
        double z_GP_partition  = 0.0;
        //we reset the coord_subdomain matrix so that we have the whole element again:
        coord_subdomain = rPoints;  
        //and we calculate its shape function values
        CalculatePosition ( coord_subdomain , x_GP_partition ,y_GP_partition ,z_GP_partition , N);
        //we check the partition sign.
        const double partition_sign = (N(0)*rDistances(0) + N(1)*rDistances(1) + N(2)*rDistances(2))/fabs(N(0)*rDistances(0) + N(1)*rDistances(1) + N(2)*rDistances(2));
        //rPartitionsSign(partition_number)=partition_sign;
        
        rGPShapeFunctionValues(partition_number,0)=N(0);
        rGPShapeFunctionValues(partition_number,1)=N(1);
        rGPShapeFunctionValues(partition_number,2)=N(2);
        
        
        //compute enriched shape function values
        double dist = 0.0;
        double abs_dist = 0.0;
        for (int j = 0; j < 3; j++)
          {
            dist += rGPShapeFunctionValues(partition_number, j) * exact_distance[j];
            abs_dist += rGPShapeFunctionValues(partition_number, j) * abs_distance[j];
          }
        
            if (partition_sign < 0.0)
              rPartitionsSign[partition_number] = -1.0;
            else
              rPartitionsSign[partition_number] = 1.0;
            
            //We use the sublement shape functions and derivatives:
            //we loop the 2 enrichment shape functions:
            for (int index_shape_function = 0; index_shape_function < 2; index_shape_function++) //enrichment shape function
              {
            if (active_node_in_enrichment_shape_function(index_shape_function) > -1) //recall some of them are inactive:
              {
                NEnriched(partition_number, index_shape_function+3 ) = one_third ; //only one gauss point. if more were to be used, it would still be simple (1/2,1/2,0);(0,1/2,1/2);(1/2,0,1/2);
                for (int j = 0; j < 2; j++) //x,y,(z)
                  {
                rGradientsValue[partition_number](index_shape_function+3, j) = DN_DX_subdomain(active_node_in_enrichment_shape_function(index_shape_function),j); 
                  }
                
              }
            else
              {
                NEnriched(partition_number, index_shape_function+3 ) = 0.0;
                for (int j = 0; j < 2; j++) //x,y,(z)
                  {
                rGradientsValue[partition_number](index_shape_function+3, j) = 0.0;
                  }
              }
            
              }

            array_1d<int,(2+1)> replacement_shape_function_nodes = ZeroVector(3);
            for (int index_shape_function = 0; index_shape_function < 3; index_shape_function++) //replacement shape function
            {
                int active_node=-1;
                for (int j = 0; j < 3; j++) //number_of_nodes in the subelement
                    if (partition_father_nodes(j,0)==index_shape_function && partition_father_nodes(j,1)==index_shape_function)
                    {
                        active_node=j;
                        break;
                    }
                if(active_node> -1)
                {
                    for (int j = 0; j < 2; j++) //x,y,(z)
                        rGradientsValue[partition_number](index_shape_function, j) = DN_DX_subdomain(active_node,j);
                     NEnriched(partition_number, index_shape_function ) = one_third;
                     replacement_shape_function_nodes(index_shape_function) = active_node;
                }
                else
                {
                  for (int j = 0; j < 2; j++) //x,y,(z)
                    rGradientsValue[partition_number](index_shape_function, j) = 0.0;
                  replacement_shape_function_nodes(index_shape_function) = -1;
                }
            }               
            

            
            //We use the sublement shape functions and derivatives:
            //we loop the 3 replacement shape functions:
            unsigned int number_of_real_nodes=0; //(each partition can have either 1 or 2);
            for (int index_shape_function = 0; index_shape_function < 3; index_shape_function++) //enrichment shape function
              {
            if (active_node_in_replacement_shape_function(index_shape_function) > -1) //recall some of them are inactive:
              number_of_real_nodes++;
              }
            
            if(useful_node_for_N0star > -1)
              { 
            for (int j = 0; j < 2; j++) //x,y,(z)
              { 
                if(partition_sign>0)
                  { 
                //first two rows are for the side where N*1 = 1
                //row 1 is for gradN*1
                rGradientpositive(3, j)=DN_DX_subdomain(useful_node_for_N0star, j);
                //row 2 is for gradN*2
                if (active_node_in_enrichment_shape_function(1) > -1)
                  rGradientpositive(4, j)=DN_DX_subdomain(active_node_in_enrichment_shape_function(1), j);
                
                for (int index_shape_function = 0; index_shape_function < 3; index_shape_function++) //replacement shape function
                  {
                    if(replacement_shape_function_nodes(index_shape_function)>-1)
                      {
                    rGradientpositive(index_shape_function, j) = DN_DX_subdomain(replacement_shape_function_nodes(index_shape_function), j);
                      }
                  }
                  }
                else
                  {
                rGradientnegative(3, j) =DN_DX_subdomain(useful_node_for_N0star, j);
                if (active_node_in_enrichment_shape_function(1) > -1)
                  rGradientnegative(4, j) =DN_DX_subdomain(active_node_in_enrichment_shape_function(1), j);
                
                for (int index_shape_function = 0; index_shape_function < 3; index_shape_function++) //replacement shape function
                  {
                    if(replacement_shape_function_nodes(index_shape_function)>-1)
                      {
                    rGradientnegative(index_shape_function, j) = DN_DX_subdomain(replacement_shape_function_nodes(index_shape_function), j);
                      }
                  }
                  } 
              } 
              }
            if(useful_node_for_N1star > -1)
              { 
            for (int j = 0; j < 2; j++) //x,y,(z)
              { 
                if(partition_sign>0)
                  { 
                //rows 3 and 4 are for the side where N*2 = 1
                //row 4 is for gradN*2
                rGradientpositive(9, j)=DN_DX_subdomain(useful_node_for_N1star, j);
                //row 3 is for gradN*1
                if (active_node_in_enrichment_shape_function(0) > -1)
                  rGradientpositive(8, j)=DN_DX_subdomain(active_node_in_enrichment_shape_function(0), j);
                
                for (int index_shape_function = 0; index_shape_function < 3; index_shape_function++) //replacement shape function
                  {
                    if(replacement_shape_function_nodes(index_shape_function)>-1)
                      {
                    rGradientpositive(index_shape_function+5, j) = DN_DX_subdomain(replacement_shape_function_nodes(index_shape_function), j);
                      }
                  }
                  }
                else
                  {
                rGradientnegative(9, j)=DN_DX_subdomain(useful_node_for_N1star, j);
                if(active_node_in_enrichment_shape_function(0) > -1)
                  rGradientnegative(8, j)=DN_DX_subdomain(active_node_in_enrichment_shape_function(0), j);
                
                        for (int index_shape_function = 0; index_shape_function < 3; index_shape_function++) //replacement shape function
                          {
                            if(replacement_shape_function_nodes(index_shape_function)>-1)
                            {
                                rGradientnegative(index_shape_function+5, j) = DN_DX_subdomain(replacement_shape_function_nodes(index_shape_function), j);
                            }
                          }
                          } 
                       }    
                }
    
            partition_number++;
                      
            
          }
        else
          found_empty_partition=true;
          }
        if (found_empty_partition==false)
          KRATOS_WATCH("WROOOONGGGGGGGGGGG");
        
        //KRATOS_WATCH(NEnriched)
        return 3;
        KRATOS_CATCH("");
        
      } 
		
		

		

		//template<class TMatrixType, class TVectorType, class TGradientType>
        static inline double CalculateVolume2D(
            const bounded_matrix<double, 3, 3 > & coordinates)
        {
            double x10 = coordinates(1,0) - coordinates(0,0);
            double y10 = coordinates(1,1) - coordinates(0,1);

            double x20 = coordinates(2,0) - coordinates(0,0);
            double y20 = coordinates(2,1) - coordinates (0,1);
            double detJ = x10 * y20-y10 * x20;
            return 0.5*detJ;
        }
        
        static inline bool CalculatePosition(const bounded_matrix<double, 3, 3 > & coordinates,
                const double xc, const double yc, const double zc,
                array_1d<double, 3 > & N
                )
        {
            double x0 = coordinates(0,0);
            double y0 = coordinates(0,1);
            double x1 = coordinates(1,0);
            double y1 = coordinates(1,1);
            double x2 = coordinates(2,0);
            double y2 = coordinates(2,1);

            double area = CalculateVol(x0, y0, x1, y1, x2, y2);
            double inv_area = 0.0;
            if (area == 0.0)
            {
                KRATOS_THROW_ERROR(std::logic_error, "element with zero area found", "");
            } else
            {
                inv_area = 1.0 / area;
            }


            N[0] = CalculateVol(x1, y1, x2, y2, xc, yc) * inv_area;
            N[1] = CalculateVol(x2, y2, x0, y0, xc, yc) * inv_area;
            N[2] = CalculateVol(x0, y0, x1, y1, xc, yc) * inv_area;
            //KRATOS_WATCH(N);

            if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0) //if the xc yc is inside the triangle return true
                return true;

            return false;
        }
        
    //3d    
    static inline void CalculateGeometryData(
        boost::numeric::ublas::bounded_matrix<double,(3+1), 3 >& coordinates,
        boost::numeric::ublas::bounded_matrix<double,4,3>& DN_DX,
        double& Volume)
    {
        
        double x10 = coordinates(1,0) - coordinates(0,0);
        double y10 = coordinates(1,1) - coordinates(0,1);
        double z10 = coordinates(1,2) - coordinates(0,2);

        double x20 = coordinates(2,0) - coordinates(0,0);
        double y20 = coordinates(2,1) - coordinates (0,1);
        double z20 = coordinates(2,2) - coordinates (0,2);
        
        double x30 = coordinates(3,0) - coordinates(0,0);
        double y30 = coordinates(3,1) - coordinates(0,1);
        double z30 = coordinates(3,2) - coordinates (0,2);

        double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;

        DN_DX(0,0) = -y20 * z30 + y30 * z20 + y10 * z30 - z10 * y30 - y10 * z20 + z10 * y20;
        DN_DX(0,1) = -z20 * x30 + x20 * z30 - x10 * z30 + z10 * x30 + x10 * z20 - z10 * x20;
        DN_DX(0,2) = -x20 * y30 + y20 * x30 + x10 * y30 - y10 * x30 - x10 * y20 + y10 * x20;
        DN_DX(1,0) = y20 * z30 - y30 * z20;
        DN_DX(1,1) = z20 * x30 - x20 * z30;
        DN_DX(1,2) = x20 * y30 - y20 * x30;
        DN_DX(2,0) = -y10 * z30 + z10 * y30;
        DN_DX(2,1) = x10 * z30 - z10 * x30;
        DN_DX(2,2) = -x10 * y30 + y10 * x30;
        DN_DX(3,0) = y10 * z20 - z10 * y20;
        DN_DX(3,1) = -x10 * z20 + z10 * x20;
        DN_DX(3,2) = x10 * y20 - y10 * x20;

        DN_DX /= detJ;
        


        Volume = detJ*0.1666666666666666666667;
    }
        static inline double CalculateVol(const double x0, const double y0,
                const double x1, const double y1,
                const double x2, const double y2
                )
        {
            return 0.5 * ((x1 - x0)*(y2 - y0)- (y1 - y0)*(x2 - x0));
        }
        
            static inline void CalculateGeometryData(
            const bounded_matrix<double, 3, 3 > & coordinates,
            boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,
            double& Area)
        {
            double x10 = coordinates(1,0) - coordinates(0,0);
            double y10 = coordinates(1,1) - coordinates(0,1);

            double x20 = coordinates(2,0) - coordinates(0,0);
            double y20 = coordinates(2,1) - coordinates (0,1);

            //Jacobian is calculated:
            //  |dx/dxi  dx/deta|   |x1-x0   x2-x0|
            //J=|               |=  |             |
            //  |dy/dxi  dy/deta|   |y1-y0   y2-y0|


            double detJ = x10 * y20-y10 * x20;

            DN_DX(0,0) = -y20 + y10;
            DN_DX(0,1) = x20 - x10;
            DN_DX(1,0) =  y20      ;
            DN_DX(1,1) = -x20     ;
            DN_DX(2,0) = -y10      ;
            DN_DX(2,1) = x10       ;

            DN_DX /= detJ;

            Area = 0.5*detJ;
        }
        
        
            

		

		
	protected:


	private:
	
	
	void Check()
	{
		if(mr_model_part.NodesBegin()->SolutionStepsDataHas(DISTANCE) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing DISTANCE variable on solution step data","");
		if(mr_model_part.NodesBegin()->SolutionStepsDataHas(PRESS_PROJ) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESS_PROJ variable on solution step data","");
		if(mr_model_part.NodesBegin()->SolutionStepsDataHas(VELOCITY) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing VELOCITY variable on solution step data","");
		if(mr_model_part.NodesBegin()->SolutionStepsDataHas(PRESSURE) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing PRESSURE variable on solution step data","");
		if(mr_model_part.NodesBegin()->SolutionStepsDataHas(PROJECTED_VELOCITY) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing PROJECTED_VELOCITY variable on solution step data","");
		if(mr_model_part.NodesBegin()->SolutionStepsDataHas(DELTA_VELOCITY) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing DELTA_VELOCITY variable on solution step data","");
		if(mr_model_part.NodesBegin()->SolutionStepsDataHas(MESH_VELOCITY) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing MESH_VELOCITY variable on solution step data","");
		if(mr_model_part.NodesBegin()->SolutionStepsDataHas(YP) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing YP variable on solution step data","");
		if(mr_model_part.NodesBegin()->SolutionStepsDataHas(NORMAL) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing NORMAL variable on solution step data","");
		if(mr_model_part.NodesBegin()->SolutionStepsDataHas(NODAL_AREA) == false) KRATOS_THROW_ERROR(std::invalid_argument,"missing NODAL_AREA variable on solution step data","");
	}
	


		
	
		

		

		
	
	ModelPart& mr_model_part;
    ModelPart& mr_model_part2;
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

#endif // KRATOS_NITSCHE_DECONDENSATION_UTILITY_INCLUDED  defined 


