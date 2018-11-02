//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Sonja Schneider
//

// From S. Schneider Implementation of a Chimera Technique TUM Master Thesis 2015

/****************************************************************************************
 * Routine Name: interpolation_utility.h
 * Content:
 * In this routine there exist 4 main-functions:
 * - DirichletBoundaryCalculation
 * - NeumannBoundaryCalculation
 * - DirichletBoundaryCalculationKreis
 * - NeumannBoundaryCalculationRechteck
 * Both functions calculate the desired values (either the Temperature on the boundary
 * of the domain or the flux (Neumann condition)) and interpolate them then to the nodes
 * onto the opposite domain segment.
 * *************************************************************************************/

/****************************************************************************************
 * Preprocessor-Directives
 * *************************************************************************************/
#if !defined KRATOS_INTERPOLATION_UTILITY_H
#define KRATOS_INTERPOLATION_UTILITY_H

// System includes
#include <math.h>
#include <string>
#include <iostream>
#include <stdlib.h>
#include "includes/define.h"
#include "includes/model_part.h"

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "geometries/triangle_2d_3.h"
#include "utilities/timer.h"
#include "includes/kratos_flags.h"
#include "utilities/binbased_fast_point_locator.h"

#include "chimera_application_variables.h"

//Database includes
#include "spatial_containers/spatial_containers.h"


namespace Kratos
{
    // Class for Calculating the Distance
    template< class T, std::size_t dim >
    class DistanceCalculator
    {
    public:
        double operator()( T const& p1, T const& p2 )
        {
            double dist = 0.0;
            for( std::size_t i = 0 ; i < dim ; i++)
            {
                double tmp = p1[i] - p2[i];
                dist += tmp*tmp;
            }
            return dist;
            //square distance because it is easier to work without the square root
        }
    };

  ///@addtogroup PureDiffusionApplication
  ///@{

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
/****************************************************************************************
 * Class InterpolationUtility
 * Annotation:
 * This class contains the four main-calculation functions:
 * - DirichletBoundaryCalculation
 * - NeumannBoundaryCalculation
 * - DirichletBoundaryCalculationKreis
 * - NeumannBoundaryCalculationRechteck
 * *************************************************************************************/
  template<std::size_t TDim >
  class InterpolationUtility
  {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of InterpolationUtility
      KRATOS_CLASS_POINTER_DEFINITION(InterpolationUtility<TDim >);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      InterpolationUtility()
      {
          //KRATOS_TRY

          //KRATOS_CATCH("")
      }


      /// Destructor.
      virtual ~InterpolationUtility() {}

      ///@}
      ///@name Operators
      ///@{

      ///@}
      ///@name Operations
      ///@{
/****************************************************************************************
 * Here all the implementation takes part
 * *************************************************************************************/      // Just a test function




  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  * FUNCTION DESCRIPTION: CreateBoundaryFaces
  * This function establishes faces between two nodes on the boundary nodes after the
  * hole elements of the background mesh are deleted.
  * +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
      void CreateBoundaryFaces(ModelPart& background, const Condition& rReferenceCondition )
      {
          //****************************************************************
          // loop over model part, store pairs of nodes on the boundary    *
          //****************************************************************

          // Create an Array to store the border Ids
          std::vector< Node<3>::Pointer > list_of_element_nodes;
          std::size_t nr_of_elements = background.NumberOfElements();
          list_of_element_nodes.reserve(nr_of_elements);
          std::vector< Condition::NodesArrayType > list_of_node_pairs; //????
          list_of_node_pairs.reserve(nr_of_elements);

          //l for filling list_of_node_pairs
          int l=0;

          // creating the list
          for(ModelPart::ElementsContainerType::iterator it = background.ElementsBegin(); it!=background.ElementsEnd(); it++)
          {
              //get the connectivity of the element into geom
              Geometry< Node<3> >& geom = it->GetGeometry();
              //Setting up counter
              int counter =0;

              for (int i = 0; i < 3; i++ )
              {
                  if(geom[i].FastGetSolutionStepValue(FLAG_VARIABLE)==1.0)
                  {
                      counter = counter +1;
                      //list_of_new_nodes.push_back( geom(i));
                  }
              }

              if (counter == 2)
              {
                  list_of_node_pairs.push_back( Condition::NodesArrayType() );
                  for(int i=0; i<3; i++)
                  {  //int j=0;
                      if(geom[i].FastGetSolutionStepValue(FLAG_VARIABLE)==1.0)
                      {
                          list_of_node_pairs[l].push_back( geom(i) );
                          //j=j+1;
                      }
                  }
                  l=l+1;
              }
          }

          int size_of_list = list_of_node_pairs.size();

          /**********************************************
          * create conditions using the list of pairs   *
          **********************************************/
          /* add conditions to ModelPart */
          //std::string ConditionName("ChimeraFluidCouplingCondition2D");
          //std::string ConditionName("ChimeraThermalCouplingCondition2D");
          //const Condition& rReferenceCondition = KratosComponents<Condition>::Get(ConditionName);
          Properties::Pointer prop = background.ConditionsBegin()->pGetProperties();

          std::size_t ConditonId = background.NumberOfConditions();
          for (int i=0 ; i< size_of_list ;i++)
          {
              Condition::Pointer pCondition = rReferenceCondition.Create(++ConditonId,list_of_node_pairs[i],prop);
              background.pConditions()->push_back(pCondition);
          }
          background.Conditions().Sort();

      }

      void Test()
        {}




      /****************************************************************************************
       * Function: Projection_Flux
       * Description:
       * First the boundary nodes of the rechteck model part are figured out and
       * saved in a list.
       *
       * With the interpolation function the interpolation of a Scalar Variable of one
       * modelpart to the other one is done...
       *
       * Variable Explanation:
       * rOrigin_ModelPart: the model part  all the variable should be taken from
       * rDestination_ModelPart: the destination model part where we want to know the values of
       *                         the variables if you want to pass only one variable
       * rOriginVariable: the name of the interpolated variable in the origin model part
       * rOriginVariable: the name of the interpolated variable in the destination model part
       * *************************************************************************************/

            template<class TDataType>

            void Projection_Flux(
                ModelPart& rOrigin_ModelPart ,
                ModelPart& rDestination_ModelPart,
                Variable<TDataType>& rOriginVariable ,
                Variable<TDataType>& rDestinationVariable,
                int counter,
                const int domain_size
                // Variable for preallocating the array 'list_of_new_nodes'
            )
            {

                KRATOS_TRY

              /**********************************************************************
               * Creating an array 'list_of_new_nodes' which contains all the nodes
               * on which an interpolation should take part
               * *******************************************************************/

                  // Create an Array to store the border Ids
                  std::vector< Node<3>::Pointer > list_of_new_nodes;
                  list_of_new_nodes.reserve(counter);

                  // creating the list
                  for(ModelPart::ElementsContainerType::iterator it = rDestination_ModelPart.ElementsBegin(); it!=rDestination_ModelPart.ElementsEnd(); it++)
                  {
                      //get the connectivity of the element
                      Geometry< Node<3> >& geom = it->GetGeometry();
/*
                      //loop on all the nodes in the geometry
                      //and get the neighbour elements to each of the nodes in the geometry
                      for(Geometry< Node<3> >::iterator kkk = geom.begin(); kkk!=geom.end(); kkk++)
                      {
                           if(kkk->FastGetSolutionStepValue(BOUNDARY_NODE)==1.0)
                               {
                                   list_of_new_nodes.push_back( *(kkk.base()));
                               }

                      }

  */

                      int size = geom.size();

                      for (int i = 0; i < size; i++ )
                      {
                          if(geom[i].FastGetSolutionStepValue(BOUNDARY_NODE)==1.0)
                                                     {

                                   list_of_new_nodes.push_back( geom(i));
                               }

                      }



                  }

                //*******************************************************************
                // properties to be used in the generation
                //*******************************************************************
                Properties::Pointer properties = rDestination_ModelPart.GetMesh().pGetProperties(1);

                //defintions for spatial search
                typedef Node<3> PointType;
                typedef Node<3>::Pointer PointTypePointer;
                typedef std::vector<PointType::Pointer>           PointVector;
                typedef std::vector<PointType::Pointer>::iterator PointIterator;
                typedef std::vector<double>               DistanceVector;
                typedef std::vector<double>::iterator     DistanceIterator;

                // ESTABLISHING SEARCH STRUCTURE
                //------------------------------
                //starting calculating time of construction of the kdtree
                boost::timer kdtree_construction;

                //*************
                // Bucket types
                typedef Bucket< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;

                //*************
                // DynamicBins;

                typedef Tree< KDTreePartition<BucketType> > tree; 		//Kdtree;

                KRATOS_INFO( "kdt construction time ") << kdtree_construction.elapsed() << std::endl;
                //finishing calculating time of construction of the kdtree

                //create a spatial database with the list of new nodes
                std::size_t bucket_size = 20;
                tree nodes_tree(list_of_new_nodes.begin(),list_of_new_nodes.end(),bucket_size);

                //work arrays
                Node<3> work_point(0,0.0,0.0,0.0);
                std::size_t MaximumNumberOfResults = 10000;
                PointVector Results(MaximumNumberOfResults);
                DistanceVector ResultsDistances(MaximumNumberOfResults);

                array_1d<double,TDim+1> N; //Shape functions vector//

                boost::timer search_and_interpolation_time;
                //loop over all of the elements in the "old" list to perform the interpolation
                for( ModelPart::ElementsContainerType::iterator el_it = rOrigin_ModelPart.ElementsBegin();
                        el_it != rOrigin_ModelPart.ElementsEnd(); el_it++)
                {
                    Geometry<Node<3> >&geom = el_it->GetGeometry();

                    //find the center and "radius" of the element
                    double xc, yc, zc,  radius;
                    CalculateCenterAndSearchRadius( geom, xc,yc,zc, radius, N);
                    work_point.X() = xc;
                    work_point.Y() = yc;
                    work_point.Z() = zc;

                    //find all of the new nodes within the radius
                    int number_of_points_in_radius;

                    //look between the new nodes which of them is inside the radius of the circumscribed circle
                    number_of_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, Results.begin(),
                                                 ResultsDistances.begin(),  MaximumNumberOfResults);


                    //check if inside
                    for( PointIterator it_found = Results.begin(); it_found != Results.begin() + number_of_points_in_radius; it_found++)
                    {
                        bool is_inside = false;
                        //once we are sure the node in inside the circle we have to see if it is inside the triangle i.e. if all the Origin element shape functions are >1

                        //double is_visited = (*it_found)->GetValue(IS_VISITED);

                        // CALCULATE THE FLUXES
                        // Vector of shape function derivatives
                        //GeometryData::ShapeFunctionsGradientsType N_der; //Shape functions vector//
                        is_inside = CalculatePosition(geom,	(*it_found)->X(),(*it_found)->Y(),(*it_found)->Z(),N);

                        //if the node falls inside the element interpolate
                        if(is_inside == true)// && is_visited != 1.0)
                        {
                            //Interpolating all the rOriginVariables to get their nodal value in the rDestination_ModelPart
                            Calculate_Flux( el_it,  N, *it_found , rOriginVariable , rDestinationVariable ,domain_size );
                            //pnode->GetValue(IS_VISITED) = 1.0;
                        }
                    }

                }

/*
                // Output of the FACE_HEAT_FLUXES
                for(
                     ModelPart::NodesContainerType::iterator inode = rDestination_ModelPart.NodesBegin();
                     inode != rDestination_ModelPart.NodesEnd();
                     inode++)
                 {
                    if(inode->X() == 0.4)
                    {
                        KRATOS_WATCH(inode->FastGetSolutionStepValue(FACE_HEAT_FLUX));
                        KRATOS_WATCH(inode->Id());
                    }

                 }
*/
                KRATOS_INFO( "search and interpolation time ") << search_and_interpolation_time.elapsed() << std::endl;
                KRATOS_CATCH("")
            }

            /**********************************************************************************
             * Projection Flux Gauss-Points
             * *******************************************************************************/

            void Projection_Flux_Gauss(
                ModelPart& rOrigin_ModelPart ,
                ModelPart& rDestination_ModelPart,
                //Variable<TDataType>& rOriginVariable ,
                //Variable<TDataType>& rDestinationVariable,
                Variable<double>& rOriginVariable,
                Variable<double>& flux,
                int counter,
                double Omega = 1.0
                // Variable for preallocating the array 'list_of_new_nodes'
            )
            {

                KRATOS_TRY

                // CREATE SEARCH STRUCTURE CONTAINING THE ORIGIN MESH
                BinBasedFastPointLocator<TDim> SearchUtil = BinBasedFastPointLocator<TDim>(rOrigin_ModelPart);
                SearchUtil.UpdateSearchDatabase();

                ProcessInfo& rProcessInfo = rDestination_ModelPart.GetProcessInfo();

                // Mark the boundary faces with a flag interface
                // Get the connectivity of the face
                std::size_t FaceTotal = 0;
                for(ModelPart::ConditionIterator itCond = rDestination_ModelPart.ConditionsBegin(); itCond != rDestination_ModelPart.ConditionsEnd(); itCond++)
                {
                    //getting the connectivity of the face
                    Geometry< Node<3> >& rGeom = itCond->GetGeometry();
                    int counter = 0;
                    for(Geometry< Node<3> >::iterator node_it = rGeom.begin(); node_it!=rGeom.end(); node_it++)
                    {
                        if(node_it->FastGetSolutionStepValue(FLAG_VARIABLE)==1.0)
                        {
                            counter = counter +1;
                        }

                    }
                    if(counter>1)
                    {
                        itCond->Set(INTERFACE,true);
                        FaceTotal++;
                    }
                }

                if (FaceTotal == 0)
                {
                    KRATOS_THROW_ERROR(std::runtime_error,"No destination nodes found for interpolation. Assign FLAG_VARIABLE=1 to destination nodes","");
                }


                // Loop over conditions in destination mesh, obtain gauss points
                for (ModelPart::ConditionIterator itCond = rDestination_ModelPart.ConditionsBegin(); itCond != rDestination_ModelPart.ConditionsEnd(); itCond++)
                {
                    // filter conditions, use only those on the boundary
                    if (itCond->Is(INTERFACE))
                    {
                        // Obtain integration points on the condition
                        Geometry< Node<3> >& rGeom = itCond->GetGeometry();
                        const GeometryData::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(itCond->GetIntegrationMethod());
                        const std::size_t NumGauss = IntegrationPoints.size();
                        Matrix NCond = rGeom.ShapeFunctionsValues(itCond->GetIntegrationMethod());

                        // Calculate elemental normal
                        array_1d<double,3> Normal;
                        itCond->Calculate(NORMAL,Normal,rProcessInfo);

                        // Normalizing the Normals
                        double length = 0.0;
                        for(std::size_t i=0;i<TDim;i++)
                        {
                            length = length + Normal[i]*Normal[i];
                        }
                        length = sqrt(length);

                        for(std::size_t i=0;i<TDim;i++)
                        {
                            Normal[i] = Normal[i]/length;
                        }

                        //Variable<double>& InterpolatedValue;
                        std::vector< double > InterpolatedValue = std::vector< double >();



                        double FluxAtGaussPoint;
                        bool Found = false;
                        for (std::size_t g = 0; g < NumGauss; g++)
                        {
                            // Obtain coordinates of integration point
                            array_1d<double,3> Coordinates(3,0.0);
                            for (std::size_t i = 0; i < rGeom.PointsNumber(); i++)
                                Coordinates += NCond(g,i)*rGeom[i].Coordinates();

                            //KRATOS_WATCH(Coordinates)
                            //KRATOS_WATCH(IntegrationPoints[g].Coordinates())

                            // Find element containing the integration point
                            Element::Pointer pElem;

                            Kratos::Vector NElem; // element shape functions evaluated at condition integration point
                            Found = SearchUtil.FindPointOnMeshSimplified(Coordinates,NElem,pElem);

                            if (Found)
                            {
                                // Interpolate a value for the integration point using pElem

                                Calculate_Flux_Gauss( Normal,  NElem, pElem , rOriginVariable , FluxAtGaussPoint);
                                InterpolatedValue.push_back(FluxAtGaussPoint);
                            }
                        }


                        std::vector< double > OldValue;


                        itCond->GetValueOnIntegrationPoints(flux,OldValue,rProcessInfo);
                        for (std::size_t g = 0; g < NumGauss; g++)
                        {
                            InterpolatedValue[g] = Omega * InterpolatedValue[g] + (1.0-Omega) * OldValue[g];
                        }

                        // Store results on condition
                        itCond->SetValueOnIntegrationPoints(flux,InterpolatedValue,rProcessInfo);
                    }

                }
                KRATOS_CATCH("");
            }



            /****************************************************************************************
             * Function: NeumannBoundaryCalculation --> Stress Tensor --> Nodal Integration
             * Description:
             * First the boundary nodes of the rechteck model part are figured out and
             * saved in a list.
             *
             * With the interpolation function the interpolation of a Scalar Variable of one
             * modelpart to the other one is done...
             *
             * Variable Explanation:
             * rOrigin_ModelPart: the model part  all the variable should be taken from
             * rDestination_ModelPart: the destination model part where we want to know the values of
             *                         the variables if you want to pass only one variable
             * rOriginVariable: the name of the interpolated variable in the origin model part
             * rOriginVariable: the name of the interpolated variable in the destination model part
             * *************************************************************************************/

                  //template<class TDataType>

                  void Projection_Traction_Nodal(
                      ModelPart& rOrigin_ModelPart ,
                      ModelPart& rDestination_ModelPart,
                      Variable< array_1d<double,3> >& velocity,
                      Variable< double >& rOriginVariablePress,
                      Variable< array_1d<double,3> >& traction,
                      int counter,
                      const int domain_size
                      // Variable for preallocating the array 'list_of_new_nodes'
                  )
                  {

                      KRATOS_TRY

                    /**********************************************************************
                     * Creating an array 'list_of_new_nodes' which contains all the nodes
                     * on which an interpolation should take part
                     * *******************************************************************/

                        // Create an Array to store the border Ids
                        std::vector< Node<3>::Pointer > list_of_new_nodes;
                        list_of_new_nodes.reserve(counter);

                        // creating the list
                        for(ModelPart::ElementsContainerType::iterator it = rDestination_ModelPart.ElementsBegin(); it!=rDestination_ModelPart.ElementsEnd(); it++)
                        {
                            //get the connectivity of the element
                            Geometry< Node<3> >& geom = it->GetGeometry();

                            int size = geom.size();

                            for (int i = 0; i < size; i++ )
                            {
                                if(geom[i].FastGetSolutionStepValue(FLAG_VARIABLE)==1.0)
                                {
                                 list_of_new_nodes.push_back( geom(i));
                                }
                            }
                        }

                      //*******************************************************************
                      // properties to be used in the generation
                      //*******************************************************************
                      Properties::Pointer properties = rDestination_ModelPart.GetMesh().pGetProperties(1);

                      //defintions for spatial search
                      typedef Node<3> PointType;
                      typedef Node<3>::Pointer PointTypePointer;
                      typedef std::vector<PointType::Pointer>           PointVector;
                      typedef std::vector<PointType::Pointer>::iterator PointIterator;
                      typedef std::vector<double>               DistanceVector;
                      typedef std::vector<double>::iterator     DistanceIterator;

                      //starting calculating time of construction of the kdtree
                      boost::timer kdtree_construction;

                      //*************
                      // Bucket types
                      typedef Bucket< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;

                      //*************
                      // DynamicBins;

                      typedef Tree< KDTreePartition<BucketType> > tree; 		//Kdtree;

                      KRATOS_INFO( "kdt construction time " )<< kdtree_construction.elapsed() << std::endl;
                      //finishing calculating time of construction of the kdtree

                      //create a spatial database with the list of new nodes
                      std::size_t bucket_size = 20;
                      tree nodes_tree(list_of_new_nodes.begin(),list_of_new_nodes.end(),bucket_size);

                      //work arrays
                      Node<3> work_point(0,0.0,0.0,0.0);
                      std::size_t MaximumNumberOfResults = 10000;
                      PointVector Results(MaximumNumberOfResults);
                      DistanceVector ResultsDistances(MaximumNumberOfResults);

                      array_1d<double,TDim+1> N; //Shape functions vector//

                      boost::timer search_and_interpolation_time;
                      //loop over all of the elements in the "old" list to perform the interpolation
                      for( ModelPart::ElementsContainerType::iterator el_it = rOrigin_ModelPart.ElementsBegin();
                              el_it != rOrigin_ModelPart.ElementsEnd(); el_it++)
                      {
                          Geometry<Node<3> >&geom = el_it->GetGeometry();

                          //find the center and "radius" of the element
                          double xc, yc, zc,  radius;
                          CalculateCenterAndSearchRadius( geom, xc,yc,zc, radius, N);
                          work_point.X() = xc;
                          work_point.Y() = yc;
                          work_point.Z() = zc;

                          //find all of the new nodes within the radius
                          int number_of_points_in_radius;

                          //look between the new nodes which of them is inside the radius of the circumscribed circle
                          number_of_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, Results.begin(),
                                                       ResultsDistances.begin(),  MaximumNumberOfResults);


                          //check if inside
                          for( PointIterator it_found = Results.begin(); it_found != Results.begin() + number_of_points_in_radius; it_found++)
                          {
                              bool is_inside = false;
                              //once we are sure the node in inside the circle we have to see if it is inside the triangle i.e. if all the Origin element shape functions are >1

                              //double is_visited = (*it_found)->GetValue(IS_VISITED);

                              // CALCULATE THE FLUXES
                              // Vector of shape function derivatives
                              //GeometryData::ShapeFunctionsGradientsType N_der; //Shape functions vector//
                              is_inside = CalculatePosition(geom,	(*it_found)->X(),(*it_found)->Y(),(*it_found)->Z(),N);

                              //if the node falls inside the element interpolate
                              if(is_inside == true)// && is_visited != 1.0)
                              {
                                  //Interpolating all the rOriginVariables to get their nodal value in the rDestination_ModelPart
                                  Calculate_Traction_Nodal( el_it,  N, *it_found , velocity, rOriginVariablePress,  traction, domain_size );
                                  //pnode->GetValue(IS_VISITED) = 1.0;
                              }
                          }

                      }

      /*
                      // Output of the FACE_HEAT_FLUXES
                      for(
                           ModelPart::NodesContainerType::iterator inode = rDestination_ModelPart.NodesBegin();
                           inode != rDestination_ModelPart.NodesEnd();
                           inode++)
                       {
                          if(inode->X() == 0.4)
                          {
                              KRATOS_WATCH(inode->FastGetSolutionStepValue(FACE_HEAT_FLUX));
                              KRATOS_WATCH(inode->Id());
                          }

                       }
      */
                      KRATOS_INFO( "search and interpolation time ") << search_and_interpolation_time.elapsed() << std::endl;
                      KRATOS_CATCH("")
                  }

            /****************************************************************************************
             * Function: Projection_Traction
             * Description:
             * First the boundary nodes of the rechteck model part are figured out and
             * saved in a list.
             *
             * With the interpolation function the interpolation of a Scalar Variable of one
             * modelpart to the other one is done...
             *
             * Variable Explanation:
             * rOrigin_ModelPart: the model part  all the variable should be taken from
             * rDestination_ModelPart: the destination model part where we want to know the values of
             *                         the variables if you want to pass only one variable
             * rOriginVariable: the name of the interpolated variable in the origin model part
             * rOriginVariable: the name of the interpolated variable in the destination model part
             * *************************************************************************************/

            //template<class TDataType>

            void Projection_Traction(
                    ModelPart& rOrigin_ModelPart ,
                    ModelPart& rDestination_ModelPart,
                    Variable< array_1d<double,3> >& velocity,
                    Variable< double >& rOriginVariablePress,
                    Variable< array_1d<double,3> >& traction,
                    int counter,
                    double Omega = 1.0
                    )
            {

                KRATOS_TRY;

                // Create search structure containing the origin mesh
                BinBasedFastPointLocator<TDim> SearchUtil = BinBasedFastPointLocator<TDim>(rOrigin_ModelPart);
                SearchUtil.UpdateSearchDatabase();

                ProcessInfo& rProcessInfo = rDestination_ModelPart.GetProcessInfo();

                // MARK THE BOUNDARY FACES WITH A FLAG INTERFACE
                //get the connectivity of the face

                for (ModelPart::ConditionIterator itCond = rDestination_ModelPart.ConditionsBegin(); itCond != rDestination_ModelPart.ConditionsEnd(); itCond++)
                {
                    //get the connectivity of the face
                    Geometry< Node<3> >& rGeom = itCond->GetGeometry();
                    int counter = 0;
                    for(Geometry< Node<3> >::iterator node_it = rGeom.begin(); node_it!=rGeom.end(); node_it++)
                    {
                        if(node_it->FastGetSolutionStepValue(FLAG_VARIABLE,0)==1.0)
                        {
                            counter = counter +1;
                        }
                    }
                    if(counter>1)
                    {
                        itCond->Set(INTERFACE,true);
                    }
                }



                // Loop over conditions in destination mesh, obtain gauss points
                for (ModelPart::ConditionIterator itCond = rDestination_ModelPart.ConditionsBegin(); itCond != rDestination_ModelPart.ConditionsEnd(); itCond++)
                {
                    // filter conditions, use only those on the boundary
                    if (itCond->Is(INTERFACE))
                    {
                        // Obtain integration points on the condition
                        Geometry< Node<3> >& rGeom = itCond->GetGeometry();
                        const GeometryData::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(itCond->GetIntegrationMethod());
                        const std::size_t NumGauss = IntegrationPoints.size();
                        Matrix NCond = rGeom.ShapeFunctionsValues(itCond->GetIntegrationMethod());

                        // Calculate elemental normal
                        array_1d<double,3> Normal;
                        itCond->Calculate(NORMAL,Normal,rProcessInfo);

                        // Normalizing the Normals
                        double length = 0.0;
                        for(std::size_t i=0;i<TDim;i++)
                        {
                            length = length + Normal[i]*Normal[i];
                        }
                        length = sqrt(length);

                        for(std::size_t i=0;i<TDim;i++)
                        {
                            Normal[i] = Normal[i]/length;
                        }

                        std::vector<array_1d<double, 3 > > InterpolatedValues;
                        array_1d<double,3> traction_result(3,0.0);
                        bool Found = false;
                        for (std::size_t g = 0; g < NumGauss; g++)
                        {
                            // Obtain coordinates of integration point
                            array_1d<double,3> Coordinates(3,0.0);
                            for (std::size_t i = 0; i < rGeom.PointsNumber(); i++)
                                Coordinates += NCond(g,i)*rGeom[i].Coordinates();

                            // Find element containing the integration point
                            Element::Pointer pElem;
                            Vector NElem; // element shape functions evaluated at condition integration point
                            Found = SearchUtil.FindPointOnMeshSimplified(Coordinates,NElem,pElem);

                            if (Found)
                            {
                                // Interpolate a value for the integration point using pElem
                                //Calculate_Traction( el_it,  N, *it_found , velocity, rOriginVariablePress,  traction, domain_size );
                                Calculate_Traction( Normal,  NElem, pElem , velocity, rOriginVariablePress,  traction_result );
                                InterpolatedValues.push_back(traction_result);
                            }
                        }

                        std::vector< array_1d<double,3> > OldValues;
                        //itCond->GetValueOnIntegrationPoints(traction,OldValues,rProcessInfo);
                        itCond->GetValueOnIntegrationPoints(traction,OldValues,rProcessInfo);
                        for (std::size_t g = 0; g < NumGauss; g++)
                        {
                            InterpolatedValues[g] = Omega * InterpolatedValues[g] + (1.0-Omega) * OldValues[g];
                        }

                        // Store results on condition
                        itCond->SetValueOnIntegrationPoints(traction,InterpolatedValues,rProcessInfo);
/* The following can be used to apply Dirchlet BC on the inflow and Neumann on the outflow
                        for (std::size_t n = 0; n < rGeom.size(); n++)
                        {
                            const array_1d<double,3>& Coordinates = rGeom[n].Coordinates();
                            Element::Pointer pElem;
                            Vector NElem;
                            bool Found = SearchUtil.FindPointOnMeshSimplified(Coordinates,NElem,pElem);

                            if (Found)
                            {
                                //Geometry element of the rOrigin_ModelPart
                                Geometry< Node<3> >& ElemGeom = pElem->GetGeometry();

                                //getting the data of the solution step
//                                array_1d<double,3>& step_data = rGeom[n].FastGetSolutionStepValue(VELOCITY); // Corresponds to vx
//                                step_data = NElem[0] * ElemGeom[0].FastGetSolutionStepValue(VELOCITY);
//                                for (std::size_t i = 1; i < ElemGeom.size(); i++)
//                                    step_data += NElem[i] * ElemGeom[i].FastGetSolutionStepValue(VELOCITY);


                                array_1d<double,3> step_data = NElem[0] * ElemGeom[0].FastGetSolutionStepValue(VELOCITY);
                                for (std::size_t i = 1; i < ElemGeom.size(); i++)
                                    step_data += NElem[i] * ElemGeom[i].FastGetSolutionStepValue(VELOCITY);

                                rGeom[n].FastGetSolutionStepValue(VELOCITY) = step_data;

                                if( inner_prod(step_data, Normal) < 0.0)
                                {
                                    rGeom[n].Fix(VELOCITY_X);
                                    rGeom[n].Fix(VELOCITY_Y);
                                }
                                else
                                {
                                    rGeom[n].Free(VELOCITY_X);
                                    rGeom[n].Free(VELOCITY_Y);

                                }
                            }
                        }*/
                    }

                }
                KRATOS_CATCH("");
            }


            /****************************************************************************************
             * Function: Projection_Traction_Cavity
             * Description:
             * First the boundary nodes of the rechteck model part are figured out and
             * saved in a list.
             *
             * With the interpolation function the interpolation of a Scalar Variable of one
             * modelpart to the other one is done...
             *
             * Variable Explanation:
             * rOrigin_ModelPart: the model part  all the variable should be taken from
             * rDestination_ModelPart: the destination model part where we want to know the values of
             *                         the variables if you want to pass only one variable
             * rOriginVariable: the name of the interpolated variable in the origin model part
             * rOriginVariable: the name of the interpolated variable in the destination model part
             * *************************************************************************************/

            //template<class TDataType>

            void Projection_Traction_Cavity(
                    ModelPart& rOrigin_ModelPart ,
                    ModelPart& rDestination_ModelPart,
                    Variable< array_1d<double,3> >& velocity,
                    Variable< double >& rOriginVariablePress,
                    Variable< array_1d<double,3> >& traction,
                    int counter,
                    double Omega = 1.0
                    )
            {

                KRATOS_TRY;

                // Create search structure containing the origin mesh
                BinBasedFastPointLocator<TDim> SearchUtil = BinBasedFastPointLocator<TDim>(rOrigin_ModelPart);
                SearchUtil.UpdateSearchDatabase();

                ProcessInfo& rProcessInfo = rDestination_ModelPart.GetProcessInfo();

                // MARK THE BOUNDARY FACES WITH A FLAG INTERFACE
                //get the connectivity of the face

                for (ModelPart::ConditionIterator itCond = rDestination_ModelPart.ConditionsBegin(); itCond != rDestination_ModelPart.ConditionsEnd(); itCond++)
                {
                    //get the connectivity of the face
                    Geometry< Node<3> >& rGeom = itCond->GetGeometry();
                    int counter = 0;
                    for(Geometry< Node<3> >::iterator node_it = rGeom.begin(); node_it!=rGeom.end(); node_it++)
                    {
                        if(node_it->FastGetSolutionStepValue(FLAG_VARIABLE,0)==1.0)
                        {
                            counter = counter +1;
                        }
                    }
                    if(counter>1)
                    {
                        itCond->Set(INTERFACE,true);
                    }
                }



                // Loop over conditions in destination mesh, obtain gauss points
                for (ModelPart::ConditionIterator itCond = rDestination_ModelPart.ConditionsBegin(); itCond != rDestination_ModelPart.ConditionsEnd(); itCond++)
                {
                    // filter conditions, use only those on the boundary
                    if (itCond->Is(INTERFACE))
                    {
                        // Obtain integration points on the condition
                        Geometry< Node<3> >& rGeom = itCond->GetGeometry();
                        const GeometryData::IntegrationPointsArrayType& IntegrationPoints = rGeom.IntegrationPoints(itCond->GetIntegrationMethod());
                        const std::size_t NumGauss = IntegrationPoints.size();
                        Matrix NCond = rGeom.ShapeFunctionsValues(itCond->GetIntegrationMethod());

                        // Calculate elemental normal
                        array_1d<double,3> Normal;
                        itCond->Calculate(NORMAL,Normal,rProcessInfo);

                        // Normalizing the Normals
                        double length = 0.0;
                        for(std::size_t i=0;i<TDim;i++)
                        {
                            length = length + Normal[i]*Normal[i];
                        }
                        length = sqrt(length);

                        for(std::size_t i=0;i<TDim;i++)
                        {
                            Normal[i] = Normal[i]/length;
                        }

                        std::vector<array_1d<double, 3 > > InterpolatedValues;
                        array_1d<double,3> traction_result(3,0.0);
                        bool Found = false;
                        for (std::size_t g = 0; g < NumGauss; g++)
                        {
                            // Obtain coordinates of integration point
                            array_1d<double,3> Coordinates(3,0.0);
                            for (std::size_t i = 0; i < rGeom.PointsNumber(); i++)
                                Coordinates += NCond(g,i)*rGeom[i].Coordinates();

                            // Find element containing the integration point
                            Element::Pointer pElem;
                            Vector NElem; // element shape functions evaluated at condition integration point
                            Found = SearchUtil.FindPointOnMeshSimplified(Coordinates,NElem,pElem);

                            if (Found)
                            {
                                // Interpolate a value for the integration point using pElem
                                //Calculate_Traction( el_it,  N, *it_found , velocity, rOriginVariablePress,  traction, domain_size );
                                Calculate_Traction( Normal,  NElem, pElem , velocity, rOriginVariablePress,  traction_result );
                                InterpolatedValues.push_back(traction_result);
                            }
                        }

                        std::vector< array_1d<double,3> > OldValues;
                        //itCond->GetValueOnIntegrationPoints(traction,OldValues,rProcessInfo);
                        itCond->GetValueOnIntegrationPoints(traction,OldValues,rProcessInfo);
                        for (std::size_t g = 0; g < NumGauss; g++)
                        {
                            InterpolatedValues[g] = Omega * InterpolatedValues[g] + (1.0-Omega) * OldValues[g];
                        }

                        // Store results on condition
                        itCond->SetValueOnIntegrationPoints(traction,InterpolatedValues,rProcessInfo);
// The following can be used to apply Dirchlet BC on the inflow and Neumann on the outflow
                        for (std::size_t n = 0; n < rGeom.size(); n++)
                        {
                            const array_1d<double,3>& Coordinates = rGeom[n].Coordinates();
                            Element::Pointer pElem;
                            Vector NElem;
                            bool Found = SearchUtil.FindPointOnMeshSimplified(Coordinates,NElem,pElem);

                            if (Found)
                            {
                                //Geometry element of the rOrigin_ModelPart
                                Geometry< Node<3> >& ElemGeom = pElem->GetGeometry();

                                //getting the data of the solution step
//                                array_1d<double,3>& step_data = rGeom[n].FastGetSolutionStepValue(VELOCITY); // Corresponds to vx
//                                step_data = NElem[0] * ElemGeom[0].FastGetSolutionStepValue(VELOCITY);
//                                for (std::size_t i = 1; i < ElemGeom.size(); i++)
//                                    step_data += NElem[i] * ElemGeom[i].FastGetSolutionStepValue(VELOCITY);


                                array_1d<double,3> step_data = NElem[0] * ElemGeom[0].FastGetSolutionStepValue(VELOCITY);
                                for (std::size_t i = 1; i < ElemGeom.size(); i++)
                                    step_data += NElem[i] * ElemGeom[i].FastGetSolutionStepValue(VELOCITY);

                                rGeom[n].FastGetSolutionStepValue(VELOCITY) = step_data;

                                if( inner_prod(step_data, Normal) < 0.0)
                                {
                                    rGeom[n].Fix(VELOCITY_X);
                                    rGeom[n].Fix(VELOCITY_Y);
                                }
                                else
                                {
                                    rGeom[n].Free(VELOCITY_X);
                                    rGeom[n].Free(VELOCITY_Y);

                                }
                            }
                        }
                    }

                }
                KRATOS_CATCH("");
            }





      /**************************************************************************************
       * Dirichlet-Boundary-Interpolation-Calculation-Function for a Geometry with Flag
       * ************************************************************************************
       * Description:
       * First the right boundary nodes of the segment1 model part are figured out and
       * saved in a list.
       *
       * With the interpolation function the interpolation of the Temperature-values of one
       * modelpart to the other one is done...
       *
       * Variable Explanation:
       * rOrigin_ModelPart: the model part  all the variable should be taken from
       * rDestination_ModelPart: the destination model part where we want to know the values of
       *                         the variables if you want to pass only one variable
       * rOriginVariable: the name of the interpolated variable in the origin model part
       * rOriginVariable: the name of the interpolated variable in the destination model part
       * *************************************************************************************/
      template<class TDataType>
      void Projection_Vector(
          ModelPart& rOrigin_ModelPart ,
          ModelPart& rDestination_ModelPart,
          Variable<TDataType>& rOriginVariable ,
          Variable<TDataType>& rDestinationVariable,
          int i,
          const int domain_size,
              const double Omega = 1.0)
      {

          KRATOS_TRY

          // Create an Array to store the border Ids
          // All Nodes with a flag in mesh 2 are border Ids
          std::vector< Node<3>::Pointer > list_of_new_nodes;
          //list_of_new_nodes.reserve(2);


          for(ModelPart::NodeIterator i_node = rDestination_ModelPart.NodesBegin();
               i_node != rDestination_ModelPart.NodesEnd();
               i_node++)
           {
              if (i_node->FastGetSolutionStepValue(FLAG_VARIABLE)==1.0)
              {
                  //KRATOS_WATCH(i_node->Id()); //print the value of the node Id which belongs to chim_mesh
                  list_of_new_nodes.push_back( *(i_node.base()));
                  i=i+1;
              }
           }

          if (list_of_new_nodes.size() == 0)
          {
              KRATOS_THROW_ERROR(std::runtime_error,"No destination nodes found for interpolation. Assign FLAG_VARIABLE=1 to destination nodes","");
          }

          //*******************************************************************
          //properties to be used in the generation
          //*******************************************************************
          Properties::Pointer properties = rDestination_ModelPart.GetMesh().pGetProperties(1);

          //defintions for spatial search
          typedef Node<3> PointType;
          typedef Node<3>::Pointer PointTypePointer;
          typedef std::vector<PointType::Pointer>           PointVector;
          typedef std::vector<PointType::Pointer>::iterator PointIterator;
          typedef std::vector<double>               DistanceVector;
          typedef std::vector<double>::iterator     DistanceIterator;

          //starting calculating time of construction of the kdtree
          boost::timer kdtree_construction;

          //*************
          // Bucket types
          typedef Bucket< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;

          //*************
          // DynamicBins;

          typedef Tree< KDTreePartition<BucketType> > tree; 		//Kdtree;

          KRATOS_INFO( "kdt construction time ") << kdtree_construction.elapsed() << std::endl;
          //finishing calculating time of construction of the kdtree

          //create a spatial database with the list of new nodes
          std::size_t bucket_size = 20;
          tree nodes_tree(list_of_new_nodes.begin(),list_of_new_nodes.end(),bucket_size);

          //work arrays
          Node<3> work_point(0,0.0,0.0,0.0);
          std::size_t MaximumNumberOfResults = 10000;
          PointVector Results(MaximumNumberOfResults);
          DistanceVector ResultsDistances(MaximumNumberOfResults);
          array_1d<double,TDim+1> N; //Shape functions vector//

          boost::timer search_and_interpolation_time;
          //loop over all of the elements in the "old" list to perform the interpolation
          for( ModelPart::ElementsContainerType::iterator el_it = rOrigin_ModelPart.ElementsBegin();
                  el_it != rOrigin_ModelPart.ElementsEnd(); el_it++)
          {
              Geometry<Node<3> >&geom = el_it->GetGeometry();

              //find the center and "radius" of the element
              double xc, yc, zc,  radius;
              CalculateCenterAndSearchRadius( geom, xc,yc,zc, radius, N);
              work_point.X() = xc;
              work_point.Y() = yc;
              work_point.Z() = zc;

              //find all of the new nodes within the radius
              int number_of_points_in_radius;

              //look between the new nodes which of them is inside the radius of the circumscribed circle
              number_of_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, Results.begin(),
                                           ResultsDistances.begin(),  MaximumNumberOfResults);


              //check if inside
              for( PointIterator it_found = Results.begin(); it_found != Results.begin() + number_of_points_in_radius; it_found++)
              {
                  bool is_inside = false;
                  //once we are sure the node in inside the circle we have to see if it is inside the triangle i.e. if all the Origin element shape functions are >1

                  //double is_visited = (*it_found)->GetValue(IS_VISITED);
                  is_inside = CalculatePosition(geom,	(*it_found)->X(),(*it_found)->Y(),(*it_found)->Z(),N);

                  //if the node falls inside the element interpolate
                  if(is_inside == true)// && is_visited != 1.0)
                  {
                      //Interpolating all the rOriginVariables to get their nodal value in the rDestination_ModelPart
                      Interpolate_Dirichlet(el_it, N, *it_found, rOriginVariable, rDestinationVariable, Omega);

                      //pnode->GetValue(IS_VISITED) = 1.0;
                      //KRATOS_WATCH( (*it_found)->FastGetSolutionStepValue(VELOCITY))
                  }
              }
          }
          KRATOS_INFO( "search and interpolation time ") << search_and_interpolation_time.elapsed() << std::endl;
          KRATOS_CATCH("");
      }


      template< class TValue >
      void ApplyRelaxation(TValue& rNewValues, const TValue& rOldValues, double Omega)
      {
          rNewValues = Omega*rNewValues + (1.0-Omega)*rOldValues;
      }


      void ExtractBoundaryFaces(ModelPart& rModelPart, ModelPart& rOutModelPart);


      ///@}
      ///@name Access
      ///@{


      ///@}
      ///@name Inquiry
      ///@{


      ///@}
      ///@name Input and output
      ///@{
      /// Turn back information as a stemplate<class T, std::size_t dim> string.
          virtual std::string Info() const
          {
              return "";
          }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const {}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const {}


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
      // Implementation of the Functions used above for Triangles
      inline void CalculateCenterAndSearchRadius(Geometry<Node<3> >&geom,
                  double& xc, double& yc, double& zc, double& R, array_1d<double,3>& N
                                                    )
          {
              double x0 = geom[0].X();
              double  y0 = geom[0].Y();
              double x1 = geom[1].X();
              double  y1 = geom[1].Y();
              double x2 = geom[2].X();
              double  y2 = geom[2].Y();


              xc = 0.3333333333333333333*(x0+x1+x2);
              yc = 0.3333333333333333333*(y0+y1+y2);
              zc = 0.0;

              double R1 = (xc-x0)*(xc-x0) + (yc-y0)*(yc-y0);
              double R2 = (xc-x1)*(xc-x1) + (yc-y1)*(yc-y1);
              double R3 = (xc-x2)*(xc-x2) + (yc-y2)*(yc-y2);

              R = R1;
              if(R2 > R) R = R2;
              if(R3 > R) R = R3;

              R = 1.01 * sqrt(R);
          }

      //***************************************
      //***************************************
      inline double CalculateVol(	const double x0, const double y0,
                                    const double x1, const double y1,
                                    const double x2, const double y2
                                 )
          {
              return 0.5*( (x1-x0)*(y2-y0)- (y1-y0)*(x2-x0) );
          }
      //***************************************
      //***************************************

      inline bool CalculatePosition(  Geometry<Node<3> >&geom,
                                      const double xc, const double yc, const double zc,
                                      array_1d<double,3>& N
                                   )
      {
          double x0 = geom[0].X();
          double  y0 = geom[0].Y();
          double x1 = geom[1].X();
          double  y1 = geom[1].Y();
          double x2 = geom[2].X();
          double  y2 = geom[2].Y();

          double area = CalculateVol(x0,y0,x1,y1,x2,y2);
          double inv_area = 0.0;
          if(area == 0.0)
          {

  // 				KRATOS_THROW_ERROR(std::logic_error,"element with zero area found","");
              //The interpolated node will not be inside an elemente with zero area
              return false;

          }
          else
          {
              inv_area = 1.0 / area;
          }


          N[0] = CalculateVol(x1,y1,x2,y2,xc,yc) * inv_area;
          N[1] = CalculateVol(x2,y2,x0,y0,xc,yc) * inv_area;
          N[2] = CalculateVol(x0,y0,x1,y1,xc,yc) * inv_area;


          if(N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <=1.0 && N[1]<= 1.0 && N[2] <= 1.0) //if the xc yc is inside the triangle return true
              return true;

          return false;
      }

      //***************************************

      inline bool CalculatePosition1(  Geometry<Node<3> >&geom,
                                      const double xc, const double yc, const double zc,
                                      array_1d<double,3>& N,
                                      GeometryData::ShapeFunctionsGradientsType& N_der
                                   )
      {
          double x0 = geom[0].X();
          double  y0 = geom[0].Y();
          double x1 = geom[1].X();
          double  y1 = geom[1].Y();
          double x2 = geom[2].X();
          double  y2 = geom[2].Y();

          double area = CalculateVol(x0,y0,x1,y1,x2,y2);
          double inv_area = 0.0;
          if(area == 0.0)
          {

  // 				KRATOS_THROW_ERROR(std::logic_error,"element with zero area found","");
              //The interpolated node will not be inside an elemente with zero area
              return false;

          }
          else
          {
              inv_area = 1.0 / area;
          }


          N[0] = CalculateVol(x1,y1,x2,y2,xc,yc) * inv_area;
          N[1] = CalculateVol(x2,y2,x0,y0,xc,yc) * inv_area;
          N[2] = CalculateVol(x0,y0,x1,y1,xc,yc) * inv_area;

          //Implementation of the Derivative of the Shape functions
          geom.ShapeFunctionsIntegrationPointsGradients(N_der,GeometryData::GI_GAUSS_1);


          if(N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <=1.0 && N[1]<= 1.0 && N[2] <= 1.0) //if the xc yc is inside the triangle return true
              return true;

          return false;
      }

      inline void CalculateCenterAndSearchRadius(Geometry<Node<3> >&geom,
              double& xc, double& yc, double& zc, double& R, array_1d<double,4>& N

                                                )
      {
          double x0 = geom[0].X();
          double  y0 = geom[0].Y();
          double  z0 = geom[0].Z();
          double x1 = geom[1].X();
          double  y1 = geom[1].Y();
          double  z1 = geom[1].Z();
          double x2 = geom[2].X();
          double  y2 = geom[2].Y();
          double  z2 = geom[2].Z();
          double x3 = geom[3].X();
          double  y3 = geom[3].Y();
          double  z3 = geom[3].Z();


          xc = 0.25*(x0+x1+x2+x3);
          yc = 0.25*(y0+y1+y2+y3);
          zc = 0.25*(z0+z1+z2+z3);

          double R1 = (xc-x0)*(xc-x0) + (yc-y0)*(yc-y0) + (zc-z0)*(zc-z0);
          double R2 = (xc-x1)*(xc-x1) + (yc-y1)*(yc-y1) + (zc-z1)*(zc-z1);
          double R3 = (xc-x2)*(xc-x2) + (yc-y2)*(yc-y2) + (zc-z2)*(zc-z2);
          double R4 = (xc-x3)*(xc-x3) + (yc-y3)*(yc-y3) + (zc-z3)*(zc-z3);

          R = R1;
          if(R2 > R) R = R2;
          if(R3 > R) R = R3;
          if(R4 > R) R = R4;

          R = sqrt(R);
      }

      //***************************************
      //***************************************
      inline double CalculateVol(	const double x0, const double y0, const double z0,
                                  const double x1, const double y1, const double z1,
                                  const double x2, const double y2, const double z2,
                                  const double x3, const double y3, const double z3
                                )
      {
          double x10 = x1 - x0;
          double y10 = y1 - y0;
          double z10 = z1 - z0;

          double x20 = x2 - x0;
          double y20 = y2 - y0;
          double z20 = z2 - z0;

          double x30 = x3 - x0;
          double y30 = y3 - y0;
          double z30 = z3 - z0;

          double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
          return  detJ*0.1666666666666666666667;

          //return 0.5*( (x1-x0)*(y2-y0)- (y1-y0)*(x2-x0) );
      }


      //***************************************
      //***************************************

      inline bool CalculatePosition(	Geometry<Node<3> >&geom,
                                      const double xc, const double yc, const double zc,
                                      array_1d<double,4>& N
                                   )
      {

          double x0 = geom[0].X();
          double  y0 = geom[0].Y();
          double  z0 = geom[0].Z();
          double x1 = geom[1].X();
          double  y1 = geom[1].Y();
          double  z1 = geom[1].Z();
          double x2 = geom[2].X();
          double  y2 = geom[2].Y();
          double  z2 = geom[2].Z();
          double x3 = geom[3].X();
          double  y3 = geom[3].Y();
          double  z3 = geom[3].Z();

          double vol = CalculateVol(x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3);

          double inv_vol = 0.0;
          if(vol < 0.0000000000001)
          {

  // 				KRATOS_THROW_ERROR(std::logic_error,"element with zero vol found","");
              //The interpolated node will not be inside an elemente with zero volume
              return false;
  // 				KRATOS_WATCH("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
          }
          else
          {
              inv_vol = 1.0 / vol;
          }

          N[0] = CalculateVol(x1,y1,z1,x3,y3,z3,x2,y2,z2,xc,yc,zc) * inv_vol;
          N[1] = CalculateVol(x3,y3,z3,x0,y0,z0,x2,y2,z2,xc,yc,zc) * inv_vol;
          N[2] = CalculateVol(x3,y3,z3,x1,y1,z1,x0,y0,z0,xc,yc,zc) * inv_vol;
          N[3] = CalculateVol(x0,y0,z0,x1,y1,z1,x2,y2,z2,xc,yc,zc) * inv_vol;


          if(N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[3] >=0.0 &&
                  N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0 && N[3] <=1.0)
              //if the xc yc zc is inside the tetrahedron return true
              return true;

          return false;
      }

      //***************************************
      //***************************************
      //projecting an array1D 2Dversion
          void Interpolate(
              ModelPart::ElementsContainerType::iterator el_it,
              const array_1d<double,3>& N,
              Node<3>::Pointer pnode,
              Variable<array_1d<double,3> >& rOriginVariable,
              Variable<array_1d<double,3> >& rDestinationVariable)
          {
              //Geometry element of the rOrigin_ModelPart
              Geometry< Node<3> >& geom = el_it->GetGeometry();

              std::size_t buffer_size = pnode->GetBufferSize();

              for(std::size_t step = 0; step<buffer_size; step++)
              {
                  //getting the data of the solution step
                  array_1d<double,3>& step_data = (pnode)->FastGetSolutionStepValue(rDestinationVariable , step);
                  //Reference or no reference???//CANCELLA
                  const array_1d<double,3>& node0_data = geom[0].FastGetSolutionStepValue(rOriginVariable , step);
                  const array_1d<double,3>& node1_data = geom[1].FastGetSolutionStepValue(rOriginVariable , step);
                  const array_1d<double,3>& node2_data = geom[2].FastGetSolutionStepValue(rOriginVariable , step);

                  //copying this data in the position of the vector we are interested in
                  for(std::size_t j= 0; j< TDim; j++)
                  {
                      step_data[j] = N[0]*node0_data[j] + N[1]*node1_data[j] + N[2]*node2_data[j];
                  }
              }

          }



// Velocity Interpolation Vx, Vy
          //projecting a scalar 2Dversion
          void Interpolate(
              ModelPart::ElementsContainerType::iterator el_it,
              const array_1d<double,3>& N,
              Node<3>::Pointer pnode,
              Variable<double>& rOriginVariable,
              Variable<double>& rDestinationVariable)
          {
              //Geometry element of the rOrigin_ModelPart
              Geometry< Node<3> >& geom = el_it->GetGeometry();

              std::size_t buffer_size = pnode->GetBufferSize();
              //facendo un loop sugli step temporali step_data come salva i dati al passo anteriore? Cioś dove passiamo l'informazione ai nodi???
              for(std::size_t step = 0; step<buffer_size; step++)
              {
                  //getting the data of the solution step
                  double& step_data = (pnode)->FastGetSolutionStepValue(rDestinationVariable , step);
                  //Reference or no reference???//CANCELLA
                  const double node0_data = geom[0].FastGetSolutionStepValue(rOriginVariable , step);
                  const double node1_data = geom[1].FastGetSolutionStepValue(rOriginVariable , step);
                  const double node2_data = geom[2].FastGetSolutionStepValue(rOriginVariable , step);

                  //copying this data in the position of the vector we are interested in

                  step_data = N[0]*node0_data + N[1]*node1_data + N[2]*node2_data;
              }
          }


          //projecting a scalar 2Dversion
          template< class TDataType >
          void Interpolate_Dirichlet(
                  ModelPart::ElementsContainerType::iterator el_it,
                  const array_1d<double,TDim+1>& N,
                  Node<3>::Pointer pnode,
                  Variable<TDataType>& rOriginVariable,
                  Variable<TDataType>& rDestinationVariable,
                  const double Omega = 1.0)
          {
              //Geometry element of the rOrigin_ModelPart
              Geometry< Node<3> >& geom = el_it->GetGeometry();

              std::size_t buffer_size = pnode->GetBufferSize();
              for(std::size_t step = 0; step<buffer_size; step++)
              {
                  //getting the data of the solution step
                  TDataType& step_data = (pnode)->FastGetSolutionStepValue(rDestinationVariable , step); // Corresponds to vx
                  step_data *= (1.0-Omega);
                  for (std::size_t i = 0; i < N.size(); i++)
                      step_data += Omega * N[i] * geom[i].FastGetSolutionStepValue(rOriginVariable , step);
              }
          }

          //projecting a scalar 2Dversion
          void Calculate_Flux(
              ModelPart::ElementsContainerType::iterator el_it,
              const array_1d<double,3>& N,
              Node<3>::Pointer pnode,
              Variable<double>& rOriginVariable,
              Variable< double >& rDestinationVariable,
              const int domain_size
              )
          {
              //Geometry element of the rOrigin_ModelPart
              Geometry< Node<3> >& geom = el_it->GetGeometry();
              //Implementation of the Derivative of the Shape functions
              GeometryData::ShapeFunctionsGradientsType N_der;
              geom.ShapeFunctionsIntegrationPointsGradients(N_der,GeometryData::GI_GAUSS_1);
              Matrix& DNDX = N_der[0];
              //KRATOS_WATCH(N_der);
              //KRATOS_WATCH(DNDX);
              const int NumNodes = TDim + 1;
              array_1d<double,3> Normal = (pnode)->FastGetSolutionStepValue(NORMAL);
              array_1d<double,3> Gradient(3,0.0);
              //getting the data of the solution step
              double& flux = (pnode)->FastGetSolutionStepValue(rDestinationVariable);

              // Normalizing the Normals
              double length = 0.0;
              for(int i=0;i<domain_size;i++)
              {
                  length = length + Normal[i]*Normal[i];
              }
              length = sqrt(length);

              for(int i=0;i<domain_size;i++)
              {
                  Normal[i] = Normal[i]/length;
              }

              // Temperature of the nodes of the triangles
              array_1d<double,NumNodes> node_data;
              for(int i=0; i<NumNodes; i++)
              {
                  node_data[i] = geom[i].FastGetSolutionStepValue(rOriginVariable);
              }

              //copying this data in the position of the vector we are interested in
              //array_1d<double,TDim> Gradient;
              for(std::size_t i=0; i<TDim; i++)
              {
                  for(int j=0; j<NumNodes; j++)
                  {
                      Gradient[i] += DNDX(j,i)*node_data[j];
                  }
              }
              /*
              Gradient[0]= DNDX(0,0)*node0_data + DNDX(1,0)*node1_data + DNDX(2,0)*node2_data;
              Gradient[1]= DNDX(0,1)*node0_data + DNDX(1,1)*node1_data + DNDX(2,1)*node2_data;
              */

              for(int i=0; i<domain_size; i++)
              {
                  flux += Gradient[i]*Normal[i];
              }
          }


          //projecting a scalar 2Dversion
          void Calculate_Flux_Gauss(
              array_1d<double,3> Normal,
              const Vector& NElem,
              Element::Pointer pElem,
              Variable<double>& rOriginVariable,
              double& flux
              )
          {
              flux = 0.0;
              //Geometry element of the rOrigin_ModelPart
              Geometry< Node<3> >& geom = pElem->GetGeometry();
              const int NumNodes = TDim +1;

              //Implementation of the Derivative of the Shape functions
              GeometryData::ShapeFunctionsGradientsType N_der;
              geom.ShapeFunctionsIntegrationPointsGradients(N_der,GeometryData::GI_GAUSS_1);
              Matrix& DNDX = N_der[0];

              array_1d<double,3> Gradient(3,0.0);

              // Temperature of the nodes of the triangles
              array_1d<double,NumNodes> node_data;
              for(int i=0; i<NumNodes; i++)
              {
                  node_data[i] = geom[i].FastGetSolutionStepValue(rOriginVariable);
              }

              //copying this data in the position of the vector we are interested in
              //array_1d<double,TDim> Gradient;
              for(std::size_t i=0; i<TDim; i++)
              {
                  for(int j=0; j<NumNodes; j++)
                  {
                      Gradient[i] += DNDX(j,i)*node_data[j];
                  }
              }

              for(std::size_t i=0; i<TDim; i++)
              {
                  flux += Gradient[i]*Normal[i];
              }
          }



          //projecting a scalar 2Dversion
        void Calculate_Traction_Nodal(
            ModelPart::ElementsContainerType::iterator el_it,
            const array_1d<double,3>& N,
            Node<3>::Pointer pnode,
            Variable< array_1d<double,3> >& velocity,      // velocity in x-dir
            //Variable<double>& rOriginVariablevy,        // velocity in y-dir
            Variable<double>& rOriginVariablePress,      // Pressure value to interpolate
            //Variable<double>& rDestinationVariable1,  // traction force in x-dir
            Variable< array_1d<double,3> >& traction,  // traction force in y-dir
            const int domain_size)
        {
            // Geometry element of the rOrigin_ModelPart
            Geometry< Node<3> >& geom = el_it->GetGeometry();
            const int NumNodes = TDim + 1; //geom.PointsNumber();

            // Implementation of the Derivative of the Shape functions
            GeometryData::ShapeFunctionsGradientsType N_der;
            geom.ShapeFunctionsIntegrationPointsGradients(N_der,GeometryData::GI_GAUSS_1);
            Matrix& DNDX = N_der[0];

            // Getting the already calculated values for the Normals
            array_1d<double,3> Normal = (pnode)->FastGetSolutionStepValue(NORMAL);

            // Creating a variable where the interpolated pressure value is stored
            double pressure = 0.0;

            // getting the data of the solution step
            // Creating pointers to store the Solution
            array_1d<double,3>& tractions = (pnode)->FastGetSolutionStepValue(traction);

            // Normalizing the Normals
            double length = 0.0;
            for(int i=0;i<domain_size;i++)
            {
                length = length + Normal[i]*Normal[i];
            }
            length = sqrt(length);

            for(int i=0;i<domain_size;i++)
            {
                Normal[i] = Normal[i]/length;
            }

            // Getting the pressure values needed for the interpolation
            array_1d<double,NumNodes> node_press;

            for(int i=0; i<NumNodes; i++)
            {
                node_press[i] = geom[i].FastGetSolutionStepValue(rOriginVariablePress);
            }

            // Gradient_ij = Du_i/Dx_j
            Matrix Gradient = ZeroMatrix(domain_size,domain_size);
            for (int n = 0; n < NumNodes; n++)
            {
                // Velocities vx of the nodes of the triangles
                const array_1d<double,3> & v_node = geom[n].FastGetSolutionStepValue(velocity);
                //KRATOS_WATCH(v_node)
                for (int i = 0; i < domain_size; i++)
                    for (int j = 0; j < domain_size; j++)
                        Gradient(i,j) += DNDX(n,j)*v_node[i];
            }

            // Creating the interpolation for the pressure
            for(int i=0;i<NumNodes;i++)
            {
                pressure = pressure + N[i]*node_press[i];
            }

            // Creating the symmetrix velocity gradients (gradv + gradv')
            Matrix matrix_der_sym = Gradient + boost::numeric::ublas::trans(Gradient);

            // Multiplication of the interpolated pressure value with the normal vector
            array_1d<double,3> press_norm(3,0.0);
            for(int i=0; i<domain_size+1; i++)
            {
                press_norm[i] = pressure * Normal[i];
            }

            // Coefficient of viscosity and density
            // Creating the interpolation for viscosity
            double viscosity = 0;
            double density = 0;

            for(int i=0; i<NumNodes; i++)
            {
                viscosity = viscosity + N[i]*geom[i].FastGetSolutionStepValue(VISCOSITY);
                density   = density  + N[i]*geom[i].FastGetSolutionStepValue(DENSITY);
            }

            // Multiplication with symmetric velocity gradient with normals and coefficient nu
            array_1d<double,3> vel_norm(domain_size,0.0);
            for(int i=0; i<domain_size; i++)
            {
                for(int j=0;j<domain_size;j++)
                {
                    vel_norm[i] = vel_norm[i] + matrix_der_sym(i,j) * Normal(j);
                }
            }

            // Calculate the traction forces

            for(int i=0; i<domain_size; i++)
            {
                tractions[i] = -press_norm[i] + viscosity * density * vel_norm[i];
            }

            //KRATOS_WATCH(tractions);
        }

          //projecting a scalar 2Dversion
          void Calculate_Traction(
              array_1d<double,3> Normal,
//              const array_1d<double,3>& NElem,
              const Vector& NElem,
              Element::Pointer pElem,
              Variable< array_1d<double,3> >& velocity,      // velocity in x-dir
              //Variable<double>& rOriginVariablevy,        // velocity in y-dir
              Variable<double>& rOriginVariablePress,      // Pressure value to interpolate
              //Variable<double>& rDestinationVariable1,  // traction force in x-dir
              array_1d<double,3>& traction  // traction force in y-dir
              )
          {
              // Geometry element of the rOrigin_ModelPart
              Geometry< Node<3> >& geom = pElem->GetGeometry();
              const int NumNodes = TDim + 1; //geom.PointsNumber();

              // Implementation of the Derivative of the Shape functions
              GeometryData::ShapeFunctionsGradientsType N_der;
              geom.ShapeFunctionsIntegrationPointsGradients(N_der,GeometryData::GI_GAUSS_1);
              Matrix& DNDX = N_der[0];

              // Getting the already calculated values for the Normals
              //array_1d<double,3> Normal = (pElem)->FastGetSolutionStepValue(NORMAL);

              // Creating a variable where the interpolated pressure value is stored
              double pressure = 0.0;

              // getting the data of the solution step
              // Creating pointers to store the Solution
              //array_1d<double,3>& tractions = (pElem)->FastGetSolutionStepValue(traction);



              // Getting the pressure values needed for the interpolation
              for(int i=0; i<NumNodes; i++)
              {
                  pressure = pressure + NElem[i]*geom[i].FastGetSolutionStepValue(rOriginVariablePress);
              }

              // Gradient_ij = Du_i/Dx_j
              Matrix Gradient = ZeroMatrix(TDim,TDim);
              for (int n = 0; n < NumNodes; n++)
              {
                  // Velocities vx of the nodes of the triangles
                  const array_1d<double,3> & v_node = geom[n].FastGetSolutionStepValue(velocity);
                  //KRATOS_WATCH(v_node)
                  for (std::size_t i = 0; i < TDim; i++)
                      for (std::size_t j = 0; j < TDim; j++)
                          Gradient(i,j) += DNDX(n,j)*v_node[i];
              }


              // Creating the symmetrix velocity gradients (gradv + gradv')
              Matrix matrix_der_sym = Gradient + boost::numeric::ublas::trans(Gradient);
              //KRATOS_WATCH(Gradient);
              //KRATOS_WATCH(matrix_der_sym);

              // Multiplication of the interpolated pressure value with the normal vector
              array_1d<double,3> press_norm(3,0.0);
              for(std::size_t i=0; i<TDim; i++)
              {
                  press_norm[i] = pressure * Normal[i];
              }

              // Coefficient of viscosity and density
              // Creating the interpolation for viscosity
              double viscosity = 0;
              double density = 0;

              for(int i=0; i<NumNodes; i++)
              {
                  viscosity = viscosity + NElem[i]*geom[i].FastGetSolutionStepValue(VISCOSITY);
                  density   = density  + NElem[i]*geom[i].FastGetSolutionStepValue(DENSITY);
              }

              // Multiplication with symmetric velocity gradient with normals and coefficient nu
              array_1d<double,3> vel_norm(3,0.0);
              for(std::size_t i=0; i<TDim; i++)
              {
                  for(std::size_t j=0;j<TDim;j++)
                  {
                      vel_norm[i] = vel_norm[i] + matrix_der_sym(i,j) * Normal[j];
                  }
              }

              // Calculate the traction forces

              for(std::size_t i=0; i<TDim; i++)
              {
                  traction[i] = -press_norm[i] + viscosity * density * vel_norm[i];
              }
/* Uncomment this to apply Neumann BC for the skew-symmetric formulation*/
              array_1d<double,3> Vgauss(3,0.0);
              for(int i=0; i<NumNodes; i++)
              {
                  Vgauss += NElem[i]*geom[i].FastGetSolutionStepValue(VELOCITY);
              }

              double prod = 0.0;
              for (std::size_t d = 0; d < TDim; d++)
                  prod += Normal[d]*Vgauss[d];

              traction -= 0.5*density*prod*Vgauss;
/**/

              //KRATOS_WATCH(traction);
          }

      //***************************************

      inline void Clear(ModelPart::NodesContainerType::iterator node_it,  int step_data_size )
          {
              std::size_t buffer_size = node_it->GetBufferSize();

              for(std::size_t step = 0; step<buffer_size; step++)
              {
                  //getting the data of the solution step
                  double* step_data = (node_it)->SolutionStepData().Data(step);

                  //copying this data in the position of the vector we are interested in
                  for(int j= 0; j< step_data_size; j++)
                  {
                      step_data[j] = 0.0;
                  }
              }
          }


      inline void ClearVariables(ModelPart::NodesContainerType::iterator node_it , Variable<array_1d<double,3> >& rVariable)
          {
              array_1d<double, 3>& Aux_var = node_it->FastGetSolutionStepValue(rVariable, 0);

              noalias(Aux_var) = ZeroVector(3);

          }


      inline void ClearVariables(ModelPart::NodesContainerType::iterator node_it,  Variable<double>& rVariable)
          {
              double& Aux_var = node_it->FastGetSolutionStepValue(rVariable, 0);

              Aux_var = 0.0;

          }



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
      InterpolationUtility& operator=(InterpolationUtility const& rOther);

      /// Copy constructor.
      InterpolationUtility(InterpolationUtility const& rOther);


      ///@}

    }; // Class InterpolationUtility

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

//  /// input stream function
//  inline std::istream& operator >> (std::istream& rIStream,
//				    InterpolationUtility& rThis);

  /// output stream function
  template<std::size_t TDim>
  inline std::ostream& operator << (std::ostream& rOStream,
                                    const InterpolationUtility<TDim>& rThis)
    {
      rThis.PrintInfo(rOStream);
      rOStream << std::endl;
      rThis.PrintData(rOStream);

      return rOStream;
    }
  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_INTERPOLATION_UTILITY_H
