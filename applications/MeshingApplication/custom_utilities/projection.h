// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Antonia Larese de Tetto
//

#if !defined(KRATOS_PROJECTION )
#define  KRATOS_PROJECTION

// External includes

// System includes
#include <string>
#include <iostream>
#include <stdlib.h>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "utilities/timer.h"
#include "geometries/triangle_2d_3.h"
#include "meshing_application_variables.h"
#include "spatial_containers/spatial_containers.h"

namespace Kratos
{

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
        return dist; //square distance because it is easier to work without the square root//
    }
};



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

/// This class allows the interpolation between non-matching meshes in 2D and 3D.
/** @author  Antonia Larese De Tetto <antoldt@cimne.upc.edu>
*
* This class allows the interpolation of a variable or of the whole model parte between non-matching meshes
* in 2D and 3D.
*
* For every node of the destination model part it is checked in which element of the origin model part it is
* contained and a linear interpolation is performed
*
* The data structure used by default is kd tree although bin, kd-tree of bins can be easily used just commenting
* and decommenting the opportune lines at the beginning of the function
*/

//class MeshTransfer
template<std::size_t TDim >
class MeshTransfer
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MeshTransfer
    KRATOS_CLASS_POINTER_DEFINITION(MeshTransfer<TDim >);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MeshTransfer() = default; //

    /// Destructor.
    virtual ~MeshTransfer() = default;


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    //If you want to pass the whole model part
    //**********************************************************************
    //**********************************************************************
    /// Interpolate the whole problem type
    /**
      * @param rOrigin_ModelPart: the model part  all the variable should be taken from
      * @param rDestination_ModelPart: the destination model part where we want to know the values of the variables
      */
    void DirectInterpolation(
        ModelPart& rOrigin_ModelPart ,
        ModelPart& rDestination_ModelPart
    )
    {
        KRATOS_TRY

        //**********************************************************************
// 			numofpts_origin = rOrigin_ModelPart.Nodes().size();
// 			numofpts_destination = rDestination_ModelPart.Nodes().size();
//
        //*******************************************************************
        //properties to be used in the generation
        Properties::Pointer properties = rDestination_ModelPart.GetMesh().pGetProperties(1);

        //defintions for spatial search
        typedef Node<3> PointType;
        typedef Node<3>::Pointer PointTypePointer;
        typedef std::vector<PointType::Pointer>           PointVector;
        typedef std::vector<PointType::Pointer>::iterator PointIterator;
        typedef std::vector<double>               DistanceVector;
        typedef std::vector<double>::iterator     DistanceIterator;


        //creating an auxiliary list for the new nodes
        PointVector list_of_new_nodes;

// 			KRATOS_WATCH("STARTING KDTREE CONSTRUCTION");
        //starting calculating time of construction of the kdtree
        auto inital_time = std::chrono::steady_clock::now();

        //*************
        // Bucket types
        typedef Bucket< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
// 			   typedef Bins< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > StaticBins;
// 			   typedef BinsDynamic< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > DynamicBins;
        //*************
        // DynamicBins;

        typedef Tree< KDTreePartition<BucketType> > tree; 		//Kdtree;
// 			   typedef Tree< OCTreePartition<BucketType> > tree; 		//Octree;
// 			   typedef Tree< StaticBins > tree; 		     		//Binstree;
// 			   typedef Tree< KDTreePartition<StaticBins> > tree; 		//KdtreeBins;
// 			   typedef typename KdtreeBins::Partitions SubPartitions;
// 			   typedef Tree< OCTreePartition<StaticBins> > tree; 		//OctreeBins;
        /*
        			   typedef Bins< TDim, PointType, stdPointVector> stdBins;
        			   typedef Tree< Bins<TDim,PointType,stdPointVector> > tree; 	//stdStaticBins;*/


        for(ModelPart::NodesContainerType::iterator node_it = rDestination_ModelPart.NodesBegin();
                node_it != rDestination_ModelPart.NodesEnd(); ++node_it)
        {
            //PointType::Pointer pnode(new PointType(*node_it));
            Node<3>::Pointer pnode = *(node_it.base());

            //putting the nodes of the destination_model part in an auxiliary list
            list_of_new_nodes.push_back( pnode );
        }

        std::cout << "kdt construction time " << Timer::ElapsedSeconds(inital_time) << std::endl;
        //finishing calculating time of construction of the kdtree
// 			KRATOS_WATCH("FINISHING KDTREE CONSTRUCTION");

        //create a spatial database with the list of new nodes
        unsigned int bucket_size = 20;
        tree nodes_tree(list_of_new_nodes.begin(),list_of_new_nodes.end(),bucket_size);


        //work arrays
        Node<3> work_point(0,0.0,0.0,0.0);
        unsigned int MaximumNumberOfResults = 10000;
        PointVector Results(MaximumNumberOfResults);
        DistanceVector ResultsDistances(MaximumNumberOfResults);
        array_1d<double,TDim+1> N; //Shape functions vector//
        int step_data_size = rDestination_ModelPart.GetNodalSolutionStepDataSize();

        for(ModelPart::NodesContainerType::iterator node_it = rDestination_ModelPart.NodesBegin();
                node_it != rDestination_ModelPart.NodesEnd(); ++node_it)
        {
            //Setting to zero the whole model part
            Clear(node_it,  step_data_size );
        }
        inital_time = std::chrono::steady_clock::now();
        //loop over all of the elements in the "old" list to perform the interpolation
        for( ModelPart::ElementsContainerType::iterator el_it = rOrigin_ModelPart.ElementsBegin();
                el_it != rOrigin_ModelPart.ElementsEnd(); el_it++)
        {
            Geometry<Node<3> >&geom = el_it->GetGeometry();

            //find the center and "radius" of the element
// 				double xc,  yc,  radius;
// 				if(TDim == 2)
// 				{
// 					 CalculateCenterAndSearchRadius( geom[0].X(), geom[0].Y(),
// 								geom[1].X(), geom[1].Y(),
// 								geom[2].X(), geom[2].Y(),
// 								xc,yc,radius);
// 	 				work_point.X() = xc; work_point.Y() = yc;
//
// 				}
// 				else
// 				{
// 					 double zc;
// 					 CalculateCenterAndSearchRadius( geom[0].X(), geom[0].Y(),geom[0].Z(),
// 								geom[1].X(), geom[1].Y(),geom[1].Z(),
// 								geom[2].X(), geom[2].Y(),geom[2].Z(),
// 								geom[3].X(), geom[3].Y(),geom[3].Z(),
// 								xc,yc,zc, radius);
// 					 work_point.X() = xc; work_point.Y() = yc; work_point.Z() = zc;
//
// 				}
            double xc, yc, zc,  radius;
            CalculateCenterAndSearchRadius( geom, xc,yc,zc, radius,N);
            work_point.X() = xc;
            work_point.Y() = yc;
            work_point.Z() = zc;

            //find all of the new nodes within the radius
            int number_of_points_in_radius;

            //look between the new nodes which of them is inside the radius of the circumscribed cyrcle
            number_of_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, Results.begin(),
                                         ResultsDistances.begin(),  MaximumNumberOfResults);


            //check if inside
            for( auto it_found = Results.begin(); it_found != Results.begin() + number_of_points_in_radius; it_found++)
            {

                bool is_inside = false;
                //once we are sure the node in inside the circle we have to see if it is inside the triangle i.e. if all the Origin element shape functions are >1
// 				        if(TDim == 2)
// 				        {
// 					       is_inside = CalculatePosition(geom[0].X(), geom[0].Y(),
// 									       geom[1].X(), geom[1].Y(),
// 									       geom[2].X(), geom[2].Y(),
// 									       (*it_found)->X(),(*it_found)->Y(),N);
//
// 				        }
// 				        else
// 				        {
// 					       is_inside = CalculatePosition(geom[0].X(), geom[0].Y(), geom[0].Z(),
// 									       geom[1].X(), geom[1].Y(), geom[1].Z(),
// 									       geom[2].X(), geom[2].Y(), geom[2].Z(),
// 									       geom[3].X(), geom[3].Y(), geom[3].Z(),
// 									       (*it_found)->X(),(*it_found)->Y(),(*it_found)->Z(),N);
//
// 				        }
                double is_visited = (*it_found)->GetValue(IS_VISITED);
                is_inside = CalculatePosition(geom,	(*it_found)->X(),(*it_found)->Y(),(*it_found)->Z(),N);

                //if the node falls inside the element interpolate
                if(is_inside == true && is_visited != 1.0)
                {
                    //Interpolating all the rVariables of the rOrigin_ModelPart to get their nodal value in the rDestination_ModelPart
                    Interpolate(  el_it,  N, step_data_size, *it_found );

                }
            }

        }
        std::cout << "search and interpolation time " << Timer::ElapsedSeconds(inital_time) << std::endl;
        KRATOS_CATCH("")
    }


    //If you want to pass only one variable
    //**********************************************************************
    //**********************************************************************
    /// Interpolate one variable
    /**
      * @param rOrigin_ModelPart: the model part  all the variable should be taken from
      * @param rDestination_ModelPart: the destination model part where we want to know the values of the variables
      * @param rOriginVariable: the name of the interpolated variable in the origin model part
      * @param rOriginVariable: the name of the interpolated variable in the destination model part
      */
    template<class TDataType>
    void DirectVariableInterpolation(
        ModelPart& rOrigin_ModelPart ,
        ModelPart& rDestination_ModelPart,
        Variable<TDataType>& rOriginVariable ,
        Variable<TDataType>& rDestinationVariable
    )
    {

        KRATOS_TRY

        //*******************************************************************
        //properties to be used in the generation
        Properties::Pointer properties = rDestination_ModelPart.GetMesh().pGetProperties(1);

        //defintions for spatial search
        typedef Node<3> PointType;
        typedef Node<3>::Pointer PointTypePointer;
        typedef std::vector<PointType::Pointer>           PointVector;
        typedef std::vector<PointType::Pointer>::iterator PointIterator;
        typedef std::vector<double>               DistanceVector;
        typedef std::vector<double>::iterator     DistanceIterator;


        //creating an auxiliary list for the new nodes
        PointVector list_of_new_nodes;

// 			KRATOS_WATCH("STARTING KDTREE CONSTRUCTION");
        //starting calculating time of construction of the kdtree
        auto inital_time = std::chrono::steady_clock::now();

        //*************
        // Bucket types
        typedef Bucket< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
// 			   typedef Bins< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > StaticBins;
// 			   typedef BinsDynamic< TDim, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > DynamicBins;
        //*************
        // DynamicBins;

        typedef Tree< KDTreePartition<BucketType> > tree; 		//Kdtree;
// 			   typedef Tree< OCTreePartition<BucketType> > tree; 		//Octree;
// 			   typedef Tree< StaticBins > tree; 		     		//Binstree;
// 			   typedef Tree< KDTreePartition<StaticBins> > tree; 		//KdtreeBins;
// 			   typedef typename KdtreeBins::Partitions SubPartitions;
// 			   typedef Tree< OCTreePartition<StaticBins> > tree; 		//OctreeBins;
        /*
        			   typedef Bins< TDim, PointType, stdPointVector> stdBins;
        			   typedef Tree< Bins<TDim,PointType,stdPointVector> > tree; 	//stdStaticBins;*/

        for(ModelPart::NodesContainerType::iterator node_it = rDestination_ModelPart.NodesBegin();
                node_it != rDestination_ModelPart.NodesEnd(); ++node_it)
        {
            node_it->GetValue(IS_VISITED) = 0.0;
        }

// KRATOS_WATCH("line 328")
        for(ModelPart::NodesContainerType::iterator node_it = rDestination_ModelPart.NodesBegin();
                node_it != rDestination_ModelPart.NodesEnd(); ++node_it)
        {

            ClearVariables(node_it, rDestinationVariable);


            //PointType::Pointer pnode(new PointType(*node_it));
            Node<3>::Pointer pnode = *(node_it.base());

            //putting the nodes of the destination_model part in an auxiliary list
            list_of_new_nodes.push_back( pnode );


        }

        std::cout << "kdt construction time " << Timer::ElapsedSeconds(inital_time) << std::endl;
        //finishing calculating time of construction of the kdtree
// 			KRATOS_WATCH("FINISHING KDTREE CONSTRUCTION");

        //create a spatial database with the list of new nodes
        unsigned int bucket_size = 20;
        tree nodes_tree(list_of_new_nodes.begin(),list_of_new_nodes.end(),bucket_size);

        //work arrays
        Node<3> work_point(0,0.0,0.0,0.0);
        unsigned int MaximumNumberOfResults = 10000;
        PointVector Results(MaximumNumberOfResults);
        DistanceVector ResultsDistances(MaximumNumberOfResults);
        array_1d<double,TDim+1> N; //Shape functions vector//
        //int step_data_size = rDestination_ModelPart.GetNodalSolutionStepDataSize();
        //unsigned int TDim = 3;
// KRATOS_WATCH("line 359")
        inital_time = std::chrono::steady_clock::now();
        //loop over all of the elements in the "old" list to perform the interpolation
        for( ModelPart::ElementsContainerType::iterator el_it = rOrigin_ModelPart.ElementsBegin();
                el_it != rOrigin_ModelPart.ElementsEnd(); el_it++)
        {
            Geometry<Node<3> >&geom = el_it->GetGeometry();

            //find the center and "radius" of the element
// 				double xc,  yc,  radius;
// 				if(TDim == 2)
// 				{
// 					 CalculateCenterAndSearchRadius( geom[0].X(), geom[0].Y(),
// 								geom[1].X(), geom[1].Y(),
// 								geom[2].X(), geom[2].Y(),
// 								xc,yc,radius);
// 	 				 work_point.X() = xc; work_point.Y() = yc;
//
// 				}
// 				else
// 				{
// 					 double zc;
// 					 CalculateCenterAndSearchRadius( geom[0].X(), geom[0].Y(),geom[0].Z(),
// 								geom[1].X(), geom[1].Y(),geom[1].Z(),
// 								geom[2].X(), geom[2].Y(),geom[2].Z(),
// 								geom[3].X(), geom[3].Y(),geom[3].Z(),
// 								xc,yc,zc, radius);
// 					 work_point.X() = xc; work_point.Y() = yc; work_point.Z() = zc;
//
// 				}
            double xc, yc, zc,  radius;
            CalculateCenterAndSearchRadius( geom, xc,yc,zc, radius, N);
            work_point.X() = xc;
            work_point.Y() = yc;
            work_point.Z() = zc;
// KRATOS_WATCH("line 391")
            //find all of the new nodes within the radius
            int number_of_points_in_radius;

            //look between the new nodes which of them is inside the radius of the circumscribed circle
            number_of_points_in_radius = nodes_tree.SearchInRadius(work_point, radius, Results.begin(),
                                         ResultsDistances.begin(),  MaximumNumberOfResults);


            //check if inside
            for( auto it_found = Results.begin(); it_found != Results.begin() + number_of_points_in_radius; it_found++)
            {
// KRATOS_WATCH("line 402")
                bool is_inside = false;
                //once we are sure the node in inside the circle we have to see if it is inside the triangle i.e. if all the Origin element shape functions are >1
// 					 if(TDim == 2)
// 					 {
// 						is_inside = CalculatePosition(geom[0].X(), geom[0].Y(),
// 										geom[1].X(), geom[1].Y(),
// 										geom[2].X(), geom[2].Y(),
// 										(*it_found)->X(),(*it_found)->Y(),N);
//
// 					 }
// 					 else
// 					 {
// 						is_inside = CalculatePosition(geom[0].X(), geom[0].Y(), geom[0].Z(),
// 										geom[1].X(), geom[1].Y(), geom[1].Z(),
// 										geom[2].X(), geom[2].Y(), geom[2].Z(),
// 										geom[3].X(), geom[3].Y(), geom[3].Z(),
// 										(*it_found)->X(),(*it_found)->Y(),(*it_found)->Z(),N);
//
// 					 }

                double is_visited = (*it_found)->GetValue(IS_VISITED);
                is_inside = CalculatePosition(geom,	(*it_found)->X(),(*it_found)->Y(),(*it_found)->Z(),N);
// KRATOS_WATCH("line 423")
// KRATOS_WATCH("IS INSIDE")
// KRATOS_WATCH((*it_found)->Id())
                //if the node falls inside the element interpolate
                if(is_inside == true && is_visited != 1.0)
                {
                    //CANCELLA insert the variable TDim
                    //Interpolating all the rOriginVariables to get their nodal value in the rDestination_ModelPart
                    Interpolate(  el_it,  N, *it_found , rOriginVariable , rDestinationVariable  );

                }
            }

        }
        std::cout << "search and interpolation time " << Timer::ElapsedSeconds(inital_time) << std::endl;
        KRATOS_CATCH("")
    }





    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a stemplate<class T, std::size_t dim> tring.
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
    ///@name Protected static Member rVariables
    ///@{


    ///@}
    ///@name Protected member rVariables
    ///@{ template<class T, std::size_t dim>


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
    ///@name Static Member rVariables
    ///@{


    ///@}
    ///@name Member rVariables
    ///@{



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
    inline double CalculateVol(	const double x0, const double y0,
                                const double x1, const double y1,
                                const double x2, const double y2
                              )
    {
        return 0.5*( (x1-x0)*(y2-y0)- (y1-y0)*(x2-x0) );
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
            KRATOS_WATCH("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
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
//el_it		     	Element iterator
//N			Shape functions
//step_data_size
//pnode			pointer to the node

    //projecting total model part 2Dversion
    void Interpolate(
        ModelPart::ElementsContainerType::iterator el_it,
        const array_1d<double,3>& N,
        int step_data_size,
        Node<3>::Pointer pnode)
    {
        //Geometry element of the rOrigin_ModelPart
        Geometry< Node<3> >& geom = el_it->GetGeometry();

        unsigned int buffer_size = pnode->GetBufferSize();

        for(unsigned int step = 0; step<buffer_size; step++)
        {
            //getting the data of the solution step
            double* step_data = (pnode)->SolutionStepData().Data(step);

            double* node0_data = geom[0].SolutionStepData().Data(step);
            double* node1_data = geom[1].SolutionStepData().Data(step);
            double* node2_data = geom[2].SolutionStepData().Data(step);

            //copying this data in the position of the vector we are interested in
            for(int j= 0; j< step_data_size; j++)
            {
                step_data[j] = N[0]*node0_data[j] + N[1]*node1_data[j] + N[2]*node2_data[j];
            }
        }
        pnode->GetValue(IS_VISITED) = 1.0;

    }
    //projecting total model part 3Dversion
    void Interpolate(
        ModelPart::ElementsContainerType::iterator el_it,
        const array_1d<double,4>& N,
        int step_data_size,
        Node<3>::Pointer pnode)//dimension or number of nodes???
    {
        //Geometry element of the rOrigin_ModelPart
        Geometry< Node<3> >& geom = el_it->GetGeometry();

        unsigned int buffer_size = pnode->GetBufferSize();

        for(unsigned int step = 0; step<buffer_size; step++)
        {
            //getting the data of the solution step
            double* step_data = (pnode)->SolutionStepData().Data(step);

            double* node0_data = geom[0].SolutionStepData().Data(step);
            double* node1_data = geom[1].SolutionStepData().Data(step);
            double* node2_data = geom[2].SolutionStepData().Data(step);
            double* node3_data = geom[3].SolutionStepData().Data(step);

            //copying this data in the position of the vector we are interested in
            for(int j= 0; j< step_data_size; j++)
            {
                step_data[j] = N[0]*node0_data[j] + N[1]*node1_data[j] + N[2]*node2_data[j] + N[3]*node3_data[j];
            }
        }
        pnode->GetValue(IS_VISITED) = 1.0;

    }


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

        unsigned int buffer_size = pnode->GetBufferSize();

        for(unsigned int step = 0; step<buffer_size; step++)
        {
            //getting the data of the solution step
            array_1d<double,3>& step_data = (pnode)->FastGetSolutionStepValue(rDestinationVariable , step);
            //Reference or no reference???//CANCELLA
            const array_1d<double,3>& node0_data = geom[0].FastGetSolutionStepValue(rOriginVariable , step);
            const array_1d<double,3>& node1_data = geom[1].FastGetSolutionStepValue(rOriginVariable , step);
            const array_1d<double,3>& node2_data = geom[2].FastGetSolutionStepValue(rOriginVariable , step);

            //copying this data in the position of the vector we are interested in
            for(unsigned int j= 0; j< TDim; j++)
            {
                step_data[j] = N[0]*node0_data[j] + N[1]*node1_data[j] + N[2]*node2_data[j];
            }
        }
        pnode->GetValue(IS_VISITED) = 1.0;

    }

    //projecting an array1D 3Dversion
    void Interpolate(
        ModelPart::ElementsContainerType::iterator el_it,
        const array_1d<double,4>& N,
        Node<3>::Pointer pnode,
        Variable<array_1d<double,3> >& rOriginVariable,
        Variable<array_1d<double,3> >& rDestinationVariable)

    {
        //Geometry element of the rOrigin_ModelPart
        Geometry< Node<3> >& geom = el_it->GetGeometry();

        unsigned int buffer_size = pnode->GetBufferSize();

        for(unsigned int step = 0; step<buffer_size; step++)
        {
            //getting the data of the solution step
            array_1d<double,3>& step_data = (pnode)->FastGetSolutionStepValue(rDestinationVariable , step);
            //Reference or no reference???//CANCELLA
            const array_1d<double,3>& node0_data = geom[0].FastGetSolutionStepValue(rOriginVariable , step);
            const array_1d<double,3>& node1_data = geom[1].FastGetSolutionStepValue(rOriginVariable , step);
            const array_1d<double,3>& node2_data = geom[2].FastGetSolutionStepValue(rOriginVariable , step);
            const array_1d<double,3>& node3_data = geom[3].FastGetSolutionStepValue(rOriginVariable , step);

            //copying this data in the position of the vector we are interested in
            for(unsigned int j= 0; j< TDim; j++)
            {
                step_data[j] = N[0]*node0_data[j] + N[1]*node1_data[j] + N[2]*node2_data[j] + N[3]*node3_data[j];
            }
        }
        pnode->GetValue(IS_VISITED) = 1.0;

    }
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

        unsigned int buffer_size = pnode->GetBufferSize();
        //facendo un loop sugli step temporali step_data come salva i dati al passo anteriore? Cioś dove passiamo l'informazione ai nodi???
        for(unsigned int step = 0; step<buffer_size; step++)
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
        pnode->GetValue(IS_VISITED) = 1.0;

    }
    //projecting a scalar 3Dversion
    void Interpolate(
        ModelPart::ElementsContainerType::iterator el_it,
        const array_1d<double,4>& N,
        Node<3>::Pointer pnode,
        Variable<double>& rOriginVariable,
        Variable<double>& rDestinationVariable)
    {
        //Geometry element of the rOrigin_ModelPart
        Geometry< Node<3> >& geom = el_it->GetGeometry();

        unsigned int buffer_size = pnode->GetBufferSize();
        //facendo un loop sugli step temporali step_data come salva i dati al passo anteriore? Cioś dove passiamo l'informazione ai nodi???
        for(unsigned int step = 0; step<buffer_size; step++)
        {
            //getting the data of the solution step
            double& step_data = (pnode)->FastGetSolutionStepValue(rDestinationVariable , step);
            //Reference or no reference???//CANCELLA
            const double node0_data = geom[0].FastGetSolutionStepValue(rOriginVariable , step);
            const double node1_data = geom[1].FastGetSolutionStepValue(rOriginVariable , step);
            const double node2_data = geom[2].FastGetSolutionStepValue(rOriginVariable , step);
            const double node3_data = geom[3].FastGetSolutionStepValue(rOriginVariable , step);

            //copying this data in the position of the vector we are interested in

            step_data = N[0]*node0_data + N[1]*node1_data + N[2]*node2_data + N[3]*node3_data;

        }
        pnode->GetValue(IS_VISITED) = 1.0;

    }
    inline void Clear(ModelPart::NodesContainerType::iterator node_it,  int step_data_size )
    {
        unsigned int buffer_size = node_it->GetBufferSize();

        for(unsigned int step = 0; step<buffer_size; step++)
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
    MeshTransfer& operator=(MeshTransfer const& rOther);


    ///@}

}; // Class MeshTransfer

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{




/// output stream function
template<std::size_t TDim>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MeshTransfer<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_PROJECTION  defined
