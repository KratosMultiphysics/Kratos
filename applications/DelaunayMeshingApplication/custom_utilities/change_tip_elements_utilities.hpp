//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined(KRATOS_CHANGE_TIP_ELEMENTS_UTILITIES_H_INCLUDED)
#define  KRATOS_CHANGE_TIP_ELEMENTS_UTILITIES_H_INCLUDED


// System includes
#include <stdlib.h>

// External includes
#include "utilities/timer.h"

// Project includes
#include "includes/variables.h"
#include "includes/model_part.h"
#include "geometries/triangle_2d_3.h"
#include "spatial_containers/spatial_containers.h"
#include "delaunay_meshing_application_variables.h"

#if !defined(KRATOS_TRIANGLE_EXTERNAL_H_INCLUDED)
#define  KRATOS_TRIANGLE_EXTERNAL_H_INCLUDED
#include "triangle.h"
#endif

namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

//TIP ELEMENTS : Elements with all its nodes in boundary

class ChangeTipElementsUtilities
{
public:
    /**@name Type Definitions */
    /*@{ */
    typedef ModelPart::ElementsContainerType                  ElementsContainerType;
    typedef ModelPart::NodesContainerType                        NodesContainerType;
    typedef ModelPart::MeshType::GeometryType::PointsArrayType      PointsArrayType;
    typedef ModelPart::MeshType::NodeType                                  NodeType;

    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */
    ChangeTipElementsUtilities(){};

    /** Destructor.
     */
    ~ChangeTipElementsUtilities(){};

    /** Operators.
     */


    //**************************************************************************
    //**************************************************************************

    bool SwapDiagonals(ModelPart& rModelPart)
    {
        KRATOS_TRY


        bool any_swap_performed = false;

	for(ModelPart::ElementsContainerType::const_iterator ie = rModelPart.ElementsBegin(); ie != rModelPart.ElementsEnd(); ++ie)
	  {

	    PointsArrayType& vertices= (ie)->GetGeometry().Points();

	    double boundary_nodes = 0;

	    int option=0;

	    unsigned int NumberOfVertices =vertices.size();
	    for(unsigned int i=0; i<NumberOfVertices; ++i)
	      {
		if(vertices[i].Is(BOUNDARY) && vertices[i].IsNot(TO_ERASE))
		  boundary_nodes+=1;
	      }

	    if(boundary_nodes == NumberOfVertices){


	      double xc1=0;
	      double yc1=0;

	      CalculateCenter( vertices[0].X(), vertices[0].Y(),
			       vertices[1].X(), vertices[1].Y(),
			       vertices[2].X(), vertices[2].Y(),
			       xc1,yc1);



	      double Area = 0;
	      double Area_Candidate = 0;

	      //1.- Set new neighbour elements to the tip element and shared element
	      WeakPointerVector<Element >& neighb_elems = ie->GetValue(NEIGHBOUR_ELEMENTS);

	      std::vector<unsigned int> new_vertices_1(3);
	      std::vector<unsigned int> new_vertices_2(3);

	      bool swap_performed = false;
	      for(WeakPointerVector< Element >::iterator ne = neighb_elems.begin(); ne!=neighb_elems.end(); ++ne)
		{

		  if (ne->Id() != ie->Id()){

		    //Geometry<Node<3> >& pGeom = (ne)->GetGeometry();
		    PointsArrayType& neighb_vertices= (ne)->GetGeometry().Points();

		    boundary_nodes=0;
		    for(unsigned int i=0; i<NumberOfVertices; ++i)
		      {
			if(neighb_vertices[i].Is(BOUNDARY))
			  boundary_nodes+=1;
		      }


		    if(boundary_nodes != NumberOfVertices){


		      double xc2=0;
		      double yc2=0;

		      CalculateCenter( neighb_vertices[0].X(), neighb_vertices[0].Y(),
				       neighb_vertices[1].X(), neighb_vertices[1].Y(),
				       neighb_vertices[2].X(), neighb_vertices[2].Y(),
				       xc2,yc2);


		      unsigned int tip_node_1 =10;
		      unsigned int tip_node_2 =10;
		      bool tip_1 = false;
		      bool tip_2 = false;
		      for(unsigned int i=0; i<NumberOfVertices; ++i)
			{
			  tip_1 = true;
			  tip_2 = true;
			  for(unsigned int j=0; j<NumberOfVertices; ++j)
			    {
			      if( vertices[i].Id() == neighb_vertices[j].Id() )
				tip_1 = false;

			      if( neighb_vertices[i].Id() == vertices[j].Id() )
				tip_2 = false;
			    }

			  if(tip_1)
			    tip_node_1=i;

			  if(tip_2)
			    tip_node_2=i;

			}

		      if(tip_node_1 ==10 || tip_node_2 ==10)
			std::cout<<" Shared Diagonal Not found in Change_Tip_Elements "<<std::endl;

		      unsigned int shared_node_1 =10;
		      unsigned int shared_node_2 =10;

		      for(unsigned int i=0; i<NumberOfVertices; ++i)
			{
			  if(tip_node_1 != i && shared_node_1 == 10)
			    shared_node_1 = i;

			  if(tip_node_1 != i && shared_node_1 != i && shared_node_2 == 10 )
			    shared_node_2 = i;
			}


		      //Compute sense of the triangle conectivities using the Area sign

		      Area = CalculateSignedArea( vertices[0].X(), vertices[0].Y(),
						  vertices[1].X(), vertices[1].Y(),
						  vertices[2].X(), vertices[2].Y() );

		      //Candidates:
		      Area_Candidate = CalculateSignedArea( vertices[tip_node_1].X(), vertices[tip_node_1].Y(),
							    neighb_vertices[tip_node_2].X(), neighb_vertices[tip_node_2].Y(),
							    vertices[shared_node_1].X(), vertices[shared_node_1].Y() );

		      double xc=0;
		      double yc=0;

		      CalculateCenter( vertices[tip_node_1].X(), vertices[tip_node_1].Y(),
				       neighb_vertices[tip_node_2].X(), neighb_vertices[tip_node_2].Y(),
				       vertices[shared_node_1].X(), vertices[shared_node_1].Y(),
				       xc,yc);


		      if(Area*Area_Candidate>0) { // A: tip_node_1, tip_node_2, shared_node_1

			new_vertices_1[0] = vertices[tip_node_1].Id();
			new_vertices_1[1] = neighb_vertices[tip_node_2].Id();
			new_vertices_1[2] = vertices[shared_node_1].Id();

			new_vertices_2[0] = neighb_vertices[tip_node_2].Id();
			new_vertices_2[1] = vertices[tip_node_1].Id();
			new_vertices_2[2] = vertices[shared_node_2].Id();

		      }
		      else{ 		 // B: tip_node_1, tip_node_2, shared_node_2

			new_vertices_1[0] = vertices[tip_node_1].Id();
			new_vertices_1[1] = neighb_vertices[tip_node_2].Id();
			new_vertices_1[2] = vertices[shared_node_2].Id();

			new_vertices_2[0] = neighb_vertices[tip_node_2].Id();
			new_vertices_2[1] = vertices[tip_node_1].Id();
			new_vertices_2[2] = vertices[shared_node_1].Id();

		      }

		      option = 0;
		      CalculateMinDistance (xc,yc,xc1,yc1,xc2,yc2,option);

		      // std::cout<<" BEFORE SWAP "<<std::endl;
		      // std::cout<<" Element i "<<ie->Id()<<" nodes ("<<vertices[0].Id()<<", "<<vertices[1].Id()<<", "<<vertices[2].Id()<<") "<<std::endl;

		      // std::cout<<" Element j "<<ne->Id()<<" nodes ("<<neighb_vertices[0].Id()<<", "<<neighb_vertices[1].Id()<<", "<<neighb_vertices[2].Id()<<") "<<std::endl;

		      //Change elements Conectivities preserving the closest center (gauss_point)
		      if(option==1){

			vertices(0) =  NodeType::Pointer(*(rModelPart.Nodes().find( new_vertices_1[0] )).base());
			vertices(1) =  NodeType::Pointer(*(rModelPart.Nodes().find( new_vertices_1[1] )).base());
			vertices(2) =  NodeType::Pointer(*(rModelPart.Nodes().find( new_vertices_1[2] )).base());

			neighb_vertices(0) = NodeType::Pointer(*(rModelPart.Nodes().find( new_vertices_2[0] )).base());
			neighb_vertices(1) = NodeType::Pointer(*(rModelPart.Nodes().find( new_vertices_2[1] )).base());
			neighb_vertices(2) = NodeType::Pointer(*(rModelPart.Nodes().find( new_vertices_2[2] )).base());

			// vertices        = new_vertices_1;
			// neighb_vertices = new_vertices_2;
			swap_performed = true;

		      // std::cout<<" AFTER SWAP 1 "<<std::endl;
		      // std::cout<<" Element i "<<ie->Id()<<" nodes ("<<vertices[0].Id()<<", "<<vertices[1].Id()<<", "<<vertices[2].Id()<<") "<<std::endl;

		      // std::cout<<" Element j "<<ne->Id()<<" nodes ("<<neighb_vertices[0].Id()<<", "<<neighb_vertices[1].Id()<<", "<<neighb_vertices[2].Id()<<") "<<std::endl;

		      }
		      else{

			vertices(0) = NodeType::Pointer(*(rModelPart.Nodes().find( new_vertices_2[0] )).base());
			vertices(1) = NodeType::Pointer(*(rModelPart.Nodes().find( new_vertices_2[1] )).base());
			vertices(2) = NodeType::Pointer(*(rModelPart.Nodes().find( new_vertices_2[2] )).base());

			neighb_vertices(0) = NodeType::Pointer(*(rModelPart.Nodes().find( new_vertices_1[0] )).base());
			neighb_vertices(1) = NodeType::Pointer(*(rModelPart.Nodes().find( new_vertices_1[1] )).base());
			neighb_vertices(2) = NodeType::Pointer(*(rModelPart.Nodes().find( new_vertices_1[2] )).base());

			// vertices        = new_vertices_2;
			// neighb_vertices = new_vertices_1;
			swap_performed  = true;

			// std::cout<<" AFTER SWAP 2 "<<std::endl;
			// std::cout<<" Element i "<<ie->Id()<<" nodes ("<<vertices[0].Id()<<", "<<vertices[1].Id()<<", "<<vertices[2].Id()<<") "<<std::endl;

			// std::cout<<" Element j "<<ne->Id()<<" nodes ("<<neighb_vertices[0].Id()<<", "<<neighb_vertices[1].Id()<<", "<<neighb_vertices[2].Id()<<") "<<std::endl;

		      }

		      if(option==0)
			std::cout<<" Something is wrong in the center calculation in Change_Tip_Elements "<<std::endl;

		    }

		  }


		}

	      if( swap_performed ){

		swap_performed = false;

		//2.- Set new neighbour elements to the tip  and shared elements
		int shared_side_1 =10;
		int shared_side_2 =10;


		for(WeakPointerVector< Element >::iterator ne = neighb_elems.begin(); ne!=neighb_elems.end(); ++ne)
		  {
		    if (ne->Id() != ie->Id() && !swap_performed){

		      PointsArrayType& neighb_vertices= (ne)->GetGeometry().Points();

		      WeakPointerVector<Element > neighb_neighb_elems = ne->GetValue(NEIGHBOUR_ELEMENTS);


		      int counter=0;
		      for(WeakPointerVector< Element >::iterator nne = neighb_neighb_elems.begin(); nne!=neighb_neighb_elems.end(); ++nne)
			{

			  if (ne->Id() != nne->Id() &&  ie->Id() != nne->Id() ){

			    PointsArrayType& neighb_neighb_vertices= (nne)->GetGeometry().Points();

			    // std::cout<<" Element k "<<ie->Id()<<" nodes ("<<neighb_neighb_vertices[0].Id()<<", "<<neighb_neighb_vertices[1].Id()<<", "<<neighb_neighb_vertices[2].Id()<<") "<<std::endl;

			    int shared_1 =0;
			    int shared_2 =0;

			    for(unsigned int i=0; i<NumberOfVertices; ++i)
			      {
				for(unsigned int j=0; j<NumberOfVertices; ++j)
				  {
				    if( vertices[i].Id() == neighb_neighb_vertices[j].Id() )
				      shared_1 +=1;

				    if( neighb_vertices[i].Id() == neighb_neighb_vertices[j].Id() )
				      shared_2 +=1;
				  }
			      }

			    if(shared_1==2)
			      shared_side_1=counter;

			    if(shared_2==2)
			      shared_side_2=counter;


			  }

			  counter++;

			}


		      if(shared_side_1 ==10 || shared_side_2 ==10 || shared_side_1 == shared_side_2 ){
			std::cout<<" Shared Elements Not found in Channge_Tip_Elements "<<std::endl;
		      }
		      else{

			// std::cout<<" Element_1 "<<ie->Id()<<" Element_2 "<<ne->Id()<<std::endl;
			// std::cout<<" shared_side_1 "<<shared_side_1<<" shared_side_2 "<<shared_side_2<<std::endl;
			// std::cout<<" Neighbours Selected ("<<neighb_neighb_elems[0].Id()<<", "<<neighb_neighb_elems[1].Id()<<", "<<neighb_neighb_elems[2].Id()<<") "<<std::endl;

			WeakPointerVector<Element >& nE_1 = ie->GetValue(NEIGHBOUR_ELEMENTS);
			// std::cout<<" Neighbours 1 pre("<<nE_1[0].Id()<<", "<<nE_1[1].Id()<<", "<<nE_1[2].Id()<<") "<<std::endl;

			if(option==2){

			  //the order is important because the iterator (ne) changes after the assignation
			  nE_1(0) = *(ie.base());
			  nE_1(1) = *(ne.base());
			  nE_1(2) = neighb_neighb_elems(shared_side_1);

			  // std::cout<<" Neighbours 1 post("<<nE_1[0].Id()<<", "<<nE_1[1].Id()<<", "<<nE_1[2].Id()<<") "<<std::endl;

			  WeakPointerVector<Element >& nE_2 = nE_1[0].GetValue(NEIGHBOUR_ELEMENTS);
			  // std::cout<<" Neighbours 2 pre("<<nE_2[0].Id()<<", "<<nE_2[1].Id()<<", "<<nE_2[2].Id()<<") "<<std::endl;

			  nE_2(0) = neighb_neighb_elems(shared_side_2);
			  nE_2(1) = nE_1(1);
			  nE_2(2) = *(ie.base());

			}
			else{

			  //the order is important because the iterator (ne) changes after the assignation
			  nE_1(0) = *(ne.base());
			  nE_1(1) = *(ie.base());
			  nE_1(2) = neighb_neighb_elems(shared_side_1);

			// std::cout<<" Neighbours 1 post("<<nE_1[0].Id()<<", "<<nE_1[1].Id()<<", "<<nE_1[2].Id()<<") "<<std::endl;

			  WeakPointerVector<Element >& nE_2 = nE_1[0].GetValue(NEIGHBOUR_ELEMENTS);
			  // std::cout<<" Neighbours 2 pre("<<nE_2[0].Id()<<", "<<nE_2[1].Id()<<", "<<nE_2[2].Id()<<") "<<std::endl;

			  nE_2(0) = nE_1(0);
			  nE_2(1) = *(ie.base());
			  nE_2(2) = neighb_neighb_elems(shared_side_2);


			}

			// std::cout<<" Neighbours 2 post("<<nE_2[0].Id()<<", "<<nE_2[1].Id()<<", "<<nE_2[2].Id()<<") "<<std::endl;
			WeakPointerVector<Element >& n_n_elems  = neighb_neighb_elems[shared_side_1].GetValue(NEIGHBOUR_ELEMENTS);
			for(WeakPointerVector< Element >::iterator nne = n_n_elems.begin(); nne!=n_n_elems.end(); ++nne)
			  {
			    if( nne->Id() == ne->Id() ){
			      n_n_elems(shared_side_1) = *(ie.base());
			    }
			  }


			swap_performed = true;
			any_swap_performed = true;

			std::cout<<" swap performed in element: "<<ie->Id()<<std::endl;

		      }

		    }

		  }

		//3.- Find condition asocied to this tip and change the Master Elements and Nodes
		if( swap_performed ){

		  for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(); ic!= rModelPart.ConditionsEnd(); ++ic)
		    {
		      Geometry< Node<3> >& rConditionGeom = ic->GetGeometry();
		      PointsArrayType& vertices_1= (ie)->GetGeometry().Points();                  //element1
		      PointsArrayType& vertices_2= ie->GetValue(NEIGHBOUR_ELEMENTS)[0].GetGeometry().Points();  //element2 in 0 position in neighbours

		      //Faces are vertices_1[0] and vertices_1[2] and vertices_2[1] and vertices_2[2] in option 1
		      //Faces are vertices_1[1] and vertices_1[2] and vertices_2[0] and vertices_2[2] in option 2

		      //option 1
		      int nda = 0;
		      int ndb = 1;

		      if(option==2)
			{
			  nda = 1;
			  ndb = 0;
			}



			if( (rConditionGeom[0].Id() == vertices_1[nda].Id()
			     && rConditionGeom[1].Id() == vertices_1[2].Id() ) ||
			    (rConditionGeom[0].Id() == vertices_1[2].Id()
			     && rConditionGeom[1].Id() == vertices_1[nda].Id() ) )
			  {
			    //usually one MasterElement and one MasterNode in 2D
			    ic->GetValue(MASTER_ELEMENTS)(0) = ( Element::WeakPointer( *(ie.base()) ) );
			    ic->GetValue(MASTER_NODES)(0) = (  Node<3>::WeakPointer( vertices_1(1) ) );

			  }


			if( (rConditionGeom[0].Id() == vertices_2[ndb].Id()
			     && rConditionGeom[1].Id() == vertices_2[2].Id() ) ||
			    (rConditionGeom[0].Id() == vertices_2[2].Id()
			     && rConditionGeom[1].Id() == vertices_2[ndb].Id() ) )
			  {
			    //usually one MasterElement and one MasterNode in 2D
			    ic->GetValue(MASTER_ELEMENTS)(0) = ( Element::WeakPointer( ie->GetValue(NEIGHBOUR_ELEMENTS)(0) ) );
			    ic->GetValue(MASTER_NODES)(0) = (  Node<3>::WeakPointer( vertices_2(0) ) );

			  }


			// if(ic->Is(CONTACT)){

			//   Element::ElementType& MasterElement = ic->GetValue(MASTER_ELEMENTS)[0];
			//   Element::NodeType&    MasterNode    = ic->GetValue(MASTER_NODES)[0];

			//   int  slave=-1;
			//   for(unsigned int i=0; i<MasterElement.GetGeometry().PointsNumber(); ++i)
			//     {
			//       if(MasterNode.Id()==MasterElement.GetGeometry()[i].Id())
			// 	{
			// 	  slave=i;
			// 	}
			//     }

			//   DenseMatrix<unsigned int> lpofa; //points that define the faces
			//   rConditionGeom.NodesInFaces(lpofa);



			//   if( (rConditionGeom[lpofa(1,slave)].Id() == vertices_1[nda].Id()
			//        && rConditionGeom[lpofa(2,slave)].Id() == vertices_1[2].Id() ) ||
			//       (rConditionGeom[lpofa(1,slave)].Id() == vertices_1[2].Id()
			//        && rConditionGeom[lpofa(2,slave)].Id() == vertices_1[nda].Id() ) )
			//     {
			//       //usually one MasterElement and one MasterNode in 2D
			//       ic->GetValue(MASTER_ELEMENTS)(0) = ( Element::WeakPointer( *(ie.base()) ) );
			//       ic->GetValue(MASTER_NODES)(0) = (  Node<3>::WeakPointer( vertices_1(1) ) );

			//     }


			//   if( (rConditionGeom[lpofa(1,slave)].Id() == vertices_2[ndb].Id()
			//        && rConditionGeom[lpofa(2,slave)].Id() == vertices_2[2].Id() ) ||
			//       (rConditionGeom[lpofa(1,slave)].Id() == vertices_2[2].Id()
			//        && rConditionGeom[lpofa(2,slave)].Id() == vertices_2[ndb].Id() ) )
			//     {
			//       //usually one MasterElement and one MasterNode in 2D
			//       ic->GetValue(MASTER_ELEMENTS)(0) = ( Element::WeakPointer( ie->GetValue(NEIGHBOUR_ELEMENTS)(0) ) );
			//       ic->GetValue(MASTER_NODES)(0) = (  Node<3>::WeakPointer( vertices_2(0) ) );

			//     }


			// }



		    }
		}

	      }
	      else{

		std::cout<<" swap NOT performed in element: "<<ie->Id()<<std::endl;

	      }



	    }

	    //Nodal Neighbours will be set later in the general search:
	    //The Diagonal swapping must be done before doing it

	  }

	return any_swap_performed;

	KRATOS_CATCH( "" )

    };


    // Swap Diagonals using the mesher output

  bool SwapDiagonals(ModelPart& rModelPart, struct triangulateio &out, std::vector<int>& PreservedElements)
    {
        KRATOS_TRY


        bool any_swap_performed = false;
	ModelPart::NodesContainerType::iterator nodes_begin = rModelPart.NodesBegin();

	for(int el = 0; el<out.numberoftriangles; ++el)
	  {

	    if(PreservedElements[el])
	    {

	    int option=0;
	    double boundary_nodes = 0;

	    int NumberOfVertices = 3;
	    int NumberOfNeighbourElements = 3;

	    array_1d<double, 3> Normal;
	    Geometry<Node<3> > vertices;
	    for(int pn=0; pn<NumberOfVertices; ++pn)
	      {
		vertices.push_back(*(nodes_begin + out.trianglelist[el*3+pn]-1).base());
		if(vertices[pn].Is(BOUNDARY) && vertices[pn].IsNot(TO_ERASE)){
		  Normal=vertices[pn].FastGetSolutionStepValue(NORMAL);
		  if(norm_2(Normal))
		    boundary_nodes+=1;
		  else
		    std::cout<<" boundary node but not Normal assigned "<<vertices[pn].Id()<<std::endl;
		}
	      }

	    if( boundary_nodes == NumberOfVertices ){

	      double xc1=0;
	      double yc1=0;

	      CalculateCenter( vertices[0].X(), vertices[0].Y(),
			       vertices[1].X(), vertices[1].Y(),
			       vertices[2].X(), vertices[2].Y(),
			       xc1,yc1);

	      double Area = 0;
	      double Area_Candidate = 0;

	      //1.- Set new neighbour elements to the tip element and shared element

	      std::vector<unsigned int> new_vertices_1(3);
	      std::vector<unsigned int> new_vertices_2(3);

	      bool swap_performed = false;

	      //number of vertices is equal to the number of elements in 2D
	      for(int ne=0; ne<NumberOfNeighbourElements; ++ne)
		{
		  int n_el = out.neighborlist[el*3+ne];

		  if( n_el>0 ){

		    if( PreservedElements[n_el-1] ){

		      boundary_nodes=0;
		      Geometry<Node<3> > neighb_vertices;
		      for(int pn=0; pn<NumberOfVertices; ++pn)
			{
			  neighb_vertices.push_back(*(nodes_begin + out.trianglelist[(n_el-1)*3+pn]-1).base());
			  if(neighb_vertices[pn].Is(BOUNDARY))
			    boundary_nodes+=1;
			}


		      if(boundary_nodes != NumberOfVertices){

			double xc2=0;
			double yc2=0;

			CalculateCenter( neighb_vertices[0].X(), neighb_vertices[0].Y(),
					 neighb_vertices[1].X(), neighb_vertices[1].Y(),
					 neighb_vertices[2].X(), neighb_vertices[2].Y(),
					 xc2,yc2);


			int tip_node_1 =10;
			int tip_node_2 =10;
			bool tip_1 = false;
			bool tip_2 = false;
			for(int i=0; i<NumberOfVertices; ++i)
			  {
			    tip_1 = true;
			    tip_2 = true;
			    for(int j=0; j<NumberOfVertices; ++j)
			      {
				if( vertices[i].Id() == neighb_vertices[j].Id() )
				  tip_1 = false;

				if( neighb_vertices[i].Id() == vertices[j].Id() )
				  tip_2 = false;
			      }

			    if(tip_1)
			      tip_node_1=i;

			    if(tip_2)
			      tip_node_2=i;

			  }

			if(tip_node_1 ==10 || tip_node_2 ==10)
			  std::cout<<" Shared Diagonal Not found in Change_Tip_Elements "<<std::endl;

			int shared_node_1 =10;
			int shared_node_2 =10;

			for(int i=0; i<NumberOfVertices; ++i)
			  {
			    if(tip_node_1 != i && shared_node_1 == 10)
			      shared_node_1 = i;

			    if(tip_node_1 != i && shared_node_1 != i && shared_node_2 == 10 )
			      shared_node_2 = i;
			  }


			//Compute sense of the triangle conectivities using the Area sign

			Area = CalculateSignedArea( vertices[0].X(), vertices[0].Y(),
						    vertices[1].X(), vertices[1].Y(),
						    vertices[2].X(), vertices[2].Y() );

			//Candidates:
			Area_Candidate = CalculateSignedArea( vertices[tip_node_1].X(), vertices[tip_node_1].Y(),
							      neighb_vertices[tip_node_2].X(), neighb_vertices[tip_node_2].Y(),
							      vertices[shared_node_1].X(), vertices[shared_node_1].Y() );

			double xc=0;
			double yc=0;

			CalculateCenter( vertices[tip_node_1].X(), vertices[tip_node_1].Y(),
					 neighb_vertices[tip_node_2].X(), neighb_vertices[tip_node_2].Y(),
					 vertices[shared_node_1].X(), vertices[shared_node_1].Y(),
					 xc,yc);


			if(Area*Area_Candidate>0) { // A: tip_node_1, tip_node_2, shared_node_1

			  new_vertices_1[0] = vertices[tip_node_1].Id();
			  new_vertices_1[1] = neighb_vertices[tip_node_2].Id();
			  new_vertices_1[2] = vertices[shared_node_1].Id();

			  new_vertices_2[0] = neighb_vertices[tip_node_2].Id();
			  new_vertices_2[1] = vertices[tip_node_1].Id();
			  new_vertices_2[2] = vertices[shared_node_2].Id();

			}
			else{ 		 // B: tip_node_1, tip_node_2, shared_node_2

			  new_vertices_1[0] = vertices[tip_node_1].Id();
			  new_vertices_1[1] = neighb_vertices[tip_node_2].Id();
			  new_vertices_1[2] = vertices[shared_node_2].Id();

			  new_vertices_2[0] = neighb_vertices[tip_node_2].Id();
			  new_vertices_2[1] = vertices[tip_node_1].Id();
			  new_vertices_2[2] = vertices[shared_node_1].Id();

			}

			option = 0;
			CalculateMinDistance (xc,yc,xc1,yc1,xc2,yc2,option);

			// std::cout<<" BEFORE SWAP "<<std::endl;
			// std::cout<<" Element i "<<ie->Id()<<" nodes ("<<vertices[0].Id()<<", "<<vertices[1].Id()<<", "<<vertices[2].Id()<<") "<<std::endl;

			// std::cout<<" Element j "<<ne->Id()<<" nodes ("<<neighb_vertices[0].Id()<<", "<<neighb_vertices[1].Id()<<", "<<neighb_vertices[2].Id()<<") "<<std::endl;

			//Change elements Conectivities preserving the closest center (gauss_point)

			  // std::cout<<" Pre triangle 1: "<<el+1<<"("<<out.trianglelist[(el)*3]<<", "<<out.trianglelist[(el)*3+1]<<", "<<out.trianglelist[(el)*3+2]<<")"<<std::endl;
			  // std::cout<<" Pre triangle 2: "<<n_el<<"("<<out.trianglelist[(n_el-1)*3]<<", "<<out.trianglelist[(n_el-1)*3+1]<<", "<<out.trianglelist[(n_el-1)*3+2]<<")"<<std::endl;

			if(option==1){

			  // std::cout<<" AFTER SWAP 1 ("<<out.trianglelist[(el)*3]<<","<<out.trianglelist[(el)*3+1]<<","<<out.trianglelist[(el)*3+2]<<") ("<<out.trianglelist[(n_el-1)*3]<<","<<out.trianglelist[(n_el-1)*3+1]<<","<<out.trianglelist[(n_el-1)*3+2]<<")"<<std::endl;

			  out.trianglelist[(el)*3]   = new_vertices_1[0];
			  out.trianglelist[(el)*3+1] = new_vertices_1[1];
			  out.trianglelist[(el)*3+2] = new_vertices_1[2];

			  out.trianglelist[(n_el-1)*3]   = new_vertices_2[0];
			  out.trianglelist[(n_el-1)*3+1] = new_vertices_2[1];
			  out.trianglelist[(n_el-1)*3+2] = new_vertices_2[2];

			  swap_performed = true;



			}
			else{


			  // std::cout<<" AFTER SWAP 1 ("<<out.trianglelist[(el)*3]<<","<<out.trianglelist[(el)*3+1]<<","<<out.trianglelist[(el)*3+2]<<") ("<<out.trianglelist[(n_el-1)*3]<<","<<out.trianglelist[(n_el-1)*3+1]<<","<<out.trianglelist[(n_el-1)*3+2]<<")"<<std::endl;

			  out.trianglelist[(el)*3]   = new_vertices_2[0];
			  out.trianglelist[(el)*3+1] = new_vertices_2[1];
			  out.trianglelist[(el)*3+2] = new_vertices_2[2];

			  out.trianglelist[(n_el-1)*3]   = new_vertices_1[0];
			  out.trianglelist[(n_el-1)*3+1] = new_vertices_1[1];
			  out.trianglelist[(n_el-1)*3+2] = new_vertices_1[2];

			  swap_performed  = true;





			}


			// std::cout<<" Tip Element: "<<el+1<<std::endl;
			// std::cout<<"(-";
			// for(int pn=0; pn<NumberOfVertices; ++pn)
			//   {
			//     std::cout<<(nodes_begin + out.trianglelist[el*3+pn]-1)->Id()<<"-";
			//   }
			// std::cout<<")"<<std::endl;

			// std::cout<<" Swapt Element: "<<n_el<<std::endl;
			// std::cout<<"(-";
			// for(int pn=0; pn<NumberOfVertices; ++pn)
			//   {
			//     std::cout<<(nodes_begin + out.trianglelist[(n_el-1)*3+pn]-1)->Id()<<"-";
			//   }
			//   std::cout<<")"<<std::endl;

			// std::cout<<" Post triangle 1: "<<el+1<<"("<<out.trianglelist[(el)*3]<<", "<<out.trianglelist[(el)*3+1]<<", "<<out.trianglelist[(el)*3+2]<<")"<<std::endl;
			// std::cout<<" Post triangle 2: "<<n_el<<"("<<out.trianglelist[(n_el-1)*3]<<", "<<out.trianglelist[(n_el-1)*3+1]<<", "<<out.trianglelist[(n_el-1)*3+2]<<")"<<std::endl;


			if(option==0)
			  std::cout<<" Something is wrong in the center calculation in Change_Tip_Elements "<<std::endl;

		      }
		    }
		  }


		}

	      if( swap_performed ){

		swap_performed = false;

		//2.- Set new neighbour elements to the tip and shared elements
		int shared_side_1 =10;
		int shared_side_2 =10;


		vertices.clear();
		for(int pn=0; pn<NumberOfVertices; ++pn)
		  {
		    vertices.push_back(*(nodes_begin + out.trianglelist[el*3+pn]-1).base());
		  }


		for(int ne=0; ne<NumberOfNeighbourElements; ++ne)
		  {
		    int n_el = out.neighborlist[el*3+ne];

		    if( n_el>0 && !swap_performed )
		      {

			if(PreservedElements[n_el-1]){

			  boundary_nodes=0;
			  Geometry<Node<3> > neighb_vertices;
			  for(int pn=0; pn<NumberOfVertices; ++pn)
			    {
			      neighb_vertices.push_back(*(nodes_begin + out.trianglelist[(n_el-1)*3+pn]-1).base());
			      if(neighb_vertices[pn].Is(BOUNDARY))
				boundary_nodes+=1;
			    }

			  if(boundary_nodes != NumberOfVertices){


			    int counter=0;
			    int sn_el=0;
			    for(int nne=0; nne<NumberOfNeighbourElements; ++nne)
			      {
				int n_n_el = out.neighborlist[(n_el-1)*3+nne];

				if( n_n_el>0 && n_n_el != n_el && (el+1) != n_n_el){

				  if(PreservedElements[n_n_el-1]){
				    Geometry<Node<3> > neighb_neighb_vertices;
				    for(int pn=0; pn<NumberOfVertices; ++pn)
				      {
					neighb_neighb_vertices.push_back(*(nodes_begin + out.trianglelist[(n_n_el-1)*3+pn]-1).base());
				      }

				    // std::cout<<" Element k "<<ie->Id()<<" nodes ("<<neighb_neighb_vertices[0].Id()<<", "<<neighb_neighb_vertices[1].Id()<<", "<<neighb_neighb_vertices[2].Id()<<") "<<std::endl;

				    int shared_1 =0;
				    int shared_2 =0;

				    for(int i=0; i<NumberOfVertices; ++i)
				      {
					for(int j=0; j<NumberOfVertices; ++j)
					  {
					    if( vertices[i].Id() == neighb_neighb_vertices[j].Id() )
					      shared_1 +=1;

					    if( neighb_vertices[i].Id() == neighb_neighb_vertices[j].Id() )
					      shared_2 +=1;
					  }
				      }

				    if(shared_1==2){
				      shared_side_1=counter;
				      sn_el = n_n_el;
				    }

				    if(shared_2==2)
				      shared_side_2=counter;

				  }
				}

				counter++;

			      }


			    if(shared_side_1 ==10 || shared_side_2 ==10 || shared_side_1 == shared_side_2 ){
			      std::cout<<" Shared Elements Not found in Channge_Tip_Elements "<<std::endl;
			    }
			    else{

			      // std::cout<<" Element_1 "<<ie->Id()<<" Element_2 "<<ne->Id()<<std::endl;
			      // std::cout<<" shared_side_1 "<<shared_side_1<<" shared_side_2 "<<shared_side_2<<std::endl;
			      // std::cout<<" Neighbours Selected ("<<neighb_neighb_elems[0].Id()<<", "<<neighb_neighb_elems[1].Id()<<", "<<neighb_neighb_elems[2].Id()<<") "<<std::endl;


			      //the order is important because the iterator (ne) changes after the assignation

			      // std::cout<<" Previous Neighbours 1: "<<el+1<<" ("<<out.neighborlist[el*3]<<", "<<out.neighborlist[el*3+1]<<", "<<out.neighborlist[el*3+2]<<")"<<std::endl;

			      // std::cout<<" Previous Neighbours 2: "<<n_el<<" ("<<out.neighborlist[(n_el-1)*3]<<", "<<out.neighborlist[(n_el-1)*3+1]<<", "<<out.neighborlist[(n_el-1)*3+2]<<")"<<std::endl;


			      int neighb_1 = out.neighborlist[(n_el-1)*3+shared_side_1];
			      int neighb_2 = out.neighborlist[(n_el-1)*3+shared_side_2];

			      if(option==2){
				//set neighbours of el
				out.neighborlist[el*3]   = -1;
				out.neighborlist[el*3+1] = n_el;
				out.neighborlist[el*3+2] = neighb_1;

				//set neighbours of n_el
				out.neighborlist[(n_el-1)*3]   = el+1;
				out.neighborlist[(n_el-1)*3+1] = -1;
				out.neighborlist[(n_el-1)*3+2] = neighb_2;
			      }
			      else{

				//set neighbours of el
				out.neighborlist[el*3]   = n_el;
				out.neighborlist[el*3+1] = -1;
				out.neighborlist[el*3+2] = neighb_1;

				//set neighbours of n_el
				out.neighborlist[(n_el-1)*3]   = -1;
				out.neighborlist[(n_el-1)*3+1] = el+1;
				out.neighborlist[(n_el-1)*3+2] = neighb_2;
			      }

			      // std::cout<<" Post Neighbours 1: ("<<out.neighborlist[el*3]<<", "<<out.neighborlist[el*3+1]<<", "<<out.neighborlist[el*3+2]<<")"<<std::endl;

			      // std::cout<<" Post Neighbours 2: ("<<out.neighborlist[(n_el-1)*3]<<", "<<out.neighborlist[(n_el-1)*3+1]<<", "<<out.neighborlist[(n_el-1)*3+2]<<")"<<std::endl;


			      for(int nne=0; nne<NumberOfNeighbourElements; ++nne)
				{
				  // std::cout<<" Pre Neigh Neighbours 2:" <<sn_el<<" ("<<out.neighborlist[(sn_el-1)*3]<<", "<<out.neighborlist[(sn_el-1)*3+1]<<", "<<out.neighborlist[(sn_el-1)*3+2]<<")"<<std::endl;
				  int i_el = out.neighborlist[(sn_el-1)*3+nne];
				  if( i_el>0 && i_el == n_el ){
				    out.neighborlist[(sn_el-1)*3+nne] = el+1;
				    // std::cout<<" Post Neigh Neighbours 2:" <<sn_el<<" ("<<out.neighborlist[(sn_el-1)*3]<<", "<<out.neighborlist[(sn_el-1)*3+1]<<", "<<out.neighborlist[(sn_el-1)*3+2]<<")"<<std::endl;
				    break;
				  }

				}

			      swap_performed = true;
			      any_swap_performed = true;

			      std::cout<<" swap performed in element: "<<el+1<<std::endl;
			      std::cout<<"(-";
			      for(int pn=0; pn<NumberOfVertices; ++pn)
			       	{
			       	  std::cout<<(nodes_begin + out.trianglelist[el*3+pn]-1)->Id();
				  std::cout<<"["<<(nodes_begin + out.trianglelist[el*3+pn]-1)->X()<<" "<<(nodes_begin + out.trianglelist[el*3+pn]-1)->Y()<<"] - ";
				}
			      std::cout<<")"<<std::endl;

			    }

			  }

			}

		      }
		  }


	      }
	      else{

		std::cout<<" swap NOT performed in element: "<<el<<std::endl;

	      }



	    }

	    //Nodal Neighbours will be set later in the general search:
	    //The Diagonal swapping must be done before doing it

	  }
        }

	return any_swap_performed;

	KRATOS_CATCH( "" )

    };

    /*@} */
    /**@name Operations */
    /*@{ */


    /*@} */
    /**@name Access */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

protected:
    /**@name Protected static Member Variables */
    /*@{ */


    /*@} */
    /**@name Protected member Variables */
    /*@{ */


    /*@} */
    /**@name Protected Operators*/
    /*@{ */


    /*@} */
    /**@name Protected Operations*/
    /*@{ */


    /*@} */
    /**@name Protected  Access */
    /*@{ */


    /*@} */
    /**@name Protected Inquiry */
    /*@{ */


    /*@} */
    /**@name Protected LifeCycle */
    /*@{ */



    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */

    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Access */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */

    inline double CalculateSignedArea(const double x0, const double y0,
				      const double x1, const double y1,
				      const double x2, const double y2)
    {
      return 0.5 * ( (x0*y1) - (x0*y2) - (x1*y0) + (x1*y2) + (x2*y0) - (x2*y1) );
    }


    inline void CalculateCenter(const double x0, const double y0,
				const double x1, const double y1,
				const double x2, const double y2,
				double& xc, double& yc)

    {
      xc = 0.3333333333333333333*(x0+x1+x2);
      yc = 0.3333333333333333333*(y0+y1+y2);
    }


    inline void CalculateMinDistance(const double xc, const double yc,
				     const double xc1, const double yc1,
				     const double xc2, const double yc2,
				     int &option)

    {
      double R1 = (xc-xc1)*(xc-xc1) + (yc-yc1)*(yc-yc1);
      double R2 = (xc-xc2)*(xc-xc2) + (yc-yc2)*(yc-yc2);

      if(R2 >= R1)
	option = 1;
      else
	option = 2;

    }


    /*@} */

}; /* Class ChangeTipElementsUtilities */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_CHANGE_TIP_ELEMENTS_UTILITIES_H_INCLUDED defined */
