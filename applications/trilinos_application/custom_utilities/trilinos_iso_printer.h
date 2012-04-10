
//   Project Name:        Kratos       
//   Last Modified by:    $Author: pbecker $
//   Date:                $Date: 2010-05-12 $
//   Revision:            $Revision: 1.0 $
//
//

#if !defined(KRATOS_TRILINOS_LOCAL_ISO_PRINTING_APP)
#define  KRATOS_TRILINOS_LOCAL_ISO_PRINTING_APP


#ifdef _OPENMP
#include <omp.h>
#endif




// System includes
#include <string>
#include <iostream> 
#include <stdlib.h>
#include <cmath>
#include <algorithm>


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/dof.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"
#include "geometries/geometry.h"
#include "geometries/geometry_data.h"


// #include "containers/array_1d.h"
// #include "processes/find_nodal_neighbours_process.h"
// #include "processes/find_elements_neighbours_process.h"
#include "containers/data_value_container.h"
#include "includes/mesh.h"
#include "utilities/math_utils.h"
//#include "utilities/split_triangle.h"
#include "utilities/split_triangle.c"
#include "utilities/split_tetrahedra.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/triangle_3d_3.h"
#include "utilities/geometry_utilities.h"
//#include "processes/node_erase_process.h"
#include "custom_utilities/parallel_fill_communicator.h"
// #include "spatial_containers/spatial_containers.h"

#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_MpiComm.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"

#include "Epetra_FEVector.h"

#include "Epetra_Import.h"
#include "Epetra_MpiComm.h"
#include "parallel_fill_communicator.h"



namespace Kratos {
    /// This Function is designed to make an array of nodes that define the elements which are cut by an isosurface
    ///this way we can later print only the results in these nodes.
    //HOW TO USE IT from python.
    /*
     #constructor:
     Iso_App = TrilinosIsosurfacePrinterApplication(model_part,Comm)
     
     #and then, in each time step in which results must be printed:
     Iso_App.ClearData()
     variable1 = PRESSURE
     variable2 = DISTANCE
     isovalue1 = 0.001
     isovalue2 = 5.0
     tolerance = 0.00000000001
     Iso_App.AddScalarVarIsosurface(variable1, isovalue1, tolerance)
     Iso_App.AddScalarVarIsosurface(variable2, isovalue2, tolerance)
     #conditions can be also added
       
     used_nodes_array=Iso_App.CreateNodesArray() #after all the isosurfaces needed, the array must be created
     size=used_nodes_array.Size()

     if size!=0: #results must be printed only when the array is not empty, otherwise it will segfault
       	gid_io.WriteNodalResults(PRESSURE,used_nodes_array,time,0)
       	gid_io.WriteNodalResults(VELOCITY,used_nodes_array,time,0)
     
          
      
    */ 
    class TrilinosIsosurfacePrinterApplication {
    public:

        typedef ModelPart::NodesContainerType NodesArrayType;
        typedef ModelPart::ElementsContainerType ElementsArrayType;
        typedef ModelPart::ConditionsContainerType ConditionsArrayType;
        typedef boost::numeric::ublas::vector<Matrix> Matrix_Order_Tensor;
        typedef boost::numeric::ublas::vector<Vector> Vector_Order_Tensor;
        typedef boost::numeric::ublas::vector<Vector_Order_Tensor> Node_Vector_Order_Tensor;
        typedef Node < 3 > PointType;
        typedef Node < 3 > ::Pointer PointPointerType;
        typedef std::vector<PointType::Pointer> PointVector;
        typedef PointVector::iterator PointIterator;

        /**constructor:
         *@param ModelPart& the model part to be refined
         *@param Epetra_MpiComm the Epetra Communicator to be used
         */

        TrilinosIsosurfacePrinterApplication(ModelPart& model_part, Epetra_MpiComm& Comm) :  mr_model_part(model_part), mrComm(Comm) {
            int nnodes = mr_model_part.Nodes().size(); //both local nodes and the ones belonging to other CPU
            m_used_nodes.resize(nnodes);
            for (int index = 0; index!=nnodes;++index) m_used_nodes[index]=false; //starting as zero (no needed nodes)
            //smallest_edge=1.0;
        } //

        ~TrilinosIsosurfacePrinterApplication() {
        }


        ///************************************************************************************************
        ///************************************************************************************************


        
        
        /// ADDSKINCONDITIONS: THIS FUNCTION ADDS TO THE NEW MODEL PART THE DATA OF THE CONDITIONS BELONGING TO THE OLD MODEL PART, THE NODES COORDINATES ALREADY EXIST. WE ONLY NEED TO COPY THEM INTO THE NEW MODEL PART.
		/** this function adds the skin condtion. 
		    WARNING: They have to be triangles and it CAN'T be empty, otherwise a segmentation fault will appear
		 * @param mr_model_part . original model part
		 **/
        
        void AddSkinConditions()
        {
            ModelPart& this_model_part = mr_model_part;
            //ModelPart& new_model_part = mr_new_model_part;
	    if (mrComm.MyPID() == 0) {	
            //KRATOS_WATCH("Adding Skin Conditions to the new model part, added in layer:")
            //KRATOS_WATCH(plane_number)
            }
            KRATOS_TRY
            
                    
            NodesArrayType& rNodes_old = this_model_part.Nodes();        //needed to find the position of the node in the node array
            NodesArrayType::iterator it_begin_node_old = rNodes_old.ptr_begin();        
                    
            //int nlocal_nodes = this_model_part.GetCommunicator().LocalMesh().Nodes().size(); //only local nodes
            //int nnodes = this_model_part.Nodes().size(); //both local nodes and the ones belonging to other CPU

             //the id of our processor:
            //double this_partition_index = double(mrComm.MyPID());
            
            //here we'll store the nodes that will be cloned. NOT ids but position in the nodes array (full, not only owned)
             //bool* used_nodes = new bool  [nnodes];
             //for (int index = 0; index!=nnodes;++index) used_nodes[index]=false; //starting as zero (no needed nodes)
            
            //now we have to search for the conditions. they can only be triangles, otherwise -> not stored
            int number_of_local_conditions=0; 
            for (ModelPart::ConditionsContainerType::iterator it = this_model_part.ConditionsBegin(); it != this_model_part.ConditionsEnd(); it++)
            {
                Geometry<Node < 3 > >& geom = it->GetGeometry();
                if (geom.size()==3){
                    ++number_of_local_conditions; //conditions are always owned
                    for (unsigned int i = 0; i < geom.size(); i++) {
                        int node_position = this_model_part.Nodes().find(geom[i].Id()) - it_begin_node_old; //probably there-s a better way to do this, i only need the position in the array, (not the ID)
                        m_used_nodes[node_position]=true; //we will have to clone this node into the new model part, no matter if owned or not.
                        }
                    }   
            }

            if (mrComm.MyPID() == 0) cout << "added conditions nodes to the list of nodes to be printed "<< endl;
            KRATOS_CATCH("")
        }
        
        
        
                
                
        //template <class TDataType>    
        //this is the function to create isosurfaces from a scalar variable. the other (a component from a vectorial variable).
        void AddScalarVarIsosurface( Variable<double>& variable, double isovalue, float tolerance) {
            KRATOS_TRY
            if (mrComm.MyPID() == 0) {
                cout <<"Printing Isosurface Nodes with the following data:"<<endl;
                KRATOS_WATCH(variable);
                KRATOS_WATCH(isovalue);
            }
		
            ModelPart& this_model_part = mr_model_part;
            //ModelPart& new_model_part = mr_new_model_part;

            vector<int> Elems_In_Plane(this_model_part.Elements().size()); //our (int) vector, where we write 1 when the element is cut by the cutting plane. when it is 2, it means we have 2 triangles (4 cutting points)




                        ElementsArrayType& rElements = this_model_part.Elements();

            //double this_partition_index = double(mrComm.MyPID());

            //first of all create a matrix with no overlap
            ElementsArrayType::iterator it_begin = rElements.ptr_begin(); //+element_partition[k];
            ElementsArrayType::iterator it_end = rElements.ptr_end(); //+element_partition[k+1];

            NodesArrayType& rNodes_old = this_model_part.Nodes();        //needed to find the position of the node in the node array
            NodesArrayType::iterator it_begin_node_old = rNodes_old.ptr_begin(); 
            
            double node_value; // node to closest point in the plane
            double neigh_value; //other node of the edge (neighbour) to closest point in the plane
            double diff_node_value;                // difference between the imposed value of the variable and its value in the node
            double diff_neigh_value;; //distance between the two nodes of the edge
            double diff_node_neigh;      
            array_1d<unsigned int, 4 > list_matching_nodes; // used to save the new nodes that match exactly old nodes  (very unlikely, but might be 4 for very plane elements)    
            unsigned int exact_nodes = 0;
            unsigned int outside_nodes = 0;
            int current_element=0;
            int number_of_cuts=0;

            for (ElementsArrayType::iterator it = it_begin; it != it_end; ++it) {
                bool print_element=false;
                ++current_element;
                number_of_cuts = 0;
                exact_nodes = 0;
                outside_nodes = 0;
                Geometry<Node < 3 > >&geom = it->GetGeometry(); //geometry of the element
                for (unsigned int i = 0; i < it->GetGeometry().size(); i++) //size = 4 ; nodes per element. NOTICE WE'LL BE LOOPING THE EDGES TWICE. THIS IS A WASTE OF TIME BUT MAKES IT EASIER TO IDENTITY ELEMENTS. LOOK BELOW.
                    //when we have a triangle inside a thetraedra, its edges (or nodes) must be cut 3 times by the plane. if we loop all 2 times we can have a counter. when it's = 6 then we have a triangle. when tetraedras are cutted 8 times then we have 2 triangles (or a cuatrilateral, the same)
                {
                    node_value= geom[i].FastGetSolutionStepValue(variable);
                    diff_node_value = isovalue - node_value; // dist = (xnode-xp)*versor closest point-plane distance
                    for (unsigned int j = 0; j < it->GetGeometry().size(); j++) //  looping on the neighbours
                    {
                        //unsigned int this_node_cutted_edges=0;
                        if (i != j) //(cant link node with itself)
                        {
                            neigh_value= geom[j].FastGetSolutionStepValue(variable);
                            diff_neigh_value = isovalue - neigh_value; 
                            diff_node_neigh = node_value - neigh_value;
                            //now that we have the two points of the edge defined we can check whether it is cut by the plane or not
                            if (fabs(diff_node_value) < (tolerance)  ) //then our node is part of the plane (this should have been done before the loop on neighbours, but this way it is easier to read .
                            {
                                number_of_cuts += 2; //since its neighbour wont take this case as a cut, we must save 2 cuts instead of one. (to reach number_of_cuts=6), 
                                ++exact_nodes;
                                list_matching_nodes[i] = geom[i].Id();
                                if ((diff_node_value * diff_neigh_value) > 0.0) //if true, it means that this node, despite being close, is actually outside the element. So if we have 2 or more ouside, we will not add this triangle in this element because we consider it property of another element
                                    ++outside_nodes;
                                break; //we exit the j loop.
                            }
                            if ((diff_node_value * diff_neigh_value) < 0.0 && (fabs(diff_neigh_value)>(tolerance))) // this means one is on top of the plane and the other on the bottom, no need to do more checks, it's in between!
                            {
                                ++number_of_cuts;
                            }
                        } //closing the i!=j if
                    } //closing the neighbour (j) loop
                } //closing the nodes (i) loop


                //now we have to save the data. we should get a list with the elements that will genereate triangles and the total number of triangles
                Elems_In_Plane[current_element - 1] = 0; //we initialize as 0
                if (exact_nodes < 4 || outside_nodes<2 ) { //this means at least one new node has to be generated
                    if (number_of_cuts == 6 || number_of_cuts == 6) //it can be 8, in that case we have 2 triangles (the cut generates a square)
                        print_element=true; 

                }
                if (print_element==true)
                {
                    for (unsigned int i = 0; i < it->GetGeometry().size(); i++) //size = 4 ; nodes per element. NOTICE WE'LL BE LOOPING THE EDGES TWICE. THIS IS A WASTE OF TIME BUT MAKES IT EASIER TO IDENTITY ELEMENTS. LOOK BELOW.
                    //when we have a triangle inside a thetraedra, its edges (or nodes) must be cut 3 times by the plane. if we loop all 2 times we can have a counter. when it's = 6 then we have a triangle. when tetraedras are cutted 8 times then we have 2 triangles (or a cuatrilateral, the same)
                    {
                        int node_position = this_model_part.Nodes().find(geom[i].Id()) - it_begin_node_old;
                        m_used_nodes[node_position]=true;
                    }
                }
            } //closing the elem loop
            
            MPI_Barrier(MPI_COMM_WORLD);
            if (mrComm.MyPID() == 0) cout << "Added nodes that are part of the isosurface to the priting list" << endl;


            KRATOS_CATCH("")
        }
        
                


        /// ************************************************************************************************
        /// ************************************************************************************************

        NodesArrayType CreateNodesArray() {
            NodesArrayType IsosurfaceNodes;
            IsosurfaceNodes.clear();
            if (mrComm.MyPID() == 0) cout <<"Creating Nodes Array"<<endl;
            //int step_data_size = mr_model_part.GetNodalSolutionStepDataSize();
            int i=0;
            //looping the nodes, no data is assigned to elements
            for (ModelPart::NodesContainerType::iterator it = mr_model_part.NodesBegin(); it != mr_model_part.NodesEnd(); it++) {
              if (m_used_nodes[i]==true)
              {
                  //Node < 3 > ::Pointer pnode
                  IsosurfaceNodes.push_back(*it);

              }
              ++i;
            }//closing node loop
            //cout <<"my ID is" << mrComm.MyPID() <<"and i have"<< total_number_of_printed_nodes <<"nodes"<<endl;
            return IsosurfaceNodes;
        }//closing subroutine
        
        void ClearData()
        {
            int nnodes = mr_model_part.Nodes().size(); //both local nodes and the ones belonging to other CPU
            //m_used_nodes.resize(nnodes);
            for (int index = 0; index!=nnodes;++index) m_used_nodes[index]=false; //starting as zero (no needed nodes)
        }
        
        



    protected:
     
        //double smallest_edge;
        boost::numeric::ublas::vector<bool> m_used_nodes;
        ModelPart& mr_model_part;
        Epetra_MpiComm& mrComm;
        



    };



} // namespace Kratos.

#endif // KRATOS_TRILINOS_LOCAL_ISO_PRINTING_APP  defined 


