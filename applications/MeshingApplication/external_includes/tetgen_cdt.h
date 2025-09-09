//
//   Project Name:        Kratos
//   Last Modified by:    $Author: anonymous $
//   Date:                $Date: 2009-01-15 14:50:34 $
//   Revision:            $Revision: 1.8 $
//




#if !defined(KRATOS_TETGEN_CDT_H_INCLUDED )
#define  KRATOS_TETGEN_CDT_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <cstdlib>
#include <boost/timer.hpp>



#include "tetgen.h" // Defined tetgenio, tetrahedralize().

// Project includes
#include "includes/define.h"
#include "utilities/geometry_utilities.h"
#include "containers/model.h"
#include "includes/model_part.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "meshing_application_variables.h"
#include "processes/entity_erase_process.h"

#include "utilities/binbased_fast_point_locator.h"


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
    class TetGenCDT
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of TetGenCDT
        KRATOS_CLASS_POINTER_DEFINITION(TetGenCDT);


        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        TetGenCDT()
        {}

        /// Destructor.
        virtual ~TetGenCDT() = default;


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{


        //*******************************************************************************************
        //*******************************************************************************************
        ///This function performs a Constrained Delaunay Triangulation once given "the skin" of a solid
        ///and eventually the list of internal nodes to be included in the triangulation.
        ///if required this will read the variable NODAL_H and use it to prescribe the output element size.
        ///
        ///@param ThisModelPart --> this variable contains the "skin" and the previous model part to be used for interpolation
        ///@param rReferenceElement --> this is a prototype element to be used in the generation of the new elements
        ///@param apply_volume_constraints --> if this parameter is set to true NODAL_H will be read and used
        ///@param property_id --> this is the id of the property to be assigned to the newly created elements
        ///
        ///variables needed
        ///NODAL_H it is used to assign the size we wish to have around a given node"
        ///IS_BOUNDARY --> this is an output flag which is assigned to all of the nodes on the boundary. Note that boundary nodes are not allowed to be removed
        ///
        ///the faces are stored in the "Conditions" of the model which will not be modified
        ///the nodes marked for erase will be removed here
        ///
        ///WARNING: nodes will be renumbered within this function
        tetgenio out , outnew;
        void addNodes(ModelPart& ThisModelPart , Element const& rReferenceElement,
            bool apply_volume_constraints,         unsigned int property_id , BinBasedFastPointLocator<3> element_finder)
        {
                //creating a new mesh
                boost::timer mesh_recreation_time;

                out.tetrahedronvolumelist = new double[out.numberoftetrahedra];

                int el_number = out.numberoftetrahedra;

                int counter = 0;
                for (int el = 0; el< el_number; el++)
                {
                    int old_base = el*4;
                    //calculate the prescribed h

                    ModelPart::NodesContainerType& rNodes  = ThisModelPart.Nodes();
                    ModelPart::NodesContainerType::iterator it1 = (rNodes).find( out.tetrahedronlist[old_base]);
                    ModelPart::NodesContainerType::iterator it2 = (rNodes).find( out.tetrahedronlist[old_base+1]);
                    ModelPart::NodesContainerType::iterator it3 = (rNodes).find( out.tetrahedronlist[old_base+2]);
                    ModelPart::NodesContainerType::iterator it4 = (rNodes).find( out.tetrahedronlist[old_base+3]);

                    if ( it1 == ThisModelPart.Nodes().end() )
                        KRATOS_THROW_ERROR(std::logic_error,"trying to use an inexisting node with id ",it1->Id());
                    if ( it2 == ThisModelPart.Nodes().end() )
                        KRATOS_THROW_ERROR(std::logic_error,"trying to use an inexisting node with id ",it2->Id());
                    if ( it3 == ThisModelPart.Nodes().end() )
                        KRATOS_THROW_ERROR(std::logic_error,"trying to use an inexisting node with id ",it3->Id());
                    if ( it4 == ThisModelPart.Nodes().end() )
                        KRATOS_THROW_ERROR(std::logic_error,"trying to use an inexisting node with id ",it4->Id());

                    Node::Pointer pn1 =  *it1.base();
                    Node::Pointer pn2 =  *it2.base();
                    Node::Pointer pn3 =  *it3.base();
                    Node::Pointer pn4 =  *it4.base();

                    double prescribed_h = (pn1)->FastGetSolutionStepValue(NODAL_H);
                    prescribed_h += (pn2)->FastGetSolutionStepValue(NODAL_H);
                    prescribed_h += (pn3)->FastGetSolutionStepValue(NODAL_H);
                    prescribed_h += (pn4)->FastGetSolutionStepValue(NODAL_H);
                    prescribed_h *= 0.25;


                    //if h is the height of a perfect tetrahedra, the edge size is edge = sqrt(3/2) h
                    //filling in the list of "IDEAL" tetrahedron volumes=1/12 * (edge)^3 * sqrt(2)~0.11785* h^3=
                    //0.2165063509*h^3
                    out.tetrahedronvolumelist[counter] = 0.217*prescribed_h*prescribed_h*prescribed_h;

                    counter += 1;
                }

                //char mesh_regen_opts[] = "u0MrqYYnSJCC";
                char mesh_regen_opts[] = "u0MrqYYnSJQ";
                tetrahedralize(mesh_regen_opts, &out, &outnew);
                KRATOS_WATCH("Adaptive remeshing executed")

                //cleaning unnecessary data
                out.deinitialize();
                out.initialize();

                std::cout << "mesh recreation time" << mesh_recreation_time.elapsed() << std::endl;

                /******************************************************************************/
                /******************************************************************************/
                /******************************************************************************/
                //here we create and interpolate the new kratos nodes
                ModelPart::NodesContainerType list_of_new_nodes;

                //Node::DofsContainerType& reference_dofs = (ThisModelPart.NodesBegin())->GetDofs();

                int n_points_before_refinement = ThisModelPart.Nodes().size();

                KRATOS_WATCH(n_points_before_refinement)
                    KRATOS_WATCH(outnew.numberofpoints)

                    //
                    InterpolateAndAddNewNodes(ThisModelPart, outnew,element_finder);

                //set the coordinates to the original value
                if (ThisModelPart.NodesBegin()->SolutionStepsDataHas(DISPLACEMENT)==true )
                {
                    for ( ModelPart::NodesContainerType::iterator it =  list_of_new_nodes.begin(); it!=list_of_new_nodes.end(); it++)
                    {
                        const array_1d<double,3>& disp = (it)->FastGetSolutionStepValue(DISPLACEMENT);
                        (it)->X0() = (it)->X() - disp[0];
                        (it)->Y0() = (it)->Y() - disp[1];
                        (it)->Z0() = (it)->Z() - disp[2];
                    }
                }

                //now add the elements
                AddElementsToModelPart(ThisModelPart,outnew,property_id,rReferenceElement);
            return  ;
        }


        void GenerateCDT(
            ModelPart& ThisModelPart ,
            Element const& rReferenceElement,
            bool apply_volume_constraints,
            unsigned int property_id)
        {

            KRATOS_TRY

             if (apply_volume_constraints==true && ThisModelPart.NodesBegin()->SolutionStepsDataHas(NODAL_H)==false )
                    KRATOS_THROW_ERROR(std::logic_error,"Add  ----NODAL_H---- variable!!!!!! ERROR","");
            if (ThisModelPart.NodesBegin()->SolutionStepsDataHas(IS_BOUNDARY)==false )
                KRATOS_THROW_ERROR(std::logic_error,"Add  ----IS_BOUNDARY---- variable!!!!!! ERROR","");

            //mark as IS_BOUNDARY the nodes on the "skin". This nodes will be maintained even if they were marked for erase by the user
            for (ModelPart::NodesContainerType::iterator inode = ThisModelPart.NodesBegin(); inode!=ThisModelPart.NodesEnd(); inode++)
                inode->FastGetSolutionStepValue(IS_BOUNDARY) = 1.0;
            for (ModelPart::ConditionsContainerType::iterator icond = ThisModelPart.ConditionsBegin(); icond!=ThisModelPart.ConditionsEnd(); icond++)
            {
                Geometry<Node >& geom = icond->GetGeometry();
                for (unsigned int i=0; i<geom.size(); i++)
                    geom[i].FastGetSolutionStepValue(IS_BOUNDARY) = 1.0;
            }

            //generate a temporary model part with the original list of nodes and the original list of elements
            ModelPart& AuxModelPart = ThisModelPart.GetModel().CreateModelPart("auxiliary");
            /*        AuxModelPart.Nodes().reserve(ThisModelPart.Nodes().size());
            for (ModelPart::NodesContainerType::iterator i_node = ThisModelPart.NodesBegin() ; i_node != ThisModelPart.NodesEnd() ; i_node++)
            (AuxModelPart.Nodes()).push_back(*(i_node.base()));*/
            AuxModelPart.Elements().reserve(ThisModelPart.Elements().size());
            for (ModelPart::ElementsContainerType::iterator i_el = ThisModelPart.ElementsBegin() ; i_el != ThisModelPart.ElementsEnd() ; i_el++)
                (AuxModelPart.Elements()).push_back(*(i_el.base()));

            //erase the elements from the model part
            ThisModelPart.Elements().clear();

            //clean up the list of nodes removing the nodes that should be erased
            //note that some of this nodes will be pointed to by some of the elements
            //and thus they will not be erased until no element points to them
            ModelPart::NodesContainerType temp_nodes_container;
            temp_nodes_container.reserve(ThisModelPart.Nodes().size());
            temp_nodes_container.swap(ThisModelPart.Nodes());
            for (ModelPart::NodesContainerType::iterator i_node = temp_nodes_container.begin() ; i_node != temp_nodes_container.end() ; i_node++)
            {
                if ( static_cast<bool>(i_node->FastGetSolutionStepValue(IS_BOUNDARY)) == true)
                    i_node->Set(TO_ERASE, false);
                if ( static_cast<bool>(i_node->Is(TO_ERASE)) == false)
                    (ThisModelPart.Nodes()).push_back(*(i_node.base()));
            }

            //reorder node Ids consecutively
            unsigned int id=1;
            for (ModelPart::NodesContainerType::iterator i_node = temp_nodes_container.begin() ; i_node != temp_nodes_container.end() ; i_node++)
                i_node->SetId(id++);

            //construct spatial structure with an auxiliary model part
            BinBasedFastPointLocator<3> element_finder(AuxModelPart);
            element_finder.UpdateSearchDatabase();

            //clearing the elements in the model part (note that some nodes will not be erased until the rInterpolationElements will be actually erased
            ThisModelPart.Elements().clear();

            /******************************************************************************/
            /******************************************************************************/
            /******************************************************************************/
            //prepare for meshing the CDT, passing the list of faces and the list of nodes
            //in this phase we simply use the nodes without adding anything
            tetgenio in;

            // All indices start from 1.
            in.firstnumber = 1;
            in.numberofpoints = ThisModelPart.Nodes().size();
            in.pointlist = new REAL[in.numberofpoints * 3];

            //writing the point coordinates in a vector
            ModelPart::NodesContainerType::iterator nodes_begin = ThisModelPart.NodesBegin();
            for (unsigned int i = 0; i<ThisModelPart.Nodes().size(); i++)
            {
                int base = i*3;
                in.pointlist[base] = (nodes_begin + i)->X();
                in.pointlist[base+1] = (nodes_begin + i)->Y();
                in.pointlist[base+2] = (nodes_begin + i)->Z();
            }


            in.numberoffacets = ThisModelPart.Conditions().size();;
            in.facetmarkerlist = new int[in.numberoffacets];
            in.facetlist = new tetgenio::facet[in.numberoffacets];
            ModelPart::ConditionsContainerType::iterator it = ThisModelPart.ConditionsBegin();
            tetgenio::facet *f;
            tetgenio::polygon *p;
            for (int ii=0; ii<in.numberoffacets ; ++ii)
            {
                f = &in.facetlist[ii];
                f->numberofpolygons = 1;
                f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
                f->numberofholes = 0;
                f->holelist = nullptr;

                Geometry<Node >& geom = (it)->GetGeometry();

                p = &f->polygonlist[0];
                p->numberofvertices = 3;
                p->vertexlist = new int[p->numberofvertices];
                p->vertexlist[0] = geom[0].Id();
                p->vertexlist[1] = geom[1].Id();
                p->vertexlist[2] = geom[2].Id();

                //increase counter
                it++;
            }

            //   char tetgen_options[] = "pqYYfnJu0MCC";
            char tetgen_options[] = "pqYYfnJu0MQ";
            tetrahedralize(tetgen_options, &in, &out);
            //freeing unnecessary memory
            in.deinitialize();
            in.initialize();

            InterpolateAndAddNewNodes(ThisModelPart, out,element_finder);

            /******************************************************************************/
            /******************************************************************************/
            /******************************************************************************/
            //in this second phase, take as input the "out" of the first phase
            //and improve it
            //HERE WE ADD THE VOLUME CONSTRAINT  -the desired volume of the equidistant tertrahedra
            //based upon average nodal_h (that is the  "a" switch
            //char regeneration_options[] = "rQJYq1.8anS";
            if (apply_volume_constraints==true)
            {
                    addNodes(ThisModelPart, rReferenceElement, true,property_id, element_finder);
            }
            else //do not add nodes
            {
                //now add the elements
                AddElementsToModelPart(ThisModelPart,out,property_id,rReferenceElement);
            }


            std::vector<unsigned int> el_per_node( outnew.numberofpoints, 0 );
            for (int iii = 0; iii< outnew.numberoftetrahedra*4; iii++)
                el_per_node[ outnew.tetrahedronlist[iii]-1 ] += 1;

            //mark for erasal nodes with no element attached
            for(unsigned int i=0; i< static_cast<unsigned int>(outnew.numberofpoints); i++)
            {
                if(el_per_node[i] == 0) //node is alone
                    ThisModelPart.Nodes().find( i+1 )->Set(TO_ERASE,true);
            }

            //do erasing
            ModelPart::NodesContainerType aux_nodes_container;
            aux_nodes_container.reserve(ThisModelPart.Nodes().size());
            aux_nodes_container.swap(ThisModelPart.Nodes());
            for (ModelPart::NodesContainerType::iterator i_node = aux_nodes_container.begin() ; i_node != aux_nodes_container.end() ; i_node++)
            {
                if ( static_cast<bool>(i_node->Is(TO_ERASE)) == false)
                    (ThisModelPart.Nodes()).push_back(*(i_node.base()));
            }

            // Remove  auxiliar model part
            ThisModelPart.GetModel().DeleteModelPart("auxiliary");

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

        /// Turn back information as a string.
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
        void AddElementsToModelPart(ModelPart& rModelPart, tetgenio& tet, unsigned int property_id,Element const& rReferenceElement)
        {
            KRATOS_TRY

                //***********************************************************************************
                //***********************************************************************************
                boost::timer adding_elems;
            //add preserved elements to the kratos
            Properties::Pointer properties = rModelPart.GetMesh().pGetProperties(property_id);
            //ModelPart::NodesContainerType::iterator nodes_begin = rModelPart.NodesBegin();
            (rModelPart.Elements()).reserve(tet.numberoftetrahedra);

            for (int iii = 0; iii< tet.numberoftetrahedra; iii++)
            {
                int id = iii + 1;
                int base = iii * 4;

                ModelPart::NodesContainerType& ModelNodes = rModelPart.Nodes();
                ModelPart::NodesContainerType::iterator it1 = (ModelNodes).find( tet.tetrahedronlist[base]);
                ModelPart::NodesContainerType::iterator it2 = (ModelNodes).find( tet.tetrahedronlist[base+1]);
                ModelPart::NodesContainerType::iterator it3 = (ModelNodes).find( tet.tetrahedronlist[base+2]);
                ModelPart::NodesContainerType::iterator it4 = (ModelNodes).find( tet.tetrahedronlist[base+3]);

                if ( it1 == rModelPart.Nodes().end() )
                    KRATOS_THROW_ERROR(std::logic_error,"trying to use an inexisting node with id ",it1->Id());
                if ( it2 == rModelPart.Nodes().end() )
                    KRATOS_THROW_ERROR(std::logic_error,"trying to use an inexisting node with id ",it2->Id());
                if ( it3 == rModelPart.Nodes().end() )
                    KRATOS_THROW_ERROR(std::logic_error,"trying to use an inexisting node with id ",it3->Id());
                if ( it4 == rModelPart.Nodes().end() )
                    KRATOS_THROW_ERROR(std::logic_error,"trying to use an inexisting node with id ",it4->Id());

                Node::Pointer pn1 =  *it1.base();
                Node::Pointer pn2 =  *it2.base();
                Node::Pointer pn3 =  *it3.base();
                Node::Pointer pn4 =  *it4.base();

                Tetrahedra3D4<Node > geom( pn1,pn2,pn3,pn4  );

                Element::Pointer p_element = rReferenceElement.Create(id, geom, properties);
                (rModelPart.Elements()).push_back(p_element);

            }
            std::cout << "time for adding elems" << adding_elems.elapsed() << std::endl;;
            rModelPart.Elements().Sort();

            KRATOS_CATCH("")
        }

        void InterpolateAndAddNewNodes(ModelPart& rModelPart, tetgenio& tet, BinBasedFastPointLocator<3>& element_finder)
        {
            unsigned int n_points_before_refinement = rModelPart.Nodes().size();

            //if the refinement was performed, we need to add it to the model part.
            if (static_cast<unsigned int>(tet.numberofpoints)>n_points_before_refinement)
            {
                //definitions for spatial search
//                 typedef Node PointType;
//                 typedef Node ::Pointer PointTypePointer;
                array_1d<double, 4 > N;
                const int max_results = 10000;
                BinBasedFastPointLocator<3>::ResultContainerType results(max_results);

                Node::DofsContainerType& reference_dofs = (rModelPart.NodesBegin())->GetDofs();

                int step_data_size = rModelPart.GetNodalSolutionStepDataSize();

                //TODO: parallelize this loop
                for (int i = n_points_before_refinement; i<tet.numberofpoints; i++)
                {
                    int id=i+1;
                    int base = i*3;
                    double& x= tet.pointlist[base];
                    double& y= tet.pointlist[base+1];
                    double& z= tet.pointlist[base+2];

                    Node::Pointer pnode = rModelPart.CreateNewNode(id,x,y,z);

                    //putting the new node also in an auxiliary list
                    //KRATOS_WATCH("adding nodes to list")
                    //list_of_new_nodes.push_back( pnode );

                    //std::cout << "new node id = " << pnode->Id() << std::endl;
                    //generating the dofs
                    for (Node::DofsContainerType::iterator iii = reference_dofs.begin();    iii != reference_dofs.end(); iii++)
                    {
                        Node::DofType &rDof = **iii;
                        Node::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );

                        (p_new_dof)->FreeDof();
                    }

                    //do interpolation
                    auto result_begin = results.begin();
                    Element::Pointer pelement;

                    bool is_found = element_finder.FindPointOnMesh(pnode->Coordinates(), N, pelement, result_begin, max_results);


                    if (is_found == true)
                    {
                        Geometry<Node >& geom = pelement->GetGeometry();

                        Interpolate( geom, N, step_data_size, pnode);
                    }


                }
            }
            std::cout << "During refinement we added " << tet.numberofpoints-n_points_before_refinement<< "nodes " <<std::endl;
        }


        //////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////
        void Interpolate( Geometry<Node >& geom, const array_1d<double,4>& N,
            unsigned int step_data_size,
            Node::Pointer pnode)
        {
            unsigned int buffer_size = pnode->GetBufferSize();


            for (unsigned int step = 0; step<buffer_size; step++)
            {

                //getting the data of the solution step
                double* step_data = (pnode)->SolutionStepData().Data(step);


                double* node0_data = geom[0].SolutionStepData().Data(step);
                double* node1_data = geom[1].SolutionStepData().Data(step);
                double* node2_data = geom[2].SolutionStepData().Data(step);
                double* node3_data = geom[3].SolutionStepData().Data(step);

                //copying this data in the position of the vector we are interested in
                for (unsigned int j= 0; j<step_data_size; j++)
                {

                    step_data[j] = N[0]*node0_data[j] + N[1]*node1_data[j] + N[2]*node2_data[j] + N[3]*node3_data[j];


                }
            }
            pnode->FastGetSolutionStepValue(IS_BOUNDARY)=0.0;
        }



        template< class T, std::size_t dim >
        class PointDistance
        {
        public:
            double operator()( T const& p1, T const& p2 )
            {
                double dist = 0.0;
                for ( std::size_t i = 0 ; i < dim ; i++)
                {
                    double tmp = p1[i] - p2[i];
                    dist += tmp*tmp;
                }
                return sqrt(dist);
            }
        };

        template< class T, std::size_t dim >
        class DistanceCalculator
        {
        public:
            double operator()( T const& p1, T const& p2 )
            {
                double dist = 0.0;
                for ( std::size_t i = 0 ; i < dim ; i++)
                {
                    double tmp = p1[i] - p2[i];
                    dist += tmp*tmp;
                }
                return dist;
            }
        };

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
        TetGenCDT& operator=(TetGenCDT const& rOther);


        ///@}

    }; // Class TetGenCDT

    ///@}

    ///@name Type Definitions
    ///@{


    ///@}
    ///@name Input and output
    ///@{


    /// input stream function
    inline std::istream& operator >> (std::istream& rIStream,
        TetGenCDT& rThis);

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream,
        const TetGenCDT& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
    ///@}


}  // namespace Kratos.

#endif // KRATOS_TETGEN_CDT_H_INCLUDED  defined
