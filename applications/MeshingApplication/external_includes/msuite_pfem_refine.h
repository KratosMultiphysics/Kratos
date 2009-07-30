//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2009-01-22 17:13:57 $
//   Revision:            $Revision: 1.5 $
//
//


#if !defined(KRATOS_MSUITE_PFEM_MODELER_H_INCLUDED )
#define  KRATOS_MSUITE_PFEM_MODELER_H_INCLUDED



// System includes
#include <string>
#include <iostream> 
#include <stdlib.h>


#include <boost/timer.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/line_2d_2.h"
#include "meshing_application.h"
#include "processes/node_erase_process.h" 
#include "spatial_containers/spatial_containers.h"

//includes of the msuite
#include "utiles.h"
#include "malla.h"



namespace Kratos {


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
    class MSuitePFEMModeler {
    public:
        ///@name Type Definitions
        ///@{
            typedef Node < 3 > PointType;
            typedef Node < 3 > ::Pointer PointPointerType;
            typedef std::vector<PointType::Pointer> PointVector;
            typedef PointVector::iterator PointIterator;
            typedef std::vector<double> DistanceVector;
            typedef std::vector<double>::iterator DistanceIterator;
            typedef Bucket < 3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > BucketType;
            typedef Tree< KDTreePartition<BucketType> > kd_tree;


        /// Pointer definition of TriGenModeler
        KRATOS_CLASS_POINTER_DEFINITION(MSuitePFEMModeler);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.

        MSuitePFEMModeler() {
        }
        //mpNodeEraseProcess(NULL){} //dimension = number of nodes

        /// Destructor.

        virtual ~MSuitePFEMModeler() {
        }


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{


        //*******************************************************************************************
        //*******************************************************************************************

        void ReGenerateMesh(
                ModelPart& ThisModelPart,
                Element const& rReferenceElement,
                Condition const& rReferenceBoundaryCondition,
                NodeEraseProcess& node_erase, bool rem_nodes = true, bool add_nodes = true,
                double my_alpha = 1.4, double h_factor = 0.5) {

            KRATOS_TRY
            if (ThisModelPart.NodesBegin()->SolutionStepsDataHas(IS_FREE_SURFACE) == false)
                KRATOS_ERROR(std::logic_error, "Add  ----IS_FREE_SURFACE---- variable!!!!!! ERROR", "");
            if (ThisModelPart.NodesBegin()->SolutionStepsDataHas(IS_STRUCTURE) == false)
                KRATOS_ERROR(std::logic_error, "Add  ----IS_STRUCTURE---- variable!!!!!! ERROR", "");
            if (ThisModelPart.NodesBegin()->SolutionStepsDataHas(IS_BOUNDARY) == false)
                KRATOS_ERROR(std::logic_error, "Add  ----IS_BOUNDARY---- variable!!!!!! ERROR", "");
            if (ThisModelPart.NodesBegin()->SolutionStepsDataHas(IS_FLUID) == false)
                KRATOS_ERROR(std::logic_error, "Add  ----IS_FLUID---- variable!!!!!! ERROR", "");

            KRATOS_WATCH("Msuite PFEM Refining Mesher")
            boost::timer auxiliary;

            //clean up the cloud of nodes prior to meshing
            CleanCloudOfNodes(ThisModelPart,node_erase,h_factor);

            KRATOS_WATCH(ThisModelPart.Nodes().size());
            KRATOS_WATCH("just after cleaning the cloud of nodes");

            //            ThisModelPart.Elements().clear();
            //            ThisModelPart.Conditions().clear();

            //interface with mesh suite
            ::malla m;

            int i, j;

            //ensure that the numeration of the nodes in Kratos is consecutive
            unsigned int index_I = 1;
            for (ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin(); in != ThisModelPart.NodesEnd(); in++)
            {
                in->SetId(index_I++);
//                in->Id() = index_I++;
                Node<3>::DofsContainerType& node_dofs = in->GetDofs();
                for(Node<3>::DofsContainerType::iterator iii = node_dofs.begin();    iii != node_dofs.end(); iii++)
		{
                    iii->SetId(in->Id());
		}
            }

            //creating array of new h values
            array1<double> new_h( ThisModelPart.Nodes().size() );


            // nodos, bounding box y h minimo
            int nlen = int(ThisModelPart.Nodes().size());
            m.n.resize(nlen);
            nodo p(&(ThisModelPart.NodesBegin()->X()));
            m.pmin = m.pmax = p;
            m.hayh = true;
            new_h[0] = (ThisModelPart.NodesBegin())->FastGetSolutionStepValue(NODAL_H);

            if (m.hayh)
                m.hmin = p.h = ThisModelPart.NodesBegin()->FastGetSolutionStepValue(NODAL_H,1); //h at the old step (before refinement)
            m.n += p;
            for (i = 1; i < nlen; i++) {
                new_h[i] = (ThisModelPart.NodesBegin() + i)->FastGetSolutionStepValue(NODAL_H);
                p.setpos(punto(&((ThisModelPart.NodesBegin() + i)->X())));
                p.set_min_max(m.pmin, m.pmax);
                if (m.hayh) {
                    double hi = (ThisModelPart.NodesBegin() + i)->FastGetSolutionStepValue(NODAL_H,1); //h at the old step (before refinement)
                    set_min(m.hmin, hi);
                    p.h = hi;
                }
                m.n += p;
            }

            // hay z?
            bool es3D = m.pmax[2] - m.pmin[2] > ERRADM;
            if (!es3D) m.tipo.set(m_planaxy);

            // epsilon (para varios usos)
            if (m.hayh)
                m.epsilon = m.hmin / 1000;
            else
                m.epsilon_bb();

            m.o = new octree(m.n, m.pmin, m.pmax, es3D ? 3 : 2); // octree
            for (i = 0; i < nlen; i++) m.o->add(i);

            // elementos y elementos de cada nodo
            bool pasa_elementos = false;
            if (pasa_elementos) {
                int nv = int( (ThisModelPart.ElementsBegin())->GetGeometry().size());
                elemento ei((nv == 2) ? e_segmento : ((nv == 3) ? e_triangulo : e_tetraedro));
                int elen = int(ThisModelPart.Elements().size());
                m.e.resize(elen);
                for (i = 0; i < elen; i++) {
                    Geometry< Node < 3 > >& geom = (ThisModelPart.ElementsBegin() + i)->GetGeometry();
                    for (j = 0; j < nv; j++) {
                        int ix = int(geom[j].Id());
                        ei[j] = ix;
                        m.n[ix].e += i;
                    }
                    m.e += ei; // presupongo que es e.len-1
                }
                m.tipo.set((nv == 2) ? m_lin : ((nv == 3) ? m_sup : m_vol));
            } else m.tipo.set(m_nodos);

            ThisModelPart.Elements().clear();
            ThisModelPart.Conditions().clear();

             KRATOS_WATCH((ThisModelPart.NodesEnd() -1 )->Id());

            //pass delaunay with alpha shape
//            bool elimina_sliver = true;
//            m.delaunay(elimina_sliver); //use "old h nodal"
//            m.alpha_shape(my_alpha);

            pline nodes_from;
            array1<double> shape_functions;
            m.delaunay_refinado(my_alpha, new_h, nodes_from, shape_functions ); //use "old h nodal"

            int added_nodes = m.n.len - ThisModelPart.Nodes().size();

            KRATOS_WATCH(ThisModelPart.Nodes().size());

            //create the new nodes
            Node<3>::DofsContainerType& reference_dofs = (ThisModelPart.NodesBegin())->GetDofs();
            int old_size = ThisModelPart.Nodes().size();
            for(int iii=old_size; iii<m.n.len; iii++)
            {
                unsigned int id=iii+1;
                Node<3>::Pointer pnode = ThisModelPart.CreateNewNode(id,m.n[iii][0],m.n[iii][1],m.n[iii][2]);
                pnode->SetBufferSize(ThisModelPart.NodesBegin()->GetBufferSize() );
                for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin();    iii != reference_dofs.end(); iii++)
		{
			Node<3>::DofType& rDof = *iii;
			Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );
			(p_new_dof)->FreeDof();
		}
            }

            //interpolate the new nodes
            int step_data_size = ThisModelPart.GetNodalSolutionStepDataSize();
            unsigned int buffer_size = ThisModelPart.NodesBegin()->GetBufferSize();
            unsigned int dim = 2;
            for(int iii=0; iii<added_nodes; iii++)
            {
                unsigned int base = iii * (dim+1);

                Node<3>::Pointer pnode = *(ThisModelPart.NodesBegin() + old_size + iii).base();

                for(unsigned int step=0; step<buffer_size; step++)
                {
                    double* step_data = (pnode)->SolutionStepData().Data(step);
//                    KRATOS_WATCH(pnode->Id());

                    double* node0_data = (ThisModelPart.NodesBegin() + nodes_from[base + 0] )->SolutionStepData().Data(step);
                    double* node1_data = (ThisModelPart.NodesBegin() + nodes_from[base + 1] )->SolutionStepData().Data(step);
                    double* node2_data = (ThisModelPart.NodesBegin() + nodes_from[base + 2] )->SolutionStepData().Data(step);
//
//                    KRATOS_WATCH((ThisModelPart.NodesBegin() + nodes_from[base + 0] )->Id());
//                    KRATOS_WATCH((ThisModelPart.NodesBegin() + nodes_from[base + 0] )->FastGetSolutionStepValue(NODAL_H));
//                    KRATOS_WATCH((ThisModelPart.NodesBegin() + nodes_from[base + 1] )->Id());
//                    KRATOS_WATCH((ThisModelPart.NodesBegin() + nodes_from[base + 1] )->FastGetSolutionStepValue(NODAL_H));
//                    KRATOS_WATCH((ThisModelPart.NodesBegin() + nodes_from[base + 2] )->Id());
//                    KRATOS_WATCH((ThisModelPart.NodesBegin() + nodes_from[base + 2] )->FastGetSolutionStepValue(NODAL_H));

                    //copying this data in the position of the vector we are interested in
                    for (int j = 0; j < step_data_size; j++)
                        step_data[j] =    shape_functions[base + 0] * node0_data[j]
                                        + shape_functions[base + 1] * node1_data[j]
                                        + shape_functions[base + 2] * node2_data[j];

//                    KRATOS_WATCH((ThisModelPart.NodesBegin() + old_size + iii)->FastGetSolutionStepValue(NODAL_H));

                }




            }




            //m.e[1000][1] //segundio nodo del elemento 1000
            //m.e[1000].n[1] //exactamente lo mismo de antes
            //m.e[1000].n //todo el array

            //loop over all elements -> copy them to kratos data structure
            Properties::Pointer properties = ThisModelPart.GetMesh().pGetProperties(1);
            ModelPart::NodesContainerType::iterator nodes_begin = ThisModelPart.NodesBegin();
            array1< elemento >& el_list = m.e;
            (ThisModelPart.Elements()).reserve(el_list.len);
            for (int i = 0; i < el_list.len; i++) {
                int id = i + 1;
                const elemento& ei = el_list[i];

                Triangle2D3<Node < 3 > > geom(
                        *((nodes_begin + ei[0]).base()),
                        *((nodes_begin + ei[1]).base()),
                        *((nodes_begin + ei[2]).base())
                        );

                Element::Pointer p_element = rReferenceElement.Create(id, geom, properties);
                (ThisModelPart.Elements()).push_back(p_element);
            }

            //generate the "faces" of the new discretization
            m.mk_frontera();
            int counter = 1;
            for (int k = 0; k < m.frontera.len; k++) //take care! there is AN ARRAY of meshes //exterior boundary k=0
            {
                const array1< elemento >& el_list_fk = m.frontera[k].e;
                for (int i = 0; i < el_list_fk.len; i++) {
                    int id = counter++;
                    const elemento& ei = el_list_fk[i];

//                    Line2D2<Node < 3 > > geom(
//                            *((nodes_begin + ei[0]).base()),
//                            *((nodes_begin + ei[1]).base())
//                            );
                    Line2D2<Node < 3 > > geom(
                            *((nodes_begin + ei[1]).base()),
                            *((nodes_begin + ei[0]).base())
                            );
                    Condition::Pointer p_cond = rReferenceBoundaryCondition.Create(id, geom, properties);
                    (ThisModelPart.Conditions()).push_back(p_cond);

                }
            }



            //put a flag on all the nodes of boundary
            for (ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin(); in != ThisModelPart.NodesEnd(); in++)
                in->FastGetSolutionStepValue(IS_BOUNDARY) = 0;
            for (int k = 0; k < m.frontera.len; k++) //take care! there is AN ARRAY of meshes //exterior boundary k=0
            {
                const pline& n_list_fk = m.frontera[k].n;
                for (int i = 0; i < n_list_fk.len; i++)
                    (nodes_begin + n_list_fk[i])->FastGetSolutionStepValue(IS_BOUNDARY) = 1;
            }

            //give to kratos the nodal neighbours
            //            m.mk_nn(); //create list of nodal neighbours
            //            for(int i=0; i<m.nn.len; i++)
            //            {
            //                //cpline es una lista cerrada
            //                const cpline& nni= m.nn[i];
            //                for(int j=0; j<nni.len; j++)
            //                    nni[j] //es el indice
            //            }

            //give to kratos the elemental neighbours
            //            m.mk_vecino();
            //            for(int i=0; i<m.e.len; i++)
            //            {
            //                pline& vi = m.vecino[i];
            //                for(int j=0; j<2; j++)
            //                    vi[j] //if there is no neighbour it returns <0
            //            }








            //now use updated information to construct the Kratos Mesh

            //            //assign new h to the nodes
            //            for(int i =0; i<m.n.len; i++)
            //                m.n[i].h = (ThisModelPart.NodesBegin() + i)->FastGetSolutionStepValue(NODAL_H);
            //
            //
            //            malla vieja(m);
            //            bool elimina_sliver = true;
            //            m.delaunay_refinado(vieja,my_alpha,elimina_sliver);




            KRATOS_WATCH("Finished remeshing with Msuite_PFEM_Refine")

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

        virtual std::string Info() const {
            return "";
        }

        /// Print information about this object.

        virtual void PrintInfo(std::ostream& rOStream) const {
        }

        /// Print object's data.

        virtual void PrintData(std::ostream& rOStream) const {
        }


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

        void CleanCloudOfNodes(
                ModelPart& ThisModelPart,
                NodeEraseProcess& node_erase,
                double h_factor)
        {
            KRATOS_TRY
            //remove nodes that are too close to the boundary in a elementwise sense

            //put all nodes into an octtree
            unsigned int bucket_size = 20;
            unsigned int max_results = 100;
            PointVector res(max_results);
            DistanceVector res_distances(max_results);

            PointVector list_of_nodes;
            list_of_nodes.reserve(ThisModelPart.Nodes().size());
            for (ModelPart::NodesContainerType::iterator i_node = ThisModelPart.NodesBegin(); i_node != ThisModelPart.NodesEnd(); i_node++) {
                (list_of_nodes).push_back(*(i_node.base()));
            }

            kd_tree nodes_tree1(list_of_nodes.begin(), list_of_nodes.end(), bucket_size);

            unsigned int n_points_in_radius;
            double radius; //radius means the distance, closer than which no node shall be allowd. if closer -> mark for erasing
            Node<3> work_point(0,0.0,0.0,0.0);
            for (ModelPart::NodesContainerType::iterator in = ThisModelPart.NodesBegin();
                    in != ThisModelPart.NodesEnd(); in++)
            {
                radius = h_factor * in->FastGetSolutionStepValue(NODAL_H,1);

                work_point[0] = in->X();
                work_point[1] = in->Y();
                work_point[2] = in->Z();

                n_points_in_radius = nodes_tree1.SearchInRadius(work_point, radius, res.begin(), res_distances.begin(), max_results);
                if (n_points_in_radius > 1) {
                    if (in->FastGetSolutionStepValue(IS_BOUNDARY,1) == 0.0 && in->FastGetSolutionStepValue(IS_STRUCTURE) == 0.0) {
                        //look if we are already erasing any of the other nodes
                        unsigned int erased_nodes = 0;
                        for (PointIterator i = res.begin(); i != res.begin() + n_points_in_radius; i++)
                            erased_nodes += (*i)->GetValue(ERASE_FLAG);

                        if (erased_nodes < 1) //we cancel the node if no other nodes are being erased
                            in->GetValue(ERASE_FLAG) = 1;

                    } else if ((in)->FastGetSolutionStepValue(IS_STRUCTURE) != 1.0) //boundary nodes will be removed if they get REALLY close to another boundary node (0.2 * h_factor)
                    {
                        //here we loop over the neighbouring nodes and if there are nodes
                        //with IS_BOUNDARY=1 which are closer than 0.2*nodal_h from our we remove the node we are considering
                        unsigned int k = 0;
                        unsigned int counter = 0;
                        for (PointIterator i = res.begin(); i != res.begin() + n_points_in_radius; i++) {
                            if ((*i)->FastGetSolutionStepValue(IS_BOUNDARY,1) == 1.0 && res_distances[k] < 0.2 * radius && res_distances[k] > 0.0) {
                                counter += 1;
                            }
                            k++;
                        }
                        if (counter > 0)
                            in->GetValue(ERASE_FLAG) = 1;
                    }
                }

                //perform the clean
            }

            //now loop on all elements, identify the nodes which are too close to the boundary and mark them for removal
            boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
            array_1d<double,3> N;
//            for (ModelPart::ElementsContainerType::iterator ie = ThisModelPart.ElementsBegin();
//                    ie != ThisModelPart.ElementsEnd(); ie++)
//            {
//                Geometry<Node<3> >& geom = ie->GetGeometry();
//                for(unsigned int i=0; i<geom.size(); i++)
//                {
//                   if(geom[i].FastGetSolutionStepValue(IS_STRUCTURE) == 0) //we identify the node that is "free" to move
//                         { geom[i].GetValue(ERASE_FLAG) = 1;}
//                }
//            }

            for (ModelPart::ElementsContainerType::iterator ie = ThisModelPart.ElementsBegin();
                    ie != ThisModelPart.ElementsEnd(); ie++)
            {
                Geometry<Node<3> >& geom = ie->GetGeometry();

                unsigned int nstructure = 0;
                for(unsigned int i=0; i<geom.size(); i++)
                    nstructure += int(geom[i].FastGetSolutionStepValue(IS_STRUCTURE) );

                if(nstructure == 2)
                {
                    //calculate shape function derivatives
                    double Volume;
                    GeometryUtils::CalculateGeometryData(geom,DN_DX,N,Volume);

                    //calculate avg h
                    double havg = geom[0].FastGetSolutionStepValue(NODAL_H) +
                                  geom[1].FastGetSolutionStepValue(NODAL_H) +
                                  geom[2].FastGetSolutionStepValue(NODAL_H);
                    havg *= 0.333333333333333;

                     for(unsigned int iii=0; iii<geom.size(); iii++)
                     {
                         if(geom[iii].FastGetSolutionStepValue(IS_STRUCTURE) != 1.0) //we identify the node that is "free" to move
                         {
                             double DNnorm = sqrt( DN_DX(iii,0)*DN_DX(iii,0) + DN_DX(iii,1)*DN_DX(iii,1) );
                             double h = 1.0 / DNnorm;

                             if(h < havg * 0.4) //cancel it if it gets too close
                             {
                                 KRATOS_WATCH(geom[iii].Id());
                                 geom[iii].GetValue(ERASE_FLAG) = 1;
                             }
                         }

                     }
                }
            }

//            ThisModelPart.Elements().clear();
//            ThisModelPart.Conditions().clear();


//             Node<3>::Pointer temp = ThisModelPart.Nodes()(807);

//             KRATOS_WATCH(temp->GetValue(ERASE_FLAG));
//             KRATOS_WATCH(temp->Id());

//             KRATOS_WATCH(ThisModelPart.Nodes().size());

            //perform the removal
            node_erase.Execute();
//             KRATOS_WATCH(temp->Id());
//
//             KRATOS_WATCH(ThisModelPart.Nodes().size());
//            KRATOS_WATCH(*temp)

            KRATOS_CATCH("");
        }

//        void Interpolate(ModelPart::NodesContainerType& old_list_of_elements, PointerVector< Node<3> >& list_of_new_nodes)
//        {
//            KRATOS_TRY
//
//            //put the new nodes in an octtree
//            PointVector aux_list;
//            aux_list.reserve(list_of_new_nodes.size());
//            for (ModelPart::NodesContainerType::iterator i_node = list_of_new_nodes.begin(); i_node != list_of_new_nodes.end(); i_node++) {
//                (aux_list).push_back(*(i_node.base()));
//            }
//
//            //loop over the existing elements to find where they fall --> if not found throw an error
//
//            KRATOS_CATCH("");
//        }



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
            MSuitePFEMModeler & operator=(MSuitePFEMModeler const& rOther);


            ///@}

        }; // Class MSuitePFEMModeler

        ///@}

        ///@name Type Definitions
        ///@{


        ///@}
        ///@name Input and output
        ///@{


        /// input stream function
        inline std::istream & operator >>(std::istream& rIStream,
                MSuitePFEMModeler& rThis);

        /// output stream function

        inline std::ostream & operator <<(std::ostream& rOStream,
                const MSuitePFEMModeler& rThis) {
            rThis.PrintInfo(rOStream);
            rOStream << std::endl;
            rThis.PrintData(rOStream);

            return rOStream;
        }
        ///@}


    } // namespace Kratos.

#endif // KRATOS_MSUITE_EXTERNAL_H_INCLUDED  defined



