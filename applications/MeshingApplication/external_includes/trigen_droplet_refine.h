// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:
//

#if !defined(KRATOS_TRIGEN_DROPLET_MODELER_H_INCLUDED )
#define  KRATOS_TRIGEN_DROPLET_MODELER_H_INCLUDED

// System includes

// External includes
#if !defined(KRATOS_TRIANGLE_EXTERNAL_H_INCLUDED)
#define  KRATOS_TRIANGLE_EXTERNAL_H_INCLUDED
#include "triangle.h"
#endif

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "geometries/triangle_2d_3.h"
#include "meshing_application_variables.h"
#include "processes/node_erase_process.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/timer.h"

namespace Kratos
{
extern "C" {
    void triangulate(char *, struct triangulateio *, struct triangulateio *,struct triangulateio *);
    //void trifree();
}


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
class TriGenDropletModeler
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TriGenModeler
    KRATOS_CLASS_POINTER_DEFINITION(TriGenDropletModeler);
    typedef Node<3> PointType;
    typedef Node<3>::Pointer PointPointerType;
    typedef std::vector<PointType::Pointer>           PointVector;
    typedef PointVector::iterator PointIterator;
    typedef std::vector<double>               DistanceVector;
    typedef std::vector<double>::iterator     DistanceIterator;
    typedef Bucket<3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KdtreeType; //Kdtree;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TriGenDropletModeler() :
        mJ(ZeroMatrix(2,2)), //local jacobian
        mJinv(ZeroMatrix(2,2)), //inverse jacobian
        mC(ZeroVector(2)), //dimension = number of nodes
        mRhs(ZeroVector(2)) {}
    //mpNodeEraseProcess(NULL){} //dimension = number of nodes

    /// Destructor.
    virtual ~TriGenDropletModeler() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    //*******************************************************************************************
    //*******************************************************************************************
    void ReGenerateMeshDroplet(
        ModelPart& ThisModelPart ,
        Element const& rReferenceElement,
        Condition const& rReferenceBoundaryCondition,
        NodeEraseProcess& node_erase, bool rem_nodes = true, bool add_nodes=true,
        double my_alpha = 1.4, double h_factor=0.5)
    {

        KRATOS_TRY
        if (ThisModelPart.NodesBegin()->SolutionStepsDataHas(IS_FREE_SURFACE)==false )
            KRATOS_THROW_ERROR(std::logic_error,"Add  ----IS_FREE_SURFACE---- variable!!!!!! ERROR","");
        if (ThisModelPart.NodesBegin()->SolutionStepsDataHas(IS_STRUCTURE)==false )
            KRATOS_THROW_ERROR(std::logic_error,"Add  ----IS_STRUCTURE---- variable!!!!!! ERROR","");
        if (ThisModelPart.NodesBegin()->SolutionStepsDataHas(IS_BOUNDARY)==false )
            KRATOS_THROW_ERROR(std::logic_error,"Add  ----IS_BOUNDARY---- variable!!!!!! ERROR","");
        if (ThisModelPart.NodesBegin()->SolutionStepsDataHas(IS_FLUID)==false )
            KRATOS_THROW_ERROR(std::logic_error,"Add  ----IS_FLUID---- variable!!!!!! ERROR","");

        KRATOS_WATCH("Trigen Droplet Refining Mesher")
        const auto inital_time = std::chrono::steady_clock::now();


//clearing elements

	int step_data_size = ThisModelPart.GetNodalSolutionStepDataSize();


        ThisModelPart.Elements().clear();

	//20150909 ajarauta
	int id = (ThisModelPart.Nodes().end() - 1)->Id() + 1;
	for(ModelPart::ConditionsContainerType::iterator ic = ThisModelPart.ConditionsBegin() ;
	    ic != ThisModelPart.ConditionsEnd() ; ic++)
	    {
		if (ic->GetGeometry().size()==2)
		{
// 		    KRATOS_WATCH("LALALALLA")
		    //Original:
		    unsigned int n_flag=ic->GetGeometry()[0].FastGetSolutionStepValue(IS_FREE_SURFACE);
		    n_flag+=ic->GetGeometry()[1].FastGetSolutionStepValue(IS_FREE_SURFACE);
		    //Enhanced for any segment at the boundary
// 		    int n_flag=ic->GetGeometry()[0].FastGetSolutionStepValue(IS_BOUNDARY);
// 		    n_flag+=ic->GetGeometry()[1].FastGetSolutionStepValue(IS_BOUNDARY);
		    unsigned int n_trip=ic->GetGeometry()[0].FastGetSolutionStepValue(TRIPLE_POINT);
		    n_trip+=ic->GetGeometry()[1].FastGetSolutionStepValue(TRIPLE_POINT);
		    //Only for IS_STRUCTURE segments
		    unsigned int n_struct=ic->GetGeometry()[0].FastGetSolutionStepValue(IS_STRUCTURE);
		    n_struct+=ic->GetGeometry()[1].FastGetSolutionStepValue(IS_STRUCTURE);
		    unsigned int n_sum = n_flag + n_trip;

		    //THIS REFINES THE NODES OF INTERNAL ELEMENTS OF THE SURFACE WHERE THE INBLOW IS: FLAG_VAR=1
		    if (n_flag==ic->GetGeometry().size() || n_struct==ic->GetGeometry().size() || n_sum==ic->GetGeometry().size())
// 		    if (n_struct==ic->GetGeometry().size())
		    {
			double x0=ic->GetGeometry()[0].X();
			double y0=ic->GetGeometry()[0].Y();
			double x1=ic->GetGeometry()[1].X();
			double y1=ic->GetGeometry()[1].Y();
			double edge01=sqrt((x0-x1)*(x0-x1)+(y0-y1)*(y0-y1));

			//ajarauta: position correction of the middle node -> CHECK WHAT HAPPENS FOR TRIPLE POINT!!
			// NORMAL is set to 0 in triple_point in the monolithic_embedded_solver right now...
			/*
			double n0x, n0y, n0norm, n0xunit,n0yunit,n1x,n1y,n1norm,n1xunit,n1yunit, x01, y01, xM, yM;
			double v01norm, scprod, costh, delt_dist, ndx, ndy, ndx_d, ndy_d, ndxnorm;
			double curv0 = ic->GetGeometry()[0].FastGetSolutionStepValue(CURVATURE);
			double curv1 = ic->GetGeometry()[1].FastGetSolutionStepValue(CURVATURE);
			double mean_curv = 0.5*(curv0 + curv1);
			double is_flat = 0.0;
			double xp = 0.0;
			double yp = 0.0;
			double n_struct = 0.0;
			n_struct += ic->GetGeometry()[0].FastGetSolutionStepValue(IS_STRUCTURE);
			n_struct += ic->GetGeometry()[1].FastGetSolutionStepValue(IS_STRUCTURE);
			if (mean_curv == 0.0 || n_struct > 1.0)
			    is_flat = 1.0;
			if (mean_curv != 0.0 && n_struct > 1.0) //avoid cases when alpha shape deletes an element with two IS_STRUCTURE nodes
			{
			    n0x = ic->GetGeometry()[0].FastGetSolutionStepValue(NORMAL_X);
			    n0y = ic->GetGeometry()[0].FastGetSolutionStepValue(NORMAL_Y);
			    n0norm = sqrt(n0x*n0x + n0y*n0y);
			    n0xunit = n0x/n0norm;
			    n0yunit = n0y/n0norm;
			    n1x = ic->GetGeometry()[1].FastGetSolutionStepValue(NORMAL_X);
			    n1y = ic->GetGeometry()[1].FastGetSolutionStepValue(NORMAL_Y);
			    n1norm = sqrt(n1x*n1x + n1y*n1y);
			    n1xunit = n1x/n1norm;
			    n1yunit = n1y/n1norm;
			    x01 = x1-x0;
			    y01 = y1-y0;
			    v01norm = sqrt(x01*x01 + y01*y01);
			    scprod = (-n1yunit)*x01 + n1xunit*y01;
			    if (scprod < 0.0)
				costh = scprod/v01norm;
			    else
				costh = -scprod/v01norm;
			    delt_dist = (1.0 + costh)/mean_curv;
			    ndx = 0.5*(n0xunit+n1xunit);
			    ndy = 0.5*(n0yunit+n1yunit);
			    ndxnorm = sqrt(ndx*ndx+ndy*ndy);
			    ndx_d = ndx*delt_dist/ndxnorm;
			    ndy_d = ndy*delt_dist/ndxnorm;
			    xM = 0.5*(x0+x1);
			    yM = 0.5*(y0+y1);
			    xp = xM + ndx_d;
			    yp = yM + ndy_d;
			}
			*/

			//////////////////////////////////////////////////////////////////////////////////////////////////////////
			// IMPORTANT!!
			// The following step only makes sense if NODAL_LENGTH is the initial value and is taken as a reference
			// for remeshing. If it is calculated every time step, it will never fulfill the condition
			// edge01 > factor*nodal_h!!!
			//////////////////////////////////////////////////////////////////////////////////////////////////////////
			double nodal_h=ic->GetGeometry()[0].FastGetSolutionStepValue(NODAL_H);
			nodal_h+=ic->GetGeometry()[1].FastGetSolutionStepValue(NODAL_H);
			nodal_h*=0.5;
			//if the edge of the segment (condition) is too long, we split it into two by adding a node in the middle

			Node<3>::DofsContainerType& reference_dofs = (ThisModelPart.NodesBegin())->GetDofs();


			double factor=2.0;
			if (edge01>factor*nodal_h)
			{
			    id++;
			    double x = 0.5*(x0+x1);
			    double y = 0.5*(y0+y1);;
			    double z = 0.0;
			    Node<3>::Pointer pnode = ThisModelPart.CreateNewNode(id,x,y,z);

			    //putting the new node also in an auxiliary list
			    //KRATOS_WATCH("adding nodes to list")

			    //std::cout << "new node id = " << pnode->Id() << std::endl;
			    //generating the dofs
			    for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin(); iii != reference_dofs.end(); iii++)
			    {
				Node<3>::DofType& rDof = **iii;
				Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );
				(p_new_dof)->FreeDof();
			    }
			    Geometry<Node<3> >& geom = ic->GetGeometry();

			    InterpolateOnEdge(geom, 0, 1, step_data_size, pnode);
			    const array_1d<double,3>& disp = pnode->FastGetSolutionStepValue(DISPLACEMENT);
			    pnode->X0() = pnode->X() - disp[0];
			    pnode->Y0() = pnode->Y() - disp[1];
			    KRATOS_WATCH("Added node at the EDGE")
			}

		    }
		}
	    }


        ThisModelPart.Conditions().clear();

        ////////////////////////////////////////////////////////////



        // bucket types
        //typedef Bucket<3, PointType, ModelPart::NodesContainerType, PointPointerType, PointIterator, DistanceIterator > BucketType;
        //typedef Bins< 3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > StaticBins;
        // bucket types

        //typedef Tree< StaticBins > Bin; 			     //Binstree;
        unsigned int bucket_size = 20;

        //performing the interpolation - all of the nodes in this list will be preserved
        unsigned int max_results = 100;
        //PointerVector<PointType> res(max_results);
        //NodeIterator res(max_results);
        PointVector res(max_results);
        DistanceVector res_distances(max_results);
        Node<3> work_point(0,0.0,0.0,0.0);
        //if the remove_node switch is activated, we check if the nodes got too close

        if (rem_nodes==true)
        {
            PointVector list_of_nodes;
            list_of_nodes.reserve(ThisModelPart.Nodes().size());
            for(ModelPart::NodesContainerType::iterator i_node = ThisModelPart.NodesBegin() ; i_node != ThisModelPart.NodesEnd() ; i_node++)
            {
                (list_of_nodes).push_back(*(i_node.base()));
            }

            KdtreeType nodes_tree1(list_of_nodes.begin(),list_of_nodes.end(), bucket_size);

            RemoveCloseNodes(ThisModelPart, nodes_tree1, node_erase, h_factor);
        }

        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        //												  	//
        //		Now we shall pass the Alpha Shape for the second time, having the "bad nodes" removed	//
        //////////////////////////////////////////////////////////////////////////////////////////////////////////
        //creating the containers for the input and output
        struct triangulateio in_mid;
        struct triangulateio out_mid;
        struct triangulateio vorout_mid;

        initialize_triangulateio(in_mid);
        initialize_triangulateio(out_mid);
        initialize_triangulateio(vorout_mid);

        //assigning the size of the input containers

        in_mid.numberofpoints = ThisModelPart.Nodes().size();
        in_mid.pointlist = (REAL *) malloc(in_mid.numberofpoints * 2 * sizeof(REAL));

        //reorder the id's and give the coordinates to the mesher
        ModelPart::NodesContainerType::iterator nodes_begin = ThisModelPart.NodesBegin();
        for(unsigned int i = 0; i<ThisModelPart.Nodes().size(); i++)
        {
            int base = i*2;
            //int base = ((nodes_begin + i)->Id()   -  1 ) * 2;

            //from now on it is consecutive
            (nodes_begin + i)->SetId(i+1);
//				(nodes_begin + i)->Id() = i+1;

            in_mid.pointlist[base] = (nodes_begin + i)->X();
            in_mid.pointlist[base+1] = (nodes_begin + i)->Y();

            auto& node_dofs = (nodes_begin + i)->GetDofs();
            for(auto iii = node_dofs.begin();    iii != node_dofs.end(); iii++)
            {
                (**iii).SetEquationId(i+1);
            }
        }
        //in_mid.numberoftriangles = ThisModelPart.Elements().size();
        //in_mid.trianglelist = (int*) malloc(in_mid.numberoftriangles * 3 * sizeof(int));

        // "P" suppresses the output .poly file. Saves disk space, but you
        //lose the ability to maintain constraining segments  on later refinements of the mesh.
        // "B" Suppresses boundary markers in the output .node, .poly, and .edge output files
        // "n" outputs a list of triangles neighboring each triangle.
        // "e" outputs edge list (i.e. all the "connectivities")
        char options1[] = "Pne";
        triangulate(options1, &in_mid, &out_mid, &vorout_mid);
        //print out the mesh generation time
        std::cout << "mesh generation time = " << Timer::ElapsedSeconds(inital_time) << std::endl;
        //number of newly generated triangles
        unsigned int el_number=out_mid.numberoftriangles;

        //PASSING THE CREATED ELEMENTS TO THE ALPHA-SHAPE
        std::vector<int> preserved_list1(el_number);
        preserved_list1.resize(el_number);

        int number_of_preserved_elems= PassAlphaShape(ThisModelPart, preserved_list1, el_number, my_alpha, out_mid);

        //freeing memory

        clean_triangulateio(in_mid);
        clean_triangulateio(vorout_mid);
        KRATOS_WATCH("ln367");
        //NOW WE SHALL PERFORM ADAPTIVE REMESHING, i.e. insert and remove nodes based upon mesh quality
        // and prescribed sizes
        struct triangulateio in2;
        struct triangulateio out2;
        struct triangulateio vorout2;

        initialize_triangulateio(in2);
        initialize_triangulateio(out2);
        initialize_triangulateio(vorout2);

//			in2.firstnumber = 1;
        in2.numberofpoints = ThisModelPart.Nodes().size();
        in2.pointlist = (REAL *) malloc(in2.numberofpoints * 2 * sizeof(REAL));

        //writing the input point list
        for(unsigned int i = 0; i<ThisModelPart.Nodes().size(); i++)
        {
            int base = i*2;
            in2.pointlist[base] = (nodes_begin + i)->X();
            in2.pointlist[base+1] = (nodes_begin + i)->Y();
        }
        in2.numberoftriangles=number_of_preserved_elems;

        in2.trianglelist = (int *) malloc(in2.numberoftriangles * 3 * sizeof(int));
        in2.trianglearealist = (REAL *) malloc(in2.numberoftriangles * sizeof(REAL));
//			in2.trianglelist = new int[in2.numberoftriangles * 3];
//			in2.trianglearealist = new double[in2.numberoftriangles];

        KRATOS_WATCH(el_number);
        int counter = 0;
        //here I will assign a huge number of NODAL_H to the free surface nodes, so that there no nodes will be added
        for(ModelPart::NodesContainerType::iterator i_node = ThisModelPart.NodesBegin() ; i_node != ThisModelPart.NodesEnd() ; i_node++)
        {
            if (i_node->FastGetSolutionStepValue(IS_FREE_SURFACE)!=0)
            {

                double& val=i_node->FastGetSolutionStepValue(NODAL_H);
                val*=2.0;
                //i_node->FastGetSolutionStepValue(NODAL_H,1)=val;
                //KRATOS_WATCH("AAAAAAAAAAAAAAAAAAAAAAAA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            }
        }


        for (unsigned int el=0; el<el_number; el++)
        {
            if( static_cast<bool>(preserved_list1[el]) == true )
            {
                //saving the list of ONLY preserved triangles, the ones that passed alpha-shape check
                int new_base = counter*3;
                int old_base = el*3;
                //copying in case it is preserved
                in2.trianglelist[new_base] = out_mid.trianglelist[old_base];
                in2.trianglelist[new_base+1] = out_mid.trianglelist[old_base+1];
                in2.trianglelist[new_base+2] = out_mid.trianglelist[old_base+2];

                //calculate the prescribed h
                double prescribed_h = (nodes_begin + out_mid.trianglelist[old_base]-1)->FastGetSolutionStepValue(NODAL_H);
                prescribed_h += (nodes_begin + out_mid.trianglelist[old_base+1]-1)->FastGetSolutionStepValue(NODAL_H);
                prescribed_h += (nodes_begin + out_mid.trianglelist[old_base+2]-1)->FastGetSolutionStepValue(NODAL_H);
                prescribed_h *= 0.3333;
                //if h is the height of a equilateral triangle, the area is 1/2*h*h
                in2.trianglearealist[counter] = 0.5*(1.5*prescribed_h*1.5*prescribed_h);
                counter += 1;
            }

        }
        //now I set back the nodal_h
        for(ModelPart::NodesContainerType::iterator i_node = ThisModelPart.NodesBegin() ; i_node != ThisModelPart.NodesEnd() ; i_node++)
        {
            if (i_node->FastGetSolutionStepValue(IS_FREE_SURFACE)!=0)
            {
                double& nodal_h=i_node->FastGetSolutionStepValue(NODAL_H);
                nodal_h/=2.0;
            }
        }


        clean_triangulateio(out_mid);
        KRATOS_WATCH("ln420");
        //here we generate a new mesh adding/removing nodes, by initializing "q"-quality mesh and "a"-area constraint switches
        //
        // MOST IMPORTANT IS "r" switch, that refines previously generated mesh!!!!!!!!!!(that is the one given inside in2)
        //char mesh_regen_opts[] = "YYJaqrn";
        //char mesh_regen_opts[] = "YJq1.4arn";
        if (add_nodes==true)
        {
            char mesh_regen_opts[] = "YJq1.4arn";
            triangulate(mesh_regen_opts, &in2, &out2, &vorout2);
            KRATOS_WATCH("Adaptive remeshing executed")
        }
        else
        {
            char mesh_regen_opts[] = "YJrn";
            triangulate(mesh_regen_opts, &in2, &out2, &vorout2);
            KRATOS_WATCH("Non-Adaptive remeshing executed")
        }

        //and now we shall find out where the new nodes belong to
        //defintions for spatial search
        typedef Node<3> PointType;
        typedef Node<3>::Pointer PointPointerType;
        typedef std::vector<PointType::Pointer>           PointVector;
        //typedef std::vector<PointType::Pointer>::iterator PointIterator;
        typedef PointVector::iterator PointIterator;
        typedef std::vector<double>               DistanceVector;
        typedef std::vector<double>::iterator     DistanceIterator;
        KRATOS_WATCH("ln449");

        typedef Bucket<3, PointType, PointVector, PointPointerType, PointIterator, DistanceIterator > BucketType;

        typedef Tree< KDTreePartition<BucketType> > kd_tree; //Kdtree;

        //int step_data_size = ThisModelPart.GetNodalSolutionStepDataSize();

        //creating an auxiliary list for the new nodes
        PointVector list_of_new_nodes;

        //node to get the DOFs from
        Node<3>::DofsContainerType& reference_dofs = (ThisModelPart.NodesBegin())->GetDofs();

        double z = 0.0;
        int n_points_before_refinement = in2.numberofpoints;
        //if points were added, we add them as nodes to the ModelPart
        if (out2.numberofpoints > n_points_before_refinement )
        {
            for(int i = n_points_before_refinement; i<out2.numberofpoints; i++)
            {
                int id=i+1;
                int base = i*2;
                double& x= out2.pointlist[base];
                double& y= out2.pointlist[base+1];

                Node<3>::Pointer pnode = ThisModelPart.CreateNewNode(id,x,y,z);

                pnode->SetBufferSize(ThisModelPart.NodesBegin()->GetBufferSize() );

                list_of_new_nodes.push_back( pnode );

                //std::cout << "new node id = " << pnode->Id() << std::endl;
                //generating the dofs
                for(Node<3>::DofsContainerType::iterator iii = reference_dofs.begin();    iii != reference_dofs.end(); iii++)
                {
                    Node<3>::DofType& rDof = **iii;
                    Node<3>::DofType::Pointer p_new_dof = pnode->pAddDof( rDof );

                    (p_new_dof)->FreeDof();
//                                                (p_new_dof)->EquationId() = -1;

                }

            }
        }

        std::cout << "During refinement we added " << out2.numberofpoints-n_points_before_refinement<< " nodes " <<std::endl;
        //unsigned int bucket_size = 20;
        //performing the interpolation - all of the nodes in this list will be preserved
        //unsigned int max_results = 50;
        //PointVector results(max_results);
        //DistanceVector results_distances(max_results);
        array_1d<double,3> N;

        array_1d<double,3> x1,x2,x3,xc;

        //int number_of_preserved_elems=0;

        int point_base;
        //WHAT ARE THOSE????
// 			Node<3> work_point(0,0.0,0.0,0.0);
        unsigned int MaximumNumberOfResults = list_of_new_nodes.size();
        PointVector Results(MaximumNumberOfResults);
        DistanceVector ResultsDistances(MaximumNumberOfResults);

        step_data_size = ThisModelPart.GetNodalSolutionStepDataSize();

        if(out2.numberofpoints-n_points_before_refinement > 0) //if we added points
        {

            kd_tree  nodes_tree2(list_of_new_nodes.begin(),list_of_new_nodes.end(),bucket_size);
            nodes_begin = ThisModelPart.NodesBegin();

            for(int el = 0; el< in2.numberoftriangles; el++)
            {
                int base = el * 3;
                //coordinates
                point_base = (in2.trianglelist[base] - 1)*2;
                x1[0] = in2.pointlist[point_base];
                x1[1] = in2.pointlist[point_base+1];

                point_base = (in2.trianglelist[base+1] - 1)*2;
                x2[0] = in2.pointlist[point_base];
                x2[1] = in2.pointlist[point_base+1];

                point_base = (in2.trianglelist[base+2] - 1)*2;
                x3[0] = in2.pointlist[point_base];
                x3[1] = in2.pointlist[point_base+1];


                //find the center and "radius" of the element
                double xc,  yc, radius;
                CalculateCenterAndSearchRadius( x1[0], x1[1],
                                                x2[0], x2[1],
                                                x3[0], x3[1],
                                                xc,yc,radius);

                //find all of the new nodes within the radius
                int number_of_points_in_radius;
                work_point.X() = xc;
                work_point.Y() = yc;
                work_point.Z() = 0.0;

                number_of_points_in_radius = nodes_tree2.SearchInRadius(work_point, radius*1.01, Results.begin(),
                                             ResultsDistances.begin(),  MaximumNumberOfResults);

                Triangle2D3<Node<3> > geom(
                    *( (nodes_begin +  in2.trianglelist[base]-1).base() 	),
                    *( (nodes_begin +  in2.trianglelist[base+1]-1).base() 	),
                    *( (nodes_begin +  in2.trianglelist[base+2]-1).base() 	)
                );

                //check if inside and eventually interpolate
                for( PointIterator it_found = Results.begin(); it_found != Results.begin() + number_of_points_in_radius; it_found++)
                {
                    bool is_inside = false;
                    is_inside = CalculatePosition(x1[0], x1[1],
                                                  x2[0], x2[1],
                                                  x3[0], x3[1],
                                                  (*it_found)->X(),(*it_found)->Y(),N);


                    if(is_inside == true)
                    {
                        Interpolate(  geom,  N, step_data_size, *(it_found ) );

                    }
                }
            }
        }

        ThisModelPart.Elements().clear();

        //set the coordinates to the original value
        for( PointVector::iterator it =  list_of_new_nodes.begin(); it!=list_of_new_nodes.end(); it++)
        {
            const array_1d<double,3>& disp = (*it)->FastGetSolutionStepValue(DISPLACEMENT);
            (*it)->X0() = (*it)->X() - disp[0];
            (*it)->Y0() = (*it)->Y() - disp[1];
            (*it)->Z0() = 0.0;
        }
        //cleaning unnecessary data
        //in2.deinitialize();
        //in2.initialize();
        //free( in2.pointlist );

        //add preserved elements to the kratos
        Properties::Pointer properties = ThisModelPart.GetMesh().pGetProperties(1);
        nodes_begin = ThisModelPart.NodesBegin();
        (ThisModelPart.Elements()).reserve(out2.numberoftriangles);

        for(int iii = 0; iii< out2.numberoftriangles; iii++)
        {
            int id = iii + 1;
            int base = iii * 3;
            Triangle2D3<Node<3> > geom(
                *( (nodes_begin +  out2.trianglelist[base]-1).base() 	),
                *( (nodes_begin +  out2.trianglelist[base+1]-1).base() 	),
                *( (nodes_begin +  out2.trianglelist[base+2]-1).base() 	)
            );


#ifdef _DEBUG
            ModelPart::NodesContainerType& ModelNodes = ThisModelPart.Nodes();
            if( *(ModelNodes).find( out2.trianglelist[base]).base() == *(ThisModelPart.Nodes().end()).base() )
                KRATOS_THROW_ERROR(std::logic_error,"trying to use an inexisting node","");
            if( *(ModelNodes).find( out2.trianglelist[base+1]).base() == *(ThisModelPart.Nodes().end()).base() )
                KRATOS_THROW_ERROR(std::logic_error,"trying to use an inexisting node","");
            if( *(ModelNodes).find( out2.trianglelist[base+2]).base() == *(ThisModelPart.Nodes().end()).base() )
                KRATOS_THROW_ERROR(std::logic_error,"trying to use an inexisting node","");

#endif

            Element::Pointer p_element = rReferenceElement.Create(id, geom, properties);
            (ThisModelPart.Elements()).push_back(p_element);

        }
        ThisModelPart.Elements().Sort();

        //filling the neighbour list
        ModelPart::ElementsContainerType::const_iterator el_begin = ThisModelPart.ElementsBegin();
        for(ModelPart::ElementsContainerType::const_iterator iii = ThisModelPart.ElementsBegin();
                iii != ThisModelPart.ElementsEnd(); iii++)
        {
            //Geometry< Node<3> >& geom = iii->GetGeometry();
            int base = ( iii->Id() - 1 )*3;

            (iii->GetValue(NEIGHBOUR_ELEMENTS)).resize(3);
            GlobalPointersVector< Element >& neighb = iii->GetValue(NEIGHBOUR_ELEMENTS);
            for(int i = 0; i<3; i++)
            {
                int index = out2.neighborlist[base+i];
                if(index > 0)
                    neighb(i) = GlobalPointer<Element>(&*(el_begin + index-1));
                else
                    neighb(i) = Element::WeakPointer();
            }
        }
        //identifying boundary, creating skin
        IdentifyBoundary(ThisModelPart, rReferenceBoundaryCondition, properties, out2);

        KRATOS_WATCH("ln749");

        clean_triangulateio(in2);
        KRATOS_WATCH("ln752");
        clean_triangulateio(out2);
        KRATOS_WATCH("ln754");
        clean_triangulateio(vorout2);
        KRATOS_WATCH("Finished remeshing with Trigen_Droplet_Refine")

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
    //KdtreeType mKdtree;

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
    boost::numeric::ublas::bounded_matrix<double,2,2> mJ; //local jacobian
    boost::numeric::ublas::bounded_matrix<double,2,2> mJinv; //inverse jacobian
    array_1d<double,2> mC; //center pos
    array_1d<double,2> mRhs; //center pos
    //NodeEraseProcess* mpNodeEraseProcess;


    ///@}
    ///@name Private Operators
    ///@{
    void RemoveCloseNodes(ModelPart& ThisModelPart, KdtreeType& nodes_tree1, NodeEraseProcess& node_erase, double& h_factor)
    {
        //unsigned int bucket_size = 20;

        //performing the interpolation - all of the nodes in this list will be preserved
        unsigned int max_results = 100;
        //PointerVector<PointType> res(max_results);
        //NodeIterator res(max_results);
        PointVector res(max_results);
        DistanceVector res_distances(max_results);
        Node<3> work_point(0,0.0,0.0,0.0);

        unsigned int n_points_in_radius;
        //radius means the distance, closer than which no node shall be allowd. if closer -> mark for erasing
        double radius;

        for(ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin(); in != ThisModelPart.NodesEnd(); in++)
        {
            radius=h_factor*in->FastGetSolutionStepValue(NODAL_H);

            work_point[0]=in->X();
            work_point[1]=in->Y();
            work_point[2]=in->Z();

            n_points_in_radius = nodes_tree1.SearchInRadius(work_point, radius, res.begin(),res_distances.begin(), max_results);
            if (n_points_in_radius>1)
            {
                if (in->FastGetSolutionStepValue(IS_BOUNDARY)==0.0 && in->FastGetSolutionStepValue(IS_STRUCTURE)==0.0)
                {
                    //look if we are already erasing any of the other nodes
                    double erased_nodes = 0;
                    for(PointIterator i=res.begin(); i!=res.begin() + n_points_in_radius ; i++)
                        erased_nodes += in->Is(TO_ERASE);

                    if( erased_nodes < 1) //we cancel the node if no other nodes are being erased
                        in->Set(TO_ERASE,true);

                }
                else if ( (in)->FastGetSolutionStepValue(IS_STRUCTURE)!=1.0) //boundary nodes will be removed if they get REALLY close to another boundary node (0.2 * h_factor)
                {
                    //here we loop over the neighbouring nodes and if there are nodes
                    //with IS_BOUNDARY=1 which are closer than 0.2*nodal_h from our we remove the node we are considering
                    unsigned int k = 0;
                    unsigned int counter = 0;
                    for(PointIterator i=res.begin(); i!=res.begin() + n_points_in_radius ; i++)
                    {
                        if ( (*i)->FastGetSolutionStepValue(IS_BOUNDARY,1)==1.0 && res_distances[k] < 0.2*radius && res_distances[k] > 0.0 )
                        {
                            // 										KRATOS_WATCH( res_distances[k] );
                            counter += 1;
                        }
                        k++;
                    }
                    if(counter > 0)
                        in->Set(TO_ERASE,true);
                }
            }

        }


        node_erase.Execute();

        KRATOS_WATCH("Number of nodes after erasing")
        KRATOS_WATCH(ThisModelPart.Nodes().size())
    }


    int PassAlphaShape(ModelPart& ThisModelPart, std::vector<int>& preserved_list1, unsigned int & el_number, double& my_alpha, struct triangulateio& out_mid)
    {
        //NOTE THAT preserved_list1 will be overwritten, only the elements that passed alpha-shaoe check will enter it

        //prepairing for alpha shape passing : creating necessary arrays
        //list of preserved elements is created: at max el_number can be preserved (all elements)


        array_1d<double,3> x1,x2,x3,xc;

        //int number_of_preserved_elems=0;
        int number_of_preserved_elems=0;
        int point_base;
        //loop for passing alpha shape
        for(unsigned int el = 0; el< el_number; el++)
        {
            int base = el * 3;

            //coordinates
            point_base = (out_mid.trianglelist[base] - 1)*2;
            x1[0] = out_mid.pointlist[point_base];
            x1[1] = out_mid.pointlist[point_base+1];

            point_base = (out_mid.trianglelist[base+1] - 1)*2;
            x2[0] = out_mid.pointlist[point_base];
            x2[1] = out_mid.pointlist[point_base+1];

            point_base = (out_mid.trianglelist[base+2] - 1)*2;
            x3[0] = out_mid.pointlist[point_base];
            x3[1] = out_mid.pointlist[point_base+1];

            //here we shall temporarily save the elements and pass them afterwards to the alpha shape
            Geometry<Node<3> > temp;

            temp.push_back( *((ThisModelPart.Nodes()).find( out_mid.trianglelist[base]).base()	) );
            temp.push_back( *((ThisModelPart.Nodes()).find( out_mid.trianglelist[base+1]).base()	) );
            temp.push_back( *((ThisModelPart.Nodes()).find( out_mid.trianglelist[base+2]).base()	) );

            int number_of_structure_nodes = int( temp[0].FastGetSolutionStepValue(IS_STRUCTURE) );
            number_of_structure_nodes += int( temp[1].FastGetSolutionStepValue(IS_STRUCTURE) );
            number_of_structure_nodes += int( temp[2].FastGetSolutionStepValue(IS_STRUCTURE) );

            //check the number of nodes of boundary
            int nfs = int( temp[0].FastGetSolutionStepValue(IS_FREE_SURFACE) );
            nfs += int( temp[1].FastGetSolutionStepValue(IS_FREE_SURFACE) );
            nfs += int( temp[2].FastGetSolutionStepValue(IS_FREE_SURFACE) );

            //check the number of nodes of boundary
            int nfluid = int( temp[0].FastGetSolutionStepValue(IS_FLUID) );
            nfluid += int( temp[1].FastGetSolutionStepValue(IS_FLUID) );
            nfluid += int( temp[2].FastGetSolutionStepValue(IS_FLUID) );

            //check the number of nodes of boundary
            int nboundary = int( temp[0].FastGetSolutionStepValue(IS_BOUNDARY) );
            nboundary += int( temp[1].FastGetSolutionStepValue(IS_BOUNDARY) );
            nboundary += int( temp[2].FastGetSolutionStepValue(IS_BOUNDARY) );
            //first check that we are working with fluid elements, otherwise throw an error
            //if (nfluid!=3)
            //	KRATOS_ERROR(std::logic_error,"THATS NOT FLUID or NOT TRIANGLE!!!!!! ERROR","");
            //otherwise perform alpha shape check


            if(number_of_structure_nodes!=3) //if it is = 3 it is a completely fixed element -> do not add it
            {

                if (nfs != 0 || nfluid != 3)  //in this case it is close to the surface so i should use alpha shape
                {

                    if( (AlphaShape(my_alpha,temp) && number_of_structure_nodes!=3)) //if alpha shape says to preserve
                    {
// 							if(nboundary==3 && number_of_structure_nodes > 1 && nfs > 0) //if it is = 3 pressure problems -> do not add it
// 							{
// 							      preserved_list1[el] = false;
// 							}
// 							else
// 							{
                        preserved_list1[el] = true;
                        number_of_preserved_elems += 1;
// 							}
                    }
//                     if (nfs == 2) //elements at boundary should ALWAYS be preserved in surface tension problems!!
// 		    {
// 			preserved_list1[el] = true;
// 			number_of_preserved_elems += 1;
// 		    }
                }
                else //internal triangle --- should be ALWAYS preserved
                {
                    double bigger_alpha = my_alpha*10.0;
                    if( AlphaShape(bigger_alpha,temp) && number_of_structure_nodes!=3)
                    {
// 						  	if(nboundary==3 && number_of_structure_nodes > 1 && nfs > 0) //if it is = 3 pressure problems -> do not add it
// 							{
// 							      preserved_list1[el] = false;
// 							}
// 							else
// 							{
                        preserved_list1[el] = true;
                        number_of_preserved_elems += 1;
// 							}
                    }
                }
            }
            else
                preserved_list1[el] = false;

        }
        return number_of_preserved_elems;
    }





    void IdentifyBoundary(ModelPart& ThisModelPart, Condition const& rReferenceBoundaryCondition,  Properties::Pointer& properties, struct triangulateio& out2)
    {

        //reset the boundary flag
        for(ModelPart::NodesContainerType::const_iterator in = ThisModelPart.NodesBegin(); in!=ThisModelPart.NodesEnd(); in++)
        {
            in->FastGetSolutionStepValue(IS_BOUNDARY) = 0;
        }
        //filling the elemental neighbours list (from now on the elements list can not change)
        ModelPart::ElementsContainerType::iterator elements_end = ThisModelPart.Elements().end();

        ThisModelPart.Elements().Unique();

        //now the boundary faces
        for(ModelPart::ElementsContainerType::iterator iii = ThisModelPart.ElementsBegin();	iii != ThisModelPart.ElementsEnd(); iii++)
        {
            int base = ( iii->Id() - 1 )*3;

            ModelPart::ElementsContainerType::iterator el_neighb;
            /*each face is opposite to the corresponding node number so
             0 ----- 1 2
             1 ----- 2 0
             2 ----- 0 1
            */

            ////finding boundaries and creating the "skin"
            //
            //********************************************************************
            //first face
            el_neighb = (ThisModelPart.Elements()).find( out2.neighborlist[base] );
            if( el_neighb == elements_end )
            {
                //std::cout << "node0" << std::endl;
                //if no neighnour is present => the face is free surface
                iii->GetGeometry()[1].FastGetSolutionStepValue(IS_BOUNDARY) = 1;
                iii->GetGeometry()[2].FastGetSolutionStepValue(IS_BOUNDARY) = 1;

                //Generate condition
                Condition::NodesArrayType temp1;
                temp1.reserve(2);
                temp1.push_back(iii->GetGeometry()(1));
                temp1.push_back(iii->GetGeometry()(2));

                Geometry< Node<3> >::Pointer cond = Geometry< Node<3> >::Pointer(new Geometry< Node<3> >(temp1) );
                int id = (iii->Id()-1)*3;
                Condition::Pointer p_cond = rReferenceBoundaryCondition.Create(id, temp1, properties);
                ThisModelPart.Conditions().push_back(p_cond);

            }

            //********************************************************************
            //second face
            el_neighb = (ThisModelPart.Elements()).find( out2.neighborlist[base+1] );
            //if( el != ThisModelPart.Elements().end() )
            if( el_neighb == elements_end )
            {
                //if no neighnour is present => the face is free surface
                iii->GetGeometry()[2].FastGetSolutionStepValue(IS_BOUNDARY) = 1;
                iii->GetGeometry()[0].FastGetSolutionStepValue(IS_BOUNDARY) = 1;

                //Generate condition
                Condition::NodesArrayType temp1;
                temp1.reserve(2);
                temp1.push_back(iii->GetGeometry()(2));
                temp1.push_back(iii->GetGeometry()(0));

                Geometry< Node<3> >::Pointer cond = Geometry< Node<3> >::Pointer(new Geometry< Node<3> >(temp1) );
                int id = (iii->Id()-1)*3+1;
                Condition::Pointer p_cond = rReferenceBoundaryCondition.Create(id, temp1, properties);
                ThisModelPart.Conditions().push_back(p_cond);


            }

            //********************************************************************
            //third face
            el_neighb = (ThisModelPart.Elements()).find( out2.neighborlist[base+2] );
            if( el_neighb == elements_end )
            {
                //if no neighnour is present => the face is free surface
                iii->GetGeometry()[0].FastGetSolutionStepValue(IS_BOUNDARY) = 1;
                iii->GetGeometry()[1].FastGetSolutionStepValue(IS_BOUNDARY) = 1;

//					Generate condition
                Condition::NodesArrayType temp1;
                temp1.reserve(2);
                temp1.push_back(iii->GetGeometry()(0));
                temp1.push_back(iii->GetGeometry()(1));

                Geometry< Node<3> >::Pointer cond = Geometry< Node<3> >::Pointer(new Geometry< Node<3> >(temp1) );
                int id = (iii->Id()-1)*3+2;
                Condition::Pointer p_cond = rReferenceBoundaryCondition.Create(id, temp1, properties);
                ThisModelPart.Conditions().push_back(p_cond);

            }


        }
    }


    //returns false if it should be removed
    bool AlphaShape(double alpha_param, Geometry<Node<3> >& pgeom)
    //bool AlphaShape(double alpha_param, Triangle2D<Node<3> >& pgeom)
    {
        KRATOS_TRY


        double x0 = pgeom[0].X();
        double x1 = pgeom[1].X();
        double x2 = pgeom[2].X();

        double y0 = pgeom[0].Y();
        double y1 = pgeom[1].Y();
        double y2 = pgeom[2].Y();

        mJ(0,0)=2.0*(x1-x0);
        mJ(0,1)=2.0*(y1-y0);
        mJ(1,0)=2.0*(x2-x0);
        mJ(1,1)=2.0*(y2-y0);


        double detJ = mJ(0,0)*mJ(1,1)-mJ(0,1)*mJ(1,0);

        mJinv(0,0) =  mJ(1,1);
        mJinv(0,1) = -mJ(0,1);
        mJinv(1,0) = -mJ(1,0);
        mJinv(1,1) =  mJ(0,0);

        bounded_matrix<double,2,2> check;

//calculate average h
        double h;
        h =  pgeom[0].FastGetSolutionStepValue(NODAL_H);
        h += pgeom[1].FastGetSolutionStepValue(NODAL_H);
        h += pgeom[2].FastGetSolutionStepValue(NODAL_H);
        h *= 0.333333333;


        if(detJ < 5e-3*h*h)
        {
            //std::cout << "detJ = " << detJ << std::endl;
            ////mark as boundary
            pgeom[0].GetSolutionStepValue(IS_BOUNDARY) = 1;
            pgeom[1].GetSolutionStepValue(IS_BOUNDARY) = 1;
            pgeom[2].GetSolutionStepValue(IS_BOUNDARY) = 1;
            return false;
        }

        else
        {

            double x0_2 = x0*x0 + y0*y0;
            double x1_2 = x1*x1 + y1*y1;
            double x2_2 = x2*x2 + y2*y2;

            //finalizing the calculation of the inverted matrix
            //std::cout<<"MATR INV"<<MatrixInverse(mJ)<<std::endl;
            mJinv /= detJ;
            //calculating the RHS
            mRhs[0] = (x1_2 - x0_2);
            mRhs[1] = (x2_2 - x0_2);

            //calculate position of the center
            noalias(mC) = prod(mJinv,mRhs);

            double radius = sqrt(pow(mC[0]-x0,2)+pow(mC[1]-y0,2));


            if (radius < h*alpha_param)
            {
                return true;
            }
            else
            {
                return false;
            }
        }


        KRATOS_CATCH("")
    }
    //AUXILLIARY FCTNS
    inline void CalculateCenterAndSearchRadius(const double x0, const double y0,
            const double x1, const double y1,
            const double x2, const double y2,
            double& xc, double& yc, double& R
                                              )
    {
        xc = 0.3333333333333333333*(x0+x1+x2);
        yc = 0.3333333333333333333*(y0+y1+y2);

        double R1 = (xc-x0)*(xc-x0) + (yc-y0)*(yc-y0);
        double R2 = (xc-x1)*(xc-x1) + (yc-y1)*(yc-y1);
        double R3 = (xc-x2)*(xc-x2) + (yc-y2)*(yc-y2);

        R = R1;
        if(R2 > R) R = R2;
        if(R3 > R) R = R3;

        R = sqrt(R);
    }

    inline double CalculateVol(	const double x0, const double y0,
                                const double x1, const double y1,
                                const double x2, const double y2
                              )
    {
        return 0.5*( (x1-x0)*(y2-y0)- (y1-y0)*(x2-x0) );
    }

    inline bool CalculatePosition(	const double x0, const double y0,
                                    const double x1, const double y1,
                                    const double x2, const double y2,
                                    const double xc, const double yc,
                                    array_1d<double,3>& N
                                 )
    {
        double area = CalculateVol(x0,y0,x1,y1,x2,y2);

        if(area < 0.000000000000001)
        {
            KRATOS_THROW_ERROR(std::logic_error,"element with zero area found","");
        }



        N[0] = CalculateVol(x1,y1,x2,y2,xc,yc)  / area;
        N[1] = CalculateVol(x2,y2,x0,y0,xc,yc)  / area;
        N[2] = CalculateVol(x0,y0,x1,y1,xc,yc)  / area;

        /*			  N[0] = CalculateVol(x0,y0,x1,y1,xc,yc) * inv_area;
        			N[1] = CalculateVol(x1,y1,x2,y2,xc,yc) * inv_area;
        			N[2] = CalculateVol(x2,y2,x0,y0,xc,yc) * inv_area;*/
        double tol = 1e-4;
        double upper_limit = 1.0+tol;
        double lower_limit = -tol;

        if(N[0] >= lower_limit && N[1] >= lower_limit && N[2] >= lower_limit && N[0] <= upper_limit && N[1] <= upper_limit && N[2] <= upper_limit) //if the xc yc is inside the triangle
            //if(N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0) //if the xc yc is inside the triangle return true
            return true;

        return false;
    }

    void Interpolate( Triangle2D3<Node<3> >& geom, const array_1d<double,3>& N,
                      unsigned int step_data_size,
                      Node<3>::Pointer pnode)
    {
        unsigned int buffer_size = pnode->GetBufferSize();
        //KRATOS_WATCH("Buffer size")
        //KRATOS_WATCH(buffer_size)

        for(unsigned int step = 0; step<buffer_size; step++)
        {

            //getting the data of the solution step
            double* step_data = (pnode)->SolutionStepData().Data(step);


            double* node0_data = geom[0].SolutionStepData().Data(step);
            double* node1_data = geom[1].SolutionStepData().Data(step);
            double* node2_data = geom[2].SolutionStepData().Data(step);

            //copying this data in the position of the vector we are interested in
            for(unsigned int j= 0; j<step_data_size; j++)
            {

                step_data[j] = N[0]*node0_data[j] + N[1]*node1_data[j] + N[2]*node2_data[j];


            }
        }
        //20150916 ajarauta
        pnode->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET) = 0.0;
	pnode->FastGetSolutionStepValue(TRIPLE_POINT) = 0.0;
	/*
	if (pnode->FastGetSolutionStepValue(IS_BOUNDARY) == 0.0)
	{
	    pnode->FastGetSolutionStepValue(TRIPLE_POINT) = 0.0;
	    pnode->FastGetSolutionStepValue(FORCE_X) = 0.0;
	    pnode->FastGetSolutionStepValue(FORCE_Y) = 0.0;
	}
	*/

        if (N[0]==0.0 && N[1]==0.0 && N[2]==0.0)
            KRATOS_THROW_ERROR(std::logic_error,"SOMETHING's wrong with the added nodes!!!!!! ERROR","");

        //if ( pnode->FastGetSolutionStepValue(BULK_MODULUS)==0.0)
        //		KRATOS_ERROR(std::logic_error,"SOMETHING's wrong with the added nodes!!!!!! ERROR","");

        //now we assure that the flag variables are set coorect!! since we add nodes inside of the fluid volume only
        //we manually reset the IS_BOUNDARY, IS_FLUID, IS_STRUCTURE, IS_FREE_SURFACE values in a right way
        //not to have values, like 0.33 0.66 resulting if we would have been interpolating them in the same way
        //as the normal variables, like Velocity etc


        pnode->FastGetSolutionStepValue(IS_BOUNDARY)=0.0;
        pnode->FastGetSolutionStepValue(IS_STRUCTURE)=0.0;
        pnode->FastGetSolutionStepValue(IS_INTERFACE)=0.0;
        pnode->Set(TO_ERASE,false);
        pnode->FastGetSolutionStepValue(IS_FREE_SURFACE)=0.0;
        pnode->FastGetSolutionStepValue(IS_FLUID)=1.0;
    }

    void initialize_triangulateio( triangulateio& tr )
    {
        tr.pointlist                  = (REAL*) NULL;
        tr.pointattributelist         = (REAL*) NULL;
        tr.pointmarkerlist            = (int*) NULL;
        tr.numberofpoints             = 0;
        tr.numberofpointattributes    = 0;
        tr.trianglelist               = (int*) NULL;
        tr.triangleattributelist      = (REAL*) NULL;
        tr.trianglearealist           = (REAL*) NULL;
        tr.neighborlist               = (int*) NULL;
        tr.numberoftriangles          = 0;
        tr.numberofcorners            = 3;
        tr.numberoftriangleattributes = 0;
        tr.segmentlist                = (int*) NULL;
        tr.segmentmarkerlist          = (int*) NULL;
        tr.numberofsegments           = 0;
        tr.holelist                   = (REAL*) NULL;
        tr.numberofholes              = 0;
        tr.regionlist                 = (REAL*) NULL;
        tr.numberofregions            = 0;
        tr.edgelist                   = (int*) NULL;
        tr.edgemarkerlist             = (int*) NULL;
        tr.normlist                   = (REAL*) NULL;
        tr.numberofedges              = 0;
    };

    void clean_triangulateio( triangulateio& tr )
    {
        if(tr.pointlist != NULL) free(tr.pointlist );
        if(tr.pointattributelist != NULL) free(tr.pointattributelist );
        if(tr.pointmarkerlist != NULL) free(tr.pointmarkerlist   );
        if(tr.trianglelist != NULL) free(tr.trianglelist  );
        if(tr.triangleattributelist != NULL) free(tr.triangleattributelist );
        if(tr.trianglearealist != NULL) free(tr.trianglearealist );
        if(tr.neighborlist != NULL) free(tr.neighborlist   );
        if(tr.segmentlist != NULL) free(tr.segmentlist    );
        if(tr.segmentmarkerlist != NULL) free(tr.segmentmarkerlist   );
        if(tr.holelist != NULL) free(tr.holelist      );
        if(tr.regionlist != NULL) free(tr.regionlist  );
        if(tr.edgelist != NULL) free(tr.edgelist   );
        if(tr.edgemarkerlist != NULL) free(tr.edgemarkerlist   );
        if(tr.normlist != NULL) free(tr.normlist  );
    };

    void InterpolateOnEdge( Geometry<Node<3> >& geom, int point1, int point2, unsigned int step_data_size, Node<3>::Pointer pnode)
    {
		unsigned int buffer_size = pnode->GetBufferSize();
		KRATOS_INFO("TriGenDropletModeler") << "buffer_size: " << buffer_size << std::endl;

		for(unsigned int step = 0; step<buffer_size; step++) {
			//getting the data of the solution step
			double* step_data = (pnode)->SolutionStepData().Data(step);

			if (point1>2 || point1<0 || point2>2 || point2<0 || (point1==point2)) {
				KRATOS_INFO("TriGenDropletModeler") << "point1: " << point1 << std::endl;
				KRATOS_INFO("TriGenDropletModeler") << "point2: " << point2 << std::endl;
				KRATOS_ERROR << "THE EDGE POINTS ARE INVALID " << std::endl;
			}

			double* node0_data = geom[point1].SolutionStepData().Data(step);
			double* node1_data = geom[point2].SolutionStepData().Data(step);

			//copying this data in the position of the vector we are interested in
			for(unsigned int j= 0; j<step_data_size; j++) {
				step_data[j] = 0.5*(node0_data[j] + node1_data[j]);
			}
		}
		//now we assure that the flag variables are set coorect!! since we add nodes inside of the fluid volume only
		//we manually reset the IS_BOUNDARY, IS_FLUID, IS_STRUCTURE, IS_FREE_SURFACE values in a right way
		//not to have values, like 0.33 0.66 resulting if we would have been interpolating them in the same way
		//as the normal variables, like Velocity etc

		//pnode->FastGetSolutionStepValue(IS_BOUNDARY)=1.0;
		//pnode->FastGetSolutionStepValue(FLAG_VARIABLE)=1.0;

		pnode->FastGetSolutionStepValue(IS_LAGRANGIAN_INLET)=0.0;
		pnode->FastGetSolutionStepValue(IS_BOUNDARY)=1.0;
		pnode->FastGetSolutionStepValue(IS_FLUID)=1.0;
		pnode->FastGetSolutionStepValue(TRIPLE_POINT)=0.0;
		pnode->FastGetSolutionStepValue(CONTACT_ANGLE)=0.0;
		//pnode->FastGetSolutionStepValue(SOLID_FRACTION_GRADIENT_X)=0.0;

		pnode->Set(TO_ERASE,false);

		if (pnode->FastGetSolutionStepValue(IS_INTERFACE)>0.9)
			pnode->FastGetSolutionStepValue(IS_INTERFACE)=1.0;
		else
			pnode->FastGetSolutionStepValue(IS_INTERFACE)=0.0;
		if (pnode->FastGetSolutionStepValue(IS_FREE_SURFACE)>0.9)
			pnode->FastGetSolutionStepValue(IS_FREE_SURFACE)=1.0;
		else
			pnode->FastGetSolutionStepValue(IS_FREE_SURFACE)=0.0;
		if (pnode->FastGetSolutionStepValue(FLAG_VARIABLE)>0.9)
			pnode->FastGetSolutionStepValue(FLAG_VARIABLE)=1.0;
		else
			pnode->FastGetSolutionStepValue(FLAG_VARIABLE)=0.0;
		if (pnode->FastGetSolutionStepValue(IS_STRUCTURE)>0.9)
		{
			pnode->FastGetSolutionStepValue(IS_STRUCTURE)=1.0;
			pnode->FastGetSolutionStepValue(IS_INTERFACE)=0.0;
			pnode->FastGetSolutionStepValue(IS_FREE_SURFACE)=0.0;
		}
		else
			pnode->FastGetSolutionStepValue(IS_STRUCTURE)=0.0;

	 };


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
    TriGenDropletModeler& operator=(TriGenDropletModeler const& rOther);


    ///@}

}; // Class TriGenPFEMModeler

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  TriGenPFEMModeler& rThis);

/// output stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  TriGenDropletModeler& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const TriGenDropletModeler& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_TRIGEN_PFEM_MODELER_H_INCLUDED  defined
