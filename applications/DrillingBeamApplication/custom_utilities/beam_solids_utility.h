#ifndef BEAM_SOLIDS_UTILITY_H
#define BEAM_SOLIDS_UTILITY_H

#include <fstream>

#include "includes/define.h"
#include "includes/variables.h"
#include "includes/model_part.h"



namespace Kratos {

    class BeamSolidsUtility{

    public:

        KRATOS_CLASS_POINTER_DEFINITION(BeamSolidsUtility);

        typedef ModelPart::NodesContainerType NodesArrayType;

        unsigned int mNodesEdge, submodel_number_of_nodes_A1, submodel_number_of_nodes_A2;

        std::vector<int> IDs_edge;

        std::vector<int> IDs_A1;
        std::vector<int> IDs_A2;

        std::vector<int> sorted_labels_A1;
        std::vector<int> sorted_labels_A2;

        NodesArrayType* mpNodesPt;       

        // std::vector<array_1d<double,3>> initial_directive_line;
        // std::vector<array_1d<double,3>> directive_line;
        // std::vector<array_1d<double,3>> normals;       

        /////   Default constructor   /////

        // BeamSolidsUtility() {}

        /////   Constructor with ModelPart   /////

        BeamSolidsUtility(ModelPart& r_model_part) {

            ModelPart& EdgePoints = r_model_part.GetSubModelPart("GENERIC_SbModel_edge");
            pNodesPt = &(EdgePoints.GetCommunicator().LocalMesh().Nodes());

        }

        /////   Destructor   /////

        virtual ~BeamSolidsUtility() {}


// -------------------------------------------------------------------------------------------------
    void ComputeEdgePointOfBeamSolids(ModelPart& r_model_part){

        // ModelPart& EdgePoints = r_model_part.GetSubModelPart("GENERIC_SbModel_edge");
        // NodesArrayType& pNodesPt = EdgePoints.GetCommunicator().LocalMesh().Nodes();                 

        array_1d<double,3> coords1, coords2, displacements, edge_coords;

        Node < 3 > Aux_Node;

        std::ofstream RES_edgePoint("RES_edgePoint.txt", std::ios_base::out | std::ios_base::app);


        nodes_edge = pNodesPt->size();
        std::cout << "submodel_number_of_nodes ..... point " << nodes_edge << std::endl;


        ModelPart::NodeIterator i_iteratorPt1 = pNodesPt->ptr_begin();
        Aux_Node = *i_iteratorPt1;

        coords1 = Aux_Node.Coordinates();
        displacements = Aux_Node.FastGetSolutionStepValue(DISPLACEMENT);

        ModelPart::NodeIterator i_iteratorPt2 = pNodesPt->ptr_begin() + 1;
        Aux_Node = *i_iteratorPt2;

        coords2 = Aux_Node.Coordinates();
        displacements = Aux_Node.FastGetSolutionStepValue(DISPLACEMENT);

        edge_coords = 0.5*(coords1 + coords2);                

        RES_edgePoint << edge_coords[0] << " " << edge_coords[1] << " " << edge_coords[2] << "\n";

    }


// ------------------------------------------------------------------------------------------------
    void ComputeDirectiveLineOfBeamSolids(ModelPart& r_model_part){

        typedef ModelPart::NodesContainerType NodesArrayType;

        Node < 3 > Aux_Node;

        array_1d<double,3> coords_A1, coords_A2, directive_coords;
        array_1d<double,3> aux_vector_rot, aux_vector_last, aux_vector_delta, aux_vector_dir;

        unsigned int aux_label;

        double modulus, directive_length=0, dlength=0;

        std::ofstream RES_directive("RES_directive.txt", std::ios_base::out | std::ios_base::app);

        ModelPart& EdgeLineA1 = r_model_part.GetSubModelPart("GENERIC_SbModel_line_A1");
        ModelPart& EdgeLineA2 = r_model_part.GetSubModelPart("GENERIC_SbModel_line_A2");

        NodesArrayType& pNodesA1 = EdgeLineA1.GetCommunicator().LocalMesh().Nodes();
        NodesArrayType& pNodesA2 = EdgeLineA2.GetCommunicator().LocalMesh().Nodes();


        submodel_number_of_nodes_A1 = pNodesA1.size();
        std::cout << "submodel_number_of_nodes A1.................." << submodel_number_of_nodes_A1 << std::endl;

        submodel_number_of_nodes_A2 = pNodesA2.size();
        std::cout << "submodel_number_of_nodes A2.................." << submodel_number_of_nodes_A2 << std::endl;

        // KRATOS_ERROR_IF(pNodesA1.size() != pNodesA2.size())<<"The number of nodes is not consistent between edges"<<std::endl;

        for (unsigned int inode = 0; inode < submodel_number_of_nodes_A1; inode++) {

            ModelPart::NodeIterator i_iteratorA1 = pNodesA1.ptr_begin() + inode;
            Aux_Node = *i_iteratorA1;

            const int& aux_id = Aux_Node.Id(); 
            IDs_A1.push_back(aux_id);

        }

        for (unsigned int inode = 0; inode < submodel_number_of_nodes_A2; inode++) {

            ModelPart::NodeIterator i_iteratorA2 = pNodesA2.ptr_begin() + inode;
            Aux_Node = *i_iteratorA2;

            const int& aux_id = Aux_Node.Id(); 
            IDs_A2.push_back(aux_id);

        }

        // unsigned int values[] = { 20, 60, 40, 25, 15, 2 };
        // int compare(const void * a, const void * b ){
        //     return ( *(int*)a - *(int*)b );
        // }

        // qsort(values, 6, sizeof(int), compare);
        // std::cout << values << std::endl;

        sorted_labels_A1 = IDs_A1;
        sorted_labels_A2 = IDs_A2; 

        for (unsigned int inode = 0; inode < submodel_number_of_nodes_A1; inode++) {

            aux_label = sorted_labels_A1[inode];
            Aux_Node = r_model_part.GetNode(aux_label);
            coords_A1 = Aux_Node.Coordinates();

            aux_label = sorted_labels_A2[inode];
            Aux_Node = r_model_part.GetNode(aux_label);
            coords_A2 = Aux_Node.Coordinates();

            directive_coords = 0.5*(coords_A1 + coords_A2);                     

            RES_directive << directive_coords[0] << " " << directive_coords[1] << " " << directive_coords[2] << "\n";                         

        }

    //         int label_B = sorted_labels_B[inode];
    //         Aux_Node = r_model_part.GetNode(label_B);
    //         coords_B = Aux_Node.Coordinates();

    //         x2 = coords_B[0];
    //         y2 = coords_B[1];
    //         z2 = coords_B[2];

    //         aux_vector_rot[0] = (y1-y2)*(z3-z2) - (z1-z2)*(y3-y2);        // calculation of the normal to the deformed plane
    //         aux_vector_rot[1] = (z1-z2)*(x3-x2) - (x1-x2)*(z3-z2);
    //         aux_vector_rot[2] = (x1-x2)*(y3-y2) - (y1-y2)*(x3-x2);

    //         modulus = sqrt( pow(aux_vector_rot[0],2) + pow(aux_vector_rot[1],2) + pow(aux_vector_rot[2],2) );

    //         aux_vector_rot = aux_vector_rot / modulus;
    //         normals.push_back(aux_vector_rot);

    //         aux_vector_last = aux_vector_dir;
    //         aux_vector_dir = (coords_A + coords_B + coords_C)/3;          // Calculation of the directive line
    //         a = aux_vector_dir[0];
    //         b = aux_vector_dir[1];
    //         c = aux_vector_dir[2];

    //         RES_directive << a << " " << b << " " << c << "\n";
    //         directive_line.push_back(aux_vector_dir);

    //         aux_vector_delta = aux_vector_dir - aux_vector_last;
    //         dlength = sqrt( pow(aux_vector_delta[0],2)+pow(aux_vector_delta[1],2)+pow(aux_vector_delta[2],2) );
    //         directive_length += dlength;

    //     directive_line.clear();
    //     normals.clear();

        }        


// -------------------------------------------------------------------------------------------------
    void ComputeReactionsOfBeamSolids(ModelPart& r_model_part){

        typedef ModelPart::NodesContainerType NodesArrayType;

        array_1d<double,3> coords,reactions;

        Node < 3 > Aux_Node;

        double feed_pressure=0, rotation_pressure=0, area=0;

        std::ofstream RES_reactions("RES_reactions.txt", std::ios_base::out | std::ios_base::app);
        area = 0.587141e-3;


        ModelPart& Base = r_model_part.GetSubModelPart("DISPLACEMENT_BC");
        NodesArrayType& pNodes = Base.GetCommunicator().LocalMesh().Nodes();

        for (unsigned int inode = 0; inode < pNodes.size(); inode++) {

            ModelPart::NodeIterator i_iteratorA = pNodes.ptr_begin() + inode;
            Aux_Node = *i_iteratorA;

            reactions = Aux_Node.FastGetSolutionStepValue(REACTION);
            coords = Aux_Node.Coordinates();

            feed_pressure += reactions[0];
            rotation_pressure += reactions[1]*coords[2] - reactions[2]*coords[1];

        }

        feed_pressure /= area;
        rotation_pressure /= area;

        RES_reactions << feed_pressure << " " << rotation_pressure << "\n";

    }

    };

} // namespace Kratos.

#endif // BEAM_SOLIDS_UTILITY_H