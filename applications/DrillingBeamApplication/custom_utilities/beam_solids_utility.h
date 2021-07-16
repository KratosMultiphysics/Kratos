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

        NodesArrayType* mpNodesPt;
        NodesArrayType* mpNodesA1;
        NodesArrayType* mpNodesA2;
        NodesArrayType* mpNodesBC;

        /////   Default constructor   /////

        BeamSolidsUtility(ModelPart& r_model_part) {

            ModelPart& EdgeNodes = r_model_part.GetSubModelPart("GENERIC_SbModel_edge");
            mpNodesPt = &(EdgeNodes.GetCommunicator().LocalMesh().Nodes());

            ModelPart& EdgeLineA1 = r_model_part.GetSubModelPart("GENERIC_SbModel_line_A1");
            mpNodesA1 = &(EdgeLineA1.GetCommunicator().LocalMesh().Nodes());

            ModelPart& EdgeLineA2 = r_model_part.GetSubModelPart("GENERIC_SbModel_line_A2");
            mpNodesA2 = &(EdgeLineA2.GetCommunicator().LocalMesh().Nodes());

            ModelPart& BaseNodes = r_model_part.GetSubModelPart("DISPLACEMENT_BC");
            mpNodesBC = &(BaseNodes.GetCommunicator().LocalMesh().Nodes());

        }

        /////   Destructor   /////

        virtual ~BeamSolidsUtility() {}


// -------------------------------------------------------------------------------------------------
    void ComputeEdgePointOfBeamSolids(ModelPart& r_model_part){

        unsigned int aux_label;

        array_1d<double,3> coords1, coords2, displacements, edge_coords,edge_coords_C;
        double edge_coords_X, edge_coords_X0;

        Node < 3 > Aux_Node;

        std::ofstream RES_edgePoint("RES_edgePoint.txt", std::ios_base::out | std::ios_base::app);


        ModelPart::NodeIterator iteratorPt1 = mpNodesPt->ptr_begin();
        Aux_Node = *iteratorPt1;

        coords1 = Aux_Node.Coordinates();
        displacements = Aux_Node.FastGetSolutionStepValue(DISPLACEMENT);

        ModelPart::NodeIterator iteratorPt2 = mpNodesPt->ptr_begin() + 1;
        Aux_Node = *iteratorPt2;

        coords2 = Aux_Node.Coordinates();
        displacements = Aux_Node.FastGetSolutionStepValue(DISPLACEMENT);

        edge_coords = 0.5*(coords1 + coords2);                       

        RES_edgePoint << edge_coords[0] << " " << edge_coords[1] << " " << edge_coords[2] << "\n";

    }


// ------------------------------------------------------------------------------------------------
    void ComputeDirectiveLineOfBeamSolids(ModelPart& r_model_part){

        unsigned int submodel_number_of_nodes_A1, submodel_number_of_nodes_A2, aux_label, aux_id;

        double modulus, directive_length=0, dlength=0;
        double coord_X_1, coord_X_2, coord_X0_1, coord_X0_2;

        std::vector<int> IDs_A1, IDs_A2, sorted_labels_A1, sorted_labels_A2;        

        array_1d<double,3> coords_A1, coords_A2, coords_directive;
        array_1d<double,3> aux_vector_rot, aux_vector_last, aux_vector_delta, aux_vector_dir;        

        Node < 3 > Aux_Node;

        std::ofstream RES_directive("RES_directiveLine.txt", std::ios_base::out | std::ios_base::app);


        submodel_number_of_nodes_A1 = mpNodesA1->size();
        submodel_number_of_nodes_A2 = mpNodesA2->size();

        std::cout << "nº nodes " << submodel_number_of_nodes_A1 << " " << submodel_number_of_nodes_A2 << std::endl;

        // KRATOS_ERROR_IF(pNodesA1.size() != pNodesA2.size())<<"The number of nodes is not consistent between edges"<<std::endl;

        for (unsigned int inode = 0; inode < submodel_number_of_nodes_A1; inode++) {

            ModelPart::NodeIterator iteratorA1 = mpNodesA1->ptr_begin() + inode;
            Aux_Node = *iteratorA1;

            aux_id = Aux_Node.Id(); 
            IDs_A1.push_back(aux_id);
            coords_A1 = Aux_Node.Coordinates();
            coord_X_1 = Aux_Node.X();
            coord_X0_1 = Aux_Node.X0();

            aux_id = Aux_Node.Id(); 
            IDs_A2.push_back(aux_id);
            coords_A2 = Aux_Node.Coordinates();
            coord_X_2 = Aux_Node.X();
            coord_X0_2 = Aux_Node.X0();            

            coords_directive = 0.5*(coords_A1 + coords_A2);

            // RES_directive << "Node: " << coords_directive[0] << " " coords_directive[1] << " " coords_directive[2] << " "<< "\n"; 

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

        double feed_pressure=0, rotation_pressure=0, area=0.587141e-3;

        array_1d<double,3> coords,reactions;

        Node < 3 > Aux_Node;

        std::ofstream RES_reactions("RES_reactions.txt", std::ios_base::out | std::ios_base::app);


        for (unsigned int inode = 0; inode < mpNodesBC->size(); inode++) {

            ModelPart::NodeIterator iteratorBC = mpNodesBC->ptr_begin() + inode;
            Aux_Node = *iteratorBC;

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