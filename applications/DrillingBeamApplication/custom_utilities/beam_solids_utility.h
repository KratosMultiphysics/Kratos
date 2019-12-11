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

        std::vector<int> sorted_labels_A;
        std::vector<int> sorted_labels_B;
        std::vector<int> sorted_labels_C; 

        std::vector<array_1d<double,3>> initial_directive_line;
        std::vector<array_1d<double,3>> directive_line; 
        std::vector<array_1d<double,3>> normals;

        /////   Default constructor   /////

        BeamSolidsUtility() {}                                            

        /////   Constructor with ModelPart   /////

        BeamSolidsUtility(ModelPart& r_model_part) {

        typedef ModelPart::NodesContainerType NodesArrayType;

        double x_A, x_B, x_C;
        // vector <double> Directive_line;        

        Node < 3 > Aux_Node;

        std::vector<int> initial_labels_A;
        std::vector<int> initial_labels_B;
        std::vector<int> initial_labels_C;

        ModelPart& EdgeLineA = r_model_part.GetSubModelPart("GENERIC_SbModel_line_A");
        ModelPart& EdgeLineB = r_model_part.GetSubModelPart("GENERIC_SbModel_line_B");
        ModelPart& EdgeLineC = r_model_part.GetSubModelPart("GENERIC_SbModel_line_C");

        NodesArrayType& pNodesA = EdgeLineA.GetCommunicator().LocalMesh().Nodes();
        NodesArrayType& pNodesB = EdgeLineB.GetCommunicator().LocalMesh().Nodes();
        NodesArrayType& pNodesC = EdgeLineC.GetCommunicator().LocalMesh().Nodes();

        unsigned int submodel_number_of_nodes = pNodesA.size();           // same number of nodes A,B,C assumed

        for (unsigned int inode = 0; inode < submodel_number_of_nodes; inode++) {

            ModelPart::NodeIterator i_iteratorA = pNodesA.ptr_begin() + inode;
            Aux_Node = *i_iteratorA;

            const int& labelA = Aux_Node.Id();
            initial_labels_A.push_back(labelA);
            x_A = Aux_Node.X();
       
            ModelPart::NodeIterator i_iteratorB = pNodesB.ptr_begin() + inode;
            Aux_Node = *i_iteratorB;

            const int& labelB = Aux_Node.Id();
            initial_labels_B.push_back(labelB);
            x_B = Aux_Node.X();

            ModelPart::NodeIterator i_iteratorC = pNodesC.ptr_begin() + inode;
            Aux_Node = *i_iteratorC;

            const int& labelC = Aux_Node.Id();
            initial_labels_C.push_back(labelC);
            x_C = Aux_Node.X();

        }

        // unsigned int values[] = { 20, 60, 40, 25, 15, 2 };

        // int compare(const void * a, const void * b ){
        //     return ( *(int*)a - *(int*)b );
        // }

        // qsort(values, 6, sizeof(int), compare);

        // std::cout << values << std::endl;

        sorted_labels_A = initial_labels_A;
        sorted_labels_B = initial_labels_B;
        sorted_labels_C = initial_labels_C;       

        }

        /////   Destructor   /////

        virtual ~BeamSolidsUtility() {}


// ------------------------------------------------------------------------------------------------

    void ComputeDirectiveLineOfBeamSolids(ModelPart& r_model_part){

        typedef ModelPart::NodesContainerType NodesArrayType;        

        double x1,y1,z1, x2,y2,z2, x3,y3,z3;
        double a,b,c;
        double modulus, directive_length=0, dlength=0;

        array_1d<double,3> coords_A, coords_B, coords_C;
        array_1d<double,3> aux_vector_rot, aux_vector_last, aux_vector_delta, aux_vector_dir;

        Node < 3 > Aux_Node;

        std::ofstream directive_RES("directive_RES.txt", std::ios_base::out | std::ios_base::app);

        const unsigned int submodel_number_of_nodes = sorted_labels_A.size();  // same number of nodes A,B,C assumed

        for (unsigned int inode = 0; inode < submodel_number_of_nodes; inode++) {

            int label_A = sorted_labels_A[inode];
            Aux_Node = r_model_part.GetNode(label_A);
            coords_A = Aux_Node.Coordinates();

            x1 = coords_A[0];
            y1 = coords_A[1];
            z1 = coords_A[2];

            int label_B = sorted_labels_B[inode];
            Aux_Node = r_model_part.GetNode(label_B);
            coords_B = Aux_Node.Coordinates();

            x2 = coords_B[0];
            y2 = coords_B[1];
            z2 = coords_B[2];

            int label_C = sorted_labels_C[inode];
            Aux_Node = r_model_part.GetNode(label_C);
            coords_C = Aux_Node.Coordinates();

            x3 = coords_C[0];
            y3 = coords_C[1];
            z3 = coords_C[2];


            aux_vector_rot[0] = (y1-y2)*(z3-z2) - (z1-z2)*(y3-y2);        // calculation of the normal to the deformed plane
            aux_vector_rot[1] = (z1-z2)*(x3-x2) - (x1-x2)*(z3-z2);
            aux_vector_rot[2] = (x1-x2)*(y3-y2) - (y1-y2)*(x3-x2);

            modulus = sqrt( pow(aux_vector_rot[0],2) + pow(aux_vector_rot[1],2) + pow(aux_vector_rot[2],2) );

            aux_vector_rot = aux_vector_rot / modulus;
            normals.push_back(aux_vector_rot);             

            aux_vector_last = aux_vector_dir;
            aux_vector_dir = (coords_A + coords_B + coords_C)/3;          // Calculation of the directive line
            a = aux_vector_dir[0];
            b = aux_vector_dir[1];
            c = aux_vector_dir[2];  

            directive_RES << a << " " << b << " " << c << "\n";            
            directive_line.push_back(aux_vector_dir);

            aux_vector_delta = aux_vector_dir - aux_vector_last;
            dlength = sqrt( pow(aux_vector_delta[0],2)+pow(aux_vector_delta[1],2)+pow(aux_vector_delta[2],2) );
            directive_length += dlength;

        }

        directive_RES << "\n Directive length ... " << directive_length << "\n";

        directive_line.clear();
        normals.clear();

    }

    void ComputeReactionsOfBeamSolids(ModelPart& r_model_part){

        typedef ModelPart::NodesContainerType NodesArrayType;
        
        double feed_pressure=0, rotation_pressure=0, area=0;
        
        array_1d<double,3> coords,reactions;

        Node < 3 > Aux_Node;


        std::ofstream reactions_RES("reactions_RES.txt", std::ios_base::out | std::ios_base::app);
        area = 5.8959e-4;

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

        reactions_RES << feed_pressure << " " << rotation_pressure << "\n";

    }

    };

} // namespace Kratos.

#endif // BEAM_SOLIDS_UTILITY_H