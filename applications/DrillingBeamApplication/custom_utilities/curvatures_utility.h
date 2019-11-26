#ifndef CURVATURES_UTILITY_H
#define CURVATURES_UTILITY_H

#include <fstream>

#include "includes/define.h"
#include "includes/variables.h"
#include "includes/model_part.h"



namespace Kratos {

    class CurvaturesUtility{

    public:

        KRATOS_CLASS_POINTER_DEFINITION(CurvaturesUtility);

        typedef ModelPart::NodesContainerType NodesArrayType;

        std::vector<int> sorted_labels_A;
        std::vector<int> sorted_labels_B;
        std::vector<int> sorted_labels_C; 

        std::vector<array_1d<double,3>> initial_directive_line;
        std::vector<array_1d<double,3>> directive_line; 
        std::vector<array_1d<double,3>> normals;

        double directive_length;      

        /// Default constructor

        CurvaturesUtility() {}

        CurvaturesUtility(ModelPart& r_model_part) {      

        typedef ModelPart::NodesContainerType NodesArrayType;

        double x_A, x_B, x_C, x0;

        array_1d<double,3> aux_vector;

        Node < 3 > Aux_Node;

        std::vector<int> initial_labels_A;
        std::vector<int> initial_labels_B;
        std::vector<int> initial_labels_C;

        unsigned int submodel_number_of_nodes = 0;


        ModelPart& EdgeLineA = r_model_part.GetSubModelPart("GENERIC_SbModel_A");
        ModelPart& EdgeLineB = r_model_part.GetSubModelPart("GENERIC_SbModel_B");
        ModelPart& EdgeLineC = r_model_part.GetSubModelPart("GENERIC_SbModel_C");

        NodesArrayType& pNodesA = EdgeLineA.GetCommunicator().LocalMesh().Nodes();
        NodesArrayType& pNodesB = EdgeLineB.GetCommunicator().LocalMesh().Nodes();
        NodesArrayType& pNodesC = EdgeLineC.GetCommunicator().LocalMesh().Nodes();

        submodel_number_of_nodes = pNodesA.size();                        // same number of nodes A,B,C assumed

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

            aux_vector[0] = ( x_A + x_B + x_C )/3;
            aux_vector[1] = 0.;
            aux_vector[2] = 0.;
            initial_directive_line.push_back(aux_vector);

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

        /// Destructor.

        virtual ~CurvaturesUtility() {}

        // vector <double> Directive_line;


// ------------------------------------------------------------------------------------------------

    void ComputeCurvatureOfBeamSolids(ModelPart& r_model_part){

        typedef ModelPart::NodesContainerType NodesArrayType;        

        double x1,y1,z1, x2,y2,z2, x3,y3,z3;
        double a,b,c;
        double modulus, dlength;

        array_1d<double,3> aux_vector,aux_vector_ini,aux_vector_dir, aux_vector_rot;      

        Node < 3 > Aux_Node;

        std::ofstream directive_RES("directive_RES.txt");
        directive_RES << "Directive line ... \n";

        std::ofstream normals_RES("normals_RES.txt");
        normals_RES << "Normals ... \n";        

        const unsigned int submodel_number_of_nodes = sorted_labels_A.size();  // same number of nodes A,B,C assumed

        for (unsigned int inode = 0; inode < submodel_number_of_nodes; inode++) {

            int label_A = sorted_labels_A[inode];
            Aux_Node = r_model_part.GetNode(label_A);
            const array_1d<double,3> & displacement_A = Aux_Node.FastGetSolutionStepValue(DISPLACEMENT);

            x1 = displacement_A[0];
            y1 = displacement_A[1];
            z1 = displacement_A[2];

            int label_B = sorted_labels_B[inode];
            Aux_Node = r_model_part.GetNode(label_B);
            const array_1d<double,3> & displacement_B = Aux_Node.FastGetSolutionStepValue(DISPLACEMENT);

            x2 = displacement_B[0];
            y2 = displacement_B[1];
            z2 = displacement_B[2];

            int label_C = sorted_labels_C[inode];
            Aux_Node = r_model_part.GetNode(label_C);
            const array_1d<double,3> & displacement_C = Aux_Node.FastGetSolutionStepValue(DISPLACEMENT);

            x3 = displacement_C[0];
            y3 = displacement_C[1];
            z3 = displacement_C[2];


            aux_vector_rot[0] = (y1-y2)*(z3-z2) - (z1-z2)*(y3-y2);        // calculation of the normal to the deformed plane
            aux_vector_rot[1] = (z1-z2)*(x3-x2) - (x1-x2)*(z3-z2);
            aux_vector_rot[2] = (x1-x2)*(y3-y2) - (y1-y2)*(x3-x2);

            modulus = sqrt( pow(aux_vector_rot[0],2) + pow(aux_vector_rot[1],2) + pow(aux_vector_rot[2],2) );

            aux_vector_rot = aux_vector_rot / modulus;
            normals.push_back(aux_vector_rot);
            a = aux_vector_rot[0];
            b = aux_vector_rot[1];
            c = aux_vector_rot[2];              
            normals_RES << a << " " << b << " " << c << "\n";            


            // aux_vector[0]  = ( x1+x2+x3 )/3;     // Definition of the deformed directive line
            // aux_vector[1]  = ( y1+y2+y3 )/3;
            // aux_vector[2]  = ( z1+z2+z3 )/3;

            // std::cout << "Pico y pala " << aux_vector << std::endl;

            aux_vector = (displacement_A + displacement_B + displacement_C)/3;
            // std::cout <<  "Directo " << aux_vector << std::endl;

            dlength = sqrt( pow(aux_vector[0],2) + pow(aux_vector[1],2) + pow(aux_vector[2],2) );
            directive_length += dlength;

            aux_vector_ini = initial_directive_line[inode];
            aux_vector_dir = aux_vector_ini + aux_vector;
            a = aux_vector_dir[0];
            b = aux_vector_dir[1];
            c = aux_vector_dir[2];                        

            directive_line.push_back(aux_vector_dir);
            directive_RES << a << " " << b << " " << c << "\n";            

        }

        directive_RES << "\n Directive length ... ";
        directive_RES << directive_length << "\n";

        directive_line.clear();
        directive_length = 0;
        normals.clear();

        directive_RES.close();        
        normals_RES.close();

    }

    };

} // namespace Kratos.

#endif // CURVATURES_UTILITY_H