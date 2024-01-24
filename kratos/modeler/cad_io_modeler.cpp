//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//


// Project includes
#include "cad_io_modeler.h"
#include "input_output/cad_json_input_sbm.h"
#include "input_output/cad_json_output.h"
#include "spatial_containers/bins_dynamic.h"
#include "utilities/math_utils.h"


namespace Kratos
{
    ///@name Stages
    ///@{

    void CadIoModeler::SetupGeometryModel()
    {
        // Read the refinements.iga.json
        const Parameters refinements_parameters = ReadParamatersFile("refinements.iga.json");
        
        KRATOS_ERROR_IF_NOT(mParameters.Has("cad_model_part_name"))
            << "Missing \"cad_model_part_name\" in CadIoModeler Parameters." << std::endl;
        const std::string cad_model_part_name = mParameters["cad_model_part_name"].GetString();
        ModelPart& cad_model_part = mpModel->HasModelPart(cad_model_part_name)
            ? mpModel->GetModelPart(cad_model_part_name)
            : mpModel->CreateModelPart(cad_model_part_name);
        
        const std::string DataFileName = mParameters.Has("geometry_file_name")
            ? mParameters["geometry_file_name"].GetString()
            : "geometry.cad.json";

        KRATOS_INFO_IF("::[CadIoModeler]::", mEchoLevel > 0) << "Importing Cad Model from: " << DataFileName << std::endl;

        // Create an object of CadJsonInputSBM in order to use its private method "ReadParamatersFile"
        CadJsonInputSBM<Node, Point> cadJsonInputSBM(DataFileName, mEchoLevel);
        // Need to read it in order know initial and final coordinates of the domain
        const Parameters geometry_parameters = cadJsonInputSBM.GetCadJsonParameters();

        //  READ GEOMETRY DIMENSIONS
        Vector knot_vector_u = geometry_parameters["breps"][0]["faces"][0]["surface"]["knot_vectors"][0].GetVector();
        Vector knot_vector_v = geometry_parameters["breps"][0]["faces"][0]["surface"]["knot_vectors"][1].GetVector();
        
        /// MODIFIED

        // IN ORDER TO CHECK IF YOU ARE USING TRIM OR SBM APPROACH
        std::ifstream file("txt_files/input_data.txt");
        std::string line;
        int SBM_technique;
        std::getline(file, line);
        std::getline(file, line); // Read the second line
        SBM_technique = std::stoi(line);
        file.close();
        
        // Create the Surrogate Model part
        ModelPart& surrogate_model_part = mpModel->CreateModelPart("surrogate_model_part");
        
        double knot_step_u; double knot_step_v;
        if (SBM_technique==0){
            // Create the snakes coordiantes
            CreateTheSnakeCoordinates(knot_vector_u, knot_vector_v, knot_step_u, knot_step_v, refinements_parameters, surrogate_model_part);
        }

        // surrogate_model_part.CreateNewProperties()
        cadJsonInputSBM.ReadModelPart(cad_model_part, surrogate_model_part);
    }

    void CadIoModeler::SetupModelPart()
    {
        if (mParameters.Has("output_geometry_file_name")) {
            std::string DataFileName = mParameters["output_geometry_file_name"].GetString();

            const std::string cad_model_part_name = mParameters["cad_model_part_name"].GetString();
            ModelPart& cad_model_part = mpModel->HasModelPart(cad_model_part_name)
                ? mpModel->GetModelPart(cad_model_part_name)
                : mpModel->CreateModelPart(cad_model_part_name);

            std::string output_file_text;
            CadJsonOutput::GetCadJsonOutput(cad_model_part, output_file_text, mEchoLevel);

            std::ofstream output_file(DataFileName);
            output_file << output_file_text;
            output_file.close();
        }
    }
    ///@}

    /// Reads in a json formatted file and returns its KratosParameters instance.
    Parameters CadIoModeler::ReadParamatersFile(
        const std::string& rDataFileName) const
    {
        // Check if rDataFileName ends with ".cad.json" and add it if needed.
        const std::string data_file_name = (rDataFileName.compare(rDataFileName.size() - 9, 9, ".iga.json") != 0)
            ? rDataFileName + ".iga.json"
            : rDataFileName;

        std::ifstream infile(data_file_name);
        KRATOS_ERROR_IF_NOT(infile.good()) << "Physics fil: "
            << data_file_name << " cannot be found." << std::endl;
        KRATOS_INFO_IF("ReadParamatersFile", mEchoLevel > 3)
            << "Reading file: \"" << data_file_name << "\"" << std::endl;

        std::stringstream buffer;
        buffer << infile.rdbuf();

        return Parameters(buffer.str());
    }









    /// MODIFIED
    void CadIoModeler::CreateTheSnakeCoordinates(Vector& knot_vector_u_complete, Vector& knot_vector_v_complete, double& knot_step_u, double& knot_step_v, const Parameters refinements_parameters, ModelPart& surrogate_model_part) {
        // const std::string cad_model_part_name = "skin_model_part";
        ModelPart& skin_model_part = mpModel->GetModelPart("skin_model_part");
            
        int insert_nb_per_span_u = refinements_parameters["refinements"][0]["parameters"]["insert_nb_per_span_u"].GetInt();
        int insert_nb_per_span_v = refinements_parameters["refinements"][0]["parameters"]["insert_nb_per_span_v"].GetInt();

        // Build the knot_vector without the repetitive knot values
        std::vector<double> knot_vector_u(1);
        knot_vector_u[0] = knot_vector_u_complete[0];

        std::vector<double> knot_vector_v(1);
        knot_vector_v[0] = knot_vector_v_complete[0];

        knot_step_u = (knot_vector_u_complete[knot_vector_u_complete.size()-1]-knot_vector_u[0]) / (insert_nb_per_span_u+1) ;
        knot_step_v = (knot_vector_v_complete[knot_vector_v_complete.size()-1]-knot_vector_v[0]) / (insert_nb_per_span_v+1) ;

        Vector meshSizes_uv(2);
        meshSizes_uv[0] = knot_step_u; 
        meshSizes_uv[1] = knot_step_v;
        surrogate_model_part.GetProcessInfo().SetValue(MARKER_MESHES, meshSizes_uv);

        for (IndexType j = 1; j < insert_nb_per_span_u + 1; ++j) {
            knot_vector_u.push_back(knot_vector_u[0]+ knot_step_u * j);
        }
        for (IndexType j = 1; j < insert_nb_per_span_v + 1; ++j) {
            knot_vector_v.push_back(knot_vector_v[0]+ knot_step_v *j);
        }
        knot_vector_u.push_back(knot_vector_u_complete[knot_vector_u_complete.size()-1]) ;
        knot_vector_v.push_back(knot_vector_v_complete[knot_vector_v_complete.size()-1]) ;


        // Create the matrix of active/inactive knot spans
        std::vector<std::vector<int>> knot_spans_available(insert_nb_per_span_v+1, std::vector<int>(insert_nb_per_span_u+1));

        PointVector points;
        for (auto &i_cond : skin_model_part.Conditions()) {
            points.push_back(PointTypePointer(new PointType(i_cond.Id(), i_cond.GetGeometry()[0].X(), i_cond.GetGeometry()[0].Y(), i_cond.GetGeometry()[0].Z())));
        }
        DynamicBins testBins(points.begin(), points.end());

        // Optimized Snake (Andrea)
        for (auto &i_cond : skin_model_part.Conditions()) {
            double x_true_boundary1 = i_cond.GetGeometry()[0].X();
            double y_true_boundary1 = i_cond.GetGeometry()[0].Y();
            double x_true_boundary2 = i_cond.GetGeometry()[1].X();
            double y_true_boundary2 = i_cond.GetGeometry()[1].Y();
    
            // Find the intersections of the skin boundary with the knot values
            int knot_span_u_1st_point = x_true_boundary1 / knot_step_u ;
            int knot_span_u_2nd_point = x_true_boundary2 / knot_step_u ;
            int knot_span_v_1st_point = y_true_boundary1 / knot_step_v ;
            int knot_span_v_2nd_point = y_true_boundary2 / knot_step_v ;

            SnakeStep(knot_spans_available, knot_span_u_1st_point, knot_span_u_2nd_point, knot_span_v_1st_point, knot_span_v_2nd_point,  x_true_boundary1, x_true_boundary2, y_true_boundary1, y_true_boundary2, knot_step_u, knot_step_v );
        }
        
        // Import lambda parameter: 0.0 -> External,  0.5 -> Optimal
        double lambda;
        std::ifstream file("txt_files/input_data.txt");
        std::string line;
        std::getline(file, line); std::getline(file, line); 
        std::getline(file, line); std::getline(file, line); // Read the 4th line
        lambda = std::stod(line);
        file.close();
        for (int i = 0; i < insert_nb_per_span_v+1; i++) {
            for (int j = 0; j < insert_nb_per_span_u+1; j++) {
                
                if (knot_spans_available[i][j] == 2) {
                    // Check the 4 neighbor knot spans -> Is there any completely inside? Note that we can just check 1 point.
                    if (i != 0 && i != insert_nb_per_span_v) {
                        if (knot_spans_available[i+1][j] == 0) {
                            Point gaussPoint = Point((j+0.5)*knot_step_u, (i+1+0.5)*knot_step_v, 0);
                            if (isPointInsideSkinBoundary(gaussPoint, testBins, skin_model_part)) {knot_spans_available[i+1][j] = 1;}
                        }
                        if (knot_spans_available[i-1][j] == 0) {
                            Point gaussPoint = Point((j+0.5)*knot_step_u, (i-1+0.5)*knot_step_v, 0);
                            if (isPointInsideSkinBoundary(gaussPoint, testBins, skin_model_part)) {knot_spans_available[i-1][j] = 1;}
                        }
                    }
                    if (j != 0 && j != insert_nb_per_span_u) {
                        if (knot_spans_available[i][j+1] == 0) {
                            Point gaussPoint = Point((j+1+0.5)*knot_step_u, (i+0.5)*knot_step_v, 0);
                            if (isPointInsideSkinBoundary(gaussPoint, testBins, skin_model_part)) {knot_spans_available[i][j+1] = 1;}
                        }
                        if (knot_spans_available[i][j-1] == 0) {
                            Point gaussPoint = Point((j-1+0.5)*knot_step_u, (i+0.5)*knot_step_v, 0);
                            if (isPointInsideSkinBoundary(gaussPoint, testBins, skin_model_part)) {knot_spans_available[i][j-1] = 1;}
                        }
                    }
                    // Create 25 "fake" GaussPoints to check if the majority are inside or outside
                    const int numFakeGaussPoints = 5;
                    int numberOfInsideGaussianPoints = 0;
                    for (int i_GPx = 0; i_GPx < numFakeGaussPoints; i_GPx++){
                        double x_coord = j*knot_step_u + knot_step_u/(numFakeGaussPoints+1)*(i_GPx+1);

                        // NOTE:: The v-knot spans are upside down in the matrix!!
                        for (int i_GPy = 0; i_GPy < numFakeGaussPoints; i_GPy++) 
                        {
                            double y_coord = i*knot_step_v + knot_step_v/(numFakeGaussPoints+1)*(i_GPy+1);
                            Point gaussPoint = Point(x_coord, y_coord, 0);  // GAUSSIAN POINT
                            if (isPointInsideSkinBoundary(gaussPoint, testBins, skin_model_part)) {
                                // Sum over the number of numFakeGaussPoints per knot span
                                numberOfInsideGaussianPoints++;
                            }
                        }
                        
                    }
                    // Mark the knot span as available or not depending on the number of Gauss Points Inside/Outside
                    if (numberOfInsideGaussianPoints < lambda*numFakeGaussPoints*numFakeGaussPoints) {
                        knot_spans_available[i][j] = -1; // Cut knot spans that have been checked
                    }
                    else{
                        knot_spans_available[i][j] = 1; // The knot span is considered DEACTIVE
                        }
                         
                }
            }
        }

        
        // Find the start of the Snake
        int start_i = 0;
        int start_j = 0;
        IndexType idSnakeNode = 1;
        for (int i = 0; i < (insert_nb_per_span_u+1); i++) {
            for (int j = 0; j < (insert_nb_per_span_v+1); j++ ) {
                if (knot_spans_available[j][i] == 1 ) {
                    // Check if one of the neighouring knot span is also == 1,
                    // otherwise the knot span is isolated I cannot be taken.
                    if ((knot_spans_available[j+1][i] == 1) || (knot_spans_available[j][i+1] == 1))
                    {
                        start_i = i;
                        start_j = j;
                        break;
                    }     
                }
            }
            // KRATOS_WATCH(knot_spans_available[start_j][start_i])
            if (knot_spans_available[start_j][start_i] == 1 ) { break; }
        }
        // KRATOS_WATCH(start_i)
        // KRATOS_WATCH(start_j)
        // KRATOS_WATCH(knot_spans_available)

        std::ofstream outputFile("txt_files/Snake_coordinates.txt");
        outputFile << std::setprecision(16);
        outputFile << knot_vector_u[start_i] << " " << knot_vector_v[start_j] << std::endl;
        
        Properties::Pointer p_cond_prop = surrogate_model_part.CreateNewProperties(0);
        Node snakeNode(1 , knot_vector_u[start_i], knot_vector_v[start_j], 0.0);
        
        // surrogate_model_part.AddNode(&snakeNode);

        surrogate_model_part.CreateNewNode(1, snakeNode);
        idSnakeNode++;
        // KRATOS_WATCH(knot_vector_u[start_i])
        // KRATOS_WATCH(knot_vector_v[start_j])
        // KRATOS_WATCH(start_i)
        // KRATOS_WATCH(start_j)
        // exit(0);
        
        // Follow the clockwise loop
        int end = 0;
        // We are going orizontally
        int direction = 0 ;      // 0 = up_vertical, 1 = right_orizontal, 2 = down_vertical, 3 = left_orizontal
        int i = start_i;         // Indici del percorso lungo i knot values
        int j = start_j;         // Indici del percorso lungo i knot values
        int I = start_i;         // Indici del percorso lungo i knot spans
        int J = start_j;         // Indici del percorso lungo i knot spans
        int steps = 0;
        while (end == 0 && steps < 100000) {
            steps ++ ;
            int is_special_case = 1; // Variable to check if we are in a super special case and we shouldn't go out of the while loop
            // KRATOS_WATCH(i)
            // KRATOS_WATCH(j)
            // KRATOS_WATCH(I)
            // KRATOS_WATCH(J)
            // Try to go left w.r.t. the current direction
            if (direction == 0) {I-- ; }
            if (direction == 1) {J++; }
            if (direction == 2) {I++ ; }
            if (direction == 3) {J-- ; }
            if (knot_spans_available[J][I] == 1) {
                // KRATOS_WATCH('trovato, sinistra') 
                direction = direction - 1;
                if (direction == -1) {direction = 3;}
            }
            else {
                // Need to try to go straight or right w.r.t. the current direction; risetto e muovo
                if (direction == 0) {I++ ; J++ ;}
                if (direction == 1) {J-- ; I++ ;}
                if (direction == 2) {I-- ; J-- ;}
                if (direction == 3) {J++ ; I-- ;}
                // KRATOS_WATCH(I)
                // KRATOS_WATCH(J)
                if (knot_spans_available[J][I] == 1) {
                    // KRATOS_WATCH('trovato, dritto')

                    // Stiamo andando a Dritti! -> Non scrivo nulla e muovo (i,j)
                    if (direction == 0) {j++ ; }
                    if (direction == 1) {i++ ; }
                    if (direction == 2) {j-- ; }
                    if (direction == 3) {i-- ; }

                    outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << std::endl; //DELETE

                    Node snakeNode(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                    surrogate_model_part.CreateNewNode(idSnakeNode, snakeNode);
                    Condition::Pointer pcond1 = surrogate_model_part.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                    idSnakeNode++;
                
                    // KRATOS_WATCH(knot_vector_u[i])
                    // KRATOS_WATCH(knot_vector_v[j])
                }
                else {
                    // Dobbiamo cercare a destra; resetto e muovo
                    if (direction == 0) {J-- ; I++ ;}
                    if (direction == 1) {I-- ; J-- ;}
                    if (direction == 2) {J++ ; I-- ;}
                    if (direction == 3) {I++ ; J++ ;}

                    if (knot_spans_available[J][I] == 1) {
                        // KRATOS_WATCH('trovato, destra')

                        // Stiamo andando a DX! -> Prima passo dritto, poi stampo, poi passo a destra (i,j), poi scrivo di nuovo
                        if (direction == 0) {j++ ; }
                        if (direction == 1) {i++ ; }
                        if (direction == 2) {j-- ; }
                        if (direction == 3) {i-- ; }
                        outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << std::endl; //DELETE
                        Node snakeNode1(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                        surrogate_model_part.CreateNewNode(idSnakeNode, snakeNode1);
                        Condition::Pointer pcond1 = surrogate_model_part.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                        idSnakeNode++;
                        // KRATOS_WATCH(knot_vector_u[i])
                        // KRATOS_WATCH(knot_vector_v[j])
                        if (direction == 0) {i++ ; }
                        if (direction == 1) {j-- ; }
                        if (direction == 2) {i-- ; }
                        if (direction == 3) {j++ ; }
                        outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << std::endl; //DELETE
                        Node snakeNode2(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                        surrogate_model_part.CreateNewNode(idSnakeNode, snakeNode2);
                        Condition::Pointer pcond2 = surrogate_model_part.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                        idSnakeNode++;
                        // KRATOS_WATCH(knot_vector_u[i])
                        // KRATOS_WATCH(knot_vector_v[j])
                        direction = direction + 1;
                        if (direction == 4) {direction = 0;}
                    }
                    else { // Super special case of "isolated" knot span to be circumnavigated
                        // KRATOS_WATCH('Super special case of "isolated" knot span')
                        is_special_case = 0;
                        // Resetto e muovo (I,J) "indietro"
                        if (direction == 0) {I-- ; J-- ;} 
                        if (direction == 1) {J++ ; I-- ;}
                        if (direction == 2) {I++ ; J++ ;}
                        if (direction == 3) {J-- ; I++ ;}
                        // Check
                        if (knot_spans_available[J][I] == 1) {
                            // First passo dritto, poi print, then move to right, then print again, then move to the right and print again
                            if (direction == 0) {j++ ; }
                            if (direction == 1) {i++ ; }
                            if (direction == 2) {j-- ; }
                            if (direction == 3) {i-- ; }
                            outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << std::endl; //DELETE
                            Node snakeNode(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                            surrogate_model_part.CreateNewNode(idSnakeNode, snakeNode);
                            Condition::Pointer pcond1 = surrogate_model_part.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                            idSnakeNode++;
                            // KRATOS_WATCH(knot_vector_u[i])
                            // KRATOS_WATCH(knot_vector_v[j])
                            if (direction == 0) {i++ ; }
                            if (direction == 1) {j-- ; }
                            if (direction == 2) {i-- ; }
                            if (direction == 3) {j++ ; }
                            outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << std::endl; //DELETE

                            Node snakeNode2(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                            surrogate_model_part.CreateNewNode(idSnakeNode, snakeNode2);
                            Condition::Pointer pcond2 = surrogate_model_part.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                            // surrogate_model_part.AddCondition(pcond2);
                            idSnakeNode++;
                            // KRATOS_WATCH(knot_vector_u[i])
                            // KRATOS_WATCH(knot_vector_v[j])
                            if (direction == 0) {j-- ; }
                            if (direction == 1) {i-- ; }
                            if (direction == 2) {j++ ; }
                            if (direction == 3) {i++ ; }
                            outputFile << knot_vector_u[i] << " " << knot_vector_v[j] << std::endl; //DELETE

                            Node snakeNode3(idSnakeNode , knot_vector_u[i], knot_vector_v[j], 0.0);
                            surrogate_model_part.CreateNewNode(idSnakeNode, snakeNode3);
                            Condition::Pointer pcond3 = surrogate_model_part.CreateNewCondition("LineCondition2D2N", idSnakeNode-1, {{idSnakeNode-1, idSnakeNode}}, p_cond_prop );
                            // surrogate_model_part.AddCondition(pcond3);
                            idSnakeNode++;
                            // KRATOS_WATCH(knot_vector_u[i])
                            // KRATOS_WATCH(knot_vector_v[j])
                            // Finally rotate the direction of 180 degrees
                            direction = direction + 2;
                            if (direction == 4) {direction = 0;}
                            if (direction == 5) {direction = 1;}
                        }
                        else{ KRATOS_WATCH('errore nello Snakes Coordinates'); 
                            // exit(0);
                            }  
                    }
                }
            }
            // Check if we have close the snake
            if (I == start_i && J == start_j && is_special_case == 1 ) {
                // End of the while loop
                end = 1;
                KRATOS_WATCH(steps)
                }
                

        }
        // Create "fictituos element" to memorize starting and ending node id for each surrogate boundary loop
        std::vector<ModelPart::IndexType> elem_nodes{1, idSnakeNode-1};
        surrogate_model_part.CreateNewElement("Element2D2N", 1, elem_nodes, p_cond_prop);
        // KRATOS_WATCH(surrogate_model_part.GetElement(1).GetGeometry()[0].Id())

        outputFile.close();
        KRATOS_WATCH('Snake process has finished')
        // exit(0);
    }


    // MODIFIED
    void CadIoModeler::SnakeStep(std::vector<std::vector<int>> &knot_spans_available, int knot_span_u_1st_point, int knot_span_u_2nd_point, int knot_span_v_1st_point,int knot_span_v_2nd_point,
                double& x_true_boundary1, double& x_true_boundary2, double& y_true_boundary1, double& y_true_boundary2, double& knot_step_u, double& knot_step_v) {
        
        if (knot_span_u_1st_point != knot_span_u_2nd_point || knot_span_v_1st_point != knot_span_v_2nd_point) { // INTERSECTION BETWEEN TRUE AND SURROGATE BOUNDARY

            // Check if we are jumping some cut knot spans. If yes we split the true segment
            if (std::abs(knot_span_v_1st_point-knot_span_v_2nd_point) > 1 || std::abs(knot_span_u_1st_point-knot_span_u_2nd_point) > 1 || 
                    (knot_span_u_1st_point != knot_span_u_2nd_point && knot_span_v_1st_point != knot_span_v_2nd_point) ) {
                KRATOS_WATCH('PARAMETER SPACE FINER THAN CAD IMPORTED GEOMETRY. You are splitting a true segment')
                double x_true_boundary_split = (x_true_boundary1+x_true_boundary2 ) / 2;
                double y_true_boundary_split = (y_true_boundary1+y_true_boundary2 ) / 2;
                int knot_span_u_point_split = x_true_boundary_split / knot_step_u ;
                int knot_span_v_point_split = y_true_boundary_split / knot_step_v ;
                SnakeStep(knot_spans_available, knot_span_u_1st_point,   knot_span_u_point_split, knot_span_v_1st_point  , knot_span_v_point_split, x_true_boundary1     , x_true_boundary_split, y_true_boundary1     , y_true_boundary_split, knot_step_u, knot_step_v );
                SnakeStep(knot_spans_available, knot_span_u_point_split, knot_span_u_2nd_point,   knot_span_v_point_split, knot_span_v_2nd_point,   x_true_boundary_split, x_true_boundary2,      y_true_boundary_split, y_true_boundary2     , knot_step_u, knot_step_v );
            }
            // Check if the true boundary crosses an u or a v knot value
            else if (knot_span_u_1st_point != knot_span_u_2nd_point) { // u knot value is crossed
                // Find the "knot_spans_available" using the intersection
                knot_spans_available[knot_span_v_1st_point][knot_span_u_1st_point] = 2;
                knot_spans_available[knot_span_v_1st_point][knot_span_u_2nd_point] = 2;

            }
            else if (knot_span_v_1st_point != knot_span_v_2nd_point) { // v knot value is crossed
                // Find the "knot_spans_available" using the intersection (Snake_coordinate classic -> External Boundary)
                knot_spans_available[knot_span_v_1st_point][knot_span_u_1st_point] = 2;
                knot_spans_available[knot_span_v_2nd_point][knot_span_u_1st_point] = 2;
            }
        }
        else {return;}


    }


    bool CadIoModeler::isPointInsideSkinBoundary(Point& point1, DynamicBins& testBins, ModelPart& skin_model_part)
    {
        // Get the nearest point of the true boundary
        PointerType pointToSearch = PointerType(new PointType(10000, point1.X(), point1.Y(), 0.0));
        PointerType nearestPoint = testBins.SearchNearestPoint(*pointToSearch);

        // SearchInRadius -> TOO SLOW
        // _________________________________________________________________________________________________________________________________
        // Vector meshSizes_uv = mpModel->GetModelPart("surrogate_model_part").GetProcessInfo().GetValue(MARKER_MESHES);
        // double meshSize = meshSizes_uv[0];
        // if (meshSizes_uv[1] > meshSize) {meshSize = meshSizes_uv[1];}
        // const double radius = sqrt(2.5)*(meshSize)+1e-15; 
        // const int numberOfResults = 1e6; 
        // ModelPart::NodesContainerType::ContainerType Results(numberOfResults);
        // std::vector<double> list_of_distances(numberOfResults);
        // int obtainedResults = testBins.SearchInRadius(*pointToSearch, radius, Results.begin(), list_of_distances.begin(), numberOfResults);
        // double minimum_distance=1e10;
        // int nearestNodeId;
        // for (uint i_distance = 0; i_distance < obtainedResults; i_distance++) {
        //     double new_distance = list_of_distances[i_distance];   
        //     if (new_distance < minimum_distance) { 
        //         minimum_distance = new_distance;
        //         nearestNodeId = i_distance;
        //         }
        // }
        // if (obtainedResults == 0) { KRATOS_WATCH('0 POINTS FOUND: EXIT');  exit(0);}
        // int id1 = Results[nearestNodeId]->Id();
        // _________________________________________________________________________________________________________________________________
        
        // Get the closest Condition the skin_model_part.Conditions
        int id1 = nearestPoint->Id();
        // KRATOS_WATCH(id2)
        auto nearestCondition1 = skin_model_part.GetCondition(id1);
        // Check if the condition is the first one and therefore the previous one does not exist
        int id2 = id1 - 1;
        if (id1 == 1) {
            int nConditions = skin_model_part.Conditions().size();
            id2 = nConditions; 
        }
        auto nearestCondition2 = skin_model_part.GetCondition(id2);
        // The two candidates nodes
        Point candidatePoint1 = nearestCondition1.GetGeometry()[1];
        Point candidatePoint2 = nearestCondition2.GetGeometry()[0];

        Point point2 = nearestCondition1.GetGeometry()[0]; // FIRST POINT IN TRUE GEOM
        Point point3 = candidatePoint1;// SECOND POINT IN TRUE GEOM
        if (MathUtils<double>::Norm(candidatePoint1-point1) > MathUtils<double>::Norm(candidatePoint2-point1)){
            // Need to invert the order to preserve the positivity of the area
            point3 = point2;
            point2 = candidatePoint2;
        }

        array_1d<double,3> v_1 = point2 - point1;
        array_1d<double,3> v_2 = point3 - point1;
        array_1d<double,3> crossProduct;
        MathUtils<double>::CrossProduct(crossProduct, v_1, v_2);

        bool isInside = false;
        if (crossProduct[2] > 0) {isInside = true;}

        return isInside;
    }

}